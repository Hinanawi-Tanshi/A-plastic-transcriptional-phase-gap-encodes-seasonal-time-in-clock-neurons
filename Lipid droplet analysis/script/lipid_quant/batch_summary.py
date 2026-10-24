"""Read, validate, summarize, export, and plot lipid-droplet batch results.

The notebook-facing API deliberately keeps all dataframe manipulation and
plotting in this module.  Result folders are expected to contain the CSV files
written by :func:`lipid_quant.pipeline.analyze_ims`.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .paths import PathInput, normalize_paths
from .summary_diagnostics import (DEFAULT_INTENSITY_METRICS, DEFAULT_LD_SIZE_METRICS,
                                  prepare_diagnostic_results, plot_ld_size_box_scatter,
                                  size_summary_metrics)


REQUIRED_TABLES = (
    "per_brain.csv",
    "per_object.csv",
    "per_plane.csv",
    "per_detection.csv",
    "config_used.yaml",
    "run_info.json",
)

DEFAULT_BRAIN_METRICS: dict[str, str] = {
    "ld_count": "LD count",
    "ld_count_density_per_1e5_um3": "LD count density (per 1e5 um3)",
    "ld_positive_volume_um3": "LD-positive core volume (um3)",
    "ld_volume_fraction_percent": "LD-positive volume fraction (%)",
    "ld_equivalent_sphere_volume_um3": "Total equivalent-sphere volume (um3)",
    "single_plane_percent": "Single-plane accepted objects (%)",
}

DEFAULT_OBJECT_METRICS: dict[str, str] = {
    "diameter_xy_um": "Selected LoG diameter (um)",
    "peak_snr": "Peak local SNR",
    "radial_contrast_snr": "Radial contrast (SNR)",
    "n_planes": "Number of z planes",
}

DEFAULT_POSTHOC_COMPARISONS: tuple[tuple[str, str], ...] = (
    ("6586 | L12D12", "6586 | L9D15"),
    ("6586 | L12D12", "17733_DbtS | L12D12"),
    ("17733_DbtS | L12D12", "17734_DbtL | L9D15"),
)

QC_BOX_METRICS: dict[str, str] = {
    "roi_volume_1e6_um3": "ROI volume (1e6 um3)",
    "broad_hysteresis_mask_percent": "Broad mask / ROI (%)",
    "single_plane_percent": "Single-plane accepted objects (%)",
}

_SAMPLE_PATTERN = re.compile(
    r"^(?P<code>[^_]+)(?:_(?P<variant>Dbt[^_]+))?.*?"
    r"_(?P<cycle>L\d+D\d+)_brain(?P<brain>\d+)(?:_|$)",
    flags=re.IGNORECASE,
)


@dataclass
class BatchResults:
    root: Path
    manifest: pd.DataFrame
    brains: pd.DataFrame
    objects: pd.DataFrame
    planes: pd.DataFrame
    detections: pd.DataFrame
    group_order: list[str]
    roots: list[Path] = field(default_factory=list)


@dataclass
class BatchReport:
    summary: pd.DataFrame
    statistics: pd.DataFrame
    validation: dict[str, pd.DataFrame]
    exported_tables: dict[str, Path]
    figures: dict[str, Path]


def discover_result_runs(
    root: PathInput,
    required_tables: Sequence[str] = REQUIRED_TABLES,
    *,
    recursive: bool = False,
) -> pd.DataFrame:
    """Discover one or multiple roots, or explicit run folders, by full path."""
    records: list[dict[str, Any]] = []
    candidates: dict[str, Path] = {}
    for batch_root in normalize_paths(root):
        if not batch_root.is_dir():
            raise NotADirectoryError(batch_root)
        if (batch_root / "per_brain.csv").is_file():
            directories = [batch_root]
        else:
            directories = batch_root.rglob("*") if recursive else batch_root.iterdir()
        for directory in directories:
            if directory.is_dir():
                candidates.setdefault(str(directory.resolve()).casefold(), directory.resolve())
    for directory in sorted(candidates.values()):
        present = {name: (directory / name).is_file() for name in required_tables}
        missing = [name for name, exists in present.items() if not exists]
        # Ignore report/output folders that are not analysis runs at all.
        if len(missing) == len(required_tables):
            continue
        record: dict[str, Any] = {
            "result_dir": directory.name,
            "result_path": str(directory.resolve()),
            "batch_root": str(directory.parent),
            "batch_id": directory.parent.name,
            "modified_time": directory.stat().st_mtime,
            "complete": not missing,
            "missing_files": "; ".join(missing),
        }
        for name, exists in present.items():
            record[f"has_{Path(name).stem}"] = bool(exists)
        records.append(record)
    if not records:
        raise FileNotFoundError(f"No analysis result folders found in: {root}")
    return pd.DataFrame.from_records(records)


def _parse_sample_metadata(
    brain_id: str,
    genotype_labels: Mapping[str, str] | None,
) -> dict[str, Any]:
    match = _SAMPLE_PATTERN.search(str(brain_id))
    if match is None:
        code = str(brain_id).split("_", 1)[0]
        variant = ""
        cycle = "unparsed"
        brain_number: Any = pd.NA
    else:
        code = match.group("code")
        variant = match.group("variant") or ""
        cycle = match.group("cycle").upper()
        brain_number = int(match.group("brain"))
    default_label = f"{code}_{variant}" if variant else code
    genotype = (genotype_labels or {}).get(code, default_label)
    return {
        "genotype_code": code,
        "variant": variant,
        "genotype": genotype,
        "light_cycle": cycle,
        "brain_number": brain_number,
        "condition": f"{genotype} | {cycle}",
    }


def _add_sample_metadata(
    frame: pd.DataFrame,
    genotype_labels: Mapping[str, str] | None,
) -> pd.DataFrame:
    result = frame.copy()
    metadata = pd.DataFrame.from_records(
        [_parse_sample_metadata(value, genotype_labels) for value in result["brain_id"]],
        index=result.index,
    )
    for column in metadata.columns:
        result[column] = metadata[column]
    result["brain_number"] = result["brain_number"].astype("Int64")
    return result


def _derived_brain_metrics(brains: pd.DataFrame) -> pd.DataFrame:
    result = brains.copy()
    result["roi_volume_1e6_um3"] = result["roi_volume_um3"] / 1e6
    result["ld_positive_volume_1e3_um3"] = result["ld_positive_volume_um3"] / 1e3
    result["ld_volume_fraction_percent"] = 100.0 * result["ld_volume_fraction"]
    result["broad_hysteresis_mask_percent"] = (
        100.0 * result["broad_hysteresis_mask_fraction"]
    )
    result["single_plane_percent"] = 100.0 * result["single_plane_object_fraction"]
    result["saturation_percent"] = 100.0 * result["saturation_fraction"]
    denominator = result["ld_count"].replace(0, np.nan)
    result["mean_equivalent_sphere_volume_um3"] = (
        result["ld_equivalent_sphere_volume_um3"] / denominator
    )
    result["mean_positive_core_volume_um3"] = (
        result["ld_positive_volume_um3"] / denominator
    )
    return result


def _normalize_group_order(
    conditions: pd.Series,
    group_order: Sequence[str] | None,
    *,
    include_unlisted_groups: bool = True,
) -> list[str]:
    if not isinstance(include_unlisted_groups, (bool, np.bool_)):
        raise ValueError("include_unlisted_groups must be True or False")
    present = [str(value) for value in conditions.dropna().unique()]
    if group_order is None:
        if not include_unlisted_groups:
            raise ValueError("Provide group_order when include_unlisted_groups=False")
        return sorted(present)
    ordered = list(dict.fromkeys(str(value) for value in group_order if str(value) in present))
    if include_unlisted_groups:
        ordered.extend(value for value in present if value not in ordered)
    elif not ordered:
        raise ValueError(
            "No groups match group_order with include_unlisted_groups=False. "
            f"Available groups: {sorted(present)}"
        )
    return ordered


def _categorize_conditions(frame: pd.DataFrame, order: Sequence[str]) -> pd.DataFrame:
    result = frame.copy()
    if "condition" in result:
        result["condition"] = pd.Categorical(
            result["condition"].astype(str),
            categories=list(order),
            ordered=True,
        )
    return result


def _read_child_table(
    run_path: Path,
    filename: str,
    brain_metadata: Mapping[str, Any],
) -> pd.DataFrame:
    frame = pd.read_csv(run_path / filename)
    frame.insert(0, "brain_id", str(brain_metadata["brain_id"]))
    frame.insert(1, "result_dir", run_path.name)
    frame["result_path"] = str(run_path.resolve())
    frame["batch_root"] = brain_metadata["batch_root"]
    frame["batch_id"] = brain_metadata["batch_id"]
    for column in ("genotype_code", "variant", "genotype", "light_cycle", "condition"):
        frame[column] = brain_metadata[column]
    frame["brain_number"] = brain_metadata["brain_number"]
    return frame


def load_batch_results(
    root: PathInput,
    genotype_labels: Mapping[str, str] | None = None,
    group_order: Sequence[str] | None = None,
    *,
    include_unlisted_groups: bool = True,
    load_objects: bool = True,
    accepted_objects_only: bool = True,
    load_planes: bool = True,
    load_detections: bool = False,
    duplicate_policy: str = "error",
    recursive: bool = False,
) -> BatchResults:
    """Load complete result folders and attach parsed biological metadata.

    ``duplicate_policy`` may be ``"error"``, ``"latest"``, or ``"all"``.
    Detections are disabled by default because they are the largest table and
    are not required for standard biological summaries.

    With ``include_unlisted_groups=True`` (default), groups absent from
    ``group_order`` are appended as before. False selects only listed groups
    for ALL child tables, plots, summaries and statistical tests; it does not
    modify individual run files. The manifest retains excluded runs with
    ``selected=False``. An explicit, nonempty matching group_order is required.
    """
    roots = normalize_paths(root)
    manifest = discover_result_runs(roots, recursive=recursive)
    incomplete = manifest.loc[~manifest["complete"]]
    if not incomplete.empty:
        details = incomplete[["result_dir", "missing_files"]].to_string(index=False)
        raise ValueError(f"Incomplete analysis result folders:\n{details}")

    brain_frames: list[pd.DataFrame] = []
    for row in manifest.itertuples(index=False):
        run_path = Path(row.result_path)
        brain = pd.read_csv(run_path / "per_brain.csv")
        if len(brain) != 1:
            raise ValueError(f"Expected one row in {run_path / 'per_brain.csv'}")
        brain["result_dir"] = run_path.name
        brain["result_path"] = str(run_path.resolve())
        brain["batch_root"] = row.batch_root
        brain["batch_id"] = row.batch_id
        info = json.loads((run_path / "run_info.json").read_text(encoding="utf-8-sig"))
        brain["source_file"] = info.get("input_file", "")
        brain["source_directory"] = str(Path(info["input_file"]).parent) if info.get("input_file") else ""
        brain["roi_method"] = info.get("roi", {}).get("method", "legacy_or_unspecified")
        roi_info=info.get('roi',{})
        scope=roi_info.get('roi_scope')
        if scope is None:
            scope='central_brain' if roi_info.get('mode') in {'auto_central_brain','guided_central_brain'} else roi_info.get('mode','not_recorded')
        brain['roi_scope']=scope
        brain["roi_status"] = info.get("roi", {}).get("roi_status", "not_recorded")
        brain["roi_review_flags"] = "; ".join(info.get("roi", {}).get("review_flags", []))
        brain["roi_algorithm_version"] = info.get("roi", {}).get("algorithm_version", "not_recorded")
        brain["roi_algorithm_source_hash"] = info.get("roi", {}).get("algorithm_source_hash", "not_recorded")
        brain["result_modified_time"] = float(row.modified_time)
        brain_frames.append(brain)
    brains = pd.concat(brain_frames, ignore_index=True)
    brains = _add_sample_metadata(brains, genotype_labels)
    order = _normalize_group_order(
        brains["condition"], group_order, include_unlisted_groups=include_unlisted_groups)
    if not include_unlisted_groups:
        # Filter before duplicate checks and before reading large child tables.
        # Hidden groups must not enter a Tukey fit or its correction family.
        brains = brains.loc[brains["condition"].isin(order)].copy().reset_index(drop=True)
    brains = _derived_brain_metrics(brains)

    duplicates = brains["brain_id"].duplicated(keep=False)
    duplicate_policy = duplicate_policy.lower()
    if duplicates.any():
        if duplicate_policy == "error":
            names = sorted(brains.loc[duplicates, "brain_id"].unique())
            raise ValueError(
                "Duplicate brain_id values were found. Set duplicate_policy='latest' "
                f"to retain the newest run. Duplicates: {names}"
            )
        if duplicate_policy == "latest":
            brains = (
                brains.sort_values("result_modified_time")
                .drop_duplicates("brain_id", keep="last")
                .reset_index(drop=True)
            )
        elif duplicate_policy != "all":
            raise ValueError("duplicate_policy must be: error, latest, or all")
    elif duplicate_policy not in {"error", "latest", "all"}:
        raise ValueError("duplicate_policy must be: error, latest, or all")

    brains = _categorize_conditions(brains, order)
    selected_runs = {str(value) for value in brains["result_path"]}
    metadata_by_run = {
        str(row.result_path): row._asdict()
        for row in brains.itertuples(index=False)
    }

    object_frames: list[pd.DataFrame] = []
    plane_frames: list[pd.DataFrame] = []
    detection_frames: list[pd.DataFrame] = []
    manifest["selected"] = manifest["result_path"].isin(selected_runs)
    for run_path_string in sorted(selected_runs):
        run_path = Path(run_path_string)
        metadata = metadata_by_run[run_path_string]
        if load_objects:
            objects = _read_child_table(run_path, "per_object.csv", metadata)
            if accepted_objects_only and "accepted" in objects:
                objects = objects.loc[objects["accepted"].astype(bool)].copy()
            object_frames.append(objects)
        if load_planes:
            planes = _read_child_table(run_path, "per_plane.csv", metadata)
            roi_area = planes["roi_area_um2"].replace(0, np.nan)
            planes["foreground_fraction"] = planes["foreground_area_um2"] / roi_area
            planes["candidate_density_per_1e3_um2"] = (
                planes["candidate_count"] / roi_area * 1e3
            )
            maximum_roi = float(planes["roi_area_um2"].max())
            planes["roi_area_relative"] = (
                planes["roi_area_um2"] / maximum_roi if maximum_roi > 0 else np.nan
            )
            active = planes["roi_area_um2"] > 0
            if active.any():
                z_min = float(planes.loc[active, "z_um"].min())
                z_max = float(planes.loc[active, "z_um"].max())
                if z_max > z_min:
                    planes["normalized_depth"] = (planes["z_um"] - z_min) / (z_max - z_min)
                else:
                    planes["normalized_depth"] = 0.5
            else:
                planes["normalized_depth"] = np.nan
            plane_frames.append(planes)
        if load_detections:
            detection_frames.append(
                _read_child_table(run_path, "per_detection.csv", metadata)
            )

    objects = pd.concat(object_frames, ignore_index=True) if object_frames else pd.DataFrame()
    planes = pd.concat(plane_frames, ignore_index=True) if plane_frames else pd.DataFrame()
    detections = (
        pd.concat(detection_frames, ignore_index=True)
        if detection_frames
        else pd.DataFrame()
    )
    objects = _categorize_conditions(objects, order)
    planes = _categorize_conditions(planes, order)
    detections = _categorize_conditions(detections, order)
    return BatchResults(
        root=roots[0],
        manifest=manifest,
        brains=brains,
        objects=objects,
        planes=planes,
        detections=detections,
        group_order=order,
        roots=roots,
    )


def validate_batch_results(results: BatchResults) -> dict[str, pd.DataFrame]:
    """Perform structural and count-consistency checks."""
    brains = results.brains
    issues: list[dict[str, str]] = []
    if brains["brain_id"].duplicated().any():
        for brain_id in brains.loc[brains["brain_id"].duplicated(False), "brain_id"].unique():
            issues.append({"severity": "error", "issue": f"duplicate brain_id: {brain_id}"})
    if (brains["roi_volume_um3"] <= 0).any():
        for brain_id in brains.loc[brains["roi_volume_um3"] <= 0, "brain_id"]:
            issues.append({"severity": "error", "issue": f"non-positive ROI volume: {brain_id}"})
    if not results.objects.empty:
        if "accepted" in results.objects:
            accepted = results.objects.loc[results.objects["accepted"].astype(bool)]
        else:
            accepted = results.objects
        observed = accepted.groupby("result_path", observed=True).size()
        expected = brains.set_index("result_path")["ld_count"].astype(int)
        mismatch = expected.ne(observed.reindex(expected.index, fill_value=0))
        for brain_id in expected.index[mismatch]:
            issues.append(
                {
                    "severity": "error",
                    "issue": (
                        f"accepted-object count mismatch: {brain_id}; "
                        f"per_brain={expected[brain_id]}, objects={observed.get(brain_id, 0)}"
                    ),
                }
            )
    if not results.planes.empty:
        observed_planes = results.planes.groupby("result_path", observed=True).size()
        expected_planes = brains.set_index("result_path")["shape_z"].astype(int)
        mismatch = expected_planes.ne(observed_planes.reindex(expected_planes.index, fill_value=0))
        for brain_id in expected_planes.index[mismatch]:
            issues.append(
                {
                    "severity": "error",
                    "issue": (
                        f"plane count mismatch: {brain_id}; expected={expected_planes[brain_id]}, "
                        f"observed={observed_planes.get(brain_id, 0)}"
                    ),
                }
            )

    group_counts = (
        brains.groupby(
            ["condition", "genotype", "light_cycle"],
            observed=True,
            dropna=False,
        )
        .size()
        .rename("n_brains")
        .reset_index()
    )
    config_hashes = (
        brains.groupby("config_hash", dropna=False)
        .size()
        .rename("n_brains")
        .reset_index()
        .sort_values("n_brains", ascending=False)
    )
    batch_counts = brains.groupby(["batch_id", "batch_root", "condition"], observed=True).size().rename("n_brains").reset_index()
    if brains["roi_method"].nunique() > 1:
        issues.append({"severity": "warning", "issue": "Multiple ROI methods loaded. Do not pool unvalidated ROI methods in one comparison"})
    if 'roi_scope' in brains and brains['roi_scope'].nunique()>1:
        issues.append({'severity':'warning','issue':'Mixed ROI scopes (whole brain / central brain / other). Counts and density denominators are not directly comparable; do not pool these results'})
    if "roi_status" in brains:
        pending = brains["roi_status"].isin(["needs_review", "failed"])
        if pending.any():
            issues.append({"severity": "warning", "issue": f"{int(pending.sum())} ROI(s) need review or failed. They remain included; explicitly review before inference"})
    if "roi_algorithm_source_hash" in brains and brains.loc[brains["roi_algorithm_source_hash"] != "not_recorded", "roi_algorithm_source_hash"].nunique() > 1:
        issues.append({"severity": "warning", "issue": "Multiple ROI implementation hashes loaded. Check algorithm-version compatibility before pooling"})
    overview = pd.DataFrame(
        [
            {
                "batch_root": "; ".join(map(str, results.roots or [results.root])),
                "n_input_roots": len(results.roots or [results.root]),
                "n_result_folders": int(len(results.manifest)),
                "n_brains_loaded": int(len(brains)),
                "n_conditions": int(brains["condition"].nunique()),
                "n_objects_loaded": int(len(results.objects)),
                "n_planes_loaded": int(len(results.planes)),
                "n_detections_loaded": int(len(results.detections)),
                "n_config_hashes": int(brains["config_hash"].nunique()),
                "n_issues": int(len(issues)),
            }
        ]
    )
    issue_frame = pd.DataFrame.from_records(issues, columns=["severity", "issue"])
    if issue_frame.empty:
        issue_frame = pd.DataFrame([{"severity": "ok", "issue": "No validation issues found"}])
    return {
        "overview": overview,
        "group_counts": group_counts,
        "batch_counts": batch_counts,
        "config_hashes": config_hashes,
        "issues": issue_frame,
    }


def summarize_brains(
    brains: pd.DataFrame,
    metrics: Mapping[str, str] | Sequence[str] = DEFAULT_BRAIN_METRICS,
    group_column: str = "condition",
) -> pd.DataFrame:
    """Return long-form descriptive statistics, treating each brain as n=1."""
    metric_names = list(metrics.keys()) if isinstance(metrics, Mapping) else list(metrics)
    missing = [metric for metric in metric_names if metric not in brains]
    if missing:
        raise KeyError(f"Brain metrics not found: {missing}")
    records: list[dict[str, Any]] = []
    for group, group_frame in brains.groupby(group_column, observed=True, sort=False):
        for metric in metric_names:
            values = pd.to_numeric(group_frame[metric], errors="coerce").dropna()
            n = int(len(values))
            mean = float(values.mean()) if n else np.nan
            standard_deviation = float(values.std(ddof=1)) if n > 1 else np.nan
            records.append(
                {
                    group_column: str(group),
                    "metric": metric,
                    "label": metrics.get(metric, metric) if isinstance(metrics, Mapping) else metric,
                    "n": n,
                    "mean": mean,
                    "std": standard_deviation,
                    "sem": standard_deviation / math.sqrt(n) if n > 1 else np.nan,
                    "median": float(values.median()) if n else np.nan,
                    "q25": float(values.quantile(0.25)) if n else np.nan,
                    "q75": float(values.quantile(0.75)) if n else np.nan,
                    "minimum": float(values.min()) if n else np.nan,
                    "maximum": float(values.max()) if n else np.nan,
                }
            )
    return pd.DataFrame.from_records(records)


def _pvalue_to_stars(p_value: float) -> str:
    if not np.isfinite(p_value):
        return ""
    if p_value < 0.0001:
        return "****"
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def _has_numerical_variation(values: np.ndarray) -> bool:
    """Do not mistake roundoff in a constant array for biological variance."""
    if len(values) < 2:
        return False
    scale = max(float(np.max(np.abs(values))), np.finfo(float).tiny)
    return bool(np.ptp(values) > 64 * np.finfo(float).eps * scale)


def calculate_anova_tukey(
    frame: pd.DataFrame,
    metrics: Mapping[str, str] | Sequence[str] = DEFAULT_BRAIN_METRICS,
    *,
    group_column: str = "condition",
    group_order: Sequence[str] | None = None,
    selected_comparisons: Sequence[tuple[str, str]] = DEFAULT_POSTHOC_COMPARISONS,
) -> pd.DataFrame:
    """Run one-way ANOVA and all-pairs Tukey HSD for each per-brain metric.

    Tukey HSD adjusted p-values control the family-wise error rate across all
    pairwise comparisons for a metric. ``selected_comparisons`` only flags the
    contrasts intended for plotting; it does not limit the fitted post-hoc test.
    """
    from scipy.stats import f_oneway
    from statsmodels.stats.multicomp import pairwise_tukeyhsd

    metric_names = list(metrics.keys()) if isinstance(metrics, Mapping) else list(metrics)
    metric_labels = metrics if isinstance(metrics, Mapping) else {}
    missing = [metric for metric in metric_names if metric not in frame]
    if missing:
        raise KeyError(f"Metrics not found for ANOVA/Tukey HSD: {missing}")
    if group_column not in frame:
        raise KeyError(f"Group column not found: {group_column}")

    present_groups = [str(value) for value in frame[group_column].dropna().unique()]
    if group_order is None:
        order = sorted(present_groups)
    else:
        order = [str(value) for value in group_order if str(value) in present_groups]
        order.extend(value for value in present_groups if value not in order)
    selected_keys = {frozenset((str(left), str(right))) for left, right in selected_comparisons}

    records: list[dict[str, Any]] = []
    for metric in metric_names:
        clean = frame[[group_column, metric]].copy()
        clean[group_column] = clean[group_column].astype(str)
        clean[metric] = pd.to_numeric(clean[metric], errors="coerce")
        clean = clean.loc[np.isfinite(clean[metric])]
        values_by_group = {
            group: clean.loc[clean[group_column] == group, metric].to_numpy(float)
            for group in order
        }
        values_by_group = {
            group: values for group, values in values_by_group.items() if values.size > 0
        }
        if len(values_by_group) < 2:
            continue

        total_n = sum(len(values) for values in values_by_group.values())
        test_note = ""
        if any(_has_numerical_variation(values) for values in values_by_group.values()):
            anova = f_oneway(*values_by_group.values())
            anova_f, anova_p = float(anova.statistic), float(anova.pvalue)
            endog = np.concatenate(list(values_by_group.values()))
            groups = np.concatenate(
                [np.repeat(group, len(values)) for group, values in values_by_group.items()]
            )
            tukey = pairwise_tukeyhsd(endog=endog, groups=groups, alpha=0.05)
            tukey_groups = [str(value) for value in tukey.groupsunique]
        else:
            # Discrete LoG medians can be exactly identical within every group.
            # A zero residual variance cannot support ANOVA/Tukey inference.
            anova_f = anova_p = np.nan
            tukey = None
            tukey_groups = sorted(values_by_group)
            test_note = "No within-group variation; ANOVA/Tukey is not estimable"
        tukey_pairs = list(combinations(tukey_groups, 2))
        for index, (group1, group2) in enumerate(tukey_pairs):
            p_adjusted = float(tukey.pvalues[index]) if tukey is not None else np.nan
            values1 = values_by_group[group1]
            values2 = values_by_group[group2]
            records.append(
                {
                    "metric": metric,
                    "label": metric_labels.get(metric, metric),
                    "test_method": "anova_tukey",
                    "omnibus_test": "one_way_anova",
                    "omnibus_statistic": anova_f,
                    "omnibus_p": anova_p,
                    "pairwise_test": "tukey_hsd",
                    "pairwise_statistic": np.nan,
                    "pairwise_df": np.nan,
                    "p_raw": np.nan,
                    "p_adjusted": p_adjusted,
                    "p_adjust_method": "tukey_hsd",
                    "anova_f": anova_f,
                    "anova_df_between": int(len(values_by_group) - 1),
                    "anova_df_within": int(total_n - len(values_by_group)),
                    "anova_p": anova_p,
                    "group1": group1,
                    "group2": group2,
                    "n_group1": int(len(values1)),
                    "n_group2": int(len(values2)),
                    "mean_group1": float(np.mean(values1)),
                    "mean_group2": float(np.mean(values2)),
                    "mean_difference_group2_minus_group1": float(tukey.meandiffs[index]) if tukey is not None else float(np.mean(values2) - np.mean(values1)),
                    "tukey_p_adjusted": p_adjusted,
                    "ci_low": float(tukey.confint[index, 0]) if tukey is not None else np.nan,
                    "ci_high": float(tukey.confint[index, 1]) if tukey is not None else np.nan,
                    "reject_0_05": bool(tukey.reject[index]) if tukey is not None else False,
                    "significance": _pvalue_to_stars(p_adjusted),
                    "test_note": test_note,
                    "display_on_plot": frozenset((group1, group2)) in selected_keys,
                }
            )
    return pd.DataFrame.from_records(records)


def calculate_pairwise_ttests(
    frame: pd.DataFrame,
    metrics: Mapping[str, str] | Sequence[str] = DEFAULT_BRAIN_METRICS,
    *,
    group_column: str = "condition",
    group_order: Sequence[str] | None = None,
    selected_comparisons: Sequence[tuple[str, str]] = DEFAULT_POSTHOC_COMPARISONS,
    equal_var: bool = False,
    p_adjust_method: str = "holm",
) -> pd.DataFrame:
    """Run all pairwise independent t-tests for each per-brain metric.

    Welch's t-test is used when ``equal_var=False`` (the default). Raw p-values
    are adjusted across all pairwise comparisons within each metric. Any method
    accepted by :func:`statsmodels.stats.multitest.multipletests` can be used;
    ``"none"`` leaves p-values unadjusted.
    """
    from scipy.stats import t as student_t
    from scipy.stats import ttest_ind_from_stats
    from statsmodels.stats.multitest import multipletests

    metric_names = list(metrics.keys()) if isinstance(metrics, Mapping) else list(metrics)
    metric_labels = metrics if isinstance(metrics, Mapping) else {}
    missing = [metric for metric in metric_names if metric not in frame]
    if missing:
        raise KeyError(f"Metrics not found for pairwise t-tests: {missing}")
    if group_column not in frame:
        raise KeyError(f"Group column not found: {group_column}")

    present_groups = [str(value) for value in frame[group_column].dropna().unique()]
    if group_order is None:
        order = sorted(present_groups)
    else:
        order = [str(value) for value in group_order if str(value) in present_groups]
        order.extend(value for value in present_groups if value not in order)
    selected_keys = {frozenset((str(left), str(right))) for left, right in selected_comparisons}
    adjustment = str(p_adjust_method).lower()

    records: list[dict[str, Any]] = []
    for metric in metric_names:
        clean = frame[[group_column, metric]].copy()
        clean[group_column] = clean[group_column].astype(str)
        clean[metric] = pd.to_numeric(clean[metric], errors="coerce")
        clean = clean.loc[np.isfinite(clean[metric])]
        values_by_group = {
            group: clean.loc[clean[group_column] == group, metric].to_numpy(float)
            for group in order
        }
        values_by_group = {
            group: values for group, values in values_by_group.items() if values.size > 0
        }
        if len(values_by_group) < 2:
            continue

        metric_records: list[dict[str, Any]] = []
        for group1, group2 in combinations(values_by_group, 2):
            values1 = values_by_group[group1]
            values2 = values_by_group[group2]
            variation1 = _has_numerical_variation(values1)
            variation2 = _has_numerical_variation(values2)
            mean1 = float(np.mean(values1)) if variation1 else float(values1[0])
            mean2 = float(np.mean(values2)) if variation2 else float(values2[0])
            test_note = ""
            if len(values1) < 2 or len(values2) < 2:
                test_note = "Fewer than two brains in a compared group"
            elif not variation1 and not variation2:
                test_note = "No within-group variation in either group; t-test is not estimable"
            if test_note:
                statistic = np.nan
                p_raw = np.nan
                degrees_freedom = np.nan
                ci_low = np.nan
                ci_high = np.nan
            else:
                # Exact zero SD for a constant group avoids catastrophic
                # cancellation inside scipy's moment calculation. For ordinary
                # groups this is the same independent t-test from sample stats.
                variance1 = float(np.var(values1, ddof=1)) if variation1 else 0.0
                variance2 = float(np.var(values2, ddof=1)) if variation2 else 0.0
                test = ttest_ind_from_stats(mean1, math.sqrt(variance1), len(values1),
                                           mean2, math.sqrt(variance2), len(values2), equal_var=equal_var)
                statistic = float(test.statistic)
                p_raw = float(test.pvalue)
                if equal_var:
                    degrees_freedom = float(len(values1) + len(values2) - 2)
                    pooled_variance = (
                        (len(values1) - 1) * variance1 + (len(values2) - 1) * variance2
                    ) / degrees_freedom
                    standard_error = math.sqrt(
                        pooled_variance * (1.0 / len(values1) + 1.0 / len(values2))
                    )
                else:
                    component1 = variance1 / len(values1)
                    component2 = variance2 / len(values2)
                    standard_error = math.sqrt(component1 + component2)
                    denominator = (
                        component1**2 / (len(values1) - 1)
                        + component2**2 / (len(values2) - 1)
                    )
                    degrees_freedom = (
                        (component1 + component2) ** 2 / denominator
                        if denominator > 0
                        else np.nan
                    )
                mean_difference = mean2 - mean1
                critical_value = (
                    float(student_t.ppf(0.975, degrees_freedom))
                    if np.isfinite(degrees_freedom)
                    else np.nan
                )
                margin = critical_value * standard_error
                ci_low = mean_difference - margin
                ci_high = mean_difference + margin
            metric_records.append(
                {
                    "metric": metric,
                    "label": metric_labels.get(metric, metric),
                    "test_method": "pairwise_ttest",
                    "omnibus_test": "",
                    "omnibus_statistic": np.nan,
                    "omnibus_p": np.nan,
                    "pairwise_test": "student_ttest" if equal_var else "welch_ttest",
                    "pairwise_statistic": statistic,
                    "pairwise_df": degrees_freedom,
                    "p_raw": p_raw,
                    "p_adjust_method": adjustment,
                    "anova_f": np.nan,
                    "anova_df_between": np.nan,
                    "anova_df_within": np.nan,
                    "anova_p": np.nan,
                    "group1": group1,
                    "group2": group2,
                    "n_group1": int(len(values1)),
                    "n_group2": int(len(values2)),
                    "mean_group1": mean1,
                    "mean_group2": mean2,
                    "mean_difference_group2_minus_group1": mean2 - mean1,
                    "tukey_p_adjusted": np.nan,
                    "ci_low": float(ci_low),
                    "ci_high": float(ci_high),
                    "test_note": test_note,
                    "display_on_plot": frozenset((group1, group2)) in selected_keys,
                }
            )

        raw_p_values = np.asarray([row["p_raw"] for row in metric_records], dtype=float)
        adjusted_p_values = np.full(raw_p_values.shape, np.nan, dtype=float)
        finite = np.isfinite(raw_p_values)
        if finite.any():
            if adjustment in {"none", "raw", "unadjusted"}:
                adjusted_p_values[finite] = raw_p_values[finite]
            else:
                adjusted_p_values[finite] = multipletests(
                    raw_p_values[finite],
                    alpha=0.05,
                    method=adjustment,
                )[1]
        for row, p_adjusted in zip(metric_records, adjusted_p_values):
            row["p_adjusted"] = float(p_adjusted)
            row["reject_0_05"] = bool(np.isfinite(p_adjusted) and p_adjusted < 0.05)
            row["significance"] = _pvalue_to_stars(float(p_adjusted))
        records.extend(metric_records)
    return pd.DataFrame.from_records(records)


def calculate_statistical_tests(
    frame: pd.DataFrame,
    metrics: Mapping[str, str] | Sequence[str] = DEFAULT_BRAIN_METRICS,
    *,
    method: str = "anova_tukey",
    group_column: str = "condition",
    group_order: Sequence[str] | None = None,
    selected_comparisons: Sequence[tuple[str, str]] = DEFAULT_POSTHOC_COMPARISONS,
    ttest_equal_var: bool = False,
    ttest_p_adjust_method: str = "holm",
) -> pd.DataFrame:
    """Dispatch to the selected per-brain group-comparison workflow."""
    normalized_method = str(method).strip().lower().replace("-", "_")
    if normalized_method in {"anova_tukey", "anova+tukey", "tukey"}:
        return calculate_anova_tukey(
            frame,
            metrics,
            group_column=group_column,
            group_order=group_order,
            selected_comparisons=selected_comparisons,
        )
    if normalized_method in {"pairwise_ttest", "multiple_ttest", "ttest", "t_test"}:
        return calculate_pairwise_ttests(
            frame,
            metrics,
            group_column=group_column,
            group_order=group_order,
            selected_comparisons=selected_comparisons,
            equal_var=ttest_equal_var,
            p_adjust_method=ttest_p_adjust_method,
        )
    raise ValueError(
        "Unknown statistical test method. Use 'anova_tukey' or 'pairwise_ttest'."
    )


def display_batch_overview(
    results: BatchResults,
    summary: pd.DataFrame | None = None,
    validation: dict[str, pd.DataFrame] | None = None,
) -> None:
    """Display compact tables in Jupyter without notebook-side manipulation."""
    from IPython.display import Markdown, display

    validation = validation or validate_batch_results(results)
    summary = summary if summary is not None else summarize_brains(results.brains)
    display(Markdown("### Batch validation"))
    display(validation["overview"])
    display(validation["issues"])
    display(Markdown("### Group sizes"))
    display(validation["group_counts"])
    display(Markdown("### Input batches (pooling does not adjust batch effects)"))
    display(validation["batch_counts"])
    display(Markdown("### Config hashes"))
    display(validation["config_hashes"])
    display(Markdown("### Per-brain table"))
    columns = [
        "batch_id",
        "roi_method",
        "roi_status",
        "roi_scope",
        "roi_review_flags",
        "condition",
        "brain_number",
        "roi_volume_um3",
        "ld_count",
        "ld_count_density_per_1e5_um3",
        "ld_positive_volume_um3",
        "ld_volume_fraction_percent",
        "single_plane_percent",
    ]
    display(results.brains[columns].sort_values(["condition", "brain_number"]))
    display(Markdown("### Group summary (long format)"))
    display(summary)


def export_batch_tables(
    results: BatchResults,
    summary: pd.DataFrame,
    validation: dict[str, pd.DataFrame],
    output_dir: str | Path,
    *,
    statistics: pd.DataFrame | None = None,
    export_objects: bool = False,
    export_planes: bool = True,
    export_detections: bool = False,
) -> dict[str, Path]:
    """Write combined tables without modifying any individual run folder."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, Path] = {}

    def write(name: str, frame: pd.DataFrame) -> None:
        path = output_dir / name
        frame.to_csv(path, index=False)
        outputs[name] = path.resolve()

    write("combined_per_brain.csv", results.brains)
    write("group_summary.csv", summary)
    if statistics is not None:
        write("statistical_test_results.csv", statistics)
        if not statistics.empty and statistics["test_method"].nunique() == 1:
            method_name = str(statistics["test_method"].iloc[0])
            write(f"{method_name}_results.csv", statistics)
    write("result_manifest.csv", results.manifest)
    for name, frame in validation.items():
        write(f"validation_{name}.csv", frame)
    if export_objects:
        if results.objects.empty:
            raise ValueError("Objects were not loaded; set load_objects=True")
        write("combined_accepted_objects.csv", results.objects)
    if export_planes:
        if results.planes.empty:
            raise ValueError("Planes were not loaded; set load_planes=True")
        write("combined_per_plane.csv", results.planes)
    if export_detections:
        if results.detections.empty:
            raise ValueError("Detections were not loaded; set load_detections=True")
        write("combined_per_detection.csv", results.detections)
    return outputs


def _condition_colors(order: Sequence[str]) -> dict[str, Any]:
    cmap = plt.get_cmap("tab10")
    return {condition: cmap(index % 10) for index, condition in enumerate(order)}


def _statistics_method_caption(statistics: pd.DataFrame | None) -> str:
    if statistics is None or statistics.empty:
        return ""
    first = statistics.iloc[0]
    if first["test_method"] == "anova_tukey":
        return "One-way ANOVA + Tukey HSD"
    test_label = "Student t-tests" if first["pairwise_test"] == "student_ttest" else "Welch t-tests"
    adjustment = str(first["p_adjust_method"])
    if adjustment in {"none", "raw", "unadjusted"}:
        return f"Pairwise {test_label}, unadjusted p-values"
    return f"Pairwise {test_label} + {adjustment} correction"


def _finish_figure(
    fig: plt.Figure,
    output_path: str | Path,
    *,
    show: bool,
    dpi: int,
) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    if show:
        try:
            from IPython.display import display

            display(fig)
        except ImportError:
            pass
    plt.close(fig)
    return output_path.resolve()


def _group_box_scatter(
    axis: plt.Axes,
    frame: pd.DataFrame,
    metric: str,
    order: Sequence[str],
    colors: Mapping[str, Any],
    *,
    seed: int,
) -> None:
    positions = np.arange(len(order), dtype=float)
    grouped_values = [
        pd.to_numeric(
            frame.loc[frame["condition"].astype(str) == condition, metric],
            errors="coerce",
        ).dropna().to_numpy()
        for condition in order
    ]
    box = axis.boxplot(
        grouped_values,
        positions=positions,
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.4},
        whiskerprops={"linewidth": 1.0},
        capprops={"linewidth": 1.0},
    )
    for patch, condition in zip(box["boxes"], order):
        patch.set_facecolor(colors[condition])
        patch.set_alpha(0.28)
        patch.set_edgecolor(colors[condition])
    rng = np.random.default_rng(seed)
    for position, condition, values in zip(positions, order, grouped_values):
        if values.size == 0:
            continue
        jitter = rng.normal(0.0, 0.055, size=values.size)
        axis.scatter(
            np.full(values.size, position) + jitter,
            values,
            s=25,
            color=colors[condition],
            edgecolor="white",
            linewidth=0.35,
            alpha=0.85,
            zorder=3,
        )
    labels = [
        f"{condition}\n(n={len(values)})"
        for condition, values in zip(order, grouped_values)
    ]
    axis.set_xticks(positions, labels=labels, rotation=25, ha="right")
    axis.grid(axis="y", alpha=0.22)


def _annotate_posthoc_comparisons(
    axis: plt.Axes,
    statistics: pd.DataFrame,
    metric: str,
    order: Sequence[str],
    comparisons: Sequence[tuple[str, str]],
    *,
    show_non_significant: bool,
) -> None:
    """Add brackets and Tukey significance stars for selected comparisons."""
    if statistics.empty:
        return
    metric_rows = statistics.loc[statistics["metric"] == metric]
    position_by_group = {str(group): index for index, group in enumerate(order)}
    annotations: list[tuple[int, int, str]] = []
    for left, right in comparisons:
        left = str(left)
        right = str(right)
        if left not in position_by_group or right not in position_by_group:
            continue
        row = metric_rows.loc[
            ((metric_rows["group1"] == left) & (metric_rows["group2"] == right))
            | ((metric_rows["group1"] == right) & (metric_rows["group2"] == left))
        ]
        if row.empty:
            continue
        label = str(row.iloc[0]["significance"])
        if not label or (label == "ns" and not show_non_significant):
            continue
        x1 = position_by_group[left]
        x2 = position_by_group[right]
        annotations.append((min(x1, x2), max(x1, x2), label))
    if not annotations:
        return

    # Short comparisons are placed closest to the data; wider comparisons are
    # stacked above them. This keeps shared endpoints and nested brackets legible.
    annotations.sort(key=lambda item: (item[1] - item[0], item[0]))
    y_min, y_max = axis.get_ylim()
    if axis.get_yscale() == "log" and y_min > 0:
        ratio = max(y_max / y_min, 1.5)
        step_factor = ratio ** 0.095
        bracket_factor = ratio ** 0.025
        y = y_max * step_factor
        for x1, x2, label in annotations:
            y_top = y * bracket_factor
            axis.plot(
                [x1, x1, x2, x2],
                [y, y_top, y_top, y],
                color="black",
                linewidth=1.0,
                clip_on=False,
            )
            axis.text(
                (x1 + x2) / 2,
                y_top,
                label,
                ha="center",
                va="bottom",
                fontsize=11,
                clip_on=False,
            )
            y = y_top * step_factor
        axis.set_ylim(y_min, y * bracket_factor)
    else:
        span = max(y_max - y_min, max(abs(y_max), 1.0) * 0.1)
        step = 0.105 * span
        bracket_height = 0.025 * span
        y = y_max + 0.055 * span
        for x1, x2, label in annotations:
            y_top = y + bracket_height
            axis.plot(
                [x1, x1, x2, x2],
                [y, y_top, y_top, y],
                color="black",
                linewidth=1.0,
                clip_on=False,
            )
            axis.text(
                (x1 + x2) / 2,
                y_top,
                label,
                ha="center",
                va="bottom",
                fontsize=11,
                clip_on=False,
            )
            y += step
        axis.set_ylim(y_min, y + 0.045 * span)


def plot_brain_metrics(
    results: BatchResults,
    output_path: str | Path,
    metrics: Mapping[str, str] = DEFAULT_BRAIN_METRICS,
    *,
    log_metrics: Sequence[str] = (),
    figure_title: str = "Per-brain lipid-droplet summary",
    statistics: pd.DataFrame | None = None,
    statistical_test_method: str = "anova_tukey",
    ttest_equal_var: bool = False,
    ttest_p_adjust_method: str = "holm",
    posthoc_comparisons: Sequence[tuple[str, str]] = DEFAULT_POSTHOC_COMPARISONS,
    annotate_statistics: bool = True,
    show_non_significant: bool = True,
    show: bool = True,
    dpi: int = 180,
) -> Path:
    """Plot per-brain distributions; every dot is one biological replicate."""
    missing = [metric for metric in metrics if metric not in results.brains]
    if missing:
        raise KeyError(f"Brain metrics not found: {missing}")
    if annotate_statistics and statistics is None:
        statistics = calculate_statistical_tests(
            results.brains,
            metrics,
            method=statistical_test_method,
            group_order=results.group_order,
            selected_comparisons=posthoc_comparisons,
            ttest_equal_var=ttest_equal_var,
            ttest_p_adjust_method=ttest_p_adjust_method,
        )
    if not metrics:
        raise ValueError("Choose at least one brain metric")
    n_columns = min(3, len(metrics))
    n_rows = int(math.ceil(len(metrics) / n_columns))
    panel_width = max(5.2, 0.70 * len(results.group_order) + 2.0)
    fig, axes = plt.subplots(n_rows, n_columns, figsize=(panel_width * n_columns, 4.6 * n_rows))
    axes_array = np.atleast_1d(axes).ravel()
    colors = _condition_colors(results.group_order)
    for index, (metric, label) in enumerate(metrics.items()):
        axis = axes_array[index]
        _group_box_scatter(
            axis,
            results.brains,
            metric,
            results.group_order,
            colors,
            seed=1931 + index,
        )
        axis.set_title(label)
        axis.set_ylabel(label)
        if metric in set(log_metrics):
            axis.set_yscale("log")
        if annotate_statistics and statistics is not None:
            _annotate_posthoc_comparisons(
                axis,
                statistics,
                metric,
                results.group_order,
                posthoc_comparisons,
                show_non_significant=show_non_significant,
            )
    for axis in axes_array[len(metrics) :]:
        axis.axis("off")
    caption = _statistics_method_caption(statistics) if annotate_statistics else ""
    title = figure_title
    if caption:
        title = f"{title}\n{caption}"
    fig.suptitle(title, fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return _finish_figure(fig, output_path, show=show, dpi=dpi)


def plot_qc_relationships(
    results: BatchResults,
    output_path: str | Path,
    *,
    statistics: pd.DataFrame | None = None,
    statistical_test_method: str = "anova_tukey",
    ttest_equal_var: bool = False,
    ttest_p_adjust_method: str = "holm",
    posthoc_comparisons: Sequence[tuple[str, str]] = DEFAULT_POSTHOC_COMPARISONS,
    annotate_statistics: bool = True,
    show_non_significant: bool = True,
    show: bool = True,
    dpi: int = 180,
) -> Path:
    """Plot ROI and segmentation diagnostics that can reveal batch artifacts."""
    brains = results.brains
    colors = _condition_colors(results.group_order)
    if annotate_statistics and statistics is None:
        statistics = calculate_statistical_tests(
            brains,
            QC_BOX_METRICS,
            method=statistical_test_method,
            group_order=results.group_order,
            selected_comparisons=posthoc_comparisons,
            ttest_equal_var=ttest_equal_var,
            ttest_p_adjust_method=ttest_p_adjust_method,
        )
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    for index, (metric, label) in enumerate(QC_BOX_METRICS.items()):
        axis = axes.ravel()[index]
        _group_box_scatter(
            axis,
            brains,
            metric,
            results.group_order,
            colors,
            seed=2718 + index,
        )
        axis.set_title(label)
        axis.set_ylabel(label)
        if annotate_statistics and statistics is not None:
            _annotate_posthoc_comparisons(
                axis,
                statistics,
                metric,
                results.group_order,
                posthoc_comparisons,
                show_non_significant=show_non_significant,
            )
    axis = axes[1, 1]
    for condition in results.group_order:
        subset = brains.loc[brains["condition"].astype(str) == condition]
        axis.scatter(
            subset["roi_volume_1e6_um3"],
            subset["ld_count_density_per_1e5_um3"],
            s=35,
            color=colors[condition],
            edgecolor="white",
            linewidth=0.4,
            alpha=0.85,
            label=condition,
        )
    axis.set_xlabel("ROI volume (1e6 um3)")
    axis.set_ylabel("LD count density (per 1e5 um3)")
    axis.set_title("Density versus analyzed volume")
    axis.grid(alpha=0.22)
    axis.legend(fontsize=8)
    caption = _statistics_method_caption(statistics) if annotate_statistics else ""
    title = "Batch-level analysis diagnostics"
    if caption:
        title = f"{title}\n{caption}"
    fig.suptitle(title, fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _finish_figure(fig, output_path, show=show, dpi=dpi)


def plot_object_distributions(
    results: BatchResults,
    output_path: str | Path,
    metrics: Mapping[str, str] = DEFAULT_OBJECT_METRICS,
    *,
    histogram_bins: int = 35,
    upper_quantile: float = 0.995,
    show: bool = True,
    dpi: int = 180,
) -> Path:
    """Plot condition-normalized distributions of accepted LD objects."""
    if results.objects.empty:
        raise ValueError("Objects were not loaded; set load_objects=True")
    if not 0.9 <= upper_quantile <= 1.0:
        raise ValueError("upper_quantile must be in [0.9, 1.0]")
    missing = [metric for metric in metrics if metric not in results.objects]
    if missing:
        raise KeyError(f"Object metrics not found: {missing}")
    n_columns = 2
    n_rows = int(math.ceil(len(metrics) / n_columns))
    fig, axes = plt.subplots(n_rows, n_columns, figsize=(7.0 * n_columns, 4.8 * n_rows))
    axes_array = np.atleast_1d(axes).ravel()
    colors = _condition_colors(results.group_order)
    for index, (metric, label) in enumerate(metrics.items()):
        axis = axes_array[index]
        all_values = pd.to_numeric(results.objects[metric], errors="coerce").dropna()
        if all_values.empty:
            axis.set_visible(False)
            continue
        if metric == "n_planes":
            lower = int(all_values.min())
            upper = int(all_values.max())
            bins = np.arange(lower - 0.5, upper + 1.5, 1.0)
        else:
            lower = float(all_values.min())
            upper = float(all_values.quantile(upper_quantile))
            if upper <= lower:
                upper = float(all_values.max()) + 1e-6
            bins = np.linspace(lower, upper, histogram_bins + 1)
        centers = 0.5 * (bins[:-1] + bins[1:])
        for condition in results.group_order:
            values = pd.to_numeric(
                results.objects.loc[
                    results.objects["condition"].astype(str) == condition,
                    metric,
                ],
                errors="coerce",
            ).dropna()
            values = values[(values >= bins[0]) & (values <= bins[-1])]
            counts, _ = np.histogram(values, bins=bins)
            proportion = counts / counts.sum() if counts.sum() else counts.astype(float)
            axis.plot(
                centers,
                proportion,
                color=colors[condition],
                linewidth=1.8,
                label=f"{condition} (n={len(values):,})",
            )
        axis.set_xlabel(label)
        axis.set_ylabel("Proportion of accepted objects")
        axis.set_title(label)
        axis.grid(alpha=0.22)
        axis.legend(fontsize=8)
    for axis in axes_array[len(metrics) :]:
        axis.axis("off")
    fig.suptitle("Accepted-object distributions", fontsize=15)
    fig.tight_layout()
    return _finish_figure(fig, output_path, show=show, dpi=dpi)


def plot_depth_profiles(
    results: BatchResults,
    output_path: str | Path,
    *,
    bins: int = 51,
    show: bool = True,
    dpi: int = 180,
) -> Path:
    """Plot mean +/- SEM profiles after normalizing each brain's active z extent."""
    if results.planes.empty:
        raise ValueError("Planes were not loaded; set load_planes=True")
    if bins < 10:
        raise ValueError("bins must be at least 10")
    frame = results.planes.loc[
        (results.planes["roi_area_um2"] > 0)
        & results.planes["normalized_depth"].between(0.0, 1.0)
    ].copy()
    frame["depth_bin"] = np.rint(frame["normalized_depth"] * (bins - 1)).astype(int)
    metrics = {
        "roi_area_relative": "ROI area / brain maximum",
        "candidate_density_per_1e3_um2": "Candidates per 1e3 um2",
        "foreground_fraction": "Foreground area / ROI area",
    }
    colors = _condition_colors(results.group_order)
    fig, axes = plt.subplots(1, len(metrics), figsize=(6.0 * len(metrics), 5.0))
    for axis, (metric, label) in zip(np.atleast_1d(axes), metrics.items()):
        per_brain = (
            frame.groupby(
                ["condition", "result_path", "depth_bin"],
                observed=True,
                sort=False,
            )[metric]
            .mean()
            .reset_index()
        )
        grouped = (
            per_brain.groupby(["condition", "depth_bin"], observed=True, sort=False)[metric]
            .agg(["mean", "std", "count"])
            .reset_index()
        )
        grouped["sem"] = grouped["std"] / np.sqrt(grouped["count"].clip(lower=1))
        for condition in results.group_order:
            subset = grouped.loc[grouped["condition"].astype(str) == condition].sort_values(
                "depth_bin"
            )
            if subset.empty:
                continue
            x = subset["depth_bin"].to_numpy(float) / (bins - 1)
            mean = subset["mean"].to_numpy(float)
            sem = subset["sem"].fillna(0.0).to_numpy(float)
            axis.plot(x, mean, color=colors[condition], linewidth=2.0, label=condition)
            axis.fill_between(x, mean - sem, mean + sem, color=colors[condition], alpha=0.16)
        axis.set_xlabel("Normalized active z depth")
        axis.set_ylabel(label)
        axis.set_title(label)
        axis.grid(alpha=0.22)
    axes[0].legend(fontsize=8)
    fig.suptitle("Normalized z-depth profiles (mean +/- SEM across brains)", fontsize=15)
    fig.tight_layout()
    return _finish_figure(fig, output_path, show=show, dpi=dpi)


def generate_standard_report(
    results: BatchResults,
    output_dir: str | Path,
    *,
    brain_metrics: Mapping[str, str] = DEFAULT_BRAIN_METRICS,
    object_metrics: Mapping[str, str] = DEFAULT_OBJECT_METRICS,
    log_brain_metrics: Sequence[str] = (),
    export_objects: bool = False,
    export_planes: bool = True,
    export_detections: bool = False,
    depth_bins: int = 51,
    object_histogram_bins: int = 35,
    object_upper_quantile: float = 0.995,
    include_intensity_diagnostics: bool = True,
    include_size_diagnostics: bool = True,
    ld_size_metrics: Mapping[str, str] = DEFAULT_LD_SIZE_METRICS,
    size_scatter_max_per_brain: int = 200,
    log_intensity_metrics: Sequence[str] = (),
    statistical_test_method: str = "anova_tukey",
    ttest_equal_var: bool = False,
    ttest_p_adjust_method: str = "holm",
    posthoc_comparisons: Sequence[tuple[str, str]] = DEFAULT_POSTHOC_COMPARISONS,
    annotate_statistics: bool = True,
    show_non_significant: bool = True,
    show: bool = True,
    dpi: int = 180,
) -> BatchReport:
    """Generate summaries plus optional signal and brain-aware size diagnostics.

    Signal reconstruction reads saved CSV columns only, not the IMS or any
    newly edited ROI. Droplet-level size boxes are descriptive; all significance
    tests (including per-brain size summaries) use one brain as one replicate.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    diagnostic_metrics = {}
    if include_intensity_diagnostics or include_size_diagnostics:
        results = prepare_diagnostic_results(
            results, include_intensity=include_intensity_diagnostics,
            size_metrics=ld_size_metrics if include_size_diagnostics else {})
        if include_intensity_diagnostics:
            diagnostic_metrics.update(DEFAULT_INTENSITY_METRICS)
        if include_size_diagnostics:
            diagnostic_metrics.update(size_summary_metrics(ld_size_metrics))
    summary = summarize_brains(results.brains, {**brain_metrics, **diagnostic_metrics})
    statistical_metrics = {**brain_metrics, **diagnostic_metrics}
    for metric, label in QC_BOX_METRICS.items():
        statistical_metrics.setdefault(metric, label)
    statistics = calculate_statistical_tests(
        results.brains,
        statistical_metrics,
        method=statistical_test_method,
        group_order=results.group_order,
        selected_comparisons=posthoc_comparisons,
        ttest_equal_var=ttest_equal_var,
        ttest_p_adjust_method=ttest_p_adjust_method,
    )
    validation = validate_batch_results(results)
    exported = export_batch_tables(
        results,
        summary,
        validation,
        output_dir,
        statistics=statistics,
        export_objects=export_objects,
        export_planes=export_planes,
        export_detections=export_detections,
    )
    figures = {
        "brain_metrics": plot_brain_metrics(
            results,
            output_dir / "brain_metrics.png",
            brain_metrics,
            log_metrics=log_brain_metrics,
            statistics=statistics,
            statistical_test_method=statistical_test_method,
            ttest_equal_var=ttest_equal_var,
            ttest_p_adjust_method=ttest_p_adjust_method,
            posthoc_comparisons=posthoc_comparisons,
            annotate_statistics=annotate_statistics,
            show_non_significant=show_non_significant,
            show=show,
            dpi=dpi,
        ),
        "qc_relationships": plot_qc_relationships(
            results,
            output_dir / "qc_relationships.png",
            statistics=statistics,
            statistical_test_method=statistical_test_method,
            ttest_equal_var=ttest_equal_var,
            ttest_p_adjust_method=ttest_p_adjust_method,
            posthoc_comparisons=posthoc_comparisons,
            annotate_statistics=annotate_statistics,
            show_non_significant=show_non_significant,
            show=show,
            dpi=dpi,
        ),
        "depth_profiles": plot_depth_profiles(
            results,
            output_dir / "depth_profiles.png",
            bins=depth_bins,
            show=show,
            dpi=dpi,
        ),
    }
    if not results.objects.empty:
        figures["object_distributions"] = plot_object_distributions(
            results,
            output_dir / "object_distributions.png",
            object_metrics,
            histogram_bins=object_histogram_bins,
            upper_quantile=object_upper_quantile,
            show=show,
            dpi=dpi,
        )
    diagnostic_plot_options = dict(
        statistics=statistics, statistical_test_method=statistical_test_method,
        ttest_equal_var=ttest_equal_var, ttest_p_adjust_method=ttest_p_adjust_method,
        posthoc_comparisons=posthoc_comparisons, annotate_statistics=annotate_statistics,
        show_non_significant=show_non_significant, show=show, dpi=dpi)
    diagnostic_info = {
        "intensity_definition": "sum(core_component_integrated_snr * noise_sigma_at_native_z)",
        "support": "Accepted tracks only; shared (z_index, core_component_label) counted once",
        "background": "Existing I_smooth minus Gaussian background; not raw intensity or integrated SNR",
        "limitations": "High-threshold core, not full LD boundary or calibrated lipid mass. Requires matched acquisition settings; no bleaching correction.",
        "size_definition": "Selected LoG diameter is a discrete-scale estimate, not a precise segmented 3D diameter",
        "size_boxes": "All accepted objects; only scatter dots sampled independently per brain",
        "statistical_unit": "One brain; no tests on pooled individual LDs",
        "non_estimable_tests": "Comparisons with no estimable within-group variance retain NaN p-values and test_note; no stars are drawn",
        "include_intensity_diagnostics": include_intensity_diagnostics,
        "include_size_diagnostics": include_size_diagnostics,
        "size_scatter_max_per_brain": size_scatter_max_per_brain,
    }
    if include_intensity_diagnostics:
        available = results.brains["ld_core_integrated_corrected_intensity_au"].notna()
        diagnostic_info["intensity_available_brains"] = int(available.sum())
        diagnostic_info["intensity_unavailable_brains"] = int((~available).sum())
        if available.any():
            figures["ld_intensity"] = plot_brain_metrics(
                results, output_dir / "ld_intensity.png", DEFAULT_INTENSITY_METRICS,
                figure_title="Background-corrected LD core signal (each dot = one brain)",
                log_metrics=log_intensity_metrics, **diagnostic_plot_options)
        if (~available).any():
            print(f"Intensity unavailable for {(~available).sum()} brain(s); see ld_intensity_note in per_brain_diagnostics.csv")
    if include_size_diagnostics and not results.objects.empty and ld_size_metrics:
        figures["ld_size_box_scatter"] = plot_ld_size_box_scatter(
            results, output_dir / "ld_size_box_scatter.png", metrics=ld_size_metrics,
            max_points_per_brain=size_scatter_max_per_brain, show=show, dpi=dpi)
        figures["ld_size_per_brain"] = plot_brain_metrics(
            results, output_dir / "ld_size_per_brain.png", size_summary_metrics(ld_size_metrics),
            figure_title="LD size summaries (each dot = one brain; tests across brains)",
            **diagnostic_plot_options)
    elif include_size_diagnostics and ld_size_metrics:
        diagnostic_info["size_diagnostics_note"] = "No accepted object rows loaded; no size figures produced"
        if results.brains["ld_count"].sum() > 0:
            print("Size diagnostics unavailable: object rows are not loaded; set load_objects=True")
    if include_intensity_diagnostics or include_size_diagnostics:
        diagnostic_columns = list(diagnostic_metrics)
        diagnostic_columns.extend(c for c in results.brains if c.startswith("ld_intensity_") or c.endswith("_valid_n_per_brain"))
        table_path = output_dir / "per_brain_diagnostics.csv"
        results.brains[["brain_id", "result_path", "condition", "ld_count", *diagnostic_columns]].to_csv(table_path, index=False)
        exported["per_brain_diagnostics.csv"] = table_path.resolve()
        method_path = output_dir / "diagnostic_methods.json"
        method_path.write_text(json.dumps(diagnostic_info, ensure_ascii=False, indent=2), encoding="utf-8")
        exported["diagnostic_methods.json"] = method_path.resolve()
    return BatchReport(
        summary=summary,
        statistics=statistics,
        validation=validation,
        exported_tables=exported,
        figures=figures,
    )


def display_report_outputs(report: BatchReport) -> None:
    """Display the output-file manifest and selected post-hoc results in Jupyter."""
    from IPython.display import Markdown, display

    records = [
        {"type": "table", "name": name, "path": str(path)}
        for name, path in report.exported_tables.items()
    ]
    records.extend(
        {"type": "figure", "name": name, "path": str(path)}
        for name, path in report.figures.items()
    )
    display(Markdown("### Saved report outputs"))
    display(pd.DataFrame.from_records(records))
    selected = report.statistics.loc[report.statistics["display_on_plot"]]
    if not selected.empty:
        method = str(selected["test_method"].iloc[0])
        display(Markdown(f"### Selected group comparisons: `{method}`"))
        columns = ["label", "group1", "group2"]
        if selected["omnibus_p"].notna().any():
            columns.append("omnibus_p")
        if selected["p_raw"].notna().any():
            columns.append("p_raw")
        columns.extend(["p_adjusted", "p_adjust_method", "significance"])
        if "test_note" in selected and selected["test_note"].ne("").any():
            columns.append("test_note")
        display(selected[columns].reset_index(drop=True))


__all__ = [
    "BatchReport",
    "BatchResults",
    "DEFAULT_BRAIN_METRICS",
    "DEFAULT_OBJECT_METRICS",
    "DEFAULT_INTENSITY_METRICS",
    "DEFAULT_LD_SIZE_METRICS",
    "DEFAULT_POSTHOC_COMPARISONS",
    "QC_BOX_METRICS",
    "calculate_anova_tukey",
    "calculate_pairwise_ttests",
    "calculate_statistical_tests",
    "discover_result_runs",
    "display_batch_overview",
    "display_report_outputs",
    "export_batch_tables",
    "generate_standard_report",
    "load_batch_results",
    "plot_brain_metrics",
    "plot_depth_profiles",
    "plot_object_distributions",
    "plot_ld_size_box_scatter",
    "prepare_diagnostic_results",
    "plot_qc_relationships",
    "summarize_brains",
    "validate_batch_results",
]
