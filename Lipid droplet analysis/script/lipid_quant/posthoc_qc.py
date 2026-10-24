"""Generate visual QC from saved lipid-droplet detections without re-detection.

The accepted tracks in per_object.csv and their per_detection.csv coordinates
are the sole source of marker positions. Raw IMS planes are read only for visual
overlays; no ROI construction, preprocessing, LoG, or z linking is rerun.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from imaris_ims_file_reader.ims import ims

from .paths import PathInput, normalize_paths
from .pipeline import _choose_qc_slices, _save_3d_detection_marker_tiff, _save_depth_qc


QC_FILES = ("projection_accepted.png", "slice_overlays_accepted.png", "depth_profiles.png")
MARKER_FILE = "detections_3d_mask.tif"
MATCH_FIELDS = {"brain_id", "input_name", "input_path", "result_name", "result_path"}


def _result_directories(inputs: PathInput, *, recursive: bool) -> list[Path]:
    found: dict[str, Path] = {}
    for root in normalize_paths(inputs):
        if not root.is_dir():
            raise ValueError(f"Expected a batch or result directory: {root}")
        if (root / "per_brain.csv").is_file():
            matches = [root]
        else:
            csv_files = root.rglob("per_brain.csv") if recursive else root.glob("*/per_brain.csv")
            matches = [path.parent for path in csv_files]
        for directory in matches:
            resolved = directory.resolve()
            found.setdefault(str(resolved).casefold(), resolved)
    if not found:
        raise FileNotFoundError("No result directories containing per_brain.csv in the selected folders")
    return sorted(found.values(), key=lambda path: str(path).casefold())


def _one_brain_row(result_dir: Path) -> pd.Series:
    table = pd.read_csv(result_dir / "per_brain.csv")
    if len(table) != 1:
        raise ValueError(f"Expected exactly one row in {result_dir / 'per_brain.csv'}")
    required = {"brain_id", "file", "ld_count", "shape_z", "shape_y", "shape_x"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Missing per_brain columns in {result_dir}: {sorted(missing)}")
    return table.iloc[0]


def plan_posthoc_qc(
    inputs: PathInput,
    *,
    regex: str | None = None,
    match_field: str = "brain_id",
    recursive: bool = False,
) -> pd.DataFrame:
    """Select saved result directories by folder scope and regex search.

    ``recursive=False`` includes only immediate result children of each batch
    root; a directly supplied result directory is always included. This avoids
    accidentally including nested ``removed`` or earlier summary folders.
    """
    if match_field not in MATCH_FIELDS:
        raise ValueError(f"match_field must be one of {sorted(MATCH_FIELDS)}")
    pattern = re.compile(regex) if regex else None
    records: list[dict[str, Any]] = []
    for directory in _result_directories(inputs, recursive=recursive):
        brain = _one_brain_row(directory)
        info_path = directory / "run_info.json"
        if not info_path.is_file():
            raise FileNotFoundError(f"Missing run_info.json in {directory}")
        info = json.loads(info_path.read_text(encoding="utf-8"))
        source_value = info.get("input_file")
        if not source_value:
            raise ValueError(f"Missing input_file in {info_path}")
        source = Path(source_value)
        fields = {
            "brain_id": str(brain["brain_id"]),
            "input_name": str(brain["file"]),
            "input_path": str(source),
            "result_name": directory.name,
            "result_path": str(directory),
        }
        if pattern and not pattern.search(fields[match_field]):
            continue
        qc_dir = directory / "qc_posthoc"
        records.append({
            "result_dir": str(directory),
            "brain_id": str(brain["brain_id"]),
            "input_file": str(source),
            "source_exists": source.is_file(),
            "ld_count": int(brain["ld_count"]),
            "qc_dir": str(qc_dir),
            "qc_complete": all((qc_dir / name).is_file() for name in QC_FILES),
            "marker_exists": (qc_dir / MARKER_FILE).is_file(),
        })
    if not records:
        raise ValueError(f"No results matched regex={regex!r} on {match_field}")
    return pd.DataFrame.from_records(records)


def _read_saved_detections(result_dir: Path) -> tuple[pd.Series, dict[str, Any], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    brain = _one_brain_row(result_dir)
    info = json.loads((result_dir / "run_info.json").read_text(encoding="utf-8"))
    config_path = result_dir / "config_used.yaml"
    if not config_path.is_file():
        raise FileNotFoundError(config_path)
    import yaml

    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    detections = pd.read_csv(result_dir / "per_detection.csv")
    objects = pd.read_csv(result_dir / "per_object.csv")
    planes = pd.read_csv(result_dir / "per_plane.csv")
    for table, required, name in (
        (detections, {"track_id", "z_index", "y_px", "x_px"}, "per_detection"),
        (objects, {"object_id", "accepted", "single_plane"}, "per_object"),
        (planes, {"z_index", "z_um", "candidate_count", "foreground_area_um2", "roi_area_um2", "raw_median", "raw_p99", "noise_sigma"}, "per_plane"),
    ):
        missing = required - set(table.columns)
        if missing:
            raise ValueError(f"Missing {name} columns in {result_dir}: {sorted(missing)}")
    if not detections.empty and "log_response" not in detections:
        raise ValueError(f"Missing log_response in {result_dir / 'per_detection.csv'}")
    if detections.empty:
        detections["log_response"] = pd.Series(dtype=float)
    accepted_text = objects["accepted"].astype(str).str.lower()
    if not accepted_text.isin(["true", "false"]).all():
        raise ValueError(f"Invalid accepted values in {result_dir / 'per_object.csv'}")
    objects = objects.copy()
    objects["accepted"] = accepted_text.eq("true")
    if objects["object_id"].duplicated().any():
        raise ValueError(f"Duplicate object_id in {result_dir / 'per_object.csv'}")
    accepted = objects.loc[objects["accepted"]]
    if len(accepted) != int(brain["ld_count"]):
        raise ValueError(f"Accepted-object count differs from per_brain.ld_count in {result_dir}")
    shape = tuple(int(brain[key]) for key in ("shape_z", "shape_y", "shape_x"))
    plane_z = pd.to_numeric(planes["z_index"], errors="coerce").to_numpy(dtype=float)
    if (not np.isfinite(plane_z).all() or
            (plane_z != np.floor(plane_z)).any() or
            sorted(plane_z.astype(int).tolist()) != list(range(shape[0]))):
        raise ValueError(f"per_plane z coverage differs from saved volume shape in {result_dir}")
    accepted_ids = set(accepted["object_id"].astype(int))
    accepted_detections = detections[detections["track_id"].isin(accepted_ids)].copy()
    if accepted_ids - set(accepted_detections["track_id"].astype(int)):
        raise ValueError(f"Accepted object without saved detection in {result_dir}")
    if not accepted_detections.empty:
        for column, upper in zip(("z_index", "y_px", "x_px"), shape):
            values = pd.to_numeric(accepted_detections[column], errors="coerce").to_numpy(dtype=float)
            if not np.isfinite(values).all() or (values < 0).any() or (values >= upper).any() or (values != np.floor(values)).any():
                raise ValueError(f"Saved accepted {column} is outside image bounds in {result_dir}")
    return brain, {"run_info": info, "config": config}, detections, objects, planes


def _raw_mip_and_slices(source: Path, brain: pd.Series, config: dict[str, Any], selected_z: set[int]) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    input_cfg = config["input"]
    reader = ims(str(source), ResolutionLevelLock=int(input_cfg["resolution_level"]), squeeze_output=True)
    try:
        shape = tuple(int(brain[key]) for key in ("shape_z", "shape_y", "shape_x"))
        if tuple(int(value) for value in reader.shape[-3:]) != shape:
            raise ValueError(f"IMS shape changed since analysis: {source}")
        if source.name != str(brain["file"]):
            raise ValueError(f"IMS filename differs from per_brain.file: {source}")
        t, c = int(input_cfg["time_index"]), int(input_cfg["channel_index"])
        mip = np.full(shape[1:], -np.inf, dtype=np.float32)
        slices: dict[int, np.ndarray] = {}
        for z in range(shape[0]):
            raw = np.asarray(reader[t, c, z, :, :], dtype=np.float32)
            mip = np.maximum(mip, raw)
            if z in selected_z:
                slices[z] = raw.copy()
        return mip, slices
    finally:
        reader.close()


def _display_limits(raw: np.ndarray) -> tuple[float, float]:
    values = raw[np.isfinite(raw)]
    if values.size == 0:
        return 0.0, 1.0
    lo, hi = np.percentile(values, [1.0, 99.8])
    return float(lo), float(max(hi, lo + 1.0))


def _save_accepted_projection(path: Path, mip: np.ndarray, detections: pd.DataFrame, objects: pd.DataFrame, title: str) -> None:
    accepted = objects[objects["accepted"]]
    points = detections[detections["track_id"].isin(accepted["object_id"])]
    if not points.empty:
        # The brightest LoG detection is a real saved XY pixel coordinate; no
        # new center or object classification is inferred from the raw image.
        points = points.sort_values("log_response", ascending=False).drop_duplicates("track_id")
        points = points.merge(accepted[["object_id", "single_plane"]], left_on="track_id", right_on="object_id", validate="one_to_one")
    lo, hi = _display_limits(mip)
    fig, ax = plt.subplots(figsize=(11, 11), constrained_layout=True)
    ax.imshow(mip, cmap="gray", vmin=lo, vmax=hi)
    if not points.empty:
        single = points["single_plane"].astype(str).str.lower().eq("true")
        for flag, color, label in ((False, "lime", "multi-z"), (True, "orange", "single-z")):
            subset = points[single == flag]
            ax.scatter(subset["x_px"], subset["y_px"], s=7, facecolors="none", edgecolors=color, linewidths=0.5, label=f"{label} ({len(subset)})")
        ax.legend(loc="lower right", fontsize=8)
    ax.set_title(f"{title}\nAccepted objects on raw maximum projection (n={len(accepted)})")
    ax.axis("off")
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _save_accepted_slices(path: Path, slices: dict[int, np.ndarray], detections: pd.DataFrame, objects: pd.DataFrame, title: str) -> None:
    selected = sorted(slices)
    if not selected:
        raise ValueError("No raw slices selected for QC")
    accepted_ids = set(objects.loc[objects["accepted"], "object_id"].astype(int))
    accepted_detections = detections[detections["track_id"].isin(accepted_ids)]
    fig, axes = plt.subplots(len(selected), 1, figsize=(11, 7 * len(selected)), constrained_layout=True)
    axes = np.atleast_1d(axes)
    for ax, z in zip(axes, selected):
        raw = slices[z]
        lo, hi = _display_limits(raw)
        ax.imshow(raw, cmap="gray", vmin=lo, vmax=hi)
        points = accepted_detections[accepted_detections["z_index"] == z]
        if not points.empty:
            ax.scatter(points["x_px"], points["y_px"], s=16, facecolors="none", edgecolors="lime", linewidths=0.7)
        ax.set_title(f"z={z}: accepted detection centers (n={len(points)})")
        ax.axis("off")
    fig.suptitle(title)
    fig.savefig(path, dpi=130)
    plt.close(fig)


def _supplement_one(result_dir: Path, *, save_qc: bool, save_marker: bool, qc_slices: int, overwrite: bool) -> str:
    brain, metadata, detections, objects, planes = _read_saved_detections(result_dir)
    qc_dir = result_dir / "qc_posthoc"
    requested = ([qc_dir / name for name in QC_FILES] if save_qc else []) + ([qc_dir / MARKER_FILE] if save_marker else [])
    missing = [path for path in requested if overwrite or not path.is_file()]
    if not missing:
        return "already_exists"
    source = Path(metadata["run_info"]["input_file"])
    raw_qc_missing = save_qc and any(path.name in QC_FILES[:2] for path in missing)
    if raw_qc_missing and not source.is_file():
        raise FileNotFoundError(f"Original IMS is required for raw-image QC: {source}")
    qc_dir.mkdir(exist_ok=True)
    if save_marker and (overwrite or not (qc_dir / MARKER_FILE).is_file()):
        shape = tuple(int(brain[key]) for key in ("shape_z", "shape_y", "shape_x"))
        diameter = int(metadata["config"].get("output", {}).get("marker_3d_diameter_px", 5))
        _save_3d_detection_marker_tiff(qc_dir / MARKER_FILE, shape, detections, objects, diameter)
    if save_qc:
        active = planes.loc[planes["roi_area_um2"] > 0, "z_index"].to_numpy(int)
        if not len(active):
            active = planes["z_index"].to_numpy(int)
        selected = _choose_qc_slices(active, qc_slices)
        if not selected:
            raise ValueError("qc_slices must be at least 1 when save_qc=True")
        raw_needed = overwrite or any(not (qc_dir / name).is_file() for name in QC_FILES[:2])
        if raw_needed:
            mip, slices = _raw_mip_and_slices(source, brain, metadata["config"], selected)
            if overwrite or not (qc_dir / QC_FILES[0]).is_file():
                _save_accepted_projection(qc_dir / QC_FILES[0], mip, detections, objects, str(brain["file"]))
            if overwrite or not (qc_dir / QC_FILES[1]).is_file():
                _save_accepted_slices(qc_dir / QC_FILES[1], slices, detections, objects, str(brain["file"]))
        if overwrite or not (qc_dir / QC_FILES[2]).is_file():
            _save_depth_qc(qc_dir / QC_FILES[2], planes, str(brain["file"]))
    provenance = {
        "method": "saved_detections_only_no_redetection",
        "result_dir": str(result_dir),
        "input_file": str(source),
        "config_hash": str(brain.get("config_hash", "")),
        "ld_count": int(brain["ld_count"]),
        "save_qc": save_qc,
        "save_3d_marker_tiff": save_marker,
        "qc_note": "Raw IMS is used only for display; markers and object identity come from saved CSVs. No SNR/foreground mask is reconstructed.",
    }
    (qc_dir / "posthoc_info.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")
    return "completed"


def run_posthoc_qc(
    inputs: PathInput,
    *,
    regex: str | None = None,
    match_field: str = "brain_id",
    recursive: bool = False,
    save_qc: bool = True,
    save_3d_marker_tiff: bool = True,
    qc_slices: int = 6,
    overwrite: bool = False,
    run: bool = False,
    raise_on_error: bool = True,
) -> pd.DataFrame:
    """Plan or generate post-hoc QC/markers in each selected result's qc_posthoc.

    Folder scope is applied first, then ``regex`` via ``re.search`` to the
    selected ``match_field``. ``run=False`` is a read-only plan. Existing
    outputs are left untouched unless ``overwrite=True``.
    """
    if not save_qc and not save_3d_marker_tiff:
        raise ValueError("Enable save_qc and/or save_3d_marker_tiff")
    if save_qc and (isinstance(qc_slices, bool) or not isinstance(qc_slices, int) or qc_slices < 1):
        raise ValueError("qc_slices must be a positive integer")
    plan = plan_posthoc_qc(inputs, regex=regex, match_field=match_field, recursive=recursive)
    plan["status"] = "planned"
    plan["error"] = ""
    if not run:
        return plan
    for index, row in plan.iterrows():
        try:
            plan.at[index, "status"] = _supplement_one(Path(row["result_dir"]), save_qc=save_qc,
                                                        save_marker=save_3d_marker_tiff,
                                                        qc_slices=qc_slices, overwrite=overwrite)
        except Exception as exc:
            plan.at[index, "status"] = "failed"
            plan.at[index, "error"] = f"{type(exc).__name__}: {exc}"
            if raise_on_error:
                raise RuntimeError(f"Post-hoc QC failed for {row['brain_id']}: {exc}") from exc
    return plan
