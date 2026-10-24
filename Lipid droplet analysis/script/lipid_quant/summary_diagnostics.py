"""Read-only signal reconstruction and brain-aware LD size diagnostics.

No IMS is reopened and no detection/segmentation is repeated. Intensity refers
to the existing high-threshold core of accepted tracks, not the broad mask.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Mapping

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


DEFAULT_INTENSITY_METRICS = {
    "ld_core_integrated_corrected_intensity_au": "Background-corrected LD core\nintegrated signal (a.u.)",
    "ld_core_signal_density_per_1e5_um3": "Corrected LD core signal\nper 1e5 um3 ROI",
}
DEFAULT_LD_SIZE_METRICS = {"diameter_xy_um": "Selected LoG diameter (um)"}


def size_summary_metrics(metrics: Mapping[str, str]) -> dict[str, str]:
    return {f"{metric}_{stat}_per_brain": f"Per-brain {stat}\n{label}"
            for metric, label in metrics.items() for stat in ("mean", "median", "p90")}


def _accepted_objects(objects: pd.DataFrame, columns=None) -> pd.DataFrame:
    frame = objects
    if columns is not None:
        # The pooled table can contain hundreds of thousands of LDs. Do not
        # copy all coordinates/shape fields just to summarize one size metric.
        wanted = list(dict.fromkeys([*columns, "accepted"]))
        frame = frame.loc[:, [c for c in wanted if c in frame]]
    if frame.empty or "accepted" not in frame:
        return frame.copy()
    return frame.loc[frame["accepted"].astype(str).str.lower().isin(["true", "1", "1.0"])].copy()


def reconstruct_core_intensity(brain: Mapping) -> dict:
    """Recover sum(I_smooth - background) = sum(core SNR * plane noise).

    Shared components are counted once per z, matching ld_positive_volume_um3.
    Negative/absent/inconsistent inputs stay unavailable, never become zeros.
    """
    output = {"ld_core_integrated_corrected_intensity_au": np.nan,
              "ld_intensity_component_count": np.nan,
              "ld_intensity_status": "unavailable", "ld_intensity_note": "",
              "ld_intensity_method": "unique_accepted_core_snr_times_plane_noise"}
    if brain.get("ld_count") == 0:
        output.update(ld_core_integrated_corrected_intensity_au=0.0,
                      ld_intensity_component_count=0, ld_intensity_status="no_accepted_objects")
        return output
    run = Path(brain["result_path"])
    required = {"track_id", "z_index", "core_component_label", "core_component_integrated_snr"}
    try:
        objects = pd.read_csv(run / "per_object.csv", usecols=lambda c: c in {"object_id", "accepted"})
        detections = pd.read_csv(run / "per_detection.csv", usecols=lambda c: c in required)
        planes = pd.read_csv(run / "per_plane.csv", usecols=lambda c: c in {"z_index", "noise_sigma"})
        if not {"object_id", "accepted"} <= set(objects) or not required <= set(detections) or not {"z_index", "noise_sigma"} <= set(planes):
            raise ValueError("Saved CSV lacks accepted flags, core integral, z index, or plane noise")
        accepted = _accepted_objects(objects)
        ids = pd.to_numeric(accepted["object_id"], errors="coerce")
        if ids.isna().any() or ids.duplicated().any() or len(ids) != brain["ld_count"]:
            raise ValueError("Accepted-object IDs/count disagree with per_brain.csv")
        for column in required:
            detections[column] = pd.to_numeric(detections[column], errors="coerce")
        relevant = detections.loc[detections["track_id"].isin(ids)].copy()
        if not set(ids) <= set(relevant["track_id"]):
            raise ValueError("Some accepted tracks are absent from per_detection.csv")
        if not np.isfinite(relevant[list(required)].to_numpy(float)).all():
            raise ValueError("Non-finite accepted core measurements")
        relevant = relevant.loc[relevant["core_component_label"] > 0]
        if (relevant["core_component_integrated_snr"] < 0).any():
            raise ValueError("Negative integrated SNR in a positive LD core")
        if not relevant.empty and (relevant.groupby(["z_index", "core_component_label"])["core_component_integrated_snr"].nunique() > 1).any():
            raise ValueError("Shared-core records disagree about their integrated SNR")
        cores = relevant.drop_duplicates(["z_index", "core_component_label"])
        reference = brain.get("ld_core_integrated_snr_voxels", np.nan)
        if pd.notna(reference) and not np.isclose(cores["core_component_integrated_snr"].sum(), reference, rtol=1e-6, atol=1e-5):
            raise ValueError("Unique accepted core integral disagrees with per_brain.csv")
        planes = planes.apply(pd.to_numeric, errors="coerce")
        if planes["z_index"].duplicated().any():
            raise ValueError("Duplicate z indices in per_plane.csv")
        noise = cores["z_index"].map(planes.set_index("z_index")["noise_sigma"])
        if not np.isfinite(noise.to_numpy(float)).all() or (noise <= 0).any():
            raise ValueError("Missing/non-positive noise for an accepted LD core plane")
        total = float(np.sum(cores["core_component_integrated_snr"].to_numpy(float) * noise.to_numpy(float), dtype=np.float64))
        if not np.isfinite(total):
            raise ValueError("Non-finite reconstructed LD core signal")
        output.update(ld_core_integrated_corrected_intensity_au=total,
                      ld_intensity_component_count=len(cores), ld_intensity_status="reconstructed")
    except (OSError, ValueError) as error:
        output["ld_intensity_note"] = str(error)
    return output


def prepare_diagnostic_results(results, *, include_intensity=True,
                               size_metrics: Mapping[str, str] = DEFAULT_LD_SIZE_METRICS):
    """Return a new BatchResults; leave source CSVs and the input object intact.

    Process one brain's selected CSV columns at a time, so the full detection
    table need not be retained in memory (load_detections may remain False).
    """
    brains = results.brains.copy().reset_index(drop=True)
    if include_intensity:
        signal = pd.DataFrame([reconstruct_core_intensity(row) for row in brains.to_dict("records")])
        for column in signal:
            brains[column] = signal[column].to_numpy()
        volume = pd.to_numeric(brains["roi_volume_um3"], errors="coerce")
        brains["ld_core_signal_density_per_1e5_um3"] = (
            brains["ld_core_integrated_corrected_intensity_au"] / volume.where(np.isfinite(volume) & (volume > 0)) * 1e5)
    objects = _accepted_objects(results.objects, ["result_path", *size_metrics])
    for metric in size_metrics:
        valid = objects.copy() if metric in objects else pd.DataFrame()
        if not valid.empty:
            valid[metric] = pd.to_numeric(valid[metric], errors="coerce")
            valid = valid.loc[np.isfinite(valid[metric]) & (valid[metric] > 0)]
        grouped = valid.groupby("result_path")[metric] if not valid.empty else None
        for stat in ("mean", "median", "p90"):
            values = (grouped.quantile(0.90) if stat == "p90" else getattr(grouped, stat)()) if grouped is not None else pd.Series(dtype=float)
            brains[f"{metric}_{stat}_per_brain"] = brains["result_path"].map(values)
        counts = grouped.size() if grouped is not None else pd.Series(dtype=int)
        brains[f"{metric}_valid_n_per_brain"] = brains["result_path"].map(counts).fillna(0).astype(int)
    return replace(results, brains=brains)


def plot_ld_size_box_scatter(results, output_path, *, metrics=DEFAULT_LD_SIZE_METRICS,
                             max_points_per_brain=200, seed=1931, show=True, dpi=180):
    """Boxes use ALL LDs; only drawn dots are sampled, independently per brain.

    This pooled-object view is descriptive. Black diamonds are brain medians;
    inferential comparisons belong in the separate per-brain summary plot.
    """
    from .batch_summary import _condition_colors, _finish_figure
    if isinstance(max_points_per_brain, bool) or int(max_points_per_brain) != max_points_per_brain or max_points_per_brain < 1:
        raise ValueError("max_points_per_brain must be a positive integer")
    if not metrics:
        raise ValueError("Choose at least one LD size metric")
    objects = _accepted_objects(results.objects, ["result_path", "condition", *metrics])
    if objects.empty:
        raise ValueError("No accepted objects loaded; set load_objects=True")
    missing = set(metrics) - set(objects)
    if missing:
        raise KeyError(f"LD size metrics not found: {sorted(missing)}")
    colors = _condition_colors(results.group_order)
    width = max(8.5, 1.35 * len(results.group_order))
    fig, axes = plt.subplots(len(metrics), 1, figsize=(width, 5.0 * len(metrics)), squeeze=False)
    rng = np.random.default_rng(seed)
    for axis, (metric, label) in zip(axes.ravel(), metrics.items()):
        frame = objects.copy()
        frame[metric] = pd.to_numeric(frame[metric], errors="coerce")
        frame = frame.loc[np.isfinite(frame[metric]) & (frame[metric] > 0)]
        groups = [frame.loc[frame["condition"].astype(str) == condition] for condition in results.group_order]
        boxes = axis.boxplot([g[metric].to_numpy(float) for g in groups], positions=np.arange(len(groups)),
                             widths=0.50, patch_artist=True, showfliers=False,
                             medianprops={"color": "black", "linewidth": 1.5})
        labels = []
        for x, (condition, group, box) in enumerate(zip(results.group_order, groups, boxes["boxes"])):
            box.set_facecolor(colors[condition]); box.set_alpha(0.25)
            box.set_edgecolor(colors[condition])
            for _, brain in group.groupby("result_path", sort=True):
                values = brain[metric].to_numpy(float)
                dots = rng.choice(values, int(max_points_per_brain), replace=False) if len(values) > max_points_per_brain else values
                axis.scatter(x + rng.uniform(-0.20, 0.20, len(dots)), dots, s=7,
                             color=colors[condition], alpha=0.20, linewidths=0, rasterized=True)
            medians = group.groupby("result_path")[metric].median().to_numpy(float)
            axis.scatter(x + rng.uniform(-0.10, 0.10, len(medians)), medians, s=23,
                         color="#222222", marker="D", edgecolors="white", linewidths=0.4, zorder=4)
            labels.append(condition.replace(" | ", "\n") + f"\nn={group.result_path.nunique()} brains\nN={len(group):,} LDs")
        axis.set_xticks(np.arange(len(groups)), labels, fontsize=9)
        axis.set_ylabel(label)
        axis.set_title(label)
        axis.grid(axis="y", alpha=0.18)
        if not frame.empty:
            # Retain the entire range even though boxplot fliers are hidden.
            low, high = float(frame[metric].min()), float(frame[metric].max())
            pad = max((high-low)*0.06, abs(high)*0.03, 0.03)
            axis.set_ylim(max(0, low-pad), high+pad)
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle("LD size by group (descriptive; no droplet-level significance tests)", fontsize=13)
    fig.legend(handles=[Line2D([], [], marker="o", ls="", color="#777777", markersize=4,
                               label=f"LD dots: at most {max_points_per_brain}/brain; boxes use all LDs"),
                        Line2D([], [], marker="D", ls="", color="#222222", markersize=4, label="Each brain's median")],
               loc="lower center", fontsize=9, ncol=2)
    fig.tight_layout(rect=(0, 0.065, 1, 0.95))
    return _finish_figure(fig, output_path, show=show, dpi=dpi)


__all__ = ["DEFAULT_INTENSITY_METRICS", "DEFAULT_LD_SIZE_METRICS", "size_summary_metrics",
           "reconstruct_core_intensity", "prepare_diagnostic_results", "plot_ld_size_box_scatter"]
