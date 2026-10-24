"""Fixed classical image-processing pipeline for lipid-droplet quantification.

The pipeline deliberately separates two endpoints:

* LD-positive burden is measured from a per-plane hysteresis mask.
* LD count is measured from 2-D multi-scale LoG peaks linked across adjacent z planes.

No neural-network model is used.  All spatial parameters are expressed in micrometres.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
import warnings
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib

# Keep Jupyter's interactive/inline backend intact.  The command-line pipeline
# still uses a headless backend so batch jobs can save QC figures without a GUI.
if "ipykernel" not in sys.modules and "matplotlib.pyplot" not in sys.modules:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py
import numpy as np
import pandas as pd
import scipy.ndimage as ndi
import yaml
from imaris_ims_file_reader.ims import ims
from scipy.optimize import linear_sum_assignment
from skimage.draw import circle_perimeter
from skimage.feature import peak_local_max
from skimage.filters import apply_hysteresis_threshold, threshold_otsu
from skimage.measure import label
from skimage.morphology import (
    binary_closing,
    binary_dilation,
    binary_erosion,
    disk,
    remove_small_holes,
    remove_small_objects,
)
from skimage.segmentation import find_boundaries
from skimage.transform import resize

from .anatomical_roi import ANATOMICAL_DEFAULTS, build_anatomical_roi, resolve_sample_guide


DEFAULT_CONFIG: dict[str, Any] = {
    "input": {
        "time_index": 0,
        "channel_index": 0,
        "resolution_level": 0,
    },
    "roi": {
        **ANATOMICAL_DEFAULTS,
        "mode": "auto_central_brain",
        "mask_path": None,
        "auto_level": 2,
        "smooth_sigma_um": 2.0,
        "contrast_high_quantile": 0.90,
        "minimum_contrast_snr": 4.0,
        "minimum_component_area_um2": 300.0,
        "closing_radius_um": 2.0,
        "erosion_margin_um": 0.0,
        "tissue_expansion_um": 2.0,
        "tissue_support_snr": 4.0,
        "optic_background_sigma_um": 8.0,
        "optic_threshold_offset_sigma": 0.0,
        "optic_seed_quantile": 0.95,
        "optic_support_fraction": 0.35,
        "optic_minimum_area_um2": 200.0,
        "optic_exclusion_mode": "z_local_central_priority",
        "optic_minimum_z_planes": 2,
        "optic_z_context_um": 6.0,
        "optic_local_closing_um": 2.0,
        "optic_local_exclusion_margin_um": 2.0,
        "optic_footprint_closing_um": 4.0,
        "optic_image_border_distance_um": 3.0,
        "optic_exclusion_margin_um": 4.0,
        "central_closing_radius_um": 3.0,
        "central_hole_maximum_area_um2": 500.0,
        "central_expansion_um": 1.2,
        "central_z_max_expansion_um": 8.0,
        "terminal_area_fraction": 0.30,
        "central_terminal_rescue_um": 12.0,
        "terminal_optic_closing_um": 0.0,
        "terminal_optic_exclusion_margin_um": 0.0,
        "terminal_maximum_area_growth": 1.5,
        "terminal_keep_multiple_components": True,
        "terminal_tissue_support_snr": 4.0,
        "terminal_minimum_component_area_um2": 100.0,
        "terminal_tissue_closing_um": 1.0,
        "retain_terminal_tissue_components": True,
        "terminal_optic_seed_quantile": 0.95,
        "terminal_optic_support_fraction": 0.0,
        "terminal_exclude_border_connected_tissue": True,
        "terminal_short_gap_maximum_um": 6.0,
        "terminal_short_gap_area_fraction": 0.25,
        "crop_yx": None,
    },
    "preprocess": {
        "smooth_sigma_xy_um": 0.18,
        "background_sigma_xy_um": 3.0,
        "noise_sample_stride": 4,
        "snr_clip_for_detection": 40.0,
    },
    "segmentation": {
        "snr_low": 2.0,
        "snr_high": 3.5,
        "minimum_area_um2": 0.07,
    },
    "detection": {
        "processing_downsample": 1,
        "minimum_diameter_um": 0.9,
        "maximum_diameter_um": 3.0,
        "number_of_scales": 6,
        "log_response_threshold": 1.1,
        "minimum_peak_snr": 3.8,
        "minimum_hessian_roundness": 0.38,
        "minimum_radial_contrast_snr": 1.2,
        "maximum_core_component_area_um2": 10.0,
        "minimum_center_distance_um": 0.45,
        "link_distance_xy_um": 0.35,
        "maximum_scale_ratio": 2.0,
        "single_plane_minimum_snr": 5.5,
        "maximum_planes_per_object": 5,
    },
    "output": {
        "save_mask_tiff": False,
        "save_roi_mask_tiff": False,
        "mask_compression": "deflate",
        "save_qc": True,
        "qc_slices": 6,
        "qc_downsample": 2,
        "full_resolution_marker_radius_px": 5,
        "save_3d_marker_tiff": True,
        "marker_3d_diameter_px": 5,
        "print_plane_progress": True,
    },
    "batch": {
        "workers": 2,
        "save_qc": False,
        "save_3d_marker_tiff": False,
        "save_mask_tiff": False,
        "save_roi_mask_tiff": False,
        "print_plane_progress": False,
    },
}


@dataclass
class PlaneResult:
    raw: np.ndarray
    background: np.ndarray
    residual: np.ndarray
    snr: np.ndarray
    roi: np.ndarray
    foreground: np.ndarray
    candidates: pd.DataFrame
    noise_sigma: float


def _deep_merge(base: dict[str, Any], update: dict[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in base.items():
        result[key] = _deep_merge(value, {}) if isinstance(value, dict) else value
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    """Load a YAML configuration and fill omitted values with stable defaults."""
    user: dict[str, Any] = {}
    if path is not None:
        with Path(path).open("r", encoding="utf-8") as handle:
            user = yaml.safe_load(handle) or {}
    old_roi = user.get("roi", {})
    if (old_roi.get("mode", "auto_central_brain") == "auto_central_brain"
            and "central_method" not in old_roi and "optic_exclusion_mode" in old_roi):
        # Do not silently reinterpret a previously saved intensity-based YAML.
        old_roi["central_method"] = "intensity_legacy"
        warnings.warn("Legacy ROI YAML detected: using intensity_legacy. Set roi.central_method='anatomical' explicitly to migrate.", UserWarning, stacklevel=2)
    from .anatomical_v2 import resolve_v2_config
    return resolve_v2_config(_deep_merge(DEFAULT_CONFIG, user))


def _config_hash(config: dict[str, Any]) -> str:
    canonical = json.dumps(config, sort_keys=True, ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()[:12]


def _robust_sigma(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float32)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 1.0
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)))
    sigma = 1.4826 * mad
    if not np.isfinite(sigma) or sigma <= np.finfo(np.float32).eps:
        sigma = float(np.std(values))
    return max(sigma, 1e-6)


def _largest_component(mask: np.ndarray) -> np.ndarray:
    connectivity = ndi.generate_binary_structure(mask.ndim, mask.ndim)
    labels, number = ndi.label(mask, structure=connectivity)
    if number == 0:
        return np.zeros_like(mask, dtype=bool)
    counts = np.bincount(labels.ravel())
    counts[0] = 0
    return labels == int(np.argmax(counts))


def _clean_optic_mask_2d(
    mask: np.ndarray,
    minimum_area_px: int,
    closing_radius_px: int,
    exclusion_margin_px: int,
) -> np.ndarray:
    """Regularize one optic-lobe mask without changing its z support."""
    result = np.asarray(mask, dtype=bool)
    if closing_radius_px > 0 and result.any():
        result = binary_closing(result, footprint=disk(closing_radius_px))
    result = remove_small_objects(result, min_size=minimum_area_px)
    if result.any():
        result = ndi.binary_fill_holes(result)
    if exclusion_margin_px > 0 and result.any():
        result = binary_dilation(result, footprint=disk(exclusion_margin_px))
    return np.asarray(result, dtype=bool)


def _border_connected_high_background_mask(
    broad_background: np.ndarray,
    candidate_domain: np.ndarray,
    image_border: np.ndarray,
    seed_quantile: float,
    support_fraction: float,
    threshold_offset_sigma: float,
    minimum_area_px: int,
) -> tuple[np.ndarray, float | None, float | None]:
    """Return only high-background pixels connected to the image border.

    Connectivity is evaluated inside an intensity-restricted support mask, not
    inside the complete tissue mask.  Consequently, a thin tissue bridge
    between central brain and optic lobe cannot cause the complete joined
    tissue component to be discarded.
    """
    domain = np.asarray(candidate_domain, dtype=bool)
    empty = np.zeros_like(domain, dtype=bool)
    if not domain.any():
        return empty, None, None
    values = np.asarray(broad_background[domain], dtype=np.float32)
    sampled = values[:: max(1, values.size // 100_000)]
    try:
        threshold = float(threshold_otsu(sampled))
    except ValueError:
        return empty, None, None
    lower_values = sampled[sampled <= np.median(sampled)]
    threshold += float(threshold_offset_sigma) * _robust_sigma(lower_values)
    seed_threshold = float(np.quantile(sampled, float(seed_quantile)))
    support_threshold = threshold + float(support_fraction) * max(
        seed_threshold - threshold,
        0.0,
    )
    optic_seed = domain & (broad_background >= seed_threshold)
    optic_support = domain & (broad_background >= support_threshold)
    optic_candidate = ndi.binary_propagation(optic_seed, mask=optic_support)
    optic_labels, optic_number = ndi.label(optic_candidate)
    optic = np.zeros_like(domain, dtype=bool)
    if optic_number:
        counts = np.bincount(optic_labels.ravel())
        border_components = set(np.unique(optic_labels[image_border]).tolist())
        border_components.discard(0)
        for component_index in border_components:
            if counts[component_index] >= minimum_area_px:
                optic |= optic_labels == component_index
    return optic, support_threshold, seed_threshold


def _select_terminal_rescue_components(
    rescue_candidate: np.ndarray,
    connection_region: np.ndarray,
    maximum_area_px: int,
    keep_multiple: bool,
) -> np.ndarray:
    """Keep small reachable tissue components that remain close to central brain."""
    labels, number = ndi.label(rescue_candidate)
    selected = np.zeros_like(rescue_candidate, dtype=bool)
    if number == 0 or maximum_area_px <= 0:
        return selected
    counts = np.bincount(labels.ravel())
    overlaps = np.bincount(labels[connection_region].ravel(), minlength=number + 1)
    component_indices = [
        index
        for index in range(1, number + 1)
        if overlaps[index] > 0 and counts[index] <= maximum_area_px
    ]
    component_indices.sort(key=lambda index: (int(overlaps[index]), -int(counts[index])), reverse=True)
    selected_area = 0
    for component_index in component_indices:
        component_area = int(counts[component_index])
        if selected_area + component_area > maximum_area_px:
            continue
        selected |= labels == component_index
        selected_area += component_area
        if not keep_multiple:
            break
    return selected


def _build_auto_tissue_roi(
    reader: Any,
    channel: int,
    config: dict[str, Any],
) -> tuple[np.ndarray, tuple[float, float, float], int]:
    roi_cfg = config["roi"]
    level = int(roi_cfg["auto_level"])
    if level >= reader.ResolutionLevels:
        level = reader.ResolutionLevels - 1
    low_volume = np.asarray(reader[level, 0, channel, :, :, :], dtype=np.float32)
    low_spacing = tuple(float(x) for x in reader.metaData[level, 0, channel, "resolution"])
    sigma_px = max(float(roi_cfg["smooth_sigma_um"]) / low_spacing[1], 0.0)
    minimum_area_px = max(
        1,
        int(round(float(roi_cfg["minimum_component_area_um2"]) / (low_spacing[1] * low_spacing[2]))),
    )
    close_radius_px = max(1, int(round(float(roi_cfg["closing_radius_um"]) / low_spacing[1])))
    expansion_px = max(
        0,
        int(round(float(roi_cfg.get("tissue_expansion_um", 0.0)) / low_spacing[1])),
    )
    contrast_high_quantile = float(roi_cfg.get("contrast_high_quantile", 0.90))
    if not 0.50 < contrast_high_quantile <= 1.0:
        raise ValueError("roi.contrast_high_quantile must be in (0.50, 1.0]")
    roi = np.zeros_like(low_volume, dtype=bool)

    for z_index, plane in enumerate(low_volume):
        smooth = ndi.gaussian_filter(plane, sigma=sigma_px, mode="nearest")
        sampled = smooth[::2, ::2]
        p50 = float(np.median(sampled))
        p_high = float(np.quantile(sampled, contrast_high_quantile))
        low_half = sampled[sampled <= p50]
        raw_noise = _robust_sigma(low_half)
        contrast_snr = float((p_high - p50) / raw_noise)
        if contrast_snr < float(roi_cfg["minimum_contrast_snr"]):
            continue
        try:
            threshold = float(threshold_otsu(sampled))
        except ValueError:
            continue
        current = smooth > threshold
        current = remove_small_objects(current, min_size=minimum_area_px)
        if expansion_px > 0 and current.any():
            background_median = float(np.median(low_half)) if low_half.size else float(p50)
            support_threshold = background_median + float(
                roi_cfg.get("tissue_support_snr", 4.0)
            ) * raw_noise
            nearby = binary_dilation(current, footprint=disk(expansion_px))
            current |= nearby & (smooth > support_threshold)
        current = binary_closing(current, footprint=disk(close_radius_px))
        current = ndi.binary_fill_holes(current)
        roi[z_index] = current

    # Retain the principal continuous tissue volume.  In central-priority mode,
    # also preserve terminal candidate pieces near its XY history: a genuine
    # small end section may be disconnected by one weak thresholded z plane.
    tissue_candidates = roi
    principal_tissue = _largest_component(tissue_candidates)
    exclusion_mode = str(
        roi_cfg.get("optic_exclusion_mode", "z_local_central_priority")
    ).lower()
    retain_terminal_components = bool(
        roi_cfg.get("retain_terminal_tissue_components", True)
    )
    if (
        exclusion_mode == "z_local_central_priority"
        and retain_terminal_components
        and principal_tissue.any()
    ):
        principal_areas = principal_tissue.reshape(principal_tissue.shape[0], -1).sum(axis=1)
        maximum_principal_area = max(1, int(principal_areas.max()))
        terminal_fraction = float(roi_cfg.get("terminal_area_fraction", 0.30))
        terminal_threshold = int(round(terminal_fraction * maximum_principal_area))
        active_principal_z = np.flatnonzero(principal_areas > 0)
        context_planes = max(
            1,
            int(
                math.ceil(
                    float(roi_cfg.get("optic_z_context_um", 6.0)) / low_spacing[0]
                )
            ),
        )
        valid_terminal_z = np.zeros(principal_tissue.shape[0], dtype=bool)
        if active_principal_z.size:
            z0 = max(0, int(active_principal_z.min()) - context_planes)
            z1 = min(principal_tissue.shape[0] - 1, int(active_principal_z.max()) + context_planes)
            valid_terminal_z[z0 : z1 + 1] = True
        terminal_z = (principal_areas <= terminal_threshold) & valid_terminal_z
        rescue_px = max(
            1,
            int(
                round(
                    float(roi_cfg.get("central_terminal_rescue_um", 12.0))
                    / low_spacing[1]
                )
            ),
        )
        principal_xy_support = binary_dilation(
            principal_tissue.any(axis=0),
            footprint=disk(rescue_px),
        )
        recovered_terminal = (
            tissue_candidates
            & terminal_z[:, None, None]
            & principal_xy_support[None, :, :]
        )
        roi = principal_tissue | recovered_terminal
    else:
        roi = principal_tissue
    erosion_px = int(round(float(roi_cfg["erosion_margin_um"]) / low_spacing[1]))
    if erosion_px > 0:
        footprint = disk(erosion_px)
        for z_index in range(roi.shape[0]):
            if roi[z_index].any():
                roi[z_index] = binary_erosion(roi[z_index], footprint=footprint)
    return roi, low_spacing, level


def _build_auto_central_brain_roi(
    reader: Any,
    channel: int,
    config: dict[str, Any],
) -> tuple[np.ndarray, tuple[float, float, float], int, dict[str, Any]]:
    """Separate lower-background central brain from bright optic lobes.

    The first stage recovers complete tissue without erosion.  A much broader
    Gaussian estimates autofluorescent background while ignoring LD-sized
    puncta.  Border-connected high-background regions form per-plane optic
    candidates.  ``fixed_xy`` reproduces the original stack-wide footprint;
    ``z_local_central_priority`` stabilizes optic masks only over a local z
    window and protects small terminal tissue that remains continuous with the
    preceding central-brain plane.
    """
    roi_cfg = config["roi"]
    tissue_config = _deep_merge(config, {"roi": {"erosion_margin_um": 0.0}})
    tissue, low_spacing, level = _build_auto_tissue_roi(reader, channel, tissue_config)
    low_volume = np.asarray(reader[level, 0, channel, :, :, :], dtype=np.float32)
    # A whole-plane Otsu threshold can miss a small terminal central section
    # when a bright optic lobe dominates the intensity histogram.  Build a
    # permissive local-background support volume for terminal rescue only; it
    # is never used to expand the confident middle-z ROI.
    tissue_sigma_px = max(float(roi_cfg["smooth_sigma_um"]) / low_spacing[1], 0.0)
    terminal_support_snr = float(roi_cfg.get("terminal_tissue_support_snr", 4.0))
    terminal_minimum_area_px = max(
        1,
        int(
            round(
                float(roi_cfg.get("terminal_minimum_component_area_um2", 100.0))
                / (low_spacing[1] * low_spacing[2])
            )
        ),
    )
    terminal_tissue_closing_px = max(
        0,
        int(round(float(roi_cfg.get("terminal_tissue_closing_um", 1.0)) / low_spacing[1])),
    )
    terminal_tissue_support = np.zeros_like(tissue, dtype=bool)
    for z_index, plane in enumerate(low_volume):
        smooth_tissue = ndi.gaussian_filter(plane, sigma=tissue_sigma_px, mode="nearest")
        sampled_tissue = smooth_tissue[::2, ::2]
        background_median = float(np.median(sampled_tissue))
        lower_half = sampled_tissue[sampled_tissue <= background_median]
        support_threshold = background_median + terminal_support_snr * _robust_sigma(lower_half)
        support = remove_small_objects(
            smooth_tissue > support_threshold,
            min_size=terminal_minimum_area_px,
        )
        if terminal_tissue_closing_px > 0 and support.any():
            support = binary_closing(
                support,
                footprint=disk(terminal_tissue_closing_px),
            )
        terminal_tissue_support[z_index] = support
    broad_sigma_px = max(float(roi_cfg["optic_background_sigma_um"]) / low_spacing[1], 0.5)
    close_radius_px = max(
        1,
        int(round(float(roi_cfg["central_closing_radius_um"]) / low_spacing[1])),
    )
    maximum_hole_area_px = max(
        1,
        int(round(float(roi_cfg["central_hole_maximum_area_um2"]) / (low_spacing[1] * low_spacing[2]))),
    )
    thresholds: list[float | None] = []
    seed_thresholds: list[float | None] = []
    terminal_thresholds: list[float | None] = []
    terminal_seed_thresholds: list[float | None] = []
    broad_background = np.zeros_like(low_volume, dtype=np.float32)
    optic_by_plane = np.zeros_like(tissue, dtype=bool)
    terminal_optic_by_plane = np.zeros_like(tissue, dtype=bool)
    border_distance_px = max(
        1,
        int(round(float(roi_cfg["optic_image_border_distance_um"]) / low_spacing[1])),
    )
    image_border = np.zeros(tissue.shape[1:], dtype=bool)
    image_border[:border_distance_px, :] = True
    image_border[-border_distance_px:, :] = True
    image_border[:, :border_distance_px] = True
    image_border[:, -border_distance_px:] = True
    minimum_optic_area_px = max(
        1,
        int(round(float(roi_cfg["optic_minimum_area_um2"]) / (low_spacing[1] * low_spacing[2]))),
    )

    for z_index, plane in enumerate(low_volume):
        broad = ndi.gaussian_filter(plane, sigma=broad_sigma_px, mode="nearest")
        broad_background[z_index] = broad
        optic, support_threshold, seed_threshold = _border_connected_high_background_mask(
            broad,
            tissue[z_index],
            image_border,
            float(roi_cfg["optic_seed_quantile"]),
            float(roi_cfg["optic_support_fraction"]),
            float(roi_cfg["optic_threshold_offset_sigma"]),
            minimum_optic_area_px,
        )
        thresholds.append(support_threshold)
        seed_thresholds.append(seed_threshold)
        optic_by_plane[z_index] = optic

        # The permissive terminal tissue support can contain a bright optic
        # lobe that was absent from the strict tissue mask.  Detect its
        # high-background pixels explicitly.  Crucially, region growing stays
        # inside the high-background support and never floods through the
        # complete tissue component into central brain.
        terminal_domain = tissue[z_index] | terminal_tissue_support[z_index]
        terminal_optic, terminal_threshold, terminal_seed_threshold = (
            _border_connected_high_background_mask(
                broad,
                terminal_domain,
                image_border,
                float(
                    roi_cfg.get(
                        "terminal_optic_seed_quantile",
                        roi_cfg["optic_seed_quantile"],
                    )
                ),
                float(
                    roi_cfg.get(
                        "terminal_optic_support_fraction",
                        roi_cfg["optic_support_fraction"],
                    )
                ),
                float(roi_cfg["optic_threshold_offset_sigma"]),
                minimum_optic_area_px,
            )
        )
        terminal_thresholds.append(terminal_threshold)
        terminal_seed_thresholds.append(terminal_seed_threshold)
        terminal_optic_by_plane[z_index] = terminal_optic

    exclusion_mode = str(
        roi_cfg.get("optic_exclusion_mode", "z_local_central_priority")
    ).lower()
    if exclusion_mode not in {"fixed_xy", "z_local_central_priority"}:
        raise ValueError(
            "roi.optic_exclusion_mode must be 'fixed_xy' or "
            "'z_local_central_priority'"
        )

    # Keep the historical stack-wide footprint available for direct
    # comparison and backward-compatible fixed_xy runs.
    minimum_z_planes = max(1, int(roi_cfg.get("optic_minimum_z_planes", 2)))
    optic_presence = optic_by_plane.sum(axis=0)
    optic_footprint = optic_presence >= minimum_z_planes
    footprint_closing_px = max(
        0,
        int(round(float(roi_cfg.get("optic_footprint_closing_um", 0.0)) / low_spacing[1])),
    )
    optic_footprint = _clean_optic_mask_2d(
        optic_footprint,
        minimum_optic_area_px,
        footprint_closing_px,
        0,
    )
    exclusion_margin_px = max(
        0,
        int(round(float(roi_cfg["optic_exclusion_margin_um"]) / low_spacing[1])),
    )
    optic_exclusion_fixed = optic_footprint.copy()
    if exclusion_margin_px > 0 and optic_exclusion_fixed.any():
        optic_exclusion_fixed = binary_dilation(
            optic_exclusion_fixed,
            footprint=disk(exclusion_margin_px),
        )

    # Local z support follows changing anatomy instead of projecting every
    # optic pixel through the entire stack.  A directly detected optic region
    # is always retained; recurrence fills ambiguous gaps only nearby in z.
    optic_z_context_um = max(0.0, float(roi_cfg.get("optic_z_context_um", 6.0)))
    optic_z_context_planes = max(0, int(math.ceil(optic_z_context_um / low_spacing[0])))
    z_kernel = np.ones(2 * optic_z_context_planes + 1, dtype=np.uint16)
    local_presence = ndi.convolve1d(
        optic_by_plane.astype(np.uint16),
        z_kernel,
        axis=0,
        mode="constant",
        cval=0,
    )
    optic_local_core = optic_by_plane | (local_presence >= minimum_z_planes)
    local_closing_px = max(
        0,
        int(round(float(roi_cfg.get("optic_local_closing_um", 2.0)) / low_spacing[1])),
    )
    local_margin_px = max(
        0,
        int(
            round(
                float(roi_cfg.get("optic_local_exclusion_margin_um", 2.0))
                / low_spacing[1]
            )
        ),
    )
    terminal_closing_px = max(
        0,
        int(round(float(roi_cfg.get("terminal_optic_closing_um", 0.0)) / low_spacing[1])),
    )
    terminal_margin_px = max(
        0,
        int(
            round(
                float(roi_cfg.get("terminal_optic_exclusion_margin_um", 0.0))
                / low_spacing[1]
            )
        ),
    )
    optic_exclusion_local = np.zeros_like(optic_by_plane, dtype=bool)
    optic_exclusion_terminal = np.zeros_like(optic_by_plane, dtype=bool)
    use_expanded_terminal_optic = bool(
        roi_cfg.get("terminal_exclude_border_connected_tissue", True)
    )
    for z_index in range(optic_by_plane.shape[0]):
        optic_exclusion_local[z_index] = _clean_optic_mask_2d(
            optic_local_core[z_index],
            minimum_optic_area_px,
            local_closing_px,
            local_margin_px,
        )
        # At terminal planes, only direct high-background evidence from the
        # current plane is allowed to remove tissue.  When enabled, the source
        # mask is detected over permissive terminal tissue support so a bright
        # edge optic lobe is still seen; connectivity itself is never grounds
        # for deleting the complete tissue component.
        terminal_optic_source = (
            terminal_optic_by_plane[z_index]
            if use_expanded_terminal_optic
            else optic_by_plane[z_index]
        )
        optic_exclusion_terminal[z_index] = _clean_optic_mask_2d(
            terminal_optic_source,
            minimum_optic_area_px,
            terminal_closing_px,
            terminal_margin_px,
        )

    if exclusion_mode == "fixed_xy":
        normal_optic_exclusion = np.broadcast_to(
            optic_exclusion_fixed[None, :, :],
            tissue.shape,
        )
    else:
        normal_optic_exclusion = optic_exclusion_local

    central_allowed = tissue & ~normal_optic_exclusion
    for z_index in range(central_allowed.shape[0]):
        if central_allowed[z_index].any():
            central_allowed[z_index] = remove_small_holes(
                binary_closing(central_allowed[z_index], footprint=disk(close_radius_px)),
                area_threshold=maximum_hole_area_px,
            )

    # Track the central brain from its largest cross-section toward both z
    # ends.  A genuine anatomical volume changes continuously; a weak depth
    # gradient misclassified as tissue can otherwise create an abrupt, very
    # large ROI in the last one or two planes.
    areas_by_z = central_allowed.reshape(central_allowed.shape[0], -1).sum(axis=1)
    anchor_z = int(np.argmax(areas_by_z))
    maximum_z_expansion_px = max(
        1,
        int(round(float(roi_cfg.get("central_z_max_expansion_um", 8.0)) / low_spacing[1])),
    )
    z_expansion_footprint = disk(maximum_z_expansion_px)
    central = np.zeros_like(central_allowed, dtype=bool)
    central[anchor_z] = _largest_component(central_allowed[anchor_z])
    anchor_area_px = max(1, int(central[anchor_z].sum()))
    terminal_area_fraction = float(roi_cfg.get("terminal_area_fraction", 0.30))
    if not 0.0 <= terminal_area_fraction <= 1.0:
        raise ValueError("roi.terminal_area_fraction must be in [0, 1]")
    terminal_area_threshold_px = int(round(terminal_area_fraction * anchor_area_px))
    terminal_rescue_px = max(
        1,
        int(round(float(roi_cfg.get("central_terminal_rescue_um", 12.0)) / low_spacing[1])),
    )
    terminal_rescue_footprint = disk(terminal_rescue_px)
    terminal_maximum_area_growth = float(
        roi_cfg.get("terminal_maximum_area_growth", 1.5)
    )
    if terminal_maximum_area_growth < 1.0:
        raise ValueError("roi.terminal_maximum_area_growth must be >= 1")
    terminal_keep_multiple = bool(
        roi_cfg.get("terminal_keep_multiple_components", True)
    )
    terminal_planes = np.zeros(central.shape[0], dtype=bool)
    terminal_reference_regions = np.zeros_like(central, dtype=bool)
    terminal_rescued_voxels = 0

    for direction in (-1, 1):
        z_index = anchor_z + direction
        terminal_active = False
        # Accumulate only confident, non-terminal central masks.  Once terminal
        # mode begins this XY reference is frozen, preventing the rescue path
        # from drifting laterally into an optic lobe.
        central_reference_footprint = central[anchor_z].copy()
        while 0 <= z_index < central.shape[0]:
            previous = central[z_index - direction]
            if not previous.any():
                break
            reachable = binary_dilation(previous, footprint=z_expansion_footprint)
            current = central_allowed[z_index] & reachable
            if current.any():
                current = remove_small_holes(
                    binary_closing(current, footprint=disk(close_radius_px)),
                    area_threshold=maximum_hole_area_px,
                )

            if exclusion_mode == "z_local_central_priority":
                if int(current.sum()) <= terminal_area_threshold_px:
                    terminal_active = True
                if terminal_active:
                    terminal_planes[z_index] = True
                    previous_area_px = max(1, int(previous.sum()))
                    maximum_terminal_area_px = max(
                        terminal_area_threshold_px,
                        int(round(previous_area_px * terminal_maximum_area_growth)),
                    )
                    rescue_corridor = binary_dilation(
                        central_reference_footprint,
                        footprint=terminal_rescue_footprint,
                    )
                    terminal_reference_regions[z_index] = rescue_corridor
                    # Current-plane border-connected optic evidence remains a
                    # hard exclusion.  Central priority overrides only optic
                    # support propagated from other z planes.
                    terminal_candidate = (
                        (tissue[z_index] | terminal_tissue_support[z_index])
                        & rescue_corridor
                        & ~optic_exclusion_terminal[z_index]
                    )
                    protected = _select_terminal_rescue_components(
                        terminal_candidate,
                        rescue_corridor,
                        maximum_terminal_area_px,
                        terminal_keep_multiple,
                    )
                    terminal_rescued_voxels += int((protected & ~current).sum())
                    current = protected

                    if current.any():
                        current = remove_small_holes(
                            binary_closing(current, footprint=disk(close_radius_px)),
                            area_threshold=maximum_hole_area_px,
                        )
                else:
                    central_reference_footprint |= current
            central[z_index] = current
            z_index += direction

    # Repair a short, catastrophic ROI collapse only when it is bracketed by
    # substantially larger central masks on both sides.  Candidate pixels must
    # be reachable from both flanks, supported by tissue in the missing plane,
    # and outside the current-plane high-background optic mask.  This guards
    # against one or two low-resolution topology errors without extrapolating
    # tissue beyond a genuine anatomical z end.
    short_gap_maximum_um = max(
        0.0,
        float(roi_cfg.get("terminal_short_gap_maximum_um", 6.0)),
    )
    short_gap_maximum_planes = int(
        math.floor(short_gap_maximum_um / low_spacing[0] + 1e-9)
    )
    short_gap_area_fraction = float(
        roi_cfg.get("terminal_short_gap_area_fraction", 0.25)
    )
    if not 0.0 <= short_gap_area_fraction < 1.0:
        raise ValueError("roi.terminal_short_gap_area_fraction must be in [0, 1)")
    short_gap_rescued_planes: list[int] = []
    short_gap_rescued_voxels = 0
    if (
        exclusion_mode == "z_local_central_priority"
        and short_gap_maximum_planes > 0
    ):
        initial_areas = central.reshape(central.shape[0], -1).sum(axis=1)
        gap_candidates: list[tuple[float, int, int]] = []
        for gap_length in range(1, short_gap_maximum_planes + 1):
            for left_z in range(0, central.shape[0] - gap_length - 1):
                right_z = left_z + gap_length + 1
                flank_area = min(
                    int(initial_areas[left_z]),
                    int(initial_areas[right_z]),
                )
                if flank_area < terminal_minimum_area_px:
                    continue
                gap_peak_area = int(initial_areas[left_z + 1 : right_z].max())
                if gap_peak_area > short_gap_area_fraction * flank_area:
                    continue
                flank_overlap = binary_dilation(
                    central[left_z],
                    footprint=z_expansion_footprint,
                ) & binary_dilation(
                    central[right_z],
                    footprint=z_expansion_footprint,
                )
                if not flank_overlap.any():
                    continue
                collapse_ratio = flank_area / max(1, gap_peak_area)
                gap_candidates.append((collapse_ratio, left_z, right_z))

        used_gap_planes: set[int] = set()
        for _, left_z, right_z in sorted(gap_candidates, reverse=True):
            gap_indices = list(range(left_z + 1, right_z))
            if any(index in used_gap_planes for index in gap_indices):
                continue
            left_reach = binary_dilation(
                central[left_z],
                footprint=z_expansion_footprint,
            )
            right_reach = binary_dilation(
                central[right_z],
                footprint=z_expansion_footprint,
            )
            interpolation_envelope = left_reach & right_reach
            maximum_gap_area_px = int(
                round(
                    max(int(initial_areas[left_z]), int(initial_areas[right_z]))
                    * terminal_maximum_area_growth
                )
            )
            replacements: list[tuple[int, np.ndarray]] = []
            for gap_z in gap_indices:
                rescue_candidate = (
                    (tissue[gap_z] | terminal_tissue_support[gap_z])
                    & interpolation_envelope
                    & ~optic_exclusion_terminal[gap_z]
                )
                rescued = _select_terminal_rescue_components(
                    rescue_candidate,
                    interpolation_envelope,
                    maximum_gap_area_px,
                    terminal_keep_multiple,
                )
                if rescued.any():
                    rescued = remove_small_holes(
                        binary_closing(rescued, footprint=disk(close_radius_px)),
                        area_threshold=maximum_hole_area_px,
                    )
                if int(rescued.sum()) <= int(initial_areas[gap_z]):
                    replacements = []
                    break
                replacements.append((gap_z, rescued))
            if not replacements:
                continue
            for gap_z, rescued in replacements:
                short_gap_rescued_voxels += int((rescued & ~central[gap_z]).sum())
                central[gap_z] = rescued
                terminal_planes[gap_z] = True
                terminal_reference_regions[gap_z] |= interpolation_envelope
                short_gap_rescued_planes.append(gap_z)
                used_gap_planes.add(gap_z)

    # Restore interior islands removed by locally bright somata/LDs.
    for z_index in range(central.shape[0]):
        if central[z_index].any():
            central[z_index] = remove_small_holes(
                binary_closing(central[z_index], footprint=disk(close_radius_px)),
                area_threshold=maximum_hole_area_px,
            )

    expansion_px = int(round(float(roi_cfg["central_expansion_um"]) / low_spacing[1]))
    if expansion_px > 0:
        footprint = disk(expansion_px)
        for z_index in range(central.shape[0]):
            if not central[z_index].any():
                continue
            expanded = binary_dilation(central[z_index], footprint=footprint)
            if exclusion_mode == "fixed_xy":
                permitted = tissue[z_index] & ~optic_exclusion_fixed
            elif terminal_planes[z_index]:
                permitted = (
                    (tissue[z_index] | terminal_tissue_support[z_index])
                    & terminal_reference_regions[z_index]
                    & ~optic_exclusion_terminal[z_index]
                )
            else:
                permitted = tissue[z_index] & ~optic_exclusion_local[z_index]
            central[z_index] |= expanded & permitted

    # fixed_xy historically enforced one strictly connected 3-D component.
    # The central-priority mode intentionally retains controlled terminal
    # pieces even when a weak intervening plane creates a one-plane gap.
    if exclusion_mode == "fixed_xy":
        central = _largest_component(central)
    diagnostics = {
        "contrast_high_quantile": float(roi_cfg.get("contrast_high_quantile", 0.90)),
        "optic_background_sigma_um": float(roi_cfg["optic_background_sigma_um"]),
        "thresholds_by_lowres_z": thresholds,
        "seed_thresholds_by_lowres_z": seed_thresholds,
        "terminal_thresholds_by_lowres_z": terminal_thresholds,
        "terminal_seed_thresholds_by_lowres_z": terminal_seed_thresholds,
        "optic_exclusion_mode": exclusion_mode,
        "optic_candidate_voxels": int(optic_by_plane.sum()),
        "terminal_optic_candidate_voxels": int(terminal_optic_by_plane.sum()),
        "terminal_expanded_border_optic_detection": use_expanded_terminal_optic,
        "optic_footprint_pixels": int(optic_footprint.sum()),
        "optic_exclusion_pixels": int(optic_exclusion_fixed.sum()),
        "optic_local_exclusion_voxels": int(optic_exclusion_local.sum()),
        "optic_terminal_exclusion_voxels": int(optic_exclusion_terminal.sum()),
        "terminal_tissue_support_voxels": int(terminal_tissue_support.sum()),
        "optic_minimum_z_planes": minimum_z_planes,
        "optic_z_context_um": optic_z_context_um,
        "optic_z_context_lowres_planes": optic_z_context_planes,
        "central_anchor_lowres_z": anchor_z,
        "central_z_max_expansion_um": float(
            roi_cfg.get("central_z_max_expansion_um", 8.0)
        ),
        "terminal_area_fraction": terminal_area_fraction,
        "terminal_lowres_z_indices": np.flatnonzero(terminal_planes).astype(int).tolist(),
        "terminal_rescued_voxels": int(terminal_rescued_voxels),
        "terminal_short_gap_maximum_um": short_gap_maximum_um,
        "terminal_short_gap_area_fraction": short_gap_area_fraction,
        "terminal_short_gap_rescued_lowres_z_indices": sorted(
            set(short_gap_rescued_planes)
        ),
        "terminal_short_gap_rescued_voxels": int(short_gap_rescued_voxels),
        "tissue_voxels": int(tissue.sum()),
        "central_voxels": int(central.sum()),
        "retained_fraction": float(central.sum() / tissue.sum()) if tissue.any() else 0.0,
    }
    return central, low_spacing, level, diagnostics


class ROIProvider:
    def __init__(
        self,
        mode: str,
        full_shape_zyx: tuple[int, int, int],
        full_spacing_zyx: tuple[float, float, float],
        low_mask: np.ndarray | None = None,
        low_spacing_zyx: tuple[float, float, float] | None = None,
        external_mask: np.ndarray | None = None,
        crop_yx: list[int] | tuple[int, int, int, int] | None = None,
    ) -> None:
        self.mode = mode
        self.full_shape = full_shape_zyx
        self.full_spacing = full_spacing_zyx
        self.low_mask = low_mask
        self.low_spacing = low_spacing_zyx
        self.external_mask = external_mask
        self.crop_yx = crop_yx
        self.roi_artifacts = None

    def plane(self, z_index: int) -> np.ndarray:
        height, width = self.full_shape[1:]
        if self.mode == "full":
            result = np.ones((height, width), dtype=bool)
        elif self.mode in {"auto_tissue", "auto_central_brain", "guided_central_brain"}:
            if self.low_mask is None or self.low_spacing is None:
                raise RuntimeError("Automatic ROI was not initialized")
            z_um = z_index * self.full_spacing[0]
            low_z = int(round(z_um / self.low_spacing[0]))
            low_z = int(np.clip(low_z, 0, self.low_mask.shape[0] - 1))
            result = resize(
                self.low_mask[low_z].astype(np.uint8),
                (height, width),
                order=0,
                preserve_range=True,
                anti_aliasing=False,
            ).astype(bool)
        elif self.mode == "mask":
            if self.external_mask is None:
                raise RuntimeError("External ROI mask was not initialized")
            if self.external_mask.ndim == 2:
                result = self.external_mask.astype(bool)
            else:
                mask_z = int(round(z_index * self.external_mask.shape[0] / self.full_shape[0]))
                mask_z = int(np.clip(mask_z, 0, self.external_mask.shape[0] - 1))
                result = self.external_mask[mask_z].astype(bool)
            if result.shape != (height, width):
                result = resize(
                    result.astype(np.uint8),
                    (height, width),
                    order=0,
                    preserve_range=True,
                    anti_aliasing=False,
                ).astype(bool)
        else:
            raise ValueError(f"Unknown ROI mode: {self.mode}")

        if self.crop_yx is not None:
            y0, y1, x0, x1 = (int(v) for v in self.crop_yx)
            crop_mask = np.zeros_like(result)
            crop_mask[max(y0, 0) : min(y1, height), max(x0, 0) : min(x1, width)] = True
            result &= crop_mask
        return result

    def active_z_indices(self) -> np.ndarray:
        if self.mode in {"auto_tissue", "auto_central_brain", "guided_central_brain"} and self.low_mask is not None and self.low_spacing is not None:
            active_low = np.flatnonzero(self.low_mask.reshape(self.low_mask.shape[0], -1).any(axis=1))
            if active_low.size:
                z0_um = active_low.min() * self.low_spacing[0]
                z1_um = active_low.max() * self.low_spacing[0]
                z0 = max(0, int(math.floor(z0_um / self.full_spacing[0])))
                z1 = min(self.full_shape[0] - 1, int(math.ceil(z1_um / self.full_spacing[0])))
                return np.arange(z0, z1 + 1)
        return np.arange(self.full_shape[0])


def _load_external_mask(path: str | Path) -> np.ndarray:
    import tifffile

    mask = np.asarray(tifffile.imread(path))
    if mask.ndim not in (2, 3):
        raise ValueError(f"ROI mask must be 2-D or 3-D; got shape {mask.shape}")
    return mask > 0


def _save_roi_mask_tiff(path: str | Path, provider: ROIProvider,
                        spacing_zyx: tuple[float, float, float], *,
                        compression: str = "deflate") -> Path:
    """Stream the same ZYX ROI planes used by analysis to an ImageJ stack."""
    import tifffile

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    z_um, y_um, x_um = spacing_zyx
    tifffile.imwrite(
        output,
        data=(provider.plane(z).astype(np.uint8) * 255 for z in range(provider.full_shape[0])),
        shape=provider.full_shape,
        dtype=np.uint8,
        imagej=True,
        compression=compression,
        resolution=(1 / x_um, 1 / y_um),
        metadata={"axes": "ZYX", "spacing": z_um, "unit": "um"},
    )
    return output


def _make_roi_provider(reader: Any, config: dict[str, Any], *, keep_roi_stages: bool = False) -> tuple[ROIProvider, dict[str, Any]]:
    roi_cfg = config["roi"]
    input_cfg = config["input"]
    channel = int(input_cfg["channel_index"])
    full_shape = tuple(int(x) for x in reader.shape[-3:])
    full_spacing = tuple(float(x) for x in reader.resolution)
    mode = str(roi_cfg["mode"])
    metadata: dict[str, Any] = {"mode": mode}
    roi_artifacts = None
    if mode in {"auto_tissue", "auto_central_brain", "guided_central_brain"}:
        method = roi_cfg.get("central_method", "anatomical")
        if mode in {"auto_central_brain", "guided_central_brain"} and method not in {"anatomical", "anatomical_v2", "intensity_legacy"}:
            raise ValueError("roi.central_method must be anatomical, anatomical_v2 or intensity_legacy")
        if mode in {"auto_central_brain", "guided_central_brain"} and method == "anatomical_v2":
            from .anatomical_v2 import build_anatomical_v2_roi
            low_mask, low_spacing, level, central_diagnostics, roi_artifacts = build_anatomical_v2_roi(reader, config, keep_stages=keep_roi_stages)
            for message in central_diagnostics["review_flags"]:
                warnings.warn(message, UserWarning, stacklevel=2)
        elif mode == "guided_central_brain" or (mode == "auto_central_brain" and method == "anatomical"):
            low_mask, low_spacing, level, central_diagnostics = build_anatomical_roi(reader, config)
        elif mode == "auto_central_brain":
            low_mask, low_spacing, level, central_diagnostics = _build_auto_central_brain_roi(
                reader,
                channel,
                config,
            )
        else:
            low_mask, low_spacing, level = _build_auto_tissue_roi(reader, channel, config)
            central_diagnostics = {}
        metadata.update(
            {
                "auto_level": level,
                "low_shape": list(low_mask.shape),
                "low_spacing_zyx_um": list(low_spacing),
                "low_roi_voxels": int(low_mask.sum()),
                **central_diagnostics,
            }
        )
        provider = ROIProvider(
            mode,
            full_shape,
            full_spacing,
            low_mask=low_mask,
            low_spacing_zyx=low_spacing,
            crop_yx=roi_cfg.get("crop_yx"),
        )
        provider.roi_artifacts = roi_artifacts
    elif mode == "mask":
        mask_path = roi_cfg.get("mask_path")
        if not mask_path:
            raise ValueError("roi.mask_path is required when roi.mode is 'mask'")
        external = _load_external_mask(mask_path)
        metadata["mask_path"] = str(Path(mask_path).resolve())
        metadata["mask_shape"] = list(external.shape)
        provider = ROIProvider(
            mode,
            full_shape,
            full_spacing,
            external_mask=external,
            crop_yx=roi_cfg.get("crop_yx"),
        )
    elif mode == "full":
        provider = ROIProvider(
            mode,
            full_shape,
            full_spacing,
            crop_yx=roi_cfg.get("crop_yx"),
        )
    else:
        raise ValueError("roi.mode must be one of: auto_central_brain, guided_central_brain, auto_tissue, mask, full")
    return provider, metadata


def _detect_blobs(
    snr: np.ndarray,
    foreground: np.ndarray,
    roi: np.ndarray,
    xy_spacing_um: float,
    config: dict[str, Any],
) -> pd.DataFrame:
    detect_cfg = config["detection"]
    minimum_diameter = float(detect_cfg["minimum_diameter_um"])
    maximum_diameter = float(detect_cfg["maximum_diameter_um"])
    number_of_scales = int(detect_cfg["number_of_scales"])
    downsample = max(1, int(detect_cfg.get("processing_downsample", 1)))
    detection_spacing_um = xy_spacing_um * downsample
    diameters_um = np.geomspace(minimum_diameter, maximum_diameter, number_of_scales)
    sigmas_px = diameters_um / (2.0 * math.sqrt(2.0) * detection_spacing_um)
    work = np.clip(snr, 0, float(config["preprocess"]["snr_clip_for_detection"])).astype(np.float32)
    valid_full = roi & foreground
    if downsample > 1:
        # Mild anti-aliasing retains a ~0.9 µm droplet as a ~3-pixel object.
        work_detection = ndi.gaussian_filter(work, sigma=0.5 * downsample, mode="nearest")[
            ::downsample, ::downsample
        ]
        valid_detection = ndi.maximum_filter(valid_full, size=downsample, mode="nearest")[
            ::downsample, ::downsample
        ]
    else:
        work_detection = work
        valid_detection = valid_full
    best_response = np.full(work_detection.shape, -np.inf, dtype=np.float32)
    best_scale_index = np.zeros(work_detection.shape, dtype=np.uint8)

    for scale_index, sigma_px in enumerate(sigmas_px):
        response = -(sigma_px**2) * ndi.gaussian_laplace(
            work_detection,
            sigma=max(float(sigma_px), 0.65),
            mode="nearest",
        )
        update = response > best_response
        best_response[update] = response[update]
        best_scale_index[update] = scale_index

    best_response[~valid_detection] = -np.inf
    minimum_distance_px = max(
        1,
        int(round(float(detect_cfg["minimum_center_distance_um"]) / detection_spacing_um)),
    )
    coordinates = peak_local_max(
        best_response,
        min_distance=minimum_distance_px,
        threshold_abs=float(detect_cfg["log_response_threshold"]),
        exclude_border=False,
    )
    if coordinates.size == 0:
        return pd.DataFrame(
            columns=[
                "y_px",
                "x_px",
                "y_um",
                "x_um",
                "peak_snr",
                "log_response",
                "diameter_um",
                "hessian_roundness",
                "radial_contrast_snr",
                "foreground_component_label",
                "component_area_um2",
                "core_component_label",
                "core_component_area_um2",
                "core_component_integrated_snr",
            ]
        )

    # Estimate curvature after Gaussian smoothing at the selected LoG scale.
    # Finite differences on the unsmoothed SNR image are dominated by pixel
    # noise and falsely label tissue edges as round; scale-matched curvature is
    # a much more stable classical blob-shape criterion.
    coordinate_roundness = np.zeros(len(coordinates), dtype=np.float32)
    coordinate_scales = best_scale_index[coordinates[:, 0], coordinates[:, 1]]
    for scale_index, sigma_px in enumerate(sigmas_px):
        selected = np.flatnonzero(coordinate_scales == scale_index)
        if selected.size == 0:
            continue
        smoothed_scale = ndi.gaussian_filter(
            work_detection,
            sigma=max(float(sigma_px), 0.65),
            mode="nearest",
        )
        ys = coordinates[selected, 0]
        xs = coordinates[selected, 1]
        interior = (
            (ys > 0)
            & (ys < work_detection.shape[0] - 1)
            & (xs > 0)
            & (xs < work_detection.shape[1] - 1)
        )
        if not interior.any():
            continue
        valid_selected = selected[interior]
        ys = ys[interior]
        xs = xs[interior]
        center = smoothed_scale[ys, xs]
        dxx = smoothed_scale[ys, xs - 1] + smoothed_scale[ys, xs + 1] - 2.0 * center
        dyy = smoothed_scale[ys - 1, xs] + smoothed_scale[ys + 1, xs] - 2.0 * center
        dxy = (
            smoothed_scale[ys + 1, xs + 1]
            - smoothed_scale[ys + 1, xs - 1]
            - smoothed_scale[ys - 1, xs + 1]
            + smoothed_scale[ys - 1, xs - 1]
        ) / 4.0
        trace = dxx + dyy
        discriminant = np.sqrt(np.maximum((dxx - dyy) ** 2 + 4.0 * dxy**2, 0.0))
        eigen_more_negative = 0.5 * (trace - discriminant)
        eigen_less_negative = 0.5 * (trace + discriminant)
        valid_peak = (eigen_more_negative < 0.0) & (eigen_less_negative < 0.0)
        roundness = np.zeros_like(trace, dtype=np.float32)
        roundness[valid_peak] = np.abs(eigen_less_negative[valid_peak]) / np.maximum(
            np.abs(eigen_more_negative[valid_peak]),
            1e-6,
        )
        coordinate_roundness[valid_selected] = roundness

    foreground_labels = label(foreground, connectivity=2)
    component_sizes = np.bincount(foreground_labels.ravel())
    core_mask = (snr >= float(config["segmentation"]["snr_high"])) & roi
    core_labels = label(core_mask, connectivity=2)
    core_component_sizes = np.bincount(core_labels.ravel())
    core_component_integrated_snr = np.bincount(
        core_labels.ravel(),
        weights=np.clip(snr, 0.0, None).ravel(),
    )
    records: list[dict[str, float]] = []
    used_full_coordinates: set[tuple[int, int]] = set()
    for coordinate_index, (y_detection, x_detection) in enumerate(coordinates):
        approximate_y = int(y_detection * downsample)
        approximate_x = int(x_detection * downsample)
        y0 = max(0, approximate_y - downsample)
        y1 = min(snr.shape[0], approximate_y + downsample + 1)
        x0 = max(0, approximate_x - downsample)
        x1 = min(snr.shape[1], approximate_x + downsample + 1)
        local_score = np.where(valid_full[y0:y1, x0:x1], snr[y0:y1, x0:x1], -np.inf)
        if not np.isfinite(local_score).any():
            continue
        local_y, local_x = np.unravel_index(int(np.nanargmax(local_score)), local_score.shape)
        y_px = int(y0 + local_y)
        x_px = int(x0 + local_x)
        if (y_px, x_px) in used_full_coordinates:
            continue
        used_full_coordinates.add((y_px, x_px))
        peak_snr = float(snr[y_px, x_px])
        if peak_snr < float(detect_cfg["minimum_peak_snr"]):
            continue
        scale_index = int(best_scale_index[y_detection, x_detection])
        hessian_roundness = float(coordinate_roundness[coordinate_index])
        if hessian_roundness < float(detect_cfg["minimum_hessian_roundness"]):
            continue

        blob_radius_px = max(
            1.0,
            float(diameters_um[scale_index]) / (2.0 * xy_spacing_um),
        )
        outer_radius_px = max(2, int(math.ceil(1.6 * blob_radius_px)))
        local_y0 = max(0, y_px - outer_radius_px)
        local_y1 = min(snr.shape[0], y_px + outer_radius_px + 1)
        local_x0 = max(0, x_px - outer_radius_px)
        local_x1 = min(snr.shape[1], x_px + outer_radius_px + 1)
        yy, xx = np.ogrid[local_y0:local_y1, local_x0:local_x1]
        radius = np.sqrt((yy - y_px) ** 2 + (xx - x_px) ** 2)
        local_roi = roi[local_y0:local_y1, local_x0:local_x1]
        core_region = (radius <= max(1.0, 0.55 * blob_radius_px)) & local_roi
        ring_region = (
            (radius >= blob_radius_px)
            & (radius <= 1.6 * blob_radius_px)
            & local_roi
        )
        if not core_region.any() or not ring_region.any():
            continue
        local_snr = snr[local_y0:local_y1, local_x0:local_x1]
        radial_contrast = float(
            np.mean(local_snr[core_region]) - np.median(local_snr[ring_region])
        )
        if radial_contrast < float(detect_cfg.get("minimum_radial_contrast_snr", 0.0)):
            continue

        core_component_label = int(core_labels[y_px, x_px])
        core_component_area_px = (
            int(core_component_sizes[core_component_label]) if core_component_label else 0
        )
        core_component_area_um2 = float(core_component_area_px * xy_spacing_um**2)
        maximum_core_area = float(detect_cfg.get("maximum_core_component_area_um2", np.inf))
        if core_component_area_um2 > maximum_core_area:
            continue
        component_label = int(foreground_labels[y_px, x_px])
        component_area_px = int(component_sizes[component_label]) if component_label else 0
        records.append(
            {
                "y_px": int(y_px),
                "x_px": int(x_px),
                "y_um": float(y_px * xy_spacing_um),
                "x_um": float(x_px * xy_spacing_um),
                "peak_snr": peak_snr,
                "log_response": float(best_response[y_detection, x_detection]),
                "diameter_um": float(diameters_um[scale_index]),
                "hessian_roundness": float(hessian_roundness),
                "radial_contrast_snr": radial_contrast,
                "foreground_component_label": component_label,
                "component_area_um2": float(component_area_px * xy_spacing_um**2),
                "core_component_label": core_component_label,
                "core_component_area_um2": core_component_area_um2,
                "core_component_integrated_snr": float(
                    core_component_integrated_snr[core_component_label]
                ),
            }
        )
    return pd.DataFrame.from_records(records)


def process_plane(
    raw: np.ndarray,
    roi: np.ndarray,
    spacing_zyx_um: tuple[float, float, float],
    config: dict[str, Any],
) -> PlaneResult:
    """Process one full-resolution XY plane using fixed physical-scale parameters."""
    preprocess_cfg = config["preprocess"]
    segmentation_cfg = config["segmentation"]
    xy_spacing = float((spacing_zyx_um[1] + spacing_zyx_um[2]) / 2.0)
    raw_float = np.asarray(raw, dtype=np.float32)
    smooth_sigma_px = float(preprocess_cfg["smooth_sigma_xy_um"]) / xy_spacing
    background_sigma_px = float(preprocess_cfg["background_sigma_xy_um"]) / xy_spacing
    if not roi.any():
        empty_float = np.zeros_like(raw_float, dtype=np.float32)
        empty_bool = np.zeros_like(roi, dtype=bool)
        return PlaneResult(
            raw=raw_float,
            background=empty_float,
            residual=empty_float.copy(),
            snr=empty_float.copy(),
            roi=roi,
            foreground=empty_bool,
            candidates=pd.DataFrame(
                columns=[
                    "y_px",
                    "x_px",
                    "y_um",
                    "x_um",
                    "peak_snr",
                    "log_response",
                    "diameter_um",
                    "hessian_roundness",
                    "radial_contrast_snr",
                    "foreground_component_label",
                    "component_area_um2",
                    "core_component_label",
                    "core_component_area_um2",
                    "core_component_integrated_snr",
                ]
            ),
            noise_sigma=1.0,
        )

    y_indices, x_indices = np.nonzero(roi)
    processing_margin_px = int(math.ceil(3.0 * background_sigma_px))
    y0 = max(0, int(y_indices.min()) - processing_margin_px)
    y1 = min(raw_float.shape[0], int(y_indices.max()) + processing_margin_px + 1)
    x0 = max(0, int(x_indices.min()) - processing_margin_px)
    x1 = min(raw_float.shape[1], int(x_indices.max()) + processing_margin_px + 1)
    raw_crop = raw_float[y0:y1, x0:x1]
    roi_crop = roi[y0:y1, x0:x1]
    smooth_crop = ndi.gaussian_filter(raw_crop, sigma=smooth_sigma_px, mode="nearest")
    background_crop = ndi.gaussian_filter(smooth_crop, sigma=background_sigma_px, mode="nearest")
    residual_crop = smooth_crop - background_crop
    stride = max(1, int(preprocess_cfg["noise_sample_stride"]))
    sampled = residual_crop[::stride, ::stride]
    sampled_roi = roi_crop[::stride, ::stride]
    if sampled_roi.any():
        noise_sigma = _robust_sigma(sampled[sampled_roi])
    else:
        noise_sigma = _robust_sigma(sampled)
    snr_crop = residual_crop / noise_sigma
    foreground_crop = apply_hysteresis_threshold(
        snr_crop,
        low=float(segmentation_cfg["snr_low"]),
        high=float(segmentation_cfg["snr_high"]),
    )
    foreground_crop &= roi_crop
    minimum_area_px = max(
        1,
        int(round(float(segmentation_cfg["minimum_area_um2"]) / (xy_spacing**2))),
    )
    foreground_crop = remove_small_objects(foreground_crop, min_size=minimum_area_px)
    candidates = _detect_blobs(snr_crop, foreground_crop, roi_crop, xy_spacing, config)
    if not candidates.empty:
        candidates["y_px"] = candidates["y_px"].astype(int) + y0
        candidates["x_px"] = candidates["x_px"].astype(int) + x0
        candidates["y_um"] = candidates["y_px"] * xy_spacing
        candidates["x_um"] = candidates["x_px"] * xy_spacing

    background = np.zeros_like(raw_float, dtype=np.float32)
    residual = np.zeros_like(raw_float, dtype=np.float32)
    snr = np.zeros_like(raw_float, dtype=np.float32)
    foreground = np.zeros_like(roi, dtype=bool)
    background[y0:y1, x0:x1] = background_crop
    residual[y0:y1, x0:x1] = residual_crop
    snr[y0:y1, x0:x1] = snr_crop
    foreground[y0:y1, x0:x1] = foreground_crop
    return PlaneResult(
        raw=raw_float,
        background=background,
        residual=residual,
        snr=snr,
        roi=roi,
        foreground=foreground,
        candidates=candidates,
        noise_sigma=noise_sigma,
    )


def _link_candidates(candidates: pd.DataFrame, config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    if candidates.empty:
        empty_objects = pd.DataFrame(
            columns=[
                "object_id",
                "z_start",
                "z_end",
                "n_planes",
                "x_um",
                "y_um",
                "z_um",
                "peak_snr",
                "log_response",
                "diameter_xy_um",
                "equivalent_sphere_volume_um3",
                "hessian_roundness",
                "radial_contrast_snr",
                "max_core_component_area_um2",
                "single_plane",
                "accepted",
            ]
        )
        return candidates.assign(track_id=pd.Series(dtype=int)), empty_objects

    detect_cfg = config["detection"]
    link_distance = float(detect_cfg["link_distance_xy_um"])
    maximum_scale_ratio = float(detect_cfg["maximum_scale_ratio"])
    candidates = candidates.sort_values(["z_index", "candidate_in_plane"]).reset_index(drop=True)
    candidates["track_id"] = -1
    next_track_id = 1
    previous_indices: list[int] = []
    previous_z: int | None = None

    for z_index, current_group in candidates.groupby("z_index", sort=True):
        current_indices = current_group.index.to_list()
        if previous_z is not None and int(z_index) == previous_z + 1 and previous_indices and current_indices:
            prev_xy = candidates.loc[previous_indices, ["y_um", "x_um"]].to_numpy(float)
            curr_xy = candidates.loc[current_indices, ["y_um", "x_um"]].to_numpy(float)
            distances = np.linalg.norm(prev_xy[:, None, :] - curr_xy[None, :, :], axis=2)
            prev_scale = candidates.loc[previous_indices, "diameter_um"].to_numpy(float)
            curr_scale = candidates.loc[current_indices, "diameter_um"].to_numpy(float)
            ratios = np.maximum(prev_scale[:, None] / curr_scale[None, :], curr_scale[None, :] / prev_scale[:, None])
            cost = distances.copy()
            invalid = (distances > link_distance) | (ratios > maximum_scale_ratio)
            cost[invalid] = 1e6
            rows, columns = linear_sum_assignment(cost)
            matched_current: set[int] = set()
            for row, column in zip(rows, columns):
                if cost[row, column] >= 1e6:
                    continue
                prev_index = previous_indices[int(row)]
                curr_index = current_indices[int(column)]
                candidates.loc[curr_index, "track_id"] = int(candidates.loc[prev_index, "track_id"])
                matched_current.add(curr_index)
            for curr_index in current_indices:
                if curr_index not in matched_current:
                    candidates.loc[curr_index, "track_id"] = next_track_id
                    next_track_id += 1
        else:
            for curr_index in current_indices:
                candidates.loc[curr_index, "track_id"] = next_track_id
                next_track_id += 1
        previous_indices = current_indices
        previous_z = int(z_index)

    object_records: list[dict[str, Any]] = []
    single_plane_min_snr = float(detect_cfg["single_plane_minimum_snr"])
    maximum_planes = max(1, int(detect_cfg.get("maximum_planes_per_object", 999999)))
    for track_id, group in candidates.groupby("track_id", sort=True):
        peak_index = group["log_response"].idxmax()
        peak = candidates.loc[peak_index]
        weights = group["log_response"].to_numpy(float)
        weights = np.maximum(weights, 1e-6)
        n_planes = int(group["z_index"].nunique())
        peak_snr = float(group["peak_snr"].max())
        accepted = (n_planes > 1 or peak_snr >= single_plane_min_snr) and n_planes <= maximum_planes
        object_records.append(
            {
                "object_id": int(track_id),
                "z_start": int(group["z_index"].min()),
                "z_end": int(group["z_index"].max()),
                "n_planes": n_planes,
                "x_um": float(np.average(group["x_um"], weights=weights)),
                "y_um": float(np.average(group["y_um"], weights=weights)),
                "z_um": float(peak["z_um"]),
                "peak_snr": peak_snr,
                "log_response": float(peak["log_response"]),
                "diameter_xy_um": float(peak["diameter_um"]),
                "equivalent_sphere_volume_um3": float(
                    math.pi / 6.0 * float(peak["diameter_um"]) ** 3
                ),
                "hessian_roundness": float(peak["hessian_roundness"]),
                "radial_contrast_snr": float(peak["radial_contrast_snr"]),
                "max_component_area_um2": float(group["component_area_um2"].max()),
                "max_core_component_area_um2": float(
                    group["core_component_area_um2"].max()
                ),
                "single_plane": bool(n_planes == 1),
                "accepted": bool(accepted),
            }
        )
    objects = pd.DataFrame.from_records(object_records)
    return candidates, objects


def _choose_qc_slices(active_z: np.ndarray, count: int) -> set[int]:
    if active_z.size == 0 or count <= 0:
        return set()
    positions = np.linspace(0, active_z.size - 1, min(count, active_z.size))
    return {int(active_z[int(round(position))]) for position in positions}


def _save_projection_qc(
    output_path: Path,
    raw_mip: np.ndarray,
    snr_mip: np.ndarray,
    mask_mip: np.ndarray,
    roi_mip: np.ndarray,
    objects: pd.DataFrame,
    spacing_zyx: tuple[float, float, float],
    title: str,
) -> None:
    raw_lo, raw_hi = np.percentile(raw_mip[roi_mip] if roi_mip.any() else raw_mip, [1, 99.8])
    snr_hi = max(5.0, float(np.percentile(snr_mip[roi_mip] if roi_mip.any() else snr_mip, 99.8)))
    fig, axes = plt.subplots(2, 2, figsize=(13, 12), constrained_layout=True)
    axes[0, 0].imshow(raw_mip, cmap="gray", vmin=raw_lo, vmax=raw_hi)
    axes[0, 0].contour(roi_mip, levels=[0.5], colors="cyan", linewidths=0.5)
    axes[0, 0].set_title("Raw max projection + ROI")
    axes[0, 1].imshow(snr_mip, cmap="magma", vmin=0, vmax=snr_hi)
    axes[0, 1].set_title("Maximum local SNR")
    axes[1, 0].imshow(raw_mip, cmap="gray", vmin=raw_lo, vmax=raw_hi)
    axes[1, 0].imshow(mask_mip, cmap="spring", alpha=np.where(mask_mip, 0.55, 0.0))
    axes[1, 0].set_title("Broad hysteresis candidate-mask projection")
    axes[1, 1].imshow(raw_mip, cmap="gray", vmin=raw_lo, vmax=raw_hi)
    if not objects.empty:
        accepted = objects[objects["accepted"]]
        multi = accepted[~accepted["single_plane"]]
        single = accepted[accepted["single_plane"]]
        axes[1, 1].scatter(
            multi["x_um"] / spacing_zyx[2],
            multi["y_um"] / spacing_zyx[1],
            s=7,
            facecolors="none",
            edgecolors="lime",
            linewidths=0.5,
            label=f"multi-z ({len(multi)})",
        )
        axes[1, 1].scatter(
            single["x_um"] / spacing_zyx[2],
            single["y_um"] / spacing_zyx[1],
            s=5,
            c="orange",
            linewidths=0,
            label=f"single-z ({len(single)})",
        )
        axes[1, 1].legend(loc="lower right", fontsize=8)
    axes[1, 1].set_title("Accepted LD centers")
    for axis in axes.flat:
        axis.axis("off")
    fig.suptitle(title)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _save_full_resolution_detection_qc(
    output_path: Path,
    raw_mip: np.ndarray,
    roi_mip: np.ndarray,
    objects: pd.DataFrame,
    spacing_zyx: tuple[float, float, float],
    marker_radius_px: int,
) -> None:
    """Save an exact-size RGB projection with every accepted center marked."""
    from PIL import Image

    display_values = raw_mip[roi_mip] if roi_mip.any() else raw_mip.ravel()
    lo, hi = np.percentile(display_values, [0.5, 99.9])
    normalized = np.clip((raw_mip - lo) / max(hi - lo, 1e-6), 0.0, 1.0)
    gray = np.round(normalized * 255.0).astype(np.uint8)
    rgb = np.repeat(gray[:, :, None], 3, axis=2)

    roi_boundary = find_boundaries(roi_mip, mode="outer")
    roi_boundary = binary_dilation(roi_boundary, footprint=disk(1))
    rgb[roi_boundary] = np.array([0, 255, 255], dtype=np.uint8)

    radius = max(2, int(marker_radius_px))
    if not objects.empty:
        accepted = objects[objects["accepted"]]
        # Draw single-plane candidates first, then multi-z objects on top.
        for single_plane, color in (
            (True, np.array([255, 165, 0], dtype=np.uint8)),
            (False, np.array([0, 255, 0], dtype=np.uint8)),
        ):
            subset = accepted[accepted["single_plane"] == single_plane]
            for row in subset.itertuples(index=False):
                center_y = int(round(float(row.y_um) / spacing_zyx[1]))
                center_x = int(round(float(row.x_um) / spacing_zyx[2]))
                for current_radius in (radius, radius + 1):
                    rr, cc = circle_perimeter(
                        center_y,
                        center_x,
                        current_radius,
                        method="andres",
                        shape=gray.shape,
                    )
                    rgb[rr, cc] = color
    Image.fromarray(rgb, mode="RGB").save(output_path)


def _save_3d_detection_marker_tiff(
    output_path: Path,
    shape_zyx: tuple[int, int, int],
    linked_candidates: pd.DataFrame,
    objects: pd.DataFrame,
    marker_diameter_px: int,
) -> None:
    """Save an exact-ZYX grayscale stack containing accepted detections only.

    Every accepted track is marked on each z plane where that object produced
    a linked 2-D detection.  Markers are filled circles with value 255 on a
    black background, and the stack contains no raw-image intensities.  Planes
    are yielded one at a time to avoid allocating the full volume.
    """
    import tifffile

    depth, height, width = (int(value) for value in shape_zyx)
    diameter = max(1, int(marker_diameter_px))
    if diameter % 2 == 0:
        raise ValueError("output.marker_3d_diameter_px must be an odd integer")
    half_width = diameter // 2
    offsets = np.arange(-half_width, half_width + 1, dtype=int)
    offset_y, offset_x = np.meshgrid(offsets, offsets, indexing="ij")
    filled_circle = offset_y**2 + offset_x**2 < (diameter / 2.0) ** 2
    circle_dy = offset_y[filled_circle]
    circle_dx = offset_x[filled_circle]

    if linked_candidates.empty or objects.empty:
        accepted_detections = linked_candidates.iloc[0:0]
    else:
        accepted_track_ids = set(
            objects.loc[objects["accepted"], "object_id"].astype(int).tolist()
        )
        accepted_detections = linked_candidates[
            linked_candidates["track_id"].astype(int).isin(accepted_track_ids)
        ]
    detections_by_z = {
        int(z_index): group[["y_px", "x_px"]].to_numpy(dtype=int)
        for z_index, group in accepted_detections.groupby("z_index", sort=False)
    }

    def marker_planes():
        for z_index in range(depth):
            plane = np.zeros((height, width), dtype=np.uint8)
            for y_px, x_px in detections_by_z.get(z_index, ()):
                marker_y = int(y_px) + circle_dy
                marker_x = int(x_px) + circle_dx
                valid = (
                    (marker_y >= 0)
                    & (marker_y < height)
                    & (marker_x >= 0)
                    & (marker_x < width)
                )
                plane[marker_y[valid], marker_x[valid]] = 255
            yield plane

    tifffile.imwrite(
        output_path,
        data=marker_planes(),
        shape=(depth, height, width),
        dtype=np.uint8,
        photometric="minisblack",
        compression="deflate",
        # Prefer classic TIFF for direct ImageJ compatibility when the full
        # uncompressed stack is below the 4-GB TIFF limit.  Larger stacks
        # transparently switch to BigTIFF and can be opened through Fiji/Bio-Formats.
        bigtiff=(depth * height * width >= 4_000_000_000),
        metadata={"axes": "ZYX"},
        rowsperstrip=256,
    )


def _save_depth_qc(output_path: Path, planes: pd.DataFrame, title: str) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True, constrained_layout=True)
    axes[0].plot(planes["z_um"], planes["raw_median"], label="raw median")
    axes[0].plot(planes["z_um"], planes["raw_p99"], label="raw p99")
    axes[0].plot(planes["z_um"], planes["noise_sigma"], label="residual noise σ")
    axes[0].set_ylabel("Intensity / residual a.u.")
    axes[0].legend()
    axes[1].plot(planes["z_um"], planes["candidate_count"], label="2-D candidates")
    axes[1].plot(planes["z_um"], planes["foreground_area_um2"], label="foreground area µm²")
    axes[1].plot(planes["z_um"], planes["roi_area_um2"], label="ROI area µm²", alpha=0.6)
    axes[1].set_xlabel("z (µm)")
    axes[1].set_ylabel("Count / area")
    axes[1].legend()
    fig.suptitle(title)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _save_slice_qc(output_path: Path, stored: list[dict[str, Any]], title: str) -> None:
    if not stored:
        return
    fig, axes = plt.subplots(len(stored), 3, figsize=(12, 3.8 * len(stored)), constrained_layout=True)
    if len(stored) == 1:
        axes = np.asarray([axes])
    for row, item in enumerate(stored):
        raw = item["raw"]
        snr = item["snr"]
        mask = item["mask"]
        roi = item["roi"]
        candidates = item["candidates"]
        raw_values = raw[roi] if roi.any() else raw
        lo, hi = np.percentile(raw_values, [1, 99.8])
        axes[row, 0].imshow(raw, cmap="gray", vmin=lo, vmax=hi)
        axes[row, 0].contour(roi, levels=[0.5], colors="cyan", linewidths=0.5)
        axes[row, 0].set_title(f"z={item['z_index']} raw + ROI")
        axes[row, 1].imshow(snr, cmap="magma", vmin=0, vmax=max(5, np.percentile(snr[roi], 99.8) if roi.any() else 5))
        axes[row, 1].set_title("local SNR")
        axes[row, 2].imshow(raw, cmap="gray", vmin=lo, vmax=hi)
        axes[row, 2].imshow(mask, cmap="spring", alpha=np.where(mask, 0.5, 0.0))
        if not candidates.empty:
            factor = int(item["downsample"])
            axes[row, 2].scatter(
                candidates["x_px"] / factor,
                candidates["y_px"] / factor,
                s=18,
                facecolors="none",
                edgecolors="lime",
                linewidths=0.7,
            )
        axes[row, 2].set_title(f"mask + LoG centers (n={len(candidates)})")
        for column in range(3):
            axes[row, column].axis("off")
    fig.suptitle(title)
    fig.savefig(output_path, dpi=170)
    plt.close(fig)


def analyze_ims(
    ims_path: str | Path,
    config: dict[str, Any],
    output_root: str | Path = "results",
    run_name: str | None = None,
) -> Path:
    """Run the complete fixed pipeline and return the output directory."""
    ims_path = Path(ims_path)
    if not ims_path.is_file():
        raise FileNotFoundError(ims_path)
    if not h5py.is_hdf5(ims_path):
        raise ValueError(f"Input is not an HDF5/IMS file: {ims_path}")
    from .anatomical_v2 import resolve_v2_config
    config = resolve_v2_config(config)
    config = resolve_sample_guide(config, ims_path)
    output_root = Path(output_root)
    config_hash = _config_hash(config)
    if run_name is None:
        run_name = f"{ims_path.stem}_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{config_hash}"
    output_dir = output_root / run_name
    output_dir.mkdir(parents=True, exist_ok=False)
    save_qc = bool(config["output"].get("save_qc", True))
    save_3d_marker = bool(config["output"].get("save_3d_marker_tiff", True))
    print_plane_progress = bool(
        config["output"].get("print_plane_progress", True)
    )
    qc_dir: Path | None = None
    if save_qc or save_3d_marker:
        qc_dir = output_dir / "qc"
        qc_dir.mkdir()

    with (output_dir / "config_used.yaml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump(config, handle, sort_keys=False, allow_unicode=True)

    reader = ims(str(ims_path), ResolutionLevelLock=int(config["input"]["resolution_level"]), squeeze_output=True)
    mask_writer = None
    try:
        channel = int(config["input"]["channel_index"])
        time_index = int(config["input"]["time_index"])
        if channel >= reader.Channels:
            raise IndexError(f"channel_index={channel}, but file has {reader.Channels} channel(s)")
        if time_index >= reader.TimePoints:
            raise IndexError(f"time_index={time_index}, but file has {reader.TimePoints} time point(s)")
        spacing = tuple(float(x) for x in reader.resolution)
        shape_zyx = tuple(int(x) for x in reader.shape[-3:])
        save_roi_qc = bool(config["output"].get("save_roi_qc", False))
        roi_provider, roi_metadata = _make_roi_provider(reader, config, keep_roi_stages=save_roi_qc)
        if config["roi"].get("central_method") == "anatomical_v2":
            with (output_dir / "roi_diagnostics.json").open("w", encoding="utf-8") as handle:
                json.dump({**roi_metadata, "input_file": str(ims_path.resolve()), "config_hash": config_hash}, handle, ensure_ascii=False, indent=2)
            if save_roi_qc and roi_provider.roi_artifacts is not None:
                from .roi_qc import save_structural_roi_qc
                pack = roi_provider.roi_artifacts
                save_structural_roi_qc(output_dir / "roi_qc", pack["volume"], roi_provider.low_mask,
                                       roi_provider.low_spacing, roi_metadata, pack["stages"], title=ims_path.name)
                roi_provider.roi_artifacts = None
                del pack
            if roi_metadata.get("roi_status") == "failed" or (roi_metadata.get("roi_status") == "needs_review" and config["roi"].get("v2_review_policy") == "error"):
                raise ValueError(f"ROI {roi_metadata['roi_status']}; see {output_dir / 'roi_diagnostics.json'}")
        active_z = roi_provider.active_z_indices()
        qc_z = (
            _choose_qc_slices(active_z, int(config["output"]["qc_slices"]))
            if save_qc
            else set()
        )
        downsample = max(1, int(config["output"]["qc_downsample"]))

        raw_mip = np.zeros(shape_zyx[1:], dtype=np.float32) if save_qc else None
        snr_mip = (
            np.full(shape_zyx[1:], -np.inf, dtype=np.float32)
            if save_qc
            else None
        )
        mask_mip = np.zeros(shape_zyx[1:], dtype=bool) if save_qc else None
        roi_mip = np.zeros(shape_zyx[1:], dtype=bool) if save_qc else None
        plane_records: list[dict[str, Any]] = []
        candidate_frames: list[pd.DataFrame] = []
        stored_slices: list[dict[str, Any]] = []
        roi_voxels = 0
        foreground_voxels = 0
        integrated_corrected_intensity = 0.0
        saturation_voxels = 0
        total_voxels = 0

        if bool(config["output"]["save_mask_tiff"]):
            import tifffile

            mask_writer = tifffile.TiffWriter(output_dir / "ld_positive_mask.tif", bigtiff=True)

        for z_index in range(shape_zyx[0]):
            raw = np.asarray(reader[time_index, channel, z_index, :, :])
            roi = roi_provider.plane(z_index)
            result = process_plane(raw, roi, spacing, config)
            candidates = result.candidates.copy()
            if not candidates.empty:
                candidates.insert(0, "candidate_in_plane", np.arange(len(candidates), dtype=int))
                candidates.insert(0, "z_um", z_index * spacing[0])
                candidates.insert(0, "z_index", z_index)
                candidate_frames.append(candidates)

            if save_qc:
                assert raw_mip is not None
                assert snr_mip is not None
                assert mask_mip is not None
                assert roi_mip is not None
                raw_mip = np.maximum(raw_mip, result.raw)
                snr_mip = np.maximum(snr_mip, result.snr)
                mask_mip |= result.foreground
                roi_mip |= result.roi
            roi_count = int(result.roi.sum())
            foreground_count = int(result.foreground.sum())
            roi_voxels += roi_count
            foreground_voxels += foreground_count
            integrated_corrected_intensity += float(result.residual[result.foreground].sum())
            dtype_max = np.iinfo(raw.dtype).max if np.issubdtype(raw.dtype, np.integer) else np.inf
            saturation_voxels += int(np.count_nonzero(raw >= dtype_max))
            total_voxels += int(raw.size)
            roi_values = result.raw[result.roi]
            plane_records.append(
                {
                    "z_index": z_index,
                    "z_um": z_index * spacing[0],
                    "roi_area_um2": roi_count * spacing[1] * spacing[2],
                    "foreground_area_um2": foreground_count * spacing[1] * spacing[2],
                    "raw_median": float(np.median(roi_values)) if roi_values.size else np.nan,
                    "raw_p99": float(np.percentile(roi_values, 99)) if roi_values.size else np.nan,
                    "noise_sigma": float(result.noise_sigma),
                    "candidate_count": int(len(candidates)),
                }
            )

            if mask_writer is not None:
                mask_writer.write(
                    result.foreground.astype(np.uint8),
                    compression=str(config["output"]["mask_compression"]),
                    photometric="minisblack",
                )

            if z_index in qc_z:
                stored_slices.append(
                    {
                        "z_index": z_index,
                        "downsample": downsample,
                        "raw": result.raw[::downsample, ::downsample].copy(),
                        "snr": result.snr[::downsample, ::downsample].copy(),
                        "mask": result.foreground[::downsample, ::downsample].copy(),
                        "roi": result.roi[::downsample, ::downsample].copy(),
                        "candidates": candidates,
                    }
                )
            if print_plane_progress:
                print(
                    f"z {z_index + 1:3d}/{shape_zyx[0]} | ROI {roi_count:8d} px | "
                    f"mask {foreground_count:7d} px | candidates {len(candidates):4d}",
                    flush=True,
                )
        if bool(config["output"].get("save_roi_mask_tiff", False)):
            _save_roi_mask_tiff(
                output_dir / "analysis_roi_used_3d.tif", roi_provider, spacing,
                compression=str(config["output"].get("mask_compression", "deflate")),
            )
    finally:
        if mask_writer is not None:
            mask_writer.close()
        reader.close()

    if candidate_frames:
        candidate_table = pd.concat(candidate_frames, ignore_index=True)
    else:
        candidate_table = pd.DataFrame(columns=["z_index", "y_px", "x_px", "track_id"])
    linked_candidates, object_table = _link_candidates(candidate_table, config)
    accepted_objects = object_table[object_table["accepted"]] if not object_table.empty else object_table
    if not linked_candidates.empty and not accepted_objects.empty:
        accepted_track_ids = set(accepted_objects["object_id"].astype(int).tolist())
        accepted_components = linked_candidates[
            linked_candidates["track_id"].astype(int).isin(accepted_track_ids)
            & (linked_candidates["core_component_label"] > 0)
        ].drop_duplicates(["z_index", "core_component_label"])
    else:
        accepted_components = pd.DataFrame()
    plane_table = pd.DataFrame.from_records(plane_records)
    voxel_volume = spacing[0] * spacing[1] * spacing[2]
    roi_volume = roi_voxels * voxel_volume
    foreground_volume = foreground_voxels * voxel_volume
    accepted_core_volume = (
        float(accepted_components["core_component_area_um2"].sum()) * spacing[0]
        if not accepted_components.empty
        else 0.0
    )
    accepted_core_integrated_snr = (
        float(accepted_components["core_component_integrated_snr"].sum())
        if not accepted_components.empty
        else 0.0
    )
    equivalent_sphere_volume = (
        float(accepted_objects["equivalent_sphere_volume_um3"].sum())
        if not accepted_objects.empty
        else 0.0
    )
    brain_record = {
        "file": ims_path.name,
        "brain_id": ims_path.stem,
        "time_index": int(config["input"]["time_index"]),
        "channel_index": int(config["input"]["channel_index"]),
        "shape_z": shape_zyx[0],
        "shape_y": shape_zyx[1],
        "shape_x": shape_zyx[2],
        "spacing_z_um": spacing[0],
        "spacing_y_um": spacing[1],
        "spacing_x_um": spacing[2],
        "roi_volume_um3": roi_volume,
        "ld_count": int(len(accepted_objects)),
        "ld_count_density_per_1e5_um3": float(len(accepted_objects) / roi_volume * 1e5) if roi_volume else np.nan,
        "ld_positive_volume_um3": accepted_core_volume,
        "ld_volume_fraction": float(accepted_core_volume / roi_volume) if roi_volume else np.nan,
        "ld_core_integrated_snr_voxels": accepted_core_integrated_snr,
        "ld_equivalent_sphere_volume_um3": equivalent_sphere_volume,
        "broad_hysteresis_mask_volume_um3": foreground_volume,
        "broad_hysteresis_mask_fraction": float(foreground_volume / roi_volume) if roi_volume else np.nan,
        "broad_hysteresis_integrated_corrected_intensity_au": integrated_corrected_intensity,
        "single_plane_object_fraction": float(accepted_objects["single_plane"].mean()) if len(accepted_objects) else np.nan,
        "saturation_fraction": float(saturation_voxels / total_voxels) if total_voxels else np.nan,
        "config_hash": config_hash,
    }
    if config["roi"].get("central_method") == "anatomical_v2":
        brain_record.update(roi_method=roi_metadata["method"], roi_status=roi_metadata["roi_status"],
                            roi_scope=roi_metadata['roi_scope'],
                            exclude_optic_lobes=roi_metadata['exclude_optic_lobes'],
                            roi_algorithm_version=roi_metadata["algorithm_version"],
                            roi_algorithm_source_hash=roi_metadata["algorithm_source_hash"],
                            roi_review_flags="; ".join(roi_metadata["review_flags"]))
    brain_table = pd.DataFrame([brain_record])
    brain_table.to_csv(output_dir / "per_brain.csv", index=False)
    plane_table.to_csv(output_dir / "per_plane.csv", index=False)
    linked_candidates.to_csv(output_dir / "per_detection.csv", index=False)
    object_table.to_csv(output_dir / "per_object.csv", index=False)

    title = f"{ims_path.name}\nconfig={config_hash}, accepted objects={len(accepted_objects)}"
    if save_qc:
        assert qc_dir is not None
        assert raw_mip is not None
        assert snr_mip is not None
        assert mask_mip is not None
        assert roi_mip is not None
        _save_projection_qc(
            qc_dir / "projection_overview.png",
            raw_mip,
            snr_mip,
            mask_mip,
            roi_mip,
            object_table,
            spacing,
            title,
        )
        _save_full_resolution_detection_qc(
            qc_dir / "detections_full_resolution.png",
            raw_mip,
            roi_mip,
            object_table,
            spacing,
            int(config["output"]["full_resolution_marker_radius_px"]),
        )
        _save_depth_qc(qc_dir / "depth_profiles.png", plane_table, title)
        _save_slice_qc(qc_dir / "slice_overlays.png", stored_slices, title)
    if save_3d_marker:
        assert qc_dir is not None
        _save_3d_detection_marker_tiff(
            qc_dir / "detections_3d_mask.tif",
            shape_zyx,
            linked_candidates,
            object_table,
            int(config["output"].get("marker_3d_diameter_px", 5)),
        )

    run_info = {
        "created_at": datetime.now().isoformat(),
        "python": sys.version,
        "input_file": str(ims_path.resolve()),
        "output_directory": str(output_dir.resolve()),
        "config_hash": config_hash,
        "roi": roi_metadata,
        "summary": brain_record,
    }
    with (output_dir / "run_info.json").open("w", encoding="utf-8") as handle:
        json.dump(run_info, handle, ensure_ascii=False, indent=2)
    return output_dir


def preview_plane(
    ims_path: str | Path,
    config: dict[str, Any],
    z_index: int,
    output_path: str | Path,
) -> Path:
    """Render a single-plane parameter preview without running the full stack."""
    ims_path = Path(ims_path)
    config = resolve_sample_guide(config, ims_path)
    output_path = Path(output_path)
    reader = ims(str(ims_path), ResolutionLevelLock=int(config["input"]["resolution_level"]), squeeze_output=True)
    try:
        channel = int(config["input"]["channel_index"])
        time_index = int(config["input"]["time_index"])
        spacing = tuple(float(x) for x in reader.resolution)
        if not 0 <= z_index < reader.shape[-3]:
            raise IndexError(f"z_index must be between 0 and {reader.shape[-3] - 1}")
        roi_provider, _ = _make_roi_provider(reader, config)
        raw = np.asarray(reader[time_index, channel, z_index, :, :])
        result = process_plane(raw, roi_provider.plane(z_index), spacing, config)
    finally:
        reader.close()

    roi_values = result.raw[result.roi] if result.roi.any() else result.raw.ravel()
    raw_lo, raw_hi = np.percentile(roi_values, [1, 99.8])
    snr_hi = max(5.0, float(np.percentile(result.snr[result.roi], 99.8)) if result.roi.any() else 5.0)
    fig, axes = plt.subplots(2, 2, figsize=(13, 12), constrained_layout=True)
    axes[0, 0].imshow(result.raw, cmap="gray", vmin=raw_lo, vmax=raw_hi)
    axes[0, 0].contour(result.roi, levels=[0.5], colors="cyan", linewidths=0.6)
    axes[0, 0].set_title("Raw + ROI")
    axes[0, 1].imshow(result.background, cmap="gray")
    axes[0, 1].set_title("Estimated background")
    axes[1, 0].imshow(result.snr, cmap="magma", vmin=0, vmax=snr_hi)
    axes[1, 0].set_title(f"Local SNR, noise σ={result.noise_sigma:.3g}")
    axes[1, 1].imshow(result.raw, cmap="gray", vmin=raw_lo, vmax=raw_hi)
    axes[1, 1].imshow(result.foreground, cmap="spring", alpha=np.where(result.foreground, 0.5, 0.0))
    if not result.candidates.empty:
        axes[1, 1].scatter(
            result.candidates["x_px"],
            result.candidates["y_px"],
            s=20,
            facecolors="none",
            edgecolors="lime",
            linewidths=0.7,
        )
    axes[1, 1].set_title(f"Mask + LoG centers (n={len(result.candidates)})")
    for axis in axes.flat:
        axis.axis("off")
    fig.suptitle(f"{ims_path.name}, z={z_index}, config={_config_hash(config)}")
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path
