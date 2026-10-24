"""Classical anatomical ROI: tissue geometry, 3-D necks, and optional 2-D guide.

Intensity is used for tissue support and weak edge costs, never as an
``optic lobe is brighter`` classifier. No lipid-droplet detections are used.
The automatic assignment of geometric lobes is a heuristic, not an atlas or
an anatomically validated classifier; diagnostic flags must be reviewed.
"""
from __future__ import annotations

import copy
import hashlib
import json
import warnings
from pathlib import Path
from typing import Any

import numpy as np
from scipy import ndimage as ndi
from skimage.draw import polygon
from skimage.filters import threshold_otsu
from skimage.morphology import binary_closing, disk, remove_small_objects
from skimage.segmentation import watershed
from skimage.transform import resize


ANATOMICAL_DEFAULTS = {
    "central_method": "anatomical",
    "anatomical_xy_um": 1.0,
    "anatomical_smooth_um": 2.0,
    "anatomical_tissue_margin_um": 2.0,
    "anatomical_minimum_area_um2": 60.0,
    "anatomical_focus_background_um": 8.0,
    "anatomical_focus_fraction": 0.08,
    "anatomical_focused_support_fraction": 0.03,
    "anatomical_peak_separation_um": 25.0,
    "anatomical_minimum_core_radius_um": 8.0,
    "anatomical_central_radius_ratio": 1.7,
    "anatomical_edge_weight": 0.12,
    "anatomical_boundary_smooth_um": 4.0,
    "anatomical_z_smooth_um": 2.0,
    "anatomical_central_margin_um": 3.0,
    "anatomical_terminal_fraction": 0.30,
    "anatomical_terminal_margin_um": 3.0,
    "anatomical_terminal_focus_fraction": 0.025,
    "anatomical_terminal_context_um": 12.0,
    "guide_mask_path": None,
    "guide_manifest_path": None,
    "guide_z_index": None,
    "guide_allow_resize": False,
    "guide_refine_band_um": 8.0,
    "guide_max_expansion_um": 20.0,
    "guide_expansion_per_z_um": 0.4,
}


class MissingGuideError(FileNotFoundError):
    """A sample has no manual guide; batch workflows may skip this sample."""


def compatible_guide_path(path: str | Path) -> Path:
    """Honor an existing path; otherwise try the same-folder ROI naming alias.

    Fiji defaults to the projection basename (sample_max_projection.roi),
    while older manifests used sample.roi. Only these .roi naming variants
    are exchanged; explicitly configured TIFF/PNG/ZIP paths are unchanged.
    """
    path = Path(path)
    if path.exists() or path.suffix.lower() != ".roi":
        return path
    suffix = "_max_projection"
    stem = path.stem
    alternate_stem = stem[:-len(suffix)] if stem.lower().endswith(suffix) else stem + suffix
    alternate = path.with_name(alternate_stem + path.suffix)
    return alternate if alternate.is_file() else path


def resolve_sample_guide(config: dict, ims_path: str | Path) -> dict:
    """Resolve a guide template or JSON mapping, and fingerprint the source.

    JSON: {"sample.ims": {"mask_path": "rois/sample.roi", "z_index": null}}.
    Use absolute IMS keys when one shared manifest covers same-named sources.
    A per-directory manifest may use basenames and the ``{source_dir}``
    template in guide_manifest_path, so the dataset remains movable.
    Relative mask paths in a manifest are relative to the manifest's directory.
    Missing .roi paths also try the same-directory name with/without
    ``_max_projection``; an existing configured path always has priority.
    ``z_index`` always refers to native level-0 z, zero based; null means MIP.
    """
    result = copy.deepcopy(config)
    roi = result["roi"]
    if roi["mode"] != "guided_central_brain":
        return result
    source = Path(ims_path).resolve()
    guide = roi.get("guide_mask_path")
    manifest_name = roi.get("guide_manifest_path")
    if manifest_name:
        manifest_path = Path(str(manifest_name).format(
            source_dir=str(source.parent), parent=source.parent.name,
            stem=source.stem,
        )).resolve()
        try:
            manifest_text = manifest_path.read_text(encoding="utf-8-sig")
        except FileNotFoundError as error:
            raise MissingGuideError(f"Guide manifest not found for {source.name}: {manifest_path}") from error
        mapping = json.loads(manifest_text)
        if not isinstance(mapping, dict):
            raise ValueError("guide_manifest_path must contain a JSON object")
        full_key = str(source).replace("\\", "/").casefold()
        exact = [v for k, v in mapping.items() if k.replace("\\", "/").casefold() == full_key]
        short = [v for k, v in mapping.items() if k.casefold() == source.name.casefold()]
        matches = exact or short
        if not matches:
            raise MissingGuideError(f"No guide entry for {source.name} in {manifest_path}")
        if len(matches) != 1:
            raise ValueError(f"Expected one unambiguous guide entry for {source}; found {len(matches)}")
        entry = matches[0]
        if not isinstance(entry, dict) or not entry.get("mask_path"):
            raise ValueError(f"Guide entry for {source.name} needs mask_path and optional z_index")
        guide_path = Path(entry["mask_path"])
        if not guide_path.is_absolute():
            guide_path = manifest_path.parent / guide_path
        guide = str(guide_path)
        roi["guide_z_index"] = entry.get("z_index", roi.get("guide_z_index"))
        roi["guide_manifest_path"] = str(manifest_path)
    elif guide:
        guide = str(guide).format(stem=source.stem, parent=source.parent.name,
                                  source_dir=str(source.parent))
    if not guide:
        raise ValueError("guided_central_brain requires guide_mask_path or guide_manifest_path")
    guide_path = compatible_guide_path(Path(guide).resolve())
    if not guide_path.exists():
        raise MissingGuideError(f"Guide not found for {source.name}: {guide_path}")
    if not guide_path.is_file():
        raise ValueError(f"Guide must be a file for {source.name}: {guide_path}")
    index = roi.get("guide_z_index")
    if index is not None and (isinstance(index, bool) or not isinstance(index, int) or index < 0):
        raise ValueError("guide_z_index must be a nonnegative integer (native zero-based z), or null")
    roi["guide_mask_path"] = str(guide_path)
    try:
        guide_bytes = guide_path.read_bytes()
    except FileNotFoundError as error:
        raise MissingGuideError(f"Guide not found for {source.name}: {guide_path}") from error
    roi["guide_sha256"] = hashlib.sha256(guide_bytes).hexdigest()
    return result


def load_guide_mask(path: str | Path, native_yx: tuple[int, int], *, allow_resize=False) -> np.ndarray:
    """Read filled grayscale TIFF/PNG or one area ImageJ .roi/.zip.

    Coordinates must use the unrotated, uncropped native IMS field of view.
    Color annotated JPEGs are deliberately NOT accepted as analysis masks.
    """
    path = Path(path)
    if path.suffix.lower() in {".roi", ".zip"}:
        from roifile import roiread
        roi = roiread(path)
        if isinstance(roi, list):
            if len(roi) != 1:
                raise ValueError("Guide ZIP must contain exactly one area ROI")
            roi = roi[0]
        if roi.roitype.name not in {"POLYGON", "RECT", "OVAL", "FREEHAND", "TRACED"}:
            raise ValueError(f"An area ROI is required, not {roi.roitype.name}")
        coords = np.asarray(roi.coordinates())
        if coords.ndim != 2 or coords.shape[0] < 3 or coords.shape[1] != 2:
            raise ValueError("Guide must have one closed area boundary (no composite ROI)")
        if (coords < 0).any() or (coords[:, 0] > native_yx[1]).any() or (coords[:, 1] > native_yx[0]).any():
            raise ValueError("ROI coordinates lie outside the native IMS XY field")
        mask = np.zeros(native_yx, bool)
        rr, cc = polygon(coords[:, 1], coords[:, 0], mask.shape)
        mask[rr, cc] = True
    else:
        if path.suffix.lower() in {".tif", ".tiff"}:
            import tifffile
            mask = np.asarray(tifffile.imread(path))
        elif path.suffix.lower() == ".png":
            from PIL import Image
            mask = np.asarray(Image.open(path))
        else:
            raise ValueError("Use a filled 2-D TIFF/PNG mask or one ImageJ .roi/.zip")
        if mask.ndim != 2 or not np.isfinite(mask).all():
            raise ValueError(f"Guide must be one finite grayscale XY plane, got {mask.shape}")
        mask = mask > 0
        if mask.shape != native_yx:
            if not allow_resize:
                raise ValueError(f"Guide XY {mask.shape} != native IMS {native_yx}; set guide_allow_resize only for the same full FOV")
            ratio = (mask.shape[0] / mask.shape[1]) / (native_yx[0] / native_yx[1])
            if not np.isclose(ratio, 1.0, rtol=0.01):
                raise ValueError("Guide aspect ratio differs from the native IMS field")
            mask = resize(mask.astype(np.uint8), native_yx, order=0, preserve_range=True, anti_aliasing=False) > 0
    if not mask.any() or mask.all():
        raise ValueError("Guide must contain a nonempty central-brain area and some exterior pixels")
    filled = ndi.binary_fill_holes(mask)
    if filled.sum() > 1.5 * mask.sum():
        raise ValueError("Guide looks like an outline: export a FILLED binary mask, or use an ImageJ .roi")
    return filled


def read_structural_volume(reader: Any, config: dict) -> tuple[np.ndarray, tuple, int]:
    """Stream an XY-downsampled volume without binning or skipping native z."""
    rc = {**ANATOMICAL_DEFAULTS, **config["roi"]}
    channel, time = int(config["input"]["channel_index"]), int(config["input"]["time_index"])
    if not 0 <= channel < reader.Channels or not 0 <= time < reader.TimePoints:
        raise IndexError("Requested channel/time is outside the IMS metadata")
    shape0 = tuple(int(x) for x in reader.metaData[0, time, channel, "shape"][-3:])
    spacing0 = tuple(float(x) for x in reader.metaData[0, time, channel, "resolution"])
    target_um = float(rc["anatomical_xy_um"])
    if target_um <= 0:
        raise ValueError("anatomical_xy_um must be positive")
    levels = [level for level in range(reader.ResolutionLevels)
              if reader.metaData[level, time, channel, "shape"][-3] == shape0[0]
              and max(reader.metaData[level, time, channel, "resolution"][1:]) <= target_um]
    level = max(levels) if levels else 0
    target_yx = tuple(max(1, int(round(n * s / max(target_um, s))))
                      for n, s in zip(shape0[1:], spacing0[1:]))
    spacing = (spacing0[0], *(n * s / target for n, s, target in zip(shape0[1:], spacing0[1:], target_yx)))
    volume = np.empty((shape0[0], *target_yx), np.float32)
    for z in range(shape0[0]):
        raw = np.asarray(reader[level, time, channel, z, :, :], dtype=np.float32)
        volume[z] = resize(raw, target_yx, order=1, preserve_range=True, anti_aliasing=True)
    return volume, spacing, level


def _sdf(mask: np.ndarray, spacing_yx) -> np.ndarray:
    # Explicit padding makes image borders real boundaries even in cropped brains.
    padded = np.pad(mask, 1, constant_values=False)
    inside = ndi.distance_transform_edt(padded, sampling=spacing_yx)[1:-1, 1:-1]
    outside = ndi.distance_transform_edt(~mask, sampling=spacing_yx)
    return (inside - outside).astype(np.float32)


def _tissue_support(volume, spacing, rc):
    sy, sx = spacing[1:]
    area = sy * sx
    smooth = ndi.gaussian_filter(volume, (0, rc["anatomical_smooth_um"] / sy,
                                          rc["anatomical_smooth_um"] / sx))
    tissue = np.zeros(volume.shape, bool)
    weak_tissue = np.zeros(volume.shape, bool)
    thresholds = []
    minimum = max(1, int(rc["anatomical_minimum_area_um2"] / area))
    margin = max(0, int(round(rc["anatomical_tissue_margin_um"] / sy)))
    for z, plane in enumerate(smooth):
        if not np.isfinite(plane).all():
            raise ValueError(f"Non-finite IMS intensity in native z={z}")
        lo, hi = np.quantile(plane, [.15, .995])
        lower_domain = plane <= np.quantile(plane, .5)
        # Large spatial variation in defocused background is not camera noise.
        # Estimating noise from that variation would erase small end sections.
        local_difference = plane - ndi.gaussian_filter(plane, 1.5 / sy)
        lower = local_difference[lower_domain]
        noise = max(1e-6, 1.4826 * np.median(np.abs(lower - np.median(lower))))
        if hi - lo <= 3 * noise:
            thresholds.append(None)
            continue
        threshold = max(float(threshold_otsu(plane)), lo + 3 * noise)
        thresholds.append(threshold)
        highpass = plane - ndi.gaussian_filter(plane, rc["anatomical_focus_background_um"] / sy)
        # An absolute fraction of structural dynamic range, not a percentile
        # of the high-pass itself: the latter always selects a fraction of a
        # completely defocused haze, even when no sharp tissue remains.
        focus_threshold = max(noise * 2, float(hi - lo) * rc["anatomical_focus_fraction"])
        focused = highpass > focus_threshold
        weak = (highpass > max(noise * 2, float(hi-lo) * rc["anatomical_terminal_focus_fraction"])) & (plane > lo + 2 * noise)
        weak = binary_closing(weak, disk(max(1, round(1 / sy))))
        weak_tissue[z] = ndi.binary_fill_holes(remove_small_objects(weak, min_size=minimum))
        initial = plane > threshold
        labels, count = ndi.label(initial)
        selected = np.zeros_like(initial)
        for index in range(1, count + 1):
            component = labels == index
            if component.sum() < minimum:
                continue
            overlap = component & focused
            if overlap.sum() >= rc["anatomical_focused_support_fraction"] * component.sum():
                selected |= component
            else:
                # Separate sharply supported terminal islands from broad defocus
                # instead of rejecting the entire connected hazy cross-section.
                selected |= overlap
        selected = remove_small_objects(selected, min_size=minimum)
        if margin and selected.any():
            nearby = ndi.distance_transform_edt(~selected, sampling=(sy, sx)) <= margin * sy
            selected |= nearby & (plane > lo + 2 * noise)
        selected = binary_closing(selected, disk(max(1, round(2 / sy))))
        tissue[z] = ndi.binary_fill_holes(selected)
    return tissue, weak_tissue, smooth, thresholds


def _geometry_markers(tissue, distance, spacing, rc):
    separation = float(rc["anatomical_peak_separation_um"])
    size = tuple(2 * max(1, int(round(separation / s))) + 1 for s in spacing)
    maxima = (distance == ndi.maximum_filter(distance, size=size, mode="constant"))
    maxima &= distance >= float(rc["anatomical_minimum_core_radius_um"])
    labels, number = ndi.label(maxima)
    candidates = [np.unravel_index(np.argmax(np.where(labels == i, distance, -1)), distance.shape)
                  for i in range(1, number + 1)]
    candidates.sort(key=lambda c: float(distance[c]), reverse=True)
    if not candidates:
        raise ValueError("No reliable anatomical core found; provide a guide ROI or reduce anatomical_minimum_core_radius_um after checking tissue support")
    primary = np.array(candidates[0])
    radius = float(distance[candidates[0]])
    markers = np.zeros(tissue.shape, np.int32)
    records, flags = [], []
    for coord in candidates:
        xy_distance = np.linalg.norm((np.array(coord)[1:] - primary[1:]) * spacing[1:])
        ratio = xy_distance / max(radius, 1e-6)
        central = ratio <= float(rc["anatomical_central_radius_ratio"])
        markers[coord] = 1 if central else 2
        records.append({"z_y_x": [int(c) for c in coord], "radius_um": float(distance[coord]),
                        "distance_from_primary_radius_units": float(ratio), "assignment": "central" if central else "peripheral"})
        if abs(ratio - rc["anatomical_central_radius_ratio"]) < .25:
            flags.append("A geometric core lies near the central/peripheral cutoff; check seed assignment")
    if not (markers == 2).any():
        flags.append("No separate peripheral core found: optic-lobe removal is uncertain; a guide ROI is recommended")
    primary_xy = primary[1:] * spacing[1:]
    edge_distance = min(*primary_xy, *((np.array(tissue.shape[1:]) - 1) * spacing[1:] - primary_xy))
    if edge_distance < radius:
        flags.append("Primary anatomical core is close to a cropped image boundary")
    return markers, records, flags


def segment_anatomical_volume(volume: np.ndarray, spacing: tuple, roi_config: dict,
                              guide_mask: np.ndarray | None = None):
    """Array API, useful for synthetic testing; guide_mask uses this volume's XY.

    Automatic: 3-D distance cores -> merge nearby central cores -> marker
    watershed at necks -> smooth signed-distance boundary + supported margin.
    Guided: central/exterior seeds on either side of the guide's search band;
    the same 3-D neck/edge partition refines only that band. True-slice guides
    allow bounded outward growth away from their native z. Both retain native
    z, supported disconnected terminal islands, and filled central interiors.
    """
    rc = {**ANATOMICAL_DEFAULTS, **roi_config}
    if np.asarray(volume).ndim != 3 or len(spacing) != 3 or min(spacing) <= 0:
        raise ValueError("Expected a ZYX volume with positive ZYX spacing")
    if guide_mask is not None:
        if guide_mask.shape != volume.shape[1:] or not guide_mask.any():
            raise ValueError("Guide dimensions must match structural XY and contain a nonempty area")
        index = rc.get("guide_z_index")
        if index is not None and (isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < volume.shape[0]):
            raise ValueError(f"guide_z_index must be native z=0..{volume.shape[0]-1}, or null")
    nonnegative = [k for k, value in ANATOMICAL_DEFAULTS.items()
                   if isinstance(value, (float, int)) and not isinstance(value, bool)]
    for key in nonnegative:
        if not np.isfinite(rc[key]) or rc[key] < 0:
            raise ValueError(f"roi.{key} must be finite and nonnegative")
    if not 0 < rc["anatomical_terminal_fraction"] < 1:
        raise ValueError("anatomical_terminal_fraction must be in (0, 1)")
    spacing = tuple(float(v) for v in spacing)
    tissue, weak_tissue, smooth, thresholds = _tissue_support(np.asarray(volume, np.float32), spacing, rc)
    if not tissue.any():
        return tissue, {"method": "anatomical_guided" if guide_mask is not None else "anatomical",
                        "native_z_preserved": True, "review_flags": ["No tissue support found"],
                        "central_area_by_z_um2": [0.0]*len(tissue), "tissue_area_by_z_um2": [0.0]*len(tissue)}
    # Treat cropped XY faces as unknown, not as proof of a narrow lobe. Z gets
    # explicit empty caps. This permits a peripheral thickness maximum at a
    # cropped face, while still assigning its identity by relative geometry.
    distance = ndi.distance_transform_edt(np.pad(tissue, ((1, 1), (0, 0), (0, 0))), sampling=spacing)[1:-1].astype(np.float32)
    distance = ndi.gaussian_filter(distance, tuple(2.0 / v for v in spacing))
    flags, cores = [], []
    allowed = None
    guide_z = rc.get("guide_z_index")
    if guide_mask is None:
        markers, cores, flags = _geometry_markers(tissue, distance, spacing, rc)
    else:
        if guide_mask.shape != tissue.shape[1:] or not guide_mask.any():
            raise ValueError("Guide dimensions must match structural XY")
        if guide_z is not None and not 0 <= guide_z < tissue.shape[0]:
            raise ValueError(f"guide_z_index must be within native z=0..{tissue.shape[0]-1}")
        sd = _sdf(guide_mask, spacing[1:])
        growth = np.zeros(tissue.shape[0]) if guide_z is None else np.minimum(
            np.abs(np.arange(tissue.shape[0]) - guide_z) * spacing[0] * rc["guide_expansion_per_z_um"],
            rc["guide_max_expansion_um"])
        guide_distance = sd[None] + growth[:, None, None]
        band = float(rc["guide_refine_band_um"])
        markers = np.zeros(tissue.shape, np.int32)
        markers[tissue & (guide_distance >= band)] = 1
        markers[tissue & (guide_distance <= -band)] = 2
        allowed = guide_distance >= -band
        if not (markers == 1).any():
            raise ValueError("Guide has no tissue-supported interior seed; check alignment or reduce guide_refine_band_um")
        if not (markers == 2).any():
            flags.append("No tissue-supported exterior guide seed; refinement uses the guide envelope only")
        if guide_z is not None and not (tissue[guide_z] & guide_mask).any():
            raise ValueError("The guide does not overlap tissue in guide_z_index; check slice index/alignment")
    # Tissue-scale gradients only weakly influence an already anatomical cut.
    gy = ndi.sobel(smooth, axis=1) / spacing[1]
    gx = ndi.sobel(smooth, axis=2) / spacing[2]
    gradient = np.hypot(gy, gx)
    del gy, gx
    scale = max(float(np.quantile(gradient[tissue], .95)), 1e-6)
    elevation = -distance / max(float(distance.max()), 1e-6) + rc["anatomical_edge_weight"] * np.clip(gradient / scale, 0, 1)
    del gradient, distance, smooth
    partition = watershed(elevation, markers, mask=tissue)
    central = partition == 1
    if allowed is not None:
        central &= allowed
    del elevation, markers, partition
    if not central.any():
        raise ValueError("Anatomical partition has no central tissue")
    raw_area = central.sum(axis=(1, 2))
    terminal = raw_area < rc["anatomical_terminal_fraction"] * raw_area.max()
    # Relax focus ONLY near a detected central section in z and XY. This
    # preserves weak terminal islands without admitting an entire defocus haze.
    # Use the original central volume, not iteratively rescued masks, so the
    # rescue cannot propagate indefinitely through empty/background z layers.
    context = int(np.ceil(rc["anatomical_terminal_context_um"] / spacing[0]))
    rescue_count = 0
    for z in np.flatnonzero(terminal):
        nearby = central[max(0, z-context):min(len(central), z+context+1)].any(axis=0)
        if not nearby.any():
            continue
        permitted = _sdf(nearby, spacing[1:]) >= -rc["anatomical_terminal_margin_um"]
        addition = weak_tissue[z] & permitted
        if allowed is not None:
            addition &= allowed[z]
        rescue_count += int((addition & ~tissue[z]).sum())
        tissue[z] |= addition
    del weak_tissue
    # SDF smoothing controls indentations and z jumps without enforcing an ellipse
    # or deleting small components just because they touch an image border.
    sdf = np.stack([_sdf(ndi.binary_fill_holes(p), spacing[1:]) if p.any()
                    else np.full(p.shape, -max(5, rc["anatomical_boundary_smooth_um"]*3), np.float32)
                    for p in central])
    smoothed = ndi.gaussian_filter(sdf, (rc["anatomical_z_smooth_um"] / spacing[0],
                                        rc["anatomical_boundary_smooth_um"] / spacing[1],
                                        rc["anatomical_boundary_smooth_um"] / spacing[2]))
    central = (smoothed >= -rc["anatomical_central_margin_um"]) & tissue
    del sdf, smoothed
    if guide_mask is not None:
        central &= allowed
        terminal_envelope = allowed
    else:
        stable_projection = central[~terminal].any(axis=0) if (~terminal).any() else central.any(axis=0)
        terminal_envelope = (_sdf(stable_projection, spacing[1:]) >= -rc["anatomical_terminal_margin_um"])[None]
    rescue = tissue & terminal[:, None, None] & terminal_envelope
    # Rescue is intersected with CURRENT tissue; it is not an extruded 2-D ROI.
    central |= rescue
    for z in range(len(central)):
        central[z] = ndi.binary_fill_holes(central[z]) & tissue[z]
        if allowed is not None:
            central[z] &= allowed[z]
    diagnostics = {
        "method": "anatomical_guided" if guide_mask is not None else "anatomical",
        "native_z_preserved": True, "geometric_cores": cores,
        "review_flags": sorted(set(flags)), "guide_kind": ("projection" if guide_z is None else "native_z_slice") if guide_mask is not None else None,
        "guide_z_index": guide_z if guide_mask is not None else None,
        "tissue_threshold_by_z": thresholds,
        "tissue_area_by_z_um2": (tissue.sum(axis=(1, 2))*spacing[1]*spacing[2]).tolist(),
        "central_area_by_z_um2": (central.sum(axis=(1, 2))*spacing[1]*spacing[2]).tolist(),
        "terminal_z_indices": np.flatnonzero(terminal).tolist(),
        "relaxed_terminal_support_voxels": rescue_count,
        "retained_fraction": float(central.sum()/tissue.sum()),
        "limitation": "Geometry-based ROI, not validated anatomical labels. Review optical slices, especially cropped/ambiguous brains."
    }
    return central, diagnostics


def build_anatomical_roi(reader: Any, config: dict):
    volume, spacing, level = read_structural_volume(reader, config)
    rc = {**ANATOMICAL_DEFAULTS, **config["roi"]}
    guide = None
    if rc["mode"] == "guided_central_brain":
        time, channel = config["input"]["time_index"], config["input"]["channel_index"]
        native_shape = reader.metaData[0, time, channel, "shape"][-3:]
        guide = load_guide_mask(rc["guide_mask_path"], tuple(native_shape[1:]), allow_resize=rc["guide_allow_resize"])
        guide = resize(guide.astype(np.uint8), volume.shape[1:], order=0, preserve_range=True, anti_aliasing=False) > 0
    mask, diagnostics = segment_anatomical_volume(volume, spacing, rc, guide)
    for message in diagnostics["review_flags"]:
        warnings.warn(message, UserWarning, stacklevel=2)
    diagnostics["guide_mask_path"] = rc["guide_mask_path"] if guide is not None else None
    return mask, spacing, level, diagnostics
