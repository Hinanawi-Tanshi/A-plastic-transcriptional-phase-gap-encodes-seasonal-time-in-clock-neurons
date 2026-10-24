"""Prepare native-XY IMS projections and per-sample central-brain guides."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from .anatomical_roi import compatible_guide_path
from .paths import PathInput, discover_ims_files


def plan_guide_projections(
    inputs: PathInput,
    output_dir: str | Path | None = None,
    *,
    recursive: bool = False,
    pattern: str = "*.ims",
) -> pd.DataFrame:
    """Plan native MIPs and guides without creating any files.

    By default each IMS parent gets its own ``central_brain_roi`` folder. The
    optional ``output_dir`` retains the earlier shared-folder layout.
    Guide filenames now follow the MIP basename, with legacy names accepted.
    """
    shared_root = Path(output_dir).resolve() if output_dir is not None else None
    mappings: dict[Path, dict] = {}
    records = []
    for source in discover_ims_files(inputs, pattern, recursive=recursive):
        if shared_root is None:
            root = source.parent / "central_brain_roi"
            projection = root / f"{source.stem}_max_projection.tif"
            guide = root / f"{source.stem}_max_projection.roi"
            manifest_key = source.name
        else:
            root = shared_root
            suffix = hashlib.sha256(str(source).casefold().encode("utf-8")).hexdigest()[:12]
            base = f"{source.stem}_{suffix}"
            projection = root / "projections" / f"{base}_max_projection.tif"
            guide = root / "guides" / f"{base}_max_projection.roi"
            manifest_key = str(source)
        manifest_path = root / "guide_manifest.json"
        if manifest_path not in mappings:
            mapping = json.loads(manifest_path.read_text(encoding="utf-8-sig")) if manifest_path.is_file() else {}
            if not isinstance(mapping, dict):
                raise ValueError(f"Guide manifest must be a JSON object: {manifest_path}")
            mappings[manifest_path] = {
                str(key).replace("\\", "/").casefold(): value for key, value in mapping.items()
            }
        normalized = mappings[manifest_path]
        entry = normalized.get(str(source).replace("\\", "/").casefold())
        if entry is None:
            entry = normalized.get(source.name.casefold())
        guide_used = guide
        if isinstance(entry, dict) and entry.get("mask_path"):
            guide_used = Path(entry["mask_path"])
            if not guide_used.is_absolute():
                guide_used = root / guide_used
        guide_used = compatible_guide_path(guide_used)
        records.append({
            "input_file": str(source),
            "manifest_path": str(manifest_path),
            "manifest_key": manifest_key,
            "projection_tif": str(projection),
            "guide_roi": str(guide),
            "guide_path_used": str(guide_used),
            "projection_exists": projection.is_file(),
            "guide_exists": guide_used.is_file(),
        })
    return pd.DataFrame.from_records(records)


def export_guide_projections(
    inputs: PathInput,
    output_dir: str | Path | None = None,
    *,
    recursive: bool = False,
    pattern: str = "*.ims",
    time_index: int = 0,
    channel_index: int = 0,
    overwrite_projections: bool = False,
    print_plane_progress: bool = False,
) -> pd.DataFrame:
    """Export raw-intensity XY MIPs and write one manifest per data directory.

    Existing projections are preserved unless explicitly overwritten. Existing
    manifest entries (including user-edited .tif mask paths) are never replaced;
    only newly discovered sources are appended. ROI files are never touched.
    """
    from export_max_projection import calculate_max_projection, save_projection_tiff

    plan = plan_guide_projections(inputs, output_dir, recursive=recursive, pattern=pattern)
    mappings: dict[Path, dict] = {}
    actions = []
    for row in plan.itertuples(index=False):
        source = Path(row.input_file)
        projection_path = Path(row.projection_tif)
        manifest_path = Path(row.manifest_path)
        if manifest_path not in mappings:
            mapping = json.loads(manifest_path.read_text(encoding="utf-8-sig")) if manifest_path.is_file() else {}
            if not isinstance(mapping, dict):
                raise ValueError(f"Guide manifest must be a JSON object: {manifest_path}")
            mappings[manifest_path] = mapping
        if projection_path.exists() and not overwrite_projections:
            action = "existing"
        else:
            projection, spacing, z_interval = calculate_max_projection(
                source, time_index=time_index, channel_index=channel_index,
                resolution_level=0, print_progress=print_plane_progress,
            )
            save_projection_tiff(
                projection_path, projection, spacing_zyx=spacing,
                metadata={
                    "axes": "YX", "source_ims": str(source),
                    "time_index": time_index, "channel_index": channel_index,
                    "resolution_level": 0, "z_interval_zero_based_half_open": list(z_interval),
                    "projection": "maximum", "intensity": "unscaled native values",
                },
                overwrite=overwrite_projections,
            )
            action = "exported"
        actions.append(action)
        mappings[manifest_path].setdefault(row.manifest_key, {
            "mask_path": str(Path(row.guide_roi).relative_to(manifest_path.parent)),
            "z_index": None,
        })
    for manifest_path, mapping in mappings.items():
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        if output_dir is not None:
            (manifest_path.parent / "guides").mkdir(exist_ok=True)
        manifest_path.write_text(json.dumps(mapping, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    result = plan.copy()
    result["projection_action"] = actions
    return result
