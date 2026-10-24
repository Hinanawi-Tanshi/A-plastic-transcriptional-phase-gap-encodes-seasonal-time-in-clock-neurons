"""ROI-only notebook previews: build the 3-D ROI once, inspect many true z planes."""
from pathlib import Path
import json
import math

import matplotlib.pyplot as plt
import numpy as np
from imaris_ims_file_reader.ims import ims
from skimage.transform import resize

from .anatomical_roi import resolve_sample_guide
from .pipeline import _make_roi_provider, _save_roi_mask_tiff


def preview_roi(ims_path, config, output_dir, z_indices=None, *, export_native_mask=False):
    """Save optical-slice overlays, area-vs-z, JSON and optional native ZYX mask.

    Slice indices refer to config.input.resolution_level (normally 0).
    Guide indices always refer to native level 0. No LD detection is performed.
    Returns a dict of paths; use IPython.display.Image for inline PNG display.
    """
    from .anatomical_v2 import resolve_v2_config
    config = resolve_v2_config(config)
    config = resolve_sample_guide(config, ims_path)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    reader = ims(str(ims_path), ResolutionLevelLock=int(config["input"]["resolution_level"]), squeeze_output=True)
    paths = {}
    try:
        shape = tuple(reader.shape[-3:])
        spacing = tuple(float(s) for s in reader.resolution)
        indices = list(dict.fromkeys(int(z) for z in (z_indices if z_indices is not None else np.linspace(0, shape[0]-1, 12).round())))
        if not indices or min(indices) < 0 or max(indices) >= shape[0]:
            raise ValueError(f"Preview indices must be within 0..{shape[0]-1}")
        provider, metadata = _make_roi_provider(reader, config, keep_roi_stages=True)
        if provider.roi_artifacts is not None:
            from .roi_qc import save_structural_roi_qc
            pack = provider.roi_artifacts
            paths.update(save_structural_roi_qc(output / "structural", pack["volume"], provider.low_mask,
                         provider.low_spacing, metadata, pack["stages"], title=Path(ims_path).name))
            paths["structural_slices"] = paths.pop("slices")
            paths["structural_depth"] = paths.pop("depth")
        ncols = min(4, len(indices))
        fig, axes = plt.subplots(math.ceil(len(indices)/ncols), ncols,
                                 figsize=(4*ncols, 4*math.ceil(len(indices)/ncols)), squeeze=False, layout="constrained")
        for axis in axes.flat:
            axis.axis("off")
        for axis, z in zip(axes.flat, indices):
            raw = np.asarray(reader[config["input"]["time_index"], config["input"]["channel_index"], z, :, :])
            target = tuple(max(1, round(n*min(1, 700/max(shape[1:])))) for n in shape[1:])
            small = resize(raw, target, preserve_range=True, anti_aliasing=True)
            mask = resize(provider.plane(z).astype(np.uint8), target, order=0, preserve_range=True, anti_aliasing=False) > 0
            axis.imshow(small, cmap="gray", vmin=np.quantile(small, .02), vmax=np.quantile(small, .998))
            if mask.any() and not mask.all():
                axis.contour(mask, levels=[.5], colors="cyan", linewidths=.7)
            axis.set_title(f"z={z} (ImageJ {z+1}); ROI cyan")
        fig.suptitle(f"{Path(ims_path).name}\n{metadata.get('method', metadata['mode'])}; independent display ranges")
        paths["slices"] = output / "roi_optical_slices.png"
        fig.savefig(paths["slices"], dpi=140)
        plt.close(fig)
        if provider.low_mask is not None:
            areas = provider.low_mask.sum(axis=(1, 2)) * np.prod(provider.low_spacing[1:])
            depths = np.arange(len(areas)) * provider.low_spacing[0]
        else:
            areas = np.array([provider.plane(z).sum() * np.prod(spacing[1:]) for z in range(shape[0])])
            depths = np.arange(len(areas)) * spacing[0]
        fig, ax = plt.subplots(figsize=(10, 3.5), layout="constrained")
        ax.plot(depths, areas, label="Final ROI")
        if "base_area_by_z_um2" in metadata:
            ax.plot(depths, metadata["base_area_by_z_um2"], label="Base ROI")
            ax.plot(depths, metadata["tissue_area_by_z_um2"], label="Tissue", alpha=.65)
            ax.legend()
        ax.set(xlabel="Depth from first slice (um)", ylabel="ROI area (um2)", title=metadata.get('roi_scope',metadata['mode']))
        ax.grid(alpha=.2)
        paths["depth"] = output / "roi_area_by_z.png"
        fig.savefig(paths["depth"], dpi=140)
        plt.close(fig)
        metadata.update({"ims_path": str(Path(ims_path).resolve()), "config": config,
                         "preview_z_indices": indices})
        paths["metadata"] = output / "roi_diagnostics.json"
        paths["metadata"].write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
        if export_native_mask:
            if config["input"]["resolution_level"] != 0:
                raise ValueError("Native ROI export requires input.resolution_level=0")
            paths["mask"] = output / ("whole_brain_roi_native.tif" if metadata.get('roi_scope')=='whole_brain' else "central_brain_roi_native.tif")
            _save_roi_mask_tiff(paths["mask"], provider, spacing)
    finally:
        reader.close()
    return paths
