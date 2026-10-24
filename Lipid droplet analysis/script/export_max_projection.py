"""Native-intensity maximum projections for the manual-ROI workflow."""
from __future__ import annotations
import json
from pathlib import Path
from typing import Sequence
import numpy as np
import tifffile
from imaris_ims_file_reader.ims import ims


def calculate_max_projection(
    ims_path: str | Path,
    *,
    time_index: int = 0,
    channel_index: int = 0,
    resolution_level: int = 0,
    z_start: int = 0,
    z_stop: int | None = None,
    print_progress: bool = True,
) -> tuple[np.ndarray, tuple[float, float, float], tuple[int, int]]:
    """Return a streaming Z maximum projection, ZYX spacing, and z interval."""
    ims_path = Path(ims_path)
    if not ims_path.is_file():
        raise FileNotFoundError(f"IMS file does not exist: {ims_path}")
    if resolution_level < 0:
        raise ValueError("resolution_level must be at least 0")

    reader = ims(
        str(ims_path),
        ResolutionLevelLock=resolution_level,
        squeeze_output=True,
    )
    try:
        if resolution_level >= int(reader.ResolutionLevels):
            raise IndexError(
                f"resolution_level={resolution_level}, but the file has "
                f"{reader.ResolutionLevels} level(s)"
            )
        if not 0 <= time_index < int(reader.TimePoints):
            raise IndexError(
                f"time_index={time_index}, but the file has {reader.TimePoints} time point(s)"
            )
        if not 0 <= channel_index < int(reader.Channels):
            raise IndexError(
                f"channel_index={channel_index}, but the file has {reader.Channels} channel(s)"
            )

        shape_zyx = tuple(int(value) for value in reader.shape[-3:])
        spacing_zyx = tuple(float(value) for value in reader.resolution)
        number_of_planes = shape_zyx[0]
        stop = number_of_planes if z_stop is None else int(z_stop)
        start = int(z_start)
        if not 0 <= start < stop <= number_of_planes:
            raise ValueError(
                f"Expected 0 <= z_start < z_stop <= {number_of_planes}; "
                f"received z_start={start}, z_stop={stop}"
            )

        projection: np.ndarray | None = None
        total = stop - start
        progress_step = max(1, total // 10)
        for completed, z_index in enumerate(range(start, stop), start=1):
            plane = np.asarray(reader[time_index, channel_index, z_index, :, :])
            if projection is None:
                projection = plane.copy()
            else:
                np.maximum(projection, plane, out=projection)
            if print_progress and (
                completed == 1 or completed == total or completed % progress_step == 0
            ):
                print(
                    f"  z {completed}/{total} (source index {z_index})",
                    flush=True,
                )
        if projection is None:  # protected by interval validation
            raise RuntimeError("No z planes were read")
        return projection, spacing_zyx, (start, stop)
    finally:
        reader.close()


def save_projection_tiff(
    output_path: str | Path,
    image: np.ndarray,
    *,
    spacing_zyx: Sequence[float],
    metadata: dict[str, object],
    compression: str | None = "zlib",
    overwrite: bool = False,
) -> Path:
    """Save a grayscale or RGB TIFF with ImageJ-readable XY calibration."""
    output_path = Path(output_path)
    if output_path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {output_path}. Set OVERWRITE=True or use --overwrite."
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    z_um, y_um, x_um = map(float, spacing_zyx)
    if min(z_um, y_um, x_um) <= 0:
        raise ValueError(f"Invalid voxel spacing: {spacing_zyx}")
    description = json.dumps(metadata, ensure_ascii=False, indent=2)
    photometric = "rgb" if image.ndim == 3 else "minisblack"
    tifffile.imwrite(
        output_path,
        image,
        photometric=photometric,
        compression=compression,
        # TIFF resolution tags are pixels per centimeter. This preserves the
        # XY pixel size when the result is opened in ImageJ/Fiji.
        resolution=(10000.0 / x_um, 10000.0 / y_um),
        resolutionunit="CENTIMETER",
        description=description,
    )
    return output_path.resolve()
