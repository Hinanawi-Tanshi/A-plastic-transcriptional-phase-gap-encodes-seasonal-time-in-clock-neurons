"""Direct 3D RNAscope quantification: LD4 and dual-channel LP.

Publication subset of the original TWX 2026 module. Detection, decomposition,
rounding, mask-1..5 counting and group statistics retain their original rules.
Batch execution with validated checkpoints is in rnascope_workflow.py.
"""
import os
import re
import gc
import logging
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import SimpleITK as sitk
from tifffile import imwrite
from tqdm import tqdm
import bigfish.detection as detection
from imaris_ims_file_reader.ims import ims


def _validate_image_mask(image, mask):
    """Reject misaligned input; preserve the original missing-label-as-zero rule."""
    if image.ndim != 3 or image.shape != mask.shape:
        raise ValueError(f"Expected matching ZYX arrays: image={image.shape}, mask={mask.shape}")
    if not np.issubdtype(mask.dtype, np.integer) or np.any(mask < 0):
        raise ValueError("Cell masks must contain non-negative integer labels")


def _save_spots(directory, channel, final_spots, raw_spots, shape, params):
    """Save coordinates from this detection, without running detection again."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(directory / f"{channel}_spots.npz",
                        final_spots=final_spots, raw_spots=raw_spots,
                        image_shape_zyx=np.asarray(shape),
                        voxel_size_nm=np.asarray(params['voxel_size_nm']))


def extract_and_process_channel(ims_path, output_dir, channel_idx=1, p_low=1, p_high=99.995):
    """
    Extract a single channel from an IMS file, apply log transform and percentile
    normalization, then save as NRRD.

    Args:
        ims_path (Path): Path to the input IMS file.
        output_dir (Path): Directory for the output NRRD file.
        channel_idx (int): Channel index (0 = first channel, 1 = second/GFP).
        p_low (float): Low percentile for normalization.
        p_high (float): High percentile for normalization.
    """
    try:
        reader = ims(str(ims_path))
        try:
            data = np.asarray(reader[0, channel_idx])
        finally:
            reader.close()

        # Log transform (log1p to handle zero-valued pixels)
        log_data = np.log1p(data.astype(np.float32))

        del data  # Free original data memory
        gc.collect()

        # Global percentile normalization
        val_low, val_high = np.percentile(log_data, (p_low, p_high))

        # Linear mapping to 0-255
        if val_high - val_low > 0:
            normalized_data = (log_data - val_low) / (val_high - val_low) * 255.0
        else:
            normalized_data = np.zeros_like(log_data)

        np.clip(normalized_data, 0, 255, out=normalized_data)

        # Convert to uint8
        final_data = normalized_data.astype(np.uint8)

        del log_data, normalized_data
        gc.collect()

        # Save file
        output_path = os.path.join(output_dir, f"{ims_path.stem}_CH{channel_idx}.nrrd")
        sitk.WriteImage(sitk.GetImageFromArray(final_data), str(output_path))
        logging.info(f"Successfully processed: {output_path}")

    except Exception as e:
        logging.error(f"Error processing file {ims_path.name}: {str(e)}")


def batch_process_ims(input_dir, output_dir, channel_idx=1, p_low=1, p_high=99.99,
                      name_filter=None):
    """
    Batch process all IMS files in a directory.

    Args:
        input_dir: Directory containing IMS files.
        output_dir: Output directory.
        channel_idx: Channel index.
        p_low: Low percentile for normalization.
        p_high: High percentile for normalization.
        name_filter: Optional filename keyword filter (e.g., "ZT18"). None = process all.
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all .ims files
    ims_files = list(input_dir.glob("*.ims"))

    # --- File filtering ---
    if name_filter:
        ims_files = [f for f in ims_files if name_filter in f.name]
        logging.info(f"Applied name filter '{name_filter}': matched {len(ims_files)} files")

    if not ims_files:
        logging.warning(f"No .ims files found in directory {input_dir}!")
        return

    logging.info(f"Found {len(ims_files)} IMS files, starting batch processing...")

    for ims_file in tqdm(ims_files, desc="Processing IMS files"):
        extract_and_process_channel(ims_file, output_dir, channel_idx, p_low, p_high)

    logging.info("Batch processing complete!")


def detect_spot(
    image,
    z_radius_nm,
    xy_radius_nm,
    threshold=None,
    voxel_size_nm=(500, 96, 96),          # (Z, Y, X) in nm
    alpha=0.7,
    beta=1.0,
    gamma=15,
    cluster_radius_nm=1000,
    nb_min_spots=10,
    decompose=True,
    cluster=True,
    verbose=True
):
    """
    Detect 3D RNA FISH signal spots using fish-quant (bigfish).

    Args:
        image: 3D array (Z, Y, X), raw image.
        z_radius_nm: Spot Z radius in nm.
        xy_radius_nm: Spot XY radius in nm.
        threshold: Detection threshold. None = auto-estimate.
        voxel_size_nm: (Z, Y, X) voxel size in nm.
        alpha, beta, gamma: Decomposition parameters.
        cluster_radius_nm: Clustering radius in nm.
        nb_min_spots: Minimum spots per cluster.
        decompose: Whether to perform dense region decomposition.
        cluster: Whether to perform clustering.
        verbose: Whether to print intermediate info.

    Returns:
        spots_final: (N, 3) array, final spot coordinates [z, y, x] (0-based, float).
        clusters: (M, 4) array or None, clusters [z, y, x, radius_nm] if cluster=True.
        spots_raw: (K, 3) array, initial spots (for debugging).
    """
    # Input validation
    if image.ndim != 3:
        raise ValueError("Input image must be 3D (Z, Y, X)")
    if not isinstance(voxel_size_nm, (tuple, list)) or len(voxel_size_nm) != 3:
        raise ValueError("voxel_size_nm must be a 3-element tuple (Z, Y, X)")

    # Uniform parameter format
    spot_radius_nm = (z_radius_nm, xy_radius_nm, xy_radius_nm)

    # 1. Initial detection
    spots_raw, estimated_threshold = detection.detect_spots(
        images=image,
        threshold=threshold,
        return_threshold=True,
        voxel_size=voxel_size_nm,
        spot_radius=spot_radius_nm
    )

    if verbose:
        print(f"Initial detection: {len(spots_raw)} spots (threshold: {estimated_threshold:.2f})")

    # 2. Dense region decomposition (optional)
    if decompose:
        spots_decomposed, _, _ = detection.decompose_dense(
            image=image,
            spots=spots_raw,
            voxel_size=voxel_size_nm,
            spot_radius=spot_radius_nm,
            alpha=alpha,
            beta=beta,
            gamma=gamma
        )
        spots_final = spots_decomposed
        if verbose:
            print(f"After decomposition: {len(spots_final)} spots")
    else:
        spots_final = spots_raw

    # 3. Clustering (optional)
    clusters_out = None
    if cluster and len(spots_final) > 0:
        spots_clustered, clusters_out = detection.detect_clusters(
            spots=spots_final,
            voxel_size=voxel_size_nm,
            radius=cluster_radius_nm,
            nb_min_spots=nb_min_spots
        )
        if verbose:
            print(f"Clusters detected: {len(clusters_out) if clusters_out is not None else 0}")

    return spots_final, clusters_out, spots_raw


def extract_zt_time(filename):
    """
    Extract ZT time from filename.
    Example: "17385_488V5_568TimCDS_DAPI_LD4_ZT17_LN1_1_2026-04-12_17.06.18.ims" -> 17

    Args:
        filename (str): IMS filename.

    Returns:
        int or None: ZT time value.
    """
    match = re.search(r'ZT(\d+)', filename)
    if match:
        return int(match.group(1))
    return None


def load_mask(mask_path):
    """
    Load mask NRRD file.

    Args:
        mask_path: Path to mask NRRD file.

    Returns:
        3D numpy array, mask labels start from 1.
    """
    mask_sitk = sitk.ReadImage(str(mask_path))
    mask_array = sitk.GetArrayFromImage(mask_sitk)  # (Z, Y, X)
    return mask_array


def load_rna_channel(ims_path, channel_idx=0):
    """
    Load a specified channel from an IMS file (raw intensity, no preprocessing).

    Args:
        ims_path: Path to IMS file.
        channel_idx: Channel index.

    Returns:
        3D numpy array (Z, Y, X).
    """
    reader = ims(str(ims_path))
    try:
        return np.asarray(reader[0, channel_idx])
    finally:
        reader.close()


def count_spots_in_mask(spots, mask_array):
    """
    Count spots within each mask label.

    Args:
        spots: (N, 3) array, spot coordinates [z, y, x] (voxel coordinates, float).
        mask_array: 3D array, mask labels (1-5).

    Returns:
        dict: {mask_label: spot_count}
    """
    counts = defaultdict(int)

    for spot in spots:
        z, y, x = spot
        # Convert to integer index
        z_int, y_int, x_int = int(round(z)), int(round(y)), int(round(x))

        # Boundary check
        if (0 <= z_int < mask_array.shape[0] and
            0 <= y_int < mask_array.shape[1] and
            0 <= x_int < mask_array.shape[2]):

            label = mask_array[z_int, y_int, x_int]
            if label > 0:  # Mask labels start from 1
                counts[label] += 1

    return dict(counts)


def process_single_sample(ims_path, mask_path, detect_params=None, verbose=True, spot_output_dir=None):
    """
    Process a single sample: load image → load mask → detect spots → count by mask.

    Args:
        ims_path: Path to IMS raw image.
        mask_path: Path to corresponding mask file.
        detect_params: Detection parameter dict.
        verbose: Whether to print detailed info.

    Returns:
        dict: Results including filename, ZT time, and per-mask spot counts.
    """
    # Default detection parameters
    if detect_params is None:
        detect_params = {
            'z_radius_nm': 500,
            'xy_radius_nm': 96,
            'threshold': None,  # Auto-estimate
            'voxel_size_nm': (500, 96, 96),
            'decompose': True,
            'cluster': False,
            'verbose': verbose
        }

    filename = Path(ims_path).stem
    zt_time = extract_zt_time(filename)

    try:
        # 1. Load RNA FISH channel (CH0)
        if verbose:
            print(f"Loading image: {filename}")
        rna_image = load_rna_channel(ims_path, channel_idx=0)

        # 2. Load mask
        if verbose:
            print(f"Loading mask: {Path(mask_path).name}")
        mask_array = load_mask(mask_path)

        # 3. Detect spots
        if verbose:
            print("Detecting spots...")
        _validate_image_mask(rna_image, mask_array)
        spots, _, raw_spots = detect_spot(rna_image, **detect_params)
        if spot_output_dir is not None:
            _save_spots(spot_output_dir, 'CH0', spots, raw_spots, rna_image.shape, detect_params)

        # Free image memory
        del rna_image
        gc.collect()

        # 4. Count spots per mask
        spot_counts = count_spots_in_mask(spots, mask_array)

        # Build result dict
        result = {
            'filename': filename,
            'zt_time': zt_time,
            'total_spots': len(spots),
        }

        # Add per-mask spot counts (mask 1-5)
        for i in range(1, 6):
            result[f'mask_{i}'] = spot_counts.get(i, 0)

        if verbose:
            print(f"ZT{zt_time}: total spots={len(spots)}, per-mask counts={spot_counts}")

        return result

    except Exception as e:
        logging.error(f"Error processing {filename}: {str(e)}")
        return None


def calculate_cell_type_statistics(df):
    """
    Calculate cell-type statistics from per-mask spot counts.

    sLNv cells: mean of mask 1-4
    5th sLNv cell: mask 5

    Args:
        df: DataFrame with columns mask_1..mask_5, zt_time.

    Returns:
        DataFrame: Statistics grouped by ZT time.
    """
    # Compute sLNv (mean of mask 1-4)
    df['sLNv'] = df[['mask_1', 'mask_2', 'mask_3', 'mask_4']].mean(axis=1)
    df['sLNv_std'] = df[['mask_1', 'mask_2', 'mask_3', 'mask_4']].std(axis=1)

    # 5th sLNv is mask 5
    df['sLNv_5th'] = df['mask_5']

    # Group by ZT time
    stats = df.groupby('zt_time').agg({
        'sLNv': ['mean', 'std', 'sem', 'count'],
        'sLNv_5th': ['mean', 'std', 'sem', 'count']
    }).reset_index()

    # Flatten column names
    stats.columns = ['zt_time',
                     'sLNv_mean', 'sLNv_std', 'sLNv_sem', 'sLNv_n',
                     'sLNv_5th_mean', 'sLNv_5th_std', 'sLNv_5th_sem', 'sLNv_5th_n']

    return stats


def plot_expression_curves(stats, output_path=None, figsize=(8, 4), lights_off=12):
    """
    Plot RNA expression time-series curves with error bars (SEM).

    Args:
        stats: Statistics DataFrame (output of calculate_cell_type_statistics).
        output_path: Path to save the figure.
        figsize: Figure size (width, height).
        lights_off: Lights-off ZT time (default 12). Gray shading from this point to right.

    Returns:
        matplotlib Figure.
    """
    import matplotlib.ticker as ticker

    fig, ax = plt.subplots(figsize=figsize)

    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.labelsize'] = 16
    plt.rcParams['legend.fontsize'] = 14

    x = stats['zt_time']

    ax.errorbar(x, stats['sLNv_mean'],
                yerr=stats['sLNv_sem'],
                fmt='o-', color="#FFA500", linewidth=2.5, markersize=8,
                capsize=5, capthick=2, label='sLNv',
                elinewidth=2)

    ax.errorbar(x, stats['sLNv_5th_mean'],
                yerr=stats['sLNv_5th_sem'],
                fmt='s-', color="#0000FF", linewidth=2.5, markersize=8,
                capsize=5, capthick=2, label='5th sLNv',
                elinewidth=2)

    ax.set_xlabel('ZT', fontweight='bold')
    ax.set_ylabel('RNA Spots per Cell', fontweight='bold')

    # X axis: major tick every 2h, minor tick every 1h
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))

    # Remove top/right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Thicker left/bottom spines
    for spine in ['left', 'bottom']:
        ax.spines[spine].set_linewidth(1.5)

    # Thicker ticks
    ax.tick_params(which='major', width=1.5, length=6, labelsize=14)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    ax.tick_params(which='minor', width=1, length=3)

    # Lights-off shading: gray from lights_off to end
    x_max = max(x.max(), 24)
    ax.axvspan(lights_off, x_max, color='gray', alpha=0.30, zorder=0)

    ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), framealpha=0.9,
              borderaxespad=0)

    fig.tight_layout()

    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {output_path}")

    plt.show()
    return fig


def visualize_spots_as_image(spots, image_shape, output_path=None):
    """
    Render detected spots as a 3D single-channel image (TIFF).

    Args:
        spots: (N, 3) array, spot coordinates [z, y, x].
        image_shape: Tuple, original image shape (Z, Y, X).
        output_path: Path to save the TIFF file.

    Returns:
        spot_image: 3D array, each voxel value = number of spots at that location.
    """
    # Create blank image
    spot_image = np.zeros(image_shape, dtype=np.uint16)

    # Count spots per voxel
    for spot in spots:
        z, y, x = spot
        z_int, y_int, x_int = int(round(z)), int(round(y)), int(round(x))

        # Boundary check
        if (0 <= z_int < image_shape[0] and
            0 <= y_int < image_shape[1] and
            0 <= x_int < image_shape[2]):
            spot_image[z_int, y_int, x_int] += 1

    print(f"Spot image statistics:")
    print(f"  - Total spots: {len(spots)}")
    print(f"  - Non-zero voxels: {np.count_nonzero(spot_image)}")
    print(f"  - Max spots per voxel: {spot_image.max()}")

    # Save as TIFF
    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        imwrite(output_path, spot_image)
        print(f"  - Saved to: {output_path}")

    return spot_image


def extract_lp_info(filename):
    """
    Extract light-pulse experiment info from LP FISH filename.

    Examples:
        ..._LP_ZT15_16_60min_LN1_9_...  → zt_range="ZT15-16", condition="60min", exp_group="Light"
        ..._LP_ZT15_16_ZT16_LN1_9_...    → zt_range="ZT15-16", condition="ZT16", exp_group="Dark"

    Args:
        filename (str): LP experiment filename.

    Returns:
        dict: {'zt_range', 'zt_start', 'zt_end', 'condition', 'exp_group'} or None.
    """
    # Match _LP_ZT{start}_{end}_{condition}_
    m = re.search(r'_LP_ZT(\d+)_(\d+)_(\d+min|ZT\d+)_', filename)
    if not m:
        return None

    zt_start = int(m.group(1))
    zt_end = int(m.group(2))
    condition = m.group(3)

    # Determine experimental group
    if condition.endswith('min'):
        exp_group = 'Light'
    elif condition.startswith('ZT'):
        exp_group = 'Dark'
    else:
        exp_group = 'Unknown'

    return {
        'zt_range': f'ZT{zt_start}-{zt_end}',
        'zt_start': zt_start,
        'zt_end': zt_end,
        'condition': condition,
        'exp_group': exp_group
    }


def detect_spots_for_channel(image, mask_array, detect_params, channel_name='CH0', spot_output_dir=None):
    """
    Detect spots for a single channel and count per mask.

    Args:
        image: 3D numpy array.
        mask_array: 3D mask array.
        detect_params: Detection parameter dict.
        channel_name: Channel name (for logging).

    Returns:
        dict: {'spots': array, 'total': int, 'mask_1'..'mask_5': int}
    """
    _validate_image_mask(image, mask_array)
    spots, _, raw_spots = detect_spot(image, **detect_params)
    if spot_output_dir is not None:
        _save_spots(spot_output_dir, channel_name, spots, raw_spots, image.shape, detect_params)
    counts = count_spots_in_mask(spots, mask_array)

    result = {'spots': spots, 'total': len(spots)}
    for i in range(1, 6):
        result[f'mask_{i}'] = counts.get(i, 0)
    return result


def process_dual_channel_sample(ims_path, mask_path, detect_params=None,
                                threshold_ch0=None, threshold_ch1=None,
                                verbose=True, spot_output_dir=None):
    """
    Process a single sample for dual-channel analysis.
    CH0 and CH1 are both RNA FISH channels; CH2 is a structural channel (mask).

    Args:
        ims_path: Path to IMS file.
        mask_path: Path to corresponding mask NRRD.
        detect_params: Common detection params (dict).
        threshold_ch0: CH0 threshold override.
        threshold_ch1: CH1 threshold override.
        verbose: Whether to print detailed info.

    Returns:
        dict: Combined CH0/CH1 results, or None on failure.
    """
    if detect_params is None:
        detect_params = {
            'z_radius_nm': 500,
            'xy_radius_nm': 96,
            'threshold': None,
            'voxel_size_nm': (500, 96, 96),
            'decompose': True,
            'cluster': False,
            'verbose': verbose
        }

    filename = Path(ims_path).stem
    lp_info = extract_lp_info(filename)

    try:
        if verbose:
            print(f"Processing: {filename}")

        # Load mask (from CH2)
        mask_array = load_mask(mask_path)

        # Process CH0
        if verbose:
            print("  Processing CH0...")
        rna_ch0 = load_rna_channel(ims_path, channel_idx=0)
        params_ch0 = detect_params.copy()
        if threshold_ch0 is not None:
            params_ch0['threshold'] = threshold_ch0
        result_ch0 = detect_spots_for_channel(rna_ch0, mask_array, params_ch0, 'CH0', spot_output_dir)
        del rna_ch0
        gc.collect()

        # Process CH1
        if verbose:
            print("  Processing CH1...")
        rna_ch1 = load_rna_channel(ims_path, channel_idx=1)
        params_ch1 = detect_params.copy()
        if threshold_ch1 is not None:
            params_ch1['threshold'] = threshold_ch1
        result_ch1 = detect_spots_for_channel(rna_ch1, mask_array, params_ch1, 'CH1', spot_output_dir)
        del rna_ch1
        gc.collect()

        # Assemble result
        result = {'filename': filename}
        if lp_info:
            result.update(lp_info)
        else:
            result.update({'zt_range': None, 'zt_start': None, 'zt_end': None,
                          'condition': None, 'exp_group': None})

        # CH0 results
        result['ch0_total_spots'] = result_ch0['total']
        for i in range(1, 6):
            result[f'ch0_mask_{i}'] = result_ch0[f'mask_{i}']

        # CH1 results
        result['ch1_total_spots'] = result_ch1['total']
        for i in range(1, 6):
            result[f'ch1_mask_{i}'] = result_ch1[f'mask_{i}']

        if verbose:
            print(f"  {lp_info['exp_group'] if lp_info else '?'} | "
                  f"CH0={result_ch0['total']} spots, CH1={result_ch1['total']} spots")

        return result

    except Exception as e:
        logging.error(f"Error processing {filename}: {str(e)}")
        return None


def calculate_dual_channel_statistics(df):
    """
    Calculate dual-channel statistics grouped by zt_range and exp_group.

    Computes sLNv (mask1-4 mean) and 5th sLNv (mask5) for both CH0 and CH1.

    Args:
        df: Dual-channel results DataFrame.

    Returns:
        (stats_ch0, stats_ch1): Tuple of DataFrames.
    """
    results = {}
    for ch_name in ['ch0', 'ch1']:
        df_ch = df.copy()
        mask_cols = [f'{ch_name}_mask_{i}' for i in range(1, 5)]
        df_ch['sLNv'] = df_ch[mask_cols].mean(axis=1)
        df_ch['sLNv_5th'] = df_ch[f'{ch_name}_mask_5']

        stats = df_ch.groupby(['zt_range', 'exp_group', 'condition']).agg(
            sLNv_mean=('sLNv', 'mean'),
            sLNv_std=('sLNv', 'std'),
            sLNv_sem=('sLNv', 'sem'),
            sLNv_n=('sLNv', 'count'),
            sLNv_5th_mean=('sLNv_5th', 'mean'),
            sLNv_5th_std=('sLNv_5th', 'std'),
            sLNv_5th_sem=('sLNv_5th', 'sem'),
            sLNv_5th_n=('sLNv_5th', 'count'),
        ).reset_index()
        results[ch_name] = stats

    return results['ch0'], results['ch1']


def plot_dual_channel_bar(stats_ch0, stats_ch1, output_path=None, figsize=(10, 5)):
    """
    Dual-channel grouped bar chart.

    X axis = condition; grouped bars for sLNv / 5th sLNv; colored by Light/Dark.

    Args:
        stats_ch0, stats_ch1: Statistics DataFrames for CH0 and CH1.
        output_path: Path to save figure.
        figsize: Figure size.

    Returns:
        matplotlib Figure.
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=False)
    channel_data = [('CH0', stats_ch0, axes[0]), ('CH1', stats_ch1, axes[1])]

    for ch_label, stats, ax in channel_data:
        # Build display label
        stats['label'] = stats['condition'] + '\n(' + stats['exp_group'] + ')'
        x = np.arange(len(stats))
        width = 0.35

        ax.bar(x - width/2, stats['sLNv_mean'], width,
               yerr=stats['sLNv_sem'], capsize=4,
               color='#2E86AB', label='sLNv', edgecolor='white')
        ax.bar(x + width/2, stats['sLNv_5th_mean'], width,
               yerr=stats['sLNv_5th_sem'], capsize=4,
               color='#E94F37', label='5th sLNv', edgecolor='white')

        ax.set_title(ch_label, fontweight='bold', fontsize=14)
        ax.set_xticks(x)
        ax.set_xticklabels(stats['label'], fontsize=10)
        ax.set_ylabel('RNA Spots per Cell', fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(which='major', width=1.5, length=5)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')

        if ch_label == 'CH0':
            ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1),
                      framealpha=0.9, prop={'weight': 'bold'})

    fig.tight_layout()
    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {output_path}")
    plt.show()
    return fig
