"""LD4 and LP batch execution around the original direct-detection functions."""
from __future__ import annotations

from datetime import datetime
import hashlib
from importlib.metadata import version
import json
import logging
import os
from pathlib import Path
import platform
import re
import uuid

import numpy as np
import pandas as pd
import rnascope_quant as core


def _resolve(base, value):
    path = Path(value).expanduser()
    return (path if path.is_absolute() else base / path).resolve()


def load_experiment(config_path, *, data_root=None, results_root=None):
    """Paths in JSON are relative to that JSON, independent of notebook cwd."""
    config_path = Path(config_path).resolve()
    cfg = json.loads(config_path.read_text(encoding="utf-8"))
    expected_mode = {"LD4": "single", "LP_dual": "dual"}
    if cfg.get("experiment") not in expected_mode or cfg.get("mode") != expected_mode[cfg["experiment"]]:
        raise ValueError("Choose the LD4 single-channel or LP_dual configuration")
    datasets = cfg.get("datasets", [])
    ids = [item["id"] for item in datasets]
    if not ids or len(set(ids)) != len(ids):
        raise ValueError("Dataset IDs must be present and unique")
    cfg["data_root"] = _resolve(config_path.parent, data_root or cfg["data_root"])
    cfg["results_root"] = _resolve(config_path.parent, results_root or cfg["results_root"])
    cfg["config_file"] = config_path
    params = cfg["detection"]
    params["voxel_size_nm"] = tuple(params["voxel_size_nm"])
    if len(params["voxel_size_nm"]) != 3 or min(params["voxel_size_nm"]) <= 0:
        raise ValueError("voxel_size_nm must contain positive Z, Y, X spacings")
    for exclusion in cfg.get("exclusions", []):
        if not exclusion.get("reason", "").strip():
            raise ValueError("Every exclusion rule needs a reason")
        re.compile(exclusion["pattern"])
    return cfg


def plan_experiment(cfg):
    """Pair masks and IMS files, keeping the original dataset-wide suffix priority.

    LD4 uses _CH1_mask. LP uses _CH3_mask if any exist in that dataset,
    otherwise _CH2_mask, as in the supplied notebook's original module.
    """
    records = []
    seen_sources = set()
    for dataset in cfg["datasets"]:
        raw = _resolve(cfg["data_root"], dataset["ims_dir"])
        masks = _resolve(cfg["data_root"], dataset["mask_dir"])
        base = {"dataset_id": dataset["id"], "filename": "", "sample_id": "",
                "ims_path": str(raw), "mask_path": str(masks)}
        if not raw.is_dir() or not masks.is_dir():
            records.append(dict(base, status="missing_directory"))
            continue
        selected = []
        for suffix in cfg["mask_suffixes"]:
            selected = sorted(masks.glob(f"*{suffix}.nrrd"))
            if selected:
                break
        if not selected:
            records.append(dict(base, status="no_masks"))
        for mask in selected:
            stem = mask.stem[:-len(suffix)]
            image = raw / f"{stem}.ims"
            key = str(image.resolve()).casefold()
            if key in seen_sources:
                raise ValueError(f"The same IMS file was selected twice: {image}")
            seen_sources.add(key)
            records.append({"dataset_id": dataset["id"], "filename": stem,
                            "sample_id": f"{dataset['id']}::{stem}",
                            "ims_path": str(image), "mask_path": str(mask),
                            "status": "ready" if image.is_file() else "missing_ims"})
    return pd.DataFrame(records)


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()


def _software():
    return {"python": platform.python_version(), **{
        name: version(name) for name in ("numpy", "scipy", "scikit-image", "scikit-learn",
                                         "big-fish", "SimpleITK", "h5py", "tifffile",
                                         "imaris-ims-file-reader", "pandas", "openpyxl")
    }}


def _parameters(cfg):
    # Output settings and exclusions do not affect detection, so changing
    # those does not require re-detecting already validated samples.
    return {"mode": cfg["mode"], "detection": cfg["detection"],
            "threshold_ch0": cfg.get("threshold_ch0"),
            "threshold_ch1": cfg.get("threshold_ch1"),
            "random_seed": cfg.get("random_seed"),
            "save_spots": cfg.get("save_spots", True),
            "core_sha256": hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest(),
            "workflow_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "software": _software()}


def _signature(filename):
    path = Path(filename)
    info = path.stat()
    return {"path": str(path.resolve()), "size": info.st_size, "mtime_ns": info.st_mtime_ns}


def _atomic_json(path, content):
    temporary = path.with_name(f"{path.name}.{uuid.uuid4().hex[:8]}.tmp")
    try:
        temporary.write_text(json.dumps(content, indent=2, default=str), encoding="utf-8")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _rng_state():
    name, state, position, gaussian, cached = np.random.get_state()
    return [name, state.tolist(), position, gaussian, cached]


def _run_sample(cfg, row, spot_dir):
    if cfg["mode"] == "single":
        if core.extract_zt_time(row.filename) is None:
            raise ValueError("Filename does not contain a ZT time")
        result = core.process_single_sample(
            row.ims_path, row.mask_path, cfg["detection"], verbose=False,
            spot_output_dir=spot_dir,
        )
    else:
        if core.extract_lp_info(row.filename) is None:
            raise ValueError("Filename does not contain LP_ZTstart_end_condition metadata")
        result = core.process_dual_channel_sample(
            row.ims_path, row.mask_path, cfg["detection"],
            threshold_ch0=cfg["threshold_ch0"], threshold_ch1=cfg["threshold_ch1"],
            verbose=False, spot_output_dir=spot_dir,
        )
    if result is None:
        raise RuntimeError("Original quantification function failed; see analysis.log for the cause")
    return dict(result, sample_id=row.sample_id, dataset_id=row.dataset_id)


def filter_samples(raw_data, exclusions):
    """Apply documented case-insensitive regex rules and retain all reasons."""
    reasons = pd.Series("", index=raw_data.index, dtype=object)
    for rule in exclusions:
        matched = raw_data["filename"].str.contains(rule["pattern"], case=False, regex=True, na=False)
        addition = f"{rule['pattern']}: {rule['reason']}"
        reasons.loc[matched] = reasons.loc[matched].map(lambda old: f"{old}; {addition}" if old else addition)
    removed = raw_data.loc[reasons != ""].copy()
    removed["exclusion_reason"] = reasons.loc[removed.index]
    return raw_data.loc[reasons == ""].copy(), removed


def _export(cfg, raw, analysis, removed, spot_index, output):
    tables = {"Raw_Data": raw, "Analysis_Data": analysis, "Removed_Samples": removed,
              "Spot_Files": spot_index}
    if cfg["mode"] == "single":
        stats = core.calculate_cell_type_statistics(analysis.copy())
        tables["Statistics"] = stats
        long = analysis.melt(id_vars=["sample_id", "dataset_id", "filename", "zt_time", "total_spots"],
                             value_vars=[f"mask_{i}" for i in range(1, 6)],
                             var_name="mask_id", value_name="spot_count")
        long["cell_type"] = long["mask_id"].map(lambda label: "sLNv_5th" if label == "mask_5" else "sLNv")
        tables["Long_Format"] = long
    else:
        stats = core.calculate_dual_channel_statistics(analysis.copy())
        tables["Statistics_CH0"], tables["Statistics_CH1"] = stats
    parameters = dict(cfg["detection"])
    if cfg["mode"] == "dual":
        parameters.update(threshold_ch0=cfg["threshold_ch0"], threshold_ch1=cfg["threshold_ch1"])
    parameters.update(random_seed=cfg.get("random_seed"),
                      mask_counting="Original labels 1-5; absent labels are counted as zero")
    tables["Detection_Params"] = pd.DataFrame([{"key": key, "value": str(value)} for key, value in parameters.items()])
    workbook = output / f"{cfg['experiment']}.xlsx"
    # Keep the original workflow's pandas/openpyxl report format.
    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        for sheet, frame in tables.items():
            frame.to_excel(writer, sheet_name=sheet, index=False)
    figure = None
    if cfg.get("make_plots", True) and not analysis.empty:
        figure = output / f"{cfg['experiment']}.png"
        if cfg["mode"] == "single":
            core.plot_expression_curves(stats, output_path=str(figure))
        else:
            core.plot_dual_channel_bar(*stats, output_path=str(figure))
    return stats, workbook, figure


def run_quantification(cfg, *, reuse_cache=True):
    """Run only ready samples, with input/parameter/code-aware checkpoints.

    Cached results require matching file size/mtime, detection parameters,
    software versions, and code. Legacy filename-only CSVs are never loaded.
    Original mask-1..5 and absent-label-as-zero semantics are preserved.
    """
    plan = plan_experiment(cfg)
    if plan.empty or not (plan["status"] == "ready").all():
        raise ValueError("Resolve missing input directories/files or masks before running; inspect plan_experiment(cfg)")
    root = cfg["results_root"] / cfg["experiment"]
    cache = root / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    output = root / datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    output.mkdir(parents=True, exist_ok=False)
    parameters = _parameters(cfg)
    snapshot = dict(cfg, effective_detection=parameters, numpy_random_state_before=_rng_state())
    _atomic_json(output / "configuration.json", snapshot)
    plan.to_csv(output / "input_manifest.csv", index=False)
    results, failures, index = [], [], []
    hits = 0
    channels = ["CH0"] if cfg["mode"] == "single" else ["CH0", "CH1"]
    handler = logging.FileHandler(output / "analysis.log", encoding="utf-8")
    logging.getLogger().addHandler(handler)
    try:
        for row in plan.itertuples(index=False):
            fingerprint = {"sample_id": row.sample_id, "parameters": parameters,
                           "ims": _signature(row.ims_path), "mask": _signature(row.mask_path)}
            key = _digest(fingerprint)
            # Short names keep deeply nested Windows publication folders usable.
            # The full digest is still validated inside each checkpoint.
            checkpoint = cache / f"{key[:32]}.json"
            spots = cache / key[:32] if cfg.get("save_spots", True) else None
            saved = None
            if reuse_cache and checkpoint.is_file():
                try:
                    saved = json.loads(checkpoint.read_text(encoding="utf-8"))
                    if saved["cache_key"] != key or not saved.get("result"):
                        saved = None
                    if spots is not None and not all((spots / f"{ch}_spots.npz").is_file() for ch in channels):
                        saved = None
                except (ValueError, KeyError, TypeError, OSError):
                    saved = None
            cache_hit = saved is not None
            try:
                if cache_hit:
                    result = saved["result"]
                    hits += 1
                else:
                    previous_state = np.random.get_state()
                    seed = cfg.get("random_seed")
                    if seed is not None:
                        # Stable per sample, including when some samples hit cache.
                        np.random.seed(int(_digest([seed, row.sample_id])[:8], 16))
                    before = _rng_state()
                    try:
                        result = _run_sample(cfg, row, spots)
                    finally:
                        if seed is not None:
                            np.random.set_state(previous_state)
                    _atomic_json(checkpoint, {"cache_key": key, "fingerprint": fingerprint,
                                              "numpy_random_state_before": before, "result": result})
                results.append(result)
                if spots is not None:
                    for channel in channels:
                        index.append({"sample_id": row.sample_id, "filename": row.filename,
                                      "channel": channel, "cache_hit": cache_hit,
                                      "spot_file": os.path.relpath(spots / f"{channel}_spots.npz", output)})
            except Exception as exception:
                logging.exception("Quantification failed: %s", row.sample_id)
                failures.append({"sample_id": row.sample_id, "filename": row.filename, "error": str(exception)})
    finally:
        logging.getLogger().removeHandler(handler)
        handler.close()
    pd.DataFrame(failures, columns=["sample_id", "filename", "error"]).to_csv(output / "failures.csv", index=False)
    if not results:
        raise RuntimeError(f"No samples succeeded. See {output / 'failures.csv'}")
    raw = pd.DataFrame(results)
    order = ["zt_time", "filename"] if cfg["mode"] == "single" else ["zt_start", "exp_group", "filename"]
    raw = raw.sort_values(order).reset_index(drop=True)
    analysis, removed = filter_samples(raw, cfg.get("exclusions", []))
    spot_index = pd.DataFrame(index, columns=["sample_id", "filename", "channel", "cache_hit", "spot_file"])
    stats, workbook, figure = _export(cfg, raw, analysis, removed, spot_index, output)
    return {"raw_data": raw, "analysis_data": analysis, "removed_samples": removed,
            "statistics": stats, "workbook": workbook, "figure": figure,
            "output_dir": output, "spot_index": spot_index,
            "cache_hits": hits, "failures": pd.DataFrame(failures)}


def export_spot_markers(spot_file, output_tiff):
    """Render the saved final coordinates from a run; do not re-detect spots."""
    with np.load(spot_file, allow_pickle=False) as saved:
        return core.visualize_spots_as_image(saved["final_spots"], tuple(saved["image_shape_zyx"]), str(output_tiff))


def preprocess_dataset(cfg, dataset_id, *, channel_index, p_low=0.5, p_high=99.995):
    """Export a normalized helper channel for external manual mask annotation.

    Spot quantification always reads the unnormalized IMS. This operation
    does not create cell labels or replace the externally annotated masks.
    """
    dataset = next(item for item in cfg["datasets"] if item["id"] == dataset_id)
    source = _resolve(cfg["data_root"], dataset["ims_dir"])
    if not source.is_dir() or not any(source.glob("*.ims")):
        raise FileNotFoundError(f"No IMS files in {source}")
    output = cfg["results_root"] / "preprocessing" / dataset_id
    core.batch_process_ims(source, output, channel_idx=channel_index, p_low=p_low, p_high=p_high)
    return output
