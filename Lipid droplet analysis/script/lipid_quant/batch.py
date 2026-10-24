"""Notebook-safe Windows process workers for multi-directory IMS analysis."""
from __future__ import annotations

import copy
import hashlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import pandas as pd

from .paths import PathInput, discover_ims_files


def batch_analysis_config(config: dict) -> dict:
    result = copy.deepcopy(config)
    if 'exclude_optic_lobes' in result.get('batch', {}):
        result['roi']['exclude_optic_lobes']=result['batch']['exclude_optic_lobes']
    for key in ("save_qc", "save_roi_qc", "save_3d_marker_tiff", "save_mask_tiff",
                "save_roi_mask_tiff", "print_plane_progress"):
        if key in result.get("batch", {}):
            result.setdefault("output", {})[key] = result["batch"][key]
    return result


def _batch_worker(path: Path, config: dict, output_root: Path, run_name: str) -> Path:
    from .pipeline import analyze_ims
    return analyze_ims(path, config, output_root, run_name=run_name)


def _prepare_batch_samples(files: list[Path], config: dict) -> list[tuple[Path, dict | None, str]]:
    """Resolve per-sample guides in input order; only missing guides are skipped."""
    from .anatomical_roi import MissingGuideError, resolve_sample_guide
    roi_cfg = config["roi"]
    if len(files) > 1 and roi_cfg["mode"] == "guided_central_brain":
        path = str(roi_cfg.get("guide_mask_path") or "")
        if not roi_cfg.get("guide_manifest_path") and "{stem}" not in path:
            raise ValueError("Guided batches require guide_manifest_path or guide_mask_path containing {stem}")
    prepared = []
    for path in files:
        try:
            prepared.append((path, resolve_sample_guide(config, path), ""))
        except MissingGuideError as error:
            prepared.append((path, None, str(error)))
    if roi_cfg["mode"] == "guided_central_brain" and len(files) > 1:
        guide_paths = [str(c["roi"]["guide_mask_path"]).casefold()
                       for _, c, _ in prepared if c is not None]
        if len(set(guide_paths)) != len(guide_paths):
            raise ValueError("Multiple IMS resolve to the same guide. Use separate per-brain masks, {parent}/{stem}, or a full-path manifest")
    return prepared


def run_batch_analysis(inputs: PathInput, config: dict, output_root: str | Path = "results", *,
                       workers: int | None = None, pattern: str = "*.ims", recursive: bool = False,
                       run: bool = True, raise_on_error: bool = True,
                       exclude_optic_lobes: bool | None = None) -> pd.DataFrame:
    """One shared worker pool across inputs. ``run=False`` returns a dry-run plan.

    Path duplicates are removed; different files with the same basename retain
    separate source hashes and output directories. Each guided sample's mask is
    resolved before dispatch. Missing manual guides are recorded as
    ``skipped_missing_roi`` (with ``skip_reason``), never replaced by auto ROI.
    Invalid guides or other configuration errors are not treated as missing.
    """
    from .anatomical_roi import MissingGuideError
    from .pipeline import _config_hash
    from .anatomical_v2 import set_optic_lobe_exclusion
    files = discover_ims_files(inputs, pattern, recursive=recursive)
    cfg = set_optic_lobe_exclusion(batch_analysis_config(config), exclude_optic_lobes)
    requested = cfg.get("batch", {}).get("workers", 1) if workers is None else workers
    if isinstance(requested, bool) or int(requested) != requested or requested < 1:
        raise ValueError("workers must be a positive integer")
    roi_cfg = cfg["roi"]
    prepared = _prepare_batch_samples(files, cfg)
    root = Path(output_root).resolve()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    jobs = []
    records = []
    for index, (path, sample_cfg, skip_reason) in enumerate(prepared):
        source_hash = hashlib.sha256(str(path).casefold().encode()).hexdigest()[:8]
        config_hash = _config_hash(sample_cfg) if sample_cfg is not None else ""
        name = f"{path.stem}_{timestamp}_{source_hash}_{config_hash}"
        record_cfg = sample_cfg if sample_cfg is not None else cfg
        records.append({"input_file": str(path), "input_directory": str(path.parent),
                        "output_directory": str(root / name) if sample_cfg is not None else "",
                        "status": "planned" if sample_cfg is not None else "skipped_missing_roi",
                        "error": "", "skip_reason": skip_reason,
                        "config_hash": config_hash,
                        "roi_method": record_cfg["roi"].get("central_method"),
                        "exclude_optic_lobes": record_cfg['roi'].get('exclude_optic_lobes',True),
                        "roi_scope": 'whole_brain' if not record_cfg['roi'].get('exclude_optic_lobes',True) else 'central_brain' if record_cfg['roi']['mode'] in {'auto_central_brain','guided_central_brain'} else record_cfg['roi']['mode'],
                        "tissue_margin_um_v1": record_cfg["roi"].get("anatomical_tissue_margin_um"),
                        "final_margin_xy_um_v2": record_cfg["roi"].get("v2_final_margin_xy_um"),
                        "save_roi_qc": record_cfg.get("output", {}).get("save_roi_qc", False),
                        "save_roi_mask_tiff": record_cfg.get("output", {}).get("save_roi_mask_tiff", False)})
        if sample_cfg is not None:
            jobs.append((index, (path, sample_cfg, root, name)))
    worker_count = min(int(requested), len(jobs))
    plan = pd.DataFrame(records)
    print(f"Batch: {len(files)} unique IMS file(s), {len(jobs)} eligible, "
          f"{len(files)-len(jobs)} skipped (missing ROI), {worker_count} worker process(es)", flush=True)
    for record in records:
        if record["status"] == "skipped_missing_roi":
            print(f"skipped_missing_roi: {Path(record['input_file']).name} ({record['skip_reason']})", flush=True)
    print(f"Effective ROI: method={roi_cfg.get('central_method')}; "
          f"exclude optic lobes={roi_cfg.get('exclude_optic_lobes',True)}; "
          f"v1 tissue margin={roi_cfg.get('anatomical_tissue_margin_um')} (ignored by v2); "
          f"v2 final XY margin={roi_cfg.get('v2_final_margin_xy_um')}; "
          f"ROI QC={cfg.get('output', {}).get('save_roi_qc', False)}", flush=True)
    if not run:
        return plan
    root.mkdir(parents=True, exist_ok=True)
    import yaml
    with (root / f"batch_config_{timestamp}.yaml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump(cfg, handle, sort_keys=False, allow_unicode=True)
    manifest = root / f"batch_manifest_{timestamp}.csv"
    plan.to_csv(manifest, index=False)

    def record_finished(index, result=None, error=None):
        if isinstance(error, MissingGuideError):
            records[index].update(status="skipped_missing_roi", skip_reason=str(error),
                                  error="", output_directory="")
        else:
            records[index]["status"] = "failed" if error is not None else "completed"
            records[index]["error"] = str(error) if error is not None else ""
        if result is not None:
            records[index]["output_directory"] = str(result)
        pd.DataFrame(records).to_csv(manifest, index=False)
        print(f"{records[index]['status']}: {files[index].name}", flush=True)

    if worker_count == 1:
        for index, job in jobs:
            try:
                record_finished(index, result=_batch_worker(*job))
            except Exception as error:
                record_finished(index, error=error)
    elif worker_count > 1:
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            pending = {executor.submit(_batch_worker, *job): index for index, job in jobs}
            for future in as_completed(pending):
                index = pending[future]
                try:
                    record_finished(index, result=future.result())
                except Exception as error:
                    record_finished(index, error=error)
    result = pd.DataFrame(records)
    failed = result.loc[result["status"] == "failed"]
    if raise_on_error and not failed.empty:
        raise RuntimeError(f"{len(failed)} batch file(s) failed. See {manifest}\n"
                           + failed[["input_file", "error"]].to_string(index=False))
    return result
