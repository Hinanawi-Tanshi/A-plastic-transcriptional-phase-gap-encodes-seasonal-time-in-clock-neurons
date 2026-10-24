"""Configuration and orchestration for the publication lipid-droplet notebook."""
from __future__ import annotations

import copy
from datetime import datetime
import json
import os
from pathlib import Path
import platform

import pandas as pd

from lipid_quant import (
    export_guide_projections, load_config, plan_guide_projections,
    run_batch_analysis, run_posthoc_qc,
)
from lipid_quant.batch_summary import (
    generate_standard_report, load_batch_results, summarize_brains,
    validate_batch_results,
)


def _resolve(base: Path, value: str | Path) -> Path:
    path = Path(value).expanduser()
    return (path if path.is_absolute() else base / path).resolve()


def load_workflow(config_path, *, data_root=None, results_root=None):
    """Resolve relative paths against the JSON file, without creating outputs.

    data_root can point to an external raw-data location. All workflow stages
    use the same input list; stage-specific batch paths are not retained.
    """
    config_path = Path(config_path).resolve()
    settings = json.loads(config_path.read_text(encoding="utf-8"))
    if not settings.get("input_paths"):
        raise ValueError("input_paths must contain at least one IMS file or directory")
    settings["data_root"] = _resolve(config_path.parent, data_root or settings["data_root"])
    settings["results_root"] = _resolve(config_path.parent, results_root or settings["results_root"])
    settings["algorithm_config"] = _resolve(config_path.parent, settings["algorithm_config"])
    settings["inputs"] = [_resolve(settings["data_root"], item) for item in settings["input_paths"]]
    settings["existing_result_dirs"] = [
        _resolve(settings["results_root"], item) for item in settings.get("existing_result_dirs", [])
    ]
    settings["config_file"] = config_path
    return settings


def analysis_config(settings):
    """Apply the original notebook's guided-CB and batch-output overrides."""
    cfg = load_config(settings["algorithm_config"])
    if settings.get("use_manual_cb_roi", True):
        if settings.get("exclude_optic_lobes", True) is not True:
            raise ValueError("Manual central-brain guides require exclude_optic_lobes=true")
        if cfg["input"]["resolution_level"] != 0:
            raise ValueError("Guided ROI analysis requires original XY resolution (level 0)")
        cfg["roi"].update(
            mode="guided_central_brain", central_method="anatomical_v2",
            guide_manifest_path="{source_dir}/central_brain_roi/guide_manifest.json",
            guide_mask_path=None, guide_z_index=None, guide_allow_resize=False,
        )
    cfg.setdefault("batch", {}).update(
        workers=settings.get("workers", 2), save_qc=False, save_roi_qc=False,
        save_3d_marker_tiff=False, save_mask_tiff=False,
        save_roi_mask_tiff=settings.get("export_final_cb_mask", False),
        print_plane_progress=False,
    )
    return cfg


def _missing_input_plan(settings):
    records = [{"input_file": str(path), "status": "missing_input",
                "detail": "Set data_root/input_paths to your IMS data before running"}
               for path in settings["inputs"] if not path.exists()]
    return pd.DataFrame(records)


def prepare_guides(settings, *, run=False):
    """Plan/export raw MIPs next to the IMS files; Fiji ROI drawing is manual.

    Existing ROI files are preserved. This stage intentionally writes
    central_brain_roi/ alongside each input when run=True. Quantification
    subsequently reads these annotations without modifying them.
    """
    kwargs = {"recursive": settings.get("recursive", False)}
    if not run:
        missing = _missing_input_plan(settings)
        if not missing.empty:
            return missing
        return plan_guide_projections(settings["inputs"], **kwargs)
    channel = load_config(settings["algorithm_config"])["input"]
    return export_guide_projections(
        settings["inputs"], **kwargs,
        time_index=channel["time_index"], channel_index=channel["channel_index"],
        overwrite_projections=False,
    )


def _new_output(root: Path, stage: str) -> Path:
    path = root / stage / datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    # Per-brain names retain the original specimen ID and provenance. Use
    # native extended paths on Windows instead of changing those identifiers.
    if os.name == "nt":
        absolute = str(path.resolve())
        if not absolute.startswith("\\\\?\\"):
            absolute = ("\\\\?\\UNC\\" + absolute[2:] if absolute.startswith("\\\\")
                        else "\\\\?\\" + absolute)
        path = Path(absolute)
    path.mkdir(parents=True, exist_ok=False)
    return path


def quantify(settings, *, run=False):
    """Plan or run a batch in a new output directory; never mix repeated runs."""
    cfg = analysis_config(settings)
    output = settings["results_root"] / "batches"
    if not run:
        missing = _missing_input_plan(settings)
        if not missing.empty:
            return missing
    if run:
        # Validate discovery and annotation configuration before creating output.
        run_batch_analysis(
            settings["inputs"], cfg, output, workers=settings.get("workers", 2),
            recursive=settings.get("recursive", False),
            exclude_optic_lobes=settings.get("exclude_optic_lobes", True), run=False,
        )
        output = _new_output(settings["results_root"], "batches")
        snapshot = dict(settings, python=platform.python_version())
        (output / "workflow_config.json").write_text(
            json.dumps(snapshot, default=str, indent=2), encoding="utf-8"
        )
    return run_batch_analysis(
        settings["inputs"], cfg, output, workers=settings.get("workers", 2),
        recursive=settings.get("recursive", False),
        exclude_optic_lobes=settings.get("exclude_optic_lobes", True), run=run,
    )


def completed_result_dirs(batch_status):
    """Select the output directories from this batch, excluding skipped samples."""
    completed = batch_status.loc[batch_status["status"] == "completed", "output_directory"]
    return [Path(value) for value in completed if str(value)]


def _selected_dirs(settings, result_dirs):
    selected = settings["existing_result_dirs"] if result_dirs is None else list(result_dirs)
    if not selected:
        raise ValueError("Choose completed result directories, or set existing_result_dirs in workflow.json")
    return selected


def supplement_qc(settings, result_dirs=None, *, run=False):
    """Generate QC from the saved detections; no droplet re-detection occurs."""
    options = settings.get("qc", {})
    return run_posthoc_qc(
        _selected_dirs(settings, result_dirs),
        regex=options.get("regex"), match_field="brain_id", recursive=False,
        save_qc=options.get("save_images", True),
        save_3d_marker_tiff=options.get("save_3d_marker", True),
        qc_slices=options.get("slices", 6), overwrite=False, run=run,
    )


def summarize(settings, result_dirs=None, *, show=True):
    """Summarize per brain and keep the original report/statistical defaults."""
    results = load_batch_results(
        _selected_dirs(settings, result_dirs),
        genotype_labels=settings["genotype_labels"], group_order=settings["group_order"],
        include_unlisted_groups=settings.get("include_unlisted_groups", False),
        load_objects=True, accepted_objects_only=True, load_planes=True,
        load_detections=False, duplicate_policy="error",
    )
    validation = validate_batch_results(results)
    summary = summarize_brains(results.brains)
    output = _new_output(settings["results_root"], "reports")
    options = copy.deepcopy(settings.get("report", {}))
    options["show"] = show
    report = generate_standard_report(results, output, **options)
    return {"results": results, "validation": validation, "summary": summary,
            "report": report, "output_dir": output}
