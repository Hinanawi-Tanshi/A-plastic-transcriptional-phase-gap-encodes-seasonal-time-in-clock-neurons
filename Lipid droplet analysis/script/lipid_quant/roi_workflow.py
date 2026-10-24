"""Notebook controls for versioned ROI experiments and explicit batch selection."""
import copy
import hashlib
import json
from datetime import datetime
from pathlib import Path
import yaml

from .anatomical_v2 import V2_DEFAULTS, resolve_v2_config, set_optic_lobe_exclusion


def make_roi_v2_config(base_config, *, final_margin_xy_um=2.0, final_margin_z_um=0.0,
                       save_roi_qc=True, review_policy='warn', roi_overrides=None,
                       exclude_optic_lobes=None):
    """Copy all detection settings unchanged; explicitly opt into v2 ROI only."""
    cfg=copy.deepcopy(base_config)
    cfg['roi'].update({**V2_DEFAULTS,**{k:v for k,v in cfg['roi'].items() if k in V2_DEFAULTS}})
    cfg['roi'].update(central_method='anatomical_v2',v2_final_margin_xy_um=final_margin_xy_um,
                      v2_final_margin_z_um=final_margin_z_um,v2_review_policy=review_policy)
    if roi_overrides: cfg['roi'].update(roi_overrides)
    cfg.setdefault('output',{})['save_roi_qc']=bool(save_roi_qc)
    cfg.setdefault('batch',{})['save_roi_qc']=bool(save_roi_qc)
    return set_optic_lobe_exclusion(cfg,exclude_optic_lobes)


def select_batch_config(source, *, base=None, trial=None, guided=None):
    """No silent fallback when a requested notebook configuration is absent."""
    choices={'base':base,'trial':trial,'guided':guided}
    if source not in choices: raise ValueError('Config source must be base, trial or guided')
    if choices[source] is None: raise ValueError(f'{source} config is unavailable; run its parameter cell first')
    cfg=resolve_v2_config(copy.deepcopy(choices[source]))
    print(f'Batch configuration source: {source}')
    print(yaml.safe_dump({'roi':cfg['roi'],'output':cfg.get('output',{}),'batch':cfg.get('batch',{})},sort_keys=False,allow_unicode=True))
    return cfg


def save_analysis_config(config, path, *, overwrite=False):
    cfg=resolve_v2_config(config)
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w' if overwrite else 'x',encoding='utf-8') as f:
        yaml.safe_dump(cfg,f,sort_keys=False,allow_unicode=True)
    return path


def preview_roi_batch(inputs, config, output_root, *, recursive=False, z_indices=None,
                      exclude_optic_lobes=None):
    """Serial ROI-only validation across directories. No LD run, no native TIFF.

    Timestamped session and source hashes prevent collisions/overwriting.
    Returns plain records, convenient for a notebook DataFrame if desired.
    Missing manual guides are recorded as skipped, just as in LD batches.
    """
    from .paths import discover_ims_files
    from .roi_preview import preview_roi
    from .anatomical_roi import MissingGuideError
    files=discover_ims_files(inputs,recursive=recursive)
    from .batch import _prepare_batch_samples, batch_analysis_config
    cfg=set_optic_lobe_exclusion(batch_analysis_config(config),exclude_optic_lobes)
    prepared=_prepare_batch_samples(files,cfg)
    root=Path(output_root)/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    root.mkdir(parents=True,exist_ok=False)
    results=[]
    for p,sample_cfg,skip_reason in prepared:
        key=hashlib.sha256(str(p).casefold().encode()).hexdigest()[:8]
        directory=root/(p.stem+'_'+key)
        record={'input_file':str(p),'output_directory':str(directory) if sample_cfg is not None else '',
                'skip_reason':skip_reason}
        try:
            if sample_cfg is None:
                record.update(status='skipped_missing_roi',roi_status='not_run',error='')
            else:
                paths=preview_roi(p,sample_cfg,directory,z_indices=z_indices,export_native_mask=False)
                info=json.loads(paths['metadata'].read_text(encoding='utf-8'))
                area=info.get('central_area_by_z_um2',[])
                dz=info.get('low_spacing_zyx_um',[0])[0]
                record.update(status='completed',roi_status=info.get('roi_status','not_recorded'),
                              roi_scope=info.get('roi_scope',info.get('mode','unknown')),
                              exclude_optic_lobes=info.get('exclude_optic_lobes',True),
                              roi_volume_um3=float(sum(area)*dz),review_flags=info.get('review_flags',[]),error='')
        except MissingGuideError as exc:
            record.update(status='skipped_missing_roi',roi_status='not_run',error='',
                          skip_reason=str(exc),output_directory='')
        except Exception as exc:
            record.update(status='failed',roi_status='failed',error=str(exc))
        results.append(record)
        (root/'roi_preview_manifest.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
        print(f'{record["status"]}: {p.name}',flush=True)
    return results


def display_roi_preview(paths):
    """Notebook rendering without plt.show()/interactive-backend requirements."""
    from IPython.display import Image, display
    for key in ('slices','structural_slices','structural_depth'):
        if key in paths: display(Image(filename=str(paths[key])))
    if 'metadata' in paths:
        info=json.loads(Path(paths['metadata']).read_text(encoding='utf-8'))
        print('ROI status:',info.get('roi_status','not_recorded'))
        for flag in info.get('review_flags',[]): print('Review:',flag)


def preview_margin_sweep(ims_path, config, output_root, *, radii_um=(0,2,4,6)):
    """Build v2 identity once and compare final XY margins on that fixed base."""
    import numpy as np
    from imaris_ims_file_reader.ims import ims
    from .anatomical_roi import resolve_sample_guide
    from .anatomical_v2 import build_anatomical_v2_roi, expand_final_roi
    from .roi_qc import save_structural_roi_qc
    cfg=resolve_v2_config(config)
    if cfg['roi'].get('central_method')!='anatomical_v2': raise ValueError('Margin sweep requires anatomical_v2')
    radii=sorted(set(float(r) for r in radii_um))
    if not radii or not np.isfinite(radii).all() or radii[0]<0: raise ValueError('Use finite nonnegative margin radii')
    cfg=resolve_sample_guide(cfg,ims_path)
    root=Path(output_root)/datetime.now().strftime('%Y%m%d_%H%M%S_%f');root.mkdir(parents=True,exist_ok=False)
    reader=ims(str(ims_path),ResolutionLevelLock=0,squeeze_output=True)
    try: _,spacing,_,diag,pack=build_anatomical_v2_roi(reader,cfg,keep_stages=True)
    finally: reader.close()
    stages=pack['stages'];records=[]
    for r in radii:
        mask=expand_final_roi(stages['base'],stages['expansion_support'],stages['protected_optic'],spacing,
                              r,cfg['roi']['v2_final_margin_z_um'],stages['allowed'])
        d=copy.deepcopy(diag)
        sweep_config=copy.deepcopy(cfg)
        sweep_config['roi']['v2_final_margin_xy_um']=r
        d.update(central_area_by_z_um2=(mask.sum((1,2))*np.prod(spacing[1:])).tolist(),
                 final_margin_xy_um=r,added_margin_volume_um3=float((mask&~stages['base']).sum()*np.prod(spacing)),
                 roi_status='needs_review',diagnostic_sweep_not_analysis=True,
                 review_flags=diag['review_flags']+['Diagnostic fixed-base margin sweep; review the added peripheral band'],
                 input_file=str(Path(ims_path).resolve()), config=sweep_config,
                 low_spacing_zyx_um=list(spacing))
        sub=root/f'margin_xy_{r:g}um'
        paths=save_structural_roi_qc(sub,pack['volume'],mask,spacing,d,stages,title=f'{Path(ims_path).name}; XY margin {r:g} um')
        (sub/'roi_diagnostics.json').write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
        records.append({'margin_xy_um':r,'roi_volume_um3':float(mask.sum()*np.prod(spacing)),
                        'added_volume_um3':d['added_margin_volume_um3'],'output_directory':str(sub)})
    (root/'margin_sweep.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
    return records
