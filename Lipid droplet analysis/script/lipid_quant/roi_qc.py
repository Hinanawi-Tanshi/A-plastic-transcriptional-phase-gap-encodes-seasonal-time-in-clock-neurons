"""Lightweight structural ROI QC, independent of lipid-detection QC."""
from pathlib import Path
import math
import numpy as np
import matplotlib.pyplot as plt


def save_structural_roi_qc(output_dir, volume, mask, spacing, diagnostics, stages, z_indices=None, title='ROI'):
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    default=np.linspace(0,len(volume)-1,8).round().astype(int).tolist()
    if z_indices is None:
        area=np.asarray(diagnostics['central_area_by_z_um2'])
        if area.any():
            default += [int(area.argmax()),int(np.argmax(abs(np.diff(area))))] if len(area)>1 else [0]
        repair=np.asarray(diagnostics.get('topology_repair_area_by_z_um2',[]),float)
        if repair.any():
            nonzero=np.flatnonzero(repair)
            default += [int(nonzero[0]),int(repair.argmax()),int(nonzero[-1])]
        interface=np.asarray(diagnostics.get('interface_area_by_z_um2',[]),float)
        if interface.any(): default.append(int(interface.argmax()))
        zs=sorted(set(default))
    else: zs=list(dict.fromkeys(int(z) for z in z_indices))
    if not zs or min(zs)<0 or max(zs)>=len(volume): raise ValueError('Invalid structural QC z indices')
    ncols=min(4,len(zs))
    fig,axes=plt.subplots(math.ceil(len(zs)/ncols),ncols,figsize=(4*ncols,4*math.ceil(len(zs)/ncols)),squeeze=False,layout='constrained')
    for ax in axes.flat: ax.axis('off')
    for ax,z in zip(axes.flat,zs):
        p=volume[z]
        ax.imshow(p,cmap='gray',vmin=np.quantile(p,.02),vmax=np.quantile(p,.998))
        overlays=[('tissue','gold',.5),('topology_repair','magenta',.7),
                  ('definite_optic','red',.55),('interface_band','darkorange',.85),
                  ('base','deepskyblue',.8),('final','lime',.8)]
        for name,color,width in overlays:
            stack=mask if name=='final' else stages.get(name)
            if stack is None: continue
            m=stack[z]
            if m.any() and not m.all(): ax.contour(m,levels=[.5],colors=[color],linewidths=width)
        ax.set_title(f'native z={z} (ImageJ {z+1})')
    scope=diagnostics.get('roi_scope','central_brain')
    fig.suptitle(f'{title}\n{scope}; tissue gold; repair magenta; definite OL red; interface orange; '
                 f'base blue; final green. {diagnostics.get("roi_status", "")}; independent contrast')
    paths={'slices':out/'roi_optical_slices.png'}
    fig.savefig(paths['slices'],dpi=130);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    p=volume.max(0)
    axes[0].imshow(p,cmap='gray',vmin=np.quantile(p,.02),vmax=np.quantile(p,.998))
    if mask.any():axes[0].contour(mask.any(0),levels=[.5],colors=['lime'],linewidths=.8)
    for core in diagnostics.get('geometric_cores',[]):
        z,y,x=core['z_y_x']
        c={'central_core':'deepskyblue','ambiguous_retained':'darkorange',
           'definite_optic':'red'}.get(core.get('retention_class'),'tomato')
        axes[0].plot(x,y,'+',color=c);axes[0].text(x+2,y+2,str(core.get('id','')),color=c,fontsize=8)
    axes[0].set_title('XY projection for localization only');axes[0].axis('off')
    depths=np.arange(len(volume))*spacing[0]
    for key,label,c in [('tissue_area_by_z_um2','Tissue','goldenrod'),('base_area_by_z_um2','Base ROI','deepskyblue'),('central_area_by_z_um2','Final ROI','seagreen')]:
        axes[1].plot(depths,np.asarray(diagnostics[key])/1000,label=label,color=c)
    if 'threshold_candidate_area_by_z_um2' in diagnostics:
        axes[1].plot(depths,np.asarray(diagnostics['threshold_candidate_area_by_z_um2'])/1000,
                     label='Per-slice threshold candidate',color='0.45',linestyle='--',linewidth=.9)
    if np.any(diagnostics.get('topology_repair_area_by_z_um2',[])):
        axes[1].plot(depths,np.asarray(diagnostics['topology_repair_area_by_z_um2'])/1000,
                     label='Topology repair area',color='magenta',linestyle=':',linewidth=1)
    if np.any(diagnostics.get('interface_area_by_z_um2',[])):
        axes[1].plot(depths,np.asarray(diagnostics['interface_area_by_z_um2'])/1000,
                     label='Retained interface',color='darkorange',linestyle=':',linewidth=1)
    axes[1].set(xlabel='Depth (um)',ylabel='Area (1000 um2)');axes[1].legend(fontsize=9)
    fig.suptitle(title+' | '+scope+' | '+diagnostics.get('roi_status',''))
    paths['depth']=out/'roi_area_by_z.png';fig.savefig(paths['depth'],dpi=140);plt.close(fig)
    return paths
