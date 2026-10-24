"""Versioned, classical ROI segmentation with deferred anatomical identity.

Geometry is a heuristic, not an atlas. Scores are NOT probabilities. All
uncertain results remain visible; no lipid detections enter this module.
"""
from __future__ import annotations

import hashlib
import json
import copy
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.morphology import binary_closing, disk, h_maxima, remove_small_objects
from skimage.segmentation import watershed

from .anatomical_roi import ANATOMICAL_DEFAULTS, _sdf

ALGORITHM_VERSION = "anatomical_v2.2"
_IMPORTED_SOURCE_HASH = hashlib.sha256(
    Path(__file__).read_bytes() + (Path(__file__).parent/'anatomical_roi.py').read_bytes()
).hexdigest()[:16]
V2_DEFAULTS = {
    "exclude_optic_lobes": True,
    "v2_whole_min_component_volume_um3": 1000.0,
    "v2_whole_min_component_fraction": 0.02,
    "v2_whole_component_gap_um": 15.0,
    "v2_whole_terminal_distance_um": 6.0,
    "v2_envelope_z_context_um": 8.0,
    "v2_envelope_boundary_tolerance_um": 3.0,
    "v2_envelope_min_repair_area_um2": 300.0,
    "v2_envelope_max_repair_area_um2": 25000.0,
    "v2_envelope_min_boundary_contact": 0.20,
    "v2_envelope_weak_contact_ratio": 0.75,
    "v2_support_margin_um": 2.0,
    "v2_edge_sharpness_per_um": 0.020,
    "v2_support_z_context_um": 6.0,
    "v2_terminal_xy_context_um": 6.0,
    "v2_peak_prominence_um": 3.0,
    "v2_core_separation_um": 12.0,
    "v2_minimum_core_radius_um": 8.0,
    "v2_boundary_smooth_xy_um": 3.0,
    "v2_boundary_smooth_z_um": 1.0,
    "v2_final_margin_xy_um": 2.0,
    "v2_final_margin_z_um": 0.0,
    "v2_optic_protection_um": 0.0,
    "v2_interface_retention_um": 8.0,
    "v2_definite_optic_crop_penalty": 0.50,
    "v2_definite_optic_max_neck_ratio": 0.55,
    "v2_definite_optic_min_radial_gap_um": 15.0,
    "v2_definite_optic_large_core_radius_um": 20.0,
    "v2_definite_optic_large_core_radial_gap_um": 40.0,
    "v2_assignment_score_gap": 0.35,
    "v2_review_policy": "warn",  # warn or error; never silently omit a sample
}


def resolved_roi_config(config):
    result = {**ANATOMICAL_DEFAULTS, **V2_DEFAULTS, **config}
    for key in V2_DEFAULTS:
        if key == "exclude_optic_lobes":
            if not isinstance(result[key], bool):
                raise ValueError('roi.exclude_optic_lobes must be true or false, not a string/number')
        elif key == "v2_review_policy":
            if result[key] not in {"warn", "error"}:
                raise ValueError("v2_review_policy must be warn or error")
        elif isinstance(result[key], bool) or not np.isfinite(result[key]) or result[key] < 0:
            raise ValueError(f"roi.{key} must be finite and nonnegative")
    if result["v2_peak_prominence_um"] <= 0 or result["v2_edge_sharpness_per_um"] <= 0:
        raise ValueError("Peak prominence and edge sharpness must be positive")
    if result['v2_whole_min_component_fraction'] > 1:
        raise ValueError('v2_whole_min_component_fraction must be between 0 and 1')
    for key in ('v2_envelope_min_boundary_contact','v2_envelope_weak_contact_ratio',
                'v2_definite_optic_crop_penalty',
                'v2_definite_optic_max_neck_ratio'):
        if result[key] > 1:
            raise ValueError(f'roi.{key} must be between 0 and 1')
    if result['v2_envelope_max_repair_area_um2'] < result['v2_envelope_min_repair_area_um2']:
        raise ValueError('v2_envelope_max_repair_area_um2 must be >= v2_envelope_min_repair_area_um2')
    for key in ("anatomical_smooth_um", "anatomical_minimum_area_um2",
                "guide_refine_band_um", "guide_max_expansion_um", "guide_expansion_per_z_um"):
        if not np.isfinite(result[key]) or result[key] < 0:
            raise ValueError(f"roi.{key} must be finite and nonnegative")
    return result


def algorithm_fingerprint():
    # A long-running notebook/process must report the code it IMPORTED, not a
    # file edited on disk while that process was already running.
    return _IMPORTED_SOURCE_HASH


def resolve_v2_config(config):
    """Materialize effective defaults and code version, only for explicit v2."""
    roi = config.get('roi', {})
    exclude = roi.get('exclude_optic_lobes', True)
    if not isinstance(exclude, bool):
        raise ValueError('roi.exclude_optic_lobes must be a boolean')
    if not exclude:
        if roi.get('mode') != 'auto_central_brain':
            raise ValueError('Whole-brain scope requires roi.mode=auto_central_brain; a central-brain guide/mask would restrict the requested whole brain')
        # Explicitly requesting whole tissue selects the v2 tissue-support path;
        # it never uses a legacy intensity-based optic classifier.
        config = copy.deepcopy(config)
        config['roi']['central_method'] = 'anatomical_v2'
    if config.get('roi',{}).get('central_method') != 'anatomical_v2':
        return config
    result=copy.deepcopy(config)
    if result['roi'].get('mode') not in {'auto_central_brain','guided_central_brain'}:
        raise ValueError('anatomical_v2 requires auto_central_brain or guided_central_brain mode')
    result['roi']=resolved_roi_config(result['roi'])
    result['_roi_implementation']={'version':ALGORITHM_VERSION,'source_hash':algorithm_fingerprint()}
    return result


def set_optic_lobe_exclusion(config, exclude=None):
    """Copy config, optionally override scope, and resolve the actual algorithm."""
    result=copy.deepcopy(config)
    if exclude is not None:
        if not isinstance(exclude, bool):
            raise ValueError('exclude_optic_lobes must be True, False, or None')
        result['roi']['exclude_optic_lobes']=exclude
        # Keep a previously present batch override consistent with this choice.
        if 'exclude_optic_lobes' in result.get('batch',{}):
            result['batch']['exclude_optic_lobes']=exclude
    return resolve_v2_config(result)


def clean_whole_brain_components(tissue, spacing, rc):
    """Conservative 3-D debris filter, without any CB/OL classification.

    Keep the dominant connected body and sizable nearby disconnected bodies.
    Small terminal islands survive only within a finite neighborhood of those
    fixed bodies. This cannot identify attached nerves/debris by anatomy.
    """
    labels,n=ndi.label(tissue,structure=np.ones((3,3,3),bool))
    volumes=np.bincount(labels.ravel()).astype(float)*np.prod(spacing)
    volumes[0]=0
    if not n:
        return tissue.copy(), {'whole_brain_components':[], 'nonbrain_removed_volume_um3':0.0}
    main=int(volumes.argmax())
    distance=ndi.distance_transform_edt(labels!=main,sampling=spacing)
    ids=np.arange(1,n+1)
    gaps=ndi.minimum(distance,labels,index=ids)
    minimum=max(rc['v2_whole_min_component_volume_um3'],rc['v2_whole_min_component_fraction']*volumes[main])
    body_ids=[i for i,gap in zip(ids,gaps) if i==main or
              (volumes[i]>=minimum and gap<=rc['v2_whole_component_gap_um'])]
    bodies=np.isin(labels,body_ids)
    if len(body_ids)>1:
        distance=ndi.distance_transform_edt(~bodies,sampling=spacing)
    # One bounded rescue, not iterative dilation from already rescued specks.
    clean=bodies | (tissue & (distance<=rc['v2_whole_terminal_distance_um']))
    kept=np.bincount(labels[clean],minlength=n+1)*np.prod(spacing)
    records=[{'id':int(i),'volume_um3':float(volumes[i]),'distance_to_main_um':float(gap),
              'retained_volume_um3':float(kept[i]),'body':bool(i in body_ids)} for i,gap in zip(ids,gaps)]
    return clean, {'whole_brain_components':records,
                   'nonbrain_removed_volume_um3':float((tissue&~clean).sum()*np.prod(spacing)),
                   'whole_brain_cleanup_limitation':'Connectivity filter; attached non-brain tissue and distant detached brain lobes require manual review.'}


def _xy_dilate(mask, radius, spacing):
    if radius <= 0 or not mask.any():
        return mask.copy()
    return ndi.distance_transform_edt(~mask, sampling=spacing) <= radius


def repair_tissue_envelope(candidate, closed_evidence, spacing, rc, first, last):
    """Repair topology flips using local-z edge evidence, not mask extrusion.

    A projection across a bounded physical z window is used only to identify
    enclosed cavities. Neighbor outer silhouettes are never copied. A cavity
    must have plausible area and boundary contact with the current plane's
    thresholded tissue, which rejects holes formed only by unrelated objects
    moving through the projection.
    """
    repair=np.zeros_like(candidate)
    if first>last or rc['v2_envelope_z_context_um']<=0:
        return repair, {'topology_repair_records':[],
                        'topology_repair_z_context_um':rc['v2_envelope_z_context_um'],
                        'topology_repair_track_count':0,
                        'topology_repair_accepted_track_count':0}
    sz,sy,sx=spacing
    half=max(1,int(np.ceil(rc['v2_envelope_z_context_um']/sz)))
    size=2*half+1
    local_union=ndi.maximum_filter1d(closed_evidence.astype(np.uint8),size=size,
                                     axis=0,mode='constant')>0
    close=disk(max(1,round(2/np.sqrt(sy*sx))))
    minimum=rc['v2_envelope_min_repair_area_um2']
    maximum=rc['v2_envelope_max_repair_area_um2']
    records=[]
    eligible=np.zeros_like(candidate)
    strong=np.zeros_like(candidate)
    strong_contact=rc['v2_envelope_min_boundary_contact']
    weak_contact=strong_contact*rc['v2_envelope_weak_contact_ratio']
    for z in range(first,last+1):
        projected=binary_closing(local_union[z],close)
        holes=ndi.binary_fill_holes(projected)&~projected
        labels,n=ndi.label(holes)
        near=_xy_dilate(closed_evidence[z],rc['v2_envelope_boundary_tolerance_um'],(sy,sx))
        for label_id in range(1,n+1):
            component=labels==label_id
            area=float(component.sum()*sy*sx)
            boundary=ndi.binary_dilation(component)&~component
            contact=float((boundary&near).sum()/max(1,boundary.sum()))
            addition=component&~candidate[z]
            added=float(addition.sum()*sy*sx)
            plausible=minimum<=area<=maximum and added>0
            weak=plausible and contact>=weak_contact
            seed=plausible and contact>=strong_contact
            if weak: eligible[z]|=component
            if seed: strong[z]|=component
            if weak or added>=minimum:
                records.append({'z':z,'cavity_area_um2':area,'added_area_um2':added,
                                'boundary_contact_fraction':contact,
                                'weak_track_candidate':bool(weak),'strong_seed':bool(seed),
                                '_representative_yx':tuple(np.argwhere(component)[0])})
    # Hysteresis in native 3-D topology: a borderline slice is repaired only
    # when its cavity continues into a slice with strong current-plane edge
    # support. This rescues a persistent open contour without accepting an
    # isolated cavity made by unrelated objects in the local-z projection.
    tracks,track_count=ndi.label(eligible,structure=np.ones((3,3,3),bool))
    accepted_tracks=np.unique(tracks[strong])
    accepted_tracks=accepted_tracks[accepted_tracks!=0]
    accepted=np.isin(tracks,accepted_tracks)
    repair=accepted&~candidate
    for record in records:
        y,x=record.pop('_representative_yx')
        track=int(tracks[record['z'],y,x])
        record['track_id']=track or None
        record['accepted']=bool(track and track in accepted_tracks)
    return repair, {'topology_repair_records':records,
                    'topology_repair_z_context_um':rc['v2_envelope_z_context_um'],
                    'topology_repair_track_count':int(track_count),
                    'topology_repair_accepted_track_count':int(len(accepted_tracks))}


def build_tissue_support(volume, spacing, rc):
    """Track intensity components; texture never deletes an interior component.

    Each XY component has a normalized boundary-sharpness score. Confident
    sections identify the supported z interval. Intermediate sections keep
    their observed intensity masks, including weak-texture interiors. Outside
    that interval, only localized current-plane structure within a finite 3-D
    neighborhood is rescued. No iterative extrusion or internal 'terminal'.
    """
    sz, sy, sx = spacing
    smooth = ndi.gaussian_filter(volume, (0, rc["anatomical_smooth_um"] / sy,
                                          rc["anatomical_smooth_um"] / sx))
    candidate = np.zeros(volume.shape, bool)
    closed_evidence = np.zeros(volume.shape, bool)
    anchors = np.zeros(volume.shape, bool)
    local_support = np.zeros(volume.shape, bool)
    relaxed = np.zeros(volume.shape, bool)
    sharpness, thresholds, records = [], [], []
    minimum = max(1, int(np.ceil(rc["anatomical_minimum_area_um2"] / (sy*sx))))
    close = disk(max(1, round(2 / np.sqrt(sy*sx))))
    for z, plane in enumerate(smooth):
        lo, hi = np.quantile(plane, [.15, .995])
        diff = plane - ndi.gaussian_filter(plane, (1.5/sy, 1.5/sx))
        low = diff[plane <= np.median(plane)]
        noise = max(1e-6, 1.4826 * np.median(np.abs(low-np.median(low))))
        dynamic = float(hi-lo)
        if dynamic <= 3*noise:
            sharpness.append(0.0); thresholds.append(None); records.append([])
            continue
        threshold = max(float(threshold_otsu(plane)), float(lo+3*noise))
        thresholds.append(threshold)
        raw = remove_small_objects(plane > threshold, min_size=minimum)
        labels, number = ndi.label(raw)
        gradient = np.hypot(ndi.sobel(plane, 0)/(8*sy), ndi.sobel(plane, 1)/(8*sx))
        plane_records = []
        for label_id, sl in enumerate(ndi.find_objects(labels), 1):
            if sl is None:
                continue
            # Use full-plane coordinates for honest treatment of cropped faces.
            comp = labels == label_id
            edge = comp ^ ndi.binary_erosion(comp)
            edge[[0, -1], :] = False; edge[:, [0, -1]] = False
            score = float(np.quantile(gradient[edge], .7)/dynamic) if edge.any() else 0.0
            area = int(comp.sum())
            plane_records.append({"area_um2": area*sy*sx, "edge_sharpness_per_um": score})
            if score >= rc["v2_edge_sharpness_per_um"]:
                anchors[z] |= comp
        records.append(plane_records)
        # The dominant area score suppresses bright tiny artifacts in a haze.
        substantial = [r for r in plane_records if r["area_um2"] >= max(60, raw.sum()*sy*sx*.1)]
        sharpness.append(max((r["edge_sharpness_per_um"] for r in substantial), default=0.0))
        closed_evidence[z] = binary_closing(raw, close)
        candidate[z] = ndi.binary_fill_holes(closed_evidence[z])
        relaxed[z] = plane > lo + max(2*noise, .015*dynamic)
        hp = plane - ndi.gaussian_filter(plane, (8/sy, 8/sx))
        weak = (hp > max(2*noise, .025*dynamic)) & relaxed[z]
        local_support[z] = ndi.binary_fill_holes(binary_closing(remove_small_objects(weak, min_size=minimum), disk(1)))

    # Stable anchors must persist in z: isolated sharp specks cannot open a
    # whole defocused slab. A one-plane real terminal island is retained below.
    anchor_area = anchors.sum((1, 2))*sy*sx
    peak = float(anchor_area.max())
    valid = (np.asarray(sharpness) >= rc["v2_edge_sharpness_per_um"]) & (anchor_area >= max(60, .03*peak))
    width = 2*max(1, int(round(1.5/sz)))+1
    density = ndi.uniform_filter1d(valid.astype(float), width, mode="constant")
    stable = valid & (density >= .5)
    flags = []
    tissue = np.zeros_like(candidate)
    indices = np.flatnonzero(stable)
    if not len(indices):
        # Do not pass a haze off as a valid whole brain.
        flags.append("No persistent sharp tissue interval; check image/guide acquisition")
        return tissue, relaxed & False, smooth, {"support_interval": None,
            "tissue_threshold_by_z": thresholds, "edge_sharpness_by_z": sharpness,
            "candidate_area_by_z_um2": (candidate.sum((1,2))*sy*sx).tolist(),
            "threshold_candidate_area_by_z_um2": (candidate.sum((1,2))*sy*sx).tolist(),
            "topology_repair_area_by_z_um2": [0.0]*len(volume),
            "maximum_topology_repair_fraction":0.0,
            "topology_repair_records":[],"topology_repair_z_context_um":rc['v2_envelope_z_context_um'],
            "topology_repair_track_count":0,"topology_repair_accepted_track_count":0,
            "_topology_repair_mask":np.zeros_like(candidate),
            "support_review_flags": flags, "terminal_z_indices": list(range(len(volume)))}
    first, last = int(indices[0]), int(indices[-1])
    threshold_candidate_area=candidate.sum((1,2))*sy*sx
    topology_repair,repair_diag=repair_tissue_envelope(candidate,closed_evidence,spacing,rc,first,last)
    candidate |= topology_repair
    tissue[first:last+1] = candidate[first:last+1]
    # Add a fixed structural support margin. Final user-facing expansion is
    # separate and never used in core extraction or tissue construction.
    for z in range(first, last+1):
        if tissue[z].any():
            tissue[z] |= _xy_dilate(tissue[z], rc["v2_support_margin_um"], (sy, sx)) & relaxed[z]
    context = int(np.ceil(rc["v2_support_z_context_um"]/sz))
    for z in list(range(max(0, first-context), first)) + list(range(last+1, min(len(volume), last+context+1))):
        near = tissue[first:min(last+1, first+context+1)].any(0) if z < first else tissue[max(first, last-context):last+1].any(0)
        permitted = _xy_dilate(near, rc["v2_terminal_xy_context_um"], (sy, sx))
        tissue[z] = local_support[z] & permitted
    # Spatial support for final margin stays independent of the margin value.
    allowed_z = np.zeros(len(volume), bool)
    allowed_z[max(0, first-context):min(len(volume), last+context+1)] = True
    relaxed &= allowed_z[:, None, None]
    outside = np.r_[np.arange(first), np.arange(last+1, len(volume))]
    relaxed[outside] &= _xy_expand_stack(local_support[outside], 2.0, (sy, sx))
    # Components entirely unrelated to the main tissue are not made brain by
    # their brightness. Keep all components supported by persistent anchors.
    labels, count = ndi.label(tissue)
    keep = np.unique(labels[anchors & tissue]); keep = keep[keep != 0]
    tissue = np.isin(labels, keep)
    repair_area=topology_repair.sum((1,2))*sy*sx
    repair_fraction=np.divide(repair_area,np.maximum(candidate.sum((1,2))*sy*sx,1))
    return tissue, relaxed, smooth, {"support_interval": [first, last],
        "tissue_threshold_by_z": thresholds, "edge_sharpness_by_z": sharpness,
        "candidate_area_by_z_um2": (candidate.sum((1,2))*sy*sx).tolist(),
        "threshold_candidate_area_by_z_um2":threshold_candidate_area.tolist(),
        "topology_repair_area_by_z_um2": repair_area.tolist(),
        "maximum_topology_repair_fraction":float(repair_fraction.max()),
        **repair_diag,
        "_topology_repair_mask":topology_repair,
        "support_review_flags": flags, "terminal_z_indices": outside.tolist()}


def _xy_expand_stack(mask, radius, spacing):
    return np.stack([_xy_dilate(p, radius, spacing) for p in mask]) if len(mask) else mask.copy()


def _core_partition(tissue, smooth, spacing, rc):
    distance = ndi.distance_transform_edt(np.pad(tissue, ((1,1),(0,0),(0,0))), sampling=spacing)[1:-1].astype(np.float32)
    distance = ndi.gaussian_filter(distance, tuple(2/s for s in spacing))
    peaks = h_maxima(distance, rc["v2_peak_prominence_um"]).astype(bool)
    peaks &= tissue & (distance >= rc["v2_minimum_core_radius_um"])
    labels, number = ndi.label(peaks)
    coords = []
    for i, sl in enumerate(ndi.find_objects(labels), 1):
        if sl is None: continue
        local = np.where(labels[sl] == i, distance[sl], -1)
        p = np.unravel_index(np.argmax(local), local.shape)
        coords.append(tuple(a+s.start for a,s in zip(p,sl)))
    coords.sort(key=lambda c: (-float(distance[c]), c))
    unique = []
    for c in coords:
        if all(np.linalg.norm((np.asarray(c)-p)*spacing) >= rc['v2_core_separation_um'] for p in unique):
            unique.append(np.asarray(c))
    coords = [tuple(int(v) for v in c) for c in unique]
    markers = np.zeros(tissue.shape, np.int32)
    for i, c in enumerate(coords, 1): markers[c] = i
    if not coords:
        return markers, distance, [], []
    gradient = np.hypot(ndi.sobel(smooth, 1)/(8*spacing[1]), ndi.sobel(smooth, 2)/(8*spacing[2]))
    gradscale = max(float(np.quantile(gradient[tissue], .95)), 1e-6)
    elevation = -distance/max(float(distance.max()), 1e-6) + .12*np.clip(gradient/gradscale, 0, 1)
    partition = watershed(elevation, markers, mask=tissue)
    del gradient, elevation, labels, peaks
    counts = np.bincount(partition.ravel(), minlength=len(coords)+1)
    records = []
    for i,c in enumerate(coords,1):
        xy = np.asarray(c[1:])*spacing[1:]
        edge = float(min(*xy, *((np.asarray(tissue.shape[1:])-1)*spacing[1:]-xy)))
        records.append({"id":i,"z_y_x":[int(v) for v in c],"radius_um":float(distance[c]),
            "edge_clearance_um":edge,"volume_um3":float(counts[i]*np.prod(spacing))})
    # Region adjacency, saddle strength, and physical contact area.
    edges = {}
    for axis in range(3):
        a=[slice(None)]*3; b=a.copy(); a[axis]=slice(None,-1); b[axis]=slice(1,None)
        pa,pb=partition[tuple(a)],partition[tuple(b)]
        border=(pa>0)&(pb>0)&(pa!=pb)
        left,right=pa[border],pb[border]
        codes=np.minimum(left,right)*(len(coords)+1)+np.maximum(left,right)
        vals=np.minimum(distance[tuple(a)][border],distance[tuple(b)][border])
        face_area=float(np.prod([s for j,s in enumerate(spacing) if j!=axis]))
        for code in np.unique(codes):
            key=(int(code//(len(coords)+1)),int(code%(len(coords)+1)))
            chosen=vals[codes==code]
            edge=edges.setdefault(key,{"a":key[0],"b":key[1],"contact_area_um2":0.,"saddle_um":0.})
            edge["contact_area_um2"]+=float(len(chosen)*face_area)
            edge["saddle_um"]=max(edge["saddle_um"],float(chosen.max()))
    for e in edges.values():
        e["neck_ratio"]=float(e["saddle_um"]/min(records[e['a']-1]['radius_um'],records[e['b']-1]['radius_um']))
    return partition, distance, records, list(edges.values())


def assign_regions(tissue, records, edges, spacing, rc):
    """Rank central-body hypotheses by radial gap, necks, and crop risk.

    Multiple central cores are allowed. An area-weighted anatomical center,
    not field-of-view center or brightest/thickest core, defines radial order.
    The selected hypothesis and runner-up are both exposed for human review.
    """
    if not records:
        return [], [], ["No geometric cores found; use a guide or inspect tissue support"]
    weights=tissue.sum(axis=0).astype(float)
    ys=weights.sum(1); xs=weights.sum(0)
    center=np.array([np.interp(.5, np.cumsum(ys)/ys.sum(), np.arange(len(ys))),
                     np.interp(.5, np.cumsum(xs)/xs.sum(), np.arange(len(xs)))])*spacing[1:]
    for r in records:
        r['distance_to_tissue_center_um']=float(np.linalg.norm(np.asarray(r['z_y_x'][1:])*spacing[1:]-center))
        r['crop_penalty']=float(max(0,1-r['edge_clearance_um']/max(r['radius_um'],1e-6)))
    order=sorted(records,key=lambda r:(r['distance_to_tissue_center_um'],r['id']))
    scale=max(float(np.median([r['radius_um'] for r in records])),1)
    hypotheses=[]
    for k in range(1,len(order)):
        members={r['id'] for r in order[:k]}
        gap=(order[k]['distance_to_tissue_center_um']-order[k-1]['distance_to_tissue_center_um'])/scale
        cut=[e for e in edges if (e['a'] in members)!=(e['b'] in members)]
        neck=float(np.mean([e['neck_ratio'] for e in cut])) if cut else 0.
        crop=float(np.average([r['crop_penalty'] for r in order[:k]],weights=[r['volume_um3'] for r in order[:k]]))
        score=3*gap-neck-2*crop
        hypotheses.append({'central_ids':sorted(members),'score':float(score),'radial_gap_radius_units':float(gap),
                           'cut_neck_ratio':neck,'central_crop_penalty':crop})
    # Absence of an identifiable periphery is uncertain, not a forced split.
    hypotheses.append({'central_ids':[r['id'] for r in records], 'score':.15,
                       'radial_gap_radius_units':0.,'cut_neck_ratio':0.,'central_crop_penalty':0.})
    hypotheses.sort(key=lambda h:(-h['score'],h['central_ids']))
    chosen=hypotheses[0]['central_ids']
    flags=[]
    if len(hypotheses)>1 and hypotheses[0]['score']-hypotheses[1]['score'] < rc['v2_assignment_score_gap']:
        flags.append('Competing central-body hypotheses have similar scores; a guide is recommended')
    if len(chosen)==len(records): flags.append('No separate peripheral region identified; optic exclusion needs review')
    if any(r['crop_penalty']>.5 for r in records if r['id'] in chosen):
        flags.append('Selected central core is strongly affected by XY cropping')
    for r in records:
        r['assignment']='central' if r['id'] in chosen else 'peripheral'
        r['tissue_center_yx_um']=center.tolist()
    return chosen,hypotheses,flags


def classify_peripheral_regions(records, edges, chosen, rc):
    """Conservatively split peripheral candidates into clear OL and uncertain.

    Clear OL requires strong crop evidence, a far narrow-neck branch, or a
    physically large core far outside the retained central set. The latter
    prevents an entire broad OL basin from being retained merely because its
    watershed neck is wide; its medial interface is restored separately.
    Small broad-neck fragments remain ambiguous because false exclusion is
    costlier there.
    """
    central=set(chosen)
    central_records=[r for r in records if r['id'] in central]
    if not central_records:
        return [],[],['No central candidates available for conservative optic classification']
    central_radius=max(r['distance_to_tissue_center_um'] for r in central_records)
    definite=[];ambiguous=[]
    for r in records:
        if r['id'] in central:
            r['retention_class']='central_core'
            continue
        direct=[e['neck_ratio'] for e in edges
                if (e['a']==r['id'] and e['b'] in central) or
                   (e['b']==r['id'] and e['a'] in central)]
        max_neck=max(direct,default=0.0)
        radial_gap=r['distance_to_tissue_center_um']-central_radius
        cropped=r['crop_penalty']>=rc['v2_definite_optic_crop_penalty']
        narrow_far=(radial_gap>=rc['v2_definite_optic_min_radial_gap_um'] and
                    max_neck<=rc['v2_definite_optic_max_neck_ratio'])
        large_far=(radial_gap>=rc['v2_definite_optic_large_core_radial_gap_um'] and
                   r['radius_um']>=rc['v2_definite_optic_large_core_radius_um'])
        r['maximum_neck_ratio_to_central']=float(max_neck)
        r['radial_gap_from_central_um']=float(radial_gap)
        evidence=[]
        if cropped: evidence.append('xy_crop')
        if narrow_far: evidence.append('far_narrow_neck')
        if large_far: evidence.append('far_large_core')
        r['optic_evidence']=evidence
        if evidence:
            definite.append(r['id']);r['retention_class']='definite_optic'
        else:
            ambiguous.append(r['id']);r['retention_class']='ambiguous_retained'
    flags=[]
    if ambiguous:
        flags.append('Ambiguous peripheral regions retained by central-tissue-priority rule: ids='+','.join(map(str,ambiguous)))
    if any(r['assignment']=='peripheral' for r in records) and not definite:
        flags.append('No peripheral candidate met definite optic criteria; optic exclusion is conservative')
    return definite,ambiguous,flags


def interface_retention_band(retained, definite_optic, spacing, width_um):
    """Z-local XY band grown only from a touching retained/definite-OL face."""
    band=np.zeros_like(retained)
    if width_um<=0:return band
    for z in range(len(retained)):
        if retained[z].any() and definite_optic[z].any():
            labels,count=ndi.label(definite_optic[z],structure=np.ones((3,3),bool))
            near=ndi.binary_dilation(retained[z],structure=np.ones((3,3),bool))
            for label_id in range(1,count+1):
                component=labels==label_id
                seed=component&near
                if seed.any():
                    distance=ndi.distance_transform_edt(~seed,sampling=spacing[1:])
                    band[z]|=component&(distance<=width_um)
    return band


def expand_final_roi(base, support, protected_optic, spacing, xy_um, z_um=0., allowed=None):
    """Monotone expansion of a fixed base; no identity or support recomputation.

    Ellipsoidal physical dilation when both radii are positive; XY-only or
    Z-only dilation when one radius is zero. Never takes pixels from protected
    optic interiors. The exact base is always retained.
    """
    if min(xy_um,z_um)<0 or not np.isfinite([xy_um,z_um]).all():
        raise ValueError('Final margins must be finite and nonnegative')
    if not base.any() or (xy_um==0 and z_um==0): return base.copy()
    if z_um==0:
        near=_xy_expand_stack(base,xy_um,spacing[1:])
    elif xy_um==0:
        near=ndi.maximum_filter1d(base,size=2*int(np.floor(z_um/spacing[0]))+1,axis=0,mode='constant')
    else:
        near=ndi.distance_transform_edt(~base,sampling=(spacing[0]/z_um,spacing[1]/xy_um,spacing[2]/xy_um))<=1
    addition=near & support & ~protected_optic
    if allowed is not None: addition &= allowed
    return base | addition


def segment_anatomical_v2(volume, spacing, roi_config, guide_mask=None, *, return_stages=False):
    rc=resolved_roi_config(roi_config)
    exclude=rc['exclude_optic_lobes']
    if not exclude and (guide_mask is not None or rc.get('mode')=='guided_central_brain'):
        raise ValueError('Whole-brain scope cannot use a central-brain guide')
    volume=np.asarray(volume,dtype=np.float32)
    spacing=tuple(float(v) for v in spacing)
    if volume.ndim!=3 or len(spacing)!=3 or not np.isfinite(spacing).all() or min(spacing)<=0 or not np.isfinite(volume).all():
        raise ValueError('Expected finite ZYX data and positive finite ZYX spacing')
    if guide_mask is not None:
        guide_mask=np.asarray(guide_mask)
        if guide_mask.shape!=volume.shape[1:] or not np.isfinite(guide_mask).all() or not guide_mask.any():
            raise ValueError('Guide must be a finite nonempty structural XY mask')
        guide_mask=guide_mask.astype(bool)
        gz=rc['guide_z_index']
        if gz is not None and (isinstance(gz,bool) or not isinstance(gz,int) or not 0<=gz<len(volume)):
            raise ValueError('guide_z_index must be a native zero-based z within the stack')
    tissue,support,smooth,diag=build_tissue_support(volume,spacing,rc)
    topology_repair=diag.pop('_topology_repair_mask')
    flags=list(diag.pop('support_review_flags'))
    if diag.get('maximum_topology_repair_fraction',0)>0.20:
        flags.append('Topology repair supplied over 20% of at least one tissue slice; inspect magenta QC contour')
    scope=('central_brain_plus_interface' if rc['v2_interface_retention_um']>0 else 'central_brain') if exclude else 'whole_brain'
    diag.update(method=('anatomical_v2_guided' if guide_mask is not None else 'anatomical_v2') if exclude else 'anatomical_v2_whole_brain',
                roi_scope=scope, exclude_optic_lobes=exclude,
                algorithm_version=ALGORITHM_VERSION, algorithm_source_hash=algorithm_fingerprint(), native_z_preserved=True,
                guide_kind=('projection' if rc['guide_z_index'] is None else 'native_z_slice') if guide_mask is not None else None,
                guide_z_index=rc['guide_z_index'] if guide_mask is not None else None)
    records=[]; hypotheses=[]; edges=[]; allowed=None
    definite_ids=[];ambiguous_ids=[]
    partition=np.zeros(volume.shape,np.int32)
    raw=np.zeros_like(tissue); optic=np.zeros_like(tissue);interface=np.zeros_like(tissue)
    definite_optic=np.zeros_like(tissue)
    if tissue.any():
        if not exclude:
            observed=tissue
            tissue,cleanup=clean_whole_brain_components(tissue,spacing,rc)
            diag.update(cleanup, input_tissue_area_by_z_um2=(observed.sum((1,2))*np.prod(spacing[1:])).tolist())
            rejected=observed&~tissue
            allowed=~rejected
            support &= allowed
            raw=tissue.copy()
            if cleanup['nonbrain_removed_volume_um3']>0.05*max(1,observed.sum()*np.prod(spacing)):
                flags.append('Whole-brain cleanup removed over 5% of observed tissue; inspect detached lobes/debris')
        elif guide_mask is None:
            partition,distance,records,edges=_core_partition(tissue,smooth,spacing,rc)
            chosen,hypotheses,identity_flags=assign_regions(tissue,records,edges,spacing,rc)
            flags.extend(identity_flags)
            definite_ids,ambiguous_ids,class_flags=classify_peripheral_regions(records,edges,chosen,rc)
            flags.extend(class_flags)
            retained_ids=list(chosen)+ambiguous_ids
            retained=np.isin(partition,retained_ids)&(partition>0)
            definite_optic=np.isin(partition,definite_ids)&(partition>0)
            interface=interface_retention_band(retained,definite_optic,spacing,rc['v2_interface_retention_um'])
            raw=retained|interface
            optic=definite_optic&~interface
            del distance
        else:
            sd=_sdf(guide_mask,spacing[1:])
            gz=rc['guide_z_index']
            if gz is not None and not (tissue[gz]&guide_mask).any():
                raise ValueError('Guide does not overlap supported tissue at guide_z_index')
            growth=np.zeros(len(volume)) if gz is None else np.minimum(np.abs(np.arange(len(volume))-gz)*spacing[0]*rc['guide_expansion_per_z_um'],rc['guide_max_expansion_um'])
            gd=sd[None]+growth[:,None,None]; band=rc['guide_refine_band_um']
            allowed=gd>=-band
            markers=np.zeros_like(partition)
            markers[tissue&(gd>=band)]=1; markers[tissue&(gd<=-band)]=2
            if not (markers==1).any(): raise ValueError('Guide has no reliable central seed; check alignment/band')
            if not (markers==2).any(): flags.append('Guide has no peripheral seed; exterior anatomy is uncertain')
            distance=ndi.distance_transform_edt(np.pad(tissue,((1,1),(0,0),(0,0))),sampling=spacing)[1:-1].astype(np.float32)
            partition=watershed(-distance,markers,mask=tissue)
            retained=(partition==1)&allowed;definite_optic=partition==2
            interface=interface_retention_band(retained,definite_optic,spacing,rc['v2_interface_retention_um'])&allowed
            raw=retained|interface; optic=definite_optic&~interface
            del distance,gd,markers
    del smooth
    base=raw.copy()
    if raw.any() and exclude:
        sdf=np.stack([_sdf(p,spacing[1:]) if p.any() else np.full(p.shape,-max(5,3*rc['v2_boundary_smooth_xy_um']),np.float32) for p in raw])
        sd=ndi.gaussian_filter(sdf,(rc['v2_boundary_smooth_z_um']/spacing[0],rc['v2_boundary_smooth_xy_um']/spacing[1],rc['v2_boundary_smooth_xy_um']/spacing[2]))
        base=(sd>=0)&tissue&~optic
        base|=interface
        del sdf,sd
        # Rescue disconnected true terminal islands only near the original CB,
        # never using iteratively rescued masks or a stack-wide XY projection.
        context=int(np.ceil(rc['v2_support_z_context_um']/spacing[0]))
        for z in diag['terminal_z_indices']:
            near=raw[max(0,z-context):min(len(raw),z+context+1)].any(0)
            permitted=_xy_dilate(near,rc['v2_terminal_xy_context_um'],spacing[1:])
            base[z] |= tissue[z]&permitted&~optic[z]
        if allowed is not None: base &= allowed
    protected=np.zeros_like(tissue)
    for z in range(len(volume)):
        if optic[z].any(): protected[z]=_sdf(optic[z],spacing[1:])>rc['v2_optic_protection_um']
    final=expand_final_roi(base,support,protected,spacing,rc['v2_final_margin_xy_um'],rc['v2_final_margin_z_um'],allowed)
    area=final.sum((1,2))*np.prod(spacing[1:])
    if area.max()>0:
        for z in range(1,len(area)-1):
            if area[z]<.4*min(area[z-1],area[z+1]) and min(area[z-1],area[z+1])>.1*area.max():
                flags.append(f'Internal ROI area collapse at native z={z}; inspect optical slices')
        if max(area[0],area[-1])>.25*area.max(): flags.append('Substantial ROI at a stack endpoint; check z coverage/defocus')
        jumps=np.flatnonzero((np.abs(np.diff(area))>.5*np.maximum(area[:-1],area[1:])) & (np.maximum(area[:-1],area[1:])>.25*area.max()))
        if len(jumps): flags.append('Large adjacent-slice ROI change near native z='+','.join(map(str,jumps.tolist()))+'; inspect terminal support/defocus')
        interval=diag.get('support_interval')
        if interval is not None and interval[1]-interval[0]>=6:
            guard=max(3,int(np.ceil(rc['v2_support_z_context_um']/spacing[0]))+1)
            start,end=interval[0]+guard,interval[1]-guard
            pairs=np.arange(start,end,dtype=int)
            relative=np.divide(np.abs(np.diff(area)),np.maximum(area[:-1],area[1:]),
                               out=np.zeros(len(area)-1,float),where=np.maximum(area[:-1],area[1:])>0)
            internal=pairs[(relative[pairs]>.25)&(np.maximum(area[pairs],area[pairs+1])>.25*area.max())]
            if len(internal):
                flags.append('Internal adjacent-slice ROI change over 25% near native z='+','.join(map(str,internal.tolist()))+'; inspect tissue continuity')
    # Classification retention is measured before the final support margin;
    # final/tissue can exceed one because support intentionally includes a
    # narrow relaxed-intensity collar outside the thresholded tissue mask.
    fraction=float(base.sum()/max(1,tissue.sum()))
    final_to_tissue=float(final.sum()/max(1,tissue.sum()))
    if exclude and guide_mask is None and fraction>.9 and len(records)>1:
        flags.append('Over 90% of tissue retained; optic exclusion needs review')
    if not final.any(): flags.append('Empty ROI')
    diag.update(geometric_cores=records, region_edges=edges, assignment_hypotheses=hypotheses,
                definite_optic_ids=definite_ids,ambiguous_retained_ids=ambiguous_ids,
                interface_retention_um=rc['v2_interface_retention_um'],
                interface_area_by_z_um2=(interface.sum((1,2))*np.prod(spacing[1:])).tolist(),
                interface_volume_um3=float(interface.sum()*np.prod(spacing)),
                review_flags=sorted(set(flags)), roi_status='failed' if not final.any() else ('needs_review' if flags else 'auto_checked'),
                tissue_area_by_z_um2=(tissue.sum((1,2))*np.prod(spacing[1:])).tolist(),
                base_area_by_z_um2=(base.sum((1,2))*np.prod(spacing[1:])).tolist(),
                 central_area_by_z_um2=area.tolist(), retained_fraction=fraction,
                 final_to_tissue_volume_ratio=final_to_tissue,
                 base_roi_volume_um3=float(base.sum()*np.prod(spacing)), added_margin_volume_um3=float((final&~base).sum()*np.prod(spacing)),
                 final_margin_xy_um=rc['v2_final_margin_xy_um'],final_margin_z_um=rc['v2_final_margin_z_um'],
                 optic_protection_um=rc['v2_optic_protection_um'],
                 limitation='Geometry/tissue heuristic, not atlas/manual ground truth. auto_checked is not anatomical validation. central_brain_plus_interface intentionally retains an interface band; whole-brain scope is not central-brain-only.',
                ignored_v1_parameters=['anatomical_tissue_margin_um','anatomical_central_radius_ratio','anatomical_focused_support_fraction','anatomical_central_margin_um'])
    if return_stages:
        return final,diag,{'tissue':tissue,'base':base,'margin':final&~base,'partition':partition,
                           'definite_optic':definite_optic,
                           'protected_optic':protected,'interface_band':interface,
                           'topology_repair':topology_repair,
                           'expansion_support':support,'allowed':allowed}
    return final,diag


def build_anatomical_v2_roi(reader, config, *, keep_stages=False):
    from .anatomical_roi import read_structural_volume, load_guide_mask
    from skimage.transform import resize
    volume, spacing, level = read_structural_volume(reader, config)
    rc = resolved_roi_config(config['roi'])
    guide = None
    if rc['mode'] == 'guided_central_brain':
        t,c = config['input']['time_index'],config['input']['channel_index']
        native = tuple(reader.metaData[0,t,c,'shape'][-2:])
        guide = load_guide_mask(rc['guide_mask_path'],native,allow_resize=rc['guide_allow_resize'])
        guide = resize(guide.astype(np.uint8),volume.shape[1:],order=0,preserve_range=True,anti_aliasing=False)>0
    pack = None
    if keep_stages:
        mask,diag,stages = segment_anatomical_v2(volume,spacing,rc,guide,return_stages=True)
        pack = {'volume':volume,'stages':stages}
    else:
        mask,diag = segment_anatomical_v2(volume,spacing,rc,guide)
    diag['guide_mask_path'] = rc['guide_mask_path'] if guide is not None else None
    if rc.get('crop_yx') is not None:
        diag['review_flags'].append('crop_yx is applied after anatomical partition; structural curves describe the uncropped ROI')
        if diag['roi_status'] != 'failed': diag['roi_status']='needs_review'
    return mask,spacing,level,diag,pack
