"""Admit current whole gap/gap-end tensors and two actual tensor joins."""
import copy
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_current_pulse_gap_background_tensor import (
    CurrentPulseGapBackgroundTensor,VIEWS,SEAMS,GATES,OPEN,NAME,RECEIPT,HERE,
    sha,pack,encode,endpoints,source_precision,canonical_tensor_groups,_verify_hashes)
from lei_ren_part1_paper_compliant_current_pulse_end_background_tensor_check import physical_rows,check_common
from lei_ren_part1_paper_compliant_current_angular_background_stress_check import finite


def check_view(view):
    if any(view[k] for k in OPEN) or not all(view[k] for k in ('actual_full_stress_not_local_difference',
            'current_selected_complete_history_graph_retained','pressure_remaining_reduction_before_enclosure',
            'effective_radial_remainder_scale_is_D1')):raise ValueError('Current full gap source scope lost')
    full=view['source_five_full_history_rows']
    for key in ('m','n','e0','J','pressure_baseline_rows','pressure_memory_rows'):
        if len(full[key])!=5 or any(row.order!=5 for row in full[key]):raise ValueError('Full axial5 ordinary history required')
        for row in full[key]:
            for value in row.coefficients:finite(value)
    if view['complete_current_future'].order!=5 or view['full_signed_absolute_Rv_pressure'].order!=5:raise ValueError('Current full C5 future/pressure source lost')
    for row in view['current_source_three_component_velocity_rows']['axial']:
        if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in row.coefficients):raise ValueError('Inactive axial input must be structurally zero')
    source=view['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']
    if len(source['theta'])!=5 or len(source['axial'])!=4:raise ValueError('Full inactive-gap history sector splitter incomplete')
    rows=physical_rows(view);nonzero=0
    if len(rows)!=290 or len(canonical_tensor_groups(view))!=71:raise ValueError('Full actual gap tensor/decomposition layout incomplete')
    for row in rows:
        for key in ('signed_coefficient','physical_lambda_exponent','physical_viscosity_exponent','radial_log_prefactor'):finite(row[key])
        for value in row['actual_source_log_parts'].values():finite(value)
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None or endpoints(row['signed_coefficient'])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Source zero carries nonzero bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
        if not row['positive_source_factors_not_materialized'] or not row['source_factors_combined_before_enclosure']:raise ValueError('Current exact physical source factors lost')
    radial=view['physical_three_component_remainder_mixed2']['radial']
    if not any(not r['exact_zero'] for sectors in radial.values() for r in sectors.values()):raise ValueError('Actual inherited radial remainder collapsed')
    return dict(rows=len(rows),nonzero=nonzero)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentPulseGapBackgroundTensor(require_checked=False)
    field.assert_graph();omitted={'current_actual_gap_gapend_tensor_views','current_actual_two_gap_tensor_interfaces'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:raise ValueError('Current gap tensor defining source/manifest differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_actual_gap_gapend_tensor_views'])!=set(VIEWS):raise ValueError('Producer gap acceptance/views exceed scope')
    counts={}
    for name,args in VIEWS.items():
        value=field.chart(*args)
        if encode(pack(value))!=raw['current_actual_gap_gapend_tensor_views'][name]:raise ValueError('Actual gap source replay differs: '+name)
        counts[name]=check_view(value)
    joins={};fresh={}
    for name in SEAMS:
        value=field.interface(name)
        if encode(pack(value))!=raw['current_actual_two_gap_tensor_interfaces'][name]:raise ValueError('Actual gap tensor function trace replay differs')
        joins[name]=check_common(value,71)
        fresh[name]=check_common(field.interface(name,Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2'),71)
    for args in (('pulse_gap_end',0,'D=1,xi=12,s=-1/mu'),('pulse_gap_end',1,'D=4mu,s=-4')):
        value=field.chart(args[0],'.537',args[1])
        if value['source_reduced_endpoint']!=args[2]:raise ValueError('Rounded coverage replaced exact reciprocal source point')
    rejected=[]
    for label,args in (('outside_gap',('pulse_gap','.2',13)),('outside_end_phase',('pulse_gap_end','.2',2)),
            ('outside_Z',('pulse_gap',2,11)),('foreign_chart',('pulse_end','.2',0)),('nonfinite_coordinate',('pulse_gap','.2','inf'))):
        try:field.chart(*args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual gap domain admitted')
    for label,kwargs in (('nonfinite_time',dict(log_tau='inf')),('nonpositive_viscosity',dict(viscosity=0)),('nonfinite_theta',dict(theta='inf'))):
        try:field.chart('pulse_gap','.2',11,**kwargs)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual gap physical sector admitted')
    clone=copy.copy(field);clone.pulse=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_or_old_native_pulse')
    else:raise ValueError('Foreign actual gap pulse owner admitted')
    clone=copy.copy(field);clone.end_tensor=copy.copy(field.end_tensor);clone.end_tensor.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_current_end_tensor')
    else:raise ValueError('Unchecked current end tensor admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_actual_gap_gapend_tensor_view_count=len(VIEWS),current_actual_gap_gapend_view_counts=counts,
        current_actual_gap_gapend_physical_rows_checked=sum(v['rows'] for v in counts.values()),
        current_actual_gap_gapend_nonzero_physical_contributions_checked=sum(v['nonzero'] for v in counts.values()),
        current_actual_two_gap_tensor_interface_rows=joins,current_actual_two_fresh_gap_tensor_interface_rows=fresh,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],actual_current_completed_tensor_adjacent_interface_count=11,
        actual_current_completed_tensor_internal_interface_count=4,current22_velocity_pressure_interface_inventory=dict(adjacent=14,internal=8),
        current_source_function_tensor_trace_theorem=field.proof,current_full_C5_future_and_selected_histories_retained=True,
        full_nonzero_radial_remainder_and_absolute_pressure_retained=True,exact_reciprocal_source_reduction_before_enclosure=True,
        invalid_domains_and_foreign_current_owners_rejected=rejected,
        scope='Actual current gap/gap-end completed physical tensors and two tensor function joins.12 regions/11 adjacent/4 internal tensor traces; remaining main/exit/entrance/core/axis/angular internal/global cone/lift/temporal remainder/energy/points/n-dependent recursion open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},all_passed=True,
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual gap/gap-end tensor and two joins PASS;12 regions/11 adjacent/4 internal;global/time open',flush=True)
    return result


if __name__=='__main__':run()
