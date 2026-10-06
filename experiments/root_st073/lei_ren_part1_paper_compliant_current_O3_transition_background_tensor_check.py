"""Admit the actual variable O3 transition tensor and its full power trace."""
import copy
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import (
    CurrentO3TransitionBackgroundTensor,VIEWS,SEAMS,GATES,OPEN,NAME,RECEIPT,HERE,
    sha,pack,encode,endpoints,source_precision,canonical_tensor_groups,_verify_hashes,read_producer)
from lei_ren_part1_paper_compliant_current_pulse_end_background_tensor_check import physical_rows,check_common
from lei_ren_part1_paper_compliant_current_angular_background_stress_check import finite


def check_view(view):
    if view['chart']!='O3_slope_mu' or any(view[k] for k in OPEN) or not all(view[k] for k in (
            'actual_full_stress_not_local_difference','source_bounds_not_resolved_physical_point_values',
            'all_variable_log_amplitude_ordinary_derivatives_retained','only_radius_and_constant_Pstar_factored_from_full_stress',
            'full_nonzero_histories_and_radial_remainder_retained_when_local_V_zero')):raise ValueError('Actual variable transition tensor scope lost')
    source=view['actual_upstream_original_pre_transition_source']
    groups=list(view['current_raw_five_history_rows'].values())+[view['current_absolute_pressure_ordinary_y_rows']]
    for group in groups:
        if len(group)!=5 or any(row.order!=5 for row in group):raise ValueError('Full actual five-history axial5 ordinary rows required')
        for row in group:
            for value in row.coefficients:finite(value)
    for label,rows in view['current_source_three_component_velocity_rows'].items():
        if len(rows)!=5 or any(row.order<(4 if label=='radial' else 5) for row in rows):raise ValueError('Actual variable source velocity rows incomplete')
        for row in rows:
            for value in row.coefficients:finite(value)
        if label=='axial' and any(endpoints(value)!=(mp.mpf(0),mp.mpf(0)) for row in rows for value in row.coefficients):raise ValueError('Original local transition V must be zero')
    for row in source['log_Utheta_ordinary_y_derivatives']:
        for value in row.coefficients:finite(value)
    grid=view['actual_upstream_physical_spatial4_time1_packet']['physical_spatial_cartesian_mixed4']
    if len(grid)!=35 or any(set(row)!= {'ux','uy','uz','p'} for row in grid.values()):raise ValueError('Actual BASE source spatial4/time1 trace incomplete')
    sectors=view['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']
    if len(sectors['theta'])!=5 or len(sectors['axial'])!=6:raise ValueError('Actual variable full meridional stress layout incomplete')
    for parts in sectors.values():
        for part in parts.values():
            if len(part['full_stress_mixed3_coefficient_enclosures'])!=10:raise ValueError('Actual full stress mixed3 lost')
            for value in part['full_stress_mixed3_coefficient_enclosures'].values():finite(value)
    rows=physical_rows(view);nonzero=0
    if len(rows)!=348 or len(canonical_tensor_groups(view))!=71:raise ValueError('Actual variable full tensor/decomposition layout incomplete')
    for row in rows:
        for key in ('signed_coefficient','physical_lambda_exponent','physical_viscosity_exponent','radial_log_prefactor'):finite(row[key])
        for value in row['actual_source_log_parts'].values():finite(value)
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None or endpoints(row['signed_coefficient'])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Zero source has nonzero physical bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
        if not row['positive_source_factors_not_materialized'] or not row['source_factors_combined_before_enclosure']:raise ValueError('Original exact source factors lost')
    remainder=view['physical_three_component_remainder_mixed2']
    if not any(not row['exact_zero'] for parts in remainder['radial'].values() for row in parts.values()):raise ValueError('Retained actual radial remainder collapsed with V')
    if any(not row['exact_zero'] for parts in remainder['axial'].values() for row in parts.values()):raise ValueError('Local V=0 axial viscosity must be structurally zero')
    return dict(rows=len(rows),nonzero=nonzero)


@source_precision
def run(field=None):
    raw=read_producer();_verify_hashes(raw)
    field=field if field is not None else CurrentO3TransitionBackgroundTensor(require_checked=False)
    field.assert_graph();omitted={'current_actual_O3_transition_tensor_views','current_actual_O3_transition_power_tensor_interface'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:raise ValueError('Actual transition source/manifest differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_actual_O3_transition_tensor_views'])!=set(VIEWS):raise ValueError('Actual transition producer views/scope differ')
    for key in ('original_arbitrary_variable_source_full_stress_theorem','current_actual_variable_source_velocity_and_remainder_theorem'):
        proof=field.proof[key]
        if not proof['passed'] or not all(proof['identities'].values()):raise ValueError('Original variable source tensor/operator theorem failed')
    source=field.proof['current_actual_transition_source_and_power_endpoint']
    if not source['passed'] or not all(source['live_original_callable_bindings'].values()):raise ValueError('Actual transition endpoint source not identified')
    counts={};variable=False
    for name,args in VIEWS.items():
        value=field.chart(*args)
        if encode(pack(value))!=raw['current_actual_O3_transition_tensor_views'][name]:raise ValueError('Actual transition full tensor replay differs: '+name)
        counts[name]=check_view(value)
        logjet=value['actual_upstream_original_pre_transition_source']['log_Utheta_ordinary_y_derivatives']
        if name in ('transition_left','transition_power'):
            if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for row in logjet[1:] for v in row.coefficients):raise ValueError('Original exact flat transition endpoint lost')
        if name=='fresh_transition':variable=any(endpoints(row[0])!=(mp.mpf(0),mp.mpf(0)) for row in logjet[1:])
    if not variable:raise ValueError('Interior source incorrectly replaced by constant pure-power rate')
    joins={};fresh={}
    for name in SEAMS:
        value=field.interface(name)
        if encode(pack(value))!=raw['current_actual_O3_transition_power_tensor_interface'][name]:raise ValueError('Actual completed transition/power trace differs')
        joins[name]=check_common(value,71)
        fresh[name]=check_common(field.interface(name,Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2'),71)
    rejected=[]
    for label,args in (('outside_transition',('O3_slope_mu','.2','1.01')),('outside_Z',('O3_slope_mu',2,'.5')),
            ('foreign_chart',('O3_power','.2','.5')),('nonfinite_coordinate',('O3_slope_mu','.2','inf'))):
        try:field.chart(*args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual transition domain admitted')
    for label,kwargs in (('nonfinite_time',dict(log_tau='inf')),('nonpositive_viscosity',dict(viscosity=0)),('nonfinite_theta',dict(theta='inf'))):
        try:field.chart('O3_slope_mu','.2','.5',**kwargs)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual transition physical sector admitted')
    clone=copy.copy(field);clone.pulse=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_pulse_owner')
    else:raise ValueError('Foreign actual pulse owner admitted')
    clone=copy.copy(field);clone.entrance_incoming_tensor=copy.copy(field.entrance_incoming_tensor);clone.entrance_incoming_tensor.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_current_parent')
    else:raise ValueError('Unchecked entrance/incoming tensor admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_actual_O3_transition_tensor_view_count=len(VIEWS),current_actual_O3_transition_view_counts=counts,
        current_actual_O3_transition_physical_rows_checked=sum(v['rows'] for v in counts.values()),
        current_actual_O3_transition_nonzero_physical_contributions_checked=sum(v['nonzero'] for v in counts.values()),
        current_actual_O3_transition_power_interface_rows=joins,current_actual_fresh_O3_transition_power_interface_rows=fresh,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],actual_current_completed_tensor_adjacent_interface_count=16,
        actual_current_completed_tensor_internal_interface_count=4,current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),
        current_actual_O3_transition_source_and_tensor_theorem=field.proof,
        actual_variable_log_amplitude_interior_and_exact_flat_endpoints_checked=True,
        full_nonzero_histories_and_actual_absolute_pressure_retained=True,
        local_V_zero_does_not_delete_radial_velocity_or_radial_remainder=True,
        invalid_domains_and_foreign_current_owners_rejected=rejected,
        scope='17 actual regions/16 adjacent/4 internal tensor traces. Whole actual variable O3 transition and one completed power join; upstream/core/axis/angular internal/global/cone/time/energy/points/n-dependent recursion open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},all_passed=True,
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual variable O3 transition full tensor and power join PASS;17 regions/16 adjacent/4 internal;global/time open',flush=True)
    return result


if __name__=='__main__':run()
