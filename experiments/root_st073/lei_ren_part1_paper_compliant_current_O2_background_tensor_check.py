"""Admit actual Rh-reference/O2 full tensors and four completed attachments."""
import copy
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_O2_background_tensor import (
    CurrentO2BackgroundTensor,VIEWS,SEAMS,DOMAINS,GATES,OPEN,NAME,RECEIPT,HERE,
    sha,pack,encode,endpoints,source_precision,canonical_tensor_groups,_verify_hashes,read_producer)
from lei_ren_part1_paper_compliant_current_pulse_end_background_tensor_check import physical_rows,check_common
from lei_ren_part1_paper_compliant_current_angular_background_stress_check import finite


def check_view(view):
    chart=view['chart']
    if chart not in DOMAINS or any(view[k] for k in OPEN) or not all(view[k] for k in (
            'actual_full_stress_not_local_difference','source_bounds_not_resolved_physical_point_values',
            'all_variable_log_amplitude_ordinary_derivatives_retained','only_radius_and_constant_Pstar_factored_from_full_stress',
            'full_nonzero_histories_and_radial_remainder_retained_when_local_V_zero')):raise ValueError('Actual full O2 source tensor scope lost')
    source=view['actual_upstream_original_pre_O2_source']
    expected=dict(Rh_reference='Rh_to_Rref',O2_slope='O2_slope',O2_axial='O2_axial_turnoff',O2_buffer='O2_11_unit_buffer')
    if source['chart']!=expected[chart] or not source['actual_five_histories_and_analytic_pressure_retained']:raise ValueError('Actual O2 source provider/history lost')
    groups=list(view['current_raw_five_history_rows'].values())+[view['current_absolute_pressure_ordinary_y_rows']]
    for group in groups:
        if len(group)!=5 or any(row.order!=5 for row in group):raise ValueError('Actual O2 five-history axial5 ordinary rows required')
        for row in group:
            for value in row.coefficients:finite(value)
    for label,rows in view['current_source_three_component_velocity_rows'].items():
        if len(rows)!=5 or any(row.order<(4 if label=='radial' else 5) for row in rows):raise ValueError('Actual O2 full variable velocity rows lost')
        for row in rows:
            for value in row.coefficients:finite(value)
        if chart=='O2_buffer' and label=='axial' and any(endpoints(value)!=(mp.mpf(0),mp.mpf(0)) for row in rows for value in row.coefficients):raise ValueError('Original buffer local V must be zero')
    grid=view['actual_upstream_physical_spatial4_time1_packet']['physical_spatial_cartesian_mixed4']
    if len(grid)!=35 or any(set(row)!= {'ux','uy','uz','p'} for row in grid.values()):raise ValueError('Actual O2 BASE spatial4/time1 packet incomplete')
    if chart in ('O2_axial','O2_buffer'):
        cutoff=source['original_cutoff_ordinary_y_derivatives']
        if len(cutoff)!=5 or not source['retained_far_tail_not_reset']:raise ValueError('Original turnoff source derivative/tail scope lost')
        for row in cutoff:finite(row)
        kernels=source['original_turnoff_kernel_enclosures']
        for key in ('B_mass','B_squared_mass','retained_far_tail'):finite(kernels[key])
        if not kernels['exponential_cell_weights_integrated_exactly'] or not kernels['integrals_normalized_at_current_radius']:raise ValueError('Original turnoff kernel construction lost')
    sectors=view['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']
    if len(sectors['theta'])!=5 or len(sectors['axial'])!=6:raise ValueError('Full O2 meridional stress layout incomplete')
    for parts in sectors.values():
        for part in parts.values():
            if len(part['full_stress_mixed3_coefficient_enclosures'])!=10:raise ValueError('Actual O2 full stress mixed3 lost')
            for value in part['full_stress_mixed3_coefficient_enclosures'].values():finite(value)
    rows=physical_rows(view);nonzero=0
    if len(rows)!=348 or len(canonical_tensor_groups(view))!=71:raise ValueError('Actual O2 full tensor/decomposition layout incomplete')
    for row in rows:
        for key in ('signed_coefficient','physical_lambda_exponent','physical_viscosity_exponent','radial_log_prefactor'):finite(row[key])
        for value in row['actual_source_log_parts'].values():finite(value)
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None or endpoints(row['signed_coefficient'])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Zero source has nonzero physical bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
        if not row['positive_source_factors_not_materialized'] or not row['source_factors_combined_before_enclosure']:raise ValueError('Exact actual O2 source factors lost')
    remainder=view['physical_three_component_remainder_mixed2']
    if not any(not row['exact_zero'] for parts in remainder['radial'].values() for row in parts.values()):raise ValueError('Actual radial O2 remainder deleted')
    if chart=='O2_buffer' and any(not row['exact_zero'] for parts in remainder['axial'].values() for row in parts.values()):raise ValueError('Buffer local V=0 axial viscosity must be structurally zero')
    return dict(rows=len(rows),nonzero=nonzero)


@source_precision
def run(field=None):
    raw=read_producer();_verify_hashes(raw)
    field=field if field is not None else CurrentO2BackgroundTensor(require_checked=False)
    field.assert_graph();omitted={'current_actual_O2_tensor_views','current_actual_four_O2_tensor_interfaces'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:raise ValueError('Actual O2 tensor manifest/source differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_actual_O2_tensor_views'])!=set(VIEWS):raise ValueError('Actual O2 producer views/scope differ')
    source=field.proof['current_actual_O2_source_and_four_endpoint_theorem']
    if not source['passed'] or not all(source['live_original_callable_bindings'].values()) or not all(source['actual_route_AST_replay'].values()) or not all(source['exact_radius_identities'].values()):raise ValueError('Actual O2 endpoint source program differs')
    if not source['ordinary_turnoff_cutoff_theorem']['passed'] or not all(source['ordinary_turnoff_cutoff_theorem']['identities'].values()):raise ValueError('Original ordinary logR cutoff theorem failed')
    counts={};slope_variable=False;axial_variable=False
    for name,args in VIEWS.items():
        value=field.chart(*args)
        if encode(pack(value))!=raw['current_actual_O2_tensor_views'][name]:raise ValueError('Actual O2 full tensor source replay differs: '+name)
        counts[name]=check_view(value);pre=value['actual_upstream_original_pre_O2_source']
        if name=='O2_slope_fresh':slope_variable=any(endpoints(row[0])!=(mp.mpf(0),mp.mpf(0)) for row in pre['log_Utheta_ordinary_y_derivatives'][1:])
        if name=='O2_axial_fresh':
            axial_variable=any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in pre['original_cutoff_ordinary_y_derivatives'][1:])
            if not any(not row['exact_zero'] for parts in value['physical_three_component_remainder_mixed2']['axial'].values() for row in parts.values()):raise ValueError('Actual active axial viscosity lost')
        if name in ('O2_axial_left','O2_axial_right'):
            cutoff=pre['original_cutoff_ordinary_y_derivatives'];expected=mp.mpf(1 if name.endswith('left') else 0)
            if endpoints(cutoff[0])!=(expected,expected) or any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in cutoff[1:]):raise ValueError('Original exact flat turnoff endpoint lost')
    if not slope_variable or not axial_variable:raise ValueError('Actual variable O2 transition source replaced by constant jets')
    joins={};fresh={}
    for name in SEAMS:
        value=field.interface(name)
        if encode(pack(value))!=raw['current_actual_four_O2_tensor_interfaces'][name]:raise ValueError('Actual O2 completed tensor trace differs')
        joins[name]=check_common(value,71)
        fresh[name]=check_common(field.interface(name,Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2'),71)
    rejected=[]
    for label,args in (('outside_reference',('Rh_reference','.2','-5.01')),('outside_slope',('O2_slope','.2','1.01')),
            ('outside_axial',('O2_axial','.2','1.01')),('outside_buffer',('O2_buffer','.2','11.01')),
            ('outside_Z',('O2_buffer',2,1)),('foreign_chart',('O3_slope_mu','.2','.5')),('nonfinite_coordinate',('O2_slope','.2','inf'))):
        try:field.chart(*args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual O2 domain admitted')
    for label,kwargs in (('nonfinite_time',dict(log_tau='inf')),('nonpositive_viscosity',dict(viscosity=0)),('nonfinite_theta',dict(theta='inf'))):
        try:field.chart('O2_buffer','.2','5.5',**kwargs)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual O2 physical sector admitted')
    clone=copy.copy(field);clone.pulse=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_pulse_owner')
    else:raise ValueError('Foreign actual pulse owner admitted')
    clone=copy.copy(field);clone.o3_tensor=copy.copy(field.o3_tensor);clone.o3_tensor.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_current_O3_parent')
    else:raise ValueError('Unchecked current O3 tensor admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_actual_O2_tensor_view_count=len(VIEWS),current_actual_O2_view_counts=counts,
        current_actual_O2_physical_rows_checked=sum(v['rows'] for v in counts.values()),
        current_actual_O2_nonzero_physical_contributions_checked=sum(v['nonzero'] for v in counts.values()),
        current_actual_four_O2_tensor_interface_rows=joins,current_actual_four_fresh_O2_tensor_interface_rows=fresh,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],actual_current_completed_tensor_adjacent_interface_count=20,
        actual_current_completed_tensor_internal_interface_count=4,current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),
        current_actual_O2_source_and_tensor_theorem=field.proof,
        actual_variable_slope_and_ordinary_logR_axial_cutoff_jets_checked=True,
        full_axial_transport_pressure_histories_and_positive_kernel_tails_retained=True,
        local_V_zero_does_not_delete_accumulated_histories_or_radial_remainder=True,
        invalid_domains_and_foreign_current_owners_rejected=rejected,
        scope='21 actual regions/20 adjacent/4 internal tensor traces. Whole actual Rh reference and O2 slope/axial/buffer, four completed joins; repaired-core/Rh tensor/axis/angular internal/global/cone/time/energy/points/n-dependent recursion open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},all_passed=True,
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual Rh reference/O2 full tensors and four joins PASS;21 regions/20 adjacent/4 internal;core/global/time open',flush=True)
    return result


if __name__=='__main__':run()
