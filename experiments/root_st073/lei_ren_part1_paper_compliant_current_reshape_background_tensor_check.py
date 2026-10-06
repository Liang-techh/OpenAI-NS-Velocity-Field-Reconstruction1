"""Admit the actual variable long-reshape full tensor and exact Rsh trace."""
import copy
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_reshape_background_tensor import (
    CurrentReshapeBackgroundTensor,VIEWS,DOMAINS,SEAMS,GATES,OPEN,NAME,RECEIPT,HERE,
    sha,pack,encode,endpoints,source_precision,canonical_tensor_groups,_verify_hashes,read_producer)
from lei_ren_part1_paper_compliant_current_pulse_end_background_tensor_check import physical_rows,check_common
from lei_ren_part1_paper_compliant_current_angular_background_stress_check import finite


def check_view(view):
    if view['chart']!='reshape' or any(view[k] for k in OPEN) or not all(view[k] for k in (
            'actual_full_stress_not_local_difference','source_bounds_not_resolved_physical_point_values',
            'full_current_R110_moments_V110_variable_B_T_kernels_and_analytic_P0_retained',
            'full_energy_baseline_and_absolute_P0_retained','derivative_coordinate_is_ordinary_logR_not_selector_phase',
            'original_variable_Bell_products_combined_before_source_caps','only_radius_and_constant_Pstar_factored_from_full_stress')):raise ValueError('Actual variable long-reshape tensor source or scope lost')
    mixed=view['actual_upstream_current_long_reshape_source'];parent=mixed['actual_inherited_axial5_packet']
    if not all(mixed[k] for k in ('current_actual_R110_feedback_in_long_reshape',
            'actual_R110_histories_and_P0_retained','full_backward_kernels_not_truncated',
            'positive_source_cap_is_enclosure_only','large_derivative_factors_combined_before_positive_source_cap')):raise ValueError('Actual original reshape histories/log source/caps reset')
    if not parent['actual_R110_histories_retained'] or not parent['original_pressure_datum_retained']:raise ValueError('Actual R110 histories or analytic pressure missing')
    if len(parent['actual_normalized_moment_shape_axial5_coefficients'])!=6 or any(len(row)!=6 for row in parent['actual_normalized_moment_shape_axial5_coefficients'].values()):raise ValueError('Full six original reshape moment shapes required')
    for key in ('full_backward_kernel_axial5_coefficients','inherited_R110_amplitude_decay_axial5_coefficients',
            'actual_inherited_R110_moment_contribution_axial5_coefficients'):
        groups=parent[key]
        expected=4 if key=='actual_inherited_R110_moment_contribution_axial5_coefficients' else 3
        if len(groups)!=expected or any(len(row)!=6 for row in groups.values()):raise ValueError('Original full positive kernels/inherited sources missing')
        for row in groups.values():
            for value in row:finite(value)
    for group in list(view['current_raw_five_history_rows'].values())+[view['current_absolute_pressure_ordinary_y_rows']]:
        if len(group)!=5 or any(row.order!=5 for row in group):raise ValueError('Actual reshape ordinary y/axial5 raw histories incomplete')
        for row in group:
            for value in row.coefficients:finite(value)
    for label,rows in view['current_source_three_component_velocity_rows'].items():
        if len(rows)!=5 or any(row.order<(4 if label=='radial' else 5) for row in rows):raise ValueError('Full actual reshape velocity rows lost')
        for row in rows:
            for value in row.coefficients:finite(value)
    amplitude=view['actual_variable_amplitude_source_enclosures']
    for key,size in (('original_variable_log_amplitude_ordinary_y_axial5',4),
            ('original_variable_exponential_Bell_rows',5),('original_variable_squared_exponential_Bell_rows',5)):
        group=amplitude[key]
        if len(group)!=size or any(row.order!=5 for row in group):raise ValueError('Variable actual log amplitude/Bell jets incomplete')
        for row in group:
            for value in row.coefficients:finite(value)
    if not amplitude['source_derivative_factors_combined_before_positive_caps'] or not amplitude['positive_caps_enclose_original_exponential_not_define_it']:raise ValueError('Actual positive source selected as cap or bounded too early')
    grid=view['actual_upstream_physical_spatial4_time1_packet']['physical_spatial_cartesian_mixed4']
    if len(grid)!=35 or any(set(row)!= {'ux','uy','uz','p'} for row in grid.values()):raise ValueError('Actual reshape BASE spatial4/time1 packet incomplete')
    sectors=view['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']
    if len(sectors['theta'])!=5 or len(sectors['axial'])!=6:raise ValueError('Full actual reshape meridional stress layout incomplete')
    for parts in sectors.values():
        for part in parts.values():
            if len(part['full_stress_mixed3_coefficient_enclosures'])!=10:raise ValueError('Actual reshape full stress mixed3 lost')
            for value in part['full_stress_mixed3_coefficient_enclosures'].values():finite(value)
    rows=physical_rows(view);nonzero=0
    if len(rows)!=348 or len(canonical_tensor_groups(view))!=71:raise ValueError('Actual reshape full tensor/decomposition layout incomplete')
    for row in rows:
        for key in ('signed_coefficient','physical_lambda_exponent','physical_viscosity_exponent','radial_log_prefactor'):finite(row[key])
        for value in row['actual_source_log_parts'].values():finite(value)
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None or endpoints(row['signed_coefficient'])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Zero actual source has nonzero physical bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
        if not row['positive_source_factors_not_materialized'] or not row['source_factors_combined_before_enclosure']:raise ValueError('Exact actual reshape source factors lost')
    if not any(not row['exact_zero'] for parts in view['physical_three_component_remainder_mixed2']['radial'].values() for row in parts.values()):raise ValueError('Actual inherited radial reshape remainder deleted')
    return dict(rows=len(rows),nonzero=nonzero)


@source_precision
def run(field=None):
    raw=read_producer();_verify_hashes(raw)
    field=field if field is not None else CurrentReshapeBackgroundTensor(require_checked=False)
    field.assert_graph();omitted={'current_actual_reshape_tensor_views','current_actual_Rsh_tensor_interfaces'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:raise ValueError('Actual reshape manifest/source differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_actual_reshape_tensor_views'])!=set(VIEWS):raise ValueError('Actual reshape producer views/scope differ')
    source=field.proof['current_actual_variable_reshape_source_pressure_and_Rsh_theorem']
    normalization=field.proof['original_actual_variable_reshape_current_radius_normalization_theorem']
    if not source['passed'] or not all(source['live_original_callable_bindings'].values()) or not source['exact_actual_Rsh_radius_function_identity']:raise ValueError('Actual current reshape/P0/Rsh source failed')
    if not normalization['passed'] or len(normalization['identities'])!=141 or not all(normalization['identities'].values()):raise ValueError('Actual variable reshape current-radius unit theorem failed')
    counts={};active={}
    for name,args in VIEWS.items():
        value=field.chart(*args)
        if encode(pack(value))!=raw['current_actual_reshape_tensor_views'][name]:raise ValueError('Actual reshape full source replay differs: '+name)
        counts[name]=check_view(value);log_y=value['actual_variable_amplitude_source_enclosures']['original_variable_log_amplitude_ordinary_y_axial5']
        if name in ('reshape_left','reshape_right'):
            if endpoints(log_y[0][0])!=endpoints(field.ctx.mpf('.1')) or any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in log_y[0].coefficients[1:]) or any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for row in log_y[1:] for v in row.coefficients):raise ValueError('Original flat reshape endpoint jets changed')
        if 'fresh' in name:
            active[name]=endpoints(log_y[0][0])!=endpoints(field.ctx.mpf('.1')) and any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for row in log_y[1:] for v in row.coefficients)
    if not active or not all(active.values()):raise ValueError('Actual interior variable B/T/cutoff log jets replaced by constant rate')
    joins={};fresh={}
    for name in SEAMS:
        value=field.interface(name)
        if encode(pack(value))!=raw['current_actual_Rsh_tensor_interfaces'][name]:raise ValueError('Actual Rsh completed tensor trace differs')
        joins[name]=check_common(value,71)
        fresh[name]=check_common(field.interface(name,Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2'),71)
    rejected=[]
    for label,args in (('before_R110',('reshape','.2','-.01')),('past_Rsh',('reshape','.2','1.01')),
            ('outside_Z',('reshape',2,'.5')),('foreign_chart',('inner_reference','.2','.5')),
            ('nonfinite_coordinate',('reshape','.2','inf'))):
        try:field.chart(*args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual reshape domain admitted')
    for label,kwargs in (('nonfinite_time',dict(log_tau='inf')),('nonpositive_viscosity',dict(viscosity=0)),('nonfinite_theta',dict(theta='inf'))):
        try:field.chart('reshape','.2','.5',**kwargs)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual reshape physical sector admitted')
    try:field.interface('R110')
    except ValueError:rejected.append('unbuilt_R110_tensor_join')
    else:raise ValueError('Unbuilt switch-power tensor join admitted')
    clone=copy.copy(field);clone.reshape=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_current_reshape_owner')
    else:raise ValueError('Foreign actual reshape owner admitted')
    clone=copy.copy(field);clone.restore_tensor=copy.copy(field.restore_tensor);clone.restore_tensor.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_current_restore_parent')
    else:raise ValueError('Unchecked current restore tensor admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_actual_reshape_tensor_view_count=len(VIEWS),current_actual_reshape_view_counts=counts,
        current_actual_reshape_physical_rows_checked=sum(v['rows'] for v in counts.values()),
        current_actual_reshape_nonzero_physical_contributions_checked=sum(v['nonzero'] for v in counts.values()),
        current_actual_Rsh_tensor_interface_rows=joins,current_actual_fresh_Rsh_tensor_interface_rows=fresh,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],
        actual_current_completed_tensor_adjacent_interface_count=25,actual_current_completed_tensor_internal_interface_count=10,
        current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),current_actual_reshape_source_units_and_tensor_theorem=field.proof,
        whole_original_variable_B_T_log_source_Bell_kernels_and_R110_histories_checked=True,
        fresh_nonconstant_interior_log_amplitude_views=active,exact_original_flat_endpoint_log_jets_checked=True,
        invalid_domains_and_foreign_current_owners_rejected=rejected,
        scope='26 actual regions/25 adjacent/10 internal tensor traces. Whole actual long reshape and completed Rsh join through reference/restore/patch to full Gamma; upstream core/bridge/switch/R110/axis/angular internal/global/cone/time/energy/points/n-dependent recursion open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},all_passed=True,
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual variable long reshape full tensor and Rsh join PASS;26 regions/25 adjacent/10 internal;core/global/time open',flush=True)
    return result


if __name__=='__main__':run()
