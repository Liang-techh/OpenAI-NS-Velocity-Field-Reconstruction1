"""Admit actual current reference/restoration tensors and three attachments."""
import copy
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_restore_background_tensor import (
    CurrentRestoreBackgroundTensor,VIEWS,DOMAINS,SEAMS,GATES,OPEN,NAME,RECEIPT,HERE,
    sha,pack,encode,endpoints,source_precision,canonical_tensor_groups,_verify_hashes,read_producer)
from lei_ren_part1_paper_compliant_current_pulse_end_background_tensor_check import physical_rows,check_common
from lei_ren_part1_paper_compliant_current_angular_background_stress_check import finite


def check_view(view):
    chart=view['chart']
    if chart not in DOMAINS or any(view[k] for k in OPEN) or not all(view[k] for k in (
            'actual_full_stress_not_local_difference','source_bounds_not_resolved_physical_point_values',
            'full_current_correlated_E_centered_histories_and_original_cutoff_jets_retained',
            'full_energy_baseline_cross_terms_and_absolute_P0_retained',
            'derivative_coordinate_is_ordinary_logR_not_selector_phase',
            'positive_amplitude_caps_enclose_original_log_source_not_define_it',
            'only_radius_and_constant_Pstar_factored_from_full_stress')):raise ValueError('Actual reference/restore tensor source or scope lost')
    restore=view['actual_upstream_current_reference_restore_source'];parent=restore['actual_inherited_axial5_packet']
    expected=dict(inner_reference='Rsh_to_Rz',axial_restore='original_axial_restoration',restore_buffer='restore_exit_to_Rm')
    if restore['chart']!=expected[chart] or not (restore['actual_reference_moment_histories_and_P0_retained'] and
            restore['positive_exponential_cap_is_enclosure_only'] and parent['actual_histories_retained'] and
            parent['original_axis_pressure_retained'] and not parent['moment_reset']):raise ValueError('Actual centered histories/log source/P0 reset')
    centered=restore['actual_centered_moment_y_derivative_axial5']
    if len(centered)!=6 or any(len(rows)!=5 or any(row.order!=5 for row in rows) for rows in centered.values()):raise ValueError('Six current centered ordinary y/axial5 histories required')
    groups=list(centered.values())+list(view['current_raw_five_history_rows'].values())+[view['current_absolute_pressure_ordinary_y_rows']]
    for group in groups:
        if len(group)!=5 or any(row.order!=5 for row in group):raise ValueError('Actual restore full ordinary y/axial5 histories required')
        for row in group:
            for value in row.coefficients:finite(value)
    for label,rows in view['current_source_three_component_velocity_rows'].items():
        if len(rows)!=5 or any(row.order<(4 if label=='radial' else 5) for row in rows):raise ValueError('Full actual restore velocity rows lost')
        for row in rows:
            for value in row.coefficients:finite(value)
    alpha=restore['actual_alpha_ordinary_y_derivatives']
    if len(alpha)!=5 or len(parent['actual_E_V110_minus_4Z_axial5_coefficients'])!=6 or len(parent['pressure_axis_axial5_coefficients'])!=6:raise ValueError('Current E/P0/cutoff derivatives incomplete')
    for value in alpha:finite(value)
    amplitude=view['actual_positive_amplitude_source_enclosures']
    for key in ('actual_positive_swirl_base_enclosure','actual_positive_squared_swirl_base_enclosure'):
        value=amplitude[key];finite(value)
        if endpoints(value)[0]<0 or endpoints(value)[1]<=0:raise ValueError('Positive exact source replaced by nonpositive cap')
    if not amplitude['positive_caps_enclose_original_exponential_not_define_it']:raise ValueError('Positive cap selected as source amplitude')
    if chart=='axial_restore':
        kernels=parent['original_restoration_kernels']
        if set(kernels)!= {'mean','mixed','square'}:raise ValueError('Full original restoration kernels missing')
        for value in kernels.values():finite(value)
    grid=view['actual_upstream_physical_spatial4_time1_packet']['physical_spatial_cartesian_mixed4']
    if len(grid)!=35 or any(set(row)!= {'ux','uy','uz','p'} for row in grid.values()):raise ValueError('Actual restore BASE spatial4/time1 packet incomplete')
    sectors=view['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']
    if len(sectors['theta'])!=5 or len(sectors['axial'])!=6:raise ValueError('Full actual restore meridional stress layout incomplete')
    for parts in sectors.values():
        for part in parts.values():
            if len(part['full_stress_mixed3_coefficient_enclosures'])!=10:raise ValueError('Actual restore full stress mixed3 lost')
            for value in part['full_stress_mixed3_coefficient_enclosures'].values():finite(value)
    rows=physical_rows(view);nonzero=0
    if len(rows)!=348 or len(canonical_tensor_groups(view))!=71:raise ValueError('Actual restore full tensor/decomposition layout incomplete')
    for row in rows:
        for key in ('signed_coefficient','physical_lambda_exponent','physical_viscosity_exponent','radial_log_prefactor'):finite(row[key])
        for value in row['actual_source_log_parts'].values():finite(value)
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None or endpoints(row['signed_coefficient'])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Zero actual source has nonzero physical bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
        if not row['positive_source_factors_not_materialized'] or not row['source_factors_combined_before_enclosure']:raise ValueError('Exact actual restore source factors lost')
    if not any(not row['exact_zero'] for parts in view['physical_three_component_remainder_mixed2']['radial'].values() for row in parts.values()):raise ValueError('Actual radial restore remainder deleted')
    return dict(rows=len(rows),nonzero=nonzero)


@source_precision
def run(field=None):
    raw=read_producer();_verify_hashes(raw)
    field=field if field is not None else CurrentRestoreBackgroundTensor(require_checked=False)
    field.assert_graph();omitted={'current_actual_restore_tensor_views','current_actual_three_restore_tensor_interfaces'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:raise ValueError('Actual restore manifest/source differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_actual_restore_tensor_views'])!=set(VIEWS):raise ValueError('Actual restore producer views/scope differ')
    source=field.proof['current_actual_restore_source_pressure_and_three_endpoint_theorem']
    normalization=field.proof['original_actual_restore_current_radius_normalization_theorem']
    if not source['passed'] or not all(source['live_original_callable_bindings'].values()) or not all(source['three_exact_source_radius_identities'].values()) or not source['actual_current_Rm_open_neighborhood_function_theorem']['passed']:raise ValueError('Actual current restore/P0/three endpoint source failed')
    if not normalization['passed'] or len(normalization['identities'])!=141 or not all(normalization['identities'].values()):raise ValueError('Full actual restore current-radius unit theorem failed')
    counts={};active=False
    for name,args in VIEWS.items():
        value=field.chart(*args)
        if encode(pack(value))!=raw['current_actual_restore_tensor_views'][name]:raise ValueError('Actual restore full source replay differs: '+name)
        counts[name]=check_view(value);alpha=value['actual_upstream_current_reference_restore_source']['actual_alpha_ordinary_y_derivatives']
        if name=='axial_restore_fresh':
            active=any(endpoints(value)!=(mp.mpf(0),mp.mpf(0)) for value in alpha[1:])
            if not any(not row['exact_zero'] for parts in value['physical_three_component_remainder_mixed2']['axial'].values() for row in parts.values()):raise ValueError('Active current restore axial viscosity deleted')
        if name in ('axial_restore_left','axial_restore_right'):
            expected=mp.mpf(1 if name.endswith('left') else 0)
            if endpoints(alpha[0])!=(expected,expected) or any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in alpha[1:]):raise ValueError('Exact original flat restore endpoint lost')
        if value['chart']=='restore_buffer' and any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in alpha):raise ValueError('Postrestore local mismatch must vanish while full histories remain')
    if not active:raise ValueError('Actual variable restoration cutoff jets replaced by constants')
    joins={};fresh={}
    for name in SEAMS:
        value=field.interface(name)
        if encode(pack(value))!=raw['current_actual_three_restore_tensor_interfaces'][name]:raise ValueError('Actual restore completed tensor trace differs')
        joins[name]=check_common(value,71)
        fresh[name]=check_common(field.interface(name,Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2'),71)
    rejected=[]
    for label,args in (('outside_reference',('inner_reference','.2','1.01')),('outside_restore',('axial_restore','.2','1.01')),
            ('past_Rm',('restore_buffer','.2','-5.99')),('before_restore_exit',('restore_buffer','.2','-7.01')),
            ('outside_Z',('restore_buffer',2,'-6.5')),('foreign_chart',('actual_patch','.2','1.25')),
            ('nonfinite_coordinate',('axial_restore','.2','inf'))):
        try:field.chart(*args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual reference/restore domain admitted')
    for label,kwargs in (('nonfinite_time',dict(log_tau='inf')),('nonpositive_viscosity',dict(viscosity=0)),('nonfinite_theta',dict(theta='inf'))):
        try:field.chart('restore_buffer','.2','-6.5',**kwargs)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual restore physical sector admitted')
    clone=copy.copy(field);clone.restore=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_current_restore_owner')
    else:raise ValueError('Foreign actual restore owner admitted')
    clone=copy.copy(field);clone.patch_tensor=copy.copy(field.patch_tensor);clone.patch_tensor.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_current_patch_parent')
    else:raise ValueError('Unchecked current patch tensor admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_actual_restore_tensor_view_count=len(VIEWS),current_actual_restore_view_counts=counts,
        current_actual_restore_physical_rows_checked=sum(v['rows'] for v in counts.values()),
        current_actual_restore_nonzero_physical_contributions_checked=sum(v['nonzero'] for v in counts.values()),
        current_actual_three_restore_tensor_interface_rows=joins,current_actual_three_fresh_restore_tensor_interface_rows=fresh,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],
        actual_current_completed_tensor_adjacent_interface_count=24,actual_current_completed_tensor_internal_interface_count=10,
        current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),current_actual_restore_source_units_and_tensor_theorem=field.proof,
        current_correlated_E_partial_histories_full_energy_baseline_same_P0_and_positive_source_logs_retained=True,
        original_active_and_flat_endpoint_cutoff_jets_checked=True,postrestore_stops_at_actual_Rm_patch_inlet=True,
        invalid_domains_and_foreign_current_owners_rejected=rejected,
        scope='25 actual regions/24 adjacent/10 internal tensor traces. Whole actual inner reference/axial restore/buffer and three completed joins through actual patch/Rh to full Gamma; upstream core/bridge/switch/reshape/Rsh/axis/angular internal/global/cone/time/energy/points/n-dependent recursion open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},all_passed=True,
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual reference/restore/buffer full tensors and three joins PASS;25 regions/24 adjacent/10 internal;core/global/time open',flush=True)
    return result


if __name__=='__main__':run()
