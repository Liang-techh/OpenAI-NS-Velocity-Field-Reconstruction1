"""Admit actual current whole O3 power, entrance tensors and both joins."""
import copy
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_pulse_entrance_incoming_background_tensor import (
    CurrentPulseEntranceIncomingBackgroundTensor,VIEWS,SEAMS,GATES,OPEN,NAME,RECEIPT,HERE,
    sha,pack,encode,endpoints,source_precision,canonical_tensor_groups,_verify_hashes,read_producer)
from lei_ren_part1_paper_compliant_current_pulse_end_background_tensor_check import physical_rows,check_common
from lei_ren_part1_paper_compliant_current_angular_background_stress_check import finite


def check_view(view):
    incoming=view['chart']=='O3_power'
    if any(view[k] for k in OPEN) or not view['actual_full_stress_not_local_difference'] or not view['source_bounds_not_resolved_physical_point_values']:
        raise ValueError('Actual entrance/incoming full tensor scope lost')
    if not view['all_nonzero_incoming_histories_and_radial_remainder_retained']:raise ValueError('Full incoming histories lost')
    if incoming:
        if not view['actual_pre_pressure_function_used_in_tensor'] or not view['actual_full_upstream_power_not_only_unit_left_neighborhood']:raise ValueError('Actual whole pre-power source required')
        full=view['current_actual_normalized_power_rows']
        groups=[full[key] for key in ('m','n','e','P')]
        raw=view['actual_upstream_original_pre_power_source']
        if not raw['actual_five_histories_and_analytic_pressure_retained']:raise ValueError('Actual pre source history lost')
        for key in ('m','h','k','e','p'):
            groups.append(raw['actual_normalized_primitive_y_derivative_axial5'][key])
        grid=view['actual_upstream_physical_spatial4_time1_packet']['physical_spatial_cartesian_mixed4']
        if len(grid)!=35 or any(set(row)!= {'ux','uy','uz','p'} for row in grid.values()):raise ValueError('Actual BASE physical source trace incomplete')
    else:
        full=view['source_full_main_exit_history_rows']
        groups=[full[key] for key in ('Bh','ml','mi','nl','ni','e0','J','pressure_baseline_rows','pressure_memory_rows')]
        groups.append(view['current_forward_anchored_full_energy_rows'])
        if not view['selected_forward_backward_energy_same_source_not_two_added_terms']:raise ValueError('Forward/backward energy source identification lost')
        for kernel in view['original_partial_linear_kernel_bounds']:
            if not kernel['omitted_history_not_deleted']:raise ValueError('Forward history deleted')
            finite(kernel['enclosure']);finite(kernel['positive_omitted_tail_bound'])
    for group in groups:
        if len(group)!=5 or any(row.order!=5 for row in group):raise ValueError('Whole axial5 ordinary history rows required')
        for row in group:
            for value in row.coefficients:finite(value)
    sectors=view['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']
    if len(sectors['theta'])!=(4 if incoming else 6) or len(sectors['axial'])!=(6 if incoming else 9):raise ValueError('Full meridional stress sector layout incomplete')
    for parts in sectors.values():
        for part in parts.values():
            if len(part['full_stress_mixed3_coefficient_enclosures'])!=10:raise ValueError('Full stress mixed3 lost')
            for value in part['full_stress_mixed3_coefficient_enclosures'].values():finite(value)
    rows=physical_rows(view);nonzero=0
    if len(rows)!=(325 if incoming else 508) or len(canonical_tensor_groups(view))!=71:raise ValueError('Actual full tensor/decomposition layout incomplete')
    for row in rows:
        for key in ('signed_coefficient','physical_lambda_exponent','physical_viscosity_exponent','radial_log_prefactor'):finite(row[key])
        for value in row['actual_source_log_parts'].values():finite(value)
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None or endpoints(row['signed_coefficient'])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Source zero carries nonzero bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
        if not row['positive_source_factors_not_materialized'] or not row['source_factors_combined_before_enclosure']:raise ValueError('Exact source factors lost')
    radial=view['physical_three_component_remainder_mixed2']['radial']
    if not any(not row['exact_zero'] for parts in radial.values() for row in parts.values()):raise ValueError('Actual radial remainder collapsed')
    return dict(rows=len(rows),nonzero=nonzero)


@source_precision
def run(field=None):
    raw=read_producer();_verify_hashes(raw)
    field=field if field is not None else CurrentPulseEntranceIncomingBackgroundTensor(require_checked=False)
    field.assert_graph();omitted={'current_actual_entrance_incoming_tensor_views','current_actual_two_entrance_incoming_tensor_interfaces'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:raise ValueError('Current actual entrance/incoming defining source differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_actual_entrance_incoming_tensor_views'])!=set(VIEWS):raise ValueError('Producer scope/views differ')
    if not all(field.proof['current_actual_raw_power_normalization_theorem']['identities'].values()):raise ValueError('Actual pre normalization/source units differ')
    if not all(field.proof['actual_both_production_O2_O3_function_equality']['identities'].values()):raise ValueError('Actual production endpoint recurrence differs')
    if not all(field.proof['actual_pre_and_full_tensor_source_velocity_theorem']['identities'].values()):raise ValueError('Actual pre/lift source velocity differs')
    counts={}
    for name,args in VIEWS.items():
        value=field.chart(*args)
        if encode(pack(value))!=raw['current_actual_entrance_incoming_tensor_views'][name]:raise ValueError('Actual entrance/incoming source replay differs: '+name)
        counts[name]=check_view(value)
        if name=='inlet':
            for key in ('Bh','ml','nl'):
                for row in value['source_full_main_exit_history_rows'][key]:
                    if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in row.coefficients):raise ValueError('Inlet local forcing/moments must be structurally zero')
            if encode(pack(value['current_forward_anchored_full_energy_rows'][0]))!=encode(pack(value['current_incoming_energy'])):raise ValueError('Actual nonzero incoming energy anchor lost')
    joins={};fresh={}
    for name in SEAMS:
        value=field.interface(name)
        if encode(pack(value))!=raw['current_actual_two_entrance_incoming_tensor_interfaces'][name]:raise ValueError('Actual full entrance/incoming tensor trace differs')
        joins[name]=check_common(value,71)
        fresh[name]=check_common(field.interface(name,Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2'),71)
    rejected=[]
    for label,args in (('outside_entrance',('pulse_entrance','.2','.03')),('outside_power',('O3_power','.2',2)),
            ('outside_early_y',('early_y','.2',2)),('outside_Z',('pulse_entrance',2,'.01')),
            ('foreign_chart',('pulse_main','.2',1)),('nonfinite_phase',('O3_power','.2','inf'))):
        try:field.chart(*args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual entrance/incoming domain admitted')
    for label,kwargs in (('nonfinite_time',dict(log_tau='inf')),('nonpositive_viscosity',dict(viscosity=0)),('nonfinite_theta',dict(theta='inf'))):
        try:field.chart('O3_power','.2','.5',**kwargs)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual power physical sector admitted')
    clone=copy.copy(field);clone.pulse=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_pulse_owner')
    else:raise ValueError('Foreign actual pulse admitted')
    clone=copy.copy(field);clone.main_tensor=copy.copy(field.main_tensor);clone.main_tensor.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_current_main_tensor')
    else:raise ValueError('Unchecked current tensor graph admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_actual_entrance_incoming_tensor_view_count=len(VIEWS),current_actual_entrance_incoming_view_counts=counts,
        current_actual_entrance_incoming_physical_rows_checked=sum(v['rows'] for v in counts.values()),
        current_actual_entrance_incoming_nonzero_physical_contributions_checked=sum(v['nonzero'] for v in counts.values()),
        current_actual_two_entrance_incoming_tensor_interface_rows=joins,current_actual_two_fresh_entrance_incoming_tensor_interface_rows=fresh,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],actual_current_completed_tensor_adjacent_interface_count=15,
        actual_current_completed_tensor_internal_interface_count=4,current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),
        current_actual_pre_power_and_entrance_source_theorem=field.proof,
        full_nonzero_incoming_energy_moments_and_actual_absolute_pressure_retained=True,
        actual_pre_tensor_not_a_synthetic_negative_pulse_extension=True,actual_BASE_power_source_trace_retained=True,
        invalid_domains_and_foreign_current_owners_rejected=rejected,
        scope='16 actual regions/15 adjacent/4 internal tensor traces.Whole actual O3 power and entrance tensors,two actual joins; remaining upstream pre/core/axis/angular internal/global/cone/time/energy/points/n-dependent recursion open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},all_passed=True,
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual O3 power/entrance tensors and both joins PASS;16 regions/15 adjacent/4 internal;global/time open',flush=True)
    return result


if __name__=='__main__':run()
