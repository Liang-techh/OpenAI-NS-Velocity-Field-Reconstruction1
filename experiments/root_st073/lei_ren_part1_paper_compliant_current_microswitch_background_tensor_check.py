"""Admit whole original micro-switch tensors and completed phase1/R2 traces."""
import copy
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_microswitch_background_tensor import (
    CurrentMicroswitchBackgroundTensor,VIEWS,DOMAINS,SEAMS,GATES,OPEN,NAME,RECEIPT,HERE,
    sha,pack,encode,endpoints,source_precision,canonical_tensor_groups,_verify_hashes,read_producer)
from lei_ren_part1_paper_compliant_current_microswitch_stress_operator import SOURCE_KEY
from lei_ren_part1_paper_compliant_current_pulse_end_background_tensor_check import physical_rows,check_common
from lei_ren_part1_paper_compliant_current_angular_background_stress_check import finite

THETA={'local_transport','retained_angular_moment','retained_mixed_moment','meridional_transport','variable_radial_shear'}
AXIAL={'local_axial_transport','nonlinear_meridional_transport','retained_linear_axial_moment','retained_full_energy','actual_absolute_pressure','axial_radial_shear'}


def check_factored_record(row,minimum=4):
    if row['axial_order']<minimum or not row['fixed_basepoint_log_factors']:raise ValueError('Actual source jet order or fixed log bases lost')
    exponents=[]
    for term in row['terms']:
        key=term['source_exponents'];coefficients=term['axial_Taylor_coefficients']
        if len(key)!=4 or len(coefficients)!=row['axial_order']+1 or key in exponents:raise ValueError('Complete collected factored source terms required')
        exponents.append(key)
        for value in key+coefficients:finite(value)
    return len(exponents)


def check_view(view,c):
    if view['chart'] not in DOMAINS or any(view[k] for k in OPEN) or not all(view[k] for k in (
            'inverse_width_source_shift_before_any_resolution','every_original_stress_and_remainder_sector_expanded_without_pruning',
            'original_radial_prefactors_and_absolute_P0_included_once',
            'full_current_signed_axial_drive_comparison_own_moments_and_actual_histories_retained',
            'derivative_coordinate_is_ordinary_logR_not_selector_phase','source_factors_retained_through_full_tensor_and_physical_operators',
            'source_locals_exposed_only_without_operator_changes','actual_full_stress_not_local_difference',
            'source_bounds_not_resolved_physical_point_values')):raise ValueError('Actual micro-switch full source scope lost')
    original=view['actual_upstream_current_microswitch_source']
    if not original['inverse_hb_not_materialized'] or not original['phase_derivatives_do_not_differentiate_width_caps'] or not original['comparison_own_moments_retained'] or not original['actual_moments_not_replaced_by_comparison']:raise ValueError('Original width/comparison/actual history source lost')
    for key in ('exact_positive_width_log','width_enclosure_is_not_source','R_enclosure_only','exact_positive_swirl_source_log'):finite(original[key])
    if not original['formal_radius_tree'] or original['exact_source_R']!='100*exp(hb*s)':raise ValueError('Original formal current source radius lost')
    for value in original['actual_parent_axial5_packet']['pressure_axis_axial5_coefficients']:finite(value)
    logs=view['fixed_current_factored_source_log_bases']
    if len(logs)!=4 or endpoints(logs[0])!=endpoints(original['exact_positive_width_log']):raise ValueError('Original four exact factor log sources lost')
    for value in logs:finite(value)
    raw=view['current_unresolved_raw_source_rows'];phase=view['original_unresolved_source_phase_rows'];terms=0
    if set(raw['histories'])!={'m','h','k','e','p'} or set(raw['velocity'])!={'radial','theta','axial'} or len(raw['absolute_pressure'])!=5:raise ValueError('Full actual raw source layout incomplete')
    for rows in list(raw['histories'].values())+list(raw['velocity'].values())+[raw['absolute_pressure']]:
        if len(rows)!=5:raise ValueError('Actual ordinary y0..4 source rows lost')
        for row in rows:terms+=check_factored_record(row)
    if set(phase['physical'])!={'Utheta_over_current_Utheta','Uz','Ur_over_current_sqrt_R_over_2','P_over_Pstar2'} or len(phase['primitives'])!=5 or len(phase['Q'])!=5:raise ValueError('Original complete unresolved phase source rows lost')
    for rows in list(phase['physical'].values())+list(phase['primitives'].values())+[phase['Q']]:
        for row in rows:check_factored_record(row)
    sectors=view['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']
    if {part['original_stress_sector'] for part in sectors['theta'].values()}!=THETA or {part['original_stress_sector'] for part in sectors['axial'].values()}!=AXIAL:raise ValueError('Original five/six full stress sectors omitted')
    for parts in sectors.values():
        for part in parts.values():
            if len(part['full_stress_mixed3_coefficient_enclosures'])!=10 or len(part['source_exponents'])!=4:raise ValueError('Complete expanded stress mixed3 required')
            for value in part['full_stress_mixed3_coefficient_enclosures'].values():finite(value)
            for value in part['exact_source_log_parts'].values():finite(value)
    grid=view['actual_upstream_physical_spatial4_time1_packet']['physical_spatial_cartesian_mixed4']
    if len(grid)!=35 or any(set(row)!={'ux','uy','uz','p'} for row in grid.values()):raise ValueError('Actual current BASE spatial4/time1 packet lost')
    rows=physical_rows(view);nonzero=0
    if len(canonical_tensor_groups(view))!=71 or len(rows)<348:raise ValueError('Actual completed tensor inventory lost')
    for row in rows:
        for key in ('signed_coefficient','physical_lambda_exponent','physical_viscosity_exponent','radial_log_prefactor'):finite(row[key])
        for value in row['actual_source_log_parts'].values():finite(value)
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None or endpoints(row['signed_coefficient'])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Zero actual source has a nonzero bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
        if not row['positive_source_factors_not_materialized'] or not row['source_factors_combined_before_enclosure']:raise ValueError('Actual source factor resolved before final bound')
    if not any(not row['exact_zero'] for parts in view['physical_three_component_remainder_mixed2']['radial'].values() for row in parts.values()):raise ValueError('Full actual radial remainder removed')
    return dict(rows=len(rows),nonzero=nonzero,unresolved_raw_factored_terms=terms,
        theta_source_sectors=len(sectors['theta']),axial_source_sectors=len(sectors['axial']))


@source_precision
def run(field=None):
    raw=read_producer();_verify_hashes(raw)
    field=field if field is not None else CurrentMicroswitchBackgroundTensor(require_checked=False)
    field.assert_graph();omitted={'current_actual_microswitch_tensor_views','current_actual_phase1_R2_tensor_interfaces'}
    if encode(pack(field.manifest()))!={key:value for key,value in raw.items() if key not in omitted}:raise ValueError('Actual micro-switch manifest/source differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_actual_microswitch_tensor_views'])!=set(VIEWS):raise ValueError('Actual micro-switch producer views/scope differ')
    proof=field.proof['current_actual_microswitch_source_pressure_units_and_two_endpoint_theorem']
    unit=proof['original_phase_to_ordinary_y_and_full_raw_history_normalization'];boundary=proof['original_phase1_and_R2_complete_source_function_theorem']
    factors=proof['exact_formal_radius_four_factor_and_current_tensor_packet_units']
    if not factors['passed'] or not all(factors['complete_current_physical_unit_identities'].values()) or not all(factors['current_tensor_packet_unit_bindings'].values()):raise ValueError('Actual formal radius/R/Pstar source factors changed or counted twice')
    if not proof['passed'] or not all(proof['live_original_callable_bindings'].values()) or len(unit['identities'])!=141 or not all(unit['identities'].values()) or boundary['actual_R2_source_rows_verified']!=135 or not all(boundary['actual_R2_original_physical_and_primitive_mixed4_identities'].values()):raise ValueError('Original full source unit/endpoint theorem failed')
    counts={};unchanged=[];variable={}
    for name,args in VIEWS.items():
        value=field.chart(*args)
        if encode(pack(value))!=raw['current_actual_microswitch_tensor_views'][name]:raise ValueError('Actual full micro-switch source replay differs: '+name)
        counts[name]=check_view(value,field.ctx)
        if name in ('first_whole','second_whole','first_fresh','second_fresh'):
            original=field.switch.evaluate(args[1],args[2],'first' if args[0]=='switch_first' else 'second')
            if encode(pack(value['actual_upstream_current_microswitch_source']))!=encode(pack(original)):raise ValueError('Original evaluate output changed beyond exposed locals')
            unchanged.append(name)
        if name in ('first_whole','second_whole'):
            lo,hi=endpoints(value['coverage_coordinate']);variable[name]=lo<hi
        print('Check actual microscopic full tensor: '+name,flush=True)
    if not all(variable.values()):raise ValueError('Whole original branches replaced by endpoint samples')
    joins={};fresh={}
    for name in SEAMS:
        value=field.interface(name)
        if encode(pack(value))!=raw['current_actual_phase1_R2_tensor_interfaces'][name]:raise ValueError('Actual phase1/R2 completed tensor trace differs')
        joins[name]=check_common(value,71);fresh[name]=check_common(field.interface(name,Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2'),71)
    rejected=[]
    for label,args in (('before_first',('switch_first','.2','-.01')),('past_first',('switch_first','.2','1.01')),
            ('before_second',('switch_second','.2','.99')),('past_second',('switch_second','.2','2.01')),
            ('outside_Z',('switch_first',2,'.5')),('foreign_chart',('switch_power','.2','.5')),
            ('nonfinite_coordinate',('switch_first','.2','inf'))):
        try:field.chart(*args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual micro-switch domain admitted')
    for label,kwargs in (('nonfinite_time',dict(log_tau='inf')),('nonpositive_viscosity',dict(viscosity=0)),('nonfinite_theta',dict(theta='inf'))):
        try:field.chart('switch_first','.2','.5',**kwargs)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual micro-switch physical sector admitted')
    try:field.interface('R100')
    except ValueError:rejected.append('unbuilt_R100_tensor_join')
    else:raise ValueError('Unbuilt bridge tensor join admitted')
    clone=copy.copy(field);clone.switch=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_current_switch_owner')
    else:raise ValueError('Foreign actual switch owner admitted')
    clone=copy.copy(field);clone.power_tensor=copy.copy(field.power_tensor);clone.power_tensor.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_current_power_parent')
    else:raise ValueError('Unchecked current power tensor admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_actual_microswitch_tensor_view_count=len(VIEWS),current_actual_microswitch_view_counts=counts,
        current_actual_microswitch_physical_rows_checked=sum(value['rows'] for value in counts.values()),
        current_actual_microswitch_nonzero_physical_contributions_checked=sum(value['nonzero'] for value in counts.values()),
        current_actual_phase1_R2_tensor_interface_rows=joins,current_actual_fresh_phase1_R2_tensor_interface_rows=fresh,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],
        actual_current_completed_tensor_adjacent_interface_count=28,actual_current_completed_tensor_internal_interface_count=10,
        current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),current_actual_microswitch_source_and_tensor_theorem=field.proof,
        original_output_unchanged_except_exposed_unresolved_source_rows_checked=unchanged,
        whole_two_original_microscopic_switches_covered=True,all_original_signed_factored_source_sectors_retained=True,
        inverse_width_not_materialized_and_caps_not_selected_as_source=True,
        exact_original_phase1_controls_and_R2_135_source_row_identities_checked=True,
        invalid_domains_and_foreign_current_owners_rejected=rejected,
        scope='29 actual regions/28 adjacent/10 internal tensor traces. Whole actual first/second microscopic switches plus completed phase1/R2 joins through power/R110/reshape/reference/restore/patch to full Gamma. Core/bridge/R100/axis/angular internal/global/cone/time/energy/points/n-dependent recursion open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},all_passed=True,
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual two micro-switch full tensors/phase1/R2 joins PASS;29 regions/28 adjacent/10 internal;core/global/time open',flush=True)
    return result


if __name__=='__main__':run()
