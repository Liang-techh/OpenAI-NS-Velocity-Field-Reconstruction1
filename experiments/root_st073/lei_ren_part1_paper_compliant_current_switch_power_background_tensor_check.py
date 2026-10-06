"""Admit actual whole R2-to-R110 tensor and exact completed R110 trace."""
import copy
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_switch_power_background_tensor import (
    CurrentSwitchPowerBackgroundTensor,VIEWS,DOMAINS,SEAMS,GATES,OPEN,NAME,RECEIPT,HERE,
    sha,pack,encode,endpoints,source_precision,canonical_tensor_groups,_verify_hashes,read_producer)
from lei_ren_part1_paper_compliant_current_pulse_end_background_tensor_check import physical_rows,check_common
from lei_ren_part1_paper_compliant_current_angular_background_stress_check import finite


def check_view(view,c):
    if view['chart']!='switch_power' or any(view[k] for k in OPEN) or not all(view[k] for k in (
            'actual_full_stress_not_local_difference','source_bounds_not_resolved_physical_point_values',
            'full_current_R2_histories_V110_F0_analytic_P0_and_formal_hb_radius_retained',
            'full_energy_baseline_and_absolute_P0_retained','derivative_coordinate_is_ordinary_logR_not_selector_fraction',
            'source_locals_exposed_only_without_operator_changes','original_Bell_products_combined_before_source_caps',
            'only_radius_and_constant_Pstar_factored_from_full_stress')):raise ValueError('Actual switch-power tensor source or scope lost')
    power=view['actual_upstream_current_switch_power_source'];parent=power['actual_inherited_axial5_packet']
    if not all(power[k] for k in ('actual_R2_histories_retained','complete_postswitch_power_installed',
            'source_caps_not_radius_or_velocity_values','inverse_hb_not_materialized',
            'positive_source_cap_is_enclosure_only','large_derivative_factors_combined_before_positive_source_cap')):raise ValueError('Actual original power histories/source/caps reset')
    if not parent['actual_R2_histories_retained'] or not parent['original_pressure_datum_retained'] or not parent['exposed_existing_locals_are_source_enclosures_not_point_values']:raise ValueError('Actual original power locals or analytic pressure lost')
    if len(parent['actual_normalized_moment_shape_axial5_coefficients'])!=6 or any(len(row)!=6 for row in parent['actual_normalized_moment_shape_axial5_coefficients'].values()):raise ValueError('Full six original power moment shapes required')
    for key in ('actual_postswitch_phi_axial5','actual_postswitch_V_axial5','original_P0_axial5'):
        if len(power[key])!=6:raise ValueError('Actual original phi/V/P0 axial5 lost')
        for value in power[key]:finite(value)
    for value in power['actual_R2_parent_axial5_packet']['pressure_axis_axial5_coefficients']:finite(value)
    for key in ('width_enclosure_is_not_source','exact_positive_width_log','zeta_enclosure_only','R_enclosure_only','theta_enclosure_only'):finite(power[key])
    if not power['formal_log_radius_tree'] or endpoints(power['width_enclosure_is_not_source'])[0]<0 or endpoints(power['theta_enclosure_only'])[0]<=0 or endpoints(power['theta_enclosure_only'])[1]>1:raise ValueError('Original formal positive width/radius source not retained')
    for group in list(view['current_raw_five_history_rows'].values())+[view['current_absolute_pressure_ordinary_y_rows']]:
        if len(group)!=5 or any(row.order!=5 for row in group):raise ValueError('Actual power ordinary y/axial5 raw histories incomplete')
        for row in group:
            for value in row.coefficients:finite(value)
    for label,rows in view['current_source_three_component_velocity_rows'].items():
        if len(rows)!=5 or any(row.order<(4 if label=='radial' else 5) for row in rows):raise ValueError('Full actual power velocity rows lost')
        for row in rows:
            for value in row.coefficients:finite(value)
    amplitude=view['actual_power_amplitude_source_enclosures']
    log_y=amplitude['original_variable_log_amplitude_ordinary_y_axial5']
    if len(log_y)!=4 or any(row.order!=5 for row in log_y):raise ValueError('Original power log jets incomplete')
    if endpoints(log_y[0][0])!=endpoints(c.mpf('.1')) or any(endpoints(value)!=(mp.mpf(0),mp.mpf(0)) for value in log_y[0].coefficients[1:]) or any(endpoints(value)!=(mp.mpf(0),mp.mpf(0)) for row in log_y[1:] for value in row.coefficients):raise ValueError('Original actual power amplitude rate changed')
    if not amplitude['source_derivative_factors_combined_before_positive_caps'] or not amplitude['positive_caps_enclose_original_exponential_not_define_it']:raise ValueError('Actual power cap selected as source or bounded too early')
    grid=view['actual_upstream_physical_spatial4_time1_packet']['physical_spatial_cartesian_mixed4']
    if len(grid)!=35 or any(set(row)!= {'ux','uy','uz','p'} for row in grid.values()):raise ValueError('Actual power BASE spatial4/time1 packet incomplete')
    sectors=view['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']
    if len(sectors['theta'])!=5 or len(sectors['axial'])!=6:raise ValueError('Full actual power meridional stress layout incomplete')
    for parts in sectors.values():
        for part in parts.values():
            if len(part['full_stress_mixed3_coefficient_enclosures'])!=10:raise ValueError('Actual power full stress mixed3 lost')
            for value in part['full_stress_mixed3_coefficient_enclosures'].values():finite(value)
    rows=physical_rows(view);nonzero=0
    if len(rows)!=348 or len(canonical_tensor_groups(view))!=71:raise ValueError('Actual power full tensor/decomposition layout incomplete')
    for row in rows:
        for key in ('signed_coefficient','physical_lambda_exponent','physical_viscosity_exponent','radial_log_prefactor'):finite(row[key])
        for value in row['actual_source_log_parts'].values():finite(value)
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None or endpoints(row['signed_coefficient'])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Zero actual source has nonzero physical bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
        if not row['positive_source_factors_not_materialized'] or not row['source_factors_combined_before_enclosure']:raise ValueError('Exact actual power source factors lost')
    if not any(not row['exact_zero'] for parts in view['physical_three_component_remainder_mixed2']['radial'].values() for row in parts.values()):raise ValueError('Actual inherited radial power remainder deleted')
    return dict(rows=len(rows),nonzero=nonzero)


@source_precision
def run(field=None):
    raw=read_producer();_verify_hashes(raw)
    field=field if field is not None else CurrentSwitchPowerBackgroundTensor(require_checked=False)
    field.assert_graph();omitted={'current_actual_switch_power_tensor_views','current_actual_R110_tensor_interfaces'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:raise ValueError('Actual switch-power manifest/source differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_actual_switch_power_tensor_views'])!=set(VIEWS):raise ValueError('Actual power producer views/scope differ')
    source=field.proof['current_actual_switch_power_source_pressure_and_R110_theorem']
    boundary=source['actual_original_R110_power_post_and_reshape_boundary_function_theorem']
    if not source['passed'] or not all(source['live_original_callable_bindings'].values()) or not boundary['passed'] or boundary['total_physical_mixed4_rows_implied']!=135 or len(boundary['physical_source_function_identities'])!=9:raise ValueError('Actual original power/P0/R110 source proof failed')
    counts={};unchanged=[];variable=False
    for name,args in VIEWS.items():
        value=field.chart(*args)
        if encode(pack(value))!=raw['current_actual_switch_power_tensor_views'][name]:raise ValueError('Actual power full source replay differs: '+name)
        counts[name]=check_view(value,field.ctx)
        if name in ('switch_power_whole','switch_power_fresh'):
            exposed=value['actual_upstream_current_switch_power_source']
            original=field.switch.postpower(args[1],args[2])
            if encode(pack({k:v for k,v in exposed.items() if k!='actual_inherited_axial5_packet'}))!=encode(pack(original)):raise ValueError('Original power output changed beyond exposed locals')
            unchanged.append(name)
        if name=='switch_power_whole':
            lo,hi=endpoints(value['actual_upstream_current_switch_power_source']['zeta_enclosure_only']);variable=lo<hi
        if name=='switch_power_right' and endpoints(value['actual_upstream_current_switch_power_source']['R_enclosure_only'])!=(mp.mpf(110),mp.mpf(110)):raise ValueError('Exact R110 endpoint offset not cancelled')
    if not variable:raise ValueError('Whole original power interval replaced by one endpoint')
    joins={};fresh={}
    for name in SEAMS:
        value=field.interface(name)
        if encode(pack(value))!=raw['current_actual_R110_tensor_interfaces'][name]:raise ValueError('Actual R110 completed tensor trace differs')
        joins[name]=check_common(value,71)
        fresh[name]=check_common(field.interface(name,Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2'),71)
    rejected=[]
    for label,args in (('before_R2',('switch_power','.2','-.01')),('past_R110',('switch_power','.2','1.01')),
            ('outside_Z',('switch_power',2,'.5')),('foreign_chart',('switch_second','.2','1.5')),
            ('nonfinite_coordinate',('switch_power','.2','inf'))):
        try:field.chart(*args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual power domain admitted')
    for label,kwargs in (('nonfinite_time',dict(log_tau='inf')),('nonpositive_viscosity',dict(viscosity=0)),('nonfinite_theta',dict(theta='inf'))):
        try:field.chart('switch_power','.2','.5',**kwargs)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual power physical sector admitted')
    try:field.interface('R2')
    except ValueError:rejected.append('unbuilt_R2_tensor_join')
    else:raise ValueError('Unbuilt microswitch tensor join admitted')
    clone=copy.copy(field);clone.switch=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_current_switch_owner')
    else:raise ValueError('Foreign actual switch owner admitted')
    clone=copy.copy(field);clone.reshape_tensor=copy.copy(field.reshape_tensor);clone.reshape_tensor.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_current_reshape_parent')
    else:raise ValueError('Unchecked current reshape tensor admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_actual_switch_power_tensor_view_count=len(VIEWS),current_actual_switch_power_view_counts=counts,
        current_actual_switch_power_physical_rows_checked=sum(v['rows'] for v in counts.values()),
        current_actual_switch_power_nonzero_physical_contributions_checked=sum(v['nonzero'] for v in counts.values()),
        current_actual_R110_tensor_interface_rows=joins,current_actual_fresh_R110_tensor_interface_rows=fresh,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],
        actual_current_completed_tensor_adjacent_interface_count=26,actual_current_completed_tensor_internal_interface_count=10,
        current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),current_actual_switch_power_source_and_tensor_theorem=field.proof,
        original_output_unchanged_except_existing_source_locals_checked=unchanged,
        whole_R2_to_R110_formal_geometry_and_five_histories_retained=True,
        exact_original_R110_geometry_moments_amplitude_P0_function_identity_checked=True,
        invalid_domains_and_foreign_current_owners_rejected=rejected,
        scope='27 actual regions/26 adjacent/10 internal tensor traces. Whole actual R2..R110 power and completed R110 join through reshape/Rsh/reference/restore/patch to full Gamma; upstream core/bridge/microswitch/R100/R2/axis/angular internal/global/cone/time/energy/points/n-dependent recursion open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},all_passed=True,
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual R2..R110 power full tensor and R110 join PASS;27 regions/26 adjacent/10 internal;core/global/time open',flush=True)
    return result


if __name__=='__main__':run()
