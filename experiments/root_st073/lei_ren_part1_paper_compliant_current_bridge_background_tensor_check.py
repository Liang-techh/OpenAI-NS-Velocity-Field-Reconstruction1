"""Admit current actual whole bridge tensors and three completed source joins."""
import copy
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_bridge_background_tensor import (
    CurrentBridgeBackgroundTensor,VIEWS,DOMAINS,SEAMS,SHARDS,GATES,OPEN,NAME,RECEIPT,HERE,
    sha,pack,encode,endpoints,source_precision,canonical_tensor_groups,_verify_hashes,read_producer,read_shard)
from lei_ren_part1_paper_compliant_current_microswitch_background_tensor_check import (
    THETA,AXIAL,check_factored_record,physical_rows,check_common,finite)


def check_view(view,c):
    if view['chart'] not in DOMAINS or any(view[k] for k in OPEN) or not all(view[k] for k in (
            'microscopic_inverse_width_or_macro_identity_before_any_resolution','every_original_stress_and_remainder_sector_expanded_without_pruning',
            'original_radial_prefactors_and_absolute_P0_included_once',
            'full_current_signed_axial_drive_axial6_comparison_and_actual_own_histories_retained',
            'derivative_coordinate_is_ordinary_logR_not_macro_fraction','source_factors_retained_through_full_tensor_and_physical_operators',
            'source_locals_exposed_only_without_operator_changes','actual_full_stress_not_local_difference','source_bounds_not_resolved_physical_point_values')):
        raise ValueError('Actual full bridge source/scope lost')
    original=view['actual_upstream_current_bridge_source'];macro=view['chart']=='bridge_macro'
    if original['chart']!=view['chart'].replace('bridge_','') or not original['current_coordinate_labelled_finite_width_history_used'] or not original['actual_moments_not_replaced_by_comparison'] or not original['all_source_width_and_amplitude_factors_retained_until_final_rows']:raise ValueError('Actual coordinate/own-history bridge source lost')
    if original['comparison_radially_frozen']!=macro:raise ValueError('True smoothing comparison was replaced by frozen values')
    parent=original['actual_parent_axial5_packet']
    if parent['comparison_moments_substituted'] or not parent['current_actual_integral_and_own_feedback_source'] or parent['exact_coordinate_labelled_source']!='upstream.packet(Z,value,chart)':raise ValueError('Same signed actual six-history prefix source required')
    for value in parent['pressure_axis_axial5_coefficients']:finite(value)
    for rows in original['comparison_parent_axial6_packet'].values():
        if len(rows)!=7:raise ValueError('Full axial-six comparison source required')
        for value in rows:finite(value)
    if not original['source_radius_tree'] or not original['exact_comparison_function_source']:raise ValueError('Original formal radius/comparison source trees lost')
    for key in ('source_width_log','width_enclosure_is_not_source','y_enclosure_only','R_enclosure_only','theta_enclosure_only'):finite(original[key])
    logs=view['fixed_current_factored_source_log_bases']
    if len(logs)!=4 or endpoints(logs[0])!=endpoints(original['source_width_log']):raise ValueError('Original four bridge factor log bases lost')
    for value in logs:finite(value)
    raw=view['current_unresolved_raw_source_rows'];source=view['original_unresolved_source_coordinate_rows'];terms=0
    if set(raw['histories'])!={'m','h','k','e','p'} or set(raw['velocity'])!={'radial','theta','axial'}:raise ValueError('Complete bridge raw units lost')
    if raw['inverse_width_source_shift_applied_before_any_resolution']==macro or raw['macro_original_ordinary_y_preserved']!=macro:raise ValueError('Wrong source-coordinate width conversion')
    expected='D_y^j=D_original_coordinate^j; fraction labels coverage only' if macro else 'D_y^j=hb^-j*D_phase^j'
    if view['ordinary_logR_conversion']!=expected:raise ValueError('Macro fraction substituted for ordinary y')
    for rows in list(raw['histories'].values())+list(raw['velocity'].values())+[raw['absolute_pressure']]:
        if len(rows)!=5:raise ValueError('Full bridge ordinary y0..4 required')
        for row in rows:terms+=check_factored_record(row)
    if len(source['physical'])!=4 or len(source['primitives'])!=5 or len(source['Q'])!=5:raise ValueError('Original full unresolved bridge source rows lost')
    for rows in list(source['physical'].values())+list(source['primitives'].values())+[source['Q']]:
        if len(rows)!=5:raise ValueError('Original full coordinate source rows lost')
        for row in rows:check_factored_record(row)
    if not raw['velocity']['axial'][1]['terms']:raise ValueError('Original nonzero signed bridge axial drive removed')
    sectors=view['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']
    if {part['original_stress_sector'] for part in sectors['theta'].values()}!=THETA or {part['original_stress_sector'] for part in sectors['axial'].values()}!=AXIAL:raise ValueError('Original full eleven stress sectors omitted')
    for parts in sectors.values():
        for part in parts.values():
            if len(part['full_stress_mixed3_coefficient_enclosures'])!=10:raise ValueError('Full bridge stress mixed3 lost')
            for value in part['full_stress_mixed3_coefficient_enclosures'].values():finite(value)
            for value in part['exact_source_log_parts'].values():finite(value)
    grid=view['actual_upstream_physical_spatial4_time1_packet']['physical_spatial_cartesian_mixed4']
    if len(grid)!=35 or any(set(row)!={'ux','uy','uz','p'} for row in grid.values()):raise ValueError('Same current BASE spatial4/time1 packet incomplete')
    rows=physical_rows(view);nonzero=0
    if len(rows)<348 or len(canonical_tensor_groups(view))!=71:raise ValueError('Complete bridge tensor/decomposition inventory lost')
    for row in rows:
        for key in ('signed_coefficient','physical_lambda_exponent','physical_viscosity_exponent','radial_log_prefactor'):finite(row[key])
        for value in row['actual_source_log_parts'].values():finite(value)
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None or endpoints(row['signed_coefficient'])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Zero bridge source has nonzero bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
        if not row['positive_source_factors_not_materialized'] or not row['source_factors_combined_before_enclosure']:raise ValueError('Original bridge source factor resolved too early')
    if not any(not row['exact_zero'] for parts in view['physical_three_component_remainder_mixed2']['radial'].values() for row in parts.values()):raise ValueError('Full nonzero radial bridge remainder removed')
    return dict(rows=len(rows),nonzero=nonzero,unresolved_raw_factored_terms=terms,
        theta_source_sectors=len(sectors['theta']),axial_source_sectors=len(sectors['axial']))


@source_precision
def run(field=None):
    raw=read_producer();_verify_hashes(raw)
    field=field if field is not None else CurrentBridgeBackgroundTensor(require_checked=False);field.assert_graph()
    omitted={'current_actual_bridge_chart_shards','current_actual_bridge_tensor_interfaces'}
    if encode(pack(field.manifest()))!={key:value for key,value in raw.items() if key not in omitted}:raise ValueError('Actual bridge source manifest differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_actual_bridge_chart_shards'])!=set(SHARDS):raise ValueError('Actual bridge producer scope/layout differs')
    proof=field.proof['current_actual_bridge_source_pressure_units_and_three_endpoint_theorem']
    unit=proof['original_current_micro_and_macro_raw_unit_theorem'];boundary=proof['original_current_full_boundary_control_and_generator_theorem']
    if not proof['passed'] or not all(proof['live_original_callable_bindings'].values()) or unit['total_identities']!=282 or not all(all(rows.values()) for rows in unit['microscopic_and_macro_original_source_unit_identities'].values()) or boundary['source_generator_identity_count']!=270 or not all(all(rows.values()) for rows in boundary['full_original_second_macro_and_R100_source_generator_identities'].values()):raise ValueError('Original complete bridge unit/boundary operator theorem failed')
    counts={};unchanged=[];variable={};datahashes={}
    for chart,path in SHARDS.items():
        descriptor=raw['current_actual_bridge_chart_shards'][chart]
        expected={name:args for name,args in VIEWS.items() if args[0]==chart}
        if descriptor['path']!=path or descriptor['sha256']!=sha(path) or descriptor['views']!=list(expected):raise ValueError('Complete bridge chart shard identity changed')
        datahashes[path]=sha(path);shard=read_shard(chart)
        if set(shard)!=set(expected):raise ValueError('Unpruned six-sector bridge chart shard required')
        for name,args in expected.items():
            value=field.chart(*args)
            if encode(pack(value))!=shard[name]:raise ValueError('Actual full bridge source replay differs: '+name)
            counts[name]=check_view(value,field.ctx)
            if name.endswith('_whole') or name.endswith('_fresh'):
                original=field.bridge.evaluate(args[1],args[2],chart.replace('bridge_',''))
                if encode(pack(value['actual_upstream_current_bridge_source']))!=encode(pack(original)):raise ValueError('Original current bridge output changed beyond exposed locals')
                unchanged.append(name)
            if name.endswith('_whole'):
                lo,hi=endpoints(value['coverage_coordinate']);variable[name]=lo<hi
            if name=='macro_right' and endpoints(value['actual_upstream_current_bridge_source']['R_enclosure_only'])!=(mp.mpf(100),mp.mpf(100)):raise ValueError('Exact original R100 endpoint changed')
            print('Check actual full bridge tensor: '+name,flush=True)
    if len(variable)!=3 or not all(variable.values()):raise ValueError('Whole original bridge charts replaced by endpoint samples')
    joins={};fresh={}
    for name in SEAMS:
        value=field.interface(name)
        if encode(pack(value))!=raw['current_actual_bridge_tensor_interfaces'][name]:raise ValueError('Actual complete bridge tensor trace differs')
        joins[name]=check_common(value,71);fresh[name]=check_common(field.interface(name,Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2'),71)
    rejected=[]
    for label,args in (('before_first',('bridge_first','.2','-.01')),('past_first',('bridge_first','.2','1.01')),
            ('before_second',('bridge_second','.2','.99')),('past_second',('bridge_second','.2','2.01')),
            ('before_macro',('bridge_macro','.2','-.01')),('past_macro',('bridge_macro','.2','1.01')),
            ('outside_Z',('bridge_first',2,'.5')),('foreign_chart',('core','.2','.5')),('nonfinite_coordinate',('bridge_first','.2','inf'))):
        try:field.chart(*args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual bridge domain admitted')
    for label,kwargs in (('nonfinite_time',dict(log_tau='inf')),('nonpositive_viscosity',dict(viscosity=0)),('nonfinite_theta',dict(theta='inf'))):
        try:field.chart('bridge_first','.2','.5',**kwargs)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid bridge physical sector admitted')
    try:field.interface('core_bridge')
    except ValueError:rejected.append('unbuilt_core_tensor_join')
    else:raise ValueError('Unbuilt core tensor interface admitted')
    clone=copy.copy(field);clone.bridge=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_current_bridge_owner')
    else:raise ValueError('Foreign actual bridge owner admitted')
    clone=copy.copy(field);clone.micro_tensor=copy.copy(field.micro_tensor);clone.micro_tensor.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_current_micro_parent')
    else:raise ValueError('Unchecked current micro tensor admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_actual_bridge_tensor_view_count=len(VIEWS),current_actual_bridge_view_counts=counts,
        current_actual_bridge_physical_rows_checked=sum(value['rows'] for value in counts.values()),
        current_actual_bridge_nonzero_physical_contributions_checked=sum(value['nonzero'] for value in counts.values()),
        current_actual_bridge_tensor_interface_rows=joins,current_actual_fresh_bridge_tensor_interface_rows=fresh,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],
        actual_current_completed_tensor_adjacent_interface_count=31,actual_current_completed_tensor_internal_interface_count=10,
        current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),current_actual_bridge_source_and_tensor_theorem=field.proof,
        original_output_unchanged_except_exposed_unresolved_source_rows_checked=unchanged,
        whole_three_original_bridge_charts_covered=True,all_original_signed_factored_source_sectors_retained=True,
        macro_original_logR_units_and_micro_true_inverse_width_retained=True,
        current_three_functional_join_receipt_and_full_boundary_source_programs_consumed=True,
        invalid_domains_and_foreign_current_owners_rejected=rejected,
        scope='32 actual regions/31 adjacent/10 internal tensor traces. Whole actual bridge-first/second/macro and completed phase1/smoothing/R100 joins through microswitch/power/R110/reshape/reference/restore/patch to full Gamma. Core/core-bridge/axis/angular internal/global/cone/time/energy/points/n-dependent recursion open.',
        input_hashes={**raw['input_hashes'],**datahashes,NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual three full bridge tensors/phase1/smoothing/R100 joins PASS;32 regions/31 adjacent/10 internal;core/global/time open',flush=True)
    return result


if __name__=='__main__':run()
