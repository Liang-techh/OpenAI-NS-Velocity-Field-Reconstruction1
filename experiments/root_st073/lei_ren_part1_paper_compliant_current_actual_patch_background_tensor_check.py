"""Admit current five-moment patch tensor, actual Rh join and six flat traces."""
import copy
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_actual_patch_background_tensor import (
    CurrentActualPatchBackgroundTensor,VIEWS,SEAMS,EDGES,GATES,OPEN,NAME,RECEIPT,HERE,
    sha,pack,encode,endpoints,source_precision,canonical_tensor_groups,_verify_hashes,read_producer)
from lei_ren_part1_paper_compliant_current_pulse_end_background_tensor_check import physical_rows,check_common
from lei_ren_part1_paper_compliant_current_angular_background_stress_check import finite


def check_view(view):
    if view['chart']!='actual_patch' or any(view[k] for k in OPEN) or not all(view[k] for k in (
            'actual_full_stress_not_local_difference','source_bounds_not_resolved_physical_point_values',
            'full_current_implicit_axial5_and_partial_integrals_retained','full_energy_cross_terms_and_absolute_P0_retained',
            'ordinary_x_converted_to_logR_before_full_physical_derivatives',
            'only_radius_and_constant_Pstar_factored_from_full_stress','rounded_support_coordinate_is_enclosure_not_exact_zero_source')):raise ValueError('Actual current patch tensor scope/source lost')
    patch=view['actual_upstream_current_implicit_patch_source'];parent=patch['actual_inherited_patch_packet']
    if not (patch['same_actual_coefficient_family_and_P0_retained'] and patch['ordinary_derivatives_not_Taylor_radial_coefficients']
            and parent['actual_Rm_history_retained'] and parent['axis_pressure_not_changed']):raise ValueError('Original actual implicit patch owner/history lost')
    if len(parent['actual_coefficient_functions'])!=5 or any(row.order!=5 for row in parent['actual_coefficient_functions']):raise ValueError('Current implicit axial5 coefficient functions missing')
    for row in parent['actual_coefficient_functions']:
        for value in row.coefficients:finite(value)
    if len(parent['original_P0_axial5'])!=6 or len(parent['actual_normalized_five_defects'])!=5:raise ValueError('Actual P0/five partial histories missing')
    for value in parent['actual_partial_weights'].values():finite(value)
    groups=list(view['current_raw_five_history_rows'].values())+[view['current_absolute_pressure_ordinary_y_rows']]
    for group in groups:
        if len(group)!=5 or any(row.order!=5 for row in group):raise ValueError('Actual patch full ordinary y/axial5 histories required')
        for row in group:
            for value in row.coefficients:finite(value)
    for label,rows in view['current_source_three_component_velocity_rows'].items():
        if len(rows)!=5 or any(row.order<(4 if label=='radial' else 5) for row in rows):raise ValueError('Full actual patch velocity rows lost')
        for row in rows:
            for value in row.coefficients:finite(value)
    grid=view['actual_upstream_physical_spatial4_time1_packet']['physical_spatial_cartesian_mixed4']
    if len(grid)!=35 or any(set(row)!= {'ux','uy','uz','p'} for row in grid.values()):raise ValueError('Actual patch BASE spatial4/time1 packet incomplete')
    sectors=view['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']
    if len(sectors['theta'])!=5 or len(sectors['axial'])!=6:raise ValueError('Full actual patch meridional stress layout incomplete')
    for parts in sectors.values():
        for part in parts.values():
            if len(part['full_stress_mixed3_coefficient_enclosures'])!=10:raise ValueError('Actual full patch stress mixed3 lost')
            for value in part['full_stress_mixed3_coefficient_enclosures'].values():finite(value)
    rows=physical_rows(view);nonzero=0
    if len(rows)!=348 or len(canonical_tensor_groups(view))!=71:raise ValueError('Actual patch full tensor/decomposition layout incomplete')
    for row in rows:
        for key in ('signed_coefficient','physical_lambda_exponent','physical_viscosity_exponent','radial_log_prefactor'):finite(row[key])
        for value in row['actual_source_log_parts'].values():finite(value)
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None or endpoints(row['signed_coefficient'])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Zero actual source has nonzero bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
        if not row['positive_source_factors_not_materialized'] or not row['source_factors_combined_before_enclosure']:raise ValueError('Exact actual patch source factors lost')
    remainder=view['physical_three_component_remainder_mixed2']
    if not any(not row['exact_zero'] for parts in remainder['radial'].values() for row in parts.values()):raise ValueError('Actual radial patch remainder deleted')
    return dict(rows=len(rows),nonzero=nonzero)


@source_precision
def run(field=None):
    raw=read_producer();_verify_hashes(raw)
    field=field if field is not None else CurrentActualPatchBackgroundTensor(require_checked=False)
    field.assert_graph();omitted={'current_actual_patch_tensor_views','current_actual_patch_Rh_and_six_support_tensor_interfaces'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:raise ValueError('Actual patch manifest/source differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_actual_patch_tensor_views'])!=set(VIEWS):raise ValueError('Actual patch producer views/scope differ')
    source=field.proof['current_actual_patch_source_pressure_support_and_Rh_theorem']
    normalization=field.proof['original_actual_patch_current_radius_normalization_theorem']
    if not source['passed'] or not all(source['live_original_callable_bindings'].values()) or not normalization['passed'] or not all(normalization['identities'].values()):raise ValueError('Actual current patch source or current-radius units failed')
    if len(normalization['identities'])!=141 or len(source['six_exact_rational_support_edges'])!=6 or not source['actual_current_pressure_defining_function_bridge']['passed']:raise ValueError('Full actual source/P0/support proof missing')
    counts={};active=[]
    for name,args in VIEWS.items():
        value=field.chart(*args)
        if encode(pack(value))!=raw['current_actual_patch_tensor_views'][name]:raise ValueError('Actual patch full source replay differs: '+name)
        counts[name]=check_view(value);parent=value['actual_upstream_current_implicit_patch_source']['actual_inherited_patch_packet']
        if name.startswith('active_'):
            if not any(endpoints(row[0])!=(mp.mpf(0),mp.mpf(0)) for row in value['actual_upstream_current_implicit_patch_source']['actual_gamma_x_derivatives']):raise ValueError('Actual active beta derivatives deleted')
            if not any(not row['exact_zero'] for parts in value['physical_three_component_remainder_mixed2']['axial'].values() for row in parts.values()):raise ValueError('Actual active implicit-patch axial viscosity deleted')
            active.append(name)
        if name=='patch_Rh' and not (parent['exact_terminal_closure_from_same_implicit_map'] and parent['terminal_refinement_is_functional_identity_not_moment_reset']):raise ValueError('Current Rh functional implicit closure missing')
        if name.startswith('patch_edge_') and not value['rounded_support_coordinate_is_enclosure_not_exact_zero_source']:raise ValueError('Rounded edge promoted to zero source')
    if len(active)!=3:raise ValueError('All three original active support families required')
    joins={};fresh={}
    for name in SEAMS:
        value=field.interface(name)
        if encode(pack(value))!=raw['current_actual_patch_Rh_and_six_support_tensor_interfaces'][name]:raise ValueError('Actual patch completed tensor trace differs')
        joins[name]=check_common(value,71)
        fresh[name]=check_common(field.interface(name,Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2'),71)
    rejected=[]
    for label,args in (('below_Rm',('actual_patch','.2','.99')),('above_Rh',('actual_patch','.2','2.72')),
            ('outside_Z',('actual_patch',2,'1.25')),('foreign_chart',('Rh_reference','.2','1.25')),
            ('nonfinite_coordinate',('actual_patch','.2','inf')),('foreign_support_edge',('actual_patch','.2','edge:50'))):
        try:field.chart(*args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual patch source domain admitted')
    for label,kwargs in (('nonfinite_time',dict(log_tau='inf')),('nonpositive_viscosity',dict(viscosity=0)),('nonfinite_theta',dict(theta='inf'))):
        try:field.chart('actual_patch','.2','1.25',**kwargs)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid actual patch physical sector admitted')
    clone=copy.copy(field);clone.patch=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_actual_patch_owner')
    else:raise ValueError('Foreign actual patch owner admitted')
    clone=copy.copy(field);clone.o2_tensor=copy.copy(field.o2_tensor);clone.o2_tensor.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_current_O2_parent')
    else:raise ValueError('Unchecked current O2 tensor admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_actual_patch_tensor_view_count=len(VIEWS),current_actual_patch_view_counts=counts,
        current_actual_patch_physical_rows_checked=sum(v['rows'] for v in counts.values()),
        current_actual_patch_nonzero_physical_contributions_checked=sum(v['nonzero'] for v in counts.values()),
        current_actual_patch_Rh_and_six_support_tensor_interface_rows=joins,
        current_actual_patch_Rh_and_six_fresh_support_tensor_interface_rows=fresh,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],
        actual_current_completed_tensor_adjacent_interface_count=21,actual_current_completed_tensor_internal_interface_count=10,
        current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),
        current_actual_patch_source_units_and_tensor_theorem=field.proof,
        current_implicit_axial5_partial_integrals_full_energy_cross_terms_and_same_P0_retained=True,
        active_original_bump_families_and_full_axial_viscosity_checked=active,
        exact_functional_flat_limits_not_rounded_zero_boxes_used_for_six_internal_traces=True,
        invalid_domains_and_foreign_current_owners_rejected=rejected,
        scope='22 actual regions/21 adjacent/10 internal tensor traces. Whole actual implicit five-moment patch, actual Rh completed tensor attachment and six support traces; upstream core/bridge/switch/reshape/restore/axis/angular internal/global/cone/time/energy/points/n-dependent recursion open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},all_passed=True,
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual implicit patch/Rh/six support tensors PASS;22 regions/21 adjacent/10 internal;core/global/time open',flush=True)
    return result


if __name__=='__main__':run()
