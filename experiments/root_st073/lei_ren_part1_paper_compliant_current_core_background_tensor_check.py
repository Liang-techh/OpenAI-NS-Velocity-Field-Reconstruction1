"""Independent scoped admission of full positive-radius core T/E and core/first trace."""
import copy
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_compliant_current_core_background_tensor import (
    CurrentCoreBackgroundTensor,VIEWS,GATES,OPEN,NAME,RECEIPT,HERE,sha,pack,encode,
    endpoints,source_precision,canonical_tensor_groups,_verify_hashes,read_producer)
from lei_ren_part1_paper_compliant_current_microswitch_background_tensor_check import (
    THETA,AXIAL,check_factored_record,physical_rows,finite)
from lei_ren_part1_paper_compliant_current_core_stress_operator import (
    core_raw_unit_theorem,core_integrated_stress_theorem)
from lei_ren_part1_paper_compliant_current_core_interior_moments_check import check_view as check_moments


def leaves(value):
    if isinstance(value,dict) and 'signed_coefficient' in value:return [value]
    if isinstance(value,dict):return sum((leaves(v) for v in value.values()),[])
    if isinstance(value,(list,tuple)):return sum((leaves(v) for v in value),[])
    return []


def check_common(view,count):
    if any(view[k] for k in OPEN) or len(view['common_actual_tensor_rows'])!=count or view['current_common_tensor_contribution_count']!=count:
        raise ValueError('Current positive-radius core/first tensor trace scope/layout differs')
    if not view['source_function_equality_precedes_common_triangle_bounds'] or not view['interval_overlap_not_used_as_function_identity']:
        raise ValueError('Common source equality must precede bounds')
    nonzero=0
    for row in view['common_actual_tensor_rows'].values():
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None:raise ValueError('Zero common trace has bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
    return dict(rows=count,nonzero=nonzero)


def check_view(view,c):
    if view['chart']!='core' or any(view[k] for k in OPEN) or not all(view[k] for k in (
        'exact_zero_total_stress_is_separate_from_original_sector_ledger',
        'every_original_NS_remainder_sector_retained','original_radial_prefactors_and_absolute_P0_included_once',
        'actual_full_stress_not_local_difference','source_bounds_not_resolved_physical_point_values',
        'physical_axis_not_in_this_chart','uniform_axis_or_required_domain_energy_not_inferred')):
        raise ValueError('Actual positive-radius core/source scope lost')
    lo,hi=endpoints(view['coverage_coordinate'])
    if lo<=0 or hi>4:raise ValueError('Axis or foreign core radius admitted')
    source=view['actual_same_fixed_point_core_source']
    check_moments(source['original_same_fixed_point_density_and_moment_packet'],c)
    if not source['no_phase_width_source_or_inverse_width_in_core'] or not source['fixed_basepoint_F0_bound_encloses_source_not_selects_value']:
        raise ValueError('Original fixed core factor source lost')
    if source['formal_radius']!='epsilon_core*rho' or source['formal_unit']!='sqrt(2*R)*F0base/Pstar':
        raise ValueError('Original exact core radius/swirl units changed')
    raw=view['current_unresolved_raw_source_rows'];terms=0
    if set(raw['velocity'])!={'radial','theta','axial'} or set(raw['histories'])!={'m','h','k','e','p'}:
        raise ValueError('Full original core raw source units lost')
    for rows in list(raw['velocity'].values())+list(raw['histories'].values())+[raw['absolute_pressure']]:
        if len(rows)!=5:raise ValueError('Full ordinary core y0..4 required')
        for row in rows:terms+=check_factored_record(row)
    for flag in ('original_core_h_and_k_unit_half_retained','original_composite_energy_A_minus_swirl_B_retained',
        'original_P0_separate_from_cumulative_C','ordinary_logR_rows_before_any_source_resolution'):
        if not raw[flag]:raise ValueError('Original raw core unit factor lost')
    logs=view['fixed_current_factored_source_log_bases']
    if len(logs)!=4 or endpoints(logs[0])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Microscopic width substituted into core')
    for value in logs:finite(value)
    original=view['original_uncancelled_signed_stress_sector_enclosures']
    if {part['original_stress_sector'] for part in original['theta'].values()}!=THETA or {part['original_stress_sector'] for part in original['axial'].values()}!=AXIAL:
        raise ValueError('Original eleven signed stress sectors lost')
    sector_nonzero=0
    for parts in original.values():
        for part in parts.values():
            if len(part['full_stress_mixed3_coefficient_enclosures'])!=10:raise ValueError('Original stress mixed3 lost')
            for value in part['full_stress_mixed3_coefficient_enclosures'].values():
                finite(value);sector_nonzero+=endpoints(value)!=(mp.mpf(0),mp.mpf(0))
            for value in part['exact_source_log_parts'].values():finite(value)
    if not sector_nonzero:raise ValueError('Individual stress sectors falsely zeroed')
    for label,parts in view['current_actual_source_stress_packet']['full_meridional_stress_log_sectors'].items():
        if label not in ('theta','axial') or len(parts)!=1:raise ValueError('Separate total stress aggregate required')
        aggregate=next(iter(parts.values()))
        if aggregate['original_stress_sector']!='integrated_core_total_stress_zero':raise ValueError('Original sector relabelled as zero')
        if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in aggregate['full_stress_mixed3_coefficient_enclosures'].values()):
            raise ValueError('Exact integrated total stress failed')
    for key in ('physical_cylindrical_stress_mixed3','physical_cylindrical_stress_divergence_mixed2',
        'physical_completed_stress_tensor_cartesian','physical_completed_stress_divergence_cartesian'):
        if any(not row['exact_zero'] for row in leaves(view[key])):raise ValueError('Exact source zero lost in linear completion')
    remainder=view['physical_three_component_remainder_mixed2']
    expected={'radial':{'time','radial_viscosity','nonlinear_transport','axial_viscosity'},
        'theta':{'axial_viscosity'},'axial':{'axial_viscosity'}}
    for label,names in expected.items():
        if {name.split('__source')[0] for name in remainder[label]}!=names:raise ValueError('Original full NS remainder source sector removed')
        for part in remainder[label].values():
            if len(part)!=6:raise ValueError('Full remainder mixed2 lost')
    if not any(not row['exact_zero'] for row in leaves(remainder['radial'])):raise ValueError('Nonzero radial remainder was removed')
    native=view['actual_upstream_physical_spatial4_time1_packet']
    if native['full_point_physical_field_evaluation'] or len(native['requested_core_spatial_log_bounds'])!=35 or set(native['requested_core_time_log_bounds'])!={'ux','uy','uz','p'}:
        raise ValueError('Same original nonsingular core spatial4/time1 source packet incomplete')
    rows=physical_rows(view);nonzero=0
    if len(rows)!=111 or len(canonical_tensor_groups(view))!=71:raise ValueError('Full core completed tensor/remainder contribution layout lost')
    for row in rows:
        for key in ('signed_coefficient','physical_lambda_exponent','physical_viscosity_exponent','radial_log_prefactor'):finite(row[key])
        for value in row['actual_source_log_parts'].values():finite(value)
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None or endpoints(row['signed_coefficient'])!=(mp.mpf(0),mp.mpf(0)):
                raise ValueError('Source zero has nonzero bound')
        else:finite(row['log_absolute_upper']);nonzero+=1
        if not row['positive_source_factors_not_materialized'] or not row['source_factors_combined_before_enclosure']:
            raise ValueError('Core source factor materialized before physical bound')
    return dict(rows=len(rows),nonzero=nonzero,unresolved_raw_factored_terms=terms,
        original_nonzero_stress_coefficients_preserved=sector_nonzero,
        exact_zero_total_stress_with_nonzero_remainder=True)


@source_precision
def run(field=None):
    raw=read_producer();_verify_hashes(raw)
    field=field if field is not None else CurrentCoreBackgroundTensor(require_checked=False);field.assert_graph()
    omitted={'current_positive_radius_core_tensor_views','current_core_bridge_completed_tensor_interface'}
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k not in omitted}:raise ValueError('Current core tensor source manifest differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_positive_radius_core_tensor_views'])!=set(VIEWS):
        raise ValueError('Positive-radius core producer scope differs')
    proof=field.proof['current_same_source_core_units_integrated_stress_and_completed_first_trace']
    unit=core_raw_unit_theorem();integrated=core_integrated_stress_theorem()
    if encode(pack(unit))!=encode(pack(proof['original_core_raw_unit_theorem'])) or encode(pack(integrated))!=encode(pack(proof['original_integrated_total_core_stress_theorem'])):
        raise ValueError('Independent core raw/source stress proof replay differs')
    if not all(proof['live_original_callable_bindings'].values()) or not proof['passed']:
        raise ValueError('Same current nonlinear core source binding failed')
    # Rebuild only this adapter's local source packets; keep the checked shared graph warm.
    field.interior.cache={};field.interior.rows_cache={};field.interior.tail_cache={}
    counts={}
    for name,args in VIEWS.items():
        value=field.chart(*args)
        if encode(pack(value))!=raw['current_positive_radius_core_tensor_views'][name]:raise ValueError('Current full core tensor replay differs: '+name)
        counts[name]=check_view(value,field.ctx)
        print('Check full same-source positive-radius core tensor: '+name,flush=True)
    core_sector=field.ctx.mpf(VIEWS['compact_core'][2])
    if endpoints(core_sector)[0]>endpoints(field.ctx.mpf(VIEWS['fresh_low'][2]))[1]:
        below_published_sector=True
    else:raise ValueError('Fresh low radius is not below published compact sector')
    join=field.interface()
    if encode(pack(join))!=raw['current_core_bridge_completed_tensor_interface']:raise ValueError('Core/first completed tensor trace differs')
    joins=check_common(join,71)
    fresh=field.interface(Z=('-0.8','0.8'),log_tau=('-5','-2'),viscosity='.2')
    fresh_joins=check_common(fresh,71)
    rejected=[]
    for label,args in (('axis_not_admitted',('core','.2',0)),('negative_rho',('core','.2','-.01')),
        ('past_core',('core','.2','4.001')),('outside_Z',('core',2,1)),('foreign_chart',('bridge_first','.2',1)),
        ('nonfinite_rho',('core','.2','inf'))):
        try:field.chart(*args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid core tensor domain admitted')
    for label,kwargs in (('nonfinite_time',dict(log_tau='inf')),('nonpositive_viscosity',dict(viscosity=0)),('nonfinite_theta',dict(theta='inf'))):
        try:field.chart('core','.2',1,**kwargs)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid physical core sector admitted')
    for label,attribute in (('foreign_original_core','core'),('foreign_raw_operator','raw_source')):
        clone=copy.copy(field);setattr(clone,attribute,object())
        try:clone.assert_graph()
        except ValueError:rejected.append(label);continue
        raise ValueError('Foreign actual core/operator accepted')
    clone=copy.copy(field);clone.interior=copy.copy(field.interior);clone.interior.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_core_interior')
    else:raise ValueError('Unchecked core moment prerequisite accepted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_positive_radius_core_tensor_view_count=len(VIEWS),current_positive_radius_core_view_counts=counts,
        current_positive_radius_core_physical_rows_checked=sum(v['rows'] for v in counts.values()),
        current_core_bridge_completed_tensor_interface_rows=joins,current_fresh_core_bridge_completed_tensor_interface_rows=fresh_joins,
        current_positive_radius_core_source_and_tensor_theorem=field.proof,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],
        actual_current_completed_tensor_adjacent_interface_count=32,actual_current_completed_tensor_internal_interface_count=10,
        current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),
        low_positive_radius_below_published_compact_sector_checked=below_published_sector,
        original_signed_core_stress_ledger_retained_and_only_total_cancelled=True,
        all_six_original_NS_remainder_source_sectors_retained=True,
        invalid_domains_and_foreign_current_owners_rejected=rejected,
        scope='33 positive-radius actual tensor regions/32 adjacent/10 internal traces. Same nonlinear core moments and original full operators yield exact zero total leading core stress, nonzero full remainder and completed core/first trace. Axis/global/cone/independent time remainder/required-domain energy/resolved points/n-dependent recursion remain open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current positive-radius full core T/E and completed core/first trace PASS;33/32/10;axis/global/time open',flush=True)
    return result


if __name__=='__main__':run()
