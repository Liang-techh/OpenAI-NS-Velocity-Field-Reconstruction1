"""Independent admission of same-source nonsingular axis T/E and Cartesian C2."""
import copy
import json
from pathlib import Path
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_core_physical_field import coefficient_value
from lei_ren_part1_paper_compliant_current_core_axis_background_tensor import (
    CurrentCoreAxisBackgroundTensor,VIEWS,GATES,OPEN,NAME,RECEIPT,HERE,sha,pack,encode,
    endpoints,source_precision,_verify_hashes,read_producer)
from lei_ren_part1_paper_compliant_current_core_axis_operator import (
    INDICES2,INDICES3,indexkey,sector_metadata,raw_axis_source_theorem,cartesian_pullback_theorem,
    physical_factor_and_cylindrical_pullback_theorem)
from lei_ren_part1_paper_compliant_current_core_background_tensor_check import leaves,finite
from lei_ren_part1_paper_compliant_current_core_interior_moments_check import check_view as check_moments


def check_row(row):
    for key in ('signed_coefficient','physical_lambda_exponent','physical_viscosity_exponent','radial_log_prefactor'):finite(row[key])
    for value in row['actual_source_log_parts'].values():finite(value)
    if row['exact_zero']:
        if row['log_absolute_upper'] is not None or endpoints(row['signed_coefficient'])!=(mp.mpf(0),mp.mpf(0)):
            raise ValueError('Nonsingular source zero has bound')
    else:finite(row['log_absolute_upper'])
    if not row['positive_source_factors_not_materialized'] or not row['source_factors_combined_before_enclosure']:
        raise ValueError('Axis source factor materialized prematurely')
    return int(not row['exact_zero'])


def check_view(view,c,field):
    if view['chart']!='core_axis_Cartesian' or any(view[k] for k in OPEN) or not all(view[k] for k in (
        'no_inverse_radius_or_angle_at_axis','all_six_original_remainder_sectors_and_nonzero_derivatives_retained',
        'pressure_FTC_and_centrifugal_cancellation_from_same_C_source',
        'original_signed_stress_ledger_retained_in_checked_positive_radius_parent',
        'only_total_stress_analytic_extension_zero','source_function_equality_not_interval_overlap',
        'source_bounds_not_resolved_physical_points')):
        raise ValueError('Whole-core axis source/scope lost')
    if endpoints(view['rho'])[0]<0 or endpoints(view['rho'])[1]>4:raise ValueError('Foreign axis/source radius admitted')
    source=view['actual_same_nonlinear_core_source'];check_moments(source['same_original_core_interior_moment_packet'],c)
    if not source['original_F0_dressing_once_before_every_derivative'] or not source['no_physical_inverse_radius_used'] or not source['source_bounds_not_point_values']:
        raise ValueError('Same source and nonsingular dressing provenance lost')
    if set(source['ordinary_mixed4_source_jets'])!={'Q','V','F'}:raise ValueError('Original source components lost')
    expected_jets={'rho'+str(i)+'_Z'+str(k) for i in range(5) for k in range(5-i)}
    for rows in source['ordinary_mixed4_source_jets'].values():
        if set(rows)!=expected_jets:raise ValueError('Complete original source mixed4 jets required')
        for value in rows.values():finite(value)
    stress=view['physical_core_Cartesian_stress_mixed3'];div=view['physical_core_Cartesian_stress_divergence_mixed2']
    error=view['physical_core_Cartesian_remainder_mixed2'];meta=sector_metadata()
    if set(stress)!={indexkey(*ik) for ik in INDICES3} or set(div)!=set(error) or set(error)!={indexkey(*ik) for ik in INDICES2}:
        raise ValueError('Full Cartesian stress3/divergence2/remainder2 inventory lost')
    n=nonzero=0
    for ik in INDICES3:
        parts=stress[indexkey(*ik)]
        if set(parts)!={'xx','xy','xz','yy','yz','zz'}:raise ValueError('Completed tensor component lost')
        for row in parts.values():
            nonzero+=check_row(row);n+=1
            if not row['exact_zero'] or endpoints(row['physical_viscosity_exponent'])!=endpoints(c.mpf(1)-c.mpf(sum(ik))/2):
                raise ValueError('Exact total stress extension/unit lost')
    for ik in INDICES2:
        parts=div[indexkey(*ik)]
        if set(parts)!={'x','y','z'}:raise ValueError('Cartesian divergence component lost')
        for row in parts.values():
            nonzero+=check_row(row);n+=1
            if not row['exact_zero'] or endpoints(row['physical_viscosity_exponent'])!=endpoints(c.mpf('.5')-c.mpf(sum(ik))/2):
                raise ValueError('Exact total stress divergence extension/unit lost')
        parts=error[indexkey(*ik)]
        if set(parts)!={'x','y','z'}:raise ValueError('Cartesian remainder vector lost')
        for component,contributions in parts.items():
            expected={name for name,data in meta.items() if component in data['seed']}
            if set(contributions)!=expected:raise ValueError('Original six remainder sectors lost')
            for name,row in contributions.items():
                nonzero+=check_row(row);n+=1
                if endpoints(row['physical_viscosity_exponent'])!=endpoints(c.mpf('.5')-c.mpf(sum(ik))/2):
                    raise ValueError('Original sqrt-nu remainder unit changed')
                factor=c.mpf(str(meta[name]['epsilon_power']))-c.mpf(ik[0]+ik[1])/2
                if endpoints(row['actual_source_log_parts']['original_epsilon_core_log'])!=endpoints(factor*field.log_epsilon):
                    raise ValueError('Original epsilon factor lost or counted twice')
                if endpoints(row['actual_source_log_parts']['original_F0_base_log'])!=endpoints(meta[name]['F0_power']*field.log_F0):
                    raise ValueError('Original F0 base factor lost or differentiated twice')
                beta=coefficient_value(c,s.sympify(meta[name]['beta']),view['X'],view['Y'],view['Z'],field.core.delta,0)
                if endpoints(row['physical_lambda_exponent'])!=endpoints(beta-ik[0]-ik[1]+ik[2]*(field.core.delta-1)):
                    raise ValueError('Original per-derivative lambda exponent lost')
                if endpoints(row['radial_log_prefactor'])!=endpoints(c.mpf(0)):
                    raise ValueError('Cancelled cylindrical radius factor reintroduced at axis')
    if n!=260:raise ValueError('Full core-axis physical row inventory differs')
    if encode(pack(view['physical_core_Cartesian_momentum_decomposition_mixed2']))!=encode(pack(error)):
        raise ValueError('Exact zero-divergence momentum decomposition differs')
    if view['exact_axis']:
        base=error[indexkey(0,0,0)]
        if any(not row['exact_zero'] for comp in ('x','y') for row in base[comp].values()):
            raise ValueError('Axis transverse vector parity failed')
        transverse=[row for ik in ((1,0,0),(0,1,0)) for parts in error[indexkey(*ik)].values() for row in parts.values()]
        if not any(not row['exact_zero'] for row in transverse):raise ValueError('Nonzero transverse derivatives removed on axis')
    return dict(rows=n,nonzero=nonzero,exact_axis=view['exact_axis'],original_source_mixed4_rows=45)


@source_precision
def run(field=None):
    raw=read_producer();_verify_hashes(raw)
    field=field if field is not None else CurrentCoreAxisBackgroundTensor(require_checked=False);field.assert_graph()
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k!='current_core_axis_and_Cartesian_views'}:
        raise ValueError('Whole core axis source manifest differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_core_axis_and_Cartesian_views'])!=set(VIEWS):
        raise ValueError('Producer core-axis scope differs')
    # Rebuild the symbolic proofs independently of the producer's cache.
    for theorem in (raw_axis_source_theorem,cartesian_pullback_theorem,physical_factor_and_cylindrical_pullback_theorem):theorem.cache_clear()
    unit=raw_axis_source_theorem();coordinate=cartesian_pullback_theorem()
    if encode(pack(unit))!=encode(pack(field.proof['original_six_sector_nonsingular_source_theorem'])) or encode(pack(coordinate))!=encode(pack(field.proof['independent_nonsingular_Cartesian_coordinate_pullback'])):
        raise ValueError('Independent axis source/pullback replay differs')
    factors=physical_factor_and_cylindrical_pullback_theorem()
    if encode(pack(factors))!=encode(pack(field.proof['original_positive_radius_to_axis_physical_C2_factor_and_coefficient_theorem'])):
        raise ValueError('Independent original cylindrical/axis physical factor replay differs')
    field.core_tensor.interior.cache={};field.core_tensor.interior.rows_cache={};field.core_tensor.interior.tail_cache={}
    counts={}
    for name,args in VIEWS.items():
        view=field.cartesian(*args)
        if encode(pack(view))!=raw['current_core_axis_and_Cartesian_views'][name]:raise ValueError('Full core-axis source replay differs: '+name)
        counts[name]=check_view(view,field.ctx,field)
        if name=='axis_center':
            axial=view['physical_core_Cartesian_remainder_mixed2'][indexkey(0,0,0)]['z']['axial_axial_viscosity']
            if endpoints(axial['signed_coefficient'])[0]<=0:raise ValueError('Original nonzero axis axial remainder removed')
        print('Check nonsingular core Cartesian T/E: '+name,flush=True)
    rejected=[]
    for label,args in (('outside_core',('.2',3,0)),('outside_Z',(2,0,0)),('nonfinite_X',('.2','inf',0)),('nonfinite_Z',('inf',0,0))):
        try:field.cartesian(*args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Foreign core-axis domain admitted')
    for label,kwargs in (('nonfinite_time',dict(log_tau='inf')),('nonpositive_viscosity',dict(viscosity=0))):
        try:field.axis('.2',**kwargs)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid physical axis sector admitted')
    clone=copy.copy(field);clone.core=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_original_core')
    else:raise ValueError('Foreign original core accepted')
    clone=copy.copy(field);clone.core_tensor=copy.copy(field.core_tensor);clone.core_tensor.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_positive_radius_core')
    else:raise ValueError('Unchecked positive-radius core prerequisite accepted')
    clone=copy.copy(field);clone.remainder=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_axis_remainder_operator')
    else:raise ValueError('Foreign nonsingular remainder operator accepted')
    try:field.core_tensor.chart('core','.2',0)
    except ValueError:rejected.append('parent_cylindrical_axis_still_rejected')
    else:raise ValueError('Axis routed through singular cylindrical parent')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        current_core_axis_Cartesian_view_counts=counts,current_core_axis_Cartesian_view_count=len(VIEWS),
        current_core_axis_physical_rows_checked=sum(v['rows'] for v in counts.values()),
        current_core_axis_source_and_nonsingular_limit_theorem=field.proof,
        actual_current_tensor_regions_available=raw['actual_current_tensor_regions_available'],
        actual_current_completed_tensor_adjacent_interface_count=32,actual_current_completed_tensor_internal_interface_count=10,
        current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),
        whole_axis_and_same_source_positive_radius_analytic_attachment_certified=True,
        exact_axis_parity_and_nonzero_transverse_derivatives_retained=True,
        original_nonzero_axis_axial_remainder_retained=True,
        invalid_domains_and_foreign_current_owners_rejected=rejected,
        scope='Same actual whole core rho[0,4], whole Z[-1,1], nonsingular Cartesian axis stress3/divergence2/remainder2 with all6original sectors, all45mixed4 source jets and full tails.33regions/32adjacent/10internal unchanged. Global/cone/independent time-flat remainder/required-domain energy/resolved points/actual n-dependent recursion remain open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current same-source nonsingular axis Cartesian T/E PASS; whole core axis/C2 admitted;33/32/10 unchanged;global/time/recursion open',flush=True)
    return result


if __name__=='__main__':run()
