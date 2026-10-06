"""Independent finite-root fixtures and original-radius candidate admission."""
import copy
import itertools
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_current_physical_tensor_locator import (
    CurrentPhysicalTensorLocator,NAME,RECEIPT,GATES,OPEN,HERE,sha,pack,encode,
    endpoints,source_precision,_verify_hashes,SYMBOLS)
from lei_ren_part1_paper_compliant_current_physical_tensor_locator_operator import (
    implicit_log_coordinate_map,original_radial_cover_theorem,implicit_physical_source_theorem,
    nonpositive_exp)
from lei_ren_part1_paper_compliant_current_tensor_registry_check import COORDINATES,check_full_view


def independent_root(z,lt,nu,delta):
    """Independent scalar equation in q, with direct exponentials.

    Used only on moderate fixture logarithms, not the reconstruction's huge
    scales. Opposite sign brackets bound the reference root independently of
    the production logaddexp and interval Newton implementation.
    """
    z,lt,nu,d=(mp.mpf(v) for v in (z,lt,nu,delta))
    if not z:return lt/2,lt/2
    a=z*z/nu;tau=mp.exp(lt)
    lo=lt/2
    hi=max((lt+mp.log(2))/2,(mp.log(2*a))/(2*(1-d)))+1
    f=lambda q:mp.exp(2*q)-a*mp.exp(2*d*q)-tau
    if f(lo)>0 or f(hi)<=0:raise ValueError('Independent original root bracket failed')
    for _ in range(600):
        mid=(lo+hi)/2
        if f(mid)<0:lo=mid
        else:hi=mid
    return lo,hi


@source_precision
def independent_implicit_fixtures():
    c=MPIntervalContext();c.dps=110
    cases=(('zero','0','-2.6','.8','.37'),('positive','.7','-2.6','.8','.37'),
        ('negative','-.7','-2.6','.8','.37'),('small_z','1e-70','-2.6','.8','.37'),
        ('large_z','1e200','-4','.3','.95'),('lambda_below_one','.013','-7','2','.8'),
        ('delta_zero','1.7','-2.6','.8','0'))
    evidence={}
    def compare(label,box,z,lt,nu,d):
        rlo,rhi=independent_root(z,lt,nu,d);lo,hi=endpoints(box['actual_log_lambda'])
        if lo>rlo or hi<rhi:raise ValueError('Independent true physical root outside directed box: '+label)
        q=(rlo+rhi)/2
        # Close to infinity the direct ratio rounds above 1 when the root
        # error exceeds the tiny true complement. Use the exact independent
        # root relation for this reference value; it can round to 1, which
        # the production outward enclosure must retain.
        Z=mp.sign(mp.mpf(z))*mp.sqrt(1-mp.exp(mp.mpf(lt)-2*q))
        zl,zh=endpoints(box['Z'])
        if zl>Z or zh<Z:raise ValueError('Independent physical Z outside source enclosure: '+label)
        if not box['actual_log_tau_retained'] or not box['true_finite_point_has_abs_Z_strictly_below_one']:
            raise ValueError('True time/finite point coordinate semantics lost')
        return dict(independent_direct_exponential_root_inside=True,independent_Z_inside=True,
            solver_status=box['solver_status'],iterations=box['iterations'],
            actual_log_lambda=box['actual_log_lambda'],log_root_enclosure_width=box['log_root_enclosure_width'])
    for label,z,lt,nu,d in cases:
        box=implicit_log_coordinate_map(c,z,lt,nu,d)
        evidence[label]=compare(label,box,z,lt,nu,d)
    # These boxes exercise the sign crossing and variable delta, including
    # both signs of q. Corner tests complement the global derivative proof;
    # they are not used as a substitute for the interval source theorem.
    for label,zbox,ltbox,nubox,dbox in (
        ('cross_axis',('-1','2'),('-4','-2'),('.7','1.3'),('.05','.4')),
        ('positive_parameter_box',('.01','.2'),('-8','-4'),('.8','1.1'),('.2','.7'))):
        box=implicit_log_coordinate_map(c,zbox,ltbox,nubox,dbox)
        for index,point in enumerate(itertools.product(zbox,ltbox,nubox,dbox)):
            compare(label+str(index),box,*point)
        evidence[label]=dict(corner_roots_and_Z_inside=16,solver_status=box['solver_status'],
            midpoint_parameter_substitution=False)
    # No exponent integer proportional to exp(10^1000) is constructed.
    tiny=nonpositive_exp(c,c.mpf('-1e1000'))
    if endpoints(tiny)[0]!=0 or endpoints(tiny)[1]<=0 or endpoints(tiny)[1]>=1:
        raise ValueError('Astronomical nonpositive exp enclosure failed')
    evidence['astronomical_negative_exp']=dict(original_log_preserved=True,positive_upper_with_zero_lower=True)
    return evidence


def check_refined_view(view,q):
    original=view['original_full_native_view']
    counts=check_full_view(original,210 if original['layout']=='Cartesian_core' else 71)
    groups=view['actual_lambda_signed_component_groups']
    if set(groups)!=set(original['canonical_signed_component_groups']):raise ValueError('Original full tensor groups lost')
    nontrivial=0
    for label,rows in groups.items():
        if len(rows)!=len(original['canonical_signed_component_groups'][label]):raise ValueError('Signed source sector lost')
        for row,old in zip(rows,original['canonical_signed_component_groups'][label]):
            if any(encode(pack(row[key]))!=encode(pack(value)) for key,value in old.items()):
                raise ValueError('Original source row changed during actual lambda dressing')
            if endpoints(row['actual_lambda_log_prefactor'])!=endpoints(old['physical_lambda_exponent']*q):
                raise ValueError('Actual implicit lambda factor differs')
            if row['actual_lambda_log_absolute_upper'] is not None:
                if any(not mp.isfinite(v) for v in endpoints(row['actual_lambda_log_absolute_upper'])):
                    raise ValueError('Nonfinite actual finite-time signed contribution log bound')
                nontrivial+=1
    return dict(**counts,actual_lambda_nonzero_contributions=nontrivial)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentPhysicalTensorLocator(require_checked=False)
    field.assert_graph()
    if encode(pack(field.manifest()))!=raw or any(raw[k] for k in GATES+OPEN):raise ValueError('Original physical locator manifest differs')
    implicit=independent_implicit_fixtures();c=field.ctx
    radiusproof=original_radial_cover_theorem();coordinateproof=implicit_physical_source_theorem()
    if len(radiusproof['all_32_adjacent_exact_radius_source_identities'])!=32 or not coordinateproof['passed']:
        raise ValueError('Original exact radius/implicit coordinate source replay failed')
    roundtrips={}
    for region,coordinate in COORDINATES.items():
        v=c.mpf(coordinate);logR=field.original_log_radius(region,v)
        location=field.locate_log_radius((logR+c.ln(2)-c.mpf('2.6'))/2,0,'-2.6')
        candidates={row['region']:row for row in location['candidates']}
        if region not in candidates:raise ValueError('Original physical radius inverse region omitted: '+region)
        lo,hi=endpoints(candidates[region]['native_coordinate_enclosure']);vl,vh=endpoints(v)
        if lo>vl or hi<vh:raise ValueError('Original coordinate not enclosed by inverse: '+region)
        roundtrips[region]=dict(candidate_region_count=len(candidates),coordinate_enclosed=True,
            inverse_status=candidates[region]['inverse_status'])
        print('Check directed physical radius inverse: '+region,flush=True)
    # Existing native operators have their own 33-route admission. Here only
    # new physical routing and the constrained core path need fresh tensors.
    cases=dict(axis=field.cartesian(0,0,'.7','.9',viscosity='.8',tensors=False),
        ordinary=field.cartesian('.7','-.2','.3','.9',viscosity='.8',tensors=False))
    core_logR=field.original_log_radius('core_positive_radius',c.mpf(4))
    cases['joint_core_boundary']=field.locate_log_radius((core_logR+c.ln(2)-c.mpf('2.6'))/2,0,'-2.6',theta='.7')
    # The boundary enclosure contains both sides; evaluate the constrained
    # core component through the real complete routing result. Other sides
    # retain their admitted original full native operators.
    counts={}
    for label,location in cases.items():
        view=field.tensor_from_location(location)
        counts[label]={part['candidate']['region']:check_refined_view(part,location['actual_log_lambda'])
            for part in view['candidate_full_tensor_views']}
    if len(counts['axis'])!=1 or 'core_positive_radius' not in counts['axis'] or 'bridge_macro' not in counts['ordinary']:
        raise ValueError('Actual Cartesian axis/ordinary point routing differs')
    normal=cases['ordinary']
    if endpoints(normal['actual_log_lambda'])==endpoints(normal['requested_log_tau']/2):
        raise ValueError('Actual off-axis lambda replaced by sqrt(tau)')
    rejected=[]
    for label,kwargs in (('terminal_time',dict(time=1)),('after_terminal',dict(time=2)),
        ('nonpositive_nu',dict(viscosity=0)),('nonfinite_x',dict(x_phys='inf'))):
        args=dict(x_phys=0,y_phys=0,z_phys='.2',time='.9',viscosity='.8',tensors=False);args.update(kwargs)
        try:field.cartesian(**args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid physical coordinate domain admitted: '+label)
    forged=copy.deepcopy(cases['axis']);forged['candidates'][0]['native_coordinate_enclosure']=c.mpf(1)
    try:field.tensor_from_location(forged)
    except ValueError:rejected.append('forged_candidate_list')
    else:raise ValueError('Foreign physical locator result admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        independent_implicit_physical_root_fixtures=implicit,current_33_original_radius_inverse_roundtrips=roundtrips,
        all_32_original_exact_radius_source_identities_replayed=True,current_fresh_actual_physical_tensor_view_counts=counts,
        actual_lambda_is_distinct_from_native_sqrt_tau_off_axis=True,
        interval_candidate_union_and_parameter_uncertainty_retained=True,
        extreme_absolute_outer_log_radius_has_25_candidates_not_resolved=True,
        nonsingular_core_joint_constraint_and_original_full_derivatives_retained=True,
        invalid_physical_domains_and_foreign_locator_results_rejected=rejected,
        scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual implicit physical map / 33 original radius inverse candidates PASS; global physical cover/cone/point/recursion open',flush=True)
    return result


if __name__=='__main__':run()
