"""Focused full-axial-shear pulse-end cone and physical-unit check."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_pulse_end_cone import (
    PREFIX,DOMAIN,TAIL_DOMAIN,ADMISSIONS,FALSE_FLAGS,current_sources,
    pulse_end_cone_identities,whole_pulse_end_bounds,source_precision)
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def independent_full_shear_fixture():
    """Original shear evaluator, both axial signs and fixed-viscosity pullbacks.

    Moderate fixtures check normalization and units. The actual cone proof
    consumes continuous source boxes, never these point fixtures.
    """
    with mp.workdps(90):
        worst=mp.mpf(0); count=0; nonzero=0; tol=mp.mpf('1e-60')
        def compare(actual,expected):
            nonlocal worst,count
            error=abs(actual-expected);worst=max(worst,error);count+=1
            if error>tol*max(1,abs(actual),abs(expected)):
                raise ArithmeticError('Independent full-shear normalization differs')
        for Z,Bhs in (('.31','.9'),('.31','-.9'),('-.4','.7'),('-.4','-.7')):
            z=mp.mpf(Z);mu=mp.mpf('.1');delta=mp.mpf('.06')
            R=mp.mpf(7);B=mp.mpf('.6');D=mp.mpf('.4');Bh=mp.mpf('.3')
            C=1/(1+z*z);Ut=B*C;Uz=Ut*D*Bh
            Uy=-(mp.mpf('.5')+mu)*Ut
            Uzy=Ut*D*(mp.mpf(Bhs)-(mp.mpf('.5')+mu)*Bh)
            zero={key:mp.mpf(0) for key in ('theta','z','theta_z','z_theta','p')}
            actual=evaluate_mp_stress(mp.log(R),z,delta,Utheta=Ut,Uz=Uz,Utheta_y=Uy,
                Utheta_Z=0,Uz_y=Uzy,Uz_Z=0,moments=zero,moments_Z=zero,P=0,P_Z=0,
                precision=90,radius_override=R,axial_override=z)
            St=actual['S_theta'];Sz=actual['S_z'];F=Ut/mp.sqrt(2*R);A=-St
            sigma=Sz/A;a=2+2*mu;Q=mp.sqrt(R/2)*B
            kappa=-(St*St+Sz*Sz)/(F*St)
            compare(sigma,D*(mp.mpf(Bhs)-(mp.mpf('.5')+mu)*Bh)/(1+mu))
            compare(kappa-2,2*mu+a*sigma*sigma)
            if Sz==0 or kappa<=a:raise ArithmeticError('Fixture omitted axial shear')
            nonzero+=1
            theta=mp.mpf('.8');axial=mp.mpf('.02')*(-1 if z<0 else 1)
            Tt=Q*theta;Tz=Q*axial
            dot=Tt*St+Tz*Sz;cross=-Tt*Sz+Tz*St
            compare(-dot/(Q*A),theta-axial*sigma)
            compare(-cross/(Q*A),theta*sigma+axial)
            margin=2*dot**2-(kappa-2)*cross**2
            compare(margin/(Q*A)**2,2*(theta-axial*sigma)**2-(kappa-2)*(theta*sigma+axial)**2)
            if dot>=0 or margin<=0:raise ArithmeticError('Full-shear fixture cone failed')
            for nu,lam in (('.01','.8'),('.7','1.2')):
                nu=mp.mpf(nu);lam=mp.mpf(lam);fac=lam**(-2-delta);common=nu*fac
                Stp=common*St;Szp=common*Sz;Fp=fac*F
                kappap=-(Stp*Stp+Szp*Szp)/(nu*Fp*Stp)
                dotp=(common*Tt)*Stp+(common*Tz)*Szp
                crossp=-(common*Tt)*Szp+(common*Tz)*Stp
                compare(kappap,kappa)
                compare(2*dotp**2-(kappap-2)*crossp**2,common**4*margin)
        return dict(comparisons=count,nonzero_axial_shear_fixtures=nonzero,
            both_axial_shear_signs_exercised=True,viscosities=['.01','.7'],
            tolerance=str(tol),maximum_absolute_error=mp.nstr(worst,30),
            actual_source_cone_proof_does_not_use_fixtures=True)


@source_precision
def run():
    records,hashes,family=current_sources()
    name=PREFIX+'pulse_end_cone.json';raw=(HERE/name).read_bytes();record=json.loads(raw)
    for source,digest in record['input_hashes'].items():
        if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:
            raise ValueError('Pulse-end cone source changed: '+source)
    if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=family:
        raise ValueError('Cone report current family differs')
    proof=pulse_end_cone_identities();bounds=whole_pulse_end_bounds(records)
    if record['exact_source_cone_proof']!=encode(pack(proof)):
        raise ValueError('Original source cone identities changed')
    if record['bounds']!=encode(pack(bounds)):
        raise ValueError('Current whole-source cone bounds changed')
    if not all(proof['identities'].values()):raise ValueError('Original cone identity missing')
    for group in ('positive_log_margins','positive_algebraic_margins'):
        for key,value in bounds[group].items():
            lo,hi=endpoints(value)
            if lo<=0 or not mp.isfinite(hi):raise ArithmeticError('Continuous cone margin failed: '+key)
    for key in ('continuous_whole_original_s_Z_domain_covered',
        'actual_source_log_recipes_consumed_without_giant_log_subtraction',
        'signed_original_memory_and_full_ten_sectors_retained','full_axial_shear_not_set_to_zero',
        'global_temporal_flatness_not_inferred'):
        if not bounds[key]:raise ValueError('Continuous original scope missing: '+key)
    if bounds['phase_samples_used_as_proof'] or bounds['source_caps_used_as_defining_field_values']:
        raise ValueError('Samples or caps substituted as original fields')
    if record['domain']!=DOMAIN or record['tail_domain']!=TAIL_DOMAIN:
        raise ValueError('Original pulse and composed tail domain changed')
    for flag in ADMISSIONS:
        if not record[flag]:raise ValueError('Regional cone admission missing: '+flag)
    for flag in FALSE_FLAGS:
        if record[flag]:raise ValueError('Regional pulse-end scope overclaimed: '+flag)
    if not record['current_selected_C5_source_check_directly_consumed']:
        raise ValueError('Direct selected C5 source gate required')
    fixture=independent_full_shear_fixture()
    hashes.update(record['input_hashes']);hashes[name]=hashlib.sha256(raw).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(all_passed=True,actual_five_defect_family_sha256=family[0],
        implicit_source_sha256=family[1],domain=DOMAIN,tail_domain=TAIL_DOMAIN,
        exact_source_cone_identities_checked=len(proof['identities']),
        positive_log_source_margins_checked=len(bounds['positive_log_margins']),
        positive_algebraic_cone_margins_checked=len(bounds['positive_algebraic_margins']),
        actual_continuous_whole_end_source_bounds_recomputed=True,
        original_full_axial_shear_preserved_without_kappa_minus2_subtraction=True,
        actual_signed_memory_meridional_radial_and_axial_sectors_bounded=True,
        current_selected_C5_source_check_directly_consumed=True,
        same_source_physical_and_support_interfaces_consumed=True,
        normalized_support_bounds_restored_via_parent_modes_and_logs=True,
        phase_samples_used_as_proof=False,source_caps_used_as_defining_field_values=False,
        independent_full_shear_fixture=fixture,input_hashes=hashes,
        **{flag:True for flag in ADMISSIONS},**{flag:False for flag in FALSE_FLAGS})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('PASS continuous whole pulse-end full-shear cone and joined tail; main/global/recursion pending',flush=True)
    print('Cone identities',len(proof['identities']),'positive margins',
        len(bounds['positive_log_margins'])+len(bounds['positive_algebraic_margins']),flush=True)
    return result


if __name__=='__main__':run()
