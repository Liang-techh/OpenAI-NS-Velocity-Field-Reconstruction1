"""Focused continuous whole-gap cone and independent original-evaluator checks."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_gap_cone import (
    PREFIX,DOMAIN,TAIL_DOMAIN,ADMISSIONS,FALSE_FLAGS,current_sources,
    gap_cone_identities,whole_gap_bounds,relative_log_envelopes,source_precision)
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def independent_zero_shear_gap_fixture():
    """Moderate unchanged gap histories, original full stress and unit oracle.

    The actual cone uses continuous source enclosures. These fixtures test
    formulas and units; their parameters are not the actual selected source.
    """
    with mp.workdps(110):
        c=MPIntervalContext();c.dps=130
        mu=mp.mpf('.1');delta=mp.mpf('.06');r=1-mu;p=1+2*mu;bp=mp.mpf('.5')+mu
        Bv=mp.exp(-40);D0=mp.exp(-40);Rv=mp.exp(100);Xp=mp.mpf('.3')
        lp=mp.log(Bv)+13/(2*mu)+13;lu=mp.mpf(0)
        lrp=mp.log(Rv)-13/mu;finite=mp.log(D0)+1/mu
        envelopes=relative_log_envelopes(c,c.mpf(mu),c.mpf(finite),c.mpf(lp),c.mpf(lu),c.mpf(lrp))
        count=0;worst=mp.mpf(0);log_checks=0;nonzero=0;tol=mp.mpf('1e-60')
        def compare(actual,expected):
            nonlocal count,worst
            error=abs(actual-expected);worst=max(worst,error);count+=1
            if error>tol*max(1,abs(actual),abs(expected)):
                raise ArithmeticError('Independent gap normalization differs')
        C=lambda z:1/(1+z*z)
        M1=lambda z:-mp.mpf('.12')*(1+z)
        M2=lambda z:mp.mpf('.07')*(1-z*z)
        future=lambda z:mp.mpf('.8')+mp.mpf('.1')*z*z
        loss=lambda z:mp.mpf('.1')+mp.mpf('.05')*z*z
        P0=lambda z:-mp.mpf('.35')*C(z)**2+mp.mpf('.02')*z
        for distance,z in ((4*mu,mp.mpf('.31')),(mp.mpf(1),mp.mpf('-.4')),(mp.mpf(2),mp.mpf('.5'))):
            R=Rv*mp.exp(-distance/mu);B=Bv*mp.exp(bp*distance/mu)
            D1=D0*mp.exp((mp.mpf('.5')-mu)*distance/mu)
            D2=D0*mp.exp((mp.mpf('.5')-2*mu)*distance/mu)
            H=mp.exp(-r*(13-distance)/mu);Qp=mp.exp(-p*distance/mu);K=mp.exp(-2*distance)
            def data(zz):
                cc=C(zz);Ut=B*cc
                pressure=B*B*(-cc*cc/(2*p)+Qp*(P0(zz)+cc*cc/(2*p)))
                e=future(zz)*K-mp.expm1(-2*distance)/(4*mu)-D0*D0*loss(zz)*K
                moments=dict(theta=mp.sqrt(2)*R**mp.mpf('1.5')*Ut*(1/r+(Xp-1/r)*H),
                    z=R*Ut*D1*M1(zz),theta_z=mp.sqrt(2)*R**mp.mpf('1.5')*Ut**2*D2*M2(zz),
                    z_theta=R*Ut**2*e,p=mp.mpf(0))
                return Ut,pressure,moments
            Ut,P,moments=data(z)
            mz={name:mp.diff(lambda zz:data(zz)[2][name],z) for name in moments}
            actual=evaluate_mp_stress(mp.log(R),z,delta,Utheta=Ut,Uz=0,
                Utheta_y=-bp*Ut,Utheta_Z=B*mp.diff(C,z),Uz_y=0,Uz_Z=0,
                moments=moments,moments_Z=mz,P=P,P_Z=mp.diff(lambda zz:data(zz)[1],z),
                precision=mp.mp.dps,radius_override=R,axial_override=z)
            St,Sz=actual['S_theta'],actual['S_z'];F=Ut/mp.sqrt(2*R);A=-St
            if Sz!=0 or St>=0:raise ArithmeticError('Actual gap source shear differs')
            expectedUr=(-(1-delta)*z*moments['z']-(1-z*z)*mz['z'])/((1-delta*z*z)*mp.sqrt(2*R))
            if expectedUr==0 or moments['z']==0 or moments['theta_z']==0:
                raise ArithmeticError('Fixture lost radial/cumulative history')
            compare(actual['U_r']/expectedUr,1);nonzero+=1
            kappa=-(St*St+Sz*Sz)/(F*St);compare(kappa-2,2*mu)
            Qbase=mp.sqrt(R/2)*B
            theta=actual['T_theta']/Qbase;axial=actual['T_z']/Qbase
            dot=actual['T_theta']*St+actual['T_z']*Sz
            cross=-actual['T_theta']*Sz+actual['T_z']*St
            normalized=2*theta**2-2*mu*axial**2
            compare(-dot/(Qbase*A),theta)
            compare(-cross/(Qbase*A),axial)
            compare((2*dot**2-(kappa-2)*cross**2)/(Qbase*A)**2,normalized)
            if theta<=0 or dot>=0 or normalized<=0:
                raise ArithmeticError('Original moderate gap fixture cone failed')
            for nu,lam in ((mp.mpf('.01'),mp.mpf('.8')),(mp.mpf('.7'),mp.mpf('1.2'))):
                fac=lam**(-2-delta);common=nu*fac
                Stp,Szp=common*St,common*Sz;Fp=fac*F
                kp=-(Stp*Stp+Szp*Szp)/(nu*Fp*Stp)
                dp=common**2*dot;cp=common**2*cross
                compare(kp,kappa)
                compare(-dp/(common**2*Qbase*A),theta)
                compare((2*dp**2-(kp-2)*cp**2)/(common**4*(Qbase*A)**2),normalized)
            direct=dict(theta_signed_original_memory=mp.log(H),
                theta_meridional_Mz_history=mp.log(B*D1),
                theta_mixed_angular_axial_history=mp.log(B*D2),theta_radial_shear=-mp.log(R),
                axial_full_unperturbed_energy_and_pressure=mp.log(B),
                axial_same_absolute_pressure_memory=mp.log(B*Qp),
                axial_linear_axial_moment_history=mp.log(D1),
                axial_selected_backward_energy_loss=mp.log(B*D0**2))
            for key,value in direct.items():
                lo,hi=endpoints(envelopes[key]);log_checks+=1
                if value-hi>tol:raise ArithmeticError('Whole-gap log envelope fails moderate source: '+key)
                endpoint=4*mu if key=='axial_same_absolute_pressure_memory' else mp.mpf(2)
                if distance==endpoint:compare((lo+hi)/2,value)
        return dict(comparisons=count,log_envelope_comparisons=log_checks,
            nonzero_radial_and_cumulative_history_fixtures=nonzero,
            zero_axial_shear_from_original_evaluator=True,viscosities=['.01','.7'],
            tolerance=str(tol),maximum_absolute_error=mp.nstr(worst,30),
            small_margins_normalized_before_comparison=True,
            actual_source_cone_proof_does_not_use_fixtures=True)


@source_precision
def run():
    records,hashes,family=current_sources()
    name=PREFIX+'pulse_gap_cone.json';raw=(HERE/name).read_bytes();record=json.loads(raw)
    for path,digest in record['input_hashes'].items():
        if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
            raise ValueError('Current gap cone input changed: '+path)
    if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=family:
        raise ValueError('Cone current source family differs')
    proof=gap_cone_identities(records);bounds=whole_gap_bounds(records)
    if record['exact_source_cone_proof']!=encode(pack(proof)) or record['bounds']!=encode(pack(bounds)):
        raise ValueError('Current continuous gap proof/bounds changed')
    for group in ('positive_parameter_margins','positive_log_margins','positive_algebraic_margins'):
        for key,value in bounds[group].items():
            lo,hi=endpoints(value)
            if lo<=0 or not mp.isfinite(hi):raise ArithmeticError('Gap strict continuous margin failed: '+key)
    if endpoints(bounds['exact_axial_shear_ratio'])!=(0,0):raise ArithmeticError('Gap source axial shear is not zero')
    for flag in ('continuous_whole_original_gap_Z_domain_covered',
        'exact_source_correlations_grouped_before_interval_enclosure',
        'all_four_angular_and_four_axial_histories_retained',
        'zero_axial_shear_proved_from_actual_gap_source',
        'radial_velocity_and_five_moment_histories_not_zeroed',
        'main_pulse_axial_shear_not_assumed_zero','global_temporal_flatness_not_inferred'):
        if not bounds[flag]:raise ValueError('Gap cone scope missing: '+flag)
    if bounds['phase_samples_used_as_proof'] or bounds['source_caps_used_as_defining_field_values']:
        raise ValueError('Samples or caps selected as cone field')
    if record['domain']!=DOMAIN or record['tail_domain']!=TAIL_DOMAIN:
        raise ValueError('Whole-gap/composed tail domain differs')
    for flag in ADMISSIONS:
        if not record[flag]:raise ValueError('Whole-gap regional admission missing: '+flag)
    for flag in FALSE_FLAGS:
        if record[flag]:raise ValueError('Gap regional scope overclaimed: '+flag)
    if not record['current_selected_C5_source_check_directly_consumed']:
        raise ValueError('Actual selected C5 source gate missing')
    fixture=independent_zero_shear_gap_fixture()
    hashes.update(record['input_hashes']);hashes[name]=hashlib.sha256(raw).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(all_passed=True,actual_five_defect_family_sha256=family[0],implicit_source_sha256=family[1],
        domain=DOMAIN,tail_domain=TAIL_DOMAIN,input_hashes=hashes,
        exact_source_cone_identities_checked=len(proof['identities']),
        positive_parameter_margins_checked=len(bounds['positive_parameter_margins']),
        positive_log_source_margins_checked=len(bounds['positive_log_margins']),
        positive_algebraic_cone_margins_checked=len(bounds['positive_algebraic_margins']),
        actual_continuous_whole_gap_source_bounds_recomputed=True,
        current_selected_C5_source_check_directly_consumed=True,
        same_source_completed_physical_gap_end_join_and_downstream_tail_consumed=True,
        original_gap_zero_axial_shear_with_nonzero_radial_histories_verified=True,
        exact_source_correlations_grouped_before_interval_enclosure=True,
        independent_zero_shear_gap_fixture=fixture,
        phase_samples_used_as_proof=False,source_caps_used_as_defining_field_values=False,
        **{flag:True for flag in ADMISSIONS},**{flag:False for flag in FALSE_FLAGS})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('PASS continuous whole original gap cone and composed end/tail; main/entrance/global/recursion pending',flush=True)
    return result


if __name__=='__main__':run()
