"""Independent original-beta whole-gap similarity companion checks."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_gap_similarity_C4 import (
    CompliantPulseGapSimilarityC4,gap_shapes,DOMAIN,FALSE_FLAGS,PREFIX,source_precision)
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def independent_gap_fixture():
    """Direct full-beta integrals and original full stress; no cone claim."""
    with mp.workdps(110):
        c=MPIntervalContext();c.dps=140
        mu,delta,xp,Bv,D0,Rv,ell=map(mp.mpf,('.2','.12','.3','.45','.07','2.8','.15'))
        p=1+2*mu;r=1-mu;bp=mp.mpf('.5')+mu
        Pstar=Bv*mp.exp(13/(2*mu)+13);Pin=mp.mpf('.12')
        normal=mp.quad(lambda t:mp.exp(-1/(1-t*t)),[-1,0,1])
        def beta(t):
            x=t/ell
            return mp.exp(-1/(1-x*x))/(ell*normal) if abs(x)<1 else mp.mpf(0)
        cuts=[-ell,0,ell]
        W=[mp.quad(lambda t,lam=mp.mpf('.5')-i*mu:mp.exp(lam*t)*beta(t),cuts) for i in (1,2)]
        gram=mp.quad(lambda t:mp.exp(-2*mu*t)*beta(t)**2,cuts)
        coeff=(lambda z:mp.mpf('.1')*(1+z),lambda z:-mp.mpf('.07')*(1-z*z))
        C=lambda z:1/(1+z*z)
        F=lambda z:mp.mpf('.8')+mp.mpf('.1')*z*z
        P0=lambda z:-mp.mpf('.35')*C(z)**2+mp.mpf('.02')*z
        M=lambda z,i:-sum(coeff[j](z)*mp.exp((mp.mpf('.5')-i*mu)*center)*W[i-1]
            for j,center in enumerate((-3,-1)))
        J0=lambda z:sum(coeff[j](z)**2*mp.exp(-2*mu*center)*gram for j,center in enumerate((-3,-1)))
        count={};worst=mp.mpf(0);tol=mp.mpf('1e-55');nonzero=0
        def compare(label,bound,value):
            nonlocal worst
            lo,hi=endpoints(bound);miss=max(lo-value,value-hi,mp.mpf(0));worst=max(worst,miss)
            if miss>tol*max(1,abs(value)):raise ArithmeticError('Independent gap fixture differs: '+label)
            count[label]=count.get(label,0)+1
        def jet(fn,z):return IntervalTaylor(c,[c.mpf(v) for v in mp.taylor(fn,z,5)])
        def reference(sv,z):
            R=Rv*mp.exp(sv);B=Bv*mp.exp(-bp*sv)
            H=mp.exp(-13*r/mu-r*sv);cc=C(z)
            m1=D0*M(z,1)*mp.exp(-(mp.mpf('.5')-mu)*sv)
            m2=D0*M(z,2)*mp.exp(-(mp.mpf('.5')-2*mu)*sv)
            K=mp.exp(2*mu*sv)
            e0=F(z)*K-mp.expm1(2*mu*sv)/(4*mu);loss=J0(z)*K
            pressure=-cc*cc/(2*p)+mp.exp(p*sv)*(P0(z)+cc*cc/(2*p))
            Ut=B*cc
            radial=( -(1-delta)*z*cc*m1-(1-z*z)*mp.diff(lambda v:C(v)*D0*M(v,1)*mp.exp(-(mp.mpf('.5')-mu)*sv),z)
                )/(1-delta*z*z)
            Ur=mp.sqrt(R/2)*B*radial
            X=1/r+(xp-1/r)*H
            moments=dict(theta=mp.sqrt(2)*R**mp.mpf('1.5')*Ut*X,
                z=R*Ut*m1,theta_z=mp.sqrt(2)*R**mp.mpf('1.5')*Ut*Ut*m2,
                z_theta=R*Ut*Ut*(e0-D0*D0*loss),
                p=Pstar*Pstar*cc*cc*(Pin+(1-mp.exp(-p*(13/mu+sv)))/(2*p)))
            return dict(R=R,B=B,H=H,Ur=Ur,Ut=Ut,Uz=mp.mpf(0),P=B*B*pressure,moments=moments)
        def direct_stress(sv,z):
            data=reference(sv,z);R=data['R'];B=data['B'];Ut=data['Ut']
            mz={name:mp.diff(lambda v:reference(sv,v)['moments'][name],z) for name in data['moments']}
            return evaluate_mp_stress(mp.log(R),z,delta,Utheta=Ut,Uz=0,Utheta_y=-bp*Ut,
                Utheta_Z=B*mp.diff(C,z),Uz_y=0,Uz_Z=0,moments=data['moments'],moments_Z=mz,
                P=data['P'],P_Z=mp.diff(lambda v:reference(sv,v)['P'],z),
                precision=mp.mp.dps,radius_override=R,axial_override=z)
        for distance,z in ((mp.mpf('.8'),mp.mpf('.31')),(mp.mpf(1),mp.mpf('-.4')),(mp.mpf(2),mp.mpf('.5'))):
            sv=-distance/mu;data=reference(sv,z);R=data['R'];B=data['B'];H=data['H']
            zj=IntervalTaylor.variable(c,c.mpf(z),5)
            cc=jet(C,z);rows=gap_shapes(c,c.mpf(mu),c.mpf(delta),zj,cc,c.mpf(xp),
                jet(lambda v:M(v,1),z),jet(lambda v:M(v,2),z),jet(F,z),jet(J0,z),jet(P0,z),c.mpf(distance))
            fixture=object.__new__(CompliantPulseGapSimilarityC4)
            fixture.ctx=c;fixture.mu=c.mpf(mu);fixture.delta=c.mpf(delta);fixture.Xp=c.mpf(xp)
            fixture.M=[jet(lambda v:M(v,i),z) for i in (1,2)]
            fixture.future=jet(F,z);fixture.J0=jet(J0,z);fixture.P0=jet(P0,z)
            fixture.Pin=c.mpf(Pin);fixture.U=c.mpf(1)
            fixture.logP=c.mpf(mp.log(Pstar));fixture.logU=c.mpf(0)
            fixture.logRp=c.mpf(mp.log(Rv)-13/mu);fixture.finite=c.mpf(mp.log(D0)+1/mu)
            source_packet=fixture.packet(c.mpf(z),c.mpf(distance),right_endpoint=(distance==4*mu))
            scales={'D0':D0,'D1':D0*mp.exp((mp.mpf('.5')-mu)*distance/mu),
                'D2':D0*mp.exp((mp.mpf('.5')-2*mu)*distance/mu)}
            Q=mp.exp(-p*distance/mu)
            for label,key in (('theta','T_theta'),('axial','T_z')):
                for j in range(4):
                    for n in range(4-j):
                        bound=c.mpf(0)
                        for part in rows['stress'][label].values():
                            rp,bpw,dp,hp=part['mode']
                            factor=R**rp*B**bpw*scales[part['selected_D_recipe']]**dp*H**hp/mp.sqrt(2)
                            if part['pressure_memory']:factor*=Q
                            bound+=part['full_derivative_rows'][j][n]*math.factorial(n)*c.mpf(factor)
                        expected=mp.diff(lambda ss:mp.diff(lambda zz:direct_stress(ss,zz)[key],z,n),sv,j)
                        compare('stress',bound,expected)
            modes={'radial':(mp.sqrt(R/2)*B*scales['D1'],'Ur'),'theta':(B,'Ut'),'axial':(B*scales['D1'],'Uz')}
            for label,(factor,key) in modes.items():
                for j in range(5):
                    for n in range(5-j):
                        bound=rows['velocity'][label][j][n]*math.factorial(n)*c.mpf(factor)
                        expected=mp.diff(lambda ss:mp.diff(lambda zz:reference(ss,zz)[key],z,n),sv,j)
                        compare('velocity',bound,expected)
            if data['Ur']!=0:nonzero+=1
            for j in range(5):
                for n in range(5-j):
                    bound=(rows['pressure_baseline_rows'][j][n]+rows['pressure_memory_rows'][j][n]*c.mpf(Q))*math.factorial(n)*c.mpf(B*B)
                    expected=mp.diff(lambda ss:mp.diff(lambda zz:reference(ss,zz)['P'],z,n),sv,j)
                    compare('pressure',bound,expected)
            # Raw cumulative linear moments must keep their nonzero history.
            for label in ('z','theta_z'):
                if data['moments'][label]==0:raise ArithmeticError('Fixture lost cumulative gap moment')
                for j in range(1,5):
                    expected=mp.diff(lambda ss:reference(ss,z)['moments'][label],sv,j)
                    compare('constant_raw_moment',c.mpf(0),expected)
            for label,sectors in source_packet['five_raw_cumulative_moment_log_sectors'].items():
                for j in range(5):
                    for n in range(5-j):
                        bound=c.mpf(0);rowkey='y'+str(j)+'_Z'+str(n)
                        for sector in sectors.values():
                            factor=c.exp(sum(sector['exact_source_log_parts'].values(),c.mpf(0)))
                            bound+=sector['full_moment_mixed4_coefficient_enclosures'][rowkey]*factor
                        expected=mp.diff(lambda ss:mp.diff(lambda zz:reference(ss,zz)['moments'][label],z,n),sv,j)
                        compare('all_five_raw_moments',bound,expected)
        return dict(comparisons=count,nonzero_radial_histories=nonzero,
            direct_complete_beta_integrals_used=True,original_full_stress_evaluator_used=True,
            ordinary_stress_order3_velocity_pressure_order4=True,tolerance=str(tol),
            maximum_positive_enclosure_miss=mp.nstr(worst,30),
            source_cone_and_corrected_residual_not_tested=True)


@source_precision
def run():
    companion=CompliantPulseGapSimilarityC4()
    actual=companion.report();name=PREFIX+'pulse_gap_similarity_C4.json'
    raw=(HERE/name).read_bytes();record=json.loads(raw)
    for path,digest in record['input_hashes'].items():
        if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
            raise ValueError('Gap similarity source changed: '+path)
    if record!=encode(pack(actual)):raise ValueError('Current whole gap source report differs')
    if record['domain']!=DOMAIN:raise ValueError('Whole original inactive gap domain changed')
    identities=companion.proof['identities']
    if not all(identities.values()):raise ValueError('Source-functional gap identity failed')
    counts={};zeros={}
    for key in ('whole_original_gap','right_end_source_packet','same_chart_switch_packet'):
        packet=actual[key]
        for group,rowkey in (
            ('full_meridional_stress_log_sectors','full_stress_mixed3_coefficient_enclosures'),
            ('full_absolute_pressure_log_sectors','full_pressure_mixed4_coefficient_enclosures'),
            ('full_velocity_log_sectors','full_velocity_mixed4_coefficient_enclosures')):
            for label,value in packet[group].items():
                sectors=value if group=='full_meridional_stress_log_sectors' else {label:value}
                for sector in sectors.values():
                    for row in sector[rowkey].values():
                        if not all(mp.isfinite(v) for v in endpoints(row)):raise ArithmeticError('Nonfinite gap row')
                        counts[group]=counts.get(group,0)+1
                        if endpoints(row)==(0,0):zeros[group]=zeros.get(group,0)+1
        for sectors in packet['five_raw_cumulative_moment_log_sectors'].values():
            for sector in sectors.values():
                for value in sector['full_moment_mixed4_coefficient_enclosures'].values():
                    if not all(mp.isfinite(v) for v in endpoints(value)):raise ArithmeticError('Nonfinite raw moment row')
                    counts['five_raw_moments']=counts.get('five_raw_moments',0)+1
        for flag in ('all_five_cumulative_moments_and_same_absolute_datum_retained',
            'zero_axial_input_does_not_zero_radial_history','original_analytic_P0_source_getter_retained',
            'P0_whole_end_box_only_bounds_the_exact_terminal_function',
            'original_G1_G2_positive_factors_not_materialized'):
            if not packet[flag]:raise ValueError('Whole gap source scope missing: '+flag)
        if packet['numerical_caps_used_as_field_values']:raise ValueError('Field replaced by cap')
        for flag in FALSE_FLAGS:
            if packet[flag] or record[flag]:raise ValueError('Similarity scope overclaimed: '+flag)
    fixture=independent_gap_fixture()
    hashes=dict(companion.hashes);hashes[name]=hashlib.sha256(raw).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(all_passed=True,actual_five_defect_family_sha256=companion.family,
        implicit_source_sha256=companion.source,domain=DOMAIN,
        source_function_identities_checked=len(identities),finite_signed_rows_checked=counts,
        structural_zero_rows_checked=zeros,current_whole_gap_report_recomputed=True,
        original_uncapped_logE_reduced_source_binding_verified=True,
        same_canonical_P0_source_bound_not_chosen_as_field=True,
        full_nonzero_raw_moment_radial_and_stress_histories_preserved=True,
        actual_original_whole_inactive_gap_similarity_companion_constructed=True,
        gap_end_similarity_velocity4_moment4_stress3_pressure4_functional_join_verified=True,
        main_gap_and_gap_end_original_coordinate_interfaces_consumed=True,
        independent_gap_fixture=fixture,input_hashes=hashes,
        source_caps_used_as_defining_field_values=False,**{flag:False for flag in FALSE_FLAGS})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('PASS original whole inactive-gap similarity companion and exact end source join; physical/cone pending',flush=True)
    return result


if __name__=='__main__':run()
