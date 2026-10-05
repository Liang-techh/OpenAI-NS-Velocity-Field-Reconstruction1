"""Focused continuous main/exit cone and independent original-integral tests."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_main_exit_cone import (
    PREFIX,DOMAIN,TAIL_DOMAIN,ADMISSIONS,FALSE_FLAGS,CORRECTIONS,current_sources,
    main_cone_identities,whole_main_exit_bounds,relative_log_recipes,relative_log_envelopes,
    source_precision)
from lei_ren_part1_paper_compliant_pulse_main_exit_fixture_integrals import run as load_fixture_integrals
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

HERE=Path(__file__).parent


def independent_original_main_exit_fixture():
    """Original gp convolution, nonzero shear and full signed stress.

    Moderate smooth fixtures test formulas/units and attained cone tests;
    continuous actual-source admission uses the separate analytic proof.
    All startup history is integrated by Fubini, without deleting a tail.
    """
    with mp.workdps(100):
        mu,delta=mp.mpf('.001'),mp.mpf('1e-10')
        r=1-mu;lam1=mp.mpf('.5')-mu;p=1+2*mu;bp=mp.mpf('.5')+mu
        Rp,Pstar,D0,Xp=mp.exp(100),mp.mpf('.45'),mp.mpf('.07'),mp.mpf('.3')
        left=mp.mpf('.02');tol=mp.mpf('1e-60')
        def sigma(x):
            if x<=0:return mp.mpf(0)
            if x>=1:return mp.mpf(1)
            return 1/(1+mp.exp(1/x**2-1/(1-x)**2))
        def sj(x):
            if x<=0:return (mp.mpf(0),)*3
            if x>=1:return (mp.mpf(1),mp.mpf(0),mp.mpf(0))
            v=sigma(x);L1=2/(1-x)**3+2/x**3;L2=6/(1-x)**4-6/x**4
            return v,v*(1-v)*L1,v*(1-v)*(L2+(1-2*v)*L1**2)
        def gj(x):
            if x>=11:return (mp.mpf(0),)*3
            v,d,dd=sj(11-x);primitive=x-mp.mpf('.01')
            return primitive*v,v-primitive*d,-2*d+primitive*dd
        startup=[];startup2=[]
        for rate in (lam1,mp.mpf('.5')-2*mu):
            k=rate/mu
            startup.append(mp.quad(lambda a:sigma(50*a)*(mp.exp(k*left)-mp.exp(k*a))/k,
                [0,left/2,left]))
            startup2.append(mp.quad(lambda a:mp.exp(k*a)*50*sj(50*a)[1],[0,left/2,left]))
        def kernel(x,i):
            rate=mp.mpf('.5')-i*mu
            head=startup[i-1]/mu-mp.exp(rate*left/mu)*((left-mp.mpf('.01'))/rate-mu/rate**2)
            main=lambda xx:(xx-mp.mpf('.01'))/rate-mu/rate**2+mp.exp(-rate*xx/mu)*head
            if x<=10:return main(x)
            return mp.exp(-rate*(x-10)/mu)*main(mp.mpf(10))+mp.quad(
                lambda a:mp.exp(rate*(a-x)/mu)*gj(a)[0],[10,(10+x)/2,x])/mu
        def second_kernel(x,i):
            rate=mp.mpf('.5')-i*mu
            result=mp.exp(-rate*x/mu)*startup2[i-1]/mu
            if x>10:result+=mp.quad(lambda a:mp.exp(rate*(a-x)/mu)*gj(a)[2],[10,(10+x)/2,x])/mu
            return result
        fixture=load_fixture_integrals();vals=fixture['values']
        energy_start=mp.mpf(vals['entrance_energy'])
        energy_primitive=lambda a:-mp.exp(-2*a)*((a-mp.mpf('.01'))**2/2
            +(a-mp.mpf('.01'))/2+mp.mpf('.25'))
        def partial_energy(x):
            result=energy_start+energy_primitive(min(x,mp.mpf(10)))-energy_primitive(left)
            if x>10:result+=mp.quad(lambda a:mp.exp(-2*a)*gj(a)[0]**2,[10,(10+x)/2,x])
            return result
        C=lambda z:1/(1+z*z)
        ap=lambda z:mp.mpf('1.02')+mu*mp.mpf('.02')*z*z
        incoming=(lambda z:mu**2*mp.mpf('.1')*(z+z**3),
            lambda z:mu**2*mp.mpf('.07')*(z+z**3))
        ein=lambda z:mp.mpf('.015')+mp.mpf('.002')*z*z
        P0=lambda z:-mp.mpf('.35')*C(z)**2+mp.mpf('.02')*z
        c=MPIntervalContext();c.dps=120
        lp=mp.log(Pstar);lu=mp.mpf(0);lrp=mp.log(Rp);finite=mp.log(D0)+1/mu
        envelopes=relative_log_envelopes(c,c.mpf(mu),c.mpf(finite),c.mpf(lp),c.mpf(lu),c.mpf(lrp))
        checks={};worst=mp.mpf(0);nonzero=falling=0;logs=0;measured=[]
        def compare(name,a,b):
            nonlocal worst
            error=abs(a-b)/max(1,abs(a),abs(b));worst=max(worst,error)
            if error>tol:raise ArithmeticError('Independent full main/exit identity differs: '+name)
            checks[name]=True
        for x,z in ((mp.mpf('.5'),mp.mpf(0)),(mp.mpf(2),mp.mpf('-.4')),
            (mp.mpf('10.4'),mp.mpf('.31')),(mp.mpf('10.8'),mp.mpf(0))):
            y=x/mu;R=Rp*mp.exp(y);B=Pstar*mp.exp(-bp*y);H=mp.exp(-r*y)
            g,gd,gdd=gj(x);kk=[kernel(x,i) for i in (1,2)]
            for i in (1,2):
                rate=mp.mpf('.5')-i*mu
                recovered=g/rate-mu*gd/rate**2+mu**2*second_kernel(x,i)/rate**2
                compare('full_original_two_IBP/'+str(x)+'/'+str(i),kk[i-1],recovered)
            energy=partial_energy(x)
            def data(zz):
                cc=C(zz);ut=B*cc;aa=ap(zz);uz=ut*aa*g
                moments_i=[aa*kk[i-1]+mp.exp(-(mp.mpf('.5')-i*mu)*y)*incoming[i-1](zz) for i in (1,2)]
                ev=mp.exp(2*x)*(ein(zz)+aa*aa*energy/mu-(1-mp.exp(-2*x))/(4*mu))
                pressure=B*B*(-cc*cc/(2*p)+mp.exp(-p*(13-x)/mu)*(P0(zz)+cc*cc/(2*p)))
                moments=dict(theta=mp.sqrt(2)*R**mp.mpf('1.5')*ut*(1/r+(Xp-1/r)*H),
                    z=R*ut*moments_i[0],theta_z=mp.sqrt(2)*R**mp.mpf('1.5')*ut**2*moments_i[1],
                    z_theta=R*ut**2*ev,p=mp.mpf(0))
                return ut,uz,pressure,moments
            Ut,Uz,P,moments=data(z)
            actual=evaluate_mp_stress(mp.log(R),z,delta,Utheta=Ut,Uz=Uz,Utheta_y=-bp*Ut,
                Utheta_Z=B*mp.diff(C,z),Uz_y=B*C(z)*ap(z)*(mu*gd-bp*g),
                Uz_Z=B*mp.diff(lambda zz:C(zz)*ap(zz),z)*g,
                moments=moments,moments_Z={key:mp.diff(lambda zz:data(zz)[3][key],z) for key in moments},
                P=P,P_Z=mp.diff(lambda zz:data(zz)[2],z),precision=mp.mp.dps,radius_override=R,axial_override=z)
            a=2+2*mu;b=2*mu*ap(z)*gd-(1+2*mu)*ap(z)*g
            F=Ut/mp.sqrt(2*R);tt,tz=actual['T_theta'],actual['T_z']
            compare('full_theta_shear/'+str(x),actual['S_theta']/F,-a)
            compare('full_axial_shear/'+str(x),actual['S_z']/F,b)
            if abs(b)>mp.mpf('1e-10'):nonzero+=1
            if gd<0:falling+=1
            w=tz/tt;kap=a+b*b/a;expr=2*b*w+b*b/a+(a-2)*w*w
            dot=tt*actual['S_theta']+tz*actual['S_z']
            cross=-tt*actual['S_z']+tz*actual['S_theta']
            compare('negative_dot/'+str(x),-dot/(tt*F),a-b*w)
            compare('full_directional_factor/'+str(x),
                (2*dot**2-(kap-2)*cross**2)/(tt*F)**2,(a*a+b*b)*(2-expr))
            if tt<=0 or dot>=0 or expr>=2 or kap<=2:raise ArithmeticError('Moderate full-shear fixture cone failed')
            q=z*z/(1+z*z);q0=mu-delta/2+(1-delta)*q
            K=r*mu*(1-delta)*(mp.mpf('.5')+q)/(lam1**2*q0)
            local_m=ap(z)*kk[0]
            local_N=-ap(z)*g+(1-delta)*(mp.mpf('.5')+q)*local_m
            local_N-=(1-delta)*z*mp.diff(ap,z)*kk[0]/2
            wlocal=local_N*r/q0;rho=kk[0]-g/lam1+mu*gd/lam1**2
            expected=2*ap(z)*g-K*ap(z)*gd+mu/lam1*ap(z)*g
            expected+=r*(1-delta)*(mp.mpf('.5')+q)*ap(z)*rho/q0
            expected-=r*(1-delta)*z*mp.diff(ap,z)*kk[0]/(2*q0)
            compare('exact_correlated_local_w/'+str(x),wlocal,expected)
            for nu,physical_lambda in ((mp.mpf('.01'),mp.mpf('.8')),(mp.mpf('.7'),mp.mpf('1.2'))):
                factor=physical_lambda**(-2-delta);S=[nu*factor*actual[key] for key in ('S_theta','S_z')]
                T=[nu*factor*v for v in (tt,tz)];Fp=factor*F
                kp=-(S[0]**2+S[1]**2)/(nu*Fp*S[0])
                compare('physical_kappa/'+str(x)+'/'+str(nu),kp,kap)
                dp=T[0]*S[0]+T[1]*S[1];cp=-T[0]*S[1]+T[1]*S[0]
                compare('physical_direction/'+str(x)+'/'+str(nu),
                    (2*dp**2-(kp-2)*cp**2)/((nu*factor)**4*(tt*F)**2),(a*a+b*b)*(2-expr))
            actual_logs=relative_log_recipes(c,c.mpf(mu),c.mpf(finite),c.mpf(lp),c.mpf(lu),c.mpf(lrp),c.mpf(x))
            for key,value in actual_logs.items():
                if endpoints(value-envelopes[key])[1]>mp.mpf('1e-90'):
                    raise ArithmeticError('Grouped source envelope differs: '+key)
                logs+=1
            measured.append(dict(xi=str(x),Z=str(z),b=mp.nstr(b,20),w=mp.nstr(w,20),
                a_minus_bw=mp.nstr(a-b*w,20),directional_bracket=mp.nstr(2-expr,20)))
            print('Independent original full-shear main/exit fixture xi='+str(x)+' checked',flush=True)
        return dict(all_passed=True,comparisons=len(checks),identities=checks,log_envelope_comparisons=logs,
            nonzero_axial_shear_fixtures=nonzero,falling_side_fixtures=falling,measured_fixture_tests=measured,
            tolerance=str(tol),maximum_normalized_error=mp.nstr(worst,30),
            original_startup_integrated_by_Fubini_and_complete_exit_integrals=True,
            actual_original_full_stress_evaluator_used=True,all_moment_pressure_and_signed_histories_retained=True,
            fixture_points_are_not_continuous_actual_source_proof=True)


@source_precision
def run():
    records,hashes,family=current_sources()
    name=PREFIX+'pulse_main_exit_cone.json';raw=(HERE/name).read_bytes();record=json.loads(raw)
    if record['domain']!=DOMAIN or record['tail_domain']!=TAIL_DOMAIN:raise ValueError('Original cone scope changed')
    for path,digest in record['input_hashes'].items():
        if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Cone source changed: '+path)
    if record['exact_source_cone_proof']!=main_cone_identities(records):raise ValueError('Current cone source proof differs')
    if record['bounds']!=encode(pack(whole_main_exit_bounds(records))):raise ValueError('Current continuous cone bounds differ')
    c=MPIntervalContext();c.dps=300;counts={}
    for group in ('positive_parameter_margins','positive_monotonicity_margins','positive_shape_and_error_margins',
        'positive_log_margins','positive_algebraic_margins'):
        for key,value in record['bounds'][group].items():
            lo,hi=endpoints(read_interval(c,value))
            if lo<=0 or not mp.isfinite(hi):raise ArithmeticError('Current continuous margin failed: '+key)
        counts[group]=len(record['bounds'][group])
    for flag in ADMISSIONS:
        if not record[flag]:raise ValueError('Missing main/exit cone admission: '+flag)
    for flag in FALSE_FLAGS:
        if record[flag]:raise ValueError('Global scope overclaimed: '+flag)
    for flag in ('continuous_whole_main_exit_Z_domain_covered','full_order_one_axial_shear_and_falling_side_retained',
        'all_signed_theta_and_axial_corrections_retained','original_energy_cancellation_and_absolute_pressure_used'):
        if not record['bounds'][flag]:raise ValueError('Actual full source scope missing: '+flag)
    if record['bounds']['phase_samples_used_as_proof'] or record['bounds']['source_caps_used_as_defining_field_values']:
        raise ValueError('Samples or caps substituted as defining proof/field')
    oracle=independent_original_main_exit_fixture()
    hashes.update(record['input_hashes']);hashes[name]=hashlib.sha256(raw).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for suffix in ('.py','.json'):
        path=PREFIX+'pulse_main_exit_fixture_integrals'+suffix
        hashes[path]=hashlib.sha256((HERE/path).read_bytes()).hexdigest()
    result=dict(all_passed=True,actual_five_defect_family_sha256=family[0],implicit_source_sha256=family[1],
        domain=DOMAIN,tail_domain=TAIL_DOMAIN,input_hashes=hashes,
        exact_source_cone_identities=len(record['exact_source_cone_proof']['identities']),positive_margin_counts=counts,
        continuous_full_shear_main_exit_proof_recomputed=True,
        admitted_main_exit_physical_Cartesian_operator_receipt_consumed_without_rerun=True,
        independent_original_main_exit_fixture=oracle,
        **{flag:True for flag in ADMISSIONS},**{flag:False for flag in FALSE_FLAGS})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('PASS continuous original main/exit full-shear cone and gap/end/tail composition; entrance/global/recursion pending',flush=True)
    return result


if __name__=='__main__':run()
