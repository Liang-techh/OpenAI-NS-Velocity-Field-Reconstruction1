"""Independent convex comparison, actual moment and width-product checks."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_inner_bridge_profiles import (
    logarithm,cumulative_enclosures,direction,IntervalTaylor,MTH,MZ,MTHZ,MZT,MP)
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_inner_bridge_profiles.json'


def sigma(x):
    if x<=0:return mp.mpf(0)
    if x>=1:return mp.mpf(1)
    a=mp.exp(-1/x**2);b=mp.exp(-1/(1-x)**2)
    return a/(a+b)


def comparison_and_chi_fixture():
    c=MPIntervalContext();c.dps=100;h=mp.mpf('.025');count=0;kernel=0
    with mp.workdps(80):
        alpha=lambda q:1-sigma(q-1)
        lo=mp.mpf(4);hi=4*mp.exp(2*h)
        # Polynomial in Z for the log core and V core. Integrate the
        # original alpha*g_y equations; compare against positive IBP bounds.
        for z in (mp.mpf(0),mp.mpf('.4')):
            rho=c.mpf([lo,hi]);zz=IntervalTaylor.variable(c,c.mpf(z),6)
            logcore=(zz*mp.mpf('.07')+zz*zz*mp.mpf('.02')+mp.mpf('-.05'))*rho
            corev=zz*2+mp.mpf('.3')+(zz*zz*mp.mpf('.02')+mp.mpf('.07'))*rho
            barphi=logcore.exp()
            for phase in (mp.mpf(0),mp.mpf('.5'),mp.mpf(1),mp.mpf('1.5'),mp.mpf(2),mp.mpf(3)):
                stop=min(phase,2)
                pieces=[mp.mpf(0)]+([mp.mpf(1)] if stop>1 else [])+[stop]
                eff=4+mp.quad(lambda q:4*h*alpha(q)*mp.exp(h*q),pieces) if stop else mp.mpf(4)
                f=lambda q:mp.exp(eff*(-mp.mpf('.05')+mp.mpf('.07')*q+mp.mpf('.02')*q*q))
                v=lambda q:mp.mpf('.3')+2*q+eff*(mp.mpf('.07')+mp.mpf('.02')*q*q)
                for k in range(7):
                    for jet,target in ((barphi,mp.diff(f,z,k)),(corev,mp.diff(v,z,k))):
                        left,right=endpoints(jet[k]*math.factorial(k))
                        if not left-mp.mpf('1e-65')<=target<=right+mp.mpf('1e-65'):
                            raise ArithmeticError('Independent smooth comparison IBP enclosure failed')
                        count+=1
        chi=lambda y:1-(1-h)*sigma(y/h)
        for y in (mp.mpf(0),h/2,h,2*h,mp.mpf('.5'),mp.mpf(1)):
            end=min(y,h)
            mass=mp.quad(chi,[0,end])+h*max(0,y-h) if y else mp.mpf(0)
            if not -mp.mpf('1e-65')<=mass<=h*(1+y)+mp.mpf('1e-65'):
                raise ArithmeticError('Exact integrated chi width bound failed')
            kernel+=1
        return dict(independent_comparison_axial_derivatives=count,independent_integrated_chi_bounds=kernel,
            finite_parameter_fixture_only=True,actual_source_admission=False,passed=True)


def variable_profile_moment_fixture():
    c=MPIntervalContext();c.dps=100;tol=mp.mpf('1e-60')
    R,x,z=s.symbols('R x z',real=True);r=s.Rational(4,5)
    phi=1+R/25+z/10+R*z*z/100;v=s.Rational(3,10)+2*z+7*R/100+R*z*z/50
    integrands=[2*R*phi,v,2*R*phi*v,v*v,R*phi*phi,phi*phi]
    raw=[s.integrate(q,(R,0,x)).subs(x,R) for q in integrands]
    shapes=[raw[0]/R**2,raw[1]/R,raw[2]/R**2,raw[3]/R,raw[4]/R**2,raw[5]/R]
    inlet={n:s.lambdify(z,q.subs(R,r),'mpmath') for n,q in zip(('H','mean','K','A','B','C'),shapes)}
    target={n:[s.lambdify((R,z),s.diff(q,z,k),'mpmath') for k in range(7)]
            for n,q in zip(('H','mean','K','A','B','C'),shapes)}
    f0=s.exp(13*z/100+z*z/50);ps=s.Rational(7,5);p0=-2-z*z/10
    d=1-z*z;L=1-s.Rational(3,100)*z*z
    physical=[f0*raw[0],raw[1],f0*raw[2],raw[3]-f0*f0*raw[4],f0*f0*raw[5]]
    Aop=lambda q:(1-s.Rational(3,100))*z*q+d*s.diff(q,z)
    Pop=lambda q:2*(1+s.Rational(3,100))*z*q-d*s.diff(q,z)
    angular=(1-s.Rational(3,200))*physical[0]-(1-s.Rational(3,100))*z*s.diff(physical[0],z)/2
    angular-=d*s.diff(physical[2],z)-(2*s.Rational(3,100)-1)*z*physical[2]
    D_R=(-R+Aop(physical[1]))/(L*R)+angular/(2*L*R**2*f0*phi)
    drive=(-R+Aop(physical[1]))*v+(1-s.Rational(3,100))*(physical[1]-z*s.diff(physical[1],z))/2
    drive+=2*s.Rational(3,100)*z*physical[3]-d*s.diff(physical[3],z)+R*Pop(ps**2*p0+physical[4])
    drive/=2*L
    dir_targets=[[s.lambdify((R,z),s.diff(q,z,k),'mpmath') for k in range(6)] for q in (D_R,drive)]
    amplitude=s.lambdify(z,f0,'mpmath');p0fn=s.lambdify(z,p0,'mpmath')
    moment_count=0;direction_count=0
    for zz in (mp.mpf(0),mp.mpf('.4')):
        initial={n:IntervalTaylor(c,[c.mpf([a-tol,a+tol])/math.factorial(k)
                 for k in range(7) for a in (mp.diff(fn,zz,k),)]) for n,fn in inlet.items()}
        for radius in (mp.mpf('.8'),mp.mpf(1),mp.mpf(4)):
            rho=c.mpf([mp.mpf('.8'),radius]);jz=IntervalTaylor.variable(c,c.mpf(zz),6)
            pc=1+jz/10+(jz*jz/100+mp.mpf('.04'))*rho
            vc=mp.mpf('.3')+jz*2+(jz*jz/50+mp.mpf('.07'))*rho
            data=cumulative_enclosures(c,c.mpf(mp.mpf('.8')/radius),initial,pc,vc)
            jets=dict(H=data[MTH],mean=data[MZ],K=data[MTHZ],A=data[MZT]['axial'],B=data[MZT]['swirl'],C=data[MP])
            for name,jet in jets.items():
                for k in range(7):
                    a,b=endpoints(jet[k]*math.factorial(k));expected=target[name][k](radius,zz)
                    if not a-tol*100<=expected<=b+tol*100:raise ArithmeticError('Actual variable-profile cumulative moment enclosure failed')
                    moment_count+=1
            pnow=1+jz/10+(jz*jz/100+mp.mpf('.04'))*radius
            vnow=mp.mpf('.3')+jz*2+(jz*jz/50+mp.mpf('.07'))*radius
            ratios=[c.mpf(mp.diff(amplitude,zz,k)/amplitude(zz)) for k in range(7)]
            ratios2=[c.mpf(mp.diff(lambda q:amplitude(q)**2,zz,k)/amplitude(zz)**2) for k in range(7)]
            pressure=IntervalTaylor(c,[c.mpf(mp.diff(p0fn,zz,k))/math.factorial(k) for k in range(7)])
            bounds=direction(c,c.mpf(zz),c.mpf('.03'),pnow,vnow,data,pressure,ratios,ratios2)
            for k in range(6):
                values=(bounds['D_over_R'][k],
                    (bounds['drive_hydro'][k]+bounds['drive_pressure'][k]*mp.mpf('1.96'))*radius
                    +bounds['drive_swirl'][k]*(radius**2*amplitude(zz)**2))
                for interval,fn in zip(values,dir_targets):
                    a,b=endpoints(interval*math.factorial(k));expected=fn[k](radius,zz)
                    if not a-tol*100<=expected<=b+tol*100:raise ArithmeticError('Independent own-moment direction axial5 failed')
                    direction_count+=1
    return dict(independently_integrated_variable_profile_moment_derivatives=moment_count,
        independent_own_moment_inertial_direction_derivatives=direction_count,
        finite_parameter_fixture_only=True,actual_source_admission=False,passed=True)


def structural_identities():
    R,z,delta=s.symbols('R z delta',real=True);Mz=s.Function('Mz')(R,z)
    V=s.diff(Mz,R);m=Mz/R;d=1-z*z;L=1-delta*z*z
    Q=(2*z*V-(1-delta)*z*m-d*s.diff(m,z))/L
    divergence=L*(Q+R*s.diff(Q,R))-(1+delta)*z*V+d*s.diff(V,z)-2*z*R*s.diff(V,R)
    if s.simplify(divergence)!=0:raise ArithmeticError('Actual cumulative mean divergence identity failed')
    phi,barphi,F0,I=s.symbols('phi barphi F0 Iz',positive=True)
    if s.simplify((F0*phi)*(I/(F0*barphi))-(phi/barphi)*I)!=0:
        raise ArithmeticError('Actual-to-comparison axial amplitude cancellation failed')
    return dict(exact_actual_mean_divergence=1,exact_axial_amplitude_cancellation=1,passed=True)


def run():
    with mp.workdps(280):
        receipt=json.loads((HERE/NAME).read_bytes())
        for name,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Bridge source changed: '+name)
        c=MPIntervalContext();c.dps=240
        sign=receipt['angular_log_correction_nonpositive_source']
        certificate=json.loads((HERE/sign['certificate']).read_bytes())
        if (not certificate['exact_implicit_exit_field_specified']
                or not certificate['actual_whole_axis_Ra_R110_relaxed_cone_analytically_certified']
                or certificate['admitted_inner_parameter_family_sha256']!=sign['admitted_inner_parameter_family_sha256']):
            raise ValueError('Angular log sign lacks accepted same-source comparison positivity')
        logh=read_interval(c,receipt['source_log_hb_enclosure']);proof_count=0;bound_count=0
        for proof in receipt['source_log_product_cap_proofs']:
            upper=read_interval(c,proof['input_absolute_upper']);additional=read_interval(c,proof['additional_source_log'])
            cap=read_interval(c,proof['cap']);computed=logh+additional+c.ln(upper)
            if endpoints(computed)[1]>endpoints(c.ln(cap))[0] or endpoints(cap)[0]<=0:
                raise ArithmeticError('Width/product cap not proved from actual source logs')
            proof_count+=1
        packets=[receipt[n] for n in ('whole_bridge','exit','terminal')]+receipt['samples']
        for packet in packets:
            for name in ('F_actual_over_F0_axial5_coefficients','F_actual_true_axial5_divided_by_F0',
                         'Uz_actual_axial5_coefficients','actual_Q_axial4_coefficients',
                         'pressure_axis_axial5_coefficients','pressure_increment_true_axial5_divided_by_R_F0_squared'):
                expected=5 if name=='actual_Q_axial4_coefficients' else 6
                if len(packet[name])!=expected:raise ValueError('Bridge axial derivative order lost')
                for value in packet[name]:
                    lo,hi=endpoints(read_interval(c,value))
                    if lo>hi or not mp.isfinite(lo) or not mp.isfinite(hi):raise ArithmeticError('Nonfinite bridge profile enclosure')
                    bound_count+=1
            if endpoints(read_interval(c,packet['F_actual_over_F0_axial5_coefficients'][0]))[0]<=0:
                raise ArithmeticError('Actual bridge positive swirl was lost')
            if packet['actual_bridge_radial_mixed4_certified']:raise ValueError('Unbuilt radial mixed derivatives were promoted')
        exit=receipt['exit']
        if not exit['core_exit_exact'] or any(endpoints(read_interval(c,v))!=(mp.mpf(0),mp.mpf(0))
                for v in exit['log_F_actual_over_exit_F_coefficients']):raise ArithmeticError('Actual bridge exit history changed')
        result=dict(actual_five_defect_family_sha256=receipt['actual_five_defect_family_sha256'],implicit_source_sha256=receipt['implicit_source_sha256'],
            comparison_fixture=comparison_and_chi_fixture(),moment_fixture=variable_profile_moment_fixture(),structural_identities=structural_identities(),
            actual_source_width_product_log_proofs_checked=proof_count,actual_bridge_profile_bounds_checked=bound_count,
            angular_log_sign_bound_to_accepted_same_source_comparison_theorem=True,
            source_bound_smoothed_comparison_axial6_enclosures_available=True,
            actual_prescribed_shear_bridge_axial5_enclosures_available=True,
            actual_bridge_radial_recovery_axial4_enclosures_available=True,
            actual_bridge_radial_mixed4_certified=False,bridge_100_110_switches_installed=False,
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,temporal_recursion=False,
            all_passed=True,input_hashes={**receipt['input_hashes'],NAME:hashlib.sha256((HERE/NAME).read_bytes()).hexdigest(),
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Actual width-source log products, convex comparison and actual cumulative moments/divergence PASS',flush=True)
    return result


if __name__=='__main__':run()
