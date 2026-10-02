"""Independent post-pulse normalization, history and Gamma-tail checks."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_corrected_outer_field import CompliantCorrectedOuterField
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_corrected_outer_field.json'


def identities():
    passed={}
    def zero(name,expression):
        if s.simplify(expression)!=0:raise ArithmeticError('Corrected outer identity failed: '+name)
        passed[name]=True
    mu,a,t,J,T,w=s.symbols('mu a t J T w',real=True)
    rate=1-mu;k=1-a;bp=s.Rational(1,2)+mu
    zero('steep_in_energy_weight',s.exp(t)*s.exp(-2*bp*t-2*rate*J)-s.exp(-2*mu*t-2*rate*J))
    zero('steep_in_pressure_weight',s.exp(-2*bp*t-2*rate*J)/2-s.exp(-(1+2*mu)*t-2*rate*J)/2)
    zero('Ns_ratio',s.exp(1)*s.exp(-2*bp-rate)-s.exp(-1-mu))
    zero('Ps_ratio',s.exp(-2*bp-rate)-s.exp(-2-mu))
    zero('steep_power_energy_ratio',s.exp(T)*s.exp(-3*T)-s.exp(-2*T))
    zero('steep_power_pressure_ratio',s.exp(-s.Rational(3,2)*T)**2-s.exp(-3*T))
    zero('steep_out_energy_weight',s.exp(t)*s.exp(-3*t+2*k*J)-s.exp(-2*t+2*k*J))
    zero('steep_out_pressure_weight',s.exp(-3*t+2*k*J)/2-s.exp(-3*t+2*k*J)/2)
    zero('Nt_ratio',s.exp(1)*s.exp(-3+k)-s.exp(-1-a))
    zero('Pt_ratio',s.exp(-3+k)-s.exp(-2-a))
    zero('waiting_energy_ratio',s.exp(w)*s.exp(-(1+2*a)*w)-s.exp(-2*a*w))
    zero('waiting_pressure_ratio',s.exp(-(s.Rational(1,2)+a)*w)**2-s.exp(-(1+2*a)*w))
    y=s.symbols('y',real=True);U=s.Function('U')(y);M=s.Function('M')(y);W=s.Function('W')(y);Q=s.Function('Q')(y)
    nx=s.exp(s.Rational(3,2)*y)*U
    zero('forward_angular_normalized_ODE',s.diff(M/nx,y).subs(s.diff(M,y),nx)-(1-(s.Rational(3,2)+s.diff(U,y)/U)*M/nx))
    e=W/(2*s.exp(y)*U**2)
    zero('remaining_energy_moment_ODE',s.diff(e,y).subs(s.diff(W,y),-s.exp(y)*U**2)-(-s.Rational(1,2)-(1+2*s.diff(U,y)/U)*e))
    p=Q/U**2
    zero('remaining_pressure_normalized_ODE',s.diff(p,y).subs(s.diff(Q,y),-U**2/2)-(-s.Rational(1,2)-2*s.diff(U,y)/U*p))
    C,K=s.symbols('C K',real=True)
    zero('angular_defect_is_retained_constant',(C*s.exp(-k*t)/K)*(K*s.exp(k*t))-C)
    p0,mrv,prv,pr,scale=s.symbols('p0 mrv prv pr scale',real=True)
    forward=p0+mrv+scale*(prv-pr);backward=-scale*pr
    zero('forward_backward_pressure_offset',forward-backward-(p0+mrv+scale*prv))
    eps,BW,D=s.symbols('epsilon BW D',real=True)
    zero('heat_square_deficit_retained',(1-eps*BW)**2-(1-eps*BW-D)**2-D*(2*(1-eps*BW)-D))
    zero('separate_epsilon_squared_atom',(1-eps*BW)**2-(1-2*eps*BW+eps**2*BW**2))
    q=s.symbols('q',positive=True);sig=s.symbols('sig',real=True)
    zero('flatten_Ev_to_Ev0_units',q**2*(q**-2*(q/2)**(2*sig))-(q/2)**(2*sig))
    # Exact Gamma angular primitive: derivative of the expression inside
    # its positive Gamma expectation gives the required renormalized tail.
    x=s.symbols('x',positive=True)
    candidate=((1+x)**(1-a)-1)/(1-a)
    zero('Gamma_angular_integrated_identity',s.diff(candidate,x)-(1+x)**(-a))
    return passed


def gamma_fixture():
    """Direct hypergeometric/Gamma evaluation, independent of Taylor bounds."""
    c=MPIntervalContext();c.dps=65
    f=object.__new__(CompliantCorrectedOuterField);f.ctx=c;f.delta=c.mpf('.2');f.a=f.delta/2
    f.k=1-f.a;f.phrate=1+f.delta;f.repair=SimpleNamespace(strong_S_cap=c.mpf('.002'))
    f.pulse=SimpleNamespace(factor=lambda logvalue:c.exp(logvalue))
    with mp.workdps(90):
        Z=mp.mpf('.5');t=mp.mpf(3);a=mp.mpf('.1');delta=2*a;S=mp.mpf('.002')*mp.exp(-t)
        xi=2*(1-Z**2)*S;result=f.heat_local(c.mpf(Z),c.mpf(t))
        def H(x):return x**(-1-a)*mp.hyperu(1+a,2,1/x)
        def Hp(x):return -a*(1+a)*x**(-2-a)*mp.hyperu(2+a,2,1/x)
        hv=H(xi);hxi=-a*(1+a)*xi**(-2-a)*mp.hyperu(2+a,2,1/xi)
        az=mp.quad(lambda v:mp.exp(-v)*v**a*(1+xi*v)**(1-a),[0,1,mp.inf])/((1-a)*mp.gamma(1+a))
        refs={'K_value':(result['K'][0],hv),'K_Z':(result['K'][1],-4*Z*S*hxi),'angular_target':(result['X_target'][0],az/hv)}
        azi=-4*Z*S*(1+a)*xi**(-2-a)*mp.hyperu(2+a,3,1/xi)
        hvz=-4*Z*S*hxi
        refs['angular_target_Z']=(result['X_target'][1],(azi*hv-az*hvz)/hv**2)
        for label,rate,key in (('pressure',1+delta,'pressure_tail_over_Utheta_squared'),('energy',delta,'energy_tail_over_R_Utheta_squared')):
            L=mp.mpf(30)
            base=mp.quad(lambda u:mp.exp(-rate*u)*H(xi*mp.exp(-u))**2,[0,1,5,15,L])+mp.exp(-rate*L)/rate
            tailerror=2*a*(1+a)*xi*mp.exp(-(rate+1)*L)/(rate+1)
            if label=='pressure':base/=2;tailerror/=2
            bz=mp.quad(lambda u:mp.exp(-(rate+1)*u)*2*H(xi*mp.exp(-u))*Hp(xi*mp.exp(-u)),[0,1,5,15,L])*(-4*Z*S)
            zerror=8*a*(1+a)*abs(Z)*S*mp.exp(-(rate+1)*L)/(rate+1)
            if label=='pressure':bz/=2;zerror/=2
            derivative=(bz*hv-2*base*hvz)/hv**3
            refs[label+'_Z']=(result[key][1],derivative)
            derivative_remainder=zerror/hv**2+2*tailerror*abs(hvz)/hv**3
            dzlo,dzhi=endpoints(result[key][1])
            if not dzlo<=derivative-derivative_remainder<=derivative+derivative_remainder<=dzhi:
                raise ArithmeticError('Independent Gamma derivative outside enclosure: '+label)
            refs[label]=(result[key][0],base/hv**2)
            lo,hi=endpoints(result[key][0])
            if not lo<=(base-tailerror)/hv**2<=base/hv**2<=hi:
                raise ArithmeticError('Independent Gamma integral outside directed enclosure: '+label)
        passed={}
        for label,(interval,value) in refs.items():
            lo,hi=endpoints(interval)
            if not lo<=value<=hi:raise ArithmeticError('Gamma fixture failed: '+label)
            passed[label]=True
        return dict(delta=str(delta),inverse_radius_source=str(S),direct_Gamma_and_tail_checks=passed,
                    finite_tail_remainder_bounded=True,passed=True)


def run():
    r=json.loads((HERE/NAME).read_bytes());hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Corrected outer source changed: '+name)
    exact=identities();c=MPIntervalContext();c.dps=240;get=lambda v:read_interval(c,v)
    samples=r['samples'];interfaces={}
    keys=('formal_log_Utheta_over_Ev0','Mtheta_over_sqrt2_R_3half_Utheta',
          'Mztheta_over_R_Utheta_squared','Mp_over_Pstar_squared','P_over_Pstar_squared')
    for li,ri,label in ((2,3,'flatten_power'),(4,5,'power_bumps'),(9,10,'bumps_steepin'),(11,12,'steepin_power'),
                        (14,15,'power_steepout'),(16,17,'steepout_waiting'),(18,19,'waiting_collar')):
        for key in keys:
            for v,w in zip(samples[li][key]['coefficients'],samples[ri][key]['coefficients']):
                a,b=endpoints(get(v));d,e=endpoints(get(w))
                if max(a,d)>min(b,e):raise ArithmeticError('Outer interface disjoint: '+label+' '+key)
        interfaces[label]=True
    allpoints=samples+r['whole_Z_examples']+[r['angular_interval_crossing_support'],r['whole_Z_Rv_interface']]
    for point in allpoints:
        for key in ('Uz','Ur','Mz_over_R_Utheta','Mtheta_z_over_sqrt2_R_3half_Utheta_squared'):
            for value in point[key]['coefficients']:
                if endpoints(get(value))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Pulse terminal zero history not inherited')
        for key in keys:
            if len(point[key]['coefficients'])!=2:raise ArithmeticError('Corrected outer lost C1: '+key)
        if endpoints(get(point['Mztheta_over_R_Utheta_squared']['coefficients'][0]))[0]<=0:
            raise ArithmeticError('Positive remaining swirl energy lost')
        if endpoints(get(point['independent_backward_pressure_over_Utheta_squared']['coefficients'][0]))[1]>=0:
            raise ArithmeticError('Negative backward pressure diagnostic lost')
        for p0,mpv,p in zip(point['P0_over_Pstar_squared']['coefficients'],point['Mp_over_Pstar_squared']['coefficients'],point['P_over_Pstar_squared']['coefficients']):
            expect=get(p0)+get(mpv);actual=get(p)
            if not endpoints(actual)[0]<=endpoints(expect)[0]<=endpoints(expect)[1]<=endpoints(actual)[1]:
                raise ArithmeticError('Original P0+Mp pressure source replaced')
        if point['full_global_pressure_terminal_identity_verified'] or point['full_outer_five_moment_match'] or point['whole_outer_cone_certified'] or point['temporal_recursion']:
            raise ValueError('Corrected outer scope overclaimed')
    for actual,expected in zip(r['post_Rv_energy_tail_in_Ev0_units']['coefficients'],r['pulse_future_energy_in_Ev0_units']['coefficients']):
        a,b=endpoints(get(actual));d,e=endpoints(get(expected))
        if max(a,d)>min(b,e):raise ArithmeticError('Complete post-Rv energy no longer agrees with selected pulse source')
    for current,pulse in ((samples[0],r['source_pulse_terminal']),(r['whole_Z_Rv_interface'],r['whole_Z_pulse_terminal'])):
        for v in current['forward_Mp_increment_in_Ev0_squared']['coefficients']:
            if endpoints(get(v))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Empty Rv pressure increment is not exactly zero')
        for value,source in zip(current['Mp_over_Pstar_squared']['coefficients'],pulse['pressure']['Mp_over_Pstar_squared']['coefficients']):
            if endpoints(get(value))!=endpoints(get(source)):raise ArithmeticError('Whole-Z incoming pressure history changed at Rv')
    if r['full_global_pressure_terminal_identity_verified'] or r['forward_angular_heat_target_equality_verified']:
        raise ValueError('Unproved absolute terminal identity claimed')
    fixture=gamma_fixture();hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],implicit_source_sha256=r['implicit_source_sha256'],
        independent_outer_identities=exact,directed_interfaces=interfaces,independent_Gamma_fixture=fixture,
        whole_Z_zero_histories_and_positive_energy_checked=True,unchanged_forward_P0_plus_Mp_checked=True,
        selected_pulse_future_energy_source_agreement_checked=True,all_passed=True,
        post_Rv_corrected_velocity_and_five_primitives_callable=True,
        full_global_pressure_terminal_identity_verified=False,forward_angular_heat_target_equality_verified=False,
        full_outer_five_moment_match=False,C4_outer_bounds_certified=False,whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Corrected post-Rv outer:',len(exact),'independent identities;',len(interfaces),'interfaces; inherited whole-Z histories, pressure source and direct Gamma fixtures PASS',flush=True)
    return result


if __name__=='__main__':run()
