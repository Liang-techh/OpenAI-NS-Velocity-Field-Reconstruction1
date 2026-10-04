"""Original sigmoid entry: same full future, stress mixed3 and pressure mixed4.

The right datum is the admitted original power phase zero. Angular recovery
retains the native positive forward history; energy and absolute pressure
use directed signed remaining integrals. No upstream stress/cone is admitted.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_steep_power_stress_C3 import CompliantSteepPowerStressC3
from lei_ren_part1_paper_compliant_steep_exit_stress_C3 import source_precision,quotient_rows
from lei_ren_part1_paper_compliant_steep_waiting_C4 import transition_kernels
from lei_ren_part1_paper_compliant_collar_stress_C3 import (
    collar_stress_rows,axial_derivative,shifted_rows)
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_heat_pressure_C4 import pressure_y_rows
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_axial_pulse_field import sigma_enclosure
from lei_ren_part1_paper_compliant_pulse_physical_bounds import P
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
DOMAIN=dict(t=[0,1],Z=[-1,1],offset='q=-wait-Ts-2+t',ordinary_derivative='d_q=d_t')


def entry_signed_kernels(c,t,KS,mu,delta,cells=128):
    """Whole remaining entry integrals; f=int_t^1 sigma=1/2-J(t).

    Traverse backwards with correlated positive length (1-t)/cells.
    Both v-t and 1-v retain their correlated affine parameter expressions.
    KS is the same original power inlet, never the later exit inlet KQ.
    """
    t=c.mpf(t); lo,hi=endpoints(t)
    if lo<0 or hi>1:raise ValueError('Original entry t in[0,1] required')
    if not isinstance(cells,int) or cells<1:raise ValueError('Positive cell count required')
    length=1-t; ds=length/cells; f=c.mpf(0); r=1-mu; a=delta/2; k=1-a; p=1+delta
    IA=c.mpf(0); IE=c.mpf(0); IP=c.mpf(0)
    if endpoints(length)[1]>0:
        for i in reversed(range(cells)):
            vleft=t+length*i/cells; vright=t+length*(i+1)/cells
            v=c.mpf([max(mp.mpf(0),endpoints(vleft)[0]),min(mp.mpf(1),endpoints(vright)[1])])
            nextf=f+ds*sigma_enclosure(c,v)
            fc=c.mpf([endpoints(f)[0],min(mp.mpf('.5'),endpoints(nextf)[1])])
            remaining=length*(1-c.mpf([i,i+1])/cells)
            exponent=(mu-a)*remaining+r*fc
            Km=(KS-1)+KS*c.expm1(exponent)
            Qm=(KS**2-1)+KS**2*c.expm1(2*exponent)
            distance=length*c.mpf([i,i+1])/cells
            IA+=ds*c.exp(k*distance)*Km
            IE+=ds*c.exp(-delta*distance)*Qm
            IP+=ds*c.exp(-p*distance)*Qm
            f=nextf
    if (lo,hi)==(mp.mpf(0),mp.mpf(0)):f=c.mpf('.5')
    else:f=c.mpf([endpoints(f)[0],min(mp.mpf('.5'),endpoints(f)[1])])
    return dict(f=f,angular=IA,energy=IE,pressure=IP,
                exact_f_definition='integral_t^1 original_sigma(v)dv = 1/2-original_J(t)',
                exact_K_definition='KS*exp((mu-delta/2)*(1-v)+(1-mu)*f(v))',
                exact_angular_definition='integral_t^1 exp(k*(v-t))*(K(v)-1)dv',
                exact_energy_definition='integral_t^1 exp(-delta*(v-t))*(K(v)^2-1)dv',
                exact_pressure_definition='integral_t^1 exp(-p*(v-t))*(K(v)^2-1)dv',
                positive_correlated_cell_lengths=True,signed_integrands=True,cells=cells)


def entry_shape(heat,terminal,t,kernels):
    c=heat.ctx; one=IntervalTaylor.constant(c,1,5); KS=terminal['KS']; r=1-heat.mu
    sig=sigma_jets(c,t)
    g=[heat.a-heat.mu-r*sig[0]]+[-r*sig[j]*math.factorial(j) for j in range(1,4)]
    exponent=(heat.mu-heat.a)*(1-t)+r*kernels['f']
    K=[one*(KS*c.exp(exponent))]
    for n in range(1,5):
        K.append(sum((K[n-1-j]*g[j]*math.comb(n-1,j) for j in range(n)),one*0))
    square=product_rows(K,K)
    Km=[one*((KS-1)+KS*c.expm1(exponent))]+K[1:]
    Qm=[one*((KS**2-1)+KS**2*c.expm1(2*exponent))]+square[1:]
    return dict(K_rows=K,K_defect_rows=Km,K_squared_defect_rows=Qm,
                log_K_rate_rows=g,original_sigma_jets=sig,remaining=1-t)


def entry_angular_rows(heat,X0,shape):
    """Actual native X history, differentiated by its original variable rate."""
    one=IntervalTaylor.constant(heat.ctx,1,5); g=shape['log_K_rate_rows']
    h=[(1-heat.mu)*(1-shape['original_sigma_jets'][0])]+g[1:]; X=[X0]
    for n in range(1,5):
        X.append((one if n==1 else one*0)-sum(
            (X[n-1-j]*h[j]*math.comb(n-1,j) for j in range(n)),one*0))
    return X


def entry_defect_rows(heat,terminal,t,kernels,shape,X0):
    """Same full moments: angular native correlation and signed quadratic FTC."""
    c=heat.ctx; one=IntervalTaylor.constant(c,1,5); ell=1-t
    Ad=shape['K_rows'][0]*X0-one/heat.k
    Ed=terminal['energy_defect_rows'][0]*c.exp(-heat.delta*ell)+one*kernels['energy']
    Pd=terminal['pressure_defect_rows'][0]*c.exp(-heat.prate*ell)+one*kernels['pressure']/2
    A=[Ad]; E=[Ed]; Pr=[Pd]; Km=shape['K_defect_rows']; Qm=shape['K_squared_defect_rows']
    for j in range(4):
        A.append(Km[j]-A[j]*heat.k)
        E.append(E[j]*heat.delta-Qm[j])
        Pr.append(Pr[j]*heat.prate-Qm[j]/2)
    return dict(angular_defect_rows=A,energy_defect_rows=E,pressure_defect_rows=Pr,
                K_defect_rows=Km,K_squared_defect_rows=Qm)


def entry_stress_rows(heat,shape,defects,Xrows,Z,q):
    """Original full-factor stress, retaining common K in angular cancellation."""
    rows=collar_stress_rows(heat,shape,defects,Z,q)
    z=IntervalTaylor.variable(heat.ctx,Z,5); one=IntervalTaylor.constant(heat.ctx,1,5)
    b=(1-heat.delta)/2; L=1-z*z*heat.delta
    H=[Xrows[j]*heat.k-z*axial_derivative(Xrows[j])*b-(one if j==0 else one*0) for j in range(4)]
    raw=[jet/L for jet in product_rows(shape['K_rows'][:4],H)]
    inertial=shifted_rows(raw,-heat.a)
    rows['theta_inertial']=inertial
    rows['theta']=[inertial[j]+rows['theta_shear'][j] for j in range(4)]
    rows['normalized_coefficient_rows']['theta'][0]=raw[0]+rows['theta_shear'][0]
    return rows


def entry_pressure_rows(heat,defects,shape,q):
    one=IntervalTaylor.constant(heat.ctx,1,5)
    current=one/(2*heat.prate)+defects['pressure_defect_rows'][0]
    return pressure_y_rows(shape['K_rows'],current,heat.prate,heat.pressure_scale,q)


def entry_transport_identities():
    a,mu,t,z=s.symbols('a mu t Z',real=True); delta=2*a; k=1-a; p=1+delta; r=1-mu
    sig=s.Function('original_sigma')(t); f=s.Function('remaining_sigma')(t); J=s.Function('original_J')(t)
    KS=s.Symbol('same_original_power_inlet_KS',positive=True)
    K=KS*s.exp((mu-a)*(1-t)+r*f); g=a-mu-r*sig
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Entry identity failed: '+name)
        proofs[name]=True
    zero('entry_remaining_f_plus_J_constant_by_FTC',s.diff(f+J,t).subs({s.diff(f,t):-sig,s.diff(J,t):sig}))
    zero('entry_remaining_f_plus_J_right_endpoint_half',s.Rational(0)+s.Rational(1,2)-s.Rational(1,2))
    zero('actual_entry_K_log_rate',s.diff(K,t).subs(s.diff(f,t),-sig)-g*K)
    zero('actual_entry_KS_right_endpoint',K.subs({t:1,f:0})-KS)
    X=s.Function('native_X')(t,z); en=s.Function('native_half_energy')(t,z)
    zero('actual_native_KX_full_angular_ODE',g*K*X+K*(1-r*(1-sig)*X)-(K-k*K*X))
    zero('actual_native_2K2_full_energy_ODE',4*g*K*K*en+2*K*K*(-s.Rational(1,2)+(2*mu+2*r*sig)*en)
         -(delta*2*K*K*en-K*K))
    C,q=s.symbols('same_exact_pressure_scale q',real=True); absolute=s.Function('native_absolute_pressure')(t,z)
    normalized=-absolute*s.exp(p*q)/C
    zero('actual_native_absolute_pressure_full_normalized_ODE',
         -s.exp(p*q)/C*(C*s.exp(-p*q)*K*K/2)+p*normalized-(p*normalized-K*K/2))
    IA=s.Function('signed_remaining_A')(t); IE=s.Function('signed_remaining_E')(t); IP=s.Function('signed_remaining_P')(t)
    AS=s.Function('same_full_AdS')(z); ES=s.Function('same_full_EdS')(z); PS=s.Function('same_full_PdS')(z)
    Ad=AS*s.exp(k*(1-t))-IA; Ed=ES*s.exp(-delta*(1-t))+IE; Pd=PS*s.exp(-p*(1-t))+IP/2
    rules={s.diff(IA,t):-(K-1)-k*IA,s.diff(IE,t):-(K*K-1)+delta*IE,s.diff(IP,t):-(K*K-1)+p*IP}
    for label,row,rhs,start in (('angular',Ad,K-1-k*Ad,AS),('energy',Ed,delta*Ed-(K*K-1),ES),
                                ('pressure',Pd,p*Pd-(K*K-1)/2,PS)):
        zero('same_full_future_'+label+'_backward_FTC',s.diff(row,t).subs(rules)-rhs)
        zero('same_full_future_'+label+'_power_inlet',row.subs({t:1,IA:0,IE:0,IP:0})-start)
    zero('source_variable_shear_margin',delta-2*g-(2*mu+2*r*sig))
    return dict(identities=proofs,entry_f_equals_half_minus_original_J_by_FTC_and_endpoint=True,
                same_full_moment_backward_FTC_and_native_normalization_verified=True,
                angular_backward_exponential_has_positive_k_sign=True,
                source_kappa_minus2_range='[2*mu,2] for 0<mu<1 and 0<=sigma<=1',
                whole_original_domain=DOMAIN,interval_overlap_used_as_functional_proof=False)


class SourceAST:
    def __init__(self):self.hashes={}; self.bindings={}
    def method(self,stem,name):
        path=HERE/(PREFIX+stem+'.py'); self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        return next(n for n in ast.walk(ast.parse(path.read_text(encoding='utf8'))) if isinstance(n,ast.FunctionDef) and n.name==name)
    def expression(self,stem,name,target,*,wanted=None,augmented=False,count=1):
        fn=self.method(stem,name)
        values=[n.value for n in ast.walk(fn) if
            (not augmented and isinstance(n,ast.Assign) and any(ast.unparse(v)==target for v in n.targets)) or
            (augmented and isinstance(n,ast.AugAssign) and isinstance(n.op,ast.Add) and ast.unparse(n.target)==target)]
        if wanted is not None:
            expected=ast.dump(ast.parse(wanted,mode='eval').body)
            values=[v for v in values if ast.dump(v)==expected]
        if len(values)!=count:raise ValueError('Actual entry source changed: '+stem+'.'+name+'.'+target)
        self.bindings[stem+'.'+name+'.'+target]=True
        return values[0]
    def evaluate(self,node,env):return eval(compile(ast.Expression(node),'<actual entry source>','eval'),env,env)
    def replay(self,stem,name,env,remove_pressure_context=False):
        fn=self.method(stem,name); fn.decorator_list=[]
        if remove_pressure_context:
            values=[n for n in fn.body if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='c'
                    and ast.unparse(n.value)=='pressure_numerator.ctx']
            if len(values)!=1:raise ValueError('Actual pressure context changed')
            fn.body.remove(values[0])
        exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual entry row '+name+'>','exec'),env)
        return env[name]


def entry_power_join_binding():
    """Actual source rows at entry t=1 versus power phase=0, arbitrary axial data."""
    a,mu,z,C,wait,Ts=s.symbols('a mu Z C wait Ts',real=True)
    c=SimpleNamespace(exp=s.exp,expm1=lambda v:s.exp(v)-1,mpf=lambda v:s.Rational(str(v)))
    heat=SimpleNamespace(ctx=c,a=a,mu=mu,S=s.Symbol('same_source_inverse_Rtail'),delta=2*a,k=1-a,prate=1+2*a,
                         pressure_scale=C,steep=SimpleNamespace(wait=wait,Ts=Ts))
    stub=SimpleNamespace(constant=lambda ctx,value,order:value,variable=lambda ctx,value,order:value)
    env=dict(math=math,c=c,IntervalTaylor=stub,axial_derivative=lambda v:s.diff(v,z),
             sigma_jets=lambda ctx,t:[s.Integer(1)]+[s.Integer(0)]*4)
    asts=SourceAST()
    for stem,name in (('collar_Gamma_C4','product_rows'),('collar_stress_C3','shifted_rows'),
        ('collar_stress_C3','collar_stress_rows'),('steep_power_stress_C3','power_angular_rows'),
        ('steep_power_stress_C3','power_shape'),('steep_power_stress_C3','power_defect_rows'),
        ('steep_power_stress_C3','power_stress_rows'),('steep_power_stress_C3','power_pressure_rows'),
        ('steep_entry_stress_C3','entry_shape'),('steep_entry_stress_C3','entry_angular_rows'),
        ('steep_entry_stress_C3','entry_defect_rows'),('steep_entry_stress_C3','entry_stress_rows'),
        ('steep_entry_stress_C3','entry_pressure_rows')):asts.replay(stem,name,env)
    asts.replay('heat_pressure_C4','pressure_y_rows',env,True)
    power={key+'_defect_rows':[s.Function('same_full_'+key+'_Q')(z)] for key in ('angular','energy','pressure')}
    power.update(KQ=s.Symbol('same_source_KQ',positive=True),original_XS=s.Function('same_native_XS')(z))
    ps=env['power_shape'](heat,power,s.Integer(0))
    pd=env['power_defect_rows'](heat,power,s.Integer(0),ps)
    terminal={key:[pd[key][0]] for key in ('angular_defect_rows','energy_defect_rows','pressure_defect_rows')}
    terminal['KS']=ps['K_rows'][0]
    kernels=dict(f=s.Integer(0),angular=s.Integer(0),energy=s.Integer(0),pressure=s.Integer(0))
    es=env['entry_shape'](heat,terminal,s.Integer(1),kernels)
    xr=env['entry_angular_rows'](heat,power['original_XS'],es)
    ed=env['entry_defect_rows'](heat,terminal,s.Integer(1),kernels,es,power['original_XS'])
    left=env['entry_stress_rows'](heat,es,ed,xr,z,-wait-Ts-1)
    right=env['power_stress_rows'](heat,power,s.Integer(0),ps,z)
    checks={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Actual entry/power join failed: '+name)
        checks[name]=True
    for j in range(5):
        zero('shape_y'+str(j),es['K_rows'][j]-ps['K_rows'][j])
        for key in ('angular_defect_rows','energy_defect_rows','pressure_defect_rows'):
            zero(key+'_y'+str(j),ed[key][j]-pd[key][j])
    for label in ('theta','axial'):
        for j in range(4):
            for n in range(4-j):zero(label+'_y'+str(j)+'_Z'+str(n),s.diff(left[label][j]-right[label][j],z,n))
    lp=env['entry_pressure_rows'](heat,ed,es,-wait-Ts-1)
    rp=env['power_pressure_rows'](heat,power,s.Integer(0),ps)
    for j in range(5):
        for n in range(5-j):zero('pressure_y'+str(j)+'_Z'+str(n),s.diff(lp[j]-rp[j],z,n))
    return dict(identities=checks,actual_entry_power_K4_stress_mixed3_pressure_mixed4_AST_join_verified=True,
                arbitrary_actual_terminal_axial_functions_used=True,source_function_equality_not_sample_overlap=True,
                input_hashes=asts.hashes)


def entry_source_bridge(power):
    """Bind original source chart, native ODEs, actual primitive and same datum."""
    asts=SourceAST(); checks={}; hashes=dict(power.hashes)
    if entry_signed_kernels.__globals__['sigma_enclosure'] is not transition_kernels.__globals__['sigma_enclosure']:
        raise ValueError('Entry remaining and original forward primitive callable differ')
    exit_source=power.exit_source
    if not exit_source.proof['original_sigmoid_reflection_gives_half_integral']:
        raise ValueError('Original reflection/half integral required')
    if not exit_source.bridge['original_sigma_between0and1_reflection_and_flat_endpoints_bound']:
        raise ValueError('Actual original sigmoid flat endpoints required')
    asts.expression('steep_waiting_C4','transition_kernels','nextJ',wanted='J+ds*sigma_enclosure(c,v)',count=2)
    asts.expression('steep_waiting_C4','transition_kernels','J',wanted='c.mpf(0)')
    asts.expression('steep_waiting_C4','transition_kernels','J',wanted="c.mpf('.5')")
    asts.expression('steep_entry_stress_C3','entry_signed_kernels','nextf',wanted='f+ds*sigma_enclosure(c,v)')
    asts.expression('steep_entry_stress_C3','entry_signed_kernels','f',wanted='c.mpf(0)')
    asts.expression('steep_entry_stress_C3','entry_signed_kernels','f',wanted="c.mpf('.5')")
    for stem,method,test,assignment in (
        ('steep_waiting_C4','transition_kernels','endpoints(t)==(mp.mpf(1),mp.mpf(1))',"J=c.mpf('.5')"),
        ('steep_entry_stress_C3','entry_signed_kernels','(lo,hi)==(mp.mpf(0),mp.mpf(0))',"f=c.mpf('.5')")):
        node=asts.method(stem,method)
        matches=[n for n in ast.walk(node) if isinstance(n,ast.If)
            and ast.dump(n.test)==ast.dump(ast.parse(test,mode='eval').body)
            and any(ast.dump(v)==ast.dump(ast.parse(assignment).body[0]) for v in n.body)]
        if len(matches)!=1:raise ValueError('Actual primitive half-endpoint guard changed')
        asts.bindings[stem+'.'+method+'.exact_half_endpoint_guard']=True
    for target,expr in {'remaining':'length*(1-c.mpf([i,i+1])/cells)',
        'distance':'length*c.mpf([i,i+1])/cells','exponent':'(mu-a)*remaining+r*fc',
        'Km':'(KS-1)+KS*c.expm1(exponent)','Qm':'(KS**2-1)+KS**2*c.expm1(2*exponent)',
        'fc':"c.mpf([endpoints(f)[0],min(mp.mpf('.5'),endpoints(nextf)[1])])"}.items():
        asts.expression('steep_entry_stress_C3','entry_signed_kernels',target,wanted=expr)
    for target,expr in {'IA':'ds*c.exp(k*distance)*Km','IE':'ds*c.exp(-delta*distance)*Qm',
                         'IP':'ds*c.exp(-p*distance)*Qm'}.items():
        asts.expression('steep_entry_stress_C3','entry_signed_kernels',target,wanted=expr,augmented=True)
    # Guard both original density branches and the entry normalization.
    for target,expr in {'angular':'ds*c.exp(rate*(v-jc))',
                       'pressure':'ds*c.exp(-(1+2*mu)*v-2*rate*jc)/2',
                       'remaining':'ds*c.exp(-2*mu*v-2*rate*jc)'}.items():
        asts.expression('steep_waiting_C4','transition_kernels',target,wanted=expr,augmented=True)
    a,mu,t,Ts,wait,eps,thetaR,Ev2=s.symbols('a mu t Ts wait eps thetaR exact_Ev0_squared',real=True)
    r=1-mu; k=1-a; p=1+2*a; bh=s.Rational(1,2)+a; bp=s.Rational(1,2)+mu
    J=s.Function('actual_J')(t); sig=s.Function('original_sigma')(t)
    IA=s.Function('native_forward_angular')(t); IE=s.Function('native_remaining_energy')(t); IP=s.Function('native_forward_pressure')(t)
    XR,after,PR=s.symbols('actual_XR actual_after_entry actual_PR',real=True)
    c=SimpleNamespace(exp=s.exp,expm1=lambda v:s.exp(v)-1,mpf=lambda v:s.Rational(str(v)))
    obj=SimpleNamespace(thetaR=thetaR,bp=bp,rate=r,mu=mu,outer=SimpleNamespace(flatten=SimpleNamespace(Ev2=Ev2)))
    env=dict(self=obj,c=c,one=s.Integer(1),t=t,J=J,data=dict(XR=XR,after_entry=after,PR=PR),
             kernels=dict(angular=IA,remaining_energy=IE,pressure=IP))
    native={key:asts.evaluate(asts.expression('steep_waiting_C4','steep_in',key),env) for key in ('theta','X','energy','pressure')}
    rules={s.diff(J,t):sig,s.diff(IA,t):s.exp(r*(t-J)),s.diff(IE,t):-s.exp(-2*mu*t-2*r*J),
           s.diff(IP,t):s.exp(-(1+2*mu)*t-2*r*J)/2}
    comparisons=dict(theta=s.diff(native['theta'],t).subs(rules)-(-bp-r*sig)*native['theta'],
        angular=s.diff(native['X'],t).subs(rules)-(1-r*(1-sig)*native['X']),
        energy=s.diff(native['energy'],t).subs(rules)-(-s.Rational(1,2)+(2*mu+2*r*sig)*native['energy']),
        pressure=s.diff(native['pressure'],t).subs(rules)-Ev2*native['theta']**2/2)
    for name,value in comparisons.items():
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Original entry ODE changed: '+name)
        checks['actual_native_'+name+'_source_ODE_replayed']=True
    asts.expression('steep_waiting_C4','steep_in','sig',wanted='sigma_jets(c,t)')
    asts.expression('steep_waiting_C4','steep_in','rates',wanted='[-self.rate*sig[j]*math.factorial(j) for j in range(4)]')
    # Replay actual amplitude definitions and reference, preserving the datum.
    obj.k=k; obj.Ts=Ts
    for target,expr in {'self.Ts':'4*(c.ln(2)-c.ln(self.delta))','self.epsilon':"c.mpf('.001')*self.delta",
        'self.thetaR':'c.exp(-100*self.bp-self.bp*self.outer.Lrel)/2','self.rate':'1-self.mu'}.items():
        asts.expression('steep_waiting_C4','__init__',target,wanted=expr)
    obj.thetaS=asts.evaluate(asts.expression('steep_waiting_C4','__init__','self.thetaS'),env)
    obj.thetaQ=asts.evaluate(asts.expression('steep_waiting_C4','__init__','self.thetaQ'),env)
    obj.thetaT=asts.evaluate(asts.expression('steep_waiting_C4','__init__','self.thetaT'),env)
    theta_base=obj.thetaT*s.exp(-bh*wait)/(1-eps)
    KS=(1-eps)*s.exp(k*(Ts+s.Rational(1,2)))
    q=-wait-Ts-2+t; K=KS*s.exp((mu-a)*(1-t)+r*(s.Rational(1,2)-J))
    if s.simplify(s.expand_power_exp(native['theta']-theta_base*s.exp(-bh*q)*K))!=0:
        raise ArithmeticError('Actual entry velocity normalization changed')
    checks['actual_native_theta_equals_same_reference_B_times_K_verified']=True
    checks['actual_KS_is_original_power_phase0_not_exit_KQ']=True
    for target,expr in {'XS':"(XR+self.infull['angular'])*c.exp(-self.rate/2)",
        'PS':"PR+self.infull['pressure']*(self.outer.flatten.Ev2*self.thetaR**2)",
        'after_entry':'(after_power*c.exp(-2*self.Ts)+decay_integral(c,2,self.Ts))*c.exp(-1-self.mu)'}.items():
        asts.expression('steep_waiting_C4','data',target,wanted=expr)
    # Actual native PS/PQ assignments and the admitted full exit pressure
    # recover the native entry-end datum, including its axial dependence.
    psym,ipfull=s.symbols('same_native_PS same_full_entry_pressure_integral',real=True)
    obj.infull=dict(pressure=ipfull)
    psnative=asts.evaluate(asts.expression('steep_waiting_C4','data','PS'),dict(env,PR=PR))
    atentry=native['pressure'].subs({t:1,J:s.Rational(1,2),IP:ipfull})
    if s.simplify(atentry-psnative)!=0:raise ArithmeticError('Actual entry absolute-pressure PS endpoint changed')
    decay=lambda ctx,rate,length:(1-s.exp(-rate*length))/rate
    pqnative=asts.evaluate(asts.expression('steep_waiting_C4','data','PQ'),dict(env,PS=psym,decay_integral=decay))
    C=Ev2*theta_base**2; qQ=-wait-1; KQ=(1-eps)*s.exp(k/2)
    # This normalized PQ is identified by the accepted absolute-pressure
    # history bridge; no free datum or independently selected cap is used.
    fullPQ=-pqnative/(C*s.exp(-p*qQ))
    pressureenv=dict(c=c,IntervalTaylor=SimpleNamespace(constant=lambda ctx,value,order:value))
    powerpressure=asts.replay('steep_power_stress_C3','power_pressure_rows',pressureenv)
    he=SimpleNamespace(ctx=c,prate=p,pressure_scale=C,steep=SimpleNamespace(wait=wait))
    leftpressure=powerpressure(he,dict(KQ=KQ,pressure_defect_rows=[fullPQ-1/(2*p)]),s.Integer(0),dict(remaining=Ts))
    if s.simplify(s.expand_power_exp(leftpressure[0]-psym))!=0:
        raise ArithmeticError('Same full exit absolute-pressure datum does not recover actual PS')
    powernative=asts.evaluate(asts.expression('steep_waiting_C4','steep_power','pressure'),
        dict(env,t=s.Integer(0),data=dict(PS=psym),decay_integral=decay))
    if s.simplify(powernative-psym)!=0:raise ArithmeticError('Actual native power inlet pressure changed')
    if not exit_source.bridge['actual_original_absolute_pressure_datum_consumed']:
        raise ValueError('Actual normalized full exit PQ datum missing')
    checks.update(actual_native_entry_pressure_right_endpoint_equals_original_PS=True,
        actual_full_exit_normalized_PQ_datum_consumed=True,
        actual_power_pressure_backward_full_PQ_recovers_native_PS_AST_verified=True,
        actual_native_entry_pressure_endpoint_equals_same_full_power_pressure_verified=True)
    asts.expression('power_angular_C4','_packet','mixed',
        wanted='flatten_mixed(c,self.mu,one,[c.mpf(0)]+rates,theta,X,energy,pressure,self.flatten.Ev2)')
    fn=asts.method('steep_waiting_C4','packet')
    route=[n for n in ast.walk(fn) if isinstance(n,ast.Call) and ast.unparse(n.func)=='self.outer._packet']
    if len(route)!=1 or [ast.unparse(v) for v in route[0].args[:5]]!=['Z','theta','X','energy','pressure']:
        raise ValueError('Original entry packet field route changed')
    asts.bindings['steep_waiting_C4.packet_to_actual_outer_packet']=True
    path=HERE/(PREFIX+'global_physical_assembly.py'); tree=ast.parse(path.read_text(encoding='utf8')); radiusenv={}
    for node in tree.body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ('MICRO','POST','PULSE'):radiusenv[target.id]=ast.literal_eval(node.value)
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='radius')
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual entry radius>','exec'),radiusenv)
    origin,Lrel=s.symbols('logRp Lrel',real=True)
    assembly=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify),logRp=origin,params=SimpleNamespace(mu=mu))
    provider=SimpleNamespace(outer=SimpleNamespace(Lrel=Lrel),Ts=Ts,wait=wait)
    actualq=radiusenv['radius'](assembly,'steep_entry',t,{},provider)[0]-radiusenv['radius'](assembly,'waiting',s.Integer(1),{},provider)[0]
    if s.simplify(actualq-q)!=0:raise ArithmeticError('Original entry radius/reference changed')
    hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    name=PREFIX+'steep_waiting_C4_check.json'; receipt=json.loads((HERE/name).read_bytes())
    if not receipt['all_passed']:raise ValueError('Original entry native C4 joins required')
    for name0,digest in receipt['input_hashes'].items():
        if hashlib.sha256((HERE/name0).read_bytes()).hexdigest()!=digest:raise ValueError('Original entry source changed: '+name0)
    hashes.update(receipt['input_hashes']); hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    flags=('angular_entry_exact_theta','angular_entry_exact_X','angular_entry_exact_energy','angular_entry_exact_pressure',
        'entry_power_exact_theta','entry_power_exact_X','entry_power_exact_energy','entry_power_exact_pressure',
        'original_entry_X_history','original_entry_backward_energy_units','new_data_XS_history','new_data_PS_history',
        'new_data_after_entry','same_original_postangular_future_decomposition','new_exact_preheat_Gamma_source',
        'new_exact_steep_power_length','new_actual_epsilon_relation','new_angular_terminal_velocity_scale','new_thetaS_velocity_history')
    consumed={}
    for flag in flags:
        if not receipt['functional_production_source_identities'][flag]:raise ValueError('Actual entry field/history identity missing: '+flag)
        consumed[flag]=True
    for flag,value in receipt['functional_production_source_identities'].items():
        if flag.startswith(('entry_power_','angular_entry_')) and ('_axial' in flag or '_y' in flag):
            if not value:raise ValueError('Actual entry high-order native field join missing: '+flag)
            consumed[flag]=True
    for flag in ('actual_source_XS_correlation_is_same_full_angular_future',
                 'actual_absolute_pressure_datum_and_exact_Ev0_units_consumed',
                 'actual_meridional_source_and_terminal_zero_primitive_route_consumed'):
        if not power.bridge[flag]:raise ValueError('Same actual power full future missing: '+flag)
        checks['consumed_'+flag]=True
    for method,target,expr in (
        ('terminal',"self.cache[key]['KS']","point['steep_power_shape']['K_rows'][0][0]"),
        ('entry_defect_rows','Ad',"shape['K_rows'][0]*X0-one/heat.k"),
        ('entry_defect_rows','Ed',"terminal['energy_defect_rows'][0]*c.exp(-heat.delta*ell)+one*kernels['energy']"),
        ('entry_defect_rows','Pd',"terminal['pressure_defect_rows'][0]*c.exp(-heat.prate*ell)+one*kernels['pressure']/2"),
        ('steep_in','mixed[P]',"{'y'+str(j)+'_Z'+str(n):jet[n]*math.factorial(n) for j,jet in enumerate(pressure) for n in range(5-j)}"),
        ('steep_in','fields[P]','[jet.truncate(4-j) for j,jet in enumerate(pressure)]')):
        asts.expression('steep_entry_stress_C3',method,target,wanted=expr)
    join=entry_power_join_binding(); hashes.update(join['input_hashes']); hashes.update(asts.hashes)
    return dict(identities=checks,actual_entry_source_AST_bindings=asts.bindings,
        original_native_C4_field_and_history_identities_consumed=consumed,
        actual_forward_J_and_backward_f_same_sigma_FTC_half_endpoint_bound=True,
        actual_entry_f_equals_half_minus_original_J_verified=True,
        actual_original_entry_production_radius_reference_verified=True,
        actual_complete_power_inlet_moments_and_absolute_datum_consumed=True,
        native_X_KX_equals_same_full_future_by_original_ODE_and_right_endpoint=True,
        native_energy_2K2_equals_same_full_future_by_original_ODE_and_right_endpoint=True,
        native_absolute_pressure_equals_same_full_future_by_original_ODE_and_right_endpoint=True,
        native_angular_entry_left_field_pressure_mixed4_join_consumed=True,
        upstream_angular_stress_companion_constructed=False,
        actual_entry_power_formula_join=join,original_coefficients_or_velocity_changed=False,input_hashes=hashes)


class CompliantSteepEntryStressC3:
    @source_precision
    def __init__(self,cells=128):
        self.power_source=CompliantSteepPowerStressC3(); self.heat=self.power_source.heat; self.steep=self.heat.steep
        self.ctx=self.heat.ctx; self.family,self.source=self.power_source.family,self.power_source.source
        self.cells=cells; self.cache={}; self.kernel_cache={}; self.hashes=dict(self.power_source.hashes)
        name=PREFIX+'steep_power_stress_C3_check.json'; receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or not receipt['actual_original_steep_power_similarity_stress_recovered']:
            raise ValueError('Accepted original power complete future required')
        for path,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Entry prerequisite changed: '+path)
        self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        if endpoints(self.heat.mu)[0]<=0 or endpoints(self.heat.mu)[1]>=1:raise ValueError('Original 0<mu<1 required')
        self.proof=entry_transport_identities(); self.bridge=entry_source_bridge(self.power_source)
        self.hashes.update(self.bridge['input_hashes']); self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def terminal(self,Z):
        key=tuple(endpoints(self.ctx.mpf(Z)))
        if key not in self.cache:
            point=self.power_source.steep_power(Z,0)
            self.cache[key]={label:[jet] for label,jet in point['steep_power_future_defect_zeroth_Taylor'].items()}
            self.cache[key]['KS']=point['steep_power_shape']['K_rows'][0][0]
        return self.cache[key]

    @source_precision
    def kernels(self,Z,t):
        key=tuple(endpoints(self.ctx.mpf(t)))
        if key not in self.kernel_cache:
            self.kernel_cache[key]=entry_signed_kernels(self.ctx,t,self.terminal(Z)['KS'],self.heat.mu,self.heat.delta,self.cells)
        return self.kernel_cache[key]

    @source_precision
    def steep_in(self,Z,t):
        c=self.ctx; Z=c.mpf(Z); t=c.mpf(t); self.steep._phase(t)
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1:raise ValueError('Original Z in[-1,1] required')
        point=dict(self.steep.steep_in(Z,t)); terminal=self.terminal(Z); kernels=self.kernels(Z,t)
        shape=entry_shape(self.heat,terminal,t,kernels); one=IntervalTaylor.constant(c,1,5)
        Xrows=entry_angular_rows(self.heat,point['angular_Taylor'],shape)
        defects=entry_defect_rows(self.heat,terminal,t,kernels,shape,Xrows[0])
        q=-self.steep.wait-self.steep.Ts-2+t; rows=entry_stress_rows(self.heat,shape,defects,Xrows,Z,q)
        pressure=entry_pressure_rows(self.heat,defects,shape,q)
        for label in ('angular','energy'):
            point['original_forward_'+label+'_Taylor']=point[label+'_Taylor']
            point['original_forward_'+label+'_y_derivative_Taylor']=point[label+'_y_derivative_Taylor']
        point['original_forward_pressure_over_Pstar_squared_Taylor']=point['pressure_over_Pstar_squared_Taylor']
        E=[one/self.heat.delta+defects['energy_defect_rows'][0]]+defects['energy_defect_rows'][1:]
        energy=quotient_rows(E,[jet*2 for jet in product_rows(shape['K_rows'],shape['K_rows'])])
        mixed=dict(point['physical_mixed_derivatives_total_order_le4']); fields=dict(point['physical_velocity_and_pressure_y_derivative_Taylor'])
        point['original_forward_pressure_mixed_bounds']=mixed[P]
        mixed[P]={'y'+str(j)+'_Z'+str(n):jet[n]*math.factorial(n) for j,jet in enumerate(pressure) for n in range(5-j)}
        fields[P]=[jet.truncate(4-j) for j,jet in enumerate(pressure)]
        margin=2*self.heat.mu+2*(1-self.heat.mu)*shape['original_sigma_jets'][0]
        point.update(angular_Taylor=Xrows[0],angular_y_derivative_Taylor=[jet.truncate(5-j) for j,jet in enumerate(Xrows)],
            energy_Taylor=energy[0],energy_y_derivative_Taylor=[jet.truncate(5-j) for j,jet in enumerate(energy)],
            pressure_over_Pstar_squared_Taylor=pressure[0],pressure_y_derivative_axial5_Taylor=pressure,
            physical_mixed_derivatives_total_order_le4=mixed,physical_velocity_and_pressure_y_derivative_Taylor=fields,
            steep_entry_phase=t,steep_entry_offset_from_Rtail=q,steep_entry_shape=shape,
            steep_entry_signed_remaining_kernels=kernels,
            steep_entry_similarity_stress_mixed3_factored={label:{'y'+str(j)+'_Z'+str(n):rows[label][j][n]*math.factorial(n)
                for j in range(4) for n in range(4-j)} for label in ('theta','axial')},
            steep_entry_similarity_stress_y_derivative_Taylor={label:[jet.truncate(3-j) for j,jet in enumerate(rows[label])]
                for label in ('theta','axial','theta_inertial','theta_shear')},
            steep_entry_future_defect_zeroth_Taylor={label:values[0] for label,values in defects.items()},
            steep_entry_shear_strength_kappa_minus2=margin,
            steep_entry_angular_shear_negative_by_exact_source=True,
            steep_entry_meridional_moments_and_velocities={label:one*0 for label in ('Mz','Mtheta_z','Ur','Uz')},
            actual_original_steep_entry_similarity_stress_recovered=True,
            actual_steep_entry_absolute_pressure_same_source_mixed4_available=True,
            steep_entry_power_stress_mixed3_join_verified=True,steep_entry_power_pressure_mixed4_join_verified=True,
            native_angular_entry_left_field_pressure_join_verified=True,upstream_angular_stress_companion_constructed=False,
            same_complete_waiting_collar_Gamma_future_used=True,original_source_histories_retained=True,
            source_defined_positive_stress_factors=dict(theta='sqrt(R/2)*B',axial='sqrt(R/2)*B^2',
                B='exact_Ev0*theta_base*exp(-bh*q)',R='Rtail*exp(q)',q='-wait-Ts-2+t'),
            steep_entry_physical_decomposition_constructed=False,steep_entry_cone_certified=False,
            global_admissible_stress_lift_constructed=False,physical_energy_integral_certified=False,temporal_recursion=False)
        return point

    def report(self):
        keys=('Z','steep_entry_phase','steep_entry_offset_from_Rtail','steep_entry_shape','steep_entry_signed_remaining_kernels',
            'steep_entry_similarity_stress_mixed3_factored','steep_entry_similarity_stress_y_derivative_Taylor',
            'steep_entry_future_defect_zeroth_Taylor','pressure_y_derivative_axial5_Taylor',
            'steep_entry_shear_strength_kappa_minus2','steep_entry_meridional_moments_and_velocities','source_defined_positive_stress_factors')
        def summary(z,t):
            point=self.steep_in(z,t); return {key:point[key] for key in keys}
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            scope='Original sigmoid entry t[0,1],Z[-1,1]; same full future similarity stress mixed3 and absolute pressure mixed4',
            samples=[summary(z,t) for z,t in (('0','0'),('.5','.5'),('-.5','1'))],whole_steep_entry=summary([-1,1],[0,1]),
            transport_identities=self.proof,steep_entry_source_bridge=self.bridge,
            actual_original_steep_entry_similarity_stress_recovered=True,
            actual_steep_entry_absolute_pressure_same_source_mixed4_available=True,
            steep_entry_power_stress_mixed3_join_verified=True,steep_entry_power_pressure_mixed4_join_verified=True,
            native_angular_entry_left_field_pressure_join_verified=True,upstream_angular_stress_companion_constructed=False,
            steep_entry_physical_decomposition_constructed=False,steep_entry_cone_certified=False,
            global_admissible_stress_lift_constructed=False,physical_energy_integral_certified=False,
            temporal_recursion=False,input_hashes=self.hashes)


def run():
    result=CompliantSteepEntryStressC3().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original steep entry same full-future stress mixed3, absolute pressure mixed4 and power joins generated',flush=True)
    return result


if __name__=='__main__':run()
