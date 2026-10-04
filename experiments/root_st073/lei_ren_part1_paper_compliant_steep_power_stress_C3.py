"""Original Ts-long steep power: SAME complete-future stress and pressure.

Backward resonant angular and quadratic moment modes are integrated exactly.
Unit baselines and pressure stationary modes cancel before enclosure. Native
velocities, original phase, datum and forward histories remain unchanged.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_steep_exit_stress_C3 import (
    CompliantSteepExitStressC3,source_precision,quotient_rows,exit_signed_kernels)
from lei_ren_part1_paper_compliant_steep_waiting_C4 import transition_kernels
from lei_ren_part1_paper_compliant_collar_stress_C3 import axial_derivative
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_pulse_physical_bounds import P
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
DOMAIN=dict(phase=[0,1],Z=[-1,1],original_offset='y=Ts*phase',
            remaining='ell=Ts*(1-phase)',offset='q=-wait-1-ell')


def power_angular_rows(heat,terminal,phase):
    """Same full-future A/K, using its source-functional inlet correlation."""
    one=IntervalTaylor.constant(heat.ctx,1,5)
    X=terminal['original_XS']+one*(heat.steep.Ts*phase)
    return [X,one]+[one*0]*3


def power_shape(heat,terminal,phase):
    c=heat.ctx; one=IntervalTaylor.constant(c,1,5)
    ell=heat.steep.Ts*(1-phase); K=terminal['KQ']*c.exp(heat.k*ell)
    rows=[one*K*(-heat.k)**j for j in range(5)]
    square=product_rows(rows,rows)
    return dict(K_rows=rows,K_defect_rows=[rows[0]-one]+rows[1:],
                K_squared_defect_rows=[square[0]-one]+square[1:],
                log_K_rate_rows=[-heat.k]+[c.mpf(0)]*3,remaining=ell)


def power_defect_rows(heat,terminal,phase,shape):
    """Exact backward full-moment FTC, with stable unit-baseline defects."""
    c=heat.ctx; one=IntervalTaylor.constant(c,1,5); ell=shape['remaining']; KQ=terminal['KQ']
    AdQ=terminal['angular_defect_rows'][0]; EdQ=terminal['energy_defect_rows'][0]
    PdQ=terminal['pressure_defect_rows'][0]
    X=power_angular_rows(heat,terminal,phase)[0]
    Ad=shape['K_rows'][0]*X-one/heat.k
    Ed=EdQ*c.exp(-heat.delta*ell)+one*(c.expm1(-heat.delta*ell)/heat.delta
         +KQ**2*c.exp(-heat.delta*ell)*c.expm1(2*ell)/2)
    Pd=PdQ*c.exp(-heat.prate*ell)+one*(c.expm1(-heat.prate*ell)/(2*heat.prate)
         +KQ**2*c.exp(-heat.prate*ell)*c.expm1(3*ell)/6)
    A=[Ad]; E=[Ed]; Pr=[Pd]; Km=shape['K_defect_rows']; Qm=shape['K_squared_defect_rows']
    for j in range(4):
        A.append(Km[j]-A[j]*heat.k)
        E.append(E[j]*heat.delta-Qm[j])
        Pr.append(Pr[j]*heat.prate-Qm[j]/2)
    return dict(angular_defect_rows=A,energy_defect_rows=E,pressure_defect_rows=Pr,
                K_defect_rows=Km,K_squared_defect_rows=Qm)


def power_stress_rows(heat,terminal,phase,shape,Z):
    """Closed modes of (3.16)-(3.18), including source factor derivatives."""
    c=heat.ctx; z=IntervalTaylor.variable(c,Z,5); one=IntervalTaylor.constant(c,1,5)
    ell=shape['remaining']; K=shape['K_rows'][0]; KQ=terminal['KQ']
    k=heat.k; delta=heat.delta; p=heat.prate; b=(1-delta)/2; L=1-z*z*delta; d=1-z*z
    AdQ=terminal['angular_defect_rows'][0]; EdQ=terminal['energy_defect_rows'][0]
    PdQ=terminal['pressure_defect_rows'][0]
    X=power_angular_rows(heat,terminal,phase)[0]
    It=K*(X*k-z*axial_derivative(X)*b-one)/L
    q=-heat.steep.wait-1-ell
    St=K*(-4*heat.S*c.exp(-q))
    # E=K^2/2+exp(-delta*ell)*(E_Q-KQ^2/2), and
    # P=K^2/6+exp(-p*ell)*(P_Q-KQ^2/6). Remove the unit
    # pair z*(exp(-delta*ell)-exp(-p*ell)) before enclosure.
    pure=-z*K*K*k/(3*L)
    Ne=z*EdQ*delta-z*(delta*KQ**2/2)-d*axial_derivative(EdQ)/2
    Np=-z*PdQ*(2*p)+z*(p*KQ**2/3)+d*axial_derivative(PdQ)
    Je=Ne/L*c.exp(-delta*ell); Jp=Np/L*c.exp(-p*ell)
    inertial=[]; shear=[]; axial=[]
    for j in range(4):
        inertial.append(It*(-1)**j+(K*k/L)*(j*(-1)**(j-1) if j else 0))
        shear.append(St*(-2)**j)
        unit_difference=(-c.expm1(-ell) if j%2==0 else 1+c.exp(-ell))
        unit=z/L*c.exp(-delta*ell)*(-c.mpf('.5'))**j*unit_difference
        axial.append(pure*(-c.mpf('2.5'))**j+Je*(-c.mpf('.5'))**j+Jp*c.mpf('.5')**j+unit)
    return dict(theta=[inertial[j]+shear[j] for j in range(4)],axial=axial,
                theta_inertial=inertial,theta_shear=shear,
                normalized_coefficient_rows=dict(theta=[It+St],axial=[axial[0]]))


def power_pressure_rows(heat,terminal,phase,shape):
    """Full absolute pressure, retaining the exit inlet datum exactly."""
    c=heat.ctx; one=IntervalTaylor.constant(c,1,5); ell=shape['remaining']; KQ=terminal['KQ']
    Pq=one/(2*heat.prate)+terminal['pressure_defect_rows'][0]
    factor=heat.pressure_scale*c.exp(heat.prate*(heat.steep.wait+1))
    rows=[-(Pq+one*(KQ**2*c.expm1(3*ell)/6))*factor]
    rows.extend(one*(KQ**2*c.exp(3*ell)*(-3)**(j-1)/2)*factor for j in range(1,5))
    return rows


def power_transport_identities():
    """Exact resonant/full-quadratic FTC and fully factored stress modes."""
    a,q,z,qQ=s.symbols('a q Z qQ',real=True); delta=2*a; k=1-a; p=1+delta; ell=qQ-q
    KQ=s.symbols('original_KQ',positive=True); K=KQ*s.exp(k*ell)
    AdQ=s.Function('full_AdQ')(z); EdQ=s.Function('full_EdQ')(z); PdQ=s.Function('full_PdQ')(z)
    A=s.exp(k*ell)*(1/k+AdQ-KQ*ell)
    E=K*K/2+s.exp(-delta*ell)*(1/delta+EdQ-KQ*KQ/2)
    Pr=K*K/6+s.exp(-p*ell)*(1/(2*p)+PdQ-KQ*KQ/6)
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Steep power identity failed: '+name)
        proofs[name]=True
    for label,row,rhs,start in (
        ('angular',A,K-k*A,1/k+AdQ),('energy',E,delta*E-K*K,1/delta+EdQ),
        ('pressure',Pr,p*Pr-K*K/2,1/(2*p)+PdQ)):
        zero('same_full_future_'+label+'_backward_FTC',s.diff(row,q)-rhs)
        zero('same_full_future_'+label+'_exit_endpoint',row.subs(q,qQ)-start)
    Edstable=EdQ*s.exp(-delta*ell)+(s.exp(-delta*ell)-1)/delta+KQ*KQ*s.exp(-delta*ell)*(s.exp(2*ell)-1)/2
    Pdstable=PdQ*s.exp(-p*ell)+(s.exp(-p*ell)-1)/(2*p)+KQ*KQ*s.exp(-p*ell)*(s.exp(3*ell)-1)/6
    zero('stable_energy_defect_with_exact_rate_two',E-1/delta-Edstable)
    zero('stable_pressure_defect_with_exact_rate_three',Pr-1/(2*p)-Pdstable)
    L=1-delta*z*z; d=1-z*z; b=(1-delta)/2
    X=A/K; It=(k*A-b*z*s.diff(A,z)-K)/L
    zero('actual_resonant_angular_X_has_slope_one',s.diff(X,q)-1)
    zero('normalized_angular_stress_resonance',It-K*(k*X-b*z*s.diff(X,z)-1)/L)
    St=-4*s.Symbol('exact_inverse_Rtail',positive=True)*s.exp(-q)*K
    zero('actual_power_shear_rate_and_strength',2*s.diff(K,q)/K-delta+2)
    Cz=(delta*z*E-d*s.diff(E,z)/2-2*p*z*Pr+d*s.diff(Pr,z))/L
    pure=-k*z*K*K/(3*L)
    Ne=delta*z*EdQ-delta*z*KQ*KQ/2-d*s.diff(EdQ,z)/2
    Np=-2*p*z*PdQ+p*z*KQ*KQ/3+d*s.diff(PdQ,z)
    Je=Ne*s.exp(-delta*ell)/L; Jp=Np*s.exp(-p*ell)/L
    zero('whole_axial_stress_unit_baselines_cancelled',
         Cz-pure-Je-Jp-z*s.exp(-delta*ell)*(1-s.exp(-ell))/L)
    theta_factor=s.exp(-a*q); axial_factor=s.exp(-(s.Rational(1,2)+delta)*q)
    for j in range(4):
        expected=It*(-1)**j+K*k/L*(j*(-1)**(j-1) if j else 0)
        zero('angular_inertial_full_factor_y'+str(j),s.diff(theta_factor*It,q,j)/theta_factor-expected)
        zero('angular_shear_full_factor_y'+str(j),s.diff(theta_factor*St,q,j)/theta_factor-St*(-2)**j)
        unit=z*s.exp(-delta*ell)/L*((-s.Rational(1,2))**j-s.exp(-ell)*s.Rational(1,2)**j)
        expected=pure*(-s.Rational(5,2))**j+Je*(-s.Rational(1,2))**j+Jp*s.Rational(1,2)**j+unit
        zero('axial_three_modes_and_units_full_factor_y'+str(j),s.diff(axial_factor*Cz,q,j)/axial_factor-expected)
    absolute=-s.exp(-p*q)*Pr; factor=s.exp(-p*qQ)
    zero('absolute_pressure_datum_and_current_density_split',
         absolute+factor*(1/(2*p)+PdQ+KQ*KQ*(s.exp(3*ell)-1)/6))
    for j in range(1,5):
        zero('absolute_pressure_original_density_y'+str(j),
             s.diff(absolute,q,j)-factor*KQ*KQ*s.exp(3*ell)*(-3)**(j-1)/2)
    return dict(identities=proofs,original_full_future_backward_FTC_verified=True,
                exact_resonant_Kq_over_K_equals_minus_k=True,exact_kappa_minus2_equals_two=True,
                rate_delta_plus_two_k_equals_two=True,rate_p_plus_two_k_equals_three=True,
                actual_source_pressure_datum_retained=True,whole_original_domain=DOMAIN,
                sampled_overlap_used_as_functional_proof=False)


def power_exit_join_binding():
    """Replay actual row ASTs on arbitrary axial terminal data at phase1."""
    a,eps,z,C,wait,Ts=s.symbols('a eps Z C wait Ts',real=True)
    c=SimpleNamespace(exp=s.exp,expm1=lambda value:s.exp(value)-1,mpf=lambda value:s.Rational(str(value)))
    heat=SimpleNamespace(ctx=c,a=a,eps=eps,S=s.Symbol('exact_inverse_Rtail'),delta=2*a,k=1-a,
                         prate=1+2*a,pressure_scale=C,steep=SimpleNamespace(wait=wait,Ts=Ts))
    stub=SimpleNamespace(constant=lambda ctx,value,order:value,variable=lambda ctx,value,order:value)
    env=dict(math=math,IntervalTaylor=stub,axial_derivative=lambda value:s.diff(value,z),
             sigma_jets=lambda ctx,t:[s.Integer(0)]*5)
    hashes={}
    def replay(stem,name,remove_context=False):
        path=HERE/(PREFIX+stem+'.py'); hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        fn=next(n for n in ast.walk(ast.parse(path.read_text(encoding='utf8'))) if isinstance(n,ast.FunctionDef) and n.name==name)
        fn.decorator_list=[]
        if remove_context:
            assignments=[n for n in fn.body if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='c'
                         and ast.unparse(n.value)=='pressure_numerator.ctx']
            if len(assignments)!=1:raise ValueError('Pressure context lookup changed')
            fn.body.remove(assignments[0]); env['c']=c
        exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual power/exit '+name+'>','exec'),env)
        return env[name]
    replay('collar_Gamma_C4','product_rows'); replay('collar_stress_C3','shifted_rows')
    generic=replay('collar_stress_C3','collar_stress_rows')
    replay('heat_pressure_C4','pressure_y_rows',True)
    replay('waiting_stress_C3','waiting_defect_rows')
    exit_shape_fn=replay('steep_exit_stress_C3','exit_shape')
    exit_defect_fn=replay('steep_exit_stress_C3','exit_defect_rows')
    exit_pressure_fn=replay('steep_exit_stress_C3','exit_pressure_rows')
    replay('steep_power_stress_C3','power_angular_rows')
    shape_fn=replay('steep_power_stress_C3','power_shape')
    defect_fn=replay('steep_power_stress_C3','power_defect_rows')
    stress_fn=replay('steep_power_stress_C3','power_stress_rows')
    pressure_fn=replay('steep_power_stress_C3','power_pressure_rows')
    kernels=dict(f=s.Rational(1,2),angular=s.Symbol('signed_exit_angular'),
                 energy=s.Symbol('signed_exit_energy'),pressure=s.Symbol('signed_exit_pressure'))
    right_shape=exit_shape_fn(heat,s.Integer(0),kernels)
    base={label+'_defect_rows':[s.Function(label+'_full_collar')(z)] for label in ('angular','energy','pressure')}
    right_defects=exit_defect_fn(heat,base,s.Integer(0),kernels,right_shape)
    terminal={label:[right_defects[label][0]] for label in ('angular_defect_rows','energy_defect_rows','pressure_defect_rows')}
    terminal['KQ']=right_shape['K_rows'][0]
    terminal['original_XS']=(1/heat.k+terminal['angular_defect_rows'][0])/terminal['KQ']-Ts
    left_shape=shape_fn(heat,terminal,s.Integer(1)); left_defects=defect_fn(heat,terminal,s.Integer(1),left_shape)
    left=stress_fn(heat,terminal,s.Integer(1),left_shape,z)
    right=generic(heat,right_shape,right_defects,z,-wait-1)
    checks={}
    for j in range(5):
        if s.simplify(left_shape['K_rows'][j]-right_shape['K_rows'][j])!=0:raise ArithmeticError('Actual power/exit K join failed')
        checks['shape_y'+str(j)]=True
        for label in ('angular','energy','pressure'):
            if j==0 and s.simplify(left_defects[label+'_defect_rows'][0]-terminal[label+'_defect_rows'][0])!=0:
                raise ArithmeticError('Actual power/exit full moment endpoint failed')
    for label in ('theta','axial'):
        for j in range(4):
            for n in range(4-j):
                if s.simplify(s.diff(left[label][j]-right[label][j],z,n))!=0:raise ArithmeticError('Actual power/exit stress join failed')
                checks[label+'_y'+str(j)+'_Z'+str(n)]=True
    lp=pressure_fn(heat,terminal,s.Integer(1),left_shape)
    rp=exit_pressure_fn(heat,right_defects,right_shape,-wait-1)
    for j in range(5):
        for n in range(5-j):
            if s.simplify(s.diff(lp[j]-rp[j],z,n))!=0:raise ArithmeticError('Actual power/exit pressure join failed')
            checks['pressure_y'+str(j)+'_Z'+str(n)]=True
    return dict(identities=checks,actual_power_exit_K_mixed4_join_verified=True,
                actual_power_exit_stress_mixed3_AST_join_verified=True,
                actual_power_exit_pressure_mixed4_AST_join_verified=True,
                arbitrary_full_terminal_axial_functions_used=True,input_hashes=hashes)


def power_angular_endpoint_binding(exit_source):
    """Bind the full companion AQ/KQ to original XQ through common history.

    Actual native XS/XQ/XT/waiting/exit expressions and actual backward
    defect functions are replayed. The signed angular integral is related
    to the native forward integral by its actual integrand, not by overlap.
    """
    collar=exit_source.waiting_source.collar_source; history=collar.history
    if not history['complete_terminal_moment_history_bridge_verified']:
        raise ValueError('Complete actual terminal history required')
    flags=('actual_C4_angular_terminal_is_original_XR',
           'same_entry_exit_angular_integrals_follow_original_sigma_J_weights',
           'actual_angular_defect_zero_from_original_waiting_and_repair_equations')
    for flag in flags:
        if not history['identities'][flag]:raise ValueError('Actual angular inlet history missing: '+flag)
    if not exit_source.proof['original_sigmoid_reflection_gives_half_integral']:
        raise ValueError('Same source sigmoid half integral required for f=J-v+1/2')
    if not exit_source.bridge['original_sigma_between0and1_reflection_and_flat_endpoints_bound']:
        raise ValueError('Actual source sigmoid reflection/flat endpoints required')
    if exit_signed_kernels.__globals__['sigma_enclosure'] is not transition_kernels.__globals__['sigma_enclosure']:
        raise ValueError('Signed remaining and native forward primitive sigmoid callables differ')
    if not collar.bridge['identities']['actual_zero_angular_history_defect_transfers_to_all_collar_offsets']:
        raise ValueError('Original full collar angular moment must equal original zero-defect history')
    k=s.Symbol('actual_k',positive=True); z,v,J,Ts,wait,mu,eps=s.symbols('Z v original_J Ts wait mu eps',real=True)
    c=SimpleNamespace(exp=s.exp,expm1=lambda value:s.exp(value)-1,mpf=lambda value:s.Rational(str(value)))
    hashes={}; bindings={}
    def expression(stem,method,target,augmented=False,wanted=None):
        path=HERE/(PREFIX+stem+'.py'); hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        fn=next(n for n in ast.walk(ast.parse(path.read_text(encoding='utf8'))) if isinstance(n,ast.FunctionDef) and n.name==method)
        values=[n.value for n in ast.walk(fn) if
                (not augmented and isinstance(n,ast.Assign) and any(ast.unparse(t)==target for t in n.targets)) or
                (augmented and isinstance(n,ast.AugAssign) and isinstance(n.op,ast.Add) and ast.unparse(n.target)==target)]
        if wanted is not None:
            expected=ast.dump(ast.parse(wanted,mode='eval').body)
            values=[node for node in values if ast.dump(node)==expected]
        if len(values)!=1:raise ValueError('Actual angular endpoint expression not unique: '+stem+'.'+method+'.'+target)
        bindings[stem+'.'+method+'.'+target]=True
        return values[0]
    def evaluate(node,environment):return eval(compile(ast.Expression(node),'<actual angular endpoint expression>','eval'),{},environment)
    def source_function(stem,name):
        path=HERE/(PREFIX+stem+'.py'); hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        return next(n for n in ast.walk(ast.parse(path.read_text(encoding='utf8'))) if isinstance(n,ast.FunctionDef) and n.name==name)
    def assignments(fn,target,wanted,count):
        name=getattr(fn,'name','endpoint_J1')
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t)==target for t in n.targets)]
        expected=ast.dump(ast.parse(wanted,mode='eval').body)
        matches=[n for n in values if ast.dump(n)==expected]
        if len(matches)!=count:raise ValueError('Actual primitive source binding failed: '+name+'.'+target)
        bindings[name+'.'+target+'.'+wanted]=True
        return matches[0]
    forward_fn=source_function('steep_waiting_C4','transition_kernels')
    remaining_fn=source_function('steep_exit_stress_C3','exit_signed_kernels')
    forward_step=assignments(forward_fn,'nextJ','J+ds*sigma_enclosure(c,v)',2)
    assignments(forward_fn,'J','c.mpf(0)',1)
    endpoint_nodes=[n for n in ast.walk(forward_fn) if isinstance(n,ast.If) and
                    ast.dump(n.test)==ast.dump(ast.parse('endpoints(t)==(mp.mpf(1),mp.mpf(1))',mode='eval').body)]
    if len(endpoint_nodes)!=1:raise ValueError('Actual primitive J(1) source guard missing')
    J1=evaluate(assignments(endpoint_nodes[0],'J',"c.mpf('.5')",1),dict(c=c))
    remaining_step=assignments(remaining_fn,'nextf','f+ds*(1-sigma_enclosure(c,v))',1)
    assignments(remaining_fn,'length','1-t',1)
    f1=evaluate(assignments(remaining_fn,'f','c.mpf(0)',1),dict(c=c))
    if sum(isinstance(n,ast.If) and ast.dump(n.test)==ast.dump(ast.parse('endpoints(length)[1]>0',mode='eval').body)
           for n in ast.walk(remaining_fn))!=1:raise ValueError('Actual remaining primitive zero-length loop guard missing')
    point_fn=source_function('axial_pulse_field','point')
    odds_node=assignments(point_fn,'odds','-1/y**2+1/(1-y)**2',1)
    exp_node=assignments(point_fn,'e','c.exp(odds)',1)
    logistic=[n.value for n in ast.walk(point_fn) if isinstance(n,ast.Return) and
              ast.dump(n.value)==ast.dump(ast.parse('e/(1+e)',mode='eval').body)]
    if len(logistic)!=1:raise ValueError('Actual sigmoid logistic return source binding missing')
    odds=evaluate(odds_node,dict(y=v)); exponential=evaluate(exp_node,dict(c=c,odds=odds))
    sigma=evaluate(logistic[0],dict(e=exponential))
    if s.simplify(odds+odds.subs(v,1-v))!=0:raise ArithmeticError('Actual source sigmoid odds reflection failed')
    u=s.Symbol('actual_sigmoid_odds',real=True)
    reflected=evaluate(logistic[0],dict(e=evaluate(exp_node,dict(c=c,odds=u))))
    if s.simplify(reflected+reflected.subs(u,-u)-1)!=0:raise ArithmeticError('Actual logistic source reflection failed')
    sigma_callable=lambda ctx,x:sigma
    Jprime=evaluate(forward_step,dict(J=s.Integer(0),ds=s.Integer(1),c=c,v=v,sigma_enclosure=sigma_callable))
    fprime=-evaluate(remaining_step,dict(f=s.Integer(0),ds=s.Integer(1),c=c,v=v,sigma_enclosure=sigma_callable))
    if s.simplify(fprime-Jprime+1)!=0 or s.simplify(f1-J1+1-s.Rational(1,2))!=0:
        raise ArithmeticError('Actual primitive FTC/endpoint proof of f-J+v-1/2 failed')
    Ientry,Ifull=s.symbols('actual_entry_integral actual_exit_forward_integral',real=True)
    original=SimpleNamespace(rate=1-mu,k=k,Ts=Ts,infull=dict(angular=Ientry),outfull=dict(angular=Ifull))
    env=dict(c=c,self=original,XR=s.Function('actual_original_XR')(z))
    XS=evaluate(expression('steep_waiting_C4','data','XS'),env); env['XS']=XS
    XQ=evaluate(expression('steep_waiting_C4','data','XQ'),env); env['XQ']=XQ
    XT=evaluate(expression('steep_waiting_C4','data','XT'),env)
    native_tail=evaluate(expression('steep_waiting_C4','waiting','X'),dict(env,data=dict(XT=XT),t=wait))
    native_exit=evaluate(expression('steep_waiting_C4','steep_out','X'),dict(env,data=dict(XQ=XQ),kernels=dict(angular=s.Integer(0)),J=s.Integer(0)))
    Ktail=1-eps; KQ=Ktail*s.exp(k/2)
    Km=evaluate(expression('steep_exit_stress_C3','exit_signed_kernels','Km'),
                dict(c=c,eps=eps,K0=Ktail,k=k,fc=J-v+s.Rational(1,2)))
    signed=evaluate(expression('steep_exit_stress_C3','exit_signed_kernels','IA',True,'ds*c.exp(k*distance)*Km'),
                    dict(c=c,ds=s.Integer(1),k=k,distance=v,Km=Km))
    forward=evaluate(expression('steep_waiting_C4','transition_kernels','angular',True,'ds*c.exp(rate*jc)'),
                     dict(c=c,ds=s.Integer(1),rate=k,jc=J))
    if s.simplify(s.expand_power_exp(signed-KQ*forward+s.exp(k*v)))!=0:
        raise ArithmeticError('Actual signed/forward exit angular integral density split failed')
    baseline=(s.exp(k)-1)/k
    if s.simplify(s.integrate(s.exp(k*v),(v,0,1))-baseline)!=0:
        raise ArithmeticError('Exact unit angular density integral failed')
    # Common zero-history closure fixes A0=Ktail*X_native_at_collar.
    # It is an already source-proven functional equality, not a numerical
    # assignment to a fresh terminal moment or interval overlap assertion.
    base={label+'_defect_rows':[s.Function('same_'+label+'_collar_defect')(z)] for label in ('angular','energy','pressure')}
    base['angular_defect_rows']=[Ktail*native_tail-1/k]
    stub=SimpleNamespace(constant=lambda ctx,value,order:value)
    heat=SimpleNamespace(ctx=c,k=k,a=1-k,eps=eps,delta=2-2*k,prate=3-2*k,
                         steep=SimpleNamespace(wait=wait))
    replay_env=dict(math=math,IntervalTaylor=stub,sigma_jets=lambda ctx,t:[s.Integer(0)]*5)
    def replay(stem,name):
        path=HERE/(PREFIX+stem+'.py'); hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        fn=next(n for n in ast.walk(ast.parse(path.read_text(encoding='utf8'))) if isinstance(n,ast.FunctionDef) and n.name==name)
        fn.decorator_list=[]
        exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual endpoint quotient '+name+'>','exec'),replay_env)
        return replay_env[name]
    replay('collar_Gamma_C4','product_rows'); replay('waiting_stress_C3','waiting_defect_rows')
    shape_fn=replay('steep_exit_stress_C3','exit_shape')
    defects_fn=replay('steep_exit_stress_C3','exit_defect_rows')
    kernels=dict(f=s.Rational(1,2),angular=KQ*Ifull-baseline,energy=s.Symbol('actual_signed_exit_energy'),
                 pressure=s.Symbol('actual_signed_exit_pressure'))
    shape=shape_fn(heat,s.Integer(0),kernels); defects=defects_fn(heat,base,s.Integer(0),kernels,shape)
    companion=(1/k+defects['angular_defect_rows'][0])/shape['K_rows'][0]
    if s.simplify(s.expand_power_exp(companion-native_exit))!=0:
        raise ArithmeticError('Actual companion full AQ/KQ differs from original steep_out XQ')
    if s.simplify(native_exit-(XS+Ts))!=0:raise ArithmeticError('Actual XQ/XS source relation lost')
    return dict(identities={flag:True for flag in flags},actual_source_expression_bindings=bindings,
                actual_signed_angular_integrand_equals_KQ_forward_integrand_minus_unit=True,
                same_original_sigma_callable_and_reflection_half_integral_consumed=True,
                actual_J_forward_primitive_and_half_endpoint_AST_bound=True,
                actual_logistic_return_and_reflection_AST_bound=True,
                actual_remaining_f_J_identity_derived_by_FTC_and_endpoint=True,
                actual_remaining_f_equals_original_J_minus_v_plus_half_verified=True,
                actual_signed_exit_integral_equals_KQ_forward_integral_minus_exponential_unit_integral=True,
                common_zero_history_identifies_full_A0_with_actual_native_collar_X=True,
                actual_child_AQ_over_KQ_equals_original_steep_out_XQ_verified=True,
                actual_original_XQ_equals_source_XS_plus_Ts_verified=True,
                source_function_equality_not_interval_overlap=True,input_hashes=hashes)


def power_source_bridge(exit_source):
    """Bind original chart, source ODEs and admitted complete exit moments."""
    hashes=dict(exit_source.hashes); checks={}; bindings={}
    def syntax(stem,method,target,expression):
        path=HERE/(PREFIX+stem+'.py'); hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        fn=next(n for n in ast.walk(ast.parse(path.read_text(encoding='utf8'))) if isinstance(n,ast.FunctionDef) and n.name==method)
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(v)==target for v in n.targets)]
        if sum(ast.dump(v)==ast.dump(ast.parse(expression,mode='eval').body) for v in values)!=1:
            raise ValueError('Power source changed: '+stem+'.'+target)
        bindings[stem+'.'+method+'.'+target]=True
    for target,expression in {'t':'self.Ts*phase','left':'self.Ts*(1-phase)',
        'theta':"one*self.thetaS*c.exp(-c.mpf('1.5')*t)",
        'X':"data['XS']+t",'energy':"(data['after_power']-c.mpf('.5'))*c.exp(-2*left)/2+c.mpf('.25')",
        'pressure':"data['PS']+decay_integral(c,3,t)*(self.outer.flatten.Ev2*self.thetaS**2/2)"}.items():
        syntax('steep_waiting_C4','steep_power',target,expression)
    syntax('steep_waiting_C4','__init__','self.Ts','4*(c.ln(2)-c.ln(self.delta))')
    syntax('steep_waiting_C4','__init__','self.thetaQ',"self.thetaS*c.exp(-c.mpf('1.5')*self.Ts)")
    syntax('steep_waiting_C4','data','XQ','XS+self.Ts')
    syntax('steep_waiting_C4','data','XS',"(XR+self.infull['angular'])*c.exp(-self.rate/2)")
    syntax('power_angular_C4','_packet','mixed','flatten_mixed(c,self.mu,one,[c.mpf(0)]+rates,theta,X,energy,pressure,self.flatten.Ev2)')
    a,t,T,wait,eps,thetaS,Ev2=s.symbols('a t Ts wait eps thetaS Ev0_squared',real=True)
    k=1-a; delta=2*a; p=1+delta; bh=s.Rational(1,2)+a
    thetaQ=thetaS*s.exp(-s.Rational(3,2)*T); thetaT=thetaQ*s.exp(-s.Rational(3,2)+k/2)
    theta_base=thetaT*s.exp(-bh*wait)/(1-eps); ell=T-t; q=-wait-1-ell
    KQ=(1-eps)*s.exp(k/2); K=KQ*s.exp(k*ell)
    theta=thetaS*s.exp(-s.Rational(3,2)*t)
    if s.simplify(s.expand_power_exp(theta-theta_base*s.exp(-bh*q)*K))!=0:
        raise ArithmeticError('Actual power original velocity normalization failed')
    checks['actual_original_power_velocity_normalization']=True
    # Replay guarded production assignments, with exact defining decay
    # integral; the actual helper returns a directed bound on that integral.
    tree=ast.parse((HERE/(PREFIX+'steep_waiting_C4.py')).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='steep_power')
    expressions={}
    for target in ('theta','X','energy','pressure'):
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(v)==target for v in n.targets)]
        if len(values)!=1:raise ValueError('Original power source assignment not unique: '+target)
        expressions[target]=values[0]
    XS,after,PS=s.symbols('source_XS source_after_power source_PS',real=True)
    c=SimpleNamespace(exp=s.exp,mpf=lambda v:s.Rational(str(v)))
    obj=SimpleNamespace(thetaS=thetaS,outer=SimpleNamespace(flatten=SimpleNamespace(Ev2=Ev2)))
    env=dict(c=c,self=obj,one=s.Integer(1),t=t,left=T-t,data=dict(XS=XS,after_power=after,PS=PS),
             decay_integral=lambda ctx,rate,length:(1-s.exp(-rate*length))/rate)
    source={target:eval(compile(ast.Expression(expr),'<actual power '+target+'>','eval'),{},env) for target,expr in expressions.items()}
    for label,value in {'theta':s.diff(source['theta'],t)+s.Rational(3,2)*source['theta'],
        'angular':s.diff(source['X'],t)-1,'energy':s.diff(source['energy'],t)-2*source['energy']+s.Rational(1,2),
        'pressure':s.diff(source['pressure'],t)-Ev2*source['theta']**2/2}.items():
        if s.simplify(value)!=0:raise ArithmeticError('Actual power source ODE changed: '+label)
        checks['actual_original_'+label+'_source_ODE']=True
    receipt=json.loads((HERE/(PREFIX+'steep_waiting_C4_check.json')).read_bytes())
    for flag in ('power_exit_exact_theta','power_exit_exact_X','power_exit_exact_energy','power_exit_exact_pressure'):
        if not receipt['all_passed'] or not receipt['functional_production_source_identities'][flag]:raise ValueError('Original power/exit source join missing')
        checks['consumed_'+flag]=True
    for flag in ('actual_same_complete_waiting_collar_Gamma_history_consumed','actual_original_absolute_pressure_datum_consumed',
        'actual_terminal_meridional_zeros_and_source_zero_density_FTC_consumed','exact_S_is_positive_source_not_enclosure_endpoint',
        'actual_production_steep_exit_radius_reference_verified','actual_original_velocity_reference_normalization_verified'):
        if not exit_source.bridge[flag]:raise ValueError('Actual full exit input missing: '+flag)
        checks['consumed_'+flag]=True
    for flag,value in exit_source.bridge['exact_Ev0_squared_pressure_amplitude_binding_consumed'].items():
        if not value:raise ValueError('Exact power pressure amplitude input missing: '+flag)
        checks['consumed_exact_amplitude_'+flag]=True
    for flag,value in exit_source.bridge['exact_source_factor_identities_consumed'].items():
        if not value:raise ValueError('Exact power stress source factor missing: '+flag)
        checks['consumed_exact_factor_'+flag]=True
    waiting=exit_source.waiting_source; collar=waiting.collar_source
    if not collar.bridge['identities']['actual_zero_angular_history_defect_transfers_to_all_collar_offsets']:
        raise ValueError('Actual original angular zero-history closure required for power inlet correlation')
    if not waiting.bridge['actual_waiting_endpoint_full_moment_history_consumed']:
        raise ValueError('Same original waiting angular endpoint required')
    if not waiting.proof['identities']['original_waiting_angular_normalization']:
        raise ValueError('Original waiting angular normalization required')
    if not exit_source.proof['identities']['original_angular_KX_normalization_ODE']:
        raise ValueError('Original exit angular normalization ODE required')
    if not exit_source.bridge['original_source_ODEs_replayed']['angular'] or not exit_source.bridge['full_future_normalized_moment_ODE_uniqueness_used']:
        raise ValueError('Original full exit angular-history transfer required')
    checks.update(consumed_actual_zero_angular_history_through_collar_waiting_exit=True,
                  source_XQ_equals_full_AQ_over_KQ_by_common_history_and_ODE_uniqueness=True)
    syntax('steep_power_stress_C3','terminal',"self.cache[key]['original_XS']","self.steep.data(Z)['XS']")
    syntax('steep_power_stress_C3','power_angular_rows','X',"terminal['original_XS']+one*(heat.steep.Ts*phase)")
    syntax('steep_power_stress_C3','power_defect_rows','X','power_angular_rows(heat,terminal,phase)[0]')
    syntax('steep_power_stress_C3','power_defect_rows','Ad',"shape['K_rows'][0]*X-one/heat.k")
    syntax('steep_power_stress_C3','power_stress_rows','X','power_angular_rows(heat,terminal,phase)[0]')
    phase=s.Symbol('source_phase',real=True); XS=s.Function('actual_source_XS')(s.Symbol('Z'))
    full_XQ=s.Symbol('same_full_AQ_over_KQ')
    difference=(full_XQ-T*(1-phase))-(XS+T*phase)
    if s.simplify(difference.subs(full_XQ,XS+T))!=0:
        raise ArithmeticError('Same full-future power angular inlet correlation failed')
    checks['source_XS_plus_Ts_phase_is_same_full_future_A_over_K_function']=True
    endpoint=power_angular_endpoint_binding(exit_source); hashes.update(endpoint['input_hashes'])
    path=HERE/(PREFIX+'global_physical_assembly.py'); tree=ast.parse(path.read_text(encoding='utf8')); radius_env={}
    for node in tree.body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ('MICRO','POST','PULSE'):radius_env[target.id]=ast.literal_eval(node.value)
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='radius')
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual power phase radius>','exec'),radius_env)
    mu,origin,length,phase=s.symbols('mu logRp Lrel phase',real=True)
    native=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify),logRp=origin,params=SimpleNamespace(mu=mu))
    provider=SimpleNamespace(outer=SimpleNamespace(Lrel=length),Ts=T,wait=wait)
    actualq=radius_env['radius'](native,'steep_power',phase,{},provider)[0]-radius_env['radius'](native,'waiting',s.Integer(1),{},provider)[0]
    if s.simplify(actualq-(-wait-1-T*(1-phase)))!=0:raise ArithmeticError('Original power radius source changed')
    checks['actual_original_power_phase_radius_reference']=True
    hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    join=power_exit_join_binding(); hashes.update(join['input_hashes'])
    return dict(identities=checks,actual_power_source_bindings=bindings,actual_power_exit_formula_join=join,
                same_full_exit_endpoint_and_normalized_moment_ODE_uniqueness_used=True,
                actual_meridional_source_and_terminal_zero_primitive_route_consumed=True,
                actual_absolute_pressure_datum_and_exact_Ev0_units_consumed=True,
                actual_source_XS_correlation_is_same_full_angular_future=True,
                actual_full_future_AQ_over_KQ_original_source_endpoint_binding=endpoint,
                huge_independent_AQ_over_KQ_minus_Ts_not_used_to_define_inlet=True,
                original_source_coefficients_or_velocities_changed=False,input_hashes=hashes)


class CompliantSteepPowerStressC3:
    @source_precision
    def __init__(self):
        self.exit_source=CompliantSteepExitStressC3(); self.heat=self.exit_source.heat; self.steep=self.heat.steep
        self.ctx=self.heat.ctx; self.family,self.source=self.exit_source.family,self.exit_source.source
        self.cache={}; self.hashes=dict(self.exit_source.hashes)
        name=PREFIX+'steep_exit_stress_C3_check.json'; receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or not receipt['actual_original_steep_exit_similarity_stress_recovered']:
            raise ValueError('Accepted actual complete exit moments required')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Power source input changed: '+source)
        self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.proof=power_transport_identities(); self.bridge=power_source_bridge(self.exit_source)
        self.hashes.update(self.bridge['input_hashes']); self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def terminal(self,Z):
        key=tuple(endpoints(self.ctx.mpf(Z)))
        if key not in self.cache:
            point=self.exit_source.steep_out(Z,0)
            self.cache[key]={label:[jet] for label,jet in point['steep_exit_future_defect_zeroth_Taylor'].items()}
            self.cache[key]['KQ']=point['steep_exit_shape']['K_rows'][0][0]
            self.cache[key]['original_XS']=self.steep.data(Z)['XS']
        return self.cache[key]

    @source_precision
    def steep_power(self,Z,phase):
        c=self.ctx; Z=c.mpf(Z); phase=c.mpf(phase); self.steep._phase(phase)
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1:raise ValueError('Source Z in[-1,1] required')
        terminal=self.terminal(Z); shape=power_shape(self.heat,terminal,phase)
        defects=power_defect_rows(self.heat,terminal,phase,shape)
        rows=power_stress_rows(self.heat,terminal,phase,shape,Z); pressure=power_pressure_rows(self.heat,terminal,phase,shape)
        point=dict(self.steep.steep_power(Z,phase)); one=IntervalTaylor.constant(c,1,5)
        for label in ('angular','energy'):
            point['original_forward_'+label+'_Taylor']=point[label+'_Taylor']
            point['original_forward_'+label+'_y_derivative_Taylor']=point[label+'_y_derivative_Taylor']
        point['original_forward_pressure_over_Pstar_squared_Taylor']=point['pressure_over_Pstar_squared_Taylor']
        point['original_forward_pressure_mixed_bounds']=point['physical_mixed_derivatives_total_order_le4'][P]
        A=[one/self.heat.k+defects['angular_defect_rows'][0]]+defects['angular_defect_rows'][1:]
        E=[one/self.heat.delta+defects['energy_defect_rows'][0]]+defects['energy_defect_rows'][1:]
        angular=power_angular_rows(self.heat,terminal,phase)
        energy=quotient_rows(E,[jet*2 for jet in product_rows(shape['K_rows'],shape['K_rows'])])
        mixed=dict(point['physical_mixed_derivatives_total_order_le4']); fields=dict(point['physical_velocity_and_pressure_y_derivative_Taylor'])
        mixed[P]={'y'+str(j)+'_Z'+str(n):jet[n]*math.factorial(n) for j,jet in enumerate(pressure) for n in range(5-j)}
        fields[P]=[jet.truncate(4-j) for j,jet in enumerate(pressure)]
        point.update(angular_Taylor=angular[0],angular_y_derivative_Taylor=[jet.truncate(5-j) for j,jet in enumerate(angular)],
                     energy_Taylor=energy[0],energy_y_derivative_Taylor=[jet.truncate(5-j) for j,jet in enumerate(energy)],
                     pressure_over_Pstar_squared_Taylor=pressure[0],pressure_y_derivative_axial5_Taylor=pressure,
                     physical_mixed_derivatives_total_order_le4=mixed,physical_velocity_and_pressure_y_derivative_Taylor=fields,
                     steep_power_phase=phase,steep_power_remaining_to_exit=shape['remaining'],
                     steep_power_offset_from_Rtail=-self.steep.wait-1-shape['remaining'],steep_power_shape=shape,
                     steep_power_similarity_stress_mixed3_factored={label:{'y'+str(j)+'_Z'+str(n):rows[label][j][n]*math.factorial(n)
                         for j in range(4) for n in range(4-j)} for label in ('theta','axial')},
                     steep_power_similarity_stress_y_derivative_Taylor={label:[jet.truncate(3-j) for j,jet in enumerate(rows[label])]
                         for label in ('theta','axial','theta_inertial','theta_shear')},
                     steep_power_future_defect_zeroth_Taylor={label:values[0] for label,values in defects.items()},
                     steep_power_correlated_source_inlet_X=terminal['original_XS'],
                     actual_same_full_angular_future_recovered_with_source_XS_identity=True,
                     steep_power_meridional_moments_and_velocities={label:one*0 for label in ('Mz','Mtheta_z','Ur','Uz')},
                     steep_power_shear_strength_kappa_minus2=c.mpf(2),
                     actual_original_steep_power_similarity_stress_recovered=True,
                     actual_steep_power_absolute_pressure_same_source_mixed4_available=True,
                     steep_power_exit_stress_mixed3_join_verified=True,steep_power_exit_pressure_mixed4_join_verified=True,
                     same_complete_waiting_collar_Gamma_future_used=True,original_source_histories_retained=True,
                     source_defined_positive_stress_factors=dict(theta='sqrt(R/2)*B',axial='sqrt(R/2)*B^2',
                         B='Ev0*theta_base*exp(-bh*q)',R='Rtail*exp(q)',q='-wait-1-Ts*(1-phase)',
                         exact_Ev0='Pstar*U*exp(-13/(2*mu)-13); Ev2 is an enclosure only'),
                     steep_power_cone_certified=False,steep_power_physical_stress_remainder_identity_verified=False,
                     global_admissible_stress_lift_constructed=False,physical_energy_integral_certified=False,temporal_recursion=False)
        return point

    def report(self):
        keys=('Z','steep_power_phase','steep_power_remaining_to_exit','steep_power_offset_from_Rtail','steep_power_shape',
              'steep_power_similarity_stress_mixed3_factored','steep_power_similarity_stress_y_derivative_Taylor',
              'steep_power_future_defect_zeroth_Taylor','pressure_y_derivative_axial5_Taylor',
              'angular_y_derivative_Taylor','steep_power_correlated_source_inlet_X',
              'steep_power_meridional_moments_and_velocities','steep_power_shear_strength_kappa_minus2','source_defined_positive_stress_factors')
        def summary(z,phase):
            point=self.steep_power(z,phase); return {key:point[key] for key in keys}
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,domain=DOMAIN,
                    scope='Original steep_power phase[0,1], y=Ts*phase, source Z[-1,1]; complete-future stress mixed3 and pressure mixed4',
                    samples=[summary(z,phase) for z,phase in (('0','0'),('.5','.5'),('-.5','1'))],whole_steep_power=summary([-1,1],[0,1]),
                    transport_identities=self.proof,steep_power_source_bridge=self.bridge,
                    actual_original_steep_power_similarity_stress_recovered=True,
                    actual_steep_power_absolute_pressure_same_source_mixed4_available=True,
                    actual_same_full_angular_future_recovered_with_source_XS_identity=True,
                    steep_power_exit_stress_mixed3_join_verified=True,steep_power_exit_pressure_mixed4_join_verified=True,
                    steep_power_cone_certified=False,steep_power_physical_stress_remainder_identity_verified=False,
                    global_admissible_stress_lift_constructed=False,physical_energy_integral_certified=False,
                    temporal_recursion=False,input_hashes=self.hashes)


def run():
    result=CompliantSteepPowerStressC3().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original full Ts steep power common-future stress mixed3 and absolute pressure mixed4 generated',flush=True)
    return result


if __name__=='__main__':run()
