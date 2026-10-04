"""Original steep exit: common full-future moments, stress and pressure.

The actual waiting/collar future is transported backward across the original
sigma transition. Unit baselines cancel before interval evaluation. Signed
finite integrals below extend those same moments; they are not replacement
terminal data or a truncation of the Gamma future.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_waiting_stress_C3 import (
    CompliantWaitingStressC3, waiting_defect_rows, source_precision)
from lei_ren_part1_paper_compliant_collar_stress_C3 import (
    collar_stress_rows, collar_moment_stress_identities)
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


def exit_signed_kernels(c,t,k,delta,eps,cells=128):
    """Directed remaining sigma primitive and three signed weighted integrals.

    f(v)=int_v^1(1-sigma), hence f=J-v+1/2 by the original symmetry.
    Traverse from the flat right endpoint. Every cell has the correlated
    positive length (1-t)/cells. v-t is (1-t)*u, never an independent
    subtraction of endpoint boxes. Signed K-1 and K^2-1 use expm1.
    """
    t=c.mpf(t); lo,hi=endpoints(t)
    if lo<0 or hi>1:raise ValueError('Original steep exit t in[0,1] required')
    if not isinstance(cells,int) or cells<1:raise ValueError('Positive cell count required')
    length=1-t; ds=length/cells; f=c.mpf(0)
    IA=c.mpf(0); IE=c.mpf(0); IP=c.mpf(0); p=1+delta; K0=1-eps
    if endpoints(length)[1]>0:
        for i in reversed(range(cells)):
            a=t+length*i/cells; b=t+length*(i+1)/cells
            v=c.mpf([max(mp.mpf(0),endpoints(a)[0]),min(mp.mpf(1),endpoints(b)[1])])
            nextf=f+ds*(1-sigma_enclosure(c,v))
            fc=c.mpf([endpoints(f)[0],min(mp.mpf('.5'),endpoints(nextf)[1])])
            Km=-eps+K0*c.expm1(k*fc)
            Qm=-2*eps+eps**2+K0**2*c.expm1(2*k*fc)
            distance=length*c.mpf([i,i+1])/cells
            IA+=ds*c.exp(k*distance)*Km
            IE+=ds*c.exp(-delta*distance)*Qm
            IP+=ds*c.exp(-p*distance)*Qm
            f=nextf
    if (lo,hi)==(mp.mpf(0),mp.mpf(0)):f=c.mpf('.5')
    else:f=c.mpf([endpoints(f)[0],min(mp.mpf('.5'),endpoints(f)[1])])
    return dict(f=f,angular=IA,energy=IE,pressure=IP,
                exact_f_definition='integral_t^1 (1-original_sigma(v))dv = J(t)-t+1/2',
                exact_angular_definition='integral_t^1 exp(k*(v-t))*(K(v)-1)dv',
                exact_energy_definition='integral_t^1 exp(-delta*(v-t))*(K(v)^2-1)dv',
                exact_pressure_definition='integral_t^1 exp(-p*(v-t))*(K(v)^2-1)dv',
                positive_correlated_cell_lengths=True,signed_integrands=True,cells=cells)


def exit_shape(heat,t,kernels):
    """Actual K=K0 exp(k*(J-t+1/2)), ordinary y rows through four."""
    c=heat.ctx; one=IntervalTaylor.constant(c,1,5); eps=heat.eps; K0=1-eps
    sig=sigma_jets(c,t)
    g=[heat.k*(sig[0]-1)]+[heat.k*sig[j]*math.factorial(j) for j in range(1,4)]
    Km0=-eps+K0*c.expm1(heat.k*kernels['f'])
    K=[one*(1+Km0)]
    for n in range(1,5):
        K.append(sum((K[n-1-j]*g[j]*math.comb(n-1,j) for j in range(n)),one*0))
    Km=[one*Km0]+K[1:]
    square=product_rows(K,K)
    Qm=[one*(-2*eps+eps**2+K0**2*c.expm1(2*heat.k*kernels['f']))]+square[1:]
    return dict(K_rows=K,K_defect_rows=Km,K_squared_defect_rows=Qm,
                log_K_rate_rows=g,original_sigma_jets=sig)


def exit_defect_rows(heat,base,t,kernels,shape):
    """Backward full-future FTC from the SAME waiting left endpoint."""
    c=heat.ctx; one=IntervalTaylor.constant(c,1,5); length=1-t
    left=waiting_defect_rows(heat,base,-heat.steep.wait)
    Ad=left['angular_defect_rows'][0]*c.exp(heat.k*length)-one*kernels['angular']
    Ed=left['energy_defect_rows'][0]*c.exp(-heat.delta*length)+one*kernels['energy']
    Pd=left['pressure_defect_rows'][0]*c.exp(-heat.prate*length)+one*kernels['pressure']/2
    A=[Ad]; E=[Ed]; Pr=[Pd]; Km=shape['K_defect_rows']; Qm=shape['K_squared_defect_rows']
    for j in range(4):
        A.append(Km[j]-A[j]*heat.k)
        E.append(E[j]*heat.delta-Qm[j])
        Pr.append(Pr[j]*heat.prate-Qm[j]/2)
    return dict(angular_defect_rows=A,energy_defect_rows=E,pressure_defect_rows=Pr,
                K_defect_rows=Km,K_squared_defect_rows=Qm)


def quotient_rows(numerator,denominator):
    """Ordinary y derivatives of axial Taylor-valued N/D."""
    out=[numerator[0]/denominator[0]]
    for n in range(1,len(numerator)):
        out.append((numerator[n]-sum((denominator[j]*out[n-j]*math.comb(n,j)
                    for j in range(1,n+1)),numerator[0]*0))/denominator[0])
    return out


def exit_pressure_rows(heat,defects,shape,q):
    one=IntervalTaylor.constant(heat.ctx,1,5)
    current=one/(2*heat.prate)+defects['pressure_defect_rows'][0]
    return pressure_y_rows(shape['K_rows'],current,heat.prate,heat.pressure_scale,q)


def exit_transport_identities():
    """Full future split, original normalized ODEs and source shear sign."""
    a,t,z=s.symbols('a t Z',real=True); delta=2*a; k=1-a; p=1+delta
    eps=s.symbols('eps',real=True); K=s.Function('K')(t); sig=s.Function('sigma')(t)
    A1=s.Function('Ad_waiting_left')(z); E1=s.Function('Ed_waiting_left')(z); P1=s.Function('Pd_waiting_left')(z)
    IA=s.Function('signed_angular_integral')(t)
    IE=s.Function('signed_energy_integral')(t); IP=s.Function('signed_pressure_integral')(t)
    Ad=s.exp(k*(1-t))*A1-IA
    Ed=s.exp(-delta*(1-t))*E1+IE; Pd=s.exp(-p*(1-t))*P1+IP/2
    rules={s.diff(IA,t):-(K-1)-k*IA,
           s.diff(IE,t):-(K*K-1)+delta*IE,s.diff(IP,t):-(K*K-1)+p*IP}
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Steep exit identity failed: '+name)
        proofs[name]=True
    for label,row,rhs in (('angular',Ad,K-1-k*Ad),('energy',Ed,delta*Ed-(K*K-1)),
                          ('pressure',Pd,p*Pd-(K*K-1)/2)):
        zero('actual_backward_'+label+'_FTC',s.diff(row,t).subs(rules)-rhs)
    for label,row,start in (('angular',Ad,A1),('energy',Ed,E1),('pressure',Pd,P1)):
        zero('same_complete_future_'+label+'_waiting_left',row.subs({t:1,IA:0,IE:0,IP:0})-start)
    X=s.Function('original_X')(t); en=s.Function('original_half_energy')(t)
    Krate=k*(sig-1)
    zero('original_angular_KX_normalization_ODE',Krate*K*X+K*(1-k*sig*X)-(K-k*K*X))
    zero('original_energy_2K2_normalization_ODE',4*K*K*Krate*en+2*K*K*(-s.Rational(1,2)+(2-2*k*sig)*en)
         -(delta*(2*K*K*en)-K*K))
    P=s.Function('absolute_pressure')(t); C,q=s.symbols('C q',real=True)
    zero('original_absolute_pressure_normalization_ODE',
         -C*K*K*s.exp(-p*q)/2*s.exp(p*q)/C+p*(-P*s.exp(p*q)/C)
         -(p*(-P*s.exp(p*q)/C)-K*K/2))
    f=s.Function('remaining_sigma_primitive')(t)
    sourceK=(1-eps)*s.exp(k*f)
    zero('actual_K_log_rate_from_remaining_sigma_primitive',
         s.diff(sourceK,t).subs(s.diff(f,t),sig-1)-Krate*sourceK)
    zero('stable_K_defect_is_actual_source_K_minus1',-eps+(1-eps)*(s.exp(k*f)-1)-(sourceK-1))
    zero('stable_quadratic_defect_is_actual_source_K_squared_minus1',
         -2*eps+eps**2+(1-eps)**2*(s.exp(2*k*f)-1)-(sourceK**2-1))
    zero('actual_shear_margin',delta-2*Krate-(delta+2*k*(1-sig)))
    odds=1/(1-t)**2-1/t**2
    zero('original_sigmoid_reflection_odds',odds.subs(t,1-t)+odds)
    positive_e=s.symbols('positive_logistic_exponential',positive=True)
    logistic=positive_e/(1+positive_e)
    zero('original_sigmoid_positive_complement',1-logistic-1/(1+positive_e))
    zero('original_sigmoid_reflection_and_half_integral',logistic+1/(1+positive_e)-1)
    return dict(identities=proofs,actual_original_normalized_moment_ODEs_verified=True,
                complete_future_endpoint_and_backward_FTC_verified=True,
                original_sigmoid_between0and1_with_flat_endpoints=True,
                original_sigmoid_reflection_gives_half_integral=True,
                source_shear_margin_at_least_delta_from_0_le_sigma_le1=True,
                interval_overlap_used_as_functional_proof=False)


def exit_join_binding():
    """Execute actual row ASTs at the flat endpoint on arbitrary axial data."""
    a,eps,z,wait,C=s.symbols('a eps Z wait C',real=True)
    c=SimpleNamespace(exp=s.exp,expm1=lambda value:s.exp(value)-1,mpf=lambda value:s.Rational(str(value)))
    heat=SimpleNamespace(ctx=c,a=a,eps=eps,S=s.Symbol('source_inverse_Rtail'),delta=2*a,
                         k=1-a,prate=1+2*a,pressure_scale=C,steep=SimpleNamespace(wait=wait))
    stub=SimpleNamespace(constant=lambda ctx,value,order:value,variable=lambda ctx,value,order:value)
    environment=dict(s=s,math=math,c=c,IntervalTaylor=stub,sigma_jets=lambda ctx,t:[s.Integer(1)]+[s.Integer(0)]*4,
                     axial_derivative=lambda value:s.diff(value,z))
    hashes={}
    def replay(stem,name,remove_context=False):
        path=HERE/(PREFIX+stem+'.py'); hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        fn=next(n for n in ast.walk(ast.parse(path.read_text(encoding='utf8'))) if isinstance(n,ast.FunctionDef) and n.name==name)
        fn.decorator_list=[]
        if remove_context:
            nodes=[n for n in fn.body if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='c'
                   and ast.unparse(n.value)=='pressure_numerator.ctx']
            if len(nodes)!=1:raise ValueError('Pressure context lookup changed')
            fn.body.remove(nodes[0])
        exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual exit join '+name+'>','exec'),environment)
        return environment[name]
    replay('collar_Gamma_C4','product_rows'); replay('collar_stress_C3','shifted_rows')
    generic=replay('collar_stress_C3','collar_stress_rows')
    replay('heat_pressure_C4','pressure_y_rows',True)
    waitingdef=replay('waiting_stress_C3','waiting_defect_rows')
    waitingstress=replay('waiting_stress_C3','waiting_stress_rows')
    waitingpressure=replay('waiting_stress_C3','waiting_pressure_rows')
    shape_fn=replay('steep_exit_stress_C3','exit_shape')
    defect_fn=replay('steep_exit_stress_C3','exit_defect_rows')
    pressure_fn=replay('steep_exit_stress_C3','exit_pressure_rows')
    base={label+'_defect_rows':[s.Function(label+'_full_terminal')(z)] for label in ('angular','energy','pressure')}
    kernels=dict(f=s.Integer(0),angular=s.Integer(0),energy=s.Integer(0),pressure=s.Integer(0))
    shape=shape_fn(heat,s.Integer(1),kernels); defects=defect_fn(heat,base,s.Integer(1),kernels,shape)
    left=generic(heat,shape,defects,z,-wait); right=waitingstress(heat,base,z,-wait)
    proofs={}
    for label in ('theta','axial'):
        for j in range(4):
            for n in range(4-j):
                if s.simplify(s.diff(left[label][j]-right[label][j],z,n))!=0:raise ArithmeticError('Actual exit stress join failed')
                proofs[label+'_y'+str(j)+'_Z'+str(n)]=True
    # The production pressure context lookup was already replayed above;
    # exit_pressure_rows uses only the remaining actual scalar expressions.
    lp=pressure_fn(heat,defects,shape,-wait); rp=waitingpressure(heat,base,-wait)
    for j in range(5):
        for n in range(5-j):
            if s.simplify(s.diff(lp[j]-rp[j],z,n))!=0:raise ArithmeticError('Actual exit pressure join failed')
            proofs['pressure_y'+str(j)+'_Z'+str(n)]=True
    return dict(identities=proofs,actual_steep_exit_waiting_stress_mixed3_AST_join_verified=True,
                actual_steep_exit_waiting_pressure_mixed4_AST_join_verified=True,
                arbitrary_full_terminal_axial_functions_not_sample_overlap=True,input_hashes=hashes)


def exit_source_bridge(waiting):
    """Original source normalization, complete history and exact joins."""
    hashes=dict(waiting.hashes); tree=ast.parse((HERE/(PREFIX+'steep_waiting_C4.py')).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='steep_out')
    expressions={}
    for target in ('theta','X','energy','pressure'):
        rows=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(v)==target for v in n.targets)]
        if len(rows)!=1:raise ValueError('Original exit expression changed: '+target)
        expressions[target]=rows[0]
    a,t,wait,J=s.symbols('a t wait J',real=True); k=1-a; delta=2*a; p=1+delta; bh=s.Rational(1,2)+a
    eps,thetaQ,Ev2=s.symbols('eps thetaQ Ev0_squared',real=True)
    XQ,WF,PQ,IA,IE,IP=s.symbols('XQ waiting_future PQ forward_angular remaining_energy forward_pressure',real=True)
    obj=SimpleNamespace(ctx=SimpleNamespace(exp=s.exp,mpf=lambda v:s.Rational(str(v))),thetaQ=thetaQ,k=k,delta=delta,
                        outer=SimpleNamespace(flatten=SimpleNamespace(Ev2=Ev2)))
    env=dict(self=obj,c=obj.ctx,one=s.Integer(1),t=t,J=J,
             data=dict(XQ=XQ,waiting_future=WF,PQ=PQ),kernels=dict(angular=IA,remaining_energy=IE,pressure=IP))
    source={name:eval(compile(ast.Expression(expr),'<original exit '+name+'>','eval'),{},env) for name,expr in expressions.items()}
    thetaT=thetaQ*s.exp(-s.Rational(3,2)+k/2)
    theta_base=thetaT*s.exp(-bh*wait)/(1-eps); q=-wait-1+t
    K=(1-eps)*s.exp(k*(J-t+s.Rational(1,2)))
    if s.simplify(s.expand_power_exp(source['theta']-theta_base*s.exp(-bh*q)*K))!=0:
        raise ArithmeticError('Actual exit velocity normalization failed')
    # Bind actual source ODEs, not just names from an old check receipt.
    tt=s.symbols('tt',real=True); JJ=s.Function('J')(tt); sig=s.Function('sigma')(tt)
    aa=s.Function('angular_kernel')(tt); ee=s.Function('remaining_energy_kernel')(tt); pp=s.Function('pressure_kernel')(tt)
    substitutions={t:tt,J:JJ,IA:aa,IE:ee,IP:pp}
    rules={s.diff(JJ,tt):sig,s.diff(aa,tt):s.exp(k*JJ),
           s.diff(ee,tt):-s.exp(-2*tt+2*k*JJ),s.diff(pp,tt):s.exp(-3*tt+2*k*JJ)/2}
    th,xx,en,pr=[source[name].subs(substitutions) for name in ('theta','X','energy','pressure')]
    comparisons=dict(theta=s.diff(th,tt).subs(rules)-(-s.Rational(3,2)+k*sig)*th,
                     angular=s.diff(xx,tt).subs(rules)-(1-k*sig*xx),
                     energy=s.diff(en,tt).subs(rules)-(-s.Rational(1,2)+(2-2*k*sig)*en),
                     pressure=s.diff(pr,tt).subs(rules)-Ev2*th*th/2)
    for name,value in comparisons.items():
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Original exit ODE changed: '+name)
    receipt=json.loads((HERE/(PREFIX+'steep_waiting_C4_check.json')).read_bytes())
    if not receipt['all_passed'] or not receipt['angular_steep_and_internal_joins_certified']:raise ValueError('Original exit/waiting field join missing')
    for name,digest in receipt['input_hashes'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Exit original source changed: '+name)
    for flag in ('exit_waiting_exact_X','exit_waiting_exact_energy','exit_waiting_exact_pressure','waiting_exit_full_infinite_Gamma_future'):
        if not receipt['functional_production_source_identities'][flag]:raise ValueError('Original full history join missing: '+flag)
    pressure_bridge=waiting.collar_source.pressure.bridge
    if not pressure_bridge['complete_defining_function_history_bridge_verified']:
        raise ValueError('Actual absolute pressure amplitude/history bridge missing')
    amplitude_flags=('production_C4_pressure_scale','same_actual_heat_velocity_amplitude',
                     'Rtail_pressure_units_equal_theta_base_squared',
                     'production_compliant_flatten_mixed_C4_self.logEv2_parts',
                     'production_compliant_corrected_outer_field_self.logEv0Pstar2_parts')
    for flag in amplitude_flags:
        if not pressure_bridge['identities'][flag]:raise ValueError('Exact pressure amplitude link missing: '+flag)
    bindings={}
    def syntax(stem,method,target,expression,augmented=False):
        path=HERE/(PREFIX+stem+'.py'); hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        method_node=next(n for n in ast.walk(ast.parse(path.read_text(encoding='utf8')))
                         if isinstance(n,ast.FunctionDef) and n.name==method)
        values=[n.value for n in ast.walk(method_node)
                if (not augmented and isinstance(n,ast.Assign) and any(ast.unparse(v)==target for v in n.targets))
                or (augmented and isinstance(n,ast.AugAssign) and isinstance(n.op,ast.Add) and ast.unparse(n.target)==target)]
        expected=ast.dump(ast.parse(expression,mode='eval').body)
        if sum(ast.dump(v)==expected for v in values)!=1:raise ValueError('Actual exit formula changed: '+stem+':'+target)
        bindings[stem+'.'+method+'.'+target]=True
    syntax('collar_Gamma_C4','__init__','self.Ev2','self.outer.flatten.Ev2')
    syntax('collar_Gamma_C4','__init__','self.outer','self.steep.outer')
    for target,expr in {
        'Ad':"left['angular_defect_rows'][0]*c.exp(heat.k*length)-one*kernels['angular']",
        'Ed':"left['energy_defect_rows'][0]*c.exp(-heat.delta*length)+one*kernels['energy']",
        'Pd':"left['pressure_defect_rows'][0]*c.exp(-heat.prate*length)+one*kernels['pressure']/2"}.items():
        syntax('steep_exit_stress_C3','exit_defect_rows',target,expr)
    syntax('steep_exit_stress_C3','exit_shape','Km0',"-eps+K0*c.expm1(heat.k*kernels['f'])")
    syntax('steep_exit_stress_C3','exit_shape','Qm',"[one*(-2*eps+eps**2+K0**2*c.expm1(2*heat.k*kernels['f']))]+square[1:]")
    syntax('steep_exit_stress_C3','exit_shape','sig','sigma_jets(c,t)')
    syntax('steep_exit_stress_C3','exit_signed_kernels','nextf','f+ds*(1-sigma_enclosure(c,v))')
    for target,expr in {
        'fc':"c.mpf([endpoints(f)[0],min(mp.mpf('.5'),endpoints(nextf)[1])])",
        'Km':'-eps+K0*c.expm1(k*fc)',
        'Qm':'-2*eps+eps**2+K0**2*c.expm1(2*k*fc)',
        'distance':'length*c.mpf([i,i+1])/cells'}.items():
        syntax('steep_exit_stress_C3','exit_signed_kernels',target,expr)
    syntax('axial_pulse_field','point','odds','-1/y**2+1/(1-y)**2')
    syntax('axial_pulse_field','point','e','c.exp(odds)')
    syntax('flat_pulse_derivatives','_sigma_left','odds','1/(1-x)**2-1/x**2')
    syntax('flat_pulse_derivatives','_sigma_left','value','e/(1+e)')
    syntax('flat_pulse_derivatives','sigma_jets','left','_sigma_left(c,1-c.mpf([a,b]))')
    collar=waiting.collar_source
    if not collar.bridge['exact_S_not_cap_endpoint']:raise ValueError('Actual exact inverse radius missing')
    source_factors=('consumed_exact_heat_radius_and_xi_binding_verified',
                    'actual_pressure_and_energy_have_same_Ev0_squared_not_runtime_cap_value')
    for flag in source_factors:
        if not collar.bridge['identities'][flag]:raise ValueError('Actual steep exit source factor missing: '+flag)
    for target,expr in {'IA':'ds*c.exp(k*distance)*Km','IE':'ds*c.exp(-delta*distance)*Qm',
                        'IP':'ds*c.exp(-p*distance)*Qm'}.items():
        syntax('steep_exit_stress_C3','exit_signed_kernels',target,expr,True)
    # Actual production chart, relative to the SAME waiting tail radius.
    path=HERE/(PREFIX+'global_physical_assembly.py'); hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    radius_tree=ast.parse(path.read_text(encoding='utf8')); radius_env={}
    for node in radius_tree.body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ('MICRO','POST','PULSE'):
                    radius_env[target.id]=ast.literal_eval(node.value)
    radius=next(n for n in ast.walk(radius_tree) if isinstance(n,ast.FunctionDef) and n.name=='radius')
    exec(compile(ast.Module(body=[radius],type_ignores=[]),'<actual steep exit radius>','exec'),radius_env)
    Ts,Lrel,mu,origin=s.symbols('Ts Lrel mu logRp',real=True)
    native=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify),logRp=origin,params=SimpleNamespace(mu=mu))
    provider=SimpleNamespace(outer=SimpleNamespace(Lrel=Lrel),Ts=Ts,wait=wait)
    actualq=radius_env['radius'](native,'steep_exit',t,{},provider)[0]-radius_env['radius'](native,'waiting',s.Integer(1),{},provider)[0]
    if s.simplify(actualq-q)!=0:raise ArithmeticError('Actual steep exit radius reference changed')
    join=exit_join_binding(); hashes.update(join['input_hashes'])
    return dict(original_source_ODEs_replayed={name:True for name in comparisons},
                source_bindings=bindings,actual_production_steep_exit_radius_reference_verified=True,
                exact_Ev0_squared_pressure_amplitude_binding_consumed={flag:True for flag in amplitude_flags},
                exact_source_factor_identities_consumed={flag:True for flag in source_factors},
                exact_S_definition='exp(-(actual logRref+13/mu+tail_finite)) = 1/Rtail > 0',
                exact_S_is_positive_source_not_enclosure_endpoint=True,
                original_sigma_between0and1_reflection_and_flat_endpoints_bound=True,
                actual_original_velocity_reference_normalization_verified=True,
                actual_same_complete_waiting_collar_Gamma_history_consumed=waiting.bridge['actual_waiting_endpoint_full_moment_history_consumed'],
                actual_original_absolute_pressure_datum_consumed=waiting.bridge['actual_waiting_absolute_pressure_endpoint_consumed'],
                actual_terminal_meridional_zeros_and_source_zero_density_FTC_consumed=waiting.bridge['actual_meridional_terminal_zeros_and_zero_density_FTC_consumed'],
                actual_original_exit_waiting_field_and_full_moment_endpoint_join_consumed=True,
                full_future_normalized_moment_ODE_uniqueness_used=True,
                actual_steep_exit_waiting_formula_join=join,
                original_datum_coefficients_or_velocity_changed=False,input_hashes=hashes)


class CompliantSteepExitStressC3:
    @source_precision
    def __init__(self,cells=128):
        self.waiting_source=CompliantWaitingStressC3(); self.heat=self.waiting_source.heat
        self.steep=self.heat.steep; self.ctx=self.heat.ctx; self.family=self.heat.family; self.source=self.heat.source
        self.cells=cells; self.cache={}; self.hashes=dict(self.waiting_source.hashes)
        name=PREFIX+'waiting_stress_C3_check.json'; receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or not receipt['actual_original_waiting_similarity_stress_recovered']:raise ValueError('Accepted waiting complete moments required')
        for path,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Exit prerequisite changed: '+path)
        self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.proof=exit_transport_identities(); self.general_proof=collar_moment_stress_identities()
        self.bridge=exit_source_bridge(self.waiting_source); self.hashes.update(self.bridge['input_hashes'])
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def kernels(self,t):
        key=tuple(endpoints(self.ctx.mpf(t)))
        if key not in self.cache:
            self.cache[key]=exit_signed_kernels(self.ctx,t,self.heat.k,self.heat.delta,self.heat.eps,self.cells)
        return self.cache[key]

    @source_precision
    def steep_out(self,Z,t):
        c=self.ctx; Z=c.mpf(Z); t=c.mpf(t); self.steep._phase(t)
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1:raise ValueError('Source Z in[-1,1] required')
        one=IntervalTaylor.constant(c,1,5); kernels=self.kernels(t); shape=exit_shape(self.heat,t,kernels)
        base=self.waiting_source.endpoint(Z); defects=exit_defect_rows(self.heat,base,t,kernels,shape)
        q=-self.steep.wait-1+t; rows=collar_stress_rows(self.heat,shape,defects,Z,q)
        pressure=exit_pressure_rows(self.heat,defects,shape,q); point=dict(self.steep.steep_out(Z,t))
        for label in ('angular','energy'):
            point['original_forward_'+label+'_Taylor']=point[label+'_Taylor']
            point['original_forward_'+label+'_y_derivative_Taylor']=point[label+'_y_derivative_Taylor']
        point['original_forward_pressure_over_Pstar_squared_Taylor']=point['pressure_over_Pstar_squared_Taylor']
        A=[one/self.heat.k+defects['angular_defect_rows'][0]]+defects['angular_defect_rows'][1:]
        E=[one/self.heat.delta+defects['energy_defect_rows'][0]]+defects['energy_defect_rows'][1:]
        angular=quotient_rows(A,shape['K_rows']); energy=quotient_rows(E,[jet*2 for jet in product_rows(shape['K_rows'],shape['K_rows'])])
        mixed=dict(point['physical_mixed_derivatives_total_order_le4']); fields=dict(point['physical_velocity_and_pressure_y_derivative_Taylor'])
        point['original_forward_pressure_mixed_bounds']=mixed[P]
        mixed[P]={'y'+str(j)+'_Z'+str(n):jet[n]*math.factorial(n) for j,jet in enumerate(pressure) for n in range(5-j)}
        fields[P]=[jet.truncate(4-j) for j,jet in enumerate(pressure)]
        margin=shape['log_K_rate_rows'][0]*(-2)+self.heat.delta
        point.update(angular_Taylor=angular[0],angular_y_derivative_Taylor=[jet.truncate(5-j) for j,jet in enumerate(angular)],
                     energy_Taylor=energy[0],energy_y_derivative_Taylor=[jet.truncate(5-j) for j,jet in enumerate(energy)],
                     pressure_over_Pstar_squared_Taylor=pressure[0],pressure_y_derivative_axial5_Taylor=pressure,
                     physical_mixed_derivatives_total_order_le4=mixed,physical_velocity_and_pressure_y_derivative_Taylor=fields,
                     steep_exit_offset_from_Rtail=q,steep_exit_phase=t,steep_exit_shape=shape,
                     steep_exit_signed_remaining_kernels=kernels,
                     steep_exit_similarity_stress_mixed3_factored={label:{'y'+str(j)+'_Z'+str(n):rows[label][j][n]*math.factorial(n)
                         for j in range(4) for n in range(4-j)} for label in ('theta','axial')},
                     steep_exit_similarity_stress_y_derivative_Taylor={label:[jet.truncate(3-j) for j,jet in enumerate(rows[label])]
                         for label in ('theta','axial','theta_inertial','theta_shear')},
                     steep_exit_future_defect_zeroth_Taylor={label:values[0] for label,values in defects.items()},
                     steep_exit_shear_strength_kappa_minus2=margin,
                     steep_exit_shear_strength_kappa_gt2_certified=endpoints(margin)[0]>0,
                     steep_exit_angular_shear_negative_by_exact_source=True,
                     steep_exit_meridional_moments_and_velocities={label:one*0 for label in ('Mz','Mtheta_z','Ur','Uz')},
                     actual_original_steep_exit_similarity_stress_recovered=True,
                     actual_steep_exit_absolute_pressure_same_source_mixed4_available=True,
                     full_sigma_phi_Gamma_future_integral_used=True,original_angular_energy_pressure_histories_retained=True,
                     steep_exit_waiting_stress_mixed3_join_verified=True,steep_exit_waiting_pressure_mixed4_join_verified=True,
                     source_defined_positive_stress_factors=dict(theta='sqrt(R/2)*B',axial='sqrt(R/2)*B^2',
                         B='Ev0*theta_base*exp(-bh*q)',R='Rtail*exp(q)',q='-wait-1+t',
                         exact_Ev0='Pstar*U*exp(-13/(2*mu)-13); Ev2 is an enclosure only'),
                     steep_exit_physical_axial_viscosity_exact_zero=False,steep_exit_cone_certified=False,
                     global_admissible_stress_lift_constructed=False,physical_energy_integral_certified=False,temporal_recursion=False)
        return point

    def report(self):
        keys=('Z','steep_exit_phase','steep_exit_offset_from_Rtail','steep_exit_shape','steep_exit_signed_remaining_kernels',
              'steep_exit_similarity_stress_mixed3_factored','steep_exit_similarity_stress_y_derivative_Taylor',
              'steep_exit_future_defect_zeroth_Taylor','pressure_y_derivative_axial5_Taylor',
              'steep_exit_shear_strength_kappa_minus2','steep_exit_meridional_moments_and_velocities','source_defined_positive_stress_factors')
        def summary(z,t):
            point=self.steep_out(z,t); return {key:point[key] for key in keys}
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                    scope='Original steep_out t[0,1], Z[-1,1]; complete-future similarity stress mixed3 and pressure mixed4',
                    samples=[summary(z,t) for z,t in (('0','0'),('.5','.5'),('-.5','1'))],whole_steep_exit=summary([-1,1],[0,1]),
                    transport_identities=self.proof,general_moment_stress_identities=self.general_proof,steep_exit_source_bridge=self.bridge,
                    actual_original_steep_exit_similarity_stress_recovered=True,
                    actual_steep_exit_absolute_pressure_same_source_mixed4_available=True,
                    steep_exit_waiting_stress_mixed3_join_verified=True,steep_exit_waiting_pressure_mixed4_join_verified=True,
                    steep_exit_physical_axial_viscosity_exact_zero=False,steep_exit_cone_certified=False,
                    global_admissible_stress_lift_constructed=False,physical_energy_integral_certified=False,
                    temporal_recursion=False,input_hashes=self.hashes)


def run():
    result=CompliantSteepExitStressC3().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual steep exit common full-future stress mixed3, absolute pressure mixed4 and waiting joins generated',flush=True)
    return result


if __name__=='__main__':run()
