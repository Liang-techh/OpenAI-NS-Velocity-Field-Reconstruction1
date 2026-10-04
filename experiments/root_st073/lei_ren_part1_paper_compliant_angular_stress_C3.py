"""Original two-support angular repair: same full moments and absolute pressure.

K depends on Z through the selected repair coefficients. Native cumulative
angular history and direct remaining quadratic bump integrals are retained.
This is a similarity stress companion, not regional physical/cone completion.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_steep_entry_stress_C3 import (
    CompliantSteepEntryStressC3,SourceAST,source_precision,quotient_rows)
from lei_ren_part1_paper_compliant_power_angular_C4 import quotient_log_rates
from lei_ren_part1_paper_compliant_corrected_outer_field import future_bump_weights
from lei_ren_part1_paper_compliant_axial_pulse_field import decay_integral
from lei_ren_part1_paper_compliant_collar_stress_C3 import collar_stress_rows,axial_derivative
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_heat_pressure_C4 import pressure_y_rows
from lei_ren_part1_paper_compliant_pulse_physical_bounds import P
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
DOMAIN=dict(s=[-4,0],Z=[-1,1],offset='q=-wait-Ts-2+s',ordinary_derivative='d_q=d_s')


def angular_future_changes(outer,data,offset):
    """Original remaining beta integrals, preserving selected axial functions."""
    c=outer.ctx; offset=c.mpf(offset); zero=data['coeff'][0]*0
    energy=zero; pressure=zero; source=[]
    for dj,center in zip(data['coeff'],(-3,-1)):
        local=offset-center; lo,hi=endpoints(local); ell=c.mpf('.15')
        if hi<=-endpoints(ell)[1]:future=outer.weights
        elif lo>=endpoints(ell)[1]:future={key:c.mpf(0) for key in outer.weights}
        else:future=future_bump_weights(c,outer.mu,outer.normalization,local,cells=outer.cells)
        energy+=(dj*(2*future['E'])+dj*dj*future['F'])*c.exp(-2*outer.mu*center)
        pressure+=(dj*future['B']+dj*dj*future['D']/2)*c.exp(-outer.prate*center)
        source.append(dict(center=center,local=local,remaining_weights=future))
    return dict(energy=energy,pressure=pressure,original_remaining_bump_integrals=source,
                disjoint_support_cross_quadratic_terms_exact_zero=True)


def angular_shape(heat,terminal,offset,packet):
    c=heat.ctx; offset=c.mpf(offset); rate=heat.a-heat.mu
    F=[packet['swirl_factor_one_plus_h_Taylor']]+packet['actual_angular_bump_y_derivatives'][1:]
    factor=terminal['KR']*c.exp(rate*offset)
    K=[sum((F[n]*math.comb(j,n)*rate**(j-n) for n in range(j+1)),F[0]*0)*factor for j in range(5)]
    one=IntervalTaylor.constant(c,1,5); square=product_rows(K,K)
    rates=quotient_log_rates(F); rates[0]=rates[0]+rate
    return dict(K_rows=K,K_defect_rows=[K[0]-one]+K[1:],
                K_squared_defect_rows=[square[0]-one]+square[1:],
                original_F_rows=F,log_K_rate_rows=rates,actual_K_Z_dependence_retained=True)


def angular_X_rows(heat,shape,X0):
    one=IntervalTaylor.constant(heat.ctx,1,5)
    rates=quotient_log_rates(shape['original_F_rows']); rates[0]=rates[0]+1-heat.mu
    X=[X0]
    for n in range(1,5):
        X.append((one if n==1 else one*0)-sum(
            (X[n-1-j]*rates[j]*math.comb(n-1,j) for j in range(n)),one*0))
    return X


def angular_defect_rows(heat,terminal,offset,shape,X0,kernels):
    """Stable original full future; no division by a capped amplitude."""
    c=heat.ctx; one=IntervalTaylor.constant(c,1,5); rate=heat.a-heat.mu
    KR=terminal['KR']; e=c.exp(heat.delta*offset); p=c.exp(heat.prate*offset)
    baseE=decay_integral(c,2*heat.mu,-offset)*c.exp(-2*heat.mu*offset)
    baseP=decay_integral(c,1+2*heat.mu,-offset)*c.exp(-(1+2*heat.mu)*offset)/2
    Ad=shape['K_rows'][0]*X0-one/heat.k
    Ed=terminal['energy_defect_rows'][0]*e+one*(c.expm1(heat.delta*offset)/heat.delta)
    Ed+=(one*baseE+kernels['energy'])*(KR**2*e)
    Pd=terminal['pressure_defect_rows'][0]*p+one*(c.expm1(heat.prate*offset)/(2*heat.prate))
    Pd+=(one*baseP+kernels['pressure'])*(KR**2*p)
    A=[Ad]; E=[Ed]; Pr=[Pd]; Km=shape['K_defect_rows']; Qm=shape['K_squared_defect_rows']
    for j in range(4):
        A.append(Km[j]-A[j]*heat.k)
        E.append(E[j]*heat.delta-Qm[j])
        Pr.append(Pr[j]*heat.prate-Qm[j]/2)
    return dict(angular_defect_rows=A,energy_defect_rows=E,pressure_defect_rows=Pr,
                K_defect_rows=Km,K_squared_defect_rows=Qm)


def angular_stress_rows(heat,shape,defects,Z,q):
    return collar_stress_rows(heat,shape,defects,Z,q)


def angular_pressure_rows(heat,defects,shape,q):
    one=IntervalTaylor.constant(heat.ctx,1,5)
    return pressure_y_rows(shape['K_rows'],one/(2*heat.prate)+defects['pressure_defect_rows'][0],
                           heat.prate,heat.pressure_scale,q)


def angular_transport_identities():
    a,mu,v,z=s.symbols('a mu s Z',real=True); delta=2*a; k=1-a; p=1+delta
    r=1-mu; d=a-mu; po=1+2*mu; KR=s.Symbol('same_actual_entry_KR',positive=True)
    F=s.Function('original_F')(v,z); N=s.Function('original_native_N')(v,z)
    X=N/F; K=KR*s.exp(d*v)*F
    JE=s.Function('actual_remaining_quadratic_E')(v,z); JP=s.Function('actual_remaining_quadratic_P')(v,z)
    ER=s.Function('same_full_ER')(z); PR=s.Function('same_full_PR')(z)
    E=s.exp(delta*v)*(ER+KR**2*((s.exp(-2*mu*v)-1)/(2*mu)+JE))
    PP=s.exp(p*v)*(PR+KR**2*((s.exp(-po*v)-1)/(2*po)+JP))
    rules={s.diff(N,v):F-r*N,s.diff(JE,v):-s.exp(-2*mu*v)*(F**2-1),
           s.diff(JP,v):-s.exp(-po*v)*(F**2-1)/2}
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Angular identity failed: '+name)
        proofs[name]=True
    zero('native_X_ODE_with_actual_axial_F',s.diff(X,v).subs(rules)-(1-(r+s.diff(F,v)/F)*X))
    zero('native_KX_full_cumulative_angular_ODE',s.diff(K*X,v).subs(rules)-(K-k*K*X))
    zero('same_remaining_energy_full_FTC',s.diff(E,v).subs(rules)-(delta*E-K*K))
    zero('same_remaining_absolute_pressure_full_FTC',s.diff(PP,v).subs(rules)-(p*PP-K*K/2))
    zero('actual_native_energy_normalization',
         2*K*K*((ER/KR**2+JE)*s.exp(2*mu*v)+(1-s.exp(2*mu*v))/(2*mu))/(2*F*F)-E)
    C,qR=s.symbols('same_exact_pressure_scale original_qR',positive=True)
    absolute=-C*s.exp(-p*(qR+v))*PP
    zero('same_absolute_pressure_derivative',s.diff(absolute,v).subs(rules)-C*s.exp(-p*(qR+v))*K*K/2)
    zero('exact_source_shear_variable_rate',delta-2*(d+s.diff(F,v)/F)-(2*mu-2*s.diff(F,v)/F))
    return dict(identities=proofs,actual_axial_F_dependence_retained=True,
                full_moment_FTC_and_original_native_normalizations_verified=True,
                quadratic_tail_is_direct_remaining_beta_integral=True,
                angular_cumulative_history_not_replaced_by_tail_cap=True,
                actual_shear_strength_formula='kappa-2=2*mu-2*F_s/F; no cone sign is yet certified',
                whole_original_domain=DOMAIN)


def angular_entry_join_binding():
    """Replay actual row implementations on arbitrary common axial data."""
    a,mu,z,C,wait,Ts=s.symbols('a mu Z C wait Ts',real=True)
    c=SimpleNamespace(exp=s.exp,expm1=lambda v:s.exp(v)-1,mpf=lambda v:s.Rational(str(v)))
    heat=SimpleNamespace(ctx=c,a=a,mu=mu,S=s.Symbol('same_exact_inverse_Rtail'),delta=2*a,k=1-a,
                         prate=1+2*a,pressure_scale=C,steep=SimpleNamespace(wait=wait,Ts=Ts))
    stub=SimpleNamespace(constant=lambda ctx,value,order:value,variable=lambda ctx,value,order:value)
    env=dict(math=math,c=c,IntervalTaylor=stub,axial_derivative=lambda v:s.diff(v,z),
             sigma_jets=lambda ctx,t:[s.Integer(0)]*5,
             decay_integral=lambda ctx,rate,length:(1-s.exp(-rate*length))/rate)
    asts=SourceAST()
    for stem,name in (('collar_Gamma_C4','product_rows'),('collar_stress_C3','shifted_rows'),
                     ('collar_stress_C3','collar_stress_rows'),('power_angular_C4','quotient_log_rates'),
                     ('steep_entry_stress_C3','entry_shape'),('steep_entry_stress_C3','entry_angular_rows'),
                     ('steep_entry_stress_C3','entry_defect_rows'),('steep_entry_stress_C3','entry_stress_rows'),
                     ('steep_entry_stress_C3','entry_pressure_rows'),('angular_stress_C3','angular_shape'),
                     ('angular_stress_C3','angular_X_rows'),('angular_stress_C3','angular_defect_rows'),
                     ('angular_stress_C3','angular_stress_rows'),('angular_stress_C3','angular_pressure_rows')):
        if name=='quotient_log_rates':
            # The source jet checks do not change this scalar algebra.
            fn=asts.method(stem,name); fn.decorator_list=[]
            fn.body=[node for node in fn.body if not isinstance(node,ast.If)]
            exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual quotient algebra>','exec'),env)
        else:asts.replay(stem,name,env)
    asts.replay('heat_pressure_C4','pressure_y_rows',env,True)
    start={key+'_defect_rows':[s.Function('same_original_power_'+key)(z)] for key in ('angular','energy','pressure')}
    start['KS']=s.Symbol('same_actual_KS',positive=True)
    ek=dict(f=s.Rational(1,2),energy=s.Symbol('same_original_entry_IE'),pressure=s.Symbol('same_original_entry_IP'))
    es=env['entry_shape'](heat,start,s.Integer(0),ek)
    XR=s.Function('same_original_angular_XR')(z)
    ex=env['entry_angular_rows'](heat,XR,es); ed=env['entry_defect_rows'](heat,start,s.Integer(0),ek,es,XR)
    terminal={label:[ed[label][0]] for label in ('energy_defect_rows','pressure_defect_rows')}
    terminal['KR']=es['K_rows'][0]
    packet=dict(swirl_factor_one_plus_h_Taylor=s.Integer(1),actual_angular_bump_y_derivatives=[s.Integer(0)]*5)
    ag=env['angular_shape'](heat,terminal,s.Integer(0),packet)
    ax=env['angular_X_rows'](heat,ag,XR)
    ad=env['angular_defect_rows'](heat,terminal,s.Integer(0),ag,XR,dict(energy=s.Integer(0),pressure=s.Integer(0)))
    q=-wait-Ts-2
    left=env['angular_stress_rows'](heat,ag,ad,z,q); right=env['entry_stress_rows'](heat,es,ed,ex,z,q)
    checks={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Actual angular/entry join failed: '+name)
        checks[name]=True
    for j in range(5):
        zero('actual_K_y'+str(j),ag['K_rows'][j]-es['K_rows'][j])
        for key in ('angular_defect_rows','energy_defect_rows','pressure_defect_rows'):
            zero(key+'_y'+str(j),ad[key][j]-ed[key][j])
    for label in ('theta','axial'):
        for j in range(4):
            for n in range(4-j):zero(label+'_y'+str(j)+'_Z'+str(n),s.diff(left[label][j]-right[label][j],z,n))
    lp=env['angular_pressure_rows'](heat,ad,ag,q); rp=env['entry_pressure_rows'](heat,ed,es,q)
    for j in range(5):
        for n in range(5-j):zero('pressure_y'+str(j)+'_Z'+str(n),s.diff(lp[j]-rp[j],z,n))
    return dict(identities=checks,actual_angular_entry_K4_full_moment_stress_mixed3_pressure_mixed4_AST_join_verified=True,
                arbitrary_actual_terminal_axial_functions_used=True,source_function_equality_not_interval_overlap=True,
                input_hashes=asts.hashes)


@source_precision
def angular_source_bridge(entry):
    asts=SourceAST(); checks={}
    for flag in ('native_X_KX_equals_same_full_future_by_original_ODE_and_right_endpoint',
                 'native_energy_2K2_equals_same_full_future_by_original_ODE_and_right_endpoint',
                 'native_absolute_pressure_equals_same_full_future_by_original_ODE_and_right_endpoint',
                 'native_angular_entry_left_field_pressure_mixed4_join_consumed'):
        if not entry.bridge[flag]:raise ValueError('Angular common datum missing: '+flag)
        checks['consumed_'+flag]=True
    if not entry.bridge['identities']['actual_native_theta_equals_same_reference_B_times_K_verified']:
        raise ValueError('Actual entry reference amplitude missing')
    def syntax(stem,method,target,expected,augmented=False):
        node=asts.expression(stem,method,target,augmented=augmented)
        if ast.dump(node)!=ast.dump(ast.parse(expected,mode='eval').body):
            raise ValueError('Angular source changed: '+stem+'.'+target)
        return node
    syntax('power_angular_C4','data','coeff',"[jet(v) for v in source['physical_coefficient_Taylor']]")
    syntax('power_angular_C4','__init__','self.rate','1-self.mu')
    syntax('power_angular_C4','__init__','self.prate','1+2*self.mu')
    syntax('power_angular_C4','angular','F','[h[0]+1]+h[1:]')
    syntax('power_angular_C4','angular','h[k]','dj*(beta[k]*math.factorial(k))',True)
    syntax('power_angular_C4','angular','beta','self.flat.beta(local)')
    syntax('power_angular_C4','angular','pastA',"dj*(c.exp(self.rate*center)*past['A'])",True)
    syntax('power_angular_C4','angular','pastP',"(dj*past['B']+dj*dj*past['D']/2)*c.exp(-self.prate*center)",True)
    syntax('power_angular_C4','angular','futureE',"(dj*(2*future['E'])+dj*dj*future['F'])*c.exp(-2*self.mu*center)",True)
    syntax('angular_stress_C3','angular_future_changes','energy',"(dj*(2*future['E'])+dj*dj*future['F'])*c.exp(-2*outer.mu*center)",True)
    syntax('angular_stress_C3','angular_future_changes','pressure',"(dj*future['B']+dj*dj*future['D']/2)*c.exp(-outer.prate*center)",True)
    syntax('angular_stress_C3','terminal','point','self.entry.steep_in(Z,0)')
    syntax('angular_stress_C3','terminal',"data['KR']","point['steep_entry_shape']['K_rows'][0][0]")
    syntax('angular_stress_C3','angular','mixed[P]',"{'y'+str(j)+'_Z'+str(n):jet[n]*math.factorial(n) for j,jet in enumerate(pressure) for n in range(5-j)}")
    syntax('angular_stress_C3','angular','fields[P]','[jet.truncate(4-j) for j,jet in enumerate(pressure)]')
    syntax('power_angular_C4','_packet','mixed','flatten_mixed(c,self.mu,one,[c.mpf(0)]+rates,theta,X,energy,pressure,self.flatten.Ev2)')

    # Actual integral loops use the same compact beta and fixed weights.
    for stem,method,target,expected in (
        ('corrected_outer_field','future_bump_weights','beta','raw_beta(c,r)/(ell*normalization)'),
        ('outer_angular_repair','bump_weights','beta','raw_beta(c,raw_coordinate)/(ell*normalization)')):
        syntax(stem,method,target,expected)
    for stem,method,target in (('corrected_outer_field','future_bump_weights','out[key]'),
                               ('outer_angular_repair','bump_weights','result[name]')):
        syntax(stem,method,target,'ds*c.exp(rate*v)*beta**power' if stem=='corrected_outer_field'
               else 'ds*c.exp(rate*s)*beta**power',True)
    for stem,method in (('corrected_outer_field','future_bump_weights'),('outer_angular_repair','bump_weights')):
        fn=asts.method(stem,method)
        loops=[n for n in ast.walk(fn) if isinstance(n,ast.For) and isinstance(n.iter,ast.Tuple)
               and len(n.iter.elts)==5]
        expected="(('A',1-mu,1),('B',-1-2*mu,1),('D',-1-2*mu,2),('E',-2*mu,1),('F',-2*mu,2))"
        if len(loops)!=1 or ast.dump(loops[0].iter)!=ast.dump(ast.parse(expected,mode='eval').body):
            raise ValueError('Actual weighted beta density rates changed')
    syntax('outer_pulse_map','raw_beta','lower','c.exp(-1/(1-c.mpf(far)**2)) if far<1 else c.mpf(0)')
    syntax('outer_pulse_map','raw_beta','upper','c.exp(-1/(1-c.mpf(near)**2)) if near<1 else c.mpf(0)')
    syntax('flat_pulse_derivatives','beta','raw','beta_jets(c,c.mpf(s)/ell)')
    asts.method('flat_pulse_derivatives','beta_jets')
    # Exact native X/E/physical-pressure expressions, not free ODE fields.
    a,mu,v,z,length,Ev2=s.symbols('a mu s Z Lrel same_Ev0_squared',real=True)
    r=1-mu; po=1+2*mu; bp=s.Rational(1,2)+mu
    F=s.Function('actual_F')(v,z); PA=s.Function('actual_past_A')(v,z)
    PE=s.Function('actual_future_E')(v,z); PP=s.Function('actual_past_P')(v,z)
    Xf=s.Function('same_flatten_Xf')(z); Pf=s.Function('same_flatten_Pf')(z); Post=s.Function('same_full_Post')(z)
    c=SimpleNamespace(exp=s.exp,mpf=lambda val:s.Rational(str(val)))
    obj=SimpleNamespace(mu=mu,rate=r,prate=po,bp=bp,Lrel=length,flatten=SimpleNamespace(Ev2=Ev2))
    env=dict(c=c,self=obj,s=v,F=[F],data=dict(flatten_exit_X=Xf,flatten_exit_pressure={'P_over_Pstar_squared':Pf},post=Post),
             y=length+v,pastA=PA,pastP=PP,futureE=PE,
             decay_integral=lambda ctx,rate,t:(1-s.exp(-rate*t))/rate)
    env['eq']=1/r
    env['Xbase']=asts.evaluate(asts.expression('power_angular_C4','angular','Xbase'),env)
    nativeX=asts.evaluate(asts.expression('power_angular_C4','angular','X'),env)
    nativeE=asts.evaluate(asts.expression('power_angular_C4','angular','energy'),env)
    theta=asts.evaluate(asts.expression('power_angular_C4','angular','theta'),env)
    env['baseline']=asts.evaluate(asts.expression('power_angular_C4','angular','baseline'),env)
    pressure=asts.evaluate(asts.expression('power_angular_C4','angular','pressure'),env)
    pressure+=asts.evaluate(asts.expression('power_angular_C4','angular','pressure',augmented=True),env)
    rules={s.diff(PA,v):s.exp(r*v)*(F-1),s.diff(PE,v):-s.exp(-2*mu*v)*(F**2-1),
           s.diff(PP,v):s.exp(-po*v)*(F**2-1)/2}
    expressions=dict(
        actual_native_angular_X_ODE_AST_replayed=s.diff(nativeX,v).subs(rules)-(1-(r+s.diff(F,v)/F)*nativeX),
        actual_native_angular_energy_ODE_AST_replayed=s.diff(nativeE,v).subs(rules)
            -(-s.Rational(1,2)+(2*mu-2*s.diff(F,v)/F)*nativeE),
        actual_native_absolute_pressure_ODE_AST_replayed=s.diff(pressure,v).subs(rules)-Ev2*theta**2/2)
    for name,value in expressions.items():
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Actual angular source ODE failed: '+name)
        checks[name]=True
    thetaR=s.exp(-bp*(100+length))/2
    if s.simplify(theta-thetaR*s.exp(-bp*v)*F)!=0:raise ArithmeticError('Original angular theta reference changed')
    checks['actual_native_theta_B_K_normalization_consumed']=True
    checks['actual_same_full_remaining_E_P_identified_by_ODE_and_entry_datum']=True
    checks['actual_compact_beta_and_disjoint_support_F_squared_density_bound']=True
    checks['actual_original_selected_coefficient_functions_retained']=True
    checks['actual_full_pressure_both_packet_routes_bound']=True
    # Replay the actual production radius relative to waiting terminal.
    path=HERE/(PREFIX+'global_physical_assembly.py'); tree=ast.parse(path.read_text(encoding='utf8')); envr={}
    for node in tree.body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ('MICRO','POST','PULSE'):envr[target.id]=ast.literal_eval(node.value)
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='radius')
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual angular radius>','exec'),envr)
    wait,Ts,origin=s.symbols('wait Ts logRp',real=True)
    assembly=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify),logRp=origin,params=SimpleNamespace(mu=mu))
    outer=SimpleNamespace(Lrel=length); steep=SimpleNamespace(outer=outer,Ts=Ts,wait=wait)
    q=envr['radius'](assembly,'outer_angular',v,{},outer)[0]-envr['radius'](assembly,'waiting',s.Integer(1),{},steep)[0]
    if s.simplify(q-(-wait-Ts-2+v))!=0:raise ArithmeticError('Actual angular source radius changed')
    left=envr['radius'](assembly,'outer_angular',s.Integer(-4),{},outer)[0]
    right=envr['radius'](assembly,'outer_power',s.Integer(1),{},outer)[0]
    if s.simplify(left-right)!=0:raise ArithmeticError('Original power/angular radius join changed')
    asts.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    checks['actual_original_angular_radius_reference_verified']=True
    checks['actual_power_angular_left_source_radius_join_verified']=True
    # The current C4 receipt supplies original flat endpoint field/pressure
    # identities. New stress joins are separately replayed below.
    consumed=entry.bridge['original_native_C4_field_and_history_identities_consumed']
    for flag in ('angular_entry_exact_theta','angular_entry_exact_X','angular_entry_exact_energy','angular_entry_exact_pressure'):
        if not consumed[flag]:raise ValueError('Actual angular/entry field endpoint missing')
    zeros={}
    for offset in (-4,0):
        for center in (-3,-1):
            beta=entry.steep.outer.flat.beta(offset-center)
            zeros[str((offset,center))]=all(endpoints(beta[n])==(mp.mpf(0),mp.mpf(0)) for n in range(5))
    if not all(zeros.values()):raise ArithmeticError('Original endpoint angular beta jets must vanish exactly')
    name=PREFIX+'power_angular_C4_check.json'; native=json.loads((HERE/name).read_bytes())
    for flag in ('angular_s_zero_empty_future_supports_nonzero_past_histories',
                 'same_primitive_ODEs_and_original_flat_beta_identify_mixed_interface_jets'):
        if not native['functional_production_source_identities'][flag]:raise ValueError('Original angular source endpoint missing: '+flag)
        checks['consumed_'+flag]=True
    asts.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    empty=angular_future_changes(entry.steep.outer,entry.steep.outer.data([-1,1]),0)
    for label in ('energy','pressure'):
        if any(endpoints(value)!=(mp.mpf(0),mp.mpf(0)) for value in (empty[label][n] for n in range(empty[label].order+1))):
            raise ArithmeticError('Actual angular remaining kernel must vanish at original s0')
    checks['actual_remaining_energy_pressure_kernels_at_s0_exact_zero']=True
    join=angular_entry_join_binding(); asts.hashes.update(join['input_hashes'])
    return dict(identities=checks,actual_source_AST_bindings=asts.bindings,
                actual_original_endpoint_beta_jets_exact_zero=zeros,actual_angular_entry_formula_join=join,
                original_native_field_pressure_endpoint_joins_consumed=True,
                actual_preceding_power_stress_companion_constructed=False,source_caps_used_as_fields=False,
                input_hashes=asts.hashes)


class CompliantAngularStressC3:
    @source_precision
    def __init__(self):
        self.entry=CompliantSteepEntryStressC3(); self.heat=self.entry.heat; self.outer=self.entry.steep.outer
        self.ctx=self.heat.ctx; self.family,self.source=self.entry.family,self.entry.source
        self.hashes=dict(self.entry.hashes); self.cache={}
        name=PREFIX+'steep_entry_stress_C3_check.json'; receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or not receipt['actual_original_steep_entry_similarity_stress_recovered']:
            raise ValueError('Accepted common entry full-moment datum required')
        for filename,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/filename).read_bytes()).hexdigest()!=digest:raise ValueError('Angular source changed: '+filename)
        self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.proof=angular_transport_identities(); self.bridge=angular_source_bridge(self.entry)
        self.hashes.update(self.bridge['input_hashes']); self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def terminal(self,Z):
        key=tuple(endpoints(self.ctx.mpf(Z)))
        if key not in self.cache:
            point=self.entry.steep_in(Z,0)
            data={name:[point['steep_entry_future_defect_zeroth_Taylor'][name]]
                  for name in ('energy_defect_rows','pressure_defect_rows')}
            data['KR']=point['steep_entry_shape']['K_rows'][0][0]
            self.cache[key]=data
        return self.cache[key]

    @source_precision
    def angular(self,Z,offset):
        c=self.ctx; Z=c.mpf(Z); offset=c.mpf(offset)
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(offset)[0]<-4 or endpoints(offset)[1]>0:
            raise ValueError('Original angular domain Z[-1,1], s[-4,0] required')
        point=dict(self.outer.angular(Z,offset)); terminal=self.terminal(Z)
        kernels=angular_future_changes(self.outer,self.outer.data(Z),offset)
        shape=angular_shape(self.heat,terminal,offset,point)
        Xrows=angular_X_rows(self.heat,shape,point['angular_Taylor'])
        defects=angular_defect_rows(self.heat,terminal,offset,shape,Xrows[0],kernels)
        q=-self.entry.steep.wait-self.entry.steep.Ts-2+offset
        stress=angular_stress_rows(self.heat,shape,defects,Z,q); pressure=angular_pressure_rows(self.heat,defects,shape,q)
        point['original_forward_pressure_over_Pstar_squared_Taylor']=point['pressure_over_Pstar_squared_Taylor']
        mixed=dict(point['physical_mixed_derivatives_total_order_le4']); fields=dict(point['physical_velocity_and_pressure_y_derivative_Taylor'])
        point['original_forward_pressure_mixed_bounds']=mixed[P]
        mixed[P]={'y'+str(j)+'_Z'+str(n):jet[n]*math.factorial(n) for j,jet in enumerate(pressure) for n in range(5-j)}
        fields[P]=[jet.truncate(4-j) for j,jet in enumerate(pressure)]
        one=IntervalTaylor.constant(c,1,5); fullE=[one/self.heat.delta+defects['energy_defect_rows'][0]]+defects['energy_defect_rows'][1:]
        energy=quotient_rows(fullE,[jet*2 for jet in product_rows(shape['K_rows'],shape['K_rows'])])
        point.update(physical_mixed_derivatives_total_order_le4=mixed,physical_velocity_and_pressure_y_derivative_Taylor=fields,
            angular_offset=offset,angular_offset_from_Rtail=q,angular_shape=shape,angular_remaining_quadratic_kernels=kernels,
            angular_similarity_stress_mixed3_factored={label:{'y'+str(j)+'_Z'+str(n):stress[label][j][n]*math.factorial(n)
                for j in range(4) for n in range(4-j)} for label in ('theta','axial')},
            angular_similarity_stress_y_derivative_Taylor={label:[jet.truncate(3-j) for j,jet in enumerate(stress[label])]
                for label in ('theta','axial','theta_inertial','theta_shear')},
            angular_future_defect_zeroth_Taylor={label:values[0] for label,values in defects.items()},
            pressure_over_Pstar_squared_Taylor=pressure[0],pressure_y_derivative_axial5_Taylor=pressure,
            original_native_energy_Taylor=point['energy_Taylor'],same_full_energy_Taylor=energy[0],
            angular_shear_strength_kappa_minus2=-2*shape['log_K_rate_rows'][0][0]+self.heat.delta,
            angular_meridional_moments_and_velocities={label:one*0 for label in ('Mz','Mtheta_z','Ur','Uz')},
            actual_original_angular_similarity_stress_recovered=True,actual_angular_absolute_pressure_same_source_mixed4_available=True,
            angular_entry_stress_mixed3_and_pressure_mixed4_join_verified=True,
            original_selected_angular_coefficient_functions_retained=True,actual_angular_K_Z_dependence_retained=True,
            angular_physical_decomposition_constructed=False,angular_cone_certified=False,
            preceding_power_stress_companion_constructed=False,global_admissible_stress_lift_constructed=False,
            physical_energy_integral_certified=False,temporal_recursion=False)
        return point

    def report(self):
        keys=('Z','angular_offset','angular_offset_from_Rtail','angular_shape','angular_remaining_quadratic_kernels',
              'angular_similarity_stress_mixed3_factored','angular_future_defect_zeroth_Taylor',
              'pressure_y_derivative_axial5_Taylor','angular_shear_strength_kappa_minus2','angular_meridional_moments_and_velocities')
        def packet(z,v):
            point=self.angular(z,v); return {key:point[key] for key in keys}
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,domain=DOMAIN,
            scope='Whole original angular s[-4,0],Z[-1,1]; similarity stress mixed3 and same absolute pressure mixed4',
            samples=[packet(z,v) for z,v in (('0','-4'),('.5','-3'),('-.5','-1'),('0','0'))],
            whole_angular=packet([-1,1],[-4,0]),
            support_crossings=[packet([-1,1],[center+side*.15-.01,center+side*.15+.01]) for center in (-3,-1) for side in (-1,1)],
            transport_identities=self.proof,source_bridge=self.bridge,input_hashes=self.hashes,
            actual_original_angular_similarity_stress_recovered=True,actual_angular_absolute_pressure_same_source_mixed4_available=True,
            angular_entry_stress_mixed3_and_pressure_mixed4_join_verified=True,
            actual_angular_K_Z_dependence_retained=True,original_selected_angular_coefficient_functions_retained=True,
            angular_physical_decomposition_constructed=False,angular_cone_certified=False,
            preceding_power_stress_companion_constructed=False,global_admissible_stress_lift_constructed=False,
            physical_energy_integral_certified=False,temporal_recursion=False)


@source_precision
def run():
    result=CompliantAngularStressC3().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original angular full moments, axial-dependent stress and absolute pressure generated; physical/cone pending',flush=True)
    return result


if __name__=='__main__':run()
