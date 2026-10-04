"""Whole original sigmoid-entry two-vector cone with a correlated ODE bound.

The native angular history is unchanged. Its small positive inlet gap is
retained before subtracting unit baselines; a continuous ODE comparison
covers the entire entry. The completed tensor/global cone is not admitted.
"""
import ast
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_steep_entry_stress_C3 import (
    CompliantSteepEntryStressC3,SourceAST,source_precision)
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_collar_physical_C2 import stress_source_log_parts
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
DOMAIN=dict(t=[0,1],Z=[-1,1],offset='q=-wait-Ts-2+t',ordinary_derivative='d_q=d_t')
TAIL_DOMAIN='Rtail*exp(-wait-Ts-2)<=R<Rtail*exp(3); Gamma stress exactly zero beyond'


def steep_entry_cone_identities():
    a,mu,t,z,wait,Ts,eps=s.symbols('a mu t Z wait Ts eps',real=True)
    delta=2*a; k=1-a; bh=s.Rational(1,2)+a; r=1-mu; b=(1-delta)/2
    sigma=s.Function('original_sigma')(t); h=r*(1-sigma); g=a-mu-r*sigma
    X=s.Function('native_X')(t,z); M=k*X-b*z*s.diff(X,z)-1
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Entry cone identity failed: '+name)
        proofs[name]=True
    zero('actual_variable_kappa_minus2',delta-2*g-(2*mu+2*r*sigma))
    zero('continuous_native_M_ODE',s.diff(M,t).subs(
        {s.diff(X,t,z):-h*s.diff(X,z),s.diff(X,t):1-h*X},simultaneous=True)+h*M-(mu-a+r*sigma))
    xf=s.Function('actual_flatten_exit_X')(z); past=s.Function('actual_full_angular_repair')(z)
    length=s.Symbol('original_Lrel',positive=True)
    XR=1/r+(xf-1/r)*s.exp(-r*length)+past
    stable=(mu-a)/r+s.exp(-r*length)*(k*(xf-1/r)-b*z*s.diff(xf,z))+k*past-b*z*s.diff(past,z)
    zero('same_original_inlet_M_without_unit_subtraction',k*XR-b*z*s.diff(XR,z)-1-stable)
    gap=s.Symbol('positive_inlet_M_lower',positive=True)
    lower=gap*s.exp(-r*t)
    zero('continuous_comparison_subsolution',s.diff(lower,t)+h*lower+r*sigma*lower)
    q=-wait-Ts-2+t; K,Sin=s.symbols('actual_K exact_positive_inverse_Rtail',positive=True)
    shear=2*Sin*s.exp(-q)*(g-(1+a))*K
    zero('strict_original_shear_negative_positive_source_factors',shear+2*Sin*s.exp(-q)*(1+mu+r*sigma)*K)
    zero('source_shear_rate_upper_equals_two',1+mu+r-2)
    base=s.Symbol('same_Bbase',positive=True); B=base*s.exp(-bh*q)
    Bmax=base*s.exp(bh*(wait+Ts+2))
    zero('whole_entry_B_decreases',B-Bmax*s.exp(-bh*t))
    zero('whole_entry_B_derivative',s.diff(B,t)+bh*B)
    Ev0,thetaT=s.symbols('same_Ev0 same_thetaT',positive=True)
    zero('actual_inlet_Bmax',Bmax.subs(base,Ev0*thetaT*s.exp(-bh*wait)/(1-eps))
         -Ev0*thetaT*s.exp(bh*(Ts+2))/(1-eps))
    zero('shared_Ts_log_cancellation',-s.Rational(3,2)*Ts+bh*(Ts+2)-(-k*Ts+2*bh))
    qt,ct,cz,st,Bcurrent,m=s.symbols('Qtheta Ctheta Cz Stheta B kappa_minus2',real=True)
    zero('normalized_original_two_vector_cone',
         2*(qt*ct*qt*st)**2-m*(qt**2*Bcurrent*cz*st)**2
         -qt**4*st**2*(2*ct**2-m*Bcurrent**2*cz**2))
    return dict(identities=proofs,whole_original_domain=DOMAIN,
                native_M_definition='k*X-b*Z*X_Z-1',
                continuous_comparison='M>=M0_lower*exp(-(1-mu)) on the whole t[0,1]',
                comparison_requires_mu_above_a_and_positive_inlet_M=True,
                actual_variable_kappa_minus2_range='[2*mu,2]',
                sufficient_directional_condition='2*Bmax^2*(Cz_abs_upper/Ctheta_lower)^2<2',
                exact_S_positive_inverse_Rtail_not_cap_endpoint=True,
                same_positive_physical_nu_lambda_factors_preserve_two_vector_cone=True,
                completed_diagonal_does_not_enter_two_vector_cone=True)


def steep_entry_Bmax_binding():
    tree=ast.parse((HERE/(PREFIX+'collar_physical_C2.py')).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='stress_source_log_parts')
    env={}; exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual entry B logs>','exec'),env)
    a,mu,wait,Ts,length,lp,lu,lrp,lone=s.symbols('a mu wait Ts Lrel logP logU logRp logone',real=True)
    c=SimpleNamespace(mpf=s.sympify,ln=s.log)
    flatten=SimpleNamespace(logEv2_parts=dict(inlet_log=2*lu))
    assembly=SimpleNamespace(ctx=c,logP=lp,logRp=lrp,dispatch=SimpleNamespace(provider=lambda chart:flatten))
    heat=SimpleNamespace(ctx=c,a=a,mu=mu,k=1-a,bh=s.Rational(1,2)+a,
                         steep=SimpleNamespace(wait=wait,Ts=Ts,logone=lone,outer=SimpleNamespace(Lrel=length)))
    original=env['stress_source_log_parts'](assembly,heat,-wait-Ts-2)['B']
    tree=ast.parse(Path(__file__).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='whole_steep_entry_bounds')
    names=('logparts',"logparts['steep']","logparts['waiting_and_current']")
    assignments=[n for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t) in names for t in n.targets)]
    if len(assignments)!=3:raise ValueError('Owned entry Bmax cancellation changed')
    env.update(assembly=assembly,h=heat,c=c)
    exec(compile(ast.Module(body=assignments,type_ignores=[]),'<owned entry B logs>','exec'),env)
    current=env['logparts']; checks={}
    if set(original)!=set(current):raise ValueError('Entry B source log keys changed')
    for key in set(original)-{'steep','waiting_and_current'}:
        if s.simplify(original[key]-current[key])!=0:raise ArithmeticError('Entry B source mismatch: '+key)
        checks[key]=True
    if s.simplify(sum(original.values())-sum(current.values()))!=0:raise ArithmeticError('Entry original Ts/wait cancellation failed')
    checks['same_original_Ts_wait_cancellation']=True
    return dict(identities=checks,actual_Bmax_is_original_source_at_q_minus_wait_minus_Ts_minus2=True,
                source_cap_values_used_as_field=False)


@source_precision
def steep_entry_cone_source_bridge(stress):
    asts=SourceAST(); checks={}
    for flag in ('actual_entry_f_equals_half_minus_original_J_verified',
                 'actual_original_entry_production_radius_reference_verified',
                 'actual_complete_power_inlet_moments_and_absolute_datum_consumed',
                 'native_X_KX_equals_same_full_future_by_original_ODE_and_right_endpoint',
                 'native_energy_2K2_equals_same_full_future_by_original_ODE_and_right_endpoint',
                 'native_absolute_pressure_equals_same_full_future_by_original_ODE_and_right_endpoint'):
        if not stress.bridge[flag]:raise ValueError('Entry cone native history missing: '+flag)
        checks['consumed_'+flag]=True
    if not stress.bridge['actual_entry_power_formula_join']['actual_entry_power_K4_stress_mixed3_pressure_mixed4_AST_join_verified']:
        raise ValueError('Entry cone actual right join missing')
    exitbridge=stress.power_source.exit_source.bridge
    if not exitbridge['exact_S_is_positive_source_not_enclosure_endpoint']:raise ValueError('Entry cone exact S missing')
    if not exitbridge['actual_original_velocity_reference_normalization_verified']:raise ValueError('Entry cone source normalization missing')
    if not stress.bridge['identities']['actual_native_angular_source_ODE_replayed']:
        raise ValueError('Actual native entry angular ODE required for continuous comparison')
    if not exitbridge['original_sigma_between0and1_reflection_and_flat_endpoints_bound']:
        raise ValueError('Original sigmoid source bounds required')
    checks['consumed_actual_native_angular_source_ODE_replayed']=True
    checks['consumed_original_sigma_between0and1_reflection_and_flat_endpoints_bound']=True
    logpath=HERE/(PREFIX+'collar_physical_C2.py')
    asts.hashes[logpath.name]=hashlib.sha256(logpath.read_bytes()).hexdigest()

    def syntax(stem,method,target,expected,augmented=False):
        node=asts.expression(stem,method,target,augmented=augmented)
        if ast.dump(node)!=ast.dump(ast.parse(expected,mode='eval').body):
            raise ValueError('Entry cone source changed: '+stem+'.'+target)
        return node
    syntax('steep_entry_stress_C3','entry_shape','sig','sigma_jets(c,t)')
    syntax('steep_entry_stress_C3','entry_angular_rows','h',"[((1-heat.mu)*(1-shape['original_sigma_jets'][0]))]+g[1:]")
    angular=asts.method('steep_entry_stress_C3','entry_angular_rows')
    append=[node.args[0] for node in ast.walk(angular) if isinstance(node,ast.Call)
            and isinstance(node.func,ast.Attribute) and ast.unparse(node.func)=='X.append']
    if len(append)!=1:raise ValueError('Actual entry angular recurrence changed')
    import math
    mu0,sig0,X0=s.symbols('source_mu original_sigma native_X0',real=True)
    derivative=asts.evaluate(append[0],dict(n=1,one=s.Integer(1),X=[X0],h=[(1-mu0)*(1-sig0)],math=math))
    if s.simplify(derivative-(1-(1-mu0)*(1-sig0)*X0))!=0:
        raise ArithmeticError('Actual entry X derivative differs from comparison ODE')
    checks['actual_entry_angular_first_derivative_AST_replayed']=True
    checks['actual_original_sigma_source_is_Z_independent']=True
    syntax('power_angular_C4','__init__','self.rate','1-self.mu')
    syntax('power_angular_C4','data','f','self.flatten.flatten(Z,100)')
    syntax('power_angular_C4','data','coeff',"[jet(v) for v in source['physical_coefficient_Taylor']]")
    syntax('power_angular_C4','__init__','self.weights',"{k:read_interval(c,v) for k,v in angular['bump_weights'].items()}")
    value=asts.expression('power_angular_C4','data','value')
    if not isinstance(value,ast.Call) or not isinstance(value.func,ast.Name) or value.func.id!='dict':
        raise ValueError('Actual angular terminal data constructor changed')
    items={keyword.arg:keyword.value for keyword in value.keywords}
    for key,expected in (('coeff','coeff'),('flatten_exit_X',"f['angular_Taylor']")):
        if ast.dump(items[key])!=ast.dump(ast.parse(expected,mode='eval').body):
            raise ValueError('Actual angular terminal source data changed: '+key)
    syntax('steep_waiting_C4','data','terminal','self.outer.angular(Z,0)')
    syntax('steep_waiting_C4','data','XR',"terminal['angular_Taylor']")
    native_entry=syntax('steep_waiting_C4','steep_in','X',"(data['XR']+kernels['angular'])*c.exp(-self.rate*(t-J))")

    syntax('power_angular_C4','angular','local','s-center')
    syntax('power_angular_C4','angular','beta','self.flat.beta(local)')
    syntax('power_angular_C4','angular','h[k]','dj*(beta[k]*math.factorial(k))',True)
    density=syntax('power_angular_C4','angular','pastA',"dj*(c.exp(self.rate*center)*past['A'])",True)
    syntax('power_angular_C4','angular','F','[h[0]+1]+h[1:]')
    syntax('power_angular_C4','angular','y','self.Lrel+s')
    eq=syntax('power_angular_C4','angular','eq','IntervalTaylor.constant(c,1/self.rate,5)')
    xbase=syntax('power_angular_C4','angular','Xbase',"eq+(data['flatten_exit_X']-1/self.rate)*c.exp(-self.rate*y)")
    xnative=syntax('power_angular_C4','angular','X','(Xbase+pastA*c.exp(-self.rate*s))/F[0]')
    fn=asts.method('power_angular_C4','angular')
    loops=[n for n in ast.walk(fn) if isinstance(n,ast.For) and ast.unparse(n.target)=='(dj, center)']
    if len(loops)!=1 or ast.dump(loops[0].iter)!=ast.dump(ast.parse("zip(data['coeff'],(-3,-1))",mode='eval').body):
        raise ValueError('Actual angular repair supports changed')
    full=[n for n in ast.walk(fn) if isinstance(n,ast.If) and ast.unparse(n.test)=='lo >= endpoints(ell)[1]']
    if len(full)!=1 or ast.unparse(full[0].body[0])!='past = self.weights':
        raise ValueError('Actual terminal angular full weights branch changed')
    beta=asts.method('flat_pulse_derivatives','beta_jets')
    guard=next(n for n in ast.walk(beta) if isinstance(n,ast.If))
    if ast.unparse(guard.test)!='hi <= -1 or lo >= 1' or ast.unparse(guard.body[0])!='return IntervalTaylor.constant(c, 0, 4)':
        raise ValueError('Actual compact beta zero-support guard changed')
    syntax('flat_pulse_derivatives','beta','raw','beta_jets(c,c.mpf(s)/ell)')
    for coordinate in (1,3):
        jets=stress.steep.outer.flat.beta(coordinate)
        if any(endpoints(jets[n])!=(mp.mpf(0),mp.mpf(0)) for n in range(5)):
            raise ArithmeticError('Original angular terminal beta must be exactly zero')
    a,mu,z,length=s.symbols('a mu Z Lrel',real=True); r=1-mu; k=1-a; b=(1-2*a)/2
    xf=s.Function('actual_flatten_exit_X')(z); dj=[s.Function('actual_d'+str(i))(z) for i in range(2)]
    weight=s.Symbol('actual_full_bump_A',positive=True); c=SimpleNamespace(exp=s.exp)
    obj=SimpleNamespace(rate=r,Lrel=length)
    env=dict(self=obj,c=c,IntervalTaylor=SimpleNamespace(constant=lambda c,v,n:v))
    equilibrium=asts.evaluate(eq,env)
    past=sum(asts.evaluate(density,dict(env,dj=d,center=center,past={'A':weight})) for d,center in zip(dj,(-3,-1)))
    base=asts.evaluate(xbase,dict(env,eq=equilibrium,data={'flatten_exit_X':xf},y=length))
    native=asts.evaluate(xnative,dict(env,Xbase=base,pastA=past,s=s.Integer(0),F=[s.Integer(1)]))
    stable=(mu-a)/r+s.exp(-r*length)*(k*(xf-1/r)-b*z*s.diff(xf,z))+k*past-b*z*s.diff(past,z)
    if s.simplify(k*native-b*z*s.diff(native,z)-1-stable)!=0:
        raise ArithmeticError('Actual angular terminal stable M source differs')
    inlet=asts.evaluate(native_entry,dict(env,data={'XR':native},kernels={'angular':s.Integer(0)},t=s.Integer(0),J=s.Integer(0)))
    if s.simplify(inlet-native)!=0:raise ArithmeticError('Original entry inlet must retain actual XR')
    checks['actual_native_entry_inlet_retains_same_angular_XR']=True
    checks['actual_angular_terminal_M_stable_formula_AST_verified']=True
    checks['actual_flatten_exit_X_and_selected_angular_coefficients_retained']=True
    checks['actual_full_bump_weights_and_compact_zero_terminal_support_verified']=True
    checks['actual_entry_power_similarity_join_consumed']=True
    checks['exact_S_positive_source_consumed']=True
    return dict(identities=checks,actual_source_AST_bindings=asts.bindings,input_hashes=asts.hashes,
                actual_inlet_M_recovered_without_unit_cancellation=True,
                actual_Bmax_source_log_binding=steep_entry_Bmax_binding(),
                source_caps_used_as_fields=False,phase_samples_used_as_proof=False)


@source_precision
def whole_steep_entry_bounds(stress,assembly):
    h=stress.heat; c=stress.ctx; out=stress.steep.outer
    packet=stress.steep_in([-1,1],[0,1]); rows=packet['steep_entry_similarity_stress_mixed3_factored']
    rawCt=rows['theta']['y0_Z0']; Cz=rows['axial']['y0_Z0']
    data=out.data([-1,1]); Xf=data['flatten_exit_X']; r=1-h.mu; b=(1-h.delta)/2
    # Actual native terminal XR is equilibrium plus the same attenuated
    # flatten history and both selected repair moments. Keep its small
    # positive equilibrium gap instead of subtracting rounded unit boxes.
    equilibrium_gap=(h.mu-h.a)/r
    flattened_correction=c.exp(-r*out.Lrel)*(h.k*abs(Xf[0]-1/r)+abs(b)*abs(Xf[1]))
    repairs=[]
    for dj,center in zip(data['coeff'],(-3,-1)):
        repairs.append(c.exp(r*center)*out.weights['A']*(h.k*abs(dj[0])+abs(b)*abs(dj[1])))
    correction=flattened_correction+sum(repairs,c.mpf(0))
    inlet_M_lower=c.mpf(endpoints(equilibrium_gap-correction)[0])
    M_lower=c.mpf(endpoints(inlet_M_lower*c.exp(-r))[0])
    # This is an upper bound on the exact shear coefficient; S>0 is
    # established from its source definition, not the zero cap endpoint.
    shear_upper=c.mpf(endpoints(4*h.S*c.exp(stress.steep.wait+stress.steep.Ts+2))[1])
    bracket_lower=c.mpf(endpoints(M_lower-shear_upper)[0])
    K=packet['steep_entry_shape']['K_rows'][0][0]
    K_lower=c.mpf(endpoints(K)[0]); theta_lower=c.mpf(endpoints(K_lower*bracket_lower)[0])
    axial_upper=c.mpf(endpoints(abs(Cz))[1])
    logparts=stress_source_log_parts(assembly,h,c.mpf(0))['B']
    logparts['steep']=-h.k*h.steep.Ts
    logparts['waiting_and_current']=2*h.bh
    logBmax=sum(logparts.values(),c.mpf(0))
    m=c.mpf([endpoints(2*h.mu)[0],mp.mpf(2)])
    sigma=packet['steep_entry_shape']['original_sigma_jets'][0]
    if endpoints(sigma)[0]<0 or endpoints(sigma)[1]>1:
        raise ArithmeticError('Original whole-entry sigmoid must lie in[0,1]')

    margins=dict(mu=h.mu,mu_minus_a=h.mu-h.a,rate=r,delta=h.delta,epsilon=h.eps,
                 one_minus_epsilon=1-h.eps,k=h.k,p=h.prate,bh=h.bh,Lmin=1-h.delta,
                 wait=h.steep.wait,Ts=h.steep.Ts,source_K_lower=K_lower,
                 inlet_M=inlet_M_lower,continuous_M=M_lower,theta_bracket=bracket_lower,theta_stress=theta_lower)
    for label,value in margins.items():
        if endpoints(value)[0]<=0 or not mp.isfinite(endpoints(value)[1]):raise ArithmeticError('Entry cone positive margin failed: '+label)
    if endpoints(h.mu)[1]>=1 or endpoints(h.eps)[1]>=1:raise ArithmeticError('Original entry parameters outside source domain')
    if endpoints(h.delta)[1]>=1:raise ArithmeticError('Original entry L must be at most one and positive')
    if endpoints(out.weights['A'])[0]<=0:raise ArithmeticError('Actual original full bump A weight positive required')
    if endpoints(axial_upper)[0]<=0 or not mp.isfinite(endpoints(axial_upper)[1]):raise ArithmeticError('Finite nonzero entry axial bound required')
    if not all(mp.isfinite(v) for v in endpoints(logBmax)):raise ArithmeticError('Finite actual inlet B logs required')
    logcone=c.ln(2)+2*logBmax+2*c.ln(axial_upper/theta_lower); threshold=c.ln(2)
    if endpoints(logcone)[1]>=endpoints(threshold)[0]:raise ArithmeticError('Whole original entry directional cone failed')
    return dict(positive_margins=margins,actual_whole_source_raw_theta_enclosure=rawCt,actual_whole_source_axial_enclosure=Cz,
                raw_whole_theta_box_used_as_positivity_proof=False,
                actual_equilibrium_M_gap=equilibrium_gap,actual_flatten_history_correction_absolute_bound=flattened_correction,
                actual_selected_angular_repair_correction_absolute_bounds=repairs,actual_total_correction_absolute_bound=correction,
                actual_inlet_M_uniform_lower=inlet_M_lower,continuous_M_uniform_lower=M_lower,
                exact_source_shear_coefficient_absolute_upper=shear_upper,theta_uniform_positive_lower=theta_lower,
                axial_uniform_absolute_upper=axial_upper,actual_kappa_minus2_enclosure=m,exact_kappa_minus2_upper=c.mpf(2),
                actual_Bmax_source_log_parts=logparts,actual_log_Bmax_enclosure=logBmax,
                log_directional_term_upper=logcone,log_threshold=threshold,
                whole_t_Z_domain_and_original_Ts_used=True,continuous_M_ODE_comparison_used=True,
                phase_grid_sampling_used_as_proof=False,source_caps_only_define_bounds_not_field_values=True,
                interval_shear_upper_zero_not_used_for_strict_source_sign=True)


def run():
    with mp.workdps(300):
        stress=CompliantSteepEntryStressC3(); assembly=CompliantGlobalPhysicalAssembly()
        if (stress.family,stress.source)!=(assembly.family,assembly.source):raise ValueError('Entry cone field family differs')
        hashes=dict(stress.hashes); hashes.update(assembly.hashes)
        for stem,gates in (
            ('steep_entry_stress_C3_check',('actual_original_steep_entry_similarity_stress_recovered',)),
            ('steep_entry_physical_C2_check',('actual_regional_physical_steep_entry_stress_remainder_identity_verified',
                'steep_entry_power_physical_stress_mixed3_and_remainder_mixed2_join_verified')),
            ('steep_power_cone_check',('steep_power_exit_waiting_collar_two_vector_cone_certified',))):
            name=PREFIX+stem+'.json'; receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or any(not receipt[gate] for gate in gates):raise ValueError('Entry cone prerequisite missing: '+name)
            if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(stress.family,stress.source):
                raise ValueError('Entry cone prerequisite family differs')
            for source,digest in receipt['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Entry cone prerequisite changed: '+source)
            hashes.update(receipt['input_hashes']); hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        print('Original entry source and physical/power prerequisites admitted',flush=True)
        proof=steep_entry_cone_identities(); bridge=steep_entry_cone_source_bridge(stress)
        bounds=whole_steep_entry_bounds(stress,assembly); hashes.update(bridge['input_hashes']); hashes.update(assembly.dispatch.hashes)
        hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(actual_five_defect_family_sha256=stress.family,implicit_source_sha256=stress.source,
                    scope='Whole original entry t[0,1], source Z[-1,1]; physical tau>0, r>0, |Z|<1',
                    domain=DOMAIN,tail_domain=TAIL_DOMAIN,exact_identities=proof,source_bridge=bridge,bounds=bounds,input_hashes=hashes,
                    steep_entry_cone_certified=True,actual_steep_entry_theta_stress_positive=True,
                    actual_steep_entry_source_shear_strictly_negative=True,actual_steep_entry_directional_cone_uniform_margin=True,
                    steep_entry_power_exit_waiting_collar_two_vector_cone_certified=True,
                    same_source_steep_entry_power_similarity_and_physical_joins_consumed=True,
                    fixed_positive_viscosity_two_vector_cone_transfer_verified=True,
                    steep_entry_regional_remainder_exact_zero=False,upstream_angular_stress_companion_constructed=False,
                    completed_full_tensor_cone_certified=False,whole_outer_cone_certified=False,
                    global_admissible_stress_lift_constructed=False,independently_bounded_global_flat_remainder=False,
                    physical_energy_integral_certified=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
        print('Whole original entry cone generated by continuous native M comparison; power tail joined, upstream/global pending',flush=True)
        return result


class CertifiedSteepEntryPhysical:
    """Attach the accepted original two-vector cone to the unchanged field."""
    def __init__(self):
        from lei_ren_part1_paper_compliant_steep_entry_physical_C2 import CompliantSteepEntryPhysicalC2
        path=HERE/(PREFIX+'steep_entry_cone_check.json'); receipt=json.loads(path.read_bytes())
        if not receipt['all_passed'] or not receipt['steep_entry_cone_certified']:raise ValueError('Accepted entry cone required')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Certified entry source changed: '+source)
        self.field=CompliantSteepEntryPhysicalC2(); self.family,self.source=self.field.family,self.field.source
        if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(self.family,self.source):
            raise ValueError('Entry physical cone family differs')
        self.receipt=receipt

    def steep_in(self,*args,**kwargs):
        point=self.field.steep_in(*args,**kwargs)
        point.update(steep_entry_cone_certified=True,steep_entry_power_exit_waiting_collar_two_vector_cone_certified=True,
                     certified_cone_scope='Original two-vector only; completed diagonal/global cone not certified',
                     completed_full_tensor_cone_certified=False)
        return point


if __name__=='__main__':run()
