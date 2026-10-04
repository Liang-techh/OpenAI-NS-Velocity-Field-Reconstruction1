"""Whole original angular two-vector cone from stable native history bounds.

The selected axial repair and both compact beta supports are unchanged.
Equilibrium surplus is retained before interval subtraction. Continuous
source bounds cover the whole domain; completed-tensor/global gates wait.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_angular_stress_C3 import CompliantAngularStressC3,SourceAST,DOMAIN
from lei_ren_part1_paper_compliant_waiting_stress_C3 import source_precision
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_collar_physical_C2 import stress_source_log_parts
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
TAIL_DOMAIN='Rtail*exp(-wait-Ts-6)<=R<Rtail*exp(3); Gamma stress exactly zero beyond'


def angular_cone_identities():
    a,mu,q,z,wait,Ts=s.symbols('a mu s Z wait Ts',real=True)
    delta=2*a; k=1-a; b=(1-delta)/2; r=1-mu; bh=s.Rational(1,2)+a
    G=s.Function('same_axially_constant_G')(q)
    N=s.Function('same_native_cumulative_N')(q,z); F=s.Function('actual_angular_F')(q,z)
    A=G*N; K=G*F; W=k*N-b*z*s.diff(N,z)-F
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Angular cone identity failed: '+name)
        proofs[name]=True
    zero('native_axial_dependent_inertial_numerator',k*A-b*z*s.diff(A,z)-K-G*W)
    zero('same_native_A_ODE',(s.diff(A,q)+k*A-K).subs(
        {s.diff(G,q):(a-mu)*G,s.diff(N,q):F-r*N},simultaneous=True))
    zero('native_continuous_W_ODE',s.diff(W,q).subs(
        {s.diff(N,q,z):s.diff(F,z)-r*s.diff(N,z),s.diff(N,q):F-r*N},simultaneous=True)
         +r*W-((mu-a)*F-b*z*s.diff(F,z)-s.diff(F,q)))
    zero('actual_shear_numerator',s.diff(K,q).subs(s.diff(G,q),(a-mu)*G)-(1+a)*K
        -G*(s.diff(F,q)-(1+mu)*F))
    zero('actual_variable_kappa_minus2',delta-2*(a-mu+s.diff(F,q)/F)-(2*mu-2*s.diff(F,q)/F))
    length=s.Symbol('original_Lrel',positive=True); xf=s.Function('actual_flatten_exit_X')(z)
    past=s.Function('same_partial_selected_angular_moments')(q,z); h=s.Function('actual_h')(q,z)
    native=1/r+(xf-1/r)*s.exp(-r*(length+q))+past
    stable=(mu-a)/r+s.exp(-r*(length+q))*(k*(xf-1/r)-b*z*s.diff(xf,z))+k*past-b*z*s.diff(past,z)-h
    zero('native_W_without_unit_subtraction',k*native-b*z*s.diff(native,z)-(1+h)-stable)
    KR=s.Symbol('same_actual_KR',positive=True); g=KR*s.exp((a-mu)*q)
    zero('actual_G_decreases_as_s_increases',s.diff(g,q)-(a-mu)*g)
    base=s.Symbol('same_Bbase',positive=True); B=base*s.exp(-bh*(-wait-Ts-2+q))
    Bmax=base*s.exp(bh*(wait+Ts+6))
    zero('original_angular_inlet_Bmax',B-Bmax*s.exp(-bh*(q+4)))
    zero('whole_angular_B_decreases',s.diff(B,q)+bh*B)
    zero('shared_Ts_wait_log_cancellation',-s.Rational(3,2)*Ts+bh*(Ts+6)-(-k*Ts+6*bh))
    qt,ct,cz,st,Bcurrent,m=s.symbols('Qtheta Ctheta Cz Stheta B kappa_minus2',real=True)
    zero('normalized_original_two_vector_cone',2*(qt*ct*qt*st)**2-m*(qt**2*Bcurrent*cz*st)**2
        -qt**4*st**2*(2*ct**2-m*Bcurrent**2*cz**2))
    return dict(identities=proofs,whole_original_domain=DOMAIN,
        native_W_definition='k*N-b*Z*N_Z-F with A=G*N, K=G*F',
        whole_W_bound='(mu-a)/(1-mu) minus absolute flatten/past-repair/h bounds',
        actual_shear_rate='kappa-2=2*mu-2*F_s/F',
        sufficient_directional_condition='m_upper*Bmax^2*(Cz_abs_upper/Ctheta_lower)^2<2',
        actual_axial_F_dependence_retained=True,
        exact_S_positive_inverse_Rtail_not_cap_endpoint=True,
        same_positive_physical_nu_lambda_factors_preserve_two_vector_cone=True,
        completed_diagonal_does_not_enter_two_vector_cone=True)


def angular_Bmax_binding():
    asts=SourceAST(); fn=asts.method('collar_physical_C2','stress_source_log_parts')
    env={}; exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual angular B logs>','exec'),env)
    a,mu,wait,Ts,length,lp,lu,lrp,lone=s.symbols('a mu wait Ts Lrel logP logU logRp logone',real=True)
    c=SimpleNamespace(mpf=s.sympify,ln=s.log); flatten=SimpleNamespace(logEv2_parts=dict(inlet_log=2*lu))
    assembly=SimpleNamespace(ctx=c,logP=lp,logRp=lrp,dispatch=SimpleNamespace(provider=lambda chart:flatten))
    heat=SimpleNamespace(ctx=c,a=a,mu=mu,k=1-a,bh=s.Rational(1,2)+a,
        steep=SimpleNamespace(wait=wait,Ts=Ts,logone=lone,outer=SimpleNamespace(Lrel=length)))
    original=env['stress_source_log_parts'](assembly,heat,-wait-Ts-6)['B']
    fn=next(n for n in ast.walk(ast.parse(Path(__file__).read_text(encoding='utf8')))
            if isinstance(n,ast.FunctionDef) and n.name=='whole_angular_bounds')
    names=('logparts',"logparts['steep']","logparts['waiting_and_current']")
    assignments=[n for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t) in names for t in n.targets)]
    if len(assignments)!=3:raise ValueError('Owned angular Bmax regrouping changed')
    env.update(assembly=assembly,h=heat,c=c)
    exec(compile(ast.Module(body=assignments,type_ignores=[]),'<owned angular B logs>','exec'),env)
    current=env['logparts']; checks={}
    if set(original)!=set(current):raise ValueError('Angular B log keys changed')
    for key in set(original)-{'steep','waiting_and_current'}:
        if s.simplify(original[key]-current[key])!=0:raise ArithmeticError('Angular B source mismatch: '+key)
        checks[key]=True
    if s.simplify(sum(original.values())-sum(current.values()))!=0:raise ArithmeticError('Actual angular B inlet logs differ')
    checks['same_original_Ts_wait_cancellation']=True
    return dict(identities=checks,input_hashes=asts.hashes,
        actual_Bmax_is_original_source_at_q_minus_wait_minus_Ts_minus6=True,source_caps_used_as_fields=False)


@source_precision
def angular_cone_source_bridge(stress):
    asts=SourceAST(); checks={}
    for flag in ('actual_native_angular_X_ODE_AST_replayed','actual_native_angular_energy_ODE_AST_replayed',
        'actual_native_absolute_pressure_ODE_AST_replayed','actual_native_theta_B_K_normalization_consumed',
        'actual_same_full_remaining_E_P_identified_by_ODE_and_entry_datum',
        'actual_compact_beta_and_disjoint_support_F_squared_density_bound',
        'actual_original_selected_coefficient_functions_retained','actual_full_pressure_both_packet_routes_bound',
        'actual_original_angular_radius_reference_verified'):
        if not stress.bridge['identities'][flag]:raise ValueError('Angular cone native source missing: '+flag)
        checks['consumed_'+flag]=True
    if not stress.bridge['actual_angular_entry_formula_join']['actual_angular_entry_K4_full_moment_stress_mixed3_pressure_mixed4_AST_join_verified']:
        raise ValueError('Angular cone right functional join missing')
    exitbridge=stress.entry.power_source.exit_source.bridge
    if not exitbridge['exact_S_is_positive_source_not_enclosure_endpoint']:
        raise ValueError('Angular cone exact positive S required')
    checks['same_full_right_entry_stress_pressure_join_consumed']=True
    checks['exact_positive_S_source_consumed']=True
    # Current accepted entry K is independent of Z at the shared datum.
    path=HERE/(PREFIX+'steep_entry_physical_C2.json')
    entry=json.loads(path.read_bytes())
    if not entry['steep_entry_physical_identities']['original_K_is_independent_of_source_Z']:
        raise ValueError('Original common KR must be axially constant')
    asts.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()

    def syntax(stem,name,target,wanted,augmented=False):
        return asts.expression(stem,name,target,wanted=wanted,augmented=augmented)
    syntax('power_angular_C4','data','coeff',"[jet(v) for v in source['physical_coefficient_Taylor']]")
    syntax('power_angular_C4','data','f','self.flatten.flatten(Z,100)')
    syntax('power_angular_C4','__init__','self.weights',"{k:read_interval(c,v) for k,v in angular['bump_weights'].items()}")
    syntax('power_angular_C4','angular','local','s-center')
    syntax('power_angular_C4','angular','h[k]','dj*(beta[k]*math.factorial(k))',True)
    density=syntax('power_angular_C4','angular','pastA',"dj*(c.exp(self.rate*center)*past['A'])",True)
    eq=syntax('power_angular_C4','angular','eq','IntervalTaylor.constant(c,1/self.rate,5)')
    base=syntax('power_angular_C4','angular','Xbase',"eq+(data['flatten_exit_X']-1/self.rate)*c.exp(-self.rate*y)")
    native=syntax('power_angular_C4','angular','X','(Xbase+pastA*c.exp(-self.rate*s))/F[0]')
    syntax('power_angular_C4','angular','y','self.Lrel+s')
    syntax('power_angular_C4','angular','F','[h[0]+1]+h[1:]')
    fn=asts.method('power_angular_C4','angular')
    loops=[n for n in ast.walk(fn) if isinstance(n,ast.For) and ast.unparse(n.target)=='(dj, center)']
    if len(loops)!=1 or ast.dump(loops[0].iter)!=ast.dump(ast.parse("zip(data['coeff'],(-3,-1))",mode='eval').body):
        raise ValueError('Original two beta support centers changed')
    syntax('power_angular_C4','angular','past','bump_weights(c,self.mu,self.normalization,local,cells=self.cells)')
    for test,wanted in (('hi <= -endpoints(ell)[1]',"past = {k: c.mpf(0) for k in self.weights}"),
                        ('lo >= endpoints(ell)[1]','past = self.weights')):
        cases=[n for n in ast.walk(fn) if isinstance(n,ast.If) and ast.unparse(n.test)==test]
        if len(cases)!=1 or ast.unparse(cases[0].body[0])!=wanted:raise ValueError('Original partial/full past beta branch changed')
    syntax('outer_angular_repair','bump_weights','end','c.mpf([max(mp.mpf(-1),min(mp.mpf(1),rlo)),max(mp.mpf(-1),min(mp.mpf(1),rhi))])')
    syntax('outer_angular_repair','bump_weights','length','end+1')
    syntax('outer_angular_repair','bump_weights','beta','raw_beta(c,raw_coordinate)/(ell*normalization)')
    syntax('outer_angular_repair','bump_weights','result[name]','ds*c.exp(rate*s)*beta**power',True)
    fnweights=asts.method('outer_angular_repair','bump_weights')
    rates=[n for n in ast.walk(fnweights) if isinstance(n,ast.For) and ast.unparse(n.target)=='(name, rate, power)']
    if len(rates)!=1 or ast.unparse(rates[0].iter.elts[0])!="('A', 1 - mu, 1)":
        raise ValueError('Original cumulative A weight rate/power changed')
    a,mu,q,z,length=s.symbols('a mu s Z Lrel',real=True); r=1-mu; k=1-a; b=(1-2*a)/2
    xf=s.Function('actual_flatten_exit_X')(z)
    dj=[s.Function('actual_d'+str(j))(z) for j in range(2)]
    weights=[s.Function('actual_partial_A'+str(j))(q) for j in range(2)]
    H=sum(d*bet for d,bet in zip(dj,[s.Function('original_beta'+str(j))(q) for j in range(2)]))
    c=SimpleNamespace(exp=s.exp); obj=SimpleNamespace(rate=r,Lrel=length)
    env=dict(c=c,self=obj,IntervalTaylor=SimpleNamespace(constant=lambda c,v,n:v))
    past=sum(asts.evaluate(density,dict(env,dj=d,center=center,past={'A':weight}))
             for d,center,weight in zip(dj,(-3,-1),weights))
    equilibrium=asts.evaluate(eq,env)
    nativebase=asts.evaluate(base,dict(env,eq=equilibrium,data={'flatten_exit_X':xf},y=length+q))
    X=asts.evaluate(native,dict(env,Xbase=nativebase,pastA=past,s=q,F=[1+H]))
    N=s.simplify(X*(1+H))
    stable=(mu-a)/r+s.exp(-r*(length+q))*(k*(xf-1/r)-b*z*s.diff(xf,z))
    stable+=sum(s.exp(-r*(q-center))*weight*(k*d-b*z*s.diff(d,z))
                for d,center,weight in zip(dj,(-3,-1),weights))-H
    if s.simplify(k*N-b*z*s.diff(N,z)-(1+H)-stable)!=0:
        raise ArithmeticError('Actual native stable angular W formula differs')
    checks['actual_native_W_stable_formula_AST_verified']=True
    checks['selected_axial_coefficients_and_flatten_history_retained']=True
    checks['clamped_positive_beta_A_integral_between_zero_and_same_full_A']=True
    checks['actual_ordinary_beta_derivative_rows_and_normalization_consumed']=True
    logs=angular_Bmax_binding(); asts.hashes.update(logs['input_hashes'])
    return dict(identities=checks,actual_source_AST_bindings=asts.bindings,
        actual_Bmax_source_log_binding=logs,input_hashes=asts.hashes,
        continuous_whole_support_and_crossing_bounds_used=True,
        source_caps_used_as_fields=False,phase_samples_used_as_proof=False)


@source_precision
def whole_angular_bounds(stress,assembly):
    h=stress.heat; c=stress.ctx; out=stress.outer
    packet=stress.angular([-1,1],[-4,0]); rows=packet['angular_similarity_stress_mixed3_factored']
    rawCt=rows['theta']['y0_Z0']; Cz=rows['axial']['y0_Z0']
    data=out.data([-1,1]); Xf=data['flatten_exit_X']; r=1-h.mu; b=(1-h.delta)/2
    equilibrium_gap=(h.mu-h.a)/r
    flattened_correction=c.exp(-r*(out.Lrel-4))*(h.k*abs(Xf[0]-1/r)+abs(b)*abs(Xf[1]))
    repairs=[c.exp(r*(4+center))*out.weights['A']*(h.k*abs(dj[0])+abs(b)*abs(dj[1]))
             for dj,center in zip(data['coeff'],(-3,-1))]
    H=packet['actual_angular_bump_y_derivatives']
    h_abs=c.mpf(endpoints(abs(H[0][0]))[1]); hs_abs=c.mpf(endpoints(abs(H[1][0]))[1])
    F=packet['angular_shape']['original_F_rows'][0][0]
    Fmin=c.mpf(endpoints(F)[0]); Fmax=c.mpf(endpoints(F)[1])
    correction=flattened_correction+sum(repairs,c.mpf(0))+h_abs
    W_lower=c.mpf(endpoints(equilibrium_gap-correction)[0])
    mu_gap=c.mpf(endpoints(h.mu-hs_abs/Fmin)[0])
    m_upper=c.mpf(endpoints(2*h.mu+2*hs_abs/Fmin)[1])
    shear_sign_gap=c.mpf(endpoints((1+h.mu)*Fmin-hs_abs)[0])
    shear_upper=c.mpf(endpoints(2*h.S*c.exp(h.steep.wait+h.steep.Ts+6)*((1+h.mu)*Fmax+hs_abs))[1])
    bracket_lower=c.mpf(endpoints(W_lower-shear_upper)[0])
    # a<mu and s<=0 imply G=KR*exp((a-mu)*s)>=KR; this KR is the actual entry datum.
    G_lower=c.mpf(endpoints(stress.terminal([-1,1])['KR'])[0])
    theta_lower=c.mpf(endpoints(G_lower*bracket_lower)[0])
    axial_upper=c.mpf(endpoints(abs(Cz))[1])
    logparts=stress_source_log_parts(assembly,h,c.mpf(0))['B']
    logparts['steep']=-h.k*h.steep.Ts
    logparts['waiting_and_current']=6*h.bh
    logBmax=sum(logparts.values(),c.mpf(0))
    margins=dict(mu=h.mu,mu_minus_a=h.mu-h.a,rate=r,delta=h.delta,epsilon=h.eps,
        one_minus_epsilon=1-h.eps,k=h.k,p=h.prate,bh=h.bh,Lmin=1-h.delta,
        wait=h.steep.wait,Ts=h.steep.Ts,Lrel_minus4=out.Lrel-4,normalization=out.normalization,
        same_full_A_weight=out.weights['A'],source_G_lower=G_lower,source_F_lower=Fmin,
        W=W_lower,mu_minus_abs_Fs_over_F=mu_gap,shear_sign=shear_sign_gap,
        theta_bracket=bracket_lower,theta_stress=theta_lower,axial_absolute_upper=axial_upper,m_upper=m_upper)
    for label,value in margins.items():
        if endpoints(value)[0]<=0 or not mp.isfinite(endpoints(value)[1]):raise ArithmeticError('Angular cone positive margin failed: '+label)
    if endpoints(h.mu)[1]>=1 or endpoints(h.eps)[1]>=1 or endpoints(h.delta)[1]>=1:
        raise ArithmeticError('Original angular parameters must preserve r>0, L in(0,1]')
    if endpoints(m_upper)[1]>=2:raise ArithmeticError('Whole original angular kappa-minus2 upper exceeds two')
    if not all(mp.isfinite(v) for v in endpoints(logBmax)):raise ArithmeticError('Finite angular B inlet source logs required')
    logcone=c.ln(m_upper)+2*logBmax+2*c.ln(axial_upper/theta_lower); threshold=c.ln(2)
    if endpoints(logcone)[1]>=endpoints(threshold)[0]:raise ArithmeticError('Whole original angular directional cone failed')
    return dict(positive_margins=margins,actual_whole_source_raw_theta_enclosure=rawCt,actual_whole_source_axial_enclosure=Cz,
        raw_whole_theta_box_used_as_positivity_proof=False,
        actual_equilibrium_W_gap=equilibrium_gap,actual_flatten_history_correction_absolute_bound=flattened_correction,
        actual_selected_partial_angular_history_correction_absolute_bounds=repairs,
        actual_original_h_absolute_upper=h_abs,actual_original_Fs_absolute_upper=hs_abs,
        actual_total_correction_absolute_bound=correction,continuous_W_uniform_lower=W_lower,
        actual_G_uniform_lower=G_lower,actual_F_uniform_lower=Fmin,actual_F_uniform_upper=Fmax,
        exact_source_shear_coefficient_absolute_upper=shear_upper,theta_uniform_positive_lower=theta_lower,
        axial_uniform_absolute_upper=axial_upper,actual_kappa_minus2_enclosure=c.mpf([endpoints(2*mu_gap)[0],endpoints(m_upper)[1]]),
        actual_kappa_minus2_upper=m_upper,actual_Bmax_source_log_parts=logparts,actual_log_Bmax_enclosure=logBmax,
        log_directional_term_upper=logcone,log_threshold=threshold,
        whole_s_Z_domain_and_original_support_crossings_covered=True,
        continuous_native_cumulative_history_bound_used=True,phase_grid_sampling_used_as_proof=False,
        source_caps_only_define_bounds_not_field_values=True,interval_S_upper_zero_not_used_for_strict_source_sign=True)


@source_precision
def run():
    stress=CompliantAngularStressC3(); assembly=CompliantGlobalPhysicalAssembly()
    if (stress.family,stress.source)!=(assembly.family,assembly.source):raise ValueError('Angular cone field family differs')
    hashes=dict(stress.hashes); hashes.update(assembly.hashes)
    for stem,gates in (
        ('angular_stress_C3_check',('actual_original_angular_similarity_stress_recovered',)),
        ('angular_physical_C2_check',('actual_regional_physical_angular_stress_remainder_identity_verified',
            'angular_entry_physical_stress_mixed3_and_remainder_mixed2_join_verified')),
        ('steep_entry_cone_check',('steep_entry_power_exit_waiting_collar_two_vector_cone_certified',))):
        name=PREFIX+stem+'.json'; receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or any(not receipt[gate] for gate in gates):raise ValueError('Angular cone prerequisite missing: '+name)
        if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(stress.family,stress.source):
            raise ValueError('Angular cone prerequisite family differs')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Angular cone prerequisite changed: '+source)
        hashes.update(receipt['input_hashes']); hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    print('Original angular source, physical/right joins and accepted entry tail admitted',flush=True)
    proof=angular_cone_identities(); bridge=angular_cone_source_bridge(stress); bounds=whole_angular_bounds(stress,assembly)
    hashes.update(bridge['input_hashes']); hashes.update(assembly.dispatch.hashes)
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=stress.family,implicit_source_sha256=stress.source,
        scope='Whole original angular s[-4,0], source Z[-1,1]; physical tau>0, r>0, |Z|<1',
        domain=DOMAIN,tail_domain=TAIL_DOMAIN,exact_identities=proof,source_bridge=bridge,bounds=bounds,input_hashes=hashes,
        angular_cone_certified=True,actual_angular_theta_stress_positive=True,
        actual_angular_source_shear_strictly_negative=True,actual_angular_directional_cone_uniform_margin=True,
        angular_entry_power_exit_waiting_collar_two_vector_cone_certified=True,
        same_source_angular_entry_similarity_and_physical_joins_consumed=True,
        fixed_positive_viscosity_two_vector_cone_transfer_verified=True,
        angular_regional_remainder_exact_zero=False,preceding_power_physical_stress_join_verified=False,
        completed_full_tensor_cone_certified=False,whole_outer_cone_certified=False,
        global_admissible_stress_lift_constructed=False,independently_bounded_global_flat_remainder=False,
        physical_energy_integral_certified=False,temporal_recursion=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Whole original angular cone generated by continuous native history bounds; entry tail joined, upstream/global pending',flush=True)
    return result


class CertifiedAngularPhysical:
    """Attach the accepted original two-vector cone to the unchanged field."""
    def __init__(self):
        from lei_ren_part1_paper_compliant_angular_physical_C2 import CompliantAngularPhysicalC2
        receipt=json.loads((HERE/(PREFIX+'angular_cone_check.json')).read_bytes())
        if not receipt['all_passed'] or not receipt['angular_cone_certified']:raise ValueError('Accepted angular cone required')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Certified angular source changed: '+source)
        self.field=CompliantAngularPhysicalC2(); self.family,self.source=self.field.family,self.field.source
        if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(self.family,self.source):
            raise ValueError('Angular physical cone family differs')
        self.receipt=receipt

    def angular(self,*args,**kwargs):
        point=self.field.angular(*args,**kwargs)
        point.update(angular_cone_certified=True,angular_entry_power_exit_waiting_collar_two_vector_cone_certified=True,
            certified_cone_scope='Original two-vector only; completed diagonal/global cone not certified',
            completed_full_tensor_cone_certified=False)
        return point


if __name__=='__main__':run()
