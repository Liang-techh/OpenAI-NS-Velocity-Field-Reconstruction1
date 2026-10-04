"""Whole original preceding-power cone from correlated flatten-exit history.

The actual long Lrel, canonical moments/pressure, velocities and physical
units remain. This regional two-vector cone does not establish global
completed-tensor admissibility, flatness, energy or coefficient recursion.
"""
import ast
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_outer_power_stress_C3 import CompliantOuterPowerStressC3,SourceAST,DOMAIN
from lei_ren_part1_paper_compliant_waiting_stress_C3 import source_precision
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_collar_physical_C2 import stress_source_log_parts
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
TAIL_DOMAIN='Rtail*exp(-wait-Ts-2-Lrel)<=R<Rtail*exp(3); Gamma stress exactly zero beyond'


def outer_power_cone_identities():
    a,mu,u,z,omega=s.symbols('a mu u Z omega',real=True)
    r=1-mu;k=1-a;b=(1-2*a)/2;f=(1+z*z)/2;j=2*z*z/(1+z*z)
    Xf=s.Function('actual_flatten_exit_X')(z)
    X=1/r+(Xf-1/r)*s.exp(-r*u)
    H=k*(Xf-1/r)-b*z*s.diff(Xf,z)
    W=(mu-a)/r+H*s.exp(-r*u)
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Outer power cone identity failed: '+name)
        proofs[name]=True
    zero('actual_native_W_stable_formula',k*X-b*z*s.diff(X,z)-1-W)
    zero('actual_continuous_W_ODE',s.diff(W,u)+r*W-(mu-a))
    xv=s.Symbol('actual_positive_Xv');t=s.Symbol('original_flatten_length',positive=True)
    I=s.Function('actual_weighted_flatten_integral')(z)
    xf=xv*s.exp(-r*t)/f+I
    past=s.exp(-r*t)*((k+b*j)*xv/f-k/r)
    zero('actual_flatten_past_axial_correlation',k*(xf-1/r)-b*z*s.diff(xf,z)
        -(past+k*I-b*z*s.diff(I,z)-k/r*(1-s.exp(-r*t))))
    density=s.exp(-omega*s.log(f))
    zero('actual_flatten_integrand_axial_correlation',
        k*(density-1)-b*z*s.diff(density,z)
        -(k*(density-1)+b*j*omega*density))
    v=s.Symbol('Z_squared',nonnegative=True)
    zero('spatial_log_lower_bound_remainder',
        k*(1-v)/2+b*2*v/(1+v)-b
        -(a*(1-v)/2+b*v*(1-v)/(1+v)))
    ff=s.Symbol('f',positive=True)
    zero('log_tangent_bound_derivative',s.diff(-s.log(ff)-(1-ff),ff)+(1-ff)/ff)
    sigx=s.Symbol('sigma_x',positive=True)
    zero('original_sigma_left_half_odds_factor',
        1/(1-sigx)**2-1/sigx**2-(2*sigx-1)/(sigx**2*(1-sigx)**2))
    zero('actual_flatten_positive_lower_regrouping',
        b*(s.exp(-50*r)-s.exp(-100*r))/(2*r)-k*s.exp(-100*r)/r
        -s.exp(-100*r)*(b*(s.exp(50*r)-1)/2-k)/r)
    q,KR,S=s.symbols('q Kright source_inverse_Rtail',real=True)
    K=KR*s.exp((a-mu)*q)
    zero('actual_source_shear',s.diff(K,q)-(1+a)*K+(1+mu)*K)
    zero('actual_constant_kappa_minus2',2*a-2*(a-mu)-2*mu)
    wait,Ts,length=s.symbols('wait Ts Lrel',real=True)
    bp=s.Rational(1,2)+mu;bh=s.Rational(1,2)+a
    zero('actual_Bmax_shared_Ts_Lrel_wait_regrouping',
        -bp*(100+length)-s.Rational(3,2)*Ts+bh*(Ts+2+length)
        -(-100*bp-k*Ts+(a-mu)*length+2*bh))
    qt,ct,cz,st,Bcurrent,m=s.symbols('Qtheta Ctheta Cz Stheta B kappa_minus2',real=True)
    zero('normalized_original_two_vector_cone',
        2*(qt*ct*qt*st)**2-m*(qt**2*Bcurrent*cz*st)**2
        -qt**4*st**2*(2*ct**2-m*Bcurrent**2*cz**2))
    return dict(identities=proofs,whole_original_domain=DOMAIN,
        actual_native_W='(mu-a)/(1-mu)+Hf(Z)*exp(-(1-mu)*u), u=(Lrel-4)*phase',
        actual_correlated_Hf_integrand='k*(f^(-omega)-1)+b*j*omega*f^(-omega)',
        positive_Hf_lemma=dict(f_range=[.5,1],omega_range=[0,1],Z_squared_range=[0,1],
            exponential_tangent='exp(x)-1>=x for x>=0',
            log_tangent='-log(f)>=1-f for 0<f<=1',
            spatial_factor='k*(1-Z^2)/2+b*2*Z^2/(1+Z^2)>=b',
            original_left_half_sigma='sigma(v/100)<=1/2 for 0<=v<=50',
            strict_parameter_condition='b*expm1(50*(1-mu))-2*k>0',
            Hf_lower='exp(-100*r)*(b*expm1(50*r)/2-k)/r'),
        sufficient_directional_condition='2*mu*Bmax^2*(Cz_abs_upper/Ctheta_lower)^2<2',
        exact_S_positive_inverse_Rtail_not_cap_endpoint=True,
        same_positive_physical_nu_lambda_factors_preserve_two_vector_cone=True,
        completed_diagonal_does_not_enter_two_vector_cone=True)


def outer_power_Bmax_binding():
    asts=SourceAST();fn=asts.method('collar_physical_C2','stress_source_log_parts')
    env={};exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual power B logs>','exec'),env)
    a,mu,wait,Ts,length,lp,lu,lrp,lone=s.symbols('a mu wait Ts Lrel logP logU logRp logone',real=True)
    c=SimpleNamespace(mpf=s.sympify,ln=s.log)
    flatten=SimpleNamespace(logEv2_parts=dict(inlet_log=2*lu))
    assembly=SimpleNamespace(ctx=c,logP=lp,logRp=lrp,dispatch=SimpleNamespace(provider=lambda chart:flatten))
    h=SimpleNamespace(ctx=c,a=a,mu=mu,k=1-a,bh=s.Rational(1,2)+a,
        steep=SimpleNamespace(wait=wait,Ts=Ts,logone=lone,outer=SimpleNamespace(Lrel=length)))
    original=env['stress_source_log_parts'](assembly,h,-wait-Ts-2-length)['B']
    fn=asts.method('outer_power_cone','whole_outer_power_bounds')
    names=('logparts',"logparts['steep']","logparts['flatten_power']","logparts['waiting_and_current']")
    assignments=[n for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t) in names for t in n.targets)]
    if len(assignments)!=4:raise ValueError('Owned power Bmax regrouping changed')
    env.update(assembly=assembly,h=h,c=c,out=SimpleNamespace(Lrel=length))
    exec(compile(ast.Module(body=assignments,type_ignores=[]),'<owned power B logs>','exec'),env)
    current=env['logparts'];checks={}
    if set(original)!=set(current):raise ValueError('Power B source keys changed')
    for key in set(original)-{'steep','flatten_power','waiting_and_current'}:
        if s.simplify(original[key]-current[key])!=0:raise ArithmeticError('Power B source mismatch: '+key)
        checks[key]=True
    if s.simplify(sum(original.values())-sum(current.values()))!=0:raise ArithmeticError('Original power inlet B logs differ')
    checks['same_original_Ts_Lrel_wait_regrouping']=True
    return dict(identities=checks,input_hashes=asts.hashes,
        actual_Bmax_is_original_source_at_q_minus_wait_minus_Ts_minus2_minusLrel=True,source_caps_used_as_fields=False)


@source_precision
def outer_power_cone_source_bridge(stress):
    asts=SourceAST();checks={}
    for flag in ('actual_native_power_X_ODE_AST_replayed','actual_native_power_energy_ODE_AST_replayed',
        'actual_native_power_absolute_pressure_ODE_AST_replayed','actual_native_power_theta_rate_AST_replayed',
        'actual_native_power_complete_history_identified_by_ODE_and_angular_datum',
        'actual_native_power_theta_equals_same_B_times_K_by_rate_and_right_datum',
        'actual_original_phase_to_ordinary_q_radius_verified','actual_selected_full_angular_future_retained'):
        if not stress.bridge['identities'][flag]:raise ValueError('Power cone actual source missing: '+flag)
        checks['consumed_'+flag]=True
    if not stress.bridge['actual_outer_power_angular_formula_join']['actual_outer_power_angular_K4_full_moment_stress_mixed3_pressure_mixed4_AST_join_verified']:
        raise ValueError('Power cone right functional join missing')
    if not stress.angular.entry.power_source.exit_source.bridge['exact_S_is_positive_source_not_enclosure_endpoint']:
        raise ValueError('Power cone exact positive S source required')
    if not stress.angular.entry.power_source.exit_source.bridge['original_sigma_between0and1_reflection_and_flat_endpoints_bound']:
        raise ValueError('Original flatten sigma source properties required')
    def syntax(stem,name,target,wanted,augmented=False):
        return asts.expression(stem,name,target,wanted=wanted,augmented=augmented)
    syntax('power_angular_C4','data','f','self.flatten.flatten(Z,100)')
    native=syntax('power_angular_C4','power','X',"one/self.rate+(data['flatten_exit_X']-1/self.rate)*c.exp(-self.rate*y)")
    syntax('power_angular_C4','power','y','(self.Lrel-4)*phase')
    syntax('flatten_mixed_C4','__init__','self.Xv',"read_interval(c,pulse['whole_Z_terminal']['Mtheta_over_sqrt2_R_3half_Utheta']['coefficients'][0])")
    qnode=syntax('flatten_mixed_C4','flatten','q','IntervalTaylor(c,[1+Z**2,2*Z,1,0,0,0])')
    rhonode=syntax('flatten_mixed_C4','flatten','rho','log_taylor(q)-c.ln(2)')
    syntax('flatten_mixed_C4','flatten','sig','sigma_jets(c,t/100)')
    Fnode=syntax('flatten_mixed_C4','flatten','F','(rho*sj[0]).exp()')
    Fcnode=syntax('flatten_mixed_C4','flatten','Fc','(rho*sigma_jets(c,v/100)[0]).exp()')
    syntax('flatten_mixed_C4','flatten','Xint','q*0')
    syntax('flatten_mixed_C4','flatten','Xint',
        'Fc*(c.exp(-self.rate*(t*(self.cells-i-1)/self.cells))*decay_integral(c,self.rate,length))',True)
    xfnode=syntax('flatten_mixed_C4','flatten','X','(Xint+self.Xv*c.exp(-self.rate*t))/F')
    syntax('flat_pulse_derivatives','_sigma_left','odds','1/(1-x)**2-1/x**2')
    syntax('flat_pulse_derivatives','_sigma_left','e','positive_exp(c,odds)')
    syntax('flat_pulse_derivatives','_sigma_left','value','e/(1+e)')
    asts.method('flat_pulse_derivatives','sigma_jets')
    a,mu,z,u=s.symbols('a mu Z u',real=True);r=1-mu;k=1-a;b=(1-2*a)/2
    xf=s.Function('same_actual_flatten_X')(z)
    env=dict(c=SimpleNamespace(exp=s.exp),self=SimpleNamespace(rate=r),one=s.Integer(1),data={'flatten_exit_X':xf},y=u)
    X=asts.evaluate(native,env)
    W=(mu-a)/r+(k*(xf-1/r)-b*z*s.diff(xf,z))*s.exp(-r*u)
    if s.simplify(k*X-b*z*s.diff(X,z)-1-W)!=0:raise ArithmeticError('Actual native power stable W source differs')
    f=(1+z*z)/2;xv=s.Symbol('same_positive_Xv');I=s.Function('actual_normalized_flatten_weighted_integral')(z)
    scalar=SimpleNamespace(exp=s.exp,ln=s.log)
    qscalar=asts.evaluate(qnode,dict(Z=z,c=scalar,IntervalTaylor=lambda ctx,values:values[0]))
    rho=asts.evaluate(rhonode,dict(q=qscalar,c=scalar,log_taylor=s.log))
    class ScalarExponential(ast.NodeTransformer):
        def visit_Call(self,node):
            node=self.generic_visit(node)
            if isinstance(node.func,ast.Attribute) and node.func.attr=='exp' and not node.args and not node.keywords:
                return ast.copy_location(ast.Call(func=ast.Name(id='source_exp',ctx=ast.Load()),args=[node.func.value],keywords=[]),node)
            return node
    def replay_exponential(node,env):
        node=ast.fix_missing_locations(ScalarExponential().visit(node))
        return asts.evaluate(node,dict(env,source_exp=s.exp))
    Fendpoint=replay_exponential(Fnode,dict(rho=rho,sj=[s.Integer(1)]))
    sig,v=s.symbols('actual_sigma source_flatten_v',real=True)
    Fc=replay_exponential(Fcnode,dict(rho=rho,c=scalar,v=v,sigma_jets=lambda ctx,value:[sig]))
    if s.simplify(Fendpoint-f)!=0 or s.simplify(Fc-s.exp(sig*rho))!=0:
        raise ArithmeticError('Actual original flatten F/Fc exponential sources differ')
    if s.simplify(s.expand_log(sig*rho-s.log(Fendpoint)-(sig-1)*s.log(f),force=True))!=0:
        raise ArithmeticError('Actual normalized flatten density differs from f^(-omega)')
    # Xint is the original unnormalized integral. Its normalized density is
    # Fc/Fendpoint=f^(sigma-1), so source Xint=f*I_norm at this endpoint.
    sourcexf=asts.evaluate(xfnode,dict(Xint=f*I,self=SimpleNamespace(Xv=xv,rate=r),c=scalar,t=s.Integer(100),F=Fendpoint))
    if s.simplify(sourcexf-(I+xv*s.exp(-100*r)/f))!=0:raise ArithmeticError('Actual normalized flatten endpoint history differs')
    c=stress.ctx
    half=sigma_jets(c,c.mpf([0,.5]))[0];whole=sigma_jets(c,c.mpf([0,1]))[0];endpoint=sigma_jets(c,c.mpf(1))[0]
    if endpoints(half)[0]<0 or endpoints(half)[1]>.5 or endpoints(whole)[0]<0 or endpoints(whole)[1]>1 or endpoints(endpoint)!=(1,1):
        raise ValueError('Actual original sigma bounds/endpoints changed')
    checks.update(actual_native_W_stable_formula_AST_verified=True,
        actual_flatten_original_Xv_q_F_and_positive_weighted_integral_source_bound=True,
        actual_flatten_endpoint_value_and_axial_integrand_correlation_verified=True,
        actual_F_endpoint_equals_original_q_over_two_AST_verified=True,
        actual_Fc_over_endpoint_F_equals_f_power_minus_omega_AST_verified=True,
        original_source_Xint_equals_f_times_same_normalized_integral_bound=True,
        actual_original_sigma_left_half_and_whole_bounds_verified=True,
        positive_Hf_proof_uses_correlated_source_integral_not_raw_Xf_boxes=True,
        same_full_right_angular_stress_pressure_join_consumed=True,exact_positive_S_source_consumed=True)
    logs=outer_power_Bmax_binding();asts.hashes.update(logs['input_hashes'])
    return dict(identities=checks,actual_source_AST_bindings=asts.bindings,actual_sigma_half_bound=pack(half),
        actual_sigma_whole_bound=pack(whole),actual_sigma_endpoint=pack(endpoint),
        actual_Bmax_source_log_binding=logs,input_hashes=asts.hashes,
        continuous_whole_original_domain_used=True,source_caps_used_as_fields=False,phase_samples_used_as_proof=False)


@source_precision
def whole_outer_power_bounds(stress,assembly):
    h=stress.heat;c=stress.ctx;out=stress.outer;r=1-h.mu;b=(1-h.delta)/2
    packet=stress.power([-1,1],[0,1]);rows=packet['outer_power_similarity_stress_mixed3_factored']
    rawCt=rows['theta']['y0_Z0'];Cz=rows['axial']['y0_Z0']
    equilibrium_gap=(h.mu-h.a)/r
    H_parameter_gap=b*c.expm1(50*r)-2*h.k
    H_lower=c.mpf(endpoints(c.exp(-100*r)*H_parameter_gap/(2*r))[0])
    W_lower=c.mpf(endpoints(equilibrium_gap)[0])
    shear_upper=c.mpf(endpoints(2*h.S*c.exp(h.steep.wait+h.steep.Ts+2+out.Lrel)*(1+h.mu))[1])
    bracket_lower=c.mpf(endpoints(W_lower-shear_upper)[0])
    # d=a-mu<0 and y<=0: original K is at least actual angular s=-4 Kright.
    K_lower=c.mpf(endpoints(stress.terminal([-1,1])['Kright'])[0])
    theta_lower=c.mpf(endpoints(K_lower*bracket_lower)[0])
    axial_upper=c.mpf(endpoints(abs(Cz))[1]);m_upper=c.mpf(endpoints(2*h.mu)[1])
    logparts=stress_source_log_parts(assembly,h,c.mpf(0))['B']
    logparts['steep']=-h.k*h.steep.Ts
    logparts['flatten_power']=-100*(c.mpf('.5')+h.mu)+(h.a-h.mu)*out.Lrel
    logparts['waiting_and_current']=2*h.bh
    logBmax=sum(logparts.values(),c.mpf(0))
    margins=dict(mu=h.mu,mu_minus_a=h.mu-h.a,rate=r,delta=h.delta,epsilon=h.eps,
        one_minus_epsilon=1-h.eps,k=h.k,p=h.prate,bh=h.bh,b=b,Lmin=1-h.delta,
        wait=h.steep.wait,Ts=h.steep.Ts,Lrel_minus4=out.Lrel-4,actual_Xv=out.flatten.Xv,
        correlated_Hf_parameter_gap=H_parameter_gap,correlated_Hf_lower=H_lower,
        W=W_lower,source_K_lower=K_lower,shear_sign=1+h.mu,
        theta_bracket=bracket_lower,theta_stress=theta_lower,axial_absolute_upper=axial_upper,m_upper=m_upper)
    for label,value in margins.items():
        if endpoints(value)[0]<=0 or not mp.isfinite(endpoints(value)[1]):raise ArithmeticError('Power cone positive margin failed: '+label)
    if endpoints(h.mu)[1]>=1 or endpoints(h.eps)[1]>=1 or endpoints(h.delta)[1]>=1 or endpoints(m_upper)[1]>=2:
        raise ArithmeticError('Original power parameter/shear range failed')
    if not all(mp.isfinite(v) for v in endpoints(logBmax)):raise ArithmeticError('Finite original power B inlet logs required')
    logcone=c.ln(m_upper)+2*logBmax+2*c.ln(axial_upper/theta_lower);threshold=c.ln(2)
    if endpoints(logcone)[1]>=endpoints(threshold)[0]:raise ArithmeticError('Whole original power directional cone failed')
    return dict(positive_margins=margins,actual_whole_source_raw_theta_enclosure=rawCt,actual_whole_source_axial_enclosure=Cz,
        raw_whole_theta_box_used_as_positivity_proof=False,actual_equilibrium_W_gap=equilibrium_gap,
        actual_correlated_flatten_Hf_uniform_lower=H_lower,actual_Hf_strict_parameter_gap=H_parameter_gap,
        negative_past_equilibrium_bound=h.k*c.exp(-100*r)/r,
        continuous_W_uniform_lower=W_lower,actual_K_uniform_lower=K_lower,
        exact_source_shear_coefficient_absolute_upper=shear_upper,theta_uniform_positive_lower=theta_lower,
        axial_uniform_absolute_upper=axial_upper,actual_kappa_minus2_enclosure=2*h.mu,
        actual_kappa_minus2_upper=m_upper,actual_Bmax_source_log_parts=logparts,actual_log_Bmax_enclosure=logBmax,
        log_directional_term_upper=logcone,log_threshold=threshold,
        whole_original_phase_Z_domain_covered=True,correlated_flatten_history_uniform_lower_used=True,
        continuous_native_cumulative_history_bound_used=True,phase_grid_sampling_used_as_proof=False,
        source_caps_only_define_bounds_not_field_values=True,interval_S_lower_zero_not_used_for_strict_source_sign=True)


@source_precision
def run():
    stress=CompliantOuterPowerStressC3();assembly=CompliantGlobalPhysicalAssembly()
    if (stress.family,stress.source)!=(assembly.family,assembly.source):raise ValueError('Power cone field family differs')
    hashes=dict(stress.hashes);hashes.update(assembly.hashes)
    for stem,gates in (
        ('outer_power_stress_C3_check',('actual_original_outer_power_similarity_stress_recovered',)),
        ('outer_power_physical_C2_check',('actual_regional_physical_outer_power_stress_remainder_identity_verified',
            'outer_power_angular_physical_stress_mixed3_and_remainder_mixed2_join_verified')),
        ('angular_cone_check',('angular_entry_power_exit_waiting_collar_two_vector_cone_certified',))):
        name=PREFIX+stem+'.json';receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or any(not receipt[gate] for gate in gates):raise ValueError('Power cone prerequisite missing: '+name)
        if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(stress.family,stress.source):
            raise ValueError('Power cone prerequisite family differs')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Power cone prerequisite changed: '+source)
        hashes.update(receipt['input_hashes']);hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    proof=outer_power_cone_identities();bridge=outer_power_cone_source_bridge(stress);bounds=whole_outer_power_bounds(stress,assembly)
    hashes.update(bridge['input_hashes']);hashes.update(assembly.dispatch.hashes)
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=stress.family,implicit_source_sha256=stress.source,
        scope='Whole original outer_power phase[0,1],Z[-1,1]; physical tau>0,r>0,|Z|<1',
        domain=DOMAIN,tail_domain=TAIL_DOMAIN,exact_identities=proof,source_bridge=bridge,bounds=bounds,input_hashes=hashes,
        outer_power_cone_certified=True,actual_outer_power_theta_stress_positive=True,
        actual_outer_power_source_shear_strictly_negative=True,actual_outer_power_directional_cone_uniform_margin=True,
        outer_power_angular_entry_exit_waiting_collar_two_vector_cone_certified=True,
        same_source_power_angular_similarity_and_physical_joins_consumed=True,
        fixed_positive_viscosity_two_vector_cone_transfer_verified=True,
        outer_power_regional_remainder_exact_zero=False,left_flatten_physical_stress_join_verified=False,
        completed_full_tensor_cone_certified=False,whole_outer_cone_certified=False,
        global_admissible_stress_lift_constructed=False,independently_bounded_global_flat_remainder=False,
        physical_energy_integral_certified=False,temporal_recursion=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Whole original preceding-power cone generated from correlated flatten history; angular tail joined, flatten/global pending',flush=True)
    return result


class CertifiedOuterPowerPhysical:
    """Attach the regional original two-vector cone to the unchanged field."""
    def __init__(self):
        from lei_ren_part1_paper_compliant_outer_power_physical_C2 import CompliantOuterPowerPhysicalC2
        receipt=json.loads((HERE/(PREFIX+'outer_power_cone_check.json')).read_bytes())
        if not receipt['all_passed'] or not receipt['outer_power_cone_certified']:raise ValueError('Accepted power cone required')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Certified power source changed: '+source)
        self.field=CompliantOuterPowerPhysicalC2();self.family,self.source=self.field.family,self.field.source
        if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(self.family,self.source):
            raise ValueError('Power physical cone family differs')
        self.receipt=receipt
    def power(self,*args,**kwargs):
        point=self.field.power(*args,**kwargs)
        point.update(outer_power_cone_certified=True,outer_power_angular_entry_exit_waiting_collar_two_vector_cone_certified=True,
            certified_cone_scope='Original two-vector only; completed diagonal/global cone not certified',
            completed_full_tensor_cone_certified=False)
        return point


if __name__=='__main__':run()
