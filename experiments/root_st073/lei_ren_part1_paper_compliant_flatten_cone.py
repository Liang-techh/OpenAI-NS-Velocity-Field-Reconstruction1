"""Whole original flatten cone from signed native pulse memory and variable F.

Same fields/full moments/absolute pressure, original100 units and physical
units remain. Regional two-vector cone is distinct from full/global gates.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_flatten_stress_C3 import CompliantFlattenStressC3,SourceAST,DOMAIN,source_precision
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_collar_physical_C2 import stress_source_log_parts
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
TAIL_DOMAIN='Rtail*exp(-wait-Ts-102-Lrel)<=R<Rtail*exp(3); Gamma stress exactly zero beyond'


def flatten_cone_identities():
    a,mu,t,z,Xp=s.symbols('a mu t Z same_original_Xp',real=True)
    r=1-mu;k=1-a;b=(1-2*a)/2;g=(mu-a)/r
    f=(1+z*z)/2;rho=s.log(f);j=2*z*z/(1+z*z)
    sig=s.Function('original_sigma')(t);F=s.exp(rho*sig)
    N=s.Function('native_N')(t,z);K=s.exp((a-mu)*(t-100))*F/f
    X=N/F;A=K*X
    W=((k+b*j)*N-b*z*s.diff(N,z))/F-1
    rules={s.diff(N,t):F-r*N,s.diff(N,t,z):s.diff(F,z)-r*s.diff(N,z)}
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Flatten cone identity failed: '+name)
        proofs[name]=True
    zero('actual_same_KX_inertial_numerator',(k*A-b*z*s.diff(A,z)-K)/K-W)
    forcing=mu-a+b*j*(1-sig)-rho*s.diff(sig,t)
    zero('actual_native_W_variable_rate_ODE',
        (s.diff(W,t)+(r+rho*s.diff(sig,t))*W-forcing).subs(rules,simultaneous=True))
    V=W-g;positive=b*j*(1-sig)-rho*s.diff(sig,t)*(1+g)
    zero('actual_shifted_V_equilibrium_ODE',
        (s.diff(V,t)+(r+rho*s.diff(sig,t))*V-positive).subs(rules,simultaneous=True))
    Iv=s.Function('same_positive_shifted_forcing_integral')(t,z)
    V0=s.Function('same_signed_initial_V')(z)
    representation=g+s.exp(-r*t)/F*(V0+Iv)
    zero('actual_continuous_integrating_factor_representation',
        (s.diff(representation,t)+(r+rho*s.diff(sig,t))*(representation-g)-positive)
        .subs(s.diff(Iv,t),s.exp(r*t)*F*positive))
    Xv=1/r+(Xp-1/r)*s.exp(-13*r/mu)
    zero('actual_stable_incoming_kXv_minus_one',k*Xv-1-g-k*(Xp-1/r)*s.exp(-13*r/mu))
    zero('actual_shifted_initial_V_with_axial_factor',(k+b*j)*Xv-1-g
        -(b*j*Xv+k*(Xp-1/r)*s.exp(-13*r/mu)))
    zero('actual_original_flatten_kappa_minus2',2*a-2*(a-mu+rho*s.diff(sig,t))
        -(2*mu-2*rho*s.diff(sig,t)))
    zero('actual_strict_source_shear',s.diff(K,t)-(1+a)*K-K*(-(1+mu)+rho*s.diff(sig,t)))
    bp=s.Rational(1,2)+mu;bh=s.Rational(1,2)+a
    Ts,L=s.symbols('Ts Lrel',real=True)
    zero('actual_inlet_B_log_correlation',-bp*(100+L)-s.Rational(3,2)*Ts+bh*(100+L+Ts+2)
        -(-(1-a)*Ts+(a-mu)*(100+L)+2*bh))
    qt,ct,cz,st,B,m=s.symbols('Qtheta Ctheta Cz Stheta B m',real=True)
    zero('original_two_vector_cone_normalization',
        2*(qt*ct*qt*st)**2-m*(qt**2*B*cz*st)**2-qt**4*st**2*(2*ct**2-m*B**2*cz**2))
    return dict(identities=proofs,whole_original_domain=DOMAIN,
        original_Xv='1/r+(Xp-1/r)*exp(-13*r/mu), Z-independent',
        actual_W='((k+b*j)*N-b*Z*N_Z)/F-1; N=F*X,j=2*Z^2/(1+Z^2)',
        actual_W_ODE='W_t+(r+rho*sigma_t)*W=mu-a+b*j*(1-sigma)-rho*sigma_t',
        actual_shifted_V_ODE='V_t+(r+rho*sigma_t)*V=b*j*(1-sigma)-rho*sigma_t*(1+g), V=W-g',
        continuous_uniform_lower='W>=g-2*k*abs(Xp-1/r)*exp(-13*r/mu), g=(mu-a)/r',
        positive_forcing_requirements=dict(rho_nonpositive=True,sigma_increasing=True,
            sigma_between0and1=True,b_and_j_nonnegative=True,actual_Xv_positive_required=True,F_at_least_one_half=True),
        source_caps_used_as_fields=False,completed_diagonal_is_outside_two_vector_cone=True)


def flatten_Bmax_binding():
    asts=SourceAST();fn=asts.method('collar_physical_C2','stress_source_log_parts');env={}
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual flatten B logs>','exec'),env)
    a,mu,wait,Ts,L,lp,lu,lrp,lone=s.symbols('a mu wait Ts Lrel logP logU logRp logone',real=True)
    c=SimpleNamespace(mpf=s.sympify,ln=s.log)
    assembly=SimpleNamespace(ctx=c,logP=lp,logRp=lrp,dispatch=SimpleNamespace(provider=lambda chart:
        SimpleNamespace(logEv2_parts=dict(inlet_log=2*lu))))
    h=SimpleNamespace(ctx=c,a=a,mu=mu,k=1-a,bh=s.Rational(1,2)+a,
        steep=SimpleNamespace(wait=wait,Ts=Ts,logone=lone,outer=SimpleNamespace(Lrel=L)))
    original=env['stress_source_log_parts'](assembly,h,-wait-Ts-102-L)['B']
    fn=asts.method('flatten_cone','whole_flatten_bounds')
    names=('logparts',"logparts['steep']","logparts['flatten_power']","logparts['waiting_and_current']")
    assignments=[node for node in ast.walk(fn) if isinstance(node,ast.Assign) and any(ast.unparse(v) in names for v in node.targets)]
    if len(assignments)!=4:raise ValueError('Owned flatten B regrouping changed')
    env.update(assembly=assembly,h=h,c=c,out=SimpleNamespace(Lrel=L))
    exec(compile(ast.Module(body=assignments,type_ignores=[]),'<owned flatten B logs>','exec'),env)
    current=env['logparts'];checks={}
    if set(original)!=set(current):raise ValueError('Flatten B source keys changed')
    for key in set(original)-{'steep','flatten_power','waiting_and_current'}:
        if s.simplify(original[key]-current[key])!=0:raise ArithmeticError('Flatten B unchanged source key differs')
        checks[key]=True
    if s.simplify(sum(original.values())-sum(current.values()))!=0:raise ArithmeticError('Flatten original inlet B differs')
    checks['same_original_Ts_Lrel_wait_regrouping']=True
    return dict(identities=checks,input_hashes=asts.hashes,
        actual_Bmax_is_original_source_at_q_minus_wait_minus_Ts_minus102_minusLrel=True,source_caps_used_as_fields=False)


def flatten_cone_source_bridge(stress):
    asts=SourceAST();checks={}
    def syntax(stem,name,target,wanted=None,augmented=False):
        return asts.expression(stem,name,target,wanted=wanted,augmented=augmented)
    syntax('axial_pulse_field','__init__','inlet',"""self.pulse.buffer.power('0',1)""")
    originalXp=syntax('axial_pulse_field','__init__','self.Xp',
        """inlet['Mtheta_over_sqrt2_R_3half_Pstar'][0]/inlet['Utheta_over_Pstar'][0]""")
    nativeXv=syntax('axial_pulse_field','end','X',
        '1/self.rate+(self.Xp-1/self.rate)*self.factor(-13*self.rate/self.mu-self.rate*s)')
    syntax('flatten_mixed_C4','__init__','self.Xv',
        """read_interval(c,pulse['whole_Z_terminal']['Mtheta_over_sqrt2_R_3half_Utheta']['coefficients'][0])""")
    syntax('flatten_mixed_C4','flatten','X','(Xint+self.Xv*c.exp(-self.rate*t))/F')
    syntax('flatten_mixed_C4','flatten','Xint',
        'Fc*(c.exp(-self.rate*(t*(self.cells-i-1)/self.cells))*decay_integral(c,self.rate,length))',True)
    syntax('flatten_mixed_C4','flatten','F','(rho*sj[0]).exp()')
    syntax('flatten_mixed_C4','flatten','Fc','(rho*sigma_jets(c,v/100)[0]).exp()')
    syntax('flatten_mixed_C4','flatten','sj','[sig[k]*math.factorial(k)/100**k for k in range(5)]')
    syntax('power_inlet_C4','__init__','self.Xp',"""self.inlet_H/self.constants['U']""")
    syntax('power_inlet_C4','__init__','self.inlet_H',"""read_interval(c,sample['Mtheta_over_sqrt2_R_3half_Pstar'][0])*q""")
    syntax('flatten_cone','whole_flatten_bounds','Xp','stress.native.inlet.Xp')
    mu,hp,u,z=s.symbols('mu same_Hp same_U Z',real=True);qpoly=1+z*z
    canonicalXp=asts.evaluate(originalXp,dict(inlet=dict(Mtheta_over_sqrt2_R_3half_Pstar=[hp/qpoly],Utheta_over_Pstar=[u/qpoly])))
    exactXv=asts.evaluate(nativeXv,dict(self=SimpleNamespace(rate=1-mu,mu=mu,Xp=canonicalXp,factor=s.exp),s=s.Integer(0)))
    if s.simplify(canonicalXp-hp/u)!=0 or s.simplify(s.diff(exactXv,z))!=0:
        raise ArithmeticError('Original canonical pulse memory is not scalar')
    checks['actual_original_pulse_Xp_ratio_and_exact_terminal_memory_AST_verified']=True
    checks['actual_original_Xv_is_Z_independent_function_not_chosen_box_value']=True
    name=PREFIX+'power_inlet_C4_check.json';inlet=json.loads((HERE/name).read_bytes())
    for flag in ('canonical_Xp_constant','actual_O2_O3_source_chain_has_exact_canonical_whole_Z_shapes',
        'single_sample_only_encloses_proved_Z_independent_Hp_Pin_constants'):
        if not inlet['exact_functional_production_and_join_identities'][flag]:raise ValueError('Actual incoming canonical shape missing')
        checks['consumed_'+flag]=True
    asts.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    for flag in ('actual_normalized_f_is_original_q_over_two_AST_verified',
        'actual_flatten_future_energy_half_of_complete_C5_source_AST_verified',
        'actual_complete_future_contains_selected_outer_angular_repairs_and_entire_heat',
        'actual_flatten_radius_and_right_power_interface_AST_verified',
        'ordinary_t_equals_logR_derivatives'):
        if not stress.bridge['identities'][flag]:raise ValueError('Same actual flatten source missing: '+flag)
        checks['consumed_'+flag]=True
    if not stress.bridge['actual_flatten_power_formula_join']['actual_flatten_power_K4_full_moment_stress_mixed3_pressure_mixed4_AST_join_verified']:
        raise ValueError('Same full flatten/power source join missing')
    checks['same_full_right_power_stress_pressure_join_consumed']=True
    # Source monotonicity: original logistic derivative, reflection and
    # exact flat endpoint branches; tail caps only enclose those functions.
    syntax('flat_pulse_derivatives','_sigma_left','odds','1/(1-x)**2-1/x**2')
    syntax('flat_pulse_derivatives','_sigma_left','e','positive_exp(c,odds)')
    syntax('flat_pulse_derivatives','_sigma_left','value','e/(1+e)')
    syntax('flat_pulse_derivatives','_sigma_left','q','value*(1-value)')
    syntax('flat_pulse_derivatives','_sigma_left','L1','2/(1-x)**3+2/x**3')
    syntax('flat_pulse_derivatives','sigma_jets','left','_sigma_left(c,1-c.mpf([a,b]))')
    derivative=syntax('flat_pulse_derivatives','_sigma_left','derivatives')
    if not isinstance(derivative,ast.List) or ast.dump(derivative.elts[1])!=ast.dump(ast.parse('q*L1',mode='eval').body):
        raise ValueError('Actual positive sigmoid derivative changed')
    fn=asts.method('flat_pulse_derivatives','sigma_jets')
    reflected=ast.parse('rows.append(IntervalTaylor(c,[1-left[0]]+[(-1)**(n+1)*left[n] for n in range(1,5)]))',mode='eval').body
    if not any(ast.dump(node)==ast.dump(reflected) for node in ast.walk(fn)):raise ValueError('Actual sigmoid reflection sign changed')
    for value,expected in ((0,[0,0,0,0,0]),(1,[1,0,0,0,0])):
        actual=sigma_jets(stress.ctx,stress.ctx.mpf(value))
        if any(endpoints(actual[n])!=(expected[n],expected[n]) for n in range(5)):raise ArithmeticError('Actual flat sigma endpoint changed')
    checks['actual_original_sigma_monotonic_positive_derivative_and_reflection_AST_verified']=True
    checks['actual_original_sigma_flat_endpoints_and100_unit_scaling_verified']=True
    if not stress.power.angular.entry.power_source.exit_source.bridge['exact_S_is_positive_source_not_enclosure_endpoint']:
        raise ValueError('Exact positive inverse-Rtail source required')
    checks['exact_positive_S_source_consumed']=True
    logs=flatten_Bmax_binding();asts.hashes.update(logs['input_hashes'])
    return dict(identities=checks,actual_source_AST_bindings=asts.bindings,input_hashes=asts.hashes,
        actual_Bmax_source_log_binding=logs,
        signed_incoming_defect_retained=True,serialized_Xv_box_used_as_defining_value=False,
        Xp_and_Xv_positive_bounds_required=True,continuous_original_domain_used=True,
        source_caps_used_as_fields=False,phase_samples_used_as_proof=False)


@source_precision
def whole_flatten_bounds(stress,assembly):
    h=stress.heat;c=stress.ctx;out=stress.power.outer;r=1-h.mu;b=(1-h.delta)/2
    packet=stress.flatten([-1,1],[0,100]);rows=packet['flatten_similarity_stress_mixed3_factored']
    rawCt=rows['theta']['y0_Z0'];Cz=rows['axial']['y0_Z0']
    Xp=stress.native.inlet.Xp
    equilibrium_gap=(h.mu-h.a)/r;gap_lower=c.mpf(endpoints(equilibrium_gap)[0])
    # A safe bound on the actual signed defect, without cancellation near1.
    defect_upper=c.mpf(endpoints(abs(Xp)+1/r)[1])
    log_actual_memory_upper=c.ln(2*h.k*defect_upper)-13*r/h.mu
    log_memory_cap=c.ln(gap_lower)-1000
    memory_log_gap=log_memory_cap-log_actual_memory_upper
    if endpoints(memory_log_gap)[0]<=0:raise ArithmeticError('Original pulse memory not below certified cap')
    memory_upper=c.mpf(endpoints(c.exp(log_memory_cap))[1])
    W_lower=c.mpf(endpoints(gap_lower-memory_upper)[0])
    cover=[];sigma_t_upper=mp.mpf(0)
    for i in range(8):
        lo=mp.mpf(i)/8;hi=mp.mpf(i+1)/8;jet=sigma_jets(c,c.mpf([lo,hi]))
        upper=endpoints(jet[1]/100)[1];sigma_t_upper=max(sigma_t_upper,upper)
        cover.append(dict(original_x_interval=[lo,hi],actual_sigma_t_enclosure=jet[1]/100))
    slope=c.mpf(sigma_t_upper)
    m_upper=c.mpf(endpoints(2*h.mu+2*c.ln(2)*slope)[1])
    shear_upper=c.mpf(endpoints(2*h.S*c.exp(h.steep.wait+h.steep.Ts+102+out.Lrel)*(1+h.mu+c.ln(2)*slope))[1])
    bracket_lower=c.mpf(endpoints(W_lower-shear_upper)[0])
    # d<0,y<=0,F/f>=1 over original domain.
    K_lower=c.mpf(endpoints(stress.terminal([-1,1])['Kright'])[0])
    theta_lower=c.mpf(endpoints(K_lower*bracket_lower)[0]);axial_upper=c.mpf(endpoints(abs(Cz))[1])
    logparts=stress_source_log_parts(assembly,h,c.mpf(0))['B']
    logparts['steep']=-h.k*h.steep.Ts
    logparts['flatten_power']=(h.a-h.mu)*(100+out.Lrel)
    logparts['waiting_and_current']=2*h.bh
    logBmax=sum(logparts.values(),c.mpf(0))
    margins=dict(mu=h.mu,mu_minus_a=h.mu-h.a,rate=r,delta=h.delta,epsilon=h.eps,
        one_minus_epsilon=1-h.eps,k=h.k,p=h.prate,bh=h.bh,b=b,Lmin=1-h.delta,
        wait=h.steep.wait,Ts=h.steep.Ts,Lrel_minus4=out.Lrel-4,actual_Xp=Xp,actual_Xv=stress.native.Xv,
        equilibrium_gap=gap_lower,signed_memory_log_gap=memory_log_gap,W=W_lower,
        source_K_lower=K_lower,shear_sign=1+h.mu,theta_bracket=bracket_lower,
        theta_stress=theta_lower,axial_absolute_upper=axial_upper,m_upper=m_upper,m_below2=2-m_upper)
    for label,value in margins.items():
        if endpoints(value)[0]<=0 or not mp.isfinite(endpoints(value)[1]):raise ArithmeticError('Flatten cone positive margin failed: '+label)
    if endpoints(h.mu)[1]>=1 or endpoints(h.eps)[1]>=1 or endpoints(h.delta)[1]>=1:
        raise ArithmeticError('Original flatten parameter range failed')
    logcone=c.ln(m_upper)+2*logBmax+2*c.ln(axial_upper/theta_lower);threshold=c.ln(2)
    if endpoints(logcone)[1]>=endpoints(threshold)[0]:raise ArithmeticError('Whole original flatten directional cone failed')
    return dict(positive_margins=margins,actual_whole_source_raw_theta_enclosure=rawCt,actual_whole_source_axial_enclosure=Cz,
        actual_equilibrium_W_gap=equilibrium_gap,original_Xp_enclosure=Xp,original_signed_defect_absolute_upper=defect_upper,
        exact_original_memory_log_upper=log_actual_memory_upper,certified_memory_log_cap=log_memory_cap,
        signed_memory_absolute_upper=memory_upper,continuous_W_uniform_lower=W_lower,
        sigma_t_continuous_source_cover=cover,actual_sigma_t_upper=slope,
        actual_K_uniform_lower=K_lower,exact_source_shear_coefficient_absolute_upper=shear_upper,
        theta_uniform_positive_lower=theta_lower,axial_uniform_absolute_upper=axial_upper,
        actual_kappa_minus2_lower=2*h.mu,actual_kappa_minus2_upper=m_upper,
        actual_Bmax_source_log_parts=logparts,actual_log_Bmax_enclosure=logBmax,
        log_directional_term_upper=logcone,log_threshold=threshold,
        whole_original_t_Z_domain_covered=True,continuous_variable_rate_W_comparison_used=True,
        actual_Xp_Xv_sign_and_signed_pulse_memory_retained=True,
        raw_whole_theta_box_used_as_positivity_proof=False,phase_grid_sampling_used_as_proof=False,
        slope_cover_uses_source_intervals_not_point_sampling=True,
        source_caps_only_define_bounds_not_field_values=True,interval_S_lower_zero_not_used_for_strict_source_sign=True)


@source_precision
def run():
    stress=CompliantFlattenStressC3();assembly=CompliantGlobalPhysicalAssembly()
    if (stress.family,stress.source)!=(assembly.family,assembly.source):raise ValueError('Flatten cone source family differs')
    hashes=dict(stress.hashes);hashes.update(assembly.hashes)
    for stem,gates in (
        ('flatten_stress_C3_check',('actual_original_flatten_similarity_stress_recovered',)),
        ('flatten_physical_C2_check',('actual_regional_physical_flatten_stress_remainder_identity_verified',
            'flatten_power_physical_stress_mixed3_and_remainder_mixed2_join_verified')),
        ('outer_power_cone_check',('outer_power_angular_entry_exit_waiting_collar_two_vector_cone_certified',))):
        name=PREFIX+stem+'.json';receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or any(not receipt[gate] for gate in gates):raise ValueError('Flatten cone prerequisite missing')
        if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(stress.family,stress.source):
            raise ValueError('Flatten cone prerequisite family differs')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Flatten cone source changed: '+source)
        hashes.update(receipt['input_hashes']);hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    proof=flatten_cone_identities();bridge=flatten_cone_source_bridge(stress);bounds=whole_flatten_bounds(stress,assembly)
    hashes.update(bridge['input_hashes']);hashes.update(assembly.dispatch.hashes)
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=stress.family,implicit_source_sha256=stress.source,
        scope='Whole original flatten t[0,100],Z[-1,1]; physical tau>0,r>0,|Z|<1',
        domain=DOMAIN,tail_domain=TAIL_DOMAIN,exact_identities=proof,source_bridge=bridge,bounds=bounds,input_hashes=hashes,
        flatten_cone_certified=True,actual_flatten_theta_stress_positive=True,
        actual_flatten_source_shear_strictly_negative=True,actual_flatten_directional_cone_uniform_margin=True,
        flatten_power_angular_entry_exit_waiting_collar_two_vector_cone_certified=True,
        same_source_flatten_power_similarity_and_physical_joins_consumed=True,
        fixed_positive_viscosity_two_vector_cone_transfer_verified=True,
        flatten_regional_remainder_exact_zero=False,left_pulse_physical_stress_join_verified=False,
        completed_full_tensor_cone_certified=False,whole_outer_cone_certified=False,
        global_admissible_stress_lift_constructed=False,independently_bounded_global_flat_remainder=False,
        physical_energy_integral_certified=False,temporal_recursion=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Whole original100-unit flatten cone and full outer tail generated; left pulse/global/recursion pending',flush=True)
    return result


class CertifiedFlattenPhysical:
    def __init__(self):
        from lei_ren_part1_paper_compliant_flatten_physical_C2 import CompliantFlattenPhysicalC2
        receipt=json.loads((HERE/(PREFIX+'flatten_cone_check.json')).read_bytes())
        if not receipt['all_passed'] or not receipt['flatten_cone_certified']:raise ValueError('Accepted flatten cone required')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Certified flatten source changed: '+source)
        self.field=CompliantFlattenPhysicalC2();self.family,self.source=self.field.family,self.field.source
        if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(self.family,self.source):
            raise ValueError('Flatten physical cone family differs')
        self.receipt=receipt
    def flatten(self,*args,**kwargs):
        point=self.field.flatten(*args,**kwargs)
        point.update(flatten_cone_certified=True,flatten_power_angular_entry_exit_waiting_collar_two_vector_cone_certified=True,
            certified_cone_scope='Original two-vector only; completed diagonal/global cone not certified',
            completed_full_tensor_cone_certified=False)
        return point


if __name__=='__main__':run()
