"""Continuous original pulse-end cone with full nonzero axial shear.

Bounds consume the accepted whole-source sectors and retain their exact
B/D/H/R log recipes. They certify the regional two-vector criterion only.
"""
import ast
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import DOMAIN, SourceAST, source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_compliant_'
TAIL_DOMAIN = 'Rtail*exp(-wait-Ts-106-Lrel)<=R<Rtail*exp(3); Gamma stress exactly zero beyond'
ADMISSIONS = (
    'pulse_end_cone_certified', 'actual_pulse_end_theta_stress_positive',
    'actual_pulse_end_source_shear_strictly_negative',
    'actual_pulse_end_full_axial_shear_kappa_verified',
    'actual_pulse_end_continuous_directional_cone_verified',
    'pulse_end_flatten_power_angular_entry_exit_waiting_collar_two_vector_cone_certified',
    'same_source_pulse_end_flatten_physical_interface_consumed',
    'fixed_positive_viscosity_two_vector_cone_transfer_verified')
FALSE_FLAGS = (
    'completed_full_tensor_cone_certified', 'whole_outer_cone_certified',
    'global_admissible_stress_lift_constructed', 'independently_bounded_global_flat_remainder',
    'physical_energy_integral_certified', 'full_background_NS_validation', 'temporal_recursion')
GATES = {
    'pulse_end_stress_C3_check': ('actual_original_pulse_end_similarity_stress_recovered',
        'whole_original_end_domain_retained', 'all_meridional_cross_terms_retained'),
    'pulse_end_flatten_join_check': ('pulse_end_flatten_full_moment_stress_pressure_functional_join_verified',),
    'pulse_end_physical_C2_check': ('actual_pulse_end_physical_decomposition_constructed',
        'pulse_end_flatten_completed_physical_interface_verified'),
    'pulse_end_support_interfaces_check': (
        'actual_pulse_end_stress3_and_physical_error2_flat_interface_bounds_available',),
    'fifth_axial_jets_check': ('actual_selected_ap_c1_c2_C5_available',),
    'flatten_cone_check': ('flatten_power_angular_entry_exit_waiting_collar_two_vector_cone_certified',)}
RELATIVE_MODES = {
    'theta': {'equilibrium': (0, 0, 0, 0), 'signed_original_memory': (0, 0, 0, 1),
        'meridional_transport': (0, 1, 1, 0), 'radial_shear': (-1, 0, 0, 0)},
    'axial': {'full_energy_and_pressure': (0, 1, 0, 0), 'axial_transport': (0, 0, 1, 0),
        'nonlinear_meridional_transport': (0, 1, 2, 0), 'linear_axial_moment': (0, 0, 1, 0),
        'selected_backward_energy_loss': (0, 1, 2, 0), 'axial_radial_shear': (-1, 0, 1, 0)}}


def current_sources():
    records = {}; hashes = {}; family = None
    for stem in ('pulse_end_stress_C3', 'pulse_mixed_C4', *GATES):
        name = PREFIX+stem+'.json'; raw = (HERE/name).read_bytes(); record = json.loads(raw)
        current = (record['actual_five_defect_family_sha256'], record['implicit_source_sha256'])
        if family is None: family = current
        if current != family: raise ValueError('Pulse cone source family differs: '+stem)
        if stem in GATES and (not record['all_passed'] or any(not record[flag] for flag in GATES[stem])):
            raise ValueError('Pulse cone prerequisite missing: '+stem)
        for source, digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest() != digest:
                raise ValueError('Pulse cone source changed: '+source)
        hashes.update(record['input_hashes']); hashes[name] = hashlib.sha256(raw).hexdigest()
        records[stem] = record
    if records['pulse_end_stress_C3']['domain'] != DOMAIN:
        raise ValueError('Whole original pulse end required')
    whole = records['pulse_end_stress_C3']['whole_original_end']
    for key, expected in (('s', (-4, 0)), ('Z', (-1, 1))):
        c = MPIntervalContext(); c.dps = 240
        if endpoints(read_interval(c, whole[key])) != expected:
            raise ValueError('Whole source domain differs: '+key)
    return records, hashes, family


def pulse_end_cone_identities():
    asts = SourceAST(); proofs = {}
    z, mu, delta, bh, bh_s = s.symbols('Z mu delta Bhat Bhat_s', real=True)
    R, B, D, H = s.symbols('R B D H', positive=True)
    C = 1/(1+z*z); r = 1-mu; b0 = (1-delta)/2
    L = 1-delta*z*z; g = (mu-delta/2)/r; a = 2+2*mu
    def zero(name, expression):
        if s.simplify(s.expand_power_exp(expression)) != 0:
            raise ArithmeticError('Pulse cone identity failed: '+name)
        proofs[name] = True
    equilibrium = asts.evaluate(asts.expression('pulse_end_stress_C3', 'pulse_coefficients', 'equilibrium'),
        dict(C=C, mu=mu, delta=delta, r=r, z=z, b=b0, L=L, axial_derivative=lambda v:s.diff(v,z)))
    zero('actual_equilibrium_positive_rewrite',
        equilibrium-(C*g+2*b0*z*z/(r*(1+z*z)**2))/L)
    fn = asts.method('pulse_end_stress_C3', 'pulse_coefficients')
    for label in ('theta', 'axial'):
        node = next(n.value for n in ast.walk(fn) if isinstance(n,ast.Assign)
            and any(ast.unparse(t)==label for t in n.targets))
        names = []
        for key, value in zip(node.keys, node.values):
            name = ast.literal_eval(key); names.append(name)
            mode = next(ast.literal_eval(k.value) for k in value.keywords if k.arg=='mode')
            relative = tuple(s.Rational(v)-w for v,w in zip(mode,(s.Rational(1,2),1,0,0)))
            if relative != RELATIVE_MODES[label][name]: raise ValueError('Actual pulse mode changed')
            zero('actual_'+label+'_'+name+'_relative_factor',
                R**s.Rational(mode[0])*B**mode[1]*D**mode[2]*H**mode[3]/s.sqrt(2)
                /(s.sqrt(R/2)*B)-R**relative[0]*B**relative[1]*D**relative[2]*H**relative[3])
        if set(names) != set(RELATIVE_MODES[label]): raise ValueError('Full ten stress sectors required')
    path = HERE/'lei_ren_part1_paper_mp_stress.py'
    asts.hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    original = next(n for n in ast.walk(ast.parse(path.read_text(encoding='utf8')))
        if isinstance(n,ast.FunctionDef) and n.name=='evaluate_mp_stress')
    env = dict(a=B*C, ay=-(s.Rational(1,2)+mu)*B*C,
        by=B*C*D*(bh_s-(s.Rational(1,2)+mu)*bh), root=s.sqrt(2*R), R=R)
    shear = {}
    for target in ('Stheta', 'Sz'):
        node = next(n.value for n in ast.walk(original) if isinstance(n,ast.Assign)
            and any(ast.unparse(t)==target for t in n.targets))
        if isinstance(node,ast.IfExp): node = node.body
        shear[target] = asts.evaluate(node,env)
    F = B*C/s.sqrt(2*R); A = a*F
    sigma = D*(bh_s-(s.Rational(1,2)+mu)*bh)/(1+mu)
    zero('actual_original_negative_theta_shear', shear['Stheta']+A)
    zero('actual_original_full_axial_shear_ratio', shear['Sz']-A*sigma)
    kappa = -(shear['Stheta']**2+shear['Sz']**2)/(F*shear['Stheta'])
    zero('actual_full_kappa_a_plus_b_squared_over_a', kappa-a*(1+sigma**2))
    zero('actual_kappa_minus2_without_large_subtraction', kappa-2-(2*mu+a*sigma**2))
    theta, axial, sig, Q, AA = s.symbols('theta axial sigma Q A', real=True)
    dot = Q*theta*(-AA)+Q*axial*AA*sig
    cross = Q*theta*(-AA*sig)+Q*axial*(-AA)
    zero('actual_negative_dot_normalization', -dot-Q*AA*(theta-axial*sig))
    zero('actual_cross_direction_normalization', -cross-Q*AA*(theta*sig+axial))
    km = s.symbols('kappa_minus2',real=True)
    zero('actual_directional_cone_normalization',
        2*dot**2-km*cross**2-(Q*AA)**2*(2*(theta-axial*sig)**2-km*(theta*sig+axial)**2))
    nu, fac, ff, st, sz = s.symbols('nu lambda_factor F Stheta Sz',positive=True)
    zero('actual_viscosity_normalized_physical_kappa',
        -((nu*fac*st)**2+(nu*fac*sz)**2)/(nu*(fac*ff)*(nu*fac*st))
        +(st**2+sz**2)/(ff*st))
    zero('actual_physical_directional_margin_common_positive_factor',
        2*((nu*fac)**2*dot)**2-km*((nu*fac)**2*cross)**2
        -(nu*fac)**4*(2*dot**2-km*cross**2))
    for target, wanted in (
        ('q', 'IntervalTaylor(c,[1+z0**2,2*z0,1,0,0,0])'),
        ('logR', 'self.assembly.logRp+13/mu+v'),
        ('logD', 'self.native.logE'),
        ('logH', '-13*(1-mu)/mu-(1-mu)*v')):
        asts.expression('pulse_end_stress_C3','end',target,wanted=wanted)
    asts.expression('pulse_end_stress_C3','end','logB',wanted=
        "dict(logPstar=self.assembly.logP,actual_log_inlet_U=c.ln(self.flatten.native.U),inverse_mu=-13/(2*mu),finite=-13-(c.mpf('.5')+mu)*v)")
    asts.expression('pulse_end_stress_C3','end','Bh[j]',
        wanted='Cj*(beta[j]*math.factorial(j))',augmented=True)
    # Exact original production-radius differences, never subtract giant log boxes.
    tree = ast.parse((HERE/(PREFIX+'global_physical_assembly.py')).read_text(encoding='utf8'))
    radius_env = {}
    for node in tree.body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ('MICRO','POST','PULSE'):
                    radius_env[target.id] = ast.literal_eval(node.value)
    asts.replay('global_physical_assembly','radius',radius_env)
    origin, v, wait, Ts, length = s.symbols('logRp s wait Ts Lrel',real=True)
    assembly = SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify),logRp=origin,params=SimpleNamespace(mu=mu))
    radial = radius_env['radius']
    end = radial(assembly,'pulse_end',v,{},None)[0]
    flatten = radial(assembly,'flatten',0,{},None)[0]
    heat = SimpleNamespace(steep=SimpleNamespace(wait=wait,Ts=Ts,outer=SimpleNamespace(Lrel=length)))
    tail = radial(assembly,'heat_collar',0,{},heat)[0]
    zero('actual_pulse_end_same_flatten_radius_plus_s',end-flatten-v)
    zero('actual_pulse_end_left_original_tail_radius_offset',end.subs(v,-4)-tail+wait+Ts+106+length)
    return dict(identities=proofs,input_hashes=asts.hashes,actual_source_AST_bindings=asts.bindings,
        equilibrium_lower_formula='Etheta >= g/2; C in[1/2,1], Z^2>=0, 0<L<=1, b0>=0',
        equilibrium_upper_formula='Etheta <= (g+2*b0/r)/(1-delta)',
        source_shear='S=A*(-1,sigma), A=(2+2mu)*C*B/sqrt(2R)>0',
        actual_kappa='kappa-2=2mu+(2+2mu)*sigma^2',
        physical_kappa='-|Sphys|^2/(nu*Fphys*Sphys_theta), fixed nu>0',
        positive_source_factors_are_exact_exponentials=True,
        normalized_coefficients_restored_with_all_parent_modes_and_logs=True,
        full_nonzero_axial_shear_retained=True)


def absolute_bound(c, value):
    return c.mpf(max(abs(v) for v in endpoints(value)))


@source_precision
def whole_pulse_end_bounds(records):
    c = MPIntervalContext(); c.dps = 240
    selected = records['pulse_mixed_C4']
    mu = read_interval(c,selected['selected_mu']); delta = read_interval(c,selected['selected_delta'])
    r = 1-mu; a = 2+2*mu; bp = c.mpf('.5')+mu; b0 = (1-delta)/2
    g = (mu-delta/2)/r
    def positive(name, value):
        lo, hi = endpoints(value)
        if lo<=0 or not mp.isfinite(hi): raise ArithmeticError('Pulse end positive margin failed: '+name)
    for name, value in dict(mu=mu,delta=delta,rate=r,one_minus_delta=1-delta,b0=b0,g=g).items():
        positive(name,value)
    g_lower = c.mpf(endpoints(g)[0])
    theta_lower = g_lower/4
    eq_upper = c.mpf(endpoints((g+2*b0/r)/(1-delta))[1])
    theta_upper = c.mpf(endpoints(eq_upper+g_lower/8)[1])
    epsilon = c.exp(-1000)
    axial_cap = theta_lower*epsilon
    sigma_cap = theta_lower/theta_upper*epsilon
    whole = records['pulse_end_stress_C3']['whole_original_end']
    sectors = whole['full_meridional_stress_log_sectors']
    eqlogs = sectors['theta']['equilibrium']['exact_source_log_parts']
    logBparts = {key:read_interval(c,eqlogs[key]) for key in
        ('logPstar','actual_log_inlet_U','inverse_mu','finite')}
    logB = sum(logBparts.values(),c.mpf(0))
    logR = 2*read_interval(c,eqlogs['source_logR'])
    logD = read_interval(c,sectors['theta']['meridional_transport']['exact_source_log_parts']['selected_log_end_scale'])
    logH = read_interval(c,sectors['theta']['signed_original_memory']['exact_source_log_parts']['signed_original_memory_log'])
    logs = (logR,logB,logD,logH); rows = {}; margins = {}
    for label, source_rows in sectors.items():
        if set(source_rows) != set(RELATIVE_MODES[label]): raise ValueError('Current full sector set changed')
        for name, sector in source_rows.items():
            expected_mode = tuple(v+w for v,w in zip(RELATIVE_MODES[label][name],(.5,1,0,0)))
            if tuple(sector['mode']) != expected_mode: raise ValueError('Current sector normalization changed')
            if label=='theta' and name=='equilibrium': continue
            coefficient = absolute_bound(c,read_interval(c,sector['full_stress_mixed3_coefficient_enclosures']['s0_Z0']))
            if endpoints(coefficient)[1]==0: raise ValueError('Expected nonzero source coefficient bound')
            relative = RELATIVE_MODES[label][name]
            actual_log = sum((power*log for power,log in zip(relative,logs) if power),c.mpf(0))+c.ln(coefficient)
            target_log = c.ln(g_lower/32) if label=='theta' else c.ln(theta_lower)-1000-c.ln(6)
            gap = target_log-actual_log; key = label+'_'+name
            positive(key,gap); margins[key]=gap
            rows[key]=dict(actual_signed_coefficient=read_interval(c,sector['full_stress_mixed3_coefficient_enclosures']['s0_Z0']),
                coefficient_absolute_upper=coefficient,relative_log_mode=relative,
                actual_relative_log_absolute_upper=actual_log,certified_log_cap=target_log,positive_log_gap=gap)
    Bh = [read_interval(c,whole['source_formal_Bhat_rows'][k]['coefficients'][0]) for k in (0,1)]
    sigma_coefficient = (absolute_bound(c,Bh[1])+bp*absolute_bound(c,Bh[0]))/(1+mu)
    actual_sigma_log = logD+c.ln(sigma_coefficient)
    sigma_target_log = c.ln(theta_lower/theta_upper)-1000
    margins['full_axial_shear_ratio_log_gap']=sigma_target_log-actual_sigma_log
    positive('full_axial_shear_ratio_log_gap',margins['full_axial_shear_ratio_log_gap'])
    dot_lower = theta_lower-axial_cap*sigma_cap
    cross_upper = theta_upper*sigma_cap+axial_cap
    km_lower = 2*mu
    km_upper = 2*mu+a*sigma_cap**2
    directional_margin = 2*dot_lower**2-km_upper*cross_upper**2
    algebraic = dict(theta_lower=theta_lower,theta_upper=theta_upper,
        equilibrium_lower_slack=g_lower/2-3*g_lower/32-theta_lower,
        equilibrium_upper_slack=g_lower/8-3*g_lower/32,
        negative_dot_lower=dot_lower,cross_absolute_upper=cross_upper,
        kappa_minus2_lower=km_lower,kappa_minus2_upper=km_upper,
        kappa_minus2_below2=2-km_upper,directional_margin=directional_margin)
    for name,value in algebraic.items(): positive(name,value)
    return dict(positive_log_margins=margins,positive_algebraic_margins=algebraic,
        mu=mu,delta=delta,equilibrium_gap=g,source_log_recipes=dict(B=logBparts,R=logR,D=logD,H=logH),
        normalized_stress_sector_bounds=rows,
        axial_uniform_absolute_upper=axial_cap,sigma_uniform_absolute_upper=sigma_cap,
        full_axial_shear_coefficient_upper=sigma_coefficient,
        actual_full_axial_shear_ratio_log_upper=actual_sigma_log,certified_shear_ratio_log_cap=sigma_target_log,
        normalized_cone='2*(theta-axial*sigma)^2-(kappa-2)*(theta*sigma+axial)^2>0',
        continuous_whole_original_s_Z_domain_covered=True,
        actual_source_log_recipes_consumed_without_giant_log_subtraction=True,
        signed_original_memory_and_full_ten_sectors_retained=True,
        full_axial_shear_not_set_to_zero=True,
        phase_samples_used_as_proof=False,source_caps_used_as_defining_field_values=False,
        support_receipt_scope='Whole-source coefficient difference bounds; absolute factors restored via parent recipes',
        global_temporal_flatness_not_inferred=True)


@source_precision
def run():
    records,hashes,family = current_sources()
    proof = pulse_end_cone_identities(); hashes.update(proof['input_hashes'])
    bounds = whole_pulse_end_bounds(records)
    hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result = dict(actual_five_defect_family_sha256=family[0],implicit_source_sha256=family[1],
        domain=DOMAIN,tail_domain=TAIL_DOMAIN,exact_source_cone_proof=proof,bounds=bounds,input_hashes=hashes,
        current_selected_C5_source_check_directly_consumed=True,
        regional_two_vector_cone_distinct_from_completed_tensor_cone=True,
        **{flag:True for flag in ADMISSIONS},**{flag:False for flag in FALSE_FLAGS})
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Whole original pulse-end cone with full axial shear generated; gap/main/global/recursion pending',flush=True)
    return result


class CertifiedPulseEndPhysical:
    """Add current regional admission to the unchanged physical source evaluator."""
    def __init__(self,cells=64):
        from lei_ren_part1_paper_compliant_pulse_end_physical_C2 import CompliantPulseEndPhysicalC2
        self.receipt = json.loads(Path(__file__).with_name(PREFIX+'pulse_end_cone_check.json').read_bytes())
        if not self.receipt['all_passed'] or not self.receipt['pulse_end_cone_certified']:
            raise ValueError('Accepted pulse-end cone required')
        for source,digest in self.receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:
                raise ValueError('Certified pulse-end source changed: '+source)
        self.physical = CompliantPulseEndPhysicalC2(cells=cells)
        if (self.physical.family,self.physical.source) != (
            self.receipt['actual_five_defect_family_sha256'],self.receipt['implicit_source_sha256']):
            raise ValueError('Certified pulse-end physical source differs')

    def end(self,*args,**kwargs):
        result = self.physical.end(*args,**kwargs)
        result.update(**{flag:True for flag in ADMISSIONS},**{flag:False for flag in FALSE_FLAGS})
        return result


if __name__=='__main__': run()
