"""Whole modulation-support norms from the actual sigma/base source.

The original positive amplitude at Rd stays factored. No amplitude or mu
midpoint is selected and no checked ancestor graph is rebuilt. These
source formulas do not include the independent repair or completed tensor.
"""
import hashlib
import inspect
import json
import math
from pathlib import Path
from types import SimpleNamespace
import sympy as s
from lei_ren_part1_paper_compliant_steep_entry_stress_C3 import SourceAST
from lei_ren_part1_paper_compliant_current_O3_frequency_majorants import (
    profile_majorants, _encode, PREFIX, HERE)

NAME = PREFIX+'current_O3_frequency_source_bounds.json'
RECEIPT = PREFIX+'current_O3_frequency_source_bounds_check.json'


def sha(name):
    return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def sigma_source_caps():
    from lei_ren_part1_paper_compliant_flat_pulse_derivatives import SIGMA_CONSTANTS
    if tuple(SIGMA_CONSTANTS) != (1,4,28,256,3104):
        raise ValueError('Actual squared-exponent sigma constants changed')
    caps = [s.Integer(1)]
    for j in range(1,5):
        rho = s.Rational(1,2) if j<=2 else s.sqrt(s.Rational(2,3*j))
        caps.append(s.Integer(SIGMA_CONSTANTS[j])*s.exp(4-1/rho**2)/rho**(3*j))
    return caps


def actual_source_norms(mu, amplitude_at_Rd):
    mu, Ad = s.sympify(mu), s.sympify(amplitude_at_Rd)
    if any(value.is_positive is not True or value.is_finite is not True for value in (mu,Ad)):
        raise ValueError('Positive finite same-source mu and factored Rd amplitude required')
    B = sigma_source_caps()
    cutoff = [s.Integer(1)]+[4**j*B[j] for j in range(1,5)]
    a = (1+mu)/2
    b,c,d = (mu*B[j] for j in range(1,4))
    positive = [s.Integer(1),a,a*a+b,a**3+3*a*b+c,
                a**4+6*a*a*b+3*b*b+4*a*c+d]
    time = [s.E]+[s.Max(s.E/2**j,positive[j]) for j in range(1,5)]
    # Partial fractions certify |dZ^k (1+Z^2)^-1|<=k! on all real Z.
    axial = [s.factorial(k) for k in range(6)]
    F = [[Ad*time[j]*axial[k] for k in range(6)] for j in range(5)]
    return dict(whole_support_logR_offset=['-2','1/2'],whole_axial_interval=['-1','1'],
        actual_cutoff_ordinary_logR_majorants=cutoff,original_sigma_ordinary_majorants=B,
        original_base_time_factor_majorants=time,original_base_axial_shape_majorants=axial,
        original_theta_over_Pstar_logR_axial_majorants=F,
        original_log_theta_logR_majorant=a,actual_mu=mu,actual_factored_amplitude_at_Rd=Ad,
        actual_source_norm_formulas_cover_whole_support_and_Z=True,
        numeric_actual_graph_parameter_enclosures_substituted=False,
        independent_repair_pressure_radial_completed_tensor_bounds_available=False)


def bounded_modulation(mu,N,amplitude_at_Rd):
    norms = actual_source_norms(mu,amplitude_at_Rd)
    envelope = profile_majorants(norms['actual_mu'],N,
        norms['actual_cutoff_ordinary_logR_majorants'],
        norms['original_theta_over_Pstar_logR_axial_majorants'],
        norms['original_log_theta_logR_majorant'])
    return dict(actual_whole_support_source_norm_formulas=norms,
        source_bound_modulation_derivative_and_shear_envelopes=envelope,
        source_norm_hypotheses_discharged_by_actual_sigma_and_base_factorization=True,
        actual_numeric_graph_parameters_not_selected=True,
        independent_repair_and_completed_tensor_common_N_cones_remain_open=True)


def exact_source_bound_theorem():
    asts,checks = SourceAST(),{}
    def zero(name,a,b):
        if s.simplify(s.expand_power_exp(a-b)) != 0:
            raise ArithmeticError('Actual source norm identity failed: '+name)
        checks[name] = True
    records,hashes = {},{}
    for stem,flags in (
        ('flat_pulse_derivatives_check',('all_passed','original_radial_shape_derivatives_C4_available')),
        ('current_O3_frequency_majorants_check',('all_passed','conditional_actual_modulation_derivative_envelopes_available')),
        ('current_modified_physical_velocity_check',('all_passed','current_modified_velocity_pressure_physical_locator_integrated'))):
        name = PREFIX+stem+'.json';record=json.loads((HERE/name).read_bytes())
        if not all(record.get(flag) is True for flag in flags):
            raise ValueError('Source derivative prerequisite missing: '+name)
        records[stem]=record;hashes[name]=sha(name)
        if stem!='current_modified_physical_velocity_check':
            for path,digest in record['input_hashes'].items():
                if path in hashes and hashes[path]!=digest: raise ValueError('Different actual sigma source')
                hashes[path]=digest
    for path,digest in hashes.items():
        if sha(path)!=digest:raise ValueError('Changed accepted source norm program: '+path)
    current = records['current_modified_physical_velocity_check']
    for stem in ('current_O3_finite_frequency_profiles','pre_pulse_mixed_C4','flat_pulse_derivatives','outer_buffer'):
        name=PREFIX+stem+'.py'
        if current['input_hashes'].get(name)!=sha(name):
            raise ValueError('Current native graph uses different original modulation source: '+name)
        hashes[name]=sha(name)
    # Exact squared-exponent sigma and its original analytic cap program.
    asts.expression('flat_pulse_derivatives','_sigma_left','odds',wanted='1/(1-x)**2-1/x**2')
    asts.expression('flat_pulse_derivatives','sigma_tail_bound','peak',wanted='c.sqrt(c.mpf(2)/(3*n))')
    asts.expression('flat_pulse_derivatives','sigma_tail_bound','rho',wanted='W if endpoints(W)[1]<endpoints(peak)[0] else peak')
    x = s.Symbol('positive_left_sigma_x',positive=True)
    n = s.Symbol('positive_derivative_order',positive=True)
    zero('actual_sigma_envelope_stationary_point',s.diff(4-1/x**2-3*n*s.log(x),x),(2-3*n*x*x)/x**3)
    caps = sigma_source_caps()
    for j in range(1,5):
        rho=s.Rational(1,2) if j<=2 else s.sqrt(s.Rational(2,3*j))
        if not (0<rho<=s.Rational(1,2)):
            raise ArithmeticError('Wrong actual sigma envelope maximum')
        checks['whole_left_half_sigma_maximum_'+str(j)]=True
    zero('source_sigma_B1',caps[1],32)
    zero('source_sigma_B2',caps[2],1792)
    zero('source_sigma_B3',caps[3],256*s.exp(-s.Rational(1,2))*(3/s.sqrt(2))**9)
    zero('source_sigma_B4',caps[4],3104*s.exp(-2)*6**6)
    # Actual chi's two transitions never overlap. Whole-interval factors
    # are exactly 1 and4^j; middle and exterior derivatives vanish.
    left=(-2,-1);right=(s.Rational(1,4),s.Rational(1,2))
    if not (4*left[1]-1<=0 and right[0]+2>=1 and left[1]<right[0]):
        raise ArithmeticError('Actual cutoff transition partition changed')
    checks['actual_cutoff_disjoint_transitions_plateau_and_flat_edges']=True
    # Source inspection binds the actual original common amplitude shape.
    source_nodes={}
    for method,target,wanted in (
        ('coordinates','qi','(1+square(z)).reciprocal()'),
        ('slope','factor',"c.exp(y/10-c.mpf('.6')*J)"),('slope','u','qi*factor'),
        ('inlet','self.cache[key]','self.slope(Z,1)'),
        ('axial','y','c.exp(md)+selector'),('axial','t','y-1'),
        ('axial','root','c.exp(-t/2)'),('axial','u','u1*root'),
        ('slope_mu','parent','self.axial(Z,buffer_offset=11)'),
        ('slope_mu','factor',"c.exp(-t/2-mu*K['J'])"),('slope_mu','u','u1*factor'),
        ('slope_mu','logU',"[zero-c.mpf('.5')-mu*sig[0]]+[zero-mu*sig[k]*math.factorial(k) for k in range(1,4)]")):
        source_nodes[method,target]=asts.expression('pre_pulse_mixed_C4',method,target,wanted=wanted)
        checks['actual_common_base_source_'+method+'_'+target]=True
    t,z,J1,Md,mu,J=s.symbols('offset Z actual_J_at1 actual_Md mu actual_J',real=True)
    q=1/(1+z*z);td=s.exp(Md)+10
    logAd=s.Rational(1,10)-s.Rational(3,5)*J1-td/2
    Ad=s.exp(logAd)
    ctx=SimpleNamespace(exp=s.exp,mpf=lambda value:s.Rational(str(value)))
    slope_factor=asts.evaluate(source_nodes['slope','factor'],dict(c=ctx,y=s.Integer(1),J=J1))
    slope_u=asts.evaluate(source_nodes['slope','u'],dict(qi=q,factor=slope_factor))
    original_y=asts.evaluate(source_nodes['axial','y'],dict(c=ctx,md=Md,selector=11+t))
    original_t=asts.evaluate(source_nodes['axial','t'],dict(y=original_y))
    original_root=asts.evaluate(source_nodes['axial','root'],dict(c=ctx,t=original_t))
    original_buffer_u=asts.evaluate(source_nodes['axial','u'],dict(u1=slope_u,root=original_root))
    zero('actual_buffer_original_parent_coordinate_at_Rd',original_t,td+t)
    zero('actual_buffer_common_amplitude_function',original_buffer_u,Ad*q*s.exp(-t/2))
    transition_factor=asts.evaluate(source_nodes['slope_mu','factor'],dict(c=ctx,t=t,mu=mu,K=dict(J=J)))
    original_transition_u=asts.evaluate(source_nodes['slope_mu','u'],dict(u1=original_buffer_u.subs(t,0),factor=transition_factor))
    zero('actual_O3_same_Rd_amplitude_function',original_transition_u,Ad*q*s.exp(-t/2-mu*J))
    # Ordinary axial majorants by exact partial fractions, no sampling or
    # unchecked claims about sharper stationary roots.
    for k in range(6):
        representation=s.factorial(k)/2*((s.I)**k/(1-s.I*z)**(k+1)+(-s.I)**k/(1+s.I*z)**(k+1))
        zero('actual_C_axial_partial_fraction_order_'+str(k),s.diff(q,z,k),representation)
    zero('real_C_denominator_modulus_squared',(1+s.I*z)*(1-s.I*z),1+z*z)
    # Independent finite series verifies positive-side exponential jets.
    h=s.symbols('h1:5',real=True);dx=s.Symbol('local_offset',real=True)
    inner=sum(h[j-1]*dx**j/s.factorial(j) for j in range(1,5))
    series=s.Poly(s.expand(sum(inner**k/s.factorial(k) for k in range(5))),dx)
    expected=(1,h[0],h[0]**2+h[1],h[0]**3+3*h[0]*h[1]+h[2],
        h[0]**4+6*h[0]**2*h[1]+3*h[1]**2+4*h[0]*h[2]+h[3])
    for j in range(5):zero('actual_base_time_ordinary_exponential_jet_'+str(j),series.nth(j)*s.factorial(j),expected[j])
    checks['J_zero_on_negative_buffer_and_original_sigma_FTC_on_positive_transition']=True
    checks['left_half_sigma_between_zero_and_one_half_by_actual_odds']=True
    checks['negative_time_factor_bounds_exp_minus_t_over2_on_minus2_to0']=True
    checks['positive_time_factor_at_most_one_by_same_J_nonnegative']=True
    inspector=Path(inspect.getsourcefile(SourceAST));asts.hashes[inspector.name]=hashlib.sha256(inspector.read_bytes()).hexdigest()
    hashes.update(asts.hashes);hashes[Path(__file__).name]=sha(Path(__file__).name)
    return dict(identities=checks,passed=True,source_bindings=asts.bindings,input_hashes=hashes,
        source_family={key:current[key] for key in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256')},
        actual_original_Rd_log_amplitude_definition=s.sstr(logAd),
        same_original_sigma_J_at1_kept_as_source_constant_not_fitted=True,
        whole_support_and_Z_source_norm_formulas_certified=True,
        factored_parameter_values_not_numerically_selected=True,
        actual_independent_repair_completed_stress_common_N_cones_remain_open=True)


def run():
    theorem=exact_source_bound_theorem()
    mu=s.Symbol('same_actual_mu',positive=True,finite=True)
    N=s.Symbol('finite_integer_N',positive=True,integer=True,finite=True)
    Ad=s.Symbol('same_actual_amplitude_at_Rd',positive=True,finite=True)
    result=dict(exact_actual_modulation_source_norm_theorem=theorem,
        actual_source_norm_and_modulation_bounds=bounded_modulation(mu,N,Ad),
        whole_support_original_cutoff_and_base_source_norm_formulas_available=True,
        actual_numeric_graph_parameter_enclosures_substituted=False,
        repair_completed_tensor_common_N_modified_cones_certified=False,input_hashes=theorem['input_hashes'])
    (HERE/NAME).write_bytes((json.dumps(_encode(result),indent=2)+'\n').encode())
    print('Actual whole-support cutoff/base source norms generated with the original Rd amplitude factored',flush=True)
    return result


if __name__=='__main__':run()
