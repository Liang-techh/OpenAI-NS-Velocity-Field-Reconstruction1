"""Original normalized preheat pressure point jets with separate error budgets.

Actual defining quadrature supplies a point approximation; independent
directed rectangles enclose its coefficient and arithmetic error. A new
source-density envelope bounds the finite-end and all late original stages.
Those stages remain exact in the source frame and nonzero in the budget.
Physical Pstar² amplification is retained; no native phase or full oracle
is installed by this normalized pressure service.
"""
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_source_parameter_frame as source
from lei_ren_part1_paper_interval_outer_slope_field import transition_integrals
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE,PREFIX,sha=source.HERE,source.PREFIX,source.sha
NAME=PREFIX+'current_original_pressure_point_jets.json'
RECEIPT=PREFIX+'current_original_pressure_point_jets_check.json'
GATE='original_normalized_P0_point_Z_jets_with_directed_coefficient_and_all_late_source_error_bounds'


def late_source_error_theorem(frame):
    """Bind every post-yd density to one integrable original envelope."""
    pressure_receipt=source.inertial.pressure.RECEIPT
    admitted=json.loads((HERE/pressure_receipt).read_bytes())
    identification=admitted['original_pressure_function_identification']
    master=identification['original_master_source']
    assert admitted['all_passed'] and identification['passed'] and master['passed']
    assert identification['same_raw_waiting_equation_and_unique_positive_root']
    assert identification['same_original_cutoff_primitive_and_stage_offsets']
    assert master['exact_cutoff']=='sigma(x)=edge(x)/(edge(x)+edge(1-x)), edge(x)=exp(-1/x^2) for x>0; otherwise0'
    assert master['exact_primitive']=='J(t)=integral_0^t sigma(v)dv; J(1)=1/2 by symmetry'
    phi_binding=source.source_assignment(PREFIX+'collar_Gamma_C4.py','phi_jets','log_upper','-4/D**2')
    assert frame.delta_branch_proof['passed'] and frame.definitions['Md']==40
    operator=frame.owner.pressure;p=operator.partition
    z,t,mu,delta,epsilon=(p[key] for key in ('z','t','mu','delta','epsilon'))
    q=p['q'];J=p['J'];sigma=p['sigma'];phi=p['phi'];r=1-mu;k=1-delta/2
    A=p['Tw']+13/mu+100+p['L'];B=-mu-2*mu*(A+1)-r-2*p['Ts']-2+k
    factor=q*q/4;K=1-epsilon*(1-sigma(t)+sigma(t)*phi(t))
    expected={
        'slope_transition_mu':s.exp(-2*mu*J(t)),
        'power_buffer':s.exp(-mu-2*mu*t),
        'pulse_reserved':s.exp(-mu-2*mu*(p['Tw']+t)),
        'z_flatten':s.exp(-mu-2*mu*(p['Tw']+13/mu+t))*(q/2)**(2*sigma(t/100)),
        'power_buffer_rel':factor*s.exp(-mu-2*mu*(p['Tw']+13/mu+100+t)),
        'steep_transition_in':factor*s.exp(-mu-2*mu*(A+t)-2*r*J(t)),
        'steep_power':factor*s.exp(-mu-2*mu*(A+1)-r-2*t),
        'steep_transition_out':factor*s.exp(-mu-2*mu*(A+1)-r-2*p['Ts']-2*t+2*k*J(t)),
        'waiting':factor*s.exp(B-delta*t),
        'heat_collar':factor*s.exp(B-delta*(p['W']+t))*K*K/(1-epsilon)**2,
        'exterior_power_tail':factor*s.exp(B-delta*(p['W']+t))/(1-epsilon)**2}
    origins=dict(zip(expected,('d','w','p','v','f','rel','s','q','t','tail','tail')))
    assert set(expected)==set(operator.original_densities)-{'reference_extension','slope_transition_ref','axial_turnoff'}
    identities={}
    for name,value in expected.items():
        global_y=p['offsets'][origins[name]]+t
        envelope=s.exp(s.Rational(3,5)-p['yd']-(global_y-p['yd']))/(2*q*q)
        actual=operator.original_densities[name]/envelope
        difference=s.expand_log(s.expand_power_exp((actual-value).rewrite(s.exp)),force=True)
        assert s.simplify(difference)==0,name
        identities[name]=True
    removed=operator.original_densities['axial_turnoff'].subs(t,p['yd']-1+t)
    removed_expected=s.exp(s.Rational(3,5)-p['yd']-t)/(2*q*q)
    assert s.simplify(removed-removed_expected)==0
    assert s.integrate(s.exp(-t),(t,0,s.oo))==1
    x=s.Symbol('nonnegative_abs_Z',real=True)
    assert s.expand((1+x*x)-2*x-(x-1)**2)==0
    assert s.expand(q*q-4*z*z-(1-z*z)**2)==0
    assert s.cancel(s.diff((-4+20*x)/(1+x)**2,x)-(28-20*x)/(1+x)**3)==0
    # q²/4<=1; 0<=J(t)<=t; r>0; 0<k<=1; W>=0.
    # Every displayed exponent is nonpositive. In the out-transition,
    # -2t+2kJ(t)<=-delta*t. K in[1-epsilon,1] on the collar.
    # Therefore the global density <= envelope/(1-epsilon)².
    # Integrating the envelope from yd to infinity gives exp(3/5-yd)/(2q²).
    # The actual fourteen stages form the same contiguous original y line.
    density_jet_checks=0
    for name,density in operator.original_densities.items():
        if name in source.inertial.pressure.BETA2:
            factors=(-4*z/q,(-4+20*z*z)/q**2)
        elif name in source.inertial.pressure.BETA0:factors=(0,0)
        else:
            assert name=='z_flatten';v=sigma(t/100)-1
            factors=(4*z*v/q,4*v/q+8*z*z*v*(2*sigma(t/100)-3)/q**2)
        for order,factor_Z in enumerate(factors,1):
            assert s.simplify(s.diff(density,z,order)/density-factor_Z)==0,(name,order)
            density_jet_checks+=1
    return dict(passed=True,original_late_density_envelope_identities=identities,
        removed_finite_axial_tail_identity_and_unit_exponential_integral_checked=True,
        source_inequality_bindings=dict(receipt=pressure_receipt,receipt_sha256=sha(pressure_receipt),
            original_cutoff_and_primitive=master,original_phi_nonpositive_exponent_binding=phi_binding,
            raw_waiting_root_is_same_unique_positive_source_root=True,
            native_mu_delta_epsilon_hypotheses_from_Md40_branch_proof=frame.delta_branch_proof,
            source_cutoff_ratio_is_in_zero_one_and_primitive_is_integral_of_that_ratio=True,
            phi_is_exp_of_nonpositive_exponent_with_flat_zero_extension=True,
            elementary_q_squared_and_beta2_factor_inequalities_checked=True),
        exact_density_Z_factor_identities=density_jet_checks,
        source_domain='real |Z|<=1; t>=0; 0<mu<1; 0<delta<1; 0<epsilon<1/2; raw W>=0',
        original_step_primitive_bounds='0<=sigma<=1, 0<=phi<=1, 0<=J(t)<=t; fixed source bounds',
        amplitude_ratio_bound='late density/envelope <= (1-epsilon)^-2 <4',
        global_integrable_envelope='exp(3/5-yd)*exp(-(y-yd))/(2*(1+Z²)^2), y>=yd',
        derivative_factor_proof=dict(beta2=[1,2,4],beta0=[1,0,0],flatten=[1,2,10],
            first='4|Z|/(1+Z²)<=2 by (|Z|-1)²>=0',
            beta2_second='(-4+20x)/(1+x)² increases from -4 to4 for x=Z² in[0,1]; derivative=(28-20x)/(1+x)^3>0',
            flatten_second='|sigma-1|<=1, |2sigma-3|<=3, Z²/(1+Z²)²<=1/4: 4+24/4=10'),
        extended_axial_removed_tail='exp(3/5-yd)/(2*(1+Z²)^2)',
        total_remainder_factors=[5,10,44],
        absolute_error_by_order='factor_j*exp(3/5-yd)/(2*(1+Z²)^2)',
        raw_waiting_root_and_fixed_stage_limits_from_checked_original_operator=True,
        original_finite_axial_end_and_all_eleven_late_stages_in_error_budget=True,
        positive_late_source_terms_not_set_to_zero=True)


def exact_Z(value):
    value=source.exact_rational(value)
    if value.is_real is not True or value.is_finite is not True or not -1<=value<=1:
        raise ValueError('Original pressure requires real Z in[-1,1]')
    return value


class OriginalNormalizedPressurePointJets:
    def __init__(self,frame=None,*,cells=4096):
        if type(cells) is not int or cells<16:raise ValueError('At least 16 directed coefficient cells required')
        self.frame=frame or source.OriginalO2SourceParameterFrame();self.family=self.frame.family
        admitted=json.loads((HERE/source.RECEIPT).read_bytes())
        if not admitted['all_passed'] or not admitted[source.GATE] or admitted['source_family']!=self.family:
            raise ValueError('Checked original correlated source parameter frame required')
        self.hashes={**admitted['input_hashes'],source.RECEIPT:sha(source.RECEIPT)}
        for filename,digest in self.hashes.items():
            if sha(filename)!=digest:raise ValueError('Accepted pressure parameter source changed: '+filename)
        self.proof=late_source_error_theorem(self.frame);self.cells=cells
        self.ctx=c=self.frame.owner.profiles.ctx
        self.interval=iv=MPIntervalContext();iv.dps=c.dps+40
        with mp.workdps(iv.dps+40):
            # True point quadrature and true directed enclosure are separate.
            mass=c.quad(lambda v:c.exp(v/5-c.mpf('1.2')*self.frame.owner.profiles.J(v)),[0,c.mpf('.5'),1])
            self.alpha=c.mpf(5)/2+mass/2+c.exp(-c.mpf(2)/5)/2
            _,masses=transition_integrals(iv,Fraction(1),cells)
            self.alpha_enclosure=iv.mpf(5)/2+masses[1]/2+iv.exp(-iv.mpf(2)/5)/2
            lo,hi=endpoints(self.alpha_enclosure)
            if not lo<=self.alpha<=hi:raise ArithmeticError('Defining pressure quadrature lies outside its independent directed enclosure')
        for filename in ('lei_ren_part1_paper_interval_outer_slope_field.py',
            'lei_ren_part1_paper_interval_long_reshape_field.py',
            'lei_ren_part1_paper_schedule_endpoint_enclosures.py'):
            if self.hashes.get(filename)!=sha(filename):raise ValueError('Unbound original directed pressure coefficient source: '+filename)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def evaluate(self,Z):
        Z=exact_Z(Z);c=self.ctx;iv=self.interval
        with mp.workdps(iv.dps+40):
            z=c.mpf(int(Z.p))/int(Z.q);q=1+z*z
            iz=iv.mpf(int(Z.p))/int(Z.q);iq=1+iz*iz
            factors=(q**-2,-4*z/q**3,(-4+20*z*z)/q**4)
            ifactors=(iq**-2,-4*iz/iq**3,(-4+20*iz*iz)/iq**4)
            jets=[]
            yd=iv.exp(40)+11
            for order,(factor,ifactor,tail_factor) in enumerate(zip(factors,ifactors,(5,10,44))):
                approximation=-self.alpha*factor
                baseline=-self.alpha_enclosure*ifactor
                coefficient_error=endpoints(abs(baseline-iv.mpf(approximation)))[1]
                tail_log=endpoints(iv.mpf(3)/5-yd+iv.ln(tail_factor)-iv.ln(2)-2*iv.ln(iq))[1]
                if order==1 and Z==0:tail_log=None
                jets.append(dict(ordinary_Z_order=order,approximate_normalized_value=approximation,
                    coefficient_quadrature_and_arithmetic_absolute_error_upper=coefficient_error,
                    finite_end_and_all_late_source_error_log_upper=tail_log,
                    finite_end_and_all_late_source_error=dict(zero=tail_log is None,
                        log_upper=tail_log,representation='exact_zero' if tail_log is None else 'positive_exp_of_log_upper'),
                    exact_symmetry_zero=order==1 and Z==0))
        return dict(source_family=self.family,original_Z_exact=str(Z),normalized_pressure_ordinary_Z_jets=jets,
            original_source_point_approximation_with_explicit_error_budgets=True,
            native_pressure_integral_or_caps_selected_as_exact_values=False,
            full_physical_pressure_error_representation='Pstar²*(coefficient_error + exp(late_source_error_log_upper)); zero tail only for exact odd symmetry at Z=0',
            physical_Pstar_squared_error_amplification_not_discarded=True,
            full_native_phase_and_global_scalar_oracle_installed=False)


def run():
    began=time.monotonic();owner=OriginalNormalizedPressurePointJets()
    samples=[owner.evaluate(z) for z in ('-.7','0','.37','1')]
    record=dict(**{GATE:True},source_family=owner.family,original_late_source_error_theorem=owner.proof,
        actual_defining_quadrature_alpha=owner.alpha,
        directed_alpha_enclosure=dict(zip(('lower','upper'),endpoints(owner.alpha_enclosure))),
        directed_coefficient_cells=owner.cells,actual_original_pressure_point_queries=samples,
        current_original_normalized_pressure_point_jet_service_installed=True,
        original_late_integrals_retained_in_exact_frame_and_nonzero_error_budget=True,
        physical_Pstar_squared_error_amplification_not_discarded=True,
        original_p1_p2_scalar_point_values_installed=False,
        numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(source.inertial.profiles.loop.OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,
        scope='Approximate actual normalized original P0/Pstar² and ordinary Z jets through2, from defining early quadrature plus independently directed coefficient/roundoff error and source-bound all-late-stage error. Native physical Pstar² amplification remains explicit; full conditioned p1/p2, phase, global numerical oracle, controls and corrected field are open.')
    (HERE/NAME).write_text(json.dumps(source.inertial.profiles.loop.encoded(record),indent=2)+'\n',encoding='utf8')
    print('Original normalized pressure point jets and all-late-stage error budgets connected',flush=True)
    return record


if __name__=='__main__':run()
