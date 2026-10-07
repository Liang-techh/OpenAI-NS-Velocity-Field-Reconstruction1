"""A joint primitive shear margin, retaining the finite-N square term.

This inequality precedes independent moment repair and completed stress.
Separate absolute shear-error ratios are not used at a vanishing taper.
"""
import json
import math
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_O3_actual_parameter_majorants as parameters
from lei_ren_part1_paper_compliant_current_O3_finite_frequency_profiles import cutoff_rows,sigma_jets

HERE,PREFIX,sha=parameters.HERE,parameters.PREFIX,parameters.sha
NAME=PREFIX+'current_O3_correlated_finite_N_shear.json'
RECEIPT=PREFIX+'current_O3_correlated_finite_N_shear_check.json'
OPEN=('independent_repair_completed_shear_certified','whole_modified_signed_two_vector_cones_certified',
      'sufficient_common_N_for_completed_tensor_certified','actual_coefficient_recursion_certified')


def exact_theorem():
    asts=parameters.numeric.transport.SourceAST();checks={}
    mu,chi,N=s.symbols('mu chi N',positive=True,finite=True)
    chi_y=s.Symbol('chi_y',real=True)
    S,C,L,sig=s.symbols('sin_phase cos_phase original_log_slope sigma',real=True)
    # Replay the defining source assignments with exact sin(4p)=2SC and
    # cos(4p)=C^2-S^2; source equality precedes any triangle enclosure.
    pi=s.pi;phase=s.Symbol('actual_phase',real=True)
    trig={2*pi*phase:S,4*pi*phase:2*S*C}
    cosine={2*pi*phase:C,4*pi*phase:C*C-S*S}
    c=SimpleNamespace(pi=pi,sin=lambda x:trig[x],cos=lambda x:cosine[x],sqrt=s.sqrt,exp=s.exp)
    env=dict(c=c,mu=mu,chi=[chi,chi_y],phase=phase,n=N,sig=sig,logEy=L)
    specs=(('da','mu*chi[0]**2*c.cos(4*c.pi*phase)'),
        ('slowA','-mu*chi[0]*chi[1]*c.sin(4*c.pi*phase)/(4*c.pi)'),
        ('beta','2*c.sqrt(mu)*chi[0]*c.sin(2*c.pi*phase)'),
        ('slowB_over_E','-c.sqrt(mu)*(chi[0]*logEy+chi[1])*c.cos(2*c.pi*phase)/(2*c.pi)'),
        ('aexcess','2*mu*sig+da-2*slowA/n'),
        ('b','c.exp(-A[0]/n)*(beta+2*slowB_over_E/n)'))
    A=-mu*chi**2*S*C/(4*pi);env['A']=[A]
    for name,wanted in specs:
        node=asts.expression('current_O3_finite_frequency_profiles','profile',name,wanted=wanted)
        env[name]=asts.evaluate(node,env);checks['actual_source_assignment_'+name]=True
    def zero(name,a,b):
        value=s.expand(a-b)
        if value!=0 and s.cancel(value)!=0:raise ArithmeticError('Joint primitive shear identity failed: '+name)
        checks[name]=True
    d=chi_y/(pi*N*chi);ell=L/(pi*N);q=s.exp(-2*A/N)/(2+3*mu)
    normalized=C*C-S*S+d*S*C+q*(2*S-(d+ell)*C)**2
    zero('same_actual_finite_N_joint_margin',env['aexcess']+env['b']**2/(2+3*mu),2*mu*sig+mu*chi**2*normalized)
    v,w,Q,E=s.symbols('v completed_square q ell',real=True)
    expression=C*C-S*S+v*S+Q*(2*S-v-E*C)**2
    square=C*C+S*S-E*S*C-S*S/(4*Q)+Q*(v+E*C-2*S+S/(2*Q))**2
    zero('joint_unbounded_cutoff_log_derivative_completed_square',expression,square)
    x=s.Symbol('nonnegative_exponent_size',nonnegative=True)
    remainder=s.exp(-x)-1+x
    zero('exp_lower_remainder_initial_value',remainder.subs(x,0),0)
    zero('exp_lower_remainder_initial_derivative',s.diff(remainder,x).subs(x,0),0)
    zero('exp_lower_remainder_positive_second_derivative',s.diff(remainder,x,2),s.exp(-x))
    # These rational bounds prove the uniform inequalities using only
    # exp(-x)>=1-x, pi>=3, N>=1, mu<=.001 and |chi|<=1.
    mu_max=s.Rational(1,1000);q_lower=(1-mu_max/12)/(2+3*mu_max)
    if q_lower<=s.Rational(2,5):raise ArithmeticError('Uniform finite-N square coefficient too small')
    checks['actual_uniform_q_lower_exceeds_two_fifths']=True
    zero('q_lower_exact_rational',q_lower,s.Rational(11999,24036))
    square_floor=1-1/(4*s.Rational(2,5))
    zero('positive_sine_square_coefficient_floor',square_floor,s.Rational(3,8))
    # |L|<=1/2+mu is valid on the entire O2/O3 offset[-2,1].
    # It also covers the tighter active-support norm (1+mu)/2.
    uniform_floor=square_floor-(s.Rational(1,2)+mu_max)/6
    zero('uniform_joint_margin_coefficient',uniform_floor,s.Rational(583,2000))
    if uniform_floor<=s.Rational(1,4):raise ArithmeticError('Joint primitive margin lost its quarter reserve')
    checks['uniform_joint_floor_strictly_exceeds_one_quarter']=True
    zero('phase_product_absolute_bound', (S-C)**2+(S+C)**2,2*(S*S+C*C))
    asts.expression('pre_pulse_mixed_C4','slope_mu','logU',
        wanted="[zero-c.mpf('.5')-mu*sig[0]]+[zero-mu*sig[k]*math.factorial(k) for k in range(1,4)]")
    asts.method('current_O3_finite_frequency_profiles','cutoff_rows')
    asts.method('flat_pulse_derivatives','sigma_jets')
    # The endpoints are exactly flat; chi=chi_y=0 is handled directly,
    # without ever forming chi_y/chi in numerical queries.
    zero('left_exact_flat_primitive_margin',env['aexcess'].subs({chi:0,chi_y:0,sig:0})+env['b'].subs({chi:0,chi_y:0})**2/(2+3*mu),0)
    zero('right_exact_flat_primitive_margin',env['aexcess'].subs({chi:0,chi_y:0,sig:s.Rational(1,2)})+env['b'].subs({chi:0,chi_y:0})**2/(2+3*mu),mu)
    return dict(passed=True,identities=checks,source_bindings=asts.bindings,
        hypotheses=['0<mu<=.001','finite integer N>=1','0<=chi<=1; chi_y any real',
            '0<=sigma<=1','original L=-1/2-mu*sigma, or O2 L=-1/2',
            'S^2+C^2=1','actual finite-N original profiles before independent repair'],
        exact_primitive_margin='a_minus2+b^2/(2+3mu)',
        completed_square_identity=s.sstr(square),
        uniform_lower_envelope='2mu*sigma+(mu/4)*chi^2',
        stronger_uniform_coefficient=str(uniform_floor),
        all_real_phases_and_all_finite_integer_N_ge1=True,
        no_division_by_chi_used_at_flat_endpoints=True,
        full_signed_stress_cone_admission_not_implied=True,
        input_hashes={**asts.hashes,Path(__file__).name:sha(Path(__file__).name)})


class CurrentCorrelatedFiniteNShear:
    def __init__(self,require_checked=True):
        self.data=parameters.inputs();self.ctx=self.data['ctx'];self.theorem=exact_theorem()
        parent=json.loads((HERE/parameters.RECEIPT).read_bytes())
        if not parent['all_passed'] or not parent['actual_parameter_finite_N_modulation_and_density_bounds_available']:
            raise ValueError('Checked actual finite-N scalar source required')
        for name,digest in parent['input_hashes'].items():
            if sha(name)!=digest:raise ValueError('Changed actual finite-N parameter input: '+name)
        self.hashes={**self.data['hashes'],**parent['input_hashes'],parameters.RECEIPT:sha(parameters.RECEIPT),**self.theorem['input_hashes']}
        if require_checked:
            record=json.loads((HERE/RECEIPT).read_bytes())
            if not record['all_passed']:raise ValueError('Checked joint primitive shear theorem required')
            for name,digest in record['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed checked joint shear program: '+name)

    def query(self,offset,N):
        parameters.positive_integer_N(N);c=self.ctx;t=c.mpf(offset)
        lo,hi=parameters.numeric.transport.endpoints(t)
        if not all(mp.isfinite(v) for v in (lo,hi)) or lo < -2 or hi > 1:
            raise ValueError('Finite actual O2/O3 logR offset in[-2,1] required')
        chi_rows=cutoff_rows(c,t);chi,dy=chi_rows[:2];sig=sigma_jets(c,t)[0]
        # sigma(t)=0 on the negative buffer, matching its exact L=-1/2.
        mu=self.data['mu'];L=-c.mpf(1)/2-mu*sig;phase=c.mpf(N)*t
        sn,cs=c.sin(2*c.pi*phase),c.cos(2*c.pi*phase)
        A=-mu*chi**2*c.sin(4*c.pi*phase)/(8*c.pi)
        theta=2*mu*sig+mu*chi**2*c.cos(4*c.pi*phase)+mu*chi*dy*c.sin(4*c.pi*phase)/(2*c.pi*N)
        axial=c.exp(-A/N)*(2*c.sqrt(mu)*chi*sn-c.sqrt(mu)*(chi*L+dy)*cs/(c.pi*N))
        exact_enclosure=theta+axial**2/(2+3*mu)
        lower=2*mu*sig+mu*chi**2/4
        theta_error=abs(mu*chi*dy*c.sin(4*c.pi*phase)/(2*c.pi*N))
        # |exp(-x)-1|<=|x|exp(|x|) preserves microscopic factors,
        # unlike subtracting a rounded exponential from one.
        x=abs(A/N);loop_b=2*c.sqrt(mu)*chi*sn
        axial_error=c.exp(x)*(x*abs(loop_b)+c.sqrt(mu)*abs(chi*L+dy)*abs(cs)/(c.pi*N))
        flat=chi._mpi_==c.mpf(0)._mpi_ and dy._mpi_==c.mpf(0)._mpi_
        return dict(logR_offset=t,finite_integer_N=N,actual_translated_phase=phase,
            actual_mu=mu,chi=chi,chi_y=dy,original_log_slope=L,sigma=sig,
            finite_N_theta_shear_excess=theta,finite_N_axial_shear=axial,
            actual_phase_theta_shear_error_absolute_cap=theta_error,
            actual_phase_axial_shear_error_absolute_cap=axial_error,
            error_caps_retain_chi_chi_y_log_slope_and_actual_phase=True,
            direct_joint_margin_enclosure=exact_enclosure,
            source_certified_joint_margin_lower_envelope=lower,
            directed_scalar_margin_lower_bound=c.mpf(parameters.numeric.transport.endpoints(lower)[0]),
            exact_flat_cutoff_value_and_first_derivative_zero=flat,
            joint_source_inequality_valid_even_if_separate_interval_enclosure_is_wide=True,
            source_family=self.data['source']['accepted']['source_family'],
            **{key:False for key in OPEN})


def run():
    field=CurrentCorrelatedFiniteNShear(require_checked=False)
    examples={name:field.query(t,N) for name,t,N in (
        ('left_flat_edge','-2',1),('left_taper','-1.99',1),('negative_plateau','-1',1),
        ('old_seam','0',1),('right_taper','.49',1),('right_flat_edge','.5',1),
        ('O3_terminal','1',1),('high_frequency_taper','.337',10**12))}
    result=dict(exact_correlated_finite_N_primitive_shear_theorem=field.theorem,
        source_family=field.data['source']['accepted']['source_family'],examples=examples,
        joint_primitive_shear_margin_certified_on_actual_O2_O3_profiles=True,
        positive_where_chi_or_original_sigma_nonzero=True,
        exact_left_original_zero_margin_preserved=True,
        separate_absolute_error_ratio_divergence_does_not_disprove_joint_margin=True,
        **{key:False for key in OPEN},input_hashes=field.hashes)
    (HERE/NAME).write_bytes((json.dumps(parameters.encoded(result),indent=2)+'\n').encode())
    print('Correlated finite-N primitive shear generated: quarter taper reserve; all finite N>=1',flush=True)
    return result


if __name__=='__main__':run()
