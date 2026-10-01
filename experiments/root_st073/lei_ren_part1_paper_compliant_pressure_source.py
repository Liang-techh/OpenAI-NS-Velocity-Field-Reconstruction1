"""Separate Section 7 epsilon source and correlated pressure perturbation.

The legacy .01 source is not mutated. Late mass boxes are regenerated for
.001, and the difference of the TRUE correlated sources is bounded after
y_t, rather than by subtracting unrelated fourteen-atom interval boxes.
No old finite core or matching receipt is relabelled as a new solution.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_logarithmic_outer_parameters import LogarithmicOuterParameters
from lei_ren_part1_paper_logarithmic_pressure_datum import LogarithmicPressureDatum, BETA2, BETA0
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent


class CompliantOuterParameters(LogarithmicOuterParameters):
    def __init__(self, Md='40', precision=160):
        super().__init__(Md, precision)
        c = self.ctx
        with mp.workdps(precision+40):
            self.log_epsilon = c.ln(c.mpf('.001')) + self.log_delta
            self.epsilon = c.exp(self.log_epsilon)
            # Re-evaluate the continuous equation at the new epsilon.
            self.waiting = self.waiting_root_enclosure()
            self.edges[9] = self.waiting['root_interval']

    def report(self):
        result = super().report()
        result['actual_epsilon_relation'] = 'epsilon=.001*delta'
        result['Section7_c_epsilon_numeric_gate'] = True
        result['all_Section7_hypotheses_certified'] = False
        return result


class CompliantPressureDatum(LogarithmicPressureDatum):
    def __init__(self, Md='40', precision=160):
        self.parameters = CompliantOuterParameters(Md, precision)
        self.ctx = c = self.parameters.ctx
        p = self.parameters
        with mp.workdps(precision+40):
            refinement_name = 'lei_ren_part1_paper_candidate_pressure_mass_refinement.json'
            refined = json.loads((HERE/refinement_name).read_bytes())
            expected = '736bbadbde99bc2f3d098d279d61ef4cb64418368263a4aba7b275e7f8892de4'
            if refined['accepted_schedule_sha256'] != expected:
                raise ValueError('Unexpected epsilon-independent early atom')
            for name, digest in refined['input_hashes'].items():
                if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
                    raise ValueError('Early-atom dependency changed: '+name)
            names = (Path(__file__).name, refinement_name,
                     'lei_ren_part1_paper_logarithmic_pressure_datum.py',
                     'lei_ren_part1_paper_logarithmic_outer_parameters.py',
                     'lei_ren_part1_paper_candidate_pressure_function.py',
                     'lei_ren_part1_paper_outer.py')
            self.input_hashes = {n: hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
            self.input_hashes.update(refined['input_hashes'])
            self.tail_upper = c.exp(c.mpf('.6')-p.yd)/(2*(1-p.epsilon)**2)
            box = c.mpf([0, endpoints(self.tail_upper)[1]])
            self.stages = {n: dict(beta=2, mass=box) for n in BETA2}
            self.stages.update({n: dict(beta=0, mass=box) for n in BETA0})
            self.stages['reference_extension']['mass'] = c.mpf('2.5')
            self.stages['slope_transition_ref']['mass'] = read_interval(
                c, refined['stages']['slope_transition_ref']['refined_mass'])
            self.stages['axial_turnoff']['mass'] = c.exp(c.mpf('.6'))*(c.exp(-1)-c.exp(-p.yd))/2
            self.stages['z_flatten'] = dict(beta='variable [0,2]', mass=box)
            self.m2 = sum((self.stages[n]['mass'] for n in BETA2), c.mpf(0))
            self.m0 = sum((self.stages[n]['mass'] for n in BETA0), c.mpf(0))
            self.rho = c.mpf('.25')
            self.flatten_complex_upper = self.tail_upper/(1-self.rho**2)**2
            self.definition = dict(
                Md=str(Md), logPstar='exp(Md)+11', c_mu='.001', c_delta='.001', c_epsilon='.001',
                delta='min(1e-200,exp(-4logPstar-30))',
                waiting='unique positive root of the continuous raw preheat waiting equation',
                angular_profile='Section 6.1 reference-plus-outer ansatz, H replaced by 1',
                Tw='-60log(mu)', Ts='4log(2/delta)', Tf=100,
                cutoff_and_schedule_python_sha256=self.input_hashes['lei_ren_part1_paper_outer.py'])
            self.source_sha = hashlib.sha256(json.dumps(self.definition, sort_keys=True).encode()).hexdigest()
            enclosure = dict(source_sha256=self.source_sha, input_hashes=self.input_hashes,
                stages=encode(_pack(self.stages)), complex_flatten_upper=encode(self.flatten_complex_upper),
                normalized_m2=encode(self.m2), normalized_m0=encode(self.m0))
            self.datum_sha = hashlib.sha256(json.dumps(enclosure, sort_keys=True).encode()).hexdigest()


def pressure_perturbation(new, old=None):
    """Uniform complex pressure difference bound for the raw preheat sources.

    Before y_t both angular profiles coincide exactly. After y_f they are
    Z-independent (H is replaced by 1), so their pressure difference is an
    entire constant in Z. y_t contains 13/mu, making a coarse post-y_d bound
    needlessly weak after multiplication by Pstar^2. A proved, finite log
    cap avoids expanding exp(-13/mu), without replacing that tail by zero.
    """
    c = new.ctx
    old = old or LogarithmicPressureDatum(new.parameters.Md, c.dps)
    p = new.parameters
    with mp.workdps(c.dps+40):
        if old.definition['c_epsilon'] != '.01' or new.definition['c_epsilon'] != '.001':
            raise ValueError('This perturbation compares the pinned .01 and .001 sources')
        for name in ('Md', 'logPstar', 'c_mu', 'c_delta', 'delta', 'Tw', 'Ts', 'Tf',
                     'cutoff_and_schedule_python_sha256'):
            if old.definition[name] != new.definition[name]:
                raise ValueError('Unchanged upstream source requirement: '+name)
        if endpoints(p.mu)[1] > mp.mpf('.001') or endpoints(p.delta)[1] > mp.mpf('1e-200'):
            raise ValueError('Positive lengths and slope envelopes require small mu/delta')
        for datum in (new, old):
            if endpoints(datum.parameters.epsilon)[0] <= 0 or endpoints(datum.parameters.epsilon)[1] >= mp.mpf('.5'):
                raise ValueError('Positive collar amplitude required')
        cap_cut = 10*p.logPstar+1000
        pulse = 13/p.mu
        if endpoints(pulse-cap_cut)[0] <= 0:
            raise ValueError('Cannot certify selected finite log cap below y_t')
        # Joint tails, not the sum of individual stage envelopes.
        late_new = c.exp(c.mpf('.6')-cap_cut)/(2*(1-p.epsilon)**2)
        late_old = c.exp(c.mpf('.6')-cap_cut)/(2*(1-c.mpf(old.parameters.epsilon))**2)
        normalized = late_new+late_old
        physical = c.exp(2*p.logPstar)*normalized
        coarse = c.mpf(old.tail_upper)+new.tail_upper
        if endpoints(normalized)[0] <= 0 or endpoints(physical)[0] <= 0:
            raise ValueError('Pressure perturbation caps must remain strictly positive')
        return dict(
            old_implicit_source_sha256=old.source_sha, new_implicit_source_sha256=new.source_sha,
            new_datum_enclosure_sha256=new.datum_sha,
            unchanged_profile_region='all y<y_t; waiting starts at common y_t',
            unchanged_Z_dependence_region='all y<=y_f, including complete variable-beta flatten',
            changed_profile_region='waiting endpoint and preheat collar/exterior, entirely Z-independent',
            source_difference_formula='Delta P0=-(Pstar^2/2)*integral_y_t^infinity (u_new^2-u_old^2) dy',
            normalized_velocity_symbol='u=Utheta/Pstar',
            c_infinity_does_not_rescale_upstream=True,
            c_infinity_relation='c_inf=Utheta(R_t)*R_t^((1+delta)/2)/(1-epsilon), independent of waiting tau',
            radius_convention='R_t is the common START of waiting; R_tail=R_t*exp(tau); c_inf also equals Utheta(R_tail)*R_tail^((1+delta)/2)/(1-epsilon)',
            comparison_uses_same_reference_radius=True,
            pressure_integral_independent_of_reference_radius='dR/R=dy and the normalized preheat profile is a function of y,Z alone',
            inherited_parameter_constructor_dependency='legacy constructor initializes epsilon-independent data; epsilon, waiting and edges[9] are explicitly rebuilt',
            Delta_P0_is_exact_constant_in_Z=True,
            Delta_P0_Z_and_all_higher_Z_derivatives_exactly_zero=True,
            complex_difference_bound_uses_constant_extension=True,
            common_after_yt_squared_velocity_envelope='u^2<=exp(.6-y)/(1-epsilon)^2',
            finite_cap_use='Integrate the envelope from actual y_t, then use exp(-y_t)<=exp(-selected_cut); no pointwise profile claim below y_t',
            exact_formal_yt_lower_term='13/mu', selected_finite_log_cut=cap_cut,
            pulse_length_exceeds_selected_cut_margin=pulse-cap_cut,
            no_exp_minus_inverse_mu_expansion_used=True,
            new_joint_after_yt_tail_cap=late_new, old_joint_after_yt_tail_cap=late_old,
            normalized_pressure_difference_abs_upper=normalized,
            physical_pressure_difference_abs_upper=physical,
            physical_pressure_difference_log_upper=c.ln(physical),
            coarse_after_yd_difference_upper=coarse,
            coarse_after_yd_bound_not_used_for_core_transfer=True,
            exact_difference_value_evaluated=False,
            nonzero_tail_not_replaced_by_zero=True,
            new_core_contraction_certified=False, old_finite_core_transferred=False,
            all_Section7_hypotheses_certified=False, temporal_recursion=False)


def run():
    new = CompliantPressureDatum()
    with mp.workdps(210):
        result = dict(compliant_source=new.report(), pressure_perturbation=pressure_perturbation(new),
                      input_hashes=new.input_hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)), indent=2)+'\n', encoding='utf-8')
    print('Separate epsilon=.001 datum: 14 positive atoms; exact constant-Z pressure perturbation bounded with physical scale', flush=True)
    return result


if __name__ == '__main__':
    run()
