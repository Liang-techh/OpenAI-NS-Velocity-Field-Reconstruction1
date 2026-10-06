"""Native affine pulse terminal -> one current angular/pressure repair branch.

The exact terminal is a source-defined function, not the old suppression
box. Shallow copies replay its waiting/radius/amplitude dependencies while
leaving the admitted physical graph unchanged. The old uniform contraction
is reused only on its certified ball; both derivative/future paths share
the new repair. Installing that branch in all physical charts is separate.
"""
import copy
import json
from pathlib import Path
from types import MethodType,SimpleNamespace

import mpmath as mp

from lei_ren_part1_paper_compliant_current_heat_pressure_stress import (
    CurrentHeatPressureStress,HERE,PREFIX,sha,function,binding,pack)
from lei_ren_part1_paper_compliant_outer_angular_repair import CompliantAngularRepair,magnitude,intersect
from lei_ren_part1_paper_compliant_outer_angular_candidate import SharedOuterAngularCandidate
from lei_ren_part1_paper_compliant_outer_buffer import decay_integral
from lei_ren_part1_paper_compliant_fifth_axial_jets import CompliantFifthAxialJets
from lei_ren_part1_paper_compliant_axial_pulse_field import CompliantAxialPulseField
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

NAME=PREFIX+'current_exact_repair_branch.json'
RECEIPT=PREFIX+'current_exact_repair_branch_check.json'
GATES=('current_native_Xv_exact_source_bound_to_repair',
       'current_exact_angular_pressure_unique_branch_recovered',
       'current_exact_repair_angular_C5_and_future_C1_available')
SCOPES=('current_exact_repair_installed_in_all_physical_charts',
    'current_heat_terminal_constants_eliminated','heat_exterior_stress_identity_certified',
    'global_completed_tensor_admissibility','admissible_stress_lift_constructed',
    'full_background_NS_validation','physical_energy_integral_certified',
    'independently_bounded_flat_remainder','full_point_physical_field_evaluation',
    'full_cartesian_vector_derivatives_certified','temporal_recursion')
VIEWS={'whole_Z':[-1,1],'axis':0,'fresh':'.431','endpoint':1}


def contained(inner,outer,label):
    il,ih=endpoints(inner);ol,oh=endpoints(outer)
    if not ol<=il<=ih<=oh:raise ArithmeticError('Uniform bound lost: '+label)
    return True


def replay_heat(heat,angular):
    """Recompute every Xv-dependent heat source; fixed Gamma atoms are reused."""
    out=copy.copy(heat);out.angular=angular;c=out.ctx;p=out.params
    out.tail_finite=p.yd+1+p.Tw+100-30*p.log_mu+2+p.Ts+angular.waiting
    out.logradius_terms=dict(selected_reference=out.logRref,pulse_term=13/p.mu,finite_offset=out.tail_finite)
    out.logS_terms={name:-value for name,value in out.logradius_terms.items()}
    out.logtail_relative=angular.waiting_field('0',1)['relative_log_profile'][0]
    out.logone=angular.waiting_logone
    if endpoints(-out.logRref-13/p.mu-out.tail_finite)[1]>=-1000:
        raise ArithmeticError('Replayed exact inverse-radius source lost its bound')
    return out


def replay_repair(original,angular,heat):
    out=copy.copy(original);out.angular=angular;out.heat=heat;c=out.ctx
    out.logtail_distance=2+out.params.Ts+angular.waiting
    out.logEtail_over_Erel=-c.mpf('1.5')*(2+out.params.Ts)+(out.rate+angular.restore_rate)/2-(c.mpf('.5')+out.delta/2)*angular.waiting
    out.log_theta_multiplier=(out.rate+angular.restore_rate)/2+angular.restore_rate*angular.waiting-angular.waiting_logone
    out.log_pressure_multiplier=2*out.logEtail_over_Erel-2*angular.waiting_logone+c.ln(out.delta/2)
    # Keep the already proved strong cap as a bound for the replayed source.
    out.inverse_radius_margin=heat.logRref+13/out.mu+heat.tail_finite+out.strong_logS_cap
    if endpoints(out.inverse_radius_margin)[0]<=0:raise ArithmeticError('Current inverse-radius bound failed')
    out.theta_heat_over_scale_cap=c.exp(out.log_theta_multiplier+out.strong_logS_cap-out.logscale)
    out.pressure_heat_over_scale_cap=c.exp(out.log_pressure_multiplier+out.strong_logS_cap-out.logscale)
    whole=out.defects([-1,1])
    out.rhs_scaled_sup=endpoints(c.mpf(magnitude(whole['r_scaled'][0]))+c.mpf(magnitude(whole['s_scaled'][0])))[1]
    out.rhs_derivative_scaled_sup=endpoints(c.mpf(magnitude(whole['r_scaled'][1]))+c.mpf(magnitude(whole['s_scaled'][1])))[1]
    out.linear_size=out.inverse_norm*max(endpoints(c.exp(out.rate)/out.weights['A'])[1],
        endpoints(c.exp(-3*out.prate)/out.weights['B'])[1])*out.rhs_scaled_sup
    # Weights, scale, old common ball, nonlinear bound and Lipschitz stay fixed.
    if endpoints(out.linear_size+out.nonlinear_size)[1]>=out.box or endpoints(out.lipschitz)[1]>=mp.mpf('.05'):
        raise ArithmeticError('Current exact-source map no longer contracts on the common ball')
    return out,whole


def replay_future(original,repair):
    out=copy.copy(original);out.repair=repair;out.angular=repair.angular;out.heat=repair.heat
    c=out.ctx
    out.Ntail=out.Nt*c.exp(-out.delta*out.angular.waiting)
    out.waiting_energy=decay_integral(c,out.delta,out.angular.waiting)
    out.tail_multiplier=out.Ntail*c.exp(-2*out.angular.waiting_logone)
    return out


class CurrentExactRepairBranch:
    @source_precision
    def __init__(self,companion=None,require_checked=True):
        self.companion=companion if companion is not None else CurrentHeatPressureStress()
        if not self.companion.acceptance_loaded:raise ValueError('Checked common current heat companion required')
        self.heat=self.companion.heat;self.pulse=self.heat.outer.pulse;self.flatten=self.heat.outer.flatten
        self.family=self.companion.family;self.source=self.companion.source;self.datum_sha=self.companion.datum_sha
        original=self.heat.future.repair;old4=self.pulse.fifth.angular4
        if type(original) is not CompliantAngularRepair or type(original.angular) is not SharedOuterAngularCandidate:
            raise ValueError('Original repair algorithms required')
        self.ctx=c=original.ctx;self.hashes=dict(self.companion.hashes)
        self.theorem=accepted(PREFIX+'outer_angular_repair_check.json',self.family,self.source,
            'actual_implicit_functional_angular_pressure_repair_independently_checked')
        _verify_hashes(self.theorem)
        self.hashes.update(self.theorem['input_hashes'])
        self.parameter_bridge=self.companion.source_owner.before.before.bindings['current_native_parameter_source_bridge']
        self.graph=dict(same_current_pulse_and_flatten=self.flatten.pulse is self.pulse,
            same_current_future=self.heat.future is self.pulse.fifth.fourth.energy.base,
            checked_current_parameter_function_bridge=self.parameter_bridge['passed'],
            both_original_repair_classes=type(old4.repair) is CompliantAngularRepair,
            both_pressure_sources_same=all((r.angular.initial.family,r.angular.initial.datum.source_sha,
                r.angular.initial.datum.datum_sha)==(self.family,self.source,self.datum_sha)
                for r in (original,old4.repair)),
            original_unique_branch_theorem=self.theorem['whole_Z_contraction_and_derivative_inverse_checked'])
        if not all(self.graph.values()):raise ValueError('Native and repair source graphs differ')
        angular=copy.copy(original.angular)
        angular.loguRp0=c.ln(c.mpf(endpoints(self.flatten.U)))
        angular.native_terminal_factor=self.pulse.factor
        # This is the native affine terminal recipe. The exponent stays formal;
        # its directed cap below is used only to enclose its nonzero exponential.
        angular.Xp=c.mpf(endpoints(self.pulse.Xp));equilibrium=1/angular.rate
        difference=angular.Xp-equilibrium;abs_difference=c.mpf(magnitude(difference))
        self.native_exponent=-13*angular.rate/angular.mu
        self.suppression_log=c.ln(1+abs_difference)+self.native_exponent
        if endpoints(self.suppression_log)[1]>=-1000:raise ArithmeticError('Native affine terminal suppression not proved')
        native_decay_bound=c.mpf(endpoints(angular.native_terminal_factor(-13*self.pulse.rate/self.pulse.mu)))
        analytic_decay_bound=c.mpf([0,endpoints(c.exp(-1000)/(1+abs_difference))[1]])
        decay_bound=intersect(c,native_decay_bound,analytic_decay_bound)
        angular.Xv=equilibrium+difference*decay_bound
        contained(angular.Xv,original.angular.Xv,'exact native Xv in old uniform suppression box')
        self.native_terminal_source=dict(definition='Xv*=1/(1-mu)+(Xp-1/(1-mu))*exp(-13*(1-mu)/mu)',
            Xp_definition='current pulse.buffer.power(0,1): Mtheta_over_sqrt2_R3half_Pstar / Utheta_over_Pstar',
            Xp_enclosure=angular.Xp,exponent=self.native_exponent,log_suppression_upper=self.suppression_log,
            Xv_enclosure=angular.Xv,exponential_enclosure=decay_bound,
            current_native_factor_callable_used=True,current_log_Utheta_Rp_origin=angular.loguRp0,
            exponential_cap_is_not_defining_value=True,old_suppression_box_is_not_defining_value=True)
        terminal=angular.before_waiting('0');Xt=terminal['X'][0];k=angular.restore_rate;eq=1/k
        if endpoints(Xt-eq)[0]<=0:raise ArithmeticError('Native prescribed waiting logarithm not positive')
        angular.waiting=(c.ln(Xt-eq)+angular.waiting_logone-angular.params.log_epsilon-c.ln(eq+angular.collarJ))/k
        contained(angular.waiting,original.angular.waiting,'recomputed native waiting root')
        selected_heat=replay_heat(original.heat,angular)
        self.repair,self.whole_rhs=replay_repair(original,angular,selected_heat)
        self.future=replay_future(self.heat.future,self.repair)
        self.angular4=copy.copy(old4);self.angular4.repair=self.repair;self.angular4.ctx=c;self.angular4.cache={}
        # Expose only the original angular fifth method. A copied complete
        # FifthAxialJets object would also expose stale future/select methods.
        self.angular5=SimpleNamespace(ctx=self.pulse.fifth.ctx,angular4=self.angular4,angular_cache={})
        self.angular5.angular=MethodType(CompliantFifthAxialJets.angular,self.angular5)
        self.graph.update(current_derivative_and_future_paths_share_one_repair=self.angular4.repair is self.future.repair,
            current_heat_angular_aliases_rebound=self.repair.angular is self.repair.heat.angular is self.future.angular,
            original_angular_fifth_method_only=self.angular5.angular.__func__ is CompliantFifthAxialJets.angular,
            same_native_terminal_factor_callable=angular.native_terminal_factor.__func__ is CompliantAxialPulseField.factor,
            current_heat_amplitude_origin_rebound=self.flatten.U is self.pulse.high.constants['U'],
            admitted_native_graph_not_mutated=self.heat.future.repair is original and self.pulse.fifth.angular4 is old4,
            both_old_default_weight_sources_reused=True)
        if not all(self.graph.values()):raise ValueError('Replayed current repair source graph differs')
        self.uniform=dict(old_common_coefficient_ball_radius=self.repair.box,
            current_scaled_rhs_sup=self.repair.rhs_scaled_sup,current_scaled_linear_map_size=self.repair.linear_size,
            nonlinear_map_size=self.repair.nonlinear_size,Lipschitz=self.repair.lipschitz,
            inverse_radius_log_margin=self.repair.inverse_radius_margin,
            native_terminal_containment_proved=True,recomputed_waiting_containment_proved=True,
            exact_inverse_radius_definition='S=exp(-(logRref+13/mu+tail_finite(waiting(Xv*))))',
            whole_Z_unique_smooth_branch_on_same_ball=True)
        for stem in ('current_pulse_flatten_source','current_power_angular_source','actual_Rp_source_join',
            'outer_angular_candidate','outer_angular_repair','outer_angular_repair_check','exact_heat_component',
            'angular_high_jets','fifth_axial_jets','future_swirl_energy','five_moment_repair','pulse_radial_C4'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.hashes[PREFIX+'outer_angular_repair_check.json']=sha(PREFIX+'outer_angular_repair_check.json')
        self.runtime={};self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in SCOPES):
                raise ValueError('Exact repair admission source or scope differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def evaluate(self,Z):
        z=self.ctx.mpf(Z);key=z._mpi_
        if key in self.runtime:return self.runtime[key]
        angular4=self.angular4.coefficients(z);angular5=self.angular5.angular(endpoints(z));future=self.future.future(z)
        out=dict(Z=z,current_native_terminal_source=self.native_terminal_source,
            current_angular_C4=angular4,current_angular_C5=angular5,current_complete_future_swirl_energy_C1=future,
            recomputed_waiting=self.repair.angular.waiting,
            current_log_inverse_radius_terms=self.repair.heat.logS_terms,
            current_log_heat_relative_amplitude=self.repair.heat.logtail_relative-self.repair.heat.logone,
            current_heat_amplitude_source_terms=dict(
                formal_Utheta_Rv0_origin=self.repair.angular.packet({})['formal_log_amplitude_origin'],
                logPstar=self.repair.params.logPstar,log_tail_relative=self.repair.heat.logtail_relative,
                minus_log_one_minus_epsilon=-self.repair.heat.logone,radial_power=(1+self.repair.delta)/2,
                logRtail_terms=self.repair.heat.logradius_terms),
            current_exact_unique_branch_equations=dict(angular='A*(exp(-3rate)*d1+exp(-rate)*d2)=r',
                pressure='B*(exp(3prate)*d1+exp(prate)*d2)+D/2*(exp(3prate)*d1^2+exp(prate)*d2^2)=sH'),
            same_current_native_source_graph=self.graph,uniform_unique_branch=self.uniform,
            old_cached_coefficients_not_relabelled=True,positive_exact_S_not_replaced_by_cap=True,
            current_prescribed_pressure_datum_unchanged=True)
        out.update({gate:self.acceptance_loaded for gate in GATES});out.update({gate:False for gate in SCOPES})
        self.runtime[key]=out;return out

    @source_precision
    def report(self):
        out=dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_native_terminal_source=self.native_terminal_source,
            same_current_native_source_graph=self.graph,uniform_unique_branch=self.uniform,
            current_exact_repair_views={name:self.evaluate(Z) for name,Z in VIEWS.items()},input_hashes=self.hashes)
        out.update({gate:False for gate in GATES+SCOPES});return out


@source_precision
def run(field=None):
    field=field if field is not None else CurrentExactRepairBranch(require_checked=False)
    result=field.report();(HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current native affine terminal replayed into one angular/pressure branch; physical installation and terminal constants pending',flush=True)
    return result


if __name__=='__main__':run()
