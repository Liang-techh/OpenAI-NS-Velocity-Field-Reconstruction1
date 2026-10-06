"""One current repair -> complete future C5 -> selected pulse functions.

Only the unchanged incoming/native pulse kernels are reused. Every future,
angular and selected coefficient cache belongs to the current exact branch.
The rebound pulse is a callable restricted view; installation in all native
physical owners, full exterior stress and temporal recursion remain separate.
"""
import ast
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp

from lei_ren_part1_paper_compliant_current_pressure_terminal_closure import (
    CurrentPressureTerminalClosure,HERE,PREFIX,sha,binding,function,pack)
from lei_ren_part1_paper_compliant_outer_buffer import decay_integral
from lei_ren_part1_paper_compliant_future_energy_high_jets import CompliantFutureEnergyHighJets
from lei_ren_part1_paper_compliant_axial_amplitude_selection import CompliantAxialAmplitude
from lei_ren_part1_paper_compliant_axial_high_jets import CompliantAxialHighJets
from lei_ren_part1_paper_compliant_fifth_axial_jets import CompliantFifthAxialJets
from lei_ren_part1_paper_compliant_pulse_mixed_C4 import CompliantPulseMixedC4
from lei_ren_part1_paper_compliant_pulse_high_jets import _SelectedSource
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

NAME=PREFIX+'current_selected_energy_source.json'
RECEIPT=PREFIX+'current_selected_energy_source_check.json'
GATES=('current_complete_future_energy_C4_C5_source_owned',
       'current_selected_ap_c1_c2_C5_source_owned',
       'current_selected_pulse_restricted_mixed4_callable')
OPEN=('current_exact_repair_installed_in_all_physical_charts',
      'heat_exterior_stress_identity_certified','global_completed_tensor_admissibility',
      'admissible_stress_lift_constructed','full_background_NS_validation',
      'physical_energy_integral_certified','independently_bounded_flat_remainder',
      'full_cartesian_vector_derivatives_certified','temporal_recursion')
VIEWS={'whole_Z':[-1,1],'axis':0,'fresh':'.709','endpoint':1}
PULSE_VIEWS={'main':('main','.709','1'), 'gap':('gap','.709','12'),
             'end_bump':('end','.709','-1'), 'terminal':('end',[-1,1],0)}


def original_algorithm_source_bindings():
    specs={
        ('compliant_future_energy_high_jets','future'):{
            'f':'self.base','high':'self.angular.coefficients(endpoints(Z))',
            'old':'f.future(Z)'},
        ('compliant_axial_amplitude_selection','select'):{
            'energy':'self.future.future(Z)','target':'future-incoming+self.base',
            'ap':'-halflinear+c.sqrt(halflinear**2-A0[0]/A2)',
            'scaled':'[uj+selected*vj for uj,vj in zip(u,self.v)]'},
        ('compliant_axial_high_jets','select'):{
            'b':'self.base',
            'future':"copy_jet(c, self.energy.future(endpoints(Z))['Section7_34_weighted_future_Taylor'])",
            'target':'future-incoming_energy+b.base','old':'b.select(Z)'},
        ('compliant_fifth_axial_jets','future'):{
            'f':'self.fourth.energy.base','high':'self.angular(Z)',
            'old':'self.fourth.energy.future(endpoints(Z))',
            'admitted':"append_fifth(c,old['complete_future_energy_Taylor'],total[5])"},
        ('compliant_fifth_axial_jets','select'):{
            'old':'self.fourth.select(Z)','b':'self.fourth.base',
            'future':"self.future(Z)['Section7_34_weighted_future_Taylor']",
            'target':'future-energy*b.mu+b.base',
            'ap':"selected_fifth(c,old['selected_ap_Taylor'],A2,A1,A0,old['positive_root_derivative_denominator'])"},
        ('compliant_pulse_radial_C4','data'):{'selected':'self.fifth.select(Z)'},
        ('compliant_axial_pulse_field','end'):{
            'future':"self.selection.future.future(Z)['complete_future_energy_Taylor']/2",
            'e':'future*c.exp(2*self.mu*s)+baseline-end_energy*(c.exp(2*self.mu*s)*self.E2cap)'},
        ('compliant_pulse_mixed_C4','_high_packet'):{
            'mixed':'transport_mixed(self,Z,point,Brows,u)'}}
    out={}
    for (stem,method),rows in specs.items():
        for target,expression in rows.items():
            binding(stem,method,target,expression)
            out[stem+'.'+method+':'+target]=True
    # This method first composes the full integral, then narrows its C1
    # prefix. Both assignments must be present in their original order.
    values=[node.value for node in ast.walk(function('compliant_future_energy_high_jets','future'))
            if isinstance(node,ast.Assign) and any(ast.unparse(t)=='total' for t in node.targets)]
    expected=['flatten+power+angular_change+postrel-heatdifference','IntervalTaylor(c, coeffs)']
    if [ast.dump(v) for v in values]!=[ast.dump(ast.parse(v,mode='eval').body) for v in expected]:
        raise ValueError('Complete C4 future composition/prefix narrowing changed')
    out['complete_future_integral_then_same_current_C1_narrowing']=True
    return out


@source_precision
def replay_complete_future(exact):
    """Use the original C1 integral recipe with current parameters and repair."""
    out=copy.copy(exact.future);out.repair=exact.repair
    out.angular=exact.repair.angular;out.heat=exact.repair.heat
    out.params=exact.repair.params;out.mu=exact.repair.mu;out.delta=exact.repair.delta
    c=out.ctx
    out.Lrel=-30*out.params.log_mu
    out.Ns=c.exp(-1-out.mu);out.Nq=out.Ns*c.exp(-2*out.params.Ts)
    out.Nt=out.Nq*c.exp(-1-out.delta/2)
    out.Ntail=out.Nt*c.exp(-out.delta*out.angular.waiting)
    out.steep_power=decay_integral(c,c.mpf(2),out.params.Ts)
    out.waiting_energy=decay_integral(c,out.delta,out.angular.waiting)
    out.tail_multiplier=out.Ntail*c.exp(-2*out.angular.waiting_logone)
    out.pulse_energy_attenuation=c.exp(-26)
    out.weighted_factor=out.mu*out.pulse_energy_attenuation/2
    # These kernels contain only unchanged mu/delta and the fixed sigma.
    # The checked defining-parameter bridge, not interval equality, permits reuse.
    out.kernels=dict(out.kernels)
    return out


class CurrentSelectedEnergySource:
    @source_precision
    def __init__(self,pressure=None,require_checked=True):
        self.pressure=pressure if pressure is not None else CurrentPressureTerminalClosure()
        if not self.pressure.acceptance_loaded:raise ValueError('Checked current pressure/repair required')
        self.exact=self.pressure.exact;self.family=self.exact.family
        self.source=self.exact.source;self.datum_sha=self.exact.datum_sha
        if not self.exact.parameter_bridge['passed']:raise ValueError('Defining parameter bridge required')
        old=self.exact.pulse.fifth;old4=old.fourth
        if (type(old),type(old4),type(old4.base),type(old4.energy),type(self.exact.pulse))!=(
                CompliantFifthAxialJets,CompliantAxialHighJets,CompliantAxialAmplitude,
                CompliantFutureEnergyHighJets,CompliantPulseMixedC4):
            raise ValueError('Unchanged typed original algorithms required')
        self.ctx=c=old.ctx;self.future=replay_complete_future(self.exact)
        self.amplitude=copy.copy(old4.base);self.amplitude.future=self.future
        self.amplitude.ctx=self.future.ctx
        self.amplitude.mu=c.mpf(endpoints(self.future.mu));self.amplitude.cache={}
        # Recompute the parameter-dependent caps/slopes. The independent
        # incoming factors and native Gram kernels retain their defining recipe.
        b=self.amplitude;b.base=(1-c.exp(-26))/4
        b.log_end_scale=b.pulse.logscale
        b.log_energy_cap=2*self.future.params.log_mu-1000;b.energy_cap=c.exp(b.log_energy_cap)
        b.energy_log_margins=[b.log_energy_cap-(self.future.params.log_mu+c.ln(Kj)+2*b.log_end_scale)
                              for Kj in b.pulse.basis['end_energy_weights']]
        if any(endpoints(v)[0]<=0 for v in b.energy_log_margins):raise ArithmeticError('End energy cap lost')
        b.nu=[c.mpf([0,endpoints(b.energy_cap)[1]]) for _ in range(2)]
        Pi=b.pulse.rows['scaled_full_rows'];b.v=b.linear_inverse(-b.mu*Pi[0],-(Pi[1]-Pi[0]))
        b.incoming_factor_log_definitions=[-13/(2*b.mu)+13*row-b.pulse.rows['common_logpref'] for row in (1,2)]
        b.incoming_factor_margins=[cap-exact for cap,exact in zip(
            b.incoming_factor_log_caps,b.incoming_factor_log_definitions)]
        if any(endpoints(v)[0]<=0 for v in b.incoming_factor_margins):raise ArithmeticError('Incoming factor cap lost')
        self.energy4=copy.copy(old4.energy);self.energy4.base=self.future
        self.energy4.ctx=self.future.ctx
        self.energy4.angular=self.exact.angular4;self.energy4.cache={}
        self.axial4=copy.copy(old4);self.axial4.energy=self.energy4
        self.axial4.ctx=self.future.ctx
        self.axial4.base=self.amplitude;self.axial4.cache={}
        self.fifth=copy.copy(old);self.fifth.fourth=self.axial4
        self.fifth.angular4=self.exact.angular4
        self.fifth.angular_cache={};self.fifth.energy_cache={};self.fifth.cache={}
        self.pulse=copy.copy(self.exact.pulse);self.pulse.fifth=self.fifth
        self.pulse.high=SimpleNamespace(ctx=c,base=self.amplitude,constants=self.axial4.constants,
            select=self.fifth.select,energy=SimpleNamespace(future=self.fifth.future))
        self.pulse.selection=_SelectedSource(self.pulse.high)
        self.pulse.pulse=self.amplitude.pulse;self.pulse.data_cache={}
        self.original_sources=(old,old4,old4.base,old4.energy,self.exact.pulse)
        self.original_cache_ids=tuple(id(cache) for cache in (
            old.angular_cache,old.energy_cache,old.cache,old4.cache,old4.energy.cache,self.exact.pulse.data_cache))
        self.source_bindings=original_algorithm_source_bindings()
        self.assert_graph()
        self.runtime={};self.pulse_runtime={};self.hashes=dict(self.pressure.hashes)
        for stem in ('future_swirl_energy','future_energy_high_jets','axial_amplitude_selection',
            'axial_high_jets','fifth_axial_jets','pulse_high_jets','pulse_radial_C4',
            'pulse_mixed_C4','axial_pulse_field','outer_pulse_map','flat_pulse_derivatives'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.amplitude.whole=self.amplitude.select([-1,1])
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current selected source scope/datum differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT)
            self.acceptance_loaded=True

    def assert_graph(self):
        old,old4,oldbase,oldenergy,oldpulse=self.original_sources
        graph=dict(
            all_current_future_repair_owners_share_one_branch=self.future.repair is self.exact.repair is self.energy4.angular.repair,
            current_future_parameters_owned=self.future.params is self.exact.repair.params,
            current_C1_C4_Taylor_contexts_match=self.amplitude.ctx is self.energy4.ctx is self.axial4.ctx is self.future.ctx is self.exact.angular4.ctx,
            current_C1_C4_C5_future_owned=self.energy4.base is self.axial4.base.future is self.fifth.fourth.energy.base is self.future,
            current_C4_C5_selected_amplitude_owned=self.fifth.fourth.base is self.amplitude,
            current_angular_C4_C5_owned=self.fifth.angular4 is self.energy4.angular is self.exact.angular4,
            native_incoming_shape_constants_preserved=self.axial4.constants is old4.constants is self.pulse.high.constants,
            native_incoming_and_fixed_kernel_source_preserved=self.amplitude.pulse is oldbase.pulse is self.pulse.pulse,
            current_pulse_selected_and_energy_callables_owned=self.pulse.fifth is self.fifth and
                self.pulse.selection.high is self.pulse.high and self.pulse.high.select.__self__ is self.fifth and
                self.pulse.high.energy.future.__self__ is self.fifth,
            defining_parameter_identity_consumed=self.exact.parameter_bridge['passed'],
            frozen_native_graph_unmutated=old.fourth is old4 and old4.base is oldbase and old4.energy is oldenergy
                and oldpulse.fifth is old,
            all_future_and_selection_caches_separate=not any(id(cache) in self.original_cache_ids for cache in (
                self.fifth.angular_cache,self.fifth.energy_cache,self.fifth.cache,self.axial4.cache,
                self.energy4.cache,self.pulse.data_cache)),
            inherited_original_C4_C5_and_pulse_algorithms=all(method.__func__ is original for method,original in (
                (self.energy4.future,CompliantFutureEnergyHighJets.future),
                (self.amplitude.select,CompliantAxialAmplitude.select),
                (self.axial4.select,CompliantAxialHighJets.select),
                (self.fifth.future,CompliantFifthAxialJets.future),(self.fifth.select,CompliantFifthAxialJets.select),
                (self.pulse._high_packet,CompliantPulseMixedC4._high_packet))))
        if not all(graph.values()):raise ValueError('Current selected source owner graph differs: '+str(graph))
        return graph

    @source_precision
    def evaluate(self,Z):
        self.assert_graph();key=tuple(endpoints(self.ctx.mpf(Z)))
        if key not in self.runtime:
            self.runtime[key]=dict(Z=self.ctx.mpf(Z),
                current_future_C4=self.energy4.future(Z),current_future_C5=self.fifth.future(Z),
                current_selected_C4=self.axial4.select(Z),current_selected_C5=self.fifth.select(Z),
                source_owner_graph=self.assert_graph(),
                exact_positive_heat_deficit_and_end_factors_retained=True,
                **dict.fromkeys(OPEN,False))
        return self.runtime[key]

    @source_precision
    def pulse_point(self,chart,Z,coordinate):
        self.assert_graph();key=(chart,tuple(endpoints(self.ctx.mpf(Z))),str(coordinate))
        if key not in self.pulse_runtime:
            if chart not in ('main','gap','gap_from_end','end'):raise ValueError('Unknown pulse chart')
            self.pulse_runtime[key]=getattr(self.pulse,chart)(Z,coordinate)
        return self.pulse_runtime[key]

    @source_precision
    def report(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,source_owner_graph=self.assert_graph(),
            original_algorithm_AST_bindings=self.source_bindings,
            defining_parameter_source_bridge=self.exact.parameter_bridge,
            current_selected_energy_views={name:self.evaluate(z) for name,z in VIEWS.items()},
            current_rebound_pulse_views={name:self.pulse_point(*args) for name,args in PULSE_VIEWS.items()},
            incoming_and_unchanged_kernel_scope='same native incoming functions, mu/delta definitions, fixed Gram kernels and formal row factors',
            pulse_installation_scope='restricted callable pulse owner; all downstream native physical owners not yet replaced',
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentSelectedEnergySource(require_checked=False)
    result=field.report()
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current complete future C5 and selected ap/c1/c2 rebuilt; restricted pulse mixed4 callable',flush=True)
    return result


if __name__=='__main__':run()
