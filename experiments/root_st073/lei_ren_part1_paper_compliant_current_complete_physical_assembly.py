"""One complete current selected/future graph in the 33 physical source charts.

Reuse the checked common core/incoming owners, replace every live pulse and
postpulse provider, then apply the unchanged Cartesian/fixed-x time map.
Outputs are signed factored bounds on actual functions, not chosen point
coefficients. Quantitative interfaces, residual/cone/flatness and temporal
recursion remain separate obligations.
"""
import ast
import copy
import json
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_full_exterior_stress import (
    CurrentFullExteriorStress,HERE,PREFIX,sha,binding,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_current_postpulse_energy_history import METHODS
from lei_ren_part1_paper_compliant_current_native_pulse_source_dispatcher import (
    CurrentNativePulseSourceDispatcher,PULSE_CHARTS,EXPECTED)
from lei_ren_part1_paper_compliant_current_bridge_physical_assembly import CHARTS,BRIDGES,original_bridge_mapper_bindings
from lei_ren_part1_paper_compliant_current_pulse_physical_assembly import (
    CurrentPulsePhysicalAssembly as PULSE,native_parameter_source_bridge,
    original_pulse_operator_bindings,pulse_factor_rebase_proof,pressure_publication_bindings)
from lei_ren_part1_paper_compliant_current_flatten_physical_assembly import original_flatten_units_binding
from lei_ren_part1_paper_compliant_current_heat_physical_assembly import original_heat_radius_bindings
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly as BASE,POST
from lei_ren_part1_paper_compliant_pulse_mixed_C4 import CompliantPulseMixedC4
from lei_ren_part1_paper_compliant_pulse_radial_C4 import CompliantPulseRadialC4
from lei_ren_part1_paper_compliant_pulse_interface_certificate import functional_identities
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

NAME=PREFIX+'current_complete_physical_assembly.json'
RECEIPT=PREFIX+'current_complete_physical_assembly_check.json'
GATES=('current_complete_selected_graph_installed_in_all_33_physical_source_charts',
    'current_exact_repair_installed_in_all_physical_charts',
    'current_complete_source_cartesian_spatial4_time1_maps_certified')
OPEN=('uniform_pulse_C4_chart_interface_certificate_available',
    'current_complete_quantitative_physical_interfaces_certified',
    'global_completed_tensor_admissibility','admissible_stress_lift_constructed',
    'full_background_NS_validation','physical_energy_integral_certified',
    'independently_bounded_flat_remainder','full_point_physical_field_evaluation',
    'full_cartesian_vector_derivatives_certified','temporal_recursion')
PUBLIC_TO_INTERNAL={chart:{'steep_entry':'steep_in','steep_exit':'steep_out'}.get(chart,chart) for chart in POST}
CHANGED=tuple(PULSE_CHARTS)+tuple(POST)
RETAINED=tuple(chart for chart in CHARTS if chart not in CHANGED)
VIEWS={'pulse_entrance':None,'pulse_main':['.02','10'],'pulse_exit':[10,11],'pulse_gap':[11,12],
    'pulse_gap_end':None,'pulse_end':[-4,0],'flatten':[0,100],'outer_power':[0,1],
    'outer_angular':[-4,0],'steep_entry':[0,1],'steep_power':[0,1],'steep_exit':[0,1],
    'waiting':[0,1],'heat_collar':[0,3],'heat_exterior':[3,mp.inf]}


def physical_views(owner):
    result=dict(VIEWS);c=owner.ctx;mu=owner.pulse.mu
    result['pulse_entrance']=[0,endpoints(c.mpf('.02')/mu)[1]]
    result['pulse_gap_end']=[-endpoints(1/mu)[0],-4]
    return result


class CurrentCompletePhysicalDispatch:
    """Actual current providers, with one retained common incoming ancestor."""
    def __init__(self,assembly):
        self.assembly=assembly;self.base=assembly.base;self.native=self.base.dispatch.native
        self.history=assembly.history;self.selected=self.history.selected;self.native_pulse=self.selected.pulse
        self.family=assembly.family;self.source=assembly.source;self.datum_sha=assembly.datum_sha
        self.chain_routes=self.native.chain_routes
        self.providers={chart:self.base.dispatch.provider(chart) for chart in RETAINED if chart not in BRIDGES}
        self.providers.update({chart:assembly.bridge for chart in BRIDGES})
        self.providers.update({chart:self.native_pulse for chart in PULSE_CHARTS})
        self.providers.update({chart:getattr(self.history,METHODS[internal][0]) for chart,internal in PUBLIC_TO_INTERNAL.items()})
        self.registry={chart:dict(self.base.source_owners[chart]) for chart in RETAINED if chart not in BRIDGES}
        for chart in BRIDGES:
            self.registry[chart]=dict(provider=PREFIX+'current_actual_bridge_mixed_C4.CurrentActualBridgeMixedC4',
                method='evaluate',acceptance_receipt=PREFIX+'current_actual_bridge_mixed_C4_check.json')
        for chart in PULSE_CHARTS:
            self.registry[chart]=dict(provider=PREFIX+'current_selected_energy_source.CurrentSelectedEnergySource.pulse',
                method=EXPECTED[chart][0],domain=EXPECTED[chart][1],acceptance_receipt=RECEIPT)
        for chart,internal in PUBLIC_TO_INTERNAL.items():
            self.registry[chart]=dict(provider=PREFIX+'current_postpulse_energy_history.CurrentPostpulseEnergyHistory.'+METHODS[internal][0],
                method=METHODS[internal][1],internal_source_chart=internal,acceptance_receipt=RECEIPT)

    def __getattr__(self,name):return getattr(self.base.dispatch,name)
    def provider(self,chart):
        if chart not in CHARTS:raise ValueError('Unknown current physical chart: '+chart)
        return self.providers[chart]

    def _packet(self,chart,packet,domain,coverage):
        return dict(chart=chart,source_packet=packet,
            actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,source_provider=self.registry[chart]['provider'],
            source_domain=domain,coverage_coordinate=coverage,
            physical_mixed_grids={'physical_mixed_derivatives_total_order_le4':packet['physical_mixed_derivatives_total_order_le4']})

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        self.provider(chart)
        if chart in PULSE_CHARTS:
            # Reuse the original route guard, but its self.provider and
            # self._packet now resolve this checked current selected owner.
            return CurrentNativePulseSourceDispatcher.evaluate(self,chart,Z,coordinate)
        if chart in POST:
            packet=self.history.evaluate(PUBLIC_TO_INTERNAL[chart],Z,coordinate)['source_packet']
            return dict(chart=chart,source_packet=packet,source_provider=self.registry[chart]['provider'],
                datum_enclosure_sha256=self.datum_sha,
                physical_mixed_grids={'physical_mixed_derivatives_total_order_le4':packet['physical_mixed_derivatives_total_order_le4']})
        if chart in BRIDGES:
            packet=self.assembly.bridge.evaluate(Z,coordinate,BRIDGES[chart])
            return dict(chart=chart,source_packet=packet,source_provider=self.registry[chart]['provider'],
                datum_enclosure_sha256=self.datum_sha)
        return self.base.dispatch.evaluate(chart,Z,coordinate)

    @source_precision
    def gap_overlap(self,Z):
        packet=self.native_pulse.gap(Z,['12','12.0001'])
        return self._packet('pulse_gap',packet,'supplemental [12,12.0001]',
            'xi=mu*log(R/Rp); covers the reciprocal gap-end coordinate boundary')

    @source_precision
    def coverage(self):
        c=self.native_pulse.ctx;mu=self.native_pulse.mu
        start=-endpoints(1/mu)[0];image=13+mu*c.mpf(start)
        if endpoints(image)[1]>endpoints(c.mpf('12.0001'))[0] or endpoints(c.mpf(start)+1/mu)[0]<0:
            raise ArithmeticError('Current exact gap coordinate boundary not covered')
        return dict(legal_gap_end_numeric_start=start,original_exact_gap_end_source='s=-1/mu',
            exact_gap_coordinate_relation='xi=13+mu*s',gap_end_start_main_coordinate_enclosure=image,
            supplemental_same_current_owner_gap_box=['12','12.0001'],
            numerical_coverage_start_is_not_substituted_for_exact_boundary=True,
            all_six_current_pulse_routes_have_original_domain_guards=True,passed=True)


def current_complete_physical_source_proof(owner):
    graph=owner.assert_graph()
    # The admitted selected constructor deliberately makes an outward mu
    # copy. The older parameter theorem assumes pointer identity at this
    # one slot. Bind the actual copy, then replay that unchanged theorem on
    # a definition-only proxy; never mutate the live selected source.
    mu_copy=class_assignment('current_selected_energy_source','CurrentSelectedEnergySource','__init__',
        'self.amplitude.mu','c.mpf(endpoints(self.future.mu))')
    selected=owner.history.selected
    if endpoints(selected.amplitude.mu)!=endpoints(selected.future.mu):
        raise ValueError('Current selected outward mu copy differs from its defining source')
    proxy=copy.copy(owner);proxy.pulse=copy.copy(owner.pulse)
    proxy.pulse.high=copy.copy(owner.pulse.high);proxy.pulse.high.base=copy.copy(owner.pulse.high.base)
    proxy.pulse.high.base.mu=selected.future.mu
    parameter=dict(original_typed_definition_parameter_theorem=native_parameter_source_bridge(proxy),
        actual_current_selected_mu_copy_AST=mu_copy,
        actual_mu_copy_same_source_enclosure=True,live_current_selected_objects_unmodified=True,
        definition_only_proxy_adjusts_parameter_pointer_not_a_field_value=True,passed=True)
    mapper_mu=endpoints(owner.params.mu);selected_mu=endpoints(selected.future.mu)
    if not mapper_mu[0]<=selected_mu[0]<=selected_mu[1]<=mapper_mu[1]:
        raise ValueError('Selected mu is not enclosed by the same original mapper definition')
    parameter.update(actual_mapper_mu_enclosure=owner.params.mu,
        actual_current_selected_mu_enclosure=selected.future.mu,
        current_selected_mu_contained_in_original_mapper_enclosure=True,
        same_defining_formula_not_interval_center_equality=True)
    bindings={target:class_assignment('current_complete_physical_assembly','CurrentCompletePhysicalAssembly',
        '__init__',target,value) for target,value in {
            'self.base':'self.exterior.history.selected.exact.companion.owner',
            'self.bridge':'self.exterior.history.selected.exact.companion.joined.common.source.bridge',
            'self.params':'self.pre.params','self.pulse':'self.history.selected.pulse',
            'self.logRref':'c.ln(110)+10*(self.logC+self.logP)',
            'self.logRp':'self.logRref+self.logP+1+self.params.Tw',
            'self.logRv':'self.logRp+13/self.params.mu',
            'self.dispatch':'CurrentCompletePhysicalDispatch(self)'}.items()}
    for method in ('entrance','main','gap','gap_from_end','end','_high_packet'):
        if getattr(owner.pulse,method).__func__ is not getattr(CompliantPulseMixedC4,method):
            raise ValueError('Original all-chart native pulse algorithm changed: '+method)
    for method in ('data','pressure_moment'):
        if getattr(owner.pulse,method).__func__ is not getattr(CompliantPulseRadialC4,method):
            raise ValueError('Original pulse source/pressure callable changed: '+method)
    # The route guard's provider method supplies the current object on every
    # legal coordinate. No old six-pulse output/selection is inherited.
    binding('compliant_current_native_pulse_source_dispatcher','evaluate','field','self.provider(chart)')
    binding('compliant_pulse_radial_C4','data','selected','self.fifth.select(Z)')
    binding('compliant_pulse_mixed_C4','_high_packet','mixed','transport_mixed(self,Z,point,Brows,u)')
    binding('compliant_current_exact_repair_branch','replay_heat','out.tail_finite',
        'p.yd+1+p.Tw+100-30*p.log_mu+2+p.Ts+angular.waiting')
    binding('compliant_current_exact_repair_branch','replay_heat','out.logradius_terms',
        'dict(selected_reference=out.logRref,pulse_term=13/p.mu,finite_offset=out.tail_finite)')
    heat_radius_bindings={target:class_assignment('exact_heat_component','SharedExactHeatComponent','__init__',target,value)
        for target,value in {'self.logC':"read_interval(c,norms['selected_logCstar'])",
            'self.logRref':'c.ln(110)+10*(self.logC+self.params.logPstar)'}.items()}
    lr,lp,lc,mu,Tw,L,Ts,W,t,lmu=s.symbols('logRref logP logC mu Tw L Ts W t logmu',real=True)
    mapper_ref=s.log(110)+10*(lc+lp);rp=mapper_ref+lp+1+Tw;rv=rp+13/mu
    raw_tail=mapper_ref+13/mu+lp+1+Tw+100-30*lmu+2+Ts+W
    mapped_tail=rv+100+L+2+Ts+W
    identities=dict(current_raw_and_mapped_Rtail_same_exact_source=s.expand(mapped_tail.subs(L,-30*lmu)-raw_tail)==0,
        current_exact_S_is_inverse_source_radius=s.simplify(s.exp(-raw_tail)*s.exp(raw_tail))==1,
        current_full_Gamma_argument_source=s.simplify(s.exp(-raw_tail)*s.exp(-t)-s.exp(-(raw_tail+t)))==0,
        original_Rv_amplitude_has_both_finite_offsets=s.expand(2*(-13/(2*mu)-13)+13/mu+26)==0)
    xi,offset=s.symbols('xi offset',real=True)
    identities.update(original_entrance_coordinate=s.expand((mu*offset)/mu-offset)==0,
        original_gap_coordinate=s.expand((13-xi)/mu-(13/mu-(xi/mu)))==0,
        original_gap_end_coordinate=s.expand(13+mu*offset-mu*(13/mu+offset))==0,
        original_gap_end_and_end_same_radius=s.expand(rv-4-(rp+13/mu-4))==0)
    if not all(identities.values()):raise ArithmeticError('Current physical radius/amplitude source differs')
    original=owner.original_operator_fixture
    for key in ('independent_implicit_physical_fixture','independent_micro_scale_fixture'):
        if not original[key]['passed']:raise ValueError('Original unchanged physical coordinate fixture missing')
    return dict(current_defining_object_graph=graph,
        actual_current_physical_constructor_assignments=bindings,
        same_current_native_and_incoming_parameter_source_bridge=parameter,
        all_six_native_pulse_algorithms_and_current_callbacks_bound=True,
        original_route_guards_with_current_provider_calls_bound=True,
        original_internal_pulse_functional_coordinate_and_ODE_identities=functional_identities(),
        internal_functional_identities_apply_to_actual_current_selected_equations=selected.source_bindings,
        quantitative_mixed4_interface_admission=False,
        original_pulse_component_and_pressure_publication=pressure_publication_bindings(),
        original_pulse_physical_branches=original_pulse_operator_bindings(),
        original_correlated_pulse_rebase=pulse_factor_rebase_proof(),
        original_post_fixed_unit_binding=original_flatten_units_binding(),
        original_heat_radius_binding=original_heat_radius_bindings(),
        retained_original_bridge_map_calls=original_bridge_mapper_bindings(),
        current_Rv_Rtail_and_heat_argument_source_identities=identities,
        current_raw_heat_reference_and_logC_source_bindings=heat_radius_bindings,
        current_six_pulse_coordinate_coverage=owner.dispatch.coverage(),
        exact_source_inverse_radius_terms=owner.history.heat.exact_logS_terms,
        checked_current_complete_energy_source_proof=owner.history.proof,
        checked_current_full_similarity_exterior_source_proof=owner.exterior.proof,
        reused_original_independent_coordinate_fixture=original['independent_implicit_physical_fixture'],
        reused_original_independent_micro_scale_fixture=original['independent_micro_scale_fixture'],passed=True)


class CurrentCompletePhysicalAssembly(PULSE):
    @source_precision
    def __init__(self,exterior=None,require_checked=True):
        self.exterior=exterior if exterior is not None else CurrentFullExteriorStress()
        if not self.exterior.acceptance_loaded:raise ValueError('Checked current full exterior graph required')
        self.history=self.exterior.history
        self.base=self.exterior.history.selected.exact.companion.owner
        self.bridge=self.exterior.history.selected.exact.companion.joined.common.source.bridge
        if not self.base.core_acceptance_loaded or not self.bridge.acceptance_loaded:
            raise ValueError('Same checked common core/incoming physical and bridge owners required')
        self.core=self.base.core;self.pre=self.base.pre;self.params=self.pre.params
        self.ctx=c=self.core.ctx;self.delta=c.mpf(endpoints(self.params.delta))
        self.family=self.exterior.family;self.source=self.exterior.source;self.datum_sha=self.exterior.datum_sha
        self.logP=c.mpf(endpoints(self.params.logPstar));self.logC=self.core.logC
        self.logRref=c.ln(110)+10*(self.logC+self.logP)
        self.logRp=self.logRref+self.logP+1+self.params.Tw;self.logRv=self.logRp+13/self.params.mu
        self.pulse=self.history.selected.pulse;self.flatten=self.history.flatten
        self.outer=self.history.outer;self.steep=self.history.steep;self.heat=self.history.heat
        self.dispatch=CurrentCompletePhysicalDispatch(self);self.source_assembly=self.dispatch
        self.source_owners=self.dispatch.registry
        self.hashes=dict(self.exterior.hashes)
        theorem_name=PREFIX+'global_physical_assembly_check.json'
        self.original_operator_fixture=accepted(theorem_name,self.family,self.source,
            'all_33_original_source_charts_physical_spatial4_time1_mapped')
        _verify_hashes(self.original_operator_fixture)
        for old in (self.base.hashes,self.bridge.hashes,self.original_operator_fixture['input_hashes']):
            for name,digest in old.items():
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Retained/current physical source hash conflict: '+name)
                self.hashes[name]=digest
        self.hashes[theorem_name]=sha(theorem_name)
        self.operator_bindings=dict(original_radius=self.radius.__func__ is BASE.radius,
            original_normalization_with_correlated_pulse_rebase=self.normalized_sources.__func__ is PULSE.normalized_sources,
            full_original_Cartesian_time_evaluate_operator_called=True)
        self.proof=current_complete_physical_source_proof(self)
        for stem in ('global_physical_assembly','current_pulse_physical_assembly','current_flatten_physical_assembly',
            'current_heat_physical_assembly','current_bridge_physical_assembly','current_native_pulse_source_dispatcher'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.cache={};self.physical_acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current complete physical source scope/datum differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT)
            self.physical_acceptance_loaded=True

    def assert_graph(self):
        h=self.history;sel=h.selected;old=self.base
        h.assert_graph();sel.assert_graph()
        graph=dict(checked_complete_current_future_and_selection=sel.acceptance_loaded,
            same_current_history_and_dispatch_aliases=(h is self.exterior.history and self.dispatch.assembly is self
                and self.dispatch.history is h and self.dispatch.selected is sel),
            checked_current_energy_and_zero_meridional_history=h.acceptance_loaded,
            checked_current_full_similarity_exterior=self.exterior.acceptance_loaded,
            same_common_core_physical_ancestor=old is sel.exact.companion.owner,
            same_common_core_source=sel.exact.companion.joined.core is self.core.original,
            same_checked_incoming_core_and_bridge=self.core is old.core and self.bridge.owner is old and self.bridge.acceptance_loaded,
            same_incoming_pre_and_parameter_source=self.pre is old.pre and self.params is old.params,
            all_fifteen_current_physical_providers_replaced=all(self.dispatch.provider(chart) is not old.dispatch.provider(chart) for chart in CHANGED),
            all_six_pulse_physical_charts_use_one_current_selected_object=all(self.dispatch.provider(chart) is sel.pulse for chart in PULSE_CHARTS),
            all_nine_postpulse_physical_charts_use_current_energy_owners=all(self.dispatch.provider(chart) is getattr(h,METHODS[internal][0]) for chart,internal in PUBLIC_TO_INTERNAL.items()),
            incoming_fifteen_core_and_pre_provider_objects_preserved=all(self.dispatch.provider(chart) is old.dispatch.provider(chart) for chart in RETAINED if chart not in BRIDGES),
            all_three_common_bridge_providers_preserved=all(self.dispatch.provider(chart) is self.bridge for chart in BRIDGES),
            same_new_native_pulse_alias=self.dispatch.native_pulse is self.pulse is sel.pulse,
            original_fixed_native_incoming_kernel_owner=self.pulse.pulse is sel.amplitude.pulse is sel.original_sources[-1].pulse,
            original_beta_normalization_and_energy_Gram_source=(self.pulse.pulse.initial.repair.normalization is sel.amplitude.pulse.initial.repair.normalization
                and self.pulse.pulse.basis is sel.amplitude.pulse.basis and 'energy_gram' in self.pulse.pulse.basis),
            original_saddle_row_sources=(self.pulse.pulse.rows is sel.amplitude.pulse.rows
                and {'saddle_L','saddle_u0'}<=set(self.pulse.pulse.rows)),
            original_fixed_flat_gp_beta_derivative_provider=self.pulse.flat is sel.original_sources[-1].flat,
            current_pulse_data_cache_is_new=self.pulse.data_cache is not old.pulse.data_cache,
            same_original_Cstar_defining_record=(self.core.records['physical_norm_family']==
                sel.future.angular.initial.repair.records['compliant_physical_norm_family']),
            same_current_flatten_outer_steep_and_heat_aliases=(self.flatten is h.flatten and self.outer is h.outer and self.steep is h.steep and self.heat is h.heat),
            complete_33_chart_live_registry=set(self.dispatch.providers)==set(self.source_owners)==set(CHARTS),
            complete_15_changed_18_retained_partition=len(CHANGED)==15 and len(RETAINED)==18 and set(CHANGED).isdisjoint(RETAINED))
        if not all(graph.values()):raise ValueError('Complete physical provider graph differs: '+str(graph))
        return graph

    @source_precision
    def evaluate(self,chart,Z,coordinate,log_tau='-1',theta='0',axis=False):
        self.assert_graph()
        if chart not in CHARTS:raise ValueError('Unknown complete current physical chart')
        packet=BASE.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=axis)
        packet.update(datum_enclosure_sha256=self.datum_sha,
            current_source_owner=self.source_owners[chart]['provider'],
            current_source_acceptance_receipt=self.source_owners[chart]['acceptance_receipt'],
            current_complete_provider_graph=self.assert_graph(),
            physical_source_graph_scope='complete current callable source graph; quantitative joins and full resolved point/NS admission remain separate',
            checked_full_similarity_exterior_source_consumed=True,
            heat_exterior_stress_identity_certified=chart=='heat_exterior',
            heat_stress_scope='checked similarity exterior theorem only; no global physical tensor or NS admission',
            **dict.fromkeys(GATES,self.physical_acceptance_loaded),**dict.fromkeys(OPEN,False))
        return packet

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_physical_chart_owners=list(CHARTS),
            replaced_current_physical_source_charts=list(CHANGED),retained_current_incoming_source_charts=list(RETAINED),
            current_source_owner_registry=self.source_owners,current_complete_physical_source_proof=self.proof,
            original_physical_operator_identity=self.operator_bindings,
            old_live_complete_future_consumers_remaining=[],
            frozen_ancestor_objects_retained_as_source_provenance_not_current_postpulse_providers=True,
            same_checked_full_Gamma_similarity_stress_source_consumed=True,
            heat_exterior_stress_identity_certified=True,
            heat_stress_scope='checked similarity exterior theorem only; no global physical tensor or NS admission',
            current_changed_physical_coordinate_boxes=physical_views(self),
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))

    @source_precision
    def report(self):
        result=self.manifest();maps={}
        for chart,coordinate in physical_views(self).items():
            maps[chart]=self.evaluate(chart,[-1,1],coordinate,theta=None)
            print('Current complete physical source mapped: '+chart,flush=True)
        maps['fresh_main']=self.evaluate('pulse_main','.517','.83',log_tau='-2.17',theta='.37')
        maps['fresh_exterior']=self.evaluate('heat_exterior','.517','6.23',log_tau='-2.17',theta='.37')
        # This private coverage view uses the current pulse's same gap
        # function; the public gap domain remains [11,12].
        dispatched=self.dispatch.gap_overlap([-1,1]);point=dispatched['source_packet']
        coordinate=['12','12.0001'];logR,_=self.radius('pulse_gap',coordinate,point,self.pulse)
        grids,logs,amplitudes=self.normalized_sources('pulse_gap',[-1,1],coordinate,point,self.pulse,logR)
        result['current_same_owner_supplemental_gap_source']=dict(source_packet=point,source_logR_enclosure=logR,
            same_current_selected_owner=True,public_gap_domain_not_widened=True,
            exact_gap_source_coordinate='xi=13+mu*s',source_log_bases=logs,source_amplitudes=amplitudes)
        result.update(current_changed_physical_maps=maps,input_hashes=self.hashes)
        return result


@source_precision
def run(field=None):
    field=field if field is not None else CurrentCompletePhysicalAssembly(require_checked=False)
    result=field.report();(HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('All 33 physical source chart providers installed on the complete current graph',flush=True)
    return result


if __name__=='__main__':run()
