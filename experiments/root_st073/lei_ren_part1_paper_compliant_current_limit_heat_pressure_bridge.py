"""Same repaired-limit pulse/flatten objects through preheat and pressure.

Transport the admitted original defining functions into a downstream-only
runtime. No old core/interface constructor is called and no global smooth
field, physical chart installation or selected future-energy claim is made.
"""
import json
import time
from pathlib import Path

import lei_ren_part1_paper_compliant_current_limit_native_pulse_bridge as inlet
import lei_ren_part1_paper_compliant_current_power_angular_source as power
import lei_ren_part1_paper_compliant_current_steep_waiting_source as steep
import lei_ren_part1_paper_compliant_current_heat_source as heat
import lei_ren_part1_paper_compliant_current_heat_pressure_stress as companion
import lei_ren_part1_paper_compliant_current_exact_repair_branch as repair
import lei_ren_part1_paper_compliant_current_collar_stress_mixed_C4 as collar
import lei_ren_part1_paper_compliant_current_angular_terminal_closure as angular
import lei_ren_part1_paper_compliant_current_pressure_terminal_balance as balance
import lei_ren_part1_paper_compliant_current_raw_preheat_pressure_operator as raw
import lei_ren_part1_paper_compliant_current_pressure_terminal_closure as pressure
from lei_ren_part1_paper_compliant_pressure_source import CompliantPressureDatum, CompliantOuterParameters
from lei_ren_part1_paper_compliant_fifth_axial_jets import CompliantFifthAxialJets
from lei_ren_part1_paper_compliant_future_swirl_energy import CompliantFutureSwirlEnergy
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE,PREFIX,sha,require=inlet.HERE,inlet.PREFIX,inlet.sha,inlet.require
NAME=PREFIX+'current_limit_heat_pressure_bridge.json'
RECEIPT=PREFIX+'current_limit_heat_pressure_bridge_check.json'
GATE='current_repaired_limit_same_object_preheat_heat_pressure_runtime_installed'
OPEN=tuple(dict.fromkeys(inlet.OPEN+pressure.OPEN+(
    'current_exact_repair_installed_in_all_physical_charts',
    'current_selected_complete_future_energy_installed_in_pulse',
    'current_heat_physical_owner_installed','global_completed_tensor_admissibility')))


def admit(stem,gate,owner):
    return inlet.checked(stem,gate,owner.family_record,owner.hashes)


def metadata(before,owner):
    owner.before=before
    owner.family_record=before.family_record
    owner.family=before.family;owner.source=before.source;owner.datum_sha=before.datum_sha
    owner.hashes=dict(before.hashes)
    owner.hashes[Path(__file__).name]=sha(Path(__file__).name)


def native_prefix_binding(before,outer,old):
    """Bind actual native objects; no current-core compatibility shim."""
    f=outer.future;p=outer.pulse;fifth=outer.fifth
    datums=(p.pulse.initial.datum,p.selection.future.angular.initial.datum,
        f.angular.initial.datum,fifth.angular4.repair.angular.initial.datum)
    require(type(f) is CompliantFutureSwirlEnergy and type(fifth) is CompliantFifthAxialJets,
        'Original native future and fifth-jet classes required')
    for datum in datums:
        require(type(datum) is CompliantPressureDatum and type(datum.parameters) is CompliantOuterParameters,
            'Original analytic datum and parameter classes required')
        require(datum.parameters.Md=='40' and datum.parameters.precision==160,
            'Pinned Md40 parameter family required')
        require(type(datum).normalized_jets is CompliantPressureDatum.normalized_jets,
            'Same inherited analytic pressure callable required')
        for key in ('definition','source_sha','datum_sha','input_hashes'):
            require(getattr(datum,key)==getattr(datums[0],key),'Native defining datum differs: '+key)
        require(datum.source_sha==before.source and datum.datum_sha==before.datum_sha,
            'Repaired-limit/native defining source differs')
    original=old['current_native_parameter_source_bridge']
    require(original['passed'] and original['exact_Md40_delta_choice_branch_proved'],
        'Unchanged original exact parameter source theorem required')
    graph=dict(same_bridge_pulse=p is before.pulse,same_bridge_flatten=outer.flatten is before.flatten,
        same_fifth_object=fifth is p.fifth,same_prefix_future=f is fifth.fourth.energy.base,
        future_parameter_aliases=f.params is f.repair.params is f.heat.params is f.angular.params is datums[2].parameters,
        angular_C4_parameter_aliases=fifth.angular4.repair.params is datums[3].parameters,
        prefix_mu_aliases=f.mu is f.params.mu is f.repair.mu,
        prefix_delta_aliases=f.delta is f.repair.delta is f.heat.delta is f.angular.delta is f.angular.initial.delta,
        original_power_callable=outer.power.__func__ is power.BASE.power,
        original_angular_callable=outer.angular.__func__ is power.BASE.angular,
        original_angular_C5_callable=fifth.angular.__func__ is CompliantFifthAxialJets.angular,
        same_native_provider_context=outer.ctx is p.ctx,
        same_native_mu=outer.mu is p.mu,same_native_delta=outer.delta is p.delta)
    require(all(graph.values()),'Repaired-limit native prefix object graph differs')
    return dict(passed=True,actual_native_only_object_graph=graph,
        original_exact_parameter_formula_theorem_reused=True,
        actual_mu_delta_defining_assignments=original['actual_mu_delta_defining_assignments'],
        current_Rp_defining_function_receipt=PREFIX+'current_limit_Rp_native_identity_check.json',
        current_core_compatibility_not_asserted=True,
        input_hashes=original['input_hashes'])


class CurrentLimitPowerSource(power.CurrentPowerAngularSourceAssembly):
    @source_precision
    def __init__(self,before):
        metadata(before,self)
        require(before.acceptance_loaded,'Checked repaired-limit inlet required')
        self.outer=power.CurrentPowerAngularC4(before);self.hashes.update(self.outer.hashes)
        admit('current_power_angular_source',power.GATE,self)
        old=json.loads((HERE/(PREFIX+'current_power_angular_source.json')).read_bytes())
        require(PREFIX+'current_power_angular_source.json' in self.hashes,
            'Original power source report must be pinned by its receipt')
        theorem=admit('power_angular_C4','flatten_power_and_power_angular_joins_certified',self)
        self.functional_proof=power.functional_source_identities()
        require(self.functional_proof==theorem['functional_production_source_identities'],
            'Original power/angular defining functions changed')
        old_binding=old['current_power_angular_source_bindings']
        self.bindings=dict(current_native_parameter_source_bridge=native_prefix_binding(before,self.outer,old_binding),
            current_prefix_and_future_decomposition=power.exact_prefix_and_decomposition_bindings())
        require(self.bindings['current_prefix_and_future_decomposition']['passed'],
            'Actual C4/C5 complete future decomposition required')
        self.registry=dict(before.chain_routes)
        self.registry['flatten']=dict(provider=PREFIX+'current_pulse_flatten_source.CurrentFlattenMixedC4',
            method='flatten',coverage_coordinate='t=log(R/Rv)',domain='[0,100]',acceptance_receipt=inlet.RECEIPT)
        self.registry.update({chart:dict(provider=PREFIX+'current_power_angular_source.CurrentPowerAngularC4',
            method='power' if chart=='outer_power' else 'angular',
            coverage_coordinate='phase=(logR-logRf)/(Lrel-4)' if chart=='outer_power' else 's=logR-logRrel',
            domain='[0,1]' if chart=='outer_power' else '[-4,0]',acceptance_receipt=RECEIPT)
            for chart in power.NEW_CHARTS})
        self.acceptance_loaded=True


class CurrentLimitSteepSource(steep.CurrentSteepWaitingSourceAssembly):
    @source_precision
    def __init__(self,before):
        metadata(before,self)
        self.steep=steep.CurrentSteepWaitingC4(before);self.hashes.update(self.steep.hashes)
        admit('current_steep_waiting_source',steep.GATE,self)
        theorem=admit('steep_waiting_C4','angular_steep_and_internal_joins_certified',self)
        self.functional_proof=steep.functional_source_identities()
        require(self.functional_proof==theorem['functional_production_source_identities'],
            'Original steep/waiting defining functions changed')
        self.bindings=steep.current_source_bindings(self)
        self.hashes.update(self.bindings['exact_current_kernel_and_parameter_source_bindings']['input_hashes'])
        self.registry=dict(before.registry)
        self.registry.update({chart:dict(provider=PREFIX+'current_steep_waiting_source.CurrentSteepWaitingC4',
            method=steep.METHODS[chart],coverage_coordinate='t=logR-logRorigin' if chart in ('steep_entry','steep_exit') else 'phase=offset/original_length',
            domain='[0,1]',acceptance_receipt=RECEIPT) for chart in steep.NEW_CHARTS})
        self.acceptance_loaded=True


class CurrentLimitHeatSource(heat.CurrentHeatSourceAssembly):
    @source_precision
    def __init__(self,before):
        metadata(before,self)
        self.heat=heat.CurrentCollarGammaC4(before);self.hashes.update(self.heat.hashes)
        admit('current_heat_source',heat.GATE,self)
        theorem=admit('collar_Gamma_C4','waiting_collar_and_collar_Gamma_joins_certified',self)
        gamma=admit('heat_pressure_C4','absolute_pressure_same_source_mixed4_available',self)
        self.functional_proof=heat.functional_source_identities()
        require(self.functional_proof==theorem['functional_production_source_identities'],
            'Original collar/Gamma defining functions changed')
        self.gamma_evidence=gamma['defining_source_bridge']['defining_function_bridge']['gamma_derivative_enclosure']
        require(self.gamma_evidence['verified'],'Original complete Gamma enclosure theorem required')
        self.bindings=heat.current_source_bindings(self)
        self.registry=dict(before.registry)
        self.registry.update({chart:dict(provider=PREFIX+'current_heat_source.CurrentCollarGammaC4',
            method=heat.METHODS[chart],coverage_coordinate='t=log(R/Rtail)',
            domain='[0,3]' if chart=='heat_collar' else '[3,infinity)',acceptance_receipt=RECEIPT)
            for chart in heat.NEW_CHARTS})
        self.acceptance_loaded=True


class CurrentLimitHeatCompanion(companion.CurrentHeatPressureStress):
    """Original downstream pressure/stress algorithms, no common-core claim."""
    @source_precision
    def __init__(self,source_owner):
        self.source_owner=source_owner;self.family_record=source_owner.family_record
        self.heat=source_owner.heat;self.ctx=self.heat.ctx
        self.family=source_owner.family;self.source=source_owner.source;self.datum_sha=source_owner.datum_sha
        self.hashes=dict(source_owner.hashes)
        original=admit('current_heat_pressure_stress',companion.GATES[0],self)
        # Reuse only source-independent formula theorems from this receipt.
        # Its old common-core object graph is deliberately not transported.
        names=('actual_retained_forward_pressure_proof','actual_constant_stress_proof',
            'current_actual_pressure_source_split_proof')
        self.formula_evidence={name:original[name] for name in names}
        require(all(row['passed'] for row in self.formula_evidence.values()),
            'Original pressure/stress formula evidence required')
        self.graph=dict(checked_repaired_limit_heat_source=source_owner.acceptance_loaded,
            same_current_heat_source=self.heat is source_owner.heat,
            same_current_waiting=self.heat.steep is source_owner.before.steep,
            same_current_future=self.heat.future is self.heat.steep.future,
            same_exact_raw_heat=self.heat.exact_heat is self.heat.future.heat,
            same_actual_native_context=self.ctx is self.heat.steep.ctx is self.heat.outer.ctx,
            original_terminal_constants_callable=self.terminal_constants.__func__ is companion.CurrentHeatPressureStress.terminal_constants,
            original_pressure_stress_callable=self.evaluate.__func__ is companion.CurrentHeatPressureStress.evaluate)
        require(all(self.graph.values()),'Actual repaired-limit heat companion graph differs')
        self.current_bindings=heat.current_source_bindings(source_owner)
        self.history_bindings=companion.history_source_bindings()
        self.constants={};self.runtime={};self.acceptance_loaded=True
        self.common_core_interface_claim=False


class CurrentLimitPressureClosure(pressure.CurrentPressureTerminalClosure):
    @source_precision
    def __init__(self,balance,raw,family_record):
        super().__init__(balance=balance,raw=raw,require_checked=False)
        self.family_record=family_record
        original=admit('current_pressure_terminal_closure',pressure.GATES[1],self)
        old=dict(original['original_pressure_function_identification'])
        current=dict(encode(pack(self.proof)))
        old_parameter=old.pop('common_exact_parameter_function_bridge')
        new_parameter=current.pop('common_exact_parameter_function_bridge')
        require(old==current,'Original pressure function proof changed beyond parameter transport')
        require(new_parameter['passed'] and new_parameter['current_core_compatibility_not_asserted'],
            'Actual native-only parameter source bridge required')
        require(old_parameter['actual_mu_delta_defining_assignments']==
            new_parameter['actual_mu_delta_defining_assignments'],
            'Original exact parameter defining functions differ')
        require(self.raw.angular is self.balance.angular and self.raw.exact is self.balance.exact,
            'Raw and repaired pressure branches must share one actual source')
        self.transport=dict(only_changed_proof_component='common_exact_parameter_function_bridge',
            original_fourteen_stage_and_P0_proof_unchanged=True,
            original_exact_mu_delta_formulas_unchanged=True,
            old_core_parameter_graph_not_imported=True,
            actual_native_parameter_graph_substituted=True,passed=True)
        self.acceptance_loaded=True


class CurrentLimitHeatPressureBridge:
    @source_precision
    def __init__(self,bridge=None,require_checked=True):
        self.inlet=bridge if bridge is not None else inlet.CurrentLimitNativePulseBridge()
        require(self.inlet.acceptance_loaded,'Checked same-object repaired-limit pulse bridge required')
        self.family_record=self.inlet.family_record
        self.power=CurrentLimitPowerSource(self.inlet)
        self.steep=CurrentLimitSteepSource(self.power)
        self.heat_source=CurrentLimitHeatSource(self.steep)
        self.companion=CurrentLimitHeatCompanion(self.heat_source)
        self.exact=repair.CurrentExactRepairBranch(companion=self.companion)
        self.collar=collar.CurrentCollarStressMixedC4(companion=self.companion)
        self.angular=angular.CurrentAngularTerminalClosure(exact=self.exact,collar=self.collar)
        self.balance=balance.CurrentPressureTerminalBalance(angular=self.angular)
        self.raw=raw.CurrentRawPreheatPressureOperator(angular=self.angular)
        self.pressure=CurrentLimitPressureClosure(balance=self.balance,raw=self.raw,family_record=self.family_record)
        self.family=self.pressure.family;self.source=self.pressure.source;self.datum_sha=self.pressure.datum_sha
        self.hashes=dict(self.pressure.hashes);self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.graph=self.object_graph()
        require(all(self.graph.values()),'Actual downstream pressure source graph differs')
        self.acceptance_loaded=False
        if require_checked:
            admit('current_limit_heat_pressure_bridge',GATE,self)
            self.acceptance_loaded=True

    def object_graph(self):
        return dict(power_consumes_repaired_limit_inlet=self.power.before is self.inlet,
            exact_same_pulse=self.exact.pulse is self.power.outer.pulse is self.inlet.pulse,
            exact_same_flatten=self.exact.flatten is self.power.outer.flatten is self.inlet.flatten,
            steep_consumes_same_power=self.steep.steep.outer is self.power.outer,
            heat_consumes_same_steep=self.heat_source.heat.steep is self.steep.steep,
            companion_consumes_same_heat=self.companion.source_owner is self.heat_source,
            exact_and_collar_share_companion=self.exact.companion is self.collar.companion is self.companion,
            angular_consumes_same_exact=self.angular.exact is self.exact,
            balance_and_raw_share_angular=self.balance.angular is self.raw.angular is self.angular,
            closure_consumes_actual_balance_and_raw=self.pressure.balance is self.balance and self.pressure.raw is self.raw,
            unchanged_independent_P0_source=self.pressure.datum_sha==self.inlet.datum_sha,
            one_replayed_unique_repair=self.exact.angular4.repair is self.exact.future.repair,
            same_replayed_Gamma_radius=self.angular.heat is self.raw.angular.heat is self.pressure.heat,
            original_P0_raw_function_identification=self.pressure.proof['passed'],
            original_pressure_balance=self.balance.proof['passed'],original_angular_closure=self.angular.proof['passed'])

    @source_precision
    def evaluate(self,chart,Z,t):
        packet=self.pressure.evaluate(chart,Z,t)
        return dict(source_packet=packet,source_family=self.family_record,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False),
            output_kind='same-object downstream source enclosures; no numeric point selection')


def run():
    began=time.monotonic();field=CurrentLimitHeatPressureBridge(require_checked=False)
    result=dict(source_family=field.family_record,**{GATE:False},all_passed=True,
        candidate_same_object_runtime_constructed=True,
        actual_same_object_graph=field.graph,actual_native_prefix_binding=field.power.bindings,
        current_heat_source_binding=field.heat_source.bindings,
        downstream_companion_graph=field.companion.graph,
        reused_downstream_only_companion_formula_evidence=field.companion.formula_evidence,
        current_pressure_parameter_transport=field.pressure.transport,
        ordered_downstream_registry=field.heat_source.registry,
        reused_original_pressure_function_identification=field.pressure.proof,
        reused_original_pressure_balance=field.balance.proof,
        reused_original_angular_closure=field.angular.proof,
        common_core_interface_not_reconstructed=True,
        old_all_N_inner_producers_not_replayed=True,
        source_enclosures_not_defining_point_values=True,
        current_exact_unique_branch_not_installed_in_old_selected_pulse=True,
        **dict.fromkeys(OPEN,False))
    result['actual_heat_handshakes']={}
    for name,chart,t in (('preheat','heat_collar',0),('exterior','heat_exterior',3)):
        result['actual_heat_handshakes'][name]=encode(pack(field.evaluate(chart,0,t)))
        print('Current limit actual pressure/heat handshake: '+name,flush=True)
    result['input_hashes']=field.hashes;result['execution_seconds']=time.monotonic()-began
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
