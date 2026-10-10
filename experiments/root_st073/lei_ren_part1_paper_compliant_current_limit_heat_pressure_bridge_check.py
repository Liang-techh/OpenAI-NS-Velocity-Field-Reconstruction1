"""Check new same-object injections; reuse admitted downstream mathematics."""
import json
import time
from pathlib import Path
import lei_ren_part1_paper_compliant_current_limit_heat_pressure_bridge as bridge


def zero(row):
    return row['lower_exact_mpf_tuple'][1]==row['upper_exact_mpf_tuple'][1]==0


def run():
    began=time.monotonic();report=json.loads((bridge.HERE/bridge.NAME).read_bytes())
    assert report['all_passed'] and report['candidate_same_object_runtime_constructed']
    assert report[bridge.GATE] is False
    hashes=dict(report['input_hashes'])
    for name,digest in hashes.items():assert bridge.sha(name)==digest,name
    family=report['source_family']
    bridge.inlet.checked('current_limit_native_pulse_bridge',bridge.inlet.GATE,family,hashes)
    old=bridge.inlet.checked('current_pressure_terminal_closure',bridge.pressure.GATES[1],family,hashes)
    current=dict(report['reused_original_pressure_function_identification'])
    accepted=dict(old['original_pressure_function_identification'])
    native_parameter=current.pop('common_exact_parameter_function_bridge')
    old_parameter=accepted.pop('common_exact_parameter_function_bridge')
    assert current==accepted
    assert native_parameter['actual_mu_delta_defining_assignments']==old_parameter['actual_mu_delta_defining_assignments']
    assert native_parameter['current_core_compatibility_not_asserted']
    assert report['current_pressure_parameter_transport']['passed']
    assert all(report['actual_same_object_graph'].values())
    native=report['actual_native_prefix_binding']['current_native_parameter_source_bridge']
    assert native['passed'] and all(native['actual_native_only_object_graph'].values())
    assert native['current_core_compatibility_not_asserted']
    assert report['current_heat_source_binding']['passed']
    assert all(report['downstream_companion_graph'].values())
    assert report['reused_original_pressure_balance']['passed'] and report['reused_original_angular_closure']['passed']
    bindings={}
    specs={
        'CurrentLimitPowerSource':{'self.outer':'power.CurrentPowerAngularC4(before)'},
        'CurrentLimitSteepSource':{'self.steep':'steep.CurrentSteepWaitingC4(before)'},
        'CurrentLimitHeatSource':{'self.heat':'heat.CurrentCollarGammaC4(before)'},
        'CurrentLimitHeatCompanion':{'self.source_owner':'source_owner','self.heat':'source_owner.heat'},
        'CurrentLimitHeatPressureBridge':{
            'self.power':'CurrentLimitPowerSource(self.inlet)',
            'self.steep':'CurrentLimitSteepSource(self.power)',
            'self.heat_source':'CurrentLimitHeatSource(self.steep)',
            'self.companion':'CurrentLimitHeatCompanion(self.heat_source)',
            'self.exact':'repair.CurrentExactRepairBranch(companion=self.companion)',
            'self.collar':'collar.CurrentCollarStressMixedC4(companion=self.companion)',
            'self.angular':'angular.CurrentAngularTerminalClosure(exact=self.exact,collar=self.collar)',
            'self.balance':'balance.CurrentPressureTerminalBalance(angular=self.angular)',
            'self.raw':'raw.CurrentRawPreheatPressureOperator(angular=self.angular)',
            'self.pressure':'CurrentLimitPressureClosure(balance=self.balance,raw=self.raw,family_record=self.family_record)'}}
    for cls,rows in specs.items():
        bindings[cls]=bridge.inlet.identity.ast_assignments('current_limit_heat_pressure_bridge',cls,'__init__',rows)
    # Independently inspect the live source chain, not its reported object_graph.
    supplied=bridge.inlet.CurrentLimitNativePulseBridge()
    field=bridge.CurrentLimitHeatPressureBridge(bridge=supplied,require_checked=False)
    assert field.inlet is supplied
    assert field.power.before is supplied
    assert field.power.outer.pulse is supplied.pulse is field.exact.pulse
    assert field.power.outer.flatten is supplied.flatten is field.exact.flatten
    assert field.steep.steep.outer is field.power.outer
    assert field.heat_source.heat.steep is field.steep.steep
    assert field.companion.source_owner is field.heat_source
    assert not field.companion.common_core_interface_claim
    assert not hasattr(field.companion,'joined') and not hasattr(field.companion,'owner')
    assert 'current_common_heat_source_graph' not in field.companion.formula_evidence
    for name,row in field.companion.formula_evidence.items():
        assert row['passed'] and row==report['reused_downstream_only_companion_formula_evidence'][name]
    assert field.exact.companion is field.collar.companion is field.companion
    assert field.angular.exact is field.exact and field.angular.collar is field.collar
    assert field.balance.angular is field.raw.angular is field.angular
    assert field.pressure.balance is field.balance and field.pressure.raw is field.raw
    assert field.raw.exact.flatten is supplied.flatten
    assert field.exact.angular4.repair is field.exact.future.repair
    assert field.angular.heat is field.raw.angular.heat is field.pressure.heat
    expected=set(bridge.inlet.PULSE_CHARTS)|{'flatten'}|set(bridge.power.NEW_CHARTS)|set(bridge.steep.NEW_CHARTS)|set(bridge.heat.NEW_CHARTS)
    assert set(report['ordered_downstream_registry'])==expected==set(field.heat_source.registry)
    assert len(expected)==15
    p0=None
    for name,view in report['actual_heat_handshakes'].items():
        packet=view['source_packet']
        assert packet['pressure_is_original_forward_function_after_terminal_identity']
        assert packet['original_axis_pressure_not_replaced_or_tail_patched']
        assert packet['whole_Z_pressure_function_identity_consumed']
        assert packet['pressure_scope_only_full_exterior_receipt_still_required']
        for key in ('actual_Dtheta_Taylor','actual_Cp_Taylor',
                    'source_proved_original_P0_plus_native_raw_integral_Taylor',
                    'source_proved_original_pressure_balance_Taylor'):
            rows=packet[key]['coefficients'];assert len(rows)==6 and all(zero(v) for v in rows)
        datum=packet['original_analytic_P0_Taylor_retained']
        assert len(datum['coefficients'])==6 and not zero(datum['coefficients'][0])
        if p0 is None:p0=datum
        else:assert datum==p0
        for key in bridge.pressure.GATES:assert packet[key]
        for key in bridge.OPEN:assert view[key] is False
    for key in bridge.OPEN:assert report[key] is False
    assert report['common_core_interface_not_reconstructed'] and report['old_all_N_inner_producers_not_replayed']
    assert report['source_enclosures_not_defining_point_values']
    assert report['current_exact_unique_branch_not_installed_in_old_selected_pulse']
    hashes[bridge.NAME]=bridge.sha(bridge.NAME);hashes[Path(__file__).name]=bridge.sha(Path(__file__).name)
    result=dict(source_family=family,all_passed=True,**{bridge.GATE:True},
        independent_live_same_object_chain_checked=True,
        original_downstream_algorithms_injected_with_actual_bridge_objects=True,
        all_fifteen_downstream_routes_retained=True,
        strict_original_pressure_function_receipt_reused=True,
        original_nonzero_independent_P0_retained_at_both_heat_handshakes=True,
        source_proved_terminal_constants_zero_not_box_overlap=True,
        common_core_shim_absent=True,scope_not_promoted=True,
        actual_constructor_AST_bindings=bindings,**dict.fromkeys(bridge.OPEN,False),
        input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (bridge.HERE/bridge.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('PASS_CURRENT_LIMIT_HEAT_PRESSURE_BRIDGE same live pulse/flatten through 15 routes and pressure closure',flush=True)
    return result


if __name__=='__main__':run()
