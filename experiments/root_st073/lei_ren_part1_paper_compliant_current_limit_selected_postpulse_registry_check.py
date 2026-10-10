"""Check this branch injection and registry; retain existing integral proofs."""
import json
import time
from pathlib import Path

import lei_ren_part1_paper_compliant_current_limit_selected_postpulse_registry as bridge


def zero(row):
    return row['lower_exact_mpf_tuple'][1]==row['upper_exact_mpf_tuple'][1]==0


@bridge.source_precision
def run(field=None):
    began=time.monotonic();raw=json.loads((bridge.HERE/bridge.NAME).read_bytes())
    assert raw['all_passed'] and raw['candidate_selected_postpulse_registry_constructed']
    assert not any(raw[key] for key in bridge.GATES+bridge.OPEN)
    hashes=dict(raw['input_hashes']);family=raw['source_family']
    for name,digest in hashes.items():assert bridge.sha(name)==digest,name
    old_history=bridge.selected.source.inlet.checked('current_postpulse_energy_history',bridge.history.GATES[0],family,hashes)
    old_exterior=bridge.selected.source.inlet.checked('current_full_exterior_stress',bridge.exterior.GATES[0],family,hashes)
    bridge.selected.source.inlet.checked('current_limit_selected_pulse_flatten',bridge.selected.GATE,family,hashes)
    field=field if field is not None else bridge.CurrentLimitSelectedPostpulseRegistry(require_checked=False)
    assert field.family_record==family and not field.acceptance_loaded
    h=field.history;sel=field.selected_owner.selected;pressure=field.selected_owner.source_field.pressure
    # Inspect the actual injected objects as well as binding their constructors.
    assert h.selected is sel and h.flatten.pulse is h.outer.pulse is sel.pulse
    assert h.outer.flatten is h.flatten and h.outer.fifth is sel.fifth
    assert h.outer.future is h.steep.future is h.heat.future is sel.future
    assert h.heat.repair is sel.exact.repair is field.selected_owner.source_field.exact.repair
    assert h.heat.exact_heat is sel.future.heat
    assert field.exterior.history is h and field.exterior.heat is h.heat
    assert field.exterior.pressure is sel.pressure is pressure
    assert h.flatten is not field.selected_owner.flatten
    assert h.flatten.inlet.datum is sel.future.angular.initial.datum is sel.exact.repair.angular.initial.datum
    assert bridge.serialized(field.object_graph())==raw['actual_same_object_graph'] and all(field.object_graph().values())
    assert bridge.serialized(h.proof)==raw['unchanged_energy_source_proof']
    assert bridge.serialized(field.exterior.proof)==raw['current_full_terminal_history_transfer']
    transport=bridge.transport_proofs(h,field.exterior,old_history,old_exterior)
    assert bridge.serialized(transport)==raw['strict_original_source_proof_transport'] and transport['passed']
    assert h.proof['selected_forward_cumulative_energy_equals_same_remaining_integral_by_FTC']
    assert h.proof['zero_meridional_histories_propagate_from_selected_terminal_by_FTC']
    assert field.exterior.proof['all_current_five_terminal_moment_functions_identified']
    bindings=bridge.selected.source.inlet.identity.ast_assignments('current_limit_selected_postpulse_registry',
        'CurrentLimitSelectedPostpulseRegistry','__init__',{
            'self.history':'history.CurrentPostpulseEnergyHistory(selected=self.selected_owner.selected)',
            'self.exterior':'exterior.CurrentFullExteriorStress(history=self.history)',
            'self.proof_transport':'transport_proofs(self.history,self.exterior,old_history,old_exterior)',
            'self.registry':'self.route_registry()'})
    expected=set(bridge.selected.source.inlet.PULSE_CHARTS)|set(bridge.DOMAINS)
    assert len(expected)==15 and set(field.registry)==expected
    assert field.registry==raw['ordered_current_source_registry']
    for chart,route in field.registry.items():
        assert route['acceptance_receipt']==bridge.RECEIPT and route['selected_source_receipt']==bridge.selected.RECEIPT
        assert 'current_limit_selected_postpulse_registry.CurrentLimitSelectedPostpulseRegistry.' in route['provider']
        assert all(route[key]==value for key,value in family.items())
        assert route['input_hash_ledger']==bridge.RECEIPT+'#input_hashes'
        if chart in bridge.selected.source.inlet.PULSE_CHARTS:
            method,domain=bridge.selected.source.inlet.EXPECTED[chart]
            assert route['method']==method and route['domain']==domain
            assert route['coverage_coordinate']==bridge.ROUTES[chart][3]
            assert field.provider(chart) is sel.pulse
        else:
            native=bridge.ALIASES.get(chart,chart);owner,method=bridge.history.METHODS[native]
            assert route['method']==method
            assert field.provider(chart) is (field.exterior if chart=='heat_exterior' else getattr(h,owner))
    rejected=[]
    for chart,bad in (('flatten',-1),('flatten',101),('outer_power',-1),('outer_power',2),
        ('outer_angular',-5),('outer_angular',1),('steep_entry',2),('steep_power',2),
        ('steep_exit',2),('waiting',-1),('heat_collar',4),('heat_exterior',2)):
        try:field.evaluate(chart,0,bad)
        except ValueError:rejected.append((chart,bad));continue
        raise AssertionError('Out-of-domain route accepted: '+chart)
    for chart,bad in (('pulse_entrance',-1),('pulse_main',0),('pulse_exit',9),
        ('pulse_gap',10),('pulse_gap_end',0),('pulse_end',1)):
        try:field.evaluate(chart,0,bad)
        except ValueError:rejected.append((chart,bad));continue
        raise AssertionError('Out-of-domain pulse route accepted: '+chart)
    for Z in (-2,2):
        try:field.evaluate('flatten',Z,0)
        except ValueError:rejected.append(('Z',Z));continue
        raise AssertionError('Out-of-domain Z accepted')
    actual={}
    for chart,Z,coordinate in (('pulse_end','.427',0),('flatten','.427',0),
        ('outer_angular','.427',0),('waiting','.427',1),('heat_collar','.427',0),
        ('heat_exterior','.427','4.23')):
        point=field.evaluate(chart,Z,coordinate);actual[chart]=bridge.serialized(point)
        assert actual[chart]==raw['actual_source_handshakes'][chart]
        assert point['acceptance_receipt']==bridge.RECEIPT
        assert point['source_provider']==field.registry[chart]['provider']
        owner=field.provider(chart);method=getattr(owner,field.registry[chart]['method']).__func__
        provider=point['underlying_source_provider']
        assert provider['actual_method']==method.__module__+'.'+method.__qualname__
        assert provider['actual_method_source_file']==method.__module__+'.py'
        assert provider['actual_method_source_sha256']==hashes[provider['actual_method_source_file']]
        assert all(provider[key]==value for key,value in family.items())
        assert point['same_repaired_limit_selected_future_in_all_postpulse_routes']
        assert not any(point[key] for key in bridge.GATES+bridge.OPEN)
    rv=actual['pulse_end']['source_packet'];flat=actual['flatten']['source_packet']
    assert rv['pressure']['P0_over_Pstar_squared']==flat['pressure']['P0_over_Pstar_squared']
    assert not zero(flat['pressure']['P0_over_Pstar_squared']['coefficients'][0])
    for packet in actual.values():
        assert packet['original_analytic_P0_Taylor_retained']==flat['pressure']['P0_over_Pstar_squared']
    for chart in ('flatten','outer_angular','waiting','heat_collar','heat_exterior'):
        packet=field.evaluate(chart,'.427',0 if chart in ('flatten','outer_angular','heat_collar') else 1 if chart=='waiting' else '4.23')['source_packet']
        assert bridge.history.endpoints(packet['energy_Taylor'][0])[0]>0,chart
    ext=actual['heat_exterior']['source_packet']
    assert ext['all_current_five_terminal_histories_consumed']
    assert ext['actual_energy_is_positive_same_full_Gamma_future_half']
    assert ext['stress_zero_is_source_identity_not_interval_overlap']
    assert ext['heat_exterior_stress_identity_certified']
    assert 'heat_exterior_stress_identity_certified' not in bridge.OPEN
    for key in ('actual_Dtheta_Taylor','actual_Cp_Taylor'):
        assert len(ext[key]['coefficients'])==6 and all(zero(row) for row in ext[key]['coefficients'])
    for jet in ext['exterior_meridional_moment_Taylor'].values():assert all(zero(row) for row in jet['coefficients'])
    indices={'y%d_Z%d'%(j,n) for j in range(5) for n in range(5-j)}
    for rows in ext['exterior_stress_mixed4'].values():
        assert set(rows)==indices and all(zero(row) for row in rows.values())
    assert set(ext['stable_current_absolute_pressure_mixed4'])==indices
    assert raw['no_second_exact_repair_constructed'] and raw['original_selected_and_native_graphs_preserved']
    hashes[bridge.NAME]=bridge.sha(bridge.NAME);hashes[Path(__file__).name]=bridge.sha(Path(__file__).name)
    receipt=dict(source_family=family,all_passed=True,**dict.fromkeys(bridge.GATES,True),**dict.fromkeys(bridge.OPEN,False),
        actual_current_constructor_AST_bound=bindings,
        strict_source_proof_transport=transport,
        actual_same_selected_branch_future_pressure_and_exact_heat_checked=True,
        all_fifteen_routes_have_actual_provider_and_new_receipt=True,
        out_of_domain_source_routes_rejected=rejected,
        fresh_current_selected_postpulse_and_exterior_handshakes_checked=True,
        independent_nonzero_P0_and_positive_complete_energy_retained=True,
        current_full_exterior_source_identity_mixed4_stress_rows=30,
        absolute_post_2Rc_histories_and_global_physical_assembly_not_claimed=True,
        input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (bridge.HERE/bridge.RECEIPT).write_text(json.dumps(bridge.serialized(receipt),indent=2)+'\n',encoding='utf8')
    print('PASS_CURRENT_LIMIT_SELECTED_POSTPULSE_REGISTRY: 15 routes, same repair, full exterior source histories/stress',flush=True)
    return receipt


if __name__=='__main__':run()
