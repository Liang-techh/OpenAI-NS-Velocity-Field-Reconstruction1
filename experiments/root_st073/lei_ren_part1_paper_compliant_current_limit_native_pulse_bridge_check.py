"""Source admission and live handshake checks without replaying inner owners."""
import json
from pathlib import Path
import time
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_limit_native_pulse_bridge as bridge
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def rows(value):return value['coefficients']


def is_zero(value):
    return value['lower_exact_mpf_tuple'][1]==value['upper_exact_mpf_tuple'][1]==0


def run():
    began=time.monotonic();report=json.loads((bridge.HERE/bridge.NAME).read_bytes())
    assert report[bridge.GATE] and report['all_passed']
    hashes=dict(report['input_hashes'])
    for name,digest in hashes.items():assert bridge.sha(name)==digest,name
    family=report['source_family']
    bridge.checked('current_limit_Rp_native_identity',bridge.identity.GATE,family,hashes)
    bridge.checked('actual_Rp_source_join','current_Rp_external_pulse_join_certified',family,hashes)
    bridge.checked('current_pulse_flatten_source','current_pulse_terminal_flatten_source_ownership_certified',family,hashes)
    bridge.checked('pulse_end_flatten_join','pulse_end_flatten_full_moment_stress_pressure_functional_join_verified',family,hashes)
    bindings=bridge.identity.ast_assignments('current_limit_native_pulse_bridge',
        'CurrentLimitNativePulseBridge','__init__',{
            'self.native_pulse':'CompliantPulseMixedC4()',
            'self.pulse':'CompliantPulseMixedC4()',
            'self.flatten':'CurrentFlattenMixedC4(self.pulse,self.family,self.source,self.hashes,cells)',
            'self.flatten_binding':'current_terminal_source_binding(self)',
            'self.current_limit_Rp_frame':"power['exact_power_to_Rp']['pulse_input_frame']"})
    assert all(report['live_native_defining_object_graph'].values())
    flat=report['actual_same_object_flatten_source_binding']
    assert flat['passed'] and all(flat['current_defining_object_graph'].values())
    assert not flat['source_caps_used_as_defining_field_values'] and not flat['interval_overlap_used_as_join_proof']
    assert flat['zero_linear_moments_are_from_empty_future_supports']
    registry=report['original_native_route_registry']
    assert set(registry)==set(bridge.PULSE_CHARTS)
    for chart,(method,domain) in bridge.EXPECTED.items():
        assert registry[chart]['method']==method and registry[chart]['domain']==domain
    rp=report['actual_native_Rp_handshake']['source_packet']
    rv=report['actual_Rv_terminal_histories'];f=report['actual_flatten_inlet_handshake']['source_packet']
    assert rp['allfive_partial_primitives_callable'] and rp['incoming_and_omitted_tails_retained']
    assert f['original_pressure_datum_and_complete_positive_future_energy_preserved']
    assert f['actual_terminal_zero_linear_and_radial_histories_inherited']
    P0=rp['pressure']['P0_over_Pstar_squared']
    assert P0==rv['P0_over_Pstar_squared']==f['pressure']['P0_over_Pstar_squared']
    assert len(rows(P0))==6
    assert all(is_zero(rows(P0)[n]) for n in (1,3,5))
    # Actual selected pulse cancellation, not a reset of the entrance rows.
    for key in ('Mz_over_R_Utheta','Mtheta_z_over_sqrt2_R_3half_Utheta_squared'):
        assert all(is_zero(row) for row in rows(rv[key])),key
    assert not is_zero(rows(rp['Mz_over_R_Utheta'])[1])
    assert not is_zero(rows(rp['Mtheta_z_over_sqrt2_R_3half_Utheta_squared'])[1])
    assert not is_zero(rows(P0)[0]) and not is_zero(rows(rv['Mp_over_Pstar_squared'])[0])
    c=MPIntervalContext();c.dps=160
    def read(row):
        return c.make_mpf(tuple(tuple(x) for x in (row['lower_exact_mpf_tuple'],row['upper_exact_mpf_tuple'])))
    for packet in (rp['pressure'],rv,f['pressure']):
        p0,mp,p=(rows(packet[k]) for k in ('P0_over_Pstar_squared','Mp_over_Pstar_squared','P_over_Pstar_squared'))
        assert len(p0)==len(mp)==len(p)==6
        for a,b,total in zip(p0,mp,p):
            lo,hi=endpoints(read(total)-read(a)-read(b));assert lo<=0<=hi
    # This pressure check verifies arithmetic of directed returns only. Its
    # analytic source identity is the separate exact defining-function proof.
    for key in bridge.OPEN:assert report[key] is False,key
    assert report['actual_native_O4_constructor_consumes_current_limit_frame']
    assert report['actual_same_pulse_terminal_flatten_source_installed']
    assert report['current_C1_limit_not_promoted_to_global_C4']
    hashes[bridge.NAME]=bridge.sha(bridge.NAME);hashes[Path(__file__).name]=bridge.sha(Path(__file__).name)
    receipt=dict(**{bridge.GATE:True},all_passed=True,source_family=family,
        actual_native_constructor_and_flatten_AST_binding=bindings,
        actual_Rp_Rv_flatten_handshakes_checked=True,
        unchanged_independent_P0_six_rows_and_even_parity=True,
        actual_nonzero_Rp_moments_retained_and_Rv_support_zeros_checked=True,
        seven_same_object_downstream_source_routes_checked=True,
        global_C4_numerics_exterior_heat_and_recursion_not_claimed=True,
        input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (bridge.HERE/bridge.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    print('PASS_CURRENT_LIMIT_NATIVE_PULSE_BRIDGE seven source routes / live Rp-Rv-flatten handshakes',flush=True)
    return receipt


if __name__=='__main__':run()
