"""Admission of the new injection and fresh selected Rv/flatten handshakes."""
import json
import time
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_limit_selected_pulse_flatten as bridge


def run():
    began=time.monotonic();report=json.loads((bridge.HERE/bridge.NAME).read_bytes())
    assert report['all_passed'] and report['candidate_selected_source_constructed']
    assert report[bridge.GATE] is False
    hashes=dict(report['input_hashes']);family=report['source_family']
    for name,digest in hashes.items():assert bridge.sha(name)==digest,name
    bridge.source.inlet.checked('current_limit_heat_pressure_bridge',bridge.source.GATE,family,hashes)
    old=bridge.source.inlet.checked('current_selected_energy_source',bridge.selected.GATES[0],family,hashes)
    assert all(old[key] for key in bridge.selected.GATES)
    assert report['reused_original_selected_source_proof']==old['current_future_selected_source_proof']
    assert all(report['actual_same_object_graph'].values()) and all(report['selected_owner_graph'].values())
    assert report['actual_flatten_terminal_source_binding']['passed']
    bindings=bridge.source.inlet.identity.ast_assignments('current_limit_selected_pulse_flatten',
        'CurrentLimitSelectedPulseFlatten','__init__',{
            'self.selected':'selected.CurrentSelectedEnergySource(pressure=self.source_field.pressure)',
            'self.pulse':'self.selected.pulse','self.native_pulse':'self.selected.pulse',
            'self.flatten':'source.inlet.CurrentFlattenMixedC4(self.pulse,self.family,self.source,self.hashes,cells)',
            'self.flatten_binding':'source.inlet.current_terminal_source_binding(self)',
            'self.proof':'source_proof(self.selected)'})
    assert set(report['original_route_registry'])==set(bridge.source.inlet.PULSE_CHARTS)
    for chart,(method,domain) in bridge.source.inlet.EXPECTED.items():
        route=report['original_route_registry'][chart]
        assert route['method']==method and route['domain']==domain
    c=MPIntervalContext();c.dps=160
    def read(row):
        return c.make_mpf(tuple(tuple(row[k]) for k in ('lower_exact_mpf_tuple','upper_exact_mpf_tuple')))
    def rows(value):return [read(row) for row in value['coefficients']]
    def equal_enclosed(a,b):
        lo,hi=bridge.source.inlet.endpoints(a-b);assert lo<=0<=hi
    def zero(row):return row['lower_exact_mpf_tuple'][1]==row['upper_exact_mpf_tuple'][1]==0
    fresh=report['fresh_current_future_and_selected_C5']
    E4=rows(fresh['current_future_C4']['complete_future_energy_Taylor'])
    E5=rows(fresh['current_future_C5']['complete_future_energy_Taylor'])
    assert len(E4)==5 and len(E5)==6 and bridge.source.inlet.endpoints(E5[0])[0]>0
    for n in range(5):assert E4[n]._mpi_==E5[n]._mpi_
    fourth=fresh['current_selected_C4'];fifth=fresh['current_selected_C5']
    for a,b in zip([fourth['selected_ap_Taylor']]+fourth['selected_scaled_end_coefficient_Taylor'],
                   [fifth['selected_ap_Taylor']]+fifth['selected_scaled_end_coefficient_Taylor']):
        aa,bb=rows(a),rows(b);assert len(aa)==5 and len(bb)==6
        for n in range(5):assert aa[n]._mpi_==bb[n]._mpi_
    ap=rows(fifth['selected_ap_Taylor'])
    lo,hi=bridge.source.inlet.endpoints(ap[0]);assert c.mpf('.9')<lo<=hi<c.mpf('1.2')
    assert bridge.source.inlet.endpoints(read(fifth['positive_root_derivative_denominator']))[0]>0
    rv=report['actual_Rv_terminal']['source_packet']
    flat=report['actual_flatten_inlet']['source_packet']
    assert rv['pulse_all_mixed_derivatives_total_order_le4_available']
    for key in ('Uz_over_Utheta','Mz_over_R_Utheta',
                'Mtheta_z_over_sqrt2_R_3half_Utheta_squared','Ur_over_sqrt_R_over_2_Utheta'):
        assert all(zero(row) for row in rv[key]['coefficients']),key
    terminal=rows(rv['Mztheta_over_R_Utheta_squared']);inlet=rows(flat['energy_Taylor'])
    for n in range(6):
        equal_enclosed(terminal[n],E5[n]/2)
        equal_enclosed(inlet[n],E5[n]/2)
    P0=rv['pressure']['P0_over_Pstar_squared']
    assert P0==flat['pressure']['P0_over_Pstar_squared']
    assert len(P0['coefficients'])==6 and not zero(P0['coefficients'][0])
    for key in bridge.OPEN:assert report[key] is False
    for key in ('exact_branch_and_independent_pressure_owned_by_current_bridge',
                'replayed_complete_future_installed_in_restricted_pulse',
                'replayed_complete_future_installed_in_new_flatten','original_admitted_graph_unmutated'):
        assert report[key]
    hashes[bridge.NAME]=bridge.sha(bridge.NAME);hashes[Path(__file__).name]=bridge.sha(Path(__file__).name)
    receipt=dict(source_family=family,all_passed=True,**{bridge.GATE:True},
        actual_selected_pressure_and_flatten_constructor_AST_bound=bindings,
        strict_original_selected_equation_function_proof_reused=True,
        actual_new_C5_future_callback_and_original_graph_retention_checked=True,
        fresh_C5_future_and_selection_preserve_C4_prefix=True,
        actual_Rv_and_flatten_use_same_complete_future_half=True,
        actual_terminal_meridional_support_zeros_and_nonzero_P0_retained=True,
        all_physical_and_later_source_owner_installation_not_claimed=True,
        **dict.fromkeys(bridge.OPEN,False),input_hashes=hashes,
        execution_seconds=time.monotonic()-began)
    (bridge.HERE/bridge.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    print('PASS_CURRENT_LIMIT_SELECTED_PULSE_FLATTEN fresh C5 selection, Rv full future/2, same flatten',flush=True)
    return receipt


if __name__=='__main__':run()
