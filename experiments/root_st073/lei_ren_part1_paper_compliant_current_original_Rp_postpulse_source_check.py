"""Independent current postpulse caller, history and seam diagnostics.

Functional identities come from pinned original algorithms. Fresh interval
overlaps are diagnostics, not proofs of global mixed regularity or closure.
"""
import copy
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_original_Rp_postpulse_source as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

PAIRS=(('flatten_exit','power_inlet'),('power_exit','angular_inlet'),
    ('angular_exit','steep_entry_inlet'),('steep_entry_exit','steep_power_inlet'),
    ('steep_power_exit','steep_exit_inlet'),('steep_exit_exit','waiting_inlet'),
    ('waiting_exit','collar_inlet'),('collar_exit','exterior_inlet'))


def overlaps(a,b,label):
    assert a.order==b.order==5,label
    for n in range(6):
        lo,hi=current.selected.inlet.endpoints(a[n]);bl,bh=current.selected.inlet.endpoints(b[n])
        assert max(lo,bl)<=min(hi,bh),(label,n)
    return 6


def packet_pressure(packet):
    return packet['pressure']['P_over_Pstar_squared'] if 'pressure' in packet else packet['pressure_over_Pstar_squared_Taylor']


def reject_stale_links(owner):
    rejected=[]
    for label in ('old_pulse','old_flatten','old_angular_owner','stale_heat_future'):
        bad=copy.copy(owner);bad.outer=copy.copy(owner.outer);bad.steep=copy.copy(owner.steep);bad.heat=copy.copy(owner.heat)
        # Keep unrelated graph links correct so the intended stale link is
        # independently rejected by the guard that reads it.
        bad.steep.outer=bad.outer;bad.heat.outer=bad.outer;bad.heat.steep=bad.steep
        bad.power=copy.copy(owner.power);bad.power.outer=bad.outer
        bad.steep_source=copy.copy(owner.steep_source);bad.steep_source.steep=bad.steep;bad.steep_source.before=bad.power
        bad.heat_source=copy.copy(owner.heat_source);bad.heat_source.heat=bad.heat;bad.heat_source.before=bad.steep_source
        bad.companion=copy.copy(owner.companion);bad.companion.heat=bad.heat;bad.companion.source_owner=bad.heat_source
        if label=='old_pulse':bad.outer.pulse=owner.before.seed.pulse
        elif label=='old_flatten':bad.outer.flatten=owner.before.seed.exact.flatten
        elif label=='old_angular_owner':bad.outer.fifth=owner.before.seed.fifth
        else:bad.heat.future=copy.copy(owner.before.future)
        try:bad.assert_graph()
        except ValueError:rejected.append(label)
        else:raise AssertionError('Stale downstream link accepted: '+label)
    return rejected


@source_precision
def run(owner=None):
    began=time.monotonic();raw=json.loads((current.HERE/current.NAME).read_bytes())
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    owner=owner if owner is not None else current.CurrentOriginalRpPostpulseSource(require_checked=False)
    assert not owner.acceptance_loaded and raw['candidate_current_postpulse_constructed']
    assert not any(raw[key] for key in current.GATES+current.OPEN)
    assert raw['source_family']==owner.family_record
    assert all(owner.assert_graph().values()) and raw['actual_same_object_graph']==owner.assert_graph()
    b=owner.before;p=owner.outer;st=owner.steep;h=owner.heat
    assert owner.power.before is b and p.pulse is b.pulse and p.flatten is b.flatten
    assert p.fifth is b.fifth and p.future is st.future is h.future is b.future
    assert h.repair is b.future.repair is b.angular4.repair is b.seed.exact.repair
    assert h.exact_heat is b.future.heat is b.future.repair.heat
    assert owner.companion.heat is h and owner.companion.source_owner is owner.heat_source
    assert p.flatten.inlet.datum is b.datum is b.future.angular.initial.datum
    assert p.ctx is st.ctx is h.ctx is b.ctx
    assert len(owner.registry)==15
    for chart,(provider,method) in current.METHODS.items():
        assert owner.provider(chart) is getattr(owner,provider)
        assert owner.registry[chart]['method']==method
    for chart in b.chain_routes:assert owner.provider(chart) is b.pulse
    assert owner.provider('flatten') is b.flatten
    ast_bindings=current.downstream.inlet.identity.ast_assignments('current_original_Rp_postpulse_source',
        'CurrentOriginalRpPostpulseSource','__init__',{
            'self.power':'downstream.CurrentLimitPowerSource(self.before)',
            'self.steep_source':'downstream.CurrentLimitSteepSource(self.power)',
            'self.heat_source':'downstream.CurrentLimitHeatSource(self.steep_source)',
            'self.companion':'downstream.CurrentLimitHeatCompanion(self.heat_source)'})
    for actual,reported,formula in (
        (owner.power.functional_proof,'retained_original_power_functional_identities',current.downstream.power.functional_source_identities),
        (owner.steep_source.functional_proof,'retained_original_steep_functional_identities',current.downstream.steep.functional_source_identities),
        (owner.heat_source.functional_proof,'retained_original_heat_functional_identities',current.downstream.heat.functional_source_identities)):
        assert actual==raw[reported]==formula() and all(actual.values())
    assert owner.companion.history_bindings==raw['current_absolute_heat_history_binding']
    assert owner.companion.history_bindings['passed']
    assert owner.power.bindings['current_prefix_and_future_decomposition']['passed']
    assert raw['current_complete_future_binding']==owner.power.bindings['current_prefix_and_future_decomposition']
    assert raw['current_heat_source_binding']==owner.heat_source.bindings
    assert not owner.heat_source.bindings['interval_overlap_used_as_join_proof']
    rejected=reject_stale_links(owner)
    Z=raw['fresh_Z'];views={};rows=0
    for name,chart,coordinate in current.VIEWS:
        # The checker reads the actual method outputs from the current live
        # providers, not the serialized candidate's truth flags.
        packet=owner.evaluate(chart,Z,coordinate)['source_packet'];views[name]=packet
        history=owner.histories(chart,Z,coordinate,packet)
        assert current.selected.inlet.encode(current.selected.inlet.pack(history))==raw['fresh_actual_source_views'][name]['five_source_histories']
        for key in ('Mz_over_R_Utheta','Mtheta_z_over_sqrt2_R_3half_Utheta_squared'):
            assert history[key].order==5
            for n in range(6):assert current.selected.inlet.endpoints(history[key][n])==(0,0)
        assert history['Mtheta_over_sqrt2_R_3half_Utheta'] is packet['angular_Taylor']
        assert history['Mztheta_over_R_Utheta_squared'] is packet['energy_Taylor']
        assert current.selected.inlet.endpoints(packet['energy_Taylor'][0])[0]>0
        assert history['P0_over_Pstar_squared'].order==5
        assert current.selected.inlet.endpoints(history['P0_over_Pstar_squared'][0])!=(0,0)
        overlaps(history['P_over_Pstar_squared'],history['Mp_over_Pstar_squared']+history['P0_over_Pstar_squared'],name+'/P0+Mp')
    for left,right in PAIRS:
        for key in ('theta_over_Ev0_Taylor','angular_Taylor','energy_Taylor'):
            rows+=overlaps(views[left][key],views[right][key],left+'->'+right+'/'+key)
        rows+=overlaps(packet_pressure(views[left]),packet_pressure(views[right]),left+'->'+right+'/pressure')
    assert rows==192
    constants=owner.heat_constants(Z)
    assert current.selected.inlet.encode(current.selected.inlet.pack(constants))==raw['retained_actual_heat_terminal_constants']
    assert constants['zero_constants_not_assumed']
    assert owner.before.assert_graph()['old_selected_cache_contents_unmutated']
    assert raw['no_repair_replay_or_reselection'] and raw['physical_radius_handle_only']
    result=dict(all_passed=True,source_family=owner.family_record,**dict.fromkeys(current.GATES,True),
        actual_independent_current_owner_graph=owner.assert_graph(),actual_constructor_AST_bindings=ast_bindings,
        stale_downstream_owner_links_rejected=rejected,all_fifteen_routes_use_current_selected_owner=True,
        same_current_unique_repair_and_complete_future_retained=True,
        original_analytic_P0_and_five_normalized_source_histories_retained=True,
        actual_eighteen_source_calls_checked=True,fresh_Z=Z,eight_seam_Taylor_diagnostic_rows=rows,
        seam_overlap_not_used_as_function_identity_or_global_contract=True,
        actual_heat_angular_and_pressure_constants_retained_not_declared_zero=True,
        old_selected_cache_contents_unmutated=True,
        scope='Fifteen current native source routes and normalized histories; physical radius, absolute heat closure and global mixed contracts remain open',
        **dict.fromkeys(current.OPEN,False),input_hashes={**owner.hashes,
            current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.selected.inlet.encode(current.selected.inlet.pack(result)),indent=2)+'\n',
        encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_POSTPULSE 15 current routes; 18 calls; 192 seam diagnostics; five source histories',flush=True)
    return result


if __name__=='__main__':run()
