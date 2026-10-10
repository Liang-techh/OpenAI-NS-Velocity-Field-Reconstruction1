"""Repaired-limit pressure branch -> selected future -> actual pulse/flatten.

Reuse the checked selected-energy adapter with explicit pressure injection.
The copied restricted pulse owns the replayed C1/C4/C5 energy and selection;
the original admitted inlet/pressure graph is never mutated.
"""
import json
import time
from pathlib import Path
import lei_ren_part1_paper_compliant_current_limit_heat_pressure_bridge as source
import lei_ren_part1_paper_compliant_current_selected_energy_source as selected
from lei_ren_part1_paper_compliant_current_selected_energy_source_check import source_proof
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE,PREFIX,sha,require=source.HERE,source.PREFIX,source.sha,source.require
NAME=PREFIX+'current_limit_selected_pulse_flatten.json'
RECEIPT=PREFIX+'current_limit_selected_pulse_flatten_check.json'
GATE='current_repaired_limit_selected_future_pulse_flatten_source_installed'
OPEN=source.OPEN+('selected_repaired_pulse_installed_in_all_downstream_physical_owners',
    'selected_repaired_flatten_installed_in_all_downstream_source_owners')


class CurrentLimitSelectedPulseFlatten(source.inlet.CurrentNativePulseSourceDispatcher):
    @source_precision
    def __init__(self,field=None,require_checked=True,cells=128):
        self.source_field=field if field is not None else source.CurrentLimitHeatPressureBridge()
        require(self.source_field.acceptance_loaded,'Checked repaired-limit heat/pressure source required')
        self.family_record=self.source_field.family_record
        self.family=self.source_field.family;self.source=self.source_field.source;self.datum_sha=self.source_field.datum_sha
        self.hashes=dict(self.source_field.hashes)
        old=source.inlet.checked('current_selected_energy_source',selected.GATES[0],self.family_record,self.hashes)
        require(all(old[key] for key in selected.GATES),'All original selected-source gates required')
        self.selected=selected.CurrentSelectedEnergySource(pressure=self.source_field.pressure)
        self.proof=source_proof(self.selected)
        require(self.proof==old['current_future_selected_source_proof'],
            'Original selected future/implicit amplitude proof differs')
        self.native_pulse=self.pulse=self.selected.pulse
        self.hashes.update(self.selected.hashes)
        self.flatten=source.inlet.CurrentFlattenMixedC4(self.pulse,self.family,self.source,self.hashes,cells)
        self.dispatch=self
        self.flatten_binding=source.inlet.current_terminal_source_binding(self)
        self.chain_routes=dict(self.source_field.inlet.chain_routes)
        self.Rp_acceptance_loaded=True;self.chain_acceptance_loaded=False
        self.graph=self.object_graph();require(all(self.graph.values()),'Current selected pulse/flatten source graph differs')
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.acceptance_loaded=False
        if require_checked:
            source.inlet.checked('current_limit_selected_pulse_flatten',GATE,self.family_record,self.hashes)
            self.acceptance_loaded=self.chain_acceptance_loaded=True

    def object_graph(self):
        old=self.source_field.inlet
        return dict(same_actual_pressure_source=self.selected.pressure is self.source_field.pressure,
            same_actual_exact_repair=self.selected.exact is self.source_field.exact,
            same_unique_repair=self.selected.future.repair is self.source_field.exact.repair,
            all_selected_owners_share_branch=all(self.selected.assert_graph().values()),
            new_restricted_pulse_used=self.pulse is self.selected.pulse and self.pulse is not old.pulse,
            new_flatten_uses_selected_pulse=self.flatten.pulse is self.pulse,
            flatten_uses_current_selected_future=self.flatten.pulse.selection.high is self.pulse.high,
            same_current_C5_callback=self.pulse.high.energy.future.__self__ is self.selected.fifth,
            original_pulse_owner_unmutated=old.pulse is self.source_field.exact.pulse,
            original_flatten_owner_unmutated=old.flatten.pulse is old.pulse,
            same_independent_analytic_datum=self.flatten.inlet.datum.source_sha==self.source and
                self.flatten.inlet.datum.datum_sha==self.datum_sha,
            actual_terminal_adapter_bound=self.flatten_binding['passed'])

    def provider(self,chart):
        if chart in source.inlet.PULSE_CHARTS:return self.pulse
        if chart=='flatten':return self.flatten
        raise ValueError('Selected repaired source owns only pulse and flatten routes')

    def _packet(self,*args,**kwargs):
        out=super()._packet(*args,**kwargs)
        out['acceptance_receipt']=RECEIPT
        out.update(**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        return out

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        if chart in source.inlet.PULSE_CHARTS:return super().evaluate(chart,Z,coordinate)
        if chart!='flatten':raise ValueError('Unknown selected repaired pulse/flatten route')
        c=self.flatten.ctx;t=c.mpf(coordinate);lo,hi=source.inlet.endpoints(t)
        if lo<0 or hi>100:raise ValueError('Original flatten domain [0,100] required')
        return dict(chart=chart,source_packet=self.flatten.flatten(Z,t),acceptance_receipt=RECEIPT,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False),
            output_kind='selected repaired pulse/flatten source enclosures; no numeric point selection')


def run():
    began=time.monotonic();field=CurrentLimitSelectedPulseFlatten(require_checked=False)
    result=dict(source_family=field.family_record,**{GATE:False},candidate_selected_source_constructed=True,
        actual_same_object_graph=field.graph,selected_owner_graph=field.selected.assert_graph(),
        reused_original_selected_source_proof=field.proof,
        actual_flatten_terminal_source_binding=field.flatten_binding,
        original_route_registry=field.chain_routes,
        exact_branch_and_independent_pressure_owned_by_current_bridge=True,
        replayed_complete_future_installed_in_restricted_pulse=True,
        replayed_complete_future_installed_in_new_flatten=True,
        original_admitted_graph_unmutated=True,**dict.fromkeys(OPEN,False))
    Z='.371'
    result['fresh_current_future_and_selected_C5']=encode(pack(field.selected.evaluate(Z)))
    print('Current limit selected future and amplitude: fresh Z',flush=True)
    result['actual_Rv_terminal']=encode(pack(field.evaluate('pulse_end',Z,0)))
    result['actual_flatten_inlet']=encode(pack(field.evaluate('flatten',Z,0)))
    print('Current limit selected Rv terminal and flatten inlet',flush=True)
    result['all_passed']=True;result['input_hashes']=field.hashes
    result['execution_seconds']=time.monotonic()-began
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
