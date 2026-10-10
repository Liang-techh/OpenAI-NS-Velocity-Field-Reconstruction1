"""Complete current Rp inlet -> copied actual selected pulse and flatten.

The original algorithms and current repair/future branch remain unchanged.
Fresh source-bound interval constants are installed in all active caches.
This is a restricted native source runtime, not a global point field.
"""
import copy
import json
from pathlib import Path
import time
from types import SimpleNamespace
import lei_ren_part1_paper_compliant_current_original_Rp_interval_inlet as inlet
import lei_ren_part1_paper_compliant_current_selected_energy_source as selected
import lei_ren_part1_paper_compliant_current_limit_heat_pressure_bridge as heat
from lei_ren_part1_paper_compliant_current_native_pulse_source_dispatcher import (
    CurrentNativePulseSourceDispatcher,PULSE_CHARTS,EXPECTED,ROUTES)
from lei_ren_part1_paper_compliant_current_pulse_flatten_source import (
    CurrentFlattenMixedC4,current_terminal_source_binding)
from lei_ren_part1_paper_compliant_pulse_high_jets import _SelectedSource
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import FlatPulseDerivatives
from lei_ren_part1_paper_compliant_current_selected_energy_source_check import source_proof
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE,PREFIX,sha=inlet.HERE,inlet.PREFIX,inlet.sha
NAME=PREFIX+'current_original_Rp_selected_pulse.json'
RECEIPT=PREFIX+'current_original_Rp_selected_pulse_check.json'
GATES=('selected_native_pulse_constructor_consumes_current_C3_frame',
       'current_original_Rp_interval_constants_consumed_by_selected_pulse',
       'current_original_Rp_same_selected_pulse_flatten_owner_installed')
PENDING=('native_interval_inlet_callback_installed','current_numeric_point_field_oracle_installed',
    'current_original_Rp_selected_pulse_installed_in_all_postpulse_owners')+inlet.exact.OPEN


def caches(owner):
    return (owner.amplitude.cache,owner.energy4.cache,owner.axial4.cache,
        owner.fifth.angular_cache,owner.fifth.energy_cache,owner.fifth.cache,owner.pulse.data_cache,
        owner.fifth.angular4.cache)


class CurrentOriginalRpSelectedPulse(CurrentNativePulseSourceDispatcher):
    @source_precision
    def __init__(self,interval=None,seed=None,require_checked=True,cells=128):
        self.inlet=interval if interval is not None else inlet.CurrentOriginalRpIntervalInlet()
        if not self.inlet.acceptance_loaded:raise ValueError('Accepted current interval inlet required')
        if seed is None:
            prefix=heat.CurrentLimitHeatPressureBridge()
            seed=selected.CurrentSelectedEnergySource(pressure=prefix.pressure)
        if type(seed) is not selected.CurrentSelectedEnergySource or not seed.acceptance_loaded:
            raise ValueError('Accepted typed current selected branch required')
        seed.assert_graph();self.seed=seed
        self.family_record=self.inlet.family_record
        self.family,self.source,self.datum_sha=seed.family,seed.source,seed.datum_sha
        if self.family_record!=dict(zip(('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256'),
                (self.family,self.source,self.datum_sha))):raise ValueError('Current inlet/selected family differs')
        self.hashes=dict(self.inlet.hashes);inlet.add_hashes(self.hashes,seed.hashes)
        for name in (inlet.NAME,inlet.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.current_selection_source_proof=source_proof(seed)
        self.ctx=c=seed.future.ctx;self.future=seed.future;self.pressure=seed.pressure;self.exact=seed.exact
        self.original_cache_snapshots=tuple((id(cache),tuple(cache.keys())) for cache in caches(seed))
        self.amplitude=copy.copy(seed.amplitude);self.amplitude.cache={}
        self.energy4=copy.copy(seed.energy4);self.energy4.cache={}
        self.angular4=copy.copy(seed.fifth.angular4);self.angular4.cache={}
        self.energy4.angular=self.angular4
        self.axial4=copy.copy(seed.axial4);self.axial4.energy=self.energy4;self.axial4.base=self.amplitude;self.axial4.cache={}
        self.fifth=copy.copy(seed.fifth);self.fifth.fourth=self.axial4
        self.fifth.angular4=self.angular4
        self.fifth.angular_cache={};self.fifth.energy_cache={};self.fifth.cache={}
        self.pulse=copy.copy(seed.pulse);self.pulse.fifth=self.fifth;self.pulse.data_cache={}
        for owner in (self.amplitude,self.energy4,self.axial4,self.fifth,self.pulse):owner.ctx=c
        self.angular4.ctx=c
        self.pulse.flat=FlatPulseDerivatives(c)
        if (self.pulse.flat.family,self.pulse.flat.source)!=(self.family,self.source):
            raise ValueError('Current beta provider source differs')
        inlet.add_hashes(self.hashes,self.pulse.flat.hashes)
        # The seed's pulse facade retains an older context while its current
        # energy/selection uses the future context. Rebox scalar enclosures;
        # copied Taylor consumers now agree on the actual branch context.
        for key,value in vars(self.pulse).copy().items():
            if hasattr(value,'_mpi_'):setattr(self.pulse,key,c.mpf(inlet.endpoints(value)))
        self.pulse.full_end_energy_weights=[c.mpf(inlet.endpoints(value)) for value in self.pulse.full_end_energy_weights]
        # Rebox full directed enclosures into this actual branch's context.
        self.constants={key:(c.mpf(inlet.endpoints(value)) if hasattr(value,'_mpi_') else value)
            for key,value in self.inlet.constants.items()}
        self.axial4.constants=self.constants
        self.pulse.high=SimpleNamespace(ctx=c,base=self.amplitude,constants=self.constants,
            select=self.fifth.select,energy=SimpleNamespace(future=self.fifth.future))
        self.pulse.selection=_SelectedSource(self.pulse.high)
        self.pulse.pulse=self.amplitude.pulse
        self.pulse.inlet_H=c.mpf(inlet.endpoints(self.inlet.inlet_H))
        self.pulse.inlet_P=c.mpf(inlet.endpoints(self.inlet.inlet_P))
        self.pulse.Xp=self.pulse.inlet_H/self.constants['U']
        # P0 is the actual current selected/future datum, not the saved inlet object.
        self.datum=self.pulse.selection.future.angular.initial.datum
        if (self.datum.definition!=self.inlet.datum.definition or self.datum.source_sha!=self.source
                or self.datum.datum_sha!=self.datum_sha):raise ValueError('Current independent P0 source differs')
        self.current_Rp_physical_log_radius=self.inlet.exact.frame.functions['original_logRp']
        self.native_pulse=self.pulse;self.dispatch=self
        self.chain_routes={}
        for chart,(method,domain) in EXPECTED.items():
            route=ROUTES[chart]
            if route[:3]!=('pulse_mixed_C4','CompliantPulseMixedC4',method) or route[4]!=domain:
                raise ValueError('Original selected pulse route changed: '+chart)
            self.chain_routes[chart]=dict(provider=PREFIX+'pulse_mixed_C4.CompliantPulseMixedC4',
                method=method,coverage_coordinate=route[3],domain=domain,acceptance_receipt=RECEIPT)
        self.flatten=CurrentFlattenMixedC4(self.pulse,self.family,self.source,self.hashes,cells)
        self.flatten_binding=current_terminal_source_binding(self)
        self.Rp_acceptance_loaded=True;self.chain_acceptance_loaded=False;self.acceptance_loaded=False
        self.assert_graph()
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or any(receipt[key] for key in PENDING):
                raise ValueError('Checked current selected consumer receipt required')
            if receipt['source_family']!=self.family_record:raise ValueError('Selected consumer family differs')
            inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=self.chain_acceptance_loaded=True

    def assert_graph(self):
        self.seed.assert_graph()
        graph=dict(fresh_shared_incoming_constant_view=self.constants is self.axial4.constants is
                self.fifth.fourth.constants is self.pulse.high.constants is self.flatten.inlet.constants,
            old_incoming_constant_dictionary_unmutated=self.seed.axial4.constants is self.seed.pulse.high.constants
                and self.constants is not self.seed.axial4.constants,
            all_active_pulse_inputs_current=self.pulse.inlet_H is self.flatten.inlet.inlet_H and
                self.pulse.inlet_P is self.flatten.inlet.inlet_P and self.pulse.Xp is self.flatten.inlet.Xp,
            same_current_C1_C4_C5_future=self.future is self.amplitude.future is self.energy4.base is
                self.fifth.fourth.energy.base,
            same_current_angular_and_unique_repair=self.fifth.angular4 is self.energy4.angular is self.angular4
                and self.angular4.repair is self.future.repair is self.seed.exact.repair,
            copied_current_angular_cache_owner=self.angular4 is not self.seed.fifth.angular4,
            shared_future_and_angle_are_cache_free=not any('cache' in key for obj in
                (self.future,self.future.angular,self.future.repair,self.future.heat) for key in vars(obj)),
            same_current_actual_P0_object=self.flatten.inlet.datum is self.datum is
                self.pulse.selection.future.angular.initial.datum,
            actual_selection_callbacks_bound=self.pulse.high.select.__self__ is self.fifth and
                self.pulse.high.energy.future.__self__ is self.fifth and self.pulse.fifth is self.fifth
                and self.pulse.selection.high is self.pulse.high,
            active_C4_C5_owner_links_current=self.fifth.fourth is self.axial4 and
                self.axial4.base is self.pulse.high.base is self.amplitude and self.axial4.energy is self.energy4,
            fixed_native_kernel_owner_preserved=self.pulse.pulse is self.amplitude.pulse is self.seed.pulse.pulse,
            same_source_context=self.pulse.ctx is self.axial4.ctx is self.energy4.ctx is self.amplitude.ctx is self.ctx,
            current_beta_provider_context=self.pulse.flat.ctx is self.angular4.ctx is self.ctx
                and self.pulse.flat is not self.seed.pulse.flat,
            copied_mutable_caches_separate=not any(id(cache) in {entry[0] for entry in self.original_cache_snapshots}
                for cache in caches(self)),
            old_selected_cache_contents_unmutated=self.original_cache_snapshots==tuple(
                (id(cache),tuple(cache.keys())) for cache in caches(self.seed)),
            copied_native_pulse_used_for_all_routes=all(self.provider(chart) is self.pulse for chart in PULSE_CHARTS),
            same_pulse_used_by_flatten=self.flatten.pulse is self.pulse,
            exact_current_constant_and_interval_receipts_consumed=self.inlet.acceptance_loaded and self.inlet.exact.acceptance_loaded)
        if not all(graph.values()):raise ValueError('Current selected consumer graph differs: '+str(graph))
        return graph

    def provider(self,chart):
        if chart in PULSE_CHARTS:return self.pulse
        if chart=='flatten':return self.flatten
        raise ValueError('Only the selected pulse and flatten are owned here')

    def _packet(self,*args,**kwargs):
        packet=super()._packet(*args,**kwargs)
        packet['acceptance_receipt']=RECEIPT
        packet.update(**dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(PENDING,False))
        return packet

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        self.assert_graph()
        if chart in PULSE_CHARTS:return super().evaluate(chart,Z,coordinate)
        if chart!='flatten':raise ValueError('Unknown current selected source chart')
        t=self.ctx.mpf(coordinate);lo,hi=inlet.endpoints(t)
        if lo<0 or hi>100:raise ValueError('Actual flatten domain [0,100] required')
        return dict(chart=chart,source_packet=self.flatten.flatten(Z,t),acceptance_receipt=RECEIPT,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(PENDING,False))


@source_precision
def run(interval=None,seed=None):
    began=time.monotonic();owner=CurrentOriginalRpSelectedPulse(interval,seed,require_checked=False)
    Z='.371';views={}
    for name,chart,coordinate in (('Rp_entrance','pulse_entrance',0),('active_end_bump','pulse_end',-3),
            ('Rv_terminal','pulse_end',0),('flatten_inlet','flatten',0)):
        views[name]=owner.evaluate(chart,Z,coordinate)
        print('Actual current Rp selected consumer:',name,flush=True)
    result=dict(source_family=owner.family_record,candidate_current_selected_consumer_constructed=True,
        actual_same_object_graph=owner.assert_graph(),actual_current_interval_constants=owner.constants,
        current_selection_source_proof=owner.current_selection_source_proof,
        actual_same_selected_pulse_flatten_binding=owner.flatten_binding,
        original_six_route_registry=owner.chain_routes,
        current_physical_logRp=dict(defining_function_node=owner.current_Rp_physical_log_radius.node,
            accepted_defining_graph=inlet.exact.frame.NAME,accepted_radius_identity_receipt=inlet.exact.frame.RECEIPT),
        fresh_actual_source_handshakes=views,**dict.fromkeys(GATES+PENDING,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_text(json.dumps(inlet.encode(inlet.pack(result)),indent=2)+'\n',encoding='utf8',newline='\n')
    return owner


if __name__=='__main__':run()
