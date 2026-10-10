"""Certified repaired-limit Rp functions -> actual native pulse and flatten.

Only downstream source owners are constructed. The repaired inlet is
identified by exact defining functions; directed bounds stay bounds.
"""
import gzip
import json
import time
from pathlib import Path
from types import SimpleNamespace
import lei_ren_part1_paper_compliant_current_limit_Rp_native_identity as identity
from lei_ren_part1_paper_compliant_current_native_pulse_source_dispatcher import (
    CurrentNativePulseSourceDispatcher, PULSE_CHARTS, EXPECTED, ROUTES)
from lei_ren_part1_paper_compliant_current_pulse_flatten_source import (
    CurrentFlattenMixedC4, current_terminal_source_binding, THEOREM)
from lei_ren_part1_paper_compliant_pulse_mixed_C4 import CompliantPulseMixedC4
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE,PREFIX,sha,require=identity.HERE,identity.PREFIX,identity.sha,identity.require
NAME=PREFIX+'current_limit_native_pulse_bridge.json'
RECEIPT=PREFIX+'current_limit_native_pulse_bridge_check.json'
GATE='current_repaired_limit_native_pulse_and_flatten_source_owners_installed'
OPEN=('actual_numeric_point_values_installed','patched_Rh_functional_join_proved',
    'physical_original_exterior_five_targets_closed','actual_temporal_scale_recursion_installed',
    'uniform_pulse_C4_chart_interface_certificate_available','current_flatten_physical_owner_installed',
    'full_current_core_to_heat_physical_assembly')


def checked(stem,gate,family,hashes):
    name=PREFIX+stem+'_check.json';row=json.loads((HERE/name).read_bytes())
    require(row['all_passed'] and row[gate],'Accepted native certificate required: '+name)
    if 'source_family' in row:require(row['source_family']==family,'Different current family: '+name)
    else:
        for key in ('implicit_source_sha256','datum_enclosure_sha256','actual_five_defect_family_sha256'):
            if key in row:require(row[key]==family[key],'Different native source: '+name+'/'+key)
    for path,digest in row['input_hashes'].items():
        require(sha(path)==digest,'Changed native dependency: '+path)
        require(path not in hashes or hashes[path]==digest,'Native dependency conflict: '+path)
        hashes[path]=digest
    hashes[name]=sha(name)
    return row


class CurrentLimitNativePulseBridge(CurrentNativePulseSourceDispatcher):
    """Strict original six pulse routes and one same-object flatten owner."""
    @source_precision
    def __init__(self,require_checked=True,cells=128):
        proof=json.loads((HERE/identity.NAME).read_bytes())
        self.identity=proof;self.identity_family=self.family_record=proof['source_family'];self.hashes={}
        checked('current_limit_Rp_native_identity',identity.GATE,self.family_record,self.hashes)
        require(proof[identity.GATE] and proof['actual_Rp_native_frame_source_function_identified'],
            'Exact repaired-limit inlet function identity required')
        self.family=self.family_record['actual_five_defect_family_sha256']
        self.source=self.family_record['implicit_source_sha256'];self.datum_sha=self.family_record['datum_enclosure_sha256']
        power=json.loads(gzip.decompress((HERE/identity.current.NAME).read_bytes()))
        self.current_limit_Rp_frame=power['exact_power_to_Rp']['pulse_input_frame']
        self.current_limit_Rp_radius=power['exact_power_to_Rp']['Rp']
        self.current_limit_relative_zero=power['exact_power_to_Rp']['relative_zero_propagation']
        checked('actual_Rp_source_join','current_Rp_external_pulse_join_certified',self.family_record,self.hashes)
        checked('current_native_pulse_source_dispatcher','current_native_pulse_source_ownership_certified',self.family_record,self.hashes)
        checked('current_pulse_flatten_source','current_pulse_terminal_flatten_source_ownership_certified',self.family_record,self.hashes)
        checked('pulse_end_flatten_join','pulse_end_flatten_full_moment_stress_pressure_functional_join_verified',self.family_record,self.hashes)
        self.native_pulse=self.pulse=CompliantPulseMixedC4()
        initial=self.pulse.pulse.initial;future=self.pulse.selection.future.angular.initial
        ps=proof['current_P0_source_binding']
        require(initial.family==future.family==self.family,'Live native family differs')
        for datum in (initial.datum,future.datum):
            require(datum.source_sha==self.source and datum.datum_sha==self.datum_sha
                and datum.definition==ps['analytic_definition'],'Live native analytic pressure datum differs')
        pre=SimpleNamespace(initial=initial,buffer=self.pulse.pulse.buffer,family=self.family,
            datum=initial.datum,params=initial.params)
        self.live_graph=identity.current.native.native_source_graph(pre,self.pulse)
        self.chain_routes={}
        for chart,(method,domain) in EXPECTED.items():
            route=ROUTES[chart]
            require(route[:3]==('pulse_mixed_C4','CompliantPulseMixedC4',method) and route[4]==domain,
                'Original native route changed: '+chart)
            self.chain_routes[chart]=dict(provider=PREFIX+'pulse_mixed_C4.CompliantPulseMixedC4',
                method=method,coverage_coordinate=route[3],domain=domain,acceptance_receipt=RECEIPT)
        self.Rp_acceptance_loaded=True;self.chain_acceptance_loaded=True
        self.flatten=CurrentFlattenMixedC4(self.pulse,self.family,self.source,self.hashes,cells)
        self.dispatch=self
        self.flatten_binding=current_terminal_source_binding(self)
        for path,digest in self.pulse.hashes.items():
            require(sha(path)==digest,'Changed live native source: '+path)
            require(path not in self.hashes or self.hashes[path]==digest,'Different live source: '+path)
            self.hashes[path]=digest
        for stem in ('current_limit_native_pulse_bridge','current_pulse_flatten_source'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            checked('current_limit_native_pulse_bridge',GATE,self.family_record,self.hashes)
            self.acceptance_loaded=True

    def provider(self,chart):
        if chart in PULSE_CHARTS:return self.pulse
        if chart=='flatten':return self.flatten
        raise ValueError('Repaired-limit bridge owns only native pulse and flatten source routes')

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        if chart in PULSE_CHARTS:packet=super().evaluate(chart,Z,coordinate)
        elif chart=='flatten':
            c=self.flatten.ctx;t=c.mpf(coordinate);lo,hi=endpoints(t)
            if lo<0 or hi>100:raise ValueError('Original flatten domain is [0,100]')
            raw=self.flatten.flatten(Z,t)
            packet=dict(chart=chart,source_packet=raw,source_coordinate_domain='[0,100]',
                source_coverage_coordinate='t=log(R/Rv)',
                output_kind='native source enclosures; no production point selection')
        else:raise ValueError('Unsupported downstream source chart: '+chart)
        packet.update(repaired_current_limit_Rp_source_function_identity_certified=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        return packet

    def manifest(self):
        return dict(**{GATE:self.acceptance_loaded},source_family=self.family_record,
            actual_Rp_native_frame_source_function_identified=True,
            actual_native_O4_constructor_consumes_current_limit_frame=True,
            actual_same_pulse_terminal_flatten_source_installed=True,
            repaired_current_Rp_frame=self.current_limit_Rp_frame,
            repaired_current_Rp_radius=self.current_limit_Rp_radius,
            current_relative_zero_propagation=self.current_limit_relative_zero,
            original_native_route_registry=self.chain_routes,
            live_native_defining_object_graph=self.live_graph,
            actual_same_object_flatten_source_binding=self.flatten_binding,
            same_single_native_pulse_object_for_six_charts=True,
            actual_Rv_terminal_five_history_adapter_available=True,
            current_C1_limit_not_promoted_to_global_C4=True,
            underlying_functions_identified_not_numeric_interval_outputs_equal=True,
            accepted_current_all_N_inner_source_producers_not_replayed=True,
            **dict.fromkeys(OPEN,False),input_hashes=dict(self.hashes))


def run():
    began=time.monotonic();owner=CurrentLimitNativePulseBridge(require_checked=False)
    result=owner.manifest();result[GATE]=True
    # One live handshake at the true inlet. Each return retains independent
    # P0 and all actual incoming moment functions; no midpoint is selected.
    result['actual_native_Rp_handshake']=encode(pack(owner.evaluate('pulse_entrance',0,0)))
    result['actual_Rv_terminal_histories']=encode(pack(owner.flatten.terminal_histories(0)))
    result['actual_flatten_inlet_handshake']=encode(pack(owner.evaluate('flatten',0,0)))
    result['input_hashes']=dict(owner.hashes);result['execution_seconds']=time.monotonic()-began
    result['all_passed']=True
    (HERE/NAME).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('CURRENT_LIMIT_NATIVE_PULSE_BRIDGE actual pulse/flatten owners; inlet and Rv handshakes',flush=True)
    return result


if __name__=='__main__':run()
