"""Live common-unit source covers for the seventeen generic-loop charts.

One genuine original seed supplies all providers. No full tensor graph, saved
source cover fallback, selected midpoint, new cutoff or numerical operator is
introduced. These are coordinate-dependent enclosures of original functions.
"""
import gzip
import importlib
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_generic_left_inlet as inlet
from lei_ren_part1_paper_compliant_current_microswitch_stress_operator import (
    compiled_original_microswitch_with_rows, raw_microswitch_rows)
from lei_ren_part1_paper_compliant_current_switch_power_source_adapter import compiled_original_postpower_with_locals
from lei_ren_part1_paper_compliant_current_reshape_stress_operator import raw_reshape_rows
from lei_ren_part1_paper_compliant_current_restore_stress_operator import raw_restore_rows
from lei_ren_part1_paper_compliant_current_patch_stress_operator import raw_patch_rows
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import raw_pre_velocity_rows
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import copy_jet
from lei_ren_part1_paper_compliant_current_pre_pulse_source_dispatcher import original_history_source_bindings

packets=inlet.packets;HERE,PREFIX,sha=inlet.HERE,inlet.PREFIX,inlet.sha
NAME=PREFIX+'current_native_generic_source_packets.json'
RECEIPT=PREFIX+'current_native_generic_source_packets_check.json'
GATE='current_seventeen_generic_original_native_coordinate_cover_backends_executed'
VIEWS=PREFIX+'current_native_generic_source_packets_views.json.gz'

DOMAINS=dict(inlet.DOMAINS,switch_first=(0,1),switch_second=(1,2),switch_power=(0,1),
    reshape=(0,1),inner_reference=(0,1),axial_restore=(0,1),restore_buffer=(-7,-6),
    actual_patch=(1,'e'),Rh_reference=(-5,0),O2_slope=(0,1),O2_axial=(0,1),O2_buffer=(0,11),
    O3_slope_mu=(0,1),O3_power=(0,1))
GROUPS={
    'micro':('current_microswitch_background_tensor','current_actual_two_microswitch_full_tensors_available'),
    'power':('current_switch_power_background_tensor','current_actual_switch_power_full_tensor_available'),
    'reshape':('current_reshape_background_tensor','current_actual_long_reshape_full_tensor_available'),
    'restore':('current_restore_background_tensor','current_actual_inner_reference_axial_restore_buffer_full_tensors_available'),
    'patch':('current_actual_patch_background_tensor','current_actual_five_moment_patch_full_tensor_available'),
    'pre':('current_O2_background_tensor','current_actual_O2_reference_slope_axial_buffer_full_tensors_available'),
    'O3':('current_O3_transition_background_tensor','current_actual_O3_transition_full_tensor_available')}
PROVIDERS={
    'switch':('actual_switch_mixed_C4','CompliantActualSwitchMixedC4',('evaluate','postpower')),
    'reshape':('actual_long_reshape_mixed_C4','CompliantActualLongReshapeMixedC4',('evaluate',)),
    'restore':('actual_reference_restore_mixed_C4','CompliantActualReferenceRestoreMixedC4',('reference_branch','restoration','postrestore')),
    'patch':('actual_feedback_patch_mixed_C4','CompliantActualFeedbackPatchMixedC4',('evaluate',)),
    'pre':('pre_pulse_mixed_C4','CompliantPrePulseMixedC4',('reference','slope','axial','slope_mu','power','packet'))}


class NativeGenericSourcePackets(inlet.NativeBridgeSourcePackets):
    def __init__(self,bridge):
        super().__init__(bridge)
        self.seed=bridge.owner
        self.switch=bridge.switch
        self.reshape=self.seed.dispatch.provider('reshape')
        self.restore=self.seed.dispatch.provider('inner_reference')
        self.patch=self.seed.dispatch.provider('actual_patch')
        self.pre=self.seed.pre
        self.native_graph=self.assert_native_graph()
        self.receipts={}
        for group,(stem,gate) in GROUPS.items():
            name=PREFIX+stem+'_check.json'
            receipt=inlet.accepted(name,bridge.family,bridge.source,gate)
            if receipt['datum_enclosure_sha256']!=bridge.datum_sha:
                raise ValueError('Same original native chart axis datum required')
            self.service.bind_hashes(receipt['input_hashes'])
            self.service.bind_hashes({name:sha(name)})
            self.receipts[group]=receipt
        for stem,gate in (('current_pre_pulse_source_dispatcher','current_Rh_to_Rp_source_chain_certified'),
                          ('pre_pulse_mixed_C4','all_pre_pulse_mixed4_available')):
            name=PREFIX+stem+'_check.json'
            receipt=inlet.accepted(name,bridge.family,bridge.source,gate)
            self.service.bind_hashes(receipt['input_hashes'])
            self.service.bind_hashes({name:sha(name)})
        self.micro_evaluator,self.micro_proof=compiled_original_microswitch_with_rows()
        proof=self.receipts['micro']['current_actual_microswitch_source_and_tensor_theorem']
        if self.micro_proof!=proof['original_output_only_compiler_theorem']:
            raise ValueError('Accepted original micro source compiler differs')
        self.power_evaluator,self.power_proof=compiled_original_postpower_with_locals()
        proof=self.receipts['power']['current_actual_switch_power_source_and_tensor_theorem'][
            'current_actual_switch_power_source_pressure_and_R110_theorem']
        if self.power_proof!=proof['original_output_only_source_locals_adapter']:
            raise ValueError('Accepted original power locals compiler differs')
        self.history_bindings=original_history_source_bindings()
        for evidence in (self.micro_proof,self.power_proof):
            self.service.bind_hashes(evidence['input_hashes'])
        for stem in ('current_native_generic_source_packets','current_microswitch_stress_operator',
            'current_reshape_stress_operator','current_restore_stress_operator','current_patch_stress_operator',
            'current_pre_pulse_stress_operator','current_pre_pulse_source_dispatcher','pre_pulse_mixed_C4'):
            name=PREFIX+stem+'.py';self.service.bind_hashes({name:sha(name)})

    def assert_native_graph(self):
        from lei_ren_part1_paper_compliant_current_core_physical_assembly import CurrentCorePhysicalAssembly
        if (type(self.seed) is not CurrentCorePhysicalAssembly or not self.seed.core_acceptance_loaded or
                not self.seed.physical_acceptance_loaded or not all(self.seed.current_provider_graph().values())):
            raise ValueError('One accepted original core/heat physical seed required')
        graph=dict(same_bridge_seed=self.bridge.owner is self.seed,
            same_context=self.ctx is self.seed.ctx,
            same_family_datum=(self.seed.family,self.seed.source,self.seed.datum_sha)==
                (self.bridge.family,self.bridge.source,self.bridge.datum_sha),
            same_switch=self.switch is self.seed.dispatch.provider('switch_first') is self.seed.dispatch.provider('switch_second') is self.seed.dispatch.provider('switch_power'),
            same_reshape=self.reshape is self.seed.dispatch.provider('reshape') is self.restore.long_mixed,
            same_restoration=self.restore is self.seed.dispatch.provider('axial_restore') is self.seed.dispatch.provider('restore_buffer') is self.patch.patch.reference_mixed,
            same_patch=self.patch is self.seed.dispatch.anchor,
            same_pre=self.pre is self.seed.dispatch.rh_reference,
            same_pre_parameters=self.pre.params is self.pre.datum.parameters,
            same_pre_seed_parameters=self.pre.params is self.seed.params,
            same_pre_datum=self.pre.datum is self.pre.initial.repair.datum,
            same_pre_source_family=(self.pre.family,self.pre.source,self.pre.datum.datum_sha)==
                (self.bridge.family,self.bridge.source,self.bridge.datum_sha),
            same_native_provider_graph=all(self.seed.dispatch.native.provider_graph_identity().values()))
        for key,(stem,class_name,methods) in PROVIDERS.items():
            provider=getattr(self,key);original=getattr(importlib.import_module(PREFIX+stem),class_name)
            graph[key+'_exact_original_class']=type(provider) is original
            graph[key+'_same_context']=provider.ctx is (self.pre.initial.ctx if key=='pre' else self.ctx)
            for method in methods:
                graph[key+'_'+method+'_unchanged']=getattr(provider,method).__func__ is getattr(original,method)
        if not all(graph.values()):
            raise ValueError('Original native generic provider graph differs: '+', '.join(k for k,v in graph.items() if not v))
        return graph

    @inlet.source_precision
    def original_raw(self,chart,Z,coordinate):
        self.assert_native_graph();c=self.ctx;z=c.mpf(Z);v=c.mpf(coordinate);ep=packets.recovery.endpoints
        if chart not in DOMAINS:
            raise ValueError('One of the17 original generic source charts required')
        lo,hi=DOMAINS[chart];hi=ep(c.exp(1))[1] if hi=='e' else hi
        if (not all(mp.isfinite(x) for x in ep(z)+ep(v)) or ep(z)[0]<-1 or ep(z)[1]>1 or ep(v)[0]<lo or ep(v)[1]>hi):
            raise ValueError('Original native chart domain and Z[-1,1] required')
        if chart in inlet.DOMAINS:
            original=self.evaluator(self.bridge,z,v,chart.replace('bridge_',''))
            raw=self.adapter(c,original)
            p0=original['actual_parent_axial5_packet']['pressure_axis_axial5_coefficients']
            provider=self.bridge;group='bridge'
        elif chart in ('switch_first','switch_second'):
            original=self.micro_evaluator(self.switch,z,v,chart.replace('switch_',''))
            raw=raw_microswitch_rows(c,original)
            p0=original['actual_parent_axial5_packet']['pressure_axis_axial5_coefficients']
            provider=self.switch;group='micro'
        elif chart in ('switch_power','reshape'):
            provider=self.switch if chart=='switch_power' else self.reshape
            original=self.power_evaluator(provider,z,v) if chart=='switch_power' else provider.evaluate(z,v)
            raw=raw_reshape_rows(c,original,c.mpf(ep(provider.invP2)))
            p0=original['actual_inherited_axial5_packet']['pressure_axis_axial5_coefficients']
            group='power' if chart=='switch_power' else 'reshape'
        elif chart in ('inner_reference','axial_restore','restore_buffer'):
            method={'inner_reference':'reference_branch','axial_restore':'restoration','restore_buffer':'postrestore'}[chart]
            original=getattr(self.restore,method)(z,v)
            raw=raw_restore_rows(c,original,c.mpf(ep(self.restore.invP2)))
            p0=original['actual_inherited_axial5_packet']['pressure_axis_axial5_coefficients']
            provider=self.restore;group='restore'
        elif chart=='actual_patch':
            original=self.patch.evaluate(v,z)
            raw=raw_patch_rows(c,original,c.mpf(ep(self.patch.invP2)))
            p0=original['actual_inherited_patch_packet']['original_P0_axial5']
            provider=self.patch;group='patch'
        else:
            if chart=='Rh_reference':original=self.pre.reference(z,v)
            elif chart=='O2_slope':original=self.pre.slope(z,v)
            elif chart=='O2_axial':original=self.pre.axial(z,phase=v)
            elif chart=='O2_buffer':original=self.pre.axial(z,buffer_offset=v)
            elif chart=='O3_slope_mu':original=self.pre.slope_mu(z,v)
            else:original=self.pre.power(z,v)
            histories={key:[copy_jet(c,row) for row in original['actual_normalized_primitive_y_derivative_axial5'][key]]
                for key in packets.recovery.RATES}
            p0=original['original_P0_axial5_coefficients'];datum=packets.jet(c,p0)
            raw=dict(velocity=raw_pre_velocity_rows(c,original),histories=histories,
                absolute_pressure=[histories['p'][0]+datum]+histories['p'][1:])
            provider=self.pre;group='O3' if chart.startswith('O3') else 'pre'
        logs=raw['algebra'].logs if 'algebra' in raw else (c.mpf(0),2*self.seed.logP,c.mpf(0),c.mpf(0))
        return dict(chart=chart,Z=z,coordinate=v,raw=raw,P0=packets.jet(c,p0[:6]),logs=logs,
            native_provider=provider,native_provider_group=group,source_packet=original,
            positive_caps_enclose_original_source_not_define_it=True,
            original_radial_half_shift_and_width_conversion_already_applied=True)

    @inlet.source_precision
    def query(self,chart,Z,coordinate):
        source=self.original_raw(chart,Z,coordinate);raw=source['raw'];c=self.ctx
        algebra=packets.FactoredAlgebra(c,source['logs'],[])
        def rows(values):
            if len(values)!=5:raise ValueError('All five ordinary logR rows required')
            result=tuple(packets.decode_row(algebra,packets.factored_rows_record(value)) for value in values)
            packets.factored_check_jets(c,result)
            return result
        native_v={key:rows(value) for key,value in raw['velocity'].items()}
        native_m={key:rows(value) for key,value in raw['histories'].items()}
        if set(native_v)!=set(('theta','axial','radial')) or set(native_m)!=set(packets.recovery.RATES):
            raise ValueError('All original velocities and five histories required')
        velocity={key:tuple(algebra.shift(row,packets.INVERSE_S) if key in ('axial','radial') else row
            for row in values) for key,values in native_v.items()}
        histories={key:tuple(algebra.shift(row,packets.INVERSE_S) if key in ('m','k') else row
            for row in values) for key,values in native_m.items()}
        pressure=rows(raw['absolute_pressure']);p0=algebra.lift(source['P0'])
        if p0.order!=5:raise ValueError('Separate original normalized P0 axial5 required')
        provenance=dict(mode='injected_whole_native_original_source_cover',chart=chart,
            Z_box=source['Z'],coordinate_box=source['coordinate'],source_family=self.family,
            owner_class=type(source['native_provider']).__name__,owner_module=type(source['native_provider']).__module__,
            native_provider_group=source['native_provider_group'],one_original_seed_graph_asserted=True,
            arbitrary_coordinates_evaluated=True,cache_cover=False,
            derivative_coordinate='ordinary y=log R',selector_derivatives_not_substituted=True,
            physical_radial_half_shift_and_original_width_conversion_applied_once=True,
            common_Pstar_unit_shift_applied_once=True,P0_read_directly_not_subtracted=True,
            pre_source_context_rebased_with_original_copy_jet=source['native_provider_group'] in ('pre','O3'),
            numerical_caps_not_used_as_defining_function_values=True,
            source_factor_resolution_performed=False,full_tensor_or_physical_point_requested=False,
            native_operator_outputs_are_covers_not_selected_point_values=True)
        self.queries.append(provenance)
        return packets.CurrentSourcePacket(chart,dict(self.family),algebra,
            packets.IntervalTaylor.variable(c,source['Z'],5),velocity,histories,pressure,p0,native_v,native_m,provenance)


def run(bridge=None,runtime_evidence=None):
    began=time.monotonic()
    if bridge is None:bridge,runtime_evidence=inlet.native_bridge_owner()
    with inlet.CheckedSourceRuntime():
        backend=NativeGenericSourcePackets(bridge)
        views={}
        points=dict(bridge_first='.1337',bridge_second='1.831',bridge_macro='.537',
            switch_first='.537',switch_second='1.337',switch_power='.537',reshape='.537',
            inner_reference='.537',axial_restore='.537',restore_buffer='-6.337',actual_patch='1.337',
            Rh_reference='-2.337',O2_slope='.537',O2_axial='.1337',O2_buffer='5.337',O3_slope_mu='.537')
        points['O3_power']=backend.ctx.mpf('.537')/backend.pre.params.Tw
        for chart in DOMAINS:
            packet=backend.query(chart,(-1,1),points[chart])
            views[chart]=packet.record()
            print('Live original generic source:',chart,flush=True)
        # Preserve the actual strict positive generic inlet through this same
        # backend. Its original histories/P0 and zero initial-condition scope
        # are inherited from the accepted native inlet helper.
        left=backend.left_inlet();packet=left.pop('original_source_packet')
        views['generic_left_inlet']=dict(left,original_source_packet=packet.record())
    (HERE/VIEWS).write_bytes(gzip.compress((json.dumps(packets.encode(views),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    hashes=dict(backend.service.hashes);hashes[VIEWS]=sha(VIEWS)
    result=dict(source_family=backend.family,**{GATE:True},
        native_generic_chart_count=len(DOMAINS),successful_coordinate_queries=len(backend.queries),
        original_native_provider_graph=backend.native_graph,
        original_native_source_chart_domains=DOMAINS,actual_query_provenance=backend.queries,
        actual_whole_Z_chart_cover_views=VIEWS,
        O3_power_sample_local_offset_encloses_nominal_point537=True,
        source_backend_kind='original live coordinate enclosures, not selected numerical point functions',
        unchanged_original_history_source_bindings=backend.history_bindings,
        unchanged_micro_source_compiler=backend.micro_proof,unchanged_power_source_locals_compiler=backend.power_proof,
        original_common_unit_theorem=backend.unit_proof,runtime_construction_evidence=runtime_evidence,
        source_ancestor_or_full_tensor_graph_rebuilt_for_chart_query=False,
        actual_changed_defect_integral_functions_installed=False,
        actual_repair_control_functions_installed=False,actual_terminal_Z_function_closure_installed=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,
        scope='All17 original generic chart providers queried live from one accepted native seed, with common factored coordinate covers and original five histories/P0. No point selection, changed loop/integral/control functions, common finite N or global cone admission.',
        input_hashes=hashes)
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('All17 native original generic source backends written',flush=True)
    return result


if __name__=='__main__':run()
