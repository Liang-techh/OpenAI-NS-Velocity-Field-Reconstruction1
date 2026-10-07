"""Live original bridge source covers without the downstream tensor graph.

Only original parent-acquisition and raw-row operators run. Their already
checked normalization is retained, followed by the existing Pstar conversion.
This is a live coordinate cover backend, not a selected point-function backend
or a replacement admission for CurrentSourcePackets' full tensor owner API.
"""
import gzip
import json
import math
import time
from pathlib import Path
import mpmath as mp

import lei_ren_part1_paper_compliant_current_generic_shear_source_packets as packets
from lei_ren_part1_paper_compliant_current_actual_bridge_mixed_C4 import CurrentActualBridgeMixedC4
from lei_ren_part1_paper_compliant_current_bridge_stress_operator import (
    compiled_current_bridge_with_rows, compiled_raw_bridge_rows)
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_current_checked_source_runtime import (
    CheckedSourceRuntime, native_bridge_owner)

HERE, PREFIX, sha = packets.HERE, packets.PREFIX, packets.sha
NAME = PREFIX + 'current_native_generic_left_inlet.json'
VIEWS = PREFIX + 'current_native_generic_left_inlet_views.json.gz'
RECEIPT = PREFIX + 'current_native_generic_left_inlet_check.json'
GATE = 'current_original_generic_left_inlet_live_source_cover_query_executed'
TENSOR_RECEIPT = PREFIX + 'current_bridge_background_tensor_check.json'
DOMAINS = {'bridge_first': (0,1), 'bridge_second': (1,2), 'bridge_macro': (0,1)}


class NativeBridgeSourcePackets:
    """Three original bridge charts, with an explicitly supplied native owner."""
    def __init__(self, bridge):
        if type(bridge) is not CurrentActualBridgeMixedC4 or not bridge.acceptance_loaded:
            raise ValueError('Exact existing accepted native bridge required')
        if not all(bridge.current_provider_graph().values()):
            raise ValueError('One original current bridge object graph required')
        self.bridge = bridge
        self.ctx = bridge.ctx
        self.family = dict(zip(packets.FAMILY_KEYS,(bridge.family,bridge.source,bridge.datum_sha)))
        self.service = packets.CurrentSourcePackets()
        if self.service.family != self.family:
            raise ValueError('Native bridge and generic source family/datum differ')
        self.service.bind_hashes(bridge.hashes)
        receipt = accepted(TENSOR_RECEIPT,bridge.family,bridge.source,
            'current_actual_three_bridge_full_tensors_available')
        if receipt['datum_enclosure_sha256'] != bridge.datum_sha:
            raise ValueError('Original raw-row normalization datum differs')
        self.service.bind_hashes(receipt['input_hashes'])
        self.service.bind_hashes({TENSOR_RECEIPT:sha(TENSOR_RECEIPT)})
        proof = receipt['current_actual_bridge_source_and_tensor_theorem']
        normalization = proof['current_actual_bridge_source_pressure_units_and_three_endpoint_theorem'][
            'original_current_micro_and_macro_raw_unit_theorem']
        if not normalization['passed'] or normalization['total_identities'] != 282:
            raise ValueError('Original microscopic/macro raw-row normalization theorem required')
        self.evaluator,self.compiler_proof = compiled_current_bridge_with_rows()
        self.adapter,self.adapter_proof = compiled_raw_bridge_rows()
        if (self.compiler_proof != proof['original_output_only_compiler_theorem'] or
                self.adapter_proof != proof['original_current_raw_unit_adapter_theorem']):
            raise ValueError('Live original raw-row compiler/adapter differs from acceptance')
        self.service.bind_hashes(self.compiler_proof['input_hashes'])
        self.service.bind_hashes(self.adapter_proof['input_hashes'])
        self.unit_proof = packets.unit_theorem()
        self.service.bind_hashes({Path(__file__).name:sha(Path(__file__).name)})
        self.queries = []

    @source_precision
    def query(self, chart, Z, coordinate):
        if chart not in DOMAINS:
            raise ValueError('Original native first/second/macro bridge chart required')
        if not self.bridge.acceptance_loaded or not all(self.bridge.current_provider_graph().values()):
            raise ValueError('Original checked native bridge graph changed')
        c=self.ctx;z=c.mpf(Z);v=c.mpf(coordinate);ep=packets.recovery.endpoints
        zl,zh=ep(z);vl,vh=ep(v);lo,hi=DOMAINS[chart]
        if not all(mp.isfinite(x) for x in (zl,zh,vl,vh)) or zl < -1 or zh > 1 or vl < lo or vh > hi:
            raise ValueError('Original chart domain and whole Z source interval required')
        original = self.evaluator(self.bridge,z,v,chart.replace('bridge_',''))
        raw = self.adapter(c,original)
        if not raw['original_radial_prefactors_included_once'] or not raw['full_absolute_P0_included_once']:
            raise ValueError('Original radial factors and absolute pressure must be retained once')
        if chart == 'bridge_macro':
            if not raw['macro_original_ordinary_y_preserved']:
                raise ValueError('Macro rows must retain their original ordinary logR coordinate')
        elif not raw['inverse_width_source_shift_applied_before_any_resolution']:
            raise ValueError('True original microscopic width inversion required')
        logs=raw['algebra'].logs
        if len(logs)!=4 or any(not mp.isfinite(x) for q in logs for x in ep(q)):
            raise ValueError('Four finite formal original source log bases required')
        algebra=packets.FactoredAlgebra(c,logs,[])
        def rows(values):
            if len(values)!=5:
                raise ValueError('Five original ordinary logR rows required')
            converted=tuple(packets.decode_row(algebra,packets.factored_rows_record(value)) for value in values)
            packets.factored_check_jets(c,converted)
            return converted
        native_v={key:rows(value) for key,value in raw['velocity'].items()}
        native_m={key:rows(value) for key,value in raw['histories'].items()}
        if set(native_v)!=set(('theta','axial','radial')) or set(native_m)!=set(packets.recovery.RATES):
            raise ValueError('Three original velocities and all five original histories required')
        velocity={key:tuple(algebra.shift(row,packets.INVERSE_S) if key in ('axial','radial') else row
                  for row in value) for key,value in native_v.items()}
        histories={key:tuple(algebra.shift(row,packets.INVERSE_S) if key in ('m','k') else row
                   for row in value) for key,value in native_m.items()}
        pressure=rows(raw['absolute_pressure'])
        p0=algebra.lift(packets.jet(c,original['actual_parent_axial5_packet']['pressure_axis_axial5_coefficients'][:6]))
        if p0.order!=5:
            raise ValueError('Independent original normalized axis pressure axial5 required')
        provenance=dict(mode='injected_native_original_source_cover',chart=chart,Z_box=z,coordinate_box=v,
            source_family=self.family,owner_class=type(self.bridge).__name__,owner_module=type(self.bridge).__module__,
            original_source_graph_asserted=True,arbitrary_coordinates_evaluated=True,cache_cover=False,
            derivative_coordinate='ordinary y=log R',
            coordinate_contract='original macro coverage selector; rows are ordinary logR' if chart=='bridge_macro'
                else 'phase; D_y^j=original_hb^-j D_phase^j before normalization',
            original_radial_and_pressure_factors_counted_once=True,
            current_swirl_factor_not_applied_twice=True,source_factor_resolution_performed=False,
            P0_read_from_original_source_not_subtracted=True,
            full_tensor_or_physical_point_evaluation_requested=False,
            native_operator_outputs_are_covers_not_selected_point_values=True)
        self.queries.append(provenance)
        return packets.CurrentSourcePacket(chart,dict(self.family),algebra,
            packets.IntervalTaylor.variable(c,z,5),velocity,histories,pressure,p0,native_v,native_m,provenance)

    def left_inlet(self, Z=(-1,1)):
        # This is the positive interior s_c/2 inlet, not the core/bridge seam0.
        import lei_ren_part1_paper_compliant_current_generic_five_defect_bounds as current
        receipt=json.loads((HERE/current.RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt[current.GATE] or receipt['source_family']!=self.family:
            raise ValueError('Same checked generic inlet/domain family required')
        self.service.bind_hashes(receipt['input_hashes'])
        self.service.bind_hashes({current.RECEIPT:sha(current.RECEIPT)})
        for stem,gate in (
            ('current_inner_exit_strict_collar','current_inner_exit_strict_collar_attached_to_current_source_graph_certified'),
            ('current_generic_shear_loop_domain','current_original_generic_loop_two_sided_collars_and_repair_geometry_certified'),
            ('current_generic_shear_uniform_inputs','current_original_whole_generic_input_margins_and_log_scales_certified'),
            ('current_generic_loop_function_sources','current_original_generic_loop_velocity_and_five_increment_function_graphs_defined')):
            name=PREFIX+stem+'_check.json';admitted=json.loads((HERE/name).read_bytes())
            family=admitted.get('source_family') or {k:admitted[k] for k in packets.FAMILY_KEYS}
            if not admitted.get('all_passed') or not admitted.get(gate) or family!=self.family:
                raise ValueError('Same-family accepted inlet/flat-branch prerequisite required: '+stem)
            self.service.bind_hashes(admitted['input_hashes'])
            self.service.bind_hashes({name:sha(name)})
        left_name=PREFIX+'current_inner_exit_strict_collar.json'
        left=json.loads((HERE/left_name).read_bytes())
        read=lambda row:packets.interval(self.ctx,row)
        phase=read(left['explicit_current_inner_exit_strict_collar']['selected_first_phase_endpoint'])/2
        ep=packets.recovery.endpoints
        if ep(phase)[0]<=0 or ep(phase)[1]>=1:
            raise ValueError('Strict positive original first-bridge interior inlet required')
        domain_name=PREFIX+'current_generic_shear_loop_domain.json'
        uniform_name=PREFIX+'current_generic_shear_uniform_inputs.json'
        domain=json.loads((HERE/domain_name).read_bytes())['current_original_generic_loop_domain']
        uniform=json.loads((HERE/uniform_name).read_bytes())['current_actual_logarithmic_loop_scales']
        eta_log=read(uniform['selected_positive_eta_log'])
        if domain['source_family']!=self.family or not domain['original_both_loop_edge_collars_strict']:
            raise ValueError('Same current two-sided strict generic loop collars required')
        left_excess=read(domain['left_collar_kappa_minus2_lower'])
        if ep(left_excess)[0]<=0 or ep(eta_log)[1]>ep(self.ctx.ln(left_excess)-self.ctx.ln(2))[0]:
            raise ValueError('Actual eta must lie below half of the positive left kappa excess')
        if ep(eta_log)[1] > ep(read(domain['allowed_eta_log_upper_for_automatic_constant_edges']))[0]:
            raise ValueError('Same actual eta and exact flat left-edge collar required')
        self.service.bind_hashes({name:sha(name) for name in (left_name,domain_name,uniform_name)})
        packet=self.query('bridge_first',Z,phase)
        state=packets.FactoredRecoveryState(packet)
        return dict(original_source_packet=packet,original_inlet_history_covers=state.original,
            original_axis_pressure_over_Pstar_squared=packet.P0,generic_inlet_defect_exact_zero=state.defect,
            source_family=self.family,source_query_coordinate=phase,
            exact_formal_radius_binding=dict(radius='Ra*exp(original_hb*s_c/2)',
                source=domain['formal_left_log_radius_offset'],positive_microscopic_offset_not_added_to_huge_log_Ra=True),
            original_P0_and_all_five_incoming_histories_retained=True,
            left_loop_q_A_B_exact_zero_by_existing_checked_collar=True,
            flat_left_branch_source_argument=dict(
                identity='Delta=H0-2=D+kappa-2, D>=0; Delta>=left_kappa_excess>=2*eta',
                source='accepted same-family loop-domain, uniform-input identity and generic flat-branch graph receipts',
                branch='q=A=B=0 when Delta>=eta; no active inverse/loop evaluated on flat collar'),
            generic_inlet_defect_zero_is_declared_Section11_initial_condition=True,
            defect_transport_or_terminal_closure_not_inferred_from_initial_zero=True,
            returned_rows_are_covers_not_selected_point_function_values=True,
            full_tensor_owner_API_gate_or_physical_point_scope_not_promoted=True,
            input_hashes=dict(self.service.hashes))


def run(bridge=None, runtime_evidence=None):
    started=time.monotonic()
    if bridge is None:
        bridge,runtime_evidence=native_bridge_owner()
        print('Genuine native bridge source owner ready',flush=True)
    with CheckedSourceRuntime():
        backend=NativeBridgeSourcePackets(bridge)
        inlet=backend.left_inlet()
        fresh=backend.query('bridge_first','.537','.1337')
        # Independent coordinates exercise the same microscopic and macro
        # original operators; all outputs remain directed coefficient covers.
        second=backend.query('bridge_second','-.317','1.831')
        macro=backend.query('bridge_macro','.731','.537')
    packet=inlet.pop('original_source_packet')
    views=dict(left_inlet=dict(inlet,original_source_packet=packet.record()),
        fresh_first=fresh.record(),fresh_second=second.record(),fresh_macro=macro.record())
    payload=(json.dumps(packets.encode(views),separators=(',',':'))+'\n').encode()
    (HERE/VIEWS).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    hashes=dict(backend.service.hashes)
    hashes.update({VIEWS:sha(VIEWS),Path(__file__).name:sha(Path(__file__).name),
        PREFIX+'current_checked_source_runtime.py':sha(PREFIX+'current_checked_source_runtime.py')})
    result=dict(source_family=backend.family,**{GATE:True},
        successful_live_coordinate_query_count=len(backend.queries),query_provenance=backend.queries,
        original_native_graph_assertions=bridge.current_provider_graph(),
        exact_native_owner_type=type(bridge).__name__,native_source_owner_acceptance_loaded=bridge.acceptance_loaded,
        runtime_construction_evidence=runtime_evidence,
        original_common_unit_theorem=backend.unit_proof,
        unchanged_raw_row_compiler_theorem=backend.compiler_proof,
        unchanged_raw_row_adapter_theorem=backend.adapter_proof,
        native_inlet_source_and_common_factored_cover_views=VIEWS,
        five_original_inlet_histories_and_separate_P0_retained=True,
        left_inlet_coordinate_is_positive_sc_over_2_not_seam_zero=True,
        numerical_caps_used_as_defining_function_values=False,
        actual_defect_integral_functions_installed=False,actual_repair_control_functions_installed=False,
        actual_terminal_Z_function_closure_installed=False,
        source_function_or_full_tensor_owner_interface_gate_promoted=False,
        **dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-started,
        scope='Actual injected native source queries on three original bridge charts, including strict interior generic inlet; formal coefficient covers only. No changed loop functions/integrals, full tensor owner, point values, repair controls, common N or global cone.',
        input_hashes=hashes)
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Native generic inlet source rows written',flush=True)
    return result


if __name__=='__main__':run()
