"""Native Rp m15 versus the corrected once-only reserved source overlay.

This creates physical ordinary-C1 function expressions. The corrected graph
has its own namespace; its integer handles never alias native graph handles.
No legacy numerical provider is reconfigured and no old error packet changes.
"""
from __future__ import annotations

from fractions import Fraction
import ast
import hashlib
import json

PREFIX = 'lei_ren_part1_paper_compliant_current_original_'
PARENT_REPORT = PREFIX + 'native_Rc_2Rc_pressure_join_lift.json'
SIGNED_REPORT = PREFIX + 'native_Rc_Rp_signed_pressure_source.json'
PULSE = PREFIX + 'Rp_native_C5_pulse_source.json.gz'
PULSE_CHECK = PREFIX + 'Rp_native_C5_pulse_source_check.json'
PULSE_SOURCE = PREFIX + 'Rp_native_C5_pulse_source.py'
FULL = PREFIX + 'n1_full_defect_linear_functions.json.gz'
CORRECTED = PREFIX + 'corrected_reserved_m15_source.json'
PROVIDER = PREFIX + 'corrected_own_swirl_provider.json'
QUIET = PREFIX + 'true_quiet_Rp_C5_pressure_C3_functions.json.gz'
PULSE5 = PREFIX + 'n1_pulse_pressure_C3_functions.json.gz'
PULSE3 = PREFIX + 'n1_pulse_C3_functions.json.gz'
SCHEMA = 'original_native_Rp_m15_source_comparison_v1'
RECORD_HASHES = {
    PARENT_REPORT: '1f4b6aedf6683299117fb3bc2cddb90241a216f2b83296dddba04f60b53df39f',
    SIGNED_REPORT: 'a8d76d71b7d527c65c4eb29cc4e8340fc6592723246a87df65e496d3bce43a59',
    PULSE: '442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe',
    PULSE_CHECK: '320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15',
    PULSE_SOURCE: '369e6404d885dd4548e4f6451f2d5a5ed4160c0abb79a4aca2849520ede03541',
    FULL: 'd6b193bf527d4e8cd4afe415b1e941c928f32288b63362783074381ef7480fc1',
    CORRECTED: 'fa57cdda6a9de6e6b32fb1b2f265d663825662be959b8a3370d56432f6ff07a5',
    PROVIDER: 'bf4ab120edd8dc02925977dc37614215a2be00a83b8393d1d5763c9ad50cba42',
}
DECODED_RECORD_HASHES = {
    PULSE: '56137195254ae71d7eedfccc82ac7b5c048c8f3acd8de87fedfee5c6725fd7b1',
    FULL: '739f8ad50179cfe75d90fcc1a9beed13bf1d76e6cfa79a64f45d8aa2d057dbf0',
}
STAGE_FLAGS = (
    'native_Rp_C5_pulse_inlet_source_identity_proved',
    'native_complete_pulse_m15_source_identity_proved',
    'registered_reserved_m15_source_bound_once',
    'native_m15_signed_comparison_defined',
)
FALSE_FLAGS = (
    'same_density_certificate', 'finite_N_history_identity', 'physical_Q_y_Z1',
    'global_five_moment_closure', 'legacy_feedback_reference_alias_proved',
    'prior_adapter_geometry_alias_proved', 'prior_adapter_P0_alias_proved',
    'time_order_recursion_claimed', 'numeric_values_computed', 'old_operators_replayed',
    'full_goal_complete', 'old212_adapter_source_identity_proved',
    'independent_terminal_function_match_complete',
    'native_vs_registered_preRp_m15_identity_proved',
)
NATIVE_GRAPH = '5cdd0c8fd6a758b9e176bf67bd29ea2f576c4ddae2c7cb63d438ea1ca8388c1e'
QUIET_GRAPH = '15d51cd83dd82f41eba990979ace934bd5c8d5c339cd516b9fc96bc20a64af56'
FULL_GRAPH = '17103716c6843d2152d0764b7662d9341592646dd3f8d73733e893a9a44c208f'
CORRECTED_GRAPH = 'bb8c47eb39b4bfce588a7ffe4311d07acf23f7098d8f5d49ee07859ea2525834'
FRAME_KEYS = ('u', 'm1', 'm2', 'X', 'energy', 'Mp', 'P0', 'pressure')


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode()).hexdigest()


class PhysicalComparisonGraph:
    """Affine combinations of new source ports, never archived graph values."""

    def __init__(self):
        self.nodes, self.forms, self.ports = [], {}, {}

    def append(self, operation, arguments, form=None, **data):
        ref = SCHEMA + ':s' + str(len(self.nodes))
        self.nodes.append(dict(id=ref, operation=operation, arguments=arguments, **data))
        self.forms[ref] = {ref: Fraction(1)} if form is None else form
        return ref

    def port(self, namespace, source_digest, handle, order):
        key = (namespace, source_digest, handle, order)
        if key not in self.ports:
            self.ports[key] = self.append('physical_source_function_port', [],
                namespace=namespace, source_graph_digest=source_digest, source_node=handle,
                ordinary_Z_order=order, physical_m15=True, bound_or_midpoint=False)
        return self.ports[key]

    def linear(self, terms):
        form = {}
        for coefficient, ref in terms:
            for atom, value in self.forms[ref].items():
                form[atom] = form.get(atom, Fraction(0)) + Fraction(coefficient) * value
        form = {key: value for key, value in form.items() if value}
        return self.append('signed_function_sum', [r for _, r in terms], form,
            coefficients=[str(Fraction(c)) for c, _ in terms])

    def normal(self, ref, aliases=None):
        aliases = aliases or {}
        result = {}
        for atom, value in self.forms[ref].items():
            key = aliases.get(atom, atom)
            result[key] = result.get(key, Fraction(0)) + value
        return {k: str(v) for k, v in sorted(result.items()) if v}


def derive(source_records, *, expected_input_digest):
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)

    check('input_digest', digest(source_records) == expected_input_digest)
    for name, expected in RECORD_HASHES.items():
        record = source_records[name]
        check('raw_and_decoded:' + name, record['raw_authenticated_before_decode'] is True
            and record['raw_sha256'] == expected and record['decoded_sha256'] == DECODED_RECORD_HASHES.get(name, expected)
            and hashlib.sha256(record['source_text'].encode()).hexdigest() == record['decoded_sha256'])
    records = {n: json.loads(v['source_text']) for n, v in source_records.items() if not n.endswith('.py')}
    parent, signed = records[PARENT_REPORT], records[SIGNED_REPORT]
    native, receipt = records[PULSE], records[PULSE_CHECK]
    full, corrected, provider = records[FULL], records[CORRECTED], records[PROVIDER]
    parent_result, signed_result = parent['component_result'], signed['component_result']
    check('accepted_native_pressure_parents', parent_result['passed'] is True
        and parent_result['native_Rc_Rp_natural_pressure_composition_proved'] is True
        and signed_result['passed'] is True and signed_result['signed_pressure_cocycle_proved'] is True)
    nodes, old_nodes = native['exact_graph_nodes'], full['exact_graph_nodes']
    check('native_pulse_graph', len(nodes) == 280182 and digest(nodes) == NATIVE_GRAPH)
    check('native_quiet_prefix', digest(nodes[:265856]) == QUIET_GRAPH)
    check('actual_old_full_common_prefix', len(old_nodes) == 53398 and digest(old_nodes) == FULL_GRAPH
        and nodes[:53398] == old_nodes)
    binding = native['Rp_native_C5_pulse_source_bindings']
    interface = signed_result['function_source_interface']
    quiet = binding['inputs'][QUIET]
    check('actual_current_Rp_frame_owner', quiet == interface['quiet']['source']
        and binding['current_frame'] == quiet['Rp_frame'] == interface['Rp_frame']
        and quiet['Rp_endpoint'] == interface['native_Rp'])
    check('eight_genuine_C5_source_inputs', set(binding['current_frame']) == set(FRAME_KEYS)
        and all(len(v) == 6 for v in binding['current_frame'].values()))
    check('native_source_receipt', receipt['all_passed'] is True
        and receipt['input_hashes'][PULSE] == RECORD_HASHES[PULSE]
        and all(value['passed'] is True for value in receipt['checks'].values())
        and receipt['current_original_actual_Rp_C5_frame_equals_canonical_native_smooth_functions'] is True)
    check('all_Z_frame_identity_and_source_callback',
        len(receipt['checks']['compact_quiet']['actual_eight_Rp_frame_C5_rows']) == 48
        and receipt['checks']['callback_and_consumer']['actual_current_eight_function_source_callback_rows'] == 48
        and receipt['checks']['callback_and_consumer']['exact_expression_substitution_not_a_second_solution_or_numerical_center'] is True)
    check('same_physical_source_parameters', binding['geometry']['N'] == 18
        and binding['geometry']['Pstar'] == 33572 and binding['geometry']['delta'] == 25220
        and binding['geometry']['Rp'] == quiet['Rp'] and binding['current_P0'] == quiet['P0'])
    check('no_second_physical_scale', binding['definition']['normalized_frame_not_multiplied_by_Pstar'] is True
        and binding['definition']['physical_velocity_Pstar_and_pressure_Pstar_squared_once'] is True)
    check('same_source_family', native['source_family'] == receipt['source_family'] == full['source_family']
        == parent_result['source_family'] == provider['corrected_own_swirl_source_inputs']['source_family'])
    source_ast = ast.parse(source_records[PULSE_SOURCE]['source_text'])
    check('source_substitution_owner_not_run', any(isinstance(n, ast.ClassDef) and n.name == 'PulseSourceSubstitution' for n in source_ast.body)
        and 'canonical_native_function_identity_required_before_acceptance' in source_records[PULSE_SOURCE]['source_text'])
    old_full = full['full_defect_linear_bindings']
    physical_pairs = {name: old_full['segments'][name]['m15'] for name in ('inner','lower','preRp','pulse')}
    old_pulse = binding['inputs'][PULSE5]['pressure']['m15']
    current_pulse = binding['pressure']['m15']
    check('independent_complete_m15_source', old_pulse[:2] == physical_pairs['pulse'] == [35556,35557]
        and len(old_pulse) == len(current_pulse) == 4
        and all(binding['rewritten_nodes'][str(a)] == b for a,b in zip(old_pulse,current_pulse)))
    check('all_three_complete_pulse_charts', set(binding['pulse_coordinates']) == {'main','gap','end'}
        and set(binding['leading_fields']) == {'main','gap','end'})
    old_native = binding['inputs'][PULSE3]['native']
    check('all_original_native_C3_inputs_substituted', binding['old_native'] == old_native
        and all(binding['input_substitution'][str(handle)] == binding['current_frame'][key][order]
            for key, rows in old_native.items() for order,handle in enumerate(rows) if handle not in (0,1)))
    source = corrected['corrected_reserved_m15_source']
    patch, roots = source['graph_patch'], source['physical_C1_new_roots']
    check('qualified_corrected_source_namespace', source['passed'] is True and source['source_only'] is True
        and source['source_executed'] is False and source['accepted_ancestors_replayed'] is False
        and patch['base_graph_node_count'] == 53444 and patch['full_node_count'] == 53488
        and patch['full_canonical_nodes_sha256'] == CORRECTED_GRAPH
        and source['full_prefix_sha256'] == 'ece7f7faf63cb7aaaaf3148c103e214ba7700b1d2fad0a73b6c7734a8ba89b66'
        and corrected['corrected_terminal_source_inputs']['full_graph_sha256'] == RECORD_HASHES[FULL])
    check('provider_reuses_same_corrected_patch', provider['corrected_reserved_m15_source'] == source
        and provider['corrected_terminal_source_inputs']['source_pairs'] == corrected['corrected_terminal_source_inputs']['source_pairs'])
    source_pairs = corrected['corrected_terminal_source_inputs']['source_pairs']
    check('corrected_inputs_use_actual_old_full_ports', all(source_pairs[name] == rows for name,rows in physical_pairs.items())
        and source_pairs['full'] == old_full['moments']['m15'])
    for order in range(2):
        full_ref, reserved_ref = roots['corrected_full_m15'][order], roots['corrected_reserved_m15'][order]
        check('actual_corrected_full_definition_Z' + str(order), patch['appended_nodes'][full_ref-53444] == dict(operation='sum',
            arguments=[physical_pairs['inner'][order],physical_pairs['lower'][order],roots['corrected_preRp_m15'][order],physical_pairs['pulse'][order]]))
        check('actual_corrected_reserved_definition_Z' + str(order), patch['appended_nodes'][reserved_ref-53444]
            == dict(operation='negative',argument=full_ref))
    overlay = provider['corrected_own_swirl_provider']
    leaf = overlay['response_source_leaf']
    ledger = provider['corrected_own_swirl_source_inputs']['ledger']
    check('accepted_once_only_source_overlay', overlay['passed'] is True and overlay['source_only_fields'] is True
        and overlay['source_slot_replaced_once'] is True and overlay['extra_own_swirl_added'] is False
        and leaf['own_swirl_registration_count'] == 1 and leaf['name'] == 'reserved_n1_own_swirl_once'
        and leaf['physical_C1_roots'] == roots['corrected_reserved_m15']
        and leaf['physical_full_C1_roots'] == roots['corrected_full_m15']
        and leaf['source_identity'] == 'reserved_corrected=-full_corrected' and leaf['exact_angular_response'] == [0,0,1])
    check('legacy_runtime_not_reconfigured', overlay['full_forward_chain_uses_new_provider'] is False
        and overlay['legacy_default_provider_chain_reconfigured'] is False
        and leaf['current_forward_germ_join_admitted'] is False)
    unit = provider['corrected_own_swirl_source_inputs']['unit_receipt']
    check('physical_unit_and_ordinary_derivative', unit['positive'] is True and unit['Z_independent'] is True
        and unit['ordinary_Z_derivative'] is True and unit['unit_node'] == leaf['unit_node'] == 53427
        and unit['external_graph_ports']['U_Z_exact_zero'] is True and leaf['unit_Z1_exact_zero'] is True
        and leaf['physical_measure'] == '-R*theta0*Utheta1 ds' and leaf['ordinary_Z_orders'] == [0,1])
    check('original_registration_and_six_channels_preserved',
        [ledger[key] for key in ('raw_count','forward_count','reduced_count','raw_reserved_index','reduced_reserved_index')] == [213,212,183,212,182]
        and ledger['six_channels_are_five_signed_plus_core_swirl_error_only'] is True
        and all(row['all_212_forward_and_error_sectors_retained'] is True for row in ledger['samples']))
    passed = all(checks.values())
    graph, outputs, identities = PhysicalComparisonGraph(), [], []
    for order in (0,1):
        native_port = lambda handle: graph.port('native_pulse_master',NATIVE_GRAPH,handle,order)
        corrected_port = lambda handle: graph.port('corrected_reserved_physical',CORRECTED_GRAPH,handle,order)
        inner, lower, pulse_old = [native_port(physical_pairs[key][order]) for key in ('inner','lower','pulse')]
        pre_current = native_port(quiet['Rp_endpoint']['pressure']['m15'][order])
        pre_corrected = corrected_port(roots['corrected_preRp_m15'][order])
        reserved = corrected_port(roots['corrected_reserved_m15'][order])
        pulse_current = native_port(current_pulse[order])
        forward = graph.linear([(1,inner),(1,lower),(1,pre_current),(1,reserved)])
        target = graph.linear([(-1,pulse_current)])
        discrepancy = graph.linear([(1,forward),(-1,target)])
        full_corrected = graph.linear([(1,inner),(1,lower),(1,pre_corrected),(1,pulse_old)])
        reserved_expanded = graph.linear([(-1,full_corrected)])
        discrepancy_expanded = graph.linear([(1,inner),(1,lower),(1,pre_current),(1,reserved_expanded),(1,pulse_current)])
        pre_difference = graph.linear([(1,pre_current),(-1,pre_corrected)])
        pulse_difference = graph.linear([(1,pulse_current),(-1,pulse_old)])
        reduced = graph.linear([(1,pre_difference),(1,pulse_difference)])
        aliases = {pulse_current:pulse_old}
        algebra_ok = graph.normal(discrepancy_expanded) == graph.normal(reduced)
        pulse_ok = graph.normal(reduced,aliases) == graph.normal(pre_difference)
        check('new_signed_comparison_Z' + str(order), algebra_ok and pulse_ok)
        proof = graph.append('qualified_signed_source_comparison', [discrepancy,reserved_expanded,reduced,pre_difference],
            admitted=passed and algebra_ok and pulse_ok, ordinary_Z_order=order,
            reserved_port=reserved, saved_reserved_expansion_proved=True,
            native_pulse_equals_old_pulse_source_proved=passed, terminal_difference_claimed_zero=False)
        outputs.append(dict(ordinary_Z_order=order, forward_m15=forward, independent_negative_complete_pulse_target=target,
            discrepancy=discrepancy, corrected_reserved_port=reserved, reserved_expansion=reserved_expanded,
            discrepancy_expanded=discrepancy_expanded, preRp_difference=pre_difference, pulse_difference=pulse_difference,
            reduced_difference=reduced, contraction=proof, canonical_difference=graph.normal(reduced,aliases),
            physical_scale='physical signed m15 source ports; Pstar and radial measures already applied once',
            registered_source_overlay_count=1, source_overlay_is_legacy_runtime=False))
        identities.append(dict(ordinary_Z_order=order, signed_expansion_proved=algebra_ok,
            native_pulse_congruence_admitted=passed, reduced_to_preRp_difference=pulse_ok, remaining_difference_zero_proved=False))
    passed = all(checks.values())
    result = dict(schema=SCHEMA, passed=passed, source_only=True, ordinary_Z_orders=[0,1],
        ordinary_Z_not_temporal_n=True, input_digest=expected_input_digest, source_family=native['source_family'],
        source_prerequisite_audit=checks, definition_graph=dict(namespace=SCHEMA,nodes=graph.nodes,node_count=len(graph.nodes),
            nodes_digest=digest(graph.nodes),outputs=outputs,proof_graph_not_numeric_evaluator=True),
        comparison_identities=identities, source_namespaces=dict(native=dict(digest=NATIVE_GRAPH,node_count=280182),
            corrected_physical=dict(digest=CORRECTED_GRAPH,node_count=53488,full_prefix_count=53444)),
        source_interface=dict(native_Rp_frame=binding['current_frame'],native_Rp_pressure=quiet['Rp_endpoint']['pressure'],
            native_complete_pulse_pressure=binding['pressure'],registered_corrected_physical_roots=roots,
            registered_unit_node=53427,raw_reserved_index=212,reduced_reserved_index=182,
            corrected_source_overlay_installed=True,legacy_numeric_chain_reconfigured=False,
            main_gap_end_coordinates=binding['pulse_coordinates'],independent_m15_not_P1=True),
        source_witness=dict(native_frame_C5_rows=48,complete_pulse_m15_substitution=list(zip(old_pulse,current_pulse)),
            full_prefix_equal=True,corrected_source_overlay_identity='reserved_corrected=-full_corrected',
            reduced_identity='Delta=(native_preRp-corrected_preRp)+(pulse_native-pulse_old)=native_preRp-corrected_preRp',
            current_vs_corrected_preRp_source_identity_unproved=True, original_common_source_domain_required=True,
            genuine_saved_source_regularities_required=True, numerical_containment_claimed=False),
        retained_construction_data=parent['retained_construction_data'],
        ledger_provenance=dict(parent_report=ledger['parent_report'],parent_sha256=ledger['parent_sha256'],
            raw_count=213,forward_count=212,reduced_count=183,raw_reserved_index=212,reduced_reserved_index=182,
            sample_digests=[{k:v for k,v in row.items() if k.endswith('_digest') or k=='saved_Z'} for row in ledger['samples']]),
        source_inspection_audit=dict(new_source_consumer_calls=1,accepted_ancestor_operator_calls=0,
            old_definition_graphs_rebuilt=False,old_operators_replayed=False,numeric_values_computed=False,
            quadrature_or_norm_calls=0,source_value_reducer_calls=0,temporal_recursion_calls=0),
        legacy_numeric_chain_reconfigured=False,old_error_terms_removed=False,extra_reserved_contribution_added=False,
        open_obligations=['native_preRp equals registered corrected preRp as physical C1 source functions',
            'legacy212 numerical packet source binding','independent terminal function match and global five moments',
            'other physical joins/exact heat/stress/temporal recursion/corrections/Cartesian NS'],
        unresolved_gates=[k for k,v in checks.items() if not v],runtime_model_metadata=None)
    result.update({flag:passed for flag in STAGE_FLAGS})
    result.update({flag:False for flag in FALSE_FLAGS})
    return result
