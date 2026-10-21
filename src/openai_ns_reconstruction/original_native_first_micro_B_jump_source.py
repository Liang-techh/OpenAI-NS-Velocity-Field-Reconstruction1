"""Bind the saved free-Z micro seam theorem to the actual native IBP jump.

Only new source-port congruence is constructed.  No archived algebra,
fingerprint, derivative, integral, or numerical evaluator is executed.
"""
import ast
import hashlib
import json

PREFIX = 'lei_ren_part1_paper_compliant_current_original_'
SCHEMA = 'original_native_first_micro_B_jump_source_v1'
PARENT = PREFIX + 'native_active_S1_weighted_ibp.json'
LEGACY = PREFIX + 'n1_micro_correlated_m15.json'
MOD = PREFIX + 'modified_Q4_pressure_C3_functions.json.gz'
LOOP = PREFIX + 'prescribed_loop_C5_functions.json.gz'
FINGERPRINT = PREFIX + 'micro_viscous_m15_numeric.json'
SEAM_CHECK = PREFIX + 'micro_active_radial4_pressure_C3_functions_check.json'
LEGACY_SOURCE = 'src/openai_ns_reconstruction/original_n1_micro_correlated_m15.py'
FINGERPRINT_SOURCE = 'src/openai_ns_reconstruction/original_micro_viscous_m15.py'
RECORD_NAMES = (PARENT, LEGACY, MOD, LOOP, FINGERPRINT, SEAM_CHECK,
                LEGACY_SOURCE, FINGERPRINT_SOURCE)
RECORD_HASHES = {
    PARENT: 'aab72bbc7a882e9bf089a43dd3407babe085f1c0ec61583578d20a5e7c9bec28',
    LEGACY: '7868f1899398a54aaeba1b799b8b1015465718445444966648295e66ae6c97aa',
    MOD: 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a',
    LOOP: '8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa',
    FINGERPRINT: '0903f9100aa46abbfd9b69c667c91b147b364cb0a33cbcc526f4549a53d7a5b8',
    SEAM_CHECK: '9cb8fd79c64ada9d447f189e6bddd21bc3ef4ccab142f7a10bd8e2cec82a24a8',
    LEGACY_SOURCE: '7c7a7022e2b059c35dff925a25b4e44b6ea45ed9139f6b2010f5d919955859f3',
    FINGERPRINT_SOURCE: 'fe514784962f2d2ab968b66c537d6e907dd5f8c1d30183123e6498ba030a33b7',
}
STAGE_FLAGS = (
    'actual_native_first_micro_B_jump_C1_source_identity_proved',
    'saved_free_Z_seam_certificate_bound_to_actual_native_ports',
    'nonzero_complete_predecessor_histories_retained',
    'native_IN_representation_retains_fifteen_other_primitive_jumps',
)
FALSE_FLAGS = (
    'all_native_B_jumps_zero_proved', 'native_active_m15_source_identity_proved',
    'native_compact_m15_source_identity_proved', 'native_old_density_identity_proved',
    'terminal_difference_zero_proved', 'global_five_moment_closure',
    'numeric_values_computed', 'old_error_terms_removed', 'old_operators_replayed',
    'legacy_numeric_chain_reconfigured', 'time_order_recursion_claimed',
    'full_goal_complete',
)
B_FORMULA = 'R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)'
B_Z_FORMULA = ('R^2*(Q_Z/L-Q*L_Z/L^2+Q*Q_Z-4*Z*U^2/L^2'
               '-4*Z^2*U*U_Z/L^2+4*Z^2*U^2*L_Z/L^3)')


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode()).hexdigest()


class Graph:
    def __init__(self):
        self.nodes = []

    def node(self, operation, arguments=(), **metadata):
        ref = len(self.nodes)
        self.nodes.append(dict(operation=operation, arguments=list(arguments), **metadata))
        return ref


def derive(source_records, *, expected_input_digest):
    checks, data, trees = {}, {}, {}

    def check(name, condition):
        checks[name] = bool(condition)

    check('frozen_input_digest', digest(source_records) == expected_input_digest)
    for name in RECORD_NAMES:
        record = source_records[name]
        text = record['source_text']
        check('authenticated:' + name,
            record['raw_sha256'] == RECORD_HASHES[name]
            and record['raw_authenticated_before_decode'] is True
            and hashlib.sha256(text.encode()).hexdigest() == record['decoded_sha256'])
        if name.endswith('.py'):
            trees[name] = ast.parse(text)
        else:
            data[name] = json.loads(text)

    parent, legacy = data[PARENT], data[LEGACY]
    parent_result = parent['component_result']
    parent_graph = parent['definition_graph']
    nodes = parent_graph['nodes']
    old = legacy['actual_micro_correlated_m15']
    seam = old['seam_certificate']
    modified, loop = data[MOD], data[LOOP]
    mb = modified['modified_Q4_pressure_C3_bindings']
    lb = loop['prescribed_loop_C5_bindings']
    a, b = 'first_micro', 'second_micro'
    first, second = lb['charts'][a], lb['charts'][b]
    fa, fb = mb['charts'][a]['physical'], mb['charts'][b]['physical']
    check('accepted_unchanged_parent_graph', parent_result['passed'] is True
        and parent_graph['node_count'] == 545
        and digest(nodes) == parent_graph['nodes_digest'] ==
            '10d77faec21c11aef2bccfeec26659e41243bd87d9c07cd941d14a93507b4729')
    check('same_source_family', all(value['source_family'] == parent['source_family']
        for value in (legacy, old, modified, loop, data[SEAM_CHECK], data[FINGERPRINT])))
    check('accepted_legacy_source_identity', legacy['checks']['passed'] is True
        and old['passed'] is True and seam['passed'] is True
        and old['internal_RQy_and_reduced_B_seams_cancel_as_functions'] is True
        and legacy['checks']['packet_invariants']['exact_same_function_internal_seam'] is True
        and old['accepted_seam_fingerprints_replayed'] is False
        and old['old_six_numeric_terms_and_all_error_caps_unchanged'] is True)
    check('legacy_native_graph_pins_identical', all(
        old['source_pins']['experiments/root_st073/' + name]['expected_sha256'] == RECORD_HASHES[name]
        for name in (MOD, LOOP, FINGERPRINT, SEAM_CHECK)))
    check('legacy_scope_free_Z_C1_not_samples',
        seam['source_graph_prefix_and_all_14_field_ports_checked'] is True
        and seam['free_Z_not_numeric_Z_substitution'] is True
        and seam['ordinary_Z_orders'] == [0, 1]
        and seam['exact_s1_seam'] is True
        and seam['accepted_fingerprint_computation_replayed'] is False
        and seam['fingerprint_scope_AST']['passed'] is True)
    check('fingerprint_report_same_accepted_function_theorem',
        seam['accepted_fingerprint_report_sha256'] == RECORD_HASHES[FINGERPRINT]
        and data[FINGERPRINT]['checks']['passed'] is True
        and data[FINGERPRINT]['checks']['source_identity_rows_per_packet'] == 14
        and data[FINGERPRINT]['checks']['direct_Qy_and_Z1_source_ids'] is True
        and data[FINGERPRINT]['checks']['exact_cancelled_interface'] == 'first_s1')
    traces = [sample['trace_identity'] for sample in data[FINGERPRINT]['samples']]
    trace = traces[0]
    check('saved_free_Z_trace_is_shared_not_boxed', all(row == trace for row in traces)
        and trace['passed'] is True and trace['signed_trace_centers_evaluated'] is False
        and trace['upper_caps_not_used_as_function_values'] is True
        and trace['common_ordinary_Z_independent_radius'] is True)
    check('raw_native_prefix_not_rebuilt',
        modified['exact_graph_nodes'][:len(loop['exact_graph_nodes'])] == loop['exact_graph_nodes'])
    check('exact_common_native_s1_geometry',
        all(first[key] == second[key] for key in
            ('coordinate', 'radius', 'physical_measure', 'actual_phase', 'roots', 'global_offset'))
        and first['domain'][1] == second['domain'][0] == 1
        and modified['exact_graph_nodes'][1] == dict(operation='exact_rational', numerator=1, denominator=1)
        and first['coordinate'] == trace['coordinate_node']
        and first['radius'] == trace['radius_node']
        and first['actual_phase'] == trace['phase_node']
        and first['physical_measure'] == trace['physical_measure_node']
        and lb['N'] == trace['selected_integer_node'] == 18
        and trace['seam_node'] == 1)
    check('actual_nonreset_predecessor_and_pressure_histories',
        second['incoming'] == lb['endpoints'][a]
        and mb['charts'][b]['incoming_integrals'] == mb['endpoints'][a]['integrals']
        and seam['complete_predecessor_history_retained'] is True
        and seam['shared_phase_N_radius_and_Pstar'] is True)
    expected = []
    for field in ('U', 'm', 'Q'):
        for radial in (0, 1):
            key = field if radial == 0 else field + '_y'
            for order in (0, 1):
                expected.append((field, radial, order, fa[key][order], fb[key][order]))
    for order in (0, 1):
        expected.append(('R_Q_y', 0, order,
            mb['charts'][a]['finite_boundary']['R_Q_y'][order],
            mb['charts'][b]['finite_boundary']['R_Q_y'][order]))
    rows = seam['source_identity_rows']
    fields = ('field', 'physical_y_order', 'ordinary_Z_order', 'left_node', 'right_node')
    check('all14_native_field_ports_match_saved_seam_theorem',
        [tuple(row[key] for key in fields) for row in rows] == expected
        and rows == trace['source_identity_rows']
        and all(len(row['canonical_endpoint_sha256']) == 64 for row in rows))
    check('saved_high_source_join_admission', data[SEAM_CHECK]['all_passed'] is True
        and data[SEAM_CHECK]['current_original_n1_first_second_micro_full_mixed_source_join_installed'] is True)
    legacy_seam_function = next(item for item in trees[LEGACY_SOURCE].body
        if isinstance(item, ast.FunctionDef) and item.name == 'seam_certificate')
    scope_function = next(item for item in ast.walk(trees[FINGERPRINT_SOURCE])
        if isinstance(item, ast.FunctionDef) and item.name == 'at')
    check('unchanged_legacy_free_Z_scope_source_bound',
        'source_identity_rows' in ast.unparse(legacy_seam_function)
        and 'self.one' in ast.unparse(scope_function)
        and 'bound_variable' in ast.unparse(scope_function))

    target = nodes[520]
    right_id, left_id = target['arguments']
    left, right = nodes[left_id], nodes[right_id]
    check('actual_parent_jump_is_first_native_seam', target['operation'] == 'actual_internal_primitive_jump'
        and target['jump_index'] == 0 and target['arguments'] == [344, 334]
        and left['chart'] == a and right['chart'] == b
        and left['endpoint_bound'] == 'domain[1]' and right['endpoint_bound'] == 'domain[0]'
        and left['formula'] == right['formula'] == B_FORMULA)

    def actual_parent_port(port_id, archive, source_values):
        row = nodes[port_id]
        return row.get('source_record') == archive and row.get('source_nodes') == source_values

    check('actual_parent_B_geometry_endpoints_bound',
        actual_parent_port(left['coordinate'], LOOP, [first['coordinate']])
        and actual_parent_port(right['coordinate'], LOOP, [second['coordinate']])
        and actual_parent_port(left['domain'], LOOP, first['domain'])
        and actual_parent_port(right['domain'], LOOP, second['domain'])
        and actual_parent_port(left['radius'], LOOP, [first['radius']])
        and actual_parent_port(right['radius'], LOOP, [second['radius']]))
    check('actual_parent_B_C1_U_Q_ports_bound', all(
        actual_parent_port(endpoint['physical_ports'][field], MOD, physical[field])
        for endpoint, physical in ((left, fa), (right, fb)) for field in ('U', 'Q')))
    mn = modified['exact_graph_nodes']
    denominators = [mn[physical[field][0]]['denominator']
        for physical in (fa, fb) for field in ('Q', 'Q_y')]
    check('same_actual_positive_L_for_both_native_B_functions',
        denominators == [25272] * 4 and nodes[18]['source_nodes'] == [25272]
        and nodes[18]['source_record'] == MOD and nodes[18]['source_positive'] is True)
    check('parent_IBP_representation_and_other15_jumps_unchanged',
        parent_graph['outputs']['native_weighted_ibp_representation'] == 537
        and parent_graph['outputs']['new_weighted_IBP_identity_certificate'] == 544
        and nodes[537]['arguments'][2:18] == list(range(520, 536))
        and nodes[544]['all16_jump_nodes'] == list(range(520, 536))
        and all(nodes[index]['jump_index'] == index - 520 for index in range(520, 536)))
    retained = parent['retained_construction_data']
    check('construction_data_frozen', retained['cumulative_terms'] == 212
        and retained['raw_signed_terms'] == 213 and retained['reduced_terms'] == 183
        and retained['six_signed_channels_retained'] is True
        and retained['all_old_error_data_retained'] is True
        and retained['reserved_contribution_registered_once'] is True
        and retained['inherited_errors_deleted'] is False)
    passed = all(checks.values())

    graph = Graph()
    parent_jump = graph.node('accepted_parent_source_function_port', source_record=PARENT,
        source_node=520, ordinary_Z_orders=[0, 1], parent_graph_digest=parent_graph['nodes_digest'])
    theorem = graph.node('accepted_free_Z_seam_theorem_port', source_record=LEGACY,
        selector='actual_micro_correlated_m15.seam_certificate',
        ordinary_Z_orders=[0, 1], accepted_fingerprint_computation_replayed=False)
    radius = graph.node('actual_common_radius_at_native_s1', source_record=LOOP,
        radius_node=first['radius'], coordinate_node=first['coordinate'],
        endpoint_graph_handle=1, endpoint_exact_value='1', ordinary_Z_independent=True)
    L = graph.node('actual_common_positive_L_port', source_record=MOD,
        source_node=25272, ordinary_Z_orders=[0, 1], formula='1-delta*Z^2')
    canonical = []
    for field, quantity in (('U', 'U'), ('Q', 'Q')):
        for order in (0, 1):
            row = next(row for row in rows if row['field'] == field
                and row['physical_y_order'] == 0 and row['ordinary_Z_order'] == order)
            canonical.append(graph.node('source_equal_native_endpoint_C1_row', (theorem,),
                quantity=quantity, ordinary_Z_order=order, native_left_node=row['left_node'],
                native_right_node=row['right_node'], canonical_endpoint_sha256=row['canonical_endpoint_sha256'],
                inherited_free_Z_identity=True, admitted=passed))
    common_B = graph.node('common_actual_B_C1_function', (radius, L, *canonical),
        formula=B_FORMULA, ordinary_Z1_formula=B_Z_FORMULA,
        L_Z='-2*delta*Z', physical_radius_Z_independent=True,
        smooth_algebra_congruence_from_U_Q_R_L=True, ordinary_Z_orders=[0, 1])
    identity = graph.node('actual_native_B_endpoint_C1_source_congruence',
        (common_B, theorem), parent_left_B_node=left_id, parent_right_B_node=right_id,
        inherited_native_U_Q_C1_ports=canonical, admitted=passed)
    zero = graph.node('actual_native_first_B_jump_zero_C1', (parent_jump, identity),
        definition='B_second_lower-B_first_upper', same_source_difference_zero=True,
        ordinary_Z_orders=[0, 1], owned_Z_cell=['0', '1/2'], admitted=passed,
        numeric_zero_substituted_for_unknown_function=False, old_packet_errors_deleted=False)
    untouched = [graph.node('unchanged_parent_internal_B_jump_port', source_record=PARENT,
        source_node=index, jump_index=index - 520, ordinary_Z_orders=[0, 1],
        source_identity_proved=False) for index in range(521, 536)]
    outer = [graph.node('unchanged_parent_IBP_source_port', source_record=PARENT,
        source_node=index) for index in (nodes[537]['arguments'][0], nodes[537]['arguments'][1])]
    K = [graph.node('unchanged_parent_chart_K_sum_port', source_record=PARENT,
        source_node=index) for index in nodes[537]['arguments'][18:]]
    representation = graph.node('native_IN_with_first_B_jump_proved_zero',
        (*outer, *untouched, *K, zero),
        definition='I_N=B_lastupper-B_firstlower-sum(J_1..J_15)+sum17(integral K_i dy)',
        parent_native_S1_node=parent_graph['outputs']['native_S1'],
        parent_representation_node=537, eliminated_parent_jump_node=520,
        elimination_basis='actual source C1 identity; not an enclosure', admitted=passed)
    result = dict(schema=SCHEMA, passed=passed, source_only=True,
        source_family=parent['source_family'], ordinary_Z_orders=[0, 1], ordinary_Z_not_temporal_n=True,
        source_prerequisite_audit=checks,
        source_interface=dict(parent_jump_node=520, parent_B_endpoints=[left_id, right_id],
            parent_graph_digest=parent_graph['nodes_digest'], actual_native_14_seam_rows=rows,
            common_native_coordinate=first['coordinate'], common_native_radius=first['radius'],
            actual_L_node=25272, endpoint_graph_handle=1, endpoint_exact_value='1',
            selected_N_graph_handle=18, selected_N_is_not_numeric_18=True,
            owned_Z_cell=['0', '1/2'], source_B_formula=B_FORMULA,
            first_jump_zero=passed, remaining_parent_jump_nodes=list(range(521, 536)),
            old_archived_parent_jump_rows_unchanged=True, old_packet_errors_deleted=False),
        definition_graph=dict(namespace=SCHEMA, nodes=graph.nodes, node_count=len(graph.nodes),
            nodes_digest=digest(graph.nodes), outputs=dict(first_B_jump_zero=zero,
                B_endpoint_source_congruence=identity, native_IN_representation=representation),
            proof_graph_not_numeric_evaluator=True, archived_graphs_rebuilt=False),
        retained_construction_data=retained,
        unresolved_gates=[name for name, value in checks.items() if not value],
        remaining_construction_gates=['fifteen_other_native_B_jump_source_bindings',
            'native_old_active_density_source_comparison', 'compact_density_source_comparison',
            'same_source_old_numerical_packet_error_theorem', 'global_five_moment_closure'],
        runtime_model_metadata=None, accepted_fingerprint_computation_replayed=False,
        accepted_IBP_algebra_reexecuted=False,
        source_inspection_audit=dict(new_source_consumer_calls=1, accepted_ancestor_operator_calls=0,
            old_operators_replayed=False, numeric_values_computed=False,
            new_actual_source_port_congruence_only=True))
    result.update({name: passed for name in STAGE_FLAGS})
    result.update({name: False for name in FALSE_FLAGS})
    return result
