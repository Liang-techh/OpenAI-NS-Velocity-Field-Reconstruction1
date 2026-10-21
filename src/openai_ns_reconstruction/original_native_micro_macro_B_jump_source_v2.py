"""Bind saved actual endpoint functions to the second native B jump.

This leaf constructs only source-function congruence. It never executes the
endpoint canonicalizer, an ancestor algebra, or a numerical evaluator.
"""
import hashlib
import json

PREFIX = 'lei_ren_part1_paper_compliant_current_original_'
SCHEMA = 'original_native_micro_macro_B_jump_source_v2'
PARENT = PREFIX + 'native_micro_macro_endpoint_semantics.json'
WEIGHTED = PREFIX + 'native_active_S1_weighted_ibp.json'
FIRST = PREFIX + 'native_first_micro_B_jump_source.json'
MOD = PREFIX + 'modified_Q4_pressure_C3_functions.json.gz'
LOOP = PREFIX + 'prescribed_loop_C5_functions.json.gz'
RECORD_NAMES = (PARENT, WEIGHTED, FIRST, MOD, LOOP)
RECORD_HASHES = {
    PARENT: 'ba2672ac070b39a0212033eeb01e80d4eb8ea3a916eb916c7c7c22e9c2d8c783',
    WEIGHTED: 'aab72bbc7a882e9bf089a43dd3407babe085f1c0ec61583578d20a5e7c9bec28',
    FIRST: '588c4bed6a73504a5b6af2be4fb8341ac6f3729db568678fed5693f36db07db9',
    MOD: 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a',
    LOOP: '8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa',
}
STAGE_FLAGS = (
    'actual_native_micro_macro_B_C1_source_identity_proved',
    'actual_native_two_B_jump_source_identities_bound',
    'native_IN_retains_fourteen_other_primitive_jumps',
    'accepted_endpoint_semantics_not_reexecuted',
)
FALSE_FLAGS = (
    'all_native_B_jumps_zero_proved', 'native_old_density_identity_proved',
    'native_active_m15_source_identity_proved', 'native_compact_m15_source_identity_proved',
    'terminal_difference_zero_proved', 'global_five_moment_closure',
    'numeric_values_computed', 'old_error_terms_removed', 'old_operators_replayed',
    'time_order_recursion_claimed', 'full_goal_complete',
)
B_FORMULA = 'R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)'
B_Z_FORMULA = ('2*R*R_Z*(Q/L+Q^2/2-2*Z^2*U^2/L^2)'
    '+R^2*(Q_Z/L-Q*L_Z/L^2+Q*Q_Z-4*Z*U^2/L^2'
    '-4*Z^2*U*U_Z/L^2+4*Z^2*U^2*L_Z/L^3)')
IN_FORMULA = 'I_N=B_lastupper-B_firstlower-sum(J_2..J_15)+sum17(integral K_i dy)'


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
    checks, data = {}, {}

    def check(name, value):
        checks[name] = bool(value)

    check('frozen_input_digest', digest(source_records) == expected_input_digest)
    check('exact_record_inventory', set(source_records) == set(RECORD_NAMES))
    for name in RECORD_NAMES:
        record = source_records[name]
        text = record['source_text']
        check('authenticated:' + name, record['raw_sha256'] == RECORD_HASHES[name]
            and record['raw_authenticated_before_decode'] is True
            and hashlib.sha256(text.encode()).hexdigest() == record['decoded_sha256'])
        data[name] = json.loads(text)
    parent, weighted, first = data[PARENT], data[WEIGHTED], data[FIRST]
    endpoint = parent['component_result']
    eg = endpoint['definition_graph']
    wg, fg = weighted['definition_graph'], first['definition_graph']
    wn, fn = wg['nodes'], fg['nodes']
    modified = data[MOD]['modified_Q4_pressure_C3_bindings']
    loop = data[LOOP]['prescribed_loop_C5_bindings']
    left_name, right_name = 'second_micro', 'frozen_macro'
    lc, rc = loop['charts'][left_name], loop['charts'][right_name]
    lp, rp = (modified['charts'][name]['physical'] for name in (left_name, right_name))
    check('same_exact_source_family', all(value['source_family'] == endpoint['source_family']
        for name, value in data.items() if name != PARENT))
    check('accepted_actual_endpoint_proof', endpoint['passed'] is True
        and endpoint['source_only'] is True
        and all(endpoint['source_prerequisite_audit'].values())
        and eg['node_count'] == 9291
        and eg['nodes_digest'] == '783c7b51eba674d445e0dd19dec20e255191a041547fa48d2a717bb7ee79388f'
        and all(endpoint[name] is True for name in (
            'actual_native_micro_macro_U_C1_identity_proved',
            'actual_native_micro_macro_m_Z2_identity_proved',
            'actual_native_micro_macro_Q_C1_identity_proved',
            'actual_native_micro_macro_R_identity_proved')))
    check('actual_free_Z_scope_not_point_values', endpoint['ordinary_Z_orders'] == [0, 1]
        and endpoint['ordinary_Z_not_temporal_n'] is True
        and endpoint['source_interface']['owned_Z_cell'] == ['0', '1/2']
        and endpoint['source_interface']['common_U_Q_binding_used'] is False
        and endpoint['source_interface']['canonical_difference_does_not_prove_function_inequality'] is True
        and endpoint['numeric_values_computed'] is False)
    check('accepted_parent_J0_identity', first['component_result']['passed'] is True
        and first['component_result']['source_interface']['first_jump_zero'] is True
        and fg['outputs']['first_B_jump_zero'] == 10 and fn[10]['admitted'] is True
        and fn[0]['source_record'] == WEIGHTED and fn[0]['source_node'] == 520)
    check('weighted_IBP_graph_unchanged', weighted['component_result']['passed'] is True
        and wg['nodes_digest'] == '10d77faec21c11aef2bccfeec26659e41243bd87d9c07cd941d14a93507b4729'
        and digest(wn) == wg['nodes_digest'] and wg['node_count'] == 545)

    rows = endpoint['endpoint_comparisons']
    expected = {}
    for field, count in (('U', 2), ('m', 3), ('Q', 2)):
        for order in range(count):
            expected[field + '_Z' + str(order)] = (field, order, lp[field][order], rp[field][order], MOD)
    expected['R_Z0'] = ('R', 0, lc['radius'], rc['radius'], LOOP)
    check('exact_eight_actual_endpoint_rows', set(rows) == set(expected))
    for label, (field, order, left_node, right_node, archive) in expected.items():
        row, left, right = rows[label], rows[label]['left'], rows[label]['right']
        check('actual_source_row:' + label, row['equal'] is True
            and row['field'] == field and row['ordinary_Z_order'] == order
            and left['source_node_ids'] == [left_node] and right['source_node_ids'] == [right_node]
            and left['coordinate_node'] == lc['coordinate'] == 66988
            and right['coordinate_node'] == rc['coordinate'] == 70265
            and left['endpoint_graph_handle'] == lc['domain'][1] == 48
            and right['endpoint_graph_handle'] == rc['domain'][0] == 0
            and left['canonical_digest'] == right['canonical_digest']
            and left['canonical_handle_digests'] == right['canonical_handle_digests']
            and archive in eg['canonical_graphs'])
    target, left, right = wn[521], wn[345], wn[355]
    check('actual_native_J1_and_B_endpoint_definitions', target['operation'] == 'actual_internal_primitive_jump'
        and target['arguments'] == [355, 345] and target['jump_index'] == 1
        and left['chart'] == left_name and left['endpoint_bound'] == 'domain[1]'
        and right['chart'] == right_name and right['endpoint_bound'] == 'domain[0]'
        and left['formula'] == right['formula'] == B_FORMULA)

    def port(index, archive, source_nodes):
        row = wn[index]
        return row.get('source_record') == archive and row.get('source_nodes') == source_nodes

    check('actual_B_U_Q_C1_source_ports', all(
        port(bound['physical_ports'][field], MOD, physical[field])
        for bound, physical in ((left, lp), (right, rp)) for field in ('U', 'Q')))
    check('actual_B_radius_coordinate_domain_ports', all(
        port(bound[field], LOOP, chart[source_key] if field == 'domain' else [chart[source_key]])
        for bound, chart in ((left, lc), (right, rc))
        for field, source_key in (('coordinate', 'coordinate'), ('domain', 'domain'), ('radius', 'radius'))))
    mn = data[MOD]['exact_graph_nodes']
    check('same_actual_positive_L', port(18, MOD, [25272]) and wn[18]['source_positive'] is True
        and all(mn[physical['Q'][0]]['denominator'] == 25272 for physical in (lp, rp)))
    def radius_is_Z_independent(root):
        # Operation-specific references: rational numerators/denominators are
        # literal integers, never graph handles.
        nodes, seen, todo = data[LOOP]['exact_graph_nodes'], set(), [root]
        while todo:
            index = todo.pop()
            if index in seen:
                continue
            seen.add(index)
            row = nodes[index]
            operation = row['operation']
            if operation == 'exact_rational':
                continue
            if operation == 'current_original_source_parameter':
                if row.get('Z_independent') is not True:
                    return False
            elif operation == 'bound_variable':
                if index not in (lc['coordinate'], rc['coordinate']):
                    return False
            elif operation in ('sum', 'product'):
                todo.extend(row['arguments'])
            elif operation == 'negative':
                todo.append(row['argument'])
            elif operation == 'positive_quotient':
                if not row.get('source_positive_certificate'):
                    return False
                todo.extend((row['numerator'], row['denominator']))
            elif operation == 'analytic_unary' and row['name'] in ('exp', 'log', 'sqrt', 'positive_sqrt', 'log1p'):
                todo.append(row['argument'])
            else:
                return False
        return True

    check('same_regular_radius_function_has_no_Z_dependent_source',
        all(radius_is_Z_independent(chart['radius']) for chart in (lc, rc)))
    check('previous_IN_source_representation_bound', fg['outputs']['native_IN_representation'] == 45
        and fn[45]['admitted'] is True and fn[45]['eliminated_parent_jump_node'] == 520
        and fn[45]['arguments'][:2] == [26, 27]
        and fn[45]['arguments'][2:17] == list(range(11, 26))
        and fn[45]['arguments'][17:34] == list(range(28, 45))
        and fn[45]['arguments'][-1] == 10
        and all(fn[i]['source_node'] == i + 510 for i in range(11, 26))
        and [fn[i]['source_node'] for i in (26, 27)] == wn[537]['arguments'][:2]
        and [fn[i]['source_node'] for i in range(28, 45)] == wn[537]['arguments'][18:])
    retained = parent['retained_construction_data']
    check('all_old_construction_data_and_errors_frozen', retained == first['retained_construction_data']
        == weighted['retained_construction_data']
        and retained['cumulative_terms'] == 212 and retained['raw_signed_terms'] == 213
        and retained['reduced_terms'] == 183 and retained['six_signed_channels_retained'] is True
        and retained['all_old_error_data_retained'] is True
        and retained['reserved_contribution_registered_once'] is True
        and retained['inherited_errors_deleted'] is False)
    passed = all(checks.values())

    graph = Graph()
    theorem = graph.node('accepted_actual_endpoint_source_identity_port', source_record=PARENT,
        selector='component_result.endpoint_comparisons', canonicalizer_reexecuted=False)
    first_zero = graph.node('accepted_J0_zero_source_port', source_record=FIRST, source_node=10,
        weighted_source_node=520, ordinary_Z_orders=[0, 1])
    jump = graph.node('actual_J1_source_port', source_record=WEIGHTED, source_node=521,
        definition='B_frozen_macro_lower-B_second_micro_upper')
    bindings = [graph.node('accepted_actual_endpoint_function_row', (theorem,),
        label=label, source_record=archive, source_nodes=[ln, rn], ordinary_Z_order=order,
        canonical_digest=rows[label]['left']['canonical_digest'], admitted=passed)
        for label, (_, order, ln, rn, archive) in expected.items()]
    common_L = graph.node('actual_common_positive_L', source_record=MOD, source_node=25272,
        formula='1-delta*Z^2', ordinary_Z1_formula='-2*delta*Z')
    common_B = graph.node('actual_B_C1_composition_congruence', (*bindings, common_L),
        formula=B_FORMULA, ordinary_Z1_formula=B_Z_FORMULA,
        same_radius_function_differentiated_with_inherited_C1_regular_source=True,
        no_radius_Z_zero_assumption_required=True, admitted=passed)
    zero = graph.node('actual_J1_zero_C1', (jump, common_B), weighted_B_ports=[345, 355],
        same_source_difference_zero=passed, ordinary_Z_orders=[0, 1], owned_Z_cell=['0', '1/2'],
        numeric_zero_substituted_for_unknown_function=False, admitted=passed)
    remaining = [graph.node('unchanged_weighted_B_jump_source_port', source_record=WEIGHTED,
        source_node=i, jump_index=i - 520, source_identity_proved=False) for i in range(522, 536)]
    outer = [graph.node('unchanged_weighted_external_B_source_port', source_record=WEIGHTED,
        source_node=i) for i in wn[537]['arguments'][:2]]
    kernels = [graph.node('unchanged_weighted_K_integral_source_port', source_record=WEIGHTED,
        source_node=i) for i in wn[537]['arguments'][18:]]
    representation = graph.node('native_IN_with_J0_J1_source_zero_identities',
        (*outer, *remaining, *kernels, first_zero, zero), formula=IN_FORMULA if passed else fn[45]['definition'],
        weighted_native_S1_node=540, original_representation_node=537,
        eliminated_weighted_jump_nodes=[520, 521] if passed else [520], admitted=passed)
    result = dict(schema=SCHEMA, passed=passed, source_only=True,
        source_family=endpoint['source_family'], ordinary_Z_orders=[0, 1], ordinary_Z_not_temporal_n=True,
        source_prerequisite_audit=checks,
        source_interface=dict(parent_jump_node=521, parent_B_endpoints=[345, 355],
            source_seam=['second_micro.upper', 'frozen_macro.lower'],
            inherited_J0_zero_source_record=FIRST, inherited_J0_zero_source_node=10,
            accepted_endpoint_source_record=PARENT, accepted_endpoint_comparisons=list(expected),
            actual_L_node=25272, owned_Z_cell=['0', '1/2'], source_B_formula=B_FORMULA,
            J1_zero=passed, proved_weighted_jump_nodes=[520, 521] if passed else [520],
            remaining_parent_jump_nodes=list(range(522 if passed else 521, 536)),
            native_IN_formula=IN_FORMULA if passed else fn[45]['definition'], old_archived_jump_rows_unchanged=True,
            old_packet_errors_deleted=False),
        definition_graph=dict(namespace=SCHEMA, nodes=graph.nodes, node_count=len(graph.nodes),
            nodes_digest=digest(graph.nodes), outputs=dict(J1_zero=zero, B_C1_congruence=common_B,
                native_IN_representation=representation), proof_graph_not_numeric_evaluator=True),
        retained_construction_data=retained,
        unresolved_gates=[name for name, value in checks.items() if not value],
        remaining_construction_gates=['fourteen_other_native_B_jump_source_bindings',
            'native_old_active_density_source_comparison', 'compact_density_source_comparison',
            'same_source_old_numerical_packet_error_theorem', 'global_five_moment_closure'],
        runtime_model_metadata=None,
        source_inspection_audit=dict(new_source_consumer_calls=1, accepted_ancestor_operator_calls=0,
            canonicalizer_reexecuted=False, old_operators_replayed=False, numeric_values_computed=False))
    result.update({name: passed for name in STAGE_FLAGS})
    result.update({name: False for name in FALSE_FLAGS})
    return result
