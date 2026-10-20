"""Lift actual Rc/2Rc physical pressure traces to ordinary C3 source joins.

Saved source theorems are premises. This leaf constructs a new physical flux
jet from actual full-m ports and discharges the saved natural jump functions.
No archived producer, graph evaluator or source-value reducer is executed.
"""
from __future__ import annotations

import ast
from math import comb
import hashlib
import json

PREFIX = 'lei_ren_part1_paper_compliant_current_original_'
PARENT_REPORT = PREFIX + 'native_Rc_Rp_signed_pressure_source.json'
OLD_JOIN = PREFIX + 'Rc_pressure_join.json'
OLD_TAIL = PREFIX + 'tail_m15_numeric.json'
QUIET = PREFIX + 'true_quiet_Rp_C5_pressure_C3_functions.json.gz'
OUTER_CHECK = PREFIX + 'true_outer_C6_C5_functions_check.json'
COMPACT_CHECK = PREFIX + 'true_compact_C5_functions_check.json'
QUIET_CHECK = PREFIX + 'true_quiet_Rp_C5_pressure_C3_functions_check.json'
PHYSICAL_SOURCE = PREFIX + 'modified_Q4_pressure_C3_functions.py'
COMPACT_SOURCE = PREFIX + 'true_compact_C5_functions.py'
QUIET_SOURCE = PREFIX + 'true_quiet_Rp_C5_pressure_C3_functions.py'
SCHEMA = 'original_native_Rc_2Rc_pressure_join_lift_v1'
RECORD_HASHES = {
    PARENT_REPORT: 'a8d76d71b7d527c65c4eb29cc4e8340fc6592723246a87df65e496d3bce43a59',
    OLD_JOIN: '5f160cfa1f251253155bdad0466c462c74c5d7349bdf5cbcb9b3b7209822e9c4',
    OLD_TAIL: 'a394d3b538e420936e382b36b63dc7405b7a098b57a557c668deba7a2cc60c9a',
    QUIET: 'bbc8ae0361bdfbb75255707e7a234eb39f531533e6f39ad718cf262e27ad42c4',
    OUTER_CHECK: 'ce5669d01621683c3be857917f206b50af067b44feefe4b000c2ba99b2785bba',
    COMPACT_CHECK: '6c933d077f9c151eda246acd4d5f694d32dc8272fc3fb65a031b9b9cf28c8e7b',
    QUIET_CHECK: '3cb8b2e7c48ee62ce2b3dba9c6b0d68b901010d6a154b18d466d3221e2c8e7c0',
    PHYSICAL_SOURCE: '7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7',
    COMPACT_SOURCE: 'f0b0333b19a0e2504d023b97e841033d3cbcc23427c57ad9672264866c01ec7b',
    QUIET_SOURCE: '0db6ea1a97184662c8e80f011844c12d143eeeb11333510f2fbbb8dcc19de068',
}
DECODED_RECORD_HASHES = {QUIET: '3bc0f9984a0aeb6f447d177c40c25116b8a59c3b90ee09b8661f029cbe032a40'}
STAGE_FLAGS = (
    'native_Rc_pressure_join_C3_proved', 'native_2Rc_pressure_join_C3_proved',
    'native_Rc_Rp_natural_pressure_composition_proved',
    'natural_physical_join_defects_zero_proved',
)
FALSE_FLAGS = (
    'same_density_certificate', 'finite_N_history_identity', 'physical_Q_y_Z1',
    'global_five_moment_closure', 'legacy_feedback_reference_alias_proved',
    'prior_adapter_geometry_alias_proved', 'prior_adapter_P0_alias_proved',
    'time_order_recursion_claimed', 'numeric_values_computed', 'old_operators_replayed',
    'full_goal_complete', 'old212_adapter_source_identity_proved',
)
GRAPH_DIGEST = '15d51cd83dd82f41eba990979ace934bd5c8d5c339cd516b9fc96bc20a64af56'
OUTER_DIGEST = '7ec4add5799967f01186af73b98d1c64ef45462b8a2c13ace2e24fc7ca3e7ad5'
COMPACT_DIGEST = '13f106b857ad6e709065810cade2d0dbbf05d62a6619d7165bfc9624bad1e8c3'


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode()).hexdigest()


class SourceJetGraph:
    """New ordinary derivative expressions; no method evaluates source atoms."""

    def __init__(self):
        self.nodes, self.cache = [], {}
        self.zero = self.node('exact_rational', [], numerator=0, denominator=1)
        self.one = self.node('exact_rational', [], numerator=1, denominator=1)

    def node(self, operation, arguments, **data):
        row = dict(operation=operation, arguments=arguments, **data)
        key = digest(row)
        if key not in self.cache:
            ref = SCHEMA + ':j' + str(len(self.nodes))
            self.cache[key] = ref
            self.nodes.append(dict(id=ref, **row))
        return self.cache[key]

    def constant(self, value):
        return self.node('exact_rational', [], numerator=value, denominator=1)

    def source(self, handle, order=0):
        return self.node('actual_source_function_port', [], source_graph_digest=GRAPH_DIGEST,
            source_node=handle, ordinary_Z_order=order)

    def add(self, *args):
        return self.node('sum', list(args))

    def scale(self, value, coefficient):
        return self.node('product', [self.constant(coefficient), value])

    def mul(self, *args):
        return self.node('product', list(args))

    def fixed(self, value, order):
        return [value] + [self.zero] * order

    def jadd(self, left, right):
        return [self.add(a, b) for a, b in zip(left, right)]

    def jscale(self, values, coefficient):
        return [self.scale(v, coefficient) for v in values]

    def jmul(self, left, right):
        return [self.add(*(self.scale(self.mul(left[j], right[k-j]), comb(k, j))
            for j in range(k+1))) for k in range(len(left))]

    def jdiv(self, numerator, denominator):
        out = []
        for k in range(len(numerator)):
            known = self.add(*(self.scale(self.mul(denominator[j], out[k-j]), comb(k, j)) for j in range(1, k+1)))
            rhs = numerator[k] if k == 0 else self.add(numerator[k], self.scale(known, -1))
            out.append(self.node('positive_quotient', [rhs, denominator[0]],
                positive_guard='original common L=1-delta*Z^2>0', ordinary_Z_order=k))
        return out

    def common_flux(self, mean_rows, radius, pstar, delta):
        z = [self.source(20), self.one, self.zero, self.zero]
        one, de = self.fixed(self.one, 3), self.fixed(self.source(delta), 3)
        zz = self.jmul(z, z)
        D = self.jadd(one, self.jscale(zz, -1))
        L = self.jadd(one, self.jscale(self.jmul(de, zz), -1))
        physical = [self.mul(self.source(pstar), self.source(ref, k)) for k, ref in enumerate(mean_rows[:5])]
        m, mz = physical[:4], physical[1:5]
        numerator = self.jadd(self.jmul(self.jmul(self.jadd(one, self.jscale(de, -1)), z), m), self.jmul(D, mz))
        qy = self.jdiv(numerator, L)
        flux = [self.mul(self.source(radius), value) for value in qy]
        return dict(Q=self.jscale(qy, -1), Q_y=qy, Q_plus_Q_y=[self.zero]*4,
            R_Q_y=flux, physical_mean_rows=physical, denominator=L,
            formula='RQ_y=R*Pstar*((1-delta)*Z*m+(1-Z^2)*D_Z(m))/(1-delta*Z^2)',
            common_flux_individually_zero_claim=False)


def ast_function(source_text, name):
    return next(node for node in ast.walk(ast.parse(source_text)) if isinstance(node, ast.FunctionDef) and node.name == name)


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
    data = {n: json.loads(v['source_text']) for n, v in source_records.items() if not n.endswith('.py')}
    parent, old, quiet = data[PARENT_REPORT], data[OLD_JOIN], data[QUIET]
    accepted = parent['component_result']
    check('accepted_pressure_cocycle', accepted['passed'] is True and accepted['signed_pressure_cocycle_proved'] is True
        and accepted['natural_physical_join_defects_zero_proved'] is False)
    nodes = quiet['exact_graph_nodes']
    check('actual_native_master_graph', len(nodes) == 265856 and digest(nodes) == GRAPH_DIGEST
        and digest(nodes[:258857]) == OUTER_DIGEST and digest(nodes[:263344]) == COMPACT_DIGEST)
    interface = accepted['function_source_interface']
    cb, qb = interface['compact']['source'], interface['quiet']['source']
    rc, ct = interface['native_Rc'], cb['terminal']
    check('actual_quiet_binding', qb == quiet['true_quiet_Rp_C5_C3_bindings'])
    check('same_source_family', old['source_family'] == quiet['source_family'] == accepted['source_family'])
    profile, join = old['Rc_profile_join'], old['Rc_pressure_join']
    ports = join['actual_physical_pressure_ports']
    check('old_join_source_proof_not_fixture', old['checks']['passed'] is True and join['passed'] is True
        and join['source_executed'] is False and join['input_source_prefix_sha256'] == OUTER_DIGEST
        and ports['passed'] is True and ports['source_only'] is True and ports['source_executed'] is False
        and all(ports['checks'].values()) and join['bounds_or_midpoints_registered_as_values'] is False)
    check('qualified_profile_zero_theorem', profile['passed'] is True and profile['source_only'] is True
        and profile['source_executed'] is False and profile['numeric_caps_used'] is False
        and profile['input_source_prefix_sha256'] == OUTER_DIGEST and profile['full_V_and_V_y_exact_zero'] is True
        and profile['compact_upper_full_V_and_V_y_exact_zero'] is True
        and profile['compact_upper_profile_V_and_V_y_exact_zero'] is True)
    check('Rc_actual_trace_ports', join['endpoint_interfaces']['outer']['boundary'] == rc['boundary']
        and join['endpoint_interfaces']['compact']['boundary'] == cb['inlet']['boundary']
        and join['endpoint_interfaces']['outer']['physical'] == rc['physical'])
    check('saved_profile_attached_to_same_Rc_proof', join['required_profile_proof_digest'] == digest(profile)
        and profile['actual_functions']['outer_O3_power_V'] == rc['modified_fields']['V']
        and profile['actual_functions']['outer_O3_power_V_y'] == rc['radial']['V_y']
        and profile['actual_functions']['compact_inlet_V'] == cb['inlet']['modified_fields']['V']
        and profile['actual_functions']['compact_inlet_V_y'] == cb['inlet']['radial']['V_y'])
    check('Rc_actual_mean_base_identity', join['full_mean_identities'][0]['passed'] is True
        and join['full_mean_identities'][0]['outer_full_mean'] == rc['modified_fields']['m'][0]
        and join['full_mean_identities'][0]['compact_full_mean'] == cb['inlet']['modified_fields']['m'][0]
        and join['deterministic_recovery_congruence']['Q_Qy_ordinary_C1_equal'] is True
        and join['deterministic_recovery_congruence']['radius_same_source'] is True)
    tail = data[OLD_TAIL]
    tail_seams = tail['source_seam_proofs']
    tail_join = tail_seams['joins']['compact_quiet']
    check('old_2Rc_source_calculus_qualified', tail_seams['passed'] is True and tail_join['passed'] is True
        and tail_seams['source_values_selected'] is False and tail_seams['numeric_caps_used'] is False
        and tail_seams['ancestor_replayed'] is False and tail['source_family'] == quiet['source_family']
        and tail['input_hashes'][QUIET] == RECORD_HASHES[QUIET]
        and tail_join['left_coordinate'] == 260181 and tail_join['right_coordinate'] == qb['s']
        and tail_join['left_endpoint_node'] == 48 and tail_join['right_endpoint_node'] == 0)
    for row in tail_join['trace_rows']:
        quantity, order = row['quantity'], row['ordinary_Z_order']
        if quantity in ('Q', 'Qy'):
            field = 'Q' if quantity == 'Q' else 'Q_y'
            check('old_2Rc_actual_pressure_port:' + quantity + str(order), row['passed'] is True
                and row['left']['recipe']['expression'] == cb['physical'][field][order]
                and row['right']['recipe']['expression'] == qb['physical'][field][order]
                and row['right']['node'] == qb['inlet']['physical'][field][order]
                and row['physical_formula']['sha256'] == RECORD_HASHES[PHYSICAL_SOURCE])
        if quantity == 'RQy':
            check('old_2Rc_actual_RQy_port:' + str(order), row['passed'] is True
                and row['left']['node'] == ct['boundary']['R_Q_y'][order]
                and row['right']['node'] == qb['inlet']['boundary']['R_Q_y'][order]
                and row['common_radius_and_Qy_source_deduction'] is True)
    check('same_physical_parameters', cb['N'] == qb['N'] == 18 and cb['Pstar'] == qb['Pstar'] == 33572
        and cb['delta'] == qb['delta'] == 25220 and cb['P0'] == qb['P0'])
    check('exact_compact_to_quiet_owner', qb['incoming'] == ct and qb['R2'] == ct['radius'])
    for name in (OUTER_CHECK, COMPACT_CHECK, QUIET_CHECK):
        receipt = data[name]
        check('genuine_source_regularity:' + name, receipt['all_passed'] is True
            and receipt['source_family'] == quiet['source_family'] and all(v['passed'] is True for v in receipt['checks'].values()))
    check('compact_C5_pressure_C3', data[COMPACT_CHECK]['current_original_true_rminus_to_2Rc_signed_pressure_C3_functions_installed'] is True)
    check('quiet_C5_pressure_C3', data[QUIET_CHECK]['current_original_true_rminus_to_Rp_signed_pressure_C3_functions_installed'] is True)
    physical_ast = ast_function(source_records[PHYSICAL_SOURCE]['source_text'], 'physical_pressure_fields')
    physical_text = ast.unparse(physical_ast)
    check('same_actual_physical_recovery_formula', all(term in physical_text for term in (
        "U = a.scale(fields['V'], Pstar)", "Uy = a.scale(radial['V_y'], Pstar)",
        "m = a.scale(fields['m'], Pstar)", 'my = a.sub(U, m)',
        'Q = div(numerator(U, m))', 'Qy = div(numerator(Uy, my))')))
    quiet_flow = ast.unparse(ast_function(source_records[QUIET_SOURCE]['source_text'], 'power_flow'))
    check('actual_quiet_mean_flow_recipe', "m=a.scale(a.trim(own['m']), decay(1))" in quiet_flow
        and 'zero = a.fixed(0)' in quiet_flow and 'V=zero' in quiet_flow)
    check('quiet_zero_V_Vy_functions', qb['fields']['V'] == qb['radial']['V_y'] == [0]*6)
    mean_flow = qb['fields']['m']
    factors = []
    for order, (ref, incoming) in enumerate(zip(mean_flow, ct['modified_fields']['m'])):
        row = nodes[ref]
        valid = row.get('operation') == 'product' and incoming in row.get('arguments', []) and len(row.get('arguments', [])) == 2
        remaining = [arg for arg in row.get('arguments', []) if arg != incoming]
        check('actual_quiet_full_mean_factor_Z' + str(order), valid and len(remaining) == 1)
        factors.extend(remaining)
    decay = factors[0]
    check('same_Z_independent_quiet_decay', len(factors) == 6 and len(set(factors)) == 1
        and nodes[decay].get('operation') == 'analytic_unary' and nodes[decay].get('name') == 'exp'
        and nodes[nodes[decay]['argument']] == dict(operation='negative', argument=qb['s']))
    check('quiet_inlet_actual_mean_substitutions', all(nodes[target] == dict(operation='function_substitution',
        expression=source, variable=qb['s'], value=0, Z_independent_substitution=True)
        for source, target in zip(mean_flow, qb['inlet']['fields']['m'])))
    check('quiet_actual_radius_recipe', nodes[qb['R']]['operation'] == 'product' and qb['R2'] in nodes[qb['R']]['arguments']
        and all(nodes[arg] == dict(operation='analytic_unary', name='exp', argument=qb['s'])
            for arg in nodes[qb['R']]['arguments'] if arg != qb['R2']))
    check('quiet_inlet_radius_substitution', nodes[qb['inlet']['radius']] == dict(operation='function_substitution',
        expression=qb['R'], variable=qb['s'], value=0, Z_independent_substitution=True))
    check('actual_compact_upper_profile_ports', all(nodes[target] == dict(operation='function_substitution',
        expression=source, variable=260181, value=48, Z_independent_substitution=True)
        for key, source_group, endpoint_group in (
            ('V', cb['modified_fields'], ct['modified_fields']), ('V_y', cb['radial'], ct['radial']))
        for source, target in zip(source_group[key], endpoint_group[key])))
    check('old_upper_profile_identifies_actual_profiles', profile['compact_flat_bump_nodes'] == cb['profiles']['bumps']
        and profile['actual_functions']['compact_profile_F_G_V']['0']['V'] == cb['profiles']['radial_y_rows'][0]['V']
        and profile['actual_functions']['compact_profile_F_G_V']['1']['V'] == cb['profiles']['radial_y_rows'][1]['V'])
    check('compact_upper_profiles_are_actual_physical_inputs',
        cb['modified_fields']['V'] == cb['profiles']['radial_y_rows'][0]['V']
        and cb['radial']['V_y'] == cb['profiles']['radial_y_rows'][1]['V'])
    compact_ast = ast.unparse(ast.parse(source_records[COMPACT_SOURCE]['source_text']))
    quiet_ast = ast.unparse(ast.parse(source_records[QUIET_SOURCE]['source_text']))
    check('actual_physical_recovery_calls',
        'self.physical = radial.physical_pressure_fields(g, self.modified_fields, self.radial, self.delta, self.Pstar)' in compact_ast
        and 'self.physical = radial.physical_pressure_fields(g, self.fields, self.radial, self.delta, self.Pstar)' in quiet_ast)
    check('ordinary_jet_widths', all(len(endpoint['modified_fields']['m']) == 6 for endpoint in (rc, ct))
        and len(qb['inlet']['fields']['m']) == 6)
    graph = SourceJetGraph()
    proof = graph.node('accepted_source_theorems', [], old_report_sha256=RECORD_HASHES[OLD_JOIN], old_2Rc_report_sha256=RECORD_HASHES[OLD_TAIL],
        Rc_C0_C1_actual_function_identity=True, compact_upper_actual_V_Vy_zero=True,
        original_source_domain_only=True, old_operator_called=False)
    source_domains = dict(common_original_positive_L_domain=True, genuine_saved_C5_histories=True,
        genuine_saved_pressure_C3=True, Z_independent_physical_scale_and_radial_geometry=True,
        numeric_high_order_norm_admission=False, extension_outside_original_domain=False)
    regularity = graph.node('accepted_genuine_source_regularity', [proof], **source_domains)
    outputs, zero_rows = {}, []
    passed = all(checks.values())
    for label, before, after, radius, mean in (
        ('Rc', rc, cb['inlet'], cb['Rc'], rc['modified_fields']['m']),
        ('2Rc', ct, qb['inlet'], ct['radius'], ct['modified_fields']['m']),
    ):
        common = graph.common_flux(mean, radius, cb['Pstar'], cb['delta'])
        base = graph.node('actual_primitive_function_congruence', [proof, regularity], seam=label,
            full_mean_left=mean[0], full_mean_right=(after['modified_fields'] if label == 'Rc' else after['fields'])['m'][0],
            actual_V_Vy_zero_on_both_sides=True, ordinary_Z_identity_for_all_original_Z=True,
            mean_identity_source=('saved_exact_Rc_same_history_identity' if label == 'Rc' else 'qualified_saved_2Rc_source_calculus_and_actual_quiet_full_mean_factor'),
            scalar_identity=('inherited_same_source_sum' if label == 'Rc' else 'exp(0)=1'),
            actual_common_radius_source=radius, admitted=passed)
        lift = graph.node('differentiate_actual_function_identity', [base, regularity], seam=label,
            primitive_mean_orders=[0,1,2,3,4], primitive_velocity_orders=[0,1,2,3],
            pressure_orders=[0,1,2,3], all_Z_base_identity_required=True,
            higher_cross_side_values_assumed=False, admitted=passed)
        rows = {}
        for moment, boundary in (('m13', 'Q_plus_Q_y'), ('m15', 'R_Q_y')):
            rows[moment] = []
            parent_rows = accepted['definition_graph']['outputs'][moment]
            for order in range(4):
                left, right = before['boundary'][boundary][order], after['boundary'][boundary][order]
                left_port, right_port = graph.source(left, order), graph.source(right, order)
                common_ref = common[boundary][order]
                edge = graph.node('actual_physical_boundary_congruence', [lift, left_port, right_port, common_ref],
                    seam=label, ordinary_Z_order=order, actual_left_handle=left, actual_right_handle=right,
                    common_expression=common_ref, physical_Q_not_loop_q=True, admitted=passed)
                old_jump = parent_rows[order]['natural_join_defects'][0 if label == 'Rc' else 1]
                zero = graph.node('natural_pressure_jump_exact_zero', [edge, graph.zero], seam=label,
                    moment=moment, ordinary_Z_order=order, parent_jump=old_jump,
                    source_identity='actual_right_boundary-actual_left_boundary', admitted=passed,
                    common_flux_individually_zero_claim=(boundary == 'Q_plus_Q_y'))
                rows[moment].append(zero)
                zero_rows.append(dict(seam=label, moment=moment, ordinary_Z_order=order, actual_left_handle=left,
                    actual_right_handle=right, parent_jump=old_jump, exact_zero=zero, admitted=passed))
        outputs[label] = dict(common_physical_flux_C3=common, primitive_identity=base, derivative_lift=lift,
            pressure_join_zero_rows=rows, common_flux_not_a_bound_or_midpoint=True)
    composition = []
    for moment, rows in accepted['definition_graph']['outputs'].items():
        for order, row in enumerate(rows):
            composition.append(graph.node('natural_signed_pressure_composition',
                [outputs['Rc']['pressure_join_zero_rows'][moment][order], outputs['2Rc']['pressure_join_zero_rows'][moment][order]],
                moment=moment, ordinary_Z_order=order, parent_source_report_sha256=RECORD_HASHES[PARENT_REPORT],
                parent_natural_segments=row['natural_segments'], parent_telescoped=row['telescoped'],
                parent_absolute_Rp=row['absolute_Rp'], parent_jump_functions=row['natural_join_defects'],
                admitted=passed, old212_error_deletion_authorized=False))
    result = dict(schema=SCHEMA, passed=passed, source_only=True, ordinary_Z_orders=[0,1,2,3],
        ordinary_Z_not_temporal_n=True, input_digest=expected_input_digest,
        pressure_join_zero_scope='Only actual native Rc and 2Rc pressure traces; sixteen scalar ordinary-Z jump functions',
        native_C2_C5_transfer_orders=[0,1,2], pressure_C3_origin='pinned genuine C5/pressure-C3 source functions',
        source_family=quiet['source_family'], source_prerequisite_audit=checks,
        source_witness=dict(master_graph_digest=GRAPH_DIGEST, old_Rc_proof_digest=digest(join),
            old_profile_proof_digest=digest(profile), old_2Rc_source_calculus_digest=digest(tail_join), physical_recovery_AST_digest=digest(ast.dump(physical_ast)),
            quiet_full_mean_decay_handle=decay, source_domains=source_domains),
        definition_graph=dict(namespace=SCHEMA, nodes=graph.nodes, node_count=len(graph.nodes), nodes_digest=digest(graph.nodes),
            outputs=outputs, natural_pressure_composition=composition, proof_graph_not_numeric_evaluator=True),
        actual_jump_zero_rows=zero_rows, retained_construction_data=parent['retained_construction_data'],
        source_inspection_audit=dict(new_source_consumer_calls=1, accepted_ancestor_operator_calls=0,
            old_definition_graphs_rebuilt=False, old_operators_replayed=False, numeric_values_computed=False,
            quadrature_or_norm_calls=0, source_value_reducer_calls=0, temporal_recursion_calls=0),
        open_obligations=['old212 same-source numerical packet binding', 'native Rp independent complete pulse and once-only reserved match',
            'other physical seams and global five moment closure', 'exact heat/stress/temporal recursion/corrections/Cartesian NS'],
        unresolved_gates=[key for key,value in checks.items() if not value], runtime_model_metadata=None)
    result.update({flag: passed for flag in STAGE_FLAGS})
    result.update({flag: False for flag in FALSE_FLAGS})
    return result
