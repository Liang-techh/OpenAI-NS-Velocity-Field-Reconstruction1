"""Native Rc-to-Rp pressure functions with explicit signed natural join defects.

Only authenticated saved definitions are consumed. The new affine graph is a
function expression interface, not an evaluator for the archived source graph.
"""
from __future__ import annotations

from fractions import Fraction
import ast
import hashlib
import json

PREFIX = 'lei_ren_part1_paper_compliant_current_original_'
PARENT_REPORT = PREFIX + 'native_Rm_m_Z_source_lift.json'
OUTER_REPORT = PREFIX + 'native_C2_outer_six_source.json'
OUTER = PREFIX + 'true_outer_C6_C5_functions.json.gz'
COMPACT = PREFIX + 'true_compact_C5_functions.json.gz'
QUIET = PREFIX + 'true_quiet_Rp_C5_pressure_C3_functions.json.gz'
OUTER_CHECK = PREFIX + 'true_outer_C6_C5_functions_check.json'
COMPACT_CHECK = PREFIX + 'true_compact_C5_functions_check.json'
QUIET_CHECK = PREFIX + 'true_quiet_Rp_C5_pressure_C3_functions_check.json'
OUTER_SOURCE = PREFIX + 'true_outer_C6_C5_functions.py'
COMPACT_SOURCE = PREFIX + 'true_compact_C5_functions.py'
QUIET_SOURCE = PREFIX + 'true_quiet_Rp_C5_pressure_C3_functions.py'
SCHEMA = 'original_native_Rc_Rp_signed_pressure_source_v1'
RECORD_HASHES = {
    PARENT_REPORT: '2e2654a8be84f2f92ae149f3ed09ad3ef85679180024626138aa5f613211e254',
    OUTER_REPORT: 'ac4900e6282dbb339708fcb202d9f21eea959f24e27d682f6f17521610216459',
    OUTER: '75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d',
    COMPACT: '3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269',
    QUIET: 'bbc8ae0361bdfbb75255707e7a234eb39f531533e6f39ad718cf262e27ad42c4',
    OUTER_CHECK: 'ce5669d01621683c3be857917f206b50af067b44feefe4b000c2ba99b2785bba',
    COMPACT_CHECK: '6c933d077f9c151eda246acd4d5f694d32dc8272fc3fb65a031b9b9cf28c8e7b',
    QUIET_CHECK: '3cb8b2e7c48ee62ce2b3dba9c6b0d68b901010d6a154b18d466d3221e2c8e7c0',
    OUTER_SOURCE: '4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3',
    COMPACT_SOURCE: 'f0b0333b19a0e2504d023b97e841033d3cbcc23427c57ad9672264866c01ec7b',
    QUIET_SOURCE: '0db6ea1a97184662c8e80f011844c12d143eeeb11333510f2fbbb8dcc19de068',
}
DECODED_RECORD_HASHES = {
    OUTER: '53b1bbf15c8c2c5a8400d7e22ef940c6a502397df8cf26bb23e51ff3869e9139',
    COMPACT: 'b6a7227da5e2427561cd22fcff614003cc14f3182537e1dbe278a6de84322510',
    QUIET: '3bc0f9984a0aeb6f447d177c40c25116b8a59c3b90ee09b8661f029cbe032a40',
}
STAGE_FLAGS = (
    'native_Rc_incoming_source_bound', 'native_compact_quiet_source_continuation_proved',
    'signed_pressure_cocycle_proved', 'same_source_internal_boundary_cancellation_proved',
)
FALSE_FLAGS = (
    'same_density_certificate', 'finite_N_history_identity', 'physical_Q_y_Z1',
    'global_five_moment_closure', 'legacy_feedback_reference_alias_proved',
    'prior_adapter_geometry_alias_proved', 'prior_adapter_P0_alias_proved',
    'time_order_recursion_claimed', 'numeric_values_computed', 'old_operators_replayed',
    'full_goal_complete', 'old212_adapter_source_identity_proved',
    'natural_physical_join_defects_zero_proved',
)
GRAPH_HASHES = (
    '7ec4add5799967f01186af73b98d1c64ef45462b8a2c13ace2e24fc7ca3e7ad5',
    '13f106b857ad6e709065810cade2d0dbbf05d62a6619d7165bfc9624bad1e8c3',
    '15d51cd83dd82f41eba990979ace934bd5c8d5c339cd516b9fc96bc20a64af56',
)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode()).hexdigest()


class AffineFunctions:
    """Exact rational coefficients of NEW expressions; source ports stay atoms."""

    def __init__(self):
        self.nodes, self.forms, self.ports = [], {}, {}

    def append(self, operation, arguments, form=None, **data):
        ref = SCHEMA + ':f' + str(len(self.nodes))
        self.nodes.append(dict(id=ref, operation=operation, arguments=arguments, **data))
        self.forms[ref] = {ref: Fraction(1)} if form is None else form
        return ref

    def port(self, node, **context):
        key = digest(dict(node=node, **context))
        if key not in self.ports:
            self.ports[key] = self.append('actual_source_function_port', [],
                source_graph_digest=GRAPH_HASHES[2], source_node=node, **context)
        return self.ports[key]

    def at(self, node, variable, value, order):
        source = self.port(node, ordinary_Z_order=order)
        return self.append('unevaluated_Z_independent_endpoint_substitution', [source],
            source_node=node, variable=variable, value=value, ordinary_Z_order=order,
            coordinate_absent_from_incoming_prefix=True)

    def combination(self, terms):
        form = {}
        for coefficient, ref in terms:
            for atom, value in self.forms[ref].items():
                form[atom] = form.get(atom, Fraction(0)) + Fraction(coefficient) * value
        form = {k: v for k, v in form.items() if v}
        return self.append('exact_signed_affine_function', [r for _, r in terms], form,
            coefficients=[str(Fraction(c)) for c, _ in terms])

    def equal(self, left, right):
        return self.forms[left] == self.forms[right]

    def normal_form(self, ref):
        return {k: str(v) for k, v in sorted(self.forms[ref].items())}


def _sum(nodes, ref, left, right):
    return nodes[ref] == dict(operation='sum', arguments=[left, right])


def _negative(nodes, ref, source):
    return nodes[ref] == dict(operation='negative', argument=source)


def _scale(nodes, ref, source, coefficient):
    row = nodes[ref]
    if row.get('operation') != 'product' or len(row.get('arguments', [])) != 2:
        return False
    first, second = row['arguments']
    other = second if first == source else first if second == source else None
    if other is None or nodes[other].get('operation') != 'exact_rational':
        return False
    return Fraction(nodes[other]['numerator'], nodes[other]['denominator']) == coefficient


def _pressure_constructors(nodes, owner, local_key, cumulative_key):
    """Inspect only defining constructor metadata, never reduce a source value."""
    out = {}
    for moment, integral, boundary, coefficient, sign in (
        ('m13', 'S0', 'Q_plus_Q_y', Fraction(-1, 2), 1),
        ('m15', 'S1', 'R_Q_y', Fraction(1, 2), -1),
    ):
        for order in range(4):
            cumulative = owner[cumulative_key][integral][order]
            local = owner[local_key][integral][order]
            incoming = owner['incoming_integrals'][integral][order]
            pressure = nodes[owner['pressure'][moment][order]]
            valid = _sum(nodes, cumulative, incoming, local)
            valid = valid and pressure.get('operation') == 'sum' and len(pressure.get('arguments', [])) == 2
            if valid:
                scaled, signed_boundary = pressure['arguments']
                valid = _scale(nodes, scaled, cumulative, coefficient)
                difference = signed_boundary
                if sign == -1:
                    valid = valid and nodes[signed_boundary].get('operation') == 'negative'
                    difference = nodes[signed_boundary].get('argument', 0)
                row = nodes[difference]
                valid = valid and row.get('operation') == 'sum' and len(row.get('arguments', [])) == 2
                if valid:
                    positive, negative = row['arguments']
                    valid = positive == owner['boundary'][boundary][order] and _negative(
                        nodes, negative, owner['initial'][boundary][order])
            out[moment + '_Z' + str(order)] = bool(valid)
    return out


def _endpoint(nodes, group, endpoint, variable, value):
    return all(nodes[target] == dict(operation='function_substitution', expression=source,
        variable=variable, value=value, Z_independent_substitution=True)
        for key in ('integrals', 'pressure', 'boundary')
        for source, target in zip(group[key].values(), endpoint[key].values())
        for source, target in zip(source, target))


def derive(source_records, *, expected_input_digest):
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)

    check('input_digest', digest(source_records) == expected_input_digest)
    for name, expected in RECORD_HASHES.items():
        record = source_records[name]
        check('authenticated:' + name, record['raw_authenticated_before_decode'] is True
            and record['raw_sha256'] == expected
            and record['decoded_sha256'] == DECODED_RECORD_HASHES.get(name, expected)
            and hashlib.sha256(record['source_text'].encode()).hexdigest() == record['decoded_sha256'])
        if name.endswith('.py'):
            ast.parse(record['source_text'])
    data = {name: json.loads(record['source_text']) for name, record in source_records.items()
        if not name.endswith('.py')}
    outer, compact, quiet = [data[name] for name in (OUTER, COMPACT, QUIET)]
    graphs = [d['exact_graph_nodes'] for d in (outer, compact, quiet)]
    check('exact_graph_counts', [len(g) for g in graphs] == [258857, 263344, 265856])
    check('canonical_graph_hashes', tuple(digest(g) for g in graphs) == GRAPH_HASHES)
    check('exact_outer_compact_prefix', graphs[1][:258857] == graphs[0])
    check('exact_compact_quiet_prefix', graphs[2][:263344] == graphs[1])
    family = outer['source_family']
    check('same_original_source_family', compact['source_family'] == quiet['source_family'] == family)
    ob, cb, qb = outer['true_outer_C6_C5_bindings'], compact['true_compact_C5_bindings'], quiet['true_quiet_Rp_C5_C3_bindings']
    rc, bg, ct, rp = ob['endpoints']['O3_power'], ob['background_endpoints']['O3_power'], cb['terminal'], qb['Rp_endpoint']
    parent, native = data[PARENT_REPORT], data[OUTER_REPORT]['component_result']
    check('accepted_local_Rm_source_parent', parent['component_result']['passed'] is True
        and parent['component_result']['Rm_physical_Qy_Z1_function_identity_proved'] is True)
    check('accepted_native_outer_source', native['passed'] is True and native['native_outer_six_source_identity'] is True
        and native['native_outer_six_nonreset_history_identity_proved'] is True and native['terminal_endpoint'] == 'Rc'
        and native['source_witness']['accepted_reports'][OUTER]['raw_sha256'] == RECORD_HASHES[OUTER])
    for archive, receipt in ((OUTER, OUTER_CHECK), (COMPACT, COMPACT_CHECK), (QUIET, QUIET_CHECK)):
        saved = data[receipt]
        check('accepted_source_receipt:' + archive, saved['all_passed'] is True and saved['source_family'] == family
            and saved['input_hashes'][archive] == RECORD_HASHES[archive]
            and all(v['passed'] is True for v in saved['checks'].values()))
    check('same_native_N_mu', cb['N'] == qb['N'] == ob['N'] == 18 and cb['mu'] == qb['mu'] == ob['mu'] == 17)
    check('same_physical_Pstar_P0_delta', cb['Pstar'] == qb['Pstar'] == 33572 and cb['P0'] == qb['P0']
        and cb['delta'] == qb['delta'] == 25220)
    check('same_initial_boundary', ob['initial'] == cb['initial'] == qb['initial'])
    check('actual_Rc_correction_incoming', cb['incoming_correction'] == rc['correction'])
    check('actual_Rc_background_incoming', cb['incoming_background'] == bg['history'] and cb['A'] == bg['fields']['E'])
    check('actual_Rc_integral_incoming', cb['incoming_integrals'] == rc['integrals'] and cb['Rc'] == rc['radius'])
    check('quiet_exact_compact_exit', qb['incoming'] == ct and qb['relative_incoming'] == ct['correction']
        and qb['leading_incoming'] == ct['background_fields'] and qb['incoming_integrals'] == ct['integrals']
        and qb['incoming_pressure'] == ct['pressure'] and qb['R2'] == ct['radius'])
    check('no_second_N_division', cb['definition']['actual_selected_N_Rc_correction_used_without_additional_N_division'] is True
        and cb['branch_transfer']['exact_source_correspondence']['physical_selected_N_incoming_not_divided_again'] is True)
    check('same_genuine_compact_branch', cb['branch_transfer']['same_original_B_Q_N_and_exact_bump_integrals'] is True
        and cb['branch_transfer']['original_C1_limit_proof']['contraction_at_most'] == '1/2')
    check('relative_only_terminal_identity', all(graphs[2][ref]['relative_correction_only'] is True
        and graphs[2][ref]['absolute_leading_histories_and_pressure_not_zeroed'] is True
        and graphs[2][ref]['actual_endpoint_rows'] == ct['correction'][key]
        for key, ref in cb['terminal_certificates'].items()))
    check('quiet_unresolved_pulse_retained', qb['definition']['selected_native_pulse_C5_consumes_this_Rp_frame'] is False)
    nodes = graphs[2]
    oc = dict(ob['charts']['O3_power'], initial=ob['initial'])
    constructors = {
        'outer': _pressure_constructors(nodes, oc, 'local_integrals', 'cumulative_integrals'),
        'compact': _pressure_constructors(nodes, cb, 'pressure_local', 'integrals'),
        'quiet': _pressure_constructors(nodes, qb, 'pressure_local', 'integrals'),
    }
    for label, rows in constructors.items():
        check('actual_pressure_constructors:' + label, all(rows.values()))
    compact_x = 260181
    check('new_coordinates_after_incoming_prefix', compact_x >= len(graphs[0]) and qb['s'] >= len(graphs[1]))
    check('compact_endpoint_substitutions', _endpoint(nodes, cb, ct, compact_x, 48))
    check('quiet_endpoint_substitutions', _endpoint(nodes, qb, rp, qb['s'], qb['length']))
    check('outer_endpoint_substitutions', _endpoint(nodes, dict(integrals=oc['cumulative_integrals'],
        pressure=oc['pressure'], boundary=oc['boundary']), rc, oc['coordinate'], 48))
    local_metadata = []
    for label, owner, variable, lower in (('compact', cb, compact_x, 1), ('quiet', qb, qb['s'], 0)):
        for key, power in (('S0', 1), ('S1', 2)):
            for order, ref in enumerate(owner['pressure_local'][key]):
                row = nodes[ref]
                valid = row.get('operation') == 'definite_integral' and row['lower'] == lower and row['upper'] == variable
                valid = valid and row['physical_radius_power'] == power and row['physical_Jacobian_applied_once'] is True
                valid = valid and row['nonzero_preceding_S0_S1_retained'] is True and row['ordinary_slow_Z_derivative_order'] == order
                check('signed_local_integral:' + label + ':' + key + ':' + str(order), valid)
                local_metadata.append(dict(chart=label, integral=key, ordinary_Z_order=order, source_node=ref, metadata=row))
    dag, outputs, identities = AffineFunctions(), {}, []
    for moment, integral, boundary, coefficient, sign in (
        ('m13', 'S0', 'Q_plus_Q_y', '-1/2', 1), ('m15', 'S1', 'R_Q_y', '1/2', -1),
    ):
        rows = []
        for order in range(4):
            source = lambda ref: dag.port(ref, ordinary_Z_order=order)
            b0, b1, b2 = [source(endpoint['boundary'][boundary][order]) for endpoint in (rc, ct, rp)]
            bi1, bi2 = [source(endpoint['boundary'][boundary][order]) for endpoint in (cb['inlet'], qb['inlet'])]
            local1 = dag.at(cb['pressure_local'][integral][order], compact_x, 48, order)
            local2 = dag.at(qb['pressure_local'][integral][order], qb['s'], qb['length'], order)
            natural1 = dag.combination([(coefficient, local1), (sign, b1), (-sign, bi1)])
            natural2 = dag.combination([(coefficient, local2), (sign, b2), (-sign, bi2)])
            jump1 = dag.combination([(1, bi1), (-1, b0)])
            jump2 = dag.combination([(1, bi2), (-1, b1)])
            correction1 = dag.combination([(sign, jump1)])
            correction2 = dag.combination([(sign, jump2)])
            increment1 = dag.combination([(coefficient, local1), (sign, b1), (-sign, b0)])
            increment2 = dag.combination([(coefficient, local2), (sign, b2), (-sign, b1)])
            expanded = dag.combination([(1, natural1), (1, natural2), (1, correction1), (1, correction2)])
            telescoped = dag.combination([(coefficient, local1), (coefficient, local2), (sign, b2), (-sign, b0)])
            natural_sum = dag.combination([(1, natural1), (1, natural2)])
            initial_boundary = source(cb['initial'][boundary][order])
            rc_integral = source(rc['integrals'][integral][order])
            rc_pressure = dag.combination([(coefficient, rc_integral), (sign, b0), (-sign, initial_boundary)])
            rp_pressure = dag.combination([(coefficient, rc_integral), (coefficient, local1), (coefficient, local2),
                (sign, b2), (-sign, initial_boundary)])
            absolute = dag.combination([(1, rc_pressure), (1, expanded)])
            proofs = dict(first_increment=dag.equal(increment1, dag.combination([(1, natural1), (1, correction1)])),
                second_increment=dag.equal(increment2, dag.combination([(1, natural2), (1, correction2)])),
                corrected_telescope=dag.equal(expanded, telescoped), absolute_pressure=dag.equal(absolute, rp_pressure),
                unproved_jumps_not_silently_dropped=not dag.equal(natural_sum, telescoped),
                identical_internal_boundary_absent=b1 not in dag.forms[telescoped])
            check(moment + ':Z' + str(order) + ':new_signed_identities', all(proofs.values()))
            identities.append(dict(moment=moment, ordinary_Z_order=order, **proofs))
            rows.append(dict(ordinary_Z_order=order, local_integrals=[local1, local2], natural_segments=[natural1, natural2],
                natural_join_defects=[jump1, jump2], signed_join_corrections=[correction1, correction2],
                predecessor_increments=[increment1, increment2], expanded=expanded, telescoped=telescoped,
                absolute_Rc=rc_pressure, absolute_Rp=rp_pressure, saved_absolute_ports=[rc['pressure'][moment][order], rp['pressure'][moment][order]],
                telescoped_normal_form=dag.normal_form(telescoped), natural_join_zero_proved=False))
        outputs[moment] = rows
    passed = all(checks.values())
    result = dict(schema=SCHEMA, source_only=True, passed=passed, ordinary_Z_orders=[0, 1, 2, 3],
        ordinary_Z_not_temporal_n=True, input_digest=expected_input_digest, source_prerequisite_audit=checks,
        source_family=family, source_witness=dict(archive_graph_counts=[len(g) for g in graphs], archive_graph_digests=list(GRAPH_HASHES),
            constructor_checks=constructors, local_integral_metadata=local_metadata, common_original_domain_required=True,
            genuine_saved_C5_and_pressure_C3_regularities_required=True, numerical_containment_claimed=False),
        function_source_interface=dict(compact=dict(coordinate=compact_x, domain=[1, 2], radius='Rc*x', source=cb),
            quiet=dict(coordinate=qb['s'], domain=[0, qb['length']], radius='R2*exp(s)', source=qb),
            native_Rc=rc, native_Rp=rp, Rp_frame=qb['Rp_frame'], source_graph_digest=GRAPH_HASHES[2],
            quiet_closed_alternative=qb['quiet_closed'], quiet_closed_added_to_generic=False),
        definition_graph=dict(namespace=SCHEMA, nodes=dag.nodes, node_count=len(dag.nodes), nodes_digest=digest(dag.nodes),
            outputs=outputs, proof_graph_not_numeric_evaluator=True), signed_identities=identities,
        retained_construction_data=parent['retained_construction_data'],
        source_inspection_audit=dict(new_source_consumer_calls=1, accepted_ancestor_operator_calls=0,
            old_definition_graphs_rebuilt=False, old_operators_replayed=False, numeric_values_computed=False,
            quadrature_or_norm_calls=0, source_value_reducer_calls=0, temporal_recursion_calls=0),
        open_obligations=['Rc natural physical boundary equality', '2Rc natural physical boundary equality',
            'old212 actual signed-source adapter equality', 'independent complete pulse and reserved terminal match',
            'global physical seams and five moments', 'heat/stress/temporal recursion/corrections/Cartesian NS'],
        unresolved_gates=[name for name, value in checks.items() if not value], runtime_model_metadata=None)
    result.update({name: passed for name in STAGE_FLAGS})
    result.update({name: False for name in FALSE_FLAGS})
    return result
