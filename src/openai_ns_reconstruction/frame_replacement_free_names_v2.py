"""Repair: conservative free-name evidence with explicit edge diagnostics.

The analysis is syntactic and does not prove functional dependence, inlet
identity, capture, or the correctness of the source graph's substitution
semantics. It never expands the protected geometry boundaries.
"""
import base64
import gzip
import hashlib
import json
from collections import Counter

P = 'lei_ren_part1_paper_compliant_current_original_'
ARCHIVE = 'experiments/root_st073/' + P + 'Rp_native_C5_pulse_source.json.gz'
ARCHIVE_SHA = '442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe'
REPORT = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_frame_selected_inlet_scope.json'
REPORT_SHA = '8b1d56bb21ba81c715a6f430e08ab9a2e66edbeebe21c68c5eec439f2a09d92d'
RESULT = 'tmp/frame_inlet_lexical_audit_v1/result.json'
RESULT_SHA = '0dc68e71d1727458b8d3a274863eeec9a077431335a05b8c0a05240d1a1cbb6d'
PROTECTED = frozenset((17, 33572, 25220, 263344, 263348, 263349, 262918, 205147, 263346))
TARGET_NAMES = ('Z', 'original_Gamma_expectation_v', 'original_Gamma_log_radius_offset')
SCHEMA = {'exact_rational': (), 'mathematical_pi': (), 'mathematical_positive_infinity': (),
          'sum': ('arguments',), 'product': ('arguments',), 'negative': ('argument',),
          'analytic_unary': ('argument',), 'positive_quotient': ('numerator', 'denominator'),
          'exact_real_comparison': ('left', 'right'), 'logical_or': ('arguments',),
          'compact_raw_beta': ('argument',), 'exact_original_source_flat_sigma': ('argument',),
          'exact_original_source_heat_phi': ('argument',)}
MAX_STATES = 100000


def inspect(nodes, root, names, protected=PROTECTED):
    """Overapproximate free names using a conservative lexical model.

    A substitution's value is traversed conditionally: the target variable
    might not occur in its expression. Unknown nodes prevent exclusion claims.
    """
    stack = [(root, (), False, ())]
    seen, unknown, ops = set(), [], Counter()
    direct = {n: [] for n in names}
    conditional = {n: [] for n in names}
    while stack:
        node, scope, maybe, path = stack.pop()
        key = (node, scope, maybe)
        if key in seen:
            continue
        if len(seen) >= MAX_STATES:
            unknown.append({'reason': 'state_limit', 'limit': MAX_STATES})
            break
        seen.add(key)
        if type(node) is not int or not 0 <= node < len(nodes):
            unknown.append({'reason': 'invalid_edge', 'node': repr(node)})
            continue
        if node in protected:
            unknown.append({'reason': 'protected_geometry_boundary', 'node': node})
            continue
        row = nodes[node]
        if not isinstance(row, dict):
            unknown.append({'reason': 'invalid_row', 'node': node})
            continue
        op = row.get('operation')
        ops[str(op)] += 1
        if op == 'bound_variable':
            name = row.get('name')
            if not isinstance(name, str) or not name:
                unknown.append({'reason': 'invalid_variable_name', 'node': node})
            elif name in direct and name not in scope:
                bucket = conditional if maybe else direct
                if len(bucket[name]) < 3:
                    bucket[name].append({'node': node, 'path': list(path), 'scope': list(scope)})
            continue
        def push(child, active, is_conditional, edge):
            if type(child) is not int or not 0 <= child < len(nodes):
                unknown.append({'reason': 'invalid_edge', 'node': node, 'edge': edge, 'child': repr(child)})
            elif child == node or any(parent == child for parent, _ in path):
                unknown.append({'reason': 'cyclic_edge', 'node': node, 'edge': edge, 'child': child})
            else:
                stack.append((child, active, is_conditional, path + ((node, edge),)))
        if op == 'definite_integral':
            var = row.get('variable')
            if not isinstance(var, str) or not var or any(k not in row for k in ('lower', 'upper', 'integrand')):
                unknown.append({'reason': 'malformed_integral', 'node': node})
                continue
            push(row['lower'], scope, maybe, 'lower')
            push(row['upper'], scope, maybe, 'upper')
            push(row['integrand'], scope + (var,), maybe, 'integrand')
        elif op == 'function_substitution':
            var = row.get('variable')
            if type(var) is not int or not 0 <= var < len(nodes) or var in protected or var == node:
                unknown.append({'reason': 'invalid_or_protected_substitution_variable', 'node': node,
                                'variable': repr(var), 'value_type': type(var).__name__,
                                'protected': type(var) is int and var in protected})
                continue
            decl = nodes[var]
            if not isinstance(decl, dict) or decl.get('operation') != 'bound_variable' or not isinstance(decl.get('name'), str) or not decl['name'] or any(k not in row for k in ('expression', 'value')):
                unknown.append({'reason': 'malformed_substitution', 'node': node})
                continue
            push(row['expression'], scope + (decl['name'],), maybe, 'expression')
            push(row['value'], scope, True, 'value_conditional')
        elif op in SCHEMA:
            for field in SCHEMA[op]:
                value = row.get(field)
                if value is None:
                    unknown.append({'reason': 'missing_operand', 'node': node, 'field': field})
                elif isinstance(value, list):
                    for index, child in enumerate(value):
                        push(child, scope, maybe, field + '[' + str(index) + ']')
                else:
                    push(value, scope, maybe, field)
        else:
            unknown.append({'reason': 'unsupported_operation', 'node': node, 'operation': op})
    results = {}
    for name in names:
        if direct[name]:
            state = 'unbound_syntactic_witness'
        elif conditional[name]:
            state = 'conditional_syntactic_witness'
        elif unknown:
            state = 'unknown_due_to_unexpanded_nodes'
        else:
            state = 'excluded_in_lexical_model'
        results[name] = {'state': state, 'direct_witnesses': direct[name],
                         'conditional_witnesses': conditional[name]}
    return {'root': root, 'states_visited': len(seen), 'operations': dict(ops),
            'names': results, 'unknown': unknown}


def recover(root):
    report_raw = (root / REPORT).read_bytes()
    if hashlib.sha256(report_raw).hexdigest() != REPORT_SHA:
        raise ValueError('saved inlet scope report hash mismatch')
    report = json.loads(report_raw)
    result_raw = base64.b64decode(report['artifacts'][RESULT]['base64'], validate=True)
    if hashlib.sha256(result_raw).hexdigest() != RESULT_SHA:
        raise ValueError('saved inlet scope result hash mismatch')
    saved = json.loads(result_raw)
    if saved.get('status') != 'returned' or saved['result'].get('status') != 'outer_contexts_unresolved':
        raise ValueError('unexpected saved scope result')
    if saved['result'].get('archive_sha256') != ARCHIVE_SHA:
        raise ValueError('saved scope archive pin mismatch')
    hits = [hit for row in saved['result']['records'] for hit in row['hits'] if hit['binders']]
    pairs = {(hit['old'], hit['replacement']) for hit in hits}
    names = {binder[2] for hit in hits for binder in hit['binders']}
    if pairs != {(31655, 265772)} or names != set(TARGET_NAMES):
        raise ValueError('saved bound-hit target or binder names changed')
    raw = (root / ARCHIVE).read_bytes()
    if hashlib.sha256(raw).hexdigest() != ARCHIVE_SHA:
        raise ValueError('raw gzip hash mismatch before decode')
    archive = json.loads(gzip.decompress(raw))
    frame = archive['Rp_native_C5_pulse_source_bindings']
    if frame['input_substitution'].get('31655') != 265772:
        raise ValueError('saved replacement differs from bound-hit receipt')
    nodes = archive['exact_graph_nodes']
    # Lexical fixtures are new and separate from all archived source graphs.
    fixture = [{'operation': 'bound_variable', 'name': 'x'},
               {'operation': 'bound_variable', 'name': 'y'},
               {'operation': 'sum', 'arguments': [0, 1]},
               {'operation': 'definite_integral', 'variable': 'x', 'integrand': 2, 'lower': 1, 'upper': 1},
               {'operation': 'function_substitution', 'variable': 0, 'expression': 2, 'value': 1},
               {'operation': 'unknown_fixture'},
               {'operation': 'function_substitution', 'variable': 0, 'expression': 0, 'value': 1},
               {'operation': 'definite_integral', 'variable': 'x', 'integrand': 6, 'lower': 1, 'upper': 1}]
    a, b, c = (inspect(fixture, i, ('x', 'y'), frozenset()) for i in (3, 4, 5))
    d = inspect(fixture, 6, ('x', 'y'), frozenset())
    e = inspect(fixture, 3, ('x',), frozenset((2,)))
    f = inspect(fixture, 7, ('x', 'y'), frozenset())
    if not (a['names']['x']['state'] == 'excluded_in_lexical_model' and
            a['names']['y']['state'] == 'unbound_syntactic_witness' and
            b['names']['x']['state'] == 'excluded_in_lexical_model' and
            b['names']['y']['state'] == 'unbound_syntactic_witness' and
            c['names']['x']['state'] == 'unknown_due_to_unexpanded_nodes' and
            d['names']['y']['state'] == 'conditional_syntactic_witness' and
            e['names']['x']['state'] == 'unknown_due_to_unexpanded_nodes' and
            f['names']['x']['state'] == 'excluded_in_lexical_model'):
        raise ValueError('new free-name fixtures failed')
    analysis = inspect(nodes, 265772, TARGET_NAMES)
    contexts = Counter(tuple(tuple(b) for b in hit['binders']) for hit in hits)
    overlap = []
    for context, count in sorted(contexts.items()):
        relevant = sorted({binder[2] for binder in context})
        overlap.append({'binder_stack': context, 'saved_state_count_across_roots': count,
                        'name_states': {name: analysis['names'][name]['state'] for name in relevant}})
    return {'schema': 'frame_replacement_free_names_v2',
            'status': 'unresolved_source_semantics' if analysis['unknown'] else 'lexical_model_inspected',
            'saved_hit_count_including_cross_root_repetition': len(hits), 'saved_pair': [31655, 265772],
            'source_archive_sha256': ARCHIVE_SHA, 'saved_result_sha256': RESULT_SHA,
            'fixture_passed': True, 'replacement_analysis': analysis, 'saved_binder_contexts': overlap,
            'syntactic_evidence_only': True, 'functional_dependence_proved': False,
            'actual_capture_proved': False, 'actual_capture_free_proved': False,
            'ordinary_Z_transport_proved': False, 'physical_endpoint_proved': False,
            'source_values_computed': False, 'old_errors_removed': False, 'runtime_model_metadata': None}
