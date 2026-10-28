"""Bounded, conservative occurrence audit of selected pulse inlet rewrites.

Only saved ap/C1/C2 roots through ordinary-Z order two are inspected.  No
replacement subtree is traversed: occurrences beneath a binder are unknown,
not certified.  This does not check a physical endpoint or Z transport.
"""
import gzip
import hashlib
import json
from collections import Counter

P = 'lei_ren_part1_paper_compliant_current_original_'
ARCHIVE = 'experiments/root_st073/' + P + 'Rp_native_C5_pulse_source.json.gz'
ARCHIVE_SHA = '442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe'
C3 = P + 'n1_pulse_pressure_C3_functions.json.gz'
C3_SHA = 'bdfcc70f57d8cea981f508074690a7e940f884382c2452805f126a3f96bdb59d'
MAX_STATES_PER_ROOT = 100000
FORBIDDEN = frozenset((17, 33572, 25220, 263344, 263348, 263349, 262918, 205147, 263346))
CHILDREN = {
    'exact_rational': (), 'bound_variable': (), 'mathematical_pi': (),
    'mathematical_positive_infinity': (), 'sum': ('arguments',),
    'product': ('arguments',), 'negative': ('argument',),
    'analytic_unary': ('argument',), 'positive_quotient': ('numerator', 'denominator'),
    'exact_real_comparison': ('left', 'right'), 'logical_or': ('arguments',),
    'compact_raw_beta': ('argument',), 'exact_original_source_flat_sigma': ('argument',),
    'exact_original_source_heat_phi': ('argument',),
}


def inspect(nodes, mapping, roots, forbidden=FORBIDDEN):
    """Record inlet occurrences in lexical context, stopping at unknowns.

    No free-variable claim is made for replacements.  An occurrence in an
    empty binder context cannot be captured by an enclosing binder.  A
    nonempty context, forbidden boundary, malformed edge or unknown operation
    is unresolved.  Memoization includes the complete binder context.
    """
    records = []
    for label, root in roots:
        stack = [(root, ())]
        seen = set()
        hits, unresolved, operations = [], [], Counter()
        while stack:
            node, scope = stack.pop()
            key = (node, scope)
            if key in seen:
                continue
            if len(seen) >= MAX_STATES_PER_ROOT:
                unresolved.append({'reason': 'state_limit', 'limit': MAX_STATES_PER_ROOT})
                break
            seen.add(key)
            if type(node) is not int or node < 0 or node >= len(nodes):
                unresolved.append({'reason': 'invalid_edge', 'node': repr(node)})
                continue
            if node in forbidden:
                unresolved.append({'reason': 'forbidden_geometry_boundary', 'node': node})
                continue
            if node in mapping:
                hits.append({'old': node, 'replacement': mapping[node],
                             'binders': list(scope), 'outer_capture_excluded': not scope})
                continue
            row = nodes[node]
            if not isinstance(row, dict):
                unresolved.append({'reason': 'invalid_row', 'node': node})
                continue
            op = row.get('operation')
            operations[str(op)] += 1
            def push(child, child_scope):
                if type(child) is not int or child < 0 or child >= node:
                    unresolved.append({'reason': 'non_prior_or_invalid_edge', 'node': node, 'child': repr(child)})
                else:
                    stack.append((child, child_scope))
            if op == 'definite_integral':
                variable = row.get('variable')
                if not isinstance(variable, str) or not variable:
                    unresolved.append({'reason': 'invalid_integral_binder', 'node': node})
                    continue
                if not all(k in row for k in ('integrand', 'lower', 'upper')):
                    unresolved.append({'reason': 'missing_integral_operand', 'node': node})
                    continue
                binder = ('integral', node, variable)
                push(row['lower'], scope)
                push(row['upper'], scope)
                push(row['integrand'], scope + (binder,))
            elif op == 'function_substitution':
                variable = row.get('variable')
                if type(variable) is not int or variable < 0 or variable >= node:
                    unresolved.append({'reason': 'invalid_substitution_variable', 'node': node})
                    continue
                if variable in mapping:
                    unresolved.append({'reason': 'rewritten_substitution_binder', 'node': node, 'variable': variable})
                    continue
                if variable in forbidden:
                    unresolved.append({'reason': 'forbidden_substitution_binder', 'node': node, 'variable': variable})
                    continue
                var = nodes[variable]
                if not isinstance(var, dict) or var.get('operation') != 'bound_variable' or not isinstance(var.get('name'), str) or not var['name']:
                    unresolved.append({'reason': 'unknown_substitution_binder', 'node': node, 'variable': variable})
                    continue
                if 'expression' not in row or 'value' not in row:
                    unresolved.append({'reason': 'missing_substitution_operand', 'node': node})
                    continue
                binder = ('substitution', node, var['name'])
                push(row['value'], scope)
                push(row['expression'], scope + (binder,))
            elif op in CHILDREN:
                for name in CHILDREN[op]:
                    value = row.get(name)
                    if value is None:
                        unresolved.append({'reason': 'missing_operand', 'node': node, 'operand': name})
                    elif isinstance(value, list):
                        for item in value:
                            push(item, scope)
                    else:
                        push(value, scope)
            else:
                unresolved.append({'reason': 'unsupported_operation', 'node': node, 'operation': op})
        contextual = [x for x in hits if x['binders']]
        records.append({'label': label, 'root': root, 'states_visited': len(seen),
                        'operations': dict(operations), 'hits': hits, 'unresolved': unresolved,
                        'outer_capture_excluded_for_all_hits': bool(hits) and not contextual and not unresolved,
                        'active_binder_hits': len(contextual)})
    return records


def recover(root):
    raw = (root / ARCHIVE).read_bytes()
    if hashlib.sha256(raw).hexdigest() != ARCHIVE_SHA:
        raise ValueError('raw gzip hash mismatch before decode')
    archive = json.loads(gzip.decompress(raw))
    nodes = archive['exact_graph_nodes']
    frame = archive['Rp_native_C5_pulse_source_bindings']
    if archive['input_hashes'].get(C3) != C3_SHA:
        raise ValueError('saved C3 archive pin mismatch')
    mapping = {}
    for old, new in frame['input_substitution'].items():
        if not isinstance(old, str) or not old.isdecimal() or str(int(old)) != old:
            raise ValueError('noncanonical inlet mapping key')
        index = int(old)
        if index in mapping or index < 2 or index >= len(nodes):
            raise ValueError('invalid inlet placeholder')
        mapping[index] = new
    if not mapping or any(type(v) is not int or v < 0 or v >= len(nodes) for v in mapping.values()):
        raise ValueError('invalid saved input mapping')
    selected = frame['inputs'][C3]['selected_pulse']
    roots = []
    for name in ('ap', 'C1', 'C2'):
        if not isinstance(selected[name], list) or len(selected[name]) < 3:
            raise ValueError('missing selected pulse jet rows')
        roots.extend((name + '[' + str(k) + ']', selected[name][k]) for k in range(3))
    # Small fixtures exercise only the new lexical algorithm, never old sources.
    fixture = [{'operation': 'bound_variable', 'name': 'x'},
               {'operation': 'bound_variable', 'name': 'v'},
               {'operation': 'exact_rational'},
               {'operation': 'definite_integral', 'variable': 'v', 'integrand': 0, 'lower': 2, 'upper': 2},
               {'operation': 'function_substitution', 'variable': 1, 'expression': 0, 'value': 2}]
    fixture_result = inspect(fixture, {0: 1}, [('open', 0), ('integral', 3), ('substitution', 4)], frozenset())
    fixture_ok = (fixture_result[0]['outer_capture_excluded_for_all_hits'] and
                  fixture_result[1]['active_binder_hits'] == 1 and
                  fixture_result[2]['active_binder_hits'] == 1 and
                  not any(r['unresolved'] for r in fixture_result))
    if not fixture_ok:
        raise ValueError('new lexical audit fixture failed')
    rows = inspect(nodes, mapping, roots)
    all_roots_clear = all(r['outer_capture_excluded_for_all_hits'] for r in rows)
    return {'schema': 'frame_selected_inlet_scope_v1', 'archive_sha256': ARCHIVE_SHA,
            'status': 'outer_contexts_clear' if all_roots_clear else 'outer_contexts_unresolved',
            'fixture_ok': fixture_ok, 'all_roots_clear': all_roots_clear,
            'root_scope': 'selected ap/C1/C2 ordinary-Z rows 0,1,2 only',
            'records': rows, 'source_values_computed': False,
            'integrals_evaluated': False, 'replacement_subtrees_traversed': False,
            'actual_capture_free_proved': False, 'ordinary_Z_transport_proved': False,
            'endpoint_identity_proved': False, 'old_errors_removed': False}
