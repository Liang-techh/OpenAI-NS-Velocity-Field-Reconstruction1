"""Expand selected previously opaque Q4 nodes one layer without re-evaluation.

This is an exact source-structure transform, not an evaluator. Child graph
references stay deferred, including both sides and the binder of substitution.
The caller must authenticate the source graph and previous result separately.
"""
from copy import deepcopy

SCHEMA = 'original_q4_endpoint_expand_v1'
TARGETS = {
    'negative': ('argument',),
    'analytic_unary': ('argument',),
    'function_substitution': ('expression', 'variable', 'value'),
}


class Mismatch(ValueError):
    pass


def expand(graph, prior_result):
    """Replace only old opaque refs; do not recurse into newly revealed rows."""
    expression = prior_result.get('expression') if isinstance(prior_result, dict) else None
    flat = prior_result.get('flat_branches') if isinstance(prior_result, dict) else None
    counts = {name: 0 for name in TARGETS}
    unique = {name: set() for name in TARGETS}
    remaining = set()
    try:
        if not isinstance(graph, dict) or not isinstance(graph.get('exact_graph_nodes'), list):
            raise Mismatch('missing source graph')
        if not isinstance(prior_result, dict) or prior_result.get('status') != 'unresolved':
            raise Mismatch('expected prior unresolved result')
        if not isinstance(expression, dict) or not isinstance(flat, list):
            raise Mismatch('missing prior expression or flat branches')
        nodes = graph['exact_graph_nodes']
        inactive = {row['inactive_body_not_visited'] for row in flat
                    if isinstance(row, dict) and type(row.get('inactive_body_not_visited')) is int}
        memo = {}

        def visit(value):
            if isinstance(value, list):
                return [visit(item) for item in value]
            if not isinstance(value, dict):
                return value
            if value.get('kind') != 'opaque_source_ref':
                return {key: visit(item) for key, item in value.items()}
            ref, operation = value.get('node'), value.get('operation')
            if type(ref) is not int or ref < 0 or ref >= len(nodes):
                raise Mismatch('opaque source reference out of range')
            if ref in inactive:
                raise Mismatch('prior inactive body cannot be expanded')
            if operation not in TARGETS:
                remaining.add((ref, str(operation)))
                return deepcopy(value)
            counts[operation] += 1
            unique[operation].add(ref)
            # The previous result may repeat a graph handle in several source
            # expressions; inspect that newly opened row at most once.
            if ref in memo:
                template = memo[ref]
                if template['operation'] != operation:
                    raise Mismatch('same node has inconsistent operations')
                result = deepcopy(template)
                result['source_ref'] = deepcopy(value)
                return result
            row = nodes[ref]
            if not isinstance(row, dict) or row.get('operation') != operation:
                raise Mismatch('opaque source operation mismatch')
            children = {}
            for key in TARGETS[operation]:
                child = row.get(key)
                if type(child) is not int or child < 0 or child >= len(nodes):
                    raise Mismatch('missing or invalid source child ' + key)
                # No child row is inspected or evaluated here.
                children[key] = {'kind': 'source_graph_ref', 'node': child}
            if operation == 'analytic_unary' and (not isinstance(row.get('name'), str) or not row['name']):
                raise Mismatch('missing analytic unary name')
            result = {'kind': 'expanded_source_node', 'node': ref, 'operation': operation,
                      'source_ref': deepcopy(value), 'source_row': deepcopy(row),
                      'source_children': children, 'substitution_applied': False}
            for key in ('source_reference_id', 'assumptions'):
                if key in row:
                    result[key] = deepcopy(row[key])
            if operation == 'function_substitution' and 'Z_independent_substitution' in row:
                result['Z_independent_substitution'] = row['Z_independent_substitution']
            memo[ref] = deepcopy(result)
            return result

        output = visit(expression)
        status = 'expanded' if any(unique.values()) else 'unresolved'
        return {
            'schema': SCHEMA, 'status': status, 'expression': output,
            'flat_branches': deepcopy(flat), 'occurrences_by_operation': counts,
            'opened_nodes_by_operation': {name: sorted(refs) for name, refs in unique.items()},
            'remaining_opaque_refs': [{'node': ref, 'operation': op} for ref, op in sorted(remaining)],
            'scope': {'source_structure_only': True, 'new_source_child_rows_evaluated': False,
                      'substitution_applied': False, 'capture_avoiding_substitution_completed': False,
                      'numeric_field_value_computed': False, 'corrected_owner_identity_proved': False,
                      'old_errors_removed': False, 'flat_branch_decisions_recomputed': False},
        }
    except Mismatch as exc:
        return {'schema': SCHEMA, 'status': 'mismatch', 'error': str(exc),
                'expression': deepcopy(expression), 'flat_branches': deepcopy(flat),
                'scope': {'numeric_field_value_computed': False, 'old_errors_removed': False}}
