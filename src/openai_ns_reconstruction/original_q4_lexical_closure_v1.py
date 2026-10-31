"""Create deferred, explicitly scoped candidates for saved Q4 substitutions.

No source expression or replacement is evaluated.  Node-ID binding is recorded
as a candidate; source occurrence-to-binder identity remains unproved.
"""
from copy import deepcopy
import hashlib
import json

SCHEMA = 'original_q4_lexical_closure_v1'


class Invalid(ValueError):
    pass


def node(nodes, ref):
    if type(ref) is not int or ref < 0 or ref >= len(nodes) or not isinstance(nodes[ref], dict):
        raise Invalid('invalid source node')
    return nodes[ref]


def node_hash(row):
    raw = json.dumps(row, sort_keys=True, separators=(',', ':'), ensure_ascii=True).encode()
    return hashlib.sha256(raw).hexdigest()


def valid_dyadic(bits):
    if not isinstance(bits, list) or len(bits) != 4 or any(type(x) is not int for x in bits):
        raise Invalid('invalid dyadic')
    sign, m, exponent, count = bits
    if sign not in (0, 1) or m < 0 or count != m.bit_length() or (m and m % 2 == 0):
        raise Invalid('noncanonical dyadic')
    if m == 0 and bits != [0, 0, 0, 0]:
        raise Invalid('noncanonical zero')
    return True


def check_env(nodes, env):
    if not isinstance(env, list):
        raise Invalid('invalid lexical environment')
    for entry in env:
        if not isinstance(entry, dict) or not isinstance(entry.get('binding'), dict):
            raise Invalid('invalid environment binding')
        source = node(nodes, entry.get('variable_node'))
        if source.get('operation') not in ('bound_variable', 'current_original_source_parameter'):
            raise Invalid('environment entry is not a typed source variable')
        binding = entry['binding']
        if binding.get('kind') != 'exact_dyadic':
            raise Invalid('outer binding must be an explicit exact dyadic')
        valid_dyadic(binding.get('tuple'))


def scope():
    return {'numeric_field_value_computed': False, 'source_value_evaluated': False,
            'old_errors_removed': False, 'occurrence_binding_identity_verified': False,
            'capture_avoidance_proved': False, 'replacement_value_evaluated': False,
            'corrected_owner_identity_proved': False}


def close(graph, expression, env):
    """Build a source closure, capturing replacement in the outer environment.

    The node-ID environment is only a candidate representation until source
    occurrence/binder identity is independently established.
    """
    try:
        if not isinstance(graph, dict) or not isinstance(graph.get('exact_graph_nodes'), list):
            raise Invalid('missing source graph')
        nodes = graph['exact_graph_nodes']
        if not isinstance(expression, dict) or expression.get('kind') != 'expanded_source_node' or expression.get('operation') != 'function_substitution':
            raise Invalid('not an expanded substitution')
        ref = expression.get('node')
        row = node(nodes, ref)
        if row != expression.get('source_row') or row.get('operation') != 'function_substitution':
            raise Invalid('saved substitution row differs from source')
        children = expression.get('source_children')
        if not isinstance(children, dict):
            raise Invalid('missing source children')
        for key in ('expression', 'variable', 'value'):
            child = children.get(key)
            if not isinstance(child, dict) or child.get('kind') != 'source_graph_ref' or child.get('node') != row.get(key):
                raise Invalid('source child mismatch: ' + key)
            index = row.get(key)
            if type(index) is not int or index < 0 or index >= len(nodes):
                raise Invalid('source child out of range: ' + key)
        binder = node(nodes, row['variable'])
        if binder.get('operation') != 'bound_variable' or not isinstance(binder.get('name'), str) or not binder['name']:
            raise Invalid('substitution binder is not a typed bound variable')
        check_env(nodes, env)
        outer = deepcopy(env)
        replacement = {'kind': 'source_thunk', 'source_node': row['value'],
                       'environment': deepcopy(outer)}
        extended = outer + [{'variable_node': row['variable'], 'binding': replacement}]
        result = {'kind': 'source_closure', 'expression_node': row['expression'],
                  'environment': extended, 'source_substitution_node': ref,
                  'binder_node': row['variable'], 'binder_name': binder['name'],
                  'node_identity_is_unverified_candidate': True}
        return {'schema': SCHEMA, 'status': 'closure_created', 'expression': result, 'scope': scope()}
    except (Invalid, KeyError, TypeError) as exc:
        return {'schema': SCHEMA, 'status': 'mismatch', 'error': str(exc), 'scope': scope()}


def bind_saved(graph, expanded_row, port, sc_digest, coordinate_digest):
    """Attach candidate closures to one saved endpoint expression.

    The caller independently authenticates the graph, prior expansion, port,
    and expected record hashes. Only already-expanded substitution rows and
    their typed binders are read; child expressions are never evaluated.
    """
    try:
        nodes = graph['exact_graph_nodes']
        if not isinstance(nodes, list) or expanded_row.get('status') != 'expanded':
            raise Invalid('missing expanded source')
        parameter, coordinate = port['parameter_node'], port['coordinate_node']
        if parameter == coordinate:
            raise Invalid('parameter and coordinate collide')
        scrow, crow = node(nodes, parameter), node(nodes, coordinate)
        if node_hash(scrow) != sc_digest or node_hash(crow) != coordinate_digest:
            raise Invalid('parameter or coordinate source record hash mismatch')
        if scrow.get('operation') != 'current_original_source_parameter' or scrow.get('name') != 'sc' or scrow.get('Z_independent') is not True:
            raise Invalid('not the Z-independent sc source parameter')
        if crow.get('operation') != 'bound_variable':
            raise Invalid('coordinate is not a typed bound variable')
        sc = port['parameter_exact_mpf_tuple']
        half = port['coordinate_exact_mpf_tuple']
        valid_dyadic(sc); valid_dyadic(half)
        if sc[0] != 0 or sc[1] <= 0 or half != [sc[0], sc[1], sc[2] - 1, sc[3]]:
            raise Invalid('invalid positive endpoint instance')
        env = [{'variable_node': parameter, 'binding': {'kind': 'exact_dyadic', 'tuple': deepcopy(sc)}},
               {'variable_node': coordinate, 'binding': {'kind': 'exact_dyadic', 'tuple': deepcopy(half)}}]
        check_env(nodes, env)
        saved = expanded_row['expression']
        if not isinstance(saved, dict):
            raise Invalid('missing saved expression')
        decisions, memo = [], {}

        def visit(value):
            if isinstance(value, list):
                return [visit(x) for x in value]
            if not isinstance(value, dict):
                return value
            if value.get('kind') == 'expanded_source_node' and value.get('operation') == 'function_substitution':
                ref = value.get('node')
                if ref not in memo:
                    memo[ref] = close(graph, value, env)
                decision = memo[ref]
                decisions.append({'node': ref, 'status': decision['status']})
                if decision['status'] != 'closure_created':
                    return deepcopy(value)
                # Preserve the entire source row as provenance; no child is traversed.
                return {'kind': 'candidate_lexical_application', 'source': deepcopy(value),
                        'closure': deepcopy(decision['expression']),
                        'occurrence_binding_identity_verified': False}
            return {key: visit(item) for key, item in value.items()}

        expression = visit(saved)
        created = sum(result['status'] == 'closure_created' for result in memo.values())
        return {'schema': SCHEMA, 'status': 'candidate_closures' if created else 'unresolved',
                'expression': expression, 'port': deepcopy(port),
                'unique_closures': [{'node': ref, 'status': result['status'],
                                     'error': result.get('error')} for ref, result in sorted(memo.items())],
                'unique_created': created, 'occurrences': len(decisions),
                'flat_branches': deepcopy(expanded_row['flat_branches']),
                'scope': {**scope(), 'source_archive_and_port_authentication_required_of_caller': True,
                          'binder_and_coordinate_record_types_checked': True,
                          'source_expression_children_evaluated': False}}
    except (Invalid, KeyError, TypeError) as exc:
        return {'schema': SCHEMA, 'status': 'unresolved', 'error': str(exc),
                'port': deepcopy(port), 'scope': scope()}
