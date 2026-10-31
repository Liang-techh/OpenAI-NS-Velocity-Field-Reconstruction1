"""Inspect the deferred first-micro integral integrands without evaluating them.

The integral's string binder is not promoted to a graph identity.  This
consumer reports only candidate bound-variable handles and opaque frontiers.
"""
SCHEMA = 'original_q4_integrand_binding_v1'
MAX_VISITS = 3000
MAX_DEPTH = 100
MAX_FRONTIER = 512
CHILDREN = {
    'sum': ('arguments',), 'product': ('arguments',),
    'negative': ('argument',), 'positive_quotient': ('numerator', 'denominator'),
    'analytic_unary': ('argument',), 'exact_real_comparison': ('left', 'right'),
    'logical_or': ('arguments',),
}


class Invalid(ValueError):
    pass


def get_node(nodes, ref):
    if type(ref) is not int or ref < 0 or ref >= len(nodes) or not isinstance(nodes[ref], dict):
        raise Invalid('invalid source reference')
    return nodes[ref]


def inspect(graph, frontier):
    """Return a bounded candidate-variable inspection, never an identity proof."""
    try:
        if not isinstance(graph, dict) or not isinstance(graph.get('exact_graph_nodes'), list):
            raise Invalid('missing source graph')
        if not isinstance(frontier, dict) or frontier.get('operation') != 'definite_integral':
            raise Invalid('not an integral frontier')
        nodes = graph['exact_graph_nodes']
        integral_ref = frontier.get('node_id')
        row = get_node(nodes, integral_ref)
        name = frontier.get('binder_name')
        if row.get('operation') != 'definite_integral' or not isinstance(name, str) or not name or row.get('variable') != name:
            raise Invalid('integral binder mismatch')
        body = frontier.get('integrand_node_id')
        if row.get('integrand') != body:
            raise Invalid('integrand reference mismatch')
        get_node(nodes, body)
        lo, hi = row.get('lower'), row.get('upper')
        if get_node(nodes, lo) != {'operation': 'exact_rational', 'numerator': 0, 'denominator': 1}:
            raise Invalid('lower bound is not exact zero')
        if get_node(nodes, hi) != {'operation': 'exact_rational', 'numerator': 1, 'denominator': 1}:
            raise Invalid('upper bound is not exact one')
        visited, active, seen = set(), set(), set()
        unshadowed, shadowed = set(), set()
        pending = []
        truncated = False

        def stop(ref, op, reason):
            nonlocal truncated
            if len(pending) >= MAX_FRONTIER:
                truncated = True
                return
            pending.append({'node': ref, 'operation': op, 'reason': reason})

        def walk(ref, under_shadow, depth):
            nonlocal truncated
            row = get_node(nodes, ref)
            op = row.get('operation')
            key = (ref, under_shadow)
            if key in active:
                stop(ref, op, 'source cycle')
                return
            if key in seen:
                return
            if len(visited) >= MAX_VISITS or depth > MAX_DEPTH:
                truncated = True
                stop(ref, op, 'inspection budget')
                return
            seen.add(key); visited.add(ref); active.add(key)
            try:
                if op == 'bound_variable':
                    varname = row.get('name')
                    if not isinstance(varname, str) or not varname:
                        raise Invalid('malformed bound variable')
                    if varname == name:
                        (shadowed if under_shadow else unshadowed).add(ref)
                    return
                if op == 'exact_rational':
                    if type(row.get('numerator')) is not int or type(row.get('denominator')) is not int or row['denominator'] <= 0:
                        raise Invalid('malformed exact rational')
                    return
                if op == 'function_substitution':
                    binder = get_node(nodes, row.get('variable'))
                    nested_name = binder.get('name')
                    if binder.get('operation') != 'bound_variable' or not isinstance(nested_name, str) or not nested_name:
                        stop(ref, op, 'nested binder is not a typed variable')
                        return
                    # Candidate lexical interpretation only; name identity is
                    # not certified, so this does not establish capture safety.
                    walk(row.get('value'), under_shadow, depth + 1)
                    walk(row.get('expression'), under_shadow or nested_name == name, depth + 1)
                    return
                if op not in CHILDREN:
                    # Includes nested integrals, lazy branches, core-source
                    # atoms and unknown operators. Do not inspect children.
                    stop(ref, op, 'source operation opaque')
                    return
                for field in CHILDREN[op]:
                    children = row.get(field)
                    if field != 'arguments':
                        children = [children]
                    if not isinstance(children, list) or not children:
                        raise Invalid('malformed source children')
                    for child in children:
                        walk(child, under_shadow, depth + 1)
            finally:
                active.remove(key)

        walk(body, False, 0)
        all_candidates = unshadowed | shadowed
        conflict = len(all_candidates) > 1
        complete = not pending and not truncated and not conflict
        return {
            'schema': SCHEMA, 'status': 'observed',
            'evidence': {'integral_node': integral_ref, 'integrand_node': body,
                         'binder_name': name, 'unshadowed_binder_nodes': sorted(unshadowed),
                         'shadowed_binder_nodes': sorted(shadowed), 'same_name_node_conflict': conflict,
                         'unit_interval_bounds_verified': True, 'visited_nodes': len(visited),
                         'inspection_truncated': truncated, 'integral_source_row': row},
            'frontier': pending,
            'scope': {'integral_binder_identity_proved': False, 'integrand_scan_complete': complete,
                      'candidate_name_scope_only': True, 'numeric_field_value_computed': False,
                      'old_errors_removed': False, 'source_values_evaluated': False},
        }
    except (Invalid, KeyError, TypeError) as exc:
        return {'schema': SCHEMA, 'status': 'mismatch', 'error': str(exc),
                'scope': {'integral_binder_identity_proved': False, 'integrand_scan_complete': False,
                          'numeric_field_value_computed': False, 'old_errors_removed': False}}


def inspect_core(graph, frontier):
    """Expose the actual analytic-core source interface, not an enclosure value."""
    try:
        if not isinstance(graph, dict) or not isinstance(graph.get('exact_graph_nodes'), list):
            raise Invalid('missing source graph')
        if not isinstance(frontier, dict) or frontier.get('operation') != 'current_original_unique_analytic_core_component':
            raise Invalid('not an analytic core frontier')
        nodes = graph['exact_graph_nodes']
        ref = frontier.get('node_id')
        row = get_node(nodes, ref)
        if row.get('operation') != frontier['operation']:
            raise Invalid('core source operation mismatch')
        quantity = row.get('quantity')
        if quantity not in ('Phi', 'Psi', 'mean_Psi'):
            raise Invalid('unknown core quantity')
        for key in ('radial_order', 'ordinary_slow_Z_derivative_order'):
            if type(row.get(key)) is not int or row[key] < 0:
                raise Invalid('invalid ' + key)
        for key in ('actual_twenty_term_operator', 'unique_analytic_Banach_solution',
                    'admitted_solution_not_a_model_or_enclosure_center'):
            if row.get(key) is not True:
                raise Invalid('missing core source declaration: ' + key)
        if not row.get('source_family') or not isinstance(row.get('component_definition'), str) or not row['component_definition']:
            raise Invalid('missing core source definition')
        for key in ('source_provider', 'radial_recurrence_provider', 'actual_source_provider', 'actual_green_inverse_provider'):
            if not row.get(key):
                raise Invalid('missing core provider provenance: ' + key)
        rho, z = row.get('rho'), row.get('Z_variable')
        for key, ref_id in (('rho', rho), ('Z_variable', z)):
            variable = get_node(nodes, ref_id)
            if variable.get('operation') != 'bound_variable' or not isinstance(variable.get('name'), str) or not variable['name']:
                raise Invalid(key + ' is not a typed source variable')
        if rho != frontier.get('expected_binder_node'):
            raise Invalid('core rho differs from saved substitution binder')
        return {
            'schema': SCHEMA, 'status': 'source_ref',
            'evidence': {'component_node': ref, 'quantity': quantity,
                         'radial_order': row['radial_order'],
                         'ordinary_Z_order': row['ordinary_slow_Z_derivative_order'],
                         'rho_node': rho, 'Z_node': z, 'rho_matches_saved_binder': True,
                         'source_row': row},
            'scope': {'source_ref_only': True, 'core_solution_admission_proved': False,
                      'numeric_field_value_computed': False, 'old_errors_removed': False},
        }
    except (Invalid, KeyError, TypeError) as exc:
        return {'schema': SCHEMA, 'status': 'mismatch', 'error': str(exc),
                'scope': {'source_ref_only': False, 'core_solution_admission_proved': False,
                          'numeric_field_value_computed': False, 'old_errors_removed': False}}
