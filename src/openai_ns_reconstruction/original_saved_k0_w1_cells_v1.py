"""Expose saved K0-only W1 target-cell sectors without recomputing them.

The supplied object is an already authenticated, decoded historical output.
No source integral, provider, graph or old helper is invoked.  The old sector
coefficients and their directed binary endpoint records are copied unchanged;
no conversion, midpoint, interval addition or cancellation is performed.
"""
from copy import deepcopy

SCHEMA = 'original_saved_k0_w1_cells_v1'
ORDERS = ('Z0', 'Z1')
FIELDS = ('F1', 'Uz1')
FAMILY = ('actual_five_defect_family_sha256', 'implicit_source_sha256', 'datum_enclosure_sha256')


def _endpoint_record(value, label):
    if not isinstance(value, dict):
        raise ValueError(label + ' must be an encoded interval')
    for key in ('lower_exact_mpf_tuple', 'upper_exact_mpf_tuple'):
        item = value.get(key)
        if not isinstance(item, list) or len(item) != 4 or any(type(x) is not int for x in item):
            raise ValueError(label + ' must retain four-integer binary endpoints')
    return value


def bind(output):
    """Copy the eight saved finite target profiles, preserving their limits.

    Authentication, true source ownership, interval validity, cell coverage,
    missing Kchi terms and all-word tails are separate caller responsibilities.
    The indexed input is the output of the earlier K0 computation, not a new
    full W1 solution and not a candidate for direct ledger accumulation.
    """
    if not isinstance(output, dict):
        raise TypeError('saved output object required')
    family, logs = output.get('source_family'), output.get('source_logs')
    if not isinstance(family, dict) or any(not isinstance(family.get(k), str) or len(family[k]) != 64 for k in FAMILY):
        raise ValueError('saved source-family identifiers required')
    if not isinstance(logs, dict):
        raise ValueError('saved formal source logs required')
    for key in ('Kchi_minus_K0_included', 'allword_tail_added', 'g4_pressure_added',
                'appended_to_cumulative_m15', 'complete_m15_computed', 'terminal_closure_proved'):
        if output.get(key) is not False:
            raise ValueError('saved K0 scope flag must be false: ' + key)
    scope = output.get('source_kernel_scope')
    if not isinstance(scope, str) or 'K0' not in scope:
        raise ValueError('saved K0 kernel scope required')
    profile = output.get('profile')
    cells = profile.get('target_cells') if isinstance(profile, dict) else None
    if not isinstance(cells, list) or len(cells) != 8:
        raise ValueError('eight saved target cells required')
    selected, counts = [], {z: {field: 0 for field in FIELDS} for z in ORDERS}
    for index, cell in enumerate(cells):
        if not isinstance(cell, dict) or cell.get('index') != index:
            raise ValueError('saved cell order or index differs')
        if cell.get('source_indices') != list(range(index + 1)) or cell.get('triangular_diagonal_source') != index:
            raise ValueError('saved source/target cell structure differs')
        result = {'index': index, 'source_indices': deepcopy(cell['source_indices']),
                  'triangular_diagonal_source': index, 'axis_endpoint_continuation': cell.get('axis_endpoint_continuation'),
                  't': deepcopy(_endpoint_record(cell.get('t'), 't')),
                  'x': deepcopy(_endpoint_record(cell.get('x'), 'x')),
                  'rho': deepcopy(_endpoint_record(cell.get('rho'), 'rho')), 'ordinary_Z': {}}
        rows = cell.get('W1_by_order')
        if not isinstance(rows, dict):
            raise ValueError('saved W1 ordinary rows required')
        for order in ORDERS:
            row = rows.get(order)
            if not isinstance(row, dict):
                raise ValueError('saved W1 ordinary Z row missing')
            result['ordinary_Z'][order] = {}
            for field in FIELDS:
                terms = row.get(field)
                if not isinstance(terms, list) or not terms:
                    raise ValueError('saved finite W1 source sectors missing')
                selected_terms = []
                for term in terms:
                    if not isinstance(term, dict):
                        raise ValueError('saved finite sector must be an object')
                    _endpoint_record(term.get('coefficient'), 'coefficient')
                    fpower, lpower = term.get('F0_power'), term.get('Lambda_power')
                    source, target = term.get('source_cell'), term.get('target_cell')
                    if type(fpower) is not int or type(lpower) is not int:
                        raise ValueError('saved integer source powers required')
                    if type(source) is not int or source < 0 or source > index or target != index:
                        raise ValueError('saved sector source/target mismatch')
                    if not isinstance(term.get('matrix_sector_identity'), str) or not isinstance(term.get('W0_sector_identity'), str):
                        raise ValueError('saved matrix and W0 sector identities required')
                    selected_terms.append({'coefficient': deepcopy(term['coefficient']),
                                           'formal_source_powers': [0, fpower, 0, lpower],
                                           'original_term': deepcopy(term)})
                counts[order][field] += len(selected_terms)
                result['ordinary_Z'][order][field] = selected_terms
        selected.append(result)
    return {'schema': SCHEMA, 'status': 'saved_finite_subterm',
            'source_family': deepcopy(family), 'source_logs': deepcopy(logs),
            'formal_source_basis': ['hb', 'F0_anchor', 'Pstar_squared', 'Lambda'],
            'four_power_mapping_scope': 'saved producer mapped each profile sector to [0,F0_power,0,Lambda_power]; this binder does not re-prove that mapping',
            'other_W1_components_omitted': ['K1', 'P1', 'd_x_F1', 'd_x_Uz1'],
            'W0_cells_and_fixed_exit_omitted': True,
            'local_F0_is_a_distinct_unbound_scale': True,
            'ordinary_Z_orders': [0, 1], 'cells': selected, 'term_counts': counts,
            'source_kernel_scope': scope, 'actual_W1_complete': False,
            'Kchi_minus_K0_included': False, 'allword_tail_added': False,
            'pressure_added': False, 'appended_to_cumulative_m15': False,
            'point_value_selected': False, 'new_source_values_computed': False,
            'coverage_recertified': False, 'binder_or_live_owner_verified': False}
