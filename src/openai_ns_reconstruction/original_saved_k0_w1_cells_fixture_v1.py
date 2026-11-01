"""Fresh structural fixture for saved K0 profile binding; no source inputs."""
from copy import deepcopy


def _interval():
    return {'lower_exact_mpf_tuple': [0, 1, -2, 1], 'upper_exact_mpf_tuple': [0, 1, -1, 1], 'directed': True}


def _synthetic():
    cells = []
    for i in range(8):
        rows = {}
        for order in ('Z0', 'Z1'):
            rows[order] = {}
            for name in ('F1', 'Uz1'):
                rows[order][name] = [{'coefficient': _interval(), 'F0_power': 2, 'Lambda_power': -3,
                                      'source_cell': i, 'target_cell': i,
                                      'matrix_sector_identity': 'synthetic-matrix-' + str(i),
                                      'W0_sector_identity': 'synthetic-w0-' + str(i),
                                      'other_provenance': {'retained': True}}]
        cells.append({'index': i, 'source_indices': list(range(i + 1)), 'triangular_diagonal_source': i,
                      'axis_endpoint_continuation': i == 0, 't': _interval(), 'x': _interval(), 'rho': _interval(),
                      'W1_by_order': rows})
    return {'source_family': {'actual_five_defect_family_sha256': 'a'*64, 'implicit_source_sha256': 'b'*64,
                              'datum_enclosure_sha256': 'c'*64},
            'source_logs': {'synthetic': 'only'}, 'source_kernel_scope': 'synthetic K0 only',
            'Kchi_minus_K0_included': False, 'allword_tail_added': False, 'g4_pressure_added': False,
            'appended_to_cumulative_m15': False, 'complete_m15_computed': False, 'terminal_closure_proved': False,
            'profile': {'target_cells': cells}}


def _rejects(leaf, obj):
    try:
        leaf.bind(obj)
    except (TypeError, ValueError, KeyError):
        return True
    return False


def run_fixture(leaf):
    data = _synthetic()
    result = leaf.bind(data)
    term = result['cells'][3]['ordinary_Z']['Z1']['F1'][0]
    checks = {
        'eight_cells_preserved': len(result['cells']) == 8,
        'powers_preserved': term['formal_source_powers'] == [0, 2, 0, -3],
        'binary_endpoints_unchanged': term['coefficient'] == _interval(),
        'provenance_unchanged': term['original_term']['other_provenance']['retained'] is True,
        'Kchi_and_tails_not_claimed': result['Kchi_minus_K0_included'] is False and result['allword_tail_added'] is False,
        'local_F0_not_identified_with_anchor': result['local_F0_is_a_distinct_unbound_scale'] is True,
        'no_field_or_ledger_claim': result['actual_W1_complete'] is False and result['appended_to_cumulative_m15'] is False,
    }
    bad = deepcopy(data); bad['profile']['target_cells'][0]['W1_by_order']['Z1'].pop('Uz1')
    checks['missing_derivative_rejected'] = _rejects(leaf, bad)
    bad = deepcopy(data); bad['profile']['target_cells'][1]['W1_by_order']['Z0']['F1'][0]['target_cell'] = 2
    checks['wrong_target_rejected'] = _rejects(leaf, bad)
    bad = deepcopy(data); bad['Kchi_minus_K0_included'] = True
    checks['unrelated_kernel_scope_rejected'] = _rejects(leaf, bad)
    bad = deepcopy(data); del bad['profile']['target_cells'][0]['W1_by_order']['Z0']['F1'][0]['Lambda_power']
    checks['missing_power_rejected'] = _rejects(leaf, bad)
    bad = deepcopy(data); bad['profile']['target_cells'][0]['W1_by_order']['Z0']['F1'][0]['coefficient'] = 1
    checks['scalar_coefficient_rejected'] = _rejects(leaf, bad)
    return {'schema': 'original_saved_k0_w1_cells_fixture_v1', 'passed': all(checks.values()),
            'checks': checks, 'synthetic_only': True, 'source_values_computed': False}
