"""Fresh synthetic archive fixture for the sc endpoint consumer."""
from __future__ import annotations

import copy


def example():
    family = {'fixture_family': 'independent_archive_sc_test_A'}
    loop_nodes = [
        {'operation': 'bound_variable', 'name': 'synthetic_native_coordinate_A'},
        {'operation': 'exact_rational', 'numerator': 1, 'denominator': 2},
        {'operation': 'current_original_source_parameter', 'name': 'sc',
         'Z_independent': True, 'source_family': family},
        {'operation': 'product', 'arguments': [2, 1]},
        {'operation': 'radius_fixture', 'tag': 'A'},
        {'operation': 'symbol', 'name': 'synthetic_upper_A'},
    ]
    q4_nodes = copy.deepcopy(loop_nodes)
    q4_nodes.extend([
        {'operation': 'fixture_Qy', 'order': 0},
        {'operation': 'fixture_Qy', 'order': 1},
        {'operation': 'product', 'arguments': [6, 4]},
        {'operation': 'product', 'arguments': [7, 4]},
        {'operation': 'function_substitution', 'expression': 8, 'variable': 0,
         'value': 3, 'Z_independent_substitution': True},
        {'operation': 'function_substitution', 'expression': 9, 'variable': 0,
         'value': 3, 'Z_independent_substitution': True},
    ])
    q4 = {'source_family': family, 'exact_graph_nodes': q4_nodes,
          'modified_Q4_pressure_C3_bindings': {
              'initial': {'R_Q_y': [10, 11]},
              'charts': {'first_micro': {
                  'finite_boundary': {'R_Q_y': [8, 9]},
                  'physical': {'Q_y': [6, 7]},
              }},
          }}
    loop = {'source_family': family, 'exact_graph_nodes': loop_nodes,
            'prescribed_loop_C5_bindings': {
                'sc': 2,
                'charts': {'first_micro': {'coordinate': 0, 'domain': [3, 5], 'radius': 4}},
            }}
    return q4, loop


def run_fixture(verify):
    checks = {}
    diagnostics = {}
    q4, loop = example()
    try:
        result = verify(q4, loop)
        checks['positive_archived_substitution'] = result['status'] == 'verified' and result['passed'] is True
        checks['both_ordinary_orders'] = result['scope']['ordinary_Z_orders'] == [0, 1]
        checks['no_cross_owner_claim'] = not result['scope']['native_corrected_function_identity_proved']
        checks['no_full_prefix_claim'] = not result['scope']['full_graph_prefix_equality_claimed']
        checks['sc_binding'] = result['evidence']['sc'] == loop['prescribed_loop_C5_bindings']['sc']
        checks['no_numeric_points'] = not result['scope']['numeric_values_computed']
    except Exception as exc:
        diagnostics['positive_exception'] = type(exc).__name__ + ': ' + str(exc)
        checks['positive_no_exception'] = False
    for case in ('wrong_half', 'wrong_sc_family', 'sc_depends_on_Z', 'wrong_lower',
                 'different_shared_record', 'wrong_expression', 'wrong_variable',
                 'wrong_value', 'wrong_radius', 'missing_domain'):
        q4, loop = example()
        if case == 'wrong_half':
            loop['exact_graph_nodes'][1]['denominator'] = 3
        elif case == 'wrong_sc_family':
            loop['exact_graph_nodes'][2]['source_family'] = {'fixture_family': 'other_B'}
        elif case == 'sc_depends_on_Z':
            loop['exact_graph_nodes'][2]['Z_independent'] = False
        elif case == 'wrong_lower':
            loop['exact_graph_nodes'][3]['arguments'] = [1, 2]
        elif case == 'different_shared_record':
            q4['exact_graph_nodes'][0]['name'] = 'different_coordinate_B'
        elif case == 'wrong_expression':
            q4['exact_graph_nodes'][10]['expression'] = 9
        elif case == 'wrong_variable':
            q4['exact_graph_nodes'][11]['variable'] = 5
        elif case == 'wrong_value':
            q4['exact_graph_nodes'][10]['value'] = 5
        elif case == 'wrong_radius':
            q4['exact_graph_nodes'][8]['arguments'] = [6, 5]
        elif case == 'missing_domain':
            del loop['prescribed_loop_C5_bindings']['charts']['first_micro']['domain']
        try:
            bad = verify(q4, loop)
            checks['reject_' + case] = bad.get('passed') is False and bad.get('status') in ('mismatch', 'unresolved')
        except Exception as exc:
            checks['reject_' + case] = False
            diagnostics[case + '_exception'] = type(exc).__name__ + ': ' + str(exc)
    return {'passed': bool(checks) and all(checks.values()), 'checks': checks,
            'diagnostics': diagnostics, 'synthetic_only': True}


__all__ = ['run_fixture']
