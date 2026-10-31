"""Fresh synthetic tests of Q4 archived physical row source binding."""
from __future__ import annotations
import copy


def example():
    family = {'fixture': 'physical_rows_distinct_family_B'}
    loop_nodes = [
        {'operation': 'exact_rational', 'numerator': 0, 'denominator': 1},
        {'operation': 'fixture_logP', 'name': 'fresh_log_scale_B'},
        {'operation': 'analytic_unary', 'name': 'exp', 'argument': 1},
    ]
    def add(row):
        loop_nodes.append(row)
        return len(loop_nodes) - 1
    v = [add({'operation': 'fixture_V', 'order': k}) for k in range(3)]
    m = [add({'operation': 'fixture_m', 'order': k}) for k in range(3)]
    h = [add({'operation': 'fixture_history', 'order': k}) for k in range(3)]
    incoming = [add({'operation': 'fixture_incoming', 'order': k}) for k in range(3)]
    qnodes = copy.deepcopy(loop_nodes)
    def qadd(row):
        qnodes.append(row)
        return len(qnodes)-1
    vy = [qadd({'operation': 'fixture_radial_Vy', 'order': k}) for k in range(3)]
    phys = {key: [qadd({'operation': 'product', 'arguments': [i, 2]}) for i in source]
            for key, source in [('U', v), ('m', m), ('U_y', vy)]}
    q4 = {'source_family': family, 'exact_graph_nodes': qnodes,
          'modified_Q4_pressure_C3_bindings': {'charts': {'first_micro': {
              'physical': phys, 'radial': {'V_y': vy}}}}}
    loop = {'source_family': family, 'exact_graph_nodes': loop_nodes,
            'prescribed_loop_C5_bindings': {'charts': {'first_micro': {
                'modified_fields': {'V': v, 'm': m}, 'history': {'m': h},
                'incoming': {'m': incoming}}}}}
    return q4, loop


def run_fixture(verify):
    checks, diagnostics = {}, {}
    q4, loop = example()
    try:
        result = verify(q4, loop)
        checks['positive_physical_rows'] = result['status'] == 'verified' and result['passed'] is True
        checks['all_ordinary_rows'] = result['scope']['ordinary_Z_orders'] == [0, 1, 2]
        checks['history_sources_retained'] = len(result['evidence']['m_history_rows']) == 3
        checks['incoming_sources_retained'] = len(result['evidence']['m_incoming_rows']) == 3
        checks['no_cross_owner_identity'] = not result['scope']['native_corrected_function_identity_proved']
        checks['no_point_values'] = not result['scope']['numeric_values_computed']
    except Exception as exc:
        checks['positive_no_exception'] = False
        diagnostics['positive_exception'] = type(exc).__name__ + ': ' + str(exc)
    for case in ('prefix_change', 'wrong_U_scale', 'wrong_m_source', 'wrong_Uy_scale',
                 'wrong_scale_kind', 'missing_history', 'short_Uy', 'family_mismatch'):
        q4, loop = example()
        physical = q4['modified_Q4_pressure_C3_bindings']['charts']['first_micro']['physical']
        if case == 'prefix_change':
            q4['exact_graph_nodes'][3]['order'] = 97
        elif case == 'wrong_U_scale':
            q4['exact_graph_nodes'][physical['U'][1]]['arguments'][1] = 1
        elif case == 'wrong_m_source':
            q4['exact_graph_nodes'][physical['m'][2]]['arguments'][0] = 3
        elif case == 'wrong_Uy_scale':
            q4['exact_graph_nodes'][physical['U_y'][0]]['arguments'][1] = 1
        elif case == 'wrong_scale_kind':
            q4['exact_graph_nodes'][2]['name'] = 'sin'
            loop['exact_graph_nodes'][2]['name'] = 'sin'
        elif case == 'missing_history':
            del loop['prescribed_loop_C5_bindings']['charts']['first_micro']['history']
        elif case == 'short_Uy':
            physical['U_y'] = physical['U_y'][:1]
        elif case == 'family_mismatch':
            loop['source_family'] = {'fixture': 'other_family_C'}
        try:
            bad = verify(q4, loop)
            checks['reject_' + case] = bad.get('passed') is False and bad.get('status') in ('mismatch','unresolved')
        except Exception as exc:
            checks['reject_' + case] = False
            diagnostics[case + '_exception'] = type(exc).__name__ + ': ' + str(exc)
    return {'passed': bool(checks) and all(checks.values()), 'checks': checks,
            'diagnostics': diagnostics, 'synthetic_only': True}


__all__ = ['run_fixture']
