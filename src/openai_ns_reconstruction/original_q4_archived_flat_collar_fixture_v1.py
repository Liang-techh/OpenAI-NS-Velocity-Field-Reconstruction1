"""Fresh synthetic fixture for source-owned Q4 flat-collar graph."""
from __future__ import annotations
import copy


def example():
    fam = {'fixture': 'flat_collar_distinct_family_C'}
    g = [{'operation': 'exact_rational', 'numerator': 0, 'denominator': 1},
         {'operation': 'exact_rational', 'numerator': 1, 'denominator': 1}]
    def add(row):
        g.append(row)
        return len(g)-1
    sc = add({'operation': 'current_original_source_parameter', 'name': 'sc',
              'Z_independent': True, 'source_family': fam})
    half = add({'operation': 'exact_rational', 'numerator': 1, 'denominator': 2})
    lower = add({'operation': 'product', 'arguments': [sc, half]})
    three = add({'operation': 'exact_rational', 'numerator': 3, 'denominator': 4})
    collar = add({'operation': 'product', 'arguments': [three, sc]})
    coord = add({'operation': 'bound_variable', 'name': 'fresh_native_coordinate_C'})
    shear = add({'operation': 'fixture_shear'})
    comp = add({'operation': 'exact_real_comparison', 'operator': '<=', 'left': coord, 'right': collar})
    pred = add({'operation': 'logical_or', 'arguments': [shear, comp]})
    N = add({'operation': 'exact_positive_integer_expression',
             'one_common_N_for_all_charts_and_control_denominators': True})
    invN = add({'operation': 'positive_quotient', 'numerator': 1, 'denominator': N})
    loop = {}
    for key in ('A', 'B_over_Pstar', 'A_phi', 'B_over_Pstar_phi'):
        loop[key] = []
        for k in range(3):
            body = add({'operation': 'fixture_active_body', 'label': key, 'order': k})
            loop[key].append(add({'operation': 'original_lazy_flat_branch', 'active_body': body,
                                  'flat_value': 0, 'flat_predicate': pred,
                                  'original_collar_coordinate': coord, 'original_collar_upper': collar,
                                  'ordinary_slow_Z_derivative_order': k,
                                  'same_original_flat_extension_through_order5': True,
                                  'active_body_not_evaluated_on_flat_branch': True}))
    base = [add({'operation': 'fixture_background_velocity_C', 'order': k}) for k in range(3)]
    dv = [add({'operation': 'product', 'arguments': [b, invN]}) for b in loop['B_over_Pstar']]
    vn = [add({'operation': 'sum', 'arguments': [v, d]}) for v, d in zip(base, dv)]
    history = [add({'operation': 'fixture_background_history_C', 'order': k}) for k in range(3)]
    incoming = [0, 0, 0]
    loop_archive = {'source_family': fam, 'exact_graph_nodes': g,
                    'prescribed_loop_C5_bindings': {
                        'sc': sc, 'N': N, 'charts': {'first_micro': {
                            'coordinate': coord, 'domain': [lower, 1], 'loop': loop,
                            'roots': {'V': base}, 'finite_N': {'V_N': vn, 'delta_V': dv},
                            'history': {'m': history}, 'incoming': {'m': incoming},
                        }}}}
    q4 = {'source_family': fam, 'exact_graph_nodes': copy.deepcopy(g),
          'modified_Q4_pressure_C3_bindings': {'charts': {'first_micro': {'radial': {}}}}}
    return q4, loop_archive


def run_fixture(verify):
    checks, diagnostics = {}, {}
    q4, loop = example()
    try:
        out = verify(q4, loop)
        checks['positive_graph_binding'] = out['status'] == 'verified' and out['passed'] is True
        checks['sc_positive_remains_condition'] = out['scope']['sc_positive_independently_certified'] is False
        checks['conditional_delta_V_zero'] = out['scope']['delta_V_zero_at_endpoint_conditional'] is True
        checks['no_physical_Uy_claim'] = out['scope']['Q4_total_Uy_correction_zero_proved'] is False
        checks['ordinary_Z0_Z1_Z2'] = out['scope']['ordinary_Z_orders'] == [0, 1, 2]
        checks['history_handles_preserved'] = len(out['evidence']['m_history_rows']) == 3
        checks['no_corrected_identity'] = out['scope']['corrected_source_identity_proved'] is False
    except Exception as exc:
        checks['positive_no_exception'] = False
        diagnostics['positive_exception'] = type(exc).__name__ + ': ' + str(exc)
    for case in ('wrong_sc_family','wrong_half','wrong_collar','wrong_coordinate',
                 'wrong_lazy_order','wrong_flat_value','changed_q4_lazy','wrong_N',
                 'wrong_deltaV','wrong_background','missing_history'):
        q4, loop = example()
        b = loop['prescribed_loop_C5_bindings']
        lc = b['charts']['first_micro']
        nodes = loop['exact_graph_nodes']
        if case == 'wrong_sc_family':
            nodes[b['sc']]['source_family'] = {'fixture':'foreign_family_D'}
        elif case == 'wrong_half':
            nodes[nodes[lc['domain'][0]]['arguments'][1]]['denominator'] = 5
        elif case == 'wrong_collar':
            nodes[nodes[lc['loop']['A'][0]]['original_collar_upper']]['arguments'][1] = 1
        elif case == 'wrong_coordinate':
            nodes[nodes[lc['loop']['A'][0]]['flat_predicate']]['arguments'] = [8]
        elif case == 'wrong_lazy_order':
            nodes[lc['loop']['B_over_Pstar'][1]]['ordinary_slow_Z_derivative_order'] = 0
        elif case == 'wrong_flat_value':
            nodes[lc['loop']['A_phi'][2]]['flat_value'] = 1
        elif case == 'changed_q4_lazy':
            q4['exact_graph_nodes'][lc['loop']['A'][0]]['flat_value'] = 1
        elif case == 'wrong_N':
            nodes[b['N']]['one_common_N_for_all_charts_and_control_denominators'] = False
        elif case == 'wrong_deltaV':
            nodes[lc['finite_N']['delta_V'][1]]['arguments'][0] = 1
        elif case == 'wrong_background':
            nodes[lc['finite_N']['V_N'][2]]['arguments'][0] = 1
        else:
            del lc['history']['m']
        try:
            bad = verify(q4, loop)
            checks['reject_' + case] = bad.get('passed') is False and bad.get('status') in ('mismatch', 'unresolved')
        except Exception as exc:
            checks['reject_' + case] = False
            diagnostics[case + '_exception'] = type(exc).__name__ + ': ' + str(exc)
    return {'passed': bool(checks) and all(checks.values()), 'checks': checks,
            'diagnostics': diagnostics, 'synthetic_only': True}


__all__ = ['run_fixture']
