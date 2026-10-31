"""Bind archived first-micro flat-loop rows to the sc/2 endpoint.

No source provider is imported or replayed.  The graph result is conditional
on sc > 0; it is not a physical point value or a cross-owner identity.
"""
from __future__ import annotations

from collections.abc import Mapping

SCHEMA = 'original_q4_archived_flat_collar_v1'


class Missing(ValueError):
    pass


class Mismatch(ValueError):
    pass


def f(obj, key, label):
    if not isinstance(obj, Mapping) or key not in obj:
        raise Missing(label)
    return obj[key]


def n(nodes, index, label):
    if not isinstance(nodes, list) or type(index) is not int or not 0 <= index < len(nodes):
        raise Missing(label)
    row = nodes[index]
    if not isinstance(row, Mapping):
        raise Missing(label)
    return row


def rows(obj, label):
    if not isinstance(obj, list) or len(obj) < 3 or any(type(k) is not int for k in obj[:3]):
        raise Missing(label)
    return obj[:3]


def verify(q4_archive, loop_archive):
    checks, evidence, errors = {}, {}, []

    def check(label, test):
        checks[label] = bool(test)
        if not test:
            raise Mismatch(label)

    try:
        qnodes = f(q4_archive, 'exact_graph_nodes', 'q4.nodes')
        lnodes = f(loop_archive, 'exact_graph_nodes', 'loop.nodes')
        family = f(loop_archive, 'source_family', 'loop.family')
        check('family_metadata_consistent', f(q4_archive, 'source_family', 'q4.family') == family)
        lb = f(loop_archive, 'prescribed_loop_C5_bindings', 'loop.bindings')
        lc = f(f(lb, 'charts', 'loop.charts'), 'first_micro', 'loop.first_micro')
        qc = f(f(f(q4_archive, 'modified_Q4_pressure_C3_bindings', 'q4.bindings'),
                 'charts', 'q4.charts'), 'first_micro', 'q4.first_micro')
        sc = f(lb, 'sc', 'loop.sc')
        N = f(lb, 'N', 'loop.N')
        coordinate = f(lc, 'coordinate', 'loop.coordinate')
        domain = f(lc, 'domain', 'loop.domain')
        if not isinstance(domain, list) or len(domain) < 1:
            raise Missing('loop.domain.lower')
        lower = domain[0]
        lower_row = n(lnodes, lower, 'lower')
        arguments = lower_row.get('arguments')
        if not isinstance(arguments, list) or len(arguments) != 2:
            raise Missing('lower.arguments')
        half = arguments[1]
        check('lower_is_sc_over_two',
              lower_row.get('operation') == 'product' and arguments[0] == sc
              and n(lnodes, half, 'half') == {'operation': 'exact_rational', 'numerator': 1, 'denominator': 2})
        sc_row = n(lnodes, sc, 'sc')
        check('sc_is_same_Z_independent_source',
              sc_row.get('operation') == 'current_original_source_parameter'
              and sc_row.get('name') == 'sc' and sc_row.get('source_family') == family
              and sc_row.get('Z_independent') is True)
        N_row = n(lnodes, N, 'N')
        check('N_is_common_positive_integer_expression',
              N_row.get('operation') == 'exact_positive_integer_expression'
              and N_row.get('one_common_N_for_all_charts_and_control_denominators') is True)
        loop = f(lc, 'loop', 'loop.functions')
        refs = {key: rows(f(loop, key, 'loop.' + key), 'loop.' + key + '.rows')
                for key in ('A', 'B_over_Pstar', 'A_phi', 'B_over_Pstar_phi')}
        first = n(lnodes, refs['A'][0], 'A[0]')
        collar = first.get('original_collar_upper')
        predicate = first.get('flat_predicate')
        collar_row = n(lnodes, collar, 'collar')
        args = collar_row.get('arguments')
        if not isinstance(args, list) or len(args) != 2:
            raise Missing('collar.arguments')
        three_quarters = args[0]
        fraction = n(lnodes, three_quarters, 'three quarters')
        check('collar_is_three_quarters_sc',
              collar_row.get('operation') == 'product' and args[1] == sc
              and fraction == {'operation': 'exact_rational', 'numerator': 3, 'denominator': 4})
        pred = n(lnodes, predicate, 'predicate')
        pred_args = pred.get('arguments')
        if pred.get('operation') != 'logical_or' or not isinstance(pred_args, list):
            raise Mismatch('flat predicate logical_or')
        comparison_found = any(
            (lambda row: row.get('operation') == 'exact_real_comparison'
             and row.get('operator') == '<=' and row.get('left') == coordinate
             and row.get('right') == collar)(n(lnodes, ref, 'predicate branch'))
            for ref in pred_args
        )
        check('predicate_contains_coordinate_le_collar', comparison_found)
        check('flat_graph_zero_exact',
              n(lnodes, 0, 'zero') == {'operation': 'exact_rational', 'numerator': 0, 'denominator': 1})
        for key, ids in refs.items():
            for k, index in enumerate(ids):
                row = n(lnodes, index, key + '[' + str(k) + ']')
                evidence['last_inspected'] = {'field': key, 'order': k, 'id': index, 'node': row}
                check(key + '_Z' + str(k) + '_same_lazy_flat_branch',
                      row.get('operation') == 'original_lazy_flat_branch'
                      and row.get('flat_value') == 0
                      and row.get('flat_predicate') == predicate
                      and row.get('original_collar_coordinate') == coordinate
                      and row.get('original_collar_upper') == collar
                      and row.get('ordinary_slow_Z_derivative_order') == k
                      and row.get('same_original_flat_extension_through_order5') is True
                      and row.get('active_body_not_evaluated_on_flat_branch') is True)
                check(key + '_Z' + str(k) + '_node_matches_q4_archive',
                      row == n(qnodes, index, 'q4.' + key))
        # The finite-N velocity correction is B_over_Pstar/N.
        finite = f(lc, 'finite_N', 'loop.finite_N')
        delta_v = rows(f(finite, 'delta_V', 'finite_N.delta_V'), 'finite_N.delta_V')
        total_v = rows(f(finite, 'V_N', 'finite_N.V_N'), 'finite_N.V_N')
        base_v = rows(f(f(lc, 'roots', 'loop.roots'), 'V', 'roots.V'), 'roots.V.rows')
        invN = None
        for k, index in enumerate(delta_v):
            row = n(lnodes, index, 'delta_V')
            args = row.get('arguments')
            evidence['last_inspected'] = {'field': 'delta_V', 'order': k, 'id': index, 'node': row}
            if row.get('operation') != 'product' or not isinstance(args, list) or len(args) != 2:
                raise Mismatch('delta_V product Z' + str(k))
            if invN is None:
                invN = args[1]
            check('delta_V_Z' + str(k) + '_is_B_times_common_invN',
                  args == [refs['B_over_Pstar'][k], invN])
            total = n(lnodes, total_v[k], 'V_N')
            check('V_N_Z' + str(k) + '_retains_background_plus_correction',
                  total.get('operation') == 'sum' and total.get('arguments') == [base_v[k], index])
        inv = n(lnodes, invN, 'inverse N')
        check('invN_references_same_N',
              inv.get('operation') == 'positive_quotient' and inv.get('numerator') == 1
              and inv.get('denominator') == N)
        # Compare every directly used condition/scale/root node, not just its
        # integer label, across the two authenticated archive namespaces.
        for index in (0, 1, sc, half, lower, coordinate, collar, three_quarters,
                      predicate, N, invN, *pred_args, *delta_v, *total_v, *base_v):
            check('selected_shared_node_' + str(index),
                  n(qnodes, index, 'q4.selected') == n(lnodes, index, 'loop.selected'))
        history = f(lc, 'history', 'loop.history')
        incoming = f(lc, 'incoming', 'loop.incoming')
        hist_m = rows(f(history, 'm', 'history.m'), 'history.m rows')
        incoming_m = rows(f(incoming, 'm', 'incoming.m'), 'incoming.m rows')
        evidence = {
            'sc': sc, 'half': half, 'lower': lower, 'coordinate': coordinate,
            'collar': collar, 'predicate': predicate, 'common_N_node': N,
            'lazy_source_rows': refs, 'delta_V_rows': delta_v,
            'background_V_rows': base_v, 'V_N_rows': total_v,
            'm_history_rows': hist_m, 'm_incoming_correction_rows': incoming_m,
            'q4_radial_available': isinstance(qc.get('radial'), Mapping),
        }
        status = 'verified'
    except Missing as exc:
        status = 'unresolved'
        errors.append(str(exc))
    except Mismatch as exc:
        status = 'mismatch'
        errors.append(str(exc))
    return {
        'schema': SCHEMA, 'status': status, 'passed': status == 'verified',
        'checks': checks, 'evidence': evidence, 'errors': errors,
        'scope': {
            'flat_at_sc_over_two_conditional_on_sc_positive': status == 'verified',
            'sc_positive_independently_certified': False,
            'delta_V_zero_at_endpoint_conditional': status == 'verified',
            'Q4_total_Uy_correction_zero_proved': False,
            'm_history_equal_across_owners_proved': False,
            'corrected_source_identity_proved': False,
            'numeric_values_computed': False,
            'old_errors_removed': False,
            'ordinary_Z_orders': [0, 1, 2],
        },
    }
