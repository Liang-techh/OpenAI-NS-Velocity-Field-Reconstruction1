"""Bind archived Q4 slow B derivatives to the admitted flat predicate.

Previous frozen receipts are decoded by the caller after raw authentication.
This source-only check does not invoke their helpers, evaluate a field, prove
sc positive, or identify any corrected source.
"""
from __future__ import annotations
from collections.abc import Mapping

SCHEMA = 'original_q4_slow_B_flat_v1'
Q4 = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz'
LOOP = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz'
RAW_SHA = {
    Q4: 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a',
    LOOP: '8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa',
}


class Missing(ValueError):
    pass


class Mismatch(ValueError):
    pass


def get(obj, key, label):
    if not isinstance(obj, Mapping) or key not in obj:
        raise Missing(label)
    return obj[key]


def nodes_row(nodes, index, label):
    if not isinstance(nodes, list) or type(index) is not int or not 0 <= index < len(nodes):
        raise Missing(label)
    value = nodes[index]
    if not isinstance(value, Mapping):
        raise Missing(label)
    return value


def first3(value, label):
    if not isinstance(value, list) or len(value) < 3 or any(type(x) is not int for x in value[:3]):
        raise Missing(label)
    return value[:3]


def receipt_data(bundle, label):
    result = get(bundle, 'result', label + '.result')
    archive = get(result, 'archive_result', label + '.archive_result')
    if get(archive, 'status', label + '.status') != 'verified':
        raise Mismatch(label + ' is not verified')
    hashes = get(result, 'input_hashes', label + '.input_hashes')
    return archive, hashes


def verify(q4, loop, flat_receipt, physical_receipt):
    checks, evidence, errors = {}, {}, []

    def check(label, value):
        checks[label] = bool(value)
        if not value:
            raise Mismatch(label)

    try:
        qnodes = get(q4, 'exact_graph_nodes', 'q4.nodes')
        lnodes = get(loop, 'exact_graph_nodes', 'loop.nodes')
        qfamily = get(q4, 'source_family', 'q4.family')
        lfamily = get(loop, 'source_family', 'loop.family')
        check('archive_families_match', qfamily == lfamily)
        flat, fh = receipt_data(flat_receipt, 'flat')
        physical, ph = receipt_data(physical_receipt, 'physical')
        # Bind both earlier archived results to these exact raw archive hashes.
        check('prior_receipt_input_hashes_match',
              all(fh.get(k) == ph.get(k) == RAW_SHA[k] for k in (Q4, LOOP)))
        flat_scope = get(flat, 'scope', 'flat.scope')
        check('prior_flat_is_conditional',
              flat_scope.get('flat_at_sc_over_two_conditional_on_sc_positive') is True
              and flat_scope.get('sc_positive_independently_certified') is False)
        fe = get(flat, 'evidence', 'flat.evidence')
        pe = get(physical, 'evidence', 'physical.evidence')
        predicate = get(fe, 'predicate', 'flat.predicate')
        coordinate = get(fe, 'coordinate', 'flat.coordinate')
        N = get(fe, 'common_N_node', 'flat.N')
        collar = get(fe, 'collar', 'flat.collar')
        bphi = first3(get(get(fe, 'lazy_source_rows', 'flat.lazy'), 'B_over_Pstar_phi', 'flat.Bphi'), 'flat.Bphi rows')
        lb = get(loop, 'prescribed_loop_C5_bindings', 'loop.bindings')
        lc = get(get(lb, 'charts', 'loop.charts'), 'first_micro', 'loop.first_micro')
        check('same_coordinate_and_N_nodes', lc.get('coordinate') == coordinate and lb.get('N') == N)
        jac = get(lc, 'physical_measure', 'loop.physical_measure')
        loop_Bphi = first3(get(get(lc, 'loop', 'loop.loop'), 'B_over_Pstar_phi', 'loop.Bphi'), 'loop.Bphi rows')
        check('Bphi_rows_are_same_archived_sources', loop_Bphi == bphi)
        qb = get(q4, 'modified_Q4_pressure_C3_bindings', 'q4.bindings')
        qc = get(get(qb, 'charts', 'q4.charts'), 'first_micro', 'q4.first_micro')
        radial = get(qc, 'radial', 'q4.radial')
        background = get(qc, 'radial_background', 'q4.radial_background')
        bg = first3(get(background, 'V', 'background.V'), 'background.V rows')
        slow = first3(get(radial, 'slow_B_y', 'slow_B_y'), 'slow_B_y rows')
        total = first3(get(radial, 'V_y', 'radial.V_y'), 'radial.V_y rows')
        out_u = first3(get(get(qc, 'physical', 'q4.physical'), 'U_y', 'physical.U_y'), 'physical.U_y rows')
        normalized = first3(get(get(pe, 'source_normalized_rows', 'physical.source'), 'U_y', 'physical.source.U_y'), 'physical.source.U_y rows')
        saved_u = first3(get(get(pe, 'physical_rows', 'physical.rows'), 'U_y', 'physical.U_y'), 'physical.U_y rows')
        check('prior_physical_receipt_binds_actual_Uy_rows', total == normalized and out_u == saved_u)
        deriv_ids, inverse_jacobian, inverse_N = [], None, None
        for k in range(3):
            row = nodes_row(qnodes, slow[k], 'slow_B_y')
            evidence['last_inspected'] = {'field': 'slow_B_y', 'order': k, 'node': slow[k], 'row': row}
            args = row.get('arguments')
            if row.get('operation') != 'product' or not isinstance(args, list) or len(args) != 2:
                raise Mismatch('slow_B_y product Z' + str(k))
            jinv, derivative = args
            if inverse_jacobian is None:
                inverse_jacobian = jinv
            check('same_inverse_jacobian_Z' + str(k), jinv == inverse_jacobian)
            drow = nodes_row(qnodes, derivative, 'native derivative lazy branch')
            deriv_ids.append(derivative)
            check('slow_derivative_flat_branch_Z' + str(k),
                  drow.get('operation') == 'original_lazy_flat_branch'
                  and drow.get('flat_value') == 0
                  and drow.get('flat_predicate') == predicate
                  and drow.get('original_collar_coordinate') == coordinate
                  and drow.get('original_collar_upper') == collar
                  and drow.get('ordinary_slow_Z_derivative_order') == k
                  and drow.get('same_original_lazy_native_derivative') is True
                  and drow.get('same_original_flat_extension_through_order5') is True
                  and drow.get('active_body_not_evaluated_on_flat_branch') is True)
            total_row = nodes_row(qnodes, total[k], 'total radial V_y')
            evidence['last_inspected'] = {'field': 'radial.V_y', 'order': k, 'node': total[k], 'row': total_row}
            targs = total_row.get('arguments')
            if total_row.get('operation') != 'sum' or not isinstance(targs, list) or len(targs) != 3:
                raise Mismatch('total radial V_y sum Z' + str(k))
            check('total_has_background_and_Bphi_Z' + str(k),
                  targs[0] == bg[k] and targs[2] == bphi[k])
            scaled = nodes_row(qnodes, targs[1], 'scaled slow B')
            sargs = scaled.get('arguments')
            if scaled.get('operation') != 'product' or not isinstance(sargs, list) or len(sargs) != 2:
                raise Mismatch('scaled slow B product Z' + str(k))
            if inverse_N is None:
                inverse_N = sargs[1]
            check('same_slowB_over_N_Z' + str(k),
                  sargs == [slow[k], inverse_N])
        jrow = nodes_row(qnodes, inverse_jacobian, 'inverse Jacobian')
        check('inverse_jacobian_matches_chart_measure',
              jrow.get('operation') == 'positive_quotient'
              and jrow.get('numerator') == 1 and jrow.get('denominator') == jac)
        nrow = nodes_row(qnodes, inverse_N, 'inverse N')
        check('inverse_N_matches_common_N',
              nrow.get('operation') == 'positive_quotient'
              and nrow.get('numerator') == 1 and nrow.get('denominator') == N)
        # Actual Bphi graph rows were checked in the prior flat receipt.
        for k, ref in enumerate(bphi):
            check('actual_Bphi_node_matches_loop_Z' + str(k),
                  nodes_row(qnodes, ref, 'q4.Bphi') == nodes_row(lnodes, ref, 'loop.Bphi'))
        evidence = {
            'coordinate': coordinate, 'collar': collar, 'predicate': predicate,
            'common_N_node': N, 'jacobian_node': jac, 'inverse_N_node': inverse_N,
            'inverse_jacobian_node': inverse_jacobian,
            'slow_B_y_rows': slow, 'slow_derivative_lazy_rows': deriv_ids,
            'background_V_y_rows': bg, 'B_phi_rows': bphi,
            'total_V_y_rows': total, 'physical_U_y_rows': out_u,
            'Pstar_node_from_prior_receipt': get(pe, 'Pstar_candidate_node', 'physical.Pstar'),
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
            'Q4_Uy_correction_zero_at_endpoint_conditional_on_sc_positive': status == 'verified',
            'sc_positive_independently_certified': False,
            'background_Uy_retained': status == 'verified',
            'corrected_source_identity_proved': False,
            'history_function_equality_claimed': False,
            'numeric_values_computed': False,
            'old_errors_removed': False,
            'ordinary_Z_orders': [0, 1, 2],
        },
    }
