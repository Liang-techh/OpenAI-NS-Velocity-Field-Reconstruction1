"""Source-only Q4/loop physical-row binding on decoded archived graphs.

Caller authenticates input bytes.  This module reads archived nodes and never
imports a producer, evaluates source values or equates corrected endpoints.
"""
from __future__ import annotations

from collections.abc import Mapping

SCHEMA = 'original_q4_archived_physical_rows_v1'


class Missing(ValueError):
    pass


class Mismatch(ValueError):
    pass


def field(obj, key, label):
    if not isinstance(obj, Mapping) or key not in obj:
        raise Missing(label)
    return obj[key]


def rows(value, label):
    if not isinstance(value, list) or len(value) < 3 or any(type(x) is not int for x in value[:3]):
        raise Missing(label)
    return value[:3]


def node(nodes, index, label):
    if not isinstance(nodes, list) or type(index) is not int or not 0 <= index < len(nodes):
        raise Missing(label)
    item = nodes[index]
    if not isinstance(item, Mapping):
        raise Missing(label)
    return item


def is_zero(row):
    return row.get('operation') == 'exact_rational' and row.get('numerator') == 0 and row.get('denominator') == 1


def verify(q4_archive, loop_archive):
    checks, evidence, errors = {}, {}, []

    def check(label, condition):
        checks[label] = bool(condition)
        if not condition:
            raise Mismatch(label)

    try:
        qnodes = field(q4_archive, 'exact_graph_nodes', 'q4.graph')
        lnodes = field(loop_archive, 'exact_graph_nodes', 'loop.graph')
        if not isinstance(qnodes, list) or not isinstance(lnodes, list):
            raise Missing('archived graph arrays')
        evidence['graph_sizes'] = {'q4': len(qnodes), 'loop': len(lnodes)}
        qfamily = field(q4_archive, 'source_family', 'q4.family')
        lfamily = field(loop_archive, 'source_family', 'loop.family')
        check('family_metadata_matches', qfamily == lfamily)
        check('loop_archive_is_exact_q4_graph_prefix',
              len(lnodes) <= len(qnodes) and all(qnodes[k] == row for k, row in enumerate(lnodes)))
        qb = field(q4_archive, 'modified_Q4_pressure_C3_bindings', 'q4.bindings')
        lb = field(loop_archive, 'prescribed_loop_C5_bindings', 'loop.bindings')
        qc = field(field(qb, 'charts', 'q4.charts'), 'first_micro', 'q4.first_micro')
        lc = field(field(lb, 'charts', 'loop.charts'), 'first_micro', 'loop.first_micro')
        physical = field(qc, 'physical', 'q4.physical')
        radial = field(qc, 'radial', 'q4.radial')
        normalized = field(lc, 'modified_fields', 'loop.modified_fields')
        history = field(lc, 'history', 'loop.history')
        incoming = field(lc, 'incoming', 'loop.incoming')
        source = {
            'U': rows(field(normalized, 'V', 'loop.V'), 'loop.V rows'),
            'm': rows(field(normalized, 'm', 'loop.m'), 'loop.m rows'),
            'U_y': rows(field(radial, 'V_y', 'q4.radial.V_y'), 'q4.radial.V_y rows'),
        }
        physical_rows = {key: rows(field(physical, key, 'q4.physical.' + key), key + ' rows')
                         for key in ('U', 'm', 'U_y')}
        history_m = rows(field(history, 'm', 'loop.history.m'), 'loop.history.m rows')
        incoming_m = rows(field(incoming, 'm', 'loop.incoming.m'), 'loop.incoming.m rows')
        evidence.update(source_normalized_rows=source, physical_rows=physical_rows,
                        m_history_rows=history_m, m_incoming_rows=incoming_m)
        # A nonzero physical product supplies the common scale candidate.
        candidate = None
        for key in ('U', 'm', 'U_y'):
            for src_id, out_id in zip(source[key], physical_rows[key]):
                src, out = node(qnodes, src_id, key + ' source'), node(qnodes, out_id, key + ' physical')
                evidence['last_inspected'] = {'field': key, 'source_id': src_id,
                                              'physical_id': out_id,
                                              'source_node': src, 'physical_node': out}
                if is_zero(src):
                    continue
                args = out.get('arguments')
                if out.get('operation') != 'product' or not isinstance(args, list) or len(args) != 2 or args[0] != src_id:
                    raise Mismatch('physical product source ' + key)
                candidate = args[1] if candidate is None else candidate
                check('common_scale_' + key + '_' + str(len(checks)), args[1] == candidate)
        if candidate is None:
            raise Missing('no nonzero source row to identify scale')
        scale = node(qnodes, candidate, 'Pstar candidate')
        check('common_scale_is_exponential_source',
              scale.get('operation') == 'analytic_unary' and scale.get('name') == 'exp'
              and type(scale.get('argument')) is int)
        for key in ('U', 'm', 'U_y'):
            for k, (src_id, out_id) in enumerate(zip(source[key], physical_rows[key])):
                src, out = node(qnodes, src_id, key + ' source'), node(qnodes, out_id, key + ' physical')
                if is_zero(src):
                    check(key + '_Z' + str(k) + '_scale', is_zero(out))
                else:
                    check(key + '_Z' + str(k) + '_scale',
                          out.get('operation') == 'product' and out.get('arguments') == [src_id, candidate])
        # Record history operands and exact source rows; do not substitute zero.
        history_links = []
        for k, (source_id, hist_id) in enumerate(zip(source['m'], history_m)):
            sr = node(qnodes, source_id, 'm normalized row')
            hr = node(qnodes, hist_id, 'm history row')
            args = sr.get('arguments', [])
            history_links.append({
                'order': k, 'modified_m_node': source_id, 'history_node': hist_id,
                'history_is_explicit_sum_operand': sr.get('operation') == 'sum' and hist_id in args,
                'history_is_exact_graph_zero': is_zero(hr),
            })
        evidence = {
            'Pstar_candidate_node': candidate, 'Pstar_exponential_argument_node': scale['argument'],
            'source_normalized_rows': source, 'physical_rows': physical_rows,
            'm_history_rows': history_m, 'm_incoming_rows': incoming_m,
            'm_history_links': history_links, 'loop_prefix_node_count': len(lnodes),
            'q4_node_count': len(qnodes),
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
            'archived_q4_loop_graph_prefix_only': True,
            'Pstar_source_identity_beyond_native_graph_proved': False,
            'native_corrected_function_identity_proved': False,
            'live_owner_identity_proved': False,
            'history_function_equality_claimed': False,
            'ordinary_Z_orders': [0, 1, 2],
            'numeric_values_computed': False,
            'old_errors_removed': False,
        },
    }
