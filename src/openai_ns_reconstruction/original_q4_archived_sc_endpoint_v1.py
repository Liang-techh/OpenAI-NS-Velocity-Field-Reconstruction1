"""Read selected archived Q4/loop source rows without replaying producers.

The caller authenticates raw archive bytes before decoding.  This verifier
compares only the named source rows, not a complete graph prefix, numerical
values, live-object identity or the corrected endpoint.
"""
from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping

SCHEMA = 'original_q4_archived_sc_endpoint_v1'


class Missing(ValueError):
    pass


class Mismatch(ValueError):
    pass


def need(value, label):
    if not value:
        raise Mismatch(label)


def field(obj, key, label):
    if not isinstance(obj, Mapping) or key not in obj:
        raise Missing(label)
    return obj[key]


def node(nodes, index, label):
    if not isinstance(nodes, list) or type(index) is not int or not 0 <= index < len(nodes):
        raise Missing(label)
    item = nodes[index]
    if not isinstance(item, Mapping):
        raise Missing(label)
    return item


def pair(obj, label):
    if not isinstance(obj, list) or len(obj) < 2 or any(type(x) is not int for x in obj[:2]):
        raise Missing(label)
    return obj[:2]


def digest(item):
    raw = json.dumps(item, sort_keys=True, separators=(',', ':'), ensure_ascii=True).encode()
    return hashlib.sha256(raw).hexdigest()


def verify(q4_archive, loop_archive):
    checks = {}
    evidence = {}
    errors = []

    def check(label, condition):
        checks[label] = bool(condition)
        need(condition, label)

    try:
        qnodes = field(q4_archive, 'exact_graph_nodes', 'q4.exact_graph_nodes')
        lnodes = field(loop_archive, 'exact_graph_nodes', 'loop.exact_graph_nodes')
        qb = field(q4_archive, 'modified_Q4_pressure_C3_bindings', 'q4.bindings')
        lb = field(loop_archive, 'prescribed_loop_C5_bindings', 'loop.bindings')
        qfamily = field(q4_archive, 'source_family', 'q4.source_family')
        lfamily = field(loop_archive, 'source_family', 'loop.source_family')
        check('archive_family_metadata_matches', qfamily == lfamily)
        lc = field(field(lb, 'charts', 'loop.charts'), 'first_micro', 'loop.first_micro')
        qc = field(field(qb, 'charts', 'q4.charts'), 'first_micro', 'q4.first_micro')
        sc = field(lb, 'sc', 'loop.sc')
        coord = field(lc, 'coordinate', 'loop.coordinate')
        dom = field(lc, 'domain', 'loop.domain')
        if not isinstance(dom, list) or len(dom) < 2:
            raise Missing('loop.domain[0]')
        lower = dom[0]
        radius = field(lc, 'radius', 'loop.radius')
        lowrow = node(lnodes, lower, 'loop.lower')
        args = lowrow.get('arguments')
        if not isinstance(args, list) or len(args) != 2:
            raise Missing('loop.lower.arguments')
        half = args[1]
        check('lower_is_exported_sc_times_half',
              lowrow.get('operation') == 'product' and args[0] == sc)
        halfrow = node(lnodes, half, 'loop.half')
        check('half_is_exact_one_half',
              halfrow.get('operation') == 'exact_rational'
              and halfrow.get('numerator') == 1 and halfrow.get('denominator') == 2)
        scrow = node(lnodes, sc, 'loop.sc')
        check('sc_is_family_Z_independent_source_parameter',
              scrow.get('operation') == 'current_original_source_parameter'
              and scrow.get('name') == 'sc'
              and scrow.get('source_family') == lfamily
              and scrow.get('Z_independent') is True)
        shared = {'sc': sc, 'half': half, 'lower': lower, 'coordinate': coord, 'radius': radius}
        for label, index in shared.items():
            left = node(lnodes, index, 'loop.' + label)
            right = node(qnodes, index, 'q4.' + label)
            check('selected_node_record_matches_' + label, left == right)
        initial = pair(field(field(qb, 'initial', 'q4.initial'), 'R_Q_y', 'q4.initial.R_Q_y'),
                       'q4.initial.R_Q_y rows')
        boundary = pair(field(field(qc, 'finite_boundary', 'q4.finite_boundary'),
                              'R_Q_y', 'q4.finite_boundary.R_Q_y'), 'q4.boundary rows')
        qy = pair(field(field(qc, 'physical', 'q4.physical'), 'Q_y', 'q4.physical.Q_y'),
                  'q4.physical.Q_y rows')
        for k in range(2):
            row = node(qnodes, initial[k], 'q4.initial Z' + str(k))
            check('initial_substitution_Z' + str(k), row == {
                'operation': 'function_substitution', 'expression': boundary[k],
                'variable': coord, 'value': lower, 'Z_independent_substitution': True})
            brow = node(qnodes, boundary[k], 'q4.boundary Z' + str(k))
            check('boundary_is_radius_times_Qy_Z' + str(k),
                  brow == {'operation': 'product', 'arguments': [qy[k], radius]})
        evidence = {
            'sc': sc, 'half': half, 'lower': lower, 'coordinate': coord,
            'radius': radius, 'initial_R_Q_y_Z0_Z1': initial,
            'finite_R_Q_y_Z0_Z1': boundary, 'Q_y_Z0_Z1': qy,
            'selected_shared_node_sha256': {key: digest(node(lnodes, index, key))
                                            for key, index in shared.items()},
            'initial_node_sha256': [digest(node(qnodes, index, 'initial')) for index in initial],
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
            'archived_selected_source_nodes_only': True,
            'full_graph_prefix_equality_claimed': False,
            'live_owner_identity_proved': False,
            'native_corrected_function_identity_proved': False,
            'numeric_values_computed': False,
            'endpoint_cancellation_claimed': False,
            'old_errors_removed': False,
            'ordinary_Z_orders': [0, 1],
        },
    }
