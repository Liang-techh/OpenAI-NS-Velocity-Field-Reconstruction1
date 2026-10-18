"""Constant-``U`` long-reshape and reference background ``I1`` arithmetic.

This module owns only the native constant-``U`` ``Q`` / pressure sector for
the long-reshape and reference charts.  It reads the current post-switch
``U`` and ``M`` rows, carries the incoming mean through both charts, and
keeps every huge radius in logarithmic form.  It does not evaluate the
accepted finite-``N`` loop or the native ``phi`` / ``A`` / ``B`` / ``C``
history providers.
"""

from __future__ import annotations

import gzip
import json
import math
from collections.abc import Mapping
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext

from .original_core_signed_direction import read, record, require, sha, source_guard


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_switch_post_finite_N_I1.json"
PARENT_SHA256 = "9d0d957bae129853030bad4cd5c484eb27e8b13b0623601e458fc4e95e93bc2b"
MODULE = "src/openai_ns_reconstruction/original_long_reference_background_I1_numeric.py"

C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
C6_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
Q4_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
Q4_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"

ORDERS = (0, 1, 2, 3, 4)
C1_ORDERS = (0, 1)
Q_ORDERS = (0, 1, 2)
I1_NAMES = ("q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z")
TAIL_BITS = 1024


def _lo(value):
    from mpmath import make_mpf
    return make_mpf(value._mpi_[0])


def _hi(value):
    from mpmath import make_mpf
    return make_mpf(value._mpi_[1])


def _point(c, value):
    if hasattr(value, "_mpi_"):
        low, high = _lo(value), _hi(value)
        require(low == high, "Point interval required")
        value = low
    value = c.mpf(value)
    return c.mpf([value, value])


def _read_record(c, value, *, directed=True):
    if isinstance(value, Mapping) and "value" in value and "lower_exact_mpf_tuple" not in value:
        value = value["value"]
    require(isinstance(value, Mapping), "Directed exact record required")
    require("lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value,
            "Directed exact tuple required")
    if directed:
        require(value.get("directed") is True, "Saved source record must be directed")
    out = read(c, value)
    require(_lo(out) <= _hi(out), "Malformed directed interval")
    return out


def _sum(c, values):
    out = c.mpf(0)
    for value in values:
        out = out + value
    return out


def _prod(c, left, right, order):
    return _sum(c, (c.mpf(math.comb(order, j)) * left[j] * right[order - j]
                    for j in range(order + 1)))


def _scale(c, scalar, value):
    return c.mpf(scalar) * value


def _neg(c, value):
    return c.mpf(0) - value


def _hull_interval(c, left, right):
    return c.mpf([min(_lo(left), _lo(right)), max(_hi(left), _hi(right))])


def _reciprocal(c, values, count):
    require(len(values) >= count, "Reciprocal rows are incomplete")
    require(_lo(values[0]) > 0, "Positive reciprocal denominator required")
    one = _point(c, 1)
    result = [one / values[0]]
    for order in range(1, count):
        tail = _sum(c, (c.mpf(math.comb(order, j)) * values[j]
                        * result[order - j] for j in range(1, order + 1)))
        result.append(_neg(c, tail / values[0]))
    return result


def _assert_lt(left, right, message):
    require(_hi(left) < _lo(right), message)


def _rows(c, rows, count, name):
    require(isinstance(rows, list) and len(rows) >= count,
            name + " must provide every ordinary row")
    found = {}
    for index, row in enumerate(rows):
        require(isinstance(row, Mapping), name + " row must be a mapping")
        order = row.get("ordinary_Z_order", index)
        require(type(order) is int and order in ORDERS,
                name + " has an invalid ordinary order")
        require(order not in found, name + " repeats an ordinary order")
        found[order] = _read_record(c, row.get("value", row))
    require(set(found) >= set(range(count)), name + " has missing ordinary rows")
    return [found[index] for index in range(count)]


def _pack_rows(rows, *, kind=None):
    result = []
    for order, value in enumerate(rows):
        item = {"ordinary_Z_order": order, "value": record(value)}
        if kind is not None:
            item["kind"] = kind
        result.append(item)
    return result


def _pin(root, path, expected):
    actual_path = root / path
    require(actual_path.is_file(), "Pinned source is missing: " + path)
    actual = sha(actual_path)
    require(actual == expected, "Pinned source changed: " + path)
    return dict(path=path, expected_sha256=expected, passed=True, source_executed=False)


def _load_parent(root, parent):
    path = root / PARENT
    require(path.is_file(), "Pinned switch-post parent is missing")
    require(sha(path) == PARENT_SHA256, "Pinned switch-post parent changed")
    disk = json.loads(path.read_bytes())
    require(parent == disk, "The exact accepted switch-post parent mapping is required")
    return disk


def _load_c6_cache(root):
    path = root / C6_CACHE
    require(path.is_file(), "Pinned C6 cache is missing")
    value = json.loads(gzip.open(path, "rb").read())
    require(isinstance(value, Mapping), "Pinned C6 cache must be a mapping")
    return value


def _source_pins(root, parent):
    expected = parent.get("input_hashes", {})
    paths = (C6_SOURCE, C6_CACHE, Q4_SOURCE, Q4_CACHE)
    pins = {PARENT: _pin(root, PARENT, PARENT_SHA256)}
    for path in paths:
        digest = expected.get(Path(path).name)
        require(digest is not None, "Parent does not pin " + path)
        pins[path] = _pin(root, path, digest)
    pins[MODULE] = _pin(root, MODULE, sha(root / MODULE))
    return pins


def _source_guards(root):
    guards = (
        source_guard(root, C6_SOURCE, "__init__", (
            "q=g.symbol('true_prescribed_long_reshape_fraction');y=g.mul(T,q);R=g.mul(c(110),g.unary('exp',y))",
            "actual=dict(phi=phi,U=incoming['actual']['U'],phi_y=a.mul(phi,slope),U_y=a.fixed(0)",
            "self.add_chart('long_reshape',q,R,T,actual,incoming['moments'])",
            "self.logref=logref=g.mul(c(10),g.add(self.core.parameters['logC'],self.core.parameters['logP']))",
            "self.reference_length=length=g.sub(g.sub(logref,T),c(8))",
            "q=g.symbol('true_prescribed_reference_fraction');y=g.add(T,g.mul(length,q));R=g.mul(c(110),g.unary('exp',y))",
            "actual=dict(phi=phi,U=incoming['actual']['U'],phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0)",
            "self.add_chart('reference',q,R,length,actual,incoming['moments'])",
        )),
        source_guard(root, C6_SOURCE, "add_chart", (
            "J=previous.histories(a,actual['phi'],actual['U'],incoming,phase,jacobian,'true_'+name)",
            "fields,radial,generic=previous.generic_fields(a,actual,J,self.F0,self.P0,self.Pstar,R,self.delta)",
        )),
        source_guard(root, Q4_SOURCE, "generic_functions", (
            "Q=div(a.sub(a.scale(a.mul(z,V),2),transport))",
        )),
        source_guard(root, Q4_SOURCE, "physical_pressure_fields", (
            "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
        )),
        source_guard(root, Q4_SOURCE, "reduced_pressure_S", (
            "linear=div(a.add(Q,a.scale(a.mul(a.sub(one,de),a.mul(z,Qz)),'1/2'),Qy))",
            "quadratic=a.mul(Q,a.add(a.scale(Q,'1/2'),Qy))",
        )),
    )
    require(isinstance(guards, tuple), "Source guards must remain a tuple")
    return guards


def _scalar_graph_guard(c6, parent):
    nodes = c6.get("exact_graph_nodes")
    bindings = c6.get("prescribed_R100_Rm_C6_bindings")
    require(isinstance(nodes, list) and isinstance(bindings, Mapping),
            "C6 exact graph and bindings are required")
    require(nodes[2] == {"operation": "exact_rational", "numerator": 40, "denominator": 1},
            "C6 node 2 must be the exact 40 scalar")
    require(nodes[3] == {"operation": "analytic_unary", "name": "exp", "argument": 2},
            "C6 node 3 must be exp(40)")
    require(nodes[4] == {"operation": "exact_rational", "numerator": 11, "denominator": 1},
            "C6 node 4 must be the exact 11 scalar")
    require(nodes[5] == {"operation": "sum", "arguments": [3, 4]},
            "C6 node 5 must be logP=exp(40)+11")
    require(nodes[2083] == {"operation": "exact_rational", "numerator": 10, "denominator": 1} and
            nodes[81513] == {"operation": "sum", "arguments": [35560, 5]},
            "C6 logref scalar factors must be exact 10 and true_core_logC+logP")
    require(nodes[6].get("operation") == "current_original_source_parameter" and
            nodes[6].get("name") == "logC" and nodes[6].get("Z_independent") is True and
            nodes[6].get("quantity_is_exact_original_parameter_not_range_endpoint") is True,
            "C6 node 6 must be the selected logC source parameter")
    require(nodes[7].get("operation") == "current_original_source_parameter" and
            nodes[7].get("name") == "T" and nodes[7].get("Z_independent") is True and
            nodes[7].get("quantity_is_exact_original_parameter_not_range_endpoint") is True,
            "C6 node 7 must be the exact T source parameter")
    native_logC = nodes[35560]
    require(native_logC.get("operation") == "current_original_source_parameter" and
            native_logC.get("name") == "true_core_logC" and
            native_logC.get("Z_independent") is True and
            native_logC.get("exact_fixed_source_parameter") is True and
            native_logC.get("source_path") == "common.core.logC" and
            native_logC.get("parameter_enclosure_not_selected_value") is True and
            isinstance(native_logC.get("source_provider"), Mapping),
            "C6 node 35560 must be the fixed true_core_logC source")
    require(nodes[6].get("source_family") == native_logC.get("source_family") and
            "selected_logCstar" in nodes[6].get("definition", "") and
            "exact singleton tuple" in nodes[6].get("definition", ""),
            "C6 native true_core_logC and selected logC source families must agree")
    selected = parent.get("actual_micro_mean_Z3_bounds", {}).get(
        "physical_scale_logs", {}).get("selected_logC")
    upstream = parent.get("n1_forward_actual_scale_jets", {}).get(
        "global_scale_records", {}).get("logCstar")
    require(isinstance(selected, Mapping) and isinstance(upstream, Mapping) and
            selected.get("directed") is True and
            upstream.get("directed_interval") is True and
            selected.get("lower_exact_mpf_tuple") ==
            selected.get("upper_exact_mpf_tuple") and
            upstream.get("lower_exact_mpf_tuple") ==
            upstream.get("upper_exact_mpf_tuple") and
            selected.get("lower_exact_mpf_tuple") ==
            upstream.get("lower_exact_mpf_tuple"),
            "C6 true_core_logC tuple must equal the parent selected_logC tuple")
    require(bindings.get("T") == 7 and bindings.get("logref") == 81514 and
            bindings.get("reference_length") == 81519,
            "C6 scalar binding indices changed")
    require(nodes[81514] == {"operation": "product", "arguments": [2083, 81513]},
            "C6 logref graph binding changed")
    require(nodes[81515] == {"operation": "negative", "argument": 7} and
            nodes[81516] == {"operation": "sum", "arguments": [81514, 81515]} and
            nodes[81517] == {"operation": "exact_rational", "numerator": 8, "denominator": 1} and
            nodes[81518] == {"operation": "negative", "argument": 81517} and
            nodes[81519] == {"operation": "sum", "arguments": [81516, 81518]},
            "C6 reference-length graph binding changed")
    definition = bindings.get("definition", {})
    require(definition.get("all_six_nonzero_incoming_and_every_chart_memory_retained") is True and
            definition.get("all_seven_source_rows_and_six_generic_output_rows_genuine") is True and
            definition.get("finite_N_modified_fields_correction_histories_and_all17_not_installed") is True,
            "C6 history/source-row contract changed")
    geometry = bindings.get("geometry_certificate", {})
    require(geometry.get("actual_native_T_and_logC_match_common_source") is True and
            geometry.get("all_four_native_T_source_tuples_equal") is True,
            "C6 native T/logC geometry identity changed")
    return dict(
        passed=True,
        node_indices=dict(logP=5, logC=6, T=7, logref=81514, reference_length=81519),
        node_operations=dict(
            logP="node 5 = node 3 + node 4 = exp(node 2) + 11",
            logC="node 35560 = current_original_source_parameter(true_core_logC); node 6 is the selected-logC alias",
            T="node 7 = current_original_source_parameter(T)",
            logref="node 81514 = node 2083 * node 81513 = 10*(logC+logP)",
            reference_length="node 81519 = logref - T - 8",
        ),
        history_contract=dict(
            all_six_nonzero_incoming_and_every_chart_memory_retained=True,
            all_seven_source_rows_and_six_generic_output_rows_genuine=True,
            finite_N_modified_fields_correction_histories_and_all17_not_installed=True,
            actual_native_T_and_logC_match_common_source=True,
            all_four_native_T_source_tuples_equal=True,
        ),
        true_core_logC_binding=dict(
            native_node=35560,
            selected_alias_node=6,
            source_path="common.core.logC",
            fixed_source_parameter=True,
            parent_selected_logC_tuple_equal=True,
            tuple_bound_from_parent_selected_logC=True,
            source_provider=dict(native_logC.get("source_provider")),
        ),
        source_executed=False,
    )


def _directions(c, parent):
    bg = parent["actual_switch_post_background_I1_numeric"]
    first = parent["actual_first_switch_background_I1_numeric"]
    require(bg.get("passed") is True and first.get("passed") is True,
            "Current post and first background packets are required")
    directions = bg["comparison_directions"]
    first_directions = first["comparison_directions"]
    require(directions == first_directions, "Current post directions must retain the first source")
    qproof = first["Q_b_proof"]
    delta = _read_record(c, qproof["delta"])
    h = _point(c, 1) - delta
    return dict(
        Z=_rows(c, directions["Z"], 3, "source Z direction"),
        d=_rows(c, directions["d"], 3, "source d direction"),
        Lz=_rows(c, directions["L"], 3, "source Lz direction"),
        delta=delta,
        h=h,
        range_Z=[_point(c, "37/100"), _point(c, "93/250")],
        source_formula=qproof.get("formula"),
    )


def _current_post_fields(c, parent):
    post = parent["actual_switch_post_background_I1_numeric"]["charts"]["post_switch"]
    require(post.get("passed") is True and
            post.get("source_proof", {}).get("source_executed") is False,
            "Current post source packet is required")
    native = post.get("native_fields")
    require(isinstance(native, Mapping), "Current post native fields are required")
    U = _rows(c, native["U"], 5, "current post U")
    M = _rows(c, native["M"], 5, "current post M")
    return dict(U=U, M=M, source_contract=post.get("source_proof", {}))


def _inherited_pressure_lemma(root, parent):
    """Validate the parent pressure identity without re-running its proof."""
    packet = parent.get("actual_switch_post_finite_N_I1_numeric")
    require(isinstance(packet, Mapping),
            "Accepted switch-post finite-N packet is required")
    lemma = packet.get("pressure_B_plus_D_identity")
    require(isinstance(lemma, Mapping) and lemma.get("passed") is True,
            "Accepted B+D pressure identity is required")
    require(lemma.get("symbolic_remainder") == "0" and
            lemma.get("Q4_reduced_pressure_source_guarded") is True and
            lemma.get("ordinary_Z_identity") is True and
            lemma.get("physical_y_numeric_bounds_claimed") is False,
            "Accepted B+D pressure identity contract changed")
    pins = packet.get("source_pins", {})
    expected = parent.get("input_hashes", {}).get(Path(Q4_SOURCE).name)
    q4_pin = pins.get(Q4_SOURCE)
    require(expected is not None and isinstance(q4_pin, Mapping) and
            q4_pin.get("expected_sha256") == expected and
            q4_pin.get("passed") is True and
            q4_pin.get("source_executed", False) is False,
            "Accepted Q4 source pin is required for inherited B+D identity")
    return dict(
        source_stage="actual_switch_post_finite_N_I1_numeric.pressure_B_plus_D_identity",
        passed=True,
        symbolic_remainder="0",
        Q4_reduced_pressure_source_guarded=True,
        ordinary_Z_identity=True,
        physical_y_numeric_bounds_claimed=False,
        Q4_source_pin_verified=True,
        reexecuted=False,
        identity=lemma.get("identity"),
        boundary_integrand=lemma.get("boundary_integrand"),
        mean_relation=lemma.get("mean_relation"),
    )


def _tail_cap(c):
    r = c.power(c.mpf(2), -TAIL_BITS)
    require(_lo(r) > 0, "Positive directed tail cap required")
    return r, c.log(r)


def _moments(c, L, r, log_r, name):
    require(_lo(L) > 0, name + " length must be positive")
    _assert_lt(-L, log_r, name + " tail must be below the source cap")
    log_L = c.log(L)
    log_L_minus_2L = log_L - c.mpf(2) * L
    _assert_lt(log_L_minus_2L, log_r,
               name + " L*exp(-2L) must be below the source cap")
    one = _point(c, 1)
    J0 = c.mpf([(c.mpf(1) - _hi(r)) / 2, c.mpf("0.5")])
    tail = c.mpf([c.mpf(0), _hi(r)])
    pointwise_decay = c.mpf([c.mpf(0), c.mpf(1)])
    J1 = tail
    J2 = tail
    return dict(
        J=[J0, J1, J2],
        length=L,
        tail_cap=r,
        pointwise_decay=pointwise_decay,
        source_proof=dict(
            exact_formulas=dict(
                J0="(1-exp(-2L))/2",
                J1="exp(-L)-exp(-2L)",
                J2="L*exp(-2L)",
            ),
            inequalities=dict(
                L_positive=True,
                minus_L_lt_log_r=True,
                log_L_minus_2L_lt_log_r=True,
            ),
            L=_record_value(L),
            log_L=record(log_L),
            log_L_minus_2L=record(log_L_minus_2L),
            log_r=record(log_r),
            pointwise_decay_interval=record(pointwise_decay),
            endpoint_decay_interval=record(tail),
            exp_minus_L_nonzero=True,
            exp_minus_2L_nonzero=True,
            tails_not_set_to_zero=True,
            no_midpoint=True,
            source_executed=False,
        ),
    )


def _record_value(value):
    return record(value)


def _mean_transport(c, U, M_in, r):
    decay = c.mpf([c.mpf(0), _hi(r)])
    decay_jet = [decay] + [_point(c, 0) for unused in ORDERS[1:]]
    difference = [M_in[index] - U[index] for index in ORDERS]
    endpoint = [U[order] + _prod(c, difference, decay_jet, order)
                for order in ORDERS]
    histories = [_hull_interval(c, M_in[order], U[order]) for order in ORDERS]
    return dict(
        endpoint=endpoint,
        histories=histories,
        decay=decay,
        source=dict(
            formula="M=U+(M_in-U)*exp(-s)",
            endpoint_decay_box="[0,r] with exp(-L)>0",
            positive_kernel_hull="hull(M_in,U)",
            incoming_memory_retained=True,
            reset=False,
            source_executed=False,
        ),
    )


def _background(c, alpha, beta, U, directions, moments, D, normalization):
    Z, Lz = directions["Z"], directions["Lz"]
    J0, J1, J2 = moments["J"]
    invL = _reciprocal(c, Lz, 2)
    invL2 = [_prod(c, invL, invL, order) for order in C1_ORDERS]
    Z2 = [_prod(c, Z, Z, order) for order in C1_ORDERS]
    U2 = [_prod(c, U, U, order) for order in C1_ORDERS]
    alpha_over_L = [_prod(c, alpha, invL, order) for order in C1_ORDERS]
    alpha_squared_half = [_scale(c, c.mpf("0.5"), _prod(c, alpha, alpha, order))
                          for order in C1_ORDERS]
    Z2U2 = [_prod(c, Z2, U2, order) for order in C1_ORDERS]
    z2_u2_over_l2 = [_prod(c, Z2U2, invL2, order) for order in C1_ORDERS]
    C0 = [_sum(c, (alpha_over_L[order], alpha_squared_half[order],
                   _scale(c, -2, z2_u2_over_l2[order])))
          for order in C1_ORDERS]
    beta_over_L = [_prod(c, beta, invL, order) for order in C1_ORDERS]
    alpha_beta = [_prod(c, alpha, beta, order) for order in C1_ORDERS]
    C1 = [_sum(c, (beta_over_L[order], alpha_beta[order]))
          for order in C1_ORDERS]
    C2 = [_scale(c, c.mpf("0.5"), _prod(c, beta, beta, order))
          for order in C1_ORDERS]
    B = [_sum(c, (_scale(c, 2, _scale(c, J0, C0[order])),
                  _scale(c, J1, C1[order])))
         for order in C1_ORDERS]
    S1 = [_sum(c, (B[order], D[order])) for order in C1_ORDERS]
    return dict(
        B=B, S1=S1, C0=C0, C1=C1, C2=C2,
        receipt=dict(
            formula="Q=alpha+beta*exp(-y); G=C0+C1*exp(-y)+C2*exp(-2y)",
            C0="alpha/Lz+alpha^2/2-2*Z^2*U^2/Lz^2",
            C1="beta/Lz+alpha*beta",
            C2="beta^2/2",
            normalized_formula="B/W=2*C0*J0+C1*J1",
            S1_formula="S1=B+D",
            C2_symbolically_cancelled=True,
            C2_cancelled_before_interval_boxing=True,
            C2_not_used_in_B_rows=True,
            no_endpoint_subtraction=True,
            normalization=normalization,
            source_executed=False,
        ),
    )


def _chart_arithmetic(c, U, M_in, directions, moments, normalization, name):
    Z, d, Lz, h = directions["Z"], directions["d"], directions["Lz"], directions["h"]
    delta = directions["delta"]
    invL = _reciprocal(c, Lz, 3)
    one_plus_delta = _point(c, 1) + delta
    U3 = U[:3]
    M3 = M_in[:3]
    U_Z = [U[1], U[2], U[3]]
    M_minus_U = [M3[index] - U3[index] for index in Q_ORDERS]
    M_Z_minus_U_Z = [M_in[index + 1] - U[index + 1] for index in Q_ORDERS]
    one_delta_Z = [one_plus_delta * Z[index] for index in Q_ORDERS]
    alpha_num = [_prod(c, one_delta_Z, U3, order)
                 - _prod(c, d, U_Z, order) for order in Q_ORDERS]
    alpha = [_prod(c, alpha_num, invL, order) for order in Q_ORDERS]
    hZ = [h * Z[index] for index in Q_ORDERS]
    beta_num = [_neg(c, _prod(c, hZ, M_minus_U, order))
                - _prod(c, d, M_Z_minus_U_Z, order)
                for order in Q_ORDERS]
    beta = [_prod(c, beta_num, invL, order) for order in Q_ORDERS]
    J0, J1, J2 = moments["J"]
    pointwise_decay = moments["pointwise_decay"]
    endpoint_decay = c.mpf([c.mpf(0), _hi(moments["tail_cap"])])
    Q = [alpha[order] + beta[order] * pointwise_decay
         for order in Q_ORDERS]
    Q_endpoint = [alpha[order] + beta[order] * endpoint_decay
                  for order in Q_ORDERS]
    integrated_q = [_scale(c, J0, alpha[order]) + _scale(c, J1, beta[order])
                    for order in Q_ORDERS]
    q_Z = [integrated_q[1], integrated_q[2]]
    q_squared = [(_scale(c, J0, _prod(c, alpha, alpha, order))
                  + _scale(c, c.mpf(2) * J1, _prod(c, alpha, beta, order))
                  + _scale(c, J2, _prod(c, beta, beta, order)))
                 for order in C1_ORDERS]
    U_q_Z = [_prod(c, U[:2], q_Z, order) for order in C1_ORDERS]
    U_squared = [_scale(c, J0, _prod(c, U[:2], U[:2], order))
                 for order in C1_ORDERS]
    U_U_Z = [_scale(c, J0, _prod(c, U[:2], U_Z[:2], order))
             for order in C1_ORDERS]
    invL2 = [_prod(c, invL, invL, order) for order in C1_ORDERS]
    hZ1 = [h * Z[index] for index in C1_ORDERS]
    Z2 = [_prod(c, Z, Z, order) for order in C1_ORDERS]
    h_Z2 = [h * Z2[index] for index in C1_ORDERS]
    Zd = [_prod(c, Z, d, order) for order in C1_ORDERS]
    coeff = dict(
        minus_inv_L=[_neg(c, value) for value in invL[:2]],
        hZ_over_2L=[_scale(c, c.mpf("0.5"), _prod(c, hZ1, invL, order))
                    for order in C1_ORDERS],
        minus_half=[_point(c, "-0.5"), _point(c, 0)],
        d_over_L=[_prod(c, d, invL, order) for order in C1_ORDERS],
        two_h_Z2_over_L2=[_scale(c, 2, _prod(c, h_Z2, invL2, order))
                          for order in C1_ORDERS],
        two_Zd_over_L2=[_scale(c, 2, _prod(c, Zd, invL2, order))
                        for order in C1_ORDERS],
    )
    integrals = dict(
        q=integrated_q[:2], q_Z=q_Z, q_squared=q_squared,
        U_q_Z=U_q_Z, U_squared=U_squared, U_U_Z=U_U_Z,
    )
    D = []
    for order in C1_ORDERS:
        D.append(_sum(c, (
            _prod(c, coeff["minus_inv_L"], integrals["q"], order),
            _prod(c, coeff["hZ_over_2L"], integrals["q_Z"], order),
            _prod(c, coeff["minus_half"], integrals["q_squared"], order),
            _prod(c, coeff["d_over_L"], integrals["U_q_Z"], order),
            _prod(c, coeff["two_h_Z2_over_L2"], integrals["U_squared"], order),
            _prod(c, coeff["two_Zd_over_L2"], integrals["U_U_Z"], order),
        )))
    background = _background(c, alpha, beta, U, directions, moments, D, normalization)
    return dict(
        Q=Q, Q_endpoint=Q_endpoint, integrated_q=integrated_q,
        alpha=alpha, beta=beta,
        integrals=integrals, D=D, coefficients=coeff,
        background=background,
        source=dict(
            chart=name,
            Q_pointwise_formula="Q=alpha+beta*[0,1] for the whole segment",
            Q_endpoint_formula="Q(L)=alpha+beta*[0,r] at the endpoint",
            alpha_formula="((1+delta)*Z*U-d*U_Z)/Lz",
            beta_formula="(-h*Z*(M_in-U)-d*(M_in_Z-U_Z))/Lz",
            integrated_q_formula="I[q]=alpha*J0+beta*J1",
            D_formula=("D=-I[q]/Lz+h*Z*I[q_Z]/(2Lz)-I[q_squared]/2+"
                       "d*I[U_q_Z]/Lz+2*h*Z^2*I[U_squared]/Lz^2+"
                       "2*Z*d*I[U_U_Z]/Lz^2"),
            U_squared_coefficient="2*h*Z^2/Lz^2",
            canonical_coefficient_key="two_h_Z2_over_L2",
            Q_and_integrated_Iq_distinct=True,
            pointwise_decay_interval="[0,1]",
            endpoint_decay_interval="[0,r]",
            decay_Z_independent=True,
            physical_Q4_source_guarded_only=True,
            source_executed=False,
        ),
    )


def _join_rows(c, long_rows, reference_rows, ratio):
    require(len(long_rows) == len(reference_rows), "Joined rows have different lengths")
    return [ratio * long_rows[index] + reference_rows[index]
            for index in range(len(long_rows))]


def _joined(c, long_packet, reference_packet, ratio, log_ratio):
    six = {}
    for name in I1_NAMES:
        six[name] = _pack_rows(_join_rows(
            c, _rows(c, long_packet["arithmetic"]["six_integrals"][name], 2,
                     "long " + name),
            _rows(c, reference_packet["arithmetic"]["six_integrals"][name], 2,
                  "reference " + name), ratio))
    long_B = _rows(c, long_packet["arithmetic"]["B_over_W"], 2, "long B")
    reference_B = _rows(c, reference_packet["arithmetic"]["B_over_W"], 2,
                        "reference B")
    long_D = _rows(c, long_packet["arithmetic"]["D"], 2, "long D")
    reference_D = _rows(c, reference_packet["arithmetic"]["D"], 2,
                        "reference D")
    long_S1 = _rows(c, long_packet["arithmetic"]["S1_over_W"], 2, "long S1")
    reference_S1 = _rows(c, reference_packet["arithmetic"]["S1_over_W"], 2,
                         "reference S1")
    return dict(
        normalization="Wref=Rz^2",
        Wlong_over_Wref=record(ratio),
        log_Wlong_over_Wref=record(log_ratio),
        six_integrals=six,
        B_over_Wref=_pack_rows(_join_rows(c, long_B, reference_B, ratio)),
        D_over_Wref=_pack_rows(_join_rows(c, long_D, reference_D, ratio)),
        S1_over_Wref=_pack_rows(_join_rows(c, long_S1, reference_S1, ratio)),
        source_proof=dict(
            exact_ratio="Wlong/Wref=exp(-2*Lref)",
            ratio_positive_nonzero=True,
            ratio_upper_is_tail_cap=True,
            no_selected_zero=True,
            no_independent_radius_added=True,
            source_executed=False,
        ),
    )


def _chart_packet(c, name, normalization, U, M_in, transport, moments,
                  arithmetic, log_R_upper, log_W):
    return dict(
        passed=True,
        chart=name,
        normalization=normalization,
        length=record(moments["length"]),
        log_R_upper=record(log_R_upper),
        log_W=record(log_W),
        W_formula="W=R_upper^2; radius not materialized",
        U_constant_rows=_pack_rows(U, kind="native_constant_U"),
        incoming_M=_pack_rows(M_in, kind="incoming_current_post_M" if name == "long_reshape"
                              else "incoming_long_endpoint_M"),
        endpoint_M=_pack_rows(transport["endpoint"], kind="tight_endpoint_M"),
        whole_segment_M_hull=_pack_rows(transport["histories"], kind="positive_kernel_M_hull"),
        mean_endpoint_decay_interval=record(transport["decay"]),
        normalized_moments=[record(value) for value in moments["J"]],
        arithmetic=dict(
            Q_pointwise=_pack_rows(arithmetic["Q"], kind="pointwise_Q"),
            Q_endpoint=_pack_rows(arithmetic["Q_endpoint"], kind="endpoint_Q"),
            integrated_q=_pack_rows(arithmetic["integrated_q"], kind="integrated_Iq"),
            alpha=_pack_rows(arithmetic["alpha"], kind="alpha_C2"),
            beta=_pack_rows(arithmetic["beta"], kind="beta_C2"),
            six_integrals={name: _pack_rows(rows, kind="ordinary_I1")
                           for name, rows in arithmetic["integrals"].items()},
            D=_pack_rows(arithmetic["D"], kind="ordinary_D_over_W"),
            B_over_W=_pack_rows(arithmetic["background"]["B"], kind="ordinary_B_over_W"),
            S1_over_W=_pack_rows(arithmetic["background"]["S1"], kind="ordinary_S1_over_W"),
            coefficients={name: _pack_rows(rows, kind="ordinary_C1_coefficient")
                          for name, rows in arithmetic["coefficients"].items()},
            C0=_pack_rows(arithmetic["background"]["C0"], kind="background_C0_C1"),
            C1=_pack_rows(arithmetic["background"]["C1"], kind="background_C1_C1"),
            C2_cancelled=_pack_rows(arithmetic["background"]["C2"], kind="background_C2_cancelled"),
            algebra=arithmetic["background"]["receipt"],
            source=arithmetic["source"],
        ),
        tail_proof=moments["source_proof"],
        source_executed=False,
    )


def build(root, parent):
    """Build native constant-``U`` long/reference ``Q`` and pressure rows."""
    root = Path(root)
    disk_parent = _load_parent(root, parent)
    pins = _source_pins(root, disk_parent)
    inherited_lemma = _inherited_pressure_lemma(root, disk_parent)
    guards = _source_guards(root)
    c6 = _load_c6_cache(root)
    graph = _scalar_graph_guard(c6, disk_parent)
    c = MPIntervalContext()
    c.prec = 1600

    scalar_source = disk_parent["actual_micro_mean_Z3_bounds"]["physical_scale_logs"]
    logC = _read_record(c, scalar_source["selected_logC"])
    node_T = c6["prescribed_R100_Rm_C6_bindings"]["geometry_certificate"]["actual_T_singleton_tuple"]
    require(node_T["lower_exact_mpf_tuple"] == node_T["upper_exact_mpf_tuple"],
            "Actual T source must be a singleton")
    T = _read_record(c, node_T, directed=False)
    logP = c.exp(c.mpf(40)) + c.mpf(11)
    logref = c.mpf(10) * (logC + logP)
    Lref = logref - T - c.mpf(8)
    r, log_r = _tail_cap(c)
    _assert_lt(-T, log_r, "Source proof exp(-T)<r failed")
    _assert_lt(-Lref, log_r, "Source proof exp(-Lref)<r failed")
    ratio = c.mpf([c.mpf(0), _hi(r)])
    log_ratio = -c.mpf(2) * Lref
    _assert_lt(log_ratio, log_r, "Source proof Wlong/Wref<r failed")

    directions = _directions(c, disk_parent)
    current = _current_post_fields(c, disk_parent)
    U, M_current = current["U"], current["M"]
    long_transport = _mean_transport(c, U, M_current, r)
    reference_transport = _mean_transport(c, U, long_transport["endpoint"], r)
    long_moments = _moments(c, T, r, log_r, "long_reshape")
    reference_moments = _moments(c, Lref, r, log_r, "reference")
    long_arithmetic = _chart_arithmetic(
        c, U, M_current, directions, long_moments, "Wlong", "long_reshape")
    reference_arithmetic = _chart_arithmetic(
        c, U, long_transport["endpoint"], directions,
        reference_moments, "Wref", "reference")
    log110 = c.log(c.mpf(110))
    long_log_R = log110 + T
    reference_log_R = log110 + logref - c.mpf(8)
    long_log_W = c.mpf(2) * long_log_R
    reference_log_W = c.mpf(2) * reference_log_R
    long_packet = _chart_packet(
        c, "long_reshape", "Wlong", U, M_current, long_transport,
        long_moments, long_arithmetic, long_log_R, long_log_W)
    reference_packet = _chart_packet(
        c, "reference", "Wref", U, long_transport["endpoint"],
        reference_transport, reference_moments, reference_arithmetic,
        reference_log_R, reference_log_W)
    long_packet["reference_incoming_M"] = _pack_rows(
        long_transport["endpoint"], kind="reference_incoming_long_endpoint_M")
    joined = _joined(c, long_packet, reference_packet, ratio, log_ratio)

    result = dict(
        passed=True,
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        parentSHA=PARENT_SHA256,
        ordinary_Z_orders=list(ORDERS),
        I1_ordinary_Z_orders=list(C1_ORDERS),
        Q_pointwise_ordinary_Z_orders=list(Q_ORDERS),
        source_pins=pins,
        source_guards=guards,
        sourceguards_tuple=True,
        source_executed=False,
        physical_y_derivatives_claimed=False,
        source_family=disk_parent.get("source_family"),
        scalar_graph=dict(
            **graph,
            logP=record(logP),
            logC=record(logC),
            T=record(T),
            logref=record(logref),
            reference_length=record(Lref),
        ),
        source_proof=dict(
            raw_actual_T_singleton_tuple=[
                list(node_T["lower_exact_mpf_tuple"]),
                list(node_T["upper_exact_mpf_tuple"]),
            ],
            current_post_U_M_source="parent.actual_switch_post_background_I1_numeric.charts.post_switch.native_fields",
            C6_cache_scalar_source="prescribed_R100_Rm_C6_functions.json.gz",
            source_executed=False,
        ),
        inherited_lemma_reference=inherited_lemma,
        radii=dict(
            long_log_R_upper=record(long_log_R),
            reference_log_R_upper=record(reference_log_R),
            long_log_W=record(long_log_W),
            reference_log_W=record(reference_log_W),
            long_W_formula="Wlong=110^2*exp(2T), represented by logWlong",
            reference_W_formula="Wref=110^2*exp(2*(logref-8)), represented by logWref",
            radius_not_materialized=True,
            exp_huge_positive_not_materialized=True,
        ),
        tail_control=dict(
            r=record(r),
            log_r=record(log_r),
            T_tail_bound="exp(-T)<r",
            reference_tail_bound="exp(-Lref)<r",
            W_ratio_bound="Wlong/Wref=exp(-2*Lref)<r",
            tail_values_are_nonzero=True,
            lower_zero_is_an_interval_enclosure_not_a_selected_value=True,
            no_midpoint=True,
        ),
        directions=dict(
            Z=_pack_rows(directions["Z"], kind="source_Z_direction"),
            d=_pack_rows(directions["d"], kind="source_d_direction"),
            Lz=_pack_rows(directions["Lz"], kind="source_Lz_direction"),
            delta=record(directions["delta"]),
            h=record(directions["h"]),
            source_formula=directions["source_formula"],
            Z_domain="[37/100,93/250]",
            Z_domain_bounds=[record(value) for value in directions["range_Z"]],
        ),
        charts={"long_reshape": long_packet, "reference": reference_packet},
        joined_long_reference_background=joined,
        source_contract=dict(
            passed=True,
            current_post_U_and_M_rows_consumed=True,
            current_post_ordinary_rows=list(ORDERS),
            c6_long_reference_U_constant=True,
            c6_histories_retained_source_guard=True,
            physical_Q4_source_guarded=True,
            long_endpoint_M_formula="M_Rsh=U+(M_R110-U)*exp(-T)",
            reference_incoming_is_long_endpoint=True,
            reference_endpoint_M_formula="M_Rz=U+(M_Rsh-U)*exp(-Lref)",
            M_not_reset_to_U=True,
            whole_segment_pointwise_Q_decay="[0,1]",
            endpoint_Q_decay="[0,r]",
            pointwise_Q_distinct_from_integrated_Iq=True,
            normalized_by_Wlong_and_Wref=True,
            common_Wref_join_included=True,
            C2_cancelled_before_interval_boxing=True,
            U_squared_coefficient="2*h*Z^2/Lz^2",
            inherited_S1_lemma_reused=True,
            no_radius_materialization=True,
            source_executed=False,
        ),
        limits=dict(
            native_phi_not_evaluated=True,
            native_A_B_C_histories_not_evaluated=True,
            finite_N_loop_not_evaluated=True,
            current_N_not_enclosed=True,
            full_pressure_global_closure=False,
            physical_y_derivatives_claimed=False,
            accepted_producer_checker_owner_norm_recurrence_quadrature_not_called=True,
            source_executed=False,
        ),
    )
    json.dumps(result)
    return result


def fixture_bound_check():
    """Run only arithmetic fixtures owned by this module."""
    c = MPIntervalContext()
    c.prec = 180
    one = _point(c, 1)
    product = _prod(c, [one] * 5, [one] * 5, 4)
    require(_lo(product) == 16 and _hi(product) == 16,
            "Fresh ordinary binomial product fixture failed")
    r = c.mpf([c.mpf(0), c.mpf("0.1")])
    U = [_point(c, 2), _point(c, 0), _point(c, 0), _point(c, 0), _point(c, 0)]
    M = [_point(c, 4), _point(c, 0), _point(c, 0), _point(c, 0), _point(c, 0)]
    transport = _mean_transport(c, U, M, r)
    require(_lo(transport["endpoint"][0]) >= 2 and
            _hi(transport["endpoint"][0]) < 2.3,
            "Fresh mean endpoint fixture failed")
    directions = dict(
        Z=[one, _point(c, 0), _point(c, 0)],
        d=[one, _point(c, 0), _point(c, 0)],
        Lz=[one, _point(c, 0), _point(c, 0)],
        delta=_point(c, 0), h=one,
    )
    moments = dict(J=[one, one, one], tail_cap=r,
                   pointwise_decay=c.mpf([c.mpf(0), c.mpf(1)]), length=one,
                   source_proof=dict())
    alpha = [one, _point(c, 0), _point(c, 0)]
    beta = [c.mpf(2), _point(c, 0), _point(c, 0)]
    background = _background(c, alpha, beta,
                             [_point(c, 1), _point(c, 0), _point(c, 0),
                              _point(c, 0), _point(c, 0)],
                             directions, moments,
                             [_point(c, 0), _point(c, 0)], "Wlong")
    require(_lo(background["B"][0]) == 3 and _hi(background["B"][0]) == 3,
            "Fresh C2 cancellation fixture failed")
    require(background["receipt"]["C2_symbolically_cancelled"] is True,
            "Fresh C2 cancellation receipt failed")
    q_moments = dict(
        J=[c.mpf([c.mpf("0.5"), c.mpf("0.5")]),
           c.mpf([c.mpf("0.1"), c.mpf("0.1")]),
           c.mpf([c.mpf("0.1"), c.mpf("0.1")])],
        tail_cap=c.mpf("0.1"),
        pointwise_decay=c.mpf([c.mpf(0), c.mpf(1)]),
        length=one,
        source_proof=dict(),
    )
    q_U = [_point(c, 1), _point(c, 0), _point(c, 0),
           _point(c, 0), _point(c, 0)]
    q_M = [_point(c, 1), _point(c, -2), _point(c, 0),
           _point(c, 0), _point(c, 0)]
    q_arithmetic = _chart_arithmetic(
        c, q_U, q_M, directions, q_moments, "Wlong", "fixture")
    require(_lo(q_arithmetic["Q"][0]) <= 1 and
            _hi(q_arithmetic["Q"][0]) >= 3,
            "Whole-segment pointwise Q fixture failed")
    require(_lo(q_arithmetic["Q_endpoint"][0]) <= 1 and
            _hi(q_arithmetic["Q_endpoint"][0]) >= 1,
            "Endpoint Q fixture failed")
    require(_hi(q_arithmetic["Q_endpoint"][0]) < c.mpf("1.3"),
            "Endpoint Q did not retain [0,r] decay")
    require(_hi(q_arithmetic["integrated_q"][0]) < 1,
            "Integrated Q fixture did not remain distinct")
    return dict(
        passed=True,
        ordinary_binomial_order4=16,
        mean_endpoint_tail_box=True,
        C2_cancelled_before_interval_boxing=True,
        whole_segment_Q_uses_0_1=True,
        endpoint_Q_uses_0_r=True,
        integrated_Q_distinct=True,
        source_executed=False,
    )


def self_check():
    return fixture_bound_check()
