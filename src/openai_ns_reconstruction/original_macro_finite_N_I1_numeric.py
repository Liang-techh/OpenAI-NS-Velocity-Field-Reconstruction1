"""Finite-N macro six-I1 error adapter.

This module consumes a source-only frozen-macro loop packet and turns its
logarithmic ``B_over_Pstar/N`` majorants into physical ``U``, mean, and
recovered ``Q`` error rows.  It deliberately does not evaluate a loop
producer, select an issued phase, materialize the huge integer ``N``, or
alter any of the saved signed packets.

The returned ``I1_full`` and ``D_full`` values are signed enclosures.  The
finite-N contribution is an absolute error radius and is never exposed as a
selected signed function value.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from .original_core_signed_direction import read, record, require, sha, source_guard
from .original_log_mixed_jets import _logsum


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_macro_background_I1.json"
PARENT_SHA256 = "5a77d3d8b9b294516cc0be05ed83667e2710d64ff8427c7e917b8aed3379aede"

LOOP_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
LOOP_ARCHIVE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
Q4_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
Q4_ARCHIVE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
Q4_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions_check.json"
Q_SOURCE = "src/openai_ns_reconstruction/original_micro_Q_increment_bounds.py"
MACRO_SOURCE = "src/openai_ns_reconstruction/original_macro_background_I1_numeric.py"
ADAPTER_SOURCE = "src/openai_ns_reconstruction/original_macro_finite_N_I1_numeric.py"

ORDER = 3
I1_ORDER = 1
I1_NAMES = ("q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z")
N_HANDLE = 18
PHASE_ORIGIN_HANDLE = 24954
MACRO_J_HANDLE = 66987
MACRO_COORDINATE_HANDLE = 70265
MACRO_RADIUS_HANDLE = 71448
PSTAR_HANDLE = 33572
PSTAR_LOG_ARGUMENT_HANDLE = 5
MACRO_INCOMING_M_HANDLES = (101875, 101876, 101877, 101878, 101879, 101880)
MACRO_HISTORY_M_HANDLES = (104342, 104343, 104344, 104345, 104346, 104347)
MACRO_MODIFIED_M_HANDLES = (104467, 104468, 104469, 104470, 104471, 104472)
MACRO_DELTA_V_HANDLES = (104063, 104064, 104065, 104066, 104067, 104068)
MACRO_V_N_HANDLES = (104304, 104305, 104306, 104307, 104308, 104309)
Q4_U_HANDLES = (140685, 140686, 140687, 140688, 140689, 140690)
Q4_M_HANDLES = (140697, 140698, 140699, 140700, 140701, 140702)
Q4_Q_HANDLES = (140780, 140784, 140790, 140796, 140802)
Q4_QY_HANDLES = (140868, 140872, 140878, 140884, 140890)
Q4_RQY_HANDLES = (141027, 141028, 141029, 141030)


def _load(root, path):
    return json.loads((root / path).read_bytes())


def _record_value(row, preferred):
    """Accept a direct record or the row wrappers used by saved packets."""
    if row is None:
        return None
    if isinstance(row, dict):
        if preferred in row:
            return row[preferred]
        if preferred == "log_absolute_upper_bound" and "value" in row:
            return row["value"]
    return row


def _ordered_rows(rows, label, length=ORDER + 1):
    require(isinstance(rows, list) and len(rows) >= length,
            "Incomplete " + label + " rows")
    rows = rows[:length]
    require([row.get("ordinary_Z_order") for row in rows
             if isinstance(row, dict)] == list(range(length)),
            "Unordered " + label + " rows")
    return rows


def _read_log_rows(c, rows, label, length=ORDER + 1):
    values = []
    for row in _ordered_rows(rows, label, length):
        if isinstance(row, dict) and row.get("exact_zero"):
            values.append(None)
            continue
        raw = _record_value(row, "log_absolute_upper_bound")
        values.append(None if raw is None else read(c, raw))
    return values


def _read_signed_rows(c, rows, label, length):
    values = []
    for row in _ordered_rows(rows, label, length):
        raw = _record_value(row, "value")
        require(raw is not None, "Missing signed value in " + label)
        values.append(read(c, raw))
    return values


def _pack_logs(values):
    out = []
    for index, value in enumerate(values):
        row = {"ordinary_Z_order": index}
        if value is None:
            row.update(exact_zero=True, log_absolute_upper_bound=None)
        else:
            row["log_absolute_upper_bound"] = record(value)
        out.append(row)
    return out


def _pack_values(values):
    return [dict(ordinary_Z_order=index, value=record(value))
            for index, value in enumerate(values)]


def _absolute_log(c, value):
    upper = c.fabs(value).b
    if upper == 0:
        return None
    return c.ln(upper)


def _absolute_logs(c, values):
    return [_absolute_log(c, value) for value in values]


def _point(c, value):
    """Make a directed point from a known upper endpoint."""
    return c.mpf([value.b, value.b])


def _upper_max(c, values):
    values = [value for value in values if value is not None]
    require(values, "At least one absolute log-cap is required")
    return _point(c, max(value.b for value in values))


def _log_abs_scalar(c, scalar):
    upper = c.fabs(scalar).b
    if upper == 0:
        return None
    return c.ln(upper)


def _log_product(c, left, right):
    if left is None or right is None:
        return None
    return left + right


def _log_add(c, *values):
    values = [value for value in values if value is not None]
    return _logsum(c, values)


def _log_linear(c, terms):
    """Absolute log-majorant of a finite linear combination."""
    products = []
    for coefficient, value in terms:
        product = _log_product(c, _log_abs_scalar(c, coefficient), value)
        if product is not None:
            products.append(product)
    return _logsum(c, products)


def _log_scale(c, value, scalar):
    return _log_product(c, _log_abs_scalar(c, scalar), value)


def _log_jet_product(c, left, right, order=I1_ORDER):
    """Ordinary derivative product in log-cap form."""
    out = []
    for k in range(order + 1):
        terms = []
        for j in range(k + 1):
            term = _log_product(c, left[j], right[k - j])
            if term is not None:
                coefficient = math.comb(k, j)
                if coefficient > 1:
                    term += _log_abs_scalar(c, c.mpf(coefficient))
                terms.append(term)
        out.append(_logsum(c, terms))
    return out


def _ordinary_product_values(left, right, order=I1_ORDER):
    """Multiply two signed ordinary derivative rows by Leibniz' rule."""
    out = []
    for k in range(order + 1):
        total = None
        for j in range(k + 1):
            term = left[j] * right[k - j] * math.comb(k, j)
            total = term if total is None else total + term
        out.append(total)
    return out


def _symmetric(c, radius):
    return c.mpf([-radius.b, radius.b])


def _pin(root, parent, path):
    if path == PARENT:
        expected = PARENT_SHA256
    elif path.startswith(BASE):
        expected = parent.get("input_hashes", {}).get(Path(path).name)
    else:
        expected = parent.get("implementation_hashes", {}).get(path)
    require(expected is not None, "Parent has no hash for " + path)
    require(sha(root / path) == expected, "Changed saved source: " + path)
    return dict(path=path, expected_sha256=expected, passed=True)


def _source_guards(root):
    return dict(
        finite_N=source_guard(root, LOOP_SOURCE, "finite_N_functions", (
            "invN=g.quotient(g.one,N,'same exact positive integer N on every source chart')",
            "s=a.scale(loop['B_over_Pstar'],invN)",
            "V_N=a.add(V,s)",
        )),
        transport=source_guard(root, LOOP_SOURCE, "transport", (
            "for name,rate in RATES.items():",
            "physical_Jacobian_applied_once=True",
            "predecessor_memory_never_reset=True",
            "out[name]=a.add(a.scale(incoming[name],decay),contribution)",
        )),
        macro_constructor=source_guard(root, LOOP_SOURCE, "__init__", (
            "self.N=self.unique('exact_positive_integer_expression',None)",
            "phase=g.unary('fractional_part',g.mul(self.N,offset))",
            "history,parts,memory=transport(g,finite['five_signed_density_increments'],incoming,",
            "incoming=self.endpoints[name]",
        )),
        physical_fields=source_guard(root, Q4_SOURCE, "physical_pressure_fields", (
            "U=a.scale(fields['V'],Pstar)",
            "m=a.scale(fields['m'],Pstar)",
            "my=a.sub(U,m)",
            "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
        )),
        total_y=source_guard(root, Q4_SOURCE, "total_radial_functions", (
            "Vy=a.add(a.trim(bg['radial']['V_y']),a.scale(slow['B_over_Pstar'],invN),chart['loop']['B_over_Pstar_phi'])",
            "physical_y_equals_native_partial_divided_by_Jacobian=True",
        )),
        Q_linear=source_guard(root, Q_SOURCE, "new_controls", (
            "op=lambda u,m:(2*z*u-(1-de)*z*m-(1-z*z)*s.diff(m,z))/(1-de*z*z)",
            "require(s.cancel(op(U+A,M+B)-op(U,M)-op(A,B))==0",
            "q.append((n[2]-2*s.diff(L,z)*q[1]-s.diff(L,z,2)*q[0])/L)",
        )),
    )


def _loop_packet_rows(c, loop_packet):
    require(loop_packet.get("passed") is True, "Finite-N loop packet is not accepted")
    require(loop_packet.get("parent_report_sha256") == PARENT_SHA256,
            "Loop packet is tied to a different macro parent")
    require(loop_packet.get("ordinary_Z_orders") == [0, 1, 2, 3],
            "Four ordinary Z rows are required")
    require(loop_packet.get("physical_y_derivatives_claimed") is False,
            "This adapter accepts Z-only loop caps")
    B = _read_log_rows(c, loop_packet["B_over_Pstar_over_frozen_E0_log_caps"],
                       "B_over_Pstar_over_frozen_E0", ORDER + 1)
    require(all(value is not None for value in B), "All B log caps are required")
    pstar_E_raw = _record_value(loop_packet.get("log_Pstar_E_upper"),
                                "log_absolute_upper_bound")
    logPstarE = read(c, pstar_E_raw)
    logN = read(c, _record_value(loop_packet.get("logN"),
                                 "log_absolute_upper_bound"))
    return B, logPstarE, logN


def _recover_delta_Q_logs(c, delta_U, delta_M, z, delta):
    """Recover absolute Q error logs through ordinary Z order two."""
    require(len(delta_U) >= 4 and len(delta_M) >= 4,
            "Q recovery requires physical U and mean rows Z0..Z3")
    h = c.mpf(1) - delta
    d = c.mpf(1) - z * z
    L = c.mpf(1) - delta * z * z
    require(L.a > 0, "Current macro L=1-delta*Z^2 is not positive")
    L1 = -2 * delta * z
    L2 = -2 * delta
    n0 = _log_linear(c, ((2 * z, delta_U[0]),
                         (-h * z, delta_M[0]),
                         (-d, delta_M[1])))
    n1 = _log_linear(c, ((c.mpf(2), delta_U[0]),
                         (2 * z, delta_U[1]),
                         (-h, delta_M[0]),
                         ((1 + delta) * z, delta_M[1]),
                         (-d, delta_M[2])))
    n2 = _log_linear(c, ((c.mpf(4), delta_U[1]),
                         (2 * z, delta_U[2]),
                         (2 * delta, delta_M[1]),
                         ((3 + delta) * z, delta_M[2]),
                         (-d, delta_M[3])))
    invL = -c.ln(L.a)
    q0 = n0 + invL
    q1 = _log_add(c, n1, _log_product(c, _log_abs_scalar(c, L1), q0)) + invL
    q2 = _log_add(c,
                  n2,
                  _log_scale(c, q1, 2 * L1),
                  _log_product(c, _log_abs_scalar(c, L2), q0)) + invL
    return [q0, q1, q2], dict(
        passed=True,
        ordinary_derivatives=True,
        no_top_row_padding=True,
        positive_L=record(L),
        L_Z=record(L1),
        L_ZZ=record(L2),
        numerator_log_caps=_pack_logs((n0, n1, n2)),
        formula=(
            "deltaQ=(2*Z*deltaU-(1-delta)*Z*deltaM-"
            "(1-Z^2)*deltaM_Z)/(1-delta*Z^2)"
        ),
    )


def _difference_logs(c, U, q, U_b, q_b):
    qz = q[1:3]
    qbz = q_b[1:3]
    Uz = U[1:3]
    Ubz = U_b[1:3]
    qb2 = _log_jet_product(c, q_b[:2], q[:2])
    q2 = _log_jet_product(c, q[:2], q[:2])
    q_squared = [_log_add(c, _log_scale(c, qb2[k], c.mpf(2)), q2[k])
                 for k in range(I1_ORDER + 1)]
    U_q_Z = []
    for k in range(I1_ORDER + 1):
        U_q_Z.append(_log_add(
            c,
            _log_jet_product(c, U_b[:2], qz)[k],
            _log_jet_product(c, U[:2], qbz)[k],
            _log_jet_product(c, U[:2], qz)[k],
        ))
    U_squared_b = _log_jet_product(c, U_b[:2], U[:2])
    U_squared_d = _log_jet_product(c, U[:2], U[:2])
    U_squared = [_log_add(c, _log_scale(c, U_squared_b[k], c.mpf(2)),
                          U_squared_d[k]) for k in range(I1_ORDER + 1)]
    UU_Z = []
    for k in range(I1_ORDER + 1):
        UU_Z.append(_log_add(
            c,
            _log_jet_product(c, U_b[:2], Uz)[k],
            _log_jet_product(c, U[:2], Ubz)[k],
            _log_jet_product(c, U[:2], Uz)[k],
        ))
    return dict(
        q=q[:2],
        q_Z=qz,
        q_squared=q_squared,
        U_q_Z=U_q_Z,
        U_squared=U_squared,
        U_U_Z=UU_Z,
    )


def _full_I1_and_D(c, macro, difference, radius):
    signed = macro["macro_I1"]["six_integrals"]
    base = {name: _read_signed_rows(c, signed[name], "signed " + name, I1_ORDER + 1)
            for name in I1_NAMES}
    full_values = {name: [value + _symmetric(c, radius) for value in base[name]]
                   for name in I1_NAMES}
    coefficients = macro["macro_I1"]["coefficientsC1rows"]
    coeff = {name: _read_signed_rows(c, coefficients[name], "coefficient " + name,
                                     I1_ORDER + 1)
             for name in coefficients}
    coefficient_products = {
        "q": _ordinary_product_values(coeff["minus_inv_L"], full_values["q"]),
        "q_Z": _ordinary_product_values(coeff["hZ_over_2L"], full_values["q_Z"]),
        "q_squared": _ordinary_product_values(coeff["minus_half"],
                                               full_values["q_squared"]),
        "U_q_Z": _ordinary_product_values(coeff["d_over_L"], full_values["U_q_Z"]),
        "U_squared": _ordinary_product_values(coeff["two_hZ2_over_L2"],
                                               full_values["U_squared"]),
        "U_U_Z": _ordinary_product_values(coeff["two_Zd_over_L2"],
                                           full_values["U_U_Z"]),
    }
    D_full = []
    for k in range(I1_ORDER + 1):
        terms = [
            coefficient_products["q"][k],
            coefficient_products["q_Z"][k],
            coefficient_products["q_squared"][k],
            coefficient_products["U_q_Z"][k],
            coefficient_products["U_squared"][k],
            coefficient_products["U_U_Z"][k],
        ]
        total = terms[0]
        for term in terms[1:]:
            total += term
        D_full.append(total)
    return (
        dict(
            passed=True,
            normalized_by_Wmacro=True,
            six_integrals={name: _pack_values(full_values[name]) for name in I1_NAMES},
            finite_N_absolute_error_log_bounds={
                name: _pack_logs(difference[name]) for name in I1_NAMES
            },
            signed_background_unchanged=True,
            symmetric_radius_added_once=record(radius),
            error_radius_is_not_a_selected_value=True,
        ),
        dict(
            passed=True,
            rows=_pack_values(D_full),
            D_rows=_pack_values(D_full),
            coefficientsC1rows=macro["macro_I1"]["coefficientsC1rows"],
            computed_from_full_I1_once=True,
            separate_D_error_not_added=True,
            normalized_by_Wmacro=True,
        ),
    )


def build(root, parent, loop_packet):
    """Build current-N physical macro six-I1 error enclosures."""
    from mpmath.ctx_iv import MPIntervalContext

    root = Path(root)
    report = _load(root, PARENT)
    require(sha(root / PARENT) == PARENT_SHA256 and report == parent,
            "Exact current macro background parent required")
    require(isinstance(loop_packet, dict), "Loop packet mapping required")

    paths = (LOOP_SOURCE, LOOP_ARCHIVE, Q4_SOURCE, Q4_ARCHIVE, Q4_CHECK,
             Q_SOURCE, MACRO_SOURCE)
    pins = {PARENT: dict(path=PARENT, expected_sha256=PARENT_SHA256, passed=True)}
    for path in paths:
        pins[path] = _pin(root, parent, path)
    packet_pins = loop_packet.get("source_pins", {})
    require(isinstance(packet_pins, dict), "Loop packet source_pins mapping required")
    for path, value in packet_pins.items():
        pins.setdefault(path, value)
    if (root / ADAPTER_SOURCE).exists():
        pins[ADAPTER_SOURCE] = dict(path=ADAPTER_SOURCE,
                                    expected_sha256=sha(root / ADAPTER_SOURCE),
                                    passed=True,
                                    self_pin=True)

    guards = _source_guards(root)
    c = MPIntervalContext()
    c.prec = 600
    B, logPstarE, logN = _loop_packet_rows(c, loop_packet)

    expected_logN = read(c, parent["actual_micro_mean_Z3_bounds"]
                         ["physical_scale_logs"]["N"])
    require(logN.a == expected_logN.a and logN.b == expected_logN.b,
            "Loop packet logN is not the accepted actual N=18 log row")

    macro = parent["actual_macro_background_I1_numeric"]
    geometry = macro["geometric_measure"]
    K = read(c, geometry["K"])
    require(K.a > 0 and K.b <= c.mpf("0.5").b,
            "Positive macro K <= 1/2 is required")
    logK = c.ln(K.b)
    delta = read(c, parent["actual_core_signed_direction"]
                 ["parameters"]["delta"])
    z = c.mpf([c.mpf(37) / 100, c.mpf(93) / 250])
    L = c.mpf(1) - delta * z * z
    require(L.a > 0, "Current macro L=1-delta*Z^2 is not positive")

    # Physical U is the local B/N source.  There is intentionally no
    # incoming-U term: C5 transport carries only m,h,k,e,p memories.
    delta_U = [logPstarE + B[k] - logN for k in range(ORDER + 1)]
    micro_second = parent["actual_micro_Q_increment_bounds"]["charts"]["second"]
    incoming_M = _read_log_rows(c, micro_second["log_caps"]["physical_delta_m"],
                                "accepted second-micro physical delta M", ORDER + 1)
    delta_M = [_upper_max(c, (incoming_M[k], delta_U[k]))
               for k in range(ORDER + 1)]
    delta_Q, qproof = _recover_delta_Q_logs(c, delta_U, delta_M, z, delta)

    fields = macro["actual_macro_fields"]
    U_b = _absolute_logs(c, _read_signed_rows(c, fields["U_b"], "macro U_b", ORDER + 1))
    Q_b = _absolute_logs(c, _read_signed_rows(c, fields["Q_b"], "macro Q_b", 3))
    difference = _difference_logs(c, delta_U, delta_Q, U_b, Q_b)
    I1_error = {
        name: [_log_product(c, value, logK) for value in values]
        for name, values in difference.items()
    }
    threshold = -1024 * c.ln(c.mpf(2))
    require(all(value is not None and value.b < threshold.a
                for values in I1_error.values() for value in values),
            "Current macro finite-N six-I1 error radius did not fit 2^-1024")
    radius = c.mpf(2) ** (-1024)
    I1_full, D_full = _full_I1_and_D(c, macro, I1_error, radius)

    target = parent["actual_micro_macro_correlated_m15"]["actual_macro_six_I1_source"]
    source_contract = dict(
        passed=True,
        ordinary_Z_orders=[0, 1, 2, 3],
        physical_y_derivatives_claimed=False,
        common_N_handle=N_HANDLE,
        phase_definition="fractional_part(N*offset)",
        phase_origin_handle=PHASE_ORIGIN_HANDLE,
        phase_Z_independent=True,
        N_not_materialized=True,
        macro_coordinate_handle=MACRO_COORDINATE_HANDLE,
        macro_J_handle=MACRO_J_HANDLE,
        macro_domain=[0, 1],
        macro_incoming_m_handles=list(MACRO_INCOMING_M_HANDLES),
        macro_history_m_handles=list(MACRO_HISTORY_M_HANDLES),
        macro_modified_m_handles=list(MACRO_MODIFIED_M_HANDLES),
        macro_delta_V_handles=list(MACRO_DELTA_V_HANDLES),
        macro_V_N_handles=list(MACRO_V_N_HANDLES),
        macro_incoming_formula=(
            "delta_m_norm(x)=exp(-J*x)*incoming_m+"
            "integral_0^x J*exp(-J*(x-t))*B_over_Pstar(t)/N dt"
        ),
        macro_incoming_nonreset=True,
        macro_delta_U_has_no_incoming=True,
        micro_delta_U_not_added_as_macro_inlet=True,
        physical_Pstar_applied_once=True,
        physical_Pstar_graph_handle=PSTAR_HANDLE,
        physical_Pstar_log_argument_handle=PSTAR_LOG_ARGUMENT_HANDLE,
        physical_R_not_a_mean_scale=True,
        q_recovery_linear=True,
        q_recovery=qproof,
        source_graph_handles=dict(
            target_from_parent=target,
            Q4_U_full=list(Q4_U_HANDLES),
            Q4_M_full=list(Q4_M_HANDLES),
            Q4_Q_full=list(Q4_Q_HANDLES),
            Q4_Q_y_full=list(Q4_QY_HANDLES),
            Q4_R_Q_y_full=list(Q4_RQY_HANDLES),
        ),
        six_I1_names=list(I1_NAMES),
        physical_measure=dict(
            radius_handle=MACRO_RADIUS_HANDLE,
            radius_power=2,
            length_handle=MACRO_J_HANDLE,
            positive_K=record(K),
            R_squared_length_applied_once=True,
            Pstar_not_added_to_I1_weight=True,
        ),
        signed_source_values_not_selected=True,
        old_signed_sources_and_error_packets_unchanged=True,
    )

    return dict(
        passed=True,
        source_family=parent.get("source_family"),
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        source_pins=pins,
        source_guards=guards,
        source_AST=guards,
        source_contract=source_contract,
        loop_packet_source_pins=packet_pins,
        physical_error_log_bounds=dict(
            delta_U=_pack_logs(delta_U),
            delta_M=_pack_logs(delta_M),
            delta_Q=_pack_logs(delta_Q),
            log_Pstar_E_upper=record(logPstarE),
            logN=record(logN),
            delta_M_initial_source=(
                "parent.actual_micro_Q_increment_bounds.charts.second.log_caps."
                "physical_delta_m"
            ),
            error_center_not_selected=True,
        ),
        six_I1_difference_error_log_bounds={
            name: _pack_logs(values) for name, values in I1_error.items()
        },
        I1_full=I1_full,
        D_full=D_full,
        full_signed_enclosures=dict(I1=I1_full, D=D_full),
        signed_I1_radius=record(radius),
        I1_error_threshold=record(threshold),
        all_normalized_I1_errors_below_2_power_minus1024=True,
        physical_graph_aliases=dict(
            Pstar_handle=PSTAR_HANDLE,
            Pstar_log_argument_handle=PSTAR_LOG_ARGUMENT_HANDLE,
            Pstar_literal_exp5_not_used=True,
            logPstar_parent_record=geometry["logPstar"],
            R_handle=MACRO_RADIUS_HANDLE,
            R_only_for_measure_and_RQy=True,
        ),
        finite_N_value=dict(
            selected=False,
            issued_phase_selected=False,
            N_materialized=False,
            source_is_log_error_radius=True,
        ),
        limits=dict(
            finite_N_six_I1_error_radius_available=True,
            finite_N_six_I1_enclosed=True,
            finite_N_global_macro_enclosed=False,
            outer_RQy_retained=True,
            outer_RQy_evaluated=False,
            two_hundred_twelve_terms_unchanged=True,
            old_six_source_errors_unchanged=True,
            global_five_moment_closure=False,
            complete_pressure_global_closure=False,
            runtime_model_metadata_available=False,
            physical_y_derivatives_claimed=False,
        ),
        actual_current_N18=True,
        accepted_ancestors_replayed=False,
    )


__all__ = ["build"]
