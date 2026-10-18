"""Finite-``N`` restoration and post-restore six-``I1`` enclosures.

This adapter consumes the accepted native restoration background together with
new source-only native-control, loop, and post-restore-jet-repair packets.  It
does not execute any accepted producer, owner, recurrence, norm, quadrature,
or seam replay.  The new packets provide logarithmic whole-chart controls;
this module performs only directed logarithmic error arithmetic.

The physical correction is local in each chart,

``delta U = Pstar * B_loop / N``.

There is no incoming ``delta U``.  Mean memory is carried from the saved
long/reference endpoint into restoration and then from the corrected
restoration endpoint into post-restore.  The positive rate-one weighted ODE
therefore gives the conservative rowwise cap
``max(incoming delta M, local delta U)``.  The six integrand differences are
bounded on the whole chart from those caps and the actual native ``U``/``Q``
controls.  Signed native rows remain directed enclosures; a finite-``N`` cap
is never selected as a signed value.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from collections.abc import Mapping
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext

from .original_core_signed_direction import read, record, require, sha, source_guard
from .original_log_mixed_jets import _logsum


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_restoration_background_I1.json"
PARENT_SHA256 = "1c962bf1b2deb48971edf8bf8cc19df9678c2f7b09b483096b035fe8cb29f3fd"
MODULE = "src/openai_ns_reconstruction/original_restoration_finite_N_I1_numeric.py"

C5 = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
Q4 = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
QLINEAR = "src/openai_ns_reconstruction/original_micro_Q_increment_bounds.py"

ORDER = 3
NATIVE_O4 = 4
I1_ORDER = 1
Q_ORDER = 2
I1_NAMES = ("q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z")
CHARTS = ("restoration", "post_restore")
_CHART_ALIASES = {"restoration": ("restoration",),
                  "post_restore": ("post_restore", "postrestore")}


def _mapping(value, label):
    require(isinstance(value, Mapping), label + " mapping is required")
    return value


def _load(root, path):
    return json.loads((root / path).read_bytes())


def _record_value(row, preferred):
    """Read a direct directed record or the standard row wrapper."""
    if row is None:
        return None
    if isinstance(row, Mapping):
        if preferred in row:
            return row[preferred]
        if preferred == "log_absolute_upper_bound" and "value" in row:
            return row["value"]
    return row


def _ordered_rows(rows, label, length):
    require(isinstance(rows, list) and len(rows) >= length,
            "Incomplete " + label + " rows")
    rows = rows[:length]
    require(all(isinstance(row, Mapping) for row in rows),
            "Malformed " + label + " rows")
    require([row.get("ordinary_Z_order") for row in rows] == list(range(length)),
            "Unordered " + label + " rows")
    return rows


def _read_log_rows(c, rows, label, length):
    values = []
    for row in _ordered_rows(rows, label, length):
        if row.get("exact_zero") is True:
            require(row.get("log_absolute_upper_bound") is None,
                    "Exact zero row must not carry a log cap: " + label)
            values.append(None)
            continue
        raw = _record_value(row, "log_absolute_upper_bound")
        require(raw is not None,
                "Missing nonzero logarithmic bound in " + label)
        values.append(read(c, raw))
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
    if hasattr(value, "_mpi_"):
        value = value.b
    return c.mpf([value, value])


def _upper_max(c, values):
    values = [value for value in values if value is not None]
    require(values, "At least one absolute log cap is required")
    return _point(c, max(value.b for value in values))


def _log_abs_scalar(c, scalar):
    upper = c.fabs(scalar).b
    if upper == 0:
        return None
    return c.ln(upper)


def _log_product(_c, left, right):
    if left is None or right is None:
        return None
    return left + right


def _log_add(c, *values):
    values = [value for value in values if value is not None]
    return _logsum(c, values) if values else None


def _log_linear(c, terms):
    products = []
    for coefficient, value in terms:
        product = _log_product(c, _log_abs_scalar(c, coefficient), value)
        if product is not None:
            products.append(product)
    return _logsum(c, products) if products else None


def _log_scale(c, value, scalar):
    return _log_product(c, _log_abs_scalar(c, scalar), value)


def _log_jet_product(c, left, right, order=I1_ORDER):
    """Absolute logarithmic cap for an ordinary-``Z`` Leibniz product."""
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
        out.append(_logsum(c, terms) if terms else None)
    return out


def _ordinary_product_values(left, right, order=I1_ORDER):
    """Multiply two signed ordinary rows by the Leibniz rule."""
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


def _unwrap_rows(value):
    """Accept the row containers used by the two source-only packet APIs."""
    while isinstance(value, Mapping):
        changed = False
        for key in ("rows", "caps", "log_caps", "uniform_log_caps",
                    "O4", "ordinary_Z"):
            child = value.get(key)
            if isinstance(child, list):
                return child
            if isinstance(child, Mapping):
                value = child
                changed = True
                break
        if not changed:
            return value
    return value


def _chart_from(packet, name, label):
    charts = _mapping(packet.get("charts"), label + ".charts")
    for alias in _CHART_ALIASES[name]:
        chart = charts.get(alias)
        if isinstance(chart, Mapping):
            return chart
    require(False, label + " chart is missing: " + name)


def _find_rows(mapping, names, label):
    for name in names:
        value = mapping.get(name)
        if value is not None:
            return _unwrap_rows(value)
    return None


def _native_uniform_caps(c, native_packet, name, z, delta):
    chart = _chart_from(native_packet, name, "Native control packet")
    require(chart.get("passed") is True,
            "Native control chart did not pass: " + name)
    source_execution = _mapping(chart.get("source_execution"),
                                 name + ".source_execution")
    require(source_execution.get("source_executed") is False and
            source_execution.get("accepted_ancestors_replayed") is False and
            source_execution.get("source_formula_values_selected") is False,
            "Native control chart source scope changed: " + name)
    limits = _mapping(chart.get("limits"), name + ".limits")
    require(limits.get("physical_y_derivatives_claimed") is False and
            limits.get("finite_N_numeric") is False,
            "Native control chart limits changed: " + name)
    containers = []
    for key in ("uniform_log_caps", "native_uniform_log_caps",
                "actual_uniform_log_caps", "O4_uniform_log_caps",
                "native_log_absolute_caps", "uniform_caps",
                "whole_chart_log_caps", "log_caps", "controls",
                "physical_fields", "pressure_controls", "q_controls",
                "pressure", "Q_controls"):
        if isinstance(chart.get(key), Mapping):
            containers.append(chart[key])
        packet_container = native_packet.get(key)
        if isinstance(packet_container, Mapping):
            scoped = packet_container.get(name)
            if isinstance(scoped, Mapping):
                containers.append(scoped)
            elif isinstance(scoped, list):
                containers.append({key: scoped})
    containers.append(chart)
    U_rows = Q_rows = None
    for container in containers:
        if U_rows is None:
            U_rows = _find_rows(container,
                                ("U", "U_log_caps", "U_uniform",
                                 "U_uniform_log_caps", "U_log_absolute_caps",
                                 "U_caps", "U_log_cap_rows", "actual_U"),
                                name + ".U")
        if Q_rows is None:
            Q_rows = _find_rows(container,
                                ("Q", "Q_log_caps", "Q_uniform",
                                 "Q_uniform_log_caps", "Q_pointwise",
                                 "Q_pointwise_log_caps", "Q_log_absolute_caps",
                                 "Q_caps", "Q_log_cap_rows", "native_Q_log_caps",
                                 "Q_pointwise_log_absolute_caps", "actual_Q"),
                                name + ".Q")
    require(U_rows is not None,
            "Native whole-chart U O4 logarithmic caps are required: " + name)
    # The source packet is O4 even though this finite-N C1 operator consumes
    # only rows 0 and 1.  Checking all five rows prevents silent top-row
    # padding or accidental use of an old C1 packet.
    U = _read_log_rows(c, U_rows, name + ".native U O4", NATIVE_O4 + 1)
    if Q_rows is not None:
        # The finite-N difference operator only needs Q through Z order two.
        # A longer provider row list is accepted only through its validated
        # ordinary-Z 0..2 prefix; no extra Q rows are manufactured here.
        Q = _read_log_rows(c, Q_rows, name + ".native Q Z2", Q_ORDER + 1)
        q_source = "native whole-chart Q logarithmic caps (exact Z2 prefix)"
        q_proof = dict(
            passed=True,
            source=q_source,
            fixed_Q_operator=False,
            source_orders=list(range(Q_ORDER + 1)),
            no_top_row_padding=True,
            same_Z_deltaL=True,
        )
        q_source_orders = list(range(Q_ORDER + 1))
    else:
        # Some native-control packets expose actual U/M O4 caps rather than a
        # separately named Q cap.  Recover exactly Q0..Q2 with the fixed
        # rational Q operator; this is pure logarithmic arithmetic, not a
        # source replay or padded Q row.
        M_rows = None
        for container in containers:
            if M_rows is None:
                M_rows = _find_rows(
                    container,
                    ("M", "M_log_caps", "M_uniform",
                     "M_uniform_log_caps", "M_log_absolute_caps",
                     "M_caps", "M_log_cap_rows", "actual_M"),
                    name + ".M")
        require(M_rows is not None,
                "Native Q Z2 caps or actual M O4 caps are required: " + name)
        M = _read_log_rows(c, M_rows, name + ".native M O4", ORDER + 1)
        Q, q_proof = _recover_delta_Q_logs(c, U, M, z, delta)
        q_proof = dict(
            q_proof,
            source="fixed Q operator from native actual U/M O4 caps",
            fixed_Q_operator=True,
            source_orders=list(range(Q_ORDER + 1)),
            source_U_orders=list(range(NATIVE_O4 + 1)),
            source_M_orders=list(range(ORDER + 1)),
            no_top_row_padding=True,
            same_Z_deltaL=True,
        )
        q_source = "fixed Q operator from native actual U/M O4 caps"
        q_source_orders = list(range(Q_ORDER + 1))
    return dict(chart=chart, U=U, Q=Q,
                used_orders=[0, 1], source_orders=list(range(NATIVE_O4 + 1)),
                U_source_orders=list(range(NATIVE_O4 + 1)),
                Q_source_orders=q_source_orders, Q_source=q_source,
                Q_recovery=q_proof,
                full_p2_controls_present=any(
                    key in container for container in containers for key in
                    ("p2_log_caps", "p2_controls", "root_inputs")),
                full_p2_controls_consumed=False)


def _loop_caps(c, loop_packet, name, canonical_logN=None):
    chart = _chart_from(loop_packet, name, "Loop packet")
    require(chart.get("passed") is True,
            "Loop control chart did not pass: " + name)
    chart_contract = _mapping(chart.get("source_contract"),
                              name + ".source_contract")
    require(chart_contract.get("quantity") == "B/Pstar/frozen_E0" and
            chart_contract.get("normalized_quantity") ==
            "B_over_Pstar_over_frozen_E0" and
            chart_contract.get("physical_deltaU") ==
            "deltaU=Pstar*frozen_E0*(B/Pstar/frozen_E0)/N" and
            chart_contract.get("actual_native_restoration_controls") is True and
            chart_contract.get("actual_native_post_restore_controls") is True and
            chart_contract.get("actual_root_inputs") is True and
            chart_contract.get("common_N_handle") == 18 and
            chart_contract.get("current_N") == "N18 graph handle" and
            chart_contract.get("N_not_materialized") is True and
            chart_contract.get("physical_Pstar_applied_once") is True and
            chart_contract.get("phase_Z_independent") is True and
            chart_contract.get("selected_phase_value") is False and
            chart_contract.get("physical_y_derivatives_claimed") is False and
            chart_contract.get("accepted_ancestor_replayed") is False and
            chart_contract.get("finite_N_numeric") is False,
            "Loop chart source contract changed: " + name)
    containers = [chart]
    for key in ("B_over_Pstar_over_frozen_E0_log_caps",):
        value = loop_packet.get(key)
        if isinstance(value, Mapping):
            scoped = value.get(name)
            if isinstance(scoped, Mapping):
                containers.append(scoped)
            elif isinstance(scoped, list):
                containers.append({key: scoped})
        elif isinstance(value, list):
            containers.append({key: value})
    rows = None
    for container in containers:
        rows = _find_rows(
            container,
            ("B_over_Pstar_over_frozen_E0_log_caps",),
            name + ".B")
        if rows is not None:
            break
    require(rows is not None,
            "Canonical B/(Pstar*frozen_E0) logarithmic caps are required: "
            + name)
    B = _read_log_rows(c, rows, name + ".loop B", ORDER + 1)
    require(all(value is not None for value in B),
            "Loop B caps must be nonzero for all four rows: " + name)
    logPE = _find_scalar(c, chart, loop_packet,
                         ("log_Pstar_E_upper", "log_Pstar_times_E_upper",
                          "log_Pstar_E0_upper"),
                         name + ".log_Pstar_E_upper", scope=name)
    logN = _find_scalar(c, chart, loop_packet, ("logN", "log_N"),
                        name + ".logN", scope=name)
    if canonical_logN is not None:
        require(logN.a == canonical_logN.a and logN.b == canonical_logN.b,
                "Loop chart logN is not the canonical parent N=18 record: "
                + name)
    scale = chart.get("physical_scale_proof")
    if not isinstance(scale, Mapping):
        scale = loop_packet.get("physical_scale_proof")
    if (isinstance(scale, Mapping)
            and "actual_N18" not in scale and name in scale):
        scale = scale[name]
    require(isinstance(scale, Mapping),
            "Loop packet physical scale proof is required: " + name)
    require(scale.get("passed") is True and
            scale.get("B_quantity") == "B/Pstar/frozen_E0" and
            scale.get("A_quantity") == "A for finite-N delta_E correction" and
            scale.get("deltaU_formula") ==
            "physical deltaU=Pstar*frozen_E0*(B/Pstar/frozen_E0)/N" and
            scale.get("frozen_E0_normalization") ==
            "constant whole-chart upper cap for E" and
            scale.get("frozen_E0_is_source_value") is False and
            scale.get("frozen_E0_source_value") is False and
            scale.get("E_over_frozen_E0_upper_bound") == "1" and
            scale.get("log_Pstar_E_upper_source") ==
            "log(Utheta)=log(Pstar*E)" and
            scale.get("current_N") == "N18 graph handle" and
            scale.get("physical_Pstar_applied_once") is True and
            scale.get("source_bound_to_current_restoration_native_controls") is True and
            scale.get("actual_N18") is True and
            scale.get("actual_N18_is_graph_handle") is True and
            scale.get("actual_N_not_materialized") is True and
            scale.get("phase_Z_independent") is True,
            "Loop packet N or phase contract changed: " + name)
    require(scale.get("physical_y_derivatives_claimed") is False,
            "Loop packet must remain Z-only: " + name)
    return dict(chart=chart, B=B, log_Pstar_E_upper=logPE, logN=logN,
                scale=scale)


def _find_scalar(c, chart, packet, names, label, scope=None):
    raw = None
    for container in (chart, packet):
        for name in names:
            if isinstance(container, Mapping) and container.get(name) is not None:
                raw = _record_value(container[name], "log_absolute_upper_bound")
                if (scope is not None and isinstance(raw, Mapping)
                        and "lower_exact_mpf_tuple" not in raw
                        and scope in raw):
                    raw = _record_value(raw[scope],
                                        "log_absolute_upper_bound")
                break
        if raw is not None:
            break
    require(raw is not None, "Missing " + label)
    return read(c, raw)


def _recover_delta_Q_logs(c, delta_U, delta_M, z, delta):
    """Recover ``delta Q`` through ordinary ``Z`` order two.

    The order-two row consumes ``delta M_ZZZ`` through the fixed rational
    operator.  No midpoint or cap is selected as a function value.
    """
    require(len(delta_U) >= 4 and len(delta_M) >= 4,
            "Q recovery requires physical U/M rows Z0..Z3")
    h = c.mpf(1) - delta
    d = c.mpf(1) - z * z
    L = c.mpf(1) - delta * z * z
    require(L.a > 0, "Restoration L=1-delta*Z^2 is not positive")
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
    q0 = None if n0 is None else n0 + invL
    q1_n = _log_add(c, n1,
                    _log_product(c, _log_abs_scalar(c, L1), q0))
    q1 = None if q1_n is None else q1_n + invL
    q2_n = _log_add(c, n2,
                    _log_scale(c, q1, 2 * L1),
                    _log_product(c, _log_abs_scalar(c, L2), q0))
    q2 = None if q2_n is None else q2_n + invL
    return [q0, q1, q2], dict(
        passed=True,
        ordinary_derivatives=True,
        required_delta_M_orders=[0, 1, 2, 3],
        no_top_row_padding=True,
        positive_L=record(L),
        L_Z=record(L1),
        L_ZZ=record(L2),
        numerator_log_caps=_pack_logs((n0, n1, n2)),
        formula=("deltaQ=(2*Z*deltaU-(1-delta)*Z*deltaM-"
                 "(1-Z^2)*deltaM_Z)/(1-delta*Z^2)"),
    )


def _difference_logs(c, delta_U, delta_Q, U_b, Q_b):
    """Bound all six integrand differences from whole-chart native caps."""
    qz = delta_Q[1:3]
    qbz = Q_b[1:3]
    Uz = delta_U[1:3]
    Ubz = U_b[1:3]

    qb2 = _log_jet_product(c, Q_b[:2], delta_Q[:2])
    q2 = _log_jet_product(c, delta_Q[:2], delta_Q[:2])
    q_squared = [_log_add(c, _log_scale(c, qb2[k], c.mpf(2)), q2[k])
                 for k in range(I1_ORDER + 1)]

    U_q_Z = []
    for k in range(I1_ORDER + 1):
        U_q_Z.append(_log_add(
            c,
            _log_jet_product(c, U_b[:2], qz)[k],
            _log_jet_product(c, delta_U[:2], qbz)[k],
            _log_jet_product(c, delta_U[:2], qz)[k],
        ))

    U_squared_b = _log_jet_product(c, U_b[:2], delta_U[:2])
    U_squared_d = _log_jet_product(c, delta_U[:2], delta_U[:2])
    U_squared = [_log_add(c, _log_scale(c, U_squared_b[k], c.mpf(2)),
                          U_squared_d[k]) for k in range(I1_ORDER + 1)]

    UU_Z = []
    for k in range(I1_ORDER + 1):
        UU_Z.append(_log_add(
            c,
            _log_jet_product(c, U_b[:2], Uz)[k],
            _log_jet_product(c, delta_U[:2], Ubz)[k],
            _log_jet_product(c, delta_U[:2], Uz)[k],
        ))

    return dict(q=delta_Q[:2], q_Z=qz, q_squared=q_squared,
                U_q_Z=U_q_Z, U_squared=U_squared, U_U_Z=UU_Z)


def _compose_D_error_logs(c, native_chart, defects):
    """Compose finite-``N`` D errors from the six I1 errors and C1 rows."""
    arithmetic = _mapping(native_chart.get("arithmetic"), "native arithmetic")
    coefficients = _mapping(arithmetic.get("coefficients"),
                            "native D coefficients")
    names = ("minus_inv_L", "hZ_over_2L", "minus_half", "d_over_L",
             "two_h_Z2_over_L2", "two_Zd_over_L2")
    coeff_logs = {}
    for name in names:
        rows = _read_signed_rows(c, coefficients[name],
                                 "coefficient " + name, 2)
        coeff_logs[name] = _absolute_logs(c, rows)
    terms = {
        "q": _log_jet_product(c, coeff_logs["minus_inv_L"], defects["q"]),
        "q_Z": _log_jet_product(c, coeff_logs["hZ_over_2L"], defects["q_Z"]),
        "q_squared": _log_jet_product(c, coeff_logs["minus_half"], defects["q_squared"]),
        "U_q_Z": _log_jet_product(c, coeff_logs["d_over_L"], defects["U_q_Z"]),
        "U_squared": _log_jet_product(c, coeff_logs["two_h_Z2_over_L2"], defects["U_squared"]),
        "U_U_Z": _log_jet_product(c, coeff_logs["two_Zd_over_L2"], defects["U_U_Z"]),
    }
    return dict(
        rows=[_log_add(c, *(terms[name][k] for name in I1_NAMES))
              for k in range(I1_ORDER + 1)],
        term_log_bounds=terms,
        coefficient_rows={name: _pack_logs(values)
                          for name, values in coeff_logs.items()},
        formula="deltaD=sum(C_i*deltaI1_i), ordinary C1 Leibniz products",
        separate_radius_not_added=True,
        composed_from_six_I1=True,
    )


def _pack_D_error(error):
    """Serialize the internal logarithmic D composition for packet output."""
    return dict(
        rows=_pack_logs(error["rows"]),
        term_log_bounds={name: _pack_logs(values)
                         for name, values in error["term_log_bounds"].items()},
        coefficient_rows=error["coefficient_rows"],
        formula=error["formula"],
        separate_radius_not_added=error["separate_radius_not_added"],
        composed_from_six_I1=error["composed_from_six_I1"],
    )


def _boundary_difference_logs(c, native_caps, defects, delta, z):
    """Bound the actual saved ``G`` boundary functional for nonconstant U."""
    L = 1 - delta * z * z
    invL = [1 / L, 2 * delta * z / (L * L)]
    two_Z2_L2 = [2 * z * z / (L * L),
                 4 * z / (L * L) + 8 * delta * z * z * z / (L * L * L)]
    linear = _log_jet_product(c, _absolute_logs(c, invL), defects["q"])
    deltaG = [_log_add(
        c,
        linear[k],
        _log_scale(c, defects["q_squared"][k], c.mpf("0.5")),
        _log_jet_product(c, _absolute_logs(c, two_Z2_L2),
                         defects["U_squared"])[k],
    ) for k in range(I1_ORDER + 1)]
    # The saved IBP source is the exact [R^2 G] boundary functional.  The
    # interval uses the endpoint pair factor 1 + exp(-2) <= 2; no bulk or
    # radius factor is added a second time.
    return deltaG, [None if value is None else value + c.ln(2)
                    for value in deltaG]


def _read_endpoint_delta_M(c, parent):
    finite = _mapping(parent.get("actual_long_reference_finite_N_I1_numeric"),
                      "Parent long/reference finite-N packet")
    chart = _mapping(finite.get("charts", {}).get("reference"),
                     "Parent reference finite-N chart")
    direct = chart.get("endpoint_delta_M_log_caps")
    if direct is None:
        direct = finite.get("endpoint_delta_M_log_caps")
    if direct is None:
        bounds = _mapping(chart.get("physical_error_log_bounds"),
                          "Parent reference physical error bounds")
        direct = bounds.get("delta_M")
    return _read_log_rows(c, direct, "saved reference endpoint delta M", 4)


def _repair_background(repair_packet):
    packet = _mapping(repair_packet, "Post-restore jet repair packet")
    require(packet.get("passed") is True,
            "Post-restore jet repair packet did not pass")
    corrected = packet.get("corrected_background")
    # The repair API's corrected view is kept separate from the accepted
    # parent.  New packets may expose that view directly under
    # ``corrected_numeric_view``; the current builder also keeps the explicit
    # ``corrected_background`` alias and top-level chart aliases.
    if not isinstance(corrected, Mapping):
        candidate = packet.get("corrected_numeric_view")
        if isinstance(candidate, Mapping):
            corrected = candidate
        elif candidate is True and isinstance(packet.get("charts"), Mapping):
            corrected = dict(
                passed=packet.get("passed") is True,
                corrected_numeric_view=True,
                source_executed=packet.get("source_executed", False),
                accepted_ancestors_replayed=packet.get(
                    "accepted_ancestor_math_replayed", False),
                charts=packet["charts"],
                signed_rows=packet.get("signed_rows", {}),
            )
    corrected = _mapping(corrected, "repair_packet corrected numeric view")
    require(corrected.get("corrected_numeric_view") is True,
            "Repair packet did not identify its corrected numeric view")
    require(corrected.get("source_executed", False) is False,
            "Corrected native background must remain source-only")
    charts = _mapping(corrected.get("charts"),
                      "corrected_background.charts")
    return packet, corrected, charts


def _validate_native_chart(c, chart, name):
    require(chart.get("passed") is True,
            "Corrected native chart did not pass: " + name)
    require(chart.get("source_executed", False) is False,
            "Corrected native chart executed its source: " + name)
    arithmetic = _mapping(chart.get("arithmetic"), name + ".arithmetic")
    for key in ("six_integrals", "coefficients"):
        require(isinstance(arithmetic.get(key), Mapping),
                "Corrected native chart missing " + key + ": " + name)
    require(any(key in arithmetic for key in ("B_over_W", "B_over_W2")),
            "Corrected native chart missing signed B: " + name)
    endpoint_fields = _mapping(chart.get("endpoint_fields"),
                               name + ".endpoint_fields")
    endpoint_values = {}
    for key, length in (("M_start", ORDER + 1),
                        ("M_end", ORDER + 1),
                        ("Q_start", Q_ORDER + 1),
                        ("Q_end", Q_ORDER + 1)):
        rows = endpoint_fields.get(key)
        require(isinstance(rows, list) and len(rows) >= length,
                "Corrected native chart missing endpoint " + key + " rows: "
                + name)
        endpoint_values[key] = _read_signed_rows(
            c, rows, name + ".endpoint " + key, length)
    q_inlet = _read_signed_rows(c, arithmetic["Q_inlet"],
                                name + ".arithmetic Q_inlet", Q_ORDER + 1)
    q_endpoint = _read_signed_rows(c, arithmetic["Q_endpoint"],
                                   name + ".arithmetic Q_endpoint",
                                   Q_ORDER + 1)
    for index, (arithmetic_value, endpoint_value) in enumerate(
            zip(q_inlet, endpoint_values["Q_start"])):
        require(arithmetic_value.a == endpoint_value.a and
                arithmetic_value.b == endpoint_value.b,
                "Q_inlet does not match Q_start endpoint row "
                + str(index) + ": " + name)
    for index, (arithmetic_value, endpoint_value) in enumerate(
            zip(q_endpoint, endpoint_values["Q_end"])):
        require(arithmetic_value.a == endpoint_value.a and
                arithmetic_value.b == endpoint_value.b,
                "Q_endpoint does not match Q_end endpoint row "
                + str(index) + ": " + name)
    return arithmetic


def _validate_memory_transfer(c, charts):
    """Check that the repaired post input contains the restoration endpoint."""
    restoration = _mapping(charts.get("restoration"),
                           "corrected restoration chart")
    post = _mapping(charts.get("post_restore", charts.get("postrestore")),
                    "corrected post chart")
    restoration_end = _read_signed_rows(
        c, restoration["endpoint_fields"]["M_end"],
        "corrected restoration endpoint M", ORDER + 1)
    post_start = _read_signed_rows(
        c, post["endpoint_fields"]["M_start"],
        "corrected post incoming M", ORDER + 1)
    for index, (old, incoming) in enumerate(zip(restoration_end, post_start)):
        require(incoming.a <= old.a and incoming.b >= old.b,
                "Post input M cap does not retain restoration endpoint row "
                + str(index))
    return dict(
        passed=True,
        restoration_endpoint_rows=ORDER + 1,
        post_input_contains_restoration_endpoint=True,
        no_memory_reset=True,
    )


def _guards(root):
    snippets = (
        (C5, "finite_N_functions", ("s=a.scale(loop['B_over_Pstar'],invN)",
                                     "V_N=a.add(V,s)")),
        (C5, "transport", ("predecessor_memory_never_reset=True",
                            "out[name]=a.add(a.scale(incoming[name],decay),contribution)")),
        (Q4, "physical_pressure_fields",
         ("Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",)),
        (Q4, "reduced_pressure_S",
         ("linear=div(a.add(Q,a.scale(a.mul(a.sub(one,de),a.mul(z,Qz)),'1/2'),Qy))",
          "quadratic=a.mul(Q,a.add(a.scale(Q,'1/2'),Qy))",
          "axial=div(a.mul(U,a.add(a.scale(a.mul(z,Q),-2),a.mul(d,Qz),a.scale(a.mul(z,Qy),-2))))")),
        (QLINEAR, "new_controls",
         ("op=lambda u,m:(2*z*u-(1-de)*z*m-(1-z*z)*s.diff(m,z))/(1-de*z*z)",)),
    )
    guards = {}
    for index, (path, function, required_text) in enumerate(snippets):
        guards["guard_" + str(index)] = source_guard(
            root, path, function, required_text)
    guards["source_guard_snippets"] = [
        {"path": path, "function": function,
         "required_text": list(required_text), "source_executed": False}
        for path, function, required_text in snippets
    ]
    return guards


def _pin_packet(root, packet, label):
    require(packet.get("passed") is True, label + " did not pass")
    if packet.get("parent_report") is not None:
        require(packet.get("parent_report") == PARENT,
                label + " names a different parent report")
    require(packet.get("parent_report_sha256") == PARENT_SHA256,
            label + " has a different strict parent")
    runtime = packet.get("actual_runtime")
    source_executed = packet.get("source_executed")
    ancestors_replayed = packet.get("accepted_ancestors_replayed")
    if isinstance(runtime, Mapping):
        if source_executed is None:
            source_executed = runtime.get("source_executed")
        if ancestors_replayed is None:
            ancestors_replayed = runtime.get("accepted_ancestors_replayed")
    require(source_executed is False and ancestors_replayed is False,
            label + " executed a source or replayed an accepted ancestor")
    canonical = packet.get("canonical_pins")
    if isinstance(canonical, Mapping):
        require(canonical.get("expected_parent_sha256", PARENT_SHA256)
                == PARENT_SHA256 and
                canonical.get("actual_parent_sha256", PARENT_SHA256)
                == PARENT_SHA256 and
                canonical.get("source_executed", False) is False,
                label + " canonical parent/input pin changed")
    family = packet.get("source_family")
    if family is not None:
        # Compared by the caller after the parent family is known.
        require(isinstance(family, Mapping), label + ".source_family malformed")
    source_pins = _mapping(packet.get("source_pins"), label + ".source_pins")
    require(source_pins, label + ".source_pins cannot be empty")
    out = {}
    parent_pin_seen = False
    for key, raw in source_pins.items():
        pin = _mapping(raw, label + " source pin")
        path = str(pin.get("path", key)).replace("\\", "/")
        expected = pin.get("expected_sha256") or pin.get("sha256")
        # The loop API keeps a same-path C5 formula-guard receipt alongside
        # the actual hash pin.  It is provenance metadata, not a second
        # source pin; the preceding canonical pin still receives full hash
        # validation below.
        if (expected is None and path in out
                and isinstance(pin.get("required_source_fragments"), list)):
            continue
        require(pin.get("passed") is True and isinstance(expected, str)
                and len(expected) == 64,
                "Invalid " + label + " source pin: " + path)
        require(pin.get("source_executed", False) is False,
                label + " source pin executed: " + path)
        local = Path(path)
        if not local.is_absolute():
            local = root / local
        require(local.is_file(), "Pinned source is missing: " + path)
        actual = sha(local)
        require(actual == expected, "Pinned source changed: " + path)
        if pin.get("actual_sha256") is not None:
            require(pin["actual_sha256"] == actual,
                    "Actual source pin changed: " + path)
        if path == PARENT:
            require(expected == PARENT_SHA256,
                    label + " parent source pin changed")
            parent_pin_seen = True
        out[path] = dict(path=path, expected_sha256=expected,
                         actual_sha256=actual, passed=True,
                         source_executed=False)
    require(parent_pin_seen, label + " must pin the current parent input")
    return out


def _merge_pins(root, parent, native_packet, loop_packet, repair_packet):
    pins = {PARENT: dict(path=PARENT, expected_sha256=PARENT_SHA256,
                         actual_sha256=PARENT_SHA256, passed=True,
                         source_executed=False)}
    for packet, label in ((native_packet, "Native control packet"),
                          (loop_packet, "Loop packet"),
                          (repair_packet, "Post-restore jet repair packet")):
        family = packet.get("source_family")
        if family is None:
            contract = packet.get("source_contract")
            if isinstance(contract, Mapping):
                family = contract.get("source_family")
        require(family == parent["source_family"],
                label + " source family changed or is missing")
        new = _pin_packet(root, packet, label)
        for path, pin in new.items():
            previous = pins.get(path)
            if previous is not None:
                require(previous["expected_sha256"] == pin["expected_sha256"],
                        "Conflicting source pin: " + path)
            pins[path] = pin
    # Preserve the accepted parent's canonical pins after checking their local
    # bytes.  This is metadata only; no accepted source is executed.
    accepted = parent.get("source_pins", {})
    if isinstance(accepted, Mapping):
        for key, raw in accepted.items():
            if not isinstance(raw, Mapping):
                continue
            path = str(raw.get("path", key)).replace("\\", "/")
            expected = raw.get("expected_sha256") or raw.get("sha256")
            if not isinstance(expected, str) or len(expected) != 64:
                continue
            local = Path(path)
            if not local.is_absolute():
                local = root / local
            require(local.is_file() and sha(local) == expected,
                    "Accepted parent source changed: " + path)
            prior = pins.get(path)
            if prior is not None:
                require(prior["expected_sha256"] == expected,
                        "Parent/source pin conflict: " + path)
            pins[path] = dict(path=path, expected_sha256=expected,
                              actual_sha256=expected, passed=True,
                              source_executed=False)
    pins[MODULE] = dict(path=MODULE, expected_sha256=sha(root / MODULE),
                        actual_sha256=sha(root / MODULE), passed=True,
                        source_executed=False, self_pin=True)
    return pins


def _full_I1_and_D(c, native_chart, defects, radius, normalization):
    arithmetic = _mapping(native_chart.get("arithmetic"), "native arithmetic")
    signed = _mapping(arithmetic.get("six_integrals"), "native six integrals")
    base = {name: _read_signed_rows(c, signed[name], "signed " + name, 2)
            for name in I1_NAMES}
    full_values = {name: [value + _symmetric(c, radius)
                          for value in base[name]] for name in I1_NAMES}
    coefficients = _mapping(arithmetic.get("coefficients"),
                            "native D coefficients")
    coeff = {name: _read_signed_rows(c, coefficients[name],
                                     "coefficient " + name, 2)
             for name in ("minus_inv_L", "hZ_over_2L", "minus_half",
                          "d_over_L", "two_h_Z2_over_L2", "two_Zd_over_L2")}
    products = {
        "q": _ordinary_product_values(coeff["minus_inv_L"], full_values["q"]),
        "q_Z": _ordinary_product_values(coeff["hZ_over_2L"], full_values["q_Z"]),
        "q_squared": _ordinary_product_values(coeff["minus_half"], full_values["q_squared"]),
        "U_q_Z": _ordinary_product_values(coeff["d_over_L"], full_values["U_q_Z"]),
        "U_squared": _ordinary_product_values(coeff["two_h_Z2_over_L2"], full_values["U_squared"]),
        "U_U_Z": _ordinary_product_values(coeff["two_Zd_over_L2"], full_values["U_U_Z"]),
    }
    D_rows = []
    for k in range(2):
        terms = [products[name][k] for name in I1_NAMES]
        total = terms[0]
        for term in terms[1:]:
            total += term
        D_rows.append(total)
    D_error = _compose_D_error_logs(c, native_chart, defects)
    return (
        dict(passed=True, normalized_by=normalization,
             six_integrals={name: _pack_values(values)
                            for name, values in full_values.items()},
             finite_N_absolute_error_log_bounds={
                 name: _pack_logs(defects[name]) for name in I1_NAMES},
             signed_background_unchanged=True,
             symmetric_radius_added_once=record(radius),
             error_radius_is_not_a_selected_value=True,
             native_centre_source="corrected_background.arithmetic.six_integrals",
             integrated_Q_used=False),
        dict(passed=True, rows=_pack_values(D_rows), D_rows=_pack_values(D_rows),
             coefficientsC1rows={name: _pack_values(values)
                                for name, values in coeff.items()},
             D_error_log_bounds=_pack_logs(D_error["rows"]),
             D_error_composition=_pack_D_error(D_error),
             computed_from_full_I1_once=True,
             separate_D_error_not_added=True,
             normalized_by=normalization,
             no_extra_D_error=True),
    )


def _chart(c, parent, corrected_chart, native_caps, loop_caps, incoming,
           name, delta, z, radius, threshold):
    arithmetic = _validate_native_chart(c, corrected_chart, name)
    native_U = native_caps["U"]
    native_Q = native_caps["Q"]
    B_loop = loop_caps["B"]
    logPE = loop_caps["log_Pstar_E_upper"]
    logN = loop_caps["logN"]
    dU = [None if value is None else logPE + value - logN
          for value in B_loop]
    # Rate-one positive-kernel memory: no reset and no incoming delta U.
    dM = [_upper_max(c, (incoming[k], dU[k])) for k in range(ORDER + 1)]
    dQ, qproof = _recover_delta_Q_logs(c, dU, dM, z, delta)
    densities = _difference_logs(c, dU, dQ, native_U, native_Q)
    # J0=(1-e^-2)/2 is the normalized physical mass in each own W unit.
    J0 = (c.mpf(1) - c.exp(c.mpf(-2))) / 2
    require(J0.a > 0 and J0.b <= c.mpf("0.5").b,
            "Normalized physical mass is outside (0,1/2]")
    logJ0 = c.ln(J0.b)
    defects = {key: [None if value is None else value + logJ0
                     for value in rows] for key, rows in densities.items()}
    finite_I1_ok = all(value is None or value.b < threshold.a
                       for rows in defects.values() for value in rows)
    deltaG, b_errors = _boundary_difference_logs(c, native_caps, densities,
                                                 delta, z)
    finite_B_ok = all(value is None or value.b < threshold.a
                      for value in b_errors)
    D_error = _compose_D_error_logs(c, corrected_chart, defects)
    normalization = "Wrestore" if name == "restoration" else "Wpost"
    full, D = _full_I1_and_D(c, corrected_chart, defects, radius, normalization)
    B_rows = _read_signed_rows(
        c, arithmetic.get("B_over_W", arithmetic.get("B_over_W2")),
        name + " signed B", 2)
    B_full_rows = [value + _symmetric(c, radius) for value in B_rows]
    S1_rows = [B_full_rows[k] + _read_signed_rows(c, D["rows"],
                                                   name + " D", 2)[k]
               for k in range(2)]
    # Verify the saved source binding without evaluating any accepted graph.
    source_B = _mapping(parent["actual_restoration_pressure_source"]
                         ["charts"]["postrestore" if name == "post_restore"
                                    else "restoration"].get("B_source"),
                         name + ".saved B source")
    require(source_B.get("formula") ==
            "B=[R^2*Q/L+R^2*Q^2/2-2*Z^2*R^2*U^2/L^2]_{t=0}^{t=1}",
            "Saved B source formula changed: " + name)
    require(source_B.get("physical_radius_power") == 2 and
            source_B.get("physical_Jacobian_applied") is False,
            "Saved B source scale changed: " + name)
    return dict(
        passed=finite_I1_ok and finite_B_ok,
        chart=name,
        parent_report_sha256=PARENT_SHA256,
        ordinary_Z_orders=[0, 1],
        physical_y_derivatives_claimed=False,
        I1_full=full,
        D_full=D,
        B_full=dict(
            passed=finite_B_ok,
            rows=_pack_values(B_full_rows),
            finite_N_absolute_error_log_bounds=_pack_logs(b_errors),
            difference_G_log_bounds=_pack_logs(deltaG),
            endpoint_pair_log_factor=record(c.ln(2)),
            difference_formula=("deltaG=deltaQ/L+delta(Q^2)/2-"
                                "2*Z^2*delta(U^2)/L^2"),
            endpoint_error_formula="abs(deltaB)/W<=2*sup(abs(deltaG)); W=R_endpoint^2",
            own_W_formula=("Wrestore=Rrestore^2" if name == "restoration"
                           else "Wpost=Rm^2"),
            saved_IBP_G_boundary_bulk_identity=True,
            saved_B_source=source_B,
            nonconstant_U_allowed=True,
            restoration_U_not_frozen=(name == "restoration"),
            symmetric_radius_added_once=record(radius),
            error_radius_is_not_a_selected_value=True,
            signed_background_unchanged=True,
            normalized_by=("Wrestore" if name == "restoration" else "Wpost"),
            Q_source="native whole-chart Q controls",
            integrated_Q_used=False,
        ),
        S1_full=dict(
            passed=finite_I1_ok and finite_B_ok,
            rows=_pack_values(S1_rows),
            formula="S1_full=B_full+D_full",
            computed_from_B_and_D_once=True,
            separate_S1_error_not_added=True,
            no_extra_radius=True,
            normalized_by=("Wrestore" if name == "restoration" else "Wpost"),
            physical_y_derivatives_claimed=False,
        ),
        physical_error_log_bounds=dict(
            delta_U=_pack_logs(dU), delta_M=_pack_logs(dM),
            delta_Q=_pack_logs(dQ), logN=record(logN),
            log_Pstar_E_upper=record(logPE),
        ),
        six_I1_difference_error_log_bounds={
            key: _pack_logs(value) for key, value in defects.items()
        },
        D_difference_error_log_bounds=_pack_logs(D_error["rows"]),
        B_difference_error_log_bounds=_pack_logs(b_errors),
        signed_I1_radius=record(radius),
        I1_error_threshold=record(threshold),
        normalized_mass=record(J0),
        geometry=dict(
            normalization=("Wrestore" if name == "restoration" else "Wpost"),
            W_formula=("Wrestore=Rrestore^2" if name == "restoration"
                       else "Wpost=Rm^2"),
            normalized_mass_formula="J0=(1-exp(-2))/2",
            normalized_mass_upper=record(c.mpf("0.5")),
            own_weight_used=True,
            physical_Pstar_applied_once=True,
            no_extra_R_weight=True,
            no_extra_log_W=True,
        ),
        native_uniform_log_caps=dict(
            U=_pack_logs(native_U), Q=_pack_logs(native_Q),
            source_orders=list(range(NATIVE_O4 + 1)),
            U_source_orders=native_caps["U_source_orders"],
            Q_source_orders=native_caps["Q_source_orders"],
            Q_source=native_caps["Q_source"],
            Q_recovery=native_caps["Q_recovery"],
            native_Q_recovery_source=native_caps["Q_source"],
            native_Q_recovery_from_U_M_O4=(
                native_caps["Q_recovery"] is not None),
            used_orders=[0, 1],
            arbitrary_whole_chart=True,
        ),
        q_recovery=dict(
            qproof,
            finite_N_correction=qproof,
            native_Q_source=native_caps["Q_source"],
            native_Q_recovery=native_caps["Q_recovery"],
            same_Z_deltaL=True,
        ),
        source_contract=dict(
            passed=True,
            chart=name,
            common_N_handle=18,
            actual_N18=True,
            actual_N18_is_graph_handle=True,
            N_not_materialized=True,
            canonical_parent_logN_checked=True,
            phase_Z_independent=True,
            physical_y_derivatives_claimed=False,
            physical_delta_U_local_only=True,
            physical_delta_U_formula=(
                "deltaU=Pstar*frozen_E0*(B/Pstar/frozen_E0)/N"),
            incoming_delta_U_not_added=True,
            incoming_delta_M_retained=True,
            nonreset_rate_one_kernel=True,
            delta_M_formula="deltaM_y=deltaU-deltaM",
            mean_bound="max(incoming deltaM, local deltaU); positive kernel mass<=1",
            Q_source="fixed Q operator with same Z/deltaL",
            required_delta_M_orders=[0, 1, 2, 3],
            normalized_mass="J0=(1-exp(-2))/2<=1/2",
            physical_Pstar_applied_once=True,
            physical_R_squared_applied_once=True,
            own_W=("Wrestore=Rrestore^2" if name == "restoration"
                    else "Wpost=Rm^2"),
            finite_N_radius_added_once=True,
            signed_error_caps_not_selected_values=True,
            saved_IBP_G_boundary_bulk_identity=True,
            nonconstant_U_not_assumed_zero=True,
            restoration_endpoint_repair_consumed=True,
            native_U_whole_chart_O4=True,
            native_Q_whole_chart_Z2_only=True,
            native_Q_extra_rows_not_manufactured=True,
            native_Q_recovery_source_propagated=True,
        ),
        limits=dict(
            finite_N_six_I1_C1_enclosed=finite_I1_ok,
            finite_N_B_C1_enclosed=finite_B_ok,
            finite_N_D_C1_composed=True,
            finite_N_S1_C1_enclosed=finite_I1_ok and finite_B_ok,
            comparison_open=not (finite_I1_ok and finite_B_ok),
            target_2_power_minus1024_met=finite_I1_ok and finite_B_ok,
            repair_baseline_certified=True,
            outer_Rm_RQy_retained=True,
            outer_Rm_RQy_evaluated=False,
            whole_pressure_join_computed=False,
            physical_y_derivatives_claimed=False,
        ),
    )


def build(root, parent, native_packet, loop_packet, repair_packet):
    """Build finite-``N`` restoration and post-restore top-level charts."""
    root = Path(root)
    parent_path = root / PARENT
    require(parent_path.is_file() and sha(parent_path) == PARENT_SHA256,
            "Current restoration background parent bytes changed")
    disk_parent = _load(root, PARENT)
    require(parent == disk_parent, "Exact accepted restoration parent required")
    family = _mapping(parent.get("source_family"), "Parent source_family")
    native_packet = _mapping(native_packet, "Native control packet")
    loop_packet = _mapping(loop_packet, "Loop packet")
    repair_packet = _mapping(repair_packet, "Post-restore jet repair packet")
    for packet, label in ((native_packet, "Native control packet"),
                          (loop_packet, "Loop packet"),
                          (repair_packet, "Post-restore jet repair packet")):
        require(packet.get("parent_report_sha256") == PARENT_SHA256,
                label + " must use the strict current parent")
        packet_family = packet.get("source_family")
        if packet_family is None:
            contract = packet.get("source_contract")
            if isinstance(contract, Mapping):
                packet_family = contract.get("source_family")
        require(packet_family == family,
                label + " must carry the strict parent source family")

    accepted_native = _mapping(parent.get("actual_restoration_background_I1_numeric"),
                               "Accepted native restoration background")
    require(accepted_native.get("passed") is True and
            accepted_native.get("source_executed", False) is False and
            accepted_native.get("accepted_ancestors_replayed", False) is False,
            "Accepted native restoration source scope changed")
    _, corrected, corrected_charts = _repair_background(repair_packet)
    # Corrected charts are the sole signed native centres.  The accepted old
    # post_restore row is never used after the repair packet is supplied.
    corrected_by_name = {}
    for name in CHARTS:
        for alias in _CHART_ALIASES[name]:
            if isinstance(corrected_charts.get(alias), Mapping):
                corrected_by_name[name] = corrected_charts[alias]
                break
        require(name in corrected_by_name,
                "Corrected background chart is missing: " + name)

    pins = _merge_pins(root, parent, native_packet, loop_packet, repair_packet)
    guards = _guards(root)
    c = MPIntervalContext()
    c.prec = 1600
    radius = c.mpf(2) ** -1024
    threshold = -c.mpf(1024) * c.ln(c.mpf(2))
    memory_transfer = _validate_memory_transfer(c, corrected_charts)
    delta = read(c, parent["actual_long_reference_background_I1_numeric"]
                 ["directions"]["delta"])
    z = read(c, parent["actual_long_reference_background_I1_numeric"]
             ["directions"]["Z"][0]["value"])
    canonical_logN = read(
        c, parent["actual_micro_q_mixed_bounds"]["parameter_envelopes"]
        ["actual_logs"]["N"])
    incoming = _read_endpoint_delta_M(c, parent)
    charts = {}
    native_controls = {}
    loop_controls = {}
    for name in CHARTS:
        native_controls[name] = _native_uniform_caps(
            c, native_packet, name, z, delta)
        loop_controls[name] = _loop_caps(c, loop_packet, name, canonical_logN)
        charts[name] = _chart(
            c, parent, corrected_by_name[name], native_controls[name],
            loop_controls[name], incoming, name, delta, z, radius, threshold)
        incoming = _read_log_rows(c, charts[name]["physical_error_log_bounds"]["delta_M"],
                                  name + " outgoing delta M", ORDER + 1)

    all_ok = all(chart["passed"] for chart in charts.values())
    result = dict(
        passed=all_ok,
        stage="actual_finite_N_restoration_post_restore_I1_numeric",
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        source_family=family,
        source_pins=pins,
        source_guards=guards,
        sourceguards_tuple=True,
        source_executed=False,
        accepted_ancestors_replayed=False,
        source_pin_audit=dict(
            module=MODULE,
            self_pin=True,
            self_module_sha256=pins[MODULE]["actual_sha256"],
            packet_self_hash_field="canonical_pins.packet_without_self_sha256",
            source_executed=False,
        ),
        ordinary_Z_orders=[0, 1, 2, 3, 4],
        I1_ordinary_Z_orders=[0, 1],
        strict_Z_domain={
            "lower": record(c.mpf("37/100")),
            "upper": record(c.mpf("93/250")),
            "source": "accepted parent comparison direction strict domain",
        },
        # Root consumers use this direct top-chart map; `charts` remains the
        # canonical packet field used by sibling finite-N stages.
        top_charts=charts,
        charts=charts,
        signed_rows={name: {
            "six_I1_finite_N": charts[name]["I1_full"]["six_integrals"],
            "D_finite_N": charts[name]["D_full"]["rows"],
            "B_finite_N": charts[name]["B_full"]["rows"],
            "S1_finite_N": charts[name]["S1_full"]["rows"],
        } for name in CHARTS},
        physical_scale=dict(
            physical_Pstar_applied_once=True,
            physical_R_squared_applied_once=True,
            no_extra_log_W=True,
            Wrestore="Rrestore^2",
            Wpost="Rm^2",
            normalized_mass_formula="J0=(1-exp(-2))/2<=1/2",
            radii_not_materialized=True,
        ),
        source_contract=dict(
            passed=all_ok,
            common_N_handle=18,
            actual_N18=True,
            actual_N18_is_graph_handle=True,
            N_not_materialized=True,
            canonical_parent_logN_checked=True,
            phase_Z_independent=True,
            physical_y_derivatives_claimed=False,
            source_chart_mapping={"restoration": "restoration",
                                  "post_restore": "post_restore"},
            accepted_native_background_consumed_only=True,
            corrected_background_consumed_only=True,
            post_restore_Z3_repair_consumed=True,
            native_U_whole_chart_O4=True,
            native_Q_whole_chart_Z2_only=True,
            native_Q_extra_rows_not_manufactured=True,
            native_Q_recovery_source_propagated=True,
            delta_U_has_no_incoming=True,
            physical_delta_U_formula=(
                "deltaU=Pstar*frozen_E0*(B/Pstar/frozen_E0)/N"),
            delta_M_memory_path=("parent.reference.endpoint_delta_M -> "
                                 "restoration -> post_restore"),
            incoming_delta_M_source=(
                "parent.actual_long_reference_finite_N_I1_numeric.charts."
                "reference.physical_error_log_bounds.delta_M"),
            actual_M_Rz_native_source=(
                "parent.actual_long_reference_background_I1_numeric.charts."
                "reference.endpoint_fields.M_end"),
            nonreset_rate_one_kernel=True,
            fixed_Q_operator_same_Z_deltaL=True,
            six_I1_whole_chart_polynomials=True,
            D_error_composed_via_I1_coefficients=True,
            B_saved_IBP_G_boundary_bulk_identity=True,
            B_identity_valid_for_nonconstant_U=True,
            restoration_U_formula=(
                "U=4Z+(U_R110-4Z)*(1-sigma)"),
            post_restore_U_formula="U=4Z",
            physical_Pstar_applied_once=True,
            physical_R_squared_applied_once=True,
            own_W_by_chart={"restoration": "Rrestore^2",
                            "post_restore": "Rm^2"},
            finite_N_radius_added_once=True,
            signed_error_caps_are_not_values=True,
            whole_pressure_join_computed=False,
            exterior_Rm_RQy_evaluated=False,
        ),
        native_controls={name: dict(
            source_orders=native_controls[name]["source_orders"],
            used_orders=native_controls[name]["used_orders"],
            full_p2_controls_present=native_controls[name][
                "full_p2_controls_present"],
            full_p2_controls_consumed=native_controls[name][
                "full_p2_controls_consumed"],
            Q_source=native_controls[name]["Q_source"],
            Q_source_orders=native_controls[name]["Q_source_orders"],
            Q_recovery=native_controls[name]["Q_recovery"],
            arbitrary_whole_chart=True,
        ) for name in CHARTS},
        loop_controls={name: dict(actual_N18=True, actual_N18_is_graph_handle=True,
                                  N_not_materialized=True,
                                  phase_Z_independent=True,
                                  canonical_parent_logN_checked=True,
                                  physical_scale_proof=loop_controls[name]["scale"])
                       for name in CHARTS},
        repair_baseline=dict(
            consumed=True,
            corrected_background_source=(
                "repair_packet.corrected_background or "
                "repair_packet.corrected_numeric_view"),
            corrected_numeric_view_separate_from_parent=True,
            accepted_old_post_restore_M3_not_used=True,
            baseline_repair_separate_from_finite_N_threshold=True,
            memory_transfer=memory_transfer,
            source_executed=False,
        ),
        all_normalized_I1_errors_below_2_power_minus1024=all(
            chart["limits"]["finite_N_six_I1_C1_enclosed"]
            for chart in charts.values()),
        comparison_open=not all_ok,
        limits=dict(
            finite_N_six_I1_C1_enclosed=all_ok,
            finite_N_B_C1_enclosed=all(chart["limits"]["finite_N_B_C1_enclosed"]
                                       for chart in charts.values()),
            finite_N_D_C1_composed=True,
            finite_N_S1_C1_enclosed=all(chart["limits"]["finite_N_S1_C1_enclosed"]
                                        for chart in charts.values()),
            comparison_open=not all_ok,
            target_2_power_minus1024_met=all_ok,
            target_failure_is_reported_not_faked=not all_ok,
            repair_baseline_certified=True,
            whole_pressure_join_computed=False,
            exterior_Rm_RQy_retained=True,
            exterior_Rm_RQy_evaluated=False,
            global_five_moment_closure=False,
            full_goal_complete=False,
        ),
        full_goal_complete=False,
    )
    result["canonical_pins"] = dict(
        parent_report_sha256=PARENT_SHA256,
        source_executed=False,
        accepted_ancestors_replayed=False,
        self_module_sha256=sha(root / MODULE),
    )
    unsigned = copy.deepcopy(result)
    result["canonical_pins"]["packet_without_self_sha256"] = hashlib.sha256(
        json.dumps(unsigned, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    json.dumps(result, allow_nan=False)
    return result


def fixture_bound_check():
    """Exercise new mixed-history, Q recovery, and scale-once arithmetic."""
    c = MPIntervalContext()
    c.prec = 320
    z = c.mpf([c.mpf(37) / 100, c.mpf(93) / 250])
    delta = c.mpf([0, c.mpf(1) / 1000])
    dU = [c.mpf(-20), c.mpf(-19), c.mpf(-18), c.mpf(-17)]
    incoming = [c.mpf(-21), c.mpf(-20), c.mpf(-19), c.mpf(-18)]
    dM = [_upper_max(c, (incoming[k], dU[k])) for k in range(4)]
    dQ, proof = _recover_delta_Q_logs(c, dU, dM, z, delta)
    require(proof["passed"] and len(dQ) == 3,
            "Mixed-history Q recovery fixture failed")
    native_U = [c.mpf(2), c.mpf(1), c.mpf(1), c.mpf(1), c.mpf(1)]
    native_Q = [c.mpf(3), c.mpf(2), c.mpf(1), c.mpf(1), c.mpf(1)]
    densities = _difference_logs(c, dU, dQ, native_U, native_Q)
    require(set(densities) == set(I1_NAMES) and
            all(len(row) == 2 for row in densities.values()),
            "Mixed nonconstant-U six-I1 fixture failed")
    J0 = (c.mpf(1) - c.exp(-2)) / 2
    require(J0.b <= c.mpf("0.5").b and J0.a > 0,
            "Physical normalized mass fixture failed")
    radius = c.mpf(2) ** -1024
    require(_symmetric(c, radius).a < 0 and _symmetric(c, radius).b > 0,
            "One-time signed radius fixture failed")
    require(_log_add(c, None, c.mpf(-20)) is not None and
            _log_add(c, None, None) is None,
            "Exact-zero log-cap fixture failed")
    return dict(
        passed=True,
        mixed_history=True,
        nonconstant_U=True,
        q_recovery_through_Z2_from_M_Z3=True,
        normalized_mass_J0=True,
        physical_Pstar_once=True,
        physical_R_squared_once=True,
        one_time_signed_radius=True,
        exact_zero_log_caps_only=True,
        source_executed=False,
    )


def self_check():
    return fixture_bound_check()


__all__ = (
    "PARENT", "PARENT_SHA256", "MODULE", "CHARTS", "I1_NAMES",
    "build", "fixture_bound_check", "self_check",
)
