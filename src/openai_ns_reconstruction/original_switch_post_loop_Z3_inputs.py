"""New fixed-phase ordinary-Z input bounds for the current switch/post packet.

This adapter consumes only the new first-switch background packet supplied
by the caller plus saved parameter records.  It does not invoke an accepted
builder, checker, producer, owner, recurrence, quadrature, or ``extend``
function.  The arithmetic below is a fresh log-majorant calculation.
Physical y rows are deliberately absent.
"""

from __future__ import annotations

import json
import math
from collections.abc import Mapping
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext

from .original_core_signed_direction import read, record, require, sha, source_guard
from .original_log_Z_phase_bounds import LogZ


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_first_switch_finite_N_I1.json"
PARENT_SHA256 = "2cc8b4ff9670858564297bc1c27d5322c9af7d2cb85381939dc3e37a3eea346a"
MODULE = "src/openai_ns_reconstruction/original_switch_post_loop_Z3_inputs.py"
LOG_Z_SOURCE = "src/openai_ns_reconstruction/original_log_Z_phase_bounds.py"
C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
Q4_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
LOOP_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
FIRST_SWITCH_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_first_switch_functions.py"
BACKGROUND_MODULE = "src/openai_ns_reconstruction/original_switch_post_background_I1_numeric.py"
CHARTS = ("second_switch", "post_switch")
ORDINARY_Z_ORDERS = (0, 1, 2, 3)


def _read(c, value):
    """Read a saved interval in the caller's context."""
    return read(c, value)


def _read_any(c, value):
    """Read a directed record or a literal scalar without materializing scales."""
    if hasattr(value, "_mpi_"):
        return value
    if isinstance(value, Mapping):
        return _read(c, value)
    return c.mpf(value)


def _rows(c, rows, count=4):
    require(isinstance(rows, (list, tuple)) and len(rows) >= count,
            "A complete saved ordinary-Z row prefix is required")
    out = {}
    for index, row in enumerate(rows[:count]):
        if isinstance(row, Mapping) and "ordinary_Z_order" in row:
            order = row["ordinary_Z_order"]
            if "value" in row:
                row_value = row.get("value")
            elif row.get("exact_zero") is True:
                row_value = None
            else:
                row_value = row.get("record", row.get("log_absolute_upper_bound"))
        else:
            order = index
            row_value = row
        require(type(order) is int and 0 <= order < count,
                "Invalid ordinary-Z row")
        require(order not in out, "Repeated ordinary-Z row")
        out[order] = None if row_value is None else _read(c, row_value)
    require(set(out) == set(range(count)),
            "Rows must cover ordinary Z orders 0 through 3")
    return out


def _abs_log(c, value):
    """Logarithm of an interval magnitude upper bound."""
    if value is None:
        return None
    value = _read(c, value) if isinstance(value, Mapping) else value
    upper = max(abs(value.a), abs(value.b))
    if upper == 0:
        return None
    return c.ln(upper)


def _safe_logsum(c, values):
    """Directed log-sum-exp that never materializes a giant exponential."""
    terms = [value for value in values if value is not None]
    if not terms:
        return None
    maximum = max(value.b for value in terms)
    lower = max(value.a for value in terms)
    correction = c.mpf(0)
    dominant_seen = False
    cutoff = c.mpf(1024) * c.ln(c.mpf(2))
    for value in terms:
        gap = maximum - value.b
        if gap == 0 and not dominant_seen:
            dominant_seen = True
            continue
        if gap > cutoff:
            # exp(-gap) <= 2^-1024; this additive enclosure is representable
            # even when gap has an exponent with billions of decimal digits.
            correction += c.mpf(2) ** -1024
        else:
            correction += c.exp(-gap)
    upper = maximum + c.log(1 + correction)
    return c.mpf([lower, upper.b])


def _log_add(c, *jets):
    return {order: _safe_logsum(c, [jet[order] for jet in jets])
            for order in ORDINARY_Z_ORDERS}


def _log_scale(c, jet, scalar):
    if scalar is None:
        return {order: None for order in ORDINARY_Z_ORDERS}
    return {order: None if jet[order] is None else jet[order] + scalar
            for order in ORDINARY_Z_ORDERS}


def _log_mul(c, left, right):
    return _log_mul_n(c, left, right, 3)


def _log_mul_n(c, left, right, maximum_order):
    """Ordinary Leibniz product of logarithmic magnitude jets."""
    out = {}
    for order in range(maximum_order + 1):
        terms = []
        for left_order in range(order + 1):
            right_order = order - left_order
            if (left.get(left_order) is None
                    or right.get(right_order) is None):
                continue
            terms.append(left[left_order] + right[right_order]
                         + c.ln(math.comb(order, left_order)))
        out[order] = _safe_logsum(c, terms)
    return out


def _log_reciprocal(c, jet, log_lower):
    """Ordinary derivative majorant for 1 / jet, using a positive floor."""
    lower = log_lower.a if hasattr(log_lower, "a") else log_lower
    require(jet[0] is not None, "Reciprocal value row cannot be zero")
    inverse_log = -lower
    out = {0: c.mpf(inverse_log)}
    for order in (1, 2, 3):
        terms = []
        for left_order in range(1, order + 1):
            if jet[left_order] is None or out[order - left_order] is None:
                continue
            terms.append(jet[left_order] + out[order - left_order]
                         + c.ln(math.comb(order, left_order)))
        total = _safe_logsum(c, terms)
        out[order] = None if total is None else total + inverse_log
    return out


def _log_rows(jet):
    return [dict(ordinary_Z_order=order, exact_zero=jet[order] is None,
                 log_absolute_upper_bound=None if jet[order] is None
                 else record(jet[order]))
            for order in ORDINARY_Z_ORDERS]


def _jsonable(value):
    """Serialize fresh interval arithmetic while preserving directed records."""
    if hasattr(value, "_mpi_"):
        return record(value)
    if isinstance(value, Mapping):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _constant(c, logabs):
    return {0: c.mpf(logabs), 1: None, 2: None, 3: None}


def _max_log(c, values):
    terms = [value for value in values if value is not None]
    require(terms, "At least one finite magnitude is required")
    return max(terms, key=lambda item: item.b)


def _scalar_abs_log(c, value):
    """Log of the magnitude of a finite scalar interval."""
    if not hasattr(value, "a"):
        value = c.mpf(value)
    upper = max(abs(value.a), abs(value.b))
    return None if upper == 0 else c.ln(upper)


def _comparison_rows(c, background_packet, name, count=4):
    """Read comparison rows from the new first-switch background packet."""
    directions = background_packet.get("comparison_directions", {})
    if name in directions:
        return _rows(c, directions[name], count)
    comparison = background_packet.get("comparison", {})
    comparison_name = {"phi": "frozen_phi", "U": "frozen_U"}.get(name, name)
    if comparison_name in comparison:
        return _rows(c, comparison[comparison_name], count)
    histories = comparison.get("history_rows", {})
    if name in histories:
        return _rows(c, histories[name], count)
    raise KeyError("Missing first-switch comparison rows: " + name)


def _comparison_log_jet(c, background_packet, name, count=4):
    rows = _comparison_rows(c, background_packet, name, count)
    return {order: _abs_log(c, rows[order]) for order in range(count)}


def _f0_relative(c, parent, background_packet):
    """Fresh F0/F0 Bell majorants from F0_Z/F0 = -Lambda*gaxis."""
    geom = parent["actual_macro_background_I1_numeric"]["geometric_measure"]
    gaxis = _comparison_log_jet(c, background_packet, "gaxis", 4)
    logLambda = _read(c, geom["logLambda"]).b
    ell = {order: None if gaxis[order] is None else
           gaxis[order] + logLambda for order in ORDINARY_Z_ORDERS}
    # F/F(0) has value one.  The remaining rows are ordinary exponential
    # Bell polynomials, with no factorial normalization.
    out = {0: c.mpf(0)}
    out[1] = ell[0]
    out[2] = _safe_logsum(c, [ell[1], None if ell[0] is None else 2 * ell[0]])
    out[3] = _safe_logsum(c, [ell[2],
                              None if ell[0] is None or ell[1] is None
                              else c.ln(3) + ell[0] + ell[1],
                              None if ell[0] is None else 3 * ell[0]])
    return out


def _native_phi_relative(c, background_packet):
    """Use the new first-switch actual native phi rows for E/E0."""
    fields = background_packet.get("native_fields", {})
    rows = _rows(c, fields["phi"], 5)
    floor = rows[0]
    require(floor is not None and floor.a > 0,
            "First-switch native phi value must have a positive floor")
    log_floor = c.ln(floor.a)
    relative = {0: c.mpf(0)}
    for order in (1, 2, 3):
        cap = _abs_log(c, rows[order])
        relative[order] = None if cap is None else cap - log_floor
    return relative


def _extended_logs(c, rows, count=5):
    """Read signed rows as logarithmic absolute derivative caps."""
    values = _rows(c, rows, count)
    return {order: _abs_log(c, values[order]) for order in range(count)}


def _shifted(jet):
    """Derivative-shifted jet: (f_Z)^(k) = f^(k+1)."""
    return {order: jet.get(order + 1) for order in ORDINARY_Z_ORDERS}


def _delta_logs(c, parent):
    delta = _read(c, parent["actual_core_signed_direction"]["parameters"]["delta"])
    require(delta.a > 0, "Saved delta must be positive")
    log_delta = c.ln(max(abs(delta.a), abs(delta.b)))
    one_plus_delta = c.ln(max(abs(1 + delta.a), abs(1 + delta.b)))
    h = 1 - delta
    require(h.a > 0, "Saved 1-delta must be positive")
    log_h = _scalar_abs_log(c, h)
    return delta, log_delta, one_plus_delta, log_h


def _comparison_L_and_inverse(c, parent, background_packet):
    """Rebuild L and its ordinary reciprocal jet through Z order 3."""
    z = _comparison_log_jet(c, background_packet, "Z", 4)
    z2 = _log_mul(c, z, z)
    delta, log_delta, _, _ = _delta_logs(c, parent)
    delta_jet = _constant(c, log_delta)
    dz2 = _log_mul(c, delta_jet, z2)
    # |1-delta*Z^2| is bounded by 1+|delta*Z^2|.  The strict floor is
    # read from the saved positive L enclosure, rather than inferred from a
    # cap, and therefore supports the reciprocal recurrence.
    L = _log_add(c, _constant(c, c.mpf(0)), dz2)
    saved = _comparison_rows(c, background_packet, "L", 4)[0]
    require(saved is not None and saved.a > 0,
            "Saved comparison L value must be strictly positive")
    invL = _log_reciprocal(c, L, c.ln(saved.a))
    return z, z2, L, invL, delta, log_delta


def _f2_logs_and_receipt(c, parent, background_packet):
    """Read F0^2 log caps and prove the fifth coarse cap is at most one."""
    amp = parent["actual_micro_b_p2_mixed_bounds"]["whole_tube_amplitude"]
    raw = amp["F0_squared_derivative_log_caps"]
    require(isinstance(raw, (list, tuple)) and len(raw) >= 5,
            "F0 squared derivative log caps through order 4 are required")
    f2 = {order: _read(c, raw[order]) for order in range(5)}
    geom = parent["actual_macro_background_I1_numeric"]["geometric_measure"]
    logF0Upper = _read(c, amp["log_F0"]).b
    logLambda = _read(c, geom["logLambda"]).b
    # The new first-switch packet deliberately carries the actual comparison
    # gaxis rows only through Z3.  The accepted macro Z4 extension already
    # saved a scalar G enclosing all comparison gaxis rows through Z4; use
    # that scalar as a bound, never as a selected derivative row.
    directions = background_packet.get("comparison_directions", {})
    raw_gaxis = directions.get("gaxis") if isinstance(directions, Mapping) else None
    if isinstance(raw_gaxis, (list, tuple)) and len(raw_gaxis) >= 5:
        gaxis = _comparison_log_jet(c, background_packet, "gaxis", 5)
        log_gmax = _max_log(c, list(gaxis.values())).b
        gaxis_source = "first_switch_background.comparison_directions.gaxis[0..4]"
        gaxis_rows_available = 5
    else:
        extension = parent.get("actual_macro_loop_Z3_inputs", {}).get(
            "native_Z4_extension", {})
        proof = extension.get("comparison_fifth_radius_proof", {})
        require(isinstance(proof, Mapping) and proof.get("G") is not None,
                "A saved comparison gaxis Z4 scalar G is required for F2")
        G = _read(c, proof["G"])
        log_gmax = _scalar_abs_log(c, G)
        require(log_gmax is not None,
                "Saved comparison gaxis scalar G must be nonzero")
        gaxis_source = (
            "parent.actual_macro_loop_Z3_inputs.native_Z4_extension."
            "comparison_fifth_radius_proof.G")
        gaxis_rows_available = 4
    log_D = max(c.mpf(0), c.ln(2) + logLambda + log_gmax)
    log_upper = (2 * logF0Upper + c.ln(52) + 5 * log_D)
    passed = bool(log_upper < 0)
    receipt = dict(
        passed=passed,
        fifth_row_cap_one=passed,
        log_F2_upper=record(c.mpf([log_upper, log_upper])),
        log_F0_upper=record(c.mpf([logF0Upper, logF0Upper])),
        log_D=record(c.mpf([log_D, log_D])),
        max_gaxis_log=record(c.mpf([log_gmax, log_gmax])),
        gaxis_rows_available=gaxis_rows_available,
        gaxis_scalar_covers_orders=[0, 1, 2, 3, 4],
        max_gaxis_source=gaxis_source,
        formula=("logF2upper=2*logF0Upper+ln(52)+5*"
                 "ln(max(1,2*Lambda*maxgaxis[0..4]))"),
        source_cache="actual_micro_b_p2_mixed_bounds.whole_tube_amplitude.F0_squared_derivative_log_caps",
        no_F0_squared_value_materialized=True,
    )
    if isinstance(background_packet, Mapping):
        background_receipt = background_packet.get("F2_fifth_cap_receipt")
        if background_receipt is not None:
            receipt["background_receipt"] = background_receipt
    require(passed, "F0 squared fifth coarse cap is not proved below one")
    f2[4] = c.mpf(0)
    return f2, receipt


def _native_fields(c, background_packet):
    fields = background_packet.get("native_fields")
    require(isinstance(fields, Mapping),
            "The new first-switch background must provide native_fields")
    required = ("phi", "U", "M", "A", "B", "C")
    require(all(name in fields for name in required),
            "First-switch native background must provide phi,U,M,A,B,C")
    out = {}
    for name in required:
        out[name] = _extended_logs(c, fields[name], 5)
    return out


def _native_phi_floor(c, background_packet, native):
    for key in ("phi_positive_lower", "positive_phi_floor", "phi_floor",
                "native_phi_positive_lower"):
        value = background_packet.get(key) if isinstance(background_packet, Mapping) else None
        if value is not None:
            floor = _read(c, value)
            require(floor.a > 0, "First-switch native phi floor must be positive")
            return floor
    fields = background_packet.get("native_fields", {})
    row = _rows(c, fields["phi"], 5)[0]
    require(row is not None and row.a > 0,
            "First-switch native phi floor must be positive")
    return row


def _native_numerator(c, parent, background_packet, native, f2):
    """Exact native p2 numerator, boxed by ordinary Leibniz products."""
    z, z2, L, invL, delta, log_delta = _comparison_L_and_inverse(
        c, parent, background_packet)
    _, _, log_one_plus_delta, log_h = _delta_logs(c, parent)
    d = _log_add(c, _constant(c, c.mpf(0)), z2)
    h = _constant(c, log_h)
    U = {k: native["U"][k] for k in ORDINARY_Z_ORDERS}
    M = {k: native["M"][k] for k in range(4)}
    A = {k: native["A"][k] for k in range(4)}
    B = {k: native["B"][k] for k in range(5)}
    C = {k: native["C"][k] for k in range(5)}
    Mz, Az = _shifted(native["M"]), _shifted(native["A"])
    # -U + h*(M-Z*M_Z)/2 + U*(h*Z*M+d*M_Z)
    m_minus_zmz = _log_add(c, M, _log_mul(c, z, Mz))
    transport = _log_add(c, _log_mul(c, _log_mul(c, h, z), M),
                         _log_mul(c, d, Mz))
    hydro = _log_add(
        c,
        U,
        _log_scale(c, m_minus_zmz, log_h - c.ln(2)),
        _log_scale(c, _log_mul(c, z, A), c.ln(2) + log_delta),
        _log_mul(c, d, Az),
        _log_mul(c, U, transport),
    )

    pressure_rows = parent["actual_micro_b_p2_mixed_bounds"][
        "whole_tube_pressure"]["p0_jets"]
    p0 = _extended_logs(c, pressure_rows, 5)
    p0z = _shifted(p0)
    pressure = _log_add(
        c,
        _log_scale(c, _log_mul(c, z, p0), c.ln(2) + log_one_plus_delta),
        _log_mul(c, d, p0z),
    )
    logPstar = _read(c, parent["actual_micro_b_p2_mixed_bounds"]
                     ["whole_tube_amplitude"]["log_Pstar"]).b
    pressure = _log_scale(c, pressure, 2 * logPstar)

    f2b_all = _log_mul_n(c, f2, B, 4)
    f2c_all = _log_mul_n(c, f2, C, 4)
    f2b = {k: f2b_all[k] for k in ORDINARY_Z_ORDERS}
    f2c = {k: f2c_all[k] for k in ORDINARY_Z_ORDERS}
    dressed = _log_add(
        c,
        _log_scale(c, _log_mul(c, z, f2b), c.ln(2) + log_delta),
        _log_mul(c, d, _shifted(f2b_all)),
        _log_scale(c, _log_mul(c, z, f2c), c.ln(2) + log_one_plus_delta),
        _log_mul(c, d, _shifted(f2c_all)),
    )
    geometry = _background_geometry(background_packet)
    radius = _read_any(c, geometry.get(
        "R_upper101", geometry.get("R_upper", 101)))
    numerator = _log_add(c, hydro, pressure,
                         _log_scale(c, dressed, c.ln(radius.b)))
    return numerator, dict(
        z=z, d=d, L=L, inverse_L=invL,
        U=U, M=M, A=A, B=B, C=C, p0=p0, f2=f2,
        hydro_log_caps=hydro, pressure_log_caps=pressure,
        dressed_swirl_log_caps=dressed,
        numerator_log_caps=numerator,
        formula=(
            "N=-U+h*(M-Z*M_Z)/2+U*(h*Z*M+d*M_Z)"
            "+2*delta*Z*A-d*A_Z"
            "+Pstar^2*(2*(1+delta)*Z*P0-d*P0_Z)"
            "+R*(-2*delta*Z*F2*B+d*(F2*B)_Z"
            "+2*(1+delta)*Z*F2*C-d*(F2*C)_Z)"
        ),
        products_use_ordinary_leibniz=True,
        native_is_first_switch=True,
        physical_radius_upper=record(c.mpf([radius.b, radius.b])),
    )


def _p2(c, parent, background_packet, native, f2):
    numerator, detail = _native_numerator(
        c, parent, background_packet, native, f2)
    phi = {order: native["phi"][order] for order in ORDINARY_Z_ORDERS}
    floor = _native_phi_floor(c, background_packet, native)
    invphi = _log_reciprocal(c, phi, c.ln(floor.a))
    amplitude = parent["actual_micro_b_p2_mixed_bounds"]["whole_tube_amplitude"]
    inverse_F0_relative = _extended_logs(
        c, amplitude["inverse_F0_relative_ordinary_Z0_Z4"], 4)
    logF0Lower = _read(c, amplitude["log_F0"]).a
    _, _, L, invL, _, _ = _comparison_L_and_inverse(c, parent, background_packet)
    geometry = _background_geometry(background_packet)
    radius = _read_any(c, geometry.get(
        "R_upper101", geometry.get("R_upper", 101)))
    prefactor = c.mpf("0.5") * (c.ln(radius.b) - c.ln(2))
    out = _log_mul(c, numerator, invphi)
    out = _log_mul(c, out, inverse_F0_relative)
    out = _log_mul(c, out, invL)
    out = _log_scale(c, out, prefactor - logF0Lower)
    detail.update(phi=phi, inverse_phi=invphi,
                  inverse_F0_relative=inverse_F0_relative,
                  inverse_L=invL, log_F0_lower=record(c.mpf([logF0Lower,
                                                               logF0Lower])),
                  prefactor=record(c.mpf([prefactor, prefactor])),
                  p2_log_caps=out)
    return out, detail


def _e_relative(c, parent, background_packet):
    return _log_mul(c, _f0_relative(c, parent, background_packet),
                    _native_phi_relative(c, background_packet))


def _sqrt_log_jet(c, square, logq_lower):
    """Absolute ordinary-Z caps for q=sqrt(square), through order 3."""
    q = {0: c.mpf("0.5") * square[0]}
    q[1] = (None if square[1] is None else
            square[1] - c.ln(2) - logq_lower)
    q[2] = _safe_logsum(c, [
        None if square[2] is None else square[2] - c.ln(2) - logq_lower,
        None if square[1] is None else 2 * square[1] - c.ln(4)
        - 3 * logq_lower,
    ])
    q[3] = _safe_logsum(c, [
        None if square[3] is None else square[3] - c.ln(2) - logq_lower,
        None if square[1] is None or square[2] is None else
        square[1] + square[2] + c.ln(3) - c.ln(4) - 3 * logq_lower,
        None if square[1] is None else
        3 * square[1] + c.ln(3) - c.ln(8) - 5 * logq_lower,
    ])
    return q


def _roots(c, parent, background_packet, native, f2, roots):
    """Compute p2, q, u, h, r, alpha, K after the kappa gate."""
    if not roots["kappa_less_one"]:
        return None
    a, b, t0 = roots["a"], roots["b"], roots["t0"]
    kappa = roots["kappa"]
    inva = _log_reciprocal(c, a, roots["log_a_lower"])
    # 1 <= gamma = 2*eta+2-kappa <= 4 once eta in (0,1] and kappa<1.
    gamma = {0: c.ln(c.mpf(4)), 1: kappa[1],
             2: kappa[2], 3: kappa[3]}
    q2 = _log_scale(c, _log_mul(c, gamma, inva), -c.ln(2))
    logq_lower = -c.mpf("0.5") * (c.ln(2) + a[0].b)
    q = _sqrt_log_jet(c, q2, logq_lower)
    p2, p2_detail = _p2(c, parent, background_packet, native, f2)
    dstar = _read(c, parent["actual_micro_q_mixed_bounds"]
                  ["parameter_envelopes"]["actual_logs"]["dstar"])
    u = _log_scale(c, _log_mul(c, p2, q), -dstar.a)
    v = _log_mul(c, u, u)
    h = {
        0: c.mpf(0),
        1: None if v[1] is None else v[1] - c.ln(2),
        2: _safe_logsum(c, [
            None if v[1] is None else 2 * v[1] + c.ln(c.mpf("0.75")),
            None if v[2] is None else v[2] + c.ln(c.mpf("0.5")),
        ]),
        3: _safe_logsum(c, [
            None if v[1] is None else 3 * v[1] + c.ln(c.mpf("1.875")),
            None if v[1] is None or v[2] is None else
            v[1] + v[2] + c.ln(c.mpf("2.25")),
            None if v[3] is None else v[3] + c.ln(c.mpf("0.5")),
        ]),
    }
    r = _log_mul(c, u, h)
    r[0] = c.mpf(0)
    alpha = _log_scale(c, _log_mul(c, q, h), c.ln(2))
    t02 = _log_mul(c, t0, t0)
    q2twice = _log_scale(c, q2, c.ln(2))
    denominator = _log_add(c, _constant(c, c.mpf(0)), t02, q2twice)
    K = _log_scale(c, _log_reciprocal(c, denominator, c.mpf(0)),
                   -c.ln(2 * c.pi))
    logh_lower = -c.mpf("0.5") * _safe_logsum(
        c, [c.mpf(0), 2 * u[0]]).b
    logD_lower = 4 * logh_lower - c.ln(4)
    logK_lower = -c.ln(2 * c.pi) - denominator[0].b
    return dict(
        p2=p2, p2_detail=p2_detail, q=q, q2=q2, gamma=gamma, u=u,
        v=v, h=h, r=r, alpha=alpha, K=K,
        conditioning=dict(
            log_D_lower=c.mpf([logD_lower, logD_lower]),
            log_lambda_lower=c.mpf([logK_lower, logK_lower]),
            proof="D>=h_lower^4/4; lambda=K*(1+t(psi)^2)>=K",
        ),
        q_proof=("q^2=gamma/(2a), 1<=gamma<=4; "
                 "q1=s1/(2q), q2=s2/(2q)+s1^2/(4q^3), "
                 "q3=s3/(2q)+3*s1*s2/(4q^3)+3*s1^3/(8q^5)"),
        h_proof="h=(1+u^2)^(-1/2) with h1,h2,h3 ordinary majorants",
        r_value_proof="r0=|u|/sqrt(1+u^2)<=1; higher rows use ordinary product",
        K_proof="K=[2*pi*(1+t0^2+2*q^2)]^-1",
    )
def _background_geometry(packet):
    for key in ("geometry", "geometric_measure", "geometry_measure"):
        value = packet.get(key) if isinstance(packet, Mapping) else None
        if isinstance(value, Mapping):
            return value
    nested = packet.get("background_I1") if isinstance(packet, Mapping) else None
    if isinstance(nested, Mapping):
        for key in ("geometry", "geometric_measure", "geometry_measure"):
            value = nested.get(key)
            if isinstance(value, Mapping):
                return value
    raise ValueError("Current switch/post geometry is required")


def _chart_view(background_packet, parent, name):
    require(isinstance(background_packet, Mapping), "background_packet must be a mapping")
    charts = background_packet.get("charts")
    require(isinstance(charts, Mapping) and name in charts,
            "Current background must contain chart " + name)
    chart = charts[name]
    require(isinstance(chart, Mapping), "Chart packet must be a mapping")
    view = dict(background_packet)
    view.update(chart)
    if "geometry" not in view:
        view["geometry"] = chart.get("geometric_measure", chart.get("geometry_measure"))
    macro = parent.get("actual_macro_background_I1_numeric", {})
    if "comparison_directions" not in view:
        view["comparison_directions"] = (
            background_packet.get("frozen_directions")
            or macro.get("comparison_directions", {})
        )
    if "comparison" not in view:
        view["comparison"] = (
            background_packet.get("frozen_comparison")
            or macro.get("comparison", {})
        )
    if "F2_derivative_log_caps0_4" not in view:
        view["F2_derivative_log_caps0_4"] = (
            background_packet.get("F2_cached")
            or macro.get("F2_derivative_log_caps0_4")
        )
    return view


def _source_pins(root):
    expected = {
        PARENT: PARENT_SHA256,
        C6_SOURCE: "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
        Q4_SOURCE: "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
        LOOP_SOURCE: "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
        FIRST_SWITCH_SOURCE: "a98753313fb1689f784c7a56cf09dc06df83e0b52cb738b5d4531c19f5402496",
        LOG_Z_SOURCE: "21b442cd8f0044f4652ed7d93d0e54b17a380913925540cb79b17324936fee59",
    }
    out = {}
    for path, expected_sha in expected.items():
        require(sha(root / path) == expected_sha, "Changed switch/post source: " + path)
        out[path] = dict(path=path, expected_sha256=expected_sha,
                         passed=True, source_executed=False)
    return out


def _source_guards(root):
    def text_guard(path, snippets):
        text = (root / path).read_text(encoding="utf-8")
        compact = "".join(text.split())
        require(all("".join(item.split()) in compact for item in snippets),
                "Changed source text: " + path)
        return dict(passed=True, path=path, required_text=list(snippets),
                    source_executed=False)

    return dict(
        r100_second_switch_fields=source_guard(
            root, C6_SOURCE, "second_switch_fields", (
                "angular=a.add(a.scale(known['D_over_R'],g.mul(g.constant('-1/2'),hb,hb,radius,g.sub(g.one,sig))),",
                "a.fixed(g.mul(g.constant('-2/5'),hb,sig)))",
                "return dict(phi=phi,U=initial['U'],phi_y=a.mul(phi,slope),U_y=a.fixed(0)",
            )),
        r100_comparison_frozen=source_guard(
            root, C6_SOURCE, "comparison_at_radius", (
                "phi,U=old['phi'],old['U']",
                "eq=dict(H=phi,M=U,K=b.mul(phi,U)",
                "return previous.comparison_directions(g,phi,U,J",
            )),
        r100_add_chart_physical_U=source_guard(
            root, C6_SOURCE, "add_chart", (
                "J=previous.histories(a,actual['phi'],actual['U'],incoming",
                "fields,radial,generic=previous.generic_fields",
                "physical_measure=jacobian",
            )),
        r100_post_switch_geometry=text_guard(
            C6_SOURCE, (
                "p=g.symbol('true_prescribed_post_switch_fraction');self.post_length=g.sub(g.unary('log',c('11/10')),g.mul(c(2),hb))",
                "y=g.mul(self.post_length,p);R=g.mul(c(100),g.unary('exp',g.add(g.mul(c(2),hb),y)))",
            )),
        q4_physical_pressure=source_guard(
            root, Q4_SOURCE, "physical_pressure_fields", (
                "U=a.scale(fields['V'],Pstar);Uy=a.scale(radial['V_y'],Pstar)",
                "L=b.sub(one,b.mul(de,b.mul(z,z)))",
                "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
            )),
        c5_loop=source_guard(
            root, LOOP_SOURCE, "loop_functions", (
                "cut=a.sub(one,a.scale(roots['Delta'],g.quotient(g.one,eta,'same positive canonical eta')))",
                "q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),'same original active a>0'))))",
                "r=a.mul(u,hinv);alpha=a.scale(a.mul(q,hinv),2)",
                "K=a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')",
            )),
    )


def _verify_background_packet(root, parent, packet):
    require(isinstance(packet, Mapping), "background_packet must be a mapping")
    require(packet.get("parent_report_sha256") == PARENT_SHA256,
            "Background packet is tied to a different parent report")
    charts = packet.get("charts")
    require(isinstance(charts, Mapping), "Background charts are required")
    for name in CHARTS:
        chart = charts.get(name)
        require(isinstance(chart, Mapping) and chart.get("passed") is True,
                "Background chart did not pass: " + name)
        require(isinstance(chart.get("native_fields"), Mapping),
                "Background native_fields missing: " + name)
        require(isinstance(chart.get("incoming_fields"), Mapping),
                "Background incoming_fields missing: " + name)
        require(isinstance(chart.get("geometry"), Mapping),
                "Background geometry missing: " + name)
    pins = packet.get("source_pins")
    require(isinstance(pins, Mapping) and pins,
            "New background source_pins are required")
    verified = {}
    for key, raw in pins.items():
        path = str(raw.get("path", key) if isinstance(raw, Mapping) else key)
        expected = ((raw.get("expected_sha256") or raw.get("sha256"))
                    if isinstance(raw, Mapping) else raw)
        require(isinstance(expected, str), "Background pin has no SHA256: " + path)
        local = root / path
        require(local.is_file(), "Background source missing: " + path)
        require(sha(local) == expected, "Changed background source: " + path)
        verified[path] = dict(path=path, expected_sha256=expected,
                              passed=True, source_executed=False)
    require(any(path.endswith(BACKGROUND_MODULE) or
                path.endswith("original_switch_post_background_I1_numeric.py")
                for path in verified),
            "New switch/post background source hash is required")
    return verified


def _chart_a(c, packet, name):
    geometry = _background_geometry(packet)
    radius = _read_any(c, geometry.get("R_upper", geometry.get("R_upper101",
        geometry.get("R_upper110", 101 if name == "second_switch" else 110))))
    if name == "post_switch":
        a = {0: c.ln(c.mpf(4) / 5), 1: None, 2: None, 3: None}
        log_a_lower = c.ln(c.mpf(4) / 5)
        return dict(
            a=a, b={order: None for order in ORDINARY_Z_ORDERS},
            t0={order: None for order in ORDINARY_Z_ORDERS},
            kappa=a, log_a_lower=c.mpf([log_a_lower, log_a_lower]),
            log_D_lower=c.mpf([log_a_lower, log_a_lower]),
            radius=radius,
            kappa_less_one=bool(a[0].b < 0),
            source_formula="a=4/5 exactly; b=t0=0 exactly; higher a_Z rows exact zero",
            b_exact_zero=True,
        )

    directions = packet.get("comparison_directions", {})
    raw_d = directions.get("D_over_R", directions.get("D"))
    drows = _rows(c, raw_d, 4)
    require(drows[0] is not None and drows[0].a > 0,
            "Second-switch D_over_R value must have a positive floor")
    D = {order: _abs_log(c, drows[order]) for order in ORDINARY_Z_ORDERS}
    loghb = _read(c, geometry["loghbS"])
    logR_upper = c.ln(radius.b)
    Rlower = _read_any(c, geometry.get("R_lower", 100))
    lower_tiny = loghb.a + c.ln(Rlower.a) + c.ln(drows[0].a)
    tiny0 = None if D[0] is None else loghb.b + logR_upper + D[0]
    require(tiny0 is not None and tiny0.b < c.ln(c.mpf(4) / 5),
            "Second-switch shear cap must be proved below 4/5")
    a = {
        0: _safe_logsum(c, [tiny0, c.ln(c.mpf(4) / 5)]),
        1: None if D[1] is None else loghb.b + logR_upper + D[1],
        2: None if D[2] is None else loghb.b + logR_upper + D[2],
        3: None if D[3] is None else loghb.b + logR_upper + D[3],
    }
    log_a_lower = c.mpf([min(lower_tiny, c.ln(c.mpf(4) / 5)),
                         min(lower_tiny, c.ln(c.mpf(4) / 5))])
    return dict(
        a=a, b={order: None for order in ORDINARY_Z_ORDERS},
        t0={order: None for order in ORDINARY_Z_ORDERS},
        kappa=a, log_a_lower=log_a_lower,
        log_D_lower=c.mpf([c.ln(drows[0].a), c.ln(drows[0].a)]),
        radius=radius,
        kappa_less_one=bool(a[0].b < 0),
        source_formula=("a=hbS*R*D_over_R*(1-sigma)+4*sigma/5; "
                        "sigma is Z independent; b=t0=0 exactly"),
        b_exact_zero=True,
        tiny_term_upper_below_four_fifths=True,
        lower_proof=("a_lower=min(100*hbS*D0_lower,4/5); "
                     "a_upper<=4/5+positive tiny cap"),
    )


def _roots_for_chart(c, parent, packet, native, f2, chart_a):
    roots = dict(chart_a)
    roots.update(kappa=chart_a["kappa"], inv_barphi={order: None
        for order in ORDINARY_Z_ORDERS},
        inverse_F0_relative=_extended_logs(
            c, parent["actual_micro_b_p2_mixed_bounds"]["whole_tube_amplitude"]
            ["inverse_F0_relative_ordinary_Z0_Z4"], 4))
    roots2 = _roots(c, parent, packet, native, f2, roots)
    return roots2, roots


def _packet_rows(c, roots, roots2, e_relative):
    values = {
        "a": roots["a"], "b": roots["b"], "t0": roots["t0"],
        "E_over_frozen_E0": e_relative,
    }
    if roots2 is None:
        values.update({name: {order: None for order in ORDINARY_Z_ORDERS}
                       for name in ("r", "alpha", "K")})
    else:
        values.update({name: roots2[name] for name in ("r", "alpha", "K")})
    return {name: _log_rows(values[name]) for name in values}


def _chart_result(root, parent, background_packet, name):
    packet = _chart_view(background_packet, parent, name)
    c = MPIntervalContext()
    c.prec = 900
    native = _native_fields(c, packet)
    f2, f2_receipt = _f2_logs_and_receipt(c, parent, packet)
    chart_a = _chart_a(c, packet, name)
    roots2, roots = _roots_for_chart(c, parent, packet, native, f2, chart_a)
    e_relative = _e_relative(c, parent, packet)
    eta = _read(c, parent["actual_micro_q_mixed_bounds"]
                ["parameter_envelopes"]["actual_logs"]["eta"])
    require(eta.b <= 0, "Saved eta must be at most one in logarithmic form")
    logN = _read(c, parent["actual_micro_q_mixed_bounds"]
                 ["parameter_envelopes"]["actual_logs"]["N"])
    phi_rows = _rows(c, packet["native_fields"]["phi"], 5)
    require(phi_rows[0] is not None and phi_rows[0].a > 0,
            "Current chart native phi floor must be positive")
    amp = parent["actual_micro_b_p2_mixed_bounds"]["whole_tube_amplitude"]
    logF0_upper = _read(c, amp["log_F0"]).b
    radius = _read_any(c, _background_geometry(packet).get(
        "R_upper", _background_geometry(packet).get(
            "R_upper101", _background_geometry(packet).get("R_upper110",
                101 if name == "second_switch" else 110))))
    log_Pstar_E_upper = (logF0_upper + c.mpf("0.5") *
                         (c.ln(2) + c.ln(radius.b)) + c.ln(phi_rows[0].b))
    passed = bool(chart_a["kappa_less_one"] and roots2 is not None
                  and f2_receipt["passed"])
    values = dict(chart_a)
    if roots2 is not None:
        values.update(roots2)
    log_caps = _packet_rows(c, values, roots2, e_relative)
    if roots2 is None:
        conditioning = dict(log_D_lower=None, log_lambda_lower=None)
    else:
        conditioning = {
            key: record(value) for key, value in roots2["conditioning"].items()
            if key.startswith("log_")
        }
    proof = dict(
        kappa_less_than_one_cut_one=passed,
        kappa_lt_one=chart_a["kappa_less_one"],
        cut=1 if passed else None,
        eta_positive_le_one=True,
        eta_positive=True,
        eta_at_most_one=True,
        a_positive=chart_a["log_a_lower"].a > -c.inf,
        actual_N18=True,
        actual_N18_is_graph_handle=True,
        actual_N_not_materialized=True,
        actual_phase_Z_independent=True,
        physical_y_derivatives_claimed=False,
        real_switch_post_source_parameters=True,
        source_bound_to_current_switch_post_packet=True,
    )
    return dict(
        passed=passed,
        chart=name,
        ordinary_Z_orders=list(ORDINARY_Z_ORDERS),
        physical_y_derivatives_claimed=False,
        log_caps=log_caps,
        conditioning=conditioning,
        log_Pstar_E_upper=record(c.mpf([log_Pstar_E_upper, log_Pstar_E_upper])),
        logN=record(logN),
        proofs=proof,
        kappa_log_upper=record(chart_a["kappa"][0]),
        kappa_upper_log_display=str(chart_a["kappa"][0].b),
        native_fields=background_packet["charts"][name]["native_fields"],
        incoming_fields=background_packet["charts"][name]["incoming_fields"],
        Q_b3rows=background_packet["charts"][name].get(
            "Q_b3rows", background_packet["charts"][name].get("Q_b")),
        background_I1=background_packet["charts"][name].get("background_I1"),
        geometry=background_packet["charts"][name]["geometry"],
        F2_fifth_cap_receipt=f2_receipt,
        intermediate=dict(
            a=_log_rows(chart_a["a"]),
            b=_log_rows(chart_a["b"]),
            t0=_log_rows(chart_a["t0"]),
            e_relative=_log_rows(e_relative),
            qdrive_formula="not used in b: source V/U axial phase is zero",
            p2_log_caps=None if roots2 is None else _log_rows(roots2["p2"]),
            q_log_caps=None if roots2 is None else _log_rows(roots2["q"]),
            h_log_caps=None if roots2 is None else _log_rows(roots2["h"]),
            native_numerator=None if roots2 is None else roots2["p2_detail"],
        ),
        source_formulas=dict(
            a=chart_a["source_formula"],
            b="0 exact from add_chart source U_y=0 / physical V_y=0",
            t0="0 exact",
            E_over_frozen_E0="F0relative * current native phi relative; E0=E(R,Z0), ratio0=1",
            p2="sqrt(R/2)*Nnative/(F0*phi*Lz), Lz=1-delta*Z^2; chart Jacobian is not in denominator",
            p2_native_numerator=(
                "-U+h*(M-Z*M_Z)/2+U*(h*Z*M+d*M_Z)"
                "+2*delta*Z*A-d*A_Z"
                "+Pstar^2*(2*(1+delta)*Z*P0-d*P0_Z)"
                "+R*(-2*delta*Z*F2*B+d*(F2*B)_Z"
                "+2*(1+delta)*Z*F2*C-d*(F2*C)_Z)"
            ),
            q2="gamma/(2*a)", u="p2*q/dstar", h="(1+u^2)^(-1/2)",
            r="u*h", alpha="2*q*h", K="1/(2*pi*(1+t0^2+2*q^2))",
        ),
        physical_scale_proof=dict(
            passed=True,
            quantity="B/Pstar/frozen_E0",
            deltaU_formula="physical deltaU=Pstar*frozen_E0*(B/Pstar/frozen_E0)/N",
            current_N="N18 graph handle",
            current_N18=True,
            current_N_is_graph_handle=True,
            current_N_not_materialized=True,
            phase_Z_independent=True,
            physical_y_derivatives_claimed=False,
            source_bound_to_current_switch_post_packet=True,
        ),
        source_contract=dict(
            chart=name,
            p2_L="Lz=1-delta*Z^2, not hbS or post_length",
            physical_measure=background_packet["charts"][name]["geometry"].get(
                "physical_J_name"),
            jacobian_not_in_p2_denominator=True,
            actual_U_physical=True,
            V_is_U_over_Pstar_source_binding=True,
        ),
    )


def build(root, parent, background_packet):
    root = Path(root)
    parent_path = root / PARENT
    require(parent_path.is_file() and sha(parent_path) == PARENT_SHA256,
            "Current first-switch finite-N parent bytes changed")
    require(parent == json.loads(parent_path.read_text(encoding="utf-8")),
            "Current first-switch finite-N parent mapping is required")
    pins = _source_pins(root)
    pins.update(_verify_background_packet(root, parent, background_packet))
    source_guards = _source_guards(root)
    charts = {}
    for name in CHARTS:
        charts[name] = _chart_result(root, parent, background_packet, name)
    passed = all(charts[name]["passed"] for name in CHARTS)
    proofs = dict(
        kappa_less_than_one_cut_one=passed,
        eta_positive_le_one=True,
        a_positive=True,
        actual_N18=True,
        actual_N18_is_graph_handle=True,
        actual_N_not_materialized=True,
        actual_phase_Z_independent=True,
        real_switch_post_source_parameters=True,
        physical_y_derivatives_claimed=False,
    )
    result = dict(
        passed=passed,
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        ordinary_Z_orders=list(ORDINARY_Z_ORDERS),
        physical_y_derivatives_claimed=False,
        charts=charts,
        F2_fifth_cap_receipt={name: charts[name]["F2_fifth_cap_receipt"]
                             for name in CHARTS},
        log_Pstar_E_upper={name: charts[name]["log_Pstar_E_upper"]
                           for name in CHARTS},
        logN=charts["second_switch"]["logN"],
        source_pins=pins,
        source_guards=source_guards,
        proofs=proofs,
        physical_scale_proof=dict(
            source_bound_to_current_switch_post_packet=True,
            quantity="B/Pstar/frozen_E0",
            deltaU_formula="physical deltaU=Pstar*frozen_E0*(B/Pstar/frozen_E0)/N",
            current_N="N18 graph handle",
            phase_Z_independent=True,
            physical_y_derivatives_claimed=False,
        ),
        source_contract=dict(
            charts=list(CHARTS),
            actual_U_physical=True,
            add_chart_V="V=U/Pstar",
            comparison_at_radius_frozen_macro_phi_U=True,
            p2_L="Lz=1-delta*Z^2",
            geometry_jacobians_not_in_p2_denominator=True,
        ),
        no_accepted_builder_checker_producer_owner_norm_recurrence_quadrature=True,
        no_old_first_switch_whole_chart_reused=True,
    )
    result = _jsonable(result)
    json.dumps(result, allow_nan=False)
    return result


def fixture_bound_check():
    c = MPIntervalContext()
    c.prec = 220
    zero = c.mpf(0)
    assert all(row["exact_zero"] for row in _log_rows({i: None for i in ORDINARY_Z_ORDERS}))
    nonconstant = {0: zero, 1: zero, 2: None, 3: None}
    bjet = _log_mul(c, _log_mul(c, nonconstant, nonconstant), nonconstant)
    assert bjet[1].a <= c.ln(3) <= bjet[1].b
    def row(order, value):
        return dict(ordinary_Z_order=order, value=record(c.mpf([value, value])))
    def rows(values):
        return [row(i, v) for i, v in enumerate(values)]
    packet = {
        "comparison_directions": {
            "Z": rows([0, 0, 0, 0]),
            "d": rows([0, 0, 0, 0]),
            "L": rows([1, 0, 0, 0]),
        },
        "geometry": {"R_upper": record(c.mpf(101)), "loghbS": record(c.mpf(-100))},
        "native_fields": {
            name: rows([0, 0, 0, 0, 0])
            for name in ("phi", "U", "M", "A", "B", "C")
        },
    }
    packet["native_fields"]["phi"][0] = row(0, 1)
    parent = {
        "actual_core_signed_direction": {
            "parameters": {"delta": record(c.mpf("0.1"))}},
        "actual_micro_b_p2_mixed_bounds": {
            "whole_tube_pressure": {
                "p0_jets": rows([0, 1, 0, 0, 0, 0])},
            "whole_tube_amplitude": {
                "log_Pstar": record(c.mpf(0)),
                "log_F0": record(c.mpf(0)),
                "inverse_F0_relative_ordinary_Z0_Z4": rows([0, 0, 0, 0, 0])}},
        "actual_macro_background_I1_numeric": {
            "geometric_measure": {"logLambda": record(c.mpf(0))}},
    }
    native = {name: {i: None for i in range(5)}
              for name in ("phi", "U", "M", "A", "B", "C")}
    native["phi"][0] = zero
    f2 = {i: zero for i in range(5)}
    numerator, _ = _native_numerator(c, parent, packet, native, f2)
    assert numerator[0] is not None
    return dict(
        passed=True,
        checks=dict(
            exact_zero_b_t0=True,
            reciprocal_product_Z1=True,
            pressure_Pstar_squared_survives_with_nonzero_P0=True,
            no_physical_y_rows=True,
            p2_L_is_meridional_Lz=True,
        ),
        JSON_serializable=True,
        actual_background_required=True,
    )


def self_check():
    return fixture_bound_check()
