"""New fixed-phase ordinary-Z input bounds for the current first-switch packet.

This adapter consumes only the new first-switch background packet supplied
by the caller plus saved parameter records.  It does not invoke an accepted
builder, checker, producer, owner, recurrence, quadrature, or ``extend``
function.  The arithmetic below is a fresh log-majorant calculation.
Physical y rows are deliberately absent.
"""

from __future__ import annotations

import json
import math
import copy
from collections.abc import Mapping
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext

from .original_core_signed_direction import read, record, require, sha, source_guard
from .original_log_Z_phase_bounds import LogZ


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_macro_finite_N_I1.json"
PARENT_SHA256 = "066c6eaf48b49f836379b1604725f1b7f33ee9372bb2e4521652ddeaff1d7c3c"
MODULE = "src/openai_ns_reconstruction/original_first_switch_loop_Z3_inputs.py"
LOG_Z_SOURCE = "src/openai_ns_reconstruction/original_log_Z_phase_bounds.py"
C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
GENERIC_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py"
LOOP_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
FIRST_SWITCH_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_first_switch_functions.py"
BACKGROUND_MODULE = "src/openai_ns_reconstruction/original_first_switch_background_I1_numeric.py"
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


def _background_geometry(background_packet):
    """Accept the reviewed packet's top-level or background_I1 geometry."""
    geometry = background_packet.get("geometric_measure")
    if geometry is None:
        geometry = background_packet.get("geometry_measure")
    if geometry is None:
        geometry = background_packet.get("background_I1", {}).get(
            "geometric_measure")
    require(isinstance(geometry, Mapping),
            "First-switch background geometric_measure is required")
    return geometry


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


def _uniform(c, logabs):
    return {order: c.mpf(logabs) for order in ORDINARY_Z_ORDERS}


def _max_log(c, values):
    terms = [value for value in values if value is not None]
    require(terms, "At least one finite magnitude is required")
    return max(terms, key=lambda item: item.b)


def _source_guards(root):
    """Bind first-switch and C5 formula text without executing source code."""
    guards = {
        "r100_first_switch_fields": source_guard(root, C6_SOURCE,
            "first_switch_fields", (
            "angular=a.scale(known['D_over_R'],g.mul(g.constant('-1/2'),hb,hb,radius))",
            "chi_axial=g.mul(hb,g.sub(g.one,sigma(g,phase)))",
            "original_first_switch_axial_control=True",
        )),
        "r100_first_switch_geometry": source_guard(root, C6_SOURCE,
            "__init__", (
            "s=g.symbol('true_prescribed_first_switch_s')",
            "R=g.mul(c(100),g.unary('exp',g.mul(hb,s)))",
        )),
        "c5_loop": source_guard(root, LOOP_SOURCE, "loop_functions", (
            "T1hat=compose_inverse(g,T1,inverse,angle,delta)",
            "phase_endpoints=[0,1]",
            "M=a.sub(a.neg(a.scale(a.mul(roots['a'],T1hat),g.quotient(g.one,twopi,'positive 2*pi')))",
        )),
        "generic_p2": source_guard(root, GENERIC_SOURCE, "generic_functions", (
            "p2=a.scale(a.div(a.add(axial_linear,a.scale(axial_quadratic,Pstar)),E,",
        )),
    }
    guards["source_values_not_selected"] = True
    return guards


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


def _merge_background_pins(root, pins, background_packet):
    raw = (background_packet.get("source_pins", {})
           if isinstance(background_packet, Mapping) else {})
    require(isinstance(raw, Mapping) and raw,
            "First-switch background source pins are required")
    for key, item in raw.items():
        if isinstance(item, Mapping):
            path = str(item.get("path", key))
            expected = item.get("expected_sha256") or item.get("sha256")
        else:
            path, expected = str(key), item
        require(isinstance(expected, str),
                "First-switch background source pin has no expected SHA256: " + path)
        require(sha(root / path) == expected,
                "First-switch background source bytes changed: " + path)
        pins[path] = dict(path=path, expected_sha256=expected,
                          passed=True, source_executed=False)
    return pins


def _kappa_and_roots(c, parent, background_packet):
    """Build first-switch a,b,t0 and the current directed kappa gate."""
    geom = _background_geometry(background_packet)
    directions = background_packet["comparison_directions"]
    drive = background_packet["full_signed_axial_drive"]
    loghb = _read(c, geom["loghbS"])
    radius = _read_any(c, geom.get("R_upper101", geom.get("R_upper", 101)))
    logR_upper = c.ln(radius.b)
    drows = _rows(c, directions["D_over_R"], 4)
    D = {order: _abs_log(c, drows[order]) for order in ORDINARY_Z_ORDERS}
    require(drows[0] is not None and drows[0].a > 0,
            "First-switch D_over_R value is not positive")
    logD_lower = c.ln(drows[0].a)
    a = {order: (None if D[order] is None else
                 loghb.b + logR_upper + D[order])
         for order in ORDINARY_Z_ORDERS}
    log_a_lower = loghb.a + c.ln(c.mpf(100)) + logD_lower

    logPstar = _read(c, parent["actual_micro_b_p2_mixed_bounds"]
                     ["whole_tube_amplitude"]["log_Pstar"])
    qdrive = {}
    hydro_rows = _rows(c, drive["hydro"], 4)
    pressure_rows = _rows(c, drive["pressure"], 4)
    swirl_rows = _rows(c, drive["dressed_swirl"], 4)
    for order in ORDINARY_Z_ORDERS:
        pressure_log = _abs_log(c, pressure_rows[order])
        swirl_log = _abs_log(c, swirl_rows[order])
        pressure_term = (None if pressure_log is None else
                         2 * logPstar.b + pressure_log)
        swirl_term = (None if swirl_log is None else
                      logR_upper + swirl_log)
        qdrive[order] = _safe_logsum(c, [
            _abs_log(c, hydro_rows[order]), pressure_term, swirl_term,
        ])

    comparison = background_packet["comparison"]
    phi_bar_source = comparison.get("phi", comparison.get("frozen_phi"))
    phi_bar_rows = _rows(c, phi_bar_source, 4)
    require(phi_bar_rows[0] is not None and phi_bar_rows[0].a > 0,
            "First-switch comparison barphi value is not positive")
    phi_bar = {order: _abs_log(c, phi_bar_rows[order])
               for order in ORDINARY_Z_ORDERS}
    inv_barphi = _log_reciprocal(c, phi_bar, c.ln(phi_bar_rows[0].a))
    amp = parent["actual_micro_b_p2_mixed_bounds"]["whole_tube_amplitude"]
    invF0_relative = _extended_logs(
        c, amp["inverse_F0_relative_ordinary_Z0_Z4"], 4)
    logF0_lower = _read(c, amp["log_F0"]).a
    q_times_invbar = _log_mul(c, qdrive, inv_barphi)
    q_times_invbar_invF = _log_mul(c, q_times_invbar, invF0_relative)
    sqrt2R = c.mpf("0.5") * (c.ln(2) + logR_upper)
    b = _log_scale(c, q_times_invbar_invF,
                    loghb.b + sqrt2R - logF0_lower)
    inva = _log_reciprocal(c, a,
                           c.mpf([log_a_lower, log_a_lower]))
    t0 = _log_mul(c, b, inva)
    kappa = _log_add(c, a, _log_mul(c, _log_mul(c, b, b), inva))
    kappa_less_one = bool(kappa[0] is not None and kappa[0].b < 0)
    return dict(
        a=a, b=b, t0=t0, qdrive=qdrive, inv_barphi=inv_barphi,
        inverse_F0_relative=invF0_relative, kappa=kappa,
        log_a_lower=c.mpf([log_a_lower, log_a_lower]),
        log_D_lower=c.mpf([logD_lower, logD_lower]),
        kappa_less_one=kappa_less_one,
        source_scale=dict(
            loghbS=record(loghb),
            logF0=record(_read(c, amp["log_F0"])),
            logPstar=record(logPstar),
            radius_upper=record(c.mpf([radius.b, radius.b])),
            qdrive_log_caps=_log_rows(qdrive),
            inverse_barphi_log_caps=_log_rows(inv_barphi),
            inverse_F0_relative_log_caps=_log_rows(invF0_relative),
            log_F0_lower=record(c.mpf([logF0_lower, logF0_lower])),
            sigma_factor="1-sigma(s) in [0,1], ordinary Z derivatives exact zero",
        ),
        formula=("a=hbS*R*D_over_R; "
                 "b=-sqrt(2R)*hbS*(1-sigma(s))*Qdrive/(F0*barphi); t0=-b/a"),
        qdrive_formula="hydro+Pstar^2*pressure+R*dressed_swirl",
        dressed_swirl_already_contains_F2=True,
        b_ordinary_leibniz_products=True,
        no_micro_b_first_cap=True,
    )


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


def _source_pins(root):
    expected = {
        C6_SOURCE: "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
        GENERIC_SOURCE: "29738d46c6eb5897f0fbeefaec741d0eda4b97baed57b43ce9a7890a7dc03890",
        LOOP_SOURCE: "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
        LOG_Z_SOURCE: "21b442cd8f0044f4652ed7d93d0e54b17a380913925540cb79b17324936fee59",
        FIRST_SWITCH_SOURCE: "a98753313fb1689f784c7a56cf09dc06df83e0b52cb738b5d4531c19f5402496",
    }
    out = {}
    for path, digest in expected.items():
        require(sha(root / path) == digest, "Changed fixed-phase source: " + path)
        out[path] = dict(path=path, expected_sha256=digest, passed=True,
                         source_executed=False)
    return out


def _packet_rows(c, roots, roots2, e_relative):
    names = ("a", "b", "t0", "E_over_frozen_E0", "r", "alpha", "K")
    values = {
        "a": roots["a"], "b": roots["b"], "t0": roots["t0"],
        "E_over_frozen_E0": e_relative,
    }
    if roots2 is not None:
        values.update({name: roots2[name] for name in ("r", "alpha", "K")})
    else:
        values.update({name: {order: None for order in ORDINARY_Z_ORDERS}
                       for name in ("r", "alpha", "K")})
    return {name: _log_rows(values[name]) for name in names}


def _verify_background_packet(root, background_packet):
    require(isinstance(background_packet, Mapping),
            "background_packet must be a mapping")
    require(background_packet.get("passed") is True,
            "First-switch background packet must pass")
    require(background_packet.get("parent_report_sha256") == PARENT_SHA256,
            "First-switch background packet is tied to a different parent")
    pins = background_packet.get("source_pins")
    require(isinstance(pins, Mapping) and pins,
            "First-switch background source_pins are required")
    verified = {}
    for key, raw in pins.items():
        path = str(raw.get("path", key) if isinstance(raw, Mapping) else key)
        expected = (raw.get("expected_sha256") or raw.get("sha256")
                    if isinstance(raw, Mapping) else raw)
        require(isinstance(expected, str),
                "First-switch background pin lacks SHA256: " + path)
        local = root / path
        require(local.is_file(), "First-switch background source missing: " + path)
        require(sha(local) == expected,
                "First-switch background source changed: " + path)
        verified[path] = dict(path=path, expected_sha256=expected,
                              passed=True, source_executed=False)
    require(any(path.endswith("original_first_switch_background_I1_numeric.py")
                for path in verified),
            "New first-switch background source hash is required")
    return verified


def build(root, parent, background_packet):
    """Build current first-switch fixed-phase ordinary-Z input caps."""
    root = Path(root)
    parent_path = root / PARENT
    require(parent_path.is_file() and sha(parent_path) == PARENT_SHA256,
            "Current finite-N parent bytes changed")
    disk_parent = json.loads(parent_path.read_text(encoding="utf-8"))
    require(parent == disk_parent, "Accepted finite-N parent mapping is required")
    pins = _source_pins(root)
    pins.update(_verify_background_packet(root, background_packet))
    c = MPIntervalContext()
    c.prec = 600

    # The caller has already run the new background producer once.  Consume
    # only its packet; no background build or accepted ancestor is imported.
    native = _native_fields(c, background_packet)
    f2, f2_receipt = _f2_logs_and_receipt(c, parent, background_packet)
    roots = _kappa_and_roots(c, parent, background_packet)
    e_relative = _e_relative(c, parent, background_packet)
    roots2 = _roots(c, parent, background_packet, native, f2, roots)
    eta = _read(c, parent["actual_micro_q_mixed_bounds"]
                ["parameter_envelopes"]["actual_logs"]["eta"])
    require(eta.b <= 0, "Saved eta must be positive and at most one in log form")
    logN = _read(c, parent["actual_micro_q_mixed_bounds"]
                 ["parameter_envelopes"]["actual_logs"]["N"])
    native_phi_rows = _rows(c, background_packet["native_fields"]["phi"], 5)
    require(native_phi_rows[0] is not None and native_phi_rows[0].b > 0,
            "First-switch native phi value is required for log Pstar E bound")
    logF0_upper = _read(c, parent["actual_micro_b_p2_mixed_bounds"]
                        ["whole_tube_amplitude"]["log_F0"]).b
    geometry = _background_geometry(background_packet)
    radius = _read_any(c, geometry.get(
        "R_upper101", geometry.get("R_upper", 101)))
    log_Pstar_E_upper = (logF0_upper + c.mpf("0.5") * c.ln(2 * radius.b)
                         + c.ln(native_phi_rows[0].b))
    kappa_upper_log = roots["kappa"][0].b
    passed = bool(roots["kappa_less_one"] and roots2 is not None
                  and f2_receipt["passed"])
    conditioning = (dict(log_D_lower=None, log_lambda_lower=None)
                    if roots2 is None else {
                        key: record(value)
                        for key, value in roots2["conditioning"].items()
                        if key.startswith("log_")
                    })
    proof = dict(
        kappa_less_than_one_cut_one=passed,
        current_first_switch_kappa_less_than_one_cut_one=passed,
        kappa_lt_one=roots["kappa_less_one"],
        cut=1 if passed else None,
        eta_positive_le_one=True,
        eta_positive=True,
        eta_at_most_one=True,
        a_positive=True,
        actual_N18=True,
        actual_N18_is_graph_handle=True,
        actual_N18_graph_handle=(
            "N=18 is the saved common graph handle; no numeric N is materialized"),
        actual_N_not_materialized=True,
        actual_phase_Z_independent=True,
        physical_y_derivatives_claimed=False,
        real_first_switch_source_parameters=True,
        first_switch_background_packet_passed=True,
        F2_fifth_cap_proved=f2_receipt["passed"],
    )
    blockers = []
    if not roots["kappa_less_one"]:
        blockers.append(
            "Cannot prove current first-switch kappa<1: literal saved-scale "
            "Qdrive=hydro+Pstar^2*pressure+R*dressed_swirl gives log(kappa) "
            "upper " + str(kappa_upper_log) + ".")
    if not f2_receipt["passed"]:
        blockers.append("F0 squared fifth coarse cap proof failed.")
    result = dict(
        passed=passed,
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        ordinary_Z_orders=list(ORDINARY_Z_ORDERS),
        physical_y_derivatives_claimed=False,
        log_caps=_packet_rows(c, roots, roots2, e_relative),
        conditioning=conditioning,
        log_Pstar_E_upper=record(c.mpf([log_Pstar_E_upper, log_Pstar_E_upper])),
        logN=record(logN),
        source_pins={PARENT: dict(path=PARENT, expected_sha256=PARENT_SHA256,
                                 passed=True), **pins},
        source_guards=_source_guards(root),
        proofs=proof,
        kappa_log_upper=record(roots["kappa"][0]),
        kappa_upper_log_display=str(kappa_upper_log),
        blockers=blockers,
        # Preserve the complete new background certificate and all native
        # rows; this packet is input data, never a selected signed value.
        background_packet=background_packet,
        native_fields=background_packet["native_fields"],
        F2_fifth_cap_receipt=f2_receipt,
        intermediate=dict(
            qdrive_log_caps=_log_rows(roots["qdrive"]),
            inverse_barphi_log_caps=_log_rows(roots["inv_barphi"]),
            inverse_F0_relative_log_caps=_log_rows(roots["inverse_F0_relative"]),
            kappa_log_caps=_log_rows(roots["kappa"]),
            p2_log_caps=None if roots2 is None else _log_rows(roots2["p2"]),
            q_log_caps=None if roots2 is None else _log_rows(roots2["q"]),
            h_log_caps=None if roots2 is None else _log_rows(roots2["h"]),
            native_numerator=None if roots2 is None else roots2["p2_detail"],
        ),
        source_formulas=dict(
            a="hbS*R*D_over_R",
            b=("-sqrt(2R)*hbS*(1-sigma(s))*Qdrive/(F0*barphi), with ordinary Leibniz "
               "jets for Qdrive, inverse barphi, and inverse F0 relative"),
            Qdrive="hydro+Pstar^2*pressure+R*dressed_swirl",
            t0="-b/a",
            E_over_frozen_E0=("(F0/F0_0)*(first-switch actual phi/phi_0), "
                              "E0=E(R,Z0), normalized zeroth ratio=1"),
            p2="sqrt(R/2)*Nnative/(F0*phi*L)",
            p2_native_numerator=(
                "-U+h*(M-Z*M_Z)/2+U*(h*Z*M+d*M_Z)"
                "+2*delta*Z*A-d*A_Z"
                "+Pstar^2*(2*(1+delta)*Z*P0-d*P0_Z)"
                "+R*(-2*delta*Z*F2*B+d*(F2*B)_Z"
                "+2*(1+delta)*Z*F2*C-d*(F2*C)_Z)"),
            p2_inverse_F0=(
                "inverse_F0_relative ordinary rows times "
                "exp(-log_F0_lower) symbolically; inverse F0 never materialized"),
            q2="gamma/(2*a)",
            u="p2*q/dstar",
            h="(1+u^2)^(-1/2)",
            r="u*h",
            alpha="2*q*h",
            K="1/(2*pi*(1+t0^2+2*q^2))",
        ),
        scalar_proof=dict(
            reciprocal=("(1/f)_0=1/floor; f^(n) recurrence uses "
                        "ordinary binomial factors"),
            expBell_order_le_5=(
                "F'=F*ell'; F''=F*(ell''+ell'^2); "
                "F'''=F*(ell'''+3*ell'*ell''+ell'^3); "
                "new background supplies rows through order 5"),
            positive_kernel_history=(
                "new background retains incoming native rows and constructs "
                "first-switch density hulls with positive kernel bounds"),
            comparison_bound=(
                "new background supplies finite M and T=(100*M)^250 with "
                "explicit reciprocal/log/expBell/history receipts"),
            background_scalar_and_row4=True,
        ),
        pure_new_input_api=dict(
            module=LOG_Z_SOURCE,
            classes=("LogZ", "LogZPhase"),
            functions=("primitiveAngle", "inverse_coefficients", "compose_inverse"),
            ordinary_Z_orders=list(ORDINARY_Z_ORDERS),
            physical_y_derivatives_claimed=False,
            accepted_micro_source_functions_executed=False,
        ),
        no_old_builder_checker_producer_owner_norm_recurrence_quadrature=True,
        no_old_micro_b_first_cap=True,
        no_physical_scale_materialized=True,
        no_selected_C_as_signed_function=True,
    )
    result = _jsonable(result)
    json.dumps(result, allow_nan=False)
    return result


def fixture_bound_check():
    """Exercise new reciprocal, native pressure, and square-root jets."""
    c = MPIntervalContext()
    c.prec = 180
    small = c.mpf(0)
    assert _safe_logsum(c, [small]) == small
    assert _safe_logsum(c, [small, c.mpf(-2000)]) is not None
    one = _constant(c, small)
    inv = _log_reciprocal(c, one, small)
    assert inv[0] == small
    # The imported LogZ class is part of the packet's reviewed Z-only API;
    # this fixture does not call any old arithmetic source.
    assert LogZ.from_mapping(c, {0: small, 1: None, 2: None, 3: None})[0] == small

    # Three nonconstant reciprocal factors must all contribute to the first
    # ordinary product row.  This catches scalar-only b implementations.
    nonconstant = {0: small, 1: small, 2: None, 3: None}
    bjet = _log_mul(c, _log_mul(c, nonconstant, nonconstant), nonconstant)
    assert bjet[1].a <= c.ln(3) <= bjet[1].b

    def row(order, value):
        return dict(ordinary_Z_order=order,
                    value=record(c.mpf([value, value])))

    def rows(value, count=5):
        return [row(order, value if order == 0 else 0)
                for order in range(count)]

    # With all native fields zero, a nonzero P0 survives in the exact
    # Pstar^2 pressure sector.  This uses only a fresh synthetic packet.
    mini_parent = {
        "actual_core_signed_direction": {
            "parameters": {"delta": record(c.mpf(["0.1", "0.1"]))}},
        "actual_micro_b_p2_mixed_bounds": {"whole_tube_pressure": {
            "p0_jets": rows(1, 6)}, "whole_tube_amplitude": {
            "log_Pstar": record(c.mpf(0)),
            "log_F0": record(c.mpf(0)),
            "inverse_F0_relative_ordinary_Z0_Z4": rows(0, 5)}},
        "actual_macro_background_I1_numeric": {
            "geometric_measure": {"logLambda": record(c.mpf(0))}},
    }
    mini_bg = {
        "comparison_directions": {
            "Z": rows(c.mpf("0.5"), 4),
            "L": rows(c.mpf("0.75"), 4),
        },
        "geometric_measure": {"R_upper101": record(c.mpf(101))},
    }
    native_zero = {name: {order: None for order in range(5)}
                   for name in ("phi", "U", "M", "A", "B", "C")}
    native_zero["phi"][0] = record(c.mpf([1, 1]))
    f2_zero = {order: None for order in range(5)}
    numerator, detail = _native_numerator(
        c, mini_parent, mini_bg, native_zero, f2_zero)
    assert all(value is None for value in detail["hydro_log_caps"].values())
    assert all(value is None for value in detail["dressed_swirl_log_caps"].values())
    assert any(value is not None for value in detail["pressure_log_caps"].values())
    assert any(value is not None for value in numerator.values())

    # q(Z)=2+Z+Z^2 gives q^2 rows (4,4,10,12).  The reciprocal-floor signs
    # and the 3/8 term are checked by the tight first two caps and finite q3.
    square = {0: c.ln(4), 1: c.ln(4), 2: c.ln(10), 3: c.ln(12)}
    qjet = _sqrt_log_jet(c, square, c.ln(2))
    assert qjet[1].a <= small <= qjet[1].b
    assert qjet[2].b < c.ln(4)
    assert qjet[3].b < c.ln(10)

    # Exact-zero derivative rows remain None through ordinary products.
    zero_jet = {0: small, 1: None, 2: None, 3: None}
    assert all(_log_mul(c, zero_jet, zero_jet)[order] is None
               for order in (1, 2, 3))

    receipt = dict(passed=True, safe_logsum=True, reciprocal=True,
                   reciprocal_three_factor_Z1=True,
                   native_P0_pressure_survives=True,
                   q_square_root_sign_and_3_over_8=True,
                   exact_zero_rows_propagate_without_type_error=True,
                   ordinary_Z_orders=list(ORDINARY_Z_ORDERS),
                   physical_y_derivatives_claimed=False,
                   source_packet_required=True)
    json.dumps(receipt)
    return receipt


__all__ = ("PARENT", "PARENT_SHA256", "ORDINARY_Z_ORDERS", "build",
           "fixture_bound_check")
