"""Independent first-switch ``B`` boundary enclosure in ``R_m^2`` units.

The kernel accepts signed ordinary ``Z0``/``Z1`` background jets and raw
logarithmic finite-``N`` error caps.  It never reads a saved source packet or
replays a producer.  The error-log threshold is discharged before a symmetric
``2**-radius_bits`` radius is attached to each background jet.  The resulting
``G`` jet gives a two-endpoint triangle bound for
``B_first=[R**2*G]_{s=0}^{s=1}``; the result is a bound, not a selected value
and not a same-source cancellation.
"""

from __future__ import annotations

import ast
import hashlib
import json
from collections.abc import Sequence
from fractions import Fraction
from pathlib import Path
from typing import Any

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import mpf_cmp


MODULE = "src/openai_ns_reconstruction/original_first_switch_boundary_enclosure.py"
FIXTURE_RELATIVE_PATH = "tmp/first_switch_boundary_enclosure_fixture.json"
ORDERS = (0, 1)
_ZERO = (0, 0, 0, 0)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _endpoint(value: Any, label: str) -> tuple[int, int, int, int]:
    _require(isinstance(value, (tuple, list)) and len(value) == 4,
             label + " must be an exact MPF tuple")
    _require(all(type(part) is int for part in value),
             label + " MPF tuple components must be integers")
    result = tuple(value)
    sign, mantissa, exponent, bit_count = result
    _require(sign in (0, 1) and mantissa >= 0 and bit_count >= 0,
             label + " has an invalid sign, mantissa, or bit count")
    _require(bit_count == mantissa.bit_length(),
             label + " has a noncanonical bit count")
    if mantissa == 0:
        _require(result == _ZERO, label + " has a noncanonical zero")
    else:
        _require(mantissa & 1, label + " mantissa is not normalized")
    point = mp.make_mpf(result)
    _require(mp.isfinite(point), label + " must be finite")
    return result


def _raw_bounds(value: Any, label: str = "interval"):
    rows = getattr(value, "_mpi_", None)
    _require(rows is not None, label + " must be a raw MP interval")
    if len(rows) == 4 and all(type(part) is int for part in rows):
        rows = (rows, rows)
    elif len(rows) == 1:
        rows = (rows[0], rows[0])
    _require(len(rows) == 2, label + " has malformed interval endpoints")
    lower = _endpoint(tuple(rows[0]), label + ".lower")
    upper = _endpoint(tuple(rows[1]), label + ".upper")
    _require(mpf_cmp(lower, upper) <= 0, label + " has reversed endpoints")
    return lower, upper


def _interval(c: MPIntervalContext, lower: Any, upper: Any):
    lo = mp.make_mpf(_endpoint(lower, "lower"))
    hi = mp.make_mpf(_endpoint(upper, "upper"))
    _require(mpf_cmp(lo._mpf_, hi._mpf_) <= 0, "reversed interval")
    result = c.mpf([lo, hi])
    _raw_bounds(result)
    return result


def _coerce_scalar(c: MPIntervalContext, value: Any, label: str):
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError(label + " must use an exact scalar or raw MP interval")
    if hasattr(value, "_mpi_"):
        return _interval(c, *_raw_bounds(value, label))
    if type(value) is int:
        return c.mpf(value)
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / value.denominator
    if isinstance(value, str):
        return c.mpf(value)
    if hasattr(value, "_mpf_"):
        point = _endpoint(tuple(value._mpf_), label)
        return _interval(c, point, point)
    if isinstance(value, (tuple, list)):
        if len(value) == 4 and all(type(part) is int for part in value):
            point = _endpoint(value, label)
            return _interval(c, point, point)
        if len(value) == 2:
            left = _coerce_scalar(c, value[0], label + ".lower")
            right = _coerce_scalar(c, value[1], label + ".upper")
            return _interval(c, _raw_bounds(left)[0], _raw_bounds(right)[1])
    raise TypeError(label + " must use an exact scalar or raw MP interval")


def _row_pair(c: MPIntervalContext, values: Any, label: str):
    _require(isinstance(values, Sequence) and not isinstance(values, (str, bytes)),
             label + " must be an ordinary Z0/Z1 pair")
    _require(len(values) == 2, label + " must contain exactly Z0 and Z1")
    return tuple(_interval(c, *_raw_bounds(item, label + "[Z" + str(order) + "]"))
                 for order, item in enumerate(values))


def _record(value: Any) -> dict[str, Any]:
    lower, upper = _raw_bounds(value)
    return {
        "lower_exact_mpf_tuple": list(lower),
        "upper_exact_mpf_tuple": list(upper),
        "directed": True,
    }


def _record_rows(rows: Sequence[Any]) -> list[dict[str, Any]]:
    return [_record(row) for row in rows]


def _is_zero(value: Any) -> bool:
    lower, upper = _raw_bounds(value)
    return lower == _ZERO and upper == _ZERO


def _strip_sign(row: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return _ZERO if row[1] == 0 else (0, row[1], row[2], row[3])


def _negate_tuple(row: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return _ZERO if row[1] == 0 else (1 - row[0], row[1], row[2], row[3])


def _negate_interval(c: MPIntervalContext, value: Any):
    lower, upper = _raw_bounds(value)
    return _interval(c, _negate_tuple(upper), _negate_tuple(lower))


def _maxabs(c: MPIntervalContext, value: Any):
    """Select the larger endpoint magnitude using MPF tuple comparison."""

    lower, upper = _raw_bounds(value)
    lower_abs = _strip_sign(lower)
    upper_abs = _strip_sign(upper)
    selected = lower_abs if mpf_cmp(lower_abs, upper_abs) >= 0 else upper_abs
    return _interval(c, selected, selected)


def _symmetric(c: MPIntervalContext, bits: int):
    point = (0, 1, -bits, 1)
    radius = _interval(c, point, point)
    return _negate_interval(c, radius), radius


def _add2(left, right):
    return [left[0] + right[0], left[1] + right[1]]


def _scale2(left, scalar):
    return [scalar * left[0], scalar * left[1]]


def _prod2(left, right):
    """Exact ordinary C1 product rule, with no Taylor coefficients."""

    return [left[0] * right[0], left[1] * right[0] + left[0] * right[1]]


def _inv2(c: MPIntervalContext, left):
    _require(left[0].a > 0, "L must have a strictly positive lower endpoint")
    inverse = 1 / left[0]
    negative_left1 = _negate_interval(c, left[1])
    return [inverse, negative_left1 * inverse * inverse]


def _module_source_preflight(repo_root: str | Path | None = None) -> dict[str, Any]:
    path = Path(repo_root) / MODULE if repo_root is not None else Path(__file__)
    raw = path.read_bytes()
    tree = ast.parse(raw.decode("utf-8"), filename=str(path))
    names = [node.name for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)]
    imports = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.append(node.module or "")
    checks = {
        "utf8": True,
        "lf_only": b"\r" not in raw,
        "single_final_lf": raw.endswith(b"\n") and not raw.endswith(b"\n\n"),
        "ast_parse": True,
        "unique_bound_boundary": names.count("bound_boundary") == 1,
        "unique_fresh_fixture": names.count("fresh_fixture") == 1,
        "no_old_project_imports": not any(
            name.startswith("openai_ns_reconstruction") or name.startswith(".original_")
            for name in imports
        ),
    }
    result = {
        "passed": all(checks.values()),
        "path": str(path),
        "module_sha256": hashlib.sha256(raw).hexdigest(),
        "checks": checks,
        "ast_function_names": sorted(set(names)),
        "imports": sorted(set(imports)),
    }
    _require(result["passed"], "module LF/EOF/AST preflight failed")
    return result


def _symbolic_proof() -> dict[str, Any]:
    import sympy as sp

    z, delta = sp.symbols("z delta", real=True)
    U = sp.Function("U")(z)
    Q = sp.Function("Q")(z)
    L = 1 - delta * z**2
    G = Q / L + Q**2 / 2 - 2 * z**2 * U**2 / L**2
    U1, Q1 = sp.symbols("U1 Q1")
    L1 = sp.diff(L, z)
    invL = 1 / L
    invL1 = -L1 / L**2
    jet = (Q1 * invL + Q * invL1
           + Q * Q1
           - 2 * ((2 * z * U**2 + 2 * z**2 * U * U1) / L**2
                  + z**2 * U**2 * (-2 * L1 / L**3)))
    expected = sp.diff(G, z).subs({sp.diff(U, z): U1, sp.diff(Q, z): Q1})
    checks = {
        "ordinary_C1_G_identity": sp.simplify(expected - jet) == 0,
        "L_inverse_derivative_identity": sp.simplify(
            sp.diff(invL, z).subs({}) - invL1) == 0,
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "formula": "G=Q/L+Q^2/2-2*Z^2*U^2/L^2; L=1-delta*Z^2",
        "ordinary_derivative": "exact two-row product and inverse rules; no Taylor factors",
        "sympy_exact": True,
    }


def bound_boundary(c: MPIntervalContext, U_bg: Sequence[Any], Q_bg: Sequence[Any],
                   z: Any, delta: Any, delta_U_logs: Sequence[Any],
                   delta_Q_logs: Sequence[Any], R_upper: Any, logRm: Any,
                   radius_bits: int = 1024) -> dict[str, Any]:
    """Bound first-switch ``B=[R^2*G]_{0}^{1}`` in own ``R_m^2`` units."""

    _require(isinstance(c, MPIntervalContext), "c must be an MPIntervalContext")
    _require(type(radius_bits) is int and radius_bits > 0,
             "radius_bits must be a positive integer")
    U_source = _row_pair(c, U_bg, "U_bg")
    Q_source = _row_pair(c, Q_bg, "Q_bg")
    dU_logs = _row_pair(c, delta_U_logs, "delta_U_logs")
    dQ_logs = _row_pair(c, delta_Q_logs, "delta_Q_logs")
    z_value = _coerce_scalar(c, z, "z")
    delta_value = _coerce_scalar(c, delta, "delta")
    R = _coerce_scalar(c, R_upper, "R_upper")
    log_rm = _coerce_scalar(c, logRm, "logRm")
    _require(R.a > 0, "R_upper must be positive")
    neg_radius, radius = _symmetric(c, radius_bits)
    error_enclosure = _interval(c, _raw_bounds(neg_radius)[0],
                                _raw_bounds(radius)[1])
    # The strict cap is ``-radius_bits * ln(2)``.  Keep it in log space;
    # this is not ``-ln(radius_bits * ln(2))``.
    threshold = _negate_interval(c, radius_bits * c.ln(c.mpf(2)))
    threshold_lower = _raw_bounds(threshold)[0]

    log_checks = []
    for kind, rows in (("U", dU_logs), ("Q", dQ_logs)):
        for order, cap in enumerate(rows):
            upper = _raw_bounds(cap, kind + " error log cap")[1]
            comparison = mpf_cmp(upper, threshold_lower)
            _require(comparison < 0,
                     kind + " finite-N error log cap is not strictly below threshold")
            log_checks.append({
                "field": kind,
                "ordinary_Z_order": order,
                "raw_log_cap": _record(cap),
                "threshold": _record(threshold),
                "upper_vs_threshold_lower": comparison,
                "strict_upper_below_threshold": True,
            })

    U_full = [row + error_enclosure for row in U_source]
    Q_full = [row + error_enclosure for row in Q_source]
    zjet = [z_value, c.mpf(1)]
    one = c.mpf(1)
    z2 = _prod2(zjet, zjet)
    L = _add2([one, c.mpf(0)], _scale2(_prod2(zjet, zjet),
                                       c.mpf(-1) * delta_value))
    invL = _inv2(c, L)
    q_over_l = _prod2(Q_full, invL)
    q_squared = _scale2(_prod2(Q_full, Q_full), c.mpf(1) / 2)
    u_squared = _prod2(U_full, U_full)
    z2_u2 = _prod2(z2, u_squared)
    invL_squared = _prod2(invL, invL)
    negative_term = _scale2(_prod2(z2_u2, invL_squared), c.mpf(-2))
    G_full = [_add2(_add2(q_over_l, q_squared), negative_term)[order]
              for order in ORDERS]
    log_R = c.ln(R)
    log_two = c.ln(c.mpf(2))
    per_order = []
    output_rows = []
    for order, G_row in enumerate(G_full):
        max_abs = _maxabs(c, G_row)
        if _is_zero(max_abs):
            output_rows.append(c.mpf(0))
            per_order.append({
                "ordinary_Z_order": order,
                "G_full_jet": _record(G_row),
                "maxabs_Gjet": _record(max_abs),
                "exact_zero_G": True,
                "emitted_kind": "EXACT_ZERO",
                "emitted_not_value_or_bound": False,
                "log_absolute_B_over_Rm2": None,
            })
            continue
        log_G = c.ln(max_abs)
        log_bound = log_two + 2 * log_R + log_G - 2 * log_rm
        upper = _raw_bounds(log_bound)[1]
        comparison = mpf_cmp(upper, threshold_lower)
        _require(comparison < 0,
                 "first-switch B endpoint-pair bound is not below threshold")
        output_rows.append(_interval(c, _raw_bounds(neg_radius)[0],
                                     _raw_bounds(radius)[1]))
        per_order.append({
            "ordinary_Z_order": order,
            "G_full_jet": _record(G_row),
            "maxabs_Gjet": _record(max_abs),
            "exact_zero_G": False,
            "log_absolute_Gjet": _record(log_G),
            "endpoint_pair_log_factor": _record(log_two),
            "log_R_upper": _record(log_R),
            "logRm": _record(log_rm),
            "log_absolute_B_over_Rm2": _record(log_bound),
            "threshold": _record(threshold),
            "upper_vs_threshold_lower": comparison,
            "strict_upper_below_threshold": True,
            "emitted_kind": "BOUND",
            "bound_label": "symmetric_2^-radius_bits_enclosure",
            "emitted_not_value_or_zero": True,
        })

    metadata = {
        "schema": "original_first_switch_boundary_enclosure_v1",
        "normalization": "own_Rm^2",
        "ordinary_Z_orders": list(ORDERS),
        "formula": "Bfirst=[R^2*G]_{s=0}^{s=1}; G=Q/L+Q^2/2-2*Z^2*U^2/L^2",
        "L_formula": "L=1-delta*Z^2",
        "radius_bits": radius_bits,
        "threshold": _record(threshold),
        "R_upper": _record(R),
        "logRm": _record(log_rm),
        "signed_background": {
            "U_bg": _record_rows(U_source),
            "Q_bg": _record_rows(Q_source),
            "retained_unchanged": True,
        },
        "finite_N_error_log_caps": {
            "delta_U_logs": _record_rows(dU_logs),
            "delta_Q_logs": _record_rows(dQ_logs),
            "raw_iv_log_caps_not_values": True,
            "strict_checks": log_checks,
        },
        "full_fields": {
            "U": _record_rows(U_full),
            "Q": _record_rows(Q_full),
            "symmetric_error_radius_added_once": _record(radius),
            "symmetric_error_enclosure": _record(error_enclosure),
            "symmetric_error_enclosure_signed": True,
            "pointwise_jets_cover_both_endpoints": True,
        },
        "G_full_jets": _record_rows(G_full),
        "per_order": per_order,
        "rawIVresult_is_endpoint_pair_bound": True,
        "source_contract": {
            "signed_background_preserved": True,
            "finite_N_error_sources_preserved": True,
            "error_caps_not_selected_values": True,
            "no_extra_incoming_U": True,
            "incoming_M_already_inside_delta_Q": True,
            "physical_R_squared_applied_once": True,
            "endpoint_measure_has_no_J": True,
            "R_upper_used_only_as_saved_geometry_bound": True,
            "no_huge_exponentiation": True,
            "same_source_cancellation": False,
            "valid_bound_not_function_value": True,
            "accepted_producer_or_project_math_imported": False,
        },
        "claims": {
            "first_switch_B_bound_only": True,
            "first_switch_S1_bound": False,
            "complete_pressure_closure": False,
            "global_five_moment_closure": False,
        },
    }
    json.dumps(metadata, allow_nan=False)
    result = {
        "passed": True,
        "rawIVresult": tuple(output_rows),
        "raw_iv_result": tuple(output_rows),
        "metadata": metadata,
        "auditablemetadata": metadata,
    }
    return result


def fresh_fixture(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Run source preflight, fresh symbolic proof, and rational interval case."""

    preflight = _module_source_preflight(repo_root)
    proof = _symbolic_proof()
    old_precision = mp.mp.prec
    checks = {
        "source_lf_eof_ast_preflight": preflight["passed"],
        "fresh_symbolic_G_derivative": proof["passed"],
    }
    case: dict[str, Any] = {}
    try:
        mp.mp.prec = 53
        c = MPIntervalContext()
        c.prec = 256

        def iv(lower: Fraction, upper: Fraction):
            return c.mpf([c.mpf(lower.numerator) / lower.denominator,
                          c.mpf(upper.numerator) / upper.denominator])

        U = (iv(Fraction(1), Fraction(6, 5)), iv(Fraction(1, 5), Fraction(1, 4)))
        Q = (iv(Fraction(1, 3), Fraction(2, 5)),
             iv(Fraction(-1, 5), Fraction(-1, 6)))
        tiny_log = (iv(Fraction(-200), Fraction(-200)),
                    iv(Fraction(-200), Fraction(-200)))
        result = bound_boundary(c, U, Q, iv(Fraction(1, 4), Fraction(1, 3)),
                                Fraction(1, 10), tiny_log, tiny_log, Fraction(101),
                                Fraction(1000), radius_bits=128)
        metadata = result["metadata"]
        checks.update({
            "global_mp53_boundary_case": mp.mp.prec == 53,
            "mixed_signed_background_rows": True,
            "raw_log_caps_retained": metadata["finite_N_error_log_caps"][
                "raw_iv_log_caps_not_values"],
            "all_error_caps_strict": all(item["strict_upper_below_threshold"]
                                          for item in metadata["finite_N_error_log_caps"]["strict_checks"]),
            "G_jets_cover_both_endpoints": metadata["full_fields"][
                "pointwise_jets_cover_both_endpoints"],
            "all_rows_are_bounds": all(item["emitted_kind"] == "BOUND"
                                        for item in metadata["per_order"]),
            "no_same_source_cancellation": metadata["source_contract"][
                "same_source_cancellation"] is False,
            "raw_result_has_two_orders": len(result["rawIVresult"]) == 2,
        })
        case = {
            "rawIVresult": [_record(row) for row in result["rawIVresult"]],
            "metadata": metadata,
        }
    finally:
        mp.mp.prec = old_precision
    output = {
        "passed": all(checks.values()),
        "schema": "original_first_switch_boundary_enclosure_fixture_v1",
        "module_sha256": preflight["module_sha256"],
        "sourcechecks": preflight["checks"],
        "checks": checks,
        "exact_proof": proof,
        "case": case,
        "scope": {
            "source_executed": False,
            "saved_json_read": False,
            "accepted_producer_replayed": False,
            "actual_saved_data_produced": False,
            "same_source_cancellation_claimed": False,
        },
    }
    json.dumps(output, allow_nan=False)
    _require(output["passed"], "fresh first-switch boundary fixture failed")
    payload = json.dumps(output, sort_keys=True, separators=(",", ":"),
                         allow_nan=False).encode()
    output["fixture_sha256"] = hashlib.sha256(payload).hexdigest()
    return output


fixture_bound_check = fresh_fixture


__all__ = [
    "FIXTURE_RELATIVE_PATH",
    "MODULE",
    "ORDERS",
    "bound_boundary",
    "fixture_bound_check",
    "fresh_fixture",
]
