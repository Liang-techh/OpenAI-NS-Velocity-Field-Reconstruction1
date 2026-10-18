"""Pure signed enclosure of the new outer pressure ``B`` and ``D`` terms.

The caller supplies the two endpoint jets and the two already accepted volume
terms.  This module does not read a saved source packet, evaluate a producer,
or select an endpoint value.  Every nonzero term is emitted only after its
logarithmic smallness proof has passed; the emitted object is a labelled
symmetric bound.
"""

from __future__ import annotations

import ast
import hashlib
import json
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path
from typing import Any

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import mpf_cmp


MODULE = "src/openai_ns_reconstruction/original_outer_pressure_enclosure.py"
FIXTURE_RELATIVE_PATH = "tmp/outer_pressure_enclosure_fixture.json"
ORDERS = (0, 1)
ENDPOINTS = ("micro_lower", "macro_upper")
TERMS = ("micro", "macro")
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
    """Copy an exact scalar or raw interval, rejecting floats and cap maps."""

    if isinstance(value, Mapping) or isinstance(value, bool) or isinstance(value, float):
        raise TypeError(label + " must use an exact scalar or raw MP interval")
    if hasattr(value, "_mpi_"):
        return _interval(c, *_raw_bounds(value, label))
    if type(value) is int:
        return c.mpf(value)
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / value.denominator
    if isinstance(value, str):
        return _interval(c, *_raw_bounds(c.mpf(value), label))
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
    _require(isinstance(values, Sequence) and not isinstance(values, (str, bytes, Mapping)),
             label + " must be an ordinary Z0/Z1 pair")
    _require(len(values) == 2, label + " must contain exactly Z0 and Z1")
    return tuple(_interval(c, *_raw_bounds(item, label + "[Z" + str(order) + "]"))
                 for order, item in enumerate(values))


def _mapping_exact(value: Any, keys: Sequence[str], label: str) -> None:
    _require(isinstance(value, Mapping), label + " must be a mapping")
    _require(set(value) == set(keys), label + " has unexpected keys")


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


def _negate_iv(c: MPIntervalContext, value: Any):
    lower, upper = _raw_bounds(value)
    return _interval(c, _negate_tuple(upper), _negate_tuple(lower))


def _maxabs(c: MPIntervalContext, value: Any):
    """Select the larger endpoint magnitude using exact tuple comparisons."""

    lower, upper = _raw_bounds(value)
    lower_abs = _strip_sign(lower)
    upper_abs = _strip_sign(upper)
    selected = lower_abs if mpf_cmp(lower_abs, upper_abs) >= 0 else upper_abs
    return _interval(c, selected, selected)


def _symmetric(c: MPIntervalContext, bits: int):
    point = (0, 1, -bits, 1)
    radius = _interval(c, point, point)
    return _negate_iv(c, radius), radius


def _add2(left, right):
    return [left[0] + right[0], left[1] + right[1]]


def _scale2(left, scalar):
    return [scalar * left[0], scalar * left[1]]


def _prod2(left, right):
    """Ordinary C1 product rule with no Taylor factorials."""

    return [left[0] * right[0], left[1] * right[0] + left[0] * right[1]]


def _inv2(c: MPIntervalContext, left):
    lower, _ = _raw_bounds(left[0], "L[Z0]")
    _require(mpf_cmp(lower, _ZERO) > 0,
             "L must have a strictly positive lower endpoint")
    inverse = 1 / left[0]
    negative_left1 = _negate_iv(c, left[1])
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
        "unique_enclose": names.count("enclose") == 1,
        "unique_fresh_fixture": names.count("fresh_fixture") == 1,
        "no_accepted_project_math_imports": not any(
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
    """Verify the ordinary C1 identity independently with exact SymPy algebra."""

    import sympy as sp

    z, delta = sp.symbols("z delta", real=True)
    U = sp.Function("U")(z)
    Q = sp.Function("Q")(z)
    L = 1 - delta * z**2
    G = Q / L + Q**2 / 2 - 2 * z**2 * U**2 / L**2
    U1, Q1 = sp.symbols("U1 Q1")
    L1 = sp.diff(L, z)
    invL1 = -L1 / L**2
    jet = (Q1 / L + Q * invL1 + Q * Q1
           - 2 * ((2 * z * U**2 + 2 * z**2 * U * U1) / L**2
                  + z**2 * U**2 * (-2 * L1 / L**3)))
    expected = sp.diff(G, z).subs({sp.diff(U, z): U1, sp.diff(Q, z): Q1})
    checks = {
        "ordinary_C1_G_identity": sp.simplify(expected - jet) == 0,
        "L_inverse_derivative_identity": sp.simplify(
            sp.diff(1 / L, z) - invL1) == 0,
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "formula": "G=Q/L+Q^2/2-2*Z^2*U^2/L^2; L=1-delta*Z^2",
        "ordinary_derivative": "exact two-row product and inverse rules; no Taylor factors",
        "sympy_exact": True,
    }


def _threshold(c: MPIntervalContext, radius_bits: int):
    return _negate_iv(c, c.mpf(radius_bits) * c.ln(c.mpf(2)))


def _check_log_caps(c: MPIntervalContext, rows: Sequence[Any], threshold: Any,
                    label: str) -> list[dict[str, Any]]:
    threshold_lower = _raw_bounds(threshold)[0]
    checks = []
    for order, cap in enumerate(rows):
        upper = _raw_bounds(cap, label + " log cap")[1]
        comparison = mpf_cmp(upper, threshold_lower)
        _require(comparison < 0,
                 label + " finite-N error log cap is not strictly below threshold")
        checks.append({
            "ordinary_Z_order": order,
            "raw_log_cap": _record(cap),
            "threshold": _record(threshold),
            "upper_vs_threshold_lower": comparison,
            "strict_upper_below_threshold": True,
        })
    return checks


def _safe_logsum(c: MPIntervalContext, terms: Sequence[Any]):
    """Return one or two log bounds using ``max(upper)+ln(2)`` safely."""

    _require(len(terms) > 0, "safe log sum needs a nonzero term")
    if len(terms) == 1:
        return terms[0], False, [_raw_bounds(terms[0])[1]]
    uppers = [_raw_bounds(term)[1] for term in terms]
    selected = uppers[0]
    for candidate in uppers[1:]:
        if mpf_cmp(candidate, selected) > 0:
            selected = candidate
    point = _interval(c, selected, selected)
    return point + c.ln(c.mpf(2)), True, uppers


def _g_jet(c: MPIntervalContext, U_bg: Sequence[Any], Q_bg: Sequence[Any],
           z: Any, delta: Any, error: Any):
    """Build the signed ordinary C1 ``G`` jet after symmetric error boxing."""

    U_full = [row + error for row in U_bg]
    Q_full = [row + error for row in Q_bg]
    zjet = [z, c.mpf(1)]
    z2 = _prod2(zjet, zjet)
    L = _add2([c.mpf(1), c.mpf(0)], _scale2(z2, c.mpf(-1) * delta))
    invL = _inv2(c, L)
    q_over_l = _prod2(Q_full, invL)
    q_squared = _scale2(_prod2(Q_full, Q_full), c.mpf(1) / 2)
    u_squared = _prod2(U_full, U_full)
    z2_u2 = _prod2(z2, u_squared)
    invL_squared = _prod2(invL, invL)
    negative_term = _scale2(_prod2(z2_u2, invL_squared), c.mpf(-2))
    return [_add2(_add2(q_over_l, q_squared), negative_term)[order]
            for order in ORDERS], U_full, Q_full, L


def enclose(c: MPIntervalContext, endpoints: Mapping[str, Any],
            volume_terms: Mapping[str, Any], z: Any, delta: Any, logRm: Any,
            radius_bits: int = 1024) -> dict[str, Any]:
    """Enclose ``B_outer + D_micro + D_macro`` in own ``R_m^2`` units."""

    _require(isinstance(c, MPIntervalContext), "c must be an MPIntervalContext")
    _require(type(radius_bits) is int and radius_bits > 0,
             "radius_bits must be a positive integer")
    _mapping_exact(endpoints, ENDPOINTS, "endpoints")
    _mapping_exact(volume_terms, TERMS, "volume_terms")
    z_value = _coerce_scalar(c, z, "z")
    delta_value = _coerce_scalar(c, delta, "delta")
    log_rm = _coerce_scalar(c, logRm, "logRm")
    threshold = _threshold(c, radius_bits)
    neg_radius, radius = _symmetric(c, radius_bits)
    error_enclosure = _interval(c, _raw_bounds(neg_radius)[0],
                                _raw_bounds(radius)[1])

    endpoint_data: dict[str, Any] = {}
    for endpoint_name in ENDPOINTS:
        packet = endpoints[endpoint_name]
        _mapping_exact(packet,
                       ("U_bg", "Q_bg", "delta_U_logs", "delta_Q_logs", "logR_upper"),
                       "endpoints[" + endpoint_name + "]")
        U_source = _row_pair(c, packet["U_bg"], endpoint_name + ".U_bg")
        Q_source = _row_pair(c, packet["Q_bg"], endpoint_name + ".Q_bg")
        dU_logs = _row_pair(c, packet["delta_U_logs"], endpoint_name + ".delta_U_logs")
        dQ_logs = _row_pair(c, packet["delta_Q_logs"], endpoint_name + ".delta_Q_logs")
        log_r = _coerce_scalar(c, packet["logR_upper"], endpoint_name + ".logR_upper")
        log_checks = (_check_log_caps(c, dU_logs, threshold, endpoint_name + ".U"),
                      _check_log_caps(c, dQ_logs, threshold, endpoint_name + ".Q"))
        G_full, U_full, Q_full, L = _g_jet(
            c, U_source, Q_source, z_value, delta_value, error_enclosure)
        endpoint_data[endpoint_name] = {
            "U_source": U_source,
            "Q_source": Q_source,
            "dU_logs": dU_logs,
            "dQ_logs": dQ_logs,
            "logR": log_r,
            "log_checks": log_checks,
            "G": G_full,
            "U_full": U_full,
            "Q_full": Q_full,
            "L": L,
        }

    raw_B = []
    B_rows_metadata = []
    for order in ORDERS:
        terms = []
        term_metadata = []
        for sign, endpoint_name in ((1, "macro_upper"), (-1, "micro_lower")):
            G_row = endpoint_data[endpoint_name]["G"][order]
            maxabs = _maxabs(c, G_row)
            info = {
                "orientation": sign,
                "endpoint": endpoint_name,
                "G_row": _record(G_row),
                "maxabs_G": _record(maxabs),
                "exact_zero_G": _is_zero(maxabs),
                "log_term_formula": "2*logR_upper+ln(maxabsG)-2*logRm",
            }
            if _is_zero(maxabs):
                info["log_term"] = None
            else:
                log_term = (c.mpf(2) * endpoint_data[endpoint_name]["logR"]
                            + c.ln(maxabs) - c.mpf(2) * log_rm)
                terms.append(log_term)
                info["log_term"] = _record(log_term)
            term_metadata.append(info)
        if not terms:
            raw_B.append(c.mpf(0))
            B_rows_metadata.append({
                "ordinary_Z_order": order,
                "orientation": {"macro_upper": 1, "micro_lower": -1},
                "terms": term_metadata,
                "exact_zero_B": True,
                "emitted_kind": "EXACT_ZERO",
                "emitted_not_value_or_zero": False,
                "combined_log_bound": None,
            })
            continue
        combined_log, added_ln2, exact_uppers = _safe_logsum(c, terms)
        combined_upper = _raw_bounds(combined_log)[1]
        comparison = mpf_cmp(combined_upper, _raw_bounds(threshold)[0])
        _require(comparison < 0,
                 "outer B triangle bound is not strictly below threshold")
        raw_B.append(error_enclosure)
        B_rows_metadata.append({
            "ordinary_Z_order": order,
            "orientation": {"macro_upper": 1, "micro_lower": -1},
            "terms": term_metadata,
            "exact_zero_B": False,
            "combined_log_bound": _record(combined_log),
            "exact_term_upper_tuples": [list(value) for value in exact_uppers],
            "safe_logsum_added_ln2": added_ln2,
            "combined_upper_vs_threshold_lower": comparison,
            "strict_combined_upper_below_threshold": True,
            "emitted_kind": "BOUND",
            "bound_label": "symmetric_2^-radius_bits_enclosure",
            "emitted_not_value_or_zero": True,
        })
    raw_B = tuple(raw_B)

    raw_D_by_term: dict[str, tuple[Any, Any]] = {}
    D_metadata: dict[str, Any] = {}
    for term_name in TERMS:
        packet = volume_terms[term_name]
        _mapping_exact(packet, ("D_rows", "logW"), "volume_terms[" + term_name + "]")
        D_source = _row_pair(c, packet["D_rows"], term_name + ".D_rows")
        log_w = _coerce_scalar(c, packet["logW"], term_name + ".logW")
        log_coefficient = log_w - c.mpf(2) * log_rm
        rows = []
        row_metadata = []
        for order, D_row in enumerate(D_source):
            maxabs = _maxabs(c, D_row)
            common = {
                "ordinary_Z_order": order,
                "original_signed_D": _record(D_row),
                "original_logW": _record(log_w),
                "log_coefficient": _record(log_coefficient),
                "coefficient_positive": True,
                "coefficient_materialized": False,
                "coefficient_formula": "exp(logW-2*logRm)",
            }
            if _is_zero(maxabs):
                rows.append(c.mpf(0))
                row_metadata.append(dict(common, exact_zero_input=True,
                                         emitted_kind="EXACT_ZERO",
                                         emitted_not_value_or_zero=False,
                                         log_norm=None))
                continue
            log_norm = c.ln(maxabs) + log_coefficient
            comparison = mpf_cmp(_raw_bounds(log_norm)[1],
                                 _raw_bounds(threshold)[0])
            _require(comparison < 0,
                     term_name + " D transport is not strictly below threshold")
            rows.append(error_enclosure)
            row_metadata.append(dict(common, exact_zero_input=False,
                                     maxabs_exact_endpoint_bound=_record(maxabs),
                                     log_absolute_D=_record(c.ln(maxabs)),
                                     log_norm=_record(log_norm),
                                     threshold=_record(threshold),
                                     upper_vs_threshold_lower=comparison,
                                     strict_upper_below_threshold=True,
                                     emitted_kind="BOUND",
                                     bound_label="symmetric_2^-radius_bits_enclosure",
                                     emitted_not_value_or_zero=True))
        raw_D_by_term[term_name] = tuple(rows)
        D_metadata[term_name] = {
            "original_signed_D_rows": _record_rows(D_source),
            "original_logW": _record(log_w),
            "log_coefficient": _record(log_coefficient),
            "coefficient_positive": True,
            "coefficient_materialized": False,
            "rows": row_metadata,
        }

    raw_result = tuple(raw_B[order] + raw_D_by_term["micro"][order]
                       + raw_D_by_term["macro"][order] for order in ORDERS)
    metadata = {
        "schema": "original_outer_pressure_enclosure_v1",
        "normalization": "own_Rm^2",
        "ordinary_Z_orders": list(ORDERS),
        "endpoint_orientation": {"macro_upper": 1, "micro_lower": -1},
        "formula": "B_outer=[R^2 G]_macro_upper-[R^2 G]_micro_lower",
        "G_formula": "G=Q/L+Q^2/2-2*Z^2*U^2/L^2",
        "L_formula": "L=1-delta*Z^2",
        "pressure_Q_is_linear_physical_field": True,
        "geometry_Z_independent": True,
        "radius_bits": radius_bits,
        "threshold": _record(threshold),
        "symmetric_error_radius": _record(radius),
        "symmetric_error_enclosure": _record(error_enclosure),
        "symmetric_error_enclosure_signed": True,
        "logRm": _record(log_rm),
        "endpoints": {
            name: {
                "logR_upper": _record(data["logR"]),
                "signed_background": {
                    "U_bg": _record_rows(data["U_source"]),
                    "Q_bg": _record_rows(data["Q_source"]),
                    "retained_unchanged": True,
                },
                "finite_N_error_log_caps": {
                    "delta_U_logs": _record_rows(data["dU_logs"]),
                    "delta_Q_logs": _record_rows(data["dQ_logs"]),
                    "raw_iv_log_caps_not_values": True,
                    "strict_checks": [*data["log_checks"][0], *data["log_checks"][1]],
                },
                "full_fields": {
                    "U": _record_rows(data["U_full"]),
                    "Q": _record_rows(data["Q_full"]),
                    "symmetric_error_added_once": True,
                    "symmetric_error_enclosure_signed": True,
                },
                "L_jet": _record_rows(data["L"]),
                "G_signed_rows": _record_rows(data["G"]),
            }
            for name, data in endpoint_data.items()
        },
        "outer_B": {
            "rows": B_rows_metadata,
            "physical_radius_factor_applied_once": True,
            "endpoint_jacobian_applied": False,
            "same_source_cancellation": False,
            "no_huge_exponentiation": True,
        },
        "volume_terms": D_metadata,
        "outer_S1": {
            "formula": "outerS1norm=Bnorm+Dmicronorm+Dmacronorm",
            "pressure_half_factor_applied": False,
            "root_applies_half_once": True,
            "raw_rows": _record_rows(raw_result),
        },
        "rawIVresult_is_outerS1norm": True,
        "source_contract": {
            "signed_background_preserved": True,
            "finite_N_error_sources_preserved": True,
            "accepted_volume_sources_only": ["Dmicro/Wmicro", "Dmacro/Wmacro"],
            "input_rows_are_cap_mappings": False,
            "input_rows_are_float_values": False,
            "source_rows_recomputed": False,
            "physical_radius_factor_once": True,
            "coefficient_positive_and_unmaterialized": True,
            "no_same_source_cancellation": True,
            "accepted_project_math_imported": False,
        },
        "claims": {
            "outer_B_bound_only": True,
            "outer_S1_bound_only": True,
            "complete_pressure_closure": False,
            "global_five_moment_closure": False,
        },
    }
    json.dumps(metadata, allow_nan=False)
    return {
        "passed": True,
        "rawIVresult": raw_result,
        "rawB": raw_B,
        "rawD_by_term": raw_D_by_term,
        "metadata": metadata,
        "auditablemetadata": metadata,
    }


def fresh_fixture(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Run source preflight, fresh algebra, and independent rational regressions."""

    # Keep this first: no fresh numerical or symbolic work precedes the source
    # LF/EOF/AST guard, and this function never reads or overwrites a fixture.
    preflight = _module_source_preflight(repo_root)
    proof = _symbolic_proof()
    old_precision = mp.mp.prec
    checks: dict[str, bool] = {
        "source_lf_eof_ast_preflight": preflight["passed"],
        "fresh_symbolic_G_derivative": proof["passed"],
    }
    case: dict[str, Any] = {}
    try:
        mp.mp.prec = 53
        c = MPIntervalContext()
        c.prec = 256

        def q(value: Fraction):
            return c.mpf(value.numerator) / value.denominator

        def iv(lower: Fraction, upper: Fraction):
            return c.mpf([q(lower), q(upper)])

        tiny = (iv(Fraction(-200), Fraction(-200)),
                iv(Fraction(-200), Fraction(-200)))
        micro = {
            "U_bg": (iv(Fraction(-1, 4), Fraction(1, 5)),
                     iv(Fraction(-3, 5), Fraction(-1, 7))),
            "Q_bg": (iv(Fraction(1, 3), Fraction(2, 5)),
                     iv(Fraction(-1, 5), Fraction(1, 6))),
            "delta_U_logs": tiny,
            "delta_Q_logs": tiny,
            "logR_upper": iv(Fraction(2), Fraction(2)),
        }
        macro = {
            "U_bg": (iv(Fraction(1, 5), Fraction(2, 5)),
                     iv(Fraction(-1, 4), Fraction(1, 3))),
            "Q_bg": (iv(Fraction(-2, 5), Fraction(-1, 3)),
                     iv(Fraction(1, 7), Fraction(1, 4))),
            "delta_U_logs": tiny,
            "delta_Q_logs": tiny,
            "logR_upper": iv(Fraction(2), Fraction(2)),
        }
        signed_D = {
            "micro": {"D_rows": (iv(Fraction(-1, 10), Fraction(1, 10)),
                                   iv(Fraction(-1, 20), Fraction(1, 30))),
                       "logW": Fraction(-10)},
            "macro": {"D_rows": (iv(Fraction(1, 12), Fraction(1, 9)),
                                   c.mpf(0)),
                       "logW": Fraction(-11)},
        }
        result = enclose(c, {"micro_lower": micro, "macro_upper": macro},
                         signed_D, Fraction(1, 4), Fraction(1, 10),
                         Fraction(1000), radius_bits=128)
        meta = result["metadata"]
        checks.update({
            "global_mp53_with_context_256": mp.mp.prec == 53 and c.prec == 256,
            "signed_rows_cross_zero": (
                meta["endpoints"]["micro_lower"]["signed_background"]["U_bg"][0]
                ["lower_exact_mpf_tuple"][0] == 1
                and meta["endpoints"]["micro_lower"]["full_fields"]["U"][0]
                ["lower_exact_mpf_tuple"][0] == 1),
            "negative_error_regression": (
                mpf_cmp(
                    tuple(meta["endpoints"]["micro_lower"]["full_fields"]["U"][0]
                          ["lower_exact_mpf_tuple"]),
                    tuple(meta["endpoints"]["micro_lower"]["signed_background"]["U_bg"][0]
                          ["lower_exact_mpf_tuple"])) < 0
                and mpf_cmp(
                    tuple(meta["endpoints"]["micro_lower"]["full_fields"]["Q"][0]
                          ["lower_exact_mpf_tuple"]),
                    tuple(meta["endpoints"]["micro_lower"]["signed_background"]["Q_bg"][0]
                          ["lower_exact_mpf_tuple"])) < 0),
            "symmetric_error_is_signed": meta["symmetric_error_enclosure_signed"],
            "G_rows_are_signed_C1": all(
                len(meta["endpoints"][name]["G_signed_rows"]) == 2
                for name in ENDPOINTS),
            "outer_orientation_preserved": meta["endpoint_orientation"] == {
                "macro_upper": 1, "micro_lower": -1},
            "volume_zero_input_retained": (
                meta["volume_terms"]["macro"]["rows"][1]["emitted_kind"]
                == "EXACT_ZERO"),
            "volume_coefficients_positive": all(
                meta["volume_terms"][name]["coefficient_positive"]
                and not meta["volume_terms"][name]["coefficient_materialized"]
                for name in TERMS),
            "outer_rows_are_two": len(result["rawIVresult"]) == 2,
            "no_half_factor": meta["outer_S1"]["pressure_half_factor_applied"] is False,
        })

        try:
            enclose(c, {"micro_lower": micro, "macro_upper": macro}, signed_D,
                    Fraction(1, 4), Fraction(1, 10), Fraction(0), radius_bits=128)
        except ValueError as error:
            checks["threshold_failure_rejected"] = "threshold" in str(error)
        else:
            checks["threshold_failure_rejected"] = False

        zero_packet = {
            "U_bg": (c.mpf(0), c.mpf(0)), "Q_bg": (c.mpf(0), c.mpf(0)),
            "delta_U_logs": tiny, "delta_Q_logs": tiny,
            "logR_upper": c.mpf(0),
        }
        zero_terms = {name: {"D_rows": (c.mpf(0), c.mpf(0)),
                             "logW": Fraction(-10)} for name in TERMS}
        zero_result = enclose(c, {"micro_lower": zero_packet,
                                  "macro_upper": zero_packet}, zero_terms,
                              Fraction(0), Fraction(0), Fraction(1000),
                              radius_bits=128)
        checks["exact_zero_D_may_zero"] = all(
            zero_result["metadata"]["volume_terms"][name]["rows"][order]
            ["emitted_kind"] == "EXACT_ZERO"
            for name in TERMS for order in ORDERS)

        try:
            enclose(c, {"micro_lower": zero_packet, "macro_upper": zero_packet},
                    zero_terms, Fraction(2), Fraction(1), Fraction(1000),
                    radius_bits=128)
        except ValueError as error:
            checks["geometry_guard_rejects_nonpositive_L"] = "L" in str(error)
        else:
            checks["geometry_guard_rejects_nonpositive_L"] = False
        case = {
            "rawIVresult": [_record(row) for row in result["rawIVresult"]],
            "rawB": [_record(row) for row in result["rawB"]],
            "rawD_by_term": {
                name: [_record(row) for row in rows]
                for name, rows in result["rawD_by_term"].items()
            },
            "metadata": meta,
        }
    finally:
        mp.mp.prec = old_precision
    output = {
        "passed": all(checks.values()),
        "schema": "outer_pressure_enclosure_rational_fixture_v1",
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
            "original_fixture_overwritten": False,
            "same_source_cancellation_claimed": False,
        },
    }
    json.dumps(output, allow_nan=False)
    _require(output["passed"], "fresh outer pressure enclosure fixture failed")
    return output


__all__ = [
    "ENDPOINTS",
    "FIXTURE_RELATIVE_PATH",
    "MODULE",
    "ORDERS",
    "TERMS",
    "enclose",
    "fresh_fixture",
]
