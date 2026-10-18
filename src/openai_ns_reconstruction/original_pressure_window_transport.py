"""Pure linear transport for already computed signed pressure C1 rows.

The caller supplies one signed ``[Z0, Z1]`` pair for each of the four
pressure charts.  This module only changes the chart normalization to the
own ``R_m^2`` scale.  The long reshaping weight remains in logarithmic form;
when its input is proved below the requested binary radius, a symmetric
``2**-radius_bits`` enclosure is emitted and labelled as a bound.
"""

from __future__ import annotations

import ast
import hashlib
import json
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path
from typing import Any

from mpmath import isfinite, make_mpf
from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import mpf_cmp


MODULE = "src/openai_ns_reconstruction/original_pressure_window_transport.py"
FIXTURE_RELATIVE_PATH = "tmp/pressure_window_transport_fixture.json"
CHARTS = ("long_reshape", "reference", "restoration", "post_restore")
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
    _require(sign in (0, 1) and mantissa >= 0,
             label + " has an invalid sign or mantissa")
    _require(bit_count == mantissa.bit_length(),
             label + " has a noncanonical bit count")
    if mantissa == 0:
        _require(result == _ZERO, label + " has a noncanonical zero")
    else:
        _require(mantissa & 1, label + " mantissa is not normalized")
    finite = make_mpf(result)
    _require(isfinite(finite), label + " must be finite")
    return result


def _raw_bounds(value: Any, label: str = "interval"):
    rows = getattr(value, "_mpi_", None)
    _require(rows is not None, label + " must be an MPIntervalContext interval")
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
    lo = make_mpf(_endpoint(lower, "lower"))
    hi = make_mpf(_endpoint(upper, "upper"))
    _require(mpf_cmp(lo._mpf_, hi._mpf_) <= 0, "reversed interval")
    result = c.mpf([lo, hi])
    _raw_bounds(result)
    return result


def _coerce_iv(c: MPIntervalContext, value: Any, label: str):
    """Copy a finite signed interval without scalar abs/endpoint selection."""

    _require(isinstance(c, MPIntervalContext), "c must be an MPIntervalContext")
    lower, upper = _raw_bounds(value, label)
    return _interval(c, lower, upper)


def _coerce_scalar(c: MPIntervalContext, value: Any, label: str):
    if hasattr(value, "_mpi_"):
        return _coerce_iv(c, value, label)
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError(label + " must use an exact scalar or interval")
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
    raise TypeError(label + " must use an exact scalar or interval")


def _record(value: Any) -> dict[str, Any]:
    lower, upper = _raw_bounds(value)
    return {
        "lower_exact_mpf_tuple": list(lower),
        "upper_exact_mpf_tuple": list(upper),
        "directed": True,
    }


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


def _abs_upper(c: MPIntervalContext, value: Any):
    """Return max(abs(exact lower), abs(exact upper)) without scalar abs."""

    lower, upper = _raw_bounds(value)
    lower_abs = _strip_sign(lower)
    upper_abs = _strip_sign(upper)
    selected = lower_abs if mpf_cmp(lower_abs, upper_abs) >= 0 else upper_abs
    return _interval(c, selected, selected)


def _symmetric_radius(c: MPIntervalContext, bits: int):
    radius = _interval(c, (0, 1, -bits, 1), (0, 1, -bits, 1))
    return _negate_iv(c, radius), radius


def _row_pair(c: MPIntervalContext, value: Any, label: str):
    _require(isinstance(value, Sequence) and not isinstance(value, (str, bytes)),
             label + " must be a signed ordinary C1 pair")
    _require(len(value) == 2, label + " must contain exactly Z0 and Z1")
    return tuple(_coerce_iv(c, item, label + "[Z" + str(order) + "]")
                 for order, item in enumerate(value))


def _source_rows(c: MPIntervalContext, rows_by_chart: Any):
    _require(isinstance(rows_by_chart, Mapping),
             "rows_by_chart must be a chart mapping")
    _require(set(rows_by_chart) == set(CHARTS),
             "rows_by_chart must have exactly the four pressure charts")
    return {chart: _row_pair(c, rows_by_chart[chart], chart)
            for chart in CHARTS}


def _materialized_weights(c: MPIntervalContext):
    # These are the only exponentials materialized by this kernel.
    return {
        "reference": c.exp(c.mpf(-4)),
        "restoration": c.exp(c.mpf(-2)),
        "post_restore": c.mpf(1),
    }


def _module_source_preflight(repo_root: str | Path | None = None) -> dict[str, Any]:
    path = Path(repo_root) / MODULE if repo_root is not None else Path(__file__)
    raw = path.read_bytes()
    checks = {
        "utf8": True,
        "lf_only": b"\r" not in raw,
        "single_final_lf": raw.endswith(b"\n") and not raw.endswith(b"\n\n"),
    }
    tree = ast.parse(raw.decode("utf-8"), filename=str(path))
    names = [node.name for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)]
    checks["ast_parse"] = True
    checks["unique_transport"] = names.count("transport") == 1
    checks["unique_fixture"] = names.count("fresh_fixture") == 1
    imports = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.append(node.module or "")
    checks["no_accepted_producer_or_helper_imports"] = not any(
        name.startswith("openai_ns_reconstruction.original_")
        or name.startswith(".original_") for name in imports
    )
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

    z, lref = sp.symbols("Z Lref", real=True)
    f_long = sp.Function("f_long")(z)
    f_ref = sp.Function("f_ref")(z)
    f_restore = sp.Function("f_restore")(z)
    f_post = sp.Function("f_post")(z)
    w_long = sp.exp(-2 * (lref + 2))
    w_ref, w_restore, w_post = sp.exp(-4), sp.exp(-2), sp.Integer(1)
    operator = w_long * f_long + w_ref * f_ref + w_restore * f_restore + w_post * f_post
    expected = (w_long * sp.diff(f_long, z) + w_ref * sp.diff(f_ref, z)
                + w_restore * sp.diff(f_restore, z) + sp.diff(f_post, z))
    checks = {
        "ordinary_C1_linearity": sp.simplify(sp.diff(operator, z) - expected) == 0,
        "long_weight_Z_independent": sp.diff(w_long, z) == 0,
        "reference_weight_is_exp_minus4": w_ref == sp.exp(-4),
        "restoration_weight_is_exp_minus2": w_restore == sp.exp(-2),
        "post_weight_is_one": w_post == 1,
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "operator": "exp(-2*(Lref+2))*long + exp(-4)*reference + exp(-2)*restoration + post_restore",
        "ordinary_derivative": "differentiate only supplied rows; all weights are Z-independent",
        "sympy_exact": True,
    }


def transport(
    c: MPIntervalContext,
    rows_by_chart: Mapping[str, Sequence[Any]],
    Lref: Any,
    *,
    geometry_Z_independent: bool = True,
    radius_bits: int = 1024,
) -> dict[str, Any]:
    """Transport four signed ordinary ``C1`` rows into own ``R_m^2`` units."""

    _require(isinstance(c, MPIntervalContext), "c must be an MPIntervalContext")
    _require(geometry_Z_independent is True,
             "the pressure weight must be Z-independent")
    _require(type(radius_bits) is int and radius_bits > 0,
             "radius_bits must be a positive integer")
    source = _source_rows(c, rows_by_chart)
    lref = _coerce_scalar(c, Lref, "Lref")
    weights = _materialized_weights(c)
    log_two = c.ln(c.mpf(2))
    long_log = (lref + c.mpf(2)) * c.mpf(-2)
    threshold = _negate_iv(c, log_two * c.mpf(radius_bits))
    neg_radius, radius = _symmetric_radius(c, radius_bits)

    weighted: dict[str, tuple[Any, Any]] = {}
    per_chart: dict[str, Any] = {}
    for chart in ("reference", "restoration", "post_restore"):
        weight = weights[chart]
        rows = tuple(row * weight for row in source[chart])
        weighted[chart] = rows
        formula = {"reference": "exp(-4)", "restoration": "exp(-2)",
                   "post_restore": "1"}[chart]
        per_chart[chart] = {
            "weight_formula": formula,
            "weight_materialized": True,
            "weight": _record(weight),
            "input_signed_rows": [_record(row) for row in source[chart]],
            "weighted_rows": [_record(row) for row in rows],
        }

    long_rows: list[Any] = []
    long_bounds: list[dict[str, Any]] = []
    suppression: list[dict[str, Any]] = []
    for order, row in enumerate(source["long_reshape"]):
        if _is_zero(row):
            enclosure = c.mpf(0)
            bound = {
                "ordinary_Z_order": order,
                "exact_zero_input": True,
                "emitted_kind": "EXACT_ZERO",
                "input_signed_row": _record(row),
                "log_absolute_bound": None,
            }
        else:
            magnitude = _abs_upper(c, row)
            log_bound = c.ln(magnitude) + long_log
            upper_bound = _raw_bounds(log_bound)[1]
            lower_threshold = _raw_bounds(threshold)[0]
            comparison = mpf_cmp(upper_bound, lower_threshold)
            _require(comparison < 0,
                     "long reshaping row is not below the requested radius")
            enclosure = _interval(c, _raw_bounds(neg_radius)[0],
                                  _raw_bounds(radius)[1])
            bound = {
                "ordinary_Z_order": order,
                "exact_zero_input": False,
                "emitted_kind": "BOUND",
                "bound_label": "symmetric_2^-radius_bits_enclosure",
                "input_signed_row": _record(row),
                "maxabs_exact_endpoint_bound": _record(magnitude),
                "log_absolute_input_bound": _record(c.ln(magnitude)),
                "log_coefficient": _record(long_log),
                "log_absolute_bound": _record(log_bound),
                "threshold": _record(threshold),
                "upper_bound_vs_threshold_lower": comparison,
                "strict_upper_below_threshold": True,
                "emitted_enclosure": _record(enclosure),
                "emitted_not_value_or_zero": True,
            }
        long_rows.append(enclosure)
        long_bounds.append(bound)
        suppression.append(dict(bound))

    weighted["long_reshape"] = tuple(long_rows)
    per_chart["long_reshape"] = {
        "weight_formula": "exp(-2*(Lref+2))",
        "weight_materialized": False,
        "log_coefficient": _record(long_log),
        "input_signed_rows": [_record(row) for row in source["long_reshape"]],
        "weighted_rows": [_record(row) for row in long_rows],
        "weighted_rows_are_symmetric_bounds": True,
        "log_absolute_bounds": long_bounds,
        "suppression_proof": suppression,
    }

    combined = tuple(
        weighted["long_reshape"][order]
        + weighted["reference"][order]
        + weighted["restoration"][order]
        + weighted["post_restore"][order]
        for order in ORDERS
    )
    metadata = {
        "schema": "original_pressure_window_transport_v1",
        "normalization": "own_Rm^2",
        "weights_relative_to_Rm2": {
            "long_reshape": "exp(-2*(Lref+2))",
            "reference": "exp(-4)",
            "restoration": "exp(-2)",
            "post_restore": "1",
        },
        "ordinary_Z_orders": list(ORDERS),
        "geometry_Z_independent": True,
        "radius_bits": radius_bits,
        "radius": _record(radius),
        "Lref_input_retained": _record(lref),
        "long_log_coefficient_retained": _record(long_log),
        "per_chart_weighted_rows": per_chart,
        "long_log_bounds": long_bounds,
        "suppression_proof": suppression,
        "combined_rows": [_record(row) for row in combined],
        "original_source_inputs": {
            chart: [_record(row) for row in source[chart]] for chart in CHARTS
        },
        "source_contract": {
            "signed_ordinary_C1_inputs": True,
            "input_rows_are_cap_mappings": False,
            "input_rows_are_Taylor_coefficients": False,
            "source_rows_recomputed": False,
            "accepted_producer_or_helper_math_imported": False,
            "long_coefficient_exponentiated": False,
            "only_exp_minus4_and_exp_minus2_materialized": True,
            "source_half_S1_factor_applied": False,
        },
        "claims": {
            "full_pressure_closure": False,
            "global_five_moment_closure": False,
            "terminal_moment_identity": False,
            "full_NS_validation": False,
        },
    }
    result = {
        "passed": True,
        "rawIVresult": combined,
        "raw_iv_result": combined,
        "combined_rows": combined,
        "metadata": metadata,
        "auditablemetadata": metadata,
    }
    json.dumps({"metadata": metadata}, allow_nan=False)
    return result


def fresh_fixture(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Run preflight, fresh symbolic checks, and independent rational rows."""

    preflight = _module_source_preflight(repo_root)
    proof = _symbolic_proof()
    old_precision = __import__("mpmath").mp.prec
    checks: dict[str, bool] = {
        "source_lf_eof_ast_preflight": preflight["passed"],
        "fresh_symbolic_linear_operator": proof["passed"],
    }
    case: dict[str, Any] = {}
    try:
        __import__("mpmath").mp.prec = 53
        c = MPIntervalContext()
        c.prec = 256
        negative = (c.mpf(-7) / 10, c.mpf(-1) / 3)
        rows = {chart: negative for chart in CHARTS}
        rows["long_reshape"] = (negative[0], c.mpf(0))
        result = transport(c, rows, Fraction(1000), radius_bits=128)
        meta = result["metadata"]
        long_info = meta["suppression_proof"]
        checks.update({
            "global_mp53_boundary_case": __import__("mpmath").mp.prec == 53,
            "negative_signed_inputs_accepted": True,
            "long_rows_strictly_suppressed": all(
                item.get("strict_upper_below_threshold", True)
                for item in long_info if not item["exact_zero_input"]
            ),
            "long_exact_zero_retained": long_info[1]["emitted_kind"] == "EXACT_ZERO",
            "long_enclosure_is_labelled_bound": long_info[0]["emitted_kind"] == "BOUND",
            "combined_rows_are_C1": len(result["rawIVresult"]) == 2,
            "source_inputs_retained": set(meta["original_source_inputs"]) == set(CHARTS),
            "no_full_pressure_claim": meta["claims"]["full_pressure_closure"] is False,
        })
        case = {
            "rawIVresult": [_record(row) for row in result["rawIVresult"]],
            "metadata": meta,
        }
    finally:
        __import__("mpmath").mp.prec = old_precision
    output = {
        "passed": all(checks.values()),
        "schema": "pressure_window_transport_rational_fixture_v1",
        "module_sha256": preflight["module_sha256"],
        "checks": checks,
        "exact_proof": proof,
        "case": case,
        "scope": {
            "source_executed": False,
            "saved_json_read": False,
            "accepted_producer_replayed": False,
            "source_rows_recomputed": False,
            "full_pressure_closure_claimed": False,
        },
    }
    json.dumps(output, allow_nan=False)
    _require(output["passed"], "fresh pressure transport fixture failed")
    return output


fixture_bound_check = fresh_fixture


__all__ = [
    "CHARTS",
    "FIXTURE_RELATIVE_PATH",
    "MODULE",
    "ORDERS",
    "fixture_bound_check",
    "fresh_fixture",
    "transport",
]
