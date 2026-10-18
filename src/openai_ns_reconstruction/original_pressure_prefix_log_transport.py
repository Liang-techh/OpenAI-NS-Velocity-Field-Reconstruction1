"""Pure logarithmic transport for the three signed pressure prefix atoms.

The caller supplies ordinary ``Z0``/``Z1`` signed enclosures and each atom's
own logarithmic weight.  This kernel keeps the coefficient in log form.  It
proves a requested binary radius before replacing a nonzero weighted row by a
labelled symmetric enclosure; an unproved row is rejected rather than capped
or silently set to zero.
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


MODULE = "src/openai_ns_reconstruction/original_pressure_prefix_log_transport.py"
FIXTURE_RELATIVE_PATH = "tmp/pressure_prefix_log_transport_fixture.json"
CHARTS = ("first_switch", "second_switch", "post_power")
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


def _raw_bounds(value: Any, label: str = "interval") -> tuple[tuple[int, int, int, int],
                                                                   tuple[int, int, int, int]]:
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
    lo = mp.make_mpf(_endpoint(lower, "lower"))
    hi = mp.make_mpf(_endpoint(upper, "upper"))
    _require(mpf_cmp(lo._mpf_, hi._mpf_) <= 0, "reversed interval")
    result = c.mpf([lo, hi])
    _raw_bounds(result)
    return result


def _coerce_scalar(c: MPIntervalContext, value: Any, label: str):
    """Accept exact scalars or directed intervals, never floats or mappings."""

    if isinstance(value, Mapping) or isinstance(value, bool) or isinstance(value, float):
        raise TypeError(label + " must use an exact scalar or MP interval")
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
    raise TypeError(label + " must use an exact scalar or MP interval")


def _coerce_row(c: MPIntervalContext, value: Any, label: str):
    _require(isinstance(value, Sequence) and not isinstance(value, (str, bytes, Mapping)),
             label + " must be a signed ordinary-Z pair")
    _require(len(value) == 2, label + " must contain exactly Z0 and Z1")
    rows = []
    for order, item in enumerate(value):
        _require(hasattr(item, "_mpi_"),
                 label + "[Z" + str(order) + "] must be a raw MP interval")
        rows.append(_interval(c, *_raw_bounds(item, label + "[Z" + str(order) + "]")))
    return tuple(rows)


def _records(value: Any) -> dict[str, Any]:
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


def _negate_interval(c: MPIntervalContext, value: Any):
    lower, upper = _raw_bounds(value)
    return _interval(c, _negate_tuple(upper), _negate_tuple(lower))


def _maxabs(c: MPIntervalContext, value: Any):
    """Select the larger endpoint magnitude using exact MPF tuple comparisons."""

    lower, upper = _raw_bounds(value)
    lower_abs = _strip_sign(lower)
    upper_abs = _strip_sign(upper)
    selected = lower_abs if mpf_cmp(lower_abs, upper_abs) >= 0 else upper_abs
    return _interval(c, selected, selected)


def _symmetric_radius(c: MPIntervalContext, bits: int):
    point = (0, 1, -bits, 1)
    radius = _interval(c, point, point)
    return _negate_interval(c, radius), radius


def _mapping_exact(mapping: Any, label: str) -> None:
    _require(isinstance(mapping, Mapping), label + " must be a chart mapping")
    _require(set(mapping) == set(CHARTS),
             label + " must have exactly first_switch, second_switch, post_power")


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
    checks["no_accepted_project_imports"] = not any(
        (node.module or "").startswith("openai_ns_reconstruction")
        or (node.module or "").startswith(".original_")
        for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)
    )
    result = {
        "passed": all(checks.values()),
        "path": str(path),
        "module_sha256": hashlib.sha256(raw).hexdigest(),
        "checks": checks,
        "ast_function_names": sorted(set(names)),
    }
    _require(result["passed"], "module LF/EOF/AST preflight failed")
    return result


def _symbolic_proof() -> dict[str, Any]:
    """Independently verify ordinary-Z differentiation with frozen weights."""

    import sympy as sp

    z, log_rm = sp.symbols("z log_rm", real=True)
    logs = sp.symbols("logW_first_switch logW_second_switch logW_post_power", real=True)
    functions = [sp.Function("f_" + chart)(z) for chart in CHARTS]
    coefficients = [sp.exp(log_w - 2 * log_rm) for log_w in logs]
    operator = sum(coefficient * function
                   for coefficient, function in zip(coefficients, functions))
    derivative = sum(coefficient * sp.diff(function, z)
                     for coefficient, function in zip(coefficients, functions))
    checks = {
        "ordinary_Z_derivative_linearity": sp.simplify(sp.diff(operator, z) - derivative) == 0,
        "all_weights_Z_independent": all(sp.diff(coefficient, z) == 0
                                          for coefficient in coefficients),
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "operator": "sum(exp(logW_chart - 2*logRm)*f_chart for chart in three charts)",
        "ordinary_derivative": "differentiate supplied rows only; weights are Z-independent",
        "sympy_exact": True,
    }


def transport(c: MPIntervalContext, signed_rows_by_chart: Mapping[str, Sequence[Any]],
              logW_by_chart: Mapping[str, Any], logRm: Any, radius_bits: int = 1024,
              geometry_Z_independent: bool = True) -> dict[str, Any]:
    """Transport three signed ordinary-Z rows using proved logarithmic weights."""

    _require(isinstance(c, MPIntervalContext), "c must be an MPIntervalContext")
    _require(geometry_Z_independent is True, "the pressure weight must be Z-independent")
    _require(type(radius_bits) is int and radius_bits > 0,
             "radius_bits must be a positive integer")
    _mapping_exact(signed_rows_by_chart, "signed_rows_by_chart")
    _mapping_exact(logW_by_chart, "logW_by_chart")
    rows = {chart: _coerce_row(c, signed_rows_by_chart[chart], chart)
            for chart in CHARTS}
    log_weights = {chart: _coerce_scalar(c, logW_by_chart[chart],
                                         "logW_by_chart[" + chart + "]")
                   for chart in CHARTS}
    log_rm = _coerce_scalar(c, logRm, "logRm")
    neg_radius, radius = _symmetric_radius(c, radius_bits)
    threshold = _negate_interval(c, c.mpf(radius_bits) * c.ln(c.mpf(2)))
    lower_threshold = _raw_bounds(threshold)[0]
    weighted_rows: dict[str, tuple[Any, Any]] = {}
    per_chart: dict[str, Any] = {}

    for chart in CHARTS:
        log_coefficient = log_weights[chart] - 2 * log_rm
        chart_rows = []
        row_metadata = []
        for order, row in enumerate(rows[chart]):
            source_atom = {
                "chart": chart,
                "ordinary_Z_order": order,
                "source_atom": chart + ".ordinaryZ" + str(order),
                "signed_function_enclosure": True,
                "geometry_Z_independent": True,
            }
            common = {
                "ordinary_Z_order": order,
                "source_atom_metadata": source_atom,
                "original_signed_row": _records(row),
                "original_logW": _records(log_weights[chart]),
                "log_coefficient": _records(log_coefficient),
            }
            if _is_zero(row):
                enclosure = c.mpf(0)
                bound = dict(common, exact_zero_input=True, emitted_kind="EXACT_ZERO",
                             log_absolute_bound=None, emitted_enclosure=_records(enclosure),
                             emitted_not_value_or_zero=False)
            else:
                magnitude = _maxabs(c, row)
                log_absolute = c.ln(magnitude)
                log_bound = log_absolute + log_coefficient
                upper_bound = _raw_bounds(log_bound)[1]
                comparison = mpf_cmp(upper_bound, lower_threshold)
                if comparison >= 0:
                    raise ValueError(
                        "suppression proof failed for " + chart + ".ordinaryZ" + str(order)
                        + ": exact upper log bound is not below the requested threshold; "
                        "no fallback cap or zero is permitted"
                    )
                enclosure = _interval(c, _raw_bounds(neg_radius)[0],
                                      _raw_bounds(radius)[1])
                bound = dict(common, exact_zero_input=False, emitted_kind="BOUND",
                             bound_label="symmetric_2^-radius_bits_enclosure",
                             maxabs_exact_endpoint_bound=_records(magnitude),
                             log_absolute_input_bound=_records(log_absolute),
                             log_absolute_bound=_records(log_bound),
                             threshold=_records(threshold),
                             upper_bound_vs_threshold_lower=comparison,
                             strict_upper_below_threshold=True,
                             emitted_enclosure=_records(enclosure),
                             emitted_not_value_or_zero=True)
            chart_rows.append(enclosure)
            row_metadata.append(bound)
        weighted_rows[chart] = tuple(chart_rows)
        per_chart[chart] = {
            "weight_formula": "exp(logW_" + chart + " - 2*logRm)",
            "weight_materialized": False,
            "original_logW": _records(log_weights[chart]),
            "log_coefficient": _records(log_coefficient),
            "rows": row_metadata,
        }

    combined = tuple(sum((weighted_rows[chart][order] for chart in CHARTS), c.mpf(0))
                     for order in ORDERS)
    metadata = {
        "schema": "original_pressure_prefix_log_transport_v1",
        "normalization": "logarithmic_weight_to_Rm2",
        "charts": list(CHARTS),
        "ordinary_Z_orders": list(ORDERS),
        "geometry_Z_independent": True,
        "radius_bits": radius_bits,
        "radius": _records(radius),
        "logRm": _records(log_rm),
        "threshold": _records(threshold),
        "per_chart": per_chart,
        "combined_rows": [_records(row) for row in combined],
        "original_source_inputs": {
            chart: [_records(row) for row in rows[chart]] for chart in CHARTS
        },
        "source_contract": {
            "signed_ordinary_Z_inputs": True,
            "input_rows_are_cap_mappings": False,
            "input_rows_are_float_values": False,
            "source_rows_recomputed": False,
            "accepted_project_math_imported": False,
            "unmaterializable_weight_exponentiated": False,
            "source_half_S1_factor_applied": False,
        },
    }
    result = {
        "passed": True,
        "rawIVresult": combined,
        "raw_iv_result": combined,
        "metadata": metadata,
        "auditablemetadata": metadata,
    }
    json.dumps(metadata, allow_nan=False)
    return result


def fresh_fixture(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Build an independent rational fixture after source preflight."""

    preflight = _module_source_preflight(repo_root)
    proof = _symbolic_proof()
    old_prec = mp.mp.prec
    checks: dict[str, bool] = {
        "source_lf_eof_ast_preflight": preflight["passed"],
        "fresh_symbolic_ordinary_Z_derivative": proof["passed"],
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
        positive = iv(Fraction(1, 5), Fraction(1, 4))
        negative = iv(Fraction(-3, 5), Fraction(-1, 3))
        zero = c.mpf(0)
        rows = {
            "first_switch": (negative, positive),
            "second_switch": (positive, negative),
            "post_power": (zero, negative),
        }
        logs = {chart: Fraction(-chart_index - 1, 10)
                for chart_index, chart in enumerate(CHARTS)}
        result = transport(c, rows, logs, Fraction(1000), radius_bits=128)
        metadata = result["metadata"]
        checks.update({
            "global_mp53_tuple_safe": all(
                item["maxabs_exact_endpoint_bound"]["upper_exact_mpf_tuple"][1] != 0
                for chart in metadata["per_chart"].values()
                for item in chart["rows"] if not item["exact_zero_input"]
            ),
            "mixed_signed_rows_accepted": True,
            "exact_zero_retained": metadata["per_chart"]["post_power"]["rows"][0]["emitted_kind"] == "EXACT_ZERO",
            "all_nonzero_rows_labelled_bounds": all(
                item["emitted_kind"] == "BOUND" and item["emitted_not_value_or_zero"]
                for chart in metadata["per_chart"].values()
                for item in chart["rows"] if not item["exact_zero_input"]
            ),
            "weights_not_materialized": all(
                chart["weight_materialized"] is False
                for chart in metadata["per_chart"].values()
            ),
            "combined_rows_are_two": len(result["rawIVresult"]) == 2,
            "source_half_S1_not_applied": metadata["source_contract"]["source_half_S1_factor_applied"] is False,
        })
        try:
            transport(c, rows, logs, Fraction(1000), radius_bits=4096)
        except ValueError as error:
            checks["guard_threshold_rejects"] = "suppression proof failed" in str(error)
        else:
            checks["guard_threshold_rejects"] = False
        case = {"rawIVresult": [_records(row) for row in result["rawIVresult"]],
                "metadata": metadata}
    finally:
        mp.mp.prec = old_prec
    output: dict[str, Any] = {
        "passed": all(checks.values()),
        "schema": "pressure_prefix_log_transport_rational_fixture_v1",
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
        },
    }
    payload = json.dumps(output, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    output["fixture_sha256"] = hashlib.sha256(payload).hexdigest()
    json.dumps(output, allow_nan=False)
    _require(output["passed"], "fresh pressure prefix log transport fixture failed")
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
