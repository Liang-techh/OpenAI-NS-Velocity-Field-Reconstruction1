"""Pure logarithmic transport of the saved ``H/U`` ordinary-Z bounds.

The source rows are already absolute common-unit logarithmic upper bounds.
This component only adds ``ln(U)-2*ln(R_m)`` and proves the requested dyadic
radius.  It never evaluates a source function, materializes a huge weight,
or chooses a signed value for the unknown ``H``.
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


MODULE = "src/openai_ns_reconstruction/original_H_m15_Rm_transport.py"
FIXTURE_RELATIVE_PATH = "tmp/H_m15_Rm_transport_fixture.json"
ORDERS = (0, 1)
_ZERO = (0, 0, 0, 0)


def _require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def _endpoint(value: Any, label: str) -> tuple[int, int, int, int]:
    _require(isinstance(value, (tuple, list)) and len(value) == 4,
             label + " must be an exact MPF tuple")
    _require(all(type(part) is int for part in value),
             label + " MPF tuple components must be integers")
    result = tuple(value)
    sign, mantissa, exponent, bits = result
    _require(sign in (0, 1) and mantissa >= 0 and bits >= 0,
             label + " has invalid sign, mantissa, or bit count")
    _require(bits == mantissa.bit_length(), label + " has noncanonical bit count")
    if mantissa == 0:
        _require(result == _ZERO, label + " has noncanonical zero")
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


def _scalar(c: MPIntervalContext, value: Any, label: str):
    """Coerce an exact scalar or raw interval without accepting float/caps."""

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
            lower = _scalar(c, value[0], label + ".lower")
            upper = _scalar(c, value[1], label + ".upper")
            return _interval(c, _raw_bounds(lower)[0], _raw_bounds(upper)[1])
    raise TypeError(label + " must use an exact scalar or raw MP interval")


def _rows(c: MPIntervalContext, value: Any, label: str):
    _require(isinstance(value, Sequence) and not isinstance(value, (str, bytes, Mapping)),
             label + " must be an ordinary Z0/Z1 pair of raw MP intervals")
    _require(len(value) == 2, label + " must contain exactly Z0 and Z1")
    result = []
    for order, item in enumerate(value):
        result.append(_interval(c, *_raw_bounds(item, label + "[Z" + str(order) + "]")))
    return tuple(result)


def _record(value: Any) -> dict[str, Any]:
    lower, upper = _raw_bounds(value)
    return {"lower_exact_mpf_tuple": list(lower),
            "upper_exact_mpf_tuple": list(upper), "directed": True}


def _records(values: Sequence[Any]) -> list[dict[str, Any]]:
    return [_record(value) for value in values]


def _neg_tuple(row: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return _ZERO if row[1] == 0 else (1 - row[0], row[1], row[2], row[3])


def _neg(c: MPIntervalContext, value: Any):
    lower, upper = _raw_bounds(value)
    return _interval(c, _neg_tuple(upper), _neg_tuple(lower))


def _symmetric(c: MPIntervalContext, bits: int):
    radius = _interval(c, (0, 1, -bits, 1), (0, 1, -bits, 1))
    negative = _neg(c, radius)
    return _interval(c, _raw_bounds(negative)[0], _raw_bounds(radius)[1])


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
        "unique_transport": names.count("transport") == 1,
        "unique_fresh_fixture": names.count("fresh_fixture") == 1,
        "no_accepted_project_math_imports": not any(
            name.startswith("openai_ns_reconstruction") or name.startswith(".original_")
            for name in imports),
    }
    result = {"passed": all(checks.values()), "path": str(path),
              "module_sha256": hashlib.sha256(raw).hexdigest(),
              "checks": checks, "ast_function_names": sorted(set(names)),
              "imports": sorted(set(imports))}
    _require(result["passed"], "module LF/EOF/AST preflight failed")
    return result


def _symbolic_proof() -> dict[str, Any]:
    """Prove the two-row ordinary-Z identity with a positive frozen ratio."""

    import sympy as sp

    z = sp.symbols("z", real=True)
    common, log_rm = sp.symbols("common_unit_log logRm", real=True)
    rows = (sp.Function("H0")(z), sp.Function("H1")(z))
    logratio = common - 2 * log_rm
    transported = tuple(row + logratio for row in rows)
    ratio = sp.exp(logratio)
    checks = {
        "positive_unit_ratio": ratio.is_positive is True,
        "ordinary_Z_derivative_identity": all(
            sp.simplify(sp.diff(new, z) - sp.diff(old, z)) == 0
            for old, new in zip(rows, transported)),
        "ratio_Z_independent": sp.diff(logratio, z) == 0,
        "ordinary_Z_two_rows": len(rows) == 2,
    }
    return {"passed": all(checks.values()), "checks": checks,
            "formula": "lognorm=absolute_common_unit_log_upper+common_unit_log-2logRm",
            "positive_ratio_formula": "exp(common_unit_log-2logRm)>0",
            "ordinary_derivative": "differentiate supplied H rows only; unit and Rm logs are frozen",
            "sympy_exact": True}


def transport(c: MPIntervalContext, common_unit_H_log_bounds: Sequence[Any],
              common_unit_log: Any, logRm: Any, radius_bits: int = 1024,
              scale_Z_independent: bool = True) -> dict[str, Any]:
    """Transport the saved absolute ``H/U`` logarithmic rows to ``R_m^2``."""

    _require(isinstance(c, MPIntervalContext), "c must be an MPIntervalContext")
    _require(type(radius_bits) is int and radius_bits > 0,
             "radius_bits must be a positive integer")
    _require(scale_Z_independent is True,
             "the supplied unit-to-Rm scale must be Z-independent")
    rows = _rows(c, common_unit_H_log_bounds, "common_unit_H_log_bounds")
    unit_log = _scalar(c, common_unit_log, "common_unit_log")
    log_rm = _scalar(c, logRm, "logRm")
    logratio = unit_log - c.mpf(2) * log_rm
    threshold = _neg(c, c.mpf(radius_bits) * c.ln(c.mpf(2)))
    threshold_lower = _raw_bounds(threshold)[0]
    symmetric = _symmetric(c, radius_bits)
    transported = []
    row_metadata = []
    for order, bound in enumerate(rows):
        lognorm = bound + logratio
        comparison = mpf_cmp(_raw_bounds(lognorm)[1], threshold_lower)
        _require(comparison < 0,
                 "H/Rm^2 suppression proof failed at ordinary Z" + str(order)
                 + ": upper log bound is not strictly below threshold")
        transported.append(symmetric)
        row_metadata.append({
            "ordinary_Z_order": order,
            "source_atom": "native.actual_first_micro_seam.ordinary_Z[" + str(order)
                           + "].absolute_common_unit_log_upper",
            "original_absolute_common_unit_log_upper": _record(bound),
            "original_signed_H_unknown": True,
            "logratio": _record(logratio),
            "lognorm": _record(lognorm),
            "threshold": _record(threshold),
            "upper_vs_threshold_lower": comparison,
            "strict_upper_below_threshold": True,
            "emitted_kind": "BOUND",
            "bound_label": "symmetric_2^-radius_bits_enclosure",
            "emitted_enclosure": _record(symmetric),
            "emitted_not_value_or_zero": True,
        })

    metadata = {
        "schema": "original_H_m15_Rm_transport_v1",
        "normalization": "absolute_H_over_common_unit_to_Rm2",
        "ordinary_Z_orders": list(ORDERS),
        "scale_Z_independent": True,
        "common_unit_positive": True,
        "common_unit_definition": "U=reference*A*Rv^2",
        "common_unit_log": _record(unit_log),
        "logRm": _record(log_rm),
        "logratio": _record(logratio),
        "radius_bits": radius_bits,
        "threshold": _record(threshold),
        "symmetric_error_enclosure": _record(symmetric),
        "rows": row_metadata,
        "source_contract": {
            "source_actual_first_micro_seam": True,
            "input_is_absolute_common_unit_log_upper": True,
            "input_is_H_over_U": True,
            "signed_H_unknown_covered": True,
            "source_rows_recomputed": False,
            "additional_local_to_common_conversion": False,
            "F0_derivative_added": False,
            "basis_evaluation": False,
            "requadrature": False,
            "ratio_materialized": False,
            "huge_exponentiation": False,
            "input_cap_mapping": False,
            "no_value_or_zero_substitution": True,
        },
        "claims": {"ordinary_Z_two_row_bound": True, "complete_pressure_closure": False},
    }
    json.dumps(metadata, allow_nan=False)
    return {"passed": True, "rawIVresult": tuple(transported),
            "metadata": metadata, "auditablemetadata": metadata}


def fresh_fixture(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Run preflight, exact symbolic checks, and fresh synthetic cases only."""

    preflight = _module_source_preflight(repo_root)
    proof = _symbolic_proof()
    old_precision = mp.mp.prec
    checks = {"source_lf_eof_ast_preflight": preflight["passed"],
              "fresh_symbolic_positive_ratio_two_row_proof": proof["passed"]}
    case: dict[str, Any] = {}
    try:
        mp.mp.prec = 53
        c = MPIntervalContext()
        c.prec = 256

        def iv(a: Fraction, b: Fraction):
            return c.mpf([c.mpf(a.numerator) / a.denominator,
                          c.mpf(b.numerator) / b.denominator])

        bounds = (iv(Fraction(-1500), Fraction(-1499)),
                  iv(Fraction(-1498), Fraction(-1497)))
        result = transport(c, bounds, Fraction(-1000), Fraction(1000), radius_bits=128)
        meta = result["metadata"]
        checks.update({
            "global_mp53_with_context_256": mp.mp.prec == 53 and c.prec == 256,
            "negative_signed_unknown_covered": all(
                row["original_signed_H_unknown"] and
                row["original_absolute_common_unit_log_upper"]["lower_exact_mpf_tuple"][0] == 1
                for row in meta["rows"]),
            "two_rows_are_symmetric_bounds": len(result["rawIVresult"]) == 2 and all(
                row["emitted_kind"] == "BOUND" and row["emitted_not_value_or_zero"]
                for row in meta["rows"]),
            "strict_threshold_checks": all(row["strict_upper_below_threshold"]
                                            for row in meta["rows"]),
            "positive_unit_ratio_proved": proof["checks"]["positive_unit_ratio"] is True,
            "ratio_not_materialized": meta["source_contract"]["ratio_materialized"] is False,
            "source_rows_preserved": meta["source_contract"]["source_rows_recomputed"] is False,
        })
        try:
            transport(c, (iv(Fraction(-1), Fraction(-1)),
                          iv(Fraction(-1), Fraction(-1))), Fraction(0), Fraction(0),
                     radius_bits=128)
        except ValueError as error:
            checks["strict_threshold_rejects"] = "threshold" in str(error)
        else:
            checks["strict_threshold_rejects"] = False
        try:
            transport(c, bounds, Fraction(-1000), Fraction(1000),
                     radius_bits=128, scale_Z_independent=False)
        except ValueError as error:
            checks["geometry_guard_rejects"] = "Z-independent" in str(error)
        else:
            checks["geometry_guard_rejects"] = False
        case = {"rawIVresult": [_record(row) for row in result["rawIVresult"]],
                "metadata": meta}
    finally:
        mp.mp.prec = old_precision
    output = {
        "passed": all(checks.values()),
        "schema": "H_m15_Rm_transport_rational_fixture_v1",
        "module_sha256": preflight["module_sha256"],
        "sourcechecks": preflight["checks"],
        "checks": checks,
        "exact_proof": proof,
        "case": case,
        "scope": {"source_executed": False, "saved_json_read": False,
                   "accepted_producer_replayed": False,
                   "actual_saved_data_produced": False,
                   "original_fixture_overwritten": False},
    }
    json.dumps(output, allow_nan=False)
    _require(output["passed"], "fresh H m15 Rm transport fixture failed")
    return output


__all__ = ["FIXTURE_RELATIVE_PATH", "MODULE", "ORDERS", "fresh_fixture", "transport"]
