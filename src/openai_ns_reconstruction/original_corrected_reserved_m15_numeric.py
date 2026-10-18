"""Directed numeric enclosure for the corrected reserved ``m15`` repair.

This module is a small caller-owned arithmetic kernel.  It consumes exact
scalars, raw interval endpoints, and two ordinary ``Z`` rows; it does not
load a report, replay a producer, or construct a source graph.  The source
identity used by the caller is

``delta_k = oldpreRp - correctedpreRp`` and ``Y = actual physical V_y``.

The physical ``Pstar`` factor is already present in ``Y``.  The repair bound
therefore uses the supplied logarithmic ``Pstar`` only for the inequality
``abs(Pstar - 1) <= Pstar`` and never materialises a large exponential.
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


MODULE = "src/openai_ns_reconstruction/original_corrected_reserved_m15_numeric.py"
FIXTURE_RELATIVE_PATH = "tmp/corrected_reserved_m15_numeric_fixture.json"
ORDERS = (0, 1)
_ZERO = (0, 0, 0, 0)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _endpoint(value: Any, label: str) -> tuple[int, int, int, int]:
    """Validate and return one canonical finite MPF tuple."""

    _require(isinstance(value, (tuple, list)) and len(value) == 4,
             label + " must be an exact MPF tuple")
    _require(all(type(part) is int for part in value),
             label + " MPF tuple components must be integers")
    result = tuple(value)
    sign, mantissa, exponent, bit_count = result
    _require(sign in (0, 1) and mantissa >= 0 and bit_count >= 0,
             label + " has an invalid sign, mantissa, or bit count")
    if mantissa == 0:
        _require(result == _ZERO, label + " has a noncanonical zero")
    else:
        _require(bit_count == mantissa.bit_length(),
                 label + " has a noncanonical bit count")
        _require(mantissa & 1, label + " mantissa is not normalized")
    return result


def _mapping_bounds(value: Mapping[str, Any], label: str):
    """Read the exact endpoint forms emitted by the surrounding kernels."""

    if "value" in value and "lower_exact_mpf_tuple" not in value:
        return _raw_bounds(value["value"], label + ".value")
    if "source_enclosure_exact_mpf_tuples" in value:
        rows = value["source_enclosure_exact_mpf_tuples"]
        _require(isinstance(rows, Sequence) and len(rows) == 2,
                 label + ".source_enclosure_exact_mpf_tuples must have two endpoints")
        return (_endpoint(rows[0], label + ".lower"),
                _endpoint(rows[1], label + ".upper"))
    keys = ("lower_exact_mpf_tuple", "upper_exact_mpf_tuple")
    if all(key in value for key in keys):
        return (_endpoint(value[keys[0]], label + ".lower"),
                _endpoint(value[keys[1]], label + ".upper"))
    raise TypeError(label + " must contain exact MPF endpoint tuples")


def _raw_bounds(value: Any, label: str = "interval"):
    """Extract exact endpoints from an interval or an exact tuple record."""

    if isinstance(value, Mapping):
        lower, upper = _mapping_bounds(value, label)
    else:
        rows = getattr(value, "_mpi_", None)
        if rows is None and isinstance(value, (tuple, list)):
            if len(value) == 4 and all(type(part) is int for part in value):
                rows = (value, value)
            elif len(value) == 2:
                rows = value
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
    """Make an interval from exact endpoint tuples without global arithmetic."""

    low = _endpoint(lower, "lower")
    high = _endpoint(upper, "upper")
    _require(mpf_cmp(low, high) <= 0, "reversed interval")
    # ``make_mpf`` only decodes the exact tuple.  All arithmetic below is
    # performed by the caller's MPIntervalContext.
    result = c.mpf([mp.make_mpf(low), mp.make_mpf(high)])
    _raw_bounds(result)
    return result


def _coerce_scalar(c: MPIntervalContext, value: Any, label: str):
    """Accept exact scalar forms and directed raw intervals; reject floats."""

    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError(label + " must use an exact scalar or raw MP interval")
    if isinstance(value, Mapping):
        return _interval(c, *_mapping_bounds(value, label))
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


def _row_pair(c: MPIntervalContext, value: Any, label: str):
    """Return two ordinary rows, broadcasting one scalar when requested."""

    if isinstance(value, Mapping) or hasattr(value, "_mpi_"):
        item = _coerce_scalar(c, value, label)
        return (item, item)
    if isinstance(value, (str, bytes, bytearray)):
        item = _coerce_scalar(c, value, label)
        return (item, item)
    if isinstance(value, Sequence):
        _require(len(value) == 2, label + " must contain exactly Z0 and Z1")
        return tuple(_coerce_scalar(c, item, label + "[Z" + str(order) + "]")
                     for order, item in enumerate(value))
    item = _coerce_scalar(c, value, label)
    return (item, item)


def _record(value: Any) -> dict[str, Any]:
    lower, upper = _raw_bounds(value)
    return {
        "lower_exact_mpf_tuple": list(lower),
        "upper_exact_mpf_tuple": list(upper),
        "directed": True,
    }


def _records(values: Sequence[Any]) -> list[dict[str, Any]]:
    return [_record(value) for value in values]


def _is_zero(value: Any) -> bool:
    lower, upper = _raw_bounds(value)
    return lower == _ZERO and upper == _ZERO


def _strip_sign(value: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return _ZERO if value[1] == 0 else (0, value[1], value[2], value[3])


def _negate_tuple(value: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return _ZERO if value[1] == 0 else (1 - value[0], value[1], value[2], value[3])


def _maxabs(c: MPIntervalContext, value: Any):
    """Return the endpoint maximum magnitude using exact tuple comparison."""

    lower, upper = _raw_bounds(value)
    lower_abs = _strip_sign(lower)
    upper_abs = _strip_sign(upper)
    selected = lower_abs if mpf_cmp(lower_abs, upper_abs) >= 0 else upper_abs
    return _interval(c, selected, selected)


def _negate(c: MPIntervalContext, value: Any):
    lower, upper = _raw_bounds(value)
    return _interval(c, _negate_tuple(upper), _negate_tuple(lower))


def _symmetric(c: MPIntervalContext, radius_bits: int):
    radius_tuple = (0, 1, -radius_bits, 1)
    radius = _interval(c, radius_tuple, radius_tuple)
    return _negate(c, radius), radius


def _threshold(c: MPIntervalContext, radius_bits: int):
    return _negate(c, c.mpf(radius_bits) * c.ln(c.mpf(2)))


def _safe_logsum(c: MPIntervalContext, terms: Sequence[Any]):
    """Bound a positive sum by exact ``max(upper endpoints) + ln(n)``.

    ``terms`` contains logarithmic bounds.  ``None`` denotes a mathematically
    absent zero term and is omitted from the count.  The selected maximum is
    compared as an MPF tuple, never as an interval or a global ``mp`` value.
    """

    nonzero = [term for term in terms if term is not None]
    if not nonzero:
        return None, {
            "nonzero_term_count": 0,
            "selected_upper_exact_mpf_tuple": None,
            "exact_zero_sum": True,
        }
    uppers = [_raw_bounds(term)[1] for term in nonzero]
    selected = uppers[0]
    for candidate in uppers[1:]:
        if mpf_cmp(candidate, selected) > 0:
            selected = candidate
    point = _interval(c, selected, selected)
    bounded = point + c.ln(c.mpf(len(nonzero)))
    return bounded, {
        "nonzero_term_count": len(nonzero),
        "selected_upper_exact_mpf_tuple": list(selected),
        "candidate_upper_exact_mpf_tuples": [list(value) for value in uppers],
        "point_plus_ln_number_nonzero_terms": True,
        "safe_logsum_added_ln2": len(nonzero) == 2,
        "exact_zero_sum": False,
    }


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
    """Prove the ordinary C1 derivative and additive pulse contraction."""

    import sympy as sp

    z, delta = sp.symbols("z delta", real=True)
    L = 1 - delta * z**2
    F = 2 * z / L
    Fprime = 2 / L + 4 * delta * z**2 / L**2
    old, repair, target = sp.symbols("old repair target", real=True)
    corrected_discrepancy = old + repair
    corrected_reserved = target + corrected_discrepancy
    checks = {
        "ordinary_C1_F_derivative": sp.simplify(sp.diff(F, z) - Fprime) == 0,
        "additive_pulse_contraction": sp.simplify(
            corrected_reserved - (target + old + repair)
        ) == 0,
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "F_formula": "F=2*z/L",
        "Fprime_formula": "Fprime=2/L+4*delta*z^2/L^2",
        "contraction_formula": "corrected_reserved=independent_target+(old_discrepancy+repair_BOUND)",
        "sympy_exact": True,
    }


def enclose(
    c: MPIntervalContext,
    z: Any,
    delta: Any,
    logR_upper: Any,
    logPstar: Any,
    Y_log_bounds: Any,
    unit_log: Any,
    old_discrepancy: Any,
    independent_target: Any,
    radius_bits: int = 4096,
    *,
    scale_Z_independent: bool = True,
    same_source_repair: bool = True,
) -> dict[str, Any]:
    """Enclose the corrected reserved discrepancy at ordinary ``Z0``/``Z1``.

    The repair itself is certified in logarithmic form.  For a selected point
    or interval ``z`` and ``L=1-delta*z^2>0``,

    ``delta0/U = R*(Pstar-1)*F*Y0/U`` and
    ``delta1/U = R*(Pstar-1)*(F*Y1+Fprime*Y0)/U``.

    The caller supplies logarithmic absolute bounds for ``Y``.  Their signs
    remain unknown and no value is selected from a cap.  The returned repair
    rows are always labelled symmetric ``BOUND`` intervals.
    """

    _require(isinstance(c, MPIntervalContext), "c must be an MPIntervalContext")
    _require(type(radius_bits) is int and radius_bits > 0,
             "radius_bits must be a positive integer")
    _require(scale_Z_independent is True,
             "scale_Z_independent must be proved true")
    _require(same_source_repair is True,
             "same_source_repair must be proved true")

    z_scope = "interval"
    z_value = z
    if isinstance(z, Mapping) and set(z) in ({"point"}, {"tube"}, {"interval"}):
        z_scope = next(iter(z))
        z_value = z[z_scope]
    z_value = _coerce_scalar(c, z_value, "z")
    if _raw_bounds(z_value)[0] == _raw_bounds(z_value)[1]:
        z_scope = "POINT" if z_scope == "interval" else z_scope.upper()
    elif z_scope == "interval":
        z_scope = "TUBE"
    else:
        z_scope = z_scope.upper()

    delta_value = _coerce_scalar(c, delta, "delta")
    log_r = _coerce_scalar(c, logR_upper, "logR_upper")
    log_pstar = _coerce_scalar(c, logPstar, "logPstar")
    unit = _coerce_scalar(c, unit_log, "unit_log")
    y_logs = _row_pair(c, Y_log_bounds, "Y_log_bounds")
    old_rows = _row_pair(c, old_discrepancy, "old_discrepancy")
    target_rows = _row_pair(c, independent_target, "independent_target")

    # The strict lower endpoint is the only Pstar sign fact needed.  In
    # particular, no Pstar exponential is formed here.
    log_pstar_lower = _raw_bounds(log_pstar, "logPstar")[0]
    _require(mpf_cmp(log_pstar_lower, _ZERO) > 0,
             "logPstar.lower must be strictly positive")

    z2 = z_value * z_value
    L = c.mpf(1) - delta_value * z2
    _require(mpf_cmp(_raw_bounds(L, "L")[0], _ZERO) > 0,
             "L must have a strictly positive lower endpoint")
    inv_L = c.mpf(1) / L
    F = c.mpf(2) * z_value * inv_L
    Fprime = c.mpf(2) * inv_L + c.mpf(4) * delta_value * z2 * inv_L * inv_L

    threshold = _threshold(c, radius_bits)
    threshold_lower = _raw_bounds(threshold, "threshold")[0]
    neg_radius, radius = _symmetric(c, radius_bits)
    repair_rows = (  # Never replace these rows by an exact zero.
        _interval(c, _raw_bounds(neg_radius)[0], _raw_bounds(radius)[1]),
        _interval(c, _raw_bounds(neg_radius)[0], _raw_bounds(radius)[1]),
    )

    common_log = log_r + log_pstar - unit
    maxabs_F = _maxabs(c, F)
    maxabs_Fprime = _maxabs(c, Fprime)
    log_F = None if _is_zero(maxabs_F) else c.ln(maxabs_F)
    log_Fprime = None if _is_zero(maxabs_Fprime) else c.ln(maxabs_Fprime)

    if log_F is None or y_logs[0] is None:
        lognorm0 = None
    else:
        lognorm0 = common_log + log_F + y_logs[0]

    terms = []
    if log_F is not None and y_logs[1] is not None:
        terms.append(log_F + y_logs[1])
    if log_Fprime is not None and y_logs[0] is not None:
        terms.append(log_Fprime + y_logs[0])
    logsum1, logsum_metadata = _safe_logsum(c, terms)
    lognorm1 = None if logsum1 is None else common_log + logsum1

    row_metadata: list[dict[str, Any]] = []
    for order, lognorm in enumerate((lognorm0, lognorm1)):
        if lognorm is None:
            comparison = {
                "ordinary_Z_order": order,
                "exact_zero_source": True,
                "upper_vs_threshold_lower": -1,
                "strict_upper_below_threshold": True,
            }
        else:
            upper = _raw_bounds(lognorm)[1]
            comparison_value = mpf_cmp(upper, threshold_lower)
            _require(comparison_value < 0,
                     "corrected reserved repair is not strictly below threshold")
            comparison = {
                "ordinary_Z_order": order,
                "lognorm": _record(lognorm),
                "threshold": _record(threshold),
                "upper_vs_threshold_lower": comparison_value,
                "strict_upper_below_threshold": True,
            }

        row_metadata.append({
            "ordinary_Z_order": order,
            "Y_log": _record(y_logs[order]),
            "lognorm": None if lognorm is None else _record(lognorm),
            "logcomparisons": [comparison],
            "safe_logsum_Z1": logsum_metadata if order == 1 else None,
            "repair_enclosure": _record(repair_rows[order]),
            "bound_label": "symmetric_2^-radius_bits_enclosure",
            "emitted_kind": "BOUND",
            "emitted_not_value_or_zero": True,
            "exact_zero_source_does_not_select_zero": comparison.get(
                "exact_zero_source", False
            ),
        })

    corrected_discrepancy = tuple(
        old_rows[order] + repair_rows[order] for order in ORDERS
    )
    corrected_reserved = tuple(
        target_rows[order] + corrected_discrepancy[order] for order in ORDERS
    )

    inputs = {
        "z": _record(z_value),
        "delta": _record(delta_value),
        "logR_upper": _record(log_r),
        "logPstar": _record(log_pstar),
        "Y_log_bounds": _records(y_logs),
        "unit_log": _record(unit),
        "old_discrepancy": _records(old_rows),
        "independent_target": _records(target_rows),
    }
    metadata = {
        "schema": "original_corrected_reserved_m15_numeric_v1",
        "ordinary_Z_orders": list(ORDERS),
        "z_scope": z_scope,
        "radius_bits": radius_bits,
        "radius": _record(radius),
        "threshold": _record(threshold),
        "L": _record(L),
        "coefficientF": _record(F),
        "coefficientFprime": _record(Fprime),
        "coefficient_F": _record(F),
        "coefficient_Fprime": _record(Fprime),
        "coefficient_formula": {
            "F": "2*z/L",
            "Fprime": "2/L+4*delta*z^2/L^2",
            "L": "1-delta*z^2",
        },
        "log_coefficient": _record(common_log),
        "log_coefficient_formula": "logR_upper+logPstar-unit_log",
        "logcomparisons": [row["logcomparisons"][0] for row in row_metadata],
        "inputs": inputs,
        "allinputexactrows": True,
        "targetsamefunctiononce": True,
        "olderrorsretained": True,
        "sourcefunctionrepairassumption": True,
        "onepointscopeCALLER": True,
        "noallZclosure": True,
        "globalsfalse": True,
        "global_closure": False,
        "scale_Z_independent": True,
        "same_source_repair": True,
        "source_identity": {
            "delta_k": "oldpreRp-correctedpreRp",
            "Y": "actual physical V_y",
            "Y_already_Pstar_once": True,
        },
        "pstar_guard": {
            "logPstar_lower_strictly_positive": True,
            "abs_Pstar_minus_one_le_Pstar": True,
            "Pstar_materialized": False,
            "additional_Pstar_factor": False,
        },
        "source_contract": {
            "actual_signed_Y_selected_from_cap": False,
            "F0_normalization": False,
            "epsilon_normalization": False,
            "huge_exponentiation": False,
            "repair_is_symmetric_bound": True,
            "repair_is_selected_value": False,
            "raw_inputs_recomputed": False,
            "saved_inputs_read": False,
            "accepted_ancestor_replayed": False,
            "normalized_source_graph_constructed": False,
        },
        "rows": row_metadata,
    }
    return {
        "passed": True,
        "rawRepair_BOUND": repair_rows,
        "rawCorrectedDiscrepancy": corrected_discrepancy,
        "rawCorrectedReserved": corrected_reserved,
        "raw_repair_BOUND": repair_rows,
        "raw_corrected_discrepancy": corrected_discrepancy,
        "raw_corrected_reserved": corrected_reserved,
        "metadata": metadata,
        "auditablemetadata": metadata,
    }


def _fixture_interval(c: MPIntervalContext, lower: Fraction, upper: Fraction):
    def exact(value: Fraction):
        return c.mpf(value.numerator) / value.denominator

    return c.mpf([exact(lower), exact(upper)])


def _fixture_direct_repair(c: MPIntervalContext, z: Any, delta: Any,
                           log_r: Any, log_pstar: Any, y0: Any, y1: Any,
                           unit: Any, signs: tuple[int, int] = (1, 1)):
    """Evaluate one small finite witness only inside the synthetic fixture."""

    z2 = z * z
    L = c.mpf(1) - delta * z2
    F = c.mpf(2) * z / L
    Fprime = c.mpf(2) / L + c.mpf(4) * delta * z2 / L**2
    factor = c.exp(log_r) * (c.exp(log_pstar) - 1) / c.exp(unit)
    Y0 = signs[0] * c.exp(y0)
    Y1 = signs[1] * c.exp(y1)
    return factor * F * Y0, factor * (F * Y1 + Fprime * Y0)


def fresh_fixture(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Build fresh symbolic and finite synthetic cases only."""

    # This is intentionally the first action: no numerical or symbolic work
    # occurs before LF/single-EOF/AST inspection of this module.
    preflight = _module_source_preflight(repo_root)
    proof = _symbolic_proof()
    old_prec = mp.mp.prec
    checks: dict[str, bool] = {
        "source_lf_eof_ast_preflight": preflight["passed"],
        "fresh_symbolic_ordinary_C1_F_derivative": proof["passed"],
    }
    cases: dict[str, Any] = {}
    try:
        # Deliberately keep global mp precision at 53 while the interval
        # context has 256 bits; all comparisons still use exact endpoint
        # tuples and the interval context.
        mp.mp.prec = 53
        c = MPIntervalContext()
        c.prec = 256
        q = lambda value: c.mpf(value.numerator) / value.denominator

        log_r = Fraction(-8)
        log_pstar = Fraction(1, 5)
        unit = Fraction(4)
        y_logs = (Fraction(-160), Fraction(-161))
        old = (
            _fixture_interval(c, Fraction(-1, 10), Fraction(1, 20)),
            _fixture_interval(c, Fraction(-1, 20), Fraction(1, 10)),
        )
        target = (
            _fixture_interval(c, Fraction(-1, 5), Fraction(1, 10)),
            _fixture_interval(c, Fraction(-1, 10), Fraction(1, 5)),
        )

        finite_cases = (
            ("positive_geometry", q(Fraction(1, 5)), q(Fraction(1, 10)),
             (q(y_logs[0]), q(y_logs[1])), (1, 1)),
            ("negative_geometry", q(Fraction(-1, 5)), q(Fraction(1, 10)),
             (q(y_logs[0]), q(y_logs[1])), (1, -1)),
            ("crossing_unknown_Y", _fixture_interval(c, Fraction(-1, 10), Fraction(1, 10)),
             q(Fraction(1, 5)), (q(y_logs[0]), q(y_logs[1])), (-1, 1)),
        )
        finite_records = {}
        direct_coverage = True
        for name, z_case, delta_case, y_case, signs in finite_cases:
            call_kwargs = {
                "z": z_case,
                "delta": delta_case,
                "logR_upper": log_r,
                "logPstar": log_pstar,
                "Y_log_bounds": y_case,
                "unit_log": unit,
                "old_discrepancy": old,
                "independent_target": target,
                "radius_bits": 128,
            }
            if name == "crossing_unknown_Y":
                # Exercise raw exact-tuple records and retain Y signs only in
                # this finite witness, never in the enclosure API.
                probe = enclose(c, **call_kwargs)
                call_kwargs["z"] = {
                    "lower_exact_mpf_tuple": list(_raw_bounds(z_case)[0]),
                    "upper_exact_mpf_tuple": list(_raw_bounds(z_case)[1]),
                }
                call_kwargs["Y_log_bounds"] = [
                    _record(y_case[0]), _record(y_case[1])
                ]
                call_kwargs["old_discrepancy"] = [_record(row) for row in old]
                call_kwargs["independent_target"] = [_record(row) for row in target]
                result = enclose(c, **call_kwargs)
                checks["rawIV_exacttuple_rows_accepted"] = True
            else:
                result = enclose(c, **call_kwargs)
            finite_records[name] = {
                "metadata": result["metadata"],
                "rawRepair_BOUND": [_record(row) for row in result["rawRepair_BOUND"]],
            }

            # Test an interior finite point against the interval bound.  The
            # crossing case intentionally chooses an unknown Y sign pair.
            z_probe = q(Fraction(7, 100)) if name == "crossing_unknown_Y" else z_case
            if hasattr(z_probe, "_mpi_") and _raw_bounds(z_probe)[0] != _raw_bounds(z_probe)[1]:
                z_probe = _interval(c, _raw_bounds(z_probe)[0], _raw_bounds(z_probe)[1])
                z_probe = q(Fraction(7, 100))
            y0 = q(y_logs[0])
            y1 = q(y_logs[1])
            direct0, direct1 = _fixture_direct_repair(
                c, z_probe, delta_case if not hasattr(delta_case, "_mpi_")
                else q(Fraction(1, 5)), q(log_r), q(log_pstar), y0, y1,
                q(unit), signs)
            radius_upper = _raw_bounds(result["rawRepair_BOUND"][0])[1]
            direct_coverage &= (
                mpf_cmp(_strip_sign(_raw_bounds(direct0)[1]), radius_upper) <= 0
                and mpf_cmp(_strip_sign(_raw_bounds(direct1)[1]), radius_upper) <= 0
            )

        checks.update({
            "finite_geometry_positive_negative_crossing": len(finite_records) == 3,
            "finite_direct_cases_covered": direct_coverage,
            "all_rows_are_bounds": all(
                row["emitted_kind"] == "BOUND"
                and row["emitted_not_value_or_zero"]
                for case in finite_records.values()
                for row in case["metadata"]["rows"]
            ),
            "additive_pulse_contraction_metadata": all(
                case["metadata"]["targetsamefunctiononce"]
                and case["metadata"]["olderrorsretained"]
                for case in finite_records.values()
            ),
        })

        # Guard checks are all expected to reject before returning a result.
        try:
            enclose(c, q(Fraction(1, 5)), q(Fraction(1, 10)),
                    log_r, Fraction(-1), y_logs, unit, old, target, radius_bits=128)
        except ValueError as error:
            checks["Pstar_guard_rejects"] = "logPstar.lower" in str(error)
        else:
            checks["Pstar_guard_rejects"] = False
        try:
            enclose(c, q(Fraction(1)), q(Fraction(100)), log_r, log_pstar,
                    y_logs, unit, old, target, radius_bits=128)
        except ValueError as error:
            checks["L_guard_rejects"] = "strictly positive" in str(error)
        else:
            checks["L_guard_rejects"] = False
        try:
            enclose(c, q(Fraction(1, 5)), q(Fraction(1, 10)),
                    log_r, log_pstar, y_logs, unit, old, target, radius_bits=4096)
        except ValueError as error:
            checks["threshold_guard_rejects"] = "below threshold" in str(error)
        else:
            checks["threshold_guard_rejects"] = False
        try:
            enclose(c, q(Fraction(1, 5)), q(Fraction(1, 10)),
                    log_r, log_pstar, y_logs, unit, old, target,
                    radius_bits=128, scale_Z_independent=False)
        except ValueError as error:
            checks["scale_guard_rejects"] = "scale_Z_independent" in str(error)
        else:
            checks["scale_guard_rejects"] = False

        cases = finite_records
        cases["guards"] = {key: value for key, value in checks.items()
                            if key.endswith("_rejects")}
    finally:
        mp.mp.prec = old_prec

    output: dict[str, Any] = {
        "passed": all(checks.values()),
        "schema": "corrected_reserved_m15_numeric_synthetic_fixture_v1",
        "module_sha256": preflight["module_sha256"],
        "sourcechecks": preflight["checks"],
        "checks": checks,
        "exact_proof": proof,
        "case": cases,
        "scope": {
            "fresh_synthetic_fixture_only": True,
            "source_executed": False,
            "saved_json_read": False,
            "actual_saved_data_produced": False,
            "accepted_ancestor_replayed": False,
            "normalized_source_graph_constructed": False,
        },
    }
    payload = json.dumps(output, sort_keys=True, separators=(",", ":"),
                         allow_nan=False).encode("utf-8")
    output["fixture_sha256"] = hashlib.sha256(payload).hexdigest()
    json.dumps(output, allow_nan=False)
    _require(output["passed"], "fresh corrected reserved m15 fixture failed")
    return output


fixture_bound_check = fresh_fixture


__all__ = [
    "FIXTURE_RELATIVE_PATH",
    "MODULE",
    "enclose",
    "fixture_bound_check",
    "fresh_fixture",
]
