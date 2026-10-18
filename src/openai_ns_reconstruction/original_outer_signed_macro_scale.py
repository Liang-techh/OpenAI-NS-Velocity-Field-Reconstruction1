"""Signed outer ``B``/``S1`` assembly at the exact macro scale.

The two ordinary C1 jets and the two signed D rows are caller supplied.  The
kernel only transports the micro terms by their logarithmic factors, forms the
signed macro-scale rows, and projects the result to ``R_m^2``.  No source
producer or saved data is read and no exponential weight is materialized.
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


MODULE = "src/openai_ns_reconstruction/original_outer_signed_macro_scale.py"
FIXTURE_RELATIVE_PATH = "tmp/outer_signed_macro_scale_fixture.json"
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
             label + " must be an ordinary Z0/Z1 pair")
    _require(len(value) == 2, label + " must contain exactly Z0 and Z1")
    return tuple(_interval(c, *_raw_bounds(item, label + "[Z" + str(k) + "]"))
                 for k, item in enumerate(value))


def _record(value: Any) -> dict[str, Any]:
    lower, upper = _raw_bounds(value)
    return {"lower_exact_mpf_tuple": list(lower),
            "upper_exact_mpf_tuple": list(upper), "directed": True}


def _records(values: Sequence[Any]) -> list[dict[str, Any]]:
    return [_record(value) for value in values]


def _zero(value: Any) -> bool:
    lower, upper = _raw_bounds(value)
    return lower == _ZERO and upper == _ZERO


def _strip(row: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return _ZERO if row[1] == 0 else (0, row[1], row[2], row[3])


def _neg_tuple(row: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return _ZERO if row[1] == 0 else (1 - row[0], row[1], row[2], row[3])


def _neg(c: MPIntervalContext, value: Any):
    lower, upper = _raw_bounds(value)
    return _interval(c, _neg_tuple(upper), _neg_tuple(lower))


def _maxabs(c: MPIntervalContext, value: Any):
    lower, upper = _raw_bounds(value)
    left, right = _strip(lower), _strip(upper)
    chosen = left if mpf_cmp(left, right) >= 0 else right
    return _interval(c, chosen, chosen)


def _symmetric(c: MPIntervalContext, bits: int):
    radius = _interval(c, (0, 1, -bits, 1), (0, 1, -bits, 1))
    return _neg(c, radius), radius


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
        "unique_assemble": names.count("assemble") == 1,
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
    import sympy as sp

    z = sp.symbols("z", real=True)
    gm, ga, dm, da = (sp.Function(name)(z) for name in ("Gm", "Ga", "Dm", "Da"))
    wb, wd, wr = sp.symbols("wb wd wr", positive=True)
    b = ga - wb * gm
    s = b + da + wd * dm
    projected = wr * s
    checks = {
        "ordinary_Z_linear_frozen_weights": sp.simplify(
            sp.diff(projected, z) - wr * (sp.diff(ga, z) - wb * sp.diff(gm, z)
                                           + sp.diff(da, z) + wd * sp.diff(dm, z))) == 0,
        "B_weight_Z_independent": all(sp.diff(weight, z) == 0 for weight in (wb, wd, wr)),
    }
    return {"passed": all(checks.values()), "checks": checks,
            "formula": "B=G_macro-exp(2logRmicro-ln10000)*G_micro; "
                        "S1=B+D_macro+exp(logWmicro-ln10000)*D_micro",
            "ordinary_derivative": "differentiate supplied C1 rows only; weights frozen",
            "sympy_exact": True}


def _threshold(c: MPIntervalContext, bits: int):
    return _neg(c, c.mpf(bits) * c.ln(c.mpf(2)))


def _bound_log(c: MPIntervalContext, magnitude: Any, factor: Any,
               threshold: Any, label: str):
    log_bound = c.ln(magnitude) + factor
    comparison = mpf_cmp(_raw_bounds(log_bound)[1], _raw_bounds(threshold)[0])
    _require(comparison < 0, label + " is not strictly below threshold")
    return log_bound, comparison


def assemble(c: MPIntervalContext, G_micro: Sequence[Any], G_macro: Sequence[Any],
             D_micro: Sequence[Any], D_macro: Sequence[Any],
             logR_micro_upper: Any, logW_micro: Any, logRm: Any,
             radius_bits: int = 1024, macro_endpoint_identity: bool = True,
             geometry_Z_independent: bool = True) -> dict[str, Any]:
    """Assemble signed macro-scale ``B``/``S1`` and project to ``R_m^2``."""

    _require(isinstance(c, MPIntervalContext), "c must be an MPIntervalContext")
    _require(type(radius_bits) is int and radius_bits > 0,
             "radius_bits must be a positive integer")
    _require(macro_endpoint_identity is True,
             "exact macro endpoint identity R(m=1)=100 is required")
    _require(geometry_Z_independent is True,
             "all supplied scale factors must be Z-independent")
    gm, ga = _rows(c, G_micro, "G_micro"), _rows(c, G_macro, "G_macro")
    dm, da = _rows(c, D_micro, "D_micro"), _rows(c, D_macro, "D_macro")
    log_r_micro = _scalar(c, logR_micro_upper, "logR_micro_upper")
    log_w_micro = _scalar(c, logW_micro, "logW_micro")
    log_rm = _scalar(c, logRm, "logRm")
    log_w_macro = c.ln(c.mpf(10000))
    threshold = _threshold(c, radius_bits)
    neg_radius, radius = _symmetric(c, radius_bits)
    epsilon = _interval(c, _raw_bounds(neg_radius)[0], _raw_bounds(radius)[1])

    micro_b, b_rows = [], []
    micro_d, d_rows = [], []
    b_factor = c.mpf(2) * log_r_micro - log_w_macro
    d_factor = log_w_micro - log_w_macro
    for order in ORDERS:
        gm_row, dm_row = gm[order], dm[order]
        gm_abs, dm_abs = _maxabs(c, gm_row), _maxabs(c, dm_row)
        b_info = {"ordinary_Z_order": order, "original_signed_G": _record(gm_row),
                  "original_logR_micro_upper": _record(log_r_micro),
                  "log_factor": _record(b_factor), "coefficient_positive": True,
                  "coefficient_materialized": False, "formula": "exp(2logRmicro-ln10000)"}
        d_info = {"ordinary_Z_order": order, "original_signed_D": _record(dm_row),
                  "original_logW_micro": _record(log_w_micro),
                  "log_factor": _record(d_factor), "coefficient_positive": True,
                  "coefficient_materialized": False, "formula": "exp(logWmicro-ln10000)"}
        if _zero(gm_abs):
            micro_b.append(c.mpf(0))
            b_rows.append(dict(b_info, exact_zero_input=True, emitted_kind="EXACT_ZERO",
                               emitted_not_value_or_zero=False, log_bound=None))
        else:
            blog, bcmp = _bound_log(c, gm_abs, b_factor, threshold,
                                    "micro G/Wmacro transport")
            micro_b.append(epsilon)
            b_rows.append(dict(b_info, exact_zero_input=False,
                               maxabs_source=_record(gm_abs), log_bound=_record(blog),
                               upper_vs_threshold_lower=bcmp,
                               strict_upper_below_threshold=True, emitted_kind="BOUND",
                               bound_label="symmetric_2^-radius_bits_enclosure",
                               emitted_not_value_or_zero=True))
        if _zero(dm_abs):
            micro_d.append(c.mpf(0))
            d_rows.append(dict(d_info, exact_zero_input=True, emitted_kind="EXACT_ZERO",
                               emitted_not_value_or_zero=False, log_bound=None))
        else:
            dlog, dcmp = _bound_log(c, dm_abs, d_factor, threshold,
                                    "micro D/Wmacro transport")
            micro_d.append(epsilon)
            d_rows.append(dict(d_info, exact_zero_input=False,
                               maxabs_source=_record(dm_abs), log_bound=_record(dlog),
                               upper_vs_threshold_lower=dcmp,
                               strict_upper_below_threshold=True, emitted_kind="BOUND",
                               bound_label="symmetric_2^-radius_bits_enclosure",
                               emitted_not_value_or_zero=True))

    raw_b = tuple(ga[k] + _neg(c, micro_b[k]) for k in ORDERS)
    raw_s1 = tuple(raw_b[k] + da[k] + micro_d[k] for k in ORDERS)
    raw_result, s1_rows = [], []
    rm_factor = log_w_macro - c.mpf(2) * log_rm
    for order, row in enumerate(raw_s1):
        magnitude = _maxabs(c, row)
        info = {"ordinary_Z_order": order, "original_signed_S1_macro": _record(row),
                "log_factor": _record(rm_factor), "formula": "exp(ln10000-2logRm)"}
        if _zero(magnitude):
            raw_result.append(c.mpf(0))
            s1_rows.append(dict(info, exact_zero_input=True, emitted_kind="EXACT_ZERO",
                                emitted_not_value_or_zero=False, log_bound=None))
        else:
            slog, scmp = _bound_log(c, magnitude, rm_factor, threshold,
                                    "outer S1/Rm^2 transport")
            raw_result.append(epsilon)
            s1_rows.append(dict(info, exact_zero_input=False,
                                maxabs_source=_record(magnitude), log_bound=_record(slog),
                                upper_vs_threshold_lower=scmp,
                                strict_upper_below_threshold=True, emitted_kind="BOUND",
                                bound_label="symmetric_2^-radius_bits_enclosure",
                                emitted_not_value_or_zero=True))

    metadata = {
        "schema": "original_outer_signed_macro_scale_v1",
        "normalization": {"macro": "Wmacro=10000", "terminal": "Rm^2"},
        "ordinary_Z_orders": list(ORDERS), "geometry_Z_independent": True,
        "macro_endpoint_identity": True,
        "macro_endpoint_source_identity": "R(m=1)=100 exactly; Wmacro=R^2=10000 exactly",
        "macro_endpoint_numeric_cap_not_selected": True,
        "radius_bits": radius_bits, "threshold": _record(threshold),
        "symmetric_error_radius": _record(radius),
        "symmetric_error_enclosure_signed": True,
        "logRm": _record(log_rm), "logWmacro": _record(log_w_macro),
        "inputs": {"G_micro": _records(gm), "G_macro": _records(ga),
                   "D_micro": _records(dm), "D_macro": _records(da),
                   "logR_micro_upper": _record(log_r_micro),
                   "logW_micro": _record(log_w_micro)},
        "micro_G_transport": {"log_formula": "2logRmicro-ln10000+ln(maxabsGmicro)",
                               "rows": b_rows, "nonzero_emitted_as_bound": True},
        "outer_B_macro_scale": {"formula": "G_macro-micro_G_BOUND",
                                 "rows": _records(raw_b), "signed": True,
                                 "orientation": {"macro": 1, "micro": -1}},
        "micro_D_transport": {"log_formula": "logWmicro-ln10000+ln(maxabsDmicro)",
                               "rows": d_rows, "nonzero_emitted_as_bound": True},
        "outer_S1_macro_scale": {"formula": "B_macro_scale+D_macro+micro_D_BOUND",
                                  "rows": _records(raw_s1), "signed": True,
                                  "pressure_half_factor_applied": False},
        "outer_S1_Rm2_transport": {"log_formula": "ln(maxabsS1)+ln10000-2logRm",
                                    "rows": s1_rows},
        "rawIVresult_is_outerS1normRm2": True,
        "source_contract": {"signed_inputs_preserved": True,
                             "micro_terms_are_bounds_not_selected_zero": True,
                             "weights_unmaterialized": True,
                             "no_source_rows_recomputed": True,
                             "no_half_factor": True,
                             "no_endpoint_jacobian": True,
                             "accepted_project_math_imported": False,
                             "no_huge_exponentiation": True},
        "claims": {"outer_S1_bound_only": True, "complete_pressure_closure": False},
    }
    json.dumps(metadata, allow_nan=False)
    return {"passed": True, "rawB_macro_scale": raw_b,
            "rawS1_macro_scale": raw_s1, "rawIVresult": tuple(raw_result),
            "metadata": metadata, "auditablemetadata": metadata}


def fresh_fixture(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Run the source guard, exact proof, and fresh synthetic signed cases."""

    preflight = _module_source_preflight(repo_root)
    proof = _symbolic_proof()
    old_precision = mp.mp.prec
    checks = {"source_lf_eof_ast_preflight": preflight["passed"],
              "fresh_symbolic_frozen_weight_linear_identity": proof["passed"]}
    case: dict[str, Any] = {}
    try:
        mp.mp.prec = 53
        c = MPIntervalContext()
        c.prec = 256

        def iv(a: Fraction, b: Fraction):
            return c.mpf([c.mpf(a.numerator) / a.denominator,
                          c.mpf(b.numerator) / b.denominator])

        gm = (iv(Fraction(-2), Fraction(1)), iv(Fraction(-1, 2), Fraction(1, 3)))
        ga = (iv(Fraction(-1, 3), Fraction(1, 4)), iv(Fraction(1, 5), Fraction(1, 2)))
        dm = (iv(Fraction(-1, 4), Fraction(1, 5)), iv(Fraction(-1, 10), Fraction(1, 10)))
        da = (iv(Fraction(1, 7), Fraction(1, 6)), iv(Fraction(-1, 8), Fraction(-1, 9)))
        result = assemble(c, gm, ga, dm, da, Fraction(-1000), Fraction(-2000),
                          Fraction(1000), radius_bits=128)
        meta = result["metadata"]
        checks.update({
            "global_mp53_with_context_256": mp.mp.prec == 53 and c.prec == 256,
            "negative_weighted_micro_source_covers_sign": (
                meta["inputs"]["G_micro"][0]["lower_exact_mpf_tuple"][0] == 1
                and meta["micro_G_transport"]["rows"][0]["emitted_kind"] == "BOUND"
                and meta["micro_G_transport"]["rows"][0]["emitted_not_value_or_zero"]),
            "signed_B_macro_rows": (
                meta["outer_B_macro_scale"]["signed"]
                and meta["outer_B_macro_scale"]["rows"][0]["lower_exact_mpf_tuple"][0] == 1),
            "signed_S1_macro_rows": (
                meta["outer_S1_macro_scale"]["signed"]
                and any(meta["outer_S1_macro_scale"]["rows"][order]
                        ["lower_exact_mpf_tuple"][0] == 1
                        for order in ORDERS)),
            "micro_nonzero_rows_are_bounds": all(
                row["emitted_kind"] == "BOUND" and row["emitted_not_value_or_zero"]
                for row in meta["micro_G_transport"]["rows"] + meta["micro_D_transport"]["rows"]),
            "rm_rows_are_two": len(result["rawIVresult"]) == 2,
            "no_half_factor": meta["outer_S1_macro_scale"]["pressure_half_factor_applied"] is False,
            "exact_macro_identity_retained": meta["macro_endpoint_numeric_cap_not_selected"]
            and meta["macro_endpoint_identity"],
        })
        try:
            assemble(c, gm, ga, dm, da, Fraction(-1000), Fraction(-2000),
                     Fraction(0), radius_bits=128)
        except ValueError as error:
            checks["threshold_guard_rejects"] = "threshold" in str(error)
        else:
            checks["threshold_guard_rejects"] = False
        zero = (c.mpf(0), c.mpf(0))
        zero_result = assemble(c, zero, zero, zero, zero, Fraction(-1000),
                               Fraction(-2000), Fraction(1000), radius_bits=128)
        checks["exact_zero_rows_may_zero"] = all(
            zero_result["metadata"][name]["rows"][order]["emitted_kind"] == "EXACT_ZERO"
            for name in ("micro_G_transport", "micro_D_transport", "outer_S1_Rm2_transport")
            for order in ORDERS)
        try:
            assemble(c, gm, ga, dm, da, Fraction(-1000), Fraction(-2000),
                     Fraction(1000), macro_endpoint_identity=False)
        except ValueError as error:
            checks["macro_identity_guard_rejects"] = "identity" in str(error)
        else:
            checks["macro_identity_guard_rejects"] = False
        try:
            assemble(c, gm, ga, dm, da, Fraction(-1000), Fraction(-2000),
                     Fraction(1000), geometry_Z_independent=False)
        except ValueError as error:
            checks["geometry_guard_rejects"] = "Z-independent" in str(error)
        else:
            checks["geometry_guard_rejects"] = False
        case = {"rawB_macro_scale": [_record(row) for row in result["rawB_macro_scale"]],
                "rawS1_macro_scale": [_record(row) for row in result["rawS1_macro_scale"]],
                "rawIVresult": [_record(row) for row in result["rawIVresult"]],
                "metadata": meta}
    finally:
        mp.mp.prec = old_precision
    output = {"passed": all(checks.values()),
              "schema": "outer_signed_macro_scale_rational_fixture_v1",
              "module_sha256": preflight["module_sha256"],
              "sourcechecks": preflight["checks"], "checks": checks,
              "exact_proof": proof, "case": case,
              "scope": {"source_executed": False, "saved_json_read": False,
                        "accepted_producer_replayed": False,
                        "actual_saved_data_produced": False,
                        "original_fixture_overwritten": False}}
    json.dumps(output, allow_nan=False)
    _require(output["passed"], "fresh outer signed macro scale fixture failed")
    return output


__all__ = ["FIXTURE_RELATIVE_PATH", "MODULE", "ORDERS", "assemble", "fresh_fixture"]
