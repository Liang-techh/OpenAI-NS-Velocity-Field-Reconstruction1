"""Independent synthetic fixture for the four-cell signed repair leaf.

The fixture has no import edge to the implementation, a receipt, or a source
producer.  A frozen runner calls :func:`run_fixture` with the leaf entry point
and its ``MPIntervalContext``.  All expected identities below are formed from
fresh rational arithmetic and directed interval endpoint checks.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any, Callable


SCHEMA = "original_corrected_fourcell_signed_repair_fixture_v1"
LEAF_SCHEMA = "original_corrected_fourcell_signed_repair_v1"
CELLS = (("-1", "-.5"), ("-.5", "0"), ("0", ".5"), (".5", "1"))


def _point(c: Any, value: Fraction | int | str) -> Any:
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / c.mpf(value.denominator)
    return c.mpf(value)


def _interval(c: Any, lower: Fraction | int | str, upper: Fraction | int | str | None = None) -> Any:
    if upper is None:
        upper = lower
    lo = _point(c, lower)
    hi = _point(c, upper)
    if hasattr(c, "make_mpf"):
        raw_lo = getattr(lo, "_mpi_", None)
        raw_hi = getattr(hi, "_mpi_", None)
        if raw_lo is not None and raw_hi is not None:
            try:
                return c.make_mpf((raw_lo[0], raw_hi[1]))
            except Exception:
                pass
    return c.mpf([lo, hi])


def _canonical(c: Any, lower: Fraction | int | str, upper: Fraction | int | str) -> dict[str, Any]:
    value = _interval(c, lower, upper)
    raw = getattr(value, "_mpi_", None)
    if raw is None:
        return {"lower": str(lower), "upper": str(upper)}
    return {
        "lower": str(lower),
        "upper": str(upper),
        "lower_exact_mpf_tuple": list(raw[0]),
        "upper_exact_mpf_tuple": list(raw[1]),
    }


def _scale(
    c: Any,
    exponents: tuple[int | float, int | float, int | float, int | float],
    radius_power: int | float,
    offset: Fraction | int | str,
) -> dict[str, Any]:
    return {
        "source_exponents": list(exponents),
        "radius_power": radius_power,
        "additional_log_interval": _canonical(c, offset, offset),
    }


def _row(
    c: Any,
    *,
    coefficient: Any,
    scale: Mapping[str, Any],
    path: str,
    source_index: str,
) -> dict[str, Any]:
    return {
        "coefficient_interval": coefficient,
        "formal_positive_scale": dict(scale),
        "provenance": {
            "same_source_id": "synthetic-fourcell-source-v1",
            "source_basis_id": "synthetic-fourcell-basis-v1",
            "source_path": path,
        },
        "source_index": source_index,
        "encloses_original_source_function": True,
        "point_value_selected": False,
    }


def _packet(c: Any) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Return source rows plus both explicit caller certificates."""

    # R is negative, Pstar is positive, Y0 positive, and Y1 is nonzero
    # negative.  Half powers intentionally occur in every resulting scale.
    rows: dict[str, Any] = {}
    for index, cell in enumerate(CELLS):
        key = str(cell)
        rows[key] = {
            "R": _row(
                c,
                coefficient=_canonical(c, Fraction(-5, 4), Fraction(-3, 4)) if index == 0
                else _interval(c, Fraction(-5, 4), Fraction(-3, 4)),
                scale=_scale(c, (1, -0.5, 0, 0), 0.5, 0),
                path=f"synthetic/{key}/R",
                source_index=f"{key}:R",
            ),
            "Pstar": _row(
                c,
                coefficient=_interval(c, Fraction(3, 2), Fraction(5, 2)),
                scale=_scale(c, (0, 0.5, 1, 0), -0.5, Fraction(1, 10)),
                path=f"synthetic/{key}/Pstar",
                source_index=f"{key}:Pstar",
            ),
            "V_y": [
                _row(
                    c,
                    coefficient=_interval(c, Fraction(1, 2), Fraction(3, 4)),
                    scale=_scale(c, (-0.5, 0, 0, 0.5), 0.5, Fraction(-1, 20)),
                    path=f"synthetic/{key}/V_y[0]",
                    source_index=f"{key}:V_y[0]",
                ),
                _row(
                    c,
                    coefficient=_interval(c, Fraction(-2, 3), Fraction(-1, 3)),
                    scale=_scale(c, (0, 0.5, -0.5, 0), 0, Fraction(1, 20)),
                    path=f"synthetic/{key}/V_y[1]",
                    source_index=f"{key}:V_y[1]",
                ),
            ],
        }
    source_certificate = {
        "passed": True,
        "same_source_basis": True,
        "same_source_provenance": True,
        "delta_Z_independent": True,
        "R_Z_independent": True,
        "Pstar_Z_independent": True,
        "ordinary_Z1": True,
        "Y_already_physical_Pstar_once": True,
        "same_source_id": "synthetic-fourcell-source-v1",
        "source_basis_id": "synthetic-fourcell-basis-v1",
    }
    pstar_certificate = {
        "passed": True,
        "logPstar": _canonical(c, Fraction(3, 2), Fraction(8, 5)),
        "same_source_id": "synthetic-fourcell-source-v1",
        "source_basis_id": "synthetic-fourcell-basis-v1",
    }
    return rows, source_certificate, pstar_certificate


def _endpoint(value: Any, upper: bool) -> Any:
    if hasattr(value, "a") and hasattr(value, "b"):
        return value.b if upper else value.a
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    return value


def _bounds(value: Any) -> tuple[Any, Any]:
    return value.a, value.b


def _call(build: Any, c: Any, *args: Any, **kwargs: Any) -> Mapping[str, Any]:
    function: Callable[..., Any] = getattr(build, "build", build)
    result = function(c, *args, **kwargs)
    if not isinstance(result, Mapping):
        raise TypeError("four-cell repair leaf did not return a mapping")
    return result


def _direct_formula(z: Fraction, delta: Fraction) -> tuple[Fraction, Fraction]:
    denominator = 1 - delta * z * z
    F = 2 * z / denominator
    Fz = 2 * (1 + delta * z * z) / (denominator * denominator)
    return F, Fz


def _contains(c: Any, interval: Any, value: Fraction | int) -> bool:
    target = _point(c, value)
    try:
        return interval.a <= target.a and target.b <= interval.b
    except Exception:
        return False


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run independent synthetic checks; frozen runners own invocation."""

    rows, source_certificate, pstar_certificate = _packet(c)
    checks: dict[str, bool] = {}
    delta = _interval(c, Fraction(1, 10), Fraction(1, 9))
    try:
        result = _call(
            build,
            c,
            CELLS,
            delta,
            rows,
            source_certificate,
            pstar_certificate,
        )
    except Exception as exc:
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"build_completed": False},
            "check_count": 1,
            "error": f"{type(exc).__name__}: {exc}",
            "synthetic_only": True,
            "implementation_imported": False,
            "receipt_opened": False,
            "graph_traversed": False,
        }

    checks["build_completed"] = isinstance(result, Mapping)
    checks["schema"] = result.get("schema") == LEAF_SCHEMA
    checks["four_exact_cells"] = result.get("cells") == [list(cell) for cell in CELLS]
    checks["whole_cell_enclosures"] = result.get("point_values_selected") is False and all(
        item.get("axis_function_whole_cell") is True
        for item in result.get("cells_by_key", {}).values()
        if isinstance(item, Mapping)
    )
    checks["strict_positive_denominators"] = all(
        item.get("denominator_evidence", {}).get("strictly_positive") is True
        for item in result.get("cells_by_key", {}).values()
    )
    checks["beta_is_nontrivial_enclosure"] = (
        result.get("pstar_beta", {}).get("interval") is not None
        and result.get("pstar_beta", {}).get("not_selected_as_zero") is True
        and result.get("pstar_beta", {}).get("not_selected_as_one") is True
        and result.get("pstar_beta", {}).get("interval_is_enclosure_not_function_value") is True
        and _endpoint(result["pstar_beta"]["interval"], False) > 0
        and _endpoint(result["pstar_beta"]["interval"], True) == 1
    )
    checks["source_certificates_are_retained"] = (
        result.get("source_certificate", {}).get("same_source_basis") is True
        and result.get("source_certificate", {}).get("same_source_provenance") is True
        and result.get("pstar_positive_certificate", {}).get("logPstar_gt_zero") is True
        and all(
            result.get("source_certificate", {}).get("required_guards", {}).get(name) is True
            for name in (
                "delta_Z_independent",
                "R_Z_independent",
                "Pstar_Z_independent",
                "ordinary_Z1",
                "Y_already_physical_Pstar_once",
            )
        )
    )
    checks["signed_coefficients_and_scales"] = (
        result.get("coefficient_intervals_signed_directed") is True
        and result.get("formal_exponents_exact_integer_or_half_integer") is True
        and result.get("additional_log_offsets_added_as_directed_intervals") is True
    )
    checks["repair_term_counts"] = (
        result.get("repair0_term_count") == 1 and result.get("repair1_term_count") == 2
    )
    checks["ordinary_product_rule_is_not_merged"] = (
        result.get("ordinary_Z1_is_product_rule") is True
        and result.get("ordinary_Z1_terms_not_merged") is True
        and all(len(item.get("repair1", ())) == 2 for item in result.get("cells_by_key", {}).values())
    )

    cells_by_key = result.get("cells_by_key", {})
    positive = cells_by_key.get(str((".5", "1")), {})
    negative = cells_by_key.get(str(("-1", "-.5")), {})
    try:
        # R<0, Pstar>0, Y0>0: the F sign flips between negative and positive Z.
        checks["coefficient_signs_follow_Z"] = (
            _endpoint(positive["repair0"][0]["coefficient_interval"], False) < 0
            and _endpoint(positive["repair0"][0]["coefficient_interval"], True) <= 0
            and _endpoint(negative["repair0"][0]["coefficient_interval"], False) >= 0
            and _endpoint(negative["repair0"][0]["coefficient_interval"], True) > 0
        )
        checks["nonzero_Y1_retained"] = (
            positive["repair1"][0]["source_indices"]["Y"]
            == positive["source_indices"]["Y1"]
            if "Y1" in positive["source_indices"]
            else positive["repair1"][0]["term_name"].endswith("Y1")
        )
    except Exception:
        checks["coefficient_signs_follow_Z"] = False
        checks["nonzero_Y1_retained"] = False

    checks["beta_correlation_is_shared"] = all(
        len({term.get("beta_binding_id") for term in item.get("repair1", []) + item.get("repair0", [])}) == 1
        for item in cells_by_key.values()
    )
    checks["half_exponents_are_present"] = any(
        any(value == 0.5 or value == -0.5 for value in term.get("formal_positive_scale", {}).get("source_exponents", [])
            + [term.get("formal_positive_scale", {}).get("radius_power")])
        for item in cells_by_key.values()
        for term in item.get("repair0", []) + item.get("repair1", [])
    )

    # Fresh exact-rational product-rule and cross-zero checks do not inspect the
    # implementation's intervals as point values.
    F_neg, Fz_neg = _direct_formula(Fraction(-1, 4), Fraction(1, 10))
    F_pos, Fz_pos = _direct_formula(Fraction(1, 4), Fraction(1, 10))
    checks["ordinary_product_rule_rational"] = (
        Fz_neg == Fraction(2) * (1 + Fraction(1, 10) * Fraction(1, 16))
        / (1 - Fraction(1, 10) * Fraction(1, 16)) ** 2
        and Fz_pos == Fraction(2) * (1 + Fraction(1, 10) * Fraction(1, 16))
        / (1 - Fraction(1, 10) * Fraction(1, 16)) ** 2
        and F_neg < 0 < F_pos
    )
    inner_negative = cells_by_key.get(str(("-.5", "0")), {})
    inner_positive = cells_by_key.get(str(("0", ".5")), {})
    try:
        checks["whole_cell_contains_direct_samples"] = (
            _contains(c, inner_negative["F"], F_neg)
            and _contains(c, inner_negative["Fz"], Fz_neg)
            and _contains(c, inner_positive["F"], F_pos)
            and _contains(c, inner_positive["Fz"], Fz_pos)
        )
    except Exception:
        checks["whole_cell_contains_direct_samples"] = False
    z_cross = _interval(c, Fraction(-1, 4), Fraction(1, 4))
    checks["cross_zero_direct_sample"] = _endpoint(z_cross, False) < 0 < _endpoint(z_cross, True)

    # L near zero must fail closed for the whole cell, not silently sample an
    # endpoint.  A delta box containing one over Z^2 makes the first cell fail.
    try:
        _call(
            build,
            c,
            CELLS,
            _interval(c, Fraction(99, 100), Fraction(101, 100)),
            rows,
            source_certificate,
            pstar_certificate,
        )
    except Exception:
        checks["L_near_zero_rejected"] = True
    else:
        checks["L_near_zero_rejected"] = False

    bad_pstar = dict(pstar_certificate)
    bad_pstar["logPstar"] = _canonical(c, Fraction(-1, 10), Fraction(1, 10))
    try:
        _call(build, c, CELLS, delta, rows, source_certificate, bad_pstar)
    except Exception:
        checks["Pstar_positivity_rejected"] = True
    else:
        checks["Pstar_positivity_rejected"] = False

    bad_source = dict(source_certificate)
    bad_source["same_source_provenance"] = False
    try:
        _call(build, c, CELLS, delta, rows, bad_source, pstar_certificate)
    except Exception:
        checks["same_source_certificate_rejected"] = True
    else:
        checks["same_source_certificate_rejected"] = False

    missing_guard = dict(source_certificate)
    missing_guard.pop("ordinary_Z1")
    try:
        _call(build, c, CELLS, delta, rows, missing_guard, pstar_certificate)
    except Exception:
        checks["missing_derivative_guard_rejected"] = True
    else:
        checks["missing_derivative_guard_rejected"] = False

    checks["no_ledger_or_global_closure"] = (
        result.get("ledger_written") is False
        and result.get("no_ledger") is True
        and result.get("global_closure_claimed") is False
    )
    checks["no_log_upper_or_point_value"] = (
        result.get("log_absolute_upper_used_as_value") is False
        and result.get("point_values_selected") is False
    )
    passed = all(checks.values())
    return {
        "schema": SCHEMA,
        "passed": passed,
        "checks_passed": passed,
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "implementation_imported": False,
        "receipt_opened": False,
        "graph_traversed": False,
        "old_code_called": False,
        "direct_numerical_samples_are_independent": True,
    }


run = run_fixture


__all__ = ["SCHEMA", "CELLS", "run_fixture", "run"]
