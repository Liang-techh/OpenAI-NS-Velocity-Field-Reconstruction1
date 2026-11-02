"""Independent rational checks for :mod:`original_collar_g5_word2_v1`.

The fixture owns fresh synthetic rows and passes ``leaf.build`` explicitly via
``run_fixture``.  It does not import a producer, read saved data, or call an
older formal-axis/helper leaf.  Its return packet contains only JSON-friendly
scalars, strings, lists, and mappings; directed intervals remain inside the
runner's interval encoder.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Mapping


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
FIELDS = ("F", "Fx")
SENTINEL_AXIS_ID = "axis:wholeZ"

__all__ = ["SOURCE_POWER_BASIS", "FIELDS", "SENTINEL_AXIS_ID", "run_fixture", "run"]


def _q(c: Any, value: Fraction | int) -> Any:
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / value.denominator
    return c.mpf(value)


def _box(c: Any, lower: Fraction | int, upper: Fraction | int) -> Any:
    lo, hi = _q(c, lower), _q(c, upper)
    return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _term(
    c: Any,
    value: Fraction | int | Any,
    source_id: str,
    powers: tuple[int, int, int, int],
    *,
    path: str,
    error_class: str,
    origin: Any,
    error_only: bool = False,
) -> dict[str, Any]:
    coefficient = value if hasattr(value, "_mpi_") else _q(c, value)
    return {
        "coefficient": coefficient,
        "powers": list(powers),
        "source_ids": [source_id],
        "source_paths": [path],
        "error_classes": [error_class],
        "source_error_classes": [error_class],
        "origin": origin,
        "role": "symmetric_error_source" if error_only else "signed_raw_physical_g5",
        "uncertainty": error_only,
        "error_only": error_only,
        "symmetric_error": error_only,
        "signed_field": not error_only,
        "frame": "raw physical g5",
    }


def _axis(c: Any, value: Fraction | int, *, path: str) -> dict[str, Any]:
    return {
        "coefficient": _q(c, value),
        "powers": [0, 0, 0, 0],
        "source_ids": [SENTINEL_AXIS_ID],
        "source_paths": [path],
        "error_classes": ["fixture_axis_wholeZ"],
        "origin": {"operation": "fixture_whole_Z_axis", "source_ids": [SENTINEL_AXIS_ID]},
        "whole_Z": True,
        "certified_whole_Z": True,
        "scope": "same-source axis Kchi whole-Z interval",
    }


def _bounds(value: Any) -> tuple[Any, Any]:
    return value.a, value.b


def _nonzero(value: Any) -> bool:
    lower, upper = _bounds(value)
    return lower != 0 or upper != 0


def _contains_zero(value: Any) -> bool:
    lower, upper = _bounds(value)
    return lower <= 0 <= upper


def _all_terms(packet: Mapping[str, Any]) -> list[dict[str, Any]]:
    # Accept either one field packet (Z0/Z1) or the complete F/Fx packet.
    if "Z0" in packet and "Z1" in packet:
        return list(packet["Z0"]) + list(packet["Z1"])
    result: list[dict[str, Any]] = []
    for field in FIELDS:
        for label in ("Z0", "Z1"):
            result.extend(packet[field][label])
    return result


def _source_powers(g5: Mapping[str, Any]) -> dict[str, tuple[int, int, int, int]]:
    result: dict[str, tuple[int, int, int, int]] = {}
    for label in ("Z0", "Z1"):
        for term in g5[label]:
            result[term["source_ids"][0]] = tuple(term["powers"])
    return result


def _rejected(call: Any) -> bool:
    try:
        call()
    except (TypeError, ValueError, KeyError):
        return True
    return False


def _factorial_int(n: int) -> int:
    result = 1
    for index in range(2, n + 1):
        result *= index
    return result


def _exp_upper_fraction(x: Fraction, terms: int = 128) -> Fraction:
    """Independent rational upper bound for exp(x), with a geometric tail."""

    if x < 0 or x >= terms + 2:
        raise ValueError("probe exp bound requires 0 <= x < terms+2")
    partial = Fraction(0)
    term = Fraction(1)
    for index in range(terms + 1):
        if index:
            term *= x / index
        partial += term
    next_term = term * x / (terms + 1)
    ratio_cap = x / (terms + 2)
    return partial + next_term / (1 - ratio_cap)


def _word_series_fraction(A: Fraction, factorial_offset: int, terms: int, *, derivative: bool = False) -> Fraction:
    """Finite positive n=3+k series with exact factorial ratios."""

    total = Fraction(0)
    denominator_ratio = 1
    for k in range(terms + 1):
        if k:
            denominator_ratio *= (factorial_offset + 2 * k - 1) * (factorial_offset + 2 * k)
        weight = Fraction(3 + k, 3) if derivative else Fraction(1)
        total += weight * (A ** k) / denominator_ratio
    return total


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run structural, provenance, formal-power, and negative-path checks."""

    if hasattr(build, "build"):
        build = build.build
    if not callable(build):
        raise TypeError("build must be callable or expose build")

    g5 = {
        "Z0": [
            _term(c, Fraction(3, 7), "fixture:raw:g5:Z0:signed", (1, 2, -1, 3), path="raw.g5.Z0.signed", error_class="fixture_signed_Z0", origin={"tag": "signed_Z0"}),
            _term(c, _box(c, Fraction(-1, 11), Fraction(1, 11)), "fixture:raw:g5:Z0:error", (-2, 0, 1, -1), path="raw.g5.Z0.error", error_class="fixture_error_Z0", origin={"tag": "error_Z0"}, error_only=True),
        ],
        "Z1": [
            _term(c, Fraction(-5, 13), "fixture:raw:g5:Z1:signed", (-1, 1, 0, 2), path="raw.g5.Z1.signed", error_class="fixture_signed_Z1", origin={"tag": "signed_Z1"}),
            _term(c, _box(c, Fraction(-2, 17), Fraction(2, 17)), "fixture:raw:g5:Z1:error", (0, -1, 2, 0), path="raw.g5.Z1.error", error_class="fixture_error_Z1", origin={"tag": "error_Z1"}, error_only=True),
        ],
    }
    chi = _axis(c, Fraction(-1, 5), path="axis.wholeZ.chi")
    chi_z = _axis(c, Fraction(2, 7), path="axis.wholeZ.chi_Z")
    chi["uncertainty"] = True
    chi_z["uncertainty"] = True
    sc = _q(c, Fraction(1, 8))
    result = build(c, g5, chi, chi_z, sc)
    source_powers = _source_powers(g5)
    checks: dict[str, bool] = {
        "schema": result.get("schema") == "original_collar_g5_word2_v1",
        "fields": tuple(result.get("fields", ())) == FIELDS,
        "formal_basis": tuple(result.get("source_power_basis", ())) == SOURCE_POWER_BASIS,
        "rows_complete": len(result.get("g5_jets", {}).get("Z0", [])) == 2 and len(result.get("g5_jets", {}).get("Z1", [])) == 2,
        "geometry_contract": result.get("geometry", {}).get("S_formula") == "S=(sc/4)*exp(sc/16)" and result.get("geometry", {}).get("L_bound") == "L <= hb*S" and result.get("geometry", {}).get("K0x_bound") == "K0x <= 1",
        "formal_hb_not_materialized": result.get("scope", {}).get("formal_hb_not_materialized") is True and result.get("conditions", {}).get("hb_materialized") is False,
        "raw_frame_no_conversion": result.get("scope", {}).get("raw_physical_g5_source_frame") is True and result.get("scope", {}).get("F0_Lambda_conversion") is False,
        "axis_source_certified": result.get("conditions", {}).get("common_axis_source_operator") is True and SENTINEL_AXIS_ID in result.get("chi", {}).get("source_ids", []),
        "nonzero_chi_and_chi_Z": _nonzero(result["chi"]["coefficient"]) and _nonzero(result["chi_z"]["coefficient"]),
        "independent_large_A_majorant": result["tail_majorant"]["large_A_safe"] is True and "exp(A/4)" in result["tail_majorant"]["F_value_formula"] and "sqrt(A)/2" not in result["tail_majorant"]["F_value_formula"],
        "tail_insertion_majorant": "(1+A/12)*exp(A/4)" in result["tail_majorant"]["F_Z1_insertion_formula"],
        "alternative_not_combined": result["alternative_decomposition"]["prior_j2_tail_combined"] is False and result["alternative_decomposition"]["old_j2_error_subtracted"] is False,
        "not_ledger": result.get("scope", {}).get("ledger_written") is False and result.get("appended_to_212_ledger") is False,
    }

    # Independent large-A rational audit.  Set A=2*M*S^2=100 with M=1/2,
    # S=10, and compare exact positive n=3+k prefixes against a rational
    # upper for exp(A/4).  The derivative prefix uses the exact n/3 weight;
    # no production tail value or old helper participates in this check.
    probe_A = Fraction(100)
    probe_M = Fraction(1, 2)
    probe_S = Fraction(10)
    probe_Mz = Fraction(1, 3)
    probe_exp_upper = _exp_upper_fraction(probe_A / 4)
    probe_value_F_base = (2 * probe_M) ** 3 * probe_S ** 8 / _factorial_int(8)
    probe_value_Fx_base = (2 * probe_M) ** 3 * probe_S ** 7 / _factorial_int(7)
    probe_insert_F_base = 3 * (2 * probe_Mz) * (2 * probe_M) ** 2 * probe_S ** 8 / _factorial_int(8)
    probe_insert_Fx_base = 3 * (2 * probe_Mz) * (2 * probe_M) ** 2 * probe_S ** 7 / _factorial_int(7)
    value_F_prefix = probe_value_F_base * _word_series_fraction(probe_A, 8, 20)
    value_Fx_prefix = probe_value_Fx_base * _word_series_fraction(probe_A, 7, 20)
    insert_F_prefix = probe_insert_F_base * _word_series_fraction(probe_A, 8, 20, derivative=True)
    insert_Fx_prefix = probe_insert_Fx_base * _word_series_fraction(probe_A, 7, 20, derivative=True)
    checks["large_A_F_upper"] = probe_value_F_base * probe_exp_upper >= value_F_prefix
    checks["large_A_Fx_upper"] = probe_value_Fx_base * probe_exp_upper >= value_Fx_prefix
    probe_derivative_upper = (1 + probe_A / 12) * probe_exp_upper
    checks["large_A_F_derivative_upper"] = probe_insert_F_base * probe_derivative_upper >= insert_F_prefix
    checks["large_A_Fx_derivative_upper"] = probe_insert_Fx_base * probe_derivative_upper >= insert_Fx_prefix

    n2 = result["n2"]
    tail = result["tail"]
    for field, add, tail_add in (("F", 6, 8), ("Fx", 5, 7)):
        n2_terms = _all_terms(n2[field])
        tail_terms = _all_terms(tail[field])
        checks[field + ":n2_nonempty"] = bool(n2_terms)
        checks[field + ":tail_nonempty"] = bool(tail_terms)
        checks[field + ":n2_hb_increment"] = all(
            term["powers"][0] == source_powers[term["g5_source_ids"][0]][0] + add
            and term["powers"][1:] == list(source_powers[term["g5_source_ids"][0]][1:])
            for term in n2_terms
        )
        checks[field + ":tail_hb_increment"] = all(
            term["powers"][0] == source_powers[term["g5_source_ids"][0]][0] + tail_add
            and term["powers"][1:] == list(source_powers[term["g5_source_ids"][0]][1:])
            for term in tail_terms
        )
        checks[field + ":tail_error_only"] = all(term["error_only"] is True and term["symmetric_error"] is True and not term["signed_field"] for term in tail_terms)
        checks[field + ":source_paths_preserved"] = all(term.get("g5_source_paths") for term in n2_terms + tail_terms)

    z1_f = n2["F"]["Z1"]
    z1_fx = n2["Fx"]["Z1"]
    checks["F:Z1_two_routes"] = any(term.get("route") == "chi^2*g5_Z1" for term in z1_f) and any("chi_outer*chi_Z_inner" in term.get("route", "") for term in z1_f) and any("chi_inner*chi_Z_outer" in term.get("route", "") for term in z1_f)
    checks["Fx:Z1_two_routes"] = any(term.get("route") == "chi^2*g5_Z1" for term in z1_fx) and any("chi_outer*chi_Z_inner" in term.get("route", "") for term in z1_fx) and any("chi_inner*chi_Z_outer" in term.get("route", "") for term in z1_fx)
    checks["Z1_split_derivative_paths"] = all(term.get("insertion_copies_total") == 2 for term in z1_f + z1_fx if term.get("word_kind") == "n2_chi_Z_insertion")
    checks["error_role_retained"] = any(
        term.get("g5_source_ids") == ["fixture:raw:g5:Z0:error"]
        and term.get("error_only") is True
        and term.get("signed_field") is False
        and term.get("coefficient_is_signed_expression") is False
        and term.get("derived_role") != "signed_expression"
        for term in _all_terms(n2["F"])
    )
    checks["axis_uncertainty_flows"] = all(term.get("uncertainty") is True for term in _all_terms(n2["F"]) + _all_terms(n2["Fx"]))
    implicit_signed = dict(g5["Z0"][0])
    implicit_signed.pop("signed_field", None)
    implicit_rows = dict(g5)
    implicit_rows["Z0"] = [implicit_signed]
    implicit_result = build(c, implicit_rows, chi, chi_z, sc)
    checks["implicit_signedness_not_invented"] = all(
        term.get("signed_field") is False and term.get("coefficient_is_signed_expression") is False
        for term in _all_terms(implicit_result["n2"]["F"])
        if term.get("g5_source_ids") == ["fixture:raw:g5:Z0:signed"]
    )
    checks["tail_is_individual"] = all(term.get("source_term_index") in (0, 1) for term in _all_terms(tail["F"]) + _all_terms(tail["Fx"]))
    checks["tail_starts_at_three"] = result["tail_metadata"]["tail_starts_at_j"] == 3 and all(term.get("word_order") == 3 for term in _all_terms(tail["F"]) + _all_terms(tail["Fx"]))

    zero_chi = dict(chi)
    zero_chi["coefficient"] = _q(c, 0)
    zero_result = build(c, g5, zero_chi, chi_z, sc)
    checks["zero_chi_edge"] = (
        all(not _nonzero(term["coefficient"]) for term in _all_terms(zero_result["n2"]["F"]) + _all_terms(zero_result["n2"]["Fx"]))
        and all(not _nonzero(term["coefficient"]) for term in _all_terms(zero_result["tail"]["F"]) + _all_terms(zero_result["tail"]["Fx"]))
        and _nonzero(zero_result["chi_z"]["coefficient"])
    )

    bad_reverse = dict(g5)
    bad_reverse["Z0"] = [dict(g5["Z0"][0], coefficient={"lower": 1, "upper": -1})]
    bad_asymmetric = dict(g5)
    bad_asymmetric["Z0"] = [dict(g5["Z0"][1], coefficient=_box(c, Fraction(-1, 11), Fraction(2, 11)))]
    bad_unknown = dict(g5)
    bad_unknown["Z0"] = [dict(g5["Z0"][0], source_ids=["unknown:g5"])]
    bad_old = dict(g5)
    bad_old["Z0"] = [dict(g5["Z0"][0], error_classes=["old_j2_tail_error"])]
    bad_shear = dict(chi, scope="lower-collar shear")
    bad_axis = dict(chi_z, source_ids=["axis:other"])
    bad_axis["origin"] = {"operation": "bad_axis", "source_ids": ["axis:other"]}
    checks["reverse_interval_rejected"] = _rejected(lambda: build(c, bad_reverse, chi, chi_z, sc))
    checks["asymmetric_error_rejected"] = _rejected(lambda: build(c, bad_asymmetric, chi, chi_z, sc))
    checks["unknown_source_rejected"] = _rejected(lambda: build(c, bad_unknown, chi, chi_z, sc))
    checks["old_error_rejected"] = _rejected(lambda: build(c, bad_old, chi, chi_z, sc))
    checks["lower_shear_rejected"] = _rejected(lambda: build(c, g5, bad_shear, chi_z, sc))
    checks["different_axis_rejected"] = _rejected(lambda: build(c, g5, chi, bad_axis, sc))

    return {
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "independent_rational_sources": True,
        "large_A_bound_checked": True,
        "zero_chi_nonzero_chi_Z_checked": True,
        "receipt_or_ancestor_read": False,
        "old_helpers_called": False,
        "midpoint_values_used": False,
        "json_safe": True,
        "cases": ["n2 Z0/Z1 routes", "split chi_Z paths", "n>=3 large-A tail", "zero chi edge", "negative input guards"],
    }


run = run_fixture
