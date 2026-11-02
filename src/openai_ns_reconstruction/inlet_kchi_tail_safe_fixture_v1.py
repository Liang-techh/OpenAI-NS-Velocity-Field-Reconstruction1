"""Fresh independent checks for :mod:`inlet_kchi_tail_safe_v1`.

The fixture owns all synthetic intervals and provenance.  It receives a build
callable explicitly, never imports a producer or an earlier leaf, and does not
read saved data.  Its rational prefix checks are independent positive-series
audits of the all-order ``exp(A/4)`` and ordinary-Z insertion bounds.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Mapping


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
FIELDS = ("F", "Fx")
BRANCHES = ("F_from_F", "F_from_Fx", "Fx_from_F", "Fx_from_Fx")
SENTINEL_AXIS_ID = "fixture:axis:wholeZ:safe-v1"

__all__ = [
    "SOURCE_POWER_BASIS",
    "FIELDS",
    "BRANCHES",
    "SENTINEL_AXIS_ID",
    "run_fixture",
    "run",
]


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
    path: Any,
    error_class: str,
    error_only: bool = False,
    uncertainty: bool = True,
    signed_field: bool | None = None,
) -> dict[str, Any]:
    coefficient = value if hasattr(value, "_mpi_") else _q(c, value)
    result = {
        "coefficient": coefficient,
        "powers": list(powers),
        "source_ids": [source_id],
        "source_path": path,
        "error_classes": [error_class],
        "source_error_classes": [error_class],
        "uncertainty": uncertainty,
        "error_only": error_only,
        "symmetric_error": error_only,
        "original_term": {"source_ids": [source_id], "source_path": path, "tag": source_id},
    }
    if signed_field is not None:
        result["signed_field"] = signed_field
    return result


def _axis(c: Any, value: Fraction | int, *, path: Any, derivative: bool = False) -> dict[str, Any]:
    return {
        "coefficient": _q(c, value),
        "powers": [0, 0, 0, 0],
        "source_ids": [SENTINEL_AXIS_ID],
        "source_path": path,
        "error_classes": ["fixture_safe_axis_Z1" if derivative else "fixture_safe_axis_Z0"],
        "uncertainty": True,
        "whole_Z": True,
        "scope": "same-source axis Kchi whole-Z interval",
    }


def _bounds(value: Any) -> tuple[Any, Any]:
    return value.a, value.b


def _upper(value: Any) -> Any:
    return _bounds(value)[1]


def _nonzero(value: Any) -> bool:
    lower, upper = _bounds(value)
    return lower != 0 or upper != 0


def _symmetric(value: Any) -> bool:
    lower, upper = _bounds(value)
    return lower == -upper


def _all_terms(result: Mapping[str, Any]) -> list[dict[str, Any]]:
    terms: list[dict[str, Any]] = []
    for field in FIELDS:
        for label in ("Z0", "Z1"):
            terms.extend(result["tail"][field][label])
    return terms


def _factorial_int(n: int) -> int:
    result = 1
    for index in range(2, n + 1):
        result *= index
    return result


def _exp_upper_fraction(x: Fraction, terms: int = 160) -> Fraction:
    """Rational exp upper bound with a geometric remainder."""

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
    """Finite positive ``j=3+k`` prefix using exact factorial ratios."""

    total = Fraction(0)
    denominator_ratio = 1
    for k in range(terms + 1):
        if k:
            denominator_ratio *= (factorial_offset + 2 * k - 1) * (factorial_offset + 2 * k)
        weight = Fraction(3 + k, 3) if derivative else Fraction(1)
        total += weight * (A**k) / denominator_ratio
    return total


def _rejected(call: Any) -> bool:
    try:
        call()
    except (TypeError, ValueError, KeyError):
        return True
    return False


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run structural, mathematical, provenance, zero-edge and bad-input checks."""

    if hasattr(build, "build"):
        build = build.build
    if not callable(build):
        raise TypeError("build must be callable or expose build")

    inlet = {
        "F": {
            "Z0": [
                _term(c, Fraction(4, 9), "fixture:safe:F:Z0:positive", (1, 2, -1, 3), path=["G0g", "F", "Z0", 0], error_class="fixture_positive_F_Z0", signed_field=True),
                _term(c, Fraction(-5, 12), "fixture:safe:F:Z0:negative", (-2, 0, 1, -1), path=["G0g", "F", "Z0", 1], error_class="fixture_negative_F_Z0", signed_field=True),
            ],
            "Z1": [
                _term(c, Fraction(-7, 15), "fixture:safe:F:Z1:negative", (-1, 1, 0, 2), path=["G0g", "F", "Z1", 0], error_class="fixture_negative_F_Z1", signed_field=True),
                _term(c, _box(c, Fraction(-3, 19), Fraction(3, 19)), "fixture:safe:F:Z1:error", (0, -1, 2, 0), path=["G0g", "F", "Z1", 1], error_class="fixture_error_F_Z1", error_only=True),
            ],
        },
        "Fx": {
            "Z0": [
                _term(c, Fraction(5, 11), "fixture:safe:Fx:Z0:positive", (2, -1, 2, 1), path=["G0g", "Fx", "Z0", 0], error_class="fixture_positive_Fx_Z0", signed_field=True),
                _term(c, Fraction(-2, 13), "fixture:safe:Fx:Z0:negative", (0, 3, -2, 0), path=["G0g", "Fx", "Z0", 1], error_class="fixture_negative_Fx_Z0", signed_field=True),
            ],
            "Z1": [
                _term(c, Fraction(1, 17), "fixture:safe:Fx:Z1:positive", (-3, 0, 1, 2), path=["G0g", "Fx", "Z1", 0], error_class="fixture_positive_Fx_Z1", signed_field=True),
                _term(c, _box(c, Fraction(-1, 23), Fraction(1, 23)), "fixture:safe:Fx:Z1:error", (1, -2, 0, -1), path=["G0g", "Fx", "Z1", 1], error_class="fixture_error_Fx_Z1", error_only=True),
            ],
        },
    }
    chi = _axis(c, Fraction(-1, 3), path=["axis", "wholeZ", "chi"])
    chi_z = _axis(c, Fraction(2, 5), path=["axis", "wholeZ", "chi_Z"], derivative=True)
    sc = _q(c, Fraction(1, 6))
    result = build(c, inlet, chi, chi_z, sc)

    source_powers = {
        term["source_ids"][0]: tuple(term["powers"])
        for field in FIELDS
        for label in ("Z0", "Z1")
        for term in inlet[field][label]
    }
    source_terms = {
        term["source_ids"][0]: term
        for field in FIELDS
        for label in ("Z0", "Z1")
        for term in inlet[field][label]
    }
    checks: dict[str, bool] = {
        "schema": result.get("schema") == "inlet_kchi_tail_safe_v1",
        "status": result.get("status") == "conditional_formal_Hb_axis_Kchi_inlet_j_ge_3_tail",
        "fields": tuple(result.get("fields", ())) == FIELDS,
        "formal_basis": tuple(result.get("source_power_basis", ())) == SOURCE_POWER_BASIS,
        "branch_names": tuple(result.get("branch_names", ())) == BRANCHES,
        "rows_complete": all(len(result["tail"][field][label]) == (4 if label == "Z0" else 8) for field in FIELDS for label in ("Z0", "Z1")),
        "geometry_contract": (
            result["geometry"]["S_formula"] == "S=(sc/4)*exp(sc/16)"
            and result["geometry"]["L_bound"] == "L <= hb*S"
            and result["geometry"]["K0_bound"] == "0 <= K0(x,t) <= x-t"
            and result["geometry"]["K0x_bound"] == "0 <= K0_x(x,t) <= 1"
        ),
        "formal_hb_not_materialized": result["geometry"]["formal_hb_not_materialized"] is True and result["conditions"]["hb_materialized"] is False,
        "only_tail": result["tail_starts_at_j"] == 3 and result["finite_prefix"] == {field: {label: [] for label in ("Z0", "Z1")} for field in FIELDS},
        "safe_exp_majorant": (
            result["tail_majorant"]["large_A_safe"] is True
            and "exp(A/4)" in result["tail_majorant"]["F_value_formula"]
            and "sqrt(A)" not in result["tail_majorant"]["F_value_formula"]
        ),
        "derivative_majorant": "(1+A/12)*exp(A/4)" in result["tail_majorant"]["F_Z1_insertion_formula"],
        "factorial_proof_present": "4^k*k!" in result["tail_majorant"]["factorial_comparison"] and "(3+k)/3" in result["tail_majorant"]["insertion_weight"],
        "alternative_not_combined": result["alternative_decomposition"]["old_j2_tail_combined"] is False and result["alternative_decomposition"]["old_j2_error_subtracted"] is False,
        "not_ledger": result["scope"]["ledger_written"] is False and result["scope"]["receipt_reads"] is False,
    }

    all_terms = _all_terms(result)
    required_term_keys = {
        "coefficient", "coefficient_interval", "powers", "formal_powers", "power_basis",
        "source_ids", "parent_source_ids", "inlet_source_ids", "chi_source_ids",
        "source_paths", "parent_source_paths", "inlet_source_paths", "error_classes",
        "source_error_classes", "error_only", "symmetric_error", "uncertainty", "signed_field",
        "tail", "tail_starts_at_j", "word_order", "field", "ordinary_Z_order",
        "formal_hb_power_added", "coefficient_contains_formal_hb", "signed_tail_value",
    }
    checks["all_term_fields"] = all(required_term_keys.issubset(term) for term in all_terms)
    checks["all_tail_error_only_symmetric"] = all(
        term["tail"] is True
        and term["error_only"] is True
        and term["symmetric_error"] is True
        and term["signed_field"] is False
        and _symmetric(term["coefficient"])
        for term in all_terms
    )
    checks["tail_word_order_three"] = all(term["word_order"] == 3 and term["tail_starts_at_j"] == 3 for term in all_terms)
    checks["formal_hb_shifts"] = all(
        term["powers"][0] == source_powers[term["inlet_source_ids"][0]][0] + term["formal_hb_power_added"]
        and term["powers"][1:] == list(source_powers[term["inlet_source_ids"][0]][1:])
        for term in all_terms
    )
    checks["branch_shift_values"] = all(
        result["branches"][term["branch"]]["formal_hb_power_added"]
        == {"F_from_F": 6, "F_from_Fx": 7, "Fx_from_F": 5, "Fx_from_Fx": 6}[term["branch"]]
        for term in all_terms
    )
    checks["source_paths_and_classes"] = all(
        term["inlet_source_paths"]
        and term["inlet_source_ids"]
        and term["inlet_source_ids"][0] in source_terms
        and source_terms[term["inlet_source_ids"][0]]["error_classes"][0] in term["error_classes"]
        for term in all_terms
    )
    checks["inherited_error_role_retained"] = any(
        term["inherited_error_only"] is True
        and term["parent_error_only"] is True
        and term["parent_signed_field"] is False
        and term["inlet_source_ids"][0] == "fixture:safe:F:Z1:error"
        for term in all_terms
    )
    checks["uncertainty_not_signed_inference"] = all(
        term["parent_error_only"] is False
        and term["parent_signed_field"] is True
        for term in all_terms
        if term["inlet_source_ids"][0] == "fixture:safe:F:Z0:positive"
    )
    z1_terms = result["tail"]["F"]["Z1"] + result["tail"]["Fx"]["Z1"]
    checks["ordinary_Z1_source_and_insertion"] = (
        any(term["word_kind"] == "j_ge_3_value_tail" and term["source_order"] == 1 for term in z1_terms)
        and any(term["word_kind"] == "j_ge_3_chi_Z_insertion_tail" and term["chi_Z_source_ids"] == [SENTINEL_AXIS_ID] for term in z1_terms)
        and all(term["ordinary_Z_order"] == 1 for term in z1_terms)
    )
    checks["chi_provenance"] = all(SENTINEL_AXIS_ID in term["chi_source_ids"] for term in all_terms)

    # Independent large-A rational prefixes: A=144, M=1/2 and S=12.
    probe_A = Fraction(144)
    probe_M = Fraction(1, 2)
    probe_S = Fraction(12)
    probe_Mz = Fraction(2, 5)
    probe_exp_upper = _exp_upper_fraction(probe_A / 4)
    probe_value_F_base = (2 * probe_M) ** 3 * probe_S**6 / _factorial_int(6)
    probe_value_Fx_base = (2 * probe_M) ** 3 * probe_S**5 / _factorial_int(5)
    probe_insert_F_base = 3 * (2 * probe_Mz) * (2 * probe_M) ** 2 * probe_S**6 / _factorial_int(6)
    probe_insert_Fx_base = 3 * (2 * probe_Mz) * (2 * probe_M) ** 2 * probe_S**5 / _factorial_int(5)
    value_F_prefix = probe_value_F_base * _word_series_fraction(probe_A, 6, 24)
    value_Fx_prefix = probe_value_Fx_base * _word_series_fraction(probe_A, 5, 24)
    insert_F_prefix = probe_insert_F_base * _word_series_fraction(probe_A, 6, 24, derivative=True)
    insert_Fx_prefix = probe_insert_Fx_base * _word_series_fraction(probe_A, 5, 24, derivative=True)
    checks["large_A_F_upper"] = probe_value_F_base * probe_exp_upper >= value_F_prefix
    checks["large_A_Fx_upper"] = probe_value_Fx_base * probe_exp_upper >= value_Fx_prefix
    probe_derivative_upper = (1 + probe_A / 12) * probe_exp_upper
    checks["large_A_F_derivative_upper"] = probe_insert_F_base * probe_derivative_upper >= insert_F_prefix
    checks["large_A_Fx_derivative_upper"] = probe_insert_Fx_base * probe_derivative_upper >= insert_Fx_prefix

    small_chi = dict(chi)
    small_chi["coefficient"] = _q(c, Fraction(1, 1000))
    small_result = build(c, inlet, small_chi, chi_z, sc)
    checks["small_chi_nonzero_tail"] = any(_nonzero(term["coefficient"]) for term in _all_terms(small_result))

    zero_chi = dict(chi)
    zero_chi["coefficient"] = _q(c, 0)
    zero_result = build(c, inlet, zero_chi, chi_z, sc)
    checks["zero_chi_nonzero_chi_Z"] = (
        not _nonzero(zero_result["chi"]["coefficient"])
        and _nonzero(zero_result["chi_z"]["coefficient"])
        and all(not _nonzero(term["coefficient"]) for term in _all_terms(zero_result))
    )

    # The source packet has a signed negative interval, an inherited symmetric
    # error, and a positive/negative ordinary-Z pair; none is merged.
    checks["positive_negative_inputs"] = (
        _upper(inlet["F"]["Z0"][0]["coefficient"]) > 0
        and _upper(inlet["F"]["Z0"][1]["coefficient"]) < 0
        and any(term["inlet_source_ids"][0] == "fixture:safe:F:Z0:negative" for term in all_terms)
    )

    missing_fx = {"F": inlet["F"]}
    bad_reverse = {"F": dict(inlet["F"], Z0=[dict(inlet["F"]["Z0"][0], coefficient={"lower": 1, "upper": -1})]), "Fx": inlet["Fx"]}
    bad_asymmetric = {"F": dict(inlet["F"], Z1=[dict(inlet["F"]["Z1"][1], coefficient=_box(c, Fraction(-1, 7), Fraction(2, 7)))]), "Fx": inlet["Fx"]}
    # Build this malformed term from scratch so no nested ``original_term`` or
    # metadata source ID can leak back into the leaf's provenance collector.
    bad_unknown_term = {
        "coefficient": inlet["F"]["Z0"][0]["coefficient"],
        "powers": list(inlet["F"]["Z0"][0]["powers"]),
        "source_ids": [],
        "source_path": [],
        "error_classes": [],
        "uncertainty": True,
        "original_term": {},
    }
    bad_unknown = {"F": dict(inlet["F"], Z0=[bad_unknown_term]), "Fx": inlet["Fx"]}
    bad_axis = dict(chi_z, source_ids=["fixture:axis:other"])
    bad_shear = dict(chi, scope="lower-collar shear")
    bad_chi_power = dict(chi, powers=[1, 0, 0, 0])
    checks["missing_fx_rejected"] = _rejected(lambda: build(c, missing_fx, chi, chi_z, sc))
    checks["reverse_interval_rejected"] = _rejected(lambda: build(c, bad_reverse, chi, chi_z, sc))
    checks["asymmetric_error_rejected"] = _rejected(lambda: build(c, bad_asymmetric, chi, chi_z, sc))
    checks["missing_source_rejected"] = _rejected(lambda: build(c, bad_unknown, chi, chi_z, sc))
    checks["different_axis_rejected"] = _rejected(lambda: build(c, inlet, chi, bad_axis, sc))
    checks["lower_shear_rejected"] = _rejected(lambda: build(c, inlet, bad_shear, chi_z, sc))
    checks["axis_power_rejected"] = _rejected(lambda: build(c, inlet, bad_chi_power, chi_z, sc))
    checks["nonpositive_sc_rejected"] = _rejected(lambda: build(c, inlet, chi, chi_z, _q(c, 0)))

    return {
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "independent_rational_sources": True,
        "large_A_bound_checked": True,
        "small_and_zero_chi_checked": True,
        "positive_negative_and_Z1_checked": True,
        "receipt_or_ancestor_read": False,
        "old_helpers_called": False,
        "midpoint_values_used": False,
        "cases": [
            "positive and negative grouped F/Fx intervals",
            "ordinary Z1 source and chi_Z insertion",
            "large-A rational j=3+k prefixes",
            "small and zero chi",
            "malformed interval/provenance/axis inputs",
        ],
    }


run = run_fixture
