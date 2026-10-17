"""Evaluate the saved uniform raw first-chart all-word tail for N = 1.

The coupled first-step archive already contains the source-bound constants for
the raw Akin-preconditioned series.  This module consumes those records as
directed majorants and performs only the small N = 1 arithmetic

    Mbase * beta_1**(N + 1) / (1 - beta_1),  N = 1,

plus the saved ordinary-Z Cauchy factor.  The exact MPF tuples from the archive
are retained; decimal/mpmath values are display and comparison values only.
The normalized S^{-1} map is deliberately left unapplied because its
Z-derivative action is not a scalar division.  Field-prefix projections belong
to a separate source-bound stage and are not copied here.

No source owner, producer, recurrence, quadrature, or accepted numerical
fixture is imported or executed here.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

CUTOFF_REPORT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_cutoff_source_binding.json"
CUTOFF_REPORT_SHA256 = "4c36c2dc2e709d7ca740a689af9c9e16802569f300dfe9997851b9bb72142856"
COUPLED_ARCHIVE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step.json.gz"
COUPLED_ARCHIVE_SHA256 = "30b807af4147658ec0978ebd253cc1b2b31bed67bfed430a965e66461450add1"
COUPLED_SOURCE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step.py"
COUPLED_SOURCE_SHA256 = "97913b82b8570b220921cd89918569a4019ff3495bdbea0548b3c53e9d8c93e7"
COUPLED_CHECK = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step_check.json"
COUPLED_CHECK_SHA256 = "6a58caf1db41697620820830cbe7f76620874a8b429b351ad3d58f637e4e1416"
REMAINDER_MODULE = "src/openai_ns_reconstruction/original_n1_first_volterra_remainder.py"

INNER_ARCHIVE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.json.gz"
INNER_ARCHIVE_SHA256 = "1b44be8c6cbfffb692c7a0430ecbc8a590f46001a233c99c59ca8c5e2c44a6fe"
PARENT_REPORT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_collar_transport.json"
PARENT_REPORT_SHA256 = "7c05f7438a27d9caa4378cb89b9ab6d33b44138086e86292b73f2fd8e3f4bf96"

ORDINARY_Z_ORDERS = (0, 1)
WORD_DEPTH = 1
MP_WORK_PRECISION = 2300


def _repo_root(root: str | Path | None) -> Path:
    if root is None:
        return Path(__file__).resolve().parents[2]
    return Path(root).resolve()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _load_gzip_json(path: Path) -> Any:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def _pin(root: Path, relative: str, expected: str | None = None) -> dict[str, Any]:
    path = root / relative
    result: dict[str, Any] = {
        "path": relative,
        "exists": path.is_file(),
        "sha256": None,
        "expected_sha256": expected,
        "sha256_match": False,
        "source_executed": False,
    }
    if path.is_file():
        result["sha256"] = _sha256(path)
        result["sha256_match"] = expected is None or result["sha256"] == expected
    return result


def _mapping_digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    ).hexdigest()


def _contains_all(text: str, fragments: Sequence[str]) -> bool:
    return all(fragment in text for fragment in fragments)


def _tuple(record: Mapping[str, Any], side: str = "upper") -> list[int] | None:
    value = record.get(f"{side}_exact_mpf_tuple")
    if not isinstance(value, list) or len(value) != 3 and len(value) != 4:
        return None
    try:
        return [int(item) for item in value]
    except (TypeError, ValueError):
        return None


def _record_has_tuples(record: Any) -> bool:
    return isinstance(record, Mapping) and _tuple(record, "lower") is not None and _tuple(record, "upper") is not None


def _mp_value(record: Mapping[str, Any], side: str = "upper") -> Any:
    """Convert a saved MPF tuple without expanding its binary exponent."""
    import mpmath as mp

    value = _tuple(record, side)
    if value is None:
        raise ValueError(f"saved {side} exact MPF tuple is missing")
    if len(value) == 3:
        sign, mantissa, exponent = value
    else:
        sign, mantissa, exponent, _bits = value
    if mantissa < 0 or exponent != int(exponent):
        raise ValueError("invalid saved MPF tuple")
    result = mp.mpf(mantissa) * mp.power(2, int(exponent))
    return -result if sign else result


def _mp_tuple_value(value: Sequence[int]) -> Any:
    import mpmath as mp

    if len(value) not in (3, 4):
        raise ValueError("invalid MPF endpoint tuple")
    sign, mantissa, exponent = (int(value[0]), int(value[1]), int(value[2]))
    result = mp.mpf(mantissa) * mp.power(2, exponent)
    return -result if sign else result


def _iv_endpoint(value: Any, side: str = "upper") -> list[int]:
    endpoint = value._mpi_[1 if side == "upper" else 0]
    return [int(item) for item in endpoint]


def _iv_from_record(record: Mapping[str, Any], side: str = "upper") -> Any:
    import mpmath as mp

    endpoint = _tuple(record, side)
    if endpoint is None:
        raise ValueError("saved MPF endpoint tuple is missing")
    scalar = _mp_tuple_value(endpoint)
    return mp.iv.mpf([scalar, scalar])


def _log10_abs(record: Mapping[str, Any], side: str = "upper") -> str:
    import mpmath as mp

    with mp.workprec(MP_WORK_PRECISION):
        value = abs(_mp_value(record, side))
        if value == 0:
            return "-inf"
        return mp.nstr(mp.log10(value), 80)


def _interval_projection(record: Mapping[str, Any]) -> dict[str, Any]:
    """Keep directed display and exact endpoint tuples from a saved record."""
    if not _record_has_tuples(record):
        raise ValueError("selected saved interval lacks exact MPF endpoints")
    result: dict[str, Any] = {
        "lower_exact_mpf_tuple": list(_tuple(record, "lower") or ()),
        "upper_exact_mpf_tuple": list(_tuple(record, "upper") or ()),
        "directed_interval": record.get("directed_interval", True),
    }
    for side in ("lower", "upper"):
        if record.get(side) is not None:
            result[f"{side}_decimal_display"] = record[side]
        result[f"{side}_log10_abs"] = _log10_abs(record, side)
    if record.get("width") is not None:
        result["width_decimal_display"] = record["width"]
    return result


def _upper_value(record: Mapping[str, Any]) -> Any:
    return _mp_value(record, "upper")


def _source_parent_gate(parent: Any, cutoff: Mapping[str, Any]) -> tuple[bool, str | None]:
    if parent is None:
        return False, None
    if isinstance(parent, (str, Path)):
        path = Path(parent)
        if not path.is_file():
            return False, None
        try:
            parent = _load_json(path)
        except (OSError, ValueError, TypeError):
            return False, None
    if not isinstance(parent, Mapping):
        return False, None
    digest = _mapping_digest(parent)
    return digest == _mapping_digest(cutoff), digest


def _raw_tail_rows(sample: Mapping[str, Any], orders: Sequence[int]) -> tuple[dict[str, Any], dict[str, Any]]:
    tail = sample.get("analytic_preconditioned_tail")
    derivative_tail = sample.get("derivative_tail")
    if not isinstance(tail, Mapping) or not isinstance(derivative_tail, Mapping):
        raise ValueError("actual_first_keep lacks saved tail records")
    required = ("initial_base_common_tube_bound", "geometric_tail_ratio_upper", "omitted_sum_bound")
    if not all(_record_has_tuples(tail.get(name)) for name in required):
        raise ValueError("saved Mbase, beta, or omitted tail record lacks exact endpoints")
    gap = derivative_tail.get("Cauchy_gap")
    if not _record_has_tuples(gap):
        raise ValueError("saved Cauchy gap lacks exact endpoints")
    mbase = tail["initial_base_common_tube_bound"]
    beta = tail["geometric_tail_ratio_upper"]
    saved_tail = tail["omitted_sum_bound"]
    saved_z1 = derivative_tail.get("remaining_first_Z_derivative_norm_upper")
    if not _record_has_tuples(saved_z1):
        raise ValueError("saved first ordinary-Z tail lacks exact endpoints")
    import mpmath as mp

    with mp.workprec(MP_WORK_PRECISION):
        mbase_u = _upper_value(mbase)
        beta_u = _upper_value(beta)
        gap_u = _upper_value(gap)
        gap_l = _mp_value(gap, "lower")
        if mbase_u < 0 or beta_u < 0 or beta_u >= 1 or gap_l <= 0 or gap_u <= 0:
            raise ValueError("saved raw tail inputs fail positivity/geometric guards")
        old_iv_precision = mp.iv.prec
        mp.iv.prec = MP_WORK_PRECISION
        try:
            mbase_iv = _iv_from_record(mbase, "upper")
            beta_iv = _iv_from_record(beta, "upper")
            gap_iv = mp.iv.mpf([_mp_tuple_value(_tuple(gap, "lower") or ()), _mp_tuple_value(_tuple(gap, "upper") or ())])
            recomputed_iv = mbase_iv * beta_iv ** (WORD_DEPTH + 1) / (1 - beta_iv)
            recomputed_z1_iv = recomputed_iv / gap_iv
            recomputed_tuple = _iv_endpoint(recomputed_iv)
            recomputed_z1_tuple = _iv_endpoint(recomputed_z1_iv)
        finally:
            mp.iv.prec = old_iv_precision
        recomputed = _mp_tuple_value(recomputed_tuple)
        recomputed_z1 = _mp_tuple_value(recomputed_z1_tuple)
        source_tail_u = _upper_value(saved_tail)
        source_z1_u = _upper_value(saved_z1)
        if recomputed > source_tail_u or recomputed_z1 > source_z1_u:
            raise ValueError("fresh directed N = 1 arithmetic exceeds its saved source endpoint")
        arithmetic = {
            "positive_Mbase_and_beta_below_one": bool(mbase_u >= 0 and 0 <= beta_u < 1),
            "Cauchy_gap_positive": bool(_mp_value(gap, "lower") > 0),
            "natural_word_depth": WORD_DEPTH,
            "retained_preconditioned_indices": list(tail.get("retained_preconditioned_indices", [])),
            "formula": "Mbase_upper * beta_N_upper^(N+1) / (1-beta_N_upper)",
            "N": WORD_DEPTH,
            "Mbase_upper": _interval_projection(mbase),
            "beta_N_upper": _interval_projection(beta),
            "recomputed_upper_exact_mpf_tuple": recomputed_tuple,
            "recomputed_upper_mpf_decimal": mpf_string(recomputed),
            "recomputed_upper_log10_abs": mpf_log10_string(recomputed),
            "saved_omitted_sum_upper": _interval_projection(saved_tail),
            "recomputed_le_saved_upper": True,
            "ordinary_Z_gap_source": "coupled.samples.actual_first_keep.derivative_tail.Cauchy_gap",
            "ordinary_Z_factor": "p! / Cauchy_gap^p",
            "fresh_recomputed_order_1_upper_exact_mpf_tuple": recomputed_z1_tuple,
            "fresh_recomputed_order_1_upper_mpf_decimal": mpf_string(recomputed_z1),
            "fresh_recomputed_order_1_log10_abs": mpf_log10_string(recomputed_z1),
            "saved_order_1_upper": _interval_projection(saved_z1),
            "recomputed_order_1_le_saved_upper": True,
            "binary_exponents_are_factored": True,
            "interval_arithmetic_precision_bits": MP_WORK_PRECISION,
            "no_giant_fraction_denominator_materialized": True,
        }
        rows = {}
        for order in orders:
            if order == 0:
                endpoint = saved_tail
                fresh = recomputed
                fresh_tuple = recomputed_tuple
                source_name = "analytic_preconditioned_tail.omitted_sum_bound"
            elif order == 1:
                endpoint = saved_z1
                fresh = recomputed_z1
                fresh_tuple = recomputed_z1_tuple
                source_name = "derivative_tail.remaining_first_Z_derivative_norm_upper"
            else:
                # The coupled source only saves the first Cauchy derivative.  A
                # higher order row is not silently synthesized from a scalar.
                rows[str(order)] = {
                    "ordinary_Z_order": order,
                    "directed_upper_available": False,
                    "reason": "saved coupled first-step tail exposes only the Cauchy first-Z gap",
                }
                continue
            rows[str(order)] = {
                "ordinary_Z_order": order,
                "directed_upper_available": True,
                "tail_upper": _interval_projection(endpoint),
                "tail_upper_source": source_name,
                "fresh_recomputed_upper_exact_mpf_tuple": fresh_tuple,
                "fresh_recomputed_upper_mpf_decimal": mpf_string(fresh),
                "fresh_recomputed_upper_log10_abs": mpf_log10_string(fresh),
                "fresh_recomputed_le_saved_upper": True,
                "cauchy_gap": _interval_projection(gap),
                "cauchy_factor_formula": "factorial(p) / Cauchy_gap^p",
            }
    return arithmetic, rows


def mpf_string(value: Any) -> str:
    import mpmath as mp

    return mp.nstr(value, 80)


def mpf_log10_string(value: Any) -> str:
    import mpmath as mp

    if value == 0:
        return "-inf"
    return mp.nstr(mp.log10(abs(value)), 80)


def hydrate_selected_word_tail(
    root: str | Path | None = None,
    parent: Any = None,
    *,
    ordinary_Z_orders: Sequence[int] = ORDINARY_Z_ORDERS,
) -> dict[str, Any]:
    """Read the selected first-collar source records and hydrate N = 1 bounds."""
    repo = _repo_root(root)
    orders = tuple(int(order) for order in ordinary_Z_orders)
    if not orders or any(order < 0 or order > 3 for order in orders):
        raise ValueError("ordinary_Z_orders must be a nonempty subset of 0..3")
    cutoff_path = repo / CUTOFF_REPORT
    coupled_path = repo / COUPLED_ARCHIVE
    cutoff = _load_json(cutoff_path)
    coupled = _load_gzip_json(coupled_path)
    first = coupled.get("samples", {}).get("actual_first_keep")
    if not isinstance(first, Mapping):
        raise ValueError("saved coupled archive has no actual_first_keep sample")
    raw_remainder = cutoff.get("actual_first_raw_word_remainder")
    norm = cutoff.get("n1_actual_canonical_norm_tail")
    if not isinstance(raw_remainder, Mapping) or not isinstance(norm, Mapping):
        raise ValueError("latest cutoff binding lacks raw remainder or normalized tail records")
    parent_ok, parent_digest = _source_parent_gate(cutoff if parent is None else parent, cutoff)
    arithmetic, tail_rows = _raw_tail_rows(first, orders)
    coupled_text = (repo / COUPLED_SOURCE).read_text(encoding="utf-8") if (repo / COUPLED_SOURCE).is_file() else ""
    source_pins = {
        "cutoff_report": _pin(repo, CUTOFF_REPORT, CUTOFF_REPORT_SHA256),
        "coupled_archive": _pin(repo, COUPLED_ARCHIVE, COUPLED_ARCHIVE_SHA256),
        "coupled_source": _pin(repo, COUPLED_SOURCE, COUPLED_SOURCE_SHA256),
        "coupled_check": _pin(repo, COUPLED_CHECK, COUPLED_CHECK_SHA256),
        "inner_volterra_archive": _pin(repo, INNER_ARCHIVE, INNER_ARCHIVE_SHA256),
        "native_collar_parent": _pin(repo, PARENT_REPORT, PARENT_REPORT_SHA256),
        "remainder_module": _pin(repo, REMAINDER_MODULE, cutoff["implementation_hashes"][REMAINDER_MODULE]),
    }
    pin_ok = all(value.get("exists") and value.get("sha256_match") for value in source_pins.values())
    source_family_ok = (
        cutoff.get("source_family") == SOURCE_FAMILY
        and coupled.get("source_family") == SOURCE_FAMILY
        and first.get("source_family") == SOURCE_FAMILY
        and raw_remainder.get("source_family") == SOURCE_FAMILY
    )
    tail = first.get("analytic_preconditioned_tail", {})
    derivative_tail = first.get("derivative_tail", {})
    uniform_scope_ok = (tail.get("all_k_derivative_words_allowed") is True and
        tail.get("scope") == "uniform analytic series tail for canonical core/actual-first-collar problem; later completed-leading supports unbound")
    tail_formula_source_ok = _contains_all(
        coupled_text,
        (
            "def preconditioned_tail",
            "tail=preconditioned_tail(c,view['initial_Gg_common_tube_bound']",
            "Mbase*[a*(1+a)*(Cp+C1*k/loss)]^k/k!",
            "all_k_derivative_words_allowed=True",
            "old_unpreconditioned_tail_reused=False",
            "derivative_margin=tail['solution_Z_radius']/2",
        ),
    )
    raw_tail_gate = (
        raw_remainder.get("passed") is True
        and raw_remainder.get("actual_all_word_cauchy_majorant_available") is True
        and raw_remainder.get("raw_preconditioned_tail_majorant_source_bound") is True
        and raw_remainder.get("normalized_tail_not_used_as_raw_tail") is True
        and raw_remainder.get("numeric_actual_first_chart_field_evaluated") is False
    )
    normalized_map = {
        "supported": False,
        "applied": False,
        "raw_tail_only": True,
        "axis_transform_source": norm.get("axis_transform"),
        "derivative_map_source": norm.get("derivative_map"),
        "reason": "S^{-1} has Z-dependent scale factors; scalar division would omit the product-rule S_Z terms",
        "source_rows_are_reference_majorants_only": norm.get("certified_values_are_upper_majorants_not_function_values") is True,
    }
    checks = {
        "source_pins": pin_ok,
        "source_family": source_family_ok,
        "parent_argument_gate": parent_ok,
        "cutoff_parent_source_file_pin": source_pins["native_collar_parent"].get("sha256_match") is True,
        "cutoff_report_passed": cutoff.get("actual_n1_cutoff_history_identity", {}).get("passed") is True,
        "coupled_source_family": coupled.get("source_family") == SOURCE_FAMILY,
        "uniform_core_first_tail_scope": uniform_scope_ok,
        "actual_N1_indices": tail.get("retained_preconditioned_indices") == [0, 1],
        "raw_tail_source_gate": raw_tail_gate,
        "tail_formula_source": tail_formula_source_ok,
        "Mbase_beta_gap_exact_tuples": all(
            _record_has_tuples(tail.get(name)) for name in ("initial_base_common_tube_bound", "geometric_tail_ratio_upper", "omitted_sum_bound")
        ) and _record_has_tuples(derivative_tail.get("Cauchy_gap")) and _record_has_tuples(derivative_tail.get("remaining_first_Z_derivative_norm_upper")),
        "beta_N_strictly_below_one": arithmetic["positive_Mbase_and_beta_below_one"],
        "Cauchy_gap_positive": arithmetic["Cauchy_gap_positive"],
        "fresh_order_0_arithmetic": tail_rows.get("0", {}).get("fresh_recomputed_le_saved_upper") is True,
        "fresh_order_1_arithmetic": tail_rows.get("1", {}).get("fresh_recomputed_le_saved_upper") is True,
        "ordinary_Z_orders_requested": list(orders) == [0, 1] or all(order in (0, 1) for order in orders),
        "normalized_scalar_division_not_used": normalized_map["applied"] is False,
        "source_execution_guard": (
            cutoff.get("actual_n1_cutoff_history_identity", {}).get("source_executed") is False
            and raw_remainder.get("source_executed") is False
        ),
        "accepted_ancestor_replay_guard": (
            raw_remainder.get("old_numeric_producer_executed") is False
            and cutoff.get("actual_cutoff_m15_tail_transport", {}).get("accepted_ancestor_or_quadrature_replayed") is False
            and cutoff.get("actual_first_raw_n1_function", {}).get("accepted_owner_builder_fixture_recurrence_quadrature_replayed") is False
        ),
        "numeric_field_value_guard": (
            raw_remainder.get("numeric_actual_first_chart_field_evaluated") is False
            and cutoff.get("actual_first_raw_n1_function", {}).get("point_values_selected") is False
            and cutoff.get("actual_first_raw_n1_function", {}).get("finite_N_numerical_calculation") is False
            and norm.get("computed_word_partial_sum_or_source_point_values") is False
        ),
    }
    passed = all(checks.values())
    return {
        "module": "original_n1_selected_word_tail",
        "passed": passed,
        "source_only": True,
        "source_executed": False,
        "accepted_ancestors_replayed": False,
        "actual_N1_tail_evaluated": passed,
        "actual_N1_tail_bound_hydrated": passed,
        "tail_is_a_directed_majorant_not_a_field_value": True,
        "numeric_actual_first_chart_field_evaluated": False,
        "signed_point_value_evaluated": False,
        "natural_word_depth": WORD_DEPTH,
        "ordinary_Z_orders": list(orders),
        "scope": "uniform N=1 raw core/first-collar tail majorant",
        "tail_source_sample": "actual_first_keep tail constants only; saved field point excluded",
        "parent_argument_digest": parent_digest,
        "source_family": dict(SOURCE_FAMILY),
        "source_pins": source_pins,
        "checks": checks,
        "arithmetic": arithmetic,
        "raw_tail_enclosures": tail_rows,
        "normalized_error_mapping": normalized_map,
        "limitations": [
            "The arithmetic hydrates the saved raw all-word majorant at N=1; it does not evaluate a signed first-chart field.",
            "No finite-word field or signed point value is formed by this tail-only module.",
            "The normalized S^{-1} map is not applied because its ordinary-Z derivative action is not scalar division.",
            "Only ordinary Z orders 0 and 1 have saved coupled Cauchy-tail data; higher requested orders remain unavailable.",
        ],
    }


def build(root: str | Path | None = None, parent: Any = None) -> dict[str, Any]:
    """Build the source-pinned selected-word tail receipt."""
    return hydrate_selected_word_tail(root, parent)


__all__ = [
    "SOURCE_FAMILY",
    "hydrate_selected_word_tail",
    "build",
]
