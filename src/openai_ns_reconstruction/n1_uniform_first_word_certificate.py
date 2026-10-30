"""Conditional uniform first-word enclosure for the saved first collar.

This is a pure algebra leaf.  It reads no archive, source file, graph, or
endpoint tuple and it executes no producer.  The caller supplies the already
saved ``q=1`` rows, the already saved raw ``N=1`` all-word majorant, an
interval lifting function, and an interval hull.  The q=1 transport row is
useful here because :func:`original_n1_selected_word_numeric.transport` puts
``y`` in the *whole* box ``[0, sc/8]`` at q=1.  Thus that row contains every
phase ``q in [0, 1]`` once the caller supplies the positive-
``sc``/``0 < hb < 1/2`` and phase-monotonicity premises.

The output is deliberately conditional.  It keeps the signed source-sector
prefix formal (coefficient interval times the same formal positive source
exponential) and adds, separately, the two positive local radii and the raw
all-word tail.  It never turns an analysis check, norm cap, or tail majorant
into a field value.  It also never identifies a core extension with a collar.

The saved selected-word packet supports ordinary Z orders 0 and 1.  The two
requested fields are kept independent.  ``eta`` is only a phase/cutoff
premise: the returned certificate records ``0 <= eta <= 1``, the original
one/transition/zero support through ``3/4``, and ``eta_Z = 0``; it does not
differentiate or evaluate eta.

Expected caller data
--------------------

``saved_q1_rows`` is either one selected-word row, a sequence of rows, or a
mapping containing ``selected_coordinate_rows``/``q1_rows``.  Rows must have
``q == "1"`` and an ``ordinary_Z`` mapping.  For every requested order and
field, the order row must contain ``prefix``, ``W0``, ``W1`` and
``local_raw_error_radii``.  Prefix terms have the saved shape
``source_exponents``, ``coefficient_interval`` and
``combined_source_log_interval``.

``saved_tail`` is the returned mapping of
``original_n1_selected_word_tail.hydrate_selected_word_tail`` (or an
equivalent source-bound mapping).  Its ``raw_tail_enclosures`` rows expose
``tail_upper`` for orders 0 and 1.  The source-bound tail flags are checked so
that a normalized tail or a selected point value cannot enter this leaf.

``interval(value)`` lifts one saved record to the caller's directed interval
object.  ``hull(values)`` returns the caller's outward interval hull.  The
caller also supplies ``nonnegative(value)``; it must return exactly ``True``
for every local/tail radius and their sums.  No endpoint attribute is
inspected here, so the interval implementation remains caller-owned.

``source_binding`` and ``phase_monotonicity`` are explicit evidence maps.
Their required flags are intentionally not inferred from names or from a
numeric overlap.  A successful result therefore means
``conditional_uniform_certificate`` rather than an unconditional source
proof.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from typing import Any


FIELDS = ("F1", "Uz1")
ORDERS = (0, 1)
ALL_SELECTED_FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")

# These are the support facts that the source packet records for the original
# eta(q) = 1 - sigma(4q - 2).  They are returned as premises, never evaluated.
ETA_SUPPORT = {
    "one": "q in [0, 1/2]",
    "transition": "q in (1/2, 3/4)",
    "zero": "q in [3/4, 1]",
}

_SOURCE_FLAGS = (
    "same_source_family",
    "same_raw_first_function",
    "same_raw_core_entry",
    "same_raw_all_word_tail",
    "raw_vector_tail_dominates_each_field",
    "ordinary_Z_rows_same_source",
    "same_saved_Z_tube",
    "prefix_is_homogeneous_not_function_value",
    "tail_is_majorant_not_function_value",
    "analysis_checks_not_function_values",
    "core_extension_not_used_as_collar",
)

_PHASE_FLAGS = (
    "q_interval_is_exact_0_to_1",
    "eta_in_unit_interval",
    "eta_support_through_3_over_4",
    "eta_Z_zero",
    "sc_positive",
    "hb_positive_less_half",
    "transport_q1_y_box_is_0_to_sc_over_8",
    "all_q_y_is_subset_of_q1_box",
    "transport_is_Z_independent",
    "local_radii_monotone_in_qsc",
    "local_radii_nonnegative",
)


IntervalObject = Any
IntervalLift = Callable[[Any], IntervalObject]
Hull = Callable[[Iterable[IntervalObject]], IntervalObject]
Nonnegative = Callable[[IntervalObject], bool]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping")
    return value


def _require_flags(evidence: Any, required: Sequence[str], name: str) -> dict[str, Any]:
    evidence_map = _mapping(evidence, name)
    missing = [flag for flag in required if evidence_map.get(flag) is not True]
    _require(not missing, f"{name} lacks explicit true flags: {', '.join(missing)}")
    return dict(evidence_map)


def _lift(interval: IntervalLift, value: Any, label: str) -> IntervalObject:
    _require(value is not None, f"{label}: saved interval record is missing")
    try:
        result = interval(value)
    except Exception as exc:  # pragma: no cover - caller arithmetic owns errors
        raise TypeError(f"{label}: caller interval lifting failed") from exc
    _require(result is not None, f"{label}: interval lifting returned None")
    return result


def _hull(hull: Hull, values: Sequence[IntervalObject], label: str) -> IntervalObject:
    _require(values, f"{label}: empty hull input")
    try:
        result = hull(values)
    except Exception as exc:  # pragma: no cover - caller arithmetic owns errors
        raise TypeError(f"{label}: caller interval hull failed") from exc
    _require(result is not None, f"{label}: interval hull returned None")
    return result


def _check_nonnegative(nonnegative: Nonnegative, value: IntervalObject, label: str) -> None:
    try:
        passed = nonnegative(value)
    except Exception as exc:  # pragma: no cover - caller arithmetic owns errors
        raise TypeError(f"{label}: caller nonnegative check failed") from exc
    _require(passed is True, f"{label}: nonnegative radius was not certified")


def _is_one(value: Any) -> bool:
    if isinstance(value, bool):
        return False
    if isinstance(value, int):
        return value == 1
    if isinstance(value, str):
        return value.strip() == "1"
    # Avoid accepting binary floats as a phase identity.  Fraction-like
    # caller objects may still expose numerator/denominator exactly.
    numerator = getattr(value, "numerator", None)
    denominator = getattr(value, "denominator", None)
    return numerator == 1 and denominator == 1


def _row_sequence(saved_rows: Any) -> list[Mapping[str, Any]]:
    if isinstance(saved_rows, Mapping):
        for key in ("selected_coordinate_rows", "q1_rows", "rows"):
            if key in saved_rows:
                saved_rows = saved_rows[key]
                break
        else:
            if "q" in saved_rows:
                return [saved_rows]
            raise TypeError(
                "saved_q1_rows must be a q=1 row or contain selected_coordinate_rows/q1_rows"
            )
    if isinstance(saved_rows, (str, bytes)):
        raise TypeError("saved_q1_rows must not be text")
    try:
        rows = list(saved_rows)
    except TypeError as exc:
        raise TypeError("saved_q1_rows must be an iterable of rows") from exc
    _require(rows, "saved_q1_rows must not be empty")
    result: list[Mapping[str, Any]] = []
    for index, row in enumerate(rows):
        result.append(_mapping(row, f"saved_q1_rows[{index}]"))
    return result


def _q1_rows(saved_rows: Any) -> list[Mapping[str, Any]]:
    rows = [row for row in _row_sequence(saved_rows) if _is_one(row.get("q"))]
    _require(rows, "saved_q1_rows contains no q=1 row")
    return rows


def _order_row(row: Mapping[str, Any], order: int, row_index: int) -> Mapping[str, Any]:
    ordinary = _mapping(row.get("ordinary_Z"), f"q=1 row {row_index}.ordinary_Z")
    value = ordinary.get(str(order), ordinary.get(order))
    return _mapping(value, f"q=1 row {row_index}.ordinary_Z[{order}]")


def _field_rows(order_row: Mapping[str, Any], component: str, field: str, label: str) -> Any:
    component_map = _mapping(order_row.get(component), f"{label}.{component}")
    _require(field in component_map, f"{label}.{component} lacks field {field}")
    return component_map[field]


def _sector_parts(
    *,
    prefix: Mapping[str, Any],
    field: str,
    order: int,
    row_index: int,
    interval: IntervalLift,
    hull: Hull,
) -> dict[tuple[int, ...], dict[str, Any]]:
    terms = prefix.get(field)
    if not isinstance(terms, Sequence) or isinstance(terms, (str, bytes)):
        raise TypeError(f"q=1 row {row_index} Z{order} prefix[{field}] must be a term sequence")
    groups: dict[tuple[int, ...], dict[str, Any]] = {}
    for term_index, raw_term in enumerate(terms):
        term = _mapping(raw_term, f"q=1 row {row_index} Z{order} prefix[{field}][{term_index}]")
        exponents = term.get("source_exponents")
        if not isinstance(exponents, Sequence) or isinstance(exponents, (str, bytes)):
            raise TypeError(f"q=1 row {row_index} Z{order} source exponents must be a sequence")
        key = tuple(exponents)
        _require(len(key) == 4, f"q=1 row {row_index} Z{order} source basis must have four exponents")
        _require(all(isinstance(value, int) and not isinstance(value, bool) for value in key),
                 f"q=1 row {row_index} Z{order} source exponents must be exact integers")
        _require(key not in groups, f"duplicate source sector {key}")
        coefficient = _lift(
            interval,
            term.get("coefficient_interval"),
            f"q=1 row {row_index} Z{order} {field} sector {key} coefficient",
        )
        source_log = _lift(
            interval,
            term.get("combined_source_log_interval"),
            f"q=1 row {row_index} Z{order} {field} sector {key} source log",
        )
        group = groups.setdefault(key, {"coefficients": [], "source_logs": [], "row_terms": []})
        group["coefficients"].append(coefficient)
        group["source_logs"].append(source_log)
        group["row_terms"].append(term_index)
    _require(groups, f"q=1 row {row_index} Z{order} prefix[{field}] is empty")
    result: dict[tuple[int, ...], dict[str, Any]] = {}
    for key, group in groups.items():
        result[key] = {
            "coefficient_interval": _hull(
                hull, group["coefficients"],
                f"q=1 Z{order} {field} sector {key} coefficient",
            ),
            "combined_source_log_interval": _hull(
                hull, group["source_logs"],
                f"q=1 Z{order} {field} sector {key} source log",
            ),
            "source_exponents": key,
            "row_terms": tuple(group["row_terms"]),
        }
    return result


def _tail_record(saved_tail: Mapping[str, Any], order: int) -> Any:
    rows = saved_tail.get("raw_tail_enclosures")
    if isinstance(rows, Mapping):
        row = rows.get(str(order), rows.get(order))
        if isinstance(row, Mapping):
            _require(row.get("directed_upper_available") is True,
                     f"saved raw N=1 tail Z{order} is not directed")
            if "tail_upper" in row:
                return row["tail_upper"]

    # A compact caller packet may expose one record per order directly.  This
    # path remains source-bound only because the tail flags below are still
    # mandatory.
    for key in (str(order), order):
        if key in saved_tail:
            value = saved_tail[key]
            if isinstance(value, Mapping) and "tail_upper" in value:
                return value["tail_upper"]
            return value
    orders = saved_tail.get("ordinary_Z")
    if isinstance(orders, Mapping):
        value = orders.get(str(order), orders.get(order))
        if isinstance(value, Mapping):
            for key in ("tail_upper", "upper", "value"):
                if key in value:
                    return value[key]
    raise KeyError(f"saved raw N=1 tail has no order {order} upper bound")


def _tail_contract(saved_tail: Mapping[str, Any]) -> dict[str, Any]:
    checks = {
        "passed": saved_tail.get("passed") is True,
        "source_only": saved_tail.get("source_only") is True,
        "N1_bound_hydrated": saved_tail.get("actual_N1_tail_bound_hydrated") is True,
        "majorant_not_value": saved_tail.get("tail_is_a_directed_majorant_not_a_field_value") is True,
        "no_numeric_field": saved_tail.get("numeric_actual_first_chart_field_evaluated") is False,
        "no_signed_point": saved_tail.get("signed_point_value_evaluated") is False,
        "raw_not_normalized": (
            saved_tail.get("raw_tail_only") is True
            or _mapping(saved_tail.get("normalized_error_mapping", {}), "saved_tail.normalized_error_mapping").get("raw_tail_only") is True
        ),
        "normalization_not_applied": (
            not isinstance(saved_tail.get("normalized_error_mapping"), Mapping)
            or saved_tail["normalized_error_mapping"].get("applied") is False
        ),
        "natural_word_depth_N1": saved_tail.get("natural_word_depth") == 1,
    }
    _require(all(checks.values()), "saved tail is not the raw source-bound N=1 majorant")
    return checks


def assemble(
    saved_q1_rows: Any,
    saved_tail: Mapping[str, Any],
    *,
    interval: IntervalLift,
    hull: Hull,
    nonnegative: Nonnegative,
    source_binding: Mapping[str, Any],
    phase_monotonicity: Mapping[str, Any],
    fields: Sequence[str] = FIELDS,
    orders: Sequence[int] = ORDERS,
) -> dict[str, Any]:
    """Assemble a conditional uniform q=0..1 raw first-word certificate.

    ``interval`` and ``hull`` are the caller's arithmetic.  The returned
    sectors are not multiplied by a materialized exponential; they are formal
    signed source sectors.  Each error row is an additive *positive radius*:

    ``local_W0 + local_W1 + raw_all_word_N1_tail``.

    The q=1 homogeneous prefix already ranges over ``y in [0, sc/8]``.  The
    phase evidence therefore justifies reusing it for all ``q in [0,1]``;
    local radii use the q=1 monotone maximum.  No q=0 midpoint, derivative
    approximation, core extension, or collar function is synthesized here.
    """

    if not callable(interval):
        raise TypeError("interval must be the caller's interval lifting function")
    if not callable(hull):
        raise TypeError("hull must be the caller's interval hull function")
    if not callable(nonnegative):
        raise TypeError("nonnegative must be the caller's radius-certification function")
    source_evidence = _require_flags(source_binding, _SOURCE_FLAGS, "source_binding")
    phase_evidence = _require_flags(phase_monotonicity, _PHASE_FLAGS, "phase_monotonicity")
    tail = _mapping(saved_tail, "saved_tail")
    tail_checks = _tail_contract(tail)

    requested_fields = tuple(fields)
    requested_orders = tuple(orders)
    _require(requested_fields == FIELDS,
             f"only independent fields {FIELDS!r} are supported")
    _require(requested_orders and all(order in ORDERS for order in requested_orders),
             "saved q=1 raw first-word data support only ordinary Z orders 0 and 1")
    _require(len(set(requested_orders)) == len(requested_orders), "duplicate ordinary Z order")

    rows = _q1_rows(saved_q1_rows)
    _require(len(rows) == 1, "one authentic q=1 row required")
    z_tubes = [row.get("selected_Z_interval") for row in rows]
    _require(all(value is not None for value in z_tubes),
             "every q=1 row must carry the saved Z tube")
    _require(all(value == z_tubes[0] for value in z_tubes[1:]),
             "q=1 rows do not share one saved Z tube")

    by_order: dict[str, dict[str, Any]] = {}
    term_checks: dict[str, dict[str, Any]] = {field: {} for field in requested_fields}
    for order in requested_orders:
        order_rows = [_order_row(row, order, row_index) for row_index, row in enumerate(rows)]
        for row_index, order_row in enumerate(order_rows):
            _field_rows(order_row, "W0", "F1", f"q=1 row {row_index} Z{order}")
            _field_rows(order_row, "W1", "F1", f"q=1 row {row_index} Z{order}")
            _field_rows(order_row, "prefix", "F1", f"q=1 row {row_index} Z{order}")
            _field_rows(order_row, "W0", "Uz1", f"q=1 row {row_index} Z{order}")
            _field_rows(order_row, "W1", "Uz1", f"q=1 row {row_index} Z{order}")
            _field_rows(order_row, "prefix", "Uz1", f"q=1 row {row_index} Z{order}")
            radii = _mapping(order_row.get("local_raw_error_radii"),
                             f"q=1 row {row_index} Z{order}.local_raw_error_radii")
            _mapping(radii.get("W0"), f"q=1 row {row_index} Z{order}.local_raw_error_radii.W0")
            _mapping(radii.get("W1"), f"q=1 row {row_index} Z{order}.local_raw_error_radii.W1")
            _require(order_row.get("allword_tail_included") is False,
                     f"q=1 row {row_index} Z{order} prefix already includes all-word tail")

        field_results: dict[str, Any] = {}
        for field in requested_fields:
            # Prefix sectors are grouped only by exact source exponent.  The
            # q=1 row's interval already ranges over y; hulling multiple saved
            # q=1 rows keeps that formal sector structure without selecting a
            # point coefficient.
            sector_groups: dict[tuple[int, ...], list[dict[str, Any]]] = {}
            for row_index, order_row in enumerate(order_rows):
                parts = _sector_parts(
                    prefix=_mapping(order_row.get("prefix"), f"q=1 row {row_index} Z{order}.prefix"),
                    field=field,
                    order=order,
                    row_index=row_index,
                    interval=interval,
                    hull=hull,
                )
                for key, part in parts.items():
                    sector_groups.setdefault(key, []).append(part)
            sectors: list[dict[str, Any]] = []
            for key in sorted(sector_groups):
                parts = sector_groups[key]
                coefficient = _hull(hull, [part["coefficient_interval"] for part in parts],
                                    f"uniform Z{order} {field} sector {key} coefficient")
                source_log = _hull(hull, [part["combined_source_log_interval"] for part in parts],
                                   f"uniform Z{order} {field} sector {key} source log")
                sectors.append({
                    "source_exponents": key,
                    "coefficient_interval": coefficient,
                    "cutoff_coefficient_interval": _hull(hull, [coefficient, coefficient * 0],
                        f"cutoff Z{order} {field} sector {key}"),
                    "combined_source_log_interval": source_log,
                    "formal_term": "coefficient_interval * exp(combined_source_log_interval)",
                    "signed_sector": True,
                    "exponential_materialized": False,
                    "source_rows_hulled": len(parts),
                })
            _require(sectors, f"uniform Z{order} {field} signed prefix has no sectors")

            local_rows = [_mapping(order_row.get("local_raw_error_radii"),
                                  f"q=1 row {index} Z{order}.local_raw_error_radii")
                          for index, order_row in enumerate(order_rows)]
            local_parts: dict[str, IntervalObject] = {}
            local_checks: dict[str, Any] = {}
            for component in ("W0", "W1"):
                values: list[IntervalObject] = []
                component_rows = []
                for row_index, radius_map in enumerate(local_rows):
                    component_map = _mapping(radius_map.get(component),
                                             f"q=1 row {row_index} Z{order} local {component}")
                    value = _lift(
                        interval,
                        component_map.get(field),
                        f"q=1 Z{order} {field} local {component}",
                    )
                    _check_nonnegative(nonnegative, value,
                                       f"q=1 Z{order} {field} local {component}")
                    values.append(value)
                    component_rows.append(value)
                radius = _hull(hull, values, f"uniform Z{order} {field} local {component}")
                _check_nonnegative(nonnegative, radius,
                                   f"uniform Z{order} {field} local {component} hull")
                local_parts[component] = radius
                local_checks[component] = {
                    "all_q1_rows_nonnegative": True,
                    "q1_hull_nonnegative": True,
                    "row_count": len(component_rows),
                }

            local_total = local_parts["W0"] + local_parts["W1"]
            _check_nonnegative(nonnegative, local_total,
                               f"uniform Z{order} {field} local W0+W1")
            tail_value = _lift(interval, _tail_record(tail, order),
                               f"saved raw N=1 tail Z{order}")
            _check_nonnegative(nonnegative, tail_value,
                               f"saved raw N=1 tail Z{order}")
            total = local_total + tail_value
            _check_nonnegative(nonnegative, total,
                               f"uniform Z{order} {field} local+tail radius")
            field_results[field] = {
                "formal_signed_sectors": sectors,
                "homogeneous_prefix_is_q1_y_box": True,
                "homogeneous_prefix_is_not_a_field_value": True,
                "cutoff_prefix_formula": "sum cutoff_coefficient_interval * exp(same source log)",
                "cutoff_error_radius": total,
                "cutoff_zero_for_q_ge_3_over_4": True,
                "error_bound": {
                    "local_W0_radius": local_parts["W0"],
                    "local_W1_radius": local_parts["W1"],
                    "local_prefix_radius": local_total,
                    "raw_all_word_N1_tail_radius": tail_value,
                    "total_additive_raw_positive_radius": total,
                    "formula": "local_W0 + local_W1 + raw_all_word_N1_tail",
                    "all_terms_nonnegative": True,
                    "tail_added_once": True,
                },
                "term_checks": local_checks,
                "independent_field": True,
            }
            term_checks[field][str(order)] = {
                "formal_sectors_present": True,
                "local_W0_nonnegative": True,
                "local_W1_nonnegative": True,
                "local_sum_nonnegative": True,
                "raw_tail_nonnegative": True,
                "additive_local_plus_tail_nonnegative": True,
                "same_saved_Z_tube": True,
            }
        by_order[str(order)] = {
            "ordinary_Z_order": order,
            "fields": field_results,
            "tail_source": "saved_tail.raw_tail_enclosures[order].tail_upper",
        }

    return {
        "passed": True,
        "status": "conditional_uniform_certificate",
        "certificate_kind": "q1_formal_signed_sector_plus_additive_raw_positive_radius",
        "conditional_on_caller_evidence": True,
        "source_binding": source_evidence,
        "phase_monotonicity": phase_evidence,
        "q_scope": {
            "q_interval": "[0,1]",
            "eta_premise": "0 <= eta <= 1",
            "eta_support": dict(ETA_SUPPORT),
            "eta_Z": "0",
            "q1_transport_y_box": "[0, sc/8]",
            "all_q_y_subset_q1_box": True,
            "local_radii_use_q1_monotone_maximum": True,
            "same_saved_Z_tube_only": True,
            "saved_Z_tube": z_tubes[0],
        },
        "fields": list(requested_fields),
        "ordinary_Z_orders": list(requested_orders),
        "orders": by_order,
        "term_checks": term_checks,
        "tail_checks": tail_checks,
        "signed_sectors_are_not_function_values": True,
        "positive_radii_are_not_signed_increments": True,
        "core_extension_not_used_as_collar": True,
        "analysis_checks_not_promoted_to_values": True,
        "limitations": [
            "The result is conditional on the caller's source-binding and phase-monotonicity evidence.",
            "It covers only the common saved Z tube carried by the q=1 rows; it is not arbitrary-Z coverage.",
            "Formal exp(source-log) factors are not materialized and no sector midpoint is selected.",
            "The raw N=1 tail is an omitted-word majorant, not a signed field value.",
            "No core-extension identity is used as a collar identity.",
        ],
    }


__all__ = ["FIELDS", "ORDERS", "ETA_SUPPORT", "assemble"]
