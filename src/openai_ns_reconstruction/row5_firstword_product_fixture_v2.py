"""Fresh synthetic checks for :mod:`row5_firstword_product_v2`.

This fixture is independent of v1.  It supplies fresh IDs, paths and rational
values, and checks directed containment by constructing directed point
intervals from exact MPF endpoint tuples.  It never compares raw MPF tuples,
which can have misleading lexicographic ordering for negative values.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Mapping


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
STATE_COMPONENTS = (1, 2, 3, 5)
OUTPUT_ROWS = ("h5", "h5_Z1")

__all__ = ["SOURCE_POWER_BASIS", "STATE_COMPONENTS", "OUTPUT_ROWS", "run_fixture", "run"]


def _clone(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clone(item) for item in value]
    return value


def _q(c: Any, value: Fraction | int) -> Any:
    rational = Fraction(value)
    return c.mpf(rational.numerator) / rational.denominator


def _raw_endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None and len(raw) == 2:
        return raw[1 if upper else 0]
    raw = getattr(value, "_mpf_", None)
    if raw is not None:
        return raw
    if hasattr(value, "a") and hasattr(value, "b"):
        endpoint = value.b if upper else value.a
        nested = getattr(endpoint, "_mpi_", None)
        if nested is not None and len(nested) == 2:
            return nested[1 if upper else 0]
        nested = getattr(endpoint, "_mpf_", None)
        if nested is not None:
            return nested
        return endpoint
    return value


def _raw_endpoints(value: Any) -> tuple[Any, Any]:
    return _raw_endpoint(value, False), _raw_endpoint(value, True)


def _interval(c: Any, lower: Fraction | int, upper: Fraction | int) -> Any:
    left = _q(c, lower)
    right = _q(c, upper)
    return c.make_mpf((_raw_endpoint(left, False), _raw_endpoint(right, True)))


def _point_interval(c: Any, raw: Any) -> Any:
    endpoint = tuple(raw) if isinstance(raw, (list, tuple)) else raw
    return c.make_mpf((endpoint, endpoint))


def _inside(c: Any, enclosure: Any, value: Fraction | int) -> bool:
    """Check enclosure using interval endpoints, never raw-tuple ordering."""

    if isinstance(enclosure, Mapping) and "lower_exact_mpf_tuple" in enclosure:
        lower_raw = tuple(enclosure["lower_exact_mpf_tuple"])
        upper_raw = tuple(enclosure["upper_exact_mpf_tuple"])
    else:
        lower_raw, upper_raw = _raw_endpoints(enclosure)
    lower = _point_interval(c, lower_raw)
    upper = _point_interval(c, upper_raw)
    point = _q(c, value)
    return bool(lower <= point.a and point.b <= upper)


def _term(
    c: Any,
    coefficient: Fraction | int,
    source_id: str,
    path: Any,
    powers: tuple[int, int, int, int],
    *,
    bounds: tuple[Fraction | int, Fraction | int] | None = None,
    scalar: int = 1,
    source_family: str | None = None,
    role: str | None = None,
    reference: tuple[int, int] | None = None,
    uncertainty: bool | None = False,
    error_only: bool | None = False,
    symmetric_error: bool | None = False,
    signed_field: bool | None = True,
    duplicate_id: str | None = None,
) -> dict[str, Any]:
    directed = _q(c, coefficient) if bounds is None else _interval(c, bounds[0], bounds[1])
    result: dict[str, Any] = {
        "coefficient": directed,
        "coefficient_interval": directed,
        "formal_powers": list(powers),
        "powers": list(powers),
        "source_ids": [duplicate_id or source_id],
        "source_path": _clone(path),
        "source_paths": [_clone(path)],
        "error_classes": [f"v2-error:{source_id}"],
        "source_error_classes": [f"v2-error:{source_id}"],
        "uncertainty": uncertainty,
        "error_only": error_only,
        "symmetric_error": symmetric_error,
        "signed_field": signed_field,
        "scalar": scalar,
    }
    if source_family is not None:
        result["source_family"] = source_family
    if role is not None:
        result["role"] = role
    if reference is None:
        result["source_only"] = source_family == "g5"
    else:
        result["state_reference"] = {"component": reference[0], "ordinary_Z_order": reference[1]}
    return result


def _state_term(
    c: Any,
    component: int,
    order: int,
    coefficient: Fraction | int,
    powers: tuple[int, int, int, int],
    *,
    bounds: tuple[Fraction | int, Fraction | int] | None = None,
    uncertainty: bool | None = False,
    error_only: bool | None = False,
    symmetric_error: bool | None = False,
    signed_field: bool | None = True,
) -> dict[str, Any]:
    return _term(
        c,
        coefficient,
        f"v2:state:W{component}:Z{order}",
        ["v2", "state", f"W{component}", f"Z{order}", 0],
        powers,
        bounds=bounds,
        uncertainty=uncertainty,
        error_only=error_only,
        symmetric_error=symmetric_error,
        signed_field=signed_field,
    )


def _expression(c: Any) -> dict[str, Any]:
    duplicate = "v2:shared:g5:occurrence"
    return {
        "schema": "row5_affine_forcing_v1",
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "rows": {
            "h5": [
                _term(c, 1, "v2:g5:error", ["v2", "receipt", "g5", "Z0", 0], (0, 0, 0, 0), bounds=(-2, 2), source_family="g5", role="source_only_g5", uncertainty=True, error_only=True, symmetric_error=True, signed_field=False, duplicate_id=duplicate),
                _term(c, -3, "v2:g5:main", ["v2", "receipt", "g5", "Z0", 1], (1, 0, 0, 0), source_family="g5", role="source_only_g5", duplicate_id=duplicate),
                _term(c, Fraction(4, 9), "v2:C0:W1", ["v2", "receipt", "C0", "Z0", "5,1"], (1, -2, 0, 1), source_family="C0", role="C0_Z0_times_state", reference=(1, 0)),
                _term(c, -5, "v2:C0:W1:dup", ["v2", "receipt", "C0", "Z0", "5,1", "dup"], (0, 0, 1, 0), source_family="C0", role="C0_Z0_times_state_duplicate_source", reference=(1, 0)),
                _term(c, -6, "v2:B1:W2Z1", ["v2", "receipt", "B1", "Z0", "5,2"], (0, 1, 0, -1), source_family="B1", role="B1_Z0_times_state_Z1", reference=(2, 1)),
                _term(c, Fraction(13, 17), "v2:axis:chi", ["v2", "receipt", "axis", "chi"], (0, 0, 0, 0), source_family="axis", role="two_chi_times_W1", reference=(1, 0), scalar=2, uncertainty=True),
            ],
            "h5_Z1": [
                _term(c, 7, "v2:g5:Z1", ["v2", "receipt", "g5", "Z1"], (0, 1, 0, 0), source_family="g5", role="source_only_g5_Z1"),
                _term(c, Fraction(-7, 11), "v2:C0:Z1:W3", ["v2", "receipt", "C0", "Z1", "5,3"], (0, 0, 0, -1), source_family="C0", role="C0_Z1_times_state", reference=(3, 0)),
                _term(c, 5, "v2:C0:Z0:W5Z1", ["v2", "receipt", "C0", "Z0", "5,5"], (0, 0, 1, 0), source_family="C0", role="C0_Z0_times_state_Z1", reference=(5, 1)),
                _term(c, -9, "v2:B1:Z1:W1Z1", ["v2", "receipt", "B1", "Z1", "5,1"], (1, 0, 0, 0), source_family="B1", role="B1_Z1_times_state_Z1", reference=(1, 1)),
                _term(c, 8, "v2:B1:Z0:W3Z2", ["v2", "receipt", "B1", "Z0", "5,3"], (0, -1, 0, 0), source_family="B1", role="B1_Z0_times_state_Z2", reference=(3, 2)),
                _term(c, -4, "v2:B1:Z0:W1Z2", ["v2", "receipt", "B1", "Z0", "5,1"], (0, 0, 0, 0), source_family="B1", role="B1_Z0_times_state_Z2", reference=(1, 2)),
                _term(c, 7, "v2:B1:Z0:W2Z2", ["v2", "receipt", "B1", "Z0", "5,2"], (0, 0, 0, 0), source_family="B1", role="B1_Z0_times_state_Z2", reference=(2, 2)),
                _term(c, -15, "v2:axis:chiZ", ["v2", "receipt", "axis", "chi_Z"], (0, 0, 0, 0), source_family="axis", role="two_chi_Z_times_W1", reference=(1, 0), scalar=2, uncertainty=True),
                _term(c, 9, "v2:axis:chi", ["v2", "receipt", "axis", "chi"], (0, 0, 0, 0), source_family="axis", role="two_chi_times_W1_Z", reference=(1, 1), scalar=2, uncertainty=True),
            ],
        },
    }


def _state(c: Any) -> dict[str, dict[str, list[dict[str, Any]]]]:
    rows: dict[str, dict[str, list[dict[str, Any]]]] = {str(component): {} for component in STATE_COMPONENTS}
    rows["1"]["Z0"] = [_state_term(c, 1, 0, 7, (0, 3, -1, -1), uncertainty=True)]
    rows["1"]["Z1"] = [_state_term(c, 1, 1, -3, (-1, 0, 0, 1))]
    rows["1"]["Z2"] = [_state_term(c, 1, 2, 5, (0, 0, 1, 0))]
    rows["2"]["Z0"] = [_state_term(c, 2, 0, 11, (0, 0, 0, 0))]
    rows["2"]["Z1"] = [_state_term(c, 2, 1, 0, (0, 0, 0, 2), bounds=(-2, 2), uncertainty=True, error_only=True, symmetric_error=True, signed_field=False)]
    rows["2"]["Z2"] = [_state_term(c, 2, 2, -11, (1, 0, 0, 0))]
    rows["3"]["Z0"] = [_state_term(c, 3, 0, 19, (0, -1, 0, 0))]
    rows["3"]["Z1"] = [_state_term(c, 3, 1, -23, (0, 0, 0, 0))]
    rows["3"]["Z2"] = [_state_term(c, 3, 2, -17, (0, 0, -1, 0))]
    rows["5"]["Z0"] = [_state_term(c, 5, 0, 29, (0, 0, 0, 0))]
    rows["5"]["Z1"] = [_state_term(c, 5, 1, 31, (0, 1, 0, 0))]
    return rows


def _callable(value: Any, name: str) -> Any:
    candidate = getattr(value, name, value)
    if not callable(candidate):
        raise TypeError(f"{name} must be callable or expose {name}")
    return candidate


def _rejected(call: Any) -> bool:
    try:
        call()
    except (TypeError, ValueError, KeyError, ZeroDivisionError):
        return True
    return False


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    build_fn = _callable(build, "build")
    expression = _expression(c)
    state = _state(c)
    try:
        result = build_fn(c, expression, state)
    except Exception as exc:
        return {"schema": "row5_firstword_product_fixture_v2", "passed": False, "checks": {"build_completed": False}, "check_count": 1, "error": f"{type(exc).__name__}: {exc}", "synthetic_only": True, "implementation_imported": False}

    checks: dict[str, bool] = {}
    checks["build_completed"] = isinstance(result, Mapping)
    checks["schema"] = result.get("schema") == "row5_firstword_product_v2"
    checks["source_schema_preserved"] = result.get("source_schema") == "row5_affine_forcing_v1"
    checks["compact_top_level_shape"] = all(key in result for key in ("rows", "source_nodes", "state_nodes")) and all(
        key not in result for key in ("products", "h5_products", "h5_Z1_products", "source_node_table", "state_node_table")
    )
    checks["products_are_compact"] = all(
        "source_node_id" in product
        and "state_node_id" in product
        and "coefficient_enclosure" in product
        and product.get("scalar_applied_to_coefficient") is False
        and product.get("term_semantics") == "scalar * coefficient_enclosure * formal_monomial; scalar deferred"
        and all(key not in product for key in ("coefficient", "coefficient_interval", "coefficient_endpoints", "source_ids", "source_paths"))
        for label in OUTPUT_ROWS
        for product in result["rows"][label]
    )
    checks["directed_enclosures"] = all(
        isinstance(product["coefficient_enclosure"], Mapping)
        and product["coefficient_enclosure"].get("directed") is True
        and "lower_exact_mpf_tuple" in product["coefficient_enclosure"]
        and "upper_exact_mpf_tuple" in product["coefficient_enclosure"]
        for label in OUTPUT_ROWS
        for product in result["rows"][label]
    )
    checks["not_function_value"] = result.get("row_totals") is None and result.get("interval_sum_performed") is False and result.get("full_h5") is False and result.get("formal_atoms_materialized") is False
    checks["residual_separate"] = result.get("unknown_residual", {}).get("symbol") == "R" and result["unknown_residual"].get("included_in_products") is False

    source_nodes = result["source_nodes"]
    state_nodes = result["state_nodes"]
    products = result["rows"]
    all_products = products["h5"] + products["h5_Z1"]
    checks["paths_present"] = all(node.get("parent_input_paths") for node in list(source_nodes.values()) + list(state_nodes.values()))
    checks["nonzero_source_present"] = any(not _inside(c, node["coefficient_interval"], 0) for node in source_nodes.values())
    checks["source_only_g5_occurrences_retained"] = sum(p["source_only"] for p in products["h5"]) == 2 and sum(p["source_only"] for p in products["h5_Z1"]) == 1
    checks["source_occurrences_not_deduped"] = len([node for node in source_nodes.values() if node.get("provenance_ids") == ["v2:shared:g5:occurrence"]]) == 2
    checks["shared_state_node_reused"] = sum(product.get("state_node_id") == "state:W1:Z0:0" for product in all_products) >= 3
    checks["error_sign_provenance_kept"] = any(source_nodes[p["source_node_id"]].get("error_only") is True for p in products["h5"]) and any(p.get("state_node_id") is not None and state_nodes[p["state_node_id"]].get("error_only") is True for p in products["h5"])
    checks["no_cancellation"] = all(p["provenance"].get("cancellation_performed") is False for p in all_products)

    def source_node(product: Mapping[str, Any]) -> Mapping[str, Any]:
        return source_nodes[product["source_node_id"]]

    axis_products = [p for p in all_products if source_node(p).get("source_family") == "axis"]
    checks["axis_scalar_separate"] = bool(axis_products) and all(p.get("scalar") == 2 and p.get("scalar_applied_to_coefficient") is False for p in axis_products)
    checks["axis_coefficient_once"] = any(_inside(c, p["coefficient_enclosure"], Fraction(91, 170)) for p in axis_products)
    checks["positive_product_contained"] = any(_inside(c, p["coefficient_enclosure"], Fraction(14, 45)) for p in products["h5"] if source_node(p).get("role") == "C0_Z0_times_state")
    checks["negative_product_contained"] = any(_inside(c, p["coefficient_enclosure"], Fraction(-133, 55)) for p in products["h5_Z1"] if source_node(p).get("role") == "C0_Z1_times_state")
    checks["formal_power_sum"] = next(p for p in products["h5"] if source_node(p).get("role") == "C0_Z0_times_state").get("formal_powers") == [1, 1, -1, 0]
    checks["B1_W_ZZ"] = next(p for p in products["h5_Z1"] if source_node(p).get("role") == "B1_Z0_times_state_Z2").get("state_node_id") == "state:W3:Z2:0"
    checks["ordinary_Z1_roles"] = {source_node(p).get("role") for p in products["h5_Z1"]} >= {"C0_Z1_times_state", "C0_Z0_times_state_Z1", "B1_Z1_times_state_Z1", "B1_Z0_times_state_Z2", "two_chi_Z_times_W1", "two_chi_times_W1_Z"}
    checks["negative_powers"] = any(any(power < 0 for power in p.get("formal_powers", ())) for p in all_products)
    conditions = result.get("conditions", {})
    checks["owner_guard"] = all(conditions.get(key) == "caller condition; not proved here" for key in ("same_frame", "same_scale", "same_gauge", "live_owner", "same_state_owner", "same_state_Z0_Z1_Z2_owner"))
    checks["scope_flags"] = result.get("scope", {}).get("source_only_g5_retention_scope") == "one output per input occurrence; no global deduplication" and result.get("scope", {}).get("no_product_deduplication") is True and result.get("scope", {}).get("no_error_cancellation") is True and result.get("scalar_application") == "deferred_to_consumer_once"

    missing_checks: dict[str, bool] = {}
    for component in (1, 2, 3):
        missing = _clone(state)
        del missing[str(component)]["Z2"]
        missing_checks[str(component)] = _rejected(lambda s=missing: build_fn(c, expression, s))
    checks["missing_W1_W2_W3_Z2_rejected"] = all(missing_checks.values())
    no_w5_z2 = _clone(state)
    checks["W5_Z2_not_required"] = _rejected(lambda: build_fn(c, expression, no_w5_z2)) is False
    with_atoms = _clone(state)
    with_atoms["formal_atoms"] = {"hb": _q(c, 2)}
    checks["atom_values_rejected"] = _rejected(lambda: build_fn(c, expression, with_atoms))
    bad_scalar = _clone(state)
    bad_scalar["1"]["Z0"][0]["scalar"] = 2
    checks["nonunit_state_scalar_rejected"] = _rejected(lambda: build_fn(c, expression, bad_scalar))
    bad_basis = _clone(state)
    bad_basis["1"]["Z0"][0]["power_basis"] = ["wrong", "basis", "only", "here"]
    checks["bad_basis_rejected"] = _rejected(lambda: build_fn(c, expression, bad_basis))
    bad_float_expression = _clone(expression)
    bad_float_expression["rows"]["h5"][2]["coefficient"] = 1.25
    bad_float_expression["rows"]["h5"][2]["coefficient_interval"] = 1.25
    checks["python_float_rejected"] = _rejected(lambda: build_fn(c, bad_float_expression, state))
    bad_float_state = _clone(state)
    bad_float_state["1"]["Z0"][0]["coefficient"] = [1.0, 2.0]
    bad_float_state["1"]["Z0"][0]["coefficient_interval"] = [1.0, 2.0]
    checks["python_float_endpoint_rejected"] = _rejected(lambda: build_fn(c, expression, bad_float_state))

    return {
        "schema": "row5_firstword_product_fixture_v2",
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "ordinary_Z1_product_rule_checked": True,
        "B1_W_ZZ_checked": True,
        "owner_guard_checked": True,
        "atom_guard_checked": True,
        "float_guard_checked": True,
        "implementation_imported": False,
    }


def run(build: Any, c: Any) -> dict[str, Any]:
    return run_fixture(build, c)
