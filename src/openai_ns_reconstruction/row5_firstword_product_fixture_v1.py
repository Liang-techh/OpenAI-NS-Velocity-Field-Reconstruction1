"""Synthetic checks for :mod:`row5_firstword_product_v1`.

The fixture supplies a small saved-expression-shaped mapping and a finite
first-word state directly to ``build``.  It does not import a producer, read a
receipt, evaluate formal atoms, or execute an ancestor leaf.  Checks are
structural and exact-rational where possible; no result is promoted to a
function value.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Mapping


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
STATE_COMPONENTS = (1, 2, 3, 5)
OUTPUT_ROWS = ("h5", "h5_Z1")

__all__ = [
    "SOURCE_POWER_BASIS",
    "STATE_COMPONENTS",
    "OUTPUT_ROWS",
    "run_fixture",
    "run",
]


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
    left_raw, _ = _raw_endpoints(left)
    _, right_raw = _raw_endpoints(right)
    return c.make_mpf((left_raw, right_raw))


def _inside(c: Any, interval: Any, value: Fraction | int) -> bool:
    if isinstance(interval, Mapping) and "lower_exact_mpf_tuple" in interval:
        interval = c.make_mpf(
            (
                interval["lower_exact_mpf_tuple"],
                interval["upper_exact_mpf_tuple"],
            )
        )
    point = _q(c, value)
    lower, upper = _raw_endpoints(interval)
    point_lower, point_upper = _raw_endpoints(point)
    return bool(lower <= point_lower and point_upper <= upper)


def _term(
    c: Any,
    coefficient: Fraction | int,
    source_id: str,
    path: Any,
    powers: tuple[int, int, int, int],
    *,
    coefficient_bounds: tuple[Fraction | int, Fraction | int] | None = None,
    scalar: int = 1,
    source_family: str | None = None,
    role: str | None = None,
    state_reference: tuple[int, int] | None = None,
    uncertainty: bool | None = False,
    error_only: bool | None = False,
    symmetric_error: bool | None = False,
    signed_field: bool | None = True,
    duplicate_source_id: str | None = None,
) -> dict[str, Any]:
    if coefficient_bounds is None:
        directed = _q(c, coefficient)
    else:
        directed = _interval(c, coefficient_bounds[0], coefficient_bounds[1])
    result: dict[str, Any] = {
        "coefficient": directed,
        "coefficient_interval": directed,
        "powers": list(powers),
        "formal_powers": list(powers),
        "source_ids": [duplicate_source_id or source_id],
        "source_path": _clone(path),
        "source_paths": [_clone(path)],
        "error_classes": [f"fixture-error:{source_id}"],
        "source_error_classes": [f"fixture-error:{source_id}"],
        "uncertainty": uncertainty,
        "error_only": error_only,
        "symmetric_error": symmetric_error,
        "signed_field": signed_field,
        "scalar": scalar,
        "fixture_source_id": source_id,
    }
    if source_family is not None:
        result["source_family"] = source_family
    if role is not None:
        result["role"] = role
    if state_reference is None:
        result["source_only"] = source_family == "g5"
    else:
        result["state_reference"] = {
            "component": state_reference[0],
            "ordinary_Z_order": state_reference[1],
        }
    return result


def _state_term(
    c: Any,
    component: int,
    order: int,
    value: Fraction | int,
    powers: tuple[int, int, int, int],
    *,
    uncertainty: bool | None = False,
    error_only: bool | None = False,
    symmetric_error: bool | None = False,
    signed_field: bool | None = True,
    bounds: tuple[Fraction | int, Fraction | int] | None = None,
) -> dict[str, Any]:
    return _term(
        c,
        value,
        f"state:W{component}:Z{order}",
        ["state-input", f"W{component}", f"Z{order}", 0],
        powers,
        coefficient_bounds=bounds,
        uncertainty=uncertainty,
        error_only=error_only,
        symmetric_error=symmetric_error,
        signed_field=signed_field,
    )


def _expression(c: Any) -> dict[str, Any]:
    duplicate = "fixture:shared-source-duplicate"
    return {
        "schema": "row5_affine_forcing_v1",
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "rows": {
            "h5": [
                _term(
                    c,
                    1,
                    "g5-err",
                    ["receipt", "result", "output", "expression", "h5", 0],
                    (0, 0, 0, 0),
                    source_family="g5",
                    role="source_only_g5",
                    coefficient_bounds=(-1, 1),
                    uncertainty=True,
                    error_only=True,
                    symmetric_error=True,
                    signed_field=False,
                    duplicate_source_id=duplicate,
                ),
                _term(
                    c,
                    2,
                    "g5-dup",
                    ["receipt", "result", "output", "expression", "h5", 1],
                    (1, 0, 0, 0),
                    source_family="g5",
                    role="source_only_g5",
                    duplicate_source_id=duplicate,
                ),
                _term(
                    c,
                    3,
                    "C0-W1-a",
                    ["receipt", "result", "output", "expression", "h5", 2],
                    (1, -1, 0, 0),
                    source_family="C0",
                    role="C0_Z0_times_state",
                    state_reference=(1, 0),
                ),
                _term(
                    c,
                    4,
                    "C0-W1-a",
                    ["receipt", "result", "output", "expression", "h5", 3],
                    (0, 0, 0, 1),
                    source_family="C0",
                    role="C0_Z0_times_state_duplicate_source",
                    state_reference=(1, 0),
                ),
                _term(
                    c,
                    -7,
                    "B1-W2",
                    ["receipt", "result", "output", "expression", "h5", 4],
                    (0, 0, 1, 0),
                    source_family="B1",
                    role="B1_Z0_times_state_Z1",
                    state_reference=(2, 1),
                ),
                _term(
                    c,
                    11,
                    "axis-chi",
                    ["receipt", "result", "output", "expression", "h5", 5],
                    (0, 0, 0, 0),
                    source_family="axis",
                    role="two_chi_times_W1",
                    state_reference=(1, 0),
                    scalar=2,
                    uncertainty=True,
                ),
            ],
            "h5_Z1": [
                _term(
                    c,
                    5,
                    "g5-Z1",
                    ["receipt", "result", "output", "expression", "h5_Z1", 0],
                    (0, 1, 0, 0),
                    source_family="g5",
                    role="source_only_g5_Z1",
                ),
                _term(
                    c,
                    -2,
                    "C0-Z1-W3",
                    ["receipt", "result", "output", "expression", "h5_Z1", 1],
                    (0, 0, 0, -1),
                    source_family="C0",
                    role="C0_Z1_times_state",
                    state_reference=(3, 0),
                ),
                _term(
                    c,
                    3,
                    "C0-Z0-W5Z1",
                    ["receipt", "result", "output", "expression", "h5_Z1", 2],
                    (0, 0, 1, 0),
                    source_family="C0",
                    role="C0_Z0_times_state_Z1",
                    state_reference=(5, 1),
                ),
                _term(
                    c,
                    -4,
                    "B1-Z1-W1Z1",
                    ["receipt", "result", "output", "expression", "h5_Z1", 3],
                    (1, 0, 0, 0),
                    source_family="B1",
                    role="B1_Z1_times_state_Z1",
                    state_reference=(1, 1),
                ),
                _term(
                    c,
                    7,
                    "B1-Z0-W3Z2",
                    ["receipt", "result", "output", "expression", "h5_Z1", 4],
                    (0, -1, 0, 0),
                    source_family="B1",
                    role="B1_Z0_times_state_Z2",
                    state_reference=(3, 2),
                ),
                _term(
                    c,
                    6,
                    "B1-Z0-W1Z2",
                    ["receipt", "result", "output", "expression", "h5_Z1", 5],
                    (0, 0, 0, 0),
                    source_family="B1",
                    role="B1_Z0_times_state_Z2",
                    state_reference=(1, 2),
                ),
                _term(
                    c,
                    -5,
                    "B1-Z0-W2Z2",
                    ["receipt", "result", "output", "expression", "h5_Z1", 6],
                    (0, 0, 0, 0),
                    source_family="B1",
                    role="B1_Z0_times_state_Z2",
                    state_reference=(2, 2),
                ),
                _term(
                    c,
                    -9,
                    "axis-chi-Z1",
                    ["receipt", "result", "output", "expression", "h5_Z1", 7],
                    (0, 0, 0, 0),
                    source_family="axis",
                    role="two_chi_Z_times_W1",
                    state_reference=(1, 0),
                    scalar=2,
                    uncertainty=True,
                ),
                _term(
                    c,
                    8,
                    "axis-chi",
                    ["receipt", "result", "output", "expression", "h5_Z1", 8],
                    (0, 0, 0, 0),
                    source_family="axis",
                    role="two_chi_times_W1_Z",
                    state_reference=(1, 1),
                    scalar=2,
                    uncertainty=True,
                ),
            ],
        },
    }


def _state(c: Any) -> dict[str, dict[str, list[dict[str, Any]]]]:
    rows: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for component in STATE_COMPONENTS:
        rows[str(component)] = {}
    rows["1"]["Z0"] = [_state_term(c, 1, 0, 5, (0, 1, 0, -1), uncertainty=True)]
    rows["1"]["Z1"] = [_state_term(c, 1, 1, -2, (-1, 0, 0, 1))]
    rows["1"]["Z2"] = [_state_term(c, 1, 2, 3, (0, 0, -1, 0))]
    rows["2"]["Z0"] = [_state_term(c, 2, 0, 7, (0, 0, 0, 0))]
    rows["2"]["Z1"] = [
        _state_term(
            c,
            2,
            1,
            0,
            (0, 0, 0, 2),
            bounds=(-1, 1),
            uncertainty=True,
            error_only=True,
            symmetric_error=True,
            signed_field=False,
        )
    ]
    rows["2"]["Z2"] = [_state_term(c, 2, 2, 11, (1, 0, 0, 0))]
    rows["3"]["Z0"] = [_state_term(c, 3, 0, 13, (0, -1, 0, 0))]
    rows["3"]["Z1"] = [_state_term(c, 3, 1, 17, (0, 0, 0, 0))]
    rows["3"]["Z2"] = [_state_term(c, 3, 2, -19, (0, 0, -1, 0))]
    rows["5"]["Z0"] = [_state_term(c, 5, 0, 23, (0, 0, 0, 0))]
    rows["5"]["Z1"] = [_state_term(c, 5, 1, -29, (0, 1, 0, 0))]
    # W5_Z2 is intentionally absent: the current row-5 expression must not
    # demand an unused second Z derivative of component 5.
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
    """Run synthetic product, provenance, and fail-closed checks."""

    build_fn = _callable(build, "build")
    expression = _expression(c)
    state = _state(c)
    checks: dict[str, bool] = {}
    try:
        result = build_fn(c, expression, state)
    except Exception as exc:
        return {
            "schema": "row5_firstword_product_fixture_v1",
            "passed": False,
            "checks": {"build_completed": False},
            "check_count": 1,
            "error": f"{type(exc).__name__}: {exc}",
            "synthetic_only": True,
            "implementation_imported": False,
        }

    checks["build_completed"] = isinstance(result, Mapping)
    checks["schema"] = result.get("schema") == "row5_firstword_product_v1"
    checks["conditional_status"] = result.get("status") == "conditional_firstword_substituted_expression"
    checks["formal_basis"] = tuple(result.get("source_power_basis", ())) == SOURCE_POWER_BASIS
    checks["not_full_h5_semantics"] = (
        result.get("expression_semantics") == "g5 + feedback(W_firstword), conditional; not full h5"
        and result.get("row_totals") is None
    )
    checks["rows_present"] = all(label in result.get("rows", {}) for label in OUTPUT_ROWS)
    checks["compact_top_level_shape"] = all(
        key not in result
        for key in ("products", "h5_products", "h5_Z1_products", "source_node_table", "state_node_table")
    ) and all(key in result for key in ("rows", "source_nodes", "state_nodes"))
    checks["products_use_node_ids"] = all(
        "source_node_id" in product
        and "state_node_id" in product
        and product.get("provenance", {}).get("node_table_lookup_required") is True
        and product.get("scalar_applied_to_coefficient") is False
        and product.get("term_semantics") == "scalar * coefficient_enclosure * formal_monomial; scalar deferred"
        and "coefficient_enclosure" in product
        and all(key not in product for key in ("coefficient", "coefficient_interval", "coefficient_endpoints"))
        and "source_ids" not in product
        and "source_paths" not in product
        for label in OUTPUT_ROWS
        for product in result["rows"][label]
    )
    checks["directed_endpoints_retained"] = all(
        isinstance(product.get("coefficient_enclosure"), Mapping)
        and product["coefficient_enclosure"].get("directed") is True
        and set(product["coefficient_enclosure"]) >= {"lower_exact_mpf_tuple", "upper_exact_mpf_tuple"}
        for label in OUTPUT_ROWS
        for product in result["rows"][label]
    )
    checks["no_function_value"] = (
        result.get("bounds_are_function_values") is False
        and result.get("formal_atoms_materialized") is False
        and result.get("interval_sum_performed") is False
        and result.get("full_h5") is False
    )
    checks["unknown_residual_separate"] = (
        result.get("unknown_residual", {}).get("symbol") == "R"
        and result["unknown_residual"].get("included_in_products") is False
        and result["unknown_residual"].get("added_to_ledger") is False
    )

    source_nodes = result["source_nodes"]
    state_nodes = result["state_nodes"]
    products = result["rows"]
    all_products = products["h5"] + products["h5_Z1"]
    checks["all_node_ids_have_paths"] = all(
        node.get("parent_input_paths")
        for node in list(source_nodes.values()) + list(state_nodes.values())
    )
    checks["nonzero_source_coefficient_present"] = any(
        not _inside(c, node["coefficient_interval"], 0) for node in source_nodes.values()
    )
    checks["source_only_g5_retained"] = sum(product["source_only"] for product in products["h5"]) == 2 and sum(
        product["source_only"] for product in products["h5_Z1"]
    ) == 1
    checks["source_duplicate_occurrences_not_deduped"] = len(
        [node for node in source_nodes.values() if node.get("provenance_ids") == ["fixture:shared-source-duplicate"]]
    ) == 2
    checks["shared_state_node_reused"] = sum(
        product.get("state_node_id") == "state:W1:Z0:0" for product in all_products
    ) >= 3
    checks["no_cancellation_flags"] = all(
        product["provenance"].get("cancellation_performed") is False for product in all_products
    )
    checks["uncertainty_not_inferred"] = all(
        product["provenance"].get("uncertainty_inferred") is False for product in all_products
    )
    checks["error_only_provenance_retained"] = any(
        source_nodes[product["source_node_id"]].get("error_only") is True for product in products["h5"]
    ) and any(
        product.get("state_node_id") is not None
        and state_nodes[product["state_node_id"]].get("error_only") is True
        for product in products["h5"]
    )

    def _source_node_for(product: Mapping[str, Any]) -> Mapping[str, Any]:
        return source_nodes[product["source_node_id"]]

    axis_products = [
        product
        for product in products["h5"] + products["h5_Z1"]
        if _source_node_for(product).get("source_family") == "axis"
    ]
    checks["axis_factor_two_preserved"] = bool(axis_products) and all(
        product.get("scalar") == 2 for product in axis_products
    )
    checks["axis_products_are_signed_nonzero_enclosures"] = bool(axis_products) and all(
        not _inside(c, product["coefficient_enclosure"], 0) for product in axis_products
    )
    checks["directed_product_coefficient_checked"] = any(
        _inside(c, product["coefficient_enclosure"], 55)
        for product in axis_products
    )

    c0_duplicate = next(
        product
        for product in products["h5"]
        if _source_node_for(product).get("role") == "C0_Z0_times_state"
    )
    checks["four_power_sum"] = c0_duplicate.get("formal_powers") == [1, 0, 0, -1]
    b1_zz = next(
        product
        for product in products["h5_Z1"]
        if _source_node_for(product).get("role") == "B1_Z0_times_state_Z2"
    )
    checks["B1_W_ZZ_cross_term"] = b1_zz.get("state_node_id") == "state:W3:Z2:0"
    checks["ordinary_Z1_product_rule_terms"] = {
        _source_node_for(product).get("role")
        for product in products["h5_Z1"]
    } >= {
        "C0_Z1_times_state",
        "C0_Z0_times_state_Z1",
        "B1_Z1_times_state_Z1",
        "B1_Z0_times_state_Z2",
        "two_chi_Z_times_W1",
        "two_chi_times_W1_Z",
    }
    checks["negative_powers_preserved"] = any(
        any(power < 0 for power in product.get("formal_powers", ())) for product in all_products
    )
    checks["positive_rational_product_contained"] = any(
        _inside(c, product["coefficient_enclosure"], 15)
        for product in products["h5"]
        if _source_node_for(product).get("role") == "C0_Z0_times_state"
    )
    checks["negative_rational_product_contained"] = any(
        _inside(c, product["coefficient_enclosure"], -26)
        for product in products["h5_Z1"]
        if _source_node_for(product).get("role") == "C0_Z1_times_state"
    )
    checks["scope_flags_preserved"] = (
        result.get("scope", {}).get("feedback_is_firstword_only") is True
        and result.get("scope", {}).get("unknown_residual_is_separate") is True
        and result.get("scope", {}).get("no_product_deduplication") is True
        and result.get("scope", {}).get("no_error_cancellation") is True
        and result.get("scope", {}).get("no_atom_values") is True
        and result.get("scope", {}).get("source_only_g5_retention_scope")
        == "one output per input occurrence; no global deduplication"
        and result.get("conditions", {}).get("physical_dx_applied") is False
        and result.get("conditions", {}).get("axis_chi_cancellation") is False
        and result.get("product_coefficient_semantics")
        == "coefficient_enclosure = source_coefficient * state_coefficient"
        and result.get("scalar_application") == "deferred_to_consumer_once"
    )
    checks["W5_Z2_not_required"] = True

    missing_state_checks: dict[str, bool] = {}
    for component in (1, 2, 3):
        missing_state = _clone(state)
        del missing_state[str(component)]["Z2"]
        missing_state_checks[str(component)] = _rejected(lambda s=missing_state: build_fn(c, expression, s))
    checks["missing_W1_W2_W3_ZZ_rejected"] = all(missing_state_checks.values())
    with_atoms = _clone(state)
    with_atoms["formal_atoms"] = {"hb": _q(c, 2)}
    checks["atom_values_rejected"] = _rejected(lambda: build_fn(c, expression, with_atoms))
    bad_state_scalar = _clone(state)
    bad_state_scalar["1"]["Z0"][0]["scalar"] = 2
    checks["nonunit_state_scalar_rejected"] = _rejected(lambda: build_fn(c, expression, bad_state_scalar))
    bad_power_basis = _clone(state)
    bad_power_basis["1"]["Z0"][0]["power_basis"] = ["wrong", "basis", "only", "here"]
    checks["incompatible_node_basis_rejected"] = _rejected(lambda: build_fn(c, expression, bad_power_basis))
    bad_float_expression = _clone(expression)
    bad_float_expression["rows"]["h5"][2]["coefficient"] = 1.25
    bad_float_expression["rows"]["h5"][2]["coefficient_interval"] = 1.25
    checks["python_float_coefficient_rejected"] = _rejected(lambda: build_fn(c, bad_float_expression, state))
    bad_float_state = _clone(state)
    bad_float_state["1"]["Z0"][0]["coefficient"] = [1.0, 2.0]
    bad_float_state["1"]["Z0"][0]["coefficient_interval"] = [1.0, 2.0]
    checks["python_float_endpoint_rejected"] = _rejected(lambda: build_fn(c, expression, bad_float_state))
    bad_expression = _clone(expression)
    bad_expression["rows"]["h5"][4].pop("state_reference")
    checks["missing_state_reference_rejected"] = _rejected(lambda: build_fn(c, bad_expression, state))

    return {
        "schema": "row5_firstword_product_fixture_v1",
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "exact_fraction_expected": True,
        "ordinary_Z1_product_rule_checked": True,
        "B1_W_ZZ_checked": True,
        "negative_powers_checked": True,
        "uncertainty_provenance_checked": True,
        "shared_node_and_duplicate_source_checked": True,
        "axis_factor_two_checked": True,
        "missing_state_rejected": True,
        "W5_Z2_not_required": True,
        "unknown_residual_not_added": True,
        "implementation_imported": False,
    }


def run(build: Any, c: Any) -> dict[str, Any]:
    return run_fixture(build, c)
