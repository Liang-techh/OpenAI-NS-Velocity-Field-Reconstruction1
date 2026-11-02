"""Independent synthetic fixture for :mod:`row5_baseline_product_v4`.

The expression intentionally contains every state reference used by the
row-5 formula, including W5 Z0/Z1 and all three W1/W2/W3 Z2 terms.  A
separate reduced-expression check verifies that an allowed but unreferenced
state column is not made mandatory.  Values and provenance IDs are fresh for
version 4 and include signed intervals, second derivatives, deferred scalar,
and error metadata.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Mapping


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
OUTPUT_ROWS = ("h5", "h5_Z1")
STATE_SECTORS = ("core", "local")
ALL_REFERENCES = (
    (1, 0), (2, 0), (3, 0), (5, 0),
    (1, 1), (2, 1), (3, 1), (5, 1),
    (1, 2), (2, 2), (3, 2),
)

__all__ = ["SOURCE_POWER_BASIS", "OUTPUT_ROWS", "STATE_SECTORS", "ALL_REFERENCES", "run_fixture", "run"]


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


def _endpoints(value: Any) -> tuple[Any, Any]:
    return _raw_endpoint(value, False), _raw_endpoint(value, True)


def _interval(c: Any, lower: Fraction | int, upper: Fraction | int) -> Any:
    return c.make_mpf((_raw_endpoint(_q(c, lower), False), _raw_endpoint(_q(c, upper), True)))


def _point(c: Any, value: Any) -> Any:
    endpoint = tuple(value) if isinstance(value, (list, tuple)) else value
    return c.make_mpf((endpoint, endpoint))


def _inside(c: Any, enclosure: Any, value: Fraction | int) -> bool:
    if isinstance(enclosure, Mapping) and "lower_exact_mpf_tuple" in enclosure:
        lower_raw = tuple(enclosure["lower_exact_mpf_tuple"])
        upper_raw = tuple(enclosure["upper_exact_mpf_tuple"])
    else:
        lower_raw, upper_raw = _endpoints(enclosure)
    lower = _point(c, lower_raw)
    upper = _point(c, upper_raw)
    candidate = _q(c, value)
    return bool(lower <= candidate.a and candidate.b <= upper)


def _origin(prefix: str, *, error_only: bool = False) -> dict[str, Any]:
    return {
        "source_ids": [prefix],
        "source_path": ["synthetic-v4", prefix],
        "error_classes": [f"v4-error:{prefix}"],
        "source_error_classes": [f"v4-error:{prefix}"],
        "uncertainty": bool(error_only),
        "error_only": bool(error_only),
        "symmetric_error": bool(error_only),
        "signed_field": not error_only,
    }


def _term(
    c: Any,
    coefficient: Fraction | int,
    prefix: str,
    powers: tuple[int, int, int, int],
    *,
    bounds: tuple[Fraction | int, Fraction | int] | None = None,
    scalar: int = 1,
    family: str | None = None,
    role: str | None = None,
    reference: tuple[int, int] | None = None,
    error_only: bool = False,
) -> dict[str, Any]:
    directed = _q(c, coefficient) if bounds is None else _interval(c, bounds[0], bounds[1])
    result: dict[str, Any] = {
        "coefficient": directed,
        "coefficient_interval": directed,
        "formal_powers": list(powers),
        "powers": list(powers),
        "source_ids": [prefix],
        "source_path": ["synthetic-v4", prefix],
        "source_paths": [["synthetic-v4", prefix]],
        "input_id": prefix + ":input",
        "origin": _origin(prefix, error_only=error_only),
        "error_classes": [f"v4-error:{prefix}"],
        "source_error_classes": [f"v4-error:{prefix}"],
        "uncertainty": bool(error_only),
        "error_only": bool(error_only),
        "symmetric_error": bool(error_only),
        "signed_field": not error_only,
        "scalar": scalar,
    }
    if family is not None:
        result["source_family"] = family
    if role is not None:
        result["role"] = role
    if reference is not None:
        result["state_reference"] = {"component": reference[0], "ordinary_Z_order": reference[1]}
    elif family == "g5":
        result["source_only_g5"] = True
    return result


def _expression(c: Any) -> dict[str, Any]:
    # Every allowed affine state reference appears at least once.  The rows
    # retain the formula's C0/B1 product-rule pieces and both axis pieces.
    h5 = [
        _term(c, 31, "v4:g5:zero", (0, 0, 0, 0), family="g5", role="source_only_g5"),
        _term(c, -3, "v4:C0:W1Z0", (1, -1, 0, 2), bounds=(-3, -2), family="C0", role="C0_Z0_times_state", reference=(1, 0)),
        _term(c, 4, "v4:C0:W2Z0", (0, 1, 0, -1), family="C0", role="C0_Z0_times_state", reference=(2, 0)),
        _term(c, -5, "v4:C0:W3Z0", (0, 0, 1, 0), family="C0", role="C0_Z0_times_state", reference=(3, 0)),
        _term(c, 6, "v4:C0:W5Z0", (0, 0, 0, 1), family="C0", role="C0_Z0_times_state", reference=(5, 0)),
        _term(c, -7, "v4:B1:W1Z1", (1, 0, 0, 0), family="B1", role="B1_Z0_times_state_Z1", reference=(1, 1)),
        _term(c, 8, "v4:B1:W2Z1", (0, -1, 0, 1), family="B1", role="B1_Z0_times_state_Z1", reference=(2, 1)),
        _term(c, -9, "v4:B1:W3Z1", (0, 0, -1, 0), family="B1", role="B1_Z0_times_state_Z1", reference=(3, 1)),
        _term(c, -11, "v4:axis:chi", (0, 0, 0, 0), scalar=2, family="axis", role="two_chi_times_W1", reference=(1, 0)),
    ]
    h5_z1 = [
        _term(c, -37, "v4:g5:one", (0, 1, 0, 0), family="g5", role="source_only_g5_Z1"),
        _term(c, -2, "v4:C0Z1:W1Z0", (0, 0, 1, 0), family="C0", role="C0_Z1_times_state", reference=(1, 0)),
        _term(c, 3, "v4:C0Z1:W2Z0", (0, 0, 0, -1), family="C0", role="C0_Z1_times_state", reference=(2, 0)),
        _term(c, -4, "v4:C0Z1:W3Z0", (1, 0, 0, 0), family="C0", role="C0_Z1_times_state", reference=(3, 0)),
        _term(c, 5, "v4:C0Z1:W5Z0", (0, 0, 1, 1), family="C0", role="C0_Z1_times_state", reference=(5, 0)),
        _term(c, -6, "v4:C0Z0:W1Z1", (0, 1, 0, 0), family="C0", role="C0_Z0_times_state_Z1", reference=(1, 1)),
        _term(c, 7, "v4:C0Z0:W2Z1", (0, 0, 0, 0), family="C0", role="C0_Z0_times_state_Z1", reference=(2, 1)),
        _term(c, -8, "v4:C0Z0:W3Z1", (0, -1, 1, 0), family="C0", role="C0_Z0_times_state_Z1", reference=(3, 1)),
        _term(c, 9, "v4:C0Z0:W5Z1", (0, 0, -1, 0), family="C0", role="C0_Z0_times_state_Z1", reference=(5, 1)),
        _term(c, -10, "v4:B1Z1:W1Z1", (1, 0, 0, 0), family="B1", role="B1_Z1_times_state_Z1", reference=(1, 1)),
        _term(c, 11, "v4:B1Z1:W2Z1", (0, 0, 1, 0), family="B1", role="B1_Z1_times_state_Z1", reference=(2, 1)),
        _term(c, -12, "v4:B1Z1:W3Z1", (0, 0, 0, -1), family="B1", role="B1_Z1_times_state_Z1", reference=(3, 1)),
        _term(c, -13, "v4:B1Z0:W1Z2", (0, 0, 0, 0), family="B1", role="B1_Z0_times_state_Z2", reference=(1, 2)),
        _term(c, 14, "v4:B1Z0:W2Z2", (0, 1, -1, 0), family="B1", role="B1_Z0_times_state_Z2", reference=(2, 2)),
        _term(c, -15, "v4:B1Z0:W3Z2", (0, 0, 1, -1), family="B1", role="B1_Z0_times_state_Z2", reference=(3, 2), error_only=True),
        _term(c, 16, "v4:axis:chiZ", (0, 0, 0, 0), scalar=2, family="axis", role="two_chi_Z_times_W1", reference=(1, 0)),
        _term(c, -18, "v4:axis:chi", (0, 0, 0, 0), scalar=2, family="axis", role="two_chi_times_W1_Z", reference=(1, 1)),
    ]
    return {"schema": "row5_affine_forcing_v1", "source_power_basis": list(SOURCE_POWER_BASIS), "rows": {"h5": h5, "h5_Z1": h5_z1}}


def _state(c: Any) -> dict[str, dict[tuple[int, int], dict[str, list[dict[str, Any]]]]]:
    state: dict[str, dict[tuple[int, int], dict[str, list[dict[str, Any]]]]] = {}
    for cell_index, cell in enumerate(("v4-cell-left", "v4-cell-right")):
        state[cell] = {}
        for component, order in ALL_REFERENCES:
            prefix = f"v4:{cell}:W{component}:Z{order}"
            if (component, order) == (1, 0):
                core, local = 5, -4
            elif (component, order) == (3, 2):
                core, local = -6, 3
            elif order == 2:
                core, local = component + 2, -(component + 3)
            else:
                core, local = component + order + 2, -(component + order + 3)
            state[cell][(component, order)] = {
                "core": [_term(c, core, prefix + ":core", (0, component - 2, order - 1, cell_index), family="W0", role="W0_core_history", error_only=(component, order) == (3, 2))],
                "local": [_term(c, local, prefix + ":local", (cell_index, order, 2 - component, 0), family="W0", role="W0_local")],
            }
    return state


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
        return {"schema": "row5_baseline_product_fixture_v4", "passed": False, "checks": {"build_completed": False}, "check_count": 1, "error": f"{type(exc).__name__}: {exc}", "synthetic_only": True, "implementation_imported": False}

    checks: dict[str, bool] = {}
    checks["build_completed"] = isinstance(result, Mapping)
    checks["schema"] = result.get("schema") == "row5_baseline_product_v4"
    checks["source_schema"] = result.get("source_schema") == "row5_affine_forcing_v1"
    checks["required_refs"] = {
        (item.get("component"), item.get("ordinary_Z_order"))
        for item in result.get("required_state_references", [])
    } == set(ALL_REFERENCES)
    checks["rows_have_core_local"] = all(
        set(result.get("rows", {}).get(row, {}).get(cell, {})) == {"core", "local"}
        for row in OUTPUT_ROWS for cell in state
    )
    source_nodes = result.get("source_nodes", {})
    products = result.get("rows", {})
    all_products = [
        product
        for row in OUTPUT_ROWS
        for cell_rows in products.get(row, {}).values()
        for sector in STATE_SECTORS
        for product in cell_rows.get(sector, [])
    ]
    checks["g5_unique"] = len(result.get("skipped_g5", [])) == 2 and all(source_nodes[p["source_node_id"]].get("source_family") != "g5" for p in all_products)
    checks["families_retained"] = {source_nodes[p["source_node_id"]].get("source_family") for p in all_products} == {"C0", "B1", "axis"}
    checks["core_local_separate"] = all(products[row][cell]["core"] and products[row][cell]["local"] for row in OUTPUT_ROWS for cell in state)
    checks["negative_interval_product"] = any(
        source_nodes[p["source_node_id"]].get("source_id") == "v4:C0:W1Z0"
        and p.get("sector") == "core"
        and _inside(c, p["coefficient_enclosure"], -15)
        for p in all_products
    )
    checks["second_derivative"] = any(p.get("state_component") == 3 and p.get("state_ordinary_Z_order") == 2 for p in all_products)
    checks["only_W1_W2_W3_Z2"] = {
        (p.get("state_component"), p.get("state_ordinary_Z_order"))
        for p in all_products if p.get("state_ordinary_Z_order") == 2
    } == {(1, 2), (2, 2), (3, 2)}
    axis_products = [p for p in all_products if p.get("source_family") == "axis"]
    checks["scalar_deferred_once"] = bool(axis_products) and all(p.get("scalar") == 2 and p.get("scalar_applied_to_coefficient") is False for p in axis_products)
    checks["scalar_not_multiplied"] = any(_inside(c, p["coefficient_enclosure"], -55) for p in axis_products)
    checks["powers_added"] = any(
        source_nodes[p["source_node_id"]].get("source_id") == "v4:C0:W1Z0"
        and p.get("formal_powers") == [1, -2, -1, 2]
        for p in all_products
    )
    checks["bidirectional_provenance"] = all(
        p.get("source_provenance", {}).get("node_id") == p.get("source_node_id")
        and p.get("state_provenance", {}).get("node_id") == p.get("state_node_id")
        and "error_provenance" in p
        for p in all_products
    )
    checks["error_separate"] = any(
        p.get("error_provenance", {}).get("source", {}).get("error_only") is True
        and p.get("error_provenance", {}).get("mixed_or_cancelled") is False
        for p in all_products
    )
    checks["no_cancellation"] = all(
        p.get("provenance", {}).get("endpoint_cancellation_performed") is False
        and p.get("provenance", {}).get("error_cancellation_performed") is False
        for p in all_products
    )
    checks["R_hi_unknown"] = result.get("R_hi", {}).get("symbol") == "R_hi" and result["R_hi"].get("included_in_products") is False
    checks["not_summed_or_materialized"] = result.get("interval_sum_performed") is False and result.get("formal_atoms_materialized") is False

    missing = _clone(state)
    del missing["v4-cell-left"][(3, 2)]
    checks["missing_required_rejected"] = _rejected(lambda: build_fn(c, expression, missing))
    bad_w5_z2 = _clone(state)
    bad_w5_z2["v4-cell-left"][(5, 2)] = bad_w5_z2["v4-cell-left"][(1, 2)]
    checks["W5_Z2_rejected"] = _rejected(lambda: build_fn(c, expression, bad_w5_z2))
    with_atoms = _clone(state)
    with_atoms["formal_atoms"] = {"hb": _q(c, 2)}
    checks["atom_guard"] = _rejected(lambda: build_fn(c, expression, with_atoms))
    bad_scalar = _clone(state)
    bad_scalar["v4-cell-left"][(1, 0)]["core"][0]["scalar"] = 2
    checks["state_scalar_rejected"] = _rejected(lambda: build_fn(c, expression, bad_scalar))

    # Dynamic-reference guard: remove the C0_Z0*W5_Z1 source from a copy of
    # the complete expression.  The full state still contains W5 Z1, but it
    # must not become a mandatory reference merely because it is allowed.
    reduced = _clone(expression)
    reduced["rows"]["h5_Z1"] = [
        term for term in reduced["rows"]["h5_Z1"]
        if term.get("state_reference") != {"component": 5, "ordinary_Z_order": 1}
        or term.get("role") != "C0_Z0_times_state_Z1"
    ]
    reduced_state = _clone(state)
    del reduced_state["v4-cell-left"][(5, 1)]
    del reduced_state["v4-cell-right"][(5, 1)]
    try:
        reduced_result = build_fn(c, reduced, reduced_state)
        reduced_refs = {
            (item.get("component"), item.get("ordinary_Z_order"))
            for item in reduced_result.get("required_state_references", [])
        }
        checks["unreferenced_state_not_required"] = (5, 1) not in reduced_refs
    except Exception:
        checks["unreferenced_state_not_required"] = False

    return {
        "schema": "row5_baseline_product_fixture_v4",
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "fresh_values_and_ids": True,
        "required_reference_coverage_checked": True,
        "dynamic_reference_guard_checked": True,
        "negative_values_checked": True,
        "ordinary_Z2_checked": True,
        "g5_unique_retention_checked": True,
        "scalar_once_checked": True,
        "error_provenance_checked": True,
        "implementation_imported": False,
    }


def run(build: Any, c: Any) -> dict[str, Any]:
    return run_fixture(build, c)

