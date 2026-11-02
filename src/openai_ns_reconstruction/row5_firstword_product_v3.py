"""Version-3 compact first-word product consumer for the saved row-5 formula.

This is a standalone consumer.  It accepts the saved v1 source expression and
a caller-owned finite first-word state, retains source/state provenance in
node tables, and emits only directed coefficient enclosures plus formal-power
products.  Formal atoms, midpoints, row sums, and the unknown residual are not
materialized.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "row5_firstword_product_v3"
COMPONENT = "conditional_row5_firstword_product"
SOURCE_SCHEMA = "row5_affine_forcing_v1"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
STATE_COMPONENTS = (1, 2, 3, 5)
STATE_ORDERS = (0, 1, 2)
OUTPUT_ROWS = ("h5", "h5_Z1")

_MISSING = object()
_COEFFICIENT_KEYS = ("coefficient_interval", "coefficient", "coefficient_enclosure", "directed_interval", "enclosure", "interval")
_POWER_KEYS = ("formal_powers", "powers", "source_powers")
_ID_KEYS = ("source_ids", "source_id", "state_ids", "state_id", "input_id", "term_id")

__all__ = [
    "SCHEMA",
    "COMPONENT",
    "SOURCE_SCHEMA",
    "SOURCE_POWER_BASIS",
    "STATE_COMPONENTS",
    "STATE_ORDERS",
    "OUTPUT_ROWS",
    "build",
    "build_firstword_product",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _clone(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if _is_sequence(value):
        return [_clone(item) for item in value]
    return value


def _validate_context(c: Any) -> None:
    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must provide mpf and make_mpf directed interval operations")


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


def _endpoint_payload(value: Any) -> dict[str, Any]:
    lower, upper = _raw_endpoints(value)
    return {"lower": _clone(lower), "upper": _clone(upper)}


def _enclosure_payload(value: Any) -> dict[str, Any]:
    lower, upper = _raw_endpoints(value)
    return {
        "lower_exact_mpf_tuple": _clone(lower),
        "upper_exact_mpf_tuple": _clone(upper),
        "directed": True,
        "function_value": False,
    }


def _reject_python_float(value: Any, label: str) -> None:
    if isinstance(value, float):
        raise TypeError(f"{label} may not contain a Python float")
    if isinstance(value, Mapping):
        for key, nested in value.items():
            _reject_python_float(nested, f"{label}.{key}")
    elif _is_sequence(value):
        for index, nested in enumerate(value):
            _reject_python_float(nested, f"{label}[{index}]")


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    return c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))


def _exact_mpf_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    if not _is_sequence(value) or len(value) != 4:
        raise TypeError(f"{label} must be a four-entry exact MPF tuple")
    if any(type(item) is not int for item in value):
        raise TypeError(f"{label} must contain integers only")
    sign, mantissa, exponent, bitcount = (int(item) for item in value)
    if sign not in (0, 1) or mantissa < 0 or bitcount < 0:
        raise ValueError(f"{label} has invalid exact MPF tuple fields")
    if mantissa == 0 and bitcount != 0:
        raise ValueError(f"{label} is non-finite")
    return sign, mantissa, exponent, bitcount


def _interval(c: Any, value: Any, label: str) -> Any:
    _reject_python_float(value, label)
    if isinstance(value, Mapping):
        nested = next((value[key] for key in _COEFFICIENT_KEYS if key in value), _MISSING)
        if nested is not _MISSING and nested is not value:
            return _interval(c, nested, label + ".coefficient")
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            if "lower_exact_mpf_tuple" not in value or "upper_exact_mpf_tuple" not in value:
                raise TypeError(f"{label} requires both exact MPF endpoints")
            lower = _exact_mpf_tuple(value["lower_exact_mpf_tuple"], label + ".lower_exact_mpf_tuple")
            upper = _exact_mpf_tuple(value["upper_exact_mpf_tuple"], label + ".upper_exact_mpf_tuple")
            result = c.make_mpf((lower, upper))
        else:
            lower = next((value[key] for key in ("lower", "lo", "left") if key in value), _MISSING)
            upper = next((value[key] for key in ("upper", "hi", "right") if key in value), _MISSING)
            if lower is _MISSING or upper is _MISSING:
                raise TypeError(f"{label} is not interval-compatible")
            result = _make_interval(c, lower, upper)
    elif _is_sequence(value) and len(value) == 2:
        result = _make_interval(c, value[0], value[1])
    else:
        try:
            result = c.mpf(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{label} is not interval-compatible") from exc
    lower, upper = _raw_endpoints(result)
    if c.make_mpf((upper, upper)) < c.make_mpf((lower, lower)):
        raise ValueError(f"{label} has reversed directed endpoints")
    return result


def _powers(raw: Mapping[Any, Any], label: str) -> list[int]:
    value = next((raw[key] for key in _POWER_KEYS if key in raw), _MISSING)
    if value is _MISSING:
        raise KeyError(f"{label}.formal_powers is required")
    if isinstance(value, Mapping):
        value = [value.get(name, _MISSING) for name in SOURCE_POWER_BASIS]
    if not _is_sequence(value) or len(value) != 4:
        raise TypeError(f"{label}.formal_powers must contain four integers")
    if any(isinstance(item, bool) or not isinstance(item, int) for item in value):
        raise TypeError(f"{label}.formal_powers must contain integers only")
    return [int(item) for item in value]


def _items(value: Any) -> list[Any]:
    if value is _MISSING or value is None:
        return []
    return list(value) if _is_sequence(value) else [value]


def _metadata_values(raw: Mapping[Any, Any], names: Sequence[str]) -> list[Any]:
    values: list[Any] = []
    for name in names:
        if name in raw:
            values.extend(_clone(item) for item in _items(raw[name]))
    for wrapper in ("origin", "metadata", "provenance", "original_term"):
        nested = raw.get(wrapper)
        if isinstance(nested, Mapping):
            for name in names:
                if name in nested:
                    values.extend(_clone(item) for item in _items(nested[name]))
    return values


def _paths(raw: Mapping[Any, Any]) -> list[Any]:
    values: list[Any] = []
    singular = ("source_path", "state_path", "input_path", "path", "parent_receipt_path")
    plural = ("source_paths", "state_paths", "input_paths", "paths")
    for name in singular:
        if name in raw:
            values.append(_clone(raw[name]))
    for name in plural:
        if name in raw:
            values.extend(_clone(item) for item in _items(raw[name]))
    for wrapper in ("origin", "metadata", "provenance", "original_term"):
        nested = raw.get(wrapper)
        if not isinstance(nested, Mapping):
            continue
        for name in singular:
            if name in nested:
                values.append(_clone(nested[name]))
        for name in plural:
            if name in nested:
                values.extend(_clone(item) for item in _items(nested[name]))
    return values


def _bool_or_none(raw: Mapping[Any, Any], name: str) -> bool | None:
    value = raw.get(name, _MISSING)
    if value is _MISSING:
        for wrapper in ("origin", "metadata", "provenance", "original_term"):
            nested = raw.get(wrapper)
            if isinstance(nested, Mapping) and name in nested:
                value = nested[name]
                break
    if value is _MISSING:
        return None
    if not isinstance(value, bool):
        raise TypeError(f"{name} must be boolean when supplied")
    return value


def _normalise_node(c: Any, raw: Any, label: str, *, node_kind: str, node_id: str) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must be a mapping")
    if {"formal_atoms", "atom_values", "formal_atom_values", "atoms"}.intersection(raw):
        raise ValueError(f"{label} may not contain materialized formal atoms")
    explicit_basis = raw.get("power_basis", raw.get("formal_power_basis", _MISSING))
    if explicit_basis is not _MISSING and tuple(explicit_basis) != SOURCE_POWER_BASIS:
        raise ValueError(f"{label}.power_basis is not the four-atom source basis")
    coefficient_raw = next((raw[key] for key in _COEFFICIENT_KEYS if key in raw), _MISSING)
    if coefficient_raw is _MISSING:
        raise KeyError(f"{label}.coefficient enclosure is required")
    coefficient = _interval(c, coefficient_raw, label + ".coefficient")
    powers = _powers(raw, label)
    identifiers = _metadata_values(raw, _ID_KEYS)
    if not identifiers:
        raise ValueError(f"{label} must retain at least one source/state node ID")
    paths = [value for value in _paths(raw) if value is not None and value != []]
    if not paths:
        raise ValueError(f"{label} must retain a parent receipt/input path")
    scalar = raw.get("scalar", 1)
    if isinstance(scalar, bool) or not isinstance(scalar, int):
        raise TypeError(f"{label}.scalar must be an integer")
    if node_kind == "finite_firstword_state" and scalar != 1:
        raise ValueError(f"{label}.scalar must be 1; state scalar is not absorbed")
    result = _clone(dict(raw))
    result.update(
        {
            "node_id": str(node_id),
            "node_kind": node_kind,
            "coefficient": coefficient,
            "coefficient_interval": coefficient,
            "coefficient_endpoints": _endpoint_payload(coefficient),
            "formal_powers": list(powers),
            "powers": list(powers),
            "power_basis": list(SOURCE_POWER_BASIS),
            "scalar": int(scalar),
            "provenance_ids": identifiers,
            "parent_input_paths": paths,
            "source_ids": identifiers,
            "source_paths": paths,
            "uncertainty": _bool_or_none(raw, "uncertainty"),
            "error_only": _bool_or_none(raw, "error_only"),
            "symmetric_error": _bool_or_none(raw, "symmetric_error"),
            "signed_field": _bool_or_none(raw, "signed_field"),
            "error_classes": _metadata_values(raw, ("error_classes", "error_class", "uncertainty_class")),
            "source_error_classes": _metadata_values(raw, ("source_error_classes", "source_error_class")),
            "original_node": _clone(dict(raw)),
            "midpoint_values_used": False,
            "function_value": False,
        }
    )
    if not result["source_error_classes"]:
        result["source_error_classes"] = list(result["error_classes"])
    return result


def _expression_rows(expression: Mapping[str, Any]) -> dict[str, Sequence[Any]]:
    rows = expression.get("rows")
    if not isinstance(rows, Mapping):
        rows = expression
    result: dict[str, Sequence[Any]] = {}
    for label in OUTPUT_ROWS:
        values = rows.get(label, _MISSING)
        if values is _MISSING or not _is_sequence(values):
            raise KeyError(f"expression is missing additive row {label}")
        result[label] = values
    return result


def _reference(raw: Mapping[str, Any], label: str) -> tuple[int, int] | None:
    value = raw.get("state_reference", raw.get("state_ref", _MISSING))
    if value is _MISSING or value is None:
        return None
    if not isinstance(value, Mapping):
        raise TypeError(f"{label}.state_reference must be a mapping or null")
    component = value.get("component", _MISSING)
    order = value.get("ordinary_Z_order", value.get("order", _MISSING))
    if isinstance(component, str) and component.isdigit():
        component = int(component)
    if isinstance(order, str) and order.isdigit():
        order = int(order)
    if isinstance(component, bool) or not isinstance(component, int) or component not in STATE_COMPONENTS:
        raise ValueError(f"{label} has invalid state component reference")
    if isinstance(order, bool) or not isinstance(order, int) or order not in STATE_ORDERS:
        raise ValueError(f"{label} has invalid state ordinary-Z reference")
    return component, order


def _state_component(state: Mapping[Any, Any], component: int) -> Any:
    for key in (str(component), component, f"W{component}", f"W_{component}"):
        if key in state:
            return state[key]
    raise KeyError(f"state is missing W{component}")


def _state_terms(state: Mapping[Any, Any], component: int, order: int) -> list[Any]:
    rows = _state_component(state, component)
    if not isinstance(rows, Mapping):
        raise TypeError(f"state W{component} must be an order mapping")
    value = _MISSING
    for key in (f"Z{order}", order, str(order), f"ordinary_Z{order}"):
        if key in rows:
            value = rows[key]
            break
    if value is _MISSING:
        raise KeyError(f"state is missing W{component}_Z{order}")
    if isinstance(value, Mapping) and "terms" in value:
        value = value["terms"]
    values = _items(value)
    if not values:
        raise ValueError(f"state W{component}_Z{order} must contain at least one term")
    return values


def _order_from_key(key: Any) -> int | None:
    if key == 0 or (isinstance(key, str) and key.lower() in {"0", "z0", "ordinary_z0"}):
        return 0
    if key == 1 or (isinstance(key, str) and key.lower() in {"1", "z1", "ordinary_z1"}):
        return 1
    if key == 2 or (isinstance(key, str) and key.lower() in {"2", "z2", "ordinary_z2"}):
        return 2
    return None


def _index_state_nodes(c: Any, state: Mapping[Any, Any]) -> dict[str, dict[str, Any]]:
    table: dict[str, dict[str, Any]] = {}
    for component in STATE_COMPONENTS:
        try:
            rows = _state_component(state, component)
        except KeyError:
            continue
        if not isinstance(rows, Mapping):
            raise TypeError(f"state W{component} must be an order mapping")
        for raw_key in rows:
            order = _order_from_key(raw_key)
            if order is None:
                continue
            for index, raw_state in enumerate(_state_terms(state, component, order)):
                node_id = f"state:W{component}:Z{order}:{index}"
                if node_id in table:
                    raise ValueError(f"state node {node_id} is supplied more than once")
                node = _normalise_node(c, raw_state, f"state.W{component}.Z{order}[{index}]", node_kind="finite_firstword_state", node_id=node_id)
                node["state_component"] = component
                node["state_ordinary_Z_order"] = order
                node["used_in_products"] = False
                table[node_id] = node
    return table


def _multiply(c: Any, left: Any, right: Any, label: str) -> Any:
    return _interval(c, _interval(c, left, label + ".source") * _interval(c, right, label + ".state"), label + ".product")


def _product_record(c: Any, row: str, source_node_id: str, state_node_id: str | None, source_node: Mapping[str, Any], state_node: Mapping[str, Any] | None, source_index: int, state_index: int | None) -> dict[str, Any]:
    coefficient = _interval(c, source_node["coefficient_interval"], f"{row}[{source_index}].source_only") if state_node is None else _multiply(c, source_node["coefficient_interval"], state_node["coefficient_interval"], f"{row}[{source_index}]x{state_index}")
    state_powers = [0, 0, 0, 0] if state_node is None else state_node["formal_powers"]
    return {
        "product_id": f"product:{row}:{source_index}:{state_node_id or 'source-only'}",
        "row": row,
        "source_node_id": source_node_id,
        "state_node_id": state_node_id,
        "source_only": state_node is None,
        "state_reference": None if state_node is None else _clone(source_node.get("state_reference")),
        "scalar": int(source_node.get("scalar", 1)),
        "formal_powers": [int(a) + int(b) for a, b in zip(source_node["formal_powers"], state_powers)],
        "coefficient_enclosure": _enclosure_payload(coefficient),
        "coefficient_is_function_value": False,
        "scalar_applied_to_coefficient": False,
        "term_semantics": "scalar * coefficient_enclosure * formal_monomial; scalar deferred",
        "directed": True,
        "source_index": int(source_index),
        "state_index": None if state_index is None else int(state_index),
        "provenance": {"source_node_id": source_node_id, "state_node_id": state_node_id, "node_table_lookup_required": True, "uncertainty_inferred": False, "sign_inferred": False, "cancellation_performed": False},
        "unknown_residual_included": False,
        "midpoint_values_used": False,
    }


def build(c: Any, expression: Mapping[str, Any], state: Mapping[Any, Any]) -> dict[str, Any]:
    _validate_context(c)
    if not isinstance(expression, Mapping):
        raise TypeError("expression must be a saved row5 expression mapping")
    if expression.get("schema") not in (None, SOURCE_SCHEMA):
        raise ValueError("expression schema is not row5_affine_forcing_v1")
    if not isinstance(state, Mapping):
        raise TypeError("state must be a direct W-component mapping")
    if {"formal_atoms", "atom_values", "formal_atom_values", "atoms"}.intersection(state):
        raise ValueError("formal atom values are forbidden in first-word substitution")
    if tuple(expression.get("source_power_basis", SOURCE_POWER_BASIS)) != SOURCE_POWER_BASIS:
        raise ValueError("expression has an incompatible four-atom formal basis")
    rows = _expression_rows(expression)
    source_nodes: dict[str, dict[str, Any]] = {}
    state_nodes = _index_state_nodes(c, state)
    products: dict[str, list[dict[str, Any]]] = {label: [] for label in OUTPUT_ROWS}
    required_refs: set[tuple[int, int]] = set()
    for row in OUTPUT_ROWS:
        for source_index, raw_source in enumerate(rows[row]):
            source_node_id = f"expression:{row}:{source_index}"
            source_node = _normalise_node(c, raw_source, f"expression.{row}[{source_index}]", node_kind="source_expression", node_id=source_node_id)
            source_nodes[source_node_id] = source_node
            reference = _reference(raw_source, f"expression.{row}[{source_index}]")
            if reference is None:
                if raw_source.get("source_family") not in (None, "g5") and not raw_source.get("source_only", False):
                    raise ValueError(f"{row}[{source_index}] has no state reference but is not source-only g5")
                products[row].append(_product_record(c, row, source_node_id, None, source_node, None, source_index, None))
                continue
            required_refs.add(reference)
            component, order = reference
            for state_index, raw_state in enumerate(_state_terms(state, component, order)):
                state_node_id = f"state:W{component}:Z{order}:{state_index}"
                state_node = state_nodes.get(state_node_id)
                if state_node is None:
                    raise KeyError(f"state node {state_node_id} was not indexed")
                state_node["used_in_products"] = True
                products[row].append(_product_record(c, row, source_node_id, state_node_id, source_node, state_node, source_index, state_index))
    for component, order in sorted(required_refs):
        _state_terms(state, component, order)
    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "conditional_firstword_substituted_expression",
        "source_schema": SOURCE_SCHEMA,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "rows": products,
        "expression_semantics": "g5 + feedback(W_firstword), conditional; not full h5",
        "row_totals": None,
        "source_nodes": source_nodes,
        "state_nodes": state_nodes,
        "required_state_references": [{"component": component, "ordinary_Z_order": order} for component, order in sorted(required_refs)],
        "source_duplicates_preserved": True,
        "same_state_node_reused_without_dedup": True,
        "formal_powers_materialized": False,
        "formal_atoms_materialized": False,
        "product_coefficient_semantics": "coefficient_enclosure = source_coefficient * state_coefficient",
        "scalar_application": "deferred_to_consumer_once",
        "interval_sum_performed": False,
        "bounds_are_function_values": False,
        "physical_h5_materialized": False,
        "full_h5": False,
        "full_state": False,
        "state_solved": False,
        "ns_validation": False,
        "unknown_residual": {"symbol": "R", "definition": "W - W_finite", "preserved": True, "included_in_products": False, "added_to_ledger": False},
        "conditions": {
            "finite_firstword_only": True,
            "fixed_x": "caller/source expression condition; not proved here",
            "rho_equals_x_squared": "caller premise; no chain rule applied",
            "geometry_Z_independent": "caller premise; not proved here",
            "same_frame": "caller condition; not proved here",
            "same_scale": "caller condition; not proved here",
            "same_gauge": "caller condition; not proved here",
            "live_owner": "caller condition; not proved here",
            "same_state_owner": "caller condition; not proved here",
            "same_state_Z0_Z1_Z2_owner": "caller condition; not proved here",
            "physical_dx_applied": False,
            "extra_F0_or_Lambda_transform": False,
            "axis_chi_cancellation": False,
        },
        "scope": {
            "source_only_g5_retained_once": True,
            "source_only_g5_retention_scope": "one output per input occurrence; no global deduplication",
            "feedback_is_firstword_only": True,
            "unknown_residual_is_separate": True,
            "no_product_deduplication": True,
            "no_error_cancellation": True,
            "no_atom_values": True,
            "no_function_value": True,
        },
        "resource_limits": {"source_node_count": len(source_nodes), "state_node_count": len(state_nodes), "h5_product_count": len(products["h5"]), "h5_Z1_product_count": len(products["h5_Z1"])},
    }


build_firstword_product = build
