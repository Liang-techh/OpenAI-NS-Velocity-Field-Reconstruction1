"""Independent row-5 baseline ``W0`` product leaf, version 4.

The input is a saved ``row5_affine_forcing_v1`` expression.  The state is a
cell-indexed mapping with the explicit shape::

    state[cell][(component, ordinary_Z_order)] = {
        "core": [term, ...], "local": [term, ...]
    }

The leaf emits source/state products without summing them.  Standalone
``g5`` records are skipped because they belong to the first-word product;
all ``C0``/``B1`` records and the explicit ``2*chi`` records remain.  Core
and local W0 history are separate output branches, directed endpoints and
the four formal powers are retained, and the scalar is deferred exactly
once.  ``R_hi`` remains an unknown caller-owned remainder.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "row5_baseline_product_v4"
COMPONENT = "conditional_row5_baseline_W0_product_v4"
SOURCE_SCHEMA = "row5_affine_forcing_v1"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
OUTPUT_ROWS = ("h5", "h5_Z1")
STATE_COMPONENTS = (1, 2, 3, 5)
STATE_ORDERS = (0, 1, 2)
STATE_SECTORS = ("core", "local")
ALLOWED_STATE_REFERENCES = frozenset(
    {
        (1, 0), (2, 0), (3, 0), (5, 0),
        (1, 1), (2, 1), (3, 1), (5, 1),
        (1, 2), (2, 2), (3, 2),
    }
)

_MISSING = object()
_COEFFICIENT_KEYS = (
    "coefficient_interval", "source_interval", "coefficient", "raw_interval",
    "directed_interval", "enclosure", "interval", "box", "value",
)
_POWER_KEYS = (
    "formal_powers", "powers", "source_powers", "formal_source_powers",
    "source_formal_powers",
)
_ID_KEYS = (
    "source_ids", "source_id", "state_ids", "state_id", "input_ids",
    "input_id", "term_id", "node_id", "id",
)
_PATH_KEYS = (
    "source_paths", "source_path", "state_paths", "state_path", "input_paths",
    "input_path", "parent_input_paths", "path", "paths",
)
_ATOM_KEYS = {"formal_atoms", "atom_values", "formal_atom_values", "atoms"}

__all__ = [
    "SCHEMA", "COMPONENT", "SOURCE_SCHEMA", "SOURCE_POWER_BASIS",
    "OUTPUT_ROWS", "STATE_COMPONENTS", "STATE_ORDERS", "STATE_SECTORS",
    "ALLOWED_STATE_REFERENCES", "build", "build_baseline_product",
    "build_row5_baseline_product", "row5_baseline_product",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _items(value: Any) -> list[Any]:
    if value is _MISSING or value is None:
        return []
    return list(value) if _is_sequence(value) else [value]


def _clone(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if _is_sequence(value):
        return [_clone(item) for item in value]
    return value


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


def _reject_float(value: Any, label: str) -> None:
    if isinstance(value, float):
        raise TypeError(f"{label} may not contain a Python float")
    if isinstance(value, Mapping):
        for key, nested in value.items():
            _reject_float(nested, f"{label}.{key}")
    elif _is_sequence(value):
        for index, nested in enumerate(value):
            _reject_float(nested, f"{label}[{index}]")


def _exact_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
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


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    return c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))


def _interval(c: Any, value: Any, label: str) -> Any:
    _reject_float(value, label)
    if isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            if "lower_exact_mpf_tuple" not in value or "upper_exact_mpf_tuple" not in value:
                raise TypeError(f"{label} requires both exact MPF endpoints")
            lower = _exact_tuple(value["lower_exact_mpf_tuple"], label + ".lower_exact_mpf_tuple")
            upper = _exact_tuple(value["upper_exact_mpf_tuple"], label + ".upper_exact_mpf_tuple")
            result = c.make_mpf((lower, upper))
        else:
            nested = next((value[key] for key in _COEFFICIENT_KEYS if key in value), _MISSING)
            if nested is not _MISSING and nested is not value:
                return _interval(c, nested, label + ".nested")
            lower = next((value[key] for key in ("lower", "lo", "left") if key in value), _MISSING)
            upper = next((value[key] for key in ("upper", "hi", "right") if key in value), _MISSING)
            if lower is _MISSING or upper is _MISSING:
                raise TypeError(f"{label} is not interval-compatible")
            result = _make_interval(c, lower, upper)
    elif _is_sequence(value) and len(value) == 2:
        result = _make_interval(c, value[0], value[1])
    else:
        try:
            point = c.mpf(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{label} is not interval-compatible") from exc
        result = point if getattr(point, "_mpi_", None) is not None else _make_interval(c, point, point)
    lower, upper = _raw_endpoints(result)
    if c.make_mpf((upper, upper)) < c.make_mpf((lower, lower)):
        raise ValueError(f"{label} has reversed directed endpoints")
    return result


def _endpoint_payload(value: Any) -> dict[str, Any]:
    lower, upper = _raw_endpoints(value)
    return {
        "lower_exact_mpf_tuple": _clone(lower),
        "upper_exact_mpf_tuple": _clone(upper),
        "directed": True,
        "function_value": False,
    }


def _powers(raw: Mapping[Any, Any], label: str) -> list[int]:
    value = next((raw[key] for key in _POWER_KEYS if key in raw), _MISSING)
    if value is _MISSING:
        for wrapper_name in ("origin", "metadata", "provenance", "original_term"):
            wrapper = raw.get(wrapper_name)
            if isinstance(wrapper, Mapping):
                value = next((wrapper[key] for key in _POWER_KEYS if key in wrapper), _MISSING)
                if value is not _MISSING:
                    break
    if value is _MISSING:
        raise KeyError(f"{label}.formal_powers is required")
    if isinstance(value, Mapping):
        value = [value.get(atom, _MISSING) for atom in SOURCE_POWER_BASIS]
    if not _is_sequence(value) or len(value) != 4:
        raise TypeError(f"{label}.formal_powers must contain four integers")
    if any(isinstance(item, bool) or not isinstance(item, int) for item in value):
        raise TypeError(f"{label}.formal_powers must contain integers only")
    return [int(item) for item in value]


def _metadata(raw: Mapping[Any, Any], names: Sequence[str]) -> list[Any]:
    values: list[Any] = []
    for name in names:
        if name in raw:
            values.extend(_clone(item) for item in _items(raw[name]))
    for wrapper_name in ("origin", "metadata", "provenance", "original_term"):
        wrapper = raw.get(wrapper_name)
        if isinstance(wrapper, Mapping):
            for name in names:
                if name in wrapper:
                    values.extend(_clone(item) for item in _items(wrapper[name]))
    return values


def _paths(raw: Mapping[Any, Any]) -> list[Any]:
    values: list[Any] = []
    for name in _PATH_KEYS:
        if name in raw:
            if name.endswith("_paths") or name in {"paths", "parent_input_paths"}:
                values.extend(_clone(item) for item in _items(raw[name]))
            else:
                values.append(_clone(raw[name]))
    for wrapper_name in ("origin", "metadata", "provenance", "original_term"):
        wrapper = raw.get(wrapper_name)
        if isinstance(wrapper, Mapping):
            for name in _PATH_KEYS:
                if name in wrapper:
                    if name.endswith("_paths") or name in {"paths", "parent_input_paths"}:
                        values.extend(_clone(item) for item in _items(wrapper[name]))
                    else:
                        values.append(_clone(wrapper[name]))
    return values


def _bool(raw: Mapping[Any, Any], name: str) -> bool | None:
    value = raw.get(name, _MISSING)
    if value is _MISSING and name == "uncertainty":
        value = raw.get("is_uncertainty", _MISSING)
    if value is _MISSING:
        for wrapper_name in ("origin", "metadata", "provenance", "original_term"):
            wrapper = raw.get(wrapper_name)
            if isinstance(wrapper, Mapping) and name in wrapper:
                value = wrapper[name]
                break
    if value is _MISSING:
        return None
    if not isinstance(value, bool):
        raise TypeError(f"{name} must be boolean when supplied")
    return value


def _origin(raw: Mapping[Any, Any], label: str, required: bool) -> Any:
    value = raw.get("origin", _MISSING)
    if value is _MISSING:
        value = raw.get("provenance", raw.get("metadata", _MISSING))
    if value is _MISSING:
        if required:
            raise KeyError(f"{label}.origin is required")
        return _clone(raw.get("original_term", raw))
    return _clone(value)


def _node(
    c: Any,
    raw: Any,
    label: str,
    *,
    node_kind: str,
    node_id: str,
    origin_required: bool,
    path_required: bool,
    state_scalar: bool,
) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must be a term mapping")
    if _ATOM_KEYS.intersection(raw):
        raise ValueError(f"{label} may not contain materialized formal atoms")
    coefficient_raw = next((raw[key] for key in _COEFFICIENT_KEYS if key in raw), _MISSING)
    if coefficient_raw is _MISSING:
        raise KeyError(f"{label}.source_interval/coefficient is required")
    coefficient = _interval(c, coefficient_raw, label + ".coefficient")
    powers = _powers(raw, label)
    ids = _metadata(raw, _ID_KEYS)
    if not ids:
        raise ValueError(f"{label} must retain at least one provenance ID")
    paths = [item for item in _paths(raw) if item is not None and item != []]
    if path_required and not paths:
        raise ValueError(f"{label} must retain a provenance/input path")
    scalar = raw.get("scalar", 1)
    if isinstance(scalar, bool) or not isinstance(scalar, int):
        raise TypeError(f"{label}.scalar must be an integer")
    if state_scalar and scalar != 1:
        raise ValueError(f"{label}.scalar must be 1; state scalar is not absorbed")
    errors = _metadata(raw, ("error_classes", "error_class", "uncertainty_class"))
    source_errors = _metadata(raw, ("source_error_classes", "source_error_class"))
    if not source_errors:
        source_errors = list(errors)
    result = _clone(dict(raw))
    result.update(
        {
            "node_id": str(node_id),
            "node_kind": node_kind,
            "coefficient": coefficient,
            "coefficient_interval": coefficient,
            "source_interval": coefficient,
            "coefficient_endpoints": _endpoint_payload(coefficient),
            "source_interval_endpoints": _endpoint_payload(coefficient),
            "formal_powers": list(powers),
            "powers": list(powers),
            "power_basis": list(SOURCE_POWER_BASIS),
            "scalar": int(scalar),
            "provenance_ids": ids,
            "source_ids": list(ids),
            "source_id": ids[0],
            "parent_input_paths": paths,
            "source_paths": list(paths),
            "origin": _origin(raw, label, origin_required),
            "uncertainty": _bool(raw, "uncertainty"),
            "error_only": _bool(raw, "error_only"),
            "symmetric_error": _bool(raw, "symmetric_error"),
            "signed_field": _bool(raw, "signed_field"),
            "error_classes": errors,
            "source_error_classes": source_errors,
            "original_term": _clone(dict(raw)),
            "midpoint_values_used": False,
            "function_value": False,
        }
    )
    return result


def _reference(raw: Mapping[str, Any], label: str) -> tuple[int, int] | None:
    value = next((raw[key] for key in ("state_reference", "state_ref", "state", "reference") if key in raw), _MISSING)
    if value is _MISSING or value is None:
        return None
    if not isinstance(value, Mapping):
        raise TypeError(f"{label}.state_reference must be a mapping or null")
    component = value.get("component", value.get("state_component", _MISSING))
    order = value.get("ordinary_Z_order", value.get("ordinary_order", value.get("order", _MISSING)))
    if isinstance(component, str) and component.isdigit():
        component = int(component)
    if isinstance(order, str) and order.isdigit():
        order = int(order)
    if isinstance(component, bool) or not isinstance(component, int) or component not in STATE_COMPONENTS:
        raise ValueError(f"{label} has invalid state component reference")
    if isinstance(order, bool) or not isinstance(order, int) or order not in STATE_ORDERS:
        raise ValueError(f"{label} has invalid state ordinary-Z reference")
    reference = (int(component), int(order))
    if reference not in ALLOWED_STATE_REFERENCES:
        raise ValueError(f"{label} requests unsupported baseline state reference {reference}")
    return reference


def _state_key(key: Any, label: str) -> tuple[int, int]:
    if not isinstance(key, tuple) or len(key) != 2:
        raise TypeError(f"{label} must use tuple (component, ordinary_order) keys")
    component, order = key
    if isinstance(component, bool) or not isinstance(component, int) or component not in STATE_COMPONENTS:
        raise ValueError(f"{label} has invalid state component")
    if isinstance(order, bool) or not isinstance(order, int) or order not in STATE_ORDERS:
        raise ValueError(f"{label} has invalid ordinary-Z order")
    reference = (int(component), int(order))
    if reference not in ALLOWED_STATE_REFERENCES:
        raise ValueError(f"{label} requests unsupported baseline state reference {reference}")
    return reference


def _is_g5(raw: Mapping[str, Any]) -> bool:
    family = raw.get("source_family", raw.get("family", ""))
    if isinstance(family, str) and family.lower() == "g5":
        return True
    if raw.get("source_only_g5") is True:
        return True
    role = raw.get("role", raw.get("term_role", ""))
    return isinstance(role, str) and "g5" in role.lower() and "source" in role.lower()


def _state_index(c: Any, state: Mapping[Any, Any]) -> tuple[list[dict[str, Any]], dict[tuple[int, tuple[int, int], str], list[dict[str, Any]]], dict[str, dict[str, Any]]]:
    if not isinstance(state, Mapping):
        raise TypeError("state must be a cell mapping")
    if _ATOM_KEYS.intersection(state):
        raise ValueError("formal atom values are forbidden in baseline state")
    if not state:
        raise ValueError("state must contain at least one cell")
    cells: list[dict[str, Any]] = []
    lookup: dict[tuple[int, tuple[int, int], str], list[dict[str, Any]]] = {}
    table: dict[str, dict[str, Any]] = {}
    for cell_index, (cell, payload) in enumerate(state.items()):
        if not isinstance(payload, Mapping):
            raise TypeError(f"state[{cell!r}] must be a mapping")
        info = {"index": cell_index, "cell": _clone(cell), "references": {}}
        cells.append(info)
        for raw_ref, sectors in payload.items():
            reference = _state_key(raw_ref, f"state[{cell!r}]")
            if not isinstance(sectors, Mapping):
                raise TypeError(f"state[{cell!r}][{reference!r}] must be a sector mapping")
            if "core" not in sectors or "local" not in sectors:
                raise KeyError(f"state[{cell!r}][{reference!r}] requires core and local")
            component, order = reference
            ref_label = f"W{component}:Z{order}"
            info["references"][ref_label] = {}
            for sector in STATE_SECTORS:
                values = _items(sectors[sector])
                if not values:
                    raise ValueError(f"state[{cell!r}][{reference!r}][{sector}] must not be empty")
                nodes: list[dict[str, Any]] = []
                for term_index, raw_term in enumerate(values):
                    node_id = f"state:{cell_index}:W{component}:Z{order}:{sector}:{term_index}"
                    node = _node(
                        c,
                        raw_term,
                        f"state[{cell!r}][{reference!r}][{sector}][{term_index}]",
                        node_kind="baseline_W0_state",
                        node_id=node_id,
                        origin_required=True,
                        path_required=False,
                        state_scalar=True,
                    )
                    node.update(
                        {
                            "cell": _clone(cell),
                            "cell_index": cell_index,
                            "state_component": component,
                            "state_ordinary_Z_order": order,
                            "state_sector": sector,
                            "state_term_index": term_index,
                            "used_in_products": False,
                        }
                    )
                    if node_id in table:
                        raise ValueError(f"state node {node_id} is duplicated")
                    table[node_id] = node
                    nodes.append(node)
                lookup[(cell_index, reference, sector)] = nodes
                info["references"][ref_label][sector] = [node["node_id"] for node in nodes]
    return cells, lookup, table


def _multiply(c: Any, left: Any, right: Any, label: str) -> Any:
    return _interval(c, _interval(c, left, label + ".source") * _interval(c, right, label + ".state"), label + ".product")


def _provenance(node: Mapping[str, Any]) -> dict[str, Any]:
    # Full source and uncertainty metadata stay in the referenced node table.
    return {"node_id": node["node_id"], "node_table_lookup_required": True}


def _product(c: Any, row: str, index: int, source: Mapping[str, Any], state: Mapping[str, Any], source_index: int, state_index: int) -> dict[str, Any]:
    coefficient = _multiply(c, source["coefficient_interval"], state["coefficient_interval"], f"{row}[{index}]")
    powers = [int(a) + int(b) for a, b in zip(source["formal_powers"], state["formal_powers"])]
    return {
        "product_id": f"product:{row}:{index}",
        "row": row,
        "cell": _clone(state["cell"]),
        "cell_index": int(state["cell_index"]),
        "sector": state["state_sector"],
        "source_node_id": source["node_id"],
        "state_node_id": state["node_id"],
        "source_index": int(source_index),
        "state_index": int(state_index),
        "source_family": source.get("source_family"),
        "source_role": source.get("role", source.get("term_role")),
        "source_only": False,
        "source_only_g5": False,
        "state_reference": _clone(source.get("state_reference")),
        "state_component": int(state["state_component"]),
        "state_ordinary_Z_order": int(state["state_ordinary_Z_order"]),
        "scalar": int(source.get("scalar", 1)),
        "scalar_applied_to_coefficient": False,
        "term_semantics": "scalar * coefficient_enclosure * formal_monomial; scalar deferred",
        "coefficient_enclosure": _endpoint_payload(coefficient),
        "coefficient_is_function_value": False,
        "formal_powers": powers,
        "powers": list(powers),
        "power_basis": list(SOURCE_POWER_BASIS),
        "directed": True,
        "function_value": False,
        "source_provenance": _provenance(source),
        "state_provenance": _provenance(state),
        "error_provenance": {
            "source": {
                "uncertainty": source.get("uncertainty"),
                "error_only": source.get("error_only"),
                "error_classes": _clone(source.get("error_classes", [])),
                "source_error_classes": _clone(source.get("source_error_classes", [])),
            },
            "state": {
                "uncertainty": state.get("uncertainty"),
                "error_only": state.get("error_only"),
                "error_classes": _clone(state.get("error_classes", [])),
                "source_error_classes": _clone(state.get("source_error_classes", [])),
            },
            "mixed_or_cancelled": False,
        },
        "provenance": {
            "source_node_id": source["node_id"],
            "state_node_id": state["node_id"],
            "source_and_state_retained": True,
            "node_table_lookup_required": True,
            "uncertainty_inferred": False,
            "sign_inferred": False,
            "endpoint_cancellation_performed": False,
            "error_cancellation_performed": False,
        },
        "unknown_residual_included": False,
        "midpoint_values_used": False,
        "formal_atoms_materialized": False,
    }


def build(c: Any, expression: Mapping[str, Any], state: Mapping[Any, Any]) -> dict[str, Any]:
    """Build un-summed baseline products for the affine expression.

    Required state references are derived only from the expression's actual
    ``state_reference`` records.  Extra allowed state columns are indexed and
    retained but never made mandatory.  Z2 is accepted only for W1, W2 and
    W3; W5 has only Z0 and Z1.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must provide mpf and make_mpf directed interval operations")
    if not isinstance(expression, Mapping):
        raise TypeError("expression must be a saved row5 affine expression mapping")
    if expression.get("schema") not in (None, SOURCE_SCHEMA):
        raise ValueError("expression schema is not row5_affine_forcing_v1")
    if tuple(expression.get("source_power_basis", SOURCE_POWER_BASIS)) != SOURCE_POWER_BASIS:
        raise ValueError("expression has an incompatible four-atom formal basis")
    if _ATOM_KEYS.intersection(expression):
        raise ValueError("formal atom values are forbidden in baseline product input")
    rows = expression.get("rows", expression)
    if not isinstance(rows, Mapping):
        raise TypeError("expression.rows must be a mapping")
    for row in OUTPUT_ROWS:
        if row not in rows or not _is_sequence(rows[row]):
            raise KeyError(f"expression is missing additive row {row}")

    cells, lookup, state_nodes = _state_index(c, state)
    source_nodes: dict[str, dict[str, Any]] = {}
    occurrences: dict[str, list[tuple[int, dict[str, Any], tuple[int, int] | None, bool]]] = {row: [] for row in OUTPUT_ROWS}
    required: set[tuple[int, int]] = set()
    skipped: list[dict[str, Any]] = []
    for row in OUTPUT_ROWS:
        for source_index, raw_source in enumerate(rows[row]):
            if not isinstance(raw_source, Mapping):
                raise TypeError(f"expression.{row}[{source_index}] must be a term mapping")
            node_id = f"expression:{row}:{source_index}"
            source = _node(
                c,
                raw_source,
                f"expression.{row}[{source_index}]",
                node_kind="row5_affine_source",
                node_id=node_id,
                origin_required=False,
                path_required=True,
                state_scalar=False,
            )
            reference = _reference(raw_source, f"expression.{row}[{source_index}]")
            source["state_reference"] = _clone(reference)
            source["state_ref"] = _clone(reference)
            source["output_row"] = row
            source["source_index"] = source_index
            source_nodes[node_id] = source
            if reference is None:
                if not _is_g5(raw_source):
                    raise ValueError(f"expression.{row}[{source_index}] has no state reference and is not g5")
                skipped.append(
                    {
                        "source_node_id": node_id,
                        "row": row,
                        "source_index": source_index,
                        "source_family": raw_source.get("source_family", raw_source.get("family")),
                        "retained_by": "row5_firstword_product_v4",
                        "product_emitted": False,
                    }
                )
                occurrences[row].append((source_index, source, None, True))
                continue
            required.add(reference)
            occurrences[row].append((source_index, source, reference, False))

    output: dict[str, dict[Any, dict[str, list[dict[str, Any]]]]] = {
        row: {info["cell"]: {sector: [] for sector in STATE_SECTORS} for info in cells}
        for row in OUTPUT_ROWS
    }
    counts = {row: {sector: 0 for sector in STATE_SECTORS} for row in OUTPUT_ROWS}
    product_index = 0
    for row in OUTPUT_ROWS:
        for source_index, source, reference, is_skipped in occurrences[row]:
            if is_skipped or reference is None:
                continue
            for info in cells:
                cell_index = int(info["index"])
                for sector in STATE_SECTORS:
                    terms = lookup.get((cell_index, reference, sector))
                    if not terms:
                        raise KeyError(f"state cell {info['cell']!r} is missing {reference!r} {sector} sector")
                    for state_index, state_node in enumerate(terms):
                        state_node["used_in_products"] = True
                        output[row][info["cell"]][sector].append(
                            _product(c, row, product_index, source, state_node, source_index, state_index)
                        )
                        counts[row][sector] += 1
                        product_index += 1

    references = [
        {"component": component, "ordinary_Z_order": order}
        for component, order in sorted(required)
    ]
    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "conditional_row5_baseline_W0_products_v4",
        "source_schema": SOURCE_SCHEMA,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "output_rows": list(OUTPUT_ROWS),
        "rows": output,
        "row_totals": None,
        "source_nodes": source_nodes,
        "state_nodes": state_nodes,
        "cells": cells,
        "state_shape": "state[cell][(component, ordinary_order)] = {'core': [terms], 'local': [terms]}",
        "state_sectors": list(STATE_SECTORS),
        "required_state_references": references,
        "allowed_state_references": [
            {"component": component, "ordinary_Z_order": order}
            for component, order in sorted(ALLOWED_STATE_REFERENCES)
        ],
        "skipped_g5": skipped,
        "g5_consumed_by_firstword": None,
        "g5_retention_in_firstword": "caller premise; prior pinned receipt documents it",
        "g5_repeated_here": False,
        "source_duplicates_preserved": True,
        "state_nodes_reused_without_dedup": True,
        "core_local_separate": True,
        "formal_powers_materialized": False,
        "formal_atoms_materialized": False,
        "product_coefficient_semantics": "coefficient_enclosure = affine_source_coefficient * W0_sector_coefficient",
        "scalar_application": "deferred_to_consumer_once",
        "interval_sum_performed": False,
        "bounds_are_function_values": False,
        "physical_h5_materialized": False,
        "full_h5": False,
        "full_state": False,
        "state_solved": False,
        "ns_validation": False,
        "R_hi": {
            "symbol": "R_hi",
            "definition": "W - (W0_baseline + W_firstword); conditional on separately binding the saved firstword",
            "preserved": True,
            "known": False,
            "included_in_products": False,
            "added_to_ledger": False,
        },
        "unknown_residual": {
            "symbol": "R_hi",
            "definition": "W - (W0_baseline + W_firstword); conditional on separately binding the saved firstword",
            "preserved": True,
            "included_in_products": False,
            "added_to_ledger": False,
        },
        "conditions": {
            "fixed_x_ordinary_Z": "caller/source expression condition; not proved here",
            "rho_equals_x_squared": "caller premise; no coordinate conversion performed here",
            "geometry_Z_independent": "caller premise; not proved here",
            "same_frame": "caller condition; not proved here",
            "same_scale": "caller condition; not proved here",
            "same_gauge": "caller condition; not proved here",
            "same_state_owner": "caller condition; not proved here",
            "same_state_Z0_Z1_Z2_owner": "caller condition; not proved here",
            "physical_dx_applied": False,
            "endpoint_cancellation": False,
            "error_cancellation": False,
            "axis_chi_cancellation": False,
            "R_hi": "unknown; retained outside products",
        },
        "scope": {
            "source_only_g5_skipped": True,
            "g5_retained_by_firstword": "caller premise; prior pinned receipt documents it",
            "all_C0_B1_products_retained": True,
            "axis_two_chi_products_retained": True,
            "core_and_local_output_separate": True,
            "no_product_deduplication": True,
            "no_endpoint_cancellation": True,
            "no_error_cancellation": True,
            "no_atom_values": True,
            "no_function_value": True,
            "R_hi_unknown": True,
        },
        "resource_limits": {
            "cell_count": len(cells),
            "source_node_count": len(source_nodes),
            "state_node_count": len(state_nodes),
            "skipped_g5_count": len(skipped),
            "product_count": product_index,
            "product_counts": counts,
        },
    }


build_baseline_product = build
build_row5_baseline_product = build
row5_baseline_product = build

