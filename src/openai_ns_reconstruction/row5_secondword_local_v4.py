"""Compact local second-feedback word from the saved first-word ledger.

This leaf consumes the *already emitted* ``row5_firstword_product_v4``
product ledger.  It selects only raw ``C0``/``B1`` feedback products and
attaches them to four fixed ``q`` cells.  The two source rows are the
ordinary ``Z0`` and ``Z1`` rows, respectively.  A cell stores one Green
factor and a term stores only references to that factor and to one product;
the full source-node and product records are retained once in their tables.

The K0 bounds used here are deliberately formal in ``hb``:

``K0 prefix <= hb**2 * S_i**2 / 2`` and
``d_x K0 prefix <= hb * S_i``.

Only the coefficient part of each factor is evaluated.  The output therefore
does not materialize ``hb`` or any source atom, does not sum signed fields,
and does not claim the core second-feedback inlet or a whole-collar ``R_hi``
bound.  ``g5`` and ``+2 chi``/axis products are excluded by source-family
selection; this is a scope filter, not an algebraic cancellation.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "row5_secondword_local_v4"
INPUT_SCHEMA = "row5_firstword_product_v4"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
SOURCE_ROWS = ("h5", "h5_Z1")
OUTPUT_FIELDS = ("F1", "d_x_F1")
CELL_COUNT = 4

_MISSING = object()
_COEFFICIENT_KEYS = (
    "coefficient_enclosure",
    "coefficient_interval",
    "coefficient",
    "directed_interval",
    "enclosure",
    "interval",
)

__all__ = [
    "SCHEMA",
    "INPUT_SCHEMA",
    "SOURCE_POWER_BASIS",
    "SOURCE_ROWS",
    "OUTPUT_FIELDS",
    "CELL_COUNT",
    "build",
    "build_secondword_local",
    "evaluate_term",
    "enclose",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _clone(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if _is_sequence(value):
        return [_clone(item) for item in value]
    return value


def _reject_python_float(value: Any, label: str) -> None:
    if isinstance(value, float):
        raise TypeError(f"{label} may not contain a Python float")
    if isinstance(value, Mapping):
        for key, nested in value.items():
            _reject_python_float(nested, f"{label}.{key}")
    elif _is_sequence(value):
        for index, nested in enumerate(value):
            _reject_python_float(nested, f"{label}[{index}]")


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


def _is_interval(value: Any) -> bool:
    return getattr(value, "_mpi_", None) is not None and hasattr(value, "a") and hasattr(value, "b")


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


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    return c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))


def _interval(c: Any, value: Any, label: str) -> Any:
    """Decode a directed interval without endpoint selection or floats."""

    _reject_python_float(value, label)
    if _is_interval(value):
        result = _make_interval(c, value, value)
    elif isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            if "lower_exact_mpf_tuple" not in value or "upper_exact_mpf_tuple" not in value:
                raise TypeError(f"{label} requires both exact MPF endpoints")
            lower = _exact_mpf_tuple(value["lower_exact_mpf_tuple"], label + ".lower_exact_mpf_tuple")
            upper = _exact_mpf_tuple(value["upper_exact_mpf_tuple"], label + ".upper_exact_mpf_tuple")
            result = c.make_mpf((lower, upper))
        else:
            nested = next((value[key] for key in _COEFFICIENT_KEYS if key in value), _MISSING)
            if nested is not _MISSING and nested is not value:
                return _interval(c, nested, label + ".coefficient")
            lower = next((value[key] for key in ("lower", "lo", "left") if key in value), _MISSING)
            upper = next((value[key] for key in ("upper", "hi", "right") if key in value), _MISSING)
            if lower is _MISSING or upper is _MISSING:
                raise TypeError(f"{label} is not interval-compatible")
            result = _make_interval(c, lower, upper)
    elif _is_sequence(value) and len(value) == 2:
        result = _make_interval(c, value[0], value[1])
    else:
        if isinstance(value, bool):
            raise TypeError(f"{label} may not be boolean")
        try:
            result = c.mpf(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{label} is not interval-compatible") from exc
    if result.a > result.b:
        raise ValueError(f"{label} has reversed directed endpoints")
    return result


def _enclosure_payload(value: Any) -> dict[str, Any]:
    lower, upper = _raw_endpoints(value)
    return {
        "lower_exact_mpf_tuple": _clone(lower),
        "upper_exact_mpf_tuple": _clone(upper),
        "directed": True,
        "function_value": False,
    }


def _nonnegative(c: Any, upper: Any) -> Any:
    zero = c.mpf(0)
    if not _is_interval(upper):
        upper = c.mpf(upper)
    return _make_interval(c, zero, upper)


def _upper(value: Any) -> Any:
    return value.b if _is_interval(value) else value


def _lower(value: Any) -> Any:
    return value.a if _is_interval(value) else value


def _abs_upper(value: Any) -> Any:
    lower = _lower(value)
    upper = _upper(value)
    lower_abs = -lower if lower < 0 else lower
    upper_abs = -upper if upper < 0 else upper
    return max(lower_abs, upper_abs)


def _validate_context(c: Any) -> None:
    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf") or not hasattr(c, "exp"):
        raise TypeError("c must provide mpmath interval mpf, make_mpf, and exp")


def _powers(raw: Mapping[Any, Any], label: str) -> list[int]:
    value = _MISSING
    for key in ("formal_powers", "powers", "source_powers"):
        if key in raw:
            value = raw[key]
            break
    if value is _MISSING:
        raise KeyError(f"{label}.formal_powers is required")
    if isinstance(value, Mapping):
        value = [value.get(name, _MISSING) for name in SOURCE_POWER_BASIS]
    if not _is_sequence(value) or len(value) != 4:
        raise TypeError(f"{label}.formal_powers must contain four integers")
    if any(isinstance(item, bool) or not isinstance(item, int) for item in value):
        raise TypeError(f"{label}.formal_powers must contain integers only")
    return [int(item) for item in value]


def _validate_no_atoms(raw: Mapping[Any, Any], label: str) -> None:
    if {"formal_atoms", "atom_values", "formal_atom_values", "atoms"}.intersection(raw):
        raise ValueError(f"{label} may not materialize formal atoms")


def _cell_intervals(c: Any, q_cells: Any) -> list[Any]:
    if q_cells is None:
        cells = [_make_interval(c, c.mpf(index) / 4, c.mpf(index + 1) / 4) for index in range(CELL_COUNT)]
    else:
        if not _is_sequence(q_cells) or len(q_cells) != CELL_COUNT:
            raise ValueError("exactly four collar q cells are required")
        cells = [_interval(c, value, f"q_cells[{index}]") for index, value in enumerate(q_cells)]
    previous = None
    for index, cell in enumerate(cells):
        if cell.a < 0 or cell.b > 1 or cell.a >= cell.b:
            raise ValueError(f"q_cells[{index}] must be a nonempty subinterval of [0,1]")
        if previous is not None and cell.a != previous.b:
            raise ValueError("four q cells must be contiguous")
        previous = cell
    if cells[0].a != 0 or cells[-1].b != 1:
        raise ValueError("four q cells must cover [0,1]")
    return cells


def _cell_factors(c: Any, sc: Any, cells: Sequence[Any]) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for index, q_interval in enumerate(cells):
        qhi = q_interval.b
        s_interval = sc * q_interval / 4
        # S_i is intentionally independent of the formal hb atom.
        S_i = (sc / 4) * qhi * c.exp(sc / 16)
        K0_coefficient = _nonnegative(c, _upper(S_i * S_i / 2))
        d_x_K0_coefficient = _nonnegative(c, _upper(S_i))
        cell_id = f"q{index}"
        result[cell_id] = {
            "cell_factor_id": cell_id,
            "q_interval": _enclosure_payload(q_interval),
            "s_interval": _enclosure_payload(s_interval),
            "S_i": _enclosure_payload(S_i),
            "geometry": {
                "x_formula": "x=2*exp(hb*s/2)",
                "s_formula": "s=sc*q/4",
                "S_i_formula": "S_i=(sc/4)*qhi*exp(sc/16)",
                "formal_hb_not_materialized": True,
            },
            "factors": {
                "F1": {
                    "factor_id": f"{cell_id}:F1",
                    "formal_hb_power": 2,
                    "coefficient_enclosure": _enclosure_payload(K0_coefficient),
                    "mass_formula": "K0_prefix in [0,hb**2*S_i**2/2]",
                    "kernel": "K0",
                },
                "d_x_F1": {
                    "factor_id": f"{cell_id}:d_x_F1",
                    "formal_hb_power": 1,
                    "coefficient_enclosure": _enclosure_payload(d_x_K0_coefficient),
                    "mass_formula": "d_x_K0_prefix in [0,hb*S_i]",
                    "kernel": "d_x_K0",
                },
            },
        }
    return result


def _error_branch(product: Mapping[str, Any], source_node: Mapping[str, Any], state_node: Mapping[str, Any] | None) -> str:
    flags = (
        product.get("error_only"),
        product.get("uncertainty"),
        product.get("symmetric_error"),
        source_node.get("error_only"),
        source_node.get("uncertainty"),
        source_node.get("symmetric_error"),
    )
    if state_node is not None:
        flags = flags + (
            state_node.get("error_only"),
            state_node.get("uncertainty"),
            state_node.get("symmetric_error"),
        )
    nodes = (product, source_node) + ((state_node,) if state_node is not None else ())
    tagged = any(node.get(key) for node in nodes for key in ("error_classes", "source_error_classes"))
    return "error" if any(value is True for value in flags) or tagged else "finite"


def _source_family(source_node: Mapping[str, Any], product: Mapping[str, Any], label: str) -> str:
    family = source_node.get("source_family", _MISSING)
    if family is _MISSING:
        family = product.get("source_family", _MISSING)
    if not isinstance(family, str) or not family:
        raise KeyError(f"{label} is missing source_family")
    product_family = product.get("source_family", _MISSING)
    if product_family is not _MISSING and product_family != family:
        raise ValueError(f"{label} source_family disagrees with source node")
    return family


def _prepare_products(c: Any, firstword_product: Mapping[str, Any]) -> tuple[dict[str, dict[str, Any]], dict[str, list[str]], list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    if firstword_product.get("schema") != INPUT_SCHEMA:
        raise ValueError("firstword_product must have row5_firstword_product_v4 schema")
    _validate_no_atoms(firstword_product, "firstword_product")
    if tuple(firstword_product.get("source_power_basis", ())) != SOURCE_POWER_BASIS:
        raise ValueError("firstword source power basis differs")
    source_nodes = firstword_product.get("source_nodes", _MISSING)
    if not isinstance(source_nodes, Mapping) or not source_nodes:
        raise KeyError("firstword_product.source_nodes is required")
    rows = firstword_product.get("rows", _MISSING)
    if not isinstance(rows, Mapping):
        raise KeyError("firstword_product.rows is required")
    state_nodes = firstword_product.get("state_nodes", {})
    if not isinstance(state_nodes, Mapping):
        raise TypeError("firstword_product.state_nodes must be a mapping")
    products: dict[str, dict[str, Any]] = {}
    selected_by_row: dict[str, list[str]] = {row: [] for row in SOURCE_ROWS}
    excluded: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    for row in SOURCE_ROWS:
        raw_products = rows.get(row, _MISSING)
        if raw_products is _MISSING or not _is_sequence(raw_products) or not raw_products:
            raise KeyError(f"firstword_product.rows[{row}] is missing or empty")
        for index, raw_product in enumerate(raw_products):
            label = f"firstword_product.rows[{row}][{index}]"
            if not isinstance(raw_product, Mapping):
                raise TypeError(f"{label} must be a product mapping")
            _validate_no_atoms(raw_product, label)
            if raw_product.get("row") != row:
                raise ValueError(f"{label} belongs to a different output row")
            product_id = raw_product.get("product_id", _MISSING)
            if not isinstance(product_id, str) or not product_id:
                raise KeyError(f"{label}.product_id is required")
            if product_id in seen_ids:
                raise ValueError(f"duplicate first-word product_id: {product_id}")
            seen_ids.add(product_id)
            source_node_id = raw_product.get("source_node_id", _MISSING)
            if not isinstance(source_node_id, str) or source_node_id not in source_nodes:
                raise KeyError(f"{label}.source_node_id is missing from source_nodes")
            source_node = source_nodes[source_node_id]
            if not isinstance(source_node, Mapping):
                raise TypeError(f"source_nodes[{source_node_id}] must be a mapping")
            family = _source_family(source_node, raw_product, label)
            if family not in ("C0", "B1", "g5", "axis"):
                raise ValueError(f"{label} has unsupported source_family {family!r}")
            if family not in ("C0", "B1"):
                reason = "source_only_g5" if family == "g5" else "axis_plus_2chi_excluded" if family == "axis" else "non_raw_feedback_family"
                excluded.append(
                    {
                        "product_id": product_id,
                        "row": row,
                        "source_node_id": source_node_id,
                        "source_family": family,
                        "reason": reason,
                        "algebraic_cancellation": False,
                    }
                )
                continue
            if raw_product.get("source_only") is not False:
                raise ValueError(f"{label} raw C0/B1 product must have source_only=False")
            if raw_product.get("scalar_applied_to_coefficient") is not False:
                raise ValueError(f"{label} scalar must remain unapplied in the input ledger")
            coefficient_raw = raw_product.get("coefficient_enclosure", _MISSING)
            if coefficient_raw is _MISSING:
                raise KeyError(f"{label}.coefficient_enclosure is required")
            coefficient = _interval(c, coefficient_raw, label + ".coefficient_enclosure")
            powers = _powers(raw_product, label)
            scalar = raw_product.get("scalar", _MISSING)
            if scalar is _MISSING:
                raise KeyError(f"{label}.scalar is required")
            if isinstance(scalar, bool) or not isinstance(scalar, int):
                raise TypeError(f"{label}.scalar must be an integer")
            state_node_id = raw_product.get("state_node_id")
            if not isinstance(state_node_id, str) or state_node_id not in state_nodes:
                raise KeyError(f"{label}.state_node_id must identify a saved state node")
            state_node = state_nodes[state_node_id]
            if not isinstance(state_node, Mapping):
                raise TypeError(f"state_nodes[{state_node_id}] must be a mapping")
            source_powers = _powers(source_node, label + ".source_node")
            state_powers = _powers(state_node, label + ".state_node")
            if powers != [a + b for a, b in zip(source_powers, state_powers)]:
                raise ValueError(f"{label} formal powers disagree with source/state nodes")
            def error_classes(node):
                result = []
                for key in ("error_classes", "source_error_classes"):
                    for item in node.get(key, []):
                        if item not in result:
                            result.append(_clone(item))
                return result
            record = {
                "product_id": product_id,
                "row": row,
                "source_node_id": source_node_id,
                "source_node_ref": source_node_id,
                "state_node_id": state_node_id,
                "state_node_ref": state_node_id,
                "source_family": family,
                "error_branch": _error_branch(raw_product, source_node, state_node),
                "source_error_classes": error_classes(source_node),
                "state_error_classes": error_classes(state_node),
                "product_error_classes": error_classes(raw_product),
                "scalar": int(scalar),
                "formal_powers": list(powers),
                "coefficient_enclosure": _enclosure_payload(coefficient),
                "original_coefficient_enclosure": _clone(coefficient_raw),
                "source_node_provenance_ref": {
                    "source_node_id": source_node_id,
                    "state_node_id": state_node_id,
                    "node_table_lookup_required": True,
                },
                "scalar_applied_to_coefficient": False,
                "no_error_cancellation": True,
            }
            products[product_id] = record
            selected_by_row[row].append(product_id)
    if not any(selected_by_row.values()):
        raise ValueError("no raw C0/B1 feedback products were selected")
    return products, selected_by_row, excluded, _clone(dict(source_nodes)), _clone(dict(state_nodes))


def evaluate_term(c: Any, product: Mapping[str, Any], cell_factor: Mapping[str, Any]) -> dict[str, Any]:
    """Evaluate one factorized signed product against one cell factor.

    The returned intervals contain only the coefficient part of the formal
    kernel mass.  ``scalar`` is applied exactly once; the ``hb`` power is
    carried in metadata and is never evaluated.
    """

    if not isinstance(product, Mapping) or not isinstance(cell_factor, Mapping):
        raise TypeError("product and cell_factor mappings are required")
    product_id = product.get("product_id", _MISSING)
    factor_id = cell_factor.get("factor_id", _MISSING)
    if not isinstance(product_id, str) or not product_id:
        raise KeyError("product.product_id is required")
    if not isinstance(factor_id, str) or not factor_id:
        raise KeyError("cell_factor.factor_id is required")
    coefficient = product.get("coefficient_enclosure", _MISSING)
    if coefficient is _MISSING:
        raise KeyError("product.coefficient_enclosure is required")
    factor_coefficient = cell_factor.get("coefficient_enclosure", _MISSING)
    if factor_coefficient is _MISSING:
        raise KeyError("cell_factor.coefficient_enclosure is required")
    powers = _powers(product, product_id)
    scalar = product.get("scalar", _MISSING)
    if scalar is _MISSING or isinstance(scalar, bool) or not isinstance(scalar, int):
        raise TypeError("product.scalar must be an integer")
    hb_power = cell_factor.get("formal_hb_power", _MISSING)
    if isinstance(hb_power, bool) or not isinstance(hb_power, int) or hb_power < 0:
        raise ValueError("cell_factor.formal_hb_power must be a nonnegative integer")
    coefficient_interval = _interval(c, coefficient, product_id + ".coefficient")
    factor_interval = _interval(c, factor_coefficient, factor_id + ".coefficient")
    signed = coefficient_interval * int(scalar) * factor_interval
    absolute = _nonnegative(c, _abs_upper(coefficient_interval) * abs(int(scalar)) * _upper(factor_interval))
    output_powers = list(powers)
    output_powers[0] += int(hb_power)
    return {
        "product_id": product_id,
        "factor_id": factor_id,
        "scalar": int(scalar),
        "scalar_applied_once": True,
        "formal_powers": list(powers),
        "output_formal_powers": output_powers,
        "kernel_formal_hb_power": int(hb_power),
        "signed_coefficient_enclosure": _enclosure_payload(_interval(c, signed, product_id + ".signed_term")),
        "absolute_coefficient_upper": _enclosure_payload(absolute),
        "formal_hb_not_materialized": True,
        "no_cancellation": True,
    }


def _product_abs_upper(c: Any, product: Mapping[str, Any]) -> Any:
    coefficient = _interval(c, product["coefficient_enclosure"], product["product_id"] + ".coefficient")
    scalar = product["scalar"]
    return _nonnegative(c, _abs_upper(coefficient) * abs(int(scalar)))


def _product_group_id(order: int, branch: str, powers: Sequence[int]) -> str:
    return "product_group:Z{}:{}:powers:{}".format(order, branch, ",".join(str(power) for power in powers))


def _bound_group_id(cell_id: str, field: str, order: int, branch: str, powers: Sequence[int]) -> str:
    return "bound:{}:{}:Z{}:{}:powers:{}".format(cell_id, field, order, branch, ",".join(str(power) for power in powers))


def build(c: Any, firstword_product: Mapping[str, Any], *, sc: Any, hb_admissibility: Any = _MISSING, q_cells: Any = None) -> dict[str, Any]:
    """Build the factorized four-cell local second-feedback word.

    ``hb_admissibility=True`` is a selected caller certificate.  No concrete
    ``hb`` value is accepted or evaluated.  ``sc`` and optional ``q_cells``
    are directed exact intervals.  The input product ledger is not re-derived
    and no source callback is run.
    """

    _validate_context(c)
    if hb_admissibility is not True:
        raise ValueError("selected hb_admissibility=True certificate is required")
    if not isinstance(firstword_product, Mapping):
        raise TypeError("firstword_product must be a saved product mapping")
    sc_interval = _interval(c, sc, "sc")
    if sc_interval.a <= 0:
        raise ValueError("sc must be strictly positive")
    cells = _cell_intervals(c, q_cells)
    products, selected_by_row, excluded, source_nodes, state_nodes = _prepare_products(c, firstword_product)
    cell_factors = _cell_factors(c, sc_interval, cells)

    product_groups: list[dict[str, Any]] = []
    product_groups_by_row: dict[str, list[dict[str, Any]]] = {row: [] for row in SOURCE_ROWS}
    for order, row in enumerate(SOURCE_ROWS):
        grouped: dict[tuple[Any, ...], dict[str, Any]] = {}
        for product_id in selected_by_row[row]:
            product = products[product_id]
            key = (product["error_branch"], tuple(product["formal_powers"]))
            if key not in grouped:
                branch, powers = key
                grouped[key] = {
                    "product_group_id": _product_group_id(order, branch, powers),
                    "ordinary_Z_order": order,
                    "source_row": row,
                    "error_branch": branch,
                    "formal_powers": list(powers),
                    "member_product_ids": [],
                    "member_count": 0,
                    "source_error_classes": [],
                    "state_error_classes": [],
                    "product_error_classes": [],
                    "no_cancellation": True,
                    "scalar_applied_once": True,
                    "formal_hb_not_materialized": True,
                    "product_abs_coefficient_upper": _nonnegative(c, 0),
                }
            group = grouped[key]
            group["member_product_ids"].append(product_id)
            group["member_count"] += 1
            for key in ("source_error_classes", "state_error_classes", "product_error_classes"):
                for error_class in product[key]:
                    if error_class not in group[key]:
                        group[key].append(_clone(error_class))
            group["product_abs_coefficient_upper"] = group["product_abs_coefficient_upper"] + _product_abs_upper(c, product)
        for group in grouped.values():
            product_groups.append(group)
            product_groups_by_row[row].append(group)

    fields: dict[str, dict[str, dict[str, Any]]] = {field: {} for field in OUTPUT_FIELDS}
    for order, row in enumerate(SOURCE_ROWS):
        group_ids = [group["product_group_id"] for group in product_groups_by_row[row]]
        for field in OUTPUT_FIELDS:
            fields[field][f"Z{order}"] = {
                "source_row": row,
                "product_group_ids": list(group_ids),
                "cell_factor_ids": list(cell_factors),
                "factor_refs": {
                    cell_id: f"cell_factors[{cell_id}].factors[{field}]" for cell_id in cell_factors
                },
                "formal_hb_power": 2 if field == "F1" else 1,
                "factorized": True,
                "term_expansion_performed": False,
                "whole_collar_source_enclosure": True,
            }

    absolute_bounds: list[dict[str, Any]] = []
    for group in product_groups:
        input_powers = group["formal_powers"]
        for field in OUTPUT_FIELDS:
            for cell_id, cell in cell_factors.items():
                factor = cell["factors"][field]
                factor_coefficient = _interval(c, factor["coefficient_enclosure"], factor["factor_id"] + ".coefficient")
                output_powers = list(input_powers)
                output_powers[0] += int(factor["formal_hb_power"])
                bound = _nonnegative(c, _upper(group["product_abs_coefficient_upper"] * factor_coefficient))
                absolute_bounds.append(
                    {
                        "group_id": _bound_group_id(cell_id, field, group["ordinary_Z_order"], group["error_branch"], output_powers),
                        "cell_factor_id": cell_id,
                        "cell_factor_ref": f"cell_factors[{cell_id}].factors[{field}]",
                        "field": field,
                        "ordinary_Z_order": group["ordinary_Z_order"],
                        "source_row": group["source_row"],
                        "error_branch": group["error_branch"],
                        "source_error_classes": _clone(group["source_error_classes"]),
                        "state_error_classes": _clone(group["state_error_classes"]),
                        "product_error_classes": _clone(group["product_error_classes"]),
                        "formal_powers": output_powers,
                        "input_formal_powers": list(input_powers),
                        "member_product_group_id": group["product_group_id"],
                        "member_product_ids_ref": f"product_groups[{group['product_group_id']}].member_product_ids",
                        "member_count": group["member_count"],
                        "product_abs_coefficient_upper": _enclosure_payload(group["product_abs_coefficient_upper"]),
                        "cell_factor_coefficient_upper": _enclosure_payload(factor_coefficient),
                        "absolute_coefficient_upper": _enclosure_payload(bound),
                        "scalar_applied_once": True,
                        "formal_hb_not_materialized": True,
                        "no_cancellation": True,
                        "same_source_reused_across_cells": True,
                        "cell_bound_is_not_independent_source_cell_data": True,
                    }
                )

    return {
        "schema": SCHEMA,
        "component": "local_second_feedback_word_v4",
        "status": "conditional_local_second_feedback_word",
        "input_schema": INPUT_SCHEMA,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "source_rows": {"Z0": "h5", "Z1": "h5_Z1"},
        "output_fields": list(OUTPUT_FIELDS),
        "source_nodes": source_nodes,
        "state_nodes": state_nodes,
        "products": products,
        "product_groups": product_groups,
        "product_group_ids_by_row": {
            row: [group["product_group_id"] for group in product_groups_by_row[row]] for row in SOURCE_ROWS
        },
        "excluded_products": excluded,
        "cell_factors": cell_factors,
        "fields": fields,
        "absolute_coefficient_bounds": absolute_bounds,
        "geometry": {
            "x_formula": "x=2*exp(hb*s/2)",
            "s_formula": "s=sc*q/4",
            "S_i_formula": "S_i=(sc/4)*qhi*exp(sc/16)",
            "q_cell_count": CELL_COUNT,
            "hb_admissibility": "selected certificate hb_admissibility=True",
            "hb_admitted_but_not_evaluated": True,
            "hb_not_materialized": True,
            "sc_interval": _enclosure_payload(sc_interval),
            "formal_hb_not_materialized": True,
        },
        "source_identity_conditions": {
            "same_source_uniform_on_each_full_collar_cell": "caller condition; not proved here",
            "common_source_owner": "caller condition; not proved here",
            "fixed_normalized_x_geometry": "caller condition; not proved here",
            "whole_collar_callback_identity": "caller condition; not proved here",
        },
        "core_feedback_inlet": {
            "status": "unknown",
            "known": False,
            "included": False,
            "reset_to_zero": False,
            "scope": "local feedback source only",
        },
        "scope": {
            "raw_K0_C0_B1_feedback_only": True,
            "source_only_g5_skipped": True,
            "axis_plus_2chi_skipped": True,
            "axis_skip_is_algebraic_cancellation": False,
            "four_q_cells_only": True,
            "whole_collar_source_enclosure": True,
            "source_rows_are_whole_collar_enclosures": True,
            "whole_collar_source_enclosure_verified": False,
            "local_seed_only": True,
            "complete_Rhi_bound": False,
            "allword_tail_added": False,
            "core_history_reset": False,
            "signed_interval_sum_performed": False,
            "absolute_bound_sum_performed": True,
            "factorized_no_product_cell_cartesian_expansion": True,
            "no_error_cancellation": True,
            "no_product_deduplication": True,
            "formal_atoms_materialized": False,
            "formal_hb_materialized": False,
        },
        "unknown_residual": {
            "symbol": "R",
            "included": False,
            "added_to_bounds": False,
        },
    }


build_secondword_local = build
build.evaluate_term = evaluate_term
enclose = build


