"""Target-cell propagation of the five ``+2 chi E51`` kernel errors.

This leaf consumes only two already-decoded output packets.  It does not
import a producer, reconstruct ``W0``, call an older kernel helper, or form a
fixed ``x=2`` exit.  The output is a target-cell enclosure for the five
omitted/error-only paths

``G0 P Rg``, ``D P Rg``, ``R P G0g``, ``R P Dg``, and ``R P Rg``.

``Rg`` is the saved ``j>=2`` Kchi potential-kernel envelope.  Its saved
source-width factor is already part of the envelope and is never multiplied
again.  The only new width used here is the outer R source-cell bound
``1/4``.  G0 uses the deliberately conservative full-triangle bounds
``x_upper**2/8`` and ``x_upper/4``; D uses the saved bare ``T K0`` masses
from each Kchi target/source-cell record.  Every emitted coefficient is a
symmetric interval and carries ``error_only=True``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any


SCHEMA = "original_kchi_compensation_target_tail_v1"
CELL_COUNT = 8
ORDERS = (0, 1)
FIELDS = ("F1", "d_x_F1")
BRANCHES = (
    "G0_P_Rg",
    "D_P_Rg",
    "R_P_G0g",
    "R_P_Dg",
    "R_P_Rg",
)
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")

# These are ordinary derivatives, not Taylor coefficients.  The JSON-friendly
# forms are part of the public report so a caller can audit the routes without
# interpreting the term stream.
P_PRODUCT_RULE_Z0 = ((0, 0, 1, "chi*V_Z0"),)
P_PRODUCT_RULE_Z1 = (
    (1, 0, 1, "chi_Z*V_Z0"),
    (0, 1, 1, "chi*V_Z1"),
)
D_PRODUCT_RULE_Z0 = ((0, 0, 1, "chi*P_Z0"),)
D_PRODUCT_RULE_Z1 = (
    (1, 0, 1, "chi_Z*P_Z0"),
    (0, 1, 1, "chi*P_Z1"),
)
R_KERNEL_RULE_Z0 = ((0, 0, 1, "R_Z0*P_Z0"),)
R_KERNEL_RULE_Z1 = (
    (1, 0, 1, "R_Z1*P_Z0"),
    (0, 1, 1, "R_Z0*P_Z1"),
)


class TargetTailInputError(ValueError):
    """Raised when a caller packet is not a usable decoded output."""


def _mapping(value: Any, label: str) -> Mapping[Any, Any]:
    if not isinstance(value, Mapping):
        raise TargetTailInputError(label + " must be a mapping")
    return value


def _sequence(value: Any, label: str, *, allow_empty: bool = True) -> list[Any]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise TargetTailInputError(label + " must be a sequence")
    result = list(value)
    if not allow_empty and not result:
        raise TargetTailInputError(label + " must not be empty")
    return result


def _copy(value: Any) -> Any:
    return deepcopy(value)


def _lookup(mapping: Mapping[Any, Any], key: Any, label: str) -> Any:
    candidates = [key]
    if isinstance(key, int):
        candidates.extend((str(key), "Z" + str(key)))
    elif isinstance(key, str) and key.isdigit():
        candidates.extend((int(key), "Z" + key))
    elif isinstance(key, str) and key.upper().startswith("Z"):
        text = key[1:]
        if text.isdigit():
            candidates.extend((text, int(text)))
    for candidate in candidates:
        if candidate in mapping:
            return mapping[candidate]
    raise TargetTailInputError(label + " is missing " + str(key))


def _order_value(mapping: Mapping[Any, Any], order: int, label: str) -> Any:
    return _lookup(mapping, order, label + ".order")


def _unwrap_output(value: Any, label: str) -> Mapping[Any, Any]:
    """Accept a decoded output or the small receipt-result wrapper around it."""

    root = _mapping(value, label)
    result = root.get("result")
    if isinstance(result, Mapping):
        nested = result.get("output")
        if isinstance(nested, Mapping):
            selected = dict(nested)
        else:
            selected = dict(result)
        # The low-word runner keeps source_family/source_logs beside its
        # ``result`` payload.  Carry those nonnumeric binding records inward
        # so conditional equality is still checked after selecting
        # ``result.outer_cells``.
        for key in ("source_family", "source_logs"):
            if key not in selected and key in root:
                selected[key] = _copy(root[key])
        return selected
    nested = root.get("output")
    if isinstance(nested, Mapping) and "outer_cells" not in root:
        selected = dict(nested)
        for key in ("source_family", "source_logs"):
            if key not in selected and key in root:
                selected[key] = _copy(root[key])
        return selected
    return root


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _interval(c: Any, value: Any, label: str = "interval") -> Any:
    """Decode directed endpoint records without converting through floats."""

    if _is_interval(value):
        return value
    if isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            if value.get("directed") is not True:
                raise TargetTailInputError(label + " must retain directed=True")
            lower = tuple(value["lower_exact_mpf_tuple"])
            upper = tuple(value["upper_exact_mpf_tuple"])
            if len(lower) != 4 or len(upper) != 4:
                raise TargetTailInputError(label + " has malformed endpoint tuples")
            return c.make_mpf((lower, upper))
        if "lower" in value and "upper" in value:
            return c.mpf([value["lower"], value["upper"]])
        if "lo" in value and "hi" in value:
            return c.mpf([value["lo"], value["hi"]])
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return c.mpf([value[0], value[1]])
    return c.mpf(value)


def _point(c: Any, numerator: int, denominator: int) -> Any:
    value = c.mpf(numerator) / denominator
    raw = value._mpi_
    return c.make_mpf((raw[0], raw[1]))


def _upper_endpoint(c: Any, value: Any) -> Any:
    interval = _interval(c, value)
    raw = interval._mpi_[1]
    return c.make_mpf((raw, raw))


def _lower_endpoint(c: Any, value: Any) -> Any:
    interval = _interval(c, value)
    raw = interval._mpi_[0]
    return c.make_mpf((raw, raw))


def _abs_upper(c: Any, value: Any) -> Any:
    return _upper_endpoint(c, abs(_interval(c, value)))


def _symmetric(c: Any, value: Any) -> Any:
    radius = _abs_upper(c, value)
    negative = -radius
    return c.make_mpf((negative._mpi_[0], radius._mpi_[1]))


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    seen: set[Any] = set()
    for value in values:
        try:
            key = (type(value).__name__, value)
            hash(key)
        except TypeError:
            key = (type(value).__name__, repr(value))
        if key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _items(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    return [value]


def _nested_origin(raw: Any) -> Mapping[Any, Any]:
    if not isinstance(raw, Mapping):
        return {}
    origin = raw.get("origin")
    return origin if isinstance(origin, Mapping) else {}


def _ids(raw: Any, default_id: Any) -> list[Any]:
    if not isinstance(raw, Mapping):
        return [default_id]
    origin = _nested_origin(raw)
    values = _items(raw.get("source_ids", raw.get("source_id")))
    values.extend(_items(origin.get("source_ids", origin.get("source_id"))))
    values = [value for value in values if value is not None and value != ""]
    return _unique(values or [default_id])


def _classes(raw: Any) -> list[Any]:
    if not isinstance(raw, Mapping):
        return []
    origin = _nested_origin(raw)
    values = _items(raw.get("error_classes", raw.get("error_class")))
    values.extend(_items(origin.get("error_classes", origin.get("error_class"))))
    return _unique(values)


def _powers(raw: Any, label: str) -> tuple[int, int, int, int]:
    if not isinstance(raw, Mapping):
        return (0, 0, 0, 0)
    value = raw.get("powers")
    if value is None:
        value = (
            raw.get("hb_power", 0),
            raw.get("F0_power", raw.get("F0_anchor_power", 0)),
            raw.get("Pstar_squared_power", raw.get("Pstar_power", 0)),
            raw.get("Lambda_power", raw.get("lambda_power", 0)),
        )
    if (
        isinstance(value, (str, bytes, bytearray))
        or not isinstance(value, Sequence)
        or len(value) != 4
        or any(type(item) is not int for item in value)
    ):
        raise TargetTailInputError(label + " must contain four integer powers")
    return tuple(int(item) for item in value)


def _normal_term(c: Any, raw: Any, *, label: str, default_id: Any) -> dict[str, Any]:
    if isinstance(raw, Mapping):
        if "coefficient" in raw:
            coefficient = raw["coefficient"]
        elif "lower_exact_mpf_tuple" in raw and "upper_exact_mpf_tuple" in raw:
            # The kchi runner encodes chi[0]/chi[1] as bare directed
            # endpoint records rather than one-term dictionaries.  Treat the
            # record as the coefficient while retaining a formal zero power
            # vector and a synthetic source id; it is still an input interval,
            # never a point value.
            coefficient = raw
        else:
            raise TargetTailInputError(label + " is missing coefficient")
    else:
        coefficient = raw
    result = {
        "coefficient": _interval(c, coefficient, label + ".coefficient"),
        "powers": _powers(raw, label),
        "source_ids": _ids(raw, default_id),
        "error_classes": _classes(raw),
        "uncertainty": bool(
            isinstance(raw, Mapping)
            and (raw.get("uncertainty") is True or raw.get("is_uncertainty") is True)
        ),
        "source_id": default_id,
    }
    if isinstance(raw, Mapping):
        result["source_id"] = raw.get("source_id", result["source_ids"][0])
        result["input_origin"] = _copy(raw.get("origin"))
        result["input_error_only"] = bool(raw.get("error_only") is True)
    return result


def _terms(c: Any, raw: Any, *, label: str, default_id: Any) -> list[dict[str, Any]]:
    if raw is None:
        return []
    if isinstance(raw, Mapping) and ("coefficient" in raw or "powers" in raw):
        values = [raw]
    elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
        values = list(raw)
    else:
        values = [raw]
    return [
        _normal_term(c, value, label=f"{label}[{index}]", default_id=f"{default_id}:{index}")
        for index, value in enumerate(values)
    ]


def _cell_list(value: Any, label: str) -> list[Mapping[Any, Any]]:
    if isinstance(value, Mapping):
        cells = [_lookup(value, index, label) for index in range(CELL_COUNT)]
    else:
        cells = _sequence(value, label, allow_empty=False)
    if len(cells) != CELL_COUNT:
        raise TargetTailInputError(label + " must contain eight canonical cells")
    return [_mapping(cell, f"{label}[{index}]") for index, cell in enumerate(cells)]


def _family_and_logs(low: Mapping[Any, Any], kchi: Mapping[Any, Any]) -> dict[str, Any]:
    low_family = low.get("source_family")
    kchi_family = kchi.get("source_family")
    low_logs = low.get("source_logs")
    kchi_logs = kchi.get("source_logs")
    family_present = low_family is not None and kchi_family is not None
    logs_present = low_logs is not None and kchi_logs is not None
    if family_present and low_family != kchi_family:
        raise TargetTailInputError("source_family differs between low and Kchi outputs")
    if logs_present and low_logs != kchi_logs:
        raise TargetTailInputError("source_logs differs between low and Kchi outputs")
    return {
        "source_family": _copy(low_family if low_family is not None else kchi_family),
        "source_logs": _copy(low_logs if low_logs is not None else kchi_logs),
        "family_equal_conditionally": bool(family_present),
        "logs_equal_conditionally": bool(logs_present),
        "live_source_owner_proved": False,
    }


def _rows_from_w0(c: Any, low: Mapping[Any, Any]) -> list[dict[int, list[dict[str, Any]]]]:
    cells = _cell_list(low.get("outer_cells"), "low_output.outer_cells")
    rows: list[dict[int, list[dict[str, Any]]]] = []
    for index, cell in enumerate(cells):
        packet = _mapping(cell.get("W0_ordinary_Z"), f"W0 cell {index}.W0_ordinary_Z")
        row: dict[int, list[dict[str, Any]]] = {}
        for order in ORDERS:
            columns = _mapping(_order_value(packet, order, f"W0 cell {index}"), f"W0 cell {index}.Z{order}")
            values = _lookup(columns, 1, f"W0 cell {index}.Z{order} column")
            row[order] = _terms(c, values, label=f"W0[{index}] Z{order}", default_id=f"W0:{index}:Z{order}")
        rows.append(row)
    return rows


def _rows_from_delta(c: Any, kchi: Mapping[Any, Any]) -> list[dict[int, list[dict[str, Any]]]]:
    first = _mapping(kchi.get("first_correction"), "kchi_output.first_correction")
    cells = _cell_list(first.get("target_cells"), "kchi_output.first_correction.target_cells")
    rows: list[dict[int, list[dict[str, Any]]]] = []
    for index, cell in enumerate(cells):
        field = _mapping(cell.get("F1"), f"Delta cell {index}.F1")
        row = {
            order: _terms(
                c,
                _order_value(field, order, f"Delta cell {index}.F1"),
                label=f"Delta[{index}] Z{order}",
                default_id=f"Delta:{index}:Z{order}",
            )
            for order in ORDERS
        }
        rows.append(row)
    return rows


def _rows_from_tail(c: Any, kchi: Mapping[Any, Any]) -> list[dict[int, list[dict[str, Any]]]]:
    tail_root = _mapping(kchi.get("W0_separate_symmetric_kernel_tail"), "Kchi W0 tail")
    field = _mapping(tail_root.get("F1"), "Kchi W0 tail.F1")
    # The saved Rg envelope is a uniform whole-core bound.  Reusing the same
    # rows at every source cell is a structural envelope conversion, not eight
    # independent source measurements.
    common = {
        order: _terms(
            c,
            _order_value(field, order, "Kchi W0 tail.F1"),
            label=f"Rg Z{order}",
            default_id=f"Rg:Z{order}",
        )
        for order in ORDERS
    }
    # P=2 chi E51 consumes column 1 (F1), including when the outer
    # inverse produces d_x_F1.  The latter uses its own *outer kernel*
    # bound; substituting the input d_x_F1 tail would change E51 to E55.
    for order in ORDERS:
        for term in common[order]:
            if term.get("uncertainty") is not True or "Kchi_j_ge_2_kernel_tail" not in term.get("error_classes", []):
                raise TargetTailInputError("Rg must retain its saved kernel-error class and uncertainty")
            lower, upper = term["coefficient"]._mpi_
            if lower[1] == 0 and upper[1] == 0:
                continue
            if lower[0] != 1 or upper[0] != 0 or lower[1:] != upper[1:]:
                raise TargetTailInputError("Rg coefficient must be a symmetric error enclosure")
    return [
        {
            order: [dict(term, source_ids=list(term["source_ids"])) for term in common[order]]
            for order in ORDERS
        }
        for _ in range(CELL_COUNT)
    ]


def _chi_rows(c: Any, kchi: Mapping[Any, Any]) -> dict[int, list[dict[str, Any]]]:
    root = _mapping(kchi.get("chi"), "Kchi chi")
    raw = root.get("chi", root)
    if isinstance(raw, Mapping):
        values = {order: _order_value(raw, order, "Kchi chi") for order in ORDERS}
    elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
        if len(raw) < 2:
            raise TargetTailInputError("Kchi chi must contain Z0 and Z1")
        values = {order: raw[order] for order in ORDERS}
    else:
        raise TargetTailInputError("Kchi chi must be an order mapping or sequence")
    return {
        order: _terms(c, values[order], label=f"chi Z{order}", default_id=f"chi:Z{order}")
        for order in ORDERS
    }


def _tail_bounds(c: Any, kchi: Mapping[Any, Any]) -> dict[str, dict[int, Any]]:
    root = _mapping(kchi.get("kernel_tail_bounds"), "Kchi kernel_tail_bounds")
    result: dict[str, dict[int, Any]] = {}
    for field in FIELDS:
        values = root.get(field)
        if isinstance(values, Mapping):
            result[field] = {
                order: _interval(c, _order_value(values, order, f"tail {field}"), f"tail {field} Z{order}")
                for order in ORDERS
            }
        else:
            values_list = _sequence(values, f"tail {field}", allow_empty=False)
            if len(values_list) < 2:
                raise TargetTailInputError(f"tail {field} must contain Z0/Z1")
            result[field] = {
                order: _interval(c, values_list[order], f"tail {field} Z{order}")
                for order in ORDERS
            }
    return result


def _saved_d_mass(c: Any, kchi_cells: list[Mapping[Any, Any]], target: int,
                  source: int, field: str) -> Any:
    cell = kchi_cells[target]
    by_source = _mapping(
        cell.get("W0_correction_by_source_cell"),
        f"Kchi target {target}.W0_correction_by_source_cell",
    )
    source_packet = _mapping(
        _lookup(by_source, source, f"Kchi target {target} source"),
        f"Kchi target {target} source {source}",
    )
    masses = _mapping(source_packet.get("kernel_masses"), f"Kchi target {target} source {source}.kernel_masses")
    return _interval(c, _lookup(masses, field, f"saved D mass {field}"),
                     f"saved D mass target {target} source {source} {field}")


def _add_powers(left: Sequence[int], right: Sequence[int]) -> tuple[int, int, int, int]:
    return tuple(a + b for a, b in zip(left, right))  # type: ignore[return-value]


def _term_ids(*terms: Mapping[str, Any]) -> list[Any]:
    return _unique([item for term in terms for item in term.get("source_ids", [])])


def _term_classes(*terms: Mapping[str, Any]) -> list[Any]:
    return _unique([item for term in terms for item in term.get("error_classes", [])])


def _p_term(
    c: Any,
    chi_term: Mapping[str, Any],
    value_term: Mapping[str, Any],
    *,
    target: int,
    source: int,
    output_order: int,
    chi_order: int,
    value_order: int,
    value_kind: str,
    symmetric: bool,
) -> dict[str, Any]:
    coefficient = 2 * chi_term["coefficient"] * value_term["coefficient"]
    if symmetric:
        coefficient = _symmetric(c, coefficient)
    ids = _term_ids(chi_term, value_term)
    classes = _term_classes(chi_term, value_term)
    route = "chi_Z*V_Z0" if chi_order == 1 else "chi*V_Z1" if value_order == 1 else "chi*V_Z0"
    powers = _add_powers(chi_term["powers"], value_term["powers"])
    return {
        "coefficient": coefficient,
        "powers": powers,
        "formal_source_powers": list(powers),
        "source_ids": ids,
        "parent_source_ids": list(ids),
        "parent_source_ids_by_factor": [
            list(chi_term.get("source_ids", [])),
            list(value_term.get("source_ids", [])),
        ],
        "parent_powers_by_factor": [
            list(chi_term["powers"]),
            list(value_term["powers"]),
        ],
        "error_classes": classes,
        "uncertainty": True,
        "error_only": symmetric,
        "symmetric_error": symmetric,
        "signed_field": not symmetric,
        "target_cell": target,
        "source_cell": source,
        "ordinary_Z_order": output_order,
        "chi_order": chi_order,
        "value_order": value_order,
        "P_route": route,
        "P_formula": "+2 chi E51",
        "value_kind": value_kind,
        "parent_input_may_be_negative": bool(value_term["coefficient"].a < 0),
        "inner_Rg_width_readded": False,
    }


def _p_rows(
    c: Any,
    chi: Mapping[int, Sequence[Mapping[str, Any]]],
    values: Mapping[int, Sequence[Mapping[str, Any]]],
    *,
    target: int,
    source: int,
    value_kind: str,
    symmetric: bool,
) -> dict[int, list[dict[str, Any]]]:
    result: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    for chi_term in chi[0]:
        for value_term in values[0]:
            result[0].append(_p_term(
                c, chi_term, value_term, target=target, source=source,
                output_order=0, chi_order=0, value_order=0,
                value_kind=value_kind, symmetric=symmetric,
            ))
    for chi_term in chi[1]:
        for value_term in values[0]:
            result[1].append(_p_term(
                c, chi_term, value_term, target=target, source=source,
                output_order=1, chi_order=1, value_order=0,
                value_kind=value_kind, symmetric=symmetric,
            ))
    for chi_term in chi[0]:
        for value_term in values[1]:
            result[1].append(_p_term(
                c, chi_term, value_term, target=target, source=source,
                output_order=1, chi_order=0, value_order=1,
                value_kind=value_kind, symmetric=symmetric,
            ))
    return result


def _error_term(
    c: Any,
    coefficient: Any,
    *,
    parent_terms: Sequence[Mapping[str, Any]],
    target: int,
    source: int,
    order: int,
    field: str,
    branch: str,
    route: str,
    metadata: Mapping[str, Any],
) -> dict[str, Any]:
    ids = _term_ids(*parent_terms)
    classes = _term_classes(*parent_terms)
    powers: tuple[int, int, int, int] = (0, 0, 0, 0)
    for term in parent_terms:
        powers = _add_powers(powers, term["powers"])
    classes = _unique(classes + [branch, "potential_compensation_target_kernel_error"])
    result = {
        "coefficient": _symmetric(c, coefficient),
        "powers": powers,
        "formal_source_powers": list(powers),
        "source_ids": ids,
        "parent_source_ids": list(ids),
        "parent_source_ids_by_factor": [
            list(term.get("source_ids", [])) for term in parent_terms
        ],
        "parent_powers_by_factor": [
            list(term.get("powers", ())) for term in parent_terms
        ],
        "parent_factor_provenance": [
            {
                "source_ids_by_factor": _copy(term.get("parent_source_ids_by_factor")),
                "powers_by_factor": _copy(term.get("parent_powers_by_factor")),
            }
            for term in parent_terms
            if term.get("parent_source_ids_by_factor") is not None
        ],
        "error_classes": classes,
        "uncertainty": True,
        "error_only": True,
        "symmetric_error": True,
        "signed_field": False,
        "target_cell": target,
        "source_cell": source,
        "ordinary_Z_order": order,
        "field": field,
        "branch": branch,
        "route": route,
        "symbolic_cancellation_claimed": False,
        "parent_term_count": len(parent_terms),
        "parent_routes": [term.get("P_route") for term in parent_terms
                           if term.get("P_route") is not None],
        "parent_input_may_be_negative": any(
            bool(term.get("parent_input_may_be_negative", False)) for term in parent_terms
        ),
    }
    result.update(dict(metadata))
    return result


def _g0_mass_upper(c: Any, target: int, field: str) -> Any:
    x_upper = _point(c, target + 1, 4)
    return x_upper * x_upper / 8 if field == "F1" else x_upper / 4


def _append_g0_p_rg(
    c: Any,
    output: list[dict[str, Any]],
    p_rows: Mapping[int, Sequence[Mapping[str, Any]]],
    *,
    target: int,
    source: int,
    field: str,
    order: int,
) -> None:
    mass = _g0_mass_upper(c, target, field)
    for p_term in p_rows[order]:
        output.append(_error_term(
            c,
            mass * p_term["coefficient"],
            parent_terms=(p_term,),
            target=target,
            source=source,
            order=order,
            field=field,
            branch="G0_P_Rg",
            route="G0_mass_after_P_Rg",
            metadata={
                "operator": "G0 P Rg",
                "kernel_mass_upper_bound": mass,
                "kernel_mass_field": field,
                "kernel_mass_source": "full_triangle_upper_x2_over_8_or_x_over_4",
                "conservative_G0_mass_upper": True,
                "saved_TK0_mass_used": False,
                "outer_source_width_bound": None,
                "inner_Rg_width_readded": False,
            },
        ))


def _append_d_p_rg(
    c: Any,
    output: list[dict[str, Any]],
    p_rows: Mapping[int, Sequence[Mapping[str, Any]]],
    chi: Mapping[int, Sequence[Mapping[str, Any]]],
    *,
    target: int,
    source: int,
    field: str,
    order: int,
    mass: Any,
) -> None:
    routes = ((0, 0, "chi*P_Z0"),) if order == 0 else (
        (1, 0, "chi_Z*P_Z0"),
        (0, 1, "chi*P_Z1"),
    )
    for chi_order, p_order, route in routes:
        for chi_term in chi[chi_order]:
            for p_term in p_rows[p_order]:
                output.append(_error_term(
                    c,
                    -2 * chi_term["coefficient"] * mass * p_term["coefficient"],
                    parent_terms=(chi_term, p_term),
                    target=target,
                    source=source,
                    order=order,
                    field=field,
                    branch="D_P_Rg",
                    route=route,
                    metadata={
                        "operator": "D P Rg",
                        "kernel_mass": mass,
                        "kernel_mass_field": field,
                        "kernel_mass_source": "saved_kchi_TK0_mass",
                        "saved_TK0_mass_used": True,
                        "saved_TK0_mass_reintegrated": False,
                        "chi_levels": ["P", "D"],
                        "outer_chi_order": chi_order,
                        "P_order": p_order,
                        "outer_source_width_bound": None,
                        "inner_Rg_width_readded": False,
                    },
                ))


def _append_r_p(
    c: Any,
    output: list[dict[str, Any]],
    p_rows: Mapping[int, Sequence[Mapping[str, Any]]],
    r_bounds: Mapping[str, Mapping[int, Any]],
    *,
    target: int,
    source: int,
    field: str,
    order: int,
    branch: str,
    value_kind: str,
) -> None:
    routes = ((0, 0, "R_Z0*P_Z0"),) if order == 0 else (
        (1, 0, "R_Z1*P_Z0"),
        (0, 1, "R_Z0*P_Z1"),
    )
    width = _point(c, 1, 4)
    for r_order, p_order, route in routes:
        bound = r_bounds[field][r_order]
        for p_term in p_rows[p_order]:
            output.append(_error_term(
                c,
                width * bound * p_term["coefficient"],
                parent_terms=(p_term,),
                target=target,
                source=source,
                order=order,
                field=field,
                branch=branch,
                route=route,
                metadata={
                    "operator": branch.replace("_", " "),
                    "value_kind": value_kind,
                    "kernel_bound": bound,
                    "kernel_bound_field": field,
                    "kernel_bound_order": r_order,
                    "kernel_bound_source": "saved_kchi_kernel_tail_bounds",
                    "outer_source_width_bound": width,
                    "outer_source_width_bound_max": "1/4",
                    "outer_source_width_applied_once": True,
                    "inner_Rg_width_readded": False,
                    "saved_TK0_mass_used": False,
                },
            ))


def _empty_fields() -> dict[str, dict[str, list[dict[str, Any]]]]:
    return {
        "Z0": {field: [] for field in FIELDS},
        "Z1": {field: [] for field in FIELDS},
    }


def build(c: Any, low_output: Any, kchi_output: Any) -> dict[str, Any]:
    """Build target-cell symmetric errors for all five P-tail paths.

    ``low_output`` and ``kchi_output`` may be direct decoded outputs or the
    small ``{"result": {"output": ...}}`` wrapper.  The returned packet has
    no fixed-exit branch and no aggregate/cancellation operation.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must provide directed mpf and make_mpf")
    low = _unwrap_output(low_output, "low_output")
    kchi = _unwrap_output(kchi_output, "kchi_output")
    binding = _family_and_logs(low, kchi)
    w0_rows = _rows_from_w0(c, low)
    delta_rows = _rows_from_delta(c, kchi)
    rg_rows = _rows_from_tail(c, kchi)
    chi = _chi_rows(c, kchi)
    r_bounds = _tail_bounds(c, kchi)
    first = _mapping(kchi.get("first_correction"), "kchi_output.first_correction")
    kchi_cells = _cell_list(first.get("target_cells"), "kchi_output.first_correction.target_cells")

    # P rows are source-cell local.  The three packets are independent parent
    # streams; no row is summed with another row.
    p_w0: list[dict[int, list[dict[str, Any]]]] = []
    p_delta: list[dict[int, list[dict[str, Any]]]] = []
    p_rg: list[dict[int, list[dict[str, Any]]]] = []
    for source in range(CELL_COUNT):
        p_w0.append(_p_rows(c, chi, w0_rows[source], target=source,
                            source=source, value_kind="G0g", symmetric=False))
        p_delta.append(_p_rows(c, chi, delta_rows[source], target=source,
                               source=source, value_kind="Dg", symmetric=False))
        p_rg.append(_p_rows(c, chi, rg_rows[source], target=source,
                            source=source, value_kind="Rg", symmetric=True))

    branches: dict[str, Any] = {}
    for branch in BRANCHES:
        branches[branch] = {
            "branch": branch,
            "error_only": True,
            "symmetric_error": True,
            "target_cell_profile": True,
            "fixed_x2_exit_computed": False,
            "complete_feedback": False,
            "complete_W1": False,
            "appended_to_ledger": False,
            "feedback_word_tail_included": False,
            "collar_included": False,
            "pressure_included": False,
            "target_cells": [],
        }

    for target in range(CELL_COUNT):
        x_lower = _point(c, target, 4)
        x_upper = _point(c, target + 1, 4)
        target_meta = {
            "target_index": target,
            "target_radius_lower": x_lower,
            "target_radius_upper_bound": x_upper,
            "target_radius_segment": [x_lower, x_upper],
            "source_indices": list(range(target + 1)),
            "source_j_le_target_i": True,
            "not_fixed_exit": True,
        }
        for branch in BRANCHES:
            branches[branch]["target_cells"].append(dict(target_meta, by_source_cell={}))

        for source in range(target + 1):
            g0_fields = _empty_fields()
            d_fields = _empty_fields()
            r_w0_fields = _empty_fields()
            r_delta_fields = _empty_fields()
            r_rg_fields = _empty_fields()
            for order in ORDERS:
                for field in FIELDS:
                    _append_g0_p_rg(
                        c, g0_fields[f"Z{order}"][field], p_rg[source],
                        target=target, source=source, field=field, order=order,
                    )
                    _append_d_p_rg(
                        c, d_fields[f"Z{order}"][field], p_rg[source], chi,
                        target=target, source=source, field=field, order=order,
                        mass=_saved_d_mass(c, kchi_cells, target, source, field),
                    )
                    _append_r_p(
                        c, r_w0_fields[f"Z{order}"][field], p_w0[source], r_bounds,
                        target=target, source=source, field=field, order=order,
                        branch="R_P_G0g", value_kind="G0g",
                    )
                    _append_r_p(
                        c, r_delta_fields[f"Z{order}"][field], p_delta[source], r_bounds,
                        target=target, source=source, field=field, order=order,
                        branch="R_P_Dg", value_kind="Dg",
                    )
                    _append_r_p(
                        c, r_rg_fields[f"Z{order}"][field], p_rg[source], r_bounds,
                        target=target, source=source, field=field, order=order,
                        branch="R_P_Rg", value_kind="Rg",
                    )

            packet_by_branch = {
                "G0_P_Rg": g0_fields,
                "D_P_Rg": d_fields,
                "R_P_G0g": r_w0_fields,
                "R_P_Dg": r_delta_fields,
                "R_P_Rg": r_rg_fields,
            }
            for branch, fields in packet_by_branch.items():
                branches[branch]["target_cells"][target]["by_source_cell"][str(source)] = {
                    "target_index": target,
                    "source_index": source,
                    "target_radius_upper_bound": x_upper,
                    "fields": fields,
                    "source_width_bound": _point(c, 1, 4) if branch.startswith("R_") else None,
                    "inner_Rg_width_readded": False,
                    "no_aggregation": True,
                    "no_cancellation": True,
                }

    return {
        "schema": SCHEMA,
        "source_family": binding["source_family"],
        "source_logs": binding["source_logs"],
        "source_binding": binding,
        "branches": branches,
        "target_cell_count": CELL_COUNT,
        "target_coordinate": "physical x",
        "target_radius_segments": "[i/4,(i+1)/4], i=0,...,7",
        "target_radius_upper_bound_formula": "x_upper=(i+1)/4",
        "ordinary_Z_orders": [0, 1],
        "fields": list(FIELDS),
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "formal_powers_not_expanded": True,
        "directed_endpoint_policy": "input endpoint tuples are decoded with c.make_mpf; output radii use upper abs and symmetric [-r,r]",
        "coefficients_are_field_values": False,
        "P_formula": "P=+2 chi E51; P_Z0=2 chi V; P_Z1=2(chi_Z V+chi V_Z)",
        "D_formula": "D=-2 chi T K0 applied after P",
        "R_formula": "R tail majorant; R_Z1 P0 + R_Z0 P1",
        "ordinary_Z0_product_rule": {
            "P": P_PRODUCT_RULE_Z0,
            "D": D_PRODUCT_RULE_Z0,
            "R": R_KERNEL_RULE_Z0,
        },
        "ordinary_Z1_product_rule": {
            "P": P_PRODUCT_RULE_Z1,
            "D": D_PRODUCT_RULE_Z1,
            "R": R_KERNEL_RULE_Z1,
        },
        "G0_mass_policy": "conservative full-triangle upper bounds x_upper^2/8 and x_upper/4",
        "D_mass_policy": "saved kchi first-correction target/source kernel_masses only",
        "R_mass_policy": "saved kernel_tail_bounds multiplied by outer source width at most 1/4 once",
        "inner_Rg_source_width_bound": 2,
        "inner_Rg_width_readded": False,
        "all_coefficients_symmetric_error_only": True,
        "fixed_x2_exit_computed": False,
        "feedback_word_tail_included": False,
        "collar_included": False,
        "pressure_included": False,
        "appended_to_ledger": False,
        "complete_W1": False,
        "complete_m15": False,
        "terminal_closure_proved": False,
        "source_owner_identity_proved": False,
        "saved_inputs_reintegrated": False,
        "old_kernel_helper_called": False,
    }


target_tail = build
build_target_tail = build


__all__ = [
    "SCHEMA", "CELL_COUNT", "ORDERS", "FIELDS", "BRANCHES",
    "SOURCE_POWER_BASIS", "P_PRODUCT_RULE_Z0", "P_PRODUCT_RULE_Z1",
    "D_PRODUCT_RULE_Z0", "D_PRODUCT_RULE_Z1", "R_KERNEL_RULE_Z0",
    "R_KERNEL_RULE_Z1", "TargetTailInputError", "build", "target_tail",
    "build_target_tail",
]
