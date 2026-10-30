"""Finite core ``W0`` physical ``m15`` source sectors.

The caller supplies the eight saved core cells and the already enclosed
leading physical fields.  This module only performs the finite cell
products and the exact ``rho = 4*t**2`` cell weights; it does not replay a
producer, read a receipt, construct ``W1``, or add an omitted-word tail.

The four formal powers are ordered as ``(hb, F0_anchor, Pstar_squared,
Lambda)``.  A physical axial product contributes one external
``Lambda**-1`` and a physical swirl product contributes two.  ``g4`` is
already the saved raw source ``2*epsilon*sqrt(rho)*N1p``; the pressure
factor adds only the one remaining ``Lambda**-1`` and the exact mass
``-8/3``.

All input coefficients are directed intervals.  Products and sums are
performed in the caller's interval context.  Sectors are grouped only by
formal powers, the uncertainty flag, and the error-class set.  Grouping is
outward interval addition and never asserts covariance, point values, or a
same-source cancellation.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
OUTER_CELL_COUNT = 8
FIELDS = ("F1", "Uz1")

# The core source certificate is deliberately small.  It authenticates the
# source identity, the Z-independent physical epsilon geometry, and the fact
# that all eight W0 outer cells came from the caller's saved construction.
REQUIRED_CERTIFICATES = (
    "g4_equals_two_epsilon_x_N1p",
    "geometry_epsilon_Z_independent",
    "actual_W0_outer_cells",
)

__all__ = [
    "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "OUTER_CELL_COUNT",
    "FIELDS",
    "REQUIRED_CERTIFICATES",
    "build",
]


def _is_interval(value: Any) -> bool:
    """Return whether *value* is a directed interval from the caller."""

    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _raw_endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    return value


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    """Construct an interval from directed endpoints without midpointing."""

    return c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))


def _same_interval(left: Any, right: Any) -> bool:
    """Compare raw directed endpoints, not interval proxy objects."""

    try:
        return left._mpi_ == right._mpi_
    except AttributeError:
        return left == right


def _as_items(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    return [value]


def _freeze(value: Any) -> Any:
    """Make arbitrary provenance values usable as deterministic keys."""

    try:
        hash(value)
    except TypeError:
        return (type(value).__name__, repr(value))
    return (type(value).__name__, value)


def _unique(values: Iterable[Any]) -> list[Any]:
    result: list[Any] = []
    seen: set[Any] = set()
    for value in values:
        key = _freeze(value)
        if key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _origin_info(value: Any, *, _seen: set[int] | None = None) -> dict[str, Any]:
    """Collect source IDs, uncertainty, and error classes from provenance."""

    if _seen is None:
        _seen = set()
    marker = id(value)
    if marker in _seen:
        return {"source_ids": [], "error_classes": [], "uncertainty": False}

    if isinstance(value, Mapping):
        _seen.add(marker)
        raw_ids = value.get("source_ids", value.get("sources"))
        if raw_ids is None:
            raw_ids = value.get("source_id", value.get("source"))
        ids = list(_as_items(raw_ids))
        raw_classes = value.get(
            "error_classes",
            value.get("uncertainty_class", value.get("error_class")),
        )
        classes = list(_as_items(raw_classes))
        operation = str(value.get("operation", "")).lower()
        error_kind = str(value.get("error_kind", value.get("kind", ""))).lower()
        uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False)))
        uncertain = uncertain or any(
            token in operation or token in error_kind
            for token in ("cauchy", "symmetric", "uncertainty")
        )
        for key in ("origin", "metadata", "parents"):
            child = value.get(key)
            if child is None or child is value:
                continue
            if isinstance(child, Sequence) and not isinstance(child, (str, bytes, bytearray)):
                nested = [_origin_info(item, _seen=_seen) for item in child]
            else:
                nested = [_origin_info(child, _seen=_seen)]
            for info in nested:
                ids.extend(info["source_ids"])
                classes.extend(info["error_classes"])
                uncertain = uncertain or bool(info["uncertainty"])
        return {
            "source_ids": _unique(ids),
            "error_classes": _unique(classes),
            "uncertainty": bool(uncertain),
        }

    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        ids: list[Any] = []
        classes: list[Any] = []
        uncertain = False
        for item in value:
            info = _origin_info(item, _seen=_seen)
            ids.extend(info["source_ids"])
            classes.extend(info["error_classes"])
            uncertain = uncertain or bool(info["uncertainty"])
        return {
            "source_ids": _unique(ids),
            "error_classes": _unique(classes),
            "uncertainty": uncertain,
        }

    return {"source_ids": [], "error_classes": [], "uncertainty": False}


def _origin_ref(value: Any) -> dict[str, Any]:
    info = _origin_info(value)
    return {
        "source_ids": list(info["source_ids"]),
        "error_classes": list(info["error_classes"]),
        "uncertainty": bool(info["uncertainty"]),
    }


def _new_origin(operation: str, *inputs: Any, **metadata: Any) -> dict[str, Any]:
    """Union provenance for a product while retaining both source paths."""

    ids: list[Any] = []
    classes: list[Any] = []
    uncertain = bool(metadata.get("uncertainty", False))
    for item in inputs:
        info = _origin_info(item)
        ids.extend(info["source_ids"])
        classes.extend(info["error_classes"])
        uncertain = uncertain or bool(info["uncertainty"])
    ids.extend(_as_items(metadata.get("source_ids")))
    classes.extend(_as_items(metadata.get("error_classes")))
    origin: dict[str, Any] = {
        "operation": operation,
        "source_ids": _unique(ids),
        "error_classes": _unique(classes),
        "uncertainty": bool(uncertain),
        "symbolic_cancellation_claimed": False,
    }
    if inputs:
        origin["parents"] = tuple(_origin_ref(item) for item in inputs)
    for key, value in metadata.items():
        if key not in {"source_ids", "error_classes", "uncertainty"}:
            origin[key] = value
    return origin


def _check_powers(value: Any, label: str) -> tuple[int, int, int, int]:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes, bytearray))
        or len(value) != 4
        or any(not isinstance(power, int) or isinstance(power, bool) for power in value)
    ):
        raise ValueError(f"{label} must contain four integer formal powers")
    return tuple(int(power) for power in value)  # type: ignore[return-value]


def _legacy_powers(value: Mapping[str, Any], label: str) -> tuple[int, int, int, int] | None:
    """Accept the older ``F0_power``/``Lambda_power`` spelling as a bridge."""

    if "powers" in value:
        return _check_powers(value["powers"], f"{label}.powers")
    if "F0_power" not in value and "Lambda_power" not in value:
        return None
    f0 = value.get("F0_power", 0)
    lam = value.get("Lambda_power", 0)
    hb = value.get("hb_power", 0)
    pstar = value.get("Pstar_squared_power", value.get("Pstar_power", 0))
    powers = (hb, f0, pstar, lam)
    return _check_powers(powers, f"{label}.formal powers")


def _check_term(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a sector mapping")
    powers = _legacy_powers(value, label)
    if powers is None:
        raise ValueError(f"{label} requires powers or F0_power/Lambda_power")
    if "coefficient" not in value:
        raise ValueError(f"{label} requires coefficient")
    if not _is_interval(value["coefficient"]):
        raise TypeError(f"{label}.coefficient must be a directed interval")
    if "origin" not in value:
        raise ValueError(f"{label} requires origin provenance")
    info = _origin_info(value["origin"])
    # A runner may keep the source identity on the sector itself as well as
    # inside origin.  Include both, without inventing an ID when neither was
    # supplied.
    ids = list(info["source_ids"])
    for key in ("source_ids", "source_id", "input_id"):
        if value.get(key) is not None:
            ids.extend(_as_items(value[key]))
    ids = _unique(ids)
    if not ids:
        raise ValueError(f"{label}.origin must retain source_ids")
    origin_value = value["origin"]
    if isinstance(origin_value, Mapping):
        origin: Any = dict(origin_value)
    else:
        origin = {}
    origin["source_ids"] = list(ids)
    origin["error_classes"] = list(info["error_classes"])
    origin["uncertainty"] = bool(info["uncertainty"])
    return {
        "powers": powers,
        "coefficient": value["coefficient"],
        "origin": origin,
        "source_ids": ids,
        "error_classes": list(info["error_classes"]),
        "uncertainty": bool(info["uncertainty"]),
    }


def _term_list(value: Any, label: str, *, allow_empty: bool = False) -> list[dict[str, Any]]:
    if isinstance(value, Mapping) and ("powers" in value or "coefficient" in value):
        values = [value]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
    else:
        raise TypeError(f"{label} must be a sector sequence")
    if not values and not allow_empty:
        raise ValueError(f"{label} must not be empty")
    return [_check_term(item, f"{label}[{index}]") for index, item in enumerate(values)]


def _order_key(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"boolean is not an ordinary-Z order in {label}")
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        try:
            return int(value)
        except ValueError as exc:
            raise ValueError(f"invalid ordinary-Z order {value!r} in {label}") from exc
    raise ValueError(f"ordinary-Z order required in {label}")


def _orders(value: Any, label: str) -> dict[int, Any]:
    if isinstance(value, Mapping):
        result: dict[int, Any] = {}
        for raw_key, item in value.items():
            key = _order_key(raw_key, label)
            if key in result:
                raise ValueError(f"duplicate ordinary-Z order {key} in {label}")
            result[key] = item
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        result = {index: item for index, item in enumerate(value)}
    else:
        raise TypeError(f"{label} must provide ordinary-Z order rows")
    if set(result) != set(ORDINARY_Z_ORDERS):
        raise ValueError(f"{label} must provide exactly Z0 and Z1")
    return result


def _normalize_leading(value: Any, label: str) -> dict[int, list[dict[str, Any]]]:
    rows = _orders(value, label)
    return {order: _term_list(rows[order], f"{label}[{order}]") for order in ORDINARY_Z_ORDERS}


def _lookup_column(row: Mapping[Any, Any], names: Sequence[Any], label: str) -> Any:
    found: list[Any] = []
    for name in names:
        if name in row:
            found.append(row[name])
    if not found:
        raise ValueError(f"{label} is missing required W0 column")
    if len(found) > 1:
        # Numeric and named aliases are accepted when they point at the same
        # sector object.  Avoid equality on interval/list payloads here: that
        # can itself invoke arithmetic or an ambiguous interval comparison.
        first = found[0]
        if any(item is not first for item in found[1:]):
            raise ValueError(f"{label} has duplicate aliases for one W0 column")
    return found[0]


def _normalize_outer_cells(c: Any, value: Any) -> list[dict[str, Any]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise TypeError("outer_cells must be a sequence of eight saved cells")
    if len(value) != OUTER_CELL_COUNT:
        raise ValueError("outer_cells must contain exactly eight cells")

    normalized: list[dict[str, Any]] = []
    for index, cell in enumerate(value):
        label = f"outer_cells[{index}]"
        if not isinstance(cell, Mapping):
            raise TypeError(f"{label} must be a cell mapping")
        if "t" not in cell or "rho" not in cell or "W0_ordinary_Z" not in cell:
            raise ValueError(f"{label} requires t, rho, and W0_ordinary_Z")
        t = cell["t"]
        rho = cell["rho"]
        if not _is_interval(t) or not _is_interval(rho):
            raise TypeError(f"{label}.t and .rho must be directed intervals")

        left = c.mpf(index) / c.mpf(8)
        right = c.mpf(index + 1) / c.mpf(8)
        expected_t = _make_interval(c, left, right)
        expected_rho = 4 * expected_t * expected_t
        if not _same_interval(t, expected_t) or not _same_interval(rho, expected_rho):
            raise ValueError(f"{label} is not the directed rho=4*t^2 outer cell")

        rows = _orders(cell["W0_ordinary_Z"], f"{label}.W0_ordinary_Z")
        fields: dict[int, dict[str, list[dict[str, Any]]]] = {}
        for order in ORDINARY_Z_ORDERS:
            row = rows[order]
            if not isinstance(row, Mapping):
                raise TypeError(f"{label}.W0_ordinary_Z[{order}] must be a column mapping")
            fields[order] = {
                "F1": _term_list(
                    _lookup_column(row, (1, "1", "F1"), f"{label}.W0[{order}].F1"),
                    f"{label}.W0[{order}].F1",
                ),
                "Uz1": _term_list(
                    _lookup_column(row, (2, "2", "Uz1", "UZ1", "U_z1"), f"{label}.W0[{order}].Uz1"),
                    f"{label}.W0[{order}].Uz1",
                ),
            }
        normalized.append({"index": index, "t": t, "rho": rho, "fields": fields, "raw": cell})
    return normalized


def _validate_source_atoms(source_atoms: Any) -> None:
    names: set[Any] = set()
    if isinstance(source_atoms, Mapping):
        names = set(source_atoms)
    elif isinstance(source_atoms, Sequence) and not isinstance(source_atoms, (str, bytes, bytearray)):
        for item in source_atoms:
            if isinstance(item, Mapping):
                name = item.get("name", item.get("atom", item.get("power")))
                if name is not None:
                    names.add(name)
            else:
                names.add(item)
    else:
        raise ValueError("source_atoms must name the four formal source atoms")
    if names != set(SOURCE_POWER_BASIS):
        raise ValueError(
            "source_atoms must contain exactly hb, F0_anchor, Pstar_squared, Lambda"
        )


def _validate_certificate(certificate: Any) -> None:
    if not isinstance(certificate, Mapping):
        raise ValueError("caller geometry and source certificates are required")
    missing = [name for name in REQUIRED_CERTIFICATES if certificate.get(name) is not True]
    if missing:
        raise ValueError("missing caller certificate: " + ", ".join(missing))


def _group_key(term: Mapping[str, Any]) -> tuple[Any, bool, tuple[Any, ...]]:
    info = _origin_info(term["origin"])
    classes = tuple(sorted((_freeze(item) for item in info["error_classes"]), key=repr))
    return tuple(term["powers"]), bool(info["uncertainty"]), classes


def _aggregate(terms: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Outward-add only equal-power/equal-error-class sectors."""

    groups: dict[tuple[Any, bool, tuple[Any, ...]], dict[str, Any]] = {}
    for term in terms:
        info = _origin_info(term["origin"])
        key = _group_key(term)
        group = groups.get(key)
        if group is None:
            group = {
                "powers": tuple(term["powers"]),
                "coefficient": term["coefficient"],
                "ids": {},
                "classes": {},
                "uncertainty": False,
                "merged_term_count": 0,
            }
            groups[key] = group
        else:
            group["coefficient"] = group["coefficient"] + term["coefficient"]
        for value in info["source_ids"]:
            group["ids"][_freeze(value)] = value
        for value in info["error_classes"]:
            group["classes"][_freeze(value)] = value
        group["uncertainty"] = group["uncertainty"] or bool(info["uncertainty"])
        group["merged_term_count"] += 1

    result: list[dict[str, Any]] = []
    for group in groups.values():
        powers = tuple(group["powers"])
        origin = {
            "operation": "aggregate",
            "source_ids": list(group["ids"].values()),
            "error_classes": list(group["classes"].values()),
            "uncertainty": bool(group["uncertainty"]),
            "merged_term_count": group["merged_term_count"],
            "symbolic_cancellation_claimed": False,
        }
        result.append(
            {
                "powers": powers,
                "F0_power": powers[1],
                "Lambda_power": powers[3],
                "coefficient": group["coefficient"],
                "origin": origin,
            }
        )
    return result


def _product(
    leading: Mapping[str, Any],
    word: Mapping[str, Any],
    factor: Any,
    *,
    powers_offset: tuple[int, int, int, int],
    channel: str,
    order: int,
    cell: int,
    derivative_piece: str,
) -> dict[str, Any]:
    powers = tuple(
        a + b + offset
        for a, b, offset in zip(leading["powers"], word["powers"], powers_offset)
    )
    return {
        "powers": powers,
        "coefficient": factor * leading["coefficient"] * word["coefficient"],
        "origin": _new_origin(
            "core_w0_m15_product",
            leading["origin"],
            word["origin"],
            channel=channel,
            ordinary_Z_order=order,
            core_cell=cell,
            derivative_piece=derivative_piece,
            external_formal_powers=powers_offset,
            symbolic_cancellation_claimed=False,
        ),
    }


def _cell_weight(c: Any, cell: int, channel: str) -> Any:
    left = c.mpf(cell) / c.mpf(8)
    right = c.mpf(cell + 1) / c.mpf(8)
    if channel == "axial":
        return c.mpf(8) * (right**2 - left**2)
    if channel == "swirl":
        return -c.mpf(16) * (right**4 - left**4)
    raise ValueError(f"unknown core W0 channel {channel!r}")


def _assemble_velocity(
    c: Any,
    cells: Sequence[Mapping[str, Any]],
    leading: Mapping[int, list[dict[str, Any]]],
    *,
    field: str,
    channel: str,
    external_powers: tuple[int, int, int, int],
) -> dict[int, dict[int, list[dict[str, Any]]]]:
    output: dict[int, dict[int, list[dict[str, Any]]]] = {
        order: {index: [] for index in range(OUTER_CELL_COUNT)}
        for order in ORDINARY_Z_ORDERS
    }
    for cell in cells:
        index = int(cell["index"])
        factor = _cell_weight(c, index, channel)
        word0 = cell["fields"][0][field]
        word1 = cell["fields"][1][field]
        raw_by_order: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
        for left in leading[0]:
            for right in word0:
                raw_by_order[0].append(
                    _product(
                        left,
                        right,
                        factor,
                        powers_offset=external_powers,
                        channel=channel,
                        order=0,
                        cell=index,
                        derivative_piece="leading0*word0",
                    )
                )
        for left in leading[1]:
            for right in word0:
                raw_by_order[1].append(
                    _product(
                        left,
                        right,
                        factor,
                        powers_offset=external_powers,
                        channel=channel,
                        order=1,
                        cell=index,
                        derivative_piece="leading1*word0",
                    )
                )
        for left in leading[0]:
            for right in word1:
                raw_by_order[1].append(
                    _product(
                        left,
                        right,
                        factor,
                        powers_offset=external_powers,
                        channel=channel,
                        order=1,
                        cell=index,
                        derivative_piece="leading0*word1",
                    )
                )
        for order in ORDINARY_Z_ORDERS:
            output[order][index] = _aggregate(raw_by_order[order])
    return output


def _pressure_factor(c: Any) -> Any:
    # -epsilon * 8 * Integral_0^1 t^2 dt = -epsilon * 8/3.
    return -c.mpf(8) / c.mpf(3)


def _assemble_pressure(
    c: Any,
    g4: Mapping[int, list[dict[str, Any]]],
) -> dict[int, list[dict[str, Any]]]:
    factor = _pressure_factor(c)
    output: dict[int, list[dict[str, Any]]] = {}
    for order in ORDINARY_Z_ORDERS:
        raw: list[dict[str, Any]] = []
        for index, source in enumerate(g4[order]):
            powers = tuple(source["powers"][axis] + (0, 0, 0, -1)[axis] for axis in range(4))
            raw.append(
                {
                    "powers": powers,
                    "coefficient": factor * source["coefficient"],
                    "origin": _new_origin(
                        "core_w0_m15_pressure_product",
                        source["origin"],
                        ordinary_Z_order=order,
                        g4_input_term_index=index,
                        external_formal_powers=(0, 0, 0, -1),
                        uniform_g4_box=True,
                        pressure_mass="8/3",
                        symbolic_cancellation_claimed=False,
                    ),
                }
            )
        output[order] = _aggregate(raw)
    return output


def build(
    c: Any,
    *,
    outer_cells: Any,
    leading_U: Any,
    leading_F: Any,
    g4: Any,
    source_family: Any,
    source_atoms: Any,
    certificate: Any,
) -> dict[str, Any]:
    """Build finite core W0 axial, swirl, and leading pressure sectors.

    ``outer_cells`` contains eight exact ``t`` cells.  Their
    ``W0_ordinary_Z`` rows use columns ``1/F1`` and ``2/Uz1`` and provide
    saved cell enclosures, not selected point values.  ``leading_U`` and
    ``leading_F`` provide physical whole-core rows for ordinary ``Z0/Z1``;
    the supplied physical F row already contains its F0 connection.  The
    supplied ``g4`` rows are consumed directly for pressure, with no extra
    epsilon or ``Pstar`` recomputation.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must be the caller's MPIntervalContext")
    if source_family is None or (isinstance(source_family, str) and not source_family.strip()):
        raise ValueError("source_family is required")
    _validate_source_atoms(source_atoms)
    _validate_certificate(certificate)

    cells = _normalize_outer_cells(c, outer_cells)
    leading_u = _normalize_leading(leading_U, "leading_U")
    leading_f = _normalize_leading(leading_F, "leading_F")
    g4_rows = _normalize_leading(g4, "g4")

    axial = _assemble_velocity(
        c,
        cells,
        leading_u,
        field="Uz1",
        channel="axial",
        external_powers=(0, 0, 0, -1),
    )
    swirl = _assemble_velocity(
        c,
        cells,
        leading_f,
        field="F1",
        channel="swirl",
        external_powers=(0, 0, 0, -2),
    )
    pressure = _assemble_pressure(c, g4_rows)

    return {
        "source_family": source_family,
        "source_atoms": source_atoms,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "outer_cells": outer_cells,
        "outer_cell_count": OUTER_CELL_COUNT,
        "axial": axial,
        "swirl": swirl,
        "pressure": pressure,
        "axial_weight_formula": "8*(t_hi^2-t_lo^2)",
        "swirl_weight_formula": "-16*(t_hi^4-t_lo^4)",
        "pressure_weight": c.mpf(8) / c.mpf(3),
        "pressure_signed_weight_formula": "-8/3",
        "axial_formula": "epsilon*integral(2*U_lead*Uz1_W0*drho)",
        "swirl_formula": "-epsilon^2*integral(2*rho*F_lead*F1_W0*drho)",
        "pressure_formula": "-epsilon*8*integral(t^2*g4*dt)",
        "formal_axial_external_powers": (0, 0, 0, -1),
        "formal_swirl_external_powers": (0, 0, 0, -2),
        "formal_pressure_external_powers": (0, 0, 0, -1),
        "rho_formula": "4*t^2",
        "g4_identity": "g4=2*epsilon*sqrt(rho)*N1p",
        "g4_consumed_directly": True,
        "ordinary_Z1_product_rule": True,
        "leading_F_is_physical_with_F0_connection": True,
        "F_Z_connection_added_here": False,
        "pressure_endpoint_added": False,
        "endpoint_added": False,
        "eta_applied": False,
        "W1_included": False,
        "w1_included": False,
        "allword_tail_added": False,
        "allword_tail_included": False,
        "tail_added": False,
        "tail_included": False,
        "ledger_consulted": False,
        "ledger_212_consulted": False,
        "no_ledger": True,
        "complete_m15_computed": False,
        "complete_m15": False,
        "m15_computed": False,
        "total_m15_computed": False,
        "source_ids_union_retained": True,
        "error_classes_union_retained": True,
        "uncertainty_retained": True,
        "aggregation_key": "(powers, uncertainty, error_classes)",
        "aggregation_is_outward_interval_addition_only": True,
        "same_source_cancellation_claimed": False,
        "symbolic_cancellation_claimed": False,
        "point_value_computed": False,
        "numeric_point_evaluation": False,
        "geometry_assumptions_caller_certified": True,
        "geometry_certificate_independently_verified": False,
        "scope": "finite core W0 physical m15 axial, swirl, and raw g4 pressure sectors only",
    }
