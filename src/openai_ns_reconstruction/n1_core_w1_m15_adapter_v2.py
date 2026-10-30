"""Independent finite core ``W1`` axial and swirl ``m15`` assembly.

This leaf consumes already-normalized, interval-valued target cells.  Each of
the eight cells carries the generic ``W0_ordinary_Z`` field name used by the
current profile handoff; semantically those fields are the saved core ``W1``
``F1`` and ``Uz1`` sectors.  The name is retained only for compatibility with
that handoff and does not mean that ``W0`` is reconstructed here.

The physical cell factors are assembled locally from ``rho = 4*t**2``:

``axial = epsilon * integral(2*U_lead*Uz1_W1*drho)``

``swirl = -epsilon**2 * integral(2*rho*F_lead*F1_W1*drho)``

The saved sector coefficient is an enclosure for a target cell, so its
coefficient is multiplied by the exact directed cell mass
``8*(t_hi**2-t_lo**2)`` for axial and
``-16*(t_hi**4-t_lo**4)`` for swirl.  Ordinary ``Z1`` uses the two-term
Leibniz rule, ``leading_Z1*W1_Z0 + leading_Z0*W1_Z1``.  Products and sums
retain all four formal powers and union all source IDs/error classes through
provenance; no ``F0`` connection or source value is recomputed.

There is intentionally no ``g4`` argument and no pressure assembly.  This
independent implementation therefore cannot duplicate a pressure contribution
that another stage may own.  It also does not import or call a previous-round
builder or any historical source producer.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any


COMPONENT = "n1_core_w1_m15_adapter"
SCHEMA = "n1_core_w1_m15_adapter_v2"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
OUTER_CELL_COUNT = 8
FIELDS = ("F1", "Uz1")
REQUIRED_CERTIFICATES = ("geometry_epsilon_Z_independent",)


__all__ = [
    "COMPONENT",
    "SCHEMA",
    "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "OUTER_CELL_COUNT",
    "FIELDS",
    "REQUIRED_CERTIFICATES",
    "build",
    "build_w1_m15",
    "adapt",
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
    """Construct a directed interval without midpoint conversion."""

    return c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))


def _same_interval(left: Any, right: Any) -> bool:
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
    """Collect source IDs, uncertainty, and error classes recursively."""

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
    """Union provenance for one product without inventing source identity."""

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
    """Accept legacy names while canonicalizing every term to four powers."""

    if "powers" in value:
        return _check_powers(value["powers"], f"{label}.powers")
    if "F0_power" not in value and "Lambda_power" not in value:
        return None
    powers = (
        value.get("hb_power", 0),
        value.get("F0_power", 0),
        value.get("Pstar_squared_power", value.get("Pstar_power", 0)),
        value.get("Lambda_power", 0),
    )
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
    ids = list(info["source_ids"])
    for key in ("source_ids", "source_id", "input_id"):
        if value.get(key) is not None:
            ids.extend(_as_items(value[key]))
    ids = _unique(ids)
    if not ids:
        raise ValueError(f"{label}.origin must retain source_ids")
    origin_value = value["origin"]
    origin = dict(origin_value) if isinstance(origin_value, Mapping) else {}
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
        raise ValueError(f"{label} is missing required W1 column")
    if len(found) > 1 and any(item is not found[0] for item in found[1:]):
        raise ValueError(f"{label} has duplicate aliases for one W1 column")
    return found[0]


def _normalize_profile(c: Any, value: Any) -> list[dict[str, Any]]:
    if isinstance(value, Mapping):
        for wrapper in ("outer_cells", "profile", "target_cells"):
            if wrapper in value:
                value = value[wrapper]
                break
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise TypeError("profile must be a sequence of eight target cells")
    if len(value) != OUTER_CELL_COUNT:
        raise ValueError("profile must contain exactly eight target cells")

    normalized: list[dict[str, Any]] = []
    for index, cell in enumerate(value):
        label = f"profile[{index}]"
        if not isinstance(cell, Mapping):
            raise TypeError(f"{label} must be a cell mapping")
        if "t" not in cell or "rho" not in cell:
            raise ValueError(f"{label} requires t and rho")
        t, rho = cell["t"], cell["rho"]
        if not _is_interval(t) or not _is_interval(rho):
            raise TypeError(f"{label}.t and .rho must be directed intervals")
        left = c.mpf(index) / c.mpf(OUTER_CELL_COUNT)
        right = c.mpf(index + 1) / c.mpf(OUTER_CELL_COUNT)
        expected_t = _make_interval(c, left, right)
        expected_rho = c.mpf(4) * expected_t * expected_t
        if not _same_interval(t, expected_t) or not _same_interval(rho, expected_rho):
            raise ValueError(f"{label} is not the directed rho=4*t^2 target cell")

        generic = cell.get("W0_ordinary_Z")
        if generic is None:
            for alias in ("W1_ordinary_Z", "W1", "ordinary_Z"):
                if alias in cell:
                    generic = cell[alias]
                    break
        if generic is None:
            raise ValueError(
                f"{label} requires W0_ordinary_Z (the generic saved W1 profile field)"
            )
        rows = _orders(generic, f"{label}.W1_ordinary_Z")
        fields: dict[int, dict[str, list[dict[str, Any]]]] = {}
        for order in ORDINARY_Z_ORDERS:
            row = rows[order]
            if not isinstance(row, Mapping):
                raise TypeError(f"{label}.W1_ordinary_Z[{order}] must be a column mapping")
            fields[order] = {
                "F1": _term_list(
                    _lookup_column(row, (1, "1", "F1"), f"{label}.W1[{order}].F1"),
                    f"{label}.W1[{order}].F1",
                ),
                "Uz1": _term_list(
                    _lookup_column(row, (2, "2", "Uz1", "UZ1", "U_z1"), f"{label}.W1[{order}].Uz1"),
                    f"{label}.W1[{order}].Uz1",
                ),
            }
        normalized.append(
            {"index": index, "t": t, "rho": rho, "fields": fields, "raw": cell}
        )
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
        raise ValueError("caller geometry certificate is required")
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
        result.append(
            {
                "powers": powers,
                "coefficient": group["coefficient"],
                "origin": {
                    "operation": "core_w1_m15_aggregate",
                    "source_ids": list(group["ids"].values()),
                    "error_classes": list(group["classes"].values()),
                    "uncertainty": bool(group["uncertainty"]),
                    "merged_term_count": group["merged_term_count"],
                    "symbolic_cancellation_claimed": False,
                },
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
            "core_w1_m15_product",
            leading["origin"],
            word["origin"],
            channel=channel,
            ordinary_Z_order=order,
            target_cell=cell,
            derivative_piece=derivative_piece,
            external_formal_powers=powers_offset,
            symbolic_cancellation_claimed=False,
        ),
    }


def _cell_weight(c: Any, cell: int, channel: str) -> Any:
    left = c.mpf(cell) / c.mpf(OUTER_CELL_COUNT)
    right = c.mpf(cell + 1) / c.mpf(OUTER_CELL_COUNT)
    if channel == "axial":
        return c.mpf(8) * (right**2 - left**2)
    if channel == "swirl":
        return -c.mpf(16) * (right**4 - left**4)
    raise ValueError(f"unknown core W1 channel {channel!r}")


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
                        derivative_piece="leading0*W1_0",
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
                        derivative_piece="leading1*W1_0",
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
                        derivative_piece="leading0*W1_1",
                    )
                )
        for order in ORDINARY_Z_ORDERS:
            output[order][index] = _aggregate(raw_by_order[order])
    return output


def _coalesce(value: Any, aliases: Mapping[str, Any], names: Sequence[str], label: str) -> Any:
    candidates = [] if value is None else [(label, value)]
    candidates.extend((name, aliases[name]) for name in names if aliases.get(name) is not None)
    if not candidates:
        return None
    first_name, first = candidates[0]
    for name, candidate in candidates[1:]:
        if candidate is not first:
            raise TypeError(f"{label} was supplied more than once ({first_name}, {name})")
    return first


def build(
    c: Any,
    *,
    profile: Any | None = None,
    leading_U: Any | None = None,
    leading_F: Any | None = None,
    source_family: Any | None = None,
    source_atoms: Any | None = None,
    certificate: Any | None = None,
    **aliases: Any,
) -> dict[str, Any]:
    """Build finite core ``W1`` axial and swirl sectors only.

    ``profile`` is the caller-normalized eight-cell structure.  For each cell,
    ``W0_ordinary_Z`` (or ``W1_ordinary_Z``) contains order ``0``/``1`` rows,
    where column ``1`` is ``F1`` and column ``2`` is ``Uz1``.  ``leading_U``
    and ``leading_F`` are physical leading rows and are forwarded only as
    interval data to the local product operation; their F0 connection is not
    reconstructed.
    """

    cells = _coalesce(
        profile,
        aliases,
        ("outer_cells", "w1_profile", "target_cells"),
        "profile",
    )
    unexpected = set(aliases) - {"outer_cells", "w1_profile", "target_cells"}
    if unexpected:
        names = ", ".join(sorted(str(name) for name in unexpected))
        if "g4" in unexpected:
            raise ValueError("g4 is not part of the independent W1 velocity adapter")
        raise TypeError(f"unexpected adapter arguments: {names}")

    if cells is None:
        raise TypeError("profile is required")
    if leading_U is None or leading_F is None:
        raise TypeError("leading_U and leading_F are required")
    if source_family is None or (isinstance(source_family, str) and not source_family.strip()):
        raise ValueError("source_family is required")
    if source_atoms is None:
        raise ValueError("source_atoms is required")
    if certificate is None:
        raise ValueError("certificate is required")
    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must be the caller's MPIntervalContext")

    _validate_source_atoms(source_atoms)
    _validate_certificate(certificate)
    normalized_cells = _normalize_profile(c, cells)
    leading_u = _normalize_leading(leading_U, "leading_U")
    leading_f = _normalize_leading(leading_F, "leading_F")

    axial = _assemble_velocity(
        c,
        normalized_cells,
        leading_u,
        field="Uz1",
        channel="axial",
        external_powers=(0, 0, 0, -1),
    )
    swirl = _assemble_velocity(
        c,
        normalized_cells,
        leading_f,
        field="F1",
        channel="swirl",
        external_powers=(0, 0, 0, -2),
    )

    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "source_family": source_family,
        "source_atoms": source_atoms,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "outer_cell_count": OUTER_CELL_COUNT,
        "axial": axial,
        "swirl": swirl,
        "axial_weight_formula": "8*(t_hi^2-t_lo^2)",
        "swirl_weight_formula": "-16*(t_hi^4-t_lo^4)",
        "axial_formula": "epsilon*integral(2*U_lead*Uz1_W1*drho)",
        "swirl_formula": "-epsilon^2*integral(2*rho*F_lead*F1_W1*drho)",
        "formal_axial_external_powers": (0, 0, 0, -1),
        "formal_swirl_external_powers": (0, 0, 0, -2),
        "rho_formula": "4*t^2",
        "profile_field_name": "W0_ordinary_Z",
        "profile_semantic_role": "saved core W1 target-cell F1/Uz1 sectors",
        "ordinary_Z1_product_rule": True,
        "leading_F_is_physical_with_F0_connection": True,
        "F_Z_connection_added_here": False,
        "F0_recomputed": False,
        "F0_recalculated": False,
        "source_ids_union_retained": True,
        "error_classes_union_retained": True,
        "uncertainty_retained": True,
        "formal_powers_retained": True,
        "aggregation_key": "(powers, uncertainty, error_classes)",
        "aggregation_is_outward_interval_addition_only": True,
        "same_source_cancellation_claimed": False,
        "symbolic_cancellation_claimed": False,
        "g4_used": False,
        "g4_pressure_added": False,
        "g4_pressure_duplicate_prevented": True,
        "pressure_included": False,
        "pressure_added": False,
        "only_axial_and_swirl_returned": True,
        "W1_included": True,
        "w1_included": True,
        "allword_tail_added": False,
        "allword_tail_included": False,
        "tail_added": False,
        "tail_included": False,
        "complete_m15_computed": False,
        "complete_m15": False,
        "m15_computed": False,
        "total_m15_computed": False,
        "geometry_assumptions_caller_certified": True,
        "geometry_certificate_independently_verified": False,
        "scope": "finite core W1 axial and swirl m15 sectors only; no pressure",
    }


build_w1_m15 = build
adapt = build
