"""Finite physical ``W0`` first-collar ``m15`` source sectors.

This module is deliberately a source-only adapter.  The caller supplies the
physical whole-phase leading ``U`` and ``F`` sectors and the ordinary-Z0/Z1
``W0`` sectors.  No producer, archive, JSON receipt, W1 ledger, or all-word
tail is consulted here.

The flat complementary cutoff is

``eta(q) = 1 - sigma(4*q - 2)``,

with exact cell masses ``(1/4, 1/4, 1/8, 0)`` and total mass ``5/8``.  The
whole-phase core is integrated once with the total mass.  Local sectors are
kept independent on active cells 0, 1, and 2; cell 3 is absent because its
cutoff mass is exactly zero.

All coefficients remain directed intervals.  Products add formal powers,
and aggregation is performed only for equal powers, uncertainty flag, and
error-class set.  Source IDs and error classes are unions, never a
same-source cancellation claim.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
FIELDS = ("Uz1", "F1")
BRANCHES = ("core", "local")
ALL_LOCAL_CELLS = (0, 1, 2, 3)
ACTIVE_LOCAL_CELLS = (0, 1, 2)

__all__ = [
    "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "FIELDS",
    "BRANCHES",
    "ALL_LOCAL_CELLS",
    "ACTIVE_LOCAL_CELLS",
    "build",
]


def _is_interval(value: Any) -> bool:
    """Return whether *value* is a caller-owned directed interval."""

    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _raw_endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    return value


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    """Construct an interval from directed endpoints without midpointing."""

    return c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))


def _scalar_bounds(value: Any) -> tuple[Any, Any]:
    """Expose raw directed endpoints for validation only.

    The interval context's ``.a`` and ``.b`` are singleton interval objects;
    the raw ``_mpi_`` endpoints are the exact directed bounds and are never
    replaced by a midpoint or a rounded scalar.
    """

    raw = getattr(value, "_mpi_", None)
    if raw is None:
        return value, value
    # ``iv.mpf._mpi_`` stores the low-level endpoint tuples, not ordinary
    # scalar MPF objects.  Convert only for comparisons; coefficients and
    # factors themselves stay in the caller's interval context.
    from mpmath import mp

    return mp.make_mpf(raw[0]), mp.make_mpf(raw[1])


def _as_items(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    return [value]


def _freeze(value: Any) -> Any:
    """Make an arbitrary provenance value usable as a grouping key."""

    try:
        hash(value)
    except TypeError:
        return (type(value).__name__, repr(value))
    return (type(value).__name__, value)


def _unique(values: Sequence[Any]) -> list[Any]:
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
        return {
            "source_ids": [],
            "error_classes": [],
            "uncertainty": False,
            "path_count": 1,
            "source_count": 1,
        }

    if isinstance(value, Mapping):
        _seen.add(marker)
        ids = list(_as_items(value.get("source_ids", value.get("sources"))))
        classes = list(_as_items(value.get(
            "error_classes",
            value.get("uncertainty_class", value.get("error_class")),
        )))
        operation = str(value.get("operation", "")).lower()
        error_kind = str(value.get("error_kind", value.get("kind", ""))).lower()
        uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False)))
        uncertain = uncertain or any(
            token in operation or token in error_kind
            for token in ("cauchy", "symmetric", "uncertainty")
        )

        nested: list[dict[str, Any]] = []
        for key in ("origin", "metadata", "parents"):
            child = value.get(key)
            if child is None or child is value:
                continue
            if isinstance(child, Sequence) and not isinstance(child, (str, bytes, bytearray)):
                nested.extend(_origin_info(item, _seen=_seen) for item in child)
            else:
                nested.append(_origin_info(child, _seen=_seen))
        for info in nested:
            ids.extend(info["source_ids"])
            classes.extend(info["error_classes"])
            uncertain = uncertain or bool(info["uncertainty"])

        raw_path = value.get("path_count", value.get("aggregate_path_count"))
        raw_source = value.get("source_count", value.get("aggregate_source_count"))
        try:
            path_count = max(1, int(raw_path)) if raw_path is not None else 1
        except (TypeError, ValueError):
            path_count = 1
        try:
            source_count = max(1, int(raw_source)) if raw_source is not None else path_count
        except (TypeError, ValueError):
            source_count = path_count
        if nested:
            path_count = max(path_count, max(int(item["path_count"]) for item in nested))
            source_count = max(source_count, max(int(item["source_count"]) for item in nested))
        return {
            "source_ids": _unique(ids),
            "error_classes": _unique(classes),
            "uncertainty": bool(uncertain),
            "path_count": path_count,
            "source_count": source_count,
        }

    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        nested = [_origin_info(item, _seen=_seen) for item in value]
        if nested:
            ids: list[Any] = []
            classes: list[Any] = []
            uncertain = False
            path_count = 1
            source_count = 1
            for info in nested:
                ids.extend(info["source_ids"])
                classes.extend(info["error_classes"])
                uncertain = uncertain or bool(info["uncertainty"])
                path_count = max(path_count, int(info["path_count"]))
                source_count = max(source_count, int(info["source_count"]))
            return {
                "source_ids": _unique(ids),
                "error_classes": _unique(classes),
                "uncertainty": uncertain,
                "path_count": path_count,
                "source_count": source_count,
            }

    return {
        "source_ids": [],
        "error_classes": [],
        "uncertainty": False,
        "path_count": 1,
        "source_count": 1,
    }


def _origin_ref(value: Any) -> dict[str, Any]:
    info = _origin_info(value)
    return {
        "source_ids": list(info["source_ids"]),
        "error_classes": list(info["error_classes"]),
        "uncertainty": bool(info["uncertainty"]),
        "path_count": int(info["path_count"]),
        "source_count": int(info["source_count"]),
    }


def _new_origin(operation: str, *inputs: Any, **metadata: Any) -> dict[str, Any]:
    """Build a provenance union for a product or aggregate."""

    infos = [_origin_info(item) for item in inputs]
    ids: list[Any] = []
    classes: list[Any] = []
    uncertain = bool(metadata.get("uncertainty", False))
    if operation.endswith("_product") or operation in {"product", "multiply"}:
        path_count = 1
        source_count = 1
        for info in infos:
            path_count *= int(info["path_count"])
            source_count *= int(info["source_count"])
    elif operation in {"aggregate", "add", "sum"}:
        path_count = 0
        source_count = 0
        for info in infos:
            path_count += int(info["path_count"])
            source_count += int(info["source_count"])
    else:
        path_count = int(infos[0]["path_count"]) if infos else 1
        source_count = int(infos[0]["source_count"]) if infos else 1

    for info in infos:
        ids.extend(info["source_ids"])
        classes.extend(info["error_classes"])
        uncertain = uncertain or bool(info["uncertainty"])
    if "source_ids" in metadata:
        ids.extend(_as_items(metadata["source_ids"]))
    if "error_classes" in metadata:
        classes.extend(_as_items(metadata["error_classes"]))

    origin: dict[str, Any] = {
        "operation": operation,
        "source_ids": _unique(ids),
        "error_classes": _unique(classes),
        "uncertainty": bool(uncertain),
        "path_count": max(1, int(path_count)),
        "source_count": max(1, int(source_count)),
        "source_id_count": len(_unique(ids)),
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


def _check_sector(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a sector mapping")
    if not {"powers", "coefficient", "origin"}.issubset(value):
        raise ValueError(f"{label} requires powers, coefficient, and origin")
    if not _is_interval(value["coefficient"]):
        raise TypeError(f"{label}.coefficient must be a directed interval")
    powers = _check_powers(value["powers"], f"{label}.powers")
    info = _origin_info(value["origin"])
    if not info["source_ids"]:
        raise ValueError(f"{label}.origin must retain source_ids")
    return {
        "powers": powers,
        "coefficient": value["coefficient"],
        "origin": value["origin"],
    }


def _term_list(
    value: Any,
    label: str,
    *,
    allow_empty: bool = False,
) -> list[dict[str, Any]]:
    if isinstance(value, Mapping) and "powers" in value:
        values = [value]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
    else:
        raise TypeError(f"{label} must be a sector sequence")
    if not values and not allow_empty:
        raise ValueError(f"{label} must not be empty")
    return [_check_sector(item, f"{label}[{index}]") for index, item in enumerate(values)]


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
        raise TypeError(f"{label} must be an ordinary-Z order mapping")
    if set(result) != set(ORDINARY_Z_ORDERS):
        raise ValueError(f"{label} must provide exactly Z0 and Z1")
    return result


def _normalize_leading(value: Any, label: str) -> dict[int, list[dict[str, Any]]]:
    raw = _orders(value, label)
    return {
        order: _term_list(raw[order], f"{label}[{order}]")
        for order in ORDINARY_Z_ORDERS
    }


def _cell_orders(value: Any, label: str) -> dict[int, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a cell mapping")
    result: dict[int, Any] = {}
    for raw_key, item in value.items():
        key = _order_key(raw_key, label)
        if key in result:
            raise ValueError(f"duplicate cell {key} in {label}")
        result[key] = item
    if set(result) != set(ALL_LOCAL_CELLS):
        raise ValueError(f"{label} must provide cells 0, 1, 2, and 3")
    return result


def _normalize_w0(value: Any) -> dict[int, dict[str, dict[str, Any]]]:
    raw = _orders(value, "w0")
    result: dict[int, dict[str, dict[str, Any]]] = {}
    for order in ORDINARY_Z_ORDERS:
        row = raw[order]
        if not isinstance(row, Mapping) or set(row) != set(FIELDS):
            raise ValueError(f"w0[{order}] must contain exactly Uz1 and F1")
        result[order] = {}
        for field in FIELDS:
            packet = row[field]
            if not isinstance(packet, Mapping) or set(packet) != set(BRANCHES):
                raise ValueError(
                    f"w0[{order}][{field}] must contain exactly core and local"
                )
            cells = _cell_orders(packet["local"], f"w0[{order}][{field}].local")
            result[order][field] = {
                "core": _term_list(packet["core"], f"w0[{order}][{field}].core"),
                "local": {
                    cell: _term_list(
                        cells[cell],
                        f"w0[{order}][{field}].local[{cell}]",
                        allow_empty=cell == 3,
                    )
                    for cell in ALL_LOCAL_CELLS
                },
            }
    return result


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
        raise ValueError("caller geometry and cutoff certificates are required")
    required = (
        "hb_strictly_between_zero_and_half",
        "sc_phase_Z_independent",
        "hb_epsilon_Z_independent",
        "original_cutoff_mass_5_over_8",
        "original_cutoff_mass_by_cell",
    )
    missing = [name for name in required if certificate.get(name) is not True]
    if missing:
        raise ValueError("missing caller certificate: " + ", ".join(missing))


def _group_key(term: Mapping[str, Any]) -> tuple[Any, bool, tuple[Any, ...]]:
    info = _origin_info(term["origin"])
    classes = tuple(sorted((_freeze(item) for item in info["error_classes"]), key=repr))
    return tuple(term["powers"]), bool(info["uncertainty"]), classes


def _aggregate(terms: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Outward-add only terms with equal powers and uncertainty classes."""

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
                "path_count": 0,
                "source_count": 0,
                "merged_term_count": 0,
            }
            groups[key] = group
        else:
            # The interval implementation performs directed outward addition.
            group["coefficient"] = group["coefficient"] + term["coefficient"]
        for value in info["source_ids"]:
            group["ids"][_freeze(value)] = value
        for value in info["error_classes"]:
            group["classes"][_freeze(value)] = value
        group["uncertainty"] = group["uncertainty"] or bool(info["uncertainty"])
        group["path_count"] += int(info["path_count"])
        group["source_count"] += int(info["source_count"])
        group["merged_term_count"] += 1

    return [
        {
            "powers": group["powers"],
            "coefficient": group["coefficient"],
            "origin": {
                "operation": "aggregate",
                "source_ids": list(group["ids"].values()),
                "error_classes": list(group["classes"].values()),
                "uncertainty": group["uncertainty"],
                "path_count": group["path_count"],
                "source_count": group["source_count"],
                "merged_term_count": group["merged_term_count"],
                "symbolic_cancellation_claimed": False,
            },
        }
        for group in groups.values()
    ]


def _physical_product(
    leading: Mapping[str, Any],
    word: Mapping[str, Any],
    coefficient_factor: Any,
    *,
    ordinary_order: int,
    field: str,
    branch: str,
    cell: int | None,
    derivative_piece: str,
    external_lambda_power: int,
    radial_bound: str,
    eta_mass: Any,
    smax: Any,
) -> dict[str, Any]:
    powers = tuple(a + b for a, b in zip(leading["powers"], word["powers"]))
    powers = (
        powers[0] + 1,
        powers[1],
        powers[2],
        powers[3] + external_lambda_power,
    )
    return {
        "powers": powers,
        "coefficient": (
            leading["coefficient"]
            * word["coefficient"]
            * coefficient_factor
        ),
        "origin": _new_origin(
            "w0_physical_m15_product",
            leading,
            word,
            ordinary_Z_order=ordinary_order,
            field=field,
            branch=branch,
            phase_cell=cell,
            derivative_piece=derivative_piece,
            external_formal_powers=(1, 0, 0, external_lambda_power),
            eta_mass=eta_mass,
            smax=smax,
            radial_bound=radial_bound,
            symbolic_cancellation_claimed=False,
        ),
    }


def _positive_radial_interval(c: Any, smax: Any, denominator_factor: int) -> Any:
    """Return ``[1, 1/(1-smax/denominator_factor)]``."""

    one = c.mpf(1)
    denominator = one - smax / c.mpf(denominator_factor)
    lower, _ = _scalar_bounds(denominator)
    if lower <= 0:
        raise ValueError("invalid positive radial denominator")
    upper = one / denominator
    return _make_interval(c, one, upper)


def _region_factor(
    c: Any,
    sc: Any,
    *,
    mass: Any,
    smax: Any,
    kind: str,
) -> tuple[Any, Any, str, int]:
    """Build one independent positive radial factor for one output region."""

    if kind == "axial":
        radial = _positive_radial_interval(c, smax, 2)
        factor = c.mpf(2) * sc * mass * radial
        return factor, radial, "[1,1/(1-smax/2)]", -1
    if kind == "swirl":
        radial = _positive_radial_interval(c, smax, 1)
        factor = c.mpf(-8) * sc * mass * radial
        return factor, radial, "[1,1/(1-smax)]", -2
    raise ValueError(f"unknown physical source kind {kind!r}")


def _region_products(
    leading: Mapping[int, list[dict[str, Any]]],
    words: Mapping[int, dict[str, Any]],
    *,
    order: int,
    branch: str,
    cell: int | None,
    coefficient_factor: Any,
    field: str,
    radial_bound: str,
    eta_mass: Any,
    smax: Any,
    external_lambda_power: int,
) -> list[dict[str, Any]]:
    word0 = words[0][field][branch] if branch == "core" else words[0][field]["local"][cell]  # type: ignore[index]
    raw: list[dict[str, Any]] = []
    if order == 0:
        for left in leading[0]:
            for right in word0:
                raw.append(_physical_product(
                    left,
                    right,
                    coefficient_factor,
                    ordinary_order=0,
                    field=field,
                    branch=branch,
                    cell=cell,
                    derivative_piece="lead0*word0",
                    external_lambda_power=external_lambda_power,
                    radial_bound=radial_bound,
                    eta_mass=eta_mass,
                    smax=smax,
                ))
        return raw

    word1 = words[1][field][branch] if branch == "core" else words[1][field]["local"][cell]  # type: ignore[index]
    for left in leading[1]:
        for right in word0:
            raw.append(_physical_product(
                left,
                right,
                coefficient_factor,
                ordinary_order=1,
                field=field,
                branch=branch,
                cell=cell,
                derivative_piece="lead1*word0",
                external_lambda_power=external_lambda_power,
                radial_bound=radial_bound,
                eta_mass=eta_mass,
                smax=smax,
            ))
    for left in leading[0]:
        for right in word1:
            raw.append(_physical_product(
                left,
                right,
                coefficient_factor,
                ordinary_order=1,
                field=field,
                branch=branch,
                cell=cell,
                derivative_piece="lead0*word1",
                external_lambda_power=external_lambda_power,
                radial_bound=radial_bound,
                eta_mass=eta_mass,
                smax=smax,
            ))
    return raw


def _assemble_kind(
    c: Any,
    sc: Any,
    leading: Mapping[int, list[dict[str, Any]]],
    words: Mapping[int, dict[str, Any]],
    *,
    field: str,
    kind: str,
) -> tuple[dict[int, dict[str, Any]], dict[str, Any]]:
    """Assemble core once and active local cells independently."""

    total_mass = c.mpf(5) / c.mpf(8)
    cell_masses = {
        0: c.mpf(1) / c.mpf(4),
        1: c.mpf(1) / c.mpf(4),
        2: c.mpf(1) / c.mpf(8),
        3: c.mpf(0),
    }
    result: dict[int, dict[str, Any]] = {}
    factors: dict[str, Any] = {"core": {}, "local": {}}

    for order in ORDINARY_Z_ORDERS:
        result[order] = {"core": [], "local": {}}

        # The core is a whole-phase sector and gets exactly one total-mass
        # factor, never one factor per q-cell.
        core_smax = sc / c.mpf(4)
        core_factor, core_radial, core_bound, lambda_power = _region_factor(
            c,
            sc,
            mass=total_mass,
            smax=core_smax,
            kind=kind,
        )
        result[order]["core"] = _aggregate(_region_products(
            leading,
            words,
            order=order,
            branch="core",
            cell=None,
            coefficient_factor=core_factor,
            field=field,
            radial_bound=core_bound,
            eta_mass=total_mass,
            smax=core_smax,
            external_lambda_power=lambda_power,
        ))
        factors["core"][order] = {
            "coefficient_factor": core_factor,
            "radial_interval": core_radial,
            "eta_mass": total_mass,
            "smax": core_smax,
            "radial_bound": core_bound,
        }

        for cell in ACTIVE_LOCAL_CELLS:
            mass = cell_masses[cell]
            smax = sc * c.mpf(cell + 1) / c.mpf(16)
            factor, radial, bound, lambda_power = _region_factor(
                c,
                sc,
                mass=mass,
                smax=smax,
                kind=kind,
            )
            result[order]["local"][cell] = _aggregate(_region_products(
                leading,
                words,
                order=order,
                branch="local",
                cell=cell,
                coefficient_factor=factor,
                field=field,
                radial_bound=bound,
                eta_mass=mass,
                smax=smax,
                external_lambda_power=lambda_power,
            ))
            factors["local"].setdefault(cell, {})[order] = {
                "coefficient_factor": factor,
                "radial_interval": radial,
                "eta_mass": mass,
                "smax": smax,
                "radial_bound": bound,
            }

    return result, factors


def build(
    c: Any,
    *,
    sc: Any,
    leading_U: Any,
    leading_F: Any,
    w0: Any,
    source_family: Any,
    source_atoms: Any,
    certificate: Any,
) -> dict[str, Any]:
    """Build physical axial and swirl ``W0`` first-word sectors through Z1.

    ``leading_U[0/1]`` and ``leading_F[0/1]`` are physical whole-phase
    leading lists.  ``w0[order]["Uz1"]`` and ``w0[order]["F1"]`` each carry a
    ``core`` list and four local cell lists.  Z1 is assembled as
    ``lead1*word0 + lead0*word1`` independently for every branch/cell.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must be the caller's MPIntervalContext")
    if source_family is None or (isinstance(source_family, str) and not source_family.strip()):
        raise ValueError("source_family is required")
    _validate_source_atoms(source_atoms)
    _validate_certificate(certificate)

    sc_iv = c.mpf(sc)
    if not _is_interval(sc_iv):
        raise TypeError("sc must be a directed interval")
    sc_lower, sc_upper = _scalar_bounds(sc_iv)
    limit = c.mpf(1) / c.mpf(1024)
    limit_lower, _ = _scalar_bounds(limit)
    if sc_lower <= 0 or sc_upper >= limit_lower:
        raise ValueError("sc must satisfy 0 < sc < 1/1024")
    if sc_lower != sc_upper:
        raise ValueError("sc must be exact (a point interval)")

    leading_u = _normalize_leading(leading_U, "leading_U")
    leading_f = _normalize_leading(leading_F, "leading_F")
    words = _normalize_w0(w0)

    axial, axial_factors = _assemble_kind(
        c,
        sc_iv,
        leading_u,
        words,
        field="Uz1",
        kind="axial",
    )
    swirl, swirl_factors = _assemble_kind(
        c,
        sc_iv,
        leading_f,
        words,
        field="F1",
        kind="swirl",
    )

    eta_cells = [
        c.mpf(1) / c.mpf(4),
        c.mpf(1) / c.mpf(4),
        c.mpf(1) / c.mpf(8),
        c.mpf(0),
    ]
    q_cells = [
        {
            "cell": cell,
            "q_lower": c.mpf(cell) / c.mpf(4),
            "q_upper": c.mpf(cell + 1) / c.mpf(4),
            "s_lower": sc_iv * c.mpf(cell) / c.mpf(16),
            "s_upper": sc_iv * c.mpf(cell + 1) / c.mpf(16),
            "eta_mass": eta_cells[cell],
            "mass": eta_cells[cell],
            "mass_exact": True,
        }
        for cell in ALL_LOCAL_CELLS
    ]

    return {
        "source_family": source_family,
        "source_atoms": source_atoms,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "axial": axial,
        "swirl": swirl,
        "q_cells": q_cells,
        "eta_profile": "1-sigma(4q-2)",
        "eta_cell_masses": eta_cells,
        "cell_masses": eta_cells,
        "eta_integral": c.mpf(5) / c.mpf(8),
        "cutoff_mass": c.mpf(5) / c.mpf(8),
        "eta_cell_masses_exact": True,
        "positive_cutoff_mass_exact": True,
        "core_uses_total_mass_once": True,
        "local_active_cells": list(ACTIVE_LOCAL_CELLS),
        "local_cell_3_omitted": True,
        "independent_region_factors": True,
        "axial_region_factors": axial_factors,
        "swirl_region_factors": swirl_factors,
        "axial_formula": "epsilon*hb*(sc/4)*integral(2*rho*U*W0*eta*dq)",
        "swirl_formula": "-epsilon^2*hb*(sc/4)*integral(2*rho^2*F*W0*eta*dq)",
        "rho_formula": "4*exp(hb*sc*q/4)",
        "formal_axial_external_powers": (1, 0, 0, -1),
        "formal_swirl_external_powers": (1, 0, 0, -2),
        "formal_hb_measure": True,
        "geometry_assumptions_caller_certified": True,
        "geometry_certificate_independently_verified": False,
        "physical_epsilon_scale_axial": "Lambda^-1",
        "physical_epsilon_scale_swirl": "Lambda^-2",
        "leading_U_role": "physical leading whole-phase U",
        "leading_F_role": "physical leading whole-phase F; F_Z connection already included",
        "F_Z_connection_added_here": False,
        "source_ids_union_retained": True,
        "error_classes_union_retained": True,
        "uncertainty_retained": True,
        "aggregation_key": "(powers, uncertainty, error_classes)",
        "aggregation_is_outward_interval_addition_only": True,
        "same_source_cancellation_claimed": False,
        "symbolic_cancellation_claimed": False,
        "allword_tail_added": False,
        "allword_tail_included": False,
        "tail_added": False,
        "tail_included": False,
        "m15_computed": False,
        "complete_m15_computed": False,
        "complete_m15": False,
        "total_m15_computed": False,
        "finite_w0_first_word_only": True,
        "scope": "finite physical W0 first-word m15 sectors; core whole-phase and active local cells remain separate",
    }

