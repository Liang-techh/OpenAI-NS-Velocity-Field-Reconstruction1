"""Cellwise physical first-collar ``m15`` projection of saved Q/P sectors.

This leaf consumes already grouped gauge-exit Q/P packets and applies only the
homogeneous ``K0`` velocity transport needed by the physical first-collar
contribution.  It adds no Q/P feedback, first potential, terminal closure, or
ledger entry.  Formal source powers are ``(hb, F0_anchor, Pstar_squared,
Lambda)`` and are never materialised.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any

SCHEMA = "n1_collar_qp_k0_m15_cells_v1"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
VELOCITY_FIELDS = ("F1", "Uz1", "d_x_F1", "d_x_Uz1")
IGNORED_FIELDS = ("K1", "P1")
ALL_CELLS = (0, 1, 2, 3)
ACTIVE_CELLS = (0, 1, 2)
ETA_CELL_MASS = ((1, 4), (1, 4), (1, 8), (0, 1))

GROUP_NAMES = (
    "Q_E01", "Q_E10", "Q_E11",
    "P_G0_P_W0", "P_G0_P_DeltaW0", "P_D_P_W0", "P_D_P_DeltaW0",
    "Q_inner_R_G0", "Q_inner_R_D",
    "Q_outer_R_Q_W0", "Q_outer_R_Q_Delta", "Q_outer_R_Q_Rg",
    "P_tail_G0_P_Rg", "P_tail_D_P_Rg", "P_tail_R_P_W0",
    "P_tail_R_P_Delta", "P_tail_R_P_Rg",
)
ERROR_GROUPS = frozenset(GROUP_NAMES[7:])
GROUP_PATHS = {
    "Q_E01": ("raw_Q_fixed_exit", "E01", "aggregate"),
    "Q_E10": ("raw_Q_fixed_exit", "E10", "aggregate"),
    "Q_E11": ("raw_Q_fixed_exit", "E11", "aggregate"),
    "P_G0_P_W0": ("saved_potential_compensation_fixed_exit", "G0_P_W0", "exit_aggregate"),
    "P_G0_P_DeltaW0": ("saved_potential_compensation_fixed_exit", "G0_P_DeltaW0", "exit_aggregate"),
    "P_D_P_W0": ("saved_potential_compensation_fixed_exit", "D_P_W0", "exit_aggregate"),
    "P_D_P_DeltaW0": ("saved_potential_compensation_fixed_exit", "D_P_DeltaW0", "exit_aggregate"),
    "Q_inner_R_G0": ("inner_W0_kernel_tail_fixed_exit", "G0", "aggregate"),
    "Q_inner_R_D": ("inner_W0_kernel_tail_fixed_exit", "D", "aggregate"),
    "Q_outer_R_Q_W0": ("outer_kernel_tail_fixed_exit", "R_Q_W0", "Z"),
    "Q_outer_R_Q_Delta": ("outer_kernel_tail_fixed_exit", "R_Q_Delta", "Z"),
    "Q_outer_R_Q_Rg": ("outer_kernel_tail_fixed_exit", "R_Q_Rg", "Z"),
    "P_tail_G0_P_Rg": ("potential_compensation_kernel_tail_fixed_exit", "G0_P_Rg", "Z"),
    "P_tail_D_P_Rg": ("potential_compensation_kernel_tail_fixed_exit", "D_P_Rg", "Z"),
    "P_tail_R_P_W0": ("potential_compensation_kernel_tail_fixed_exit", "R_P_W0", "Z"),
    "P_tail_R_P_Delta": ("potential_compensation_kernel_tail_fixed_exit", "R_P_Delta", "Z"),
    "P_tail_R_P_Rg": ("potential_compensation_kernel_tail_fixed_exit", "R_P_Rg", "Z"),
}

__all__ = [
    "SCHEMA", "SOURCE_POWER_BASIS", "ORDINARY_Z_ORDERS", "FIELDS",
    "VELOCITY_FIELDS", "IGNORED_FIELDS", "ALL_CELLS", "ACTIVE_CELLS",
    "GROUP_NAMES", "ERROR_GROUPS", "build", "build_qp_k0_m15",
]


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _raw_endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    return raw[1 if upper else 0] if raw is not None else value


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    """Construct a directed interval without midpointing."""

    return c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))


def _scalar_bounds(value: Any) -> tuple[Any, Any]:
    raw = getattr(value, "_mpi_", None)
    if raw is None:
        return value, value
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
    """Collect provenance recursively, including top-level term metadata."""

    if _seen is None:
        _seen = set()
    marker = id(value)
    if marker in _seen:
        return {"source_ids": [], "error_classes": [], "uncertainty": False,
                "error_only": False, "symmetric_error": False, "path_count": 1,
                "source_count": 1}
    if isinstance(value, Mapping):
        _seen.add(marker)
        ids = list(_as_items(value.get("source_ids", value.get("sources"))))
        ids.extend(_as_items(value.get("source_id")))
        classes = list(_as_items(value.get("error_classes", value.get("error_class"))))
        operation = str(value.get("operation", "")).lower()
        kind = str(value.get("error_kind", value.get("kind", ""))).lower()
        uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False)))
        error_only = bool(value.get("error_only", value.get("is_error_only", False)))
        symmetric = bool(value.get("symmetric_error", value.get("is_symmetric_error", False)))
        uncertain = uncertain or any(t in operation or t in kind for t in ("error", "uncertainty", "cauchy", "symmetric"))
        error_only = error_only or "error" in operation or "error" in kind
        symmetric = symmetric or "symmetric" in operation or "symmetric" in kind
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
            uncertain = uncertain or info["uncertainty"]
            error_only = error_only or info["error_only"]
            symmetric = symmetric or info["symmetric_error"]
        try:
            paths = max(1, int(value.get("path_count", value.get("aggregate_path_count", 1))))
        except (TypeError, ValueError):
            paths = 1
        try:
            sources = max(1, int(value.get("source_count", value.get("aggregate_source_count", paths))))
        except (TypeError, ValueError):
            sources = paths
        if nested:
            paths = max(paths, max(int(item["path_count"]) for item in nested))
            sources = max(sources, max(int(item["source_count"]) for item in nested))
        return {"source_ids": _unique(ids), "error_classes": _unique(classes),
                "uncertainty": bool(uncertain), "error_only": bool(error_only),
                "symmetric_error": bool(symmetric), "path_count": paths,
                "source_count": sources}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        nested = [_origin_info(item, _seen=_seen) for item in value]
        ids: list[Any] = []
        classes: list[Any] = []
        uncertain = error_only = symmetric = False
        paths = sources = 1
        for info in nested:
            ids.extend(info["source_ids"])
            classes.extend(info["error_classes"])
            uncertain = uncertain or info["uncertainty"]
            error_only = error_only or info["error_only"]
            symmetric = symmetric or info["symmetric_error"]
            paths = max(paths, int(info["path_count"]))
            sources = max(sources, int(info["source_count"]))
        return {"source_ids": _unique(ids), "error_classes": _unique(classes),
                "uncertainty": uncertain, "error_only": error_only,
                "symmetric_error": symmetric, "path_count": paths,
                "source_count": sources}
    return {"source_ids": [], "error_classes": [], "uncertainty": False,
            "error_only": False, "symmetric_error": False, "path_count": 1,
            "source_count": 1}


def _origin_ref(value: Any) -> dict[str, Any]:
    info = _origin_info(value)
    return {key: info[key] for key in ("source_ids", "error_classes", "uncertainty", "error_only", "symmetric_error", "path_count", "source_count")}


def _new_origin(operation: str, *inputs: Any, **metadata: Any) -> dict[str, Any]:
    infos = [_origin_info(item) for item in inputs]
    ids: list[Any] = []
    classes: list[Any] = []
    uncertain = bool(metadata.get("uncertainty", False))
    error_only = bool(metadata.get("error_only", False))
    symmetric = bool(metadata.get("symmetric_error", False))
    if operation.endswith("_product") or operation in {"product", "multiply"}:
        paths = sources = 1
        for info in infos:
            paths *= int(info["path_count"])
            sources *= int(info["source_count"])
    elif operation in {"aggregate", "add", "sum"}:
        paths = sources = 0
        for info in infos:
            paths += int(info["path_count"])
            sources += int(info["source_count"])
    else:
        paths = int(infos[0]["path_count"]) if infos else 1
        sources = int(infos[0]["source_count"]) if infos else 1
    for info in infos:
        ids.extend(info["source_ids"])
        classes.extend(info["error_classes"])
        uncertain = uncertain or info["uncertainty"]
        error_only = error_only or info["error_only"]
        symmetric = symmetric or info["symmetric_error"]
    ids.extend(_as_items(metadata.get("source_ids")))
    classes.extend(_as_items(metadata.get("error_classes")))
    origin = {
        "operation": operation, "source_ids": _unique(ids),
        "error_classes": _unique(classes), "uncertainty": bool(uncertain),
        "error_only": bool(error_only), "symmetric_error": bool(symmetric),
        "path_count": max(1, int(paths)), "source_count": max(1, int(sources)),
        "source_id_count": len(_unique(ids)),
        "symbolic_cancellation_claimed": False,
    }
    if inputs:
        origin["parents"] = tuple(_origin_ref(item) for item in inputs)
    origin.update({key: value for key, value in metadata.items() if key not in {"source_ids", "error_classes", "uncertainty", "error_only", "symmetric_error"}})
    return origin


def _check_powers(value: Any, label: str) -> tuple[int, int, int, int]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)) or len(value) != 4 or any(not isinstance(item, int) or isinstance(item, bool) for item in value):
        raise ValueError(f"{label} must contain four integer formal powers")
    return tuple(int(item) for item in value)  # type: ignore[return-value]


def _check_sector(value: Any, label: str, *, group_error: bool = False) -> dict[str, Any]:
    if not isinstance(value, Mapping) or "powers" not in value or "coefficient" not in value:
        raise ValueError(f"{label} requires powers and coefficient")
    if not _is_interval(value["coefficient"]):
        raise TypeError(f"{label}.coefficient must be a directed interval")
    powers = _check_powers(value["powers"], f"{label}.powers")
    info = _origin_info(value)
    if not info["source_ids"]:
        raise ValueError(f"{label} must retain source_ids")
    origin = value.get("origin") or {"source_ids": list(info["source_ids"]), "error_classes": list(info["error_classes"])}
    result = dict(value)
    result.update({"powers": powers, "coefficient": value["coefficient"], "origin": origin,
                   "source_ids": list(info["source_ids"]), "error_classes": list(info["error_classes"]),
                   "uncertainty": bool(info["uncertainty"] or group_error),
                   "error_only": bool(info["error_only"] or group_error),
                   "symmetric_error": bool(info["symmetric_error"] or group_error)})
    return result


def _term_list(value: Any, label: str, *, allow_empty: bool, group_error: bool) -> list[dict[str, Any]]:
    if isinstance(value, Mapping) and "powers" in value:
        values = [value]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
    else:
        raise TypeError(f"{label} must be a sector sequence")
    if not values and not allow_empty:
        raise ValueError(f"{label} must not be empty")
    return [_check_sector(item, f"{label}[{index}]", group_error=group_error) for index, item in enumerate(values)]


def _order_key(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"boolean is not an ordinary-Z order in {label}")
    if isinstance(value, int):
        return value
    text = str(value).strip().lower()
    if text in {"z0", "0"}:
        return 0
    if text in {"z1", "1"}:
        return 1
    raise ValueError(f"invalid ordinary-Z order {value!r} in {label}")


def _orders(value: Any, label: str) -> dict[int, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must provide explicit Z0 and Z1 rows")
    result: dict[int, Any] = {}
    for raw_key, item in value.items():
        key = _order_key(raw_key, label)
        if key in result:
            raise ValueError(f"duplicate ordinary-Z order {key} in {label}")
        result[key] = item
    if set(result) != set(ORDINARY_Z_ORDERS):
        raise ValueError(f"{label} must provide exactly Z0 and Z1")
    return result


def _normalize_leading(value: Any, label: str) -> dict[int, list[dict[str, Any]]]:
    rows = _orders(value, label)
    return {order: _term_list(rows[order], f"{label}[Z{order}]", allow_empty=False, group_error=False) for order in ORDINARY_Z_ORDERS}


def _normalize_groups(value: Any) -> dict[str, dict[str, dict[int, list[dict[str, Any]]]]]:
    if not isinstance(value, Mapping):
        raise TypeError("inlet_groups must be a mapping")
    for wrapper in ("groups", "inlet_groups", "incoming"):
        if wrapper in value and isinstance(value[wrapper], Mapping):
            value = value[wrapper]
            break
    if set(value) != set(GROUP_NAMES):
        missing = [name for name in GROUP_NAMES if name not in value]
        extra = [name for name in value if name not in GROUP_NAMES]
        raise ValueError(f"inlet_groups must contain exactly 17 named groups; missing={missing}, extra={extra}")
    result: dict[str, dict[str, dict[int, list[dict[str, Any]]]]] = {}
    for group in GROUP_NAMES:
        packet = value[group]
        if not isinstance(packet, Mapping) or set(packet) != set(FIELDS):
            raise ValueError(f"inlet_groups[{group}] must contain exactly {FIELDS}")
        error = group in ERROR_GROUPS
        result[group] = {}
        for field in FIELDS:
            rows = _orders(packet[field], f"inlet_groups[{group}][{field}]")
            result[group][field] = {order: _term_list(rows[order], f"inlet_groups[{group}][{field}][Z{order}]", allow_empty=True, group_error=error) for order in ORDINARY_Z_ORDERS}
    return result


def _validate_source_atoms(source_atoms: Any) -> None:
    if isinstance(source_atoms, Mapping):
        names = set(source_atoms)
    elif isinstance(source_atoms, Sequence) and not isinstance(source_atoms, (str, bytes, bytearray)):
        names = {item.get("name", item.get("atom", item.get("power"))) if isinstance(item, Mapping) else item for item in source_atoms}
    else:
        raise ValueError("source_atoms must name the four formal source atoms")
    if names != set(SOURCE_POWER_BASIS):
        raise ValueError("source_atoms must contain exactly hb, F0_anchor, Pstar_squared, Lambda")


def _certified(certificate: Mapping[str, Any], *names: str) -> bool:
    return any(certificate.get(name) is True for name in names)


def _validate_certificate(certificate: Any) -> None:
    if not isinstance(certificate, Mapping):
        raise ValueError("caller geometry and cutoff certificates are required")
    requirements = (
        ("hb_strictly_between_zero_and_half", "strict_hb_admitted"),
        ("sc_phase_Z_independent", "hb_sc_Z_independent"),
        ("hb_epsilon_Z_independent", "rho_radius_Z_independent"),
        ("original_cutoff_mass_5_over_8", "cutoff_mass_5_over_8_exact", "actual_cutoff_complete_support_measure"),
        ("original_cutoff_mass_by_cell", "cutoff_mass_by_cell_exact", "eta_cell_masses_exact"),
    )
    missing = [names[0] for names in requirements if not _certified(certificate, *names)]
    if missing:
        raise ValueError("missing caller certificate: " + ", ".join(missing))


def _group_key(term: Mapping[str, Any]) -> tuple[Any, bool, bool, tuple[Any, ...], Any]:
    info = _origin_info(term["origin"])
    classes = tuple(sorted((_freeze(item) for item in info["error_classes"]), key=repr))
    origin = term.get("origin")
    branch: Any = None
    if isinstance(origin, Mapping) and isinstance(origin.get("path"), Mapping):
        branch = origin["path"].get("branch")
    elif isinstance(origin, Mapping):
        branch = origin.get("branch")
    return tuple(term["powers"]), bool(info["uncertainty"]), bool(info["error_only"]), classes, branch


def _aggregate(terms: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Outward-add equal formal/error classes within one group/cell."""

    groups: dict[tuple[Any, bool, bool, tuple[Any, ...], Any], dict[str, Any]] = {}
    for term in terms:
        info = _origin_info(term["origin"])
        key = _group_key(term)
        group = groups.get(key)
        if group is None:
            group = {
                "powers": tuple(term["powers"]), "coefficient": term["coefficient"],
                "ids": {}, "classes": {}, "uncertainty": bool(info["uncertainty"]),
                "error_only": bool(info["error_only"]), "symmetric_error": bool(info["symmetric_error"]),
                "path_count": 0, "source_count": 0, "merged_term_count": 0,
                "path_counts": {},
            }
            groups[key] = group
        else:
            group["coefficient"] = group["coefficient"] + term["coefficient"]
        for value in info["source_ids"]:
            group["ids"][_freeze(value)] = value
        for value in info["error_classes"]:
            group["classes"][_freeze(value)] = value
        group["uncertainty"] = group["uncertainty"] or bool(info["uncertainty"])
        group["error_only"] = group["error_only"] or bool(info["error_only"])
        group["symmetric_error"] = group["symmetric_error"] and bool(info["symmetric_error"])
        group["path_count"] += int(info["path_count"])
        group["source_count"] += int(info["source_count"])
        group["merged_term_count"] += 1
        origin = term.get("origin")
        if isinstance(origin, Mapping) and origin.get("paths") is not None:
            # A second aggregation pass receives branch aggregates.  Preserve
            # their compact path multiplicities instead of expanding products.
            for item in origin.get("paths", ()):
                if not isinstance(item, Mapping):
                    continue
                raw_path = item.get("path", ())
                path_key = tuple(raw_path) if isinstance(raw_path, (list, tuple)) else (repr(raw_path),)
                count = int(item.get("count", 1))
                group["path_counts"][path_key] = group["path_counts"].get(path_key, 0) + count
        elif isinstance(origin, Mapping) and origin.get("path") is not None:
            path = origin["path"]
            path_key = tuple(path.get(name) for name in ("group", "channel", "field", "branch", "phase_cell", "ordinary_Z_order", "ordinary_Z_piece")) if isinstance(path, Mapping) else (repr(path),)
            group["path_counts"][path_key] = group["path_counts"].get(path_key, 0) + 1
    result: list[dict[str, Any]] = []
    for group in groups.values():
        path_entries = [{"path": list(path), "count": count} for path, count in group["path_counts"].items()]
        origin = {
            "operation": "aggregate", "source_ids": list(group["ids"].values()),
            "error_classes": list(group["classes"].values()), "uncertainty": group["uncertainty"],
            "error_only": group["error_only"], "symmetric_error": group["symmetric_error"],
            "path_count": group["path_count"], "source_count": group["source_count"],
            "merged_term_count": group["merged_term_count"],
            "paths": path_entries,
            "symbolic_cancellation_claimed": False,
        }
        if len(path_entries) == 1 and len(path_entries[0]["path"]) == 7:
            names = ("group", "channel", "field", "branch", "phase_cell", "ordinary_Z_order", "ordinary_Z_piece")
            origin["path"] = dict(zip(names, path_entries[0]["path"]))
        branch_values = {entry["path"][3] for entry in path_entries if len(entry["path"]) == 7}
        if len(branch_values) == 1:
            origin["branch"] = next(iter(branch_values))
        result.append({
            "powers": group["powers"], "coefficient": group["coefficient"], "origin": origin,
            "source_ids": list(origin["source_ids"]), "error_classes": list(origin["error_classes"]),
            "uncertainty": bool(origin["uncertainty"]), "error_only": bool(origin["error_only"]),
            "symmetric_error": bool(origin["symmetric_error"]),
        })
    return result


def _cell_factors(c: Any, sc: Any, cell: int) -> dict[str, Any]:
    q_lower = c.mpf(cell) / c.mpf(4)
    q_upper = c.mpf(cell + 1) / c.mpf(4)
    q_interval = _make_interval(c, q_lower, q_upper)
    mass_num, mass_den = ETA_CELL_MASS[cell]
    mass = c.mpf(mass_num) / c.mpf(mass_den)

    # 0 <= hb*sc*q/4 <= sc*q_upper/8.  The closure denominator is positive
    # under the caller's sc bound; hb itself is kept formal.
    x_upper = sc * q_upper / c.mpf(8)
    denominator = c.mpf(1) - x_upper
    denominator_lower, _ = _scalar_bounds(denominator)
    if denominator_lower <= 0:
        raise ValueError("invalid certified radial closure denominator")
    rho_candidate = c.mpf(4) / denominator
    rho_squared_candidate = rho_candidate * rho_candidate
    rho_interval = _make_interval(c, c.mpf(4), rho_candidate)
    rho_squared_interval = _make_interval(c, c.mpf(16), rho_squared_candidate)
    # For d = hb*sc*q/4, phi(d)=(1-exp(-d))/d satisfies
    # exp(-d_max) <= phi(d) <= 1.  The caller bound exp(-d_max) >= 1-d_max
    # supplies a directed positive lower endpoint without materialising hb.
    phase_interval = _make_interval(c, denominator, c.mpf(1))

    # The factors match +2*sc*M_i*hb*Lambda^-1*(rho/4) and
    # -8*sc*M_i*hb*Lambda^-2*(rho/4)^2, with hb/Lambda powers carried below.
    axial_Uc = c.mpf(2) * sc * mass * (rho_interval / c.mpf(4))
    axial_Uxc = c.mpf(2) * sc * mass * (sc * q_interval / c.mpf(4)) * (rho_interval / c.mpf(4))
    swirl_Fc = -c.mpf(8) * sc * mass * (rho_squared_interval / c.mpf(16))
    swirl_Fxc = -c.mpf(8) * sc * mass * (sc * q_interval / c.mpf(4)) * phase_interval * (rho_squared_interval / c.mpf(16))
    return {
        "cell": cell, "q_interval": q_interval, "q_lower": q_lower, "q_upper": q_upper,
        "eta_mass": mass, "eta_mass_exact": True,
        "rho_interval": rho_interval, "rho_squared_interval": rho_squared_interval,
        "phase_interval": phase_interval,
        "closure_denominator": denominator,
        "closure": "exp(x) <= 1/(1-x); exp(-d_max) <= (1-exp(-d))/d <= 1; 0 <= hb < 1/2",
        "factors": {"axial_Uc": axial_Uc, "axial_Uxc": axial_Uxc,
                    "swirl_Fc": swirl_Fc, "swirl_Fxc": swirl_Fxc},
    }


def _product(
    leading: Mapping[str, Any], inlet: Mapping[str, Any], factor: Any, *,
    group: str, channel: str, field: str, branch: str, cell: int, order: int,
    derivative_piece: str, external_hb_power: int, external_lambda_power: int,
    factor_formula: str,
) -> dict[str, Any]:
    base = tuple(a + b for a, b in zip(leading["powers"], inlet["powers"]))
    powers = (base[0] + external_hb_power, base[1], base[2], base[3] + external_lambda_power)
    left_info = _origin_info(leading)
    right_info = _origin_info(inlet)
    error_only = group in ERROR_GROUPS or left_info["error_only"] or right_info["error_only"]
    # A symmetric error interval remains symmetric after multiplication by
    # either a signed point/interval factor; do not require the leading sector
    # itself to be tagged as an error.
    symmetric = error_only and (left_info["symmetric_error"] or right_info["symmetric_error"])
    path = {
        "group": group, "channel": channel, "field": field, "branch": branch,
        "phase_cell": cell, "ordinary_Z_order": order, "ordinary_Z_piece": derivative_piece,
        "leading_source_ids": list(left_info["source_ids"]),
        "inlet_source_ids": list(right_info["source_ids"]),
    }
    origin = _new_origin(
        "qp_k0_m15_product", leading, inlet,
        error_classes=["symmetric_error"] if error_only else [],
        uncertainty=bool(left_info["uncertainty"] or right_info["uncertainty"] or error_only),
        error_only=error_only, symmetric_error=symmetric,
        group=group, group_path=GROUP_PATHS[group], channel=channel, field=field,
        branch=branch, phase_cell=cell, ordinary_Z_order=order,
        ordinary_Z_piece=derivative_piece, path=path, factor_formula=factor_formula,
        external_formal_powers=(external_hb_power, 0, 0, external_lambda_power),
        symbolic_cancellation_claimed=False,
    )
    return {
        "powers": powers, "coefficient": leading["coefficient"] * inlet["coefficient"] * factor,
        "origin": origin, "source_ids": list(origin["source_ids"]),
        "error_classes": list(origin["error_classes"]), "uncertainty": bool(origin["uncertainty"]),
        "error_only": bool(origin["error_only"]), "symmetric_error": bool(origin["symmetric_error"]),
    }


def _branch_products(
    leading: Mapping[int, list[dict[str, Any]]],
    packet: Mapping[str, Mapping[int, list[dict[str, Any]]]], *,
    group: str, channel: str, field: str, branch: str, cell: int,
    factor: Any, lambda_power: int, hb_power: int, factor_formula: str,
) -> dict[int, list[dict[str, Any]]]:
    raw: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    for left in leading[0]:
        for right in packet[field][0]:
            raw[0].append(_product(
                left, right, factor, group=group, channel=channel, field=field,
                branch=branch, cell=cell, order=0, derivative_piece="leading*inlet_Z0",
                external_hb_power=hb_power, external_lambda_power=lambda_power,
                factor_formula=factor_formula,
            ))
    for left in leading[1]:
        for right in packet[field][0]:
            raw[1].append(_product(
                left, right, factor, group=group, channel=channel, field=field,
                branch=branch, cell=cell, order=1, derivative_piece="leading_Z*inlet_Z0",
                external_hb_power=hb_power, external_lambda_power=lambda_power,
                factor_formula=factor_formula,
            ))
    for left in leading[0]:
        for right in packet[field][1]:
            raw[1].append(_product(
                left, right, factor, group=group, channel=channel, field=field,
                branch=branch, cell=cell, order=1, derivative_piece="leading*inlet_Z1",
                external_hb_power=hb_power, external_lambda_power=lambda_power,
                factor_formula=factor_formula,
            ))
    return raw


def _channel_result(
    c: Any, sc: Any, leading: Mapping[int, list[dict[str, Any]]],
    packet: Mapping[str, Mapping[int, list[dict[str, Any]]]], *, group: str, channel: str,
) -> tuple[dict[int, dict[int, list[dict[str, Any]]]], dict[int, dict[int, dict[str, list[dict[str, Any]]]]], dict[int, dict[str, Any]]]:
    if channel == "axial":
        branches = (
            ("Uc", "Uz1", 1, "axial_Uc", "+2*sc*M_i*hb*Lambda^-1*(rho/4)"),
            ("Uxc", "d_x_Uz1", 2, "axial_Uxc", "+2*sc*M_i*hb^2*Lambda^-1*(sc*q/4)*(rho/4)"),
        )
        lambda_power = -1
    elif channel == "swirl":
        branches = (
            ("Fc", "F1", 1, "swirl_Fc", "-8*sc*M_i*hb*Lambda^-2*(rho/4)^2"),
            ("Fxc", "d_x_F1", 2, "swirl_Fxc", "-8*sc*M_i*hb^2*Lambda^-2*(sc*q/4)*(rho/4)^2"),
        )
        lambda_power = -2
    else:
        raise ValueError(f"unknown physical channel {channel!r}")
    terms: dict[int, dict[int, list[dict[str, Any]]]] = {order: {cell: [] for cell in ALL_CELLS} for order in ORDINARY_Z_ORDERS}
    paths: dict[int, dict[int, dict[str, list[dict[str, Any]]]]] = {order: {cell: {} for cell in ALL_CELLS} for order in ORDINARY_Z_ORDERS}
    geometries: dict[int, dict[str, Any]] = {}
    for cell in ALL_CELLS:
        geometry = _cell_factors(c, sc, cell)
        geometries[cell] = geometry
        for order in ORDINARY_Z_ORDERS:
            paths[order][cell] = {branch[0]: [] for branch in branches}
        if cell not in ACTIVE_CELLS:
            continue
        combined: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
        for branch, field, hb_power, factor_name, formula in branches:
            raw = _branch_products(
                leading, packet, group=group, channel=channel, field=field, branch=branch,
                cell=cell, factor=geometry["factors"][factor_name], lambda_power=lambda_power,
                hb_power=hb_power, factor_formula=formula,
            )
            # Aggregate within a branch/cell/order to keep the output bounded;
            # aggregate origins retain source IDs and every contributing path.
            for order in ORDINARY_Z_ORDERS:
                branch_terms = _aggregate(raw[order])
                paths[order][cell][branch] = branch_terms
                combined[order].extend(branch_terms)
        for order in ORDINARY_Z_ORDERS:
            terms[order][cell] = _aggregate(combined[order])
    return terms, paths, geometries


def _normalize_certificate_alias(certificate: Any, geometry_certificate: Any) -> Any:
    if certificate is not None and geometry_certificate is not None and certificate is not geometry_certificate:
        raise ValueError("certificate and geometry_certificate cannot both be supplied")
    return geometry_certificate if geometry_certificate is not None else certificate


def build(
    c: Any, *, sc: Any, leading_U: Any, leading_F: Any, inlet_groups: Any = None,
    source_family: Any = None, source_atoms: Any = None, geometry_certificate: Any = None,
    groups: Any = None, incoming: Any = None, certificate: Any = None,
) -> dict[str, Any]:
    """Build the finite physical Q/P first-collar contribution.

    Leading U/F sectors are whole-collar Z0/Z1 lists.  ``inlet_groups`` is the
    exact 17-group mapping emitted by the Q/P inlet stage, with six fields and
    explicit Z0/Z1 lists per field.  Results keep each group independent and
    expose cellwise aggregate sectors plus branch paths.
    """
    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must be the caller's MPIntervalContext")
    if source_family is None or (isinstance(source_family, str) and not source_family.strip()):
        raise ValueError("source_family is required")
    if source_atoms is None:
        raise ValueError("source_atoms is required")
    _validate_source_atoms(source_atoms)
    cert = _normalize_certificate_alias(geometry_certificate, certificate)
    _validate_certificate(cert)

    sc_iv = c.mpf(sc)
    if not _is_interval(sc_iv):
        raise TypeError("sc must be a directed interval")
    sc_lower, sc_upper = _scalar_bounds(sc_iv)
    limit_lower, _ = _scalar_bounds(c.mpf(1) / c.mpf(1024))
    if sc_lower <= 0 or sc_upper >= limit_lower:
        raise ValueError("sc must satisfy 0 < sc < 1/1024")
    if sc_lower != sc_upper:
        raise ValueError("sc must be exact (a point interval)")

    if inlet_groups is not None and (groups is not None or incoming is not None):
        raise ValueError("inlet_groups cannot be mixed with groups or incoming")
    raw_groups = inlet_groups if inlet_groups is not None else groups if groups is not None else incoming
    if raw_groups is None:
        raise ValueError("inlet_groups is required")
    leading_u = _normalize_leading(leading_U, "leading_U")
    leading_f = _normalize_leading(leading_F, "leading_F")
    normalized_groups = _normalize_groups(raw_groups)

    cell_geometry = {cell: _cell_factors(c, sc_iv, cell) for cell in ALL_CELLS}
    group_results: dict[str, Any] = {}
    for group in GROUP_NAMES:
        packet = normalized_groups[group]
        axial, axial_paths, axial_geometry = _channel_result(c, sc_iv, leading_u, packet, group=group, channel="axial")
        swirl, swirl_paths, swirl_geometry = _channel_result(c, sc_iv, leading_f, packet, group=group, channel="swirl")
        group_results[group] = {
            "group": group, "group_path": GROUP_PATHS[group],
            "group_kind": "error_only" if group in ERROR_GROUPS else "finite_QP",
            "error_only_group": group in ERROR_GROUPS,
            "input_error_only": group in ERROR_GROUPS,
            "symmetric_error_input": group in ERROR_GROUPS,
            "axial": axial, "swirl": swirl,
            "paths": {"axial": axial_paths, "swirl": swirl_paths},
            "input_fields": list(FIELDS),
            "input_sector_counts": {field: {f"Z{order}": len(packet[field][order]) for order in ORDINARY_Z_ORDERS} for field in FIELDS},
            "input_sector_paths_retained_in_runner_index": True,
            "velocity_fields_used": {"axial": ["Uz1", "d_x_Uz1"], "swirl": ["F1", "d_x_F1"]},
            "ignored_fields": list(IGNORED_FIELDS), "P1_is_not_m15": True, "K1_is_not_m15": True,
            "ordinary_Z1_product_rule": "leading_Z*inlet_Z0 + leading*inlet_Z1",
            "Z1_product_rule": "leading_Z*inlet_Z0 + leading*inlet_Z1",
            "source_ids_retained": True, "error_provenance_retained": True,
        }

    eta_masses = [c.mpf(num) / c.mpf(den) for num, den in ETA_CELL_MASS]
    q_cells = [{
        "cell": cell, "q_lower": cell_geometry[cell]["q_lower"], "q_upper": cell_geometry[cell]["q_upper"],
        "q_interval": cell_geometry[cell]["q_interval"], "eta_mass": eta_masses[cell],
        "eta_mass_exact": True, "velocity_cutoff_active": cell in ACTIVE_CELLS,
    } for cell in ALL_CELLS]
    return {
        "schema": SCHEMA, "source_family": source_family, "source_atoms": source_atoms,
        "source_power_basis": list(SOURCE_POWER_BASIS), "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "fields": list(FIELDS), "velocity_fields": list(VELOCITY_FIELDS), "ignored_fields": list(IGNORED_FIELDS),
        "inlet_groups": list(GROUP_NAMES), "group_results": group_results,
        "group_paths": {group: list(GROUP_PATHS[group]) for group in GROUP_NAMES}, "q_cells": q_cells,
        "eta_profile": "1-sigma(4q-2)", "eta_cell_masses": eta_masses, "eta_cell_masses_exact": True,
        "eta_integral": c.mpf(5) / c.mpf(8), "cutoff_mass": c.mpf(5) / c.mpf(8), "cell_masses": eta_masses,
        "cutoff_mass_by_cell": eta_masses, "positive_cutoff_mass_exact": True, "cutoff_mass_exact": True,
        "active_cells": list(ACTIVE_CELLS), "zero_mass_cell": 3, "zero_mass_applies_to_velocity_only": True,
        "radial_phase_enclosures": cell_geometry,
        "rho_formula": "4*exp(hb*sc*q/4)", "d_rho_formula": "hb*(sc/4)*rho*dq",
        "H0U_formula": "Uc + hb*sc*q/4*Uxc", "H0F_formula": "Fc + (1-exp(-hb*sc*q/4))*Fxc",
        "axial_formula": "+2*Lambda^-1*Ulead*H0U*d_rho",
        "swirl_formula": "-2*Lambda^-2*rho*Flead*H0F*d_rho",
        "ordinary_Z1_product_rule": "leading_Z*inlet_Z0 + leading*inlet_Z1",
        "Z1_product_rule": "leading_Z*inlet_Z0 + leading*inlet_Z1",
        "formal_external_powers": {"axial_Uc": (1, 0, 0, -1), "axial_Uxc": (2, 0, 0, -1), "swirl_Fc": (1, 0, 0, -2), "swirl_Fxc": (2, 0, 0, -2)},
        "formal_source_powers_not_materialized": True, "formal_hb_not_materialized": True,
        "formal_hb_measure": True, "physical_epsilon_scale_axial": "Lambda^-1", "physical_epsilon_scale_swirl": "Lambda^-2",
        "geometry_certificate_caller_owned": True, "geometry_certificate_independently_verified": False,
        "certified_radial_closure": True,
        "radial_closure_inequality": "exp(x) <= 1/(1-x), x <= sc*q_upper/8",
        "H0F_phase_closure": "exp(-d_max) <= (1-exp(-d))/d <= 1, d=hb*sc*q/4",
        "signed_axial_coefficients": True, "signed_swirl_coefficients": True,
        "source_ids_union_retained": True, "error_classes_union_retained": True, "uncertainty_retained": True,
        "independent_group_results": True, "aggregation_key": "(group, cell, channel, Z order, inlet branch, powers, uncertainty, error_classes)",
        "aggregation_is_outward_interval_addition_only": True, "P1_conflated_with_m15": False,
        "new_QP_feedback_added": False, "first_potential_added": False, "terminal_closure_proved": False,
        "appended_to_212_ledger": False, "appended_to_ledger": False,
        "allword_tail_added": False, "allword_tail_included": False,
        "finite_physical_qp_k0_m15_computed": True,
        "physical_m15_finite_contribution_computed": True,
        "m15_computed": False,
        "complete_m15_computed": False, "complete_m15": False,
        "scope": "standalone finite physical first-collar Q/P K0 source-sector projection; no feedback, first potential, or terminal closure",
    }


build_qp_k0_m15 = build


