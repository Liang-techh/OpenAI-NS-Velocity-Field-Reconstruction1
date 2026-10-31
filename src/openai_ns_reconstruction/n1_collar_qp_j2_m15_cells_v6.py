"""Cellwise physical ``m15`` enclosure for the Q/P Kchi ``j>=2`` tail.

This leaf consumes the already grouped raw ``Fc/Fxc`` Q/P inlet sectors and
returns only symmetric error radii for the omitted homogeneous potential
words.  It does not re-enter K0, the first potential word, local forcing, or
any feedback word.  The physical output is the swirl velocity m15 error
projection of the tail; no signed m15 field or ledger term is produced.

The formal source basis is ``(hb, F0_anchor, Pstar_squared, Lambda)``.  The
whole-collar majorant keeps ``hb`` formal.  Its F branches add ``+4`` and
``+5`` formal hb powers for ``Fc`` and ``Fxc`` respectively.  The physical
``d-rho`` factor adds one further hb and ``Lambda^-2``.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any


SCHEMA = "n1_collar_qp_j2_m15_cells_v6"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
INLET_COMPONENTS = ("Fc", "Fxc")
BRANCHES = ("F_from_Fc", "F_from_Fxc")
ALL_CELLS = (0, 1, 2, 3)
ACTIVE_CELLS = (0, 1, 2)
ETA_CELL_MASS = ((1, 4), (1, 4), (1, 8), (0, 1))
FORMAL_HB_POWER_ADDED = {"F_from_Fc": 4, "F_from_Fxc": 5}
FORMAL_M15_HB_ADDED = {branch: power + 1 for branch, power in FORMAL_HB_POWER_ADDED.items()}

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
    "SCHEMA", "SOURCE_POWER_BASIS", "ORDINARY_Z_ORDERS", "INLET_COMPONENTS",
    "BRANCHES", "GROUP_NAMES", "ERROR_GROUPS", "FORMAL_HB_POWER_ADDED",
    "FORMAL_M15_HB_ADDED", "build", "build_qp_j2_m15", "build_tail", "evaluate",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _items(value: Any) -> list[Any]:
    if value is None:
        return []
    return list(value) if _is_sequence(value) else [value]


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


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _raw(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    raw = getattr(value, "_mpf_", None)
    return raw if raw is not None else value


def _make(c: Any, lower: Any, upper: Any) -> Any:
    return c.make_mpf((_raw(lower, False), _raw(upper, True)))


def _iv(c: Any, value: Any, label: str = "interval") -> Any:
    if _is_interval(value):
        return _make(c, value, value)
    if isinstance(value, Mapping):
        for key in ("value", "interval", "box", "enclosure", "data"):
            if key in value:
                return _iv(c, value[key], label)
        for left, right in (("lower", "upper"), ("lo", "hi"), ("left", "right")):
            if left in value and right in value:
                return _make(c, value[left], value[right])
        raise TypeError(f"{label} requires an interval-compatible value")
    if _is_sequence(value) and len(value) == 2:
        return _make(c, value[0], value[1])
    try:
        return c.mpf(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{label} is not interval-compatible") from exc


def _bounds(c: Any, value: Any) -> tuple[Any, Any]:
    value = _iv(c, value)
    raw = value._mpi_
    return c.make_mpf((raw[0], raw[0])), c.make_mpf((raw[1], raw[1]))


def _lo(c: Any, value: Any) -> Any:
    return _bounds(c, value)[0]


def _hi(c: Any, value: Any) -> Any:
    return _bounds(c, value)[1]


def _zero(c: Any) -> Any:
    return c.mpf(0)


def _abs_upper(c: Any, value: Any) -> Any:
    lo, hi = _bounds(c, value)
    if lo >= 0:
        return hi
    if hi <= 0:
        return -lo
    return max(-lo, hi)


def _symmetric(c: Any, radius: Any) -> Any:
    radius = _hi(c, radius)
    if radius < 0:
        raise ValueError("symmetric radius must be nonnegative")
    return _make(c, -radius, radius)


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _iv(c, value, label)
    lo, hi = _bounds(c, value)
    if hi < 0:
        raise ValueError(f"{label} is negative")
    if lo < 0:
        lo = _zero(c)
    return _make(c, lo, hi)


def _positive(c: Any, value: Any, label: str) -> Any:
    value = _iv(c, value, label)
    if _lo(c, value) <= 0:
        raise ValueError(f"{label} must be strictly positive")
    return value


def _origin_info(value: Any, *, seen: set[int] | None = None) -> dict[str, Any]:
    result = {
        "source_ids": [], "error_classes": [], "uncertainty": False,
        "error_only": False, "symmetric_error": False, "path_count": 1,
        "source_count": 1,
    }
    if seen is None:
        seen = set()
    marker = id(value)
    if marker in seen:
        return result
    if isinstance(value, Mapping):
        seen.add(marker)
        ids: list[Any] = []
        classes: list[Any] = []
        for key in ("source_ids", "sources", "source_id", "source", "input_id"):
            ids.extend(_items(value.get(key)))
        for key in ("error_classes", "error_class", "uncertainty_class"):
            classes.extend(_items(value.get(key)))
        operation = str(value.get("operation", "")).lower()
        kind = str(value.get("kind", value.get("error_kind", ""))).lower()
        uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False)))
        error_only = bool(value.get("error_only", value.get("is_error_only", False)))
        symmetric = bool(value.get("symmetric_error", value.get("is_symmetric_error", False)))
        uncertain = uncertain or any(token in operation or token in kind for token in ("error", "uncertainty", "symmetric"))
        error_only = error_only or "error" in operation or "error" in kind
        symmetric = symmetric or "symmetric" in operation or "symmetric" in kind
        nested: list[dict[str, Any]] = []
        for key in ("origin", "metadata", "parents"):
            child = value.get(key)
            if child is None or child is value:
                continue
            for item in (_items(child) if _is_sequence(child) else [child]):
                nested.append(_origin_info(item, seen=seen))
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
        return {
            "source_ids": _unique(ids), "error_classes": _unique(classes),
            "uncertainty": bool(uncertain), "error_only": bool(error_only),
            "symmetric_error": bool(symmetric), "path_count": paths,
            "source_count": sources,
        }
    if _is_sequence(value):
        nested = [_origin_info(item, seen=seen) for item in value]
        ids = [item for info in nested for item in info["source_ids"]]
        classes = [item for info in nested for item in info["error_classes"]]
        return {
            "source_ids": _unique(ids), "error_classes": _unique(classes),
            "uncertainty": any(info["uncertainty"] for info in nested),
            "error_only": any(info["error_only"] for info in nested),
            "symmetric_error": any(info["symmetric_error"] for info in nested),
            "path_count": max([1] + [int(info["path_count"]) for info in nested]),
            "source_count": max([1] + [int(info["source_count"]) for info in nested]),
        }
    return result


def _origin_ref(value: Any) -> dict[str, Any]:
    info = _origin_info(value)
    return {key: info[key] for key in ("source_ids", "error_classes", "uncertainty", "error_only", "symmetric_error", "path_count", "source_count")}


def _new_origin(operation: str, *inputs: Any, **metadata: Any) -> dict[str, Any]:
    infos = [_origin_info(item) for item in inputs]
    ids = [item for info in infos for item in info["source_ids"]]
    classes = [item for info in infos for item in info["error_classes"]]
    if operation.endswith("_product") or operation in {"product", "multiply"}:
        paths = sources = 1
        for info in infos:
            paths *= int(info["path_count"])
            sources *= int(info["source_count"])
    elif operation in {"aggregate", "add", "sum"}:
        paths = sum(int(info["path_count"]) for info in infos)
        sources = sum(int(info["source_count"]) for info in infos)
    else:
        paths = int(infos[0]["path_count"]) if infos else 1
        sources = int(infos[0]["source_count"]) if infos else 1
    ids.extend(_items(metadata.get("source_ids")))
    classes.extend(_items(metadata.get("error_classes")))
    origin = {
        "operation": operation, "source_ids": _unique(ids),
        "error_classes": _unique(classes), "uncertainty": True,
        "error_only": True, "symmetric_error": True,
        "path_count": max(1, int(paths)), "source_count": max(1, int(sources)),
        "source_id_count": len(_unique(ids)),
        "symbolic_cancellation_claimed": False,
    }
    if inputs:
        origin["parents"] = tuple(_origin_ref(item) for item in inputs)
    origin.update({key: value for key, value in metadata.items() if key not in {"source_ids", "error_classes"}})
    return origin


def _order_key(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"boolean is not an ordinary-Z order in {label}")
    if isinstance(value, int):
        return value
    text = str(value).strip().lower()
    if text in {"0", "z0", "value", "ordinary_z0"}:
        return 0
    if text in {"1", "z1", "derivative", "ordinary_z1", "ordinary_z"}:
        return 1
    raise ValueError(f"invalid ordinary-Z order {value!r} in {label}")


def _term_like(value: Any) -> bool:
    return isinstance(value, Mapping) and any(key in value for key in ("coefficient", "coefficient_interval", "value", "interval", "powers"))


def _term(c: Any, raw: Any, label: str, default_id: str, *, group_error: bool = False) -> dict[str, Any]:
    item = dict(raw) if isinstance(raw, Mapping) else {"coefficient": raw}
    coefficient = next((item[key] for key in ("coefficient", "coefficient_interval", "value", "interval") if key in item), None)
    if coefficient is None:
        raise ValueError(f"{label} is missing coefficient")
    coefficient = _iv(c, coefficient, label + ".coefficient")
    info = _origin_info(item)
    ids = _unique(info["source_ids"])
    if not ids:
        ids = [default_id]
    classes = _unique(info["error_classes"])
    error_only = bool(info["error_only"] or group_error)
    symmetric = bool(info["symmetric_error"] or group_error)
    if error_only and "symmetric_error" not in classes:
        classes.append("symmetric_error")
    powers = item.get("powers")
    if powers is None:
        powers = (
            item.get("hb_power", 0), item.get("F0_power", item.get("F0_anchor_power", 0)),
            item.get("Pstar_squared_power", item.get("Pstar_power", 0)),
            item.get("Lambda_power", item.get("lambda_power", 0)),
        )
    if not _is_sequence(powers) or len(powers) != 4 or any(isinstance(power, bool) or not isinstance(power, int) for power in powers):
        raise ValueError(f"{label}.powers must contain four integer formal powers")
    origin_raw = item.get("origin")
    if isinstance(origin_raw, Mapping):
        origin = dict(origin_raw)
        origin.setdefault("source_ids", ids)
        origin.setdefault("error_classes", classes)
        origin.setdefault("uncertainty", bool(info["uncertainty"] or error_only))
    else:
        origin = {
            "source_ids": ids, "error_classes": classes,
            "uncertainty": bool(info["uncertainty"] or error_only),
            "error_only": error_only, "symmetric_error": symmetric,
        }
    result = dict(item)
    result.update({
        "powers": tuple(int(power) for power in powers), "coefficient": coefficient,
        "source_ids": ids, "error_classes": classes,
        "uncertainty": bool(info["uncertainty"] or error_only), "error_only": error_only,
        "symmetric_error": symmetric, "origin": origin,
    })
    if isinstance(origin_raw, str):
        result["origin_text"] = origin_raw
    return result


def _term_list(c: Any, value: Any, label: str, *, allow_empty: bool, group_error: bool) -> list[dict[str, Any]]:
    if isinstance(value, Mapping) and _term_like(value):
        values = [value]
    elif _is_sequence(value):
        values = list(value)
    else:
        values = [value]
    if not values and not allow_empty:
        raise ValueError(f"{label} must not be empty")
    return [_term(c, item, f"{label}[{index}]", f"{label}:{index}", group_error=group_error) for index, item in enumerate(values)]


def _rows(c: Any, value: Any, label: str, *, allow_empty: bool, group_error: bool) -> dict[int, list[dict[str, Any]]]:
    raw = value
    if isinstance(raw, Mapping):
        for wrapper in ("ordinary_Z", "source_Z", "source_jets", "jets", "rows"):
            if wrapper in raw and isinstance(raw[wrapper], (Mapping, Sequence)):
                raw = raw[wrapper]
                break
    if _is_sequence(raw) and len(raw) == 2 and not all(_term_like(item) for item in raw):
        return {order: _term_list(c, raw[order], f"{label}.Z{order}", allow_empty=allow_empty, group_error=group_error) for order in ORDINARY_Z_ORDERS}
    if not isinstance(raw, Mapping):
        raise ValueError(f"{label} requires explicit ordinary Z0 and Z1 rows")
    found: dict[int, Any] = {}
    for key, item in raw.items():
        order = _order_key(key, label)
        if order in found:
            raise ValueError(f"{label} supplies Z{order} twice")
        found[order] = item
    if set(found) != set(ORDINARY_Z_ORDERS):
        raise ValueError(f"{label} requires explicit ordinary Z0 and Z1 rows")
    return {order: _term_list(c, found[order], f"{label}.Z{order}", allow_empty=allow_empty, group_error=group_error) for order in ORDINARY_Z_ORDERS}


def _lookup(mapping: Mapping[str, Any], names: Sequence[str]) -> Any:
    for name in names:
        if name in mapping:
            return mapping[name]
    lowered = {str(key).lower(): value for key, value in mapping.items()}
    for name in names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    return None


def _normalize_packet(c: Any, packet: Any, label: str, *, group_error: bool) -> dict[str, dict[int, list[dict[str, Any]]]]:
    if not isinstance(packet, Mapping):
        raise TypeError(f"{label} must map Fc and Fxc")
    fc = _lookup(packet, ("Fc", "F_c", "incoming_Fc", "F"))
    fxc = _lookup(packet, ("Fxc", "F_xc", "Fx_c", "incoming_Fxc", "Fx"))
    if fc is None or fxc is None:
        raise ValueError(f"{label} requires both Fc and Fxc")
    return {
        "Fc": _rows(c, fc, label + ".Fc", allow_empty=False, group_error=group_error),
        "Fxc": _rows(c, fxc, label + ".Fxc", allow_empty=False, group_error=group_error),
    }


def _normalize_groups(c: Any, value: Any) -> dict[str, dict[str, dict[int, list[dict[str, Any]]]]]:
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
    return {
        group: _normalize_packet(c, value[group], f"inlet_groups[{group}]", group_error=group in ERROR_GROUPS)
        for group in GROUP_NAMES
    }


def _normalize_leading(c: Any, value: Any) -> dict[int, list[dict[str, Any]]]:
    rows = _rows(c, value, "leading_F", allow_empty=False, group_error=False)
    return rows


def _normalize_chi(c: Any, value: Any) -> dict[int, list[dict[str, Any]]]:
    if not isinstance(value, Mapping):
        raise TypeError("chi must provide ordinary Z0 and Z1 rows")
    raw = value
    for wrapper in ("ordinary_Z", "chi_jets", "jets", "rows"):
        if wrapper in raw and isinstance(raw[wrapper], Mapping):
            raw = raw[wrapper]
            break
    rows = _rows(c, raw, "chi", allow_empty=False, group_error=False)
    for order in ORDINARY_Z_ORDERS:
        if len(rows[order]) != 1:
            raise ValueError(f"chi Z{order} requires exactly one source-bound sector")
        if tuple(rows[order][0]["powers"]) != (0, 0, 0, 0):
            raise ValueError("chi formal powers must be zero when folded into the nonlinear majorant")
        if not rows[order][0]["source_ids"]:
            raise ValueError(f"chi Z{order} must retain a source ID")
    return rows


def _validate_atoms(source_atoms: Any) -> None:
    if isinstance(source_atoms, Mapping):
        names = set(source_atoms)
    elif _is_sequence(source_atoms):
        names = {
            item.get("name", item.get("atom", item.get("power"))) if isinstance(item, Mapping) else item
            for item in source_atoms
        }
    else:
        raise ValueError("source_atoms must name the four formal source atoms")
    if names != set(SOURCE_POWER_BASIS):
        raise ValueError("source_atoms must contain exactly hb, F0_anchor, Pstar_squared, Lambda")


def _validate_certificate(certificate: Any) -> None:
    if not isinstance(certificate, Mapping):
        raise ValueError("caller geometry and cutoff certificates are required")
    evidence = certificate.get("geometry", certificate)
    if not isinstance(evidence, Mapping):
        raise ValueError("geometry certificate must be a mapping")
    aliases = (
        ("hb_strictly_between_zero_and_half", "strict_hb_admitted"),
        ("sc_phase_Z_independent", "hb_sc_Z_independent"),
        ("hb_epsilon_Z_independent", "rho_radius_Z_independent"),
        ("original_cutoff_mass_by_cell", "cutoff_mass_by_cell_exact", "eta_cell_masses_exact"),
    )
    missing = [
        names[0]
        for names in aliases
        if not any(certificate.get(name) is True or evidence.get(name) is True for name in names)
    ]
    if missing:
        raise ValueError("missing caller certificate: " + ", ".join(missing))


def _kernel_bounds(c: Any, sc: Any, chi: Mapping[int, Sequence[Mapping[str, Any]]]) -> dict[str, Any]:
    phase = sc / c.mpf(8)
    A = sc / c.mpf(16)
    E = _nonnegative(c, c.exp(A), "E")
    C0 = _nonnegative(c, _abs_upper(c, chi[0][0]["coefficient"]), "C0")
    C1 = _nonnegative(c, _abs_upper(c, chi[1][0]["coefficient"]), "C1")
    k0 = _nonnegative(c, c.mpf(8) * phase**2 * E**2 * C0, "k0")
    kz = _nonnegative(c, c.mpf(8) * phase**2 * E**2 * C1, "kz")
    T = _nonnegative(c, c.exp(k0 / c.mpf(4)), "T")
    return {
        "c": phase, "A": A, "E": E, "C0": C0, "C1": C1,
        "k0": k0, "kz": kz, "T": T,
        "formulas": {
            "c": "sc/8", "A": "sc/16", "E": "exp(A)",
            "C0": "sup|chi|", "C1": "sup|chi_Z|",
            "k0": "8*c**2*E**2*C0", "kz": "8*c**2*E**2*C1",
            "T": "exp(k0/4)",
        },
    }


def _cell_factors(c: Any, sc: Any, cell: int) -> dict[str, Any]:
    q_lower = c.mpf(cell) / c.mpf(4)
    q_upper = c.mpf(cell + 1) / c.mpf(4)
    mass_num, mass_den = ETA_CELL_MASS[cell]
    mass = c.mpf(mass_num) / c.mpf(mass_den)
    denominator = c.mpf(1) - sc * q_upper / c.mpf(8)
    if _lo(c, denominator) <= 0:
        raise ValueError("invalid positive rho closure denominator")
    rho_upper = c.mpf(4) / denominator
    rho_squared = _make(c, c.mpf(16), rho_upper**2)
    physical_factor = c.mpf(8) * sc * mass * (rho_squared / c.mpf(16))
    return {
        "cell": cell, "q_lower": q_lower, "q_upper": q_upper,
        "q_interval": _make(c, q_lower, q_upper), "eta_mass": mass,
        "eta_mass_exact": True, "rho_squared_interval": rho_squared,
        "closure_denominator": denominator,
        "physical_factor": physical_factor,
        "physical_factor_formula": "8*sc*M_i*(rho^2/16)",
        "rho_formula": "4*exp(hb*sc*q/4)",
        "d_rho_formula": "hb*(sc/4)*rho*dq",
        "velocity_cutoff_active": cell in ACTIVE_CELLS,
    }


def _group_key(term: Mapping[str, Any]) -> tuple[Any, ...]:
    origin = term["origin"]
    classes = tuple(sorted((_freeze(item) for item in origin.get("error_classes", ())), key=repr))
    path = origin.get("path", {}) if isinstance(origin, Mapping) else {}
    return tuple(term["powers"]), classes, path.get("ordinary_Z_piece") if isinstance(path, Mapping) else None


def _aggregate(terms: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[Any, ...], dict[str, Any]] = {}
    for term in terms:
        info = _origin_info(term["origin"])
        key = _group_key(term)
        group = groups.get(key)
        if group is None:
            group = {
                "powers": tuple(term["powers"]), "coefficient": term["coefficient"],
                "ids": {}, "classes": {}, "paths": {}, "merged_term_count": 0,
                "path_count": 0, "source_count": 0,
                "incoming_ids": {}, "leading_ids": {}, "chi_ids": {},
                "chi_classes": {},
            }
            groups[key] = group
        else:
            group["coefficient"] = group["coefficient"] + term["coefficient"]
        for value in info["source_ids"]:
            group["ids"][_freeze(value)] = value
        for value in info["error_classes"]:
            group["classes"][_freeze(value)] = value
        for value in term.get("incoming_source_ids", ()):
            group["incoming_ids"][_freeze(value)] = value
        for value in term.get("leading_source_ids", ()):
            group["leading_ids"][_freeze(value)] = value
        for value in term.get("chi_source_ids", ()):
            group["chi_ids"][_freeze(value)] = value
        for value in term.get("chi_source_classes", ()):
            group["chi_classes"][_freeze(value)] = value
        group["merged_term_count"] += 1
        group["path_count"] += int(info["path_count"])
        group["source_count"] += int(info["source_count"])
        path = term.get("origin", {}).get("path", {})
        path_key = repr(path)
        entry = group["paths"].setdefault(path_key, {"path": dict(path), "count": 0})
        entry["count"] += 1
    result: list[dict[str, Any]] = []
    for group in groups.values():
        paths = list(group["paths"].values())
        origin = {
            "operation": "aggregate", "source_ids": list(group["ids"].values()),
            "error_classes": list(group["classes"].values()), "uncertainty": True,
            "error_only": True, "symmetric_error": True,
            "path_count": group["path_count"], "source_count": group["source_count"],
            "merged_term_count": group["merged_term_count"], "paths": paths,
            "incoming_source_ids": list(group["incoming_ids"].values()),
            "leading_source_ids": list(group["leading_ids"].values()),
            "chi_source_ids": list(group["chi_ids"].values()),
            "chi_source_classes": list(group["chi_classes"].values()),
            "symbolic_cancellation_claimed": False,
        }
        if len(paths) == 1:
            origin["path"] = paths[0]["path"]
        item = {
            "powers": group["powers"], "coefficient": group["coefficient"],
            "origin": origin, "source_ids": list(origin["source_ids"]),
            "error_classes": list(origin["error_classes"]), "uncertainty": True,
            "error_only": True, "symmetric_error": True,
            "incoming_source_ids": list(group["incoming_ids"].values()),
            "leading_source_ids": list(group["leading_ids"].values()),
            "chi_source_ids": list(group["chi_ids"].values()),
            "chi_source_classes": list(group["chi_classes"].values()),
        }
        if paths:
            path = paths[0]["path"]
            item.update({
                "group": path.get("group"), "branch": path.get("branch"),
                "phase_cell": path.get("phase_cell"),
                "ordinary_Z_order": path.get("ordinary_Z_order"),
                "ordinary_Z_piece": path.get("ordinary_Z_piece"),
                "tail_starts_at_j": 2,
                "tail_formal_hb_power_added": FORMAL_HB_POWER_ADDED.get(path.get("branch")),
                "formal_hb_power_added": FORMAL_M15_HB_ADDED.get(path.get("branch")),
                "physical_hb_power_added": 1,
                "formal_external_powers": (1, 0, 0, -2),
                "coefficient_is_symmetric_radius": True,
            })
        result.append(item)
    return result


def _tail_term(
    c: Any, *, leading: Mapping[str, Any], incoming: Mapping[str, Any],
    chi_terms: Sequence[Mapping[str, Any]], factor: Any, group: str,
    branch: str, cell: int, order: int, piece: str,
    leading_index: int, incoming_index: int, chi_index: int,
    chi_z_index: int | None, tail_hb_power: int,
) -> dict[str, Any]:
    leading_info = _origin_info(leading)
    incoming_info = _origin_info(incoming)
    chi_ids = _unique(source_id for item in chi_terms for source_id in item.get("source_ids", []))
    chi_classes = _unique(error_class for item in chi_terms for error_class in item.get("error_classes", []))
    ids = _unique(list(leading_info["source_ids"]) + list(incoming_info["source_ids"]) + chi_ids)
    classes = _unique(list(leading_info["error_classes"]) + list(incoming_info["error_classes"]) + chi_classes + ["j2_kernel_remainder", "symmetric_error"])
    radius = (
        _abs_upper(c, leading["coefficient"])
        * _abs_upper(c, incoming["coefficient"])
        * factor
    )
    powers = tuple(a + b for a, b in zip(leading["powers"], incoming["powers"]))
    powers = (powers[0] + tail_hb_power + 1, powers[1], powers[2], powers[3] - 2)
    path = {
        "group": group, "group_path": list(GROUP_PATHS[group]), "branch": branch,
        "phase_cell": cell, "ordinary_Z_order": order, "ordinary_Z_piece": piece,
        "leading_index": leading_index, "incoming_index": incoming_index,
        "chi_index": chi_index, "chi_Z_index": chi_z_index,
    }
    origin = _new_origin(
        "qp_j2_m15_symmetric_error_product", leading, incoming, *chi_terms,
        source_ids=ids, error_classes=classes, group=group,
        group_path=list(GROUP_PATHS[group]), branch=branch, field="F",
        incoming_component="Fxc" if branch.endswith("Fxc") else "Fc",
        phase_cell=cell, ordinary_Z_order=order, ordinary_Z_piece=piece,
        leading_index=leading_index, incoming_index=incoming_index,
        chi_index=chi_index, chi_Z_index=chi_z_index, path=path,
        chi_source_ids=chi_ids, chi_source_classes=chi_classes,
        leading_source_ids=list(leading_info["source_ids"]),
        incoming_source_ids=list(incoming_info["source_ids"]),
        tail_formal_hb_power_added=tail_hb_power,
        formal_hb_power_added=tail_hb_power + 1, physical_hb_power_added=1,
        formal_external_powers=(1, 0, 0, -2),
        coefficient_is_symmetric_radius=True, coefficient_contains_formal_hb=False,
        signed_source_term=False, signed_leading_term=False,
        tail_starts_at_j=2, kernel_tail_order="j>=2",
        first_potential_reexecuted=False, feedback_tail_included=False,
        symbolic_cancellation_claimed=False,
    )
    return {
        "powers": powers, "coefficient": _symmetric(c, radius), "origin": origin,
        "source_ids": ids, "error_classes": classes, "uncertainty": True,
        "error_only": True, "symmetric_error": True,
        "chi_source_ids": chi_ids, "chi_source_classes": chi_classes,
        "incoming_source_ids": list(incoming_info["source_ids"]),
        "leading_source_ids": list(leading_info["source_ids"]),
        "group": group, "branch": branch, "phase_cell": cell,
        "ordinary_Z_order": order, "ordinary_Z_piece": piece,
        "tail_starts_at_j": 2, "tail_formal_hb_power_added": tail_hb_power,
        "formal_hb_power_added": tail_hb_power + 1,
        "physical_hb_power_added": 1, "formal_external_powers": (1, 0, 0, -2),
        "coefficient_is_symmetric_radius": True,
    }


def _branch_paths(
    c: Any, *, group: str, branch: str, packet: Mapping[str, Mapping[int, Sequence[Mapping[str, Any]]]],
    leading: Mapping[int, Sequence[Mapping[str, Any]]], chi: Mapping[int, Sequence[Mapping[str, Any]]],
    cell: int, physical_factor: Any, constants: Mapping[str, Any],
) -> dict[int, list[dict[str, Any]]]:
    component = "Fxc" if branch.endswith("Fxc") else "Fc"
    source_scale = c.mpf(1) if component == "Fc" else c.mpf(1) * constants["sc"] / c.mpf(4)
    tail_hb_power = FORMAL_HB_POWER_ADDED[branch]
    base_factor = source_scale * constants["k0"]**2 * constants["T"] / c.mpf(2)
    chi_factor = source_scale * constants["kz"] * constants["k0"] * constants["T"]
    terms: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    if cell not in ACTIVE_CELLS:
        return terms
    for leading_index, left in enumerate(leading[0]):
        for incoming_index, right in enumerate(packet[component][0]):
            for chi_index, chi0 in enumerate(chi[0]):
                terms[0].append(_tail_term(
                    c, leading=left, incoming=right, chi_terms=(chi0,), factor=physical_factor * base_factor,
                    group=group, branch=branch, cell=cell, order=0, piece="leading*tail_Z0",
                    leading_index=leading_index, incoming_index=incoming_index, chi_index=chi_index,
                    chi_z_index=None, tail_hb_power=tail_hb_power,
                ))
    for leading_index, left in enumerate(leading[1]):
        for incoming_index, right in enumerate(packet[component][0]):
            for chi_index, chi0 in enumerate(chi[0]):
                terms[1].append(_tail_term(
                    c, leading=left, incoming=right, chi_terms=(chi0,), factor=physical_factor * base_factor,
                    group=group, branch=branch, cell=cell, order=1, piece="leading_Z*tail_Z0",
                    leading_index=leading_index, incoming_index=incoming_index, chi_index=chi_index,
                    chi_z_index=None, tail_hb_power=tail_hb_power,
                ))
    for leading_index, left in enumerate(leading[0]):
        for incoming_index, right in enumerate(packet[component][1]):
            for chi_index, chi0 in enumerate(chi[0]):
                terms[1].append(_tail_term(
                    c, leading=left, incoming=right, chi_terms=(chi0,), factor=physical_factor * base_factor,
                    group=group, branch=branch, cell=cell, order=1, piece="leading*tail_inlet_Z1",
                    leading_index=leading_index, incoming_index=incoming_index, chi_index=chi_index,
                    chi_z_index=None, tail_hb_power=tail_hb_power,
                ))
    for leading_index, left in enumerate(leading[0]):
        for incoming_index, right in enumerate(packet[component][0]):
            for chi_index, chi0 in enumerate(chi[0]):
                for chi_z_index, chi1 in enumerate(chi[1]):
                    terms[1].append(_tail_term(
                        c, leading=left, incoming=right, chi_terms=(chi0, chi1), factor=physical_factor * chi_factor,
                        group=group, branch=branch, cell=cell, order=1, piece="leading*tail_chi_Z1",
                        leading_index=leading_index, incoming_index=incoming_index, chi_index=chi_index,
                        chi_z_index=chi_z_index, tail_hb_power=tail_hb_power,
                    ))
    return terms


def build(
    c: Any, *, sc: Any, inlet_groups: Any, leading_F: Any, chi: Any,
    source_family: Any, source_atoms: Any, certificate: Any,
) -> dict[str, Any]:
    """Build cellwise symmetric Q/P inlet Kchi-tail swirl m15 errors."""

    if not all(hasattr(c, name) for name in ("mpf", "make_mpf", "exp")):
        raise TypeError("c must be a directed interval context")
    if source_family is None or (isinstance(source_family, str) and not source_family.strip()):
        raise ValueError("source_family is required")
    _validate_atoms(source_atoms)
    _validate_certificate(certificate)
    sc_iv = _iv(c, sc, "sc")
    if _lo(c, sc_iv) <= 0:
        raise ValueError("sc must be positive")
    sc_lo, sc_hi = _bounds(c, sc_iv)
    if sc_lo != sc_hi:
        raise ValueError("sc must be exact (a point interval)")
    groups = _normalize_groups(c, inlet_groups)
    leading = _normalize_leading(c, leading_F)
    chi_rows = _normalize_chi(c, chi)
    constants = _kernel_bounds(c, sc_iv, chi_rows)
    constants = dict(constants)
    constants["sc"] = sc_iv
    cell_geometry = {cell: _cell_factors(c, sc_iv, cell) for cell in ALL_CELLS}

    group_results: dict[str, Any] = {}
    for group in GROUP_NAMES:
        packet = groups[group]
        paths: dict[str, dict[int, dict[int, dict[str, list[dict[str, Any]]]]]] = {
            "swirl": {order: {cell: {branch: [] for branch in BRANCHES} for cell in ALL_CELLS} for order in ORDINARY_Z_ORDERS}
        }
        swirl: dict[int, dict[int, list[dict[str, Any]]]] = {order: {cell: [] for cell in ALL_CELLS} for order in ORDINARY_Z_ORDERS}
        radius_factors: dict[str, Any] = {}
        for cell in ALL_CELLS:
            physical_factor = cell_geometry[cell]["physical_factor"]
            for branch in BRANCHES:
                raw = _branch_paths(
                    c, group=group, branch=branch, packet=packet, leading=leading,
                    chi=chi_rows, cell=cell, physical_factor=physical_factor, constants=constants,
                )
                for order in ORDINARY_Z_ORDERS:
                    paths["swirl"][order][cell][branch] = _aggregate(raw[order])
                    swirl[order][cell].extend(paths["swirl"][order][cell][branch])
                component = "Fxc" if branch.endswith("Fxc") else "Fc"
                source_scale = c.mpf(1) if component == "Fc" else sc_iv / c.mpf(4)
                radius_factors[branch] = {
                    "Z0": source_scale * constants["k0"]**2 * constants["T"] / c.mpf(2),
                    "Z1_inlet": source_scale * constants["k0"]**2 * constants["T"] / c.mpf(2),
                    "Z1_chi": source_scale * constants["kz"] * constants["k0"] * constants["T"],
                    "formal_hb_power_added": FORMAL_HB_POWER_ADDED[branch],
                }
        group_results[group] = {
            "group": group, "group_path": list(GROUP_PATHS[group]),
            "group_kind": "error_only" if group in ERROR_GROUPS else "finite_QP",
            "error_only_group": group in ERROR_GROUPS,
            "input_error_only": group in ERROR_GROUPS,
            "input_components": list(INLET_COMPONENTS),
            "input_sector_counts": {component: {f"Z{order}": len(packet[component][order]) for order in ORDINARY_Z_ORDERS} for component in INLET_COMPONENTS},
            "swirl": swirl, "paths": paths,
            "radius_factors": radius_factors,
            "ordinary_Z1_product_rule": "leading_Z*tail_Z0 + leading*tail_inlet_Z1 + leading*tail_chi_Z1",
            "source_ids_retained": True, "chi_ids_retained": True,
            "error_provenance_retained": True, "all_outputs_error_only": True,
            "all_outputs_symmetric": True, "first_potential_reexecuted": False,
            "feedback_tail_included": False,
        }

    q_cells = [
        {
            "cell": cell, "q_lower": cell_geometry[cell]["q_lower"], "q_upper": cell_geometry[cell]["q_upper"],
            "q_interval": cell_geometry[cell]["q_interval"], "eta_mass": cell_geometry[cell]["eta_mass"],
            "eta_mass_exact": True, "velocity_cutoff_active": cell in ACTIVE_CELLS,
            "zero_mass_velocity_only": cell == 3,
        }
        for cell in ALL_CELLS
    ]
    return {
        "schema": SCHEMA, "status": "conditional_qp_j2_m15_symmetric_tail", "source_family": source_family, "source_atoms": source_atoms,
        "source_power_basis": list(SOURCE_POWER_BASIS), "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "inlet_components": list(INLET_COMPONENTS), "branches": list(BRANCHES),
        "inlet_groups": list(GROUP_NAMES), "group_results": group_results,
        "group_paths": {group: list(GROUP_PATHS[group]) for group in GROUP_NAMES},
        "chi_jets": {f"Z{order}": list(chi_rows[order]) for order in ORDINARY_Z_ORDERS},
        "leading_F_sector_counts": {f"Z{order}": len(leading[order]) for order in ORDINARY_Z_ORDERS},
        "q_cells": q_cells, "active_cells": list(ACTIVE_CELLS), "zero_mass_cell": 3,
        "eta_profile": "1-sigma(4q-2)", "eta_cell_masses": [cell_geometry[cell]["eta_mass"] for cell in ALL_CELLS],
        "eta_cell_masses_exact": True, "eta_integral": c.mpf(5) / c.mpf(8),
        "physical_channel": "swirl", "physical_swirl_formula": "-2*Lambda^-2*rho*Flead*F_tail*d_rho",
        "physical_coefficient_formula": "8*sc*M_i*(rho^2/16)",
        "rho_formula": "4*exp(hb*sc*q/4)", "d_rho_formula": "hb*(sc/4)*rho*dq",
        "radial_phase_enclosures": cell_geometry, "kernel_bounds": constants,
        "geometry_certificate": certificate, "geometry_certificate_caller_owned": True,
        "component_error_radii": {
            "F": {"B_formula": "|Fc| + hb*(sc/4)*|Fxc|", "BZ_formula": "|Fc_Z| + hb*(sc/4)*|Fxc_Z|", "R_F": "B*hb^4*k0^2*T/2", "R_FZ": "BZ*hb^4*k0^2*T/2 + B*hb^4*kz*k0*T"},
        },
        "ordinary_Z1_product_rule": "leading_Z*tail_Z0 + leading*tail_inlet_Z1 + leading*tail_chi_Z1",
        "formal_hb_power_added": dict(FORMAL_HB_POWER_ADDED),
        "formal_m15_hb_power_added": dict(FORMAL_M15_HB_ADDED),
        "formal_external_powers": (1, 0, 0, -2), "formal_powers_not_materialized": True,
        "source_ids_preserved": True, "chi_ids_preserved": True, "error_classes_preserved": True,
        "symmetric_error_outputs": True, "error_only_outputs_explicit": True,
        "signed_source_terms_used_as_radii": True, "signed_values_returned": False,
        "same_source_cancellation_claimed": False, "appended_to_212_ledger": False,
        "first_potential_added": False, "baseline_K0_reentered": False,
        "feedback_tail_included": False, "allword_tail_included": False,
        "m15_tail_error_computed": True, "m15_computed": False,
        "complete_m15_computed": False, "complete_m15": False,
        "terminal_closure_proved": False, "source_receipts_read": False,
        "old_builder_called": False, "ancestor_imports": False,
        "scope": "Q/P inlet homogeneous Kchi j>=2 symmetric swirl m15 error only; no first word or full feedback tail",
    }


build_qp_j2_m15 = build
build_tail = build
evaluate = build












