"""Cellwise physical m15 projection of the first Q/P collar potential word.

The caller owns the fixed-``x=2`` six-component Q/P inlet packets, the
whole-collar leading ``F`` sectors, and the two source-bound ``chi`` sectors.
This leaf constructs only the first potential word in the ``F`` channel and
then integrates its physical swirl contribution over the three nonzero
cutoff cells.  The no-potential endpoint result is deliberately not used as
a profile, and no old leaf, receipt, or producer is imported.

For a target ``x(q) = 2 exp(hb*sc*q/8)``, ``y = hb*sc*q/8`` and
``K0(x,t) = t*(1-t**2/x**2)/2``, the two positive kernel masses are bounded
after their formal ``hb`` powers are removed by
``M_Fc / hb**2 >= 2*(sc/8)**2*q**2*exp(-2*A_q)`` and
``M_Fxc / hb**3 >= (4/3)*(sc/8)**3*q**3*exp(-4*A_q)``,
with corresponding upper exponents ``2*A_q`` and ``2*A_q`` and
``A_q=sc*q_upper/16`` on a cell.  The second estimate retains the physical
``H0F=(1-4/t**2) Fxc`` factor.  The physical measure contributes one further
formal ``hb`` and the swirl channel contributes ``Lambda**-2``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Iterable

SCHEMA = "n1_collar_qp_j1_m15_cells_v1"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
VELOCITY_FIELDS = ("F1", "d_x_F1")
IGNORED_FIELDS = ("Uz1", "K1", "P1", "d_x_Uz1")
BRANCHES = ("F_from_Fc", "F_from_Fxc")
ALL_CELLS = (0, 1, 2, 3)
ACTIVE_CELLS = (0, 1, 2)
ETA_CELL_MASS = ((1, 4), (1, 4), (1, 8), (0, 1))
FORMAL_HB_POWER_ADDED = {"F_from_Fc": 2, "F_from_Fxc": 3}
PHYSICAL_EXTERNAL_POWERS = (1, 0, 0, -2)
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
MASS_FORMULAS = {
    "F_from_Fc": "2*(sc/8)^2*q^2*[exp(-2*A_q),exp(2*A_q)]",
    "F_from_Fxc": "(4/3)*(sc/8)^3*q^3*[exp(-4*A_q),exp(2*A_q)] with H0F=(1-4/t^2)*Fxc",
}

__all__ = [
    "SCHEMA", "SOURCE_POWER_BASIS", "ORDINARY_Z_ORDERS", "FIELDS",
    "VELOCITY_FIELDS", "IGNORED_FIELDS", "BRANCHES", "ALL_CELLS",
    "ACTIVE_CELLS", "GROUP_NAMES", "ERROR_GROUPS", "FORMAL_HB_POWER_ADDED",
    "MASS_FORMULAS", "build", "build_qp_j1_m15", "evaluate",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _items(value: Any) -> list[Any]:
    if value is None:
        return []
    return list(value) if _is_sequence(value) else [value]


def _unique(values: Iterable[Any]) -> list[Any]:
    result: list[Any] = []
    for value in values:
        if not any(value == old for old in result):
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


def _interval(c: Any, value: Any, label: str = "interval") -> Any:
    if _is_interval(value):
        return _make(c, value, value)
    if isinstance(value, Mapping):
        for key in ("value", "interval", "box", "enclosure", "data"):
            if key in value:
                return _interval(c, value[key], label)
        if "lower" in value and "upper" in value:
            return _make(c, value["lower"], value["upper"])
        if "lo" in value and "hi" in value:
            return _make(c, value["lo"], value["hi"])
        raise TypeError(f"{label} requires an interval-compatible value")
    if _is_sequence(value) and len(value) == 2:
        return _make(c, value[0], value[1])
    try:
        return c.mpf(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{label} is not interval-compatible") from exc


def _bounds(c: Any, value: Any) -> tuple[Any, Any]:
    raw = _interval(c, value)._mpi_
    return c.make_mpf((raw[0], raw[0])), c.make_mpf((raw[1], raw[1]))


def _scalar_bounds(c: Any, value: Any) -> tuple[Any, Any]:
    return _bounds(c, value)


def _zero(c: Any) -> Any:
    return c.mpf(0)


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _interval(c, value, label)
    lower, upper = _bounds(c, value)
    if upper < 0:
        raise ValueError(f"{label} is negative")
    if lower < 0:
        lower = _zero(c)
    return _make(c, lower, upper)


def _positive(c: Any, value: Any, label: str) -> Any:
    value = _interval(c, value, label)
    lower, _ = _bounds(c, value)
    if lower <= 0:
        raise ValueError(f"{label} must be strictly positive")
    return value


def _symmetric(c: Any, value: Any) -> Any:
    lower, upper = _bounds(c, value)
    radius = max(-lower, upper)
    return _make(c, -radius, radius)


def _expbox(c: Any, left: Any, right: Any) -> Any:
    lower, _ = _bounds(c, left)
    _, upper = _bounds(c, right)
    return _make(c, c.exp(lower), c.exp(upper))


def _freeze(value: Any) -> Any:
    try:
        hash(value)
    except TypeError:
        return (type(value).__name__, repr(value))
    return (type(value).__name__, value)


def _metadata(value: Any, *, seen: set[int] | None = None) -> dict[str, Any]:
    empty = {
        "source_ids": [], "error_classes": [], "uncertainty": False,
        "error_only": False, "symmetric_error": False,
    }
    if seen is None:
        seen = set()
    if not isinstance(value, Mapping):
        return empty
    marker = id(value)
    if marker in seen:
        return empty
    seen.add(marker)
    ids: list[Any] = []
    for key in ("source_ids", "sources", "source_id", "source", "input_id"):
        if key in value:
            ids.extend(_items(value[key]))
    classes: list[Any] = []
    for key in ("error_classes", "error_class", "uncertainty_class"):
        if key in value:
            classes.extend(_items(value[key]))
    operation = str(value.get("operation", "")).lower()
    kind = str(value.get("kind", value.get("error_kind", ""))).lower()
    uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False)))
    error_only = bool(value.get("error_only", value.get("is_error_only", False)))
    symmetric = bool(value.get("symmetric_error", value.get("is_symmetric_error", False)))
    uncertain = uncertain or any(token in operation or token in kind for token in ("error", "uncertainty", "symmetric"))
    symmetric = symmetric or "symmetric" in operation or "symmetric" in kind
    error_only = error_only or "error_only" in operation or "error_only" in kind
    for key in ("origin", "metadata", "parents"):
        child = value.get(key)
        if child is None or child is value:
            continue
        children = _items(child) if _is_sequence(child) else [child]
        for item in children:
            nested = _metadata(item, seen=seen)
            ids.extend(nested["source_ids"])
            classes.extend(nested["error_classes"])
            uncertain = uncertain or nested["uncertainty"]
            error_only = error_only or nested["error_only"]
            symmetric = symmetric or nested["symmetric_error"]
    return {
        "source_ids": _unique(ids), "error_classes": _unique(classes),
        "uncertainty": bool(uncertain), "error_only": bool(error_only),
        "symmetric_error": bool(symmetric),
    }


def _origin_descriptor(value: Any) -> Any:
    if isinstance(value, str):
        return value
    if isinstance(value, Mapping):
        return dict(value)
    return None


def _powers(value: Any, label: str) -> tuple[int, int, int, int]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} requires four formal powers")
    powers = value.get("powers", value.get("formal_powers"))
    if powers is None:
        powers = (
            value.get("hb_power", 0),
            value.get("F0_power", value.get("F0_anchor_power", 0)),
            value.get("Pstar_squared_power", value.get("Pstar_power", 0)),
            value.get("Lambda_power", value.get("lambda_power", 0)),
        )
    if not _is_sequence(powers) or len(powers) != 4 or any(isinstance(item, bool) or not isinstance(item, int) for item in powers):
        raise ValueError(f"{label} must contain four integer formal powers")
    return tuple(int(item) for item in powers)


def _term(c: Any, raw: Any, label: str, *, group_error: bool = False) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must be a sector mapping")
    coefficient = next((raw[key] for key in ("coefficient", "coefficient_interval", "value", "interval") if key in raw), None)
    if coefficient is None:
        raise ValueError(f"{label} is missing coefficient")
    powers = _powers(raw, label + ".powers")
    info = _metadata(raw)
    if not info["source_ids"]:
        raise ValueError(f"{label} must retain source_ids")
    classes = list(info["error_classes"])
    error_only = bool(group_error or info["error_only"] or info["symmetric_error"])
    symmetric = bool(group_error or info["symmetric_error"] or error_only)
    uncertain = bool(info["uncertainty"] or error_only)
    if symmetric and "symmetric_error" not in classes:
        classes.append("symmetric_error")
    origin_raw = raw.get("origin")
    origin = dict(origin_raw) if isinstance(origin_raw, Mapping) else {"source_ids": list(info["source_ids"])}
    origin.update({
        "source_ids": _unique(info["source_ids"]),
        "error_classes": _unique(classes),
        "uncertainty": uncertain,
        "error_only": error_only,
        "symmetric_error": symmetric,
    })
    result = dict(raw)
    result.update({
        "powers": powers,
        "coefficient": _interval(c, coefficient, label + ".coefficient"),
        "origin": origin,
        "source_ids": _unique(info["source_ids"]),
        "error_classes": _unique(classes),
        "uncertainty": uncertain,
        "error_only": error_only,
        "symmetric_error": symmetric,
    })
    if isinstance(origin_raw, str):
        result["origin_text"] = origin_raw
    return result


def _term_list(c: Any, raw: Any, label: str, *, allow_empty: bool, group_error: bool = False) -> list[dict[str, Any]]:
    if isinstance(raw, Mapping) and any(key in raw for key in ("coefficient", "coefficient_interval", "value", "interval", "powers")):
        values = [raw]
    elif _is_sequence(raw):
        values = list(raw)
    else:
        raise TypeError(f"{label} must be a sector sequence")
    if not values and not allow_empty:
        raise ValueError(f"{label} must not be empty")
    return [_term(c, value, f"{label}[{index}]", group_error=group_error) for index, value in enumerate(values)]


def _order(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"boolean is not an ordinary-Z order in {label}")
    if isinstance(value, int):
        result = value
    else:
        text = str(value).strip().lower()
        result = 0 if text in {"z0", "0", "value"} else 1 if text in {"z1", "1", "derivative", "ordinary_z1"} else -1
    if result not in ORDINARY_Z_ORDERS:
        raise ValueError(f"invalid ordinary-Z order {value!r} in {label}")
    return result


def _rows(c: Any, raw: Any, label: str, *, allow_empty: bool, group_error: bool = False) -> dict[int, list[dict[str, Any]]]:
    if isinstance(raw, Mapping):
        for wrapper in ("ordinary_Z", "source_Z", "jets", "rows", "incoming", "state"):
            if wrapper in raw and isinstance(raw[wrapper], (Mapping, Sequence)):
                raw = raw[wrapper]
                break
    if not isinstance(raw, Mapping):
        raise ValueError(f"{label} requires explicit ordinary Z0 and Z1 rows")
    found: dict[int, Any] = {}
    for key, value in raw.items():
        index = _order(key, label)
        if index in found:
            raise ValueError(f"{label} supplies Z{index} twice")
        found[index] = value
    if set(found) != set(ORDINARY_Z_ORDERS):
        raise ValueError(f"{label} requires explicit ordinary Z0 and Z1 rows")
    return {
        order: _term_list(c, found[order], f"{label}.Z{order}", allow_empty=allow_empty, group_error=group_error)
        for order in ORDINARY_Z_ORDERS
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
    result: dict[str, dict[str, dict[int, list[dict[str, Any]]]]] = {}
    for group in GROUP_NAMES:
        packet = value[group]
        if not isinstance(packet, Mapping) or set(packet) != set(FIELDS):
            raise ValueError(f"inlet_groups[{group}] must contain exactly {FIELDS}")
        group_error = group in ERROR_GROUPS
        result[group] = {
            field: _rows(
                c, packet[field], f"inlet_groups[{group}][{field}]",
                allow_empty=True, group_error=group_error,
            )
            for field in FIELDS
        }
    return result


def _normalize_leading(c: Any, value: Any, label: str) -> dict[int, list[dict[str, Any]]]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must provide explicit Z0 and Z1 rows")
    return _rows(c, value, label, allow_empty=False)


def _normalize_chi(c: Any, value: Any) -> dict[int, dict[str, Any]]:
    if value is None:
        raise ValueError("chi is required")
    raw = value
    if isinstance(raw, Mapping):
        for wrapper in ("ordinary_Z", "chi_jets", "jets", "rows"):
            if wrapper in raw and isinstance(raw[wrapper], (Mapping, Sequence)):
                raw = raw[wrapper]
                break
    rows = _rows(c, raw, "chi", allow_empty=False)
    result: dict[int, dict[str, Any]] = {}
    for order in ORDINARY_Z_ORDERS:
        if len(rows[order]) != 1:
            raise ValueError(f"chi Z{order} must contain exactly one source sector")
        term = rows[order][0]
        if tuple(term["powers"]) != (0, 0, 0, 0):
            raise ValueError(f"chi Z{order} must have zero formal powers")
        result[order] = term
    return result


def _validate_source_atoms(value: Any) -> None:
    if isinstance(value, Mapping):
        names = set(value)
    elif _is_sequence(value):
        names = {
            item.get("name", item.get("atom", item.get("power"))) if isinstance(item, Mapping) else item
            for item in value
        }
    else:
        raise ValueError("source_atoms must name the four formal source atoms")
    if names != set(SOURCE_POWER_BASIS):
        raise ValueError("source_atoms must contain exactly hb, F0_anchor, Pstar_squared, Lambda")


def _validate_certificate(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError("caller geometry and cutoff certificates are required")
    evidence = value.get("geometry", value)
    if not isinstance(evidence, Mapping):
        raise ValueError("geometry certificate must be a mapping")
    strict = any(evidence.get(name) is True or value.get(name) is True for name in (
        "hb_strictly_between_zero_and_half", "strict_hb_admitted", "strict_hb_admissibility",
    ))
    hb_sc = any(evidence.get(name) is True or value.get(name) is True for name in (
        "sc_phase_Z_independent", "hb_sc_Z_independent",
    ))
    rho = any(evidence.get(name) is True or value.get(name) is True for name in (
        "hb_epsilon_Z_independent", "rho_radius_Z_independent",
    ))
    cutoff = any(evidence.get(name) is True or value.get(name) is True for name in (
        "original_cutoff_mass_5_over_8", "cutoff_mass_5_over_8_exact", "actual_cutoff_complete_support_measure",
    ))
    cell_mass = any(evidence.get(name) is True or value.get(name) is True for name in (
        "original_cutoff_mass_by_cell", "cutoff_mass_by_cell_exact", "eta_cell_masses_exact",
    ))
    missing = []
    if not strict:
        missing.append("strict_hb_admitted")
    if not hb_sc:
        missing.append("hb_sc_Z_independent")
    if not rho:
        missing.append("rho_radius_Z_independent")
    if not cutoff:
        missing.append("original_cutoff_mass_5_over_8")
    if not cell_mass:
        missing.append("original_cutoff_mass_by_cell")
    if missing:
        raise ValueError("missing caller certificate: " + ", ".join(missing))
    return {
        "strict_hb_admitted": True, "hb_sc_Z_independent": True,
        "rho_radius_Z_independent": True, "original_cutoff_mass_5_over_8": True,
        "original_cutoff_mass_by_cell": True, "raw": dict(value),
    }


def _cell_factors(c: Any, sc: Any, cell: int) -> dict[str, Any]:
    q_lower = c.mpf(cell) / c.mpf(4)
    q_upper = c.mpf(cell + 1) / c.mpf(4)
    q_interval = _make(c, q_lower, q_upper)
    eta_num, eta_den = ETA_CELL_MASS[cell]
    eta_mass = c.mpf(eta_num) / c.mpf(eta_den)
    phase = sc / c.mpf(8)
    # y=hb*(sc/8)*q and 0<hb<1/2, so y<=A_q=sc*q_upper/16.
    A_q = sc * q_upper / c.mpf(16)
    exp_fc_lower = _expbox(c, -2 * A_q, -2 * A_q)
    exp_fc_upper = _expbox(c, 2 * A_q, 2 * A_q)
    exp_fxc_lower = _expbox(c, -4 * A_q, -4 * A_q)
    exp_fxc_upper = _expbox(c, 2 * A_q, 2 * A_q)
    mass_fc = _nonnegative(
        c,
        _make(
            c,
            c.mpf(2) * phase**2 * q_lower**2 * exp_fc_lower,
            c.mpf(2) * phase**2 * q_upper**2 * exp_fc_upper,
        ),
        "F<-Fc positive K0 mass",
    )
    mass_fxc = _nonnegative(
        c,
        _make(
            c,
            c.mpf(4) / c.mpf(3) * phase**3 * q_lower**3 * exp_fxc_lower,
            c.mpf(4) / c.mpf(3) * phase**3 * q_upper**3 * exp_fxc_upper,
        ),
        "F<-Fxc positive K0 mass",
    )
    # rho^2=16 exp(hb*sc*q/2), and hb<1/2 gives the cellwise upper
    # exponent sc*q_upper/4.  The lower endpoint is the exact positive 16.
    rho_squared = _make(c, c.mpf(16), c.mpf(16) * _expbox(c, sc * q_upper / c.mpf(4), sc * q_upper / c.mpf(4)))
    factors = {
        "F_from_Fc": _nonnegative(c, sc * eta_mass * rho_squared * mass_fc, "Fc physical factor"),
        "F_from_Fxc": _nonnegative(c, sc * eta_mass * rho_squared * mass_fxc, "Fxc physical factor"),
    }
    return {
        "cell": cell,
        "q_lower": q_lower,
        "q_upper": q_upper,
        "q_interval": q_interval,
        "eta_mass": eta_mass,
        "eta_mass_exact": True,
        "velocity_cutoff_active": cell in ACTIVE_CELLS,
        "A_q": A_q,
        "rho_squared_interval": rho_squared,
        "rho_formula": "4*exp(hb*sc*q/4)",
        "d_rho_formula": "hb*sc*rho/4*dq",
        "positive_kernel_masses": {"F_from_Fc": mass_fc, "F_from_Fxc": mass_fxc},
        "mass_formulas": dict(MASS_FORMULAS),
        "physical_factors": factors,
        "physical_factor_formulas": {
            "F_from_Fc": "(-2*Lambda^-2*rho*Flead)*(-2*chi*K0*Fc)*(hb*sc*rho/4*dq), integrated once",
            "F_from_Fxc": "(-2*Lambda^-2*rho*Flead)*(-2*chi*K0*(1-4/t^2)*Fxc)*(hb*sc*rho/4*dq), integrated once",
        },
        "kernel": "K0(x,t)=t*(1-t^2/x^2)/2",
        "target_x": "2*exp(hb*sc*q/8)",
        "target_q_is_variable": True,
        "endpoint_profile_reused": False,
    }


def _product(
    c: Any,
    leading: Mapping[str, Any],
    chi: Mapping[str, Any],
    inlet: Mapping[str, Any],
    factor: Any,
    *,
    group: str,
    branch: str,
    cell: int,
    order: int,
    piece: str,
    potential_piece: str,
    mass: Any,
) -> dict[str, Any]:
    powers0 = tuple(a + b + d for a, b, d in zip(leading["powers"], chi["powers"], inlet["powers"]))
    hb_added = FORMAL_HB_POWER_ADDED[branch]
    powers = (powers0[0] + hb_added + PHYSICAL_EXTERNAL_POWERS[0], powers0[1], powers0[2], powers0[3] + PHYSICAL_EXTERNAL_POWERS[3])
    left_info = _metadata(leading)
    chi_info = _metadata(chi)
    inlet_info = _metadata(inlet)
    group_error = group in ERROR_GROUPS
    error_only = bool(group_error or left_info["error_only"] or chi_info["error_only"] or inlet_info["error_only"])
    symmetric = bool(error_only or left_info["symmetric_error"] or chi_info["symmetric_error"] or inlet_info["symmetric_error"])
    classes = _unique(left_info["error_classes"] + chi_info["error_classes"] + inlet_info["error_classes"])
    if error_only and "symmetric_error" not in classes:
        classes.append("symmetric_error")
    uncertainty = bool(error_only or left_info["uncertainty"] or chi_info["uncertainty"] or inlet_info["uncertainty"])
    coefficient = leading["coefficient"] * chi["coefficient"] * inlet["coefficient"] * factor
    if error_only:
        coefficient = _symmetric(c, coefficient)
    source_ids = _unique(left_info["source_ids"] + chi_info["source_ids"] + inlet_info["source_ids"])
    path = {
        "group": group, "channel": "swirl", "field": "F1", "branch": branch,
        "phase_cell": cell, "ordinary_Z_order": order,
        "ordinary_Z_piece": piece, "potential_piece": potential_piece,
        "leading_source_ids": list(left_info["source_ids"]),
        "chi_source_ids": list(chi_info["source_ids"]),
        "inlet_source_ids": list(inlet_info["source_ids"]),
    }
    origin = {
        "operation": "qp_j1_physical_m15_product",
        "source_ids": source_ids,
        "leading_source_ids": list(left_info["source_ids"]),
        "chi_source_ids": list(chi_info["source_ids"]),
        "inlet_source_ids": list(inlet_info["source_ids"]),
        "error_classes": classes,
        "uncertainty": uncertainty,
        "error_only": error_only,
        "symmetric_error": symmetric,
        "group": group,
        "group_path": GROUP_PATHS[group],
        "channel": "swirl",
        "field": "F1",
        "branch": branch,
        "phase_cell": cell,
        "ordinary_Z_order": order,
        "ordinary_Z_piece": piece,
        "potential_piece": potential_piece,
        "potential_operator": "-2*chi*K0*H0F",
        "physical_operator": "-2*Lambda^-2*rho*Flead*d_rho",
        "kernel_mass": mass,
        "kernel_mass_is_positive": True,
        "kernel_mass_is_coefficient_only": True,
        "potential_formal_hb_power_added": hb_added,
        "physical_measure_formal_hb_power_added": 1,
        "formal_external_powers": PHYSICAL_EXTERNAL_POWERS,
        "coefficient_contains_formal_hb": False,
        "H0F_formula": "Fc+(1-4/t^2)*Fxc",
        "leading_origin": _origin_descriptor(leading.get("origin")),
        "chi_origin": _origin_descriptor(chi.get("origin")),
        "inlet_origin": _origin_descriptor(inlet.get("origin")),
        "path": path,
        "symbolic_cancellation_claimed": False,
    }
    return {
        "powers": powers,
        "coefficient": coefficient,
        "origin": origin,
        "source_ids": source_ids,
        "leading_source_ids": list(left_info["source_ids"]),
        "chi_source_ids": list(chi_info["source_ids"]),
        "inlet_source_ids": list(inlet_info["source_ids"]),
        "error_classes": classes,
        "uncertainty": uncertainty,
        "error_only": error_only,
        "symmetric_error": symmetric,
        "group": group,
        "channel": "swirl",
        "field": "F1",
        "branch": branch,
        "phase_cell": cell,
        "ordinary_Z_order": order,
        "ordinary_Z_piece": piece,
        "potential_piece": potential_piece,
        "kernel_mass": mass,
        "physical_factor": factor,
        "potential_formal_hb_power_added": hb_added,
        "physical_measure_formal_hb_power_added": 1,
        "formal_external_powers": PHYSICAL_EXTERNAL_POWERS,
        "signed_source_term": not error_only,
        "baseline_K0_reentered": False,
        "endpoint_profile_reused": False,
    }


def _group_key(term: Mapping[str, Any]) -> tuple[Any, ...]:
    classes = tuple(sorted((_freeze(value) for value in term.get("error_classes", ())), key=repr))
    return tuple(term["powers"]), bool(term.get("error_only")), bool(term.get("symmetric_error")), classes


def _aggregate(c: Any, terms: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[Any, ...], dict[str, Any]] = {}
    for term in terms:
        key = _group_key(term)
        info = groups.get(key)
        if info is None:
            info = {
                "powers": tuple(term["powers"]), "coefficient": term["coefficient"],
                "source_ids": [], "leading_source_ids": [], "chi_source_ids": [], "inlet_source_ids": [],
                "error_classes": list(term.get("error_classes", ())),
                "uncertainty": bool(term.get("uncertainty")),
                "error_only": bool(term.get("error_only")),
                "symmetric_error": bool(term.get("symmetric_error")),
                "path_count": 0, "source_count": 0, "merged_term_count": 0,
                "paths": {},
            }
            groups[key] = info
        else:
            info["coefficient"] = info["coefficient"] + term["coefficient"]
        for target, source in (
            ("source_ids", term.get("source_ids", ())),
            ("leading_source_ids", term.get("leading_source_ids", ())),
            ("chi_source_ids", term.get("chi_source_ids", ())),
            ("inlet_source_ids", term.get("inlet_source_ids", ())),
        ):
            info[target] = _unique(info[target] + list(source))
        info["error_classes"] = _unique(info["error_classes"] + list(term.get("error_classes", ())))
        info["uncertainty"] = bool(info["uncertainty"] or term.get("uncertainty"))
        info["error_only"] = bool(info["error_only"] or term.get("error_only"))
        info["symmetric_error"] = bool(info["symmetric_error"] or term.get("symmetric_error"))
        info["path_count"] += 1
        info["source_count"] += max(1, len(term.get("source_ids", ())))
        info["merged_term_count"] += 1
        origin = term.get("origin")
        path = origin.get("path") if isinstance(origin, Mapping) else None
        if not isinstance(path, Mapping):
            path = {
                "group": term.get("group"), "channel": term.get("channel", "swirl"),
                "field": term.get("field", "F1"), "branch": term.get("branch"),
                "phase_cell": term.get("phase_cell"), "ordinary_Z_order": term.get("ordinary_Z_order"),
                "ordinary_Z_piece": term.get("ordinary_Z_piece"),
            }
        path_key = tuple(path.get(name) for name in ("group", "channel", "field", "branch", "phase_cell", "ordinary_Z_order", "ordinary_Z_piece"))
        info["paths"][path_key] = info["paths"].get(path_key, 0) + 1
    result: list[dict[str, Any]] = []
    for info in groups.values():
        path_entries = [{"path": list(path), "count": count} for path, count in info["paths"].items()]
        origin = {
            "operation": "aggregate",
            "source_ids": list(info["source_ids"]),
            "leading_source_ids": list(info["leading_source_ids"]),
            "chi_source_ids": list(info["chi_source_ids"]),
            "inlet_source_ids": list(info["inlet_source_ids"]),
            "error_classes": list(info["error_classes"]),
            "uncertainty": bool(info["uncertainty"]),
            "error_only": bool(info["error_only"]),
            "symmetric_error": bool(info["symmetric_error"]),
            "path_count": info["path_count"],
            "source_count": info["source_count"],
            "merged_term_count": info["merged_term_count"],
            "paths": path_entries,
            "symbolic_cancellation_claimed": False,
        }
        coefficient = info["coefficient"]
        if info["error_only"]:
            coefficient = _symmetric(c, coefficient)
        result.append({
            "powers": info["powers"], "coefficient": coefficient, "origin": origin,
            "source_ids": list(info["source_ids"]),
            "leading_source_ids": list(info["leading_source_ids"]),
            "chi_source_ids": list(info["chi_source_ids"]),
            "inlet_source_ids": list(info["inlet_source_ids"]),
            "error_classes": list(info["error_classes"]),
            "uncertainty": bool(info["uncertainty"]),
            "error_only": bool(info["error_only"]),
            "symmetric_error": bool(info["symmetric_error"]),
            "merged_term_count": info["merged_term_count"],
            "path_count": info["path_count"],
        })
    return result


def _branch_terms(
    c: Any,
    leading: Mapping[int, list[dict[str, Any]]],
    chi: Mapping[int, Mapping[str, Any]],
    packet: Mapping[str, Mapping[int, list[dict[str, Any]]]],
    *,
    group: str,
    branch: str,
    cell: int,
    factor: Any,
    mass: Any,
    order: int,
) -> list[dict[str, Any]]:
    source_field = "F1" if branch == "F_from_Fc" else "d_x_F1"
    if order == 0:
        pieces = [("leading*chi*inlet_Z0", 0, 0, 0)]
    else:
        pieces = [
            ("leading_Z*chi*inlet_Z0", 1, 0, 0),
            ("leading*chi_Z*inlet_Z0", 0, 1, 0),
            ("leading*chi*inlet_Z1", 0, 0, 1),
        ]
    result: list[dict[str, Any]] = []
    for piece, leading_order, chi_order, inlet_order in pieces:
        potential_piece = "chi_Z*inlet_Z0" if chi_order == 1 else "chi*inlet_Z1" if inlet_order == 1 else "chi*inlet_Z0"
        for left in leading[leading_order]:
            for inlet in packet[source_field][inlet_order]:
                result.append(_product(
                    c, left, chi[chi_order], inlet, factor,
                    group=group, branch=branch, cell=cell, order=order,
                    piece=piece, potential_piece=potential_piece, mass=mass,
                ))
    return result


def build(
    c: Any,
    *,
    sc: Any,
    leading_F: Any,
    chi: Any,
    inlet_groups: Any = None,
    source_family: Any = None,
    source_atoms: Any = None,
    certificate: Any = None,
    groups: Any = None,
    incoming: Any = None,
) -> dict[str, Any]:
    """Build the standalone first-potential physical swirl sectors.

    The only accepted physical channel is ``F``.  ``Z1`` is assembled from
    the three Leibniz pieces
    ``leading_Z*chi*inlet_Z0 + leading*chi_Z*inlet_Z0 + leading*chi*inlet_Z1``.
    Every active cell receives one positive, cell-specific integrated factor;
    no endpoint value is reused as a q profile.
    """

    if not all(hasattr(c, key) for key in ("mpf", "make_mpf", "exp")):
        raise TypeError("c must be a directed interval context")
    if source_family is None or (isinstance(source_family, str) and not source_family.strip()):
        raise ValueError("source_family is required")
    if source_atoms is None:
        raise ValueError("source_atoms is required")
    _validate_source_atoms(source_atoms)
    cert = _validate_certificate(certificate)
    if inlet_groups is not None and (groups is not None or incoming is not None):
        raise ValueError("inlet_groups cannot be mixed with groups or incoming")
    raw_groups = inlet_groups if inlet_groups is not None else groups if groups is not None else incoming
    if raw_groups is None:
        raise ValueError("inlet_groups is required")
    sc_iv = _interval(c, sc, "sc")
    if not _is_interval(sc_iv):
        raise TypeError("sc must be a directed interval")
    sc_lower, sc_upper = _scalar_bounds(c, sc_iv)
    limit_lower, _ = _scalar_bounds(c, c.mpf(1) / c.mpf(1024))
    if sc_lower <= 0 or sc_upper >= limit_lower:
        raise ValueError("sc must satisfy 0 < sc < 1/1024")
    if sc_lower != sc_upper:
        raise ValueError("sc must be exact (a point interval)")
    leading = _normalize_leading(c, leading_F, "leading_F")
    chi_rows = _normalize_chi(c, chi)
    packets = _normalize_groups(c, raw_groups)
    cell_geometry = {cell: _cell_factors(c, sc_iv, cell) for cell in ALL_CELLS}

    group_results: dict[str, Any] = {}
    for group in GROUP_NAMES:
        packet = packets[group]
        paths: dict[int, dict[int, dict[str, list[dict[str, Any]]]]] = {
            order: {cell: {branch: [] for branch in BRANCHES} for cell in ALL_CELLS}
            for order in ORDINARY_Z_ORDERS
        }
        for order in ORDINARY_Z_ORDERS:
            for cell in ACTIVE_CELLS:
                geometry = cell_geometry[cell]
                for branch in BRANCHES:
                    raw_terms = _branch_terms(
                        c, leading, chi_rows, packet,
                        group=group, branch=branch, cell=cell,
                        factor=geometry["physical_factors"][branch],
                        mass=geometry["positive_kernel_masses"][branch],
                        order=order,
                    )
                    paths[order][cell][branch] = _aggregate(c, raw_terms)
        group_results[group] = {
            "group": group,
            "group_path": GROUP_PATHS[group],
            "group_kind": "error_only" if group in ERROR_GROUPS else "finite_QP",
            "error_only_group": group in ERROR_GROUPS,
            "input_error_only": group in ERROR_GROUPS,
            "paths": {"swirl": paths},
            "potential_channel": "F_only",
            "potential_branches": list(BRANCHES),
            "potential_operator": "-2*chi*K0*H0F",
            "physical_operator": "-2*Lambda^-2*rho*Flead*d_rho",
            "velocity_fields_used": list(VELOCITY_FIELDS),
            "ignored_fields": list(IGNORED_FIELDS),
            "P1_is_not_m15": True,
            "K1_is_not_m15": True,
            "axial_not_added": True,
            "pressure_not_added": True,
            "ordinary_Z1_product_rule": "leading_Z*chi*inlet_Z0 + leading*chi_Z*inlet_Z0 + leading*chi*inlet_Z1",
            "source_ids_retained": True,
            "chi_source_ids_retained": True,
            "error_provenance_retained": True,
            "input_sector_counts": {
                field: {f"Z{order}": len(packet[field][order]) for order in ORDINARY_Z_ORDERS}
                for field in FIELDS
            },
        }

    eta_masses = [c.mpf(num) / c.mpf(den) for num, den in ETA_CELL_MASS]
    q_cells = [
        {
            "cell": cell,
            "q_lower": cell_geometry[cell]["q_lower"],
            "q_upper": cell_geometry[cell]["q_upper"],
            "q_interval": cell_geometry[cell]["q_interval"],
            "eta_mass": eta_masses[cell],
            "eta_mass_exact": True,
            "velocity_cutoff_active": cell in ACTIVE_CELLS,
        }
        for cell in ALL_CELLS
    ]
    return {
        "schema": SCHEMA,
        "status": "conditional_partial_qp_first_potential_physical_swirl_m15",
        "source_family": source_family,
        "source_atoms": source_atoms,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "fields": list(FIELDS),
        "velocity_fields": list(VELOCITY_FIELDS),
        "ignored_fields": list(IGNORED_FIELDS),
        "inlet_groups": list(GROUP_NAMES),
        "group_results": group_results,
        "group_paths": {group: list(GROUP_PATHS[group]) for group in GROUP_NAMES},
        "q_cells": q_cells,
        "eta_profile": "1-sigma(4q-2)",
        "eta_cell_masses": eta_masses,
        "eta_cell_masses_exact": True,
        "eta_integral": c.mpf(5) / c.mpf(8),
        "cutoff_mass": c.mpf(5) / c.mpf(8),
        "cutoff_mass_by_cell": eta_masses,
        "positive_cutoff_mass_exact": True,
        "active_cells": list(ACTIVE_CELLS),
        "zero_mass_cell": 3,
        "zero_mass_applies_to_velocity_only": True,
        "radial_phase_enclosures": cell_geometry,
        "rho_formula": "4*exp(hb*sc*q/4)",
        "d_rho_formula": "hb*sc*rho/4*dq",
        "target_x_formula": "2*exp(hb*sc*q/8)",
        "target_q_is_variable": True,
        "H0F_formula": "Fc+(1-4/t^2)*Fxc",
        "K0_formula": "K0(x,t)=t*(1-t^2/x^2)/2",
        "potential_formula": "F_j1(q)=-2*chi*integral_2^x(q) K0(x(q),t)*H0F(t)dt",
        "swirl_formula": "-2*Lambda^-2*rho*Flead*F_j1*d_rho",
        "ordinary_Z1_product_rule": "leading_Z*chi*inlet_Z0 + leading*chi_Z*inlet_Z0 + leading*chi*inlet_Z1",
        "formal_potential_hb_powers": dict(FORMAL_HB_POWER_ADDED),
        "formal_external_powers": PHYSICAL_EXTERNAL_POWERS,
        "formal_hb_measure": True,
        "formal_source_powers_not_materialized": True,
        "formal_hb_not_materialized": True,
        "physical_epsilon_scale_swirl": "Lambda^-2",
        "geometry_certificate_caller_owned": True,
        "geometry_certificate_independently_verified": False,
        "geometry_certificate": cert,
        "certificate_passthrough": cert["raw"],
        "strict_hb_admitted": True,
        "admissible_hb": "0<hb<1/2",
        "positive_kernel_masses": True,
        "physical_measure_integrated_once": True,
        "source_ids_union_retained": True,
        "error_classes_union_retained": True,
        "uncertainty_retained": True,
        "independent_group_results": True,
        "aggregation_key": "(group, branch, cell, ordinary-Z order, powers, uncertainty, error classes)",
        "aggregation_is_outward_interval_addition_only": True,
        "symmetric_error_outputs_separate": True,
        "new_QP_feedback_added": False,
        "first_potential_added": True,
        "only_F_swirl_added": True,
        "finite_physical_qp_j1_m15_computed": True,
        "physical_m15_finite_contribution_computed": True,
        "axial_added": False,
        "pressure_added": False,
        "old_endpoint_correction_reused": False,
        "old_leaf_called": False,
        "endpoint_profile_reused": False,
        "allword_tail_added": False,
        "appended_to_ledger": False,
        "terminal_closure_proved": False,
        "m15_computed": False,
        "complete_m15_computed": False,
        "complete_collar_field": False,
        "scope": "standalone finite physical first-collar Q/P first-potential F-to-swirl m15 cell projection; no axial, pressure, feedback, or terminal closure",
    }


build_qp_j1_m15 = build
evaluate = build
