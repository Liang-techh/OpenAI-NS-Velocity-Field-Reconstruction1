"""Independent fixed-E00 collar seam leaf, version 3.

This leaf keeps integer Z indices internally and exposes the caller-facing
``inlet`` and ``chi_jets`` with the explicit ``Z0``/``Z1`` keys.  It contains
only the fixed H0 seam, algebraic K0 transfer, and the first homogeneous Kchi
word.  The j>=2 remainder and q>0 local forcing are deliberately omitted.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any

SCHEMA = "n1_e00_collar_fixed_seam_v4"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDERS = (0, 1)
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
CANONICAL_FIELDS = ("F", "U", "K", "P", "Fx", "Ux")
BRANCHES = ("F_from_Fc", "F_from_Fxc", "Fx_from_Fc", "Fx_from_Fxc")
FIRST_HB_POWER_ADDED = {"F_from_Fc": 2, "F_from_Fxc": 3, "Fx_from_Fc": 1, "Fx_from_Fxc": 2}

_ALIASES = {
    "F1": ("F1", "F", "Fc", "F_c"),
    "Uz1": ("Uz1", "U", "Uc", "U_c", "Uz"),
    "K1": ("K1", "K", "Kc", "K_c"),
    "P1": ("P1", "P", "Pc", "P_c"),
    "d_x_F1": ("d_x_F1", "Fx", "Fxc", "F_xc", "d_x_F"),
    "d_x_Uz1": ("d_x_Uz1", "Ux", "Uxc", "U_xc", "d_x_Uz"),
}
_K0 = {
    "F_from_Fc": ("F1", "F1", "F", 0),
    "F_from_Fxc": ("F1", "d_x_F1", "Fx", 1),
    "U_from_Uc": ("Uz1", "Uz1", "U", 0),
    "U_from_Uxc": ("Uz1", "d_x_Uz1", "Ux", 1),
    "K_from_Kc": ("K1", "K1", "K", 0),
    "K_from_Uxc": ("K1", "d_x_Uz1", "Ux", 1),
    "P_from_Pc": ("P1", "P1", "P", 0),
    "Fx_from_Fxc": ("d_x_F1", "d_x_F1", "Fx", 0),
    "Ux_from_Uxc": ("d_x_Uz1", "d_x_Uz1", "Ux", 0),
}
_POTENTIAL = {
    "F_from_Fc": ("F1", "F1", "F"),
    "F_from_Fxc": ("d_x_F1", "F1", "Fx"),
    "Fx_from_Fc": ("F1", "d_x_F1", "F"),
    "Fx_from_Fxc": ("d_x_F1", "d_x_F1", "Fx"),
}


def _seq(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _items(value: Any) -> list[Any]:
    return [] if value is None else list(value) if _seq(value) else [value]


def _uniq(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    for value in values:
        if not any(value == old for old in result):
            result.append(value)
    return result


def _is_iv(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _raw(value: Any, upper: bool) -> Any:
    bits = getattr(value, "_mpi_", None)
    if bits is not None:
        return bits[1 if upper else 0]
    bits = getattr(value, "_mpf_", None)
    return bits if bits is not None else value


def _make(c: Any, lower: Any, upper: Any) -> Any:
    if isinstance(lower, Fraction):
        lower = c.mpf(lower.numerator) / lower.denominator
    if isinstance(upper, Fraction):
        upper = c.mpf(upper.numerator) / upper.denominator
    return c.make_mpf((_raw(lower, False), _raw(upper, True)))


def _iv(c: Any, value: Any, label: str = "interval") -> Any:
    if _is_iv(value):
        return _make(c, value, value)
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / value.denominator
    if isinstance(value, Mapping):
        for key in ("value", "interval", "box", "enclosure", "data"):
            if key in value:
                return _iv(c, value[key], label)
        if "lower" in value and "upper" in value:
            return _make(c, value["lower"], value["upper"])
        if "lo" in value and "hi" in value:
            return _make(c, value["lo"], value["hi"])
        raise TypeError(f"{label} requires an interval-compatible value")
    if _seq(value) and len(value) == 2:
        return _make(c, value[0], value[1])
    try:
        return c.mpf(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{label} is not interval-compatible") from exc


def _bounds(c: Any, value: Any) -> tuple[Any, Any]:
    bits = _iv(c, value)._mpi_
    return c.make_mpf((bits[0], bits[0])), c.make_mpf((bits[1], bits[1]))


def _lo(c: Any, value: Any) -> Any:
    return _bounds(c, value)[0]


def _hi(c: Any, value: Any) -> Any:
    return _bounds(c, value)[1]


def _zero(c: Any) -> Any:
    return c.mpf(0)


def _one(c: Any) -> Any:
    return c.mpf(1)


def _exact_zero(c: Any, value: Any) -> bool:
    return _lo(c, value) == 0 and _hi(c, value) == 0


def _hull(c: Any, *values: Any) -> Any:
    return _zero(c) if not values else _make(c, min(_lo(c, value) for value in values), max(_hi(c, value) for value in values))


def _abs_upper(c: Any, value: Any) -> Any:
    low, high = _bounds(c, value)
    return max(-low, high)


def _symmetric(c: Any, value: Any) -> Any:
    radius = _abs_upper(c, value)
    return _make(c, -radius, radius)


def _refs(value: Any) -> list[Any]:
    if not isinstance(value, Mapping):
        return []
    refs = [value[key] for key in ("term_path", "path", "source_path", "operator_path", "old_full_prefix_ref", "old_full_prefix", "old_prefix_ref", "full_prefix_ref", "prefix_ref") if key in value]
    for key in ("origin", "provenance", "metadata"):
        child = value.get(key)
        if isinstance(child, Mapping):
            refs.extend(_refs(child))
    return _uniq(refs)


def _metadata(value: Mapping[str, Any]) -> tuple[list[Any], list[Any], bool, bool]:
    ids = list(_items(value.get("source_ids", value.get("source_id", value.get("input_id")))))
    classes = list(_items(value.get("error_classes", value.get("error_class", value.get("uncertainty_class")))))
    classes.extend(_items(value.get("source_error_classes")))
    error_only = bool(value.get("error_only", value.get("symmetric_error", value.get("is_error_only", False))))
    uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False))) or error_only
    for key in ("origin", "provenance", "metadata"):
        child = value.get(key)
        if isinstance(child, Mapping):
            classes.extend(_items(child.get("error_classes", child.get("error_class"))))
            classes.extend(_items(child.get("source_error_classes")))
            uncertain = uncertain or bool(child.get("uncertainty", False))
            error_only = error_only or bool(child.get("error_only", child.get("symmetric_error", False)))
    return _uniq(ids), _uniq(classes), uncertain, error_only


def _powers(value: Mapping[str, Any], label: str) -> tuple[int, int, int, int]:
    powers = value.get("powers", value.get("formal_powers", value.get("source_powers")))
    if powers is None:
        powers = (value.get("hb_power"), value.get("F0_power", value.get("F0_anchor_power")), value.get("Pstar_squared_power", value.get("Pstar_power")), value.get("Lambda_power", value.get("lambda_power")))
    if not _seq(powers) or len(powers) != 4 or any(isinstance(item, bool) or not isinstance(item, int) for item in powers):
        raise ValueError(f"{label} must carry four integer formal powers")
    return tuple(int(item) for item in powers)


def _term(c: Any, raw: Any, label: str) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must be a mapping sector")
    item = dict(raw)
    coefficient = next((item[key] for key in ("coefficient", "coefficient_interval", "value", "interval") if key in item), None)
    if coefficient is None:
        raise ValueError(f"{label} is missing coefficient")
    ids, classes, uncertain, error_only = _metadata(item)
    if not ids:
        raise ValueError(f"{label} must carry explicit source_ids")
    if error_only and "symmetric_error" not in classes:
        classes.append("symmetric_error")
    result = dict(item)
    result.update({"coefficient": _iv(c, coefficient, label + ".coefficient"), "powers": _powers(item, label), "source_ids": ids, "source_error_classes": list(classes), "error_classes": classes, "uncertainty": uncertain, "error_only": error_only, "symmetric_error": error_only, "term_path": item.get("term_path", item.get("path")), "path_refs": _refs(item)})
    return result


def _order(value: Any) -> int | None:
    if isinstance(value, int) and not isinstance(value, bool):
        return value if value in ORDERS else None
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"0", "z0", "value"}:
            return 0
        if text in {"1", "z1", "derivative", "ordinary_z1"}:
            return 1
    return None


def _term_like(value: Any) -> bool:
    return isinstance(value, Mapping) and any(key in value for key in ("coefficient", "coefficient_interval", "value", "interval", "powers"))


def _terms(c: Any, raw: Any, label: str) -> list[dict[str, Any]]:
    if raw is None:
        raise ValueError(f"{label} is required")
    values = [raw] if _term_like(raw) else list(raw) if _seq(raw) else [raw]
    return [_term(c, value, f"{label}[{index}]") for index, value in enumerate(values)]


def _rows(c: Any, raw: Any, label: str) -> dict[int, list[dict[str, Any]]]:
    if isinstance(raw, Mapping):
        for wrapper in ("ordinary_Z", "source_Z", "source_jets", "jets", "rows", "incoming"):
            if wrapper in raw and isinstance(raw[wrapper], (Mapping, Sequence)):
                raw = raw[wrapper]
                break
    if isinstance(raw, Mapping):
        found: dict[int, Any] = {}
        for key, value in raw.items():
            index = _order(key)
            if index is not None:
                if index in found:
                    raise ValueError(f"{label} supplies Z{index} twice")
                found[index] = value
        if set(found) != {0, 1}:
            raise ValueError(f"{label} requires explicit Z0 and Z1")
        return {index: _terms(c, found[index], f"{label}.Z{index}") for index in ORDERS}
    if _seq(raw) and len(raw) == 2:
        return {index: _terms(c, raw[index], f"{label}.Z{index}") for index in ORDERS}
    raise ValueError(f"{label} requires explicit Z0 and Z1")


def _lookup(raw: Mapping[Any, Any], names: Sequence[str]) -> Any:
    for name in names:
        if name in raw:
            return raw[name]
    lowered = {str(key).lower(): value for key, value in raw.items()}
    return next((lowered[name.lower()] for name in names if name.lower() in lowered), None)


def _unwrap(raw: Any) -> Any:
    if not isinstance(raw, Mapping):
        return raw
    for wrapper in ("inlet", "E00", "six_fields", "fields", "state", "components", "incoming"):
        if wrapper in raw and isinstance(raw[wrapper], Mapping):
            return raw[wrapper]
    return raw


def _packet(c: Any, raw: Any) -> dict[str, dict[int, list[dict[str, Any]]]]:
    raw = _unwrap(raw)
    if not isinstance(raw, Mapping):
        raise TypeError("inlet must be a six-field mapping")
    return {field: _rows(c, _lookup(raw, _ALIASES[field]), "inlet." + field) for field in FIELDS}


def _chi(c: Any, raw: Any) -> dict[int, list[dict[str, Any]]]:
    if raw is None:
        raise ValueError("chi is required")
    return _rows(c, raw, "chi")


def _atoms(raw: Any) -> list[Any]:
    if raw is None:
        raise ValueError("source_atoms is required")
    result = list(raw.keys()) if isinstance(raw, Mapping) else list(raw)
    if set(result) != set(SOURCE_POWER_BASIS):
        raise ValueError("source_atoms must contain exactly four formal atoms")
    return result


def _certificate(raw: Any) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise ValueError("ownership certificate is required")
    owners = raw.get("owners", raw.get("owner", {}))
    owners = owners if isinstance(owners, Mapping) else {}
    found: dict[str, Any] = {}
    for name in ("q", "x", "rho", "chi"):
        owner = next((raw[key] for key in (f"actual_{name}_owner", f"{name}_owner", f"{name}_owner_ref", f"{name}_owner_proved") if key in raw), None)
        if owner is None:
            owner = owners.get(name, owners.get("actual_" + name))
        if owner is None or owner is False or owner == "":
            raise ValueError(f"certificate lacks actual {name} owner")
        found[name] = owner
    if not (raw.get("fixed_geometry") is True or raw.get("ordinary_Z_geometry_fixed") is True or raw.get("geometry_Z_independent") is True):
        raise ValueError("certificate must prove fixed geometry")
    if not (raw.get("rho_radius_Z_independent") is True or raw.get("rho_Z_independent") is True):
        raise ValueError("certificate must prove rho/radius Z independence")
    return {"owners": found, "fixed_geometry": True, "rho_radius_Z_independent": True, "raw": dict(raw)}


def _q_cells(c: Any, raw: Any) -> list[dict[str, Any]]:
    if raw is None:
        raw = ((0, Fraction(1, 4)), (Fraction(1, 4), Fraction(1, 2)), (Fraction(1, 2), 1))
    if not _seq(raw) or len(raw) != 3:
        raise ValueError("exactly three q cells are required")
    result: list[dict[str, Any]] = []
    previous = None
    for index, cell in enumerate(raw):
        if isinstance(cell, Mapping):
            left, right = cell.get("left", cell.get("lower")), cell.get("right", cell.get("upper"))
        elif _seq(cell) and len(cell) == 2:
            left, right = cell
        else:
            raise ValueError("q cells require endpoint pairs")
        low, high = _iv(c, left, "q.left"), _iv(c, right, "q.right")
        low_value, high_value = _lo(c, low), _hi(c, high)
        if low_value < 0 or high_value > 1 or high_value <= low_value:
            raise ValueError("q cells must increase inside [0,1]")
        if previous is not None and low_value != previous:
            raise ValueError("q cells must meet exactly")
        previous = high_value
        result.append({"index": index, "q_interval": _hull(c, low, high), "left": low, "right": high})
    if previous != 1:
        raise ValueError("q cells must end at q=1")
    return result


def _k0_factors(c: Any, hb: Any, sc: Any, q: Any) -> dict[str, dict[str, Any]]:
    q = _iv(c, q, "q")
    zero = _zero(c)
    d_bound = hb * sc * q / 4
    d_coeff = sc * q / 4
    added = {"F_from_Fc": 0, "F_from_Fxc": 1, "U_from_Uc": 0, "U_from_Uxc": 1, "K_from_Kc": 0, "K_from_Uxc": 1, "P_from_Pc": 0, "Fx_from_Fxc": 0, "Ux_from_Uxc": 0}
    if _exact_zero(c, q):
        factors = {"F_from_Fc": _one(c), "F_from_Fxc": zero, "U_from_Uc": _one(c), "U_from_Uxc": zero, "K_from_Kc": _one(c), "K_from_Uxc": zero, "P_from_Pc": _one(c), "Fx_from_Fxc": _one(c), "Ux_from_Uxc": _one(c)}
        return {name: {"factor": factor, "formal_hb_power_added": added[name], "formula": "H0 identity at d=0"} for name, factor in factors.items()}
    phi = _make(c, c.exp(-_hi(c, d_bound)), _one(c))
    factors = {"F_from_Fc": _one(c), "F_from_Fxc": d_coeff * phi, "U_from_Uc": _one(c), "U_from_Uxc": d_coeff, "K_from_Kc": c.exp(-d_bound), "K_from_Uxc": -d_coeff * phi, "P_from_Pc": _one(c), "Fx_from_Fxc": c.exp(-c.mpf(3) * d_bound / 2), "Ux_from_Uxc": c.exp(-d_bound / 2)}
    formulas = {"F_from_Fc": "1", "F_from_Fxc": "hb*(sc*q/4)*phi(d)", "U_from_Uc": "1", "U_from_Uxc": "hb*sc*q/4", "K_from_Kc": "exp(-d)", "K_from_Uxc": "-hb*(sc*q/4)*phi(d)", "P_from_Pc": "1", "Fx_from_Fxc": "exp(-3*d/2)", "Ux_from_Uxc": "exp(-d/2)"}
    return {name: {"factor": factors[name], "formal_hb_power_added": added[name], "formula": formulas[name]} for name in factors}


def _route(c: Any, source: Mapping[str, Any], factor: Any, *, field: str, component: str, branch: str, source_order: int, hb_added: int, seam: bool) -> dict[str, Any]:
    identity = not _exact_zero(c, factor) and _lo(c, factor) == 1 and _hi(c, factor) == 1
    coefficient = source["coefficient"] if seam and identity else _zero(c) if seam else source["coefficient"] * factor
    error_only = bool(source.get("error_only") or source.get("symmetric_error"))
    if error_only and not seam:
        coefficient = _symmetric(c, coefficient)
    powers = tuple(source["powers"])
    source_classes = list(source.get("source_error_classes", source.get("error_classes", [])))
    return {"coefficient": coefficient, "powers": (powers[0] + hb_added, powers[1], powers[2], powers[3]), "incoming_powers": powers, "source_ids": list(source["source_ids"]), "incoming_source_ids": list(source["source_ids"]), "source_error_classes": source_classes, "error_classes": list(source.get("error_classes", source_classes)), "uncertainty": bool(source.get("uncertainty") or error_only), "error_only": error_only, "symmetric_error": error_only, "field": field, "source_component": component, "incoming_component": component, "branch": branch, "source_order": source_order, "term_path": source.get("term_path"), "path_refs": list(source.get("path_refs", [])), "old_full_prefix_refs": list(source.get("path_refs", [])), "factor": factor, "factor_is_identity": identity, "formal_hb_power_added": hb_added, "formal_powers_preserved": True, "signed_source_term": not error_only, "identity_difference": _zero(c) if seam else None, "identity_difference_exact_zero": bool(seam), "physical_radial_derivative_preserved": field in {"d_x_F1", "d_x_Uz1"}, "ordinary_Z_geometry_differentiated": False}


def _k0_stage(c: Any, packet: Mapping[str, Any], *, hb: Any, sc: Any, q: Any, seam: bool, index: int | None) -> dict[str, Any]:
    factors = _k0_factors(c, hb, sc, q)
    branches = {name: {"branch": name, "Z0": [], "Z1": []} for name in _K0}
    for branch, (field, source_field, component, _) in _K0.items():
        data = factors[branch]
        branches[branch].update({"field": field, "source_component": component, "factor": data["factor"], "factor_formula": data["formula"], "formal_hb_power_added": data["formal_hb_power_added"], "operator": "fixed-E00 K0 algebraic transfer", "algebraic_transfer_only": True, "append_to_old_coarse_prefix": False})
        for order in ORDERS:
            for source in packet[source_field][order]:
                branches[branch][f"Z{order}"].append(_route(c, source, data["factor"], field=field, component=component, branch=branch, source_order=order, hb_added=data["formal_hb_power_added"], seam=seam))
    return {"stage": "K0_transfer", "branches": branches, "q_interval": _iv(c, q), "q_index": index, "seam": seam, "six_component": True, "algebraic_transfer_only": True, "ordinary_Z1_product_rule": "factor*inlet_Z1 with fixed Z-independent factor", "ordinary_Z_geometry_differentiated": False, "physical_radial_derivative_preserved": True, "baseline_K0_reentered": False, "append_to_old_coarse_prefix": False}


def _potential_term(c: Any, chi: Mapping[str, Any], source: Mapping[str, Any], *, branch: str, field: str, component: str, chi_order: int, source_order: int, output_order: int, mass: Any, index: int | None, seam: bool) -> dict[str, Any]:
    base = tuple(a + b for a, b in zip(chi["powers"], source["powers"]))
    powers = (base[0] + FIRST_HB_POWER_ADDED[branch], base[1], base[2], base[3])
    ids = _uniq(list(chi["source_ids"]) + list(source["source_ids"]))
    source_classes = list(source.get("source_error_classes", source.get("error_classes", [])))
    classes = _uniq(list(chi.get("error_classes", [])) + list(source.get("error_classes", source_classes)))
    error_only = bool(chi.get("error_only") or source.get("error_only") or chi.get("symmetric_error") or source.get("symmetric_error"))
    if error_only and "symmetric_error" not in classes:
        classes.append("symmetric_error")
    coefficient = _zero(c) if seam else -2 * chi["coefficient"] * source["coefficient"] * mass
    if error_only and not seam:
        coefficient = _symmetric(c, coefficient)
    piece = "chi_Z*inlet_Z0" if (chi_order, source_order) == (1, 0) else "chi*inlet_Z1" if (chi_order, source_order) == (0, 1) else "chi*inlet_Z0"
    refs = _uniq(list(chi.get("path_refs", [])) + list(source.get("path_refs", [])))
    return {"coefficient": coefficient, "powers": powers, "source_ids": ids, "chi_source_ids": list(chi["source_ids"]), "incoming_source_ids": list(source["source_ids"]), "source_error_classes": source_classes, "error_classes": classes, "uncertainty": bool(error_only or chi.get("uncertainty") or source.get("uncertainty")), "error_only": error_only, "symmetric_error": error_only, "sentinel_error_propagated": any("sentinel_error" in str(item) for item in source_classes), "origin": {"operation": "first_Kchi_E00_fixed_seam_v4", "source_ids": ids, "source_error_classes": source_classes, "error_classes": classes, "leibniz_piece": piece, "old_full_prefix_refs": refs}, "branch": branch, "field": field, "incoming_component": component, "chi_order": chi_order, "incoming_order": source_order, "output_order": output_order, "source_order": source_order, "leibniz_piece": piece, "kernel_mass": mass, "kernel_mass_is_coefficient_only": True, "formal_hb_power_added": FIRST_HB_POWER_ADDED[branch], "coefficient_contains_formal_hb": False, "physical_dx_count": 1, "physical_radial_derivative_preserved": field == "d_x_F1", "signed_source_term": not error_only, "baseline_K0_included": False, "baseline_K0_reentered": False, "q_index": index, "source_term_path": source.get("term_path"), "chi_term_path": chi.get("term_path"), "path_refs": refs, "old_full_prefix_refs": refs, "exact_zero_at_q0": seam}


def _masses(c: Any, sc: Any, q: Any) -> dict[str, dict[str, Any]]:
    q = _iv(c, q, "q")
    zero = _zero(c)
    formulas = {"F_from_Fc": "2*(sc/8)^2*q^2*exp([-2Aq,+2Aq])", "F_from_Fxc": "(4/3)*(sc/8)^3*q^3*exp([-4Aq,+2Aq])", "Fx_from_Fc": "2*(sc/8)*q*exp([-3Aq,+Aq])", "Fx_from_Fxc": "2*(sc/8)^2*q^2*exp([-5Aq,+Aq])"}
    if _exact_zero(c, q):
        return {branch: {"mass": zero, "formula": formulas[branch]} for branch in BRANCHES}
    c0, A = sc / 8, sc / 16
    q2, q3 = q * q, q * q * q
    values = {"F_from_Fc": _hull(c, 2 * c0**2 * q2 * c.exp(-2 * A * q), 2 * c0**2 * q2 * c.exp(2 * A * q)), "F_from_Fxc": _hull(c, c.mpf(4) / 3 * c0**3 * q3 * c.exp(-4 * A * q), c.mpf(4) / 3 * c0**3 * q3 * c.exp(2 * A * q)), "Fx_from_Fc": _hull(c, 2 * c0 * q * c.exp(-3 * A * q), 2 * c0 * q * c.exp(A * q)), "Fx_from_Fxc": _hull(c, 2 * c0**2 * q2 * c.exp(-5 * A * q), 2 * c0**2 * q2 * c.exp(A * q))}
    return {branch: {"mass": values[branch], "formula": formulas[branch]} for branch in BRANCHES}


def _first_stage(c: Any, packet: Mapping[str, Any], chi: Mapping[int, Sequence[Mapping[str, Any]]], masses: Mapping[str, Any], *, index: int | None, seam: bool) -> dict[str, Any]:
    branches = {name: {"branch": name, "Z0": [], "Z1": []} for name in BRANCHES}
    for branch, (source_field, _, component) in _POTENTIAL.items():
        field = "F1" if branch.startswith("F_") else "d_x_F1"
        data = masses[branch]
        branches[branch].update({"field": field, "source_component": component, "operator": "-2*chi*K0*H0", "kernel_mass": data["mass"], "kernel_mass_formula": data["formula"], "formal_hb_power_added": FIRST_HB_POWER_ADDED[branch], "signed_sectors": True, "baseline_K0_included": False, "local_forcing_included": False})
        for chi_order, source_order, output_order in ((0, 0, 0), (1, 0, 1), (0, 1, 1)):
            for chi_term in chi[chi_order]:
                for source in packet[source_field][source_order]:
                    branches[branch][f"Z{output_order}"].append(_potential_term(c, chi_term, source, branch=branch, field=field, component=component, chi_order=chi_order, source_order=source_order, output_order=output_order, mass=data["mass"], index=index, seam=seam))
    return {"stage": "first_potential", "branches": branches, "first_potential_only": True, "baseline_K0_included": False, "j2_included": False, "local_forcing_included": False, "ordinary_Z1_product_rule": "chi_Z*inlet_Z0 + chi*inlet_Z1", "ordinary_Z_geometry_differentiated": False, "q_index": index, "exact_zero_at_q0": seam}


def _j2_stage(index: int | None, seam: bool) -> dict[str, Any]:
    branches = {name: {"branch": name, "Z0": [], "Z1": [], "field": "F1" if name.startswith("F_") else "d_x_F1", "formal_hb_power_added": FIRST_HB_POWER_ADDED[name] + 2, "j2_bound_status": "not_computed_in_this_leaf"} for name in BRANCHES}
    return {"stage": "j2_tail", "branches": branches, "structural_zero_fields": ["Uz1", "K1", "P1", "d_x_Uz1"], "j2_kernel_remainder_included": False, "j2_remainder_included": False, "j2_bound_status": "not_computed_in_this_leaf", "symmetric_error_outputs": True, "baseline_K0_included": False, "j1_included": False, "ordinary_Z1_product_rule": "chi_Z*inlet_Z0 + chi*inlet_Z1", "q_index": index, "exact_zero_at_q0": seam}


def _local(index: int | None, seam: bool) -> dict[str, Any]:
    return {"stage": "local_forcing", "branches": {}, "local_source_computed": False, "local_forcing_status": "exact_zero_at_q0" if seam else "not_computed", "contributes_to_inlet": False, "q_index": index, "exact_zero_at_q0": seam}


def _ledger(packet: Mapping[str, Any]) -> dict[str, dict[str, list[dict[str, Any]]]]:
    return {field: {f"Z{order}": [{"source_ids": list(term["source_ids"]), "powers": tuple(term["powers"]), "term_path": term.get("term_path"), "path_refs": list(term.get("path_refs", [])), "old_full_prefix_refs": list(term.get("path_refs", [])), "source_error_classes": list(term.get("source_error_classes", term.get("error_classes", []))), "error_classes": list(term.get("error_classes", []))} for term in packet[field][order]] for order in ORDERS} for field in FIELDS}


def _public_rows(rows: Mapping[int, Sequence[Mapping[str, Any]]]) -> dict[str, list[dict[str, Any]]]:
    return {f"Z{order}": [dict(term) for term in rows[order]] for order in ORDERS}


def _public_packet(packet: Mapping[str, Mapping[int, Sequence[Mapping[str, Any]]]]) -> dict[str, dict[str, list[dict[str, Any]]]]:
    return {field: _public_rows(rows) for field, rows in packet.items()}


def _seam(c: Any, packet: Mapping[str, Any], chi: Mapping[int, Sequence[Mapping[str, Any]]], *, hb: Any, sc: Any) -> dict[str, Any]:
    zero = _zero(c)
    masses = _masses(c, sc, zero)
    k0 = _k0_stage(c, packet, hb=hb, sc=sc, q=zero, seam=True, index=None)
    first = _first_stage(c, packet, chi, masses, index=None, seam=True)
    return {"q": zero, "x": c.mpf(2), "rho": c.mpf(4), "d": zero, "h0_identity": {"d": zero, "factors": {name: {"factor": data["factor"], "formal_hb_power_added": data["formal_hb_power_added"], "formula": data["formula"]} for name, data in _k0_factors(c, hb, sc, zero).items()}, "coefficient_differences": {field: {f"Z{order}": [zero for _ in packet[field][order]] for order in ORDERS} for field in FIELDS}, "componentwise_zero": True, "literal_directed_arithmetic": True, "no_midpoint": True, "ordinary_Z_geometry_differentiated": False, "physical_radial_derivative_preserved": True}, "k0_transfer": k0, "first_potential": first, "j2_tail": _j2_stage(None, True), "local_forcing": _local(None, True), "potential_kernel_masses": {name: data["mass"] for name, data in masses.items()}, "all_potential_masses_zero": all(_exact_zero(c, data["mass"]) for data in masses.values()), "local_forcing_exact_zero": True, "higher_potential_words_exact_zero": True, "source_term_ledger": _ledger(packet), "old_full_prefix_references_retained": True, "conditional_algebraic_seam_only": True, "same_source_function_identity_proved": False, "operator_equality_to_old_W1_proved": False, "no_dup_across_unrelated_owners_proved": False}


def build(c: Any, inlet: Any = None, chi: Any = None, *, source_family: Any = None, source_atoms: Any = None, certificate: Any = None, geometry_certificate: Any = None, q_cells: Any = None, hb: Any = None, sc: Any = None, **unexpected: Any) -> dict[str, Any]:
    """Build the V4 fixed-E00 seam and homogeneous first-word q cells."""
    if not all(hasattr(c, key) for key in ("mpf", "make_mpf", "exp")):
        raise TypeError("c must be a directed interval context")
    if unexpected:
        raise TypeError("unexpected fixed-seam options: " + ", ".join(sorted(unexpected)))
    if inlet is None or chi is None:
        raise TypeError("caller-owned inlet and chi are required")
    if source_family is None or source_family == "":
        raise ValueError("source_family is required")
    if geometry_certificate is not None:
        if certificate is not None and geometry_certificate != certificate:
            raise ValueError("certificates disagree")
        certificate = geometry_certificate
    cert = _certificate(certificate)
    atoms = _atoms(source_atoms)
    packet, chi_rows = _packet(c, inlet), _chi(c, chi)
    raw = cert["raw"]
    hb_input = _iv(c, raw.get("hb", (0, Fraction(1, 2))) if hb is None else hb, "hb")
    if _lo(c, hb_input) < 0 or _hi(c, hb_input) > c.mpf(1) / 2:
        raise ValueError("hb must lie in [0,1/2]")
    hb_interval = _make(c, _zero(c), c.mpf(1) / 2)
    sc_value = raw.get("sc") if sc is None else sc
    if sc_value is None:
        raise TypeError("explicit caller-owned sc is required")
    sc_interval = _iv(c, sc_value, "sc")
    if _lo(c, sc_interval) <= 0:
        raise ValueError("sc must be strictly positive")
    cells = _q_cells(c, q_cells if q_cells is not None else raw.get("q_cells"))
    seam = _seam(c, packet, chi_rows, hb=hb_interval, sc=sc_interval)
    cell_results: list[dict[str, Any]] = []
    for cell in cells:
        q = cell["q_interval"]
        masses = _masses(c, sc_interval, q)
        cell_results.append({"index": cell["index"], "q_interval": q, "q_left": cell["left"], "q_right": cell["right"], "k0_transfer": _k0_stage(c, packet, hb=hb_interval, sc=sc_interval, q=q, seam=False, index=cell["index"]), "first_potential": _first_stage(c, packet, chi_rows, masses, index=cell["index"], seam=False), "j2_tail": _j2_stage(cell["index"], False), "local_forcing": _local(cell["index"], False), "potential_kernel_masses": masses, "ordinary_Z1_product_rule": "chi_Z*inlet_Z0 + chi*inlet_Z1", "ordinary_Z_geometry_differentiated": False, "homogeneous": True, "append_to_old_coarse_prefix": False})
    geometry = {"fixed_target": "x=2,rho=4", "x_endpoint": c.mpf(2), "rho": c.mpf(4), "hb_interval": hb_interval, "hb_closure_only": True, "sc": sc_interval, "q_cells": [cell["q_interval"] for cell in cells], "ordinary_Z_geometry_fixed": True, "ordinary_Z_geometry_differentiated": False, "rho_radius_Z_independent": True, "physical_radial_derivative_preserved": True, "formal_hb_not_materialized": True, "formal_source_powers_not_materialized": True, "midpoint_values_used": False}
    return {"schema": SCHEMA, "status": "conditional_fixed_e00_collar_seam_and_q_transport_v4", "source_family": source_family, "source_atoms": atoms, "source_power_basis": list(SOURCE_POWER_BASIS), "fields": list(FIELDS), "canonical_fields": list(CANONICAL_FIELDS), "ordinary_Z_orders": [0, 1], "inlet": _public_packet(packet), "chi_jets": _public_rows(chi_rows), "source_term_ledger": _ledger(packet), "geometry": geometry, "geometry_certificate": cert, "certificate_passthrough": cert["raw"], "q_cells": cell_results, "seam": seam, "ordinary_Z1_product_rule": "chi_Z*inlet_Z0 + chi*inlet_Z1", "ordinary_Z_geometry_differentiated": False, "fixed_geometry_Z_independent": True, "physical_radial_derivative_preserved": True, "source_ids_preserved": True, "formal_powers_preserved": True, "error_classes_preserved": True, "signed_coefficients_retained": True, "uncertainty_separate": True, "error_only_outputs_explicit": True, "first_potential_hb_powers_added": dict(FIRST_HB_POWER_ADDED), "new_fixed_E00_source_distinct": True, "newly_fixed_missing_four_fields_consumed": True, "K0_j1_disjoint": True, "K0_transfer_algebraic_only": True, "legacy_coarse_W1_transfer_reused": False, "append_to_old_coarse_prefix": False, "old_W1_function_equal": False, "operator_equality_to_old_W1_proved": False, "no_dup_across_unrelated_owners_proved": False, "conditional_algebraic_seam_only": True, "same_source_function_identity_proved": False, "j2_remainder_included": False, "j2_bound_status": "not_computed_in_this_leaf", "local_forcing_exact_zero_at_q0": True, "higher_potential_words_exact_zero_at_q0": True, "m15_computed": False, "complete_collar": False, "complete_inlet": False, "old_builder_called": False, "ancestor_imports": False, "receipt_reads": False, "quadrature_used": False, "midpoint_values_used": False, "automatic_code_execution": False}


build_seam = build
build_transport = build
evaluate = build

__all__ = ["SCHEMA", "SOURCE_POWER_BASIS", "FIELDS", "CANONICAL_FIELDS", "BRANCHES", "FIRST_HB_POWER_ADDED", "build", "build_seam", "build_transport", "evaluate"]
