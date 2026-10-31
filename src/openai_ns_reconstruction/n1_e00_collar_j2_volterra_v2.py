"""Fixed-E00 Kchi Volterra ``j >= 2`` tail.

The leaf consumes caller-owned interval rows only.  It returns symmetric
error enclosures and never evaluates a source function or writes a ledger.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any

SCHEMA = "n1_e00_collar_j2_volterra_v2"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDERS = (0, 1)
ORDINARY_Z_ORDERS = ORDERS
INLET_COMPONENTS = ("F1", "Fx")
FIELDS = INLET_COMPONENTS
CANONICAL_FIELDS = ("F", "Fx")
LEGACY_FIELDS = ("F1", "d_x_F1")
BRANCHES = ("F_from_Fc", "F_from_Fxc", "Fx_from_Fc", "Fx_from_Fxc")
FORMAL_HB_POWER_ADDED = {
    "F_from_Fc": 4,
    "F_from_Fxc": 5,
    "Fx_from_Fc": 3,
    "Fx_from_Fxc": 4,
}
J2_HB_POWER_ADDED = dict(FORMAL_HB_POWER_ADDED)
_SOURCE_ALIASES = {
    "F1": ("F1", "Fc", "F_c", "F"),
    "Fx": ("Fx", "Fxc", "F_xc", "d_x_F1", "d_x_F"),
}
_COEFFICIENT_KEYS = (
    "coefficient", "coefficient_interval", "value", "interval", "box", "enclosure",
)
_POWER_KEYS = ("powers", "formal_powers", "source_powers")

__all__ = [
    "SCHEMA", "SOURCE_POWER_BASIS", "ORDERS", "ORDINARY_Z_ORDERS",
    "INLET_COMPONENTS", "FIELDS", "CANONICAL_FIELDS", "LEGACY_FIELDS",
    "BRANCHES", "FORMAL_HB_POWER_ADDED", "J2_HB_POWER_ADDED",
    "kernel_bounds", "tail_kernel_bounds", "build", "build_tail",
    "build_j2_tail", "evaluate",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _items(value: Any) -> list[Any]:
    return [] if value is None else list(value) if _is_sequence(value) else [value]


def _unique(values: Sequence[Any]) -> list[Any]:
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
    if isinstance(lower, Fraction):
        lower = c.mpf(lower.numerator) / lower.denominator
    if isinstance(upper, Fraction):
        upper = c.mpf(upper.numerator) / upper.denominator
    if not hasattr(lower, "_mpi_") and not hasattr(lower, "_mpf_"):
        lower = c.mpf(lower)
    if not hasattr(upper, "_mpi_") and not hasattr(upper, "_mpf_"):
        upper = c.mpf(upper)
    return c.make_mpf((_raw(lower, False), _raw(upper, True)))


def _iv(c: Any, value: Any, label: str = "interval") -> Any:
    """Coerce an interval without choosing a midpoint."""
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
    raw = _iv(c, value)._mpi_
    return c.make_mpf((raw[0], raw[0])), c.make_mpf((raw[1], raw[1]))


def _lo(c: Any, value: Any) -> Any:
    return _bounds(c, value)[0]


def _hi(c: Any, value: Any) -> Any:
    return _bounds(c, value)[1]


def _zero(c: Any) -> Any:
    return c.mpf(0)


def _exact_zero(c: Any, value: Any) -> bool:
    return _lo(c, value) == 0 and _hi(c, value) == 0


def _abs_upper(c: Any, value: Any) -> Any:
    lo, hi = _bounds(c, value)
    if lo >= 0:
        return hi
    if hi <= 0:
        return -lo
    return max(-lo, hi)


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _iv(c, value, label)
    lo, hi = _bounds(c, value)
    if hi < 0:
        raise ValueError(f"{label} must be nonnegative")
    if lo < 0:
        lo = _zero(c)
    return _make(c, lo, hi)


def _positive(c: Any, value: Any, label: str) -> Any:
    value = _iv(c, value, label)
    if _lo(c, value) <= 0:
        raise ValueError(f"{label} must be strictly positive")
    return value


def _symmetric(c: Any, radius: Any) -> Any:
    radius = _hi(c, _nonnegative(c, radius, "symmetric radius"))
    return _make(c, -radius, radius)


def _lookup(raw: Mapping[Any, Any], names: Sequence[str], label: str) -> Any:
    for name in names:
        if name in raw:
            return raw[name]
    lowered = {str(key).lower(): value for key, value in raw.items()}
    for name in names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    raise KeyError(f"{label} is required")


def _order(key: Any) -> int | None:
    if isinstance(key, int) and not isinstance(key, bool):
        return key if key in ORDERS else None
    if isinstance(key, str):
        text = key.lower()
        if text in {"0", "z0", "value", "ordinary_z0"}:
            return 0
        if text in {"1", "z1", "derivative", "ordinary_z1", "ordinary_z"}:
            return 1
    return None


def _term_like(value: Any) -> bool:
    return isinstance(value, Mapping) and any(
        key in value for key in _COEFFICIENT_KEYS + _POWER_KEYS
    )


def _metadata(raw: Mapping[str, Any]) -> tuple[list[Any], list[Any], list[Any], bool, bool]:
    ids: list[Any] = []
    source_classes: list[Any] = []
    error_classes: list[Any] = []
    for key in ("source_ids", "sources", "source_id", "source", "input_id"):
        if key in raw:
            ids.extend(_items(raw[key]))
    for key in ("source_error_classes", "source_error_class"):
        if key in raw:
            source_classes.extend(_items(raw[key]))
    for key in ("error_classes", "error_class", "uncertainty_class"):
        if key in raw:
            error_classes.extend(_items(raw[key]))
    for key in ("origin", "provenance", "metadata"):
        child = raw.get(key)
        if isinstance(child, Mapping):
            nested = _metadata(child)
            # A provenance/origin ID does not satisfy the explicit source-ID
            # requirement of this leaf; callers must bind source_ids on the
            # term itself.
            source_classes.extend(nested[1])
            error_classes.extend(nested[2])
    error_only = bool(raw.get("error_only", raw.get("symmetric_error", False)))
    uncertain = bool(raw.get("uncertainty", raw.get("is_uncertainty", False))) or error_only
    return _unique(ids), _unique(source_classes), _unique(error_classes), uncertain, error_only


def _powers(raw: Mapping[str, Any], label: str) -> tuple[int, int, int, int]:
    powers: Any = next((raw[key] for key in _POWER_KEYS if key in raw), None)
    if powers is None:
        powers = (
            raw.get("hb_power"),
            raw.get("F0_power", raw.get("F0_anchor_power")),
            raw.get("Pstar_squared_power", raw.get("Pstar_power")),
            raw.get("Lambda_power", raw.get("lambda_power")),
        )
    if not _is_sequence(powers) or len(powers) != 4 or any(
        isinstance(item, bool) or not isinstance(item, int) for item in powers
    ):
        raise ValueError(f"{label}.powers must contain four integer formal powers")
    return tuple(int(item) for item in powers)


def _term(c: Any, raw: Any, label: str, default_id: str, *, require_powers: bool = True) -> dict[str, Any]:
    item = dict(raw) if isinstance(raw, Mapping) else {"coefficient": raw}
    coefficient = next((item[key] for key in _COEFFICIENT_KEYS if key in item), None)
    if coefficient is None:
        raise ValueError(f"{label} is missing coefficient interval")
    ids, source_classes, error_classes, uncertain, error_only = _metadata(item)
    if not ids:
        raise ValueError(f"{label} must retain at least one source ID")
    if require_powers:
        powers = _powers(item, label)
    else:
        try:
            powers = _powers(item, label)
        except ValueError:
            powers = (0, 0, 0, 0)
    if not error_classes:
        error_classes = list(source_classes)
    return {
        **item,
        "coefficient": _iv(c, coefficient, label + ".coefficient"),
        "powers": powers,
        "source_ids": _unique(ids),
        "source_error_classes": _unique(source_classes),
        "error_classes": _unique(error_classes),
        "uncertainty": bool(uncertain),
        "error_only": bool(error_only),
        "symmetric_error": bool(error_only),
    }


def _term_list(c: Any, raw: Any, label: str, default_id: str, *, require_powers: bool = True) -> list[dict[str, Any]]:
    if isinstance(raw, Mapping) and _term_like(raw):
        values = [raw]
    elif _is_sequence(raw):
        values = list(raw)
    else:
        values = [raw]
    if not values:
        raise ValueError(f"{label} must not be empty")
    return [
        _term(c, value, f"{label}[{index}]", f"{default_id}:{index}", require_powers=require_powers)
        for index, value in enumerate(values)
    ]


def _rows(c: Any, raw: Any, label: str, *, require_powers: bool = True) -> dict[int, list[dict[str, Any]]]:
    if isinstance(raw, Mapping):
        for wrapper in ("ordinary_Z", "source_Z", "jets", "rows", "inlet"):
            if wrapper in raw and isinstance(raw[wrapper], (Mapping, Sequence)) and not _term_like(raw):
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
            raise ValueError(f"{label} requires explicit Z0 and Z1 rows")
        return {
            index: _term_list(c, found[index], f"{label}.Z{index}", f"{label}:Z{index}", require_powers=require_powers)
            for index in ORDERS
        }
    if _is_sequence(raw) and len(raw) == 2:
        return {
            index: _term_list(c, raw[index], f"{label}.Z{index}", f"{label}:Z{index}", require_powers=require_powers)
            for index in ORDERS
        }
    raise ValueError(f"{label} requires explicit Z0 and Z1 rows")


def _normalize_sources(c: Any, raw: Any) -> dict[str, dict[int, list[dict[str, Any]]]]:
    if raw is None:
        raise TypeError("inlet/source rows are required")
    if not isinstance(raw, Mapping):
        raise TypeError("inlet/source rows must be a mapping")
    for wrapper in ("inlet", "sources", "input", "branches"):
        if wrapper in raw and isinstance(raw[wrapper], Mapping) and not any(
            name in raw for name in _SOURCE_ALIASES["F1"] + _SOURCE_ALIASES["Fx"]
        ):
            raw = raw[wrapper]
            break
    result: dict[str, dict[int, list[dict[str, Any]]]] = {}
    for component in INLET_COMPONENTS:
        try:
            value = _lookup(raw, _SOURCE_ALIASES[component], f"inlet.{component}")
        except KeyError:
            branch_name = "F_from_Fc" if component == "F1" else "F_from_Fxc"
            value = _lookup(raw, (branch_name,), f"inlet.{component}")
        result[component] = _rows(c, value, f"inlet.{component}")
    return result


def _normalize_chi(c: Any, chi: Any, chi_Z: Any) -> dict[int, list[dict[str, Any]]]:
    if chi is None:
        raise TypeError("chi is required")
    raw = chi
    if isinstance(raw, Mapping):
        for wrapper in ("chi_jets", "ordinary_Z", "jets", "rows"):
            if wrapper in raw and isinstance(raw[wrapper], (Mapping, Sequence)) and not _term_like(raw):
                raw = raw[wrapper]
                break
    if chi_Z is not None:
        if isinstance(raw, Mapping) and any(_order(key) == 1 for key in raw):
            raise ValueError("chi_Z supplied twice")
        if isinstance(raw, Mapping) and not _term_like(raw):
            z0 = next((raw[key] for key in ("Z0", "z0", 0, "0") if key in raw), None)
            if z0 is None:
                raise ValueError("chi Z0 is required")
        else:
            z0 = raw
        z1 = chi_Z
    else:
        rows = _rows(c, raw, "chi", require_powers=True)
        z0, z1 = rows[0], rows[1]
    rows = {
        0: _term_list(c, z0, "chi.Z0", "chi:Z0", require_powers=True),
        1: _term_list(c, z1, "chi.Z1", "chi:Z1", require_powers=True),
    }
    if len(rows[0]) != 1 or len(rows[1]) != 1:
        raise ValueError("chi and chi_Z each require exactly one source-bound interval")
    if any(tuple(term["powers"]) != (0, 0, 0, 0) for terms in rows.values() for term in terms):
        raise ValueError("chi and chi_Z formal powers must be zero in the tail majorant")
    return rows


def _validate_context(c: Any) -> None:
    if not all(hasattr(c, name) for name in ("mpf", "make_mpf", "exp")):
        raise TypeError("c must be a caller-owned directed interval context")


def _validate_hb(c: Any, hb: Any) -> Any:
    hb = _iv(c, hb, "hb")
    if _lo(c, hb) < 0 or _hi(c, hb) > c.mpf(1) / 2:
        raise ValueError("hb must be enclosed in [0,1/2]")
    return hb


def _sum_powers(*powers: Sequence[int]) -> tuple[int, int, int, int]:
    return tuple(sum(int(item[index]) for item in powers) for index in range(4))  # type: ignore[return-value]


def _scaled_powers(power: Sequence[int], factor: int) -> tuple[int, int, int, int]:
    return tuple(int(item) * factor for item in power)  # type: ignore[return-value]


def _source_radius(c: Any, source: Mapping[str, Any]) -> Any:
    return _nonnegative(c, _abs_upper(c, source["coefficient"]), "input source radius")


def _chi_radius(c: Any, row: Mapping[str, Any]) -> Any:
    return _nonnegative(c, _abs_upper(c, row["coefficient"]), "chi radius")


def _q_interval(c: Any, q: Any, label: str = "q") -> Any:
    q = _iv(c, q, label)
    if _lo(c, q) < 0:
        raise ValueError(f"{label} must be nonnegative")
    if _hi(c, q) > 1:
        raise ValueError(f"{label} must be enclosed in [0,1]")
    if _lo(c, q) > _hi(c, q):
        raise ValueError(f"{label} endpoints are reversed")
    return q


def _normalise_q_cells(c: Any, q_cells: Any) -> list[Any]:
    if isinstance(q_cells, Mapping):
        values = [value for _, value in sorted(q_cells.items(), key=lambda item: int(item[0]))]
    elif _is_sequence(q_cells):
        values = list(q_cells)
    else:
        raise TypeError("q_cells must be an ordered sequence of interval pairs")
    if not values:
        raise ValueError("q_cells must not be empty")
    result = [_q_interval(c, value, f"q_cells[{index}]") for index, value in enumerate(values)]
    for index in range(1, len(result)):
        if _lo(c, result[index - 1]) > _lo(c, result[index]):
            raise ValueError("q_cells must be ordered")
    return result


def _geometry(c: Any, *, hb: Any, sc: Any, q: Any) -> dict[str, Any]:
    a = c.mpf(2)
    x = a * c.exp(hb * sc * q / 8)
    L = _nonnegative(c, x - a, "L")
    scq = sc * q
    B = _nonnegative(c, (scq / 4) * c.exp(scq / 16), "B")
    d = _nonnegative(c, hb * scq / 4, "d")
    return {
        "a": a,
        "x": x,
        "L": L,
        "B": B,
        "d": d,
        "q_interval": q,
        "x_formula": "2*exp(hb*sc*q/8)",
        "B_formula": "(sc*q/4)*exp(sc*q/16)",
        "L_bound": "L <= hb*B",
        "kernel_bounds": {
            "K0": "0 <= K0(x,t) <= x-t",
            "d_x_K0": "0 <= d_x K0(x,t) <= 1",
        },
        "formal_hb_not_materialized": True,
        "hb_upper_used_for_exponent": "1/2",
        "ordinary_Z_geometry_differentiated": False,
        "physical_radial_derivative_preserved": True,
        "midpoint_values_used": False,
    }


def kernel_bounds(
    c: Any,
    *,
    hb: Any,
    sc: Any,
    q: Any,
    chi: Any,
    chi_Z: Any = None,
) -> dict[str, Any]:
    """Return coefficient-only directed majorants for one q interval."""
    _validate_context(c)
    hb_i = _validate_hb(c, hb)
    sc_i = _positive(c, sc, "sc")
    q_i = _q_interval(c, q)
    chi_rows = _normalize_chi(c, chi, chi_Z)
    geometry = _geometry(c, hb=hb_i, sc=sc_i, q=q_i)
    M = _chi_radius(c, chi_rows[0][0])
    Mz = _chi_radius(c, chi_rows[1][0])
    B_upper = _hi(c, geometry["B"])
    M_upper = _hi(c, M)
    Mz_upper = _hi(c, Mz)
    half = c.mpf(1) / 2
    exponent_F = _hi(c, c.exp(M_upper * B_upper * B_upper * half))
    exponent_Z1 = _hi(c, c.exp(M_upper * B_upper * B_upper))
    two_M_sq = (2 * M_upper) ** 2
    four_M_Mz = 4 * M_upper * Mz_upper
    base_F = two_M_sq * (B_upper ** 4) * exponent_F / 24
    base_Fx = two_M_sq * (B_upper ** 3) * exponent_F / 6
    base_F_Z1_chi = four_M_Mz * (B_upper ** 4) * 2 * exponent_Z1 / 24
    base_Fx_Z1_chi = four_M_Mz * (B_upper ** 3) * 2 * exponent_Z1 / 6
    factors = {
        "F_from_Fc": _nonnegative(c, base_F, "F tail factor"),
        "F_from_Fxc": _nonnegative(c, base_F * B_upper, "Fxc-to-F tail factor"),
        "Fx_from_Fc": _nonnegative(c, base_Fx, "Fx tail factor"),
        "Fx_from_Fxc": _nonnegative(c, base_Fx * B_upper, "Fxc-to-Fx tail factor"),
    }
    derivative_factors = {
        "F_from_Fc": _nonnegative(c, base_F_Z1_chi, "F chi_Z tail factor"),
        "F_from_Fxc": _nonnegative(c, base_F_Z1_chi * B_upper, "Fxc chi_Z tail factor"),
        "Fx_from_Fc": _nonnegative(c, base_Fx_Z1_chi, "Fx chi_Z tail factor"),
        "Fx_from_Fxc": _nonnegative(c, base_Fx_Z1_chi * B_upper, "Fxc chi_Z tail factor"),
    }
    if _exact_zero(c, q_i):
        factors = {name: _zero(c) for name in BRANCHES}
        derivative_factors = {name: _zero(c) for name in BRANCHES}
    return {
        "hb": hb_i,
        "sc": sc_i,
        "q": q_i,
        "M": M,
        "Mz": Mz,
        "B": geometry["B"],
        "L": geometry["L"],
        "x": geometry["x"],
        "d": geometry["d"],
        "B_upper": B_upper,
        "M_upper": M_upper,
        "Mz_upper": Mz_upper,
        "exponent_F": exponent_F,
        "exponent_Z1_chi": exponent_Z1,
        "factors": factors,
        "derivative_factors": derivative_factors,
        "formulas": {
            "F": "(2*M)^2*B^4*exp(M*B^2/2)/24",
            "F_from_Fxc": "B*(2*M)^2*B^4*exp(M*B^2/2)/24",
            "Fx": "(2*M)^2*B^3*exp(M*B^2/2)/6",
            "Fx_from_Fxc": "B*(2*M)^2*B^3*exp(M*B^2/2)/6",
            "Z1_chi_F": "(4*M*Mz)*B^4*2*exp(M*B^2)/24",
            "Z1_chi_Fx": "(4*M*Mz)*B^3*2*exp(M*B^2)/6",
        },
        "tail_starts_at_j": 2,
        "hb_upper_for_exponent": half,
        "formal_hb_not_materialized": True,
        "q_zero_exact": _exact_zero(c, q_i),
    }


tail_kernel_bounds = kernel_bounds


def _tail_term(
    c: Any,
    *,
    source: Mapping[str, Any],
    chi0: Mapping[str, Any],
    chi1: Mapping[str, Any],
    branch: str,
    source_order: int,
    chi_order: int,
    factor: Any,
    kernel: Mapping[str, Any],
) -> dict[str, Any]:
    derivative_piece = chi_order == 1
    if derivative_piece:
        chi_terms = (chi1, chi0)
        chi_powers = _sum_powers(chi1["powers"], chi0["powers"])
        leibniz_piece = "chi_Z*inlet_Z0"
    else:
        chi_terms = (chi0, chi0)
        chi_powers = _scaled_powers(chi0["powers"], 2)
        leibniz_piece = "chi*inlet_Z1" if source_order == 1 else "chi*inlet_Z0"
    powers = list(_sum_powers(source["powers"], chi_powers))
    powers[0] += FORMAL_HB_POWER_ADDED[branch]
    radius = _nonnegative(c, _source_radius(c, source) * factor, "j2 tail radius")
    coefficient = _symmetric(c, radius)
    input_ids = list(source.get("source_ids", ()))
    chi_ids_raw: list[Any] = []
    chi_source_classes: list[Any] = []
    chi_error_classes: list[Any] = []
    for term in chi_terms:
        chi_ids_raw.extend(_items(term.get("source_ids")))
        chi_source_classes.extend(_items(term.get("source_error_classes")))
        chi_error_classes.extend(_items(term.get("error_classes")))
    source_ids = _unique(input_ids + chi_ids_raw)
    source_classes = _unique(_items(source.get("source_error_classes")) + chi_source_classes)
    error_classes = _unique(
        _items(source.get("error_classes"))
        + chi_error_classes
        + ["j2_kernel_remainder", "symmetric_error"]
    )
    field = "F1" if branch.startswith("F_") else "d_x_F1"
    input_component = "F1" if branch.endswith("Fc") else "Fx"
    return {
        "coefficient": coefficient,
        "coefficient_is_symmetric_radius": True,
        "radius": radius,
        "radius_factor": factor,
        "powers": tuple(powers),
        "input_powers": tuple(source["powers"]),
        "chi_powers": tuple(chi_powers),
        "formal_hb_power_added": FORMAL_HB_POWER_ADDED[branch],
        "coefficient_contains_formal_hb": False,
        "source_ids": source_ids,
        "input_source_ids": input_ids,
        "incoming_source_ids": input_ids,
        "chi_source_ids": chi_ids_raw,
        "chi_source_ids_unique": _unique(chi_ids_raw),
        "source_error_classes": source_classes,
        "input_error_classes": list(source.get("error_classes", ())),
        "chi_source_error_classes": _unique(chi_source_classes),
        "error_classes": error_classes,
        "uncertainty": True,
        "error_only": True,
        "symmetric_error": True,
        "signed_source_term": False,
        "signed_tail_value": False,
        "signed_input_coefficient": source["coefficient"],
        "signed_input_coefficient_preserved": True,
        "branch": branch,
        "field": field,
        "canonical_field": "F" if branch.startswith("F_") else "Fx",
        "input_component": input_component,
        "source_component": input_component,
        "incoming_component": input_component,
        "source_order": source_order,
        "input_order": source_order,
        "chi_order": chi_order,
        "output_order": 1 if derivative_piece or source_order == 1 else 0,
        "leibniz_piece": leibniz_piece,
        "tail_starts_at_j": 2,
        "kernel_tail_order": "j>=2",
        "operator_word_coefficient": "(-2)^j",
        "absolute_majorant_drops_alternating_sign": True,
        "initial_K0_attached_once": True,
        "baseline_K0_included": False,
        "j1_included": False,
        "physical_dx_count_minimum": 4 if branch.startswith("F_") else 3,
        "formal_hb_power_source": "L^4 <= hb^4 B^4" if branch.startswith("F_") else "L^3 <= hb^3 B^3",
        "h_fxc_bound": "(1-exp(-d))*|Fxc| <= hb*B*|Fxc|" if branch.endswith("Fxc") else None,
        "chi_factor_count": 2,
        "chi_derivative_j_factor_included": derivative_piece,
        "chi_derivative_extra_factor": 2 if derivative_piece else None,
        "q_zero_exact": bool(kernel.get("q_zero_exact")),
        "origin": {
            "operation": "n1_e00_collar_j2_volterra_v2",
            "branch": branch,
            "tail_starts_at_j": 2,
            "leibniz_piece": leibniz_piece,
            "input_source_ids": input_ids,
            "chi_source_ids": chi_ids_raw,
            "source_ids": source_ids,
            "error_classes": error_classes,
            "formal_hb_power_added": FORMAL_HB_POWER_ADDED[branch],
            "symmetric_error": True,
            "no_ledger": True,
        },
    }


def _branch(
    c: Any,
    *,
    branch: str,
    sources: Mapping[str, Mapping[int, Sequence[Mapping[str, Any]]]],
    chi_rows: Mapping[int, Sequence[Mapping[str, Any]]],
    kernel: Mapping[str, Any],
) -> dict[str, Any]:
    component = "F1" if branch.endswith("Fc") else "Fx"
    rows = {"Z0": [], "Z1": []}
    source_rows = sources[component]
    factors = kernel["factors"]
    derivative_factors = kernel["derivative_factors"]
    for source in source_rows[0]:
        rows["Z0"].append(
            _tail_term(
                c, source=source, chi0=chi_rows[0][0], chi1=chi_rows[1][0],
                branch=branch, source_order=0, chi_order=0,
                factor=factors[branch], kernel=kernel,
            )
        )
    for source in source_rows[1]:
        rows["Z1"].append(
            _tail_term(
                c, source=source, chi0=chi_rows[0][0], chi1=chi_rows[1][0],
                branch=branch, source_order=1, chi_order=0,
                factor=factors[branch], kernel=kernel,
            )
        )
    for source in source_rows[0]:
        rows["Z1"].append(
            _tail_term(
                c, source=source, chi0=chi_rows[0][0], chi1=chi_rows[1][0],
                branch=branch, source_order=0, chi_order=1,
                factor=derivative_factors[branch], kernel=kernel,
            )
        )
    return {
        "branch": branch,
        "field": "F1" if branch.startswith("F_") else "d_x_F1",
        "canonical_field": "F" if branch.startswith("F_") else "Fx",
        "input_component": component,
        "source_component": component,
        "Z0": rows["Z0"],
        "Z1": rows["Z1"],
        "response": rows,
        "tail_starts_at_j": 2,
        "ordinary_Z1_product_rule": "chi*inlet_Z1 + chi_Z*inlet_Z0",
        "chi_derivative_j_factor_included": True,
        "signed_sectors": False,
        "symmetric_error_outputs": True,
        "baseline_K0_included": False,
        "j1_included": False,
        "no_ledger": True,
        "formal_hb_power_added": FORMAL_HB_POWER_ADDED[branch],
    }


def _atoms(source_atoms: Any) -> list[Any]:
    if source_atoms is None:
        return list(SOURCE_POWER_BASIS)
    values = list(source_atoms.keys()) if isinstance(source_atoms, Mapping) else list(source_atoms)
    if set(values) != set(SOURCE_POWER_BASIS):
        raise ValueError("source_atoms must contain exactly the four formal source atoms")
    return values


def _build_one(
    c: Any,
    *,
    hb: Any,
    sc: Any,
    q: Any,
    sources: Mapping[str, Mapping[int, Sequence[Mapping[str, Any]]]],
    chi_rows: Mapping[int, Sequence[Mapping[str, Any]]],
    source_family: Any,
    source_atoms: Any,
    certificate: Any,
    q_index: int | None,
) -> dict[str, Any]:
    kernel = kernel_bounds(c, hb=hb, sc=sc, q=q, chi=chi_rows)
    branches = {
        branch: _branch(c, branch=branch, sources=sources, chi_rows=chi_rows, kernel=kernel)
        for branch in BRANCHES
    }
    seam = bool(kernel["q_zero_exact"])
    geometry = _geometry(c, hb=kernel["hb"], sc=kernel["sc"], q=kernel["q"])
    return {
        "schema": SCHEMA,
        "status": "conditional_fixed_e00_j2_volterra_tail",
        "source_family": source_family,
        "source_atoms": _atoms(source_atoms),
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "fields": list(FIELDS),
        "canonical_fields": list(CANONICAL_FIELDS),
        "ordinary_Z_orders": [0, 1],
        "branches": branches,
        "branch_names": list(BRANCHES),
        "response": {name: branches[name]["response"] for name in BRANCHES},
        "tail": {name: branches[name]["response"] for name in BRANCHES},
        "inlet": {
            component: {f"Z{order}": list(sources[component][order]) for order in ORDERS}
            for component in INLET_COMPONENTS
        },
        "chi_jets": {f"Z{order}": list(chi_rows[order]) for order in ORDERS},
        "geometry": geometry,
        "kernel": kernel,
        "kernel_bounds": kernel,
        "chi_Z_required": True,
        "ordinary_Z1_product_rule": "chi*inlet_Z1 + chi_Z*inlet_Z0",
        "ordinary_Z_geometry_differentiated": False,
        "fixed_geometry_Z_independent": True,
        "physical_radial_derivative_preserved": True,
        "formal_hb_powers_added": dict(FORMAL_HB_POWER_ADDED),
        "source_ids_preserved": True,
        "chi_ids_preserved": True,
        "source_error_classes_preserved": True,
        "error_classes_preserved": True,
        "formal_powers_preserved": True,
        "symmetric_error_outputs": True,
        "error_only_outputs_explicit": True,
        "j2_remainder_included": True,
        "j2_kernel_remainder_included": True,
        "allword_tail_added": False,
        "first_potential_included": False,
        "baseline_K0_included": False,
        "j1_included": False,
        "appended_to_ledger": False,
        "ledger_written": False,
        "no_ledger": True,
        "q_interval": kernel["q"],
        "q_index": q_index,
        "seam": {
            "q": kernel["q"],
            "exact_zero": seam,
            "componentwise_zero": seam,
            "all_tail_coefficients_zero": seam,
            "higher_potential_words_exact_zero": seam,
            "literal_directed_arithmetic": True,
        },
        "certificate_passthrough": certificate,
        "old_builder_called": False,
        "ancestor_imports": False,
        "receipt_reads": False,
        "midpoint_values_used": False,
        "automatic_code_execution": False,
    }

def build(
    c: Any,
    inlet: Any = None,
    chi: Any = None,
    chi_Z: Any = None,
    *,
    hb: Any = None,
    sc: Any = None,
    q: Any = None,
    q_interval: Any = None,
    q_cells: Any = None,
    source: Any = None,
    sources: Any = None,
    source_family: Any = None,
    source_atoms: Any = None,
    certificate: Any = None,
    geometry_certificate: Any = None,
) -> dict[str, Any]:
    """Build a fixed-E00 ``j >= 2`` tail enclosure."""
    _validate_context(c)
    if hb is None or sc is None:
        raise TypeError("hb and sc interval inputs are required")
    if geometry_certificate is not None:
        if certificate is not None and geometry_certificate != certificate:
            raise ValueError("certificate and geometry_certificate disagree")
        certificate = geometry_certificate
    aliases = [value is not None for value in (inlet, source, sources)]
    if sum(aliases) > 1:
        raise ValueError("inlet, source and sources are aliases; supply only one")
    raw_source = inlet if inlet is not None else source if source is not None else sources
    source_rows = _normalize_sources(c, raw_source)
    chi_rows = _normalize_chi(c, chi, chi_Z)
    _atoms(source_atoms)
    _validate_hb(c, hb)
    _positive(c, sc, "sc")
    if q is not None and q_interval is not None:
        raise ValueError("q and q_interval are aliases; supply only one")
    if q is None:
        q = q_interval
    if q_cells is not None and q is not None:
        raise ValueError("q and q_cells are mutually exclusive")
    if q_cells is not None:
        intervals = _normalise_q_cells(c, q_cells)
    elif q is not None:
        intervals = [_q_interval(c, q)]
    else:
        raise TypeError("q or q_cells interval input is required")
    cells = [
        _build_one(
            c,
            hb=hb,
            sc=sc,
            q=interval,
            sources=source_rows,
            chi_rows=chi_rows,
            source_family=source_family,
            source_atoms=source_atoms,
            certificate=certificate,
            q_index=index if q_cells is not None else None,
        )
        for index, interval in enumerate(intervals)
    ]
    if len(cells) == 1:
        result = dict(cells[0])
        result["q_cells"] = cells
        result["cell_count"] = 1
        return result
    return {
        "schema": SCHEMA,
        "status": "conditional_fixed_e00_j2_volterra_tail_cells",
        "source_family": source_family,
        "source_atoms": _atoms(source_atoms),
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "fields": list(FIELDS),
        "canonical_fields": list(CANONICAL_FIELDS),
        "ordinary_Z_orders": [0, 1],
        "q_cells": cells,
        "cell_count": len(cells),
        "q_interval": None,
        "branches": {},
        "response": {},
        "tail": {},
        "inlet": {
            component: {f"Z{order}": list(source_rows[component][order]) for order in ORDERS}
            for component in INLET_COMPONENTS
        },
        "chi_jets": {f"Z{order}": list(chi_rows[order]) for order in ORDERS},
        "formal_hb_powers_added": dict(FORMAL_HB_POWER_ADDED),
        "j2_remainder_included": True,
        "symmetric_error_outputs": True,
        "appended_to_ledger": False,
        "ledger_written": False,
        "no_ledger": True,
        "midpoint_values_used": False,
        "automatic_code_execution": False,
    }


build_tail = build
build_j2_tail = build
evaluate = build


