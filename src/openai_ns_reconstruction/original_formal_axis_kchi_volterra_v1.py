"""Formal-Hb axis-Kchi Volterra words on the complete collar.

This leaf is intentionally source-shaped rather than a numerical solver.  It
keeps ``hb`` as a formal power and only evaluates caller-owned directed
intervals for the coefficient factors.  The two input ordinary-Z rows are
copied independently, every finite word is signed, and omitted words are
returned as separate symmetric error records.

The scope is deliberately narrow: the formal ``F``/``Fx`` sector only.  It
does not construct the other four fields, recurse in time, read a receipt or
saved data, append a ledger, claim an m15 result, or prove a same-owner
identity.  The caller supplies ``sc`` and whole-Z axis ``chi``/``chi_z``
intervals; no lower-collar shear is accepted as a substitute.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "original_formal_axis_kchi_volterra_v1"
COMPONENT = "formal_Hb_axis_Kchi_volterra"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
ORDINARY_Z_LABELS = ("Z0", "Z1")
FIELDS = ("F", "Fx")
BASE_NAMES = ("Fc", "Fxc", "forcing")
MAX_ORDER = 32
MAX_INPUT_TERMS = 256
MAX_OUTPUT_TERMS = 4096

_MISSING = object()
_COEFFICIENT_KEYS = (
    "coefficient",
    "coefficient_interval",
    "interval",
    "enclosure",
    "box",
    "value",
)
_POWER_KEYS = ("powers", "formal_powers", "source_powers")

__all__ = [
    "SCHEMA", "COMPONENT", "SOURCE_POWER_BASIS", "ORDINARY_Z_ORDERS",
    "ORDINARY_Z_LABELS", "FIELDS", "BASE_NAMES", "build",
    "build_formal_axis_kchi_volterra", "build_tail", "enclose", "evaluate",
]


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _items(value: Any) -> list[Any]:
    if value is None or value is _MISSING:
        return []
    return list(value) if _is_sequence(value) else [value]


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    seen: set[tuple[str, str]] = set()
    for value in values:
        key = (type(value).__name__, repr(value))
        if key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _lookup(mapping: Mapping[Any, Any], names: Sequence[Any], label: str) -> Any:
    for name in names:
        try:
            if name in mapping:
                return mapping[name]
        except (TypeError, ValueError):
            continue
    lowered = {str(key).lower(): value for key, value in mapping.items()}
    for name in names:
        if str(name).lower() in lowered:
            return lowered[str(name).lower()]
    if label:
        raise KeyError(f"{label} is required")
    return _MISSING


def _raw(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    raw = getattr(value, "_mpf_", None)
    if raw is not None:
        return raw
    if hasattr(value, "a") and hasattr(value, "b"):
        endpoint = value.b if upper else value.a
        nested = getattr(endpoint, "_mpi_", None)
        return nested[1 if upper else 0] if nested is not None else endpoint
    return value


def _make(c: Any, lower: Any, upper: Any) -> Any:
    return c.make_mpf((_raw(lower, False), _raw(upper, True)))


def _iv(c: Any, value: Any, label: str = "interval") -> Any:
    """Coerce a value to a directed interval without midpointing."""

    if hasattr(value, "_mpi_"):
        return c.make_mpf((_raw(value, False), _raw(value, True)))
    if _is_mapping(value):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            if "lower_exact_mpf_tuple" not in value or "upper_exact_mpf_tuple" not in value:
                raise TypeError(f"{label} requires both exact MPF endpoint tuples")
            lower = _exact_mpf_tuple(c, value["lower_exact_mpf_tuple"], label + ".lower_exact_mpf_tuple")
            upper = _exact_mpf_tuple(c, value["upper_exact_mpf_tuple"], label + ".upper_exact_mpf_tuple")
            try:
                interval = c.make_mpf((lower, upper))
            except Exception as exc:
                raise TypeError(f"{label} exact MPF endpoints are malformed") from exc
            if _lo(c, interval) > _hi(c, interval):
                raise ValueError(f"{label} exact MPF endpoints are reversed")
            return interval
        nested = _lookup(value, ("interval", "enclosure", "box", "coefficient", "coefficient_interval", "value"), "")
        if nested is not _MISSING and nested is not value:
            return _iv(c, nested, label)
        lower = _lookup(value, ("lower", "lo", "left"), "")
        upper = _lookup(value, ("upper", "hi", "right"), "")
        if lower is not _MISSING and upper is not _MISSING:
            return _make(c, lower, upper)
        raise TypeError(f"{label} is not interval-compatible")
    if _is_sequence(value) and len(value) == 2:
        return _make(c, value[0], value[1])
    try:
        return c.mpf(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{label} is not interval-compatible") from exc


def _exact_mpf_tuple(c: Any, value: Any, label: str) -> tuple[int, int, int, int]:
    if not isinstance(value, (tuple, list)) or len(value) != 4:
        raise TypeError(f"{label} must be a four-integer MPF tuple")
    if any(type(item) is not int for item in value):
        raise TypeError(f"{label} may not contain decimal strings or floats")
    result = tuple(int(item) for item in value)
    sign, mantissa, _exponent, bitcount = result
    if sign not in (0, 1) or mantissa < 0 or bitcount < 0:
        raise ValueError(f"{label} has invalid MPF sign/mantissa/bitcount")
    if mantissa == 0 and bitcount != 0:
        raise ValueError(f"{label} is a non-finite MPF special value")
    try:
        point = c.make_mpf((result, result))
        isfinite = getattr(c, "isfinite", None)
        if callable(isfinite) and not bool(isfinite(point)):
            raise ValueError(f"{label} must be finite")
        isinf = getattr(c, "isinf", None)
        isnan = getattr(c, "isnan", None)
        if callable(isinf) and bool(isinf(point)):
            raise ValueError(f"{label} must be finite")
        if callable(isnan) and bool(isnan(point)):
            raise ValueError(f"{label} must be finite")
    except ValueError:
        raise
    except Exception as exc:
        raise TypeError(f"{label} is not a valid finite MPF tuple") from exc
    return result


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


def _one(c: Any) -> Any:
    return c.mpf(1)


def _range(c: Any, lower: Any, upper: Any) -> Any:
    lower, upper = _lo(c, lower), _hi(c, upper)
    if upper < lower:
        raise ValueError("directed interval endpoints are reversed")
    return c.make_mpf((_raw(lower, False), _raw(upper, True)))


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _iv(c, value, label)
    lower, upper = _lo(c, value), _hi(c, value)
    if upper < 0:
        raise ValueError(f"{label} must have nonnegative upper endpoint")
    if lower < 0:
        lower = _zero(c)
    return _range(c, lower, upper)


def _positive(c: Any, value: Any, label: str) -> Any:
    value = _iv(c, value, label)
    if _lo(c, value) <= 0:
        raise ValueError(f"{label} must be strictly positive")
    return value


def _abs_upper(c: Any, value: Any) -> Any:
    lower, upper = _lo(c, value), _hi(c, value)
    if lower >= 0:
        return upper
    if upper <= 0:
        return -lower
    return max(-lower, upper)


def _symmetric(c: Any, radius: Any) -> Any:
    upper = _hi(c, _nonnegative(c, radius, "symmetric radius"))
    return c.make_mpf((_raw(-upper, False), _raw(upper, True)))


def _factorial(c: Any, n: int) -> Any:
    result = _one(c)
    for index in range(2, n + 1):
        result = result * index
    return result


def _power_interval(c: Any, value: Any, exponent: int) -> Any:
    if exponent < 0:
        raise ValueError("interval powers must be nonnegative")
    result = _one(c)
    for _ in range(exponent):
        result = result * value
    return result


def _add_powers(*powers: Sequence[int]) -> list[int]:
    return [sum(int(item[index]) for item in powers) for index in range(4)]


def _scaled_powers(power: Sequence[int], multiplier: int) -> list[int]:
    return [int(item) * multiplier for item in power]


def _rows(c: Any, raw: Any, label: str, *, allow_empty: bool = False) -> dict[int, list[dict[str, Any]]]:
    current = raw
    if _is_mapping(current):
        for wrapper in ("ordinary_Z", "source_Z", "rows", "jets", "terms"):
            candidate = current.get(wrapper, _MISSING)
            if candidate is not _MISSING and candidate is not current:
                current = candidate
                break
    if _is_mapping(current):
        found: dict[int, Any] = {}
        for key, value in current.items():
            text = str(key).lower()
            index = (
                0 if key == 0 or text in {"0", "z0", "ordinary_z0", "value"}
                else 1 if key == 1 or text in {"1", "z1", "ordinary_z1", "derivative"}
                else None
            )
            if index is not None:
                if index in found:
                    raise ValueError(f"{label} supplies Z{index} twice")
                found[index] = value
        if set(found) != {0, 1}:
            raise ValueError(f"{label} requires explicit Z0 and Z1 rows")
        return {index: _term_list(c, found[index], f"{label}.Z{index}", allow_empty=allow_empty) for index in ORDINARY_Z_ORDERS}
    if _is_sequence(current) and len(current) == 2:
        return {index: _term_list(c, current[index], f"{label}.Z{index}", allow_empty=allow_empty) for index in ORDINARY_Z_ORDERS}
    raise TypeError(f"{label} requires ordinary Z0 and Z1 term rows")


def _term_list(c: Any, raw: Any, label: str, *, allow_empty: bool = False) -> list[dict[str, Any]]:
    if _is_mapping(raw) and any(key in raw for key in _COEFFICIENT_KEYS):
        values = [raw]
    elif _is_sequence(raw):
        values = list(raw)
    else:
        raise TypeError(f"{label} must be a term list")
    if not values and not allow_empty:
        raise ValueError(f"{label} must not be empty")
    return [_term(c, value, f"{label}[{index}]") for index, value in enumerate(values)]


def _normalise_powers(raw: Mapping[Any, Any], label: str) -> list[int]:
    powers = next((raw[key] for key in _POWER_KEYS if key in raw), _MISSING)
    if powers is _MISSING:
        powers = (
            raw.get("hb_power", _MISSING),
            raw.get("F0_power", raw.get("F0_anchor_power", _MISSING)),
            raw.get("Pstar_squared_power", raw.get("Pstar_power", _MISSING)),
            raw.get("Lambda_power", raw.get("lambda_power", _MISSING)),
        )
    if _is_mapping(powers):
        powers = [powers.get(atom, _MISSING) for atom in SOURCE_POWER_BASIS]
    if not _is_sequence(powers) or len(powers) != 4:
        raise TypeError(f"{label}.powers must contain four integer formal powers")
    if any(item is _MISSING or isinstance(item, bool) or not isinstance(item, int) for item in powers):
        raise TypeError(f"{label}.powers must contain four integer formal powers")
    return [int(item) for item in powers]


def _metadata(raw: Mapping[Any, Any], label: str) -> dict[str, Any]:
    ids: list[Any] = []
    for key in ("source_ids", "source_id", "sources", "source"):
        if key in raw:
            ids.extend(_items(raw[key]))
    if not ids and _is_mapping(raw.get("origin")):
        origin = raw["origin"]
        for key in ("source_ids", "source_id", "sources", "source"):
            if key in origin:
                ids.extend(_items(origin[key]))
    if not ids:
        raise ValueError(f"{label} must retain at least one source_id")
    source_paths = _items(raw.get("source_paths", raw.get("source_path", [])))
    source_error_classes = _items(raw.get("source_error_classes", raw.get("source_error_class", [])))
    error_classes = _items(raw.get("error_classes", raw.get("error_class", [])))
    if not error_classes:
        error_classes = list(source_error_classes)
    error_only = raw.get("error_only", False)
    symmetric = raw.get("symmetric_error", error_only)
    if not isinstance(error_only, bool) or not isinstance(symmetric, bool):
        raise TypeError(f"{label}.error_only and symmetric_error must be boolean")
    if error_only and not symmetric:
        raise ValueError(f"{label} error_only requires symmetric_error")
    return {
        "source_ids": _unique(ids),
        "source_paths": _unique(source_paths),
        "source_error_classes": _unique(source_error_classes),
        "error_classes": _unique(error_classes),
        "error_only": error_only,
        "symmetric_error": symmetric,
        "uncertainty": bool(raw.get("uncertainty", raw.get("is_uncertainty", False)) or error_only),
    }


def _term(c: Any, raw: Any, label: str) -> dict[str, Any]:
    if not _is_mapping(raw):
        raise TypeError(f"{label} must be a term mapping")
    coefficient = next((raw[key] for key in _COEFFICIENT_KEYS if key in raw), _MISSING)
    if coefficient is _MISSING:
        raise KeyError(f"{label} is missing coefficient")
    metadata = _metadata(raw, label)
    coefficient_interval = _iv(c, coefficient, f"{label}.coefficient")
    if metadata["error_only"]:
        lower, upper = coefficient_interval._mpi_
        if not (lower[1] == 0 and upper[1] == 0):
            if lower[0] != 1 or upper[0] != 0 or lower[1:] != upper[1:]:
                raise ValueError(f"{label} error_only requires symmetric binary endpoints")
    powers = _normalise_powers(raw, label)
    result = dict(raw)
    result.update(
        {
            "coefficient": coefficient_interval,
            "coefficient_interval": coefficient_interval,
            "powers": list(powers),
            "formal_powers": list(powers),
            "power_basis": list(SOURCE_POWER_BASIS),
            **metadata,
        }
    )
    return result


def _unwrap(raw: Any, names: Sequence[str]) -> Any:
    current = raw
    for _ in range(4):
        if not _is_mapping(current):
            break
        nested = _MISSING
        for name in names:
            if name in current and current[name] is not current:
                nested = current[name]
                break
        if nested is _MISSING:
            break
        current = nested
    return current


def _normalise_inlet(c: Any, inlet: Any) -> dict[str, dict[int, list[dict[str, Any]]]]:
    if not _is_mapping(inlet):
        raise TypeError("inlet must be a mapping with F and Fx ordinary rows")
    current = _unwrap(inlet, ("inlet", "sources", "input"))
    if not _is_mapping(current):
        raise TypeError("inlet must be a mapping with F and Fx ordinary rows")
    f_raw = next((current[key] for key in ("F", "Fc", "F1") if key in current), _MISSING)
    fx_raw = next((current[key] for key in ("Fx", "Fxc", "d_x_F", "d_x_F1") if key in current), _MISSING)
    if f_raw is _MISSING or fx_raw is _MISSING:
        raise ValueError("inlet must contain both F and Fx rows")
    return {"F": _rows(c, f_raw, "inlet.F"), "Fx": _rows(c, fx_raw, "inlet.Fx")}


def _normalise_forcing(c: Any, forcing: Any) -> dict[int, list[dict[str, Any]]]:
    if not _is_mapping(forcing):
        raise TypeError("forcing must be a mapping with uniform h5 Z0 and Z1 rows")
    current = _unwrap(forcing, ("forcing", "h5", "source", "rows"))
    return _rows(c, current, "forcing", allow_empty=True)


def _normalise_axis_interval(c: Any, raw: Any, label: str) -> dict[str, Any]:
    if not _is_mapping(raw):
        raise TypeError(f"{label} must be a source-backed whole-Z interval mapping")
    axis_raw = dict(raw)
    if not any(key in axis_raw for key in _POWER_KEYS + ("hb_power", "F0_power", "F0_anchor_power", "Pstar_squared_power", "Pstar_power", "Lambda_power", "lambda_power")):
        axis_raw["powers"] = [0, 0, 0, 0]
    term = _term(c, axis_raw, label)
    if term["error_only"]:
        raise ValueError(f"{label} must be an axis-potential enclosure, not an error-only packet")
    if term["powers"] != [0, 0, 0, 0]:
        raise ValueError(f"{label} formal powers must be zero")
    scope = " ".join(str(raw.get(key, "")) for key in ("scope", "domain", "role", "support")).lower()
    if raw.get("whole_Z", True) is False or "lower-collar" in scope or "lower collar" in scope or "shear" in scope:
        raise ValueError(f"{label} must be the whole-Z axis Kchi interval")
    return term


def _validate_context(c: Any) -> None:
    if not all(hasattr(c, name) for name in ("mpf", "make_mpf", "exp", "sqrt")):
        raise TypeError("c must provide directed mpf/make_mpf/exp/sqrt arithmetic")


def _validate_order(max_order: Any) -> int:
    if isinstance(max_order, bool) or not isinstance(max_order, int):
        raise TypeError("max_order must be a nonnegative integer")
    if max_order < 0 or max_order > MAX_ORDER:
        raise ValueError(f"max_order must lie in [0,{MAX_ORDER}]")
    return max_order


def _validate_min_order(min_order: Any, max_order: int) -> int:
    if isinstance(min_order, bool) or not isinstance(min_order, int):
        raise TypeError("min_order must be a nonnegative integer")
    if min_order < 0 or min_order > max_order + 1:
        raise ValueError("min_order must lie between zero and max_order+1")
    return min_order


def _validate_sc(c: Any, sc: Any) -> Any:
    if sc is None:
        raise TypeError("caller-certified positive sc interval is required")
    if _is_mapping(sc):
        for key in ("positive", "certified_positive", "Z_independent", "z_independent"):
            if key in sc and not isinstance(sc[key], bool):
                raise TypeError(f"sc.{key} must be boolean")
            if key in sc and sc[key] is False:
                raise ValueError(f"sc.{key} certificate is false")
    return _positive(c, _iv(c, sc, "sc"), "sc")


def _merge_ids(*groups: Sequence[Any]) -> list[Any]:
    values: list[Any] = []
    for group in groups:
        values.extend(group)
    return _unique(values)


def _merge_classes(*groups: Sequence[Any]) -> list[Any]:
    values: list[Any] = []
    for group in groups:
        values.extend(group)
    return _unique(values)


def _has_shared_id(left: Sequence[Any], right: Sequence[Any]) -> bool:
    return any(any(candidate == other for other in right) for candidate in left)


def _provenance(parent: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "source_ids": list(parent.get("source_ids", [])),
        "source_paths": list(parent.get("source_paths", [])),
        "source_error_classes": list(parent.get("source_error_classes", [])),
        "error_classes": list(parent.get("error_classes", [])),
        "term_path": parent.get("term_path"),
        "origin": parent.get("origin"),
    }


def _output_term(
    parent: Mapping[str, Any],
    *,
    coefficient: Any,
    powers: Sequence[int],
    base_name: str,
    output_field: str,
    z_order: int,
    word_order: int,
    word_kind: str,
    leibniz_piece: str,
    chi_terms: Sequence[Mapping[str, Any]],
    formal_hb_power_added: int,
    tail: bool = False,
    radius_formula: str | None = None,
) -> dict[str, Any]:
    chi_ids = _merge_ids(*(term.get("source_ids", []) for term in chi_terms))
    chi_z_ids = _merge_ids(*(term.get("source_ids", []) for term in chi_terms[1:])) if len(chi_terms) > 1 else []
    chi_classes = _merge_classes(*(term.get("error_classes", []) for term in chi_terms))
    inherited_error = bool(parent.get("error_only", False))
    if tail:
        result_error_only = True
        result_symmetric = True
        kind = "volterra_tail"
    else:
        result_error_only = inherited_error
        result_symmetric = bool(parent.get("symmetric_error", False)) if inherited_error else False
        kind = "inherited_error" if inherited_error else "signed_finite"
    result = {
        "coefficient": coefficient,
        "coefficient_interval": coefficient,
        "coefficient_is_symmetric_radius": bool(tail),
        "signed_parent_coefficient": parent.get("coefficient"),
        "signed_source_coefficient_preserved": True,
        "powers": list(powers),
        "formal_powers": list(powers),
        "power_basis": list(SOURCE_POWER_BASIS),
        "source_ids": _merge_ids(parent.get("source_ids", []), chi_ids),
        "parent_source_ids": list(parent.get("source_ids", [])),
        "chi_source_ids": chi_ids,
        "chi_Z_source_ids": chi_z_ids,
        "source_paths": list(parent.get("source_paths", [])),
        "source_error_classes": _merge_classes(parent.get("source_error_classes", []), chi_classes),
        "error_classes": _merge_classes(parent.get("error_classes", []), chi_classes, ["formal_volterra_tail"] if tail else []),
        "parent_source_error_classes": list(parent.get("source_error_classes", [])),
        "parent_error_classes": list(parent.get("error_classes", [])),
        "uncertainty": bool(tail or parent.get("uncertainty", False) or inherited_error),
        "error_only": result_error_only,
        "symmetric_error": result_symmetric,
        "inherited_error_only": inherited_error,
        "signed_finite": bool(not tail and not inherited_error),
        "signed_tail_value": False,
        "point_value_selected": False,
        "tail": bool(tail),
        "term_kind": kind,
        "base": base_name,
        "field": output_field,
        "ordinary_Z": f"Z{z_order}",
        "ordinary_Z_order": z_order,
        "word_order": word_order,
        "word_kind": word_kind,
        "operator_word_coefficient": "(-2)^n" if word_order else "1",
        "leibniz_piece": leibniz_piece,
        "formal_hb_power_added": formal_hb_power_added,
        "coefficient_contains_formal_hb": False,
        "parent_provenance": _provenance(parent),
        "no_midpoint": True,
    }
    if radius_formula is not None:
        result["radius_formula"] = radius_formula
    return result


def _geometry(c: Any, sc: Any) -> dict[str, Any]:
    q_interval = c.make_mpf((_raw(_zero(c), False), _raw(_one(c), True)))
    S = _nonnegative(c, (sc / 4) * c.exp(sc / 16), "S")
    return {
        "q_interval": q_interval,
        "x_formula": "x(q)=2*exp(hb*sc*q/8), 0<=q<=1",
        "x_end_formula": "x_end=2*exp(hb*sc/8)",
        "L_formula": "L=x_end-2",
        "S": S,
        "S_formula": "S=(sc/4)*exp(sc/16)",
        "L_bound": "L <= hb*S",
        "K0_formula": "K0(x,t)=t*(1-t^2/x^2)/2",
        "K0_bound": "0 <= K0(x,t) <= x-t",
        "Kx_formula": "Kx(x,t)=(t/x)^3",
        "Kx_bound": "0 <= Kx(x,t) <= 1",
        "kernel_bounds_proved_under_formal_condition": True,
        "L_bound_proved_under_formal_condition": True,
        "formal_hb_not_materialized": True,
        "ordinary_Z_geometry_differentiated": False,
        "hb_sc_geometry_Z_independent": True,
        "midpoint_values_used": False,
    }


def _base_specs(c: Any, sc: Any, geometry: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    S = geometry["S"]
    zero, one = _zero(c), _one(c)
    Fxc_F_factor = _range(c, zero, sc / 4)
    Fxc_Fx_factor = _range(c, c.exp(-3 * sc / 16), one)
    forcing_F_factor = _range(c, zero, S ** 2 / 2)
    forcing_Fx_factor = _range(c, zero, S)
    return {
        "Fc": {
            "input": "F", "F_factor": one, "F_hb": 0,
            "Fx_direct_factor": None, "Fx_direct_hb": 0,
            "factor_formula_F": "1", "factor_formula_Fx": "0 at n=0",
        },
        "Fxc": {
            "input": "Fx", "F_factor": Fxc_F_factor, "F_hb": 1,
            "Fx_direct_factor": Fxc_Fx_factor, "Fx_direct_hb": 0,
            "factor_formula_F": "hb*[0,sc/4]", "factor_formula_Fx": "[exp(-3*sc/16),1]",
        },
        "forcing": {
            "input": "forcing", "F_factor": forcing_F_factor, "F_hb": 2,
            "Fx_direct_factor": forcing_Fx_factor, "Fx_direct_hb": 1,
            "factor_formula_F": "hb^2*[0,S^2/2]", "factor_formula_Fx": "hb*[0,S]",
        },
    }


def _word_multiplier(c: Any, chi: Mapping[str, Any], order: int) -> Any:
    return ((-2) ** order) * _power_interval(c, chi["coefficient"], order)


def _insertion_multiplier(c: Any, chi: Mapping[str, Any], chi_z: Mapping[str, Any], order: int) -> Any:
    return ((-2) ** order) * order * _power_interval(c, chi["coefficient"], order - 1) * chi_z["coefficient"]


def _finite_term(
    c: Any,
    parent: Mapping[str, Any],
    *,
    factor: Any,
    factor_hb: int,
    output_field: str,
    base_name: str,
    z_order: int,
    order: int,
    chi: Mapping[str, Any],
    chi_z: Mapping[str, Any],
    S: Any,
    insertion: bool,
) -> dict[str, Any]:
    if insertion:
        multiplier = _insertion_multiplier(c, chi, chi_z, order)
        chi_power = _add_powers(_scaled_powers(chi["powers"], order - 1), chi_z["powers"])
        chi_terms = (chi, chi_z)
        piece = "n*chi^(n-1)*chi_Z*base_Z0"
    else:
        multiplier = _word_multiplier(c, chi, order)
        chi_power = _scaled_powers(chi["powers"], order)
        chi_terms = (chi,)
        piece = "chi^n*base_Z1" if z_order == 1 else "chi^n*base_Z0"
    hb_add = factor_hb + (2 * order if output_field == "F" else 2 * order - 1)
    word_order_power = 2 * order if output_field == "F" else 2 * order - 1
    word_factor = _range(c, _zero(c), S ** word_order_power / _factorial(c, word_order_power))
    powers = _add_powers(parent["powers"], chi_power)
    powers[0] += hb_add
    result = _output_term(
        parent,
        coefficient=parent["coefficient"] * factor * word_factor * multiplier,
        powers=powers,
        base_name=base_name,
        output_field=output_field,
        z_order=z_order,
        word_order=order,
        word_kind="chi_Z_insertion" if insertion else "finite_word",
        leibniz_piece=piece,
        chi_terms=chi_terms,
        formal_hb_power_added=hb_add,
    )
    result["word_factor"] = word_factor
    result["word_factor_formula"] = "[0,S^(2n)/(2n)!]" if output_field == "F" else "[0,S^(2n-1)/(2n-1)!]"
    return result


def _base_term(
    parent: Mapping[str, Any],
    *,
    factor: Any,
    factor_hb: int,
    output_field: str,
    base_name: str,
    z_order: int,
    direct: bool = False,
) -> dict[str, Any]:
    powers = list(parent["powers"])
    powers[0] += factor_hb
    return _output_term(
        parent,
        coefficient=parent["coefficient"] * factor,
        powers=powers,
        base_name=base_name,
        output_field=output_field,
        z_order=z_order,
        word_order=0,
        word_kind="direct_base" if direct else "base",
        leibniz_piece="base_Z1" if z_order == 1 else "base_Z0",
        chi_terms=(),
        formal_hb_power_added=factor_hb,
    )


def _tail_radius(c: Any, *, output_field: str, insertion: bool, m: int, M: Any, Mz: Any, S: Any, z: Any) -> tuple[Any, str]:
    A = 2 * M * S ** 2
    Az = 2 * Mz * S ** 2
    if output_field == "F" and not insertion:
        radius = A ** m / _factorial(c, 2 * m) * c.exp(z)
        formula = "hb^(2m)*A^m/(2m)!*exp(z)"
    elif output_field == "F" and insertion:
        radius = A ** (m - 1) * (m * Az) / _factorial(c, 2 * m) * (1 + z) * c.exp(z)
        formula = "hb^(2m)*m*Az*A^(m-1)/(2m)!*(1+z)*exp(z)"
    elif output_field == "Fx" and not insertion:
        radius = (2 * M) ** m * S ** (2 * m - 1) / _factorial(c, 2 * m - 1) * c.exp(z)
        formula = "hb^(2m-1)*(2M)^m*S^(2m-1)/(2m-1)!*exp(z)"
    else:
        radius = (m * 2 * Mz) * (2 * M) ** (m - 1) * S ** (2 * m - 1) / _factorial(c, 2 * m - 1) * (1 + z) * c.exp(z)
        formula = "hb^(2m-1)*m*2Mz*(2M)^(m-1)*S^(2m-1)/(2m-1)!*(1+z)*exp(z)"
    return _nonnegative(c, radius, "Volterra tail majorant"), formula


def _tail_term(
    c: Any,
    parent: Mapping[str, Any],
    *,
    base_name: str,
    output_field: str,
    z_order: int,
    insertion: bool,
    factor: Any,
    factor_hb: int,
    m: int,
    M: Any,
    Mz: Any,
    S: Any,
    z: Any,
    chi: Mapping[str, Any],
    chi_z: Mapping[str, Any],
) -> dict[str, Any]:
    majorant, formula = _tail_radius(c, output_field=output_field, insertion=insertion, m=m, M=M, Mz=Mz, S=S, z=z)
    radius = _abs_upper(c, parent["coefficient"]) * _abs_upper(c, factor) * _hi(c, majorant)
    if insertion:
        chi_power = _add_powers(_scaled_powers(chi["powers"], m - 1), chi_z["powers"])
        chi_terms = (chi, chi_z)
        piece, word_kind = "tail: base_Z0 chi_Z insertion", "chi_Z_insertion_tail"
    else:
        chi_power = _scaled_powers(chi["powers"], m)
        chi_terms = (chi,)
        piece = "tail: base_Z1 value" if z_order == 1 else "tail: base_Z0 value"
        word_kind = "value_tail"
    order_hb = 2 * m if output_field == "F" else 2 * m - 1
    powers = _add_powers(parent["powers"], chi_power)
    powers[0] += factor_hb + order_hb
    return _output_term(
        parent,
        coefficient=_symmetric(c, radius),
        powers=powers,
        base_name=base_name,
        output_field=output_field,
        z_order=z_order,
        word_order=m,
        word_kind=word_kind,
        leibniz_piece=piece,
        chi_terms=chi_terms,
        formal_hb_power_added=factor_hb + order_hb,
        tail=True,
        radius_formula=formula,
    )


def _branch(
    c: Any,
    *,
    base_name: str,
    spec: Mapping[str, Any],
    rows: Mapping[int, Sequence[Mapping[str, Any]]],
    chi: Mapping[str, Any],
    chi_z: Mapping[str, Any],
    max_order: int,
    min_order: int,
    M: Any,
    Mz: Any,
    S: Any,
    z: Any,
) -> dict[str, Any]:
    finite = {field: {"Z0": [], "Z1": []} for field in FIELDS}
    tail = {field: {"Z0": [], "Z1": []} for field in FIELDS}
    factor_F, factor_F_hb = spec["F_factor"], spec["F_hb"]

    # min_order=0 is the standalone synthetic envelope and includes its n=0
    # bases.  min_order=2 is the real-source deduplicated path: old n=0/n=1
    # records are already owned elsewhere and are not reconstructed here.
    if min_order == 0:
        for parent in rows[0]:
            finite["F"]["Z0"].append(_base_term(parent, factor=factor_F, factor_hb=factor_F_hb, output_field="F", base_name=base_name, z_order=0))
        for parent in rows[1]:
            finite["F"]["Z1"].append(_base_term(parent, factor=factor_F, factor_hb=factor_F_hb, output_field="F", base_name=base_name, z_order=1))
        direct_factor = spec["Fx_direct_factor"]
        if direct_factor is not None:
            for parent in rows[0]:
                finite["Fx"]["Z0"].append(_base_term(parent, factor=direct_factor, factor_hb=spec["Fx_direct_hb"], output_field="Fx", base_name=base_name, z_order=0, direct=True))
            for parent in rows[1]:
                finite["Fx"]["Z1"].append(_base_term(parent, factor=direct_factor, factor_hb=spec["Fx_direct_hb"], output_field="Fx", base_name=base_name, z_order=1, direct=True))

    for order in range(max(1, min_order), max_order + 1):
        for parent in rows[0]:
            finite["F"]["Z0"].append(_finite_term(c, parent, factor=factor_F, factor_hb=factor_F_hb, output_field="F", base_name=base_name, z_order=0, order=order, chi=chi, chi_z=chi_z, S=S, insertion=False))
            finite["Fx"]["Z0"].append(_finite_term(c, parent, factor=factor_F, factor_hb=factor_F_hb, output_field="Fx", base_name=base_name, z_order=0, order=order, chi=chi, chi_z=chi_z, S=S, insertion=False))
            finite["F"]["Z1"].append(_finite_term(c, parent, factor=factor_F, factor_hb=factor_F_hb, output_field="F", base_name=base_name, z_order=1, order=order, chi=chi, chi_z=chi_z, S=S, insertion=True))
            finite["Fx"]["Z1"].append(_finite_term(c, parent, factor=factor_F, factor_hb=factor_F_hb, output_field="Fx", base_name=base_name, z_order=1, order=order, chi=chi, chi_z=chi_z, S=S, insertion=True))
        for parent in rows[1]:
            finite["F"]["Z1"].append(_finite_term(c, parent, factor=factor_F, factor_hb=factor_F_hb, output_field="F", base_name=base_name, z_order=1, order=order, chi=chi, chi_z=chi_z, S=S, insertion=False))
            finite["Fx"]["Z1"].append(_finite_term(c, parent, factor=factor_F, factor_hb=factor_F_hb, output_field="Fx", base_name=base_name, z_order=1, order=order, chi=chi, chi_z=chi_z, S=S, insertion=False))

    m = max_order + 1
    for parent in rows[0]:
        tail["F"]["Z0"].append(_tail_term(c, parent, base_name=base_name, output_field="F", z_order=0, insertion=False, factor=factor_F, factor_hb=factor_F_hb, m=m, M=M, Mz=Mz, S=S, z=z, chi=chi, chi_z=chi_z))
        tail["Fx"]["Z0"].append(_tail_term(c, parent, base_name=base_name, output_field="Fx", z_order=0, insertion=False, factor=factor_F, factor_hb=factor_F_hb, m=m, M=M, Mz=Mz, S=S, z=z, chi=chi, chi_z=chi_z))
        tail["F"]["Z1"].append(_tail_term(c, parent, base_name=base_name, output_field="F", z_order=1, insertion=True, factor=factor_F, factor_hb=factor_F_hb, m=m, M=M, Mz=Mz, S=S, z=z, chi=chi, chi_z=chi_z))
        tail["Fx"]["Z1"].append(_tail_term(c, parent, base_name=base_name, output_field="Fx", z_order=1, insertion=True, factor=factor_F, factor_hb=factor_F_hb, m=m, M=M, Mz=Mz, S=S, z=z, chi=chi, chi_z=chi_z))
    for parent in rows[1]:
        tail["F"]["Z1"].append(_tail_term(c, parent, base_name=base_name, output_field="F", z_order=1, insertion=False, factor=factor_F, factor_hb=factor_F_hb, m=m, M=M, Mz=Mz, S=S, z=z, chi=chi, chi_z=chi_z))
        tail["Fx"]["Z1"].append(_tail_term(c, parent, base_name=base_name, output_field="Fx", z_order=1, insertion=False, factor=factor_F, factor_hb=factor_F_hb, m=m, M=M, Mz=Mz, S=S, z=z, chi=chi, chi_z=chi_z))

    return {
        "base": base_name,
        "input_component": spec["input"],
        "finite_prefix": finite,
        "finite": finite,
        "ordinary_Z": finite,
        "response": finite,
        "tail": tail,
        "tail_ordinary_Z": tail,
        "F": finite["F"],
        "Fx": finite["Fx"],
        "tail_F": tail["F"],
        "tail_Fx": tail["Fx"],
        "tail_starts_at_j": m,
        "ordinary_Z1_product_rule": "chi^n*base_Z1 + n*chi^(n-1)*chi_Z*base_Z0",
        "source_component_separate": True,
        "signed_finite_terms": True,
        "inherited_error_terms_separate": True,
        "formal_hb_not_materialized": True,
        "baseline_word_included": min_order == 0,
        "min_order": min_order,
        "prior_orders_not_recomputed": [0, 1] if min_order >= 2 else [0] if min_order == 1 else [],
        "new_tail_is_alternative_to_prior_j2": min_order >= 2,
        "complete_W1": False,
    }


def _merge_fields(branches: Mapping[str, Mapping[str, Any]], key: str) -> dict[str, dict[str, list[Any]]]:
    result = {field: {"Z0": [], "Z1": []} for field in FIELDS}
    for branch in branches.values():
        rows = branch[key]
        for field in FIELDS:
            for label in ORDINARY_Z_LABELS:
                result[field][label].extend(rows[field][label])
    return result


def build(
    c: Any,
    inlet: Any,
    forcing: Any,
    chi: Any,
    chi_z: Any,
    sc: Any,
    *,
    max_order: int = 0,
    min_order: int = 0,
) -> dict[str, Any]:
    """Build formal-Hb F/Fx Volterra words and independent tails.

    ``hb`` is a formal source atom, not a numerical argument.  The explicit
    condition ``0 < hb < 1/2`` is carried as a caller assumption and no
    coefficient evaluates it.  ``min_order=2`` is the real-source
    deduplicated interface: it emits only signed n=2..N words and the new
    j>=N+1 tail, as an alternative to an already-owned j=0/j=1 and old j>=2
    envelope.  The alternatives are never added together here.
    """

    _validate_context(c)
    order = _validate_order(max_order)
    start = _validate_min_order(min_order, order)
    inlet_rows = _normalise_inlet(c, inlet)
    forcing_rows = _normalise_forcing(c, forcing)
    chi_term = _normalise_axis_interval(c, chi, "chi")
    chi_z_term = _normalise_axis_interval(c, chi_z, "chi_z")
    if not _has_shared_id(chi_term["source_ids"], chi_z_term["source_ids"]):
        raise ValueError("chi and chi_z must retain a common whole-Z axis-Kchi source ID")
    sc_interval = _validate_sc(c, sc)
    term_count = sum(len(inlet_rows[field][z]) for field in ("F", "Fx") for z in ORDINARY_Z_ORDERS)
    term_count += sum(len(forcing_rows[z]) for z in ORDINARY_Z_ORDERS)
    if term_count > MAX_INPUT_TERMS:
        raise ValueError(f"input term count exceeds resource limit {MAX_INPUT_TERMS}")
    geometry = _geometry(c, sc_interval)
    specs = _base_specs(c, sc_interval, geometry)
    M = _nonnegative(c, _abs_upper(c, chi_term["coefficient"]), "M=sup|chi|")
    Mz = _nonnegative(c, _abs_upper(c, chi_z_term["coefficient"]), "Mz=sup|chi_Z|")
    S = geometry["S"]
    A = _nonnegative(c, 2 * M * S ** 2, "A")
    Az = _nonnegative(c, 2 * Mz * S ** 2, "Az")
    z = _nonnegative(c, c.sqrt(A) / 2, "z")
    branches: dict[str, dict[str, Any]] = {}
    for name in BASE_NAMES:
        input_name = specs[name]["input"]
        rows = forcing_rows if name == "forcing" else inlet_rows[input_name]
        branches[name] = _branch(c, base_name=name, spec=specs[name], rows=rows, chi=chi_term, chi_z=chi_z_term, max_order=order, min_order=start, M=M, Mz=Mz, S=S, z=z)
    finite = _merge_fields(branches, "finite_prefix")
    tail = _merge_fields(branches, "tail")
    output_count = sum(len(finite[field][z]) + len(tail[field][z]) for field in FIELDS for z in ORDINARY_Z_LABELS)
    if output_count > MAX_OUTPUT_TERMS:
        raise ValueError(f"output term count exceeds resource limit {MAX_OUTPUT_TERMS}")
    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "conditional_formal_Hb_axis_Kchi_volterra_finite_tail",
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "ordinary_Z_labels": list(ORDINARY_Z_LABELS),
        "fields": list(FIELDS),
        "base_names": list(BASE_NAMES),
        "branch_names": list(BASE_NAMES),
        "inlet": inlet_rows,
        "forcing": forcing_rows,
        "chi": chi_term,
        "chi_z": chi_z_term,
        "sc": sc_interval,
        "branches": branches,
        "finite_prefix": finite,
        "finite": finite,
        "response": finite,
        "ordinary_Z": finite,
        "tail": tail,
        "tail_ordinary_Z": tail,
        "F": finite["F"],
        "Fx": finite["Fx"],
        "tail_F": tail["F"],
        "tail_Fx": tail["Fx"],
        "ordinary_Z1_product_rule": "chi^n*base_Z1 + n*chi^(n-1)*chi_Z*base_Z0",
        "max_order": order,
        "min_order": start,
        "finite_orders_emitted": ([0] if start == 0 else []) + list(range(max(1, start), order + 1)),
        "tail_starts_at_j": order + 1,
        "geometry": geometry,
        "tail_majorant": {
            "M": M, "Mz": Mz, "A": A, "Az": Az, "z": z, "m": order + 1,
            "zero_power_convention": "0^0=1 when M=0 and N=0",
            "F_value": "hb^(2m)*A^m/(2m)!*exp(z)",
            "F_Z_insertion": "hb^(2m)*m*Az*A^(m-1)/(2m)!*(1+z)*exp(z)",
            "Fx_value": "hb^(2m-1)*(2M)^m*S^(2m-1)/(2m-1)!*exp(z)",
            "Fx_Z_insertion": "hb^(2m-1)*m*2Mz*(2M)^(m-1)*S^(2m-1)/(2m-1)!*(1+z)*exp(z)",
        },
        "conditions": {
            "hb_strictly_between_zero_and_half": True,
            "hb_condition_statement": "0 < hb < 1/2 (caller-certified formal condition)",
            "hb_condition_is_caller_assumption": True,
            "hb_materialized": False,
            "sc_positive": True,
            "hb_Z_independent": True,
            "sc_Z_independent": True,
            "geometry_Z_independent": True,
            "chi_scope": "same-source axis Kchi whole-Z interval",
            "chi_source_common": True,
            "inlet_rows_complete": True,
            "forcing_rows_complete": True,
            "forcing_rows_supplied": bool(forcing_rows[0] or forcing_rows[1]),
            "forcing_empty_means_no_local_forcing": not bool(forcing_rows[0] or forcing_rows[1]),
            "fake_zero_forcing_filled": False,
            "ordinary_Z_rows_independent": True,
        },
        "alternative_envelopes": {
            "new_finite_and_tail": True,
            "prior_j0_j1_recomputed": start == 0,
            "prior_j2_tail_supplied": False,
            "new_tail_is_alternative_to_prior_j2": start >= 2,
            "combined": False,
            "never_add_old_and_new_tail": True,
        },
        "scope": {
            "whole_collar_conservative_interval": True,
            "complete_W1": False,
            "full_six_field": False,
            "time_recursion": False,
            "m15": False,
            "ledger_written": False,
            "entered_212_ledger": False,
            "same_owner_proof": False,
            "numeric_function_evaluation": False,
            "midpoint_values_used": False,
        },
        "resource_limits": {
            "max_order": MAX_ORDER,
            "max_input_terms": MAX_INPUT_TERMS,
            "max_output_terms": MAX_OUTPUT_TERMS,
            "input_term_count": term_count,
            "output_term_count": output_count,
        },
        "source_ids_preserved": True,
        "parent_provenance_copied": True,
        "input_error_only_preserved": True,
        "inherited_error_terms_separate": True,
        "tail_terms_independent": True,
        "signed_finite_terms": True,
        "formal_hb_powers_not_materialized": True,
        "no_ancestor_helper_called": True,
        "receipt_reads": False,
        "saved_data_read": False,
    }


build_formal_axis_kchi_volterra = build
build_tail = build
enclose = build
evaluate = build


