"""Safe whole-collar ``j >= 3`` Kchi tails for a saved ``G0g`` inlet.

This is a deliberately independent consumer of already grouped ``F`` and
``Fx`` ordinary-Z inlet intervals.  It does not import or call an earlier
formal-axis or inlet-tail leaf, does not reconstruct a finite word, and does
not consume a receipt.  The only returned sectors are symmetric error-only
majorants for the omitted words starting at ``j = 3``.

The chart assumptions are the caller's source contract

``L <= hb*S``, ``S = (sc/4)*exp(sc/16)``, ``0 <= K0(x,t) <= x-t`` and
``0 <= K0_x(x,t) <= 1``.

Put ``M = sup|chi|``, ``Mz = sup|chi_Z|`` and ``A = 2*M*S**2``.  For
``m = 3`` the positive F series is bounded by

``(2*M)**m*S**(2*m)/(2*m)! * exp(A/4)``

and the Fx series by the same expression with ``S**(2*m-1)/(2*m-1)!``.
The ordinary-Z insertion has the extra factor
``(1 + A/(4*m)) = (1 + A/12)``.  Indeed, for ``r = 5`` or ``6`` and
``k >= 0``, ``(r+2*k)!/r! >= 4**k*k!``.  After writing ``j=m+k``, the
insertion weight is ``(m+k)/m``; summing the resulting exponential series
gives ``(1+A/(4*m))*exp(A/4)``.  This is an all-A factorial majorant, so no
smallness condition on A is silently used.

The input F and Fx routes have different formal width powers.  The four
routes are F<-F, F<-Fx, Fx<-F and Fx<-Fx, with formal hb shifts
``6, 7, 5, 6`` respectively.  The hb atom is never materialized: it is
recorded only in the output power vector and formula metadata.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "inlet_kchi_tail_safe_v1"
COMPONENT = "formal_Hb_axis_Kchi_inlet_j_ge_3_tail"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
ORDINARY_Z_LABELS = ("Z0", "Z1")
FIELDS = ("F", "Fx")
INLET_COMPONENTS = ("F", "Fx")
BRANCHES = ("F_from_F", "F_from_Fx", "Fx_from_F", "Fx_from_Fx")
TAIL_START = 3
FORMAL_HB_POWER_ADDED = {
    "F_from_F": 6,
    "F_from_Fx": 7,
    "Fx_from_F": 5,
    "Fx_from_Fx": 6,
}
FIELD_FORMAL_HB_POWER_ADDED = {"F": 6, "Fx": 5}
TAIL_ERROR_CLASS = "formal_kchi_j_ge_3_tail"
MAX_INPUT_TERMS = 1024
MAX_OUTPUT_TERMS = 8192

_COEFFICIENT_KEYS = (
    "coefficient",
    "coefficient_interval",
    "interval",
    "enclosure",
    "box",
    "value",
)
_POWER_KEYS = ("powers", "formal_powers", "source_powers")
_MISSING = object()

__all__ = [
    "SCHEMA",
    "COMPONENT",
    "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "ORDINARY_Z_LABELS",
    "FIELDS",
    "INLET_COMPONENTS",
    "BRANCHES",
    "TAIL_START",
    "FORMAL_HB_POWER_ADDED",
    "FIELD_FORMAL_HB_POWER_ADDED",
    "TAIL_ERROR_CLASS",
    "build",
    "build_inlet_tail",
    "build_tail",
    "inlet_tail",
    "evaluate",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _items(value: Any) -> list[Any]:
    if value is None or value is _MISSING:
        return []
    return list(value) if _is_sequence(value) else [value]


def _clone(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if _is_sequence(value):
        return [_clone(item) for item in value]
    return value


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    seen: set[tuple[str, str]] = set()
    for value in values:
        key = (type(value).__name__, repr(value))
        if key not in seen:
            seen.add(key)
            result.append(value)
    return result


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
        if nested is not None:
            return nested[1 if upper else 0]
        nested = getattr(endpoint, "_mpf_", None)
        if nested is not None:
            return nested
        return endpoint
    return value


def _make(c: Any, lower: Any, upper: Any) -> Any:
    return c.make_mpf((_raw(lower, False), _raw(upper, True)))


def _exact_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    if not isinstance(value, (tuple, list)) or len(value) != 4:
        raise TypeError(f"{label} must be a four-integer MPF tuple")
    if any(type(item) is not int for item in value):
        raise TypeError(f"{label} may not contain decimal strings or floats")
    sign, mantissa, exponent, bitcount = tuple(int(item) for item in value)
    if sign not in (0, 1) or mantissa < 0 or bitcount < 0:
        raise ValueError(f"{label} has invalid MPF sign/mantissa/bitcount")
    if mantissa == 0 and bitcount != 0:
        raise ValueError(f"{label} is a non-finite MPF special value")
    return sign, mantissa, exponent, bitcount


def _iv(c: Any, value: Any, label: str = "interval") -> Any:
    """Decode a directed interval without a midpoint or float conversion."""

    if hasattr(value, "_mpi_"):
        result = _make(c, value, value)
    elif isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            if "lower_exact_mpf_tuple" not in value or "upper_exact_mpf_tuple" not in value:
                raise TypeError(f"{label} requires both exact MPF endpoints")
            lower = _exact_tuple(value["lower_exact_mpf_tuple"], label + ".lower_exact_mpf_tuple")
            upper = _exact_tuple(value["upper_exact_mpf_tuple"], label + ".upper_exact_mpf_tuple")
            result = c.make_mpf((lower, upper))
        else:
            nested = next((value[key] for key in _COEFFICIENT_KEYS if key in value), _MISSING)
            if nested is not _MISSING and nested is not value:
                result = _iv(c, nested, label)
            else:
                lower = next((value[key] for key in ("lower", "lo", "left") if key in value), _MISSING)
                upper = next((value[key] for key in ("upper", "hi", "right") if key in value), _MISSING)
                if lower is _MISSING or upper is _MISSING:
                    raise TypeError(f"{label} is not interval-compatible")
                if isinstance(lower, (tuple, list)) and isinstance(upper, (tuple, list)):
                    lower = _exact_tuple(lower, label + ".lower")
                    upper = _exact_tuple(upper, label + ".upper")
                result = _make(c, lower, upper)
    elif _is_sequence(value) and len(value) == 2:
        result = _make(c, value[0], value[1])
    else:
        try:
            point = c.mpf(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{label} is not interval-compatible") from exc
        result = _make(c, point, point)
    raw = getattr(result, "_mpi_", None)
    if raw is None:
        raise TypeError(f"{label} did not produce a directed interval")
    lower = c.make_mpf((raw[0], raw[0]))
    upper = c.make_mpf((raw[1], raw[1]))
    if upper < lower:
        raise ValueError(f"{label} has reversed directed endpoints")
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


def _range(c: Any, lower: Any, upper: Any, label: str) -> Any:
    lower, upper = _lo(c, lower), _hi(c, upper)
    if upper < lower:
        raise ValueError(f"{label} has reversed endpoints")
    return _make(c, lower, upper)


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _iv(c, value, label)
    lower, upper = _lo(c, value), _hi(c, value)
    if upper < 0:
        raise ValueError(f"{label} has a negative upper endpoint")
    if lower < 0:
        lower = _zero(c)
    return _range(c, lower, upper, label)


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
    return _make(c, -upper, upper)


def _factorial(c: Any, n: int) -> Any:
    result = _one(c)
    for index in range(2, n + 1):
        result *= index
    return result


def _pow_int(c: Any, value: Any, exponent: int) -> Any:
    if exponent < 0:
        raise ValueError("integer interval powers must be nonnegative")
    result = _one(c)
    for _ in range(exponent):
        result *= value
    return result


def _is_bool(value: Any) -> bool:
    return isinstance(value, bool)


def _normalise_powers(raw: Mapping[Any, Any], label: str) -> list[int]:
    powers = next((raw[key] for key in _POWER_KEYS if key in raw), _MISSING)
    if powers is _MISSING:
        powers = (
            raw.get("hb_power", 0),
            raw.get("F0_power", raw.get("F0_anchor_power", 0)),
            raw.get("Pstar_squared_power", raw.get("Pstar_power", 0)),
            raw.get("Lambda_power", raw.get("lambda_power", 0)),
        )
    if isinstance(powers, Mapping):
        powers = [powers.get(atom, _MISSING) for atom in SOURCE_POWER_BASIS]
    if not _is_sequence(powers) or len(powers) != 4:
        raise TypeError(f"{label}.powers must contain four integer formal powers")
    if any(_is_bool(item) or not isinstance(item, int) for item in powers):
        raise TypeError(f"{label}.powers must contain four integer formal powers")
    return [int(item) for item in powers]


def _meta_values(raw: Mapping[Any, Any], names: Sequence[str]) -> list[Any]:
    values: list[Any] = []
    for name in names:
        if name in raw:
            values.extend(_items(raw[name]))
    for container_name in ("origin", "metadata", "original_term"):
        container = raw.get(container_name)
        if isinstance(container, Mapping):
            for name in names:
                if name in container:
                    values.extend(_items(container[name]))
    return _unique(values)


def _path_values(raw: Mapping[Any, Any]) -> list[Any]:
    values: list[Any] = []
    if "source_path" in raw:
        values.append(_clone(raw["source_path"]))
    if "source_paths" in raw:
        paths = raw["source_paths"]
        if _is_sequence(paths) and any(_is_sequence(item) or isinstance(item, Mapping) for item in paths):
            values.extend(_clone(item) for item in paths)
        else:
            values.append(_clone(paths))
    for container_name in ("origin", "metadata", "original_term"):
        container = raw.get(container_name)
        if isinstance(container, Mapping):
            if "source_path" in container:
                values.append(_clone(container["source_path"]))
            if "source_paths" in container:
                values.extend(_items(container["source_paths"]))
    return _unique(values)


def _normalise_term(c: Any, raw: Any, label: str) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must be a term mapping")
    coefficient = next((raw[key] for key in _COEFFICIENT_KEYS if key in raw), _MISSING)
    if coefficient is _MISSING:
        raise KeyError(f"{label} is missing coefficient")
    coefficient_interval = _iv(c, coefficient, label + ".coefficient")
    powers = _normalise_powers(raw, label)
    source_ids = _meta_values(raw, ("source_ids", "source_id", "sources", "source", "input_id"))
    if not source_ids:
        raise ValueError(f"{label} must retain at least one source_id")
    source_paths = _path_values(raw)
    error_classes = _meta_values(raw, ("error_classes", "error_class", "uncertainty_class"))
    source_error_classes = _meta_values(raw, ("source_error_classes", "source_error_class"))
    uncertainty = raw.get("uncertainty", raw.get("is_uncertainty", False))
    error_only = raw.get("error_only", False)
    symmetric_error = raw.get("symmetric_error", error_only)
    signed_field = raw.get("signed_field", False)
    for key, value in (
        ("uncertainty", uncertainty),
        ("error_only", error_only),
        ("symmetric_error", symmetric_error),
        ("signed_field", signed_field),
    ):
        if not isinstance(value, bool):
            raise TypeError(f"{label}.{key} must be boolean")
    if error_only and not symmetric_error:
        raise ValueError(f"{label} error_only requires symmetric_error")
    if error_only and signed_field:
        raise ValueError(f"{label} error_only may not be signed_field")
    if error_only and _lo(c, coefficient_interval) != -_hi(c, coefficient_interval):
        raise ValueError(f"{label} error_only coefficient must be symmetric")
    if not source_error_classes:
        source_error_classes = list(error_classes)
    result = dict(raw)
    result.update(
        {
            "coefficient": coefficient_interval,
            "coefficient_interval": coefficient_interval,
            "powers": list(powers),
            "formal_powers": list(powers),
            "power_basis": list(SOURCE_POWER_BASIS),
            "source_ids": list(source_ids),
            "source_paths": list(source_paths),
            "source_path": _clone(raw.get("source_path")),
            "source_error_classes": list(source_error_classes),
            "error_classes": list(error_classes),
            "uncertainty": uncertainty,
            "error_only": error_only,
            "symmetric_error": symmetric_error,
            "signed_field": signed_field,
            "original_term": _clone(raw.get("original_term")),
        }
    )
    return result


def _order_key(key: Any) -> int | None:
    if key == 0 or (isinstance(key, str) and key.lower() in {"0", "z0", "ordinary_z0", "value"}):
        return 0
    if key == 1 or (isinstance(key, str) and key.lower() in {"1", "z1", "ordinary_z1", "derivative"}):
        return 1
    return None


def _term_list(c: Any, raw: Any, label: str) -> list[dict[str, Any]]:
    if isinstance(raw, Mapping) and any(key in raw for key in _COEFFICIENT_KEYS):
        values = [raw]
    elif _is_sequence(raw):
        values = list(raw)
    else:
        raise TypeError(f"{label} must be a term sequence")
    if not values:
        raise ValueError(f"{label} must not be empty")
    return [_normalise_term(c, item, f"{label}[{index}]") for index, item in enumerate(values)]


def _rows(c: Any, raw: Any, label: str) -> dict[int, list[dict[str, Any]]]:
    current = raw
    if isinstance(current, Mapping):
        for wrapper in ("ordinary_Z", "source_Z", "rows", "jets", "terms"):
            nested = current.get(wrapper, _MISSING)
            if nested is not _MISSING and nested is not current:
                current = nested
                break
    if not isinstance(current, Mapping):
        raise TypeError(f"{label} requires explicit Z0 and Z1 rows")
    found: dict[int, Any] = {}
    for key, value in current.items():
        index = _order_key(key)
        if index is not None:
            if index in found:
                raise ValueError(f"{label} supplies Z{index} twice")
            found[index] = value
    if set(found) != {0, 1}:
        raise ValueError(f"{label} requires explicit ordinary Z0 and Z1 rows")
    return {index: _term_list(c, found[index], f"{label}.Z{index}") for index in ORDINARY_Z_ORDERS}


def _lookup(mapping: Mapping[Any, Any], names: Sequence[str]) -> Any:
    for name in names:
        if name in mapping:
            return mapping[name]
    lowered = {str(key).lower(): value for key, value in mapping.items()}
    for name in names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    return _MISSING


def _normalise_inlet(c: Any, inlet: Any) -> dict[str, dict[int, list[dict[str, Any]]]]:
    if not isinstance(inlet, Mapping):
        raise TypeError("inlet must contain F and Fx ordinary rows")
    current = inlet
    for wrapper in ("inlet", "G0g", "g0g", "aggregate", "output", "fields"):
        nested = current.get(wrapper, _MISSING)
        if isinstance(nested, Mapping) and nested is not current:
            current = nested
            break
    f_raw = _lookup(current, ("F", "Fc", "F1"))
    fx_raw = _lookup(current, ("Fx", "Fxc", "d_x_F", "d_x_F1"))
    # Saved aggregate packets sometimes put the field names inside each Z row.
    if f_raw is _MISSING or fx_raw is _MISSING:
        row_values: dict[int, Mapping[Any, Any]] = {}
        for key, value in current.items():
            order = _order_key(key)
            if order is not None and isinstance(value, Mapping):
                row_values[order] = value
        if set(row_values) == {0, 1}:
            f_raw = {f"Z{order}": _lookup(row_values[order], ("F", "Fc", "F1")) for order in ORDINARY_Z_ORDERS}
            fx_raw = {f"Z{order}": _lookup(row_values[order], ("Fx", "Fxc", "d_x_F", "d_x_F1")) for order in ORDINARY_Z_ORDERS}
    if f_raw is _MISSING or fx_raw is _MISSING or f_raw is None or fx_raw is None:
        raise ValueError("inlet must contain both F and Fx rows")
    if f_raw is _MISSING or fx_raw is _MISSING:
        raise ValueError("inlet must contain both F and Fx rows")
    return {"F": _rows(c, f_raw, "inlet.F"), "Fx": _rows(c, fx_raw, "inlet.Fx")}


def _axis_root(value: Any) -> str | None:
    text = str(value).replace("_", "").replace("-", "").lower()
    marker = "axis:wholez"
    if marker not in text:
        return None
    return text.split(marker, 1)[0] + marker


def _normalise_axis(c: Any, raw: Any, label: str) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must be a whole-Z axis interval mapping")
    current = raw
    for wrapper in ("chi", "chi_z", "chi_Z", "axis", "value"):
        nested = current.get(wrapper, _MISSING)
        if isinstance(nested, Mapping) and nested is not current and not any(key in current for key in _COEFFICIENT_KEYS):
            current = nested
            break
    item = dict(current)
    item.setdefault("powers", [0, 0, 0, 0])
    term = _normalise_term(c, item, label)
    if term["powers"] != [0, 0, 0, 0]:
        raise ValueError(f"{label} formal powers must be zero")
    if term["error_only"]:
        raise ValueError(f"{label} must be an axis value enclosure, not error-only")
    if raw.get("whole_Z", raw.get("whole_z", True)) is False or current.get("whole_Z", current.get("whole_z", True)) is False:
        raise ValueError(f"{label} must cover whole Z")
    scope = " ".join(str(current.get(key, raw.get(key, ""))) for key in ("scope", "domain", "role", "support")).lower()
    if any(token in scope for token in ("lower-collar", "lower collar", "shear", "lower_band", "lower band")):
        raise ValueError(f"{label} may not use lower-collar shear chi")
    term["axis_whole_Z_certified"] = True
    term["scope"] = current.get("scope", "whole-Z axis Kchi")
    return term


def _axis_pair(chi: Mapping[str, Any], chi_z: Mapping[str, Any]) -> list[Any]:
    left = set(chi.get("source_ids", []))
    right = set(chi_z.get("source_ids", []))
    common = [value for value in chi.get("source_ids", []) if value in right]
    if common:
        return _unique(common)
    left_roots = {_axis_root(value) for value in left} - {None}
    right_roots = {_axis_root(value) for value in right} - {None}
    if left_roots.intersection(right_roots):
        return sorted(left_roots.intersection(right_roots))
    raise ValueError("chi and chi_z must retain a common whole-Z source operator")


def _validate_context(c: Any) -> None:
    if not all(hasattr(c, name) for name in ("mpf", "make_mpf", "exp")):
        raise TypeError("c must provide directed mpf, make_mpf and exp arithmetic")


def _validate_sc(c: Any, sc: Any) -> Any:
    if sc is None:
        raise TypeError("positive sc is required")
    if isinstance(sc, Mapping):
        for key in ("positive", "certified_positive", "Z_independent", "z_independent", "sc_Z_independent"):
            if key in sc and not isinstance(sc[key], bool):
                raise TypeError(f"sc.{key} must be boolean")
            if key in sc and sc[key] is False:
                raise ValueError(f"sc.{key} certificate is false")
    return _positive(c, _iv(c, sc, "sc"), "sc")


def _geometry(c: Any, sc: Any) -> dict[str, Any]:
    S = _nonnegative(c, (sc / 4) * c.exp(sc / 16), "S")
    return {
        "S": S,
        "S_formula": "S=(sc/4)*exp(sc/16)",
        "x_formula": "x(q)=2*exp(hb*sc*q/8), 0<=q<=1",
        "L_bound": "L <= hb*S",
        "K0_bound": "0 <= K0(x,t) <= x-t",
        "K0x_bound": "0 <= K0_x(x,t) <= 1",
        "hb_condition": "0 < hb < 1/2 (caller-certified formal condition)",
        "formal_hb_not_materialized": True,
        "geometry_Z_independent": True,
        "midpoint_values_used": False,
    }


def _tail_constants(c: Any, chi: Mapping[str, Any], chi_z: Mapping[str, Any], S: Any) -> dict[str, Any]:
    M = _range(c, _zero(c), _abs_upper(c, chi["coefficient"]), "M")
    Mz = _range(c, _zero(c), _abs_upper(c, chi_z["coefficient"]), "Mz")
    A = _nonnegative(c, 2 * M * _pow_int(c, S, 2), "A=2*M*S^2")
    Az = _nonnegative(c, 2 * Mz * _pow_int(c, S, 2), "Az=2*Mz*S^2")
    m = TAIL_START
    exp_A_over_4 = _nonnegative(c, c.exp(A / 4), "exp(A/4)")
    two_M = 2 * M
    two_Mz = 2 * Mz
    value_F = _nonnegative(c, _pow_int(c, two_M, m) * _pow_int(c, S, 2 * m) / _factorial(c, 2 * m) * exp_A_over_4, "F value tail radius")
    value_Fx = _nonnegative(c, _pow_int(c, two_M, m) * _pow_int(c, S, 2 * m - 1) / _factorial(c, 2 * m - 1) * exp_A_over_4, "Fx value tail radius")
    derivative_factor = (1 + A / (4 * m)) * exp_A_over_4
    # Keep the algebra in a directed-interval form while exposing the exact
    # simplified factor in the public metadata.
    insertion_F = _nonnegative(c, m * two_Mz * _pow_int(c, two_M, m - 1) * _pow_int(c, S, 2 * m) / _factorial(c, 2 * m) * derivative_factor, "F Z1 insertion tail radius")
    insertion_Fx = _nonnegative(c, m * two_Mz * _pow_int(c, two_M, m - 1) * _pow_int(c, S, 2 * m - 1) / _factorial(c, 2 * m - 1) * derivative_factor, "Fx Z1 insertion tail radius")
    return {
        "M": M,
        "Mz": Mz,
        "A": A,
        "Az": Az,
        "m": m,
        "exp_A_over_4": exp_A_over_4,
        "value_F": value_F,
        "value_Fx": value_Fx,
        "insertion_F": insertion_F,
        "insertion_Fx": insertion_Fx,
        "F_value_formula": "(2*M)^3*S^6/6!*exp(A/4)",
        "Fx_value_formula": "(2*M)^3*S^5/5!*exp(A/4)",
        "F_Z1_insertion_formula": "3*(2*Mz)*(2*M)^2*S^6/6!*(1+A/12)*exp(A/4)",
        "Fx_Z1_insertion_formula": "3*(2*Mz)*(2*M)^2*S^5/5!*(1+A/12)*exp(A/4)",
        "factorial_offsets": {"F": 6, "Fx": 5},
        "factorial_comparison": "(r+2*k)!/r! >= 4^k*k! for r=5,6",
        "insertion_weight": "(3+k)/3",
        "bound_proof": "write j=3+k; factorial comparison gives exp(A/4), and (3+k)/3 gives (1+A/12)*exp(A/4)",
        "large_A_safe": True,
        "zero_power_convention": "0^0=1; here m=3 so a zero chi majorant vanishes",
    }


def _merge(values: Sequence[Any]) -> list[Any]:
    return _unique([item for value in values for item in value])


def _branch_factor(c: Any, sc: Any, branch: str) -> tuple[Any, str]:
    if branch in ("F_from_F", "Fx_from_F"):
        return _one(c), "1"
    return _range(c, _zero(c), sc / 4, "Fx-to-field source factor"), "[0,sc/4]"


def _output_term(
    c: Any,
    parent: Mapping[str, Any],
    *,
    chi_terms: Sequence[Mapping[str, Any]],
    branch: str,
    group: str,
    field: str,
    input_component: str,
    source_order: int,
    chi_order: int,
    insertion: bool,
    majorant: Any,
    base_factor: Any,
    base_factor_formula: str,
    source_index: int,
) -> dict[str, Any]:
    hb_added = FORMAL_HB_POWER_ADDED[branch]
    radius_factor = _nonnegative(c, base_factor * majorant, branch + " radius factor")
    radius = _nonnegative(c, _abs_upper(c, parent["coefficient"]) * _hi(c, radius_factor), "tail radius")
    chi_ids = _merge([term.get("source_ids", []) for term in chi_terms])
    chi_paths = _merge([term.get("source_paths", []) for term in chi_terms])
    chi_classes = _merge([term.get("error_classes", []) for term in chi_terms])
    source_ids = _unique(list(parent.get("source_ids", [])) + chi_ids)
    # Keep the inlet's source path in the primary path field, matching the
    # saved grouped-sector record.  Axis paths remain available separately in
    # ``chi_source_paths``; they are not folded into the inlet path.
    source_paths = list(parent.get("source_paths", []))
    source_error_classes = _unique(list(parent.get("source_error_classes", [])) + chi_classes)
    error_classes = _unique(list(parent.get("error_classes", [])) + chi_classes + [TAIL_ERROR_CLASS, "symmetric_error"])
    powers = list(parent["powers"])
    powers[0] += hb_added
    if insertion:
        route = "3*chi^2*chi_Z*inlet_Z0"
        word_kind = "j_ge_3_chi_Z_insertion_tail"
        leibniz_piece = "j*chi^(j-1)*chi_Z*inlet_Z0"
    else:
        route = f"chi^3*inlet_Z{source_order}"
        word_kind = "j_ge_3_value_tail"
        leibniz_piece = f"chi^j*inlet_Z{source_order}"
    origin = {
        "operation": "inlet_kchi_j_ge_3_safe_tail",
        "group": group,
        "branch": branch,
        "field": field,
        "input_component": input_component,
        "source_order": source_order,
        "chi_order": chi_order,
        "leibniz_piece": leibniz_piece,
        "source_ids": source_ids,
        "source_paths": source_paths,
        "error_classes": error_classes,
        "error_only": True,
        "symmetric_error": True,
    }
    coefficient = _symmetric(c, radius)
    return {
        "coefficient": coefficient,
        "coefficient_interval": coefficient,
        "coefficient_is_symmetric_radius": True,
        "coefficient_is_function_value": False,
        "known_error_upper_is_function_value": False,
        "coefficient_contains_formal_hb": False,
        "signed_tail_value": False,
        "signed_source_coefficient_preserved": True,
        "parent_coefficient": parent["coefficient"],
        "signed_parent_coefficient": parent["coefficient"],
        "powers": powers,
        "formal_powers": list(powers),
        "formal_power_shift": [hb_added, 0, 0, 0],
        "power_basis": list(SOURCE_POWER_BASIS),
        "source_ids": source_ids,
        "parent_source_ids": list(parent.get("source_ids", [])),
        "inlet_source_ids": list(parent.get("source_ids", [])),
        "chi_source_ids": chi_ids,
        "chi_Z_source_ids": _merge([chi_terms[1].get("source_ids", [])]) if insertion else [],
        "source_paths": source_paths,
        "source_path": _clone(parent.get("source_path")),
        "parent_source_paths": list(parent.get("source_paths", [])),
        "inlet_source_paths": list(parent.get("source_paths", [])),
        "chi_source_paths": chi_paths,
        "source_error_classes": source_error_classes,
        "parent_source_error_classes": list(parent.get("source_error_classes", [])),
        "parent_error_classes": list(parent.get("error_classes", [])),
        "error_classes": error_classes,
        "uncertainty": True,
        "parent_uncertainty": bool(parent.get("uncertainty", False)),
        "error_only": True,
        "inherited_error_only": bool(parent.get("error_only", False)),
        "parent_error_only": bool(parent.get("error_only", False)),
        "symmetric_error": True,
        "parent_symmetric_error": bool(parent.get("symmetric_error", False)),
        "signed_field": False,
        "parent_signed_field": parent.get("signed_field", False),
        "derived_role": "symmetric_error_only",
        "tail": True,
        "tail_starts_at_j": TAIL_START,
        "word_order": TAIL_START,
        "word_kind": word_kind,
        "term_kind": "volterra_tail",
        "field": field,
        "ordinary_Z": f"Z{1 if insertion or source_order == 1 else 0}",
        "ordinary_Z_order": 1 if insertion or source_order == 1 else 0,
        "source_order": source_order,
        "chi_order": chi_order,
        "leibniz_piece": leibniz_piece,
        "route": route,
        "branch": branch,
        "group": group,
        "input_component": input_component,
        "incoming_component": input_component,
        "source_term_index": source_index,
        "formal_hb_power_added": hb_added,
        "formal_hb_not_materialized": True,
        "radius": coefficient,
        "radius_factor": radius_factor,
        "base_factor": base_factor,
        "base_factor_formula": base_factor_formula,
        "tail_majorant_role": "j>=3 positive factorial majorant",
        "tail_radius_is_symmetric": True,
        "error_only_tail": True,
        "signed_expression_interval": False,
        "inherited_error_preserved": True,
        "source_term_metadata_preserved": True,
        "original_term": _clone(parent.get("original_term")),
        "parent_origin": _clone(parent.get("origin")),
        "origin": origin,
        "midpoint_values_used": False,
        "old_j2_tail_combined": False,
        "old_j2_error_subtracted": False,
        "prior_word_recomputed": False,
    }


def _branch_set(
    c: Any,
    inlet: Mapping[str, Mapping[int, Sequence[Mapping[str, Any]]]],
    chi: Mapping[str, Any],
    chi_z: Mapping[str, Any],
    constants: Mapping[str, Any],
    sc: Any,
    branch: str,
    group: str,
) -> dict[str, Any]:
    field = "F" if branch.startswith("F_") else "Fx"
    input_component = "Fx" if branch.endswith("_Fx") else "F"
    rows = {"Z0": [], "Z1": []}
    base_factor, base_factor_formula = _branch_factor(c, sc, branch)
    value_majorant = constants["value_F"] if field == "F" else constants["value_Fx"]
    insertion_majorant = constants["insertion_F"] if field == "F" else constants["insertion_Fx"]
    for index, parent in enumerate(inlet[input_component][0]):
        rows["Z0"].append(_output_term(
            c, parent, chi_terms=(chi,), branch=branch, group=group, field=field,
            input_component=input_component, source_order=0, chi_order=0,
            insertion=False, majorant=value_majorant, base_factor=base_factor,
            base_factor_formula=base_factor_formula, source_index=index,
        ))
    for index, parent in enumerate(inlet[input_component][1]):
        rows["Z1"].append(_output_term(
            c, parent, chi_terms=(chi,), branch=branch, group=group, field=field,
            input_component=input_component, source_order=1, chi_order=0,
            insertion=False, majorant=value_majorant, base_factor=base_factor,
            base_factor_formula=base_factor_formula, source_index=index,
        ))
    for index, parent in enumerate(inlet[input_component][0]):
        rows["Z1"].append(_output_term(
            c, parent, chi_terms=(chi, chi_z), branch=branch, group=group, field=field,
            input_component=input_component, source_order=0, chi_order=1,
            insertion=True, majorant=insertion_majorant, base_factor=base_factor,
            base_factor_formula=base_factor_formula, source_index=index,
        ))
    return {
        "branch": branch,
        "field": field,
        "input_component": input_component,
        "response": rows,
        "tail": rows,
        "remainder": rows,
        "Z0": rows["Z0"],
        "Z1": rows["Z1"],
        "tail_starts_at_j": TAIL_START,
        "kernel_tail_order": "j>=3",
        "formal_hb_power_added": FORMAL_HB_POWER_ADDED[branch],
        "formal_hb_not_materialized": True,
        "radius_factors": {
            "Z0": base_factor * value_majorant,
            "Z1_source": base_factor * value_majorant,
            "Z1_chi": base_factor * insertion_majorant,
        },
        "signed_sectors": False,
        "all_terms_error_only": True,
        "error_only_sectors_separate": True,
        "source_terms_individual": True,
        "old_j2_tail_combined": False,
    }


def _aggregate(branches: Mapping[str, Mapping[str, Any]]) -> dict[str, dict[str, list[Any]]]:
    return {
        "F": {
            "Z0": list(branches["F_from_F"]["Z0"]) + list(branches["F_from_Fx"]["Z0"]),
            "Z1": list(branches["F_from_F"]["Z1"]) + list(branches["F_from_Fx"]["Z1"]),
        },
        "Fx": {
            "Z0": list(branches["Fx_from_F"]["Z0"]) + list(branches["Fx_from_Fx"]["Z0"]),
            "Z1": list(branches["Fx_from_F"]["Z1"]) + list(branches["Fx_from_Fx"]["Z1"]),
        },
    }


def _empty_packet() -> dict[str, dict[str, list[Any]]]:
    return {field: {label: [] for label in ORDINARY_Z_LABELS} for field in FIELDS}


def build(c: Any, inlet: Any, chi: Any, chi_z: Any, sc: Any) -> dict[str, Any]:
    """Build only the safe independent inlet ``j >= 3`` tail.

    The positional interface is intentional: ``inlet`` is the saved grouped
    F/Fx packet, ``chi`` and ``chi_z`` are whole-Z axis intervals, and ``sc``
    is the positive geometry scale.  No numerical ``hb`` argument is accepted.
    """

    _validate_context(c)
    inlet_rows = _normalise_inlet(c, inlet)
    chi_term = _normalise_axis(c, chi, "chi")
    chi_z_term = _normalise_axis(c, chi_z, "chi_z")
    common_axis_ids = _axis_pair(chi_term, chi_z_term)
    sc_interval = _validate_sc(c, sc)
    input_count = sum(len(inlet_rows[field][order]) for field in INLET_COMPONENTS for order in ORDINARY_Z_ORDERS)
    if input_count > MAX_INPUT_TERMS:
        raise ValueError(f"input term count exceeds cap {MAX_INPUT_TERMS}")
    geometry = _geometry(c, sc_interval)
    constants = _tail_constants(c, chi_term, chi_z_term, geometry["S"])
    branches = {
        branch: _branch_set(c, inlet_rows, chi_term, chi_z_term, constants, sc_interval, branch, "G0g")
        for branch in BRANCHES
    }
    tail = _aggregate(branches)
    output_count = sum(len(tail[field][label]) for field in FIELDS for label in ORDINARY_Z_LABELS)
    if output_count > MAX_OUTPUT_TERMS:
        raise ValueError(f"output term count exceeds cap {MAX_OUTPUT_TERMS}")
    empty = _empty_packet()
    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "conditional_formal_Hb_axis_Kchi_inlet_j_ge_3_tail",
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "ordinary_Z_labels": list(ORDINARY_Z_LABELS),
        "fields": list(FIELDS),
        "inlet_components": list(INLET_COMPONENTS),
        "inlet": inlet_rows,
        "G0g": inlet_rows,
        "G0g_group": "G0g",
        "inlet_groups": ["G0g"],
        "grouped_inlet_packets_separate": True,
        "grouped_inlet_packets_merged": False,
        "chi": chi_term,
        "chi_z": chi_z_term,
        "chi_Z": chi_z_term,
        "chi_source_operator_ids": common_axis_ids,
        "sc": sc_interval,
        "geometry": geometry,
        "branches": branches,
        "branch_names": list(BRANCHES),
        "all_branches": {branch: value for branch, value in branches.items()},
        "tail": tail,
        "tail_ordinary_Z": tail,
        "tail_response": tail,
        "response": tail,
        "F": tail["F"],
        "Fx": tail["Fx"],
        "finite_prefix": empty,
        "finite": empty,
        "finite_orders_emitted": [],
        "tail_starts_at_j": TAIL_START,
        "ordinary_Z1_product_rule": "chi^j*inlet_Z1 + j*chi^(j-1)*chi_Z*inlet_Z0",
        "tail_majorant": constants,
        "conditions": {
            "hb_strictly_between_zero_and_half": True,
            "hb_condition_is_caller_assumption": True,
            "hb_materialized": False,
            "sc_positive": True,
            "sc_Z_independent": True,
            "geometry_Z_independent": True,
            "chi_scope": "same-source axis Kchi whole-Z interval",
            "chi_source_common": True,
            "chi_whole_Z": True,
            "chi_z_whole_Z": True,
            "inlet_rows_complete": True,
            "ordinary_Z_rows_independent": True,
        },
        "alternative_decomposition": {
            "new_tail_only": True,
            "tail_starts_at_j": TAIL_START,
            "prior_j0_j1_j2_recomputed": False,
            "old_j2_tail_recomputed": False,
            "old_j2_tail_combined": False,
            "old_j2_error_subtracted": False,
            "prior_errors_deleted": False,
            "combined": False,
            "alternative_to_prior_j2": True,
        },
        "scope": {
            "whole_collar_conservative_interval": True,
            "complete_W1": False,
            "full_six_field": False,
            "time_recursion": False,
            "m15": False,
            "ledger_written": False,
            "numeric_function_evaluation": False,
            "midpoint_values_used": False,
            "old_builder_called": False,
            "ancestor_imports": False,
            "receipt_reads": False,
            "saved_data_read": False,
            "formal_hb_not_materialized": True,
        },
        "resource_limits": {
            "max_input_terms": MAX_INPUT_TERMS,
            "max_output_terms": MAX_OUTPUT_TERMS,
            "input_term_count": input_count,
            "output_term_count": output_count,
        },
        "source_ids_preserved": True,
        "source_paths_preserved": True,
        "error_classes_preserved": True,
        "formal_powers_preserved": True,
        "inherited_error_terms_separate": True,
        "tail_terms_individual": True,
        "tail_terms_symmetric_error_only": True,
        "all_tail_terms_error_only": True,
        "all_tail_terms_symmetric": True,
        "signed_tail_value": False,
        "formal_hb_powers_not_materialized": True,
        "no_old_leaf_called": True,
        "no_old_words_recomputed": True,
        "no_old_errors_removed": True,
    }


build_inlet_tail = build
build_tail = build
inlet_tail = build
evaluate = build
