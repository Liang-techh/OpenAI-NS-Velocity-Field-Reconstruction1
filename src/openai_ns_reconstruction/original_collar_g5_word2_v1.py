"""Whole-collar raw-physical ``g5`` word ``n=2`` and independent tail.

This source leaf consumes only caller-owned whole-collar raw physical ``g5``
sectors and certified whole-``Z`` axis potential intervals.  It emits the
previously omitted ``n=2`` local forcing word and a separate ``n>=3``
majorant.  Existing local ``n=0``/``n=1`` rows and an earlier ``j>=2``
envelope are not reconstructed, read, or combined.

``hb`` is a formal source atom.  Coefficient arithmetic evaluates only
directed interval factors in the caller's context ``c``; formal powers stay
in the four-component source power vector.  No ``F0``/``Lambda`` frame
conversion is made for a raw physical ``g5`` source.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any


SCHEMA = "original_collar_g5_word2_v1"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
ORDINARY_Z_LABELS = ("Z0", "Z1")
FIELDS = ("F", "Fx")
WORD_ORDER = 2
TAIL_START = 3
MAX_G5_TERMS = 256
MAX_OUTPUT_TERMS = 4096
MAX_PROVENANCE_ITEMS = 64

FORMAL_HB_POWER_ADDED = {"F": 6, "Fx": 5}
FORMAL_HB_TAIL_POWER_ADDED = {"F": 8, "Fx": 7}

_COEFFICIENT_KEYS = (
    "coefficient", "coefficient_interval", "raw_interval", "interval",
    "enclosure", "box", "value",
)
_FORBIDDEN_SOURCE_MARKERS = (
    "unknown", "old_j2", "old-j2", "prior_j2", "prior-j2",
    "legacy_j2", "legacy-j2", "j>=2", "j>= 2", "reentry", "re-entry",
    "ancestor",
)
_FORBIDDEN_ERROR_MARKERS = (
    "old_j2", "old-j2", "prior_j2", "prior-j2", "legacy_j2", "legacy-j2",
    "j>=2", "j>= 2", "reentry", "re-entry", "previous_error", "old_error",
)

__all__ = [
    "SCHEMA", "SOURCE_POWER_BASIS", "ORDINARY_Z_ORDERS", "ORDINARY_Z_LABELS",
    "FIELDS", "WORD_ORDER", "TAIL_START", "MAX_G5_TERMS", "MAX_OUTPUT_TERMS",
    "FORMAL_HB_POWER_ADDED", "FORMAL_HB_TAIL_POWER_ADDED", "build",
    "build_word2", "build_tail", "evaluate",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _items(value: Any) -> list[Any]:
    if value is None:
        return []
    return list(value) if _is_sequence(value) else [value]


def _clone(value: Any) -> Any:
    try:
        return deepcopy(value)
    except Exception:
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
    """Construct a directed interval from exact caller-owned endpoints."""

    return c.make_mpf((_raw(lower, False), _raw(upper, True)))


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _iv(c: Any, value: Any, label: str = "interval") -> Any:
    """Decode an interval without converting through float or a midpoint."""

    if _is_interval(value):
        result = _make(c, value, value)
    elif isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            if "lower_exact_mpf_tuple" not in value or "upper_exact_mpf_tuple" not in value:
                raise TypeError(f"{label} requires both exact MPF endpoint tuples")
            lower, upper = value["lower_exact_mpf_tuple"], value["upper_exact_mpf_tuple"]
            if not isinstance(lower, (tuple, list)) or not isinstance(upper, (tuple, list)):
                raise TypeError(f"{label} exact endpoints must be MPF tuples")
            if len(lower) != 4 or len(upper) != 4:
                raise TypeError(f"{label} exact endpoints must have four entries")
            result = c.make_mpf((tuple(lower), tuple(upper)))
        else:
            nested = next((value[key] for key in _COEFFICIENT_KEYS if key in value), None)
            if nested is not None and nested is not value:
                result = _iv(c, nested, label)
            else:
                lower = next((value[key] for key in ("lower", "lo", "left") if key in value), None)
                upper = next((value[key] for key in ("upper", "hi", "right") if key in value), None)
                if lower is None or upper is None:
                    raise TypeError(f"{label} is not interval-compatible")
                result = _make(c, lower, upper)
    elif _is_sequence(value) and len(value) == 2:
        result = _make(c, value[0], value[1])
    else:
        try:
            result = c.mpf(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{label} is not interval-compatible") from exc
        if not hasattr(result, "_mpi_"):
            result = _make(c, result, result)
    raw = getattr(result, "_mpi_", None)
    if raw is None:
        raise TypeError(f"{label} did not produce a directed interval")
    lower = c.make_mpf((raw[0], raw[0]))
    upper = c.make_mpf((raw[1], raw[1]))
    if upper < lower:
        raise ValueError(f"{label} has reversed directed endpoints")
    return result


def _bounds(c: Any, value: Any) -> tuple[Any, Any]:
    raw = getattr(value, "_mpi_", None)
    if raw is None:
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


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _iv(c, value, label)
    lower, upper = _lo(c, value), _hi(c, value)
    if upper < 0:
        raise ValueError(f"{label} has a negative upper endpoint")
    if lower < 0:
        lower = _zero(c)
    return _make(c, lower, upper)


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
        result = result * index
    return result


def _pow_int(c: Any, value: Any, exponent: int) -> Any:
    if exponent < 0:
        raise ValueError("integer interval powers must be nonnegative")
    result = _one(c)
    for _ in range(exponent):
        result = result * value
    return result


def _range(c: Any, lower: Any, upper: Any, label: str) -> Any:
    lower, upper = _lo(c, lower), _hi(c, upper)
    if upper < lower:
        raise ValueError(f"{label} has reversed endpoints")
    return _make(c, lower, upper)


def _normalise_powers(raw: Mapping[Any, Any], label: str) -> list[int]:
    powers = raw.get("powers", raw.get("formal_powers", raw.get("source_powers")))
    if powers is None:
        powers = (
            raw.get("hb_power", 0),
            raw.get("F0_power", raw.get("F0_anchor_power", 0)),
            raw.get("Pstar_squared_power", raw.get("Pstar_power", 0)),
            raw.get("Lambda_power", raw.get("lambda_power", 0)),
        )
    if isinstance(powers, Mapping):
        powers = [powers.get(atom, 0) for atom in SOURCE_POWER_BASIS]
    if not _is_sequence(powers) or len(powers) != 4:
        raise TypeError(f"{label}.powers must contain four integer formal powers")
    if any(isinstance(item, bool) or not isinstance(item, int) for item in powers):
        raise TypeError(f"{label}.powers must contain four integer formal powers")
    return [int(item) for item in powers]


def _metadata_items(raw: Mapping[Any, Any], names: Sequence[str]) -> list[Any]:
    values: list[Any] = []
    for name in names:
        if name in raw:
            values.extend(_items(raw[name]))
    for container_name in ("origin", "metadata"):
        container = raw.get(container_name)
        if isinstance(container, Mapping):
            for name in names:
                if name in container:
                    values.extend(_items(container[name]))
    return _unique(values)


def _forbidden(values: Sequence[Any], markers: Sequence[str], label: str) -> None:
    for value in values:
        text = str(value).strip().lower()
        if not text:
            raise ValueError(f"{label} contains an empty source/provenance value")
        if any(marker in text for marker in markers):
            raise ValueError(f"{label} contains a forbidden old/unknown source marker")


def _source_ids(raw: Mapping[Any, Any], label: str) -> list[Any]:
    values = _metadata_items(raw, ("source_ids", "source_id", "sources", "input_id"))
    if not values:
        raise ValueError(f"{label} must retain at least one source_id")
    if len(values) > MAX_PROVENANCE_ITEMS:
        raise ValueError(f"{label} exceeds source ID cap {MAX_PROVENANCE_ITEMS}")
    _forbidden(values, _FORBIDDEN_SOURCE_MARKERS, label + ".source_ids")
    return values


def _source_paths(raw: Mapping[Any, Any], label: str) -> list[Any]:
    values = _metadata_items(raw, ("source_paths", "source_path"))
    if len(values) > MAX_PROVENANCE_ITEMS:
        raise ValueError(f"{label} exceeds source path cap {MAX_PROVENANCE_ITEMS}")
    return values


def _error_classes(raw: Mapping[Any, Any], label: str) -> list[Any]:
    values = _metadata_items(raw, ("error_classes", "error_class", "source_error_classes", "source_error_class"))
    if len(values) > MAX_PROVENANCE_ITEMS:
        raise ValueError(f"{label} exceeds error-class cap {MAX_PROVENANCE_ITEMS}")
    _forbidden(values, _FORBIDDEN_ERROR_MARKERS, label + ".error_classes")
    return values


def _normalise_term(c: Any, raw: Any, label: str) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must be a term mapping")
    coefficient = next((raw[key] for key in _COEFFICIENT_KEYS if key in raw), None)
    if coefficient is None:
        raise ValueError(f"{label} is missing coefficient")
    coefficient_interval = _iv(c, coefficient, label + ".coefficient")
    powers = _normalise_powers(raw, label)
    source_ids = _source_ids(raw, label)
    source_paths = _source_paths(raw, label)
    error_classes = _error_classes(raw, label)
    uncertainty = raw.get("uncertainty", raw.get("is_uncertainty", False))
    error_only = raw.get("error_only", False)
    symmetric_error = raw.get("symmetric_error", error_only)
    # A raw grouped sector is not a signed field merely because it is not
    # marked error-only.  Signedness must be an explicit source assertion.
    signed_field = raw.get("signed_field", False)
    for key, value in (("uncertainty", uncertainty), ("error_only", error_only), ("symmetric_error", symmetric_error), ("signed_field", signed_field)):
        if not isinstance(value, bool):
            raise TypeError(f"{label}.{key} must be boolean")
    if error_only and not symmetric_error:
        raise ValueError(f"{label} error_only requires symmetric_error")
    if error_only and signed_field:
        raise ValueError(f"{label} error_only may not be labelled signed_field")
    if error_only and _lo(c, coefficient_interval) != -_hi(c, coefficient_interval):
        raise ValueError(f"{label} error_only coefficient must be symmetric")
    frame_values = [raw[key] for key in ("frame", "source_frame", "coordinate_frame") if key in raw]
    frame_text = " ".join(str(value).lower() for value in frame_values)
    if any(token in frame_text for token in ("f0 conversion", "lambda conversion", "converted", "normalized")):
        raise ValueError(f"{label} is not in the raw physical g5 source frame")
    result = _clone(dict(raw))
    result.update({
        "coefficient": coefficient_interval, "coefficient_interval": coefficient_interval,
        "powers": list(powers), "formal_powers": list(powers),
        "power_basis": list(SOURCE_POWER_BASIS), "source_ids": list(source_ids),
        "source_paths": list(source_paths), "error_classes": list(error_classes),
        "uncertainty": uncertainty, "error_only": error_only,
        "symmetric_error": symmetric_error, "signed_field": signed_field,
        "input_origin": _clone(raw.get("origin")),
        "input_role": _clone(raw.get("role", raw.get("source_role"))),
    })
    return result


def _order_key(key: Any) -> str | None:
    if key == 0 or (isinstance(key, str) and key.lower() in {"0", "z0", "ordinary_z0", "value"}):
        return "Z0"
    if key == 1 or (isinstance(key, str) and key.lower() in {"1", "z1", "ordinary_z1", "derivative"}):
        return "Z1"
    return None


def _term_list(c: Any, raw: Any, label: str) -> list[dict[str, Any]]:
    if isinstance(raw, Mapping) and any(key in raw for key in _COEFFICIENT_KEYS + ("powers", "formal_powers")):
        values = [raw]
    elif _is_sequence(raw):
        values = list(raw)
    else:
        raise TypeError(f"{label} must be a term sequence")
    if not values:
        raise ValueError(f"{label} must not be empty")
    return [_normalise_term(c, item, f"{label}[{index}]") for index, item in enumerate(values)]


def _normalise_g5(c: Any, g5: Any) -> dict[str, list[dict[str, Any]]]:
    if not isinstance(g5, Mapping):
        raise TypeError("g5 must be a mapping with explicit Z0 and Z1 rows")
    current = g5
    for wrapper in ("g5", "ordinary_Z", "source_Z", "rows"):
        nested = current.get(wrapper) if isinstance(current, Mapping) else None
        if isinstance(nested, Mapping) and nested is not current:
            current = nested
            break
    found: dict[str, Any] = {}
    for key, value in current.items():
        order = _order_key(key)
        if order is not None:
            if order in found:
                raise ValueError(f"g5 supplies {order} twice")
            found[order] = value
    if set(found) != set(ORDINARY_Z_LABELS):
        raise ValueError("g5 requires explicit Z0 and Z1 rows")
    rows = {label: _term_list(c, found[label], "g5." + label) for label in ORDINARY_Z_LABELS}
    if sum(len(rows[label]) for label in ORDINARY_Z_LABELS) > MAX_G5_TERMS:
        raise ValueError(f"g5 term count exceeds cap {MAX_G5_TERMS}")
    return rows


def _axis_marker(value: Any) -> bool:
    text = str(value).replace("_", "").replace("-", "").lower()
    return "axis:wholez" in text or ("axis" in text and "wholez" in text)


def _axis_root(value: Any) -> str | None:
    if not _axis_marker(value):
        return None
    text = str(value).replace("_", "").replace("-", "").lower()
    prefix = text.split("axis:wholez", 1)[0]
    return prefix + "axis:wholez"


def _normalise_axis(c: Any, raw: Any, label: str) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must be a whole-Z certified interval mapping")
    if not any(key in raw for key in _COEFFICIENT_KEYS):
        raise ValueError(f"{label} is missing coefficient")
    term = _normalise_term(c, raw, label)
    if term["powers"] != [0, 0, 0, 0]:
        raise ValueError(f"{label} formal powers must be zero")
    if term["error_only"]:
        raise ValueError(f"{label} must be an axis potential interval, not an error-only packet")
    if raw.get("whole_Z", True) is False or raw.get("whole_z", True) is False:
        raise ValueError(f"{label} must cover whole Z")
    if any(raw.get(key) is False for key in ("certified_whole_Z", "whole_Z_certified", "Z_independent")):
        raise ValueError(f"{label} whole-Z certificate is false")
    scope_text = " ".join(str(raw.get(key, "")).lower() for key in ("scope", "domain", "role", "support"))
    if any(token in scope_text for token in ("lower-collar", "lower collar", "shear", "lower_band", "lower band")):
        raise ValueError(f"{label} may not use a lower-collar shear chi")
    if not any(_axis_marker(value) for value in term["source_ids"]):
        raise ValueError(f"{label} must retain an axis:wholeZ source ID")
    term["axis_whole_Z_certified"] = True
    term["scope"] = raw.get("scope", "whole-Z axis Kchi")
    return term


def _validate_axis_pair(chi: Mapping[str, Any], chi_z: Mapping[str, Any]) -> list[Any]:
    chi_roots = {root for root in (_axis_root(value) for value in chi["source_ids"]) if root is not None}
    chi_z_roots = {root for root in (_axis_root(value) for value in chi_z["source_ids"]) if root is not None}
    common = sorted(chi_roots.intersection(chi_z_roots))
    if not common:
        raise ValueError("chi and chi_z must retain a common whole-Z axis source operator")
    return common


def _validate_context(c: Any) -> None:
    if not all(hasattr(c, name) for name in ("mpf", "make_mpf", "exp", "sqrt")):
        raise TypeError("c must provide directed mpf, make_mpf, exp and sqrt arithmetic")


def _add_powers(*powers: Sequence[int]) -> list[int]:
    return [sum(int(power[index]) for power in powers) for index in range(4)]


def _scaled_powers(power: Sequence[int], multiplier: int) -> list[int]:
    return [int(item) * multiplier for item in power]


def _merge_ids(*terms: Mapping[str, Any]) -> list[Any]:
    return _unique([item for term in terms for item in term.get("source_ids", [])])


def _merge_paths(*terms: Mapping[str, Any]) -> list[Any]:
    return _unique([item for term in terms for item in term.get("source_paths", [])])


def _merge_classes(*terms: Mapping[str, Any]) -> list[Any]:
    return _unique([item for term in terms for item in term.get("error_classes", [])])


def _normalise_sc(c: Any, sc: Any) -> tuple[Any, Any]:
    if sc is None:
        raise TypeError("sc is required")
    raw_record = _clone(sc)
    if isinstance(sc, Mapping):
        for key in ("Z_independent", "z_independent", "geometry_Z_independent", "sc_Z_independent"):
            if key in sc and sc[key] is False:
                raise ValueError(f"sc.{key} certificate is false")
    interval = _positive(c, _iv(c, sc, "sc"), "sc")
    return interval, raw_record


def _geometry(c: Any, sc: Any) -> dict[str, Any]:
    four, sixteen = c.mpf(4), c.mpf(16)
    S = _nonnegative(c, (sc / four) * c.exp(sc / sixteen), "S")
    return {
        "S": S, "S_formula": "S=(sc/4)*exp(sc/16)",
        "x_formula": "x(q)=2*exp(hb*sc*q/8)", "x_end_formula": "x(1)=2*exp(hb*sc/8)",
        "L_bound": "L <= hb*S", "K0_bound": "0 <= K0(x,t) <= x-t",
        "K0x_bound": "K0x <= 1", "Kx_bound": "K0x <= 1",
        "hb_condition": "0 < hb < 1/2", "hb_sc_geometry_Z_independent": True,
        "sc_geometry_Z_independent": True, "geometry_Z_independent": True,
        "formal_hb_not_materialized": True, "midpoint_values_used": False,
    }


def _tail_constants(c: Any, chi: Mapping[str, Any], chi_z: Mapping[str, Any], S: Any) -> dict[str, Any]:
    M = _make(c, _zero(c), _abs_upper(c, chi["coefficient"]))
    Mz = _make(c, _zero(c), _abs_upper(c, chi_z["coefficient"]))
    A = _nonnegative(c, c.mpf(2) * M * _pow_int(c, S, 2), "A=2*M*S^2")
    z = _nonnegative(c, c.sqrt(A) / c.mpf(2), "z=sqrt(A)/2")
    # The exp(sqrt(A)/2) shortcut is not a valid all-A majorant.  The
    # factorial comparison uses exp(A/4), and the m=3 derivative insertion
    # uses its corresponding (1+A/12) factor.
    exp_A_over_4 = c.exp(A / c.mpf(4))
    two_M, two_Mz = c.mpf(2) * M, c.mpf(2) * Mz
    value_F = _nonnegative(c, _pow_int(c, two_M, 3) * _pow_int(c, S, 8) / _factorial(c, 8) * exp_A_over_4, "F tail radius")
    value_Fx = _nonnegative(c, _pow_int(c, two_M, 3) * _pow_int(c, S, 7) / _factorial(c, 7) * exp_A_over_4, "Fx tail radius")
    insertion_factor = (c.mpf(1) + A / c.mpf(12)) * exp_A_over_4
    insertion_F = _nonnegative(c, c.mpf(3) * two_Mz * _pow_int(c, two_M, 2) * _pow_int(c, S, 8) / _factorial(c, 8) * insertion_factor, "F Z1 insertion tail radius")
    insertion_Fx = _nonnegative(c, c.mpf(3) * two_Mz * _pow_int(c, two_M, 2) * _pow_int(c, S, 7) / _factorial(c, 7) * insertion_factor, "Fx Z1 insertion tail radius")
    return {
        "M": M, "Mz": Mz, "A": A, "z": z, "m": 3,
        "exp_A_over_4": exp_A_over_4, "value_F": value_F, "value_Fx": value_Fx,
        "insertion_F": insertion_F, "insertion_Fx": insertion_Fx,
        "F_value_formula": "hb^8*(2M)^3*S^8/8!*exp(A/4)",
        "Fx_value_formula": "hb^7*(2M)^3*S^7/7!*exp(A/4)",
        "F_Z1_insertion_formula": "hb^8*3*(2Mz)*(2M)^2*S^8/8!*(1+A/12)*exp(A/4)",
        "Fx_Z1_insertion_formula": "hb^7*3*(2Mz)*(2M)^2*S^7/7!*(1+A/12)*exp(A/4)",
        "large_A_safe": True,
        "bound_proof": "factorial denominator comparison with exp(A/4), not exp(sqrt(A)/2)",
        "source_radius_policy": "individual source absolute radius; no aggregation",
        "old_j2_bound_reused": False, "old_j2_error_subtracted": False,
    }


def _n2_factor(c: Any, chi: Mapping[str, Any], chi_z: Mapping[str, Any], S: Any, field: str, insertion: bool) -> tuple[Any, str]:
    degree = 6 if field == "F" else 5 if field == "Fx" else None
    if degree is None:
        raise ValueError(f"unknown field {field}")
    taylor = _range(c, _zero(c), _pow_int(c, S, degree) / _factorial(c, degree), "n2 Taylor factor")
    if insertion:
        scale = c.mpf(4) * chi["coefficient"] * chi_z["coefficient"] * taylor
        formula = f"4*chi*chi_Z*[0,S^{degree}/{degree}!]"
    else:
        scale = c.mpf(4) * _pow_int(c, chi["coefficient"], 2) * taylor
        formula = f"4*chi^2*[0,S^{degree}/{degree}!]"
    return scale, formula


def _derived_term(
    source: Mapping[str, Any], *, coefficient: Any, powers: Sequence[int], field: str,
    output_order: int, route: str, word_kind: str, hb_added: int,
    chi_terms: Sequence[Mapping[str, Any]], tail: bool, source_index: int,
    radius_formula: str | None = None, factor_formula: str | None = None,
    insertion_copy: int | None = None,
) -> dict[str, Any]:
    chi_ids = _merge_ids(*chi_terms) if chi_terms else []
    chi_paths = _merge_paths(*chi_terms) if chi_terms else []
    chi_classes = _merge_classes(*chi_terms) if chi_terms else []
    source_ids = _unique(list(source.get("source_ids", [])) + chi_ids)
    classes = _unique(list(source.get("error_classes", [])) + chi_classes)
    derived_class = "n_ge_3_independent_tail" if tail else "n2_signed_expression"
    if derived_class not in classes:
        classes.append(derived_class)
    if tail and "symmetric_error" not in classes:
        classes.append("symmetric_error")
    source_error_only = bool(source.get("error_only", False))
    source_symmetric = bool(source.get("symmetric_error", source_error_only))
    source_uncertainty = bool(source.get("uncertainty", False))
    axis_uncertainty = any(bool(term.get("uncertainty", False)) for term in chi_terms)
    finite_signed = (not tail) and (not source_error_only) and bool(source.get("signed_field", False))
    source_role = _clone(source.get("role", source.get("source_role")))
    if tail:
        derived_role = "symmetric_error_only"
    elif source_error_only:
        # Keep the caller's role visible; an inherited error packet is never
        # relabelled as a signed finite source expression.
        derived_role = source_role if source_role is not None else "inherited_error_only"
    elif finite_signed:
        derived_role = "signed_expression"
    else:
        derived_role = source_role if source_role is not None else "uncertainty_unverified"
    result = _clone(dict(source))
    result.update({
        "coefficient": coefficient, "coefficient_interval": coefficient,
        "powers": list(powers), "formal_powers": list(powers),
        "power_basis": list(SOURCE_POWER_BASIS), "formal_power_shift": [hb_added, 0, 0, 0],
        "source_ids": source_ids, "g5_source_ids": list(source.get("source_ids", [])),
        "g5_source_paths": list(source.get("source_paths", [])),
        "chi_source_ids": chi_ids, "chi_source_paths": chi_paths,
        "chi_source_classes": chi_classes,
        "parent_source_ids": list(source.get("source_ids", [])),
        "parent_source_paths": list(source.get("source_paths", [])),
        "parent_origin": _clone(source.get("origin")),
        "parent_error_classes": list(source.get("error_classes", [])),
        "source_uncertainty": source_uncertainty,
        "axis_uncertainty": axis_uncertainty,
        "source_role": source_role,
        "error_classes": classes,
        "uncertainty": True if tail else (source_uncertainty or axis_uncertainty),
        "error_only": True if tail else source_error_only,
        "symmetric_error": True if tail else source_symmetric,
        "signed_field": False if tail else finite_signed,
        "derived_role": derived_role,
        "source_role_preserved": True,
        "coefficient_is_signed_expression": finite_signed,
        "coefficient_is_function_value": False,
        "known_error_upper_is_function_value": False,
        "midpoint_values_used": False, "coefficient_contains_formal_hb": False,
        "word_order": TAIL_START if tail else WORD_ORDER,
        "tail": tail, "tail_starts_at_j": TAIL_START,
        "field": field, "ordinary_Z_order": output_order,
        "source_order": 0 if "g5_Z0" in route else 1,
        "route": route, "word_kind": word_kind,
        "operator": "sum(n>=3)(-2*chi)^n*G0^(n+1)*g5" if tail else "(-2*chi)^2*G0^3*g5",
        "signed_source_coefficient_preserved": True, "signed_tail_value": False,
        "old_j2_tail_combined": False, "old_j2_error_subtracted": False,
        "entered_212_ledger": False, "source_frame": "raw_physical_g5",
        "F0_Lambda_conversion": False, "source_term_index": source_index,
        "chi_formal_powers": [list(term.get("powers", (0, 0, 0, 0))) for term in chi_terms],
        "chi_factor_count": 3 if tail else 2, "source_paths_preserved": True,
    })
    if insertion_copy is not None:
        result["insertion_copy"] = insertion_copy
        result["insertion_copies_total"] = 2
        result["derivative_factor_is_split"] = True
        result["derivative_factor"] = 1
    if tail:
        result.update({"signed_expression_interval": False, "inherited_error_only": source_error_only, "tail_radius_is_symmetric": True, "error_only_tail": True})
    else:
        result.update({"signed_expression_interval": finite_signed, "inherited_error_only": source_error_only})
    if radius_formula is not None:
        result["radius_formula"] = radius_formula
    if factor_formula is not None:
        result["factor_formula"] = factor_formula
    return result


def _build_n2(c: Any, rows: Mapping[str, Sequence[Mapping[str, Any]]], chi: Mapping[str, Any], chi_z: Mapping[str, Any], S: Any) -> dict[str, dict[str, list[dict[str, Any]]]]:
    result = {field: {label: [] for label in ORDINARY_Z_LABELS} for field in FIELDS}
    for index, source in enumerate(rows["Z0"]):
        for field in FIELDS:
            scale, formula = _n2_factor(c, chi, chi_z, S, field, False)
            powers = _add_powers(source["powers"], _scaled_powers(chi["powers"], 2))
            powers[0] += FORMAL_HB_POWER_ADDED[field]
            result[field]["Z0"].append(_derived_term(
                source, coefficient=source["coefficient"] * scale, powers=powers, field=field,
                output_order=0, route="chi^2*g5_Z0", word_kind="n2_value_Z0",
                hb_added=FORMAL_HB_POWER_ADDED[field], chi_terms=(chi,), tail=False, source_index=index,
                factor_formula=formula,
            ))
    for index, source in enumerate(rows["Z1"]):
        for field in FIELDS:
            scale, formula = _n2_factor(c, chi, chi_z, S, field, False)
            powers = _add_powers(source["powers"], _scaled_powers(chi["powers"], 2))
            powers[0] += FORMAL_HB_POWER_ADDED[field]
            result[field]["Z1"].append(_derived_term(
                source, coefficient=source["coefficient"] * scale, powers=powers, field=field,
                output_order=1, route="chi^2*g5_Z1", word_kind="n2_value_Z1",
                hb_added=FORMAL_HB_POWER_ADDED[field], chi_terms=(chi,), tail=False, source_index=index,
                factor_formula=formula,
            ))
    # d_Z(chi^2) is represented as two provenance-preserving factors.  Each
    # record carries 4*chi*chi_Z; their sum is the required factor 2, but no
    # provenance or signed source interval is merged.
    for index, source in enumerate(rows["Z0"]):
        for field in FIELDS:
            scale, formula = _n2_factor(c, chi, chi_z, S, field, True)
            powers = _add_powers(source["powers"], chi["powers"], chi_z["powers"])
            powers[0] += FORMAL_HB_POWER_ADDED[field]
            for copy_index, route in enumerate((
                "chi_outer*chi_Z_inner*g5_Z0", "chi_inner*chi_Z_outer*g5_Z0",
            ), start=1):
                result[field]["Z1"].append(_derived_term(
                    source, coefficient=source["coefficient"] * scale, powers=powers, field=field,
                    output_order=1, route=route, word_kind="n2_chi_Z_insertion",
                    hb_added=FORMAL_HB_POWER_ADDED[field], chi_terms=(chi, chi_z), tail=False,
                    source_index=index, insertion_copy=copy_index,
                    factor_formula=formula,
                ))
    return result


def _build_tail(c: Any, rows: Mapping[str, Sequence[Mapping[str, Any]]], chi: Mapping[str, Any], chi_z: Mapping[str, Any], constants: Mapping[str, Any]) -> dict[str, dict[str, list[dict[str, Any]]]]:
    result = {field: {label: [] for label in ORDINARY_Z_LABELS} for field in FIELDS}
    for index, source in enumerate(rows["Z0"]):
        for field in FIELDS:
            majorant = constants["value_F"] if field == "F" else constants["value_Fx"]
            powers = _add_powers(source["powers"], _scaled_powers(chi["powers"], 3))
            powers[0] += FORMAL_HB_POWER_TAIL_POWER_ADDED[field]
            result[field]["Z0"].append(_derived_term(
                source, coefficient=_symmetric(c, _abs_upper(c, source["coefficient"]) * majorant), powers=powers,
                field=field, output_order=0, route="chi^3*g5_Z0_value_tail", word_kind="n_ge_3_value_tail",
                hb_added=FORMAL_HB_POWER_TAIL_ADDED[field], chi_terms=(chi,), tail=True, source_index=index,
                radius_formula=constants["F_value_formula"] if field == "F" else constants["Fx_value_formula"],
            ))
    for index, source in enumerate(rows["Z1"]):
        for field in FIELDS:
            majorant = constants["value_F"] if field == "F" else constants["value_Fx"]
            powers = _add_powers(source["powers"], _scaled_powers(chi["powers"], 3))
            powers[0] += FORMAL_HB_POWER_TAIL_ADDED[field]
            result[field]["Z1"].append(_derived_term(
                source, coefficient=_symmetric(c, _abs_upper(c, source["coefficient"]) * majorant), powers=powers,
                field=field, output_order=1, route="chi^3*g5_Z1_value_tail", word_kind="n_ge_3_value_tail",
                hb_added=FORMAL_HB_POWER_TAIL_ADDED[field], chi_terms=(chi,), tail=True, source_index=index,
                radius_formula=constants["F_value_formula"] if field == "F" else constants["Fx_value_formula"],
            ))
    for index, source in enumerate(rows["Z0"]):
        for field in FIELDS:
            majorant = constants["insertion_F"] if field == "F" else constants["insertion_Fx"]
            powers = _add_powers(source["powers"], _scaled_powers(chi["powers"], 2), chi_z["powers"])
            powers[0] += FORMAL_HB_POWER_TAIL_ADDED[field]
            result[field]["Z1"].append(_derived_term(
                source, coefficient=_symmetric(c, _abs_upper(c, source["coefficient"]) * majorant), powers=powers,
                field=field, output_order=1, route="3*(2Mz)*(2M)^2*g5_Z0_insertion_tail", word_kind="n_ge_3_chi_Z_insertion_tail",
                hb_added=FORMAL_HB_POWER_TAIL_ADDED[field], chi_terms=(chi, chi_z), tail=True, source_index=index,
                radius_formula=constants["F_Z1_insertion_formula"] if field == "F" else constants["Fx_Z1_insertion_formula"],
            ))
    return result


def _term_count(packet: Mapping[str, Mapping[str, Sequence[Any]]]) -> int:
    return sum(len(packet[field][label]) for field in FIELDS for label in ORDINARY_Z_LABELS)


def build(c: Any, g5: Any, chi: Any, chi_z: Any, sc: Any) -> dict[str, Any]:
    """Build the ``n=2`` word and independent ``n>=3`` whole-collar tail.

    The positional interface is intentional: ``build(c, g5, chi, chi_z,
    sc)`` makes the raw source frame and both whole-Z axis inputs explicit.
    No numerical ``hb`` argument is accepted or materialized.
    """

    _validate_context(c)
    rows = _normalise_g5(c, g5)
    chi_term = _normalise_axis(c, chi, "chi")
    chi_z_term = _normalise_axis(c, chi_z, "chi_z")
    common_axis_roots = _validate_axis_pair(chi_term, chi_z_term)
    sc_interval, sc_record = _normalise_sc(c, sc)
    geometry = _geometry(c, sc_interval)
    constants = _tail_constants(c, chi_term, chi_z_term, geometry["S"])
    n2 = _build_n2(c, rows, chi_term, chi_z_term, geometry["S"])
    tail = _build_tail(c, rows, chi_term, chi_z_term, constants)
    output_count = _term_count(n2) + _term_count(tail)
    if output_count > MAX_OUTPUT_TERMS:
        raise ValueError(f"output term count exceeds cap {MAX_OUTPUT_TERMS}")
    n2_terms = [term for field in FIELDS for label in ORDINARY_Z_LABELS for term in n2[field][label]]
    contains_signed_expression = any(term.get("coefficient_is_signed_expression") is True for term in n2_terms)
    branch_meta = {
        "n2": {
            "word_order": WORD_ORDER, "tail_starts_at_j": TAIL_START,
            "operator": "(-2*chi)^2*G0^3*g5", "source_frame": "raw_physical_g5",
            "formal_hb_powers": dict(FORMAL_HB_POWER_ADDED),
            "F_factor": "hb^6*4*chi^2*[0,S^6/6!]", "Fx_factor": "hb^5*4*chi^2*[0,S^5/5!]",
            "ordinary_Z1": "4*(chi^2*g5_Z1 + 2*chi*chi_Z*g5_Z0)",
            "ordinary_Z1_value_route": "chi^2*g5_Z1",
            "ordinary_Z1_insertion_route": "2*chi*chi_Z*g5_Z0",
            "Taylor_factor_not_divided_again": True, "source_terms_individual": True,
            "chi_Z_insertion_paths_separate": True, "response": n2,
        },
        "n_ge_3_tail": {
            "tail_starts_at_j": TAIL_START,
            "operator": "sum(n>=3)(-2*chi)^n*G0^(n+1)*g5",
            "formal_hb_powers": dict(FORMAL_HB_TAIL_POWER_ADDED), "response": tail,
            "all_terms_error_only": True, "all_terms_symmetric": True,
            "source_terms_individual": True,
        },
    }
    alternative = {
        "n2_plus_n_ge_3": True,
        "new_n2_and_tail_are_alternative_to_prior_j2": True,
        "prior_local_g5_n0_n1_recomputed": False, "prior_j2_tail_recomputed": False,
        "prior_j2_tail_combined": False, "old_j2_bound_reused": False,
        "old_j2_error_subtracted": False, "entered_212_ledger": False,
    }
    scope = {
        "whole_collar": True, "whole_Z_axis_potential": True,
        "lower_collar_shear_chi_accepted": False, "raw_physical_g5_source_frame": True,
        "F0_Lambda_conversion": False, "formal_hb_not_materialized": True,
        "midpoint_values_used": False, "same_source_operator_conditioned": True,
        "source_owner_proved": False, "old_helper_called": False, "ancestor_imports": False,
        "saved_data_read": False, "ledger_written": False, "complete_W1": False,
        "complete_six_field": False, "m15": False, "time_recursion": False,
    }
    return {
        "schema": SCHEMA,
        "status": "conditional_raw_physical_g5_word2_and_independent_n_ge_3_tail",
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS), "ordinary_Z_labels": list(ORDINARY_Z_LABELS),
        "fields": list(FIELDS), "g5": rows, "g5_jets": rows,
        "chi": chi_term, "chi_z": chi_z_term, "chi_Z": chi_z_term,
        "chi_source_operator_roots": common_axis_roots, "sc": sc_interval, "sc_record": sc_record,
        "geometry": geometry,
        "kernel": {
            "G0_formula": "G0", "G0_bound": "0 <= K0(x,t) <= x-t", "K0x_bound": "K0x <= 1",
            "physical_x_formula": "x(q)=2*exp(hb*sc*q/8)", "L_bound": "L <= hb*S",
        },
        "n2": n2, "word2": n2, "response": n2, "finite_prefix": n2,
        "tail": tail, "n_ge_3_tail": tail, "tail_response": tail,
        "branches": branch_meta, "n2_metadata": branch_meta["n2"],
        "tail_metadata": branch_meta["n_ge_3_tail"], "tail_majorant": constants,
        "alternative_decomposition": alternative, "scope": scope,
        "conditions": {
            "hb_strictly_between_zero_and_half": True, "hb_condition_is_caller_assumption": True,
            "hb_materialized": False, "sc_positive": True, "sc_geometry_Z_independent": True,
            "chi_whole_Z_certified": True, "chi_z_whole_Z_certified": True,
            "common_axis_source_operator": True, "g5_rows_complete": True,
        },
        "resource_limits": {
            "max_g5_terms": MAX_G5_TERMS, "max_output_terms": MAX_OUTPUT_TERMS,
            "input_g5_term_count": sum(len(rows[label]) for label in ORDINARY_Z_LABELS),
            "output_term_count": output_count,
        },
        "source_ids_preserved": True, "source_paths_preserved": True, "origin_preserved": True,
        "error_classes_preserved": True, "uncertainty_preserved_for_n2": True,
        "error_only_input_preserved_for_n2": True, "error_only_tail_outputs": True,
        "symmetric_tail_outputs": True, "signed_expression_coefficients": contains_signed_expression,
        "contains_signed_expression_terms": contains_signed_expression,
        "coefficients_are_function_values": False, "no_old_error_mixing": True,
        "no_unknown_source": True, "no_hb_argument": True, "no_ancestor_helper_called": True,
        "receipt_reads": False, "saved_data_read": False, "appended_to_212_ledger": False,
    }


build_word2 = build
build_tail = build
evaluate = build

