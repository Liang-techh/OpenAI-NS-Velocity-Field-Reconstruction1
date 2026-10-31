"""Scale-preserving, one-cell ``m15`` ``W1`` density algebra.

This is deliberately a small algebra leaf.  It does not integrate a density,
construct an original source, or import another reconstruction module.  The
caller supplies exact directed ranges for one fixed ``rho`` cell and for the
ordinary ``Z0``/``Z1`` rows of ``U0``, ``Delta``, ``Phi``, ``F1``, ``Uz1`` and
``L = F0_Z/F0``.  A range backend is supplied by the caller; the bundled
``FractionIntervalBackend`` is only a convenient exact fallback.

The two signed densities are kept in coefficient space:

``axial = epsilon * 2 * (U0 + Delta) * Uz1``

``swirl = epsilon**2 * F0 * (-2 * rho * Phi * F1)``

The ordinary ``Z1`` rows use the product rule, with ``rho_Z = epsilon_Z = 0``
as explicit caller premises.  In particular, this module never divides by
``epsilon`` and never materializes either ``epsilon`` or ``F0``.  Every input
finite range and additive tail remains a separately labelled product term;
terms are not aggregated into a value or an error cap.

Each supplied Z1 row denotes the ordinary derivative of the full supplied
field, already expressed in its own formal powers.  In particular, if an
input carries an F0 power, its Z1 row must already include that input's F0
derivative.  The additional L term below differentiates the leading F0 only.
Optional original N1p rows remain a separate pressure channel.  A caller may
express an internal Lambda factor as a formal epsilon power -1; powers are
combined symbolically and no epsilon inverse is numerically evaluated.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import Any


SCHEMA = "original_core_scaled_m15_density_v1"
COMPONENT = "original_core_scaled_m15_density"
ORDINARY_Z_ORDERS = (0, 1)
ORDINARY_Z_LABELS = ("Z0", "Z1")
FIELDS = ("U0", "Delta", "Phi", "F1", "Uz1", "L")
FORMAL_POWER_BASIS = ("epsilon", "F0")


class DensityInputError(ValueError):
    """Raised when a density packet is not an exact, admissible range packet."""


class ExactRangeBackendError(TypeError):
    """Raised when the caller's range backend cannot provide exact intervals."""


@dataclass(frozen=True)
class RationalInterval:
    """An exact directed interval used by :class:`FractionIntervalBackend`."""

    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        if type(self.lower) is not Fraction or type(self.upper) is not Fraction:
            raise TypeError("RationalInterval endpoints must be Fraction values")
        if self.lower > self.upper:
            raise ValueError("RationalInterval lower endpoint exceeds upper endpoint")

    @property
    def lo(self) -> Fraction:
        return self.lower

    @property
    def hi(self) -> Fraction:
        return self.upper

    def __iter__(self):
        yield self.lower
        yield self.upper

    def __len__(self) -> int:
        return 2

    def __getitem__(self, index: int) -> Fraction:
        if index == 0:
            return self.lower
        if index == 1:
            return self.upper
        raise IndexError(index)

    def as_tuple(self) -> tuple[Fraction, Fraction]:
        return self.lower, self.upper


def _exact(value: Any, label: str) -> Fraction:
    """Decode an exact scalar without accepting floating point values."""

    if isinstance(value, bool):
        raise DensityInputError(f"{label} must not be bool")
    if isinstance(value, float):
        raise DensityInputError(f"{label} must be an exact rational, not float")
    if type(value) is Fraction:
        return value
    if type(value) is int:
        return Fraction(value)
    if isinstance(value, str):
        text = value.strip()
        if not text:
            raise DensityInputError(f"{label} is empty")
        try:
            return Fraction(text)
        except (TypeError, ValueError, ZeroDivisionError) as exc:
            raise DensityInputError(f"{label} must be an exact rational") from exc
    raise DensityInputError(f"{label} must be an exact rational")


class FractionIntervalBackend:
    """Minimal exact range backend for synthetic callers and fixtures.

    A production caller can provide an object with the same ``interval``,
    ``add``, ``mul`` and ``scale`` methods.  The leaf only asks the backend for
    directed range arithmetic; it does not choose a midpoint or a precision.
    """

    exact = True

    def interval(self, lower: Any, upper: Any | None = None) -> RationalInterval:
        lo = _exact(lower, "interval.lower")
        hi = lo if upper is None else _exact(upper, "interval.upper")
        return RationalInterval(lo, hi)

    make_interval = interval

    def endpoints(self, value: RationalInterval) -> tuple[Fraction, Fraction]:
        return value.lower, value.upper

    def add(self, left: RationalInterval, right: RationalInterval) -> RationalInterval:
        return RationalInterval(left.lower + right.lower, left.upper + right.upper)

    def mul(self, left: RationalInterval, right: RationalInterval) -> RationalInterval:
        values = (
            left.lower * right.lower,
            left.lower * right.upper,
            left.upper * right.lower,
            left.upper * right.upper,
        )
        return RationalInterval(min(values), max(values))

    def scale(self, factor: Any, value: RationalInterval) -> RationalInterval:
        f = _exact(factor, "interval scale")
        if f >= 0:
            return RationalInterval(value.lower * f, value.upper * f)
        return RationalInterval(value.upper * f, value.lower * f)

    def neg(self, value: RationalInterval) -> RationalInterval:
        return self.scale(-1, value)


_MISSING = object()


def _key_name(value: Any) -> str:
    return str(value).strip().lower().replace("-", "_").replace(" ", "_")


def _first(mapping: Mapping[Any, Any], *names: str) -> Any:
    for name in names:
        if name in mapping:
            return mapping[name]
    normalised = {_key_name(key): key for key in mapping}
    for name in names:
        key = normalised.get(_key_name(name))
        if key is not None:
            return mapping[key]
    return _MISSING


def _has(mapping: Mapping[Any, Any], *names: str) -> bool:
    return _first(mapping, *names) is not _MISSING


def _items(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    return [value]


def _stable(value: Any) -> tuple[str, str]:
    try:
        hash(value)
    except TypeError:
        return type(value).__name__, repr(value)
    return type(value).__name__, repr(value)


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    seen: set[tuple[str, str]] = set()
    for value in values:
        marker = _stable(value)
        if marker not in seen:
            seen.add(marker)
            result.append(value)
    return result


def _reject_claims(value: Any, path: str = "input", seen: set[int] | None = None) -> None:
    """Reject point/midpoint assertions and already-integrated inputs."""

    if seen is None:
        seen = set()
    marker = id(value)
    if marker in seen:
        return
    if isinstance(value, Mapping):
        seen.add(marker)
        for raw_key, child in value.items():
            key = _key_name(raw_key)
            if key in {
                "point",
                "is_point",
                "point_value",
                "midpoint",
                "is_midpoint",
                "midpoint_claim",
                "point_claim",
            } and child not in (False, None):
                raise DensityInputError(f"{path}.{raw_key} point or midpoint claims are not accepted")
            if key in {
                "invalid",
                "invalid_flag",
                "is_invalid",
                "valid",
                "is_valid",
            } and child not in (False, None, True):
                raise DensityInputError(f"{path}.{raw_key} has an invalid flag")
            if key in {"invalid", "invalid_flag", "is_invalid"} and child is True:
                raise DensityInputError(f"{path}.{raw_key} is invalid")
            if key in {"valid", "is_valid"} and child is False:
                raise DensityInputError(f"{path}.{raw_key} is false")
            if key in {
                "integration_performed",
                "integrated",
            } and child is True:
                raise DensityInputError(f"{path}.{raw_key} makes an unverified claim")
            _reject_claims(child, f"{path}.{raw_key}", seen)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        seen.add(marker)
        for index, child in enumerate(value):
            _reject_claims(child, f"{path}[{index}]", seen)


def _extract_endpoints(value: Any, backend: Any) -> tuple[Any, Any] | None:
    endpoint_method = getattr(backend, "endpoints", None)
    if callable(endpoint_method):
        try:
            pair = endpoint_method(value)
            if isinstance(pair, Sequence) and not isinstance(pair, (str, bytes, bytearray)) and len(pair) == 2:
                return pair[0], pair[1]
        except Exception:
            pass
    if isinstance(value, RationalInterval):
        return value.lower, value.upper
    if isinstance(value, Mapping):
        lower = _first(value, "lower", "lo", "left", "a")
        upper = _first(value, "upper", "hi", "right", "b")
        if lower is not _MISSING and upper is not _MISSING:
            return lower, upper
    for lower_name, upper_name in (("lower", "upper"), ("lo", "hi"), ("a", "b")):
        if hasattr(value, lower_name) and hasattr(value, upper_name):
            return getattr(value, lower_name), getattr(value, upper_name)
    raw = getattr(value, "_mpi_", None)
    if isinstance(raw, Sequence) and len(raw) == 2:
        return raw[0], raw[1]
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)) and len(value) == 2:
        return value[0], value[1]
    return None


class _Backend:
    """Small adapter around the caller's exact-range backend."""

    def __init__(self, backend: Any | None) -> None:
        self.raw = FractionIntervalBackend() if backend is None else backend
        if not any(callable(getattr(self.raw, name, None)) for name in ("interval", "make_interval")):
            raise ExactRangeBackendError("backend must provide interval or make_interval")
        if getattr(self.raw, "exact", True) is False:
            raise ExactRangeBackendError("an exact-range backend is required")

    def endpoints(self, value: Any, label: str) -> tuple[Any, Any]:
        pair = _extract_endpoints(value, self.raw)
        if pair is None:
            raise DensityInputError(f"{label} is not a directed interval")
        lo, hi = pair
        if isinstance(lo, float) or isinstance(hi, float) or isinstance(lo, bool) or isinstance(hi, bool):
            raise DensityInputError(f"{label} endpoints must be exact, not float or bool")
        try:
            if lo > hi:
                raise DensityInputError(f"{label} interval is reversed")
        except TypeError as exc:
            raise ExactRangeBackendError(f"backend endpoints for {label} are not ordered") from exc
        return lo, hi

    def interval(self, value: Any, label: str) -> Any:
        if isinstance(value, Mapping):
            nested = _first(value, "interval", "bounds", "range", "finite", "value")
            if nested is not _MISSING and nested is not value:
                value = nested
        pair = _extract_endpoints(value, self.raw)
        if pair is None:
            raise DensityInputError(f"{label} must supply lower and upper exact endpoints")
        lo, hi = pair
        # Input endpoint scalars are normalized through Fraction so a caller
        # cannot accidentally smuggle a binary float into an exact backend.
        lo = _exact(lo, f"{label}.lower")
        hi = _exact(hi, f"{label}.upper")
        if lo > hi:
            raise DensityInputError(f"{label} interval is reversed")
        make = getattr(self.raw, "interval", None) or getattr(self.raw, "make_interval", None)
        try:
            result = make(lo, hi)
        except Exception as exc:
            raise ExactRangeBackendError(f"backend could not construct {label}") from exc
        self.endpoints(result, label)
        return result

    def add(self, left: Any, right: Any, label: str) -> Any:
        operation = getattr(self.raw, "add", None)
        try:
            result = operation(left, right) if callable(operation) else left + right
        except Exception as exc:
            raise ExactRangeBackendError(f"backend could not add {label}") from exc
        self.endpoints(result, label)
        return result

    def mul(self, left: Any, right: Any, label: str) -> Any:
        operation = getattr(self.raw, "mul", None)
        try:
            result = operation(left, right) if callable(operation) else left * right
        except Exception as exc:
            raise ExactRangeBackendError(f"backend could not multiply {label}") from exc
        self.endpoints(result, label)
        return result

    def scale(self, factor: int | Fraction, value: Any, label: str) -> Any:
        operation = getattr(self.raw, "scale", None)
        try:
            result = operation(factor, value) if callable(operation) else value * factor
        except Exception as exc:
            raise ExactRangeBackendError(f"backend could not scale {label}") from exc
        self.endpoints(result, label)
        return result

    def serial(self, value: Any) -> dict[str, Any]:
        serializer = getattr(self.raw, "serialize", None)
        if callable(serializer):
            result = serializer(value)
            if isinstance(result, Mapping):
                return {
                    "lower": result.get("lower", result.get("lo", result.get("left"))),
                    "upper": result.get("upper", result.get("hi", result.get("right"))),
                }
        lo, hi = self.endpoints(value, "result interval")
        return {"lower": lo, "upper": hi}


def _provenance(value: Any, *, seen: set[int] | None = None) -> dict[str, Any]:
    """Collect source/term/error provenance without dropping old metadata."""

    if seen is None:
        seen = set()
    marker = id(value)
    if marker in seen:
        return {"source_ids": [], "term_ids": [], "source_paths": [], "error_classes": [], "uncertainty": False, "scale": []}
    if isinstance(value, Mapping):
        seen.add(marker)
        source_ids: list[Any] = []
        term_ids: list[Any] = []
        source_paths: list[Any] = []
        error_classes: list[Any] = []
        scales: list[Any] = []
        for key in ("source_ids", "sources", "source_id", "source", "input_source_ids"):
            raw = _first(value, key)
            if raw is not _MISSING:
                source_ids.extend(_items(raw))
        for key in ("term_ids", "term_id", "parent_id", "input_term_ids", "id"):
            raw = _first(value, key)
            if raw is not _MISSING:
                term_ids.extend(_items(raw))
        raw_paths = _first(value, "source_paths", "source_path")
        if raw_paths is not _MISSING:
            source_paths.extend(_items(raw_paths))
        for key in ("error_classes", "error_class", "uncertainty_class"):
            raw = _first(value, key)
            if raw is not _MISSING:
                error_classes.extend(_items(raw))
        for key in ("original_scale", "original_scale_metadata", "scale_metadata", "scale", "scale_powers"):
            raw = _first(value, key)
            if raw is not _MISSING:
                scales.append({"key": key, "value": raw})
        uncertain = bool(_first(value, "uncertainty", "is_uncertainty") not in (_MISSING, False, None))
        for key in ("origin", "metadata", "provenance", "parents"):
            child = _first(value, key)
            if child is _MISSING or child is value:
                continue
            nested = _provenance(child, seen=seen)
            source_ids.extend(nested["source_ids"])
            term_ids.extend(nested["term_ids"])
            source_paths.extend(nested["source_paths"])
            error_classes.extend(nested["error_classes"])
            scales.extend(nested["scale"])
            uncertain = uncertain or nested["uncertainty"]
        return {
            "source_ids": _unique(source_ids),
            "term_ids": _unique(term_ids),
            "source_paths": _unique(source_paths),
            "error_classes": _unique(error_classes),
            "uncertainty": uncertain,
            "scale": scales,
        }
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        seen.add(marker)
        result = {"source_ids": [], "term_ids": [], "source_paths": [], "error_classes": [], "uncertainty": False, "scale": []}
        for child in value:
            nested = _provenance(child, seen=seen)
            for name in ("source_ids", "term_ids", "source_paths", "error_classes", "scale"):
                result[name].extend(nested[name])
            result["uncertainty"] = result["uncertainty"] or nested["uncertainty"]
        for name in ("source_ids", "term_ids", "source_paths", "error_classes"):
            result[name] = _unique(result[name])
        return result
    return {"source_ids": [], "term_ids": [], "source_paths": [], "error_classes": [], "uncertainty": False, "scale": []}


@dataclass(frozen=True)
class _Component:
    interval: Any
    field: str
    ordinary_order: int
    kind: str
    label: str
    tail_id: Any
    source_ids: tuple[Any, ...]
    term_ids: tuple[Any, ...]
    source_paths: tuple[Any, ...]
    error_classes: tuple[Any, ...]
    uncertainty: bool
    scale_metadata: tuple[Any, ...]
    formal_powers: tuple[int, int] | None
    raw_metadata: Any


def _component_metadata(value: Any, inherited: Any, field: str, order: int) -> dict[str, Any]:
    info = _provenance(value)
    parent = _provenance(inherited)
    return {
        "source_ids": tuple(_unique(parent["source_ids"] + info["source_ids"])),
        "term_ids": tuple(_unique(parent["term_ids"] + info["term_ids"])),
        "source_paths": tuple(_unique(parent["source_paths"] + info["source_paths"])),
        "error_classes": tuple(_unique(parent["error_classes"] + info["error_classes"])),
        "uncertainty": bool(parent["uncertainty"] or info["uncertainty"]),
        "scale": tuple(parent["scale"] + info["scale"]),
        "field": field,
        "ordinary_Z_order": order,
    }


def _formal_powers(value: Any, label: str) -> tuple[int, int] | None:
    """Read caller-supplied epsilon/F0 powers for an optional pressure row."""

    if not isinstance(value, Mapping):
        return None
    raw = _first(value, "formal_powers", "scale_powers", "powers")
    if raw is _MISSING:
        return None
    if isinstance(raw, Mapping):
        epsilon = _first(raw, "epsilon", "epsilon_power", "e")
        f0 = _first(raw, "F0", "F0_power", "f0", "f")
        if epsilon is _MISSING or f0 is _MISSING:
            raise DensityInputError(f"{label}.formal_powers must name epsilon and F0")
        raw = (epsilon, f0)
    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes, bytearray)) or len(raw) != 2:
        raise DensityInputError(f"{label}.formal_powers must contain two integer powers")
    powers: list[int] = []
    for index, power in enumerate(raw):
        if isinstance(power, bool) or not isinstance(power, int):
            raise DensityInputError(f"{label}.formal_powers[{index}] must be an integer")
        powers.append(int(power))
    return powers[0], powers[1]


def _component_interval(raw: Any, backend: _Backend, label: str) -> Any:
    if isinstance(raw, Mapping):
        nested = _first(raw, "interval", "bounds", "range", "coefficient", "value")
        if nested is not _MISSING and nested is not raw:
            raw = nested
    return backend.interval(raw, label)


def _component(raw: Any, *, kind: str, field: str, order: int, index: int, inherited: Any, backend: _Backend, label: str) -> _Component:
    if isinstance(raw, Mapping):
        nested = _first(raw, "interval", "bounds", "range", "coefficient", "value")
        interval_raw = raw if nested is _MISSING else nested
    else:
        interval_raw = raw
    interval = _component_interval(interval_raw, backend, label)
    info = _component_metadata(raw, inherited, field, order)
    own_id = _first(raw, "tail_id", "component_id", "term_id", "parent_id", "id") if isinstance(raw, Mapping) else _MISSING
    if own_id is _MISSING:
        own_id = index
    if kind == "finite":
        component_label = "finite" if index == 0 else f"finite:{own_id}"
        tail_id = None
    else:
        component_label = f"tail:{own_id}"
        tail_id = own_id
    powers = _formal_powers(raw, label)
    if powers is None:
        powers = _formal_powers(inherited, label)
    for scale in info["scale"]:
        value = scale.get("value") if isinstance(scale, Mapping) else None
        if isinstance(value, Mapping) and any(value.get(key, 0) != 0 for key in ("epsilon_power", "F0_power")) and powers is None:
            raise DensityInputError(f"{label} has factored scale metadata but no explicit formal_powers")
    return _Component(
        interval=interval,
        field=field,
        ordinary_order=order,
        kind=kind,
        label=component_label,
        tail_id=tail_id,
        source_ids=info["source_ids"],
        term_ids=info["term_ids"],
        source_paths=info["source_paths"],
        error_classes=info["error_classes"],
        uncertainty=info["uncertainty"] or kind == "tail",
        scale_metadata=info["scale"],
        formal_powers=powers,
        raw_metadata=raw,
    )


def _as_component_rows(raw: Any, *, kind: str, field: str, order: int, inherited: Any, backend: _Backend, label: str, allow_empty: bool = False) -> list[_Component]:
    if raw is _MISSING or raw is None:
        if allow_empty:
            return []
        raise DensityInputError(f"{label} is missing")
    if isinstance(raw, Mapping):
        nested = _first(raw, "terms", "components")
        if nested is not _MISSING:
            raw = nested
    if isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
        # A pair of exact endpoint scalars is one interval, not two terms.
        if len(raw) == 2:
            try:
                _exact(raw[0], label + ".lower")
                _exact(raw[1], label + ".upper")
            except DensityInputError:
                values = list(raw)
            else:
                values = [raw]
        else:
            values = list(raw)
    else:
        values = [raw]
    if not values and not allow_empty:
        raise DensityInputError(f"{label} must contain at least one interval")
    return [
        _component(value, kind=kind, field=field, order=order, index=index, inherited=inherited, backend=backend, label=f"{label}[{index}]")
        for index, value in enumerate(values)
    ]


def _field_rows(raw: Any, field: str) -> tuple[Any, Any, Any]:
    inherited = raw
    if isinstance(raw, Mapping):
        nested = _first(raw, "ordinary", "ordinary_Z", "ordinary_z", "rows", "jet")
        if nested is not _MISSING:
            inherited = raw
            raw = nested
    if isinstance(raw, Mapping):
        z0 = _first(raw, "Z0", "z0", "0", 0, "finite_Z0", "finite_z0")
        z1 = _first(raw, "Z1", "z1", "1", 1, "finite_Z1", "finite_z1")
    elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)) and len(raw) == 2:
        z0, z1 = raw[0], raw[1]
    else:
        z0 = z1 = _MISSING
    if z0 is _MISSING or z1 is _MISSING:
        raise DensityInputError(f"{field} requires explicit ordinary Z0 and Z1 rows")
    return z0, z1, inherited


def _normalise_field(raw: Any, field: str, backend: _Backend) -> dict[int, list[_Component]]:
    z0, z1, inherited = _field_rows(raw, field)
    rows: dict[int, list[_Component]] = {}
    for order, row in ((0, z0), (1, z1)):
        row_inherited = inherited
        if isinstance(row, Mapping):
            finite = _first(row, "finite", "finite_terms", "finite_interval")
            tails = _first(row, "tails", "tail", "additive_tails", "additive_tail", "tail_terms", "error_tails")
            if finite is _MISSING:
                interval_keys = ("lower", "upper", "lo", "hi", "left", "right", "interval", "bounds", "coefficient", "value")
                finite = row if any(_has(row, key) for key in interval_keys) else _MISSING
            if finite is _MISSING and _has(row, "terms", "components"):
                finite = _first(row, "terms", "components")
            if finite is _MISSING:
                raise DensityInputError(f"{field}.Z{order} requires a finite interval")
        else:
            finite, tails = row, _MISSING
        if tails is _MISSING:
            raise DensityInputError(f"{field}.Z{order} requires explicit additive tails (possibly an empty list)")
        finite_components = _as_component_rows(
            finite,
            kind="finite",
            field=field,
            order=order,
            inherited=row_inherited,
            backend=backend,
            label=f"{field}.Z{order}.finite",
        )
        tail_components = _as_component_rows(
            tails,
            kind="tail",
            field=field,
            order=order,
            inherited=row_inherited,
            backend=backend,
            label=f"{field}.Z{order}.tail",
            allow_empty=True,
        ) if tails is not _MISSING else []
        rows[order] = finite_components + tail_components
    return rows


def _lookup_field(fields: Mapping[Any, Any], field: str) -> Any:
    aliases = {
        "U0": ("U0", "u0", "U", "u"),
        "Delta": ("Delta", "delta", "DELTA", "Uz_minus_U0", "Uz-U0"),
        "Phi": ("Phi", "phi", "PHI"),
        "F1": ("F1", "f1", "F_1"),
        "Uz1": ("Uz1", "uz1", "UZ1", "U_z1", "u_z1"),
        "L": ("L", "l", "log_derivative", "logDerivative"),
    }
    value = _first(fields, *aliases[field])
    if value is _MISSING:
        raise DensityInputError(f"fields is missing {field}")
    return value


def _zero_derivative(inputs: Mapping[Any, Any], name: str) -> dict[str, Any]:
    aliases = {
        "rho_Z": ("rho_Z", "rho_z", "rho_derivative", "rho_Z_zero"),
        "epsilon_Z": ("epsilon_Z", "epsilon_z", "epsilon_derivative", "epsilon_Z_zero"),
        "eta_Z": ("eta_Z", "eta_z", "eta_derivative", "eta_Z_zero"),
    }[name]
    raw = _first(inputs, *aliases)
    if raw is _MISSING:
        geometry = _first(inputs, "geometry", "geometry_certificate", "premises")
        if isinstance(geometry, Mapping):
            raw = _first(geometry, *aliases)
    if raw is _MISSING:
        raise DensityInputError(f"explicit {name}=0 premise is required")
    if isinstance(raw, bool):
        if raw is not True:
            raise DensityInputError(f"{name} must be explicitly zero")
        return {"value": 0, "verified": True, "source": "caller_zero_flag"}
    if isinstance(raw, Mapping):
        if _has(raw, "zero", "is_zero") and _first(raw, "zero", "is_zero") is True:
            return {"value": 0, "verified": True, "source": "caller_zero_flag"}
        nested = _first(raw, "interval", "bounds", "range", "value", "derivative")
        if nested is not _MISSING:
            raw = nested
    pair = _extract_endpoints(raw, FractionIntervalBackend())
    if pair is None:
        value = _exact(raw, name)
        if value != 0:
            raise DensityInputError(f"{name} must be zero")
        return {"value": 0, "verified": True, "source": "caller_exact_zero"}
    lo, hi = _exact(pair[0], f"{name}.lower"), _exact(pair[1], f"{name}.upper")
    if lo != 0 or hi != 0:
        raise DensityInputError(f"{name} must be the exact zero interval")
    return {"value": {"lower": lo, "upper": hi}, "verified": True, "source": "caller_zero_interval"}


def _rho(inputs: Mapping[Any, Any], backend: _Backend) -> tuple[Any, dict[str, Any]]:
    raw = _first(inputs, "rho", "rho_interval", "fixed_rho")
    if raw is _MISSING:
        cell = _first(inputs, "cell")
        if isinstance(cell, Mapping):
            raw = _first(cell, "rho", "rho_interval", "fixed_rho")
    if raw is _MISSING:
        raise DensityInputError("fixed scalar rho interval is required")
    value = backend.interval(raw, "rho")
    lower, _ = backend.endpoints(value, "rho")
    if lower < 0:
        raise DensityInputError("rho must be nonnegative")
    return value, _provenance(raw)


def _field_packet(inputs: Mapping[Any, Any], backend: _Backend) -> tuple[dict[str, dict[int, list[_Component]]], Any, dict[str, Any]]:
    fields = _first(inputs, "fields", "ordinary_fields", "jet_fields")
    if fields is _MISSING:
        fields = inputs
    if not isinstance(fields, Mapping):
        raise DensityInputError("fields must be a mapping")
    result = {field: _normalise_field(_lookup_field(fields, field), field, backend) for field in FIELDS}
    rho, rho_info = _rho(inputs, backend)
    return result, rho, rho_info


def _merge_metadata(components: Sequence[_Component]) -> dict[str, Any]:
    source_ids: list[Any] = []
    term_ids: list[Any] = []
    source_paths: list[Any] = []
    error_classes: list[Any] = []
    scales: list[Any] = []
    uncertain = False
    for component in components:
        source_ids.extend(component.source_ids)
        term_ids.extend(component.term_ids)
        source_paths.extend(component.source_paths)
        error_classes.extend(component.error_classes)
        scales.extend(component.scale_metadata)
        uncertain = uncertain or component.uncertainty
    return {
        "source_ids": _unique(source_ids),
        "term_ids": _unique(term_ids),
        "source_paths": _unique(source_paths),
        "error_classes": _unique(error_classes),
        "uncertainty": uncertain,
        "scale_metadata": scales,
    }


def _term(
    backend: _Backend,
    components: Sequence[_Component],
    factor: int,
    *,
    channel: str,
    order: int,
    piece: str,
    index: int,
    powers: tuple[int, int],
) -> dict[str, Any]:
    coefficient = components[0].interval
    for component in components[1:]:
        coefficient = backend.mul(coefficient, component.interval, f"{channel}.Z{order}.{piece}")
    coefficient = backend.scale(factor, coefficient, f"{channel}.Z{order}.{piece}")
    info = _merge_metadata(components)
    tail_components = [
        {"field": component.field, "label": component.label, "tail_id": component.tail_id}
        for component in components
        if component.kind == "tail"
    ]
    component_labels = [f"{component.field}:{component.label}" for component in components]
    finite_tail_signature = {component.field: component.kind for component in components}
    original_scale = {
        component.field: list(component.scale_metadata)
        for component in components
        if component.scale_metadata
    }
    ordinary_label = ORDINARY_Z_LABELS[order]
    sign = 1 if factor >= 0 else -1
    term_id = f"{channel}.{ordinary_label}.{piece}.{index}"
    origin = {
        "operation": "scaled_m15_density_distributive_product",
        "term_id": term_id,
        "parents": list(info["term_ids"]),
        "source_ids": list(info["source_ids"]),
        "error_classes": list(info["error_classes"]),
        "uncertainty": bool(info["uncertainty"]),
        "symbolic_cancellation_claimed": False,
        "aggregation_performed": False,
    }
    serial = backend.serial(coefficient)
    # The formula contributes the outer powers; caller-factored input powers
    # must also be composed, not merely copied into provenance metadata.
    outer_powers = powers
    powers = tuple(outer_powers[i] + sum((component.formal_powers or (0, 0))[i]
                                        for component in components) for i in range(2))
    return {
        "term_id": term_id,
        "channel": channel,
        "ordinary_Z_order": order,
        "ordinary_Z": ordinary_label,
        "piece": piece,
        "coefficient": serial,
        "coefficient_interval": serial,
        "signed_factor": factor,
        "sign": sign,
        "formal_powers": {"epsilon": powers[0], "F0": powers[1]},
        "outer_formal_powers": {"epsilon": outer_powers[0], "F0": outer_powers[1]},
        "powers": [powers[0], powers[1]],
        "epsilon_power": powers[0],
        "F0_power": powers[1],
        "component_labels": component_labels,
        "finite_tail_signature": finite_tail_signature,
        "tail_components": tail_components,
        "has_additive_tail": bool(tail_components),
        "source_ids": list(info["source_ids"]),
        "input_source_ids": list(info["source_ids"]),
        "term_ids": list(info["term_ids"]),
        "input_term_ids": list(info["term_ids"]),
        "source_paths": list(info["source_paths"]),
        "error_classes": list(info["error_classes"]),
        "uncertainty": bool(info["uncertainty"]),
        "error_only": bool(tail_components),
        "original_scale_metadata": original_scale,
        "scale_metadata": original_scale,
        "origin": origin,
    }


def _make_terms(backend: _Backend, fields: Mapping[str, dict[int, list[_Component]]], rho: Any, rho_info: Mapping[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    axial: list[dict[str, Any]] = []
    swirl: list[dict[str, Any]] = []
    counter = 0
    # Z0 axial: 2*U0*Uz1 + 2*Delta*Uz1.
    for left_name in ("U0", "Delta"):
        for left in fields[left_name][0]:
            for right in fields["Uz1"][0]:
                axial.append(_term(backend, (left, right), 2, channel="axial", order=0, piece=f"{left_name}_0*Uz1_0", index=counter, powers=(1, 0)))
                counter += 1
    counter = 0
    # Z1 axial product rule, epsilon_Z=0.
    for left_name in ("U0", "Delta"):
        for left in fields[left_name][1]:
            for right in fields["Uz1"][0]:
                axial.append(_term(backend, (left, right), 2, channel="axial", order=1, piece=f"{left_name}_Z*Uz1_0", index=counter, powers=(1, 0)))
                counter += 1
    for left_name in ("U0", "Delta"):
        for left in fields[left_name][0]:
            for right in fields["Uz1"][1]:
                axial.append(_term(backend, (left, right), 2, channel="axial", order=1, piece=f"{left_name}_0*Uz1_Z", index=counter, powers=(1, 0)))
                counter += 1

    # rho is a fixed scalar factor; it is intentionally not differentiated.
    rho_component = _Component(
        interval=rho,
        field="rho",
        ordinary_order=0,
        kind="finite",
        label="finite",
        tail_id=None,
        source_ids=tuple(rho_info["source_ids"]),
        term_ids=tuple(rho_info["term_ids"]),
        source_paths=tuple(rho_info["source_paths"]),
        error_classes=tuple(rho_info["error_classes"]),
        uncertainty=rho_info["uncertainty"],
        scale_metadata=tuple(rho_info["scale"]),
        formal_powers=(0, 0),
        raw_metadata=None,
    )
    counter = 0
    for phi in fields["Phi"][0]:
        for f1 in fields["F1"][0]:
            swirl.append(_term(backend, (rho_component, phi, f1), -2, channel="swirl", order=0, piece="rho*Phi_0*F1_0", index=counter, powers=(2, 1)))
            counter += 1
    counter = 0
    # Z1 swirl: -2*rho*(L*Phi*F1 + Phi_Z*F1 + Phi*F1_Z).
    for log_derivative in fields["L"][0]:
        for phi in fields["Phi"][0]:
            for f1 in fields["F1"][0]:
                swirl.append(_term(backend, (rho_component, log_derivative, phi, f1), -2, channel="swirl", order=1, piece="rho*L_0*Phi_0*F1_0", index=counter, powers=(2, 1)))
                counter += 1
    for phi in fields["Phi"][1]:
        for f1 in fields["F1"][0]:
            swirl.append(_term(backend, (rho_component, phi, f1), -2, channel="swirl", order=1, piece="rho*Phi_Z*F1_0", index=counter, powers=(2, 1)))
            counter += 1
    for phi in fields["Phi"][0]:
        for f1 in fields["F1"][1]:
            swirl.append(_term(backend, (rho_component, phi, f1), -2, channel="swirl", order=1, piece="rho*Phi_0*F1_Z", index=counter, powers=(2, 1)))
            counter += 1
    return axial, swirl


def _pressure_input(inputs: Mapping[Any, Any], backend: _Backend) -> dict[int, list[_Component]] | None:
    """Normalize the optional *original* ``N1p`` row without merging it."""

    include = _first(inputs, "include_pressure")
    raw = _first(inputs, "N1p", "n1p", "pressure_N1p", "pressure")
    if raw is not _MISSING and isinstance(raw, Mapping):
        nested = _first(raw, "N1p", "n1p", "original_N1p")
        if nested is not _MISSING:
            raw = nested
    if raw is not _MISSING:
        if include is _MISSING:
            raise DensityInputError("N1p requires an explicit include_pressure flag")
        if include is not True:
            raise DensityInputError("N1p is present while include_pressure is not true")
    if include is _MISSING or include is False:
        return None
    if include is not True:
        raise DensityInputError("include_pressure must be the exact boolean True or False")
    if raw is _MISSING:
        raise DensityInputError("include_pressure=True requires original N1p Z0/Z1 rows")
    if isinstance(raw, Mapping):
        normalized = _first(raw, "normalized", "normalized_N1p", "N1p_normalized", "normalized_pressure")
        if normalized not in (_MISSING, False, None):
            raise DensityInputError("normalized N1p cannot be supplied to the original pressure channel")
    pressure = _normalise_field(raw, "N1p", backend)
    for order in ORDINARY_Z_ORDERS:
        for component in pressure[order]:
            if component.formal_powers is None:
                # The pressure channel remains disconnected when its original
                # scale ledger is absent; no guessed epsilon/F0 power is
                # allowed to enter the velocity-density result.
                raise DensityInputError("original N1p requires explicit formal powers")
    return pressure


def _make_pressure_terms(
    backend: _Backend,
    pressure: dict[int, list[_Component]],
    rho: Any,
    rho_info: Mapping[str, Any],
) -> list[dict[str, Any]]:
    rho_component = _Component(
        interval=rho,
        field="rho",
        ordinary_order=0,
        kind="finite",
        label="finite",
        tail_id=None,
        source_ids=tuple(rho_info["source_ids"]),
        term_ids=tuple(rho_info["term_ids"]),
        source_paths=tuple(rho_info["source_paths"]),
        error_classes=tuple(rho_info["error_classes"]),
        uncertainty=rho_info["uncertainty"],
        scale_metadata=tuple(rho_info["scale"]),
        formal_powers=(0, 0),
        raw_metadata=None,
    )
    terms: list[dict[str, Any]] = []
    counter = 0
    for order in ORDINARY_Z_ORDERS:
        for component in pressure[order]:
            assert component.formal_powers is not None
            term = _term(
                backend,
                (rho_component, component),
                -1,
                channel="pressure_N1p",
                order=order,
                piece="rho*N1p" if order == 0 else "rho*N1p_Z",
                index=counter,
                powers=(2, 0),
            )
            term["pressure_input_formal_powers"] = list(component.formal_powers)
            term["pressure_channel"] = "original N1p only"
            terms.append(term)
            counter += 1
    return terms


def _certificate(inputs: Mapping[Any, Any], *, cell_count: int) -> dict[str, Any]:
    scope = "one fixed-rho density cell; no integration or whole-cell enclosure"
    common = {
        "verified": False,
        "unverified": True,
        "scope": scope,
        "reason": "density algebra has no source owner, binder, or whole-cell proof",
    }
    return {
        "whole_cell": dict(common, claim="whole-cell coverage not established", cell_count=cell_count),
        "owner": dict(common, claim="source owner not established"),
        "binder": dict(common, claim="source/binder identity not established"),
        "whole_cell_verified": False,
        "owner_verified": False,
        "binder_verified": False,
        "input_certificate_received": _first(inputs, "certificate", "geometry_certificate", "binder_certificate") is not _MISSING,
        "unverified": True,
    }


def _cell_input(inputs: Mapping[Any, Any]) -> Mapping[Any, Any]:
    cell = _first(inputs, "cell")
    if isinstance(cell, Mapping):
        merged = dict(inputs)
        merged.pop("cell", None)
        merged.update(cell)
        return merged
    return inputs


def _eta_status(inputs: Mapping[Any, Any]) -> dict[str, Any]:
    """Record whether saved ``F1``/``Uz1`` rows already include ``eta``."""

    raw = _first(inputs, "responses_include_eta", "eta_included", "F1_Uz1_include_eta")
    if raw is _MISSING:
        raise DensityInputError("responses_include_eta=True is required for the cutoff density")
    if raw is not True and raw is not False:
        raise DensityInputError("responses_include_eta must be boolean")
    if raw is True:
        eta_zero = _zero_derivative(inputs, "eta_Z")
        return {
            "included": True,
            "verified": True,
            "status": "F1 and Uz1 are caller-certified eta-included responses",
            "eta_Z": eta_zero,
        }
    raise DensityInputError("eta-free responses cannot be used for the cutoff density")


def _one_cell(inputs: Mapping[Any, Any], backend: _Backend, cell_index: int) -> dict[str, Any]:
    inputs = _cell_input(inputs)
    _reject_claims(inputs)
    for forbidden in ("epsilon", "epsilon_value", "epsilon_scale", "F0", "F0_value", "F0_scale", "materialized_F0"):
        if _first(inputs, forbidden) is not _MISSING:
            raise DensityInputError(f"{forbidden} must not be supplied or materialized")
    fields, rho, rho_info = _field_packet(inputs, backend)
    zero_rho = _zero_derivative(inputs, "rho_Z")
    zero_epsilon = _zero_derivative(inputs, "epsilon_Z")
    eta_status = _eta_status(inputs)
    axial, swirl = _make_terms(backend, fields, rho, rho_info)
    pressure = _pressure_input(inputs, backend)
    pressure_terms = _make_pressure_terms(backend, pressure, rho, rho_info) if pressure is not None else []
    cell_id = _first(inputs, "cell_id", "id")
    if cell_id is _MISSING:
        cell_id = cell_index
    for term in axial + swirl + pressure_terms:
        term["cell_index"] = cell_index
        term["cell_id"] = cell_id
        term["term_id"] = f"cell{cell_index}:" + term["term_id"]
        term["origin"]["term_id"] = term["term_id"]
    axial_by_order = {str(order): [term for term in axial if term["ordinary_Z_order"] == order] for order in ORDINARY_Z_ORDERS}
    swirl_by_order = {str(order): [term for term in swirl if term["ordinary_Z_order"] == order] for order in ORDINARY_Z_ORDERS}
    pressure_by_order = {str(order): [term for term in pressure_terms if term["ordinary_Z_order"] == order] for order in ORDINARY_Z_ORDERS}
    return {
        "cell_index": cell_index,
        "cell_id": cell_id,
        "rho": backend.serial(rho),
        "rho_source_ids": list(rho_info["source_ids"]),
        "axial": axial_by_order,
        "swirl": swirl_by_order,
        "pressure": {"included": pressure is not None, "Z0": pressure_by_order["0"], "Z1": pressure_by_order["1"]},
        "terms": axial + swirl + pressure_terms,
        "zero_derivative_premises": {"rho_Z": zero_rho, "epsilon_Z": zero_epsilon, "eta_Z": eta_status["eta_Z"]},
        "eta_status": eta_status,
        "fields_present": list(FIELDS),
    }


def density(inputs: Any = None, backend: Any | None = None, **kwargs: Any) -> dict[str, Any]:
    """Build signed per-cell axial and swirl density terms.

    Canonical use is ``density(packet, backend=backend)``.  ``packet`` has a
    scalar ``rho`` interval, explicit zero premises ``rho_Z`` and
    ``epsilon_Z``, and ``fields`` containing the six fields with ``Z0`` and
    ``Z1`` rows.  A row is ``{"finite": interval, "tails": [interval, ...]}``;
    tails are additive and remain distinct in the returned distributive terms.
    ``cells`` may be supplied for a list of independent fixed-rho cells.
    """

    if inputs is None and "packet" in kwargs:
        inputs = kwargs.pop("packet")
    if backend is None:
        backend = kwargs.pop("interval_backend", kwargs.pop("exact_backend", None))
    if kwargs:
        names = ", ".join(sorted(str(key) for key in kwargs))
        raise TypeError(f"unexpected density arguments: {names}")
    # Also accept density(backend, packet) for small adapter callers.
    if not isinstance(inputs, Mapping) and isinstance(backend, Mapping):
        inputs, backend = backend, inputs
    if not isinstance(inputs, Mapping):
        raise DensityInputError("inputs must be a mapping")
    range_backend = _Backend(backend)
    cells_raw = _first(inputs, "cells", "target_cells")
    if cells_raw is _MISSING:
        cell_inputs = [inputs]
    else:
        if not isinstance(cells_raw, Sequence) or isinstance(cells_raw, (str, bytes, bytearray)) or not cells_raw:
            raise DensityInputError("cells must be a nonempty sequence")
        cell_inputs = []
        shared = dict(inputs)
        shared.pop("cells", None)
        shared.pop("target_cells", None)
        for index, raw_cell in enumerate(cells_raw):
            if not isinstance(raw_cell, Mapping):
                raise DensityInputError(f"cells[{index}] must be a mapping")
            merged = dict(shared)
            merged.update(raw_cell)
            cell_inputs.append(merged)
    cells = [_one_cell(cell, range_backend, index) for index, cell in enumerate(cell_inputs)]
    certificate = _certificate(inputs, cell_count=len(cells))
    all_axial = [term for cell in cells for term in cell["axial"]["0"] + cell["axial"]["1"]]
    all_swirl = [term for cell in cells for term in cell["swirl"]["0"] + cell["swirl"]["1"]]
    all_pressure = [term for cell in cells for term in cell["pressure"]["Z0"] + cell["pressure"]["Z1"]]
    all_terms = all_axial + all_swirl + all_pressure
    pressure_included = bool(all_pressure)
    eta_statuses = [cell["eta_status"] for cell in cells]
    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "ordinary_Z_orders": [0, 1],
        "ordinary_Z_labels": list(ORDINARY_Z_LABELS),
        "formal_power_basis": list(FORMAL_POWER_BASIS),
        "formal_powers_are_separate": True,
        "cells": cells,
        "cell_count": len(cells),
        "axial": {"Z0": [term for term in all_axial if term["ordinary_Z_order"] == 0], "Z1": [term for term in all_axial if term["ordinary_Z_order"] == 1]},
        "swirl": {"Z0": [term for term in all_swirl if term["ordinary_Z_order"] == 0], "Z1": [term for term in all_swirl if term["ordinary_Z_order"] == 1]},
        "pressure": {"included": pressure_included, "Z0": [term for term in all_pressure if term["ordinary_Z_order"] == 0], "Z1": [term for term in all_pressure if term["ordinary_Z_order"] == 1]},
        "terms": all_terms,
        "term_count": len(all_terms),
        "finite_term_count": sum(1 for term in all_terms if not term["has_additive_tail"]),
        "tail_term_count": sum(1 for term in all_terms if term["has_additive_tail"]),
        "formulas": {
            "axial_Z0": "epsilon * 2*(U0+Delta)*Uz1",
            "swirl_Z0": "epsilon^2*F0 * (-2*rho*Phi*F1)",
            "axial_Z1": "epsilon * 2*((U0_Z+Delta_Z)*Uz1 + (U0+Delta)*Uz1_Z)",
            "swirl_Z1": "epsilon^2*F0 * (-2*rho*(L*Phi*F1 + Phi_Z*F1 + Phi*F1_Z))",
            "pressure_N1p_Z0": "-epsilon^2*rho*N1p (optional original channel; not merged)",
            "pressure_N1p_Z1": "-epsilon^2*rho*N1p_Z (optional original channel; not merged)",
            "rho_Z": "0",
            "epsilon_Z": "0",
        },
        "rho_is_fixed": True,
        "rho_Z_is_zero": True,
        "epsilon_Z_is_zero": True,
        "eta_status": eta_statuses[0] if len(eta_statuses) == 1 else eta_statuses,
        "geometry_and_eta_premises_are_caller_claims": True,
        "input_Z1_is_full_ordinary_derivative_before_factoring": True,
        "eta_Z_is_zero_when_eta_included": all(status["eta_Z"].get("verified") is True for status in eta_statuses if status["included"]),
        "responses_include_eta": all(status["included"] for status in eta_statuses),
        "tails_distributed": True,
        "tails_remain_separate": True,
        "aggregation_performed": False,
        "error_tails_collapsed": False,
        "old_errors_erased": False,
        "epsilon_materialized": False,
        "F0_materialized": False,
        "epsilon_divided": False,
        "integration_performed": False,
        "source_imports_used": False,
        "actual_source_imported": False,
        "pressure_channel_included": pressure_included,
        "pressure_channel_is_separate": True,
        "pressure_channel_is_original_N1p": pressure_included,
        "N1p_normalized": False,
        "W1_axial_swirl_only_by_default": not pressure_included,
        "N1p_pressure_not_merged_into_velocity": True,
        "pressure_endpoint_term_not_summed": True,
        "full_m15_claimed": False,
        "complete_m15_computed": False,
        "endpoint_sum_performed": False,
        "certificate": certificate,
        "whole_cell_certificate": certificate["whole_cell"],
        "owner_certificate": certificate["owner"],
        "binder_certificate": certificate["binder"],
        "scope": "signed coefficient-space per-cell W1 m15 density only",
    }


build = density
evaluate = density
build_density = density


def serialize(value: Any) -> Any:
    """Convert result-like values to JSON-shaped containers, retaining Fraction."""

    if isinstance(value, RationalInterval):
        return {"lower": value.lower, "upper": value.upper}
    if isinstance(value, Mapping):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    return value


__all__ = [
    "SCHEMA",
    "COMPONENT",
    "ORDINARY_Z_ORDERS",
    "ORDINARY_Z_LABELS",
    "FIELDS",
    "FORMAL_POWER_BASIS",
    "DensityInputError",
    "ExactRangeBackendError",
    "RationalInterval",
    "FractionIntervalBackend",
    "density",
    "build",
    "evaluate",
    "build_density",
    "serialize",
]
