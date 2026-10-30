"""Compact source-factored first-collar Laurent forcing leaf.

This leaf is deliberately self contained.  It accepts the actual radial
receipt (or ordinary derivative grids) and never imports a producer, executes
an older implementation, or reads a saved JSON packet.  The physical
equations and the ordinary ``Z0`` through ``Z3`` interface are the first
collar equations used by the Laurent forcing leaf.

The difference from the unbounded Laurent representation is important:
every Taylor coefficient is reduced immediately by the key
``(four formal powers, uncertainty)``.  Reduction is directed interval
addition only; it is not a functional identity, a source cancellation, or a
symmetry argument.  The reduced origin retains the union of parent
``source_ids`` and an exact aggregate ``path_count``.  A parent marked
``uncertainty=True`` can therefore never be added to the ordinary source
class.  Origins are small references rather than recursive copies of input
records.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
_ZERO_POWERS = (0, 0, 0, 0)
_MAX_SOURCE_ORDER = 5
_OUTPUT_ORDERS = (0, 1, 2, 3)
_RHO_MIN = "4"
_RHO_MAX = "4.1"
_Z_MIN = "-1"
_Z_MAX = "1"

__all__ = [
    "SOURCE_POWER_BASIS",
    "IntervalTaylorJet",
    "TaylorJet",
    "LaurentTaylorJet",
    "FormalLaurentJet",
    "build_first_collar_laurent_forcing_jets",
    "first_collar_laurent_forcing_jets",
    "derive_first_collar_laurent_forcing_jets",
    "build_first_collar_laurent_forcing",
    "first_collar_laurent_forcing",
    "build_collar_laurent_forcing",
    "build_first_collar_laurent_forcing_z3",
    "first_collar_laurent_forcing_z3",
    "build_first_collar_forcing_laurent_jets",
    "first_collar_forcing_laurent_jets",
    "build_first_collar_forcing_laurent",
    "first_collar_forcing_laurent",
]


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") or (hasattr(value, "a") and hasattr(value, "b"))


def _bounds(value: Any) -> tuple[Any, Any]:
    if not _is_interval(value):
        raise TypeError("an mpmath interval is required")
    try:
        return value.a, value.b
    except AttributeError:  # pragma: no cover - defensive for iv variants
        rows = value._mpi_
        return rows[0], rows[1]


def _interval(context: Any, value: Any, label: str = "interval") -> Any:
    """Coerce without selecting a midpoint or an endpoint."""

    if _is_interval(value):
        lo, hi = _bounds(value)
        return context.mpf([lo, hi])
    if isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            return context.make_mpf(
                (tuple(value["lower_exact_mpf_tuple"]), tuple(value["upper_exact_mpf_tuple"]))
            )
        if "lower" in value and "upper" in value:
            return context.mpf([value["lower"], value["upper"]])
        # Receipts sometimes wrap an interval in one named field.
        for name in ("interval", "enclosure", "value", "coefficient", "coefficient_interval"):
            if name in value:
                return _interval(context, value[name], label + "." + name)
        raise TypeError(f"{label} needs lower/upper endpoints")
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)) and len(value) == 2:
        return context.mpf([value[0], value[1]])
    try:
        return context.mpf(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"unsupported {label}") from exc


def _require_nonzero(value: Any, label: str) -> None:
    lo, hi = _bounds(value)
    if lo <= 0 <= hi:
        raise ZeroDivisionError(f"{label} interval contains zero")


def _require_positive(value: Any, label: str) -> None:
    lo, _ = _bounds(value)
    if lo <= 0:
        raise ValueError(f"{label} must be strictly positive")


def _require_domain(value: Any, lower: Any, upper: Any, context: Any, label: str) -> None:
    lo, hi = _bounds(value)
    if lo < context.mpf(lower) or hi > context.mpf(upper):
        raise ValueError(f"{label} is outside [{lower}, {upper}]")


def _factorial_div(value: Any, order: int, context: Any) -> Any:
    return value / context.mpf(math.factorial(order))


def _as_ids(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if isinstance(value, Mapping):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    return [value]


def _as_labels(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    return [value]


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    seen = set()
    for value in values:
        try:
            key = (type(value).__name__, value)
            hash(key)
        except TypeError:
            key = (type(value).__name__, repr(value))
        if key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _origin_info(origin: Any) -> dict[str, Any]:
    """Read only the provenance fields needed for compact aggregation."""

    if isinstance(origin, Mapping):
        source_ids = list(_as_ids(origin.get("source_ids", origin.get("sources"))))
        error_classes = _unique(_as_labels(origin.get("error_classes", origin.get("uncertainty_class"))))
        operation = str(origin.get("operation", "")).lower()
        error_kind = str(origin.get("error_kind", origin.get("kind", ""))).lower()
        uncertain = bool(
            origin.get("uncertainty", origin.get("is_uncertainty", False))
            or "cauchy" in operation
            or "symmetric" in operation
            or "cauchy" in error_kind
            or "symmetric" in error_kind
        )
        try:
            path_count = int(origin.get("path_count", origin.get("aggregate_path_count", 1)))
        except (TypeError, ValueError):
            path_count = 1
        if path_count < 1:
            path_count = 1
        return {
            "source_ids": _unique(source_ids),
            "uncertainty": uncertain,
            "error_classes": error_classes,
            "path_count": path_count,
            "operation": origin.get("operation"),
        }
    # A prior radial receipt may have wrapped its dictionary origin in a
    # small tuple.  Inspect nested provenance records only; scalar tuple
    # labels are not source ids and are intentionally ignored.
    if isinstance(origin, Sequence) and not isinstance(origin, (str, bytes, bytearray)):
        nested = [_origin_info(item) for item in origin if isinstance(item, (Mapping, Sequence)) and not isinstance(item, (str, bytes, bytearray))]
        source_ids: list[Any] = []
        error_classes: list[Any] = []
        uncertain = False
        path_count = 1
        for info in nested:
            source_ids.extend(info["source_ids"])
            error_classes.extend(info.get("error_classes", []))
            uncertain = uncertain or bool(info["uncertainty"])
            path_count = max(path_count, int(info["path_count"]))
        return {"source_ids": _unique(source_ids), "uncertainty": uncertain, "error_classes": _unique(error_classes), "path_count": path_count, "operation": None}
    return {"source_ids": [], "uncertainty": False, "error_classes": [], "path_count": 1, "operation": None}


def _origin_ref(origin: Any) -> dict[str, Any]:
    info = _origin_info(origin)
    result: dict[str, Any] = {
        "source_ids": list(info["source_ids"]),
        "uncertainty": bool(info["uncertainty"]),
        "error_classes": list(info.get("error_classes", [])),
        "path_count": int(info["path_count"]),
    }
    if info.get("operation") is not None:
        result["operation"] = info["operation"]
    return result


def _origin(operation: str, *inputs: Any, **metadata: Any) -> dict[str, Any]:
    """Create a bounded provenance reference, never a recursive input copy."""

    infos = [_origin_info(item) for item in inputs]
    ids: list[Any] = []
    error_classes: list[Any] = []
    uncertain = bool(metadata.get("uncertainty", False)) or bool(
        metadata.get("symmetric_error", False) or metadata.get("cauchy_error", False)
    )
    if operation in {"multiply", "product"}:
        path_count = 1
        for info in infos:
            path_count *= int(info["path_count"])
    elif operation in {"aggregate", "add", "sum"}:
        path_count = 0
        for info in infos:
            path_count += int(info["path_count"])
    else:
        path_count = int(infos[0]["path_count"]) if infos else 1
    for info in infos:
        ids.extend(info["source_ids"])
        error_classes.extend(info.get("error_classes", []))
        uncertain = uncertain or bool(info["uncertainty"])
    if "source_ids" in metadata:
        ids.extend(_as_ids(metadata["source_ids"]))
    if "error_classes" in metadata:
        error_classes.extend(_as_labels(metadata["error_classes"]))
    result: dict[str, Any] = {
        "operation": operation,
        "source_ids": _unique(ids),
        "uncertainty": bool(uncertain),
        "error_classes": _unique(error_classes),
        "path_count": max(1, int(path_count)),
        "aggregate_path_count": max(1, int(path_count)),
    }
    # Parent references are intentionally shallow: they preserve the source
    # identity without retaining an entire product tree or input record.
    if inputs:
        result["parents"] = tuple(_origin_ref(item) for item in inputs)
    for key, value in metadata.items():
        if key not in {"source_ids", "uncertainty", "error_classes", "symmetric_error", "cauchy_error"}:
            result[key] = value
    return result


class IntervalTaylorJet:
    """Factorial-normalised interval Taylor algebra for numeric geometry."""

    __slots__ = ("context", "coefficients")

    def __init__(self, context: Any, coefficients: Sequence[Any]):
        values = tuple(coefficients)
        if not values:
            raise ValueError("a Taylor jet needs at least its constant coefficient")
        self.context = context
        self.coefficients = values

    @property
    def order(self) -> int:
        return len(self.coefficients) - 1

    @classmethod
    def constant(cls, context: Any, value: Any, order: int) -> "IntervalTaylorJet":
        if order < 0:
            raise ValueError("nonnegative Taylor order required")
        value = _interval(context, value, "constant")
        return cls(context, [value] + [context.mpf(0)] * order)

    @classmethod
    def variable(cls, context: Any, value: Any, order: int) -> "IntervalTaylorJet":
        if order < 1:
            return cls.constant(context, value, order)
        value = _interval(context, value, "Taylor variable")
        return cls(context, [value, context.mpf(1)] + [context.mpf(0)] * (order - 1))

    def truncate(self, order: int) -> "IntervalTaylorJet":
        if order < 0 or order > self.order:
            raise ValueError("requested Taylor order is unavailable")
        return type(self)(self.context, self.coefficients[: order + 1])

    def differentiate(self) -> "IntervalTaylorJet":
        if self.order == 0:
            raise ValueError("cannot differentiate a jet with no known derivative")
        return type(self)(self.context, [(index + 1) * self[index + 1] for index in range(self.order)])

    dz = differentiate

    def __getitem__(self, index: int) -> Any:
        return self.coefficients[index]

    def _coerce(self, other: Any) -> "IntervalTaylorJet":
        if isinstance(other, IntervalTaylorJet):
            if other.context is not self.context:
                raise TypeError("Taylor jets must use the caller's one interval context")
            return other
        return type(self).constant(self.context, other, self.order)

    def __add__(self, other: Any) -> "IntervalTaylorJet":
        other = self._coerce(other)
        order = min(self.order, other.order)
        return type(self)(self.context, [self[index] + other[index] for index in range(order + 1)])

    __radd__ = __add__

    def __neg__(self) -> "IntervalTaylorJet":
        return type(self)(self.context, [-value for value in self.coefficients])

    def __sub__(self, other: Any) -> "IntervalTaylorJet":
        return self + (-self._coerce(other))

    def __rsub__(self, other: Any) -> "IntervalTaylorJet":
        return self._coerce(other) - self

    def __mul__(self, other: Any) -> "IntervalTaylorJet":
        other = self._coerce(other)
        order = min(self.order, other.order)
        coefficients = []
        for index in range(order + 1):
            value = self.context.mpf(0)
            for left in range(index + 1):
                value += self[left] * other[index - left]
            coefficients.append(value)
        return type(self)(self.context, coefficients)

    __rmul__ = __mul__

    def reciprocal(self) -> "IntervalTaylorJet":
        _require_nonzero(self[0], "Taylor denominator")
        coefficients = [self.context.mpf(1) / self[0]]
        for index in range(1, self.order + 1):
            value = self.context.mpf(0)
            for left in range(1, index + 1):
                value += self[left] * coefficients[index - left]
            coefficients.append(-value / self[0])
        return type(self)(self.context, coefficients)

    def __truediv__(self, other: Any) -> "IntervalTaylorJet":
        return self * self._coerce(other).reciprocal()

    def __rtruediv__(self, other: Any) -> "IntervalTaylorJet":
        return self._coerce(other) * self.reciprocal()

    def __pow__(self, power: int) -> "IntervalTaylorJet":
        if not isinstance(power, int) or power < 0:
            raise ValueError("only nonnegative integer Taylor powers are supported")
        result = type(self).constant(self.context, 1, self.order)
        for _ in range(power):
            result = result * self
        return result


TaylorJet = IntervalTaylorJet


def _mapping_get(mapping: Mapping[Any, Any], names: Sequence[Any]) -> Any | None:
    for name in names:
        try:
            if name in mapping:
                return mapping[name]
        except (TypeError, ValueError):
            continue
    return None


def _source_value(grid: Any, radial_order: int, axial_order: int, label: str) -> Any:
    """Read one required ordinary derivative; never synthesize a row."""

    if grid is None:
        raise TypeError(f"{label} source grid is required")
    # A live jet-like object is accepted, but only for its explicit radial-0
    # row.  The receipt itself is otherwise treated as a plain mapping.
    if hasattr(grid, "coefficients") and not isinstance(grid, Mapping):
        if radial_order != 0 or axial_order >= len(grid.coefficients):
            raise KeyError(f"missing {label} ordinary derivative rho^{radial_order} Z^{axial_order}")
        return grid.coefficients[axial_order]
    if isinstance(grid, Mapping):
        nested = _mapping_get(grid, ("grid", "derivatives", "radial_rows", "actual"))
        if nested is not None and nested is not grid:
            try:
                return _source_value(nested, radial_order, axial_order, label)
            except (KeyError, TypeError):
                pass
        direct = _mapping_get(
            grid,
            (
                (radial_order, axial_order),
                f"rho{radial_order}_Z{axial_order}",
                f"rho{radial_order}_z{axial_order}",
                f"rho_{radial_order}_Z_{axial_order}",
                f"rho{radial_order}Z{axial_order}",
                f"{radial_order},{axial_order}",
                f"{radial_order}:{axial_order}",
            ),
        )
        if direct is not None:
            return direct
        if radial_order == 0:
            direct_z = _mapping_get(grid, (axial_order, str(axial_order), f"Z{axial_order}", f"z{axial_order}"))
            if direct_z is not None:
                return direct_z
        row = _mapping_get(grid, (radial_order, str(radial_order), f"rho{radial_order}", f"rho_{radial_order}", f"radial_{radial_order}"))
        if isinstance(row, Mapping):
            nested_row = _mapping_get(row, (axial_order, str(axial_order), f"Z{axial_order}", f"z{axial_order}"))
            if nested_row is not None:
                return nested_row
        elif isinstance(row, Sequence) and not isinstance(row, (str, bytes, bytearray)) and axial_order < len(row):
            return row[axial_order]
    if isinstance(grid, Sequence) and not isinstance(grid, (str, bytes, bytearray)):
        if radial_order < len(grid):
            row = grid[radial_order]
            if isinstance(row, Sequence) and not isinstance(row, (str, bytes, bytearray)):
                if axial_order < len(row):
                    return row[axial_order]
        # one-dimensional Z row
        if radial_order == 0 and axial_order < len(grid):
            return grid[axial_order]
    raise KeyError(f"missing {label} ordinary derivative rho^{radial_order} Z^{axial_order}")


def _pick(mapping: Mapping[Any, Any], names: Sequence[str]) -> Any | None:
    return _mapping_get(mapping, names)


def _receipt_containers(*values: Any) -> list[Mapping[Any, Any]]:
    """Expose common actual-radial receipt wrappers without reading files."""

    result: list[Mapping[Any, Any]] = []
    pending = [value for value in values if isinstance(value, Mapping)]
    seen: set[int] = set()
    while pending:
        current = pending.pop(0)
        if id(current) in seen:
            continue
        seen.add(id(current))
        result.append(current)
        for name in ("radial_rows", "rows", "actual", "fields", "receipt", "source"):
            child = current.get(name)
            if isinstance(child, Mapping):
                pending.append(child)
    return result


def _resolve_inputs(
    *,
    Phi_derivatives: Any,
    Uz_derivatives: Any,
    radial_average_Uz_derivatives: Any,
    mean_derivatives: Any,
    log_F0: Any,
    F0_log: Any,
    logF0: Any,
    F0_log_interval: Any,
    containers: Sequence[Mapping[Any, Any]],
) -> tuple[Any, Any, Any, Any]:
    def first(value: Any, names: Sequence[str]) -> Any:
        if value is not None:
            return value
        for container in containers:
            found = _pick(container, names)
            if found is not None:
                return found
        return None

    phi = first(Phi_derivatives, ("Phi_derivatives", "actual_Phi", "Phi", "phi_derivatives", "phi"))
    uz = first(Uz_derivatives, ("Uz_derivatives", "actual_U", "U", "U_derivatives", "uz_derivatives", "u"))
    mean = first(
        radial_average_Uz_derivatives,
        ("radial_average_Uz_derivatives", "actual_mean", "mean_derivatives", "M_derivatives", "mean", "M"),
    )
    if mean is None:
        mean = first(mean_derivatives, ("mean_derivatives", "actual_mean", "M_derivatives", "mean", "M"))
    log_value = first(log_F0, ("log_F0", "F0_log", "logF0", "F0_log_interval"))
    if log_value is None:
        log_value = first(F0_log, ("log_F0", "F0_log", "logF0", "F0_log_interval"))
    if log_value is None:
        log_value = first(logF0, ("log_F0", "F0_log", "logF0", "F0_log_interval"))
    if log_value is None:
        log_value = first(F0_log_interval, ("log_F0", "F0_log", "logF0", "F0_log_interval"))
    if log_value is None:
        log_value = first(None, ("F0_exact_positive_log_enclosure", "F0_log_enclosure", "log_F0_interval"))
    if phi is None or uz is None or mean is None:
        raise TypeError("Phi, Uz, and radial-average Uz derivative grids are all required")
    if log_value is None:
        raise TypeError("the caller's positive F0 log interval is required")
    return phi, uz, mean, log_value


def _looks_like_source_grid(value: Any) -> bool:
    if isinstance(value, Mapping):
        return not ("lower" in value and "upper" in value)
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) < 3:
            return False
        return any(isinstance(row, Mapping) or (isinstance(row, Sequence) and not isinstance(row, (str, bytes, bytearray))) for row in value[:3])
    return False


def _parse_positional_inputs(
    args: Sequence[Any], *, x: Any, Phi_derivatives: Any, Uz_derivatives: Any,
    radial_average_Uz_derivatives: Any, log_F0: Any,
) -> tuple[Any, Any, Any, Any, Any]:
    remaining = list(args)
    if remaining and x is None and not _looks_like_source_grid(remaining[0]):
        x = remaining.pop(0)
    slots = [Phi_derivatives, Uz_derivatives, radial_average_Uz_derivatives]
    for index in range(3):
        if not remaining:
            break
        if slots[index] is not None:
            raise TypeError("a source grid was supplied both positionally and by keyword")
        if not _looks_like_source_grid(remaining[0]):
            break
        slots[index] = remaining.pop(0)
    if remaining:
        if x is None and len(remaining) >= 2:
            x = remaining.pop(0)
        if remaining:
            if log_F0 is not None:
                raise TypeError("log F0 was supplied both positionally and by keyword")
            log_F0 = remaining.pop(0)
    if remaining:
        raise TypeError("too many positional first-collar forcing inputs")
    return x, slots[0], slots[1], slots[2], log_F0


def _power_from_names(term: Mapping[Any, Any]) -> tuple[int, int, int, int] | None:
    powers = term.get("powers")
    if powers is not None:
        if isinstance(powers, Mapping):
            values = tuple(powers.get(name, 0) for name in SOURCE_POWER_BASIS)
        else:
            try:
                values = tuple(powers)
            except TypeError as exc:
                raise TypeError("Laurent powers must be a four-item sequence") from exc
        if len(values) != 4:
            raise ValueError("Laurent powers must have four entries: hb,F0_anchor,Pstar_squared,Lambda")
    else:
        aliases = (
            ("hb_power", "h_b_power", "hb"),
            ("F0_anchor_power", "F0_power", "f0_power", "F0_anchor", "F0"),
            ("Pstar_squared_power", "Pstar_power", "pstar_squared_power", "Pstar_squared", "Pstar2"),
            ("Lambda_power", "lambda_power", "lambda", "Lambda"),
        )
        values = []
        any_found = False
        for keys in aliases:
            value = _mapping_get(term, keys)
            any_found = any_found or value is not None
            values.append(0 if value is None else value)
        if not any_found:
            return None
    result: list[int] = []
    for index, value in enumerate(values):
        if isinstance(value, bool):
            raise TypeError(f"Laurent power {SOURCE_POWER_BASIS[index]} must be an integer")
        try:
            integer = int(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"Laurent power {SOURCE_POWER_BASIS[index]} must be an integer") from exc
        if integer != value:
            raise TypeError(f"Laurent power {SOURCE_POWER_BASIS[index]} must be an integer")
        result.append(integer)
    return tuple(result)  # type: ignore[return-value]


def _term_record(context: Any, value: Any, *, label: str, radial_order: int | None = None,
                 axial_order: int | None = None, index: int = 0) -> dict[str, Any]:
    original = value
    if isinstance(value, Mapping):
        coefficient = value.get("coefficient", value.get("coefficient_interval", value.get("value")))
        powers = _power_from_names(value)
        if coefficient is None and powers is None and "lower" in value and "upper" in value:
            coefficient, powers = value, _ZERO_POWERS
        if coefficient is None and "coefficients" in value:
            coefficients = value["coefficients"]
            if not isinstance(coefficients, Sequence):
                raise TypeError(f"{label} term coefficients must be a sequence")
            coefficient = coefficients[axial_order or 0]
        if coefficient is None:
            raise KeyError(f"{label} Laurent term needs coefficient")
        if powers is None:
            powers = _ZERO_POWERS
        input_origin = value.get("origin")
        if input_origin is None:
            input_origin = {
                "source_ids": value.get("source_ids", []),
                "uncertainty": bool(value.get("uncertainty", False)),
                "error_classes": value.get("error_classes", value.get("uncertainty_class", [])),
            }
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
        if len(values) == 2 and isinstance(values[0], Sequence) and not isinstance(values[0], (str, bytes, bytearray)):
            powers, coefficient, input_origin = _power_from_names({"powers": values[0]}), values[1], None
        elif len(values) == 2:
            powers, coefficient, input_origin = _ZERO_POWERS, value, None
        elif len(values) == 5:
            powers, coefficient, input_origin = _power_from_names({"powers": values[:4]}), values[4], None
        elif len(values) == 3 and isinstance(values[0], Sequence) and not isinstance(values[0], (str, bytes, bytearray)):
            powers, coefficient, input_origin = _power_from_names({"powers": values[0]}), values[1], values[2]
        elif len(values) == 3:
            powers, coefficient, input_origin = (0, int(values[0]), 0, int(values[1])), values[2], None
        else:
            raise TypeError(f"unsupported {label} Laurent term sequence")
    else:
        powers, coefficient, input_origin = _ZERO_POWERS, value, None
    coefficient_iv = _interval(context, coefficient, f"{label} coefficient")
    origin = _origin(
        "source_term",
        input_origin,
        source_ids=(input_origin.get("source_ids", []) if isinstance(input_origin, Mapping) else []),
        label=label,
        rho_order=radial_order,
        Z_order=axial_order,
        term_index=index,
    )
    # Do not retain ``original``: the source receipt remains the caller's
    # object and an origin is intentionally only a lightweight reference.
    del original
    return {"powers": tuple(powers), "coefficient": coefficient_iv, "origin": origin}


def _term_like(value: Any) -> bool:
    if isinstance(value, Mapping):
        return any(key in value for key in ("powers", "coefficient", "coefficient_interval", "value", "terms", "sectors", "source_terms", "hb_power", "F0_power", "Lambda_power"))
    if _is_interval(value):
        return False
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
        return len(values) in (3, 5) or (len(values) == 2 and isinstance(values[0], Sequence))
    return False


def _term_list(context: Any, value: Any, *, label: str, radial_order: int, axial_order: int) -> list[dict[str, Any]]:
    if isinstance(value, Mapping):
        for key in ("terms", "sectors", "sector_terms", "laurent_terms", "source_terms"):
            if key in value:
                entries = value[key]
                if not isinstance(entries, Sequence) or isinstance(entries, (str, bytes, bytearray)):
                    raise TypeError(f"{label} {key} must be a sequence")
                if not entries:
                    return []
                return [_term_record(context, entry, label=label, radial_order=radial_order, axial_order=axial_order, index=index) for index, entry in enumerate(entries)]
        return [_term_record(context, value, label=label, radial_order=radial_order, axial_order=axial_order)]
    if _is_interval(value) or not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        return [_term_record(context, value, label=label, radial_order=radial_order, axial_order=axial_order)]
    values = list(value)
    if not values:
        return []
    if len(values) == 2 and not any(_term_like(item) for item in values):
        return [_term_record(context, value, label=label, radial_order=radial_order, axial_order=axial_order)]
    if not all(isinstance(item, Mapping) and _term_like(item) for item in values) and _term_like(value):
        return [_term_record(context, value, label=label, radial_order=radial_order, axial_order=axial_order)]
    return [_term_record(context, entry, label=label, radial_order=radial_order, axial_order=axial_order, index=index) for index, entry in enumerate(values)]


def _aggregate_terms(context: Any, terms: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Directed-add terms by powers and uncertainty, preserving provenance."""

    buckets: dict[tuple[tuple[int, int, int, int], bool, tuple[Any, ...]], dict[str, Any]] = {}
    for raw in terms:
        powers = tuple(int(value) for value in raw["powers"])
        info = _origin_info(raw.get("origin"))
        error_classes = tuple(sorted(str(value) for value in info.get("error_classes", [])))
        key = (powers, bool(info["uncertainty"]), error_classes)
        if key not in buckets:
            buckets[key] = {
                "powers": powers,
                "coefficient": raw["coefficient"],
                "origin": _origin("sector", raw.get("origin"), source_ids=info["source_ids"], uncertainty=info["uncertainty"], error_classes=error_classes),
            }
        else:
            current = buckets[key]
            current["coefficient"] = current["coefficient"] + raw["coefficient"]
            current["origin"] = _origin("aggregate", current["origin"], raw.get("origin"), power_key=powers, uncertainty=bool(info["uncertainty"]), error_classes=error_classes)
    return list(buckets.values())


class LaurentTaylorJet:
    """Taylor coefficients compacted by formal powers and uncertainty."""

    __slots__ = ("context", "coefficients")

    def __init__(self, context: Any, coefficients: Sequence[Sequence[Mapping[str, Any]]]):
        values = tuple(tuple(dict(term) for term in _aggregate_terms(context, row)) for row in coefficients)
        if not values:
            raise ValueError("a Laurent Taylor jet needs at least its constant coefficient")
        self.context = context
        self.coefficients = values

    @property
    def order(self) -> int:
        return len(self.coefficients) - 1

    def __getitem__(self, index: int) -> tuple[dict[str, Any], ...]:
        return self.coefficients[index]

    @classmethod
    def zero(cls, context: Any, order: int, origin: Any = "zero") -> "LaurentTaylorJet":
        return cls(context, [[] for _ in range(order + 1)])

    @classmethod
    def constant(cls, context: Any, value: Any, order: int, *, powers: Sequence[int] = _ZERO_POWERS, origin: Any = "constant") -> "LaurentTaylorJet":
        return cls(context, [[{"powers": tuple(int(power) for power in powers), "coefficient": _interval(context, value, "Laurent constant"), "origin": _origin("constant", origin)}]] + [[] for _ in range(order)])

    @classmethod
    def variable(cls, context: Any, value: Any, order: int, *, origin: Any = "variable") -> "LaurentTaylorJet":
        if order < 1:
            return cls.constant(context, value, order, origin=origin)
        term0 = {"powers": _ZERO_POWERS, "coefficient": _interval(context, value, "Laurent variable"), "origin": _origin("variable", origin)}
        term1 = {"powers": _ZERO_POWERS, "coefficient": context.mpf(1), "origin": _origin("variable_derivative", origin)}
        return cls(context, [[term0], [term1]] + [[] for _ in range(order - 1)])

    @classmethod
    def from_numeric(cls, jet: IntervalTaylorJet, *, label: str) -> "LaurentTaylorJet":
        rows = []
        for index, coefficient in enumerate(jet.coefficients):
            if coefficient.a == 0 and coefficient.b == 0:
                rows.append([])
            else:
                rows.append([{"powers": _ZERO_POWERS, "coefficient": coefficient, "origin": _origin("numeric_geometry", label=label, Taylor_order=index)}])
        return cls(jet.context, rows)

    def truncate(self, order: int) -> "LaurentTaylorJet":
        if order < 0 or order > self.order:
            raise ValueError("requested Laurent Taylor order is unavailable")
        return type(self)(self.context, self.coefficients[: order + 1])

    def differentiate(self) -> "LaurentTaylorJet":
        if self.order == 0:
            raise ValueError("cannot differentiate a jet with no known derivative")
        return type(self)(self.context, [[{"powers": term["powers"], "coefficient": (index + 1) * term["coefficient"], "origin": _origin("Z_derivative", term["origin"], derivative_order=index + 1)} for term in self[index + 1]] for index in range(self.order)])

    dz = differentiate

    def _coerce(self, other: Any) -> "LaurentTaylorJet":
        if isinstance(other, LaurentTaylorJet):
            if other.context is not self.context:
                raise TypeError("Laurent Taylor jets must use the caller's one interval context")
            return other
        return type(self).constant(self.context, other, self.order, origin="scalar_coercion")

    def __add__(self, other: Any) -> "LaurentTaylorJet":
        other = self._coerce(other)
        order = min(self.order, other.order)
        return type(self)(self.context, [list(self[index]) + list(other[index]) for index in range(order + 1)])

    __radd__ = __add__

    def __neg__(self) -> "LaurentTaylorJet":
        return type(self)(self.context, [[{"powers": term["powers"], "coefficient": -term["coefficient"], "origin": _origin("negate", term["origin"])} for term in row] for row in self.coefficients])

    def __sub__(self, other: Any) -> "LaurentTaylorJet":
        return self + (-self._coerce(other))

    def __rsub__(self, other: Any) -> "LaurentTaylorJet":
        return self._coerce(other) - self

    def __mul__(self, other: Any) -> "LaurentTaylorJet":
        other = self._coerce(other)
        order = min(self.order, other.order)
        rows: list[list[dict[str, Any]]] = []
        for index in range(order + 1):
            row: list[dict[str, Any]] = []
            for left_index in range(index + 1):
                for first in self[left_index]:
                    for second in other[index - left_index]:
                        row.append({"powers": tuple(a + b for a, b in zip(first["powers"], second["powers"])), "coefficient": first["coefficient"] * second["coefficient"], "origin": _origin("multiply", first["origin"], second["origin"], Taylor_order=index)})
            rows.append(row)
        return type(self)(self.context, rows)

    __rmul__ = __mul__

    def __truediv__(self, other: Any) -> "LaurentTaylorJet":
        if isinstance(other, LaurentTaylorJet):
            raise ValueError("formal Laurent denominators are not inverted")
        if isinstance(other, IntervalTaylorJet):
            return self * type(self).from_numeric(other.reciprocal(), label="numeric_denominator")
        denominator = _interval(self.context, other, "numeric denominator")
        _require_nonzero(denominator, "numeric denominator")
        return self * (self.context.mpf(1) / denominator)


FormalLaurentJet = LaurentTaylorJet


def _lift_numeric(jet: IntervalTaylorJet, label: str) -> LaurentTaylorJet:
    return LaurentTaylorJet.from_numeric(jet, label=label)


def _formal_divide_by_numeric(value: LaurentTaylorJet, denominator: IntervalTaylorJet, label: str) -> LaurentTaylorJet:
    _require_nonzero(denominator[0], label)
    return value * _lift_numeric(denominator.reciprocal(), f"reciprocal({label})")


def _source_jet(context: Any, grid: Any, radial_order: int, axial_order: int, label: str) -> LaurentTaylorJet:
    rows = []
    for order in range(axial_order + 1):
        value = _source_value(grid, radial_order, order, label)
        terms = _term_list(context, value, label=label, radial_order=radial_order, axial_order=order)
        rows.append([{"powers": term["powers"], "coefficient": _factorial_div(term["coefficient"], order, context), "origin": _origin("factorial_normalize_source_derivative", term["origin"], factorial=math.factorial(order))} for term in terms])
    return LaurentTaylorJet(context, rows)


def _check_grid_presence(grid: Any, label: str, radial_orders: Sequence[int]) -> None:
    for radial_order in radial_orders:
        for axial_order in range(_MAX_SOURCE_ORDER - radial_order + 1):
            _source_value(grid, radial_order, axial_order, label)


def _axial_second(field: LaurentTaylorJet, field_r: LaurentTaylorJet, field_rr: LaurentTaylorJet, z: LaurentTaylorJet, rho: LaurentTaylorJet, delta: LaurentTaylorJet, weight: LaurentTaylorJet, L_numeric: IntervalTaylorJet) -> LaurentTaylorJet:
    one = LaurentTaylorJet.constant(z.context, 1, z.order, origin="one")
    d = one - z * z
    field_z = field.differentiate()
    field_r_z = field_r.differentiate()
    field_zz = field_z.differentiate()
    A = weight * z * field + d * field_z - 2 * z * rho * field_r
    A_Z = weight * field + (weight - 2) * z * field_z + d * field_zz - 2 * rho * field_r - 2 * z * rho * field_r_z
    A_R = (weight - 2) * z * field_r + d * field_r_z - 2 * z * rho * field_rr
    first = _formal_divide_by_numeric((weight - 1 + delta) * z * A, L_numeric, "L")
    second_inside = _formal_divide_by_numeric(A_Z, L_numeric, "L")
    second_inside = second_inside + _formal_divide_by_numeric(2 * delta * z * A, L_numeric * L_numeric, "L^2")
    second = d * second_inside
    third = _formal_divide_by_numeric(2 * z * rho * A_R, L_numeric, "L")
    return _formal_divide_by_numeric(first + second - third, L_numeric, "L")


def _relative_f0_jet(context: Any, z: IntervalTaylorJet, delta: IntervalTaylorJet, Lambda: IntervalTaylorJet, required_j: Any, order: int) -> IntervalTaylorJet:
    one = IntervalTaylorJet.constant(context, 1, order)
    d = one - z * z
    L = one - delta * z * z
    j = IntervalTaylorJet.constant(context, required_j, order)
    H = z * (one - delta) / 2 + d * (4 * z + j)
    sigma = required_j / context.mpf(500)
    h_square = H * H
    h_square_coefficients = list(h_square.coefficients)
    h_square_coefficients[0] = H[0] ** 2
    denominator = IntervalTaylorJet(context, h_square_coefficients) + sigma * sigma
    _require_nonzero(denominator[0], "H^2+sigma^2")
    ell = -(Lambda * L * H) / denominator
    coefficients = [context.mpf(1)]
    for n in range(1, order + 1):
        value = context.mpf(0)
        for k in range(n):
            value += ell[k] * coefficients[n - 1 - k]
        coefficients.append(value / n)
    return IntervalTaylorJet(context, coefficients)


def _ordinary_terms(jet: LaurentTaylorJet, order: int, label: str) -> list[dict[str, Any]]:
    if jet.order < order:
        raise ValueError(f"{label} only has Taylor order {jet.order}; Z^{order} is unavailable")
    multiplier = math.factorial(order)
    return [{"powers": tuple(term["powers"]), "coefficient": term["coefficient"] * multiplier, "origin": _origin("ordinary_Z_derivative", term["origin"], Z_order=order, family=label)} for term in jet[order]]


def _shift_powers(jet: LaurentTaylorJet, shift: Sequence[int], *, operation: str, family: str) -> LaurentTaylorJet:
    if len(shift) != 4:
        raise ValueError("a Laurent power shift needs four entries")
    shift_tuple = tuple(int(value) for value in shift)
    return LaurentTaylorJet(jet.context, [[{"powers": tuple(a + b for a, b in zip(term["powers"], shift_tuple)), "coefficient": term["coefficient"], "origin": _origin(operation, term["origin"], family=family, power_shift=shift_tuple)} for term in row] for row in jet.coefficients])


def build_first_collar_laurent_forcing_jets(
    c: Any, rho: Any, Z: Any, delta: Any, Lambda: Any, required_j: Any, *positional: Any,
    Phi_derivatives: Any = None, Uz_derivatives: Any = None, radial_average_Uz_derivatives: Any = None,
    log_F0: Any = None, x: Any = None, F0_log: Any = None, logF0: Any = None, F0_log_interval: Any = None,
    mean_derivatives: Any = None, whole_core: Any = None, source: Any = None,
    actual_receipt: Any = None, radial_receipt: Any = None, actual_radial_receipt: Any = None, receipt: Any = None,
    actual: Any = None, actual_radial: Any = None,
    hb: Any = None, Pstar_squared: Any = None, source_power_basis: Any = None, **unknown: Any,
) -> dict[str, Any]:
    """Return compact source-factored ``g4/g5/g6`` ordinary ``Z0..Z3`` rows."""

    if unknown:
        names = ", ".join(sorted(unknown))
        raise TypeError(f"unknown first-collar Laurent forcing input(s): {names}")
    if not hasattr(c, "mpf"):
        raise TypeError("c must be the caller's MPIntervalContext")
    if source_power_basis is not None and tuple(source_power_basis) != SOURCE_POWER_BASIS:
        raise ValueError(f"source_power_basis must be {SOURCE_POWER_BASIS!r}")
    x, Phi_derivatives, Uz_derivatives, radial_average_Uz_derivatives, log_F0 = _parse_positional_inputs(
        positional, x=x, Phi_derivatives=Phi_derivatives, Uz_derivatives=Uz_derivatives,
        radial_average_Uz_derivatives=radial_average_Uz_derivatives, log_F0=log_F0,
    )
    containers = _receipt_containers(whole_core, source, actual_receipt, radial_receipt, actual_radial_receipt, receipt, actual, actual_radial)
    if x is None:
        for container in containers:
            supplied_x = _pick(container, ("x", "radial_x", "sqrt_rho"))
            if supplied_x is not None:
                x = supplied_x
                break
    if x is None:
        raise TypeError("caller-supplied x is required; x is never inferred from rho")
    phi_grid, uz_grid, mean_grid, log_f0 = _resolve_inputs(
        Phi_derivatives=Phi_derivatives, Uz_derivatives=Uz_derivatives,
        radial_average_Uz_derivatives=radial_average_Uz_derivatives, mean_derivatives=mean_derivatives,
        log_F0=log_F0, F0_log=F0_log, logF0=logF0, F0_log_interval=F0_log_interval, containers=containers,
    )

    rho_iv, z_iv, delta_iv = _interval(c, rho, "rho"), _interval(c, Z, "Z"), _interval(c, delta, "delta")
    lambda_iv, j_iv, x_iv, log_f0_iv = _interval(c, Lambda, "Lambda"), _interval(c, required_j, "required_j"), _interval(c, x, "x"), _interval(c, log_f0, "log F0")
    _require_domain(rho_iv, _RHO_MIN, _RHO_MAX, c, "rho")
    _require_domain(z_iv, _Z_MIN, _Z_MAX, c, "Z")
    _require_positive(lambda_iv, "Lambda")
    _require_positive(j_iv, "required_j")
    _require_positive(x_iv, "x")
    x_squared = x_iv ** 2
    if x_squared.a > rho_iv.a or x_squared.b < rho_iv.b:
        raise ValueError("x must enclose the positive square root of the whole rho interval")

    # Complete derivative admission is checked before any algebra.
    _check_grid_presence(phi_grid, "Phi_derivatives", (0, 1, 2))
    _check_grid_presence(uz_grid, "Uz_derivatives", (0, 1, 2))
    _check_grid_presence(mean_grid, "radial_average_Uz_derivatives", (0, 1))
    degree = _MAX_SOURCE_ORDER
    rho_num = IntervalTaylorJet.constant(c, rho_iv, degree)
    z_num = IntervalTaylorJet.variable(c, z_iv, degree)
    delta_num = IntervalTaylorJet.constant(c, delta_iv, degree)
    lambda_num = IntervalTaylorJet.constant(c, lambda_iv, degree)
    one_num = IntervalTaylorJet.constant(c, 1, degree)
    d_num = one_num - z_num * z_num
    L_num = one_num - delta_num * z_num * z_num
    _require_nonzero(L_num[0], "L=1-delta*Z^2")
    rho_s, z_s, delta_s = _lift_numeric(rho_num, "rho"), _lift_numeric(z_num, "Z"), _lift_numeric(delta_num, "delta")
    one_s, d_s = _lift_numeric(one_num, "one"), _lift_numeric(d_num, "1-Z^2")
    U0, U1, U2 = _source_jet(c, uz_grid, 0, 5, "Uz_derivatives"), _source_jet(c, uz_grid, 1, 4, "Uz_derivatives"), _source_jet(c, uz_grid, 2, 3, "Uz_derivatives")
    M0, M1 = _source_jet(c, mean_grid, 0, 5, "radial_average_Uz_derivatives"), _source_jet(c, mean_grid, 1, 4, "radial_average_Uz_derivatives")
    Phi0, Phi1, Phi2 = _source_jet(c, phi_grid, 0, 5, "Phi_derivatives"), _source_jet(c, phi_grid, 1, 4, "Phi_derivatives"), _source_jet(c, phi_grid, 2, 3, "Phi_derivatives")
    M0_Z, M1_Z, U1_Z = M0.differentiate(), M1.differentiate(), U1.differentiate()
    Q0 = _formal_divide_by_numeric(2 * z_s * U0 - (1 - delta_s) * z_s * M0 - d_s * M0_Z, L_num, "L")
    Qr = _formal_divide_by_numeric(2 * z_s * U1 - (1 - delta_s) * z_s * M1 - d_s * M1_Z, L_num, "L")
    rhoQrr = _formal_divide_by_numeric(2 * z_s * rho_s * U2 - (1 - delta_s) * z_s * (U1 - 2 * M1) - d_s * (U1_Z - 2 * M1_Z), L_num, "L")
    Q0_Z = Q0.differentiate()
    S = _formal_divide_by_numeric(Q0 + (1 - delta_s) * z_s * Q0_Z / 2 + rho_s * Qr, L_num, "L")
    S = S + Q0 * (Q0 / 2 + rho_s * Qr)
    S = S + _formal_divide_by_numeric(U0 * (-2 * z_s * Q0 + d_s * Q0_Z - 2 * z_s * rho_s * Qr), L_num, "L")
    P = 2 * Qr + rhoQrr
    x_formal = _lift_numeric(IntervalTaylorJet.constant(c, x_iv, degree), "caller_x")
    g4_power0, g4_power1 = 2 * x_formal * P, -x_formal * S
    relative_f0 = _lift_numeric(_relative_f0_jet(c, z_num, delta_num, lambda_num, j_iv, degree), "relative_F0")
    F0_ratio, F1_ratio, F2_ratio = (Phi0 * relative_f0).truncate(5), (Phi1 * relative_f0).truncate(4), (Phi2 * relative_f0).truncate(3)
    weight_F = _lift_numeric(IntervalTaylorJet.constant(c, -2, degree) - delta_num, "weight_F")
    weight_U = _lift_numeric(IntervalTaylorJet.constant(c, -1, degree) - delta_num, "weight_U")
    Nt_over_F0_anchor = -_axial_second(F0_ratio, F1_ratio, F2_ratio, z_s, rho_s, delta_s, weight_F, L_num)
    Nz = -_axial_second(U0, U1, U2, z_s, rho_s, delta_s, weight_U, L_num)
    g6_power1 = 2 * Nz - _formal_divide_by_numeric(4 * rho_s * z_s * P, L_num, "L")
    g6_power2 = _formal_divide_by_numeric(2 * rho_s * z_s * S, L_num, "L")

    families = {
        "4": ((g4_power0, (0, 0, 0, 0), "g4_Lambda0"), (g4_power1, (0, 0, 0, -1), "g4_Lambda-1")),
        "5": ((2 * Nt_over_F0_anchor, (0, 1, 0, -1), "g5_F0_anchor_Lambda-1"),),
        "6": ((g6_power1, (0, 0, 0, -1), "g6_Lambda-1"), (g6_power2, (0, 0, 0, -2), "g6_Lambda-2")),
    }
    sectors: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for order in _OUTPUT_ORDERS:
        rows: dict[str, list[dict[str, Any]]] = {}
        for family, entries in families.items():
            emitted: list[dict[str, Any]] = []
            for jet, shift, sector_label in entries:
                emitted.extend(_ordinary_terms(_shift_powers(jet, shift, operation="external_formal_power", family=sector_label), order, f"g{family}"))
            rows[family] = emitted
        sectors[str(order)] = rows

    # Dynamic raw B0 and B1 use the same actual input jets and the very same
    # relative F0 jet as the forcing. Only the external epsilon and F0 atoms
    # are shifted; no absolute F0 or hb value is materialized.
    F = _shift_powers(F0_ratio, (0,1,0,0), operation="formal_F0", family="F")
    Fr = _shift_powers(F1_ratio, (0,1,0,0), operation="formal_F0", family="Fr")
    Fz = F.differentiate()
    Uz = U0.differentiate()
    H = z_s*(1-delta_s)/2 + d_s*U0
    B = Q0 + _formal_divide_by_numeric(1-2*z_s*U0,L_num,"L")
    E = z_s*U0 - c.mpf(1)/2
    Cf = F + rho_s*Fr
    Zf = _formal_divide_by_numeric((-2-delta_s)*z_s*F+d_s*Fz-2*z_s*rho_s*Fr,L_num,"L")
    Zw = _formal_divide_by_numeric((-1-delta_s)*z_s*U0+d_s*Uz-2*z_s*rho_s*U1,L_num,"L")
    C0 = {
        "4,1": (4*x_formal*F,1),
        "5,1": (2*(Q0+_formal_divide_by_numeric(E*(-2+delta_s),L_num,"L")),1),
        "5,2": (2*Zf-_formal_divide_by_numeric(2*(-1+delta_s)*z_s*Cf,L_num,"L"),1),
        "5,3": (-_formal_divide_by_numeric(2*(1+delta_s)*z_s*Cf,L_num,"L"),1),
        "5,5": (x_formal*B,1),
        "6,1": (-_formal_divide_by_numeric(8*z_s*rho_s*F,L_num,"L"),2),
        "6,2": (2*(_formal_divide_by_numeric(E*(-1+delta_s),L_num,"L")+Zw
                   -_formal_divide_by_numeric(U1*z_s*(-1+delta_s)*rho_s,L_num,"L")),1),
        "6,3": (-_formal_divide_by_numeric(2*(1+delta_s)*z_s*rho_s*U1,L_num,"L"),1),
        "6,4": (-_formal_divide_by_numeric(4*z_s,L_num,"L"),1),
        "6,6": (x_formal*B,1),
    }
    B1 = {
        "5,1": (_formal_divide_by_numeric(2*H,L_num,"L"),1),
        "5,2": (-_formal_divide_by_numeric(2*d_s*Cf,L_num,"L"),1),
        "5,3": (-_formal_divide_by_numeric(2*d_s*Cf,L_num,"L"),1),
        "6,2": (_formal_divide_by_numeric(2*(H-U1*d_s*rho_s),L_num,"L"),1),
        "6,3": (-_formal_divide_by_numeric(2*d_s*rho_s*U1,L_num,"L"),1),
        "6,4": (_formal_divide_by_numeric(2*d_s,L_num,"L"),1),
    }
    matrix={}
    for label,entries in (("C0",C0),("B1",B1)):
        matrix[label]={str(k):{key:_ordinary_terms(_shift_powers(value,(0,0,0,-ep),
            operation="external_epsilon",family=label+"_"+key),k,label+"_"+key)
            for key,(value,ep) in entries.items()} for k in range(3)}

    proof = {
        "kind": "actual_collar_forcing_compact_v1",
        "same_source_dynamic_C0_and_B1_Z0_to_Z2": True,
        "F0_relative_derivative_applied_once": True,
        "C0_61_has_two_external_epsilon_factors": True,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(_OUTPUT_ORDERS),
        "ordinary_Z_derivative_order": 3,
        "source_derivatives_are_ordinary": True,
        "Taylor_coefficients_divided_by_factorial": True,
        "unknown_high_derivatives_padded": False,
        "source_grid_total_order_cap": _MAX_SOURCE_ORDER,
        "required_source_rows": {"Phi": "rho^i Z^k for i=0,1,2 and 0<=k<=5-i", "U": "rho^i Z^k for i=0,1,2 and 0<=k<=5-i", "mean": "rho^i Z^k for i=0,1 and 0<=k<=5-i"},
        "compact_aggregation_key": "(powers, uncertainty, error_classes)",
        "aggregation_is_directed_interval_addition_only": True,
        "source_ids_union_retained": True,
        "aggregate_path_count_retained": True,
        "uncertainty_error_classes_retained": True,
        "uncertainty_never_mixed_with_source": True,
        "cauchy_or_symmetric_error_requires_uncertainty": True,
        "no_source_function_identity_or_cancellation": True,
        "rho_domain": "[4,4.1]",
        "Z_domain": "[-1,1]",
        "rho_zero_safe": True,
        "rho_zero_rule": "rho*m_rr = Ur-2*mr; no division by rho",
        "interval_outward_rounding": True,
        "denominator_checks": ["L=1-delta*Z^2 excludes zero", "H^2+sigma^2 excludes zero"],
        "formal_laurent_terms_compacted": True,
        "hb_negative_powers_preserved": True,
        "input_source_provenance_retained": True,
        "input_scalar_intervals_promoted_to_zero_sector": True,
        "F0_power_is_formal": True,
        "F0_current_Z_factor_formal": True,
        "F0_log_interval": log_f0_iv,
        "F0_upper_not_used_as_value": True,
        "F0_absolute_exp_materialized": False,
        "Pstar_squared_materialized": False,
        "hb_materialized": False,
        "x_caller_supplied": True,
        "x_inferred_from_rho": False,
        "positive_x_squared_encloses_rho": True,
        "g5_anchor_rule": "F(Z+h)/F0(Z) jet at each current Z; no extra D_F",
        "Lambda_inside_relative_F0_derivative_coefficients_is_enclosed_numerically": True,
        "outer_Lambda_power_remains_formal": True,
        "g4_formula": "2*epsilon*x*(-S/2+Lambda*P), epsilon=Lambda^-1",
        "g5_formula": "2*epsilon*(-axial_second(F))",
        "g6_formula": "2*epsilon*Nz-4*z*rho*epsilon^2*Np/L",
        "Q_formula": "Q0=(2*z*U0-(1-delta)*z*M0-d*M0_Z)/L",
        "Qr_formula": "Qr=(2*z*U1-(1-delta)*z*M1-d*M1_Z)/L",
        "rhoQrr_formula": "rho*Qrr=(2*z*rho*U2-(1-delta)*z*(U1-2*M1)-d*(U1_Z-2*M1_Z))/L",
        "P_formula": "P=2*Qr+rho*Qrr",
        "S_formula": "(Q0+(1-delta)*z*Q0_Z/2+rho*Qr)/L+Q0*(Q0/2+rho*Qr)+U0*(-2*z*Q0+d*Q0_Z-2*z*rho*Qr)/L",
        "axial_second_formula": "[(a-1+delta)*z*A/L+d*(AZ/L+2*delta*z*A/L^2)-2*z*rho*AR/L]/L",
        "weights": {"F": "-2-delta", "U": "-1-delta"},
        "x_formula": "caller-supplied x interval; no sqrt(rho) inference",
        "actual_radial_receipt_consumed": bool(actual_receipt is not None or radial_receipt is not None or actual_radial_receipt is not None or receipt is not None or actual is not None or actual_radial is not None),
        "hb_argument_supplied": hb is not None,
        "Pstar_squared_argument_supplied": Pstar_squared is not None,
    }
    jets = {"Q0": Q0, "Qr": Qr, "rhoQrr": rhoQrr, "S": S, "P": P, "Np_formal": -S / 2 + _lift_numeric(lambda_num, "internal_Lambda") * P, "Nt_over_F0_anchor": Nt_over_F0_anchor, "Nz": Nz, "relative_F0": relative_f0, "x": x_formal}
    by_family = {f"g{family}": {f"Z{order}": sectors[str(order)][family] for order in _OUTPUT_ORDERS} for family in families}
    flat_rows = {f"g{family}_Z{order}": sectors[str(order)][family] for family in families for order in _OUTPUT_ORDERS}
    result: dict[str, Any] = {"matrix": matrix, "sectors": sectors, "rows": sectors, "by_order": sectors, "source_power_basis": list(SOURCE_POWER_BASIS), "proof": proof, "metadata": proof, "jets": jets, "log_F0": log_f0_iv, "F0_log": log_f0_iv, "F0_log_interval": log_f0_iv, **by_family, **flat_rows}
    result.update(sectors)
    return result


def first_collar_laurent_forcing_jets(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def derive_first_collar_laurent_forcing_jets(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def build_first_collar_laurent_forcing(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def first_collar_laurent_forcing(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def build_collar_laurent_forcing(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def build_first_collar_laurent_forcing_z3(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def first_collar_laurent_forcing_z3(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def build_first_collar_forcing_laurent_jets(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def first_collar_forcing_laurent_jets(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def build_first_collar_forcing_laurent(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def first_collar_forcing_laurent(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build_first_collar_laurent_forcing_jets(*args, **kwargs)

