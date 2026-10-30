"""Laurent source sectors for the first-collar ``g4/g5/g6`` forcing.

This module is an independent source leaf.  It repeats the finite Taylor
algebra needed by the first-collar forcing equations, but it does not import
or execute an older forcing implementation.  The caller supplies ordinary
mixed derivative grids for ``Phi``, ``U`` and the radial mean.  A grid entry
may be a scalar interval (which is promoted to the zero source sector) or a
list of signed Laurent terms.  A Laurent term has the form::

    {"powers": (hb, F0_anchor, Pstar_squared, Lambda),
     "coefficient": interval,
     "origin": ...}

The four powers are formal integer exponents.  In particular, negative
``hb`` exponents are accepted and are never discarded.  Products form a
Cartesian product of terms and additions concatenate terms; equal powers are
*not* combined.  This leaves every source contribution and its provenance
visible to the downstream factored Green adapter.

The returned forcing sectors use the order-major shape::

    result["sectors"]["0"]["4"] -> list of terms

where the row labels are ``"4"``, ``"5"`` and ``"6"`` and the orders are
ordinary ``Z0`` through ``Z3`` derivatives.  ``result["rows"]`` is the same
mapping.  Compatibility views ``g4``, ``g5`` and ``g6`` are also provided.

The relative positive ``F0`` jet follows the same local log-derivative recipe
as the source equations.  ``Lambda`` is used as a numerical interval only
inside that relative jet.  The external Lambda factors of ``g4/g5/g6`` are
formal exponent shifts in the fourth power slot.  No exponentiation of the
caller log-``F0`` interval, no positive amplitude bound, and no ``Pstar`` or
``hb`` value is used as a field value.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


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


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
_ZERO_POWERS = (0, 0, 0, 0)
_MAX_SOURCE_ORDER = 5
_OUTPUT_ORDERS = (0, 1, 2, 3)
_RHO_MIN = "4"
_RHO_MAX = "4.1"
_Z_MIN = "-1"
_Z_MAX = "1"


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") or (hasattr(value, "a") and hasattr(value, "b"))


def _bounds(value: Any) -> tuple[Any, Any]:
    """Return directed endpoints without selecting a representative."""

    if not _is_interval(value):
        raise TypeError("an mpmath interval is required")
    try:
        return value.a, value.b
    except AttributeError:  # pragma: no cover - defensive for iv variants
        rows = value._mpi_
        return rows[0], rows[1]


def _interval(context: Any, value: Any, label: str = "interval") -> Any:
    """Coerce a scalar, endpoint pair, or interval with directed rounding."""

    if _is_interval(value):
        lo, hi = _bounds(value)
        return context.mpf([lo, hi])
    if isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            return context.make_mpf(
                (
                    tuple(value["lower_exact_mpf_tuple"]),
                    tuple(value["upper_exact_mpf_tuple"]),
                )
            )
        if "lower" in value and "upper" in value:
            return context.mpf([value["lower"], value["upper"]])
        raise TypeError(f"{label} needs lower/upper endpoints")
    if isinstance(value, (list, tuple)) and len(value) == 2:
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


def _origin(operation: str, *inputs: Any, **metadata: Any) -> dict[str, Any]:
    """Build small nested provenance records without evaluating source atoms."""

    result: dict[str, Any] = {"operation": operation}
    if inputs:
        result["inputs"] = inputs
    result.update(metadata)
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
        zero = context.mpf(0)
        return cls(context, [value] + [zero] * order)

    @classmethod
    def variable(cls, context: Any, value: Any, order: int) -> "IntervalTaylorJet":
        if order < 1:
            return cls.constant(context, value, order)
        value = _interval(context, value, "Taylor variable")
        zero = context.mpf(0)
        return cls(context, [value, context.mpf(1)] + [zero] * (order - 1))

    def truncate(self, order: int) -> "IntervalTaylorJet":
        if order < 0:
            raise ValueError("nonnegative Taylor order required")
        if order > self.order:
            raise ValueError("cannot pad a Taylor jet with unknown derivatives")
        return IntervalTaylorJet(self.context, self.coefficients[: order + 1])

    def differentiate(self) -> "IntervalTaylorJet":
        if self.order == 0:
            raise ValueError("cannot differentiate a jet with no known derivative")
        return IntervalTaylorJet(
            self.context,
            [(index + 1) * self.coefficients[index + 1] for index in range(self.order)],
        )

    dz = differentiate

    def __getitem__(self, index: int) -> Any:
        return self.coefficients[index]

    def _coerce(self, other: Any) -> "IntervalTaylorJet":
        if isinstance(other, IntervalTaylorJet):
            if other.context is not self.context:
                raise TypeError("Taylor jets must use the caller's one interval context")
            return other
        return IntervalTaylorJet.constant(self.context, other, self.order)

    def __add__(self, other: Any) -> "IntervalTaylorJet":
        other = self._coerce(other)
        order = min(self.order, other.order)
        return IntervalTaylorJet(
            self.context,
            [self[index] + other[index] for index in range(order + 1)],
        )

    __radd__ = __add__

    def __neg__(self) -> "IntervalTaylorJet":
        return IntervalTaylorJet(self.context, [-value for value in self.coefficients])

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
        return IntervalTaylorJet(self.context, coefficients)

    __rmul__ = __mul__

    def reciprocal(self) -> "IntervalTaylorJet":
        _require_nonzero(self[0], "Taylor denominator")
        context = self.context
        coefficients = [context.mpf(1) / self[0]]
        for index in range(1, self.order + 1):
            value = context.mpf(0)
            for left in range(1, index + 1):
                value += self[left] * coefficients[index - left]
            coefficients.append(-value / self[0])
        return IntervalTaylorJet(context, coefficients)

    def __truediv__(self, other: Any) -> "IntervalTaylorJet":
        return self * self._coerce(other).reciprocal()

    def __rtruediv__(self, other: Any) -> "IntervalTaylorJet":
        return self._coerce(other) * self.reciprocal()

    def __pow__(self, power: int) -> "IntervalTaylorJet":
        if not isinstance(power, int) or power < 0:
            raise ValueError("only nonnegative integer Taylor powers are supported")
        result = IntervalTaylorJet.constant(self.context, 1, self.order)
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

    if isinstance(grid, Mapping):
        names = (
            (radial_order, axial_order),
            f"rho{radial_order}_Z{axial_order}",
            f"rho{radial_order}_z{axial_order}",
            f"rho_{radial_order}_Z_{axial_order}",
            f"rho{radial_order}Z{axial_order}",
            f"{radial_order},{axial_order}",
            f"{radial_order}:{axial_order}",
        )
        direct = _mapping_get(grid, names)
        if direct is not None:
            return direct

        radial_names = (
            radial_order,
            str(radial_order),
            f"rho{radial_order}",
            f"rho_{radial_order}",
        )
        row = _mapping_get(grid, radial_names)
        if isinstance(row, Mapping):
            axial_names = (axial_order, str(axial_order), f"Z{axial_order}", f"z{axial_order}")
            nested = _mapping_get(row, axial_names)
            if nested is not None:
                return nested

    if isinstance(grid, Sequence) and not isinstance(grid, (str, bytes, bytearray)):
        if radial_order < len(grid):
            row = grid[radial_order]
            if isinstance(row, Sequence) and not isinstance(row, (str, bytes, bytearray)):
                if axial_order < len(row):
                    return row[axial_order]

    raise KeyError(f"missing {label} ordinary derivative rho^{radial_order} Z^{axial_order}")


def _pick(mapping: Mapping[Any, Any], names: Sequence[str]) -> Any | None:
    return _mapping_get(mapping, names)


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
    whole_core: Any,
    source: Any,
) -> tuple[Any, Any, Any, Any]:
    containers = [value for value in (whole_core, source) if isinstance(value, Mapping)]

    def first(value: Any, names: Sequence[str]) -> Any:
        if value is not None:
            return value
        for container in containers:
            found = _pick(container, names)
            if found is not None:
                return found
        return None

    phi = first(Phi_derivatives, ("Phi_derivatives", "Phi", "phi_derivatives"))
    uz = first(Uz_derivatives, ("Uz_derivatives", "Uz", "U_derivatives", "uz_derivatives"))
    mean = first(
        radial_average_Uz_derivatives,
        ("radial_average_Uz_derivatives", "mean_derivatives", "M_derivatives", "mean"),
    )
    if mean is None:
        mean = first(mean_derivatives, ("mean_derivatives", "M_derivatives", "mean"))

    log_value = first(log_F0, ("log_F0", "F0_log", "logF0", "F0_log_interval"))
    if log_value is None:
        log_value = first(F0_log, ("log_F0", "F0_log", "logF0", "F0_log_interval"))
    if log_value is None:
        log_value = first(logF0, ("log_F0", "F0_log", "logF0", "F0_log_interval"))
    if log_value is None:
        log_value = first(
            F0_log_interval,
            ("log_F0", "F0_log", "logF0", "F0_log_interval"),
        )
    if log_value is None:
        for container in containers:
            log_value = _pick(
                container,
                ("F0_exact_positive_log_enclosure", "F0_log_enclosure", "log_F0_interval"),
            )
            if log_value is not None:
                break

    if phi is None or uz is None or mean is None:
        raise TypeError("Phi, Uz, and radial-average Uz derivative grids are all required")
    if log_value is None:
        raise TypeError("the caller's positive F0 log interval is required")
    return phi, uz, mean, log_value


def _looks_like_source_grid(value: Any) -> bool:
    """Distinguish a positional grid from an endpoint pair used for ``x``."""

    if isinstance(value, Mapping):
        return not ("lower" in value and "upper" in value)
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) < 3:
            return False
        return any(
            isinstance(row, Mapping)
            or (isinstance(row, Sequence) and not isinstance(row, (str, bytes, bytearray)))
            for row in value[:3]
        )
    return False


def _parse_positional_inputs(
    args: Sequence[Any],
    *,
    x: Any,
    Phi_derivatives: Any,
    Uz_derivatives: Any,
    radial_average_Uz_derivatives: Any,
    log_F0: Any,
) -> tuple[Any, Any, Any, Any, Any]:
    """Accept old natural positional source forms without changing semantics."""

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
        found = []
        any_found = False
        for keys in aliases:
            value = _mapping_get(term, keys)
            if value is not None:
                any_found = True
            found.append(0 if value is None else value)
        if not any_found:
            return None
        values = tuple(found)
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


def _term_record(
    context: Any,
    value: Any,
    *,
    label: str,
    radial_order: int | None = None,
    axial_order: int | None = None,
    index: int = 0,
) -> dict[str, Any]:
    """Normalize one source term while retaining the complete input record."""

    original: Any = value
    if isinstance(value, Mapping):
        coefficient = value.get("coefficient", value.get("value"))
        powers = _power_from_names(value)
        if coefficient is None and powers is None and "lower" in value and "upper" in value:
            # The ordinary scalar interval spelling is also accepted for a
            # grid cell.  It carries no formal source powers.
            coefficient = value
            powers = _ZERO_POWERS
        if coefficient is None and "coefficients" in value:
            coefficients = value["coefficients"]
            if not isinstance(coefficients, Sequence):
                raise TypeError(f"{label} term coefficients must be a sequence")
            # A source grid ordinarily has one coefficient per source cell;
            # accepting a sequence here makes the term usable by callers that
            # provide all axial entries in one record, without padding a
            # missing derivative row.
            coefficient = coefficients[axial_order or 0]
        if coefficient is None:
            raise KeyError(f"{label} Laurent term needs coefficient")
        if powers is None:
            # A mapping containing only an interval-like lower/upper pair is
            # handled by _term_list, not here.  A plain coefficient mapping is
            # intentionally accepted as the zero formal sector for backwards
            # compatibility with the scalar source leaf.
            powers = _ZERO_POWERS
        input_origin = value.get("origin")
        if input_origin is None:
            input_origin = {
                key: item
                for key, item in value.items()
                if key not in {"coefficient", "value", "coefficients", "powers"}
            }
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
        if len(values) == 2 and isinstance(values[0], Sequence) and not isinstance(values[0], (str, bytes, bytearray)):
            powers = _power_from_names({"powers": values[0]})
            coefficient = values[1]
            input_origin = None
        elif len(values) == 2:
            # Endpoint pairs are scalar interval input, not an implicit
            # two-term source list.  The caller can still provide two
            # Laurent terms by wrapping each term record in a list whose
            # entries are themselves term-like.
            powers = _ZERO_POWERS
            coefficient = value
            input_origin = None
        elif len(values) == 5:
            powers = _power_from_names({"powers": values[:4]})
            coefficient = values[4]
            input_origin = None
        elif len(values) == 3 and isinstance(values[0], Sequence) and not isinstance(values[0], (str, bytes, bytearray)):
            powers = _power_from_names({"powers": values[0]})
            coefficient = values[1]
            input_origin = values[2]
        elif len(values) == 3:
            # Compatibility with old (F0_power, Lambda_power, coefficient)
            # rows.  The omitted hb and Pstar powers are exactly zero.
            powers = (0, int(values[0]), 0, int(values[1]))
            coefficient = values[2]
            input_origin = None
        else:
            raise TypeError(f"unsupported {label} Laurent term sequence")
    else:
        powers = _ZERO_POWERS
        coefficient = value
        input_origin = None

    coefficient_iv = _interval(context, coefficient, f"{label} coefficient")
    source_origin = _origin(
        "source_term",
        input_origin,
        label=label,
        rho_order=radial_order,
        Z_order=axial_order,
        term_index=index,
        original_record=original,
    )
    return {"powers": tuple(powers), "coefficient": coefficient_iv, "origin": source_origin}


def _term_like(value: Any) -> bool:
    if isinstance(value, Mapping):
        return any(key in value for key in ("powers", "coefficient", "value", "coefficients", "hb_power", "F0_power", "Lambda_power"))
    if _is_interval(value):
        return False
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
        return len(values) in (3, 5) or (len(values) == 2 and isinstance(values[0], Sequence))
    return False


def _term_list(
    context: Any,
    value: Any,
    *,
    label: str,
    radial_order: int,
    axial_order: int,
) -> list[dict[str, Any]]:
    """Return a nonempty list of Laurent terms for one grid cell."""

    if isinstance(value, Mapping):
        for key in ("terms", "sectors", "laurent_terms", "source_terms"):
            if key in value:
                entries = value[key]
                if not isinstance(entries, Sequence) or isinstance(entries, (str, bytes, bytearray)):
                    raise TypeError(f"{label} {key} must be a sequence")
                if not entries:
                    raise ValueError(f"{label} {key} may not be empty")
                return [
                    _term_record(
                        context,
                        entry,
                        label=label,
                        radial_order=radial_order,
                        axial_order=axial_order,
                        index=index,
                    )
                    for index, entry in enumerate(entries)
                ]
        return [
            _term_record(
                context,
                value,
                label=label,
                radial_order=radial_order,
                axial_order=axial_order,
            )
        ]
    if _is_interval(value) or not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        return [
            _term_record(
                context,
                value,
                label=label,
                radial_order=radial_order,
                axial_order=axial_order,
            )
        ]

    values = list(value)
    if not values:
        # An explicitly present empty sector list denotes exact zero; a
        # missing derivative key is rejected by _source_value before here.
        return []
    # A numeric pair is an endpoint interval, not a two-term Laurent list.
    if len(values) == 2 and not any(_term_like(item) for item in values):
        return [
            _term_record(
                context,
                value,
                label=label,
                radial_order=radial_order,
                axial_order=axial_order,
            )
        ]
    if not all(isinstance(item, Mapping) and _term_like(item) for item in values) and _term_like(value):
        return [
            _term_record(
                context,
                value,
                label=label,
                radial_order=radial_order,
                axial_order=axial_order,
            )
        ]
    return [
        _term_record(
            context,
            entry,
            label=label,
            radial_order=radial_order,
            axial_order=axial_order,
            index=index,
        )
        for index, entry in enumerate(values)
    ]


class LaurentTaylorJet:
    """Factorial-normalised Taylor coefficients carrying unmerged sectors."""

    __slots__ = ("context", "coefficients")

    def __init__(self, context: Any, coefficients: Sequence[Sequence[Mapping[str, Any]]]):
        values = tuple(tuple(dict(term) for term in row) for row in coefficients)
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
    def constant(
        cls,
        context: Any,
        value: Any,
        order: int,
        *,
        powers: Sequence[int] = _ZERO_POWERS,
        origin: Any = "constant",
    ) -> "LaurentTaylorJet":
        coefficient = _interval(context, value, "Laurent constant")
        term = {"powers": tuple(int(power) for power in powers), "coefficient": coefficient, "origin": origin}
        return cls(context, [[term]] + [[] for _ in range(order)])

    @classmethod
    def variable(cls, context: Any, value: Any, order: int, *, origin: Any = "variable") -> "LaurentTaylorJet":
        if order < 1:
            return cls.constant(context, value, order, origin=origin)
        coefficient = _interval(context, value, "Laurent variable")
        one = context.mpf(1)
        term0 = {"powers": _ZERO_POWERS, "coefficient": coefficient, "origin": origin}
        term1 = {"powers": _ZERO_POWERS, "coefficient": one, "origin": _origin("variable_derivative", origin)}
        return cls(context, [[term0], [term1]] + [[] for _ in range(order - 1)])

    @classmethod
    def from_numeric(cls, jet: IntervalTaylorJet, *, label: str) -> "LaurentTaylorJet":
        coefficients = []
        for index, coefficient in enumerate(jet.coefficients):
            if coefficient.a == 0 and coefficient.b == 0:
                coefficients.append([])
                continue
            coefficients.append(
                [
                    {
                        "powers": _ZERO_POWERS,
                        "coefficient": coefficient,
                        "origin": _origin("numeric_geometry", label=label, Taylor_order=index),
                    }
                ]
            )
        return cls(jet.context, coefficients)

    def truncate(self, order: int) -> "LaurentTaylorJet":
        if order < 0:
            raise ValueError("nonnegative Taylor order required")
        if order > self.order:
            raise ValueError("cannot pad a Laurent Taylor jet with unknown derivatives")
        return LaurentTaylorJet(self.context, self.coefficients[: order + 1])

    def differentiate(self) -> "LaurentTaylorJet":
        if self.order == 0:
            raise ValueError("cannot differentiate a jet with no known derivative")
        coefficients: list[list[dict[str, Any]]] = []
        for index in range(self.order):
            row = []
            for term in self[index + 1]:
                row.append(
                    {
                        "powers": term["powers"],
                        "coefficient": (index + 1) * term["coefficient"],
                        "origin": _origin("Z_derivative", term["origin"], derivative_order=index + 1),
                    }
                )
            coefficients.append(row)
        return LaurentTaylorJet(self.context, coefficients)

    dz = differentiate

    def _coerce(self, other: Any) -> "LaurentTaylorJet":
        if isinstance(other, LaurentTaylorJet):
            if other.context is not self.context:
                raise TypeError("Laurent Taylor jets must use the caller's one interval context")
            return other
        return LaurentTaylorJet.constant(self.context, other, self.order, origin="scalar_coercion")

    def __add__(self, other: Any) -> "LaurentTaylorJet":
        other = self._coerce(other)
        order = min(self.order, other.order)
        coefficients = []
        for index in range(order + 1):
            # Deliberately concatenate.  Equal powers are distinct signed
            # source sectors and are not algebraically collapsed here.
            coefficients.append(
                [
                    {
                        "powers": term["powers"],
                        "coefficient": term["coefficient"],
                        "origin": _origin("add", term["origin"]),
                    }
                    for term in self[index]
                ]
                + [
                    {
                        "powers": term["powers"],
                        "coefficient": term["coefficient"],
                        "origin": _origin("add", term["origin"]),
                    }
                    for term in other[index]
                ]
            )
        return LaurentTaylorJet(self.context, coefficients)

    __radd__ = __add__

    def __neg__(self) -> "LaurentTaylorJet":
        return LaurentTaylorJet(
            self.context,
            [
                [
                    {
                        "powers": term["powers"],
                        "coefficient": -term["coefficient"],
                        "origin": _origin("negate", term["origin"]),
                    }
                    for term in row
                ]
                for row in self.coefficients
            ],
        )

    def __sub__(self, other: Any) -> "LaurentTaylorJet":
        return self + (-self._coerce(other))

    def __rsub__(self, other: Any) -> "LaurentTaylorJet":
        return self._coerce(other) - self

    def __mul__(self, other: Any) -> "LaurentTaylorJet":
        other = self._coerce(other)
        order = min(self.order, other.order)
        coefficients: list[list[dict[str, Any]]] = []
        for index in range(order + 1):
            row: list[dict[str, Any]] = []
            for left in range(index + 1):
                for first in self[left]:
                    for second in other[index - left]:
                        powers = tuple(a + b for a, b in zip(first["powers"], second["powers"]))
                        row.append(
                            {
                                "powers": powers,
                                "coefficient": first["coefficient"] * second["coefficient"],
                                "origin": _origin("multiply", first["origin"], second["origin"], Taylor_order=index),
                            }
                        )
            coefficients.append(row)
        return LaurentTaylorJet(self.context, coefficients)

    __rmul__ = __mul__

    def __truediv__(self, other: Any) -> "LaurentTaylorJet":
        if isinstance(other, LaurentTaylorJet):
            raise ValueError("formal Laurent denominators are not inverted")
        if isinstance(other, IntervalTaylorJet):
            return self * LaurentTaylorJet.from_numeric(other.reciprocal(), label="numeric_denominator")
        denominator = _interval(self.context, other, "numeric denominator")
        _require_nonzero(denominator, "numeric denominator")
        return self * (self.context.mpf(1) / denominator)


FormalLaurentJet = LaurentTaylorJet


def _lift_numeric(jet: IntervalTaylorJet, label: str) -> LaurentTaylorJet:
    return LaurentTaylorJet.from_numeric(jet, label=label)


def _formal_divide_by_numeric(value: LaurentTaylorJet, denominator: IntervalTaylorJet, label: str) -> LaurentTaylorJet:
    _require_nonzero(denominator[0], label)
    return value * _lift_numeric(denominator.reciprocal(), f"reciprocal({label})")


def _source_jet(
    context: Any,
    grid: Any,
    radial_order: int,
    axial_order: int,
    label: str,
) -> LaurentTaylorJet:
    coefficients = []
    for order in range(axial_order + 1):
        value = _source_value(grid, radial_order, order, label)
        terms = _term_list(
            context,
            value,
            label=label,
            radial_order=radial_order,
            axial_order=order,
        )
        coefficients.append(
            [
                {
                    "powers": term["powers"],
                    "coefficient": _factorial_div(term["coefficient"], order, context),
                    "origin": _origin(
                        "factorial_normalize_source_derivative",
                        term["origin"],
                        factorial=math.factorial(order),
                    ),
                }
                for term in terms
                if not (term['coefficient'].a == 0 and term['coefficient'].b == 0)
            ]
        )
    return LaurentTaylorJet(context, coefficients)


def _axial_second(
    field: LaurentTaylorJet,
    field_r: LaurentTaylorJet,
    field_rr: LaurentTaylorJet,
    z: LaurentTaylorJet,
    rho: LaurentTaylorJet,
    delta: LaurentTaylorJet,
    weight: LaurentTaylorJet,
    L_numeric: IntervalTaylorJet,
) -> LaurentTaylorJet:
    """The inner-Volterra axial second operator in Laurent Taylor form."""

    one = LaurentTaylorJet.constant(z.context, 1, z.order, origin="one")
    d = one - z * z
    field_z = field.differentiate()
    field_r_z = field_r.differentiate()
    field_zz = field_z.differentiate()
    A = weight * z * field + d * field_z - 2 * z * rho * field_r
    A_Z = (
        weight * field
        + (weight - 2) * z * field_z
        + d * field_zz
        - 2 * rho * field_r
        - 2 * z * rho * field_r_z
    )
    A_R = (weight - 2) * z * field_r + d * field_r_z - 2 * z * rho * field_rr
    first = _formal_divide_by_numeric((weight - 1 + delta) * z * A, L_numeric, "L")
    second_inside = _formal_divide_by_numeric(A_Z, L_numeric, "L")
    second_inside = second_inside + _formal_divide_by_numeric(2 * delta * z * A, L_numeric * L_numeric, "L^2")
    second = d * second_inside
    third = _formal_divide_by_numeric(2 * z * rho * A_R, L_numeric, "L")
    return _formal_divide_by_numeric(first + second - third, L_numeric, "L")


def _relative_f0_jet(
    context: Any,
    z: IntervalTaylorJet,
    delta: IntervalTaylorJet,
    Lambda: IntervalTaylorJet,
    required_j: Any,
    order: int,
) -> IntervalTaylorJet:
    """Build ``F0(Z+h)/F0(Z)`` from the local same-source log derivative."""

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
    result = []
    for index, term in enumerate(jet[order]):
        result.append(
            {
                "powers": tuple(term["powers"]),
                "coefficient": term["coefficient"] * multiplier,
                "origin": _origin("ordinary_Z_derivative", term["origin"], Z_order=order, term_index=index, family=label),
            }
        )
    return result


def _shift_powers(
    jet: LaurentTaylorJet,
    shift: Sequence[int],
    *,
    operation: str,
    family: str,
) -> LaurentTaylorJet:
    if len(shift) != 4:
        raise ValueError("a Laurent power shift needs four entries")
    shift_tuple = tuple(int(value) for value in shift)
    return LaurentTaylorJet(
        jet.context,
        [
            [
                {
                    "powers": tuple(a + b for a, b in zip(term["powers"], shift_tuple)),
                    "coefficient": term["coefficient"],
                    "origin": _origin(operation, term["origin"], family=family, power_shift=shift_tuple),
                }
                for term in row
            ]
            for row in jet.coefficients
        ],
    )


def build_first_collar_laurent_forcing_jets(
    c: Any,
    rho: Any,
    Z: Any,
    delta: Any,
    Lambda: Any,
    required_j: Any,
    *positional: Any,
    Phi_derivatives: Any = None,
    Uz_derivatives: Any = None,
    radial_average_Uz_derivatives: Any = None,
    log_F0: Any = None,
    x: Any = None,
    F0_log: Any = None,
    logF0: Any = None,
    F0_log_interval: Any = None,
    mean_derivatives: Any = None,
    whole_core: Any = None,
    source: Any = None,
    hb: Any = None,
    Pstar_squared: Any = None,
    source_power_basis: Any = None,
    **unknown: Any,
) -> dict[str, Any]:
    """Return unmerged Laurent sectors for ``g4/g5/g6`` ordinary ``Z0..Z3``.

    Required source rows are ``rho^i Z^k`` with ``i=0,1,2`` and
    ``0 <= k <= 5-i`` for ``Phi`` and ``U``; the mean uses ``i=0,1`` with the
    same axial limits.  A missing row is an immediate error.  Source terms
    are never evaluated at ``hb``, ``F0_anchor`` or ``Pstar_squared`` values.
    """

    if unknown:
        names = ", ".join(sorted(unknown))
        raise TypeError(f"unknown first-collar Laurent forcing input(s): {names}")
    if not hasattr(c, "mpf"):
        raise TypeError("c must be the caller's MPIntervalContext")
    if source_power_basis is not None and tuple(source_power_basis) != SOURCE_POWER_BASIS:
        raise ValueError(f"source_power_basis must be {SOURCE_POWER_BASIS!r}")

    x, Phi_derivatives, Uz_derivatives, radial_average_Uz_derivatives, log_F0 = _parse_positional_inputs(
        positional,
        x=x,
        Phi_derivatives=Phi_derivatives,
        Uz_derivatives=Uz_derivatives,
        radial_average_Uz_derivatives=radial_average_Uz_derivatives,
        log_F0=log_F0,
    )
    containers = [value for value in (whole_core, source) if isinstance(value, Mapping)]
    if x is None:
        for container in containers:
            supplied_x = _pick(container, ("x", "radial_x", "sqrt_rho"))
            if supplied_x is not None:
                x = supplied_x
                break
    if x is None:
        raise TypeError("caller-supplied x is required; x is never inferred from rho")

    phi_grid, uz_grid, mean_grid, log_f0 = _resolve_inputs(
        Phi_derivatives=Phi_derivatives,
        Uz_derivatives=Uz_derivatives,
        radial_average_Uz_derivatives=radial_average_Uz_derivatives,
        mean_derivatives=mean_derivatives,
        log_F0=log_F0,
        F0_log=F0_log,
        logF0=logF0,
        F0_log_interval=F0_log_interval,
        whole_core=whole_core,
        source=source,
    )

    rho_iv = _interval(c, rho, "rho")
    z_iv = _interval(c, Z, "Z")
    delta_iv = _interval(c, delta, "delta")
    lambda_iv = _interval(c, Lambda, "Lambda")
    j_iv = _interval(c, required_j, "required_j")
    x_iv = _interval(c, x, "x")
    log_f0_iv = _interval(c, log_f0, "log F0")
    _require_domain(rho_iv, _RHO_MIN, _RHO_MAX, c, "rho")
    _require_domain(z_iv, _Z_MIN, _Z_MAX, c, "Z")
    _require_positive(lambda_iv, "Lambda")
    _require_positive(j_iv, "required_j")
    _require_positive(x_iv, "x")
    x_squared = x_iv ** 2
    if x_squared.a > rho_iv.a or x_squared.b < rho_iv.b:
        raise ValueError("x must enclose the positive square root of the whole rho interval")

    degree = _MAX_SOURCE_ORDER
    rho_num = IntervalTaylorJet.constant(c, rho_iv, degree)
    z_num = IntervalTaylorJet.variable(c, z_iv, degree)
    delta_num = IntervalTaylorJet.constant(c, delta_iv, degree)
    lambda_num = IntervalTaylorJet.constant(c, lambda_iv, degree)
    one_num = IntervalTaylorJet.constant(c, 1, degree)
    d_num = one_num - z_num * z_num
    L_num = one_num - delta_num * z_num * z_num
    _require_nonzero(L_num[0], "L=1-delta*Z^2")

    rho = _lift_numeric(rho_num, "rho")
    z = _lift_numeric(z_num, "Z")
    delta = _lift_numeric(delta_num, "delta")
    one = _lift_numeric(one_num, "one")
    d = _lift_numeric(d_num, "1-Z^2")
    U0 = _source_jet(c, uz_grid, 0, 5, "Uz_derivatives")
    U1 = _source_jet(c, uz_grid, 1, 4, "Uz_derivatives")
    U2 = _source_jet(c, uz_grid, 2, 3, "Uz_derivatives")
    M0 = _source_jet(c, mean_grid, 0, 5, "radial_average_Uz_derivatives")
    M1 = _source_jet(c, mean_grid, 1, 4, "radial_average_Uz_derivatives")
    Phi0 = _source_jet(c, phi_grid, 0, 5, "Phi_derivatives")
    Phi1 = _source_jet(c, phi_grid, 1, 4, "Phi_derivatives")
    Phi2 = _source_jet(c, phi_grid, 2, 3, "Phi_derivatives")

    M0_Z = M0.differentiate()
    M1_Z = M1.differentiate()
    U1_Z = U1.differentiate()
    Q0 = _formal_divide_by_numeric(2 * z * U0 - (1 - delta) * z * M0 - d * M0_Z, L_num, "L")
    Qr = _formal_divide_by_numeric(2 * z * U1 - (1 - delta) * z * M1 - d * M1_Z, L_num, "L")
    rhoQrr = _formal_divide_by_numeric(
        2 * z * rho * U2
        - (1 - delta) * z * (U1 - 2 * M1)
        - d * (U1_Z - 2 * M1_Z),
        L_num,
        "L",
    )
    Q0_Z = Q0.differentiate()
    S = _formal_divide_by_numeric(Q0 + (1 - delta) * z * Q0_Z / 2 + rho * Qr, L_num, "L")
    S = S + Q0 * (Q0 / 2 + rho * Qr)
    S = S + _formal_divide_by_numeric(U0 * (-2 * z * Q0 + d * Q0_Z - 2 * z * rho * Qr), L_num, "L")
    P = 2 * Qr + rhoQrr

    x_num = IntervalTaylorJet.constant(c, x_iv, degree)
    x_formal = _lift_numeric(x_num, "caller_x")
    g4_power0 = 2 * x_formal * P
    g4_power1 = -x_formal * S

    relative_f0_num = _relative_f0_jet(c, z_num, delta_num, lambda_num, j_iv, degree)
    relative_f0 = _lift_numeric(relative_f0_num, "relative_F0")
    F0_ratio = (Phi0 * relative_f0).truncate(5)
    F1_ratio = (Phi1 * relative_f0).truncate(4)
    F2_ratio = (Phi2 * relative_f0).truncate(3)
    weight_F = _lift_numeric(IntervalTaylorJet.constant(c, -2, degree) - delta_num, "weight_F")
    weight_U = _lift_numeric(IntervalTaylorJet.constant(c, -1, degree) - delta_num, "weight_U")
    Nt_over_F0_anchor = -_axial_second(F0_ratio, F1_ratio, F2_ratio, z, rho, delta, weight_F, L_num)
    Nz = -_axial_second(U0, U1, U2, z, rho, delta, weight_U, L_num)

    g6_power1 = 2 * Nz - _formal_divide_by_numeric(4 * rho * z * P, L_num, "L")
    g6_power2 = _formal_divide_by_numeric(2 * rho * z * S, L_num, "L")

    # The external powers are metadata only.  Internal Lambda use above is
    # numerical interval enclosure for the relative-F0 coefficient recipe.
    families = {
        "4": (
            (g4_power0, (0, 0, 0, 0), "g4_Lambda0"),
            (g4_power1, (0, 0, 0, -1), "g4_Lambda-1"),
        ),
        "5": ((2 * Nt_over_F0_anchor, (0, 1, 0, -1), "g5_F0_anchor_Lambda-1"),),
        "6": (
            (g6_power1, (0, 0, 0, -1), "g6_Lambda-1"),
            (g6_power2, (0, 0, 0, -2), "g6_Lambda-2"),
        ),
    }
    sectors: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for order in _OUTPUT_ORDERS:
        rows: dict[str, list[dict[str, Any]]] = {}
        for family, entries in families.items():
            emitted: list[dict[str, Any]] = []
            for jet, shift, sector_label in entries:
                shifted = _shift_powers(jet, shift, operation="external_formal_power", family=sector_label)
                emitted.extend(_ordinary_terms(shifted, order, f"g{family}"))
            rows[family] = emitted
        sectors[str(order)] = rows

    proof = {
        "kind": "first_collar_forcing_laurent_jets_v1",
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(_OUTPUT_ORDERS),
        "ordinary_Z_derivative_order": 3,
        "source_derivatives_are_ordinary": True,
        "Taylor_coefficients_divided_by_factorial": True,
        "unknown_high_derivatives_padded": False,
        "source_grid_total_order_cap": _MAX_SOURCE_ORDER,
        "required_source_rows": {
            "Phi": "rho^i Z^k for i=0,1,2 and 0<=k<=5-i",
            "U": "rho^i Z^k for i=0,1,2 and 0<=k<=5-i",
            "mean": "rho^i Z^k for i=0,1 and 0<=k<=5-i",
        },
        "rho_domain": "[4,4.1]",
        "Z_domain": "[-1,1]",
        "rho_zero_safe": True,
        "rho_zero_rule": "rho*m_rr = Ur-2*mr; no division by rho",
        "interval_outward_rounding": True,
        "denominator_checks": ["L=1-delta*Z^2 excludes zero", "H^2+sigma^2 excludes zero"],
        "formal_laurent_terms_unmerged": True,
        "hb_negative_powers_preserved": True,
        "input_source_provenance_retained": True,
        "input_scalar_intervals_promoted_to_zero_sector": True,
        "F0_power_is_formal": True,
        "F0_current_Z_factor_formal": True,
        "F0_log_interval_supplied": True,
        "F0_log_interval": log_f0_iv,
        "F0_upper_not_used_as_value": True,
        "F0_absolute_exp_materialized": False,
        "Pstar_squared_materialized": False,
        "hb_materialized": False,
        "no_input_midpoint_or_error_bound_chosen_by_leaf": True,
        "no_sector_combination": True,
        "x_caller_supplied": True,
        "x_inferred_from_rho": False,
        "positive_x_squared_encloses_rho": True,
        "g5_anchor_rule": "F(Z+h)/F0(Z) jet at each current Z; no extra D_F",
        "Lambda_inside_relative_F0_derivative_coefficients_is_enclosed_numerically": True,
        "outer_Lambda_power_remains_formal": True,
        "lambda_role_source": (
            "caller Lambda interval is used in the relative F0 jet; "
            "row Lambda powers are formal and are not numerically materialized"
        ),
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
        "hb_argument_supplied": hb is not None,
        "Pstar_squared_argument_supplied": Pstar_squared is not None,
    }
    jets = {
        "Q0": Q0,
        "Qr": Qr,
        "rhoQrr": rhoQrr,
        "S": S,
        "P": P,
        "Np_formal": -S / 2 + _lift_numeric(lambda_num, "internal_Lambda") * P,
        "Nt_over_F0_anchor": Nt_over_F0_anchor,
        "Nz": Nz,
        "relative_F0": relative_f0,
        "x": x_formal,
    }
    by_family = {
        f"g{family}": {
            f"Z{order}": sectors[str(order)][family]
            for order in _OUTPUT_ORDERS
        }
        for family in families
    }
    flat_rows = {
        f"g{family}_Z{order}": sectors[str(order)][family]
        for family in families
        for order in _OUTPUT_ORDERS
    }
    result: dict[str, Any] = {
        "sectors": sectors,
        "rows": sectors,
        "by_order": sectors,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "proof": proof,
        "metadata": proof,
        "jets": jets,
        "log_F0": log_f0_iv,
        "F0_log": log_f0_iv,
        "F0_log_interval": log_f0_iv,
        **by_family,
        **flat_rows,
    }
    # Direct order keys make the order-major callback convenient while the
    # named ``sectors`` mapping remains the canonical result field.
    result.update(sectors)
    return result


def first_collar_laurent_forcing_jets(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Alias for :func:`build_first_collar_laurent_forcing_jets`."""

    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def derive_first_collar_laurent_forcing_jets(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compatibility alias for the Laurent source leaf."""

    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def build_first_collar_laurent_forcing(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Short name for the formal-sector first-collar forcing evaluator."""

    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def first_collar_laurent_forcing(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Short alias for :func:`build_first_collar_laurent_forcing_jets`."""

    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def build_collar_laurent_forcing(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compatibility alias used by collar source dispatchers."""

    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def build_first_collar_laurent_forcing_z3(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compatibility alias retaining the ordinary ``z3`` operation name."""

    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def first_collar_laurent_forcing_z3(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Short ``z3`` alias for the Laurent source leaf."""

    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def build_first_collar_forcing_laurent_jets(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Word-order compatibility alias for downstream source dispatchers."""

    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def first_collar_forcing_laurent_jets(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Word-order compatibility alias for downstream source dispatchers."""

    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def build_first_collar_forcing_laurent(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compact compatibility alias for Laurent forcing dispatchers."""

    return build_first_collar_laurent_forcing_jets(*args, **kwargs)


def first_collar_forcing_laurent(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compact compatibility alias for Laurent forcing dispatchers."""

    return build_first_collar_laurent_forcing_jets(*args, **kwargs)

