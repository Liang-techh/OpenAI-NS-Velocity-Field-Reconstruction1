"""Conditional first-collar radial Laurent source leaf.

This module is deliberately a small, dependency-free adapter.  It does not
import a producer, read a saved packet, or select a point from an interval.
The caller supplies the actual zeroth-radial ``Phi``, ``U`` and own mean
histories, the comparison-core ``Phi``/``Psi`` radial grid, and the ordinary
``Z`` derivatives of the original cutoff.  The adapter then transcribes the
radial identities used by the original collar fields and by ``first_rows``.

The returned entries are signed source sectors with the common shape::

    {"powers": (hb, F0_anchor, Pstar_squared, Lambda),
     "coefficient": directed_mpmath_interval,
     "origin": provenance}

``powers[0]`` is allowed to be positive or negative.  In particular, the
``chi_r`` identity is represented as a power-zero term minus a power-minus-one
term; the inverse width is never evaluated.  Source sectors are not merged,
so caller provenance survives products and sums.

This is a conditional source leaf, not an assertion that the comparison core
is the actual collar.  The caller must prove the supplied ``rho`` and
same-source identities.  Those contracts are recorded, not independently
verified here.
"""

from __future__ import annotations

import copy
import math
from collections.abc import Mapping, Sequence
from typing import Any


__all__ = [
    "IntervalTaylorJet",
    "TaylorJet",
    "SectorJet",
    "build",
    "build_collar_laurent_radial",
    "collar_laurent_radial",
    "derive_collar_laurent_radial",
    "build_radial",
]


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
_RADIAL_ORDERS = (0, 1, 2)
_TOTAL_Z_ORDER = 5


def _is_interval(value: Any) -> bool:
    """Return whether *value* has the mpmath interval shape."""

    return hasattr(value, "_mpi_") or (hasattr(value, "a") and hasattr(value, "b"))


def _bounds(value: Any) -> tuple[Any, Any]:
    if not _is_interval(value):
        raise TypeError("an mpmath interval is required")
    try:
        return value.a, value.b
    except AttributeError:  # pragma: no cover - defensive for alternate iv APIs
        rows = value._mpi_
        return rows[0], rows[1]


def _interval(context: Any, value: Any, label: str = "interval") -> Any:
    """Coerce an interval without choosing a midpoint or endpoint."""

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


def _is_term_mapping(value: Any) -> bool:
    return isinstance(value, Mapping) and (
        "powers" in value or "coefficient" in value or "coefficient_interval" in value
    )


def _is_term_list(value: Any) -> bool:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        return False
    return bool(value) and all(_is_term_mapping(item) for item in value)


def _mapping_get(mapping: Mapping[Any, Any], names: Sequence[Any]) -> Any | None:
    for name in names:
        try:
            if name in mapping:
                return mapping[name]
        except (TypeError, ValueError):
            continue
    return None


def _source_value(grid: Any, radial_order: int, axial_order: int, label: str) -> Any:
    """Read one grid cell, rejecting missing high ordinary derivatives."""

    if grid is None:
        raise TypeError(f"{label} source grid is required")

    # Live jet objects are accepted as a convenience for callers that already
    # normalized a single radial row.  Their cells are still read explicitly;
    # no absent derivative is manufactured here.
    if isinstance(grid, (SectorJet, IntervalTaylorJet)):
        if radial_order != 0 or axial_order > grid.order:
            raise KeyError(f"missing {label} ordinary derivative rho^{radial_order} Z^{axial_order}")
        if isinstance(grid, SectorJet):
            return list(grid[axial_order])
        return grid[axial_order]

    # A wrapper is convenient for callers that keep a named derivative row.
    if isinstance(grid, Mapping):
        nested_grid = _mapping_get(grid, ("grid", "derivatives", "rows", "values"))
        if nested_grid is not None:
            return _source_value(nested_grid, radial_order, axial_order, label)

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

        # A one-dimensional row may be keyed by Z order.
        if radial_order == 0:
            direct_z = _mapping_get(
                grid,
                (axial_order, str(axial_order), f"Z{axial_order}", f"z{axial_order}"),
            )
            if direct_z is not None:
                return direct_z

        row = _mapping_get(
            grid,
            (
                radial_order,
                str(radial_order),
                f"rho{radial_order}",
                f"rho_{radial_order}",
                f"radial_{radial_order}",
            ),
        )
        if row is not None:
            if isinstance(row, Mapping):
                nested = _mapping_get(
                    row,
                    (axial_order, str(axial_order), f"Z{axial_order}", f"z{axial_order}"),
                )
                if nested is not None:
                    return nested
            elif isinstance(row, Sequence) and not isinstance(row, (str, bytes, bytearray)):
                if axial_order < len(row):
                    return row[axial_order]

    if isinstance(grid, Sequence) and not isinstance(grid, (str, bytes, bytearray)):
        # A plain row is the natural representation for the caller-supplied
        # actual fields.  Sector cells are sequences themselves, so detect
        # them before treating the outer sequence as radial rows.
        one_dimensional = bool(grid) and (
            _is_interval(grid[0])
            or _is_term_mapping(grid[0])
            or _is_term_list(grid[0])
            or not isinstance(grid[0], Sequence)
            or isinstance(grid[0], (str, bytes, bytearray))
        )
        if one_dimensional:
            if radial_order != 0:
                raise KeyError(
                    f"missing {label} ordinary derivative rho^{radial_order} Z^{axial_order}"
                )
            if axial_order < len(grid):
                return grid[axial_order]
        elif radial_order < len(grid):
            row = grid[radial_order]
            if isinstance(row, Sequence) and not isinstance(row, (str, bytes, bytearray)):
                if axial_order < len(row):
                    return row[axial_order]

    raise KeyError(f"missing {label} ordinary derivative rho^{radial_order} Z^{axial_order}")


def _term_powers(raw: Mapping[Any, Any], label: str) -> tuple[int, int, int, int]:
    powers = raw.get("powers")
    if powers is None:
        # Compatibility with the earlier first-collar sector spelling.  The
        # canonical output always uses one four-entry tuple.
        powers = (
            raw.get("hb_power", 0),
            raw.get("F0_power", 0),
            raw.get("Pstar_squared_power", raw.get("Pstar_power", 0)),
            raw.get("Lambda_power", 0),
        )
    if not isinstance(powers, (tuple, list)) or len(powers) != 4:
        raise ValueError(f"{label} sector needs four source powers")
    if any(not isinstance(power, int) or isinstance(power, bool) for power in powers):
        raise TypeError(f"{label} source powers must be integers")
    return tuple(int(power) for power in powers)


def _sector_terms(
    context: Any,
    value: Any,
    label: str,
    *,
    source_family: Any,
    default_origin: Any,
) -> list[dict[str, Any]]:
    """Normalize one source cell to canonical, unmerged sector records."""

    if isinstance(value, Mapping) and not _is_term_mapping(value):
        nested = _mapping_get(value, ("sectors", "sector_terms", "terms"))
        if nested is not None:
            return _sector_terms(
                context,
                nested,
                label,
                source_family=source_family,
                default_origin=default_origin,
            )

    if _is_term_mapping(value):
        raw_terms: Sequence[Any] = (value,)
    elif _is_term_list(value):
        raw_terms = value
    else:
        coefficient = _interval(context, value, label)
        return [
            {
                "powers": (0, 0, 0, 0),
                "coefficient": coefficient,
                "origin": ("caller_scalar", source_family, default_origin),
            }
        ]

    result: list[dict[str, Any]] = []
    for index, raw in enumerate(raw_terms):
        powers = _term_powers(raw, f"{label}[{index}]")
        coefficient_raw = raw.get("coefficient", raw.get("coefficient_interval"))
        if coefficient_raw is None:
            raise ValueError(f"{label}[{index}] sector coefficient is missing")
        coefficient = _interval(context, coefficient_raw, f"{label}[{index}] coefficient")
        origin = ("caller_sector", raw.get("origin", default_origin), raw.get("semantics", "unclassified"), source_family, index)
        result.append({"powers": powers, "coefficient": coefficient, "origin": origin})
    return result


class IntervalTaylorJet:
    """Factorial-normalized ordinary-Z Taylor jet over directed intervals."""

    __slots__ = ("context", "coefficients")

    def __init__(self, context: Any, coefficients: Sequence[Any]):
        values = tuple(coefficients)
        if not values:
            raise ValueError("a Taylor jet needs at least one coefficient")
        self.context = context
        self.coefficients = values

    @property
    def order(self) -> int:
        return len(self.coefficients) - 1

    @classmethod
    def constant(cls, context: Any, value: Any, order: int) -> "IntervalTaylorJet":
        if order < 0:
            raise ValueError("nonnegative Taylor order required")
        value = _interval(context, value, "Taylor constant")
        zero = context.mpf(0)
        return cls(context, [value] + [zero] * order)

    def truncate(self, order: int) -> "IntervalTaylorJet":
        if order < 0 or order > self.order:
            raise ValueError("requested Taylor order is unavailable")
        return type(self)(self.context, self.coefficients[: order + 1])

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
        return type(self)(
            self.context,
            [self[index] + other[index] for index in range(order + 1)],
        )

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


TaylorJet = IntervalTaylorJet


class SectorJet:
    """Taylor jet whose coefficients are unmerged formal source sectors."""

    __slots__ = ("context", "coefficients")

    def __init__(self, context: Any, coefficients: Sequence[Sequence[Mapping[str, Any]]]):
        self.context = context
        # Discard only exactly-zero newly generated terms; interval errors and
        # all nonzero source sectors remain separate and retain provenance.
        self.coefficients = tuple(tuple(dict(term) for term in row if not (term["coefficient"].a == 0 and term["coefficient"].b == 0)) for row in coefficients)
        if not self.coefficients:
            raise ValueError("a sector jet needs at least one coefficient")
        for row in self.coefficients:
            for term in row:
                _validate_sector(term, context)

    @property
    def order(self) -> int:
        return len(self.coefficients) - 1

    @classmethod
    def from_cells(
        cls,
        context: Any,
        cells: Sequence[Sequence[Mapping[str, Any]]],
    ) -> "SectorJet":
        return cls(context, cells)

    @classmethod
    def constant(
        cls,
        context: Any,
        value: Any,
        order: int,
        *,
        powers: tuple[int, int, int, int] = (0, 0, 0, 0),
        origin: Any = ("constant",),
    ) -> "SectorJet":
        coefficient = _interval(context, value, "sector constant")
        zero = context.mpf(0)
        rows = [[{"powers": powers, "coefficient": coefficient, "origin": origin}]]
        rows.extend(
            [{"powers": powers, "coefficient": zero, "origin": ("constant_zero", origin)}]
            for _ in range(order)
        )
        return cls(context, rows)

    def truncate(self, order: int) -> "SectorJet":
        if order < 0 or order > self.order:
            raise ValueError("requested sector Taylor order is unavailable")
        return type(self)(self.context, self.coefficients[: order + 1])

    def __getitem__(self, index: int) -> tuple[dict[str, Any], ...]:
        return self.coefficients[index]

    def _coerce(self, other: Any) -> "SectorJet":
        if isinstance(other, SectorJet):
            if other.context is not self.context:
                raise TypeError("sector jets must use the caller's one interval context")
            return other
        if isinstance(other, IntervalTaylorJet):
            if other.context is not self.context:
                raise TypeError("sector and Taylor jets must share one interval context")
            rows = []
            for index in range(other.order + 1):
                rows.append(
                    [
                        {
                            "powers": (0, 0, 0, 0),
                            "coefficient": other[index],
                            "origin": ("numeric_jet", index),
                        }
                    ]
                )
            return type(self)(self.context, rows)
        return type(self).constant(self.context, other, self.order, origin=("scalar",))

    def __add__(self, other: Any) -> "SectorJet":
        other = self._coerce(other)
        order = min(self.order, other.order)
        # Deliberately concatenate rather than combine equal powers.  This is
        # what keeps the original signed source provenance visible.
        return type(self)(
            self.context,
            [list(self[index]) + list(other[index]) for index in range(order + 1)],
        )

    __radd__ = __add__

    def __neg__(self) -> "SectorJet":
        return type(self)(
            self.context,
            [
                [
                    {
                        "powers": term["powers"],
                        "coefficient": -term["coefficient"],
                        "origin": ("negate", term["origin"]),
                    }
                    for term in row
                ]
                for row in self.coefficients
            ],
        )

    def __sub__(self, other: Any) -> "SectorJet":
        return self + (-self._coerce(other))

    def __rsub__(self, other: Any) -> "SectorJet":
        return self._coerce(other) - self

    def __mul__(self, other: Any) -> "SectorJet":
        other = self._coerce(other)
        order = min(self.order, other.order)
        rows: list[list[dict[str, Any]]] = []
        for index in range(order + 1):
            terms: list[dict[str, Any]] = []
            for left_index in range(index + 1):
                right_index = index - left_index
                for left in self[left_index]:
                    for right in other[right_index]:
                        terms.append(
                            {
                                "powers": tuple(
                                    a + b
                                    for a, b in zip(left["powers"], right["powers"])
                                ),
                                "coefficient": left["coefficient"] * right["coefficient"],
                                "origin": (
                                    "product",
                                    left["origin"],
                                    right["origin"],
                                ),
                            }
                        )
            rows.append(terms)
        return type(self)(self.context, rows)

    __rmul__ = __mul__

    def __truediv__(self, other: Any) -> "SectorJet":
        if isinstance(other, SectorJet):
            raise ValueError("Laurent sector denominators are not inverted")
        return self * self._coerce(other).to_numeric().reciprocal()

    def shift_hb(self, amount: int, *, origin: Any = "hb_shift") -> "SectorJet":
        if not isinstance(amount, int) or isinstance(amount, bool):
            raise TypeError("hb shift must be an integer")
        return type(self)(
            self.context,
            [
                [
                    {
                        "powers": (
                            term["powers"][0] + amount,
                            *term["powers"][1:],
                        ),
                        "coefficient": term["coefficient"],
                        "origin": (origin, amount, term["origin"]),
                    }
                    for term in row
                ]
                for row in self.coefficients
            ],
        )

    def to_numeric(self) -> IntervalTaylorJet:
        """Collapse only for a width-free denominator, with a strict check."""

        rows: list[Any] = []
        for index, terms in enumerate(self.coefficients):
            if any(term["powers"] != (0, 0, 0, 0) for term in terms):
                raise ValueError("a formal sector cannot be used as a denominator")
            value = self.context.mpf(0)
            for term in terms:
                value += term["coefficient"]
            rows.append(value)
        return IntervalTaylorJet(self.context, rows)


def _validate_sector(term: Mapping[str, Any], context: Any) -> None:
    if set(("powers", "coefficient", "origin")) - set(term):
        raise ValueError("sector metadata missing powers, coefficient, or origin")
    powers = term["powers"]
    if not isinstance(powers, tuple) or len(powers) != 4:
        raise ValueError("sector powers must be a four-entry tuple")
    if any(not isinstance(power, int) or isinstance(power, bool) for power in powers):
        raise TypeError("sector powers must be integers")
    if not _is_interval(term["coefficient"]):
        raise TypeError("sector coefficient must be a directed interval")
    # This check is intentionally identity based: a source coefficient built
    # in another interval context must not silently cross the caller boundary.
    try:
        _ = term["coefficient"].a + context.mpf(0)
    except (AttributeError, TypeError, ValueError):
        raise TypeError("sector coefficient is not compatible with caller context")


def _ordinary_sector_rows(jet: SectorJet, max_order: int) -> dict[tuple[int, int], list[dict[str, Any]]]:
    rows: dict[tuple[int, int], list[dict[str, Any]]] = {}
    for order in range(max_order + 1):
        factor = math.factorial(order)
        rows[(0, order)] = [
            {
                "powers": term["powers"],
                "coefficient": term["coefficient"] * factor,
                "origin": term["origin"],
            }
            for term in jet[order]
        ]
    return rows


def _ordinary_sector_grid(jets: Sequence[SectorJet]) -> dict[tuple[int, int], list[dict[str, Any]]]:
    rows: dict[tuple[int, int], list[dict[str, Any]]] = {}
    for radial_order, jet in enumerate(jets):
        for axial_order in range(jet.order + 1):
            factor = math.factorial(axial_order)
            rows[(radial_order, axial_order)] = [
                {
                    "powers": term["powers"],
                    "coefficient": term["coefficient"] * factor,
                    "origin": term["origin"],
                }
                for term in jet[axial_order]
            ]
    return rows


def _numeric_source_jet(
    context: Any,
    grid: Any,
    radial_order: int,
    axial_order: int,
    label: str,
) -> IntervalTaylorJet:
    values = []
    for order in range(axial_order + 1):
        raw = _source_value(grid, radial_order, order, label)
        # Comparison rows are ordinary interval data, not formal actual
        # sectors.  A zero-power one-term record is accepted as a convenience
        # but its role remains comparison-only in the output metadata.
        terms = _sector_terms(
            context,
            raw,
            f"{label} rho^{radial_order} Z^{order}",
            source_family="comparison",
            default_origin=(label, radial_order, order),
        )
        if any(term["powers"] != (0, 0, 0, 0) for term in terms):
            raise ValueError(f"{label} comparison source must be width-free")
        value = context.mpf(0)
        for term in terms:
            value += term["coefficient"]
        values.append(value / math.factorial(order))
    return IntervalTaylorJet(context, values)


def _sector_source_jet(
    context: Any,
    grid: Any,
    radial_order: int,
    axial_order: int,
    label: str,
    source_family: Any,
) -> SectorJet:
    rows = []
    for order in range(axial_order + 1):
        raw = _source_value(grid, radial_order, order, label)
        terms = _sector_terms(
            context,
            raw,
            f"{label} rho^{radial_order} Z^{order}",
            source_family=source_family,
            default_origin=(label, radial_order, order),
        )
        rows.append(
            [
                {
                    "powers": term["powers"],
                    "coefficient": term["coefficient"] / math.factorial(order),
                    "origin": term["origin"],
                }
                for term in terms
            ]
        )
    return SectorJet(context, rows)


def _numeric_row_jet(
    context: Any,
    value: Any,
    max_order: int,
    label: str,
) -> IntervalTaylorJet:
    if isinstance(value, IntervalTaylorJet):
        if value.context is not context:
            raise TypeError(f"{label} jet uses a different interval context")
        if value.order < max_order:
            raise KeyError(f"missing {label} ordinary Z^{max_order} derivative")
        return value.truncate(max_order)
    values = []
    for order in range(max_order + 1):
        raw = _source_value(value, 0, order, label)
        values.append(_interval(context, raw, f"{label} Z^{order}") / math.factorial(order))
    return IntervalTaylorJet(context, values)


def _sector_row_jet(
    context: Any,
    value: Any,
    max_order: int,
    label: str,
    source_family: Any,
) -> SectorJet:
    if isinstance(value, SectorJet):
        if value.context is not context:
            raise TypeError(f"{label} jet uses a different interval context")
        if value.order < max_order:
            raise KeyError(f"missing {label} ordinary Z^{max_order} derivative")
        return value.truncate(max_order)
    rows = []
    for order in range(max_order + 1):
        raw = _source_value(value, 0, order, label)
        terms = _sector_terms(
            context,
            raw,
            f"{label} Z^{order}",
            source_family=source_family,
            default_origin=(label, 0, order),
        )
        rows.append(
            [
                {
                    "powers": term["powers"],
                    "coefficient": term["coefficient"] / math.factorial(order),
                    "origin": term["origin"],
                }
                for term in terms
            ]
        )
    return SectorJet(context, rows)


def _sector_from_numeric(
    context: Any,
    value: IntervalTaylorJet,
    *,
    origin: Any,
) -> SectorJet:
    return SectorJet(
        context,
        [
            [
                {
                    "powers": (0, 0, 0, 0),
                    "coefficient": coefficient,
                    "origin": (origin, index),
                }
            ]
            for index, coefficient in enumerate(value.coefficients)
        ],
    )


def _numeric_constants(context: Any, value: Any, order: int, label: str) -> IntervalTaylorJet:
    return IntervalTaylorJet.constant(context, _interval(context, value, label), order)


def _check_grid_presence(grid: Any, label: str, *, radial_orders: Sequence[int]) -> None:
    """Force every requested source cell to be present before algebra starts."""

    for radial_order in radial_orders:
        for axial_order in range(_TOTAL_Z_ORDER - radial_order + 1):
            _source_value(grid, radial_order, axial_order, label)


def _resolve_alias(primary: Any, aliases: Sequence[tuple[str, Any]], label: str) -> Any:
    if primary is not None:
        return primary
    for alias, value in aliases:
        if value is not None:
            return value
    raise TypeError(f"caller-supplied {label} is required")


def build_collar_laurent_radial(
    c: Any,
    rho: Any,
    Z: Any,
    log_hb: Any,
    actual_Phi: Any = None,
    actual_U: Any = None,
    own_mean: Any = None,
    comparison_Phi: Any = None,
    comparison_U: Any = None,
    sigma: Any = None,
    sigma_prime: Any = None,
    *,
    actual_phi: Any = None,
    actual_u: Any = None,
    actual_mean: Any = None,
    actual_mean_derivatives: Any = None,
    Phi: Any = None,
    U: Any = None,
    mean: Any = None,
    phi: Any = None,
    u: Any = None,
    M: Any = None,
    Phi_derivatives: Any = None,
    U_derivatives: Any = None,
    mean_derivatives: Any = None,
    core_Phi: Any = None,
    core_U: Any = None,
    core_phi: Any = None,
    core_u: Any = None,
    comparison_core_Phi: Any = None,
    comparison_core_U: Any = None,
    comparison_phi: Any = None,
    comparison_u: Any = None,
    chi: Any = None,
    chi_r: Any = None,
    source_family: Any = None,
    source_atoms: Mapping[str, Any] | None = None,
    geometry_certificate: Mapping[str, Any] | None = None,
    source_provenance: Any = None,
    **unknown: Any,
) -> dict[str, Any]:
    """Build first-collar radial ``Phi/U/mean`` sector jets from actual inputs.

    The comparison core U derivatives already include epsilon. Their source
    identity and all inherited interval errors remain caller responsibilities.

    Inputs are ordinary derivative grids.  The actual fields need only be
    supplied at radial order zero (``Z0`` through ``Z5``); radial orders one
    and two are derived from the caller's comparison-core radial grid.  For a
    derivative of radial order ``i`` the output contains ordinary ``Z`` rows
    through ``Z(5-i)``.  A source cell may be a scalar directed interval or a
    list of canonical sector records.

    ``sigma`` and ``sigma_prime`` are the ordinary ``Z`` grids of the original
    cutoff and its phase derivative.  If ``chi``/``chi_r`` are supplied
    instead, they must already be sector-factored grids; this escape hatch is
    useful when a caller has an independently retained cutoff proof.
    """

    if unknown:
        names = ", ".join(sorted(unknown))
        raise TypeError(f"unknown collar Laurent radial input(s): {names}")
    if not hasattr(c, "mpf"):
        raise TypeError("c must be the caller's mpmath interval context")

    actual_Phi = _resolve_alias(
        actual_Phi,
        (("actual_phi", actual_phi), ("Phi", Phi), ("phi", phi), ("Phi_derivatives", Phi_derivatives)),
        "actual Phi",
    )
    actual_U = _resolve_alias(
        actual_U,
        (("actual_u", actual_u), ("U", U), ("u", u), ("U_derivatives", U_derivatives)),
        "actual U",
    )
    own_mean = _resolve_alias(
        own_mean,
        (
            ("actual_mean", actual_mean),
            ("actual_mean_derivatives", actual_mean_derivatives),
            ("mean", mean),
            ("M", M),
            ("mean_derivatives", mean_derivatives),
        ),
        "own mean",
    )
    comparison_Phi = _resolve_alias(
        comparison_Phi,
        (
            ("core_Phi", core_Phi),
            ("core_phi", core_phi),
            ("comparison_core_Phi", comparison_core_Phi),
            ("comparison_phi", comparison_phi),
        ),
        "comparison core Phi",
    )
    comparison_U = _resolve_alias(
        comparison_U,
        (
            ("core_U", core_U),
            ("core_u", core_u),
            ("comparison_core_U", comparison_core_U),
            ("comparison_u", comparison_u),
        ),
        "comparison core Psi",
    )

    if source_family is None or source_family == "":
        raise ValueError("source_family is required and must identify the caller's source")
    if (
        geometry_certificate is None
        or not isinstance(geometry_certificate, Mapping)
        or not geometry_certificate
    ):
        raise ValueError("caller geometry/rho contract is required as geometry_certificate")
    required_geometry = ("hb_strictly_between_zero_and_half", "sc_phase_Z_independent", "hb_epsilon_Z_independent", "sigma_Z_independent")
    if any(geometry_certificate.get(key) is not True for key in required_geometry):
        raise ValueError("caller geometry and fixed-Z cutoff certificates required")
    if source_atoms is None:
        raise ValueError("four original caller source atom records are required")
    else:
        if not isinstance(source_atoms, Mapping) or set(source_atoms) != set(SOURCE_POWER_BASIS):
            raise ValueError(
                "source_atoms must identify exactly hb, F0_anchor, Pstar_squared, and Lambda"
            )
        source_atoms = dict(source_atoms)
        source_atoms_supplied = True

    rho_iv = _interval(c, rho, "rho")
    log_hb_iv = _interval(c, log_hb, "log hb")
    z_iv = _interval(c, Z, "Z")
    _require_positive(rho_iv, "rho")
    if log_hb_iv.b >= c.ln(c.mpf("0.5")).a or rho_iv.a < 4 or rho_iv.b > c.mpf("4.1") or z_iv.a < -1 or z_iv.b > 1:
        raise ValueError("original collar rho/Z/hb domain required")

    # Validate all ordinary source cells before performing any algebra.  A
    # missing comparison Z5, for example, must be reported rather than
    # silently padded by a zero Taylor coefficient.
    _check_grid_presence(actual_Phi, "actual_Phi", radial_orders=(0,))
    _check_grid_presence(actual_U, "actual_U", radial_orders=(0,))
    _check_grid_presence(own_mean, "own_mean", radial_orders=(0,))
    _check_grid_presence(comparison_Phi, "comparison_Phi", radial_orders=_RADIAL_ORDERS)
    _check_grid_presence(comparison_U, "comparison_U", radial_orders=_RADIAL_ORDERS)

    # Actual zeroth-radial fields are source sectors.  The source-cell parser
    # accepts scalar intervals for simple callers but never turns a missing
    # derivative into zero.
    phi0 = _sector_source_jet(c, actual_Phi, 0, _TOTAL_Z_ORDER, "actual_Phi", source_family)
    u0 = _sector_source_jet(c, actual_U, 0, _TOTAL_Z_ORDER, "actual_U", source_family)
    m0 = _sector_source_jet(c, own_mean, 0, _TOTAL_Z_ORDER, "own_mean", source_family)

    pc = [
        _numeric_source_jet(c, comparison_Phi, radial_order, _TOTAL_Z_ORDER - radial_order, "comparison_Phi")
        for radial_order in _RADIAL_ORDERS
    ]
    ps = [
        _sector_source_jet(c, comparison_U, radial_order, _TOTAL_Z_ORDER - radial_order, "comparison_U", source_family)
        for radial_order in _RADIAL_ORDERS
    ]

    # Comparison radial direction.  The conversion to sectors happens only
    # after the ratio has been formed, so it cannot be mistaken for an actual
    # Phi or Psi field in the returned packet.
    A = pc[1] / pc[0]
    Ar = pc[2] / pc[0] - A * A
    A_sector = _sector_from_numeric(c, A, origin=("comparison_core_Phi_ratio", source_family))
    Ar_sector = _sector_from_numeric(c, Ar, origin=("comparison_core_Phi_radial_curvature", source_family))
    psi1_sector = ps[1]
    psi2_sector = ps[2]

    if chi is not None and sigma is not None:
        raise ValueError("supply either sigma or an independently factored chi, not both")
    if chi_r is not None and sigma_prime is not None:
        raise ValueError("supply either sigma_prime or an independently factored chi_r, not both")
    if (chi is None) != (chi_r is None) and (sigma is None or sigma_prime is None):
        raise ValueError("chi and chi_r must be supplied together when sigma data are absent")

    if chi is None:
        if sigma is None or sigma_prime is None:
            raise TypeError("sigma and sigma_prime ordinary Z grids are required")
        sigma_j = _numeric_row_jet(c, sigma, _TOTAL_Z_ORDER, "sigma")
        sigma_p_j = _numeric_row_jet(c, sigma_prime, _TOTAL_Z_ORDER, "sigma_prime")
        for k in range(1, _TOTAL_Z_ORDER + 1):
            if sigma_j[k].a != 0 or sigma_j[k].b != 0 or sigma_p_j[k].a != 0 or sigma_p_j[k].b != 0:
                raise ValueError("cutoff jets must be Z-independent at fixed phase")
        sigma_s = _sector_from_numeric(c, sigma_j, origin=("original_sigma", source_family))
        sigma_p_s = _sector_from_numeric(c, sigma_p_j, origin=("original_sigma_prime", source_family))
        one = SectorJet.constant(c, 1, _TOTAL_Z_ORDER, origin=("unit",))
        chi_s = one - sigma_s + sigma_s.shift_hb(1, origin="chi_hb_positive")
        rho_inv = _numeric_constants(c, 1 / rho_iv, _TOTAL_Z_ORDER, "rho^-1")
        sigma_p_over_rho = sigma_p_s * rho_inv
        chi_r_s = sigma_p_over_rho - sigma_p_over_rho.shift_hb(-1, origin="chi_r_hb_inverse")
        chi_source = "derived from caller sigma and sigma_prime"
    else:
        chi_s = _sector_row_jet(c, chi, _TOTAL_Z_ORDER, "chi", source_family)
        chi_r_s = _sector_row_jet(c, chi_r, _TOTAL_Z_ORDER, "chi_r", source_family)
        chi_source = "caller supplied pre-factored chi and chi_r"

    # Exact radial identities from collar_fields/first_rows.  All quantities
    # on the right retain their source sectors; no comparison row is returned
    # under an actual field name.
    phi1 = phi0 * chi_s * A_sector
    phi2 = phi0 * (
        chi_r_s * A_sector
        + chi_s * chi_s * A_sector * A_sector
        + chi_s * Ar_sector
    )

    ratio = phi0 / pc[0]
    # The comparison input is physical core U, whose radial derivatives
    # already include epsilon=Lambda^-1. Do not apply it a second time.
    u1 = ratio * chi_s * psi1_sector
    u2 = ratio * (
        chi_r_s * psi1_sector
        + chi_s * (chi_s - 1) * A_sector * psi1_sector
        + chi_s * psi2_sector
    )

    rho_inv = _numeric_constants(c, 1 / rho_iv, _TOTAL_Z_ORDER, "rho^-1")
    mean1 = (u0 - m0) * rho_inv
    mean2 = (u1 - 2 * mean1) * rho_inv

    jets = {
        "Phi": (phi0, phi1.truncate(4), phi2.truncate(3)),
        "U": (u0, u1.truncate(4), u2.truncate(3)),
        "mean": (m0, mean1.truncate(4), mean2.truncate(3)),
    }
    radial_rows = {name: _ordinary_sector_grid(rows) for name, rows in jets.items()}
    chi_grid = _ordinary_sector_grid((chi_s,))
    chi_r_grid = _ordinary_sector_grid((chi_r_s,))

    metadata = {
        "schema": "actual_collar_radial_core_U_v1",
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "sector_schema": {
            "powers": "(hb, F0_anchor, Pstar_squared, Lambda)",
            "coefficient": "directed mpmath interval",
            "origin": "caller/source operation provenance",
        },
        "source_family": source_family,
        "source_atoms": copy.deepcopy(source_atoms),
        "source_atoms_supplied": source_atoms_supplied,
        "source_provenance": copy.deepcopy(source_provenance),
        "rho_interval": rho_iv,
        "Z_interval": z_iv,
        "log_hb_interval_contract": log_hb_iv,
        "actual_inputs": {
            "Phi": "caller-supplied actual Phi zeroth-radial history",
            "U": "caller-supplied actual U zeroth-radial history",
            "mean": "caller-supplied actual own mean zeroth-radial history",
        },
        "comparison_inputs": {
            "Phi": "caller-supplied comparison-core Phi radial grid",
            "U": "caller-supplied comparison-core physical U radial grid",
        },
        "comparison_source_is_not_actual": True,
        "actual_source_identity_proved": False,
        "same_source_contract": "caller-proved; not independently verified",
        "rho_contract": "caller-supplied rho; no radius or source map inferred",
        "geometry_certificate": copy.deepcopy(dict(geometry_certificate)),
        "geometry_certificate_independently_verified": False,
        "hb_formal": True,
        "hb_source_log_used_only_for_contract": True,
        "comparison_U_radial_derivatives_already_include_epsilon": True,
        "comparison_U_source_powers_preserved": True,
        "epsilon_not_applied_twice": True,
        "hb_positive_and_negative_powers_preserved": True,
        "chi_formula": "chi=1+(hb-1)*sigma=(1-sigma)+hb*sigma",
        "chi_r_formula": "chi_r=(hb-1)*sigma_prime/(hb*rho)=sigma_prime/rho-hb^-1*sigma_prime/rho",
        "radial_formulas": {
            "Phi_rho": "Phi*chi*(Phi_core_rho/Phi_core)",
            "Phi_rhorho": "Phi*(chi_r*A+chi^2*A^2+chi*(Phi_core_rhorho/Phi_core-A^2))",
            "U_rho": "(Phi/Phi_core)*chi*U_core_rho",
            "U_rhorho": "(Phi/Phi_core)*(chi_r*U_core_rho+chi*(chi-1)*A*U_core_rho+chi*U_core_rhorho)",
            "mean_rho": "(U-mean)/rho",
            "mean_rhorho": "(U_rho-2*mean_rho)/rho",
        },
        "input_derivatives_are_ordinary": True,
        "internal_taylor_coefficients_divided_by_factorial": True,
        "output_derivatives_are_ordinary": True,
        "output_orders": {str(i): list(range(_TOTAL_Z_ORDER - i + 1)) for i in _RADIAL_ORDERS},
        "required_actual_rows": {"rho0": "Z0..Z5"},
        "required_comparison_rows": {
            "rho0": "Z0..Z5",
            "rho1": "Z0..Z4",
            "rho2": "Z0..Z3",
        },
        "missing_derivatives_rejected": True,
        "unknown_inputs_rejected": True,
        "chi_source": chi_source,
        "actual_point_value_selected": False,
        "comparison_promoted_to_actual": False,
        "old_program_imported": False,
        "saved_packet_read": False,
    }

    return {
        "source_family": source_family,
        "source_atoms": copy.deepcopy(source_atoms),
        "Phi": radial_rows["Phi"],
        "U": radial_rows["U"],
        "mean": radial_rows["mean"],
        "M": radial_rows["mean"],
        "actual_Phi": radial_rows["Phi"],
        "actual_U": radial_rows["U"],
        "actual_mean": radial_rows["mean"],
        "actual": radial_rows,
        "radial_rows": radial_rows,
        "rows": radial_rows,
        "chi": chi_grid,
        "chi_r": chi_r_grid,
        "comparison_roles": metadata["comparison_inputs"],
        "comparison": {
            "role": "comparison-only source used in radial directions",
            "is_actual": False,
            "Phi": "comparison_Phi radial grid supplied by caller",
            "U": "comparison_U radial grid supplied by caller",
        },
        "metadata": metadata,
        "proof": metadata,
    }


def build(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Short public alias for :func:`build_collar_laurent_radial`."""

    return build_collar_laurent_radial(*args, **kwargs)


def collar_laurent_radial(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compatibility alias retaining the descriptive operation name."""

    return build_collar_laurent_radial(*args, **kwargs)


def derive_collar_laurent_radial(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compatibility alias for callers using a derivative-oriented name."""

    return build_collar_laurent_radial(*args, **kwargs)


def build_radial(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compatibility alias for the radial adapter."""

    return build_collar_laurent_radial(*args, **kwargs)
