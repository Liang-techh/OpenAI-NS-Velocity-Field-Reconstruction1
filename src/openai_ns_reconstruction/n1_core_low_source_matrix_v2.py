"""Independent numerical low-order core source and raw matrix leaf.

The caller owns the source provenance.  ``build`` only consumes directed
intervals for the ordinary mixed derivatives of ``Phi``, ``U`` and ``M`` and
returns formal ``F0``/``Lambda`` sectors.  Taylor coefficients are divided by
factorials internally; every returned coefficient is converted back to an
ordinary ``Z`` derivative.

No producer module, persisted artifact, or floating-point representative is
used here.  The only required context operation is ``c.mpf`` (plus the usual
directed interval arithmetic and ``sqrt``).
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


__all__ = [
    "IntervalTaylorJet",
    "TaylorJet",
    "build",
    "build_source",
    "build_matrix",
]


_TOTAL_ORDER = 5
_OUTPUT_ORDER = 2
_RADIAL_ROWS = ((0, 5), (1, 4), (2, 3))


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") or (hasattr(value, "a") and hasattr(value, "b"))


def _bounds(value: Any) -> tuple[Any, Any]:
    if not _is_interval(value):
        raise TypeError("a directed interval is required")
    try:
        return value.a, value.b
    except AttributeError:
        rows = value._mpi_
        return rows[0], rows[1]


def _interval(context: Any, value: Any, label: str = "interval") -> Any:
    """Coerce without selecting a point from an interval."""

    if _is_interval(value):
        lo, hi = _bounds(value)
        return context.mpf([lo, hi])
    if isinstance(value, Mapping):
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


def _mapping_get(mapping: Mapping[Any, Any], names: Sequence[Any]) -> Any | None:
    for name in names:
        try:
            if name in mapping:
                return mapping[name]
        except (TypeError, ValueError):
            continue
    return None


def _source_value(grid: Any, radial_order: int, axial_order: int, label: str) -> Any:
    """Read one ordinary derivative; missing entries are never padded."""

    if grid is None:
        raise TypeError(f"{label} source grid is required")
    if isinstance(grid, Mapping):
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
        row = _mapping_get(
            grid,
            (radial_order, str(radial_order), f"rho{radial_order}", f"rho_{radial_order}"),
        )
        if isinstance(row, Mapping):
            nested = _mapping_get(
                row,
                (axial_order, str(axial_order), f"Z{axial_order}", f"z{axial_order}"),
            )
            if nested is not None:
                return nested
    if isinstance(grid, Sequence) and not isinstance(grid, (str, bytes, bytearray)):
        if radial_order < len(grid):
            row = grid[radial_order]
            if isinstance(row, Sequence) and not isinstance(row, (str, bytes, bytearray)):
                if axial_order < len(row):
                    return row[axial_order]
    raise KeyError(f"missing {label} ordinary derivative rho^{radial_order} Z^{axial_order}")


def _validate_grid(context: Any, grid: Any, label: str) -> None:
    for radial_order, axial_max in _RADIAL_ROWS:
        for axial_order in range(axial_max + 1):
            _interval(
                context,
                _source_value(grid, radial_order, axial_order, label),
                f"{label} rho^{radial_order} Z^{axial_order}",
            )


def _source_jet(
    context: Any,
    grid: Any,
    radial_order: int,
    axial_order: int,
    label: str,
) -> "IntervalTaylorJet":
    coefficients = []
    for order in range(axial_order + 1):
        value = _interval(
            context,
            _source_value(grid, radial_order, order, label),
            f"{label} rho^{radial_order} Z^{order}",
        )
        coefficients.append(value / context.mpf(math.factorial(order)))
    return IntervalTaylorJet(context, coefficients)


class IntervalTaylorJet:
    """Truncated directed-interval Taylor algebra.

    ``coefficients[k]`` encloses ``f^(k)(Z)/k!``.  Operations truncate to the
    smallest known order.  Only constants are extended with known zero
    derivatives; source rows are read explicitly and never zero-filled.
    """

    __slots__ = ("context", "coefficients")

    def __init__(self, context: Any, coefficients: Sequence[Any]):
        values = tuple(coefficients)
        if not values:
            raise ValueError("a Taylor jet needs a constant coefficient")
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
            raise ValueError("requested Taylor order is not available")
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
                raise TypeError("Taylor jets must use one interval context")
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
        coefficients = [self.context.mpf(1) / self[0]]
        for index in range(1, self.order + 1):
            value = self.context.mpf(0)
            for left in range(1, index + 1):
                value += self[left] * coefficients[index - left]
            coefficients.append(-value / self[0])
        return IntervalTaylorJet(self.context, coefficients)

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


def _square(jet: IntervalTaylorJet) -> IntervalTaylorJet:
    """Square a jet, using interval ``** 2`` for its base coefficient."""

    product = jet * jet
    coefficients = list(product.coefficients)
    coefficients[0] = jet[0] ** 2
    return IntervalTaylorJet(jet.context, coefficients)


def _axial_second(
    field: IntervalTaylorJet,
    field_r: IntervalTaylorJet,
    field_rr: IntervalTaylorJet,
    z: IntervalTaylorJet,
    rho: IntervalTaylorJet,
    delta: IntervalTaylorJet,
    weight: IntervalTaylorJet,
) -> IntervalTaylorJet:
    """Raw ``Z_(a-1+delta) Z_a`` operator from the original source."""

    one = IntervalTaylorJet.constant(z.context, 1, z.order)
    d = one - z * z
    L = one - delta * z * z
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
    return (
        (weight - 1 + delta) * z * A / L
        + d * (A_Z / L + 2 * delta * z * A / (L * L))
        - 2 * z * rho * A_R / L
    ) / L


def _relative_f0_jet(
    context: Any,
    z: IntervalTaylorJet,
    delta: IntervalTaylorJet,
    Lambda: IntervalTaylorJet,
    required_j: Any,
    order: int,
) -> IntervalTaylorJet:
    """Build the current-point jet ``F0(Z+h)/F0(Z)``."""

    one = IntervalTaylorJet.constant(context, 1, order)
    d = one - _square(z)
    L = one - delta * _square(z)
    j = IntervalTaylorJet.constant(context, required_j, order)
    H = z * (one - delta) / 2 + d * (4 * z + j)
    sigma = _interval(context, required_j, "j") / context.mpf(500)
    denominator = _square(H) + sigma ** 2
    _require_nonzero(denominator[0], "H^2+(j/500)^2")
    ell = -(Lambda * L * H) / denominator
    coefficients = [context.mpf(1)]
    for n in range(1, order + 1):
        value = context.mpf(0)
        for k in range(n):
            value += ell[k] * coefficients[n - 1 - k]
        coefficients.append(value / context.mpf(n))
    return IntervalTaylorJet(context, coefficients)


def _ordinary(jet: IntervalTaylorJet, order: int, label: str) -> Any:
    if jet.order < order:
        raise ValueError(f"{label} has no ordinary Z^{order} derivative")
    return jet[order] * math.factorial(order)


def _term(f0_power: int, Lambda_power: int, jet: IntervalTaylorJet, order: int, label: str) -> dict[str, Any]:
    if f0_power not in (0, 1) or Lambda_power > 0:
        raise ValueError("invalid formal source exponents")
    return {
        "F0_power": f0_power,
        "Lambda_power": Lambda_power,
        "coefficient": _ordinary(jet, order, label),
    }


def _assemble(
    context: Any,
    rho: Any,
    Z: Any,
    delta: Any,
    Lambda: Any,
    required_j: Any,
    Phi_grid: Any,
    U_grid: Any,
    M_grid: Any,
    *,
    want_g: bool,
    want_matrix: bool,
) -> dict[str, Any]:
    if not hasattr(context, "mpf"):
        raise TypeError("c must provide directed interval mpf arithmetic")
    if want_g is False and want_matrix is False:
        raise ValueError("at least one of source_only or matrix_only must be false")

    # Validate the complete supplied total-order-five grids before any algebra.
    # This makes a missing high derivative an explicit input error rather than
    # an implicit zero.
    _validate_grid(context, Phi_grid, "Phi_grid")
    _validate_grid(context, U_grid, "U_grid")
    _validate_grid(context, M_grid, "M_grid")

    rho_iv = _interval(context, rho, "rho")
    z_iv = _interval(context, Z, "Z")
    delta_iv = _interval(context, delta, "delta")
    lambda_iv = _interval(context, Lambda, "Lambda")
    j_iv = _interval(context, required_j, "j")
    degree = _TOTAL_ORDER
    rho_j = IntervalTaylorJet.constant(context, rho_iv, degree)
    z_j = IntervalTaylorJet.variable(context, z_iv, degree)
    delta_j = IntervalTaylorJet.constant(context, delta_iv, degree)
    lambda_j = IntervalTaylorJet.constant(context, lambda_iv, degree)
    one = IntervalTaylorJet.constant(context, 1, degree)
    d_j = one - _square(z_j)
    L_j = one - delta_j * _square(z_j)
    _require_nonzero(L_j[0], "L=1-delta*Z^2")

    # Whole-core source rows use ordinary derivatives and are converted once.
    U0 = _source_jet(context, U_grid, 0, 5, "U_grid")
    U1 = _source_jet(context, U_grid, 1, 4, "U_grid")
    U2 = _source_jet(context, U_grid, 2, 3, "U_grid")
    M0 = _source_jet(context, M_grid, 0, 5, "M_grid")
    M1 = _source_jet(context, M_grid, 1, 4, "M_grid")
    Phi0 = _source_jet(context, Phi_grid, 0, 5, "Phi_grid")
    Phi1 = _source_jet(context, Phi_grid, 1, 4, "Phi_grid")
    Phi2 = _source_jet(context, Phi_grid, 2, 3, "Phi_grid")

    M0_Z = M0.differentiate()
    M1_Z = M1.differentiate()
    U1_Z = U1.differentiate()
    Q = (2 * z_j * U0 - (1 - delta_j) * z_j * M0 - d_j * M0_Z) / L_j
    Qr = (2 * z_j * U1 - (1 - delta_j) * z_j * M1 - d_j * M1_Z) / L_j
    rhoQrr = (
        2 * z_j * rho_j * U2
        - (1 - delta_j) * z_j * (U1 - 2 * M1)
        - d_j * (U1_Z - 2 * M1_Z)
    ) / L_j
    Q_Z = Q.differentiate()
    S = (
        (Q + (1 - delta_j) * z_j * Q_Z / 2 + rho_j * Qr) / L_j
        + Q * (Q / 2 + rho_j * Qr)
        + U0 * (-2 * z_j * Q + d_j * Q_Z - 2 * z_j * rho_j * Qr) / L_j
    )
    P = 2 * Qr + rhoQrr

    # Compute the current-point F0-normalized radial jets once.  Both source
    # and matrix sectors consume this same product; rebuilding it in the two
    # branches would repeat the F0 connection algebra.
    relative_f0 = _relative_f0_jet(
        context, z_j, delta_j, lambda_j, j_iv, degree
    )
    F0_ratio = (Phi0 * relative_f0).truncate(5)
    F1_ratio = (Phi1 * relative_f0).truncate(4)
    F2_ratio = (Phi2 * relative_f0).truncate(3)

    result: dict[str, Any] = {}
    if want_g:
        # The relative jet already includes the F0 connection exactly once.
        weight_F = IntervalTaylorJet.constant(context, -2, degree) - delta_j
        weight_U = IntervalTaylorJet.constant(context, -1, degree) - delta_j
        Nt_over_F0 = -_axial_second(
            F0_ratio, F1_ratio, F2_ratio, z_j, rho_j, delta_j, weight_F
        )
        Nz = -_axial_second(U0, U1, U2, z_j, rho_j, delta_j, weight_U)
        x_j = IntervalTaylorJet.constant(context, context.sqrt(rho_iv), degree)
        g_rows: dict[str, dict[str, list[dict[str, Any]]]] = {}
        for order in range(_OUTPUT_ORDER + 1):
            g_rows[str(order)] = {
                "4": [
                    _term(0, 0, 2 * x_j * P, order, "g4 Lambda^0"),
                    _term(0, -1, -x_j * S, order, "g4 Lambda^-1"),
                ],
                "5": [_term(1, -1, 2 * Nt_over_F0, order, "g5 F0 coefficient")],
                "6": [
                    _term(0, -1, 2 * Nz - 4 * rho_j * z_j * P / L_j, order, "g6 Lambda^-1"),
                    _term(0, -2, 2 * rho_j * z_j * S / L_j, order, "g6 Lambda^-2"),
                ],
            }
        result["g"] = g_rows

    if want_matrix:
        # F/F0, Fr/F0 and their Z derivatives are used directly.  No extra
        # F0 connection is applied after this product.
        f = (Phi0 * relative_f0).truncate(5)
        fr = (Phi1 * relative_f0).truncate(4)
        frr = (Phi2 * relative_f0).truncate(3)
        fz = f.differentiate()
        w = U0
        wr = U1
        wz = w.differentiate()
        H = z_j * (1 - delta_j) / 2 + d_j * w
        B = Q + (1 - 2 * z_j * w) / L_j
        E = z_j * w - context.mpf("0.5")
        Cf = f + rho_j * fr
        Zf = ((-2 - delta_j) * z_j * f + d_j * fz - 2 * z_j * rho_j * fr) / L_j
        Zw = ((-1 - delta_j) * z_j * w + d_j * wz - 2 * z_j * rho_j * wr) / L_j

        c0_jets = {
            "4,1": (1, -1, f * (4 * context.sqrt(rho_iv))),
            "5,1": (0, -1, 2 * (Q + E * (-2 + delta_j) / L_j)),
            "5,2": (1, -1, 2 * Zf - 2 * (-1 + delta_j) * z_j * Cf / L_j),
            "5,3": (1, -1, -2 * (1 + delta_j) * z_j * Cf / L_j),
            "5,5": (0, -1, B * context.sqrt(rho_iv)),
            "6,1": (1, -2, -8 * z_j * rho_j * f / L_j),
            "6,2": (0, -1, 2 * (E * (-1 + delta_j) / L_j + Zw - wr * z_j * (-1 + delta_j) * rho_j / L_j)),
            "6,3": (0, -1, -2 * (1 + delta_j) * z_j * rho_j * wr / L_j),
            "6,4": (0, -1, -4 * z_j / L_j),
            "6,6": (0, -1, B * context.sqrt(rho_iv)),
        }
        b1_jets = {
            "5,1": (0, -1, 2 * H / L_j),
            "5,2": (1, -1, -2 * d_j * Cf / L_j),
            "5,3": (1, -1, -2 * d_j * Cf / L_j),
            "6,2": (0, -1, 2 * (H - wr * d_j * rho_j) / L_j),
            "6,3": (0, -1, -2 * d_j * rho_j * wr / L_j),
            "6,4": (0, -1, 2 * d_j / L_j),
        }
        c0_rows: dict[str, dict[str, list[dict[str, Any]]]] = {}
        b1_rows: dict[str, dict[str, list[dict[str, Any]]]] = {}
        for order in range(_OUTPUT_ORDER + 1):
            c0_rows[str(order)] = {
                key: [_term(f0_power, lambda_power, jet, order, f"C0 {key}")]
                for key, (f0_power, lambda_power, jet) in c0_jets.items()
            }
            b1_rows[str(order)] = {
                key: [_term(f0_power, lambda_power, jet, order, f"B1 {key}")]
                for key, (f0_power, lambda_power, jet) in b1_jets.items()
            }
        result["C0"] = c0_rows
        result["B1"] = b1_rows

    return result


def build(
    c: Any,
    rho: Any,
    Z: Any,
    delta: Any,
    Lambda: Any,
    j: Any,
    Phi_grid: Any,
    U_grid: Any,
    M_grid: Any,
    *,
    source_only: bool = False,
    matrix_only: bool = False,
) -> dict[str, Any]:
    """Return ordinary ``Z0/Z1/Z2`` source and raw dynamic matrix rows.

    The grids contain ordinary derivatives keyed by ``rho{i}_Z{k}``, by
    ``(i,k)``, or as nested row sequences.  Each requested field is required
    through total order five: ``i=0..2`` with ``k<=5-i``.  ``source_only``
    returns only ``g``; ``matrix_only`` returns only ``C0`` and ``B1``.
    """

    if source_only and matrix_only:
        raise ValueError("source_only and matrix_only are mutually exclusive")
    return _assemble(
        c,
        rho,
        Z,
        delta,
        Lambda,
        j,
        Phi_grid,
        U_grid,
        M_grid,
        want_g=not matrix_only,
        want_matrix=not source_only,
    )


def build_source(*args: Any, **kwargs: Any) -> dict[str, Any]:
    kwargs["source_only"] = True
    return build(*args, **kwargs)


def build_matrix(*args: Any, **kwargs: Any) -> dict[str, Any]:
    kwargs["matrix_only"] = True
    return build(*args, **kwargs)
