"""Independent first-order core matrix leaf.

This module is deliberately smaller than the source/matrix leaf used by the
older all-order program.  It computes only the ordinary ``Z0`` and ``Z1``
coefficients of the ten ``C0`` and six ``B1`` raw matrix rows.  The source
grids are read directly, and the formal ``F0`` and ``Lambda`` powers stay on
each returned term; no source or matrix values are selected from an interval.

The local ``Dual`` algebra stores an interval value and its ordinary first
``Z`` derivative.  The ``F0`` connection is represented formally by a value
one and the local logarithmic slope ``ell``.  Thus the derivative of
``Phi/F0`` is ``Phi_Z + ell*Phi`` while ``ell_Z`` is used for its second
ordinary derivative, exactly as required by the raw source formulas.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


__all__ = ["Dual", "build"]


def _is_interval(value: Any) -> bool:
    """Return whether *value* is an interval-like object."""

    return hasattr(value, "_mpi_") or (hasattr(value, "a") and hasattr(value, "b"))


def _bounds(value: Any) -> tuple[Any, Any]:
    """Get interval endpoints without taking a midpoint."""

    if hasattr(value, "a") and hasattr(value, "b"):
        return value.a, value.b
    if hasattr(value, "_mpi_"):
        rows = value._mpi_
        return rows[0], rows[1]
    # A scalar is a degenerate interval for the benefit of a scalar context.
    return value, value


def _interval(context: Any, value: Any, label: str = "interval") -> Any:
    """Coerce a value to the caller's arithmetic context, preserving bounds."""

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
    """Reject a denominator whose directed interval contains zero."""

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
    """Read one supplied ordinary mixed derivative.

    The low matrix uses only the rows needed by its first-order formulas; no
    complete-grid validation or zero padding is performed here.
    """

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


class Dual:
    """A first-order directed-interval dual ``(value, d/dZ value)``."""

    __slots__ = ("context", "value", "derivative")

    def __init__(self, context: Any, value: Any, derivative: Any = 0) -> None:
        if not hasattr(context, "mpf"):
            raise TypeError("c must provide directed interval mpf arithmetic")
        self.context = context
        self.value = _interval(context, value, "dual value")
        self.derivative = _interval(context, derivative, "dual derivative")

    @classmethod
    def constant(cls, context: Any, value: Any) -> "Dual":
        return cls(context, value, context.mpf(0))

    @classmethod
    def variable(cls, context: Any, value: Any) -> "Dual":
        return cls(context, value, context.mpf(1))

    def _coerce(self, other: Any) -> "Dual":
        if isinstance(other, Dual):
            if other.context is not self.context:
                raise TypeError("dual values must use one interval context")
            return other
        return Dual.constant(self.context, other)

    def __add__(self, other: Any) -> "Dual":
        other = self._coerce(other)
        return Dual(self.context, self.value + other.value, self.derivative + other.derivative)

    __radd__ = __add__

    def __neg__(self) -> "Dual":
        return Dual(self.context, -self.value, -self.derivative)

    def __sub__(self, other: Any) -> "Dual":
        return self + (-self._coerce(other))

    def __rsub__(self, other: Any) -> "Dual":
        return self._coerce(other) - self

    def __mul__(self, other: Any) -> "Dual":
        other = self._coerce(other)
        return Dual(
            self.context,
            self.value * other.value,
            self.derivative * other.value + self.value * other.derivative,
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "Dual":
        _require_nonzero(self.value, "dual denominator")
        inverse = self.context.mpf(1) / self.value
        return Dual(
            self.context,
            inverse,
            -self.derivative / (self.value ** 2),
        )

    def __truediv__(self, other: Any) -> "Dual":
        return self * self._coerce(other).reciprocal()

    def __rtruediv__(self, other: Any) -> "Dual":
        return self._coerce(other) * self.reciprocal()

    def __pow__(self, power: int) -> "Dual":
        if not isinstance(power, int) or power < 0:
            raise ValueError("only nonnegative integer dual powers are supported")
        result = Dual.constant(self.context, 1)
        for _ in range(power):
            result = result * self
        return result


def _square(value: Dual) -> Dual:
    """Square a dual while retaining a tight interval base square."""

    return Dual(
        value.context,
        value.value ** 2,
        2 * value.value * value.derivative,
    )


def _term(f0_power: int, lambda_power: int, value: Dual, order: int, label: str) -> dict[str, Any]:
    if f0_power not in (0, 1) or lambda_power > 0:
        raise ValueError(f"invalid formal powers for {label}")
    if order == 0:
        coefficient = value.value
    elif order == 1:
        coefficient = value.derivative
    else:
        raise ValueError("the low matrix exports only ordinary Z0 and Z1")
    return {
        "F0_power": f0_power,
        "Lambda_power": lambda_power,
        "coefficient": coefficient,
    }


def _read(context: Any, grid: Any, radial_order: int, axial_order: int, label: str) -> Any:
    return _interval(
        context,
        _source_value(grid, radial_order, axial_order, label),
        f"{label} rho^{radial_order} Z^{axial_order}",
    )


def build(
    c: Any,
    rho: Any,
    z: Any,
    delta: Any,
    Lambda: Any,
    j: Any,
    Phi_grid: Any,
    U_grid: Any,
    M_grid: Any,
) -> dict[str, Any]:
    """Build the ordinary ``Z0``/``Z1`` raw ``C0`` and ``B1`` matrix rows.

    ``Phi_grid``, ``U_grid`` and ``M_grid`` provide ordinary mixed source
    derivatives.  At minimum this leaf reads ``rho0_Z0..2`` and
    ``rho1_Z0..1`` for ``Phi``/``U``, and ``rho0_Z0..2`` for ``M``.  Keys of
    the form ``rho{i}_Z{k}`` are accepted directly, along with the equivalent
    tuple, nested-row, and legacy spelling forms handled locally here.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "sqrt"):
        raise TypeError("c must provide directed interval mpf arithmetic and sqrt")

    rho_value = _interval(c, rho, "rho")
    z_value = _interval(c, z, "z")
    delta_value = _interval(c, delta, "delta")
    lambda_value = _interval(c, Lambda, "Lambda")
    j_value = _interval(c, j, "j")

    # Only these source rows are needed for the ten/six first-order matrix
    # formulas.  In particular, no radial-second rows and no Z3+ rows are
    # touched by this leaf.
    phi0_0 = _read(c, Phi_grid, 0, 0, "Phi_grid")
    phi0_1 = _read(c, Phi_grid, 0, 1, "Phi_grid")
    phi0_2 = _read(c, Phi_grid, 0, 2, "Phi_grid")
    phi1_0 = _read(c, Phi_grid, 1, 0, "Phi_grid")
    phi1_1 = _read(c, Phi_grid, 1, 1, "Phi_grid")

    u0_0 = _read(c, U_grid, 0, 0, "U_grid")
    u0_1 = _read(c, U_grid, 0, 1, "U_grid")
    u0_2 = _read(c, U_grid, 0, 2, "U_grid")
    u1_0 = _read(c, U_grid, 1, 0, "U_grid")
    u1_1 = _read(c, U_grid, 1, 1, "U_grid")

    m0_0 = _read(c, M_grid, 0, 0, "M_grid")
    m0_1 = _read(c, M_grid, 0, 1, "M_grid")
    m0_2 = _read(c, M_grid, 0, 2, "M_grid")

    one = Dual.constant(c, 1)
    rho_dual = Dual.constant(c, rho_value)
    z_dual = Dual.variable(c, z_value)
    delta_dual = Dual.constant(c, delta_value)
    lambda_dual = Dual.constant(c, lambda_value)
    j_dual = Dual.constant(c, j_value)
    d = one - _square(z_dual)
    L = one - delta_dual * _square(z_dual)
    _require_nonzero(L.value, "L=1-delta*z^2")

    # The F0 value is formal one.  ell is itself a dual, so its derivative is
    # the required ell_Z; the relative F0 factor has pair [1, ell].
    H0 = z_dual * (one - delta_dual) / 2 + d * (4 * z_dual + j_dual)
    sigma = j_value / c.mpf(500)
    denominator = _square(H0) + _square(Dual.constant(c, sigma))
    _require_nonzero(denominator.value, "H^2+(j/500)^2")
    ell = -(lambda_dual * L * H0) / denominator
    relative_f0 = Dual(c, c.mpf(1), ell.value)

    # f and fr are Phi/F0 duals.  f_z carries the second ordinary derivative
    # explicitly, including ell_Z and ell^2, rather than asking for an extra
    # source row or a second-order output.
    phi0 = Dual(c, phi0_0, phi0_1)
    phi1 = Dual(c, phi1_0, phi1_1)
    f = phi0 * relative_f0
    fr = phi1 * relative_f0
    fz = Dual(
        c,
        f.derivative,
        phi0_2
        + 2 * ell.value * phi0_1
        + (ell.derivative + ell.value ** 2) * phi0_0,
    )

    # U/W, radial U/W, and M/Mz are all first-order duals.  The explicit
    # second axial rows only supply the derivative slots of wz and Mz.
    w = Dual(c, u0_0, u0_1)
    wr = Dual(c, u1_0, u1_1)
    wz = Dual(c, u0_1, u0_2)
    m = Dual(c, m0_0, m0_1)
    mz = Dual(c, m0_1, m0_2)

    Q = (2 * z_dual * w - (one - delta_dual) * z_dual * m - d * mz) / L
    H = z_dual * (one - delta_dual) / 2 + d * w
    B = Q + (one - 2 * z_dual * w) / L
    E = z_dual * w - c.mpf("0.5")
    Cf = f + rho_dual * fr
    Zf = (
        (-2 - delta_dual) * z_dual * f
        + d * fz
        - 2 * z_dual * rho_dual * fr
    ) / L
    Zw = (
        (-1 - delta_dual) * z_dual * w
        + d * wz
        - 2 * z_dual * rho_dual * wr
    ) / L

    x = c.sqrt(rho_value)
    c0_jets = {
        "4,1": (1, -1, f * (4 * x)),
        "5,1": (0, -1, 2 * (Q + E * (-2 + delta_dual) / L)),
        "5,2": (1, -1, 2 * Zf - 2 * (-1 + delta_dual) * z_dual * Cf / L),
        "5,3": (1, -1, -2 * (1 + delta_dual) * z_dual * Cf / L),
        "5,5": (0, -1, B * x),
        "6,1": (1, -2, -8 * z_dual * rho_dual * f / L),
        "6,2": (
            0,
            -1,
            2
            * (
                E * (-1 + delta_dual) / L
                + Zw
                - wr * z_dual * (-1 + delta_dual) * rho_dual / L
            ),
        ),
        "6,3": (0, -1, -2 * (1 + delta_dual) * z_dual * rho_dual * wr / L),
        "6,4": (0, -1, -4 * z_dual / L),
        "6,6": (0, -1, B * x),
    }
    b1_jets = {
        "5,1": (0, -1, 2 * H / L),
        "5,2": (1, -1, -2 * d * Cf / L),
        "5,3": (1, -1, -2 * d * Cf / L),
        "6,2": (0, -1, 2 * (H - wr * d * rho_dual) / L),
        "6,3": (0, -1, -2 * d * rho_dual * wr / L),
        "6,4": (0, -1, 2 * d / L),
    }

    return {
        "C0": {
            str(order): {
                key: [_term(f0_power, lambda_power, value, order, f"C0 {key}")]
                for key, (f0_power, lambda_power, value) in c0_jets.items()
            }
            for order in (0, 1)
        },
        "B1": {
            str(order): {
                key: [_term(f0_power, lambda_power, value, order, f"B1 {key}")]
                for key, (f0_power, lambda_power, value) in b1_jets.items()
            }
            for order in (0, 1)
        },
    }
