"""Source-factored first-collar forcing ordinary ``Z`` jets.

This is a dependency-free source leaf.  The caller supplies the actual,
same-source ``Phi``, ``U`` and radial-mean ordinary mixed derivative grids;
there is no import of the core producer and no reconstruction from a stored
bound or result packet.  The algebra is the static ``n1_core_forcing_z3``
algebra extended to the first collar ``4 <= rho <= 4.1`` and to ordinary
``Z`` orders zero through three.

Source rows are ordinary derivatives.  Internally they become factorial-
normalised Taylor coefficients, and are converted back to ordinary
derivatives only when a sector row is emitted.  Every emitted row retains its
formal ``F0_power`` and ``Lambda_power`` and its signed directed coefficient;
sectors with equal powers are deliberately not combined. Input derivative
grids are scalar intervals: this leaf does not preserve formal hb powers
within them, and it does not certify caller source identity.

``Lambda`` is enclosed numerically only while differentiating the relative
positive ``F0`` jet.  The outer Lambda powers in the returned sector rows are
formal metadata.  The supplied log ``F0`` interval is retained as metadata;
no positive ``F0`` value, error radius, midpoint, or upper endpoint is used
as a field value.  ``x`` is required from the caller and is never inferred
from ``rho``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


__all__ = [
    "IntervalTaylorJet",
    "TaylorJet",
    "build_first_collar_forcing_jets",
    "first_collar_forcing_jets",
    "derive_first_collar_forcing_jets",
    "build_first_collar_forcing_z3",
    "first_collar_forcing_z3",
]


_MAX_SOURCE_ORDER = 5
_OUTPUT_ORDERS = (0, 1, 2, 3)
_RHO_MIN = "4"
_RHO_MAX = "4.1"
_Z_MIN = "-1"
_Z_MAX = "1"


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") or (hasattr(value, "a") and hasattr(value, "b"))


def _bounds(value: Any) -> tuple[Any, Any]:
    """Return directed endpoints without choosing an interval representative."""

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


class IntervalTaylorJet:
    """Truncated interval Taylor algebra with factorial-normalised terms."""

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
    """Read one required ordinary derivative; never synthesize a missing row."""

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


def _source_jet(
    context: Any,
    grid: Any,
    radial_order: int,
    axial_order: int,
    label: str,
) -> IntervalTaylorJet:
    coefficients = []
    for order in range(axial_order + 1):
        value = _interval(
            context,
            _source_value(grid, radial_order, order, label),
            f"{label} rho^{radial_order} Z^{order}",
        )
        coefficients.append(_factorial_div(value, order, context))
    return IntervalTaylorJet(context, coefficients)


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
    """Accept both keyword-oriented and natural positional source calls."""

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

    # New positional form is (x, Phi, Uz, mean, logF0), while the convenient
    # source-first form is (Phi, Uz, mean, x, logF0).  A lone value after
    # three grids is interpreted as logF0 when x was already named.
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


def _axial_second(
    field: IntervalTaylorJet,
    field_r: IntervalTaylorJet,
    field_rr: IntervalTaylorJet,
    z: IntervalTaylorJet,
    rho: IntervalTaylorJet,
    delta: IntervalTaylorJet,
    weight: IntervalTaylorJet,
) -> IntervalTaylorJet:
    """The inner-Volterra axial second operator, in Taylor form."""

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
    """Build ``F0(Z+h)/F0(Z)`` from its supplied log-derivative formula."""

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


def _ordinary_derivative(jet: IntervalTaylorJet, order: int, label: str) -> Any:
    if jet.order < order:
        raise ValueError(f"{label} only has Taylor order {jet.order}; Z^{order} is unavailable")
    return jet[order] * math.factorial(order)


def _row(F0_power: int, Lambda_power: int, coefficient: Any) -> dict[str, Any]:
    if F0_power not in (0, 1):
        raise ValueError("only F0 powers zero and one are admitted")
    if Lambda_power > 0:
        raise ValueError("the external Lambda exponent may not be positive")
    return {
        "F0_power": F0_power,
        "Lambda_power": Lambda_power,
        "coefficient": coefficient,
    }


def build_first_collar_forcing_jets(
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
    **unknown: Any,
) -> dict[str, Any]:
    """Return unmerged signed ``g4/g5/g6`` source sectors for ``Z0..Z3``.

    The required source grids contain ordinary mixed derivatives indexed by
    ``(rho_order, Z_order)`` or equivalent nested/row forms.  The required
    rows are ``i=0,1,2`` through ``k <= 5-i``.  Positional source arguments
    are accepted in either ``(x, Phi, Uz, mean, logF0)`` or
    ``(Phi, Uz, mean, x, logF0)`` form; keyword arguments are preferred.
    """

    if unknown:
        names = ", ".join(sorted(unknown))
        raise TypeError(f"unknown first-collar forcing input(s): {names}")
    if not hasattr(c, "mpf"):
        raise TypeError("c must be the caller's MPIntervalContext")

    x, Phi_derivatives, Uz_derivatives, radial_average_Uz_derivatives, log_F0 = (
        _parse_positional_inputs(
            positional,
            x=x,
            Phi_derivatives=Phi_derivatives,
            Uz_derivatives=Uz_derivatives,
            radial_average_Uz_derivatives=radial_average_Uz_derivatives,
            log_F0=log_F0,
        )
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
    rho_j = IntervalTaylorJet.constant(c, rho_iv, degree)
    z_j = IntervalTaylorJet.variable(c, z_iv, degree)
    delta_j = IntervalTaylorJet.constant(c, delta_iv, degree)
    lambda_j = IntervalTaylorJet.constant(c, lambda_iv, degree)
    one = IntervalTaylorJet.constant(c, 1, degree)
    d_j = one - z_j * z_j
    L_j = one - delta_j * z_j * z_j
    _require_nonzero(L_j[0], "L=1-delta*Z^2")

    # These exact source orders are the minimum needed by an ordinary Z^3
    # output after the axial second operator and Q/M_Z differentiation.
    U0 = _source_jet(c, uz_grid, 0, 5, "Uz_derivatives")
    U1 = _source_jet(c, uz_grid, 1, 4, "Uz_derivatives")
    U2 = _source_jet(c, uz_grid, 2, 3, "Uz_derivatives")
    M0 = _source_jet(c, mean_grid, 0, 5, "radial_average_Uz_derivatives")
    M1 = _source_jet(c, mean_grid, 1, 4, "radial_average_Uz_derivatives")
    Phi0 = _source_jet(c, phi_grid, 0, 5, "Phi_derivatives")
    Phi1 = _source_jet(c, phi_grid, 1, 4, "Phi_derivatives")
    Phi2 = _source_jet(c, phi_grid, 2, 3, "Phi_derivatives")

    # rho*m_rr = U_r - 2*m_r avoids division by rho, including at the axis.
    M0_Z = M0.differentiate()
    M1_Z = M1.differentiate()
    U1_Z = U1.differentiate()
    Q0 = (2 * z_j * U0 - (1 - delta_j) * z_j * M0 - d_j * M0_Z) / L_j
    Qr = (2 * z_j * U1 - (1 - delta_j) * z_j * M1 - d_j * M1_Z) / L_j
    rhoQrr = (
        2 * z_j * rho_j * U2
        - (1 - delta_j) * z_j * (U1 - 2 * M1)
        - d_j * (U1_Z - 2 * M1_Z)
    ) / L_j
    Q0_Z = Q0.differentiate()
    S = (
        (Q0 + (1 - delta_j) * z_j * Q0_Z / 2 + rho_j * Qr) / L_j
        + Q0 * (Q0 / 2 + rho_j * Qr)
        + U0 * (-2 * z_j * Q0 + d_j * Q0_Z - 2 * z_j * rho_j * Qr) / L_j
    )
    P = 2 * Qr + rhoQrr

    # x is deliberately caller-bound; no sqrt(rho) fallback is permitted.
    x_j = IntervalTaylorJet.constant(c, x_iv, degree)
    g4_power0 = 2 * x_j * P
    g4_power1 = -x_j * S

    relative_f0 = _relative_f0_jet(c, z_j, delta_j, lambda_j, j_iv, degree)
    F0_ratio = (Phi0 * relative_f0).truncate(5)
    F1_ratio = (Phi1 * relative_f0).truncate(4)
    F2_ratio = (Phi2 * relative_f0).truncate(3)
    weight_F = IntervalTaylorJet.constant(c, -2, degree) - delta_j
    weight_U = IntervalTaylorJet.constant(c, -1, degree) - delta_j
    Nt_over_F0_anchor = -_axial_second(
        F0_ratio, F1_ratio, F2_ratio, z_j, rho_j, delta_j, weight_F
    )
    Nz = -_axial_second(U0, U1, U2, z_j, rho_j, delta_j, weight_U)

    powers = {
        "g4": ((0, 0, g4_power0), (0, -1, g4_power1)),
        "g5": ((1, -1, 2 * Nt_over_F0_anchor),),
        "g6": (
            (0, -1, 2 * Nz - 4 * rho_j * z_j * P / L_j),
            (0, -2, 2 * rho_j * z_j * S / L_j),
        ),
    }
    rows: dict[str, list[dict[str, Any]]] = {}
    for family, terms in powers.items():
        for order in _OUTPUT_ORDERS:
            rows[f"{family}_Z{order}"] = [
                _row(f0_power, lambda_power, _ordinary_derivative(jet, order, family))
                for f0_power, lambda_power, jet in terms
            ]

    proof = {
        "kind": "first_collar_forcing_jets_z0_z3_v1",
        "same_source_identity_proved": False,
        "hb_factors_inside_input_grids_preserved": False,
        "input_grids_certified_by_caller": True,
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
        "denominator_checks": [
            "L=1-delta*Z^2 excludes zero",
            "H^2+sigma^2 excludes zero",
        ],
        "epsilon_source_factor": "epsilon=Lambda^-1 retained only as Lambda_power",
        "F0_power_is_formal": True,
        "F0_current_Z_factor_formal": True,
        "F0_log_interval_supplied": True,
        "F0_log_interval": log_f0_iv,
        "F0_upper_not_used_as_value": True,
        "F0_absolute_exp_materialized": False,
        "no_input_midpoint_or_error_bound_chosen_by_leaf": True,
        "no_sector_combination": True,
        "x_caller_supplied": True,
        "x_inferred_from_rho": False,
        "positive_x_squared_encloses_rho": True,
        "g5_anchor_rule": "F(Z+h)/F0(Z) jet at each current Z; no extra D_F",
        "Lambda_inside_relative_F0_derivative_coefficients_is_enclosed_numerically": True,
        "outer_Lambda_power_remains_formal": True,
        "lambda_role_source": (
            "caller Lambda interval is used in ell for the relative F0 jet; "
            "row Lambda_power is not numerically materialized"
        ),
        "formula_source": "static n1_core_forcing_z3 algebra, independently transcribed",
        "g4_formula": "2*epsilon*x*(-S/2+Lambda*P)",
        "g5_formula": "2*epsilon*(-axial_second(F))",
        "g6_formula": "2*epsilon*Nz-4*z*rho*epsilon^2*Np/L",
        "Q_formula": "Q0=(2*z*U0-(1-delta)*z*M0-d*M0_Z)/L",
        "Qr_formula": "Qr=(2*z*U1-(1-delta)*z*M1-d*M1_Z)/L",
        "rhoQrr_formula": (
            "rho*Qrr=(2*z*rho*U2-(1-delta)*z*(U1-2*M1)-d*(U1_Z-2*M1_Z))/L"
        ),
        "P_formula": "P=2*Qr+rho*Qrr",
        "S_formula": (
            "(Q0+(1-delta)*z*Q0_Z/2+rho*Qr)/L+Q0*(Q0/2+rho*Qr)+"
            "U0*(-2*z*Q0+d*Q0_Z-2*z*rho*Qr)/L"
        ),
        "axial_second_formula": "[(a-1+delta)*z*A/L+d*(AZ/L+2*delta*z*A/L^2)-2*z*rho*AR/L]/L",
        "weights": {"F": "-2-delta", "U": "-1-delta"},
        "x_formula": "caller-supplied x interval; no sqrt(rho) inference",
    }
    jets = {
        "Q0": Q0,
        "Qr": Qr,
        "rhoQrr": rhoQrr,
        "S": S,
        "P": P,
        "Np_formal": -S / 2 + lambda_j * P,
        "Nt_over_F0_anchor": Nt_over_F0_anchor,
        "Nz": Nz,
        "relative_F0": relative_f0,
        "x": x_j,
    }
    by_family = {
        family: {f"Z{order}": rows[f"{family}_Z{order}"] for order in _OUTPUT_ORDERS}
        for family in powers
    }
    return {
        **rows,
        "rows": rows,
        "sectors": by_family,
        "g4": by_family["g4"],
        "g5": by_family["g5"],
        "g6": by_family["g6"],
        "proof": proof,
        "metadata": proof,
        "jets": jets,
        "log_F0": log_f0_iv,
        "F0_log": log_f0_iv,
        "F0_log_interval": log_f0_iv,
    }


def first_collar_forcing_jets(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Alias for :func:`build_first_collar_forcing_jets`."""

    return build_first_collar_forcing_jets(*args, **kwargs)


def derive_first_collar_forcing_jets(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compatibility alias for source-leaf callers."""

    return build_first_collar_forcing_jets(*args, **kwargs)


def build_first_collar_forcing_z3(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compatibility alias retaining the historical ``z3`` operation name."""

    return build_first_collar_forcing_jets(*args, **kwargs)


def first_collar_forcing_z3(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Short alias for :func:`build_first_collar_forcing_jets`."""

    return build_first_collar_forcing_jets(*args, **kwargs)
