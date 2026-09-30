"""Executable Section 6 outer schedule from Lei--Ren Part I.

This module is a source-pinned velocity-schedule layer.  It implements the
reference branch and the uncorrected outer candidate in Section 6 of
Z. Lei and X. Ren, arXiv:2609.35406v1, equations (4.3), (4.6), (4.9),
(4.14), (4.18)--(4.22), and (6.1)--(6.5).  The regular swirl variable is

    F = U_theta / sqrt(2 R),

and the radial velocity is deliberately not synthesized here.  Pressure,
the five cumulative moments, the waiting-length matching, the Section 7
corrections, the Section 10 inner correction, stress cones, PDE residuals,
and scale recursion are separate tasks.

The paper leaves several absolute constants existential.  Callers therefore
provide numerical values for them, and the metadata records those choices as
uncertified assumptions.  Stage positions are stored as ``Decimal`` offsets
in ``log(R)``.  This is essential when ``13 / mu`` or ``30 log(1 / mu)`` is
larger than the spacing of binary64 numbers at the corresponding checkpoint.
The implementation does not materialize all physical radii as floats.

Source: https://arxiv.org/html/2609.35406v1
"""

from __future__ import annotations

from decimal import Decimal, InvalidOperation, ROUND_CEILING, localcontext
import math
import operator
from pathlib import Path
import sys
from typing import Any, Iterable

import numpy as np
from numpy.polynomial.legendre import leggauss


_ROOT = Path(__file__).resolve().parents[2]
_SRC = _ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from openai_ns_reconstruction.schedule_pressure import (  # noqa: E402
    outgoing_sigma,
    outgoing_sigma_derivative,
)
from openai_ns_reconstruction.heat_exterior import heat_factor  # noqa: E402


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_SECTION = "Sections 4.1, 6.1, equations (4.3), (4.6), (4.9), (4.14), (4.18)-(4.22), (6.1)-(6.5)"
DECIMAL_PRECISION = 128
MAX_DECIMAL_PRECISION = 4096
SIGMA_QUADRATURE_ORDER = 96
_LOG_FLOAT_MAX = math.log(np.finfo(float).max)
_D0 = Decimal(0)
_D1 = Decimal(1)
_D2 = Decimal(2)
_D_HALF = Decimal("0.5")


def _decimal(value: Any, name: str) -> Decimal:
    """Parse one finite Decimal without routing through binary64."""

    if isinstance(value, bool):
        raise TypeError(f"{name} must be a finite real number")
    try:
        result = value if isinstance(value, Decimal) else Decimal(str(value))
    except (InvalidOperation, TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not result.is_finite():
        raise ValueError(f"{name} must be finite")
    return result


def _decimal_float(value: Decimal, name: str) -> float:
    """Convert a Decimal to float only when it is representable."""

    result = float(value)
    if not math.isfinite(result):
        raise OverflowError(f"{name} is outside binary64 range")
    return result


def _dexp(value: Decimal, precision: int = DECIMAL_PRECISION) -> Decimal:
    with localcontext() as context:
        context.prec = int(precision)
        context.Emax = 999_999_999
        context.Emin = -999_999_999
        return value.exp()


def _dln(value: Decimal, precision: int = DECIMAL_PRECISION) -> Decimal:
    if value <= 0:
        raise ValueError("logarithm requires a positive Decimal")
    with localcontext() as context:
        context.prec = int(precision)
        context.Emax = 999_999_999
        context.Emin = -999_999_999
        return value.ln()


def _dlog1p(value: Decimal) -> Decimal:
    if value <= -1:
        raise ValueError("log1p requires an argument greater than -1")
    return _dln(_D1 + value)


def _decimal_add_exact(left: Decimal, right: Decimal) -> Decimal:
    """Add finite Decimal checkpoints without rounding the stored relation."""

    precision = len(left.as_tuple().digits) + len(right.as_tuple().digits) + 8
    with localcontext() as context:
        context.prec = precision
        context.Emax = 999_999_999
        context.Emin = -999_999_999
        return left + right


def _as_z(value: Any) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("Z must be finite with |Z| <= 1") from exc
    if not math.isfinite(result) or abs(result) > 1.0:
        raise ValueError("Z must be finite with |Z| <= 1")
    return result


def _sigma(value: Decimal | float) -> float:
    """Paper's flat step, exp(-1/t^2), evaluated from existing source code."""

    if isinstance(value, Decimal):
        if value <= 0:
            return 0.0
        if value >= 1:
            return 1.0
        value = float(value)
    return float(outgoing_sigma(float(value)))


def _sigma_prime(value: Decimal | float) -> float:
    if isinstance(value, Decimal):
        if value <= 0 or value >= 1:
            return 0.0
        value = float(value)
    return float(outgoing_sigma_derivative(float(value)))


def _flat_edge(value: Decimal | float) -> float:
    """The paper's ``mathfrak f(t) = exp(-1/t^2)`` for t > 0."""

    if isinstance(value, Decimal):
        if value <= 0:
            return 0.0
        value = float(value)
    value = float(value)
    if value <= 0.0:
        return 0.0
    if not math.isfinite(value):
        return 1.0
    exponent = -1.0 / (value * value)
    if exponent < -745.0:
        return 0.0
    return float(math.exp(exponent))


def _flat_edge_prime(value: Decimal | float) -> float:
    if isinstance(value, Decimal):
        if value <= 0:
            return 0.0
        value = float(value)
    value = float(value)
    if value <= 0.0 or not math.isfinite(value):
        return 0.0
    f = _flat_edge(value)
    if f == 0.0:
        return 0.0
    return float(2.0 * f / value**3)


def _exp_to_float(log_value: Decimal) -> float:
    """Return exp(log_value), with an explicit zero for representable underflow."""

    log_float = float(log_value)
    if log_float < math.log(np.nextafter(0.0, 1.0)):
        return 0.0
    if log_float > _LOG_FLOAT_MAX:
        raise OverflowError("profile value is outside binary64 range")
    return float(math.exp(log_float))


class PaperOuterSchedule:
    """Decimal-stable source schedule and its derived logarithmic checkpoints.

    ``logPstar`` and ``logRref`` are natural logarithms.  ``waiting_length``
    is the dimensionless ``tau`` in the paper.  The class accepts numerical
    values even when the paper's existential smallness conditions do not hold;
    ``metadata()`` reports those conditions instead of promoting them to a
    theorem claim.
    """

    _STAGE_ORDER = (
        "reference",
        "slope_transition_ref",
        "axial_turnoff",
        "slope_transition_mu",
        "power_buffer",
        "pulse_reserved",
        "z_flatten",
        "power_buffer_rel",
        "steep_transition_in",
        "steep_power",
        "steep_transition_out",
        "waiting",
        "heat_connection",
        "exact_heat",
    )

    def __init__(
        self,
        *,
        logPstar: Any,
        logRref: Any,
        delta: Any,
        Md: Any,
        c_mu: Any,
        c_delta: Any,
        c_epsilon: Any,
        waiting_length: Any = 0,
        sigma_quadrature_order: Any = SIGMA_QUADRATURE_ORDER,
    ) -> None:
        self.logPstar = _decimal(logPstar, "logPstar")
        self.logRref = _decimal(logRref, "logRref")
        self.delta = _decimal(delta, "delta")
        self.Md = _decimal(Md, "Md")
        self.c_mu = _decimal(c_mu, "c_mu")
        self.c_delta = _decimal(c_delta, "c_delta")
        self.c_epsilon = _decimal(c_epsilon, "c_epsilon")
        self.waiting_length = _decimal(waiting_length, "waiting_length")
        try:
            order = operator.index(sigma_quadrature_order)
        except TypeError as exc:
            raise TypeError("sigma_quadrature_order must be an integer") from exc
        if order <= 0:
            raise ValueError("sigma_quadrature_order must be positive")
        self.sigma_quadrature_order = int(order)

        if self.delta <= 0 or self.delta >= Decimal("0.005"):
            raise ValueError("delta must satisfy 0 < delta < 1/200")
        if self.Md <= 0:
            raise ValueError("Md must be positive")
        if self.c_mu <= 0 or self.c_delta <= 0 or self.c_epsilon <= 0:
            raise ValueError("c_mu, c_delta and c_epsilon must be positive")
        if self.waiting_length < 0:
            raise ValueError("waiting_length must be nonnegative")

        # The stage ``13 / mu`` can contain more decimal places than the
        # module default can retain.  Choose precision from the number of
        # decimal digits in ``exp(-log_mu)`` before constructing checkpoints.
        # This keeps additions such as ``y_v + 100`` exact for every accepted
        # input while retaining a finite resource bound for pathological
        # unbounded inputs.
        with localcontext() as context:
            context.prec = DECIMAL_PRECISION
            context.Emax = 999_999_999
            context.Emin = -999_999_999
            preliminary_log_mu = _dln(self.c_mu) - Decimal(4) * self.logPstar
            if preliminary_log_mu < 0:
                log10 = _dln(Decimal(10))
                exponent_digits = int(
                    ((-preliminary_log_mu) / log10).to_integral_value(
                        rounding=ROUND_CEILING
                    )
                )
            else:
                exponent_digits = 0
            self.decimal_precision = max(
                DECIMAL_PRECISION, exponent_digits + 80
            )
            if self.decimal_precision > MAX_DECIMAL_PRECISION:
                raise ValueError(
                    "log_mu requires more than the supported Decimal precision "
                    f"cap ({MAX_DECIMAL_PRECISION} digits)"
                )
            # Short alias for downstream schedule helpers that refer to the
            # active Decimal precision as ``schedule.precision``.
            self.precision = self.decimal_precision

        with localcontext() as context:
            context.prec = self.decimal_precision
            context.Emax = 999_999_999
            context.Emin = -999_999_999
            self.log_mu = _dln(self.c_mu, self.decimal_precision) - Decimal(4) * self.logPstar
            self.mu = _dexp(self.log_mu, self.decimal_precision)
            self.Td = _dexp(self.Md, self.decimal_precision) + Decimal(10)
            self.Tw = Decimal(60) * (-self.log_mu)
            self.Ts = Decimal(4) * _dln(Decimal(2) / self.delta, self.decimal_precision)
            self.Tf = Decimal(100)
            self.epsilon = self.c_epsilon * self.delta
            if self.epsilon >= 1:
                raise ValueError("epsilon=c_epsilon*delta must be less than one")

            self.y_h = Decimal(-5)
            self.y_d = Decimal(1) + self.Td
            self.y_w = self.y_d + Decimal(1)
            self.y_p = self.y_w + self.Tw
            self.y_v = self.y_p + Decimal(13) / self.mu
            self.y_f = _decimal_add_exact(self.y_v, self.Tf)
            self.y_rel = self.y_f - Decimal(30) * self.log_mu
            self.y_s = self.y_rel + Decimal(1)
            self.y_q = self.y_s + self.Ts
            self.y_t = self.y_q + Decimal(1)
            self.y_tail = self.y_t + self.waiting_length
            self.y_b = _decimal_add_exact(self.y_tail, Decimal(3))
            if self.y_f - self.y_v != Decimal(100):
                raise ArithmeticError("Decimal stage precision lost y_f - y_v = 100")
            if self.y_b - self.y_tail != Decimal(3):
                raise ArithmeticError("Decimal stage precision lost y_b - y_tail = 3")
            self.logR_h = self.logRref + self.y_h
            self.logR_d = self.logRref + self.y_d
            self.logR_w = self.logRref + self.y_w
            self.logR_p = self.logRref + self.y_p
            self.logR_v = self.logRref + self.y_v
            self.logR_f = self.logRref + self.y_f
            self.logR_rel = self.logRref + self.y_rel
            self.logR_s = self.logRref + self.y_s
            self.logR_q = self.logRref + self.y_q
            self.logR_t = self.logRref + self.y_t
            self.logR_tail = self.logRref + self.y_tail
            self.logR_b = self.logRref + self.y_b

            self._sigma_J1 = self._compute_sigma_primitive()
            self._sigma_J1_decimal = Decimal(str(self._sigma_J1))
            self._log_c_inf = (
                self._log_A(self.y_tail)
                + (Decimal(1) + self.delta) * self.logR_tail / Decimal(2)
                - _dln(
                    Decimal(2) * (Decimal(1) - self.epsilon),
                    self.decimal_precision,
                )
            )
            # Keep the logarithm as the authoritative representation.  A
            # source-invalid demo can make c_inf astronomically large (for
            # example when delta <= c_delta * mu is violated), and asking
            # Decimal.exp() to materialize it would request an impossible
            # number of digits.  Expose a Decimal value only in binary64's
            # finite logarithmic range, where callers can use it directly.
            log_min_float = Decimal(str(math.log(np.nextafter(0.0, 1.0))))
            log_max_float = Decimal(str(_LOG_FLOAT_MAX))
            if log_min_float <= self._log_c_inf <= log_max_float:
                self.c_inf_decimal = _dexp(
                    self._log_c_inf, self.decimal_precision
                )
            else:
                self.c_inf_decimal = None

        self._stage_bounds = self._make_stage_bounds()
        self._profile = PaperOuterProfile(self)

    def _compute_sigma_primitive(self) -> float:
        nodes, weights = leggauss(self.sigma_quadrature_order)
        x = (np.asarray(nodes, dtype=float) + 1.0) / 2.0
        w = np.asarray(weights, dtype=float) / 2.0
        return float(np.dot(w, np.asarray([outgoing_sigma(float(v)) for v in x])))

    def _J(self, argument: Decimal) -> Decimal:
        """Primitive of sigma using fixed [0,1] quadrature plus linear tail."""

        if argument <= 0:
            return _D0
        if argument >= 1:
            return self._sigma_J1_decimal + argument - _D1
        return Decimal(str(self._integrate_sigma_partial(float(argument))))

    def _integrate_sigma_partial(self, endpoint: float) -> float:
        nodes, weights = leggauss(self.sigma_quadrature_order)
        x = 0.5 * endpoint * (np.asarray(nodes, dtype=float) + 1.0)
        w = 0.5 * endpoint * np.asarray(weights, dtype=float)
        return float(np.dot(w, np.asarray([outgoing_sigma(float(v)) for v in x])))

    def _log_A(self, y: Decimal) -> Decimal:
        with localcontext() as context:
            context.prec = self.decimal_precision
            context.Emax = 999_999_999
            context.Emin = -999_999_999
            return (
                self.logPstar
                + Decimal("0.1") * y
                - Decimal("0.6") * self._J(y)
                - self.mu * self._J(y - self.y_d)
                - (Decimal(1) - self.mu) * self._J(y - self.y_rel)
                + (Decimal(1) - self.delta / Decimal(2))
                * self._J(y - self.y_rel - Decimal(1) - self.Ts)
            )

    def _slope_s(self, y: Decimal) -> Decimal:
        return (
            Decimal("0.1")
            - Decimal("0.6") * Decimal(str(_sigma(y)))
            - self.mu * Decimal(str(_sigma(y - self.y_d)))
            - (Decimal(1) - self.mu) * Decimal(str(_sigma(y - self.y_rel)))
            + (Decimal(1) - self.delta / Decimal(2))
            * Decimal(str(_sigma(y - self.y_rel - Decimal(1) - self.Ts)))
        )

    def _make_stage_bounds(self) -> dict[str, tuple[Decimal, Decimal | None]]:
        return {
            "reference": (Decimal("-Infinity"), _D0),
            "slope_transition_ref": (_D0, Decimal(1)),
            "axial_turnoff": (Decimal(1), self.y_d),
            "slope_transition_mu": (self.y_d, self.y_w),
            "power_buffer": (self.y_w, self.y_p),
            "pulse_reserved": (self.y_p, self.y_v),
            "z_flatten": (self.y_v, self.y_f),
            "power_buffer_rel": (self.y_f, self.y_rel),
            "steep_transition_in": (self.y_rel, self.y_s),
            "steep_power": (self.y_s, self.y_q),
            "steep_transition_out": (self.y_q, self.y_t),
            "waiting": (self.y_t, self.y_tail),
            "heat_connection": (self.y_tail, self.y_b),
            "exact_heat": (self.y_b, None),
        }

    @property
    def profile(self) -> "PaperOuterProfile":
        return self._profile

    @property
    def checkpoint_offsets(self) -> dict[str, Decimal]:
        return {
            "reference": _D0,
            "R_h": self.y_h,
            "R_d": self.y_d,
            "R_w": self.y_w,
            "R_p": self.y_p,
            "R_v": self.y_v,
            "R_f": self.y_f,
            "R_rel": self.y_rel,
            "R_s": self.y_s,
            "R_q": self.y_q,
            "R_t": self.y_t,
            "R_tail": self.y_tail,
            "R_b": self.y_b,
        }

    @property
    def checkpoint_log_radii(self) -> dict[str, Decimal]:
        return {
            name: self.logRref + offset
            for name, offset in self.checkpoint_offsets.items()
        }

    def _coerce_log_radius(self, value: Any) -> Decimal:
        return _decimal(value, "log_radius")

    def _coerce_radius_log(self, radius: Any) -> Decimal:
        value = _decimal(radius, "R")
        if value <= 0:
            raise ValueError("R must be positive")
        return _dln(value, self.decimal_precision)

    def stage_local_coordinate(self, log_radius: Any) -> dict[str, Any]:
        """Return the exact stage checkpoint and local log offset."""

        y = self._coerce_log_radius(log_radius) - self.logRref
        for name in self._STAGE_ORDER:
            start, end = self._stage_bounds[name]
            if y >= start and (end is None or y <= end):
                return {
                    "stage": name,
                    "checkpoint_offset": start,
                    "checkpoint_log_radius": self.logRref + start,
                    "offset": y - start,
                    "offset_string": str(y - start),
                }
        raise ArithmeticError("stage partition failed")

    def _axial_cutoff(self, y: Decimal) -> tuple[float, float]:
        t = y - Decimal(1)
        if t <= 0:
            return 1.0, 0.0
        cutoff_end = _dexp(self.Md, self.decimal_precision) - Decimal(1)
        if t >= cutoff_end:
            return 0.0, 0.0
        argument = _dln(_D1 + t, self.decimal_precision) / self.Md
        value = 1.0 - _sigma(argument)
        # d/dy [1 - sigma(log(1+t)/Md)].
        derivative = -_sigma_prime(argument) / (float(self.Md) * float(_D1 + t))
        return float(value), float(derivative)

    def _heat_info(self, log_radius: Decimal, z: float) -> dict[str, float | str]:
        d = max(0.0, 1.0 - z * z)
        if d == 0.0:
            inv_r = (
                math.exp(-float(log_radius))
                if float(log_radius) < 745.0
                else 0.0
            )
            return {
                "xi": 0.0,
                "H": 1.0,
                "H_prime": float(heat_factor(0.0, float(self.delta / 2), derivative=1)),
                "bound": 0.0,
                "method": "checked_heat_factor",
                "inv_R": inv_r,
            }
        if float(log_radius) <= _LOG_FLOAT_MAX:
            log_radius_float = float(log_radius)
            inv_R = math.exp(-log_radius_float) if log_radius_float < 745.0 else 0.0
            xi = 2.0 * d * inv_R
            if math.isfinite(xi):
                h = float(self.delta / 2)
                return {
                    "xi": xi,
                    "H": float(heat_factor(xi, h)),
                    "H_prime": float(heat_factor(xi, h, derivative=1)),
                    "bound": float(h * (1.0 + h) * xi),
                    "method": "checked_heat_factor",
                    "inv_R": inv_R,
                }
        with localcontext() as context:
            context.prec = self.decimal_precision
            context.Emax = 999_999_999
            context.Emin = -999_999_999
            inv_r_decimal = _dexp(-log_radius, self.decimal_precision)
            xi_decimal = Decimal(str(2.0 * d)) * inv_r_decimal
            h_decimal = self.delta / Decimal(2)
            bound_decimal = h_decimal * (Decimal(1) + h_decimal) * xi_decimal
            bound = float(bound_decimal) if bound_decimal < Decimal(str(np.finfo(float).max)) else float("inf")
        return {
            "xi": float(xi_decimal) if xi_decimal < Decimal(str(np.finfo(float).max)) else 0.0,
            "H": 1.0,
            "H_prime": 0.0,
            "bound": bound,
            "method": "H_to_1_declared_bound",
            "inv_R": float(inv_r_decimal) if inv_r_decimal < Decimal(str(np.finfo(float).max)) else 0.0,
        }

    def _tail_log_values(self, log_radius: Decimal, y: Decimal, z: float) -> dict[str, Any]:
        t = y - self.y_tail
        sigma_t = _sigma(t)
        sigma_t_prime = _sigma_prime(t)
        heat = self._heat_info(log_radius, z)
        heat_value = float(heat["H"])
        heat_prime = float(heat["H_prime"])
        flat_argument = (Decimal(3) - t) / Decimal(2)
        flat_value = _flat_edge(flat_argument)
        flat_prime = _flat_edge_prime(flat_argument)
        flat_factor = 1.0 - float(self.epsilon) * flat_value
        flat_factor_t = float(self.epsilon) * flat_prime / 2.0
        one_minus_epsilon = 1.0 - float(self.epsilon)
        K = ((1.0 - sigma_t) * one_minus_epsilon
             + sigma_t * heat_value * flat_factor)
        if not math.isfinite(K) or K <= 0.0:
            raise ArithmeticError("terminal heat interpolation lost positivity")
        xi = float(heat["xi"])
        heat_t = -xi * heat_prime
        K_t = (
            -sigma_t_prime * one_minus_epsilon
            + sigma_t_prime * heat_value * flat_factor
            + sigma_t * (heat_t * flat_factor + heat_value * flat_factor_t)
        )
        inv_r = float(heat["inv_R"])
        heat_z = heat_prime * (-4.0 * z * inv_r)
        logK_z = sigma_t * flat_factor * heat_z / K
        a = (Decimal(1) + self.delta) / Decimal(2)
        log_u = self._log_c_inf - a * log_radius + Decimal(str(math.log(K)))
        return {
            "log_u": log_u,
            "log_u_z": Decimal(str(logK_z)),
            "slope": -a + Decimal(str(K_t / K)),
            "heat_method": heat["method"],
            "heat_limit_bound": float(heat["bound"]),
            "heat_xi": xi,
            "K": K,
        }

    def _evaluate_log_radius(self, log_radius: Decimal, z: float) -> dict[str, Any]:
        z = _as_z(z)
        with localcontext() as context:
            context.prec = self.decimal_precision
            context.Emax = 999_999_999
            context.Emin = -999_999_999
            y = log_radius - self.logRref
            stage = self.stage_local_coordinate(log_radius)
            if y < 0:
                log_u = self.logPstar - _dln(
                    Decimal(str(1.0 + z * z)), self.decimal_precision
                ) + Decimal("0.1") * y
                log_u_z = Decimal(str(-2.0 * z / (1.0 + z * z)))
                slope = Decimal("0.1")
                uz = 4.0 * z
                uz_z = 4.0
                heat_method = None
                heat_limit_bound = 0.0
            elif y < self.y_tail:
                log_a = self._log_A(y)
                z_factor_log = math.log1p(z * z)
                flat_sigma = _sigma((y - self.y_v) / self.Tf)
                log_u = (
                    log_a
                    - Decimal(str(z_factor_log))
                    + Decimal(str(flat_sigma * (z_factor_log - math.log(2.0))))
                )
                log_u_z = Decimal(str(
                    -2.0 * z / (1.0 + z * z) * (1.0 - flat_sigma)
                ))
                slope = self._slope_s(y) + Decimal(str(
                    _sigma_prime((y - self.y_v) / self.Tf)
                    / float(self.Tf)
                    * (z_factor_log - math.log(2.0))
                ))
                uz_cutoff, _uz_cutoff_y = self._axial_cutoff(y)
                uz = 4.0 * z * uz_cutoff
                uz_z = 4.0 * uz_cutoff
                heat_method = None
                heat_limit_bound = 0.0
            else:
                tail = self._tail_log_values(log_radius, y, z)
                log_u = tail["log_u"]
                log_u_z = tail["log_u_z"]
                slope = tail["slope"]
                uz = 0.0
                uz_z = 0.0
                heat_method = tail["heat_method"]
                heat_limit_bound = tail["heat_limit_bound"]

            log_f = log_u - (Decimal(str(math.log(2.0))) + log_radius) / Decimal(2)
            f_value = _exp_to_float(log_f)
            u_value = _exp_to_float(log_u)
            log_u_z_float = float(log_u_z)
            return {
                "log_radius": log_radius,
                "y": y,
                "stage": stage["stage"],
                "checkpoint_offset": stage["checkpoint_offset"],
                "checkpoint_log_radius": stage["checkpoint_log_radius"],
                "stage_offset": stage["offset"],
                "stage_offset_string": stage["offset_string"],
                "log_angular_amplitude": log_u,
                "logF": log_f,
                "Utheta": u_value,
                "F": f_value,
                "Uz": float(uz),
                "Uz_Z": float(uz_z),
                "logarithmic_slope": slope,
                "logF_slope": slope - Decimal("0.5"),
                "dlogU_dZ": log_u_z,
                "dlogF_dZ": log_u_z,
                "Utheta_Z": float(u_value * log_u_z_float),
                "F_Z": float(f_value * log_u_z_float),
                "heat_method": heat_method,
                "heat_limit_bound": heat_limit_bound,
            }

    def at_log_radius(self, log_radius: Any, Z: Any) -> dict[str, Any]:
        """Evaluate the schedule from an exact/logarithmic radius."""

        return self._evaluate_log_radius(self._coerce_log_radius(log_radius), _as_z(Z))

    def at_radius(self, radius: Any, Z: Any) -> dict[str, Any]:
        """Evaluate at a positive physical radius without storing stage radii."""

        return self._evaluate_log_radius(self._coerce_radius_log(radius), _as_z(Z))

    def _metadata_conditions(self) -> dict[str, Any]:
        with localcontext() as context:
            context.prec = self.decimal_precision
            context.Emax = 999_999_999
            context.Emin = -999_999_999
            return {
                "Pstar_gt_exp_Td": bool(self.logPstar > self.Td),
                "mu_le_one_over_60": bool(self.mu <= Decimal(1) / Decimal(60)),
                "delta_le_c_delta_mu": bool(self.delta <= self.c_delta * self.mu),
                "epsilon_lt_one": bool(self.epsilon < 1),
                "waiting_length_nonnegative": bool(self.waiting_length >= 0),
                "Rref_threshold": "R_ref >= R_corr is not certifiable because R_corr is existential",
            }

    def metadata(self) -> dict[str, Any]:
        """Return JSON-ready source, schedule, and limitation metadata."""

        return {
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "source_section": SOURCE_SECTION,
            "api": "PaperOuterSchedule + PaperOuterProfile",
            "inputs": {
                "logPstar": str(self.logPstar),
                "logRref": str(self.logRref),
                "delta": str(self.delta),
                "Md": str(self.Md),
                "c_mu": str(self.c_mu),
                "c_delta": str(self.c_delta),
                "c_epsilon": str(self.c_epsilon),
                "waiting_length": str(self.waiting_length),
            },
            "derived": {
                "log_mu": str(self.log_mu),
                "mu": str(self.mu),
                "Td": str(self.Td),
                "Tw": str(self.Tw),
                "Ts": str(self.Ts),
                "Tf": str(self.Tf),
                "epsilon": str(self.epsilon),
                "log_c_inf": str(self._log_c_inf),
                "c_inf": (
                    str(self.c_inf_decimal)
                    if self.c_inf_decimal is not None
                    else None
                ),
                "c_inf_numeric_available": self.c_inf_decimal is not None,
                "sigma_primitive_J1": self._sigma_J1,
                "sigma_quadrature_order": self.sigma_quadrature_order,
            },
            "checkpoint_offsets": {name: str(value) for name, value in self.checkpoint_offsets.items()},
            "checkpoint_log_radii": {name: str(value) for name, value in self.checkpoint_log_radii.items()},
            "stage_order": list(self._STAGE_ORDER),
            "decimal_precision": self.decimal_precision,
            "radii_representation": "Decimal logarithmic checkpoints; no full float radius materialization",
            "heat_unrepresentable_policy": "H=1 with declared bound h*(1+h)*xi",
            "paper_absolute_constants_numeric": True,
            "theorem_conditions": self._metadata_conditions(),
            "claims": {
                "pressure_implemented": False,
                "five_moments_closed": False,
                "waiting_length_matched": False,
                "section_7_corrections_implemented": False,
                "section_10_inner_correction_implemented": False,
                "relaxed_cone_validated": False,
                "admissible_cone_validated": False,
                "pde_validated": False,
                "scale_recursion_validated": False,
            },
        }


class PaperOuterProfile:
    """Callable view of :class:`PaperOuterSchedule` at physical radii."""

    def __init__(self, schedule: PaperOuterSchedule) -> None:
        self.schedule = schedule

    def _evaluate(self, R: Any, Z: Any) -> dict[str, Any]:
        return self.schedule.at_radius(R, Z)

    def at_log_radius(self, log_radius: Any, Z: Any) -> dict[str, Any]:
        return self.schedule.at_log_radius(log_radius, Z)

    def stage_local_coordinate(self, log_radius: Any) -> dict[str, Any]:
        return self.schedule.stage_local_coordinate(log_radius)

    def Utheta(self, R: Any, Z: Any) -> float:
        return float(self._evaluate(R, Z)["Utheta"])

    U = Utheta

    def F(self, R: Any, Z: Any) -> float:
        return float(self._evaluate(R, Z)["F"])

    def Uz(self, R: Any, Z: Any) -> float:
        return float(self._evaluate(R, Z)["Uz"])

    def log_angular_amplitude(self, R: Any, Z: Any) -> Decimal:
        return self._evaluate(R, Z)["log_angular_amplitude"]

    def logF(self, R: Any, Z: Any) -> Decimal:
        return self._evaluate(R, Z)["logF"]

    def logarithmic_slope(self, R: Any, Z: Any) -> Decimal:
        return self._evaluate(R, Z)["logarithmic_slope"]

    def logF_slope(self, R: Any, Z: Any) -> Decimal:
        return self._evaluate(R, Z)["logF_slope"]

    def dlogU_dZ(self, R: Any, Z: Any) -> Decimal:
        return self._evaluate(R, Z)["dlogU_dZ"]

    def dlogF_dZ(self, R: Any, Z: Any) -> Decimal:
        return self._evaluate(R, Z)["dlogF_dZ"]

    def Utheta_Z(self, R: Any, Z: Any) -> float:
        return float(self._evaluate(R, Z)["Utheta_Z"])

    def F_Z(self, R: Any, Z: Any) -> float:
        return float(self._evaluate(R, Z)["F_Z"])

    def Uz_Z(self, R: Any, Z: Any) -> float:
        return float(self._evaluate(R, Z)["Uz_Z"])

    def metadata(self) -> dict[str, Any]:
        return self.schedule.metadata()


__all__ = [
    "PaperOuterSchedule",
    "PaperOuterProfile",
    "SOURCE",
    "SOURCE_VERSION",
    "SOURCE_SECTION",
]
