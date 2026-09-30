"""Signed-log flat-shape defect kernels.

For positive ``k``, ``m``, and ``T`` this module evaluates

.. math::

   I(k,m,B,T)=\int_0^T e^{-k\ell}\operatorname{expm1}
       (mB\,\sigma(\ell/T))\,d\ell,

where ``sigma`` is the flat logistic bump used by the finite reshape.  The
integral is signed when ``B`` is signed.  A saddle-centered finite quadrature
keeps the huge ``T`` and extremely small flat tail visible.  The derivative
with respect to ``B`` is evaluated as a separate positive signed-log integral;
it is never reconstructed by subtracting two nearby values.

This is a numerical finite-kernel adapter.  Its quadrature window and tails
are explicitly unenclosed, and it makes no field, cone, or global claim.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import mpmath as mp


def _mp(value: Any) -> mp.mpf:
    return value if isinstance(value, mp.mpf) else mp.mpf(str(value))


def _work_precision(T: mp.mpf, precision: int) -> int:
    """Retain enough decimal places to locate a ``T**(2/3)`` saddle."""

    digits = 0
    if T > 1:
        digits = max(0, int(mp.floor(mp.log10(T))) + 1)
    return max(80, int(precision), 320 + digits)


def _logistic_data(s: mp.mpf) -> tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    """Return ``q``, ``1-q``, logit ``D``, and ``D'`` for ``q=sigma(s)``."""

    if s <= 0:
        return mp.mpf(0), mp.mpf(1), mp.ninf, mp.inf
    if s >= 1:
        # At the exact right endpoint q is already saturated. Returning a
        # zero slope keeps (1-q) D' from becoming the indeterminate 0*inf;
        # the limiting log-density derivative is simply -k there.
        return mp.mpf(1), mp.mpf(0), mp.inf, mp.mpf(0)
    D = -1 / s**2 + 1 / (1 - s) ** 2
    if D <= 0:
        log_q = D - mp.log1p(mp.exp(D))
        log_one_minus_q = -mp.log1p(mp.exp(D))
    else:
        log_q = -mp.log1p(mp.exp(-D))
        log_one_minus_q = -D - mp.log1p(mp.exp(-D))
    q = mp.exp(log_q)
    one_minus_q = mp.exp(log_one_minus_q)
    D_prime = 2 / s**3 + 2 / (1 - s) ** 3
    return q, one_minus_q, D, D_prime


def _log_abs_expm1(value: mp.mpf) -> tuple[int, mp.mpf]:
    """Return sign and ``log(abs(expm1(value)))`` without cancellation."""

    value = _mp(value)
    if value == 0:
        return 0, mp.ninf
    sign = 1 if value > 0 else -1
    magnitude = abs(value)
    if magnitude < mp.mpf(".5"):
        # expm1(value)/value is positive on both sides of zero.
        ratio = mp.expm1(value) / value
        return sign, mp.log(magnitude) + mp.log(ratio)
    if value > 50:
        return sign, value + mp.log1p(-mp.exp(-value))
    if value < -50:
        return sign, mp.log1p(-mp.exp(value))
    return sign, mp.log(abs(mp.expm1(value)))


def _expm1_log_derivative_factor(value: mp.mpf) -> mp.mpf:
    """Return ``value*exp(value)/expm1(value)`` stably."""

    value = _mp(value)
    if value == 0:
        return mp.mpf(1)
    if abs(value) < mp.mpf(".5"):
        return value * mp.exp(value) / mp.expm1(value)
    if value > 50:
        return value / (1 - mp.exp(-value))
    if value < -50:
        return -value * mp.exp(value) / (1 - mp.exp(value))
    return value * mp.exp(value) / mp.expm1(value)


def _grid(T: mp.mpf, *, points_per_decade: int = 8) -> list[mp.mpf]:
    """Build a logarithmic candidate grid, retaining both endpoints."""

    if T <= 1:
        decades = 10
    else:
        decades = max(10, int(mp.ceil(mp.log10(T))) + 8)
    count = decades * points_per_decade
    values = {T}
    for j in range(count + 1):
        values.add(T * mp.power(10, -mp.mpf(j) / points_per_decade))
    # A linear set is useful for resolved T and catches a mode near the origin.
    if T <= 1000:
        for j in range(1, 33):
            values.add(T * mp.mpf(j) / 32)
    return sorted(value for value in values if 0 < value <= T)


def _bisect_mode(
    left: mp.mpf,
    right: mp.mpf,
    derivative: Callable[[mp.mpf], mp.mpf],
    *,
    iterations: int = 220,
) -> mp.mpf:
    """Bisect a positive-to-negative derivative crossing."""

    lo, hi = left, right
    for _ in range(iterations):
        mid = (lo + hi) / 2
        value = derivative(mid)
        if value > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def _locate_mode(
    T: mp.mpf,
    phi: Callable[[mp.mpf], mp.mpf],
    derivative: Callable[[mp.mpf], mp.mpf],
) -> tuple[mp.mpf, mp.mpf, str]:
    """Locate the dominant finite saddle or boundary maximum."""

    candidates = _grid(T)
    values: list[tuple[mp.mpf, mp.mpf]] = []
    derivatives: list[tuple[mp.mpf, mp.mpf]] = []
    for point in candidates:
        value = phi(point)
        values.append((value, point))
        derivatives.append((derivative(point), point))

    roots: list[mp.mpf] = []
    previous_d, previous_x = derivatives[0]
    for current_d, current_x in derivatives[1:]:
        if previous_d > 0 and current_d < 0:
            roots.append(_bisect_mode(previous_x, current_x, derivative))
        previous_d, previous_x = current_d, current_x

    options: list[tuple[mp.mpf, mp.mpf, str]] = [
        (value, point, "grid") for value, point in values if mp.isfinite(value)
    ]
    options.extend((phi(root), root, "saddle") for root in roots)
    if not options:
        raise ArithmeticError("flat-shape kernel has no finite log-density")
    best_value, best_point, method = max(options, key=lambda item: item[0])
    # Include the exact right endpoint when a narrow positive-B layer sits
    # close to T; the derivative formula is singular only at the left end.
    endpoint_value = phi(T)
    if endpoint_value > best_value:
        return T, endpoint_value, "right_endpoint"
    return best_point, best_value, method


def _curvature(
    mode: mp.mpf,
    T: mp.mpf,
    k: mp.mpf,
    derivative: Callable[[mp.mpf], mp.mpf],
) -> tuple[mp.mpf, mp.mpf]:
    """Return a negative curvature and its Gaussian width."""

    if mode <= 0 or mode >= T:
        width = mp.sqrt(max(mode, T * mp.mpf("1e-30")) / max(k, mp.mpf("1")))
        return -1 / width**2, width
    h = mode * mp.mpf("1e-6")
    h = min(h, mode / 3, (T - mode) / 3)
    if h <= 0:
        width = mp.sqrt(mode / max(k, mp.mpf("1")))
        return -1 / width**2, width
    second = (derivative(mode + h) - derivative(mode - h)) / (2 * h)
    if not mp.isfinite(second) or second >= 0:
        width = mp.sqrt(mode / max(k, mp.mpf("1")))
        return -1 / width**2, width
    return second, 1 / mp.sqrt(-second)


def _window_integral(
    *,
    T: mp.mpf,
    mode: mp.mpf,
    log_mode: mp.mpf,
    width: mp.mpf,
    phi: Callable[[mp.mpf], mp.mpf],
    order: int,
    window: mp.mpf,
) -> mp.mpf:
    """Integrate a normalized log-density over a finite saddle window.

    The quadrature coordinate is ``v = (ell - mode) / width``.  Splitting
    into unit-sized ``v`` panels avoids asking one Gauss rule to resolve the
    skewed tails of a broad resolved-scale saddle.  When the whole physical
    interval is at most 100 widths wide, it is inexpensive to integrate the
    whole interval and avoid an arbitrary tail cut.  The returned quantity is
    the physical ``d ell`` integral, with the common ``exp(log_mode)`` factor
    removed.
    """

    if width <= 0:
        return mp.mpf(0)
    total_span = T / width
    if total_span <= 100:
        lower = -mode / width
        upper = (T - mode) / width
    else:
        lower = max(-window, -mode / width)
        upper = min(window, (T - mode) / width)
    if upper <= lower:
        return mp.mpf(0)
    nodes, weights = mp.gauss_quadrature(int(order), "legendre")
    panel_count = max(1, int(mp.ceil(upper - lower)))
    total = mp.mpf(0)
    for panel in range(panel_count):
        panel_lower = lower + (upper - lower) * panel / panel_count
        panel_upper = lower + (upper - lower) * (panel + 1) / panel_count
        midpoint = (panel_lower + panel_upper) / 2
        half = (panel_upper - panel_lower) / 2
        for node, weight in zip(nodes, weights):
            v = midpoint + half * node
            point = mode + width * v
            delta = phi(point) - log_mode
            total += half * weight * mp.exp(delta)
    return width * total


@dataclass(frozen=True)
class FlatShapeKernel:
    """Reusable evaluator for one ``(k,m,T)`` flat-shape kernel."""

    k: mp.mpf
    m: mp.mpf
    T: mp.mpf
    precision: int = 260
    order: int = 32
    window: mp.mpf = mp.mpf(12)

    def __init__(
        self,
        k: Any,
        m: Any,
        T: Any,
        *,
        precision: int = 260,
        order: int = 32,
        window: Any = 12,
    ) -> None:
        if int(order) != order or int(order) < 8:
            raise ValueError("order must be an integer at least 8")
        if int(precision) != precision or int(precision) < 80:
            raise ValueError("precision must be an integer at least 80")
        precision_value = int(precision)
        with mp.workdps(max(80, precision_value)):
            k_value, m_value, T_value = _mp(k), _mp(m), _mp(T)
            window_value = _mp(window)
        if not all(mp.isfinite(value) for value in (k_value, m_value, T_value, window_value)):
            raise ValueError("k, m, T, and window must be finite")
        if not (k_value > 0 and m_value > 0 and T_value > 0):
            raise ValueError("k, m, and T must be positive")
        if window_value <= 0:
            raise ValueError("window must be positive")
        object.__setattr__(self, "k", k_value)
        object.__setattr__(self, "m", m_value)
        object.__setattr__(self, "T", T_value)
        object.__setattr__(self, "precision", precision_value)
        object.__setattr__(self, "order", int(order))
        object.__setattr__(self, "window", window_value)

    def evaluate(self, B: Any) -> dict[str, Any]:
        """Return signed-log value and independent first ``B`` sensitivity."""

        workdps = _work_precision(self.T, self.precision)
        with mp.workdps(workdps):
            k, m, T = self.k, self.m, self.T
            B_value = _mp(B)
            if not mp.isfinite(B_value):
                raise ValueError("B must be finite")

            def q_data(ell: mp.mpf):
                return _logistic_data(ell / T)

            def log_phi(ell: mp.mpf) -> mp.mpf:
                q, _, _, _ = q_data(ell)
                if q == 0:
                    return mp.ninf
                _, log_abs = _log_abs_expm1(m * B_value * q)
                return -k * ell + log_abs

            def phi_prime(ell: mp.mpf) -> mp.mpf:
                q, one_minus_q, _, D_prime = q_data(ell)
                if q == 0:
                    return mp.inf
                r = m * B_value * q
                factor = _expm1_log_derivative_factor(r)
                return -k + factor * one_minus_q * D_prime / T

            def sensitivity_phi(ell: mp.mpf) -> mp.mpf:
                q, _, _, _ = q_data(ell)
                if q == 0:
                    return mp.ninf
                return -k * ell + mp.log(m) + mp.log(q) + m * B_value * q

            def sensitivity_prime(ell: mp.mpf) -> mp.mpf:
                q, one_minus_q, _, D_prime = q_data(ell)
                if q == 0:
                    return mp.inf
                r = m * B_value * q
                return -k + one_minus_q * D_prime / T * (1 + r)

            sensitivity_mode, sensitivity_log, sensitivity_method = _locate_mode(
                T, sensitivity_phi, sensitivity_prime
            )
            _, sensitivity_width = _curvature(
                sensitivity_mode, T, k, sensitivity_prime
            )
            sensitivity_mass = _window_integral(
                T=T,
                mode=sensitivity_mode,
                log_mode=sensitivity_log,
                width=sensitivity_width,
                phi=sensitivity_phi,
                order=self.order,
                window=self.window,
            )
            sensitivity_log_abs = sensitivity_log + mp.log(sensitivity_mass)

            if B_value == 0:
                value_sign = 0
                value_log_abs = mp.ninf
                mode = sensitivity_mode
                log_mode = mp.ninf
                width = sensitivity_width
                method = "zero_B"
            else:
                value_sign, _ = _log_abs_expm1(m * B_value)
                mode, log_mode, method = _locate_mode(T, log_phi, phi_prime)
                _, width = _curvature(mode, T, k, phi_prime)
                mass = _window_integral(
                    T=T,
                    mode=mode,
                    log_mode=log_mode,
                    width=width,
                    phi=log_phi,
                    order=self.order,
                    window=self.window,
                )
                value_log_abs = log_mode + mp.log(mass)

            value = (
                mp.mpf(0)
                if value_sign == 0
                else value_sign * mp.exp(value_log_abs)
            )
            sensitivity = mp.exp(sensitivity_log_abs)
            return {
                "k": k,
                "m": m,
                "B": B_value,
                "T": T,
                "sign": value_sign,
                "log_abs": value_log_abs,
                "value": value,
                "derivative_sign": 1,
                "derivative_log_abs": sensitivity_log_abs,
                "derivative_value": sensitivity,
                "mode": mode,
                "width": width,
                "derivative_mode": sensitivity_mode,
                "derivative_width": sensitivity_width,
                "mode_method": method,
                "derivative_mode_method": sensitivity_method,
                "work_precision": workdps,
                "quadrature_order": self.order,
                "saddle_window": self.window,
                "quadrature_window_unenclosed": True,
                "tail_unenclosed": True,
                "signed_log_representation": True,
                "first_B_sensitivity_separate": True,
                "global_field_installed": False,
                "limitations": (
                    "Finite saddle-window quadrature; tails and quadrature remainder "
                    "are not enclosed."
                ),
            }


def evaluate_flat_shape_defect(
    k: Any,
    m: Any,
    B: Any,
    T: Any,
    precision: int = 260,
    order: int = 32,
    window: Any = 12,
) -> dict[str, Any]:
    """Evaluate one flat-shape kernel and its first ``B`` derivative."""

    return FlatShapeKernel(
        k,
        m,
        T,
        precision=precision,
        order=order,
        window=window,
    ).evaluate(B)


evaluate = evaluate_flat_shape_defect


__all__ = ["FlatShapeKernel", "evaluate_flat_shape_defect", "evaluate"]
