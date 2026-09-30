"""Continuous MP angular correction installed on one :class:`AngularCorrection`.

The legacy angular correction uses a float quadrature approximation to the
compact bump.  This adapter keeps the same signed-log receipt shape while
replacing the bump and its three weighted atoms by one shared
``ContinuousAxialBump``.  The small angular solve is performed after dividing
the two coefficients by the positive angular discrepancy.  This is essential
for the source-scale values, whose absolute coefficients can have logarithms
far below ordinary floating-point range.

The heat and pre-heat inputs intentionally retain the scope of the existing
``heat_defects`` module.  Heat constants are its central first-Taylor values;
the pre-heat term is a fixed-node MP re-evaluation of its positive difference
ODE.  Neither is an exact heat or pressure cancellation certificate.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import json
import math
import sys
import types
from types import SimpleNamespace
from typing import Any

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp


_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from lei_ren_part1_paper_angular_correction import AngularCorrection  # noqa: E402
from lei_ren_part1_paper_continuous_axial_basis import ContinuousAxialBump  # noqa: E402
from lei_ren_part1_paper_continuous_axial_pulse import ContinuousAxialPulse  # noqa: E402
from lei_ren_part1_paper_heat_defects import heat_defects  # noqa: E402


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_EQUATIONS = "(7.13), (7.17), (7.20)--(7.23)"
DEFAULT_BASIS_PRECISION = 100


def _mp(value: Any) -> mp.mpf:
    """Parse a value without routing an MP value through binary64."""

    if isinstance(value, mp.mpf):
        return value
    return mp.mpf(str(value))


def _z_string(value: Any, precision: int) -> str:
    """Return the cache key and receipt representation for an MP axial Z."""

    with mp.workdps(max(int(precision), 32)):
        z = _mp(value)
        if not mp.isfinite(z) or abs(z) > 1:
            raise ValueError("Z must be finite and lie in [-1,1]")
        return mp.nstr(z, max(int(precision), 32))


def _signed(value: Any, precision: int) -> dict[str, Any]:
    """Encode a possibly tiny signed MP value in the established receipt form."""

    with mp.workdps(max(int(precision), 32)):
        x = _mp(value)
        if x == 0:
            return {
                "sign": 0,
                "log_abs": None,
                "arbitrary_exponent_value": "0",
            }
        return {
            "sign": 1 if x > 0 else -1,
            "log_abs": mp.nstr(mp.log(abs(x)), int(precision)),
            "arbitrary_exponent_value": mp.nstr(x, int(precision)),
        }


def _signed_value(value: Any, precision: int) -> mp.mpf:
    """Decode both current full-value and older signed-log receipts."""

    if isinstance(value, mp.mpf):
        return value
    if not isinstance(value, dict):
        return _mp(value)
    exact = value.get("arbitrary_exponent_value")
    if exact is not None:
        return _mp(exact)
    sign = int(value.get("sign", 0))
    log_abs = value.get("log_abs")
    if sign == 0 or log_abs is None:
        return mp.mpf(0)
    with mp.workdps(max(int(precision), 32)):
        return (1 if sign > 0 else -1) * mp.exp(_mp(log_abs))


def _sigma_jet(x: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    """MP version of the shared flat step and its first derivative."""

    return ContinuousAxialPulse.sigma_pair(_mp(x))


class ContinuousAngularCorrectionProvider:
    """Live continuous bump atoms and coefficient derivatives.

    The provider owns the one ``ContinuousAxialBump`` instance used for
    pointwise values, coefficient weights, and optional correction integrals.
    ``correction`` is kept by identity; installation mutates that object so
    existing consumers do not silently acquire a second angular correction.
    """

    def __init__(
        self,
        correction: AngularCorrection,
        *,
        precision: int | None = None,
        basis_precision: int = DEFAULT_BASIS_PRECISION,
    ) -> None:
        if not isinstance(correction, AngularCorrection):
            raise TypeError("correction must be an AngularCorrection")
        self.correction = correction
        self.schedule = correction.schedule
        self.incoming = correction.incoming
        self.precision = max(
            int(precision if precision is not None else correction.precision),
            int(getattr(self.schedule, "decimal_precision", 0)),
            64,
        )
        self.basis_precision = max(int(basis_precision), 32)
        # The shared basis is only resolved at basis_precision.  Integrating
        # it at a 400 digit outer arithmetic precision would spend time on
        # digits that the normalizer and pointwise atom do not contain.
        self.integral_precision = self.basis_precision
        self.order = max(int(getattr(correction, "order", 128)), 16)
        self.basis = ContinuousAxialBump(
            precision=self.basis_precision,
            ell=".15",
        )
        self.ell = self.basis.ell
        # Decimal schedule parameters can have hundreds of significant
        # digits.  Parse them only inside the provider work context; parsing
        # an incoming Decimal at mpmath's default 15 digits would destroy the
        # small-mu separation between A and B.
        with mp.workdps(self.precision):
            self.mu = _mp(self.schedule.mu)

        # These are the only retained heat constants.  Their logs are copied
        # as strings first so a Decimal/MP exponent never becomes a float.
        central_heat = heat_defects(self.schedule, 0, order=self.order)
        if central_heat.get("log_r_H") is None or central_heat.get("log_s_H") is None:
            raise ArithmeticError("central heat defects did not provide positive constants")
        with mp.workdps(self.precision):
            self.C_r = mp.exp(_mp(central_heat["log_r_H"]))
            self.C_s = mp.exp(_mp(central_heat["log_s_H"]))

        self._heat0 = central_heat
        self._preheat_cache: dict[str, dict[str, Any]] = {}
        self._coefficient_cache: dict[str, dict[str, Any]] = {}
        self._weighted_cache: dict[tuple[str, int, str, str], mp.mpf] = {}
        self._prepare_preheat_baseline()
        self.weights = self._make_weights()
        # A few legacy callers ask the weight object for the energy atom.  A
        # SimpleNamespace is deliberate: all values remain MP and the old
        # attribute names remain available.
        self.weights.weighted_atom = self.weighted_atom
        self.full_outer_closed = False
        self.angular_pressure_and_moments_complete = False

    # ------------------------------------------------------------------
    # Shared continuous basis and weighted atoms
    # ------------------------------------------------------------------
    def _make_weights(self) -> SimpleNamespace:
        with mp.workdps(self.precision):
            normalizer = _mp(self.basis.normalizer)
            beta_inf = mp.exp(-1) / (self.ell * normalizer)
            beta_l2_squared = self.weighted_atom(0, power=2)
            A = self.weighted_atom(1 - self.mu, power=1)
            B = self.weighted_atom(-(1 + 2 * self.mu), power=1)
            D = self.weighted_atom(-(1 + 2 * self.mu), power=2)
            return SimpleNamespace(
                mu=self.mu,
                ell=self.ell,
                quadrature_order=self.order,
                normalization_integral=normalizer,
                normalizer=normalizer,
                beta_inf=beta_inf,
                beta_l2_squared=beta_l2_squared,
                A_mu=A,
                B_mu=B,
                D_mu=D,
                basis=self.basis,
                continuous_definition=(
                    "beta_ell(s)=exp(-1/(1-(s/ell)^2))/(ell*I0) "
                    "on |s|<ell"
                ),
            )

    def _integrate_weighted_interval(
        self,
        lam: mp.mpf,
        power: int,
        lower: mp.mpf,
        upper: mp.mpf,
    ) -> mp.mpf:
        """Directly integrate one clipped local interval.

        Flat-endpoint scaling is used when a requested interval lies wholly
        near either support endpoint.  This avoids asking adaptive quadrature
        to discover an exponentially small positive atom from an absolute
        stopping criterion.  No full-minus-partial subtraction is used.
        """

        if upper <= lower:
            return mp.mpf(0)
        ell = self.ell
        lo = max(_mp(lower), -ell)
        hi = min(_mp(upper), ell)
        if hi <= lo:
            return mp.mpf(0)
        norm = _mp(self.basis.normalizer)
        p = int(power)
        # Work in x=s/ell, where the raw flat bump has a fixed support.
        xlo = lo / ell
        xhi = hi / ell
        lamell = _mp(lam) * ell

        def log_raw(x: mp.mpf) -> mp.mpf:
            if abs(x) >= 1:
                return mp.ninf
            return -1 / (1 - x * x)

        def log_integrand(x: mp.mpf) -> mp.mpf:
            return lamell * x + p * log_raw(x)

        # Near a flat endpoint the maximum is at the endpoint for the
        # supported source parameters.  Keep the guard conservative; direct
        # quadrature is used for intervals reaching the central region.
        scaled = False
        if xhi <= mp.mpf("-.5") and abs(lamell) <= p:
            scale = log_integrand(xhi)
            scaled = True
        elif xlo >= mp.mpf(".5") and abs(lamell) <= p:
            scale = log_integrand(xlo)
            scaled = True
        else:
            scale = mp.mpf(0)

        points = [xlo]
        for point in (mp.mpf("-.5"), mp.mpf("0"), mp.mpf(".5")):
            if xlo < point < xhi:
                points.append(point)
        points.append(xhi)
        if scaled:
            integral = mp.quad(
                lambda x: mp.exp(log_integrand(x) - scale), points
            )
            integral *= mp.exp(scale)
        else:
            integral = mp.quad(lambda x: mp.exp(log_integrand(x)), points)
        return ell * integral / (ell * norm) ** p

    def weighted_atom(
        self,
        lam: Any,
        power: int = 1,
        lower: Any | None = None,
        upper: Any | None = None,
    ) -> mp.mpf:
        """Return ``∫ exp(lam*s) beta(s)**power ds`` on a clipped support.

        ``lower`` and ``upper`` are local coordinates relative to one bump
        center.  Omitting both means the complete support ``[-ell,ell]``.
        Bounds are MP values and the interval is integrated directly, which
        preserves tiny endpoint atoms without cancellation.
        """

        p = int(power)
        if p <= 0:
            raise ValueError("power must be a positive integer")
        with mp.workdps(self.integral_precision):
            if lower is None:
                lo = -self.ell
            else:
                lo = _mp(lower)
            if upper is None:
                hi = self.ell
            else:
                hi = _mp(upper)
            lam_mp = _mp(lam)
            key = (
                mp.nstr(lam_mp, self.integral_precision),
                p,
                mp.nstr(lo, self.integral_precision),
                mp.nstr(hi, self.integral_precision),
            )
            cached = self._weighted_cache.get(key)
            if cached is not None:
                return cached
            value = self._integrate_weighted_interval(lam_mp, p, lo, hi)
            self._weighted_cache[key] = value
            return value

    def _shifted_atom(
        self,
        lam: mp.mpf,
        center: mp.mpf,
        power: int = 1,
        lower_t: Any | None = None,
        upper_t: Any | None = None,
    ) -> mp.mpf:
        """Weighted atom for a bump centered at ``center`` in global ``t``."""

        # Parse global bounds in the outer context before the local integral
        # intentionally rounds to the declared basis precision.
        with mp.workdps(self.precision):
            lam_mp = _mp(lam)
            center_mp = _mp(center)
            lo = None if lower_t is None else _mp(lower_t) - center_mp
            hi = None if upper_t is None else _mp(upper_t) - center_mp
            return mp.exp(lam_mp * center_mp) * self.weighted_atom(
                lam_mp, power=power, lower=lo, upper=hi
            )

    # ------------------------------------------------------------------
    # Fixed baseline pre-heat model
    # ------------------------------------------------------------------
    def _prepare_preheat_baseline(self) -> None:
        """Cache one float ODE baseline and fixed quadrature node values."""

        a = float(1 - self.mu)
        duration = float(self.schedule.Tf)
        try:
            initial = float(self.incoming["incoming_flattening_X"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                "incoming must provide incoming_flattening_X for continuous preheat"
            ) from exc
        if not math.isfinite(initial) or duration <= 0:
            raise ValueError("incoming_flattening_X and Tf must define a finite positive baseline")

        base_log = -math.log(2.0)

        def rhs(t: float, state: list[float]) -> list[float]:
            with mp.workdps(self.precision):
                sigma, sigma_prime = _sigma_jet(mp.mpf(str(t / duration)))
                derivative = float(sigma_prime) / duration
            base_rate = a + base_log * derivative
            return [1.0 - base_rate * state[0]]

        result = solve_ivp(
            rhs,
            (0.0, duration),
            [initial],
            dense_output=True,
            rtol=1e-11,
            atol=1e-13,
            max_step=duration / 200.0,
        )
        if not result.success or result.sol is None:
            raise ArithmeticError("continuous preheat baseline ODE failed")

        nodes, quadrature_weights = leggauss(self.order)
        cached_nodes: list[tuple[mp.mpf, mp.mpf, mp.mpf]] = []
        with mp.workdps(self.precision):
            for node, weight in zip(nodes, quadrature_weights):
                x = (mp.mpf(str(float(node))) + 1) / 2
                t = _mp(self.schedule.Tf) * x
                baseline_float = float(result.sol(float(t))[0])
                if not math.isfinite(baseline_float):
                    raise ArithmeticError("preheat baseline ODE returned a nonfinite value")
                cached_nodes.append(
                    (
                        x,
                        mp.mpf(str(float(weight))) / 2,
                        mp.mpf(str(baseline_float)),
                    )
                )
        self._preheat_nodes = tuple(cached_nodes)
        self._preheat_ode = result
        self._preheat_initial = initial
        self._preheat_duration = duration

    def _preheat_at(self, z_exact: str) -> dict[str, Any]:
        cached = self._preheat_cache.get(z_exact)
        if cached is not None:
            return cached
        with mp.workdps(self.precision):
            z = _mp(z_exact)
            L = mp.log1p(z * z)
            dz = 2 * z / (1 + z * z)
            a = 1 - self.mu
            T = _mp(self.schedule.Tf)
            k = -mp.log(2) + L
            Q0 = mp.mpf(0)
            Q1 = mp.mpf(0)
            for x, weight, baseline in self._preheat_nodes:
                sigma, sigma_prime = _sigma_jet(x)
                W = mp.exp(-a * T * (1 - x) - k * (1 - sigma))
                atom = W * sigma_prime * baseline
                Q0 += weight * atom
                Q1 += weight * (1 - sigma) * atom
            contraction = mp.power(self.mu, 30 * (1 - self.mu))
            rpre = contraction * L * Q0
            rpreZ = contraction * dz * (Q0 - L * Q1)
            result = {
                "Z": float(z),
                "Z_exact": z_exact,
                "exact_zero": z == 0,
                "r_pre": _signed(rpre, self.precision),
                "r_pre_Z": _signed(rpreZ, self.precision),
                "log_r_pre": (
                    mp.nstr(mp.log(rpre), self.precision) if rpre > 0 else None
                ),
                "log_r_pre_Z": (
                    mp.nstr(mp.log(abs(rpreZ)), self.precision) if rpreZ else None
                ),
                "Q0": mp.nstr(Q0, self.precision),
                "Q1": mp.nstr(Q1, self.precision),
                "preheat_contraction": mp.nstr(contraction, self.precision),
                "baseline_ode_rtol": 1e-11,
                "baseline_quadrature_order": self.order,
                "baseline_evaluated_once_as_float": True,
                "sigma_and_Z_weights_MP": True,
                "scope": (
                    "Fixed-node MP re-evaluation of the positive pre-heat "
                    "difference ODE; baseline solve_ivp and quadrature are "
                    "numerical inputs, not an exact preheat identity"
                ),
            }
            self._preheat_cache[z_exact] = result
            return result

    def _heat_at(self, z: mp.mpf) -> dict[str, Any]:
        with mp.workdps(self.precision):
            one_minus_z2 = 1 - z * z
            rH = self.C_r * one_minus_z2
            sH = self.C_s * one_minus_z2
            rHZ = -2 * z * self.C_r
            sHZ = -2 * z * self.C_s
            return {
                "Z": float(z),
                "Z_exact": mp.nstr(z, self.precision),
                "exact_zero": one_minus_z2 == 0,
                "r_H": _signed(rH, self.precision),
                "s_H": _signed(sH, self.precision),
                "r_H_Z": _signed(rHZ, self.precision),
                "s_H_Z": _signed(sHZ, self.precision),
                "log_r_H": mp.nstr(mp.log(rH), self.precision) if rH > 0 else None,
                "log_s_H": mp.nstr(mp.log(sH), self.precision) if sH > 0 else None,
                "central_log_r_H": str(self._heat0["log_r_H"]),
                "central_log_s_H": str(self._heat0["log_s_H"]),
                "input_model": "r_H=C_r*(1-Z^2), s_H=C_s*(1-Z^2)",
                "Z_derivatives": "r_H_Z=-2*Z*C_r, s_H_Z=-2*Z*C_s",
                "source_heat_first_Taylor_constants": True,
                "quadrature_order": self.order,
                "scope": (
                    "Positive central first-Taylor heat constants inherited "
                    "from heat_defects; exact heat closure is open"
                ),
            }

    # ------------------------------------------------------------------
    # Coefficient solve and pointwise jets
    # ------------------------------------------------------------------
    def _compute_coefficients(self, z_exact: str) -> dict[str, Any]:
        with mp.workdps(self.precision):
            z = _mp(z_exact)
            pre = self._preheat_at(z_exact)
            heat = self._heat_at(z)
            rpre = _signed_value(pre["r_pre"], self.precision)
            rpreZ = _signed_value(pre["r_pre_Z"], self.precision)
            rH = _signed_value(heat["r_H"], self.precision)
            rHZ = _signed_value(heat["r_H_Z"], self.precision)
            s = _signed_value(heat["s_H"], self.precision)
            sZ = _signed_value(heat["s_H_Z"], self.precision)
            r = rpre + rH
            rZ = rpreZ + rHZ
            if r <= 0:
                raise ArithmeticError("actual angular target must remain positive")

            A = _mp(self.weights.A_mu)
            B = _mp(self.weights.B_mu)
            D = _mp(self.weights.D_mu)
            alpha = mp.exp(-2 * (1 - self.mu))
            beta = mp.exp(-2 * (1 + 2 * self.mu))
            rho0 = mp.exp(1 - self.mu) / A
            s_factor = mp.exp(-3 * (1 + 2 * self.mu)) / B
            q = D / (2 * B)

            # Solve for e_i=d_i/r.  The quadratic is the same normalized
            # small branch as the legacy implementation, but every operation
            # remains MP and the pressure target is never converted to float.
            sigma = s_factor * s / r
            c2 = q * r * (1 + beta * alpha * alpha)
            c1 = 1 - beta * alpha - 2 * q * beta * alpha * rho0 * r
            c0 = beta * rho0 + q * beta * rho0 * rho0 * r - sigma
            discriminant = c1 * c1 - 4 * c2 * c0
            scale = max(abs(c1 * c1), abs(4 * c2 * c0), mp.mpf(1))
            if c1 <= 0 or discriminant <= 0:
                raise ArithmeticError(
                    "source small-branch assumptions failed in continuous MP solve"
                )
            # This form avoids cancellation when c0 is tiny relative to c1.
            e1 = -2 * c0 / (c1 + mp.sqrt(discriminant))
            e2 = rho0 - alpha * e1
            d1 = r * e1
            d2 = r * e2

            # Differentiate the unscaled equations directly.  This keeps the
            # derivative of a signed tiny coefficient visible in its receipt.
            jacobian = mp.matrix(
                [
                    [alpha, mp.mpf(1)],
                    [1 + 2 * q * d1, beta * (1 + 2 * q * d2)],
                ]
            )
            rhs_Z = mp.matrix([rho0 * rZ, s_factor * sZ])
            dZ = mp.lu_solve(jacobian, rhs_Z)
            d1Z, d2Z = dZ[0], dZ[1]
            normalized_e1Z = (d1Z * r - d1 * rZ) / (r * r)
            normalized_e2Z = (d2Z * r - d2 * rZ) / (r * r)

            angular_residual = alpha * d1 + d2 - rho0 * r
            pressure_residual = (
                d1
                + beta * d2
                + q * (d1 * d1 + beta * d2 * d2)
                - s_factor * s
            )
            derivative_residual_1 = alpha * d1Z + d2Z - rho0 * rZ
            derivative_residual_2 = (
                (1 + 2 * q * d1) * d1Z
                + beta * (1 + 2 * q * d2) * d2Z
                - s_factor * sZ
            )
            # Match the legacy receipt's cancellation-free normalized
            # residuals.  Keep the direct residuals separately for diagnosis;
            # neither field is a claim of exact tiny-pressure cancellation.
            scaled_angular_residual = angular_residual / r
            scaled_pressure_residual = pressure_residual / r

            negative = max(-min(d1, d2), mp.mpf(0))
            if negative * _mp(self.weights.beta_inf) >= 1:
                raise ArithmeticError("continuous bump multiplier lost positivity")
            deficit = negative * _mp(self.weights.beta_inf)
            row: dict[str, Any] = {
                "Z": float(z),
                "Z_exact": z_exact,
                "r_pre": pre,
                "heat": heat,
                "r": _signed(r, self.precision),
                "r_Z": _signed(rZ, self.precision),
                "s": _signed(s, self.precision),
                "s_Z": _signed(sZ, self.precision),
                "log_r": mp.nstr(mp.log(r), self.precision),
                "log_s": mp.nstr(mp.log(s), self.precision) if s > 0 else None,
                "d1": _signed(d1, self.precision),
                "d2": _signed(d2, self.precision),
                "d1_Z": _signed(d1Z, self.precision),
                "d2_Z": _signed(d2Z, self.precision),
                "d1_over_r": mp.nstr(e1, self.precision),
                "d2_over_r": mp.nstr(e2, self.precision),
                "d1_over_r_Z": mp.nstr(normalized_e1Z, self.precision),
                "d2_over_r_Z": mp.nstr(normalized_e2Z, self.precision),
                "A_mu": mp.nstr(A, self.precision),
                "B_mu": mp.nstr(B, self.precision),
                "D_mu": mp.nstr(D, self.precision),
                "alpha": mp.nstr(alpha, self.precision),
                "beta": mp.nstr(beta, self.precision),
                "rho_over_r": mp.nstr(rho0, self.precision),
                "sigma": mp.nstr(sigma, self.precision),
                "q": mp.nstr(q, self.precision),
                "discriminant": mp.nstr(discriminant, self.precision),
                "scaled_algebraic_residuals": [
                    mp.nstr(scaled_angular_residual, self.precision),
                    mp.nstr(scaled_pressure_residual, self.precision),
                ],
                "absolute_algebraic_residuals": [
                    mp.nstr(angular_residual, self.precision),
                    mp.nstr(pressure_residual, self.precision),
                ],
                "derivative_algebraic_residuals": [
                    mp.nstr(derivative_residual_1, self.precision),
                    mp.nstr(derivative_residual_2, self.precision),
                ],
                "bump_multiplier_positive": True,
                "log_maximum_negative_multiplier_deficit": (
                    mp.nstr(mp.log(deficit), self.precision) if deficit else None
                ),
                "decimal_precision": self.precision,
                "basis_precision": self.basis_precision,
                "bump_quadrature_order": self.order,
                "continuous_bump_definition_shared": True,
                "input_status": (
                    "MP continuous bump weights with preheat baseline and "
                    "central first-Taylor heat inputs; exact pressure cancellation "
                    "is not claimed"
                ),
                "full_outer_closed": False,
            }
            self._coefficient_cache[z_exact] = row
            return row

    def coefficients(self, Z: Any) -> dict[str, Any]:
        """Return a cached coefficient receipt keyed by the full MP Z string."""

        key = _z_string(Z, self.precision)
        row = self._coefficient_cache.get(key)
        return row if row is not None else self._compute_coefficients(key)

    def _jet_coefficients(self, Z: Any, coefficients: dict[str, Any] | None) -> dict[str, Any]:
        if coefficients is None:
            return self.coefficients(Z)
        zkey = _z_string(Z, self.precision)
        supplied = coefficients.get("Z_exact") if isinstance(coefficients, dict) else None
        if supplied is None:
            supplied_z = coefficients.get("Z") if isinstance(coefficients, dict) else None
            if supplied_z is not None:
                supplied = _z_string(supplied_z, self.precision)
        if supplied is not None and supplied != zkey:
            raise ValueError("coefficient receipt uses a different exact Z")
        return coefficients

    def value_jet(
        self,
        t: Any,
        Z: Any,
        coefficients: dict[str, Any] | None = None,
    ) -> dict[str, mp.mpf]:
        """Return ``h,h_t,h_tt,h_Z,h_tZ`` from the shared continuous bumps."""

        with mp.workdps(self.precision):
            t_mp = _mp(t)
            # Check support before solving coefficients.  The outer field calls
            # this at many radii where the correction is identically zero.
            centers = (mp.mpf(-3), mp.mpf(-1))
            if all(abs(t_mp - center) >= self.ell for center in centers):
                zero = mp.mpf(0)
                return {"h": zero, "h_t": zero, "h_tt": zero, "h_Z": zero, "h_tZ": zero}
            row = self._jet_coefficients(Z, coefficients)
            d1 = _signed_value(row.get("d1"), self.precision)
            d2 = _signed_value(row.get("d2"), self.precision)
            d1Z = _signed_value(row.get("d1_Z", 0), self.precision)
            d2Z = _signed_value(row.get("d2_Z", 0), self.precision)
            values = [self.basis.values(t_mp - center) for center in centers]
            h = d1 * values[0]["beta"] + d2 * values[1]["beta"]
            ht = d1 * values[0]["beta_s"] + d2 * values[1]["beta_s"]
            htt = d1 * values[0]["beta_ss"] + d2 * values[1]["beta_ss"]
            hZ = d1Z * values[0]["beta"] + d2Z * values[1]["beta"]
            htZ = d1Z * values[0]["beta_s"] + d2Z * values[1]["beta_s"]
            return {"h": h, "h_t": ht, "h_tt": htt, "h_Z": hZ, "h_tZ": htZ}

    def relative_bump(
        self,
        t: Any,
        Z: Any,
        *,
        coefficients: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Return a signed-log receipt for ``h(t,Z)``.

        ``relative_correction`` remains separate from ``log1p`` because an
        ordinary ``1+h`` expression can round to exactly one for a tiny h.
        """

        with mp.workdps(self.precision):
            jet = self.value_jet(t, Z, coefficients=coefficients)
            h = jet["h"]
            receipt = _signed(h, self.precision)
            receipt.update(
                relative_correction=mp.nstr(h, self.precision),
                log1p_relative_correction=mp.nstr(mp.log1p(h), self.precision),
                t_exact=mp.nstr(_mp(t), self.precision),
                Z_exact=_z_string(Z, self.precision),
                continuous_bump_definition_shared=True,
            )
            return receipt

    # ------------------------------------------------------------------
    # Optional correction integrals used by pressure/energy consumers
    # ------------------------------------------------------------------
    def correction_integral(
        self,
        lam: Any,
        Z: Any,
        *,
        power: int = 1,
        lower: Any | None = None,
        upper: Any | None = None,
        coefficients: dict[str, Any] | None = None,
    ) -> mp.mpf:
        """Integrate ``exp(lam*t) h(t)^power`` over the requested t interval."""

        with mp.workdps(self.precision):
            row = self._jet_coefficients(Z, coefficients)
            p = int(power)
            if p == 1:
                d = [_signed_value(row["d1"], self.precision), _signed_value(row["d2"], self.precision)]
                return sum(
                    d_i * self._shifted_atom(_mp(lam), center, 1, lower, upper)
                    for d_i, center in zip(d, (mp.mpf(-3), mp.mpf(-1)))
                )
            if p == 2:
                d = [_signed_value(row["d1"], self.precision), _signed_value(row["d2"], self.precision)]
                return sum(
                    d_i * d_i * self._shifted_atom(_mp(lam), center, 2, lower, upper)
                    for d_i, center in zip(d, (mp.mpf(-3), mp.mpf(-1)))
                )
            raise ValueError("correction_integral power must be 1 or 2")

    def pressure_correction_integral(
        self,
        Z: Any,
        *,
        lower: Any | None = None,
        upper: Any | None = None,
        coefficients: dict[str, Any] | None = None,
    ) -> mp.mpf:
        """Return ``∫ exp(-(1+2mu)t) (h+h²/2) dt`` on a t interval."""

        with mp.workdps(self.precision):
            lam = -(1 + 2 * self.mu)
            return self.correction_integral(
                lam, Z, lower=lower, upper=upper, coefficients=coefficients
            ) + mp.mpf(".5") * self.correction_integral(
                lam,
                Z,
                power=2,
                lower=lower,
                upper=upper,
                coefficients=coefficients,
            )

    def independent_bump_replay(
        self,
        row: dict[str, Any],
        *,
        order: int = 192,
    ) -> dict[str, Any]:
        """Replay the two normalized equations with the continuous atoms.

        This keeps the old method's receipt keys while avoiding its float
        ``numpy.exp`` path, which cannot consume MP support coordinates.
        """

        # Retain ``order`` for receipt compatibility; MP atoms define the
        # actual replay and do not route through float quadrature.
        with mp.workdps(self.precision):
            e1 = _mp(row["d1_over_r"])
            e2 = _mp(row["d2_over_r"])
            lam_angular = 1 - self.mu
            lam_pressure = -(1 + 2 * self.mu)
            angular = sum(
                e * self._shifted_atom(lam_angular, center, 1)
                for e, center in zip((e1, e2), (mp.mpf(-3), mp.mpf(-1)))
            )
            pressure = sum(
                e * self._shifted_atom(lam_pressure, center, 1)
                for e, center in zip((e1, e2), (mp.mpf(-3), mp.mpf(-1)))
            )
            r = mp.exp(_mp(row["log_r"]))
            target_ratio = (
                mp.exp(_mp(row["log_s"]) - _mp(row["log_r"]))
                if row.get("log_s") is not None
                else mp.mpf(0)
            )
            return {
                "angular_increment_over_r": angular,
                "angular_scaled_difference_from_one": abs(angular - 1),
                "pressure_linear_increment_over_r": pressure,
                "log_pressure_target_over_r": (
                    mp.nstr(mp.log(target_ratio), 30) if target_ratio else None
                ),
                "log_quadratic_term_scale_over_r": row["log_r"],
                "quadrature_order": order,
                "continuous_atoms_replayed": True,
                "scope": (
                    "Independent continuous-atom replay; quadratic pressure "
                    "term remains at the signed tiny-target scale"
                ),
            }


def _install_method(correction: AngularCorrection, provider: ContinuousAngularCorrectionProvider) -> None:
    """Patch the existing object while retaining direct legacy methods."""

    correction._legacy_continuous_coefficients = correction.coefficients
    correction._legacy_continuous_relative_bump = correction.relative_bump

    @lru_cache(maxsize=128)
    def cached_coefficients(owner: AngularCorrection, Z: Any) -> dict[str, Any]:
        return provider.coefficients(Z)

    def relative_bump(
        owner: AngularCorrection,
        t: Any,
        Z: Any,
        *,
        coefficients: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        return provider.relative_bump(t, Z, coefficients=coefficients)

    def value_jet(
        owner: AngularCorrection,
        t: Any,
        Z: Any,
        coefficients: dict[str, Any] | None = None,
    ) -> dict[str, mp.mpf]:
        return provider.value_jet(t, Z, coefficients=coefficients)

    correction.coefficients = types.MethodType(cached_coefficients, correction)
    correction.relative_bump = types.MethodType(relative_bump, correction)
    correction.value_jet = types.MethodType(value_jet, correction)
    correction.independent_bump_replay = provider.independent_bump_replay
    # Replace only the shared weight view; the correction object itself keeps
    # its identity for every existing profile/tail/pressure consumer.
    correction.weights = provider.weights
    correction.weighted_atom = provider.weighted_atom
    correction._continuous_angular_provider = provider
    # CorrectedPressureAdapter's existing integration seam uses this shorter
    # legacy attribute name; point it at the same provider by identity.
    correction._continuous_provider = provider
    correction.full_outer_closed = False


def install_continuous_angular_correction(
    correction: AngularCorrection,
    precision: int | None = None,
    basis_precision: int = DEFAULT_BASIS_PRECISION,
) -> ContinuousAngularCorrectionProvider:
    """Install and return a continuous provider on the same correction object."""

    existing = getattr(correction, "_continuous_angular_provider", None)
    if existing is not None:
        return existing
    provider = ContinuousAngularCorrectionProvider(
        correction,
        precision=precision,
        basis_precision=basis_precision,
    )
    _install_method(correction, provider)
    return provider


def _fixture() -> tuple[AngularCorrection, ContinuousAngularCorrectionProvider]:
    """Build the small PaperOuterSchedule fixture used by the CLI checks."""

    from lei_ren_part1_paper_waiting_length import (
        apply_waiting_root,
        incoming_angular_ratio,
        solve_waiting_length,
    )
    from lei_ren_part1_paper_outer import PaperOuterSchedule

    base = PaperOuterSchedule(
        logPstar=14,
        logRref=10,
        delta="1e-32",
        Md=".5",
        c_mu=".001",
        c_delta=".001",
        c_epsilon=".01",
    )
    incoming = incoming_angular_ratio(base)
    matched = apply_waiting_root(base, solve_waiting_length(base, incoming["incoming_X"]))
    correction = AngularCorrection(matched, incoming, precision=160, order=64)
    provider = install_continuous_angular_correction(correction, precision=160, basis_precision=100)
    return correction, provider


def run() -> dict[str, Any]:
    """Run derivative, parity, support, and endpoint checks on the fixture."""

    correction, provider = _fixture()
    with mp.workdps(provider.precision):
        z = mp.mpf(".37")
        hstep = mp.mpf("1e-8")
        t = mp.mpf("-3") + mp.mpf(".041")
        row = provider.coefficients(z)
        center = provider.value_jet(t, z, coefficients=row)
        left = provider.value_jet(t - hstep, z, coefficients=row)
        right = provider.value_jet(t + hstep, z, coefficients=row)
        h_t_fd = (right["h"] - left["h"]) / (2 * hstep)
        h_tt_fd = (right["h_t"] - left["h_t"]) / (2 * hstep)
        zleft = provider.value_jet(t, z - hstep, coefficients=provider.coefficients(z - hstep))
        zright = provider.value_jet(t, z + hstep, coefficients=provider.coefficients(z + hstep))
        h_z_fd = (zright["h"] - zleft["h"]) / (2 * hstep)
        h_tz_fd = (zright["h_t"] - zleft["h_t"]) / (2 * hstep)

        normalized_derivative_checks = {
            "d1_over_r": mp.nstr(
                abs(
                    (
                        _signed_value(row["d1_Z"], provider.precision) * _signed_value(row["r"], provider.precision)
                        - _signed_value(row["d1"], provider.precision) * _signed_value(row["r_Z"], provider.precision)
                    )
                    / (_signed_value(row["r"], provider.precision) ** 2)
                    - _mp(row["d1_over_r_Z"])
                ),
                30,
            ),
            "d2_over_r": mp.nstr(
                abs(
                    (
                        _signed_value(row["d2_Z"], provider.precision) * _signed_value(row["r"], provider.precision)
                        - _signed_value(row["d2"], provider.precision) * _signed_value(row["r_Z"], provider.precision)
                    )
                    / (_signed_value(row["r"], provider.precision) ** 2)
                    - _mp(row["d2_over_r_Z"])
                ),
                30,
            ),
        }
        parity_plus = provider.coefficients(z)
        parity_minus = provider.coefficients(-z)
        parity = {
            "d1_even_difference": mp.nstr(
                abs(_signed_value(parity_plus["d1"], provider.precision) - _signed_value(parity_minus["d1"], provider.precision)),
                30,
            ),
            "d2_even_difference": mp.nstr(
                abs(_signed_value(parity_plus["d2"], provider.precision) - _signed_value(parity_minus["d2"], provider.precision)),
                30,
            ),
            "d1_Z_odd_difference": mp.nstr(
                abs(_signed_value(parity_plus["d1_Z"], provider.precision) + _signed_value(parity_minus["d1_Z"], provider.precision)),
                30,
            ),
            "d2_Z_odd_difference": mp.nstr(
                abs(_signed_value(parity_plus["d2_Z"], provider.precision) + _signed_value(parity_minus["d2_Z"], provider.precision)),
                30,
            ),
        }
        outside = provider.value_jet(mp.mpf("0"), z)
        endpoint_rows = {str(endpoint): provider.coefficients(endpoint) for endpoint in (-1, 0, 1)}
        endpoint_rpre_positive = all(
            _signed_value(endpoint_rows[str(endpoint)]["r_pre" ]["r_pre"], provider.precision) > 0
            for endpoint in (-1, 1)
        )
        result: dict[str, Any] = {
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "source_equations": SOURCE_EQUATIONS,
            "provider": {
                "precision": provider.precision,
                "basis_precision": provider.basis_precision,
                "integral_working_precision": provider.integral_precision,
                "continuous_bump_shared_by_identity": correction.weights.basis is provider.basis,
                "weights": {
                    key: mp.nstr(getattr(provider.weights, key), provider.precision)
                    for key in ("A_mu", "B_mu", "D_mu", "beta_l2_squared")
                },
                "full_outer_closed": False,
            },
            "jet_checks": {
                "h_t_relative_error": mp.nstr(abs(h_t_fd / center["h_t"] - 1), 30),
                "h_tt_relative_error": mp.nstr(abs(h_tt_fd / center["h_tt"] - 1), 30),
                "h_Z_relative_error": mp.nstr(abs(h_z_fd / center["h_Z"] - 1), 30),
                "h_tZ_relative_error": mp.nstr(abs(h_tz_fd / center["h_tZ"] - 1), 30),
                "normalized_coefficient_derivative_absolute_errors": normalized_derivative_checks,
            },
            "parity": parity,
            "support_check": {key: mp.nstr(value, 20) for key, value in outside.items()},
            "endpoint_check": {
                "r_pre_positive_at_pm1": endpoint_rpre_positive,
                "heat_r_zero_at_pm1": all(
                    endpoint_rows[str(endpoint)]["heat"]["exact_zero"] for endpoint in (-1, 1)
                ),
                "heat_r_positive_at_zero": _signed_value(endpoint_rows["0"]["heat"]["r_H"], provider.precision) > 0,
            },
            "coefficients_cache_keys_are_exact_strings": all(
                isinstance(endpoint_rows[str(endpoint)]["Z_exact"], str) for endpoint in (-1, 0, 1)
            ),
            "signed_tiny_correction_preserved": True,
            "pressure_cancellation_exact": False,
            "quadrature_enclosure_certified": False,
            "scope": (
                "Continuous angular bump values, jets, weighted atoms and MP "
                "small-branch coefficient receipts; heat Taylor and fixed-node "
                "preheat inputs remain numerical inherited approximations."
            ),
        }
    output = _HERE / "lei_ren_part1_paper_continuous_angular_correction.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"jet_checks": result["jet_checks"], "endpoint_check": result["endpoint_check"]}), flush=True)
    return result


if __name__ == "__main__":
    run()
