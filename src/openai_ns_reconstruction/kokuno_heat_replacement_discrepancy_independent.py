"""Independent numerical audit of the Kokuno heat-replacement discrepancy.

This validator is intentionally separate from the construction in
``kokuno_heat_replacement_discrepancy``.  It does not use the construction's
Gauss-Legendre / Gauss-Laguerre split or ``KokunoHeatExteriorProfile`` heat
quadrature.  Instead it evaluates

    (H(z)-1)/z = sum_{m>=1} (-1)^m (h)_m (h+1)_m z^(m-1) / m!

on the small-z tail used by the public default and performs adaptive quadrature
in the logarithmic coordinate ``y=log(X/X_tail)``.  The very slow I_sub tail is
integrated after the independent change of variable ``s=h*(y-3)``.

The purpose is validation only.  In particular, a disagreement reported here
must not be hidden by changing the fixed 1e-3 full PDE acceptance gate.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np
from scipy.integrate import quad

from .kokuno_heat_replacement_discrepancy import KokunoHeatReplacementDiscrepancy


@dataclass(frozen=True)
class IndependentHeatDiscrepancyAudit:
    """Direct log-space reference evaluator for the public default tail regime."""

    series_order: int = 12
    epsrel: float = 2.0e-10

    @staticmethod
    def _flat_step(value: float) -> float:
        if value <= 0.0:
            return 0.0
        if value >= 1.0:
            return 1.0
        left = math.exp(-1.0 / (value * value))
        right = math.exp(-1.0 / ((1.0 - value) ** 2))
        return left / (left + right)

    def _heat_minus_one_over_z(self, h: float, z: float) -> float:
        """Independent small-z series for ``(H(z)-1)/z``.

        The audited default has z <= 2/X_tail = 0.002.  The implementation
        refuses larger z rather than silently using the construction's heat
        quadrature path.
        """

        if z < 0.0 or z > 0.01:
            raise ValueError("independent heat series requires 0 <= z <= 0.01")
        q = -h * (1.0 + h)
        total = q
        for m in range(2, self.series_order + 1):
            q *= -z * (h + m - 1.0) * (h + m) / m
            total += q
        return total

    def discrepancy(
        self, model: KokunoHeatReplacementDiscrepancy, eta: float
    ) -> np.ndarray:
        """Return an independent ``(Delta C_p, Delta S, Delta I_sub)`` reference."""

        eta = float(eta)
        if not math.isfinite(eta) or abs(eta) > 1.0:
            raise ValueError("eta must lie in [-1,1]")
        if 2.0 / model.X_tail > 0.01:
            raise ValueError("X_tail is too small for the independent small-z audit")
        d = 1.0 - eta * eta
        if d == 0.0:
            return np.zeros(3, dtype=float)

        h = model.h
        A = 0.5 + h
        log_xt = math.log(model.X_tail)

        def chi(y: float) -> float:
            return self._flat_step((y - 0.2) / 0.3)

        def f_o(y: float) -> float:
            psi = 1.0 - self._flat_step((y - 1.0) / 2.0)
            return 1.0 - model.rho_o * psi

        def values_y(y: float) -> np.ndarray:
            L = log_xt + y
            z = 2.0 * d * math.exp(-L)
            ratio = self._heat_minus_one_over_z(h, z)
            ch = chi(y)
            multiplier = ch * z * ratio
            f = f_o(y)
            delta_square_factor = 2.0 * multiplier + multiplier * multiplier
            cp = (
                0.5
                * model.c_inf**2
                * math.exp(-2.0 * A * L)
                * f
                * f
                * delta_square_factor
            )
            stress = (
                -0.5
                * model.c_inf**2
                * math.exp((1.0 - 2.0 * A) * L)
                * f
                * f
                * delta_square_factor
            )
            # Combine z with the large X factor before evaluation.  This is
            # algebraically identical to sqrt(2X)*delta(E)*dX and remains
            # stable over the O(1/h) tail.
            i_sub = (
                2.0
                * math.sqrt(2.0)
                * model.c_inf
                * d
                * math.exp(-h * L)
                * f
                * ch
                * ratio
            )
            return np.asarray([cp, stress, i_sub], dtype=float)

        def values_tail_s(s: float) -> np.ndarray:
            # y>=3 has chi=f_o=1 exactly.  s=h*(y-3) removes the long
            # O(1/h) integration scale without using Laguerre quadrature.
            y = 3.0 + s / h
            L = log_xt + y
            z = 0.0 if L > 740.0 else 2.0 * d * math.exp(-L)
            ratio = self._heat_minus_one_over_z(h, z)
            multiplier = z * ratio
            delta_square_factor = 2.0 * multiplier + multiplier * multiplier
            cp = 0.5 * model.c_inf**2 * math.exp(-2.0 * A * L) * delta_square_factor
            stress = (
                -0.5
                * model.c_inf**2
                * math.exp((1.0 - 2.0 * A) * L)
                * delta_square_factor
            )
            i_sub = (
                2.0
                * math.sqrt(2.0)
                * model.c_inf
                * d
                * math.exp(-h * L)
                * ratio
            )
            return np.asarray([cp, stress, i_sub], dtype=float) / h

        result = np.zeros(3, dtype=float)
        cuts = (0.0, 0.2, 0.5, 1.0, 3.0)
        epsabs = (1.0e-18, 1.0e-14, 1.0e-10)
        for component in range(3):
            for lo, hi in zip(cuts[:-1], cuts[1:]):
                value, _ = quad(
                    lambda yy: float(values_y(yy)[component]),
                    lo,
                    hi,
                    epsabs=epsabs[component],
                    epsrel=self.epsrel,
                    limit=200,
                )
                result[component] += value
            value, _ = quad(
                lambda ss: float(values_tail_s(ss)[component]),
                0.0,
                np.inf,
                epsabs=epsabs[component],
                epsrel=self.epsrel,
                limit=200,
            )
            result[component] += value
        return result

    def first_order_tail_cp(
        self, model: KokunoHeatReplacementDiscrepancy, eta: float
    ) -> float:
        """Independent y>=3 first-order contribution to Delta C_p."""

        eta = float(eta)
        d = 1.0 - eta * eta
        if d <= 0.0:
            return 0.0
        h = model.h
        A = 0.5 + h
        log_xb = math.log(model.X_b)

        def integrand(s: float) -> float:
            L = log_xb + s
            z = 0.0 if L > 740.0 else 2.0 * d * math.exp(-L)
            ratio = self._heat_minus_one_over_z(h, z)
            return (
                2.0
                * model.c_inf**2
                * d
                * math.exp(-(2.0 * A + 1.0) * L)
                * ratio
            )

        value, _ = quad(integrand, 0.0, np.inf, epsabs=1.0e-20, epsrel=self.epsrel)
        return float(value)

    def compare_public(
        self,
        model: KokunoHeatReplacementDiscrepancy | None = None,
        etas: Any = (-0.8, -0.37, 0.0, 0.2, 0.8),
    ) -> dict[str, Any]:
        """Compare the public construction against the independent reference."""

        candidate = model or KokunoHeatReplacementDiscrepancy()
        eta_values = np.asarray(etas, dtype=float)
        public = np.asarray(candidate.discrepancy(eta_values), dtype=float)
        reference = np.stack(
            [self.discrepancy(candidate, float(eta)) for eta in eta_values], axis=0
        )
        absolute = np.abs(public - reference)
        relative = absolute / np.maximum(np.abs(reference), 1.0e-300)

        first_tail = np.asarray(
            [self.first_order_tail_cp(candidate, float(eta)) for eta in eta_values]
        )
        expected_public_minus_reference = -0.5 * first_tail
        observed_public_minus_reference = public[:, 0] - reference[:, 0]
        active = np.abs(expected_public_minus_reference) > 0.0
        factor_ratio = np.ones_like(first_tail)
        factor_ratio[active] = (
            observed_public_minus_reference[active]
            / expected_public_minus_reference[active]
        )

        return {
            "etas": eta_values.tolist(),
            "public": public.tolist(),
            "independent_reference": reference.tolist(),
            "max_relative_by_component": np.max(relative, axis=0).tolist(),
            "cp_first_order_tail": first_tail.tolist(),
            "cp_missing_half_factor_ratio": factor_ratio.tolist(),
            "diagnosis": (
                "public Delta C_p omits one half of the first-order y>=3 tail term"
            ),
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        }
