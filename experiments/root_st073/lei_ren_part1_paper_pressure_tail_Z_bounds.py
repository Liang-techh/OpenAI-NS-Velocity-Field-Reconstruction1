"""Z-derivative receipts for the same-profile pressure tail.

The source angular factor on the finite flattening interval is

    U_theta = E(y) (1 + Z**2)**(-1) * ((1 + Z**2)/2)**sigma(y).

Consequently the uncorrected pressure integrand has the pointwise envelope
``abs(d_Z(U_theta**2/2)) <= 4*abs(Z)/(1+Z**2) * U_theta**2 <= 2*U_theta**2``.
After flattening, the nominal source is Z-independent until the heat
connection.  The heat factor gives a separate analytic envelope using
``|H'(xi)| <= h*(1+h)``.  The angular bump correction is retained as a
separate signed term.  Its coefficient derivatives are not silently set to
zero, so the returned corrected-tail derivative remains explicitly
unresolved.

All pressure values are reported in both normalized ``P/Pstar**2`` and
physical units.  Gauss quadrature refinement is reported separately from the
pointwise analytic envelopes; this module does not claim a global PDE,
regular-core, cone, or theorem certification.

Source: Lei--Ren, arXiv:2609.35406v1, Sections 4.14--4.25 and 7.17--7.23.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys
from typing import Any, Callable

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SRC = ROOT / "src"
for path in (SRC, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_paper_axis_pressure_jets import (  # noqa: E402
    AxisPressureJets,
    _from_signed_log,
    _signed_log,
)
from lei_ren_part1_paper_corrected_profile import CorrectedSourceProfile  # noqa: E402


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"


def _relative(value: mp.mpf, scale: mp.mpf, precision: int) -> dict[str, Any]:
    return {
        "absolute": _signed_log(value, precision),
        "relative": mp.nstr(abs(value / scale), precision) if scale else None,
    }


def _finite_z(value: Any) -> tuple[mp.mpf, float]:
    try:
        z = mp.mpf(str(value))
    except Exception as exc:  # pragma: no cover - defensive input boundary
        raise ValueError("Z must be a finite real number in [-1,1]") from exc
    if not mp.isfinite(z) or abs(z) > 1:
        raise ValueError("Z must be a finite real number in [-1,1]")
    return z, float(z)


class PressureTailZBounds:
    """Analytic baseline-tail derivative envelope plus unresolved bump term.

    ``jets`` must be an ``AxisPressureJets`` instance constructed from the
    same ``CorrectedSourceProfile`` used by the pressure adapter.  The bound
    integrates a pointwise analytic envelope in logarithmic radius, so no
    astronomical physical radius is materialized.
    """

    def __init__(self, jets: AxisPressureJets, *, quadrature_order: int = 128) -> None:
        if not isinstance(jets, AxisPressureJets):
            raise TypeError("jets must be an AxisPressureJets instance")
        if int(quadrature_order) != quadrature_order or quadrature_order < 16:
            raise ValueError("quadrature_order must be an integer at least 16")
        self.jets = jets
        self.profile = jets.profile
        self.schedule = jets.schedule
        self.pressure = jets.pressure
        self.precision = int(jets.precision)
        self.quadrature_order = int(quadrature_order)
        self._bound_cache: dict[tuple[str, int], dict[str, Any]] = {}

    def _integrate_y(
        self,
        left: Decimal,
        right: Decimal,
        function: Callable[[Decimal], mp.mpf],
        order: int,
    ) -> mp.mpf:
        """Integrate a local y-coordinate without converting checkpoints to float."""

        if right <= left:
            return mp.mpf("0")
        nodes, weights = leggauss(int(order))
        with localcontext() as context:
            context.prec = max(int(self.schedule.decimal_precision), self.precision)
            half = (right - left) / Decimal(2)
            midpoint = (right + left) / Decimal(2)
            total = mp.mpf("0")
            for node, weight in zip(nodes, weights):
                y = midpoint + half * Decimal(str(float(node)))
                total += mp.mpf(str(float(weight))) * function(y)
        return mp.mpf(str(half)) * total

    def _log_u_sq_over_pstar_sq(self, log_radius: Decimal, z_for_schedule: float) -> mp.mpf:
        row = self.schedule.at_log_radius(log_radius, z_for_schedule)
        with mp.workdps(self.precision):
            return mp.exp(
                2
                * (
                    mp.mpf(str(row["log_angular_amplitude"]))
                    - mp.mpf(str(self.schedule.logPstar))
                )
            )

    def _flatten_bound(self, z: mp.mpf, *, order: int) -> mp.mpf:
        """Bound the finite flattening interval in normalized pressure units."""

        z_factor = 4 * abs(z) / (1 + z * z)
        if z_factor == 0:
            return mp.mpf("0")

        def integrand(y: Decimal) -> mp.mpf:
            with localcontext() as context:
                context.prec = max(int(self.schedule.decimal_precision), self.precision)
                log_radius = self.schedule.logRref + y
            # U(Z)^2 <= U(0)^2 on 0 <= sigma <= 1.  The requested envelope
            # is 4|Z|/(1+Z^2) times U^2, so using U(0)^2 is conservative.
            return z_factor * self._log_u_sq_over_pstar_sq(log_radius, 0.0)

        return self._integrate_y(self.schedule.y_v, self.schedule.y_f, integrand, order)

    def _heat_connection_bound(self, z: mp.mpf, *, order: int) -> mp.mpf:
        """Bound Z dependence of K_delta on [y_tail,y_b]."""

        h = mp.mpf(str(self.schedule.delta)) / 2
        a1 = h * (1 + h)
        z_factor = 4 * abs(z) * a1
        if z_factor == 0:
            return mp.mpf("0")

        def integrand(y: Decimal) -> mp.mpf:
            with localcontext() as context:
                context.prec = max(int(self.schedule.decimal_precision), self.precision)
                log_radius = self.schedule.logRref + y
                log_A = self.schedule._log_A(y)
            # K <= 1 and |d_Z H(2(1-Z^2)/R)| <= 4|Z| a1/R.
            exponent = (
                2 * (mp.mpf(str(log_A)) - mp.mpf(str(self.schedule.logPstar)))
                - mp.mpf(str(log_radius))
            )
            return z_factor * mp.exp(exponent)

        return self._integrate_y(
            self.schedule.y_tail,
            self.schedule.y_b,
            integrand,
            order,
        )

    def _exact_heat_bound(self, z: mp.mpf) -> mp.mpf:
        """Analytic heat-exterior derivative envelope in normalized units."""

        h = mp.mpf(str(self.schedule.delta)) / 2
        a1 = h * (1 + h)
        z_factor = 4 * abs(z) * a1 / (2 + mp.mpf(str(self.schedule.delta)))
        if z_factor == 0:
            return mp.mpf("0")
        exponent = (
            2
            * (
                mp.mpf(str(self.schedule._log_c_inf))
                - mp.mpf(str(self.schedule.logPstar))
            )
            - (2 + mp.mpf(str(self.schedule.delta)))
            * mp.mpf(str(self.schedule.logR_b))
        )
        return z_factor * mp.exp(exponent)

    def baseline_derivative_bound(
        self,
        Z: Any,
        *,
        quadrature_order: int | None = None,
    ) -> dict[str, Any]:
        """Return the analytic envelope for ``|d_Z P_baseline(Rv,Z)|``.

        The returned finite Gauss value is accompanied by 64/128-style order
        refinement in the executable receipt.  The pointwise heat and
        flattening envelopes are analytic; quadrature truncation is not
        silently called a rigorous global error certificate.
        """

        z, z_float = _finite_z(Z)
        order = self.quadrature_order if quadrature_order is None else int(quadrature_order)
        if order < 16:
            raise ValueError("quadrature_order must be at least 16")
        key = (mp.nstr(z, self.precision), order)
        if key in self._bound_cache:
            return self._bound_cache[key]
        with mp.workdps(self.precision):
            flatten = self._flatten_bound(z, order=order)
            heat_connection = self._heat_connection_bound(z, order=order)
            exact_heat = self._exact_heat_bound(z)
            total = flatten + heat_connection + exact_heat
            pstar_squared = mp.exp(2 * mp.mpf(str(self.schedule.logPstar)))
            result = {
                "Z": z_float,
                "quadrature_order": order,
                "pressure_units": "P/Pstar^2",
                "flattening_bound_normalized": _signed_log(flatten, self.precision),
                "heat_connection_bound_normalized": _signed_log(
                    heat_connection, self.precision
                ),
                "exact_heat_bound_normalized": _signed_log(exact_heat, self.precision),
                "total_bound_normalized": _signed_log(total, self.precision),
                "flattening_bound_physical": _signed_log(
                    flatten * pstar_squared, self.precision
                ),
                "heat_connection_bound_physical": _signed_log(
                    heat_connection * pstar_squared, self.precision
                ),
                "exact_heat_bound_physical": _signed_log(
                    exact_heat * pstar_squared, self.precision
                ),
                "total_bound_physical": _signed_log(
                    total * pstar_squared, self.precision
                ),
                "Pstar_squared": mp.nstr(pstar_squared, self.precision),
                "pointwise_formula": (
                    "flatten: |d_Z(U^2/2)| <= 4|Z|/(1+Z^2) U^2 <= 2U^2; "
                    "heat: |H'| <= h(1+h), K<=1"
                ),
                "nominal_Z_independent_interval": "[y_f,y_tail]",
                "quadrature_status": (
                    "analytic pointwise envelope; finite Gauss evaluation with "
                    "order refinement, no global quadrature proof"
                ),
            }
        self._bound_cache[key] = result
        return result

    def correction_derivative_receipt(
        self,
        Z: Any,
        *,
        quadrature_order: int | None = None,
    ) -> dict[str, Any]:
        """Expose the corrected bump term and its unresolved coefficient derivative.

        At ``Rv`` the source amplitude ``E_rel`` is already flattened, so its
        Z derivative is exactly zero.  If ``L_j`` and ``Q_j`` are the bump
        linear and quadratic pressure weights, then

        ``d_Z I_b = d1_Z (L1+d1 Q1) + d2_Z (L2+d2 Q2)``.

        The current coefficient receipt does not provide certified ``d_j_Z``;
        this method returns the exact algebraic form and leaves that term
        unresolved instead of estimating it by float finite differences.
        """

        z, z_float = _finite_z(Z)
        order = self.quadrature_order if quadrature_order is None else int(quadrature_order)
        if order < 16:
            raise ValueError("quadrature_order must be at least 16")
        with mp.workdps(self.precision):
            coefficient = self.profile.angular.coefficients(z_float)
            increment = self.pressure.pressure_increment(
                z_float, quadrature_order=order
            )
            d1 = _from_signed_log(coefficient["d1"])
            d2 = _from_signed_log(coefficient["d2"])
            rows = increment["pressure_term_receipt"]["rows"]
            linear = [mp.mpf(str(row["linear_integral"])) for row in rows]
            quadratic = [mp.mpf(str(row["quadratic_integral"])) for row in rows]
            basis_1 = linear[0] + d1 * quadratic[0]
            basis_2 = linear[1] + d2 * quadratic[1]
            bump_integral = _from_signed_log(
                increment["delta_p_bump_over_Erel_squared"]
            )
            log_u_rel = self.schedule.at_log_radius(
                self.schedule.logR_rel, 0.0
            )["log_angular_amplitude"]
            scale = mp.exp(
                2
                * (
                    mp.mpf(str(log_u_rel))
                    - mp.mpf(str(self.schedule.logPstar))
                )
            )
            correction_normalized = -scale * bump_integral
            pstar_squared = mp.exp(2 * mp.mpf(str(self.schedule.logPstar)))
            return {
                "Z": z_float,
                "quadrature_order": order,
                "pressure_units": "P/Pstar^2",
                "correction_P_Rv_normalized": _signed_log(
                    correction_normalized, self.precision
                ),
                "correction_P_Rv_physical": _signed_log(
                    correction_normalized * pstar_squared, self.precision
                ),
                "d1": coefficient["d1"],
                "d2": coefficient["d2"],
                "linear_weights": [mp.nstr(value, self.precision) for value in linear],
                "quadratic_weights": [
                    mp.nstr(value, self.precision) for value in quadratic
                ],
                "formal_dI_b_dZ": (
                    "d1_Z*(L1+d1*Q1) + d2_Z*(L2+d2*Q2)"
                ),
                "basis_for_unknown_d1_Z": _signed_log(basis_1, self.precision),
                "basis_for_unknown_d2_Z": _signed_log(basis_2, self.precision),
                "formal_dP_b_dZ": (
                    "-E_rel^2/Pstar^2 * [d1_Z*(L1+d1 Q1)+d2_Z*(L2+d2 Q2)]"
                ),
                "E_rel_Z_status": "exactly zero after source flattening",
                "coefficient_derivative_status": (
                    "unresolved: current AngularCorrection receipt has d1,d2 "
                    "values but no certified Z derivatives"
                ),
                "finite_difference_used": False,
                "full_corrected_tail_C1_bound": False,
                "pressure_increment_receipt": increment,
            }

    def receipt(self, Z: Any, *, quadrature_order: int | None = None) -> dict[str, Any]:
        """Combine baseline bound, same-profile tail value, and bump status."""

        order = self.quadrature_order if quadrature_order is None else int(quadrature_order)
        bound = self.baseline_derivative_bound(Z, quadrature_order=order)
        correction = self.correction_derivative_receipt(
            Z, quadrature_order=order
        )
        tail = self.jets.pressure_receipt(
            Z, quadrature_order=order, include_direct_reference=False
        )
        with mp.workdps(self.precision):
            baseline_tail = _from_signed_log(
                tail["future_tail"]["baseline_pressure_over_Pstar_squared"]
            )
            bound_value = _from_signed_log(bound["total_bound_normalized"])
            ratio = abs(bound_value / baseline_tail) if baseline_tail else mp.mpf("0")
        return {
            "Z": float(Z),
            "quadrature_order": order,
            "same_profile_tail": tail,
            "uncorrected_baseline_derivative_bound": bound,
            "corrected_bump_derivative": correction,
            "baseline_tail_bound_to_value_ratio": mp.nstr(ratio, self.precision),
            "total_derivative_status": (
                "baseline uncorrected C1 envelope supplied; corrected angular "
                "bump coefficient derivative unresolved"
            ),
            "finite_difference_used": False,
            "scope": (
                "Future-pressure Z derivative diagnostic for the same source "
                "profile. It does not certify the corrected tail derivative, "
                "PDE, cone, global energy, regular core, or theorem conditions."
            ),
        }


def _make_true_source_profile(
    *, precision: int = 160, order: int = 128
) -> tuple[CorrectedSourceProfile, mp.mpf, mp.mpf]:
    """Build the explicit Lambda=10^36 source candidate without materializing R."""

    with mp.workdps(precision):
        Lambda = mp.mpf("1e36")
        log_cstar = 2 * mp.log(Lambda)
        log_pstar = mp.mpf("14")
        log_rref = mp.log(110) + 10 * (log_cstar + log_pstar)
    from lei_ren_part1_paper_outer import PaperOuterSchedule
    from lei_ren_part1_paper_axial_energy_tail import build_default_tail

    schedule = PaperOuterSchedule(
        logPstar=mp.nstr(log_pstar, precision),
        logRref=mp.nstr(log_rref, precision),
        delta="1e-32",
        Md=".5",
        c_mu=".001",
        c_delta=".001",
        c_epsilon=".01",
    )
    tail = build_default_tail(
        precision=precision,
        quadrature_order=order,
        schedule=schedule,
        match_waiting=True,
    )
    profile = CorrectedSourceProfile(
        precision=precision,
        order=order,
        schedule=tail.schedule,
        angular=tail.correction,
        tail=tail,
    )
    return profile, Lambda, log_cstar


def run() -> dict[str, Any]:
    precision = 160
    order = 128
    profile, Lambda, log_cstar = _make_true_source_profile(
        precision=precision, order=order
    )
    jets = AxisPressureJets(profile, quadrature_order=order, R_a=4 / Lambda)
    bounds = PressureTailZBounds(jets, quadrature_order=order)
    samples = []
    for z in (0.0, 0.3, 0.7, 1.0):
        samples.append(bounds.receipt(z, quadrature_order=order))

    refinement = []
    for z in (0.3, 0.7):
        coarse = bounds.baseline_derivative_bound(z, quadrature_order=64)
        fine = bounds.baseline_derivative_bound(z, quadrature_order=128)
        with mp.workdps(precision):
            c = _from_signed_log(coarse["total_bound_normalized"])
            f = _from_signed_log(fine["total_bound_normalized"])
            difference = abs(f - c)
            relative = abs(difference / f) if f else mp.mpf("0")
        refinement.append(
            {
                "Z": z,
                "order_64": coarse["total_bound_normalized"],
                "order_128": fine["total_bound_normalized"],
                "absolute_difference": _signed_log(difference, precision),
                "relative_difference": mp.nstr(relative, precision),
            }
        )

    report = {
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_equations": "4.14--4.25, 7.17--7.23",
        "parameters": {
            "Lambda": mp.nstr(Lambda, precision),
            "logCstar": mp.nstr(log_cstar, precision),
            "logPstar": "14",
            "delta": "1e-32",
            "R_a": mp.nstr(4 / Lambda, precision),
            "quadrature_order": order,
        },
        "profile_construction": profile.metadata()["construction"],
        "schedule_conditions": profile.schedule.metadata().get(
            "theorem_conditions", {}
        ),
        "samples": samples,
        "quadrature_refinement": refinement,
        "checks": {
            "true_source_profile_Lambda_1e36": True,
            "flattening_prefactor_bound_used": True,
            "heat_factor_derivative_bound_used": True,
            "nominal_independent_interval_marked": True,
            "angular_bump_derivative_explicitly_unresolved": True,
            "finite_difference_used": False,
            "full_corrected_tail_C1_bound": False,
            "PDE_or_theorem_certified": False,
        },
        "scope": (
            "The baseline uncorrected future-pressure derivative has an "
            "analytic pointwise envelope evaluated by finite Gauss quadrature. "
            "The same-profile angular bump value is retained, but d1_Z and d2_Z "
            "are unresolved; no corrected-tail C1 bound is claimed."
        ),
    }
    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "baseline_physical_bounds": [
                    row["uncorrected_baseline_derivative_bound"][
                        "total_bound_physical"
                    ]
                    for row in samples
                ],
                "correction_derivative_status": [
                    row["corrected_bump_derivative"]["coefficient_derivative_status"]
                    for row in samples
                ],
                "quadrature_refinement": refinement,
            },
            default=str,
        )
    )
    return report


if __name__ == "__main__":
    run()
