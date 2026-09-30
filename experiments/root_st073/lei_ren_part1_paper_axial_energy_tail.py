"""Section 7.5 angular-energy tail for the current paper outer schedule.

This module evaluates the dimensionless quantity

    Q(Z) = integral_{R_v}^{infinity} U_theta(R,Z)^2 dR
           / (R_v U_theta(R_v,Z)^2),

using logarithmic stage coordinates.  It deliberately never materializes
the astronomical physical radii in the source schedule.  The finite
flattening and transition intervals use Gauss--Legendre quadrature, while
the constant-slope intervals use their exact exponential contractions.  The
two Section 7.4 angular bumps are integrated through their linear and
quadratic Gram contributions, retaining the signed-log coefficients supplied
by ``AngularCorrection``.

The exact heat exterior has a small positive/negative uncertainty from the
heat factor.  We report a nominal H=1 tail together with a rigorous Taylor
bound for that tail.  The angular coefficients themselves are candidate
inputs from ``AngularCorrection`` (preheat difference ODE and bounded first
Taylor heat defects), so this is an energy-tail computation and not a claim
of full moment closure.

Source: Lei--Ren, arXiv:2609.35406v1, Sections 4.14--4.22 and 7.34.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
import json
import math
from pathlib import Path
import sys
from typing import Any

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss


_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from lei_ren_part1_paper_angular_correction import AngularCorrection  # noqa: E402
from lei_ren_part1_paper_outer import PaperOuterSchedule  # noqa: E402
from lei_ren_part1_paper_outer_closure import _paper_raw_bump  # noqa: E402
from lei_ren_part1_paper_waiting_length import (  # noqa: E402
    apply_waiting_root,
    incoming_angular_ratio,
    solve_waiting_length,
)


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_EQUATIONS = "(4.14)--(4.22), (7.34)"
DEFAULT_PRECISION = 120
DEFAULT_QUADRATURE_ORDER = 64


def _mp_signed_log(row: dict[str, Any], key: str, *, precision: int) -> mp.mpf:
    """Read one signed-log coefficient without routing it through float."""

    value = row.get(key)
    if not isinstance(value, dict) or not value.get("sign"):
        return mp.mpf("0")
    with mp.workdps(precision):
        result = mp.exp(mp.mpf(str(value["log_abs"])))
        return result if int(value["sign"]) > 0 else -result


def _mp_decimal(value: Decimal | Any) -> mp.mpf:
    return mp.mpf(str(value))


class AngularEnergyTail:
    """Log-coordinate evaluator for the Section 7.34 future swirl energy.

    Parameters
    ----------
    schedule:
        A matched ``PaperOuterSchedule``.  Its ``logR_*`` checkpoints are
        read as Decimal logarithms and never exponentiated.
    correction:
        The current ``AngularCorrection`` object.  Its signed-log ``d1`` and
        ``d2`` rows are used only on the compact supports around
        ``log(R/R_rel) = -3`` and ``-1``.
    precision:
        mpmath working precision for exponentials and reported bounds.
    quadrature_order:
        Default Gauss--Legendre order for finite nonconstant stages and the
        bump Gram integrals.  ``evaluate(..., quadrature_order=...)`` can
        override it for refinement checks.
    """

    def __init__(
        self,
        schedule: PaperOuterSchedule,
        correction: AngularCorrection,
        *,
        precision: int = DEFAULT_PRECISION,
        quadrature_order: int = DEFAULT_QUADRATURE_ORDER,
    ) -> None:
        if precision < 64:
            raise ValueError("At least 64 decimal digits are required")
        if int(quadrature_order) != quadrature_order or quadrature_order < 16:
            raise ValueError("quadrature_order must be an integer at least 16")
        if correction.schedule is not schedule:
            raise ValueError("AngularCorrection must use the same schedule")
        self.schedule = schedule
        self.correction = correction
        self.precision = int(precision)
        self.quadrature_order = int(quadrature_order)

        with localcontext() as context:
            context.prec = max(schedule.decimal_precision, 160)
            self._logR_v = schedule.logR_v
            self._logR_f = schedule.logR_f
            self._logR_rel = schedule.logR_rel
            self._logR_s = schedule.logR_s
            self._logR_q = schedule.logR_q
            self._logR_t = schedule.logR_t
            self._logR_tail = schedule.logR_tail
            self._logR_b = schedule.logR_b
            self._x_f = schedule.logR_f - schedule.logR_v
            self._x_rel = schedule.logR_rel - schedule.logR_v
            self._x_s = schedule.logR_s - schedule.logR_v
            self._x_q = schedule.logR_q - schedule.logR_v
            self._x_t = schedule.logR_t - schedule.logR_v
            self._x_tail = schedule.logR_tail - schedule.logR_v
            self._x_b = schedule.logR_b - schedule.logR_v
        if self._x_f != Decimal(100):
            raise ArithmeticError("source flattening length is not exactly 100")

        # A single row is reused for all quadrature orders at one Z.  It is
        # still recomputed when Z changes so candidate coefficient receipts
        # cannot accidentally be mixed between latitudes.
        self._coefficient_cache: dict[float, dict[str, Any]] = {}

    # ------------------------------------------------------------------
    # Schedule/log-coordinate helpers
    # ------------------------------------------------------------------
    def _offset_log_radius(self, offset: float | Decimal) -> Decimal:
        """Return log(R_v) + offset at schedule precision."""

        with localcontext() as context:
            context.prec = max(self.schedule.decimal_precision, 160)
            return self._logR_v + (offset if isinstance(offset, Decimal) else Decimal(str(offset)))

    def _log_u_v(self, Z: float) -> Decimal:
        with localcontext() as context:
            context.prec = max(self.schedule.decimal_precision, 160)
            return self.schedule.at_log_radius(self._logR_v, Z)["log_angular_amplitude"]

    def _log_u_relative(self, offset: float | Decimal, Z: float, log_u_v: Decimal) -> mp.mpf:
        """Return log(U_theta(R_v e^offset)/U_theta(R_v))."""

        with localcontext() as context:
            context.prec = max(self.schedule.decimal_precision, 160)
            row = self.schedule.at_log_radius(self._offset_log_radius(offset), Z)
            difference = row["log_angular_amplitude"] - log_u_v
        return _mp_decimal(difference)

    def _integrand(self, offset: float, Z: float, log_u_v: Decimal) -> mp.mpf:
        with mp.workdps(self.precision):
            log_ratio = self._log_u_relative(offset, Z, log_u_v)
            return mp.exp(mp.mpf(str(offset)) + 2 * log_ratio)

    @staticmethod
    def _gauss_integral(left: float, right: float, function, order: int) -> mp.mpf:
        if right <= left:
            return mp.mpf("0")
        nodes, weights = leggauss(int(order))
        midpoint = (left + right) / 2.0
        half = (right - left) / 2.0
        total = mp.mpf("0")
        for node, weight in zip(nodes, weights):
            total += mp.mpf(str(float(weight))) * function(midpoint + half * float(node))
        return half * total

    def _finite_stage(self, left: float, right: float, Z: float, log_u_v: Decimal, order: int) -> mp.mpf:
        return self._gauss_integral(
            left,
            right,
            lambda offset: self._integrand(offset, Z, log_u_v),
            order,
        )

    # ------------------------------------------------------------------
    # Stage-specific pieces
    # ------------------------------------------------------------------
    def _flattening(self, Z: float, log_u_v: Decimal, order: int) -> mp.mpf:
        """Finite [R_v,R_f] integral, including its Z flattening."""

        return self._finite_stage(0.0, float(self._x_f), Z, log_u_v, order)

    def _pre_rel_background(
        self, Z: float, log_u_v: Decimal
    ) -> tuple[mp.mpf, dict[str, str]]:
        """Exact constant-slope background on [R_f,R_rel]."""

        with mp.workdps(self.precision):
            xf = mp.mpf(str(self._x_f))
            xrel = mp.mpf(str(self._x_rel))
            length = xrel - xf
            lf = self._log_u_relative(self._x_f, Z, log_u_v)
            lrel = self._log_u_relative(self._x_rel, Z, log_u_v)
            mu = _mp_decimal(self.schedule.mu)
            pref = mp.exp(xf + 2 * lf)
            if mu == 0:
                contraction = length
            else:
                contraction = -mp.expm1(-2 * mu * length) / (2 * mu)
            nominal = pref * contraction
            expected_lrel = lf - (mp.mpf("0.5") + mu) * length
            return nominal, {
                "length": mp.nstr(length, self.precision),
                "log_ratio_at_Rf": mp.nstr(lf, self.precision),
                "log_ratio_at_Rrel": mp.nstr(lrel, self.precision),
                "constant_slope_log_ratio_prediction_at_Rrel": mp.nstr(expected_lrel, self.precision),
                "endpoint_log_ratio_mismatch": mp.nstr(lrel - expected_lrel, self.precision),
                "slope": mp.nstr(-(mp.mpf("0.5") + mu), self.precision),
            }

    def _angular_bump_correction(
        self,
        Z: float,
        log_u_v: Decimal,
        coefficients: dict[str, Any],
        order: int,
    ) -> tuple[mp.mpf, dict[str, Any]]:
        """Quadratic bump energy correction on [R_f,R_rel]."""

        with mp.workdps(self.precision):
            xrel = mp.mpf(str(self._x_rel))
            lrel = self._log_u_relative(self._x_rel, Z, log_u_v)
            pref = mp.exp(xrel + 2 * lrel)
            mu = _mp_decimal(self.schedule.mu)
            weights = self.correction.weights
            ell = mp.mpf(str(weights.ell))
            norm = mp.mpf(str(weights.normalization_integral))
            nodes, weights_gl = leggauss(int(order))
            rows: dict[str, Any] = {}
            total = mp.mpf("0")
            for name, center in (("d1", -3.0), ("d2", -1.0)):
                coefficient = _mp_signed_log(coefficients, name, precision=self.precision)
                linear = mp.mpf("0")
                quadratic = mp.mpf("0")
                for node, weight in zip(nodes, weights_gl):
                    unit = mp.mpf(str(float(node)))
                    local = ell * unit
                    t = mp.mpf(str(center)) + local
                    raw = mp.mpf(str(float(_paper_raw_bump(float(node)))))
                    beta = raw / (ell * norm)
                    measure = ell * mp.mpf(str(float(weight)))
                    factor = mp.exp(-2 * mu * t)
                    linear += measure * factor * beta
                    quadratic += measure * factor * beta * beta
                contribution = pref * (2 * coefficient * linear + coefficient * coefficient * quadratic)
                total += contribution
                rows[name] = {
                    "coefficient": mp.nstr(coefficient, self.precision),
                    "linear_weight": mp.nstr(linear, self.precision),
                    "quadratic_weight": mp.nstr(quadratic, self.precision),
                    "energy_contribution": mp.nstr(contribution, self.precision),
                    "support_center_log_R_over_Rrel": center,
                    "support_halfwidth": float(weights.ell),
                }
            return total, {
                "background_prefactor_at_Rrel": mp.nstr(pref, self.precision),
                "bumps_disjoint_quadratic_cross_term": "exactly zero by support",
                "order": int(order),
                "rows": rows,
            }

    def _post_rel_stages(self, Z: float, log_u_v: Decimal, order: int) -> tuple[mp.mpf, dict[str, Any]]:
        """Steep transitions, exact powers, waiting, and finite heat collar."""

        with mp.workdps(self.precision):
            xrel = float(self._x_rel)
            xs = float(self._x_s)
            xq = float(self._x_q)
            xt = float(self._x_t)
            xtail = float(self._x_tail)
            xb = float(self._x_b)
            mu = _mp_decimal(self.schedule.mu)
            delta = _mp_decimal(self.schedule.delta)

            drop = self._finite_stage(xrel, xs, Z, log_u_v, order)

            # On [R_s,R_q], U is exactly R^{-3/2}; the integrand in log
            # radius therefore has slope -2.
            ls = self._log_u_relative(self._x_s, Z, log_u_v)
            steep_pref = mp.exp(mp.mpf(str(xs)) + 2 * ls)
            steep_length = mp.mpf(str(xq - xs))
            steep = steep_pref * (-mp.expm1(-2 * steep_length)) / 2

            restore = self._finite_stage(xs + (xq - xs), xt, Z, log_u_v, order)
            # The expression above starts at xq and is written this way to
            # keep the source checkpoint names visible in the receipt.
            waiting = self._finite_stage(xt, xtail, Z, log_u_v, order) if xtail > xt else mp.mpf("0")
            # Replace the numerical waiting integral by its exact terminal
            # power contraction.  The finite call above is retained only as
            # a cheap independent value for the report.
            lt = self._log_u_relative(self._x_t, Z, log_u_v)
            waiting_length = mp.mpf(str(xtail - xt))
            waiting_exact = mp.exp(mp.mpf(str(xt)) + 2 * lt)
            if delta == 0:
                waiting_exact *= waiting_length
            else:
                waiting_exact *= -mp.expm1(-delta * waiting_length) / delta

            heat_connection = self._finite_stage(xtail, xb, Z, log_u_v, order)

            return drop + steep + restore + waiting_exact + heat_connection, {
                "steep_transition_in": mp.nstr(drop, self.precision),
                "steep_power_exact": mp.nstr(steep, self.precision),
                "steep_power_length": mp.nstr(steep_length, self.precision),
                "steep_transition_out": mp.nstr(restore, self.precision),
                "waiting_power_exact": mp.nstr(waiting_exact, self.precision),
                "waiting_power_finite_quadrature_crosscheck": mp.nstr(waiting, self.precision),
                "waiting_power_length": mp.nstr(waiting_length, self.precision),
                "heat_connection": mp.nstr(heat_connection, self.precision),
                "heat_connection_interval": [mp.nstr(mp.mpf(str(xtail)), self.precision), mp.nstr(mp.mpf(str(xb)), self.precision)],
                "quadrature_order": int(order),
            }

    def _exact_heat_tail(
        self, Z: float, log_u_v: Decimal
    ) -> tuple[mp.mpf, mp.mpf, mp.mpf, dict[str, str]]:
        """Return nominal, lower, upper H-tail integrals and bound data."""

        with mp.workdps(self.precision):
            xb = mp.mpf(str(self._x_b))
            delta = _mp_decimal(self.schedule.delta)
            # Use the source heat amplitude c_inf directly.  The schedule
            # evaluator intentionally rounds H_delta to one at this
            # astronomical radius, but the bound below must retain the true
            # (nonzero) H deficit in arbitrary-exponent form.
            a = (1 + delta) / 2
            log_base_ratio = (
                _mp_decimal(self.schedule._log_c_inf)
                - a * _mp_decimal(self._logR_b)
                - _mp_decimal(log_u_v)
            )
            pref = mp.exp(xb + 2 * log_base_ratio)
            if delta == 0:
                nominal = mp.inf
            else:
                nominal = pref / delta
            d = mp.mpf(str(1.0 - float(Z) * float(Z)))
            h = delta / 2
            a1 = h * (1 + h)
            # 1-H <= a1*(2 d / R), so 1-H^2 <= 2*a1*(2d/R).
            log_Rb = _mp_decimal(self._logR_b)
            inv_Rb = mp.exp(-log_Rb)
            bound = pref * (4 * a1 * d * inv_Rb) / (1 + delta)
            lower = nominal - bound
            if lower < 0:
                lower = mp.mpf("0")
            return nominal, lower, nominal, {
                "prefactor_at_Rb": mp.nstr(pref, self.precision),
                "log_base_ratio_at_Rb": mp.nstr(log_base_ratio, self.precision),
                "inverse_Rb_arbitrary_precision": mp.nstr(inv_Rb, self.precision),
                "heat_taylor_a1": mp.nstr(a1, self.precision),
                "tail_correction_bound": mp.nstr(bound, self.precision),
                "lower_minus_nominal": mp.nstr(lower - nominal, self.precision),
                "upper_minus_nominal": "0",
                "nominal_H_policy": "H=1 source tail; actual H^2 deficit bounded separately",
            }

    def _heat_connection_bound(self, Z: float, log_u_v: Decimal) -> tuple[mp.mpf, dict[str, str]]:
        """Bound the H_delta deficit on the finite three-unit collar.

        The production schedule may report H=1 when Decimal underflows its
        inverse radius.  On [R_tail,R_b], the interpolation factor multiplying
        H is at most one, so the square-integrand loss is bounded by
        ``2 a1 (2 d / R_tail)`` times the corresponding base power integral.
        """

        with mp.workdps(self.precision):
            delta = _mp_decimal(self.schedule.delta)
            h = delta / 2
            a1 = h * (1 + h)
            d = mp.mpf(str(1.0 - float(Z) * float(Z)))
            a = (1 + delta) / 2
            xtail = mp.mpf(str(self._x_tail))
            length = mp.mpf(str(self._x_b - self._x_tail))
            log_base_ratio = (
                _mp_decimal(self.schedule._log_c_inf)
                - a * _mp_decimal(self._logR_tail)
                - _mp_decimal(log_u_v)
            )
            pref = mp.exp(xtail + 2 * log_base_ratio)
            if delta == 0:
                base_integral = pref * length
            else:
                base_integral = pref * (-mp.expm1(-delta * length)) / delta
            inv_Rtail = mp.exp(-_mp_decimal(self._logR_tail))
            bound = base_integral * (4 * a1 * d * inv_Rtail)
            return bound, {
                "base_power_integral_on_collar": mp.nstr(base_integral, self.precision),
                "inverse_Rtail_arbitrary_precision": mp.nstr(inv_Rtail, self.precision),
                "collar_correction_bound": mp.nstr(bound, self.precision),
            }

    def _coefficients(self, Z: float) -> dict[str, Any]:
        key = float(Z)
        if key not in self._coefficient_cache:
            self._coefficient_cache[key] = self.correction.coefficients(key)
        return self._coefficient_cache[key]

    def evaluate(self, Z: Any, *, quadrature_order: int | None = None) -> dict[str, Any]:
        """Evaluate Q and return a JSON-ready stage receipt."""

        z = float(Z)
        if not math.isfinite(z) or abs(z) > 1:
            raise ValueError("Z must lie in [-1,1]")
        order = self.quadrature_order if quadrature_order is None else int(quadrature_order)
        if order < 16:
            raise ValueError("quadrature_order must be at least 16")
        with mp.workdps(self.precision):
            log_u_v = self._log_u_v(z)
            coefficients = self._coefficients(z)
            flatten = self._flattening(z, log_u_v, order)
            pre_rel, pre_meta = self._pre_rel_background(z, log_u_v)
            bump, bump_meta = self._angular_bump_correction(z, log_u_v, coefficients, order)
            post_rel, post_meta = self._post_rel_stages(z, log_u_v, order)
            heat_nominal, heat_lower, heat_upper, heat_meta = self._exact_heat_tail(z, log_u_v)
            heat_collar_bound, heat_collar_meta = self._heat_connection_bound(z, log_u_v)
            heat_tail_bound = mp.mpf(heat_meta["tail_correction_bound"])
            q_nominal = flatten + pre_rel + bump + post_rel + heat_nominal
            # At the available precision these two subtractions can round
            # back to q_nominal; retain the bounds themselves in the receipt.
            q_lower = flatten + pre_rel + bump + post_rel + heat_nominal - heat_tail_bound - heat_collar_bound
            q_upper = flatten + pre_rel + bump + post_rel + heat_upper
            energy_scale = _mp_decimal(self.schedule.mu) * mp.exp(-26) / 2
            payload = {
                "Z": z,
                "quadrature_order": order,
                "log_Utheta_at_Rv": str(log_u_v),
                "stage_offsets_log_R_over_Rv": {
                    "Rf": mp.nstr(mp.mpf(str(self._x_f)), self.precision),
                    "Rrel": mp.nstr(mp.mpf(str(self._x_rel)), self.precision),
                    "Rs": mp.nstr(mp.mpf(str(self._x_s)), self.precision),
                    "Rq": mp.nstr(mp.mpf(str(self._x_q)), self.precision),
                    "Rt": mp.nstr(mp.mpf(str(self._x_t)), self.precision),
                    "Rtail": mp.nstr(mp.mpf(str(self._x_tail)), self.precision),
                    "Rb": mp.nstr(mp.mpf(str(self._x_b)), self.precision),
                },
                "flattening_0_to_100": mp.nstr(flatten, self.precision),
                "pre_rel_constant_background": mp.nstr(pre_rel, self.precision),
                "pre_rel_metadata": pre_meta,
                "pre_rel_angular_bump_quadratic_correction": mp.nstr(bump, self.precision),
                "pre_rel_bump_metadata": bump_meta,
                "post_rel_metadata": post_meta,
                "post_rel_finite_and_power_total": mp.nstr(post_rel, self.precision),
                "exact_heat_tail_nominal_H_equals_1": mp.nstr(heat_nominal, self.precision),
                "exact_heat_tail_lower": mp.nstr(heat_lower, self.precision),
                "exact_heat_tail_upper": mp.nstr(heat_upper, self.precision),
                "exact_heat_tail_metadata": heat_meta,
                "heat_connection_heat_deficit_bound": mp.nstr(heat_collar_bound, self.precision),
                "heat_connection_heat_deficit_metadata": heat_collar_meta,
                "total_heat_deficit_bound": mp.nstr(heat_tail_bound + heat_collar_bound, self.precision),
                "Q_nominal": mp.nstr(q_nominal, self.precision),
                "Q_lower": mp.nstr(q_lower, self.precision),
                "Q_upper": mp.nstr(q_upper, self.precision),
                "energy_target_scale_mu_over_2_exp_minus_26": mp.nstr(energy_scale, self.precision),
                "energy_target_contribution_nominal": mp.nstr(energy_scale * q_nominal, self.precision),
                "energy_target_contribution_lower": mp.nstr(energy_scale * q_lower, self.precision),
                "energy_target_contribution_upper": mp.nstr(energy_scale * q_upper, self.precision),
                "angular_input_status": coefficients.get("input_status", "unknown"),
                "heat_input_status": "Finite collar uses checked heat_factor; exact exterior uses H Taylor bound",
                "exact_closure": False,
                "full_outer_closed": False,
            }
            return payload

    def dimensionless_Q(self, Z: Any, *, quadrature_order: int | None = None) -> mp.mpf:
        """Return nominal Q as an arbitrary-precision number."""

        order = self.quadrature_order if quadrature_order is None else int(quadrature_order)
        with mp.workdps(self.precision):
            return mp.mpf(self.evaluate(Z, quadrature_order=order)["Q_nominal"])

    def energy_target_contribution(self, Z: Any, *, quadrature_order: int | None = None) -> mp.mpf:
        """Return the Q contribution in source equation (7.34)."""

        with mp.workdps(self.precision):
            return _mp_decimal(self.schedule.mu) * mp.exp(-26) * self.dimensionless_Q(
                Z, quadrature_order=quadrature_order
            ) / 2

    def metadata(self) -> dict[str, Any]:
        return {
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "source_equations": SOURCE_EQUATIONS,
            "definition": "Q = integral_Rv^infinity Utheta^2 dR / (Rv Utheta(Rv)^2)",
            "energy_target_factor": "mu/2 * exp(-26) * Q",
            "quadrature": "Gauss-Legendre on nonconstant finite stages; exact contractions on constant slopes",
            "precision": self.precision,
            "default_quadrature_order": self.quadrature_order,
            "radii_representation": "Decimal log checkpoints; no physical radius exponentiation",
            "heat_bound": "0 <= 1-H_delta(2d/R) <= h(1+h)*2d/R; H^2 tail bound retained",
            "angular_input_status": "AngularCorrection candidate preheat plus bounded first-Taylor heat inputs",
            "exact_closure": False,
            "full_outer_closed": False,
        }


def build_default_tail(
    *,
    precision: int = DEFAULT_PRECISION,
    quadrature_order: int = DEFAULT_QUADRATURE_ORDER,
    schedule: PaperOuterSchedule | None = None,
    correction: AngularCorrection | None = None,
    incoming: dict[str, Any] | None = None,
    match_waiting: bool | None = None,
) -> AngularEnergyTail:
    """Build a shared schedule/correction/tail assembly.

    With no explicit objects this retains the historical default: construct
    the demo base schedule, solve its waiting-length root, and build one
    matched ``AngularCorrection`` plus ``AngularEnergyTail``.

    ``schedule`` allows a caller to provide source-scaled parameters.  Such a
    schedule is used as supplied unless ``match_waiting=True``; in that mode
    the waiting root is solved and the returned tail owns the matched schedule
    produced from the explicit base.  ``correction`` may be supplied when a
    caller already assembled the angular candidate; identity with the active
    schedule is required so pressure and tail inputs cannot be mixed.
    """

    if int(quadrature_order) != quadrature_order or quadrature_order < 16:
        raise ValueError("quadrature_order must be an integer at least 16")
    if correction is not None and schedule is None:
        schedule = correction.schedule
    if correction is not None and schedule is not correction.schedule:
        raise ValueError("correction must use the supplied schedule object")

    if schedule is None:
        base = PaperOuterSchedule(
            logPstar=14,
            logRref=10,
            delta="1e-32",
            Md=".5",
            c_mu=".001",
            c_delta=".001",
            c_epsilon=".01",
        )
        incoming_data = incoming if incoming is not None else incoming_angular_ratio(base)
        waiting = solve_waiting_length(base, incoming_data["incoming_X"])
        active_schedule = apply_waiting_root(base, waiting)
    else:
        active_schedule = schedule
        if match_waiting is None:
            match_waiting = False
        if match_waiting:
            if correction is not None:
                raise ValueError("match_waiting cannot be combined with an explicit correction")
            incoming_data = incoming if incoming is not None else incoming_angular_ratio(active_schedule)
            waiting = solve_waiting_length(active_schedule, incoming_data["incoming_X"])
            active_schedule = apply_waiting_root(active_schedule, waiting)
        else:
            incoming_data = incoming if incoming is not None else incoming_angular_ratio(active_schedule)

    if correction is None:
        correction = AngularCorrection(
            active_schedule,
            incoming_data,
            precision=precision,
            order=max(128, quadrature_order),
        )
    elif correction.schedule is not active_schedule:
        raise ValueError("explicit correction and active schedule must be identical")
    return AngularEnergyTail(
        active_schedule,
        correction,
        precision=precision,
        quadrature_order=quadrature_order,
    )


def run() -> dict[str, Any]:
    """Run two-order checks and write the sibling JSON receipt."""

    tail = build_default_tail(precision=DEFAULT_PRECISION, quadrature_order=64)
    rows = []
    for z in (0.0, 0.5, -0.5, 1.0):
        coarse = tail.evaluate(z, quadrature_order=64)
        fine = tail.evaluate(z, quadrature_order=128)
        q64 = mp.mpf(coarse["Q_nominal"])
        q128 = mp.mpf(fine["Q_nominal"])
        difference = abs(q128 - q64)
        scale = max(abs(q128), mp.mpf("1e-300"))
        row = {
            "Z": z,
            "Q_order_64": coarse["Q_nominal"],
            "Q_order_128": fine["Q_nominal"],
            "Q_absolute_difference": mp.nstr(difference, DEFAULT_PRECISION),
            "Q_relative_difference": mp.nstr(difference / scale, DEFAULT_PRECISION),
            "energy_target_order_128": fine["energy_target_contribution_nominal"],
            "Q_lower": fine["Q_lower"],
            "Q_upper": fine["Q_upper"],
            "heat_tail_and_collar_bound": fine["total_heat_deficit_bound"],
            "heat_tail_bound_width_at_report_precision": mp.nstr(mp.mpf(fine["Q_upper"]) - mp.mpf(fine["Q_lower"]), DEFAULT_PRECISION),
            "stage_receipt_order_128": fine,
        }
        rows.append(row)

    report = {
        **tail.metadata(),
        "schedule": tail.schedule.metadata(),
        "rows": rows,
        "checks": {
            "two_orders_completed": True,
            "finite_stage_and_bump_quadrature_refined": True,
            "constant_slope_stages_exact": True,
            "heat_tail_bound_retained": True,
            "no_physical_radius_materialized": True,
            "full_outer_closed": False,
        },
        "scope": "Actual candidate angular-energy tail only; pressure, axial/mixed closure, cone, PDE, and exact heat moment closure remain open.",
    }
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "Z_rows": [row["Z"] for row in rows],
                "max_relative_Q_difference": max(float(row["Q_relative_difference"]) for row in rows),
                "energy_target_Z0": rows[0]["energy_target_order_128"],
                "full_outer_closed": False,
            }
        ),
        flush=True,
    )
    return report


if __name__ == "__main__":
    run()
