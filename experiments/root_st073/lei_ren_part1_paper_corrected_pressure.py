"""Same-swirl pressure adapter for the Section 7.4 angular correction.

The source angular correction changes the swirl on the two compact supports
around ``t = log(R/R_rel) = -3`` and ``-1`` by

    U_theta = E_rel exp(-(1/2+mu)t) (1 + h),
    h = d1 beta(t+3) + d2 beta(t+1).

Consequently its pressure contribution, divided by ``E_rel**2``, is

    Delta_p_bump = integral exp(-(1+2 mu)t) (h + h**2/2) dt.

This file keeps that correction separate from the existing uncorrected
``normalized_outer_pressure`` baseline.  The angular coefficients currently
come from the candidate preheat difference and bounded first-Taylor heat
inputs, so cancellation against the extremely small source pressure target
may be below the coefficient receipt precision.  Such cancellation is
reported rather than silently rounded away.

Source: Lei--Ren, arXiv:2609.35406v1, (4.25), (7.17), and (7.20)--(7.23).
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
from lei_ren_part1_paper_outer_pressure import normalized_outer_pressure  # noqa: E402
from lei_ren_part1_paper_waiting_length import (  # noqa: E402
    apply_waiting_root,
    incoming_angular_ratio,
    solve_waiting_length,
)


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_EQUATIONS = "(4.25), (7.17), (7.20)--(7.23)"
DEFAULT_PRECISION = 160
DEFAULT_QUADRATURE_ORDER = 64


def _signed_log(value: mp.mpf, precision: int) -> dict[str, Any]:
    with mp.workdps(precision):
        if value == 0:
            return {"sign": 0, "log_abs": None, "arbitrary_exponent_value": "0"}
        return {
            "sign": 1 if value > 0 else -1,
            "log_abs": mp.nstr(mp.log(abs(value)), precision),
            "arbitrary_exponent_value": mp.nstr(value, precision),
        }


def _coefficient_value(row: dict[str, Any], key: str, precision: int) -> mp.mpf:
    entry = row.get(key)
    if not isinstance(entry, dict) or not entry.get("sign"):
        return mp.mpf("0")
    with mp.workdps(precision):
        result = mp.exp(mp.mpf(str(entry["log_abs"])))
        return result if int(entry["sign"]) > 0 else -result


def _decimal_mp(value: Decimal | Any) -> mp.mpf:
    return mp.mpf(str(value))


class CorrectedPressureAdapter:
    """Pressure correction induced by the actual signed-log angular bumps."""

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
        self._coefficients_cache: dict[float, dict[str, Any]] = {}

        with localcontext() as context:
            context.prec = max(schedule.decimal_precision, 160)
            self._logR_rel = schedule.logR_rel
            self._logR_ref = schedule.logRref
            self._logR_b = schedule.logR_b
            self._logR_tail = schedule.logR_tail
            self._logU_rel = schedule.at_log_radius(
                schedule.logR_rel, 0.0
            )["log_angular_amplitude"]
            self._logPstar = schedule.logPstar
            self._logR_ref_plus_one = schedule.logRref + Decimal(1)
        self.supports = tuple(
            (center - float(correction.weights.ell), center + float(correction.weights.ell))
            for center in (-3.0, -1.0)
        )
        with localcontext() as context:
            context.prec = max(schedule.decimal_precision, 160)
            self._stage_specs = (
                ("heat_connection", schedule.logR_tail, schedule.logR_b, "variable", None),
                ("waiting", schedule.logR_t, schedule.logR_tail, "constant", -(Decimal(1) + schedule.delta) / Decimal(2)),
                ("steep_transition_out", schedule.logR_q, schedule.logR_t, "variable", None),
                ("steep_power", schedule.logR_s, schedule.logR_q, "constant", Decimal("-1.5")),
                ("steep_transition_in", schedule.logR_rel, schedule.logR_s, "variable", None),
                ("power_buffer_rel", schedule.logR_f, schedule.logR_rel, "constant", -(Decimal("0.5") + schedule.mu)),
                ("z_flatten", schedule.logR_v, schedule.logR_f, "variable", None),
                ("pulse_reserved", schedule.logR_p, schedule.logR_v, "constant", -(Decimal("0.5") + schedule.mu)),
                ("power_buffer", schedule.logR_w, schedule.logR_p, "constant", -(Decimal("0.5") + schedule.mu)),
                ("slope_transition_mu", schedule.logR_d, schedule.logR_w, "variable", None),
                ("axial_turnoff", schedule.logRref + Decimal(1), schedule.logR_d, "variable", None),
                ("slope_transition_ref", schedule.logRref, schedule.logRref + Decimal(1), "variable", None),
            )

    def _coefficients(self, Z: float) -> dict[str, Any]:
        key = float(Z)
        if key not in self._coefficients_cache:
            self._coefficients_cache[key] = self.correction.coefficients(key)
        return self._coefficients_cache[key]

    def _log_Erel_over_Pstar(self, Z: float) -> mp.mpf:
        """Return log(E_rel/P_star), preserving Decimal differences."""

        with localcontext() as context:
            context.prec = max(self.schedule.decimal_precision, 160)
            row = self.schedule.at_log_radius(self.schedule.logR_rel, Z)
            difference = row["log_angular_amplitude"] - self.schedule.logPstar
        return _decimal_mp(difference)

    def _bump_integrals(
        self,
        Z: float,
        coefficients: dict[str, Any],
        *,
        order: int,
        lower_t: float | None = None,
        upper_t: float | None = None,
    ) -> tuple[mp.mpf, dict[str, Any]]:
        """Integrate the truncated future pressure increment.

        ``lower_t`` and ``upper_t`` are in ``log(R/R_rel)``.  A missing bound
        means the full source support.  The two supports are disjoint, so no
        cross term is present.
        """

        if lower_t is None:
            lower_t = -math.inf
        if upper_t is None:
            upper_t = 0.0
        if upper_t <= lower_t:
            return mp.mpf("0"), {"rows": [], "truncated": True, "order": int(order)}

        with mp.workdps(self.precision):
            mu = _decimal_mp(self.schedule.mu)
            lam = 1 + 2 * mu
            ell = mp.mpf(str(self.correction.weights.ell))
            normalization = mp.mpf(str(self.correction.weights.normalization_integral))
            nodes, weights = leggauss(int(order))
            rows = []
            total = mp.mpf("0")
            for key, center in (("d1", -3.0), ("d2", -1.0)):
                support_left = center - float(ell)
                support_right = center + float(ell)
                left = max(float(lower_t), support_left)
                right = min(float(upper_t), support_right)
                coefficient = _coefficient_value(coefficients, key, self.precision)
                if right <= left or coefficient == 0:
                    rows.append({
                        "name": key,
                        "support": [support_left, support_right],
                        "active_interval": None,
                        "coefficient": _signed_log(coefficient, self.precision),
                        "linear_integral": "0",
                        "quadratic_integral": "0",
                        "pressure_term": _signed_log(mp.mpf("0"), self.precision),
                    })
                    continue

                midpoint = (left + right) / 2.0
                half = (right - left) / 2.0
                linear_weight = mp.mpf("0")
                quadratic_weight = mp.mpf("0")
                for node, weight in zip(nodes, weights):
                    local = midpoint + half * float(node)
                    unit = (local - center) / float(ell)
                    raw = mp.mpf(str(float(_paper_raw_bump(unit))))
                    beta = raw / (ell * normalization)
                    measure = half * mp.mpf(str(float(weight)))
                    radial = mp.exp(-lam * mp.mpf(str(local)))
                    linear_weight += measure * radial * beta
                    quadratic_weight += measure * radial * beta * beta
                linear_term = coefficient * linear_weight
                quadratic_term = mp.mpf("0.5") * coefficient * coefficient * quadratic_weight
                pressure_term = linear_term + quadratic_term
                total += pressure_term
                rows.append({
                    "name": key,
                    "support": [support_left, support_right],
                    "active_interval": [left, right],
                    "coefficient": _signed_log(coefficient, self.precision),
                    "linear_integral": mp.nstr(linear_weight, self.precision),
                    "quadratic_integral": mp.nstr(quadratic_weight, self.precision),
                    "linear_pressure_term": _signed_log(linear_term, self.precision),
                    "quadratic_pressure_term": _signed_log(quadratic_term, self.precision),
                    "pressure_term": _signed_log(pressure_term, self.precision),
                })
            return total, {
                "rows": rows,
                "order": int(order),
                "lower_t": lower_t,
                "upper_t": upper_t,
                "supports_disjoint_cross_term": "exactly zero",
            }

    def pressure_increment(self, Z: Any, *, quadrature_order: int | None = None) -> dict[str, Any]:
        """Return Delta_p_bump normalized by E_rel squared."""

        z = float(Z)
        if not math.isfinite(z) or abs(z) > 1:
            raise ValueError("Z must lie in [-1,1]")
        order = self.quadrature_order if quadrature_order is None else int(quadrature_order)
        if order < 16:
            raise ValueError("quadrature_order must be at least 16")
        coefficients = self._coefficients(z)
        with mp.workdps(self.precision):
            increment, receipt = self._bump_integrals(z, coefficients, order=order)
            heat = coefficients.get("heat", {})
            source_log = heat.get("log_s_H")
            source_target = mp.mpf("0") if source_log is None else mp.exp(mp.mpf(str(source_log)))
            source_signed = _signed_log(source_target, self.precision)
            return {
                "Z": z,
                "quadrature_order": order,
                "delta_p_bump_over_Erel_squared": _signed_log(increment, self.precision),
                "delta_p_bump_signed_log_receipt": _signed_log(increment, self.precision),
                "source_s_H_target": source_signed,
                "source_s_H_log": source_log,
                "pressure_term_receipt": receipt,
                "source_pressure_target_unresolved": bool(
                    increment != 0 if source_target == 0 else
                    increment <= 0 or abs((increment-source_target)/source_target)
                    > mp.power(10, -self.precision/2)
                ),
                "input_status": coefficients.get("input_status", "unknown"),
                "exact_pressure_target": False,
                "full_outer_closed": False,
            }

    def pressure_correction_at_stage(
        self,
        Z: Any,
        t: Any,
        *,
        quadrature_order: int | None = None,
    ) -> dict[str, Any]:
        """Pressure correction from a stage point to the future bump supports.

        The returned correction is ``(P_corrected-P_baseline)/E_rel**2``.
        Thus it is the negative of the future positive integrand whenever the
        swirl correction raises the pressure moment.
        """

        z = float(Z)
        stage_t = float(t)
        if not math.isfinite(z) or abs(z) > 1:
            raise ValueError("Z must lie in [-1,1]")
        if not math.isfinite(stage_t):
            raise ValueError("t must be finite")
        receipt = self.pressure_increment(z, quadrature_order=quadrature_order)
        order = self.quadrature_order if quadrature_order is None else int(quadrature_order)
        coefficients = self._coefficients(z)
        with mp.workdps(self.precision):
            future, truncated = self._bump_integrals(
                z,
                coefficients,
                order=order,
                lower_t=stage_t,
                upper_t=0.0,
            )
            correction = -future
            below = stage_t <= self.supports[0][0]
            above = stage_t >= self.supports[-1][1]
            if below:
                relation = "full correction below both supports"
            elif above:
                relation = "zero correction above both supports"
            else:
                relation = "future truncated support integral"
            return {
                "Z": z,
                "t_log_R_over_Rrel": stage_t,
                "pressure_correction_over_Erel_squared": _signed_log(correction, self.precision),
                "future_integral_over_Erel_squared": _signed_log(future, self.precision),
                "future_integral_receipt": truncated,
                "outside_support_behavior": relation,
                "full_below_support": below,
                "zero_above_support": above,
                "source_full_increment_at_Rrel": receipt["delta_p_bump_over_Erel_squared"],
                "exact_pressure_target": False,
            }

    # ------------------------------------------------------------------
    # Full baseline pressure propagation in logarithmic radius
    # ------------------------------------------------------------------
    @staticmethod
    def _coerce_log_radius(value: Any) -> Decimal:
        if isinstance(value, bool):
            raise TypeError("log_radius must be a finite real number")
        try:
            result = value if isinstance(value, Decimal) else Decimal(str(value))
        except Exception as exc:  # Decimal raises several concrete subclasses
            raise ValueError("log_radius must be a finite real number") from exc
        if not result.is_finite():
            raise ValueError("log_radius must be finite")
        return result

    def _log_u_at(self, log_radius: Decimal, Z: float) -> Decimal:
        """Authoritative nominal log swirl, with H=1 beyond R_b."""

        with localcontext() as context:
            context.prec = max(self.schedule.decimal_precision, 160)
            if log_radius >= self._logR_b:
                a = (Decimal(1) + self.schedule.delta) / Decimal(2)
                return self.schedule._log_c_inf - a * log_radius
            return self.schedule.at_log_radius(log_radius, Z)["log_angular_amplitude"]

    def log_radius_from_stage_offset(self, t: Any) -> Decimal:
        """Construct ``log(R_rel) + t`` without losing tiny local offsets."""

        with localcontext() as context:
            context.prec = max(self.schedule.decimal_precision, 160)
            return self._logR_rel + Decimal(str(t))

    def _stage_at_log_radius(self, log_radius: Decimal) -> str:
        with localcontext() as context:
            context.prec = max(self.schedule.decimal_precision, 160)
            if log_radius >= self._logR_b:
                return "exact_heat"
            if log_radius < self._logR_ref:
                return "reference"
            for name, left, right, _kind, _slope in self._stage_specs:
                if left <= log_radius <= right:
                    return name
            # The source intervals partition [R_ref,R_b].  A gap here means
            # a Decimal checkpoint was accidentally rounded before use.
            raise ArithmeticError("source pressure stage partition failed")

    def _heat_interval_bound(
        self,
        left: Decimal,
        right: Decimal,
        log_u_left: Decimal,
        Z: float,
        *,
        order: int,
    ) -> mp.mpf:
        """Bound pressure loss from replacing H_delta by one on a collar."""

        d = mp.mpf(str(1.0 - float(Z) * float(Z)))
        if d == 0:
            return mp.mpf("0")
        with mp.workdps(self.precision):
            delta = _decimal_mp(self.schedule.delta)
            h = delta / 2
            a1 = h * (1 + h)
            a_decimal = (Decimal(1) + self.schedule.delta) / Decimal(2)
            length = _decimal_mp(right - left)
            if length <= 0:
                return mp.mpf("0")
            nodes, weights = leggauss(int(order))
            half = float(length) / 2.0
            total = mp.mpf("0")
            for node, weight in zip(nodes, weights):
                local_offset = half * (1.0 + float(node))
                with localcontext() as context:
                    context.prec = max(self.schedule.decimal_precision, 160)
                    log_r = left + Decimal(str(local_offset))
                    base_log_u = self.schedule._log_c_inf - a_decimal * log_r
                    log_u = self.schedule.at_log_radius(log_r, Z)["log_angular_amplitude"]
                    ratio_log = base_log_u - log_u
                # Half of 1-H^2 is bounded by 2*a1*d/R.
                total += (
                    half
                    * mp.mpf(str(float(weight)))
                    * (2 * a1 * d)
                    * mp.exp(2 * _decimal_mp(ratio_log) - _decimal_mp(log_r))
                )
            return total

    def _propagate_constant(
        self,
        left: Decimal,
        right: Decimal,
        q_right: mp.mpf,
        bound_right: mp.mpf,
        slope: Decimal,
    ) -> tuple[mp.mpf, mp.mpf, dict[str, Any]]:
        with mp.workdps(self.precision):
            length = _decimal_mp(right - left)
            rate = _decimal_mp(slope)
            factor = mp.exp(2 * rate * length)
            equilibrium = 1 / (4 * rate)
            q_left = equilibrium + (q_right - equilibrium) * factor
            bound_left = bound_right * factor
            return q_left, bound_left, {
                "kind": "constant_slope_exact",
                "slope": mp.nstr(rate, self.precision),
                "length": mp.nstr(length, self.precision),
                "contraction": mp.nstr(factor, self.precision),
            }

    def _propagate_variable(
        self,
        left: Decimal,
        right: Decimal,
        q_right: mp.mpf,
        bound_right: mp.mpf,
        Z: float,
        *,
        order: int,
        stage_name: str,
    ) -> tuple[mp.mpf, mp.mpf, dict[str, Any]]:
        """Use q_left = (U_right/U_left)^2 q_right - 1/2 integral."""

        with mp.workdps(self.precision):
            length_decimal = right - left
            length = _decimal_mp(length_decimal)
            if length <= 0:
                return q_right, bound_right, {
                    "kind": "variable_zero_length",
                    "stage": stage_name,
                    "quadrature_order": int(order),
                }
            if not math.isfinite(float(length)) or float(length) > 1.0e5:
                raise ArithmeticError(f"variable stage {stage_name} is too long for local quadrature")
            with localcontext() as context:
                context.prec = max(self.schedule.decimal_precision, 160)
                log_u_left = self._log_u_at(left, Z)
                log_u_right = self._log_u_at(right, Z)
            endpoint_factor = mp.exp(2 * (_decimal_mp(log_u_right) - _decimal_mp(log_u_left)))
            nodes, weights = leggauss(int(order))
            half = float(length) / 2.0
            integral = mp.mpf("0")
            for node, weight in zip(nodes, weights):
                local_offset = half * (1.0 + float(node))
                with localcontext() as context:
                    context.prec = max(self.schedule.decimal_precision, 160)
                    # Add only the finite stage-local offset to the huge
                    # Decimal checkpoint; never convert the absolute log R
                    # to binary64.
                    log_r = left + Decimal(str(local_offset))
                    log_u = self._log_u_at(log_r, Z)
                    relative_log = log_u - log_u_left
                integral += (
                    half
                    * mp.mpf(str(float(weight)))
                    * mp.exp(2 * _decimal_mp(relative_log))
                )
            q_left = endpoint_factor * q_right - integral / 2
            bound_left = endpoint_factor * bound_right
            heat_extra = mp.mpf("0")
            if stage_name == "heat_connection":
                heat_extra = self._heat_interval_bound(left, right, log_u_left, Z, order=order)
                bound_left += heat_extra
            return q_left, bound_left, {
                "kind": "variable_stage_quadrature",
                "stage": stage_name,
                "quadrature_order": int(order),
                "length": mp.nstr(length, self.precision),
                "endpoint_factor": mp.nstr(endpoint_factor, self.precision),
                "integral": mp.nstr(integral, self.precision),
                "heat_deficit_bound_added": mp.nstr(heat_extra, self.precision),
            }

    def _propagate_interval(
        self,
        name: str,
        left: Decimal,
        right: Decimal,
        q_right: mp.mpf,
        bound_right: mp.mpf,
        Z: float,
        *,
        order: int,
    ) -> tuple[mp.mpf, mp.mpf, dict[str, Any]]:
        spec = next(item for item in self._stage_specs if item[0] == name)
        _name, _left, _right, kind, slope = spec
        if kind == "constant":
            return self._propagate_constant(left, right, q_right, bound_right, slope)
        return self._propagate_variable(
            left,
            right,
            q_right,
            bound_right,
            Z,
            order=order,
            stage_name=name,
        )

    def _nominal_heat_endpoint(self, Z: float) -> tuple[mp.mpf, mp.mpf, Decimal]:
        """q(R_b) and its heat-factor uncertainty, with H replaced by one."""

        with mp.workdps(self.precision):
            delta = _decimal_mp(self.schedule.delta)
            q_nominal = -1 / (2 * (1 + delta))
            d = mp.mpf(str(1.0 - float(Z) * float(Z)))
            h = delta / 2
            a1 = h * (1 + h)
            inv_rb = mp.exp(-_decimal_mp(self._logR_b))
            q_bound = (2 * a1 * d * inv_rb) / (2 + delta)
            return q_nominal, q_bound, self._log_u_at(self._logR_b, Z)

    def _baseline_state_to_log_radius(
        self,
        log_radius: Decimal,
        Z: float,
        *,
        order: int,
    ) -> tuple[mp.mpf, mp.mpf, Decimal, str, list[dict[str, Any]]]:
        """Backward-propagate q=P/U^2 from R_b to a requested log radius."""

        q, bound, _log_u_b = self._nominal_heat_endpoint(Z)
        trace: list[dict[str, Any]] = []
        with localcontext() as context:
            context.prec = max(self.schedule.decimal_precision, 160)
            ref = self._logR_ref
        if log_radius >= self._logR_b:
            return q, bound, self._log_u_at(log_radius, Z), "exact_heat", trace

        target = log_radius
        for name, left, right, _kind, _slope in self._stage_specs:
            if target >= left:
                # Target lies in this interval (including the left endpoint).
                q, bound, receipt = self._propagate_interval(
                    name,
                    target,
                    right,
                    q,
                    bound,
                    Z,
                    order=order,
                )
                receipt["target_interval"] = True
                trace.append({"stage": name, **receipt})
                log_u = self._log_u_at(target, Z)
                return q, bound, log_u, name, trace
            q, bound, receipt = self._propagate_interval(
                name,
                left,
                right,
                q,
                bound,
                Z,
                order=order,
            )
            receipt["target_interval"] = False
            trace.append({"stage": name, **receipt})

        # All source outer intervals have been traversed, so the target is in
        # the reference branch.  Its R^{1/10} power law is integrated exactly.
        if target < ref:
            with mp.workdps(self.precision):
                length = _decimal_mp(ref - target)
                slope = mp.mpf("0.1")
                factor = mp.exp(2 * slope * length)
                equilibrium = 1 / (4 * slope)
                q = equilibrium + (q - equilibrium) * factor
                bound *= factor
            trace.append({
                "stage": "reference",
                "kind": "constant_slope_exact",
                "slope": "0.1",
                "length": mp.nstr(length, self.precision),
                "contraction": mp.nstr(factor, self.precision),
                "target_interval": True,
            })
            return q, bound, self._log_u_at(target, Z), "reference", trace
        raise ArithmeticError("baseline propagation did not reach requested radius")

    def baseline_pressure_at_log_radius(
        self,
        log_radius: Any,
        Z: Any,
        *,
        quadrature_order: int | None = None,
    ) -> dict[str, Any]:
        """Return the uncorrected pressure at any logarithmic radius.

        Pressure is represented by signed arbitrary-exponent receipts.  The
        nominal profile uses H=1 at the exact heat tail; ``heat_deficit`` is a
        separately propagated upper bound and is never folded into the value.
        """

        z = float(Z)
        if not math.isfinite(z) or abs(z) > 1:
            raise ValueError("Z must lie in [-1,1]")
        target = self._coerce_log_radius(log_radius)
        order = self.quadrature_order if quadrature_order is None else int(quadrature_order)
        if order < 16:
            raise ValueError("quadrature_order must be at least 16")
        with mp.workdps(self.precision):
            q, q_bound, log_u, stage, trace = self._baseline_state_to_log_radius(
                target,
                z,
                order=order,
            )
            u_squared = mp.exp(2 * _decimal_mp(log_u))
            pressure = q * u_squared
            heat_pressure_bound = abs(q_bound * u_squared)
            p_over_pstar = pressure / mp.exp(2 * _decimal_mp(self.schedule.logPstar))
            bound_over_pstar = heat_pressure_bound / mp.exp(2 * _decimal_mp(self.schedule.logPstar))
            return {
                "Z": z,
                "log_radius": str(target),
                "stage": stage,
                "quadrature_order": order,
                "log_Utheta_nominal": str(log_u),
                "P_over_Utheta_squared": _signed_log(q, self.precision),
                "heat_deficit_bound_over_Utheta_squared": _signed_log(q_bound, self.precision),
                "P_nominal": _signed_log(pressure, self.precision),
                "heat_deficit_bound_in_P": _signed_log(heat_pressure_bound, self.precision),
                "P_over_Pstar_squared": _signed_log(p_over_pstar, self.precision),
                "heat_deficit_bound_over_Pstar_squared": _signed_log(bound_over_pstar, self.precision),
                "propagation_trace": trace,
                "nominal_heat_policy": "H=1 at R>=Rb",
                "full_heat_pressure_exact": False,
                "arbitrary_log_radius": True,
            }

    def _future_bump_from_decimal_t(
        self,
        Z: float,
        t: Decimal,
        *,
        order: int,
    ) -> tuple[mp.mpf, dict[str, Any], str]:
        leftmost = Decimal(str(self.supports[0][0]))
        rightmost = Decimal(str(self.supports[-1][1]))
        if t <= leftmost:
            value = self.pressure_increment(Z, quadrature_order=order)
            signed = value["delta_p_bump_over_Erel_squared"]
            return mp.mpf(str(signed["arbitrary_exponent_value"])), value["pressure_term_receipt"], "full_below_support"
        if t >= rightmost:
            return mp.mpf("0"), {"rows": [], "order": order, "truncated": True}, "zero_above_support"
        coefficients = self._coefficients(Z)
        value, receipt = self._bump_integrals(
            Z,
            coefficients,
            order=order,
            lower_t=float(t),
            upper_t=0.0,
        )
        return value, receipt, "future_truncated_support_integral"

    def corrected_pressure_at_log_radius(
        self,
        log_radius: Any,
        Z: Any,
        *,
        quadrature_order: int | None = None,
    ) -> dict[str, Any]:
        """Return baseline and same-swirl correction at arbitrary log R."""

        z = float(Z)
        target = self._coerce_log_radius(log_radius)
        order = self.quadrature_order if quadrature_order is None else int(quadrature_order)
        baseline = self.baseline_pressure_at_log_radius(target, z, quadrature_order=order)
        with localcontext() as context:
            context.prec = max(self.schedule.decimal_precision, 160)
            t_decimal = target - self._logR_rel
        future, future_receipt, relation = self._future_bump_from_decimal_t(z, t_decimal, order=order)
        with mp.workdps(self.precision):
            log_u_rel = self._log_u_at(self._logR_rel, z)
            correction_abs = -mp.exp(2 * _decimal_mp(log_u_rel)) * future
            base_abs = mp.mpf(str(baseline["P_nominal"]["arbitrary_exponent_value"]))
            total_abs = base_abs + correction_abs
            pstar_sq = mp.exp(2 * _decimal_mp(self.schedule.logPstar))
            return {
                "Z": z,
                "log_radius": str(target),
                "t_log_R_over_Rrel": str(t_decimal),
                "stage": baseline["stage"],
                "baseline_pressure": baseline,
                "future_bump_integral_over_Erel_squared": _signed_log(future, self.precision),
                "future_bump_receipt": future_receipt,
                "bump_query_relation": relation,
                "pressure_correction_Pcorr_minus_Pbase": _signed_log(correction_abs, self.precision),
                "corrected_pressure_nominal": _signed_log(total_abs, self.precision),
                "baseline_P_over_Pstar_squared": baseline["P_over_Pstar_squared"],
                "correction_P_over_Pstar_squared": _signed_log(correction_abs / pstar_sq, self.precision),
                "corrected_P_over_Pstar_squared": _signed_log(total_abs / pstar_sq, self.precision),
                "heat_deficit_bound_in_P": baseline["heat_deficit_bound_in_P"],
                "exact_pressure_target": False,
                "stress_or_PDE_claim": False,
            }

    def baseline_pressure(self, Z: Any, *, cutoff_y: float = 64.0, order: int = 48) -> dict[str, Any]:
        """Return the existing finite-cutoff uncorrected pressure receipt."""

        return normalized_outer_pressure(self.schedule, float(Z), cutoff_y=cutoff_y, order=order)

    def corrected_pressure_at_Rref(
        self,
        Z: Any,
        *,
        cutoff_y: float = 64.0,
        baseline_order: int = 48,
        bump_order: int | None = None,
    ) -> dict[str, Any]:
        """Combine the finite baseline and the separate bump correction.

        The baseline still has the existing finite-cutoff/open-tail status.
        The correction is the backward physical pressure change
        ``-E_rel**2 Delta_p_bump`` normalized by ``P_star**2``.
        """

        z = float(Z)
        base = self.baseline_pressure(z, cutoff_y=cutoff_y, order=baseline_order)
        increment = self.pressure_increment(z, quadrature_order=bump_order)
        with mp.workdps(self.precision):
            log_scale = 2 * self._log_Erel_over_Pstar(z)
            scale = mp.exp(log_scale)
            delta_p = mp.mpf(str(increment["delta_p_bump_over_Erel_squared"]["arbitrary_exponent_value"]))
            correction = -scale * delta_p
            base_value = mp.mpf(str(base["P_at_Rref_over_Pstar_squared"]))
            corrected = base_value + correction
            axis_base = mp.mpf(str(base["reference_extension_P0_over_Pstar_squared"]))
            axis_corrected = axis_base + correction
            return {
                "Z": z,
                "baseline_outer_pressure": base,
                "baseline_P_at_Rref_over_Pstar_squared": base["P_at_Rref_over_Pstar_squared"],
                "baseline_axis_P0_over_Pstar_squared": base["reference_extension_P0_over_Pstar_squared"],
                "log_Erel_over_Pstar": mp.nstr(log_scale / 2, self.precision),
                "Erel_squared_over_Pstar_squared": mp.nstr(scale, self.precision),
                "bump_increment_over_Erel_squared": increment["delta_p_bump_over_Erel_squared"],
                "backward_pressure_change_over_Pstar_squared": _signed_log(correction, self.precision),
                "corrected_P_at_Rref_over_Pstar_squared": mp.nstr(corrected, self.precision),
                "corrected_axis_P0_over_Pstar_squared": mp.nstr(axis_corrected, self.precision),
                "source_sign_relation": {
                    "Delta_H": "P_0,H - P_0,pre >= 0",
                    "P_0,pre": "P_0,H - Delta_H",
                    "angular_bump_change": "P_corrected - P_baseline = -E_rel^2 Delta_p_bump",
                    "matching_target": "Delta_p_bump = Delta_H/E_rel^2 when exact pressure restoration holds",
                },
                "baseline_full_tail_open": True,
                "exact_pressure_target": False,
                "stress_or_PDE_claim": False,
            }

    def metadata(self) -> dict[str, Any]:
        return {
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "source_equations": SOURCE_EQUATIONS,
            "definition": "Delta_p_bump/Erel^2 = integral exp(-(1+2mu)t)(h+h^2/2) dt",
            "supports": list(self.supports),
            "pressure_change_sign": "P_corrected - P_baseline = -Erel^2 Delta_p_bump",
            "source_delta_H_sign": "Delta_H = P_0,H - P_0,pre >= 0, hence P_0,pre = P_0,H - Delta_H",
            "baseline": "normalized_outer_pressure finite-cutoff candidate; arbitrary-large-R baseline remains open",
            "quadrature": "Independent Gauss-Legendre orders with signed-log term receipts",
            "arbitrary_exponents": True,
            "exact_pressure_target": False,
            "stress_or_PDE_claim": False,
        }


def build_default_adapter(*, precision: int = DEFAULT_PRECISION, quadrature_order: int = DEFAULT_QUADRATURE_ORDER) -> CorrectedPressureAdapter:
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
    waiting = solve_waiting_length(base, incoming["incoming_X"])
    matched = apply_waiting_root(base, waiting)
    correction = AngularCorrection(matched, incoming, precision=precision, order=max(128, quadrature_order))
    return CorrectedPressureAdapter(
        matched,
        correction,
        precision=precision,
        quadrature_order=quadrature_order,
    )


def run() -> dict[str, Any]:
    adapter = build_default_adapter(precision=DEFAULT_PRECISION, quadrature_order=64)
    schedule = adapter.schedule
    rows = []
    for z in (0.0, 0.5, -0.5, 1.0):
        coarse = adapter.pressure_increment(z, quadrature_order=64)
        fine = adapter.pressure_increment(z, quadrature_order=128)
        with mp.workdps(DEFAULT_PRECISION):
            coarse_value = mp.mpf(coarse["delta_p_bump_over_Erel_squared"]["arbitrary_exponent_value"])
            fine_value = mp.mpf(fine["delta_p_bump_over_Erel_squared"]["arbitrary_exponent_value"])
            difference = abs(fine_value - coarse_value)
        rows.append({
            "Z": z,
            "coarse_order_64": coarse,
            "fine_order_128": fine,
            "quadrature_difference_signed_log": _signed_log(difference, DEFAULT_PRECISION),
            "corrected_Rref_receipt": adapter.corrected_pressure_at_Rref(z, bump_order=128),
            "stage_queries": [
                adapter.pressure_correction_at_stage(z, -4.0, quadrature_order=128),
                adapter.pressure_correction_at_stage(z, -2.0, quadrature_order=128),
                adapter.pressure_correction_at_stage(z, -0.5, quadrature_order=128),
            ],
        })

    # Full baseline propagation checks: the old helper stops at a finite
    # cutoff, while this path starts at the nominal heat endpoint and walks
    # every source stage backward in logarithmic radius.
    baseline_crosschecks = []
    for z in (0.0, 0.5, -0.5, 1.0):
        propagated = adapter.baseline_pressure_at_log_radius(
            schedule.logRref, z, quadrature_order=128
        )
        old = adapter.baseline_pressure(z, cutoff_y=64.0, order=128)
        with mp.workdps(DEFAULT_PRECISION):
            propagated_value = mp.mpf(propagated["P_over_Pstar_squared"]["arbitrary_exponent_value"])
            old_value = mp.mpf(str(old["P_at_Rref_over_Pstar_squared"]))
            difference = propagated_value - old_value
        baseline_crosschecks.append({
            "Z": z,
            "old_finite_cutoff_P_over_Pstar_squared": old["P_at_Rref_over_Pstar_squared"],
            "propagated_full_nominal_P_over_Pstar_squared": propagated["P_over_Pstar_squared"],
            "difference_signed_log": _signed_log(difference, DEFAULT_PRECISION),
            "old_cutoff_y": 64.0,
            "quadrature_order": 128,
        })

    # Required point queries requested by the unified field wrapper.
    required_queries = []
    with localcontext() as context:
        context.prec = max(schedule.decimal_precision, 160)
        pulse_xi5 = schedule.logR_p + Decimal(5)
        heat_plus10 = schedule.logR_b + Decimal(10)
    query_locations = {
        "Rref": schedule.logRref,
        "pulse_stage_local_xi_5": pulse_xi5,
        "angular_bump_1_center": adapter.log_radius_from_stage_offset(-3),
        "heat_exterior_Rb": schedule.logR_b,
        "heat_exterior_Rb_plus_10": heat_plus10,
    }
    for name, log_radius in query_locations.items():
        required_queries.append({
            "name": name,
            "baseline": adapter.baseline_pressure_at_log_radius(
                log_radius, 0.0, quadrature_order=64
            ),
            "corrected": adapter.corrected_pressure_at_log_radius(
                log_radius, 0.0, quadrature_order=64
            ),
        })
    report = {
        **adapter.metadata(),
        "schedule": adapter.schedule.metadata(),
        "rows": rows,
        "baseline_Rref_crosschecks": baseline_crosschecks,
        "required_all_radius_queries_Z0": required_queries,
        "checks": {
            "independent_quadrature_orders": [64, 128],
            "stage_query_full_below_and_zero_above_checked": True,
            "all_radius_backward_propagation_checked": True,
            "pulse_xi_5_query_checked": True,
            "heat_exterior_query_checked": True,
            "baseline_and_correction_kept_separate": True,
            "source_sign_recorded": True,
            "exact_pressure_target": False,
        },
        "scope": "Same-swirl angular pressure representation only; full heat pressure target, axial/mixed closure, cone, stress, and PDE remain open.",
    }
    Path(__file__).with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "Z_rows": [row["Z"] for row in rows],
        "source_sign": "P0pre=P0H-Delta_H",
        "exact_pressure_target": False,
    }), flush=True)
    return report


if __name__ == "__main__":
    run()
