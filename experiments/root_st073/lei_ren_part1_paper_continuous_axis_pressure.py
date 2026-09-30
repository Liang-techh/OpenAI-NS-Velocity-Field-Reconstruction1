"""Continuous pre-``Rv`` axis-pressure datum for the nonlinear core.

This adapter is deliberately narrower than :class:`AxisPressureJets`.  It
requires the shared ``PaperOuterSchedule`` to have the continuous angular
provider installed already, then integrates the common ``Z=0`` pre-``Rv``
angular amplitude.  Variable transition stages use MP Gauss nodes and the
constant stages use their exact exponential primitives.  The post-``Rv``
pressure tail is left out of the axis jet and is reported as uncertified.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
from pathlib import Path
import sys
from typing import Any

import mpmath as mp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SRC = ROOT / "src"
for path in (SRC, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_paper_axis_pressure_jets import (  # noqa: E402
    AxisPressureJets,
    REFERENCE_EXTENSION_INTEGRAL,
)


class ContinuousAxisPressureJets(AxisPressureJets):
    """Axis pressure jets from the shared continuous angular pre-``Rv`` data.

    The returned anchor is normalized by ``Pstar**2``.  With
    ``F = Utheta / sqrt(2 R)``, the exact ``Z=0`` separated prefix is

    ``K = -2.5 - integral(F**2 dR / Pstar**2)``.

    The ``-2.5`` term is the declared reference-to-axis extension.  No
    post-``Rv`` value or derivative is folded into this dominant jet.
    """

    _PRE_RV_STAGES = (
        ("slope_transition_ref", None),
        ("axial_turnoff", "-0.5"),
        ("slope_transition_mu", None),
        ("power_buffer", "power"),
        ("pulse_reserved", "power"),
    )

    def __init__(
        self,
        profile: Any,
        *,
        quadrature_order: int = 192,
        R_a: Any | None = None,
    ) -> None:
        schedule = profile.schedule
        if getattr(schedule, "_continuous_angular_provider", None) is None:
            raise ValueError(
                "ContinuousAxisPressureJets requires an installed continuous "
                "angular schedule"
            )
        super().__init__(
            profile,
            quadrature_order=quadrature_order,
            R_a=R_a,
        )
        with mp.workdps(self.precision):
            self.mu = mp.mpf(str(self.schedule.mu))
        self._mp_nodes_cache: dict[int, tuple[tuple[mp.mpf, mp.mpf], ...]] = {}

    def _mp_gauss_nodes(
        self, order: int
    ) -> tuple[tuple[mp.mpf, mp.mpf], ...]:
        order = int(order)
        cached = self._mp_nodes_cache.get(order)
        if cached is not None:
            return cached
        with mp.workdps(self.precision):
            nodes, weights = mp.gauss_quadrature(order, "legendre")
            cached = tuple(
                (
                    (mp.mpf(node) + 1) / 2,
                    mp.mpf(weight) / 2,
                )
                for node, weight in zip(nodes, weights)
            )
        self._mp_nodes_cache[order] = cached
        return cached

    def _log_u_over_pstar(self, offset: Any) -> mp.mpf:
        """Return ``log(Utheta/Pstar)`` at an exact schedule offset, ``Z=0``."""

        log_radius = self._log_at_offset(offset)
        row = self.schedule.at_log_radius(log_radius, 0)
        return mp.mpf(str(row["log_angular_amplitude"])) - mp.mpf(
            str(self.schedule.logPstar)
        )

    def _stage_integral(
        self,
        name: str,
        slope_kind: str | None,
        left: Decimal,
        right: Decimal,
        order: int,
    ) -> tuple[mp.mpf, dict[str, Any]]:
        with mp.workdps(self.precision):
            left_mp = mp.mpf(str(left))
            right_mp = mp.mpf(str(right))
            length = right_mp - left_mp
            if length <= 0:
                return mp.mpf("0"), {
                    "stage": name,
                    "left": str(left),
                    "right": str(right),
                    "method": "empty",
                    "integral": "0",
                }

            if slope_kind is None:
                nodes = self._mp_gauss_nodes(order)
                total = mp.fsum(
                    weight
                    * mp.exp(
                        2
                        * self._log_u_over_pstar(
                            left_mp + length * node
                        )
                    )
                    for node, weight in nodes
                )
                value = length * total / 2
                receipt = {
                    "stage": name,
                    "left": str(left),
                    "right": str(right),
                    "method": "MP Gauss Legendre",
                    "quadrature_order": order,
                    "integral": mp.nstr(value, self.precision),
                }
                return value, receipt

            if slope_kind == "-0.5":
                slope = mp.mpf("-0.5")
            elif slope_kind == "power":
                slope = -mp.mpf("0.5") - self.mu
            else:
                raise ValueError(f"unknown pre-Rv slope kind: {slope_kind}")

            rate = 2 * slope
            log_ratio_left = self._log_u_over_pstar(left)
            if rate == 0:
                atom = length
            else:
                atom = mp.expm1(rate * length) / rate
            value = mp.mpf("0.5") * mp.exp(2 * log_ratio_left) * atom
            receipt = {
                "stage": name,
                "left": str(left),
                "right": str(right),
                "method": "exact exponential atom",
                "slope": mp.nstr(slope, self.precision),
                "integral": mp.nstr(value, self.precision),
            }
            return value, receipt

    def _preflatten_integral(self, order: int) -> dict[str, Any]:
        """Integrate the shared continuous ``Z=0`` prefix from ``Rref`` to ``Rv``."""

        order = int(order)
        if order in self._integral_cache:
            return self._integral_cache[order]
        bounds = self.schedule._make_stage_bounds()
        with mp.workdps(self.precision):
            intervals: list[dict[str, Any]] = []
            values: list[mp.mpf] = []
            for name, slope_kind in self._PRE_RV_STAGES:
                left, right = bounds[name]
                if right is None or right > self.schedule.y_v:
                    raise ValueError(f"pre-Rv stage {name} has an invalid endpoint")
                value, receipt = self._stage_integral(
                    name, slope_kind, left, right, order
                )
                values.append(value)
                intervals.append(receipt)

            early = mp.fsum(values[:3])
            power = mp.fsum(values[3:])
            total = early + power
            result = {
                "quadrature_order": order,
                "early_0_to_yw": mp.nstr(early, self.precision),
                "power_yw_to_yv": mp.nstr(power, self.precision),
                "total_Rref_to_Rv": mp.nstr(total, self.precision),
                "intervals": intervals,
                "slope_transition_ref": "finite MP quadrature",
                "axial_turnoff_slope": "-0.5 exact",
                "slope_transition_mu": "finite MP quadrature",
                "power_slope": mp.nstr(-mp.mpf("0.5") - self.mu, self.precision),
                "Z_separation": "exact factor (1+Z^2)^-2 before Rv",
                "normalization": (
                    "integral Utheta^2/(2 Pstar^2) dy = "
                    "integral F^2 dR/Pstar^2"
                ),
                "continuous_angular_schedule_required": True,
                "z0_prefix": True,
                "post_Rv_tail_omitted": True,
                "post_Rv_tail_certified": False,
            }
        self._integral_cache[order] = result
        return result

    def _anchor(self, order: int) -> dict[str, Any]:
        order = int(order)
        if order in self._anchor_cache:
            return self._anchor_cache[order]
        with mp.workdps(self.precision):
            integral = self._preflatten_integral(order)
            K = -REFERENCE_EXTENSION_INTEGRAL - mp.mpf(
                integral["total_Rref_to_Rv"]
            )
            result = {
                "quadrature_order": order,
                "pressure_units": "P0_over_Pstar_squared",
                "normalization": "all scalar pressure values in this receipt are P/Pstar^2",
                "anchor_K_actual_Z0": mp.nstr(K, self.precision),
                "anchor_K_integral_split": mp.nstr(K, self.precision),
                "anchor_difference": "0",
                "full_P0_Z0": None,
                "P_Rref_Z0": None,
                "tail_P_Rv_Z0": None,
                "full_reference_receipt_Z0": None,
                "tail_receipt_Z0": None,
                "preflatten_integral": integral,
                "anchor_definition": (
                    "-2.5 minus continuous pre-Rv pressure integral; exact Z0 prefix"
                ),
                "post_Rv_pressure_tail_omitted": True,
                "post_Rv_pressure_tail_certified": False,
                "tail_derivative_status": "post-Rv tail omitted; no derivative claim",
                "full_pressure_exact": False,
            }
        self._anchor_cache[order] = result
        return result


__all__ = ["ContinuousAxisPressureJets"]
