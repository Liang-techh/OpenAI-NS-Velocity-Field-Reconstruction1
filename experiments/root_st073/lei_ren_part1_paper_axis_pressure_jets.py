"""Analytic dominant axis-pressure jets with an explicit future-tail receipt.

For the current source schedule, before ``R_v`` the angular amplitude has the
exact separated form ``Utheta(R,Z)=E(R)/(1+Z^2)``.  Therefore the pressure
primitive from the reference branch through ``R_v`` is a single constant
times ``(1+Z^2)^-2``.  This module computes that constant by splitting at
``R_v`` and anchors it at the actual same-profile pressure at ``Z=0``.  The
post-``R_v`` pressure is retained as a separate signed-log tail value and
bound; it is not differentiated or folded into the analytic Taylor jet.

This is an axis-pressure representation for continuation inputs.  It does
not certify a regular core, pressure-target equality, PDE, cone, global
energy, or theorem conditions.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
from functools import lru_cache
import json
from pathlib import Path
import sys
from typing import Any

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SRC = ROOT / "src"
for path in (SRC, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_paper_corrected_profile import CorrectedSourceProfile  # noqa: E402


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
REFERENCE_EXTENSION_INTEGRAL = mp.mpf("2.5")


def _signed_log(value: Any, precision: int) -> dict[str, Any]:
    value_mp = mp.mpf(value)
    if value_mp == 0:
        return {"sign": 0, "log_abs": None, "arbitrary_exponent_value": "0"}
    return {
        "sign": 1 if value_mp > 0 else -1,
        "log_abs": mp.nstr(mp.log(abs(value_mp)), precision),
        "arbitrary_exponent_value": mp.nstr(value_mp, precision),
    }


def _from_signed_log(row: dict[str, Any]) -> mp.mpf:
    if not row or not row.get("sign") or row.get("log_abs") is None:
        return mp.mpf("0")
    value = mp.exp(mp.mpf(str(row["log_abs"])))
    return value if int(row["sign"]) > 0 else -value


def _relative(defect: mp.mpf, scale: mp.mpf, precision: int) -> dict[str, Any]:
    return {
        "absolute": _signed_log(defect, precision),
        "relative": mp.nstr(abs(defect / scale), precision) if scale else None,
    }


class AxisPressureJets:
    """Shared-profile axis pressure and analytic dominant Taylor jets.

    ``profile`` must be a ``CorrectedSourceProfile`` whose pressure adapter,
    schedule, and angular correction are shared.  All returned pressure
    values are normalized by ``Pstar**2`` unless a receipt explicitly says
    otherwise.  ``R_a`` is optional continuation metadata and is never used
    to invent a core pressure.
    """

    def __init__(
        self,
        profile: CorrectedSourceProfile,
        *,
        quadrature_order: int = 128,
        R_a: Any | None = None,
    ) -> None:
        if int(quadrature_order) != quadrature_order or quadrature_order < 16:
            raise ValueError("quadrature_order must be an integer at least 16")
        self.profile = profile
        self.schedule = profile.schedule
        self.pressure = profile.pressure
        if self.pressure.schedule is not self.schedule or self.pressure.correction is not profile.angular:
            raise ValueError("profile pressure, schedule, and angular correction must be shared")
        self.precision = int(profile.precision)
        self.quadrature_order = int(quadrature_order)
        self.R_a = None if R_a is None else mp.mpf(str(R_a))
        self._integral_cache: dict[int, dict[str, Any]] = {}
        self._anchor_cache: dict[int, dict[str, Any]] = {}

    def _log_at_offset(self, offset: Any) -> Decimal:
        with localcontext() as context:
            context.prec = max(int(self.schedule.decimal_precision), self.precision)
            return self.schedule.logRref + Decimal(str(offset))

    def _preflatten_integral(self, order: int) -> dict[str, Any]:
        """Compute ``int_Rref^Rv F^2 dR/Pstar^2`` using exact Z separation."""

        order = int(order)
        if order in self._integral_cache:
            return self._integral_cache[order]
        nodes, weights = leggauss(order)
        y_d = self.schedule.y_d
        y_w = self.schedule.y_w
        edges = (Decimal(0), Decimal(1), y_d, y_w)
        with mp.workdps(self.precision):
            early = mp.mpf("0")
            interval_receipts = []
            for left, right in zip(edges[:-1], edges[1:]):
                length = right - left
                half = mp.mpf(str(length)) / 2
                midpoint = (mp.mpf(str(left)) + mp.mpf(str(right))) / 2
                total = mp.mpf("0")
                for node, weight in zip(nodes, weights):
                    local = midpoint + half * mp.mpf(str(float(node)))
                    log_radius = self._log_at_offset(local)
                    row = self.schedule.at_log_radius(log_radius, 0.0)
                    log_ratio = mp.mpf(str(row["log_angular_amplitude"])) - mp.mpf(
                        str(self.schedule.logPstar)
                    )
                    total += mp.mpf(str(float(weight))) * mp.exp(2 * log_ratio)
                interval = half * total / 2
                early += interval
                interval_receipts.append(
                    {
                        "left": str(left),
                        "right": str(right),
                        "integral": mp.nstr(interval, self.precision),
                    }
                )

            # On [y_w,y_v] the source has the exact slope
            # Utheta/Pstar = const * exp(-(1/2+mu)(y-y_w)).
            log_u_w = self.schedule.at_log_radius(self._log_at_offset(y_w), 0.0)[
                "log_angular_amplitude"
            ]
            u_w_squared = mp.exp(
                2
                * (
                    mp.mpf(str(log_u_w))
                    - mp.mpf(str(self.schedule.logPstar))
                )
            )
            slope = mp.mpf("0.5") + mp.mpf(str(self.schedule.mu))
            length = mp.mpf(str(self.schedule.y_v - y_w))
            contraction = -mp.expm1(-2 * slope * length)
            power = u_w_squared * contraction / (4 * slope)
            total = early + power
            result = {
                "quadrature_order": order,
                "early_0_to_yw": mp.nstr(early, self.precision),
                "power_yw_to_yv": mp.nstr(power, self.precision),
                "total_Rref_to_Rv": mp.nstr(total, self.precision),
                "intervals": interval_receipts,
                "power_slope": mp.nstr(slope, self.precision),
                "power_contraction": mp.nstr(contraction, self.precision),
                "Z_separation": "exact factor (1+Z^2)^-2 before Rv",
                "no_physical_radius_materialized": True,
            }
        self._integral_cache[order] = result
        return result

    def _tail_pressure(self, Z: Any, order: int) -> dict[str, Any]:
        return self.pressure.corrected_pressure_at_log_radius(
            self.schedule.logR_v,
            Z,
            quadrature_order=order,
        )

    def preflatten_integral(self, *, quadrature_order: int | None = None) -> dict[str, Any]:
        """Return the exact-Z-separated reference-to-``Rv`` integral receipt."""

        order = self.quadrature_order if quadrature_order is None else int(quadrature_order)
        if order < 16:
            raise ValueError("quadrature_order must be at least 16")
        return self._preflatten_integral(order)

    def anchor_receipt(self, *, quadrature_order: int | None = None) -> dict[str, Any]:
        """Return the actual ``Z=0`` full-axis-minus-tail anchor receipt."""

        order = self.quadrature_order if quadrature_order is None else int(quadrature_order)
        if order < 16:
            raise ValueError("quadrature_order must be at least 16")
        return self._anchor(order)

    def _full_reference_pressure(self, Z: Any, order: int) -> dict[str, Any]:
        return self.pressure.corrected_pressure_at_log_radius(
            self.schedule.logRref,
            Z,
            quadrature_order=order,
        )

    def _anchor(self, order: int) -> dict[str, Any]:
        order = int(order)
        if order in self._anchor_cache:
            return self._anchor_cache[order]
        with mp.workdps(self.precision):
            tail = self._tail_pressure(0.0, order)
            full = self._full_reference_pressure(0.0, order)
            tail_value = _from_signed_log(tail["corrected_P_over_Pstar_squared"])
            # The adapter value at Rref is P(Rref).  The regular-axis
            # extension contributes the exact reference integral 2.5 at
            # Z=0, so the axis value is P(Rref)-2.5.
            reference_value = _from_signed_log(full["corrected_P_over_Pstar_squared"])
            full_value = reference_value - REFERENCE_EXTENSION_INTEGRAL
            anchor = full_value - tail_value
            integral = self._preflatten_integral(order)
            integral_anchor = -REFERENCE_EXTENSION_INTEGRAL - mp.mpf(
                integral["total_Rref_to_Rv"]
            )
            result = {
                "quadrature_order": order,
                "pressure_units": "P0_over_Pstar_squared",
                "normalization": "all scalar pressure values in this receipt are P/Pstar^2",
                "anchor_K_actual_Z0": mp.nstr(anchor, self.precision),
                "anchor_K_integral_split": mp.nstr(integral_anchor, self.precision),
                "anchor_difference": mp.nstr(anchor - integral_anchor, self.precision),
                "full_P0_Z0": _signed_log(full_value, self.precision),
                "P_Rref_Z0": _signed_log(reference_value, self.precision),
                "tail_P_Rv_Z0": _signed_log(tail_value, self.precision),
                "full_reference_receipt_Z0": full,
                "tail_receipt_Z0": tail,
                "preflatten_integral": integral,
                "anchor_definition": "actual full P0(Z=0) minus same-profile tail P(Rv,Z=0)",
                "tail_derivative_status": "not expanded or derivative-certified",
            }
        self._anchor_cache[order] = result
        return result

    def dominant_taylor(self, degree: int = 8, *, quadrature_order: int | None = None) -> dict[str, Any]:
        """Return arbitrary-degree coefficients of ``K/(1+Z^2)^2``."""

        if int(degree) != degree or degree < 0:
            raise ValueError("degree must be a nonnegative integer")
        order = self.quadrature_order if quadrature_order is None else int(quadrature_order)
        anchor = self._anchor(order)
        with mp.workdps(self.precision):
            K = mp.mpf(anchor["anchor_K_actual_Z0"])
            coefficients = []
            for n in range(int(degree) + 1):
                if n % 2:
                    value = mp.mpf("0")
                else:
                    k = n // 2
                    value = K * ((-1) ** k) * (k + 1)
                coefficients.append(mp.nstr(value, self.precision))
            return {
                "degree": int(degree),
                "coefficients": coefficients,
                "anchor_K_actual_Z0": mp.nstr(K, self.precision),
                "pressure_units": "P0_over_Pstar_squared",
                "normalization": "all coefficients and anchor values are P/Pstar^2",
                "Pstar_squared": mp.nstr(
                    mp.exp(2 * mp.mpf(str(self.schedule.logPstar))), self.precision
                ),
                "formula": "coefficient[2k]=K*(-1)^k*(k+1), coefficient[2k+1]=0",
                "future_tail_omitted_from_jet": True,
                "future_tail_derivative_bound": None,
                "anchor": anchor,
            }

    def physical_taylor(
        self, degree: int = 8, *, quadrature_order: int | None = None
    ) -> dict[str, Any]:
        """Return the same dominant jet rescaled to physical pressure units.

        ``dominant_taylor`` is deliberately normalized by ``Pstar**2`` for
        stable source-scale arithmetic.  This helper performs the explicit
        conversion so callers cannot accidentally treat the normalized
        anchor as physical pressure.
        """

        normalized = self.dominant_taylor(
            degree, quadrature_order=quadrature_order
        )
        with mp.workdps(self.precision):
            pstar_squared = mp.exp(2 * mp.mpf(str(self.schedule.logPstar)))
            physical_coefficients = [
                mp.nstr(mp.mpf(value) * pstar_squared, self.precision)
                for value in normalized["coefficients"]
            ]
            physical_anchor = mp.mpf(normalized["anchor_K_actual_Z0"]) * pstar_squared
        return {
            "degree": normalized["degree"],
            "coefficients": physical_coefficients,
            "anchor_K_actual_Z0": mp.nstr(physical_anchor, self.precision),
            "pressure_units": "physical_P",
            "source_normalization": "physical coefficients = normalized coefficients * Pstar^2",
            "Pstar_squared": mp.nstr(pstar_squared, self.precision),
            "formula": normalized["formula"],
            "future_tail_omitted_from_jet": normalized["future_tail_omitted_from_jet"],
            "future_tail_derivative_bound": normalized["future_tail_derivative_bound"],
            "normalized_jet": normalized,
        }

    def pressure_receipt(
        self,
        Z: Any,
        *,
        quadrature_order: int | None = None,
        include_direct_reference: bool = True,
    ) -> dict[str, Any]:
        """Return split pressure in ``P/Pstar**2`` units and a tail receipt."""

        z = float(Z)
        if not np.isfinite(z) or abs(z) > 1:
            raise ValueError("Z must lie in [-1,1]")
        order = self.quadrature_order if quadrature_order is None else int(quadrature_order)
        anchor = self._anchor(order)
        with mp.workdps(self.precision):
            K = mp.mpf(anchor["anchor_K_actual_Z0"])
            tail = self._tail_pressure(z, order)
            tail_value = _from_signed_log(tail["corrected_P_over_Pstar_squared"])
            dominant = K / (1 + mp.mpf(str(z)) ** 2) ** 2
            full_value = dominant + tail_value
            direct = self._full_reference_pressure(z, order) if include_direct_reference else None
            direct_reference_value = (
                _from_signed_log(direct["corrected_P_over_Pstar_squared"])
                if direct is not None
                else mp.mpf("0")
            )
            reference_extension = REFERENCE_EXTENSION_INTEGRAL / (
                1 + mp.mpf(str(z)) ** 2
            ) ** 2
            direct_value = direct_reference_value - reference_extension
            pstar_squared = mp.exp(2 * mp.mpf(str(self.schedule.logPstar)))
            heat_bound = _from_signed_log(tail["heat_deficit_bound_in_P"]) / pstar_squared
            bump_change = _from_signed_log(tail["correction_P_over_Pstar_squared"])
            return {
                "Z": z,
                "quadrature_order": order,
                "pressure_units": "P0_over_Pstar_squared",
                "Pstar_squared": mp.nstr(pstar_squared, self.precision),
                "dominant_P0": _signed_log(dominant, self.precision),
                "future_tail_P_Rv": _signed_log(tail_value, self.precision),
                "full_P0": _signed_log(full_value, self.precision),
                "direct_full_reference_P": _signed_log(direct_value, self.precision)
                if direct is not None
                else None,
                "direct_P_at_Rref": _signed_log(direct_reference_value, self.precision)
                if direct is not None
                else None,
                "decomposition_error_vs_direct": _relative(
                    full_value - direct_value,
                    direct_value,
                    self.precision,
                )
                if direct is not None
                else None,
                "future_tail": {
                    "corrected_pressure_over_Pstar_squared": tail[
                        "corrected_P_over_Pstar_squared"
                    ],
                    "baseline_pressure_over_Pstar_squared": tail[
                        "baseline_P_over_Pstar_squared"
                    ],
                    "angular_bump_correction_over_Pstar_squared": tail[
                        "correction_P_over_Pstar_squared"
                    ],
                    "heat_deficit_bound_over_Pstar_squared": _signed_log(
                        heat_bound,
                        self.precision,
                    ),
                    "same_profile_adapter_receipt": tail,
                    "derivative_status": "value/bound only; no derivative bound claimed",
                },
                "dominant_formula": "K/(1+Z^2)^2 with K anchored at actual Z=0 full P0 minus tail",
                "preflatten_split": anchor["preflatten_integral"],
                "full_pressure_exact": False,
                "tail_derivative_bound": None,
            }

    def taylor_coefficients(self, degree: int = 8, *, quadrature_order: int | None = None) -> dict[str, Any]:
        """Alias for callers requesting analytic dominant pressure jets."""

        return self.dominant_taylor(degree, quadrature_order=quadrature_order)


def _make_source_profile(*, precision: int = 120, order: int = 128) -> tuple[CorrectedSourceProfile, mp.mpf, mp.mpf]:
    with mp.workdps(precision):
        Lambda = mp.mpf("1e18")
        logCstar = 2 * mp.log(Lambda)
        logPstar = mp.mpf("14")
        logRref = mp.log(110) + 10 * (logCstar + logPstar)
    from lei_ren_part1_paper_outer import PaperOuterSchedule
    from lei_ren_part1_paper_axial_energy_tail import build_default_tail

    base = PaperOuterSchedule(
        logPstar=mp.nstr(logPstar, precision),
        logRref=mp.nstr(logRref, precision),
        delta="1e-32",
        Md=".5",
        c_mu=".001",
        c_delta=".001",
        c_epsilon=".01",
    )
    tail = build_default_tail(
        precision=precision,
        quadrature_order=order,
        schedule=base,
        match_waiting=True,
    )
    profile = CorrectedSourceProfile(
        precision=precision,
        order=order,
        schedule=tail.schedule,
        angular=tail.correction,
        tail=tail,
    )
    return profile, Lambda, logCstar


def run() -> dict[str, Any]:
    precision = 120
    order = 128
    profile, Lambda, logCstar = _make_source_profile(precision=precision, order=order)
    jets = AxisPressureJets(profile, quadrature_order=order, R_a=4 / Lambda)
    anchor = jets.anchor_receipt(quadrature_order=order)
    anchor_coarse = jets.anchor_receipt(quadrature_order=64)
    refinement = {
        "order_64": jets._preflatten_integral(64),
        "order_128": jets._preflatten_integral(128),
    }
    with mp.workdps(precision):
        i64 = mp.mpf(refinement["order_64"]["total_Rref_to_Rv"])
        i128 = mp.mpf(refinement["order_128"]["total_Rref_to_Rv"])
        refinement_difference = abs(i128 - i64)
    rows = [jets.pressure_receipt(z, quadrature_order=order) for z in (0.0, 0.3, 0.7)]
    tail_refinement = []
    for z in (0.0, 0.3, 0.7):
        tail_64 = jets.pressure_receipt(z, quadrature_order=64)
        tail_128 = jets.pressure_receipt(z, quadrature_order=128)
        with mp.workdps(precision):
            tail_value_64 = _from_signed_log(tail_64["future_tail_P_Rv"])
            tail_value_128 = _from_signed_log(tail_128["future_tail_P_Rv"])
            tail_abs_difference = abs(tail_value_128 - tail_value_64)
            tail_rel_difference = (
                abs((tail_value_128 - tail_value_64) / tail_value_128)
                if tail_value_128
                else mp.mpf("0")
            )
            full_value_64 = _from_signed_log(tail_64["full_P0"])
            full_value_128 = _from_signed_log(tail_128["full_P0"])
            full_abs_difference = abs(full_value_128 - full_value_64)
            full_rel_difference = (
                abs((full_value_128 - full_value_64) / full_value_128)
                if full_value_128
                else mp.mpf("0")
            )
        tail_refinement.append(
            {
                "Z": z,
                "tail_order_64": tail_64["future_tail_P_Rv"],
                "tail_order_128": tail_128["future_tail_P_Rv"],
                "tail_absolute_difference": _signed_log(
                    tail_abs_difference, precision
                ),
                "tail_relative_difference": mp.nstr(tail_rel_difference, precision),
                "full_P0_absolute_difference": _signed_log(
                    full_abs_difference, precision
                ),
                "full_P0_relative_difference": mp.nstr(full_rel_difference, precision),
            }
        )
    report = {
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "parameters": {
            "j": "0.02",
            "Lambda": mp.nstr(Lambda, precision),
            "logCstar": mp.nstr(logCstar, precision),
            "logCstar_lower_bound": "2 log Lambda; A_Omega lower bound set to zero only",
            "R_a": mp.nstr(4 / Lambda, precision),
            "logRref": str(profile.schedule.logRref),
        },
        "profile_construction": profile.metadata()["construction"],
        "anchor": anchor,
        "dominant_taylor_degree_12": jets.dominant_taylor(12),
        "preflatten_quadrature_refinement": {
            "order_64": refinement["order_64"],
            "order_128": refinement["order_128"],
            "total_integral_absolute_difference": mp.nstr(refinement_difference, precision),
            "anchor_order_64": anchor_coarse["anchor_K_actual_Z0"],
            "anchor_order_128": anchor["anchor_K_actual_Z0"],
            "anchor_absolute_difference": mp.nstr(
                abs(
                    mp.mpf(anchor["anchor_K_actual_Z0"])
                    - mp.mpf(anchor_coarse["anchor_K_actual_Z0"])
                ),
                precision,
            ),
        },
        "rows": rows,
        "tail_quadrature_refinement": tail_refinement,
        "checks": {
            "decomposition_samples": [row["Z"] for row in rows],
            "future_tail_value_and_bound_retained": True,
            "future_tail_derivative_bound_claimed": False,
            "dominant_jet_is_analytic": True,
            "anchor_comes_from_actual_Z0_full_minus_tail": True,
            "full_pressure_exact": False,
            "core_or_theorem_certified": False,
        },
        "scope": (
            "Same-profile axis-pressure split and dominant analytic Taylor jet. "
            "The post-Rv tail remains a signed numerical value with a heat bound; "
            "no derivative bound is asserted for that tail, and no regular-core, "
            "PDE, cone, five-moment, or theorem certification is made."
        ),
    }
    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "anchor_K": anchor["anchor_K_actual_Z0"],
                "integral_anchor_difference": anchor["anchor_difference"],
                "preflatten_refinement": mp.nstr(refinement_difference, precision),
                "decomposition_relative": [
                    (row["Z"], row["decomposition_error_vs_direct"]["relative"])
                    for row in rows
                ],
            },
            default=str,
        )
    )
    return report


if __name__ == "__main__":
    run()
