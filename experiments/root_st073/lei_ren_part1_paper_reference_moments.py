"""Actual five cumulative moments on the temporary source reference.

This module is deliberately a small adapter around ``CorrectedSourceProfile``.
It evaluates the five moments from the same ``Utheta`` and ``Uz`` returned by
that profile; it never replaces an integral by an incoming or boundary target.
The branch ``log(R) <= log(R_ref)`` is evaluated analytically from the paper's
reference profile (4.3).  The finite source interval from ``R_ref`` through
the axial turnoff is integrated in logarithmic radius with arbitrary
precision.  No claim is made about the later outer stages.

The source definitions are (2.21)--(2.22) of Lei--Ren Part I:

    Mtheta   = integral sqrt(2 R) Utheta dR
    Mz       = integral Uz dR
    Mtheta_z = integral sqrt(2 R) Utheta Uz dR
    Mz_theta = integral (Uz^2 - Utheta^2/2) dR
    Mp       = integral Utheta^2/(2 R) dR.

Source: https://arxiv.org/html/2609.35406v1, version 2609.35406v1.
"""

from __future__ import annotations

from decimal import Decimal, InvalidOperation, localcontext
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lei_ren_part1_paper_corrected_profile import CorrectedSourceProfile  # noqa: E402


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
MOMENT_KEYS = ("theta", "z", "theta_z", "z_theta", "p")


def _decimal(value: Any, name: str) -> Decimal:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a finite real number")
    try:
        out = value if isinstance(value, Decimal) else Decimal(str(value))
    except (InvalidOperation, TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not out.is_finite():
        raise ValueError(f"{name} must be finite")
    return out


def _mp(value: Any) -> mp.mpf:
    return mp.mpf(str(value))


def _z(value: Any) -> mp.mpf:
    out = _mp(value)
    if not mp.isfinite(out) or abs(out) > 1:
        raise ValueError("Z must be finite with |Z| <= 1")
    return out


def signed_log(value: Any, *, digits: int = 100) -> dict[str, Any]:
    """Represent an arbitrary-exponent value without binary64 conversion."""

    with mp.workdps(max(80, int(digits))):
        value = _mp(value)
        if value == 0:
            return {
                "sign": 0,
                "log_abs": None,
                "arbitrary_exponent_value": "0",
            }
        return {
            "sign": 1 if value > 0 else -1,
            "log_abs": mp.nstr(mp.log(abs(value)), digits),
            "arbitrary_exponent_value": mp.nstr(value, digits),
        }


def _relative_error(value: mp.mpf, reference: mp.mpf) -> mp.mpf:
    with mp.workdps(max(160, mp.mp.dps)):
        scale = max(abs(reference), mp.mpf("1e-120"))
        return abs(value - reference) / scale


class PaperReferenceMoments:
    """Five actual cumulative moments for one corrected source profile.

    Parameters
    ----------
    profile:
        A ``CorrectedSourceProfile`` instance.  Its ``values`` method is the
        only velocity source used on the finite outer interval.
    precision:
        Decimal digits used by mpmath arithmetic.  The profile precision is
        used as a lower bound.
    quadrature_order:
        Gauss--Legendre order for each finite logarithmic stage.

    The callable interface is ``moments(logR, Z)`` and returns raw mpmath
    values under the five keys expected by the source stress evaluator.
    ``moments_Z`` and ``moments_with_Z`` provide an explicit fourth-order
    centered stencil from the same callable.
    """

    def __init__(
        self,
        profile: CorrectedSourceProfile,
        *,
        precision: int | None = None,
        quadrature_order: int = 64,
    ) -> None:
        if not hasattr(profile, "schedule") or not hasattr(profile, "values"):
            raise TypeError("profile must expose schedule and values(logR, Z)")
        if quadrature_order < 8:
            raise ValueError("quadrature_order must be at least 8")
        self.profile = profile
        self.schedule = profile.schedule
        self.precision = max(
            80,
            int(precision if precision is not None else getattr(profile, "precision", 160)),
            int(getattr(self.schedule, "decimal_precision", 80)),
        )
        self.quadrature_order = int(quadrature_order)
        with mp.workdps(self.precision):
            self._log_rref = _mp(self.schedule.logRref)
            self._log_pstar = _mp(self.schedule.logPstar)
        with mp.workdps(self.precision):
            self._rref = mp.exp(self._log_rref)
            self._axial_turnoff_y = mp.exp(_mp(self.schedule.Md))
        self._y_d_decimal = _decimal(self.schedule.y_d, "schedule.y_d")
        with mp.workdps(self.precision):
            self._y_d = _mp(self._y_d_decimal)

    def __call__(
        self,
        logR: Any,
        Z: Any,
        *,
        quadrature_order: int | None = None,
    ) -> dict[str, mp.mpf]:
        return self.moments_at_log_radius(logR, Z, quadrature_order=quadrature_order)

    def _target_y(self, logR: Any) -> tuple[Decimal, mp.mpf]:
        log_r = _decimal(logR, "logR")
        with localcontext() as ctx:
            ctx.prec = max(self.precision, int(getattr(self.schedule, "decimal_precision", 80)))
            y_decimal = log_r - _decimal(self.schedule.logRref, "schedule.logRref")
        with mp.workdps(self.precision):
            return y_decimal, _mp(y_decimal)

    def _reference_moments(self, y: mp.mpf, z: mp.mpf) -> dict[str, mp.mpf]:
        """Exact (4.3) moments from zero to ``R_ref exp(y)``.

        With ``x=R/R_ref``, ``Utheta=Uref*x^(1/10)`` and ``Uz=4Z``:

        ``Mtheta=(5/8)sqrt(2) Rref^(3/2) Uref*x^(8/5)``;
        ``Mz=4 Z Rref*x``;
        ``Mtheta_z=4 Z Mtheta``;
        ``Mztheta=16 Z^2 Rref*x-(5/12)Rref Uref^2*x^(6/5)``;
        ``Mp=(5/2)Uref^2*x^(1/5)``.
        """

        one = mp.mpf(1)
        x = mp.exp(y)
        uref = mp.exp(self._log_pstar) / (one + z * z)
        rref = self._rref
        mtheta = (
            mp.mpf(5) / 8
            * mp.sqrt(2)
            * rref ** (mp.mpf(3) / 2)
            * uref
            * x ** (mp.mpf(8) / 5)
        )
        mz = 4 * z * rref * x
        return {
            "theta": mtheta,
            "z": mz,
            "theta_z": 4 * z * mtheta,
            "z_theta": 16 * z * z * rref * x
            - mp.mpf(5) / 12 * rref * uref * uref * x ** (mp.mpf(6) / 5),
            "p": mp.mpf(5) / 2 * uref * uref * x ** (mp.mpf(1) / 5),
        }

    def _profile_at_y(self, y: mp.mpf, z: mp.mpf) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
        """Return ``R``, ``Utheta``, ``Uz`` at a finite source coordinate."""

        log_r = self.profile.log_at(
            self.schedule.logRref,
            mp.nstr(y, self.precision),
        )
        values = self.profile.values(log_r, float(z))
        radius = self._rref * mp.exp(y)
        return radius, _mp(values["Utheta"]), _mp(values["Uz"])

    def _integrand_at_y(self, y: mp.mpf, z: mp.mpf) -> dict[str, mp.mpf]:
        radius, utheta, uz = self._profile_at_y(y, z)
        # dR = R dy.  The p integrand consequently has no R factor.
        theta_density = mp.sqrt(2 * radius) * utheta * radius
        return {
            "theta": theta_density,
            "z": uz * radius,
            "theta_z": theta_density * uz,
            "z_theta": (uz * uz - utheta * utheta / 2) * radius,
            "p": utheta * utheta / 2,
        }

    def _finite_interval(
        self,
        left: mp.mpf,
        right: mp.mpf,
        z: mp.mpf,
        *,
        order: int,
    ) -> dict[str, mp.mpf]:
        if right <= left:
            return {key: mp.mpf(0) for key in MOMENT_KEYS}
        nodes, weights = leggauss(order)
        half = (right - left) / 2
        center = (right + left) / 2
        total = {key: mp.mpf(0) for key in MOMENT_KEYS}
        for node, weight in zip(nodes, weights):
            y = center + half * mp.mpf(str(float(node)))
            row = self._integrand_at_y(y, z)
            weight_mp = mp.mpf(str(float(weight)))
            for key in MOMENT_KEYS:
                total[key] += weight_mp * row[key]
        return {key: half * total[key] for key in MOMENT_KEYS}

    def _finite_increment(
        self,
        target_y: mp.mpf,
        z: mp.mpf,
        *,
        order: int,
    ) -> dict[str, mp.mpf]:
        # The source changes formula at y=1, reaches the exact axial turnoff
        # at y=exp(Md), and ends the requested finite source interval at y_d.
        # Splitting at all three locations prevents quadrature from straddling
        # the flat axial cutoff.
        cuts = [mp.mpf(0)]
        if target_y > 1:
            cuts.append(mp.mpf(1))
        if target_y > self._axial_turnoff_y:
            cuts.append(self._axial_turnoff_y)
        cuts.append(target_y)
        total = {key: mp.mpf(0) for key in MOMENT_KEYS}
        for left, right in zip(cuts, cuts[1:]):
            increment = self._finite_interval(left, right, z, order=order)
            for key in MOMENT_KEYS:
                total[key] += increment[key]
        return total

    def moments_at_log_radius(
        self,
        logR: Any,
        Z: Any,
        *,
        quadrature_order: int | None = None,
    ) -> dict[str, mp.mpf]:
        """Return the five actual cumulative moments through ``logR``.

        Supported domain is the reference branch and the finite initial source
        interval ``logR_ref <= logR <= logR_d``.  Later source stages are
        intentionally rejected because they require their own tail integrals.
        """

        with mp.workdps(self.precision):
            y_decimal, y = self._target_y(logR)
            z = _z(Z)
            if y_decimal > self._y_d_decimal:
                raise ValueError(
                    "reference moments stop at the axial turnoff logR_d; "
                    "later outer stages are not implemented"
                )
            if y <= 0:
                return self._reference_moments(y, z)
            order = self.quadrature_order if quadrature_order is None else int(quadrature_order)
            if order < 8:
                raise ValueError("quadrature_order must be at least 8")
            total = self._reference_moments(mp.mpf(0), z)
            increment = self._finite_increment(y, z, order=order)
            for key in MOMENT_KEYS:
                total[key] += increment[key]
            return total

    # Short aliases used by stress adapters.
    moments = moments_at_log_radius
    moment_values = moments_at_log_radius

    def moments_Z(
        self,
        logR: Any,
        Z: Any,
        *,
        step: Any = "1e-3",
        quadrature_order: int | None = None,
    ) -> dict[str, mp.mpf]:
        """Fourth-order centered Z derivatives of the same five moments."""

        with mp.workdps(self.precision):
            z = _z(Z)
            h = _mp(step)
            if h <= 0 or abs(z) + 2 * h >= 1:
                raise ValueError("Z derivative stencil must remain inside |Z| < 1")
            samples = {
                -2: self.moments_at_log_radius(logR, z - 2 * h, quadrature_order=quadrature_order),
                -1: self.moments_at_log_radius(logR, z - h, quadrature_order=quadrature_order),
                1: self.moments_at_log_radius(logR, z + h, quadrature_order=quadrature_order),
                2: self.moments_at_log_radius(logR, z + 2 * h, quadrature_order=quadrature_order),
            }
            return {
                key: (
                    samples[-2][key]
                    - 8 * samples[-1][key]
                    + 8 * samples[1][key]
                    - samples[2][key]
                )
                / (12 * h)
                for key in MOMENT_KEYS
            }

    def moments_with_Z(
        self,
        logR: Any,
        Z: Any,
        *,
        step: Any = "1e-3",
        quadrature_order: int | None = None,
    ) -> dict[str, Any]:
        """Return raw moments and their explicit fourth-order Z stencil."""

        moments = self.moments_at_log_radius(logR, Z, quadrature_order=quadrature_order)
        return {
            "moments": moments,
            "moments_Z": self.moments_Z(
                logR, Z, step=step, quadrature_order=quadrature_order
            ),
            "Z_step": _mp(step),
            "Z_stencil": "centered five-point fourth-order",
        }

    def moment_receipt(
        self,
        logR: Any,
        Z: Any,
        *,
        quadrature_order: int | None = None,
        digits: int | None = None,
    ) -> dict[str, Any]:
        """JSON-ready signed-log receipt while retaining raw API separately."""

        digits = self.precision if digits is None else int(digits)
        values = self.moments_at_log_radius(
            logR, Z, quadrature_order=quadrature_order
        )
        return {
            key: signed_log(value, digits=digits) for key, value in values.items()
        }

    def z_derivative_receipt(
        self,
        logR: Any,
        Z: Any,
        *,
        coarse_step: Any = "1e-3",
        fine_step: Any = "5e-4",
        quadrature_order: int | None = None,
        digits: int | None = None,
    ) -> dict[str, Any]:
        """Signed-log moments_Z values and explicit stencil refinement."""

        digits = self.precision if digits is None else int(digits)
        with mp.workdps(self.precision):
            coarse = self.moments_Z(
                logR, Z, step=coarse_step, quadrature_order=quadrature_order
            )
            fine = self.moments_Z(
                logR, Z, step=fine_step, quadrature_order=quadrature_order
            )
            refinement = {
                key: {
                    "absolute": mp.nstr(abs(fine[key] - coarse[key]), digits),
                    "relative_to_fine": mp.nstr(
                        _relative_error(coarse[key], fine[key]), digits
                    ),
                }
                for key in MOMENT_KEYS
            }
            return {
                "coarse_step": mp.nstr(_mp(coarse_step), digits),
                "fine_step": mp.nstr(_mp(fine_step), digits),
                "coarse": {key: signed_log(value, digits=digits) for key, value in coarse.items()},
                "fine": {key: signed_log(value, digits=digits) for key, value in fine.items()},
                "refinement": refinement,
                "stencil": "centered five-point fourth-order",
            }

    def radial_derivative_check(
        self,
        y: Any,
        Z: Any,
        *,
        step: Any = "1e-4",
        quadrature_order: int | None = None,
        digits: int | None = None,
    ) -> dict[str, Any]:
        """Check ``dM/dR`` against each defining integrand independently."""

        digits = self.precision if digits is None else int(digits)
        with mp.workdps(self.precision):
            y0 = _mp(y)
            h = _mp(step)
            if h <= 0 or y0 - h <= 0 or y0 + h >= self._y_d:
                raise ValueError("radial derivative check requires h < y < y_d-h")
            z = _z(Z)
            log_rm = self.profile.log_at(
                self.schedule.logRref, mp.nstr(y0 - h, self.precision)
            )
            log_rp = self.profile.log_at(
                self.schedule.logRref, mp.nstr(y0 + h, self.precision)
            )
            lower = self.moments_at_log_radius(log_rm, z, quadrature_order=quadrature_order)
            upper = self.moments_at_log_radius(log_rp, z, quadrature_order=quadrature_order)
            denominator = self._rref * mp.exp(y0) * (mp.exp(h) - mp.exp(-h))
            finite = {key: (upper[key] - lower[key]) / denominator for key in MOMENT_KEYS}
            # The p integrand in y is Utheta^2/2, while dMp/dR is that divided by R.
            # Compute all five expected d/dR values explicitly for the receipt.
            radius, utheta, uz = self._profile_at_y(y0, z)
            expected_R = {
                "theta": mp.sqrt(2 * radius) * utheta,
                "z": uz,
                "theta_z": mp.sqrt(2 * radius) * utheta * uz,
                "z_theta": uz * uz - utheta * utheta / 2,
                "p": utheta * utheta / (2 * radius),
            }
            rows = {}
            for key in MOMENT_KEYS:
                rows[key] = {
                    "finite_difference": signed_log(finite[key], digits=digits),
                    "defining_integrand": signed_log(expected_R[key], digits=digits),
                    "absolute_error": mp.nstr(abs(finite[key] - expected_R[key]), digits),
                    "relative_error": mp.nstr(_relative_error(finite[key], expected_R[key]), digits),
                }
            return {
                "y": mp.nstr(y0, digits),
                "Z": mp.nstr(z, digits),
                "step": mp.nstr(h, digits),
                "rows": rows,
                "derivative_variable": "R (centered finite difference in log R)",
            }

    def metadata(self) -> dict[str, Any]:
        return {
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "source_equations": "(2.21)-(2.22), reference profile (4.3)",
            "api": "PaperReferenceMoments(profile).moments(logR,Z)",
            "moment_keys": list(MOMENT_KEYS),
            "reference_branch": "exact analytic integrals from R=0 to R=Rref*exp(y), y<=0",
            "finite_branch": "Gauss-Legendre in y=log(R/Rref), split at y=1 and y=exp(Md), 0<=y<=y_d",
            "axial_turnoff_log_offset": mp.nstr(self._axial_turnoff_y, self.precision),
            "finite_upper_log_offset": str(self._y_d_decimal),
            "arbitrary_exponent": True,
            "boundary_target_substitution": False,
            "profile_source": "CorrectedSourceProfile.values(logR,Z)",
            "pressure_integrability": {
                "reference_Utheta_exponent": "1/10",
                "reference_Mp_integrand_in_x": "x^(-4/5) dx",
                "integrable_at_axis": True,
                "reference_Mp_formula": "(5/2) Uref^2 x^(1/5)",
            },
            "scope": "Temporary reference and initial source through axial turnoff; no later outer closure claim",
        }


def build_reference_moment_callable(
    profile: CorrectedSourceProfile,
    *,
    precision: int | None = None,
    quadrature_order: int = 64,
) -> PaperReferenceMoments:
    """Construct the callable object used by the actual stress adapter."""

    return PaperReferenceMoments(
        profile, precision=precision, quadrature_order=quadrature_order
    )


def compute_reference_moments(
    profile: CorrectedSourceProfile,
    logR: Any,
    Z: Any,
    *,
    precision: int | None = None,
    quadrature_order: int = 64,
) -> dict[str, mp.mpf]:
    """One-shot convenience API returning raw mpmath moment values."""

    return PaperReferenceMoments(
        profile, precision=precision, quadrature_order=quadrature_order
    ).moments(logR, Z)


def _difference_receipt(
    coarse: dict[str, mp.mpf], fine: dict[str, mp.mpf], *, digits: int
) -> dict[str, Any]:
    with mp.workdps(max(80, int(digits))):
        return {
            key: {
                "absolute": mp.nstr(abs(fine[key] - coarse[key]), digits),
                "relative_to_fine": mp.nstr(_relative_error(coarse[key], fine[key]), digits),
            }
            for key in MOMENT_KEYS
        }


def run() -> dict[str, Any]:
    """Run focused reference, finite-stage, Z, and M_R checks."""

    profile = CorrectedSourceProfile(precision=160, order=64)
    engine = PaperReferenceMoments(profile, precision=160, quadrature_order=64)
    schedule = profile.schedule
    y_d = _mp(schedule.y_d)
    # Keep the sample coordinate as a decimal string until inside an mp
    # context; constructing mp.mpf at the process default precision would
    # otherwise bake a binary-sized rounding into a 160-digit receipt.
    z = "0.3"
    samples_y = [mp.mpf("-1"), mp.mpf("0"), mp.mpf("0.5"), mp.mpf("2"), y_d]
    samples = []
    for y in samples_y:
        log_r = profile.log_at(schedule.logRref, mp.nstr(y, engine.precision))
        values = engine.moments(log_r, z)
        samples.append(
            {
                "y": mp.nstr(y, engine.precision),
                "logR": str(log_r),
                "moments": {key: signed_log(value, digits=engine.precision) for key, value in values.items()},
            }
        )

    y_test = mp.mpf("2")
    log_test = profile.log_at(schedule.logRref, mp.nstr(y_test, engine.precision))
    coarse = engine.moments(log_test, z, quadrature_order=32)
    fine = engine.moments(log_test, z, quadrature_order=64)
    z_refinement = engine.z_derivative_receipt(
        log_test, z, coarse_step="1e-3", fine_step="5e-4", digits=engine.precision
    )
    radial_checks = [
        engine.radial_derivative_check(y, z, step="1e-4", digits=engine.precision)
        for y in (mp.mpf("0.5"), mp.mpf("2"), mp.mpf("10"))
    ]
    report = {
        "metadata": engine.metadata(),
        "schedule_inputs": {
            "logPstar": str(schedule.logPstar),
            "logRref": str(schedule.logRref),
            "delta": str(schedule.delta),
            "Md": str(schedule.Md),
            "y_d": str(schedule.y_d),
        },
        "reference_formula_receipt": {
            "Utheta": "Uref*x^(1/10)",
            "Uz": "4 Z",
            "Mtheta": "(5/8)*sqrt(2)*Rref^(3/2)*Uref*x^(8/5)",
            "Mz": "4*Z*Rref*x",
            "Mtheta_z": "4*Z*Mtheta",
            "Mztheta": "16*Z^2*Rref*x-(5/12)*Rref*Uref^2*x^(6/5)",
            "Mp": "(5/2)*Uref^2*x^(1/5)",
        },
        "samples": samples,
        "quadrature_refinement_at_y_2": _difference_receipt(coarse, fine, digits=engine.precision),
        "Z_derivative_refinement_at_y_2": z_refinement,
        "independent_M_R_checks": radial_checks,
        "claims": {
            "actual_profile_integrals": True,
            "boundary_targets_substituted": False,
            "reference_pressure_integrable": True,
            "later_outer_moments_closed": False,
            "stress_cone_or_PDE_certified": False,
        },
    }
    return report


if __name__ == "__main__":
    receipt = run()
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(receipt, indent=2, default=str) + "\n", encoding="utf-8")
    summary = {
        "sample_count": len(receipt["samples"]),
        "reference_pressure_integrable": receipt["claims"]["reference_pressure_integrable"],
        "radial_checks": len(receipt["independent_M_R_checks"]),
    }
    print(json.dumps(summary, indent=2))
