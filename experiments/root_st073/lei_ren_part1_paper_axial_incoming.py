"""Incoming axial moment data for Lei--Ren Part I, Section 7.5.

This receipt evaluates the temporary reference-plus-outer candidate up to
``R_p``.  It supplies the signed logarithms of the two linear inputs in
(7.31),

    m_1 = M^z(R_p)/(R_p E_p),
    m_2 = M^{theta z}(R_p)/(R_p sqrt(2 R_p) E_p^2),

and the prior axial-energy quantity

    E_prior = mu M^{z theta}(R_p)/(R_p E_p^2),

used on the right of (7.34).  The source moment definitions are (2.21)--
(2.22), lines 2159--2216 of the cached TeX; the normalized inputs and the
energy equation are (7.30)--(7.34), lines 10793--11075.

The long ``R^{-1/2-mu}`` stage is integrated analytically in logarithmic
radius.  Finite transition stages use Gauss--Legendre quadrature, and all
large exponential ratios are retained as arbitrary-precision signed logs.
This is a temporary source candidate: no regular axis core, radial velocity,
axial pulse, pressure closure, cone, PDE, or global admissibility claim is
made here.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
import json
import math
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

from lei_ren_part1_paper_outer import (  # noqa: E402
    PaperOuterSchedule,
    SOURCE,
    SOURCE_SECTION,
    SOURCE_VERSION,
)


DEFAULT_INPUTS = {
    "logPstar": "14",
    "logRref": "10",
    "delta": "1e-32",
    "Md": ".5",
    "c_mu": ".001",
    "c_delta": ".001",
    "c_epsilon": ".01",
}


def _mp_decimal(value: Any) -> mp.mpf:
    return mp.mpf(str(value))


def _row(schedule: PaperOuterSchedule, y: Any, z: float) -> dict[str, Any]:
    """Evaluate at a finite stage-local y without materializing R."""

    if isinstance(y, Decimal):
        y_decimal = y
    else:
        y_decimal = Decimal(str(y))
    return schedule.at_log_radius(schedule.logRref + y_decimal, z)


def _integrate_finite(
    schedule: PaperOuterSchedule,
    z: float,
    left: Any,
    right: Any,
    function: Callable[[dict[str, Any], mp.mpf], mp.mpf],
    *,
    order: int,
) -> mp.mpf:
    """Gauss--Legendre integral over a finite log-radius stage."""

    a = mp.mpf(str(left))
    b = mp.mpf(str(right))
    if b <= a:
        return mp.mpf("0")
    nodes, weights = leggauss(order)
    half = (b - a) / 2
    center = (a + b) / 2
    total = mp.mpf("0")
    for node, weight in zip(nodes, weights):
        y = center + half * mp.mpf(str(float(node)))
        row = _row(schedule, y, z)
        total += mp.mpf(str(float(weight))) * function(row, y)
    return half * total


def _signed_log(value: mp.mpf, *, digits: int) -> dict[str, Any]:
    """Represent a signed arbitrary-exponent value without float underflow."""

    if value == 0:
        return {"sign": 0, "log_abs": None, "arbitrary_exponent_value": "0"}
    return {
        "sign": 1 if value > 0 else -1,
        "log_abs": mp.nstr(mp.log(abs(value)), digits),
        "arbitrary_exponent_value": mp.nstr(value, digits),
    }


def _signed_log_scaled(
    value: mp.mpf,
    log_scale: mp.mpf,
    *,
    digits: int,
    negate: bool = False,
) -> dict[str, Any]:
    """Scale a signed value in log space without materializing tiny factors."""

    if value == 0:
        return {"sign": 0, "log_abs": None, "arbitrary_exponent_value": "0"}
    sign = 1 if value > 0 else -1
    if negate:
        sign = -sign
    return {
        "sign": sign,
        "log_abs": mp.nstr(mp.log(abs(value)) + log_scale, digits),
        "arbitrary_exponent_value": mp.nstr(
            sign * mp.exp(mp.log(abs(value)) + log_scale), digits
        ),
    }


def _stage_integrals(
    schedule: PaperOuterSchedule,
    z: float,
    *,
    order: int,
) -> dict[str, mp.mpf]:
    """Compute dimensionless incoming moment integrals.

    The returned quantities have the common ``R_ref`` factors removed:
    ``I_z=M^z/R_ref``, ``I_theta_z=M^{theta z}/(sqrt(2)R_ref^(3/2))``,
    ``I_uz2=integral(Uz^2 dR)/R_ref``, and
    ``I_swirl=integral(Utheta^2 dR/2)/R_ref``.
    """

    y_cut = _mp_decimal(schedule.Md.exp())
    y_d = _mp_decimal(schedule.y_d)
    y_w = _mp_decimal(schedule.y_w)
    y_p = _mp_decimal(schedule.y_p)
    z_mp = mp.mpf(str(z))

    def uz_integrand(row: dict[str, Any], y: mp.mpf) -> mp.mpf:
        return mp.exp(y) * mp.mpf(str(row["Uz"]))

    def theta_z_integrand(row: dict[str, Any], y: mp.mpf) -> mp.mpf:
        return (
            mp.exp(mp.mpf("1.5") * y)
            * _mp_decimal(row["Utheta"])
            * mp.mpf(str(row["Uz"]))
        )

    def uz2_integrand(row: dict[str, Any], y: mp.mpf) -> mp.mpf:
        uz = mp.mpf(str(row["Uz"]))
        return mp.exp(y) * uz * uz

    # On y<0, U^z=4Z and U^theta=Pstar*exp(y/10)/(1+Z^2).
    i_z = mp.mpf(4) * z_mp
    i_theta_z = (
        mp.mpf(4) * z_mp / (1 + z_mp * z_mp)
        * _mp_decimal(schedule.logPstar.exp())
        / mp.mpf("1.6")
    )
    i_uz2 = mp.mpf(16) * z_mp * z_mp
    pstar = _mp_decimal(schedule.logPstar.exp())
    i_swirl = pstar * pstar / (mp.mpf("2.4") * (1 + z_mp * z_mp) ** 2)

    # Uz vanishes at y >= exp(Md), so its incoming integrals have only one
    # finite outer transition stage.
    # The source cutoff changes formula at y=1.  Splitting there keeps the
    # finite-stage refinement meaningful instead of asking one high-order
    # rule to straddle the flat transition.
    axial_stages = (
        (Decimal(0), min(Decimal(1), schedule.Md.exp())),
        (Decimal(1), y_cut),
    )
    for left, right in axial_stages:
        i_z += _integrate_finite(schedule, z, left, right, uz_integrand, order=order)
        i_theta_z += _integrate_finite(
            schedule, z, left, right, theta_z_integrand, order=order
        )
        i_uz2 += _integrate_finite(schedule, z, left, right, uz2_integrand, order=order)

    def swirl_integrand(row: dict[str, Any], y: mp.mpf) -> mp.mpf:
        u = _mp_decimal(row["Utheta"])
        return mp.exp(y) * u * u / 2

    # Resolve finite slope transitions separately, then integrate the exact
    # power stage [y_w,y_p] in closed form.  This avoids e^(13/mu).
    for left, right in (
        (Decimal(0), Decimal(1)),
        (Decimal(1), schedule.y_d),
        (schedule.y_d, schedule.y_w),
    ):
        i_swirl += _integrate_finite(
            schedule, z, left, right, swirl_integrand, order=order
        )
    row_w = _row(schedule, schedule.y_w, z)
    u_w = _mp_decimal(row_w["Utheta"])
    mu = _mp_decimal(schedule.mu)
    length = y_p - y_w
    if length > 0:
        decay_integral = -mp.expm1(-2 * mu * length) / (2 * mu)
        i_swirl += mp.exp(y_w) * u_w * u_w / 2 * decay_integral

    return {
        "I_z": i_z,
        "I_theta_z": i_theta_z,
        "I_uz2": i_uz2,
        "I_swirl": i_swirl,
    }


def incoming_axial_moments(
    schedule: PaperOuterSchedule,
    Z: Any,
    *,
    order: int = 64,
    precision: int = 160,
) -> dict[str, Any]:
    """Return signed-log m1, m2, and E_prior at the source ``R_p``."""

    z = float(Z)
    if not math.isfinite(z) or abs(z) > 1:
        raise ValueError("Z must lie in [-1,1]")
    if order < 16:
        raise ValueError("order must be at least 16")
    if precision < 80:
        raise ValueError("precision must be at least 80 decimal digits")

    with mp.workdps(precision):
        integrals = _stage_integrals(schedule, z, order=order)
        y_p = _mp_decimal(schedule.y_p)
        row_p = _row(schedule, schedule.y_p, z)
        log_ep = _mp_decimal(row_p["log_angular_amplitude"])
        log_rp_ep = y_p + log_ep
        log_rp_sqrt2rp_ep2 = mp.mpf("1.5") * y_p + 2 * log_ep
        m1 = integrals["I_z"] / mp.exp(log_rp_ep)
        m2 = integrals["I_theta_z"] / mp.exp(log_rp_sqrt2rp_ep2)
        eprior = (
            _mp_decimal(schedule.mu)
            * (integrals["I_uz2"] - integrals["I_swirl"])
            / mp.exp(y_p + 2 * log_ep)
        )
        mu = _mp_decimal(schedule.mu)
        lambda1 = mp.mpf("0.5") - mu
        lambda2 = mp.mpf("0.5") - 2 * mu
        log_row_scale1 = -mp.mpf(13) * lambda1 / mu
        log_row_scale2 = -mp.mpf(13) * lambda2 / mu
        # Equation (7.31) has RHS -m_i minus the pulse term.  Keep both
        # the scaled base m_i and its signed RHS contribution explicit;
        # the corresponding pulse module reports its positive integral.
        base_m1 = m1
        base_m2 = m2
        return {
            "Z": z,
            "order": int(order),
            "precision": int(precision),
            "log_Ep": mp.nstr(log_ep, precision),
            "dimensionless_integrals": {
                name: mp.nstr(value, precision)
                for name, value in integrals.items()
            },
            "m1_Mz_over_RpEp": _signed_log(m1, digits=precision),
            "m2_Mtheta_z_over_Rp_sqrt2RpEp2": _signed_log(m2, digits=precision),
            "E_prior_mu_Mztheta_over_RpEp2": _signed_log(eprior, digits=precision),
            "row_normalization": {
                "lambda1": mp.nstr(lambda1, precision),
                "lambda2": mp.nstr(lambda2, precision),
                "log_scale1": mp.nstr(log_row_scale1, precision),
                "log_scale2": mp.nstr(log_row_scale2, precision),
                "scaled_base_m1": _signed_log_scaled(
                    base_m1, log_row_scale1, digits=precision
                ),
                "scaled_base_m2": _signed_log_scaled(
                    base_m2, log_row_scale2, digits=precision
                ),
                "scaled_base_rhs1_minus_m1": _signed_log_scaled(
                    base_m1, log_row_scale1, digits=precision, negate=True
                ),
                "scaled_base_rhs2_minus_m2": _signed_log_scaled(
                    base_m2, log_row_scale2, digits=precision, negate=True
                ),
                "definition": (
                    "row i of (7.31) is multiplied by exp(-13*lambda_i/mu); "
                    "these are the base -m_i contributions, before the pulse term"
                ),
            },
            "moment_sign_convention": {
                "Mz": "integral Uz dR",
                "Mtheta_z": "integral sqrt(2R) Utheta Uz dR",
                "Mztheta": "integral (Uz^2 - Utheta^2/2) dR",
            },
            "scope": (
                "Temporary source reference-plus-outer candidate through Rp; "
                "no regular axis core, axial pulse, pressure closure, cone, "
                "PDE, or global admissibility claim."
            ),
        }


def _difference(a: dict[str, Any], b: dict[str, Any], key: str, precision: int) -> dict[str, Any]:
    va = mp.mpf(a[key]["arbitrary_exponent_value"])
    vb = mp.mpf(b[key]["arbitrary_exponent_value"])
    diff = abs(vb - va)
    relative = diff / abs(vb) if vb else mp.mpf("0")
    return {
        "absolute": mp.nstr(diff, precision),
        "log_absolute": mp.nstr(mp.log(diff), precision) if diff else None,
        "relative_to_fine": mp.nstr(relative, precision) if vb else None,
        "log_relative_to_fine": mp.nstr(mp.log(relative), precision)
        if relative
        else None,
    }


def run() -> dict[str, Any]:
    schedule = PaperOuterSchedule(**DEFAULT_INPUTS)
    zs = (0.0, 0.5, -0.5)
    coarse = [incoming_axial_moments(schedule, z, order=32) for z in zs]
    fine = [incoming_axial_moments(schedule, z, order=64) for z in zs]
    keys = (
        "m1_Mz_over_RpEp",
        "m2_Mtheta_z_over_Rp_sqrt2RpEp2",
        "E_prior_mu_Mztheta_over_RpEp2",
    )
    refinement = [
        {
            "Z": z,
            "m1": _difference(c, f, keys[0], 80),
            "m2": _difference(c, f, keys[1], 80),
            "E_prior": _difference(c, f, keys[2], 80),
        }
        for z, c, f in zip(zs, coarse, fine)
    ]
    report = {
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_section": SOURCE_SECTION,
        "equations": {
            "moment_definitions": "(2.21)-(2.22), lines 2159-2216",
            "incoming_linear_inputs": "(7.30)-(7.33), lines 10793-10959",
            "prior_energy_input": "(7.34), lines 10968-11075",
        },
        "api": "incoming_axial_moments(schedule, Z, order, precision)",
        "inputs": DEFAULT_INPUTS,
        "schedule": schedule.metadata(),
        "R_p_log_offset": str(schedule.y_p),
        "coarse_order_32": coarse,
        "fine_order_64": fine,
        "quadrature_refinement": refinement,
        "parity_observation": (
            "m1 and m2 are odd in Z; E_prior is even. This is a numerical "
            "symmetry observation for the supplied candidate."
        ),
        "normalization": {
            "m1": "Mz/(Rp*Ep)",
            "m2": "Mtheta_z/(Rp*sqrt(2*Rp)*Ep^2)",
            "E_prior": "mu*Mztheta/(Rp*Ep^2)",
            "linear_row_rescaling": "multiply row i by exp(-13*lambda_i/mu) before source (7.31) solve",
            "base_rhs_logs": "Each returned row_normalization scaled_base_rhs field is -m_i*exp(-13*lambda_i/mu), retained as a signed log.",
        },
        "claims": {
            "incoming_data_only": True,
            "axial_pulse_solved": False,
            "axial_closure": False,
            "regular_axis": False,
            "pressure_closed": False,
            "cone_validated": False,
            "pde_validated": False,
        },
    }
    return report


if __name__ == "__main__":
    receipt = run()
    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(receipt, indent=2, default=str) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "Z": [row["Z"] for row in receipt["fine_order_64"]],
                "quadrature_refinement": receipt["quadrature_refinement"],
            },
            indent=2,
        )
    )
