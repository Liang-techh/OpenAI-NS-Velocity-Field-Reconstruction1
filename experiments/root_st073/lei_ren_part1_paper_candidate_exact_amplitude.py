"""Exact fixed-data amplitude from the six complex poles of ``g``.

This module keeps the prescribed tiny amplitude scale.  It does not fit or
renormalize the amplitude.  Nominal roots are computed at high precision and
each root is enclosed by a directed Rouché disk before interval evaluation of
the logarithmic primitive.  The public ``exact_amplitude`` function is the
small interface used by the candidate core driver.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


HERE = Path(__file__).resolve().parent
PRECISION = 260
ROOT_PRECISION = 720
RADIUS_TEXT = "1e-240"
J_TEXT = "1e-14"
DELTA_TEXT = "1e-200"
SIGMA_DENOMINATOR_TEXT = "500"
LOGC_TEXT = "5e151"
Z0_INTERVAL_TEXT = ("-1e-14", "0")


def _lo(value: Any) -> mp.mpf:
    if hasattr(value, "_mpi_"):
        return mp.make_mpf(value._mpi_[0])
    return mp.mpf(value)


def _hi(value: Any) -> mp.mpf:
    if hasattr(value, "_mpi_"):
        return mp.make_mpf(value._mpi_[1])
    return mp.mpf(value)


def _contains_zero(value: Any) -> bool:
    return _lo(value) <= 0 <= _hi(value)


def _iv(ctx: MPIntervalContext, value: Any):
    if hasattr(value, "_mpi_"):
        return value
    if isinstance(value, dict) and "lower_exact_mpf_tuple" in value:
        return ctx.mpf(
            [
                mp.make_mpf(tuple(value["lower_exact_mpf_tuple"])),
                mp.make_mpf(tuple(value["upper_exact_mpf_tuple"])),
            ]
        )
    return ctx.mpf(str(value))


def _point_scalar(ctx: MPIntervalContext, value: mp.mpf):
    return ctx.mpf([mp.make_mpf(value._mpf_), mp.make_mpf(value._mpf_)])


def _point_complex(ctx: MPIntervalContext, value: mp.mpc):
    return ctx.mpc(_point_scalar(ctx, value.real), _point_scalar(ctx, value.imag))


def _disk_complex(
    ctx: MPIntervalContext, value: mp.mpc, radius: mp.mpf
):
    # Do the endpoint subtraction at precision above the nominal-root
    # precision.  This keeps the requested radius from being narrowed when
    # a caller enters with mpmath's default 15-digit global context.
    with mp.workdps(max(ROOT_PRECISION + 40, int(ctx.dps) + 40)):
        center = _point_complex(ctx, value)
        radius_mp = mp.mpf(radius)
        real_lower = _lo(center.real) - radius_mp
        real_upper = _hi(center.real) + radius_mp
        imag_lower = _lo(center.imag) - radius_mp
        imag_upper = _hi(center.imag) + radius_mp
    return ctx.mpc(
        ctx.mpf([real_lower, real_upper]),
        ctx.mpf([imag_lower, imag_upper]),
    )


def _encode(value: Any):
    """JSON encoding with exact interval endpoint tuples, including ivmpc."""

    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if hasattr(value, "_mpi_"):
        lower = mp.make_mpf(value._mpi_[0])
        upper = mp.make_mpf(value._mpi_[1])
        return {
            "lower": mp.nstr(lower, 90),
            "upper": mp.nstr(upper, 90),
            "width": mp.nstr(upper - lower, 60),
            "lower_exact_mpf_tuple": list(value._mpi_[0]),
            "upper_exact_mpf_tuple": list(value._mpi_[1]),
        }
    if isinstance(value, mp.mpf):
        return {
            "value": mp.nstr(value, 90),
            "log_abs": None if value == 0 else mp.nstr(mp.log(abs(value)), 70),
            "exact_mpf_tuple": list(value._mpf_),
        }
    if isinstance(value, mp.mpc):
        return {"real": _encode(value.real), "imag": _encode(value.imag)}
    if hasattr(value, "real") and hasattr(value, "imag"):
        return {"real": _encode(value.real), "imag": _encode(value.imag)}
    if isinstance(value, dict):
        return {key: _encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_encode(item) for item in value]
    return value


def _parameters(ctx: MPIntervalContext) -> dict[str, Any]:
    j = ctx.mpf(J_TEXT)
    delta = ctx.mpf(DELTA_TEXT)
    sigma = j / ctx.mpf(SIGMA_DENOMINATOR_TEXT)
    a = (9 - delta) / 2
    return {
        "j": j,
        "delta": delta,
        "sigma": sigma,
        "a": a,
        "logC": ctx.mpf(LOGC_TEXT),
        "radius": ctx.mpf(RADIUS_TEXT),
    }


def _H(ctx: MPIntervalContext, z, params: dict[str, Any], target=None):
    if target is None:
        target = ctx.mpc(0, 0)
    return (
        -4 * z**3
        - params["j"] * z**2
        + params["a"] * z
        + params["j"]
        - target
    )


def _Hprime(ctx: MPIntervalContext, z, params: dict[str, Any]):
    return -12 * z**2 - 2 * params["j"] * z + params["a"]


def _L(ctx: MPIntervalContext, z, params: dict[str, Any]):
    return 1 - params["delta"] * z**2


def _nominal_roots(params: dict[str, Any]):
    """Return high-precision nominal roots and the real H root."""

    with mp.workdps(ROOT_PRECISION):
        j = mp.mpf(J_TEXT)
        delta = mp.mpf(DELTA_TEXT)
        sigma = j / mp.mpf(SIGMA_DENOMINATOR_TEXT)
        a = (9 - delta) / 2

        def roots(target):
            return mp.polyroots(
                [-4, -j, a, j - target],
                maxsteps=2000,
                error=False,
            )

        zero_roots = roots(mp.mpc(0, 0))
        z0 = min(zero_roots, key=lambda value: abs(value))
        poles = []
        for sign in (1, -1):
            target = mp.mpc(0, sign * sigma)
            for index, value in enumerate(roots(target)):
                poles.append(
                    {
                        "target_sign": sign,
                        "nominal_index": index,
                        "nominal": mp.mpc(value.real, value.imag),
                    }
                )
        return {"z0": mp.mpc(z0.real, z0.imag), "poles": poles}


def _rouche_disk(
    ctx: MPIntervalContext,
    center: mp.mpc,
    target,
    params: dict[str, Any],
):
    radius = params["radius"]
    center_point = _point_complex(ctx, center)
    disk = _disk_complex(ctx, center, mp.mpf(RADIUS_TEXT))
    residual = _H(ctx, center_point, params, target)
    derivative = _Hprime(ctx, disk, params)
    lhs = abs(residual) + 12 * abs(disk) * radius**2 + 4 * radius**3
    rhs = abs(derivative) * radius
    margin = _lo(rhs) - _hi(lhs)
    return {
        "center": center_point,
        "disk": disk,
        "residual_abs_upper": _hi(abs(residual)),
        "derivative_abs_lower": _lo(abs(derivative)),
        "rouche_lhs_upper": _hi(lhs),
        "rouche_rhs_lower": _lo(rhs),
        "rouche_margin_lower": margin,
        "strict_rouche": margin > 0,
    }


def _root_data(ctx: MPIntervalContext):
    params = _parameters(ctx)
    nominal = _nominal_roots(params)
    sigma = params["sigma"]
    zero_target = ctx.mpc(0, 0)
    z0_rouche = _rouche_disk(ctx, nominal["z0"], zero_target, params)
    pole_records = []
    for row in nominal["poles"]:
        target = ctx.mpc(0, row["target_sign"] * sigma)
        proof = _rouche_disk(ctx, row["nominal"], target, params)
        pole_records.append(
            {
                **row,
                **proof,
                "target": target,
                "imaginary_sign_certified": (
                    _lo(proof["disk"].imag) > 0
                    or _hi(proof["disk"].imag) < 0
                ),
            }
        )

    all_centers = [nominal["z0"]] + [row["nominal"] for row in pole_records]
    separation_margins = []
    for i, left in enumerate(all_centers):
        for right in all_centers[i + 1 :]:
            left_iv = _point_complex(ctx, left)
            right_iv = _point_complex(ctx, right)
            separation_margins.append(
                _lo(abs(left_iv - right_iv)) - 2 * _hi(params["radius"])
            )
    min_separation = min(separation_margins)

    # H is strictly increasing on [-j, 0], and its endpoint signs bracket
    # the small root.  Conjugation plus uniqueness in the Rouché disk makes
    # the enclosed H root real.
    j = params["j"]
    interval_left = -j
    hprime_small_lower = params["a"] - 14 * j**2
    h_left = _H(ctx, ctx.mpc(interval_left, 0), params)
    h_right = _H(ctx, ctx.mpc(0, 0), params)
    z0_disk = z0_rouche["disk"]
    z0_real = ctx.mpf(
        [_lo(z0_disk.real), _hi(z0_disk.real)]
    )
    z0_real_root_checks = {
        "disk_in_minus_j_zero": _lo(z0_real) >= _lo(-j)
        and _hi(z0_real) <= _hi(ctx.mpf(0)),
        "Hprime_positive_on_interval": _lo(hprime_small_lower) > 0,
        "H_left_negative": _hi(h_left.real) < 0,
        "H_right_positive": _lo(h_right.real) > 0,
        "root_real_by_conjugate_uniqueness": z0_rouche["strict_rouche"],
        "nominal_z0_real": nominal["z0"].imag == 0,
    }

    return {
        "params": params,
        "z0": z0_rouche,
        "z0_real": z0_real,
        "poles": pole_records,
        "minimum_pairwise_separation_margin": min_separation,
        "root_separation_passed": min_separation > 0,
        "all_rouche_disks_strict": z0_rouche["strict_rouche"]
        and all(row["strict_rouche"] for row in pole_records),
        "z0_interval_checks": z0_real_root_checks,
        "z0_unique_in_minus_j_zero": all(z0_real_root_checks.values()),
    }


def _conjugate_pairing(ctx: MPIntervalContext, roots: list[dict[str, Any]]):
    pairs = []
    unused = set(range(len(roots)))
    with mp.workdps(ROOT_PRECISION):
        tolerance = mp.mpf("1e-600")
        while unused:
            index = min(unused)
            unused.remove(index)
            source = roots[index]
            target = min(
                unused,
                key=lambda candidate: abs(
                    roots[candidate]["nominal"]
                    - source["nominal"].conjugate()
                ),
            )
            distance = abs(
                roots[target]["nominal"]
                - source["nominal"].conjugate()
            )
            if distance > tolerance:
                raise AssertionError("complex pole conjugate pairing failed")
            unused.remove(target)
            pairs.append((index, target, distance))
    return pairs


def _real_z(ctx: MPIntervalContext, value: Any):
    if hasattr(value, "_mpi_"):
        return value
    return ctx.mpf(str(value))


def _log_ratio(
    ctx: MPIntervalContext,
    z,
    z0,
    pole: dict[str, Any],
):
    ratio = (ctx.mpc(z, 0) - pole["disk"]) / (
        ctx.mpc(z0, 0) - pole["disk"]
    )
    real_positive = _lo(ratio.real) > 0
    imag_contains_zero = _contains_zero(ratio.imag)
    branch_cut_rectangle_possible = (
        _lo(ratio.real) <= 0 and imag_contains_zero
    )
    continuous_branch = pole["imaginary_sign_certified"] and _lo(
        abs(ctx.mpc(z0, 0) - pole["disk"])
    ) > 0
    if branch_cut_rectangle_possible and not continuous_branch:
        raise ValueError("principal logarithm rectangle meets an uncertified branch cut")
    return {
        "ratio": ratio,
        "log": ctx.log(ratio),
        "ratio_real_positive": real_positive,
        "ratio_imag_contains_zero": imag_contains_zero,
        "principal_rectangle_safe": not branch_cut_rectangle_possible,
        "continuous_real_path_branch": continuous_branch,
        "path_argument": "Im((t-r)/(Z0-r))=Im(r)*(t-Z0)/|Z0-r|^2; nonzero off the anchor",
    }


def _conjugate(value):
    return value.real - value.imag * 1j


def _g_direct(ctx: MPIntervalContext, z, params):
    h = _H(ctx, z, params)
    return _L(ctx, z, params) * h / (h**2 + params["sigma"]**2)


def exact_amplitude(
    ctx: MPIntervalContext | None = None,
    Lambda: Any = "1e120",
    Z: Any = ".3",
) -> dict[str, Any]:
    """Return directed fixed-data amplitude intervals at one real Z.

    The returned ``F0_interval`` is a real positive interval.  ``G_interval``
    and ``logF0_interval`` are real intervals; complex diagnostic intervals
    are included separately.  ``Lambda`` and ``Z`` may be decimal strings or
    scalar/interval values accepted by ``MPIntervalContext``.
    """

    if ctx is None:
        ctx = MPIntervalContext()
        ctx.dps = PRECISION
    params = _parameters(ctx)
    roots = _root_data(ctx)
    z = _real_z(ctx, Z)
    z0 = roots["z0_real"]
    pole_records = roots["poles"]
    pairs = _conjugate_pairing(ctx, pole_records)

    log_terms = {}
    G_complex = ctx.mpc(0, 0)
    Gprime_from_poles = ctx.mpc(0, 0)
    for left_index, right_index, _ in pairs:
        left = pole_records[left_index]
        right = pole_records[right_index]
        left_residue = _L(ctx, left["disk"], params) / (
            2 * _Hprime(ctx, left["disk"], params)
        )
        left_log = _log_ratio(ctx, z, z0, left)
        right_log = _log_ratio(ctx, z, z0, right)
        left_term = left_residue * left_log["log"]
        right_residue = _L(ctx, right["disk"], params) / (
            2 * _Hprime(ctx, right["disk"], params)
        )
        right_term = right_residue * right_log["log"]
        G_complex += left_term + right_term
        Gprime_from_poles += left_residue / (ctx.mpc(z, 0) - left["disk"])
        Gprime_from_poles += right_residue / (ctx.mpc(z, 0) - right["disk"])
        log_terms[str(left_index)] = left_log
        log_terms[str(right_index)] = right_log

    G_interval = G_complex.real
    G_imag = G_complex.imag
    Lambda_iv = _iv(ctx, Lambda)
    logF0 = -params["logC"] - Lambda_iv * G_interval
    F0 = ctx.exp(logF0)
    g_value = _g_direct(ctx, ctx.mpc(z, 0), params)
    derivative_difference = Gprime_from_poles - g_value
    anchor_ratio = _log_ratio(ctx, z0, z0, pole_records[0])["log"]
    all_principal_safe = all(
        item["principal_rectangle_safe"]
        or item["continuous_real_path_branch"]
        for item in log_terms.values()
    )
    all_continuous = all(
        item["continuous_real_path_branch"] for item in log_terms.values()
    )
    denominator_lower = params["sigma"] ** 2
    checks = {
        "all_rouche_disks_strict": roots["all_rouche_disks_strict"],
        "root_separation_passed": roots["root_separation_passed"],
        "z0_unique_in_minus_j_zero": roots["z0_unique_in_minus_j_zero"],
        "all_pole_imaginary_signs_certified": all(
            row["imaginary_sign_certified"] for row in pole_records
        ),
        "log_path_principal_or_continuous": all_principal_safe,
        "G_imaginary_part_contains_zero": _contains_zero(G_imag),
        "Gprime_minus_g_contains_zero": _contains_zero(
            derivative_difference.real
        )
        and _contains_zero(derivative_difference.imag),
        "anchor_contains_zero": _contains_zero(anchor_ratio.real)
        and _contains_zero(anchor_ratio.imag),
        "real_denominator_positive": _lo(denominator_lower) > 0,
        "F0_positive_interval": _lo(F0) > 0,
    }
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise ValueError("exact amplitude validation failed: " + ", ".join(failed))
    return {
        "Lambda": Lambda_iv,
        "Z": z,
        "Z0_interval": z0,
        "G_interval": G_interval,
        "G_complex_interval": G_complex,
        "G_imag_interval": G_imag,
        "G_imaginary_part_contains_zero": _contains_zero(G_imag),
        "logF0_interval": logF0,
        "F0_interval": F0,
        "F0_positive_interval": _lo(F0) > 0,
        "log_terms": log_terms,
        "log_path_principal_or_continuous": all_principal_safe,
        "log_path_continuous_branch_checks": all_continuous,
        "g_direct_interval": g_value,
        "Gprime_from_poles_interval": Gprime_from_poles,
        "Gprime_minus_g_contains_zero": _contains_zero(
            derivative_difference.real
        )
        and _contains_zero(derivative_difference.imag),
        "Gprime_minus_g_interval": derivative_difference,
        "anchor_log_interval": anchor_ratio,
        "anchor_contains_zero": _contains_zero(anchor_ratio.real)
        and _contains_zero(anchor_ratio.imag),
        "real_denominator_lower_bound": denominator_lower,
        "real_denominator_positive": _lo(denominator_lower) > 0,
        "checks": checks,
        "root_data": roots,
        "conjugate_pairs": pairs,
    }


def _independent_integral_check(g_interval, ctx: MPIntervalContext):
    """Compare the pole primitive with an independent high-precision integral."""

    with mp.workdps(340):
        j = mp.mpf(J_TEXT)
        delta = mp.mpf(DELTA_TEXT)
        sigma = j / mp.mpf(SIGMA_DENOMINATOR_TEXT)
        a = (9 - delta) / 2
        z0 = min(
            mp.polyroots([-4, -j, a, j], maxsteps=2000, error=False),
            key=lambda value: abs(value),
        ).real
        z = mp.mpf(".3")

        def g(value):
            h = -4 * value**3 - j * value**2 + a * value + j
            return (1 - delta * value**2) * h / (h**2 + sigma**2)

        integral = mp.quad(g, [z0, z0 + mp.mpf("1e-14"), z])
        contained = _lo(g_interval) <= integral <= _hi(g_interval)
        return {
            "integral_value": integral,
            "contained_in_G_interval": contained,
            "split_points": [z0, z0 + mp.mpf("1e-14"), z],
            "non_singular_real_denominator_bound": (
                ctx.mpf(J_TEXT) / ctx.mpf(SIGMA_DENOMINATOR_TEXT)
            )
            ** 2,
        }


def run() -> dict[str, Any]:
    ctx = MPIntervalContext()
    ctx.dps = PRECISION
    with mp.workdps(ROOT_PRECISION + 40):
        roots = _root_data(ctx)
        values = {
            "Lambda48": exact_amplitude(ctx, Lambda="1e48", Z=".3"),
            "Lambda120": exact_amplitude(ctx, Lambda="1e120", Z=".3"),
        }
        integral_check = _independent_integral_check(
            values["Lambda120"]["G_interval"], ctx
        )
        if not integral_check["contained_in_G_interval"]:
            raise ValueError("independent g integral is outside the G interval")
        report = {
            "precision": PRECISION,
            "root_nominal_precision": ROOT_PRECISION,
            "rouche_radius": ctx.mpf(RADIUS_TEXT),
            "parameters": _parameters(ctx),
            "root_data": roots,
            "Lambda_values_at_Z_.3": values,
            "independent_g_integral_check": integral_check,
            "call_signature": "exact_amplitude(ctx, Lambda='1e120', Z='.3')",
            "analytic_primitive": "G(Z)=sum_r L(r)/(2 H'(r))*log((Z-r)/(Z0-r))",
            "F0_definition": "F0=exp(-logC-Lambda*G)",
            "nominal_fit_or_amplitude_reset": False,
            "tiny_amplitude_preserved": True,
            "source_parameter_errors_enclosed": False,
            "shared_core_regenerated": False,
            "global_matching_certified": False,
            "stress_cone_certified": False,
        }
        path = Path(__file__).with_suffix(".json")
        path.write_text(json.dumps(_encode(report), indent=2) + "\n", encoding="utf-8")
        print(
            "root disks strict",
            roots["all_rouche_disks_strict"],
            "separation",
            roots["root_separation_passed"],
            "G(.3)",
            mp.nstr(_hi(values["Lambda120"]["G_interval"]), 24),
            "F0 Lambda120 positive",
            values["Lambda120"]["F0_positive_interval"],
            flush=True,
        )
        return report


if __name__ == "__main__":
    run()
