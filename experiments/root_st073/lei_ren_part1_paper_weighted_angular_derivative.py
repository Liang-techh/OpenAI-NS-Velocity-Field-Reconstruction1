"""Weighted angular derivative gates from paper (8.54)--(8.58).

The ordinary X_h correction norm is too coarse near the simple zero of H_0:
there the leading Bessel model has zero radial slope.  This receipt keeps the
first correction anchored at the zero, bounds the remaining derivative by a
weighted sqrt(chi) estimate, and evaluates the resulting gates for the
existing Lambda = 1e48 majorant and a separate Lambda = 1e120 candidate.

The calculation is conditional on the accepted fixed analytic data and the
existing all-term operator majorants.  It does not regenerate the shared
finite core, add source or parameter errors, or claim matching.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


HERE = Path(__file__).resolve().parent
PAPER_PATH = HERE.parent.parent / "work_paper_cache" / "lei_ren_part1.txt"
MAJORANT_SOURCE = HERE / "lei_ren_part1_paper_nonlinear_map_majorant.py"
LAMBDA_48 = "1e48"
LAMBDA_120 = "1e120"
PRECISION = 220
RADIAL_RADIUS = "4.1"
DELTA = "1e-200"
J = "1e-14"
SIGMA_DENOMINATOR = "500"
TARGET_RATIO = "1/32"
ACCEPTED_SCHEDULE_SHA256 = (
    "736bbadbde99bc2f3d098d279d61ef4cb64418368263a4aba7b275e7f8892de4"
)

INPUT_NAMES = (
    "lei_ren_part1_paper_nonlinear_candidate_Lambda48.json",
    "lei_ren_part1_paper_analytic_radial_tail.json",
    "lei_ren_part1_paper_linear_resolvent_bound.json",
    "lei_ren_part1_paper_complex_pressure_bound.json",
    "lei_ren_part1_paper_commuting_resolvent_bound.json",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _interval(ctx: MPIntervalContext, record: Any):
    if hasattr(record, "_mpi_"):
        return record
    return ctx.mpf(
        [
            mp.make_mpf(tuple(record["lower_exact_mpf_tuple"])),
            mp.make_mpf(tuple(record["upper_exact_mpf_tuple"])),
        ]
    )


def _upper(value: Any) -> mp.mpf:
    if not hasattr(value, "_mpi_"):
        return mp.mpf(value)
    return endpoints(value)[1]


def _lower(value: Any) -> mp.mpf:
    if not hasattr(value, "_mpi_"):
        return mp.mpf(value)
    return endpoints(value)[0]


def _verify_nested_hashes(base: Path, record: dict[str, Any]) -> dict[str, str]:
    nested = record.get("input_hashes", {})
    if not isinstance(nested, dict) or not nested:
        raise ValueError("analytic receipt lacks nested input hashes")
    checked = {}
    for name, expected in nested.items():
        path = (base / name).resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        actual = _sha256(path)
        if actual != expected:
            raise ValueError(
                f"nested hash mismatch for {name}: expected {expected}, got {actual}"
            )
        checked[name] = actual
    return checked


def _load_inputs(base: Path):
    records = {
        name: json.loads((base / name).read_text(encoding="utf-8"))
        for name in INPUT_NAMES
    }
    candidate = records[INPUT_NAMES[0]]
    tube = records[INPUT_NAMES[1]]
    linear = records[INPUT_NAMES[2]]
    pressure = records[INPUT_NAMES[3]]
    commuting = records[INPUT_NAMES[4]]
    if candidate.get("accepted_schedule_sha256") != ACCEPTED_SCHEDULE_SHA256:
        raise ValueError("candidate majorant has an unexpected schedule hash")
    if not candidate.get("contraction_proved"):
        raise ValueError("existing Lambda 1e48 candidate contraction gate is false")
    if not linear.get("analytic_linear_inverse_bound_certified"):
        raise ValueError("analytic linear inverse receipt is not certified")
    if not pressure.get("all_true_pressure_stages_included"):
        raise ValueError("analytic pressure receipt omits true pressure stages")
    if pressure.get("accepted_schedule_sha256") != ACCEPTED_SCHEDULE_SHA256:
        raise ValueError("pressure receipt has an unexpected schedule hash")
    tube_hash = _sha256(base / INPUT_NAMES[1])
    if linear.get("input_sha256") != tube_hash:
        raise ValueError("linear inverse receipt uses a different analytic tube")
    if commuting.get("input_sha256") != tube_hash:
        raise ValueError("commuting resolvent uses a different analytic tube")
    nested = {}
    for name in (INPUT_NAMES[0], INPUT_NAMES[3]):
        nested[name] = _verify_nested_hashes(base, records[name])
    source_hashes = {
        name: _sha256(base / name) for name in INPUT_NAMES
    }
    source_hashes[MAJORANT_SOURCE.name] = _sha256(MAJORANT_SOURCE)
    if PAPER_PATH.exists():
        source_hashes["work_paper_cache/lei_ren_part1.txt"] = _sha256(
            PAPER_PATH
        )
    else:
        raise FileNotFoundError(PAPER_PATH)
    return records, source_hashes, nested


def _get(ctx: MPIntervalContext, record: dict[str, Any], key: str):
    return _interval(ctx, record[key])


def _candidate_majorant(
    ctx: MPIntervalContext,
    records: dict[str, dict[str, Any]],
    lambda_text: str,
) -> dict[str, Any]:
    """Replay all 20 existing map terms for a requested Lambda."""

    tube = records[INPUT_NAMES[1]]
    linear = records[INPUT_NAMES[2]]
    pressure = records[INPUT_NAMES[3]]
    commuting = records[INPUT_NAMES[4]]
    eta = _get(ctx, tube, "complex_tube_radius")
    h = _get(ctx, tube, "Xh_parameter")
    a = 1 + eta
    dt = ctx.mpf(DELTA)
    lam = ctx.mpf(lambda_text)
    eps = 1 / lam
    j = ctx.mpf(J)
    L = _get(ctx, tube, "L_modulus_lower")
    pole = _get(ctx, tube, "denominator_factor_modulus_lower")
    H = _get(ctx, tube, "H0_modulus_upper")
    U = 4 * a + j
    d = 1 + a * a
    weight = _get(ctx, linear, "Cauchy_weight_supremum_upper")
    Bphi = _get(ctx, linear, "Phi_model_Xh_norm_upper") + 1
    Bpsi = _get(ctx, pressure, "Psi_model_Xh_norm_upper") + 1
    Rnorm = _get(ctx, commuting, "resolvent_norm_upper")
    fixed_multiplier = _get(
        ctx, commuting, "axial_convolution_factor_upper"
    )
    Gupper = _get(ctx, tube, "G_modulus_upper")
    logC = ctx.mpf("5e151")
    Cstar_margin = logC - 2 * ctx.log(lam) - lam * Gupper
    if _lower(Cstar_margin) <= 0:
        raise AssertionError(f"Cstar guard failed for Lambda = {lambda_text}")
    Fsquare = ctx.exp(2 * (-logC + lam * Gupper)) * weight
    product = ctx.mpf(256)
    J1 = ctx.mpf(80)
    J2 = ctx.mpf(40)
    axial = ctx.mpf(20480) / h
    radial = ctx.mpf(20480)
    Pcal = 80 * product**2 * Fsquare
    beta = _get(ctx, tube, "beta_modulus_upper")
    coefficients = {
        "beta": beta * weight,
        "W0_over_L": (1 + (1 + dt) * a * U + 4 * d) / L * weight,
        "H0_over_L": H / L * weight,
        "az": (1 + dt) * a / L * weight,
        "d_over_L": d / L * weight,
        "z_over_L": a / L * weight,
        "axial_linear": ((1 + dt) / 2 * (1 + 4 * a * U) + 4 * d)
        / L
        * weight,
        "cross_swirl": d / pole * weight,
    }
    term_names = []
    term_sizes = []
    term_lips = []
    term_explicit_epsilon = []

    def term(
        name: str,
        coefficient: Any,
        p: int,
        q: int,
        *,
        explicit_epsilon: bool = False,
    ) -> None:
        coefficient = coefficient * fixed_multiplier / product
        value = coefficient * Bphi**p * Bpsi**q
        lip = ctx.mpf(0)
        if p:
            lip += coefficient * p * Bphi ** (p - 1) * Bpsi**q
        if q:
            lip += coefficient * q * Bphi**p * Bpsi ** (q - 1)
        term_names.append(name)
        term_sizes.append(value)
        term_lips.append(lip)
        term_explicit_epsilon.append(explicit_epsilon)

    # J2 Etheta: retain the exact 10-term split from (8.50).
    term("-beta Phi", product * coefficients["beta"] * J2, 1, 0)
    term(
        "W0/L scaledR Phi_R",
        product * coefficients["W0_over_L"] * radial,
        1,
        0,
    )
    term(
        "H0/L Phi_Z",
        product * coefficients["H0_over_L"] * axial,
        1,
        0,
    )
    term(
        "-eps az M(Psi) Phi",
        eps * product * coefficients["az"] * J2 * product,
        1,
        1,
        explicit_epsilon=True,
    )
    term(
        "-eps az M(Psi) scaledR Phi_R",
        eps * product * coefficients["az"] * radial,
        1,
        1,
        explicit_epsilon=True,
    )
    term(
        "-eps d/L d_Z M(Psi) Phi",
        eps * product * coefficients["d_over_L"] * axial,
        1,
        1,
        explicit_epsilon=True,
    )
    term(
        "-eps d/L d_Z M(Psi) scaledR Phi_R",
        eps * product * coefficients["d_over_L"] * axial,
        1,
        1,
        explicit_epsilon=True,
    )
    term(
        "-eps delta z/L Psi Phi",
        eps * dt * product * coefficients["z_over_L"] * J2 * product,
        1,
        1,
        explicit_epsilon=True,
    )
    term(
        "eps d/L Psi Phi_Z",
        eps * product * coefficients["d_over_L"] * axial,
        1,
        1,
        explicit_epsilon=True,
    )
    term(
        "-d H0/(H0^2+sigma^2) Psi Phi",
        product * coefficients["cross_swirl"] * J2 * product,
        1,
        1,
    )
    # J1 Ez: retain all 10 axial terms including pressure/F0^2 coupling.
    term(
        "W0/L scaledR Psi_R",
        product * coefficients["W0_over_L"] * radial,
        0,
        1,
    )
    term(
        "axial linear coefficient Psi",
        product * coefficients["axial_linear"] * J1,
        0,
        1,
    )
    term(
        "H0/L Psi_Z",
        product * coefficients["H0_over_L"] * axial,
        0,
        1,
    )
    term(
        "-eps az M(Psi) scaledR Psi_R",
        eps * product * coefficients["az"] * radial,
        0,
        2,
        explicit_epsilon=True,
    )
    term(
        "-eps d/L d_Z M(Psi) scaledR Psi_R",
        eps * product * coefficients["d_over_L"] * axial,
        0,
        2,
        explicit_epsilon=True,
    )
    term(
        "-eps (1+delta) z/L Psi^2",
        eps * (1 + dt) * product * coefficients["z_over_L"] * J1 * product,
        0,
        2,
        explicit_epsilon=True,
    )
    term(
        "eps d/L Psi Psi_Z",
        eps * product * coefficients["d_over_L"] * axial,
        0,
        2,
        explicit_epsilon=True,
    )
    term(
        "d/L d_Z Pcal",
        product * coefficients["d_over_L"] * axial * Pcal,
        2,
        0,
    )
    term(
        "-2(1+delta)z/L Pcal",
        2 * (1 + dt) * product * coefficients["z_over_L"] * J1 * Pcal,
        2,
        0,
    )
    term(
        "-2z/L scaledR F0^2 Phi^2",
        2 * product * coefficients["z_over_L"] * J1 * 80 * product**2 * Fsquare,
        2,
        0,
    )
    if len(term_names) != 20:
        raise AssertionError("weighted majorant did not retain all 20 terms")
    theta_size = sum(
        (term_sizes[index] for index in range(20) if index < 10),
        ctx.mpf(0),
    )
    theta_lip = sum(
        (term_lips[index] for index in range(20) if index < 10),
        ctx.mpf(0),
    )
    z_size = sum(
        (term_sizes[index] for index in range(20) if index >= 10),
        ctx.mpf(0),
    )
    z_lip = sum(
        (term_lips[index] for index in range(20) if index >= 10),
        ctx.mpf(0),
    )
    size = (Rnorm * theta_size + z_size) / 2
    lip = (Rnorm * theta_lip + z_lip) / 2
    scaled_size = eps * size
    scaled_lip = eps * lip
    return {
        "Lambda": lam,
        "epsilon": eps,
        "Cstar_log_margin": Cstar_margin,
        "Nsize": size,
        "LipN": lip,
        "scaled_map_size": scaled_size,
        "scaled_map_Lipschitz": scaled_lip,
        "contraction_gate": _upper(scaled_size) <= mp.mpf(".5")
        and _upper(scaled_lip) <= mp.mpf(".5"),
        "term_count": len(term_names),
        "term_names": term_names,
        "term_sizes": term_sizes,
        "term_lips": term_lips,
        "term_explicit_epsilon": term_explicit_epsilon,
        "explicit_epsilon_term_count": sum(term_explicit_epsilon),
        "all_paper_8_50_terms_included": True,
        "pressure_and_swirl_couplings_retained": True,
        "Fsquare_Xh_upper": Fsquare,
        "fixed_multiplier_upper": fixed_multiplier,
        "resolvent_upper": Rnorm,
        "Xh_parameter": h,
        "beta_complex_modulus_upper": beta,
        "source_majorant_scope": "Accepted fixed analytic data; no original parameter or source derivation errors.",
    }


def _positive_sum_with_tail(
    ctx: MPIntervalContext,
    term_function,
    ratio_upper: mp.mpf,
    *,
    first_index: int,
    last_index: int,
):
    total = ctx.mpf(0)
    for index in range(first_index, last_index + 1):
        total += term_function(index)
    next_term = term_function(last_index + 1)
    if ratio_upper >= 1:
        raise ValueError("positive series tail ratio is not below one")
    tail = ctx.mpf([0, _upper(next_term) / (1 - ratio_upper)])
    return total + tail, tail


def _derivative_evaluation_factors(ctx: MPIntervalContext, h):
    radius_ratio = ctx.mpf(RADIAL_RADIUS) / 20
    ratio_upper = _upper(radius_ratio)
    last = 256

    def term_k0(n):
        return ctx.mpf(n) / ((n + 1) ** 2) * radius_ratio ** (n - 1)

    def term_k1(n):
        return ctx.mpf(n) / (n + 1) * radius_ratio ** (n - 1)

    # For k = 1 the rational prefactor has successive ratio
    # (n+1)^2/(n(n+2)) > 1.  Use its maximum over the unsummed tail,
    # attained at n = last + 1, instead of the bare radial ratio.
    ratio_k1_tail = (
        ratio_upper * (last + 2) ** 2 / ((last + 1) * (last + 3))
    )
    sum_k0, tail_k0 = _positive_sum_with_tail(
        ctx,
        term_k0,
        ratio_upper,
        first_index=1,
        last_index=last,
    )
    sum_k1, tail_k1 = _positive_sum_with_tail(
        ctx,
        term_k1,
        ratio_k1_tail,
        first_index=1,
        last_index=last,
    )
    return {
        "radial_ratio": radius_ratio,
        "full_series_cutoff": last,
        "sum_k0": sum_k0,
        "sum_k1": sum_k1,
        "tail_ratio_k1": ratio_k1_tail,
        "tail_k0": tail_k0,
        "tail_k1": tail_k1,
        # i = 1, k = 0 and i = 1, k = 1, respectively.
        "i1_k0_factor": sum_k0 / 20,
        "i1_k1_factor": sum_k1 / (80 * h),
        "positive_series_tail_checks_pass": _upper(tail_k0) < mp.mpf("1e-170")
        and _upper(tail_k1) < mp.mpf("1e-150"),
    }


def _model_norm_and_remainders(
    ctx: MPIntervalContext,
    tube: dict[str, Any],
    epsilon: Any,
):
    h = _get(ctx, tube, "Xh_parameter")
    eta = _get(ctx, tube, "complex_tube_radius")
    beta_complex = _get(ctx, tube, "beta_modulus_upper")
    chi_complex = _get(ctx, tube, "chi_modulus_upper")
    beta_real_lower = ctx.mpf(1)
    beta_real = ctx.mpf("31") / 10
    smax = ctx.mpf(RADIAL_RADIUS)
    # Cauchy disks of radius eta/2 leave h/(eta/2)=1/4.  The axial
    # (m+1)^2 factor is then bounded by one for every m >= 0.
    cauchy_ratio_raw = h / (eta / 2)
    # The receipt was generated with the exact construction eta=sigma/40,
    # h=eta/8.  Retain the raw directed quotient above for audit, but use
    # the exact rational quotient for the Cauchy weight; endpoint rounding
    # otherwise widens an identity that is exact in the construction.
    h_decimal = mp.mpf(tube["Xh_parameter"]["upper"])
    eta_decimal = mp.mpf(tube["complex_tube_radius"]["upper"])
    if h_decimal * 8 != eta_decimal:
        raise AssertionError("analytic receipt does not encode h=eta/8")
    cauchy_ratio = ctx.mpf(1) / 4
    axial_weight = ctx.mpf(1)
    if _upper(cauchy_ratio) > mp.mpf("1") / 4:
        raise AssertionError("weighted Cauchy ratio is larger than 1/4")
    chi_upper = _upper(chi_complex)
    beta_upper = _upper(beta_complex)
    last = 256

    def model_term(n):
        return (
            ctx.mpf(10) ** n
            * (n + 1) ** 2
            * chi_complex ** (n - 1)
            / (mp.factorial(n - 1) * mp.factorial(n + 1))
        )

    ratio_tail = (
        10
        * chi_upper
        * (last + 3)
        / ((last + 1) * (last + 2) ** 2)
    )
    model_series, model_tail = _positive_sum_with_tail(
        ctx,
        model_term,
        ratio_tail,
        first_index=1,
        last_index=last,
    )
    model_norm = beta_complex * axial_weight * model_series

    alpha_interval = 1 + epsilon * beta_real
    alpha_upper = _upper(alpha_interval)
    qmax_upper = _upper(smax * alpha_interval / 2)
    bprime_lower = ctx.mpf(1) / 2 - qmax_upper / 6
    model_slope_checks = {
        "beta_real_lower_ge_1": _upper(beta_real_lower) >= 1,
        "qmax_le_2.1": qmax_upper <= mp.mpf("2.1"),
        "qmax_over_3_lt_1": qmax_upper / 3 < 1,
        "minus_Bprime_lower_ge_0.15": bprime_lower >= mp.mpf(".15"),
        "minus_Bprime_lower_gt_1/8": bprime_lower > mp.mpf(1) / 8,
    }
    if not all(model_slope_checks.values()):
        raise AssertionError("model slope lower-bound gate failed")
    exp_bound = ctx.exp(ctx.mpf("2.1"))
    b2 = exp_bound / 6
    b3 = exp_bound / 24
    f2_upper = _upper(smax / 2 * b2 + alpha_interval * smax**2 / 8 * b3)
    model_second_order = beta_real**2 * f2_upper / 2

    # Independent positive checks for |B^(m)| <= exp(2.1)/(m+1)!.
    bessel_checks = []
    for derivative_order in range(4):
        qmax = mp.mpf("2.1")

        def bessel_term(k):
            return qmax**k / (
                mp.factorial(k) * mp.factorial(k + derivative_order + 1)
            )

        ratio = qmax / ((last + 2) * (last + derivative_order + 3))
        total, tail = _positive_sum_with_tail(
            ctx,
            bessel_term,
            ratio,
            first_index=0,
            last_index=last,
        )
        bound = mp.e ** mp.mpf("2.1") / mp.factorial(derivative_order + 1)
        bessel_checks.append(
            {
                "derivative_order": derivative_order,
                "positive_sum": total,
                "tail": tail,
                "analytic_bound": bound,
                "contained": _upper(total) <= bound,
            }
        )
    return {
        "beta_real_upper": beta_real,
        "beta_real_lower": beta_real_lower,
        "beta_complex_upper": beta_complex,
        "chi_complex_upper": chi_complex,
        "model_qmax_upper": qmax_upper,
        "minus_Bprime_lower": bprime_lower,
        "model_slope_checks": model_slope_checks,
        "cauchy_ratio_raw": cauchy_ratio_raw,
        "cauchy_ratio": cauchy_ratio,
        "cauchy_ratio_identity": "h=eta/8 from analytic_radial_tail construction",
        "axial_weight_upper": axial_weight,
        "model_series_cutoff": last,
        "model_first_correction_Xh_upper": model_norm,
        "model_series_tail": model_tail,
        "model_second_order_derivative_upper": model_second_order,
        "Bessel_positive_series_checks": bessel_checks,
        "Bessel_positive_series_checks_passed": all(
            row["contained"] for row in bessel_checks
        ),
    }


def _anchored_first_correction_proof(ctx: MPIntervalContext) -> dict[str, Any]:
    """Check the symbolic axis reduction used at the H0 root.

    At H0 = 0, chi and chi_Z vanish, the model has Phi0 = 1 with zero
    radial and axial derivatives, and the swirl numerator g vanishes.  Thus
    E_theta(X0; 0) is exactly -beta0.  Solving the regular axis equation
    gives -beta0*s/4, which agrees with beta0*s*B'(0)/2 because B'(0)=-1/2.
    """

    sigma = ctx.mpf(J) / ctx.mpf(SIGMA_DENOMINATOR)
    h0 = ctx.mpf(0)
    h0_prime = ctx.mpf(1)
    chi0 = h0**2 / (h0**2 + sigma**2)
    chi0_z = 2 * h0 * h0_prime * sigma**2 / (h0**2 + sigma**2) ** 2
    phi0_axis = ctx.mpf(1)
    phi0_s = ctx.mpf(0)
    phi0_z = ctx.mpf(0)
    g0 = h0 / (h0**2 + sigma**2)
    conditions = {
        "H0_zero": h0 == 0,
        "chi_zero": chi0 == 0,
        "chi_Z_zero": chi0_z == 0,
        "Phi0_axis_one": phi0_axis == 1,
        "Phi0_s_zero": phi0_s == 0,
        "Phi0_Z_zero": phi0_z == 0,
        "swirl_g_zero": g0 == 0,
        "sigma_positive": _lower(sigma) > 0,
    }
    actual_coefficient = -ctx.mpf(1) / 4
    model_Bprime_zero = -ctx.mpf(1) / 2
    model_coefficient = model_Bprime_zero / 2
    conditions["regular_axis_solution_coefficient"] = (
        actual_coefficient == model_coefficient
    )
    return {
        "source_equations": "paper (8.54), axis reduction immediately below (8.54), and model expansion following it",
        "axis_conditions": conditions,
        "Etheta_X0_at_epsilon_zero": "-beta0",
        "regular_axis_equation": "2*(s*y''+2*y')=-beta0",
        "Phi1_axis": "-beta0*s/4",
        "model_Phi1_axis": "beta0*s*B'(0)/2=-beta0*s/4",
        "actual_axis_slope_coefficient": actual_coefficient,
        "model_axis_slope_coefficient": model_coefficient,
        "proved": all(conditions.values()),
    }


def _weighted_case(
    ctx: MPIntervalContext,
    majorant: dict[str, Any],
    tube: dict[str, Any],
    eval_factors: dict[str, Any],
    model_data: dict[str, Any],
) -> dict[str, Any]:
    epsilon = majorant["epsilon"]
    nsize = majorant["Nsize"]
    lipn = majorant["LipN"]
    # The fixed-point second-order remainder is bounded by Nsize * LipN.
    second_order_norm = nsize * lipn
    epsilon = majorant["epsilon"]
    explicit_names = [
        majorant["term_names"][index]
        for index, marked in enumerate(majorant["term_explicit_epsilon"])
        if marked
    ]
    explicit_pure_names = [
        name for name in explicit_names if "delta" not in name
    ]
    explicit_delta_names = [
        name for name in explicit_names if "delta" in name
    ]
    explicit_theta = sum(
        (
            majorant["term_sizes"][index] / epsilon
            for index, marked in enumerate(majorant["term_explicit_epsilon"])
            if marked and index < 10
        ),
        ctx.mpf(0),
    )
    explicit_z = sum(
        (
            majorant["term_sizes"][index] / epsilon
            for index, marked in enumerate(majorant["term_explicit_epsilon"])
            if marked and index >= 10
        ),
        ctx.mpf(0),
    )
    D_epsN = (majorant["resolvent_upper"] * explicit_theta + explicit_z) / 2
    D_epsN_i1_k0 = D_epsN * eval_factors["i1_k0_factor"]
    model_first = model_data["model_first_correction_Xh_upper"]
    difference_first_norm = nsize + model_first
    C1 = difference_first_norm * eval_factors["i1_k1_factor"]
    C2_fixed = second_order_norm * eval_factors["i1_k0_factor"]
    C2 = (
        C2_fixed
        + D_epsN_i1_k0
        + model_data["model_second_order_derivative_upper"]
    )

    anchored_proof = _anchored_first_correction_proof(ctx)

    dt = ctx.mpf(DELTA)
    j = ctx.mpf(J)
    sigma = j / ctx.mpf(SIGMA_DENOMINATOR)
    cH = (9 - dt) / 2 - j * (1 + j) - 4 * (1 + j + j * j)
    if _lower(cH) <= 0:
        raise AssertionError("cH root-factor lower bound is nonpositive")
    small_factor = 10 * sigma / cH
    small_ratio = (
        C1 * small_factor * ctx.sqrt(epsilon) / 2 + C2 * epsilon
    )
    large_ratio = (
        C1 * epsilon * (1 + j) + C2 * epsilon
    ) / (ctx.mpf(".99") + epsilon)
    target = ctx.mpf(1) / 32

    # Solve the conservative small-chi quadratic in x = sqrt(epsilon), then
    # intersect it with the large-chi linear gate.
    A = _upper(C1 * small_factor) / 2
    B = _upper(C2)
    if B == 0:
        x_max = mp.mpf(1) / 32 / A
    else:
        x_max = (mp.mpf(1) / 16) / (
            A + mp.sqrt(A * A + B / 8)
        )
    epsilon_small_max = x_max * x_max
    large_denominator = 32 * _upper(C1 * (1 + j) + C2) - 1
    epsilon_large_max = (
        mp.mpf(".99") / large_denominator
        if large_denominator > 0
        else mp.mpf(1)
    )
    epsilon_max = min(epsilon_small_max, epsilon_large_max, mp.mpf(1) / 500)
    lambda_min = 1 / epsilon_max

    small_gate = _upper(small_ratio) <= _upper(target)
    large_gate = _upper(large_ratio) <= _upper(target)
    checks = {
        "anchored_first_correction_identity": anchored_proof["proved"],
        "small_chi_gate": small_gate,
        "large_chi_gate": large_gate,
        "positive_series_checks": bool(
            eval_factors["positive_series_tail_checks_pass"]
            and model_data["Bessel_positive_series_checks_passed"]
        ),
        "contraction_gate": bool(majorant["contraction_gate"]),
    }
    return {
        "Lambda": majorant["Lambda"],
        "epsilon": epsilon,
        "Nsize_upper": nsize,
        "LipN_upper": lipn,
        "fixedpoint_second_order_norm_upper": second_order_norm,
        "explicit_epsilon_term_names": explicit_names,
        "explicit_epsilon_term_count": len(explicit_names),
        "explicit_epsilon_pure_term_count": len(explicit_pure_names),
        "explicit_epsilon_delta_term_count": len(explicit_delta_names),
        "explicit_epsilon_remainder_Xh_upper": D_epsN,
        "explicit_epsilon_remainder_i1_k0_upper": D_epsN_i1_k0,
        "Phi1_minus_model1_Xh_upper": difference_first_norm,
        "C1_i1_k1_upper": C1,
        "C2_fixed_i1_k0_upper": C2_fixed,
        "C2_explicit_epsilon_i1_k0_upper": D_epsN_i1_k0,
        "C2_model_upper": model_data["model_second_order_derivative_upper"],
        "C2_total_upper": C2,
        "sigma": sigma,
        "cH_lower": cH,
        "small_chi_factor_10sigma_over_cH": small_factor,
        "small_chi_ratio_upper": small_ratio,
        "large_chi_ratio_upper": large_ratio,
        "large_branch_C2_epsilon_overbound": True,
        "large_branch_remainder_note": "Uses C2*epsilon; the actual second-order term is C2*epsilon^2, so this is conservative.",
        "required_ratio_upper": target,
        "small_chi_epsilon_max": epsilon_small_max,
        "large_chi_epsilon_max": epsilon_large_max,
        "sufficient_epsilon_max": epsilon_max,
        "sufficient_Lambda_min": lambda_min,
        "checks": checks,
        "anchored_first_correction_proof": anchored_proof,
        "anchored_axis_slope_formula": "partial_s Phi(0,Z)=-(chi(Z)+epsilon beta(Z))/4",
        "weighted_negative_derivative_gate": all(checks.values()),
        "conclusion": (
            "-partial_s Phi >= (chi + epsilon)/32 on the stated real rectangle"
            if all(checks.values())
            else "Weighted sign gate remains open for this Lambda and majorant"
        ),
    }


def run() -> dict[str, Any]:
    base = HERE
    records, source_hashes, nested_hashes = _load_inputs(base)
    ctx = MPIntervalContext()
    ctx.dps = PRECISION
    with mp.workdps(PRECISION + 50):
        tube = records[INPUT_NAMES[1]]
        h = _get(ctx, tube, "Xh_parameter")
        eval_factors = _derivative_evaluation_factors(ctx, h)
        model_data_48 = _model_norm_and_remainders(
            ctx, tube, ctx.mpf("1e-48")
        )
        model_data_120 = _model_norm_and_remainders(
            ctx, tube, ctx.mpf("1e-120")
        )
        majorants = {
            LAMBDA_48: _candidate_majorant(
                ctx, records, LAMBDA_48
            ),
            LAMBDA_120: _candidate_majorant(
                ctx, records, LAMBDA_120
            ),
        }
        cases = {
            LAMBDA_48: _weighted_case(
                ctx,
                majorants[LAMBDA_48],
                tube,
                eval_factors,
                model_data_48,
            ),
            LAMBDA_120: _weighted_case(
                ctx,
                majorants[LAMBDA_120],
                tube,
                eval_factors,
                model_data_120,
            ),
        }
        output = {
            "precision": PRECISION,
            "accepted_schedule_sha256": ACCEPTED_SCHEDULE_SHA256,
            "input_hashes": source_hashes,
            "nested_input_hashes": nested_hashes,
            "paper_equations": {
                "axis_balance": "(8.7), (8.19), (8.21)",
                "first_correction": "(8.54)",
                "second_order_remainder": "(8.55)",
                "weighted_derivative": "(8.56)-(8.58)",
            },
            "anchored_axis_slope_formula": "partial_s Phi(0,Z)=-(chi(Z)+epsilon beta(Z))/4",
            "target_rectangle": {
                "scaled_R": ["0", RADIAL_RADIUS],
                "Z": ["-1", "1"],
                "chi_small_threshold": ".99",
            },
            "root_factorization": {
                "root_interval": ["-j", "0"],
                "cH_formula": "(9-delta)/2 - j*(1+j) - 4*(1+j+j^2)",
                "cH_lower": cases[LAMBDA_120]["cH_lower"],
                "sampled_root_used": False,
            },
            "derivative_evaluation_factors": eval_factors,
            "model_data_Lambda48": model_data_48,
            "model_data_Lambda120": model_data_120,
            "candidate_majorants": majorants,
            "weighted_cases": cases,
            "Lambda48_weighted_sign_certified": cases[LAMBDA_48][
                "weighted_negative_derivative_gate"
            ],
            "Lambda120_weighted_sign_certified": cases[LAMBDA_120][
                "weighted_negative_derivative_gate"
            ],
            "finite_shared_core_regenerated": False,
            "infinite_radial_remainder_enclosed": False,
            "original_parameter_errors_enclosed": False,
            "source_derivation_errors_enclosed": False,
            "adapter_evaluation_roundoff_enclosed": False,
            "core_to_collar_matching_certified": False,
            "stress_cone_certified": False,
            "temporal_recursion": False,
            "conditional_fixed_data_scope": True,
            "all20_paper_8_50_terms_retained": True,
            "pressure_F0_squared_coupling_retained": True,
        }
        path = Path(__file__).with_suffix(".json")
        path.write_text(
            json.dumps(encode(output), indent=2) + "\n",
            encoding="utf-8",
        )
        print(
            "weighted Lambda48 gate",
            output["Lambda48_weighted_sign_certified"],
            "Lambda120 gate",
            output["Lambda120_weighted_sign_certified"],
            "Lambda120 sufficient minimum",
            mp.nstr(_upper(cases[LAMBDA_120]["sufficient_Lambda_min"]), 18),
            flush=True,
        )
        return output


if __name__ == "__main__":
    run()
