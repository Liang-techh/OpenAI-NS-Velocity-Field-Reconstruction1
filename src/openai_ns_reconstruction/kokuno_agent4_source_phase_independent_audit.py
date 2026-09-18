"""Independent black-box audit for the Kokuno displayed source phase/frame adapter.

This module intentionally consumes only the public ``phase_frame`` value contract from
``KokunoSourcePhaseBoundSignedCurlFamily``.  It does not call Agent-2's tests, reuse
its finite-difference helper, inspect complete-curl internals, or treat this local
source-formula check as a Navier--Stokes residual.

The independent oracle:
* rebuilds the displayed cone/carrier/phase formulas directly from the public source;
* uses a separate nearest-nonzero-integer search implementation;
* differentiates only public phase values with centered FD6 on three frozen steps;
* covers off-grid, axis-near, h/ell perturbations, carrier k=2/3/4 transitions;
* calibrates two deliberately wrong mutations (tilt sign and axial normalization).
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_source_phase_bound_signed_curl import KokunoSourcePhaseBoundSignedCurlFamily

SEED = 9173201
STEPS = (1.0e-2, 5.0e-3, 2.5e-3)
H_VALUES = (0.0035, 0.0075, 0.0095)
ELL_VALUES = (17, 149, 301, 421)
N_CASES = 36

FROZEN_GUARDS = {
    "value_relative_max": 5.0e-13,
    "frame_relative_rms": 5.0e-13,
    "source_identity_relative_max": 5.0e-13,
    "angular_winding_relative_max": 5.0e-13,
    "fd_finest_relative_rms": 1.0e-9,
    "fd_min_refinement_ratio": 20.0,
    "wrong_tilt_mutation_relative_rms_min": 0.30,
    "wrong_axial_normalization_relative_rms_min": 0.20,
}


def _background(R: float, Z: float, coeffs: tuple[float, ...]) -> tuple[float, ...]:
    """Independent smooth nonlinear background used only by this audit."""
    fbase, gbase, a1, a2, a3, b1, b2, b3 = coeffs
    F = (
        fbase
        + a1 * math.sin(17.0 * R)
        + a2 * math.cos(13.0 * Z)
        + a3 * R * Z
        + 0.012 * R * R
    )
    G = (
        gbase
        + b1 * math.cos(11.0 * R)
        + b2 * math.sin(19.0 * Z)
        + b3 * R * Z
        - 0.009 * Z * Z
    )
    F_R = a1 * 17.0 * math.cos(17.0 * R) + a3 * Z + 0.024 * R
    F_Z = -a2 * 13.0 * math.sin(13.0 * Z) + a3 * R
    G_R = -b1 * 11.0 * math.sin(11.0 * R) + b3 * Z
    G_Z = b2 * 19.0 * math.cos(19.0 * Z) + b3 * R - 0.018 * Z
    return F, G, F_R, G_R, F_Z, G_Z


def _nearest_nonzero_reference(x: float) -> int:
    """Independent exhaustive nearest-nonzero search with the declared tie policy."""
    lo = math.floor(x) - 3
    hi = math.ceil(x) + 3
    candidates = [n for n in range(lo, hi + 1) if n != 0]
    distances = [abs(float(n) - x) for n in candidates]
    best = min(distances)
    ties = [n for n, d in zip(candidates, distances) if abs(d - best) <= 2.0e-15]
    if len(ties) == 1:
        return ties[0]
    if x > 0.0:
        return max(ties)
    if x < 0.0:
        return min(ties)
    return 1


def _case_table() -> list[dict[str, Any]]:
    rng = np.random.default_rng(SEED)
    cases: list[dict[str, Any]] = []
    for i in range(N_CASES):
        h = H_VALUES[i % len(H_VALUES)]
        ell = ELL_VALUES[i % len(ELL_VALUES)]
        axis_near = i < 9
        R = (0.036 + 0.018 * rng.random()) if axis_near else (0.18 + rng.random())
        Z = -0.7 + 1.4 * rng.random()
        theta = -2.5 + 5.0 * rng.random()
        R0 = 0.55 + 0.7 * rng.random()
        F0 = 1.0 + 0.5 * rng.random()
        a = 2.6 + 2.2 * rng.random()
        b_s = -0.75 + 1.5 * rng.random()
        u_star = 1.3 + 2.0 * rng.random()
        L_s = 0.55 + 0.65 * rng.random()
        v_plus = L_s * (0.08 + 0.40 * rng.random())
        v_minus = L_s * (0.08 + 0.40 * rng.random())
        coeffs = (
            1.15 + 0.20 * rng.random(),
            0.20 + 0.15 * rng.random(),
            0.11,
            0.08,
            0.025,
            0.09,
            0.07,
            -0.018,
        )
        cases.append(
            {
                "index": i,
                "h": h,
                "ell": ell,
                "axis_near": axis_near,
                "R": R,
                "Z": Z,
                "theta": theta,
                "R0": R0,
                "F0": F0,
                "a": a,
                "b_s": b_s,
                "u_star": u_star,
                "L_s": L_s,
                "v": (v_plus, v_minus),
                "coeffs": coeffs,
            }
        )
    return cases


def _reference(case: dict[str, Any], sigma: float, R: float, Z: float, theta: float, v: float) -> dict[str, Any]:
    h = float(case["h"])
    ell = int(case["ell"])
    R0 = float(case["R0"])
    F0 = float(case["F0"])
    a = float(case["a"])
    b_s = float(case["b_s"])
    u_star = float(case["u_star"])
    L_s = float(case["L_s"])

    epsilon = 2.0 ** (-h * ell)
    k = int(math.ceil(epsilon ** -0.5))
    t_s = -b_s / a
    denom = math.sqrt(1.0 + t_s * t_s)
    N = np.array((-1.0 / denom, -t_s / denom), dtype=float)
    K = np.array((t_s / denom, -1.0 / denom), dtype=float)
    v_s = a * (1.0 + t_s * t_s)
    g0_abs = a * F0 * denom
    g0 = g0_abs * N
    lambda0 = math.sqrt(2.0 * a * F0 * F0 * (1.0 - 2.0 / v_s))
    c0 = lambda0 / (2.0 * F0 * N[0])
    B_s = math.sqrt(
        lambda0 / (epsilon * k * k * (1.0 + u_star * u_star) ** 1.5)
    )

    tilt = K - sigma * u_star * g0 / (L_s * g0_abs * g0_abs)
    pair = B_s * tilt
    tilde_p = R0 * pair[0]
    p_z = pair[1]
    kp = _nearest_nonzero_reference(k * tilde_p)
    p = float(kp) / float(k)
    x0 = sigma * B_s * u_star / 2.0

    F, G, F_R, G_R, F_Z, G_Z = _background(R, Z, tuple(case["coeffs"]))
    H = p * F + p_z * G
    phase = p * theta + p_z * Z / epsilon + x0 * R - v * H
    n_phi = np.array(
        (
            x0 - v * (p * F_R + p_z * G_R),
            p / R,
            p_z - epsilon * v * (p * F_Z + p_z * G_Z),
        ),
        dtype=float,
    )
    return {
        "epsilon": epsilon,
        "k": k,
        "N": N,
        "K": K,
        "v_s": v_s,
        "g0_abs": g0_abs,
        "g0": g0,
        "lambda0": lambda0,
        "c0": c0,
        "B_s": B_s,
        "tilde_p": tilde_p,
        "p_z": p_z,
        "kp": kp,
        "p": p,
        "x0": x0,
        "phase": phase,
        "n_phi": n_phi,
    }


def _public_frame(
    family: KokunoSourcePhaseBoundSignedCurlFamily,
    case: dict[str, Any],
    R: float,
    Z: float,
    theta: float,
) -> dict[str, Any]:
    F, G, F_R, G_R, F_Z, G_Z = _background(R, Z, tuple(case["coeffs"]))
    return family.phase_frame(
        R=R,
        Z=Z,
        theta=theta,
        v=np.asarray((case["v"],), dtype=float),
        beta_labels=((int(case["ell"]), f"heldout-{int(case['index'])}"),),
        R0=float(case["R0"]),
        F0=float(case["F0"]),
        a=float(case["a"]),
        b_s=float(case["b_s"]),
        u_star=float(case["u_star"]),
        L_s=float(case["L_s"]),
        F=F,
        G=G,
        F_R=F_R,
        G_R=G_R,
        F_Z=F_Z,
        G_Z=G_Z,
    )


def _fd6(fun, x: float, step: float) -> float:
    return (
        -fun(x - 3.0 * step)
        + 9.0 * fun(x - 2.0 * step)
        - 45.0 * fun(x - step)
        + 45.0 * fun(x + step)
        - 9.0 * fun(x + 2.0 * step)
        + fun(x + 3.0 * step)
    ) / (60.0 * step)


def _wrong_tilt_n(case: dict[str, Any], sigma: float, R: float, Z: float, v: float) -> np.ndarray:
    """Mutation: reverse the signed source tilt while leaving x0 unchanged."""
    h = float(case["h"])
    ell = int(case["ell"])
    R0 = float(case["R0"])
    F0 = float(case["F0"])
    a = float(case["a"])
    b_s = float(case["b_s"])
    u_star = float(case["u_star"])
    L_s = float(case["L_s"])
    epsilon = 2.0 ** (-h * ell)
    k = int(math.ceil(epsilon ** -0.5))
    t_s = -b_s / a
    denom = math.sqrt(1.0 + t_s * t_s)
    N = np.array((-1.0 / denom, -t_s / denom), dtype=float)
    K = np.array((t_s / denom, -1.0 / denom), dtype=float)
    v_s = a * (1.0 + t_s * t_s)
    g0_abs = a * F0 * denom
    g0 = g0_abs * N
    lambda0 = math.sqrt(2.0 * a * F0 * F0 * (1.0 - 2.0 / v_s))
    B_s = math.sqrt(
        lambda0 / (epsilon * k * k * (1.0 + u_star * u_star) ** 1.5)
    )
    pair = B_s * (K + sigma * u_star * g0 / (L_s * g0_abs * g0_abs))
    tilde_p = R0 * pair[0]
    p_z = pair[1]
    p = _nearest_nonzero_reference(k * tilde_p) / float(k)
    x0 = sigma * B_s * u_star / 2.0
    _, _, F_R, G_R, F_Z, G_Z = _background(R, Z, tuple(case["coeffs"]))
    return np.array(
        (
            x0 - v * (p * F_R + p_z * G_R),
            p / R,
            p_z - epsilon * v * (p * F_Z + p_z * G_Z),
        ),
        dtype=float,
    )


def _rel_scalar(a: float, b: float) -> float:
    return abs(a - b) / max(1.0, abs(b))


def _rel_vector(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b) / max(1.0, float(np.linalg.norm(b))))


def run_audit() -> dict[str, Any]:
    cases = _case_table()
    families = {h: KokunoSourcePhaseBoundSignedCurlFamily(h=h) for h in H_VALUES}

    value_errors: list[float] = []
    frame_errors: list[float] = []
    source_identity_errors: list[float] = []
    angular_winding_errors: list[float] = []
    wrong_tilt_errors: list[float] = []
    wrong_axial_errors: list[float] = []
    fd_errors: dict[float, list[float]] = {step: [] for step in STEPS}
    carrier_values: set[int] = set()
    axis_near_sign_cases = 0
    off_grid_sign_cases = 0
    positive_stencil_cases = 0

    for case in cases:
        family = families[float(case["h"])]
        R = float(case["R"])
        Z = float(case["Z"])
        theta = float(case["theta"])
        base = _public_frame(family, case, R, Z, theta)
        turn = _public_frame(family, case, R, Z, theta + 2.0 * math.pi)

        if R - 3.0 * STEPS[0] > 0.0:
            positive_stencil_cases += 2

        for sidx, sigma in enumerate((1.0, -1.0)):
            v = float(case["v"][sidx])
            ref = _reference(case, sigma, R, Z, theta, v)
            carrier_values.add(int(ref["k"]))
            if bool(case["axis_near"]):
                axis_near_sign_cases += 1
            else:
                off_grid_sign_cases += 1

            public_phase = float(base["phase"][0, sidx])
            public_n = np.asarray(base["n_phi"][0, sidx], dtype=float)
            value_errors.append(_rel_scalar(public_phase, float(ref["phase"])))
            frame_errors.append(_rel_vector(public_n, np.asarray(ref["n_phi"], dtype=float)))

            for public_key, ref_key in (
                ("B_s_by_beta", "B_s"),
                ("lambda0_by_beta", "lambda0"),
                ("c0_by_beta", "c0"),
            ):
                value_errors.append(
                    _rel_scalar(float(base[public_key][0]), float(ref[ref_key]))
                )
            value_errors.append(
                _rel_scalar(float(base["tilde_p_by_beta_sign"][0, sidx]), float(ref["tilde_p"]))
            )
            value_errors.append(
                _rel_scalar(float(base["p_z_by_beta_sign"][0, sidx]), float(ref["p_z"]))
            )
            value_errors.append(
                _rel_scalar(float(base["p_by_beta_sign"][0, sidx]), float(ref["p"]))
            )
            if int(base["k_by_beta"][0]) != int(ref["k"]):
                value_errors.append(1.0)
            if int(base["kp_by_beta_sign"][0, sidx]) != int(ref["kp"]):
                value_errors.append(1.0)

            lhs = float(
                np.dot(
                    np.array((float(ref["tilde_p"]) / float(case["R0"]), float(ref["p_z"]))),
                    np.asarray(ref["g0"], dtype=float),
                )
            )
            rhs = -sigma * float(ref["B_s"]) * float(case["u_star"]) / float(case["L_s"])
            source_identity_errors.append(abs(lhs - rhs) / max(abs(rhs), 1.0e-15))
            source_identity_errors.append(
                abs(float(ref["c0"]) ** 2 - (float(ref["v_s"]) - 2.0) / 2.0)
                / max(1.0, abs((float(ref["v_s"]) - 2.0) / 2.0))
            )

            phase_turn = float(turn["phase"][0, sidx])
            winding = (phase_turn - public_phase) * int(ref["k"]) / (2.0 * math.pi)
            angular_winding_errors.append(
                abs(winding - int(ref["kp"])) / max(1.0, abs(int(ref["kp"])))
            )
            frame_turn = np.asarray(turn["n_phi"][0, sidx], dtype=float)
            angular_winding_errors.append(_rel_vector(frame_turn, public_n))

            wrong_tilt = _wrong_tilt_n(case, sigma, R, Z, v)
            wrong_tilt_errors.append(_rel_vector(wrong_tilt, public_n))
            wrong_axial = np.array(public_n, copy=True)
            wrong_axial[2] = public_n[2] / float(ref["epsilon"])
            wrong_axial_errors.append(_rel_vector(wrong_axial, public_n))

            for step in STEPS:
                def phase_r(rr: float) -> float:
                    return float(_public_frame(family, case, rr, Z, theta)["phase"][0, sidx])

                def phase_theta(tt: float) -> float:
                    return float(_public_frame(family, case, R, Z, tt)["phase"][0, sidx])

                def phase_z(zz: float) -> float:
                    return float(_public_frame(family, case, R, zz, theta)["phase"][0, sidx])

                fd_vec = np.array(
                    (
                        _fd6(phase_r, R, step),
                        _fd6(phase_theta, theta, step) / R,
                        float(ref["epsilon"]) * _fd6(phase_z, Z, step),
                    ),
                    dtype=float,
                )
                fd_errors[step].append(_rel_vector(fd_vec, public_n))

    fd_rms = {
        str(step): float(np.sqrt(np.mean(np.square(fd_errors[step]))))
        for step in STEPS
    }
    refinements = (
        fd_rms[str(STEPS[0])] / fd_rms[str(STEPS[1])],
        fd_rms[str(STEPS[1])] / fd_rms[str(STEPS[2])],
    )
    wrong_tilt_rms = float(np.sqrt(np.mean(np.square(wrong_tilt_errors))))
    wrong_axial_rms = float(np.sqrt(np.mean(np.square(wrong_axial_errors))))

    metrics = {
        "seed": SEED,
        "case_count": len(cases),
        "sign_case_count": 2 * len(cases),
        "axis_near_sign_case_count": axis_near_sign_cases,
        "off_grid_sign_case_count": off_grid_sign_cases,
        "positive_coarse_radial_stencil_sign_case_count": positive_stencil_cases,
        "h_values": list(H_VALUES),
        "ell_values": list(ELL_VALUES),
        "carrier_k_values_observed": sorted(carrier_values),
        "value_relative_max": float(max(value_errors)),
        "frame_relative_rms": float(np.sqrt(np.mean(np.square(frame_errors)))),
        "source_identity_relative_max": float(max(source_identity_errors)),
        "angular_winding_relative_max": float(max(angular_winding_errors)),
        "fd_relative_rms_by_step": fd_rms,
        "fd_refinement_ratios": list(refinements),
        "wrong_tilt_mutation_relative_rms": wrong_tilt_rms,
        "wrong_axial_normalization_relative_rms": wrong_axial_rms,
    }

    guards = {
        "value_reference": metrics["value_relative_max"] <= FROZEN_GUARDS["value_relative_max"],
        "frame_reference": metrics["frame_relative_rms"] <= FROZEN_GUARDS["frame_relative_rms"],
        "source_identities": metrics["source_identity_relative_max"]
        <= FROZEN_GUARDS["source_identity_relative_max"],
        "angular_winding": metrics["angular_winding_relative_max"]
        <= FROZEN_GUARDS["angular_winding_relative_max"],
        "fd_finest": fd_rms[str(STEPS[-1])] <= FROZEN_GUARDS["fd_finest_relative_rms"],
        "fd_refinement": min(refinements) >= FROZEN_GUARDS["fd_min_refinement_ratio"],
        "wrong_tilt_detected": wrong_tilt_rms
        >= FROZEN_GUARDS["wrong_tilt_mutation_relative_rms_min"],
        "wrong_axial_normalization_detected": wrong_axial_rms
        >= FROZEN_GUARDS["wrong_axial_normalization_relative_rms_min"],
        "carrier_transitions_covered": set((2, 3, 4)).issubset(carrier_values),
        "axis_near_covered": axis_near_sign_cases >= 18,
        "off_grid_covered": off_grid_sign_cases >= 54,
        "coarse_radial_stencils_positive": positive_stencil_cases == 2 * len(cases),
    }

    return {
        "schema": "kokuno-agent4-source-phase-independent-audit-v1",
        "task": "KOKUNO-A4-SOURCE-PHASE-INDEPENDENT-AUDIT-031",
        "parent_contract": "KokunoSourcePhaseBoundSignedCurlFamily.phase_frame",
        "independence": {
            "public_values_only": True,
            "agent2_fd_helper_reused": False,
            "complete_curl_helper_reused": False,
            "training_tensor_or_loss_used": False,
            "pressure_or_forcing_fit_used": False,
            "source_schedule_class_reused_by_reference": False,
        },
        "frozen_steps": list(STEPS),
        "frozen_guards": dict(FROZEN_GUARDS),
        "metrics": metrics,
        "guards": guards,
        "local_structural_preflight_passed": bool(all(guards.values())),
        "truth_boundary": {
            "displayed_source_phase_frame_independently_audited": bool(all(guards.values())),
            "actual_positive_order_background_bound": False,
            "actual_source_h_sigma_pulse_integrals_bound": False,
            "actual_signed_auxiliary_rectangles_bound": False,
            "actual_auxiliary_torus_mode_family_bound": False,
            "source_actual_partition_labels_instantiated": False,
            "public_source_bound_velocity_osc_ready": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        },
        "formal_project_gates_unchanged": {
            "normalized_momentum_max_l2": 1.0e-3,
            "divergence_max_l2": 1.0e-5,
            "assessed_here": False,
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args(argv)
    report = run_audit()
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.out is None:
        print(text, end="")
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
        print(args.out)
    return 0 if report["local_structural_preflight_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
