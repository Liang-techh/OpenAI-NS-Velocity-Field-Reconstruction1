"""Independent Agent-4 audit of the supplied signed complete-curl physical family.

This audit deliberately treats ``KokunoSourceSignedCompleteCurlFamily.physical_family``
as a value-only black box.  It does not read its returned amplitudes or derivative
arrays and does not reuse the complete-curl helper in the covariance oracle.

The supplied field below is a deterministic, source-compatible manufactured
family.  It is *not* the unreleased Kokuno positive-order background, pulse data,
signed rectangles, or source partition labels.  The purpose of this increment is
narrow: determine whether the new public physical sign-resolved complete-curl
handoff can retain a genuinely two-dimensional phase-mean covariance response,
and whether an intentionally duplicated sign geometry is detected as rank-deficient.

These local structural guards were frozen before exact-head execution.  They do
not alter or assess the project-wide held-out NS gates (normalized momentum
max/L2 <=1e-3; divergence max/L2 <=1e-5).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_source_signed_complete_curl_family import KokunoSourceSignedCompleteCurlFamily

SEED = 9173191
RESOLUTIONS = (12, 24, 48)
RADII = (0.08, 0.35, 0.90)
TARGET_STEP = 2.0e-4
H = 0.005
BETA_LABELS = ((8, "agent4-heldout-audit"),)

# Frozen before execution.
MIN_PHYSICAL_RANK_RATIO = 2.0e-2
MAX_RESOLUTION_RELATIVE_DRIFT = 1.0e-8
MIN_VELOCITY_RMS = 1.0e-4
MAX_DUPLICATED_SIGN_RANK_RATIO = 1.0e-8

BASE_A_C = 2.0
BASE_U_STAR = 3.0
BASE_H_PLUS = 1.1
BASE_H_MINUS = 0.9
BASE_P = 0.40
BASE_Q = 0.08
BASE_T_N = -BASE_A_C * BASE_P
BASE_T_K = BASE_U_STAR * BASE_Q
DIRECTION_GAP_ETA = 0.10

PLUS_PROTO = np.asarray([1.0 + 0.2j, 0.55 + 0.1j, -1.0 - 0.2j], dtype=np.complex128)
MINUS_PROTO = np.asarray([0.8 - 0.15j, -0.45 + 0.05j, -0.8 + 0.15j], dtype=np.complex128)
N_PHI = np.asarray([1.0, 0.0, 1.0], dtype=float)


def _frozen_offsets() -> tuple[float, float, float]:
    rng = np.random.default_rng(SEED)
    return tuple(float(x) for x in rng.uniform(0.0, 2.0 * np.pi, size=3))


def _manufactured_inputs(n: int, radius: float, *, duplicate_sign_geometry: bool) -> dict[str, Any]:
    if n not in RESOLUTIONS:
        raise ValueError(f"n must be one of {RESOLUTIONS}")
    if radius not in RADII:
        raise ValueError(f"radius must be one of {RADII}")

    off_plus, off_minus, off_theta = _frozen_offsets()
    grid = 2.0 * np.pi * np.arange(n, dtype=float) / float(n)
    plus_phase = np.broadcast_to(grid[:, None] + off_plus, (n, n))
    minus_phase = np.broadcast_to(grid[None, :] + off_minus, (n, n))
    theta = np.mod(0.31 * plus_phase + 0.17 * minus_phase + off_theta, 2.0 * np.pi)

    phase = np.empty((n, n, 1, 2), dtype=float)
    phase[..., 0, 0] = plus_phase
    phase[..., 0, 1] = minus_phase

    n_phi = np.broadcast_to(N_PHI, (n, n, 1, 2, 3)).copy()
    t_proto = np.empty((n, n, 1, 2, 3), dtype=np.complex128)
    t_proto[..., 0, 0, :] = PLUS_PROTO
    t_proto[..., 0, 1, :] = PLUS_PROTO if duplicate_sign_geometry else MINUS_PROTO

    zeros_vec = np.zeros_like(t_proto)
    eta = np.ones((n, n, 1), dtype=float)
    zeros_beta = np.zeros_like(eta)
    R = np.full((n, n, 1), float(radius), dtype=float)

    return {
        "R": R,
        "theta": theta,
        "phase": phase,
        "n_phi": n_phi,
        "t_plus_prototype": t_proto,
        "D_r_C_plus_prototype": zeros_vec,
        "D_z_C_plus_prototype": zeros_vec,
        "eta": eta,
        "D_r_eta": zeros_beta,
        "D_z_eta": zeros_beta,
        "beta_labels": BETA_LABELS,
        "A_c": BASE_A_C,
        "u_star": BASE_U_STAR,
        "h_plus": BASE_H_PLUS,
        "h_minus": BASE_H_MINUS,
        "T_N": BASE_T_N,
        "T_K": BASE_T_K,
        "D_r_A_c": 0.0,
        "D_r_u_star": 0.0,
        "D_r_h_plus": 0.0,
        "D_r_h_minus": 0.0,
        "D_r_T_N": 0.0,
        "D_r_T_K": 0.0,
        "D_z_A_c": 0.0,
        "D_z_u_star": 0.0,
        "D_z_h_plus": 0.0,
        "D_z_h_minus": 0.0,
        "D_z_T_N": 0.0,
        "D_z_T_K": 0.0,
        "direction_gap_eta": DIRECTION_GAP_ETA,
    }


def _public_velocity(
    family: KokunoSourceSignedCompleteCurlFamily,
    args: dict[str, Any],
    T_N: float,
    T_K: float,
) -> tuple[np.ndarray, np.ndarray]:
    call = dict(args)
    call["T_N"] = float(T_N)
    call["T_K"] = float(T_K)
    # Value-only black-box consumption: no returned amplitudes/tangents enter the oracle.
    out = family.physical_family(**call)
    cart = np.asarray(out["velocity_physical_cartesian_total"], dtype=float)
    theta = np.asarray(call["theta"], dtype=float)
    if cart.shape != theta.shape + (3,):
        raise RuntimeError("public Cartesian total has an unexpected shape")
    return cart, theta


def _covariance_from_public_cartesian(cart: np.ndarray, theta: np.ndarray) -> tuple[np.ndarray, float]:
    c = np.cos(theta)
    s = np.sin(theta)
    ux, uy, uz = cart[..., 0], cart[..., 1], cart[..., 2]
    ur = c * ux + s * uy
    utheta = -s * ux + c * uy
    covariance = np.asarray(
        [np.mean(ur * utheta, dtype=np.float64), np.mean(ur * uz, dtype=np.float64)],
        dtype=float,
    )
    velocity_rms = float(np.sqrt(np.mean(np.sum(cart * cart, axis=-1), dtype=np.float64)))
    return covariance, velocity_rms


def _covariance(
    family: KokunoSourceSignedCompleteCurlFamily,
    args: dict[str, Any],
    T_N: float,
    T_K: float,
) -> tuple[np.ndarray, float]:
    cart, theta = _public_velocity(family, args, T_N, T_K)
    return _covariance_from_public_cartesian(cart, theta)


def _target_jacobian(
    family: KokunoSourceSignedCompleteCurlFamily,
    args: dict[str, Any],
) -> tuple[np.ndarray, float]:
    d = TARGET_STEP
    c_np, _ = _covariance(family, args, BASE_T_N + d, BASE_T_K)
    c_nm, _ = _covariance(family, args, BASE_T_N - d, BASE_T_K)
    c_kp, _ = _covariance(family, args, BASE_T_N, BASE_T_K + d)
    c_km, _ = _covariance(family, args, BASE_T_N, BASE_T_K - d)
    c0, velocity_rms = _covariance(family, args, BASE_T_N, BASE_T_K)
    jac = np.column_stack(((c_np - c_nm) / (2.0 * d), (c_kp - c_km) / (2.0 * d)))
    if not np.all(np.isfinite(jac)) or not np.all(np.isfinite(c0)):
        raise RuntimeError("non-finite covariance audit result")
    return jac, velocity_rms


def _svd_metrics(jac: np.ndarray) -> tuple[np.ndarray, float]:
    singular_values = np.linalg.svd(jac, compute_uv=False)
    if singular_values.shape != (2,):
        raise RuntimeError("expected a 2x2 covariance Jacobian")
    smax = float(singular_values[0])
    ratio = 0.0 if smax == 0.0 else float(singular_values[1] / smax)
    return singular_values, ratio


def build_report() -> dict[str, Any]:
    family = KokunoSourceSignedCompleteCurlFamily(h=H)
    cases: list[dict[str, Any]] = []
    all_rank_ratios: list[float] = []
    all_velocity_rms: list[float] = []
    all_mutation_ratios: list[float] = []
    all_drifts: list[float] = []

    for radius in RADII:
        rows: list[dict[str, Any]] = []
        baseline_jacobians: dict[int, np.ndarray] = {}
        for n in RESOLUTIONS:
            base_args = _manufactured_inputs(n, radius, duplicate_sign_geometry=False)
            mutation_args = _manufactured_inputs(n, radius, duplicate_sign_geometry=True)

            jac, velocity_rms = _target_jacobian(family, base_args)
            mutation_jac, _ = _target_jacobian(family, mutation_args)
            singular_values, rank_ratio = _svd_metrics(jac)
            mutation_singular_values, mutation_rank_ratio = _svd_metrics(mutation_jac)

            baseline_jacobians[n] = jac
            all_rank_ratios.append(rank_ratio)
            all_velocity_rms.append(velocity_rms)
            all_mutation_ratios.append(mutation_rank_ratio)
            rows.append(
                {
                    "resolution": int(n),
                    "jacobian": jac.tolist(),
                    "singular_values": singular_values.tolist(),
                    "rank_ratio": rank_ratio,
                    "velocity_rms": velocity_rms,
                    "duplicated_sign_mutation_jacobian": mutation_jac.tolist(),
                    "duplicated_sign_mutation_singular_values": mutation_singular_values.tolist(),
                    "duplicated_sign_mutation_rank_ratio": mutation_rank_ratio,
                }
            )

        finest = baseline_jacobians[RESOLUTIONS[-1]]
        finest_norm = max(float(np.linalg.norm(finest)), np.finfo(float).tiny)
        drift_rows = []
        for n in RESOLUTIONS[:-1]:
            drift = float(np.linalg.norm(baseline_jacobians[n] - finest) / finest_norm)
            all_drifts.append(drift)
            drift_rows.append({"resolution": int(n), "relative_frobenius_drift_to_finest": drift})
        cases.append({"radius": float(radius), "resolutions": rows, "resolution_stability": drift_rows})

    metrics = {
        "minimum_physical_rank_ratio": float(min(all_rank_ratios)),
        "maximum_resolution_relative_drift": float(max(all_drifts)),
        "minimum_velocity_rms": float(min(all_velocity_rms)),
        "maximum_duplicated_sign_mutation_rank_ratio": float(max(all_mutation_ratios)),
    }
    guards = {
        "minimum_physical_rank_ratio": MIN_PHYSICAL_RANK_RATIO,
        "maximum_resolution_relative_drift": MAX_RESOLUTION_RELATIVE_DRIFT,
        "minimum_velocity_rms": MIN_VELOCITY_RMS,
        "maximum_duplicated_sign_mutation_rank_ratio": MAX_DUPLICATED_SIGN_RANK_RATIO,
    }
    passed = bool(
        metrics["minimum_physical_rank_ratio"] >= MIN_PHYSICAL_RANK_RATIO
        and metrics["maximum_resolution_relative_drift"] <= MAX_RESOLUTION_RELATIVE_DRIFT
        and metrics["minimum_velocity_rms"] >= MIN_VELOCITY_RMS
        and metrics["maximum_duplicated_sign_mutation_rank_ratio"] <= MAX_DUPLICATED_SIGN_RANK_RATIO
    )
    return {
        "schema": "kokuno-agent4-signed-complete-curl-covariance-audit-v1",
        "seed": SEED,
        "parent_exact_head": "bcf983ff0b8e059b9f5d2b051b46e839a9afb2c5",
        "resolutions": list(RESOLUTIONS),
        "radii": list(RADII),
        "target_step": TARGET_STEP,
        "source_compatible_manufactured_input": True,
        "actual_source_input": False,
        "oracle": {
            "candidate_consumption": "public velocity_physical_cartesian_total values only",
            "covariance": "independently rotate Cartesian total back to cylindrical and phase-average ur*utheta, ur*uz",
            "target_derivative": "centered finite difference in public T_N/T_K inputs",
            "complete_curl_helper_reused_by_oracle": False,
            "returned_candidate_amplitudes_or_tangents_read_by_oracle": False,
        },
        "guards_frozen_before_execution": guards,
        "metrics": metrics,
        "cases": cases,
        "local_physical_covariance_preflight_passed": passed,
        "truth_boundary": {
            "supplied_manufactured_physical_signed_covariance_rank_two_assessed": True,
            "actual_source_positive_order_background_bound": False,
            "actual_source_h_sigma_pulse_integrals_bound": False,
            "actual_source_signed_rectangles_or_modes_bound": False,
            "actual_source_partition_labels_instantiated": False,
            "public_source_bound_xyz_t_oscillatory_velocity_ready": False,
            "actual_source_physical_covariance_rank_two_assessed": False,
            "genuinely_independent_second_covariance_column_ready": False,
            "correction_ready": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "formal_project_gates_unchanged": {
            "normalized_momentum_max_and_L2": 1.0e-3,
            "divergence_max_and_L2": 1.0e-5,
            "assessed_here": False,
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    report = build_report()
    text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    # Scientific PASS/REJECT is encoded in the report instead of suppressing artifact publication.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
