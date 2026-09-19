"""Independent Agent-4 audit of the selected PA.10 pressure primitive seam.

This validator is intentionally black-box with respect to Agent-1 construction
internals. It consumes only ``KokunoPA10SelectedPressurePrimitive.evaluate``
(value and declared derivative outputs) plus immutable public candidate
identity metadata. Independent finite-difference and composite-Simpson paths
are used to audit the new pressure/mixed-derivative interface.

This is NOT a full Navier--Stokes validation. It cannot promote global leading
pressure, full same-cycle correction, held-out momentum residual, or
``pde_validated``.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_pa10_selected_pressure_primitive import KokunoPA10SelectedPressurePrimitive


SCHEMA = "kokuno-agent4-pa10-pressure-primitive-independent-audit-v1"
AGENT1_PR = 606
AGENT1_HEAD = "56d7e47e4033d869f802ee51ed380b1daa4a34a1"
AGENT1_CANDIDATE_SHA256 = "8b5cdfb749d395f25d7c5f09f9200a3b281dd346fcb4ead6c8c5b1583d9ff675"

SEED = 9173271
SAMPLE_COUNT = 24
Y_STEPS = (0.04, 0.02, 0.01)
PHI_ETA_STEPS = (0.02, 0.01, 0.005)
P_ETA_WIDTH_FRACTIONS = (0.20, 0.10, 0.05)
SIMPSON_COUNTS = (65, 129, 257)

# Frozen independently from the candidate implementation, from the exact #606
# receipt. They are used only to place pressure-eta probes on the extremely
# narrow non-underflow scale; they are not derivative oracles.
SELECTED_LAMBDA = 2.503192875997364e27
PHASE_STATIONARY_ETA = -0.004449378508888404

# Local interface gates frozen before the exact-head CI run. These are NOT the
# project PDE gates.
LOCAL_GATES = {
    "radial_relative_rms": 2.0e-5,
    "radial_relative_max": 1.0e-4,
    "radial_min_refinement_ratio": 2.0,
    "phi_eta_relative_rms": 2.0e-5,
    "phi_eta_relative_max": 1.0e-4,
    "phi_eta_min_refinement_ratio": 2.0,
    "p_eta_relative_rms": 2.0e-3,
    "p_eta_relative_max": 5.0e-3,
    "p_eta_min_refinement_ratio": 2.0,
    "fundamental_theorem_relative_rms": 1.0e-8,
    "fundamental_theorem_relative_max": 5.0e-8,
    "origin_zero_abs_max": 1.0e-12,
    "mutation_relative_rms_min": 0.5,
    "nontrivial_derivative_rms_min": 1.0e-8,
}

FINAL_PROJECT_GATES = {
    "normalized_momentum_max": 1.0e-3,
    "normalized_momentum_L2": 1.0e-3,
    "divergence_max": 1.0e-5,
    "divergence_L2": 1.0e-5,
}


def _rms(value: np.ndarray) -> float:
    array = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(array * array)))


def _metric(estimate: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    estimate = np.asarray(estimate, dtype=float)
    reference = np.asarray(reference, dtype=float)
    error = estimate - reference
    rms_ref = max(_rms(reference), 1.0e-300)
    max_ref = max(float(np.max(np.abs(reference))), 1.0e-300)
    return {
        "absolute_rms": _rms(error),
        "absolute_max": float(np.max(np.abs(error))),
        "reference_rms": _rms(reference),
        "reference_max": float(np.max(np.abs(reference))),
        "relative_rms": _rms(error) / rms_ref,
        "relative_max": float(np.max(np.abs(error))) / max_ref,
    }


def _refinement_ratios(levels: list[dict[str, float]]) -> list[float]:
    errors = [float(level["relative_rms"]) for level in levels]
    out: list[float] = []
    for coarse, fine in zip(errors[:-1], errors[1:], strict=True):
        out.append(float(coarse / max(fine, 1.0e-300)))
    return out


def _fd4_y(
    primitive: KokunoPA10SelectedPressurePrimitive,
    field: str,
    y: np.ndarray,
    eta: np.ndarray,
    h: float,
) -> np.ndarray:
    def values(offset: float) -> np.ndarray:
        return np.asarray(primitive.evaluate(y + offset * h, eta)[field], dtype=float)

    return (-values(2.0) + 8.0 * values(1.0) - 8.0 * values(-1.0) + values(-2.0)) / (
        12.0 * h
    )


def _fd4_eta(
    primitive: KokunoPA10SelectedPressurePrimitive,
    field: str,
    y: np.ndarray,
    eta: np.ndarray,
    h: float,
) -> np.ndarray:
    def values(offset: float) -> np.ndarray:
        return np.asarray(primitive.evaluate(y, eta + offset * h)[field], dtype=float)

    return (-values(2.0) + 8.0 * values(1.0) - 8.0 * values(-1.0) + values(-2.0)) / (
        12.0 * h
    )


def _simpson(values: np.ndarray, upper: float) -> float:
    values = np.asarray(values, dtype=float)
    n = int(values.size)
    if n < 3 or n % 2 != 1:
        raise ValueError("composite Simpson requires an odd node count >=3")
    h = float(upper) / float(n - 1)
    return float(
        h
        / 3.0
        * (
            values[0]
            + values[-1]
            + 4.0 * np.sum(values[1:-1:2])
            + 2.0 * np.sum(values[2:-1:2])
        )
    )


def _fundamental_theorem_level(
    primitive: KokunoPA10SelectedPressurePrimitive,
    y_targets: np.ndarray,
    eta_targets: np.ndarray,
    node_count: int,
) -> dict[str, float]:
    estimates: list[float] = []
    references: list[float] = []
    for upper, eta in zip(y_targets, eta_targets, strict=True):
        nodes = np.linspace(0.0, float(upper), int(node_count), dtype=float)
        eta_nodes = np.full(nodes.shape, float(eta), dtype=float)
        p_y = np.asarray(primitive.evaluate(nodes, eta_nodes)["p_Y"], dtype=float)
        estimates.append(_simpson(p_y, float(upper)))
        endpoint = float(primitive.evaluate(float(upper), float(eta))["p"])
        origin = float(primitive.evaluate(0.0, float(eta))["p"])
        references.append(endpoint - origin)
    return _metric(np.asarray(estimates), np.asarray(references))


def run_audit() -> dict[str, Any]:
    primitive = KokunoPA10SelectedPressurePrimitive(quadrature_points=64)
    if primitive.sha256 != AGENT1_CANDIDATE_SHA256:
        raise RuntimeError(
            "Agent-1 candidate identity moved: "
            f"{primitive.sha256} != {AGENT1_CANDIDATE_SHA256}"
        )

    rng = np.random.default_rng(SEED)
    width = float(SELECTED_LAMBDA ** -0.5)
    eta_ulp = float(
        abs(
            np.nextafter(PHASE_STATIONARY_ETA, math.inf)
            - PHASE_STATIONARY_ETA
        )
    )

    # Public pressure radial derivative and Y*Phi_Y at fresh off-grid points.
    y_radial = rng.uniform(0.35, 3.75, size=SAMPLE_COUNT)
    eta_radial = PHASE_STATIONARY_ETA + rng.uniform(
        -0.60, 0.60, size=SAMPLE_COUNT
    ) * width
    radial_public = primitive.evaluate(y_radial, eta_radial)
    ref_p_y = np.asarray(radial_public["p_Y"], dtype=float)
    ref_y_phi_y = np.asarray(radial_public["Y_Phi_Y"], dtype=float)

    p_y_levels: list[dict[str, float]] = []
    y_phi_y_levels: list[dict[str, float]] = []
    finest_p_y_estimate: np.ndarray | None = None
    for h in Y_STEPS:
        p_y_fd = _fd4_y(primitive, "p", y_radial, eta_radial, h)
        phi_y_fd = _fd4_y(primitive, "Phi", y_radial, eta_radial, h)
        y_phi_y_fd = y_radial * phi_y_fd
        p_y_levels.append({"step": h, **_metric(p_y_fd, ref_p_y)})
        y_phi_y_levels.append({"step": h, **_metric(y_phi_y_fd, ref_y_phi_y)})
        finest_p_y_estimate = p_y_fd

    p_y_ratios = _refinement_ratios(p_y_levels)
    y_phi_y_ratios = _refinement_ratios(y_phi_y_levels)

    # Public Phi_eta on a generic eta cloud. This avoids using the pressure
    # layer's Lambda-narrow scale to make a slow profile derivative look small.
    y_phi_eta = rng.uniform(0.35, 3.75, size=SAMPLE_COUNT)
    eta_phi_eta = rng.uniform(-0.72, 0.72, size=SAMPLE_COUNT)
    ref_phi_eta = np.asarray(
        primitive.evaluate(y_phi_eta, eta_phi_eta)["Phi_eta"], dtype=float
    )
    phi_eta_levels: list[dict[str, float]] = []
    for h in PHI_ETA_STEPS:
        estimate = _fd4_eta(primitive, "Phi", y_phi_eta, eta_phi_eta, h)
        phi_eta_levels.append({"step": h, **_metric(estimate, ref_phi_eta)})
    phi_eta_ratios = _refinement_ratios(phi_eta_levels)

    # Public p_eta on the independent Lambda^{-1/2} probe scale. Direct
    # finite differences use only p-values from the public evaluator.
    y_p_eta = rng.uniform(0.45, 3.70, size=SAMPLE_COUNT)
    eta_p_eta = PHASE_STATIONARY_ETA + rng.uniform(
        -0.65, 0.65, size=SAMPLE_COUNT
    ) * width
    ref_p_eta = np.asarray(
        primitive.evaluate(y_p_eta, eta_p_eta)["p_eta"], dtype=float
    )
    p_eta_levels: list[dict[str, float]] = []
    finest_p_eta_estimate: np.ndarray | None = None
    eta_steps: list[float] = []
    for fraction in P_ETA_WIDTH_FRACTIONS:
        h = float(fraction * width)
        eta_steps.append(h)
        estimate = _fd4_eta(primitive, "p", y_p_eta, eta_p_eta, h)
        p_eta_levels.append(
            {
                "step": h,
                "width_fraction": fraction,
                "step_in_local_ulps": h / max(eta_ulp, 1.0e-300),
                **_metric(estimate, ref_p_eta),
            }
        )
        finest_p_eta_estimate = estimate
    p_eta_ratios = _refinement_ratios(p_eta_levels)

    # Independent integral identity: integral_0^Y p_Y dY = p(Y)-p(0).
    y_base = np.asarray([0.73, 1.61, 2.84, 3.97], dtype=float)
    eta_offsets = np.asarray([-0.45, 0.0, 0.45], dtype=float)
    y_ft = np.repeat(y_base, eta_offsets.size)
    eta_ft = np.tile(
        PHASE_STATIONARY_ETA + eta_offsets * width,
        y_base.size,
    )
    fto_levels = [
        {
            "node_count": int(count),
            **_fundamental_theorem_level(primitive, y_ft, eta_ft, int(count)),
        }
        for count in SIMPSON_COUNTS
    ]

    origin_eta = np.asarray(
        [
            PHASE_STATIONARY_ETA,
            PHASE_STATIONARY_ETA - 0.5 * width,
            PHASE_STATIONARY_ETA + 0.5 * width,
            -0.41,
            0.37,
        ],
        dtype=float,
    )
    origin_y = np.zeros_like(origin_eta)
    origin = primitive.evaluate(origin_y, origin_eta)
    origin_zero_abs_max = float(
        max(
            np.max(np.abs(np.asarray(origin["p"], dtype=float))),
            np.max(np.abs(np.asarray(origin["Y_p_Y"], dtype=float))),
            np.max(np.abs(np.asarray(origin["Y_Phi_Y"], dtype=float))),
        )
    )

    assert finest_p_y_estimate is not None
    assert finest_p_eta_estimate is not None
    p_y_sign_flip = _metric(finest_p_y_estimate, -ref_p_y)["relative_rms"]
    p_eta_sign_flip = _metric(finest_p_eta_estimate, -ref_p_eta)["relative_rms"]

    guards: dict[str, bool] = {
        "candidate_identity_bound": primitive.sha256 == AGENT1_CANDIDATE_SHA256,
        "p_y_finest_relative_rms": p_y_levels[-1]["relative_rms"]
        <= LOCAL_GATES["radial_relative_rms"],
        "p_y_finest_relative_max": p_y_levels[-1]["relative_max"]
        <= LOCAL_GATES["radial_relative_max"],
        "p_y_refinement": min(p_y_ratios)
        >= LOCAL_GATES["radial_min_refinement_ratio"],
        "Y_Phi_Y_finest_relative_rms": y_phi_y_levels[-1]["relative_rms"]
        <= LOCAL_GATES["radial_relative_rms"],
        "Y_Phi_Y_finest_relative_max": y_phi_y_levels[-1]["relative_max"]
        <= LOCAL_GATES["radial_relative_max"],
        "Y_Phi_Y_refinement": min(y_phi_y_ratios)
        >= LOCAL_GATES["radial_min_refinement_ratio"],
        "Phi_eta_finest_relative_rms": phi_eta_levels[-1]["relative_rms"]
        <= LOCAL_GATES["phi_eta_relative_rms"],
        "Phi_eta_finest_relative_max": phi_eta_levels[-1]["relative_max"]
        <= LOCAL_GATES["phi_eta_relative_max"],
        "Phi_eta_refinement": min(phi_eta_ratios)
        >= LOCAL_GATES["phi_eta_min_refinement_ratio"],
        "p_eta_finest_relative_rms": p_eta_levels[-1]["relative_rms"]
        <= LOCAL_GATES["p_eta_relative_rms"],
        "p_eta_finest_relative_max": p_eta_levels[-1]["relative_max"]
        <= LOCAL_GATES["p_eta_relative_max"],
        "p_eta_refinement": min(p_eta_ratios)
        >= LOCAL_GATES["p_eta_min_refinement_ratio"],
        "fundamental_theorem_finest_relative_rms": fto_levels[-1]["relative_rms"]
        <= LOCAL_GATES["fundamental_theorem_relative_rms"],
        "fundamental_theorem_finest_relative_max": fto_levels[-1]["relative_max"]
        <= LOCAL_GATES["fundamental_theorem_relative_max"],
        "origin_zero": origin_zero_abs_max <= LOCAL_GATES["origin_zero_abs_max"],
        "p_y_nontrivial": _rms(ref_p_y)
        >= LOCAL_GATES["nontrivial_derivative_rms_min"],
        "p_eta_nontrivial": _rms(ref_p_eta)
        >= LOCAL_GATES["nontrivial_derivative_rms_min"],
        "p_y_sign_flip_detected": p_y_sign_flip
        >= LOCAL_GATES["mutation_relative_rms_min"],
        "p_eta_sign_flip_detected": p_eta_sign_flip
        >= LOCAL_GATES["mutation_relative_rms_min"],
    }
    failed_guards = [name for name, passed in guards.items() if not passed]
    passed = not failed_guards

    return {
        "schema": SCHEMA,
        "task": "KOKUNO-A4-PA10-PRESSURE-PRIMITIVE-INDEPENDENT-AUDIT-045",
        "provenance": {
            "agent1_pr": AGENT1_PR,
            "agent1_exact_head": AGENT1_HEAD,
            "agent1_candidate_sha256": AGENT1_CANDIDATE_SHA256,
        },
        "independence_contract": {
            "public_api_consumed": (
                "KokunoPA10SelectedPressurePrimitive.evaluate(Y,eta)"
            ),
            "agent1_private_integral_helper_called": False,
            "agent1_axis_state_called_by_validator": False,
            "agent1_f0_prime_called_by_validator": False,
            "agent1_quadrature_used_as_validation_operator": False,
            "independent_operators": [
                "centered_FD4_in_Y_on_public_p_and_Phi",
                "centered_FD4_in_eta_on_public_p_and_Phi",
                "composite_Simpson_in_Y_on_public_p_Y",
            ],
        },
        "frozen_protocol": {
            "seed": SEED,
            "sample_count_per_fd_cloud": SAMPLE_COUNT,
            "Y_steps": list(Y_STEPS),
            "Phi_eta_steps": list(PHI_ETA_STEPS),
            "p_eta_width_fractions": list(P_ETA_WIDTH_FRACTIONS),
            "pressure_eta_probe_width": width,
            "pressure_eta_steps": eta_steps,
            "phase_stationary_eta": PHASE_STATIONARY_ETA,
            "selected_lambda": SELECTED_LAMBDA,
            "simpson_node_counts": list(SIMPSON_COUNTS),
            "local_gates": dict(LOCAL_GATES),
            "final_project_gates_unchanged": dict(FINAL_PROJECT_GATES),
        },
        "results": {
            "p_Y": {
                "levels": p_y_levels,
                "refinement_ratios": p_y_ratios,
            },
            "Y_Phi_Y": {
                "levels": y_phi_y_levels,
                "refinement_ratios": y_phi_y_ratios,
            },
            "Phi_eta": {
                "levels": phi_eta_levels,
                "refinement_ratios": phi_eta_ratios,
            },
            "p_eta": {
                "levels": p_eta_levels,
                "refinement_ratios": p_eta_ratios,
                "eta_ulp_at_center": eta_ulp,
            },
            "fundamental_theorem_Y": {
                "levels": fto_levels,
            },
            "origin_zero_abs_max": origin_zero_abs_max,
            "nontrivial_reference_rms": {
                "p_Y": _rms(ref_p_y),
                "Y_Phi_Y": _rms(ref_y_phi_y),
                "Phi_eta": _rms(ref_phi_eta),
                "p_eta": _rms(ref_p_eta),
            },
            "mutation_relative_rms": {
                "p_Y_sign_flip": p_y_sign_flip,
                "p_eta_sign_flip": p_eta_sign_flip,
            },
        },
        "guards": guards,
        "failed_guards": failed_guards,
        "selected_pressure_primitive_independent_preflight_passed": passed,
        "truth_boundary": {
            "selected_center_pressure_primitive_independently_audited": passed,
            "source_pressure_radius_one_ball_bound_assessed": False,
            "global_pressure_matched": False,
            "global_leading_profile_reconstructed": False,
            "full_same_cycle_requested_stress_materialized": False,
            "public_velocity_correction_materialized": False,
            "leading_only_ns_residual_assessed": False,
            "leading_plus_oscillatory_ns_residual_assessed": False,
            "after_correction_ns_residual_assessed": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)
    receipt = run_audit()
    text = json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
