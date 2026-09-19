"""Independent Agent-4 re-audit of the repaired PA.10 public pressure value path.

This validator is deliberately black-box with respect to Agent-1 construction
internals.  It consumes only ``KokunoPA10SelectedPressurePrimitive.evaluate``
and immutable candidate identity metadata.  In particular it does not call the
candidate's local-phase quadrature, stationary-factorization helpers, radial
quadrature helper, or any training/construction tensor.

The narrow question is whether the repaired public *value* ``p(Y, eta)`` has a
numerical eta derivative consistent with the declared public ``p_eta`` on the
natural ``Lambda**(-1/2)`` layer that Agent 4 #608 previously rejected.  A
centered Cartesian-style FD6 operator and an independent Gauss--Legendre
fundamental-theorem check are used here, rather than Agent 1's FD4 self-check.

This is a selected-center interface audit only.  It is NOT a full
Navier--Stokes residual evaluation and cannot promote source radius-one-ball
bounds, global matched pressure/leading velocity, correction readiness, or
``pde_validated``.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.polynomial.legendre import leggauss

from .kokuno_pa10_selected_pressure_primitive import KokunoPA10SelectedPressurePrimitive


SCHEMA = "kokuno-agent4-pa10-local-phase-pressure-independent-reaudit-v1"
AGENT1_PR = 615
AGENT1_HEAD = "7bd597398aef8b4ef1671f79e39eb352adacdd0a"
AGENT1_CANDIDATE_SHA256 = "57ad627b55d9efd4c3669796dc80b9f67bc2af5dc42ce65d5a019d19392930e6"

# Fresh held-out cloud.  None of these samples appear in Agent 1's repair
# receipt or Agent 4 #608's prior seed.
SEED = 9173281
SAMPLE_COUNT = 32

# Keep the exact #608 public-value layer ladder so the repair is judged against
# the same numerical scale, but use a different sixth-order derivative stencil.
P_ETA_WIDTH_FRACTIONS = (0.20, 0.10, 0.05)
ETA_INTEGRAL_ORDERS = (16, 32, 64)

# Frozen independently from implementation details, from the public #615
# candidate receipt.  These locate the physical layer; they are not derivative
# oracles.
SELECTED_LAMBDA = 2.503192875997364e27
PHASE_STATIONARY_ETA = -0.004449378508888404

# Preserve #608's p_eta acceptance thresholds unchanged.  No result-dependent
# relaxation is permitted.  The integral identity is an additional independent
# guard, not a replacement for the finite-difference ladder.
LOCAL_GATES = {
    "p_eta_relative_rms": 2.0e-3,
    "p_eta_relative_max": 5.0e-3,
    "p_eta_min_refinement_ratio": 2.0,
    "eta_fundamental_theorem_relative_rms": 1.0e-8,
    "eta_fundamental_theorem_relative_max": 5.0e-8,
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
    return [
        float(coarse / max(fine, 1.0e-300))
        for coarse, fine in zip(errors[:-1], errors[1:], strict=True)
    ]


def _fd6_eta(
    primitive: KokunoPA10SelectedPressurePrimitive,
    field: str,
    y: np.ndarray,
    eta: np.ndarray,
    h: float,
) -> np.ndarray:
    """Centered sixth-order first derivative using public values only."""

    def values(offset: float) -> np.ndarray:
        return np.asarray(primitive.evaluate(y, eta + offset * h)[field], dtype=float)

    return (
        -values(-3.0)
        + 9.0 * values(-2.0)
        - 45.0 * values(-1.0)
        + 45.0 * values(1.0)
        - 9.0 * values(2.0)
        + values(3.0)
    ) / (60.0 * h)


def _eta_integral_level(
    primitive: KokunoPA10SelectedPressurePrimitive,
    y: np.ndarray,
    eta_a: np.ndarray,
    eta_b: np.ndarray,
    order: int,
) -> dict[str, float]:
    """Check int_a^b p_eta d eta = p(b)-p(a) with independent GL quadrature."""

    nodes, weights = leggauss(int(order))
    estimates: list[float] = []
    references: list[float] = []
    for local_y, lower, upper in zip(y, eta_a, eta_b, strict=True):
        half = 0.5 * (float(upper) - float(lower))
        midpoint = 0.5 * (float(upper) + float(lower))
        eta_nodes = midpoint + half * nodes
        y_nodes = np.full(eta_nodes.shape, float(local_y), dtype=float)
        derivative = np.asarray(
            primitive.evaluate(y_nodes, eta_nodes)["p_eta"], dtype=float
        )
        estimates.append(half * float(np.dot(weights, derivative)))
        p_upper = float(primitive.evaluate(float(local_y), float(upper))["p"])
        p_lower = float(primitive.evaluate(float(local_y), float(lower))["p"])
        references.append(p_upper - p_lower)
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
        abs(np.nextafter(PHASE_STATIONARY_ETA, math.inf) - PHASE_STATIONARY_ETA)
    )

    # Keep the held-out center offsets away from the exact stationary root and
    # leave room for the widest +/-3h FD6 stencil entirely inside the narrow
    # layer neighborhood.
    y = rng.uniform(0.45, 3.70, size=SAMPLE_COUNT)
    eta_offsets = rng.uniform(-0.28, 0.28, size=SAMPLE_COUNT)
    eta = PHASE_STATIONARY_ETA + eta_offsets * width
    public = primitive.evaluate(y, eta)
    ref_p_eta = np.asarray(public["p_eta"], dtype=float)

    p_eta_levels: list[dict[str, float]] = []
    finest_estimate: np.ndarray | None = None
    for fraction in P_ETA_WIDTH_FRACTIONS:
        h = float(fraction * width)
        estimate = _fd6_eta(primitive, "p", y, eta, h)
        p_eta_levels.append(
            {
                "step": h,
                "width_fraction": float(fraction),
                "step_in_local_ulps": h / max(eta_ulp, 1.0e-300),
                **_metric(estimate, ref_p_eta),
            }
        )
        finest_estimate = estimate
    p_eta_ratios = _refinement_ratios(p_eta_levels)

    # Disjoint deterministic intervals.  They are intentionally asymmetric so
    # p(b)-p(a) is nonzero and the relative fundamental-theorem metric remains
    # informative rather than being normalized by near-cancellation.
    interval_width_pairs = np.asarray(
        [
            (-0.70, -0.24),
            (-0.62, -0.12),
            (-0.48, -0.04),
            (0.04, 0.48),
            (0.12, 0.62),
            (0.24, 0.70),
        ],
        dtype=float,
    )
    y_integral = np.asarray([0.61, 1.07, 1.83, 2.41, 3.02, 3.61], dtype=float)
    eta_a = PHASE_STATIONARY_ETA + interval_width_pairs[:, 0] * width
    eta_b = PHASE_STATIONARY_ETA + interval_width_pairs[:, 1] * width
    integral_levels = [
        {
            "quadrature_order": int(order),
            **_eta_integral_level(primitive, y_integral, eta_a, eta_b, int(order)),
        }
        for order in ETA_INTEGRAL_ORDERS
    ]

    assert finest_estimate is not None
    sign_flip_response = _metric(finest_estimate, -ref_p_eta)["relative_rms"]

    guards: dict[str, bool] = {
        "candidate_identity_bound": primitive.sha256 == AGENT1_CANDIDATE_SHA256,
        "p_eta_finest_relative_rms": p_eta_levels[-1]["relative_rms"]
        <= LOCAL_GATES["p_eta_relative_rms"],
        "p_eta_finest_relative_max": p_eta_levels[-1]["relative_max"]
        <= LOCAL_GATES["p_eta_relative_max"],
        "p_eta_refinement": min(p_eta_ratios)
        >= LOCAL_GATES["p_eta_min_refinement_ratio"],
        "eta_fundamental_theorem_finest_relative_rms": integral_levels[-1][
            "relative_rms"
        ]
        <= LOCAL_GATES["eta_fundamental_theorem_relative_rms"],
        "eta_fundamental_theorem_finest_relative_max": integral_levels[-1][
            "relative_max"
        ]
        <= LOCAL_GATES["eta_fundamental_theorem_relative_max"],
        "p_eta_nontrivial": _rms(ref_p_eta)
        >= LOCAL_GATES["nontrivial_derivative_rms_min"],
        "p_eta_sign_flip_detected": sign_flip_response
        >= LOCAL_GATES["mutation_relative_rms_min"],
    }
    failed_guards = sorted(name for name, passed in guards.items() if not passed)

    truth_boundary = {
        "selected_center_public_p_eta_value_path_assessed": True,
        "source_pressure_radius_one_ball_bound_assessed": False,
        "source_pressure_radius_one_ball_lipschitz_assessed": False,
        "source_mixed_Y_Phi_Y_radius_one_ball_bound_assessed": False,
        "source_R1_R2_bounds_assessed": False,
        "source_MK_assessed": False,
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
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }

    return {
        "schema": SCHEMA,
        "agent1_pr": AGENT1_PR,
        "agent1_head": AGENT1_HEAD,
        "candidate_sha256": primitive.sha256,
        "independence_contract": {
            "public_api_only": "KokunoPA10SelectedPressurePrimitive.evaluate",
            "candidate_local_phase_quadrature_used_as_oracle": False,
            "candidate_stationary_factorization_used_as_oracle": False,
            "candidate_internal_radial_integrals_used_as_oracle": False,
            "training_or_construction_tensor_read": False,
            "validation_derivative": "centered FD6 on public p values",
            "global_identity": "independent Gauss-Legendre integral of public p_eta",
        },
        "frozen_protocol": {
            "seed": SEED,
            "sample_count": SAMPLE_COUNT,
            "selected_lambda": SELECTED_LAMBDA,
            "phase_stationary_eta": PHASE_STATIONARY_ETA,
            "lambda_minus_half_width": width,
            "p_eta_width_fractions": list(P_ETA_WIDTH_FRACTIONS),
            "eta_integral_orders": list(ETA_INTEGRAL_ORDERS),
            "local_gates": dict(LOCAL_GATES),
            "final_project_gates_unchanged": dict(FINAL_PROJECT_GATES),
        },
        "p_eta_public_reference": {
            "rms": _rms(ref_p_eta),
            "sampled_max": float(np.max(np.abs(ref_p_eta))),
        },
        "p_eta_fd6_levels": p_eta_levels,
        "p_eta_fd6_refinement_ratios": p_eta_ratios,
        "eta_fundamental_theorem_levels": integral_levels,
        "p_eta_sign_flip_relative_rms": sign_flip_response,
        "guards": guards,
        "failed_guards": failed_guards,
        "repaired_public_p_eta_independent_preflight_passed": not failed_guards,
        "truth_boundary": truth_boundary,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = run_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
