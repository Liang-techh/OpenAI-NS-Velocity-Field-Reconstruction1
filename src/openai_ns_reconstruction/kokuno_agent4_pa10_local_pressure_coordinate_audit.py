"""Independent Agent-4 audit of the PA.10 cancellation-safe pressure coordinate.

This validator consumes only the public
``KokunoPA10LocalPressureCoordinate.evaluate_local(Y, s)`` interface and its
immutable candidate identity.  It deliberately does not call Agent-1 local
phase, stationary-factorization, radial-integral, or self-check helpers.

The audit uses a disjoint held-out cloud and a centered eighth-order finite
-difference derivative in the O(1) local coordinate ``s``.  A separate
composite-Simpson fundamental-theorem check is used instead of Agent 1's
Gauss--Legendre self-check.  The physical chain rule is checked against the
value-only FD8 derivative, and the two-float eta mapping is checked with
``math.fsum``.

This is a selected-center pressure-interface audit only.  It is not a
Navier--Stokes momentum residual evaluation, a source radius-one-ball bound,
or evidence that a global leading field exists.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_pa10_local_pressure_coordinate import KokunoPA10LocalPressureCoordinate


SCHEMA = "kokuno-agent4-pa10-local-pressure-coordinate-independent-audit-v1"
AGENT1_PR = 624
AGENT1_HEAD = "e036b057d9c9fd320e2bb87aeb0fcdb1395f2387"
AGENT1_CANDIDATE_SHA256 = "2c41504eeb448b840280dbb80f40f1209f0efff9d0db2ca0e8826820e73944f6"
PRIOR_ETA_ONLY_AGENT4_PR = 617

SEED = 9173291
SAMPLE_COUNT = 40
FD8_STEPS = (0.08, 0.04, 0.02)
SIMPSON_NODE_COUNTS = (65, 129, 257)

LOCAL_GATES = {
    "p_s_relative_rms": 5.0e-6,
    "p_s_relative_max": 1.0e-5,
    "p_s_min_refinement_ratio": 8.0,
    "local_fundamental_theorem_relative_rms": 1.0e-7,
    "local_fundamental_theorem_relative_max": 5.0e-7,
    "local_fundamental_theorem_min_refinement_ratio": 4.0,
    "physical_chain_rule_relative_rms": 5.0e-6,
    "physical_chain_rule_relative_max": 1.0e-5,
    "two_term_eta_mapping_scaled_max": 5.0e-12,
    "mutation_relative_rms_min": 0.5,
    "nontrivial_p_s_rms_min": 1.0e-8,
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


def _fd8_s(
    interface: KokunoPA10LocalPressureCoordinate,
    y: np.ndarray,
    s: np.ndarray,
    h: float,
) -> np.ndarray:
    """Centered eighth-order first derivative from public pressure values."""

    def values(offset: float) -> np.ndarray:
        return np.asarray(interface.evaluate_local(y, s + offset * h)["p"], dtype=float)

    return (
        (4.0 / 5.0) * (values(1.0) - values(-1.0))
        - (1.0 / 5.0) * (values(2.0) - values(-2.0))
        + (4.0 / 105.0) * (values(3.0) - values(-3.0))
        - (1.0 / 280.0) * (values(4.0) - values(-4.0))
    ) / h


def _simpson_integral(
    interface: KokunoPA10LocalPressureCoordinate,
    y: float,
    left: float,
    right: float,
    node_count: int,
) -> float:
    if node_count < 3 or node_count % 2 == 0:
        raise ValueError("composite Simpson requires an odd node count >= 3")
    points = np.linspace(float(left), float(right), int(node_count), dtype=float)
    values = np.asarray(interface.evaluate_local(float(y), points)["p_s"], dtype=float)
    weights = np.ones(int(node_count), dtype=float)
    weights[1:-1:2] = 4.0
    weights[2:-1:2] = 2.0
    step = (float(right) - float(left)) / float(node_count - 1)
    return float(step * np.dot(weights, values) / 3.0)


def _fundamental_theorem_level(
    interface: KokunoPA10LocalPressureCoordinate,
    y_values: np.ndarray,
    intervals: np.ndarray,
    node_count: int,
) -> dict[str, float]:
    estimates: list[float] = []
    references: list[float] = []
    for y, (left, right) in zip(y_values, intervals, strict=True):
        estimates.append(
            _simpson_integral(interface, float(y), float(left), float(right), node_count)
        )
        p_left = float(interface.evaluate_local(float(y), float(left))["p"])
        p_right = float(interface.evaluate_local(float(y), float(right))["p"])
        references.append(p_right - p_left)
    return _metric(np.asarray(estimates), np.asarray(references))


def _two_term_mapping_metrics(
    interface: KokunoPA10LocalPressureCoordinate,
    s: np.ndarray,
) -> dict[str, float]:
    public = interface.evaluate_local(np.ones_like(s), s)
    hi = np.asarray(public["eta_hi"], dtype=float)
    lo = np.asarray(public["eta_lo"], dtype=float)
    target_delta = np.asarray(public["eta_offset"], dtype=float)
    recovered = np.asarray(
        [
            math.fsum((float(h), float(l), -float(interface.eta_star)))
            for h, l in zip(hi, lo, strict=True)
        ],
        dtype=float,
    )
    error = recovered - target_delta
    scaled = error / float(interface.layer_width)
    return {
        "absolute_rms": _rms(error),
        "absolute_max": float(np.max(np.abs(error))),
        "scaled_rms_in_layer_widths": _rms(scaled),
        "scaled_max_in_layer_widths": float(np.max(np.abs(scaled))),
    }


def run_audit() -> dict[str, Any]:
    interface = KokunoPA10LocalPressureCoordinate(radial_quadrature_points=64)
    if interface.sha256 != AGENT1_CANDIDATE_SHA256:
        raise RuntimeError(
            "Agent-1 local-coordinate identity moved: "
            f"{interface.sha256} != {AGENT1_CANDIDATE_SHA256}"
        )

    rng = np.random.default_rng(SEED)
    y = rng.uniform(0.37, 3.83, size=SAMPLE_COUNT)
    s = rng.uniform(-0.57, 0.57, size=SAMPLE_COUNT)
    public = interface.evaluate_local(y, s)
    ref_p_s = np.asarray(public["p_s"], dtype=float)
    ref_p_eta = np.asarray(public["p_eta"], dtype=float)

    p_s_levels: list[dict[str, float]] = []
    finest_fd: np.ndarray | None = None
    for step in FD8_STEPS:
        estimate = _fd8_s(interface, y, s, float(step))
        p_s_levels.append({"step": float(step), **_metric(estimate, ref_p_s)})
        finest_fd = estimate
    p_s_ratios = _refinement_ratios(p_s_levels)

    intervals = np.asarray(
        [
            (-0.93, -0.17),
            (-0.81, 0.13),
            (-0.66, 0.49),
            (-0.41, 0.88),
            (0.07, 0.79),
            (0.22, 0.96),
        ],
        dtype=float,
    )
    y_integral = np.asarray([0.53, 0.98, 1.57, 2.26, 3.09, 3.74], dtype=float)
    integral_levels = [
        {
            "node_count": int(count),
            **_fundamental_theorem_level(
                interface, y_integral, intervals, int(count)
            ),
        }
        for count in SIMPSON_NODE_COUNTS
    ]
    integral_ratios = _refinement_ratios(integral_levels)

    if finest_fd is None:
        raise RuntimeError("FD8 ladder did not execute")
    chain_estimate = float(interface.sqrt_lambda) * finest_fd
    chain_metric = _metric(chain_estimate, ref_p_eta)
    mapping_metric = _two_term_mapping_metrics(interface, s)
    sign_flip_response = _metric(finest_fd, -ref_p_s)["relative_rms"]

    guards: dict[str, bool] = {
        "candidate_identity_bound": interface.sha256 == AGENT1_CANDIDATE_SHA256,
        "p_s_finest_relative_rms": p_s_levels[-1]["relative_rms"]
        <= LOCAL_GATES["p_s_relative_rms"],
        "p_s_finest_relative_max": p_s_levels[-1]["relative_max"]
        <= LOCAL_GATES["p_s_relative_max"],
        "p_s_refinement": min(p_s_ratios)
        >= LOCAL_GATES["p_s_min_refinement_ratio"],
        "local_fundamental_theorem_finest_relative_rms": integral_levels[-1][
            "relative_rms"
        ]
        <= LOCAL_GATES["local_fundamental_theorem_relative_rms"],
        "local_fundamental_theorem_finest_relative_max": integral_levels[-1][
            "relative_max"
        ]
        <= LOCAL_GATES["local_fundamental_theorem_relative_max"],
        "local_fundamental_theorem_refinement": min(integral_ratios)
        >= LOCAL_GATES["local_fundamental_theorem_min_refinement_ratio"],
        "physical_chain_rule_relative_rms": chain_metric["relative_rms"]
        <= LOCAL_GATES["physical_chain_rule_relative_rms"],
        "physical_chain_rule_relative_max": chain_metric["relative_max"]
        <= LOCAL_GATES["physical_chain_rule_relative_max"],
        "two_term_eta_mapping_scaled_max": mapping_metric[
            "scaled_max_in_layer_widths"
        ]
        <= LOCAL_GATES["two_term_eta_mapping_scaled_max"],
        "p_s_nontrivial": _rms(ref_p_s) >= LOCAL_GATES["nontrivial_p_s_rms_min"],
        "p_s_sign_flip_detected": sign_flip_response
        >= LOCAL_GATES["mutation_relative_rms_min"],
    }
    failed_guards = sorted(name for name, passed in guards.items() if not passed)

    truth_boundary = {
        "prior_eta_only_agent4_pr": PRIOR_ETA_ONLY_AGENT4_PR,
        "prior_eta_only_preflight_passed": False,
        "local_coordinate_interface_assessed": True,
        "source_pressure_radius_one_ball_bound_assessed": False,
        "source_pressure_radius_one_ball_lipschitz_assessed": False,
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
        "candidate_sha256": interface.sha256,
        "independence_contract": {
            "public_api_only": "KokunoPA10LocalPressureCoordinate.evaluate_local",
            "candidate_local_phase_quadrature_used_as_oracle": False,
            "candidate_stationary_factorization_used_as_oracle": False,
            "candidate_radial_integral_helper_used_as_oracle": False,
            "candidate_diagnostic_receipt_used_as_oracle": False,
            "training_or_construction_tensor_read": False,
            "validation_derivative": "centered FD8 on public p values",
            "global_identity": "composite Simpson integral of public p_s",
            "physical_chain_rule_reference": "public p_eta versus sqrt(Lambda)*value-only FD8(p)",
        },
        "frozen_protocol": {
            "seed": SEED,
            "sample_count": SAMPLE_COUNT,
            "fd8_steps": list(FD8_STEPS),
            "simpson_node_counts": list(SIMPSON_NODE_COUNTS),
            "local_gates": dict(LOCAL_GATES),
            "final_project_gates_unchanged": dict(FINAL_PROJECT_GATES),
        },
        "selected": {
            "Lambda": float(interface.core.rescaling_lambda),
            "sqrt_Lambda": float(interface.sqrt_lambda),
            "eta_star": float(interface.eta_star),
            "layer_width": float(interface.layer_width),
        },
        "p_s_public_reference": {
            "rms": _rms(ref_p_s),
            "sampled_max": float(np.max(np.abs(ref_p_s))),
        },
        "p_s_fd8_levels": p_s_levels,
        "p_s_fd8_refinement_ratios": p_s_ratios,
        "local_fundamental_theorem_levels": integral_levels,
        "local_fundamental_theorem_refinement_ratios": integral_ratios,
        "physical_chain_rule_from_value_fd8": chain_metric,
        "two_term_eta_mapping": mapping_metric,
        "p_s_sign_flip_relative_rms": sign_flip_response,
        "guards": guards,
        "failed_guards": failed_guards,
        "local_pressure_coordinate_independent_preflight_passed": not failed_guards,
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
