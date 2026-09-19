"""Materialize the finite-head mean debt of the admitted oscillatory component.

Kokuno Agent 3 already has two separately verified ingredients:

* PR #581 materializes a real candidate-generated oscillatory self-defect mean
  and its compact radial requested-stress profiles from the independently
  admitted Agent-2 public velocity; and
* PR #514 defines an explicit repository-autonomous finite-head partition
  factor.  That factor is deliberately *not* the pinned formal theorem's
  nonconstructive Mathlib ``missingWeight`` realization.

This module performs exactly one new operation: it combines those frozen
ingredients through

    Delta C_osc,aut(r) = -w_aut * requestedStress_osc(r).

No caller may supply a defect, stress, target, or gain.  The source stress is
recomputed internally from the admitted public candidate using PR #581's fixed
operator and sampling contract, while ``w_aut`` is recomputed from the exact
PR #514 autonomous profile and frozen parameters ``N=5, n=5, q=1.3``.

The result is only an *oscillatory-component* finite-head mean debt.  Agent-1
leading/cross terms, matched pressure and restricted forcing are not yet in the
same-cycle composite defect, so this module intentionally does not open the
signed mean inverse or finite correction cycle.
"""
from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_actual_oscillatory_mean_stress import (
    ADMITTED_AGENT2_HEAD,
    ADMITTED_AGENT2_PR,
    INDEPENDENT_AGENT4_HEAD,
    INDEPENDENT_AGENT4_PR,
    materialize_actual_oscillatory_mean_stress,
)

TASK = "KOKUNO-A3-AUTONOMOUS-OSCILLATORY-MEAN-DEBT-040"
SCHEMA = "kokuno-a3-autonomous-oscillatory-mean-debt-v1"

PARENT_AGENT3_STRESS_PR = 581
PARENT_AGENT3_STRESS_HEAD = "c95dab915ae25e689e7c82bf1b8614e205f9274c"
AUTONOMOUS_FACTOR_ORIGIN_PR = 514
AUTONOMOUS_FACTOR_ORIGIN_HEAD = "50d1b8ca965b35d5e0a7210c22e2302cd30c20e8"
AUTONOMOUS_PROFILE_VERSION = "agent3-cinf-logdyadic-l2-v1"
AUTONOMOUS_PROFILE_FORMULA = (
    "b(s)=exp(1-1/(1-s^2)) for |s|<1 else 0; "
    "raw_m(q)=b(log2(q)+m); mask_m=raw_m/sqrt(sum_j raw_j^2)"
)
AUTONOMOUS_PROFILE_SUPPORT = "|log2(q)+m|<1"
PREPARED_N = 5
BAND = 5
COORDINATE_Q = 1.3
EXPECTED_AUTONOMOUS_MISSING_WEIGHT = 0.28410624176179694

# Freeze the same bounded-cost actual-candidate stress receipt parameters used by
# PR #581's dedicated CI.  These are evaluation coordinates/resolution, not
# tunable residual targets.
RADIAL_COUNT = 25
ANGULAR_COUNT = 32
TIME = 0.50
Z = 0.08
FD4_STEP = 0.005
NU = 0.01


def _autonomous_profile(s: float) -> float:
    value = float(s)
    if not math.isfinite(value):
        raise ValueError("profile coordinate must be finite")
    if abs(value) >= 1.0:
        return 0.0
    return math.exp(1.0 - 1.0 / (1.0 - value * value))


def _replay_autonomous_factor() -> dict[str, Any]:
    """Recompute, rather than hand-enter, the exact autonomous PR #514 factor."""
    q_phys = math.ldexp(COORDINATE_Q, -BAND)
    log_q = math.log2(q_phys)
    center = -log_q
    lo = math.floor(center) - 2
    hi = math.ceil(center) + 2
    raw = {
        m: _autonomous_profile(log_q + m)
        for m in range(lo, hi + 1)
    }
    raw = {m: value for m, value in raw.items() if value > 0.0}
    norm_sq = math.fsum(value * value for value in raw.values())
    if not math.isfinite(norm_sq) or norm_sq <= 0.0:
        raise ArithmeticError("autonomous partition failed to cover physical q")
    norm = math.sqrt(norm_sq)
    masks = {m: value / norm for m, value in raw.items()}
    closure = math.fsum(value * value for value in masks.values())
    if abs(closure - 1.0) > 128.0 * math.ulp(1.0):
        raise ArithmeticError("autonomous squared partition lost normalization")
    finite_head = set(range(-1, PREPARED_N))
    weight = math.fsum(value * value for m, value in masks.items() if m in finite_head)
    if abs(weight - EXPECTED_AUTONOMOUS_MISSING_WEIGHT) > 5.0e-15:
        raise RuntimeError("autonomous PR #514 factor replay drifted")
    profile_payload = json.dumps(
        {
            "version": AUTONOMOUS_PROFILE_VERSION,
            "formula": AUTONOMOUS_PROFILE_FORMULA,
            "support": AUTONOMOUS_PROFILE_SUPPORT,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return {
        "origin_pr": AUTONOMOUS_FACTOR_ORIGIN_PR,
        "origin_head": AUTONOMOUS_FACTOR_ORIGIN_HEAD,
        "profile_version": AUTONOMOUS_PROFILE_VERSION,
        "profile_sha256": hashlib.sha256(profile_payload).hexdigest(),
        "profile_formula": AUTONOMOUS_PROFILE_FORMULA,
        "profile_support": AUTONOMOUS_PROFILE_SUPPORT,
        "prepared_N": PREPARED_N,
        "band": BAND,
        "coordinate_q": COORDINATE_Q,
        "physical_q": q_phys,
        "finite_head_indices": tuple(range(-1, PREPARED_N)),
        "active_mask_indices": tuple(sorted(masks)),
        "active_mask_values": {str(m): masks[m] for m in sorted(masks)},
        "squared_partition_sum": closure,
        "autonomous_missing_weight": weight,
        "profile_is_repository_autonomous": True,
        "formal_theorem_machine_bump_identity_claimed": False,
        "formal_missing_weight_equality_claimed": False,
        "theorem_missing_weight_replaced": False,
        "theorem_missing_weight_materialized": False,
    }


def materialize_autonomous_oscillatory_mean_debt() -> dict[str, Any]:
    """Build the frozen autonomous finite-head debt from the real oscillatory component."""
    stress_signature = inspect.signature(materialize_actual_oscillatory_mean_stress)
    forbidden = {"defect", "residual", "target", "stress", "requested_stress", "gain"}
    exposed = forbidden.intersection(stress_signature.parameters)
    if exposed:
        raise RuntimeError(
            "upstream actual-stress API unexpectedly accepts caller target data: "
            + ", ".join(sorted(exposed))
        )

    stress_report = materialize_actual_oscillatory_mean_stress(
        radial_count=RADIAL_COUNT,
        angular_count=ANGULAR_COUNT,
        time=TIME,
        z=Z,
        fd4_step=FD4_STEP,
        nu=NU,
    )
    truth = stress_report["truth_boundary"]
    anti = stress_report["anti_surrogate_contract"]
    if stress_report["provenance"]["admitted_agent2_head"] != ADMITTED_AGENT2_HEAD:
        raise RuntimeError("actual oscillatory stress drifted from admitted Agent-2 identity")
    if stress_report["provenance"]["independent_agent4_head"] != INDEPENDENT_AGENT4_HEAD:
        raise RuntimeError("actual oscillatory stress drifted from independent Agent-4 identity")
    if anti["surrogate_defect_used"] is not False:
        raise RuntimeError("surrogate defect is forbidden")
    if anti["caller_supplied_target_parameters"] != []:
        raise RuntimeError("caller-supplied target parameters are forbidden")
    required_true = (
        "admitted_agent2_velocity_consumed",
        "exact_candidate_time_derivative_consumed",
        "real_oscillatory_self_defect_component_consumed",
        "real_oscillatory_quadratic_mean_component_measured",
        "oscillatory_requested_stress_component_materialized",
    )
    if not all(truth[name] is True for name in required_true):
        raise RuntimeError("upstream oscillatory mean-stress truth boundary is not admitted")
    if truth["full_same_cycle_composite_requested_stress_materialized"] is not False:
        raise RuntimeError("this component-only lane must not masquerade as the full composite defect")

    stress = stress_report["requested_stress_component"]
    if stress["ordering"] != ["theta_e2", "axial_e1"]:
        raise RuntimeError("requested-stress channel ordering drifted")
    theta = np.asarray(stress["theta_stress_values"], dtype=float)
    axial = np.asarray(stress["axial_stress_values"], dtype=float)
    radii = np.asarray(stress["radial_nodes"], dtype=float)
    if theta.ndim != 1 or axial.shape != theta.shape or radii.shape != theta.shape:
        raise RuntimeError("requested-stress profile shape contract changed")
    if not np.all(np.isfinite(theta)) or not np.all(np.isfinite(axial)):
        raise RuntimeError("requested-stress profile contains non-finite values")

    factor = _replay_autonomous_factor()
    weight = float(factor["autonomous_missing_weight"])
    debt_theta = -weight * theta
    debt_axial = -weight * axial
    debt = np.stack((debt_theta, debt_axial), axis=-1)
    stress_matrix = np.stack((theta, axial), axis=-1)
    relation_error = float(np.max(np.abs(debt + weight * stress_matrix)))
    debt_rms = float(np.sqrt(np.mean(np.sum(debt * debt, axis=-1))))
    expected_rms = weight * float(stress["requested_stress_rms"])
    if relation_error > 64.0 * math.ulp(max(1.0, float(np.max(np.abs(debt))))):
        raise ArithmeticError("finite-head debt identity lost floating-point closure")
    if not math.isclose(debt_rms, expected_rms, rel_tol=2.0e-15, abs_tol=1.0e-15):
        raise ArithmeticError("finite-head debt RMS lost scalar-weight consistency")

    return {
        "task": TASK,
        "schema": SCHEMA,
        "provenance": {
            "parent_agent3_stress_pr": PARENT_AGENT3_STRESS_PR,
            "parent_agent3_stress_head": PARENT_AGENT3_STRESS_HEAD,
            "admitted_agent2_pr": ADMITTED_AGENT2_PR,
            "admitted_agent2_head": ADMITTED_AGENT2_HEAD,
            "independent_agent4_pr": INDEPENDENT_AGENT4_PR,
            "independent_agent4_head": INDEPENDENT_AGENT4_HEAD,
            "autonomous_factor_origin_pr": AUTONOMOUS_FACTOR_ORIGIN_PR,
            "autonomous_factor_origin_head": AUTONOMOUS_FACTOR_ORIGIN_HEAD,
        },
        "evaluation_state": {
            "radial_count": RADIAL_COUNT,
            "angular_count": ANGULAR_COUNT,
            "time": TIME,
            "z": Z,
            "fd4_step": FD4_STEP,
            "nu": NU,
        },
        "autonomous_factor": factor,
        "oscillatory_requested_stress": {
            "ordering": stress["ordering"],
            "requested_stress_rms": float(stress["requested_stress_rms"]),
            "requested_stress_max_abs": float(stress["requested_stress_max_abs"]),
            "radial_nodes": radii.tolist(),
            "theta_e2_values": theta.tolist(),
            "axial_e1_values": axial.tolist(),
        },
        "oscillatory_component_mean_debt": {
            "formula": "DeltaC_osc_aut=-autonomous_missing_weight*requestedStress_osc",
            "ordering": ["theta_e2", "axial_e1"],
            "theta_e2_values": debt_theta.tolist(),
            "axial_e1_values": debt_axial.tolist(),
            "debt_rms": debt_rms,
            "debt_max_abs": float(np.max(np.abs(debt))),
            "scalar_weight_rms_prediction": expected_rms,
            "identity_closure_max_abs": relation_error,
        },
        "anti_surrogate_contract": {
            "caller_supplied_defect_or_stress_parameters": [],
            "surrogate_defect_used": False,
            "actual_oscillatory_self_defect_recomputed_internally": True,
            "autonomous_factor_frozen_before_this_debt_evaluation": True,
            "threshold_or_gain_tuned_after_observation": False,
        },
        "truth_boundary": {
            "real_oscillatory_self_defect_component_consumed": True,
            "oscillatory_requested_stress_component_materialized": True,
            "autonomous_finite_head_factor_applied": True,
            "oscillatory_component_finite_head_mean_debt_materialized": True,
            "formal_missing_weight_materialized": False,
            "formal_missing_weight_equality_claimed": False,
            "full_same_cycle_composite_requested_stress_materialized": False,
            "agent1_leading_cross_terms_included": False,
            "matched_pressure_included": False,
            "restricted_forcing_included": False,
            "candidate_finite_head_mean_debt_materialized": False,
            "real_full_candidate_defect_consumed": False,
            "signed_mean_inverse_input_ready": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_rerun_allowed": False,
            "heldout_ns_residual_assessed": False,
            "residual_reduction_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "blowup_proved": False,
        },
        "blockers": (
            "Agent-1 leading/cross terms are not yet present in the same-cycle composite defect",
            "global matched pressure and restricted-forcing terms are not yet present",
            "repository-autonomous finite-head factor is not the formal theorem-machine missingWeight",
            "signed mean inverse must remain closed until a full same-cycle target is materialized",
        ),
    }


def write_receipt(output: str | Path) -> dict[str, Any]:
    report = materialize_autonomous_oscillatory_mean_debt()
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3_autonomous_oscillatory_mean_debt_v1.json",
    )
    args = parser.parse_args()
    print(json.dumps(write_receipt(args.output), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
