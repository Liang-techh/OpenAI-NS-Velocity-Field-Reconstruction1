"""Agent-5 checkpoint for the next three Kokuno numeric prerequisite seams.

This is one minimal integration/routing increment on top of the v34 Agent-5
checkpoint.  It records three newly executable sibling deliveries without
copying their mathematics into the integration lane or promoting a local
engineering quantity into source truth:

* Agent 1 exposes the strict pointwise PA.10 *necessary* admissible band for a
  newly propagated coupled normalization.  The old selected shared-C path is
  still excluded; passing the new screen would only mean "not excluded".
* Agent 2 can numerically materialize the source-displayed signed covariance
  mass law h_sigma on explicitly labelled candidate pulse/cutoff samples.  The
  numerical samples remain repository/candidate choices, not recovered hidden
  Kokuno data, and no public physical oscillatory velocity is produced yet.
* Agent 3 supplies a deterministic repository-autonomous finite-head mean
  factor for engineering closure.  It is deliberately distinct from the
  nonconstructively selected theorem-machine missingWeight and cannot become a
  physical Delta C until same-cycle requestedStress is materialized.

The result is a narrower shortest path, not a 3D candidate and not an NS gate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from .kokuno_audited_geometry_leading_routing_checkpoint import (
    AUTONOMOUS_REALIZATION_POLICY,
    FORMAL_GATES,
    SCHEMA as PARENT_SCHEMA,
    TASK_ID as PARENT_TASK_ID,
    build_checkpoint as build_parent_checkpoint,
    checkpoint_sha256,
    validate_checkpoint as validate_parent_checkpoint,
)

SCHEMA = "kokuno-agent5-numeric-prerequisites-routing-checkpoint-v35"
TASK_ID = "KOKUNO-A5-NUMERIC-PREREQUISITES-ROUTING-035"

AGENT1_PA10_SCREEN_RECEIPT = {
    "pr": 512,
    "head": "864263ac2c9c1a030105dcf209a71b5602976d23",
    "evidence_class": "candidate_derived_necessary_screen",
    "dedicated_run": 35405746620,
    "standard_run": 35405782809,
    "dedicated_status": "success",
    "standard_status": "success",
    "artifact": "kokuno-agent1-coupled-normalization-pa10-screen-v1",
    "artifact_id": 10571888989,
    "artifact_digest": "sha256:0c943566c861d89e41a0ec878d1d802fdf2a20a4c03d18a93a559706856bf0da",
    "coupled_shared_C_dependency_recorded": True,
    "pointwise_PA10_necessary_screen_executable": True,
    "required_log_E_i_open_band_executable": True,
    "selected_realization_excluded": True,
    "posthoc_C_multiplier_applied": False,
    "passing_pointwise_screen_is_source_T_sh_certificate": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": False,
}

AGENT2_SIGNED_MASS_RECEIPT = {
    "pr": 513,
    "head": "47f918e482981e7b124b48cef1781188302c4446",
    "evidence_class": "candidate_derived_source_formula_evaluation",
    "dedicated_run": 35406299226,
    "standard_run": 35406299071,
    "dedicated_status": "success",
    "standard_status": "success",
    "focused_passed": 19,
    "focused_seconds": 0.84,
    "source_h_sigma_integral_formula_executable": True,
    "source_displayed_D_g_formula_executable": True,
    "candidate_h_sigma_numerically_materializable": True,
    "candidate_h_sigma_can_feed_reference_covariance_pair": True,
    "actual_source_h_sigma_pulse_integrals_bound": False,
    "actual_source_pulse_samples_recovered": False,
    "actual_positive_order_background_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "public_provenance_labelled_xyz_t_oscillatory_velocity_ready": False,
    "actual_source_physical_covariance_rank_two_assessed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": False,
}

AGENT3_AUTONOMOUS_MEAN_RECEIPT = {
    "pr": 514,
    "head": "50d1b8ca965b35d5e0a7210c22e2302cd30c20e8",
    "evidence_class": "repository_autonomous_engineering_factor",
    "dedicated_run": 35406792841,
    "standard_run": 35406792929,
    "dedicated_status": "success",
    "standard_status": "success",
    "artifact": "kokuno-agent3-autonomous-mean-factor",
    "artifact_id": 10572741978,
    "artifact_digest": "sha256:a9bb1e4ffad4a7be0c9aebb70ef4d53bded3bdd064cf1c21dc4c49f4aa6fdc73",
    "prepared_N": 5,
    "band": 5,
    "coordinate_q": 1.3,
    "physical_q": 0.040625,
    "active_mask_indices": [4, 5],
    "autonomous_missing_weight": 0.28410624176179694,
    "autonomous_finite_head_mean_factor_executable": True,
    "formal_theorem_machine_bump_identity_claimed": False,
    "formal_missing_weight_equality_claimed": False,
    "theorem_missing_weight_replaced": False,
    "theorem_missing_weight_materialized": False,
    "requested_stress_actual_state_values_materialized": False,
    "finite_head_mean_debt_materialized": False,
    "real_candidate_defect_consumed": False,
    "signed_mean_inverse_input_ready": False,
    "finite_correction_cycle_rerun_allowed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": False,
}

PREVIOUS_AGENT5_RECEIPT = {
    "pr": 511,
    "head": "d244804a4b2322127f42239803d7023557be2325",
    "schema": "kokuno-agent5-audited-geometry-leading-routing-checkpoint-v34",
    "dedicated_run": 35405175156,
    "standard_run": 35405175254,
    "dedicated_status": "success",
    "standard_status": "success",
    "artifact": "kokuno-agent5-audited-geometry-leading-routing-checkpoint-v34",
    "artifact_id": 10571992884,
    "artifact_digest": "sha256:80b5a1b3d5587aace739cefacd4b46984e6c0db4e7b0947ed61658b9f038f376",
    "consumed_in_executable_ancestry": True,
}


def _copy(value: Any) -> Any:
    return json.loads(json.dumps(value))


def _set_pipeline_status(payload: dict[str, Any], stage: str, status: str) -> None:
    for item in payload["pipeline"]:
        if item["stage"] == stage:
            item["ready"] = False
            item["status"] = status
            return
    raise ValueError(f"missing pipeline stage: {stage}")


def build_checkpoint() -> dict[str, Any]:
    payload = _copy(build_parent_checkpoint())
    payload.pop("checkpoint_sha256", None)
    payload["schema"] = SCHEMA
    payload["task_id"] = TASK_ID

    upstream = payload["upstream"]
    upstream["agent1_pa10_screen_sibling"] = _copy(AGENT1_PA10_SCREEN_RECEIPT)
    upstream["agent2_signed_mass_sibling"] = _copy(AGENT2_SIGNED_MASS_RECEIPT)
    upstream["agent3_autonomous_mean_sibling"] = _copy(AGENT3_AUTONOMOUS_MEAN_RECEIPT)
    upstream["previous_agent5_v34_ancestry"] = _copy(PREVIOUS_AGENT5_RECEIPT)

    states = payload["states"]
    states["leading_pa10_required_band_executable"] = True
    states["candidate_signed_h_sigma_mass_executable"] = True
    states["autonomous_finite_head_mean_factor_executable"] = True
    states["same_cycle_requested_stress_materialized"] = False

    _set_pipeline_status(
        payload,
        "leading_candidate",
        "Agent 1 now exposes the strict coupled-normalization PA.10 necessary target band; the previously selected shared-C realization remains excluded, so propagate a fresh provenance-labelled upstream Appendix-B normalization into that band and rerun the same screen before any PA.16/I3/I4/matched-pressure handoff",
    )
    _set_pipeline_status(
        payload,
        "oscillatory_augmentation",
        "Agent 2 now evaluates the displayed h_sigma mass law on labelled numerical candidate pulses over the independently audited autonomous geometry; next bind provenance-labelled positive-order/background and auxiliary modes into the audited phase/complete-curl path and emit Q-scaled by-sign/by-beta/total velocity_osc(x,y,z,t), without promoting candidate pulse samples to recovered source data",
    )
    _set_pipeline_status(
        payload,
        "mean_radial_corrections",
        "Agent 3 now has a deterministic autonomous finite-head mean factor, explicitly not theorem missingWeight; next materialize same-cycle physical requestedStress from the candidate state, form only a candidate-specific Delta C, then apply the existing signed inverse/budget/spacetime/radial/quadratic guards after the physical oscillatory covariance interface is independently audited",
    )
    _set_pipeline_status(
        payload,
        "independent_validation",
        "Agent 4 has already independently closed the autonomous rectangle geometry seam; its next highest-value target remains the first fully numerical provenance-labelled public physical oscillatory velocity and phase-mean covariance response, not the local h_sigma quadrature or autonomous mean-factor bookkeeping",
    )

    payload["routing"] = {
        "new_fact": (
            "Three upstream numerical prerequisites are now executable with explicit truth boundaries: "
            "a necessary PA.10 leading target band, candidate signed h_sigma mass evaluation, and a "
            "repository-autonomous finite-head mean factor. None is yet the missing physical candidate seam."
        ),
        "shortest_next_closure": [
            "Agent 1: propagate one fresh provenance-labelled normalization through Appendix B into the #512 required log E_i band and rerun the frozen PA.10 screen; no post-hoc C retuning and no PA.16 handoff merely because a pointwise screen passes.",
            "Agent 2: stop at no more standalone h_sigma arithmetic; bind numerical positive-order/background plus signed auxiliary modes to the audited phase/complete-curl stack and expose Q-scaled by-sign/by-beta/total velocity_osc(x,y,z,t) with explicit source-vs-autonomous provenance.",
            "Agent 3: compute theta/axial requestedStress from the same candidate cycle state and combine it with the autonomous factor only as a candidate-specific finite-head mean debt; do not name it formal missingWeight or run the finite correction cycle before the physical covariance interface passes independent audit.",
            "Agent 4: black-box audit the first public fully numerical physical oscillatory velocity, complete-curl/divergence behavior and phase-mean covariance rank using only the public velocity interface and frozen guards.",
            "Agent 5: when a non-obstructed leading candidate, independently audited public oscillatory velocity and guarded same-cycle correction coexist, instantiate the deterministic candidate artifact with velocity/pressure/restricted-forcing APIs, save/load plus Python/MATLAB smoke, and freeze it for held-out NS validation.",
        ],
    }

    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    return payload


def validate_checkpoint(payload: dict[str, Any]) -> None:
    if payload.get("schema") != SCHEMA or payload.get("task_id") != TASK_ID:
        raise ValueError("checkpoint schema/task mismatch")
    if payload.get("checkpoint_sha256") != checkpoint_sha256(payload):
        raise ValueError("checkpoint sha256 mismatch")
    if payload.get("formal_gates") != FORMAL_GATES:
        raise ValueError("formal PDE gates changed")
    if payload.get("autonomous_realization_policy") != AUTONOMOUS_REALIZATION_POLICY:
        raise ValueError("autonomous realization policy changed")

    # Reuse every v34 fail-closed check.  New sibling receipts are additive and
    # intentionally do not overwrite v34's frozen Agent-1/4 evidence.
    parent_view = _copy(payload)
    parent_view["schema"] = PARENT_SCHEMA
    parent_view["task_id"] = PARENT_TASK_ID
    parent_view["checkpoint_sha256"] = checkpoint_sha256(parent_view)
    validate_parent_checkpoint(parent_view)

    states = payload["states"]
    if states.get("leading_pa10_required_band_executable") is not True:
        raise ValueError("leading PA.10 target band was removed")
    if states.get("candidate_signed_h_sigma_mass_executable") is not True:
        raise ValueError("candidate signed h_sigma mass seam was removed")
    if states.get("autonomous_finite_head_mean_factor_executable") is not True:
        raise ValueError("autonomous finite-head mean factor was removed")
    if states.get("same_cycle_requested_stress_materialized") is not False:
        raise ValueError("same-cycle requestedStress was promoted without evidence")

    # These local prerequisites are deliberately insufficient for any readiness
    # or PDE gate promotion.
    for name in (
        "leading_ready",
        "public_provenance_labelled_xyz_t_oscillatory_velocity_ready",
        "oscillatory_ready",
        "genuinely_independent_second_covariance_column_ready",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        if states.get(name) is not False:
            raise ValueError(f"local numeric prerequisite incorrectly promoted {name}")

    a1 = payload["upstream"].get("agent1_pa10_screen_sibling")
    if a1 != AGENT1_PA10_SCREEN_RECEIPT:
        raise ValueError("Agent 1 PA.10 screen receipt changed")
    if not a1["pointwise_PA10_necessary_screen_executable"] or not a1["selected_realization_excluded"]:
        raise ValueError("Agent 1 necessary-screen semantics changed")
    if a1["passing_pointwise_screen_is_source_T_sh_certificate"] or a1["source_T_sh_lower_bound_verified"]:
        raise ValueError("Agent 1 necessary screen was laundered into source T_sh truth")
    if a1["posthoc_C_multiplier_applied"]:
        raise ValueError("post-hoc leading C retuning was introduced")

    a2 = payload["upstream"].get("agent2_signed_mass_sibling")
    if a2 != AGENT2_SIGNED_MASS_RECEIPT:
        raise ValueError("Agent 2 signed-mass receipt changed")
    if not a2["candidate_h_sigma_numerically_materializable"]:
        raise ValueError("Agent 2 candidate h_sigma seam lost")
    if a2["actual_source_h_sigma_pulse_integrals_bound"] or a2["actual_source_pulse_samples_recovered"]:
        raise ValueError("candidate h_sigma samples were promoted to recovered source data")
    if a2["public_provenance_labelled_xyz_t_oscillatory_velocity_ready"]:
        raise ValueError("local h_sigma mass was promoted to public oscillatory velocity")

    a3 = payload["upstream"].get("agent3_autonomous_mean_sibling")
    if a3 != AGENT3_AUTONOMOUS_MEAN_RECEIPT:
        raise ValueError("Agent 3 autonomous-mean receipt changed")
    if not a3["autonomous_finite_head_mean_factor_executable"]:
        raise ValueError("Agent 3 autonomous mean factor lost")
    if a3["formal_missing_weight_equality_claimed"] or a3["theorem_missing_weight_materialized"]:
        raise ValueError("autonomous factor was promoted to theorem missingWeight")
    if a3["requested_stress_actual_state_values_materialized"] or a3["finite_head_mean_debt_materialized"]:
        raise ValueError("autonomous factor was promoted to same-cycle physical mean debt")
    if a3["finite_correction_cycle_rerun_allowed"]:
        raise ValueError("finite correction cycle was opened before same-cycle defect materialization")

    previous = payload["upstream"].get("previous_agent5_v34_ancestry")
    if previous != PREVIOUS_AGENT5_RECEIPT:
        raise ValueError("previous Agent 5 v34 receipt changed")
    if payload["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is not None:
        raise ValueError("no comparable full-domain Kokuno receipt exists yet")


def write_checkpoint(path: str | Path) -> dict[str, Any]:
    payload = build_checkpoint()
    validate_checkpoint(payload)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = write_checkpoint(args.output)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
