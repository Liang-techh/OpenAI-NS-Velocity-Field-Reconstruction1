from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_independent_oscillatory_admission import (
    EXPECTED_AGENT2_HEAD,
    FROZEN_GUARDS,
    evaluate_independent_oscillatory_admission,
)


def _rejected_agent4_audit() -> dict:
    """Decision-relevant values from Agent-4 PR #531 exact-head artifact."""

    return {
        "schema": "kokuno-agent4-public-oscillatory-independent-audit-v1",
        "parent_agent2_head": EXPECTED_AGENT2_HEAD,
        "seed": 9173241,
        "protocol": {
            "primary_interface": "module-level velocity_osc(x,y,z,t) only",
            "heldout_offgrid_points": 24,
            "fd4_steps": [0.02, 0.01, 0.005],
            "phase_resolutions": [32, 64, 128],
            "pressure_or_forcing_fit": False,
            "candidate_internal_derivative_oracle_used": False,
            "training_tensor_or_loss_used": False,
            "project_registered_support": "r<2 and |z|<2, smooth zero extension",
        },
        "guards_frozen_before_actions": dict(FROZEN_GUARDS),
        "divergence": {
            "rows": [
                {"step": 0.02, "relative_rms": 0.6916422603413279, "relative_max": 0.8442456766582342},
                {"step": 0.01, "relative_rms": 0.4329981107094264, "relative_max": 0.6071654826117929},
                {"step": 0.005, "relative_rms": 0.10460795155383176, "relative_max": 0.09811443995058577},
            ],
            "relative_rms_refinement_ratios": [1.5973332059304777, 4.139246627791994],
            "local_divergence_passed": False,
        },
        "divergence_mutation": {
            "relative_rms": 0.23755538605232182,
            "caught": True,
        },
        "phase_mean_covariance": {
            "rows": [
                {"phase_resolution": 32, "rank_ratio": 0.018748193463796137},
                {"phase_resolution": 64, "rank_ratio": 0.01874819346346177},
                {"phase_resolution": 128, "rank_ratio": 0.018748193463887446},
            ],
            "relative_drift_to_finest": [5.718781655269777e-12, 2.060561561211433e-11],
            "duplicated_first_column_rank_ratio": 9.669379904398983e-19,
            "covariance_rank_preflight_passed": False,
        },
        "support": {
            "axis_near_absolute_max": 0.0,
            "radial_exterior_absolute_max": 0.0,
            "project_axial_exterior_absolute_max": 769255.3191046526,
            "radial_axis_support_passed": True,
            "project_axial_support_passed": False,
        },
        "nontriviality": {"heldout_velocity_rms": 317659.12117004656, "passed": True},
        "parameter_perturbation": {
            "h_minus_10pct_relative_change": 1.2614067667595157,
            "h_plus_10pct_relative_change": 0.009101848910306524,
            "passed": True,
        },
        "local_curl_covariance_preflight_passed": False,
        "project_support_preflight_passed": False,
        "public_oscillatory_preflight_passed": False,
        "truth_boundary": {
            "oscillatory_component_only": True,
            "global_leading_velocity_pressure_available": False,
            "materialized_correction_available": False,
            "formal_full_domain_pde_gate_assessed": False,
            "heldout_ns_residual_assessed": False,
            "pde_validated": False,
            "formal_momentum_normalized_max_l2_gate": 1.0e-3,
            "formal_divergence_max_l2_gate": 1.0e-5,
        },
    }


def _passing_audit() -> dict:
    audit = copy.deepcopy(_rejected_agent4_audit())
    audit["divergence"]["rows"][-1]["relative_rms"] = 1.0e-5
    audit["divergence"]["rows"][-1]["relative_max"] = 5.0e-5
    audit["divergence"]["relative_rms_refinement_ratios"] = [4.0, 4.0]
    audit["divergence"]["local_divergence_passed"] = True
    for row in audit["phase_mean_covariance"]["rows"]:
        row["rank_ratio"] = 0.03
    audit["phase_mean_covariance"]["covariance_rank_preflight_passed"] = True
    audit["support"]["project_axial_exterior_absolute_max"] = 0.0
    audit["support"]["project_axial_support_passed"] = True
    audit["local_curl_covariance_preflight_passed"] = True
    audit["project_support_preflight_passed"] = True
    audit["public_oscillatory_preflight_passed"] = True
    return audit


def test_current_independent_rejection_blocks_agent3_mean_debt() -> None:
    report = evaluate_independent_oscillatory_admission(_rejected_agent4_audit())
    assert report["independent_public_oscillatory_preflight_passed"] is False
    assert report["failed_guards"] == (
        "finest_relative_divergence_rms",
        "finest_relative_divergence_max",
        "minimum_divergence_refinement_ratio",
        "minimum_covariance_rank_ratio",
        "project_axial_exterior_absolute_max",
    )
    assert report["correction_ingest_allowed"] is False
    assert report["same_cycle_requested_stress_materialization_allowed"] is False
    assert report["candidate_finite_head_mean_debt_materialization_allowed"] is False
    assert report["finite_correction_cycle_rerun_allowed"] is False
    assert report["surrogate_defect_used"] is False
    assert report["pde_validated"] is False


def test_future_independent_pass_only_opens_materialization_stage() -> None:
    report = evaluate_independent_oscillatory_admission(_passing_audit())
    assert report["failed_guards"] == ()
    assert report["independent_public_oscillatory_preflight_passed"] is True
    assert report["same_cycle_requested_stress_materialization_allowed"] is True
    assert report["candidate_finite_head_mean_debt_materialization_allowed"] is True
    assert report["requested_stress_actual_state_values_materialized"] is False
    assert report["finite_head_mean_debt_materialized"] is False
    assert report["real_candidate_defect_consumed"] is False
    assert report["signed_mean_inverse_input_ready"] is False
    assert report["finite_correction_cycle_rerun_allowed"] is False
    assert report["heldout_ns_residual_assessed"] is False
    assert report["residual_reduction_claimed"] is False


def test_threshold_laundering_is_rejected() -> None:
    audit = _rejected_agent4_audit()
    audit["guards_frozen_before_actions"]["minimum_covariance_rank_ratio"] = 0.018
    with pytest.raises(ValueError, match="changed from the frozen"):
        evaluate_independent_oscillatory_admission(audit)


def test_reported_pass_cannot_override_measured_failure() -> None:
    audit = _rejected_agent4_audit()
    audit["public_oscillatory_preflight_passed"] = True
    with pytest.raises(ValueError, match="reported public oscillatory verdict"):
        evaluate_independent_oscillatory_admission(audit)


def test_truth_boundary_promotion_is_rejected() -> None:
    audit = _rejected_agent4_audit()
    audit["truth_boundary"]["heldout_ns_residual_assessed"] = True
    with pytest.raises(ValueError, match="truth boundary"):
        evaluate_independent_oscillatory_admission(audit)


def test_wrong_provider_head_is_rejected() -> None:
    audit = _rejected_agent4_audit()
    audit["parent_agent2_head"] = "0" * 40
    with pytest.raises(ValueError, match="provider head"):
        evaluate_independent_oscillatory_admission(audit)
