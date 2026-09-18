import pytest

from openai_ns_reconstruction.kokuno_agent4_source_phase_independent_audit import (
    FROZEN_GUARDS,
    STEPS,
    run_audit,
)


@pytest.fixture(scope="module")
def report():
    return run_audit()


def test_frozen_independent_phase_audit_passes(report):
    assert report["frozen_steps"] == list(STEPS)
    assert report["frozen_guards"] == FROZEN_GUARDS
    assert all(report["guards"].values())
    assert report["local_structural_preflight_passed"] is True


def test_audit_is_value_only_and_separate_from_agent2_derivatives(report):
    independence = report["independence"]
    assert independence["public_values_only"] is True
    assert independence["agent2_fd_helper_reused"] is False
    assert independence["complete_curl_helper_reused"] is False
    assert independence["training_tensor_or_loss_used"] is False
    assert independence["pressure_or_forcing_fit_used"] is False
    assert independence["source_schedule_class_reused_by_reference"] is False


def test_heldout_coverage_and_mutations_are_nontrivial(report):
    metrics = report["metrics"]
    assert metrics["sign_case_count"] == 72
    assert metrics["axis_near_sign_case_count"] >= 18
    assert metrics["off_grid_sign_case_count"] >= 54
    assert metrics["positive_coarse_radial_stencil_sign_case_count"] == 72
    assert set(metrics["carrier_k_values_observed"]) >= {2, 3, 4}
    assert metrics["wrong_tilt_mutation_relative_rms"] >= FROZEN_GUARDS[
        "wrong_tilt_mutation_relative_rms_min"
    ]
    assert metrics["wrong_axial_normalization_relative_rms"] >= FROZEN_GUARDS[
        "wrong_axial_normalization_relative_rms_min"
    ]


def test_three_level_fd6_refines_and_truth_boundary_stays_closed(report):
    metrics = report["metrics"]
    assert len(metrics["fd_relative_rms_by_step"]) == 3
    assert min(metrics["fd_refinement_ratios"]) >= FROZEN_GUARDS[
        "fd_min_refinement_ratio"
    ]
    assert metrics["fd_relative_rms_by_step"][str(STEPS[-1])] <= FROZEN_GUARDS[
        "fd_finest_relative_rms"
    ]

    truth = report["truth_boundary"]
    assert truth["displayed_source_phase_frame_independently_audited"] is True
    for key in (
        "actual_positive_order_background_bound",
        "actual_source_h_sigma_pulse_integrals_bound",
        "actual_signed_auxiliary_rectangles_bound",
        "actual_auxiliary_torus_mode_family_bound",
        "source_actual_partition_labels_instantiated",
        "public_source_bound_velocity_osc_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        assert truth[key] is False

    gates = report["formal_project_gates_unchanged"]
    assert gates["normalized_momentum_max_l2"] == 1.0e-3
    assert gates["divergence_max_l2"] == 1.0e-5
    assert gates["assessed_here"] is False
