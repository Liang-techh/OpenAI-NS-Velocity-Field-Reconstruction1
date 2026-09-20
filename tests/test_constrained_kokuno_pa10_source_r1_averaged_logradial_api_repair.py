from __future__ import annotations

import inspect

from openai_ns_reconstruction import (
    kokuno_pa10_source_r1_averaged_logradial_api_repair as repair,
)
from openai_ns_reconstruction.kokuno_pa10_source_r1_averaged_logradial_api_repair import (
    FAILED_A4_HEAD,
    FAILED_A4_RUN,
    run_repaired_independent_audit,
)


def test_repaired_averaged_logradial_audit_passes_frozen_protocol() -> None:
    payload = run_repaired_independent_audit()
    assert payload["failed_guards"] == []
    assert (
        payload[
            "source_R1_averaged_logradial_ordinary_slot_independent_preflight_passed"
        ]
        is True
    )

    assert payload["superseded_failed_a4"]["exact_head"] == FAILED_A4_HEAD
    assert payload["superseded_failed_a4"]["workflow_run"] == FAILED_A4_RUN
    assert payload["superseded_failed_a4"]["receipt_generated"] is False
    assert payload["superseded_failed_a4"]["scientific_reject"] is False

    frozen = payload["frozen_protocol"]
    assert frozen["offgrid_seed"] == 9173411
    assert frozen["offgrid_count"] == 8192
    assert frozen["identity_seed"] == 9173412
    assert frozen["identity_cases"] == 512
    assert frozen["max_public_ratio"] == 1.01

    ratios = payload["public_over_independent_ratios"]
    assert min(ratios.values()) >= 1.0
    assert max(ratios.values()) <= 1.01
    assert all(payload["negative_controls"].values())

    offgrid = payload["fresh_offgrid_eta_stress"]
    assert offgrid["fresh_random_offgrid_points"] == 8192
    assert offgrid["finite"] is True
    assert offgrid["dominated"] is True
    assert offgrid["axis_nontrivial"] is True
    assert offgrid["h_pm_0p1_percent_2D_response"] > 0.0

    identity = payload["manufactured_averaging_identity_stress"]
    assert identity["cases"] == 512
    assert identity["identity_passed"] is True
    assert identity["wrong_sign_mutation_detected"] is True
    assert identity["nontrivial"] is True


def test_repair_uses_public_product_api_not_missing_or_private_attribute() -> None:
    source = inspect.getsource(repair._observable_public_product_constant)
    assert ".product(" in source
    assert ".product_constant" not in source
    assert "._product" not in source
    assert "Fraction" not in source
    assert "bridge_post_j_replay" not in source


def test_repair_does_not_move_scientific_gates_or_truth_boundary() -> None:
    payload = run_repaired_independent_audit()
    repair_meta = payload["repair"]
    assert repair_meta["changed_scientific_formula"] is False
    assert repair_meta["changed_seed"] is False
    assert repair_meta["changed_mutation_factor"] is False
    assert repair_meta["changed_ratio_gate"] is False
    assert repair_meta["changed_final_pde_gate"] is False

    gates = payload["immutable_project_gates"]
    assert gates["normalized_momentum_max"] == 1e-3
    assert gates["normalized_momentum_L2"] == 1e-3
    assert gates["divergence_max"] == 1e-5
    assert gates["divergence_L2"] == 1e-5
    assert gates["free_residual_defined_forcing_allowed"] is False

    truth = payload["truth_boundary"]
    assert truth["failed_a4_head_remains_failed_without_receipt"] is True
    assert truth["repair_is_representation_only"] is True
    assert truth["all_R1_ordinary_slots_independently_audited"] is False
    assert truth["source_full_post_J2_R1_radius_one_ball_bound_machine_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["complete_kokuno_composite_velocity"] is False
    assert truth["leading_only_ns_residual_assessed"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
