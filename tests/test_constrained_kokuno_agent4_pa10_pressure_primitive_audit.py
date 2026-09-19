from __future__ import annotations

import math

from openai_ns_reconstruction.kokuno_agent4_pa10_pressure_primitive_audit import (
    AGENT1_CANDIDATE_SHA256,
    AGENT1_HEAD,
    FINAL_PROJECT_GATES,
    LOCAL_GATES,
    SCHEMA,
    SEED,
    run_audit,
)


def test_independent_pressure_audit_is_deterministic_and_fail_closed() -> None:
    receipt = run_audit()

    assert receipt["schema"] == SCHEMA
    assert receipt["provenance"]["agent1_exact_head"] == AGENT1_HEAD
    assert receipt["provenance"]["agent1_candidate_sha256"] == AGENT1_CANDIDATE_SHA256
    assert receipt["frozen_protocol"]["seed"] == SEED
    assert receipt["frozen_protocol"]["local_gates"] == LOCAL_GATES
    assert receipt["frozen_protocol"]["final_project_gates_unchanged"] == FINAL_PROJECT_GATES

    independence = receipt["independence_contract"]
    assert independence["agent1_private_integral_helper_called"] is False
    assert independence["agent1_axis_state_called_by_validator"] is False
    assert independence["agent1_f0_prime_called_by_validator"] is False
    assert independence["agent1_quadrature_used_as_validation_operator"] is False

    truth = receipt["truth_boundary"]
    assert truth["source_pressure_radius_one_ball_bound_assessed"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["full_same_cycle_requested_stress_materialized"] is False
    assert truth["public_velocity_correction_materialized"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False

    for channel in ("p_Y", "Y_Phi_Y", "Phi_eta", "p_eta"):
        levels = receipt["results"][channel]["levels"]
        assert len(levels) == 3
        for level in levels:
            assert math.isfinite(level["relative_rms"])
            assert math.isfinite(level["relative_max"])

    fto_levels = receipt["results"]["fundamental_theorem_Y"]["levels"]
    assert [level["node_count"] for level in fto_levels] == [65, 129, 257]
    assert all(math.isfinite(level["relative_rms"]) for level in fto_levels)


def test_final_gates_are_not_interface_gates() -> None:
    assert FINAL_PROJECT_GATES["normalized_momentum_max"] == 1.0e-3
    assert FINAL_PROJECT_GATES["normalized_momentum_L2"] == 1.0e-3
    assert FINAL_PROJECT_GATES["divergence_max"] == 1.0e-5
    assert FINAL_PROJECT_GATES["divergence_L2"] == 1.0e-5
    assert LOCAL_GATES["p_eta_relative_rms"] != FINAL_PROJECT_GATES["normalized_momentum_L2"]
