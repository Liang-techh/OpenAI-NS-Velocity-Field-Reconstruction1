from __future__ import annotations

import copy
import inspect
import math

import numpy as np
import pytest

from openai_ns_reconstruction import (
    kokuno_current_rf40_first_turn_leading_oscillatory_identity as mod,
)


def _synthetic_receipt() -> dict:
    h = "1" * 64
    return {
        "schema": mod.RECEIPT_SCHEMA,
        "semantic_sha256": h,
        "saved_semantic_sha256": h,
        "oscillatory_runtime_sha256": "2" * 64,
        "agent1_semantic_sha256": "3" * 64,
        "agent1_configuration_sha256": "4" * 64,
        "parent_agent2_source_blob_sha1": mod.PARENT_AGENT2_SOURCE_BLOB_SHA1,
        "public_contract": mod.public_contract(),
        "diagnostic": {
            "inner_probe_count": 6,
            "post_XR_probe_count": 3,
            "post_XR_X_over_XR": [math.exp(0.18), math.exp(0.57), math.e],
            "post_XR_X_over_X1": [math.exp(-0.82), math.exp(-0.43), 1.0],
            "inner_composition_closure_max_abs": 0.0,
            "post_XR_composition_closure_max_abs": 0.0,
            "inner_oscillatory_vector_rms": 1.0e-6,
            "post_XR_oscillatory_abs_max": 0.0,
            "velocity_save_load_replay_max_abs": 0.0,
            "beyond_current_X1_fail_closed": True,
            "X_h": 1.0,
            "X_R": math.exp(5.0),
            "X_1": math.exp(6.0),
        },
        "truth_boundary": mod._truth_boundary(),
        "frozen_mechanical_gates": {
            "inner_composition_closure_max_abs": mod.COMPOSITION_ATOL,
            "post_XR_composition_closure_max_abs": mod.COMPOSITION_ATOL,
            "inner_oscillatory_vector_rms_min": mod.OSCILLATORY_SIGNAL_FLOOR,
            "post_XR_oscillatory_abs_max": 0.0,
            "velocity_save_load_replay_max_abs": mod.VELOCITY_REPLAY_ATOL,
            "beyond_current_X1_fail_closed": True,
        },
    }


def test_parent_agent2_source_is_exactly_the_stacked_pr987_blob() -> None:
    assert mod.PARENT_AGENT2_PR == 987
    assert mod.PARENT_AGENT2_HEAD == "4be2c9ee898c44dd1ad2217e90601b161fe81964"
    assert mod._validate_parent_a2_source() == mod.PARENT_AGENT2_SOURCE_BLOB_SHA1


def test_agent1_target_is_exact_rf40_first_turn_sibling() -> None:
    assert mod.AGENT1_PR == 986
    assert mod.AGENT1_HEAD == "23c00b98526e187ff04d432745300a084ec859f2"
    assert mod.AGENT1_CLASS == "KokunoPA16CurrentCartesianRF40FirstTurn"
    truth = mod._required_agent1_truth()
    assert truth["current_lineage_RF40_first_turn_materialized"] is True
    assert truth["full_post_XR_RF40_current_lineage_materialized"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_public_velocity_contract_has_no_tuning_escape_hatch() -> None:
    assert list(inspect.signature(mod.velocity).parameters) == ["x", "y", "z", "t"]
    contract = mod.public_contract()
    assert contract["forbidden_velocity_inputs_present"] == []
    assert contract["leading_materialized_through_RF40_first_turn"] is True
    assert contract["full_oscillatory_runtime_payload_bound"] is True
    assert contract["identity_preserving_save_load"] is True
    assert contract["after_RF40_first_turn_fail_closed"] is True
    assert contract["new_oscillatory_parameters"] is False
    assert contract["leading_profile_reimplemented"] is False
    assert contract["mean_projection_performed"] is False
    assert contract["radial_inverse_performed"] is False
    assert contract["correction_velocity_constructed"] is False
    assert contract["pressure_or_forcing_added"] is False
    assert contract["complete_ns_residual"] is False
    assert contract["velocity_export_ready"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False


def test_truth_boundary_advances_only_first_turn_composition() -> None:
    truth = mod._truth_boundary()
    assert truth["current_leading_plus_oscillatory_velocity_through_RF40_first_turn_materialized"] is True
    assert truth["post_XR_RF40_first_turn_leading_plus_oscillatory_materialized"] is True
    assert truth["full_concrete_oscillatory_runtime_digest_bound"] is True
    assert truth["identity_preserving_RF40_first_turn_composite_save_load_available"] is True
    for key in (
        "velocity_after_RF40_first_turn_materialized",
        "full_post_XR_RF40_current_lineage_materialized",
        "RF40_axial_shutdown_current_lineage_materialized",
        "RF40_lambda_turn_current_lineage_materialized",
        "RF40_power_law_current_lineage_materialized",
        "current_lineage_cone_I1_I2_I3_I4_outer_overlays_completed",
        "outer_global_leading_velocity_materialized",
        "global_compact_support_completed",
        "agent3_correction_velocity_composed",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_velocity_pressure_forcing_api",
        "heldout_ns_residual_assessed",
        "same_protocol_comparable_to_st006",
        "residual_reduction_claimed",
        "velocity_export_ready",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "pde_validated",
    ):
        assert truth[key] is False


def test_cartesian_probe_replays_requested_similarity_X() -> None:
    class Dummy:
        D = 0.31

    X = np.asarray((0.7, 1.4, 2.1))
    eta = np.asarray((-0.3, 0.1, 0.4))
    t = np.asarray((0.37, 0.51, 0.69))
    theta = np.asarray((0.2, 0.8, 1.4))
    x, y, z, tt = mod._cartesian_probe(Dummy(), X, eta, t, theta)
    tau = 1.0 - tt
    q = tau / (1.0 - eta * eta)
    replay_X = (x * x + y * y) / (2.0 * q)
    replay_eta = z / np.power(q, Dummy.D)
    assert np.allclose(replay_X, X, rtol=0.0, atol=2.0e-15)
    assert np.allclose(replay_eta, eta, rtol=0.0, atol=2.0e-15)


def test_synthetic_receipt_gate_accepts_only_frozen_scope() -> None:
    receipt = _synthetic_receipt()
    mod.enforce_receipt(receipt)

    leaked = copy.deepcopy(receipt)
    leaked["diagnostic"]["post_XR_oscillatory_abs_max"] = 1.0e-12
    with pytest.raises(AssertionError, match="leaked"):
        mod.enforce_receipt(leaked)

    promoted = copy.deepcopy(receipt)
    promoted["truth_boundary"]["pde_validated"] = True
    with pytest.raises(AssertionError, match="truth boundary"):
        mod.enforce_receipt(promoted)

    moved_parent = copy.deepcopy(receipt)
    moved_parent["parent_agent2_source_blob_sha1"] = "0" * 40
    with pytest.raises(AssertionError, match="parent source identity"):
        mod.enforce_receipt(moved_parent)


def test_sha256_guard_rejects_malformed_identity() -> None:
    with pytest.raises(ValueError):
        mod._require_sha256("abc", "bad")
    assert mod._require_sha256("a" * 64, "good") == "a" * 64
