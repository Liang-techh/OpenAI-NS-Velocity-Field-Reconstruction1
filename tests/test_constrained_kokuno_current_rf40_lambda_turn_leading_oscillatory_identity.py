from __future__ import annotations

import copy
import inspect
import math

import numpy as np
import pytest

from openai_ns_reconstruction import (
    kokuno_current_rf40_lambda_turn_leading_oscillatory_identity as mod,
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
        "source_provenance_sha256": mod._identity._sha256(mod._source_provenance()),
        "public_contract": mod.public_contract(),
        "diagnostic": {
            "inner_probe_count": 6,
            "lambda_turn_probe_count": 3,
            "lambda_turn_stage_fraction": [0.18, 0.57, 1.0],
            "lambda_turn_X_over_X2": [math.exp(0.18), math.exp(0.57), math.e],
            "lambda_turn_X_over_X3": [math.exp(-0.82), math.exp(-0.43), 1.0],
            "inner_composition_closure_max_abs": 0.0,
            "lambda_turn_composition_closure_max_abs": 0.0,
            "inner_oscillatory_vector_rms": 1.0e-6,
            "lambda_turn_oscillatory_abs_max": 0.0,
            "velocity_save_load_replay_max_abs": 0.0,
            "beyond_current_X3_fail_closed": True,
            "X_R": 1.0,
            "X_2": math.e * math.e,
            "X_3": math.e * math.e * math.e,
        },
        "truth_boundary": mod._truth_boundary(),
        "frozen_mechanical_gates": {
            "inner_composition_closure_max_abs": mod.COMPOSITION_ATOL,
            "lambda_turn_composition_closure_max_abs": mod.COMPOSITION_ATOL,
            "inner_oscillatory_vector_rms_min": mod.OSCILLATORY_SIGNAL_FLOOR,
            "lambda_turn_oscillatory_abs_max": 0.0,
            "velocity_save_load_replay_max_abs": mod.VELOCITY_REPLAY_ATOL,
            "beyond_current_X3_fail_closed": True,
        },
    }


def test_parent_agent2_source_is_exactly_the_stacked_pr999_blob() -> None:
    assert mod.PARENT_AGENT2_PR == 999
    assert mod.PARENT_AGENT2_HEAD == "93b99292fcf141f24b7c6d7e4fbf95775e6a07e5"
    assert mod._validate_parent_a2_source() == mod.PARENT_AGENT2_SOURCE_BLOB_SHA1


def test_agent1_target_is_exact_rf40_lambda_turn_sibling() -> None:
    assert mod.AGENT1_PR == 998
    assert mod.AGENT1_HEAD == "43b295444b1e9558222d757cb551e04385eddcb5"
    assert mod.AGENT1_CLASS == "KokunoPA16CurrentCartesianRF40LambdaTurn"
    truth = mod._required_agent1_truth()
    assert truth["current_lineage_RF40_lambda_turn_materialized"] is True
    assert truth["RF40_power_law_current_lineage_materialized"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_source_provenance_is_explicit_and_non_paper_exact() -> None:
    source = mod._source_provenance()
    assert source["repository"] == "KokunoYumeto/yang-mills-interacting-workbench"
    assert source["commit"] == "143f6773feb424ad9ed3a8d116653200f20346b7"
    assert source["path"] == "navier-stokes/navier_stokes_workbench.tex"
    assert source["blob_sha1"] == "205a99807302e21a51c5eaf223390c0dfc42bcd0"
    assert source["corrected_release_date"] == "2026-09-09"
    assert source["license"] is None
    assert source["classification"] == "reimplement_math_only"
    assert source["migration_scope"] == "provenance_only"
    assert source["paper_exact_claim"] is False


def test_public_velocity_contract_has_no_tuning_escape_hatch() -> None:
    assert list(inspect.signature(mod.velocity).parameters) == ["x", "y", "z", "t"]
    contract = mod.public_contract()
    assert contract["forbidden_velocity_inputs_present"] == []
    assert contract["leading_materialized_through_RF40_lambda_turn"] is True
    assert contract["full_oscillatory_runtime_payload_bound"] is True
    assert contract["identity_preserving_save_load"] is True
    assert contract["after_RF40_lambda_turn_fail_closed"] is True
    assert contract["a2_source_provenance_bound"] is True
    assert contract["repairs_agent1_internal_source_semantic_binding"] is False
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


def test_truth_boundary_advances_only_lambda_turn_composition() -> None:
    truth = mod._truth_boundary()
    assert truth["current_leading_plus_oscillatory_velocity_through_RF40_lambda_turn_materialized"] is True
    assert truth["RF40_lambda_turn_leading_plus_oscillatory_materialized"] is True
    assert truth["full_concrete_oscillatory_runtime_digest_bound"] is True
    assert truth["identity_preserving_RF40_lambda_turn_composite_save_load_available"] is True
    assert truth["a2_semantic_identity_binds_corrected_source_provenance_block"] is True
    assert truth["agent1_internal_external_source_semantic_binding_repaired_by_this_increment"] is False
    for key in (
        "velocity_after_RF40_lambda_turn_materialized",
        "full_post_XR_RF40_current_lineage_materialized",
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
    leaked["diagnostic"]["lambda_turn_oscillatory_abs_max"] = 1.0e-12
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

    moved_source = copy.deepcopy(receipt)
    moved_source["source_provenance_sha256"] = "0" * 64
    with pytest.raises(AssertionError, match="source provenance"):
        mod.enforce_receipt(moved_source)


def test_sha256_guard_rejects_malformed_identity() -> None:
    with pytest.raises(ValueError):
        mod._require_sha256("abc", "bad")
    assert mod._require_sha256("a" * 64, "good") == "a" * 64
