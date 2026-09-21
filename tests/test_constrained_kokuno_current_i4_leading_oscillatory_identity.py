from __future__ import annotations

import copy
import inspect
import math

import numpy as np
import pytest

from openai_ns_reconstruction import kokuno_current_i4_leading_oscillatory_identity as mod


def _synthetic_receipt() -> dict:
    h = "1" * 64
    x_start = 100.0
    x_end = 1000.0
    fractions = [float(v) for v in mod._I4_STAGE_FRACTIONS]
    xs = [
        math.exp(math.log(x_start) + f * (math.log(x_end) - math.log(x_start)))
        for f in fractions
    ]
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
            "i4_probe_count": 3,
            "i4_stage_fraction": fractions,
            "i4_similarity_X": xs,
            "X_I4_start": x_start,
            "X_I4_end": x_end,
            "inner_composition_closure_max_abs": 0.0,
            "i4_composition_closure_max_abs": 0.0,
            "inner_oscillatory_vector_rms": 1.0e-6,
            "i4_oscillatory_vector_rms_observed": 0.0,
            "axis_oscillatory_abs_max": 0.0,
            "outside_oscillatory_abs_max": 0.0,
            "velocity_save_load_replay_max_abs": 0.0,
        },
        "truth_boundary": mod._truth_boundary(),
        "frozen_mechanical_gates": {
            "inner_composition_closure_max_abs": mod.COMPOSITION_ATOL,
            "i4_composition_closure_max_abs": mod.COMPOSITION_ATOL,
            "inner_oscillatory_vector_rms_min": mod.OSCILLATORY_SIGNAL_FLOOR,
            "axis_oscillatory_abs_max": 0.0,
            "outside_oscillatory_abs_max": 0.0,
            "velocity_save_load_replay_max_abs": mod.VELOCITY_REPLAY_ATOL,
        },
    }


def test_parent_a2_is_exact_1071_source_blob() -> None:
    assert mod.PARENT_AGENT2_PR == 1071
    assert mod.PARENT_AGENT2_HEAD == "48d69e37e78e7f7f0e4e9936f28ff7974719288d"
    assert mod._validate_parent_a2_source() == mod.PARENT_AGENT2_SOURCE_BLOB_SHA1


def test_target_is_exact_current_i4_leading_sibling() -> None:
    assert mod.AGENT1_PR == 1079
    assert mod.AGENT1_HEAD == "b06742ca6e189499192ede3cce40f62cdc1e35ca"
    assert mod.AGENT1_CLASS == "KokunoPA16CurrentCartesianI4LeadingPreservation"
    assert mod.AGENT1_SOURCE_BLOB_SHA1 == "6f04ce0a856b44430402576dad88438da90d1ebb"
    truth = mod._required_agent1_truth()
    assert truth["current_leading_preserved_through_I4"] is True
    assert truth["source_I4_is_reserved_mean_correction_interval"] is True
    assert truth["source_leading_E0_U0_unchanged_on_I4"] is True
    assert truth["source_positive_order_I3_profiles_materialized"] is False
    assert truth["I3_positive_order_correction_materialized"] is False
    assert truth["source_I4_mean_correction_materialized"] is False
    assert truth["current_I4_mean_correction_materialized"] is False
    assert truth["leading_I4_overlay_invented"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_corrected_source_is_structural_provenance_not_exact_claim() -> None:
    source = mod._source_provenance()
    assert source["repository"] == "KokunoYumeto/yang-mills-interacting-workbench"
    assert source["commit"] == "143f6773feb424ad9ed3a8d116653200f20346b7"
    assert source["path"] == "navier-stokes/navier_stokes_workbench.tex"
    assert source["blob_sha1"] == "205a99807302e21a51c5eaf223390c0dfc42bcd0"
    assert source["corrected_release_date"] == "2026-09-09"
    assert source["classification"] == "structural_math_provenance_only"
    assert source["paper_exact_claim"] is False


def test_public_velocity_contract_has_no_tuning_or_correction_escape_hatch() -> None:
    assert list(inspect.signature(mod.velocity).parameters) == ["x", "y", "z", "t"]
    contract = mod.public_contract()
    assert contract["forbidden_velocity_inputs_present"] == []
    assert contract["current_leading_through_I4_materialized"] is True
    assert contract["current_I4_leading_plus_oscillatory_materialized"] is True
    assert contract["full_oscillatory_runtime_payload_bound"] is True
    assert contract["identity_preserving_save_load"] is True
    assert contract["source_I4_mean_correction_materialized"] is False
    assert contract["field_summands_changed"] is False
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


def test_truth_boundary_advances_only_current_i4_composition() -> None:
    truth = mod._truth_boundary()
    assert truth["current_leading_through_I4_consumed"] is True
    assert truth["current_I4_leading_plus_frozen_complete_curl_oscillation_materialized"] is True
    assert truth["full_concrete_oscillatory_runtime_digest_bound"] is True
    assert truth["identity_preserving_current_I4_composite_save_load_available"] is True
    assert truth["field_summands_changed_by_this_increment"] is False
    for key in (
        "source_positive_order_I3_profiles_materialized",
        "I3_positive_order_correction_materialized",
        "source_I4_mean_correction_materialized",
        "current_I4_mean_correction_materialized",
        "leading_I4_overlay_invented",
        "source_terminal_tail_schedule_bound_into_current_velocity",
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
        "source_forcing_provider_materialized",
        "source_background_provider_materialized",
        "paper_exact",
        "openai_field_identified",
        "pde_validated",
    ):
        assert truth[key] is False


def test_i4_probe_fractions_are_fixed_strictly_interior_and_replay_similarity_x() -> None:
    np.testing.assert_array_equal(mod._I4_STAGE_FRACTIONS, np.asarray((0.18, 0.57, 0.92)))

    class Dummy:
        D = 0.31
        X_I4_start = 50.0
        X_I4_end = 500.0

    fractions = mod._I4_STAGE_FRACTIONS
    eta = np.asarray((-0.3, 0.1, 0.4))
    t = np.asarray((0.37, 0.51, 0.69))
    theta = np.asarray((0.2, 0.8, 1.4))
    x, y, z, tt, X = mod._cartesian_i4_probe(Dummy(), fractions, eta, t, theta)
    q = (1.0 - tt) / (1.0 - eta * eta)
    replay_X = (x * x + y * y) / (2.0 * q)
    replay_eta = z / np.power(q, Dummy.D)
    assert np.allclose(replay_X, X, rtol=0.0, atol=2.0e-12)
    assert np.allclose(replay_eta, eta, rtol=0.0, atol=2.0e-15)
    assert np.all((X > Dummy.X_I4_start) & (X < Dummy.X_I4_end))


def test_synthetic_receipt_gate_accepts_only_frozen_scope() -> None:
    receipt = _synthetic_receipt()
    mod.enforce_receipt(receipt)

    promoted = copy.deepcopy(receipt)
    promoted["truth_boundary"]["pde_validated"] = True
    with pytest.raises(AssertionError, match="truth boundary"):
        mod.enforce_receipt(promoted)

    leaked = copy.deepcopy(receipt)
    leaked["diagnostic"]["outside_oscillatory_abs_max"] = 1.0e-12
    with pytest.raises(AssertionError, match="outside support"):
        mod.enforce_receipt(leaked)

    moved_parent = copy.deepcopy(receipt)
    moved_parent["parent_agent2_source_blob_sha1"] = "0" * 40
    with pytest.raises(AssertionError, match="parent source identity"):
        mod.enforce_receipt(moved_parent)

    fake_mean = copy.deepcopy(receipt)
    fake_mean["truth_boundary"]["source_I4_mean_correction_materialized"] = True
    with pytest.raises(AssertionError, match="truth boundary"):
        mod.enforce_receipt(fake_mean)


def test_sha256_guard_rejects_malformed_identity() -> None:
    with pytest.raises(ValueError):
        mod._require_sha256("abc", "bad")
    assert mod._require_sha256("a" * 64, "good") == "a" * 64
