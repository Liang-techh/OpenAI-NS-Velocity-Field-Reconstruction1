from __future__ import annotations

import copy
import inspect

import pytest

import openai_ns_reconstruction.kokuno_current_exterior_xr_leading_oscillatory_identity as mod


def _valid_receipt():
    return {
        "schema": mod.RECEIPT_SCHEMA,
        "semantic_sha256": "1" * 64,
        "saved_semantic_sha256": "1" * 64,
        "oscillatory_runtime_sha256": "2" * 64,
        "agent1_semantic_sha256": "3" * 64,
        "agent1_configuration_sha256": "4" * 64,
        "public_contract": mod.public_contract(),
        "diagnostic": {
            "inner_probe_count": 6,
            "exterior_probe_count": 3,
            "exterior_X_over_Xh": [2.0, 10.0, 148.0],
            "exterior_X_over_XR": [0.01, 0.1, 1.0],
            "inner_composition_closure_max_abs": 0.0,
            "exterior_composition_closure_max_abs": 0.0,
            "inner_oscillatory_vector_rms": 1.0e-3,
            "exterior_oscillatory_abs_max": 0.0,
            "velocity_save_load_replay_max_abs": 0.0,
            "beyond_current_XR_fail_closed": True,
            "X_h": 1.0,
            "X_R": 10.0,
        },
        "truth_boundary": mod._truth_boundary(),
        "frozen_mechanical_gates": {
            "inner_composition_closure_max_abs": mod.COMPOSITION_ATOL,
            "exterior_composition_closure_max_abs": mod.COMPOSITION_ATOL,
            "inner_oscillatory_vector_rms_min": mod.OSCILLATORY_SIGNAL_FLOOR,
            "exterior_oscillatory_abs_max": 0.0,
            "velocity_save_load_replay_max_abs": mod.VELOCITY_REPLAY_ATOL,
            "beyond_current_XR_fail_closed": True,
        },
    }


def test_public_velocity_surface_has_no_scientific_tuning_knobs() -> None:
    assert list(inspect.signature(mod.velocity).parameters) == ["x", "y", "z", "t"]
    contract = mod.public_contract()
    assert contract["forbidden_velocity_inputs_present"] == []
    assert contract["exact_agent1_pr"] == 980
    assert contract["leading_materialized_through_XR"] is True
    assert contract["full_oscillatory_runtime_payload_bound"] is True
    assert contract["identity_preserving_save_load"] is True
    assert contract["post_XR_fail_closed"] is True
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


def test_truth_boundary_advances_only_current_composition_through_XR() -> None:
    truth = mod._truth_boundary()
    assert truth["current_leading_plus_oscillatory_velocity_through_XR_materialized"] is True
    assert truth["velocity_beyond_Xh_through_XR_materialized"] is True
    assert truth["full_concrete_oscillatory_runtime_digest_bound"] is True
    assert truth["identity_preserving_through_XR_composite_save_load_available"] is True
    assert truth["reload_recomputes_leading_and_oscillatory_identities"] is True
    for key in (
        "post_XR_velocity_materialized",
        "post_XR_RF40_current_lineage_materialized",
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


def test_exact_parent_and_agent1_identities_are_pinned() -> None:
    assert mod.PARENT_AGENT2_PR == 981
    assert mod.PARENT_AGENT2_HEAD == "9997fc55455d126f935643da36bf17eaa0491aa4"
    assert mod.AGENT1_PR == 980
    assert mod.AGENT1_HEAD == "d3c971f2c62e272333e124e532212d23cca4908d"
    assert mod.AGENT1_MODULE == (
        "openai_ns_reconstruction.kokuno_pa16_current_cartesian_exterior_to_xr"
    )
    assert mod.AGENT1_CLASS == "KokunoPA16CurrentCartesianExteriorToXR"
    assert mod.AGENT1_SOURCE_BLOB_SHA1 == "2ba28636a56b252fa485719e7b3e8da754d7e588"
    assert mod.AGENT1_TEST_BLOB_SHA1 == "ecbab9822c54f440f5555f25cd5e7e1648d84bd2"
    assert mod.AGENT1_WORKFLOW_BLOB_SHA1 == "5c2d68211e8fea526e06087fc4475f7bab9eb563"
    assert mod.PUBLIC_Z_OSCILLATORY_SOURCE_BLOB_SHA1 == (
        "4a3abc1a11f7e054f651dd7bace8ada6d5a6b653"
    )


def test_enforcer_accepts_frozen_through_XR_mechanics() -> None:
    mod.enforce_receipt(_valid_receipt())


@pytest.mark.parametrize(
    "mutate,match",
    [
        (
            lambda r: r["diagnostic"].__setitem__(
                "inner_composition_closure_max_abs", 10.0 * mod.COMPOSITION_ATOL
            ),
            "inner leading",
        ),
        (
            lambda r: r["diagnostic"].__setitem__(
                "exterior_composition_closure_max_abs", 10.0 * mod.COMPOSITION_ATOL
            ),
            "exterior leading",
        ),
        (
            lambda r: r["diagnostic"].__setitem__("inner_oscillatory_vector_rms", 0.0),
            "vacuous",
        ),
        (
            lambda r: r["diagnostic"].__setitem__("exterior_oscillatory_abs_max", 1.0e-9),
            "leaked",
        ),
        (
            lambda r: r["diagnostic"].__setitem__(
                "velocity_save_load_replay_max_abs", 10.0 * mod.VELOCITY_REPLAY_ATOL
            ),
            "save/load",
        ),
        (
            lambda r: r["diagnostic"].__setitem__("beyond_current_XR_fail_closed", False),
            "fails closed",
        ),
        (lambda r: r.__setitem__("saved_semantic_sha256", "5" * 64), "semantic identity"),
        (
            lambda r: r["truth_boundary"].__setitem__("velocity_export_ready", True),
            "illegally promoted",
        ),
    ],
)
def test_enforcer_rejects_mechanical_or_truth_mutations(mutate, match) -> None:
    receipt = copy.deepcopy(_valid_receipt())
    mutate(receipt)
    with pytest.raises((AssertionError, ValueError), match=match):
        mod.enforce_receipt(receipt)
