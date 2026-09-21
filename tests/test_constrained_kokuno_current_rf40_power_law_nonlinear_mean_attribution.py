from __future__ import annotations

import inspect

import pytest

from openai_ns_reconstruction.kokuno_current_rf40_power_law_nonlinear_mean_attribution import (
    AGENT1_LEADING_HEAD,
    AGENT1_LEADING_PR,
    AGENT1_LEADING_SOURCE_BLOB,
    AGENT2_COMPOSITE_PARENT_HEAD,
    AGENT2_COMPOSITE_PARENT_PR,
    AGENT2_COMPOSITE_PARENT_SOURCE_BLOB,
    FIXED_LEADING_SPATIAL_STEP,
    ExactCurrentRF40PowerLawNonlinearBackend,
    _expected_source_provenance,
    _validate_power_law_configuration,
    materialize_current_rf40_power_law_nonlinear_mean_attribution,
    truth_boundary,
)


def _sha(ch: str) -> str:
    return ch * 64


def _configuration() -> dict[str, object]:
    return {
        "parent_agent2": {
            "pr": AGENT2_COMPOSITE_PARENT_PR,
            "head": AGENT2_COMPOSITE_PARENT_HEAD,
            "source_blob_sha1": AGENT2_COMPOSITE_PARENT_SOURCE_BLOB,
        },
        "agent1_leading": {
            "pr": AGENT1_LEADING_PR,
            "head": AGENT1_LEADING_HEAD,
            "source_blob_sha1": AGENT1_LEADING_SOURCE_BLOB,
            "semantic_sha256": _sha("a"),
            "configuration_sha256": _sha("b"),
        },
        "oscillatory_runtime": {
            "payload_sha256": _sha("c"),
        },
        "source_provenance": _expected_source_provenance(),
        "truth_boundary": {
            "current_leading_plus_oscillatory_velocity_through_RF40_power_law_materialized": True,
            "RF40_power_law_leading_plus_oscillatory_materialized": True,
            "full_concrete_oscillatory_runtime_digest_bound": True,
            "identity_preserving_RF40_power_law_composite_save_load_available": True,
            "reload_recomputes_leading_and_oscillatory_identities": True,
            "a2_semantic_identity_binds_corrected_source_provenance_block": True,
            "agent1_internal_external_source_semantic_binding_repaired_by_this_increment": False,
            "velocity_after_RF40_power_law_materialized": False,
            "full_post_XR_RF40_current_lineage_materialized": False,
            "current_lineage_cone_I1_I2_I3_I4_outer_overlays_completed": False,
            "outer_global_leading_velocity_materialized": False,
            "global_compact_support_completed": False,
            "agent3_correction_velocity_composed": False,
            "matched_global_pressure_materialized": False,
            "restricted_forcing_materialized": False,
            "complete_velocity_pressure_forcing_api": False,
            "heldout_ns_residual_assessed": False,
            "same_protocol_comparable_to_st006": False,
            "residual_reduction_claimed": False,
            "velocity_export_ready": False,
            "visual_correspondence_verified": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "pde_validated": False,
        },
    }


def test_power_law_configuration_accepts_only_pinned_truth_runtime_and_source():
    _validate_power_law_configuration(
        _configuration(),
        runtime_semantic_sha256=_sha("d"),
        runtime_oscillatory_sha256=_sha("c"),
    )


@pytest.mark.parametrize(
    ("path", "value", "message"),
    [
        (("parent_agent2", "head"), "bad", "parent head"),
        (("parent_agent2", "source_blob_sha1"), "bad", "parent source"),
        (("agent1_leading", "head"), "bad", "leading lineage"),
        (("oscillatory_runtime", "payload_sha256"), _sha("e"), "disagrees"),
        (("source_provenance", "commit"), "bad", "provenance"),
        (
            ("truth_boundary", "velocity_after_RF40_power_law_materialized"),
            True,
            "truth boundary",
        ),
        (
            ("truth_boundary", "matched_global_pressure_materialized"),
            True,
            "truth boundary",
        ),
    ],
)
def test_power_law_configuration_fails_closed_on_identity_source_or_truth_drift(
    path, value, message
):
    config = _configuration()
    cursor = config
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = value
    with pytest.raises(ValueError, match=message):
        _validate_power_law_configuration(
            config,
            runtime_semantic_sha256=_sha("d"),
            runtime_oscillatory_sha256=_sha("c"),
        )


def test_exact_binder_rejects_unpinned_backend_before_scientific_use():
    class FakeField:
        leading_backend = object()
        semantic_sha256 = _sha("d")
        oscillatory_runtime_sha256 = _sha("c")

        def velocity(self, x, y, z, t):
            return None

        def configuration(self):
            return _configuration()

    with pytest.raises(ValueError, match="module identity"):
        ExactCurrentRF40PowerLawNonlinearBackend.bind(
            FakeField(), lambda p, t: None
        )


def test_public_contract_has_no_scientific_tuning_knobs():
    params = inspect.signature(
        materialize_current_rf40_power_law_nonlinear_mean_attribution
    ).parameters
    assert list(params) == ["backend", "radius", "z", "t"]
    for forbidden in (
        "residual", "defect", "mean", "stress", "pressure", "forcing", "gain",
        "damping", "spatial_step", "derivative_step", "angular_order", "viscosity",
        "scientific_threshold", "correction", "stage_budget",
    ):
        assert forbidden not in params
    assert FIXED_LEADING_SPATIAL_STEP == 1.0e-3


def test_truth_boundary_advances_only_power_law_mean_scope():
    boundary = truth_boundary()
    assert boundary["current_leading_plus_oscillation_through_RF40_power_law_consumed"] is True
    assert boundary["RF40_power_law_velocity_materialized"] is True
    assert boundary["current_RF40_power_law_mixed_nonlinear_mean_materialized"] is True
    assert boundary["current_RF40_power_law_quadratic_nonlinear_mean_materialized"] is True
    assert boundary["current_RF40_power_law_aggregate_nonlinear_mean_materialized"] is True
    assert boundary["agent2_oscillatory_jacobian_reimplemented_by_agent3"] is False
    assert boundary["caller_tunable_leading_spatial_step"] is False
    assert boundary["caller_supplied_surrogate_defect_allowed"] is False
    assert boundary["forbidden_public_parameters_absent"] is True
    for key in (
        "velocity_after_RF40_power_law_materialized",
        "full_post_XR_RF40_current_lineage_materialized",
        "outer_global_leading_velocity_materialized",
        "radial_inverse_performed_in_this_increment",
        "pressure_gradient_included",
        "restricted_forcing_included",
        "complete_ns_defect",
        "scoped_nonlinear_mean_authorized_as_correction_target",
        "mean_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "paper_exact",
        "pde_validated",
        "blowup_proved",
    ):
        assert boundary[key] is False, key
    assert boundary["final_normalized_momentum_gate"] == 1.0e-3
    assert boundary["final_normalized_divergence_gate"] == 1.0e-5
