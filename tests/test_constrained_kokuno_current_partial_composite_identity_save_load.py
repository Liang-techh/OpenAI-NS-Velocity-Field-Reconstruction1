from __future__ import annotations

import copy
import inspect

import pytest

import openai_ns_reconstruction.kokuno_current_partial_composite_identity_save_load as mod


class _FakeOscillatoryRuntime:
    def __init__(self, amplitude: float) -> None:
        self.amplitude = float(amplitude)

    def to_payload(self):
        return {
            "schema": "fake-osc-runtime",
            "parameters": {"amplitude": self.amplitude},
            "frozen_realization": {
                "modes": [[1, 2], [3, 5]],
                "support": {"radial": [0.15, 1.35], "axial": [-2.0, 2.0]},
            },
            "truth_boundary": {"paper_exact": False, "pde_validated": False},
        }


def _valid_receipt():
    semantic = "1" * 64
    oscillatory = "2" * 64
    parent = "3" * 64
    return {
        "schema": mod.RECEIPT_SCHEMA,
        "semantic_sha256": semantic,
        "saved_semantic_sha256": semantic,
        "reloaded_semantic_sha256": semantic,
        "oscillatory_runtime_sha256": oscillatory,
        "saved_oscillatory_runtime_sha256": oscillatory,
        "reloaded_oscillatory_runtime_sha256": oscillatory,
        "parent_composite_semantic_sha256": parent,
        "reloaded_parent_composite_semantic_sha256": parent,
        "diagnostic": {
            "velocity_save_load_replay_max_abs": 0.0,
            "oscillatory_increment_vector_rms": 1.0e-3,
            "beyond_current_Xh_fail_closed_after_reload": True,
        },
        "truth_boundary": mod._truth_boundary(),
        "frozen_mechanical_gates": {
            "velocity_save_load_replay_max_abs": mod.VELOCITY_REPLAY_ATOL,
            "oscillatory_increment_vector_rms_min": mod.OSCILLATORY_SIGNAL_FLOOR,
            "semantic_identity_exact_replay": True,
            "oscillatory_runtime_identity_exact_replay": True,
            "parent_composite_identity_exact_replay": True,
            "beyond_current_Xh_fail_closed_after_reload": True,
        },
    }


def test_complete_runtime_digest_sees_concrete_realization_drift() -> None:
    a = mod.oscillatory_runtime_sha256(_FakeOscillatoryRuntime(1.0))
    b = mod.oscillatory_runtime_sha256(_FakeOscillatoryRuntime(2.0))
    assert len(a) == len(b) == 64
    assert a != b


def test_complete_runtime_digest_is_canonical_under_mapping_order() -> None:
    payload_a = {
        "truth_boundary": {"paper_exact": False, "pde_validated": False},
        "frozen_realization": {"b": 2, "a": 1},
    }
    payload_b = {
        "frozen_realization": {"a": 1, "b": 2},
        "truth_boundary": {"pde_validated": False, "paper_exact": False},
    }

    class Runtime:
        def __init__(self, payload):
            self.payload = payload

        def to_payload(self):
            return self.payload

    assert mod.oscillatory_runtime_sha256(Runtime(payload_a)) == mod.oscillatory_runtime_sha256(
        Runtime(payload_b)
    )


def test_runtime_payload_rejects_nonfinite_or_promoted_truth() -> None:
    class BadNaN:
        def to_payload(self):
            return {
                "value": float("nan"),
                "truth_boundary": {"paper_exact": False, "pde_validated": False},
            }

    class BadTruth:
        def to_payload(self):
            return {"truth_boundary": {"paper_exact": True, "pde_validated": False}}

    with pytest.raises(ValueError, match="finite and JSON serializable"):
        mod.oscillatory_runtime_sha256(BadNaN())
    with pytest.raises(RuntimeError, match="truth boundary"):
        mod.oscillatory_runtime_sha256(BadTruth())


def test_public_velocity_surface_has_no_scientific_tuning_knobs() -> None:
    assert list(inspect.signature(mod.velocity).parameters) == ["x", "y", "z", "t"]
    contract = mod.public_contract()
    assert contract["forbidden_velocity_inputs_present"] == []
    assert contract["full_oscillatory_runtime_payload_bound"] is True
    assert contract["leading_configuration_and_semantic_bound"] is True
    assert contract["identity_preserving_partial_composite_save_load"] is True
    assert contract["reload_recomputes_runtime_digests"] is True
    assert contract["field_formula_changed"] is False
    assert contract["new_oscillatory_parameters"] is False
    assert contract["mean_projection_performed"] is False
    assert contract["radial_inverse_performed"] is False
    assert contract["correction_velocity_constructed"] is False
    assert contract["pressure_or_forcing_added"] is False
    assert contract["complete_ns_residual"] is False
    assert contract["velocity_export_ready"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False


def test_truth_boundary_closes_only_identity_save_load_gap() -> None:
    truth = mod._truth_boundary()
    assert truth["current_partial_leading_plus_oscillatory_velocity_materialized"] is True
    assert truth["full_concrete_oscillatory_runtime_digest_bound"] is True
    assert truth["identity_preserving_partial_composite_save_load_available"] is True
    assert truth["reload_recomputes_leading_and_oscillatory_identities"] is True
    for key in (
        "velocity_beyond_Xh_materialized",
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


def test_enforcer_accepts_exact_identity_replay_mechanics() -> None:
    mod.enforce_receipt(_valid_receipt())


@pytest.mark.parametrize(
    "mutate,match",
    [
        (lambda r: r.__setitem__("reloaded_semantic_sha256", "4" * 64), "semantic identity"),
        (
            lambda r: r.__setitem__("reloaded_oscillatory_runtime_sha256", "4" * 64),
            "oscillatory runtime identity",
        ),
        (
            lambda r: r.__setitem__("reloaded_parent_composite_semantic_sha256", "4" * 64),
            "parent partial-composite identity",
        ),
        (
            lambda r: r["diagnostic"].__setitem__(
                "velocity_save_load_replay_max_abs", 10.0 * mod.VELOCITY_REPLAY_ATOL
            ),
            "velocity replay",
        ),
        (
            lambda r: r["diagnostic"].__setitem__("oscillatory_increment_vector_rms", 0.0),
            "vacuous",
        ),
        (
            lambda r: r["diagnostic"].__setitem__(
                "beyond_current_Xh_fail_closed_after_reload", False
            ),
            "fails closed",
        ),
        (
            lambda r: r["truth_boundary"].__setitem__("velocity_export_ready", True),
            "illegally promoted",
        ),
    ],
)
def test_enforcer_rejects_identity_or_truth_mutations(mutate, match) -> None:
    receipt = copy.deepcopy(_valid_receipt())
    mutate(receipt)
    with pytest.raises(AssertionError, match=match):
        mod.enforce_receipt(receipt)


def test_parent_and_source_identities_are_exactly_pinned() -> None:
    assert mod.PARENT_AGENT2_PR == 975
    assert mod.PARENT_AGENT2_HEAD == "b716dfe495763439ac51f6e866211337910017cb"
    assert mod.COMPOSITION_AGENT2_PR == 970
    assert mod.COMPOSITION_AGENT2_HEAD == "3a6405bbd3d10b8c3b38078f0c989c45e5d407b4"
    assert mod.AGENT1_PR == 965
    assert mod.AGENT1_HEAD == "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1"
    assert mod.AGENT1_SOURCE_BLOB_SHA1 == "ca8b80b0451be1a8f31deaf620f542d8c2e92c0c"
    assert mod.PUBLIC_Z_OSCILLATORY_SOURCE_BLOB_SHA1 == "4a3abc1a11f7e054f651dd7bace8ada6d5a6b653"
