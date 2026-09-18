from __future__ import annotations

import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.kokuno_bounded_coordinate_independent_routing_checkpoint import (
    AGENT4_COORDINATE_AUDIT_RECEIPT,
    SCHEMA,
    build_checkpoint,
    write_checkpoint,
)
from openai_ns_reconstruction.kokuno_bounded_coordinate_routing_checkpoint import checkpoint_sha256


def test_agent4_generic_bounded_coordinate_audit_is_bound_as_pass() -> None:
    receipt = AGENT4_COORDINATE_AUDIT_RECEIPT
    assert receipt["source_pr"] == 402
    assert receipt["repository_bounded_coordinate_independent_preflight_passed"] is True
    assert receipt["actual_source_multiband_family_audited"] is False
    assert receipt["common_finest_covariance_response_relative_rms"] <= 1.0e-2
    assert receipt["band_finest_covariance_response_relative_rms"] <= 1.0e-2
    assert min(receipt["common_refinement_ratios"]) >= 1.8
    assert min(receipt["band_refinement_ratios"]) >= 1.8
    assert receipt["common_quadratic_identity_relative_error"] <= 5.0e-13
    assert receipt["simultaneous_label_velocity_permutation_relative_rms"] <= 5.0e-13
    assert receipt["minimum_label_multiplier_at_l1_0125"] >= 0.875 - 5.0e-13
    assert receipt["misaligned_outer_band_response_relative_change"] >= 0.1


def test_independent_coordinate_pass_does_not_promote_missing_physical_family() -> None:
    checkpoint = build_checkpoint()
    assert checkpoint["schema"] == SCHEMA
    assert checkpoint["upstream"]["agent4"]["latest_pr"] == 402
    assert checkpoint["states"]["repository_bounded_coordinate_independently_audited"] is True
    assert checkpoint["states"]["actual_multiband_physical_family_independently_audited"] is False
    assert checkpoint["states"]["genuinely_independent_second_covariance_column_ready"] is False
    assert checkpoint["states"]["correction_ready"] is False
    assert checkpoint["states"]["candidate_artifact_instantiated"] is False
    assert checkpoint["states"]["velocity_export_ready"] is False
    assert checkpoint["states"]["formal_full_domain_pde_gate_assessed"] is False
    assert checkpoint["states"]["pde_validated"] is False


def test_shortest_route_moves_past_generic_coordinate_audit() -> None:
    checkpoint = build_checkpoint()
    routing = checkpoint["routing"]
    joined = " ".join(routing["shortest_next_closure"]).lower()
    forbidden = " ".join(routing["do_not_do"]).lower()
    assert "at least two active slow bands" in joined
    assert "actual-family tangent/rank receipt" in joined
    assert "do not rerun the generic bounded-coordinate audit" in forbidden
    assert "0/25" in forbidden


def test_checkpoint_v23_is_deterministic_hash_bound_and_serializable(tmp_path: Path) -> None:
    first = build_checkpoint()
    second = build_checkpoint()
    assert first == second
    stored_hash = first["checkpoint_sha256"]
    unhashed = dict(first)
    del unhashed["checkpoint_sha256"]
    assert stored_hash == checkpoint_sha256(unhashed)

    output = tmp_path / "checkpoint.json"
    written = write_checkpoint(output)
    loaded = json.loads(output.read_text(encoding="utf-8"))
    assert loaded == written == first


def test_formal_gates_remain_unchanged() -> None:
    gates = build_checkpoint()["fixed_gates"]
    assert gates["held_out_normalized_full_momentum_max"] == pytest.approx(1.0e-3)
    assert gates["held_out_normalized_full_momentum_l2"] == pytest.approx(1.0e-3)
    assert gates["divergence_max"] == pytest.approx(1.0e-5)
    assert gates["divergence_l2"] == pytest.approx(1.0e-5)
    assert gates["changed"] is False
