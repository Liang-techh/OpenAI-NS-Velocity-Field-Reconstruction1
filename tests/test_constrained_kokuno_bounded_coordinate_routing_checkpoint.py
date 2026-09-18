from __future__ import annotations

import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.kokuno_bounded_coordinate_routing_checkpoint import (
    AGENT2_COORDINATE_RECEIPT,
    AGENT3_BRIDGE_RECEIPT,
    FIXED_GATES,
    SCHEMA,
    STATES,
    build_checkpoint,
    checkpoint_sha256,
    write_checkpoint,
)


def test_coordinate_budget_conversion_is_unit_safe_and_unchanged() -> None:
    frozen = AGENT3_BRIDGE_RECEIPT["frozen_physical_amplitude_budget"]
    amplitude = AGENT3_BRIDGE_RECEIPT["reference_physical_amplitude"]
    converted = AGENT3_BRIDGE_RECEIPT["converted_fractional_budget"]
    assert converted == pytest.approx(frozen / abs(amplitude), rel=1e-15, abs=0.0)
    assert converted == pytest.approx(0.10010444711945854, rel=0.0, abs=1e-16)
    assert AGENT3_BRIDGE_RECEIPT["budget_changed"] is False
    assert AGENT3_BRIDGE_RECEIPT["budget_unit_matches_jacobian"] is True
    assert AGENT2_COORDINATE_RECEIPT["parameter_count"] == 2
    assert AGENT2_COORDINATE_RECEIPT["repository_coefficient_unit_mapping_available"] is True
    assert AGENT2_COORDINATE_RECEIPT["source_coefficient_unit_mapping_available"] is False


def test_one_band_control_stays_rank_deficient_and_finite_cycle_closed() -> None:
    receipt = AGENT3_BRIDGE_RECEIPT
    assert receipt["one_band_band_coordinate_inactive"] is True
    assert receipt["band_coordinate_response_vector_rms"] == 0.0
    assert receipt["nodes_requiring_second_direction"] == 25
    assert receipt["rank2_required_nodes"] == 0
    assert receipt["rank_deficient_required_nodes"] == 25
    assert receipt["algebraic_relative_stress_residual_rms"] == pytest.approx(
        0.2101945021227428, rel=0.0, abs=1e-15
    )
    assert receipt["family_bounded_inverse_preflight_passed"] is False
    assert receipt["finite_correction_cycle_authorized"] is False
    assert STATES["genuinely_independent_second_covariance_column_ready"] is False
    assert STATES["correction_ready"] is False


def test_global_truth_states_and_fixed_gates_remain_fail_closed() -> None:
    assert FIXED_GATES == {
        "held_out_normalized_full_momentum_max": 1.0e-3,
        "held_out_normalized_full_momentum_l2": 1.0e-3,
        "divergence_max": 1.0e-5,
        "divergence_l2": 1.0e-5,
        "changed": False,
    }
    assert STATES["leading_ready"] is False
    assert STATES["oscillatory_ready"] is False
    assert STATES["candidate_artifact_instantiated"] is False
    assert STATES["velocity_export_ready"] is False
    assert STATES["formal_full_domain_pde_gate_assessed"] is False
    assert STATES["pde_validated"] is False


def test_checkpoint_is_deterministic_and_hash_bound(tmp_path: Path) -> None:
    first = build_checkpoint()
    second = build_checkpoint()
    assert first == second
    assert first["schema"] == SCHEMA
    stored_hash = first["checkpoint_sha256"]
    unhashed = dict(first)
    del unhashed["checkpoint_sha256"]
    assert stored_hash == checkpoint_sha256(unhashed)

    output = tmp_path / "checkpoint.json"
    written = write_checkpoint(output)
    loaded = json.loads(output.read_text(encoding="utf-8"))
    assert loaded == written == first


def test_shortest_route_targets_real_multiband_family_not_duplicate_columns() -> None:
    checkpoint = build_checkpoint()
    routing = checkpoint["routing"]
    joined = " ".join(routing["shortest_next_closure"]).lower()
    forbidden = " ".join(routing["do_not_do"]).lower()
    assert "at least two active slow bands" in joined
    assert "rank-two coverage" in joined
    assert "agent 4" in joined
    assert "do not fabricate extra per-label free coefficients" in forbidden
    assert "do not call delta_band a second direction" in forbidden
