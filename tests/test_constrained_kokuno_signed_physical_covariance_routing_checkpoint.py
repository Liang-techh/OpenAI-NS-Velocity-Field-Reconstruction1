from __future__ import annotations

import json

import pytest

from openai_ns_reconstruction.kokuno_signed_physical_covariance_routing_checkpoint import (
    FORMAL_GATES,
    build_checkpoint,
    checkpoint_sha256,
    validate_checkpoint,
    write_checkpoint,
)


def test_checkpoint_promotes_only_supplied_signed_physical_preflight() -> None:
    payload = build_checkpoint()
    validate_checkpoint(payload)
    states = payload["states"]

    assert states["displayed_source_signed_covariance_reference_executable"]
    assert states["signed_covariance_reference_independently_audited"]
    assert states["signed_covariance_physical_to_reference_unit_bridge_ready"]
    assert states["supplied_signed_complete_curl_physical_family_executable"]
    assert states["supplied_signed_complete_curl_physical_covariance_rank_two"]
    assert states["supplied_signed_complete_curl_physical_covariance_independently_audited"]

    assert not states["actual_source_physical_covariance_rank_two_assessed"]
    assert not states["genuinely_independent_second_covariance_column_ready"]
    assert not states["oscillatory_ready"]
    assert not states["correction_ready"]
    assert not states["candidate_artifact_instantiated"]
    assert not states["pde_validated"]


def test_checkpoint_preserves_agent4_black_box_receipt_and_frozen_guards() -> None:
    payload = build_checkpoint()
    a4 = payload["upstream"]["agent4"]

    assert a4["head"] == "44793b7f632bb0344bbd24d2a6074828b263b04b"
    assert a4["oracle_consumes_public_cartesian_total_only"]
    assert not a4["oracle_reuses_complete_curl_helper"]
    assert not a4["oracle_reads_returned_tangents_or_amplitudes"]
    assert a4["local_physical_covariance_preflight_passed"]
    assert a4["minimum_physical_rank_ratio"] == pytest.approx(0.35573170032291024)
    assert a4["maximum_resolution_relative_drift"] == pytest.approx(4.5397547533405476e-13)
    assert a4["minimum_velocity_rms"] == pytest.approx(21.242258025207185)
    assert a4["maximum_duplicated_sign_mutation_rank_ratio"] == pytest.approx(4.565665180524497e-13)
    assert a4["minimum_physical_rank_ratio"] >= a4["frozen_minimum_rank_ratio_guard"]
    assert a4["maximum_resolution_relative_drift"] <= a4["frozen_maximum_resolution_drift_guard"]
    assert a4["maximum_duplicated_sign_mutation_rank_ratio"] <= a4["frozen_maximum_duplicated_sign_rank_ratio_guard"]
    assert not a4["actual_source_physical_covariance_rank_two_assessed"]


def test_checkpoint_records_agent3_cross_check_without_promoting_real_candidate() -> None:
    payload = build_checkpoint()
    a3 = payload["upstream"]["agent3"]

    assert a3["rank_two_cells"] == 5
    assert a3["radial_cells"] == 5
    assert a3["minimum_covariance_rank_ratio"] == pytest.approx(0.3203750745181774)
    assert a3["maximum_quadratic_homogeneity_error"] < 1.0e-12
    assert a3["physical_to_reference_rule"] == "H_ref * y = DeltaC / epsilon"
    assert not a3["real_candidate_defect_consumed"]
    assert not a3["genuinely_independent_second_covariance_column_ready"]


def test_checkpoint_keeps_source_binding_and_formal_pde_gate_fail_closed() -> None:
    payload = build_checkpoint()
    states = payload["states"]
    for name in (
        "actual_positive_order_background_bound",
        "actual_source_h_sigma_pulse_integrals_bound",
        "actual_signed_auxiliary_rectangles_or_modes_bound",
        "actual_source_partition_labels_instantiated",
        "public_source_bound_xyz_t_velocity_ready",
        "real_candidate_defect_consumed",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        assert not states[name]

    assert payload["formal_gates"] == FORMAL_GATES
    assert payload["formal_gates"]["held_out_normalized_momentum_max"] == 1.0e-3
    assert payload["formal_gates"]["held_out_divergence_max"] == 1.0e-5
    baseline = payload["baseline_vs_kokuno"]["st006"]
    assert baseline["momentum_sampled_max"] == pytest.approx(0.1082289305112118)
    assert baseline["momentum_volume_l2"] == pytest.approx(0.10758432876230622)
    assert payload["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is None


def test_checkpoint_fails_closed_on_source_or_real_candidate_promotion() -> None:
    payload = build_checkpoint()
    payload["states"]["genuinely_independent_second_covariance_column_ready"] = True
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    with pytest.raises(ValueError, match="fail-closed state promoted"):
        validate_checkpoint(payload)

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["upstream"]["agent4"]["actual_source_physical_covariance_rank_two_assessed"] = True
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    with pytest.raises(ValueError, match="promoted to actual-source evidence"):
        validate_checkpoint(payload)

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["upstream"]["agent2"]["actual_source_h_sigma_pulse_integrals_bound"] = True
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    with pytest.raises(ValueError, match="manufactured/source truth boundary changed"):
        validate_checkpoint(payload)


def test_checkpoint_sha_roundtrip_and_receipt_copy_isolation(tmp_path) -> None:
    first = build_checkpoint()
    first["upstream"]["agent4"]["minimum_physical_rank_ratio"] = -1.0
    second = build_checkpoint()
    assert second["upstream"]["agent4"]["minimum_physical_rank_ratio"] > 0.0

    path = tmp_path / "checkpoint.json"
    payload = write_checkpoint(path)
    loaded = json.loads(path.read_text(encoding="utf-8"))
    assert loaded == payload
    assert loaded["checkpoint_sha256"] == checkpoint_sha256(loaded)

    loaded["formal_gates"]["held_out_normalized_momentum_max"] = 2.0e-3
    with pytest.raises(ValueError, match="checkpoint sha256 mismatch"):
        validate_checkpoint(loaded)
