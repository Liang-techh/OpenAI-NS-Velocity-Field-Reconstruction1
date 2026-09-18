from __future__ import annotations

import json

import pytest

from openai_ns_reconstruction.kokuno_audited_geometry_leading_routing_checkpoint import (
    AGENT1_RECEIPT,
    AGENT4_RECEIPT,
    FORMAL_GATES,
    PREVIOUS_AGENT5_RECEIPT,
    build_checkpoint,
    checkpoint_sha256,
    validate_checkpoint,
    write_checkpoint,
)


def _rehash(payload):
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    return payload


def test_v34_routes_audited_geometry_and_leading_obstruction_without_promotion() -> None:
    payload = build_checkpoint()
    validate_checkpoint(payload)
    states = payload["states"]

    assert states["autonomous_signed_rectangle_geometry_ready"]
    assert states["autonomous_signed_rectangle_geometry_independently_audited"]
    assert states["selected_shared_C_pa10_path_obstructed"]
    assert not states["selected_shared_C_pa16_handoff_allowed"]
    assert states["leading_reparameterization_required"]
    assert states["candidate_autonomous_realization_route_open"]

    for name in (
        "leading_ready",
        "public_provenance_labelled_xyz_t_oscillatory_velocity_ready",
        "oscillatory_ready",
        "genuinely_independent_second_covariance_column_ready",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        assert not states[name]


def test_agent1_obstruction_receipt_is_exact_and_fail_closed() -> None:
    payload = build_checkpoint()
    a1 = payload["upstream"]["agent1_latest_sibling"]

    assert a1 == AGENT1_RECEIPT
    assert a1["head"] == "a0037e5cff60be97aa01cc7273eb7af18389ed7a"
    assert a1["dedicated_run"] == 35404218836
    assert a1["dedicated_status"] == "success"
    assert a1["standard_status_at_freeze"] == "in_progress"
    assert a1["selected_shared_C_pa10_obstructed"]
    assert not a1["selected_pa16_handoff_allowed"]
    assert a1["necessary_T_sh_lower_bound"] > a1["geometric_T_sh_ceiling"]
    assert a1["separation_margin_upper_bound"] < 0.0
    assert not a1["source_B0_analytic_bound_proved"]
    assert not a1["global_pressure_matched"]


def test_agent4_autonomous_geometry_audit_is_independent_not_source_exact() -> None:
    payload = build_checkpoint()
    a4 = payload["upstream"]["agent4_sibling"]

    assert a4 == AGENT4_RECEIPT
    assert a4["head"] == "aa9d2bb0a6c196ea86daffd5a6f539785a9b0363"
    assert a4["dedicated_run"] == 35403144495
    assert a4["standard_run"] == 35403144460
    assert a4["dedicated_status"] == "success"
    assert a4["standard_status"] == "success"
    assert a4["evidence_class"] == "independent_validation"
    assert a4["autonomous_signed_rectangle_geometry_independently_audited"]
    assert a4["local_structural_preflight_passed"]
    assert a4["schedule_relative_max"] == 0.0
    assert a4["permutation_relative_max"] == 0.0
    assert a4["strict_geometry_guard_violations"] == 0
    assert a4["caller_order_assignment_mutation_mismatch_fraction"] == pytest.approx(0.8333333333333334)
    assert a4["wrong_Ls_formula_mutation_relative_rms"] == pytest.approx(0.999998125541105)
    assert not a4["source_rectangle_centers_recovered"]
    assert not a4["source_rectangle_radius_r0_recovered"]
    assert not a4["formal_full_domain_pde_gate_assessed"]
    assert not a4["pde_validated"]


def test_previous_agent5_and_formal_gates_remain_frozen() -> None:
    payload = build_checkpoint()
    previous = payload["upstream"]["previous_agent5_ancestry"]

    assert previous == PREVIOUS_AGENT5_RECEIPT
    assert previous["head"] == "96e9514766ca5c512dc30409fbfa70b732fb56e7"
    assert previous["dedicated_status"] == "success"
    assert previous["standard_status"] == "success"

    assert payload["formal_gates"] == FORMAL_GATES
    assert payload["formal_gates"]["held_out_normalized_momentum_max"] == 1.0e-3
    assert payload["formal_gates"]["held_out_normalized_momentum_l2"] == 1.0e-3
    assert payload["formal_gates"]["held_out_divergence_max"] == 1.0e-5
    baseline = payload["baseline_vs_kokuno"]["st006"]
    assert baseline["momentum_sampled_max"] == pytest.approx(0.1082289305112118)
    assert baseline["momentum_volume_l2"] == pytest.approx(0.10758432876230622)
    assert payload["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is None


def test_v34_fails_closed_on_obstruction_audit_or_pde_overpromotion() -> None:
    payload = json.loads(json.dumps(build_checkpoint()))
    payload["states"]["selected_shared_C_pa16_handoff_allowed"] = True
    with pytest.raises(ValueError, match="PA.16"):
        validate_checkpoint(_rehash(payload))

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["states"]["leading_ready"] = True
    with pytest.raises(ValueError):
        validate_checkpoint(_rehash(payload))

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["upstream"]["agent4_sibling"]["source_rectangle_centers_recovered"] = True
    with pytest.raises(ValueError):
        validate_checkpoint(_rehash(payload))

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["upstream"]["agent1_latest_sibling"]["source_B0_analytic_bound_proved"] = True
    with pytest.raises(ValueError):
        validate_checkpoint(_rehash(payload))

    payload = json.loads(json.dumps(build_checkpoint()))
    payload["states"]["pde_validated"] = True
    with pytest.raises(ValueError):
        validate_checkpoint(_rehash(payload))


def test_v34_roundtrip_is_deterministic(tmp_path) -> None:
    payload = build_checkpoint()
    path = tmp_path / "checkpoint.json"
    written = write_checkpoint(path)
    loaded = json.loads(path.read_text(encoding="utf-8"))

    assert loaded == written
    assert loaded["checkpoint_sha256"] == checkpoint_sha256(loaded)

    loaded["formal_gates"]["held_out_normalized_momentum_max"] = 2.0e-3
    with pytest.raises(ValueError, match="checkpoint sha256 mismatch"):
        validate_checkpoint(loaded)
