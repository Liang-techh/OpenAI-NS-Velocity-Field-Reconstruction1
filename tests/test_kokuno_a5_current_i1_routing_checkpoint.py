from __future__ import annotations

import copy
import inspect

import pytest

from openai_ns_reconstruction import kokuno_a4_current_i1_leading_divergence_independent_audit as a4
from openai_ns_reconstruction import kokuno_a5_current_i1_routing_checkpoint as a5


def _valid_a4_receipt(*, audit_pass: bool = True) -> dict[str, object]:
    resolutions = []
    for h in a4.SPATIAL_STEPS:
        resolutions.append(
            {
                "step": float(h),
                "sampled_max": 4.0e-7,
                "estimated_volume": 1.0,
                "pooled_volume_l2_estimate": 3.0e-7,
                "pooled_weighted_rms": 3.0e-7,
                "normalized_sampled_max": 4.0e-7,
                "normalized_weighted_rms": 3.0e-7,
                "i1_entry_seam_sampled_max": 5.0e-7,
                "axis_axis_near_sampled_max": 5.0e-7,
                "speed_rms": 1.0e-2,
                "per_region": {},
                "per_time": [],
                "worst_witness": {},
            }
        )
    return {
        "schema": a4.SCHEMA,
        "upstream_pr": 1051,
        "upstream_head": a5.AGENT1_CURRENT_I1["head"],
        "candidate_semantic_sha256": "a" * 64,
        "protocol": {
            "seed": a4.SEED,
            "spatial_steps": list(a4.SPATIAL_STEPS),
            "derivative_operator": a5.AGENT4_CURRENT_I1_AUDIT["operator"],
            "frozen_i1_closure_tolerance": 5.0e-7,
        },
        "gates": {
            "final_project_momentum_gate_unchanged": 1.0e-3,
            "final_project_divergence_gate_unchanged": 1.0e-5,
        },
        "resolutions": resolutions,
        "checks": {},
        "truth_boundary": {
            "current_i1_leading_velocity_consumed": True,
            "leading_only_ns_residual_assessed": False,
            "leading_plus_oscillatory_ns_residual_assessed": False,
            "after_correction_ns_residual_assessed": False,
            "matched_global_pressure_materialized": False,
            "restricted_forcing_materialized": False,
            "canonical_24_48_96_volume_admission_assessed": False,
            "same_protocol_comparable_to_st006": False,
            "pde_validated": False,
        },
        "audit_pass": audit_pass,
    }


def _redigest(checkpoint: dict[str, object]) -> None:
    payload = {
        k: copy.deepcopy(v)
        for k, v in checkpoint.items()
        if k not in {"schema", "task", "digest"}
    }
    checkpoint["digest"] = a5._sha256(payload)


def test_build_checkpoint_binds_exact_current_i1_and_preserves_truth_boundary() -> None:
    checkpoint = a5.build_checkpoint(_valid_a4_receipt())
    a5.validate_checkpoint(checkpoint)
    a5.enforce_scoped_a4_gate(checkpoint)

    assert checkpoint["agent1_current_i1"]["head"] == "ea59dc305b4265fcb0bb0f948c2a45f4477d1ae5"
    assert checkpoint["agent4_current_i1_audit"]["head"] == "ad9f8722b6ec89df975f9b2402e63066dc14422d"
    assert checkpoint["agent2_latest_sibling"]["consumed_by_this_checkpoint"] is False
    assert checkpoint["agent3_latest_sibling"]["consumed_by_current_i1_identity"] is False
    assert checkpoint["prior_a5_unmodulated_x4"]["evidence_transfer_to_current_i1_allowed"] is False
    assert checkpoint["readiness"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert checkpoint["truth_boundary"]["heldout_normalized_full_ns_residual_assessed"] is False
    assert checkpoint["truth_boundary"]["pde_validated"] is False
    assert checkpoint["final_gates"]["normalized_momentum_sampled_max"] == 1.0e-3
    assert checkpoint["final_gates"]["normalized_divergence_sampled_max"] == 1.0e-5
    assert checkpoint["final_gates"]["residual_defined_free_forcing_allowed"] is False


def test_scoped_negative_a4_receipt_is_saved_but_not_admitted() -> None:
    checkpoint = a5.build_checkpoint(_valid_a4_receipt(audit_pass=False))
    a5.validate_checkpoint(checkpoint)
    with pytest.raises(AssertionError, match="scoped divergence audit did not pass"):
        a5.enforce_scoped_a4_gate(checkpoint)
    assert checkpoint["truth_boundary"]["scientific_admission"] is False
    assert checkpoint["truth_boundary"]["pde_validated"] is False


def test_cross_lineage_and_project_gate_mutations_fail_closed() -> None:
    base = a5.build_checkpoint(_valid_a4_receipt())

    mutated = copy.deepcopy(base)
    mutated["prior_a5_unmodulated_x4"]["evidence_transfer_to_current_i1_allowed"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="prior A5 X4 lineage"):
        a5.validate_checkpoint(mutated)

    mutated = copy.deepcopy(base)
    mutated["agent2_latest_sibling"]["self_contained_cartesian_composite_on_agent1_1051"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="A2 sibling"):
        a5.validate_checkpoint(mutated)

    mutated = copy.deepcopy(base)
    mutated["truth_boundary"]["pde_validated"] = True
    _redigest(mutated)
    with pytest.raises(ValueError, match="truth boundary drifted"):
        a5.validate_checkpoint(mutated)

    mutated = copy.deepcopy(base)
    mutated["final_gates"]["normalized_momentum_sampled_max"] = 2.0e-3
    _redigest(mutated)
    with pytest.raises(ValueError, match="fixed project gates drifted"):
        a5.validate_checkpoint(mutated)


def test_upstream_receipt_rejects_pde_laundering_and_public_surface_has_no_knobs() -> None:
    bad = _valid_a4_receipt()
    bad["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="illegally promoted pde_validated"):
        a5.validate_upstream_a4_receipt(bad)

    params = inspect.signature(a5.build_checkpoint).parameters
    forbidden = {"threshold", "momentum_gate", "divergence_gate", "forcing", "residual", "correction"}
    assert forbidden.isdisjoint(params)
