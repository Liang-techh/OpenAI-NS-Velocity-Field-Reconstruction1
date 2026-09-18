"""Agent-5 routing checkpoint for the bounded physical-family coordinate bridge.

This module adds no Kokuno mathematics.  It binds Agent 2 PR #400's explicit
repository coefficient coordinates to Agent 3 PR #401's unit-safe covariance
bridge and records the current Agent 1/4 sibling frontiers.  The one-band
calibration remains rank deficient, so this checkpoint deliberately does not
promote a correction, candidate artifact, export, or PDE-validation state.

The formal held-out gates remain normalized momentum max/L2 <= 1e-3 and
divergence max/L2 <= 1e-5.  No local covariance or coefficient-unit diagnostic
is interchangeable with that full-domain PDE gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-bounded-coordinate-routing-checkpoint-v22"
TASK_ID = "KOKUNO-A5-BOUNDED-COORDINATE-ROUTING-022"
BASE_AGENT3_PR = 401
BASE_AGENT3_HEAD = "c8db7fc707e74beef510cdd90da06b5cf2295682"

FIXED_GATES = {
    "held_out_normalized_full_momentum_max": 1.0e-3,
    "held_out_normalized_full_momentum_l2": 1.0e-3,
    "divergence_max": 1.0e-5,
    "divergence_l2": 1.0e-5,
    "changed": False,
}

ST006_REFERENCE = {
    "candidate_sha256": "6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3",
    "validation_seed": 9172801,
    "held_out_cartesian_points": 4096,
    "validation_times": 6,
    "finest_spatial_step": 0.005,
    "momentum_sampled_max": 0.1082289305112118,
    "volume_l2": 0.10758432876230622,
    "pde_validated": False,
}

AGENT2_COORDINATE_RECEIPT = {
    "task_id": "K2-OSC-023",
    "source_pr": 400,
    "source_head": "d4464acc3e7072f473746df9b02df823cdbcb8c4",
    "dedicated_run": 35345330212,
    "standard_run": 35345330333,
    "parameter_names": ["delta_common", "delta_band"],
    "parameter_unit": "dimensionless_fractional_multiplier_of_Q_scaled_physical_velocity",
    "parameter_count": 2,
    "aggregate_l1_bound_max": 0.5,
    "repository_coefficient_unit_mapping_available": True,
    "source_coefficient_unit_mapping_available": False,
    "one_band_band_coordinate_must_be_inactive": True,
    "actual_positive_order_background_binding_ready": False,
    "public_actual_source_xyz_t_velocity_ready": False,
}

AGENT3_BRIDGE_RECEIPT = {
    "task_id": "KOKUNO-A3-BOUNDED-COORDINATE-UNIT-BRIDGE-020",
    "source_pr": 401,
    "source_head": BASE_AGENT3_HEAD,
    "dedicated_run": 35346154634,
    "artifact_id": 10547465914,
    "artifact_digest": "sha256:84298de509d06cd6bd11a9e701069c3e5c6a3578be0ee0faf3726f5b17450dd6",
    "frozen_physical_amplitude_budget": 0.012513055889932317,
    "reference_physical_amplitude": 0.125,
    "converted_fractional_budget": 0.10010444711945854,
    "budget_conversion": "B_fractional=B_physical_amplitude/abs(reference_physical_amplitude)",
    "budget_changed": False,
    "budget_unit_matches_jacobian": True,
    "common_coordinate_relative_error_to_scaled_existing_response": 4.78894088298019e-16,
    "band_coordinate_response_vector_rms": 0.0,
    "one_band_band_coordinate_inactive": True,
    "nodes_requiring_second_direction": 25,
    "rank2_required_nodes": 0,
    "rank_deficient_required_nodes": 25,
    "algebraic_relative_stress_residual_rms": 0.2101945021227428,
    "max_additional_coefficient_l1_on_required": 0.0015400684172224392,
    "additional_family_l1_budget": 0.10010444711945854,
    "family_bounded_inverse_preflight_passed": False,
    "finite_correction_cycle_authorized": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
}

UPSTREAM = {
    "agent1": {
        "latest_pr": 399,
        "head_sha": "b0f68c45de12500533297fb36d9006fbba00ea94",
        "standard_run": 35344681320,
        "eta_smooth_repair_coefficient_family_reconstructed": True,
        "repair_applied_to_selected_autonomous_modulation": True,
        "source_admissible_loop_reconstructed": False,
        "global_pressure_matched": False,
        "global_leading_profile_reconstructed": False,
        "leading_ready": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent2": {
        "latest_pr": 400,
        "head_sha": AGENT2_COORDINATE_RECEIPT["source_head"],
        "dedicated_run": AGENT2_COORDINATE_RECEIPT["dedicated_run"],
        "standard_run": AGENT2_COORDINATE_RECEIPT["standard_run"],
        "repository_coefficient_unit_mapping_available": True,
        "source_coefficient_unit_mapping_available": False,
        "bounded_coordinate_contract_ready": True,
        "public_actual_source_xyz_t_velocity_ready": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "consumed_in_executable_ancestry": True,
    },
    "agent3": {
        "latest_pr": BASE_AGENT3_PR,
        "head_sha": BASE_AGENT3_HEAD,
        "dedicated_run": AGENT3_BRIDGE_RECEIPT["dedicated_run"],
        "artifact_id": AGENT3_BRIDGE_RECEIPT["artifact_id"],
        "artifact_digest": AGENT3_BRIDGE_RECEIPT["artifact_digest"],
        "bounded_coordinate_unit_bridge_ready": True,
        "nodes_requiring_second_direction": 25,
        "rank_two_coverage": "0/25",
        "missing_relative_stress_response": 0.2101945021227428,
        "actual_multiband_family_consumed": False,
        "family_bounded_inverse_preflight_passed": False,
        "finite_correction_cycle_rerun_allowed": False,
        "consumed_in_executable_ancestry": True,
    },
    "agent4": {
        "latest_pr": 394,
        "head_sha": "69a7ca724db377510815eddcc2b152b05d7193ac",
        "dedicated_run": 35341494077,
        "standard_run": 35341494053,
        "artifact_id": 10545247621,
        "artifact_digest": "sha256:6b2cc78a2f0ba6977b938fdafb6aa62a4059964727e14ad9fa681c7246c09765",
        "source_compatible_partition_independent_preflight_passed": True,
        "bounded_coordinate_multiband_independent_audit_run": False,
        "consumed_in_executable_ancestry": False,
    },
    "previous_agent5": {
        "latest_pr": 395,
        "head_sha": "63de571a12e8fd80a47e0bf0924ce1412c1f76ce",
        "source_compatible_partition_routed": True,
        "consumed_in_executable_ancestry": False,
    },
}

AUTONOMOUS_ROUTE_POLICY = {
    "kokuno_derived_autonomous_choices_allowed": True,
    "repository_coefficient_coordinates_are_autonomous": True,
    "repository_coefficient_coordinates_may_be_called_source_recovered": False,
    "source_compatible_partition_may_be_used_as_scaffold": True,
    "selected_agent1_modulation_may_be_used_as_scaffold": True,
    "paper_exact": False,
}

FUTURE_CANDIDATE_ARTIFACT_CONTRACT = {
    "schema_name": "kokuno-candidate-artifact-v1-reserved",
    "status": "reserved_not_instantiated",
    "required_provenance": [
        "source_repository",
        "source_commit_or_version",
        "source_document",
        "source_derived_formulas",
        "autonomous_choices",
        "stage_artifact_hashes",
    ],
    "required_public_api": {
        "velocity": "velocity(points, t) -> (..., 3)",
        "pressure": "pressure(points, t) -> (...)",
        "forcing": "forcing(points, t) -> (..., 3) under the fixed/restricted contract only",
        "save_load": "deterministic metadata plus hashed numeric payloads",
    },
    "required_stage_receipts": [
        "source_profile_ingest",
        "leading",
        "oscillatory",
        "mean_radial_correction",
        "independent_validation",
    ],
}

STATES = {
    "source_compatible_partition_realization_ready": True,
    "source_compatible_partition_independently_audited": True,
    "repository_coefficient_unit_mapping_available": True,
    "source_coefficient_unit_mapping_available": False,
    "leading_ready": False,
    "oscillatory_machinery_ready": True,
    "oscillatory_ready": False,
    "public_kokuno_derived_oscillatory_xyz_t_velocity_ready": False,
    "genuinely_independent_second_covariance_column_ready": False,
    "correction_ready": False,
    "candidate_artifact_instantiated": False,
    "complete_kokuno_composite_velocity_ready": False,
    "velocity_export_ready": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
}

PIPELINE_FRONTIER = [
    {
        "stage": "source_profile_ingest",
        "ready": False,
        "status": "autonomous_partition_and_repository_coefficient_units_ready",
        "blockers": [
            "actual positive-order/background path or an explicitly provenance-labelled autonomous substitute",
            "actual multi-band family with at least two active slow bands",
        ],
    },
    {
        "stage": "leading_candidate",
        "ready": False,
        "status": "eta_smooth_i1_repair_ready_global_assembly_missing",
        "evidence": [
            "Agent-1 #399 supplies an eta-smooth PA.17 repair family for the selected autonomous modulation"
        ],
        "blockers": [
            "complete radial/heat overlay composition",
            "global matched pressure",
            "admissible fixed/restricted forcing contract on the assembled leading field",
        ],
    },
    {
        "stage": "oscillatory_augmentation",
        "ready": False,
        "status": "coordinate_units_ready_waiting_for_multiband_physical_family",
        "evidence": [
            "Agent-2 #400 exposes delta_common/delta_band with explicit dimensionless fractional units"
        ],
        "blockers": [
            "positive-order/background and mode binding",
            "public Q-scaled xyz,t velocity by beta and total",
            "two or more active slow bands so delta_band is nonzero",
            "genuinely independent covariance direction",
        ],
    },
    {
        "stage": "mean_radial_corrections",
        "ready": False,
        "status": "unit_bridge_cleared_but_one_band_rank_deficient",
        "evidence": [
            "Agent-3 #401 converts the frozen physical-amplitude budget to fractional units without threshold laundering",
            "common-coordinate calibration relative error is below 5e-16",
            "the one-band delta_band response is exactly inactive",
            "rank-two coverage remains 0/25 and missing relative stress remains 0.2101945021",
        ],
        "blockers": [
            "actual multi-band physical family",
            "rank-two coverage under the frozen bounded inverse",
            "spacetime/radial guards including retained radial d_z sigma_1",
            "Agent-4 independent audit before materialized correction promotion",
        ],
    },
    {"stage": "candidate_artifact", "ready": False, "status": "waiting_for_leading_oscillatory_correction"},
    {"stage": "independent_full_pde_validation", "ready": False, "status": "waiting_for_frozen_complete_candidate"},
    {"stage": "report_and_export", "ready": False, "status": "waiting_for_validation_receipt"},
]


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def checkpoint_sha256(payload: dict[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def build_checkpoint() -> dict[str, Any]:
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base": {"agent3_pr": BASE_AGENT3_PR, "agent3_head": BASE_AGENT3_HEAD},
        "fixed_gates": FIXED_GATES,
        "st006_reference": ST006_REFERENCE,
        "agent2_coordinate_receipt": AGENT2_COORDINATE_RECEIPT,
        "agent3_bridge_receipt": AGENT3_BRIDGE_RECEIPT,
        "upstream": UPSTREAM,
        "autonomous_route_policy": AUTONOMOUS_ROUTE_POLICY,
        "future_candidate_artifact_contract": FUTURE_CANDIDATE_ARTIFACT_CONTRACT,
        "pipeline_frontier": PIPELINE_FRONTIER,
        "states": STATES,
        "routing": {
            "new_fact": (
                "Agent-2 #400 and Agent-3 #401 close the repository coefficient-unit seam: the frozen physical-amplitude budget is converted explicitly into fractional coordinates, but the real one-band control remains rank deficient."
            ),
            "interpretation": (
                "Unit compatibility is no longer the correction-lane blocker.  The next missing mathematical object is a real multi-band physical oscillatory family whose band coordinate produces a genuinely independent covariance direction."
            ),
            "do_not_do": [
                "do not reuse 0.012513055889932317 as a dimensionless fractional budget",
                "do not call delta_band a second direction when only one slow band is active",
                "do not fabricate extra per-label free coefficients to create rank",
                "do not rerun the finite correction cycle while rank-two coverage is 0/25",
                "do not instantiate/export the final Kokuno candidate before leading, oscillatory and correction stages coexist",
                "do not weaken the fixed 1e-3 momentum or 1e-5 divergence gates",
                "do not use residual-defined free forcing",
            ],
            "shortest_next_closure": [
                "Agent 2: instantiate the positive-order/background path on the independently cleared source-compatible partition and produce a public Q-scaled by-beta/total xyz,t family with at least two active slow bands; keep source/autonomous provenance explicit.",
                "Agent 3: consume that exact multi-band family through the #401 unit bridge and #393 bounded inverse; require rank-two coverage, algebraic fit, aggregate budget, spacetime and radial guards before any correction materialization.",
                "Agent 4: independently audit the first passing multi-band tangent/rank receipt and later the frozen full candidate; only Agent 4 may promote the formal held-out PDE gate.",
                "Agent 1: compose #399's eta-smooth I1 repair with the remaining overlays and publish global leading velocity, matched pressure and admissible fixed/restricted forcing.",
                "Agent 5: instantiate the deterministic candidate artifact and Python/MATLAB smoke only after leading + oscillatory + correction coexist in one executable ancestry.",
            ],
        },
    }
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    return payload


def write_checkpoint(path: str | Path) -> dict[str, Any]:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = build_checkpoint()
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent5/bounded_coordinate_routing_checkpoint_v22.json",
    )
    args = parser.parse_args(argv)
    payload = write_checkpoint(args.output)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
