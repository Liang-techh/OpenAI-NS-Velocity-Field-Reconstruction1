"""Kokuno Agent-5 integration contract for the current RF40 lambda-turn seam.

Integration/provenance glue only.  This increment registers the first matched
current lineage through the RF40 lambda-turn endpoint ``X_3=X_w``:

    A1 #998 leading through X_3
      -> A2 #1004 identity-bound leading + frozen complete-curl oscillation
      -> A4 #1007 implementation-distinct public-Cartesian divergence audit.

A1 #1005 (leading-only RF40 power law through X_4) and A3 #1006
(first-turn radial force through X_1) are recorded only as non-consumed sibling
state.  Scoped divergence is not a complete NS residual, unresolved CI is not
scientific admission, and no correction/PDE state is promoted here.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-current-rf40-lambda-turn-composite-ingest-v1"
TASK_ID = "KOKUNO-A5-CURRENT-RF40-LAMBDA-TURN-COMPOSITE-INGEST-099"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 1002,
    "head": "bb3d943461820fe43c0bdb2063cb11b428c4ea1c",
    "branch": "codex/kokuno-a5-rf40-axial-shutdown-composite-ingest-098",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_current_rf40_axial_shutdown_composite_ingest_contract.py",
    "source_blob": "6ac6d7c41f921b52e8bc32f6ffcde5cd37602360",
    "test_path": "tests/test_kokuno_a5_current_rf40_axial_shutdown_composite_ingest_contract.py",
    "test_blob": "1efd89b81ca785da555bfdcfc36b4f67ede7285c",
    "workflow_path": ".github/workflows/kokuno-agent5-current-rf40-axial-shutdown-composite-ingest.yml",
    "workflow_blob": "228104b1ec697c8cfff0ebc8a3b8759cfb08c1b8",
}

AGENT2_LAMBDA_TURN_COMPOSITE = {
    "pr": 1004,
    "head": "e9697cdfe2249e1b8bc2f9241bd16d456c60418f",
    "branch": "codex/kokuno-a2-rf40-lambda-turn-composite-086",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_rf40_lambda_turn_leading_oscillatory_identity.py",
    "source_blob": "def2998256b2339115c2400c75b6db6ab58dfe95",
    "test_path": "tests/test_constrained_kokuno_current_rf40_lambda_turn_leading_oscillatory_identity.py",
    "test_blob": "f8880d7ead60a89be60803435c106c7b075c0f2b",
    "workflow_path": ".github/workflows/kokuno-agent2-current-rf40-lambda-turn-leading-oscillatory.yml",
    "workflow_blob": "3dcf132cd9fe8b9101a7e5ca6058c3b44b4a8a34",
    "leading_pr": 998,
    "leading_head": "43b295444b1e9558222d757cb551e04385eddcb5",
    "velocity_api": "velocity(x,y,z,t)->[...,3]",
    "materialized_domain": "current leading + frozen complete-curl oscillation through RF40 lambda turn X_3=X_w",
    "identity_preserving_save_load": True,
    "full_concrete_oscillatory_runtime_digest_bound": True,
    "corrected_source_provenance_bound_in_agent2_identity": True,
    "fails_closed_after_X_3": True,
    "power_law_composite_materialized": False,
    "pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "complete_ns_residual_assessed": False,
}

AGENT4_LAMBDA_TURN_AUDIT = {
    "pr": 1007,
    "head": "8b5b337940edd568bf30ddba7afe68815bce1fd9",
    "branch": "codex/kokuno-a4-rf40-lambda-turn-composite-divergence-audit-101",
    "audited_agent2_pr": 1004,
    "audited_agent2_head": "e9697cdfe2249e1b8bc2f9241bd16d456c60418f",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_rf40_lambda_turn_leading_oscillatory_divergence_independent_audit.py",
    "source_blob": "2ae59ab9bc8480979f0b0660b1db8b801fe847c1",
    "test_path": "tests/test_constrained_kokuno_a4_current_rf40_lambda_turn_leading_oscillatory_divergence_independent_audit.py",
    "test_blob": "813dfc007f248e5297829330df243df43db59ed7",
    "workflow_path": ".github/workflows/kokuno-agent4-current-rf40-lambda-turn-leading-oscillatory-divergence-audit.yml",
    "workflow_blob": "0fc520a698d94026ca82738dc4ab5f57252c52da",
    "reference_path": "save/reload A2 #1004, then centered FD2 of public Cartesian velocity only",
    "seed": 9173731,
    "times": [0.31, 0.47, 0.63, 0.71],
    "fd2_step_ladder": [0.02, 0.01, 0.005],
    "lambda_turn_rho_zones": [[0.08, 0.28], [0.38, 0.62], [0.72, 0.92]],
    "strict_lambda_turn_heldout_points": 72,
    "scoped_divergence_gate": 1.0e-5,
    "uses_agent2_production_jacobian_or_divergence": False,
    "uses_agent1_source_derivative_helpers": False,
    "uses_pressure_forcing_or_residual_helper": False,
    "canonical_24_48_96_whole_domain_admission": False,
    "momentum_or_full_ns_evidence": False,
    "scientific_admission": False,
}

AGENT1_POWER_LAW_SIBLING = {
    "pr": 1005,
    "head": "2c76ebdc41d6c566f43a1305034ba2ff9dce410b",
    "branch": "agent-kokuno-1/current-cartesian-rf40-power-law-100",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_rf40_power_law.py",
    "source_blob": "36a0183352a2005bd0b5f477594e05a2e4fe3868",
    "test_path": "tests/test_kokuno_pa16_current_cartesian_rf40_power_law.py",
    "test_blob": "9c88ab44c3a195b42f6f600b54563491869be64f",
    "workflow_path": ".github/workflows/kokuno-agent1-current-cartesian-rf40-power-law.yml",
    "workflow_blob": "2f871a63034ca0284761e97cdd014766c9820d94",
    "materialized_domain": "leading-only through RF40 constant-lambda power-law endpoint X_4",
    "consumed_by_agent2_1004": False,
    "covered_by_agent4_1007": False,
    "consumed_by_this_increment": False,
}

AGENT3_FIRST_TURN_RADIAL_FORCE_SIBLING = {
    "pr": 1006,
    "head": "2c10faab7c229b801cb868355e0a1de0a729445c",
    "branch": "codex/kokuno-a3-rf40-first-turn-radial-force-108",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_rf40_first_turn_nonlinear_radial_force.py",
    "source_blob": "d73a1d6e6b7c0cc341dd99f92d6f20d6eb53acff",
    "test_path": "tests/test_constrained_kokuno_current_rf40_first_turn_nonlinear_radial_force.py",
    "test_blob": "a52ae206c588f89690fed5aa2bd2d5dfc416feca",
    "workflow_path": ".github/workflows/kokuno-agent3-current-rf40-first-turn-nonlinear-radial-force.yml",
    "workflow_blob": "c974546fc915b731801d620ae4f00d49b8cb4c1c",
    "materialized_domain": "nonlinear radial force from axial stress through RF40 first turn X_1 only",
    "rf40_axial_shutdown_current_lineage_consumed": False,
    "rf40_lambda_turn_current_lineage_consumed": False,
    "rf40_power_law_current_lineage_consumed": False,
    "authorized_as_correction_target": False,
    "real_cartesian_correction_velocity": False,
    "consumed_by_this_increment": False,
}

OBSERVED_CI_AT_FREEZE = {
    "parent_a5_1002": {
        "repository_tests_run": 35603144952,
        "dedicated_run": 35603145112,
        "status": "queued",
        "conclusion": None,
    },
    "agent2_1004_lambda_turn_composite": {
        "repository_tests_run": 35606610366,
        "dedicated_run": 35606610315,
        "status": "queued",
        "conclusion": None,
    },
    "agent1_1005_power_law_sibling": {
        "repository_tests_run": 35607090902,
        "dedicated_run": 35607090935,
        "status": "queued",
        "conclusion": None,
    },
    "agent3_1006_first_turn_radial_force_sibling": {
        "repository_tests_run": 35607681349,
        "dedicated_run": 35607681863,
        "status": "queued",
        "conclusion": None,
    },
    "agent4_1007_lambda_turn_audit": {
        "repository_tests_run": 35608063788,
        "dedicated_run": 35608063570,
        "status": "queued",
        "conclusion": None,
    },
}

FROZEN_SCIENCE = {
    "viscosity": 0.01,
    "evaluation_box": [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
    "time_interval": [0.25, 0.75],
    "forcing_family": "restricted two-parameter curl forcing",
    "residual_defined_free_forcing_forbidden": True,
}

FINAL_GATE = {
    "normalized_momentum_sampled_max": 1.0e-3,
    "normalized_momentum_volume_l2": 1.0e-3,
    "divergence_sampled_max": 1.0e-5,
    "divergence_volume_l2": 1.0e-5,
    "canonical_volume_quadrature_ladder": [24, 48, 96],
}

ST006_BASELINE = {
    "same_protocol_full_pde_baseline": True,
    "momentum_sampled_max": 0.1082289305112118,
    "momentum_volume_l2": 0.10758432876230622,
}

READINESS = {
    "leading_ready": False,
    "oscillatory_ready": True,
    "correction_ready": False,
    "velocity_export_ready": False,
    "pde_validated": False,
}

TRUTH_BOUNDARY = {
    "current_leading_plus_oscillatory_velocity_through_rf40_first_turn_materialized": True,
    "current_leading_plus_oscillatory_velocity_through_rf40_axial_shutdown_materialized": True,
    "current_leading_plus_oscillatory_velocity_through_rf40_lambda_turn_materialized": True,
    "identity_preserving_lambda_turn_composite_save_load_available": True,
    "agent4_1007_independent_lambda_turn_composite_divergence_audit_registered": True,
    "agent4_1007_independent_lambda_turn_composite_divergence_audit_admitted": False,
    "agent1_1005_power_law_leading_registered_as_nonconsumed_sibling": True,
    "agent3_1006_first_turn_radial_force_registered_as_nonconsumed_sibling": True,
    "current_nonlinear_mean_through_axial_shutdown_materialized": False,
    "current_nonlinear_mean_through_lambda_turn_materialized": False,
    "current_nonlinear_radial_stress_through_lambda_turn_materialized": False,
    "current_nonlinear_radial_force_through_lambda_turn_materialized": False,
    "leading_plus_oscillatory_velocity_through_power_law_materialized": False,
    "independent_power_law_composite_audit_available": False,
    "complete_post_xr_rf40_current_lineage_materialized": False,
    "cone_i1_i2_i3_i4_outer_overlays_completed": False,
    "outer_global_leading_velocity_materialized": False,
    "global_compact_support_completed": False,
    "matched_cartesian_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
    "restricted_forcing_materialized": False,
    "complete_ns_defect_materialized": False,
    "current_nonlinear_mean_authorized_as_correction_target": False,
    "real_agent3_ns_correction_velocity_materialized": False,
    "real_candidate_finite_correction_cycle_run": False,
    "complete_candidate_api_ready": False,
    "canonical_whole_domain_divergence_l2_assessed": False,
    "heldout_normalized_ns_residual_assessed": False,
    "same_protocol_st006_comparison_available_now": False,
    "residual_reduction_claimed": False,
    "scientific_admission": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}

PIPELINE_POSITION = {
    "stage": "identity-bound leading+oscillatory composite through RF40 lambda turn X_3 with matching independent scoped divergence audit registered",
    "input": "A5 #1002 + A2 #1004 + A4 #1007; A1 #1005 and A3 #1006 retained as non-consumed siblings",
    "new_output": "checksum-bound lambda-turn composite/save-load -> independent public-velocity staged divergence registration",
    "not_output": "A3 correction-side consumption through X_3, power-law composite/audit, cone/global completion, pressure, forcing, complete defect, full residual, or PDE validation",
    "next_shortest_blocker": "recompose the frozen oscillatory runtime onto A1 #1005 power-law leading identity, then carry A3 mean/radial plumbing and A4 independent validation to that same X_4 identity before cone/global completion",
}

_EXPECTED = {
    "parent_a5": PARENT_A5,
    "agent2_lambda_turn_composite": AGENT2_LAMBDA_TURN_COMPOSITE,
    "agent4_lambda_turn_audit": AGENT4_LAMBDA_TURN_AUDIT,
    "agent1_power_law_sibling": AGENT1_POWER_LAW_SIBLING,
    "agent3_first_turn_radial_force_sibling": AGENT3_FIRST_TURN_RADIAL_FORCE_SIBLING,
    "observed_ci_at_freeze": OBSERVED_CI_AT_FREEZE,
    "frozen_science": FROZEN_SCIENCE,
    "final_gate": FINAL_GATE,
    "st006_baseline": ST006_BASELINE,
    "readiness": READINESS,
    "truth_boundary": TRUTH_BOUNDARY,
    "pipeline_position": PIPELINE_POSITION,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _require_git_heads(payload: Mapping[str, Any]) -> None:
    heads = (
        payload["parent_a5"]["head"],
        payload["agent2_lambda_turn_composite"]["head"],
        payload["agent2_lambda_turn_composite"]["leading_head"],
        payload["agent4_lambda_turn_audit"]["head"],
        payload["agent4_lambda_turn_audit"]["audited_agent2_head"],
        payload["agent1_power_law_sibling"]["head"],
        payload["agent3_first_turn_radial_force_sibling"]["head"],
    )
    if any(not isinstance(head, str) or _HEX40.fullmatch(head) is None for head in heads):
        raise ValueError("all registered git heads must be lowercase 40-hex identities")


def registration_payload() -> dict[str, Any]:
    return copy.deepcopy(_EXPECTED)


def materialize_registration_receipt() -> dict[str, Any]:
    registration = registration_payload()
    _require_git_heads(registration)
    if registration["agent4_lambda_turn_audit"]["audited_agent2_head"] != registration[
        "agent2_lambda_turn_composite"
    ]["head"]:
        raise ValueError("A4 lambda-turn audit must bind the exact A2 lambda-turn composite")
    return {
        "schema": SCHEMA_NAME,
        "task_id": TASK_ID,
        "registration": registration,
        "registration_sha256": _sha256(registration),
    }


def enforce_registration_receipt(receipt: Mapping[str, Any]) -> None:
    if receipt.get("schema") != SCHEMA_NAME:
        raise ValueError("schema drifted")
    if receipt.get("task_id") != TASK_ID:
        raise ValueError("task id drifted")
    registration = receipt.get("registration")
    if not isinstance(registration, Mapping):
        raise ValueError("registration must be a mapping")
    digest = receipt.get("registration_sha256")
    if not isinstance(digest, str) or digest != _sha256(registration):
        raise ValueError("registration SHA-256 mismatch")
    _require_git_heads(registration)
    if dict(registration) != _EXPECTED:
        raise ValueError("registration payload drifted from the frozen Agent-5 registration")
    if registration["agent4_lambda_turn_audit"]["audited_agent2_head"] != registration[
        "agent2_lambda_turn_composite"
    ]["head"]:
        raise ValueError("A4/A2 exact-head lineage mismatch")


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    receipt = materialize_registration_receipt()
    enforce_registration_receipt(receipt)
    text = json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
