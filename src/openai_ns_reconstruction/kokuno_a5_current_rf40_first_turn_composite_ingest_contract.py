"""Kokuno Agent-5 integration contract for the current RF40 first-turn seam.

This module is provenance/integration glue only.  It checksum-binds the first
matched current-lineage chain that exists strictly beyond ``X_R``:

    A1 #986 leading velocity through X_1=e X_R
      -> A2 #992 identity-bound leading + complete-curl oscillatory velocity
      -> A3 #994 nonlinear cylindrical m=0 mean attribution
      -> A4 #995 implementation-distinct public-velocity divergence audit

A1 #993 is one leading-only stage farther, through RF40 axial shutdown ``X_2``.
It is deliberately registered only as a non-consumed sibling because no matching
A2/A3/A4 axial-shutdown composite exists at this freeze.

Registration is not scientific admission.  In particular, A4 #995 supplies only
scoped divergence evidence, not momentum/full-NS residual evidence and not the
canonical [24,48,96] whole-domain admission.  The fixed 1e-3 momentum gate,
1e-5 divergence gate, and ban on residual-defined free forcing remain unchanged.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-current-rf40-first-turn-composite-ingest-v1"
TASK_ID = "KOKUNO-A5-CURRENT-RF40-FIRST-TURN-COMPOSITE-INGEST-097"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")
_HEX64 = re.compile(r"^[0-9a-f]{64}$")

PARENT_A5 = {
    "pr": 990,
    "head": "92def17ccab4192a7bd270a4f4d637138e953266",
    "branch": "codex/kokuno-a5-current-exterior-xr-composite-ingest-096",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_current_exterior_xr_composite_ingest_contract.py",
    "source_blob": "c9e09bbbf41e6524bb0c3db62a7e2b8c417fa09e",
    "role": "register matched leading+oscillatory/nonlinear-mean/divergence chain through X_R",
}

AGENT2_FIRST_TURN_COMPOSITE = {
    "pr": 992,
    "head": "2a2fad307a674ec1eb5a2bc426c0d01c3a0e7377",
    "branch": "codex/kokuno-a2-rf40-first-turn-composite-084",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_rf40_first_turn_leading_oscillatory_identity.py",
    "source_blob": "c126dace161e5dae1745185c79152736817994b8",
    "test_path": "tests/test_constrained_kokuno_current_rf40_first_turn_leading_oscillatory_identity.py",
    "test_blob": "173ddce0ce9e1b5e89a663fd450e0fc135e156b3",
    "workflow_path": ".github/workflows/kokuno-agent2-current-rf40-first-turn-leading-oscillatory.yml",
    "workflow_blob": "6931746518d2e54cf178550f21f0ef0d5b531ead",
    "leading_pr": 986,
    "leading_head": "23c00b98526e187ff04d432745300a084ec859f2",
    "velocity_api": "velocity(x,y,z,t)->[...,3]",
    "materialized_domain": "current leading + frozen complete-curl oscillation through X_1=e X_R",
    "identity_preserving_save_load": True,
    "complete_oscillatory_runtime_digest_bound": True,
    "fails_closed_after_X_1": True,
    "axial_shutdown_composite_materialized": False,
    "matched_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "complete_ns_residual_assessed": False,
}

AGENT3_FIRST_TURN_NONLINEAR_MEAN = {
    "pr": 994,
    "head": "6088e3d3055e2792ba20fbd703fcc238ce1b27b5",
    "branch": "codex/kokuno-a3-rf40-first-turn-nonlinear-mean-106",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_rf40_first_turn_nonlinear_mean_attribution.py",
    "source_blob": "0f0e1f5c31cf33e28fb6f2851c2df7250bc582d1",
    "test_path": "tests/test_constrained_kokuno_current_rf40_first_turn_nonlinear_mean_attribution.py",
    "test_blob": "580692ba7045f9c22419b756b4a736253b852edf",
    "workflow_path": ".github/workflows/kokuno-agent3-current-rf40-first-turn-nonlinear-mean.yml",
    "workflow_blob": "88e565296814e284ff03c36578049ecfe71f66c8",
    "consumed_agent2_pr": 992,
    "consumed_agent2_head": "2a2fad307a674ec1eb5a2bc426c0d01c3a0e7377",
    "decomposition": "A=(u_lead·grad)u_osc; B=(u_osc·grad)u_lead; Q=(u_osc·grad)u_osc; N=A+B+Q",
    "projector": "existing rotating cylindrical m=0 projector",
    "materialized_domain": "through RF40 first turn X_1=e X_R",
    "radial_inverse_in_this_increment": False,
    "complete_ns_defect": False,
    "authorized_as_correction_target": False,
    "real_cartesian_correction_velocity": False,
}

AGENT4_FIRST_TURN_AUDIT = {
    "pr": 995,
    "head": "70cd44ed0b5cabbcf5103df8959ff008377d17fa",
    "branch": "codex/kokuno-a4-rf40-first-turn-composite-divergence-audit-099",
    "audited_agent2_pr": 992,
    "audited_agent2_head": "2a2fad307a674ec1eb5a2bc426c0d01c3a0e7377",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_rf40_first_turn_leading_oscillatory_divergence_independent_audit.py",
    "source_blob": "3881a0bb98e6369cf8fde85556c9d7ac0c2db7fc",
    "test_path": "tests/test_constrained_kokuno_a4_current_rf40_first_turn_leading_oscillatory_divergence_independent_audit.py",
    "test_blob": "edc7bb9add03adb3babddd9961d1c52e614fd679",
    "workflow_path": ".github/workflows/kokuno-agent4-current-rf40-first-turn-leading-oscillatory-divergence-audit.yml",
    "workflow_blob": "58fc18ff4fb1410603a7d45a5a07f6e6b83ade8a",
    "reference_path": "save/reloaded public Cartesian velocity only; centered FD2 Jacobians",
    "uses_agent2_production_jacobian_or_divergence": False,
    "uses_agent1_source_derivative_helpers": False,
    "uses_pressure_forcing_or_residual_helper": False,
    "seed": 9173711,
    "times": [0.31, 0.47, 0.63, 0.71],
    "fd2_step_ladder": [0.02, 0.01, 0.005],
    "post_xr_log_strata": [[0.08, 0.28], [0.38, 0.62], [0.72, 0.92]],
    "post_xr_heldout_points": 72,
    "scoped_divergence_gate": 1.0e-5,
    "nontrivial_total_speed_rms_floor": 1.0e-10,
    "nontrivial_oscillatory_speed_rms_floor": 1.0e-12,
    "medium_to_fine_degradation_ratio_cap": 1.25,
    "numerical_floor": 2.0e-8,
    "scoped_divergence_evidence_only": True,
    "canonical_whole_domain_admission": False,
    "momentum_or_full_ns_evidence": False,
    "scientific_admission": False,
}

AGENT1_AXIAL_SHUTDOWN_SIBLING = {
    "pr": 993,
    "head": "2ac6460b483efb1c07f2fa65e7fed781a32f2718",
    "branch": "agent-kokuno-1/current-cartesian-rf40-axial-shutdown-098",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_rf40_axial_shutdown.py",
    "source_blob": "057a0514c8a941c3e59922b8158f480434b4441e",
    "test_path": "tests/test_kokuno_pa16_current_cartesian_rf40_axial_shutdown.py",
    "test_blob": "b1703748f567d7d7d505bf7f2e07aa287a8c0559",
    "workflow_path": ".github/workflows/kokuno-agent1-current-cartesian-rf40-axial-shutdown.yml",
    "workflow_blob": "a5f485f89456c7171bf3789d11f677c4d800658a",
    "materialized_domain": "leading-only through RF40 axial shutdown X_2",
    "preserves_current_incompressibility_memory": True,
    "lambda_turn_materialized": False,
    "power_law_continuation_materialized": False,
    "cone_i1_i2_i3_i4_completed": False,
    "consumed_by_agent2_992": False,
    "consumed_by_agent3_994": False,
    "covered_by_agent4_995": False,
    "consumed_by_this_increment": False,
}

OBSERVED_CI_AT_FREEZE = {
    "parent_a5_990": {
        "repository_tests_run": 35591673799,
        "dedicated_run": 35591673815,
        "status": "queued",
        "conclusion": None,
    },
    "agent2_992": {
        "repository_tests_run": 35594694349,
        "dedicated_run": 35594694416,
        "status": "queued",
        "conclusion": None,
    },
    "agent1_993_sibling": {
        "repository_tests_run": 35594713221,
        "dedicated_run": 35594713191,
        "status": "queued",
        "conclusion": None,
    },
    "agent3_994": {
        "repository_tests_run": 35595384695,
        "dedicated_run": 35595384759,
        "status": "queued",
        "conclusion": None,
    },
    "agent4_995": {
        "repository_tests_run": 35596166026,
        "dedicated_run": 35596166074,
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
    "identity_preserving_first_turn_composite_save_load_available": True,
    "current_nonlinear_m0_mean_through_rf40_first_turn_materialized": True,
    "agent4_995_independent_first_turn_composite_divergence_audit_registered": True,
    "agent4_995_independent_first_turn_composite_divergence_audit_admitted": False,
    "agent1_993_axial_shutdown_leading_only_materialized_as_sibling": True,
    "current_leading_plus_oscillatory_velocity_through_axial_shutdown_materialized": False,
    "current_nonlinear_mean_through_axial_shutdown_materialized": False,
    "independent_axial_shutdown_composite_audit_available": False,
    "complete_post_xr_rf40_current_lineage_materialized": False,
    "rf40_lambda_turn_current_lineage_materialized": False,
    "rf40_power_law_current_lineage_materialized": False,
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
    "stage": "identity-bound current leading+oscillatory composite and nonlinear mean through RF40 first turn, with matching independent staged divergence audit registered",
    "input": "A5 #990 + A2 #992 + A3 #994 + A4 #995; A1 #993 retained as newer non-consumed axial-shutdown leading-only sibling",
    "new_output": "checksum-bound first-turn composite/save-load -> nonlinear m=0 mean -> independent public-velocity staged divergence registration",
    "not_output": "axial-shutdown composite/mean/audit, remaining RF40 and cone/I1-I4 global completion, pressure, forcing, complete defect, correction velocity, full residual, or PDE validation",
    "next_shortest_blocker": "recompose exact oscillatory runtime onto A1 #993 axial-shutdown leading identity, then carry A3 mean/radial plumbing and A4 independent validation forward before matched pressure/restricted forcing and complete-defect correction",
}

_EXPECTED_TOP_LEVEL = {
    "parent_a5": PARENT_A5,
    "agent2_first_turn_composite": AGENT2_FIRST_TURN_COMPOSITE,
    "agent3_first_turn_nonlinear_mean": AGENT3_FIRST_TURN_NONLINEAR_MEAN,
    "agent4_first_turn_audit": AGENT4_FIRST_TURN_AUDIT,
    "agent1_axial_shutdown_sibling": AGENT1_AXIAL_SHUTDOWN_SIBLING,
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


def _require_hex40(value: Any, label: str) -> None:
    if not isinstance(value, str) or _HEX40.fullmatch(value) is None:
        raise ValueError(f"{label} must be one lowercase 40-hex SHA")


def _require_hex64(value: Any, label: str) -> None:
    if not isinstance(value, str) or _HEX64.fullmatch(value) is None:
        raise ValueError(f"{label} must be one lowercase 64-hex SHA-256")


def build_contract(exact_head: str) -> dict[str, Any]:
    """Build the immutable registration receipt for this bounded A5 increment."""
    _require_hex40(exact_head, "exact_head")
    payload: dict[str, Any] = {
        "schema": SCHEMA_NAME,
        "task_id": TASK_ID,
        "agent5_exact_head": exact_head,
    }
    for key, value in _EXPECTED_TOP_LEVEL.items():
        payload[key] = copy.deepcopy(value)
    payload["registration_sha256"] = _sha256(payload)
    return payload


def verify_contract(contract: Mapping[str, Any]) -> None:
    """Fail closed if provenance, gates, readiness, or truth boundaries drift."""
    if contract.get("schema") != SCHEMA_NAME or contract.get("task_id") != TASK_ID:
        raise ValueError("Agent-5 registration schema/task drifted")
    _require_hex40(contract.get("agent5_exact_head"), "agent5_exact_head")

    for key, expected in _EXPECTED_TOP_LEVEL.items():
        if contract.get(key) != expected:
            raise ValueError(f"Agent-5 registration drifted at {key}")

    recorded_digest = contract.get("registration_sha256")
    _require_hex64(recorded_digest, "registration_sha256")
    unsigned = dict(contract)
    unsigned.pop("registration_sha256", None)
    if _sha256(unsigned) != recorded_digest:
        raise ValueError("Agent-5 registration digest mismatch")

    truth = contract["truth_boundary"]
    forbidden_true = (
        "agent4_995_independent_first_turn_composite_divergence_audit_admitted",
        "current_leading_plus_oscillatory_velocity_through_axial_shutdown_materialized",
        "current_nonlinear_mean_through_axial_shutdown_materialized",
        "independent_axial_shutdown_composite_audit_available",
        "complete_post_xr_rf40_current_lineage_materialized",
        "rf40_lambda_turn_current_lineage_materialized",
        "rf40_power_law_current_lineage_materialized",
        "cone_i1_i2_i3_i4_outer_overlays_completed",
        "outer_global_leading_velocity_materialized",
        "global_compact_support_completed",
        "matched_cartesian_pressure_materialized",
        "cartesian_matched_pressure_gradient_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect_materialized",
        "current_nonlinear_mean_authorized_as_correction_target",
        "real_agent3_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "complete_candidate_api_ready",
        "canonical_whole_domain_divergence_l2_assessed",
        "heldout_normalized_ns_residual_assessed",
        "same_protocol_st006_comparison_available_now",
        "residual_reduction_claimed",
        "scientific_admission",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    )
    promoted = [key for key in forbidden_true if truth.get(key) is True]
    if promoted:
        raise ValueError(f"forbidden truth promotion: {promoted}")

    if contract["agent4_first_turn_audit"]["scientific_admission"] is not False:
        raise ValueError("queued A4 evidence cannot be scientifically admitted")
    sibling = contract["agent1_axial_shutdown_sibling"]
    if any(
        sibling[key] is not False
        for key in (
            "consumed_by_agent2_992",
            "consumed_by_agent3_994",
            "covered_by_agent4_995",
            "consumed_by_this_increment",
        )
    ):
        raise ValueError("A1 #993 axial-shutdown sibling was silently laundered into matched lineage")


def write_contract(path: str | Path, exact_head: str) -> dict[str, Any]:
    contract = build_contract(exact_head)
    verify_contract(contract)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return contract


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-head", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    contract = write_contract(args.output, args.exact_head)
    print(contract["registration_sha256"])


if __name__ == "__main__":
    main()
