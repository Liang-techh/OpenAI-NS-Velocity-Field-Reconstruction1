"""Kokuno Agent-5 integration contract for the current RF40 power-law seam.

Integration/provenance glue only. This increment registers the first matched
current lineage through the RF40 constant-lambda power-law endpoint ``X_4``:

    A1 #1005 leading through X_4
      -> A2 #1010 identity-bound leading + frozen complete-curl oscillation
      -> A4 #1012 implementation-distinct public-Cartesian divergence audit.

A3 #1011 is recorded as a lagging correction-side sibling: it materializes
nonlinear m=0 mean attribution only through RF40 axial shutdown ``X_2`` and is
not consumed by the X_4 composite/audit seam. Scoped divergence is not a
complete NS residual, unresolved CI is not scientific admission, and no
correction/PDE state is promoted here.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-current-rf40-power-law-composite-ingest-v1"
TASK_ID = "KOKUNO-A5-CURRENT-RF40-POWER-LAW-COMPOSITE-INGEST-100"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 1008,
    "head": "28860ee1f4655670c62b24880f08aa3d50459623",
    "branch": "codex/kokuno-a5-rf40-lambda-turn-composite-ingest-099",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_current_rf40_lambda_turn_composite_ingest_contract.py",
    "source_blob": "c71499e4e0ae38c9f6eb1032c2e63be7dde58894",
    "test_path": "tests/test_kokuno_a5_current_rf40_lambda_turn_composite_ingest_contract.py",
    "test_blob": "39d54fee64394ab08c90d344098d823e166735f4",
    "workflow_path": ".github/workflows/kokuno-agent5-current-rf40-lambda-turn-composite-ingest.yml",
    "workflow_blob": "5a9e60633cab21096d03a45f76c19bb1f2582303",
}

AGENT1_POWER_LAW_LEADING = {
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
    "fails_closed_after_X_4": True,
}

AGENT2_POWER_LAW_COMPOSITE = {
    "pr": 1010,
    "head": "e36d9da4b7f037e998f5b1658f8c0ea291a76b80",
    "branch": "codex/kokuno-a2-rf40-power-law-composite-087",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_rf40_power_law_leading_oscillatory_identity.py",
    "source_blob": "4f1dc0566c041b9de20f18b4ae9c49d9fd93c58d",
    "test_path": "tests/test_constrained_kokuno_current_rf40_power_law_leading_oscillatory_identity.py",
    "test_blob": "d4858abcce6e5767a7d64c96fdcea18e044eea16",
    "workflow_path": ".github/workflows/kokuno-agent2-current-rf40-power-law-leading-oscillatory.yml",
    "workflow_blob": "f38b4642b886a03eeae229b798202f60ddc3895c",
    "leading_pr": 1005,
    "leading_head": "2c76ebdc41d6c566f43a1305034ba2ff9dce410b",
    "velocity_api": "velocity(x,y,z,t)->[...,3]",
    "materialized_domain": "current leading + frozen complete-curl oscillation through RF40 constant-lambda power-law endpoint X_4",
    "identity_preserving_save_load": True,
    "full_concrete_oscillatory_runtime_digest_bound": True,
    "fails_closed_after_X_4": True,
    "velocity_after_X_4_materialized": False,
    "pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "complete_ns_residual_assessed": False,
}

AGENT4_POWER_LAW_AUDIT = {
    "pr": 1012,
    "head": "9b2be0ffecc14479dbd941eebd663d63f89f2a2f",
    "branch": "codex/kokuno-a4-rf40-power-law-composite-divergence-audit-102",
    "audited_agent2_pr": 1010,
    "audited_agent2_head": "e36d9da4b7f037e998f5b1658f8c0ea291a76b80",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_rf40_power_law_leading_oscillatory_divergence_independent_audit.py",
    "source_blob": "f5f7acf0c02f5f0712943ba841b1bc04478bcdd9",
    "test_path": "tests/test_constrained_kokuno_a4_current_rf40_power_law_leading_oscillatory_divergence_independent_audit.py",
    "test_blob": "3336c8a3524ebc762d8e7a49b06510662b361fae",
    "workflow_path": ".github/workflows/kokuno-agent4-current-rf40-power-law-leading-oscillatory-divergence-audit.yml",
    "workflow_blob": "78d5cfb380664f764387dc0c138d46d51766071c",
    "reference_path": "save/reload A2 #1010, then centered FD2 of public Cartesian velocity only",
    "seed": 9173741,
    "times": [0.31, 0.47, 0.63, 0.71],
    "fd2_step_ladder": [0.02, 0.01, 0.005],
    "strict_power_law_heldout_points": 72,
    "scoped_divergence_gate": 1.0e-5,
    "uses_agent2_production_jacobian_or_divergence": False,
    "uses_agent1_source_derivative_helpers": False,
    "uses_pressure_forcing_or_residual_helper": False,
    "canonical_24_48_96_whole_domain_admission": False,
    "momentum_or_full_ns_evidence": False,
    "scientific_admission": False,
}

AGENT3_AXIAL_SHUTDOWN_MEAN_SIBLING = {
    "pr": 1011,
    "head": "6aabf7b6dd28ae683a3776d518e9a977507911e9",
    "branch": "codex/kokuno-a3-rf40-axial-shutdown-nonlinear-mean-109",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_rf40_axial_shutdown_nonlinear_mean_attribution.py",
    "source_blob": "1ccef44607487d7c3cbdfa6dbce0067fdd302ea5",
    "test_path": "tests/test_constrained_kokuno_current_rf40_axial_shutdown_nonlinear_mean_attribution.py",
    "test_blob": "6f173d8662ce59d235a8d87c92f390626cb3d047",
    "workflow_path": ".github/workflows/kokuno-agent3-current-rf40-axial-shutdown-nonlinear-mean.yml",
    "workflow_blob": "a682c3ba00132d00a33481f1ecddd9ec32a5f824",
    "materialized_domain": "nonlinear m=0 mean attribution through RF40 axial shutdown X_2 only",
    "rf40_axial_shutdown_current_lineage_consumed": True,
    "rf40_lambda_turn_current_lineage_consumed": False,
    "rf40_power_law_current_lineage_consumed": False,
    "radial_inverse_through_X_2_materialized": False,
    "authorized_as_correction_target": False,
    "real_cartesian_correction_velocity": False,
    "consumed_by_this_power_law_seam": False,
}

OBSERVED_CI_AT_FREEZE = {
    "parent_a5_1008": {
        "repository_tests_run": 35608923025,
        "dedicated_run": 35608922952,
        "status": "queued",
        "conclusion": None,
    },
    "agent2_1010_power_law_composite": {
        "repository_tests_run": 35613259366,
        "dedicated_run": 35613258890,
        "status": "queued",
        "conclusion": None,
    },
    "agent3_1011_axial_shutdown_mean_sibling": {
        "repository_tests_run": 35613775429,
        "dedicated_run": 35613775518,
        "status": "queued",
        "conclusion": None,
    },
    "agent4_1012_power_law_audit": {
        "repository_tests_run": 35615457973,
        "dedicated_run": 35615458225,
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
    "current_leading_plus_oscillatory_velocity_through_rf40_lambda_turn_materialized": True,
    "current_leading_plus_oscillatory_velocity_through_rf40_power_law_materialized": True,
    "identity_preserving_power_law_composite_save_load_available": True,
    "agent4_1012_independent_power_law_composite_divergence_audit_registered": True,
    "agent4_1012_independent_power_law_composite_divergence_audit_admitted": False,
    "current_nonlinear_mean_through_axial_shutdown_materialized": True,
    "agent3_1011_registered_as_lagging_nonconsumed_sibling": True,
    "current_nonlinear_radial_stress_through_axial_shutdown_materialized": False,
    "current_nonlinear_radial_force_through_axial_shutdown_materialized": False,
    "current_nonlinear_mean_through_lambda_turn_materialized": False,
    "current_nonlinear_mean_through_power_law_materialized": False,
    "current_nonlinear_radial_stress_through_power_law_materialized": False,
    "current_nonlinear_radial_force_through_power_law_materialized": False,
    "velocity_after_X_4_materialized": False,
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
    "stage": "identity-bound leading+oscillatory composite through RF40 power law X_4 with matching independent scoped divergence audit registered",
    "input": "A5 #1008 + A1 #1005 + A2 #1010 + A4 #1012; A3 #1011 retained as lagging non-consumed correction-side sibling",
    "new_output": "checksum-bound X_4 composite/save-load -> independent public-velocity power-law divergence-audit registration",
    "not_output": "A3 correction-side consumption through X_4, cone/global completion, pressure, forcing, complete defect, full residual, or PDE validation",
    "next_shortest_blocker": "carry A3 nonlinear mean/radial plumbing from X_2 to the exact X_4 composite identity, while A1/A2 continue cone/I1-I4 global completion; only after global velocity completion bind matched pressure and preregistered restricted forcing",
}

_EXPECTED = {
    "parent_a5": PARENT_A5,
    "agent1_power_law_leading": AGENT1_POWER_LAW_LEADING,
    "agent2_power_law_composite": AGENT2_POWER_LAW_COMPOSITE,
    "agent4_power_law_audit": AGENT4_POWER_LAW_AUDIT,
    "agent3_axial_shutdown_mean_sibling": AGENT3_AXIAL_SHUTDOWN_MEAN_SIBLING,
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


def _require_hex40(value: str, label: str) -> None:
    if not _HEX40.fullmatch(value):
        raise ValueError(f"{label} is not an exact 40-hex git identity")


def _enforce_internal_relations(registration: Mapping[str, Any]) -> None:
    for label, obj in (
        ("parent_a5", registration["parent_a5"]),
        ("agent1", registration["agent1_power_law_leading"]),
        ("agent2", registration["agent2_power_law_composite"]),
        ("agent4", registration["agent4_power_law_audit"]),
        ("agent3", registration["agent3_axial_shutdown_mean_sibling"]),
    ):
        _require_hex40(obj["head"], f"{label}.head")
        for blob_key in ("source_blob", "test_blob", "workflow_blob"):
            _require_hex40(obj[blob_key], f"{label}.{blob_key}")

    a1 = registration["agent1_power_law_leading"]
    a2 = registration["agent2_power_law_composite"]
    a4 = registration["agent4_power_law_audit"]
    a3 = registration["agent3_axial_shutdown_mean_sibling"]

    if a2["leading_head"] != a1["head"] or a2["leading_pr"] != a1["pr"]:
        raise ValueError("A2/A1 power-law lineage mismatch")
    if a4["audited_agent2_head"] != a2["head"] or a4["audited_agent2_pr"] != a2["pr"]:
        raise ValueError("A4/A2 power-law audit lineage mismatch")
    if a3["rf40_power_law_current_lineage_consumed"]:
        raise ValueError("A3 lagging sibling cannot be laundered into X_4 consumption")
    if a3["consumed_by_this_power_law_seam"]:
        raise ValueError("A3 lagging sibling cannot be marked consumed by this seam")

    truth = registration["truth_boundary"]
    if not truth["current_leading_plus_oscillatory_velocity_through_rf40_power_law_materialized"]:
        raise ValueError("registered X_4 composite must remain materialized")
    if not truth["identity_preserving_power_law_composite_save_load_available"]:
        raise ValueError("registered X_4 save/load identity must remain available")
    if not truth["agent4_1012_independent_power_law_composite_divergence_audit_registered"]:
        raise ValueError("A4 #1012 audit registration cannot be dropped")
    if truth["agent4_1012_independent_power_law_composite_divergence_audit_admitted"]:
        raise ValueError("unresolved A4 exact-head CI cannot be scientific admission")

    forbidden_true = (
        "current_nonlinear_mean_through_power_law_materialized",
        "current_nonlinear_radial_stress_through_power_law_materialized",
        "current_nonlinear_radial_force_through_power_law_materialized",
        "velocity_after_X_4_materialized",
        "complete_post_xr_rf40_current_lineage_materialized",
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
    for key in forbidden_true:
        if truth[key]:
            raise ValueError(f"truth boundary illegally promoted: {key}")

    readiness = registration["readiness"]
    if readiness != READINESS:
        raise ValueError("readiness state drifted")
    if registration["final_gate"] != FINAL_GATE:
        raise ValueError("final scientific gate drifted")
    if registration["st006_baseline"] != ST006_BASELINE:
        raise ValueError("ST006 baseline drifted")
    if not registration["frozen_science"]["residual_defined_free_forcing_forbidden"]:
        raise ValueError("residual-defined free forcing must remain forbidden")


def materialize_registration_receipt() -> dict[str, Any]:
    registration = copy.deepcopy(_EXPECTED)
    _enforce_internal_relations(registration)
    return {
        "schema": SCHEMA_NAME,
        "task_id": TASK_ID,
        "registration": registration,
        "registration_sha256": _sha256(registration),
    }


def enforce_registration_receipt(receipt: Mapping[str, Any]) -> None:
    if receipt.get("schema") != SCHEMA_NAME:
        raise ValueError("registration schema drifted")
    if receipt.get("task_id") != TASK_ID:
        raise ValueError("registration task id drifted")
    registration = receipt.get("registration")
    if not isinstance(registration, Mapping):
        raise ValueError("registration payload missing")
    digest = receipt.get("registration_sha256")
    if digest != _sha256(registration):
        raise ValueError("registration SHA-256 mismatch")
    if registration != _EXPECTED:
        raise ValueError("registration payload drifted from frozen contract")
    _enforce_internal_relations(registration)


def write_registration_receipt(path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    receipt = materialize_registration_receipt()
    enforce_registration_receipt(receipt)
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return output


def _main() -> int:
    parser = argparse.ArgumentParser(description="Materialize the frozen Kokuno A5 RF40 power-law integration receipt.")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    path = write_registration_receipt(args.output)
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
