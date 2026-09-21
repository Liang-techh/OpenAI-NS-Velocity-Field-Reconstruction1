"""Kokuno Agent-5 integration contract for the current RF40 axial-shutdown seam.

This is integration/provenance glue only.  It registers the first matched current
lineage through the RF40 axial-shutdown endpoint ``X_2``:

    A1 #993 leading velocity through X_2
      -> A2 #999 identity-bound leading + frozen complete-curl oscillation
      -> A4 #1001 implementation-distinct public-Cartesian divergence audit.

The newest A1 #998 lambda-turn leading field and A3 #1000 first-turn radial-stress
work are recorded only as non-consumed siblings.  In particular, this module does
not claim an A3 nonlinear mean/correction through X_2, does not launder scoped
divergence into a complete NS residual, and does not promote unresolved CI to
scientific admission.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-current-rf40-axial-shutdown-composite-ingest-v1"
TASK_ID = "KOKUNO-A5-CURRENT-RF40-AXIAL-SHUTDOWN-COMPOSITE-INGEST-098"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")
_HEX64 = re.compile(r"^[0-9a-f]{64}$")

PARENT_A5 = {
    "pr": 996,
    "head": "88fdafc03ab71e3f8ffe7f39add900e548627daf",
    "branch": "codex/kokuno-a5-rf40-first-turn-composite-ingest-097",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_current_rf40_first_turn_composite_ingest_contract.py",
    "source_blob": "07cd42a7f393bf34b918090630a2a263c2610289",
    "test_path": "tests/test_kokuno_a5_current_rf40_first_turn_composite_ingest_contract.py",
    "test_blob": "89d9874509ad780426de99a673a5a46e79dd4eaf",
    "workflow_path": ".github/workflows/kokuno-agent5-current-rf40-first-turn-composite-ingest.yml",
    "workflow_blob": "22103cad6a1b5d8558cbe9371d93dc9ec8f901b8",
}

AGENT2_AXIAL_SHUTDOWN_COMPOSITE = {
    "pr": 999,
    "head": "93b99292fcf141f24b7c6d7e4fbf95775e6a07e5",
    "branch": "codex/kokuno-a2-rf40-axial-shutdown-composite-085",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_rf40_axial_shutdown_leading_oscillatory_identity.py",
    "source_blob": "255e65ce0c37652858509c33e7c9fad40b73ca97",
    "test_path": "tests/test_constrained_kokuno_current_rf40_axial_shutdown_leading_oscillatory_identity.py",
    "test_blob": "39663c3d7cb13303da08e0479b017e4cb093ad29",
    "workflow_path": ".github/workflows/kokuno-agent2-current-rf40-axial-shutdown-leading-oscillatory.yml",
    "workflow_blob": "f1811d51b247250751d84615150de6680c60739b",
    "leading_pr": 993,
    "leading_head": "2ac6460b483efb1c07f2fa65e7fed781a32f2718",
    "velocity_api": "velocity(x,y,z,t)->[...,3]",
    "materialized_domain": "current leading + frozen complete-curl oscillation through RF40 axial shutdown X_2",
    "identity_preserving_save_load": True,
    "full_concrete_oscillatory_runtime_digest_bound": True,
    "fails_closed_after_X_2": True,
    "pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "complete_ns_residual_assessed": False,
}

AGENT4_AXIAL_SHUTDOWN_AUDIT = {
    "pr": 1001,
    "head": "a4775ffbe89e6106322567a05d2f547173e39b0b",
    "branch": "codex/kokuno-a4-rf40-axial-shutdown-composite-divergence-audit-100",
    "audited_agent2_pr": 999,
    "audited_agent2_head": "93b99292fcf141f24b7c6d7e4fbf95775e6a07e5",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_rf40_axial_shutdown_leading_oscillatory_divergence_independent_audit.py",
    "source_blob": "c72b833e81c7ecff56e0b78d596434e90e13f795",
    "test_path": "tests/test_constrained_kokuno_a4_current_rf40_axial_shutdown_leading_oscillatory_divergence_independent_audit.py",
    "test_blob": "20aa826d6e05de28c948ef2f97262afa1fd1b461",
    "workflow_path": ".github/workflows/kokuno-agent4-current-rf40-axial-shutdown-leading-oscillatory-divergence-audit.yml",
    "workflow_blob": "ec698850df8840a078c7748fd01093bbadaa7e4e",
    "reference_path": "save/reload A2 #999, then centered FD2 of public Cartesian velocity only",
    "seed": 9173721,
    "times": [0.31, 0.47, 0.63, 0.71],
    "fd2_step_ladder": [0.02, 0.01, 0.005],
    "axial_shutdown_rho_zones": [[0.08, 0.28], [0.38, 0.62], [0.72, 0.92]],
    "strict_axial_shutdown_heldout_points": 72,
    "scoped_divergence_gate": 1.0e-5,
    "uses_agent2_production_jacobian_or_divergence": False,
    "uses_agent1_source_derivative_helpers": False,
    "uses_pressure_forcing_or_residual_helper": False,
    "canonical_24_48_96_whole_domain_admission": False,
    "momentum_or_full_ns_evidence": False,
    "scientific_admission": False,
}

AGENT1_LAMBDA_TURN_SIBLING = {
    "pr": 998,
    "head": "43b295444b1e9558222d757cb551e04385eddcb5",
    "branch": "agent-kokuno-1/current-cartesian-rf40-lambda-turn-099",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_rf40_lambda_turn.py",
    "source_blob": "07743c30978360e305a8863e05e7ed322d818b32",
    "test_path": "tests/test_kokuno_pa16_current_cartesian_rf40_lambda_turn.py",
    "test_blob": "efa89ba23365d3a1b72d55d7aaa7837108cbc9a6",
    "workflow_path": ".github/workflows/kokuno-agent1-current-cartesian-rf40-lambda-turn.yml",
    "workflow_blob": "97674e6f567e600300ca0e7169c1cbd02fd6bace",
    "materialized_domain": "leading-only through RF40 lambda turn X_3",
    "consumed_by_agent2_999": False,
    "covered_by_agent4_1001": False,
    "consumed_by_this_increment": False,
}

AGENT3_FIRST_TURN_RADIAL_STRESS_SIBLING = {
    "pr": 1000,
    "head": "1dde810a0a3f2f93c268ca10a86b3443dea14e9f",
    "branch": "codex/kokuno-a3-rf40-first-turn-radial-stress-107",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_rf40_first_turn_nonlinear_radial_stress.py",
    "source_blob": "6334b2fa37d49b60548147e6ae5beca6da652d80",
    "test_path": "tests/test_constrained_kokuno_current_rf40_first_turn_nonlinear_radial_stress.py",
    "test_blob": "dba29f3d13f789707bc4948e00df170cfa50627c",
    "workflow_path": ".github/workflows/kokuno-agent3-current-rf40-first-turn-nonlinear-radial-stress.yml",
    "workflow_blob": "afeb7fb8f703d298af0f7149466477a524c332a6",
    "materialized_domain": "nonlinear radial stress through RF40 first turn X_1 only",
    "rf40_axial_shutdown_current_lineage_consumed": False,
    "authorized_as_correction_target": False,
    "real_cartesian_correction_velocity": False,
    "consumed_by_this_increment": False,
}

OBSERVED_CI_AT_FREEZE = {
    "parent_a5_996": {
        "repository_tests_run": 35596837793,
        "dedicated_run": 35596837676,
        "status": "queued",
        "conclusion": None,
    },
    "agent1_998_lambda_turn_sibling": {
        "repository_tests_run": 35600348853,
        "dedicated_run": 35600348878,
        "status": "queued",
        "conclusion": None,
    },
    "agent2_999_axial_shutdown_composite": {
        "repository_tests_run": 35600785060,
        "dedicated_run": 35600785056,
        "status": "queued",
        "conclusion": None,
    },
    "agent3_1000_first_turn_radial_stress_sibling": {
        "repository_tests_run": 35600886895,
        "dedicated_run": 35600886701,
        "status": "queued",
        "conclusion": None,
    },
    "agent4_1001_axial_shutdown_audit": {
        "repository_tests_run": 35602077548,
        "dedicated_run": 35602077388,
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
    "identity_preserving_axial_shutdown_composite_save_load_available": True,
    "agent4_1001_independent_axial_shutdown_composite_divergence_audit_registered": True,
    "agent4_1001_independent_axial_shutdown_composite_divergence_audit_admitted": False,
    "current_nonlinear_mean_through_axial_shutdown_materialized": False,
    "current_nonlinear_radial_stress_through_axial_shutdown_materialized": False,
    "agent3_1000_first_turn_radial_stress_registered_as_nonconsumed_sibling": True,
    "agent1_998_lambda_turn_leading_registered_as_nonconsumed_sibling": True,
    "leading_plus_oscillatory_velocity_through_lambda_turn_materialized": False,
    "independent_lambda_turn_composite_audit_available": False,
    "complete_post_xr_rf40_current_lineage_materialized": False,
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
    "stage": "identity-bound leading+oscillatory composite through RF40 axial shutdown X_2 with matching independent scoped divergence audit registered",
    "input": "A5 #996 + A2 #999 + A4 #1001; A1 #998 and A3 #1000 retained as non-consumed siblings",
    "new_output": "checksum-bound axial-shutdown composite/save-load -> independent public-velocity staged divergence registration",
    "not_output": "A3 mean/stress/correction through X_2, lambda-turn composite/audit, remaining RF40/global completion, pressure, forcing, complete defect, full residual, or PDE validation",
    "next_shortest_blocker": "recompose frozen oscillatory runtime onto A1 #998 lambda-turn leading identity, then carry A3 mean/radial plumbing and A4 independent validation to that same X_3 identity before remaining global completion",
}

_EXPECTED = {
    "parent_a5": PARENT_A5,
    "agent2_axial_shutdown_composite": AGENT2_AXIAL_SHUTDOWN_COMPOSITE,
    "agent4_axial_shutdown_audit": AGENT4_AXIAL_SHUTDOWN_AUDIT,
    "agent1_lambda_turn_sibling": AGENT1_LAMBDA_TURN_SIBLING,
    "agent3_first_turn_radial_stress_sibling": AGENT3_FIRST_TURN_RADIAL_STRESS_SIBLING,
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
        raise ValueError(f"{label} must be a lowercase 40-hex git identity")


def _require_exact(payload: Mapping[str, Any], key: str, expected: Any) -> None:
    if payload.get(key) != expected:
        raise ValueError(f"{key} drifted from the frozen Agent-5 registration")


def registration_payload() -> dict[str, Any]:
    return copy.deepcopy(_EXPECTED)


def materialize_registration_receipt() -> dict[str, Any]:
    payload = registration_payload()
    for group in (
        payload["parent_a5"],
        payload["agent2_axial_shutdown_composite"],
        payload["agent4_axial_shutdown_audit"],
        payload["agent1_lambda_turn_sibling"],
        payload["agent3_first_turn_radial_stress_sibling"],
    ):
        _require_hex40(group["head"], "head")
        for key, value in group.items():
            if key.endswith("_blob") or key == "source_blob" or key == "test_blob" or key == "workflow_blob":
                _require_hex40(value, key)
    receipt = {
        "schema": SCHEMA_NAME,
        "task": TASK_ID,
        "registration": payload,
    }
    receipt["registration_sha256"] = _sha256(payload)
    return receipt


def enforce_registration_receipt(receipt: Mapping[str, Any]) -> None:
    _require_exact(receipt, "schema", SCHEMA_NAME)
    _require_exact(receipt, "task", TASK_ID)
    registration = receipt.get("registration")
    if registration != _EXPECTED:
        raise ValueError("registration payload drifted from the exact frozen seam")
    digest = receipt.get("registration_sha256")
    if not isinstance(digest, str) or _HEX64.fullmatch(digest) is None:
        raise ValueError("registration SHA-256 is malformed")
    if digest != _sha256(_EXPECTED):
        raise ValueError("registration SHA-256 mismatch")

    truth = registration["truth_boundary"]
    readiness = registration["readiness"]
    gate = registration["final_gate"]
    a2 = registration["agent2_axial_shutdown_composite"]
    a4 = registration["agent4_axial_shutdown_audit"]
    a1 = registration["agent1_lambda_turn_sibling"]
    a3 = registration["agent3_first_turn_radial_stress_sibling"]

    if a4["audited_agent2_head"] != a2["head"]:
        raise ValueError("A4 audit is not bound to the registered A2 composite")
    if a2["leading_head"] != "2ac6460b483efb1c07f2fa65e7fed781a32f2718":
        raise ValueError("A2 axial-shutdown composite detached from A1 #993")
    if a1["consumed_by_this_increment"] or a3["consumed_by_this_increment"]:
        raise ValueError("non-consumed sibling laundering detected")
    if a3["rf40_axial_shutdown_current_lineage_consumed"]:
        raise ValueError("A3 first-turn stress cannot be promoted to axial shutdown")
    if a4["scientific_admission"] or truth["agent4_1001_independent_axial_shutdown_composite_divergence_audit_admitted"]:
        raise ValueError("unresolved A4 audit cannot be scientifically admitted")
    if not truth["current_leading_plus_oscillatory_velocity_through_rf40_axial_shutdown_materialized"]:
        raise ValueError("registered axial-shutdown composite truth was lost")

    hard_false = (
        "current_nonlinear_mean_through_axial_shutdown_materialized",
        "current_nonlinear_radial_stress_through_axial_shutdown_materialized",
        "leading_plus_oscillatory_velocity_through_lambda_turn_materialized",
        "independent_lambda_turn_composite_audit_available",
        "complete_post_xr_rf40_current_lineage_materialized",
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
    for key in hard_false:
        if truth[key] is not False:
            raise ValueError(f"truth boundary laundering detected at {key}")

    if readiness != {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }:
        raise ValueError("readiness was promoted outside the registered scope")
    if gate != {
        "normalized_momentum_sampled_max": 1.0e-3,
        "normalized_momentum_volume_l2": 1.0e-3,
        "divergence_sampled_max": 1.0e-5,
        "divergence_volume_l2": 1.0e-5,
        "canonical_volume_quadrature_ladder": [24, 48, 96],
    }:
        raise ValueError("fixed scientific gate drifted")
    if registration["frozen_science"]["residual_defined_free_forcing_forbidden"] is not True:
        raise ValueError("free residual-defined forcing firewall was disabled")


def write_registration_receipt(path: str | Path) -> dict[str, Any]:
    receipt = materialize_registration_receipt()
    enforce_registration_receipt(receipt)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    write_registration_receipt(args.output)


if __name__ == "__main__":
    main()
