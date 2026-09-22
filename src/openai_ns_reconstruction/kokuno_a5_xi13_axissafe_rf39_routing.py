"""Fail-closed Kokuno Agent-5 routing after the xi=13/A2-axis-safe/A3-RF39 advances.

Integration/provenance glue only.  This module deliberately does not rewrite any
Agent-1--4 mathematics and does not transfer evidence across candidate identities.

Fresh frontiers frozen here:

* Agent 1 #1133 composes the current Cartesian leading field through the public
  pulse endpoint ``xi=13`` using the already-declared repository-autonomous
  principal amplitude and pointwise end-bump realization.  It is callable and
  serializable, but terminal/exterior/global leading completion is still absent.
* Agent 2 #1132 exposes an axis-safe, provider-driven physical Cartesian
  ``velocity(x,y,z,t)`` surface for one localized real conjugate harmonic pair.
  The actual source/background/mode/support provider is still external, so this
  is not a new self-contained project candidate.  The latest self-contained
  leading+oscillatory project candidate remains A2 #1117 / A1 #1107 through
  ``xi=11``.
* Agent 3 #1134 makes the RF34--RF39 compact mean-correction mechanics executable
  on top of #1126.  Its fixture is mechanics-only: the exact current-I4 fixed-Q
  RF30 covariance/derivative backend is still missing, so no current-I4 RF30
  defect or applied correction is admitted.
* Agent 4 #1127 remains the matching implementation-distinct audit authority for
  the older exact current-I4 radial-force path A3 #1118.  Its exact-head Actions
  were still queued at this checkpoint; it is scoped operator-consistency
  evidence only and cannot authorize a complete NS correction.

The final project gates remain fixed: normalized momentum sampled-max and
volume-L2 <= 1e-3; normalized divergence sampled-max and volume-L2 <= 1e-5;
canonical quadrature [24,48,96]; residual-defined/free forcing forbidden.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a5_pulse_mj_physical_harmonic_rf31_routing as parent_a5

SCHEMA = "kokuno-a5-xi13-axissafe-rf39-routing-v1"
TASK = "KOKUNO-A5-XI13-AXISSAFE-RF39-ROUTING-115"

PARENT_A5 = {
    "pr": 1128,
    "head": "f69f5842ea77d0837fc47ed9036b2e7cf799d101",
    "branch": "kokuno-agent5/pulse-mj-physical-harmonic-rf31-routing-114",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_pulse_mj_physical_harmonic_rf31_routing.py",
    "source_blob": "4a0d14b06301b8d2db6f5925676f29eea38e3c52",
    "task": "KOKUNO-A5-PULSE-MJ-PHYSICAL-HARMONIC-RF31-ROUTING-114",
    "schema": "kokuno-a5-pulse-mj-physical-harmonic-rf31-routing-v1",
}

AGENT1_XI13_LEADING = {
    "pr": 1133,
    "head": "3b4af71c4ab547bc11f9f3e320afa4b62c25b930",
    "branch": "agent1/current-cartesian-pulse-end-20260922",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_pulse_end_compensated.py",
    "source_blob": "4b8fad0781d74c934b6cc6fe7805a11d16807a92",
    "parent_pr": 1124,
    "stage": "current-cartesian-leading-through-public-pulse-end-xi13",
    "cartesian_end_compensation_composed": True,
    "public_pulse_endpoint_xi": 13.0,
    "save_load_materialized": True,
    "source_exact_amplitude_materialized": False,
    "source_exact_hidden_bump_recovered": False,
    "terminal_global_leading_velocity_materialized": False,
    "repository_tests_run": 35691163240,
    "dedicated_run": 35691163289,
    "actions_status_at_registration": "queued",
}

AGENT2_AXIS_SAFE_PHYSICAL_HARMONIC = {
    "pr": 1132,
    "head": "1693c6b0e0e337e7fab3914bbcfcab8ad54bbf3a",
    "branch": "codex/k2-osc-103-axis-safe-physical-velocity",
    "source_path": "src/openai_ns_reconstruction/kokuno_source_axis_safe_physical_velocity.py",
    "source_blob": "9067cb68dc49e08a3da7aaa429366f11d2073cdd",
    "parent_pr": 1125,
    "stage": "provider-driven-axis-safe-source-physical-cartesian-harmonic",
    "provider_driven_velocity_xyzt_materialized": True,
    "axis_zero_core_fail_closed": True,
    "project_domain_source_input_provider_materialized": False,
    "global_axis_safe_source_velocity_materialized": False,
    "self_contained_velocity_xyzt_provider": False,
    "repository_tests_run": 35691160766,
    "dedicated_run": 35691160815,
    "actions_status_at_registration": "queued",
}

LATEST_SELF_CONTAINED_PROJECT_COMPOSITE = {
    "pr": 1117,
    "head": "27741d9c0a27262f7fabf61eebaa2fbd507e9f03",
    "consumed_agent1_pr": 1107,
    "consumed_agent1_head": "45da043dd2b4cd067f005a72c7e21fd0f2bcf309",
    "stage": "project-domain-logX-leading-plus-frozen-complete-curl-through-xi11",
    "public_pulse_endpoint_xi": 11.0,
    "save_load_materialized": True,
    "matching_xi13_agent2_composite": False,
    "matching_agent4_audit_present": False,
}

AGENT3_RF34_RF39_MECHANICS = {
    "pr": 1134,
    "head": "ecbc7e3702e366a5f355800dd662b11766273a41",
    "branch": "agent3/rf34-rf39-compact-mean-correction-20260922",
    "source_path": "src/openai_ns_reconstruction/kokuno_rf34_rf39_compact_mean_correction.py",
    "source_blob": "392975eafd3787c91e2dc82ad9b735567513271d",
    "parent_pr": 1126,
    "stage": "RF34-RF39-source-specific-compact-mean-correction-mechanics",
    "rf34_rf39_compact_realization_executable": True,
    "autonomous_pointwise_bump_realization": True,
    "source_exact_pointwise_bump_recovered": False,
    "current_i4_source_chart_backend_materialized": False,
    "current_i4_rf30_defect_materialized": False,
    "current_i4_rf34_rf39_correction_materialized": False,
    "candidate_residual_evidence": False,
    "cartesian_correction_velocity_materialized": False,
    "repository_tests_run": 35691928795,
    "dedicated_run": 35691928711,
    "actions_status_at_registration": "queued",
}

CURRENT_I4_REAL_RADIAL_FORCE = {
    "pr": 1118,
    "head": "922f7aa10460ded44af212d313eceff33a2ac647",
    "consumed_agent2_pr": 1080,
    "consumed_agent1_pr": 1079,
    "stage": "current-I4-nonlinear-radial-force-partial_z_sigma1",
    "radial_force_materialized": True,
    "complete_ns_correction_authorized": False,
}

AGENT4_CURRENT_I4_RADIAL_FORCE_AUDIT = {
    "task": "K4-VAL-116",
    "pr": 1127,
    "head": "add922982c2f3a4cf5a72d84910817e63f05c487",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_i4_radial_force_independent_audit.py",
    "source_blob": "c254b248e2703eaec1c4f757657ca40b9e211c7d",
    "audited_agent3_pr": 1118,
    "audited_agent3_head": "922f7aa10460ded44af212d313eceff33a2ac647",
    "stage": "current-I4-radial-force-independent-FD4-audit",
    "present": True,
    "registered": True,
    "scoped_gate_passed": None,
    "scientifically_admitted": False,
    "authorizes_complete_ns_correction": False,
    "repository_tests_run": 35688548029,
    "dedicated_run": 35688548012,
    "actions_status_at_registration": "queued",
}

CORRECTED_SOURCE = {
    "repository": "KokunoYumeto/yang-mills-interacting-workbench",
    "commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
    "path": "navier-stokes/navier_stokes_workbench.tex",
    "blob": "205a99807302e21a51c5eaf223390c0dfc42bcd0",
    "release_date": "2026-09-09",
    "role": "structural provenance only; autonomous realization remains explicitly labeled",
}

FROZEN_GATES = {
    "normalized_momentum_sampled_max": 1.0e-3,
    "normalized_momentum_volume_l2": 1.0e-3,
    "normalized_divergence_sampled_max": 1.0e-5,
    "normalized_divergence_volume_l2": 1.0e-5,
    "canonical_quadrature": [24, 48, 96],
    "residual_defined_free_forcing_allowed": False,
    "st006_momentum_sampled_max": 0.1082289305112118,
    "st006_momentum_volume_l2": 0.10758432876230622,
}

CORE_STATE = {
    "leading_ready": False,
    "oscillatory_ready": True,
    "correction_ready": False,
    "velocity_export_ready": False,
    "pde_validated": False,
}


def _truth_boundary() -> dict[str, Any]:
    return {
        "current_cartesian_leading_through_xi13_materialized": True,
        "current_cartesian_end_compensation_composed": True,
        "source_exact_main_pulse_amplitude_materialized": False,
        "source_exact_hidden_end_bump_recovered": False,
        "terminal_global_leading_velocity_materialized": False,
        "matching_xi13_agent2_project_composite_materialized": False,
        "matching_xi13_agent4_audit_present": False,
        "provider_driven_axis_safe_source_harmonic_velocity_materialized": True,
        "project_domain_source_input_provider_materialized": False,
        "global_axis_safe_source_velocity_materialized": False,
        "source_harmonic_self_contained_velocity_xyzt_provider": False,
        "latest_self_contained_project_composite_is_agent2_1117": True,
        "rf34_rf39_compact_mechanics_materialized": True,
        "current_i4_source_chart_backend_materialized": False,
        "current_i4_rf30_defect_materialized": False,
        "current_i4_rf34_rf39_correction_materialized": False,
        "cartesian_correction_velocity_materialized": False,
        "current_i4_radial_force_materialized": True,
        "agent4_matching_current_i4_radial_force_audit_present": True,
        "agent4_matching_current_i4_radial_force_audit_registered": True,
        "agent4_matching_current_i4_radial_force_audit_scoped_gate_passed": None,
        "agent4_matching_current_i4_radial_force_audit_admitted": False,
        "scoped_current_i4_force_authorized_as_complete_ns_correction_target": False,
        "matched_cartesian_pressure_gradient_materialized": False,
        "preregistered_restricted_forcing_materialized": False,
        "restricted_forcing_proved_not_residual_defined": False,
        "complete_candidate_api_ready": False,
        "complete_identity_bound_ns_defect_materialized": False,
        "finite_correction_cycle_run": False,
        "heldout_complete_ns_residual_assessed": False,
        "canonical_whole_domain_admission_run": False,
        "same_protocol_comparable_to_st006": False,
        "scientific_admission": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
        "pde_validated": False,
    }


def build_identity_firewall() -> dict[str, bool]:
    return {
        "agent1_1133_xi13_not_consumed_by_agent2_1117_xi11": True,
        "agent2_1132_provider_adapter_not_self_contained_project_composite": True,
        "agent3_1134_mechanics_not_current_i4_candidate_residual_evidence": True,
        "agent4_1127_current_i4_force_audit_not_xi13_or_rf39_evidence": True,
        "cross_identity_evidence_transfer_allowed": False,
    }


def _canonical_bytes(value: Mapping[str, Any]) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _without_digest(registration: Mapping[str, Any]) -> dict[str, Any]:
    payload = copy.deepcopy(dict(registration))
    payload.pop("registration_sha256", None)
    return payload


def registration_sha256(registration: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_bytes(_without_digest(registration))).hexdigest()


def build_registration() -> dict[str, Any]:
    registration: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "role": "Agent-5 integration/provenance/artifact/CI routing only",
        "parent_a5": copy.deepcopy(PARENT_A5),
        "corrected_source": copy.deepcopy(CORRECTED_SOURCE),
        "frontiers": {
            "leading": copy.deepcopy(AGENT1_XI13_LEADING),
            "source_axis_safe_harmonic_sibling": copy.deepcopy(AGENT2_AXIS_SAFE_PHYSICAL_HARMONIC),
            "latest_self_contained_project_composite": copy.deepcopy(LATEST_SELF_CONTAINED_PROJECT_COMPOSITE),
            "rf34_rf39_mechanics_sibling": copy.deepcopy(AGENT3_RF34_RF39_MECHANICS),
            "current_i4_real_radial_force": copy.deepcopy(CURRENT_I4_REAL_RADIAL_FORCE),
            "current_i4_radial_force_validator": copy.deepcopy(AGENT4_CURRENT_I4_RADIAL_FORCE_AUDIT),
        },
        "identity_firewall": build_identity_firewall(),
        "truth_boundary": _truth_boundary(),
        "core_state": copy.deepcopy(CORE_STATE),
        "frozen_gates": copy.deepcopy(FROZEN_GATES),
        "shortest_closure": [
            "recompose A2 on exact A1 #1133 xi13 identity and obtain matching A4 scoped audit",
            "materialize the real fixed-Q RF30 covariance/derivative backend before A3 #1134 can become candidate-specific correction evidence",
            "finish terminal/global leading and bind matched pressure plus preregistered restricted forcing on one identity",
            "only then form complete NS defect, run real finite correction cycle, and let A4 run the held-out full-NS 1e-3 gate",
        ],
    }
    registration["registration_sha256"] = registration_sha256(registration)
    validate_registration(registration)
    return registration


def validate_registration(registration: Mapping[str, Any]) -> None:
    reg = dict(registration)
    expected_digest = registration_sha256(reg)
    if reg.get("registration_sha256") != expected_digest:
        raise ValueError("registration digest mismatch")
    if reg.get("schema") != SCHEMA or reg.get("task") != TASK:
        raise ValueError("schema/task identity drift")
    if reg.get("parent_a5") != PARENT_A5:
        raise ValueError("parent A5 identity drift")
    if reg.get("corrected_source") != CORRECTED_SOURCE:
        raise ValueError("corrected-source provenance drift")

    frontiers = reg.get("frontiers")
    if not isinstance(frontiers, Mapping):
        raise ValueError("frontiers missing")
    expected_frontiers = {
        "leading": AGENT1_XI13_LEADING,
        "source_axis_safe_harmonic_sibling": AGENT2_AXIS_SAFE_PHYSICAL_HARMONIC,
        "latest_self_contained_project_composite": LATEST_SELF_CONTAINED_PROJECT_COMPOSITE,
        "rf34_rf39_mechanics_sibling": AGENT3_RF34_RF39_MECHANICS,
        "current_i4_real_radial_force": CURRENT_I4_REAL_RADIAL_FORCE,
        "current_i4_radial_force_validator": AGENT4_CURRENT_I4_RADIAL_FORCE_AUDIT,
    }
    if dict(frontiers) != expected_frontiers:
        raise ValueError("frontier identity drift")

    if reg.get("identity_firewall") != build_identity_firewall():
        raise ValueError("identity firewall drift")
    if reg.get("core_state") != CORE_STATE:
        raise ValueError("core readiness promotion/drift")
    if reg.get("frozen_gates") != FROZEN_GATES:
        raise ValueError("frozen project gates changed")
    if reg.get("truth_boundary") != _truth_boundary():
        raise ValueError("truth boundary promotion/drift")


def write_registration(path: str | Path) -> dict[str, Any]:
    registration = build_registration()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(registration, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return registration


def load_registration(path: str | Path) -> dict[str, Any]:
    registration = json.loads(Path(path).read_text(encoding="utf-8"))
    validate_registration(registration)
    return registration


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    registration = write_registration(args.output)
    print(json.dumps({"registration_sha256": registration["registration_sha256"], "output": str(args.output)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
