"""Fail-closed Kokuno Agent-5 routing checkpoint after A1/A2/A3 advances.

This module is integration/provenance glue only.  It records the exact scientific
frontiers without transferring evidence between semantic candidate identities:

* A1 #1124 materializes current-lineage pulse-entry M/J bookkeeping and can feed
  the already-public #1116 two-row c1/c2 end solve.  It does not compose that
  end correction into the Cartesian leading field and does not recover the
  source-exact amplitude or hidden bump realization.
* A2 #1125 physicalizes the corrected-source localized complete-curl harmonic.
  It remains caller-dependent and is not a self-contained project-domain
  ``velocity(x,y,z,t)`` provider, so the latest self-contained project composite
  remains exact A2 #1117 / A1 #1107 through xi=11.
* A3 #1126 materializes a typed RF30/RF31 formula executor and five-row mechanics
  contract.  The exact current-I4 source fixed-Q covariance backend is absent,
  so #1126 is not current-I4 candidate defect evidence and solves/applies no
  correction.  The actual correction-side current-I4 frontier remains A3 #1118
  radial force ``partial_z sigma_1`` on the older A2 #1080 / A1 #1079 identity.
* A4 #1119 remains the delivered implementation-distinct current-I4 stress audit.
  K4-VAL-116 has been preregistered for the #1118 radial force, but at this
  checkpoint it is a claim only: there is no delivered exact head to register or
  admit.  A pending claim is never treated as scientific evidence.

No pressure, forcing, residual, optimization sample, correction coefficient or
scientific threshold is accepted by this module.  The final held-out full-NS
admission gate remains normalized momentum max/L2 <= 1e-3 and divergence
max/L2 <= 1e-5, with canonical quadrature [24,48,96] and residual-defined/free
forcing forbidden.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a5_xi11_composite_i4_force_routing as parent_a5

SCHEMA = "kokuno-a5-pulse-mj-physical-harmonic-rf31-routing-v1"
TASK = "KOKUNO-A5-PULSE-MJ-PHYSICAL-HARMONIC-RF31-ROUTING-114"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 1120,
    "head": "337ad149f6227fd054e496f4f38e48216b2ff2f4",
    "branch": "kokuno-agent5/xi11-composite-i4-force-routing-113",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_xi11_composite_i4_force_routing.py",
    "source_blob": "9559ea30dcc6bea7c4d8179b555e4884108c240a",
    "task": "KOKUNO-A5-XI11-COMPOSITE-I4-FORCE-ROUTING-113",
    "schema": "kokuno-a5-xi11-composite-i4-force-routing-v1",
}

AGENT1_CURRENT_PULSE_ENTRY_MJ = {
    "pr": 1124,
    "head": "4be97cbeb3b78b4d6a162447ae8fea73d0c2c477",
    "branch": "agent1/current-pulse-entry-mj-20260922",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_pulse_entry_moments.py",
    "source_blob": "d4cad0001552a88f682230a6b7fc13c993797497",
    "parent_pr": 1116,
    "parent_head": "c144d00fd81a938e497ca5d841fed6d1fe448a8b",
    "stage": "current-lineage-pulse-entry-MJ-bookkeeping",
    "current_lineage_J_entry_materialized": True,
    "current_J_assumed_zero": False,
    "source_exact_amplitude_root_materialized": False,
    "cartesian_end_compensation_composed": False,
    "terminal_global_leading_velocity_materialized": False,
    "repository_tests_run": 35687207886,
    "dedicated_run": 35687207951,
    "actions_status_at_registration": "queued",
}

AGENT2_SOURCE_PHYSICALIZED_HARMONIC = {
    "pr": 1125,
    "head": "bab607c1a4f93a586dc8a1a9fa3c0d8b5badeb58",
    "branch": "codex/k2-osc-102-physical-scaling",
    "source_path": "src/openai_ns_reconstruction/kokuno_source_physicalized_localized_harmonic.py",
    "source_blob": "b19805fc2e78d6a51b6ead6c141125a662d3ec4c",
    "parent_pr": 1117,
    "parent_head": "27741d9c0a27262f7fabf61eebaa2fbd507e9f03",
    "stage": "corrected-source-localized-complete-curl-physical-scaling",
    "source_physical_Q_scaling_applied": True,
    "self_contained_velocity_xyzt_provider": False,
    "deterministic_source_input_provider_materialized": False,
    "global_axis_safe_source_velocity_materialized": False,
    "repository_tests_run": 35687379497,
    "dedicated_run": 35687379546,
    "actions_status_at_registration": "queued",
}

LATEST_SELF_CONTAINED_PROJECT_COMPOSITE = {
    "pr": 1117,
    "head": "27741d9c0a27262f7fabf61eebaa2fbd507e9f03",
    "consumed_agent1_pr": 1107,
    "consumed_agent1_head": "45da043dd2b4cd067f005a72c7e21fd0f2bcf309",
    "stage": "project-domain-logX-leading-plus-frozen-complete-curl-through-xi11",
    "save_load_materialized": True,
    "terminal_global_velocity_materialized": False,
    "matching_agent4_audit_present": False,
}

CURRENT_I4_RADIAL_FORCE = {
    "pr": 1118,
    "head": "922f7aa10460ded44af212d313eceff33a2ac647",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_i4_nonlinear_radial_force.py",
    "source_blob": "e28eac6f7fb7f3191a024b052bcf677cff8d74d2",
    "consumed_agent2_pr": 1080,
    "consumed_agent1_pr": 1079,
    "stage": "current-I4-nonlinear-radial-force-partial_z_sigma1",
    "radial_force_materialized": True,
    "complete_ns_correction_authorized": False,
}

AGENT3_RF30_RF31_TYPED_EXECUTOR = {
    "pr": 1126,
    "head": "edcf3a3a2e87f335c268460ca0d4b3b23f4a9847",
    "branch": "agent3/rf30-rf31-typed-mean-correction-20260922",
    "source_path": "src/openai_ns_reconstruction/kokuno_rf30_rf31_typed_mean_correction.py",
    "source_blob": "fa31fdf12adceaa5b295db502cf0ccaaad778307",
    "parent_pr": 1118,
    "parent_head": "922f7aa10460ded44af212d313eceff33a2ac647",
    "stage": "RF30-RF31-typed-formula-executor-mechanics-only",
    "rf30_rf31_formula_executor_materialized": True,
    "current_i4_source_fixed_q_backend_materialized": False,
    "current_i4_candidate_defect_evidence": False,
    "correction_coefficients_solved": False,
    "correction_applied": False,
    "nonlinear_remainder_recomputed": False,
    "repository_tests_run": 35688277879,
    "dedicated_run": 35688277940,
    "actions_status_at_registration": "queued",
}

AGENT4_CURRENT_I4_STRESS_AUDIT = {
    "pr": 1119,
    "head": "3f3fd169ce79a2e44b90b50353415d057119e63f",
    "branch": "kokuno-agent4/current-i4-radial-stress-independent-audit-115",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_i4_radial_stress_independent_audit.py",
    "source_blob": "97387bcab05098ecc622372ec3aa08150d79b89d",
    "audited_agent3_pr": 1109,
    "stage": "current-I4-compact-stress-independent-audit",
    "present": True,
    "registered": True,
    "scientifically_admitted": False,
    "audits_radial_force": False,
    "repository_tests_run": 35684778994,
    "dedicated_run": 35684779061,
    "actions_status_at_registration": "queued",
}

AGENT4_RADIAL_FORCE_AUDIT_PENDING = {
    "task": "K4-VAL-116",
    "status": "claimed-pending-delivery",
    "exact_base_pr": 1118,
    "exact_base_head": "922f7aa10460ded44af212d313eceff33a2ac647",
    "claimed_operator": "implementation-distinct centered five-point FD4 reconstruction of partial_z sigma_1",
    "seed": 9173881,
    "time": 0.39,
    "heldout_z_centers": [-0.137, -0.071, 0.061],
    "physical_z_steps": [0.02, 0.01, 0.005],
    "positive_radius_count": 169,
    "positive_radius_interval": [0.005, 0.44],
    "delivered_pr": None,
    "delivered_head": None,
    "present": False,
    "registered": False,
    "scientifically_admitted": False,
}

CORRECTED_SOURCE = {
    "repository": "KokunoYumeto/yang-mills-interacting-workbench",
    "commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
    "path": "navier-stokes/navier_stokes_workbench.tex",
    "blob": "205a99807302e21a51c5eaf223390c0dfc42bcd0",
    "release_date": "2026-09-09",
    "role": "structural provenance only; not paper-exact candidate evidence",
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
        "current_lineage_J_entry_materialized": True,
        "current_cartesian_end_compensation_composed": False,
        "source_exact_main_pulse_amplitude_materialized": False,
        "terminal_global_leading_velocity_materialized": False,
        "source_physicalized_localized_harmonic_materialized": True,
        "source_physicalized_harmonic_self_contained_velocity_xyzt_provider": False,
        "latest_self_contained_project_composite_is_agent2_1117": True,
        "matching_agent4_xi11_composite_audit_present": False,
        "current_i4_radial_force_materialized": True,
        "agent4_matching_current_i4_radial_stress_audit_present": True,
        "agent4_matching_current_i4_radial_stress_audit_registered": True,
        "agent4_matching_current_i4_radial_stress_audit_admitted": False,
        "agent4_current_i4_radial_force_audit_claimed": True,
        "agent4_matching_current_i4_radial_force_audit_present": False,
        "agent4_matching_current_i4_radial_force_audit_registered": False,
        "agent4_matching_current_i4_radial_force_audit_admitted": False,
        "rf30_rf31_formula_executor_materialized": True,
        "current_i4_source_fixed_q_backend_materialized": False,
        "current_i4_rf30_defect_materialized": False,
        "current_i4_rf31_five_row_system_materialized": False,
        "correction_coefficients_solved": False,
        "correction_applied": False,
        "rf44_rf49_nonlinear_remainder_recomputed": False,
        "scoped_current_i4_force_authorized_as_complete_ns_correction_target": False,
        "matched_cartesian_pressure_gradient_materialized": False,
        "preregistered_restricted_forcing_materialized": False,
        "restricted_forcing_proved_not_residual_defined": False,
        "complete_identity_bound_ns_defect_materialized": False,
        "cartesian_correction_velocity_materialized": False,
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


def _canonical_bytes(value: Mapping[str, Any]) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _payload_without_digest(registration: Mapping[str, Any]) -> dict[str, Any]:
    payload = copy.deepcopy(dict(registration))
    payload.pop("registration_sha256", None)
    return payload


def registration_sha256(registration: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_bytes(_payload_without_digest(registration))).hexdigest()


def build_registration() -> dict[str, Any]:
    registration: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "role": "Agent-5 integration/provenance/artifact/CI routing only",
        "parent_a5": copy.deepcopy(PARENT_A5),
        "corrected_source": copy.deepcopy(CORRECTED_SOURCE),
        "frontiers": {
            "leading_end_bookkeeping": copy.deepcopy(AGENT1_CURRENT_PULSE_ENTRY_MJ),
            "latest_self_contained_project_composite": copy.deepcopy(LATEST_SELF_CONTAINED_PROJECT_COMPOSITE),
            "source_oscillatory_physicalization_sibling": copy.deepcopy(AGENT2_SOURCE_PHYSICALIZED_HARMONIC),
            "current_i4_radial_force": copy.deepcopy(CURRENT_I4_RADIAL_FORCE),
            "rf30_rf31_typed_executor_sibling": copy.deepcopy(AGENT3_RF30_RF31_TYPED_EXECUTOR),
            "current_i4_stress_validator": copy.deepcopy(AGENT4_CURRENT_I4_STRESS_AUDIT),
            "current_i4_radial_force_validator_pending": copy.deepcopy(AGENT4_RADIAL_FORCE_AUDIT_PENDING),
        },
        "identity_firewall": {
            "agent1_1124_not_consumed_by_agent2_1117": True,
            "agent2_1125_not_a_replacement_for_agent2_1117_project_composite": True,
            "agent3_1126_mechanics_not_current_i4_candidate_evidence": True,
            "agent4_1119_stress_evidence_not_radial_force_evidence": True,
            "agent4_pending_claim_not_evidence": True,
            "cross_identity_evidence_transfer_allowed": False,
        },
        "core_state": copy.deepcopy(CORE_STATE),
        "frozen_gates": copy.deepcopy(FROZEN_GATES),
        "truth_boundary": _truth_boundary(),
        "shortest_next_closure": [
            "A1 compose exact #1124 M/J inputs through #1116 c1/c2 into the same Cartesian leading identity, then continue terminal/global leading",
            "A4 deliver K4-VAL-116 on exact A3 #1118 radial force; only an exact delivered head may become registered scoped evidence",
            "A2 or source backend supply exact current-I4 fixed-Q covariance/derivative state before A3 #1126 can produce current-I4 RF30/RF31 candidate evidence",
            "global velocity + matched pressure/grad-p + preregistered non-residual-defined restricted forcing before complete NS defect",
            "only then real Cartesian correction/finite cycle and independent A4 held-out full-NS 1e-3 gate",
        ],
    }
    registration["registration_sha256"] = registration_sha256(registration)
    validate_registration(registration)
    return registration


def validate_registration(registration: Mapping[str, Any]) -> None:
    reg = dict(registration)
    if reg.get("schema") != SCHEMA or reg.get("task") != TASK:
        raise ValueError("Agent-5 schema/task drift")
    if registration_sha256(reg) != reg.get("registration_sha256"):
        raise ValueError("Agent-5 registration digest mismatch")

    for name, upstream in (
        ("parent", reg["parent_a5"]),
        ("a1", reg["frontiers"]["leading_end_bookkeeping"]),
        ("composite", reg["frontiers"]["latest_self_contained_project_composite"]),
        ("a2-source", reg["frontiers"]["source_oscillatory_physicalization_sibling"]),
        ("i4-force", reg["frontiers"]["current_i4_radial_force"]),
        ("a3-rf31", reg["frontiers"]["rf30_rf31_typed_executor_sibling"]),
        ("a4-stress", reg["frontiers"]["current_i4_stress_validator"]),
    ):
        head = upstream.get("head")
        if not isinstance(head, str) or not _HEX40.fullmatch(head):
            raise ValueError(f"{name} exact head is not frozen")

    pending = reg["frontiers"]["current_i4_radial_force_validator_pending"]
    if pending != AGENT4_RADIAL_FORCE_AUDIT_PENDING:
        raise ValueError("pending A4 K4-VAL-116 claim drift")
    if pending["delivered_head"] is not None or pending["registered"] or pending["present"]:
        raise ValueError("a pending A4 claim cannot be promoted to delivered evidence")

    if reg["identity_firewall"] != build_identity_firewall():
        raise ValueError("identity firewall drift")
    if reg["core_state"] != CORE_STATE:
        raise ValueError("core readiness drift")
    if reg["frozen_gates"] != FROZEN_GATES:
        raise ValueError("frozen project gates drift")
    if reg["truth_boundary"] != _truth_boundary():
        raise ValueError("scientific truth boundary drift")


def build_identity_firewall() -> dict[str, bool]:
    return {
        "agent1_1124_not_consumed_by_agent2_1117": True,
        "agent2_1125_not_a_replacement_for_agent2_1117_project_composite": True,
        "agent3_1126_mechanics_not_current_i4_candidate_evidence": True,
        "agent4_1119_stress_evidence_not_radial_force_evidence": True,
        "agent4_pending_claim_not_evidence": True,
        "cross_identity_evidence_transfer_allowed": False,
    }


def write_registration(path: str | Path) -> dict[str, Any]:
    reg = build_registration()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(reg, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return reg


def load_registration(path: str | Path) -> dict[str, Any]:
    reg = json.loads(Path(path).read_text(encoding="utf-8"))
    validate_registration(reg)
    return reg


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    reg = write_registration(args.output)
    print(json.dumps({"output": args.output, "registration_sha256": reg["registration_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    _main()
