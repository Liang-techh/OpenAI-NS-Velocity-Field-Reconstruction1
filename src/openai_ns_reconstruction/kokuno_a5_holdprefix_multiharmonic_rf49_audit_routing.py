"""Fail-closed Kokuno Agent-5 routing after hold-prefix/A2-diagnostic/RF49/A4-prefix advances.

This module is integration/provenance glue only. It does not rewrite Agent-1--4
mathematics, and it never transfers evidence between different candidate
identities.

Fresh frontiers frozen here:

* Agent 1 #1148 exposes the full *unedited* eta-independent exponential hold as
  a low-dimensional helper, while the unified Cartesian leading candidate is
  extended only through a conservative pre-relative-swirl prefix. The two
  public relative-swirl bumps, the complete terminal hold, exterior heat, and
  global leading field remain unmaterialized.
* Agent 2 #1149 adds an implementation-distinct three-resolution Cartesian
  divergence/vorticity diagnostic for the provider-driven bounded
  multi-harmonic complete-curl family of #1141. The source provider remains
  external, so this does not create a self-contained project candidate or a
  complete Navier--Stokes residual.
* Agent 3 #1150 adds a source-exponent-gain handoff firewall on the RF44--RF49
  mechanics route. Its 0.39999/0.17/0.89996 numbers are asymptotic
  class-exponent gains, not multiplicative raw residual contractions. The
  real same-identity current-I4 fixed-Q RF30 backend is still missing.
* Agent 4 #1152 independently audits exact Agent-1 #1148 pre-relative-swirl
  prefix divergence after save/load. Its exact-head CI is still unresolved at
  registration, so this records presence/registration only, never PASS or
  scientific admission. Older Agent-4 #1143 remains bound only to Agent-1
  #1140 and is retained as parent-stage evidence, not transferred forward.

The project gates remain frozen: normalized momentum sampled-max and volume-L2
<= 1e-3; normalized divergence sampled-max and volume-L2 <= 1e-5; canonical
quadrature [24,48,96]; residual-defined/free forcing forbidden.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a5_postpulse_multiharmonic_xi11_audit_routing as parent_a5

SCHEMA = "kokuno-a5-holdprefix-multiharmonic-rf49-audit-routing-v1"
TASK = "KOKUNO-A5-HOLDPREFIX-MULTIHARMONIC-RF49-AUDIT-ROUTING-117"

PARENT_A5 = {
    "pr": 1144,
    "head": "69be0e10c53cb6e2a369be09246aabffc8969773",
    "branch": "kokuno-agent5/postpulse-multiharmonic-xi11-audit-routing-116",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_postpulse_multiharmonic_xi11_audit_routing.py",
    "source_blob": "2f61417d9088c8ef1ffb763345a1ab1146b6c756",
    "task": "KOKUNO-A5-POSTPULSE-MULTIHARMONIC-XI11-AUDIT-ROUTING-116",
    "schema": "kokuno-a5-postpulse-multiharmonic-xi11-audit-routing-v1",
}

AGENT1_HOLD_PREFIX = {
    "pr": 1148,
    "head": "5bafa4199ef127bf50bce30a58cc69ded99e1aa9",
    "branch": "codex/kokuno-a1-postpulse-eta-independent-hold",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_postpulse_eta_independent_hold_prefix.py",
    "source_blob": "b2fd8154fc5242b273879e060ad9a19d8d4c35b8",
    "parent_pr": 1140,
    "parent_head": "c5442c11165d7f0976889b826bf58dce38a3d3dd",
    "stage": "current-cartesian-leading-postpulse-pre-relative-swirl-hold-prefix",
    "public_hold_formula": "30*log(1/lambda)",
    "inherited_lambda": 0.05,
    "unedited_hold_baseline_materialized": True,
    "bounded_velocity_xyzt_materialized": True,
    "pre_relative_swirl_guard_logx_units": 3.5,
    "relative_swirl_bumps_materialized": False,
    "source_terminal_hold_after_eta_flattening_materialized": False,
    "exterior_heat_global_leading_materialized": False,
    "terminal_global_leading_velocity_materialized": False,
    "repository_tests_run": 35701267149,
    "dedicated_run": 35701267153,
    "actions_status_at_registration": "queued",
}

AGENT2_MULTIHARMONIC_DIAGNOSTIC = {
    "task": "K2-OSC-105",
    "pr": 1149,
    "head": "db7a9a4efb4faed2fa46d321fdf5a66b8686de1f",
    "branch": "kokuno-agent2-source-multiharmonic-diagnostics-105",
    "source_path": "src/openai_ns_reconstruction/kokuno_source_bounded_multiharmonic_diagnostics.py",
    "source_blob": "1a6d07535471d164ae53e22cac560091fe5b7f85",
    "audited_agent2_pr": 1141,
    "audited_agent2_head": "4241b0b4fcdc36e750798133101baaf803ba8c44",
    "stage": "provider-driven-bounded-multiharmonic-three-resolution-diagnostic",
    "source_chart_equivalent_fd_steps": [0.02, 0.01, 0.005],
    "divergence_reported": True,
    "vorticity_morphology_reported": True,
    "offgrid_velocity_reported": True,
    "axis_core_zero_checked": True,
    "outer_support_boundary_assessed": False,
    "source_provider_self_contained": False,
    "complete_ns_residual_assessed": False,
    "oscillation_before_after_ns_residual_compared": False,
    "repository_tests_run": 35700978600,
    "dedicated_run": 35700978573,
    "actions_status_at_registration": "queued",
}

AGENT3_RF49_GAIN_FIREWALL = {
    "task": "KOKUNO-A3-RF49-SOURCE-EXPONENT-GAIN-FIREWALL-125",
    "pr": 1150,
    "head": "e0b48ccfa964db0caa9d657961d451a21e1d53bf",
    "branch": "agent3/rf49-source-exponent-gain-firewall-20260922",
    "source_path": "src/openai_ns_reconstruction/kokuno_rf49_source_exponent_gain_firewall.py",
    "source_blob": "baa5cab26ebdeaf8ae6644f452f4f53f6a204881",
    "parent_pr": 1142,
    "parent_head": "d4ff4357b1619cc13eeacfc88d453f73976fc17c",
    "stage": "rf49-source-exponent-gain-handoff-firewall",
    "kappa_s": 1.0e-5,
    "source_required_exponent_increment": 0.1,
    "wave_exponent_gain": 0.39999,
    "retained_tangential_mean_exponent_gain": 0.17,
    "defect_exponent_gain": 0.89996,
    "gain_role": "asymptotic-class-exponent-gain-not-raw-residual-contraction",
    "autonomous_defect_norm_nonexpansion_ratio": 1.0,
    "autonomous_defect_norm_roundoff_slack": 1.0e-12,
    "current_i4_source_chart_backend_materialized": False,
    "current_i4_candidate_gain_evidence": False,
    "cartesian_correction_velocity_materialized": False,
    "complete_ns_defect_materialized": False,
    "repository_tests_run": 35701295350,
    "dedicated_run": 35701295266,
    "actions_status_at_registration": "queued",
}

AGENT4_POSTPULSE_PARENT_LEADING_AUDIT = {
    "task": "K4-VAL-118",
    "pr": 1143,
    "head": "70897e18df4b58b2e023a8bfc82dfdd552db62b9",
    "branch": "kokuno-agent4/postpulse-leading-divergence-118",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_postpulse_eta_flattening_leading_divergence_independent_audit.py",
    "source_blob": "5b90b13eab506e6fe63fe360352d570a5c2ac3c1",
    "audited_agent1_pr": 1140,
    "audited_agent1_head": "c5442c11165d7f0976889b826bf58dce38a3d3dd",
    "stage": "postpulse-eta-flattening-leading-only-independent-divergence-audit",
    "fd_steps": [0.02, 0.01, 0.005],
    "seed": 9173901,
    "save_load_roundtrip": True,
    "implementation_distinct": True,
    "binary64_fd_representability_firewall": True,
    "present": True,
    "registered": True,
    "scoped_gate_passed": None,
    "scientifically_admitted": False,
    "complete_ns_momentum_residual_assessed": False,
    "pde_validated": False,
    "repository_tests_run": 35696887163,
    "dedicated_run": 35696886910,
    "actions_status_at_registration": "queued",
}

AGENT4_HOLD_PREFIX_LEADING_AUDIT = {
    "task": "K4-VAL-119",
    "pr": 1152,
    "head": "caaf740af357fd6741da5118a1254c2b945caad0",
    "branch": "kokuno-agent4/postpulse-hold-prefix-divergence-119",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_postpulse_hold_prefix_leading_divergence_independent_audit.py",
    "source_blob": "89da19c5d1fcd3335ba97b6e5e21f9c8bb472721",
    "audited_agent1_pr": 1148,
    "audited_agent1_head": "5bafa4199ef127bf50bce30a58cc69ded99e1aa9",
    "stage": "pre-relative-swirl-hold-prefix-leading-only-independent-divergence-audit",
    "fd_steps": [0.02, 0.01, 0.005],
    "seed": 9173911,
    "save_load_roundtrip": True,
    "implementation_distinct": True,
    "binary64_fd_representability_firewall": True,
    "unedited_full_hold_helper_excluded_from_scientific_derivatives": True,
    "present": True,
    "registered": True,
    "scoped_gate_passed": None,
    "scientifically_admitted": False,
    "leading_only_ns_residual_assessed": False,
    "complete_ns_momentum_residual_assessed": False,
    "pde_validated": False,
    "repository_tests_run": 35701896410,
    "dedicated_run": 35701896363,
    "actions_status_at_registration": "queued",
}

AGENT2_BOUNDED_MULTIHARMONIC = copy.deepcopy(parent_a5.AGENT2_BOUNDED_MULTIHARMONIC)
LATEST_SELF_CONTAINED_PROJECT_COMPOSITE = copy.deepcopy(parent_a5.LATEST_SELF_CONTAINED_PROJECT_COMPOSITE)
AGENT4_XI11_COMPOSITE_DIVERGENCE_AUDIT = copy.deepcopy(parent_a5.AGENT4_XI11_COMPOSITE_DIVERGENCE_AUDIT)
CURRENT_I4_REAL_RADIAL_FORCE = copy.deepcopy(parent_a5.CURRENT_I4_REAL_RADIAL_FORCE)
AGENT4_CURRENT_I4_RADIAL_FORCE_AUDIT = copy.deepcopy(parent_a5.AGENT4_CURRENT_I4_RADIAL_FORCE_AUDIT)
CORRECTED_SOURCE = copy.deepcopy(parent_a5.CORRECTED_SOURCE)
FROZEN_GATES = copy.deepcopy(parent_a5.FROZEN_GATES)
CORE_STATE = copy.deepcopy(parent_a5.CORE_STATE)


def _truth_boundary() -> dict[str, Any]:
    return {
        "postpulse_unedited_hold_baseline_materialized": True,
        "postpulse_pre_relative_swirl_cartesian_prefix_materialized": True,
        "source_relative_swirl_bumps_materialized": False,
        "source_terminal_hold_after_eta_flattening_materialized": False,
        "terminal_global_leading_velocity_materialized": False,
        "matching_latest_holdprefix_agent2_project_composite_materialized": False,
        "matching_latest_holdprefix_agent4_audit_present": True,
        "matching_latest_holdprefix_agent4_audit_registered": True,
        "matching_latest_holdprefix_agent4_audit_scoped_gate_passed": None,
        "matching_latest_holdprefix_agent4_audit_scientifically_admitted": False,
        "agent4_1143_parent_stage_leading_audit_present": True,
        "agent4_1143_scoped_gate_passed": None,
        "bounded_multiharmonic_complete_curl_family_materialized": True,
        "bounded_multiharmonic_three_resolution_diagnostic_materialized": True,
        "bounded_multiharmonic_source_provider_self_contained": False,
        "bounded_multiharmonic_complete_ns_residual_assessed": False,
        "latest_self_contained_project_composite_is_agent2_1117": True,
        "agent4_matching_xi11_composite_audit_present": True,
        "agent4_matching_xi11_composite_audit_scoped_gate_passed": None,
        "rf49_source_exponent_gain_firewall_materialized": True,
        "rf49_gains_are_raw_residual_contraction_factors": False,
        "current_i4_source_chart_backend_materialized": False,
        "current_i4_rf30_defect_materialized": False,
        "current_i4_source_specific_gain_evidence": False,
        "current_i4_radial_force_materialized": True,
        "agent4_matching_current_i4_radial_force_audit_present": True,
        "cartesian_correction_velocity_materialized": False,
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
        "agent1_1148_hold_prefix_not_consumed_by_xi11_agent2_1117": True,
        "agent4_1143_audit_of_agent1_1140_not_evidence_for_agent1_1148": True,
        "agent4_1152_audit_matches_agent1_1148_only": True,
        "agent4_1152_scoped_divergence_not_full_ns_evidence": True,
        "agent2_1149_provider_diagnostic_not_self_contained_or_full_ns_evidence": True,
        "agent3_1150_source_exponent_gain_not_raw_residual_contraction": True,
        "agent3_1150_mechanics_not_current_i4_candidate_gain_evidence": True,
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
            "leading": copy.deepcopy(AGENT1_HOLD_PREFIX),
            "bounded_multiharmonic_family": copy.deepcopy(AGENT2_BOUNDED_MULTIHARMONIC),
            "bounded_multiharmonic_diagnostic": copy.deepcopy(AGENT2_MULTIHARMONIC_DIAGNOSTIC),
            "rf49_gain_firewall": copy.deepcopy(AGENT3_RF49_GAIN_FIREWALL),
            "postpulse_parent_leading_validator": copy.deepcopy(AGENT4_POSTPULSE_PARENT_LEADING_AUDIT),
            "holdprefix_leading_validator": copy.deepcopy(AGENT4_HOLD_PREFIX_LEADING_AUDIT),
            "latest_self_contained_project_composite": copy.deepcopy(LATEST_SELF_CONTAINED_PROJECT_COMPOSITE),
            "xi11_composite_validator": copy.deepcopy(AGENT4_XI11_COMPOSITE_DIVERGENCE_AUDIT),
            "current_i4_real_radial_force": copy.deepcopy(CURRENT_I4_REAL_RADIAL_FORCE),
            "current_i4_radial_force_validator": copy.deepcopy(AGENT4_CURRENT_I4_RADIAL_FORCE_AUDIT),
        },
        "identity_firewall": build_identity_firewall(),
        "truth_boundary": _truth_boundary(),
        "core_state": copy.deepcopy(CORE_STATE),
        "frozen_gates": copy.deepcopy(FROZEN_GATES),
        "shortest_closure": [
            "materialize the two public relative-swirl bumps on the exact Agent-1 hold-prefix identity before calling the terminal hold complete",
            "finish exterior/global leading and matched pressure without relabeling the bounded prefix as global",
            "bind that exact global Agent-1 identity into one self-contained Agent-2 leading-plus-oscillatory candidate; provider-only diagnostics do not close this seam",
            "materialize the real same-identity fixed-Q RF30 covariance/derivative backend before RF49 mechanics can become candidate-specific correction evidence",
            "bind preregistered non-residual-defined restricted forcing before forming a complete NS defect",
            "only then run the real finite correction cycle and Agent-4 held-out full-NS 1e-3 gate",
        ],
    }
    registration["registration_sha256"] = registration_sha256(registration)
    validate_registration(registration)
    return registration


def validate_registration(registration: Mapping[str, Any]) -> None:
    reg = dict(registration)
    if reg.get("registration_sha256") != registration_sha256(reg):
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
        "leading": AGENT1_HOLD_PREFIX,
        "bounded_multiharmonic_family": AGENT2_BOUNDED_MULTIHARMONIC,
        "bounded_multiharmonic_diagnostic": AGENT2_MULTIHARMONIC_DIAGNOSTIC,
        "rf49_gain_firewall": AGENT3_RF49_GAIN_FIREWALL,
        "postpulse_parent_leading_validator": AGENT4_POSTPULSE_PARENT_LEADING_AUDIT,
        "holdprefix_leading_validator": AGENT4_HOLD_PREFIX_LEADING_AUDIT,
        "latest_self_contained_project_composite": LATEST_SELF_CONTAINED_PROJECT_COMPOSITE,
        "xi11_composite_validator": AGENT4_XI11_COMPOSITE_DIVERGENCE_AUDIT,
        "current_i4_real_radial_force": CURRENT_I4_REAL_RADIAL_FORCE,
        "current_i4_radial_force_validator": AGENT4_CURRENT_I4_RADIAL_FORCE_AUDIT,
    }
    if dict(frontiers) != expected_frontiers:
        raise ValueError("frontier identity drift")
    if reg.get("identity_firewall") != build_identity_firewall():
        raise ValueError("identity firewall drift")
    if reg.get("truth_boundary") != _truth_boundary():
        raise ValueError("truth boundary promotion/drift")
    if reg.get("core_state") != CORE_STATE:
        raise ValueError("core readiness promotion/drift")
    if reg.get("frozen_gates") != FROZEN_GATES:
        raise ValueError("frozen project gates changed")


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
    print(json.dumps({
        "task": registration["task"],
        "registration_sha256": registration["registration_sha256"],
        "pde_validated": registration["truth_boundary"]["pde_validated"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
