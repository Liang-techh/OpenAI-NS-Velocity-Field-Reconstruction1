"""Fail-closed CR002 audit for RF40 power-law oscillatory-support/evidence scope.

Governance only. This module does not validate Navier--Stokes, visual
correspondence, paper exactness, OpenAI-field identity, or blow-up.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

SCHEMA = "cr002-kokuno-rf40-power-law-oscillatory-support-scope-v1"
TASK_ID = "CR002-KOKUNO-RF40-POWER-LAW-OSCILLATORY-SUPPORT-SCOPE-105"
CONTRACT_REL = Path("configs/kokuno_rf40_power_law_oscillatory_support_scope.json")
CANONICAL_CONSTRAINTS_REL = Path("configs/constraints.json")
CANONICAL_CONSTRAINTS_GIT_BLOB_SHA1 = "6c559e42895a606e2ef025ade4cb448966d75814"

EXPECTED_BASE = {
    "pr": 1012,
    "head": "9b2be0ffecc14479dbd941eebd663d63f89f2a2f",
    "branch": "codex/kokuno-a4-rf40-power-law-composite-divergence-audit-102",
    "source_blob": "f5f7acf0c02f5f0712943ba841b1bc04478bcdd9",
}
EXPECTED_A2 = {
    "pr": 1010,
    "head": "e36d9da4b7f037e998f5b1658f8c0ea291a76b80",
    "source_blob": "4f1dc0566c041b9de20f18b4ae9c49d9fd93c58d",
    "task": "K2-OSC-087",
}
EXPECTED_A4 = {
    "pr": 1012,
    "head": "9b2be0ffecc14479dbd941eebd663d63f89f2a2f",
    "source_blob": "f5f7acf0c02f5f0712943ba841b1bc04478bcdd9",
    "task": "K4-VAL-102",
}
EXPECTED_CR001 = {
    "nu": 0.01,
    "physical_domain": "R^3",
    "evaluation_box": [[-2, 2], [-2, 2], [-2, 2]],
    "support": "r < 2 and abs(z) < 2",
    "time_interval": [0.25, 0.75],
    "forcing_mode": "restricted_two_parameter_family",
    "reference_energy": 1.0,
    "reference_energy_abs_tolerance": 0.001,
    "validation_seed": 914027,
    "held_out_points": 4096,
    "derivative_steps": [0.02, 0.01, 0.005],
    "quadrature_orders_per_axis": [24, 48, 96],
    "momentum_max": 0.001,
    "momentum_l2": 0.001,
    "divergence_max": 1e-5,
    "divergence_l2": 1e-5,
    "residual_defined_free_forcing_forbidden": True,
    "candidate_collapse_forbidden": True,
    "post_hoc_threshold_relaxation_forbidden": True,
}
REQUIRED_FALSE = (
    "composite_callable_means_each_summand_nonzero_everywhere",
    "power_law_composite_label_implies_nontrivial_oscillatory_contribution_on_power_law_interval",
    "power_law_independent_divergence_audit_implies_nontrivial_oscillatory_structure_validated_there",
    "zero_power_law_oscillatory_increment_is_public_source_support_fact",
    "zero_power_law_oscillatory_increment_establishes_visual_correspondence",
    "scoped_divergence_audit_implies_canonical_whole_domain_admission",
    "scoped_divergence_audit_implies_momentum_or_full_ns_validation",
    "callable_save_load_implies_global_project_domain_totality",
)
REQUIRED_TRUE = (
    "power_law_composite_callable_materialized",
    "power_law_composite_identity_preserving_save_load_materialized",
    "power_law_independent_cartesian_divergence_audit_materialized",
    "inner_oscillatory_signal_required_nontrivial",
    "registered_power_law_oscillatory_velocity_required_zero",
)

def _git_blob_sha1(path: Path) -> str:
    raw = path.read_bytes()
    return hashlib.sha1(f"blob {len(raw)}\0".encode("ascii") + raw).hexdigest()

def _load(path: Path) -> Mapping[str, Any]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, Mapping):
        raise AssertionError(f"{path} is not a JSON object")
    return raw

def _assert_equal(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise AssertionError(f"{label} drifted: expected {expected!r}, got {actual!r}")

def _replay_cr001(root: Path, contract: Mapping[str, Any]) -> None:
    path = root / CANONICAL_CONSTRAINTS_REL
    if _git_blob_sha1(path) != CANONICAL_CONSTRAINTS_GIT_BLOB_SHA1:
        raise AssertionError("canonical constraints blob drifted")
    raw = _load(path)
    derived = {
        "nu": raw["nu"],
        "physical_domain": raw["domain"]["physical"],
        "evaluation_box": raw["domain"]["evaluation_box"],
        "support": raw["domain"]["support"],
        "time_interval": raw["domain"]["time_interval"],
        "forcing_mode": raw["forcing"]["mode"],
        "reference_energy": raw["nontriviality"]["reference_energy"],
        "reference_energy_abs_tolerance": raw["nontriviality"]["reference_energy_abs_tolerance"],
        "validation_seed": raw["validation"]["seed"],
        "held_out_points": raw["validation"]["held_out_points"],
        "derivative_steps": raw["validation"]["derivative_steps"],
        "quadrature_orders_per_axis": raw["validation"]["quadrature_orders_per_axis"],
        "momentum_max": raw["validation"]["thresholds"]["pde_residual_max"],
        "momentum_l2": raw["validation"]["thresholds"]["pde_residual_L2"],
        "divergence_max": raw["validation"]["thresholds"]["divergence_max"],
        "divergence_l2": raw["validation"]["thresholds"]["divergence_L2"],
        "residual_defined_free_forcing_forbidden":
            "No residual-dependent basis or pointwise free force" in raw["forcing"]["restriction"],
        "candidate_collapse_forbidden":
            "reject collapsed candidates" in raw["nontriviality"]["enforcement"],
        "post_hoc_threshold_relaxation_forbidden":
            "changing thresholds requires a new experiment version" in raw["validation"]["failure_policy"],
    }
    _assert_equal(derived, EXPECTED_CR001, "canonical CR001 replay")
    snap = dict(contract["cr001_snapshot"])
    snap.pop("constraints_blob_sha1", None)
    _assert_equal(snap, EXPECTED_CR001, "contract CR001 snapshot")
    _assert_equal(contract["cr001_snapshot"]["constraints_blob_sha1"],
                  CANONICAL_CONSTRAINTS_GIT_BLOB_SHA1, "contract constraints blob")

def audit(root: Path) -> None:
    contract = _load(root / CONTRACT_REL)
    _assert_equal(contract.get("schema"), SCHEMA, "schema")
    _assert_equal(contract.get("task_id"), TASK_ID, "task id")
    _assert_equal(contract.get("exact_base"), EXPECTED_BASE, "exact base")

    prov = contract.get("four_way_provenance")
    if not isinstance(prov, Mapping):
        raise AssertionError("four_way_provenance missing")
    _assert_equal(set(prov), {
        "user_requirements", "public_source_facts",
        "autonomous_repository_facts", "pending_or_unknown"
    }, "four-way provenance keys")
    for key, value in prov.items():
        if not isinstance(value, list) or not value:
            raise AssertionError(f"provenance bucket {key} must be nonempty")

    upstream = contract.get("upstream_scope")
    if not isinstance(upstream, Mapping):
        raise AssertionError("upstream_scope missing")
    a2 = upstream["agent2_1010"]
    a4 = upstream["agent4_1012"]
    for k, v in EXPECTED_A2.items():
        _assert_equal(a2.get(k), v, f"A2 {k}")
    for k, v in EXPECTED_A4.items():
        _assert_equal(a4.get(k), v, f"A4 {k}")
    _assert_equal(a2.get("registered_power_law_oscillatory_abs_max_gate"), 0.0,
                  "A2 power-law oscillatory zero gate")
    _assert_equal(a2.get("inner_oscillatory_nontrivial_gate"), True,
                  "A2 inner oscillatory nontrivial gate")
    _assert_equal(a4.get("audited_agent2_pr"), 1010, "A4 audited A2 PR")
    _assert_equal(a4.get("audited_agent2_head"), EXPECTED_A2["head"], "A4 audited A2 head")
    _assert_equal(a4.get("scoped_divergence_evidence_only"), True, "A4 evidence scope")
    _assert_equal(a4.get("momentum_or_full_ns_evidence"), False, "A4 full-NS scope")

    locked = contract.get("machine_locked_distinctions")
    if not isinstance(locked, Mapping):
        raise AssertionError("machine_locked_distinctions missing")
    for key in REQUIRED_TRUE:
        _assert_equal(locked.get(key), True, key)
    for key in REQUIRED_FALSE:
        _assert_equal(locked.get(key), False, key)

    witness = contract.get("scope_logic_witness")
    if not isinstance(witness, Mapping):
        raise AssertionError("scope witness missing")
    _assert_equal(witness.get("classification"), "autonomous_mechanics_only", "witness classification")
    toy = witness.get("toy_values")
    _assert_equal(toy["u_total_at_probe"],
                  [a + b for a, b in zip(toy["u_lead_at_probe"], toy["u_osc_at_probe"])],
                  "toy additive closure")
    if any(float(v) != 0.0 for v in toy["u_osc_at_probe"]):
        raise AssertionError("toy oscillatory witness must be exactly zero")
    _assert_equal(witness.get("not_a_public_source_fact"), True, "witness source boundary")
    _assert_equal(witness.get("not_candidate_numerical_evidence"), True, "witness numerical boundary")

    canonical = contract.get("canonical_delivery_independence")
    expected_canonical = {
        "candidate_family": "eq45_supported_velocity_candidate_v1",
        "velocity_api": "openai_ns_reconstruction.eq45_supported_delivery:velocity",
        "velocity_export_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }
    _assert_equal(canonical, expected_canonical, "canonical Eq45 independence")

    kokuno = contract.get("kokuno_route_state")
    for key in ("velocity_export_ready", "visual_correspondence_verified",
                "pde_validated", "paper_exact", "openai_field_identified",
                "complete_candidate_api_ready"):
        _assert_equal(kokuno.get(key), False, f"Kokuno {key}")

    changes = contract.get("scientific_changes")
    if not isinstance(changes, Mapping) or any(changes.values()):
        raise AssertionError("governance increment must not change scientific bytes or gates")

    forbidden = contract.get("forbidden_promotions")
    if not isinstance(forbidden, list) or len(forbidden) < 10:
        raise AssertionError("forbidden promotion firewall is incomplete")
    needed_tokens = (
        "nontrivial oscillatory contribution",
        "public-source support fact",
        "canonical whole-domain admission",
        "momentum/full-NS evidence",
        "unrestricted pointwise forcing",
        "post-hoc threshold relaxation",
    )
    joined = "\n".join(str(v) for v in forbidden)
    for token in needed_tokens:
        if token not in joined:
            raise AssertionError(f"forbidden-promotion token missing: {token}")

    _replay_cr001(root, contract)

def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args(argv)
    audit(args.root)
    print("CR002 RF40 power-law oscillatory-support scope: PASS")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
