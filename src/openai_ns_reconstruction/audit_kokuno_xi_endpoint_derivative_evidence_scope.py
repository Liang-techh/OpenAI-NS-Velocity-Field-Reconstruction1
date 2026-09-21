"""Fail-closed CR002 audit for Xi endpoint derivative evidence semantics.

Agent 4 PR #936 independently reconstructs radial derivatives from save/reloaded
public values on interior points of the reference-stress interval.  It also
checks the exact Xi handoff by comparing the handoff payload with direct public
``values(Xi)`` / production ``radial_derivatives(Xi)`` calls.  Those are useful
and different checks, but the derivative part of the Xi replay is not itself an
implementation-distinct endpoint derivative reconstruction.

This module changes no candidate mathematics.  It only prevents an exact
production-to-production handoff replay from being relabeled as independent Xi
endpoint derivative validation.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_reference_stress_xi110_ingest_contract import (
    AGENT1_REFERENCE_STRESS_XI110,
    AGENT4_REFERENCE_STRESS_XI110_AUDIT,
    FINAL_GATE,
    READINESS,
    TRUTH_BOUNDARY,
)

SCHEMA = "cr002-kokuno-xi-endpoint-derivative-evidence-scope-v1"
TASK_ID = "CR002-KOKUNO-XI-ENDPOINT-DERIVATIVE-EVIDENCE-089"
SCOPE_PATH = Path("configs/kokuno_xi_endpoint_derivative_evidence_scope.json")
CONSTRAINTS_PATH = Path("configs/constraints.json")
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"
EXPECTED_PARENT = {
    "pr": 937,
    "head": "44b9349ef07d947d0da7eff053ae2409348374b8",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_reference_stress_xi110_ingest_contract.py",
    "source_blob": "531106d200b4163018f0331537b7bf6706042bd0",
}
EXPECTED_INTEGRATION = {
    "branch": "codex/cr001-constraints",
    "head": "7f6658cac1a8ffea178bec5accf2fa5a784b8249",
}
EXPECTED_PROVENANCE_KEYS = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
EXPECTED_DISTINCTIONS = {
    "interior_implementation_distinct_radial_derivative_audit_registered": True,
    "Xi_public_handoff_replay_registered": True,
    "Xi_handoff_replay_compares_public_values_to_public_values": True,
    "Xi_handoff_replay_compares_derivatives_to_same_production_radial_derivative_surface": True,
    "Xi_handoff_derivative_replay_is_independent_endpoint_derivative_validation": False,
    "Xi_radial_derivatives_independently_verified": False,
    "Xi_value_replay_implies_Xi_derivative_accuracy": False,
    "interior_derivative_accuracy_implies_endpoint_derivative_accuracy": False,
    "reference_Xi_handoff_implies_actual_candidate_Xi_handoff": False,
    "actual_candidate_Xi_derivative_handoff_verified": False,
    "actual_candidate_bridge_materialized": False,
    "kokuno_global_cartesian_velocity_materialized": False,
    "kokuno_velocity_export_ready": False,
    "pde_validated": False,
}
EXPECTED_DELIVERY_STATE = {
    "canonical_eq45_velocity_delivery_remains_independent": True,
    "canonical_eq45_velocity_export_ready": True,
    "kokuno_velocity_export_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _git_blob_sha(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()  # noqa: S324 - Git blob identity.


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload


def _constraints_snapshot(payload: Mapping[str, Any]) -> dict[str, Any]:
    domain = payload["domain"]
    forcing = payload["forcing"]
    nontriviality = payload["nontriviality"]
    validation = payload["validation"]
    thresholds = validation["thresholds"]
    return {
        "nu": payload["nu"],
        "physical_domain": domain["physical"],
        "evaluation_box": domain["evaluation_box"],
        "support": domain["support"],
        "time_interval": domain["time_interval"],
        "forcing_mode": forcing["mode"],
        "reference_energy": nontriviality["reference_energy"],
        "reference_energy_abs_tolerance": nontriviality["reference_energy_abs_tolerance"],
        "validation_seed": validation["seed"],
        "held_out_points": validation["held_out_points"],
        "derivative_steps": validation["derivative_steps"],
        "quadrature_orders_per_axis": validation["quadrature_orders_per_axis"],
        "momentum_max_gate": thresholds["pde_residual_max"],
        "momentum_l2_gate": thresholds["pde_residual_L2"],
        "divergence_max_gate": thresholds["divergence_max"],
        "divergence_l2_gate": thresholds["divergence_L2"],
    }


def mechanics_witness(scope: Mapping[str, Any]) -> dict[str, Any]:
    """Replay the logic-only shared-derivative counterexample."""
    witness = scope["mechanics_witness"]
    _require(
        witness.get("classification") == "autonomous_mechanics_only",
        "mechanics witness provenance was promoted",
    )
    Xi = float(witness["Xi"])
    h = float(witness["backward_step"])
    production = float(witness["production_derivative_at_Xi"])
    handoff = float(witness["handoff_derivative_at_Xi"])
    _require(math.isfinite(Xi) and math.isfinite(h) and h > 0.0, "invalid witness geometry")
    _require(math.isfinite(production) and math.isfinite(handoff), "invalid witness derivative")

    def f(x: float) -> float:
        return x * x

    independent = (3.0 * f(Xi) - 4.0 * f(Xi - h) + f(Xi - 2.0 * h)) / (2.0 * h)
    replay_error = abs(production - handoff)
    independent_error = abs(production - independent)
    _require(replay_error <= 1.0e-15, "shared production/handoff replay must agree exactly")
    _require(abs(independent - 2.0 * Xi) <= 1.0e-12, "three-point backward witness drift")
    _require(independent_error >= 5.0e-2, "witness no longer separates replay from accuracy")
    _require(
        abs(replay_error - float(witness["production_handoff_replay_abs_error"])) <= 1.0e-15,
        "registered replay error drift",
    )
    _require(
        abs(independent - float(witness["independent_endpoint_derivative"])) <= 1.0e-12,
        "registered independent endpoint derivative drift",
    )
    _require(
        abs(independent_error - float(witness["production_vs_independent_abs_error"])) <= 1.0e-12,
        "registered production-vs-independent error drift",
    )
    return {
        "production_derivative": production,
        "handoff_derivative": handoff,
        "production_handoff_replay_abs_error": replay_error,
        "independent_three_point_backward_derivative": independent,
        "production_vs_independent_abs_error": independent_error,
        "replay_consistency_implies_endpoint_accuracy": False,
    }


def audit_scope(scope: Mapping[str, Any], *, repository_root: str | Path) -> dict[str, Any]:
    root = Path(repository_root)
    c = copy.deepcopy(dict(scope))
    _require(c.get("schema") == SCHEMA, "scope schema drift")
    _require(c.get("task_id") == TASK_ID, "task id drift")
    _require(c.get("active_integration_reference") == EXPECTED_INTEGRATION, "active integration reference drift")
    _require(c.get("audited_parent") == EXPECTED_PARENT, "audited A5 parent identity drift")

    parent_bytes = (root / EXPECTED_PARENT["source_path"]).read_bytes()
    _require(_git_blob_sha(parent_bytes) == EXPECTED_PARENT["source_blob"], "A5 #937 source blob drift")

    constraints_bytes = (root / CONSTRAINTS_PATH).read_bytes()
    _require(_git_blob_sha(constraints_bytes) == EXPECTED_CONSTRAINTS_BLOB, "canonical CR001 constraints blob drift")
    constraints = json.loads(constraints_bytes.decode("utf-8"))
    frozen = c.get("canonical_constraints", {})
    _require(frozen.get("path") == str(CONSTRAINTS_PATH), "constraints path drift")
    _require(frozen.get("blob") == EXPECTED_CONSTRAINTS_BLOB, "constraints blob identity drift")
    for key, value in _constraints_snapshot(constraints).items():
        _require(frozen.get(key) == value, f"CR001 constraint drift: {key}")
    _require(frozen.get("residual_defined_pointwise_force_forbidden") is True, "free residual-defined forcing was enabled")
    _require(frozen.get("collapsed_candidate_forbidden") is True, "candidate collapse was enabled")
    _require(frozen.get("post_hoc_threshold_relaxation_forbidden") is True, "post-hoc threshold relaxation was enabled")

    evidence = c.get("registered_upstream_evidence", {})
    a1 = evidence.get("agent1_reference_stress", {})
    _require(a1.get("pr") == AGENT1_REFERENCE_STRESS_XI110["pr"], "A1 PR drift")
    _require(a1.get("head") == AGENT1_REFERENCE_STRESS_XI110["head"], "A1 head drift")
    _require(a1.get("source_blob") == AGENT1_REFERENCE_STRESS_XI110["source_blob"], "A1 source blob drift")
    _require(a1.get("reference_range") == AGENT1_REFERENCE_STRESS_XI110["reference_range"], "A1 reference range drift")
    _require(a1.get("Xi") == AGENT1_REFERENCE_STRESS_XI110["X_i"], "A1 Xi drift")
    _require(a1.get("public_radial_derivatives_registered") is True, "A1 radial derivative surface disappeared")
    _require(a1.get("Xi_handoff_registered") is True, "A1 Xi handoff disappeared")

    a4 = evidence.get("agent4_independent_audit", {})
    _require(a4.get("pr") == AGENT4_REFERENCE_STRESS_XI110_AUDIT["pr"], "A4 PR drift")
    _require(a4.get("head") == AGENT4_REFERENCE_STRESS_XI110_AUDIT["head"], "A4 head drift")
    _require(a4.get("source_blob") == AGENT4_REFERENCE_STRESS_XI110_AUDIT["source_blob"], "A4 source blob drift")
    protocol = AGENT4_REFERENCE_STRESS_XI110_AUDIT["protocol"]
    _require(a4.get("independent_operator") == AGENT4_REFERENCE_STRESS_XI110_AUDIT["independent_operator"], "A4 independent operator drift")
    _require(a4.get("independent_sample_x_range") == protocol["x_range"], "A4 independent X sampling range drift")
    _require(a4.get("physical_halfwidths") == protocol["physical_halfwidths"], "A4 independent derivative ladder drift")
    _require(a4.get("Xi_handoff_abs_gate") == protocol["Xi_handoff_abs_gate"], "A4 Xi handoff replay gate drift")
    _require(a4.get("Xi_direct_derivative_reference") == "field.radial_derivatives(Xi, endpoint_eta)", "Xi direct derivative reference was rewritten")
    _require(a4.get("Xi_handoff_derivative_payload") == "field.handoff_at_Xi(endpoint_eta)", "Xi handoff payload reference was rewritten")
    _require(a4.get("independent_endpoint_derivative_operator_applied_at_Xi") is False, "unregistered independent Xi endpoint operator was asserted")

    provenance = c.get("provenance_classes", {})
    _require(set(provenance) == EXPECTED_PROVENANCE_KEYS, "four-way provenance classes drift")
    _require(provenance.get("public_source_fact") == [], "repository evidence semantics were laundered into public-source fact")
    _require(
        any("endpoint-safe" in item for item in provenance.get("pending_unknown", [])),
        "pending independent Xi endpoint derivative check was dropped",
    )
    _require(
        not any("mechanics_witness" in item or "shared-bug" in item for item in provenance.get("public_source_fact", [])),
        "mechanics witness was laundered into public-source facts",
    )

    _require(c.get("governed_distinctions") == EXPECTED_DISTINCTIONS, "unsupported Xi endpoint evidence promotion or scope drift")
    _require(c.get("delivery_and_scientific_state") == EXPECTED_DELIVERY_STATE, "delivery/PDE/visual/exactness state coupling")

    # Cross-check the exact A5 registration remains fail-closed beyond the scoped reference seam.
    _require(READINESS.get("leading_ready") is False, "A5 leading readiness unexpectedly promoted")
    _require(READINESS.get("velocity_export_ready") is False, "A5 Kokuno velocity export unexpectedly promoted")
    _require(READINESS.get("pde_validated") is False, "A5 PDE state unexpectedly promoted")
    _require(TRUTH_BOUNDARY.get("actual_final_interpolation_100_to_Xi_materialized") is False, "actual candidate bridge unexpectedly promoted")
    _require(TRUTH_BOUNDARY.get("actual_G_i_at_Xi_materialized") is False, "actual G_i unexpectedly promoted")
    _require(TRUTH_BOUNDARY.get("actual_ell_i_at_Xi_materialized") is False, "actual ell_i unexpectedly promoted")
    _require(TRUTH_BOUNDARY.get("corrected_global_cartesian_leading_velocity_materialized") is False, "global Cartesian leading velocity unexpectedly promoted")
    _require(TRUTH_BOUNDARY.get("complete_candidate_api_ready") is False, "complete candidate API unexpectedly promoted")
    _require(FINAL_GATE.get("normalized_momentum_sampled_max") == 1.0e-3, "momentum max gate drift")
    _require(FINAL_GATE.get("normalized_momentum_volume_l2") == 1.0e-3, "momentum L2 gate drift")
    _require(FINAL_GATE.get("divergence_sampled_max") == 1.0e-5, "divergence max gate drift")
    _require(FINAL_GATE.get("divergence_volume_l2") == 1.0e-5, "divergence L2 gate drift")

    witness = mechanics_witness(c)
    receipt = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "audited_parent_head": EXPECTED_PARENT["head"],
        "audited_parent_source_blob": EXPECTED_PARENT["source_blob"],
        "agent4_source_blob": AGENT4_REFERENCE_STRESS_XI110_AUDIT["source_blob"],
        "canonical_constraints_blob": EXPECTED_CONSTRAINTS_BLOB,
        "mechanics_witness": witness,
        "interior_independent_derivative_audit_registered": True,
        "Xi_public_handoff_replay_registered": True,
        "Xi_independent_endpoint_derivative_verified": False,
        "actual_candidate_Xi_derivative_handoff_verified": False,
        "kokuno_velocity_export_ready": False,
        "pde_validated": False,
        "canonical_eq45_velocity_delivery_remains_independent": True,
    }
    encoded = json.dumps(receipt, sort_keys=True, separators=(",", ":"), allow_nan=False)
    receipt["receipt_sha256"] = hashlib.sha256(encoded.encode("utf-8")).hexdigest()
    return receipt


def audit_repository(repository_root: str | Path) -> dict[str, Any]:
    root = Path(repository_root)
    return audit_scope(_load_json(root / SCOPE_PATH), repository_root=root)


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = audit_repository(args.repository_root)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    _main()
