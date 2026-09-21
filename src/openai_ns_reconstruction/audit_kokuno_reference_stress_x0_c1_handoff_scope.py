"""Fail-closed CR002 audit for the Kokuno X0 reference-stress handoff scope.

This audit governs one narrow representation boundary.  Agent 1 #918 registers
an exact *value* handoff at X0 from the natural stress-free segment and then
propagates a repository-autonomous pressure-driven right branch.  Agent 4 #921
independently checks radial derivatives of that right branch and exact-X0 value
composition.  Neither operation, by itself, proves equality of a natural-side
left derivative and the autonomous continuation-side right derivative.

The audit therefore prevents a C0/value handoff from being relabeled as a C1
or source-prepared stress splice.  It changes no candidate mathematics.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_autonomous_reference_stress_ingest_contract import (
    AGENT1_AUTONOMOUS_REFERENCE_STRESS,
    AGENT4_AUTONOMOUS_REFERENCE_STRESS_AUDIT,
    FINAL_GATE,
    READINESS,
)

SCHEMA = "cr002-kokuno-reference-stress-x0-c1-handoff-scope-v1"
TASK_ID = "CR002-KOKUNO-REFERENCE-STRESS-X0-C1-HANDOFF-087"
SCOPE_PATH = Path("configs/kokuno_reference_stress_x0_c1_handoff_scope.json")
CANONICAL_CONSTRAINTS_PATH = Path("configs/constraints.json")

EXPECTED_PARENT = {
    "pr": 922,
    "head": "8bb1139f2b2d9824fd3bbac693ffa4bf2e71dea2",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_autonomous_reference_stress_ingest_contract.py",
    "source_blob": "733f76ab40bd48d916301ea9ed2d4c8f9e1eff09",
}
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"
EXPECTED_INTEGRATION = {
    "branch": "codex/cr001-constraints",
    "head": "ca37d25d19b97d893030194ebd6364160ae4355e",
}
EXPECTED_PROVENANCE_KEYS = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
EXPECTED_PUBLIC_FACTS = [
    "As registered by Agent 1 #918, the regular axial stress satisfies D_X N_s+N_s=S_n with n_s=N_s/L.",
    "As registered by Agent 1 #918, the natural stress-free segment supplies the value relation n_s(X_0,eta)=-2*U_X(X_0,eta).",
]
EXPECTED_GOVERNED_DISTINCTIONS = {
    "x0_value_handoff_registered": True,
    "right_branch_radial_derivative_audit_registered": True,
    "x0_value_equality_implies_c1_splice": False,
    "right_branch_derivative_consistency_implies_left_right_derivative_match": False,
    "autonomous_pressure_ode_implies_source_prepared_stress_derivative": False,
    "x0_c1_reference_stress_handoff_verified": False,
    "source_stress_free_derivative_compatibility_verified": False,
    "source_prepared_reference_stress_identity_verified": False,
    "autonomous_reference_stress_authorized_as_actual_fixed_kappa_stress": False,
    "autonomous_reference_stress_authorized_as_agent3_correction_target": False,
}
EXPECTED_SCIENTIFIC_STATE = {
    "canonical_eq45_velocity_delivery_remains_independent": True,
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
    return hashlib.sha1(header + data).hexdigest()  # noqa: S324 - Git object identity is SHA-1.


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload


def _canonical_constraints_snapshot(constraints: Mapping[str, Any]) -> dict[str, Any]:
    validation = constraints["validation"]
    thresholds = validation["thresholds"]
    return {
        "nu": constraints["nu"],
        "physical_domain": constraints["domain"]["physical"],
        "evaluation_box": constraints["domain"]["evaluation_box"],
        "support": constraints["domain"]["support"],
        "time_interval": constraints["domain"]["time_interval"],
        "forcing_mode": constraints["forcing"]["mode"],
        "reference_energy": constraints["nontriviality"]["reference_energy"],
        "reference_energy_abs_tolerance": constraints["nontriviality"][
            "reference_energy_abs_tolerance"
        ],
        "validation_seed": validation["seed"],
        "held_out_points": validation["held_out_points"],
        "derivative_steps": validation["derivative_steps"],
        "quadrature_orders_per_axis": validation["quadrature_orders_per_axis"],
        "momentum_max_gate": thresholds["pde_residual_max"],
        "momentum_l2_gate": thresholds["pde_residual_L2"],
        "divergence_max_gate": thresholds["divergence_max"],
        "divergence_l2_gate": thresholds["divergence_L2"],
    }


def mechanics_witness(scope: Mapping[str, Any]) -> dict[str, float | bool]:
    """Replay the mechanics-only equal-value/different-slope counterexample."""
    witness = scope["mechanics_witness"]
    _require(
        witness.get("classification") == "autonomous_mechanics_only",
        "mechanics witness provenance was promoted",
    )
    x0 = float(witness["x0"])
    common = float(witness["common_value"])
    left_slope = float(witness["left_slope"])
    right_slope = float(witness["right_slope"])

    def left(x: float) -> float:
        return common + left_slope * (x - x0)

    def right(x: float) -> float:
        return common + right_slope * (x - x0)

    value_jump = right(x0) - left(x0)
    derivative_jump = right_slope - left_slope
    _require(abs(value_jump) <= 1.0e-15, "mechanics witness lost exact value handoff")
    _require(abs(derivative_jump) > 1.0e-12, "mechanics witness lost derivative mismatch")
    _require(
        abs(derivative_jump - float(witness["derivative_jump"])) <= 2.0e-15,
        "mechanics witness derivative jump drift",
    )
    _require(witness.get("equal_value_at_x0") is True, "witness equality flag drift")
    return {
        "left_value_at_x0": left(x0),
        "right_value_at_x0": right(x0),
        "value_jump": value_jump,
        "left_derivative": left_slope,
        "right_derivative": right_slope,
        "derivative_jump": derivative_jump,
        "equal_value": True,
        "different_derivative": True,
    }


def audit_scope(
    scope: Mapping[str, Any], *, repository_root: str | Path
) -> dict[str, Any]:
    """Audit one supplied scope object against exact upstream and CR001 identities."""
    root = Path(repository_root)
    c = copy.deepcopy(dict(scope))
    _require(c.get("schema") == SCHEMA, "scope schema drift")
    _require(c.get("task_id") == TASK_ID, "scope task id drift")
    _require(c.get("active_integration_reference") == EXPECTED_INTEGRATION, "integration reference drift")
    _require(c.get("audited_parent") == EXPECTED_PARENT, "audited parent identity drift")

    parent_path = root / EXPECTED_PARENT["source_path"]
    parent_bytes = parent_path.read_bytes()
    _require(
        _git_blob_sha(parent_bytes) == EXPECTED_PARENT["source_blob"],
        "A5 #922 source blob drift",
    )

    constraints_path = root / CANONICAL_CONSTRAINTS_PATH
    constraints_bytes = constraints_path.read_bytes()
    _require(
        _git_blob_sha(constraints_bytes) == EXPECTED_CONSTRAINTS_BLOB,
        "canonical CR001 constraints blob drift",
    )
    constraints = json.loads(constraints_bytes.decode("utf-8"))
    frozen = c.get("canonical_constraints", {})
    _require(frozen.get("path") == str(CANONICAL_CONSTRAINTS_PATH), "constraints path drift")
    _require(frozen.get("blob") == EXPECTED_CONSTRAINTS_BLOB, "constraints identity drift")
    observed_constraints = _canonical_constraints_snapshot(constraints)
    for key, value in observed_constraints.items():
        _require(frozen.get(key) == value, f"CR001 constraint drift: {key}")
    _require(
        frozen.get("residual_defined_pointwise_force_forbidden") is True,
        "residual-defined pointwise forcing must remain forbidden",
    )
    _require(frozen.get("collapsed_candidate_forbidden") is True, "candidate collapse must remain forbidden")
    _require(
        frozen.get("post_hoc_threshold_relaxation_forbidden") is True,
        "post-hoc threshold relaxation must remain forbidden",
    )

    pins = c.get("upstream_pins", {})
    a1 = pins.get("agent1_reference_stress", {})
    _require(a1.get("pr") == 918, "Agent-1 PR drift")
    _require(a1.get("head") == AGENT1_AUTONOMOUS_REFERENCE_STRESS["head"], "Agent-1 head drift")
    _require(
        a1.get("source_blob") == AGENT1_AUTONOMOUS_REFERENCE_STRESS["source_blob"],
        "Agent-1 source blob drift",
    )
    _require(
        a1.get("registered_value_handoff")
        == AGENT1_AUTONOMOUS_REFERENCE_STRESS["stress_free_value_handoff"],
        "Agent-1 X0 value-handoff semantics drift",
    )
    _require(
        a1.get("right_branch_ode") == AGENT1_AUTONOMOUS_REFERENCE_STRESS["source_ode"],
        "Agent-1 right-branch ODE semantics drift",
    )
    _require(
        a1.get("pressure_provenance") == AGENT1_AUTONOMOUS_REFERENCE_STRESS["pressure_provenance"],
        "Agent-1 autonomous pressure provenance drift",
    )

    a4 = pins.get("agent4_reference_stress_audit", {})
    _require(a4.get("pr") == 921, "Agent-4 PR drift")
    _require(a4.get("head") == AGENT4_AUTONOMOUS_REFERENCE_STRESS_AUDIT["head"], "Agent-4 head drift")
    _require(
        a4.get("source_blob") == AGENT4_AUTONOMOUS_REFERENCE_STRESS_AUDIT["source_blob"],
        "Agent-4 source blob drift",
    )
    _require(
        a4.get("registered_reference_surface")
        == AGENT4_AUTONOMOUS_REFERENCE_STRESS_AUDIT["scientific_input_surface"],
        "Agent-4 public-values boundary drift",
    )
    _require("FD6" in a4.get("registered_derivative_audit", ""), "Agent-4 derivative-audit identity drift")
    _require(a4.get("exact_x0_check") == "value composition only", "X0 audit scope drift")

    provenance = c.get("provenance_classes", {})
    _require(set(provenance) == EXPECTED_PROVENANCE_KEYS, "four-way provenance classes drift")
    _require(provenance.get("public_source_fact") == EXPECTED_PUBLIC_FACTS, "public-source fact laundering/drift")
    _require(
        any("value handoff" in item.lower() for item in provenance.get("user_requirement", [])) is False,
        "source value handoff was laundered into a user requirement",
    )
    _require(
        any("witness" in item.lower() for item in provenance.get("public_source_fact", [])) is False,
        "autonomous mechanics witness was laundered into public-source fact",
    )
    _require(
        any("left radial derivative" in item for item in provenance.get("pending_unknown", [])),
        "pending left/right derivative compatibility was dropped",
    )

    _require(
        c.get("governed_distinctions") == EXPECTED_GOVERNED_DISTINCTIONS,
        "unsupported C1/source-stress promotion or scope drift",
    )
    _require(
        c.get("delivery_and_scientific_state") == EXPECTED_SCIENTIFIC_STATE,
        "delivery/PDE/visual/exactness state coupling or promotion",
    )

    # Cross-check A5's current fail-closed registration rather than trusting only
    # this new JSON scope.
    _require(READINESS.get("leading_ready") is False, "A5 leading readiness unexpectedly promoted")
    _require(READINESS.get("velocity_export_ready") is False, "A5 Kokuno velocity export unexpectedly promoted")
    _require(READINESS.get("pde_validated") is False, "A5 PDE state unexpectedly promoted")
    _require(FINAL_GATE.get("momentum_sampled_max") == 1.0e-3, "A5 momentum max gate drift")
    _require(FINAL_GATE.get("momentum_volume_l2") == 1.0e-3, "A5 momentum L2 gate drift")
    _require(FINAL_GATE.get("divergence_sampled_max") == 1.0e-5, "A5 divergence max gate drift")
    _require(FINAL_GATE.get("divergence_volume_l2") == 1.0e-5, "A5 divergence L2 gate drift")

    witness_receipt = mechanics_witness(c)
    receipt = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "audited_parent_head": EXPECTED_PARENT["head"],
        "audited_parent_source_blob": EXPECTED_PARENT["source_blob"],
        "canonical_constraints_blob": EXPECTED_CONSTRAINTS_BLOB,
        "mechanics_witness": witness_receipt,
        "x0_value_handoff_registered": True,
        "right_branch_derivative_audit_registered": True,
        "x0_c1_reference_stress_handoff_verified": False,
        "source_stress_free_derivative_compatibility_verified": False,
        "source_prepared_reference_stress_identity_verified": False,
        "pde_validated": False,
        "canonical_eq45_velocity_delivery_remains_independent": True,
    }
    encoded = json.dumps(receipt, sort_keys=True, separators=(",", ":"), allow_nan=False)
    receipt["receipt_sha256"] = hashlib.sha256(encoded.encode("utf-8")).hexdigest()
    return receipt


def audit_repository(repository_root: str | Path) -> dict[str, Any]:
    root = Path(repository_root)
    scope = _load_json(root / SCOPE_PATH)
    return audit_scope(scope, repository_root=root)


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
