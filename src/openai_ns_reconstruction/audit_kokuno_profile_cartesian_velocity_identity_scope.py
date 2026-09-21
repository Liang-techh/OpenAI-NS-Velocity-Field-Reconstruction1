"""Fail-closed CR002 audit for Kokuno profile-to-Cartesian velocity identity.

Agent 1 #925 materializes and serializes the fixed-kappa source-coordinate
profile surfaces F(X,eta), U(X,eta), E(X,eta). Agent 4 #928 independently
audits those public profile values and their radial derivatives. Agent 5 #929
registers that seam while keeping corrected/global Cartesian leading velocity
and the complete candidate API false.

This audit locks the representation boundary: profile identity and profile
derivative consistency are not, by themselves, a physical
velocity(x,y,z,t)->[u,v,w] artifact identity. It changes no candidate math.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_fixed_kappa_x100_ingest_contract import (
    AGENT1_FIXED_KAPPA_X100,
    AGENT4_FIXED_KAPPA_X100_AUDIT,
    FINAL_GATE,
    READINESS,
    TRUTH_BOUNDARY,
)

SCHEMA = "cr002-kokuno-profile-cartesian-velocity-identity-scope-v1"
TASK_ID = "CR002-KOKUNO-PROFILE-CARTESIAN-IDENTITY-088"
SCOPE_PATH = Path("configs/kokuno_profile_cartesian_velocity_identity_scope.json")
CANONICAL_CONSTRAINTS_PATH = Path("configs/constraints.json")

EXPECTED_PARENT = {
    "pr": 929,
    "head": "744a4fbc642b861f61da91de713dd06adc357f09",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_fixed_kappa_x100_ingest_contract.py",
    "source_blob": "502e64e97903853a16025815eabe31bf8b2111d0",
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
    "As registered by Agent 1 #925 from the pinned corrected public reader, the fixed-kappa continuation uses partial_X log(F)=-kappa_0*p1r/(2X), partial_X U=-kappa_0*nsr/2, and E=sqrt(2X)*F on its stated source-coordinate interval.",
    "As registered by Agent 1 #925, the displayed fixed-kappa shear schedule is kappa_0*(p1r,p2r*E_r/E).",
]
EXPECTED_GOVERNED_DISTINCTIONS = {
    "source_coordinate_profile_executable": True,
    "profile_save_load_semantic_identity_registered": True,
    "profile_radial_derivative_audit_registered": True,
    "profile_identity_implies_cartesian_velocity_identity": False,
    "profile_semantic_sha_is_velocity_candidate_sha": False,
    "profile_radial_derivatives_imply_cartesian_velocity_jacobian": False,
    "profile_values_alone_fix_unique_cartesian_component_map": False,
    "kokuno_cartesian_coordinate_time_pullback_materialized": False,
    "kokuno_cartesian_component_map_bound_to_profile_artifact": False,
    "kokuno_axis_regular_cartesian_velocity_materialized": False,
    "kokuno_global_outer_join_materialized": False,
    "kokuno_velocity_provider_save_load_identity_materialized": False,
    "kokuno_complete_candidate_api_ready": False,
    "kokuno_velocity_export_ready": False,
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
EXPECTED_DELIVERY_REQUIREMENTS = [
    "A callable physical-input API velocity(x,y,z,t)->[u,v,w] with explicit units/domain/time semantics.",
    "An executable coordinate/time pullback from physical (x,y,z,t) into every source-coordinate variable consumed by the exact profile artifact.",
    "An explicit component/scaling map from profile quantities to cylindrical and then Cartesian velocity components, bound to the same profile semantic identity.",
    "A finite axis branch or mathematically defined r->0 limit wherever cylindrical-to-Cartesian conversion is singular.",
    "A complete/global join and support/extension rule covering the declared physical delivery domain rather than only X_1<=X<=100.",
    "Save/load identity and candidate SHA bound to the actual Cartesian velocity provider/artifact, not merely to F/U/E profile configuration.",
    "Independent state flags so exportability cannot promote visual correspondence, PDE validation, paper exactness, or OpenAI-field identity.",
]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _git_blob_sha(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()  # noqa: S324 - Git blob identity is SHA-1.


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload


def _canonical_constraints_snapshot(constraints: Mapping[str, Any]) -> dict[str, Any]:
    domain = constraints["domain"]
    forcing = constraints["forcing"]
    nontriviality = constraints["nontriviality"]
    validation = constraints["validation"]
    thresholds = validation["thresholds"]
    return {
        "nu": constraints["nu"],
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
        "boundary": domain["boundary"],
        "pressure_gauge": domain["pressure_gauge"],
        "forcing_a_bounds": forcing["parameters"]["a"],
        "forcing_c_bounds": forcing["parameters"]["c"],
        "minimum_energy_each_validation_time": nontriviality["minimum_energy_each_validation_time"],
        "maximum_energy_each_validation_time": nontriviality["maximum_energy_each_validation_time"],
        "validation_times": validation["times"],
        "residual_norms": validation["norms"],
    }


def mechanics_witness(scope: Mapping[str, Any]) -> dict[str, Any]:
    """Replay a mechanics-only profile-packet/component-map counterexample."""
    witness = scope["mechanics_witness"]
    _require(
        witness.get("classification") == "autonomous_mechanics_only",
        "mechanics witness provenance was promoted",
    )
    profile = witness["profile_packet"]
    F = float(profile["F"])
    U = float(profile["U"])
    E = float(profile["E"])
    phi = float(witness["azimuth_phi"])
    _require(abs(phi) <= 1.0e-15, "mechanics witness azimuth drifted")

    # Standard cylindrical-to-Cartesian conversion, applied to two deliberately
    # different autonomous assignments of the same scalar packet.
    def cartesian(ur: float, utheta: float, uz: float) -> tuple[float, float, float]:
        return (
            ur * math.cos(phi) - utheta * math.sin(phi),
            ur * math.sin(phi) + utheta * math.cos(phi),
            uz,
        )

    a = cartesian(U, E, F)
    b = cartesian(U, F, E)
    expected_a = tuple(float(v) for v in witness["cartesian_a_at_phi0"])
    expected_b = tuple(float(v) for v in witness["cartesian_b_at_phi0"])
    _require(all(abs(x - y) <= 1.0e-15 for x, y in zip(a, expected_a)), "map-A witness drift")
    _require(all(abs(x - y) <= 1.0e-15 for x, y in zip(b, expected_b)), "map-B witness drift")
    difference = math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))
    _require(difference > 1.0e-12, "mechanics witness lost Cartesian non-uniqueness")
    _require(witness.get("same_profile_packet") is True, "same-profile witness flag drift")
    _require(witness.get("different_cartesian_velocity") is True, "different-velocity witness flag drift")
    return {
        "profile_packet": {"F": F, "U": U, "E": E},
        "cartesian_a": list(a),
        "cartesian_b": list(b),
        "cartesian_difference_l2": difference,
        "same_profile_packet": True,
        "different_cartesian_velocity": True,
    }


def audit_scope(scope: Mapping[str, Any], *, repository_root: str | Path) -> dict[str, Any]:
    """Audit one supplied scope against exact A5 profile-ingest and CR001 identities."""
    root = Path(repository_root)
    c = copy.deepcopy(dict(scope))
    _require(c.get("schema") == SCHEMA, "scope schema drift")
    _require(c.get("task_id") == TASK_ID, "scope task id drift")
    _require(c.get("active_integration_reference") == EXPECTED_INTEGRATION, "integration reference drift")
    _require(c.get("audited_parent") == EXPECTED_PARENT, "audited parent identity drift")

    parent_bytes = (root / EXPECTED_PARENT["source_path"]).read_bytes()
    _require(_git_blob_sha(parent_bytes) == EXPECTED_PARENT["source_blob"], "A5 #929 source blob drift")

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
    for key, value in _canonical_constraints_snapshot(constraints).items():
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
    a1 = pins.get("agent1_fixed_kappa_profile", {})
    _require(a1.get("pr") == 925, "Agent-1 PR drift")
    _require(a1.get("head") == AGENT1_FIXED_KAPPA_X100["head"], "Agent-1 head drift")
    _require(a1.get("source_blob") == AGENT1_FIXED_KAPPA_X100["source_blob"], "Agent-1 source blob drift")
    _require(
        a1.get("registered_public_values") == AGENT1_FIXED_KAPPA_X100["public_values"],
        "Agent-1 public profile-value surface drift",
    )
    _require(
        a1.get("registered_public_radial_derivatives")
        == AGENT1_FIXED_KAPPA_X100["public_radial_derivatives"],
        "Agent-1 radial-derivative surface drift",
    )
    _require(
        a1.get("semantic_identity_and_save_load_registered")
        == AGENT1_FIXED_KAPPA_X100["semantic_identity_and_save_load_registered"],
        "Agent-1 profile serialization identity drift",
    )
    _require(
        a1.get("source_prepared_appendixA_pressure")
        == AGENT1_FIXED_KAPPA_X100["source_prepared_appendixA_pressure"],
        "Agent-1 pressure provenance drift",
    )
    _require(
        a1.get("source_exact_fixed_kappa_continuation")
        == AGENT1_FIXED_KAPPA_X100["source_exact_fixed_kappa_continuation"],
        "Agent-1 source-exactness boundary drift",
    )

    a4 = pins.get("agent4_profile_audit", {})
    _require(a4.get("pr") == 928, "Agent-4 PR drift")
    _require(a4.get("head") == AGENT4_FIXED_KAPPA_X100_AUDIT["head"], "Agent-4 head drift")
    _require(a4.get("source_blob") == AGENT4_FIXED_KAPPA_X100_AUDIT["source_blob"], "Agent-4 source blob drift")
    for key in (
        "saved_reloaded_public_values_only",
        "production_derivatives_used_as_reference",
        "source_coordinate_profile_audit_only",
        "complete_ns_residual_audit",
    ):
        _require(a4.get(key) == AGENT4_FIXED_KAPPA_X100_AUDIT[key], f"Agent-4 audit boundary drift: {key}")

    provenance = c.get("provenance_classes", {})
    _require(set(provenance) == EXPECTED_PROVENANCE_KEYS, "four-way provenance classes drift")
    _require(provenance.get("public_source_fact") == EXPECTED_PUBLIC_FACTS, "public-source fact laundering/drift")
    _require(
        not any(
            "counterexample" in item.lower() or "witness" in item.lower()
            for item in provenance.get("public_source_fact", [])
        ),
        "autonomous mechanics witness was laundered into public-source fact",
    )
    _require(
        any("physical inputs (x,y,z,t)" in item for item in provenance.get("pending_unknown", [])),
        "pending physical coordinate/time pullback was dropped",
    )
    _require(
        any(
            "actual Kokuno Cartesian velocity provider" in item
            for item in provenance.get("pending_unknown", [])
        ),
        "pending Cartesian velocity save/load identity was dropped",
    )

    _require(
        c.get("governed_distinctions") == EXPECTED_GOVERNED_DISTINCTIONS,
        "unsupported profile-to-Cartesian promotion or scope drift",
    )
    _require(
        c.get("cartesian_velocity_delivery_requirements") == EXPECTED_DELIVERY_REQUIREMENTS,
        "Cartesian velocity delivery requirements drift",
    )
    _require(
        c.get("delivery_and_scientific_state") == EXPECTED_DELIVERY_STATE,
        "delivery/PDE/visual/exactness state coupling or promotion",
    )

    # Cross-check the parent ingest contract itself, not only the new JSON scope.
    _require(READINESS.get("leading_ready") is False, "A5 leading readiness unexpectedly promoted")
    _require(READINESS.get("velocity_export_ready") is False, "A5 Kokuno velocity export unexpectedly promoted")
    _require(READINESS.get("pde_validated") is False, "A5 PDE state unexpectedly promoted")
    _require(
        TRUTH_BOUNDARY.get("corrected_global_cartesian_leading_velocity_materialized") is False,
        "A5 corrected/global Cartesian leading velocity unexpectedly materialized",
    )
    _require(
        TRUTH_BOUNDARY.get("complete_candidate_api_ready") is False,
        "A5 complete candidate API unexpectedly promoted",
    )
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
        "source_coordinate_profile_executable": True,
        "profile_save_load_semantic_identity_registered": True,
        "profile_radial_derivative_audit_registered": True,
        "kokuno_cartesian_velocity_identity_materialized": False,
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
