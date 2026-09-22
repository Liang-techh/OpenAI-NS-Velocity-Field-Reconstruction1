"""Fail-closed CR002 audit for the #1124 pulse-entry residual role.

This module does not evaluate or modify a velocity field. It only guards the
semantic boundary between an already-frozen construction closure residual used
for current-candidate M/J bookkeeping and the preregistered Navier--Stokes
residual/forcing contract.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

EXPECTED_SCHEMA = "cr002-kokuno-current-pulse-entry-residual-role-scope-v1"
EXPECTED_TASK = "CR002-KOKUNO-PULSE-ENTRY-RESIDUAL-ROLE-138"
EXPECTED_PARENT_HEAD = "4be97cbeb3b78b4d6a162447ae8fea73d0c2c477"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload


def _git_blob_sha(path: Path) -> str:
    data = path.read_bytes()
    header = f"blob {len(data)}\0".encode()
    return hashlib.sha1(header + data).hexdigest()


def _eq(actual: Any, expected: Any, label: str, errors: list[str]) -> None:
    if actual != expected:
        errors.append(f"{label}: expected {expected!r}, got {actual!r}")


def _false(mapping: Mapping[str, Any], key: str, errors: list[str]) -> None:
    if mapping.get(key) is not False:
        errors.append(f"{key} must remain false")


def _true(mapping: Mapping[str, Any], key: str, errors: list[str]) -> None:
    if mapping.get(key) is not True:
        errors.append(f"{key} must remain true")


def bookkeeping_projection_witness(
    residual_a: Sequence[float], residual_b: Sequence[float]
) -> dict[str, Any]:
    """Show that identical M/J bookkeeping does not identify a five-row residual.

    This is an autonomous mechanics witness only. It does not model the real
    Kokuno construction residual and it is not Navier--Stokes evidence.
    """

    a = tuple(float(v) for v in residual_a)
    b = tuple(float(v) for v in residual_b)
    if len(a) != 5 or len(b) != 5:
        raise ValueError("mechanics witness residuals must each have five rows")
    if not all(math.isfinite(v) for v in a + b):
        raise ValueError("mechanics witness residuals must be finite")
    return {
        "same_mj_projection": a[:2] == b[:2],
        "same_full_residual": a == b,
        "norm_a": math.sqrt(sum(v * v for v in a)),
        "norm_b": math.sqrt(sum(v * v for v in b)),
    }


def audit_contract(
    contract_path: Path,
    constraints_path: Path,
    candidate_module_path: Path,
    project_status_path: Path,
) -> list[str]:
    errors: list[str] = []
    contract = _load_json(contract_path)

    _eq(contract.get("schema"), EXPECTED_SCHEMA, "schema", errors)
    _eq(contract.get("task_id"), EXPECTED_TASK, "task_id", errors)
    exact_parent = contract.get("exact_parent", {})
    if not isinstance(exact_parent, dict):
        errors.append("exact_parent must be an object")
        exact_parent = {}
    _eq(exact_parent.get("pr"), 1124, "exact_parent.pr", errors)
    _eq(exact_parent.get("head_sha"), EXPECTED_PARENT_HEAD, "exact_parent.head_sha", errors)
    _eq(
        exact_parent.get("candidate_module"),
        "src/openai_ns_reconstruction/kokuno_current_pulse_entry_moments.py",
        "exact_parent.candidate_module",
        errors,
    )

    provenance = contract.get("four_way_provenance", {})
    if not isinstance(provenance, dict):
        errors.append("four_way_provenance must be an object")
        provenance = {}
    expected_classes = {
        "user_requirements",
        "public_source_facts",
        "autonomous_design",
        "pending_unknown",
    }
    _eq(set(provenance), expected_classes, "four_way_provenance classes", errors)
    for key in expected_classes:
        value = provenance.get(key)
        if not isinstance(value, list) or not value:
            errors.append(f"four_way_provenance.{key} must be a nonempty list")

    role = contract.get("role_separation", {})
    if not isinstance(role, dict):
        errors.append("role_separation must be an object")
        role = {}
    _eq(
        role.get("construction_closure_residual_classification"),
        "autonomous_current_candidate_bookkeeping",
        "construction closure residual classification",
        errors,
    )
    for key in (
        "construction_closure_residual_is_public_source_data",
        "construction_closure_residual_is_forcing",
        "construction_closure_residual_is_heldout_ns_residual",
        "construction_closure_residual_may_authorize_force_fitting",
        "residual_defined_free_forcing_allowed",
        "current_cartesian_end_compensation_composed",
        "unified_global_kokuno_velocity_export_ready",
    ):
        _false(role, key, errors)
    for key in (
        "current_lineage_M_entry_materialized",
        "current_lineage_J_entry_materialized",
        "public_end_compensator_inputs_ready",
    ):
        _true(role, key, errors)

    firewall = contract.get("evidence_transfer_firewall", {})
    if not isinstance(firewall, dict):
        errors.append("evidence_transfer_firewall must be an object")
        firewall = {}
    forbidden = firewall.get("forbidden_implications")
    if not isinstance(forbidden, list) or len(forbidden) < 7:
        errors.append("evidence_transfer_firewall.forbidden_implications is incomplete")
    else:
        required_fragments = (
            "body-force",
            "held-out complete NS residual",
            "Cartesian velocity composition",
            "global velocity export readiness",
            "visual correspondence",
            "PDE validation",
            "paper exactness or OpenAI-field identity",
        )
        joined = "\n".join(str(v) for v in forbidden)
        for fragment in required_fragments:
            if fragment not in joined:
                errors.append(f"missing forbidden implication fragment: {fragment}")

    cr001 = contract.get("cr001_invariants", {})
    if not isinstance(cr001, dict):
        errors.append("cr001_invariants must be an object")
        cr001 = {}
    _eq(cr001.get("constraints_blob_sha"), EXPECTED_CONSTRAINTS_BLOB, "constraints blob", errors)
    _eq(cr001.get("nu"), 0.01, "nu", errors)
    _eq(cr001.get("physical_domain"), "R^3", "physical domain", errors)
    _eq(cr001.get("evaluation_box"), [[-2, 2], [-2, 2], [-2, 2]], "evaluation box", errors)
    _eq(cr001.get("support"), "r < 2 and abs(z) < 2", "support", errors)
    _eq(cr001.get("time_interval"), [0.25, 0.75], "time interval", errors)
    _eq(cr001.get("forcing_mode"), "restricted_two_parameter_family", "forcing mode", errors)
    _eq(cr001.get("validation_seed"), 914027, "validation seed", errors)
    _eq(cr001.get("held_out_points"), 4096, "held-out points", errors)
    _eq(cr001.get("derivative_steps"), [0.02, 0.01, 0.005], "derivative steps", errors)
    _eq(cr001.get("quadrature_orders_per_axis"), [24, 48, 96], "quadrature orders", errors)
    _false(cr001, "threshold_relaxation_this_increment", errors)
    _false(cr001, "forcing_family_change_this_increment", errors)
    _eq(
        cr001.get("thresholds"),
        {
            "divergence_max": 1e-5,
            "divergence_L2": 1e-5,
            "pde_residual_max": 1e-3,
            "pde_residual_L2": 1e-3,
        },
        "CR001 thresholds",
        errors,
    )

    if _git_blob_sha(constraints_path) != EXPECTED_CONSTRAINTS_BLOB:
        errors.append("canonical configs/constraints.json blob changed")
    constraints = _load_json(constraints_path)
    _eq(constraints.get("nu"), 0.01, "live constraints.nu", errors)
    domain = constraints.get("domain", {})
    _eq(domain.get("physical"), "R^3", "live physical domain", errors)
    _eq(domain.get("evaluation_box"), [[-2, 2], [-2, 2], [-2, 2]], "live evaluation box", errors)
    _eq(domain.get("support"), "r < 2 and abs(z) < 2", "live support", errors)
    _eq(domain.get("time_interval"), [0.25, 0.75], "live time interval", errors)
    forcing = constraints.get("forcing", {})
    _eq(forcing.get("mode"), "restricted_two_parameter_family", "live forcing mode", errors)
    restriction = str(forcing.get("restriction", ""))
    if "No residual-dependent basis or pointwise free force" not in restriction:
        errors.append("live forcing restriction no longer forbids residual-dependent free force")
    nontriviality = constraints.get("nontriviality", {})
    _eq(nontriviality.get("reference_energy"), 1.0, "live reference energy", errors)
    _eq(
        nontriviality.get("reference_energy_abs_tolerance"),
        0.001,
        "live reference energy tolerance",
        errors,
    )
    validation = constraints.get("validation", {})
    _eq(validation.get("seed"), 914027, "live validation seed", errors)
    _eq(validation.get("held_out_points"), 4096, "live held-out points", errors)
    _eq(validation.get("derivative_steps"), [0.02, 0.01, 0.005], "live derivative steps", errors)
    _eq(
        validation.get("quadrature_orders_per_axis"),
        [24, 48, 96],
        "live quadrature orders",
        errors,
    )
    _eq(validation.get("thresholds", {}).get("pde_residual_max"), 0.001, "live momentum max gate", errors)
    _eq(validation.get("thresholds", {}).get("pde_residual_L2"), 0.001, "live momentum L2 gate", errors)
    _eq(validation.get("thresholds", {}).get("divergence_max"), 1e-5, "live divergence max gate", errors)
    _eq(validation.get("thresholds", {}).get("divergence_L2"), 1e-5, "live divergence L2 gate", errors)

    module_text = candidate_module_path.read_text()
    required_module_fragments = (
        '"current_lineage_J_entry_materialized": True',
        '"current_J_from_frozen_I1_closure_residual": True',
        '"current_J_assumed_zero": False',
        '"current_cartesian_end_compensation_composed": False',
        '"restricted_forcing_materialized": False',
        '"heldout_ns_residual_assessed": False',
        '"pde_validated": False',
        "This is a current-candidate bookkeeping seam.",
    )
    for fragment in required_module_fragments:
        if fragment not in module_text:
            errors.append(f"#1124 module role/truth fragment missing: {fragment}")

    states = contract.get("independent_delivery_states", {})
    if not isinstance(states, dict):
        errors.append("independent_delivery_states must be an object")
        states = {}
    _true(states, "canonical_eq45_velocity_export_ready", errors)
    for key in (
        "kokuno_pulse_entry_route_velocity_export_ready",
        "kokuno_visual_correspondence_verified",
        "kokuno_pde_validated",
        "kokuno_paper_exact",
        "kokuno_openai_field_identified",
    ):
        _false(states, key, errors)

    project_status = _load_json(project_status_path)
    live_states = project_status.get("states", {})
    _true(live_states, "velocity_export_ready", errors)
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        _false(live_states, key, errors)

    witness = contract.get("mechanics_witness", {})
    if not isinstance(witness, dict):
        errors.append("mechanics_witness must be an object")
        witness = {}
    _eq(witness.get("classification"), "autonomous_mechanics_only", "witness classification", errors)
    try:
        receipt = bookkeeping_projection_witness(witness.get("residual_a", []), witness.get("residual_b", []))
    except (TypeError, ValueError) as exc:
        errors.append(f"invalid mechanics witness: {exc}")
    else:
        if receipt["same_mj_projection"] is not True:
            errors.append("mechanics witness must preserve the M/J projection")
        if receipt["same_full_residual"] is not False:
            errors.append("mechanics witness must differ outside the M/J projection")
        if math.isclose(receipt["norm_a"], receipt["norm_b"], rel_tol=0.0, abs_tol=0.0):
            errors.append("mechanics witness full residual norms must differ")

    return errors


def assert_contract(
    contract_path: Path,
    constraints_path: Path,
    candidate_module_path: Path,
    project_status_path: Path,
) -> None:
    errors = audit_contract(contract_path, constraints_path, candidate_module_path, project_status_path)
    if errors:
        raise ValueError("CR002 pulse-entry residual-role audit failed:\n- " + "\n- ".join(errors))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, default=Path("configs/kokuno_current_pulse_entry_residual_role_scope.json"))
    parser.add_argument("--constraints", type=Path, default=Path("configs/constraints.json"))
    parser.add_argument(
        "--candidate-module",
        type=Path,
        default=Path("src/openai_ns_reconstruction/kokuno_current_pulse_entry_moments.py"),
    )
    parser.add_argument("--project-status", type=Path, default=Path("project_status.json"))
    args = parser.parse_args()
    errors = audit_contract(args.contract, args.constraints, args.candidate_module, args.project_status)
    if errors:
        for error in errors:
            print(f"FAIL: {error}")
        return 1
    print("PASS: CR002 pulse-entry construction-residual role remains fail-closed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
