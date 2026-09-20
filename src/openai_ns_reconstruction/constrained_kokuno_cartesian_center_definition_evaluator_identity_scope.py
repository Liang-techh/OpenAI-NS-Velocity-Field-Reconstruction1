"""Fail-closed CR002 audit of source-field definition versus evaluator identity.

Agent-1 PR #803 is the first Kokuno PA.10 contraction-center object in this
lane with a directly callable and save/load-capable Cartesian
``velocity(x,y,z,t)->[...,3]``.  It correctly separates its field digest from
the expensive source-C evidence receipt.  The remaining representation seam is
one level finer: ``field_configuration()`` also hashes the numerical q-inverse
realization, including the autonomous bisection iteration count.

That is useful for exact executable replay, but it is not the same claim as a
digest of the public mathematical field definition.  This audit keeps those
roles separate without changing the evaluator, source profiles, C, any
scientific threshold, or any delivery/scientific state.
"""

from __future__ import annotations

import ast
import copy
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_pa10_cartesian_center_velocity import (
    DEFAULT_BISECTION_ITERATIONS,
    KokunoPA10CartesianCenterVelocity,
    SCHEMA as SOURCE_SCHEMA,
)

SCHEMA = "cr002-kokuno-center-definition-evaluator-identity-scope-v1"
TASK_ID = "CR002-KOKUNO-CENTER-DEFINITION-EVALUATOR-IDENTITY-072"
UPSTREAM_PR = 803
UPSTREAM_HEAD = "0666fff748d7c2659774b0aacba70042d045aab5"
UPSTREAM_SOURCE_BLOB = "f404eff536b91248f76095c111d788482a19f6a3"

_ROOT = Path(__file__).resolve().parents[2]
_CONTRACT_PATH = (
    _ROOT / "configs" / "kokuno_cartesian_center_definition_evaluator_identity_scope.json"
)
_CONSTRAINTS_PATH = _ROOT / "configs" / "constraints.json"
_SOURCE_PATH = (
    _ROOT
    / "src"
    / "openai_ns_reconstruction"
    / "kokuno_pa10_cartesian_center_velocity.py"
)


class GovernanceError(ValueError):
    """Raised when the frozen CR002 representation boundary drifts."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise GovernanceError(message)


def load_contract() -> dict[str, Any]:
    return json.loads(_CONTRACT_PATH.read_text())


def load_constraints() -> dict[str, Any]:
    return json.loads(_CONSTRAINTS_PATH.read_text())


def _class_method(tree: ast.AST, class_name: str, method_name: str) -> ast.FunctionDef:
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and item.name == method_name:
                    return item
    raise GovernanceError(f"missing {class_name}.{method_name}")


def _module_assignment(tree: ast.Module, name: str) -> Any:
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(isinstance(target, ast.Name) and target.id == name for target in node.targets):
            return ast.literal_eval(node.value)
    raise GovernanceError(f"missing module assignment {name}")


def validate_source_structure(source_text: str) -> dict[str, Any]:
    """Bind the distinction to the exact #803 evaluator representation."""
    tree = ast.parse(source_text)
    default_iterations = _module_assignment(tree, "DEFAULT_BISECTION_ITERATIONS")
    _require(default_iterations == 80, "default bisection iteration count drift")

    solve_q = _class_method(tree, "KokunoPA10CartesianCenterVelocity", "solve_q")
    solve_q_source = ast.get_source_segment(source_text, solve_q) or ""
    for token in (
        "q_star = np.power(np.abs(z_arr), 1.0 / self.D)",
        "high = q_star + tau / derivative_floor",
        "for _ in range(self.bisection_iterations):",
        "mid = low + 0.5 * (high - low)",
    ):
        _require(token in solve_q_source, f"q-inverse realization missing {token}")

    field_configuration = _class_method(
        tree, "KokunoPA10CartesianCenterVelocity", "field_configuration"
    )
    field_configuration_source = ast.get_source_segment(
        source_text, field_configuration
    ) or ""
    for token in (
        '"physical_profile_configuration": self.physical_profiles.configuration()',
        '"time_interval": [float(v) for v in self.time_interval]',
        '"q_inverse": {',
        '"method": "source-monotone-certified-bracket-bisection"',
        '"iterations": int(self.bisection_iterations)',
    ):
        _require(
            token in field_configuration_source,
            f"field_configuration identity binding missing {token}",
        )

    field_sha = _class_method(
        tree, "KokunoPA10CartesianCenterVelocity", "field_sha256"
    )
    field_sha_source = ast.get_source_segment(source_text, field_sha) or ""
    _require(
        "_canonical_json(self.field_configuration())" in field_sha_source,
        "field_sha256 no longer hashes complete field_configuration",
    )

    from_configuration = _class_method(
        tree, "KokunoPA10CartesianCenterVelocity", "from_configuration"
    )
    from_configuration_source = ast.get_source_segment(
        source_text, from_configuration
    ) or ""
    _require(
        'q_inverse.get("method") != "source-monotone-certified-bracket-bisection"'
        in from_configuration_source,
        "save/load no longer binds q-inverse method",
    )
    _require(
        'bisection_iterations=int(q_inverse["iterations"])' in from_configuration_source,
        "save/load no longer binds q-inverse iterations",
    )

    return {
        "default_bisection_iterations": int(default_iterations),
        "q_inverse_uses_configured_iteration_count": True,
        "field_configuration_binds_physical_profile_configuration": True,
        "field_configuration_binds_time_interval": True,
        "field_configuration_binds_q_inverse_method": True,
        "field_configuration_binds_q_inverse_iterations": True,
        "field_sha256_hashes_complete_field_configuration": True,
        "save_load_preserves_evaluator_realization": True,
    }


def validate_contract(contract: Mapping[str, Any]) -> None:
    _require(contract.get("schema_version") == 1, "contract schema_version drift")
    _require(contract.get("task_id") == TASK_ID, "task id drift")

    upstream = contract.get("upstream_binding", {})
    _require(upstream.get("pr") == UPSTREAM_PR, "upstream PR drift")
    _require(upstream.get("head") == UPSTREAM_HEAD, "upstream head drift")
    _require(
        upstream.get("source_blob_sha") == UPSTREAM_SOURCE_BLOB,
        "upstream source blob drift",
    )
    _require(upstream.get("schema") == SOURCE_SCHEMA, "upstream source schema drift")

    classes = contract.get("source_classification", {})
    _require(
        set(classes)
        == {
            "user_requirement",
            "public_source_fact",
            "autonomous_design",
            "pending_unknown",
        },
        "source classes must remain exactly the canonical four",
    )
    _require(all(classes[key] for key in classes), "every source class must be nonempty")
    public = "\n".join(classes["public_source_fact"]).lower()
    autonomous = "\n".join(classes["autonomous_design"]).lower()
    _require(
        "bisection_iterations" not in public
        and "default_bisection_iterations" not in public,
        "autonomous evaluator iteration count laundered into public-source facts",
    )
    _require(
        "default_bisection_iterations=80" in autonomous,
        "autonomous bisection iteration choice must stay explicit",
    )

    facts = contract.get("current_representation_facts", {})
    for key in (
        "inner_velocity_callable",
        "velocity_vectorized",
        "field_configuration_save_load_roundtrip",
        "inner_source_X_only_fail_closed",
        "field_configuration_contains_physical_profile_configuration",
        "field_configuration_contains_registered_time_interval",
        "field_configuration_contains_q_inverse_method",
        "field_configuration_contains_bisection_iterations",
        "field_sha256_hashes_complete_field_configuration",
        "bisection_iterations_change_field_sha256",
    ):
        _require(facts.get(key) is True, f"representation fact lost: {key}")
    for key in (
        "changing_only_bisection_iterations_changes_public_source_formula",
        "changing_only_bisection_iterations_changes_physical_profile_configuration",
    ):
        _require(facts.get(key) is False, f"representation scope laundering: {key}")

    roles = contract.get("identity_roles", {})
    _require(
        roles.get("current_field_sha256_role")
        == "executable_evaluator_realization_identity",
        "current field_sha256 role drift",
    )
    for key in (
        "current_field_sha_is_valid_for_exact_numeric_replay_binding",
        "bisection_iterations_are_autonomous_evaluator_parameters",
    ):
        _require(roles.get(key) is True, f"missing identity rule: {key}")
    for key in (
        "current_field_sha_is_source_mathematical_definition_identity",
        "dedicated_source_mathematical_definition_digest_exposed_by_803",
        "bisection_iterations_are_public_source_parameters",
        "field_sha_change_alone_proves_public_source_field_definition_changed",
        "mathematical_definition_equivalence_implies_bitwise_velocity_equality",
    ):
        _require(roles.get(key) is False, f"identity laundering: {key}")

    cross = contract.get("cross_implementation_validation_rules", {})
    for key in (
        "independent_validator_may_use_different_q_solver_realization",
        "source_field_correspondence_must_bind_source_formula_profile_configuration_and_domain",
        "exact_numeric_replay_must_bind_evaluator_realization_identity",
        "different_evaluator_sha_may_still_target_same_source_mathematical_field",
    ):
        _require(cross.get(key) is True, f"cross-implementation rule lost: {key}")
    for key in (
        "different_evaluator_sha_may_be_called_same_exact_executable_artifact",
        "numeric_equality_across_realizations_is_assumed_without_test",
    ):
        _require(cross.get(key) is False, f"cross-implementation laundering: {key}")

    delivery = contract.get("delivery_claim_boundary", {})
    _require(delivery.get("inner_velocity_callable") is True, "inner velocity lost")
    _require(
        delivery.get("inner_velocity_save_load_capable") is True,
        "save/load capability lost",
    )
    for key in (
        "unified_kokuno_velocity_export_ready",
        "global_leading_velocity_materialized",
        "complete_candidate_api_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(delivery.get(key) is False, f"premature delivery/scientific promotion: {key}")

    claims = contract.get("claim_boundaries", {})
    _require(claims, "missing claim boundaries")
    for key, value in claims.items():
        _require(value is False, f"claim boundary promoted: {key}")


def validate_cr001(constraints: Mapping[str, Any], contract: Mapping[str, Any]) -> None:
    lock = contract["CR001_lock"]
    _require(constraints["nu"] == lock["nu"] == 0.01, "nu drift")
    domain = constraints["domain"]
    _require(
        domain["physical"] == lock["physical_domain"] == "R^3",
        "physical domain drift",
    )
    _require(domain["evaluation_box"] == lock["evaluation_box"], "evaluation box drift")
    _require(domain["support"] == lock["support"], "support drift")
    _require(domain["time_interval"] == lock["time_interval"], "time interval drift")

    forcing = constraints["forcing"]
    _require(forcing["mode"] == lock["forcing_mode"], "forcing mode drift")
    _require(
        forcing["parameters"] == lock["forcing_parameter_bounds"],
        "forcing parameter bounds drift",
    )
    restriction = forcing["restriction"].lower()
    _require("no residual-dependent" in restriction, "residual-defined force allowed")
    _require("pointwise free force" in restriction, "pointwise free force allowed")

    nontriviality = constraints["nontriviality"]
    _require(
        nontriviality["reference_energy"] == lock["reference_energy"],
        "reference energy drift",
    )
    _require(
        nontriviality["reference_energy_abs_tolerance"]
        == lock["reference_energy_abs_tolerance"],
        "reference-energy tolerance drift",
    )
    _require(
        "reject collapsed candidates" in nontriviality["enforcement"].lower(),
        "amplitude-collapse rejection drift",
    )

    validation = constraints["validation"]
    _require(validation["seed"] == lock["validation_seed"], "validation seed drift")
    _require(
        validation["held_out_points"] == lock["held_out_points"],
        "held-out point count drift",
    )
    _require(validation["times"] == lock["validation_times"], "validation times drift")
    _require(
        validation["derivative_steps"] == lock["derivative_steps"],
        "derivative ladder drift",
    )
    _require(
        validation["quadrature_orders_per_axis"]
        == lock["quadrature_orders_per_axis"],
        "quadrature ladder drift",
    )
    thresholds = validation["thresholds"]
    for key in (
        "divergence_max",
        "divergence_L2",
        "pde_residual_max",
        "pde_residual_L2",
    ):
        _require(thresholds[key] == lock[key], f"threshold drift: {key}")
    _require(
        "changing thresholds requires a new experiment version"
        in validation["failure_policy"].lower(),
        "post-hoc threshold relaxation guard drift",
    )


def evaluator_identity_witness() -> dict[str, Any]:
    """Show that the current digest intentionally binds evaluator realization.

    No equality of velocity values across realization settings is assumed here.
    """
    left = KokunoPA10CartesianCenterVelocity(
        bisection_iterations=DEFAULT_BISECTION_ITERATIONS
    )
    right = KokunoPA10CartesianCenterVelocity(
        source_center=left.source_center,
        time_interval=left.time_interval,
        bisection_iterations=96,
    )
    left_config = left.field_configuration()
    right_config = right.field_configuration()

    _require(
        left_config["physical_profile_configuration"]
        == right_config["physical_profile_configuration"],
        "physical profile changed while varying only evaluator iterations",
    )
    _require(
        left_config["time_interval"] == right_config["time_interval"],
        "registered time interval changed while varying only evaluator iterations",
    )
    _require(
        left_config["q_inverse"]["method"] == right_config["q_inverse"]["method"],
        "q-inverse method changed unexpectedly",
    )
    _require(
        left_config["q_inverse"]["iterations"] == 80
        and right_config["q_inverse"]["iterations"] == 96,
        "witness iteration counts drifted",
    )
    _require(left.field_sha256 != right.field_sha256, "field SHA must bind iterations")

    return {
        "same_physical_profile_configuration": True,
        "same_registered_time_interval": True,
        "same_source_q_inverse_equation": True,
        "same_q_inverse_method": True,
        "different_evaluator_iteration_count": True,
        "different_current_field_sha256": True,
        "velocity_bitwise_equality_claimed": False,
        "public_source_definition_difference_claimed": False,
        "left_iterations": 80,
        "right_iterations": 96,
    }


def audit(
    *,
    contract: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    source_text: str | None = None,
) -> dict[str, Any]:
    contract_payload = (
        copy.deepcopy(dict(contract)) if contract is not None else load_contract()
    )
    constraints_payload = (
        copy.deepcopy(dict(constraints)) if constraints is not None else load_constraints()
    )
    source_payload = source_text if source_text is not None else _SOURCE_PATH.read_text()

    validate_contract(contract_payload)
    validate_cr001(constraints_payload, contract_payload)
    source_structure = validate_source_structure(source_payload)
    witness = evaluator_identity_witness()

    _require(DEFAULT_BISECTION_ITERATIONS == 80, "runtime iteration default drift")
    _require(
        SOURCE_SCHEMA == contract_payload["upstream_binding"]["schema"],
        "runtime source schema drift",
    )

    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "upstream_pr": UPSTREAM_PR,
        "upstream_head": UPSTREAM_HEAD,
        "source_structure": source_structure,
        "evaluator_identity_witness": witness,
        "identity_state": {
            "inner_velocity_callable": True,
            "current_field_sha_is_exact_evaluator_replay_identity": True,
            "source_mathematical_definition_digest_exposed": False,
            "solver_iterations_are_public_source_data": False,
            "unified_kokuno_velocity_export_ready_promoted": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
        "CR001_unchanged": True,
    }


def _main() -> None:
    print(json.dumps(audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    _main()
