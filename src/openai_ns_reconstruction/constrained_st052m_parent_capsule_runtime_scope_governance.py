"""Fail-closed CR002 governance for the live ST052-M replay capsule.

The integrated parent-capsule workflow is useful exact-source replay evidence, but
it is not yet a package-level ST052-M velocity loader or a whole-child save/load
capsule.  This auditor keeps those delivery stages distinct without changing any
velocity, forcing, residual, or acceptance threshold.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

CONTRACT_REL = Path("configs/st052m_parent_capsule_runtime_scope_contract.json")
CONSTRAINTS_REL = Path("configs/constraints.json")
STATUS_REL = Path("project_status.json")
CAPSULE_REL = Path("src/openai_ns_reconstruction/st052_parent_capsule.py")
WORKFLOW_REL = Path(".github/workflows/agent8-st052-parent-capsule-live.yml")


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_false(mapping: dict[str, Any], key: str, where: str) -> None:
    _require(mapping.get(key) is False, f"{where}.{key} must be exactly false")


def _require_true(mapping: dict[str, Any], key: str, where: str) -> None:
    _require(mapping.get(key) is True, f"{where}.{key} must be exactly true")


def validate_contract_data(
    contract: dict[str, Any],
    constraints: dict[str, Any],
    status: dict[str, Any],
    capsule_source: str,
    workflow_source: str,
) -> dict[str, Any]:
    _require(contract.get("schema") == "st052m-parent-capsule-runtime-scope-governance/v1", "unexpected governance schema")
    _require(contract.get("task") == "CR002", "governance task must remain CR002")

    subject = contract.get("subject")
    _require(isinstance(subject, dict), "subject missing")
    _require(subject.get("candidate_id") == "ST052-M", "subject candidate drift")
    _require(subject.get("source_pr") == 508, "source PR drift")
    _require(subject.get("source_head") == "b3b8bfdbe1077f9ec967d158602951997d81e17d", "source head drift")
    _require(subject.get("source_recipe_git_blob_sha1") == "e30c769052379f72afeee46ca264482884cc5ac7", "source recipe identity drift")

    classification = contract.get("classification")
    _require(isinstance(classification, dict), "classification missing")
    _require(set(classification) == {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"}, "source classification categories drift")
    for key, items in classification.items():
        _require(isinstance(items, list) and items, f"classification.{key} must be nonempty")

    observed = contract.get("observed_live_scope")
    _require(isinstance(observed, dict), "observed_live_scope missing")
    _require_true(observed, "source_recipe_blob_verified_inside_builder", "observed_live_scope")
    _require_false(observed, "source_head_verified_inside_builder", "observed_live_scope")
    _require_true(observed, "source_head_enforced_by_exact_checkout_workflow", "observed_live_scope")
    _require_true(observed, "candidate_file_emitted_by_workflow", "observed_live_scope")
    _require_false(observed, "candidate_file_committed_to_live_repository", "observed_live_scope")
    _require_false(observed, "package_level_parent_velocity_loader_exposed", "observed_live_scope")
    _require_false(observed, "package_level_parent_velocity_api_exposed", "observed_live_scope")
    _require_true(observed, "parent_replay_manifest_builder_exposed", "observed_live_scope")
    _require_false(observed, "whole_child_save_load_ready", "observed_live_scope")
    _require_false(observed, "experimental_st052_velocity_export_ready", "observed_live_scope")

    rules = contract.get("evidence_scope_rules")
    _require(isinstance(rules, dict), "evidence_scope_rules missing")
    _require("open PR #622" in str(rules.get("candidate_validation_linkage", "")), "#622 sibling seam must remain explicit")
    _require("does not authenticate" in str(rules.get("builder_manifest_alone", "")), "builder source-head limitation must remain explicit")

    truth = contract.get("truth_states")
    _require(isinstance(truth, dict), "truth_states missing")
    _require_true(truth, "parent_exact_source_replay_workflow_integrated", "truth_states")
    _require_true(truth, "parent_manifest_builder_integrated", "truth_states")
    for key in (
        "parent_runtime_ready",
        "whole_child_save_load_ready",
        "production_candidate_selected",
        "experimental_st052_velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require_false(truth, key, "truth_states")

    # The builder verifies the exact recipe blob, but the module itself does not
    # inspect source_root's Git commit.  Exact source-head binding is supplied by
    # the workflow checkout and therefore must not be laundered into a standalone
    # builder-level certificate.
    _require("_git_blob_sha1(recipe_path)" in capsule_source, "capsule must verify recipe Git blob")
    _require("SOURCE_HEAD = \"b3b8bfdbe1077f9ec967d158602951997d81e17d\"" in capsule_source, "capsule source-head declaration drift")
    for forbidden in ("subprocess.run", "subprocess.check_output", "git rev-parse", "rev-parse HEAD"):
        _require(forbidden not in capsule_source, f"builder source-head verification semantics changed; update governance: found {forbidden}")
    _require("def build_manifest(" in capsule_source, "parent manifest builder missing")
    _require("def velocity(" not in capsule_source, "parent capsule unexpectedly exposes velocity; governance requires review")
    _require("def load_candidate(" not in capsule_source, "parent capsule unexpectedly exposes candidate loader; governance requires review")
    _require("\"whole_child_save_load_ready\": False" in capsule_source, "capsule truth boundary must keep whole-child save/load false")

    _require("repository: Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1" in workflow_source, "exact-source checkout repository missing")
    _require("ref: b3b8bfdbe1077f9ec967d158602951997d81e17d" in workflow_source, "exact-source checkout ref drift")
    _require("candidate.json" in workflow_source and "validation.json" in workflow_source, "workflow must emit candidate and validation")
    _require("actions/upload-artifact@v4" in workflow_source, "workflow artifact upload missing")

    cr = contract.get("cr001_nonmutation")
    _require(isinstance(cr, dict), "cr001_nonmutation missing")
    _require(constraints.get("nu") == cr.get("nu") == 0.01, "nu drift")
    domain = constraints.get("domain")
    _require(isinstance(domain, dict), "constraints.domain missing")
    _require(domain.get("physical") == cr.get("physical_domain") == "R^3", "physical domain drift")
    _require(domain.get("evaluation_box") == cr.get("evaluation_box") == [[-2, 2], [-2, 2], [-2, 2]], "evaluation box drift")
    _require(domain.get("support") == cr.get("support") == "r < 2 and abs(z) < 2", "support drift")
    _require(domain.get("time_interval") == cr.get("time_interval") == [0.25, 0.75], "time interval drift")

    forcing = constraints.get("forcing")
    _require(isinstance(forcing, dict), "constraints.forcing missing")
    _require(forcing.get("mode") == cr.get("forcing_mode") == "restricted_two_parameter_family", "forcing mode drift")
    _require(forcing.get("parameters") == cr.get("forcing_parameter_bounds") == {"a": [0.0, 10.0], "c": [0.0, 10.0]}, "forcing bounds drift")

    nontriviality = constraints.get("nontriviality")
    _require(isinstance(nontriviality, dict), "constraints.nontriviality missing")
    _require(nontriviality.get("reference_energy") == cr.get("reference_energy") == 1.0, "reference energy drift")
    _require(nontriviality.get("reference_energy_abs_tolerance") == cr.get("reference_energy_abs_tolerance") == 0.001, "energy tolerance drift")

    validation = constraints.get("validation")
    _require(isinstance(validation, dict), "constraints.validation missing")
    _require(validation.get("seed") == cr.get("validation_seed") == 914027, "validation seed drift")
    _require(validation.get("held_out_points") == cr.get("held_out_points") == 4096, "held-out count drift")
    _require(validation.get("times") == cr.get("validation_times"), "validation times drift")
    _require(validation.get("derivative_steps") == cr.get("derivative_steps") == [0.02, 0.01, 0.005], "derivative ladder drift")
    _require(validation.get("quadrature_orders_per_axis") == cr.get("quadrature_orders_per_axis") == [24, 48, 96], "quadrature ladder drift")
    thresholds = validation.get("thresholds")
    _require(isinstance(thresholds, dict), "validation thresholds missing")
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _require(thresholds.get(key) == cr.get(key), f"{key} drift")
    _require_false(cr, "free_residual_defined_forcing_allowed", "cr001_nonmutation")
    _require_false(cr, "amplitude_collapse_success_allowed", "cr001_nonmutation")

    latest = status.get("latest_st052m_visual_candidate_evidence")
    _require(isinstance(latest, dict), "project_status.latest_st052m_visual_candidate_evidence missing")
    for key in ("production_candidate_selected", "velocity_export_ready", "visualization_ready", "visual_correspondence_verified", "pde_validated"):
        _require_false(latest, key, "project_status.latest_st052m_visual_candidate_evidence")
    states = status.get("states")
    _require(isinstance(states, dict), "project_status.states missing")
    _require_true(states, "velocity_export_ready", "project_status.states")
    _require(status.get("candidate_family") == "eq45_supported_velocity_candidate_v1", "canonical candidate relabelled")

    return {
        "schema": contract["schema"],
        "candidate_id": subject["candidate_id"],
        "builder_source_head_verified": False,
        "workflow_exact_source_head_bound": True,
        "parent_runtime_ready": False,
        "whole_child_save_load_ready": False,
        "experimental_st052_velocity_export_ready": False,
        "canonical_velocity_export_ready": True,
        "pde_validated": False,
    }


def audit_repository(repo_root: Path | None = None) -> dict[str, Any]:
    root = repo_root or Path(__file__).resolve().parents[2]
    return validate_contract_data(
        _read_json(root / CONTRACT_REL),
        _read_json(root / CONSTRAINTS_REL),
        _read_json(root / STATUS_REL),
        (root / CAPSULE_REL).read_text(),
        (root / WORKFLOW_REL).read_text(),
    )


if __name__ == "__main__":
    print(json.dumps(audit_repository(), indent=2, sort_keys=True))
