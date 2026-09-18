"""Fail-closed CR002 governance for render-readiness and visual-evidence scope.

An immutable fixed-protocol render is useful delivery evidence, but it does not
by itself satisfy repository-wide visualization readiness, public visual
correspondence, material-path evidence, or PDE acceptance.  This module changes
no velocity, candidate, route, rendering implementation, or scientific gate.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-RENDER-READINESS-SCOPE-082"
_INTEGRATION_HEAD = "09682c9a4c70148113639220387344ea6c9b9538"
_INTEGRATED_PR = 490
_SOURCE_RENDER_PR = 487
_SOURCE_RENDER_RUN = 35395433367
_SOURCE_RENDER_ARTIFACT_ID = 10568105023
_SOURCE_RENDER_ARTIFACT_DIGEST = (
    "sha256:55b69f275f2c160dfecc6e1fc25436dee5860e8e6ea35c0a0c10bec1e9e89d3f"
)
_INTEGRATION_TEST_RUN = 35396522245
_CANONICAL_SHA = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
_CANONICAL_API = "openai_ns_reconstruction.eq45_supported_delivery:velocity"
_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected object in {path}")
    return data


def load_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "render_readiness_scope_contract.json")


def load_delivery_state_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "delivery_state_contract.json")


def load_delivery_readiness_scope_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "delivery_readiness_scope_contract.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def audit_render_readiness_scope(
    contract: Mapping[str, Any] | None = None,
    delivery_state_contract: Mapping[str, Any] | None = None,
    delivery_readiness_scope_contract: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit render-evidence semantics against live delivery and CR001 contracts."""

    contract = deepcopy(dict(load_contract() if contract is None else contract))
    delivery_state = deepcopy(
        dict(load_delivery_state_contract() if delivery_state_contract is None else delivery_state_contract)
    )
    delivery_scope = deepcopy(
        dict(
            load_delivery_readiness_scope_contract()
            if delivery_readiness_scope_contract is None
            else delivery_readiness_scope_contract
        )
    )
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))

    _require_equal(contract.get("schema_version"), 1, "contract schema")
    _require_equal(contract.get("task_id"), TASK_ID, "contract task id")

    snapshot = contract.get("snapshot", {})
    _require_equal(snapshot.get("active_integration_branch"), "codex/cr001-constraints", "integration branch")
    _require_equal(snapshot.get("active_integration_head"), _INTEGRATION_HEAD, "integration head")
    _require_equal(snapshot.get("integrated_render_pr"), _INTEGRATED_PR, "integrated render PR")
    _require_equal(snapshot.get("integrated_render_merge_commit"), _INTEGRATION_HEAD, "render merge commit")
    _require_equal(snapshot.get("source_render_pr"), _SOURCE_RENDER_PR, "source render PR")
    _require_equal(snapshot.get("source_render_workflow_run"), _SOURCE_RENDER_RUN, "source render run")
    _require_equal(snapshot.get("source_render_artifact_id"), _SOURCE_RENDER_ARTIFACT_ID, "render artifact id")
    _require_equal(
        snapshot.get("source_render_artifact_digest"),
        _SOURCE_RENDER_ARTIFACT_DIGEST,
        "render artifact digest",
    )
    _require_equal(snapshot.get("integration_tests_run"), _INTEGRATION_TEST_RUN, "integration test run")
    _require_equal(snapshot.get("integration_tests_conclusion"), "success", "integration test conclusion")

    _require_equal(set(contract.get("allowed_source_classes", ())), _ALLOWED_CLASSES, "source-class vocabulary")
    source = contract.get("source_classification", {})
    _require(source and set(source.values()) <= _ALLOWED_CLASSES, "invalid source classification")
    for key in ("callable_velocity_x_y_z_t", "reproducible_visual_comparison"):
        _require_equal(source.get(key), "user_requirement", f"source class {key}")
    _require_equal(
        source.get("public_openai_visualization_observables"),
        "public_source_fact",
        "public visualization source class",
    )
    for key in (
        "render_readiness_state_names_and_scope_rules",
        "fixed_seed_camera_times_resolution_quantile_streamline_protocol",
        "sampled_grid_interpolation_and_vorticity_diagnostic",
        "cr001_validation_protocol",
    ):
        _require_equal(source.get(key), "autonomous_design", f"source class {key}")
    for key in (
        "openai_hidden_numerical_velocity",
        "openai_hidden_parameters",
        "openai_hidden_frame_time_camera_seed_mapping",
    ):
        _require_equal(source.get(key), "pending_unknown", f"source class {key}")

    readiness = contract.get("readiness_semantics", {})
    artifact = readiness.get("render_artifact_ready", {})
    visualization = readiness.get("unqualified_visualization_ready", {})
    correspondence = readiness.get("visual_correspondence_verified", {})
    repository_visualization = delivery_state.get("states", {}).get("visualization_ready", {})
    repository_correspondence = delivery_state.get("states", {}).get("visual_correspondence_verified", {})

    _require_equal(
        visualization.get("required_evidence"),
        repository_visualization.get("positive_evidence"),
        "visualization-ready evidence",
    )
    _require_equal(
        correspondence.get("required_evidence"),
        repository_correspondence.get("positive_evidence"),
        "visual-correspondence evidence",
    )
    for key in (
        "does_not_imply_visualization_ready",
        "does_not_imply_visual_correspondence_verified",
        "does_not_imply_pde_validated",
        "does_not_imply_paper_exact",
    ):
        _require_equal(artifact.get(key), True, f"render artifact scope {key}")
    for key in (
        "single_sampled_grid_render_alone_is_sufficient",
        "target_free_render_alone_is_sufficient",
        "green_render_ci_alone_is_sufficient",
    ):
        _require_equal(visualization.get(key), False, f"visualization scope {key}")
    for key in (
        "target_free_render_is_sufficient",
        "single_resolution_render_is_sufficient",
        "render_artifact_ready_is_sufficient",
    ):
        _require_equal(correspondence.get(key), False, f"correspondence scope {key}")

    previous_sampled = delivery_scope.get("readiness_semantics", {}).get("sampled_grid_export_ready", {})
    _require_equal(
        previous_sampled.get("does_not_imply_visualization_ready"),
        True,
        "sampled-grid/visualization separation",
    )
    _require_equal(
        previous_sampled.get("does_not_imply_visual_correspondence_verified"),
        True,
        "sampled-grid/correspondence separation",
    )

    observed = contract.get("integrated_pr490", {})
    _require_equal(observed.get("candidate_role"), "noncanonical_frozen_st051b_redistribution_child", "PR490 role")
    _require_equal(observed.get("redistribution_gain"), 0.025, "PR490 gain")
    _require_equal(observed.get("grid_resolution"), 33, "PR490 grid resolution")
    _require_equal(observed.get("times"), [0.25, 0.5, 0.75], "PR490 times")
    _require_equal(observed.get("streamline_count"), 48, "PR490 streamline count")
    _require_equal(
        observed.get("streamline_kind"),
        "instantaneous_integral_curve_of_sampled_velocity_direction",
        "PR490 streamline semantics",
    )
    _require_equal(observed.get("interpolation"), "trilinear_on_sampled_grid", "PR490 interpolation")
    _require_equal(observed.get("integrator"), "project_local_vectorized_rk4", "PR490 integrator")
    _require_equal(
        observed.get("vorticity_operator"),
        "second_order_cartesian_finite_difference_on_same_sampled_grid",
        "PR490 vorticity operator",
    )
    for key in (
        "streamline_time_evolution_claim",
        "streamline_may_be_claimed_as_material_path",
        "sampled_vorticity_may_be_claimed_as_pde_acceptance_evidence",
        "sampled_vorticity_may_be_claimed_as_continuum_certificate",
        "render_resolution_stability_assessed",
        "public_openai_frame_comparison_performed",
        "openai_numeric_target_used",
        "camera_or_seed_fitted_to_openai",
        "candidate_registered_as_canonical",
        "unqualified_visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        _require_equal(observed.get(key), False, f"PR490 truth {key}")
    _require_equal(observed.get("render_artifact_ready"), True, "PR490 render artifact readiness")

    forbidden = contract.get("forbidden_inferences", [])
    _require_equal(len(forbidden), 5, "forbidden inference count")
    by_target = {row.get("to"): row for row in forbidden}
    _require_equal(
        by_target["visualization_ready"].get("allowed_without_additional_evidence"),
        False,
        "render artifact cannot promote visualization readiness",
    )
    _require_equal(
        by_target["visual_correspondence_verified"].get("allowed_without_public_reference_comparison_and_resolution_stability"),
        False,
        "target-free render cannot promote visual correspondence",
    )
    _require_equal(by_target["material_path_or_time_evolution_evidence"].get("allowed"), False, "streamline/pathline separation")
    _require_equal(by_target["pde_or_divergence_acceptance_evidence"].get("allowed"), False, "vorticity/PDE separation")
    _require_equal(by_target["paper_exact_or_openai_field_identified"].get("allowed"), False, "render/exactness separation")

    canonical = contract.get("canonical_delivery", {})
    _require_equal(canonical.get("candidate_family"), project_status.get("candidate_family"), "canonical family")
    _require_equal(canonical.get("candidate_sha256"), project_status.get("candidate_sha256"), "canonical SHA")
    _require_equal(canonical.get("candidate_sha256"), _CANONICAL_SHA, "known canonical SHA")
    _require_equal(canonical.get("velocity_api"), project_status.get("velocity_api"), "canonical API")
    _require_equal(canonical.get("velocity_api"), _CANONICAL_API, "known canonical API")
    states = project_status.get("states", {})
    for key, expected in (
        ("velocity_export_ready", True),
        ("visualization_ready", False),
        ("visual_correspondence_verified", False),
        ("pde_validated", False),
        ("paper_exact", False),
        ("openai_field_identified", False),
    ):
        _require_equal(canonical.get(key), expected, f"contract canonical {key}")
        _require_equal(states.get(key), expected, f"project status {key}")
    _require_equal(canonical.get("integrated_pr490_changes_canonical_candidate"), False, "PR490 canonical candidate separation")
    _require_equal(canonical.get("integrated_pr490_changes_canonical_velocity_api"), False, "PR490 canonical API separation")

    domain = constraints.get("domain", {})
    forcing = constraints.get("forcing", {})
    nontriviality = constraints.get("nontriviality", {})
    validation = constraints.get("validation", {})
    thresholds = validation.get("thresholds", {})
    cr001 = contract.get("canonical_cr001", {})
    expected_pairs = (
        (cr001.get("nu"), constraints.get("nu"), "nu"),
        (cr001.get("physical_domain"), domain.get("physical"), "physical domain"),
        (cr001.get("evaluation_box"), domain.get("evaluation_box"), "evaluation box"),
        (cr001.get("support"), domain.get("support"), "support"),
        (cr001.get("time_interval"), domain.get("time_interval"), "time interval"),
        (cr001.get("force_mode"), forcing.get("mode"), "force mode"),
        (cr001.get("force_bounds"), forcing.get("parameters"), "force bounds"),
        (cr001.get("reference_energy"), nontriviality.get("reference_energy"), "reference energy"),
        (cr001.get("reference_energy_abs_tolerance"), nontriviality.get("reference_energy_abs_tolerance"), "energy tolerance"),
        (cr001.get("validation_seed"), validation.get("seed"), "validation seed"),
        (cr001.get("held_out_points"), validation.get("held_out_points"), "held-out points"),
        (cr001.get("validation_times"), validation.get("times"), "validation times"),
        (cr001.get("derivative_steps"), validation.get("derivative_steps"), "derivative ladder"),
        (cr001.get("divergence_max"), thresholds.get("divergence_max"), "divergence max"),
        (cr001.get("divergence_L2"), thresholds.get("divergence_L2"), "divergence L2"),
        (cr001.get("pde_residual_max"), thresholds.get("pde_residual_max"), "PDE max"),
        (cr001.get("pde_residual_L2"), thresholds.get("pde_residual_L2"), "PDE L2"),
    )
    for expected, actual, label in expected_pairs:
        _require_equal(actual, expected, label)
    _require("No residual-dependent basis" in forcing.get("restriction", ""), "free-force restriction weakened")

    truth = contract.get("truth_boundary", {})
    _require_equal(truth.get("canonical_velocity_export_ready"), True, "canonical velocity readiness")
    _require_equal(truth.get("integrated_pr490_render_artifact_ready"), True, "integrated render artifact readiness")
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require_equal(truth.get(key), False, f"truth boundary {key}")
        _require_equal(states.get(key), False, f"project truth {key}")

    mutation = contract.get("mutation_scope", {})
    for key in (
        "velocity_changed",
        "candidate_changed",
        "pressure_changed",
        "forcing_changed",
        "optimizer_changed",
        "validation_sample_changed",
        "norm_definition_changed",
        "threshold_changed",
        "active_route_changed",
        "canonical_velocity_api_changed",
        "render_implementation_changed",
    ):
        _require_equal(mutation.get(key), False, key)

    return {
        "task_id": TASK_ID,
        "status": "governance_pass",
        "integrated_render_pr": _INTEGRATED_PR,
        "render_artifact_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "instantaneous_streamlines_are_material_paths": False,
        "sampled_vorticity_is_pde_acceptance_evidence": False,
        "canonical_velocity_export_ready": True,
        "canonical_velocity_api": project_status.get("velocity_api"),
        "canonical_thresholds_unchanged": True,
        "pde_validated": False,
    }


def main() -> None:
    print(json.dumps(audit_render_readiness_scope(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
