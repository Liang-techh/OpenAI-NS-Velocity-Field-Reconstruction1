"""Fail-closed CR002 governance for sampled velocity-grid delivery identity.

A sampled visualization grid is useful delivery evidence, but it is not the same
object as the analytic candidate, a continuum ``velocity(x,y,z,t)`` API, or a
PDE-validation receipt.  This module changes no field or scientific threshold.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

TASK_ID = "CR002-VELOCITY-GRID-DELIVERY-IDENTITY-062"
_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_ST048_S_SHA = "6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09"
_ST048_HEAD = "97695a86f85ce68fb4ae70c41fc81c904d655183"
_CANONICAL_DELIVERY_SHA = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
_CANONICAL_VELOCITY_API = "openai_ns_reconstruction.eq45_supported_delivery:velocity"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"expected object in {path}")
    return data


def load_contract() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "velocity_grid_delivery_identity_contract.json")


def load_constraints() -> dict[str, Any]:
    return _load_json(_repo_root() / "configs" / "constraints.json")


def load_project_status() -> dict[str, Any]:
    return _load_json(_repo_root() / "project_status.json")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    _require(actual == expected, f"{label} drifted: {actual!r} != {expected!r}")


def audit_velocity_grid_delivery_identity(
    contract: Mapping[str, Any] | None = None,
    constraints: Mapping[str, Any] | None = None,
    project_status: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit candidate/grid/file/API identity separation and CR001 invariants."""

    contract = deepcopy(dict(load_contract() if contract is None else contract))
    constraints = deepcopy(dict(load_constraints() if constraints is None else constraints))
    project_status = deepcopy(dict(load_project_status() if project_status is None else project_status))

    _require_equal(contract.get("schema_version"), 1, "contract schema")
    _require_equal(contract.get("task_id"), TASK_ID, "contract task id")

    snapshot = contract.get("snapshot", {})
    _require_equal(snapshot.get("active_integration_branch"), "codex/cr001-constraints", "integration branch")
    _require_equal(snapshot.get("observed_delivery_pr"), 404, "observed delivery PR")
    _require_equal(snapshot.get("observed_delivery_pr_state"), "open_unintegrated", "delivery PR state")
    _require_equal(snapshot.get("observed_delivery_pr_base"), "main", "delivery PR base")

    _require_equal(set(contract.get("allowed_source_classes", ())), _ALLOWED_CLASSES, "source-class vocabulary")
    source = contract.get("source_classification", {})
    _require(source and set(source.values()) <= _ALLOWED_CLASSES, "invalid source classification")
    for key in ("callable_velocity_x_y_z_t", "save_load_delivery"):
        _require_equal(source.get(key), "user_requirement", f"source class {key}")
    _require_equal(
        source.get("public_openai_visualization_observables"),
        "public_source_fact",
        "public visualization source class",
    )
    for key in (
        "st048_candidate_recipe_and_repository_measurements",
        "grid_resolution_times_format_and_checksum_scheme",
        "cr001_validation_protocol",
    ):
        _require_equal(source.get(key), "autonomous_design", f"source class {key}")
    for key in (
        "openai_hidden_numerical_velocity",
        "openai_hidden_frame_time_camera_mapping",
        "openai_hidden_parameters",
    ):
        _require_equal(source.get(key), "pending_unknown", f"source class {key}")

    layers = contract.get("identity_layers", {})
    analytic = layers.get("analytic_candidate", {})
    _require_equal(analytic.get("candidate_id"), "ST048-S", "analytic candidate id")
    _require_equal(analytic.get("source_pr"), 390, "analytic candidate source PR")
    _require_equal(analytic.get("source_head_sha"), _ST048_HEAD, "analytic candidate source head")
    _require_equal(analytic.get("archived_raw_candidate_sha256"), _ST048_S_SHA, "analytic candidate SHA")
    _require_equal(analytic.get("role"), "protocol_scoped_open_challenger", "analytic candidate role")
    _require_equal(analytic.get("canonical_velocity_replaced"), False, "challenger canonical promotion")
    _require_equal(analytic.get("pde_validated"), False, "challenger PDE state")

    grid = layers.get("sampled_grid", {})
    _require_equal(grid.get("producer_pr"), 404, "grid producer PR")
    _require_equal(grid.get("resolution"), 33, "grid resolution")
    _require_equal(grid.get("times"), [0.25, 0.5, 0.75], "grid times")
    _require_equal(grid.get("evaluation_box"), [[-2, 2], [-2, 2], [-2, 2]], "grid box")
    _require_equal(grid.get("component_order"), ["u", "v", "w"], "grid component order")
    _require_equal(grid.get("grid_sha_role"), "derived_sample_content_identity", "grid SHA role")
    _require_equal(grid.get("grid_sha_is_candidate_sha"), False, "grid/candidate SHA separation")
    _require_equal(
        grid.get("changing_resolution_or_times_changes_grid_identity_not_analytic_candidate_identity"),
        True,
        "grid/candidate identity invariance",
    )
    _require_equal(grid.get("sampled_grid_is_continuum_velocity_api"), False, "grid/continuum separation")
    _require_equal(
        grid.get("off_grid_or_unsampled_time_interpolation_requires_explicit_new_numerical_contract"),
        True,
        "interpolation contract requirement",
    )

    serialized = layers.get("serialized_file", {})
    _require_equal(serialized.get("format"), "NetCDF", "serialized format")
    _require_equal(
        serialized.get("file_bytes_identity_is_grid_content_identity"),
        False,
        "file/grid identity separation",
    )
    _require_equal(
        serialized.get("serialization_metadata_may_change_without_changing_sampled_values"),
        True,
        "metadata/value identity separation",
    )
    _require_equal(
        serialized.get("serialized_file_hash_may_replace_candidate_sha"),
        False,
        "file/candidate SHA separation",
    )

    delivery = layers.get("canonical_callable_delivery", {})
    _require_equal(delivery.get("candidate_family"), project_status.get("candidate_family"), "canonical delivery family")
    _require_equal(delivery.get("candidate_sha256"), project_status.get("candidate_sha256"), "canonical delivery SHA")
    _require_equal(delivery.get("candidate_sha256"), _CANONICAL_DELIVERY_SHA, "known canonical delivery SHA")
    _require_equal(delivery.get("velocity_api"), project_status.get("velocity_api"), "canonical velocity API")
    _require_equal(delivery.get("velocity_api"), _CANONICAL_VELOCITY_API, "known canonical velocity API")
    _require_equal(delivery.get("velocity_export_ready"), True, "canonical export state")
    _require_equal(
        delivery.get("open_unintegrated_grid_export_changes_canonical_velocity_api"),
        False,
        "grid/API promotion separation",
    )
    _require_equal(
        delivery.get("open_unintegrated_grid_export_changes_canonical_candidate"),
        False,
        "grid/canonical candidate promotion separation",
    )
    _require_equal(delivery.get("pde_failure_blocks_callable_velocity_delivery"), False, "delivery/PDE separation")

    evidence = contract.get("grid_evidence_scope", {})
    allowed = set(evidence.get("allowed_claims", ()))
    _require(
        "the frozen analytic candidate was sampled at the declared nodes and times" in allowed,
        "missing sampled-node claim",
    )
    _require(
        "the resulting grid is a cross-language visualization artifact when its source identity is preserved" in allowed,
        "missing cross-language artifact claim",
    )
    for key in (
        "grid_ci_success_implies_candidate_promotion",
        "grid_ci_success_implies_visual_correspondence",
        "grid_ci_success_implies_pde_validated",
        "grid_derivatives_are_cr001_pde_acceptance_evidence",
        "grid_nodes_may_replace_cr001_held_out_sample",
        "three_grid_times_may_replace_six_cr001_validation_times",
        "visual_interpolation_may_be_called_exact_candidate_without_contract",
        "sampled_grid_may_be_called_exact_openai_field",
    ):
        _require_equal(evidence.get(key), False, key)

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
        (
            cr001.get("reference_energy_abs_tolerance"),
            nontriviality.get("reference_energy_abs_tolerance"),
            "reference energy tolerance",
        ),
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
    _require(grid.get("times") != validation.get("times"), "render-grid times accidentally equal validation protocol")
    _require(grid.get("resolution") ** 3 != validation.get("held_out_points"), "render grid accidentally aliases held-out sample count")

    truth = contract.get("truth_boundary", {})
    states = project_status.get("states", {})
    _require_equal(truth.get("canonical_velocity_export_ready"), True, "canonical export truth")
    _require_equal(states.get("velocity_export_ready"), True, "project export truth")
    _require_equal(
        truth.get("st048_sampled_grid_integrated_as_canonical_delivery"),
        False,
        "ST048 sampled-grid integration truth",
    )
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require_equal(truth.get(key), False, f"contract truth {key}")
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
    ):
        _require_equal(mutation.get(key), False, key)

    return {
        "task_id": TASK_ID,
        "status": "governance_pass",
        "analytic_candidate": "ST048-S",
        "sampled_grid_is_separate_identity": True,
        "sampled_grid_is_not_continuum_velocity_api": True,
        "canonical_velocity_family": project_status.get("candidate_family"),
        "canonical_velocity_api": project_status.get("velocity_api"),
        "canonical_velocity_export_ready": True,
        "st048_canonical_promotion": False,
        "grid_derivatives_are_pde_evidence": False,
        "canonical_thresholds_unchanged": True,
        "pde_validated": False,
        "visual_correspondence_verified": False,
    }


def main() -> None:
    print(json.dumps(audit_velocity_grid_delivery_identity(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
