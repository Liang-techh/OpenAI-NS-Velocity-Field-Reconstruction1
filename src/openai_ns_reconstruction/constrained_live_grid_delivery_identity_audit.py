"""CR002 audit for the live ST052-M grid-delivery identity contract.

This audit reconciles the already-integrated ST052-M whole-child capsule,
NPZ/MAT exporter, and GNU Octave software smoke without destructively replacing
the legacy sampled-grid governance contract still consumed by older fail-closed
audits. It deliberately does not promote finite-grid or rendering evidence into
continuous-function, PDE, visual-correspondence, or paper/OpenAI-field claims.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from . import st052_grid_export as grid_export
from . import st052_linear_temporal_capsule as whole_capsule


CONTRACT_FILENAME = "live_st052_grid_delivery_identity_contract.json"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise AssertionError(f"expected JSON object: {path}")
    return obj


def _assert_equal(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise AssertionError(f"{label} drifted: {actual!r} != {expected!r}")


def audit_live_grid_delivery_identity(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Fail closed on live delivery, representation, truth-state, or CR001 drift."""
    root = Path(repo_root) if repo_root is not None else _repo_root()
    contract = _load_json(root / "configs" / CONTRACT_FILENAME)
    constraints = _load_json(root / "configs" / "constraints.json")

    if contract.get("schema_version") != 2:
        raise AssertionError("live ST052 velocity-grid contract must use reconciled schema v2")
    if contract.get("task_id") != "CR002-VELOCITY-GRID-DELIVERY-IDENTITY-073":
        raise AssertionError("unexpected live ST052 velocity-grid governance task id")

    delivery = contract["current_st052_grid_delivery"]
    runtime = delivery["whole_child_runtime"]
    truth = whole_capsule.TRUTH_BOUNDARY
    _assert_equal(runtime["save_load_ready_with_exact_source_runtime"],
                  truth["whole_child_save_load_ready_with_exact_source_runtime"],
                  "ST052 source-runtime save/load state")
    _assert_equal(runtime["exact_source_runtime_identity_closed"],
                  truth["exact_source_runtime_identity_closed"],
                  "ST052 exact-source runtime identity")
    _assert_equal(runtime["historical_module_cache_isolated"],
                  truth["historical_module_cache_isolated"],
                  "ST052 historical module cache isolation")
    _assert_equal(runtime["standalone_package_parent_runtime_ready"],
                  truth["standalone_package_parent_runtime_ready"],
                  "ST052 standalone package parent runtime")
    _assert_equal(runtime["velocity_export_ready"],
                  truth["velocity_export_ready"],
                  "ST052 velocity-export readiness")

    _assert_equal(delivery["export_schema"], grid_export.SCHEMA, "grid export schema")
    _assert_equal(delivery["producer_task_id"], grid_export.TASK_ID, "grid export task id")
    _assert_equal(delivery["filenames"]["npz"], grid_export.NPZ_FILENAME, "NPZ filename")
    _assert_equal(delivery["filenames"]["mat"], grid_export.MAT_FILENAME, "MAT filename")
    _assert_equal(delivery["filenames"]["manifest"], grid_export.MANIFEST_FILENAME, "manifest filename")
    _assert_equal(delivery["spatial_points_per_axis"], grid_export.DEFAULT_GRID_POINTS, "grid points")
    _assert_equal(delivery["times"], list(grid_export.DEFAULT_TIMES), "grid times")
    _assert_equal(delivery["domain"], list(grid_export.DEFAULT_DOMAIN), "grid domain")
    _assert_equal(delivery["velocity_component_shape"],
                  [len(grid_export.DEFAULT_TIMES),
                   grid_export.DEFAULT_GRID_POINTS,
                   grid_export.DEFAULT_GRID_POINTS,
                   grid_export.DEFAULT_GRID_POINTS],
                  "velocity component shape")
    _assert_equal(delivery["array_layout"], "u,v,w: [time,x,y,z]", "array layout")
    _assert_equal(delivery["meshgrid_indexing"], "ij", "meshgrid indexing")
    _assert_equal(delivery["callable_node_replay_tolerance"],
                  grid_export.CALLABLE_GRID_PARITY_ATOL,
                  "callable node replay tolerance")
    _assert_equal(delivery["callable_node_replay_tolerance_scope"],
                  "engineering_serialization_replay_only",
                  "callable node replay tolerance scope")
    if delivery["payload_checksum_bound_to_whole_candidate_identity"] is not True:
        raise AssertionError("grid payload must remain checksum-bound to whole-candidate identity")
    if delivery["npz_mat_numerical_payload_exact_equal_required"] is not True:
        raise AssertionError("NPZ/MAT numerical payload equality requirement drifted")

    exporter_text = (
        root / "src" / "openai_ns_reconstruction" / "st052_grid_export.py"
    ).read_text(encoding="utf-8")
    for token in (
        '"array_layout": "u,v,w: [time,x,y,z]"',
        '"meshgrid_indexing": "ij"',
        '"grid_export_payload_checksum_bound_to_whole_candidate": True',
        "It is not a PDE, visualization, or source-correspondence acceptance gate.",
    ):
        if token not in exporter_text:
            raise AssertionError(f"live exporter evidence missing: {token}")

    representation = contract["representation_scope"]
    if representation.get("sampled_grid_is_finite_representation") is not True:
        raise AssertionError("finite sampled-grid representation fact must remain true")
    for key in (
        "sampled_grid_is_continuum_velocity_api",
        "node_equality_implies_continuous_function_identity",
        "off_grid_interpolant_specified",
        "off_grid_error_bound_established",
        "continuous_derivative_equivalence_established",
        "grid_derivatives_are_cr001_pde_acceptance_evidence",
        "grid_nodes_may_replace_cr001_held_out_sample",
        "five_grid_times_may_replace_six_cr001_validation_times",
        "visual_interpolation_may_be_called_exact_candidate_without_contract",
        "sampled_grid_may_be_called_exact_openai_field",
    ):
        if representation.get(key) is not False:
            raise AssertionError(f"unverified finite-grid inference promoted: {key}")

    octave = contract["octave_consumer"]
    octave_text = (root / "scripts" / "smoke_st052_grid_octave.m").read_text(encoding="utf-8")
    for token in (
        "size(S.u), [5, 33, 33, 33]",
        "expected_t = [0.25; 0.375; 0.5; 0.625; 0.75]",
        "actual_octave_runtime_executed=true",
        "visualization_ready=false",
        "visual_correspondence_verified=false",
        "pde_validated=false",
    ):
        if token not in octave_text:
            raise AssertionError(f"integrated Octave smoke evidence missing: {token}")
    for key in (
        "actual_octave_runtime_executed",
        "mat_file_load_smoke_verified",
        "fixed_png_render_smoke_verified",
        "shape_time_domain_finite_checks_executed",
        "centered_difference_omega_z_is_software_diagnostic_only",
    ):
        if octave.get(key) is not True:
            raise AssertionError(f"integrated Octave software fact must remain true: {key}")
    for key in (
        "actual_matlab_runtime_executed",
        "layout_metadata_field_enforced",
        "whole_candidate_identity_metadata_enforced",
        "grid_payload_checksum_metadata_enforced",
        "component_basis_machine_bound_from_metadata",
        "derived_grid_vorticity_scientific_semantics_verified",
    ):
        if octave.get(key) is not False:
            raise AssertionError(f"unverified Octave/metadata claim promoted: {key}")
    for forbidden_current_token in (
        "array_layout",
        "whole_candidate_identity_sha256",
        "grid_payload_sha256",
    ):
        if forbidden_current_token in octave_text:
            raise AssertionError(
                "Octave consumer metadata behavior changed; re-audit live ST052 grid semantics"
            )

    boundary = contract["truth_boundary"]
    for key in (
        "canonical_velocity_export_ready",
        "st052_source_runtime_backed_callable_save_load_ready",
        "st052_grid_export_materialized",
        "st052_octave_software_smoke_passed",
    ):
        if boundary.get(key) is not True:
            raise AssertionError(f"integrated delivery fact regressed: {key}")
    for key in (
        "st052_velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        if boundary.get(key) is not False:
            raise AssertionError(f"scientific/delivery truth state promoted without evidence: {key}")

    canonical = contract["canonical_callable_delivery"]
    if canonical.get("velocity_export_ready") is not True:
        raise AssertionError("canonical Eq45 callable delivery readiness regressed")
    if canonical.get("pde_failure_blocks_callable_velocity_delivery") is not False:
        raise AssertionError("PDE failure must not be turned into a callable-delivery blocker")
    if canonical.get("st052_grid_or_octave_delivery_replaces_canonical_eq45") is not False:
        raise AssertionError("ST052 grid/software delivery cannot silently replace canonical Eq45")

    frozen = contract["canonical_cr001"]
    _assert_equal(constraints["nu"], frozen["nu"], "nu")
    _assert_equal(constraints["domain"]["physical"], frozen["physical_domain"], "physical domain")
    _assert_equal(constraints["domain"]["evaluation_box"], frozen["evaluation_box"], "evaluation box")
    _assert_equal(constraints["domain"]["support"], frozen["support"], "support")
    _assert_equal(constraints["domain"]["time_interval"], frozen["time_interval"], "time interval")
    _assert_equal(constraints["forcing"]["mode"], frozen["force_mode"], "forcing mode")
    _assert_equal(constraints["forcing"]["parameters"], frozen["force_bounds"], "forcing bounds")
    if "No residual-dependent basis or pointwise free force" not in constraints["forcing"]["restriction"]:
        raise AssertionError("residual-defined/free-force prohibition drifted")
    _assert_equal(constraints["nontriviality"]["reference_energy"],
                  frozen["reference_energy"], "reference energy")
    _assert_equal(constraints["nontriviality"]["reference_energy_abs_tolerance"],
                  frozen["reference_energy_abs_tolerance"], "reference energy tolerance")
    if "reject collapsed candidates" not in constraints["nontriviality"]["enforcement"]:
        raise AssertionError("amplitude-collapse rejection drifted")
    validation = constraints["validation"]
    _assert_equal(validation["seed"], frozen["validation_seed"], "validation seed")
    _assert_equal(validation["held_out_points"], frozen["held_out_points"], "held-out points")
    _assert_equal(validation["times"], frozen["validation_times"], "validation times")
    _assert_equal(validation["derivative_steps"], frozen["derivative_steps"], "derivative ladder")
    _assert_equal(validation["quadrature_orders_per_axis"],
                  frozen["quadrature_orders_per_axis"], "quadrature ladder")
    thresholds = validation["thresholds"]
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _assert_equal(thresholds[key], frozen[key], key)

    return {
        "task_id": contract["task_id"],
        "live_st052_grid_delivery_checked": True,
        "octave_software_scope_checked": True,
        "continuous_representation_boundary_checked": True,
        "cr001_checked_unchanged": True,
        "truth_boundary": boundary,
    }
