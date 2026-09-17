"""Govern representation semantics for the quartic Eq45 visualization smoke.

The PR #177 smoke saves all three Cartesian velocity components on a fixed y=0
plane, but its line overlay integrates only the two-component (u,w) meridional
field.  For an axisymmetric field with swirl, that 2-D line overlay is not a
three-dimensional streamline and is not evidence of public visual correspondence.

This module is claim-boundary governance only.  It changes no velocity, support
transform, force, threshold, candidate routing, or scientific acceptance state.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .constrained_eq45_supported_phi10_quartic_visualization_smoke import (
    DEFAULT_EXTENT,
    DEFAULT_GRID_SIZE,
    DEFAULT_TIMES,
    SCHEMA as SMOKE_SCHEMA,
    governed_quartic_candidate,
    meridional_velocity_slice,
)

EXPECTED_SCHEMA = "eq45_phi10_quartic_visualization_smoke_scope_v1"
EXPECTED_SOURCE_VOCABULARY = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
EXPECTED_CANDIDATE_SHA256 = (
    "03fdae73170469ae1160297b489b1b83a27adddfc941119e116a21b98bd15049"
)
EXPECTED_PARENT_SHA256 = (
    "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
)
EXPECTED_IMPLEMENTATION_HEAD = "e7ab61b4485be1779e048608605f92c53dcd334a"


def _load(path: str | Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _same_float(left: Any, right: Any, *, atol: float = 1.0e-12) -> bool:
    try:
        return abs(float(left) - float(right)) <= atol
    except (TypeError, ValueError):
        return False


def audit_quartic_visualization_smoke_scope(
    *,
    contract_path: str | Path | None = None,
    constraints_path: str | Path | None = None,
    delivery_state_path: str | Path | None = None,
) -> dict[str, Any]:
    """Fail closed if the 2-D smoke is promoted beyond its representation scope."""
    root = Path(__file__).resolve().parents[2]
    contract = _load(
        contract_path
        or root / "configs/eq45_phi10_quartic_visualization_smoke_scope.json"
    )
    constraints = _load(constraints_path or root / "configs/constraints.json")
    delivery = _load(delivery_state_path or root / "configs/delivery_state_contract.json")

    _require(contract.get("schema") == EXPECTED_SCHEMA, "governance schema drifted")

    canonical_vocab = set(delivery["classification_vocabulary"])
    local_vocab = set(contract.get("source_vocabulary", []))
    _require(
        canonical_vocab == EXPECTED_SOURCE_VOCABULARY,
        "delivery-state source vocabulary drifted",
    )
    _require(
        local_vocab == EXPECTED_SOURCE_VOCABULARY,
        "smoke-governance source vocabulary drifted",
    )
    classifications = contract.get("source_classification", {})
    _require(
        set(classifications.values()) <= EXPECTED_SOURCE_VOCABULARY,
        "unknown source classification",
    )
    expected_classifications = {
        "callable_3d_time_varying_velocity_delivery": "user_requirement",
        "separate_velocity_visual_pde_and_exactness_states": "user_requirement",
        "fixed_y0_meridional_smoke_plane": "autonomous_design",
        "default_smoke_grid_extent_and_times": "autonomous_design",
        "speed_background": "autonomous_design",
        "two_component_u_w_line_overlay": "autonomous_design",
        "public_visual_geometry_time_camera_registration": "pending_unknown",
        "paper_exact_or_hidden_openai_identity": "pending_unknown",
    }
    _require(
        classifications == expected_classifications,
        "visualization-smoke source classification drifted",
    )

    binding = contract["cr001_binding"]
    domain = constraints["domain"]
    forcing = constraints["forcing"]
    nontriviality = constraints["nontriviality"]
    validation = constraints["validation"]
    thresholds = validation["thresholds"]
    _require(_same_float(binding["nu"], constraints["nu"]), "nu drifted")
    _require(binding["physical_domain"] == domain["physical"], "physical domain drifted")
    _require(binding["evaluation_box"] == domain["evaluation_box"], "evaluation box drifted")
    _require(binding["support"] == domain["support"], "support contract drifted")
    _require(binding["time_interval"] == domain["time_interval"], "time interval drifted")
    _require(binding["forcing_mode"] == forcing["mode"], "forcing mode drifted")
    _require(
        binding["forcing_parameters"] == sorted(forcing["parameters"]),
        "restricted forcing parameter set drifted",
    )
    _require(
        "No residual-dependent basis or pointwise free force" in forcing["restriction"],
        "restricted forcing no-free-force guard drifted",
    )
    _require(
        _same_float(binding["reference_energy_time"], nontriviality["reference_time"]),
        "reference energy time drifted",
    )
    _require(
        _same_float(binding["reference_energy"], nontriviality["reference_energy"]),
        "reference energy target drifted",
    )
    _require(
        _same_float(
            binding["reference_energy_abs_tolerance"],
            nontriviality["reference_energy_abs_tolerance"],
        ),
        "reference energy tolerance drifted",
    )
    _require(
        int(binding["held_out_points"]) == int(validation["held_out_points"]),
        "held-out sample count drifted",
    )
    _require(
        [float(x) for x in binding["derivative_steps"]]
        == [float(x) for x in validation["derivative_steps"]],
        "derivative ladder drifted",
    )
    threshold_pairs = (
        ("divergence_max_threshold", "divergence_max"),
        ("divergence_L2_threshold", "divergence_L2"),
        ("pde_residual_max_threshold", "pde_residual_max"),
        ("pde_residual_L2_threshold", "pde_residual_L2"),
    )
    for local_key, registered_key in threshold_pairs:
        _require(
            _same_float(binding[local_key], thresholds[registered_key]),
            f"{registered_key} threshold drifted",
        )

    smoke = contract["smoke_binding"]
    _require(smoke["implementation_pr"] == 177, "smoke implementation PR drifted")
    _require(
        smoke["implementation_head"] == EXPECTED_IMPLEMENTATION_HEAD,
        "smoke implementation head drifted",
    )
    _require(smoke["manifest_schema"] == SMOKE_SCHEMA, "smoke manifest schema drifted")
    _require(
        smoke["candidate_family"]
        == "Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate",
        "smoke candidate family drifted",
    )
    _require(
        smoke["candidate_sha256"] == EXPECTED_CANDIDATE_SHA256,
        "smoke candidate identity drifted",
    )
    _require(
        smoke["supported_parent_sha256"] == EXPECTED_PARENT_SHA256,
        "smoke parent identity drifted",
    )
    _require(
        smoke["observable"] == "fixed_physical_y0_meridional_velocity_slice",
        "smoke observable drifted",
    )
    _require(smoke["saved_component_order"] == ["u", "v", "w"], "saved component order drifted")
    _require(smoke["plane"] == {"axis": "y", "value": 0.0}, "smoke plane drifted")
    _require(
        [float(x) for x in smoke["default_times"]] == [float(x) for x in DEFAULT_TIMES],
        "smoke default times drifted",
    )
    _require(_same_float(smoke["default_extent"], DEFAULT_EXTENT), "smoke extent drifted")
    _require(int(smoke["default_grid_size"]) == int(DEFAULT_GRID_SIZE), "smoke grid size drifted")
    _require(
        smoke["render"] == "speed_background_plus_projected_u_w_streamlines",
        "smoke render mode drifted",
    )
    _require(smoke["line_overlay_components"] == ["u", "w"], "line overlay components drifted")
    _require(smoke["omitted_line_direction_component"] == "v", "omitted line component drifted")
    for key in (
        "public_image_used",
        "camera_or_registration_fitted",
        "visual_pass_threshold_defined",
    ):
        _require(smoke[key] is False, f"{key} must remain false")

    field = governed_quartic_candidate()
    _require(field.sha256 == EXPECTED_CANDIDATE_SHA256, "live quartic candidate identity drifted")
    _require(field.base_sha256 == EXPECTED_PARENT_SHA256, "live supported parent identity drifted")
    sampled = meridional_velocity_slice(field, 0.25, grid_size=9, extent=1.25)
    max_abs_v = float(np.max(np.abs(np.asarray(sampled["v"], dtype=float))))
    _require(max_abs_v > 1.0e-12, "swirl component unexpectedly vanished on the governed y=0 smoke sample")

    semantics = contract["representation_semantics"]
    expected_semantics = {
        "saved_slice_contains_full_cartesian_velocity_at_each_sample": True,
        "line_overlay_is_2d_integral_curve_of_y0_u_w_component_field": True,
        "line_overlay_is_true_3d_streamline": False,
        "line_overlay_is_full_3d_trajectory_projection": False,
        "omitted_v_may_be_nonzero_for_swirl_candidate": True,
        "line_overlay_may_be_called_meridional_or_poloidal_streamline_smoke": True,
        "line_overlay_may_be_called_true_3d_streamline": False,
    }
    _require(semantics == expected_semantics, "visualization-smoke representation semantics drifted")

    conclusions = contract["governed_conclusions"]
    _require(conclusions["smoke_may_establish_candidate_renderability"] is True, "renderability evidence drifted")
    _require(
        conclusions["candidate_local_rendering_may_proceed_while_pde_unvalidated"] is True,
        "PDE status may not become a candidate-local rendering blocker",
    )
    for key in (
        "smoke_may_establish_visualization_ready",
        "smoke_may_establish_visual_correspondence_verified",
        "smoke_may_establish_pde_validated",
        "smoke_may_establish_paper_exact",
        "smoke_may_establish_openai_field_identified",
        "pde_failure_or_pending_may_block_velocity_export",
    ):
        _require(conclusions[key] is False, f"{key} must remain false")

    states = contract["truth_states"]
    expected_states = {
        "velocity_export_ready": True,
        "visualization_candidate_only": True,
        "visualization_smoke_generated": True,
        "physical_support_validated": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    _require(states == expected_states, "visualization-smoke truth-state boundary drifted")

    required_forbidden = {
        "u_w_y0_line_overlay=>true_3d_streamline",
        "u_w_y0_line_overlay=>full_3d_trajectory_projection",
        "render_smoke_generated=>visualization_ready",
        "render_smoke_generated=>visual_correspondence_verified",
        "render_smoke_generated=>pde_validated",
        "visual_slice_similarity=>paper_exact_or_openai_field_identified",
        "pde_failed_or_pending=>velocity_export_not_allowed",
    }
    _require(set(contract["forbidden_inferences"]) == required_forbidden, "forbidden-inference set drifted")

    return {
        "candidate_sha256": field.sha256,
        "supported_parent_sha256": field.base_sha256,
        "sampled_y0_max_abs_v": max_abs_v,
        "saved_slice_contains_full_cartesian_velocity_at_each_sample": semantics[
            "saved_slice_contains_full_cartesian_velocity_at_each_sample"
        ],
        "line_overlay_is_true_3d_streamline": semantics["line_overlay_is_true_3d_streamline"],
        "visualization_smoke_generated": states["visualization_smoke_generated"],
        "visualization_ready": states["visualization_ready"],
        "visual_correspondence_verified": states["visual_correspondence_verified"],
        "pde_validated": states["pde_validated"],
        "velocity_export_ready": states["velocity_export_ready"],
    }


def main() -> int:
    print(json.dumps(audit_quartic_visualization_smoke_scope(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
