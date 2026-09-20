"""Fail-closed CR002 audit for the ST052-M reservoir visual-target scope.

PR #775 uses a useful but project-chosen numerical development screen.  This
module prevents those z/r/time coordinates and comparator gates from being
relabelled as public OpenAI numerical observations, independent visual
validation, or PDE evidence.

The audit is representation/claim governance only: it does not import or run
research fitting code and it changes no velocity field.
"""
from __future__ import annotations

import ast
import copy
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = REPO_ROOT / "configs/st052m_reservoir_visual_target_provenance_contract.json"
SOURCE_PATH = REPO_ROOT / "experiments/root_st052/agent7_st052m_outer_reservoir_nonlinear_replay.py"
CONSTRAINTS_PATH = REPO_ROOT / "configs/constraints.json"

EXPECTED_TASK_ID = "CR002-ST052M-RESERVOIR-VISUAL-TARGET-PROVENANCE-111"
EXPECTED_SOURCE_HEAD = "93887a59729d22113badf2ae4dac2f7d868e2703"
EXPECTED_PUBLIC_AUTHORITY_HEAD = "a6d713e3679024d58f2650e0a75ef2092a8dddff"
EXPECTED_PUBLIC_CONFIG_BLOB = "a222dbf88c4cf61e1dcc21e33cecefb7bb1aa4d1"
EXPECTED_CLASSES = [
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
]


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _literal_assignments(source_text: str) -> dict[str, Any]:
    tree = ast.parse(source_text)
    out: dict[str, Any] = {}
    for node in tree.body:
        if not isinstance(node, (ast.Assign, ast.AnnAssign)):
            continue
        if isinstance(node, ast.Assign):
            if len(node.targets) != 1 or not isinstance(node.targets[0], ast.Name):
                continue
            name = node.targets[0].id
            value_node = node.value
        else:
            if not isinstance(node.target, ast.Name) or node.value is None:
                continue
            name = node.target.id
            value_node = node.value
        try:
            out[name] = ast.literal_eval(value_node)
        except (ValueError, TypeError):
            continue
    return out


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def load_default_inputs() -> tuple[dict[str, Any], str, dict[str, Any]]:
    return (
        _load_json(CONTRACT_PATH),
        SOURCE_PATH.read_text(encoding="utf-8"),
        _load_json(CONSTRAINTS_PATH),
    )


def audit(
    contract: dict[str, Any] | None = None,
    source_text: str | None = None,
    constraints: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Validate the frozen source/autonomous/pending claim boundary.

    Optional in-memory arguments exist so mutation regressions can demonstrate
    fail-closed behavior without changing repository files.
    """
    if contract is None or source_text is None or constraints is None:
        default_contract, default_source, default_constraints = load_default_inputs()
        contract = default_contract if contract is None else contract
        source_text = default_source if source_text is None else source_text
        constraints = default_constraints if constraints is None else constraints

    contract = copy.deepcopy(contract)
    assignments = _literal_assignments(source_text)

    _require(contract.get("task_id") == EXPECTED_TASK_ID, "Agent-6 task identity drifted")
    audited = contract.get("audited_increment", {})
    _require(audited.get("pr") == 775, "audited PR drifted")
    _require(audited.get("head") == EXPECTED_SOURCE_HEAD, "audited #775 exact head drifted")
    _require(
        audited.get("task_id") == "CR003-ST052M-OUTER-RESERVOIR-NONLINEAR-REPLAY-110",
        "audited Agent-7 task drifted",
    )

    authority = contract.get("public_observable_authority", {})
    _require(authority.get("pr") == 773, "public-observable authority PR drifted")
    _require(authority.get("head") == EXPECTED_PUBLIC_AUTHORITY_HEAD, "public-observable authority head drifted")
    _require(authority.get("config_blob_sha") == EXPECTED_PUBLIC_CONFIG_BLOB, "public-observable contract blob drifted")
    _require(
        authority.get("integration_state_at_freeze") == "open_external_pr_not_live_on_constrained_integration",
        "queued/open external public-observable evidence may not be laundered into live integration",
    )
    _require(authority.get("public_numeric_targets_available") is False, "public numeric targets were invented")
    _require(
        authority.get("relevant_qualitative_observable_ids")
        == ["inward_spiraling_trajectories", "axial_stretching_trajectories"],
        "qualitative observable routing drifted",
    )

    _require(contract.get("classification_vocabulary") == EXPECTED_CLASSES, "four-class vocabulary drifted")
    classification = contract.get("classification", {})
    _require(set(classification) == set(EXPECTED_CLASSES), "classification buckets drifted")
    autonomous_text = "\n".join(classification["autonomous_design"])
    pending_text = "\n".join(classification["pending_unknown"])
    public_text = "\n".join(classification["public_source_fact"])
    _require("abs(z)=1.55 and 1.75" in autonomous_text, "visible z probes must remain autonomous design")
    _require("radii 0.6, 0.9, 1.2" in autonomous_text, "probe radii must remain autonomous design")
    _require("times 0.375, 0.50, 0.625, 0.75" in autonomous_text, "probe times must remain autonomous design")
    _require("no numerical target" in public_text, "public qualitative evidence acquired a numeric target")
    _require("OpenAI-render location" in pending_text, "coordinate registration must remain pending/unknown")
    _require("frame-to-physical-time mapping" in pending_text, "frame/time registration must remain pending/unknown")

    _require(assignments.get("TASK_ID") == audited.get("task_id"), "#775 source task ID drifted")
    _require(assignments.get("SOURCE_PARENT_PR") == 767, "#775 structural parent drifted")
    _require(assignments.get("RELATED_FRESH_HOLDOUT_PR") == 714, "fresh-holdout routing drifted")
    _require(tuple(assignments.get("VISIBLE_ABS_Z", ())) == (1.55, 1.75), "#775 visible z probes drifted")
    _require(tuple(assignments.get("RESERVOIR_ABS_Z", ())) == (1.85, 1.95), "#775 reservoir z probes drifted")
    _require(tuple(assignments.get("PROBE_RADII", ())) == (0.6, 0.9, 1.2), "#775 radial probes drifted")
    _require(tuple(assignments.get("PROBE_TIMES", ())) == (0.375, 0.5, 0.625, 0.75), "#775 time probes drifted")

    frozen = contract.get("frozen_numeric_screen", {})
    _require(frozen.get("visible_abs_z") == [1.55, 1.75], "contract/source visible z mismatch")
    _require(frozen.get("outer_reservoir_abs_z") == [1.85, 1.95], "contract/source reservoir z mismatch")
    _require(frozen.get("probe_radii") == [0.6, 0.9, 1.2], "contract/source radius mismatch")
    _require(frozen.get("probe_times") == [0.375, 0.5, 0.625, 0.75], "contract/source time mismatch")
    _require(frozen.get("visible_lobe_abs_z_interval") == [1.0, 1.8], "visible-lobe design interval drifted")
    _require(frozen.get("return_reservoir_abs_z_interval") == [1.8, 2.0], "return-reservoir design interval drifted")
    _require(frozen.get("numeric_screen_provenance") == "autonomous_design", "numeric visual screen was laundered into source fact")
    _require(frozen.get("public_numeric_target") is None, "public numerical target must remain null")
    _require(frozen.get("independent_visual_validation") is False, "development probe set was relabelled independent")
    _require(frozen.get("global_visual_validation") is False, "targeted probe set was relabelled global visual validation")

    _require(
        '"scope": "targeted development morphology fingerprint; not blind/global visual validation"' in source_text,
        "#775 targeted-development scope marker drifted",
    )
    _require(
        "outer return-flow location is a project representation choice; no OpenAI/public numeric threshold exists"
        in source_text,
        "#775 outer-reservoir public-target disclaimer drifted",
    )
    _require('"visible_abs_z_interval": [1.0, 1.8]' in source_text, "#775 visible-lobe interval drifted")
    _require('"return_reservoir_abs_z_interval": [1.8, 2.0]' in source_text, "#775 return-reservoir interval drifted")

    truth = assignments.get("TRUTH")
    _require(isinstance(truth, dict), "#775 TRUTH dictionary unavailable")
    required_false = [
        "canonical_velocity_changed",
        "saved_velocity_changed",
        "production_candidate_selected",
        "fresh_714_path_data_used",
        "parameter_grid_scan_performed",
        "optimization_performed",
        "post_result_damping_or_retuning_performed",
        "pressure_or_force_changed",
        "held_out_pde_residual_evaluated",
        "parent_pde_receipt_transferred",
        "public_image_numeric_target_used",
        "pixel_similarity_objective_used",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ]
    for key in required_false:
        _require(truth.get(key) is False, f"#775 truth boundary promoted {key}")
    _require(truth.get("targeted_visualization_direction_fingerprint_used") is True, "targeted development fingerprint marker lost")
    _require(truth.get("historical_development_paths_used_for_coefficient") is True, "selection-data provenance drifted")

    boundary = contract.get("truth_boundary", {})
    for key in (
        "canonical_velocity_changed_by_this_audit",
        "velocity_export_ready_promoted_by_this_audit",
        "visualization_ready",
        "visual_correspondence_verified",
        "source_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(boundary.get(key) is False, f"Agent-6 audit improperly promoted {key}")

    shortcut = contract.get("anti_shortcut", {})
    for key, value in shortcut.items():
        _require(value is False, f"forbidden shortcut enabled: {key}")

    # CR001 nonmutation guard: exact frozen scientific contract remains authoritative.
    _require(constraints.get("nu") == 0.01, "CR001 viscosity drifted")
    domain = constraints.get("domain", {})
    _require(domain.get("physical") == "R^3", "CR001 physical domain drifted")
    _require(domain.get("evaluation_box") == [[-2, 2], [-2, 2], [-2, 2]], "CR001 evaluation box drifted")
    _require(domain.get("support") == "r < 2 and abs(z) < 2", "CR001 support drifted")
    _require(domain.get("time_interval") == [0.25, 0.75], "CR001 time interval drifted")
    forcing = constraints.get("forcing", {})
    _require(forcing.get("mode") == "restricted_two_parameter_family", "CR001 forcing family drifted")
    _require(forcing.get("parameters") == {"a": [0.0, 10.0], "c": [0.0, 10.0]}, "CR001 force bounds drifted")
    restriction = forcing.get("restriction", "")
    _require("No residual-dependent basis or pointwise free force" in restriction, "free-force shortcut appeared")
    nontrivial = constraints.get("nontriviality", {})
    _require(nontrivial.get("reference_energy") == 1.0, "CR001 reference energy drifted")
    _require(nontrivial.get("reference_energy_abs_tolerance") == 0.001, "CR001 energy tolerance drifted")
    validation = constraints.get("validation", {})
    _require(validation.get("seed") == 914027, "CR001 validation seed drifted")
    _require(validation.get("held_out_points") == 4096, "CR001 held-out count drifted")
    _require(validation.get("derivative_steps") == [0.02, 0.01, 0.005], "CR001 derivative ladder drifted")
    _require(validation.get("quadrature_orders_per_axis") == [24, 48, 96], "CR001 quadrature ladder drifted")
    thresholds = validation.get("thresholds", {})
    _require(thresholds.get("divergence_max") == 1e-5 and thresholds.get("divergence_L2") == 1e-5, "CR001 divergence gates drifted")
    _require(thresholds.get("pde_residual_max") == 1e-3 and thresholds.get("pde_residual_L2") == 1e-3, "CR001 momentum gates drifted")
    _require("changing thresholds requires a new experiment version" in validation.get("failure_policy", ""), "post-hoc threshold relaxation became possible")

    return {
        "task_id": EXPECTED_TASK_ID,
        "audited_pr": 775,
        "audited_head": EXPECTED_SOURCE_HEAD,
        "numeric_probe_provenance": "autonomous_design",
        "public_numeric_target": None,
        "targeted_development_evidence_only": True,
        "independent_visual_validation": False,
        "global_visual_validation": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "cr001_unchanged": True,
    }


def main() -> None:
    print(json.dumps(audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
