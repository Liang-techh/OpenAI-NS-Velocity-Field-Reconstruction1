"""Rank the frozen ST052-M #587/#652 pair by a public-source directional fingerprint.

Issue #666 freezes this comparison before this module is executed.  The objective
uses only qualitative statements on OpenAI's public 2026-09-08 Navier--Stokes
page: the vortex spirals inward, becomes increasingly elongated, the central
region shrinks, and the visualization trajectories show inward spiraling plus
axial stretching.

No image-derived numeric target, hidden coefficient, camera, frame-time map,
parameter fit, basis growth, pressure/forcing change, or common rescaling is
introduced.  The older #652 target-free Pareto rejection remains true even if
this narrower public-source directional comparison prefers the same exact child.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import agent7_st052m_combined_witness_grid_robustness as grid

TASK_ID = "CR003-ST052M-PUBLIC-SOURCE-DIRECTIONAL-FINGERPRINT-094"
PREREG_ISSUE = 666
PUBLIC_SOURCE_URL = "https://openai.com/index/navier-stokes-solution/"
PUBLIC_SOURCE_DATE = "2026-09-08"
SOURCE_GRID_PR = 659
SOURCE_GRID_HEAD = "afe2e4ba92fbeb8ef82d2a6fcd3d3a706d2a2c90"
SOURCE_WITNESS_PR = 652
SOURCE_WITNESS_HEAD = "b0d98d0d5fdb3ca6a438fdf6ac0dbaaef7781faf"
SOURCE_TEMPORAL_PR = 587
SOURCE_TEMPORAL_HEAD = "0b93819095f6c8576a7571bdc2d2fbef4154944d"
DIVERGENCE_MAX = 1.0e-5
CONTRACTION_RETENTION_TOL = 1.0e-10

PUBLIC_SOURCE_MAPPING = {
    "axial_stretching": {
        "public_statement": "trajectories show inward spiraling and axial stretching; the vortex becomes increasingly elongated",
        "observable": "child-minus-#587 enstrophy aspect increment at t=.50",
        "frozen_grids": [33, 41, 49],
        "rule": "strictly positive on every frozen grid",
    },
    "inward_spiraling": {
        "public_statement": "the vortex spirals inward; trajectories show inward spiraling",
        "observable": "48-path DOP853 mean-absolute-turns fidelity and inward-path count",
        "rule": "turns fidelity strictly improves and inward-path count does not decrease",
    },
    "central_contraction": {
        "public_statement": "the central region shrinks",
        "observable": "material-path contraction magnitude relative to redistributed control",
        "rule": "challenger contraction stays positive and within 1e-10 absolute of #587",
    },
}

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "witness_retuned": False,
    "parameter_grid_scan_performed": False,
    "optimization_performed": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
    "public_image_numeric_target_used": False,
    "pixel_similarity_objective_used": False,
    "hidden_openai_parameters_inferred": False,
    "visual_acceptance_threshold_defined": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def directional_decision(grid_report: dict, witness_report: dict) -> dict:
    """Apply exactly the qualitative partial order frozen in issue #666."""
    aspect_by_grid = {
        int(k): float(v)
        for k, v in grid_report["grid_robustness_decision"][
            "aspect_increment_by_resolution"
        ].items()
    }
    expected_grids = tuple(grid.GRID_RESOLUTIONS)
    if tuple(sorted(aspect_by_grid)) != tuple(sorted(expected_grids)):
        raise RuntimeError("frozen morphology grid set drifted")
    axial_stretching = bool(all(aspect_by_grid[n] > 0.0 for n in expected_grids))

    paths = witness_report["path_relative_to_control"]
    linear_paths = paths["linear"]
    child_paths = paths["child"]
    turns_increment = float(
        child_paths["mean_absolute_turns_relative"]
        - linear_paths["mean_absolute_turns_relative"]
    )
    inward_path_count_change = int(
        child_paths["inward_path_count_delta"] - linear_paths["inward_path_count_delta"]
    )
    inward_spiraling = bool(turns_increment > 0.0 and inward_path_count_change >= 0)

    linear_contraction = float(linear_paths["contraction_magnitude_relative"])
    child_contraction = float(child_paths["contraction_magnitude_relative"])
    contraction_difference = float(child_contraction - linear_contraction)
    central_contraction_retained = bool(
        child_contraction > 0.0
        and abs(contraction_difference) <= CONTRACTION_RETENTION_TOL
    )

    structure = witness_report["structure_preflight"]
    support_exact = bool(float(structure["support_max_abs"]) == 0.0)
    divergence_ok = bool(float(structure["divergence_fd_max"]) <= DIVERGENCE_MAX)
    structural_guard = bool(support_exact and divergence_ok)

    preferred = bool(
        axial_stretching
        and inward_spiraling
        and central_contraction_retained
        and structural_guard
    )
    return {
        "aspect_increment_by_resolution": aspect_by_grid,
        "axial_stretching_guard": axial_stretching,
        "turns_fidelity_increment": turns_increment,
        "inward_path_count_change": inward_path_count_change,
        "inward_spiraling_guard": inward_spiraling,
        "linear_contraction_magnitude_relative": linear_contraction,
        "child_contraction_magnitude_relative": child_contraction,
        "contraction_difference": contraction_difference,
        "contraction_retention_tolerance": CONTRACTION_RETENTION_TOL,
        "central_contraction_retained_guard": central_contraction_retained,
        "support_max_abs": float(structure["support_max_abs"]),
        "support_exact_zero_guard": support_exact,
        "divergence_fd_max": float(structure["divergence_fd_max"]),
        "divergence_max": DIVERGENCE_MAX,
        "divergence_guard": divergence_ok,
        "structural_guard": structural_guard,
        "public_source_directionally_preferred": preferred,
    }


def run(out: Path) -> dict:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Exact #659 rerun also reruns the exact #652 evaluator once and writes the
    # corresponding source_witness_reference.json beside its own report.
    grid_path = out.with_name("grid_robustness_reference.json")
    grid_report = grid.run(grid_path)
    witness_path = grid_path.with_name("source_witness_reference.json")
    witness_report = json.loads(witness_path.read_text())

    if grid_report["source_652_formal_rejection_remains_binding"] is not True:
        raise RuntimeError("#652 formal target-free rejection was not preserved")
    if witness_report["combined_child_target_free_pareto_verified"] is not False:
        raise RuntimeError("#652 target-free Pareto verdict drifted")
    if witness_report["frozen_witness"]["swirl_a"] != grid.witness.SWIRL_A:
        raise RuntimeError("#652 frozen compact-swirl coordinate drifted")
    if abs(
        float(witness_report["frozen_witness"]["shoulder_lambda"])
        - grid.witness.SHOULDER_LAMBDA
    ) > 1.0e-15:
        raise RuntimeError("#652 frozen shoulder timing coordinate drifted")

    decision = directional_decision(grid_report, witness_report)
    preferred = bool(decision["public_source_directionally_preferred"])

    secondary = {
        "target_free_tip_thinning_increment": float(
            witness_report["nonlinear_oriented_desirability_increments_vs_linear"][
                "tip_thinning_increment"
            ]
        ),
        "target_free_axial_pair_fidelity_increment": float(
            witness_report["nonlinear_oriented_desirability_increments_vs_linear"][
                "axial_pair_fidelity_increment"
            ]
        ),
        "midpoint_energy_child_relative_to_linear": float(
            witness_report["midpoint_energy"]["child_relative_to_linear"]
        ),
        "late_morphology_relative_to_control": witness_report[
            "late_morphology_relative_to_control_descriptive_only"
        ],
        "note": (
            "These are autonomous diagnostics only. Tip thinning is not a separate "
            "requirement stated by the frozen public OpenAI prose/caption and cannot "
            "veto or promote this public-source directional comparison."
        ),
    }

    routing = (
        "exact #652 is preferred over #587 only under the frozen public-source qualitative directional fingerprint; route this exact child to independent fixed-render/save-load comparison without retuning, while preserving the #652 target-free rejection and all visual/PDE/exact-field false states"
        if preferred
        else "retain #587 for the frozen public-source directional fingerprint; do not retune #652 or add a basis from this comparison"
    )

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "public_source": {
            "url": PUBLIC_SOURCE_URL,
            "date": PUBLIC_SOURCE_DATE,
            "mapping": PUBLIC_SOURCE_MAPPING,
            "qualitative_only": True,
        },
        "candidate_pair": {
            "baseline": {"pr": SOURCE_TEMPORAL_PR, "head": SOURCE_TEMPORAL_HEAD},
            "challenger": {"pr": SOURCE_WITNESS_PR, "head": SOURCE_WITNESS_HEAD},
            "grid_audit": {"pr": SOURCE_GRID_PR, "head": SOURCE_GRID_HEAD},
        },
        "source_652_formal_target_free_rejection_remains_binding": True,
        "grid_659_tip_margin_not_grid_resolved": bool(
            not grid_report["grid_robustness_decision"][
                "pareto_rejection_grid_resolved"
            ]
        ),
        "public_source_directional_decision": decision,
        "secondary_autonomous_diagnostics": secondary,
        "public_source_directionally_preferred": preferred,
        "routing": routing,
        **TRUTH,
    }
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    r = run(args.out)
    print("public_source_directionally_preferred=", r["public_source_directionally_preferred"])
    print("decision=", r["public_source_directional_decision"])
    print("secondary=", r["secondary_autonomous_diagnostics"])
    print("routing=", r["routing"])


if __name__ == "__main__":
    main()
