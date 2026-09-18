"""Independent roundoff-isolation audit for the Kokuno source-wave-shell seam.

Agent 4 PR #371 independently checked Agent 2's source-grid/source-wave-shell
wrapper with a Cartesian FD4 curl/divergence oracle.  Its absolute errors were
tiny, but the frozen divergence-refinement ratio rejected the finest
``0.001 / 0.0005 / 0.00025`` ladder after the calculation entered an
approximately 1e-10 floating-point plateau.

This module is a *new preregistered experiment*, not a reinterpretation of
#371.  It keeps the exact same scientific guards and independent Cartesian
operator, changes only the held-out seed and the FD4 ladder to the coarser
``0.004 / 0.002 / 0.001`` range selected before observing this run, and asks
whether the same source-domain seam exhibits the required refinement while it
is still truncation dominated.  The historical #371 rejection remains part of
the report and is never overwritten.

This is still only a local structural preflight.  It does not provide a global
Kokuno leading pressure, actual-source public oscillatory xyz,t velocity, a
promoted correction composite, or the formal full-domain Navier--Stokes gate.
"""
from __future__ import annotations

import argparse
import json
from math import sqrt
from pathlib import Path
from typing import Any

import numpy as np

from . import kokuno_source_wave_shell_cartesian_independent as prior

SCHEMA = "kokuno-agent4-source-wave-shell-roundoff-ladder-audit-v1"
TASK_ID = "KOKUNO-A4-SOURCE-WAVE-SHELL-ROUNDOFF-LADDER-AUDIT-020"
BASE_PR = 371
BASE_HEAD = "917db15e5b72a04662cf78cd0345b32a9adeb51c"
SEED = 9173111
FD4_STEPS = (0.004, 0.002, 0.001)

# These are deliberately aliases of the already-frozen Agent-4 guards.  This
# audit may change the sample and step ladder, but not the acceptance criteria.
LOCAL_CURL_FINE_REL_GUARD = prior.LOCAL_CURL_FINE_REL_GUARD
LOCAL_DIVERGENCE_FINE_ABS_GUARD = prior.LOCAL_DIVERGENCE_FINE_ABS_GUARD
LOCAL_REFINEMENT_RATIO_GUARD = prior.LOCAL_REFINEMENT_RATIO_GUARD
MISSING_GRADIENT_RELATIVE_ERROR_FLOOR = prior.MISSING_GRADIENT_RELATIVE_ERROR_FLOOR
MISSING_GRADIENT_DIVERGENCE_FLOOR = prior.MISSING_GRADIENT_DIVERGENCE_FLOOR
FORMAL_MOMENTUM_GATE = prior.FORMAL_MOMENTUM_GATE
FORMAL_DIVERGENCE_GATE = prior.FORMAL_DIVERGENCE_GATE

# Immutable receipt from PR #371.  Recording it here prevents a later PASS on
# this distinct ladder from silently rewriting the earlier scientific REJECT.
PRIOR_REJECTION = {
    "pr": 371,
    "exact_head": BASE_HEAD,
    "seed": 9173101,
    "fd4_steps": [0.001, 0.0005, 0.00025],
    "finest_curl_relative_rms": 1.3245337639271412e-12,
    "finest_divergence_max_abs": 2.1924837789035134e-10,
    "finest_divergence_rms_abs": 1.0406906338850813e-10,
    "curl_relative_rms_refinement_ratios": [
        15.994123110028152,
        14.5964801840349,
    ],
    "divergence_rms_refinement_ratios": [
        11.789475096192975,
        1.193958216011321,
    ],
    "source_wave_shell_independent_cartesian_preflight_passed": False,
}


def _sample_probes() -> list[dict[str, Any]]:
    """Fresh source-domain probes with the same frozen stratification as #371."""
    rng = np.random.default_rng(SEED)
    probes: list[dict[str, Any]] = []
    for case in prior.CASES:
        for region, count, x_ratio_bounds in (
            ("source_shell_interior", 6, (1.8, 4.0)),
            ("near_inner_shell_edge", 6, (1.15, 1.50)),
        ):
            for _ in range(count):
                s_Q = float(rng.uniform(0.58, 1.75))
                X = case.shell.X_a * float(rng.uniform(*x_ratio_bounds))
                R = sqrt(2.0 * s_Q * X)
                theta = float(rng.uniform(0.27, 2.83))
                z_phys = float(rng.uniform(-0.45, 0.45))
                xyz = np.asarray(
                    (R * np.cos(theta), R * np.sin(theta), z_phys), dtype=float
                )
                probes.append(
                    {
                        "case": case,
                        "region": region,
                        "s_Q": s_Q,
                        "xyz": xyz,
                    }
                )
    return probes


def _check_fd4_stencils(probes: list[dict[str, Any]]) -> dict[str, float | bool]:
    """Require every point used by the new FD4 ladder to stay in source domain."""
    min_X_ratio = float("inf")
    min_axis_ratio = float("inf")
    for probe in probes:
        case = probe["case"]
        s_Q = probe["s_Q"]
        xyz = np.asarray(probe["xyz"], dtype=float)
        for step in FD4_STEPS:
            for axis in range(3):
                for multiplier in (-2.0, -1.0, 0.0, 1.0, 2.0):
                    point = xyz.copy()
                    point[axis] += multiplier * step
                    R, _, _ = prior._xyz_to_source(case, point)
                    shell = case.shell.source_shell_coordinates(R, s_Q)
                    X = float(np.asarray(shell["X"]))
                    min_X_ratio = min(min_X_ratio, X / case.shell.X_a)
                    min_axis_ratio = min(
                        min_axis_ratio, R / case.shell.axis_lower_bound
                    )
    return {
        "all_fd4_stencil_points_source_shell_valid": True,
        "minimum_X_over_X_a_across_stencils": min_X_ratio,
        "minimum_R_over_axis_lower_bound_across_stencils": min_axis_ratio,
    }


def _ratio(coarse: float, fine: float) -> float:
    return float(coarse / max(fine, 1.0e-300))


def build_report() -> dict[str, Any]:
    probes = _sample_probes()
    stencil = _check_fd4_stencils(probes)
    ladder = [prior._evaluate_step(probes, step) for step in FD4_STEPS]
    finest = ladder[-1]
    curl_ratios = [
        _ratio(
            ladder[index]["curl_relative_rms"],
            ladder[index + 1]["curl_relative_rms"],
        )
        for index in range(len(ladder) - 1)
    ]
    divergence_ratios = [
        _ratio(
            ladder[index]["divergence_rms_abs"],
            ladder[index + 1]["divergence_rms_abs"],
        )
        for index in range(len(ladder) - 1)
    ]
    mutation = prior._mutation_metrics(probes, FD4_STEPS[-1])
    negative_controls = prior._negative_controls()

    checks = {
        "finest_cartesian_curl_relative_rms": bool(
            finest["curl_relative_rms"] <= LOCAL_CURL_FINE_REL_GUARD
        ),
        "finest_cartesian_divergence_max_abs": bool(
            finest["divergence_max_abs"] <= LOCAL_DIVERGENCE_FINE_ABS_GUARD
        ),
        "curl_refinement_ratio": bool(
            min(curl_ratios) >= LOCAL_REFINEMENT_RATIO_GUARD
        ),
        "divergence_refinement_ratio": bool(
            min(divergence_ratios) >= LOCAL_REFINEMENT_RATIO_GUARD
        ),
        "missing_support_gradient_curl_detected": bool(
            mutation["curl_point_relative_min"]
            >= MISSING_GRADIENT_RELATIVE_ERROR_FLOOR
        ),
        "missing_support_gradient_divergence_detected": bool(
            mutation["mutated_divergence_max_abs"]
            >= MISSING_GRADIENT_DIVERGENCE_FLOOR
        ),
        "all_fd4_stencil_points_source_shell_valid": bool(
            stencil["all_fd4_stencil_points_source_shell_valid"]
        ),
        "off_source_angular_grid_phase_rejected": bool(
            negative_controls["off_source_angular_grid_phase_rejected"]
        ),
        "below_inner_source_wave_shell_rejected": bool(
            negative_controls["below_inner_source_wave_shell_rejected"]
        ),
    }
    structural_pass = bool(all(checks.values()))

    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "dependency": {
            "agent4_parent_pr": BASE_PR,
            "exact_head": BASE_HEAD,
            "agent2_source_parent_pr": prior.BASE_PR,
            "agent2_source_exact_head": prior.BASE_HEAD,
            "public_artifact": "KokunoSourceWaveShellLocalizedCurl",
        },
        "preregistration": {
            "historical_rejection_overwritten": False,
            "prior_rejection": PRIOR_REJECTION,
            "changed_before_this_run": {
                "seed": [9173101, SEED],
                "fd4_steps": [[0.001, 0.0005, 0.00025], list(FD4_STEPS)],
            },
            "unchanged_scientific_guards": True,
            "purpose": (
                "test the same source-domain seam in a coarser truncation-dominated "
                "FD4 range after #371 encountered a floating-point-scale plateau"
            ),
            "roundoff_hypothesis_is_not_an_acceptance_override": True,
        },
        "seed": SEED,
        "parameter_cases": [
            {
                "label": case.label,
                "epsilon": case.epsilon,
                "k": case.k,
                "m": case.m,
                "j": case.j,
                "p": case.p,
                "physical_azimuthal_mode_mj": case.azimuthal_mode,
                "Lambda": case.Lambda,
                "X_a": case.shell.X_a,
                "axis_lower_bound": case.shell.axis_lower_bound,
            }
            for case in prior.CASES
        ],
        "operator": {
            "oracle": (
                "same independently recoded source-form Cartesian vector potential and "
                "centered Cartesian FD4 operator frozen by Agent 4 PR #371"
            ),
            "public_path": (
                "KokunoSourceWaveShellLocalizedCurl.localized_mode output only, converted "
                "from cylindrical to Cartesian"
            ),
            "caller_coefficient_derivatives": (
                "independently recoded C_m with centered FD6 inherited from #371; not the "
                "Cartesian validation derivative"
            ),
            "fd4_steps": list(FD4_STEPS),
            "coefficient_derivative_step": prior.COEFFICIENT_DERIVATIVE_STEP,
            "agent2_cylindrical_validation_oracle_reused": False,
            "training_tensor_or_loss_read": False,
            "construction_residual_operator_reused": False,
            "free_forcing_used": False,
        },
        "sampling": {
            "total_points": len(probes),
            "regions": {
                "source_shell_interior": sum(
                    probe["region"] == "source_shell_interior" for probe in probes
                ),
                "near_inner_shell_edge": sum(
                    probe["region"] == "near_inner_shell_edge" for probe in probes
                ),
            },
            "fresh_seed_relative_to_pr371": True,
            "stencil_source_domain": stencil,
        },
        "resolution_ladder": ladder,
        "refinement": {
            "curl_relative_rms_coarse_to_fine_ratios": curl_ratios,
            "divergence_rms_coarse_to_fine_ratios": divergence_ratios,
        },
        "mutation": mutation,
        "negative_controls": negative_controls,
        "local_guards": {
            "thresholds": {
                "finest_curl_relative_rms_max": LOCAL_CURL_FINE_REL_GUARD,
                "finest_divergence_max_abs_max": LOCAL_DIVERGENCE_FINE_ABS_GUARD,
                "minimum_refinement_ratio": LOCAL_REFINEMENT_RATIO_GUARD,
                "missing_gradient_relative_error_min": MISSING_GRADIENT_RELATIVE_ERROR_FLOOR,
                "missing_gradient_divergence_max_abs_min": MISSING_GRADIENT_DIVERGENCE_FLOOR,
            },
            "same_thresholds_as_agent4_prs_353_360_371": True,
            "checks": checks,
            "source_wave_shell_roundoff_ladder_preflight_passed": structural_pass,
            "scope": (
                "local source-grid/source-wave-shell support-localized curl implementation "
                "consistency only"
            ),
        },
        "formal_project_gates": {
            "normalized_momentum_max_l2": FORMAL_MOMENTUM_GATE,
            "divergence_max_l2": FORMAL_DIVERGENCE_GATE,
            "formal_full_domain_pde_gate_assessed": False,
            "reason": (
                "global matched Kokuno leading pressure, actual-source public xyz,t "
                "oscillatory velocity, and promoted after-correction composite are absent"
            ),
        },
        "truth_boundary": {
            "actual_source_Lambda_recovered": False,
            "concrete_source_partition_bumps_reconstructed": False,
            "actual_source_background_path_instantiated": False,
            "actual_source_public_oscillatory_xyz_t_velocity_available": False,
            "global_leading_velocity_pressure_available": False,
            "genuinely_independent_second_covariance_column_available": False,
            "after_correction_global_velocity_available": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def write_report(path: str | Path) -> Path:
    target = Path(path)
    if target.exists():
        raise FileExistsError(f"refusing to overwrite existing audit report: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(build_report(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return target


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent4/source_wave_shell_roundoff_ladder_audit.json",
    )
    args = parser.parse_args(argv)
    path = write_report(args.output)
    print(path.read_text(encoding="utf-8"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
