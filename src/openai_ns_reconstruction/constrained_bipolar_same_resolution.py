"""Same-resolution CR005 replay for the bipolar six-mode and odd-Phi(0,3) families."""
from __future__ import annotations

import argparse
import contextlib
import io
import json
from pathlib import Path

import numpy as np

from .constrained_bipolar_joint_energy import run as run_joint_energy
from .eq45_supported_delivery import Eq45SupportedDeliveryField


DEFAULT_ODD_REPORT = Path("artifacts/bipolar_joint_odd3/report.json")
DEFAULT_OUTPUT = Path("artifacts/bipolar_joint_same_resolution")


def _best_physical_objective(report: dict) -> float:
    values = [
        float(row["objective"])
        for row in report["multistart_trials"]
        if row["physical_constraints_satisfied"]
    ]
    if not values:
        raise RuntimeError("report contains no energy/core-feasible multistart result")
    return min(values)


def _summary(report: dict) -> dict:
    holdout = np.asarray([row["balance_defect"] for row in report["holdout"]], dtype=float)
    modes = report["modes"]
    parameters = np.asarray(report["parameters"], dtype=float)
    velocity = parameters[: len(modes)]
    force = parameters[len(modes) :]
    if force.shape != (2,):
        raise RuntimeError("expected exactly the preregistered two force coefficients")
    training_defects = np.asarray(report["training_equality_defects"][1:], dtype=float)
    return {
        "candidate_sha256": report["candidate_sha256"],
        "training_quadrature_order": int(report["training_quadrature_order"]),
        "mode_count": len(modes),
        "best_physical_training_objective": _best_physical_objective(report),
        "training_balance_rms": float(np.sqrt(np.mean(training_defects**2))),
        "training_balance_max_abs": float(np.max(np.abs(training_defects))),
        "holdout_balance_rms": float(np.sqrt(np.mean(holdout**2))),
        "holdout_balance_max_abs": float(np.max(np.abs(holdout))),
        "holdout_balance_defects": holdout.tolist(),
        "core_ratios": [float(value) for value in report["core_ratios"]],
        "velocity_parameters": velocity.tolist(),
        "velocity_parameters_at_bound": int(np.count_nonzero(np.abs(velocity) >= 4.0 - 1e-8)),
        "force_parameters": {"a": float(force[0]), "c": float(force[1])},
        "optimizer_success": bool(report["optimizer_success"]),
    }


def run(output: str | Path = DEFAULT_OUTPUT, odd_report: str | Path = DEFAULT_ODD_REPORT) -> dict:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    baseline_dir = output / "baseline_order96"

    # Silence the nested report dump; this comparator emits one compact deterministic receipt.
    with contextlib.redirect_stdout(io.StringIO()):
        run_joint_energy(
            multistart=True,
            odd_extension=False,
            training_order=96,
            output=baseline_dir,
        )

    baseline = json.loads((baseline_dir / "report.json").read_text())
    odd = json.loads(Path(odd_report).read_text())
    if baseline["odd_phi03_extension"]:
        raise RuntimeError("baseline replay unexpectedly enabled the odd extension")
    if not odd["odd_phi03_extension"]:
        raise RuntimeError("comparison report is not the odd-Phi(0,3) extension")
    if baseline["training_quadrature_order"] != 96 or odd["training_quadrature_order"] != 96:
        raise RuntimeError("same-resolution comparison requires order-96 training on both sides")

    baseline_summary = _summary(baseline)
    odd_summary = _summary(odd)
    field = Eq45SupportedDeliveryField.load_candidate(baseline_dir / "candidate.json")
    grid = field.grid(
        np.asarray([-0.5, 0.0, 0.5]),
        np.asarray([-0.5, 0.0, 0.5]),
        np.asarray([-0.5, 0.0, 0.5]),
        np.asarray([0.25, 0.5, 0.75]),
    )
    grid_finite = bool(np.all(np.isfinite(grid)))
    grid_nonzero = bool(np.max(np.abs(grid)) > 0.0)
    if not grid_finite or not grid_nonzero:
        raise RuntimeError("baseline candidate failed direct finite/nonzero [u,v,w] grid export")

    baseline_rms = baseline_summary["holdout_balance_rms"]
    odd_rms = odd_summary["holdout_balance_rms"]
    receipt = {
        "task": "CR005-BIPOLAR-ODD3-SAME-RESOLUTION-BASELINE-036",
        "comparison_contract": {
            "training_quadrature_order": 96,
            "optimizer_seed": 9172609,
            "multistart_count": 12,
            "velocity_bounds": [-4.0, 4.0],
            "restricted_force_bounds": {"a": [0.0, 10.0], "c": [0.0, 10.0]},
            "optimizer": "SLSQP",
            "max_iterations_per_start": 500,
            "ftol": 1e-11,
            "training_times": [0.3125, 0.5, 0.6875],
            "training_derivative_step": 0.005,
            "holdout_times": [0.34375, 0.46875, 0.59375, 0.71875],
            "holdout_quadrature_order": 96,
            "holdout_derivative_step": 0.0025,
            "initial_energy_equality": 1.0,
            "central_component_ratio_bounds": [0.95, 1.05],
            "force_family": "RestrictedForce(a,c) only; no residual-defined force",
        },
        "baseline_six_mode": baseline_summary,
        "odd_phi03_extension": odd_summary,
        "odd_vs_baseline": {
            "holdout_rms_fractional_change": float(odd_rms / baseline_rms - 1.0),
            "holdout_rms_reduction_fraction": float(1.0 - odd_rms / baseline_rms),
            "training_objective_fractional_change": float(
                odd_summary["best_physical_training_objective"]
                / baseline_summary["best_physical_training_objective"]
                - 1.0
            ),
        },
        "baseline_velocity_export": {
            "grid_shape": list(grid.shape),
            "finite": grid_finite,
            "nonzero": grid_nonzero,
        },
        "interpretation_guard": {
            "optimizer_convergence_is_pde_validation": False,
            "energy_balance_is_full_momentum_validation": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "candidate_selection_resolved": False,
        },
    }
    (output / "report.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--odd-report", default=str(DEFAULT_ODD_REPORT))
    args = parser.parse_args()
    run(args.output, args.odd_report)


if __name__ == "__main__":
    main()
