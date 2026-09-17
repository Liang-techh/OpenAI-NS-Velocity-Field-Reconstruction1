"""Cross-locate F(2,0) swirl sensitivity with baseline theta-momentum error.

This is a target-free routing diagnostic.  It does not fit a coefficient and it
never evaluates the momentum residual of a perturbed candidate; Agent 3 owns
that validation lane.  The only PDE quantity consumed here is the frozen
baseline candidate's pressure-independent cylindrical-theta residual, used as
an independent spatial error map for deciding whether the already-screened
F(2,0) direction is worth promoting to a bounded candidate family.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from .constrained_eq45_bipolar_f20_capacity import F10, F20, _embedded_field, _with_swirl_delta
from .constrained_force import RestrictedForce
from .constrained_validation import residual

TASK_ID = "CR003-BIPOLAR-F20-THETA-OVERLAP-036"
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "candidate_artifact_changed": False,
    "coefficient_fitted": False,
    "perturbed_candidate_pde_residual_evaluated": False,
    "force_or_pressure_changed": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
}


def _rings(radii, z_values=(-0.3, 0.0, 0.3), n_theta: int = 12) -> np.ndarray:
    theta = 2.0 * np.pi * np.arange(n_theta) / n_theta + np.pi / (2.0 * n_theta)
    rows = []
    for radius in radii:
        for z in z_values:
            rows.append(
                np.column_stack(
                    (
                        radius * np.cos(theta),
                        radius * np.sin(theta),
                        np.full(n_theta, z, dtype=float),
                    )
                )
            )
    return np.concatenate(rows)


def _theta_unit(points: np.ndarray) -> np.ndarray:
    r = np.hypot(points[:, 0], points[:, 1])
    if np.any(r <= 1e-14):
        raise ValueError("theta direction undefined on the axis")
    return np.column_stack((-points[:, 1] / r, points[:, 0] / r, np.zeros(len(points))))


def _swirl_response(field, mode, points, *, time: float, coefficient_step: float) -> tuple[np.ndarray, float]:
    plus = _with_swirl_delta(field, mode, coefficient_step)
    minus = _with_swirl_delta(field, mode, -coefficient_step)
    response = (plus.at_points(points, time) - minus.at_points(points, time)) / (2.0 * coefficient_step)
    theta = _theta_unit(points)
    swirl = np.sum(response * theta, axis=1)
    total_rms = float(np.sqrt(np.mean(np.sum(response * response, axis=1))))
    swirl_rms = float(np.sqrt(np.mean(swirl * swirl)))
    return swirl, float(swirl_rms / max(total_rms, 1e-300))


def _theta_residual(field, force, points, *, time: float, step: float, nu: float) -> np.ndarray:
    zero_pressure = lambda x, t: np.zeros(len(x), dtype=float)
    out = residual(
        field.at_points,
        zero_pressure,
        force,
        points,
        time,
        nu=nu,
        step=step,
        time_bounds=(0.25, 0.75),
    )
    return np.sum(out["momentum"] * _theta_unit(points), axis=1)


def _rms(values) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def _cosine_nonnegative(a, b) -> float:
    a = np.abs(np.asarray(a, dtype=float))
    b = np.abs(np.asarray(b, dtype=float))
    return float(np.dot(a, b) / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-300))


def _overlap(error: np.ndarray, response: np.ndarray) -> dict[str, float]:
    err = np.abs(np.asarray(error, dtype=float))
    weight = np.asarray(response, dtype=float) ** 2
    weighted_error = float(np.sum(weight * err) / max(np.sum(weight), 1e-300))
    mean_error = float(np.mean(err))
    cutoff = float(np.quantile(err, 0.75))
    high = err >= cutoff
    top_quartile_energy_share = float(np.sum(weight[high]) / max(np.sum(weight), 1e-300))
    return {
        "abs_theta_error_mean": mean_error,
        "response_energy_weighted_abs_theta_error": weighted_error,
        "error_enrichment": weighted_error / max(mean_error, 1e-300),
        "abs_response_abs_error_cosine": _cosine_nonnegative(response, error),
        "top_quartile_error_response_energy_share": top_quartile_energy_share,
    }


def audit_f20_theta_overlap(
    *,
    time: float = 0.5,
    coefficient_step: float = 0.01,
    derivative_step: float = 0.005,
    config_path: str | Path = "configs/constraints.json",
    training_report_path: str | Path = "artifacts/bipolar_joint_capped/report.json",
) -> dict[str, Any]:
    cfg = json.loads(Path(config_path).read_text())
    training = json.loads(Path(training_report_path).read_text())
    field = _embedded_field()
    force_parameters = training["parameters"][-2:]
    force = RestrictedForce(a=float(force_parameters[0]), c=float(force_parameters[1]))

    regions = {
        "interior": _rings((0.45, 0.75, 1.05)),
        "radial_flank": _rings((1.62, 1.72, 1.82)),
    }
    points = np.concatenate(tuple(regions.values()))
    labels = np.concatenate(
        [np.full(len(values), name, dtype=object) for name, values in regions.items()]
    )

    theta_error = _theta_residual(
        field,
        force,
        points,
        time=float(time),
        step=float(derivative_step),
        nu=float(cfg["nu"]),
    )

    modes = {}
    for mode, name in ((F10, "F10"), (F20, "F20")):
        swirl_response, swirl_fraction = _swirl_response(
            field,
            mode,
            points,
            time=float(time),
            coefficient_step=float(coefficient_step),
        )
        region_rows = {}
        for region_name in regions:
            mask = labels == region_name
            region_rows[region_name] = {
                "theta_error_rms": _rms(theta_error[mask]),
                "swirl_response_rms": _rms(swirl_response[mask]),
                **_overlap(theta_error[mask], swirl_response[mask]),
            }
        region_rows["radial_flank_to_interior_response_rms"] = (
            region_rows["radial_flank"]["swirl_response_rms"]
            / max(region_rows["interior"]["swirl_response_rms"], 1e-300)
        )
        modes[name] = {
            "swirl_fraction": swirl_fraction,
            "global_overlap": _overlap(theta_error, swirl_response),
            "regions": region_rows,
        }

    f10_selectivity = modes["F10"]["regions"]["radial_flank_to_interior_response_rms"]
    f20_selectivity = modes["F20"]["regions"]["radial_flank_to_interior_response_rms"]
    f10_enrichment = modes["F10"]["global_overlap"]["error_enrichment"]
    f20_enrichment = modes["F20"]["global_overlap"]["error_enrichment"]
    return {
        "task_id": TASK_ID,
        "source_candidate_sha256": field.sha256,
        "time": float(time),
        "coefficient_step": float(coefficient_step),
        "derivative_step": float(derivative_step),
        "sample_count": int(len(points)),
        "sampling": {
            "interior_radii": [0.45, 0.75, 1.05],
            "radial_flank_radii": [1.62, 1.72, 1.82],
            "z_values": [-0.3, 0.0, 0.3],
            "angles_per_ring": 12,
        },
        "baseline_theta_error": {
            "rms_all": _rms(theta_error),
            "rms_interior": _rms(theta_error[labels == "interior"]),
            "rms_radial_flank": _rms(theta_error[labels == "radial_flank"]),
            "radial_flank_to_interior_rms": _rms(theta_error[labels == "radial_flank"])
            / max(_rms(theta_error[labels == "interior"]), 1e-300),
        },
        "modes": modes,
        "F20_vs_F10_radial_flank_selectivity_ratio": f20_selectivity / max(f10_selectivity, 1e-300),
        "F20_vs_F10_error_enrichment_ratio": f20_enrichment / max(f10_enrichment, 1e-300),
        "interpretation_contract": (
            "Overlap is a routing diagnostic only. It asks whether a public-velocity sensitivity "
            "direction places response energy where the frozen baseline has large pressure-independent "
            "theta momentum error. It does not show that moving along that direction reduces residual."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="artifacts/bipolar_f20_theta_overlap/report.json")
    args = parser.parse_args()
    report = audit_f20_theta_overlap()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
