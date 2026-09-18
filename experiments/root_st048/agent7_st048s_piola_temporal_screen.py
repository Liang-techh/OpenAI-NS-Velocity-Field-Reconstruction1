"""Screen one bounded temporal coefficient on the frozen ST048-S Piola warp.

This Constrained Agent 7 increment adds no new spatial basis.  It starts from the
fixed beta=.075 diagnostic crossing from PR #407 and asks whether one even-in-time
coefficient can improve endpoint axial morphology without radial/collar inflation:

    beta_gamma(t) = .075 + gamma * (4 t - 2)^2,
    gamma in {0,.025,.05,.075,.125}.

The contravariant Piola pullback is unchanged, so spatial divergence preservation
holds at each frozen time.  A time-dependent beta changes u_t, however, so parent
pressure/forcing and held-out momentum evidence are deliberately not transferred.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import numpy as np

import agent7_st048s_piola_warp_screen as fixed

TASK_ID = "CR003-ST048S-PIOLA-TEMPORAL-COEFFICIENT-062"
PARENT_ID = fixed.PARENT_ID
BASE_BETA = 0.075
DEFAULT_GAMMAS = (0.0, 0.025, 0.05, 0.075, 0.125)
DEFAULT_TIMES = fixed.DEFAULT_TIMES
DEFAULT_ENERGY_TIMES = fixed.DEFAULT_ENERGY_TIMES
DEFAULT_GRID_SIZE = fixed.DEFAULT_GRID_SIZE
MORPHOLOGY_ZERO_TOLERANCE = fixed.MORPHOLOGY_ZERO_TOLERANCE
DIVERGENCE_PREFLIGHT_MAX = fixed.DIVERGENCE_PREFLIGHT_MAX
TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_spatial_basis_added": False,
    "one_temporal_geometry_coefficient_screened": True,
    "production_gamma_selected": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "held_out_pde_residual_transferred_from_parent": False,
    "public_image_used": False,
    "visual_score_defined": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
}


def _validate_gamma(gamma: float) -> float:
    gamma = float(gamma)
    if not np.isfinite(gamma) or not 0.0 <= gamma <= 0.125:
        raise ValueError("gamma must lie in the preregistered [0,.125] screen")
    return gamma


def _validate_gammas(values: Iterable[float]) -> tuple[float, ...]:
    values = tuple(_validate_gamma(value) for value in values)
    if len(values) < 2 or values[0] != 0.0:
        raise ValueError("gamma grid must start at zero and contain at least two values")
    if any(a >= b for a, b in zip(values, values[1:])):
        raise ValueError("gamma grid must be strictly increasing")
    return values


def beta_at_time(time: float, gamma: float, *, base_beta: float = BASE_BETA) -> float:
    time = float(time)
    gamma = _validate_gamma(gamma)
    if not np.isfinite(time) or not 0.25 <= time <= 0.75:
        raise ValueError("time must lie in [0.25,0.75]")
    beta = float(base_beta) + gamma * (4.0 * time - 2.0) ** 2
    return fixed._validate_beta(beta)


class TemporalPiolaWarpField:
    def __init__(self, parent, gamma: float, *, base_beta: float = BASE_BETA):
        self.parent = parent
        self.gamma = _validate_gamma(gamma)
        self.base_beta = fixed._validate_beta(base_beta)
        if beta_at_time(0.25, self.gamma, base_beta=self.base_beta) > 0.20 + 1.0e-15:
            raise ValueError("temporal schedule exceeds the already-screened beta range")

    def at_points(self, points, time: float):
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
            raise ValueError("points must be finite shape (n,3)")
        beta = beta_at_time(float(time), self.gamma, base_beta=self.base_beta)
        mapped_z, jac = fixed._warp_z_and_jacobian(points[:, 2], beta)
        mapped = np.array(points, copy=True)
        mapped[:, 2] = mapped_z
        velocity = np.asarray(self.parent.at_points(mapped, float(time)), dtype=float)
        if velocity.shape != points.shape or not np.all(np.isfinite(velocity)):
            raise RuntimeError("unexpected parent velocity")
        out = np.array(velocity, copy=True)
        out[:, 0] *= jac
        out[:, 1] *= jac
        return out


def _scaled_temporal_field(parent, gamma: float, *, energy_order: int):
    raw = TemporalPiolaWarpField(parent, gamma)
    energy = fixed._energy(raw, time=0.25, order=energy_order)
    scale = fixed._normalization_scale(energy)
    return raw, fixed._ScaledField(raw, scale), scale, energy


def _relative(child: float, parent: float) -> float:
    return float(child / max(abs(parent), 1.0e-300) - 1.0)


def _structure_checks(field, *, gamma: float, times: tuple[float, ...]):
    exterior = np.array(
        [[2.01, 0.0, 0.0], [-2.01, 0.0, 0.0], [0.4, 0.2, 2.01], [0.4, -0.2, -2.01]],
        dtype=float,
    )
    support_max = max(float(np.max(np.abs(field.at_points(exterior, time)))) for time in times)

    sign_checks = []
    divergence_max = 0.0
    fd_points = np.array(
        [[0.31, 0.22, 0.17], [0.62, -0.27, -0.41], [0.88, 0.31, 0.67], [1.22, -0.36, -1.22]],
        dtype=float,
    )
    h = 1.0e-5
    coordinate_rows = []
    for time in times:
        for radius in (0.05, 0.10, 0.20):
            for z in (-0.10, 0.10):
                ux, uy, uz = field.at_points(np.array([[radius, 0.0, z]], dtype=float), time)[0]
                sign_checks.append(bool(ux < 0.0 and uy > 0.0 and z * uz > 0.0))

        divergence = np.zeros(len(fd_points), dtype=float)
        for axis in range(3):
            offset = np.zeros_like(fd_points)
            offset[:, axis] = h
            plus = field.at_points(fd_points + offset, time)
            minus = field.at_points(fd_points - offset, time)
            divergence += (plus[:, axis] - minus[:, axis]) / (2.0 * h)
        divergence_max = max(divergence_max, float(np.max(np.abs(divergence))))

        beta = beta_at_time(time, gamma)
        z_grid = np.linspace(-fixed.PHYSICAL_Z_CUT, fixed.PHYSICAL_Z_CUT, 4001)
        mapped, jac = fixed._warp_z_and_jacobian(z_grid, beta)
        coordinate_rows.append({
            "time": time,
            "beta": beta,
            "minimum_coordinate_jacobian": float(np.min(jac)),
            "maximum_coordinate_jacobian": float(np.max(jac)),
            "map_monotone_on_support": bool(np.all(np.diff(mapped) > 0.0)),
        })

    return {
        "outside_support_max_abs_velocity": support_max,
        "representative_inward_swirl_bipolar_signs_all_pass": bool(all(sign_checks)),
        "representative_sign_checks": len(sign_checks),
        "cartesian_fd_divergence_max_abs": divergence_max,
        "coordinate_rows": coordinate_rows,
        "all_coordinate_maps_monotone": bool(all(row["map_monotone_on_support"] for row in coordinate_rows)),
    }


def _normalized_field_for_parameters(parent, *, base_beta: float, gamma: float, energy_order: int):
    class _LocalTemporalField:
        def at_points(self, points, time: float):
            points_array = np.asarray(points, dtype=float)
            beta = float(base_beta) + float(gamma) * (4.0 * float(time) - 2.0) ** 2
            beta = fixed._validate_beta(beta)
            mapped_z, jac = fixed._warp_z_and_jacobian(points_array[:, 2], beta)
            mapped = np.array(points_array, copy=True)
            mapped[:, 2] = mapped_z
            velocity = np.asarray(parent.at_points(mapped, float(time)), dtype=float)
            out = np.array(velocity, copy=True)
            out[:, 0] *= jac
            out[:, 1] *= jac
            return out

    raw = _LocalTemporalField()
    scale = fixed._normalization_scale(fixed._energy(raw, time=0.25, order=energy_order))
    return fixed._ScaledField(raw, scale)


def _two_column_sensitivity(parent, *, energy_order: int = 32):
    rng = np.random.default_rng(9174161)
    points = rng.uniform(-1.35, 1.35, size=(36, 3))
    times = np.resize(np.array([0.25, 0.375, 0.50, 0.625, 0.75], dtype=float), len(points))
    eps = 1.0e-4

    def values(base_beta: float, gamma: float):
        field = _normalized_field_for_parameters(
            parent, base_beta=base_beta, gamma=gamma, energy_order=energy_order
        )
        rows = []
        for point, time in zip(points, times):
            rows.append(field.at_points(point[None, :], float(time))[0])
        return np.asarray(rows, dtype=float).ravel()

    beta_plus = values(BASE_BETA + eps, 0.0)
    beta_minus = values(BASE_BETA - eps, 0.0)
    constant_column = (beta_plus - beta_minus) / (2.0 * eps)
    gamma_zero = values(BASE_BETA, 0.0)
    gamma_plus = values(BASE_BETA, eps)
    temporal_column = (gamma_plus - gamma_zero) / eps

    norms = np.array([np.linalg.norm(constant_column), np.linalg.norm(temporal_column)], dtype=float)
    if np.any(norms <= 1.0e-14):
        raise RuntimeError("degenerate sensitivity column")
    matrix = np.column_stack((constant_column / norms[0], temporal_column / norms[1]))
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix))
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0.0 else float("inf")
    cosine = float(np.dot(matrix[:, 0], matrix[:, 1]))
    return {
        "parameter_count": 2,
        "parameters": ["constant_beta", "endpoint_even_time_gamma"],
        "sample_count": len(points),
        "finite_difference_step": eps,
        "normalized_rank": rank,
        "normalized_singular_values": singular.tolist(),
        "normalized_condition_number": condition,
        "normalized_column_cosine": cosine,
        "raw_column_norms": norms.tolist(),
    }


def audit_temporal_piola_screen(
    *,
    gammas: Iterable[float] = DEFAULT_GAMMAS,
    times: Iterable[float] = DEFAULT_TIMES,
    energy_times: Iterable[float] = DEFAULT_ENERGY_TIMES,
    energy_orders: tuple[int, int] = (48, 64),
    grid_size: int = DEFAULT_GRID_SIZE,
):
    gammas = _validate_gammas(gammas)
    times = tuple(float(time) for time in times)
    energy_times = tuple(float(time) for time in energy_times)
    if times != tuple(DEFAULT_TIMES):
        raise ValueError("this preregistered screen uses t=.25/.50/.75")
    if not energy_times or any(time < 0.25 or time > 0.75 for time in energy_times):
        raise ValueError("energy times must lie in [0.25,0.75]")
    if len(energy_orders) != 2 or not 12 <= energy_orders[0] < energy_orders[1]:
        raise ValueError("energy_orders must be two increasing orders >=12")
    if grid_size < 17 or grid_size % 2 == 0:
        raise ValueError("grid_size must be odd and >=17")

    _ignored_parent, children = fixed.transfer._load_frozen_fields()
    parent = children[PARENT_ID]

    raw_energy = {}
    refinement_max = 0.0
    for gamma in gammas:
        raw = TemporalPiolaWarpField(parent, gamma)
        raw_energy[gamma] = {}
        for time in energy_times:
            coarse = fixed._energy(raw, time=time, order=energy_orders[0])
            fine = fixed._energy(raw, time=time, order=energy_orders[1])
            raw_energy[gamma][time] = fine
            refinement_max = max(refinement_max, abs(coarse - fine) / max(abs(fine), 1.0e-12))

    baseline_raw = TemporalPiolaWarpField(parent, 0.0)
    baseline_scale = fixed._normalization_scale(raw_energy[0.0][0.25])
    baseline = fixed._ScaledField(baseline_raw, baseline_scale)
    baseline_metrics = {
        time: fixed.transfer._morphology(baseline, time=time, grid_size=grid_size)
        for time in times
    }

    rows = []
    crossings = []
    endpoint_times = (times[0], times[-1])
    for gamma in gammas:
        raw = TemporalPiolaWarpField(parent, gamma)
        scale = fixed._normalization_scale(raw_energy[gamma][0.25])
        field = fixed._ScaledField(raw, scale)
        normalized_energy = {
            time: scale * scale * raw_energy[gamma][time]
            for time in energy_times
        }
        reference_energy_pass = abs(normalized_energy[0.25] - 1.0) <= fixed.REFERENCE_ENERGY_TOLERANCE
        energy_range_pass = all(
            fixed.VALIDATION_ENERGY_RANGE[0] <= value <= fixed.VALIDATION_ENERGY_RANGE[1]
            for value in normalized_energy.values()
        )
        metrics = {
            time: fixed.transfer._morphology(field, time=time, grid_size=grid_size)
            for time in times
        }
        axial = [
            _relative(metrics[time]["axial_rms_over_support"], baseline_metrics[time]["axial_rms_over_support"])
            for time in times
        ]
        radial = [
            _relative(metrics[time]["radial_rms_over_support"], baseline_metrics[time]["radial_rms_over_support"])
            for time in times
        ]
        q90 = [
            float(metrics[time]["axial_q90_over_support"] - baseline_metrics[time]["axial_q90_over_support"])
            for time in times
        ]
        q99 = [
            float(metrics[time]["axial_q99_over_support"] - baseline_metrics[time]["axial_q99_over_support"])
            for time in times
        ]
        collar = [
            metrics[time]["outer_support_collar_enstrophy_fraction"]
            / max(baseline_metrics[time]["outer_support_collar_enstrophy_fraction"], 1.0e-300)
            for time in times
        ]
        tail065 = [
            _relative(metrics[time]["outer_065_enstrophy_fraction"], baseline_metrics[time]["outer_065_enstrophy_fraction"])
            for time in times
        ]
        tail075 = [
            _relative(metrics[time]["outer_075_enstrophy_fraction"], baseline_metrics[time]["outer_075_enstrophy_fraction"])
            for time in times
        ]
        structure = _structure_checks(field, gamma=gamma, times=times)
        preflight = bool(
            reference_energy_pass
            and energy_range_pass
            and structure["outside_support_max_abs_velocity"] <= 1.0e-12
            and structure["representative_inward_swirl_bipolar_signs_all_pass"]
            and structure["cartesian_fd_divergence_max_abs"] <= DIVERGENCE_PREFLIGHT_MAX
            and structure["all_coordinate_maps_monotone"]
        )
        endpoint_indices = (0, len(times) - 1)
        extension = bool(
            all(value >= -MORPHOLOGY_ZERO_TOLERANCE for value in axial)
            and any(
                q90[index] > MORPHOLOGY_ZERO_TOLERANCE or q99[index] > MORPHOLOGY_ZERO_TOLERANCE
                for index in endpoint_indices
            )
        )
        low_cost = bool(max(radial) <= 0.05 and max(collar) <= 2.0)
        clean = bool(gamma > 0.0 and preflight and extension and low_cost)
        if clean:
            crossings.append(gamma)
        rows.append({
            "gamma": gamma,
            "beta_by_time": {str(time): beta_at_time(time, gamma) for time in times},
            "common_energy_normalization_scale": scale,
            "normalized_energy_by_time": {str(time): value for time, value in normalized_energy.items()},
            "reference_energy_pass": reference_energy_pass,
            "validation_energy_range_pass": energy_range_pass,
            "preflight_pass": preflight,
            "temporal_endpoint_extension_rule_pass": clean,
            "axial_rms_relative_changes_vs_fixed_beta075": axial,
            "radial_rms_relative_changes_vs_fixed_beta075": radial,
            "q90_deltas_over_support_vs_fixed_beta075": q90,
            "q99_deltas_over_support_vs_fixed_beta075": q99,
            "collar_enstrophy_ratios_vs_fixed_beta075": collar,
            "outer_065_enstrophy_relative_changes_vs_fixed_beta075": tail065,
            "outer_075_enstrophy_relative_changes_vs_fixed_beta075": tail075,
            "structure": structure,
        })

    sensitivity = _two_column_sensitivity(parent)
    smallest = min(crossings) if crossings else None
    routing = (
        "A clean endpoint-morphology crossing exists with only one bounded temporal geometry coefficient; "
        "do not add another spatial basis. Any governed child must treat gamma as a new time-dependent degree, "
        "rebuild pressure/restricted forcing, and rerun fresh held-out full momentum because u_t changes."
        if smallest is not None else
        "No preregistered temporal coefficient cleanly moves endpoint q90/q99 relative to the fixed beta=.075 field. "
        "Do not grow the spatial basis from this negative result; retain fixed beta=.075 as the simpler geometry route "
        "and route the remaining mismatch to actual 3-D render/trajectory evidence before another representation change."
    )
    return {
        "task_id": TASK_ID,
        "issue": 416,
        "parent_id": PARENT_ID,
        "base_beta": BASE_BETA,
        "gamma_schedule": "beta(t)=0.075+gamma*(4*t-2)^2",
        "preregistered_gammas": list(gammas),
        "morphology_times": list(times),
        "energy_times": list(energy_times),
        "grid_size": grid_size,
        "smallest_preregistered_temporal_endpoint_crossing_gamma": smallest,
        "energy_quadrature_orders": list(energy_orders),
        "energy_quadrature_max_relative_refinement_change": refinement_max,
        "two_column_velocity_sensitivity": sensitivity,
        "rows": rows,
        "routing": routing,
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    report = audit_temporal_piola_screen()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "task_id": report["task_id"],
        "smallest_crossing_gamma": report["smallest_preregistered_temporal_endpoint_crossing_gamma"],
        "sensitivity": report["two_column_velocity_sensitivity"],
        "routing": report["routing"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
