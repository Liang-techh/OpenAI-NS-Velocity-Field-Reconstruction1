"""Screen one bounded axisymmetric swirl-gain coordinate on the ST048-S temporal Piola field.

This Constrained Agent 7 increment follows Agent-9 PR #425: the gamma=.025
temporal Piola schedule improves Eulerian axial morphology but slightly weakens
material-path winding.  Rather than add another spatial basis or more axial
geometry, screen exactly one azimuthal capacity coordinate

    u_theta -> (1 + kappa) u_theta

before one common positive E(.25)=1 normalization.

For an axisymmetric field, a constant azimuthal gain preserves the support and
adds no cylindrical divergence term.  It still changes vorticity, nonlinear
momentum and energy normalization, so parent pressure/forcing/PDE evidence is
not transferred.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import numpy as np

import agent7_st048s_piola_temporal_screen as temporal
import agent7_st048s_piola_warp_screen as fixed

TASK_ID = "CR003-ST048S-TEMPORAL-PIOLA-SWIRL-GAIN-063"
PARENT_ID = fixed.PARENT_ID
BASE_BETA = temporal.BASE_BETA
BASE_GAMMA = 0.025
DEFAULT_KAPPAS = (0.0, 0.025, 0.05, 0.075, 0.10, 0.15)
DEFAULT_TIMES = temporal.DEFAULT_TIMES
DEFAULT_ENERGY_TIMES = temporal.DEFAULT_ENERGY_TIMES
DEFAULT_GRID_SIZE = temporal.DEFAULT_GRID_SIZE
WINDING_GAIN_FLOOR = 0.02
AXIAL_RMS_ABS_CHANGE_MAX = 0.02
RADIAL_RMS_GROWTH_MAX = 0.03
COLLAR_RATIO_MAX = 1.25

TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_spatial_basis_added": False,
    "one_swirl_gain_coordinate_screened": True,
    "production_kappa_selected": False,
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


def _validate_kappa(kappa: float) -> float:
    kappa = float(kappa)
    if not np.isfinite(kappa) or not 0.0 <= kappa <= 0.15:
        raise ValueError("kappa must lie in the preregistered [0,.15] screen")
    return kappa


def _validate_kappas(values: Iterable[float]) -> tuple[float, ...]:
    values = tuple(_validate_kappa(value) for value in values)
    if len(values) < 2 or values[0] != 0.0:
        raise ValueError("kappa grid must start at zero and contain at least two values")
    if any(a >= b for a, b in zip(values, values[1:])):
        raise ValueError("kappa grid must be strictly increasing")
    return values


class SwirlGainField:
    """Multiply only the cylindrical azimuthal component by 1+kappa."""

    def __init__(self, parent, kappa: float, *, validate: bool = True):
        self.parent = parent
        self.kappa = _validate_kappa(kappa) if validate else float(kappa)
        if not np.isfinite(self.kappa) or 1.0 + self.kappa <= 0.0:
            raise ValueError("swirl multiplier must stay positive")

    def at_points(self, points, time: float):
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
            raise ValueError("points must be finite shape (n,3)")
        velocity = np.asarray(self.parent.at_points(points, float(time)), dtype=float)
        if velocity.shape != points.shape or not np.all(np.isfinite(velocity)):
            raise RuntimeError("unexpected parent velocity")

        x = points[:, 0]
        y = points[:, 1]
        radius = np.hypot(x, y)
        out = np.array(velocity, copy=True)
        mask = radius > 1.0e-14
        if np.any(mask):
            rx = x[mask] / radius[mask]
            ry = y[mask] / radius[mask]
            ux = velocity[mask, 0]
            uy = velocity[mask, 1]
            radial = rx * ux + ry * uy
            azimuthal = -ry * ux + rx * uy
            azimuthal *= 1.0 + self.kappa
            out[mask, 0] = rx * radial - ry * azimuthal
            out[mask, 1] = ry * radial + rx * azimuthal
        return out


def _relative(child: float, parent: float) -> float:
    return float(child / max(abs(parent), 1.0e-300) - 1.0)


def _raw_field(parent, *, base_beta: float = BASE_BETA, gamma: float = BASE_GAMMA, kappa: float = 0.0, validate_kappa: bool = True):
    temporal_field = temporal.TemporalPiolaWarpField(parent, gamma, base_beta=base_beta)
    return SwirlGainField(temporal_field, kappa, validate=validate_kappa)


def _scaled_field(parent, *, kappa: float, energy_order: int):
    raw = _raw_field(parent, kappa=kappa)
    energy = fixed._energy(raw, time=0.25, order=energy_order)
    scale = fixed._normalization_scale(energy)
    return raw, fixed._ScaledField(raw, scale), scale, energy


def _probe_points() -> np.ndarray:
    rows = []
    for radius in (0.6, 0.9, 1.2):
        for z in (-0.3, 0.3):
            for theta in np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False):
                rows.append([radius * np.cos(theta), radius * np.sin(theta), z])
    return np.asarray(rows, dtype=float)


def _winding_proxies(field, *, time: float):
    points = _probe_points()
    velocity = np.asarray(field.at_points(points, float(time)), dtype=float)
    x = points[:, 0]
    y = points[:, 1]
    radius = np.hypot(x, y)
    rx = x / radius
    ry = y / radius
    radial = rx * velocity[:, 0] + ry * velocity[:, 1]
    azimuthal = -ry * velocity[:, 0] + rx * velocity[:, 1]
    axial = velocity[:, 2]
    angular_rate = np.abs(azimuthal) / radius
    poloidal = np.hypot(radial, axial)
    swirl_rms = float(np.sqrt(np.mean(azimuthal * azimuthal)))
    poloidal_rms = float(np.sqrt(np.mean(poloidal * poloidal)))
    pitch_numerator = float(np.sqrt(np.mean((radius * axial) ** 2)))
    return {
        "probe_count": int(len(points)),
        "mean_abs_angular_rate": float(np.mean(angular_rate)),
        "mean_inward_radial_speed": float(np.mean(-radial)),
        "swirl_to_poloidal_rms_ratio": swirl_rms / max(poloidal_rms, 1.0e-300),
        "helical_pitch_proxy_rms": pitch_numerator / max(swirl_rms, 1.0e-300),
    }


def _normalized_parameter_field(parent, *, base_beta: float, gamma: float, kappa: float, energy_order: int):
    raw = _raw_field(
        parent,
        base_beta=base_beta,
        gamma=gamma,
        kappa=kappa,
        validate_kappa=False,
    )
    scale = fixed._normalization_scale(fixed._energy(raw, time=0.25, order=energy_order))
    return fixed._ScaledField(raw, scale)


def _three_column_sensitivity(parent, *, energy_order: int = 32):
    rng = np.random.default_rng(9174271)
    points = rng.uniform(-1.35, 1.35, size=(36, 3))
    times = np.resize(np.array([0.25, 0.375, 0.50, 0.625, 0.75], dtype=float), len(points))
    eps = 1.0e-4

    def values(base_beta: float, gamma: float, kappa: float):
        field = _normalized_parameter_field(
            parent,
            base_beta=base_beta,
            gamma=gamma,
            kappa=kappa,
            energy_order=energy_order,
        )
        rows = [
            field.at_points(point[None, :], float(time))[0]
            for point, time in zip(points, times)
        ]
        return np.asarray(rows, dtype=float).ravel()

    beta_column = (
        values(BASE_BETA + eps, BASE_GAMMA, 0.0)
        - values(BASE_BETA - eps, BASE_GAMMA, 0.0)
    ) / (2.0 * eps)
    gamma_column = (
        values(BASE_BETA, BASE_GAMMA + eps, 0.0)
        - values(BASE_BETA, BASE_GAMMA - eps, 0.0)
    ) / (2.0 * eps)
    zero = values(BASE_BETA, BASE_GAMMA, 0.0)
    swirl_column = (values(BASE_BETA, BASE_GAMMA, eps) - zero) / eps

    columns = [beta_column, gamma_column, swirl_column]
    norms = np.asarray([np.linalg.norm(column) for column in columns], dtype=float)
    if np.any(norms <= 1.0e-14):
        raise RuntimeError("degenerate sensitivity column")
    matrix = np.column_stack([column / norm for column, norm in zip(columns, norms)])
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix))
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0.0 else float("inf")
    cosine = matrix.T @ matrix
    return {
        "parameter_count": 3,
        "parameters": ["constant_beta", "endpoint_even_time_gamma", "swirl_gain_kappa"],
        "sample_count": int(len(points)),
        "finite_difference_step": eps,
        "normalized_rank": rank,
        "normalized_singular_values": singular.tolist(),
        "normalized_condition_number": condition,
        "normalized_column_cosines": cosine.tolist(),
        "raw_column_norms": norms.tolist(),
    }


def audit_swirl_gain_screen(
    *,
    kappas: Iterable[float] = DEFAULT_KAPPAS,
    times: Iterable[float] = DEFAULT_TIMES,
    energy_times: Iterable[float] = DEFAULT_ENERGY_TIMES,
    energy_orders: tuple[int, int] = (48, 64),
    grid_size: int = DEFAULT_GRID_SIZE,
):
    kappas = _validate_kappas(kappas)
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
    for kappa in kappas:
        raw = _raw_field(parent, kappa=kappa)
        raw_energy[kappa] = {}
        for time in energy_times:
            coarse = fixed._energy(raw, time=time, order=energy_orders[0])
            fine = fixed._energy(raw, time=time, order=energy_orders[1])
            raw_energy[kappa][time] = fine
            refinement_max = max(
                refinement_max,
                abs(coarse - fine) / max(abs(fine), 1.0e-12),
            )

    base_raw = _raw_field(parent, kappa=0.0)
    base_scale = fixed._normalization_scale(raw_energy[0.0][0.25])
    baseline = fixed._ScaledField(base_raw, base_scale)
    baseline_morphology = {
        time: fixed.transfer._morphology(baseline, time=time, grid_size=grid_size)
        for time in times
    }
    baseline_winding = {time: _winding_proxies(baseline, time=time) for time in times}

    rows = []
    crossings = []
    for kappa in kappas:
        raw = _raw_field(parent, kappa=kappa)
        scale = fixed._normalization_scale(raw_energy[kappa][0.25])
        field = fixed._ScaledField(raw, scale)
        normalized_energy = {
            time: scale * scale * raw_energy[kappa][time]
            for time in energy_times
        }
        reference_energy_pass = abs(normalized_energy[0.25] - 1.0) <= fixed.REFERENCE_ENERGY_TOLERANCE
        energy_range_pass = all(
            fixed.VALIDATION_ENERGY_RANGE[0] <= value <= fixed.VALIDATION_ENERGY_RANGE[1]
            for value in normalized_energy.values()
        )
        morphology = {
            time: fixed.transfer._morphology(field, time=time, grid_size=grid_size)
            for time in times
        }
        winding = {time: _winding_proxies(field, time=time) for time in times}
        angular_changes = [
            _relative(winding[time]["mean_abs_angular_rate"], baseline_winding[time]["mean_abs_angular_rate"])
            for time in times
        ]
        swirl_ratio_changes = [
            _relative(winding[time]["swirl_to_poloidal_rms_ratio"], baseline_winding[time]["swirl_to_poloidal_rms_ratio"])
            for time in times
        ]
        pitch_changes = [
            _relative(winding[time]["helical_pitch_proxy_rms"], baseline_winding[time]["helical_pitch_proxy_rms"])
            for time in times
        ]
        axial = [
            _relative(morphology[time]["axial_rms_over_support"], baseline_morphology[time]["axial_rms_over_support"])
            for time in times
        ]
        radial = [
            _relative(morphology[time]["radial_rms_over_support"], baseline_morphology[time]["radial_rms_over_support"])
            for time in times
        ]
        q90 = [
            float(morphology[time]["axial_q90_over_support"] - baseline_morphology[time]["axial_q90_over_support"])
            for time in times
        ]
        q99 = [
            float(morphology[time]["axial_q99_over_support"] - baseline_morphology[time]["axial_q99_over_support"])
            for time in times
        ]
        collar = [
            morphology[time]["outer_support_collar_enstrophy_fraction"]
            / max(baseline_morphology[time]["outer_support_collar_enstrophy_fraction"], 1.0e-300)
            for time in times
        ]
        structure = temporal._structure_checks(field, gamma=BASE_GAMMA, times=times)
        preflight = bool(
            reference_energy_pass
            and energy_range_pass
            and structure["outside_support_max_abs_velocity"] <= 1.0e-12
            and structure["representative_inward_swirl_bipolar_signs_all_pass"]
            and structure["cartesian_fd_divergence_max_abs"] <= temporal.DIVERGENCE_PREFLIGHT_MAX
            and structure["all_coordinate_maps_monotone"]
        )
        clean = bool(
            kappa > 0.0
            and preflight
            and min(angular_changes) >= WINDING_GAIN_FLOOR - 1.0e-12
            and max(abs(value) for value in axial) <= AXIAL_RMS_ABS_CHANGE_MAX + 1.0e-12
            and max(radial) <= RADIAL_RMS_GROWTH_MAX + 1.0e-12
            and max(collar) <= COLLAR_RATIO_MAX + 1.0e-12
        )
        if clean:
            crossings.append(kappa)
        rows.append({
            "kappa": kappa,
            "azimuthal_multiplier": 1.0 + kappa,
            "common_energy_normalization_scale": scale,
            "normalized_energy": {str(key): value for key, value in normalized_energy.items()},
            "reference_energy_pass": reference_energy_pass,
            "energy_range_pass": energy_range_pass,
            "winding_proxy_by_time": {str(key): value for key, value in winding.items()},
            "mean_angular_rate_relative_changes_vs_kappa0": angular_changes,
            "swirl_to_poloidal_relative_changes_vs_kappa0": swirl_ratio_changes,
            "helical_pitch_proxy_relative_changes_vs_kappa0": pitch_changes,
            "axial_rms_relative_changes_vs_kappa0": axial,
            "radial_rms_relative_changes_vs_kappa0": radial,
            "q90_deltas_over_support_vs_kappa0": q90,
            "q99_deltas_over_support_vs_kappa0": q99,
            "collar_enstrophy_ratios_vs_kappa0": collar,
            "structure": structure,
            "preflight_pass": preflight,
            "clean_winding_capacity_rule_pass": clean,
        })

    sensitivity = _three_column_sensitivity(parent)
    report = {
        "task_id": TASK_ID,
        "parent_id": PARENT_ID,
        "fixed_temporal_geometry": {
            "base_beta": BASE_BETA,
            "gamma": BASE_GAMMA,
            "beta_by_time": {
                str(time): temporal.beta_at_time(time, BASE_GAMMA)
                for time in times
            },
        },
        "kappa_grid": list(kappas),
        "times": list(times),
        "energy_times": list(energy_times),
        "energy_orders": list(energy_orders),
        "grid_size": grid_size,
        "preregistered_clean_rule": {
            "minimum_mean_angular_rate_relative_gain_each_time": WINDING_GAIN_FLOOR,
            "maximum_absolute_axial_vorticity_rms_relative_change": AXIAL_RMS_ABS_CHANGE_MAX,
            "maximum_radial_vorticity_rms_growth": RADIAL_RMS_GROWTH_MAX,
            "maximum_support_collar_enstrophy_ratio": COLLAR_RATIO_MAX,
            "requires_energy_support_core_sign_divergence_preflight": True,
        },
        "energy_quadrature_max_relative_refinement_change": refinement_max,
        "baseline_winding_proxy_by_time": {
            str(key): value for key, value in baseline_winding.items()
        },
        "rows": rows,
        "smallest_preregistered_clean_winding_crossing_kappa": min(crossings) if crossings else None,
        "three_column_velocity_sensitivity": sensitivity,
        "routing": (
            "If a nonzero clean crossing exists, treat it only as evidence that one bounded "
            "azimuthal degree can repair Eulerian winding capacity without another spatial basis. "
            "Do not promote it from this screen. A governed child must choose an explicit kappa "
            "bound/value, keep the temporal Piola chain rule, rebuild compatible pressure/restricted "
            "forcing, and rerun fresh held-out full momentum plus independent material-path/3-D "
            "streamline comparison. If no crossing exists, do not widen kappa post hoc; route back "
            "to a localized swirl-shape degree rather than more axial geometry."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--grid-size", type=int, default=DEFAULT_GRID_SIZE)
    args = parser.parse_args()
    report = audit_swirl_gain_screen(grid_size=args.grid_size)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)


if __name__ == "__main__":
    main()
