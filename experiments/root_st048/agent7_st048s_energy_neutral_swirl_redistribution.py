"""Screen one first-order energy-neutral inner/mid swirl redistribution coordinate.

Constrained Agent 7 follows #437 and Agent-9 #445.  The frozen material-path
winding deficit is concentrated at seed radii r=.6 and r=.9, while r=1.2 is
already a surplus relative to the retained ST006 benchmark.  The previous
outer-peaked ring replacement was independent but inefficient because common
energy normalization weakened inward/poloidal motion before useful winding gain.

This increment keeps frozen ST048-S and the temporal Piola schedule

    beta(t) = .075 + .025 (4t-2)^2

and replaces global kappa with exactly one azimuthal redistribution coordinate

    u_theta -> (1 + kappa_redist * h(r)) u_theta,

where h = g_inner - alpha g_outer.  Both g functions are C-infinity compact
radial-window bumps.  The coefficient alpha is computed once from the frozen
candidate at t=.25 so that the *first derivative* of kinetic energy with respect
to kappa_redist is numerically zero before the final exact E(.25)=1 common
normalization:

    integral u_theta^2 h(r) r dr dz = 0.

This is candidate-side autonomous representation design informed by the frozen
trajectory attribution in #445.  ST006 is a benchmark, not OpenAI truth.  No
pressure/forcing fit, hidden OpenAI number, free f=R, u->0 route, or held-out PDE
transfer is used here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
from numpy.polynomial.legendre import leggauss

import agent7_st048s_piola_swirl_gain_screen as global_screen
import agent7_st048s_piola_temporal_screen as temporal
import agent7_st048s_piola_warp_screen as fixed

TASK_ID = "CR003-ST048S-ENERGY-NEUTRAL-SWIRL-REDISTRIBUTION-065"
PARENT_ID = global_screen.PARENT_ID
BASE_BETA = global_screen.BASE_BETA
BASE_GAMMA = global_screen.BASE_GAMMA
DEFAULT_GAINS = (0.0, 0.015, 0.025, 0.035, 0.05)
DEFAULT_TIMES = global_screen.DEFAULT_TIMES
DEFAULT_ENERGY_TIMES = global_screen.DEFAULT_ENERGY_TIMES
DEFAULT_GRID_SIZE = global_screen.DEFAULT_GRID_SIZE
SEED_RADII = (0.6, 0.9, 1.2)
INNER_WINDOW = (0.30, 1.05)
OUTER_WINDOW = (0.95, 1.85)
INNER_GAIN_FLOOR = 0.02
MID_GAIN_FLOOR = 0.01
OUTER_GAIN_CEILING = 0.005
NORMALIZATION_DEVIATION_MAX = 0.005
INWARD_SPEED_LOSS_MAX = 0.005
AXIAL_RMS_ABS_CHANGE_MAX = 0.02
RADIAL_RMS_GROWTH_MAX = 0.03
COLLAR_RATIO_MAX = 1.25

TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_spatial_basis_added": False,
    "one_swirl_redistribution_coordinate_screened": True,
    "global_kappa_stacked": False,
    "balance_uses_candidate_energy_only": True,
    "openai_numeric_target_inferred": False,
    "production_redistribution_gain_selected": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "held_out_pde_residual_transferred_from_parent": False,
    "material_path_integration_performed": False,
    "public_image_used": False,
    "visual_score_defined": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
}


def compact_interval_bump(radius, lo: float, hi: float) -> np.ndarray:
    """C-infinity bump on lo<r<hi, normalized to one at the midpoint."""
    radius = np.asarray(radius, dtype=float)
    if not np.all(np.isfinite(radius)) or np.any(radius < 0.0):
        raise ValueError("radius must be finite and nonnegative")
    lo = float(lo)
    hi = float(hi)
    if not (0.0 <= lo < hi <= 2.0):
        raise ValueError("radial bump window must lie in [0,2]")
    mid = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    x = (radius - mid) / half
    out = np.zeros_like(radius, dtype=float)
    mask = np.abs(x) < 1.0
    xm = x[mask]
    out[mask] = np.exp(1.0 - 1.0 / (1.0 - xm * xm))
    return out


def inner_profile(radius) -> np.ndarray:
    return compact_interval_bump(radius, *INNER_WINDOW)


def outer_profile(radius) -> np.ndarray:
    return compact_interval_bump(radius, *OUTER_WINDOW)


def _validate_gain(value: float) -> float:
    value = float(value)
    if not np.isfinite(value) or not 0.0 <= value <= 0.05:
        raise ValueError("redistribution gain must lie in the preregistered [0,.05] screen")
    return value


def _validate_gains(values: Iterable[float]) -> tuple[float, ...]:
    values = tuple(_validate_gain(value) for value in values)
    if len(values) < 2 or values[0] != 0.0:
        raise ValueError("gain grid must start at zero and contain at least two values")
    if any(a >= b for a, b in zip(values, values[1:])):
        raise ValueError("gain grid must be strictly increasing")
    return values


def _swirl_energy_moment(
    field,
    profile: Callable[[np.ndarray], np.ndarray],
    *,
    time: float,
    order: int,
) -> float:
    """Return pi*int u_theta^2 profile(r) r dr dz on r<2, |z|<2."""
    if order < 12:
        raise ValueError("moment quadrature order must be >=12")
    node, weight = leggauss(int(order))
    radius = node + 1.0
    z = 2.0 * node
    wr = weight
    wz = 2.0 * weight
    rr, zz = np.meshgrid(radius, z, indexing="ij")
    points = np.column_stack((rr.ravel(), np.zeros(rr.size), zz.ravel()))
    velocity = np.asarray(field.at_points(points, float(time)), dtype=float)
    if velocity.shape != points.shape or not np.all(np.isfinite(velocity)):
        raise RuntimeError("unexpected field velocity")
    # On the x-positive meridian y=0, u_theta is exactly Cartesian uy.
    utheta_sq = (velocity[:, 1] ** 2).reshape(rr.shape)
    radial_weight = profile(rr)
    return float(np.sum((wr[:, None] * wz[None, :]) * (np.pi * rr * utheta_sq * radial_weight)))


def _balance_coefficient(field, *, time: float = 0.25, order: int = 64) -> dict[str, float]:
    inner = _swirl_energy_moment(field, inner_profile, time=time, order=order)
    outer = _swirl_energy_moment(field, outer_profile, time=time, order=order)
    if not np.isfinite(inner) or not np.isfinite(outer) or inner <= 0.0 or outer <= 0.0:
        raise RuntimeError("redistribution energy moments must be finite and positive")
    alpha = float(inner / outer)
    defect = float(inner - alpha * outer)
    relative = float(abs(defect) / max(inner + abs(alpha * outer), 1.0e-300))
    return {
        "inner_swirl_energy_moment": inner,
        "outer_swirl_energy_moment": outer,
        "outer_balance_coefficient_alpha": alpha,
        "balanced_moment_defect": defect,
        "balanced_moment_relative_defect": relative,
        "kinetic_energy_first_derivative_at_zero": 2.0 * defect,
    }


def redistribution_profile(radius, alpha: float) -> np.ndarray:
    alpha = float(alpha)
    if not np.isfinite(alpha) or alpha <= 0.0:
        raise ValueError("alpha must be positive and finite")
    return inner_profile(radius) - alpha * outer_profile(radius)


class EnergyNeutralSwirlRedistributionField:
    """Apply one balanced radial multiplier only to the cylindrical swirl component."""

    def __init__(self, parent, gain: float, alpha: float, *, validate: bool = True):
        self.parent = parent
        self.gain = _validate_gain(gain) if validate else float(gain)
        self.alpha = float(alpha)
        if not np.isfinite(self.gain) or not np.isfinite(self.alpha) or self.alpha <= 0.0:
            raise ValueError("gain/alpha must be finite and alpha positive")
        radius = np.linspace(0.0, 2.0, 4001)
        multiplier = 1.0 + self.gain * redistribution_profile(radius, self.alpha)
        if np.min(multiplier) <= 0.0:
            raise ValueError("redistribution multiplier must stay positive on support")

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
            multiplier = 1.0 + self.gain * redistribution_profile(radius[mask], self.alpha)
            if np.any(multiplier <= 0.0):
                raise ValueError("redistribution multiplier must stay positive")
            azimuthal *= multiplier
            out[mask, 0] = rx * radial - ry * azimuthal
            out[mask, 1] = ry * radial + rx * azimuthal
        return out


def _relative(child: float, parent: float) -> float:
    return float(child / max(abs(parent), 1.0e-300) - 1.0)


def _base_temporal_field(parent, *, base_beta: float = BASE_BETA, gamma: float = BASE_GAMMA):
    return temporal.TemporalPiolaWarpField(parent, gamma, base_beta=base_beta)


def _raw_field(parent, *, gain: float, alpha: float):
    return EnergyNeutralSwirlRedistributionField(_base_temporal_field(parent), gain, alpha)


def _seed_radius_proxies(field, *, time: float):
    per_radius = {}
    all_radial = []
    all_angular = []
    for radius in SEED_RADII:
        points = []
        for z in (-0.3, 0.3):
            for theta in np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False):
                points.append([radius * np.cos(theta), radius * np.sin(theta), z])
        points = np.asarray(points, dtype=float)
        velocity = np.asarray(field.at_points(points, float(time)), dtype=float)
        x, y = points[:, 0], points[:, 1]
        rr = np.hypot(x, y)
        rx, ry = x / rr, y / rr
        radial = rx * velocity[:, 0] + ry * velocity[:, 1]
        azimuthal = -ry * velocity[:, 0] + rx * velocity[:, 1]
        angular = np.abs(azimuthal) / rr
        per_radius[f"{radius:.1f}"] = {
            "path_proxy_count": int(len(points)),
            "mean_abs_angular_rate": float(np.mean(angular)),
            "mean_inward_radial_speed": float(np.mean(-radial)),
        }
        all_radial.append(radial)
        all_angular.append(angular)
    radial = np.concatenate(all_radial)
    angular = np.concatenate(all_angular)
    return {
        "per_seed_radius": per_radius,
        "aggregate_mean_abs_angular_rate": float(np.mean(angular)),
        "aggregate_mean_inward_radial_speed": float(np.mean(-radial)),
    }


def _normalized_parameter_field(
    parent,
    *,
    base_beta: float,
    gamma: float,
    global_kappa: float,
    redistribution_gain: float,
    alpha: float,
    energy_order: int,
):
    field = temporal.TemporalPiolaWarpField(parent, gamma, base_beta=base_beta)
    if global_kappa != 0.0:
        field = global_screen.SwirlGainField(field, global_kappa, validate=False)
    if redistribution_gain != 0.0:
        field = EnergyNeutralSwirlRedistributionField(
            field,
            redistribution_gain,
            alpha,
            validate=False,
        )
    scale = fixed._normalization_scale(fixed._energy(field, time=0.25, order=energy_order))
    return fixed._ScaledField(field, scale)


def _four_column_sensitivity(parent, *, alpha: float, energy_order: int = 32):
    rng = np.random.default_rng(9174291)
    points = rng.uniform(-1.35, 1.35, size=(36, 3))
    times = np.resize(np.array([0.25, 0.375, 0.50, 0.625, 0.75], dtype=float), len(points))
    eps = 1.0e-4

    def values(base_beta: float, gamma: float, global_kappa: float, redistribution_gain: float):
        field = _normalized_parameter_field(
            parent,
            base_beta=base_beta,
            gamma=gamma,
            global_kappa=global_kappa,
            redistribution_gain=redistribution_gain,
            alpha=alpha,
            energy_order=energy_order,
        )
        rows = [field.at_points(point[None, :], float(time))[0] for point, time in zip(points, times)]
        return np.asarray(rows, dtype=float).ravel()

    beta_column = (
        values(BASE_BETA + eps, BASE_GAMMA, 0.0, 0.0)
        - values(BASE_BETA - eps, BASE_GAMMA, 0.0, 0.0)
    ) / (2.0 * eps)
    gamma_column = (
        values(BASE_BETA, BASE_GAMMA + eps, 0.0, 0.0)
        - values(BASE_BETA, BASE_GAMMA - eps, 0.0, 0.0)
    ) / (2.0 * eps)
    zero = values(BASE_BETA, BASE_GAMMA, 0.0, 0.0)
    global_column = (values(BASE_BETA, BASE_GAMMA, eps, 0.0) - zero) / eps
    redistribution_column = (values(BASE_BETA, BASE_GAMMA, 0.0, eps) - zero) / eps

    columns = [beta_column, gamma_column, global_column, redistribution_column]
    norms = np.asarray([np.linalg.norm(column) for column in columns], dtype=float)
    if np.any(norms <= 1.0e-14):
        raise RuntimeError("degenerate sensitivity column")
    matrix = np.column_stack([column / norm for column, norm in zip(columns, norms)])
    singular = np.linalg.svd(matrix, compute_uv=False)
    existing = np.column_stack((beta_column, gamma_column, global_column))
    coeff, *_ = np.linalg.lstsq(existing, redistribution_column, rcond=None)
    residual = redistribution_column - existing @ coeff
    return {
        "parameter_count_in_diagnostic": 4,
        "production_recommendation_parameter_count": 3,
        "production_interpretation": "redistribution replaces global kappa; it is not an added fourth production degree",
        "parameters": ["constant_beta", "endpoint_even_time_gamma", "global_swirl_kappa_reference", "energy_neutral_swirl_redistribution"],
        "sample_count": int(len(points)),
        "finite_difference_step": eps,
        "normalized_rank": int(np.linalg.matrix_rank(matrix)),
        "normalized_singular_values": singular.tolist(),
        "normalized_condition_number": float(singular[0] / singular[-1]),
        "normalized_column_cosines": (matrix.T @ matrix).tolist(),
        "raw_column_norms": norms.tolist(),
        "redistribution_column_novelty_fraction_outside_beta_gamma_global_span": float(
            np.linalg.norm(residual) / max(np.linalg.norm(redistribution_column), 1.0e-300)
        ),
    }


def audit_energy_neutral_swirl_redistribution(
    *,
    gains: Iterable[float] = DEFAULT_GAINS,
    times: Iterable[float] = DEFAULT_TIMES,
    energy_times: Iterable[float] = DEFAULT_ENERGY_TIMES,
    energy_orders: tuple[int, int] = (48, 64),
    grid_size: int = DEFAULT_GRID_SIZE,
):
    gains = _validate_gains(gains)
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
    base_raw = _base_temporal_field(parent)
    balance_coarse = _balance_coefficient(base_raw, order=energy_orders[0])
    balance_fine = _balance_coefficient(base_raw, order=energy_orders[1])
    alpha = balance_fine["outer_balance_coefficient_alpha"]
    balance_refinement = abs(
        balance_coarse["outer_balance_coefficient_alpha"] - alpha
    ) / max(abs(alpha), 1.0e-300)

    profile_seed_values = {
        f"{radius:.1f}": float(redistribution_profile(np.array([radius]), alpha)[0])
        for radius in SEED_RADII
    }
    support_grid = np.linspace(0.0, 2.0, 4001)
    profile_grid = redistribution_profile(support_grid, alpha)

    raw_energy = {}
    energy_refinement_max = 0.0
    for gain in gains:
        raw = _raw_field(parent, gain=gain, alpha=alpha)
        raw_energy[gain] = {}
        for time in energy_times:
            coarse = fixed._energy(raw, time=time, order=energy_orders[0])
            fine = fixed._energy(raw, time=time, order=energy_orders[1])
            raw_energy[gain][time] = fine
            energy_refinement_max = max(
                energy_refinement_max,
                abs(coarse - fine) / max(abs(fine), 1.0e-12),
            )

    base_scale = fixed._normalization_scale(raw_energy[0.0][0.25])
    baseline = fixed._ScaledField(_raw_field(parent, gain=0.0, alpha=alpha), base_scale)
    baseline_morphology = {
        time: fixed.transfer._morphology(baseline, time=time, grid_size=grid_size)
        for time in times
    }
    baseline_proxy = {time: _seed_radius_proxies(baseline, time=time) for time in times}

    rows = []
    crossings = []
    for gain in gains:
        raw = _raw_field(parent, gain=gain, alpha=alpha)
        scale = fixed._normalization_scale(raw_energy[gain][0.25])
        field = fixed._ScaledField(raw, scale)
        normalized_energy = {time: scale * scale * raw_energy[gain][time] for time in energy_times}
        morphology = {
            time: fixed.transfer._morphology(field, time=time, grid_size=grid_size)
            for time in times
        }
        proxy = {time: _seed_radius_proxies(field, time=time) for time in times}

        radius_angular_changes = {}
        radius_inward_changes = {}
        for radius in SEED_RADII:
            key = f"{radius:.1f}"
            radius_angular_changes[key] = [
                _relative(
                    proxy[time]["per_seed_radius"][key]["mean_abs_angular_rate"],
                    baseline_proxy[time]["per_seed_radius"][key]["mean_abs_angular_rate"],
                )
                for time in times
            ]
            radius_inward_changes[key] = [
                _relative(
                    proxy[time]["per_seed_radius"][key]["mean_inward_radial_speed"],
                    baseline_proxy[time]["per_seed_radius"][key]["mean_inward_radial_speed"],
                )
                for time in times
            ]

        aggregate_angular = [
            _relative(proxy[time]["aggregate_mean_abs_angular_rate"], baseline_proxy[time]["aggregate_mean_abs_angular_rate"])
            for time in times
        ]
        aggregate_inward = [
            _relative(proxy[time]["aggregate_mean_inward_radial_speed"], baseline_proxy[time]["aggregate_mean_inward_radial_speed"])
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
        multiplier_min = float(np.min(1.0 + gain * profile_grid))
        multiplier_max = float(np.max(1.0 + gain * profile_grid))
        reference_energy_pass = abs(normalized_energy[0.25] - 1.0) <= fixed.REFERENCE_ENERGY_TOLERANCE
        energy_range_pass = all(
            fixed.VALIDATION_ENERGY_RANGE[0] <= value <= fixed.VALIDATION_ENERGY_RANGE[1]
            for value in normalized_energy.values()
        )
        preflight = bool(
            reference_energy_pass
            and energy_range_pass
            and multiplier_min > 0.0
            and structure["outside_support_max_abs_velocity"] <= 1.0e-12
            and structure["representative_inward_swirl_bipolar_signs_all_pass"]
            and structure["cartesian_fd_divergence_max_abs"] <= temporal.DIVERGENCE_PREFLIGHT_MAX
            and structure["all_coordinate_maps_monotone"]
        )
        clean = bool(
            gain > 0.0
            and preflight
            and min(radius_angular_changes["0.6"]) >= INNER_GAIN_FLOOR - 1.0e-12
            and min(radius_angular_changes["0.9"]) >= MID_GAIN_FLOOR - 1.0e-12
            and max(radius_angular_changes["1.2"]) <= OUTER_GAIN_CEILING + 1.0e-12
            and abs(scale - 1.0) <= NORMALIZATION_DEVIATION_MAX + 1.0e-12
            and min(aggregate_inward) >= -INWARD_SPEED_LOSS_MAX - 1.0e-12
            and max(abs(value) for value in axial) <= AXIAL_RMS_ABS_CHANGE_MAX + 1.0e-12
            and max(radial) <= RADIAL_RMS_GROWTH_MAX + 1.0e-12
            and max(collar) <= COLLAR_RATIO_MAX + 1.0e-12
        )
        if clean:
            crossings.append(gain)
        rows.append({
            "redistribution_gain": gain,
            "common_energy_normalization_scale": scale,
            "scale_relative_to_gain0": float(scale / base_scale - 1.0),
            "normalized_energy": {str(time): value for time, value in normalized_energy.items()},
            "minimum_swirl_multiplier_on_support": multiplier_min,
            "maximum_swirl_multiplier_on_support": multiplier_max,
            "seed_radius_angular_rate_relative_changes": radius_angular_changes,
            "seed_radius_inward_speed_relative_changes": radius_inward_changes,
            "aggregate_angular_rate_relative_changes": aggregate_angular,
            "aggregate_inward_speed_relative_changes": aggregate_inward,
            "axial_rms_relative_changes_vs_gain0": axial,
            "radial_rms_relative_changes_vs_gain0": radial,
            "q90_deltas_over_support_vs_gain0": q90,
            "q99_deltas_over_support_vs_gain0": q99,
            "collar_enstrophy_ratios_vs_gain0": collar,
            "preflight_pass": preflight,
            "clean_energy_neutral_redistribution_rule_pass": clean,
            "structure": structure,
        })

    crossing = min(crossings) if crossings else None
    routing = (
        "A clean first-order energy-neutral inner/mid redistribution crossing exists. Keep this as a replacement candidate for global kappa, not an added fourth production degree. Before promotion, materialize one governed child, rerun actual material paths/3D streamlines, and rebuild pressure/restricted forcing plus fresh held-out full momentum; do not transfer parent PDE evidence."
        if crossing is not None
        else
        "No preregistered clean crossing. Stop this redistribution profile and do not widen the gain grid post hoc. Use the reported radius/energy/morphology diagnostics to redesign at most one next bounded redistribution shape rather than adding arbitrary swirl basis dimensions."
    )
    return {
        "task_id": TASK_ID,
        "issue": 446,
        "parent_id": PARENT_ID,
        "frozen_temporal_piola": {"base_beta": BASE_BETA, "gamma": BASE_GAMMA},
        "agent9_routing_evidence": {
            "pr": 445,
            "meaning": "candidate-side frozen path attribution: remaining winding deficit is mainly inner r=.6, then r=.9, while r=1.2 is surplus versus ST006 benchmark",
            "st006_is_openai_truth": False,
        },
        "profile_design": {
            "inner_window": list(INNER_WINDOW),
            "outer_window": list(OUTER_WINDOW),
            "formula": "h(r)=g_inner(r)-alpha*g_outer(r); each g is exp(1-1/(1-x^2)) on its open radial window and zero outside",
            "seed_radius_profile_values": profile_seed_values,
            "profile_min": float(np.min(profile_grid)),
            "profile_max": float(np.max(profile_grid)),
        },
        "energy_balance_coarse": balance_coarse,
        "energy_balance_fine": balance_fine,
        "balance_coefficient_relative_refinement_change": float(balance_refinement),
        "maximum_energy_quadrature_relative_refinement_change": float(energy_refinement_max),
        "clean_rule": {
            "inner_r06_angular_gain_floor": INNER_GAIN_FLOOR,
            "mid_r09_angular_gain_floor": MID_GAIN_FLOOR,
            "outer_r12_angular_gain_ceiling": OUTER_GAIN_CEILING,
            "absolute_common_scale_deviation_max": NORMALIZATION_DEVIATION_MAX,
            "aggregate_inward_speed_loss_max": INWARD_SPEED_LOSS_MAX,
            "axial_rms_absolute_change_max": AXIAL_RMS_ABS_CHANGE_MAX,
            "radial_rms_growth_max": RADIAL_RMS_GROWTH_MAX,
            "collar_ratio_max": COLLAR_RATIO_MAX,
        },
        "rows": rows,
        "smallest_preregistered_clean_energy_neutral_redistribution_crossing": crossing,
        "four_column_velocity_sensitivity": _four_column_sensitivity(parent, alpha=alpha),
        "truth_boundary": dict(TRUTH_BOUNDARY),
        "routing": routing,
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("artifacts/st048_agent7_energy_neutral_swirl_redistribution/report.json"))
    args = parser.parse_args(argv)
    report = audit_energy_neutral_swirl_redistribution()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "task_id": report["task_id"],
        "crossing": report["smallest_preregistered_clean_energy_neutral_redistribution_crossing"],
        "alpha": report["energy_balance_fine"]["outer_balance_coefficient_alpha"],
        "balance_relative_defect": report["energy_balance_fine"]["balanced_moment_relative_defect"],
        "sensitivity_rank": report["four_column_velocity_sensitivity"]["normalized_rank"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
