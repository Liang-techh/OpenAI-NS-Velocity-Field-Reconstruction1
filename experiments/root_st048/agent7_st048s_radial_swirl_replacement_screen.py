"""Screen one radial-localized swirl coordinate as a replacement for global kappa.

Constrained Agent 7 follows Agent-9 PR #435.  The global post-Piola
``kappa=.05`` crossing raises cumulative material-path winding only modestly while
common energy renormalization weakens inward contraction and axial-pair
separation.  Do not add another parameter on top of beta/gamma/kappa here.
Instead replace the global azimuthal multiplier by exactly one compact radial
shape coordinate

    u_theta -> (1 + kappa_ring * g(r^2/4)) u_theta,

where ``g`` is a C-infinity bump on 0<r<2, zero on the axis and at the radial
support edge, with unit peak at r=sqrt(2).  The temporal Piola schedule remains
frozen at beta(t)=.075+.025*(4t-2)^2 and every row gets one positive common
normalization restoring E(.25)=1.

Any axisymmetric azimuthal-only multiplier contributes no cylindrical divergence
term.  It still changes vorticity, nonlinear momentum and energy, so pressure,
restricted forcing and held-out PDE evidence are not transferred.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import numpy as np

import agent7_st048s_piola_swirl_gain_screen as global_screen
import agent7_st048s_piola_temporal_screen as temporal
import agent7_st048s_piola_warp_screen as fixed

TASK_ID = "CR003-ST048S-RADIAL-SWIRL-REPLACEMENT-064"
PARENT_ID = global_screen.PARENT_ID
BASE_BETA = global_screen.BASE_BETA
BASE_GAMMA = global_screen.BASE_GAMMA
GLOBAL_REFERENCE_KAPPA = 0.05
DEFAULT_RING_GAINS = (0.0, 0.05, 0.075, 0.10, 0.15, 0.20)
DEFAULT_TIMES = global_screen.DEFAULT_TIMES
DEFAULT_ENERGY_TIMES = global_screen.DEFAULT_ENERGY_TIMES
DEFAULT_GRID_SIZE = global_screen.DEFAULT_GRID_SIZE
ANGULAR_GAIN_FLOOR = 0.02
INWARD_SPEED_LOSS_MAX = 0.015
AXIAL_RMS_ABS_CHANGE_MAX = 0.02
RADIAL_RMS_GROWTH_MAX = 0.03
COLLAR_RATIO_MAX = 1.25

TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "candidate_artifact_changed": False,
    "global_kappa_promoted": False,
    "one_radial_swirl_shape_coordinate_screened": True,
    "localized_coordinate_is_replacement_not_additional_production_parameter": True,
    "production_ring_gain_selected": False,
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


def radial_ring_profile(radius) -> np.ndarray:
    """C-infinity ring profile in q=r^2/4, normalized to g(1/2)=1."""
    radius = np.asarray(radius, dtype=float)
    if not np.all(np.isfinite(radius)) or np.any(radius < 0.0):
        raise ValueError("radius must be finite and nonnegative")
    q = 0.25 * radius * radius
    out = np.zeros_like(q, dtype=float)
    mask = (q > 0.0) & (q < 1.0)
    qm = q[mask]
    # exp(1 - 1/(4q) - 1/(4(1-q))) is C-infinity after zero extension,
    # is positive on (0,1), and equals one at q=1/2.
    out[mask] = np.exp(1.0 - 0.25 / qm - 0.25 / (1.0 - qm))
    return out


def _validate_ring_gain(value: float) -> float:
    value = float(value)
    if not np.isfinite(value) or not 0.0 <= value <= 0.20:
        raise ValueError("ring gain must lie in the preregistered [0,.20] screen")
    return value


def _validate_ring_gains(values: Iterable[float]) -> tuple[float, ...]:
    values = tuple(_validate_ring_gain(value) for value in values)
    if len(values) < 2 or values[0] != 0.0:
        raise ValueError("ring-gain grid must start at zero and contain at least two values")
    if any(a >= b for a, b in zip(values, values[1:])):
        raise ValueError("ring-gain grid must be strictly increasing")
    return values


class RadialSwirlGainField:
    """Multiply only u_theta by 1+kappa_ring*g(r), leaving u_r/u_z fixed."""

    def __init__(self, parent, ring_gain: float, *, validate: bool = True):
        self.parent = parent
        self.ring_gain = _validate_ring_gain(ring_gain) if validate else float(ring_gain)
        if not np.isfinite(self.ring_gain):
            raise ValueError("ring gain must be finite")

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
            multiplier = 1.0 + self.ring_gain * radial_ring_profile(radius[mask])
            if np.any(multiplier <= 0.0):
                raise ValueError("radial swirl multiplier must stay positive")
            azimuthal *= multiplier
            out[mask, 0] = rx * radial - ry * azimuthal
            out[mask, 1] = ry * radial + rx * azimuthal
        return out


def _raw_ring_field(parent, *, ring_gain: float):
    temporal_field = temporal.TemporalPiolaWarpField(parent, BASE_GAMMA, base_beta=BASE_BETA)
    return RadialSwirlGainField(temporal_field, ring_gain)


def _scaled(raw, *, energy_order: int):
    energy = fixed._energy(raw, time=0.25, order=energy_order)
    scale = fixed._normalization_scale(energy)
    return fixed._ScaledField(raw, scale), scale, energy


def _relative(child: float, parent: float) -> float:
    return float(child / max(abs(parent), 1.0e-300) - 1.0)


def _normalized_combined_field(
    parent,
    *,
    base_beta: float,
    gamma: float,
    global_kappa: float,
    ring_gain: float,
    energy_order: int,
):
    field = temporal.TemporalPiolaWarpField(parent, gamma, base_beta=base_beta)
    if global_kappa != 0.0:
        field = global_screen.SwirlGainField(field, global_kappa, validate=False)
    if ring_gain != 0.0:
        field = RadialSwirlGainField(field, ring_gain, validate=False)
    scale = fixed._normalization_scale(fixed._energy(field, time=0.25, order=energy_order))
    return fixed._ScaledField(field, scale)


def _four_column_sensitivity(parent, *, energy_order: int = 32):
    """Compare ring shape with beta/gamma/global-swirl directions on fixed probes."""
    rng = np.random.default_rng(9174281)
    points = rng.uniform(-1.35, 1.35, size=(36, 3))
    times = np.resize(np.array([0.25, 0.375, 0.50, 0.625, 0.75], dtype=float), len(points))
    eps = 1.0e-4

    def values(base_beta: float, gamma: float, global_kappa: float, ring_gain: float):
        field = _normalized_combined_field(
            parent,
            base_beta=base_beta,
            gamma=gamma,
            global_kappa=global_kappa,
            ring_gain=ring_gain,
            energy_order=energy_order,
        )
        rows = [
            field.at_points(point[None, :], float(time))[0]
            for point, time in zip(points, times)
        ]
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
    ring_column = (values(BASE_BETA, BASE_GAMMA, 0.0, eps) - zero) / eps

    columns = [beta_column, gamma_column, global_column, ring_column]
    norms = np.asarray([np.linalg.norm(column) for column in columns], dtype=float)
    if np.any(norms <= 1.0e-14):
        raise RuntimeError("degenerate sensitivity column")
    matrix = np.column_stack([column / norm for column, norm in zip(columns, norms)])
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix))
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0.0 else float("inf")
    cosine = matrix.T @ matrix

    existing = np.column_stack((beta_column, gamma_column, global_column))
    coeff, *_ = np.linalg.lstsq(existing, ring_column, rcond=None)
    residual = ring_column - existing @ coeff
    novelty = float(np.linalg.norm(residual) / max(np.linalg.norm(ring_column), 1.0e-300))
    return {
        "parameter_count_in_diagnostic": 4,
        "production_recommendation_parameter_count": 3,
        "production_interpretation": "ring gain replaces global kappa; it is not an added fourth production degree",
        "parameters": ["constant_beta", "endpoint_even_time_gamma", "global_swirl_kappa_reference", "radial_ring_swirl_replacement"],
        "sample_count": int(len(points)),
        "finite_difference_step": eps,
        "normalized_rank": rank,
        "normalized_singular_values": singular.tolist(),
        "normalized_condition_number": condition,
        "normalized_column_cosines": cosine.tolist(),
        "raw_column_norms": norms.tolist(),
        "ring_column_novelty_fraction_outside_beta_gamma_global_span": novelty,
    }


def audit_radial_swirl_replacement_screen(
    *,
    ring_gains: Iterable[float] = DEFAULT_RING_GAINS,
    times: Iterable[float] = DEFAULT_TIMES,
    energy_times: Iterable[float] = DEFAULT_ENERGY_TIMES,
    energy_orders: tuple[int, int] = (48, 64),
    grid_size: int = DEFAULT_GRID_SIZE,
):
    ring_gains = _validate_ring_gains(ring_gains)
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
    for gain in ring_gains:
        raw = _raw_ring_field(parent, ring_gain=gain)
        raw_energy[gain] = {}
        for time in energy_times:
            coarse = fixed._energy(raw, time=time, order=energy_orders[0])
            fine = fixed._energy(raw, time=time, order=energy_orders[1])
            raw_energy[gain][time] = fine
            refinement_max = max(
                refinement_max,
                abs(coarse - fine) / max(abs(fine), 1.0e-12),
            )

    base_raw = _raw_ring_field(parent, ring_gain=0.0)
    base_scale = fixed._normalization_scale(raw_energy[0.0][0.25])
    baseline = fixed._ScaledField(base_raw, base_scale)
    baseline_morphology = {
        time: fixed.transfer._morphology(baseline, time=time, grid_size=grid_size)
        for time in times
    }
    baseline_proxy = {time: global_screen._winding_proxies(baseline, time=time) for time in times}

    global_raw = global_screen._raw_field(parent, kappa=GLOBAL_REFERENCE_KAPPA)
    global_scale = fixed._normalization_scale(fixed._energy(global_raw, time=0.25, order=energy_orders[1]))
    global_reference = fixed._ScaledField(global_raw, global_scale)
    global_proxy = {time: global_screen._winding_proxies(global_reference, time=time) for time in times}
    global_proxy_changes = {
        "mean_angular_rate_relative_changes_vs_kappa0": [
            _relative(global_proxy[time]["mean_abs_angular_rate"], baseline_proxy[time]["mean_abs_angular_rate"])
            for time in times
        ],
        "mean_inward_speed_relative_changes_vs_kappa0": [
            _relative(global_proxy[time]["mean_inward_radial_speed"], baseline_proxy[time]["mean_inward_radial_speed"])
            for time in times
        ],
    }

    rows = []
    crossings = []
    for gain in ring_gains:
        raw = _raw_ring_field(parent, ring_gain=gain)
        scale = fixed._normalization_scale(raw_energy[gain][0.25])
        field = fixed._ScaledField(raw, scale)
        normalized_energy = {
            time: scale * scale * raw_energy[gain][time]
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
        proxy = {time: global_screen._winding_proxies(field, time=time) for time in times}
        angular_changes = [
            _relative(proxy[time]["mean_abs_angular_rate"], baseline_proxy[time]["mean_abs_angular_rate"])
            for time in times
        ]
        inward_changes = [
            _relative(proxy[time]["mean_inward_radial_speed"], baseline_proxy[time]["mean_inward_radial_speed"])
            for time in times
        ]
        swirl_ratio_changes = [
            _relative(proxy[time]["swirl_to_poloidal_rms_ratio"], baseline_proxy[time]["swirl_to_poloidal_rms_ratio"])
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
            gain > 0.0
            and preflight
            and min(angular_changes) >= ANGULAR_GAIN_FLOOR - 1.0e-12
            and min(inward_changes) >= -INWARD_SPEED_LOSS_MAX - 1.0e-12
            and max(abs(value) for value in axial) <= AXIAL_RMS_ABS_CHANGE_MAX + 1.0e-12
            and max(radial) <= RADIAL_RMS_GROWTH_MAX + 1.0e-12
            and max(collar) <= COLLAR_RATIO_MAX + 1.0e-12
        )
        if clean:
            crossings.append(gain)
        rows.append({
            "ring_gain": gain,
            "ring_profile_values_at_seed_radii": {
                str(radius): float(radial_ring_profile(np.asarray([radius]))[0])
                for radius in (0.6, 0.9, 1.2)
            },
            "common_energy_normalization_scale": scale,
            "normalized_energy": {str(key): value for key, value in normalized_energy.items()},
            "reference_energy_pass": reference_energy_pass,
            "energy_range_pass": energy_range_pass,
            "winding_proxy_by_time": {str(key): value for key, value in proxy.items()},
            "mean_angular_rate_relative_changes_vs_kappa0": angular_changes,
            "mean_inward_speed_relative_changes_vs_kappa0": inward_changes,
            "swirl_to_poloidal_relative_changes_vs_kappa0": swirl_ratio_changes,
            "axial_rms_relative_changes_vs_kappa0": axial,
            "radial_rms_relative_changes_vs_kappa0": radial,
            "q90_deltas_over_support_vs_kappa0": q90,
            "q99_deltas_over_support_vs_kappa0": q99,
            "collar_enstrophy_ratios_vs_kappa0": collar,
            "structure": structure,
            "preflight_pass": preflight,
            "clean_radial_swirl_replacement_rule_pass": clean,
        })

    sensitivity = _four_column_sensitivity(parent)
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
        "replacement_semantics": "ring gain replaces the global kappa production coordinate; diagnostic four-column rank includes both only to measure novelty",
        "ring_profile": "g(q)=exp(1-1/(4q)-1/(4(1-q))) on 0<q<1, q=r^2/4; zero otherwise",
        "ring_gain_grid": list(ring_gains),
        "times": list(times),
        "energy_times": list(energy_times),
        "energy_orders": list(energy_orders),
        "grid_size": grid_size,
        "preregistered_clean_rule": {
            "minimum_mean_angular_rate_relative_gain_each_time": ANGULAR_GAIN_FLOOR,
            "maximum_fixed_probe_inward_speed_relative_loss_each_time": INWARD_SPEED_LOSS_MAX,
            "maximum_absolute_axial_vorticity_rms_relative_change": AXIAL_RMS_ABS_CHANGE_MAX,
            "maximum_radial_vorticity_rms_growth": RADIAL_RMS_GROWTH_MAX,
            "maximum_support_collar_enstrophy_ratio": COLLAR_RATIO_MAX,
            "requires_energy_support_core_sign_divergence_preflight": True,
        },
        "energy_quadrature_max_relative_refinement_change": refinement_max,
        "baseline_winding_proxy_by_time": {str(key): value for key, value in baseline_proxy.items()},
        "global_kappa_reference": {
            "kappa": GLOBAL_REFERENCE_KAPPA,
            "common_energy_normalization_scale": global_scale,
            "proxy_by_time": {str(key): value for key, value in global_proxy.items()},
            **global_proxy_changes,
            "agent9_material_path_receipt_is_sibling_evidence_not_recomputed_here": True,
        },
        "rows": rows,
        "smallest_preregistered_clean_radial_swirl_replacement_crossing": min(crossings) if crossings else None,
        "four_column_velocity_sensitivity": sensitivity,
        "routing": (
            "If a nonzero clean crossing exists, prefer the radial ring coordinate as a replacement "
            "candidate for the global kappa direction, not as an added fourth production parameter. "
            "Route the first crossing to an independent Agent-9 material-path replay before any "
            "materialization. Only if trajectory winding/contraction improve together should Agent 1/2 "
            "consider a governed child with a rebuilt pressure/restricted force and fresh held-out full "
            "momentum. If no crossing exists, retain the global kappa result as capacity evidence and do "
            "not widen the radial-gain grid post hoc."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--grid-size", type=int, default=DEFAULT_GRID_SIZE)
    args = parser.parse_args()
    report = audit_radial_swirl_replacement_screen(grid_size=args.grid_size)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)


if __name__ == "__main__":
    main()
