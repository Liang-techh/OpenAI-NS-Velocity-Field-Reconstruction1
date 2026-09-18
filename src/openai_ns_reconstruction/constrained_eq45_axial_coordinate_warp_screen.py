"""Screen one support-preserving axial coordinate-geometry degree of freedom.

The preceding axial-envelope line can move the q90/high-enstrophy core only at
large alpha, where the support collar is strongly loaded.  This follow-up does
not add another basis.  It instead tests one smooth axial coordinate warp that
is identity at the registered physical support boundary and is constructed with
the Piola transform so an exactly divergence-free parent remains divergence
free under the coordinate change.

For |z| < Z, Z=2, write s=z/Z and

    h_beta(z) = z * (1 - beta * (1-s^2)^4),

with h_beta(z)=z outside.  The diagnostic velocity is

    u_beta(x,y,z,t) = [h'(z) u_x(x,y,h(z),t),
                       h'(z) u_y(x,y,h(z),t),
                       u_z(x,y,h(z),t)].

This is the contravariant Piola pullback for y=(x,y,h(z)); consequently
``div u_beta = h'(z) (div u)(x,y,h(z),t)`` wherever the parent is smooth.  The
fourth power makes the warp and its low derivatives approach the identity in
the support collar.  beta is an autonomous capacity parameter, not recovered
from an OpenAI image or source coefficient.

As in the preceding envelope screen, each nonzero beta receives one positive
common velocity scale to restore the registered reference energy.  That common
scale preserves instantaneous streamline directions and normalized
vorticity-location fingerprints but is not a Navier--Stokes invariance.  This
module therefore does not create a production candidate, fit pressure/forcing,
or claim held-out PDE evidence or OpenAI correspondence.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_axial_envelope_small_alpha_screen import (
    _axisymmetric_grid,
    _constraints,
)
from .constrained_eq45_bipolar_axial_shoulder_poloidal_capacity import (
    _vorticity_morphology,
)
from .constrained_eq45_bipolar_compact_poloidal_capacity import (
    _cylindrical_components,
    _probe_groups,
    _source_field,
)

TASK_ID = "CR003-AXIAL-COORDINATE-WARP-SCREEN-056"
PHYSICAL_Z_CUT = 2.0
DEFAULT_BETAS = (0.0, 0.10, 0.20, 0.30, 0.40, 0.50)
DEFAULT_MORPHOLOGY_TIMES = (0.25, 0.50, 0.75)
MORPHOLOGY_ZERO_TOLERANCE = 1.0e-12
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "canonical_velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_spatial_basis_added": False,
    "diagnostic_coordinate_family_added": True,
    "joint_renormalized_child_materialized": False,
    "production_geometry_parameter_selected": False,
    "production_bound_selected": False,
    "candidate_sha_created": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "public_image_used": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
}


def _validate_beta(beta: float) -> float:
    beta = float(beta)
    if not np.isfinite(beta) or not 0.0 <= beta <= 0.50:
        raise ValueError("beta must lie in [0,0.50] for this preregistered screen")
    return beta


def _validate_beta_grid(betas: Iterable[float]) -> tuple[float, ...]:
    values = tuple(_validate_beta(value) for value in betas)
    if len(values) < 2 or values[0] != 0.0:
        raise ValueError("beta grid must start at 0 and contain at least two values")
    if any(a >= b for a, b in zip(values, values[1:])):
        raise ValueError("beta grid must be strictly increasing")
    return values


def _warp_z_and_jacobian(z, beta: float) -> tuple[np.ndarray, np.ndarray]:
    """Return h_beta(z) and h_beta'(z) for the fixed-support warp."""
    beta = _validate_beta(beta)
    z_arr = np.asarray(z, dtype=float)
    if not np.all(np.isfinite(z_arr)):
        raise ValueError("z must be finite")
    s = z_arr / PHYSICAL_Z_CUT
    inside = np.abs(s) < 1.0
    one_minus = np.maximum(1.0 - s * s, 0.0)
    bump = one_minus**4
    mapped = np.where(inside, z_arr * (1.0 - beta * bump), z_arr)
    jacobian_inside = 1.0 - beta * bump + 8.0 * beta * s * s * one_minus**3
    jacobian = np.where(inside, jacobian_inside, 1.0)
    if not np.all(np.isfinite(mapped)) or not np.all(np.isfinite(jacobian)):
        raise RuntimeError("coordinate warp became nonfinite")
    if np.any(jacobian <= 0.0):
        raise RuntimeError("coordinate warp lost orientation")
    return mapped, jacobian


class _AxialCoordinateWarpProxy:
    """Evaluate the Piola-transformed diagnostic field without creating a candidate."""

    def __init__(self, field, beta: float):
        self._field = field
        self.beta = _validate_beta(beta)

    def __getattr__(self, name):
        return getattr(self._field, name)

    def at_points(self, points, time: float):
        points_arr = np.asarray(points, dtype=float)
        if points_arr.ndim < 1 or points_arr.shape[-1] != 3:
            raise ValueError("points must end in Cartesian component dimension 3")
        if not np.all(np.isfinite(points_arr)):
            raise ValueError("points must be finite")
        mapped_z, jacobian = _warp_z_and_jacobian(points_arr[..., 2], self.beta)
        mapped_points = np.array(points_arr, copy=True)
        mapped_points[..., 2] = mapped_z
        parent_velocity = np.asarray(self._field.at_points(mapped_points, float(time)), dtype=float)
        if parent_velocity.shape != points_arr.shape or not np.all(np.isfinite(parent_velocity)):
            raise RuntimeError("unexpected/nonfinite parent velocity")
        out = np.array(parent_velocity, copy=True)
        out[..., 0] *= jacobian
        out[..., 1] *= jacobian
        return out


def _energy(proxy, *, time: float, order: int) -> float:
    points, volume = _axisymmetric_grid(order)
    velocity = np.asarray(proxy.at_points(points, float(time)), dtype=float)
    if velocity.shape != (len(points), 3) or not np.all(np.isfinite(velocity)):
        raise RuntimeError("unexpected/nonfinite velocity in energy quadrature")
    return 0.5 * float(np.sum(volume * np.sum(velocity * velocity, axis=1)))


def _normalization_scale(raw_reference_energy: float, target: float) -> float:
    raw_reference_energy = float(raw_reference_energy)
    target = float(target)
    if not np.isfinite(raw_reference_energy) or raw_reference_energy <= 0.0:
        raise ValueError("raw reference energy must be positive and finite")
    if not np.isfinite(target) or target <= 0.0:
        raise ValueError("target energy must be positive and finite")
    return float(np.sqrt(target / raw_reference_energy))


def _morphology(proxy, *, time: float, grid_size: int) -> dict[str, float]:
    return _vorticity_morphology(proxy, time=float(time), grid_size=int(grid_size))


def _structure_checks(field, *, beta: float) -> dict[str, Any]:
    proxy = _AxialCoordinateWarpProxy(field, beta)
    replay_points = np.concatenate(tuple(_probe_groups().values()))
    replay_error = max(
        float(
            np.max(
                np.abs(
                    _AxialCoordinateWarpProxy(field, 0.0).at_points(replay_points, time)
                    - field.at_points(replay_points, time)
                )
            )
        )
        for time in (0.25, 0.50, 0.75)
    )

    support_points = np.array(
        [
            [2.01, 0.0, 0.0],
            [-2.01, 0.0, 0.0],
            [0.4, 0.2, 2.01],
            [0.4, -0.2, -2.01],
        ],
        dtype=float,
    )
    support_max = float(np.max(np.abs(proxy.at_points(support_points, 0.5))))

    check = np.array(
        [
            [0.31, 0.22, 0.17],
            [0.62, -0.27, -0.41],
            [0.88, 0.31, 0.67],
            [1.22, -0.36, -1.22],
        ],
        dtype=float,
    )
    h = 1.0e-5
    divergence = np.zeros(len(check), dtype=float)
    for axis in range(3):
        offset = np.zeros_like(check)
        offset[:, axis] = h
        plus = proxy.at_points(check + offset, 0.5)
        minus = proxy.at_points(check - offset, 0.5)
        divergence += (plus[:, axis] - minus[:, axis]) / (2.0 * h)

    core_rows = []
    for time in np.linspace(field.time_start, field.time_end, 6):
        tau = 1.0 - float(time)
        point = np.array([[0.1 * np.sqrt(tau), 0.0, 0.1 * tau**0.495]])
        velocity = proxy.at_points(point, float(time))[0]
        ur, ut, uz = _cylindrical_components(point, velocity.reshape(1, 3))
        core_rows.append(bool(ur[0] < 0.0 and ut[0] > 0.0 and uz[0] > 0.0))

    z_grid = np.linspace(-PHYSICAL_Z_CUT, PHYSICAL_Z_CUT, 4001)
    mapped, jacobian = _warp_z_and_jacobian(z_grid, beta)
    return {
        "beta_zero_parent_replay_max_abs_velocity_error": replay_error,
        "outside_physical_support_max_abs_velocity": support_max,
        "cartesian_fd_divergence_max_abs": float(np.max(np.abs(divergence))),
        "representative_core_signs_all_pass": bool(all(core_rows)),
        "minimum_coordinate_jacobian": float(np.min(jacobian)),
        "maximum_coordinate_jacobian": float(np.max(jacobian)),
        "map_monotone_on_support": bool(np.all(np.diff(mapped) > 0.0)),
        "support_endpoint_map_error": float(
            max(abs(mapped[0] + PHYSICAL_Z_CUT), abs(mapped[-1] - PHYSICAL_Z_CUT))
        ),
    }


def audit_axial_coordinate_warp_screen(
    *,
    betas: Iterable[float] = DEFAULT_BETAS,
    quadrature_orders: Iterable[int] = (64, 96),
    morphology_times: Iterable[float] = DEFAULT_MORPHOLOGY_TIMES,
    morphology_grid_size: int = 33,
) -> dict[str, Any]:
    betas = _validate_beta_grid(betas)
    orders = tuple(int(value) for value in quadrature_orders)
    if len(orders) != 2 or not 16 <= orders[0] < orders[1]:
        raise ValueError("quadrature_orders must contain two increasing orders >=16")
    morph_times = tuple(float(value) for value in morphology_times)
    if not morph_times or any(not np.isfinite(value) for value in morph_times):
        raise ValueError("morphology_times must be finite and nonempty")
    if morphology_grid_size < 17 or morphology_grid_size % 2 == 0:
        raise ValueError("morphology_grid_size must be odd and >=17")

    cfg = _constraints()
    nontriviality = cfg["nontriviality"]
    validation_times = tuple(float(value) for value in cfg["validation"]["times"])
    reference_time = float(nontriviality["reference_time"])
    target = float(nontriviality["reference_energy"])
    tolerance = float(nontriviality["reference_energy_abs_tolerance"])
    validation_lower = float(nontriviality["minimum_energy_each_validation_time"])
    validation_upper = float(nontriviality["maximum_energy_each_validation_time"])

    field = _source_field()
    proxies = {beta: _AxialCoordinateWarpProxy(field, beta) for beta in betas}
    energy_by_order: dict[int, dict[float, dict[float, float]]] = {}
    for order in orders:
        energy_by_order[order] = {
            beta: {
                time: _energy(proxies[beta], time=time, order=order)
                for time in validation_times
            }
            for beta in betas
        }
    coarse = energy_by_order[orders[0]]
    fine = energy_by_order[orders[1]]
    energy_refinement = max(
        abs(coarse[beta][time] - fine[beta][time]) / max(abs(fine[beta][time]), 1.0e-12)
        for beta in betas
        for time in validation_times
    )

    baseline_by_time = {
        time: _morphology(proxies[0.0], time=time, grid_size=morphology_grid_size)
        for time in morph_times
    }

    rows = []
    robust_rows = []
    for beta in betas:
        raw_reference = float(fine[beta][reference_time])
        baseline_parent_replay = bool(beta == 0.0 and abs(raw_reference - target) <= tolerance)
        scale = 1.0 if baseline_parent_replay else _normalization_scale(raw_reference, target)
        normalized_reference = scale * scale * raw_reference
        validation_rows = []
        validation_ok = True
        for time in validation_times:
            raw = float(fine[beta][time])
            normalized = scale * scale * raw
            passed = bool(validation_lower <= normalized <= validation_upper)
            validation_ok &= passed
            validation_rows.append(
                {
                    "time": time,
                    "raw_energy": raw,
                    "jointly_renormalized_energy": normalized,
                    "validation_energy_range_pass": passed,
                }
            )

        morphology_rows = []
        q90_gains = []
        q99_gains = []
        axial_rel = []
        radial_rel = []
        collar_ratio = []
        for time in morph_times:
            base = baseline_by_time[time]
            trial = base if beta == 0.0 else _morphology(
                proxies[beta], time=time, grid_size=morphology_grid_size
            )
            q90_gain = float(trial["axial_q90_over_Zp"] - base["axial_q90_over_Zp"])
            q99_gain = float(trial["axial_q99_over_Zp"] - base["axial_q99_over_Zp"])
            q90_gains.append(q90_gain)
            q99_gains.append(q99_gain)
            axial_rel.append(float(trial["axial_rms_over_Zp"] / base["axial_rms_over_Zp"] - 1.0))
            radial_rel.append(float(trial["radial_rms_over_Rp"] / base["radial_rms_over_Rp"] - 1.0))
            collar_ratio.append(
                float(
                    trial["physical_plateau_collar_enstrophy_fraction"]
                    / max(base["physical_plateau_collar_enstrophy_fraction"], 1.0e-300)
                )
            )
            morphology_rows.append(
                {
                    "time": time,
                    "metrics": trial,
                    "q90_gain_over_Zp": q90_gain,
                    "q99_gain_over_Zp": q99_gain,
                }
            )

        reference_ok = bool(abs(normalized_reference - target) <= tolerance)
        preflight = bool(reference_ok and validation_ok)
        robust_q90 = bool(all(value > MORPHOLOGY_ZERO_TOLERANCE for value in q90_gains))
        robust_q99 = bool(all(value > MORPHOLOGY_ZERO_TOLERANCE for value in q99_gains))
        structure = _structure_checks(field, beta=beta)
        structure_ok = bool(
            structure["outside_physical_support_max_abs_velocity"] == 0.0
            and structure["representative_core_signs_all_pass"]
            and structure["map_monotone_on_support"]
            and structure["support_endpoint_map_error"] <= 1.0e-14
        )
        row = {
            "beta": beta,
            "baseline_parent_replay_without_common_scale": baseline_parent_replay,
            "raw_reference_energy": raw_reference,
            "common_velocity_scale": scale,
            "common_scale_fractional_change": float(scale - 1.0),
            "jointly_renormalized_reference_energy": float(normalized_reference),
            "reference_energy_gate_pass": reference_ok,
            "validation_energy_range_all_pass": bool(validation_ok),
            "validation_time_energy": validation_rows,
            "representation_energy_preflight_passed": preflight,
            "structure_checks": structure,
            "structure_preflight_passed": structure_ok,
            "morphology_by_time": morphology_rows,
            "q90_gain_over_Zp_by_time": q90_gains,
            "q99_gain_over_Zp_by_time": q99_gains,
            "q90_gain_positive_at_all_morphology_times": robust_q90,
            "q99_gain_positive_at_all_morphology_times": robust_q99,
            "minimum_axial_rms_relative_change_across_morphology_times": float(min(axial_rel)),
            "maximum_axial_rms_relative_change_across_morphology_times": float(max(axial_rel)),
            "maximum_radial_rms_relative_change_across_morphology_times": float(max(radial_rel)),
            "maximum_collar_enstrophy_fraction_ratio_across_morphology_times": float(max(collar_ratio)),
            "positive_common_scale_morphology_invariance": True,
            "navier_stokes_balance_invariant_under_common_scale": False,
            "piola_divergence_identity_used": True,
        }
        if beta > 0.0 and preflight and structure_ok and robust_q90:
            robust_rows.append(row)
        rows.append(row)

    smallest_robust = robust_rows[0]["beta"] if robust_rows else None
    return {
        "task_id": TASK_ID,
        "source_candidate_sha256": str(field.sha256),
        "parent_screen_task": "CR003-AXIAL-ENVELOPE-JOINT-RENORM-SCREEN-055",
        "screen_protocol": {
            "betas": list(betas),
            "quadrature_orders": list(orders),
            "morphology_times": list(morph_times),
            "morphology_grid_size": int(morphology_grid_size),
            "physical_z_cut": PHYSICAL_Z_CUT,
            "warp_power": 4,
            "parameter_count_added_for_screen": 0,
            "new_basis_shapes_added": 0,
            "production_beta_selected": False,
        },
        "energy_quadrature_max_relative_refinement_change": float(energy_refinement),
        "baseline_morphology_by_time": baseline_by_time,
        "rows": rows,
        "smallest_sampled_preflight_passing_multi_time_q90_mover": smallest_robust,
        "smallest_mover_is_production_selection": False,
        "held_out_pde_residual": {
            "evaluated": False,
            "reason": (
                "This coordinate-warp proxy changes velocity geometry and has no governed pressure/force child. "
                "The Piola map preserves divergence structure but is not a Navier--Stokes invariance, and the "
                "subsequent common energy scale is not an NS invariance either. Fresh candidate-level pressure/force "
                "and held-out momentum checks are required after materialization."
            ),
        },
        "routing_contract": (
            "Do not add another swirl/cap/odd-q basis. This is the first larger coordinate-geometry test after the "
            "collar-heavy axial-envelope screen. If a moderate beta produces multi-time q90 gain with materially lower "
            "collar/radial cost than the alpha=1 envelope capacity result, hand at most this one Piola-warp degree to "
            "Agent 1 for an explicitly bounded, newly identified child; Agent 2/3 must then rebuild pressure/force and "
            "run fresh full momentum plus actual 3-D streamline/vorticity comparison. If q90 still requires an extreme "
            "warp, stop interior coordinate growth and move to explicit support/aspect-ratio geometry."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--morphology-grid-size", type=int, default=33)
    args = parser.parse_args(argv)
    report = audit_axial_coordinate_warp_screen(morphology_grid_size=args.morphology_grid_size)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
