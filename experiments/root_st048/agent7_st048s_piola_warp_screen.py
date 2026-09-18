"""Screen one divergence-preserving axial geometry degree on frozen ST048-S.

This is a Constrained Agent 7 representation-capacity experiment.  It does not
add a basis, fit pressure/forcing, change a scientific threshold, or promote a
candidate.  The parent is the exact frozen ST048-S recipe from PR #390, selected
by the preceding target-free morphology transfer audit as the cleaner
visualization-oriented residual challenger.

For |z|<2 define

    h_beta(z) = z * [1 - beta (1-(z/2)^2)^4]

and apply the contravariant Piola pullback

    u_beta = [h'(z) u_x(x,y,h(z),t),
              h'(z) u_y(x,y,h(z),t),
              u_z(x,y,h(z),t)].

Thus div(u_beta)=h'(z) div(u)(x,y,h(z),t) for a smooth parent.  Each nonzero
beta receives one positive common velocity scale restoring E(t=.25)=1.  That
scale is part of this capacity screen only; it is not a Navier--Stokes
invariance, so parent pressure/forcing and held-out residual evidence are not
transferred.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import numpy as np
from numpy.polynomial.legendre import leggauss

import agent7_morphology_transfer as transfer

TASK_ID = "CR003-ST048S-PIOLA-WARP-SCREEN-061"
PARENT_ID = "ST048-S"
PHYSICAL_Z_CUT = 2.0
DEFAULT_BETAS = (0.0, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20)
DEFAULT_TIMES = (0.25, 0.50, 0.75)
DEFAULT_ENERGY_TIMES = (0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75)
DEFAULT_GRID_SIZE = 33
REFERENCE_ENERGY = 1.0
REFERENCE_ENERGY_TOLERANCE = 1.0e-3
VALIDATION_ENERGY_RANGE = (0.1, 10.0)
DIVERGENCE_PREFLIGHT_MAX = 1.0e-6
MORPHOLOGY_ZERO_TOLERANCE = 1.0e-12
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "canonical_velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_spatial_basis_added": False,
    "diagnostic_coordinate_degree_added": True,
    "production_beta_selected": False,
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


def _validate_beta(beta: float) -> float:
    beta = float(beta)
    if not np.isfinite(beta) or not 0.0 <= beta <= 0.20:
        raise ValueError("beta must lie in the preregistered [0,0.20] screen")
    return beta


def _validate_betas(values: Iterable[float]) -> tuple[float, ...]:
    values = tuple(_validate_beta(v) for v in values)
    if len(values) < 2 or values[0] != 0.0:
        raise ValueError("beta grid must start at zero and contain at least two values")
    if any(a >= b for a, b in zip(values, values[1:])):
        raise ValueError("beta grid must be strictly increasing")
    return values


def _warp_z_and_jacobian(z, beta: float):
    beta = _validate_beta(beta)
    z = np.asarray(z, dtype=float)
    if not np.all(np.isfinite(z)):
        raise ValueError("z must be finite")
    s = z / PHYSICAL_Z_CUT
    inside = np.abs(s) < 1.0
    one_minus = np.maximum(1.0 - s * s, 0.0)
    bump = one_minus**4
    mapped = np.where(inside, z * (1.0 - beta * bump), z)
    jac_inside = 1.0 - beta * bump + 8.0 * beta * s * s * one_minus**3
    jac = np.where(inside, jac_inside, 1.0)
    if np.any(jac <= 0.0) or not np.all(np.isfinite(mapped)) or not np.all(np.isfinite(jac)):
        raise RuntimeError("coordinate warp lost orientation or became nonfinite")
    return mapped, jac


class _PiolaWarpField:
    def __init__(self, parent, beta: float):
        self.parent = parent
        self.beta = _validate_beta(beta)

    def at_points(self, points, time: float):
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
            raise ValueError("points must be finite shape (n,3)")
        mapped_z, jac = _warp_z_and_jacobian(points[:, 2], self.beta)
        mapped = np.array(points, copy=True)
        mapped[:, 2] = mapped_z
        velocity = np.asarray(self.parent.at_points(mapped, float(time)), dtype=float)
        if velocity.shape != points.shape or not np.all(np.isfinite(velocity)):
            raise RuntimeError("unexpected parent velocity")
        out = np.array(velocity, copy=True)
        out[:, 0] *= jac
        out[:, 1] *= jac
        return out


class _ScaledField:
    def __init__(self, parent, scale: float):
        scale = float(scale)
        if not np.isfinite(scale) or scale <= 0.0:
            raise ValueError("scale must be positive and finite")
        self.parent = parent
        self.scale = scale

    def at_points(self, points, time: float):
        return self.scale * np.asarray(self.parent.at_points(points, float(time)), dtype=float)


def _energy(field, *, time: float, order: int) -> float:
    """Axisymmetric kinetic-energy quadrature over r<2, |z|<2."""
    if order < 12:
        raise ValueError("energy quadrature order must be >=12")
    node, weight = leggauss(int(order))
    radius = node + 1.0
    z = 2.0 * node
    wr = weight
    wz = 2.0 * weight
    rr, zz = np.meshgrid(radius, z, indexing="ij")
    points = np.column_stack((rr.ravel(), np.zeros(rr.size), zz.ravel()))
    velocity = np.asarray(field.at_points(points, float(time)), dtype=float)
    speed_sq = np.sum(velocity * velocity, axis=1).reshape(rr.shape)
    # E = 1/2 int |u|^2 dV = pi int |u|^2 r dr dz for an axisymmetric field.
    return float(np.sum((wr[:, None] * wz[None, :]) * (np.pi * rr * speed_sq)))


def _normalization_scale(raw_reference_energy: float) -> float:
    raw_reference_energy = float(raw_reference_energy)
    if not np.isfinite(raw_reference_energy) or raw_reference_energy <= 0.0:
        raise ValueError("raw reference energy must be positive")
    return float(np.sqrt(REFERENCE_ENERGY / raw_reference_energy))


def _structure_checks(parent, warped, *, beta: float, scale: float, times: tuple[float, ...]):
    scaled = _ScaledField(warped, scale)
    exterior = np.array(
        [[2.01, 0.0, 0.0], [-2.01, 0.0, 0.0], [0.4, 0.2, 2.01], [0.4, -0.2, -2.01]],
        dtype=float,
    )
    support_max = max(float(np.max(np.abs(scaled.at_points(exterior, t)))) for t in times)

    sign_checks = []
    for time in times:
        for radius in (0.05, 0.10, 0.20):
            for z in (-0.10, 0.10):
                ux, uy, uz = scaled.at_points(np.array([[radius, 0.0, z]]), time)[0]
                sign_checks.append(bool(ux < 0.0 and uy > 0.0 and z * uz > 0.0))

    points = np.array(
        [[0.31, 0.22, 0.17], [0.62, -0.27, -0.41], [0.88, 0.31, 0.67], [1.22, -0.36, -1.22]],
        dtype=float,
    )
    h = 1.0e-5
    divergence = np.zeros(len(points), dtype=float)
    for axis in range(3):
        offset = np.zeros_like(points)
        offset[:, axis] = h
        plus = scaled.at_points(points + offset, 0.5)
        minus = scaled.at_points(points - offset, 0.5)
        divergence += (plus[:, axis] - minus[:, axis]) / (2.0 * h)

    z_grid = np.linspace(-PHYSICAL_Z_CUT, PHYSICAL_Z_CUT, 4001)
    mapped, jac = _warp_z_and_jacobian(z_grid, beta)
    replay_points = np.array(
        [[0.1, 0.0, 0.1], [0.5, -0.2, -0.7], [1.2, 0.3, 1.1], [1.8, 0.0, -1.7]],
        dtype=float,
    )
    beta_zero_replay = 0.0
    if beta == 0.0 and abs(scale - 1.0) <= 1.0e-14:
        beta_zero_replay = max(
            float(np.max(np.abs(scaled.at_points(replay_points, t) - parent.at_points(replay_points, t))))
            for t in times
        )
    return {
        "outside_support_max_abs_velocity": support_max,
        "representative_inward_swirl_bipolar_signs_all_pass": bool(all(sign_checks)),
        "representative_sign_checks": len(sign_checks),
        "cartesian_fd_divergence_max_abs": float(np.max(np.abs(divergence))),
        "minimum_coordinate_jacobian": float(np.min(jac)),
        "maximum_coordinate_jacobian": float(np.max(jac)),
        "map_monotone_on_support": bool(np.all(np.diff(mapped) > 0.0)),
        "support_endpoint_map_error": float(max(abs(mapped[0] + 2.0), abs(mapped[-1] - 2.0))),
        "beta_zero_parent_replay_max_abs_velocity_error": beta_zero_replay,
    }


def _delta(child: dict[str, float], parent: dict[str, float]):
    return {
        key + "_delta": float(child[key] - parent[key])
        for key in parent
    } | {
        key + "_relative": float(child[key] / max(abs(parent[key]), 1.0e-300) - 1.0)
        for key in parent
    }


def audit_st048s_piola_warp_screen(
    *,
    betas: Iterable[float] = DEFAULT_BETAS,
    times: Iterable[float] = DEFAULT_TIMES,
    energy_times: Iterable[float] = DEFAULT_ENERGY_TIMES,
    energy_orders: tuple[int, int] = (48, 64),
    grid_size: int = DEFAULT_GRID_SIZE,
):
    betas = _validate_betas(betas)
    times = tuple(float(t) for t in times)
    energy_times = tuple(float(t) for t in energy_times)
    if not times or any(t < 0.25 or t > 0.75 or not np.isfinite(t) for t in times):
        raise ValueError("morphology times must lie in [0.25,0.75]")
    if not energy_times or any(t < 0.25 or t > 0.75 or not np.isfinite(t) for t in energy_times):
        raise ValueError("energy times must lie in [0.25,0.75]")
    if len(energy_orders) != 2 or not 12 <= energy_orders[0] < energy_orders[1]:
        raise ValueError("energy_orders must be two increasing orders >=12")
    if grid_size < 17 or grid_size % 2 == 0:
        raise ValueError("grid_size must be odd and >=17")

    _parent, children = transfer._load_frozen_fields()
    parent = children[PARENT_ID]
    parent_metrics = {t: transfer._morphology(parent, time=t, grid_size=grid_size) for t in times}

    raw_energy = {}
    energy_refinement_max = 0.0
    for beta in betas:
        warped = _PiolaWarpField(parent, beta)
        raw_energy[beta] = {}
        for time in energy_times:
            coarse = _energy(warped, time=time, order=energy_orders[0])
            fine = _energy(warped, time=time, order=energy_orders[1])
            raw_energy[beta][time] = fine
            energy_refinement_max = max(
                energy_refinement_max,
                abs(coarse - fine) / max(abs(fine), 1.0e-12),
            )

    rows = []
    clean_candidates = []
    for beta in betas:
        warped = _PiolaWarpField(parent, beta)
        reference_raw = raw_energy[beta][0.25]
        scale = 1.0 if beta == 0.0 else _normalization_scale(reference_raw)
        scaled = _ScaledField(warped, scale)
        normalized_energy = {
            time: scale * scale * raw_energy[beta][time]
            for time in energy_times
        }
        reference_pass = abs(normalized_energy[0.25] - REFERENCE_ENERGY) <= REFERENCE_ENERGY_TOLERANCE
        energy_range_pass = all(VALIDATION_ENERGY_RANGE[0] <= value <= VALIDATION_ENERGY_RANGE[1]
                                for value in normalized_energy.values())
        metrics = {t: transfer._morphology(scaled, time=t, grid_size=grid_size) for t in times}
        deltas = {t: _delta(metrics[t], parent_metrics[t]) for t in times}
        axial = [deltas[t]["axial_rms_over_support_relative"] for t in times]
        radial = [deltas[t]["radial_rms_over_support_relative"] for t in times]
        q90 = [deltas[t]["axial_q90_over_support_delta"] for t in times]
        q99 = [deltas[t]["axial_q99_over_support_delta"] for t in times]
        collar = [metrics[t]["outer_support_collar_enstrophy_fraction"] /
                  max(parent_metrics[t]["outer_support_collar_enstrophy_fraction"], 1.0e-300)
                  for t in times]
        tail065 = [metrics[t]["outer_065_enstrophy_fraction"] /
                   max(parent_metrics[t]["outer_065_enstrophy_fraction"], 1.0e-300) - 1.0
                   for t in times]
        tail075 = [metrics[t]["outer_075_enstrophy_fraction"] /
                   max(parent_metrics[t]["outer_075_enstrophy_fraction"], 1.0e-300) - 1.0
                   for t in times]
        structure = _structure_checks(parent, warped, beta=beta, scale=scale, times=times)
        preflight = bool(
            reference_pass
            and energy_range_pass
            and structure["outside_support_max_abs_velocity"] <= 1.0e-12
            and structure["representative_inward_swirl_bipolar_signs_all_pass"]
            and structure["cartesian_fd_divergence_max_abs"] <= DIVERGENCE_PREFLIGHT_MAX
            and structure["map_monotone_on_support"]
        )
        extension = bool(
            all(value >= -MORPHOLOGY_ZERO_TOLERANCE for value in axial)
            and any(value > MORPHOLOGY_ZERO_TOLERANCE for value in q90 + q99)
        )
        low_cost = bool(max(radial) <= 0.05 and max(collar) <= 2.0)
        clean = bool(beta > 0.0 and preflight and extension and low_cost)
        if clean:
            clean_candidates.append(beta)
        rows.append({
            "beta": beta,
            "parameter_count": 0 if beta == 0.0 else 1,
            "raw_reference_energy": reference_raw,
            "common_energy_normalization_scale": scale,
            "normalized_energy_by_time": normalized_energy,
            "reference_energy_pass": bool(reference_pass),
            "validation_energy_range_pass": bool(energy_range_pass),
            "morphology_by_time": metrics,
            "child_minus_parent_by_time": deltas,
            "axial_rms_relative_changes": axial,
            "radial_rms_relative_changes": radial,
            "q90_deltas_over_support": q90,
            "q99_deltas_over_support": q99,
            "outer_065_enstrophy_relative_changes": tail065,
            "outer_075_enstrophy_relative_changes": tail075,
            "collar_enstrophy_ratios": collar,
            "structure": structure,
            "preflight_pass": preflight,
            "reused_clean_extension_rule_pass": clean,
        })

    smallest_clean = min(clean_candidates) if clean_candidates else None
    first_nonzero = rows[1]
    local_sensitivity = {
        "beta": first_nonzero["beta"],
        "axial_rms_relative_derivative_estimates": [v / first_nonzero["beta"] for v in first_nonzero["axial_rms_relative_changes"]],
        "radial_rms_relative_derivative_estimates": [v / first_nonzero["beta"] for v in first_nonzero["radial_rms_relative_changes"]],
    }
    if smallest_clean is None:
        routing = (
            "No beta<=0.20 in the preregistered one-degree Piola screen satisfies the reused clean-extension rule on frozen ST048-S. "
            "Do not widen the warp, add another basis, or alter pressure/forcing on this evidence. Keep ST048-S as the visualization-oriented residual backbone and consume the fixed 3-D render/grid lane first; only revisit geometry if that direct visualization shows a material axial-reach defect."
        )
    else:
        routing = (
            f"The smallest preregistered clean Piola crossing on frozen ST048-S is beta={smallest_clean:.3f}. "
            "This is capacity evidence, not a production coefficient. Do not add another spatial basis. Agent 1 may materialize at most this single bounded geometry degree with a new candidate identity and explicit energy normalization; Agent 2/3 must then rebuild compatible pressure/forcing and run fresh held-out full momentum before any promotion."
        )

    recipes = json.loads((Path(__file__).resolve().parent / "recipes.json").read_text())
    return {
        "task_id": TASK_ID,
        "issue": 405,
        "parent_id": PARENT_ID,
        "parent_raw_candidate_sha256": recipes[PARENT_ID].get("original_candidate_sha256"),
        "contract": {
            "betas": list(betas),
            "times": list(times),
            "energy_times": list(energy_times),
            "energy_orders": list(energy_orders),
            "grid_size": grid_size,
            "reference_energy": REFERENCE_ENERGY,
            "reference_energy_abs_tolerance": REFERENCE_ENERGY_TOLERANCE,
            "validation_energy_range": list(VALIDATION_ENERGY_RANGE),
            "clean_extension_rule_reused_from_pr_398": True,
            "classification": "autonomous target-free representation-capacity diagnostic",
        },
        "upstream_st048_residual_receipt": transfer._residual_receipt(),
        "parent_morphology_by_time": parent_metrics,
        "rows": rows,
        "energy_quadrature_max_relative_refinement_change": energy_refinement_max,
        "local_beta_sensitivity": local_sensitivity,
        "smallest_preregistered_clean_beta": smallest_clean,
        "held_out_full_momentum": "not evaluated: the Piola velocity warp has no governed transformed pressure/forcing child, so ST048-S PDE evidence is not transferable",
        "routing": routing,
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="../../artifacts/st048_agent7_piola_warp_screen/report.json")
    args = parser.parse_args()
    report = audit_st048s_piola_warp_screen()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "task_id": report["task_id"],
        "smallest_preregistered_clean_beta": report["smallest_preregistered_clean_beta"],
        "energy_quadrature_max_relative_refinement_change": report["energy_quadrature_max_relative_refinement_change"],
        "routing": report["routing"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
