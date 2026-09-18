"""Target-free fixed-support axial-envelope flattening capacity audit.

The existing cap-localized correction can move the far q99 tail, but the main
q90/high-enstrophy core remains too short.  This Agent-7 increment therefore
stops basis growth and probes one minimal *base-geometry* direction instead.

The current Eq45 compact profile uses the axial envelope

    B_0(s) = (1-s^2)^p_+,  s = eta/eta_cut, p>=5.

This audit introduces one diagnostic interpolation, with exactly the same
support boundary and the same C4 vanishing order,

    g_alpha(s) = (1-alpha)s^2 + alpha s^4,
    B_alpha(s) = (1-g_alpha(s))^p_+.

alpha=0 exactly replays the current field.  alpha=1 is the fixed quartic
flattening probe: it leaves the support boundary at |eta|=eta_cut while making
the interior envelope less center-weighted.  No public image, residual map, or
fitted shape parameter selects alpha=1.

This module is representation/visualization-capacity evidence only.  It does
not materialize a new candidate, select a production alpha/bound, change
forcing/pressure, evaluate a held-out PDE residual, or establish OpenAI visual
correspondence.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_bipolar_axial_cap_poloidal_capacity import _cap_response
from .constrained_eq45_bipolar_axial_shoulder_poloidal_capacity import (
    _morphology_vector,
    _vorticity_morphology,
)
from .constrained_eq45_bipolar_compact_poloidal_capacity import (
    COMPACT_POLOIDAL,
    PHI01,
    PHI03,
    _cylindrical_components,
    _novelty,
    _probe_groups,
    _rank_condition,
    _response,
    _source_field,
)
from .constrained_eq45_poloidal_streamfunction import eq45_streamfunction_poloidal_jet
from .constrained_eq45_profile_basis import Eq45CompactProfileBasis, Eq45ProfileJets
from .constrained_eq45_velocity import Eq45VelocityBackbone

TASK_ID = "CR003-AXIAL-ENVELOPE-FLATTENING-CAPACITY-053"
DEFAULT_TRIAL_ALPHA = 1.0
ALPHA_LOCAL_MIN = -0.25
ALPHA_LOCAL_MAX = 1.0
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "canonical_velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_spatial_basis_added": False,
    "diagnostic_geometry_family_added": True,
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


def _validate_alpha(alpha: float) -> float:
    alpha = float(alpha)
    if not np.isfinite(alpha) or not ALPHA_LOCAL_MIN <= alpha <= ALPHA_LOCAL_MAX:
        raise ValueError(
            f"alpha must lie in [{ALPHA_LOCAL_MIN},{ALPHA_LOCAL_MAX}] for this diagnostic"
        )
    return alpha


def _flattened_envelope_jets(
    basis: Eq45CompactProfileBasis,
    x,
    eta,
    alpha: float,
) -> Eq45ProfileJets:
    """Evaluate the existing Phi/F coefficients with only the axial envelope changed."""
    if not isinstance(basis, Eq45CompactProfileBasis):
        raise TypeError("basis must be Eq45CompactProfileBasis")
    alpha = _validate_alpha(alpha)
    x_arr = np.asarray(x, dtype=float)
    eta_arr = np.asarray(eta, dtype=float)
    try:
        x_arr, eta_arr = np.broadcast_arrays(x_arr, eta_arr)
    except ValueError as exc:
        raise ValueError("X and eta must be broadcast-compatible") from exc
    if not np.all(np.isfinite(x_arr)) or not np.all(np.isfinite(eta_arr)):
        raise ValueError("X and eta must be finite")
    if np.any(x_arr < 0.0):
        raise ValueError("X must be nonnegative")

    # Reuse the governed radial envelope exactly; replace only the eta envelope.
    bx, dbx, d2bx, _, _ = basis._envelopes(x_arr, eta_arr)
    sx = x_arr / basis.x_cut
    se = eta_arr / basis.eta_cut
    inside_eta = np.abs(se) < 1.0
    blended = (1.0 - alpha) * se**2 + alpha * se**4
    one_minus = np.where(inside_eta, 1.0 - blended, 0.0)
    # The guarded alpha interval keeps one_minus nonnegative for |se|<1.
    if np.any(one_minus < -64.0 * np.finfo(float).eps):
        raise RuntimeError("diagnostic axial envelope became negative inside support")
    one_minus = np.maximum(one_minus, 0.0)
    power = basis.cutoff_power
    be = one_minus**power
    dblended_deta = (
        2.0 * (1.0 - alpha) * se + 4.0 * alpha * se**3
    ) / basis.eta_cut
    dbe = np.where(
        inside_eta,
        -power * one_minus ** (power - 1) * dblended_deta,
        0.0,
    )

    phi = np.zeros_like(x_arr, dtype=float)
    phi_x = np.zeros_like(x_arr, dtype=float)
    phi_eta = np.zeros_like(x_arr, dtype=float)
    phi_xx = np.zeros_like(x_arr, dtype=float)
    phi_xeta = np.zeros_like(x_arr, dtype=float)
    swirl = np.zeros_like(x_arr, dtype=float)

    for index, (i, j) in enumerate(basis.mode_indices):
        sx_i = sx**i
        se_j = se**j
        mono_x = np.zeros_like(x_arr) if i == 0 else (i / basis.x_cut) * sx ** (i - 1)
        mono_xx = (
            np.zeros_like(x_arr)
            if i <= 1
            else (i * (i - 1) / basis.x_cut**2) * sx ** (i - 2)
        )
        mono_eta = (
            np.zeros_like(x_arr)
            if j == 0
            else (j / basis.eta_cut) * se ** (j - 1)
        )
        x_part = bx * sx_i
        x_part_x = dbx * sx_i + bx * mono_x
        x_part_xx = d2bx * sx_i + 2.0 * dbx * mono_x + bx * mono_xx
        eta_part = be * se_j
        eta_part_eta = dbe * se_j + be * mono_eta

        basis_value = x_part * eta_part
        basis_x = x_part_x * eta_part
        basis_eta = x_part * eta_part_eta
        basis_xx = x_part_xx * eta_part
        basis_xeta = x_part_x * eta_part_eta
        phi += basis.phi_coefficients[index] * basis_value
        phi_x += basis.phi_coefficients[index] * basis_x
        phi_eta += basis.phi_coefficients[index] * basis_eta
        phi_xx += basis.phi_coefficients[index] * basis_xx
        phi_xeta += basis.phi_coefficients[index] * basis_xeta
        swirl += basis.swirl_coefficients[index] * basis_value

    values = (phi, phi_x, phi_eta, phi_xx, phi_xeta, swirl)
    if not all(np.all(np.isfinite(value)) for value in values):
        raise RuntimeError("flattened axial-envelope profile jets became nonfinite")
    return Eq45ProfileJets(
        phi=phi,
        phi_x=phi_x,
        phi_eta=phi_eta,
        phi_xx=phi_xx,
        phi_xeta=phi_xeta,
        swirl=swirl,
    )


class _EnvelopeFieldProxy:
    """Evaluate one diagnostic envelope alpha without creating a candidate identity."""

    def __init__(self, field, alpha: float):
        self._field = field
        self.alpha = _validate_alpha(alpha)

    def __getattr__(self, name):
        return getattr(self._field, name)

    def at_points(self, points, time: float):
        points_arr = np.asarray(points, dtype=float)
        parent = self._field.parent
        time_arr = parent._validated_time(points_arr, time)
        backbone = Eq45VelocityBackbone(parent.profile_values, h=parent.h)
        coordinates = backbone.coordinates(points_arr, time_arr)
        jets = _flattened_envelope_jets(
            parent.profile_basis, coordinates.X, coordinates.eta, self.alpha
        )
        poloidal = eq45_streamfunction_poloidal_jet(
            coordinates.X,
            coordinates.eta,
            parent.h,
            jets.phi,
            jets.phi_x,
            jets.phi_eta,
            jets.phi_xx,
            jets.phi_xeta,
        )
        x = points_arr[..., 0]
        y = points_arr[..., 1]
        radius = np.hypot(x, y)
        q = coordinates.q
        psi = np.power(q, 0.5 - parent.h) * coordinates.X * jets.phi
        u_r = radius * poloidal.v0 / (2.0 * q)
        u_theta = radius * np.power(q, -1.0 - parent.h) * jets.swirl
        u_z = np.power(q, -0.5 - parent.h) * poloidal.U
        return self._field.taper.apply_cylindrical(
            points_arr, psi=psi, u_r=u_r, u_theta=u_theta, u_z=u_z
        )



def _envelope_response(field, step: float, points: np.ndarray, time: float) -> np.ndarray:
    step = float(step)
    if not np.isfinite(step) or not 0.0 < step <= 0.20:
        raise ValueError("step must lie in (0,0.20]")
    plus = _EnvelopeFieldProxy(field, step).at_points(points, float(time))
    minus = _EnvelopeFieldProxy(field, -step).at_points(points, float(time))
    return (plus - minus) / (2.0 * step)


def _response_matrix(field, *, time: float, step: float) -> tuple[np.ndarray, np.ndarray]:
    points = np.concatenate(tuple(_probe_groups().values()))
    columns = [
        _response(field, PHI01, step, points, time).reshape(-1),
        _response(field, PHI03, step, points, time).reshape(-1),
        _response(field, COMPACT_POLOIDAL, step, points, time).reshape(-1),
        _cap_response(field, step, points, time).reshape(-1),
        _envelope_response(field, step, points, time).reshape(-1),
    ]
    return points, np.column_stack(columns)


def _morphology(field, *, time: float, alpha: float, grid_size: int) -> dict[str, float]:
    return _vorticity_morphology(
        _EnvelopeFieldProxy(field, alpha), time=float(time), grid_size=grid_size
    )


def _morphology_delta(base: dict[str, float], trial: dict[str, float]) -> dict[str, float]:
    return {f"{key}_delta": float(trial[key] - base[key]) for key in base}


def _structure_checks(field, *, trial_alpha: float) -> dict[str, float]:
    # alpha=0 must be a literal replay of the current supported field.
    replay_points = np.concatenate(tuple(_probe_groups().values()))
    replay_error = max(
        float(
            np.max(
                np.abs(
                    _EnvelopeFieldProxy(field, 0.0).at_points(replay_points, t)
                    - field.at_points(replay_points, t)
                )
            )
        )
        for t in (0.25, 0.50, 0.75)
    )

    # Physical support remains exactly the same because the existing taper is reused.
    support_points = np.array(
        [
            [2.01, 0.0, 0.0],
            [-2.01, 0.0, 0.0],
            [0.4, 0.2, 2.01],
            [0.4, -0.2, -2.01],
        ],
        dtype=float,
    )
    support_max = float(
        np.max(np.abs(_EnvelopeFieldProxy(field, trial_alpha).at_points(support_points, 0.5)))
    )

    # Independent Cartesian finite-difference divergence check away from the axis.
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
    proxy = _EnvelopeFieldProxy(field, trial_alpha)
    divergence = np.zeros(len(check), dtype=float)
    for axis in range(3):
        offset = np.zeros_like(check)
        offset[:, axis] = h
        plus = proxy.at_points(check + offset, 0.5)
        minus = proxy.at_points(check - offset, 0.5)
        divergence += (plus[:, axis] - minus[:, axis]) / (2.0 * h)

    # Record whether the center-sign convention survives at representative times.
    core_rows = []
    for time in np.linspace(field.time_start, field.time_end, 6):
        tau = 1.0 - float(time)
        point = np.array([[0.1 * np.sqrt(tau), 0.0, 0.1 * tau**0.495]])
        velocity = proxy.at_points(point, float(time))[0]
        ur, ut, uz = _cylindrical_components(point, velocity.reshape(1, 3))
        core_rows.append(bool(ur[0] < 0.0 and ut[0] > 0.0 and uz[0] > 0.0))

    return {
        "alpha_zero_parent_replay_max_abs_velocity_error": replay_error,
        "outside_physical_support_max_abs_velocity": support_max,
        "cartesian_fd_divergence_max_abs": float(np.max(np.abs(divergence))),
        "representative_core_signs_all_pass": bool(all(core_rows)),
    }


def audit_axial_envelope_flattening_capacity(
    *,
    time: float = 0.5,
    alpha_steps: Iterable[float] = (0.02, 0.01),
    trial_alpha: float = DEFAULT_TRIAL_ALPHA,
    morphology_grid_size: int = 33,
) -> dict[str, Any]:
    steps = tuple(float(v) for v in alpha_steps)
    if len(steps) != 2 or not steps[0] > steps[1] > 0.0 or steps[0] > 0.20:
        raise ValueError("alpha_steps must contain two decreasing values in (0,0.20]")
    trial_alpha = _validate_alpha(trial_alpha)
    if morphology_grid_size < 17 or morphology_grid_size % 2 == 0:
        raise ValueError("morphology_grid_size must be odd and >=17")

    field = _source_field()
    points, coarse = _response_matrix(field, time=float(time), step=steps[0])
    _, fine = _response_matrix(field, time=float(time), step=steps[1])
    response_refinement = float(
        np.linalg.norm(coarse - fine) / max(np.linalg.norm(fine), 1e-300)
    )
    response_rank = _rank_condition(fine)
    envelope_novelty = _novelty(fine[:, :4], fine[:, 4])
    envelope_response_rms = float(np.linalg.norm(fine[:, 4]) / np.sqrt(fine.shape[0]))

    base_morph = _morphology(field, time=float(time), alpha=0.0, grid_size=morphology_grid_size)
    trial_morph = _morphology(
        field, time=float(time), alpha=trial_alpha, grid_size=morphology_grid_size
    )
    morph_plus = _morphology(
        field, time=float(time), alpha=steps[-1], grid_size=morphology_grid_size
    )
    morph_minus = _morphology(
        field, time=float(time), alpha=-steps[-1], grid_size=morphology_grid_size
    )
    local_morphology_response = (
        _morphology_vector(morph_plus) - _morphology_vector(morph_minus)
    ) / (2.0 * steps[-1])

    time_rows = []
    for sample_time in (field.time_start, 0.5 * (field.time_start + field.time_end), field.time_end):
        baseline = _morphology(
            field, time=float(sample_time), alpha=0.0, grid_size=morphology_grid_size
        )
        trial = _morphology(
            field, time=float(sample_time), alpha=trial_alpha, grid_size=morphology_grid_size
        )
        time_rows.append(
            {
                "time": float(sample_time),
                "baseline": baseline,
                "trial": trial,
                "delta": _morphology_delta(baseline, trial),
            }
        )

    return {
        "task_id": TASK_ID,
        "source_candidate_sha256": str(field.sha256),
        "time": float(time),
        "geometry_definition": {
            "baseline": "B0(s)=(1-s^2)^p_+",
            "diagnostic_family": "B_alpha(s)=(1-[(1-alpha)s^2+alpha*s^4])^p_+",
            "support_boundary_unchanged": True,
            "cutoff_power_unchanged": int(field.parent.profile_basis.cutoff_power),
            "eta_cut_unchanged": float(field.parent.profile_basis.eta_cut),
            "trial_alpha": trial_alpha,
            "trial_alpha_interpretation": "alpha=1 is fixed quartic interior flattening, not a fitted production value",
        },
        "representation_increment": {
            "new_spatial_basis_shapes_added": 0,
            "diagnostic_geometry_parameters_added": 1,
            "new_production_fit_parameters_added": 0,
            "production_parameter_selected": False,
        },
        "velocity_space_diagnostics": {
            "modes": ["Phi01", "Phi03", COMPACT_POLOIDAL, "AXIAL_CAP_BANDED_C4_ODD_Z_POLOIDAL", "AXIAL_ENVELOPE_ALPHA"],
            "probe_point_count": int(len(points)),
            "rank_condition": response_rank,
            "finite_difference_refinement_relative_change": response_refinement,
            "envelope_novelty_outside_existing_four_span": envelope_novelty,
            "envelope_response_rms_per_public_velocity_component": envelope_response_rms,
        },
        "local_morphology_response_at_alpha_zero": {
            "alpha_step": steps[-1],
            "metric_order": [
                "radial_rms_over_Rp",
                "axial_rms_over_Zp",
                "aspect_z_over_r",
                "outer_half_enstrophy_fraction",
                "outer_065_enstrophy_fraction",
                "outer_075_enstrophy_fraction",
            ],
            "derivative": [float(v) for v in local_morphology_response],
        },
        "baseline_morphology": base_morph,
        "trial_morphology": trial_morph,
        "trial_delta": _morphology_delta(base_morph, trial_morph),
        "time_slice_morphology": time_rows,
        "structure_checks": _structure_checks(field, trial_alpha=trial_alpha),
        "routing_contract": (
            "This screen asks whether the existing quadratic compact eta envelope, rather than missing basis rank, "
            "is obstructing q90/high-enstrophy axial reach. If fixed-support quartic flattening materially moves q90 "
            "with acceptable velocity conditioning and structural guards, Agent 1 should materialize at most one "
            "explicitly bounded envelope-shape child and rerun fresh energy/pressure/force/full-momentum validation. "
            "If q90 still does not move materially, stop compact-envelope shape growth and investigate a larger-scale "
            "support/coordinate geometry change instead."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--trial-alpha", type=float, default=DEFAULT_TRIAL_ALPHA)
    parser.add_argument("--morphology-grid-size", type=int, default=33)
    args = parser.parse_args(argv)
    report = audit_axial_envelope_flattening_capacity(
        trial_alpha=args.trial_alpha,
        morphology_grid_size=args.morphology_grid_size,
    )
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
