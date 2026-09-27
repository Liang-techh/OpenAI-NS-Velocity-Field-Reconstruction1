"""Mean compatibility replay for the locked co-designed wave candidate.

The saved constrained mean is augmented by only the selected mode-0 tangent
controls from ``wave_dynamics_codesign.json``.  The mode-0 velocity starts at
zero at the reference time, while its physical-time derivative and constant
mode-0 pressure gradient remain in the mean momentum.  The nonzero mode-1
wave contributes its complete angular-mean advection through
``single_mode_force``.  Mode-1 and mode-2 tangent controls have zero angular
mean at the instantaneous initial state and are intentionally omitted here.

This is a compatibility diagnostic for one locked candidate.  It does not
optimize, integrate in time, or validate the full oscillatory PDE.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from affine_momentum import jets, momentum  # noqa: E402
from broad_meridional_constrained import (  # noqa: E402
    load_saved_field as load_constrained_mean,
)
from broad_shear_dynamic_control import load_saved_field as load_dynamic_field  # noqa: E402
from full_wave_tangent import _unpack_full  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from joined_field import coordinates  # noqa: E402
from midplane_integrated_moment_balance import evaluate as replay_moments  # noqa: E402
from supported_fourier_basis import basis_data  # noqa: E402
from wave_mean_flux import single_mode_force  # noqa: E402
from joint_wave_mean_fit import (  # noqa: E402
    _grouped_cone_panels,
    _radial_panels,
    _wave_force_cartesian,
    _wave_metadata,
    _wave_moment_column,
)


MEAN_PATH = ROOT / "broad_meridional_constrained.json"
DYNAMIC_PATH = ROOT / "broad_shear_dynamic_control.json"
CANDIDATE_PATH = ROOT / "wave_dynamics_codesign.json"
REPLAY_PATH = ROOT / "wave_dynamics_replay.json"
MEAN_REPLAY_PATH = ROOT / "broad_meridional_constrained_replay.json"
OUTPUT_PATH = ROOT / "wave_dynamics_mean_compatibility.json"

MOMENT_ORDER = 96
MOMENT_ZFACTOR = 0.002
CONE_ORDER = 64
MARGIN = np.array([0.02, 0.0, 0.0], dtype=float)


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _save(report, path=OUTPUT_PATH):
    Path(path).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _pack(value):
    return np.asarray(value).tolist()


class MeanMode0Field:
    """Saved constrained mean plus an affine mode-0 tangent correction."""

    def __init__(self, mean, center, widths, degree, derivative, pressure, tau0):
        self.mean = mean
        self.inner = mean.inner
        self.nu = float(mean.nu)
        self.center = np.asarray(center, dtype=float)
        self.widths = np.asarray(widths, dtype=float)
        self.degree = int(degree)
        self.derivative = np.asarray(derivative, dtype=complex)
        self.pressure = np.asarray(pressure, dtype=complex)
        self.tau0 = float(tau0)
        self.join_X = getattr(mean, "join_X", None)
        self.ratio = getattr(mean, "ratio", None)

    def fields(self, points, tau):
        points = np.asarray(points, dtype=float)
        tau_value = float(np.asarray(tau, dtype=float).ravel()[0])
        velocity, pressure = self.mean.fields(points, tau_value)
        velocity = velocity.copy()
        pressure = pressure.copy()
        mode_velocity, mode_pressure, _ = basis_data(
            points, self.center, self.widths, 0, self.degree, (0.0, 0.0)
        )
        # t = -tau, hence t-t0 = tau0-tau.  The pressure controls are
        # constant in time, as in full_wave_tangent.LocalPotentialField.
        delta_t = self.tau0 - tau_value
        velocity += delta_t * np.einsum(
            "niq,q->ni", mode_velocity, self.derivative
        ).real
        pressure += np.einsum("nq,q->n", mode_pressure, self.pressure).real
        return velocity, pressure


def _wave_forcing(points, wave, nu, hspace):
    """Complete angular-mean mode-1 force in Cartesian coordinates."""

    return _wave_force_cartesian(points, wave, nu, hspace)


def _radial_integrals(inner, R, z, tau, breaks, order, center, widths):
    """Mode-0 velocity/pressure integrals for reusable moment rows."""

    Xmax = float(inner.p.X_max)
    q = float(coordinates(0.0, z / np.sqrt(inner.nu), tau, inner.h)["q"])
    ri = np.sqrt(2.0 * inner.nu * q * Xmax)
    edges = [0.0, ri, R]
    edges.extend(ri + (R - ri) * float(value) for value in breaks)
    edges = sorted(set(np.clip(edges, 0.0, R)))
    nodes, node_weights = leggauss(int(order))
    radii = np.concatenate([
        0.5 * (lo + hi) + 0.5 * (hi - lo) * nodes
        for lo, hi in zip(edges[:-1], edges[1:])
    ])
    weights = np.concatenate([
        0.5 * (hi - lo) * node_weights
        for lo, hi in zip(edges[:-1], edges[1:])
    ])
    points = np.column_stack((radii, np.zeros_like(radii), np.full_like(radii, z)))
    velocity, pressure, _ = basis_data(
        points, center, widths, 0, 2, (0.0, 0.0)
    )
    velocity = velocity.real
    pressure = pressure.real
    return {
        "velocity_theta": np.einsum("n,n,nq->q", weights, radii**2, velocity[:, 1, :]),
        "velocity_axial": np.einsum("n,n,nq->q", weights, radii, velocity[:, 2, :]),
        "pressure": np.einsum("n,n,nq->q", weights, radii, pressure),
    }


def _mode0_moment_rows(inner, tau, breaks, center, widths, order=MOMENT_ORDER,
                       zfactor=MOMENT_ZFACTOR, htime=None):
    """Return the 4x180 linear response of integrated moments to tangent x."""

    if htime is None:
        htime = 1.0e-4 * tau
    response = np.zeros((4, 180), dtype=float)
    metadata = []
    for row_index, eta in enumerate((-0.2, 0.2)):
        outer = inner.from_similarity(
            inner.p.X_max * np.array([16.0**2]), [eta], tau
        )[0]
        R, z = float(outer[0]), float(outer[2])
        zp = inner.from_similarity(
            inner.p.X_max * np.array([1.0]), [eta + 0.01], tau
        )[0, 2]
        zm = inner.from_similarity(
            inner.p.X_max * np.array([1.0]), [eta - 0.01], tau
        )[0, 2]
        hz = zfactor * abs(zp - zm) / 0.02
        tau_offsets = (-2, -1, 1, 2)
        theta_states = []
        axial_states = []
        for offset in tau_offsets:
            tau_j = tau + offset * htime
            state = _radial_integrals(
                inner, R, z, tau_j, breaks, order, center, widths
            )
            # The mode-0 coefficient is tau0-tau_j = -offset*htime.
            theta_states.append((-offset * htime) * state["velocity_theta"])
            axial_states.append((-offset * htime) * state["velocity_axial"])
        dtau_theta = (
            theta_states[0] - 8.0 * theta_states[1]
            + 8.0 * theta_states[2] - theta_states[3]
        ) / (12.0 * htime)
        dtau_axial = (
            axial_states[0] - 8.0 * axial_states[1]
            + 8.0 * axial_states[2] - axial_states[3]
        ) / (12.0 * htime)
        z_states = []
        for offset in (-2, -1, 1, 2):
            state = _radial_integrals(
                inner, R, z + offset * hz, tau, breaks, order, center, widths
            )
            z_states.append(state["pressure"])
        dz_pressure = (
            z_states[0] - 8.0 * z_states[1]
            + 8.0 * z_states[2] - z_states[3]
        ) / (12.0 * hz)
        response[row_index * 2, :27] = dtau_theta / R**2
        response[row_index * 2 + 1, :27] = dtau_axial / R
        response[row_index * 2 + 1, 27:36] = -dz_pressure / R
        metadata.append({
            "eta": float(eta),
            "outer_radius": R,
            "z": z,
            "z_step": hz,
            "time_step": htime,
            "row_order": ["theta", "axial"],
        })
    return response, metadata


def _wave_moment_forms(inner, tau, k, breaks, wave, nu, hspace, order=MOMENT_ORDER):
    """Build symmetric 4x54x54 wave-moment quadratic forms."""

    panels = _radial_panels(inner, k, (-0.2, 0.2), breaks, order)
    forms = np.zeros((4, 54, 54), dtype=float)
    for eta_index, panel in enumerate(panels):
        points = np.column_stack((
            panel["radii"],
            np.zeros(len(panel["radii"])),
            np.full(len(panel["radii"]), panel["z"]),
        ))
        from fourier_patch_evolution import basis_jets

        velocity, gradient, _, _, _ = basis_jets(
            points, wave["center"], wave["widths"], wave["mode"],
            wave["degree"], wave["carrier"], nu, hspace,
        )
        # Complex columns for x=[Re(c), Im(c)] are V and iV respectively.
        velocity_real = np.concatenate((velocity, 1j * velocity), axis=-1)
        gradient_real = np.concatenate((gradient, 1j * gradient), axis=-1)
        bilinear = 0.5 * np.einsum(
            "ncdi,ndj->ncij", gradient_real, velocity_real.conj()
        ).real
        radial = panel["radii"]
        weights = panel["weights"]
        R = panel["outer_radius"]
        # The integrated moment convention is -integral(r^2 R_theta)/R^2
        # and -integral(r R_z)/R, matching _wave_moment_column.
        forms[2 * eta_index] = -np.tensordot(
            weights * radial**2 / R**2, bilinear[:, 1], axes=(0, 0)
        )
        forms[2 * eta_index + 1] = -np.tensordot(
            weights * radial / R, bilinear[:, 2], axes=(0, 0)
        )
    forms = 0.5 * (forms + np.swapaxes(forms, 1, 2))
    return forms, panels


def _cone_replay(saved, inner, k, breaks, field, wave, nu, hspace, htime):
    """Replay 27 saved cone locations with the current mean plus wave force."""

    points, panels = _grouped_cone_panels(saved, inner, k, breaks, CONE_ORDER)
    mean_jet = jets(field, points, field.tau0, hspace, htime)
    mean_residual = momentum(mean_jet)
    wave_force = _wave_forcing(points, wave, nu, hspace)
    total_residual = mean_residual + wave_force
    diagnostics = saved["cone_problem"]["diagnostics"]
    saved_replay = None
    if MEAN_REPLAY_PATH.exists():
        saved_replay = json.loads(MEAN_REPLAY_PATH.read_text(encoding="utf-8"))
        saved_replay = saved_replay.get("cone_replay", {}).get("rows")
    rows = []
    all_margins = []
    for index, (sl, radii, weights, radius) in enumerate(panels):
        def target(residual):
            return np.array([
                -np.dot(weights * radii**2, residual[sl, 1]) / radius**2,
                -np.dot(weights * radii, residual[sl, 2]) / radius,
            ])

        mean_target = target(mean_residual)
        wave_target = target(wave_force)
        total_target = target(total_residual)
        H = np.asarray(diagnostics[index]["H"], dtype=float)
        mean_margin = H @ mean_target - MARGIN
        wave_margin = H @ wave_target
        total_margin = H @ total_target - MARGIN
        all_margins.extend(total_margin.tolist())
        row = {
            "index": index,
            "kind": diagnostics[index]["kind"],
            "label": diagnostics[index]["label"],
            "physical_center": diagnostics[index]["physical_center"],
            "H": H.tolist(),
            "mean_target": mean_target.tolist(),
            "wave_target": wave_target.tolist(),
            "total_target": total_target.tolist(),
            "mean_margin": mean_margin.tolist(),
            "wave_margin_increment": wave_margin.tolist(),
            "total_margin": total_margin.tolist(),
            "total_pass": bool(np.all(total_margin >= -1.0e-7)),
            "mean_residual_center_norm": float(np.linalg.norm(mean_residual[index])),
            "total_residual_center_norm": float(np.linalg.norm(total_residual[index])),
        }
        if saved_replay is not None and index < len(saved_replay):
            row["saved_corrected_cone_reference"] = saved_replay[index]
        rows.append(row)
    row_pass_count = int(sum(bool(row["total_pass"]) for row in rows))
    inequality_pass_count = int(np.sum(np.asarray(all_margins) >= -1.0e-7))
    return {
        "order": CONE_ORDER,
        "node_count": len(rows),
        "rows": rows,
        "minimum_margin": float(np.min(all_margins)),
        "pass_count": row_pass_count,
        "row_pass_count": row_pass_count,
        "row_count": len(rows),
        "inequality_pass_count": inequality_pass_count,
        "inequality_count": len(all_margins),
        "point_count": int(len(points)),
    }


def run(candidate_path=CANDIDATE_PATH, replay_path=REPLAY_PATH,
        output_path=OUTPUT_PATH):
    started = time.perf_counter()
    saved = json.loads(MEAN_PATH.read_text(encoding="utf-8"))
    candidate_path = Path(candidate_path)
    replay_path = Path(replay_path)
    output_path = Path(output_path)
    candidate = json.loads(candidate_path.read_text(encoding="utf-8"))
    replay = json.loads(replay_path.read_text(encoding="utf-8"))
    mean, mean_report = load_constrained_mean(MEAN_PATH)
    grouped_replacements = int(install_in_field(mean))
    dynamic, dynamic_report = load_dynamic_field(DYNAMIC_PATH)
    k = float(saved["k"])
    tau = float(saved["tau"])
    nu = float(saved["nu"])
    hspace = 5.0e-4 * np.sqrt(nu * tau)
    htime = 1.0e-4 * tau
    tangent = np.asarray(candidate["selected"]["tangent_coefficients"], dtype=float)
    derivatives, pressures = _unpack_full(tangent, 9)
    # The dynamics report stores the selected coefficients but keeps the
    # geometry/carrier in its immutable source wave report.
    if all(key in candidate for key in ("center", "widths", "carrier", "mode", "degree")):
        wave = _wave_metadata(candidate_path)
    else:
        source_wave_path = ROOT / candidate["source"]
        wave = _wave_metadata(source_wave_path)
        wave["source"] = f"{candidate_path.name}:selected + {source_wave_path.name}:geometry"
        wave["potential_coefficients"] = candidate["selected"]["coefficients_original"]
    wave["coefficients"] = np.asarray(
        wave["potential_coefficients"], dtype=float
    )[:, 0] + 1j * np.asarray(wave["potential_coefficients"], dtype=float)[:, 1]
    if len(tangent) != 180 or tangent.shape != (180,):
        raise ValueError("selected tangent vector must have 180 real entries")
    mode0_derivative = np.asarray(derivatives[0], dtype=complex)
    mode0_pressure = np.asarray(pressures[0], dtype=complex)
    if np.max(np.abs(mode0_derivative.imag)) > 1.0e-12 or np.max(np.abs(mode0_pressure.imag)) > 1.0e-12:
        raise ValueError("mode-0 tangent controls should be real")
    breaks = [float(v) for v in saved["quadrature"]["radial_split_breaks"]]
    field = MeanMode0Field(
        mean, wave["center"], wave["widths"], wave["degree"],
        mode0_derivative, mode0_pressure, tau,
    )
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "constraints_maintained": False,
        "grouped_backend_replacements": grouped_replacements,
        "k": k,
        "tau": tau,
        "physical_time": "t=-tau",
        "nu": nu,
        "finite_difference_steps": {"hspace": hspace, "htime": htime},
        "sources": {
            "mean": {"path": MEAN_PATH.name, "sha256": _sha256(MEAN_PATH)},
            "candidate": {"path": candidate_path.name, "sha256": _sha256(candidate_path),
                          "source_sha256_field": candidate.get("source_sha256")},
            "replay": {"path": replay_path.name, "sha256": _sha256(replay_path),
                        "status": replay.get("status")},
            "dynamic": {"path": DYNAMIC_PATH.name, "sha256": _sha256(DYNAMIC_PATH)},
        },
        "candidate_selected_metrics": {
            key: candidate["selected"].get(key)
            for key in ("objective", "corrected_volume_L2", "corrected_max",
                        "growth_lambda", "training_wave_energy", "fitted_wave_energy_rate")
        },
        "wave": {
            key: (value.tolist() if isinstance(value, np.ndarray) else value)
            for key, value in wave.items() if key != "coefficients"
        },
        "tangent_source": {
            "layout": "full_wave_tangent._unpack_full; mode-0 first 36 controls, modes 1/2 angular mean omitted",
            "selected_count": int(len(tangent)),
            "mode0_derivative": mode0_derivative.real.tolist(),
            "mode0_pressure": mode0_pressure.real.tolist(),
            "mode0_derivative_l2": float(np.linalg.norm(mode0_derivative)),
            "mode0_pressure_l2": float(np.linalg.norm(mode0_pressure)),
            "ignored_mode1_derivative_l2": float(np.linalg.norm(derivatives[1])),
            "ignored_mode2_derivative_l2": float(np.linalg.norm(derivatives[2])),
        },
        "quadrature": {
            "moment_order": MOMENT_ORDER,
            "moment_zfactor": MOMENT_ZFACTOR,
            "cone_order": CONE_ORDER,
            "radial_split_breaks": breaks,
            "cone_geometry": "saved 27 locations and H transforms; grouped radial panels",
        },
        "scope": (
            "Instantaneous angular-mean compatibility replay for the locked "
            "co-designed wave. Mode-0 tangent derivative and pressure are added "
            "to the saved mean; complete mode-1 angular-mean force is added. "
            "No optimization, oscillatory full-PDE validation, finite-time "
            "trajectory, continuum cone, or scale-recursion claim."
        ),
    }
    _save(report, output_path)
    print(json.dumps({"stage": "initialized", "grouped_backend_replacements": grouped_replacements,
                      "mode0_derivative_l2": report["tangent_source"]["mode0_derivative_l2"]}), flush=True)

    # Build the reusable moment linearization and wave quadratic forms before
    # the actual replay.  This is cheap relative to the order-96 field replay.
    moment_rows, moment_geometry = _mode0_moment_rows(
        dynamic.inner, tau, breaks, wave["center"], wave["widths"],
        MOMENT_ORDER, MOMENT_ZFACTOR, htime,
    )
    wave_forms, wave_panels = _wave_moment_forms(
        dynamic.inner, tau, k, breaks, wave, nu, hspace, MOMENT_ORDER,
    )
    report["reusable_moment_linearization"] = {
        "moment_rows": moment_rows.tolist(),
        "baseline_moments_source": MEAN_REPLAY_PATH.name,
        "wave_moment_forms": wave_forms.tolist(),
        "mode0_rows_nonzero_columns": [0, 36],
        "row_order": ["eta=-0.2 theta", "eta=-0.2 axial",
                       "eta=+0.2 theta", "eta=+0.2 axial"],
        "geometry": moment_geometry,
        "wave_form_definition": "x.T @ wave_moment_forms[row] @ x with x=[Re(c), Im(c)]",
        "wave_panel_count": len(wave_panels),
    }
    report["status"] = "forms_assembled"
    report["elapsed_forms_seconds"] = float(time.perf_counter() - started)
    _save(report, output_path)
    print(json.dumps({"stage": "forms_assembled", "moment_rows": list(moment_rows.shape),
                      "wave_forms": list(wave_forms.shape),
                      "elapsed_seconds": report["elapsed_forms_seconds"]}), flush=True)

    # Actual order-96 mean replay.  This is the only expensive field call in
    # the experiment; wave forcing is integrated from the same radial panels.
    mean_moments = np.asarray(replay_moments(
        field, dynamic.inner, k, MOMENT_ORDER, MOMENT_ZFACTOR,
        radial_breaks=breaks,
    ), dtype=float)
    wave_moment, _ = _wave_moment_column(
        dynamic.inner, k, breaks, wave, nu, hspace, MOMENT_ORDER,
    )
    joint_moments = mean_moments + wave_moment
    baseline_reference = None
    if MEAN_REPLAY_PATH.exists():
        baseline_reference = json.loads(MEAN_REPLAY_PATH.read_text(encoding="utf-8"))
        baseline_reference = baseline_reference.get("moment_replay")
    report["moment_replay"] = {
        "order": MOMENT_ORDER,
        "mean_mode0_corrected": mean_moments.tolist(),
        "wave_force_column": wave_moment.tolist(),
        "joint_angular_mean": joint_moments.tolist(),
        "mean_mode0_max_abs": float(np.max(np.abs(mean_moments))),
        "joint_max_abs": float(np.max(np.abs(joint_moments))),
        "wave_form_selected_check": (
            np.einsum("i,rij,j->r", np.r_[wave["coefficients"].real,
                                             wave["coefficients"].imag],
                      wave_forms,
                      np.r_[wave["coefficients"].real,
                            wave["coefficients"].imag]).tolist()
        ),
        "reusable_rows_selected_check": (moment_rows @ tangent[:180]).tolist(),
        "saved_mean_reference": baseline_reference,
        "source": "midplane_integrated_moment_balance.evaluate plus wave_mean_flux.single_mode_force radial panels",
    }
    selected_mode0_vector = np.r_[tangent[:36], np.zeros(144)]
    baseline_moments = np.asarray(
        baseline_reference["values"] if baseline_reference is not None
        else mean_moments - moment_rows @ selected_mode0_vector,
        dtype=float,
    )
    linearized_mean = baseline_moments + moment_rows @ selected_mode0_vector
    report["reusable_moment_linearization"]["baseline_moments"] = baseline_moments.tolist()
    report["reusable_moment_linearization"]["baseline_source"] = (
        f"{MEAN_REPLAY_PATH.name}:moment_replay.values"
        if baseline_reference is not None else "mean_mode0_corrected - moment_rows @ selected_mode0_vector"
    )
    report["reusable_moment_linearization"]["selected_mode0_reconstruction"] = {
        "linearized_mean_moments": linearized_mean.tolist(),
        "actual_mean_moments": mean_moments.tolist(),
        "maximum_absolute_error": float(np.max(np.abs(linearized_mean - mean_moments))),
        "note": "mode-0 controls have zero instantaneous velocity; this checks the linear moment response against the actual order-96 replay",
    }
    report["reusable_moment_linearization"]["wave_selected_reconstruction"] = {
        "quadratic_wave_moments": report["moment_replay"]["wave_form_selected_check"],
        "direct_wave_moments": wave_moment.tolist(),
        "maximum_absolute_error": float(np.max(np.abs(
            np.asarray(report["moment_replay"]["wave_form_selected_check"])
            - wave_moment
        ))),
    }
    report["status"] = "moment_complete"
    report["elapsed_moment_seconds"] = float(time.perf_counter() - started)
    _save(report, output_path)
    print(json.dumps({"stage": "moment_complete", "joint_max_abs": report["moment_replay"]["joint_max_abs"],
                      "elapsed_seconds": report["elapsed_moment_seconds"]}), flush=True)

    cone = _cone_replay(
        saved, dynamic.inner, k, breaks, field, wave, nu, hspace, htime,
    )
    report["cone_replay"] = cone
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report, output_path)
    print(json.dumps({"stage": "completed", "joint_moment_max": report["moment_replay"]["joint_max_abs"],
                      "cone_min_margin": cone["minimum_margin"],
                      "cone_pass_count": cone["pass_count"],
                      "cone_row_count": cone["row_count"],
                      "cone_inequality_pass_count": cone["inequality_pass_count"],
                      "cone_point_count": cone["point_count"],
                      "elapsed_seconds": report["elapsed_seconds"]}), flush=True)
    return report


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", default=str(CANDIDATE_PATH))
    parser.add_argument("--replay", default=str(REPLAY_PATH))
    parser.add_argument("--output", default=str(OUTPUT_PATH))
    args = parser.parse_args()
    run(args.candidate, args.replay, args.output)
