"""Bounded full-potential local-time tangent fit for the locked wave candidate.

The initial field is the constrained mean plus the saved mode-1 exact-curl
potential.  At one physical time, fit potential time derivatives and constant
pressure coefficients for modes 0, 1, and 2 against the complete Cartesian
momentum residual.  The fit is a local tangent diagnostic only: mode-0
derivatives change the mean state and no cone or moment constraint is imposed.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from affine_momentum import jets, momentum  # noqa: E402
from broad_meridional_constrained import load_saved_field as load_constrained_mean  # noqa: E402
from broad_wave_mean_fit import patch_nodes, norms  # noqa: E402
from grouped_joined_field import install_in_field  # noqa: E402
from supported_fourier_basis import basis_data  # noqa: E402


MEAN_PATH = ROOT / "broad_meridional_constrained.json"
WAVE_PATH = ROOT / "wave_stress_growth_codesign.json"
OUTPUT_PATH = ROOT / "full_wave_tangent.json"

MODE_LIST = (0, 1, 2)
FIT_ORDER = 6
HOLDOUT_ORDER = 9
FIT_ANGLES = 8
HOLDOUT_ANGLES = 12
FIT_SHIFT = 0.13
HOLDOUT_SHIFT = 0.31
DEGREE = 2


def _save(report):
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _pack_complex(value):
    value = np.asarray(value)
    return np.stack((value.real, value.imag), axis=-1).tolist()


def _decode_complex(value):
    value = np.asarray(value, dtype=float)
    return value[..., 0] + 1j * value[..., 1]


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _grid(field, center, widths, tau, breaks, order, angles, shift):
    meridional, weights, metadata = patch_nodes(
        field, np.asarray(center), np.asarray(widths), tau, order, order, breaks
    )
    theta = shift + np.arange(angles) * 2.0 * np.pi / angles
    points = np.stack(
        (
            meridional[:, 0, None] * np.cos(theta),
            meridional[:, 0, None] * np.sin(theta),
            np.broadcast_to(meridional[:, 2, None], (len(meridional), angles)),
        ),
        axis=-1,
    )
    return points.reshape(-1, 3), np.repeat(weights / angles, angles), metadata


class LocalPotentialField:
    """Mean plus modes 0, 1, and 2 with affine physical-time amplitudes."""

    def __init__(self, mean, center, widths, degree, carriers, initial,
                 derivatives, pressures, tau0):
        self.mean = mean
        self.center = tuple(float(v) for v in center)
        self.widths = tuple(float(v) for v in widths)
        self.degree = int(degree)
        self.carriers = {int(m): np.asarray(v, dtype=float) for m, v in carriers.items()}
        self.initial = tuple(np.asarray(v, dtype=complex) for v in initial)
        self.derivatives = tuple(np.asarray(v, dtype=complex) for v in derivatives)
        self.pressures = tuple(np.asarray(v, dtype=complex) for v in pressures)
        self.tau0 = float(tau0)
        self.nu = mean.nu

    def fields(self, points, tau):
        points = np.asarray(points, dtype=float)
        tau = float(np.asarray(tau).ravel()[0])
        velocity, pressure = self.mean.fields(points, tau)
        velocity = velocity.copy()
        pressure = pressure.copy()
        # t = -tau, so t - t0 = tau0 - tau.
        physical_delta = self.tau0 - tau
        for mode in MODE_LIST:
            carrier = self.carriers[mode]
            coefficient = self.initial[mode] + physical_delta * self.derivatives[mode]
            basis_velocity, basis_pressure, _ = basis_data(
                points, self.center, self.widths, mode, self.degree, carrier
            )
            velocity += np.einsum("niq,q->ni", basis_velocity, coefficient).real
            pressure += np.einsum("nq,q->n", basis_pressure, self.pressures[mode]).real
        return velocity, pressure


def _metric(residual, weights):
    result = norms(residual, weights)
    result["point_count"] = int(len(residual))
    return result


def _basis_columns(points, center, widths, carriers, degree):
    """Full 180-real-column derivative/pressure response matrix."""
    blocks = []
    layout = []
    q = (degree + 1) ** 2
    velocity_q = 3 * q
    for mode in MODE_LIST:
        velocity, _, pressure_gradient = basis_data(
            points, center, widths, mode, degree, carriers[mode]
        )
        if mode == 0:
            for j in range(velocity_q):
                blocks.append(velocity[:, :, j].real.reshape(-1))
                layout.append({"mode": mode, "kind": "velocity", "component": "real", "index": j})
            for j in range(q):
                blocks.append(pressure_gradient[:, :, j].real.reshape(-1))
                layout.append({"mode": mode, "kind": "pressure", "component": "real", "index": j})
        else:
            for kind, tensor, count in (
                ("velocity", velocity, velocity_q),
                ("pressure", pressure_gradient, q),
            ):
                for j in range(count):
                    blocks.append(tensor[:, :, j].real.reshape(-1))
                    layout.append({"mode": mode, "kind": kind, "component": "real", "index": j})
                    blocks.append(-tensor[:, :, j].imag.reshape(-1))
                    layout.append({"mode": mode, "kind": kind, "component": "imag", "index": j})
    return np.stack(blocks, axis=1), layout


def _restricted_columns(points, center, widths, carriers, degree, initial_mode1):
    """Existing 38-real-column tangent: scalar mode-1 amplitude plus p1/p2."""
    q = (degree + 1) ** 2
    velocity, _, _ = basis_data(points, center, widths, 1, degree, carriers[1])
    mode_velocity = np.einsum("niq,q->ni", velocity, initial_mode1)
    blocks = [mode_velocity.real.reshape(-1), -mode_velocity.imag.reshape(-1)]
    layout = [
        {"mode": 1, "kind": "velocity_amplitude", "component": "real"},
        {"mode": 1, "kind": "velocity_amplitude", "component": "imag"},
    ]
    for mode in (1, 2):
        _, _, pressure_gradient = basis_data(
            points, center, widths, mode, degree, carriers[mode]
        )
        for j in range(q):
            blocks.append(pressure_gradient[:, :, j].real.reshape(-1))
            layout.append({"mode": mode, "kind": "pressure", "component": "real", "index": j})
            blocks.append(-pressure_gradient[:, :, j].imag.reshape(-1))
            layout.append({"mode": mode, "kind": "pressure", "component": "imag", "index": j})
    return np.stack(blocks, axis=1), layout


def _weighted_fit(design, residual, weights):
    sqrt_weights = np.sqrt(np.asarray(weights, dtype=float))
    weighted_design = design * np.repeat(sqrt_weights, 3)[:, None]
    weighted_rhs = (-residual * sqrt_weights[:, None]).reshape(-1)
    scaled_columns = np.maximum(np.linalg.norm(weighted_design, axis=0), 1.0e-30)
    normalized = weighted_design / scaled_columns[None, :]
    coefficients_scaled, _, rank, singular = np.linalg.lstsq(
        normalized, weighted_rhs, rcond=1.0e-10
    )
    coefficients = coefficients_scaled / scaled_columns
    return coefficients, {
        "rank": int(rank),
        "column_count": int(design.shape[1]),
        "singular_values": singular.tolist(),
        "column_scales": scaled_columns.tolist(),
    }


def _unpack_full(x, q):
    x = np.asarray(x, dtype=float)
    velocity_q = 3 * q
    cursor = 0
    derivatives = []
    pressures = []
    for mode in MODE_LIST:
        if mode == 0:
            derivatives.append(x[cursor:cursor + velocity_q].astype(complex))
            cursor += velocity_q
            pressures.append(x[cursor:cursor + q].astype(complex))
            cursor += q
        else:
            derivatives.append(
                x[cursor:cursor + 2 * velocity_q:2]
                + 1j * x[cursor + 1:cursor + 2 * velocity_q:2]
            )
            cursor += 2 * velocity_q
            pressures.append(
                x[cursor:cursor + 2 * q:2]
                + 1j * x[cursor + 1:cursor + 2 * q:2]
            )
            cursor += 2 * q
    if cursor != len(x):
        raise ValueError(f"Full coefficient unpack consumed {cursor} of {len(x)}")
    return tuple(derivatives), tuple(pressures)


def _unpack_restricted(x, q, initial_mode1):
    x = np.asarray(x, dtype=float)
    velocity_q = len(initial_mode1)
    amplitude = x[0] + 1j * x[1]
    derivatives = (
        np.zeros(velocity_q, dtype=complex),
        amplitude * initial_mode1,
        np.zeros(velocity_q, dtype=complex),
    )
    cursor = 2
    pressures = [np.zeros(q, dtype=complex)]
    for _ in (1, 2):
        pressures.append(
            x[cursor:cursor + 2 * q:2]
            + 1j * x[cursor + 1:cursor + 2 * q:2]
        )
        cursor += 2 * q
    if cursor != len(x):
        raise ValueError(f"Restricted coefficient unpack consumed {cursor} of {len(x)}")
    return tuple(derivatives), tuple(pressures)


def _mode_energy_data(points, weights, center, widths, carrier, coefficient,
                      derivative):
    velocity = basis_data(points, center, widths, 1, DEGREE, carrier)[0]
    weighted = np.asarray(weights)[:, None]
    matrix = velocity.reshape(-1, velocity.shape[-1]).conj().T @ (
        np.repeat(weighted, 3, axis=0) * velocity.reshape(-1, velocity.shape[-1])
    )
    energy = 0.25 * float(np.real(np.vdot(coefficient, matrix @ coefficient)))
    rate = 0.5 * float(np.real(np.vdot(coefficient, matrix @ derivative)))
    return {
        "physical_mode1_energy": energy,
        "physical_mode1_energy_rate": rate,
        "complex_amplitude_norm": float(np.linalg.norm(coefficient)),
        "complex_derivative_norm": float(np.linalg.norm(derivative)),
    }


def _source_snapshot(mean_report, wave_report, center, widths, carriers, breaks):
    selected = wave_report["selected"]
    return {
        "mean": {
            "source_report": MEAN_PATH.name,
            "source_sha256": _sha256(MEAN_PATH),
            "k": float(mean_report["k"]),
            "tau": float(mean_report["tau"]),
            "nu": float(mean_report["nu"]),
            "coefficients": mean_report["coefficients"],
            "mode_names": mean_report.get("mode_names", []),
            "radial_split_breaks": list(breaks),
        },
        "wave": {
            "source_report": WAVE_PATH.name,
            "source_sha256": _sha256(WAVE_PATH),
            "selected_trial_index": int(selected.get("trial_index", -1)),
            "center": list(center),
            "widths": list(widths),
            "degree": int(wave_report["degree"]),
            "mode": int(wave_report["mode"]),
            "carrier": np.asarray(wave_report["carrier"], dtype=float).tolist(),
            "coefficients_original": selected["coefficients_original"],
        },
        "geometry": {
            "radial_breaks": list(breaks),
            "fit_order": FIT_ORDER,
            "holdout_order": HOLDOUT_ORDER,
            "fit_angles": FIT_ANGLES,
            "holdout_angles": HOLDOUT_ANGLES,
        },
        "convention": "physical time t=-tau; u=mean+Re(sum V_m c_m), pressure=mean+Re(sum P_m p_m)",
    }


def run():
    started = time.perf_counter()
    mean_report = json.loads(MEAN_PATH.read_text(encoding="utf-8"))
    wave_report = json.loads(WAVE_PATH.read_text(encoding="utf-8"))
    if not wave_report.get("selected"):
        raise ValueError("Locked wave co-design report has no selected candidate")
    mean, _ = load_constrained_mean()
    replacements = install_in_field(mean)
    center = tuple(float(v) for v in wave_report["center"])
    widths = tuple(float(v) for v in wave_report["widths"])
    degree = int(wave_report["degree"])
    if degree != DEGREE or int(wave_report["mode"]) != 1:
        raise ValueError("This bounded diagnostic expects the locked degree-2 mode-1 candidate")
    carrier1 = np.asarray(wave_report["carrier"], dtype=float)
    carriers = {0: np.zeros(2), 1: carrier1, 2: 2.0 * carrier1}
    initial_mode1 = _decode_complex(wave_report["selected"]["coefficients_original"])
    q = (degree + 1) ** 2
    if initial_mode1.shape != (3 * q,):
        raise ValueError(f"Expected 27 packed wave coefficients, got {initial_mode1.shape}")
    breaks = sorted(set(float(v) for v in mean_report["quadrature"]["radial_split_breaks"]) | {0.62})
    tau0 = float(mean_report["tau"])
    hspace = 5.0e-4 * np.sqrt(mean.nu * tau0)
    htime = 1.0e-4 * tau0
    initial = (np.zeros(3 * q, complex), initial_mode1, np.zeros(3 * q, complex))
    zero_derivatives = tuple(np.zeros(3 * q, complex) for _ in MODE_LIST)
    zero_pressures = tuple(np.zeros(q, complex) for _ in MODE_LIST)
    frozen = LocalPotentialField(
        mean, center, widths, degree, carriers, initial,
        zero_derivatives, zero_pressures, tau0,
    )
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "wave_integrated": False,
        "constraints_maintained": False,
        "grouped_backend_replacements": int(replacements),
        "scope": (
            "Local full-potential tangent fit for the locked mode-1 candidate. "
            "Modes 0, 1, and 2 potential time derivatives and constant pressures "
            "are fitted on a split physical patch. Mode-0 derivatives change the "
            "mean state; no moments, cones, finite-time integration, PDE, or "
            "scale-recursion acceptance is claimed."
        ),
        "inputs": _source_snapshot(mean_report, wave_report, center, widths, carriers, breaks),
        "initial_k": float(mean_report["k"]),
        "initial_tau": tau0,
        "viscosity": float(mean.nu),
        "degree": degree,
        "mode_list": list(MODE_LIST),
        "basis_q": q,
        "timesteps": {"hspace": hspace, "htime": htime},
        "fit_geometry": {
            "order": FIT_ORDER,
            "angles": FIT_ANGLES,
            "angle_shift": FIT_SHIFT,
            "radial_breaks": breaks,
        },
        "holdout_geometry": {
            "order": HOLDOUT_ORDER,
            "angles": HOLDOUT_ANGLES,
            "angle_shift": HOLDOUT_SHIFT,
            "radial_breaks": breaks,
        },
    }
    previous = None
    if OUTPUT_PATH.exists():
        try:
            candidate = json.loads(OUTPUT_PATH.read_text(encoding="utf-8"))
            candidate_inputs = candidate.get("inputs", {})
            if (
                candidate.get("frozen_cache")
                and candidate_inputs.get("mean", {}).get("source_sha256") == _sha256(MEAN_PATH)
                and candidate_inputs.get("wave", {}).get("source_sha256") == _sha256(WAVE_PATH)
            ):
                previous = candidate
        except (OSError, ValueError, KeyError):
            previous = None
    _save(report)
    print(json.dumps({"stage": "initializing", "grouped_backend_replacements": replacements,
                      "reuse_frozen_cache": previous is not None}), flush=True)

    if previous is not None:
        cache = previous["frozen_cache"]
        fit_points = np.asarray(cache["fit_points"], dtype=float)
        fit_weights = np.asarray(cache["fit_weights"], dtype=float)
        frozen_fit_residual = np.asarray(cache["fit_residual"], dtype=float)
        hold_points = np.asarray(cache["holdout_points"], dtype=float)
        hold_weights = np.asarray(cache["holdout_weights"], dtype=float)
        frozen_hold_residual = np.asarray(cache["holdout_residual"], dtype=float)
        fit_meta = previous["geometry"]["fit"]
        hold_meta = previous["geometry"]["holdout"]
        report["geometry"] = {"fit": fit_meta, "holdout": hold_meta}
        report["baseline"] = previous["baseline"]
        report["frozen_cache"] = cache
        report["status"] = "baseline_reused"
        _save(report)
        print(json.dumps({"stage": "baseline_reused", "fit_points": len(fit_points),
                          "holdout_points": len(hold_points)}), flush=True)
    else:
        fit_points, fit_weights, fit_meta = _grid(
            mean, center, widths, tau0, breaks, FIT_ORDER, FIT_ANGLES, FIT_SHIFT
        )
        hold_points, hold_weights, hold_meta = _grid(
            mean, center, widths, tau0, breaks, HOLDOUT_ORDER, HOLDOUT_ANGLES, HOLDOUT_SHIFT
        )
        print(json.dumps({
            "stage": "baseline_fit_started",
            "fit_points": len(fit_points),
            "holdout_points": len(hold_points),
        }), flush=True)
        mean_fit_residual = momentum(jets(mean, fit_points, tau0, hspace, htime))
        frozen_fit_residual = momentum(jets(frozen, fit_points, tau0, hspace, htime))
        mean_hold_residual = momentum(jets(mean, hold_points, tau0, hspace, htime))
        frozen_hold_residual = momentum(jets(frozen, hold_points, tau0, hspace, htime))
        report["geometry"] = {"fit": fit_meta, "holdout": hold_meta}
        report["baseline"] = {
            "mean_only_fit": _metric(mean_fit_residual, fit_weights),
            "frozen_wave_fit": _metric(frozen_fit_residual, fit_weights),
            "mean_only_holdout": _metric(mean_hold_residual, hold_weights),
            "frozen_wave_holdout": _metric(frozen_hold_residual, hold_weights),
        }
        # Keep the expensive Cartesian FD cache in the JSON checkpoint so later
        # basis or constraint experiments can reuse the identical residual/grid.
        report["frozen_cache"] = {
            "fit_points": fit_points.tolist(),
            "fit_weights": fit_weights.tolist(),
            "fit_residual": frozen_fit_residual.tolist(),
            "holdout_points": hold_points.tolist(),
            "holdout_weights": hold_weights.tolist(),
            "holdout_residual": frozen_hold_residual.tolist(),
            "weighting": "physical cylindrical volume weights repeated over angles",
        }
        report["status"] = "baseline_assembled"
        _save(report)
        print(json.dumps({
            "stage": "baseline_assembled",
            "mean_fit_rms": report["baseline"]["mean_only_fit"]["volume_RMS"],
            "frozen_fit_rms": report["baseline"]["frozen_wave_fit"]["volume_RMS"],
        }), flush=True)

    full_design, full_layout = _basis_columns(
        fit_points, center, widths, carriers, degree
    )
    restricted_design, restricted_layout = _restricted_columns(
        fit_points, center, widths, carriers, degree, initial_mode1
    )
    full_x, full_svd = _weighted_fit(full_design, frozen_fit_residual, fit_weights)
    restricted_x, restricted_svd = _weighted_fit(
        restricted_design, frozen_fit_residual, fit_weights
    )
    report["design"] = {
        "full_layout": full_layout,
        "full_svd": full_svd,
        "restricted_layout": restricted_layout,
        "restricted_svd": restricted_svd,
        "definition": "weighted residual columns; full columns are mode0/1/2 velocity derivatives and pressure gradients; restricted columns are scalar mode1 amplitude derivative plus mode1/2 pressures",
    }
    report["status"] = "design_assembled"
    _save(report)
    print(json.dumps({
        "stage": "design_assembled",
        "full_rank": full_svd["rank"],
        "full_columns": full_svd["column_count"],
        "restricted_rank": restricted_svd["rank"],
        "restricted_columns": restricted_svd["column_count"],
    }), flush=True)

    full_derivatives, full_pressures = _unpack_full(full_x, q)
    restricted_derivatives, restricted_pressures = _unpack_restricted(
        restricted_x, q, initial_mode1
    )
    frozen_fit_flat = frozen_fit_residual.reshape(-1)
    full_fit_linear = frozen_fit_flat + full_design @ full_x
    restricted_fit_linear = frozen_fit_flat + restricted_design @ restricted_x
    full_field = LocalPotentialField(
        mean, center, widths, degree, carriers, initial,
        full_derivatives, full_pressures, tau0,
    )
    restricted_field = LocalPotentialField(
        mean, center, widths, degree, carriers, initial,
        restricted_derivatives, restricted_pressures, tau0,
    )
    report["fits"] = {
        "restricted38": {
            "fit_linear": _metric(restricted_fit_linear.reshape(-1, 3), fit_weights),
            "holdout_actual_fd": _metric(
                momentum(jets(restricted_field, hold_points, tau0, hspace, htime)), hold_weights
            ),
            "coefficients": restricted_x.tolist(),
            "derivative_norms": [float(np.linalg.norm(v)) for v in restricted_derivatives],
            "pressure_norms": [float(np.linalg.norm(v)) for v in restricted_pressures],
            "mode1_energy": _mode_energy_data(
                fit_points, fit_weights, center, widths, carriers[1],
                initial_mode1, restricted_derivatives[1],
            ),
        },
        "full180": {
            "fit_linear": _metric(full_fit_linear.reshape(-1, 3), fit_weights),
            "holdout_actual_fd": _metric(
                momentum(jets(full_field, hold_points, tau0, hspace, htime)), hold_weights
            ),
            "coefficients": full_x.tolist(),
            "derivative_norms": [float(np.linalg.norm(v)) for v in full_derivatives],
            "pressure_norms": [float(np.linalg.norm(v)) for v in full_pressures],
            "mode1_energy": _mode_energy_data(
                fit_points, fit_weights, center, widths, carriers[1],
                initial_mode1, full_derivatives[1],
            ),
        },
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report)
    print(json.dumps({
        "stage": "completed",
        "status": report["status"],
        "full_fit_rms": report["fits"]["full180"]["fit_linear"]["volume_RMS"],
        "full_holdout_rms": report["fits"]["full180"]["holdout_actual_fd"]["volume_RMS"],
        "restricted_holdout_rms": report["fits"]["restricted38"]["holdout_actual_fd"]["volume_RMS"],
        "full_mode1_energy_rate": report["fits"]["full180"]["mode1_energy"]["physical_mode1_energy_rate"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    run()
