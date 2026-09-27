"""Affine instantaneous shape-trend rows for the frozen wave tangent.

The rows in this module describe sampled fixed-cylinder enstrophy-shape rates
with respect to the 180 real tangent controls used by ``full_wave_tangent``.
They are a local diagnostic for a later constrained fit.  The velocity and
vorticity of the frozen candidate are evaluated with the same finite spatial
stencil as ``wave_tangent_observables``; the tangent curls use the analytic
supported-Fourier jets.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from broad_meridional_constrained import load_saved_field
from broad_shear_dynamic_control import load_saved_field as load_dynamic
from full_wave_tangent import LocalPotentialField, _basis_columns, _unpack_full
from grouped_joined_field import install_in_field
from meridional_state_cache import value_and_pressure_modes
from supported_fourier_analytic_jets import basis_jets
from vortex_state_observables import fixed_cylinder


ROOT = Path(__file__).resolve().parent
DEFAULT_CANDIDATE = ROOT / "wave_moment_cone_tangent.json"
DEFAULT_OBSERVABLES = ROOT / "wave_tangent_observables.json"
DEFAULT_OUTPUT = ROOT / "wave_shape_tangent_rows.json"
FROZEN_CACHE = ROOT / "full_wave_frozen_cache.json"

MODE_LIST = (0, 1, 2)
DEGREE = 2
NU = 0.01


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(report: dict, path: Path) -> None:
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _curl(gradient: np.ndarray) -> np.ndarray:
    """Curl columns from J[..., component, derivative, column]."""
    return np.stack(
        (
            gradient[:, 2, 1, :] - gradient[:, 1, 2, :],
            gradient[:, 0, 2, :] - gradient[:, 2, 0, :],
            gradient[:, 1, 0, :] - gradient[:, 0, 1, :],
        ),
        axis=1,
    )


def _state_with_fd_jets(field, points: np.ndarray, weights: np.ndarray, tau: float):
    """Return the frozen state using vortex_state_observables' FD stencil."""
    h = 5.0e-4 * np.sqrt(field.nu * tau)
    velocity = field.fields(points, tau)[0]
    gradient = np.empty((len(points), 3, 3), dtype=float)
    for axis, direction in enumerate(np.eye(3)):
        um2, um, up, up2 = (
            field.fields(points + multiple * h * direction, tau)[0]
            for multiple in (-2, -1, 1, 2)
        )
        gradient[:, :, axis] = (um2 - 8.0 * um + 8.0 * up - up2) / (12.0 * h)
    omega = np.column_stack(
        (
            gradient[:, 2, 1] - gradient[:, 1, 2],
            gradient[:, 0, 2] - gradient[:, 2, 0],
            gradient[:, 1, 0] - gradient[:, 0, 1],
        )
    )
    density = np.sum(omega * omega, axis=1)
    energy = float(weights @ density)
    if not np.isfinite(energy) or energy <= 0.0:
        raise ValueError("Nonpositive or nonfinite sampled enstrophy")
    radius = np.linalg.norm(points[:, :2], axis=1)
    theta_velocity = (-points[:, 1] * velocity[:, 0] + points[:, 0] * velocity[:, 1]) / radius
    angular_speed = theta_velocity / radius
    zmean = float(weights @ (density * points[:, 2]) / energy)
    radial_squared = float(weights @ (density * radius * radius) / energy)
    axial_offset = points[:, 2] - zmean
    axial_squared = float(weights @ (density * axial_offset * axial_offset) / energy)
    radial_rms = float(np.sqrt(radial_squared))
    axial_rms = float(np.sqrt(axial_squared))
    spin = float(weights @ (density * angular_speed) / energy)
    reference = {
        "cylinder_enstrophy": energy,
        "enstrophy_radial_rms": radial_rms,
        "enstrophy_axial_rms": axial_rms,
        "enstrophy_aspect_ratio": axial_rms / radial_rms,
        "enstrophy_weighted_angular_speed": spin,
        "zmean": zmean,
        "radius_squared": radial_squared,
        "axial_squared": axial_squared,
    }
    return velocity, omega, density, angular_speed, reference


def _metric_rate(
    velocity_t: np.ndarray,
    omega_t: np.ndarray,
    points: np.ndarray,
    weights: np.ndarray,
    reference: dict,
    dt_dk: float,
) -> np.ndarray:
    """Return d/dk of (radial RMS, aspect, signed angular speed)."""
    radius = np.linalg.norm(points[:, :2], axis=1)
    density = reference["density"]
    omega = reference["omega"]
    energy = reference["cylinder_enstrophy"]
    radial_rms = reference["enstrophy_radial_rms"]
    axial_rms = reference["enstrophy_axial_rms"]
    radial_squared = reference["radius_squared"]
    axial_squared = reference["axial_squared"]
    spin = reference["enstrophy_weighted_angular_speed"]
    zmean = reference["zmean"]

    density_t = 2.0 * np.sum(omega * omega_t, axis=1)
    d_radial = float(
        weights @ ((radius * radius - radial_squared) * density_t)
        / (2.0 * radial_rms * energy)
    )
    axial_offset = points[:, 2] - zmean
    d_axial = float(
        weights @ ((axial_offset * axial_offset - axial_squared) * density_t)
        / (2.0 * axial_rms * energy)
    )
    d_aspect = d_axial / radial_rms - axial_rms * d_radial / radial_rms**2
    theta_velocity_t = (
        -points[:, 1] * velocity_t[:, 0] + points[:, 0] * velocity_t[:, 1]
    ) / radius
    angular_speed_t = theta_velocity_t / radius
    d_spin = float(
        weights
        @ (density * angular_speed_t + (reference["angular_speed"] - spin) * density_t)
        / energy
    )
    return dt_dk * np.array((d_radial, d_aspect, d_spin), dtype=float)


def _build_field(candidate: dict, snapshot: dict):
    mean, mean_report = load_saved_field()
    install_in_field(mean)
    dynamic, dynamic_report = load_dynamic()
    install_in_field(dynamic)
    base, *_ = value_and_pressure_modes(dynamic, dynamic_report)
    k0 = float(dynamic_report["k"])
    tau0 = float(snapshot["inputs"]["mean"]["tau"])
    if not np.isclose(0.5 * 2.0 ** (-k0), tau0, rtol=1.0e-12, atol=0.0):
        raise ValueError("Reference times differ")

    geometry = snapshot["inputs"]["wave"]
    selected = candidate["selected"]
    packed = np.asarray(selected["coefficients_original"], dtype=float)
    wave = packed[:, 0] + 1j * packed[:, 1]
    tangent = np.asarray(selected["tangent_coefficients"], dtype=float)
    derivatives, pressures = _unpack_full(tangent, (DEGREE + 1) ** 2)
    # Pressure does not affect shape observables, so omit it exactly.
    pressures = tuple(np.zeros_like(value) for value in pressures)
    carrier = np.asarray(geometry["carrier"], dtype=float)
    field = LocalPotentialField(
        mean,
        geometry["center"],
        geometry["widths"],
        geometry["degree"],
        {0: np.zeros(2), 1: carrier, 2: 2.0 * carrier},
        (np.zeros(3 * (DEGREE + 1) ** 2, complex), wave,
         np.zeros(3 * (DEGREE + 1) ** 2, complex)),
        derivatives,
        pressures,
        tau0,
    )
    return field, base, dynamic_report, geometry, tangent, mean_report


def run(candidate_path=DEFAULT_CANDIDATE, observables_path=DEFAULT_OBSERVABLES,
        output_path=DEFAULT_OUTPUT):
    started = time.perf_counter()
    candidate_path = Path(candidate_path)
    observables_path = Path(observables_path)
    output_path = Path(output_path)
    candidate_raw = candidate_path.read_bytes()
    candidate = json.loads(candidate_raw)
    if candidate.get("status") != "completed":
        raise ValueError("Candidate must be frozen before shape-row construction")
    observables_raw = observables_path.read_bytes()
    observables = json.loads(observables_raw)
    if observables.get("status") != "completed":
        raise ValueError("Shape observables must be frozen before row construction")
    snapshot = json.loads(FROZEN_CACHE.read_text(encoding="utf-8"))
    field, base, dynamic_report, geometry, tangent, mean_report = _build_field(
        candidate, snapshot
    )
    k0 = float(dynamic_report["k"])
    tau0 = float(snapshot["inputs"]["mean"]["tau"])
    points, weights, domain = fixed_cylinder(base.inner, k0, angles=12)

    report = {
        "status": "running",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "source": candidate_path.name,
        "source_sha256": hashlib.sha256(candidate_raw).hexdigest(),
        "observables_source": observables_path.name,
        "observables_sha256": hashlib.sha256(observables_raw).hexdigest(),
        "reference_k": k0,
        "physical_time": -tau0,
        "dt_dk": float(np.log(2.0) * tau0),
        "domain": domain,
        "point_count": int(len(points)),
        "control_count": int(len(tangent)),
        "scope": (
            "Instantaneous fixed-cylinder affine shape rows only. Baseline is "
            "inferred from the frozen central finite-difference observable slope "
            "minus the selected tangent response; it is not an independent "
            "analytic dynamics model. No trajectory, PDE, or scale-recursion claim."
        ),
    }
    _save(report, output_path)

    velocity, omega, density, angular_speed, reference = _state_with_fd_jets(
        field, points, weights, tau0
    )
    reference["density"] = density
    reference["omega"] = omega
    reference["angular_speed"] = angular_speed
    reference_detail = {
        key: float(value)
        for key, value in reference.items()
        if key not in {"density", "omega", "angular_speed"}
    }
    # Keep the compact ordered vector explicit for downstream constrained-fit
    # code, while retaining all sampled state quantities for auditability.
    report["reference_observables"] = [
        reference["enstrophy_radial_rms"],
        reference["enstrophy_aspect_ratio"],
        reference["enstrophy_weighted_angular_speed"],
    ]
    report["reference_observable_names"] = [
        "enstrophy_radial_rms",
        "enstrophy_aspect_ratio",
        "enstrophy_weighted_angular_speed",
    ]
    report["reference_observables_detail"] = reference_detail
    report["reference_observable_errors_vs_frozen_center"] = {
        key: float(
            reference[key]
            - observables["rows"][1][key]
        )
        for key in (
            "cylinder_enstrophy",
            "enstrophy_radial_rms",
            "enstrophy_aspect_ratio",
            "enstrophy_weighted_angular_speed",
        )
    }
    _save(report, output_path)

    carriers = {
        0: np.zeros(2),
        1: np.asarray(geometry["carrier"], dtype=float),
        2: 2.0 * np.asarray(geometry["carrier"], dtype=float),
    }
    row_blocks = []
    analytic_scales = []
    for mode in MODE_LIST:
        values, gradients, _, _, _ = basis_jets(
            points,
            geometry["center"],
            geometry["widths"],
            mode,
            DEGREE,
            carriers[mode],
            NU,
        )
        curls = _curl(gradients)
        if mode == 0:
            for index in range(values.shape[-1]):
                row_blocks.append(
                    _metric_rate(
                        values[:, :, index].real,
                        curls[:, :, index].real,
                        points,
                        weights,
                        reference,
                        report["dt_dk"],
                    )
                )
                analytic_scales.append({"mode": mode, "kind": "velocity", "component": "real", "index": index})
            for index in range((DEGREE + 1) ** 2):
                row_blocks.append(np.zeros(3, dtype=float))
                analytic_scales.append({"mode": mode, "kind": "pressure", "component": "real", "index": index})
        else:
            for kind, tensor in (("velocity", values), ("pressure", None)):
                count = values.shape[-1] if kind == "velocity" else (DEGREE + 1) ** 2
                if kind == "pressure":
                    for index in range(count):
                        row_blocks.extend((np.zeros(3), np.zeros(3)))
                        analytic_scales.extend(
                            (
                                {"mode": mode, "kind": kind, "component": "real", "index": index},
                                {"mode": mode, "kind": kind, "component": "imag", "index": index},
                            )
                        )
                    continue
                for index in range(count):
                    row_blocks.append(
                        _metric_rate(
                            tensor[:, :, index].real,
                            curls[:, :, index].real,
                            points,
                            weights,
                            reference,
                            report["dt_dk"],
                        )
                    )
                    row_blocks.append(
                        _metric_rate(
                            -tensor[:, :, index].imag,
                            -curls[:, :, index].imag,
                            points,
                            weights,
                            reference,
                            report["dt_dk"],
                        )
                    )
                    analytic_scales.extend(
                        (
                            {"mode": mode, "kind": kind, "component": "real", "index": index},
                            {"mode": mode, "kind": kind, "component": "imag", "index": index},
                        )
                    )
    shape_rows = np.asarray(row_blocks, dtype=float).T
    if shape_rows.shape != (3, len(tangent)):
        raise ValueError(f"Unexpected shape rows {shape_rows.shape}, tangent {len(tangent)}")

    # Ask full_wave_tangent for its authoritative 180-column ordering.  The
    # matrix itself is intentionally discarded; pressure columns are zero here.
    _, layout = _basis_columns(
        points, geometry["center"], geometry["widths"], carriers, DEGREE
    )
    if len(layout) != len(tangent):
        raise ValueError("Full tangent layout length differs from candidate")
    if [
        (a.get("mode"), a.get("kind"), a.get("component"), a.get("index"))
        for a in layout
    ] != [
        (a.get("mode"), a.get("kind"), a.get("component"), a.get("index"))
        for a in analytic_scales
    ]:
        raise ValueError("Analytic row order differs from full_wave_tangent layout")

    finite_slopes = np.asarray(
        [
            observables["derivative_per_k"]["enstrophy_radial_rms"],
            observables["derivative_per_k"]["enstrophy_aspect_ratio"],
            observables["derivative_per_k"]["enstrophy_weighted_angular_speed"],
        ],
        dtype=float,
    )
    selected_response = shape_rows @ tangent
    baseline = finite_slopes - selected_response
    predicted_selected = baseline + selected_response
    reference_spin = reference["enstrophy_weighted_angular_speed"]
    desired_signs = np.array((-1.0, 1.0, np.sign(reference_spin)), dtype=float)
    report.update(
        {
            "shape_metrics": [
                "enstrophy_radial_rms",
                "enstrophy_aspect_ratio",
                "enstrophy_weighted_angular_speed",
            ],
            "shape_control_rows": shape_rows.tolist(),
            "shape_baseline": baseline.tolist(),
            "finite_difference_derivative_per_k": finite_slopes.tolist(),
            "selected_tangent_response_per_k": selected_response.tolist(),
            "selected_reconstructed_derivative_per_k": predicted_selected.tolist(),
            "selected_reconstruction_max_abs_error": float(
                np.max(np.abs(predicted_selected - finite_slopes))
            ),
            "desired_signs": desired_signs.tolist(),
            "inequality_convention": (
                "shape_baseline + shape_control_rows @ x has raw components "
                "[radial_rms_rate, aspect_rate, signed_weighted_angular_speed_rate]. "
                "The desired signed rates are desired_signs multiplied by those "
                "components and constrained nonnegative."
            ),
            "control_layout": layout,
            "pressure_column_indices": [
                index for index, item in enumerate(layout) if item["kind"] == "pressure"
            ],
            "selected_coefficients": tangent.tolist(),
            "row_norms": np.linalg.norm(shape_rows, axis=1).tolist(),
            "row_max_abs": np.max(np.abs(shape_rows), axis=1).tolist(),
            "sampled_direction_flags": observables["sampled_direction_flags"],
            "source_snapshots": {
                "mean_report": mean_report.get("source", "broad_meridional_constrained.json"),
                "candidate_source_sha256": hashlib.sha256(candidate_raw).hexdigest(),
                "observable_source_sha256": hashlib.sha256(observables_raw).hexdigest(),
            },
        }
    )
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(report, output_path)
    print(
        json.dumps(
            {
                "status": report["status"],
                "rows": list(shape_rows.shape),
                "selected_reconstruction_max_abs_error": report[
                    "selected_reconstruction_max_abs_error"
                ],
                "finite_difference_derivative_per_k": finite_slopes.tolist(),
                "shape_baseline": baseline.tolist(),
                "elapsed_seconds": report["elapsed_seconds"],
            }
        ),
        flush=True,
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--observables", type=Path, default=DEFAULT_OBSERVABLES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.candidate, args.observables, args.output)
