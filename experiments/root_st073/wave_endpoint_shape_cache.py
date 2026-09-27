"""Finite-endpoint geometry oracle for the frozen wave tangent.

The cache separates the expensive mean endpoint field (velocity and its
five-point spatial Jacobian) from the analytic compact Fourier tangent.  The
resulting endpoint state is affine in the 180 real tangent controls used by
``full_wave_tangent``.  ``EndpointShape.evaluate`` computes the three sampled
enstrophy shape observables and their exact Jacobian with respect to those
controls.

This is a fixed-cylinder endpoint diagnostic.  It does not advance an NS
solution or establish PDE validity or scale recursion.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path

for _name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np

from broad_meridional_constrained import load_saved_field
from broad_shear_dynamic_control import load_saved_field as load_dynamic
from full_wave_tangent import LocalPotentialField, _unpack_full
from grouped_joined_field import install_in_field
from meridional_state_cache import value_and_pressure_modes
from supported_fourier_analytic_jets import basis_jets
from vortex_state_observables import fixed_cylinder


ROOT = Path(__file__).resolve().parent
DEFAULT_CACHE = ROOT / "wave_endpoint_shape_cache.npz"
DEFAULT_REPORT = ROOT / "wave_endpoint_shape_cache.json"
DEFAULT_CANDIDATE = ROOT / "shape_direction_margin_tangent.json"
DEFAULT_SAVED_ORIGINAL = ROOT / "wave_tangent_observables.json"
DEFAULT_SAVED_SHAPE = ROOT / "shape_direction_margin_observables.json"
FROZEN_GEOMETRY = ROOT / "full_wave_frozen_cache.json"
SHAPE_ROWS = ROOT / "wave_shape_tangent_rows.json"

MODE_LIST = (0, 1, 2)
DEGREE = 2
NU = 0.01
DELTA_K = 1.0e-6


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _decode(value):
    value = np.asarray(value, dtype=float)
    return value[..., 0] + 1j * value[..., 1]


def _pack(value):
    value = np.asarray(value)
    return np.stack((value.real, value.imag), axis=-1)


def _curl(jacobian):
    """Curl columns from J[..., component, derivative, column]."""
    return np.stack(
        (
            jacobian[:, 2, 1, :] - jacobian[:, 1, 2, :],
            jacobian[:, 0, 2, :] - jacobian[:, 2, 0, :],
            jacobian[:, 1, 0, :] - jacobian[:, 0, 1, :],
        ),
        axis=1,
    )


def _mean_endpoint_fd(mean, points, tau):
    """Mean-only endpoint u,J using the same centered five-point stencil."""
    h = 5.0e-4 * np.sqrt(float(mean.nu) * float(tau))
    velocity = np.asarray(mean.fields(points, tau)[0], dtype=float)
    gradient = np.empty((len(points), 3, 3), dtype=float)
    for axis, direction in enumerate(np.eye(3)):
        samples = [
            np.asarray(mean.fields(points + multiple * h * direction, tau)[0], dtype=float)
            for multiple in (-2, -1, 1, 2)
        ]
        um2, um, up, up2 = samples
        gradient[:, :, axis] = (um2 - 8.0 * um + 8.0 * up - up2) / (12.0 * h)
    return velocity, gradient, h


def _basis_response(points, geometry, dtau):
    """Analytic V,J response for all 180 real tangent controls."""
    center = geometry["center"]
    widths = geometry["widths"]
    degree = int(geometry["degree"])
    carrier = np.asarray(geometry["carrier"], dtype=float)
    carriers = {0: np.zeros(2), 1: carrier, 2: 2.0 * carrier}
    q = (degree + 1) ** 2
    blocks_v = []
    blocks_j = []
    layout = []
    precision = []
    for mode in MODE_LIST:
        values, gradients, _, _, _ = basis_jets(
            points, center, widths, mode, degree, carriers[mode], NU
        )
        if mode == 0:
            for index in range(3 * q):
                blocks_v.append(dtau * values[:, :, index].real)
                blocks_j.append(dtau * gradients[:, :, :, index].real)
                layout.append((mode, "velocity", "real", index))
                precision.append("analytic_basis_jets")
            for index in range(q):
                blocks_v.append(np.zeros((len(points), 3), dtype=float))
                blocks_j.append(np.zeros((len(points), 3, 3), dtype=float))
                layout.append((mode, "pressure", "real", index))
                precision.append("zero_velocity_pressure")
        else:
            for index in range(3 * q):
                blocks_v.append(dtau * values[:, :, index].real)
                blocks_j.append(dtau * gradients[:, :, :, index].real)
                layout.append((mode, "velocity", "real", index))
                precision.append("analytic_basis_jets")
                blocks_v.append(-dtau * values[:, :, index].imag)
                blocks_j.append(-dtau * gradients[:, :, :, index].imag)
                layout.append((mode, "velocity", "imag", index))
                precision.append("analytic_basis_jets")
            for index in range(q):
                blocks_v.append(np.zeros((len(points), 3), dtype=float))
                blocks_j.append(np.zeros((len(points), 3, 3), dtype=float))
                layout.append((mode, "pressure", "real", index))
                precision.append("zero_velocity_pressure")
                blocks_v.append(np.zeros((len(points), 3), dtype=float))
                blocks_j.append(np.zeros((len(points), 3, 3), dtype=float))
                layout.append((mode, "pressure", "imag", index))
                precision.append("zero_velocity_pressure")
    velocity = np.stack(blocks_v, axis=-1)
    gradient = np.stack(blocks_j, axis=-1)
    if velocity.shape != (len(points), 3, 180) or gradient.shape != (len(points), 3, 3, 180):
        raise ValueError(f"Unexpected analytic response shapes {velocity.shape} {gradient.shape}")
    return velocity, gradient, layout, precision


def _observables_and_jacobian(velocity, gradient, weights, points, control_jacobian_v=None,
                              control_jacobian_j=None):
    """Compute [radial RMS, aspect, signed weighted spin] and derivatives."""
    velocity = np.asarray(velocity, dtype=float)
    gradient = np.asarray(gradient, dtype=float)
    weights = np.asarray(weights, dtype=float)
    points = np.asarray(points, dtype=float)
    radius = np.linalg.norm(points[:, :2], axis=1)
    radius2 = radius * radius
    if np.any(radius <= 0.0):
        raise ValueError("Fixed cylinder must not include the axis")
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
        raise ValueError("Nonpositive or nonfinite endpoint enstrophy")
    angular_speed = (-points[:, 1] * velocity[:, 0] + points[:, 0] * velocity[:, 1]) / radius2
    zmean = float(weights @ (density * points[:, 2]) / energy)
    radius_squared = float(weights @ (density * radius2) / energy)
    axial_offset = points[:, 2] - zmean
    axial_squared = float(weights @ (density * axial_offset * axial_offset) / energy)
    radial_rms = float(np.sqrt(radius_squared))
    axial_rms = float(np.sqrt(axial_squared))
    spin = float(weights @ (density * angular_speed) / energy)
    values = np.array((radial_rms, axial_rms / radial_rms, spin), dtype=float)
    if control_jacobian_v is None:
        return values, None, dict(
            energy=energy, omega=omega, density=density,
            angular_speed=angular_speed, zmean=zmean,
            radius_squared=radius_squared, axial_squared=axial_squared,
            radial_rms=radial_rms, axial_rms=axial_rms, spin=spin,
        )

    dv = np.asarray(control_jacobian_v, dtype=float)
    dgrad = np.asarray(control_jacobian_j, dtype=float)
    domega = _curl(dgrad)
    ddensity = 2.0 * np.sum(omega[:, :, None] * domega, axis=1)
    denergy = weights @ ddensity
    d_radius_squared = weights @ ((radius2 - radius_squared)[:, None] * ddensity) / energy
    d_radial = d_radius_squared / (2.0 * radial_rms)
    d_zmean = weights @ ((points[:, 2] - zmean)[:, None] * ddensity) / energy
    # The d zmean term cancels because the density-weighted centered first
    # moment is zero, leaving the compact centered formula below.
    d_axial_squared = (
        weights @ ((axial_offset * axial_offset - axial_squared)[:, None] * ddensity)
        / energy
    )
    d_axial = d_axial_squared / (2.0 * axial_rms)
    d_aspect = d_axial / radial_rms - axial_rms * d_radial / radial_rms**2
    d_angular = (-points[:, 1, None] * dv[:, 0, :] + points[:, 0, None] * dv[:, 1, :]) / radius2[:, None]
    d_spin = weights @ (
        (angular_speed - spin)[:, None] * ddensity + density[:, None] * d_angular
    ) / energy
    jacobian = np.vstack((d_radial, d_aspect, d_spin))
    detail = dict(
        energy=energy, omega=omega, density=density,
        angular_speed=angular_speed, zmean=zmean,
        radius_squared=radius_squared, axial_squared=axial_squared,
        radial_rms=radial_rms, axial_rms=axial_rms, spin=spin,
        d_energy=denergy, d_zmean=d_zmean,
    )
    return values, jacobian, detail


class EndpointShape:
    """Load a frozen endpoint cache and evaluate affine tangent observables."""

    def __init__(self, cache_path=DEFAULT_CACHE):
        self.cache_path = Path(cache_path)
        with np.load(self.cache_path, allow_pickle=False) as loaded:
            data = {key: loaded[key] for key in loaded.files}
        self.points = np.asarray(data["points"], dtype=float)
        self.weights = np.asarray(data["weights"], dtype=float)
        self.velocity_offset = np.asarray(data["velocity_offset"], dtype=float)
        self.gradient_offset = np.asarray(data["gradient_offset"], dtype=float)
        self.velocity_control = np.asarray(data["velocity_control"], dtype=float)
        self.gradient_control = np.asarray(data["gradient_control"], dtype=float)
        self.reference = np.asarray(data["reference"], dtype=float)
        self.coefficients_original = np.asarray(data["coefficients_original"], dtype=float)
        self.k0 = float(np.asarray(data["k0"]).ravel()[0])
        self.k_endpoint = float(np.asarray(data["k_endpoint"]).ravel()[0])
        self.tau0 = float(np.asarray(data["tau0"]).ravel()[0])
        self.tau_endpoint = float(np.asarray(data["tau_endpoint"]).ravel()[0])
        self.dtau = float(np.asarray(data["dtau"]).ravel()[0])
        self.point_count = int(len(self.points))
        self.control_count = int(self.velocity_control.shape[-1])
        if self.velocity_control.shape != (self.point_count, 3, 180):
            raise ValueError("Endpoint velocity control shape is not (N,3,180)")
        if self.gradient_control.shape != (self.point_count, 3, 3, 180):
            raise ValueError("Endpoint gradient control shape is not (N,3,3,180)")

    def evaluate(self, control):
        """Return endpoint values and 3 by 180 Jacobian for a tangent control."""
        control = np.asarray(control, dtype=float)
        if control.shape != (180,):
            raise ValueError(f"control must have shape (180,), got {control.shape}")
        velocity = self.velocity_offset + np.einsum("ncq,q->nc", self.velocity_control, control)
        gradient = self.gradient_offset + np.einsum("ncdq,q->ncd", self.gradient_control, control)
        values, jacobian, _ = _observables_and_jacobian(
            velocity, gradient, self.weights, self.points,
            self.velocity_control, self.gradient_control
        )
        return values, jacobian


def _saved_endpoint(path, k_endpoint):
    if not path.exists():
        return None
    report = json.loads(path.read_text(encoding="utf-8"))
    if report.get("status") != "completed":
        return None
    rows = report.get("rows", [])
    for row in rows:
        if abs(float(row.get("k", np.nan)) - k_endpoint) <= 2.0e-12:
            return np.array((row["enstrophy_radial_rms"], row["enstrophy_aspect_ratio"],
                             row["enstrophy_weighted_angular_speed"]), dtype=float)
    return None


def build_cache(cache_path=DEFAULT_CACHE, report_path=DEFAULT_REPORT,
                candidate_path=DEFAULT_CANDIDATE):
    started = time.perf_counter()
    candidate_path = Path(candidate_path)
    candidate_raw = candidate_path.read_bytes()
    candidate = json.loads(candidate_raw)
    if candidate.get("status") != "completed":
        raise ValueError("Tangent candidate must be completed")
    shape_raw = SHAPE_ROWS.read_bytes()
    shape_rows = json.loads(shape_raw)
    frozen_raw = FROZEN_GEOMETRY.read_bytes()
    frozen = json.loads(frozen_raw)
    geometry = frozen["inputs"]["wave"]
    mean_input = frozen["inputs"]["mean"]
    mean, mean_report = load_saved_field()
    install_in_field(mean)
    if mean_report.get("coefficients") != mean_input.get("coefficients"):
        raise ValueError("Loaded mean coefficients differ from frozen geometry")
    dynamic, dynamic_report = load_dynamic()
    install_in_field(dynamic)
    base, *_ = value_and_pressure_modes(dynamic, dynamic_report)
    k0 = float(mean_input["k"])
    tau0 = float(mean_input["tau"])
    if abs(k0 - float(dynamic_report["k"])) > 1.0e-12:
        raise ValueError("Frozen mean and dynamic k differ")
    points, weights, domain = fixed_cylinder(base.inner, k0, angles=12)
    delta_k = DELTA_K
    k_endpoint = k0 + delta_k
    tau_endpoint = 0.5 * 2.0 ** (-k_endpoint)
    dtau = tau0 - tau_endpoint
    mean_velocity, mean_gradient, hspace = _mean_endpoint_fd(mean, points, tau_endpoint)

    packed_wave = np.asarray(candidate["selected"]["coefficients_original"], dtype=float)
    if packed_wave.shape != (27, 2):
        raise ValueError("Expected packed 27-column initial wave")
    wave = packed_wave[:, 0] + 1j * packed_wave[:, 1]
    carriers = {0: np.zeros(2), 1: np.asarray(geometry["carrier"], dtype=float),
                2: 2.0 * np.asarray(geometry["carrier"], dtype=float)}
    initial_velocity = np.zeros((len(points), 3), dtype=float)
    initial_gradient = np.zeros((len(points), 3, 3), dtype=float)
    for mode, coefficient in ((1, wave),):
        values, gradients, _, _, _ = basis_jets(
            points, geometry["center"], geometry["widths"], mode,
            int(geometry["degree"]), carriers[mode], NU
        )
        initial_velocity += np.einsum("ncq,q->nc", values, coefficient).real
        initial_gradient += np.einsum("ncdq,q->ncd", gradients, coefficient).real
    velocity_control, gradient_control, layout, precision = _basis_response(
        points, geometry, dtau
    )
    velocity_offset = mean_velocity + initial_velocity
    gradient_offset = mean_gradient + initial_gradient
    reference = np.asarray(shape_rows["reference_observables"], dtype=float)
    zero_values, zero_jacobian, zero_detail = _observables_and_jacobian(
        velocity_offset, gradient_offset, weights, points,
        velocity_control, gradient_control
    )
    endpoint_shape = EndpointShape.__new__(EndpointShape)
    endpoint_shape.points = points
    endpoint_shape.weights = weights
    endpoint_shape.velocity_offset = velocity_offset
    endpoint_shape.gradient_offset = gradient_offset
    endpoint_shape.velocity_control = velocity_control
    endpoint_shape.gradient_control = gradient_control
    endpoint_shape.reference = reference
    endpoint_shape.coefficients_original = packed_wave
    endpoint_shape.k0 = k0
    endpoint_shape.k_endpoint = k_endpoint
    endpoint_shape.tau0 = tau0
    endpoint_shape.tau_endpoint = tau_endpoint
    endpoint_shape.dtau = dtau
    endpoint_shape.point_count = len(points)
    endpoint_shape.control_count = 180

    original_candidate = json.loads((ROOT / "wave_moment_cone_tangent.json").read_text(encoding="utf-8"))
    original_control = np.asarray(original_candidate["selected"]["tangent_coefficients"], dtype=float)
    shape_control = np.asarray(candidate["selected"]["tangent_coefficients"], dtype=float)
    original_values, original_jacobian = endpoint_shape.evaluate(original_control)
    shape_values, shape_jacobian = endpoint_shape.evaluate(shape_control)
    saved_original = _saved_endpoint(DEFAULT_SAVED_ORIGINAL, k_endpoint)
    saved_shape = _saved_endpoint(DEFAULT_SAVED_SHAPE, k_endpoint)

    # One directional derivative check in scaled control units. The endpoint
    # map itself is affine, so this isolates the nonlinear observable formulas.
    scales = np.maximum(np.abs(shape_control), 1.0)
    direction = np.linspace(1.0, 2.0, 180) * scales
    direction /= np.linalg.norm(direction)
    # The spin directional derivative is about 1e-6 in these normalized
    # control units; a 1e-3 step loses several digits to subtraction from a
    # 7e5 baseline.  This larger relative-small step keeps the central
    # difference in the affine-control regime while resolving the signal.
    fd_step = 1.0e3
    plus, _ = endpoint_shape.evaluate(shape_control + fd_step * direction)
    minus, _ = endpoint_shape.evaluate(shape_control - fd_step * direction)
    fd_directional = (plus - minus) / (2.0 * fd_step)
    analytic_directional = shape_jacobian @ direction
    gradient_error = float(np.max(np.abs(fd_directional - analytic_directional)))
    gradient_rel_error = gradient_error / max(float(np.max(np.abs(fd_directional))), 1.0e-30)

    cache_path = Path(cache_path)
    np.savez_compressed(
        cache_path,
        points=points,
        weights=weights,
        velocity_offset=velocity_offset,
        gradient_offset=gradient_offset,
        velocity_control=velocity_control,
        gradient_control=gradient_control,
        reference=reference,
        coefficients_original=packed_wave,
        k0=np.array([k0]), k_endpoint=np.array([k_endpoint]),
        tau0=np.array([tau0]), tau_endpoint=np.array([tau_endpoint]),
        dtau=np.array([dtau]),
        source_hash_candidate=np.array([_sha256(candidate_path)]),
        source_hash_geometry=np.array([_sha256(FROZEN_GEOMETRY)]),
        source_hash_mean=np.array([mean_input["source_sha256"]]),
        source_hash_shape_rows=np.array([_sha256(SHAPE_ROWS)]),
    )
    report = dict(
        status="completed", accepted=False, pde_validated=False,
        scale_recursion_established=False, wave_integrated=False,
        source=candidate_path.name, source_sha256=_sha256(candidate_path),
        geometry_source=FROZEN_GEOMETRY.name, geometry_sha256=_sha256(FROZEN_GEOMETRY),
        shape_rows_source=SHAPE_ROWS.name, shape_rows_sha256=_sha256(SHAPE_ROWS),
        mean_source=mean_input["source_report"], mean_source_sha256=mean_input["source_sha256"],
        cache_path=cache_path.name, cache_sha256=_sha256(cache_path),
        reference_k=k0, endpoint_k=k_endpoint,
        physical_time_reference=-tau0, physical_time_endpoint=-tau_endpoint,
        tau_reference=tau0, tau_endpoint=tau_endpoint, dtau=dtau,
        finite_difference_hspace=hspace, point_count=int(len(points)),
        domain=domain, control_count=180, pressure_column_count=45,
        precision_convention={
            "mean_offset": "double precision mean.fields with centered five-point spatial FD",
            "initial_wave_offset": "double precision supported_fourier_analytic_jets",
            "tangent_velocity_gradient": "double precision supported_fourier_analytic_jets, multiplied by tau_reference-tau_endpoint",
            "pressure_columns": "exact zero velocity and gradient",
        },
        control_layout=[dict(mode=a, kind=b, component=c, index=d, precision=p)
                        for (a, b, c, d), p in zip(layout, precision)],
        reference_observables=reference.tolist(),
        zero_control_endpoint_observables=zero_values.tolist(),
        zero_control_detail={k: float(v) for k, v in zero_detail.items()
                             if np.isscalar(v)},
        candidate_comparisons={
            "original_wave_moment_cone": {
                "control_source": "wave_moment_cone_tangent.json",
                "endpoint_observables": original_values.tolist(),
                "endpoint_jacobian_max_abs": float(np.max(np.abs(original_jacobian))),
                "saved_fd_endpoint_observables": saved_original.tolist() if saved_original is not None else None,
                "analytic_minus_saved_fd": (original_values - saved_original).tolist() if saved_original is not None else None,
            },
            "shape_direction_margin": {
                "control_source": candidate_path.name,
                "endpoint_observables": shape_values.tolist(),
                "endpoint_jacobian_max_abs": float(np.max(np.abs(shape_jacobian))),
                "saved_fd_endpoint_observables": saved_shape.tolist() if saved_shape is not None else None,
                "analytic_minus_saved_fd": (shape_values - saved_shape).tolist() if saved_shape is not None else None,
            },
        },
        directional_gradient_check={
            "control_source": candidate_path.name, "step": fd_step,
            "analytic": analytic_directional.tolist(),
            "finite_difference": fd_directional.tolist(),
            "max_abs_error": gradient_error,
            "relative_error": gradient_rel_error,
        },
        scope=(
            "Fixed physical cylinder endpoint oracle. Mean endpoint velocity and "
            "Jacobian use the same five-point spatial stencil as saved observables; "
            "wave offset and all tangent velocity/J columns use analytic supported "
            "Fourier jets. Endpoint shape only; no PDE, trajectory, or scale-recursion claim."
        ),
        elapsed_seconds=float(time.perf_counter() - started),
    )
    report_path = Path(report_path)
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"], "point_count": report["point_count"],
        "endpoint_k": report["endpoint_k"],
        "original_endpoint": original_values.tolist(),
        "shape_endpoint": shape_values.tolist(),
        "gradient_relative_error": gradient_rel_error,
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    args = parser.parse_args()
    build_cache(args.cache, args.output, args.candidate)
