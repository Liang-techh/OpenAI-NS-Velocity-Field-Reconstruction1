"""Divergence-preserving anisotropic scale transport at one frozen time.

For ``S = diag(s**.5, s**.5, s**(.5-h))`` this module transports a reference
velocity by the Piola map

    u_P(x) = s**(.5-h) * S / det(S) * u_ref(S**(-1) x),

then restores the paper's stronger axisymmetric swirl scaling with a pure
axisymmetric swirl increment.  The latter has zero divergence independently
of its radial and axial profile.  This is a kinematic mapping diagnostic only;
pressure, time evolution, and the Navier--Stokes residual are deliberately
omitted.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

for _key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from global_two_patch_candidate import load as load_two_patch  # noqa: E402


GRID_PATH = ROOT / "enriched_endpoint_shape_cache.npz"
TWO_PATCH_MODULE = ROOT / "global_two_patch_candidate.py"
TWO_PATCH_SOURCE = ROOT / "localized_two_patch_constrained.json"
PARENT_MODULE = ROOT / "global_localized_candidate.py"
DELTA_K = 1.0e-6
S_HALF = 0.5
MEAN_ANGLES = 5


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(path: Path, report: dict) -> None:
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def cylindrical(velocity, points):
    points = np.asarray(points, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    theta = np.arctan2(points[:, 1], points[:, 0])
    c, s = np.cos(theta), np.sin(theta)
    return np.column_stack((
        c * velocity[:, 0] + s * velocity[:, 1],
        -s * velocity[:, 0] + c * velocity[:, 1],
        velocity[:, 2],
    ))


class VelocityAdapter:
    """Expose a saved field or transport as a common velocity-only object."""

    def __init__(self, field, tau):
        self.field = field
        self.tau = float(tau)

    def velocity(self, points):
        points = np.asarray(points, dtype=float)
        return np.asarray(self.field.fields(points, self.tau)[0], dtype=float)


class SolenoidalScaleTransport:
    """Piola transport plus mean-swirl correction for one scale ``s``."""

    def __init__(self, reference, scale, h, mean_angles=MEAN_ANGLES):
        self.reference = reference
        self.scale = float(scale)
        self.h = float(h)
        self.mean_angles = int(mean_angles)
        if not 0.0 < self.scale <= 1.0:
            raise ValueError("scale must lie in (0,1]")
        if self.mean_angles < 3 or self.mean_angles % 2 == 0:
            raise ValueError("mean_angles must be odd and at least 3")
        self.spatial_scale = np.array(
            [self.scale ** 0.5, self.scale ** 0.5, self.scale ** (0.5 - self.h)],
            dtype=float,
        )
        self.det_scale = float(self.scale ** (1.5 - self.h))
        self.piola_amplitude = float(self.scale ** (0.5 - self.h))
        self.horizontal_factor = float(self.scale ** -0.5)
        self.axial_factor = float(self.scale ** (-0.5 - self.h))
        self.swirl_increment_factor = self.axial_factor - self.horizontal_factor

    def _axisymmetric_swirl(self, source_points):
        source_points = np.asarray(source_points, dtype=float)
        radius = np.linalg.norm(source_points[:, :2], axis=1)
        result = np.zeros(len(source_points), dtype=float)
        active = radius > 0.0
        if not np.any(active):
            return result
        active_points = source_points[active]
        # Use fixed absolute angles rather than a theta-offset stencil.  The
        # discrete average is then a function of (r,z) only even when the
        # sampled field contains unresolved angular modes; for the current
        # |m|<=2 field this agrees with every rotated stencil to roundoff.
        angles = np.broadcast_to(
            2.0 * np.pi * np.arange(self.mean_angles)[None, :] / self.mean_angles,
            (len(active_points), self.mean_angles),
        )
        rotated = np.empty((len(active_points), self.mean_angles, 3), dtype=float)
        rotated[:, :, 0] = radius[active, None] * np.cos(angles)
        rotated[:, :, 1] = radius[active, None] * np.sin(angles)
        rotated[:, :, 2] = active_points[:, 2, None]
        values = self.reference.velocity(rotated.reshape(-1, 3)).reshape(
            len(active_points), self.mean_angles, 3
        )
        c, s = np.cos(angles), np.sin(angles)
        swirl = -s * values[:, :, 0] + c * values[:, :, 1]
        result[active] = np.mean(swirl, axis=1)
        return result

    def velocity(self, points):
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("points must have shape (n,3)")
        source_points = points / self.spatial_scale[None, :]
        reference_velocity = self.reference.velocity(source_points)
        piola = self.piola_amplitude * self.spatial_scale[None, :] / self.det_scale * reference_velocity
        mean_swirl = self._axisymmetric_swirl(source_points)
        radius = np.linalg.norm(points[:, :2], axis=1)
        correction = np.zeros_like(points)
        active = radius > 0.0
        correction[active, 0] = -points[active, 1] / radius[active] * self.swirl_increment_factor * mean_swirl[active]
        correction[active, 1] = points[active, 0] / radius[active] * self.swirl_increment_factor * mean_swirl[active]
        return piola + correction


def _load_grid():
    with np.load(GRID_PATH, allow_pickle=False) as loaded:
        points = np.asarray(loaded["points"], dtype=float)
        weights = np.asarray(loaded["weights"], dtype=float)
    if points.shape != (7776, 3) or weights.shape != (7776,):
        raise ValueError(f"Expected (7776,3)/(7776,), got {points.shape}/{weights.shape}")
    return points, weights


def _select_probes(points, count=18):
    radius = np.linalg.norm(points[:, :2], axis=1)
    mask = (radius > 0.0025) & (radius < 0.0075) & (np.abs(points[:, 2]) < 0.00055)
    candidates = points[mask]
    if len(candidates) < count:
        raise ValueError(f"Only {len(candidates)} interior probe points available")
    # Keep one angular sector and deterministic z/r spread. The fixed cache
    # already carries the physical quadrature weights; this is only a bounded
    # pointwise derivative check.
    selected = candidates[:: max(1, len(candidates) // count)][:count]
    if len(selected) != count:
        selected = candidates[:count]
    return np.asarray(selected, dtype=float)


def _fd_divergence(evaluator, points, steps):
    points = np.asarray(points, dtype=float)
    steps = np.broadcast_to(np.asarray(steps, dtype=float), (3,))
    values = np.zeros(len(points), dtype=float)
    for axis in range(3):
        offset = np.zeros(3, dtype=float)
        offset[axis] = steps[axis]
        plus = evaluator(np.concatenate((points + offset, points + 2.0 * offset), axis=0))
        minus = evaluator(np.concatenate((points - offset, points - 2.0 * offset), axis=0))
        n = len(points)
        values += (8.0 * (plus[:n, axis] - minus[:n, axis]) - (plus[n:, axis] - minus[n:, axis])) / (12.0 * steps[axis])
    return values


def _weighted_norm(values, weights):
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.sum(np.asarray(weights) * np.sum(values * values, axis=-1))))


def _finite_summary(values):
    values = np.asarray(values, dtype=float)
    norms = np.linalg.norm(values, axis=1)
    return {
        "finite": bool(np.all(np.isfinite(values))),
        "nonzero_count": int(np.count_nonzero(norms > 0.0)),
        "max_norm": float(np.max(norms)),
        "min_norm": float(np.min(norms)),
    }


def _divergence_rows(reference, transport, source_points, target_points, scale, base_weights, nu):
    h0 = 5.0e-4 * np.sqrt(float(nu * 0.5 * 2.0 ** (-11.0003)))
    rows = []
    spatial = transport.spatial_scale
    for divisor in (1.0, 2.0):
        h = h0 / divisor
        output_steps = np.full(3, h, dtype=float)
        source_steps = h / spatial
        reference_divergence = _fd_divergence(reference.velocity, source_points, source_steps)
        transported_divergence = _fd_divergence(transport.velocity, target_points, output_steps)
        predicted = reference_divergence / scale
        difference = transported_divergence - predicted
        rows.append({
            "step": h,
            "reference_divergence_max_abs": float(np.max(np.abs(reference_divergence))),
            "transported_divergence_max_abs": float(np.max(np.abs(transported_divergence))),
            "predicted_transport_divergence_max_abs": float(np.max(np.abs(predicted))),
            "fd_minus_analytic_prediction_max_abs": float(np.max(np.abs(difference))),
            "fd_minus_analytic_prediction_weighted_rms": float(np.sqrt(np.average(difference * difference, weights=base_weights))),
        })
    rows[0]["error_reduction_h_to_h2"] = rows[0]["fd_minus_analytic_prediction_max_abs"] / max(rows[1]["fd_minus_analytic_prediction_max_abs"], 1.0e-300)
    return rows


def run(output_path=ROOT / "solenoidal_scale_transport.json"):
    started = time.perf_counter()
    points, weights = _load_grid()
    full_field, localized, snapshot, candidate, hashes = load_two_patch()
    mean = VelocityAdapter(full_field, snapshot["inputs"]["mean"]["tau"])
    h = float(localized.inner.h)
    tau0 = float(snapshot["inputs"]["mean"]["tau"])
    reference = mean
    transport = SolenoidalScaleTransport(reference, S_HALF, h, MEAN_ANGLES)
    source_probes = _select_probes(points)
    target_probes = source_probes * transport.spatial_scale[None, :]
    source_cyl = cylindrical(reference.velocity(source_probes), source_probes)
    target_cyl = cylindrical(transport.velocity(target_probes), target_probes)
    raw_target_cyl = source_cyl * np.array(
        [transport.horizontal_factor, transport.axial_factor, transport.axial_factor]
    )[None, :]
    outside = np.array([
        [0.02, 0.0, 0.0],
        [0.0, 0.02, 0.0],
        [0.0, 0.0, 0.003],
    ], dtype=float)
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "scope": (
            "Frozen-time kinematic anisotropic transport of the current global "
            "two-patch velocity. The Piola map preserves divergence under the "
            "anisotropic coordinate map, and a sampled axisymmetric swirl increment "
            "restores mean swirl scaling. Pressure and time dynamics are omitted; "
            "this is not PDE or recursive-closure acceptance."
        ),
        "sources": {
            "two_patch_module": {"path": TWO_PATCH_MODULE.name, "sha256": _sha256(TWO_PATCH_MODULE)},
            "two_patch_source": {"path": TWO_PATCH_SOURCE.name, "sha256": _sha256(TWO_PATCH_SOURCE)},
            "parent_module": {"path": PARENT_MODULE.name, "sha256": _sha256(PARENT_MODULE)},
            "fixed_cylinder_cache": {"path": GRID_PATH.name, "sha256": _sha256(GRID_PATH)},
            "loader_hashes": hashes,
        },
        "inputs": {
            "reference_k": 11.0003,
            "reference_tau": tau0,
            "h": h,
            "scale": S_HALF,
            "mean_angles": MEAN_ANGLES,
            "spatial_matrix": transport.spatial_scale.tolist(),
            "determinant": transport.det_scale,
            "piola_amplitude": transport.piola_amplitude,
            "horizontal_velocity_factor": transport.horizontal_factor,
            "axial_velocity_factor": transport.axial_factor,
            "swirl_increment_factor": transport.swirl_increment_factor,
            "pressure": "omitted by design",
        },
        "analytic_divergence_argument": {
            "piola_formula": "uP(x)=s^(0.5-h)*S/det(S)*u_ref(S^-1 x)",
            "piola_divergence": "div(uP)(x)=s^(-1)*div(u_ref)(S^-1 x)",
            "swirl_increment": "delta_u=delta_u_theta(r,z)*e_theta is axisymmetric, hence div(delta_u)=0 for r>0",
            "axis_mean_extraction": "five equispaced rotations; exact for angular Fourier modes |m|<=2, verified against twelve-angle samples",
        },
        "s_half_metric": {
            "probe_count": int(len(source_probes)),
            "source_velocity": _finite_summary(reference.velocity(source_probes)),
            "transport_velocity": _finite_summary(target_cyl),
            "reference_weighted_probe_norm": _weighted_norm(source_cyl, np.ones(len(source_cyl))),
            "raw_per_component_target_norm": _weighted_norm(raw_target_cyl, np.ones(len(raw_target_cyl))),
            "transport_target_norm": _weighted_norm(target_cyl, np.ones(len(target_cyl))),
            "transport_minus_raw_per_component_relative": float(
                _weighted_norm(target_cyl - raw_target_cyl, np.ones(len(target_cyl)))
                / max(_weighted_norm(raw_target_cyl, np.ones(len(raw_target_cyl))), 1.0e-300)
            ),
            "transported_support_probe_points": target_probes.tolist(),
        },
        "fields": {},
    }
    _save(Path(output_path), report)
    print(json.dumps({"stage": "sources_loaded", "probe_count": len(source_probes)}), flush=True)

    # Verify that five rotations recover the axisymmetric swirl mean at the
    # current angular bandwidth. A twelve-angle reference is an independent
    # check, not an input to the transport itself.
    alias_points = source_probes[: min(8, len(source_probes))]
    alias_r = np.linalg.norm(alias_points[:, :2], axis=1)
    alias_theta = np.arctan2(alias_points[:, 1], alias_points[:, 0])
    alias_rows = []
    for angle_count in (5, 12):
        angles = alias_theta[:, None] + 2.0 * np.pi * np.arange(angle_count)[None, :] / angle_count
        rotated = np.empty((len(alias_points), angle_count, 3), dtype=float)
        rotated[:, :, 0] = alias_r[:, None] * np.cos(angles)
        rotated[:, :, 1] = alias_r[:, None] * np.sin(angles)
        rotated[:, :, 2] = alias_points[:, 2, None]
        values = reference.velocity(rotated.reshape(-1, 3)).reshape(len(alias_points), angle_count, 3)
        swirl = -np.sin(angles) * values[:, :, 0] + np.cos(angles) * values[:, :, 1]
        alias_rows.append(np.mean(swirl, axis=1))
    alias_difference = alias_rows[0] - alias_rows[1]
    report["axisymmetric_mean_alias_check"] = {
        "probe_count": int(len(alias_points)),
        "angles_compared": [5, 12],
        "max_abs_difference": float(np.max(np.abs(alias_difference))),
        "rms_difference": float(np.sqrt(np.mean(alias_difference * alias_difference))),
        "interpretation": "Five-angle mean is exact for current |m|<=2 angular content up to evaluation roundoff; this is a sampled bandwidth check.",
    }
    _save(Path(output_path), report)

    # Divergence convergence at the mapped points. The source derivative steps
    # are divided by S so that both finite differences represent the same
    # physical output step.
    divergence_rows = _divergence_rows(reference, transport, source_probes, target_probes, S_HALF, np.ones(len(source_probes)), 0.01)
    report["divergence_fd_check"] = {
        "rows": divergence_rows,
        "convergence_scope": "bounded 4th-order Cartesian finite differences on interior probes; not a continuum or 1e-3 certificate",
    }
    _save(Path(output_path), report)
    print(json.dumps({"stage": "divergence_checked", "rows": divergence_rows}), flush=True)

    transported_outside = transport.velocity(outside)
    report["support_and_finiteness"] = {
        "reference_probe": _finite_summary(reference.velocity(source_probes)),
        "transport_probe": _finite_summary(transport.velocity(target_probes)),
        "outside_points": outside.tolist(),
        "outside_transport_max_abs": float(np.max(np.abs(transported_outside))),
        "interpretation": "The finite sampled reference is nonzero; the mapped support remains compact under S and these outside points evaluate to zero.",
    }
    _save(Path(output_path), report)

    # Composition check: T_s1(T_s2(u)) versus T_(s1*s2)(u) at a few mapped
    # points. The five-angle mean is recomputed at each nested stage.
    s1, s2 = 0.5, 0.8
    stage2 = SolenoidalScaleTransport(reference, s2, h, MEAN_ANGLES)
    sequential = SolenoidalScaleTransport(stage2, s1, h, MEAN_ANGLES)
    direct = SolenoidalScaleTransport(reference, s1 * s2, h, MEAN_ANGLES)
    composition_source = source_probes[: min(6, len(source_probes))]
    composition_points = composition_source * direct.spatial_scale[None, :]
    sequential_values = sequential.velocity(composition_points)
    direct_values = direct.velocity(composition_points)
    composition_difference = sequential_values - direct_values
    report["composition_check"] = {
        "s1": s1,
        "s2": s2,
        "product": s1 * s2,
        "probe_count": int(len(composition_points)),
        "max_abs_difference": float(np.max(np.abs(composition_difference))),
        "relative_l2_difference": float(
            np.linalg.norm(composition_difference)
            / max(np.linalg.norm(direct_values), 1.0e-300)
        ),
        "interpretation": "Numerical composition check for the kinematic map, including nested five-angle mean extraction.",
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save(Path(output_path), report)
    print(json.dumps({
        "stage": "completed",
        "axis_mean_alias_max": report["axisymmetric_mean_alias_check"]["max_abs_difference"],
        "divergence": report["divergence_fd_check"]["rows"],
        "composition_relative": report["composition_check"]["relative_l2_difference"],
        "elapsed_seconds": report["elapsed_seconds"],
    }), flush=True)
    return report


if __name__ == "__main__":
    run()
