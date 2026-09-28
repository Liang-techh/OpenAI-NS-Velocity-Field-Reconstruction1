"""Infinitesimal physical-time generator for the solenoidal scale transport.

The frozen transport in ``solenoidal_scale_transport.py`` uses
``s = tau / tau0`` and physical time ``t = -tau``.  Therefore
``ds/dt = -1/tau0`` at ``s=1``.  If ``J = grad u_ref`` and ``m`` is the
axisymmetric reference swirl, differentiation of the Piola map gives

    u_t = [diag(1/2, 1/2, 1/2+h) u
           + J @ (1/2*x, 1/2*y, (1/2-h)*z)
           + h*m*e_theta] / tau0.

This module only supplies and checks that kinematic generator.  It does not
add pressure, a forcing term, or a Navier--Stokes residual.
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
from solenoidal_scale_transport import (  # noqa: E402
    MEAN_ANGLES,
    SolenoidalScaleTransport,
    VelocityAdapter,
)


GRID_PATH = ROOT / "enriched_endpoint_shape_cache.npz"
TRANSPORT_MODULE = ROOT / "solenoidal_scale_transport.py"
TWO_PATCH_MODULE = ROOT / "global_two_patch_candidate.py"
TWO_PATCH_SOURCE = ROOT / "localized_two_patch_constrained.json"
PARENT_MODULE = ROOT / "global_localized_candidate.py"


def _as_points(points):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    return points


def _as_vectors(name, values, count):
    values = np.asarray(values, dtype=float)
    if values.shape != (count, 3):
        raise ValueError(f"{name} must have shape ({count},3), got {values.shape}")
    return values


def _as_jacobian(jacobian, count):
    jacobian = np.asarray(jacobian, dtype=float)
    if jacobian.shape == (3, 3):
        jacobian = np.broadcast_to(jacobian, (count, 3, 3)).copy()
    if jacobian.shape != (count, 3, 3):
        raise ValueError(f"jacobian must have shape ({count},3,3) or (3,3), got {jacobian.shape}")
    return jacobian


def _axis_direction(points):
    points = _as_points(points)
    radius = np.linalg.norm(points[:, :2], axis=1)
    result = np.zeros_like(points)
    active = radius > 0.0
    result[active, 0] = -points[active, 1] / radius[active]
    result[active, 1] = points[active, 0] / radius[active]
    return result


def axisymmetric_swirl_mean(reference, points, mean_angles=MEAN_ANGLES):
    """Return the discrete absolute-angle mean of reference ``u_theta``.

    The fixed angle stencil makes this a function of ``(r,z)`` alone, even if
    a future reference contains angular modes above the resolved bandwidth.
    At the axis the cylindrical direction is undefined, so the regular value
    used by the correction is zero.
    """

    points = _as_points(points)
    mean_angles = int(mean_angles)
    if mean_angles < 3 or mean_angles % 2 == 0:
        raise ValueError("mean_angles must be odd and at least 3")
    radius = np.linalg.norm(points[:, :2], axis=1)
    result = np.zeros(len(points), dtype=float)
    active = radius > 0.0
    if not np.any(active):
        return result
    active_points = points[active]
    angles = np.broadcast_to(
        2.0 * np.pi * np.arange(mean_angles)[None, :] / mean_angles,
        (len(active_points), mean_angles),
    )
    rotated = np.empty((len(active_points), mean_angles, 3), dtype=float)
    rotated[:, :, 0] = radius[active, None] * np.cos(angles)
    rotated[:, :, 1] = radius[active, None] * np.sin(angles)
    rotated[:, :, 2] = active_points[:, 2, None]
    values = np.asarray(reference.velocity(rotated.reshape(-1, 3)), dtype=float).reshape(
        len(active_points), mean_angles, 3
    )
    swirl = -np.sin(angles) * values[:, :, 0] + np.cos(angles) * values[:, :, 1]
    result[active] = np.mean(swirl, axis=1)
    return result


def solenoidal_scale_generator(
    points,
    velocity,
    jacobian,
    tau0,
    h,
    mean_swirl,
):
    """Evaluate the physical-time derivative at ``s=1``.

    ``jacobian`` is the Cartesian spatial Jacobian with rows indexed by output
    component and columns by input coordinate.  Supplying it explicitly keeps
    the generator usable with analytic or cached spatial derivatives.
    """

    points = _as_points(points)
    count = len(points)
    velocity = _as_vectors("velocity", velocity, count)
    jacobian = _as_jacobian(jacobian, count)
    mean_swirl = np.asarray(mean_swirl, dtype=float)
    if mean_swirl.shape != (count,):
        raise ValueError(f"mean_swirl must have shape ({count},), got {mean_swirl.shape}")
    tau0 = float(tau0)
    h = float(h)
    if not np.isfinite(tau0) or tau0 <= 0.0:
        raise ValueError("tau0 must be positive and finite")
    if not np.isfinite(h):
        raise ValueError("h must be finite")
    source_direction = points * np.array([0.5, 0.5, 0.5 - h], dtype=float)[None, :]
    advective = np.einsum("nij,nj->ni", jacobian, source_direction)
    amplitude = velocity * np.array([0.5, 0.5, 0.5 + h], dtype=float)[None, :]
    correction = h * mean_swirl[:, None] * _axis_direction(points)
    return (amplitude + advective + correction) / tau0


def generator_from_reference(reference, points, tau0, h, jacobian, mean_angles=MEAN_ANGLES):
    """Convenience wrapper using a velocity reference and supplied ``J``."""

    points = _as_points(points)
    velocity = np.asarray(reference.velocity(points), dtype=float)
    mean_swirl = axisymmetric_swirl_mean(reference, points, mean_angles)
    return solenoidal_scale_generator(points, velocity, jacobian, tau0, h, mean_swirl)


def _finite_difference_jacobian(reference, points, step):
    """Fourth-order Cartesian Jacobian for a small bounded validation set."""

    points = _as_points(points)
    count = len(points)
    step = np.broadcast_to(np.asarray(step, dtype=float), (3,))
    if np.any(step <= 0.0) or not np.all(np.isfinite(step)):
        raise ValueError("step must be positive and finite")
    jacobian = np.empty((count, 3, 3), dtype=float)
    for axis in range(3):
        offset = np.zeros(3, dtype=float)
        offset[axis] = step[axis]
        plus = reference.velocity(np.concatenate((points + offset, points + 2.0 * offset), axis=0))
        minus = reference.velocity(np.concatenate((points - offset, points - 2.0 * offset), axis=0))
        n = count
        jacobian[:, :, axis] = (
            8.0 * (plus[:n] - minus[:n]) - (plus[n:] - minus[n:])
        ) / (12.0 * step[axis])
    return jacobian


def _norm(values):
    values = np.asarray(values, dtype=float)
    return float(np.linalg.norm(values))


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _load_validation_points(localized, tau0):
    inner_radius, outer_radius = localized.radial_radii(tau0)
    transition_radius = 0.5 * (inner_radius + outer_radius)
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 2.0e-4],
            [3.0e-3, 0.0, 0.0],
            [4.0e-3, 1.0e-3, 1.0e-4],
            [transition_radius, 0.0, 0.0],
            [transition_radius, 0.0, 1.0e-4],
            [outer_radius * 1.02, 0.0, 0.0],
            [1.2e-2, 0.0, 0.0],
            [0.0, 0.0, 1.0e-3],
        ],
        dtype=float,
    )
    labels = ["axis", "axis", "interior", "interior", "transition", "transition", "exterior", "exterior", "exterior"]
    return points, labels, float(inner_radius), float(outer_radius)


def _category_errors(labels, errors):
    output = {}
    for label in sorted(set(labels)):
        values = np.asarray([error for item, error in zip(labels, errors) if item == label], dtype=float)
        output[label] = {"count": int(len(values)), "max_abs": float(np.max(values)), "rms": float(np.sqrt(np.mean(values * values)))}
    return output


def run(output_path=ROOT / "scale_transport_generator.json"):
    started = time.perf_counter()
    full_field, localized, snapshot, candidate, hashes = load_two_patch()
    tau0 = float(snapshot["inputs"]["mean"]["tau"])
    h = float(localized.inner.h)
    reference = VelocityAdapter(full_field, tau0)
    points, labels, inner_radius, outer_radius = _load_validation_points(localized, tau0)
    nu = float(snapshot["inputs"]["mean"]["nu"])
    fd_step = 5.0e-4 * np.sqrt(nu * tau0)
    jacobian_coarse = _finite_difference_jacobian(reference, points, fd_step)
    jacobian_fine = _finite_difference_jacobian(reference, points, fd_step / 2.0)
    velocity = reference.velocity(points)
    mean_swirl = axisymmetric_swirl_mean(reference, points, MEAN_ANGLES)
    generator = solenoidal_scale_generator(points, velocity, jacobian_fine, tau0, h, mean_swirl)
    generator_coarse = solenoidal_scale_generator(points, velocity, jacobian_coarse, tau0, h, mean_swirl)
    report = {
        "status": "initializing",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "scope": (
            "Infinitesimal physical-time derivative of the frozen solenoidal scale map at s=1. "
            "The reference is held at tau0; pressure, Navier--Stokes dynamics, and recursive closure are omitted."
        ),
        "sources": {
            "two_patch_module": {"path": TWO_PATCH_MODULE.name, "sha256": _sha256(TWO_PATCH_MODULE)},
            "two_patch_source": {"path": TWO_PATCH_SOURCE.name, "sha256": _sha256(TWO_PATCH_SOURCE)},
            "parent_module": {"path": PARENT_MODULE.name, "sha256": _sha256(PARENT_MODULE)},
            "transport_module": {"path": TRANSPORT_MODULE.name, "sha256": _sha256(TRANSPORT_MODULE)},
            "fixed_cylinder_cache": {"path": GRID_PATH.name, "sha256": _sha256(GRID_PATH)},
            "loader_hashes": hashes,
        },
        "inputs": {
            "reference_k": float(snapshot["inputs"]["mean"]["k"]),
            "reference_tau": tau0,
            "physical_time_convention": "t=-tau",
            "h": h,
            "nu": nu,
            "mean_angles": MEAN_ANGLES,
            "inner_radius": inner_radius,
            "outer_radius": outer_radius,
            "finite_difference_jacobian_step": fd_step,
        },
        "derivation": {
            "scale_time_relation": "s=tau/tau0, ds/dt=-1/tau0 at s=1",
            "formula": "ut=(diag(1/2,1/2,1/2+h)*u + J@[x/2,y/2,(1/2-h)z] + h*m*e_theta)/tau0",
            "piola_sign": "d/ds [s^(-1/2) u_h(S^-1x)] = -u_h/2 - J_h@[x/2,y/2,(1/2-h)z]; ds/dt is negative",
            "swirl_sign": "d/ds [s^(-1/2-h)-s^(-1/2)] at s=1 is -h; ds/dt gives +h/tau0",
            "divergence_scope": "The Piola part transports divergence by s^(-1); the added axisymmetric swirl has zero divergence for r>0.",
        },
        "sample_points": {"labels": labels, "points": points.tolist()},
    }
    _save = lambda: Path(output_path).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    _save()
    print(json.dumps({"stage": "generator_assembled", "point_count": len(points)}), flush=True)

    jacobian_difference = jacobian_fine - jacobian_coarse
    generator_difference = generator - generator_coarse
    report["jacobian_refinement"] = {
        "coarse_step": fd_step,
        "fine_step": fd_step / 2.0,
        "max_abs_difference": float(np.max(np.abs(jacobian_difference))),
        "relative_l2_difference": _norm(jacobian_difference) / max(_norm(jacobian_fine), 1.0e-300),
        "category": _category_errors(labels, np.linalg.norm(jacobian_difference, axis=(1, 2))),
    }
    report["generator_jacobian_refinement"] = {
        "max_abs_difference": float(np.max(np.abs(generator_difference))),
        "relative_l2_difference": _norm(generator_difference) / max(_norm(generator), 1.0e-300),
        "category": _category_errors(labels, np.linalg.norm(generator_difference, axis=1)),
    }
    _save()

    # A physical-time forward difference uses s=1-dt/tau0 because advancing
    # t=-tau decreases tau.  Several dimensionless dt/tau0 values expose the
    # expected first-order truncation without changing the reference field.
    forward_rows = []
    for epsilon in (2.0e-3, 1.0e-3, 5.0e-4, 2.5e-4):
        dt = epsilon * tau0
        transport = SolenoidalScaleTransport(reference, 1.0 - epsilon, h, MEAN_ANGLES)
        finite_difference = (transport.velocity(points) - velocity) / dt
        error = finite_difference - generator
        row = {
            "dt_over_tau0": epsilon,
            "dt": dt,
            "max_abs_error": float(np.max(np.abs(error))),
            "weighted_rms_error": float(np.sqrt(np.mean(np.sum(error * error, axis=1)))),
            "relative_l2_error": _norm(error) / max(_norm(generator), 1.0e-300),
            "category": _category_errors(labels, np.linalg.norm(error, axis=1)),
        }
        forward_rows.append(row)
        print(json.dumps({"stage": "forward_difference", **row}), flush=True)
    for coarse, fine in zip(forward_rows, forward_rows[1:]):
        coarse["error_ratio_to_next_finer"] = coarse["max_abs_error"] / max(fine["max_abs_error"], 1.0e-300)
    report["forward_physical_time_check"] = {
        "rows": forward_rows,
        "interpretation": (
            "Forward physical-time differences of T_(1-dt/tau0) at fixed reference tau0. "
            "The expected first-order error decreases with dt; this validates the sign and amplitudes of the generator only."
        ),
    }
    report["status"] = "completed"
    report["elapsed_seconds"] = float(time.perf_counter() - started)
    _save()
    print(json.dumps({"stage": "completed", "elapsed_seconds": report["elapsed_seconds"], "forward": forward_rows}), flush=True)
    return report


if __name__ == "__main__":
    run()
