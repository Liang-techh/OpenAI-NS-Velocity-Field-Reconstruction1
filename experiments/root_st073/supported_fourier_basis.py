"""Vectorized compact polynomial Fourier potential basis.

The columns returned by :func:`basis_data` are exact curls of compactly
supported cylindrical potentials.  The velocity columns use ``-curl(A)``;
this is the sign used by ``curl_wave_patch_collocation.basis`` for a time
derivative of a potential.  A scalar Fourier polynomial and its full
Cartesian pressure gradient are returned alongside the velocity columns.

The polynomial index is ``j = a * (degree + 1) + b`` for
``xi**a * eta**b``.  Potential components are packed in radial, theta,
axial blocks, each containing all ``j`` in that order.  Thus the velocity
column for component ``c`` and polynomial ``(a, b)`` is
``c * (degree + 1)**2 + j``.
"""

from __future__ import annotations

import json
import operator
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np


COMPONENT_ORDER = ("radial", "theta", "axial")
POLYNOMIAL_ORDER = "a-major then b: j = a * (degree + 1) + b"


def _integer(value, name):
    """Return an integer argument while rejecting non-integral values."""

    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc
    return int(result)


def basis_layout(degree=2):
    """Describe polynomial ordering and component packing for ``degree``."""

    degree = _integer(degree, "degree")
    if degree < 0:
        raise ValueError("degree must be nonnegative")
    n = degree + 1
    count = n * n
    return {
        "degree": degree,
        "polynomial_count": count,
        "polynomial_order": POLYNOMIAL_ORDER,
        "component_order": COMPONENT_ORDER,
        "component_blocks": {
            name: [offset * count, (offset + 1) * count]
            for offset, name in enumerate(COMPONENT_ORDER)
        },
        "pressure_block": [0, count],
    }


def _pair(value, name):
    try:
        pair = tuple(value)
    except TypeError as exc:
        raise TypeError(f"{name} must contain two values") from exc
    if len(pair) != 2:
        raise ValueError(f"{name} must contain two values")
    result = np.asarray(pair, dtype=float)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return float(result[0]), float(result[1])


def _points(points):
    result = np.asarray(points, dtype=float)
    if result.ndim == 1 and result.shape == (3,):
        result = result[None, :]
    if result.ndim != 2 or result.shape[1] != 3:
        raise ValueError("points must have shape (N, 3)")
    if not np.all(np.isfinite(result)):
        raise ValueError("points must be finite")
    return result


def _bump(value, center, halfwidth):
    """Vectorized form of ``curl_wave_prototype.bump`` and its derivative."""

    scaled = (value - center + halfwidth) / (2.0 * halfwidth)
    inside = (scaled > 0.0) & (scaled < 1.0)
    bump_value = np.zeros_like(value, dtype=float)
    bump_derivative = np.zeros_like(value, dtype=float)
    if np.any(inside):
        s = scaled[inside]
        one_minus = 1.0 - s
        bump_value[inside] = 1024.0 * s**5 * one_minus**5
        bump_derivative[inside] = (
            5120.0 * s**4 * one_minus**4 * (1.0 - 2.0 * s)
            / (2.0 * halfwidth)
        )
    return bump_value, bump_derivative


def basis_data(
    points,
    center=(0.004, 0.0),
    widths=(0.001, 0.001),
    mode=0,
    degree=2,
    carrier=(0.0, 0.0),
):
    """Return compact velocity, pressure, and pressure-gradient basis data.

    Parameters
    ----------
    points : array_like, shape (N, 3)
        Cartesian points ``(x, y, z)``.  A single length-three point is also
        accepted and is treated as one row.
    center : pair
        Cylindrical support center ``(r0, z0)``.
    widths : pair
        Radial and axial half-widths ``(dr, dz)``.  ``dr < r0`` is required so
        the compact support stays away from the cylindrical axis.
    mode : int
        Integer angular Fourier mode ``m``.
    degree : int
        Maximum power in each of ``xi`` and ``eta``.
    carrier : pair
        Radial and axial carrier wave numbers ``(kr, kz)``.

    Returns
    -------
    velocity : complex ndarray, shape (N, 3, 3 * M)
        Cartesian velocity columns, where ``M = (degree + 1)**2``.
    pressure : complex ndarray, shape (N, M)
        Scalar Fourier polynomial columns.
    pressure_gradient : complex ndarray, shape (N, 3, M)
        Full Cartesian gradients of the pressure columns.
    """

    points = _points(points)
    r0, z0 = _pair(center, "center")
    dr, dz = _pair(widths, "widths")
    if r0 <= 0.0:
        raise ValueError("center radius r0 must be positive")
    if dr <= 0.0 or dz <= 0.0:
        raise ValueError("widths must be positive")
    if dr >= r0:
        raise ValueError("radial half-width dr must satisfy dr < r0")
    mode = _integer(mode, "mode")
    degree = _integer(degree, "degree")
    if degree < 0:
        raise ValueError("degree must be nonnegative")
    kr, kz = _pair(carrier, "carrier")

    n = degree + 1
    count = n * n
    x, y, z = points.T
    radius = np.hypot(x, y)
    theta = np.arctan2(y, x)
    safe_radius = np.where(radius > 0.0, radius, 1.0)
    xi = (radius - r0) / dr
    eta = (z - z0) / dz

    br, brr = _bump(radius, r0, dr)
    bz, bzz = _bump(z, z0, dz)
    envelope = br * bz
    envelope_r = brr * bz
    envelope_z = br * bzz

    xi_powers = [np.ones_like(xi)]
    eta_powers = [np.ones_like(eta)]
    for _ in range(degree):
        xi_powers.append(xi_powers[-1] * xi)
        eta_powers.append(eta_powers[-1] * eta)

    phase = np.exp(1j * (mode * theta + kr * (radius - r0) + kz * (z - z0)))
    velocity_cyl = np.zeros((len(points), 3, 3 * count), dtype=complex)
    pressure = np.zeros((len(points), count), dtype=complex)
    pressure_gradient_cyl = np.zeros((len(points), 3, count), dtype=complex)

    for a in range(n):
        for b in range(n):
            index = a * n + b
            polynomial = xi_powers[a] * eta_powers[b]
            q0 = envelope * polynomial
            q0_r = envelope_r * polynomial
            q0_z = envelope_z * polynomial
            if a:
                q0_r = q0_r + (a / dr) * envelope * xi_powers[a - 1] * eta_powers[b]
            if b:
                q0_z = q0_z + (b / dz) * envelope * xi_powers[a] * eta_powers[b - 1]

            q = q0 * phase
            qr = (q0_r + 1j * kr * q0) * phase
            qz = (q0_z + 1j * kz * q0) * phase
            angular = 1j * mode * q / safe_radius

            pressure[:, index] = q
            pressure_gradient_cyl[:, 0, index] = qr
            pressure_gradient_cyl[:, 1, index] = angular
            pressure_gradient_cyl[:, 2, index] = qz

            radial = index
            azimuthal = count + index
            axial = 2 * count + index

            # The first three blocks are -curl(A_r), -curl(A_theta),
            # and -curl(A_z), including the cylindrical 1/r connection terms.
            velocity_cyl[:, 1, radial] = -qz
            velocity_cyl[:, 2, radial] = angular
            velocity_cyl[:, 0, azimuthal] = qz
            velocity_cyl[:, 2, azimuthal] = -(qr + q / safe_radius)
            velocity_cyl[:, 0, axial] = -angular
            velocity_cyl[:, 1, axial] = qr

    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    velocity = np.empty_like(velocity_cyl)
    velocity[:, 0, :] = cos_theta[:, None] * velocity_cyl[:, 0, :] - sin_theta[:, None] * velocity_cyl[:, 1, :]
    velocity[:, 1, :] = sin_theta[:, None] * velocity_cyl[:, 0, :] + cos_theta[:, None] * velocity_cyl[:, 1, :]
    velocity[:, 2, :] = velocity_cyl[:, 2, :]

    pressure_gradient = np.empty_like(pressure_gradient_cyl)
    pressure_gradient[:, 0, :] = cos_theta[:, None] * pressure_gradient_cyl[:, 0, :] - sin_theta[:, None] * pressure_gradient_cyl[:, 1, :]
    pressure_gradient[:, 1, :] = sin_theta[:, None] * pressure_gradient_cyl[:, 0, :] + cos_theta[:, None] * pressure_gradient_cyl[:, 1, :]
    pressure_gradient[:, 2, :] = pressure_gradient_cyl[:, 2, :]

    # This also makes the zero-support and axis behavior explicit if a caller
    # supplies a point at a boundary where floating point cancellation leaves
    # a tiny intermediate value.
    valid = (
        (radius > 0.0)
        & (radius > r0 - dr)
        & (radius < r0 + dr)
        & (z > z0 - dz)
        & (z < z0 + dz)
    )
    velocity[~valid, :, :] = 0.0
    pressure[~valid, :] = 0.0
    pressure_gradient[~valid, :, :] = 0.0
    return velocity, pressure, pressure_gradient


def _finite_difference_checks():
    """Run inexpensive independent checks and return JSON-serializable data."""

    center = (0.004, -0.001)
    widths = (0.0008, 0.0012)
    mode = 2
    degree = 2
    carrier = (130.0, -90.0)
    radius = center[0] + widths[0] * np.array([-0.45, -0.15, 0.20, 0.50])
    angle = np.array([0.20, 0.80, 1.40, 2.20])
    z = center[1] + widths[1] * np.array([-0.40, 0.10, 0.35, 0.60])
    points = np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))
    # A small centered step keeps the finite-difference truncation below the
    # compact bump's fourth-derivative scale while avoiding roundoff from the
    # Cartesian rotation.
    h = 2.0e-9

    velocity, pressure, pressure_gradient = basis_data(
        points, center, widths, mode, degree, carrier
    )
    divergence = np.zeros((len(points), velocity.shape[2]), dtype=complex)
    pressure_fd = np.zeros_like(pressure_gradient)
    derivative_scale = 0.0
    for coordinate in range(3):
        step = np.zeros(3)
        step[coordinate] = h
        vp = basis_data(points + step, center, widths, mode, degree, carrier)[0]
        vm = basis_data(points - step, center, widths, mode, degree, carrier)[0]
        derivative = (vp - vm) / (2.0 * h)
        divergence += derivative[:, coordinate, :]
        derivative_scale = max(derivative_scale, float(np.max(np.abs(derivative))))
        pp = basis_data(points + step, center, widths, mode, degree, carrier)[1]
        pm = basis_data(points - step, center, widths, mode, degree, carrier)[1]
        pressure_fd[:, coordinate, :] = (pp - pm) / (2.0 * h)

    collocation_velocity_error = 0.0
    collocation_gradient_error = 0.0
    try:
        from .curl_wave_patch_collocation import basis as collocation_basis
    except ImportError:
        # The neighboring collocation module uses historical absolute imports
        # and therefore needs its directory on sys.path when this file is run
        # as ``python -m experiments.root_st073.supported_fourier_basis``.
        module_dir = str(Path(__file__).resolve().parent)
        if module_dir not in sys.path:
            sys.path.insert(0, module_dir)
        from curl_wave_patch_collocation import basis as collocation_basis
    wave = SimpleNamespace(
        radius=center[0], zcenter=center[1], radial_halfwidth=widths[0], axial_halfwidth=widths[1]
    )
    collocation_mode = {"m": mode, "normal": (carrier[0], 0.0, carrier[1])}
    theta0 = np.column_stack((radius, np.zeros_like(radius), z))
    velocity0, _, gradient0 = basis_data(theta0, center, widths, mode, degree, carrier)
    for row, (r_value, z_value) in enumerate(zip(radius, z)):
        target = collocation_basis(wave, collocation_mode, r_value, z_value)
        collocation_velocity_error = max(
            collocation_velocity_error,
            float(np.max(np.abs(velocity0[row] - target[:, :27]))),
        )
        collocation_gradient_error = max(
            collocation_gradient_error,
            float(np.max(np.abs(gradient0[row] - target[:, 27:]))),
        )

    outside = np.array(
        [
            [center[0] - widths[0] - 1.0e-5, 0.0, center[1]],
            [center[0] + widths[0] + 1.0e-5, 0.0, center[1]],
            [center[0], 0.0, center[1] - widths[1] - 1.0e-5],
            [center[0], 0.0, center[1] + widths[1] + 1.0e-5],
            [0.0, 0.0, center[1]],
        ]
    )
    outside_velocity, outside_pressure, outside_gradient = basis_data(
        outside, center, widths, mode, degree, carrier
    )

    return {
        "layout": basis_layout(degree),
        "degree2_shapes": {
            "velocity": list(velocity.shape),
            "pressure": list(pressure.shape),
            "pressure_gradient": list(pressure_gradient.shape),
        },
        "generic_point_count": int(len(points)),
        "finite_difference_step": h,
        "divergence_fd_max": float(np.max(np.abs(divergence))),
        "divergence_fd_relative_to_max_derivative": float(
            np.max(np.abs(divergence)) / max(derivative_scale, np.finfo(float).tiny)
        ),
        "velocity_derivative_max": derivative_scale,
        "pressure_gradient_fd_max": float(np.max(np.abs(pressure_gradient - pressure_fd))),
        "collocation_velocity_max_abs_diff": collocation_velocity_error,
        "collocation_pressure_gradient_max_abs_diff": collocation_gradient_error,
        "outside_support_max": float(
            max(
                np.max(np.abs(outside_velocity)),
                np.max(np.abs(outside_pressure)),
                np.max(np.abs(outside_gradient)),
            )
        ),
        "checks_passed": bool(
            np.max(np.abs(divergence)) < 2.0e-3
            and np.max(np.abs(pressure_gradient - pressure_fd)) < 2.0e-6
            and collocation_velocity_error < 2.0e-12
            and collocation_gradient_error < 2.0e-12
            and max(
                np.max(np.abs(outside_velocity)),
                np.max(np.abs(outside_pressure)),
                np.max(np.abs(outside_gradient)),
            )
            == 0.0
        ),
        "accepted": False,
        "pde_validated": False,
        "scope": "manufactured basis checks only; no physical NS acceptance",
    }


if __name__ == "__main__":
    result = _finite_difference_checks()
    output = Path(__file__).with_name("supported_fourier_basis_check.json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
