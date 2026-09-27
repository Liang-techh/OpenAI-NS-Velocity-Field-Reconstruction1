"""Analytic spatial jets for the compact supported Fourier potential basis.

This module is a drop-in diagnostic companion to ``fourier_patch_evolution``.
It preserves that module's ``basis_jets`` return convention while replacing
finite-difference spatial derivatives with exact polynomial/carrier
derivatives in cylindrical coordinates.  The existing basis and evolution
modules remain unchanged.
"""

from __future__ import annotations

import json
from math import comb
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from supported_fourier_basis import basis_data


def _points(points):
    result = np.asarray(points, dtype=float)
    if result.ndim == 1 and result.shape == (3,):
        result = result[None, :]
    if result.ndim != 2 or result.shape[1] != 3:
        raise ValueError("points must have shape (N, 3)")
    if not np.all(np.isfinite(result)):
        raise ValueError("points must be finite")
    return result


def _pair(value, name):
    result = np.asarray(tuple(value), dtype=float)
    if result.shape != (2,) or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain two finite values")
    return float(result[0]), float(result[1])


def _factor_derivatives(value, center, halfwidth, power, order=3):
    """Derivatives of ``xi**power * (1 - xi**2)**5`` in ``value``."""
    xi = (value - center) / halfwidth
    inside = (xi > -1.0) & (xi < 1.0)
    values = []
    for derivative_order in range(order + 1):
        value_derivative = np.zeros_like(xi, dtype=float)
        if np.any(inside):
            x = xi[inside]
            derivative = np.zeros_like(x)
            # Leibniz rule for xi**power * (1-xi)**5 * (1+xi)**5.
            # Keeping the three factors explicit avoids cancellation near a
            # compact-support edge in the expanded degree-10 polynomial.
            for first in range(derivative_order + 1):
                for second in range(derivative_order - first + 1):
                    third = derivative_order - first - second
                    if first > power or second > 5 or third > 5:
                        continue
                    falling_power = 1
                    for step in range(first):
                        falling_power *= power - step
                    falling_left = 1
                    for step in range(second):
                        falling_left *= 5 - step
                    falling_right = 1
                    for step in range(third):
                        falling_right *= 5 - step
                    derivative += (
                        comb(derivative_order, first)
                        * comb(derivative_order - first, second)
                        * falling_power * falling_left * falling_right
                        * x**(power - first)
                        * (-1.0)**second * (1.0 - x)**(5 - second)
                        * (1.0 + x)**(5 - third)
                    )
            value_derivative[inside] = derivative / halfwidth**derivative_order
        values.append(value_derivative)
    return values


def _q_derivatives(points, center, widths, mode, degree, carrier):
    """Return q derivatives ``d_r**i d_z**j q`` for 0 <= i,j <= 3."""
    points = _points(points)
    r0, z0 = _pair(center, "center")
    dr, dz = _pair(widths, "widths")
    if r0 <= 0.0 or dr <= 0.0 or dz <= 0.0 or dr >= r0:
        raise ValueError("need r0 > dr > 0 and dz > 0")
    mode = int(mode)
    degree = int(degree)
    if degree < 0:
        raise ValueError("degree must be nonnegative")
    kr, kz = _pair(carrier, "carrier")

    x, y, z = points.T
    radius = np.hypot(x, y)
    safe_radius = np.where(radius > 0.0, radius, 1.0)
    theta = np.arctan2(y, x)
    radial_inside = (radius > r0 - dr) & (radius < r0 + dr)
    axial_inside = (z > z0 - dz) & (z < z0 + dz)
    valid = (radius > 0.0) & radial_inside & axial_inside

    n = degree + 1
    count = n * n
    radial_factors = [
        _factor_derivatives(radius, r0, dr, power, order=3)
        for power in range(n)
    ]
    axial_factors = [
        _factor_derivatives(z, z0, dz, power, order=3)
        for power in range(n)
    ]
    radial_carrier = np.zeros((4, len(points)), dtype=complex)
    axial_carrier = np.zeros((4, len(points)), dtype=complex)
    for derivative_order in range(4):
        for lower in range(derivative_order + 1):
            radial_carrier[derivative_order] += (
                comb(derivative_order, lower)
                * radial_factors[0][lower]
                * (1j * kr)**(derivative_order - lower)
            )
            axial_carrier[derivative_order] += (
                comb(derivative_order, lower)
                * axial_factors[0][lower]
                * (1j * kz)**(derivative_order - lower)
            )

    phase = np.exp(1j * (mode * theta + kr * (radius - r0)
                         + kz * (z - z0)))
    qder = np.zeros((4, 4, len(points), count), dtype=complex)
    for a in range(n):
        for b in range(n):
            index = a * n + b
            radial_carrier = np.zeros((4, len(points)), dtype=complex)
            axial_carrier = np.zeros((4, len(points)), dtype=complex)
            for derivative_order in range(4):
                for lower in range(derivative_order + 1):
                    radial_carrier[derivative_order] += (
                        comb(derivative_order, lower)
                        * radial_factors[a][lower]
                        * (1j * kr)**(derivative_order - lower)
                    )
                    axial_carrier[derivative_order] += (
                        comb(derivative_order, lower)
                        * axial_factors[b][lower]
                        * (1j * kz)**(derivative_order - lower)
                    )
            for radial_order in range(4):
                for axial_order in range(4):
                    qder[radial_order, axial_order, :, index] = (
                        radial_carrier[radial_order]
                        * axial_carrier[axial_order] * phase
                    )

    # The factors already vanish outside their individual supports.  Applying
    # the complete validity mask makes axis, boundary, and exterior behavior
    # identical to supported_fourier_basis.basis_data.
    qder[:, :, ~valid, :] = 0.0
    return qder, radius, safe_radius, valid


def _rotate_cylindrical(values, cosine, sine):
    return np.einsum("nia,naj->nij", _rotation(cosine, sine), values)


def _analytic_data(points, center, widths, mode, degree, carrier):
    points = _points(points)
    qder, radius, safe_radius, valid = _q_derivatives(
        points, center, widths, mode, degree, carrier)
    npoints = len(points)
    count = (int(degree) + 1)**2
    mode_factor = 1j * int(mode)
    radial = qder[0, 0]
    cosine = np.divide(points[:, 0], safe_radius,
                       out=np.zeros(npoints), where=safe_radius != 0.0)
    sine = np.divide(points[:, 1], safe_radius,
                     out=np.zeros(npoints), where=safe_radius != 0.0)

    velocity_cyl = np.zeros((npoints, 3, 3 * count), dtype=complex)
    gradient_cyl = np.zeros((npoints, 3, 3, 3 * count), dtype=complex)
    laplacian_cyl = np.zeros_like(velocity_cyl)
    pressure_gradient_cyl = np.zeros((npoints, 3, count), dtype=complex)

    def scalar_laplacian(value, radial_derivative, radial_second,
                         axial_second):
        return (radial_second + radial_derivative / safe_radius
                - (int(mode)**2) * value / safe_radius**2 + axial_second)

    for index in range(count):
        q = qder[0, 0, :, index]
        qr = qder[1, 0, :, index]
        qz = qder[0, 1, :, index]
        qrr = qder[2, 0, :, index]
        qrz = qder[1, 1, :, index]
        qzz = qder[0, 2, :, index]
        qrrr = qder[3, 0, :, index]
        qrrz = qder[2, 1, :, index]
        qrzz = qder[1, 2, :, index]
        qzzz = qder[0, 3, :, index]
        r = safe_radius

        pressure_gradient_cyl[:, 0, index] = qr
        pressure_gradient_cyl[:, 1, index] = mode_factor * q / r
        pressure_gradient_cyl[:, 2, index] = qz

        for component, block in enumerate((0, count, 2 * count)):
            if component == 0:  # A_r = q.
                vr, vt, vz = np.zeros_like(q), -qz, mode_factor * q / r
                vr_r, vt_r, vz_r = np.zeros_like(q), -qrz, mode_factor * (qr / r - q / r**2)
                vr_z, vt_z, vz_z = np.zeros_like(q), -qzz, mode_factor * qz / r
                vr_rr, vt_rr, vz_rr = np.zeros_like(q), -qrrz, mode_factor * (qrr / r - 2 * qr / r**2 + 2 * q / r**3)
                vr_zz, vt_zz, vz_zz = np.zeros_like(q), -qzzz, mode_factor * qzz / r
            elif component == 1:  # A_theta = q.
                vr, vt, vz = qz, np.zeros_like(q), -(qr + q / r)
                vr_r, vt_r, vz_r = qrz, np.zeros_like(q), -(qrr + qr / r - q / r**2)
                vr_z, vt_z, vz_z = qzz, np.zeros_like(q), -(qrz + qz / r)
                vr_rr, vt_rr, vz_rr = qrrz, np.zeros_like(q), -(qrrr + qrr / r - 2 * qr / r**2 + 2 * q / r**3)
                vr_zz, vt_zz, vz_zz = qzzz, np.zeros_like(q), -(qrzz + qzz / r)
            else:  # A_z = q.
                vr, vt, vz = -mode_factor * q / r, qr, np.zeros_like(q)
                vr_r, vt_r, vz_r = -mode_factor * (qr / r - q / r**2), qrr, np.zeros_like(q)
                vr_z, vt_z, vz_z = -mode_factor * qz / r, qrz, np.zeros_like(q)
                vr_rr, vt_rr, vz_rr = -mode_factor * (qrr / r - 2 * qr / r**2 + 2 * q / r**3), qrrr, np.zeros_like(q)
                vr_zz, vt_zz, vz_zz = -mode_factor * qzz / r, qrzz, np.zeros_like(q)

            velocity_cyl[:, 0, block + index] = vr
            velocity_cyl[:, 1, block + index] = vt
            velocity_cyl[:, 2, block + index] = vz

            gradient_cyl[:, 0, 0, block + index] = vr_r
            gradient_cyl[:, 0, 1, block + index] = (mode_factor * vr - vt) / r
            gradient_cyl[:, 0, 2, block + index] = vr_z
            gradient_cyl[:, 1, 0, block + index] = vt_r
            gradient_cyl[:, 1, 1, block + index] = (mode_factor * vt + vr) / r
            gradient_cyl[:, 1, 2, block + index] = vt_z
            gradient_cyl[:, 2, 0, block + index] = vz_r
            gradient_cyl[:, 2, 1, block + index] = mode_factor * vz / r
            gradient_cyl[:, 2, 2, block + index] = vz_z

            scalar_r = scalar_laplacian(vr, vr_r, vr_rr, vr_zz)
            scalar_t = scalar_laplacian(vt, vt_r, vt_rr, vt_zz)
            scalar_z = scalar_laplacian(vz, vz_r, vz_rr, vz_zz)
            laplacian_cyl[:, 0, block + index] = (
                scalar_r - vr / r**2 - 2 * mode_factor * vt / r**2)
            laplacian_cyl[:, 1, block + index] = (
                scalar_t - vt / r**2 + 2 * mode_factor * vr / r**2)
            laplacian_cyl[:, 2, block + index] = scalar_z

    velocity = _rotate_cylindrical(velocity_cyl, cosine, sine)
    gradient = np.einsum("nia,nabq,njb->nijq", _rotation(cosine, sine),
                         gradient_cyl, _rotation(cosine, sine))
    laplacian = _rotate_cylindrical(laplacian_cyl, cosine, sine)
    pressure_gradient = _rotate_cylindrical(
        pressure_gradient_cyl, cosine, sine)
    velocity[~valid] = 0.0
    gradient[~valid] = 0.0
    laplacian[~valid] = 0.0
    pressure_gradient[~valid] = 0.0
    return velocity, gradient, laplacian, radial, pressure_gradient


def _rotation(cosine, sine):
    result = np.zeros((len(cosine), 3, 3), dtype=float)
    result[:, 0, 0] = cosine
    result[:, 0, 1] = -sine
    result[:, 1, 0] = sine
    result[:, 1, 1] = cosine
    result[:, 2, 2] = 1.0
    return result


def basis_jets(points, center, widths, mode, degree, carrier, nu, h=None):
    """Return analytic ``(V, J, -nu*lapV, P, gradP)`` basis jets.

    The optional ``h`` parameter is accepted for drop-in compatibility with
    ``fourier_patch_evolution.basis_jets`` and is intentionally ignored.
    """
    velocity, gradient, laplacian, pressure, pressure_gradient = _analytic_data(
        points, center, widths, mode, degree, carrier)
    return velocity, gradient, -float(nu) * laplacian, pressure, pressure_gradient


def _finite_difference_jets(points, center, widths, mode, degree, carrier,
                            nu, h):
    velocity, pressure, pressure_gradient = basis_data(
        points, center, widths, mode, degree, carrier)
    gradient = np.empty((len(points), 3, 3, velocity.shape[-1]), complex)
    laplacian = np.zeros_like(velocity)
    for axis in range(3):
        step = np.eye(3)[axis] * h
        vm2, vm, vp, vp2 = [basis_data(
            points + multiple * step, center, widths, mode, degree, carrier)[0]
            for multiple in (-2, -1, 1, 2)]
        gradient[:, :, axis] = (vm2 - 8 * vm + 8 * vp - vp2) / (12 * h)
        laplacian += (-vp2 + 16 * vp - 30 * velocity + 16 * vm - vm2) / (
            12 * h * h)
    return velocity, gradient, -float(nu) * laplacian, pressure, pressure_gradient


def _max_abs(value):
    return float(np.max(np.abs(value))) if np.size(value) else 0.0


def _reference_values(points, center, widths, mode, degree, carrier):
    analytic = basis_jets(points, center, widths, mode, degree, carrier, 0.01)
    exact = basis_data(points, center, widths, mode, degree, carrier)
    return {
        "velocity_max_abs_error": _max_abs(analytic[0] - exact[0]),
        "pressure_max_abs_error": _max_abs(analytic[3] - exact[1]),
        "pressure_gradient_max_abs_error": _max_abs(analytic[4] - exact[2]),
        "velocity_scale": _max_abs(exact[0]),
        "pressure_scale": _max_abs(exact[1]),
        "pressure_gradient_scale": _max_abs(exact[2]),
    }


def _fd_comparison(points, center, widths, mode, degree, carrier, nu):
    analytic = basis_jets(points, center, widths, mode, degree, carrier, nu)
    h_values = tuple(float(min(widths) * factor)
                     for factor in (2.0e-4, 1.0e-4, 5.0e-5))
    rows = []
    for h in h_values:
        finite_difference = _finite_difference_jets(
            points, center, widths, mode, degree, carrier, nu, h)
        rows.append({
            "h": h,
            "gradient_max_abs_error": _max_abs(analytic[1] - finite_difference[1]),
            "viscous_laplacian_max_abs_error": _max_abs(analytic[2] - finite_difference[2]),
            "gradient_scale": _max_abs(analytic[1]),
            "viscous_laplacian_scale": _max_abs(analytic[2]),
            "gradient_relative_error": _max_abs(analytic[1] - finite_difference[1]) / max(_max_abs(analytic[1]), np.finfo(float).tiny),
            "viscous_laplacian_relative_error": _max_abs(analytic[2] - finite_difference[2]) / max(_max_abs(analytic[2]), np.finfo(float).tiny),
        })
    return rows


def _saved_points(report):
    center = np.asarray(report["center"], dtype=float)
    widths = np.asarray(report["widths"], dtype=float)
    nodes, _ = leggauss(7)
    axis = 0.95 * nodes
    train = np.asarray([[center[0] + x * widths[0], 0.0,
                         center[1] + z * widths[1]]
                        for x in axis for z in axis])
    held_axis = np.linspace(-0.8, 0.8, 6)
    offsets = ([(x, z) for x in held_axis for z in held_axis]
               + [(x, z) for x in (-0.97, 0.97) for z in (-0.97, 0.97)])
    held = np.asarray([[center[0] + x * widths[0], 0.0,
                        center[1] + z * widths[1]] for x, z in offsets])
    # Use interior saved points for FD so a finite-difference stencil does not
    # cross a compact support boundary. Boundary behavior is checked separately.
    interior_held = held[np.all(np.abs((held[:, [0, 2]]
                                        - np.array([center[0], center[1]]))
                                       / widths) <= 0.6, axis=1)]
    points = np.vstack((train[:8], interior_held[:8]))
    return train, held, points


def run():
    output = Path(__file__).with_name("supported_fourier_analytic_jets.json")
    saved = json.loads((Path(__file__).with_name(
        "fourier_patch_evolution.json")).read_text(encoding="utf-8"))
    nu = 0.01

    normal_cases = [
        dict(name="normal_mode2_carrier", center=(0.004, -0.001),
             widths=(0.0008, 0.0012), mode=2, degree=2,
             carrier=(130.0, -90.0)),
        dict(name="simple_mode0", center=(0.005, 0.002), widths=(0.001, 0.001),
             mode=0, degree=2, carrier=(0.0, 0.0)),
    ]
    normal_results = []
    for case in normal_cases:
        r0, z0 = case["center"]
        dr, dz = case["widths"]
        offsets = np.asarray([[-0.45, -0.4], [-0.15, 0.1],
                              [0.2, 0.35], [0.5, 0.6]])
        points = np.column_stack((
            (r0 + dr * offsets[:, 0]), np.zeros(len(offsets)),
            z0 + dz * offsets[:, 1]))
        result = dict(case=case, point_count=len(points),
                      reference=_reference_values(
                          points, case["center"], case["widths"],
                          case["mode"], case["degree"], case["carrier"]),
                      finite_difference=_fd_comparison(
                          points, case["center"], case["widths"], case["mode"],
                          case["degree"], case["carrier"], nu))
        analytic = basis_jets(points, case["center"], case["widths"],
                              case["mode"], case["degree"], case["carrier"], nu)
        result["divergence_max_abs"] = _max_abs(
            np.einsum("niiq->nq", analytic[1]))
        result["divergence_scale"] = _max_abs(analytic[1])
        normal_results.append(result)

    train, held, saved_points = _saved_points(saved)
    saved_center = tuple(saved["center"])
    saved_widths = tuple(saved["widths"])
    saved_degree = int(saved["state_data"]["degree"])
    saved_carriers = {int(key): tuple(value)
                      for key, value in saved["carriers"].items()}
    saved_fd = []
    saved_reference = []
    saved_divergence = []
    for mode in sorted(saved_carriers):
        carrier = saved_carriers[mode]
        saved_reference.append(dict(
            mode=mode,
            values=_reference_values(saved_points, saved_center, saved_widths,
                                     mode, saved_degree, carrier)))
        saved_fd.append(dict(
            mode=mode,
            carrier=list(carrier),
            finite_difference=_fd_comparison(
                saved_points, saved_center, saved_widths, mode, saved_degree,
                carrier, nu)))
        analytic = basis_jets(saved_points, saved_center, saved_widths, mode,
                              saved_degree, carrier, nu)
        saved_divergence.append(dict(
            mode=mode,
            divergence_max_abs=_max_abs(np.einsum("niiq->nq", analytic[1])),
            gradient_scale=_max_abs(analytic[1]),
            viscous_laplacian_scale=_max_abs(analytic[2])))

    r0, z0 = saved_center
    dr, dz = saved_widths
    representative = np.asarray([
        [r0, 0.0, z0],
        [r0 - dr * (1.0 - 1.0e-6), 0.0, z0],
        [r0 + dr * (1.0 - 1.0e-10), 0.0, z0],
        [r0, 0.0, z0 - dz * (1.0 - 1.0e-6)],
        [r0, 0.0, z0 + dz * (1.0 - 1.0e-10)],
        [r0 - dr, 0.0, z0],
        [r0 + dr, 0.0, z0],
        [r0, 0.0, z0 - dz],
        [r0, 0.0, z0 + dz],
        [r0 - dr - 1.0e-5, 0.0, z0],
        [r0 + dr + 1.0e-5, 0.0, z0],
        [r0, 0.0, z0 - dz - 1.0e-5],
        [r0, 0.0, z0 + dz + 1.0e-5],
        [0.0, 0.0, z0],
    ])
    support_checks = []
    for name, point in zip(("interior", "radial_inner_near_boundary",
                            "radial_outer_near_boundary",
                            "axial_lower_near_boundary",
                            "axial_upper_near_boundary",
                            "radial_inner_boundary", "radial_outer_boundary",
                            "axial_lower_boundary", "axial_upper_boundary",
                            "radial_exterior_low", "radial_exterior_high",
                            "axial_exterior_low", "axial_exterior_high", "axis"), representative):
        values = _reference_values(
            point[None, :], saved_center, saved_widths, 4, saved_degree,
            saved_carriers[4])
        support_checks.append(dict(name=name, point=point.tolist(), **values))

    all_divergence = [row["divergence_max_abs"] for row in normal_results]
    all_divergence.extend(row["divergence_max_abs"] for row in saved_divergence)
    report = {
        "source": "Analytic cylindrical polynomial/carrier derivatives for the supported Fourier potential basis",
        "function": "basis_jets(points, center, widths, mode, degree, carrier, nu, h=None)",
        "return_contract": ["V", "J", "-nu*lapV", "P", "gradP"],
        "shapes_degree2": {
            "velocity": ["N", 3, 27], "gradient": ["N", 3, 3, 27],
            "viscous_laplacian": ["N", 3, 27], "pressure": ["N", 9],
            "pressure_gradient": ["N", 3, 9],
        },
        "derivative_method": (
            "q=xi^a eta^b (1-xi^2)^5 (1-eta^2)^5 exp(i(m theta+kr(r-r0)+kz(z-z0))); "
            "polynomial derivatives through total order 3, cylindrical connection terms, "
            "and Cartesian rotation R J_cyl R^T"),
        "normal_cases": normal_results,
        "saved_candidate": {
            "source": "fourier_patch_evolution.json",
            "center": list(saved_center), "widths": list(saved_widths),
            "degree": saved_degree, "train_point_count": len(train),
            "held_point_count": len(held),
            "fd_point_count": len(saved_points),
            "reference_values": saved_reference,
            "finite_difference": saved_fd,
            "divergence": saved_divergence,
        },
        "support_and_axis_reference_checks": support_checks,
        "divergence_summary": {
            "max_abs_over_cases": float(max(all_divergence)),
            "cases_are_sampled": True,
        },
        "scope": (
            "Analytic basis-jet diagnostic only. Finite-difference comparisons are local "
            "and step-dependent; support-boundary behavior is checked by exact basis replay. "
            "This does not change the default backend, certify global accuracy, validate NS, "
            "or establish scale recursion."),
        "accepted": False,
        "scale_recursion_established": False,
        "pde_validated": False,
    }
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "normal_cases": len(normal_results),
        "saved_modes": len(saved_carriers),
        "max_divergence": report["divergence_summary"]["max_abs_over_cases"],
        "saved_reference_max": max(
            row["values"]["velocity_max_abs_error"]
            for row in saved_reference),
    }, indent=2), flush=True)
    return report


if __name__ == "__main__":
    run()
