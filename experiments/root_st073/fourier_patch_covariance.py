"""Cheap covariance and oscillatory-energy diagnostic for the saved BDF path.

This uses only the saved Hermite states and ``supported_fourier_basis``.  It
does not construct the evolving mean field or evaluate background jets.  The
reported energy is an angular-mode proxy, and the covariance is sampled at
the fixed patch center; neither is a total-energy or evolving stress test.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import CubicHermiteSpline


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from supported_fourier_basis import basis_data  # noqa: E402


def run(
    input_path=ROOT / "fourier_patch_implicit.json",
    reference_path=ROOT / "feedback_moving_wave.json",
    output_path=ROOT / "fourier_patch_covariance.json",
):
    """Recompute the saved covariance and energy-proxy report."""

    report = json.loads(Path(input_path).read_text(encoding="utf-8"))
    if not report.get("state_data"):
        print("state_data absent; no diagnostic written")
        return None
    source = json.loads(Path(reference_path).read_text(encoding="utf-8"))
    data = report["state_data"]
    decode = lambda a: np.asarray(a)[..., 0] + 1j * np.asarray(a)[..., 1]
    times = np.asarray(data["physical_times"], dtype=float)
    states = decode(data["states"])
    slopes = decode(data["slopes"])
    state_path = CubicHermiteSpline(times, states, slopes, extrapolate=False)
    center = tuple(float(v) for v in report["center"])
    widths = tuple(float(v) for v in report["widths"])
    degree = int(report["degree"])
    modes = [int(m) for m in report["modes"]]
    carriers = {
        int(m): np.asarray(v, dtype=float) for m, v in report["carriers"].items()
    }

    angles = np.arange(64, dtype=float) * 2.0 * np.pi / 64.0
    r0, z0 = center
    center_points = np.column_stack(
        (r0 * np.cos(angles), r0 * np.sin(angles), np.full(64, z0))
    )

    quad_nodes, quad_weights = leggauss(16)
    quad_points = np.array(
        [
            [r0 + xr * widths[0], 0.0, z0 + zr * widths[1]]
            for xr in quad_nodes
            for zr in quad_nodes
        ]
    )
    quad_weight = np.outer(quad_weights, quad_weights).ravel() * quad_points[:, 0]
    quad_weight /= np.sum(quad_weight)

    def mode_velocity(points, mode, state):
        velocity, _, _ = basis_data(
            points,
            center=center,
            widths=widths,
            mode=mode,
            degree=degree,
            carrier=carriers[mode],
        )
        return np.einsum("niq,q->ni", velocity, state[mode])

    def covariance_and_energy(state):
        oscillatory = np.zeros((len(center_points), 3), dtype=complex)
        energy_density = np.zeros(len(quad_points), dtype=float)
        mode_energy = {}
        for mode in modes:
            if mode == 0:
                continue
            center_velocity = mode_velocity(center_points, mode, state)
            oscillatory += center_velocity
            quadrature_velocity = mode_velocity(quad_points, mode, state)
            contribution = 0.5 * np.sum(np.abs(quadrature_velocity) ** 2, axis=1)
            energy_density += contribution
            mode_energy[str(mode)] = float(np.sum(quad_weight * contribution))

        ca, sa = np.cos(angles), np.sin(angles)
        cylindrical = np.column_stack(
            (
                ca * oscillatory[:, 0].real + sa * oscillatory[:, 1].real,
                -sa * oscillatory[:, 0].real + ca * oscillatory[:, 1].real,
                oscillatory[:, 2].real,
            )
        )
        covariance = np.mean(cylindrical[:, 0, None] * cylindrical[:, 1:], axis=0)
        return covariance, float(np.sum(quad_weight * energy_density)), mode_energy

    target = np.asarray(source["covariance_target"], dtype=float)
    initial_covariance, initial_energy, initial_mode_energy = covariance_and_energy(
        states[0]
    )
    rows = []
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        physical_time = float(times[0] + fraction * (times[-1] - times[0]))
        state = np.asarray(state_path(physical_time))
        covariance, energy, mode_energy = covariance_and_energy(state)
        rows.append(
            {
                "fraction": fraction,
                "physical_time": physical_time,
                "state_node_count": int(len(times)),
                "energy_proxy": energy,
                "energy_proxy_ratio_to_initial": float(energy / initial_energy)
                if initial_energy
                else None,
                "covariance_radial_times_theta_z": covariance.tolist(),
                "covariance_minus_initial_reference_target": (
                    covariance - target
                ).tolist(),
                "covariance_relative_error_to_initial_reference_target": float(
                    np.linalg.norm(covariance - target)
                    / max(np.linalg.norm(target), np.finfo(float).tiny)
                ),
                "mode_energy_proxy": mode_energy,
            }
        )

    result = {
        "accepted": False,
        "pde_validated": False,
        "scope": "Wave energy proxy and fixed-center oscillatory covariance diagnostic only; no total finite-energy or evolving corrected-mean stress matching.",
        "source_report": "fourier_patch_implicit.json",
        "reference_report": "feedback_moving_wave.json",
        "degree": degree,
        "modes_included_in_oscillatory_proxy": list(range(1, max(modes) + 1)),
        "center": list(center),
        "widths": list(widths),
        "angular_sample_count": int(len(angles)),
        "quadrature_order_per_axis": 16,
        "quadrature_point_count": int(len(quad_points)),
        "quadrature_weight_sum": float(np.sum(quad_weight)),
        "state_data_node_count": int(len(times)),
        "physical_times": times.tolist(),
        "initial_reference_covariance_target": target.tolist(),
        "initial_covariance": initial_covariance.tolist(),
        "initial_covariance_minus_reference_target": (
            initial_covariance - target
        ).tolist(),
        "initial_energy_proxy": initial_energy,
        "initial_mode_energy_proxy": initial_mode_energy,
        "definition": {
            "covariance": "mean over 64 angles of v_osc,r * [v_osc,theta, v_osc,z], with v_osc = Re(sum_{m>=1} V_m c_m).",
            "energy_proxy": "normalized r-weighted 16x16 full-support Gauss average of 0.5 * sum_{m>=1} |V_m c_m|^2.",
            "trajectory": "Cubic Hermite interpolation of saved BDF states and slopes; no mean.fields/jets calls.",
        },
        "rows": rows,
    }
    output_path = Path(output_path)
    output_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(output_path),
                "initial_energy_proxy": initial_energy,
                "initial_covariance": initial_covariance.tolist(),
                "rows": rows,
            },
            indent=2,
        )
    )
    return result


if __name__ == "__main__":
    run()
