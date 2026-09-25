"""Finite moving-patch transport of projected Kelvin amplitudes.

This is a bounded physical-coordinate experiment between the registered times
``tau=0.0084`` and ``tau=0.00846``.  The two transported harmonic modes use
the existing projected Kelvin ODE at a 3 by 3 set of moving radial/axial
nodes.  Their potentials are interpolated over the patch and an analytic
cylindrical curl is evaluated, so the added velocity is divergence-free by
construction.  A scalar pressure proxy is obtained from the same local
pressure projection used by the Kelvin ODE.

The calculation is deliberately smaller than a supported Section 7--9
construction: the phase is common to the interpolated patch mode, the
projection-derived pressure is not a separately solved pressure PDE, and no
global moment or energy claim is made.  Full Cartesian finite-difference
momentum is nevertheless evaluated on disjoint spatial holdouts at two times.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from curl_wave_prototype import LocalizedCurlWave, WavePerturbedField, bump
from joint_collar_fit import kinematics
from kelvin_amplitude_transport import derivative as source_kelvin_derivative
from local_poloidal_basis_screen import load_robust_candidate
from transported_phase_screen import coefficients


ROOT = Path(__file__).resolve().parent
COMPACT = ROOT / "compact_potential"
SOURCE_NAME = "local_poloidal_10pct_source.json"
TAU0 = 0.0084
TAU1 = 0.00846
RADIAL_HALF_WIDTH = 0.00275
AXIAL_HALF_WIDTH = 0.000075
TIME_HALF_WIDTH = 0.00015
PATCH_NODES = np.array([-0.6, 0.0, 0.6], dtype=float)
HOLDOUT_NODES = np.array([-0.45, 0.45], dtype=float)


def lagrange_values(value: float, nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return quadratic Lagrange values and derivatives at ``value``."""
    values = np.ones(len(nodes), dtype=float)
    derivatives = np.zeros(len(nodes), dtype=float)
    for i, node in enumerate(nodes):
        for j, other in enumerate(nodes):
            if i == j:
                continue
            values[i] *= (value - other) / (node - other)
        for j, other in enumerate(nodes):
            if i == j:
                continue
            term = 1.0 / (node - other)
            for k, third in enumerate(nodes):
                if k == i or k == j:
                    continue
                term *= (value - third) / (node - third)
            derivatives[i] += term
    return values, derivatives


def pressure_projection_from_coefficients(mode: int, state: np.ndarray, values: np.ndarray) -> float:
    """Compute the projected-pressure scalar from cached physical jets."""
    r, z, kr, kz, _ = state[:5]
    amplitude = np.asarray(state[5:8], dtype=float)
    ur, uz, swirl, ur_r, ur_z, uz_r, uz_z, swirl_r, swirl_z = values
    jacobian = np.array(
        [
            [ur_r, -swirl, ur_z],
            [swirl + r * swirl_r, ur / r, r * swirl_z],
            [uz_r, 0.0, uz_z],
        ],
        dtype=float,
    )
    normal = np.array([kr, mode / r, kz], dtype=float)
    return float(2.0 * np.dot(normal, jacobian @ amplitude) / np.dot(normal, normal))


class TransportedKelvinPatch:
    """Exact-curl wave whose potential coefficients evolve on a 3 by 3 patch."""

    def __init__(
        self,
        source: dict,
        base,
        radial_halfwidth: float = RADIAL_HALF_WIDTH,
        axial_halfwidth: float = AXIAL_HALF_WIDTH,
        time_halfwidth: float = TIME_HALF_WIDTH,
    ):
        self.base = base
        self.nu = base.nu
        self.tau0 = float(source["tau"])
        self.radial_halfwidth = float(radial_halfwidth)
        self.axial_halfwidth = float(axial_halfwidth)
        self.time_halfwidth = float(time_halfwidth)
        self.prototype = LocalizedCurlWave(
            source,
            radial_halfwidth=self.radial_halfwidth,
            axial_halfwidth=self.axial_halfwidth,
            time_halfwidth=self.time_halfwidth,
        )
        self.waves = self.prototype.waves
        self.weights = self.prototype.weights
        self.patch_nodes = PATCH_NODES.copy()
        # The full finite-difference coefficient call is expensive for this
        # candidate.  Cache a 3 by 3 by 9 physical jet table and integrate
        # each Kelvin state against its local interpolated table.  This keeps
        # the experiment bounded while making the approximation explicit.
        # Nine cached times are enough for this bounded diagnostic and avoid
        # rerunning the expensive five-point physical jet at every RK stage.
        self.integration_times = np.linspace(
            self.tau0 - 0.00012, self.tau0 + 0.00012, 9
        )
        self._coefficient_table = np.zeros((9, len(self.integration_times), 9))
        self._trajectories: dict[tuple[int, int, int], object] = {}
        self._snapshot_cache: dict[tuple[int, float], dict] = {}
        self.source_rhs_error = 0.0
        self._build_trajectories()

    def _rhs(self, mode: int, state: np.ndarray, values: np.ndarray) -> np.ndarray:
        """Algebraic RHS copied from kelvin_amplitude_transport.derivative."""
        r, z, kr, kz, _, ar, at, az = state
        ur, uz, swirl, ur_r, ur_z, uz_r, uz_z, swirl_r, swirl_z = values
        jacobian = np.array(
            [
                [ur_r, -swirl, ur_z],
                [swirl + r * swirl_r, ur / r, r * swirl_z],
                [uz_r, 0.0, uz_z],
            ],
            dtype=float,
        )
        amplitude = np.array([ar, at, az], dtype=float)
        normal = np.array([kr, mode / r, kz], dtype=float)
        ja = jacobian @ amplitude
        rotation = swirl * np.array([-at, ar, 0.0], dtype=float)
        projection = -2.0 * normal * np.dot(normal, ja) / np.dot(normal, normal)
        amp_tau = rotation + ja + self.nu * np.dot(normal, normal) * amplitude + projection
        return np.r_[-ur, -uz,
                     ur_r * kr + uz_r * kz + mode * swirl_r,
                     ur_z * kr + uz_z * kz + mode * swirl_z,
                     mode * swirl, amp_tau]

    def _table_values(self, node_index: int, tau: float) -> np.ndarray:
        values = self._coefficient_table[node_index]
        return np.array(
            [np.interp(float(tau), self.integration_times, values[:, column])
             for column in range(values.shape[1])],
            dtype=float,
        )

    def _integrate_state(self, mode: int, initial: np.ndarray, node_index: int):
        """Fixed-step RK4 on the cached jet table, with dense linear output."""
        times = self.integration_times
        states = np.zeros((len(times), 8), dtype=float)
        center = len(times) // 2
        states[center] = initial

        def step(state: np.ndarray, time: float, h: float) -> np.ndarray:
            def rhs(at: np.ndarray, tt: float) -> np.ndarray:
                return self._rhs(mode, at, self._table_values(node_index, tt))
            k1 = rhs(state, time)
            k2 = rhs(state + 0.5 * h * k1, time + 0.5 * h)
            k3 = rhs(state + 0.5 * h * k2, time + 0.5 * h)
            k4 = rhs(state + h * k3, time + h)
            return state + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

        for index in range(center + 1, len(times)):
            states[index] = step(states[index - 1], times[index - 1],
                                 times[index] - times[index - 1])
        for index in range(center - 1, -1, -1):
            states[index] = step(states[index + 1], times[index + 1],
                                 times[index] - times[index + 1])

        def trajectory(tau: float) -> np.ndarray:
            return np.array(
                [np.interp(float(tau), times, states[:, column])
                 for column in range(states.shape[1])],
                dtype=float,
            )

        return trajectory

    def _build_trajectories(self) -> None:
        r0, _, z0 = self.prototype.radius, 0.0, self.prototype.zcenter
        for i, xi in enumerate(self.patch_nodes):
            for j, eta in enumerate(self.patch_nodes):
                node_index = 3 * i + j
                radius = r0 + self.radial_halfwidth * xi
                height = z0 + self.axial_halfwidth * eta
                for index, tau in enumerate(self.integration_times):
                    self._coefficient_table[node_index, index] = coefficients(
                        self.base, radius, height, float(tau)
                    )
        for mode_index, mode_data in enumerate(self.waves):
            m = int(mode_data["m"])
            kr, _, kz = np.asarray(mode_data["normal"], dtype=float)
            for i, xi in enumerate(self.patch_nodes):
                for j, eta in enumerate(self.patch_nodes):
                    node_index = 3 * i + j
                    radius = r0 + self.radial_halfwidth * xi
                    height = z0 + self.axial_halfwidth * eta
                    normal = np.array([kr, m / radius, kz], dtype=float)
                    amplitude = np.asarray(mode_data["amplitude"], dtype=float).copy()
                    amplitude -= normal * np.dot(normal, amplitude) / np.dot(normal, normal)
                    initial = np.r_[radius, height, kr, kz, 0.0, amplitude]
                    self._trajectories[(mode_index, i, j)] = self._integrate_state(
                        m, initial, node_index
                    )
                    if i == 1 and j == 1:
                        reference = source_kelvin_derivative(
                            self.base, m, self.tau0, initial
                        )
                        cached = self._rhs(
                            m, initial,
                            self._table_values(node_index, self.tau0),
                        )
                        self.source_rhs_error = max(
                            self.source_rhs_error,
                            float(np.max(np.abs(reference - cached))),
                        )

    def _snapshot(self, mode_index: int, tau: float) -> dict:
        key = (mode_index, float(tau))
        if key in self._snapshot_cache:
            return self._snapshot_cache[key]
        mode_data = self.waves[mode_index]
        m = int(mode_data["m"])
        center_trajectory = self._trajectories[(mode_index, 1, 1)]
        center_state = np.asarray(center_trajectory(float(tau)), dtype=float)
        rc, zc, kr, kz, phi0 = center_state[:5]
        center_normal = np.array([kr, m / rc, kz], dtype=float)
        potentials = np.zeros((3, 3, 3), dtype=float)
        pressures = np.zeros((3, 3), dtype=float)
        transverse_errors = []
        for i in range(3):
            for j in range(3):
                state = np.asarray(self._trajectories[(mode_index, i, j)](float(tau)), dtype=float)
                amplitude = state[5:8].copy()
                local_normal = np.array([state[2], m / state[0], state[3]], dtype=float)
                transverse_errors.append(
                    float(abs(np.dot(local_normal, amplitude))
                          / (np.linalg.norm(local_normal) * np.linalg.norm(amplitude)))
                    if np.linalg.norm(amplitude) else 0.0
                )
                amplitude -= center_normal * np.dot(center_normal, amplitude) / np.dot(center_normal, center_normal)
                potentials[i, j] = np.cross(center_normal, amplitude) / np.dot(center_normal, center_normal)
                pressures[i, j] = pressure_projection_from_coefficients(
                    m, state, self._table_values(3 * i + j, float(tau))
                )
        snapshot = {
            "center_state": center_state,
            "center_normal": center_normal,
            "potentials": potentials,
            "pressures": pressures,
            "max_transversality_error": float(max(transverse_errors)),
            "amplitude_norm_ratio": float(
                np.linalg.norm(center_state[5:8]) / np.linalg.norm(mode_data["amplitude"])
            ),
        }
        self._snapshot_cache[key] = snapshot
        return snapshot

    def fields(self, points, tau):
        pts = np.asarray(points, dtype=float)
        ts = np.broadcast_to(tau, (len(pts),))
        velocity = np.zeros_like(pts)
        pressure = np.zeros(len(pts), dtype=float)
        for index, (point, time) in enumerate(zip(pts, ts)):
            time = float(time)
            if not self.tau0 - self.time_halfwidth < time < self.tau0 + self.time_halfwidth:
                continue
            x, ycart, z = point
            radius = float(np.hypot(x, ycart))
            if radius == 0.0:
                continue
            theta = float(np.arctan2(ycart, x))
            total_cyl = np.zeros(3, dtype=float)
            for mode_index, (weight, mode_data) in enumerate(zip(self.weights, self.waves)):
                m = int(mode_data["m"])
                snapshot = self._snapshot(mode_index, time)
                state = snapshot["center_state"]
                rc, zc, kr, kz, phi0 = state[:5]
                xi = (radius - rc) / self.radial_halfwidth
                eta = (z - zc) / self.axial_halfwidth
                if abs(xi) >= 1.0 or abs(eta) >= 1.0:
                    continue
                radial_bump, radial_bump_r = bump(radius, rc, self.radial_halfwidth)
                axial_bump, axial_bump_z = bump(z, zc, self.axial_halfwidth)
                time_bump, _ = bump(time, self.tau0, self.time_halfwidth)
                envelope = radial_bump * axial_bump * time_bump
                envelope_r = radial_bump_r * axial_bump * time_bump
                envelope_z = radial_bump * axial_bump_z * time_bump
                Lr, dLr = lagrange_values(xi, self.patch_nodes)
                Lz, dLz = lagrange_values(eta, self.patch_nodes)
                potential = np.zeros(3, dtype=float)
                potential_r = np.zeros(3, dtype=float)
                potential_z = np.zeros(3, dtype=float)
                pressure_amp = 0.0
                for i in range(3):
                    for j in range(3):
                        basis = Lr[i] * Lz[j]
                        radial_basis = dLr[i] * Lz[j] / self.radial_halfwidth
                        axial_basis = Lr[i] * dLz[j] / self.axial_halfwidth
                        potential += basis * snapshot["potentials"][i, j]
                        potential_r += radial_basis * snapshot["potentials"][i, j]
                        potential_z += axial_basis * snapshot["potentials"][i, j]
                        pressure_amp += basis * snapshot["pressures"][i, j]
                phase = m * theta + kr * (radius - rc) + kz * (z - zc) + phi0
                sine, cosine = np.sin(phase), np.cos(phase)
                pr, pt, pz = potential
                pr_r, pt_r, pz_r = potential_r
                pr_z, pt_z, pz_z = potential_z
                curl_r = (
                    (-m * pz / radius + kz * pt) * envelope * sine
                    - (envelope_z * pt + envelope * pt_z) * cosine
                )
                curl_theta = (
                    (-kz * pr + kr * pz) * envelope * sine
                    + (envelope_z * pr + envelope * pr_z
                       - envelope_r * pz - envelope * pz_r) * cosine
                )
                curl_z = (
                    (-kr * pt + m * pr / radius) * envelope * sine
                    + (envelope_r * pt + envelope * pt_r
                       + envelope * pt / radius) * cosine
                )
                factor = float(np.sqrt(weight))
                total_cyl += factor * np.array([curl_r, curl_theta, curl_z])
                pressure[index] += factor * envelope * pressure_amp * cosine
            ca, sa = x / radius, ycart / radius
            velocity[index] = [
                ca * total_cyl[0] - sa * total_cyl[1],
                sa * total_cyl[0] + ca * total_cyl[1],
                total_cyl[2],
            ]
        return velocity, pressure

    def diagnostics(self, tau: float) -> dict:
        snapshots = [self._snapshot(index, float(tau)) for index in range(len(self.waves))]
        return {
            "tau": float(tau),
            "mode_amplitude_norm_ratios_at_center": [
                row["amplitude_norm_ratio"] for row in snapshots
            ],
            "max_node_transversality_error": max(
                row["max_transversality_error"] for row in snapshots
            ),
            "source_kelvin_rhs_max_abs_difference": self.source_rhs_error,
            "modes": [
                {
                    "m": int(mode["m"]),
                    "center": snapshots[index]["center_state"][:2].tolist(),
                    "center_normal": snapshots[index]["center_normal"].tolist(),
                }
                for index, mode in enumerate(self.waves)
            ],
        }


class AddedWaveField:
    """Add a wave's velocity and pressure to a common base field."""

    def __init__(self, base, wave):
        self.base = base
        self.wave = wave
        self.nu = base.nu

    def fields(self, points, tau):
        velocity, pressure = self.base.fields(points, tau)
        wave_velocity, wave_pressure = self.wave.fields(points, tau)
        return velocity + wave_velocity, pressure + wave_pressure


def holdout_points(wave: TransportedKelvinPatch, axis: np.ndarray, angles: np.ndarray) -> np.ndarray:
    center = wave._snapshot(0, wave.tau0)["center_state"]
    rc, zc = center[:2]
    rows = []
    for xi in axis:
        for eta in axis:
            radius = rc + wave.radial_halfwidth * xi
            height = zc + wave.axial_halfwidth * eta
            rows.extend(
                [
                    [radius * np.cos(angle), radius * np.sin(angle), height]
                    for angle in angles
                ]
            )
    return np.asarray(rows, dtype=float)


def metrics(field, points: np.ndarray, tau: float) -> dict:
    hs = 0.0005 * np.sqrt(field.nu * tau)
    ht = 0.0001 * tau
    velocity, gradient, temporal = kinematics(field, points, tau, hs, ht)
    residual = temporal + np.einsum("nij,nj->ni", gradient, velocity)
    norms = np.linalg.norm(residual, axis=1)
    return {
        "max_momentum": float(np.max(norms)),
        "rms_momentum": float(np.sqrt(np.mean(norms**2))),
        "max_fd_divergence": float(np.max(np.abs(np.trace(gradient, axis1=1, axis2=2)))),
        "max_speed": float(np.max(np.linalg.norm(velocity, axis=1))),
        "points": int(len(points)),
        "spatial_step": float(hs),
        "temporal_step": float(ht),
    }


def run(output_name: str = "transported_kelvin_patch.json") -> dict:
    source = json.loads((COMPACT / SOURCE_NAME).read_text())
    if abs(float(source["tau"]) - TAU0) > 1e-14:
        raise ValueError(f"registered source tau changed: {source['tau']}")
    base = load_robust_candidate(0.1)
    frozen_wave = LocalizedCurlWave(
        source,
        radial_halfwidth=RADIAL_HALF_WIDTH,
        axial_halfwidth=AXIAL_HALF_WIDTH,
        time_halfwidth=TIME_HALF_WIDTH,
    )
    dynamic_wave = TransportedKelvinPatch(source, base)
    frozen_field = WavePerturbedField(base, frozen_wave)
    dynamic_field = AddedWaveField(base, dynamic_wave)
    # Two Cartesian samples (two spatial holdout locations and one angle)
    # keep the independent full-momentum check bounded on this expensive base.
    angles = np.array([0.0], dtype=float)
    points = holdout_points(dynamic_wave, HOLDOUT_NODES, angles)
    rows = []
    for tau in (TAU0, TAU1):
        rows.append(
            {
                "tau": tau,
                "base": metrics(base, points, tau),
                "frozen": metrics(frozen_field, points, tau),
                "transported_kelvin_patch": metrics(dynamic_field, points, tau),
                "patch_diagnostics": dynamic_wave.diagnostics(tau),
            }
        )
    report = {
        "source": SOURCE_NAME,
        "registered_times": [TAU0, TAU1],
        "patch": {
            "radial_halfwidth": RADIAL_HALF_WIDTH,
            "axial_halfwidth": AXIAL_HALF_WIDTH,
            "time_halfwidth": TIME_HALF_WIDTH,
            "node_axis": PATCH_NODES.tolist(),
            "holdout_axis": HOLDOUT_NODES.tolist(),
            "angles_per_node": len(angles),
            "mode_count": len(dynamic_wave.waves),
            "angular_modes": [int(mode["m"]) for mode in dynamic_wave.waves],
        },
        "rows": rows,
        "construction": {
            "phase": "transported local cylindrical phase kinematics with backward-time center and wavevector signs, using cached physical jets",
            "amplitude": "same projected Kelvin physical-coordinate RHS algebra as kelvin_amplitude_transport.derivative, integrated by bounded fixed-step RK4 on cached jets; no solve_ivp trajectory is claimed",
            "velocity": "analytic cylindrical curl of an interpolated vector potential over the moving 3 by 3 patch",
            "pressure": "interpolated p_proj = 2*n dot (J*a)/|n|^2, with the sign chosen to match the -2*n*(n dot J*a)/|n|^2 amplitude projection; experimental proxy, not a pressure PDE",
            "divergence": "analytic curl identity; max_fd_divergence is an independent finite-difference diagnostic",
            "coefficient_table": "nine times at the initial 3 by 3 node locations; node motion is evolved against these cached jets to bound runtime",
        },
        "scope": "Finite 3 by 3 moving patch in the registered tau neighborhood. Complete Cartesian finite-difference momentum is evaluated on disjoint spatial holdouts at two times. No global moment closure, finite-energy extension, critical-time continuation, normalized-chart estimate, or PDE acceptance is claimed.",
        "paper_status": "Sections 7--9 are only an experimental physical-coordinate analogue here; this is not a paper-equivalent construction.",
        "accepted": False,
    }
    output = COMPACT / output_name
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"output": str(output), "rows": rows}, indent=2), flush=True)
    return report


if __name__ == "__main__":
    run()
