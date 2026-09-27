"""Grouped/vectorized evaluator for the radial ``JoinedField`` bridge.

``JoinedField.fields`` historically evaluates source coordinates and bridge
coefficients once per point.  For radial panels, many points share exactly the
same ``(tau, z)``.  ``GroupedJoinedField`` evaluates those shared quantities
once per group and vectorizes the radial polynomial and Cartesian rotation.
The wrapper is opt-in; the global ``joined_field`` implementation is left
unchanged.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
UPSTREAM = ROOT / "NS_ST073_Full_Local_Recurrence" / "upstream"
if str(UPSTREAM) not in sys.path:
    sys.path.insert(0, str(UPSTREAM))

from heat_exterior import physical  # noqa: E402
from joined_field import JoinedField, coordinates  # noqa: E402
from poloidal_bridge import coefficients, septic  # noqa: E402
from swirl_bridge import quintic, traces  # noqa: E402


class GroupedJoinedField:
    """Drop-in ``fields`` wrapper that groups points by exact ``(tau, z)``."""

    def __init__(self, joined):
        self.joined = joined
        for name in (
            "inner", "nu", "join_X", "outer_ratio", "swirl_bubble_amplitude",
            "outer_swirl_bubble_amplitude", "slope_swirl_bubble_amplitude",
            "shear_swirl_amplitude", "shear_swirl_cycles", "c",
        ):
            if hasattr(joined, name):
                setattr(self, name, getattr(joined, name))

    def _group_fields(self, points, t):
        """Evaluate one exact-z/time group, preserving original branch rules."""

        points = np.asarray(points, dtype=float)
        count = len(points)
        if count == 0:
            return np.empty((0, 3), dtype=float), np.empty(0, dtype=float)
        sn = np.sqrt(self.nu)
        radii = np.hypot(points[:, 0], points[:, 1])
        z = float(points[0, 2])

        # q and eta do not depend on r.  Calling the source inversion at the
        # axis gives the same q/eta as the original per-point calls while
        # reducing the 60 bisection loops to one vectorized scalar call.
        coord = coordinates(0.0, z / sn, float(t), self.inner.h)
        eta = float(coord["eta"])
        q = float(coord["q"])
        if (abs(eta) > self.inner.p.eta_max
                or not 0.5 * 2.0 ** (-self.inner.p.k_max) <= t <= 0.5):
            raise ValueError("Only registered axial/time slab supported")

        ri, ro, left, a = coefficients(
            self.inner, eta, float(t), self.join_X, self.outer_ratio
        )
        inner_mask = radii <= ri
        outer_mask = radii >= ro
        bridge_mask = ~(inner_mask | outer_mask)
        velocity = np.empty((count, 3), dtype=float)
        pressure = np.empty(count, dtype=float)

        if np.any(inner_mask):
            inner_data = self.inner.evaluate(points[inner_mask], float(t))
            velocity[inner_mask] = inner_data["velocity"]
            pressure[inner_mask] = inner_data["pressure"]
        if np.any(outer_mask):
            outer_data = physical(
                points[outer_mask], float(t), c=self.c
            )
            velocity[outer_mask] = outer_data["velocity"]
            pressure[outer_mask] = outer_data["pressure"]

        if np.any(bridge_mask):
            bridge_points = points[bridge_mask]
            bridge_radii = radii[bridge_mask]
            width = ro - ri
            X = self.join_X
            rho = ri / sn
            co = self.inner.coefficients(eta, q)

            def val(comp, derivative=0, column=0):
                polynomial = np.polynomial.polynomial.polyder(
                    co[comp, :, column], derivative
                )
                return float(np.polynomial.polynomial.polyval(X, polynomial)
                             / (2.0 * q) ** derivative)

            A = val(0)
            C_s, C_ss, C_sss = val(2, 1), val(2, 2), val(2, 3)
            C_z, C_sz, C_ssz = val(2, 0, 1), val(2, 1, 1), val(2, 2, 1)
            fourth = (6.0 * C_s + 24.0 * rho**2 * C_ss
                      + 8.0 * rho**4 * C_sss) / sn
            fixed_z = np.array([
                -ri * sn * rho * A,
                ri * C_z,
                C_z + 2.0 * rho**2 * C_sz,
                (6.0 * rho * C_sz + 4.0 * rho**3 * C_ssz) / sn,
            ])
            ri_z = ri * float(coord["q_z"]) / (2.0 * q * sn)
            total_z = fixed_z + np.r_[left[1:], fourth] * ri_z
            modified = np.array([
                total_z[j] + j * left[j] * ri_z / ri for j in range(4)
            ])
            a_z = septic(modified, width)
            y = (bridge_radii - ri) / width
            y_prime = -(1.0 + (self.outer_ratio - 1.0) * y) * ri_z / width
            psi_r = np.polynomial.polynomial.polyval(
                y, np.polynomial.polynomial.polyder(a)
            ) / width
            psi_z = (
                np.polynomial.polynomial.polyval(y, a_z)
                + np.polynomial.polynomial.polyval(
                    y, np.polynomial.polynomial.polyder(a)
                ) * y_prime
            )

            _, _, swirl_left, swirl_right = traces(
                self.inner, eta, float(t), self.c,
                self.join_X, self.outer_ratio
            )
            swirl_coefficients = quintic(swirl_left, swirl_right, width)
            tangential = np.polynomial.polynomial.polyval(y, swirl_coefficients)
            tangential += (
                sn * q ** (-self.inner.A) * self.swirl_bubble_amplitude
                * 64.0 * y**3 * (1.0 - y)**3
            )
            tangential += (
                sn * q ** (-self.inner.A)
                * self.outer_swirl_bubble_amplitude
                * (64.0 * y**3 * (1.0 - y)**3 / 0.421875)
                * (y / 0.75)**8
            )
            if self.slope_swirl_bubble_amplitude:
                y0 = (np.sqrt(1.0 / self.join_X) - 1.0) / (self.outer_ratio - 1.0)
                tangential += (
                    sn * q ** (-self.inner.A)
                    * self.slope_swirl_bubble_amplitude
                    * 64.0 * y**3 * (1.0 - y)**3 * (y - y0)
                )
            if self.shear_swirl_amplitude:
                y0 = (np.sqrt(1.0 / self.join_X) - 1.0) / (self.outer_ratio - 1.0)
                phase = 2.0 * np.pi * self.shear_swirl_cycles * (y - y0)
                tangential += (
                    sn * q ** (-self.inner.A)
                    * self.shear_swirl_amplitude
                    * 64.0 * y**3 * (1.0 - y)**3
                    * np.sin(phase) / (2.0 * np.pi * self.shear_swirl_cycles)
                )

            outer_point = np.array([[ro, 0.0, z]])
            outer_data = physical(outer_point, float(t), c=self.c)
            pressure_left = np.array([
                self.nu * val(3),
                2.0 * sn * rho * val(3, 1),
                2.0 * val(3, 1) + 4.0 * rho**2 * val(3, 2),
            ])
            pressure_right = np.array([
                float(outer_data["pressure"][0]),
                swirl_right[0] ** 2 / ro,
                2.0 * swirl_right[0] * swirl_right[1] / ro
                - swirl_right[0] ** 2 / ro**2,
            ])
            pressure_coefficients = quintic(pressure_left, pressure_right, width)
            bridge_pressure = np.polynomial.polynomial.polyval(
                y, pressure_coefficients
            )
            safe = np.where(bridge_radii > 0.0, bridge_radii, 1.0)
            ca = bridge_points[:, 0] / safe
            sa = bridge_points[:, 1] / safe
            radial = -psi_z / safe
            axial = psi_r / safe
            velocity[bridge_mask] = np.column_stack((
                radial * ca - tangential * sa,
                radial * sa + tangential * ca,
                axial,
            ))
            pressure[bridge_mask] = bridge_pressure
        return velocity, pressure

    def fields(self, points, tau):
        points = np.asarray(points, dtype=float)
        times = np.broadcast_to(np.asarray(tau, dtype=float), (len(points),))
        if len(points) == 0:
            return np.empty((0, 3), dtype=float), np.empty(0, dtype=float)
        groups = {}
        for index, (point, t) in enumerate(zip(points, times)):
            groups.setdefault((float(t), float(point[2])), []).append(index)
        velocity = np.empty((len(points), 3), dtype=float)
        pressure = np.empty(len(points), dtype=float)
        for (t, _), indices in groups.items():
            values, pressures = self._group_fields(points[indices], t)
            velocity[indices] = values
            pressure[indices] = pressures
        return velocity, pressure


def install_in_field(root):
    """Wrap known ``CachedAdaptiveBridge`` nodes under a field hierarchy.

    The traversal is intentionally limited to common wrapper attributes and
    containers used by the repository.  It clears a bridge cache after
    replacement so old point-array results cannot mask the grouped evaluator.
    Returns the number of replacements.
    """

    from adaptive_join_multiscale_fit import CachedAdaptiveBridge

    seen = set()
    replacements = 0

    def visit(node):
        nonlocal replacements
        if node is None or id(node) in seen:
            return
        seen.add(id(node))
        if isinstance(node, CachedAdaptiveBridge):
            if not isinstance(node.joined, GroupedJoinedField):
                node.joined = GroupedJoinedField(node.joined)
                if hasattr(node, "cache"):
                    node.cache.clear()
                replacements += 1
            return
        if isinstance(node, dict):
            for value in node.values():
                visit(value)
            return
        if isinstance(node, (list, tuple, set)):
            for value in node:
                visit(value)
            return
        for name in (
            "base", "current", "modified", "reference", "unit", "compact",
            "baseline", "values", "pressures", "modes",
        ):
            if hasattr(node, name):
                visit(getattr(node, name))

    visit(root)
    return replacements


def _make_joined(inner):
    from heat_exterior import physical

    join_X = 1.0 / 64.0
    ratio = 16.0
    tau = 0.5 * 2.0 ** (-6.0)
    point = inner.from_similarity([join_X], [0.0], tau)
    heat_amplitude = float(
        inner.evaluate(point, tau)["velocity"][0, 1]
        / physical(point, tau, c=1.0)["velocity"][0, 1]
    )
    return JoinedField(
        inner=inner, join_X=join_X, heat_amplitude=heat_amplitude,
        outer_ratio=ratio, swirl_bubble_amplitude=0.13,
        outer_swirl_bubble_amplitude=-0.07,
        slope_swirl_bubble_amplitude=0.05,
        shear_swirl_amplitude=0.03, shear_swirl_cycles=16,
    )


def _comparison(old, new, points, tau):
    started = time.perf_counter()
    old_velocity, old_pressure = old.fields(points, tau)
    old_seconds = time.perf_counter() - started
    started = time.perf_counter()
    new_velocity, new_pressure = new.fields(points, tau)
    new_seconds = time.perf_counter() - started
    velocity_diff = np.abs(old_velocity - new_velocity)
    pressure_diff = np.abs(old_pressure - new_pressure)
    velocity_scale = max(float(np.max(np.abs(old_velocity))), 1.0e-30)
    pressure_scale = max(float(np.max(np.abs(old_pressure))), 1.0e-30)
    return {
        "point_count": int(len(points)),
        "old_seconds": float(old_seconds),
        "grouped_seconds": float(new_seconds),
        "speedup": float(old_seconds / new_seconds) if new_seconds else None,
        "velocity_max_abs_diff": float(np.max(velocity_diff)),
        "velocity_relative_max_diff": float(np.max(velocity_diff) / velocity_scale),
        "pressure_max_abs_diff": float(np.max(pressure_diff)),
        "pressure_relative_max_diff": float(np.max(pressure_diff) / pressure_scale),
    }


def run():
    from adaptive_core_join_screen import AdaptiveRadialAdapter

    inner = AdaptiveRadialAdapter()
    joined = _make_joined(inner)
    grouped = GroupedJoinedField(joined)
    sn = np.sqrt(inner.nu)

    # Correctness grid: core, both sides of each join, bridge, exterior,
    # nonzero angle, several times, and vector-valued tau.
    test_rows = []
    for k, eta in ((6.4, -0.2), (11.0, 0.0), (19.0, 0.2)):
        tau = 0.5 * 2.0 ** (-k)
        co = coordinates(0.0, eta * 0.0, tau, inner.h)
        ri, ro, _, _ = coefficients(inner, eta, tau, joined.join_X, joined.outer_ratio)
        radii = np.array([
            0.5 * ri, ri - 1.0e-12, ri + 1.0e-12,
            ri + 0.37 * (ro - ri), ro - 1.0e-12, ro + 1.0e-12,
            1.2 * ro,
        ])
        angles = np.linspace(0.2, 1.7, len(radii))
        points = np.column_stack((radii * np.cos(angles),
                                  radii * np.sin(angles),
                                  np.full(len(radii),
                                          inner.from_similarity(
                                              [inner.p.X_max * 2.0], [eta], tau
                                          )[0, 2])))
        row = _comparison(joined, grouped, points, tau)
        row.update({"k": float(k), "eta": float(eta), "tau_mode": "scalar"})
        test_rows.append(row)
        vector_tau = np.full(len(points), tau)
        row = _comparison(joined, grouped, points, vector_tau)
        row.update({"k": float(k), "eta": float(eta), "tau_mode": "vector"})
        test_rows.append(row)

    # Compare full spatial/time finite-difference jets at three off-axis
    # bridge points.  This catches small field differences that can be
    # amplified by the second spatial stencil in the momentum residual.
    k = 11.0
    tau = 0.5 * 2.0 ** (-k)
    eta = 0.0
    z = float(inner.from_similarity([inner.p.X_max * 2.0], [eta], tau)[0, 2])
    ri, ro, _, _ = coefficients(inner, eta, tau, joined.join_X, joined.outer_ratio)
    radii = ri + np.array([0.2, 0.5, 0.8]) * (ro - ri)
    angles = np.array([0.31, 1.07, 2.11])
    jet_points = np.column_stack((radii * np.cos(angles),
                                  radii * np.sin(angles),
                                  np.full(3, z)))
    from affine_momentum import jets, momentum  # noqa: E402
    hspace = 0.0005 * np.sqrt(inner.nu * tau)
    htime = 0.00025 * tau
    old_jets = jets(joined, jet_points, tau, hspace, htime)
    grouped_jets = jets(grouped, jet_points, tau, hspace, htime)
    old_residual = momentum(old_jets)
    grouped_residual = momentum(grouped_jets)
    jet_difference = np.concatenate([
        (old_jets[index] - grouped_jets[index]).ravel() for index in range(3)
    ])
    residual_difference = old_residual - grouped_residual
    residual_scale = max(float(np.max(np.abs(old_residual))), 1.0e-30)
    full_jets = {
        "point_count": 3,
        "k": k,
        "velocity_pressure_points": jet_points.tolist(),
        "field_jet_max_abs_diff": float(np.max(np.abs(jet_difference))),
        "momentum_max_abs_diff": float(np.max(np.abs(residual_difference))),
        "momentum_relative_max_diff": float(np.max(np.abs(residual_difference)) / residual_scale),
        "old_momentum_max_abs": float(np.max(np.abs(old_residual))),
    }

    from adaptive_join_multiscale_fit import CachedAdaptiveBridge  # noqa: E402
    cached_probe = CachedAdaptiveBridge(_make_joined(inner))
    holder = {"baseline": {"values": [cached_probe], "modes": []}}
    replacement_count = install_in_field(holder)
    install_test = {
        "replacement_count": int(replacement_count),
        "wrapped_type": type(cached_probe.joined).__name__,
        "cache_cleared": bool(cached_probe.cache == {}),
    }

    # Domain checks should remain errors for both implementations.
    errors = []
    for label, points, tau in (
        ("time_above_registered", np.array([[0.001, 0.0, 0.0]]), 0.5001),
        ("eta_outside_registered", np.array([[0.001, 0.0, 0.4]]), 0.5 * 2.0 ** (-5.5)),
    ):
        outcomes = []
        for evaluator in (joined, grouped):
            try:
                evaluator.fields(points, tau)
            except Exception as error:  # noqa: BLE001 - compare public errors
                outcomes.append(type(error).__name__ + ":" + str(error))
            else:
                outcomes.append("no_error")
        errors.append({"case": label, "old": outcomes[0], "grouped": outcomes[1],
                       "same_error_type": outcomes[0].split(":", 1)[0] == outcomes[1].split(":", 1)[0]})

    # Cold benchmark: four exact z groups and 500 radial/angle points per
    # group.  Each evaluator is a fresh instance; no old array cache is used.
    benchmark_points = []
    benchmark_tau = []
    tau = 0.5 * 2.0 ** (-11.0)
    for eta in (-0.35, -0.1, 0.15, 0.4):
        z = float(inner.from_similarity([inner.p.X_max * 2.0], [eta], tau)[0, 2])
        co = coordinates(0.0, z / sn, tau, inner.h)
        ri, ro, _, _ = coefficients(inner, eta, tau, joined.join_X, joined.outer_ratio)
        radial = np.linspace(0.35 * ri, 1.25 * ro, 500)
        angles = np.linspace(0.13, 2.7, 500)
        benchmark_points.extend(np.column_stack((radial * np.cos(angles),
                                                   radial * np.sin(angles),
                                                   np.full(500, z))))
        benchmark_tau.extend([tau] * 500)
    benchmark_points = np.asarray(benchmark_points)
    benchmark_tau = np.asarray(benchmark_tau)
    old_cold = _make_joined(inner)
    grouped_cold = GroupedJoinedField(_make_joined(inner))
    benchmark = _comparison(old_cold, grouped_cold, benchmark_points, benchmark_tau)
    benchmark.update({"group_count": 4, "points_per_group": 500,
                      "cold_instances": True})

    report = {
        "accepted": False,
        "pde_validated": False,
        "scope": "Mathematically equivalent grouped JoinedField evaluator and finite correctness/performance checks; no PDE, recursion, or global acceptance claim.",
        "wrapper": "GroupedJoinedField",
        "parameters": {
            "join_X": joined.join_X,
            "outer_ratio": joined.outer_ratio,
            "swirl_bubble_amplitude": joined.swirl_bubble_amplitude,
            "outer_swirl_bubble_amplitude": joined.outer_swirl_bubble_amplitude,
            "slope_swirl_bubble_amplitude": joined.slope_swirl_bubble_amplitude,
            "shear_swirl_amplitude": joined.shear_swirl_amplitude,
            "shear_swirl_cycles": joined.shear_swirl_cycles,
            "heat_amplitude": joined.c,
        },
        "correctness": test_rows,
        "domain_errors": errors,
        "full_jets": full_jets,
        "benchmark": benchmark,
        "install_helper_test": install_test,
        "install_helper": "install_in_field(root)",
    }
    output = ROOT / "grouped_joined_field.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "benchmark": benchmark,
                      "max_velocity_diff": max(row["velocity_max_abs_diff"] for row in test_rows),
                      "max_pressure_diff": max(row["pressure_max_abs_diff"] for row in test_rows)}, indent=2), flush=True)
    return report


if __name__ == "__main__":
    run()
