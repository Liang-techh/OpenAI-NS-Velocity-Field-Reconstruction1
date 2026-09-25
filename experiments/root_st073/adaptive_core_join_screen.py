"""Attach the scale-adaptive ST073 core to the existing heat bridge.

This probes the actual dynamic cost of a finite-energy-oriented radial
transition; it does not assert acceptance or global finite energy.
"""

from dataclasses import replace
import json

import numpy as np

from joined_field import JoinedField, independent_fd, coordinates
from order_blend_scale import AdaptiveOrderCore, weight
from radial_continuation import ROOT
from heat_exterior import physical


class AdaptiveRadialAdapter:
    """Expose the scheduled core's coefficient jets to JoinedField."""

    def __init__(self):
        self.core = AdaptiveOrderCore()
        self.fields_by_order = self.core.fields
        base = self.fields_by_order[12]
        self.p = replace(base.p, X_max=1.0/64.0, k_max=20)
        self.nu, self.h, self.A, self.D = base.nu, base.h, base.A, base.D

    def active(self, tau):
        k = -np.log2(2.0*tau)
        if k <= 10.0:
            return 8, 8, 0.0
        if k < 12.0:
            w, _ = weight(k, tau, 10.0, 12.0)
            return 8, 10, w
        if k <= 18.0:
            return 10, 10, 0.0
        if k < 20.0:
            w, _ = weight(k, tau, 18.0, 20.0)
            return 10, 12, w
        return 12, 12, 0.0

    def coefficients(self, eta, q):
        tau = float(q)*(1.0-float(eta)**2)
        left, right, w = self.active(tau)
        l = self.fields_by_order[left].coefficients(float(eta), float(q))
        if left == right:
            return l
        r = self.fields_by_order[right].coefficients(float(eta), float(q))
        l = np.pad(l, ((0, 0), (0, r.shape[1]-l.shape[1]), (0, 0)))
        return (1.0-w)*l + w*r

    def evaluate(self, points, tau):
        points = np.asarray(points, float)
        src = points/np.sqrt(self.nu)
        coord = coordinates(np.hypot(src[:, 0], src[:, 1]), src[:, 2],
                            np.broadcast_to(np.asarray(tau, float), (len(points),)), self.h)
        # The bridge calls this adapter with one tau per point.
        times = np.broadcast_to(np.asarray(tau, float), (len(points),))
        if not np.all(times == times[0]):
            raise ValueError("Adapter expects a shared remaining time")
        return self.core.evaluate_similarity(coord["X"], coord["eta"],
                                             float(times[0]), np.arctan2(src[:, 1], src[:, 0]))

    def from_similarity(self, X, eta, tau, angle=0):
        return self.fields_by_order[8].from_similarity(X, eta, tau, angle)


def run():
    inner = AdaptiveRadialAdapter()
    tau0 = 0.5*2.0**(-6)
    match_point = inner.from_similarity([1.0/64.0], [0.0], tau0)
    core_swirl = inner.evaluate(match_point, tau0)["velocity"][0, 1]
    heat_swirl = physical(match_point, tau0, c=1.0)["velocity"][0, 1]
    amplitude = float(core_swirl/heat_swirl)
    field = JoinedField(inner=inner, join_X=1.0/64.0,
                        heat_amplitude=amplitude, outer_ratio=2.0)
    rows = []
    for k in (6.5, 11.0, 19.0):
        tau = 0.5*2.0**(-k)
        radius_ratio = np.array([1.25, 1.75, 1.25, 1.75])
        eta = np.array([0.0, 0.0, 0.25, 0.25])
        X = (1.0/64.0)*radius_ratio**2
        points = inner.from_similarity(X, eta, tau)
        step_space = 0.0005*np.sqrt(inner.nu*tau)
        step_time = 0.0001*tau
        residual, div = independent_fd(field, points, tau, step_space, step_time)
        core_points = inner.from_similarity(np.full(len(eta), 0.01), eta, tau)
        core_data = inner.evaluate(core_points, tau)
        rows.append(dict(
            k=k, tau=tau, X=X.tolist(), eta=eta.tolist(),
            active_orders=list(inner.active(tau)[:2]),
            core_momentum_max=float(np.max(np.linalg.norm(core_data["residual"], axis=1))),
            bridge_momentum_norms=np.linalg.norm(residual, axis=1).tolist(),
            bridge_momentum_max=float(np.max(np.linalg.norm(residual, axis=1))),
            bridge_divergence_max=float(np.max(np.abs(div))),
            spatial_step=step_space, remaining_time_step=step_time,
        ))
    report = dict(
        source="AdaptiveOrderCore -> JoinedField -> radial heat exterior",
        heat_amplitude=amplitude,
        join_X=1.0/64.0, outer_X=1.0/16.0,
        rows=rows,
        scope="Four bridge locations at three times, one finite-difference resolution. The heat exterior is axially unlocalized and the join has no momentum acceptance. No whole-space finite energy or critical-time force claim.",
        accepted=False, pde_validated=False, scale_recursion_established=False,
    )
    output = ROOT / "adaptive_core_join_screen.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output), rows=rows), indent=2))
    return report


if __name__ == "__main__":
    run()
