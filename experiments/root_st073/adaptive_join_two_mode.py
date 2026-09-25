"""Coupled swirl and solenoidal poloidal bridge repair at one scale."""

import json

import numpy as np
from scipy.optimize import least_squares, minimize

from adaptive_core_join_screen import AdaptiveRadialAdapter
from adaptive_join_swirl_capacity import sample_points
from affine_momentum import jets, momentum, combine
from joined_field import JoinedField, independent_fd, coordinates
from radial_continuation import ROOT
from heat_exterior import physical


class PoloidalBubble:
    def __init__(self, base, amplitude):
        self.base = base
        self.amplitude = float(amplitude)
        self.nu = base.nu

    def fields(self, points, tau):
        points = np.asarray(points, float)
        u, p = self.base.fields(points, tau)
        if self.amplitude == 0.0:
            return u, p
        r = np.hypot(points[:, 0], points[:, 1])
        safe = np.where(r > 0.0, r, 1.0)
        rootnu = np.sqrt(self.nu)
        coord = coordinates(r/rootnu, points[:, 2]/rootnu,
                            np.broadcast_to(tau, r.shape), self.base.inner.h)
        q = np.asarray(coord["q"])
        qz = np.asarray(coord["q_z"])/rootnu
        ri = np.sqrt(2.0*self.nu*q*self.base.join_X)
        width = (self.base.outer_ratio-1.0)*ri
        y_raw = (r-ri)/width
        active = (y_raw > 0.0) & (y_raw < 1.0)
        y = np.clip(y_raw, 0.0, 1.0)
        bubble = 256.0*y**4*(1.0-y)**4
        bubble_y = 1024.0*y**3*(1.0-y)**3*(1.0-2.0*y)
        y_z = (-(1.0+(self.base.outer_ratio-1.0)*y)*qz
               /(2.0*q*(self.base.outer_ratio-1.0)))
        scale = self.amplitude*self.nu**1.5*(2.0*self.base.join_X)*q**(1.0-self.base.inner.A)
        psi_r = scale*bubble_y/width
        psi_z = scale*((1.0-self.base.inner.A)*qz/q*bubble + bubble_y*y_z)
        ur = np.where(active, -psi_z/safe, 0.0)
        uz = np.where(active, psi_r/safe, 0.0)
        u[:, 0] += ur*points[:, 0]/safe
        u[:, 1] += ur*points[:, 1]/safe
        u[:, 2] += uz
        return u, p


def run():
    inner = AdaptiveRadialAdapter()
    tau0 = 0.5*2.0**(-6)
    point = inner.from_similarity([1.0/64.0], [0.0], tau0)
    amplitude = float(inner.evaluate(point, tau0)["velocity"][0, 1]
                      / physical(point, tau0, c=1.0)["velocity"][0, 1])
    ratio = 16.0
    def candidate(a):
        bridge = JoinedField(inner=inner, join_X=1.0/64.0,
                             heat_amplitude=amplitude, outer_ratio=ratio,
                             swirl_bubble_amplitude=float(a[0]))
        return PoloidalBubble(bridge, float(a[1]))
    points = np.concatenate([sample_points(inner, ratio, eta, 11.0)[0]
                             for eta in (-0.2, 0.0, 0.2)])
    tau = 0.5*2.0**(-11)
    args = (points, tau, 0.0005*np.sqrt(inner.nu*tau), 0.0001*tau)
    baseline = jets(candidate([0.0, 0.0]), *args)
    unit_jets = [jets(candidate(a), *args) for a in ([1.0, 0.0], [0.0, 1.0])]
    modes = tuple(np.stack([unit[j]-baseline[j] for unit in unit_jets])
                  for j in range(3))
    scale = float(np.max(np.linalg.norm(momentum(baseline), axis=1)))
    def residual(a):
        return (momentum(combine(baseline, modes, a))/scale).ravel()
    def jacobian(a):
        u, grad, _ = combine(baseline, modes, a)
        du, dg, dl = modes
        derivative = dl+np.einsum("knij,nj->kni", dg, u)+np.einsum("nij,knj->kni", grad, du)
        return derivative.reshape(2, -1).T/scale
    fit = least_squares(residual, np.array([-2.35, 0.0]), jac=jacobian,
                        bounds=(-40.0, 40.0), max_nfev=100,
                        ftol=1e-11, xtol=1e-11, gtol=1e-11)
    def peak_objective(a):
        return float(np.max(np.linalg.norm(momentum(combine(baseline, modes, a)), axis=1)))
    peak_fit = minimize(peak_objective, fit.x, method="Powell",
                        bounds=[(-40.0, 40.0), (-40.0, 40.0)],
                        options=dict(xtol=1e-6, ftol=1e-8, maxiter=200))
    selected = peak_fit.x if peak_objective(peak_fit.x) < peak_objective(fit.x) else fit.x
    fitted = candidate(selected)
    direct, div = independent_fd(fitted, *args)
    predicted = momentum(combine(baseline, modes, selected))
    holdouts = []
    for eta, k in ((-0.25, 11.0), (0.25, 11.0), (0.0, 10.5), (0.0, 11.5)):
        pts, t = sample_points(inner, ratio, eta, k)
        h = 0.0005*np.sqrt(inner.nu*t)
        ht = 0.0001*t
        old, _ = independent_fd(candidate([0.0, 0.0]), pts, t, h, ht)
        new, d = independent_fd(fitted, pts, t, h, ht)
        holdouts.append(dict(eta=eta, k=k,
                             baseline_max=float(np.max(np.linalg.norm(old, axis=1))),
                             fitted_max=float(np.max(np.linalg.norm(new, axis=1))),
                             divergence_max=float(np.max(np.abs(d)))))
    report = dict(
        source="AdaptiveOrderCore -> wide JoinedField + swirl and poloidal bubbles",
        ratio=ratio, training_k=11.0, training_eta=[-0.2, 0.0, 0.2],
        amplitudes=selected.tolist(), least_squares_amplitudes=fit.x.tolist(),
        fit_nfev=fit.nfev, fit_success=bool(fit.success),
        peak_fit_success=bool(peak_fit.success),
        training_baseline_max=float(np.max(np.linalg.norm(momentum(baseline), axis=1))),
        training_fitted_direct_max=float(np.max(np.linalg.norm(direct, axis=1))),
        predicted_direct_vector_max_difference=float(np.max(np.abs(predicted-direct))),
        training_divergence_max=float(np.max(np.abs(div))),
        holdouts=holdouts,
        scope="Two endpoint-preserving axisymmetric modes at one scale; nine train and twelve holdout locations. Pressure and heat exterior unchanged. No volume L2 or continuum/domain acceptance.",
        accepted=False, pde_validated=False,
    )
    output = ROOT / "adaptive_join_two_mode.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    run()
