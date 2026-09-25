"""Fit one shared 24-mode bridge correction at two dyadic scales."""

import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares

from adaptive_core_join_screen import AdaptiveRadialAdapter
from affine_momentum import jets, momentum, combine
from joined_field import JoinedField, independent_fd
from radial_continuation import ROOT
from heat_exterior import physical
from wide_modes import WideJointModes


class CachedAdaptiveBridge:
    def __init__(self, joined):
        self.joined = joined
        self.inner = joined.inner
        self.nu = joined.nu
        self.ratio = joined.outer_ratio
        self.join_X = joined.join_X
        self.cache = {}

    def fields(self, points, tau):
        points = np.asarray(points, float)
        times = np.broadcast_to(np.asarray(tau, float), (len(points),))
        key = (points.shape, points.tobytes(), times.tobytes())
        if key not in self.cache:
            self.cache[key] = self.joined.fields(points, times)
        u, p = self.cache[key]
        return u.copy(), p.copy()


def sample_points(inner, ratio, k, etas, radial_order):
    tau = 0.5*2.0**(-k)
    roots, _ = leggauss(radial_order)
    y = (roots+1.0)/2.0
    all_X = []
    all_eta = []
    for eta in etas:
        all_X.extend((inner.p.X_max*(1.0+(ratio-1.0)*y)**2).tolist())
        all_eta.extend([eta]*len(y))
    points = inner.from_similarity(np.asarray(all_X), np.asarray(all_eta), tau)
    return points, tau


def run():
    inner = AdaptiveRadialAdapter()
    tau0 = 0.5*2.0**(-6)
    point = inner.from_similarity([1.0/64.0], [0.0], tau0)
    heat_amplitude = float(inner.evaluate(point, tau0)["velocity"][0, 1]
                           / physical(point, tau0, c=1.0)["velocity"][0, 1])
    ratio = 16.0
    joined = JoinedField(inner=inner, join_X=1.0/64.0,
                         heat_amplitude=heat_amplitude, outer_ratio=ratio)
    cached = CachedAdaptiveBridge(joined)
    zero = WideJointModes(cached, np.zeros(24))
    slices = []
    for k in (11.0, 19.0):
        points, tau = sample_points(inner, ratio, k, (-0.2, 0.0, 0.2), 5)
        args = (points, tau, 0.0005*np.sqrt(inner.nu*tau), 0.0001*tau)
        baseline = jets(zero, *args)
        mode_jets = []
        for j in range(24):
            a = np.zeros(24)
            a[j] = 1.0
            unit = jets(WideJointModes(cached, a), *args)
            mode_jets.append(tuple(u-b for u, b in zip(unit, baseline)))
        modes = tuple(np.stack([m[i] for m in mode_jets]) for i in range(3))
        scale = float(np.max(np.linalg.norm(momentum(baseline), axis=1)))
        slices.append(dict(k=k, args=args, baseline=baseline,
                           modes=modes, scale=scale))

    def objective(a):
        pieces = [(momentum(combine(s["baseline"], s["modes"], a))/s["scale"]).ravel()
                  for s in slices]
        return np.concatenate([*pieces, 0.002*a])

    def jacobian(a):
        blocks = []
        for s in slices:
            u, grad, _ = combine(s["baseline"], s["modes"], a)
            du, dg, dl = s["modes"]
            derivative = (dl+np.einsum("knij,nj->kni", dg, u)
                          +np.einsum("nij,knj->kni", grad, du))
            blocks.append(derivative.reshape(24, -1).T/s["scale"])
        return np.vstack([*blocks, 0.002*np.eye(24)])

    fit = least_squares(objective, np.zeros(24), jac=jacobian,
                        bounds=(-10.0, 10.0), max_nfev=100,
                        ftol=1e-10, xtol=1e-10, gtol=1e-10)
    fitted = WideJointModes(cached, fit.x)
    training = []
    for s in slices:
        predicted = momentum(combine(s["baseline"], s["modes"], fit.x))
        direct, divergence = independent_fd(fitted, *s["args"])
        training.append(dict(
            k=s["k"], baseline_max=s["scale"],
            fitted_predicted_max=float(np.max(np.linalg.norm(predicted, axis=1))),
            fitted_direct_max=float(np.max(np.linalg.norm(direct, axis=1))),
            surrogate_direct_component_difference=float(np.max(np.abs(predicted-direct))),
            divergence_max=float(np.max(np.abs(divergence))),
        ))
    holdouts = []
    for k in (10.5, 11.5, 18.5, 19.5):
        points, tau = sample_points(inner, ratio, k, (-0.1, 0.1), 6)
        h = 0.0005*np.sqrt(inner.nu*tau)
        ht = 0.0001*tau
        old, _ = independent_fd(zero, points, tau, h, ht)
        new, divergence = independent_fd(fitted, points, tau, h, ht)
        holdouts.append(dict(
            k=k, baseline_max=float(np.max(np.linalg.norm(old, axis=1))),
            fitted_max=float(np.max(np.linalg.norm(new, axis=1))),
            divergence_max=float(np.max(np.abs(divergence))),
        ))
    report = dict(
        source="AdaptiveOrderCore -> JoinedField ratio16 + shared 24 compact modes",
        coefficient_order="12 poloidal/pressure modes, then 12 swirl modes",
        train_scales=[11.0, 19.0], train_eta=[-0.2, 0.0, 0.2],
        amplitudes=fit.x.tolist(), fit_nfev=fit.nfev,
        fit_success=bool(fit.success), initial_cost=float(np.dot(objective(np.zeros(24)), objective(np.zeros(24)))/2),
        final_cost=float(fit.cost), training=training, holdouts=holdouts,
        baseline_growth_exponent_per_k=float(np.log2(training[1]["baseline_max"] / training[0]["baseline_max"])
                                             /(training[1]["k"]-training[0]["k"])),
        fitted_growth_exponent_per_k=float(np.log2(training[1]["fitted_direct_max"] / training[0]["fitted_direct_max"])
                                           /(training[1]["k"]-training[0]["k"])),
        scope="24 shared axisymmetric endpoint-preserving modes; two scale slices, 15 train points each and 12 points on each of four time holdouts. Full Cartesian FD momentum. No moment/stress cone constraint, spatial-volume L2 or whole-field acceptance.",
        accepted=False, pde_validated=False,
    )
    output = ROOT / "adaptive_join_multiscale_fit.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output), training=training,
                          holdouts=holdouts, fit_nfev=fit.nfev), indent=2))
    return report


if __name__ == "__main__":
    run()
