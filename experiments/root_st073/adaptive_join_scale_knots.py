"""Fit independent bridge coefficients at two scales and C2-transfer them."""

import json

import numpy as np
from scipy.optimize import least_squares

from adaptive_core_join_screen import AdaptiveRadialAdapter
from adaptive_join_multiscale_fit import CachedAdaptiveBridge, sample_points
from affine_momentum import jets, momentum, combine
from joined_field import JoinedField, independent_fd
from radial_continuation import ROOT
from heat_exterior import physical
from wide_modes import WideJointModes


K_LEFT, K_RIGHT = 11.0, 19.0


def knot_weight(k):
    s = np.clip((k-K_LEFT)/(K_RIGHT-K_LEFT), 0.0, 1.0)
    return float(s**3*(10.0-15.0*s+6.0*s*s))


class InterpolatedWideModes:
    def __init__(self, base, left, right):
        self.base = base
        self.left = np.asarray(left, float)
        self.right = np.asarray(right, float)
        self.nu = base.nu

    def coefficients(self, tau):
        k = -np.log2(2.0*float(tau))
        w = knot_weight(k)
        return (1.0-w)*self.left+w*self.right

    def fields(self, points, tau):
        times = np.broadcast_to(np.asarray(tau, float), (len(points),))
        if not np.all(times == times[0]):
            raise ValueError("Expected one remaining time per field call")
        return WideJointModes(self.base, self.coefficients(float(times[0]))).fields(points, times)


def assemble_slice(base, zero, inner, ratio, k):
    points, tau = sample_points(inner, ratio, k, (-0.2, 0.0, 0.2), 5)
    args = (points, tau, 0.0005*np.sqrt(inner.nu*tau), 0.0001*tau)
    baseline = jets(zero, *args)
    mode_jets = []
    for j in range(24):
        unit_amplitude = np.zeros(24)
        unit_amplitude[j] = 1.0
        unit = jets(WideJointModes(base, unit_amplitude), *args)
        mode_jets.append(tuple(u-b for u, b in zip(unit, baseline)))
    modes = tuple(np.stack([m[i] for m in mode_jets]) for i in range(3))
    scale = float(np.max(np.linalg.norm(momentum(baseline), axis=1)))
    return dict(k=k, args=args, baseline=baseline, modes=modes, scale=scale)


def fit_slice(s, start):
    def objective(a):
        return np.r_[(momentum(combine(s["baseline"], s["modes"], a))/s["scale"]).ravel(),
                     0.002*a]

    def jacobian(a):
        u, grad, _ = combine(s["baseline"], s["modes"], a)
        du, dg, dl = s["modes"]
        derivative = dl+np.einsum("knij,nj->kni", dg, u)+np.einsum("nij,knj->kni", grad, du)
        return np.vstack([derivative.reshape(24, -1).T/s["scale"], 0.002*np.eye(24)])

    return least_squares(objective, start, jac=jacobian,
                         bounds=(-10.0, 10.0), max_nfev=100,
                         ftol=1e-10, xtol=1e-10, gtol=1e-10)


def sampled_max(field, points, tau, nu):
    h = 0.0005*np.sqrt(nu*tau)
    ht = 0.0001*tau
    residual, divergence = independent_fd(field, points, tau, h, ht)
    return float(np.max(np.linalg.norm(residual, axis=1))), float(np.max(np.abs(divergence))), residual


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
    shared_amplitude = np.asarray(json.loads(
        (ROOT/"adaptive_join_multiscale_fit.json").read_text())["amplitudes"], float)
    shared = WideJointModes(cached, shared_amplitude)
    slices = [assemble_slice(cached, zero, inner, ratio, k) for k in (K_LEFT, K_RIGHT)]
    fits = [fit_slice(s, shared_amplitude) for s in slices]
    dynamic = InterpolatedWideModes(cached, fits[0].x, fits[1].x)
    training = []
    for s, fit in zip(slices, fits):
        baseline_max = s["scale"]
        shared_max, _, _ = sampled_max(shared, s["args"][0], s["args"][1], inner.nu)
        dynamic_max, divergence, _ = sampled_max(dynamic, s["args"][0], s["args"][1], inner.nu)
        training.append(dict(k=s["k"], baseline_max=baseline_max,
                             shared_max=shared_max, knot_max=dynamic_max,
                             knot_divergence_max=divergence,
                             fit_nfev=fit.nfev, fit_success=bool(fit.success)))
    holdouts = []
    for k in (10.5, 11.5, 14.0, 15.0, 17.0, 18.5, 19.5):
        points, tau = sample_points(inner, ratio, k, (-0.1, 0.1), 6)
        baseline_max, _, _ = sampled_max(zero, points, tau, inner.nu)
        shared_max, _, _ = sampled_max(shared, points, tau, inner.nu)
        dynamic_max, divergence, dynamic_residual = sampled_max(dynamic, points, tau, inner.nu)
        frozen = WideJointModes(cached, dynamic.coefficients(tau))
        frozen_max, _, frozen_residual = sampled_max(frozen, points, tau, inner.nu)
        holdouts.append(dict(k=k, baseline_max=baseline_max,
                             shared_max=shared_max, dynamic_max=dynamic_max,
                             frozen_same_coefficients_max=frozen_max,
                             time_switch_residual_difference_max=float(np.max(
                                 np.linalg.norm(dynamic_residual-frozen_residual, axis=1))),
                             divergence_max=divergence))
    report = dict(
        source="AdaptiveOrderCore -> JoinedField ratio16 + C2 two-knot 24-mode bridge",
        knot_scales=[K_LEFT, K_RIGHT],
        left_amplitudes=fits[0].x.tolist(), right_amplitudes=fits[1].x.tolist(),
        coefficient_max_difference=float(np.max(np.abs(fits[1].x-fits[0].x))),
        coefficient_relative_l2_difference=float(np.linalg.norm(fits[1].x-fits[0].x)/np.linalg.norm(fits[0].x)),
        training=training, holdouts=holdouts,
        knot_growth_exponent_per_k=float(np.log2(training[1]["knot_max"]/training[0]["knot_max"])
                                          /(training[1]["k"]-training[0]["k"])),
        max_time_switch_to_dynamic_ratio=float(max(
            row["time_switch_residual_difference_max"]/row["dynamic_max"] for row in holdouts)),
        scope="Independent coefficients fitted at k=11 and 19, C2 quintic interpolation in k, full Cartesian momentum sampled on knot and off-knot grids. No radial moment, stress-cone, volume-L2 or whole-field acceptance.",
        accepted=False, pde_validated=False, scale_recursion_established=False,
    )
    output = ROOT/"adaptive_join_scale_knots.json"
    output.write_bytes((json.dumps(report, indent=2)+"\n").encode())
    print(json.dumps(dict(output=str(output), training=training,
                          holdouts=holdouts), indent=2))
    return report


if __name__ == "__main__":
    run()
