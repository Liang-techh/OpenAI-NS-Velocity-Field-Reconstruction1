"""One-mode azimuthal momentum repair for the adaptive wide bridge."""

import json

import numpy as np
from scipy.optimize import minimize_scalar

from adaptive_core_join_screen import AdaptiveRadialAdapter
from affine_momentum import jets, momentum, combine
from joined_field import JoinedField, independent_fd
from radial_continuation import ROOT
from heat_exterior import physical


def sample_points(inner, ratio, eta, k):
    tau = 0.5*2.0**(-k)
    y = np.array([0.25, 0.5, 0.75])
    radius_ratio = 1.0+(ratio-1.0)*y
    X = (1.0/64.0)*radius_ratio**2
    points = inner.from_similarity(X, np.full(len(y), eta), tau)
    return points, tau


def jet_pair(base, mode1, points, tau):
    args = (points, tau, 0.0005*np.sqrt(base.nu*tau), 0.0001*tau)
    zero = jets(base, *args)
    unit = jets(mode1, *args)
    modes = tuple(np.stack([u-b]) for b, u in zip(zero, unit))
    return zero, modes, args


def peak(base, modes, a):
    return float(np.max(np.linalg.norm(momentum(combine(base, modes, np.array([a]))), axis=1)))


def run():
    inner = AdaptiveRadialAdapter()
    tau0 = 0.5*2.0**(-6)
    point = inner.from_similarity([1.0/64.0], [0.0], tau0)
    amplitude = float(inner.evaluate(point, tau0)["velocity"][0, 1]
                      / physical(point, tau0, c=1.0)["velocity"][0, 1])
    ratio = 16.0
    def candidate(a):
        return JoinedField(inner=inner, join_X=1.0/64.0,
                           heat_amplitude=amplitude, outer_ratio=ratio,
                           swirl_bubble_amplitude=a)
    base, mode1 = candidate(0.0), candidate(1.0)
    points, tau = sample_points(inner, ratio, 0.0, 11.0)
    train, train_modes, args = jet_pair(base, mode1, points, tau)
    grid = np.linspace(-40.0, 40.0, 161)
    grid_values = np.array([peak(train, train_modes, a) for a in grid])
    index = int(np.argmin(grid_values))
    low, high = grid[max(index-1, 0)], grid[min(index+1, len(grid)-1)]
    optimum = minimize_scalar(lambda a: peak(train, train_modes, a),
                              bounds=(low, high), method="bounded",
                              options=dict(xatol=1e-8))
    fitted = candidate(float(optimum.x))
    direct, div = independent_fd(fitted, *args)
    predicted = momentum(combine(train, train_modes, np.array([optimum.x])))
    holdouts = []
    for eta, k in ((-0.2, 11.0), (0.2, 11.0), (0.0, 11.5)):
        pts, t = sample_points(inner, ratio, eta, k)
        old, _ = independent_fd(base, pts, t, 0.0005*np.sqrt(inner.nu*t), 0.0001*t)
        new, d = independent_fd(fitted, pts, t, 0.0005*np.sqrt(inner.nu*t), 0.0001*t)
        holdouts.append(dict(eta=eta, k=k,
                             baseline_max=float(np.max(np.linalg.norm(old, axis=1))),
                             fitted_max=float(np.max(np.linalg.norm(new, axis=1))),
                             divergence_max=float(np.max(np.abs(d)))))
    report = dict(
        source="AdaptiveOrderCore -> JoinedField with endpoint-preserving swirl bubble",
        ratio=ratio, training_k=11.0, training_eta=0.0,
        amplitude=float(optimum.x),
        train_baseline_max=peak(train, train_modes, 0.0),
        train_fitted_predicted_max=peak(train, train_modes, float(optimum.x)),
        train_fitted_direct_max=float(np.max(np.linalg.norm(direct, axis=1))),
        predicted_direct_vector_max_difference=float(np.max(np.abs(predicted-direct))),
        train_divergence_max=float(np.max(np.abs(div))),
        holdouts=holdouts,
        scope="One endpoint-preserving axisymmetric swirl-bubble coefficient; three train locations and nine holdout locations. Pressure and poloidal bridge unchanged. This is a bounded capacity screen, not complete-field acceptance.",
        accepted=False, pde_validated=False,
    )
    output = ROOT / "adaptive_join_swirl_capacity.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    run()
