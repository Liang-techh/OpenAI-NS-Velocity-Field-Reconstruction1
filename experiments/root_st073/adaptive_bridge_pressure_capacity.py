"""Fit extra compact pressure modes on the corrected adaptive bridge."""

import json

import numpy as np
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.legendre import Legendre
from scipy.optimize import lsq_linear

from adaptive_core_join_screen import AdaptiveRadialAdapter
from adaptive_join_multiscale_fit import CachedAdaptiveBridge, sample_points
from joined_field import JoinedField, independent_fd, coordinates
from radial_continuation import ROOT
from heat_exterior import physical
from wide_modes import WideJointModes


MODE_PAIRS = [(j, m) for j in range(6) for m in range(4)]
ENVELOPE = 64*Poly([0, 1])**3*Poly([1, -1])**3
BASIS = [ENVELOPE*Legendre.basis(j).convert(kind=Poly)(Poly([-1, 2]))
         for j in range(6)]


def pressure_and_gradient(points, tau, inner, ratio, coefficients):
    points = np.asarray(points, float)
    rootnu = np.sqrt(inner.nu)
    r = np.hypot(points[:, 0], points[:, 1])
    coord = coordinates(r/rootnu, points[:, 2]/rootnu,
                        np.broadcast_to(np.asarray(tau, float), r.shape), inner.h)
    q, eta = np.asarray(coord["q"]), np.asarray(coord["eta"])
    qz = np.asarray(coord["q_z"])/rootnu
    etaz = np.asarray(coord["eta_z"])/rootnu
    ri = np.sqrt(2*inner.nu*q*inner.p.X_max)
    width = (ratio-1.0)*ri
    raw_y = (r-ri)/width
    active = (raw_y > 0.0) & (raw_y < 1.0)
    y = np.clip(raw_y, 0.0, 1.0)
    yz = -(1.0+(ratio-1.0)*y)*qz/(2.0*q*(ratio-1.0))
    scale = inner.nu*q**(-2.0*inner.A)
    pressure = np.zeros(len(points))
    grad_r = np.zeros(len(points))
    grad_z = np.zeros(len(points))
    for amplitude, (j, m) in zip(coefficients, MODE_PAIRS):
        if amplitude == 0:
            continue
        B, By = BASIS[j](y), BASIS[j].deriv()(y)
        ax = (eta/.3)**m
        ax_z = np.zeros_like(eta) if m == 0 else m*(eta/.3)**(m-1)*etaz/.3
        pressure += amplitude*scale*B*ax
        grad_r += amplitude*scale*By*ax/width
        grad_z += amplitude*scale*((-2.0*inner.A)*qz/q*B*ax
                                    + By*yz*ax+B*ax_z)
    pressure = np.where(active, pressure, 0.0)
    grad_r = np.where(active, grad_r, 0.0)
    grad_z = np.where(active, grad_z, 0.0)
    safe = np.where(r > 0, r, 1.0)
    gradient = np.column_stack((grad_r*points[:, 0]/safe,
                                grad_r*points[:, 1]/safe, grad_z))
    return pressure, gradient


class PressureExtendedField:
    def __init__(self, base, inner, ratio, amplitudes):
        self.base, self.inner, self.ratio = base, inner, ratio
        self.amplitudes = np.asarray(amplitudes, float)
        self.nu = inner.nu

    def fields(self, points, tau):
        u, p = self.base.fields(points, tau)
        addition, _ = pressure_and_gradient(points, tau, self.inner,
                                            self.ratio, self.amplitudes)
        return u, p+addition


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
    prior = np.asarray(json.loads((ROOT/"adaptive_join_multiscale_fit.json").read_text())["amplitudes"])
    base = WideJointModes(cached, prior)
    training_data = []
    for k in (11.0, 19.0):
        points, tau = sample_points(inner, ratio, k, (-0.2, 0.0, 0.2), 5)
        h, ht = 0.0005*np.sqrt(inner.nu*tau), 0.0001*tau
        residual, _ = independent_fd(base, points, tau, h, ht)
        columns = []
        for j in range(len(MODE_PAIRS)):
            amplitude = np.zeros(len(MODE_PAIRS))
            amplitude[j] = 1.0
            _, gradient = pressure_and_gradient(points, tau, inner, ratio, amplitude)
            columns.append(gradient)
        matrix = np.stack(columns, axis=2)
        scale = float(np.max(np.linalg.norm(residual, axis=1)))
        training_data.append(dict(k=k, points=points, tau=tau, residual=residual,
                                  matrix=matrix, scale=scale))
    A = np.vstack([d["matrix"].reshape(-1, len(MODE_PAIRS))/d["scale"]
                   for d in training_data]+[0.002*np.eye(len(MODE_PAIRS))])
    b = np.concatenate([-d["residual"].ravel()/d["scale"]
                        for d in training_data]+[np.zeros(len(MODE_PAIRS))])
    fit = lsq_linear(A, b, bounds=(-20.0, 20.0), max_iter=500, tol=1e-11)
    fitted = PressureExtendedField(base, inner, ratio, fit.x)
    training = []
    for d in training_data:
        predicted = d["residual"]+np.einsum("nij,j->ni", d["matrix"], fit.x)
        direct, divergence = independent_fd(fitted, d["points"], d["tau"],
                                            0.0005*np.sqrt(inner.nu*d["tau"]),
                                            0.0001*d["tau"])
        training.append(dict(
            k=d["k"], baseline_max=d["scale"],
            pressure_fit_predicted_max=float(np.max(np.linalg.norm(predicted, axis=1))),
            pressure_fit_direct_max=float(np.max(np.linalg.norm(direct, axis=1))),
            pressure_invariant_angular_lower_bound=float(np.max(np.abs(d["residual"][:, 1]))),
            analytic_direct_component_difference=float(np.max(np.abs(predicted-direct))),
            divergence_max=float(np.max(np.abs(divergence))),
        ))
    holdouts = []
    for k in (10.5, 11.5, 18.5, 19.5):
        points, tau = sample_points(inner, ratio, k, (-0.1, 0.1), 6)
        h, ht = 0.0005*np.sqrt(inner.nu*tau), 0.0001*tau
        old, _ = independent_fd(base, points, tau, h, ht)
        new, divergence = independent_fd(fitted, points, tau, h, ht)
        holdouts.append(dict(
            k=k, baseline_max=float(np.max(np.linalg.norm(old, axis=1))),
            pressure_fit_max=float(np.max(np.linalg.norm(new, axis=1))),
            angular_lower_bound=float(np.max(np.abs(old[:, 1]))),
            divergence_max=float(np.max(np.abs(divergence))),
        ))
    report = dict(
        source="Shared 24-mode adaptive bridge plus 24 extra compact pressure modes",
        pressure_mode_pairs=MODE_PAIRS,
        pressure_amplitudes=fit.x.tolist(), fit_success=bool(fit.success),
        fit_cost=float(fit.cost), active_bounds=int(np.sum(np.abs(fit.x)>19.99)),
        training=training, holdouts=holdouts,
        scope="Pressure-only correction, same coefficients at k=11 and 19; scalar potential leaves velocity, divergence, and angular momentum residual invariant. Three axial by five radial train points per scale. No continuous pressure PDE, moment constraints, cone, volume L2 or whole-field acceptance.",
        accepted=False, pde_validated=False,
    )
    output = ROOT/"adaptive_bridge_pressure_capacity.json"
    output.write_bytes((json.dumps(report, indent=2)+"\n").encode())
    print(json.dumps(dict(output=str(output), training=training,
                          holdouts=holdouts), indent=2))
    return report


if __name__ == "__main__":
    run()
