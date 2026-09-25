"""C2 time blends of ST073 radial orders with complete local momentum.

The field u=(1-w)u_left+w*u_right, with similarly blended pressure,
remains solenoidal.
Its residual includes the exact time-switch and nonlinear cross terms.
"""

from dataclasses import replace
import json

import numpy as np

from radial_continuation import FullRadialField, ROOT
from full_radial import coordinates


TRANSFERS = ((8, 10, 10.0, 12.0), (10, 12, 18.0, 20.0))


class AdaptiveOrderCore:
    """Callable local C2 order schedule; physical time is t=0.5-tau."""

    def __init__(self, k_max=20):
        source = ROOT / "NS_ST073_Full_Local_Recurrence/data/ST073-V.json"
        base = FullRadialField.load(source)
        self.k_max = int(k_max)
        self.fields = {order: FullRadialField(replace(base.p, order=order, k_max=self.k_max))
                       for order in (8, 10, 12)}
        self.nu = base.nu
        self.h = base.h

    def evaluate_similarity(self, X, eta, tau, angle=0.0):
        if not (0.5*2.0**(-self.k_max) <= tau <= 0.5*2.0**(-6)):
            raise ValueError(f"Outside exploratory k=6..{self.k_max} time window")
        X, eta, angle = np.broadcast_arrays(np.asarray(X, float),
                                             np.asarray(eta, float),
                                             np.asarray(angle, float))
        shape = X.shape
        x, e, a = X.ravel(), eta.ravel(), angle.ravel()
        k = -np.log2(2.0*tau)
        if k <= 10.0:
            left_order = right_order = 8
            w = dw_dt = 0.0
        elif k < 12.0:
            left_order, right_order = 8, 10
            w, dw_dt = weight(k, tau, 10.0, 12.0)
        elif k <= 18.0:
            left_order = right_order = 10
            w = dw_dt = 0.0
        elif k < 20.0:
            left_order, right_order = 10, 12
            w, dw_dt = weight(k, tau, 18.0, 20.0)
        else:
            left_order = right_order = 12
            w = dw_dt = 0.0
        left = self.fields[left_order].evaluate_similarity(x, e, tau, a)
        if left_order == right_order:
            data = {key: left[key] for key in ("velocity", "pressure", "residual", "divergence")}
        else:
            right = self.fields[right_order].evaluate_similarity(x, e, tau, a)
            cross_cyl = np.array([
                difference_advection(self.fields[left_order], self.fields[right_order],
                                     float(xi), float(ei), tau)
                for xi, ei in zip(x, e)
            ])
            ca, sa = np.cos(a), np.sin(a)
            cross = np.stack([ca*cross_cyl[:, 0] - sa*cross_cyl[:, 1],
                              sa*cross_cyl[:, 0] + ca*cross_cyl[:, 1],
                              cross_cyl[:, 2]], axis=1)
            delta = right["velocity"] - left["velocity"]
            data = dict(
                velocity=(1.0-w)*left["velocity"] + w*right["velocity"],
                pressure=(1.0-w)*left["pressure"] + w*right["pressure"],
                residual=((1.0-w)*left["residual"] + w*right["residual"]
                          + dw_dt*delta - w*(1.0-w)*cross),
                divergence=(1.0-w)*left["divergence"] + w*right["divergence"],
            )
        data["velocity"] = data["velocity"].reshape(shape + (3,))
        data["pressure"] = data["pressure"].reshape(shape)
        data["residual"] = data["residual"].reshape(shape + (3,))
        data["divergence"] = data["divergence"].reshape(shape)
        data.update(k=float(k), orders=[left_order, right_order], weight=w)
        return data

    def evaluate(self, points, t):
        points = np.asarray(points, float)
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("Expected (n,3) Cartesian points")
        tau = 0.5 - float(t)
        src = points / np.sqrt(self.nu)
        coord = coordinates(np.hypot(src[:, 0], src[:, 1]), src[:, 2],
                            np.full(len(points), tau), self.h)
        return self.evaluate_similarity(coord["X"], coord["eta"], tau,
                                        np.arctan2(src[:, 1], src[:, 0]))

    def velocity_xyz(self, x, y, z, t):
        x, y, z = np.broadcast_arrays(x, y, z)
        points = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
        return self.evaluate(points, t)["velocity"].reshape(x.shape + (3,))


def weight(k, tau, left_k, right_k):
    s = np.clip((k - left_k) / (right_k - left_k), 0.0, 1.0)
    w = s**3 * (10.0 - 15.0*s + 6.0*s*s)
    dw_ds = 30.0*s*s*(1.0 - s)**2
    dw_dt = dw_ds / ((right_k - left_k) * tau * np.log(2.0))
    return float(w), float(dw_dt)


def local_coefficients(field, X, eta, tau):
    q = tau / (1.0 - eta*eta)
    S = 2.0*q
    coeff = field.coefficients(float(eta), float(q))
    values = []
    for i in (0, 1, 2):
        values.append([
            np.polynomial.polynomial.polyval(X, coeff[i, :, 0]),
            np.polynomial.polynomial.polyval(X, coeff[i, :, 1]),
            np.polynomial.polynomial.polyval(
                X, np.polynomial.polynomial.polyder(coeff[i, :, 0])) / S,
        ])
    return np.array(values, float), S


def difference_advection(field_left, field_right, X, eta, tau):
    c_left, S = local_coefficients(field_left, X, eta, tau)
    c_right, _ = local_coefficients(field_right, X, eta, tau)
    d = c_right - c_left
    (a, az, ass), (b, bz, bs), (c, cz, cs) = d
    r = np.sqrt(S*X)
    adv = np.array([
        r*(a*a + 2.0*S*X*a*ass + c*az - b*b),
        r*(2.0*a*b + 2.0*S*X*a*bs + c*bz),
        2.0*S*X*a*cs + c*cz,
    ])
    return np.sqrt(field_left.nu)*adv


def row(field_left, field_right, left_k, right_k, k, x_nodes, eta_nodes):
    tau = 0.5*2.0**(-k)
    w, dw_dt = weight(k, tau, left_k, right_k)
    xx, ee = np.meshgrid(x_nodes, eta_nodes, indexing="ij")
    x, eta = xx.ravel(), ee.ravel()
    d_left = field_left.evaluate_similarity(x, eta, tau)
    d_right = field_right.evaluate_similarity(x, eta, tau)
    delta_u = d_right["velocity"] - d_left["velocity"]
    cross = np.array([difference_advection(field_left, field_right, float(xi), float(ei), tau)
                      for xi, ei in zip(x, eta)])
    switch_term = dw_dt*delta_u
    nonlinear_term = -w*(1.0-w)*cross
    residual = ((1.0-w)*d_left["residual"] + w*d_right["residual"]
                + switch_term + nonlinear_term)
    magnitude = np.linalg.norm(residual, axis=1)
    peak = int(np.argmax(magnitude))
    return dict(
        k=k, tau=tau, weight=w, weight_time_derivative=dw_dt,
        left_order=field_left.p.order, right_order=field_right.p.order,
        left_momentum_max=float(np.max(np.linalg.norm(d_left["residual"], axis=1))),
        right_momentum_max=float(np.max(np.linalg.norm(d_right["residual"], axis=1))),
        blended_momentum_max=float(magnitude[peak]),
        peak_X=float(x[peak]), peak_eta=float(eta[peak]),
        time_switch_term_max=float(np.max(np.linalg.norm(switch_term, axis=1))),
        nonlinear_cross_term_max=float(np.max(np.linalg.norm(nonlinear_term, axis=1))),
        order_velocity_difference_max=float(np.max(np.linalg.norm(delta_u, axis=1))),
        blended_divergence_max=float(np.max(np.abs((1.0-w)*d_left["divergence"] + w*d_right["divergence"]))),
    )


def run():
    source = ROOT / "NS_ST073_Full_Local_Recurrence/data/ST073-V.json"
    base = FullRadialField.load(source)
    fields = {order: FullRadialField(replace(base.p, order=order, k_max=20))
              for order in (8, 10, 12)}
    x_nodes = np.linspace(0.001, 1.0/64.0, 8)
    eta_nodes = np.array([-0.3, -0.15, 0.0, 0.15, 0.3])
    transfers = []
    for left_order, right_order, left_k, right_k in TRANSFERS:
        ks = np.linspace(left_k, right_k, 9)
        rows = [row(fields[left_order], fields[right_order], left_k, right_k,
                    float(k), x_nodes, eta_nodes) for k in ks]
        transfers.append(dict(
            orders=[left_order, right_order], k_interval=[left_k, right_k],
            rows=rows,
            sampled_gate_1e_minus_3=all(r["blended_momentum_max"] < 1e-3 for r in rows),
        ))
    report = dict(
        source=str(source.relative_to(ROOT)),
        construction="C2 quintic time blends of adjacent local radial orders and pressures; exact residual identity includes w_t*(u_right-u_left)-w*(1-w)*((u_right-u_left) dot grad)(u_right-u_left).",
        x_nodes=x_nodes.tolist(), eta_nodes=eta_nodes.tolist(),
        transfers=transfers,
        sampled_gate_1e_minus_3=all(t["sampled_gate_1e_minus_3"] for t in transfers),
        arithmetic_note="On this host longdouble is float64; order-12 absolute residuals near k=18..20 approach cancellation noise. The blended samples are exploratory, not precision-certified bounds.",
        scope="40 local similarity locations per time at nine times per transfer; both k-bands exceed ST073-V registered k<=6. No continuous supremum, spatial-volume L2, outer matching, global energy or critical-time limit.",
        pde_validated=False, scale_recursion_established=False,
    )
    output = ROOT / "order_blend_scale.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output), transfers=transfers), indent=2))
    return report


if __name__ == "__main__":
    run()
