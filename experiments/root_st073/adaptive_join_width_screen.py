"""Screen radial heat-bridge width for the scale-adaptive ST073 core."""

import json

import numpy as np

from adaptive_core_join_screen import AdaptiveRadialAdapter
from joined_field import JoinedField
from radial_continuation import ROOT
from heat_exterior import physical


def momentum_terms(field, points, tau, hspace, htime):
    u, _ = field.fields(points, tau)
    gradient = np.empty((len(points), 3, 3))
    grad_pressure = np.empty_like(u)
    laplacian = np.zeros_like(u)
    for axis in range(3):
        step = np.zeros(3)
        step[axis] = hspace
        um2, pm2 = field.fields(points-2*step, tau)
        um, pm = field.fields(points-step, tau)
        up, pp = field.fields(points+step, tau)
        up2, pp2 = field.fields(points+2*step, tau)
        gradient[:, :, axis] = (um2-8*um+8*up-up2)/(12*hspace)
        grad_pressure[:, axis] = (pm2-8*pm+8*pp-pp2)/(12*hspace)
        laplacian += (-up2+16*up-30*u+16*um-um2)/(12*hspace*hspace)
    du_dt = -(field.fields(points, tau-2*htime)[0]
              -8*field.fields(points, tau-htime)[0]
              +8*field.fields(points, tau+htime)[0]
              -field.fields(points, tau+2*htime)[0])/(12*htime)
    advection = np.einsum("nij,nj->ni", gradient, u)
    terms = dict(time=du_dt, advection=advection,
                 pressure=grad_pressure, viscosity=-field.nu*laplacian)
    residual = sum(terms.values())
    divergence = np.trace(gradient, axis1=1, axis2=2)
    return residual, divergence, terms


def run():
    inner = AdaptiveRadialAdapter()
    tau0 = 0.5*2.0**(-6)
    point = inner.from_similarity([1.0/64.0], [0.0], tau0)
    amplitude = float(inner.evaluate(point, tau0)["velocity"][0, 1]
                      / physical(point, tau0, c=1.0)["velocity"][0, 1])
    rows = []
    for ratio in (2.0, 4.0, 8.0, 16.0, 32.0, 64.0):
        field = JoinedField(inner=inner, join_X=1.0/64.0,
                            heat_amplitude=amplitude, outer_ratio=ratio)
        k = 11.0
        tau = 0.5*2.0**(-k)
        y = np.array([0.25, 0.5, 0.75])
        radius_ratio = 1.0+(ratio-1.0)*y
        X = (1.0/64.0)*radius_ratio**2
        eta = np.zeros(len(y))
        points = inner.from_similarity(X, eta, tau)
        residual, divergence, terms = momentum_terms(
            field, points, tau, 0.0005*np.sqrt(inner.nu*tau), 0.0001*tau)
        peak = int(np.argmax(np.linalg.norm(residual, axis=1)))
        rows.append(dict(
            outer_radius_ratio=ratio, outer_X=(1.0/64.0)*ratio**2,
            y=y.tolist(), X=X.tolist(),
            momentum_norms=np.linalg.norm(residual, axis=1).tolist(),
            momentum_max=float(np.max(np.linalg.norm(residual, axis=1))),
            divergence_max=float(np.max(np.abs(divergence))),
            peak_term_vectors={name: value[peak].tolist() for name, value in terms.items()},
            peak_residual_vector=residual[peak].tolist(),
        ))
    report = dict(
        source="AdaptiveOrderCore -> JoinedField -> radial heat exterior",
        k=11.0, tau=0.5*2.0**(-11), eta=0.0,
        heat_amplitude=amplitude,
        rows=rows,
        best_sampled_ratio=min(rows, key=lambda r:r["momentum_max"])["outer_radius_ratio"],
        scope="Three interior bridge points per width on one axial/time slice. Changing width changes the sampled physical locations. No optimized pressure, wave stress, volume L2, or continuous domain bound.",
        accepted=False, pde_validated=False,
    )
    output = ROOT / "adaptive_join_width_screen.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output), rows=rows), indent=2))
    return report


if __name__ == "__main__":
    run()
