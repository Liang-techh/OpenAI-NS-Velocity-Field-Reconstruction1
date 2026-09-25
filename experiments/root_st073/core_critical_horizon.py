"""Explore when an unchanged ST073-V core exceeds the momentum gate.

This extends the registered time window for diagnosis only. It samples the
local polynomial field; it gives no bound on the complete matched field.
"""

from dataclasses import replace
import json

import numpy as np

from radial_continuation import FullRadialField, ROOT


def sampled_row(field, k, x_nodes, eta_nodes):
    tau = 0.5 * 2.0 ** (-k)
    xx, ee = np.meshgrid(x_nodes, eta_nodes, indexing="ij")
    x, eta = xx.ravel(), ee.ravel()
    data = field.evaluate_similarity(x, eta, tau)
    magnitude = np.linalg.norm(data["residual"], axis=1)
    peak = int(np.argmax(magnitude))
    q = tau / (1.0 - eta * eta)
    scaled = data["residual_cylindrical"] / np.sqrt(field.nu)
    scaled[:, 0] *= q ** 1.5
    scaled[:, 1:] *= (q ** (field.A + 1.0))[:, None]
    term_sum = sum(np.linalg.norm(data[name], axis=1) for name in
                   ("time_derivative", "advection", "pressure_gradient", "viscous_term"))
    roundoff_reference = np.finfo(np.longdouble).eps * term_sum
    return dict(
        k=k,
        tau=tau,
        sampled_momentum_max=float(magnitude[peak]),
        peak_X=float(x[peak]),
        peak_eta=float(eta[peak]),
        normalized_momentum_max=float(np.max(np.linalg.norm(scaled, axis=1))),
        term_sum_at_peak=float(term_sum[peak]),
        roundoff_reference_at_peak=float(roundoff_reference[peak]),
        residual_to_roundoff_reference_at_peak=float(magnitude[peak] / roundoff_reference[peak]),
        sampled_divergence_max=float(np.max(np.abs(data["divergence"]))),
    )


def run():
    source = ROOT / "NS_ST073_Full_Local_Recurrence/data/ST073-V.json"
    baseline = FullRadialField.load(source)
    orders = (8, 10, 12, 14)
    ks = (6, 8, 10, 12, 14, 15, 16, 18)
    x_nodes = np.array([0.001, 0.004, 0.008, 0.012, 1.0 / 64.0])
    eta_nodes = np.array([-0.3, -0.15, 0.0, 0.15, 0.3])
    rows = {}
    for order in orders:
        field = FullRadialField(replace(baseline.p, order=order, k_max=max(ks)))
        rows[str(order)] = [sampled_row(field, k, x_nodes, eta_nodes) for k in ks]
    peak8 = np.array([r["sampled_momentum_max"] for r in rows["8"]])
    peak10 = np.array([r["sampled_momentum_max"] for r in rows["10"]])
    report = dict(
        source=str(source.relative_to(ROOT)),
        modification="k_max is extended to 18 for diagnosis; orders 8, 10, 12 and 14 compare radial truncations with all other parameters fixed.",
        viscosity=baseline.nu,
        target_momentum_max=1e-3,
        x_nodes=x_nodes.tolist(),
        eta_nodes=eta_nodes.tolist(),
        k_values=list(ks),
        rows_by_radial_order=rows,
        relative_order_8_to_10_peak_difference=(np.abs(peak8 - peak10) / np.maximum(peak10, 1e-300)).tolist(),
        arithmetic_note="On this Windows host np.longdouble has float64 precision. eps times the sum of PDE-term norms is only a roundoff scale indicator, not an error bound. Residuals comparable to it are unresolved.",
        first_sampled_failure_k_by_order={str(order): next((r["k"] for r in rows[str(order)] if r["sampled_momentum_max"] >= 1e-3), None) for order in orders},
        scope="25 local similarity points per time, X<=1/64, |eta|<=.3; k>6 is an exploratory extension of registered time. A sampled failure disproves the unchanged field's 1e-3 gate at that point, but sampled success is not a bound. No outer matching, forcing construction, or global energy claim.",
        pde_validated=False,
        scale_recursion_established=False,
    )
    output = ROOT / "core_critical_horizon.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(dict(output=str(output), first_sampled_failure_k_by_order=report["first_sampled_failure_k_by_order"], rows_by_radial_order=rows), indent=2))
    return report


if __name__ == "__main__":
    run()
