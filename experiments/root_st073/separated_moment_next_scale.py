"""Test one genuine new-knot transfer of the finite moment-repair ansatz.

Earlier knot coefficients remain fixed. Only the 12 coefficients at k=21
are fitted; k=20 and k=22 are independent transfer checks. This is a
finite-dimensional diagnostic, not the paper's residual-improvement cycle.
"""

import json

import numpy as np
from scipy.optimize import minimize

from adaptive_bridge_moment_fit import moment_slices, outer_moments
from adaptive_bridge_recursive_defect import build_fields
from affine_momentum import combine, momentum
from joined_field import independent_fd
from radial_continuation import ROOT
from separated_moment_fit import grid_slice
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


KNOTS = (11.0, 15.0, 19.0, 21.0)


def run():
    # Extend the exploratory registered slab without altering prior reports.
    inner, fields = build_fields(k_max=24)
    base = fields["two_sided_cone"]
    previous = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    fixed = np.asarray(previous["constrained_l4"]["amplitudes"], float)
    old_last = fixed[-12:].copy()

    def field(last):
        return SeparatedMomentModes(base, np.r_[fixed, last],
                                    windows=RADIAL_WINDOWS_THREE, knots=KNOTS)

    zero = field(np.zeros(12))
    units = [field(np.eye(12)[j]) for j in range(12)]
    moments21 = moment_slices(inner, base, zero, orders=(21.0,),
                              unit_fields=units)[0]
    grid21 = grid_slice(inner, zero, units, 21.0)
    baseline_moments = outer_moments(moments21, np.zeros(12))
    moment_scales = np.empty(4)
    for component in (0, 1):
        moment_scales[component::2] = max(
            np.max(np.abs(baseline_moments[component::2])), 1e-12)

    def normalized_moments(last):
        return outer_moments(moments21, last) / moment_scales

    def grid_norms(last):
        return np.linalg.norm(momentum(combine(
            grid21["baseline"], grid21["modes"], last)), axis=1)

    def objective(last):
        values = grid_norms(last) / grid21["scale"]
        return float(np.mean(values**4) + 1e-7 * np.dot(last - old_last,
                                                          last - old_last))

    fitted = minimize(objective, old_last, method="SLSQP",
                      bounds=[(-100.0, 100.0)] * 12,
                      constraints=[dict(type="eq", fun=normalized_moments)],
                      options=dict(maxiter=300, ftol=1e-12))

    def describe(last):
        return dict(coefficients=last.tolist(),
                    normalized_moment_max=float(np.max(np.abs(
                        normalized_moments(last)))),
                    momentum_grid_max=float(np.max(grid_norms(last))))

    transfer = []
    for k in (19.0, 20.0, 21.0, 22.0):
        tau = 0.5 * 2.0**-k
        # Direct full Cartesian finite differences on a fresh spatial grid.
        from adaptive_join_multiscale_fit import sample_points
        points, _ = sample_points(inner, 16.0, k,
                                  (-0.2, 0.0, 0.2), 5)
        rows = {}
        for name, coefficients in (("frozen", old_last),
                                   ("new_knot", fitted.x)):
            # Frozen is the four-knot extension that repeats the k=19 row.
            residual, divergence = independent_fd(
                field(coefficients), points, tau,
                0.0005 * np.sqrt(inner.nu * tau), 0.0001 * tau)
            rows[name] = dict(momentum_max=float(np.max(np.linalg.norm(
                residual, axis=1))),
                divergence_max=float(np.max(np.abs(divergence))))
        transfer.append(dict(k=k, tau=tau, **rows))

    report = dict(source="One-knot physical moment and momentum transfer",
                  fixed_knots=list(KNOTS[:-1]), new_knot=21.0,
                  fixed_coefficient_count=len(fixed),
                  original_at_new_scale=describe(old_last),
                  fitted_at_new_scale=describe(fitted.x),
                  optimizer_success=bool(fitted.success),
                  optimizer_message=str(fitted.message),
                  optimizer_iterations=int(fitted.nit),
                  transfer=transfer,
                  scope="One 12-dimensional new-knot fit, four sampled physical moment equalities at k=21, 15 momentum grid points at k=21, direct FD checks at k=19,20,21,22. Earlier coefficients fixed. No uniform cone, full-domain or iterative residual-improvement claim.",
                  accepted=False, scale_recursion_established=False)
    output = ROOT / "separated_moment_next_scale.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    run()
