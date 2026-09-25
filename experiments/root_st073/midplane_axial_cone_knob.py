"""Screen a local mean-profile knob for wider two-scale axial cone support.

Perturbs the first radial-window odd poloidal coefficient at k=11 and k=19.
Reports physical cone nodes and the exact sampled moment damage; this is a
sensitivity experiment, not a moment-preserving redesigned background.
"""

import json

import numpy as np

from adaptive_bridge_moment_fit import moment_slices, outer_moments
from adaptive_bridge_recursive_defect import build_fields
from radial_continuation import ROOT
from radial_peak_cone import operator, stress_primitive
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


KNOTS = (11., 15., 19.)
ETA_NODES = (-.075, -.05, -.0125, .025, .05)
Y = .325


def cone_row(field, r, z, tau, eta):
    velocity, gradient, residual = operator(field, np.array([[r, 0., z]]), tau)
    F = velocity[0, 1] / r
    shear = np.array([gradient[0, 1, 0] - F, gradient[0, 2, 0]])
    magnitude = float(np.linalg.norm(shear))
    row = dict(eta=eta, momentum_norm=float(np.linalg.norm(residual[0])),
               lambda_squared=None, cone_pass=False)
    if magnitude == 0:
        return row
    N = shear / magnitude
    K = np.array([-N[1], N[0]])
    lam2 = float(-2 * F * N[0] * (2 * F * N[0] + magnitude))
    row["lambda_squared"] = lam2
    if lam2 <= 0 or abs(2 * F * N[0]) <= 1e-14:
        return row
    target = stress_primitive(field, float(r), float(z), tau, order=12)
    dot_n, dot_k = float(target @ N), float(target @ K)
    row.update(target_dot_N=dot_n, target_dot_K=dot_k)
    if abs(dot_n) <= 1e-14:
        return row
    ratio = abs(np.sqrt(lam2) * dot_k / (2 * F * N[0] * dot_n))
    row.update(cone_ratio=float(ratio),
               cone_pass=bool(dot_n < 0 and ratio < 1))
    return row


def run():
    inner, fields = build_fields()
    base = fields["two_sided_cone"]
    old = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    amplitudes = np.asarray(old["constrained_l4"]["amplitudes"], float)
    shaped = amplitudes.reshape(3, 2, 3, 2)
    indices = (int(np.ravel_multi_index((0, 1, 0, 1), shaped.shape)),
               int(np.ravel_multi_index((2, 1, 0, 1), shaped.shape)))

    def field(a):
        return SeparatedMomentModes(base, a,
                                    windows=RADIAL_WINDOWS_THREE,
                                    knots=KNOTS)

    zero = field(amplitudes)
    units = [field(amplitudes + np.eye(len(amplitudes))[index])
             for index in indices]
    moment_data = moment_slices(inner, base, zero, orders=KNOTS,
                                unit_fields=units)
    baseline_moments = np.array([outer_moments(s, np.zeros(2))
                                 for s in moment_data])
    cases = []
    for multiplier in (.5, 1., 1.5):
        new = amplitudes.copy()
        deltas = [(multiplier - 1) * amplitudes[index] for index in indices]
        for index, delta in zip(indices, deltas):
            new[index] += delta
        candidate = field(new)
        moments = np.array([outer_moments(s, deltas) for s in moment_data])
        scales = []
        for k in (11, 19):
            tau = .5 * 2.**-k
            X = np.full(len(ETA_NODES),
                        inner.p.X_max * (1 + 15*Y)**2)
            points = inner.from_similarity(X, np.asarray(ETA_NODES), tau)
            rows = [cone_row(candidate, float(point[0]), float(point[2]),
                             tau, eta)
                    for point, eta in zip(points, ETA_NODES)]
            scales.append(dict(k=k, rows=rows,
                               passing_eta=[row["eta"] for row in rows
                                            if row["cone_pass"]]))
        cases.append(dict(multiplier=multiplier,
                          coefficient_deltas=deltas,
                          max_absolute_sampled_moment=float(np.max(np.abs(moments))),
                          max_absolute_sampled_moment_drift=float(
                              np.max(np.abs(moments-baseline_moments))),
                          moments=moments.tolist(), scales=scales))
        print(f"finished multiplier={multiplier}: passing={[s['passing_eta'] for s in scales]}",
              flush=True)
    report = dict(source="First-window odd-poloidal axial cone sensitivity",
                  adjusted_flat_indices=list(indices),
                  baseline_coefficients=[float(amplitudes[i]) for i in indices],
                  y=Y, eta_nodes=ETA_NODES,
                  baseline_max_absolute_sampled_moment=float(
                      np.max(np.abs(baseline_moments))),
                  cases=cases,
                  scope="Two-scale physical cone analogue at five eta nodes plus exact sampled physical outer-moment values at three knot scales. Candidate coefficient perturbations are not moment-preserving or full residual improvements; passing nodes do not certify continuous support.",
                  accepted=False, scale_recursion_established=False)
    (ROOT / "midplane_axial_cone_knob.json").write_bytes(
        (json.dumps(report, indent=2) + "\n").encode())
    return report


if __name__ == "__main__":
    run()
