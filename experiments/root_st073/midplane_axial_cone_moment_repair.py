"""Repair sampled outer moments after an axial-cone widening mean change.

Fixes half-size first-window odd poloidal coefficients at k=11,19 and
optimizes the remaining 34 mean coefficients to restore all twelve
sampled physical moments. Then directly rechecks the midplane cone.
"""

import argparse
import json

import numpy as np
from scipy.optimize import minimize

from adaptive_bridge_moment_fit import moment_slices, outer_moments
from adaptive_bridge_recursive_defect import build_fields
from midplane_axial_cone_knob import ETA_NODES, KNOTS, Y, cone_row
from radial_continuation import ROOT
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


def run(all_knots=False):
    inner, fields = build_fields()
    base = fields["two_sided_cone"]
    old = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    amplitudes = np.asarray(old["constrained_l4"]["amplitudes"], float)
    shape = (3, 2, 3, 2)
    knot_indices = range(3) if all_knots else (0, 2)
    fixed_indices = [int(np.ravel_multi_index((index, 1, 0, 1), shape))
                     for index in knot_indices]
    free_indices = np.array([i for i in range(len(amplitudes))
                             if i not in fixed_indices])
    fixed_delta = np.zeros_like(amplitudes)
    fixed_delta[fixed_indices] = -.5 * amplitudes[fixed_indices]

    def field(a):
        return SeparatedMomentModes(base, a, windows=RADIAL_WINDOWS_THREE,
                                    knots=KNOTS)

    current = field(amplitudes)
    units = [field(amplitudes + np.eye(len(amplitudes))[i])
             for i in range(len(amplitudes))]
    slices = moment_slices(inner, base, current, orders=KNOTS,
                           unit_fields=units)
    moments0 = np.array([outer_moments(s, np.zeros(len(amplitudes)))
                         for s in slices]).ravel()
    unit_changes = np.array([
        (np.array([outer_moments(s, np.eye(len(amplitudes))[i])
                   for s in slices]).ravel() - moments0)
        for i in range(len(amplitudes))])
    component_scales = np.maximum(np.max(np.abs(unit_changes), axis=0),
                                  1e-8)

    def vector(z):
        delta = fixed_delta.copy()
        delta[free_indices] = z
        return delta

    def moment_defect(z):
        delta = vector(z)
        return (np.array([outer_moments(s, delta) for s in slices]).ravel()
                / component_scales)

    scales = np.maximum(np.abs(amplitudes[free_indices]), .01)

    def objective(z):
        return float(np.dot(z / scales, z / scales))

    bound = np.maximum(.05, .5 * np.abs(amplitudes[free_indices]))
    solved = minimize(objective, np.zeros(len(free_indices)),
                      method="SLSQP",
                      bounds=list(zip(-bound, bound)),
                      constraints=[dict(type="eq", fun=moment_defect)],
                      options=dict(maxiter=300, ftol=1e-12))
    delta = vector(solved.x)
    repaired = field(amplitudes + delta)
    physical_moments = np.array([outer_moments(s, delta) for s in slices])
    cone_scales = []
    for k in (11, 19):
        tau = .5 * 2.**-k
        X = np.full(len(ETA_NODES),
                    inner.p.X_max * (1 + 15 * Y)**2)
        points = inner.from_similarity(X, np.asarray(ETA_NODES), tau)
        rows = [cone_row(repaired, float(point[0]), float(point[2]),
                         tau, eta)
                for point, eta in zip(points, ETA_NODES)]
        cone_scales.append(dict(k=k, rows=rows,
                                passing_eta=[row["eta"] for row in rows
                                             if row["cone_pass"]]))
    report = dict(source="Moment-repaired first-window odd-poloidal axial-cone change",
                  all_knots=all_knots,
                  fixed_indices=fixed_indices,
                  fixed_delta=fixed_delta[fixed_indices].tolist(),
                  optimizer_success=bool(solved.success),
                  optimizer_message=str(solved.message),
                  optimizer_iterations=int(solved.nit),
                  weighted_coefficient_distance=objective(solved.x),
                  max_relative_other_coefficient_change=float(np.max(
                      np.abs(solved.x) / scales)),
                  normalized_moment_max=float(np.max(np.abs(
                      moment_defect(solved.x)))),
                  max_absolute_moment=float(np.max(np.abs(physical_moments))),
                  baseline_max_absolute_moment=float(np.max(np.abs(moments0))),
                  amplitudes=(amplitudes + delta).tolist(),
                  cone_scales=cone_scales,
                  scope="SLSQP restores twelve sampled physical outer moments at three scale knots while holding the two selected odd-poloidal coefficients at half baseline. Cone rows are a sampled physical analogue at five eta values on k=11,19, not a continuous support, full leading-profile construction, or momentum improvement.",
                  accepted=False, scale_recursion_established=False)
    output_name = ("midplane_axial_cone_all_knots_repair.json" if all_knots
                   else "midplane_axial_cone_moment_repair.json")
    (ROOT / output_name).write_bytes(
        (json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(optimizer_success=report["optimizer_success"],
                          normalized_moment_max=report["normalized_moment_max"],
                          passing=[s["passing_eta"] for s in cone_scales]),
                     indent=2))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-knots", action="store_true")
    run(all_knots=parser.parse_args().all_knots)
