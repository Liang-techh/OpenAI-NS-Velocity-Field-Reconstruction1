"""Classify which stress component must change to enter the cone.

Uses saved physical-cone geometry from delayed_pressure_cone_window.json;
no new field sampling is needed. This is a pointwise necessary-direction
diagnostic, not a realizable stress or wave construction.
"""

import json

import numpy as np

from radial_continuation import ROOT


def line_interval(directions, target, margin, component):
    lower, upper = -np.inf, np.inf
    for direction in directions:
        rhs = -margin-direction@target
        slope = direction[component]
        if slope > 1e-14:
            upper = min(upper, rhs/slope)
        elif slope < -1e-14:
            lower = max(lower, rhs/slope)
        elif rhs < 0:
            return None
    return (lower, upper) if lower <= upper else None


def minimum_abs(interval):
    if interval is None:
        return None
    return float(np.clip(0., *interval))


def closest_two_component_target(directions, target, margin):
    matrix = np.stack(directions)
    candidates = []
    if np.max(matrix@target+margin) <= 1e-10:
        candidates.append(target)
    for direction in directions:
        projected = target-((direction@target+margin)
                            /(direction@direction))*direction
        if np.max(matrix@projected+margin) <= 1e-8:
            candidates.append(projected)
    if abs(np.linalg.det(matrix)) > 1e-12:
        candidates.append(np.linalg.solve(matrix,
                                          -margin*np.ones(2)))
    if not candidates:
        return None
    return min(candidates, key=lambda item: np.linalg.norm(item-target))


def run():
    source = json.loads((ROOT/'delayed_pressure_cone_window.json').read_text())
    ratio = source['cone_ratio_margin']
    margin = source['stress_margin']
    rows = []
    for row in source['raw_geometry']:
        target = np.asarray(row['target'])
        N, K = np.asarray(row['N']), np.asarray(row['K'])
        directions = [ratio*N+sign*row['c']*K for sign in (-1., 1.)]
        relaxed_directions = [.999*N+sign*row['c']*K
                              for sign in (-1., 1.)]
        theta = line_interval(directions, target, margin, 0)
        axial = line_interval(directions, target, margin, 1)
        relaxed_theta = line_interval(relaxed_directions, target, 0., 0)
        relaxed_axial = line_interval(relaxed_directions, target, 0., 1)
        closest = closest_two_component_target(directions, target, margin)
        delta = closest-target if closest is not None else None
        rows.append(dict(X=row['X'], eta=row['eta'],
                         target=target.tolist(),
                         free_theta_feasible=theta is not None,
                         free_theta_relaxed_0999=relaxed_theta is not None,
                         min_abs_theta_change=minimum_abs(theta),
                         free_axial_feasible=axial is not None,
                         free_axial_relaxed_0999=relaxed_axial is not None,
                         min_abs_axial_change=minimum_abs(axial),
                         closest_two_component_delta=(
                             delta.tolist() if delta is not None else None),
                         closest_two_component_norm=(
                             float(np.linalg.norm(delta))
                             if delta is not None else None)))
    report = dict(source='delayed_pressure_cone_window.json',
                  ratio_margin=ratio, stress_margin=margin,
                  rows=rows,
                  scope='Pointwise geometry with one stress component '
                        'free at a time. Does not show that any divergence '
                        'stress or supported wave can realize the change.',
                  accepted=False)
    (ROOT/'delayed_pressure_cone_geometry.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(rows), flush=True)


if __name__ == '__main__':
    run()
