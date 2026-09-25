"""Select positive mode pairs from actual exact-curl center covariances."""

import itertools
import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from curl_wave_prototype import LocalizedCurlWave, cylindrical_residual
from radial_continuation import ROOT


def support(inner, source, nu):
    tau = float(source["tau"])
    X = inner.p.X_max * (1 + 15 * 0.325)**2
    ri = float(inner.from_similarity([inner.p.X_max], [-0.0125], tau)[0, 0])
    radial_halfwidth = 0.035 * 15 * ri
    below = inner.from_similarity([X], [-0.0375], tau)[0, 2]
    above = inner.from_similarity([X], [0.0125], tau)[0, 2]
    zcenter = float(source["point"][2])
    axial_halfwidth = 0.99 * min(zcenter - below, above - zcenter)
    time_halfwidth = 0.75 * axial_halfwidth**2 / nu
    return dict(radial_halfwidth=radial_halfwidth,
                axial_halfwidth=axial_halfwidth,
                time_halfwidth=time_halfwidth,
                min_tau=0.5 * 2.0**-20)


def sampled_covariance(wave, points, tau):
    values, _ = wave.fields(points, tau)
    cyl = cylindrical_residual(values, points)
    return np.mean(cyl[:, 0, None] * cyl[:, 1:], axis=0)


def run():
    inner, _ = build_fields()
    scales = []
    for k in (11, 19):
        source = json.loads((ROOT / "compact_potential" /
                             f"midplane_wave_source_k{k}.json").read_text())
        source["nu"] = inner.nu
        tau = float(source["tau"])
        args = support(inner, source, inner.nu)
        angles = np.arange(16) * 2 * np.pi / 16
        radius, _, z = source["point"]
        points = np.array([(radius * np.cos(a), radius * np.sin(a), z)
                           for a in angles])
        pulses = source["kelvin_pulses"]
        columns = []
        mode_norms = []
        for i, pulse in enumerate(pulses):
            partner = next(j for j, other in enumerate(pulses)
                           if other["mode"] != pulse["mode"])
            wave = LocalizedCurlWave(source, pulse_indices=(i, partner),
                                     **args)
            wave.weights = np.array([1.0, 0.0])
            columns.append(sampled_covariance(wave, points, tau))
            mode_norms.append(float(np.dot(wave.waves[0]["normal"],
                                            wave.waves[0]["normal"])))
        target = np.asarray(source["local_tangential_stress_primitive"])
        candidates = []
        current = tuple(LocalizedCurlWave(source, **args).pulse_indices)
        current_unconstrained = None
        current_condition = None
        for i, j in itertools.combinations(range(len(pulses)), 2):
            if pulses[i]["mode"] == pulses[j]["mode"]:
                continue
            matrix = np.column_stack((columns[i], columns[j]))
            condition = float(np.linalg.cond(matrix))
            if (i, j) == current:
                current_condition = condition
                current_unconstrained = np.linalg.solve(
                    matrix, target).tolist()
            if condition > 1e8:
                continue
            weights = np.linalg.solve(matrix, target)
            if np.all(weights > 0):
                cost = float(inner.nu * (weights[0] * mode_norms[i]
                                          + weights[1] * mode_norms[j]))
                candidates.append(dict(indices=[i, j],
                                       modes=[pulses[i]["mode"],
                                              pulses[j]["mode"]],
                                       positive_weights=weights.tolist(),
                                       carrier_viscous_cost=cost,
                                       matrix_condition=float(
                                           np.linalg.cond(matrix))))
        candidates.sort(key=lambda row: row["carrier_viscous_cost"])
        scales.append(dict(k=k, tau=tau, support=args,
                           source_ideal_nnls_indices=[
                               item["index"] for item in
                               source["covariance_nnls"]["selected"]],
                           ideal_selected_indices=list(current),
                           ideal_pair_unconstrained_physical_weights=(
                               current_unconstrained),
                           ideal_pair_physical_matrix_condition=(
                               current_condition),
                           positive_pair_count=len(candidates),
                           selected=candidates[0] if candidates else None,
                           top_pairs=candidates[:10]))
    report = dict(source="Positive covariance pairs for actual center exact-curl fields",
                  scales=scales,
                  scope="Actual exact-curl angular covariance at one point/time per scale, 16 angles, fixed proposed support. Positive two-mode pairs ranked by weighted carrier viscous cost. No open-cone, transported amplitude, full residual, or scale recursion claim.",
                  accepted=False, scale_recursion_established=False)
    output = ROOT / "midplane_physical_covariance_pairs.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output), summary=[
        {key: item[key] for key in (
            "k", "ideal_pair_unconstrained_physical_weights",
            "positive_pair_count", "selected")}
        for item in scales]), indent=2))
    return report


if __name__ == "__main__":
    run()
