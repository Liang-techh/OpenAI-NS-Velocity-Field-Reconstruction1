"""Check the repaired mean on every integer scale between the fit knots."""

import argparse
import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from midplane_axial_cone_knob import KNOTS, Y, cone_row
from radial_continuation import ROOT
from radial_peak_cone import operator
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


ETA_NODES = (-.05, -.0125, .025, .05)


def run(all_knots=False):
    inner, fields = build_fields()
    base = fields["two_sided_cone"]
    original_data = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    repair_name = ("midplane_axial_cone_all_knots_repair.json" if all_knots
                   else "midplane_axial_cone_moment_repair.json")
    repaired_data = json.loads((ROOT / repair_name).read_text())
    original = SeparatedMomentModes(
        base, np.asarray(original_data["constrained_l4"]["amplitudes"]),
        windows=RADIAL_WINDOWS_THREE, knots=KNOTS)
    repaired = SeparatedMomentModes(
        base, np.asarray(repaired_data["amplitudes"]),
        windows=RADIAL_WINDOWS_THREE, knots=KNOTS)
    scales = []
    for k in range(11, 20):
        tau = .5 * 2.**-k
        X = np.full(len(ETA_NODES), inner.p.X_max * (1 + 15 * Y)**2)
        coordinates = inner.from_similarity(X, np.asarray(ETA_NODES), tau)
        points = np.column_stack((coordinates[:, 0],
                                  np.zeros(len(ETA_NODES)),
                                  coordinates[:, 2]))
        before = np.linalg.norm(operator(original, points, tau)[2], axis=1)
        after = np.linalg.norm(operator(repaired, points, tau)[2], axis=1)
        cone = cone_row(repaired, float(points[-1, 0]),
                        float(points[-1, 2]), tau, ETA_NODES[-1])
        row = dict(k=k, original_momentum_norm=before.tolist(),
                   repaired_momentum_norm=after.tolist(),
                   repaired_to_original_max=float(np.max(after) / np.max(before)),
                   cone_at_positive_edge=cone)
        scales.append(row)
        print(f"k={k}: residual ratio={row['repaired_to_original_max']:.6g}, "
              f"cone ratio={cone.get('cone_ratio', float('nan')):.6g}",
              flush=True)
    output = dict(source="Mean-only integer-scale holdout between moment knots",
                  all_knots=all_knots,
                  y=Y, eta_nodes=ETA_NODES, scales=scales,
                  scope="Physical mean-field momentum on four local axial nodes "
                        "and one cone row per scale; no wave, interior-time "
                        "budget, spatial maximum, or recursive contraction claim.")
    output_name = ("midplane_axial_cone_all_knots_scale_scan.json" if all_knots
                   else "midplane_axial_cone_scale_scan.json")
    (ROOT / output_name).write_bytes(
        (json.dumps(output, indent=2) + "\n").encode())
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-knots", action="store_true")
    run(all_knots=parser.parse_args().all_knots)
