"""Screen the actual rectangular wave support for a connected stress cone."""

import json

import numpy as np

from midplane_axial_cone_knob import cone_row
from midplane_wider_cone_slope_projection import build_case
from radial_continuation import ROOT


def run():
    angles = np.arange(16) * 2 * np.pi / 16
    offsets = (-.8, 0., .8)
    scales = []
    for k in (11, 19):
        mean, wave, _, _ = build_case(k, angles)
        rows = []
        for xi in offsets:
            for ez in offsets:
                r = wave.radius + xi * wave.radial_halfwidth
                z = wave.zcenter + ez * wave.axial_halfwidth
                cone = cone_row(mean, float(r), float(z), wave.tau0, ez)
                rows.append(dict(radial_offset=xi, axial_offset=ez,
                                 cone_pass=cone["cone_pass"],
                                 cone_ratio=cone.get("cone_ratio"),
                                 lambda_squared=cone.get("lambda_squared"),
                                 target_dot_N=cone.get("target_dot_N")))
        scales.append(dict(k=k, rows=rows,
                           passing=sum(row["cone_pass"] for row in rows)))
        print(f"k={k}: {scales[-1]['passing']}/{len(rows)} cone samples pass",
              flush=True)
    report = dict(source="Wider-wave physical rectangular support cone screen",
                  offsets=offsets, scales=scales,
                  scope="Nine physical radial-axial samples per scale at pulse "
                        "center; no continuous spatial or temporal cone proof.",
                  full_support_cone_established=False)
    (ROOT / "midplane_wider_cone_spatial_screen.json").write_bytes(
        (json.dumps(report, indent=2) + "\n").encode())
    return report


if __name__ == "__main__":
    run()
