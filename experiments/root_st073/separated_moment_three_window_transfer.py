"""Check moment closure between the two fitted scale knots."""

import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from radial_continuation import ROOT
from radial_peak_cone import stress_primitive
from separated_moment_modes import SeparatedMomentModes, RADIAL_WINDOWS_THREE


def run():
    inner, fields = build_fields()
    base = fields["two_sided_cone"]
    result = json.loads((ROOT / "separated_moment_three_window.json").read_text())
    field = SeparatedMomentModes(
        base, result["constrained_l4"]["amplitudes"],
        windows=RADIAL_WINDOWS_THREE)
    rows = []
    for k in (11.0, 13.0, 15.0, 17.0, 19.0):
        tau = 0.5 * 2.0**-k
        for eta in (-0.2, 0.2):
            end = inner.from_similarity([inner.p.X_max * 16.0**2],
                                        [eta], tau)[0]
            args = (float(end[0]), float(end[2]), tau)
            current = stress_primitive(field, *args)
            original = stress_primitive(base, *args)
            rows.append(dict(k=k, eta=eta,
                             fitted_outer_stress=current.tolist(),
                             original_outer_stress=original.tolist(),
                             norm_ratio=float(np.linalg.norm(current)
                                              / np.linalg.norm(original))))
    report = dict(
        source="Five-scale direct physical outer-moment transfer of three-window correction",
        rows=rows,
        scope="Two axial sides on five sampled scale slices only; no uniform-in-scale or paper-profile moment theorem.",
        accepted=False, scale_recursion_established=False,
    )
    output = ROOT / "separated_moment_three_window_transfer.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output), rows=[
        {key: row[key] for key in ("k", "eta", "norm_ratio")}
        for row in rows]), indent=2))
    return report


if __name__ == "__main__":
    run()
