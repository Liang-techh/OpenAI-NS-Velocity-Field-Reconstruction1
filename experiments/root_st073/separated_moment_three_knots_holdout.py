"""Compare exact-moment and L4-constrained three-knot time holdouts."""

import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from adaptive_join_multiscale_fit import sample_points
from joined_field import independent_fd
from radial_continuation import ROOT
from radial_peak_cone import stress_primitive
from separated_moment_modes import SeparatedMomentModes, RADIAL_WINDOWS_THREE


def run():
    inner, fields = build_fields()
    base = fields["two_sided_cone"]
    source = json.loads((ROOT / "separated_moment_three_knots.json").read_text())
    rows = []
    for candidate in ("exact", "constrained_l4"):
        field = SeparatedMomentModes(
            base, source[candidate]["amplitudes"],
            windows=RADIAL_WINDOWS_THREE, knots=(11.0, 15.0, 19.0))
        for k in (13.0, 17.0):
            tau = 0.5 * 2.0**-k
            points, _ = sample_points(inner, 16.0, k, (-0.2, 0.0, 0.2), 5)
            residual, divergence = independent_fd(
                field, points, tau,
                0.0005 * np.sqrt(inner.nu * tau), 0.0001 * tau)
            ratios = []
            for eta in (-0.2, 0.2):
                end = inner.from_similarity([inner.p.X_max * 16.0**2],
                                            [eta], tau)[0]
                args = (float(end[0]), float(end[2]), tau)
                current = stress_primitive(field, *args)
                old = stress_primitive(base, *args)
                ratios.append(float(np.linalg.norm(current)
                                    / np.linalg.norm(old)))
            rows.append(dict(candidate=candidate, k=k,
                             momentum_max=float(np.max(np.linalg.norm(residual, axis=1))),
                             divergence_max=float(np.max(np.abs(divergence))),
                             outer_moment_ratios=ratios))
    report = dict(
        source="Off-knot comparison for three-knot separated moment field",
        rows=rows,
        scope="Two time holdouts and two axial slices only; no continuous scale or PDE acceptance.",
        accepted=False,
    )
    output = ROOT / "separated_moment_three_knots_holdout.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output), rows=rows), indent=2))
    return report


if __name__ == "__main__":
    run()
