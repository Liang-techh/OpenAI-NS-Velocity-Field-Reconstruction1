"""Direct finite-difference validation of the separated exact-moment field."""

import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from adaptive_join_multiscale_fit import sample_points
from joined_field import independent_fd
from radial_continuation import ROOT
from radial_peak_cone import stress_primitive
from separated_moment_modes import SeparatedMomentModes


def run():
    inner, fields = build_fields()
    base = fields["two_sided_cone"]
    fitted = json.loads((ROOT / "separated_moment_fit.json").read_text())
    field = SeparatedMomentModes(base, fitted["exact_moment"]["amplitudes"])
    rows = []
    for k in (11.0, 15.0, 19.0):
        tau = 0.5 * 2.0**-k
        points, _ = sample_points(inner, 16.0, k, (-0.2, 0.0, 0.2), 5)
        h, ht = 0.0005 * np.sqrt(inner.nu * tau), 0.0001 * tau
        residual, divergence = independent_fd(field, points, tau, h, ht)
        edge_rows = []
        hotspot_difference = []
        for eta in (-0.2, 0.2):
            end = inner.from_similarity([inner.p.X_max * 16.0**2],
                                        [eta], tau)[0]
            outer = stress_primitive(field, float(end[0]), float(end[2]), tau)
            old = stress_primitive(base, float(end[0]), float(end[2]), tau)
            edge_rows.append(dict(eta=eta, fitted_outer_stress=outer.tolist(),
                                  original_outer_stress=old.tolist(),
                                  norm_ratio=float(np.linalg.norm(outer)
                                                   / np.linalg.norm(old))))
            hotspot = inner.from_similarity(
                [inner.p.X_max * (1.0 + 15.0 * 0.05)**2], [eta], tau)
            new_u, new_p = field.fields(hotspot, tau)
            old_u, old_p = base.fields(hotspot, tau)
            hotspot_difference.append(float(max(
                np.max(np.abs(new_u - old_u)),
                np.max(np.abs(new_p - old_p)))))
        rows.append(dict(k=k,
                         direct_momentum_max=float(np.max(np.linalg.norm(residual, axis=1))),
                         direct_divergence_max=float(np.max(np.abs(divergence))),
                         edges=edge_rows,
                         hotspot_field_difference_max=max(hotspot_difference)))
    report = dict(
        source="Independent direct FD and radial integration of separated exact-moment field",
        rows=rows,
        scope="Three finite scale slices, two axial sides. Exact field equality at sampled y=.05 hotspots preserves their previously reported local physical cone analogue, but no open support cone, full-space closure, global energy, or PDE acceptance is claimed.",
        accepted=False, scale_recursion_established=False,
    )
    output = ROOT / "separated_moment_validate.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output), rows=rows), indent=2))
    return report


if __name__ == "__main__":
    run()
