"""Audit physical tangential residual moments at the outer bridge edge."""

import json

import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from radial_peak_cone import stress_primitive
from radial_continuation import ROOT


def run():
    inner, fields = build_fields()
    rows = []
    for k in (11.0, 15.0, 19.0):
        tau = 0.5 * 2.0**-k
        for eta in (-0.2, 0.2):
            X = inner.p.X_max * np.array([(1.0 + 15.0 * y)**2
                                           for y in (0.05, 1.0)])
            points = inner.from_similarity(X, [eta, eta], tau)
            hotspot_r, outer_r = map(float, points[:, 0])
            z = float(points[0, 2])
            field_rows = {}
            for name, field in fields.items():
                local = stress_primitive(field, hotspot_r, z, tau, order=12)
                outer = stress_primitive(field, outer_r, z, tau, order=12)
                field_rows[name] = dict(
                    local_stress=local.tolist(),
                    outer_stress=outer.tolist(),
                    outer_to_local_norm=float(np.linalg.norm(outer)
                                              / max(np.linalg.norm(local), 1e-30)),
                )
            rows.append(dict(k=k, eta=eta, tau=tau, hotspot_radius=hotspot_r,
                             outer_radius=outer_r, z=z, fields=field_rows))
    report = dict(
        source="Physical full-residual radial stress primitive of adaptive bridge",
        rows=rows,
        scope="Axis-to-outer-edge weighted radial integrals of full tangential Cartesian residual at six time/axial slices. Nonzero outer primitive obstructs a compact stress lift of this physical residual. This is not the paper's leading normalized moment system, and no continuum or PDE acceptance is claimed.",
        accepted=False, paper_radial_moments_satisfied=False,
    )
    output = ROOT / "adaptive_bridge_outer_moments.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output), summary=[
        dict(k=row["k"], eta=row["eta"], fields={
            name: value["outer_to_local_norm"]
            for name, value in row["fields"].items()})
        for row in rows]), indent=2))
    return report


if __name__ == "__main__":
    run()
