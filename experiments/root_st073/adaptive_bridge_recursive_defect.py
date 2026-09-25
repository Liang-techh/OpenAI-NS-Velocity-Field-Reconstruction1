"""Measure the dyadic transfer of the complete bridge momentum defect.

This is a finite-grid diagnostic of a fixed bridge ansatz, not a proof of
asymptotic recursion or a volume norm.
"""

import json

import numpy as np

from adaptive_core_join_screen import AdaptiveRadialAdapter
from adaptive_join_multiscale_fit import CachedAdaptiveBridge, sample_points
from affine_momentum import jets, momentum
from heat_exterior import physical
from joined_field import JoinedField
from radial_continuation import ROOT
from wide_modes import WideJointModes


SCALES = (11.0, 13.0, 15.0, 17.0, 19.0)


def build_fields(k_max=20):
    inner = AdaptiveRadialAdapter(k_max=k_max)
    tau0 = 0.5 * 2.0**-6
    point = inner.from_similarity([1.0 / 64.0], [0.0], tau0)
    amplitude = float(inner.evaluate(point, tau0)["velocity"][0, 1]
                      / physical(point, tau0, c=1.0)["velocity"][0, 1])
    joined = JoinedField(inner=inner, join_X=1.0 / 64.0,
                         heat_amplitude=amplitude, outer_ratio=16.0)
    cached = CachedAdaptiveBridge(joined)
    shared = json.loads((ROOT / "adaptive_join_multiscale_fit.json").read_text())
    biside = json.loads((ROOT / "adaptive_bridge_cone_biside_fit.json").read_text())
    selected = next(row for row in biside["candidates"]
                    if row["cone_weight"] == biside["selected_cone_weight"])
    return inner, {
        "unmodified": WideJointModes(cached, np.zeros(24)),
        "shared_24_mode": WideJointModes(cached, shared["amplitudes"]),
        "two_sided_cone": WideJointModes(cached, selected["amplitudes"]),
    }


def run():
    inner, fields = build_fields()
    results = {}
    for name, field in fields.items():
        rows = []
        reference = None
        for k in SCALES:
            points, tau = sample_points(inner, 16.0, k,
                                        (-0.2, 0.0, 0.2), 5)
            h, ht = 0.0005 * np.sqrt(inner.nu * tau), 0.0001 * tau
            residual = momentum(jets(field, points, tau, h, ht))
            normalized = tau**1.5 * residual
            if reference is None:
                reference = normalized.copy()
            norms = np.linalg.norm(residual, axis=1)
            rows.append(dict(
                k=k, tau=tau,
                maximum=float(np.max(norms)),
                normalized_maximum=float(np.max(np.linalg.norm(normalized, axis=1))),
                normalized_vector_l2=float(np.linalg.norm(normalized)),
                transfer_drift_to_k11=float(
                    np.linalg.norm(normalized - reference)
                    / np.linalg.norm(reference)),
                hotspot_index=int(np.argmax(norms)),
            ))
        ratios = [rows[i + 1]["normalized_vector_l2"]
                  / rows[i]["normalized_vector_l2"]
                  for i in range(len(rows) - 1)]
        results[name] = dict(
            rows=rows,
            normalized_two_halving_ratios=ratios,
            maximum_growth_exponent_per_halving=float(
                np.log2(rows[-1]["maximum"] / rows[0]["maximum"])
                / (SCALES[-1] - SCALES[0])),
        )
    report = dict(
        source="Fixed-coefficient adaptive core plus ratio-16 bridge",
        scales=list(SCALES), results=results,
        scope="Same 15 similarity-grid points on five scales; complete Cartesian finite-difference momentum, multiplied by tau^1.5. This is not a continuum, volume-L2, recursive-fixed-point, or global-PDE certificate.",
        accepted=False, scale_recursion_established=False,
    )
    output = ROOT / "adaptive_bridge_recursive_defect.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps(dict(output=str(output),
                          summary={name: dict(
                              normalized_two_halving_ratios=data["normalized_two_halving_ratios"],
                              growth_exponent=data["maximum_growth_exponent_per_halving"],
                              first=data["rows"][0], last=data["rows"][-1])
                              for name, data in results.items()}), indent=2))
    return report


if __name__ == "__main__":
    run()
