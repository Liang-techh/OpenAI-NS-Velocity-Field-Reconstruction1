import json

import numpy as np

from openai_ns_reconstruction.bipolar_spatial_channel_redistribution import REGION_NAMES, run


def _calibration_summary(report):
    return {
        "phi_delta": report["phi_coefficient_max_abs_delta"],
        "swirl_delta": report["swirl_coefficient_max_abs_delta"],
        "fine": report["finest_comparison"],
        "center": report["central_signed_components"],
        "stability": report["medium_to_fine_stability"],
        "fine_channels": [
            {
                "time": row["time"],
                "odd_fraction": row["odd3"]["channel_fraction"],
                "capped_fraction": row["capped"]["channel_fraction"],
                "odd_region_fraction": row["odd3"]["region_fraction_of_total"],
                "capped_region_fraction": row["capped"]["region_fraction_of_total"],
                "odd_swirl_concentration": row["odd3"]["swirl_concentration"],
                "capped_swirl_concentration": row["capped"]["swirl_concentration"],
            }
            for row in report["rows"]
            if row["resolution"] == report["resolutions"][-1]
        ],
    }


def test_bipolar_capped_spatial_channel_audit_calibration():
    report = run()
    assert report["resolutions"] == [24, 32, 40]
    assert report["times"] == [0.3125, 0.5, 0.75]
    assert report["odd3_sha256"] != report["capped_sha256"]
    assert report["phi_coefficient_max_abs_delta"] < 1e-10
    assert report["swirl_coefficient_max_abs_delta"] > 1.0
    for row in report["rows"]:
        for candidate in ("odd3", "capped"):
            metrics = row[candidate]
            assert np.isfinite(metrics["energy"]["total"])
            assert metrics["energy"]["total"] > 0.0
            assert abs(sum(metrics["channel_fraction"].values()) - 1.0) < 1e-12
            assert abs(sum(metrics["region_fraction_of_total"].values()) - 1.0) < 1e-12
            assert set(metrics["regions"]) == set(REGION_NAMES)
    # Calibration head: expose deterministic scientific values once, then replace
    # this fail with pinned fail-closed regressions on the final head.
    raise AssertionError(json.dumps(_calibration_summary(report), sort_keys=True))
