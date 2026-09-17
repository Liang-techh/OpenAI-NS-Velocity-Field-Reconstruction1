import json

import numpy as np

from openai_ns_reconstruction.constrained_eq45_supported_phi10_cubic_late_morphology import (
    DEFAULT_RESOLUTIONS,
    DEFAULT_TIME,
    audit_cubic_late_offkeyframe_morphology,
)


def test_cubic_late_offkeyframe_morphology_calibration():
    report = audit_cubic_late_offkeyframe_morphology()

    assert report["schema"] == "eq45_supported_phi10_cubic_late_offkeyframe_morphology_v1"
    assert report["resolutions"] == list(DEFAULT_RESOLUTIONS)
    assert report["time"] == DEFAULT_TIME
    np.testing.assert_allclose(report["tau"], 0.75, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(report["coefficients"]["static"], -0.3, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(report["coefficients"]["quadratic"], -0.16875, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(report["coefficients"]["cubic"], -0.321875, rtol=0.0, atol=1.0e-15)

    fine = report["rows"][-1]
    summary = report["finest_summary"]
    payload = {
        "fine": {
            name: {
                key: fine[name][key]
                for key in (
                    "radial_q99",
                    "axial_q99",
                    "aspect_q99",
                    "vorticity2_weighted_radial_rms",
                    "vorticity2_weighted_axial_rms",
                    "vorticity2_weighted_aspect",
                    "radial_outer_vorticity2_fraction",
                    "whole_grid_collar_vorticity2_fraction",
                    "superlevel_collar_voxel_fraction",
                )
            }
            for name in ("baseline", "quadratic", "cubic")
        },
        "summary": summary,
    }
    raise AssertionError("CR008_CUBIC_LATE_CALIBRATION=" + json.dumps(payload, sort_keys=True))
