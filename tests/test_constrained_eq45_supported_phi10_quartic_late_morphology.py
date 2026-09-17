import json

import numpy as np

from openai_ns_reconstruction.constrained_eq45_supported_phi10_quartic_late_morphology import (
    DEFAULT_RESOLUTIONS,
    DEFAULT_TIME,
    audit_quartic_late_offkeyframe_morphology,
)


def test_quartic_late_offkeyframe_morphology_calibration():
    report = audit_quartic_late_offkeyframe_morphology()

    assert report["schema"] == "eq45_supported_phi10_quartic_late_offkeyframe_morphology_v1"
    assert report["resolutions"] == list(DEFAULT_RESOLUTIONS)
    assert report["time"] == DEFAULT_TIME
    np.testing.assert_allclose(report["tau"], 0.75, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(report["coefficients"]["static"], -0.3, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(report["coefficients"]["cubic"], -0.321875, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(
        report["quartic_nullspace_coefficient"],
        -0.2831460674157303,
        rtol=0.0,
        atol=1.0e-15,
    )

    fine = report["rows"][-1]
    middle = report["rows"][-2]
    payload = {
        "coefficients": report["coefficients"],
        "fine_baseline": fine["baseline"],
        "fine_cubic": fine["cubic"],
        "fine_quartic": fine["quartic"],
        "middle_quartic": middle["quartic"],
        "summary": report["finest_summary"],
    }
    raise AssertionError("AGENT4_QUARTIC_MORPHOLOGY_CALIBRATION=" + json.dumps(payload, sort_keys=True))
