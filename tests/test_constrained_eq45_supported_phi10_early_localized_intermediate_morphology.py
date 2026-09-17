import json

import numpy as np

from openai_ns_reconstruction.constrained_eq45_supported_phi10_early_localized_intermediate_morphology import (
    DEFAULT_RESOLUTIONS,
    DEFAULT_TIME,
    audit_early_localized_intermediate_morphology,
)


def test_early_localized_intermediate_morphology_calibration():
    report = audit_early_localized_intermediate_morphology()

    assert report["schema"] == "eq45_supported_phi10_early_localized_intermediate_morphology_v1"
    assert report["resolutions"] == list(DEFAULT_RESOLUTIONS)
    assert report["time"] == DEFAULT_TIME
    np.testing.assert_allclose(report["tau"], 0.5, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(report["coefficients"]["static"], -0.3, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(report["coefficients"]["quadratic"], -0.125, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(report["coefficients"]["affine"], 0.4, rtol=0.0, atol=1.0e-15)

    payload = {
        "finest_row": report["rows"][-1],
        "finest_summary": report["finest_summary"],
    }
    raise AssertionError("CALIBRATION=" + json.dumps(payload, sort_keys=True))
