import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_c2_compact_early_morphology import (
    audit_c2_compact_early_morphology,
)


def test_c2_compact_early_morphology_calibration_receipt():
    report = audit_c2_compact_early_morphology()
    payload = {
        "time": report["time"],
        "tau": report["tau"],
        "compact_return_time": report["compact_return_time"],
        "coefficients": report["coefficients"],
        "rows": report["rows"],
        "finest_summary": report["finest_summary"],
    }
    raise AssertionError("CALIBRATION_RECEIPT=" + json.dumps(payload, sort_keys=True))


def test_c2_compact_early_morphology_fails_closed_outside_active_window():
    with pytest.raises(ValueError, match="active compact C2 transition"):
        audit_c2_compact_early_morphology(time=0.4375)
    with pytest.raises(ValueError, match="strictly increasing"):
        audit_c2_compact_early_morphology(resolutions=(49, 49, 81))
    with pytest.raises(ValueError, match="superlevel_fraction"):
        audit_c2_compact_early_morphology(superlevel_fraction=1.0)


def test_c2_compact_early_morphology_truth_boundary_stays_false():
    report = audit_c2_compact_early_morphology(resolutions=(17, 21, 25))
    truth = report["truth_boundary"]
    assert truth["visualization_candidate_only"] is True
    for key in (
        "canonical_velocity_changed",
        "diagnostic_velocity_changed",
        "compact_temporal_shape_promoted",
        "public_image_fitted",
        "physical_support_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False
    assert np.isfinite(report["finest_summary"]["compact_baseline_vorticity2_weighted_radial_rms_ratio"])
