import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_blend_phi12_axial_capacity import (
    audit_blend_phi12_axial_capacity,
)


def test_blend_phi12_axial_capacity_calibration_receipt():
    report = audit_blend_phi12_axial_capacity()
    # Calibration sentinel: expose deterministic hosted metrics before replacing
    # this with fail-closed scientific regressions on the exact same audit.
    pytest.fail(json.dumps(report, sort_keys=True))


def test_blend_phi12_axial_capacity_truth_boundary_and_guards():
    report = audit_blend_phi12_axial_capacity()
    truth = report["truth_boundary"]
    for key in (
        "velocity_changed",
        "canonical_velocity_changed",
        "production_coefficients_changed",
        "new_basis_added",
        "blend_weight_selected",
        "phi12_coefficient_selected",
        "force_or_pressure_fitted",
        "pde_objective_used",
        "public_image_fitted",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False

    with pytest.raises(ValueError, match="strictly inside"):
        audit_blend_phi12_axial_capacity(blend_weight=0.0)
    with pytest.raises(ValueError, match="delivery interval"):
        audit_blend_phi12_axial_capacity(time=0.8)
    with pytest.raises(ValueError, match="equal length"):
        audit_blend_phi12_axial_capacity(blend_steps=(0.1, 0.05), phi12_steps=(0.04,))
    with pytest.raises(ValueError, match=r"\[0,1\]"):
        audit_blend_phi12_axial_capacity(blend_weight=0.05, blend_steps=(0.1, 0.05))

    assert np.all(np.isfinite(report["singular_values"]))
