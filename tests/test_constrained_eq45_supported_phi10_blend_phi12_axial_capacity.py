import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_blend_phi12_axial_capacity import (
    audit_blend_phi12_axial_capacity,
)


def test_blend_phi12_capacity_is_independent_but_extreme_axial_collar_is_suppressed():
    report = audit_blend_phi12_axial_capacity()

    assert report["numerical_rank"] == 2
    assert np.isfinite(report["condition_number"])
    assert report["condition_number"] < 5.0
    assert report["phi12_novelty_outside_blend_span"] > 0.90
    assert abs(report["response_column_cosine"]) < 0.5
    assert max(report["step_refinement_relative_changes"]) < 1.0e-10

    np.testing.assert_allclose(
        report["singular_values"],
        [0.12329296391791501, 0.03508142769485742],
        rtol=1.0e-8,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        report["combined_response_rms"],
        [0.005449713852512864, 0.0017039714694015845],
        rtol=1.0e-8,
        atol=1.0e-12,
    )

    radial = np.asarray(report["radial_response_rms"], dtype=float)
    axial = np.asarray(report["axial_response_rms"], dtype=float)
    ratios = np.asarray(report["axial_to_radial_response_rms_ratio"], dtype=float)
    assert np.all(radial > 1.0e-3)
    assert np.all(axial < 1.0e-5)
    assert np.all(ratios < 0.01)
    # Phi(1,2) is more axially selective than the blend direction, but the
    # fixed physical support transform still suppresses both directions very
    # strongly at the extreme top/bottom collar probes.
    assert ratios[1] > 5.0 * ratios[0]


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
