import pytest

from openai_ns_reconstruction.constrained_eq45_supported_collar_vorticity import (
    governed_supported_seed,
)
from openai_ns_reconstruction.constrained_eq45_supported_early_radial_capacity import (
    audit_supported_early_radial_capacity,
)


def test_supported_early_radial_capacity_calibration():
    report = audit_supported_early_radial_capacity(governed_supported_seed())
    pytest.fail(
        repr(
            {
                "baseline_feature": report["baseline_feature"],
                "response_matrix": report["response_matrix"],
                "singular_values": report["singular_values"],
                "rank": report["numerical_rank"],
                "condition": report["condition_number"],
                "refinement": report["response_refinement_relative_change"],
                "projection": report["early_only_projection"],
                "trial": report["bounded_local_trial"],
            }
        )
    )


def test_supported_early_radial_capacity_rejects_bad_inputs():
    child = governed_supported_seed()
    with pytest.raises(ValueError):
        audit_supported_early_radial_capacity(child, coefficient_steps=(0.02, 0.04))
    with pytest.raises(ValueError):
        audit_supported_early_radial_capacity(child, times=(0.25, 0.50, 0.80))
    with pytest.raises(ValueError):
        audit_supported_early_radial_capacity(child, resolution=19)


def test_supported_early_radial_capacity_does_not_mutate_or_promote_candidate():
    child = governed_supported_seed()
    before = child.sha256
    report = audit_supported_early_radial_capacity(child, resolution=25)

    assert child.sha256 == before
    assert report["parameter_count_added"] == 0
    assert report["existing_parameters_audited"] == 2
    assert report["interpretation"]["new_spatial_basis_added"] is False
    assert report["interpretation"]["public_visual_target_fitted"] is False

    truth = report["truth_boundary"]
    assert truth["velocity_changed"] is False
    assert truth["production_coefficients_changed"] is False
    assert truth["new_basis_added"] is False
    assert truth["forcing_or_pressure_refit"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
    assert truth["blowup_proved"] is False
