import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_bipolar_axial_cap_energy_envelope import (
    PREVIOUS_CAPACITY_TRIAL,
    TASK_ID,
    _direct_energy,
    _energy,
    _energy_quadratic,
    _source_field,
    audit_axial_cap_energy_envelope,
)


def test_axial_cap_energy_quadratic_matches_direct_replay():
    field = _source_field()
    coefficients = _energy_quadratic(field, 0.25, 72)
    for value in (-0.25, -0.07, 0.0, 0.09, 0.25):
        assert _energy(coefficients, value) == pytest.approx(
            _direct_energy(field, 0.25, value, 72), rel=2e-12, abs=2e-12
        )


def test_axial_cap_fixed_parent_envelope_is_truth_bounded_and_replayable():
    report = audit_axial_cap_energy_envelope(
        quadrature_orders=(48, 96), morphology_grid_size=17
    )
    assert report["task_id"] == TASK_ID
    envelope = report["zero_centered_fixed_parent_energy_envelope"]
    assert envelope["exists"] is True
    assert 0.0 < envelope["half_width"] <= 4.0
    assert envelope["minimum_energy"] >= 0.999 - 5e-13
    assert envelope["maximum_energy"] <= 1.001 + 5e-13
    assert report["validation_time_energy_range_satisfied_over_envelope"] is True
    assert np.isfinite(report["quadrature_max_relative_coefficient_change"])
    assert report["quadrature_max_relative_coefficient_change"] < 0.2
    assert report["previous_capacity_trial_abs_coefficient"] == PREVIOUS_CAPACITY_TRIAL
    assert report["capacity_trial_to_energy_envelope_amplitude_ratio"] > 0.0
    assert max(row["absolute_difference"] for row in report["direct_energy_replay_checks"]) < 2e-11

    morphology = report["morphology"]
    for key in (
        "baseline_q90_over_Zp",
        "energy_envelope_max_q90_over_Zp",
        "capacity_trial_max_q90_over_Zp",
        "baseline_q99_over_Zp",
        "energy_envelope_max_q99_over_Zp",
        "capacity_trial_max_q99_over_Zp",
    ):
        assert np.isfinite(morphology[key])
        assert morphology[key] >= 0.0

    truth = report["truth_boundary"]
    assert truth["velocity_changed"] is False
    assert truth["axial_cap_poloidal_coefficient_selected"] is False
    assert truth["materialization_bound_selected"] is False
    assert truth["held_out_pde_residual_evaluated"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_axial_cap_energy_audit_rejects_bad_resolution_contracts():
    with pytest.raises(ValueError):
        audit_axial_cap_energy_envelope(quadrature_orders=(96, 48), morphology_grid_size=17)
    with pytest.raises(ValueError):
        audit_axial_cap_energy_envelope(quadrature_orders=(8, 16), morphology_grid_size=17)
    with pytest.raises(ValueError):
        audit_axial_cap_energy_envelope(quadrature_orders=(48, 96), morphology_grid_size=18)
