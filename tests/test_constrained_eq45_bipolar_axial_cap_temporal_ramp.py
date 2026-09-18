import math

import pytest

from openai_ns_reconstruction.constrained_eq45_bipolar_axial_cap_temporal_ramp import (
    _linear_ramp_coefficient,
    audit_axial_cap_temporal_ramp_capacity,
)


def test_linear_ramp_replays_reference_and_reaches_terminal_value():
    kwargs = dict(start_time=0.25, end_time=0.75, terminal_coefficient=0.25)
    assert _linear_ramp_coefficient(0.25, **kwargs) == 0.0
    assert _linear_ramp_coefficient(0.50, **kwargs) == pytest.approx(0.125)
    assert _linear_ramp_coefficient(0.75, **kwargs) == pytest.approx(0.25)


def test_axial_cap_temporal_ramp_preserves_registered_guards_and_reports_morphology():
    report = audit_axial_cap_temporal_ramp_capacity(
        quadrature_orders=(24, 48), morphology_grid_size=17
    )

    growth = report["representation_increment"]
    assert growth["new_spatial_basis_shapes_added"] == 0
    assert growth["diagnostic_temporal_shapes_added"] == 1
    assert growth["new_production_fit_parameters_added"] == 0
    assert growth["production_coefficient_selected"] is False
    assert report["both_signs_reference_time_parent_replay"] is True
    assert report["both_signs_broad_energy_range_passed"] is True
    assert report["both_signs_core_signs_passed"] is True
    assert report["truth_boundary"]["held_out_pde_residual_evaluated"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False
    assert report["truth_boundary"]["pde_validated"] is False

    for name, sign in (("negative_schedule", -1), ("positive_schedule", 1)):
        row = report[name]
        assert row["sign"] == sign
        assert row["reference_time_coefficient"] == 0.0
        assert row["reference_time_morphology_replay_max_abs_error"] == 0.0
        assert row["terminal_coefficient"] == pytest.approx(sign * 0.25)
        assert row["coefficient_time_derivative"] == pytest.approx(sign * 0.5)
        coefficients = [entry["coefficient"] for entry in row["morphology_rows"]]
        assert coefficients == pytest.approx([0.0, sign * 0.125, sign * 0.25])
        for entry in row["morphology_rows"]:
            for value in entry["trial_morphology"].values():
                if isinstance(value, (int, float)):
                    assert math.isfinite(value)


def test_axial_cap_temporal_ramp_rejects_invalid_inputs():
    with pytest.raises(ValueError, match="terminal_magnitude"):
        audit_axial_cap_temporal_ramp_capacity(terminal_magnitude=0.0)
    with pytest.raises(ValueError, match="morphology_grid_size"):
        audit_axial_cap_temporal_ramp_capacity(morphology_grid_size=24)
    with pytest.raises(ValueError, match="morphology_times"):
        audit_axial_cap_temporal_ramp_capacity(morphology_times=(0.5, 0.4))
    with pytest.raises(ValueError, match="outside"):
        _linear_ramp_coefficient(
            0.8, start_time=0.25, end_time=0.75, terminal_coefficient=0.25
        )
