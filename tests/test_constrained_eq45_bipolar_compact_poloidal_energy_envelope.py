import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_bipolar_compact_poloidal_energy_envelope import (
    PREVIOUS_CAPACITY_TRIAL,
    _direct_energy,
    _energy_quadratic,
    _largest_symmetric_half_width,
    _quadratic_value,
    _source_field,
    audit_compact_poloidal_energy_envelope,
)


def test_energy_quadratic_matches_direct_compact_velocity_energy():
    field = _source_field()
    coefficients = _energy_quadratic(field, 0.25, 48)
    for coefficient in (-0.05, 0.0, 0.05):
        predicted = _quadratic_value(coefficients, coefficient)
        direct = _direct_energy(field, 0.25, coefficient, 48)
        assert predicted == pytest.approx(direct, rel=2e-12, abs=2e-12)
    assert coefficients["quadratic"] > 0.0


def test_zero_centered_energy_envelope_respects_unchanged_reference_gate():
    field = _source_field()
    coefficients = _energy_quadratic(field, 0.25, 96)
    envelope = _largest_symmetric_half_width(
        coefficients,
        lower=0.999,
        upper=1.001,
        diagnostic_search_guard=field.parent.profile_basis.coefficient_limit,
    )
    assert envelope["exists"] is True
    assert envelope["half_width"] > 0.0
    assert envelope["minimum_energy"] >= 0.999 - 2e-12
    assert envelope["maximum_energy"] <= 1.001 + 2e-12
    assert envelope["guard_limited"] is False


def test_compact_poloidal_energy_envelope_is_target_free_and_non_promoting():
    report = audit_compact_poloidal_energy_envelope(
        quadrature_orders=(48, 96),
        morphology_grid_size=17,
    )
    gate = report["reference_energy_gate"]
    envelope = gate["largest_symmetric_zero_centered_envelope"]
    assert gate["zero_coefficient_gate_satisfied"] is True
    assert abs(gate["baseline_energy"] - 1.0) <= 0.001
    assert envelope["exists"] is True
    assert envelope["half_width"] > 0.0
    assert gate["diagnostic_search_guard_is_materialization_bound"] is False
    assert gate["normalization_roots_select_coefficient"] is False
    assert gate["previous_capacity_trial_abs_coefficient"] == PREVIOUS_CAPACITY_TRIAL

    assert report["validation_time_energy_envelope"]["all_times_satisfied"] is True
    assert report["core_sign_screen"]["all_envelope_endpoint_and_baseline_checks_satisfied"] is True
    morphology = report["vorticity_morphology_at_energy_envelope"]
    for key in (
        "axial_rms_span_over_Zp",
        "radial_rms_span_over_Rp",
        "axial_q90_span_over_Zp",
        "axial_q99_span_over_Zp",
        "outer_axial_enstrophy_fraction_span",
    ):
        assert np.isfinite(morphology[key])
        assert morphology[key] >= 0.0

    assert report["basis_growth"]["new_basis_shapes_added"] == 0
    assert report["truth_boundary"]["velocity_changed"] is False
    assert report["truth_boundary"]["candidate_artifact_changed"] is False
    assert report["truth_boundary"]["compact_poloidal_coefficient_selected"] is False
    assert report["truth_boundary"]["materialization_bound_selected"] is False
    assert report["truth_boundary"]["held_out_pde_residual_evaluated"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False
    assert report["truth_boundary"]["pde_validated"] is False


def test_malformed_energy_envelope_inputs_fail_closed():
    field = _source_field()
    coefficients = _energy_quadratic(field, 0.25, 32)
    with pytest.raises(ValueError):
        _largest_symmetric_half_width(
            coefficients, lower=1.001, upper=0.999, diagnostic_search_guard=1.0
        )
    with pytest.raises(ValueError):
        audit_compact_poloidal_energy_envelope(quadrature_orders=(96, 48), morphology_grid_size=17)
