from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_partial_nonlinear_radial_force import (
    FINE_PAIR_RELATIVE_STABILITY_GATE,
    RADIAL_FORCE_PIECE_CLOSURE_RELATIVE_GATE,
    SOURCE_RADIAL_FORCE_FORMULA,
    Z_DERIVATIVE_STEP_LADDER,
    _differentiate_axial_stress_ladder,
    centered_axial_stress_derivative,
    materialize_current_partial_nonlinear_radial_force,
    truth_boundary,
)
from openai_ns_reconstruction.kokuno_strict_inner_transport_radial_stress import (
    StrictInnerTransportRadialGeometry,
)


def _geometry() -> StrictInnerTransportRadialGeometry:
    radii = np.linspace(0.02, 0.44, 22)
    return StrictInnerTransportRadialGeometry(
        time=0.39,
        axial_z=-0.08,
        radii=tuple(float(v) for v in radii),
        bump_center=0.30,
        bump_halfwidth=0.10,
    )


def test_centered_axial_stress_derivative_is_exact_for_affine_z_profile():
    radii = np.asarray(_geometry().radii)
    step = 0.007
    z0 = -0.13
    slope = 0.4 + 0.7 * radii - 0.2 * radii**2
    intercept = np.sin(radii)
    plus = intercept + (z0 + step) * slope
    minus = intercept + (z0 - step) * slope
    derivative = centered_axial_stress_derivative(plus, minus, step)
    assert np.allclose(derivative, slope, rtol=0.0, atol=2.0e-14)


def test_source_radial_force_uses_same_axial_stress_for_all_nonlinear_pieces():
    geometry = _geometry()
    radii = np.asarray(geometry.radii)

    quadratic_slope = 0.3 * radii + 0.2 * radii**2
    mixed_slope = -0.1 + 0.4 * radii
    aggregate_slope = quadratic_slope + mixed_slope

    def sampler(z_value: float):
        base = np.cos(0.7 * radii)
        quadratic = 0.2 * base + z_value * quadratic_slope
        mixed = -0.1 * base + z_value * mixed_slope
        return {
            "quadratic": quadratic,
            "mixed": mixed,
            "aggregate": quadratic + mixed,
        }

    result = _differentiate_axial_stress_ladder(geometry, sampler)
    assert tuple(result["steps"]) == Z_DERIVATIVE_STEP_LADDER
    assert result["aggregate_derivative_stability_preflight_passed"] is True
    assert result["aggregate_finest_pair_relative_rms_difference"] <= 1.0e-12

    for level in result["levels"]:
        assert level["piece_closure_relative_max"] <= RADIAL_FORCE_PIECE_CLOSURE_RELATIVE_GATE
        assert np.allclose(level["quadratic"], quadratic_slope, atol=2.0e-13, rtol=0.0)
        assert np.allclose(level["mixed"], mixed_slope, atol=2.0e-13, rtol=0.0)
        assert np.allclose(level["aggregate"], aggregate_slope, atol=3.0e-13, rtol=0.0)


def test_piece_attribution_mutation_fails_closed():
    geometry = _geometry()
    radii = np.asarray(geometry.radii)

    def sampler(z_value: float):
        quadratic = z_value * radii
        mixed = 2.0 * z_value * radii
        aggregate = quadratic + mixed
        if z_value > geometry.axial_z:
            aggregate = aggregate + 1.0e-5 * radii
        return {
            "quadratic": quadratic,
            "mixed": mixed,
            "aggregate": aggregate,
        }

    with pytest.raises(ValueError, match="piece attribution closure"):
        _differentiate_axial_stress_ladder(geometry, sampler)


def test_malformed_stress_sampler_fails_closed():
    geometry = _geometry()

    def missing_piece(_z_value: float):
        values = np.ones(len(geometry.radii))
        return {"quadratic": values, "aggregate": values}

    with pytest.raises(ValueError, match="omitted mixed"):
        _differentiate_axial_stress_ladder(geometry, missing_piece)

    def wrong_shape(_z_value: float):
        values = np.ones(len(geometry.radii) - 1)
        return {"quadratic": values, "mixed": values, "aggregate": values}

    with pytest.raises(ValueError, match="wrong radial shape"):
        _differentiate_axial_stress_ladder(geometry, wrong_shape)


def test_public_contract_has_no_stress_force_or_derivative_tuning_inputs():
    params = inspect.signature(materialize_current_partial_nonlinear_radial_force).parameters
    assert list(params) == ["backend", "geometry"]
    for forbidden in (
        "residual",
        "defect",
        "mean",
        "source",
        "stress",
        "radial_force",
        "pressure",
        "forcing",
        "gain",
        "damping",
        "spatial_step",
        "derivative_step",
        "z_step",
        "angular_order",
        "viscosity",
        "correction",
        "stage_budget",
        "scientific_threshold",
    ):
        assert forbidden not in params


def test_truth_boundary_answers_radial_channel_completeness_without_promotion():
    boundary = truth_boundary()
    assert SOURCE_RADIAL_FORCE_FORMULA == "(div T)_r = partial_z sigma_1"
    assert boundary["current_partial_axial_e1_stress_consumed"] is True
    assert boundary["current_partial_radial_force_from_dz_sigma1_materialized"] is True
    assert boundary["quadratic_mixed_aggregate_radial_force_attribution_materialized"] is True
    assert boundary["independent_third_radial_moment_inverse_introduced"] is False
    assert boundary["recorded_radial_mean_equated_to_radial_force"] is False
    assert boundary["agent2_curl_or_jacobian_reimplemented_by_agent3"] is False
    assert boundary["forbidden_public_parameters_absent"] is True
    assert boundary["fine_pair_relative_stability_gate"] == FINE_PAIR_RELATIVE_STABILITY_GATE
    assert boundary["radial_force_piece_closure_relative_gate"] == (
        RADIAL_FORCE_PIECE_CLOSURE_RELATIVE_GATE
    )

    for key in (
        "velocity_beyond_Xh_materialized",
        "outer_global_leading_velocity_materialized",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect",
        "full_three_component_ns_correction_object_materialized",
        "scoped_current_partial_stress_authorized_as_correction_target",
        "mean_correction_velocity_materialized",
        "current_real_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert boundary[key] is False, key
    assert boundary["final_normalized_momentum_gate"] == 1.0e-3
    assert boundary["final_normalized_divergence_gate"] == 1.0e-5
