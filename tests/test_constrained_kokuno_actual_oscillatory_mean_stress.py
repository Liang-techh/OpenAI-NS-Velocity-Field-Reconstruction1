from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_actual_oscillatory_mean_stress import (
    ADMITTED_AGENT2_HEAD,
    INDEPENDENT_AGENT4_HEAD,
    _compact_radial_stress,
    materialize_actual_oscillatory_mean_stress,
)


def test_compact_radial_stress_enforces_weighted_moment_complement() -> None:
    radii = np.linspace(0.15, 1.35, 65)
    source = (radii - 0.15) ** 2 * (1.35 - radii) ** 2 * (1.0 + 0.3 * radii)
    for exponent in (1, 2):
        receipt = _compact_radial_stress(
            radii,
            source,
            exponent=exponent,
            bump_center=0.75,
            bump_halfwidth=0.60,
        )
        assert receipt["bump_weighted_integral"] == pytest.approx(1.0, abs=2.0e-14)
        assert abs(receipt["moment_complement_weighted_moment"]) <= 2.0e-14
        assert receipt["stress_inner_edge"] == 0.0
        assert abs(receipt["stress_outer_edge"]) <= 2.0e-13
        assert receipt["stress_rms"] > 0.0


def test_actual_materializer_has_no_caller_supplied_defect_target() -> None:
    names = set(inspect.signature(materialize_actual_oscillatory_mean_stress).parameters)
    assert names == {"radial_count", "angular_count", "time", "z", "fd4_step", "nu"}
    assert not ({"defect", "residual", "target", "stress"} & names)


def test_actual_admitted_field_materializes_nonzero_mean_stress_fail_closed() -> None:
    receipt = materialize_actual_oscillatory_mean_stress(
        radial_count=17,
        angular_count=16,
        time=0.5,
        z=0.08,
        fd4_step=0.005,
        nu=0.01,
    )
    assert receipt["provenance"]["admitted_agent2_head"] == ADMITTED_AGENT2_HEAD
    assert receipt["provenance"]["independent_agent4_head"] == INDEPENDENT_AGENT4_HEAD
    assert receipt["anti_surrogate_contract"]["surrogate_defect_used"] is False
    assert receipt["anti_surrogate_contract"]["caller_supplied_target_parameters"] == []

    mean = receipt["mean_operator"]
    stress = receipt["requested_stress_component"]
    assert np.isfinite(mean["raw_mean_vector_rms"])
    assert mean["raw_mean_vector_rms"] > 0.0
    assert np.isfinite(mean["quadratic_mean_vector_rms"])
    assert mean["quadratic_mean_vector_rms"] > 0.0
    assert mean["mean_decomposition_closure_max_abs"] <= 1.0e-9 * max(
        mean["raw_mean_sampled_max"], 1.0
    )
    assert stress["requested_stress_rms"] > 0.0
    assert abs(stress["theta_e2"]["moment_complement_weighted_moment"]) <= 1.0e-8
    assert abs(stress["axial_e1"]["moment_complement_weighted_moment"]) <= 1.0e-8
    assert abs(stress["theta_e2"]["stress_outer_edge"]) <= 1.0e-8
    assert abs(stress["axial_e1"]["stress_outer_edge"]) <= 1.0e-8

    truth = receipt["truth_boundary"]
    assert truth["admitted_agent2_velocity_consumed"] is True
    assert truth["exact_candidate_time_derivative_consumed"] is True
    assert truth["real_oscillatory_self_defect_component_consumed"] is True
    assert truth["real_oscillatory_quadratic_mean_component_measured"] is True
    assert truth["oscillatory_requested_stress_component_materialized"] is True
    assert truth["full_same_cycle_composite_requested_stress_materialized"] is False
    assert truth["candidate_finite_head_mean_debt_materialized"] is False
    assert truth["signed_mean_inverse_input_ready"] is False
    assert truth["finite_correction_cycle_rerun_allowed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["pde_validated"] is False


def test_invalid_state_and_resolution_fail_closed() -> None:
    with pytest.raises(ValueError):
        materialize_actual_oscillatory_mean_stress(radial_count=15)
    with pytest.raises(ValueError):
        materialize_actual_oscillatory_mean_stress(angular_count=8)
    with pytest.raises(ValueError):
        materialize_actual_oscillatory_mean_stress(time=0.1)
    with pytest.raises(ValueError):
        materialize_actual_oscillatory_mean_stress(z=2.0)
    with pytest.raises(ValueError):
        materialize_actual_oscillatory_mean_stress(fd4_step=0.0)
