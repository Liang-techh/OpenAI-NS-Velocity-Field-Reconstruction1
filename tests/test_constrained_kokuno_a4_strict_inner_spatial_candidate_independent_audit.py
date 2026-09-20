from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_strict_inner_spatial_candidate_independent_audit import (
    DIVERGENCE_SAMPLED_MAX_GATE,
    DIVERGENCE_SAMPLED_RMS_GATE,
    FINAL_DIVERGENCE_GATE,
    FINAL_MOMENTUM_NORMALIZED_GATE,
    HELDOUT_SEED,
    RICHARDSON_SPATIAL_STEPS,
    component_axis_error_metrics,
    divergence_from_jacobian,
    divergence_gate_failures,
    frozen_gate_failures,
    public_contract,
    richardson_axis_derivative,
    richardson_jacobian_level,
    tensor_rms,
    vector_rms,
    vorticity_from_jacobian,
)


def _manufactured_velocity(x, y, z, t):
    x, y, z, t = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    # Divergence-free cubic cyclic field.  t is intentionally present only as
    # an evaluation coordinate; this test audits spatial derivatives.
    return np.stack((y**3 + z + 0.0 * t, z**3 + x, x**3 + y), axis=-1)


def _manufactured_jacobian(x, y, z, t):
    x, y, z, t = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    shape = x.shape + (3, 3)
    out = np.zeros(shape, dtype=float)
    out[..., 0, 1] = 3.0 * y * y
    out[..., 0, 2] = 1.0
    out[..., 1, 0] = 1.0
    out[..., 1, 2] = 3.0 * z * z
    out[..., 2, 0] = 3.0 * x * x
    out[..., 2, 1] = 1.0
    return out


def _points():
    return (
        np.asarray([0.17, -0.21, 0.08, 0.31]),
        np.asarray([-0.13, 0.24, 0.05, -0.19]),
        np.asarray([0.22, 0.11, -0.18, -0.07]),
        np.asarray([0.47, 0.50, 0.53, 0.49]),
    )


def test_richardson_jacobian_locks_component_axis_convention():
    x, y, z, t = _points()
    expected = _manufactured_jacobian(x, y, z, t)
    levels = tuple(
        richardson_jacobian_level(
            _manufactured_velocity, x, y, z, t, step=step
        )
        for step in RICHARDSON_SPATIAL_STEPS
    )

    for level in levels:
        np.testing.assert_allclose(level.jacobian, expected, atol=2.0e-11, rtol=2.0e-11)

    np.testing.assert_allclose(divergence_from_jacobian(levels[-1].jacobian), 0.0, atol=2.0e-11)
    expected_vorticity = np.stack(
        (
            np.ones_like(x) - 3.0 * z * z,
            np.ones_like(x) - 3.0 * x * x,
            np.ones_like(x) - 3.0 * y * y,
        ),
        axis=-1,
    )
    np.testing.assert_allclose(
        vorticity_from_jacobian(levels[-1].jacobian),
        expected_vorticity,
        atol=2.0e-11,
        rtol=2.0e-11,
    )
    assert frozen_gate_failures(expected, levels) == []


def test_mutation_firewall_detects_scale_sign_and_divergence_injection():
    x, y, z, t = _points()
    production = _manufactured_jacobian(x, y, z, t)
    levels = tuple(
        richardson_jacobian_level(
            _manufactured_velocity, x, y, z, t, step=step
        )
        for step in RICHARDSON_SPATIAL_STEPS
    )

    scaled = 0.99 * production
    assert "fine_relative_rms" in frozen_gate_failures(scaled, levels)

    sign_flip = production.copy()
    sign_flip[..., 0, 2] *= -1.0
    sign_failures = frozen_gate_failures(sign_flip, levels)
    assert "fine_relative_rms" in sign_failures
    assert "fine_relative_sampled_max" in sign_failures

    injected = levels[-1].jacobian.copy()
    injected[..., 0, 0] += 1.0e-3
    failures = divergence_gate_failures(injected)
    assert "divergence_sampled_max" in failures
    assert "divergence_sampled_rms" in failures
    divergence = divergence_from_jacobian(injected)
    assert float(np.max(np.abs(divergence))) > DIVERGENCE_SAMPLED_MAX_GATE
    assert float(np.sqrt(np.mean(divergence * divergence))) > DIVERGENCE_SAMPLED_RMS_GATE


def test_metrics_are_tensor_and_vector_scoped():
    x, y, z, t = _points()
    jac = _manufactured_jacobian(x, y, z, t)
    perturbed = jac.copy()
    perturbed[..., 2, 0] += 1.0e-6
    metrics = component_axis_error_metrics(jac, perturbed)
    assert len(metrics) == 9
    wx = next(item for item in metrics if item["component"] == "w" and item["axis"] == "x")
    assert wx["absolute_rms"] == pytest.approx(1.0e-6)
    assert tensor_rms(jac) > 0.0
    assert vector_rms(_manufactured_velocity(x, y, z, t)) > 0.0


def test_axis_derivative_rejects_invalid_axis_and_bad_provider_shape():
    x, y, z, t = _points()
    with pytest.raises(ValueError, match="axis"):
        richardson_axis_derivative(
            _manufactured_velocity, x, y, z, t, axis=3, step=1.0e-3
        )

    def bad_velocity(x, y, z, t):
        x = np.asarray(x, dtype=float)
        return np.zeros(x.shape + (2,), dtype=float)

    with pytest.raises(ValueError, match="shape"):
        richardson_axis_derivative(
            bad_velocity, x, y, z, t, axis=0, step=1.0e-3
        )


def test_public_contract_freezes_independence_and_final_gates():
    contract = public_contract()
    assert contract["heldout_seed"] == HELDOUT_SEED == 9173541
    assert contract["richardson_spatial_steps"] == [2.4e-3, 1.2e-3, 6.0e-4]
    independent = contract["independent_reference"]
    assert independent["artifact_must_be_saved_and_reloaded"] is True
    assert independent["numerical_reference_uses"] == "reloaded public velocity(x,y,z,t) only"
    assert independent["uses_agent1_analytic_jacobian_as_reference"] is False
    assert independent["uses_agent2_oscillatory_fd6_as_reference"] is False
    assert independent["uses_agent2_fd4_verifier_as_reference"] is False
    assert independent["uses_candidate_component_tensors"] is False

    gates = contract["frozen_gates"]
    assert gates["final_momentum_normalized_max_l2"] == FINAL_MOMENTUM_NORMALIZED_GATE == 1.0e-3
    assert gates["final_divergence_max_l2"] == FINAL_DIVERGENCE_GATE == 1.0e-5
    assert gates["independent_divergence_sampled_max"] == 1.0e-5
    assert gates["independent_divergence_sampled_rms"] == 1.0e-5

    scope = contract["norm_scope"]
    assert scope["jacobian_consistency_is_pde_residual"] is False
    assert scope["divergence_rms_is_heldout_sampled_rms"] is True
    assert scope["whole_domain_volume_l2_assessed"] is False
    assert scope["complete_ns_residual_assessed"] is False
    assert scope["same_protocol_st006_comparison_valid"] is False

    truth = contract["truth_boundary"]
    assert truth["strict_inner_spatial_candidate_artifact_independently_consumed"] is True
    assert truth["strict_inner_velocity_jacobian_independently_assessed"] is True
    for key in (
        "global_leading_velocity_materialized",
        "outer_join_materialized",
        "matched_pressure_included",
        "restricted_forcing_included",
        "correction_velocity_included",
        "complete_kokuno_candidate_assembled",
        "leading_only_ns_residual_assessed",
        "leading_plus_oscillatory_ns_residual_assessed",
        "after_correction_ns_residual_assessed",
        "heldout_ns_momentum_residual_assessed",
        "formal_full_domain_pde_gate_assessed",
        "whole_domain_boundary_support_gate_assessed",
        "pde_validated",
        "paper_exact",
    ):
        assert truth[key] is False
