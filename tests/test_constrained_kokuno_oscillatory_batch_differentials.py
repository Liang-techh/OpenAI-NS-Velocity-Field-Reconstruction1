from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_oscillatory_batch_axis_safety import (
    _support,
    _time_interval,
)
from openai_ns_reconstruction.kokuno_oscillatory_batch_differentials import (
    FIXED_SPATIAL_STEP,
    SCALAR_REPLAY_ATOL,
    _curl_from_jacobian,
    batch_differentials_sha256,
    evaluate_oscillatory_batch_differentials,
    materialize_batch_differentials_receipt,
    public_contract,
    semantic_payload,
)
from openai_ns_reconstruction.kokuno_public_oscillatory_vorticity_diagnostic import (
    evaluate_vorticity_osc_fd6,
)


POINTS = np.asarray(
    [
        (0.42, 0.19, -1.10),
        (0.63, -0.27, -0.45),
        (-0.51, 0.38, 0.15),
        (-0.74, -0.22, 0.65),
    ],
    dtype=float,
)
TIMES = np.asarray((0.31, 0.41, 0.53, 0.67), dtype=float)


def test_batch_differentials_replay_existing_fd6_scalar_surface() -> None:
    got = evaluate_oscillatory_batch_differentials(POINTS, TIMES)
    assert got.velocity.shape == (4, 3)
    assert got.velocity_jacobian.shape == (4, 3, 3)
    assert got.divergence.shape == (4,)
    assert got.vorticity.shape == (4, 3)
    assert got.interior_mask.shape == (4,)
    assert np.all(got.interior_mask)

    scalar = [
        evaluate_vorticity_osc_fd6(x, y, z, t, spatial_step=FIXED_SPATIAL_STEP)
        for (x, y, z), t in zip(POINTS, TIMES, strict=True)
    ]
    expected_j = np.stack(
        [np.asarray(item["velocity_gradient_fd6"], dtype=float).reshape(3, 3) for item in scalar]
    )
    expected_w = np.stack(
        [np.asarray(item["vorticity_fd6"], dtype=float).reshape(3) for item in scalar]
    )
    expected_d = np.asarray(
        [np.asarray(item["divergence_fd6"], dtype=float).reshape(()).item() for item in scalar]
    )
    np.testing.assert_allclose(got.velocity_jacobian, expected_j, rtol=0.0, atol=SCALAR_REPLAY_ATOL)
    np.testing.assert_allclose(got.vorticity, expected_w, rtol=0.0, atol=SCALAR_REPLAY_ATOL)
    np.testing.assert_allclose(got.divergence, expected_d, rtol=0.0, atol=SCALAR_REPLAY_ATOL)
    np.testing.assert_array_equal(got.vorticity, _curl_from_jacobian(got.velocity_jacobian))
    np.testing.assert_array_equal(
        got.divergence, np.trace(got.velocity_jacobian, axis1=-2, axis2=-1)
    )


def test_higher_rank_shape_time_broadcast_and_order_invariance() -> None:
    cloud = POINTS.reshape(2, 2, 3)
    by_scalar_time = evaluate_oscillatory_batch_differentials(cloud, 0.47)
    by_full_time = evaluate_oscillatory_batch_differentials(
        cloud, np.full((2, 2), 0.47, dtype=float)
    )
    np.testing.assert_array_equal(by_scalar_time.velocity, by_full_time.velocity)
    np.testing.assert_array_equal(by_scalar_time.velocity_jacobian, by_full_time.velocity_jacobian)
    np.testing.assert_array_equal(by_scalar_time.vorticity, by_full_time.vorticity)
    assert by_scalar_time.velocity.shape == (2, 2, 3)
    assert by_scalar_time.velocity_jacobian.shape == (2, 2, 3, 3)
    assert by_scalar_time.divergence.shape == (2, 2)

    perm = np.asarray((2, 0, 3, 1))
    shuffled = evaluate_oscillatory_batch_differentials(POINTS[perm], TIMES[perm])
    inverse = np.argsort(perm)
    np.testing.assert_array_equal(shuffled.velocity[inverse], evaluate_oscillatory_batch_differentials(POINTS, TIMES).velocity)
    np.testing.assert_array_equal(shuffled.velocity_jacobian[inverse], evaluate_oscillatory_batch_differentials(POINTS, TIMES).velocity_jacobian)
    np.testing.assert_array_equal(shuffled.vorticity[inverse], evaluate_oscillatory_batch_differentials(POINTS, TIMES).vorticity)


def test_axis_support_faces_and_exterior_are_exact_zero_for_all_fields() -> None:
    r0, r1, z0, z1 = _support()
    tiny = np.nextafter(0.0, 1.0)
    points = np.asarray(
        [
            (0.0, 0.0, 0.0),
            (tiny, -tiny, 0.25),
            (0.5 * r0, 0.0, -0.4),
            (r0, 0.0, 0.0),
            (r1, 0.0, 0.0),
            (0.5 * (r0 + r1), 0.0, z0),
            (0.5 * (r0 + r1), 0.0, z1),
            (1.01 * r1, 0.0, 0.3),
        ],
        dtype=float,
    )
    got = evaluate_oscillatory_batch_differentials(points, 0.47)
    assert not np.any(got.interior_mask)
    np.testing.assert_array_equal(got.velocity, np.zeros_like(got.velocity))
    np.testing.assert_array_equal(got.velocity_jacobian, np.zeros_like(got.velocity_jacobian))
    np.testing.assert_array_equal(got.divergence, np.zeros_like(got.divergence))
    np.testing.assert_array_equal(got.vorticity, np.zeros_like(got.vorticity))
    assert np.all(np.isfinite(got.velocity_jacobian))


def test_input_and_parent_time_domain_fail_closed() -> None:
    with pytest.raises(ValueError, match="shape"):
        evaluate_oscillatory_batch_differentials(np.zeros((3, 2)), 0.4)
    bad = POINTS.copy()
    bad[0, 1] = np.nan
    with pytest.raises(ValueError, match="finite"):
        evaluate_oscillatory_batch_differentials(bad, TIMES)
    with pytest.raises(ValueError, match="broadcast"):
        evaluate_oscillatory_batch_differentials(POINTS, np.zeros(3))
    with pytest.raises(ValueError, match="finite"):
        evaluate_oscillatory_batch_differentials(POINTS, np.nan)

    t0, t1 = _time_interval()
    # Time validity is checked before support masking; even the axis cannot bypass it.
    with pytest.raises(ValueError, match="registered candidate interval"):
        evaluate_oscillatory_batch_differentials(np.asarray([(0.0, 0.0, 0.0)]), t0 - 1.0e-6)
    with pytest.raises(ValueError, match="registered candidate interval"):
        evaluate_oscillatory_batch_differentials(np.asarray([(0.0, 0.0, 0.0)]), t1 + 1.0e-6)


def test_public_api_has_no_derivative_or_residual_tuning_escape_hatch() -> None:
    params = set(inspect.signature(evaluate_oscillatory_batch_differentials).parameters)
    forbidden = {
        "amplitude", "phase", "phase_offset", "scale", "orientation", "support",
        "spatial_step", "derivative_step", "residual", "target", "forcing", "pressure",
        "viscosity", "nu", "gain", "threshold", "mean", "stress", "inverse",
    }
    assert not (params & forbidden)
    contract = public_contract()
    assert contract["forbidden_inputs_present"] == []
    assert contract["fixed_spatial_step"] == 1.0e-3
    assert contract["caller_tunable_step"] is False
    assert contract["axis_safe_by_parent_strict_support_mask"] is True
    assert contract["complete_curl_reimplemented"] is False
    assert contract["new_derivative_realization_introduced"] is False
    assert contract["mean_projection_performed"] is False
    assert contract["radial_inverse_performed"] is False
    assert contract["correction_velocity_constructed"] is False
    assert contract["complete_ns_residual"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False


def test_semantic_identity_and_receipt_truth_boundary() -> None:
    payload = semantic_payload()
    assert payload["parent_agent2_pr"] == 954
    assert payload["parent_agent2_head"] == "e43c32c3d0258f635e8b1f999d89c38524355856"
    assert payload["differential_realization"]["fixed_spatial_step"] == 1.0e-3
    assert payload["differential_realization"]["caller_tunable_step"] is False
    assert payload["truth_boundary"]["global_leading_velocity_materialized"] is False
    assert payload["truth_boundary"]["paper_exact"] is False
    assert payload["truth_boundary"]["pde_validated"] is False
    assert len(batch_differentials_sha256()) == 64

    receipt = materialize_batch_differentials_receipt()
    diagnostic = receipt["diagnostic"]
    assert diagnostic["jacobian_scalar_replay_max_abs"] <= SCALAR_REPLAY_ATOL
    assert diagnostic["vorticity_scalar_replay_max_abs"] <= SCALAR_REPLAY_ATOL
    assert diagnostic["divergence_scalar_replay_max_abs"] <= SCALAR_REPLAY_ATOL
    assert diagnostic["vorticity_curl_closure_max_abs"] == 0.0
    assert diagnostic["divergence_trace_closure_max_abs"] == 0.0
    assert diagnostic["masked_all_fields_abs_max"] == 0.0
    assert diagnostic["interior_vorticity_vector_rms"] > 1.0e-12
