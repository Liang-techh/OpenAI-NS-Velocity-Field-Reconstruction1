import numpy as np

from openai_ns_reconstruction.kokuno_independent_leading_core import (
    REGISTERED_PDE_THRESHOLD,
    evaluate_fd4,
    run_leading_core_preflight,
    sample_held_out_core_points,
)


def test_fd4_manufactured_polynomial_momentum_is_exact() -> None:
    def velocity(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(x, y, z, t)
        return np.stack((t * x + y, t * y + z, t * z + x), axis=-1)

    def pressure(x, y, z, t):
        del t
        return 0.5 * (x * x + y * y + z * z)

    points = np.array(
        [
            [0.13, -0.07, 0.19],
            [-0.11, 0.09, -0.17],
            [0.021, 0.031, -0.041],
        ],
        dtype=float,
    )
    time = 0.47
    values = evaluate_fd4(velocity, pressure, points, time, 0.013, nu=0.01)
    x, y, z = points.T
    expected = np.column_stack(
        (
            2.0 * x + time * time * x + 2.0 * time * y + z,
            2.0 * y + time * time * y + 2.0 * time * z + x,
            2.0 * z + time * time * z + 2.0 * time * x + y,
        )
    )
    np.testing.assert_allclose(values["divergence"], 3.0 * time, rtol=0.0, atol=2e-12)
    np.testing.assert_allclose(values["residual"], expected, rtol=0.0, atol=2e-11)


def test_held_out_sampling_is_reproducible_and_axis_near() -> None:
    first = sample_held_out_core_points(interior_count=8, axis_near_count=4, seed=1234)
    second = sample_held_out_core_points(interior_count=8, axis_near_count=4, seed=1234)
    np.testing.assert_array_equal(first["interior"], second["interior"])
    np.testing.assert_array_equal(first["axis_near"], second["axis_near"])
    radii = np.hypot(first["axis_near"][:, 0], first["axis_near"][:, 1])
    assert np.all((radii >= 2e-4) & (radii <= 2e-2))
    assert not np.any(first["interior"] == 0.0)


def test_public_leading_core_preflight_is_truth_bounded_and_mutation_sensitive() -> None:
    report = run_leading_core_preflight(
        interior_count=8,
        axis_near_count=4,
        seed=9172841,
        times=(0.5,),
        steps=(0.02, 0.01, 0.005),
    )
    assert report["schema"] == "kokuno-independent-leading-core-preflight-v1"
    assert report["operator"]["candidate_derivative_helpers_used"] is False
    assert report["operator"]["training_tensors_used"] is False
    assert report["operator"]["public_contract_only"] == ["velocity", "pressure"]
    assert report["formal_full_pde_gate"]["normalized_residual_threshold"] == REGISTERED_PDE_THRESHOLD
    assert report["formal_full_pde_gate"]["assessed"] is False
    assert report["formal_full_pde_gate"]["passed"] is False
    assert report["truth_boundary"]["global_support_validated"] is False
    assert report["truth_boundary"]["complete_kokuno_composite_velocity"] is False
    assert report["truth_boundary"]["pde_validated"] is False

    finest = report["finest_summary"]
    assert np.isfinite(finest["worst_local_normalized_residual_max"])
    assert np.isfinite(finest["worst_local_normalized_residual_rms"])
    assert finest["worst_local_normalized_residual_max"] > 0.0
    assert finest["worst_divergence_max"] >= 0.0

    for row in report["rows"]:
        assert row["interior"]["velocity_rms"] > 0.0
        assert row["axis_near"]["velocity_rms"] > 0.0
        assert np.isfinite(row["interior"]["volume_l2_residual_estimate"])

    mutation = report["mutation_calibration"]
    assert abs(mutation["divergence_shift_mean"] - 0.02) < 2e-7
    assert mutation["divergence_shift_max_abs_error_from_expected"] < 2e-7
    assert abs(mutation["pressure_gradient_x_shift_mean"] - 0.02) < 2e-7
    assert mutation["pressure_gradient_x_shift_max_abs_error_from_expected"] < 2e-7
    assert mutation["pressure_gradient_transverse_shift_max"] < 2e-7
