import numpy as np
import pytest

from openai_ns_reconstruction.constrained_basis_rank import (
    screen_linearized_residual_capacity,
)


def test_rrqr_detects_redundant_and_zero_directions():
    jacobian = np.array(
        [
            [3.0, 0.0, 6.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]
    )
    result = screen_linearized_residual_capacity(
        [1.0, -2.0, 0.5],
        jacobian,
        labels=["radial", "axial", "radial_duplicate", "zero"],
    )

    assert result.numerical_rank == 2
    assert 1 in result.independent_indices
    assert len(set(result.independent_indices) & {0, 2}) == 1
    assert 3 in result.dependent_indices
    assert len(set(result.dependent_indices) & {0, 2}) == 1
    assert result.column_norms[3] == 0.0


def test_projection_reports_exact_first_order_capacity_when_residual_is_in_span():
    jacobian = np.array(
        [
            [1.0, 0.0],
            [0.0, 2.0],
            [1.0, 1.0],
            [-1.0, 0.5],
        ]
    )
    coefficients = np.array([0.75, -1.25])
    residual = jacobian @ coefficients

    result = screen_linearized_residual_capacity(residual, jacobian)

    assert result.numerical_rank == 2
    assert result.projected_l2 < 1e-12
    assert result.projected_reduction_fraction == pytest.approx(1.0, abs=1e-12)


def test_projection_preserves_orthogonal_floor_and_validation_fails_closed():
    jacobian = np.array(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [0.0, 0.0],
        ]
    )
    residual = np.array([2.0, -1.0, 3.0])

    result = screen_linearized_residual_capacity(
        residual,
        jacobian,
        labels=["mode_a", "mode_b"],
    )

    assert result.projected_l2 == pytest.approx(3.0, abs=1e-12)
    assert result.projected_l2 < result.baseline_l2
    assert result.independent_labels == ("mode_a", "mode_b")

    with pytest.raises(ValueError, match="row count"):
        screen_linearized_residual_capacity([1.0, 2.0], np.eye(3))
    with pytest.raises(ValueError, match="non-finite"):
        screen_linearized_residual_capacity([1.0, np.nan], np.eye(2))
    with pytest.raises(ValueError, match="unique"):
        screen_linearized_residual_capacity([1.0, 2.0], np.eye(2), labels=["x", "x"])
