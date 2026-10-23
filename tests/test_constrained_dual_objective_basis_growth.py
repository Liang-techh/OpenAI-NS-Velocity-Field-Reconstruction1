import numpy as np
import pytest

from openai_ns_reconstruction.constrained_dual_objective_basis_growth import (
    diagnose_nested_dual_objective_basis_growth,
)


def test_nested_growth_identifies_first_shared_residual_visual_capacity():
    residual_jacobian = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    visual_jacobian = np.array([
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ])

    report = diagnose_nested_dual_objective_basis_growth(
        residual_jacobian,
        np.array([1.0, 1.0]),
        visual_jacobian,
        np.array([1.0, 1.0]),
        residual_scales=np.ones(2),
        visual_scales=np.ones(2),
        lower_bounds=np.zeros(3),
        upper_bounds=np.ones(3),
        nested_dimensions=(1, 2, 3),
        parameter_labels=("residual_only", "visual_only", "shared"),
    )

    assert report.levels[0].residual_remaining_ratio == pytest.approx(1 / np.sqrt(2))
    assert report.levels[0].visual_remaining_ratio == pytest.approx(1.0)
    assert report.levels[0].jointly_improves_both is False

    assert report.levels[1].residual_remaining_ratio == pytest.approx(1 / np.sqrt(2))
    assert report.levels[1].visual_remaining_ratio == pytest.approx(1 / np.sqrt(2))
    assert report.levels[1].jointly_improves_both is True
    assert report.smallest_jointly_improving_dimension == 2

    assert report.levels[2].residual_remaining_ratio < 1e-10
    assert report.levels[2].visual_remaining_ratio < 1e-10
    assert report.levels[2].joint_remaining_ratio < 1e-10


def test_shared_coefficient_exposes_residual_visual_sign_conflict():
    report = diagnose_nested_dual_objective_basis_growth(
        np.array([[1.0]]),
        np.array([1.0]),
        np.array([[1.0]]),
        np.array([-1.0]),
        residual_scales=np.ones(1),
        visual_scales=np.ones(1),
        lower_bounds=np.array([-1.0]),
        upper_bounds=np.array([1.0]),
        nested_dimensions=(1,),
        parameter_labels=("conflicted",),
    )
    level = report.levels[0]
    assert level.coefficients[0] == pytest.approx(0.0, abs=1e-12)
    assert level.residual_remaining_ratio == pytest.approx(1.0)
    assert level.visual_remaining_ratio == pytest.approx(1.0)
    assert level.jointly_improves_both is False
    assert report.smallest_jointly_improving_dimension is None


def test_declared_bounds_limit_capacity_and_are_reported():
    report = diagnose_nested_dual_objective_basis_growth(
        np.array([[1.0]]),
        np.array([2.0]),
        np.array([[1.0]]),
        np.array([2.0]),
        residual_scales=np.ones(1),
        visual_scales=np.ones(1),
        lower_bounds=np.array([-0.5]),
        upper_bounds=np.array([0.5]),
        nested_dimensions=(1,),
        parameter_labels=("bounded_mode",),
    )
    level = report.levels[0]
    assert level.coefficients[0] == pytest.approx(0.5, abs=1e-10)
    assert level.residual_remaining_ratio == pytest.approx(0.75, abs=1e-10)
    assert level.visual_remaining_ratio == pytest.approx(0.75, abs=1e-10)
    assert level.active_upper_bounds == ("bounded_mode",)


def test_fail_closed_on_mismatched_or_posthoc_inputs():
    common = dict(
        residual_correction_target=np.ones(2),
        visual_jacobian=np.ones((2, 2)),
        visual_correction_target=np.ones(2),
        residual_scales=np.ones(2),
        visual_scales=np.ones(2),
        lower_bounds=-np.ones(2),
        upper_bounds=np.ones(2),
        nested_dimensions=(1, 2),
    )

    with pytest.raises(ValueError, match="share parameter columns"):
        diagnose_nested_dual_objective_basis_growth(
            np.ones((2, 3)), **common
        )

    bad = dict(common)
    bad["nested_dimensions"] = (2, 1)
    with pytest.raises(ValueError, match="strictly increasing"):
        diagnose_nested_dual_objective_basis_growth(np.ones((2, 2)), **bad)

    bad = dict(common)
    bad["visual_scales"] = np.array([1.0, 0.0])
    with pytest.raises(ValueError, match="strictly positive"):
        diagnose_nested_dual_objective_basis_growth(np.ones((2, 2)), **bad)

    bad = dict(common)
    bad["visual_correction_target"] = np.zeros(2)
    with pytest.raises(ValueError, match="must be nonzero"):
        diagnose_nested_dual_objective_basis_growth(np.ones((2, 2)), **bad)
