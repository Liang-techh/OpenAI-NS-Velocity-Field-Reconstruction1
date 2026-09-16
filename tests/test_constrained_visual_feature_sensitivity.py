import numpy as np
import pytest

from openai_ns_reconstruction.constrained_visual_feature_sensitivity import (
    diagnose_visual_feature_sensitivity,
)


def test_detects_duplicate_and_dead_visual_directions():
    # p0 and p1 produce the same normalized feature direction; p2 is locally dead.
    def features(theta):
        p0, p1, p2 = theta
        return np.array([
            p0 + 2.0*p1 + p2*p2,
            3.0*p0 + 6.0*p1,
            5.0*p2*p2,
        ])

    report = diagnose_visual_feature_sensitivity(
        features,
        np.array([0.2, -0.1, 0.0]),
        steps=np.array([1e-4, 1e-4, 1e-4]),
        lower_bounds=np.array([-1.0, -1.0, -1.0]),
        upper_bounds=np.array([1.0, 1.0, 1.0]),
        feature_scales=np.ones(3),
        parameter_labels=("phi_base", "phi_duplicate", "unused_local"),
        feature_labels=("core_radius", "axial_extent", "pitch"),
    )

    assert report.numerical_rank == 1
    assert report.dead_parameter_labels == ("unused_local",)
    assert report.stencil_by_parameter == ("central", "central", "central")
    assert report.condition_number == pytest.approx(1.0)


def test_uses_bound_safe_one_sided_stencils():
    def features(theta):
        return np.array([2.0*theta[0] - theta[1], theta[0] + 4.0*theta[1]])

    report = diagnose_visual_feature_sensitivity(
        features,
        np.array([1.0, -1.0]),
        steps=np.array([0.1, 0.1]),
        lower_bounds=np.array([-1.0, -1.0]),
        upper_bounds=np.array([1.0, 1.0]),
        feature_scales=np.array([1.0, 1.0]),
    )
    assert report.stencil_by_parameter == ("backward", "forward")
    jac = np.asarray(report.raw_jacobian)
    np.testing.assert_allclose(jac, np.array([[2.0, -1.0], [1.0, 4.0]]), atol=1e-12)
    assert report.numerical_rank == 2


def test_rank_is_scale_aware_and_repeatable():
    def features(theta):
        return np.array([
            1000.0*theta[0] + theta[1],
            0.001*theta[0] + 2.0*theta[1],
        ])

    kwargs = dict(
        steps=np.array([1e-3, 1e-3]),
        lower_bounds=np.array([-1.0, -1.0]),
        upper_bounds=np.array([1.0, 1.0]),
        feature_scales=np.array([1000.0, 1.0]),
        parameter_labels=("phi_taper", "swirl_pitch"),
        feature_labels=("axial_extent", "pitch"),
    )
    a = diagnose_visual_feature_sensitivity(features, np.array([0.2, 0.3]), **kwargs)
    b = diagnose_visual_feature_sensitivity(features, np.array([0.2, 0.3]), **kwargs)
    assert a.to_dict() == b.to_dict()
    assert a.numerical_rank == 2
    assert all(v > 0 for v in a.parameter_response_norms)
    assert all(v > 0 for v in a.feature_reachability_norms)


def test_fail_closed_inputs_and_nonfinite_features():
    def ok(theta):
        return np.array([theta[0]])

    with pytest.raises(ValueError):
        diagnose_visual_feature_sensitivity(
            ok, np.array([2.0]), steps=np.array([0.1]),
            lower_bounds=np.array([-1.0]), upper_bounds=np.array([1.0]),
            feature_scales=np.array([1.0]),
        )
    with pytest.raises(ValueError):
        diagnose_visual_feature_sensitivity(
            ok, np.array([0.0]), steps=np.array([0.1]),
            lower_bounds=np.array([-1.0]), upper_bounds=np.array([1.0]),
            feature_scales=np.array([0.0]),
        )

    def bad(theta):
        return np.array([np.nan])

    with pytest.raises(ValueError):
        diagnose_visual_feature_sensitivity(
            bad, np.array([0.0]), steps=np.array([0.1]),
            lower_bounds=np.array([-1.0]), upper_bounds=np.array([1.0]),
            feature_scales=np.array([1.0]),
        )
