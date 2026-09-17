import numpy as np
import pytest

from openai_ns_reconstruction.constrained_public_morphology_trend import (
    compare_ordered_morphology_trends,
)


def test_affine_feature_scaling_preserves_order_and_increment_profile():
    candidate = np.array([[1.0, 8.0], [2.0, 6.0], [4.0, 4.0], [7.0, 1.0]])
    reference = np.column_stack((10.0 + 3.0 * candidate[:, 0], 50.0 + 2.0 * candidate[:, 1]))
    report = compare_ordered_morphology_trends(
        candidate,
        reference,
        ["axial_extent", "radial_extent"],
        candidate_id="candidate-a",
        reference_source="public-frames-v1",
        pairing_provenance="same published frame order; no hidden-time inference",
        candidate_times=[0.25, 0.4, 0.6, 0.75],
    )
    assert report.pairing_mode == "ordinal_equal_length"
    assert report.mean_spearman_rho == pytest.approx(1.0)
    assert report.mean_normalized_increment_rms == pytest.approx(0.0)
    assert report.all_endpoint_directions_agree
    assert not report.time_registration_performed
    assert not report.visual_correspondence_verified
    assert not report.pde_validated


def test_reversed_public_trend_is_not_hidden_by_scaling():
    candidate = np.array([[1.0], [2.0], [3.0], [4.0]])
    reference = np.array([[8.0], [6.0], [4.0], [2.0]])
    report = compare_ordered_morphology_trends(
        candidate,
        reference,
        ["axial_extent"],
        candidate_id="candidate-a",
        reference_source="public-frames-v1",
        pairing_provenance="declared ordinal comparison",
    )
    feature = report.features[0]
    assert feature.spearman_rho == pytest.approx(-1.0)
    assert not feature.endpoint_direction_agrees
    assert feature.normalized_increment_rms > 0.0


def test_monotone_nonlinear_pacing_has_rho_one_but_nonzero_increment_error():
    candidate = np.array([[0.0], [1.0], [2.0], [3.0]])
    reference = np.array([[0.0], [0.2], [0.8], [3.0]])
    report = compare_ordered_morphology_trends(
        candidate,
        reference,
        ["aspect_ratio"],
        candidate_id="candidate-a",
        reference_source="public-frames-v1",
        pairing_provenance="declared ordinal comparison",
    )
    assert report.features[0].spearman_rho == pytest.approx(1.0)
    assert report.features[0].normalized_increment_rms > 0.2


def test_explicit_monotone_pairs_allow_different_frame_counts_without_registration():
    candidate = np.arange(6.0)[:, None]
    reference = (20.0 + np.arange(8.0))[:, None]
    report = compare_ordered_morphology_trends(
        candidate,
        reference,
        ["axial_extent"],
        candidate_id="candidate-a",
        reference_source="public-video-keyframes",
        pairing_provenance="manual keyframe pairs declared before scoring",
        candidate_times=[0.25, 0.35, 0.45, 0.55, 0.65, 0.75],
        frame_pairs=[[0, 0], [2, 3], [5, 7]],
    )
    assert report.pairing_mode == "explicit_monotone_pairs"
    assert report.candidate_indices == (0, 2, 5)
    assert report.reference_indices == (0, 3, 7)
    assert report.features[0].spearman_rho == pytest.approx(1.0)
    assert not report.hidden_frame_times_inferred


def test_different_frame_counts_require_explicit_pairs():
    with pytest.raises(ValueError, match="automatic time alignment is forbidden"):
        compare_ordered_morphology_trends(
            np.arange(4.0)[:, None],
            np.arange(5.0)[:, None],
            ["extent"],
            candidate_id="c",
            reference_source="r",
            pairing_provenance="p",
        )


@pytest.mark.parametrize(
    "pairs",
    [
        [[0, 0], [2, 2], [1, 3]],
        [[0, 0], [1, 2], [2, 1]],
        [[0, 0], [0, 1], [2, 2]],
    ],
)
def test_frame_pairs_must_preserve_order(pairs):
    with pytest.raises(ValueError, match="strictly increasing"):
        compare_ordered_morphology_trends(
            np.arange(4.0)[:, None],
            np.arange(4.0)[:, None],
            ["extent"],
            candidate_id="c",
            reference_source="r",
            pairing_provenance="p",
            frame_pairs=pairs,
        )


def test_constant_feature_fails_closed():
    with pytest.raises(ValueError, match="constant"):
        compare_ordered_morphology_trends(
            np.ones((4, 1)),
            np.arange(4.0)[:, None],
            ["extent"],
            candidate_id="c",
            reference_source="r",
            pairing_provenance="p",
        )


def test_nonfinite_and_bad_times_fail_closed():
    good = np.arange(4.0)[:, None]
    bad = good.copy()
    bad[2, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        compare_ordered_morphology_trends(
            bad,
            good,
            ["extent"],
            candidate_id="c",
            reference_source="r",
            pairing_provenance="p",
        )
    with pytest.raises(ValueError, match="strictly increasing"):
        compare_ordered_morphology_trends(
            good,
            good,
            ["extent"],
            candidate_id="c",
            reference_source="r",
            pairing_provenance="p",
            candidate_times=[0.25, 0.5, 0.4, 0.75],
        )


def test_names_and_provenance_are_required():
    good = np.column_stack((np.arange(4.0), np.arange(4.0) ** 2))
    with pytest.raises(ValueError, match="unique"):
        compare_ordered_morphology_trends(
            good,
            good,
            ["extent", "extent"],
            candidate_id="c",
            reference_source="r",
            pairing_provenance="p",
        )
    with pytest.raises(ValueError, match="pairing_provenance"):
        compare_ordered_morphology_trends(
            good,
            good,
            ["x", "y"],
            candidate_id="c",
            reference_source="r",
            pairing_provenance="",
        )
