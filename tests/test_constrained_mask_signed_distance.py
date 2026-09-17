import numpy as np
import pytest

from openai_ns_reconstruction.constrained_mask_signed_distance import (
    build_signed_distance_visual_residual,
)


def _mask():
    mask = np.zeros((31, 41), dtype=bool)
    mask[8:24, 11:30] = True
    mask[13:18, 17:23] = False
    return mask


def _run(candidate, reference, spacing=(1.0, 1.0)):
    return build_signed_distance_visual_residual(
        candidate,
        reference,
        pixel_spacing=spacing,
        candidate_provenance="candidate render fixture",
        reference_provenance="public mask fixture",
    )


def test_identical_masks_are_zero_and_arrays_are_read_only():
    mask = _mask()
    result = _run(mask, mask)
    assert result.balanced_rms == pytest.approx(0.0)
    assert result.max_abs == pytest.approx(0.0)
    assert np.linalg.norm(result.weighted_residual_vector) == pytest.approx(0.0)
    assert not result.normalized_delta.flags.writeable
    assert result.truth_boundary["visual_correspondence_verified"] is False
    assert result.truth_boundary["pde_validated"] is False


def test_weighted_vector_norm_equals_balanced_rms():
    reference = _mask()
    candidate = np.roll(reference, 2, axis=1)
    result = _run(candidate, reference)
    assert np.linalg.norm(result.weighted_residual_vector) == pytest.approx(result.balanced_rms)
    assert result.foreground_rms > 0.0
    assert result.background_rms > 0.0


def test_no_hidden_translation_alignment():
    reference = _mask()
    candidate = np.roll(reference, 3, axis=0)
    result = _run(candidate, reference)
    assert result.balanced_rms > 0.01
    assert result.truth_boundary["alignment_performed"] is False


def test_missing_hole_is_visible_to_signed_distance():
    reference = _mask()
    candidate = reference.copy()
    candidate[13:18, 17:23] = True
    result = _run(candidate, reference)
    assert result.balanced_rms > 0.0
    assert result.max_abs > 0.0


def test_anisotropic_spacing_changes_physical_residual():
    reference = _mask()
    candidate = np.roll(reference, 2, axis=0)
    isotropic = _run(candidate, reference, (1.0, 1.0))
    stretched = _run(candidate, reference, (3.0, 1.0))
    assert stretched.balanced_rms != pytest.approx(isotropic.balanced_rms)


@pytest.mark.parametrize(
    "bad",
    [
        np.zeros((9, 9), dtype=bool),
        np.ones((9, 9), dtype=bool),
        np.full((9, 9), 0.5),
    ],
)
def test_invalid_or_collapsed_candidate_fails_closed(bad):
    reference = np.zeros((9, 9), dtype=bool)
    reference[2:7, 2:7] = True
    with pytest.raises(ValueError):
        _run(bad, reference)


def test_shape_nan_spacing_and_provenance_fail_closed():
    reference = _mask()
    with pytest.raises(ValueError):
        _run(reference[:-1], reference)

    bad = reference.astype(float)
    bad[0, 0] = np.nan
    with pytest.raises(ValueError):
        _run(bad, reference)

    with pytest.raises(ValueError):
        _run(reference, reference, (0.0, 1.0))

    with pytest.raises(ValueError):
        build_signed_distance_visual_residual(
            reference,
            reference,
            candidate_provenance="",
            reference_provenance="public fixture",
        )
