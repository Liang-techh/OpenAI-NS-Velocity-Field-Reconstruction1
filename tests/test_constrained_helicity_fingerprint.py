import numpy as np
import pytest

from openai_ns_reconstruction.constrained_helicity_fingerprint import (
    audit_helicity_resolution,
    diagnose_helicity_fingerprint,
)


def positive_beltrami(points, time):
    z = np.asarray(points)[:, 2]
    return np.column_stack((np.sin(z), np.cos(z), np.zeros_like(z)))


def negative_beltrami(points, time):
    z = np.asarray(points)[:, 2]
    return np.column_stack((np.sin(z), -np.cos(z), np.zeros_like(z)))


def shear(points, time):
    y = np.asarray(points)[:, 1]
    return np.column_stack((y, np.zeros_like(y), np.zeros_like(y)))


def test_positive_and_negative_beltrami_handedness():
    pos = diagnose_helicity_fingerprint(positive_beltrami, 0.5, box_half_width=1.0, grid_size=17, provenance="analytic +curl eigenfield")
    neg = diagnose_helicity_fingerprint(negative_beltrami, 0.5, box_half_width=1.0, grid_size=17, provenance="analytic -curl eigenfield")
    assert pos.signed_alignment == pytest.approx(1.0, abs=1e-12)
    assert pos.absolute_alignment == pytest.approx(1.0, abs=1e-12)
    assert pos.positive_alignment_weight_fraction == pytest.approx(1.0)
    assert pos.negative_alignment_weight_fraction == pytest.approx(0.0)
    assert neg.signed_alignment == pytest.approx(-1.0, abs=1e-12)
    assert neg.absolute_alignment == pytest.approx(1.0, abs=1e-12)
    assert neg.positive_alignment_weight_fraction == pytest.approx(0.0)
    assert neg.negative_alignment_weight_fraction == pytest.approx(1.0)
    assert pos.pde_validated is False and neg.visual_correspondence_verified is False


def test_orthogonal_shear_has_zero_helicity_alignment():
    fp = diagnose_helicity_fingerprint(shear, 0.5, box_half_width=1.0, grid_size=17, provenance="analytic shear")
    assert fp.helicity_integral == pytest.approx(0.0, abs=1e-14)
    assert fp.helicity_rms == pytest.approx(0.0, abs=1e-14)
    assert fp.signed_alignment == pytest.approx(0.0, abs=1e-14)
    assert fp.absolute_alignment == pytest.approx(0.0, abs=1e-14)
    assert fp.alignment_q10 == pytest.approx(0.0, abs=1e-14)
    assert fp.alignment_q90 == pytest.approx(0.0, abs=1e-14)


def test_resolution_audit_has_no_hidden_registration_or_threshold():
    audit = audit_helicity_resolution(positive_beltrami, 0.5, grid_sizes=(9, 13, 17), box_half_width=1.0, provenance="analytic resolution fixture")
    assert audit.grid_sizes == (9, 13, 17)
    assert audit.finest_grid_size == 17
    assert max(audit.signed_alignment_abs_delta_to_finest) < 1e-12
    assert max(audit.absolute_alignment_abs_delta_to_finest) < 1e-12
    assert audit.visualization_ready is False
    assert audit.pde_validated is False


def test_fail_closed_inputs():
    zero = lambda points, time: np.zeros_like(points)
    bad_shape = lambda points, time: np.zeros((len(points), 2))
    nonfinite = lambda points, time: np.full_like(points, np.nan)
    with pytest.raises(ValueError, match="identically zero"):
        diagnose_helicity_fingerprint(zero, 0.5, provenance="zero fixture")
    with pytest.raises(ValueError, match="shape"):
        diagnose_helicity_fingerprint(bad_shape, 0.5, provenance="bad shape")
    with pytest.raises(ValueError, match="finite"):
        diagnose_helicity_fingerprint(nonfinite, 0.5, provenance="nan fixture")
    with pytest.raises(ValueError, match="odd integer"):
        diagnose_helicity_fingerprint(positive_beltrami, 0.5, grid_size=16, provenance="bad grid")
    with pytest.raises(ValueError, match="provenance"):
        diagnose_helicity_fingerprint(positive_beltrami, 0.5, provenance="")
    with pytest.raises(ValueError, match="three grid"):
        audit_helicity_resolution(positive_beltrami, 0.5, grid_sizes=(9, 13), provenance="too short")
    with pytest.raises(ValueError, match="strictly increasing"):
        audit_helicity_resolution(positive_beltrami, 0.5, grid_sizes=(9, 17, 13), provenance="bad order")
