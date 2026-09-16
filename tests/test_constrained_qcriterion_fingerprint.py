import numpy as np
import pytest

from openai_ns_reconstruction.constrained_qcriterion_fingerprint import (
    audit_qcriterion_resolution,
    diagnose_qcriterion_fingerprint,
)


def test_solid_body_rotation_matches_exact_q_and_interface_contract():
    omega = 2.0

    def velocity(points, time):
        del time
        out = np.zeros_like(points, dtype=float)
        out[..., 0] = -omega * points[..., 1]
        out[..., 1] = omega * points[..., 0]
        return out

    result = diagnose_qcriterion_fingerprint(
        velocity,
        time=0.5,
        half_width=1.0,
        grid_size=17,
    )

    assert result.q_peak == pytest.approx(omega**2, abs=2e-13)
    assert result.q_rms == pytest.approx(omega**2, abs=2e-13)
    assert result.positive_q_volume_fraction == pytest.approx(1.0)
    assert result.sampled_divergence_rms < 2e-14
    assert result.positive_core_present is True
    assert result.claim_scope == "visualization_morphology_only"
    assert result.pde_validated is False
    assert result.openai_field_identified is False


def test_pure_strain_has_no_false_rotation_dominated_core():
    rate = 1.5

    def velocity(points, time):
        del time
        out = np.zeros_like(points, dtype=float)
        out[..., 0] = rate * points[..., 0]
        out[..., 1] = -rate * points[..., 1]
        return out

    result = diagnose_qcriterion_fingerprint(
        velocity,
        time=0.5,
        half_width=1.0,
        grid_size=17,
    )

    assert result.q_peak == pytest.approx(-(rate**2), abs=2e-13)
    assert result.positive_q_volume_fraction == 0.0
    assert result.positive_q_weight == 0.0
    assert result.positive_core_present is False
    assert result.positive_q_r50 is None
    assert result.positive_q_abs_z90 is None


def test_localized_swirl_resolution_audit_and_fail_closed_inputs():
    def velocity(points, time):
        del time
        x = points[..., 0]
        y = points[..., 1]
        z = points[..., 2]
        a = np.exp(-1.3 * (x*x + y*y) - 0.8 * z*z)
        out = np.empty_like(points, dtype=float)
        out[..., 0] = -y * a
        out[..., 1] = x * a
        out[..., 2] = 0.2 * z * a
        return out

    audit = audit_qcriterion_resolution(
        velocity,
        time=0.5,
        half_width=1.5,
        grid_sizes=(17, 25, 33),
    )
    assert len(audit.fingerprints) == 3
    assert audit.fingerprints[-1].positive_core_present is True
    assert audit.fingerprints[-1].positive_q_r90 is not None
    assert audit.relative_changes_to_finest[-1]["q_peak"] == 0.0
    assert audit.pde_validated is False

    with pytest.raises(ValueError, match="odd"):
        diagnose_qcriterion_fingerprint(velocity, time=0.5, half_width=1.0, grid_size=16)
    with pytest.raises(ValueError, match="three"):
        audit_qcriterion_resolution(
            velocity,
            time=0.5,
            half_width=1.0,
            grid_sizes=(17, 25),
        )

    def bad_shape(points, time):
        del time
        return np.zeros(points.shape[:-1], dtype=float)

    with pytest.raises(ValueError, match="shape"):
        diagnose_qcriterion_fingerprint(bad_shape, time=0.5, half_width=1.0, grid_size=17)

    def nonfinite(points, time):
        del time
        out = np.zeros_like(points, dtype=float)
        out[..., 0] = np.nan
        return out

    with pytest.raises(ValueError, match="non-finite"):
        diagnose_qcriterion_fingerprint(nonfinite, time=0.5, half_width=1.0, grid_size=17)
