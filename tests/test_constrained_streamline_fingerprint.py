import math

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_streamline_fingerprint import (
    audit_streamline_step_resolution,
    diagnose_streamline_fingerprint,
)


def helical_velocity(points, time):
    points = np.asarray(points, dtype=float)
    x = points[..., 0]
    y = points[..., 1]
    return np.stack((-2.0 * y, 2.0 * x, np.ones_like(x)), axis=-1)


def axial_velocity(points, time):
    points = np.asarray(points, dtype=float)
    out = np.zeros_like(points)
    out[..., 2] = 1.0
    return out


def ring_seeds(n=8, radius=1.0):
    phi = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return np.stack((radius * np.cos(phi), radius * np.sin(phi), np.zeros_like(phi)), axis=-1)


def test_helical_streamline_geometry_matches_analytic_pitch_and_extent():
    report = diagnose_streamline_fingerprint(
        helical_velocity, 0.5, ring_seeds(), half_arclength=5.0, step=0.025
    )
    expected_z_span = 10.0 / math.sqrt(5.0)
    expected_turns = 10.0 / (math.pi * math.sqrt(5.0))
    assert report["complete_fraction"] == 1.0
    assert report["speed_floor_fraction"] == 0.0
    assert report["pitch_per_turn"]["q50"] == pytest.approx(math.pi, rel=2e-7)
    assert report["z_span"]["q50"] == pytest.approx(expected_z_span, rel=2e-7)
    assert report["abs_turns"]["q50"] == pytest.approx(expected_turns, rel=2e-7)
    assert report["radial_drift"]["q90"] < 1e-8
    assert report["signed_turn_coherence"] == pytest.approx(1.0, abs=1e-12)
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False


def test_three_step_resolution_is_stable_for_smooth_helix():
    audit = audit_streamline_step_resolution(
        helical_velocity,
        0.5,
        ring_seeds(),
        half_arclength=5.0,
        steps=(0.1, 0.05, 0.025),
    )
    pitch_entries = audit["comparisons"]["pitch_per_turn.q50"]
    turns_entries = audit["comparisons"]["abs_turns.q50"]
    assert len(audit["reports"]) == 3
    assert pitch_entries[0]["relative_delta_to_finest"] < 2e-6
    assert turns_entries[0]["relative_delta_to_finest"] < 2e-6
    assert audit["comparisons"]["radial_drift.q90"][0]["absolute_delta_to_finest"] < 2e-6


def test_axial_field_has_zero_turns_and_null_pitch_without_false_helix():
    report = diagnose_streamline_fingerprint(
        axial_velocity,
        0.5,
        np.array([[1.0, 0.0, 0.0], [0.5, 0.0, 0.25]]),
        half_arclength=2.0,
        step=0.05,
    )
    assert report["abs_turns"]["q90"] == pytest.approx(0.0, abs=1e-14)
    assert math.isnan(report["pitch_per_turn"]["q50"])
    assert report["signed_turn_coherence"] is None
    assert report["z_span"]["q50"] == pytest.approx(4.0, rel=1e-12)


def test_fail_closed_on_bad_inputs_and_nonfinite_velocity():
    seeds = ring_seeds(4)
    with pytest.raises(ValueError, match="three"):
        audit_streamline_step_resolution(helical_velocity, 0.5, seeds, steps=(0.1, 0.05))
    with pytest.raises(ValueError, match="strictly decreasing"):
        audit_streamline_step_resolution(
            helical_velocity, 0.5, seeds, steps=(0.05, 0.1, 0.025)
        )
    with pytest.raises(ValueError, match="shape"):
        diagnose_streamline_fingerprint(helical_velocity, 0.5, np.zeros((3, 2)))

    def bad_velocity(points, time):
        out = np.zeros_like(points, dtype=float)
        out[..., 0] = np.nan
        return out

    with pytest.raises(ValueError, match="non-finite"):
        diagnose_streamline_fingerprint(
            bad_velocity, 0.5, seeds, half_arclength=1.0, step=0.1
        )
