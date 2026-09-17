import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_independent_radial_stress import (
    IndependentPchipRadialInverse,
    _window,
    run_real_candidate_preflight,
    validate_profile_independently,
)
from openai_ns_reconstruction.kokuno_radial_stress import (
    CompactRadialStressInverse,
    RadialMeanDefectProfile,
)


def _synthetic_profile(count=65):
    radii = np.linspace(0.12, 0.72, count)
    raw = np.sin(6.0 * radii) + 0.25 * np.cos(11.0 * radii) + 0.1 * radii
    values = raw * _window(radii, 0.12, 0.72, 5)
    return RadialMeanDefectProfile(
        radii=radii,
        raw_values=raw,
        values=values,
        exponent=2,
        component="theta",
        time=0.5,
        z=0.0,
        angular_count=16,
        annulus=(0.12, 0.72),
    )


def test_independent_pchip_inverse_agrees_with_public_stress_off_grid():
    profile = _synthetic_profile(65)
    reference = CompactRadialStressInverse(profile)
    independent = IndependentPchipRadialInverse(profile)
    query = np.linspace(0.137, 0.703, 211)

    public = reference.stress(query)
    alternate = independent.stress(query)
    normalized = np.sqrt(np.mean((public - alternate) ** 2)) / np.sqrt(np.mean(public**2))

    assert normalized < 5.0e-3
    assert abs(reference.moment - independent.moment) / abs(reference.moment) < 5.0e-3
    assert np.max(np.abs(reference.stress(np.array([0.02, 0.12, 0.72, 0.82])))) == 0.0
    assert np.max(np.abs(independent.stress(np.array([0.02, 0.12, 0.72, 0.82])))) == 0.0


def test_independent_validator_detects_mutated_public_stress():
    metrics = validate_profile_independently(_synthetic_profile(65), seed=9172761, held_out_count=193)

    assert metrics["stress_disagreement_normalized_rms"] < 5.0e-3
    # The independent RHS uses shape-preserving PCHIP while the public stress
    # uses Agent-3's cubic-spline antiderivative.  This synthetic check only
    # calibrates that the cross-operator discrepancy is small enough to expose
    # the 5% mutation; the real-artifact report keeps its own stricter status.
    assert metrics["identity_rms"] < 1.0e-4
    assert metrics["mutation_to_baseline_identity_ratio"] > 20.0
    assert metrics["raw_inverse_outer_tail_without_moment_subtraction"] > 1.0e-8
    assert metrics["public_support_edge_outside_max_abs"] == 0.0
    assert metrics["independent_support_edge_outside_max_abs"] == 0.0


def test_independent_reconstruction_is_not_the_agent3_cubic_spline_path():
    profile = _synthetic_profile(49)
    independent = IndependentPchipRadialInverse(profile)
    assert independent._source.__class__.__name__ == "PchipInterpolator"
    assert independent._bump.__class__.__name__ == "PchipInterpolator"


def test_real_report_requires_three_strictly_increasing_resolutions_before_io():
    with pytest.raises(ValueError, match="at least three"):
        run_real_candidate_preflight(radial_counts=(33, 49))
    with pytest.raises(ValueError, match="strictly increasing"):
        run_real_candidate_preflight(radial_counts=(33, 65, 49))
