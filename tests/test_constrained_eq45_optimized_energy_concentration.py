import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_optimized_energy_concentration import (
    _weighted_quantile,
    audit_optimized_eq45_energy_concentration,
)


EXPECTED_CANONICAL_SHA = "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"
EXPECTED_OPTIMIZED_SHA = "b34902ee28b3965f245eceb5867df7409498063cb0ea80f1758931723b303c98"


def test_weighted_quantile_is_weighted_and_fail_closed():
    values = np.array([0.0, 1.0, 2.0, 3.0])
    weights = np.array([1.0, 1.0, 8.0, 1.0])
    assert _weighted_quantile(values, weights, 0.5) == 2.0
    assert _weighted_quantile(values, weights, 0.9) == 2.0
    with pytest.raises(ValueError):
        _weighted_quantile(values, -weights, 0.5)
    with pytest.raises(ValueError):
        _weighted_quantile(values, weights, 1.0)


def test_optimized_eq45_energy_concentration_is_axially_rebalanced_and_resolved():
    report = audit_optimized_eq45_energy_concentration()

    assert report["canonical_candidate_sha256"] == EXPECTED_CANONICAL_SHA
    assert report["optimized_candidate_sha256"] == EXPECTED_OPTIMIZED_SHA
    assert report["resolutions"] == [24, 32, 48]
    assert report["times"] == [0.25, 0.5, 0.75]

    expected_totals = {
        0.25: (0.8942276137, 1.1763595243),
        0.50: (0.7335698129, 0.9666628448),
        0.75: (0.5225981760, 0.6893732976),
    }
    for row in report["finest_comparison"]:
        time = row["time"]
        canonical_total, optimized_total = expected_totals[time]
        assert row["canonical"]["total_energy"] == pytest.approx(canonical_total, rel=2e-8)
        assert row["optimized"]["total_energy"] == pytest.approx(optimized_total, rel=2e-8)

        # The #104 update is only mildly broader radially but substantially
        # longer axially in energy support.  This is a morphology diagnostic,
        # not a public-target acceptance condition.
        assert 0.02 < row["radial_rms_relative_change"] < 0.04
        assert 0.35 < row["axial_rms_relative_change"] < 0.37
        assert 0.30 < row["rms_aspect_relative_change"] < 0.34

    # Smooth RMS moments are well resolved from 32^3 -> 48^3.  Quantile metrics
    # are deliberately not gated this tightly because |z| levels are plane
    # quantized; their larger sensitivity remains exposed in the report.
    assert report["max_bulk_rms_relative_delta_penultimate_to_finest"] < 5e-4
    assert report["max_quantile_relative_delta_penultimate_to_finest"] > 1e-2

    assert report["velocity_changed_by_this_audit"] is False
    assert report["physical_support_validated"] is False
    assert report["visualization_ready"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
    assert report["blowup_proved"] is False


def test_energy_concentration_requires_three_resolution_levels():
    with pytest.raises(ValueError, match="at least three"):
        audit_optimized_eq45_energy_concentration(resolutions=(24, 32))
