import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_early_localized_temporal_mode import (
    Eq45SupportedPhi10EarlyLocalizedTemporalCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_temporal_mode import (
    Eq45SupportedAffineTemporalModeCandidate,
)


def _base():
    return Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())


def _points():
    return np.array(
        [
            [0.35, 0.10, -0.30],
            [0.60, -0.20, 0.15],
            [1.10, 0.25, 0.45],
            [1.72, 0.00, 0.20],
        ],
        dtype=float,
    )


def test_selected_quadratic_schedule_preserves_screened_endpoint_identities():
    base = _base()
    candidate = Eq45SupportedPhi10EarlyLocalizedTemporalCandidate(base=base)
    affine = Eq45SupportedAffineTemporalModeCandidate(
        base=base, family="phi", mode_i=1, mode_j=0, slope=1.4
    )

    assert candidate.tau(0.25) == pytest.approx(-1.0)
    assert candidate.tau(0.50) == pytest.approx(0.0)
    assert candidate.tau(0.75) == pytest.approx(1.0)
    assert candidate.midpoint_coefficient == pytest.approx(-0.3)
    assert candidate.coefficient_at(0.25) == pytest.approx(-1.7)
    assert candidate.coefficient_at(0.50) == pytest.approx(-0.3)
    assert candidate.coefficient_at(0.625) == pytest.approx(-0.125)
    assert candidate.coefficient_at(0.75) == pytest.approx(-0.3)

    points = _points()
    np.testing.assert_array_equal(
        candidate.at_points(points, 0.25),
        affine.at_points(points, 0.25),
    )
    np.testing.assert_array_equal(
        candidate.at_points(points, 0.50),
        base.at_points(points, 0.50),
    )
    np.testing.assert_array_equal(
        candidate.at_points(points, 0.75),
        base.at_points(points, 0.75),
    )

    # The intermediate late-half excursion is one quarter of the affine one:
    # +0.175 versus +0.7 relative to the static midpoint coefficient.
    quadratic_excursion = float(candidate.coefficient_at(0.625)) - candidate.midpoint_coefficient
    affine_excursion = float(affine.coefficient_at(0.625)) - affine.midpoint_coefficient
    assert quadratic_excursion / affine_excursion == pytest.approx(0.25)


def test_candidate_uses_only_existing_phi10_mode_and_supports_mixed_time_grid():
    base = _base()
    candidate = Eq45SupportedPhi10EarlyLocalizedTemporalCandidate(base=base)

    for time in (0.25, 0.50, 0.625, 0.75):
        snapshot = candidate.snapshot(time)
        source = base.parent.profile_basis
        changed = snapshot.parent.profile_basis
        assert changed.mode_indices == source.mode_indices
        assert changed.radial_degree == source.radial_degree
        assert changed.eta_degree == source.eta_degree
        assert changed.swirl_coefficients == source.swirl_coefficients
        differences = [
            index
            for index, (left, right) in enumerate(
                zip(changed.phi_coefficients, source.phi_coefficients)
            )
            if left != right
        ]
        if time in (0.50, 0.75):
            assert differences == []
        else:
            assert differences == [candidate.mode_index]
        assert snapshot.taper == base.taper

    points = _points()
    times = np.array([0.25, 0.50, 0.625, 0.75])
    values = candidate.at_points(points, times)
    assert values.shape == (4, 3)
    assert np.all(np.isfinite(values))
    for index, time in enumerate(times):
        np.testing.assert_array_equal(
            values[index],
            candidate.snapshot(float(time)).at_points(points[index], float(time)),
        )

    grid = candidate.grid(
        np.array([-0.5, 0.5]),
        np.array([-0.4, 0.4]),
        np.array([-0.3, 0.3]),
        times,
    )
    assert grid.shape == (4, 2, 2, 2, 3)
    assert np.all(np.isfinite(grid))


def test_candidate_roundtrips_with_distinct_identity_and_truth_boundary(tmp_path):
    base = _base()
    candidate = Eq45SupportedPhi10EarlyLocalizedTemporalCandidate(base=base)
    assert candidate.sha256 != base.sha256

    path = tmp_path / "phi10_early_localized.json"
    candidate.save_json(path)
    loaded = Eq45SupportedPhi10EarlyLocalizedTemporalCandidate.load_json(path)
    assert loaded.to_dict() == candidate.to_dict()
    assert loaded.sha256 == candidate.sha256

    times = np.array([0.25, 0.50, 0.625, 0.75])
    points = _points()
    np.testing.assert_array_equal(
        loaded.at_points(points, times),
        candidate.at_points(points, times),
    )

    payload = json.loads(path.read_text(encoding="utf-8"))
    truth = payload["truth_boundary"]
    assert truth["velocity_export_ready"] is True
    assert truth["visualization_candidate_only"] is True
    assert truth["screened_temporal_shape_materialized"] is True
    assert truth["spatial_basis_grown"] is False
    assert truth["canonical_velocity_changed"] is False
    assert truth["production_temporal_shape_promoted"] is False
    for key in (
        "physical_support_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False


def test_candidate_fails_closed_on_zero_nonfinite_bounds_time_and_metadata():
    base = _base()

    with pytest.raises(ValueError, match="nonzero"):
        Eq45SupportedPhi10EarlyLocalizedTemporalCandidate(base=base, early_delta=0.0)
    with pytest.raises(ValueError, match="finite"):
        Eq45SupportedPhi10EarlyLocalizedTemporalCandidate(base=base, early_delta=np.nan)
    with pytest.raises(ValueError, match="coefficient bound"):
        Eq45SupportedPhi10EarlyLocalizedTemporalCandidate(base=base, early_delta=-3.8)
    with pytest.raises(TypeError, match="base"):
        Eq45SupportedPhi10EarlyLocalizedTemporalCandidate(base=object())

    candidate = Eq45SupportedPhi10EarlyLocalizedTemporalCandidate(base=base)
    with pytest.raises(ValueError, match="declared interval"):
        candidate.coefficient_at(0.751)
    with pytest.raises(ValueError, match="snapshot time"):
        candidate.snapshot(np.array([0.25, 0.75]))

    payload = candidate.to_dict()
    payload["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        Eq45SupportedPhi10EarlyLocalizedTemporalCandidate.from_dict(payload)

    payload = candidate.to_dict()
    payload["temporal_mode"]["schedule"] = "different"
    with pytest.raises(ValueError, match="schedule"):
        Eq45SupportedPhi10EarlyLocalizedTemporalCandidate.from_dict(payload)
