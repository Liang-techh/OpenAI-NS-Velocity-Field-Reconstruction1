import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
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


def test_midpoint_is_exact_base_and_endpoints_change_one_existing_mode():
    base = _base()
    candidate = Eq45SupportedAffineTemporalModeCandidate(
        base=base, family="phi", mode_i=1, mode_j=0, slope=0.25
    )

    assert candidate.time_midpoint == 0.5
    assert candidate.tau(0.25) == pytest.approx(-1.0)
    assert candidate.tau(0.50) == pytest.approx(0.0)
    assert candidate.tau(0.75) == pytest.approx(1.0)

    original = base.parent.profile_basis.phi_coefficients[
        base.parent.profile_basis.mode_indices.index((1, 0))
    ]
    assert candidate.coefficient_at(0.25) == pytest.approx(original - 0.25)
    assert candidate.coefficient_at(0.50) == pytest.approx(original)
    assert candidate.coefficient_at(0.75) == pytest.approx(original + 0.25)

    points = _points()
    np.testing.assert_array_equal(
        candidate.at_points(points, 0.50),
        base.at_points(points, 0.50),
    )

    for time in (0.25, 0.75):
        snapshot = candidate.snapshot(time)
        assert snapshot.parent.profile_basis.radial_degree == base.parent.profile_basis.radial_degree
        assert snapshot.parent.profile_basis.eta_degree == base.parent.profile_basis.eta_degree
        assert snapshot.taper == base.taper
        changed = snapshot.parent.profile_basis.phi_coefficients
        source = base.parent.profile_basis.phi_coefficients
        differences = [
            index for index, (left, right) in enumerate(zip(changed, source)) if left != right
        ]
        assert differences == [candidate.mode_index]
        assert not np.allclose(
            candidate.at_points(points, time),
            base.at_points(points, time),
            rtol=0.0,
            atol=1e-13,
        )


def test_swirl_temporal_mode_supports_mixed_time_batches_and_grid():
    base = _base()
    candidate = Eq45SupportedAffineTemporalModeCandidate(
        base=base, family="swirl", mode_i=1, mode_j=2, slope=0.20
    )
    points = _points()[:3]
    times = np.array([0.25, 0.50, 0.75])
    values = candidate.at_points(points, times)
    assert values.shape == (3, 3)
    assert np.all(np.isfinite(values))

    for index, time in enumerate(times):
        np.testing.assert_allclose(
            values[index],
            candidate.snapshot(float(time)).at_points(points[index], float(time)),
            rtol=0.0,
            atol=0.0,
        )

    grid = candidate.grid(
        np.array([-0.5, 0.5]),
        np.array([-0.4, 0.4]),
        np.array([-0.3, 0.3]),
        times,
    )
    assert grid.shape == (3, 2, 2, 2, 3)
    assert np.all(np.isfinite(grid))


def test_temporal_candidate_roundtrips_with_distinct_identity(tmp_path):
    base = _base()
    candidate = Eq45SupportedAffineTemporalModeCandidate(
        base=base, family="phi", mode_i=1, mode_j=2, slope=-0.30
    )
    assert candidate.sha256 != base.sha256

    path = tmp_path / "temporal.json"
    candidate.save_json(path)
    loaded = Eq45SupportedAffineTemporalModeCandidate.load_json(path)
    assert loaded.to_dict() == candidate.to_dict()
    assert loaded.sha256 == candidate.sha256

    times = np.array([0.25, 0.5, 0.75])
    points = np.repeat(_points()[:1], 3, axis=0)
    np.testing.assert_array_equal(
        loaded.at_points(points, times),
        candidate.at_points(points, times),
    )

    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["truth_boundary"]["velocity_export_ready"] is True
    for key in (
        "physical_support_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert payload["truth_boundary"][key] is False
    assert payload["truth_boundary"]["spatial_basis_grown"] is False
    assert payload["truth_boundary"]["affine_temporal_mode_active"] is True


def test_temporal_candidate_fails_closed_on_invalid_or_out_of_bound_requests():
    base = _base()

    with pytest.raises(ValueError, match="nonzero"):
        Eq45SupportedAffineTemporalModeCandidate(
            base=base, family="phi", mode_i=1, mode_j=0, slope=0.0
        )
    with pytest.raises(ValueError, match="family"):
        Eq45SupportedAffineTemporalModeCandidate(
            base=base, family="bad", mode_i=1, mode_j=0, slope=0.1
        )
    with pytest.raises(ValueError, match="already exist"):
        Eq45SupportedAffineTemporalModeCandidate(
            base=base, family="phi", mode_i=0, mode_j=4, slope=0.1
        )
    with pytest.raises(ValueError, match="finite"):
        Eq45SupportedAffineTemporalModeCandidate(
            base=base, family="phi", mode_i=1, mode_j=0, slope=np.nan
        )

    # Phi(0,0) starts at +1 and the existing hard limit is 4, so slope 3.01
    # would reach +4.01 at one registered endpoint.
    with pytest.raises(ValueError, match="coefficient bound"):
        Eq45SupportedAffineTemporalModeCandidate(
            base=base, family="phi", mode_i=0, mode_j=0, slope=3.01
        )

    candidate = Eq45SupportedAffineTemporalModeCandidate(
        base=base, family="phi", mode_i=0, mode_j=0, slope=3.0
    )
    assert candidate.coefficient_at(0.75) == pytest.approx(4.0)
    with pytest.raises(ValueError, match="declared interval"):
        candidate.coefficient_at(0.751)
    with pytest.raises(ValueError, match="snapshot time"):
        candidate.snapshot(np.array([0.25, 0.75]))

    payload = candidate.to_dict()
    payload["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        Eq45SupportedAffineTemporalModeCandidate.from_dict(payload)
