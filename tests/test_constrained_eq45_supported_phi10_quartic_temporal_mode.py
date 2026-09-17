import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_cubic_temporal_mode import (
    Eq45SupportedPhi10CubicLocalizedTemporalCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_quartic_derivative_capacity import (
    balanced_nullspace_coefficient,
    quartic_phi10_delta,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_quartic_temporal_mode import (
    Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
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
            [1.68, 0.18, 1.72],
            [-0.75, 0.55, -1.70],
        ],
        dtype=float,
    )


def test_quartic_candidate_materializes_screened_anchor_identities():
    base = _base()
    candidate = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(base=base)
    cubic = Eq45SupportedPhi10CubicLocalizedTemporalCandidate(base=base)

    assert candidate.nullspace_coefficient == pytest.approx(
        -0.2831460674157303, rel=0.0, abs=1.0e-14
    )
    assert candidate.tau(0.25) == pytest.approx(-1.0)
    assert candidate.tau(0.50) == pytest.approx(0.0)
    assert candidate.tau(0.625) == pytest.approx(0.5)
    assert candidate.tau(0.75) == pytest.approx(1.0)
    assert candidate.midpoint_coefficient == pytest.approx(-0.3)
    assert candidate.coefficient_at(0.25) == pytest.approx(-1.7)
    assert candidate.coefficient_at(0.50) == pytest.approx(-0.3)
    assert candidate.coefficient_at(0.625) == pytest.approx(-0.3)
    assert candidate.coefficient_at(0.75) == pytest.approx(-0.3)

    points = _points()
    np.testing.assert_array_equal(
        candidate.at_points(points, 0.25),
        cubic.at_points(points, 0.25),
    )
    for time in (0.50, 0.625, 0.75):
        np.testing.assert_array_equal(
            candidate.at_points(points, time),
            base.at_points(points, time),
        )

    for time in (0.375, 0.5625, 0.6875):
        assert not np.array_equal(
            candidate.at_points(points, time),
            base.at_points(points, time),
        )
        tau = candidate.tau(time)
        assert candidate.coefficient_at(time) == pytest.approx(
            candidate.midpoint_coefficient
            + float(
                quartic_phi10_delta(
                    tau,
                    early_delta=candidate.early_delta,
                    nullspace_coefficient=candidate.nullspace_coefficient,
                )
            )
        )


def test_quartic_candidate_changes_only_existing_phi10_and_supports_mixed_time_grid():
    base = _base()
    candidate = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(base=base)

    for time in (0.25, 0.375, 0.50, 0.5625, 0.625, 0.6875, 0.75):
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
        if time in (0.50, 0.625, 0.75):
            assert differences == []
        else:
            assert differences == [candidate.mode_index]
        assert snapshot.taper == base.taper

    points = _points()
    times = np.array([0.25, 0.375, 0.50, 0.5625, 0.625, 0.6875])
    values = candidate.at_points(points, times)
    assert values.shape == (6, 3)
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
        np.array([0.25, 0.50, 0.625, 0.75]),
    )
    assert grid.shape == (4, 2, 2, 2, 3)
    assert np.all(np.isfinite(grid))


def test_quartic_candidate_roundtrips_with_pinned_balancing_rule_and_truth_boundary(tmp_path):
    base = _base()
    candidate = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(base=base)
    assert candidate.sha256 != base.sha256
    assert candidate.nullspace_coefficient == balanced_nullspace_coefficient(
        early_delta=candidate.early_delta
    )

    path = tmp_path / "phi10_quartic_balanced.json"
    candidate.save_json(path)
    loaded = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate.load_json(path)
    assert loaded.to_dict() == candidate.to_dict()
    assert loaded.sha256 == candidate.sha256

    times = np.array([0.25, 0.375, 0.50, 0.5625, 0.625, 0.6875])
    points = _points()
    np.testing.assert_array_equal(
        loaded.at_points(points, times),
        candidate.at_points(points, times),
    )

    payload = json.loads(path.read_text(encoding="utf-8"))
    temporal = payload["temporal_mode"]
    assert temporal["nullspace_coefficient"] == candidate.nullspace_coefficient
    assert temporal["nullspace_rule"] == "least_squares_minimize_return_node_coefficient_slope"

    truth = payload["truth_boundary"]
    assert truth["velocity_export_ready"] is True
    assert truth["visualization_candidate_only"] is True
    assert truth["screened_temporal_shape_materialized"] is True
    assert truth["derivative_balanced_temporal_shape"] is True
    assert truth["spatial_basis_grown"] is False
    assert truth["force_or_pressure_fitted"] is False
    assert truth["pde_objective_used_to_choose_nullspace_coefficient"] is False
    assert truth["public_image_fitted"] is False
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


def test_quartic_candidate_fails_closed_on_bounds_time_and_metadata():
    base = _base()

    with pytest.raises(ValueError, match="nonzero"):
        Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(
            base=base, early_delta=0.0
        )
    with pytest.raises(ValueError, match="finite"):
        Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(
            base=base, early_delta=np.nan
        )
    with pytest.raises(ValueError, match="coefficient bound"):
        Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(
            base=base, early_delta=-3.8
        )
    with pytest.raises(TypeError, match="base"):
        Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(base=object())

    candidate = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(base=base)
    with pytest.raises(ValueError, match="declared interval"):
        candidate.coefficient_at(0.751)
    with pytest.raises(ValueError, match="snapshot time"):
        candidate.snapshot(np.array([0.25, 0.75]))

    payload = candidate.to_dict()
    payload["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate.from_dict(payload)

    payload = candidate.to_dict()
    payload["temporal_mode"]["nullspace_coefficient"] += 1.0e-6
    with pytest.raises(ValueError, match="nullspace coefficient"):
        Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate.from_dict(payload)

    payload = candidate.to_dict()
    payload["temporal_mode"]["schedule"] = "different"
    with pytest.raises(ValueError, match="schedule"):
        Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate.from_dict(payload)
