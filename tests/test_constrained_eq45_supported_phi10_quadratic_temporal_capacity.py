import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_collar_vorticity import (
    governed_supported_seed,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_quadratic_temporal_capacity import (
    audit_supported_phi10_early_localized_temporal_capacity,
    quadratic_phi10_coefficient,
    quadratic_phi10_delta,
    quadratic_phi10_snapshot,
)


def test_quadratic_schedule_keeps_early_and_returns_to_static_endpoints():
    base = governed_supported_seed()
    basis = base.parent.profile_basis
    index = basis.mode_indices.index((1, 0))
    midpoint = basis.phi_coefficients[index]

    assert quadratic_phi10_delta(-1.0) == pytest.approx(-1.4, abs=1.0e-15)
    assert quadratic_phi10_delta(0.0) == 0.0
    assert quadratic_phi10_delta(1.0) == 0.0
    assert quadratic_phi10_coefficient(base, 0.25) == pytest.approx(midpoint - 1.4)
    assert quadratic_phi10_coefficient(base, 0.50) == pytest.approx(midpoint)
    assert quadratic_phi10_coefficient(base, 0.75) == pytest.approx(midpoint)

    probes = np.array(
        [[0.37, 0.11, -0.22], [1.10, -0.25, 0.42], [1.72, 0.0, 0.30]],
        dtype=float,
    )
    for time in (0.50, 0.75):
        trial = quadratic_phi10_snapshot(base, time)
        assert np.array_equal(trial.at_points(probes, time), base.at_points(probes, time))


def test_early_localized_temporal_capacity_screen_is_resolved_and_truth_bounded():
    report = audit_supported_phi10_early_localized_temporal_capacity(
        resolutions=(25, 33, 41),
    )

    assert report["schema"] == "eq45_supported_phi10_early_localized_temporal_capacity_v1"
    assert report["time_design"]["rank"] == 2
    assert report["time_design"]["basis_dimension_increment"] == 1
    assert report["time_design"]["condition_number"] < 1.2

    schedule = report["quadratic_schedule"]
    assert schedule["start_delta"] == pytest.approx(-1.4, abs=1.0e-15)
    assert schedule["midpoint_delta"] == 0.0
    assert schedule["end_delta"] == 0.0
    assert schedule["extrema"]["minimum"] >= -schedule["extrema"]["coefficient_limit"]
    assert schedule["extrema"]["maximum"] <= schedule["extrema"]["coefficient_limit"]

    summaries = {row["time"]: row for row in report["fine_reference_summaries"]}
    early = summaries[0.25]
    midpoint = summaries[0.50]
    late = summaries[0.75]

    # The quadratic trial is exactly the already-selected affine trial at t=.25,
    # so it must retain the resolved early radial-outer suppression.
    assert early["quadratic_baseline_radial_outer_vorticity2_ratio"] < 0.6
    assert early["quadratic_affine_radial_q99_ratio"] == pytest.approx(1.0, abs=1.0e-12)
    assert early["quadratic_affine_axial_q99_ratio"] == pytest.approx(1.0, abs=1.0e-12)
    assert early["quadratic_affine_radial_outer_vorticity2_ratio"] == pytest.approx(
        1.0, abs=1.0e-12
    )

    # The extra temporal degree returns exactly to the static supported field at
    # both midpoint and late endpoint rather than inheriting affine late drift.
    for row in (midpoint, late):
        assert row["quadratic_baseline_radial_q99_ratio"] == pytest.approx(1.0, abs=1.0e-12)
        assert row["quadratic_baseline_axial_q99_ratio"] == pytest.approx(1.0, abs=1.0e-12)
        assert row["quadratic_baseline_weighted_aspect_ratio"] == pytest.approx(
            1.0, abs=1.0e-12
        )
        assert row["quadratic_baseline_radial_outer_vorticity2_ratio"] == pytest.approx(
            1.0, abs=1.0e-12
        )

    assert all(
        item["exact_public_velocity_match"]
        for item in report["endpoint_equivalence"].values()
    )
    assert max(
        item["max_abs_difference"] for item in report["endpoint_equivalence"].values()
    ) == 0.0

    # At tau=.5, the quadratic coefficient excursion is one quarter of the
    # affine excursion.  Public velocity should reflect the same strong reduction
    # in temporal collateral on deterministic probes.
    intermediate = report["intermediate_collateral"]
    assert intermediate["tau"] == pytest.approx(0.5)
    assert 0.20 < intermediate["quadratic_to_affine_public_velocity_delta_rms_ratio"] < 0.30

    truth = report["truth_boundary"]
    assert truth["screen_trial_velocity_changed"] is True
    assert truth["new_spatial_basis_added"] is False
    assert truth["temporal_basis_dimension_increment"] == 1
    assert truth["production_temporal_shape_promoted"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert report["pde_validation_rerun"] is False


def test_early_localized_screen_fails_closed_on_mismatched_endpoint_amplitude():
    with pytest.raises(ValueError, match="early_delta must equal -affine_slope"):
        audit_supported_phi10_early_localized_temporal_capacity(
            early_delta=-1.0,
            affine_slope=1.4,
            resolutions=(9, 11, 13),
        )

    with pytest.raises(ValueError, match="tau must be finite"):
        quadratic_phi10_delta(np.array([0.0, np.nan]))
