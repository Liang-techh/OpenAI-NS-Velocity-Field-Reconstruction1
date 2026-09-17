import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_collar_vorticity import (
    governed_supported_seed,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_cubic_temporal_capacity import (
    audit_supported_phi10_cubic_temporal_capacity,
    cubic_phi10_coefficient,
    cubic_phi10_delta,
    cubic_phi10_delta_dtau,
    cubic_phi10_snapshot,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_quadratic_temporal_capacity import (
    quadratic_phi10_snapshot,
)


def test_cubic_schedule_preserves_early_and_returns_static_at_intermediate():
    base = governed_supported_seed()
    basis = base.parent.profile_basis
    index = basis.mode_indices.index((1, 0))
    midpoint = basis.phi_coefficients[index]

    assert cubic_phi10_delta(-1.0) == pytest.approx(-1.4, abs=1.0e-15)
    assert cubic_phi10_delta(0.0) == 0.0
    assert cubic_phi10_delta(0.5) == 0.0
    assert cubic_phi10_delta(1.0) == 0.0
    assert cubic_phi10_coefficient(base, 0.25) == pytest.approx(midpoint - 1.4)
    assert cubic_phi10_coefficient(base, 0.50) == pytest.approx(midpoint)
    assert cubic_phi10_coefficient(base, 0.625) == pytest.approx(midpoint)
    assert cubic_phi10_coefficient(base, 0.75) == pytest.approx(midpoint)

    probes = np.array(
        [[0.37, 0.11, -0.22], [1.10, -0.25, 0.42], [1.72, 0.0, 0.30]],
        dtype=float,
    )
    early_cubic = cubic_phi10_snapshot(base, 0.25)
    early_quadratic = quadratic_phi10_snapshot(base, 0.25)
    assert np.allclose(
        early_cubic.at_points(probes, 0.25),
        early_quadratic.at_points(probes, 0.25),
        rtol=0.0,
        atol=1.0e-13,
    )
    for time in (0.50, 0.625, 0.75):
        trial = cubic_phi10_snapshot(base, time)
        assert np.array_equal(trial.at_points(probes, time), base.at_points(probes, time))


def test_cubic_capacity_screen_adds_one_conditioned_temporal_degree_and_reduces_collateral():
    report = audit_supported_phi10_cubic_temporal_capacity()

    assert report["schema"] == "eq45_supported_phi10_cubic_temporal_capacity_v1"
    design = report["time_design"]
    assert design["rank"] == 3
    assert design["basis_dimension_increment"] == 1
    assert design["condition_number"] < 6.0

    late = report["late_half_localization"]
    assert late["cubic_to_quadratic_ratio"] < 0.13
    assert late["cubic_max_abs_coefficient_delta"] < 0.023

    rows = {row["time"]: row for row in report["probe_rows"]}
    assert rows[0.25]["cubic_quadratic_public_velocity_max_abs_difference"] < 1.0e-12
    assert rows[0.375]["cubic_to_quadratic_public_velocity_delta_rms_ratio"] < 0.75
    assert rows[0.625]["cubic_public_velocity_delta"]["rms"] == 0.0
    assert rows[0.625]["quadratic_public_velocity_delta"]["rms"] > 0.0
    assert rows[0.6875]["cubic_to_quadratic_public_velocity_delta_rms_ratio"] < 0.20

    for time in ("0.5", "0.625", "0.75"):
        assert report["static_return"][time]["cubic_public_velocity_delta_rms"] == 0.0
        assert report["static_return"][time]["cubic_public_velocity_delta_max_abs"] == 0.0

    sharpness = report["temporal_sharpness"]
    assert 1.4 < sharpness["cubic_to_quadratic_start_derivative_ratio"] < 1.5
    assert report["visualization_fingerprint_implication"][
        "quadratic_known_intermediate_collateral_removed_by_identity"
    ] is True

    truth = report["truth_boundary"]
    assert truth["new_spatial_basis_added"] is False
    assert truth["temporal_basis_dimension_increment"] == 1
    assert truth["production_temporal_shape_promoted"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert report["pde_validation_rerun"] is False


def test_cubic_screen_fails_closed_on_invalid_inputs():
    with pytest.raises(ValueError, match="tau must be finite"):
        cubic_phi10_delta(np.array([0.0, np.nan]))

    with pytest.raises(ValueError, match="probe_times must include"):
        audit_supported_phi10_cubic_temporal_capacity(
            probe_times=(0.25, 0.30, 0.40, 0.50, 0.70, 0.75)
        )

    assert cubic_phi10_delta_dtau(-1.0) == pytest.approx(3.033333333333333, rel=1.0e-14)
