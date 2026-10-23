import numpy as np
import pytest

from openai_ns_reconstruction.constrained_residual_localization import (
    audit_residual_localization_resolution,
    diagnose_residual_localization,
)


def _localized_residual(points: np.ndarray, time: float) -> np.ndarray:
    x, y, z = points.T
    r = np.hypot(x, y)
    theta_x = np.divide(-y, r, out=np.zeros_like(r), where=r > 0)
    theta_y = np.divide(x, r, out=np.zeros_like(r), where=r > 0)
    swirl = 0.12 * np.exp(-2.0 * r * r) * (1.0 - 0.1 * time)
    tip = 1.7 * np.exp(-18.0 * (np.abs(z) - 0.86) ** 2) * np.exp(-5.0 * r * r)
    return np.column_stack((swirl * theta_x, swirl * theta_y, tip))


def test_localization_identifies_axial_tip_dominated_residual():
    snap = diagnose_residual_localization(
        _localized_residual,
        time=0.5,
        radial_half_width=1.0,
        axial_half_height=1.0,
        grid_size=40,
    )

    assert snap.cylindrical_component_energy_fraction["axial"] > 0.95
    assert snap.region_energy_fraction["axial_tips"] > 0.75
    assert snap.region_component_energy_fraction["axial_tips"]["axial"] > 0.80
    assert snap.region_energy_fraction["radial_collar"] < 0.01
    assert sum(snap.region_energy_fraction.values()) == pytest.approx(1.0, abs=1e-12)
    assert sum(snap.cylindrical_component_energy_fraction.values()) == pytest.approx(
        1.0, abs=1e-12
    )
    assert snap.pde_validated is False
    assert snap.paper_exact is False


def test_three_level_resolution_audit_is_stable_for_smooth_fixture():
    audit = audit_residual_localization_resolution(
        _localized_residual,
        time=0.5,
        radial_half_width=1.0,
        axial_half_height=1.0,
        grid_sizes=(16, 24, 40),
    )

    assert audit.max_region_fraction_delta_to_finest[0] < 0.05
    assert audit.max_region_fraction_delta_to_finest[1] < 0.03
    assert audit.max_component_fraction_delta_to_finest[0] < 0.01
    assert audit.relative_total_rms_delta_to_finest[0] < 0.03
    assert audit.max_region_fraction_delta_to_finest[-1] == 0.0
    assert audit.max_component_fraction_delta_to_finest[-1] == 0.0
    assert audit.relative_total_rms_delta_to_finest[-1] == 0.0


def test_region_mutation_moves_energy_to_radial_collar():
    def collar_residual(points: np.ndarray, time: float) -> np.ndarray:
        x, y, z = points.T
        r = np.hypot(x, y)
        radial_x = np.divide(x, r, out=np.zeros_like(r), where=r > 0)
        radial_y = np.divide(y, r, out=np.zeros_like(r), where=r > 0)
        amp = np.exp(-120.0 * (r - 0.88) ** 2) * np.exp(-3.0 * z * z)
        return np.column_stack((amp * radial_x, amp * radial_y, np.zeros_like(amp)))

    snap = diagnose_residual_localization(
        collar_residual,
        time=0.5,
        radial_half_width=1.0,
        axial_half_height=1.0,
        grid_size=48,
    )
    assert snap.region_energy_fraction["radial_collar"] > 0.65
    assert snap.cylindrical_component_energy_fraction["radial"] > 0.999999
    assert snap.cylindrical_component_energy_fraction["axial"] == 0.0


def test_fail_closed_on_bad_inputs_and_zero_residual():
    zero = lambda points, time: np.zeros_like(points)
    bad_shape = lambda points, time: np.zeros((points.shape[0], 2))
    nan_residual = lambda points, time: np.full_like(points, np.nan)

    with pytest.raises(ValueError, match="numerically zero"):
        diagnose_residual_localization(
            zero,
            time=0.5,
            radial_half_width=1.0,
            axial_half_height=1.0,
            grid_size=16,
        )
    with pytest.raises(ValueError, match="shape"):
        diagnose_residual_localization(
            bad_shape,
            time=0.5,
            radial_half_width=1.0,
            axial_half_height=1.0,
            grid_size=16,
        )
    with pytest.raises(ValueError, match="non-finite"):
        diagnose_residual_localization(
            nan_residual,
            time=0.5,
            radial_half_width=1.0,
            axial_half_height=1.0,
            grid_size=16,
        )
    with pytest.raises(ValueError, match="at least three"):
        audit_residual_localization_resolution(
            _localized_residual,
            time=0.5,
            radial_half_width=1.0,
            axial_half_height=1.0,
            grid_sizes=(16, 24),
        )
    with pytest.raises(ValueError, match="strictly increasing"):
        audit_residual_localization_resolution(
            _localized_residual,
            time=0.5,
            radial_half_width=1.0,
            axial_half_height=1.0,
            grid_sizes=(16, 16, 32),
        )
    with pytest.raises(ValueError, match="core_radius_fraction"):
        diagnose_residual_localization(
            _localized_residual,
            time=0.5,
            radial_half_width=1.0,
            axial_half_height=1.0,
            grid_size=16,
            core_radius_fraction=0.8,
            collar_start_fraction=0.75,
        )
