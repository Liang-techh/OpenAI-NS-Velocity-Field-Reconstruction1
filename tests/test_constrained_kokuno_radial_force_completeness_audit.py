import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_radial_force_completeness_audit import (
    FINE_PAIR_RELATIVE_STABILITY_TOLERANCE,
    Z_DERIVATIVE_STEP_LADDER,
    _validate_z_steps,
    audit_radial_force_derivative_ladder,
    centered_radial_force,
)


def test_centered_radial_force_is_exact_for_affine_z_stress():
    radii = np.linspace(0.08, 0.36, 33)
    slope = 0.3 + 0.7 * radii
    base = np.sin(radii)
    step = 0.01
    sigma_plus = base + step * slope
    sigma_minus = base - step * slope
    measured = centered_radial_force(sigma_plus, sigma_minus, step)
    assert np.allclose(measured, slope, rtol=0.0, atol=2.0e-14)


def test_centered_derivative_ladder_detects_stable_second_order_limit():
    radii = np.linspace(0.08, 0.36, 33)
    truth = 1.0 + radii
    steps = Z_DERIVATIVE_STEP_LADDER
    derivatives = [truth + step * step * (0.2 + radii) for step in steps]
    audit = audit_radial_force_derivative_ladder(steps, derivatives)
    assert audit["z_steps"] == list(steps)
    assert audit["radial_force_derivative_stability_preflight_passed"] is True
    assert (
        audit["finest_pair_relative_rms_difference"]
        < FINE_PAIR_RELATIVE_STABILITY_TOLERANCE
    )


def test_derivative_ladder_preserves_a_real_numerical_rejection():
    radii = np.linspace(0.08, 0.36, 33)
    base = 1.0 + radii
    derivatives = [base, base, 1.2 * base]
    audit = audit_radial_force_derivative_ladder(
        Z_DERIVATIVE_STEP_LADDER,
        derivatives,
    )
    assert audit["radial_force_derivative_stability_preflight_passed"] is False
    assert (
        audit["finest_pair_relative_rms_difference"]
        > FINE_PAIR_RELATIVE_STABILITY_TOLERANCE
    )


def test_active_mask_is_applied_before_stability_metric():
    truth = np.ones(7)
    coarse = truth.copy()
    fine = truth.copy()
    coarse[[0, -1]] = 100.0
    fine[[0, -1]] = -100.0
    mask = np.array([False, True, True, True, True, True, False])
    audit = audit_radial_force_derivative_ladder(
        (0.02, 0.01),
        (coarse, fine),
        active_mask=mask,
    )
    assert audit["finest_pair_relative_rms_difference"] == 0.0
    assert audit["radial_force_derivative_stability_preflight_passed"] is True


def test_guards_fail_closed_on_bad_derivative_inputs():
    with pytest.raises(ValueError, match="strictly decreasing"):
        _validate_z_steps((0.01, 0.02))
    with pytest.raises(ValueError, match="positive"):
        centered_radial_force(np.ones(3), np.ones(3), 0.0)
    with pytest.raises(ValueError, match="matching 1-D"):
        centered_radial_force(np.ones(3), np.ones(4), 0.01)
    with pytest.raises(ValueError, match="active_mask"):
        audit_radial_force_derivative_ladder(
            (0.02, 0.01),
            (np.ones(3), np.ones(3)),
            active_mask=np.zeros(3, dtype=bool),
        )


def test_default_radial_force_audit_contract_is_frozen():
    assert Z_DERIVATIVE_STEP_LADDER == (0.02, 0.01, 0.005)
    assert FINE_PAIR_RELATIVE_STABILITY_TOLERANCE == 5.0e-2
