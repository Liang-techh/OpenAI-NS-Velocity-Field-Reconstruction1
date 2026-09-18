import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_missing_covariance_column_target import (
    build_missing_covariance_column_target,
)
from openai_ns_reconstruction.kokuno_radial_resolution_covariance_audit import (
    FINE_PAIR_RELATIVE_STABILITY_TOLERANCE,
    RADIAL_RESOLUTION_LADDER,
    REFERENCE_RADIAL_COUNT,
    _validate_radial_counts,
    audit_radial_resolution_ladder,
    compare_nested_receipts,
)


def _synthetic_receipt(count, *, transverse_scale=1.0):
    radii = np.linspace(0.08, 0.36, count)
    current = np.column_stack((radii, 0.25 * radii))
    transverse = np.column_stack((-0.25 * radii, radii))
    target = 0.01 * current + 0.02 * transverse_scale * transverse
    return build_missing_covariance_column_target(radii, target, current)


def test_nested_linear_receipts_are_resolution_identical():
    coarse = _synthetic_receipt(33)
    fine = _synthetic_receipt(65)
    comparison = compare_nested_receipts(coarse, fine, coefficient_budget=0.02)
    assert comparison["coarse_radial_count"] == 33
    assert comparison["fine_radial_count"] == 65
    assert comparison["max_guarded_relative_change"] < 1.0e-12
    assert comparison["current_coefficient_relative_rms_difference"] < 1.0e-12


def test_stable_nested_ladder_passes_frozen_budget_guard():
    receipts = [_synthetic_receipt(count) for count in RADIAL_RESOLUTION_LADDER]
    audit = audit_radial_resolution_ladder(
        receipts,
        coefficient_budget=0.02,
    )
    assert audit["radial_counts"] == list(RADIAL_RESOLUTION_LADDER)
    assert audit["coefficient_budget_refit_on_fine_levels"] is False
    assert audit["radial_resolution_stability_preflight_passed"] is True
    assert audit["finest_pair_max_guarded_relative_change"] < 1.0e-12


def test_fine_target_drift_fails_stability_guard_without_budget_refit():
    receipts = [
        _synthetic_receipt(33),
        _synthetic_receipt(65),
        _synthetic_receipt(129, transverse_scale=1.20),
    ]
    audit = audit_radial_resolution_ladder(
        receipts,
        coefficient_budget=0.02,
    )
    assert audit["coefficient_budget"] == 0.02
    assert audit["coefficient_budget_refit_on_fine_levels"] is False
    assert audit["radial_resolution_stability_preflight_passed"] is False
    assert (
        audit["finest_pair_max_guarded_relative_change"]
        > FINE_PAIR_RELATIVE_STABILITY_TOLERANCE
    )


def test_radial_ladder_must_be_nested_and_increasing():
    assert _validate_radial_counts(RADIAL_RESOLUTION_LADDER) == RADIAL_RESOLUTION_LADDER
    with pytest.raises(ValueError, match="strictly increasing"):
        _validate_radial_counts((33, 33, 65))
    with pytest.raises(ValueError, match="nested"):
        _validate_radial_counts((33, 64, 129))
    with pytest.raises(ValueError, match="at least 9"):
        _validate_radial_counts((7, 13))


def test_default_resolution_contract_is_frozen_before_measurement():
    assert RADIAL_RESOLUTION_LADDER == (33, 65, 129)
    assert REFERENCE_RADIAL_COUNT == 33
    assert FINE_PAIR_RELATIVE_STABILITY_TOLERANCE == 2.0e-2
