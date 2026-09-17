from functools import lru_cache

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_cubic_force_crosscheck import (
    compare_cubic_localized_phi10_force,
)


@lru_cache(maxsize=1)
def _report():
    return compare_cubic_localized_phi10_force()


def test_cubic_crosscheck_reuses_frozen_velocity_and_preregistered_force_contract():
    report = _report()

    assert report["velocity_schedule_frozen_before_pde_fit"] is True
    assert report["early_delta"] == pytest.approx(-1.4)
    assert report["mode"] == {"family": "phi", "index": [1, 0]}
    assert report["training_randomness_used"] is False
    assert report["training_random_seed"] is None
    assert report["fit_and_holdout_separate"] is True
    assert report["force_bounds"] == [0.0, 10.0]
    assert report["pressure_fitted"] is False
    assert report["new_force_direction_added"] is False
    assert report["residual_defined_force_allowed"] is False
    assert report["temporal_schedule_refit_on_pde"] is False
    assert all(report["public_identity_checks"].values())

    for name in ("static_projection", "quadratic_projection", "cubic_projection"):
        fit = report[name]["fit"]
        assert 0.0 <= fit["a"] <= 10.0
        assert 0.0 <= fit["c"] <= 10.0
        assert fit["max_iterations"] == 200
        assert fit["solver_tolerance"] == pytest.approx(1.0e-12)
        assert np.isfinite(fit["before"]["rms"])
        assert np.isfinite(fit["after"]["rms"])


def test_cubic_crosscheck_preserves_checked_static_and_quadratic_regression():
    report = _report()
    finest = report["finest_holdout_comparison"]

    # These are already-established #135/#152 results.  Keeping them fixed makes
    # this new report a cross-check of the cubic candidate, not a silent operator
    # or probe-set change.
    assert finest["spatial_step"] == pytest.approx(0.005)
    assert finest["static_zero_force_rms"] == pytest.approx(
        3.5872881930415295, rel=0.0, abs=5.0e-10
    )
    assert finest["static_projected_force_rms"] == pytest.approx(
        3.5860556463766144, rel=0.0, abs=5.0e-10
    )
    assert finest["quadratic_zero_force_rms"] == pytest.approx(
        3.706399725216429, rel=0.0, abs=5.0e-10
    )
    assert finest["quadratic_projected_force_rms"] == pytest.approx(
        3.7050946339151385, rel=0.0, abs=5.0e-10
    )

    assert np.isfinite(finest["cubic_zero_force_rms"])
    assert np.isfinite(finest["cubic_projected_force_rms"])
    assert np.isfinite(finest["cubic_vs_static_zero_force_fractional_change"])
    assert np.isfinite(finest["cubic_vs_quadratic_zero_force_fractional_change"])


def test_cubic_crosscheck_records_sharper_time_derivative_without_promotion():
    report = _report()
    sharpness = report["temporal_sharpness_context"]

    assert sharpness["quadratic_abs_d_delta_dtau_at_start"] == pytest.approx(2.1)
    assert sharpness["cubic_abs_d_delta_dtau_at_start"] == pytest.approx(3.033333333333333)
    assert sharpness["cubic_to_quadratic_start_derivative_ratio"] == pytest.approx(13.0 / 9.0)

    assert report["formal_pde_gate_assessed"] is False
    assert report["visualization_candidate_only"] is True
    assert report["canonical_velocity_changed"] is False
    for key in (
        "visual_correspondence_verified",
        "visualization_ready",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert report[key] is False
