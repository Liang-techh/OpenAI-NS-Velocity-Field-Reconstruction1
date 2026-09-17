from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_parameter_perturbation import (
    _perturbed_candidate,
    run_eq45_parameter_perturbation_audit,
)


ROOT = Path(__file__).resolve().parents[1]
SEED_ARTIFACT = ROOT / "artifacts/constrained/eq45_velocity_candidate_seed.json"


def test_audit_is_deterministic_held_out_and_keeps_truth_boundary_false():
    first = run_eq45_parameter_perturbation_audit(SEED_ARTIFACT)
    second = run_eq45_parameter_perturbation_audit(SEED_ARTIFACT)

    assert first == second
    assert first["probe_seed"] == 914131
    assert first["probe_count"] == 128
    assert first["times"] == [0.25, 0.5, 0.75]
    assert first["perturbation_seeds"] == [271828, 314159, 161803]
    assert len(first["responses"]) == 18
    assert len(first["perturbed_candidate_sha256"]) == 6

    truth = first["truth_boundary"]
    assert truth["velocity_changed_in_repository"] is False
    assert truth["generalization_diagnostic_only"] is True
    assert truth["physical_support_validated"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_small_unseen_perturbations_have_finite_bounded_public_velocity_response():
    report = run_eq45_parameter_perturbation_audit(SEED_ARTIFACT)
    summary = report["low_scale_summary"]

    assert 0.004 < summary["rms_relative_change_min"] < 0.006
    assert 0.015 < summary["rms_relative_change_max"] < 0.018
    assert 0.007 < summary["rms_relative_change_mean"] < 0.009
    assert 0.30 < summary["response_gain_min"] < 0.40
    assert 1.10 < summary["response_gain_max"] < 1.20

    for row in report["responses"]:
        assert np.isfinite(row["rms_relative_velocity_change"])
        assert np.isfinite(row["max_change_over_base_rms"])
        assert np.isfinite(row["response_gain"])
        assert row["rms_relative_velocity_change"] > 0.0
        assert row["max_change_over_base_rms"] > 0.0
        assert row["response_gain"] > 0.0


def test_same_direction_double_scale_doubles_response_and_stays_in_bounds():
    report = run_eq45_parameter_perturbation_audit(SEED_ARTIFACT)
    for calibration in report["scale_calibration"].values():
        assert calibration["high_over_low_response_ratio_min"] == pytest.approx(
            2.0, rel=1e-10, abs=1e-12
        )
        assert calibration["high_over_low_response_ratio_max"] == pytest.approx(
            2.0, rel=1e-10, abs=1e-12
        )

    base = Eq45VelocityCandidate.load_json(SEED_ARTIFACT)
    for seed in report["perturbation_seeds"]:
        changed, relative_l2 = _perturbed_candidate(base, seed=seed, scale=0.010)
        coefficients = np.asarray(
            changed.profile_basis.phi_coefficients
            + changed.profile_basis.swirl_coefficients
        )
        assert np.max(np.abs(coefficients)) <= changed.profile_basis.coefficient_limit
        assert relative_l2 == pytest.approx(0.02851692023198076)


def test_audit_fails_closed_on_invalid_sampling_or_out_of_bounds_perturbation():
    with pytest.raises(ValueError, match="probe_count"):
        run_eq45_parameter_perturbation_audit(SEED_ARTIFACT, probe_count=8)
    with pytest.raises(ValueError, match="three distinct"):
        run_eq45_parameter_perturbation_audit(
            SEED_ARTIFACT, perturbation_seeds=(1, 1, 2)
        )
    with pytest.raises(ValueError, match="candidate delivery interval"):
        run_eq45_parameter_perturbation_audit(
            SEED_ARTIFACT, times=(0.25, 0.5, 0.8)
        )

    base = Eq45VelocityCandidate.load_json(SEED_ARTIFACT)
    with pytest.raises(ValueError, match="coefficient bounds"):
        _perturbed_candidate(base, seed=271828, scale=2.0)
