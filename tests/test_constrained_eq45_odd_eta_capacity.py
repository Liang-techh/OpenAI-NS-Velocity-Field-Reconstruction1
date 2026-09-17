import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_odd_eta_capacity import (
    audit_eq45_odd_eta_capacity,
)


ROOT = Path(__file__).resolve().parents[1]
SEED = ROOT / "artifacts/constrained/eq45_velocity_candidate_seed.json"
CALIBRATION = ROOT / "artifacts/constrained/eq45_odd_eta_parity_capacity_seed.json"


def test_existing_odd_eta_modes_give_two_stable_independent_parity_controls():
    candidate = Eq45VelocityCandidate.load_json(SEED)
    original_sha = candidate.sha256
    checked = json.loads(CALIBRATION.read_text(encoding="utf-8"))

    report = audit_eq45_odd_eta_capacity(candidate)

    assert candidate.sha256 == original_sha == checked["candidate_sha256"]
    assert report["parameter_labels"] == ["Phi(0,1)", "F(0,1)"]
    assert report["metric_labels"] == checked["metric_labels"]
    assert np.allclose(
        report["baseline_normalized_parity_fingerprint"],
        checked["baseline_normalized_parity_fingerprint"],
        rtol=0.0,
        atol=5e-15,
    )
    response = np.asarray(report["finest_normalized_response_per_unit_coefficient"])
    checked_response = np.asarray(checked["finest_normalized_response_per_unit_coefficient"])
    assert np.allclose(response, checked_response, rtol=2e-12, atol=2e-13)

    assert report["numerical_rank"] == checked["numerical_rank"] == 2
    assert report["condition_number"] == pytest.approx(checked["condition_number"], rel=2e-12)
    assert report["step_refinement_relative_change"] < 1e-11

    # Phi(0,1) breaks poloidal parity but does not create odd swirl.
    assert abs(response[0, 0]) > 0.5
    assert abs(response[1, 0]) > 0.1
    assert abs(response[2, 0]) < 1e-12
    # F(0,1) is the independent odd-swirl channel and leaves poloidal metrics alone.
    assert abs(response[0, 1]) < 1e-12
    assert abs(response[1, 1]) < 1e-12
    assert abs(response[2, 1]) > 0.1
    assert report["phi_to_odd_swirl_cross_talk_fraction"] < 1e-12
    assert report["swirl_to_poloidal_cross_talk_fraction"] < 1e-12

    truth = report["truth_boundary"]
    assert truth["velocity_changed"] is False
    assert truth["odd_eta_mode_activated_in_production"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_odd_eta_capacity_fails_closed_on_invalid_or_already_active_inputs():
    candidate = Eq45VelocityCandidate.seed()

    with pytest.raises(ValueError, match="strictly decreasing"):
        audit_eq45_odd_eta_capacity(candidate, coefficient_steps=(0.02, 0.04))
    with pytest.raises(ValueError, match="integer >= 4"):
        audit_eq45_odd_eta_capacity(candidate, azimuth_count=3)
    with pytest.raises(ValueError, match="delivery interval"):
        audit_eq45_odd_eta_capacity(candidate, times=(0.24, 0.5, 0.75))
    with pytest.raises(ValueError, match="positive"):
        audit_eq45_odd_eta_capacity(candidate, radii=(0.0, 0.7))

    basis = candidate.profile_basis
    phi = list(basis.phi_coefficients)
    phi[basis.mode_indices.index((0, 1))] = 0.1
    from dataclasses import replace

    active = replace(candidate, profile_basis=replace(basis, phi_coefficients=tuple(phi)))
    with pytest.raises(ValueError, match="dormant"):
        audit_eq45_odd_eta_capacity(active)
