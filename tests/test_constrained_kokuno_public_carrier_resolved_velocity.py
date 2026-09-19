import json
import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_candidate_velocity import (
    KokunoPublicCandidateOscillatoryVelocity,
)
from openai_ns_reconstruction.kokuno_public_carrier_resolved_velocity import (
    KokunoCarrierResolvedCandidateOscillatoryVelocity,
    carrier_resolution_diagnostics,
    velocity_osc,
)


def test_carrier_resolution_is_source_formula_derived_and_materially_lower_than_parent():
    legacy = KokunoPublicCandidateOscillatoryVelocity()
    repaired = KokunoCarrierResolvedCandidateOscillatoryVelocity()

    old = carrier_resolution_diagnostics(legacy)
    new = carrier_resolution_diagnostics(repaired)

    assert old["source_phase_formula_bound"] is True
    assert new["source_phase_formula_bound"] is True
    assert old["independent_fd4_divergence_assessed"] is False
    assert new["independent_fd4_divergence_assessed"] is False

    # The rejected #539 autonomous background generated a very high angular
    # winding.  This increment changes only autonomous background constants and
    # deliberately lowers that carrier without altering the source phase rule.
    assert old["max_abs_kp"] >= 300
    assert new["max_abs_kp"] <= 8
    assert new["max_local_wavenumber"] < old["max_local_wavenumber"] / 40.0
    assert new["minimum_local_wavelength"] > old["minimum_local_wavelength"] * 40.0

    finest = new["phase_increments_by_step"]["0.005"]
    assert finest["max_local_phase_increment"] < 0.09
    assert finest["max_tangential_phase_increment"] < 0.09


def test_carrier_resolved_public_provider_remains_nontrivial_finite_and_compact():
    points = np.asarray(
        (
            (0.62, 0.13, 0.17, 0.37),
            (0.71, -0.28, -0.21, 0.49),
            (0.88, 0.19, 0.31, 0.61),
        ),
        dtype=float,
    )
    vel = velocity_osc(points[:, 0], points[:, 1], points[:, 2], points[:, 3])
    assert vel.shape == (3, 3)
    assert np.all(np.isfinite(vel))
    assert np.linalg.norm(vel) > 1.0e-8

    exterior = velocity_osc(
        np.asarray((0.0, 1.50, 0.72, 0.72)),
        np.asarray((0.0, 0.0, 0.0, 0.0)),
        np.asarray((0.0, 0.0, 2.05, -2.05)),
        np.asarray((0.41, 0.43, 0.47, 0.53)),
    )
    assert np.array_equal(exterior, np.zeros_like(exterior))


def test_carrier_resolved_payload_records_autonomous_background_without_source_laundering(tmp_path):
    field = KokunoCarrierResolvedCandidateOscillatoryVelocity()
    payload = field.to_payload()
    background = payload["frozen_realization"]["background"]
    receipt = payload["frozen_realization"]["carrier_resolution"]

    assert payload["schema"] == "kokuno-public-carrier-resolved-oscillatory-velocity-v1"
    assert background["R0"] == pytest.approx(0.82)
    assert background["F0"] == pytest.approx(16.0)
    assert background["a"] == pytest.approx(64.0)
    assert background["u_star"] == pytest.approx(8.0)
    assert background["F"] == pytest.approx(16.0)
    assert receipt["background_binding"] == "repository_autonomous_bounded_carrier_resolution_reparameterization"
    assert receipt["source_formulas_changed"] is False
    assert receipt["vector_potential_complete_curl_path_changed"] is False
    assert receipt["independent_agent4_reaudit_required"] is True

    truth = payload["truth_boundary"]
    assert truth["candidate_carrier_resolution_reparameterized"] is True
    assert truth["candidate_carrier_scales_quantified"] is True
    assert truth["independent_fd4_divergence_reaudit_required"] is True
    assert truth["independent_fd4_divergence_passed"] is False
    assert truth["actual_positive_order_background_bound"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False

    path = field.save_json(tmp_path / "carrier-resolved.json")
    raw = json.loads(path.read_text())
    loaded = KokunoCarrierResolvedCandidateOscillatoryVelocity.load_json(path)
    assert raw["sha256"] == field.sha256
    assert loaded.sha256 == field.sha256
    assert loaded.to_payload() == field.to_payload()


def test_diagnostic_rejects_invalid_resolution_inputs():
    field = KokunoCarrierResolvedCandidateOscillatoryVelocity()
    with pytest.raises(ValueError, match="reference_radius"):
        carrier_resolution_diagnostics(field, reference_radius=0.0)
    with pytest.raises(ValueError, match="fd4_steps"):
        carrier_resolution_diagnostics(field, fd4_steps=(0.01, -0.005))
