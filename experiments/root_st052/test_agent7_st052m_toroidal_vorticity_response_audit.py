from __future__ import annotations

import numpy as np

import agent7_st052m_toroidal_vorticity_response_audit as audit


def test_cartesian_curl_operator_on_solid_rotation() -> None:
    points = audit._probe_points()[:32]

    def solid_rotation(pts: np.ndarray) -> np.ndarray:
        pts = np.asarray(pts, dtype=float)
        return np.column_stack((-pts[:, 1], pts[:, 0], np.zeros(len(pts))))

    omega = audit.curl_fd(solid_rotation, points)
    expected = np.zeros_like(omega)
    expected[:, 2] = 2.0
    np.testing.assert_allclose(omega, expected, atol=2.0e-10, rtol=0.0)


def test_symbolic_identity_and_complementary_unit_vorticity_channels() -> None:
    symbolic = audit.symbolic_toroidal_curl_identity()
    assert symbolic["exact_identity"] is True
    assert symbolic["residual"] == ["0", "0", "0"]

    result = audit.unit_vorticity_selectivity()
    assert result["passes"] is True
    assert result["selectivity_passes"] is True
    assert result["toroidal"]["omega_theta_max_abs"] <= audit.LEAKAGE_GATE
    assert result["toroidal"]["omega_rz_rms"] > audit.NONTRIVIAL_RMS_MIN
    assert result["poloidal"]["omega_rz_max_abs"] <= audit.LEAKAGE_GATE
    assert result["poloidal"]["omega_theta_rms"] > audit.NONTRIVIAL_RMS_MIN
    ident = result["normalized_response"]
    assert ident["rank"] == 2
    assert abs(ident["cosine"]) <= audit.COSINE_MAX_ABS
    assert ident["condition"] <= audit.CONDITION_MAX


def test_actual_frozen_846_child_realizes_the_toroidal_vorticity_increment() -> None:
    result = audit.actual_child_vorticity_response()
    assert result["passes"] is True
    assert result["beta_calibration"]["nontrivial"] is True
    assert result["beta_calibration"]["fresh_714_data_used"] is False
    consistency = result["relative_delta_consistency"]
    assert consistency["passes"] is True
    assert consistency["relative_rms"] <= audit.DELTA_REL_RMS_MAX
    assert consistency["relative_sampled_max"] <= audit.DELTA_REL_SAMPLED_MAX
    assert result["source_numeric_vorticity_target_used"] is False
    assert result["held_out_pde_residual_evaluated"] is False

    assert audit.TRUTH["candidate_velocity_changed"] is False
    assert audit.TRUTH["new_basis_added"] is False
    assert audit.TRUTH["new_coefficient_selected"] is False
    assert audit.TRUTH["held_out_pde_residual_evaluated"] is False
    assert audit.TRUTH["visual_correspondence_verified"] is False
    assert audit.TRUTH["pde_validated"] is False
    assert audit.TRUTH["paper_exact"] is False
    assert audit.TRUTH["openai_field_identified"] is False
