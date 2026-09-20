from __future__ import annotations

import math
from pathlib import Path

import agent7_st052m_reservoir_axial_stretch_audit as audit


def test_exact_physical_axial_sign_root() -> None:
    result = audit.analytic_root_geometry()
    expected_q = (13.0 + 5.0 * math.sqrt(17.0)) / 30.0
    assert result["polynomial_coefficients"] == [-225, 195, 64]
    assert abs(result["physical_q_root"] - expected_q) < 5.0e-14
    assert abs(result["physical_r_root"] - math.sqrt(expected_q)) < 5.0e-14
    assert 0.0 < result["core_cross_section_area_fraction"] < 1.0


def test_live_reservoir_has_core_stretch_and_annular_return() -> None:
    result = audit.live_probe_signs()
    assert result["all_gates_pass"] is True
    assert result["by_radius"]["r=0.60"]["inside_axial_stretch_core"] is True
    assert result["by_radius"]["r=0.90"]["inside_axial_stretch_core"] is True
    assert result["by_radius"]["r=1.20"]["inside_axial_stretch_core"] is False
    assert result["by_radius"]["r=1.50"]["inside_axial_stretch_core"] is False
    for row in result["by_radius"].values():
        assert row["axial_sign_gate"] is True
        assert row["radial_inward_both_z_signs"] is True


def test_axial_shape_changes_sign_once_and_flux_balances() -> None:
    root = audit.EXPECTED_R_ROOT
    eps = 1.0e-5
    assert float(audit.axial_shape_factor(root - eps)) > 0.0
    assert float(audit.axial_shape_factor(root + eps)) < 0.0
    flux = audit.axial_flux_refinement()
    assert flux["passes"] is True
    assert flux["max_abs_finest_balance"] <= audit.FLUX_TOL
    assert flux["max_medium_to_fine_change"] <= audit.REFINEMENT_TOL


def test_peak_ordering() -> None:
    peaks = audit.radial_peak_geometry()
    assert peaks["positive_core_peak_radius"] < audit.EXPECTED_R_ROOT
    assert peaks["negative_annular_peak_radius"] > audit.EXPECTED_R_ROOT
    assert peaks["positive_core_peak_Cz_over_Z"] > 0.0
    assert peaks["negative_annular_peak_Cz_over_Z"] < 0.0


def test_truth_boundary_and_report(tmp_path: Path) -> None:
    report = audit.run(tmp_path / "report.json")
    assert report["task_id"] == "CR003-ST052M-RESERVOIR-AXIAL-STRETCH-AUDIT-111"
    assert report["prereg_issue"] == 784
    assert report["source_parent"] == {
        "pr": 775,
        "head": "93887a59729d22113badf2ae4dac2f7d868e2703",
    }
    decision = report["capacity_conclusion"]
    assert decision["reservoir_channel_has_visible_inward_radial_and_core_axial_stretch_capacity"] is True
    assert decision["second_poloidal_basis_justified_by_this_audit"] is False
    assert decision["retune_775_from_this_audit_allowed"] is False
    for key in (
        "candidate_velocity_changed",
        "basis_dimension_changed",
        "candidate_coefficient_changed",
        "time_law_changed",
        "pressure_or_force_changed",
        "held_out_pde_residual_evaluated",
        "public_image_numeric_target_used",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert report[key] is False
