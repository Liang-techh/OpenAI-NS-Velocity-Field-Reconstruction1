from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_physical_center_profile_contract import (
    KokunoPA10PhysicalCenterProfileContract,
)
from openai_ns_reconstruction.kokuno_pa10_source_c_normalized_physical_center import (
    DEFAULT_CERTIFICATE_INTERVALS,
    SELECTED_SOURCE_C,
    KokunoPA10SourceCNormalizedPhysicalCenter,
)


def test_source_complex_C_certificate_closes_with_same_outer_C() -> None:
    profile = KokunoPA10SourceCNormalizedPhysicalCenter()
    certificate = profile.source_C_certificate

    assert profile.C == SELECTED_SOURCE_C == 1000.0
    assert profile.certificate_intervals == DEFAULT_CERTIFICATE_INTERVALS == 2**21
    assert certificate["upstream_complex_neighborhood_machine_certified"] is True
    assert certificate["source_complex_C_normalization_agent1_machine_certified"] is True
    assert certificate["source_complex_C_normalization_independently_admitted"] is False
    assert 0.0 < certificate["phi_star_complex_sup_upper"] < profile.C
    assert certificate["C_over_phi_star_upper"] > 1.0
    assert 0.0 < certificate["normalized_g_complex_sup_upper"] < 1.0
    assert profile.physical_profiles.C == profile.C
    assert profile.axis_domain.binding.outer_schedule.C == profile.C

    truth = profile.truth_boundary
    assert truth["source_complex_C_normalization_agent1_machine_certified"] is True
    assert truth["source_complex_C_normalization_independently_admitted"] is False
    assert truth["selected_C_chosen_from_ns_residual"] is False
    assert truth["selected_C_is_CR001_energy_normalizer"] is False
    assert truth["source_center_is_final_corrected_fixed_point"] is False
    assert truth["cartesian_spacetime_velocity_materialized"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_selected_C_only_rescales_F_center_and_preserves_U_v0() -> None:
    base = KokunoPA10PhysicalCenterProfileContract()
    normalized = KokunoPA10SourceCNormalizedPhysicalCenter().physical_profiles
    assert base.C == 2.0
    assert normalized.C == 1000.0

    X = np.asarray([0.01, 0.08, 0.16, 0.31], dtype=float)
    eta = np.asarray([-0.8, -0.2, 0.35, 0.9], dtype=float)
    before = base.values(X, eta)
    after = normalized.values(X, eta)

    ratio = base.C / normalized.C
    np.testing.assert_allclose(after["F_0"], ratio * before["F_0"], rtol=3e-13, atol=0.0)
    np.testing.assert_allclose(after["E_0"], ratio * before["E_0"], rtol=3e-13, atol=0.0)
    np.testing.assert_allclose(after["U_0"], before["U_0"], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(after["v_0"], before["v_0"], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(after["V_0"], before["V_0"], rtol=0.0, atol=0.0)
    assert np.any(np.abs(after["F_0"]) > 0.0)
    assert np.any(np.abs(after["U_0"]) > 0.0)


def test_source_normalized_configuration_roundtrip(tmp_path) -> None:
    profile = KokunoPA10SourceCNormalizedPhysicalCenter()
    path = tmp_path / "source_c_normalized_profile.json"
    payload = profile.save_configuration(path)
    loaded = KokunoPA10SourceCNormalizedPhysicalCenter.load_configuration(path)

    assert json.loads(path.read_text()) == payload
    assert loaded.configuration() == payload
    assert loaded.C == profile.C == 1000.0
    X = np.asarray([0.02, 0.2], dtype=float)
    eta = np.asarray([-0.5, 0.5], dtype=float)
    before = profile.values(X, eta)
    after = loaded.values(X, eta)
    for key in ("F_0", "E_0", "U_0", "v_0", "V_0"):
        np.testing.assert_allclose(after[key], before[key], rtol=0.0, atol=0.0)


def test_contract_rejects_uncertified_old_C() -> None:
    with pytest.raises(ValueError, match="requires outer/profile C=1000"):
        KokunoPA10SourceCNormalizedPhysicalCenter(
            physical_profiles=KokunoPA10PhysicalCenterProfileContract()
        )


def test_report_keeps_downstream_scientific_gates_fail_closed() -> None:
    profile = KokunoPA10SourceCNormalizedPhysicalCenter()
    report = profile.report()
    checks = report["machine_checks"]
    truth = report["truth_boundary"]

    assert checks["source_C_condition_closed_at_agent1_self_certificate_level"] is True
    assert checks["normalized_g_complex_sup_upper_lt_one"] is True
    assert checks["same_C_propagated_into_outer_schedule"] is True
    assert checks["physical_center_F_nontrivial_on_probe"] is True
    assert checks["physical_center_U_nontrivial_on_probe"] is True
    assert checks["all_physical_center_values_finite"] is True
    assert report["normalization_contract"]["CR001_energy_normalization_used_as_source_C_proof"] is False
    assert report["normalization_contract"]["source_C_used_as_CR001_energy_normalizer"] is False
    assert report["scientific_gates"]["momentum_max_l2"] == 1.0e-3
    assert report["scientific_gates"]["divergence_max_l2"] == 1.0e-5
    assert report["scientific_gates"]["free_residual_defined_forcing_forbidden"] is True
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["complete_kokuno_composite_velocity"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
