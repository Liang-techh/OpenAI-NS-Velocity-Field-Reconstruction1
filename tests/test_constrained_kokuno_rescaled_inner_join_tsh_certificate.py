from __future__ import annotations

import copy
import math

import pytest

from openai_ns_reconstruction.kokuno_rescaled_inner_join_tsh_certificate import (
    LOG_F_SUP,
    LOG_X_RESTORE_START,
    SIGMA_PRIME_SUP,
    KokunoRescaledInnerJoinTshCertificate,
)


@pytest.fixture(scope="module")
def certificate() -> KokunoRescaledInnerJoinTshCertificate:
    return KokunoRescaledInnerJoinTshCertificate(eta_nodes=5)


def test_rescaled_tsh_uses_shared_c_and_disjoint_eta_envelope(certificate):
    envelope = certificate.envelope_report()
    assert envelope["eta_construction_nodes"] == 5
    assert envelope["eta_holdout_nodes"] == 4
    assert envelope["disjoint_holdout_within_envelope"] is True
    assert envelope["analytic_source_B0_bound_proved"] is False

    expected = 20.0 * SIGMA_PRIME_SUP * (certificate.selected_B0 + LOG_F_SUP)
    assert certificate.T_sh_lower_bound == pytest.approx(expected, rel=0.0, abs=0.0)
    assert certificate.selected_T_sh >= certificate.T_sh_lower_bound

    binding = certificate.incoming.binding
    assert certificate.selected_log_x_sep == pytest.approx(
        binding.log_x_sep(certificate.selected_T_sh), rel=0.0, abs=0.0
    )
    assert certificate.geometry_report()["selected_log_C"] == binding.selected_log_C
    assert certificate.geometry_report()["selected_log_X_R"] == binding.log_X_R


def test_current_selected_rescaled_normalization_fails_pa10_separation(certificate):
    assert certificate.selected_B0 > 0.0
    assert math.isfinite(certificate.selected_T_sh)
    assert math.isfinite(certificate.max_T_sh_for_separation)
    assert certificate.separation_geometry_margin < 0.0
    assert certificate.separation_geometry_feasible is False
    assert certificate.selected_log_x_sep > LOG_X_RESTORE_START
    assert certificate.selected_pa16_handoff_allowed is False
    with pytest.raises(ValueError, match="PA.16 handoff remains closed"):
        certificate.pa16_inputs_at_eta(0.0)


def test_report_does_not_launder_fixed_b0_c_repair_or_source_proof(certificate):
    report = certificate.report()
    geometry = report["geometry"]
    truth = report["truth_boundary"]
    assert geometry["fixed_B0_C_multiplier_repair_reported"] is False
    assert "recomputing" in geometry["reason_no_fixed_B0_C_repair"]
    assert report["selected_pa16_handoff_allowed"] is False
    assert report["source_route_ready"] is False
    assert truth["coupled_C_repair_not_inferred_from_fixed_B0"] is True
    assert truth["source_B0_analytic_bound_proved"] is False
    assert truth["source_T_sh_lower_bound_verified"] is False
    assert truth["inner_to_outer_join_completed"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_roundtrip_and_truth_metadata_fail_closed():
    original = KokunoRescaledInnerJoinTshCertificate(eta_nodes=7)
    payload = original.to_payload()
    replay = KokunoRescaledInnerJoinTshCertificate.from_payload(payload)
    assert replay.to_payload() == payload

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["inner_to_outer_join_completed"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoRescaledInnerJoinTshCertificate.from_payload(tampered)
