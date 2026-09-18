from __future__ import annotations

import copy
import math

import pytest

from openai_ns_reconstruction.kokuno_rescaled_inner_join_tsh_certificate import (
    LOG_F_SUP,
    LOG_X_RESTORE_START,
    SIGMA_PRIME_SUP,
    WITNESS_POLICY,
    KokunoRescaledInnerJoinTshCertificate,
)


@pytest.fixture(scope="module")
def certificate() -> KokunoRescaledInnerJoinTshCertificate:
    return KokunoRescaledInnerJoinTshCertificate()


def test_pointwise_witness_is_shared_c_and_gives_valid_B0_lower_bound(certificate):
    witness = certificate.pointwise_report()
    expected_eta = certificate.incoming.boundary.reference.seed.phase_stationary_eta
    assert witness["witness_policy"] == WITNESS_POLICY
    assert witness["witness_eta"] == pytest.approx(expected_eta, rel=0.0, abs=0.0)
    assert math.isfinite(witness["ell_i_at_witness"])
    assert witness["abs_ell_i_at_witness"] > 0.0
    assert witness["pointwise_B0_lower_bound"] == pytest.approx(
        abs(witness["ell_i_at_witness"]), rel=0.0, abs=0.0
    )
    assert witness["analytic_source_B0_bound_proved"] is False

    geometry = certificate.geometry_report()
    binding = certificate.incoming.binding
    assert geometry["selected_log_C"] == binding.selected_log_C
    assert geometry["selected_log_X_R"] == binding.log_X_R
    expected_lower = 20.0 * SIGMA_PRIME_SUP * (
        witness["pointwise_B0_lower_bound"] + LOG_F_SUP
    )
    assert geometry["pointwise_T_sh_lower_bound"] == pytest.approx(
        expected_lower, rel=0.0, abs=0.0
    )


def test_current_selected_scale_is_obstructed_without_needing_full_B0(certificate):
    geometry = certificate.geometry_report()
    assert math.isfinite(geometry["pointwise_T_sh_lower_bound"])
    assert math.isfinite(geometry["max_T_sh_for_current_shared_C_scale"])
    assert geometry["separation_margin_upper_bound"] < 0.0
    assert geometry["minimum_log_x_sep_from_witness"] >= LOG_X_RESTORE_START
    assert geometry["selected_path_pointwise_obstructed"] is True
    assert certificate.selected_path_pointwise_obstructed is True
    assert certificate.selected_pa16_handoff_allowed is False
    with pytest.raises(ValueError, match="pointwise B0 lower bound"):
        certificate.pa16_inputs_at_eta(0.0)


def test_report_keeps_source_and_pde_claims_fail_closed(certificate):
    report = certificate.report()
    geometry = report["geometry"]
    truth = report["truth_boundary"]
    assert geometry["fixed_B0_C_multiplier_repair_reported"] is False
    assert "recomputing upstream" in geometry["reason_no_fixed_B0_C_repair"]
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
    original = KokunoRescaledInnerJoinTshCertificate()
    payload = original.to_payload()
    replay = KokunoRescaledInnerJoinTshCertificate.from_payload(payload)
    assert replay.to_payload() == payload

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["inner_to_outer_join_completed"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoRescaledInnerJoinTshCertificate.from_payload(tampered)
