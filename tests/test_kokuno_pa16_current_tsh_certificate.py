import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_tsh_certificate import (
    LOG_F_SUP,
    LOG_X_RESTORE_START,
    SIGMA_PRIME_SUP,
    KokunoPA16CurrentTshCertificate,
)


@pytest.fixture(scope="module")
def certificate():
    return KokunoPA16CurrentTshCertificate(eta_nodes=17)


def test_current_lineage_covers_source_eta_and_envelope_has_disjoint_holdout(certificate):
    lo, hi = certificate.moments.eta_interval
    assert lo <= -1.0
    assert hi >= 1.0

    report = certificate.envelope_report()
    assert report["eta_interval"] == [-1.0, 1.0]
    assert report["eta_construction_nodes"] == 17
    assert report["eta_holdout_nodes"] == 16
    assert report["construction_max_abs_ell_i"] > 0.0
    assert report["holdout_max_abs_ell_i"] > 0.0
    assert report["selected_numerical_B0_envelope"] > report["construction_max_abs_ell_i"]
    assert report["disjoint_holdout_within_envelope"] is True
    assert report["analytic_source_B0_bound_proved"] is False


def test_source_tsh_formula_and_separation_geometry_are_replayed_exactly(certificate):
    expected_lower = 20.0 * SIGMA_PRIME_SUP * (certificate.selected_B0 + LOG_F_SUP)
    assert certificate.T_sh_lower_bound == pytest.approx(expected_lower, rel=0, abs=0)
    assert certificate.selected_T_sh > certificate.T_sh_lower_bound

    schedule = certificate.moments.outer_schedule
    expected_log_x_i = math.log(certificate.moments.X_i) - schedule.log_X_R
    expected_max = LOG_X_RESTORE_START - expected_log_x_i
    assert certificate.log_x_i == pytest.approx(expected_log_x_i, rel=0, abs=2e-14)
    assert certificate.max_T_sh_for_separation == pytest.approx(expected_max, rel=0, abs=2e-14)
    assert certificate.geometry_margin == pytest.approx(
        certificate.max_T_sh_for_separation - certificate.selected_T_sh,
        rel=0,
        abs=2e-14,
    )
    assert certificate.separation_geometry_feasible is (
        certificate.log_x_i + certificate.selected_T_sh < LOG_X_RESTORE_START
    )

    # Since X_R=110(CP_*)^10 and X_i=110, this is an independent algebraic replay.
    independent_max = 10.0 * (math.log(schedule.C) + schedule.log_P_star) - 8.0
    assert certificate.max_T_sh_for_separation == pytest.approx(
        independent_max, rel=0, abs=5e-13
    )


def test_current_xi_ell_values_are_not_a_frozen_or_zero_surrogate(certificate):
    eta = np.asarray([-0.75, -0.25, 0.0, 0.25, 0.75])
    direct = np.asarray(certificate.moments.bridge.handoff_at_Xi(eta)["ell_i"])
    sampled = certificate._ell_values(eta)
    np.testing.assert_allclose(sampled, direct, rtol=0, atol=0)
    assert np.all(np.isfinite(sampled))
    assert float(np.max(sampled) - np.min(sampled)) > 1e-8


def test_selected_join_is_fail_closed_or_binds_exact_current_outer_schedule(certificate):
    if certificate.separation_geometry_feasible:
        join = certificate.build_selected_inner_join(quadrature_points=64)
        assert join.T_sh == certificate.selected_T_sh
        assert join.outer_schedule.to_payload() == certificate.moments.outer_schedule.to_payload()
        assert join.log_x_sep < -8.0
    else:
        with pytest.raises(ValueError, match="does not fit"):
            certificate.build_selected_inner_join(quadrature_points=64)
        assert certificate.required_delta_log_X_R > 0.0
        assert certificate.required_log_C_multiplier_if_only_X_R_is_enlarged > 0.0


def test_report_and_configuration_roundtrip_preserve_truth_boundary(certificate, tmp_path):
    report = certificate.report()
    truth = report["truth_boundary"]
    assert report["source_route_ready"] is False
    assert truth["PA10_T_sh_formula_executable"] is True
    assert truth["selected_B0_C0_envelope_numerically_checked"] is True
    assert truth["source_B0_analytic_bound_proved"] is False
    assert truth["source_T_sh_lower_bound_verified"] is False
    assert truth["five_moment_repair_applied"] is False
    assert truth["inner_to_outer_join_completed"] is False
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    path = tmp_path / "current_tsh_certificate.json"
    certificate.save_configuration(path)
    replay = KokunoPA16CurrentTshCertificate.load_configuration(path)
    assert replay.configuration() == certificate.configuration()
    assert replay.semantic_sha256 == certificate.semantic_sha256

    payload = copy.deepcopy(certificate.configuration())
    payload["schema"] = "wrong-schema"
    with pytest.raises(ValueError, match="schema"):
        KokunoPA16CurrentTshCertificate.from_configuration(payload)


def test_parameter_and_source_domain_guards():
    with pytest.raises(ValueError):
        KokunoPA16CurrentTshCertificate(eta_nodes=10)
    with pytest.raises(ValueError):
        KokunoPA16CurrentTshCertificate(envelope_relative_padding=0.0)
    with pytest.raises(ValueError):
        KokunoPA16CurrentTshCertificate(envelope_absolute_padding=0.0)
