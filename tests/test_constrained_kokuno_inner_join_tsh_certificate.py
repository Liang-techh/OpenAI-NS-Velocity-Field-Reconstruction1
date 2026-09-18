import copy
import hashlib
import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_inner_join_tsh_certificate import (
    LOG_F_SUP,
    SIGMA_PRIME_SUP,
    KokunoInnerJoinTshCertificate,
)


@pytest.fixture(scope="module")
def selected_certificate():
    return KokunoInnerJoinTshCertificate(eta_nodes=9, envelope_relative_padding=0.05)


def test_selected_B0_envelope_and_PA10_formula_are_executable_and_fail_closed(selected_certificate):
    cert = selected_certificate
    envelope = cert.envelope_report()
    assert envelope["disjoint_holdout_within_envelope"] is True
    assert envelope["selected_numerical_B0_envelope"] > 0.0
    assert envelope["holdout_max_abs_ell_i"] <= envelope["selected_numerical_B0_envelope"]

    expected = 20.0 * SIGMA_PRIME_SUP * (
        envelope["selected_numerical_B0_envelope"] + LOG_F_SUP
    )
    assert cert.T_sh_lower_bound == pytest.approx(expected, rel=2e-15)
    assert cert.selected_T_sh >= cert.T_sh_lower_bound

    truth = cert.report()["truth_boundary"]
    assert truth["PA10_T_sh_formula_executable"] is True
    assert truth["selected_B0_C0_envelope_numerically_checked"] is True
    assert truth["selected_T_sh_lower_bound_numerically_instantiated"] is True
    assert truth["source_B0_analytic_bound_proved"] is False
    assert truth["source_T_sh_lower_bound_verified"] is False
    assert truth["source_outer_pressure_datum_bound"] is False
    assert truth["inner_to_outer_join_completed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_selected_separation_geometry_is_reported_without_relaxing_source_boundary(selected_certificate):
    cert = selected_certificate
    report = cert.geometry_report()
    assert report["required_log_x_sep_upper_bound"] == -8.0
    assert report["selected_log_x_sep"] == pytest.approx(
        report["log_x_i"] + report["selected_T_sh"], rel=0.0, abs=1e-12
    )
    assert report["max_T_sh_for_current_outer_schedule"] == pytest.approx(
        -8.0 - report["log_x_i"], rel=0.0, abs=1e-12
    )
    if report["separation_geometry_feasible"]:
        join = cert.build_selected_inner_join(quadrature_points=32)
        assert join.T_sh == cert.selected_T_sh
        assert join.log_x_sep < -8.0
        assert cert.required_C_multiplier_if_only_X_R_is_enlarged == 1.0
    else:
        assert report["selected_log_x_sep"] >= -8.0
        assert cert.required_C_multiplier_if_only_X_R_is_enlarged > 1.0
        with pytest.raises(ValueError, match="does not fit"):
            cert.build_selected_inner_join(quadrature_points=32)


def test_selected_Tsh_report_preserves_independent_pressure_blocker(selected_certificate):
    report = selected_certificate.report()
    pressure = report["pressure_datum"]
    assert pressure["source_outer_pressure_datum_bound"] is False
    assert pressure["actual_source_incoming_five_moment_discrepancy_bound"] is False
    assert pressure["selected_reference_pressure_scale"] == pytest.approx(1.0)
    assert pressure["required_source_lower_bound_for_reference_ansatz"] > 1e10
    assert report["source_route_ready"] is False


def test_payload_roundtrip_and_truth_tamper_fail_closed(tmp_path, selected_certificate):
    cert = selected_certificate
    path = cert.save_json(tmp_path / "tsh.json")
    replay = KokunoInnerJoinTshCertificate.load_json(path)
    assert replay.to_payload() == cert.to_payload()
    assert replay.sha256 == cert.sha256

    payload = copy.deepcopy(cert.to_payload())
    payload["truth_boundary"]["source_T_sh_lower_bound_verified"] = True
    unsigned = {key: value for key, value in payload.items() if key != "sha256"}
    payload["sha256"] = hashlib.sha256(
        json.dumps(unsigned, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoInnerJoinTshCertificate.from_payload(payload)


def test_source_constants_match_closed_form_norms():
    assert SIGMA_PRIME_SUP == 8.0
    assert LOG_F_SUP == pytest.approx(math.log(2.0), rel=0.0, abs=0.0)
    eta = np.linspace(-1.0, 1.0, 1001)
    f = 1.0 / (1.0 + eta * eta)
    assert np.max(np.abs(np.log(f))) == pytest.approx(LOG_F_SUP, rel=1e-15)
