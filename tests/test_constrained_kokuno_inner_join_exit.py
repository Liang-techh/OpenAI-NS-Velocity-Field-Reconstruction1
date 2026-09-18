import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_inner_join_exit import (
    KokunoInnerJoinExit,
    LOG_JOIN_EXIT,
    LOG_REPAIR_START,
    LOG_RESTORE_END,
    LOG_RESTORE_START,
)


def _ideal_boundary(join: KokunoInnerJoinExit, eta: float) -> tuple[float, float]:
    f = join.source_f(eta)
    return math.log(f), 4.0 * eta


def test_source_geometry_binds_XR_Xi_Xsep_and_Xh():
    join = KokunoInnerJoinExit(T_sh=20.0)
    schedule = join.outer_schedule
    assert join.log_x_i == pytest.approx(
        -10.0 * (math.log(schedule.C) + schedule.log_P_star), abs=1e-14
    )
    assert join.log_X_i == pytest.approx(math.log(110.0), abs=1e-15)
    assert join.log_X_sep == pytest.approx(math.log(110.0) + 20.0, abs=1e-14)
    assert join.log_x_sep < LOG_RESTORE_START
    assert join.log_X_h == pytest.approx(schedule.log_X_R - 5.0, abs=1e-14)
    report = join.geometry_report()
    assert report["source_T_sh_lower_bound_verified"] is False
    assert report["PA16_repair_log_x"] == [LOG_REPAIR_START, LOG_JOIN_EXIT]


def test_source_transition_and_axial_restore_match_ideal_with_flat_endpoints():
    join = KokunoInnerJoinExit(T_sh=20.0)
    eta = 0.31
    ell_i, G_i = _ideal_boundary(join, eta)
    logs = np.asarray(
        [
            join.log_x_i,
            join.log_x_i + 0.37 * join.T_sh,
            join.log_x_sep,
            LOG_RESTORE_START,
            -7.5,
            LOG_RESTORE_END,
            LOG_REPAIR_START,
            LOG_JOIN_EXIT,
        ]
    )
    x = np.exp(logs)
    values = join.profile_values_scaled(
        x,
        eta=eta,
        ell_i=ell_i,
        G_i=G_i,
        coefficients=np.zeros(5),
    )
    U0, E0 = join._ideal(logs, eta)
    np.testing.assert_allclose(values["U"], U0, rtol=0.0, atol=2e-15)
    np.testing.assert_allclose(values["E"], E0, rtol=3e-14, atol=0.0)
    np.testing.assert_allclose(values["U_x"], 0.0, rtol=0.0, atol=2e-14)
    np.testing.assert_allclose(values["E_x"], 0.1 * E0 / x, rtol=2e-12, atol=0.0)


def test_nonideal_transition_radial_derivative_matches_centered_log_difference():
    join = KokunoInnerJoinExit(T_sh=20.0)
    eta = -0.27
    f = join.source_f(eta)
    ell_i = math.log(f) + 2.5e-4
    G_i = 4.0 * eta + 1.0e-5
    log_x = join.log_x_i + 0.41 * join.T_sh
    x = math.exp(log_x)
    analytic = join.profile_values_scaled(
        np.asarray(x), eta=eta, ell_i=ell_i, G_i=G_i
    )
    step = 2.0e-6
    plus = join.profile_values_scaled(
        np.asarray(math.exp(log_x + step)), eta=eta, ell_i=ell_i, G_i=G_i
    )["E"]
    minus = join.profile_values_scaled(
        np.asarray(math.exp(log_x - step)), eta=eta, ell_i=ell_i, G_i=G_i
    )["E"]
    dE_dlogx = float((plus - minus) / (2.0 * step))
    assert x * float(analytic["E_x"]) == pytest.approx(dE_dlogx, rel=2e-7, abs=1e-11)


def test_ideal_upstream_data_produce_zero_discrepancy_and_zero_repair():
    join = KokunoInnerJoinExit(T_sh=20.0)
    eta = 0.2
    ell_i, G_i = _ideal_boundary(join, eta)
    pre = join.pre_repair_scaled_discrepancy(
        eta=eta,
        ell_i=ell_i,
        G_i=G_i,
        incoming_scaled_discrepancy=np.zeros(5),
    )
    assert np.max(np.abs(pre)) < 5e-13
    solved = join.solve_at_eta(
        eta=eta,
        ell_i=ell_i,
        G_i=G_i,
        incoming_scaled_discrepancy=np.zeros(5),
    )
    assert solved.max_abs_exit_discrepancy < 5e-13
    assert np.max(np.abs(np.asarray(solved.repair.coefficients))) < 1e-9


def test_real_transition_discrepancy_routes_into_PA16_and_closes_at_exit():
    join = KokunoInnerJoinExit(T_sh=20.0, quadrature_points=128)
    eta = 0.2
    f = join.source_f(eta)
    # Deliberately small nonideal caller data.  These are a regression probe,
    # not claimed source Appendix-B outputs.
    ell_i = math.log(f) + 1.0e-5
    G_i = 4.0 * eta + 1.0e-7
    incoming = np.asarray([2e-10, -1e-10, 5e-11, -2e-10, 1e-9])
    solved = join.solve_at_eta(
        eta=eta,
        ell_i=ell_i,
        G_i=G_i,
        incoming_scaled_discrepancy=incoming,
    )
    pre = np.asarray(solved.pre_repair_scaled_discrepancy)
    assert np.linalg.norm(pre) > 1e-9
    assert solved.repair.success is True
    assert max(abs(v) for v in solved.repair.coefficients) < join.repair.coefficient_limit
    assert solved.max_abs_exit_discrepancy <= 2.0e-11

    # The corrected profile is executable inside PA.16 and is exactly ideal at
    # the source exit because the compact bumps vanish before log x=-5.
    coeff = np.asarray(solved.repair.coefficients)
    interior_x = math.exp(-5.5)
    interior = join.profile_values_scaled(
        np.asarray(interior_x), eta=eta, ell_i=ell_i, G_i=G_i, coefficients=coeff
    )
    assert np.isfinite(float(interior["E"])) and float(interior["E"]) > 0.0
    exit_x = math.exp(LOG_JOIN_EXIT)
    exit_values = join.profile_values_scaled(
        np.asarray(exit_x), eta=eta, ell_i=ell_i, G_i=G_i, coefficients=coeff
    )
    U0, E0 = join._ideal(np.asarray(LOG_JOIN_EXIT), eta)
    assert float(exit_values["U"]) == pytest.approx(float(U0), abs=2e-15)
    assert float(exit_values["E"]) == pytest.approx(float(E0), rel=2e-14)


def test_repair_interval_fails_closed_without_coefficients():
    join = KokunoInnerJoinExit(T_sh=20.0)
    eta = 0.1
    ell_i, G_i = _ideal_boundary(join, eta)
    with pytest.raises(ValueError, match="PA.16 coefficients"):
        join.profile_values_scaled(
            np.asarray(math.exp(-5.5)), eta=eta, ell_i=ell_i, G_i=G_i
        )


def test_geometry_fails_closed_when_caller_Tsh_overruns_restore_stage():
    with pytest.raises(ValueError, match="x_sep<exp\(-8\)"):
        KokunoInnerJoinExit(T_sh=400.0)


def test_payload_roundtrip_and_truth_tamper(tmp_path):
    join = KokunoInnerJoinExit(T_sh=20.0, quadrature_points=80)
    path = join.save_json(tmp_path / "inner_join.json")
    replay = KokunoInnerJoinExit.load_json(path)
    assert replay.to_payload() == join.to_payload()
    assert replay.sha256 == join.sha256

    payload = copy.deepcopy(join.to_payload())
    payload["truth_boundary"]["actual_upstream_appendix_B_boundary_data_bound"] = True
    unsigned = {key: value for key, value in payload.items() if key != "sha256"}
    import hashlib, json

    payload["sha256"] = hashlib.sha256(
        json.dumps(unsigned, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoInnerJoinExit.from_payload(payload)
