import copy
import hashlib
import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_appendix_b_boundary import (
    A_FINAL,
    X_FINAL_START,
    X_I,
    KokunoAppendixBBoundary,
)


def test_source_geometry_keeps_two_final_widths_inside_100_to_110():
    boundary = KokunoAppendixBBoundary()
    report = boundary.geometry_report()
    assert report["X_final_start"] == X_FINAL_START
    assert report["X_i"] == X_I
    assert boundary.y100 < boundary.y_axial_end < boundary.y_angular_end < boundary.y_i
    assert report["final_constant_a_log_width"] > 0.0
    assert report["source_free_choices_are_autonomous"] is True


def test_Xi_boundary_values_are_executable_vectorized_and_nontrivial():
    boundary = KokunoAppendixBBoundary(max_step=0.08)
    eta = np.asarray([-0.2, 0.2])
    values = boundary.boundary_values(eta)
    for value in values.values():
        assert value.shape == eta.shape
        assert np.all(np.isfinite(value))
    assert np.all(values["F_i"] > 0.0)
    assert np.all(values["E_i"] > 0.0)
    np.testing.assert_allclose(
        values["ell_i"],
        np.log(boundary.reference.C * values["E_i"]),
        rtol=0.0,
        atol=2e-14,
    )

    # The selected Appendix-B activation is executable rather than silently
    # reusing the frozen reference value at X_i.
    frozen_F = boundary.reference.F(X_I, eta)
    frozen_U = boundary.reference.U(X_I, eta)
    assert np.max(np.abs(values["F_i"] - frozen_F)) > 1e-8
    assert np.max(np.abs(values["G_i"] - frozen_U)) > 1e-10


def test_boundary_solver_refinement_is_stable_at_representative_eta():
    coarse = KokunoAppendixBBoundary(rtol=5e-9, atol=5e-11, max_step=0.08)
    fine = KokunoAppendixBBoundary(rtol=1e-9, atol=1e-11, max_step=0.04)
    a = coarse._solve_eta(0.2)
    b = fine._solve_eta(0.2)
    assert a.ell_i == pytest.approx(b.ell_i, rel=2e-6, abs=2e-8)
    assert a.G_i == pytest.approx(b.G_i, rel=2e-6, abs=2e-8)
    assert b.radial_log_slope_F_i == pytest.approx(-0.5 * A_FINAL, abs=0.0)
    assert b.radial_log_slope_U_i == 0.0


def test_eta_derivative_api_is_finite_and_does_not_claim_source_exactness():
    boundary = KokunoAppendixBBoundary(max_step=0.08)
    deriv = boundary.boundary_eta_derivatives(np.asarray([0.0]))
    assert np.isfinite(float(deriv["ell_i_eta"][0]))
    assert np.isfinite(float(deriv["G_i_eta"][0]))
    truth = boundary.to_payload()["truth_boundary"]
    assert truth["actual_upstream_appendix_B_boundary_values_executable_for_selected_realization"] is True
    assert truth["incoming_five_moment_discrepancy_bound"] is False
    assert truth["source_T_sh_lower_bound_verified"] is False
    assert truth["inner_to_outer_join_completed"] is False
    assert truth["paper_exact"] is False


def test_payload_roundtrip_and_truth_tamper_fail_closed(tmp_path):
    boundary = KokunoAppendixBBoundary(max_step=0.08)
    path = boundary.save_json(tmp_path / "appendix_b_boundary.json")
    replay = KokunoAppendixBBoundary.load_json(path)
    assert replay.to_payload() == boundary.to_payload()
    assert replay.sha256 == boundary.sha256

    payload = copy.deepcopy(boundary.to_payload())
    payload["truth_boundary"]["paper_exact"] = True
    payload["sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoAppendixBBoundary.from_payload(payload)


def test_invalid_final_widths_fail_closed():
    with pytest.raises(ValueError, match="fit strictly"):
        KokunoAppendixBBoundary(
            axial_shutdown_log_width=0.06,
            angular_settle_log_width=0.06,
        )
