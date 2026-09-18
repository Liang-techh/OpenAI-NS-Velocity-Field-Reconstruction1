import copy
import hashlib
import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_rescaled_appendix_b_incoming_moments import (
    KokunoRescaledAppendixBIncomingMoments,
)


@pytest.fixture(scope="module")
def selected_receipt():
    moments = KokunoRescaledAppendixBIncomingMoments(axis_quadrature_points=48)
    eta = float(moments.boundary.reference.seed.phase_stationary_eta)
    return moments, eta, moments.result_at_eta(eta)


def test_rescaled_pa15_target_is_finite_nontrivial_and_matches_Xi_endpoint(selected_receipt):
    moments, eta, receipt = selected_receipt
    boundary = moments.boundary.boundary_values(np.asarray(eta))

    assert receipt.ell_i == pytest.approx(float(boundary["ell_i"]), rel=2e-8, abs=2e-6)
    # The 10-state moment augmentation and the authoritative 6-state boundary
    # solve use the same ODE/tolerances but adaptive DOP853 controls a different
    # error norm.  Their endpoint agreement is therefore a numerical-consistency
    # check at the solver scale, not a bitwise identity requirement.
    assert receipt.G_i == pytest.approx(float(boundary["G_i"]), rel=5e-7, abs=5e-7)
    assert np.all(np.isfinite(receipt.actual_physical_moments))
    assert np.all(np.isfinite(receipt.actual_scaled_moments))
    assert np.all(np.isfinite(receipt.ideal_scaled_moments))
    assert np.all(np.isfinite(receipt.incoming_scaled_discrepancy))
    assert receipt.actual_physical_moments[4] > 0.0
    assert receipt.max_abs_incoming_discrepancy > 1e-8

    truth = moments.to_payload()["truth_boundary"]
    assert truth["selected_source_rescaled_appendix_B_prefix_moments_executable"] is True
    assert truth["shared_C_PA15_scaling_applied"] is True
    assert truth["rescaled_PA15_incoming_discrepancy_executable"] is True
    assert truth["source_T_sh_lower_bound_verified"] is False
    assert truth["inner_to_outer_join_completed"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_shared_C_scale_is_not_old_finite_template_and_underflow_is_explicit(selected_receipt):
    moments, eta, receipt = selected_receipt

    assert moments.binding.finite_template_C_matches_selected_C is False
    assert moments.binding.delta_log_X_R_from_finite_template > 1e10
    assert moments.log_x_i < -1e10

    # On this selected shared-C scale, normalized M/I/J/S are smaller than the
    # binary64 materialization range.  That fact is recorded, not relabeled as
    # an exact source zero.  C_p is scale invariant and remains nontrivial.
    assert all(receipt.actual_scaled_underflow_rows[:4])
    assert receipt.actual_scaled_underflow_rows[4] is False
    assert all(receipt.ideal_scaled_underflow_rows)
    np.testing.assert_array_equal(np.asarray(receipt.actual_scaled_moments[:4]), 0.0)
    np.testing.assert_array_equal(np.asarray(receipt.ideal_scaled_moments), 0.0)
    assert receipt.incoming_scaled_discrepancy[4] == pytest.approx(
        receipt.actual_scaled_moments[4], rel=0.0, abs=0.0
    )
    assert receipt.incoming_scaled_discrepancy[4] > 0.0

    handoff = moments.pa16_inputs_at_eta(eta)
    assert handoff["row_order"] == ["M", "I", "J", "S", "C_p"]
    assert handoff["source_T_sh_lower_bound_verified"] is False
    np.testing.assert_array_equal(
        np.asarray(handoff["incoming_scaled_discrepancy"]),
        np.asarray(receipt.incoming_scaled_discrepancy),
    )


def test_vectorized_selected_discrepancy_has_pa15_row_shape(selected_receipt):
    moments, eta, receipt = selected_receipt
    vectorized = moments.incoming_scaled_discrepancy(np.asarray([eta]))
    assert vectorized.shape == (1, 5)
    np.testing.assert_allclose(
        vectorized[0],
        np.asarray(receipt.incoming_scaled_discrepancy),
        rtol=3e-8,
        atol=2e-8,
    )


def test_axis_quadrature_refinement_stabilizes_scale_invariant_Cp():
    coarse = KokunoRescaledAppendixBIncomingMoments(axis_quadrature_points=32)
    fine = KokunoRescaledAppendixBIncomingMoments(axis_quadrature_points=64)
    eta = float(coarse.boundary.reference.seed.phase_stationary_eta)
    a = coarse.result_at_eta(eta)
    b = fine.result_at_eta(eta)

    assert a.incoming_scaled_discrepancy[4] == pytest.approx(
        b.incoming_scaled_discrepancy[4], rel=2e-6, abs=2e-8
    )
    assert a.G_i == pytest.approx(b.G_i, rel=5e-7, abs=5e-7)
    assert a.ell_i == pytest.approx(b.ell_i, rel=2e-12, abs=2e-4)


def test_payload_roundtrip_and_truth_tamper_fail_closed(tmp_path):
    moments = KokunoRescaledAppendixBIncomingMoments(axis_quadrature_points=48)
    path = moments.save_json(tmp_path / "rescaled_appendix_b_incoming_moments.json")
    replay = KokunoRescaledAppendixBIncomingMoments.load_json(path)
    assert replay.to_payload() == moments.to_payload()
    assert replay.sha256 == moments.sha256

    payload = copy.deepcopy(moments.to_payload())
    payload["truth_boundary"]["paper_exact"] = True
    unsigned = {key: value for key, value in payload.items() if key != "sha256"}
    payload["sha256"] = hashlib.sha256(
        json.dumps(unsigned, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoRescaledAppendixBIncomingMoments.from_payload(payload)
