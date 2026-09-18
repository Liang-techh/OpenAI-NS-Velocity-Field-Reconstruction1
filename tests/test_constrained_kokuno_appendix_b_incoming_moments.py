import copy
import hashlib
import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_appendix_b_incoming_moments import (
    KokunoAppendixBIncomingMoments,
)
from openai_ns_reconstruction.kokuno_inner_join_exit import KokunoInnerJoinExit
from openai_ns_reconstruction.kokuno_outer_reserved_patch_schedule import (
    KokunoOuterReservedPatchSchedule,
)


def test_selected_appendix_b_prefix_moments_are_finite_nontrivial_and_match_Xi_boundary():
    moments = KokunoAppendixBIncomingMoments(axis_quadrature_points=64)
    receipt = moments.result_at_eta(0.2)
    boundary = moments.boundary._solve_eta(0.2)

    assert receipt.ell_i == pytest.approx(boundary.ell_i, rel=3e-7, abs=3e-9)
    assert receipt.G_i == pytest.approx(boundary.G_i, rel=3e-7, abs=3e-9)
    assert np.all(np.isfinite(receipt.actual_scaled_moments))
    assert np.all(np.isfinite(receipt.ideal_scaled_moments))
    assert np.all(np.isfinite(receipt.incoming_scaled_discrepancy))
    assert receipt.max_abs_incoming_discrepancy > 1e-10

    truth = moments.to_payload()["truth_boundary"]
    assert truth["selected_appendix_B_prefix_moments_executable"] is True
    assert truth["source_PA15_scaled_incoming_discrepancy_executable"] is True
    assert truth["source_T_sh_lower_bound_verified"] is False
    assert truth["inner_to_outer_join_completed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_axis_quadrature_refinement_stabilizes_incoming_target():
    coarse = KokunoAppendixBIncomingMoments(axis_quadrature_points=64)
    fine = KokunoAppendixBIncomingMoments(axis_quadrature_points=128)
    a = np.asarray(coarse.result_at_eta(-0.2).incoming_scaled_discrepancy)
    b = np.asarray(fine.result_at_eta(-0.2).incoming_scaled_discrepancy)
    np.testing.assert_allclose(a, b, rtol=3e-6, atol=3e-10)


def test_ideal_prefix_is_analytic_and_vectorized_discrepancy_has_PA15_row_shape():
    moments = KokunoAppendixBIncomingMoments(axis_quadrature_points=64)
    eta = 0.0
    ideal = moments.ideal_scaled_moments(eta)
    x_i = moments.x_i
    P = moments.P_star

    assert ideal[0] == 0.0
    assert ideal[2] == 0.0
    assert ideal[4] == pytest.approx(2.5 * P * P * x_i**0.2, rel=2e-14)

    grid = np.asarray([-0.2, 0.0, 0.2])
    discrepancy = moments.incoming_scaled_discrepancy(grid)
    assert discrepancy.shape == (3, 5)
    assert np.all(np.isfinite(discrepancy))
    scalar = np.asarray(moments.result_at_eta(0.2).incoming_scaled_discrepancy)
    np.testing.assert_allclose(discrepancy[2], scalar, rtol=0.0, atol=0.0)


def test_real_incoming_target_routes_additively_into_existing_PA16_pre_repair_map():
    moments = KokunoAppendixBIncomingMoments(axis_quadrature_points=64)
    receipt = moments.result_at_eta(0.2)
    join = KokunoInnerJoinExit(T_sh=1.0, outer_schedule=moments.outer_schedule)
    incoming = np.asarray(receipt.incoming_scaled_discrepancy)
    pre_real = join.pre_repair_scaled_discrepancy(
        eta=receipt.eta,
        ell_i=receipt.ell_i,
        G_i=receipt.G_i,
        incoming_scaled_discrepancy=incoming,
    )
    pre_zero = join.pre_repair_scaled_discrepancy(
        eta=receipt.eta,
        ell_i=receipt.ell_i,
        G_i=receipt.G_i,
        incoming_scaled_discrepancy=np.zeros(5),
    )
    np.testing.assert_allclose(pre_real - pre_zero, incoming, rtol=3e-8, atol=3e-10)


def test_payload_roundtrip_and_truth_tamper_fail_closed(tmp_path):
    moments = KokunoAppendixBIncomingMoments(axis_quadrature_points=64)
    path = moments.save_json(tmp_path / "appendix_b_incoming_moments.json")
    replay = KokunoAppendixBIncomingMoments.load_json(path)
    assert replay.to_payload() == moments.to_payload()
    assert replay.sha256 == moments.sha256

    payload = copy.deepcopy(moments.to_payload())
    payload["truth_boundary"]["paper_exact"] = True
    unsigned = {key: value for key, value in payload.items() if key != "sha256"}
    payload["sha256"] = hashlib.sha256(
        json.dumps(unsigned, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoAppendixBIncomingMoments.from_payload(payload)


def test_mismatched_outer_scale_binding_fails_closed():
    schedule = KokunoOuterReservedPatchSchedule(C=3.0)
    with pytest.raises(ValueError, match="boundary C"):
        KokunoAppendixBIncomingMoments(outer_schedule=schedule)
