import copy

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_outer_reserved_patch_schedule import (
    KokunoOuterReservedPatchSchedule,
)
from openai_ns_reconstruction.kokuno_pa10_actual_xi_prefix_moments import (
    KokunoPA10ActualXiPrefixMoments,
)


@pytest.fixture(scope="module")
def moments():
    return KokunoPA10ActualXiPrefixMoments(quadrature_order=24)


def test_unified_profile_replays_all_four_current_lineage_segments(moments):
    eta = 0.2
    X = np.asarray(
        [
            0.5 * moments.X_0,
            0.5 * (moments.X_0 + moments.X_1),
            0.5 * (moments.X_1 + moments.X_b_ref),
            105.0,
        ]
    )
    unified = moments.profile_values(X, eta)

    natural = moments.bridge.actual_x100.activation.reference.source_normalized.values(
        X[0], eta
    )
    activation = moments.bridge.actual_x100.activation.values(X[1], eta)
    fixed = moments.bridge.actual_x100.values(X[2], eta)
    final = moments.bridge.values(X[3], eta)

    expected_F = np.asarray(
        [
            natural["F_0"],
            activation["F_activation"],
            fixed["F_fixed_kappa"],
            final["F_final_bridge"],
        ],
        dtype=float,
    )
    expected_U = np.asarray(
        [
            natural["U_0"],
            activation["U_activation"],
            fixed["U_fixed_kappa"],
            final["U_final_bridge"],
        ],
        dtype=float,
    )
    np.testing.assert_allclose(unified["F_actual_prefix"], expected_F, rtol=0, atol=0)
    np.testing.assert_allclose(unified["U_actual_prefix"], expected_U, rtol=0, atol=0)
    np.testing.assert_allclose(
        unified["E_actual_prefix"],
        np.sqrt(2.0 * X) * expected_F,
        rtol=2e-15,
        atol=0,
    )


def test_unified_radial_derivatives_replay_public_value_evaluator(moments):
    eta = 0.15
    for X, step in ((0.5 * moments.X_0, 2e-6), (20.0, 2e-4), (105.0, 2e-4)):
        h = step * max(1.0, X)
        plus = moments.profile_values(X + h, eta)
        minus = moments.profile_values(X - h, eta)
        analytic = moments.radial_derivatives(X, eta)
        for value_key, derivative_key in (
            ("F_actual_prefix", "F_actual_prefix_X"),
            ("U_actual_prefix", "U_actual_prefix_X"),
            ("E_actual_prefix", "E_actual_prefix_X"),
        ):
            fd = (plus[value_key] - minus[value_key]) / (2.0 * h)
            scale = max(1e-9, abs(float(fd)), abs(float(analytic[derivative_key])))
            assert abs(float(fd - analytic[derivative_key])) / scale < 3e-3


def test_public_moment_density_identities_and_ideal_Cp_cancellation(moments):
    eta = -0.25
    X = np.asarray([0.2 * moments.X_0, 8.0, 104.0])
    values = moments.profile_values(X, eta)
    d = moments.moment_densities(X, eta)
    F = values["F_actual_prefix"]
    U = values["U_actual_prefix"]
    H = 2.0 * X * F
    np.testing.assert_allclose(d["M_density"], U, rtol=0, atol=0)
    np.testing.assert_allclose(d["I_density"], H, rtol=0, atol=0)
    np.testing.assert_allclose(d["J_density"], U * H, rtol=0, atol=0)
    np.testing.assert_allclose(d["S_density"], U * U - X * F * F, rtol=0, atol=0)
    np.testing.assert_allclose(d["C_p_density"], F * F, rtol=0, atol=0)

    ideal = moments.ideal_PA15_scaled_moments(eta)
    f = 1.0 / (1.0 + eta * eta)
    assert ideal[4] == pytest.approx(2.5 * f * f / moments.C**2, rel=2e-15)


def test_current_candidate_side_prefix_discrepancy_is_finite_nontrivial_and_refines(moments):
    eta = 0.2
    physical_24 = moments.physical_prefix_moments_at_Xi(eta, order=24)
    physical_32 = moments.physical_prefix_moments_at_Xi(eta, order=32)
    assert np.all(np.isfinite(physical_24))
    assert np.all(np.isfinite(physical_32))
    np.testing.assert_allclose(physical_24, physical_32, rtol=8e-4, atol=3e-8)

    actual = moments.actual_PA15_scaled_moments(eta, order=32)
    ideal = moments.ideal_PA15_scaled_moments(eta)
    discrepancy = moments.incoming_PA15_discrepancy(np.asarray(eta), order=32)
    np.testing.assert_allclose(discrepancy, actual - ideal, rtol=0, atol=0)
    assert np.all(np.isfinite(discrepancy))
    assert float(np.max(np.abs(discrepancy))) > 1e-10


def test_vectorized_discrepancy_and_pa16_input_bind_current_Xi_handoff(moments):
    eta = np.asarray([-0.2, 0.0, 0.2])
    discrepancy = moments.incoming_PA15_discrepancy(eta, order=20)
    assert discrepancy.shape == (3, 5)
    assert np.all(np.isfinite(discrepancy))

    receipt = moments.pa16_input_at_eta(0.2, order=20)
    handoff = moments.bridge.handoff_at_Xi(np.asarray(0.2))
    assert receipt["ell_i"] == pytest.approx(float(handoff["ell_i"]), rel=0, abs=0)
    assert receipt["G_i"] == pytest.approx(float(handoff["G_i"]), rel=0, abs=0)
    np.testing.assert_allclose(
        receipt["incoming_scaled_discrepancy"], discrepancy[2], rtol=0, atol=0
    )
    assert receipt["outer_schedule_sha256"] == moments.outer_schedule.sha256
    assert receipt["source_T_sh_lower_bound_verified"] is False


def test_configuration_roundtrip_and_outer_scale_tamper_fail_closed(moments, tmp_path):
    path = tmp_path / "actual_xi_prefix_moments.json"
    moments.save_configuration(path)
    replay = KokunoPA10ActualXiPrefixMoments.load_configuration(path)
    assert replay.configuration() == moments.configuration()
    assert replay.semantic_sha256 == moments.semantic_sha256

    payload = copy.deepcopy(moments.configuration())
    payload["outer_schedule"]["parameters"]["C"] = 999.0
    with pytest.raises(ValueError):
        KokunoPA10ActualXiPrefixMoments.from_configuration(payload)

    with pytest.raises(ValueError, match="outer schedule C"):
        KokunoPA10ActualXiPrefixMoments(
            outer_schedule=KokunoOuterReservedPatchSchedule(C=999.0),
            quadrature_order=20,
        )


def test_truth_boundary_remains_fail_closed_beyond_candidate_side_Xi_moments(moments):
    truth = moments.truth_boundary
    assert truth["candidate_side_upstream_five_moment_discrepancy_at_Xi_materialized"] is True
    assert truth["candidate_side_PA16_input_tuple_materialized"] is True
    assert truth["source_prepared_upstream_five_moment_discrepancy_materialized"] is False
    assert truth["source_T_sh_lower_bound_verified"] is False
    assert truth["five_moment_repair_applied"] is False
    assert truth["inner_to_outer_join_completed"] is False
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_domain_guards(moments):
    with pytest.raises(ValueError):
        moments.profile_values(-1e-6, 0.0)
    with pytest.raises(ValueError):
        moments.profile_values(moments.X_i + 1e-6, 0.0)
    with pytest.raises(ValueError, match="X>0"):
        moments.radial_derivatives(0.0, 0.0)
