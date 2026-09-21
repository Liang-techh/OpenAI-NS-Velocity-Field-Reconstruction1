from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_fixed_kappa_x100 import (
    CONTINUATION_QUADRATURE_ORDER,
    KokunoPA10FixedKappaContinuationToX100,
)


def _field() -> KokunoPA10FixedKappaContinuationToX100:
    return KokunoPA10FixedKappaContinuationToX100()


def test_activation_endpoint_value_handoff_is_exact() -> None:
    field = _field()
    eta = np.asarray([-0.6, -0.2, 0.0, 0.25, 0.6])
    X = np.full_like(eta, field.X_1)

    continued = field.values(X, eta)
    endpoint = field.activation.values(X, eta)

    np.testing.assert_allclose(
        continued["F_fixed_kappa"], endpoint["F_activation"], rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        continued["U_fixed_kappa"], endpoint["U_activation"], rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        continued["delta_log_F_from_activation_endpoint"], 0.0, rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        continued["delta_U_from_activation_endpoint"], 0.0, rtol=0.0, atol=0.0
    )


def test_fixed_kappa_public_odes_replay_with_independent_centered_difference() -> None:
    field = _field()
    span = field.X_b_ref - field.X_1
    X = field.X_1 + span * np.asarray([0.22, 0.57, 0.83])
    eta = np.asarray([-0.35, 0.1, 0.45])
    h = 2.0e-4 * np.maximum(1.0, np.abs(X))

    plus = field.values(X + h, eta)
    minus = field.values(X - h, eta)
    F_X_fd = (plus["F_fixed_kappa"] - minus["F_fixed_kappa"]) / (2.0 * h)
    U_X_fd = (plus["U_fixed_kappa"] - minus["U_fixed_kappa"]) / (2.0 * h)
    analytic = field.radial_derivatives(X, eta)

    np.testing.assert_allclose(
        F_X_fd, analytic["F_fixed_kappa_X"], rtol=6.0e-3, atol=2.0e-8
    )
    np.testing.assert_allclose(
        U_X_fd, analytic["U_fixed_kappa_X"], rtol=6.0e-3, atol=2.0e-8
    )


def test_x100_handoff_is_vectorized_nontrivial_and_shear_is_composed() -> None:
    field = _field()
    eta = np.asarray([-0.55, -0.15, 0.0, 0.3, 0.55])
    handoff = field.handoff_at_X100(eta)

    assert handoff["F_fixed_kappa"].shape == eta.shape
    assert handoff["U_fixed_kappa"].shape == eta.shape
    assert handoff["E_fixed_kappa"].shape == eta.shape
    assert np.all(np.isfinite(handoff["F_fixed_kappa"]))
    assert np.all(np.isfinite(handoff["U_fixed_kappa"]))
    assert np.all(np.isfinite(handoff["E_fixed_kappa"]))
    assert np.max(np.abs(handoff["U_fixed_kappa"])) > 1.0e-12
    assert np.max(np.abs(handoff["shear_first"])) > 1.0e-12

    np.testing.assert_allclose(
        handoff["shear_first"],
        field.kappa_0 * handoff["p1_reference"],
        rtol=0.0,
        atol=0.0,
    )
    expected_second = (
        field.kappa_0
        * handoff["p2_reference_autonomous"]
        * handoff["E_reference"]
        / handoff["E_fixed_kappa"]
    )
    np.testing.assert_allclose(
        handoff["shear_second"], expected_second, rtol=2.0e-14, atol=2.0e-14
    )


def test_quadrature_refinement_is_engineering_replay_not_pde_evidence() -> None:
    field = _field()
    eta = np.asarray([-0.4, 0.2, 0.5])
    X = np.full_like(eta, field.X_b_ref)
    coarse = field.values(X, eta, order=CONTINUATION_QUADRATURE_ORDER)
    fine = field.values(X, eta, order=64)

    np.testing.assert_allclose(
        coarse["F_fixed_kappa"], fine["F_fixed_kappa"], rtol=3.0e-3, atol=2.0e-9
    )
    np.testing.assert_allclose(
        coarse["U_fixed_kappa"], fine["U_fixed_kappa"], rtol=3.0e-3, atol=2.0e-9
    )


def test_configuration_roundtrip_and_truth_boundary(tmp_path) -> None:
    field = _field()
    path = tmp_path / "fixed_kappa_x100.json"
    saved = field.save_configuration(path)
    rebound = KokunoPA10FixedKappaContinuationToX100.load_configuration(path)

    assert rebound.configuration() == saved
    assert rebound.semantic_sha256 == field.semantic_sha256
    assert field.kappa_0 == 0.1

    truth = field.truth_boundary
    assert truth["public_fixed_kappa_continuation_formula_executable"] is True
    assert truth["candidate_side_fixed_kappa_continuation_to_X100_executable"] is True
    assert truth["selected_kappa0_chosen_from_ns_residual"] is False
    assert truth["selected_kappa0_global_source_smallness_admitted"] is False
    assert truth["source_prepared_reference_nsr_materialized"] is False
    assert truth["source_exact_fixed_kappa_continuation_to_X100_materialized"] is False
    assert truth["final_interpolation_to_Xi_materialized"] is False
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["same_protocol_comparable_to_st006"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    mutated = copy.deepcopy(saved)
    mutated["selected_kappa0"] = 0.1001
    with pytest.raises(ValueError, match="configuration is frozen"):
        KokunoPA10FixedKappaContinuationToX100.from_configuration(mutated)


def test_domain_is_fail_closed() -> None:
    field = _field()
    with pytest.raises(ValueError, match="fixed-kappa continuation domain"):
        field.values(field.X_1 * 0.99, 0.0)
    with pytest.raises(ValueError, match="fixed-kappa continuation domain"):
        field.values(field.X_b_ref + 0.1, 0.0)
    eta_lo, eta_hi = field.eta_interval
    with pytest.raises(ValueError, match="activation interval"):
        field.values(field.X_b_ref, eta_hi + 1.0e-4)
