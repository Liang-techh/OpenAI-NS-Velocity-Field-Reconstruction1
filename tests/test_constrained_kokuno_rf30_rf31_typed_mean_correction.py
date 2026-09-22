from __future__ import annotations

from dataclasses import replace
import inspect
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_rf30_rf31_typed_mean_correction import (
    CandidateIdentity,
    RF30PreMeanState,
    RF31CorrectionBasis,
    RF30RF31ContractError,
    materialize_rf30_rf31_mean_correction_system,
    truth_boundary,
)


def _trapz(y: np.ndarray, x: np.ndarray) -> float:
    total = 0.0
    for i in range(x.size - 1):
        total += 0.5 * (float(y[i]) + float(y[i + 1])) * float(x[i + 1] - x[i])
    return total


class MechanicsBackend:
    def __init__(self, *, repository: bool = False) -> None:
        self.R = np.array([1.0, 1.2, 1.5, 1.8, 2.0])
        kind = "repository-candidate" if repository else "mechanics-only"
        self.identity = CandidateIdentity(
            candidate_id="fixture-current-i4",
            candidate_sha256="a" * 64,
            evidence_kind=kind,
        )

        R = self.R
        self.state = RF30PreMeanState(
            source_candidate_id=self.identity.candidate_id,
            source_candidate_sha256=self.identity.candidate_sha256,
            source_chart_id="fixed-q-chart-fixture",
            source_chart_sha256="b" * 64,
            radius_R=tuple(R),
            base_V=tuple(0.8 + 0.15 * R),
            base_G=tuple(-0.25 + 0.10 * R),
            mean_W_rr=tuple(0.04 * R**2),
            mean_W_zr=tuple(0.03 * R),
            mean_W_thetatheta=tuple(0.02 + 0.01 * R),
            mean_W_ztheta=tuple(0.015 * R**2),
            mean_W_zz=tuple(0.05 - 0.01 * R),
            dR_mean_W_rr=tuple(0.08 * R),
            dZ_mean_W_zr=tuple(0.01 + 0.002 * R),
            actual_candidate_recomputed=repository,
            oscillatory_covariance_recomputed=repository,
            normalized_haar_mean_used=repository,
            source_fixed_q_chart_used=repository,
        )

        delta_v = []
        gamma_d = []
        for j in range(5):
            delta_v.append(tuple((R - 0.9) ** (j + 1)))
            gamma_d.append(tuple((0.2 + 0.05 * j) * (2.1 - R) ** (5 - j)))
        self.basis = RF31CorrectionBasis(
            source_candidate_id=self.identity.candidate_id,
            source_candidate_sha256=self.identity.candidate_sha256,
            source_chart_id=self.state.source_chart_id,
            source_chart_sha256=self.state.source_chart_sha256,
            radius_R=tuple(R),
            delta_v_columns=tuple(delta_v),
            gamma_d_columns=tuple(gamma_d),
            basis_selected_before_defect_evaluation=True,
            compact_support_verified=True,
            torus_independent_verified=True,
        )

    def candidate_identity(self, candidate: object) -> CandidateIdentity:
        assert candidate is self
        return self.identity

    def rf30_pre_mean_state(self, candidate: object) -> RF30PreMeanState:
        assert candidate is self
        return self.state

    def rf31_correction_basis(
        self, candidate: object, state: RF30PreMeanState
    ) -> RF31CorrectionBasis:
        assert candidate is self
        assert state is self.state
        return self.basis


def test_public_surface_exposes_no_manual_defect_or_tuning_knob() -> None:
    signature = inspect.signature(materialize_rf30_rf31_mean_correction_system)
    assert tuple(signature.parameters) == ("backend", "candidate")
    forbidden = {
        "P",
        "J_theta",
        "J_z",
        "residual",
        "defect",
        "forcing",
        "pressure",
        "correction",
        "gain",
        "damping",
        "threshold",
        "heldout",
    }
    assert forbidden.isdisjoint(signature.parameters)


def test_rf30_initial_defects_match_independent_formula_replay() -> None:
    backend = MechanicsBackend()
    receipt = materialize_rf30_rf31_mean_correction_system(backend, backend)

    R = backend.R
    s = backend.state
    Wrr = np.asarray(s.mean_W_rr)
    Wtt = np.asarray(s.mean_W_thetatheta)
    dR_Wrr = np.asarray(s.dR_mean_W_rr)
    dZ_Wzr = np.asarray(s.dZ_mean_W_zr)
    Wzt = np.asarray(s.mean_W_ztheta)
    Wzz = np.asarray(s.mean_W_zz)

    g_r = -(dR_Wrr + Wrr / R) - dZ_Wzr + Wtt / R
    expected_P = _trapz(g_r, R)
    expected_Jtheta = _trapz(R**2 * Wzt, R)
    expected_Jz = _trapz(R * Wzz, R) - 0.5 * _trapz(R**2 * g_r, R)

    assert receipt.defect.P == pytest.approx(expected_P, abs=1e-15)
    assert receipt.defect.J_theta == pytest.approx(expected_Jtheta, abs=1e-15)
    assert receipt.defect.J_z == pytest.approx(expected_Jz, abs=1e-15)
    assert receipt.defect.target_rows == pytest.approx(
        (0.0, 0.0, -expected_P, -expected_Jtheta, -expected_Jz),
        abs=1e-15,
    )
    assert np.asarray(receipt.defect.radial_source_g_r) == pytest.approx(g_r, abs=1e-15)
    assert receipt.repository_candidate_defect_evidence is False
    assert receipt.linear_five_row_system_materialized is True


def test_rf31_five_rows_match_independent_column_replay() -> None:
    backend = MechanicsBackend()
    receipt = materialize_rf30_rf31_mean_correction_system(backend, backend)
    R = backend.R
    V = np.asarray(backend.state.base_V)
    G = np.asarray(backend.state.base_G)
    dv = np.asarray(backend.basis.delta_v_columns)
    gd = np.asarray(backend.basis.gamma_d_columns)

    expected = np.zeros((5, 5))
    for j in range(5):
        expected[0, j] = _trapz(R**2 * dv[j], R)
        expected[1, j] = _trapz(R * gd[j], R)
        expected[2, j] = _trapz((2.0 * V / R) * dv[j], R)
        expected[3, j] = _trapz(R**2 * (G * dv[j] + V * gd[j]), R)
        expected[4, j] = _trapz(2.0 * R * G * gd[j] - R * V * dv[j], R)

    actual = np.asarray(receipt.matrix_B)
    assert actual == pytest.approx(expected, abs=1e-14)
    assert receipt.matrix_rank == int(np.linalg.matrix_rank(expected))
    assert math.isfinite(receipt.matrix_condition_2) or math.isinf(
        receipt.matrix_condition_2
    )


def test_repository_candidate_requires_actual_source_chart_recomputation() -> None:
    backend = MechanicsBackend(repository=True)
    receipt = materialize_rf30_rf31_mean_correction_system(backend, backend)
    assert receipt.repository_candidate_defect_evidence is True

    backend.state = replace(backend.state, actual_candidate_recomputed=False)
    with pytest.raises(
        RF30RF31ContractError,
        match="was not recomputed from candidate",
    ):
        materialize_rf30_rf31_mean_correction_system(backend, backend)

    backend = MechanicsBackend(repository=True)
    backend.state = replace(backend.state, source_fixed_q_chart_used=False)
    with pytest.raises(RF30RF31ContractError, match="not in source fixed-Q chart"):
        materialize_rf30_rf31_mean_correction_system(backend, backend)


@pytest.mark.parametrize(
    ("state_change", "basis_change", "match"),
    [
        ({"surrogate_defect_used": True}, {}, "surrogate RF30 defect"),
        ({"residual_as_forcing_shortcut_used": True}, {}, "residual-as-forcing"),
        ({"heldout_samples_used_to_construct_state": True}, {}, "held-out samples"),
        ({}, {"surrogate_target_used_to_choose_basis": True}, "surrogate defect"),
        ({}, {"heldout_samples_used_to_choose_basis": True}, "held-out samples"),
        ({}, {"basis_selected_before_defect_evaluation": False}, "frozen before defect"),
        ({}, {"compact_support_verified": False}, "compact support"),
        ({}, {"torus_independent_verified": False}, "torus independent"),
    ],
)
def test_surrogate_leakage_and_unregistered_basis_fail_closed(
    state_change: dict[str, object],
    basis_change: dict[str, object],
    match: str,
) -> None:
    backend = MechanicsBackend()
    backend.state = replace(backend.state, **state_change)
    backend.basis = replace(backend.basis, **basis_change)
    with pytest.raises(RF30RF31ContractError, match=match):
        materialize_rf30_rf31_mean_correction_system(backend, backend)


def test_candidate_and_chart_identity_mismatch_fail_closed() -> None:
    backend = MechanicsBackend()
    backend.state = replace(backend.state, source_candidate_sha256="c" * 64)
    with pytest.raises(RF30RF31ContractError, match="candidate hash mismatch"):
        materialize_rf30_rf31_mean_correction_system(backend, backend)

    backend = MechanicsBackend()
    backend.basis = replace(backend.basis, source_chart_id="other-chart")
    with pytest.raises(RF30RF31ContractError, match="source chart id mismatch"):
        materialize_rf30_rf31_mean_correction_system(backend, backend)


def test_truth_boundary_does_not_promote_current_i4_or_ns_cycle() -> None:
    boundary = truth_boundary()
    assert boundary["rf30_rf31_formula_executor_materialized"] is True
    assert boundary["current_i4_source_chart_backend_materialized"] is False
    assert boundary["current_i4_rf30_defect_materialized"] is False
    assert boundary["current_i4_rf31_five_row_system_materialized"] is False
    assert boundary["correction_coefficients_solved"] is False
    assert boundary["correction_applied"] is False
    assert boundary["rf44_rf49_nonlinear_remainder_recomputed"] is False
    assert boundary["complete_ns_defect"] is False
    assert boundary["finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["pde_validated"] is False
    assert boundary["final_normalized_momentum_gate"] == 1e-3
    assert boundary["final_normalized_divergence_gate"] == 1e-5
    assert boundary["residual_defined_free_forcing_allowed"] is False
    assert boundary["finite_stage_small_residual_is_blowup_proof"] is False
