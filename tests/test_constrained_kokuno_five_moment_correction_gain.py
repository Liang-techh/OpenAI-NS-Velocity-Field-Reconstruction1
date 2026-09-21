from __future__ import annotations

from dataclasses import replace
import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_five_moment_correction_gain import (
    FIXED_POINT_MAX_ITERATIONS,
    MOMENT_DIMENSION,
    SOURCE_LIPSCHITZ_LIMIT,
    SOURCE_READER_HEAD,
    SOURCE_READER_TEX_BLOB,
    SOURCE_SMALLNESS_LIMIT,
    CandidateIdentity,
    FiveMomentCorrectionSystem,
    MomentCorrectionGainError,
    assess_five_moment_correction_gain,
    build_mechanics_report,
    solve_five_moment_correction,
    truth_boundary,
)


class Candidate:
    candidate_id = "candidate-a"
    candidate_sha256 = "a" * 64


class Backend:
    def __init__(
        self,
        *,
        quadratic_scale: float = 0.10,
        discrepancy_scale: float = 1.0,
        candidate_hash: str | None = None,
        residual_as_forcing: bool = False,
        matrix: np.ndarray | None = None,
    ) -> None:
        self.quadratic_scale = float(quadratic_scale)
        self.discrepancy_scale = float(discrepancy_scale)
        self.candidate_hash = candidate_hash
        self.residual_as_forcing = bool(residual_as_forcing)
        self.matrix = np.eye(MOMENT_DIMENSION) if matrix is None else np.asarray(matrix, dtype=float)

    def candidate_identity(self, candidate: Candidate) -> CandidateIdentity:
        return CandidateIdentity(
            candidate_id=candidate.candidate_id,
            candidate_sha256=candidate.candidate_sha256,
            evidence_kind="mechanics-only",
        )

    def five_moment_correction_system(self, candidate: Candidate) -> FiveMomentCorrectionSystem:
        Q = np.zeros((MOMENT_DIMENSION, MOMENT_DIMENSION, MOMENT_DIMENSION))
        for i in range(MOMENT_DIMENSION):
            Q[i, i, i] = self.quadratic_scale
        d = self.discrepancy_scale * np.array([0.05, -0.04, 0.03, -0.02, 0.01])
        return FiveMomentCorrectionSystem(
            source_candidate_id=candidate.candidate_id,
            source_candidate_sha256=(
                candidate.candidate_sha256 if self.candidate_hash is None else self.candidate_hash
            ),
            fit_partition_id="held-in-moments",
            fit_sample_ids=("s0", "s1", "s2", "s3", "s4"),
            matrix_B=tuple(tuple(float(v) for v in row) for row in self.matrix),
            bilinear_Q=tuple(
                tuple(tuple(float(v) for v in row) for row in plane) for plane in Q
            ),
            discrepancy_d=tuple(float(v) for v in d),
            discrepancy_kind="mechanics-only",
            discrepancy_from_complete_ns_defect=False,
            residual_as_forcing_shortcut_used=self.residual_as_forcing,
        )


def test_source_bound_is_recomputed_and_positive_case_is_certified() -> None:
    diagnostic = assess_five_moment_correction_gain(Backend(), Candidate())
    assert diagnostic.inverse_norm_nu == pytest.approx(1.0)
    assert diagnostic.matrix_condition_2 == pytest.approx(1.0)
    assert diagnostic.bilinear_norm_upper_q == pytest.approx(np.sqrt(5.0) * 0.10)
    assert diagnostic.registered_radius == pytest.approx(
        2.0 * diagnostic.inverse_norm_nu * diagnostic.discrepancy_norm
    )
    assert diagnostic.source_smallness_parameter == pytest.approx(
        8.0
        * diagnostic.inverse_norm_nu**2
        * diagnostic.bilinear_norm_upper_q
        * diagnostic.discrepancy_norm
    )
    assert diagnostic.image_norm_upper <= diagnostic.registered_radius
    assert diagnostic.lipschitz_upper <= SOURCE_LIPSCHITZ_LIMIT
    assert diagnostic.source_smallness_parameter <= SOURCE_SMALLNESS_LIMIT
    assert diagnostic.source_smallness_condition_met
    assert diagnostic.self_map_certified
    assert diagnostic.contraction_certified
    assert not diagnostic.rejection_reasons


def test_bounded_fixed_point_solves_the_exact_quadratic_moment_equation() -> None:
    backend = Backend()
    candidate = Candidate()
    receipt = solve_five_moment_correction(backend, candidate)
    assert 1 <= receipt.iterations <= FIXED_POINT_MAX_ITERATIONS
    assert receipt.local_moment_fixed_point_converged
    assert receipt.local_moment_correction_accepted
    assert receipt.correction_nontrivial
    assert receipt.equation_residual_norm <= 1.0e-12
    assert receipt.coefficient_norm <= receipt.diagnostic.registered_radius * (1.0 + 1.0e-12)
    assert not receipt.complete_ns_correction_authorized
    assert not receipt.analytic_full_ns_cycle_gain_available
    assert not receipt.heldout_ns_residual_assessed
    assert not receipt.pde_validated

    system = backend.five_moment_correction_system(candidate)
    B = np.asarray(system.matrix_B)
    Q = np.asarray(system.bilinear_Q)
    d = np.asarray(system.discrepancy_d)
    c = np.asarray(receipt.coefficient_vector)
    residual = B @ c + np.einsum("ijk,j,k->i", Q, c, c) - d
    assert np.linalg.norm(residual) <= 1.0e-12


def test_large_quadratic_map_is_rejected_before_iteration() -> None:
    backend = Backend(quadratic_scale=20.0)
    diagnostic = assess_five_moment_correction_gain(backend, Candidate())
    assert not diagnostic.source_smallness_condition_met
    assert not diagnostic.contraction_certified
    assert "kokuno_source_smallness_condition_failed" in diagnostic.rejection_reasons
    assert "fixed_point_lipschitz_bound_exceeds_one_half" in diagnostic.rejection_reasons
    with pytest.raises(MomentCorrectionGainError, match="rejected before iteration"):
        solve_five_moment_correction(backend, Candidate())


def test_zero_discrepancy_does_not_materialize_a_trivial_correction() -> None:
    backend = Backend(discrepancy_scale=0.0)
    diagnostic = assess_five_moment_correction_gain(backend, Candidate())
    assert not diagnostic.correction_needed
    assert "zero_moment_discrepancy" in diagnostic.rejection_reasons
    with pytest.raises(MomentCorrectionGainError, match="zero_moment_discrepancy"):
        solve_five_moment_correction(backend, Candidate())


def test_candidate_identity_and_forcing_shortcut_fail_closed() -> None:
    with pytest.raises(MomentCorrectionGainError, match="candidate hash mismatch"):
        assess_five_moment_correction_gain(Backend(candidate_hash="b" * 64), Candidate())
    with pytest.raises(MomentCorrectionGainError, match="residual-as-forcing"):
        assess_five_moment_correction_gain(Backend(residual_as_forcing=True), Candidate())


def test_singular_or_malformed_moment_system_fails_closed() -> None:
    singular = np.eye(MOMENT_DIMENSION)
    singular[-1] = singular[0]
    with pytest.raises(MomentCorrectionGainError, match="singular"):
        assess_five_moment_correction_gain(Backend(matrix=singular), Candidate())

    base = Backend()

    class BadShapeBackend(Backend):
        def five_moment_correction_system(self, candidate: Candidate) -> FiveMomentCorrectionSystem:
            system = super().five_moment_correction_system(candidate)
            return replace(system, matrix_B=((1.0, 0.0), (0.0, 1.0)))

    with pytest.raises(MomentCorrectionGainError, match="5x5"):
        assess_five_moment_correction_gain(BadShapeBackend(), Candidate())

    assert base.five_moment_correction_system(Candidate()).discrepancy_from_complete_ns_defect is False


def test_public_api_has_no_caller_gain_norm_or_threshold_knobs() -> None:
    forbidden = {
        "B",
        "Q",
        "d",
        "nu",
        "q",
        "gain",
        "radius",
        "lipschitz",
        "smallness",
        "success",
        "residual",
        "defect",
        "pressure",
        "forcing",
        "damping",
        "alpha",
        "threshold",
        "scientific_threshold",
        "viscosity",
    }
    for fn in (assess_five_moment_correction_gain, solve_five_moment_correction):
        assert forbidden.isdisjoint(inspect.signature(fn).parameters)
        assert tuple(inspect.signature(fn).parameters) == ("backend", "candidate")


def test_truth_boundary_keeps_local_gain_separate_from_full_ns_cycle() -> None:
    truth = truth_boundary()
    assert truth["source_reader_head"] == SOURCE_READER_HEAD
    assert truth["source_reader_tex_blob"] == SOURCE_READER_TEX_BLOB
    assert truth["source_equation"] == "B c + Q(c,c) = d"
    assert truth["source_smallness_rule"] == "8*nu^2*q*||d||<=1"
    assert truth["local_five_moment_analytic_contraction_gain_available"]
    assert truth["divergent_local_moment_correction_rejected_before_iteration"]
    assert not truth["analytic_full_ns_cycle_gain_available"]
    assert not truth["current_strict_inner_radial_stress_authorized_as_ns_correction"]
    assert not truth["complete_ns_correction_velocity_materialized"]
    assert not truth["real_candidate_finite_correction_cycle_run"]
    assert not truth["heldout_normalized_ns_residual_assessed"]
    assert not truth["residual_reduction_claimed"]
    assert not truth["same_protocol_comparable_to_st006"]
    assert not truth["pde_validated"]
    assert truth["final_normalized_momentum_gate"] == 1.0e-3
    assert truth["final_normalized_divergence_gate"] == 1.0e-5


def test_mechanics_report_contains_both_accept_and_reject_paths() -> None:
    report = build_mechanics_report()
    assert report["candidate_residual_evidence"] is False
    assert report["accepted_mechanics"]["local_moment_fixed_point_converged"] is True
    assert report["rejected_mechanics_would_be_blocked"] is True
    assert report["rejected_mechanics"]["rejection_reasons"]
