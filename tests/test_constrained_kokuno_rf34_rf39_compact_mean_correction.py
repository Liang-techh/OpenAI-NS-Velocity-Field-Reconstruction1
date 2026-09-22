from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_rf34_rf39_compact_mean_correction import (
    COEFFICIENT_SOLVE_REL_TOL,
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    FIVE_ROW_REL_TOL,
    RF34RF39ContractError,
    _MechanicsBackend,
    _MechanicsCandidate,
    materialize_rf34_rf39_compact_mean_correction,
    truth_boundary,
)


def test_mechanics_source_closed_form_closes_all_five_rows() -> None:
    receipt = materialize_rf34_rf39_compact_mean_correction(
        _MechanicsBackend(), _MechanicsCandidate()
    )
    assert receipt.source_rf34_rf39_correction_materialized is True
    assert receipt.repository_candidate_correction_evidence is False
    assert receipt.compact_support_preserved is True
    assert receipt.two_zero_moments_preserved is True
    assert receipt.five_row_closure_max_relative <= FIVE_ROW_REL_TOL
    assert receipt.coefficient_direct_solve_max_relative <= COEFFICIENT_SOLVE_REL_TOL
    assert receipt.basis_match_max_abs <= 2.0e-11
    assert receipt.reserved_base_G_max_abs <= 2.0e-11
    assert receipt.reserved_base_V_max_relative_error <= 2.0e-8
    assert receipt.correction_chart_l2 > 0.0
    assert receipt.correction_chart_max > 0.0
    assert np.all(np.isfinite(receipt.coefficients_u0_u1_u2_s0_s1))
    assert receipt.correction_applied_to_candidate is False
    assert receipt.rf44_rf49_nonlinear_remainder_recomputed is False
    assert receipt.complete_ns_defect is False
    assert receipt.heldout_ns_residual_assessed is False
    assert receipt.pde_validated is False


def test_lambda_zero_source_degeneracy_fails_closed() -> None:
    with pytest.raises(RF34RF39ContractError, match="lambda > 0"):
        materialize_rf34_rf39_compact_mean_correction(
            _MechanicsBackend(lambda_value=0.0), _MechanicsCandidate()
        )


def test_canonical_basis_drift_fails_closed() -> None:
    with pytest.raises(RF34RF39ContractError, match="not the frozen RF34"):
        materialize_rf34_rf39_compact_mean_correction(
            _MechanicsBackend(mutate_basis=True), _MechanicsCandidate()
        )


def test_heldout_selected_context_fails_closed() -> None:
    class BadBackend(_MechanicsBackend):
        def rf34_rf39_patch_context(self, candidate, state):
            from dataclasses import replace

            return replace(
                super().rf34_rf39_patch_context(candidate, state),
                heldout_samples_used_to_choose_context=True,
            )

    with pytest.raises(RF34RF39ContractError, match="held-out"):
        materialize_rf34_rf39_compact_mean_correction(BadBackend(), _MechanicsCandidate())


def test_source_exact_bump_laundering_fails_closed() -> None:
    class BadBackend(_MechanicsBackend):
        def rf34_rf39_patch_context(self, candidate, state):
            from dataclasses import replace

            return replace(
                super().rf34_rf39_patch_context(candidate, state),
                source_exact_bump_claimed=True,
            )

    with pytest.raises(RF34RF39ContractError, match="autonomous"):
        materialize_rf34_rf39_compact_mean_correction(BadBackend(), _MechanicsCandidate())


def test_patch_must_be_frozen_before_defect() -> None:
    class BadBackend(_MechanicsBackend):
        def rf34_rf39_patch_context(self, candidate, state):
            from dataclasses import replace

            return replace(
                super().rf34_rf39_patch_context(candidate, state),
                context_frozen_before_defect_evaluation=False,
            )

    with pytest.raises(RF34RF39ContractError, match="frozen before defect"):
        materialize_rf34_rf39_compact_mean_correction(BadBackend(), _MechanicsCandidate())


def test_public_entry_point_exposes_no_scientific_tuning_knobs() -> None:
    params = tuple(inspect.signature(materialize_rf34_rf39_compact_mean_correction).parameters)
    assert params == ("backend", "candidate")


def test_truth_boundary_keeps_real_ns_claims_false() -> None:
    truth = truth_boundary()
    assert truth["rf34_rf39_compact_correction_realization_materialized"] is True
    assert truth["concrete_bump_realization_classification"] == "autonomous_repository_choice"
    assert truth["source_exact_bump_recovered"] is False
    for key in (
        "current_i4_source_chart_backend_materialized",
        "current_i4_rf30_defect_materialized",
        "current_i4_rf34_rf39_correction_materialized",
        "correction_applied_to_candidate",
        "rf44_rf49_nonlinear_remainder_recomputed",
        "cartesian_correction_velocity_materialized",
        "complete_ns_defect",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
    ):
        assert truth[key] is False, key
    assert truth["final_normalized_momentum_gate"] == FINAL_NORMALIZED_MOMENTUM_GATE == 1.0e-3
    assert truth["final_normalized_divergence_gate"] == FINAL_NORMALIZED_DIVERGENCE_GATE == 1.0e-5
    assert truth["residual_defined_free_forcing_allowed"] is False
    assert truth["finite_stage_small_residual_is_blowup_proof"] is False
