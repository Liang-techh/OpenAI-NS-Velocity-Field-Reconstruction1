from __future__ import annotations

from dataclasses import replace
import inspect
import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_rf44_rf49_postupdate_recompute import (
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    RF44RF49ContractError,
    RF47_REL_TOL,
    _RF44MechanicsBackend,
    materialize_rf44_rf49_postupdate_recompute,
    truth_boundary,
)


def test_mechanics_recomputes_nonzero_postupdate_defect_and_rf47_closes() -> None:
    r = materialize_rf44_rf49_postupdate_recompute(_RF44MechanicsBackend(), None)
    assert r.source_rf44_rf49_recompute_materialized is True
    assert r.repository_candidate_recompute_evidence is False
    assert r.current_i4_candidate_evidence is False
    assert r.rf47_closure_max_relative <= RF47_REL_TOL
    assert np.isfinite(r.defect_gain_ratio_2) and r.defect_gain_ratio_2 >= 0.0
    assert r.nonlinear_remainder_nonzero is True
    assert r.post_defect.defect_norm_2 > 0.0
    assert np.max(np.abs(r.radial_remainder_Rg)) > 0.0
    assert r.cartesian_correction_velocity_materialized is False
    assert r.complete_ns_defect is False
    assert r.heldout_ns_residual_assessed is False
    assert r.pde_validated is False


def test_heldout_constructed_post_state_fails_closed() -> None:
    class Bad(_RF44MechanicsBackend):
        def rf44_post_update_state(self, candidate, correction):
            return replace(super().rf44_post_update_state(candidate, correction),
                           heldout_samples_used_to_construct_state=True)
    with pytest.raises(RF44RF49ContractError, match="held-out"):
        materialize_rf44_rf49_postupdate_recompute(Bad(), None)


def test_surrogate_increment_fails_closed() -> None:
    class Bad(_RF44MechanicsBackend):
        def rf44_correction_increment(self, candidate, correction, pre_state, post_state):
            return replace(super().rf44_correction_increment(candidate, correction, pre_state, post_state),
                           surrogate_increment_used=True)
    with pytest.raises(RF44RF49ContractError, match="held-out/surrogate"):
        materialize_rf44_rf49_postupdate_recompute(Bad(), None)


def test_increment_must_match_exact_rf34_rf39_correction() -> None:
    class Bad(_RF44MechanicsBackend):
        def rf44_correction_increment(self, candidate, correction, pre_state, post_state):
            inc = super().rf44_correction_increment(candidate, correction, pre_state, post_state)
            u = np.asarray(inc.delta_v, dtype=float).copy(); u[len(u)//2] += 1e-5
            return replace(inc, delta_v=tuple(float(x) for x in u))
    with pytest.raises(RF44RF49ContractError, match="does not match exact RF34"):
        materialize_rf44_rf49_postupdate_recompute(Bad(), None)


def test_wave_or_base_drift_fails_closed() -> None:
    class Bad(_RF44MechanicsBackend):
        def rf44_post_update_state(self, candidate, correction):
            s = super().rf44_post_update_state(candidate, correction)
            V = np.asarray(s.base_V, dtype=float).copy(); V[len(V)//2] += 1e-6
            return replace(s, base_V=tuple(float(x) for x in V))
    with pytest.raises(RF44RF49ContractError, match="fixed wave/base"):
        materialize_rf44_rf49_postupdate_recompute(Bad(), None)


def test_post_state_must_equal_authenticated_increment() -> None:
    class Bad(_RF44MechanicsBackend):
        def rf44_post_update_state(self, candidate, correction):
            s = super().rf44_post_update_state(candidate, correction)
            v = np.asarray(s.v, dtype=float).copy(); v[len(v)//3] += 1e-6
            return replace(s, v=tuple(float(x) for x in v))
    with pytest.raises(RF44RF49ContractError, match="post-pre v"):
        materialize_rf44_rf49_postupdate_recompute(Bad(), None)


def test_operator_terms_must_be_recomputed_from_state() -> None:
    class Bad(_RF44MechanicsBackend):
        def rf44_post_update_state(self, candidate, correction):
            return replace(super().rf44_post_update_state(candidate, correction),
                           operator_terms_recomputed_from_state=False)
    with pytest.raises(RF44RF49ContractError, match="provenance is incomplete"):
        materialize_rf44_rf49_postupdate_recompute(Bad(), None)


def test_public_entry_point_has_no_defect_gain_or_threshold_knobs() -> None:
    assert tuple(inspect.signature(materialize_rf44_rf49_postupdate_recompute).parameters) == ("backend", "candidate")


def test_truth_boundary_keeps_real_ns_claims_false_and_gates_fixed() -> None:
    t = truth_boundary()
    assert t["rf44_rf49_postupdate_recompute_contract_materialized"] is True
    assert t["source_specific_defect_gain_diagnostic_materialized"] is True
    assert t["new_gain_acceptance_threshold_introduced"] is False
    for key in (
        "current_i4_source_chart_backend_materialized", "current_i4_rf30_defect_materialized",
        "current_i4_rf34_rf39_correction_materialized", "current_i4_rf44_rf49_remainder_recomputed",
        "cartesian_correction_velocity_materialized", "complete_ns_defect",
        "real_candidate_finite_correction_cycle_run", "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed", "same_protocol_comparable_to_st006", "pde_validated",
    ):
        assert t[key] is False, key
    assert t["final_normalized_momentum_gate"] == FINAL_NORMALIZED_MOMENTUM_GATE == 1e-3
    assert t["final_normalized_divergence_gate"] == FINAL_NORMALIZED_DIVERGENCE_GATE == 1e-5
    assert t["residual_defined_free_forcing_allowed"] is False
    assert t["finite_stage_small_residual_is_blowup_proof"] is False
