from __future__ import annotations

import inspect

import pytest

from openai_ns_reconstruction.kokuno_finite_correction_cycle_contract import (
    CorrectionFieldProvider,
    _analytic_cycle,
    admit_finite_correction_transition,
    audit_finite_correction_transition,
    deterministic_receipt,
)
from openai_ns_reconstruction.kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    VectorFieldProvider,
)


HELD_IN = [
    (-0.31, 0.27, 0.14, 0.37),
    (0.22, -0.19, -0.33, 0.63),
    (0.17, 0.41, -0.12, 0.46),
]
HELD_OUT = [
    (0.11, 0.36, -0.21, 0.41),
    (-0.28, -0.17, 0.29, 0.59),
    (0.33, -0.44, 0.18, 0.52),
]
UPDATE_POINTS = [
    (0.09, 0.24, -0.16, 0.35),
    (-0.21, -0.32, 0.27, 0.55),
    (0.25, 0.38, 0.11, 0.67),
]


def _audit(scale_after: float):
    return audit_finite_correction_transition(
        *_analytic_cycle(scale_after),
        HELD_IN,
        HELD_OUT,
        UPDATE_POINTS,
        viscosity=0.01,
        spatial_step=1.0e-3,
    )


def test_actual_raw_defect_gain_is_recomputed_and_half_step_passes() -> None:
    report = _audit(0.5)
    assert report.update_velocity_closure_max <= 2.0e-11
    assert report.update_velocity_dt_closure_max <= 2.0e-11
    assert report.correction_nontrivial is True
    assert report.correction_vector_rms > 0.0
    assert report.held_in_rms_ratio == pytest.approx(0.5, rel=0.0, abs=2.0e-10)
    assert report.held_out_rms_ratio == pytest.approx(0.5, rel=0.0, abs=2.0e-10)
    assert report.held_in_max_ratio == pytest.approx(0.5, rel=0.0, abs=2.0e-10)
    assert report.held_out_max_ratio == pytest.approx(0.5, rel=0.0, abs=2.0e-10)
    assert report.finite_step_gain_guard_passed is True
    assert report.correction_divergence_max <= 2.0e-11


def test_admission_api_rejects_actual_defect_expansion() -> None:
    report = _audit(1.5)
    assert report.held_in_rms_ratio == pytest.approx(1.5, rel=0.0, abs=3.0e-10)
    assert report.held_out_rms_ratio == pytest.approx(1.5, rel=0.0, abs=3.0e-10)
    assert report.finite_step_gain_guard_passed is False
    with pytest.raises(ValueError, match="gain guard rejected"):
        admit_finite_correction_transition(
            *_analytic_cycle(1.5),
            HELD_IN,
            HELD_OUT,
            UPDATE_POINTS,
            viscosity=0.01,
            spatial_step=1.0e-3,
        )


def test_rejects_transition_identity_skip() -> None:
    providers = list(_analytic_cycle(0.5))
    correction = providers[-1]
    assert isinstance(correction, CorrectionFieldProvider)
    wrong_after = CycleIdentity(
        correction.to_identity.cycle_id,
        correction.to_identity.cycle_index + 1,
        "skipped-after",
    )
    providers[4] = VectorFieldProvider(
        wrong_after,
        providers[4].source_ref,
        providers[4].evaluator,
    )
    with pytest.raises(ValueError, match="one cycle identity|advance one cycle index"):
        audit_finite_correction_transition(
            *providers,
            HELD_IN,
            HELD_OUT,
            UPDATE_POINTS,
            viscosity=0.01,
            spatial_step=1.0e-3,
        )


def test_rejects_restricted_forcing_receipt_drift() -> None:
    providers = list(_analytic_cycle(0.5))
    after_force = providers[7]
    assert isinstance(after_force, RestrictedForcingProvider)
    providers[7] = RestrictedForcingProvider(
        after_force.identity,
        after_force.source_ref,
        "different restriction receipt",
        after_force.evaluator,
    )
    with pytest.raises(ValueError, match="forcing receipt changed"):
        audit_finite_correction_transition(
            *providers,
            HELD_IN,
            HELD_OUT,
            UPDATE_POINTS,
            viscosity=0.01,
            spatial_step=1.0e-3,
        )


def test_rejects_velocity_update_not_equal_to_before_plus_correction() -> None:
    providers = list(_analytic_cycle(0.5))
    after_velocity = providers[4]
    assert isinstance(after_velocity, VectorFieldProvider)
    providers[4] = VectorFieldProvider(
        after_velocity.identity,
        after_velocity.source_ref + ":mutated",
        lambda x, y, z, t: (0.6 * t * y, 0.0, 0.0),
    )
    with pytest.raises(ValueError, match=r"u\^\(k\+1\)=u\^k\+delta_u_k closure failed"):
        audit_finite_correction_transition(
            *providers,
            HELD_IN,
            HELD_OUT,
            UPDATE_POINTS,
            viscosity=0.01,
            spatial_step=1.0e-3,
        )


def test_public_audit_has_no_surrogate_or_gain_threshold_inputs() -> None:
    names = set(inspect.signature(audit_finite_correction_transition).parameters)
    forbidden = {
        "residual",
        "defect",
        "stress",
        "target",
        "gain",
        "gain_threshold",
        "normalized_residual",
        "momentum_threshold",
        "divergence_threshold",
    }
    assert names.isdisjoint(forbidden)
    assert "before_velocity" in names
    assert "after_velocity" in names
    assert "correction" in names


def test_deterministic_receipt_keeps_real_candidate_claims_false() -> None:
    receipt = deterministic_receipt()
    report = receipt["analytic_mechanics_regression"]
    assert report["finite_step_gain_guard_passed"] is True
    assert report["held_in_rms_ratio"] == pytest.approx(0.5, abs=2.0e-10)
    assert report["held_out_rms_ratio"] == pytest.approx(0.5, abs=2.0e-10)
    boundary = receipt["truth_boundary"]
    assert boundary["finite_correction_cycle_mechanics_executable"] is True
    assert boundary["before_after_actual_defect_recomputed_from_raw_providers"] is True
    assert boundary["divergent_correction_rejected_by_admission_api"] is True
    assert boundary["formal_kokuno_correction_cycle_gain_bound_claimed"] is False
    assert boundary["restricted_forcing_semantics_independently_validated_here"] is False
    assert boundary["real_full_candidate_defect_consumed"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed_for_real_candidate"] is False
    assert boundary["pde_validated"] is False
