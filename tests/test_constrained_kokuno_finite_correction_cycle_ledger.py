from __future__ import annotations

import inspect

import pytest

from openai_ns_reconstruction.kokuno_finite_correction_cycle_ledger import (
    _analytic_cycle,
    _regression_points,
    admit_finite_correction_cycle,
    audit_finite_correction_cycle,
    deterministic_receipt,
)


def test_deterministic_receipt_records_two_step_actual_defect_contraction() -> None:
    receipt = deterministic_receipt()
    assert receipt["schema"] == "kokuno-a3-finite-correction-cycle-ledger-v1"
    assert receipt["task"] == "KOKUNO-A3-FINITE-CYCLE-LEDGER-048"
    assert (
        receipt["provenance"]["parent_agent3_head"]
        == "f7830a34784127a9ad7894047039b0ea1e8795eb"
    )

    regression = receipt["analytic_mechanics_regression"]
    report = regression["report"]
    assert report["state_count"] == 3
    assert report["step_count"] == 2
    assert report["every_step_nontrivial"] is True
    assert report["every_step_actual_defect_gain_guard_passed"] is True
    assert report["empirical_cycle_contraction_passed"] is True
    assert report["held_in_rms_contraction_factors"] == pytest.approx([0.5, 0.5], abs=2e-10)
    assert report["held_out_rms_contraction_factors"] == pytest.approx([0.5, 0.5], abs=2e-10)
    assert report["cumulative_held_in_rms_ratio"] == pytest.approx(0.25, abs=2e-10)
    assert report["cumulative_held_out_rms_ratio"] == pytest.approx(0.25, abs=2e-10)
    assert report["worst_step_held_out_rms_ratio"] == pytest.approx(0.5, abs=2e-10)
    assert regression["expansion_mutation_rejected"] is True


def test_public_cycle_admission_recomputes_every_adjacent_transition() -> None:
    held_in, held_out, update = _regression_points()
    states, corrections = _analytic_cycle((1.0, 0.5, 0.25))
    report = admit_finite_correction_cycle(
        states,
        corrections,
        held_in,
        held_out,
        update,
        viscosity=0.01,
        spatial_step=1.0e-3,
    )
    assert report.step_count == 2
    assert len(report.step_reports) == 2
    assert report.step_reports[0].to_identity == report.step_reports[1].from_identity
    assert report.cumulative_held_out_rms_ratio == pytest.approx(0.25, abs=2e-10)
    assert report.maximum_correction_divergence_max <= 2e-11


def test_divergent_intermediate_or_later_step_is_rejected() -> None:
    held_in, held_out, update = _regression_points()
    states, corrections = _analytic_cycle((1.0, 0.5, 0.75))
    audit = audit_finite_correction_cycle(
        states,
        corrections,
        held_in,
        held_out,
        update,
        viscosity=0.01,
        spatial_step=1.0e-3,
    )
    assert audit.step_reports[0].finite_step_gain_guard_passed is True
    assert audit.step_reports[1].finite_step_gain_guard_passed is False
    assert audit.every_step_actual_defect_gain_guard_passed is False
    assert audit.empirical_cycle_contraction_passed is False

    with pytest.raises(ValueError, match="actual-defect gain guard"):
        admit_finite_correction_cycle(
            states,
            corrections,
            held_in,
            held_out,
            update,
            viscosity=0.01,
            spatial_step=1.0e-3,
        )


def test_cycle_shape_rejects_missing_or_mismatched_corrections() -> None:
    held_in, held_out, update = _regression_points()
    states, corrections = _analytic_cycle((1.0, 0.5, 0.25))
    with pytest.raises(ValueError, match="exactly one correction"):
        audit_finite_correction_cycle(
            states,
            corrections[:1],
            held_in,
            held_out,
            update,
            viscosity=0.01,
            spatial_step=1.0e-3,
        )


def test_public_api_has_no_surrogate_or_threshold_injection_parameters() -> None:
    forbidden = (
        "residual",
        "defect",
        "stress",
        "target",
        "gain",
        "normalized",
        "threshold",
    )
    for fn in (audit_finite_correction_cycle, admit_finite_correction_cycle):
        names = tuple(inspect.signature(fn).parameters)
        assert all(not any(word in name for word in forbidden) for name in names)


def test_truth_boundary_does_not_promote_real_candidate_or_pde_claims() -> None:
    boundary = deterministic_receipt()["truth_boundary"]
    assert boundary["finite_correction_cycle_ledger_executable"] is True
    assert boundary["each_step_actual_defect_recomputed_from_raw_providers"] is True
    assert boundary["one_disjoint_heldin_heldout_split_preserved_across_cycle"] is True
    assert boundary["empirical_actual_defect_contraction_diagnostic_executable"] is True
    assert boundary["divergent_intermediate_step_rejected"] is True
    assert boundary["formal_kokuno_correction_cycle_gain_bound_claimed"] is False
    assert boundary["real_full_candidate_defect_consumed"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed_for_real_candidate"] is False
    assert boundary["pde_validated"] is False
    assert boundary["final_normalized_momentum_gate"] == 1e-3
    assert boundary["final_normalized_divergence_gate"] == 1e-5
