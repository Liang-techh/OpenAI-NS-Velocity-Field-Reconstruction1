import math

import pytest

from openai_ns_reconstruction.actual_signed_mean_defect import ActualSignedMeanDefectAdmission
from openai_ns_reconstruction.kokuno_actual_mean_factor_dependency import (
    BUMP_NUMERIC_DEPENDENCY,
    MISSING_WEIGHT_DEFINITION,
    PHYSICAL_SCALE_DEFINITION,
    PINNED_MATHLIB_COMMIT,
    KokunoActualMeanFactorDependencyGate,
    dyadic_q,
    missing_weight_indices,
    physical_scale_from_coordinate_q,
)
from test_actual_signed_canonical_scope import _scope
from test_actual_signed_mean_defect import _witness


def _admission():
    scope = _scope()
    return ActualSignedMeanDefectAdmission(scope, _witness(scope))


def test_missing_weight_index_set_matches_closed_integer_interval():
    assert missing_weight_indices(0) == (-1,)
    assert missing_weight_indices(1) == (-1, 0)
    assert missing_weight_indices(4) == (-1, 0, 1, 2, 3)


@pytest.mark.parametrize("value", [-1, True, 1.5])
def test_missing_weight_index_set_rejects_non_naturals(value):
    with pytest.raises(ValueError, match="nonnegative integer"):
        missing_weight_indices(value)


def test_physical_scale_partial_evaluator_uses_exact_dyadic_band_law():
    assert dyadic_q(1) == 0.5
    assert dyadic_q(3) == 0.125
    assert physical_scale_from_coordinate_q(3, 1.6) == pytest.approx(0.2, rel=0, abs=1e-16)


@pytest.mark.parametrize("coordinate_q", [0.0, -1.0, math.inf, math.nan])
def test_physical_scale_rejects_nonpositive_or_nonfinite_coordinate_q(coordinate_q):
    with pytest.raises(ValueError, match="strictly positive and finite"):
        physical_scale_from_coordinate_q(3, coordinate_q)


def test_current_formal_source_factor_is_pinned_but_numeric_debt_stays_closed():
    admission = _admission()
    receipt = KokunoActualMeanFactorDependencyGate().evaluate(admission)

    assert receipt["formal_commit"] == admission.witness.lean_commit
    assert receipt["mathlib_commit"] == PINNED_MATHLIB_COMMIT
    assert receipt["requested_cross_defect_theorem"] == admission.witness.theorem_symbol
    assert receipt["missing_weight_definition"] == MISSING_WEIGHT_DEFINITION
    assert receipt["physical_scale_definition"] == PHYSICAL_SCALE_DEFINITION
    assert receipt["bump_numeric_dependency"] == BUMP_NUMERIC_DEPENDENCY
    assert receipt["prepared_N"] == admission.witness.prepared_N
    assert receipt["band"] == admission.witness.band
    assert receipt["missing_weight_indices"] == missing_weight_indices(admission.witness.prepared_N)
    assert receipt["missing_weight_term_count"] == admission.witness.prepared_N + 1

    assert receipt["source_algebra_pinned"] is True
    assert receipt["physical_scale_algebra_executable_given_coordinate_q"] is True
    assert receipt["missing_weight_depends_on_pinned_mathlib_contdiff_bump"] is True
    assert receipt["mathlib_bump_base_selected_by_nonempty_some"] is True
    assert receipt["numeric_bump_profile_export_present"] is False
    assert receipt["missing_weight_values_materialized"] is False
    assert receipt["requested_stress_actual_state_values_materialized"] is False
    assert receipt["finite_head_mean_debt_materialized"] is False
    assert receipt["arbitrary_python_bump_may_be_promoted_to_actual"] is False
    assert receipt["surrogate_defect_used"] is False
    assert receipt["real_candidate_defect_consumed"] is False
    assert receipt["signed_mean_inverse_input_ready"] is False
    assert receipt["finite_correction_cycle_rerun_allowed"] is False
    assert receipt["heldout_ns_residual_assessed"] is False
    assert receipt["residual_reduction_claimed"] is False
    assert receipt["pde_validated"] is False

    assert set(receipt["blockers"]) == {
        "no machine-linked numeric realization of the pinned dyadic bump/missingWeight is exported",
        "actual cycle-state theta/axial residual requestedStress values remain opaque",
        "finite-head physical mean debt is not numerically materialized",
    }


def test_dependency_gate_rejects_non_admission_input():
    with pytest.raises(TypeError, match="ActualSignedMeanDefectAdmission"):
        KokunoActualMeanFactorDependencyGate().evaluate(object())
