import pytest

from openai_ns_reconstruction.actual_signed_mean_defect import ActualSignedMeanDefectAdmission
from openai_ns_reconstruction.kokuno_actual_signed_mean_correction_gate import (
    FORMAL_NUMERIC_EVIDENCE_KIND,
    PHYSICAL_TARGET_UNITS,
    ActualSignedMeanNumericTarget,
    KokunoActualSignedMeanCorrectionGate,
)
from test_actual_signed_canonical_scope import _scope
from test_actual_signed_mean_defect import _witness


def _admission():
    scope = _scope()
    return ActualSignedMeanDefectAdmission(scope, _witness(scope))


def _signed_rank_receipt(**changes):
    receipt = {
        "local_supplied_signed_physical_covariance_rank_two": True,
        "physical_complete_curl_covariance_rank_two_assessed": True,
        "actual_positive_order_background_bound": False,
        "actual_source_h_sigma_pulse_integrals_bound": False,
        "actual_signed_auxiliary_rectangles_bound": False,
        "actual_auxiliary_torus_mode_family_bound": False,
        "source_actual_partition_labels_instantiated": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "real_candidate_defect_consumed": False,
    }
    receipt.update(changes)
    return receipt


def _target(admission, **changes):
    kwargs = dict(
        application_id=admission.witness.application_id,
        delta_c_theta=(0.01, -0.02),
        delta_c_z=(0.03, -0.04),
        epsilon=(0.04, 0.01),
        provenance=admission.witness.provenance,
    )
    kwargs.update(changes)
    return ActualSignedMeanNumericTarget(**kwargs)


def test_current_formal_mean_debt_and_supplied_rank_do_not_unlock_numeric_inverse():
    admission = _admission()
    result = KokunoActualSignedMeanCorrectionGate().evaluate(
        admission, _signed_rank_receipt()
    )

    assert result["requested_cross_defect_identity_admitted"] is True
    assert result["finite_head_band_admitted"] is True
    assert result["requested_cross_defect_theorem_machine_replayed"] is False
    assert result["theorem_values_materialized"] is False
    assert result["supplied_signed_physical_covariance_rank_two"] is True
    assert result["actual_signed_complete_curl_source_family_bound"] is False
    assert result["actual_source_physical_covariance_rank_two"] is False
    assert result["real_numeric_target_bound"] is False
    assert result["signed_mean_inverse_input_ready"] is False
    assert result["real_candidate_defect_consumed"] is False
    assert result["finite_correction_cycle_rerun_allowed"] is False
    assert result["heldout_ns_residual_assessed"] is False
    assert result["residual_reduction_claimed"] is False
    assert result["pde_validated"] is False
    assert result["physical_to_reference_conversion"] == "H_ref y = Delta C / epsilon"
    assert set(result["blockers"]) == {
        "actual finite-head mean debt values remain opaque",
        "actual signed complete-curl source family is not bound",
        "actual-source physical covariance rank two is not certified",
        "no theorem-bound numeric physical mean-defect target is available",
    }


def test_hand_filled_numeric_target_is_rejected_while_formal_values_are_opaque():
    admission = _admission()
    target = _target(admission)
    with pytest.raises(ValueError, match="still marks the actual .* values opaque"):
        KokunoActualSignedMeanCorrectionGate().evaluate(
            admission,
            _signed_rank_receipt(
                actual_positive_order_background_bound=True,
                actual_source_h_sigma_pulse_integrals_bound=True,
                actual_signed_auxiliary_rectangles_bound=True,
                actual_auxiliary_torus_mode_family_bound=True,
                source_actual_partition_labels_instantiated=True,
                genuinely_independent_second_covariance_column_ready=True,
            ),
            numeric_target=target,
        )


@pytest.mark.parametrize(
    "changes, message",
    [
        ({"surrogate_defect_used": True}, "surrogate mean-defect targets are forbidden"),
        ({"epsilon": (0.04, 0.0)}, "epsilon must be strictly positive"),
        ({"evidence_kind": "numeric-scan"}, "formal-machine-replay"),
        ({"units": "reference covariance"}, "physical phase-mean covariance/stress"),
    ],
)
def test_numeric_target_contract_rejects_surrogate_or_unit_provenance_drift(changes, message):
    admission = _admission()
    with pytest.raises(ValueError, match=message):
        _target(admission, **changes)


def test_numeric_target_contract_records_exact_expected_units_and_evidence_kind():
    admission = _admission()
    target = _target(admission)
    assert target.units == PHYSICAL_TARGET_UNITS
    assert target.evidence_kind == FORMAL_NUMERIC_EVIDENCE_KIND
    assert target.surrogate_defect_used is False


def test_missing_signed_rank_boundary_fact_fails_closed():
    admission = _admission()
    receipt = _signed_rank_receipt()
    receipt.pop("actual_signed_auxiliary_rectangles_bound")
    with pytest.raises(ValueError, match="missing required keys"):
        KokunoActualSignedMeanCorrectionGate().evaluate(admission, receipt)
