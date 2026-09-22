from __future__ import annotations

import inspect

import pytest

from openai_ns_reconstruction import kokuno_a2_exact_current_i4_runtime_rebind_contract as mod


def test_truth_boundary_is_runtime_only_and_fail_closed() -> None:
    truth = mod.truth_boundary()
    assert truth["a2_1080_runtime_is_pinned_external_repository_snapshot"] is True
    assert truth["a1_1079_runtime_is_pinned_external_repository_snapshot"] is True
    assert truth["a2_960_differential_runtime_is_local_exact_blob"] is True
    assert truth["exact_backend_bind_is_the_authenticator"] is True
    assert truth["oscillatory_complete_curl_reimplemented"] is False
    assert truth["leading_profile_reimplemented"] is False
    assert truth["agent3_mean_correction_reimplemented"] is False
    assert truth["rf30_rf39_correction_executed_here"] is False
    assert truth["cartesian_delta_u_materialized"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["paper_exact"] is False
    assert truth["pde_validated"] is False
    assert truth["forbidden_scientific_controls_exposed"] is False
    assert truth["public_parameters"] == ("composite_field", "differential_function")


def test_nonexact_objects_cannot_cross_rebind_boundary() -> None:
    def fake_differential(x, t):
        return x, t

    with pytest.raises(mod.ExactRuntimeRebindError, match="runtime rebind failed"):
        mod.bind_exact_current_i4_runtime(object(), fake_differential)


def test_receipt_enforcement_rejects_truth_promotion() -> None:
    base = {
        "schema": mod.SCHEMA,
        "task": mod.TASK,
        "composite_source_blob": mod.AGENT2_COMPOSITE_SOURCE_BLOB,
        "leading_source_blob": mod.AGENT1_LEADING_SOURCE_BLOB,
        "differential_source_blob": mod.AGENT2_DIFFERENTIAL_SOURCE_BLOB,
        "exact_backend_rebind_executed": True,
        "exact_external_a2_1080_runtime_authenticated": True,
        "exact_external_a1_1079_runtime_authenticated": True,
        "local_exact_a2_960_runtime_authenticated": True,
        "rf30_rf39_correction_executed_here": False,
        "rf44_postupdate_state_materialized": False,
        "cartesian_delta_u_materialized": False,
        "heldout_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "paper_exact": False,
        "pde_validated": False,
    }
    mod.enforce_receipt(dict(base))
    for key in (
        "rf30_rf39_correction_executed_here",
        "rf44_postupdate_state_materialized",
        "cartesian_delta_u_materialized",
        "heldout_ns_residual_assessed",
        "residual_reduction_claimed",
        "paper_exact",
        "pde_validated",
    ):
        mutated = dict(base)
        mutated[key] = True
        with pytest.raises(AssertionError):
            mod.enforce_receipt(mutated)


def test_contract_exposes_no_scientific_tuning_escape_hatch() -> None:
    params = set(inspect.signature(mod.bind_exact_current_i4_runtime).parameters)
    forbidden = {
        "residual", "defect", "forcing", "pressure", "gain", "optimizer",
        "threshold", "viscosity", "nu", "coefficient", "correction", "target",
    }
    assert not (params & forbidden)
