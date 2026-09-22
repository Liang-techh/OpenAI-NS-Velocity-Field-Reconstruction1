import copy
import inspect
import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_relative_swirl_angular_correction import (
    PARENT_EXACT_HEAD,
    SCHEMA,
    SOURCE_BLOB,
    SOURCE_COMMIT,
    KokunoCurrentRelativeSwirlAngularCorrection,
)


def test_source_parent_and_truth_boundary_are_fail_closed():
    adapter = KokunoCurrentRelativeSwirlAngularCorrection()
    cfg = adapter.configuration()
    assert cfg["schema"] == SCHEMA
    assert cfg["parent_exact_head"] == PARENT_EXACT_HEAD
    assert cfg["source_commit"] == SOURCE_COMMIT
    assert cfg["source_blob"] == SOURCE_BLOB
    assert cfg["moment_order"] == ["M", "J", "I", "S", "C_p"]
    assert "delta_target_current=-delta_I" in cfg["source_formulas"]["stable_current_correction"]

    truth = adapter.truth_boundary
    assert truth["current_PA17_I_residual_consumed"] is True
    assert truth["current_lineage_angular_target_correction_materialized"] is True
    assert truth["late_hold_angular_cancellation_avoided"] is True
    assert truth["imported_base_absolute_angular_target_materialized"] is False
    assert truth["current_lineage_angular_entry_I_materialized"] is False
    assert truth["current_relative_swirl_total_target_materialized"] is False
    assert truth["current_cartesian_relative_swirl_composed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_consumes_actual_frozen_PA17_I_row_and_eta_jet():
    adapter = KokunoCurrentRelativeSwirlAngularCorrection()
    eta = np.array([-0.8, -0.2, 0.0, 0.35, 0.9])
    row, row_eta = adapter.post_i1_I_residual_with_eta(eta)
    full, full_eta = adapter.moments.post_i1_residual_normalized_with_eta(eta)
    np.testing.assert_allclose(row, full[..., 2], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(row_eta, full_eta[..., 2], rtol=0.0, atol=0.0)
    assert np.all(np.isfinite(row))
    assert np.all(np.isfinite(row_eta))
    # The frozen 5e-7 gate is not silently promoted to exact restoration.
    assert np.max(np.abs(row)) > 0.0


def test_log_scale_is_exact_PA17_physical_I_scale():
    adapter = KokunoCurrentRelativeSwirlAngularCorrection()
    repair = adapter.moments.repair
    expected = 1.5 * repair.log_X_1 + repair.log_e_1
    assert adapter.log_delta_I_scale == pytest.approx(expected, rel=0.0, abs=0.0)


def test_first_bump_denominator_matches_unedited_hold_identity():
    adapter = KokunoCurrentRelativeSwirlAngularCorrection()
    eta = np.array([-1.0, -0.4, 0.0, 0.55, 1.0])
    log_xh = adapter.log_XH_first_bump_center(eta)

    state = adapter.hold._flatten_endpoint_state(eta)
    E_flat = np.asarray(state["E"], dtype=float)
    expected = (
        1.5 * adapter.hold.log_X_flatten_end
        + 0.5 * math.log(2.0)
        + np.log(E_flat)
        + (1.0 - adapter.lambda_value) * adapter.compensator.y1
    )
    np.testing.assert_allclose(log_xh, expected, rtol=0.0, atol=0.0)
    # The source eta-flattening endpoint is mathematically eta-independent.
    np.testing.assert_allclose(E_flat, E_flat[2], rtol=5.0e-13, atol=0.0)


def test_stable_target_equals_direct_physical_scale_ratio_without_late_subtraction():
    adapter = KokunoCurrentRelativeSwirlAngularCorrection()
    eta = np.array([-0.75, 0.0, 0.6])
    row, row_eta = adapter.post_i1_I_residual_with_eta(eta)
    log_ratio = adapter.log_delta_I_scale - adapter.log_XH_first_bump_center(eta)
    expected = -row * np.exp(log_ratio)
    expected_eta = -row_eta * np.exp(log_ratio)

    target, target_eta = adapter.correction_target_with_eta(eta)
    np.testing.assert_allclose(target, expected, rtol=2.0e-15, atol=0.0)
    np.testing.assert_allclose(target_eta, expected_eta, rtol=2.0e-15, atol=0.0)
    assert np.max(np.abs(target)) > 0.0


def test_current_correction_survives_when_late_rI_subtraction_rounds_to_zero():
    adapter = KokunoCurrentRelativeSwirlAngularCorrection()
    target, _ = adapter.correction_target_with_eta(0.0)
    target = float(target)
    assert target != 0.0

    alpha = 1.0 - adapter.lambda_value
    delta_r_end = -target * math.exp(
        -alpha * (adapter.compensator.hold_length - adapter.compensator.y1)
    )
    r_star = 1.0 / alpha
    # This is exactly why the adapter carries the row-scaled defect instead of
    # reconstructing it by subtracting two O(1) doubles at the late hold.
    assert r_star + delta_r_end == r_star


def test_eta_jet_matches_centered_replay_of_the_executable_correction():
    adapter = KokunoCurrentRelativeSwirlAngularCorrection()
    eta = 0.31
    target, target_eta = adapter.correction_target_with_eta(eta)
    h = 2.0e-6
    plus, _ = adapter.correction_target_with_eta(eta + h)
    minus, _ = adapter.correction_target_with_eta(eta - h)
    fd = float((plus - minus) / (2.0 * h))
    analytic = float(target_eta)
    scale = max(abs(fd), abs(analytic), 1.0e-300)
    assert abs(fd - analytic) / scale < 2.0e-5
    assert math.isfinite(float(target))


def test_base_target_must_be_explicit_and_is_only_combined_not_invented():
    adapter = KokunoCurrentRelativeSwirlAngularCorrection()
    signature = inspect.signature(adapter.apply_to_materialized_base_target)
    assert signature.parameters["base_target"].default is inspect.Parameter.empty
    assert signature.parameters["base_target_eta"].default is inspect.Parameter.empty

    eta = np.array([-0.4, 0.2])
    base = np.array([2.0e-8, -3.0e-8])
    base_eta = np.array([1.0e-9, 4.0e-9])
    correction, correction_eta = adapter.correction_target_with_eta(eta)
    total, total_eta = adapter.apply_to_materialized_base_target(
        eta, base_target=base, base_target_eta=base_eta
    )
    np.testing.assert_allclose(total, base + correction, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(total_eta, base_eta + correction_eta, rtol=0.0, atol=0.0)


def test_save_load_semantic_identity_and_truth_mutation_fail_closed(tmp_path):
    adapter = KokunoCurrentRelativeSwirlAngularCorrection()
    path = tmp_path / "current_angular_correction.json"
    adapter.save_configuration(path)
    loaded = KokunoCurrentRelativeSwirlAngularCorrection.load_configuration(path)
    assert loaded.semantic_sha256() == adapter.semantic_sha256()
    assert loaded.configuration() == adapter.configuration()

    payload = json.loads(path.read_text())
    mutated = copy.deepcopy(payload)
    mutated["truth_boundary"]["current_lineage_angular_entry_I_materialized"] = True
    with pytest.raises(ValueError):
        KokunoCurrentRelativeSwirlAngularCorrection.from_configuration(mutated)

    mutated = copy.deepcopy(payload)
    mutated["numerical_realization"]["missing_absolute_anchor_policy"] = "assume zero"
    with pytest.raises(ValueError):
        KokunoCurrentRelativeSwirlAngularCorrection.from_configuration(mutated)


def test_public_api_has_no_residual_forcing_or_optimizer_tuning_surface():
    adapter = KokunoCurrentRelativeSwirlAngularCorrection()
    forbidden = {
        "residual",
        "forcing",
        "pressure",
        "viscosity",
        "optimizer",
        "heldout",
        "threshold",
        "target_residual",
    }
    for method_name in (
        "correction_target_with_eta",
        "apply_to_materialized_base_target",
        "log_XH_first_bump_center",
    ):
        params = set(inspect.signature(getattr(adapter, method_name)).parameters)
        assert not (params & forbidden)
