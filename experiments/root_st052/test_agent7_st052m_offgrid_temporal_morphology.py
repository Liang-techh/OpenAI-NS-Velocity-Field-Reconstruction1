from __future__ import annotations

import numpy as np

import agent7_st052m_offgrid_temporal_morphology as audit


def test_frozen_protocol_is_exactly_one_offgrid_coordinate() -> None:
    assert audit.TASK_ID == "CR003-ST052M-OFFGRID-TEMPORAL-MORPHOLOGY-139"
    assert audit.PREREG_ISSUE == 1177
    assert audit.STACK_BASE_PR == 1167
    assert audit.STACK_BASE_HEAD == "3fc6286c4b338126ce60eb167ed082a6af953335"
    assert audit.COORDINATE_NAME == "core_width_skew_0375_to_0625"
    assert (audit.T_LEFT, audit.T_RIGHT) == (0.375, 0.625)
    assert audit.SUPPORT_RADIUS_SCALE == 2.0
    assert audit.EPSILONS == (1.0e-3, 5.0e-4)
    assert audit.PROBE_NAME == "time_skewed_poloidal_extension"


def test_parent_time_skew_activation_is_visible_offgrid() -> None:
    activation = audit.alias.rep._activation_skew
    assert activation(0.25) == 0.0
    assert activation(0.50) == 0.0
    assert activation(0.75) == 0.0
    assert np.isclose(activation(0.375), -0.375, rtol=0.0, atol=1.0e-15)
    assert np.isclose(activation(0.625), 0.375, rtol=0.0, atol=1.0e-15)


def test_centered_coordinate_derivative_on_synthetic_time_skew(monkeypatch) -> None:
    point = np.asarray([[0.2, 0.0, 0.1]], dtype=float)

    def synthetic_coordinate(field_fn):
        right = float(field_fn(point, audit.T_RIGHT)[0, 0])
        left = float(field_fn(point, audit.T_LEFT)[0, 0])
        return right - left

    monkeypatch.setattr(audit, "core_width_skew_coordinate", synthetic_coordinate)

    def base(points, time):
        pts = np.asarray(points, dtype=float)
        return np.zeros((pts.shape[0], 3), dtype=float)

    def tangent(points, time):
        pts = np.asarray(points, dtype=float)
        out = np.zeros((pts.shape[0], 3), dtype=float)
        out[:, 0] = float(time)
        return out

    derivative = audit.coordinate_derivative(base, tangent, 5.0e-4)
    assert np.isclose(derivative, audit.T_RIGHT - audit.T_LEFT, rtol=0.0, atol=1.0e-12)


def test_derivative_gate_requires_stable_nonzero_same_sign_response() -> None:
    good = audit._derivative_gate(0.2005, 0.2)
    assert good["response_nonzero"] is True
    assert good["coarse_fine_same_nonzero_sign"] is True
    assert good["stable_nonzero_response"] is True

    zero = audit._derivative_gate(0.0, 0.0)
    assert zero["response_nonzero"] is False
    assert zero["stable_nonzero_response"] is False

    sign_flip = audit._derivative_gate(-0.2, 0.2)
    assert sign_flip["coarse_fine_same_nonzero_sign"] is False
    assert sign_flip["stable_nonzero_response"] is False

    drifting = audit._derivative_gate(0.3, 0.2)
    assert drifting["relative_derivative_drift"] > audit.DERIVATIVE_DRIFT_MAX
    assert drifting["stable_nonzero_response"] is False


def test_rank_helper_detects_one_new_direction() -> None:
    existing = np.zeros((9, 5), dtype=float)
    existing[:5, :] = np.eye(5)
    independent = np.zeros(9, dtype=float)
    independent[5] = 1.0
    duplicate = existing[:, 0].copy()

    rank5, condition5 = audit._column_normalized_rank(existing)
    rank6, condition6 = audit._column_normalized_rank(
        np.column_stack((existing, independent))
    )
    rank_dup, _condition_dup = audit._column_normalized_rank(
        np.column_stack((existing, duplicate))
    )

    assert rank5 == 5
    assert rank6 == 6
    assert rank_dup == 5
    assert np.isclose(condition5, 1.0)
    assert np.isclose(condition6, 1.0)


def test_truth_boundary_forbids_candidate_promotion() -> None:
    truth = audit.TRUTH
    assert truth["candidate_velocity_changed"] is False
    assert truth["basis_dimension_changed"] is False
    assert truth["new_temporal_basis_added_to_candidate"] is False
    assert truth["coefficient_selected"] is False
    assert truth["pressure_or_force_changed"] is False
    assert truth["free_residual_force_used"] is False
    assert truth["held_out_pde_residual_evaluated"] is False
    assert truth["public_image_numeric_target_used"] is False
    assert truth["direct_visualization_fingerprint_improvement"] == 0.0
    assert truth["visualization_ready"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
