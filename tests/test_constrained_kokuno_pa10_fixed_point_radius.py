from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_pa10_fixed_point_radius import (
    KokunoPA10FixedPointRadiusGate,
)


def test_selected_lambda_is_executable_and_radius_scales_as_source_formula() -> None:
    gate = KokunoPA10FixedPointRadiusGate(operator_constant_M=3.0)
    assert gate.selected_lambda > 1.0e20
    assert gate.conditional_radius == pytest.approx(3.0 / gate.selected_lambda)
    assert gate.radius_for_M(6.0) == pytest.approx(2.0 * gate.conditional_radius)


def test_missing_source_M_fails_closed_without_inventing_number() -> None:
    gate = KokunoPA10FixedPointRadiusGate()
    assert gate.conditional_radius is None
    report = gate.report()
    assert report["selected_execution"]["operator_constant_M_input"] is None
    assert report["selected_execution"]["conditional_fixed_point_radius_M_over_Lambda"] is None


def test_M_budget_is_exact_inverse_planning_relation() -> None:
    gate = KokunoPA10FixedPointRadiusGate()
    target = 1.0e-6
    budget = gate.M_budget_for_radius(target)
    assert budget == pytest.approx(target * gate.selected_lambda)
    assert gate.radius_for_M(budget) == pytest.approx(target)


def test_invalid_M_and_target_rejected() -> None:
    with pytest.raises(ValueError, match="operator_constant_M"):
        KokunoPA10FixedPointRadiusGate(operator_constant_M=0.0)
    gate = KokunoPA10FixedPointRadiusGate()
    with pytest.raises(ValueError, match="operator_constant_M"):
        gate.radius_for_M(float("nan"))
    with pytest.raises(ValueError, match="target_radius"):
        gate.M_budget_for_radius(-1.0)


def test_truth_boundary_keeps_source_and_pa16_gates_closed() -> None:
    gate = KokunoPA10FixedPointRadiusGate(operator_constant_M=1.0)
    truth = gate.report()["truth_boundary"]
    assert truth["source_fixed_point_radius_formula_executable"] is True
    assert truth["selected_lambda_is_source_existential_threshold"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_contraction_invariant_ball_machine_verified"] is False
    assert truth["source_contraction_factor_machine_verified"] is False
    assert truth["source_fixed_point_distance_machine_bound"] is False
    assert truth["source_fixed_point_solved"] is False
    assert truth["source_B0_dependencies_machine_bound"] is False
    assert truth["source_T_sh_lower_bound_verified"] is False
    assert truth["selected_pa16_handoff_allowed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_payload_roundtrip_and_truth_tamper_fail_closed(tmp_path) -> None:
    gate = KokunoPA10FixedPointRadiusGate(operator_constant_M=2.5)
    payload = gate.to_payload()
    replay = KokunoPA10FixedPointRadiusGate.from_payload(payload)
    assert replay.sha256 == gate.sha256
    assert replay.to_payload() == payload

    path = tmp_path / "fixed_point_radius.json"
    gate.save_json(path)
    loaded = KokunoPA10FixedPointRadiusGate.load_json(path)
    assert loaded.sha256 == gate.sha256
    assert json.loads(path.read_text())["sha256"] == gate.sha256

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["source_fixed_point_distance_machine_bound"] = True
    # Recompute the outer digest so the fail-close check reaches truth metadata.
    body = copy.deepcopy(tampered)
    body.pop("sha256")
    import hashlib

    canonical = json.dumps(body, sort_keys=True, separators=(",", ":"), allow_nan=False)
    tampered["sha256"] = hashlib.sha256(canonical.encode()).hexdigest()
    with pytest.raises(ValueError, match="truth boundary"):
        KokunoPA10FixedPointRadiusGate.from_payload(tampered)
