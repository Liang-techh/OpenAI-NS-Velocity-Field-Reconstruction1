from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.audit_kokuno_axis_completed_first_cell_quadrature_scope import (
    CONTRACT_PATH,
    _AUTONOMOUS_MARKER,
    audit_contract,
)


def _contract() -> dict:
    return json.loads(CONTRACT_PATH.read_text())


def test_axis_completed_first_cell_scope_audit_passes() -> None:
    report = audit_contract()
    assert report["formal_axis_based_inverse_verified"] is False
    assert report["formal_first_cell_accuracy_verified"] is False
    assert report["pde_validated"] is False
    assert report["conclusion"] == (
        "axis_completed_trapezoid_is_scoped_discrete_realization_not_formal_first_cell_accuracy_certificate"
    )


def test_constant_source_witness_distinguishes_e2_first_cell() -> None:
    report = audit_contract()
    assert report["e1"]["relative_primitive_overestimate"] == pytest.approx(0.0, abs=2e-15)
    assert report["e1"]["formal_stress_at_r_min"] == pytest.approx(-0.1)
    assert report["e1"]["trapezoid_stress_at_r_min"] == pytest.approx(-0.1)

    assert report["e2"]["formal_primitive"] == pytest.approx(0.2**3 / 3.0)
    assert report["e2"]["trapezoid_primitive"] == pytest.approx(0.5 * 0.2**3)
    assert report["e2"]["relative_primitive_overestimate"] == pytest.approx(0.5)
    assert report["e2"]["formal_stress_at_r_min"] == pytest.approx(-0.2 / 3.0)
    assert report["e2"]["trapezoid_stress_at_r_min"] == pytest.approx(-0.1)


def test_rejects_formal_first_cell_accuracy_promotion() -> None:
    data = _contract()
    data["representation_scope"]["formal_first_cell_accuracy_verified"] = True
    with pytest.raises(ValueError, match="unsupported promotion"):
        audit_contract(data)


def test_rejects_axis_inverse_promotion() -> None:
    data = _contract()
    data["representation_scope"]["formal_axis_based_inverse_verified"] = True
    with pytest.raises(ValueError, match="unsupported promotion"):
        audit_contract(data)


def test_rejects_e2_quadrature_witness_laundering() -> None:
    data = _contract()
    data["mechanics_witness"]["e2"]["relative_primitive_overestimate"] = 0.0
    with pytest.raises(ValueError, match="e2 mechanics witness drift"):
        audit_contract(data)


def test_rejects_mechanics_witness_as_candidate_evidence() -> None:
    data = _contract()
    data["mechanics_witness"]["candidate_evidence"] = True
    with pytest.raises(ValueError, match="candidate evidence"):
        audit_contract(data)


def test_rejects_autonomous_trapezoid_as_public_source_fact() -> None:
    data = _contract()
    data["provenance"]["autonomous_design"].remove(_AUTONOMOUS_MARKER)
    data["provenance"]["public_source_fact"].append(_AUTONOMOUS_MARKER)
    with pytest.raises(ValueError, match="autonomous design|public-source fact"):
        audit_contract(data)


def test_rejects_cr001_momentum_threshold_relaxation() -> None:
    data = _contract()
    data["cr001_snapshot"]["pde_residual_max"] = 0.002
    with pytest.raises(ValueError, match="threshold drift"):
        audit_contract(data)


def test_rejects_residual_defined_free_force_permission() -> None:
    data = _contract()
    data["cr001_snapshot"]["residual_defined_pointwise_free_force_allowed"] = True
    with pytest.raises(ValueError, match="free residual-defined forcing"):
        audit_contract(data)


def test_rejects_amplitude_collapse_permission() -> None:
    data = _contract()
    data["cr001_snapshot"]["amplitude_collapse_allowed"] = True
    with pytest.raises(ValueError, match="amplitude collapse"):
        audit_contract(data)


def test_rejects_pde_promotion() -> None:
    data = _contract()
    data["truth_state"]["pde_validated"] = True
    with pytest.raises(ValueError, match="unsupported truth-state promotion"):
        audit_contract(data)


def test_rejects_pending_pde_as_callable_delivery_blocker() -> None:
    data = _contract()
    data["truth_state"]["pde_pending_blocks_existing_callable_velocity_delivery"] = True
    with pytest.raises(ValueError, match="unsupported truth-state promotion"):
        audit_contract(data)
