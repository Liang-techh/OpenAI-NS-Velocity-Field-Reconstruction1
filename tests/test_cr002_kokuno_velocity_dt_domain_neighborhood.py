from __future__ import annotations

import copy
import json

import numpy as np
import pytest

import openai_ns_reconstruction.audit_kokuno_velocity_dt_domain_neighborhood as audit


def _scope() -> dict:
    return audit.load_scope()


def test_scope_contract_and_four_way_classification_are_fail_closed() -> None:
    payload = _scope()
    audit.validate_scope(payload)

    assert payload["schema"] == audit.SCHEMA
    assert payload["task"] == audit.TASK
    assert set(payload["classification"]) == {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }

    domain = payload["domain_scope"]
    assert domain["inner_source_X_only"] is True
    assert domain["outside_inner_X_raises"] is True
    assert domain["analytic_velocity_dt_is_fixed_cartesian"] is True
    assert (
        domain[
            "pointwise_velocity_dt_evaluable_does_not_imply_frozen_fd_stencil_fits"
        ]
        is True
    )
    assert domain["agent4_fd4_evidence_is_interior_subdomain_evidence"] is True
    assert domain["agent4_fd4_evidence_establishes_entire_inner_X_domain"] is False
    assert domain["agent4_fd4_evidence_establishes_inner_boundary_behavior"] is False
    assert domain["global_pde_derivative_validated"] is False

    states = payload["independent_states"]
    assert states["inner_velocity_dt_interface_available"] is True
    assert states["inner_velocity_dt_agent4_interior_audit_registered"] is True
    assert states["velocity_export_ready_promoted_by_this_contract"] is False
    assert states["visual_correspondence_verified"] is False
    assert states["pde_validated"] is False
    assert states["paper_exact"] is False
    assert states["openai_field_identified"] is False


def test_pointwise_velocity_dt_can_be_valid_while_fine_fd4_leaves_inner_X() -> None:
    witness = audit.boundary_witness()

    assert witness["pointwise_velocity_dt_finite"] is True
    assert np.isfinite(witness["pointwise_velocity_dt_norm"])
    assert witness["fd4_step"] == 1.0e-4
    assert witness["positive_outer_offset"] == 2.0e-4
    assert np.isclose(
        witness["crossing_delta_at_eta_zero"], 1.0e-4, rtol=0.0, atol=5.0e-16
    )
    assert witness["shifted_source_X_fraction"] > 1.0
    assert witness["velocity_at_t_plus_2h_rejected_by_inner_domain"] is True
    assert "not a physical singularity" in witness["interpretation"]


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        (
            "domain_scope",
            "agent4_fd4_evidence_establishes_entire_inner_X_domain",
            True,
        ),
        ("domain_scope", "full_inner_domain_pde_derivative_validated", True),
        ("domain_scope", "global_pde_derivative_validated", True),
        ("independent_states", "velocity_export_ready_promoted_by_this_contract", True),
        ("independent_states", "visual_correspondence_verified", True),
        ("independent_states", "pde_validated", True),
        ("independent_states", "paper_exact", True),
        ("independent_states", "openai_field_identified", True),
    ],
)
def test_rejects_scope_or_truth_state_promotion(
    section: str, key: str, value: bool
) -> None:
    payload = copy.deepcopy(_scope())
    payload[section][key] = value
    with pytest.raises(ValueError):
        audit.validate_scope(payload)


def test_rejects_provenance_laundering() -> None:
    payload = copy.deepcopy(_scope())
    payload["classification"]["public_source_fact"].append(
        "Agent-4 FD4 sampling locations and finite-difference steps are paper parameters."
    )
    payload["classification"]["autonomous_design"] = []
    with pytest.raises(ValueError):
        audit.validate_scope(payload)


def test_rejects_cr001_threshold_relaxation(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    constraints = json.loads(audit.CR001_PATH.read_text(encoding="utf-8"))
    constraints["validation"]["thresholds"]["pde_residual_max"] = 2.0e-3
    mutated = tmp_path / "constraints.json"
    mutated.write_text(json.dumps(constraints), encoding="utf-8")
    monkeypatch.setattr(audit, "CR001_PATH", mutated)

    with pytest.raises(ValueError, match="threshold drift"):
        audit.validate_scope(_scope())


def test_rejects_cr001_free_force_or_amplitude_collapse(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    constraints = json.loads(audit.CR001_PATH.read_text(encoding="utf-8"))
    constraints["forcing"]["restriction"] = "pointwise free force allowed"
    mutated = tmp_path / "constraints.json"
    mutated.write_text(json.dumps(constraints), encoding="utf-8")
    monkeypatch.setattr(audit, "CR001_PATH", mutated)

    with pytest.raises(ValueError, match="free-force prohibition drift"):
        audit.validate_scope(_scope())


def test_rejects_agent4_protocol_drift(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(audit.agent4, "TIME_STEPS", (4.0e-4, 2.0e-4, 2.0e-4))
    with pytest.raises(ValueError, match="FD4 step drift"):
        audit.validate_scope(_scope())


def test_audit_report_keeps_scientific_states_independent() -> None:
    report = audit.audit()
    assert report["scope_validated"] is True
    assert report["boundary_witness"][
        "velocity_at_t_plus_2h_rejected_by_inner_domain"
    ] is True
    truth = report["truth_boundary"]
    assert truth["inner_velocity_dt_interface_available"] is True
    assert truth["inner_velocity_dt_agent4_interior_audit_registered"] is True
    assert truth["velocity_export_ready_promoted_by_this_contract"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
