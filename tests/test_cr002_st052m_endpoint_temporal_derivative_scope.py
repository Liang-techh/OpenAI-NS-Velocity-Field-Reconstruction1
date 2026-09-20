from __future__ import annotations

import copy
import json

import pytest

import openai_ns_reconstruction.audit_st052m_endpoint_temporal_derivative_scope as audit


def _scope() -> dict:
    return audit.load_scope()


def test_contract_and_derivative_witness_are_fail_closed() -> None:
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
    assert payload["classification"]["public_source_fact"] == []

    witness = audit.derivative_witness()
    assert witness == {
        "g2_start": 0.0,
        "g2_mid": 1.0,
        "g2_end": 0.0,
        "g2_prime_start": 8.0,
        "g2_prime_mid": 0.0,
        "g2_prime_end": -8.0,
    }


def test_endpoint_velocity_value_and_velocity_dt_semantics_stay_separate() -> None:
    scope = _scope()["scope"]
    assert scope["g2_zero_endpoints_preserve_additive_velocity_values_at_endpoints"] is True
    assert scope["g2_zero_endpoints_imply_velocity_dt_endpoint_preservation"] is False
    assert scope["g2_zero_endpoints_imply_momentum_residual_endpoint_neutrality"] is False
    assert (
        scope[
            "nonzero_g2_prime_endpoints_require_explicit_velocity_dt_accounting_if_beta_times_C_is_nonzero"
        ]
        is True
    )
    assert scope["held_out_pde_residual_evaluated_by_this_increment"] is False


def test_initial_energy_value_preservation_does_not_promote_pde_or_nontriviality() -> None:
    payload = _scope()
    scope = payload["scope"]
    assert (
        scope[
            "g2_zero_at_t025_preserves_initial_kinetic_energy_value_if_parent_and_domain_are_unchanged"
        ]
        is True
    )
    assert scope["initial_energy_value_preservation_is_not_pde_validation"] is True
    assert scope["initial_energy_value_preservation_is_not_new_nontriviality_evidence"] is True
    assert payload["independent_states"]["pde_validated"] is False


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("scope", "g2_zero_endpoints_imply_velocity_dt_endpoint_preservation", True),
        ("scope", "g2_zero_endpoints_imply_momentum_residual_endpoint_neutrality", True),
        ("scope", "held_out_pde_residual_evaluated_by_this_increment", True),
        ("independent_states", "temporal_child_materialized", True),
        ("independent_states", "temporal_coefficient_selected", True),
        ("independent_states", "visual_correspondence_verified", True),
        ("independent_states", "pde_validated", True),
        ("independent_states", "paper_exact", True),
        ("independent_states", "openai_field_identified", True),
    ],
)
def test_rejects_semantic_or_scientific_promotion(
    section: str, key: str, value: bool
) -> None:
    payload = copy.deepcopy(_scope())
    payload[section][key] = value
    with pytest.raises(ValueError):
        audit.validate_scope(payload)


def test_rejects_autonomous_time_law_provenance_laundering() -> None:
    payload = copy.deepcopy(_scope())
    payload["classification"]["public_source_fact"].append(
        "OpenAI published g2(t)=16*(t-.25)*(.75-t) as its numerical time law."
    )
    with pytest.raises(ValueError):
        audit.validate_scope(payload)


def test_rejects_held_out_pde_tuning_permission() -> None:
    payload = copy.deepcopy(_scope())
    payload["time_basis"][
        "coefficient_may_not_be_selected_from_cr001_held_out_pde_samples"
    ] = False
    with pytest.raises(ValueError):
        audit.validate_scope(payload)


def test_rejects_cr001_threshold_relaxation(tmp_path, monkeypatch: pytest.MonkeyPatch) -> None:
    constraints = json.loads(audit.CR001_PATH.read_text(encoding="utf-8"))
    constraints["validation"]["thresholds"]["pde_residual_max"] = 2.0e-3
    mutated = tmp_path / "constraints.json"
    mutated.write_text(json.dumps(constraints), encoding="utf-8")
    monkeypatch.setattr(audit, "CR001_PATH", mutated)
    with pytest.raises(ValueError):
        audit.validate_scope(_scope())


def test_rejects_cr001_free_force_drift(tmp_path, monkeypatch: pytest.MonkeyPatch) -> None:
    constraints = json.loads(audit.CR001_PATH.read_text(encoding="utf-8"))
    constraints["forcing"]["restriction"] = "pointwise residual-defined free force allowed"
    mutated = tmp_path / "constraints.json"
    mutated.write_text(json.dumps(constraints), encoding="utf-8")
    monkeypatch.setattr(audit, "CR001_PATH", mutated)
    with pytest.raises(ValueError):
        audit.validate_scope(_scope())


def test_rejects_upstream_g2_source_drift(tmp_path, monkeypatch: pytest.MonkeyPatch) -> None:
    source = audit.UPSTREAM_PATH.read_text(encoding="utf-8")
    mutated = tmp_path / "agent7_st052m_endpoint_temporal_curvature_preflight.py"
    mutated.write_text(source.replace("16.0 *", "15.0 *", 1), encoding="utf-8")
    monkeypatch.setattr(audit, "UPSTREAM_PATH", mutated)
    with pytest.raises(ValueError, match="source blob drift"):
        audit.validate_scope(_scope())


def test_audit_keeps_delivery_and_scientific_states_independent() -> None:
    report = audit.audit()
    assert report["scope_validated"] is True
    assert report["derivative_witness"]["g2_prime_start"] == 8.0
    assert report["derivative_witness"]["g2_prime_end"] == -8.0
    truth = report["truth_boundary"]
    assert truth["existing_callable_velocity_delivery_invalidated"] is False
    assert truth["velocity_export_ready_promoted_by_this_contract"] is False
    assert truth["temporal_child_materialized"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
