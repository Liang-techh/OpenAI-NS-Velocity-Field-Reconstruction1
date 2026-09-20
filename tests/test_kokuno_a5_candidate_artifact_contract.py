from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_candidate_artifact_contract import (
    FINAL_PROJECT_GATES,
    PIPELINE_ORDER,
    PUBLIC_API_CONTRACT,
    deterministic_contract,
    load_contract,
    validate_contract,
    write_contract,
)


def test_contract_is_fail_closed_and_round_trips(tmp_path):
    exact_head = "a" * 40
    payload = deterministic_contract(exact_head=exact_head)
    validate_contract(payload)

    assert payload["pipeline_order"] == PIPELINE_ORDER
    assert payload["public_api_contract"] == PUBLIC_API_CONTRACT
    assert payload["exact_head"] == exact_head
    assert payload["stage_state"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert payload["candidate_payload"]["materialized"] is False
    assert payload["candidate_payload"]["evaluator_bindings"] is None
    assert payload["final_project_gates_unchanged"] == FINAL_PROJECT_GATES
    assert payload["truth_boundary"]["contract_only_not_candidate"] is True
    assert payload["truth_boundary"]["pde_validated"] is False

    path = tmp_path / "contract.json"
    written = write_contract(path, exact_head=exact_head)
    loaded = load_contract(path)
    assert loaded == written


def test_contract_rejects_candidate_or_scientific_promotion():
    base = deterministic_contract()

    for mutate in (
        lambda p: p["candidate_payload"].update(
            {"materialized": True, "candidate_id": "surrogate"}
        ),
        lambda p: p["truth_boundary"].update({"pde_validated": True}),
        lambda p: p["truth_boundary"].update({"velocity_export_ready": True}),
        lambda p: p["baseline_vs_kokuno"].update(
            {"kokuno_same_protocol_full_candidate_residual_available": True}
        ),
        lambda p: p["baseline_vs_kokuno"].update({"comparison_performed": True}),
        lambda p: p["final_project_gates_unchanged"].update(
            {"normalized_momentum_max": 2e-3}
        ),
    ):
        payload = copy.deepcopy(base)
        mutate(payload)
        with pytest.raises(RuntimeError):
            validate_contract(payload)


def test_contract_rejects_forcing_escape_hatch_and_digest_tamper():
    base = deterministic_contract()

    payload = copy.deepcopy(base)
    payload["public_api_contract"]["forcing"]["may_not_be_defined_from_residual"] = False
    with pytest.raises(RuntimeError):
        validate_contract(payload)

    payload = copy.deepcopy(base)
    payload["contract_sha256"] = "0" * 64
    with pytest.raises(RuntimeError):
        validate_contract(payload)


def test_contract_rejects_pipeline_reordering():
    payload = deterministic_contract()
    payload["pipeline_order"] = list(reversed(PIPELINE_ORDER))
    with pytest.raises(RuntimeError):
        validate_contract(payload)
