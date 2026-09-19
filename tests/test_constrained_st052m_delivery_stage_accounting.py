from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_st052m_delivery_stage_accounting import (
    CONTRACT_PATH,
    ROOT,
    GovernanceError,
    audit,
)


def _contract() -> dict:
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def _write_contract(tmp_path: Path, payload: dict) -> Path:
    path = tmp_path / "contract.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _reject(tmp_path: Path, mutate) -> None:
    payload = deepcopy(_contract())
    mutate(payload)
    with pytest.raises(GovernanceError):
        audit(ROOT, contract_path=_write_contract(tmp_path, payload))


def test_current_stage_accounting_passes() -> None:
    result = audit()
    assert result["status"] == "pass"
    assert result["render_receipt_available"] is True
    assert result["live_materialization_complete"] is False
    assert result["experimental_child_velocity_export_ready"] is False
    assert result["delivery_route_complete"] is False
    assert result["canonical_candidate_unchanged"] is True
    assert result["cr001_unchanged"] is True


def test_rejects_render_as_live_materialization(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["render_evidence"].__setitem__("counts_as_live_ancestry_materialization", True),
    )


def test_rejects_render_as_post_materialization_stage(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["delivery_stage_accounting"].__setitem__(
            "render_before_materialization_satisfies_post_materialization_render_stage", True
        ),
    )


def test_rejects_experimental_child_export_promotion(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["exact_source_child"].__setitem__("velocity_export_ready", True),
    )


def test_rejects_automatic_render_evidence_transfer(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["later_materialized_child_evidence_transfer"].__setitem__(
            "automatic_transfer_allowed", True
        ),
    )


def test_rejects_parameter_only_identity_transfer(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["later_materialized_child_evidence_transfer"].__setitem__(
            "same_parameters_alone_are_sufficient", True
        ),
    )


def test_rejects_visual_correspondence_promotion(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["render_evidence"].__setitem__("counts_as_visual_correspondence_verified", True),
    )


def test_rejects_invalid_source_class(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["source_classification"].__setitem__("delivery_stage_accounting_policy", "internal_fact"),
    )


def test_rejects_cr001_threshold_drift(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["canonical_cr001"].__setitem__("pde_residual_max", 0.01),
    )


def test_rejects_delivery_route_reordering(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["snapshot"].__setitem__(
            "project_status_delivery_shortest_route", "render_then_materialize_exact_st052m_child"
        ),
    )
