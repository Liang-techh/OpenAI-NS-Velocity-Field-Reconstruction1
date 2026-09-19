from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_st052m_linear_temporal_governance import (
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


def test_current_linear_temporal_governance_passes() -> None:
    result = audit()
    assert result["status"] == "pass"
    assert result["activation_law_bound_to_identity"] is True
    assert result["endpoint_identity_scoped"] is True
    assert result["fresh_time_derivative_required_for_pde"] is True
    assert result["static_evidence_transfer_fail_closed"] is True
    assert result["temporal_child_velocity_export_ready"] is False
    assert result["canonical_candidate_unchanged"] is True
    assert result["cr001_unchanged"] is True


def test_rejects_endpoints_as_interior_identity(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["endpoint_identities"].__setitem__("endpoint_identity_implies_interior_identity", True),
    )


def test_rejects_endpoint_velocity_interpolation_claim(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["temporal_representation"].__setitem__(
            "velocity_is_linear_interpolation_of_endpoint_velocities", True
        ),
    )


def test_rejects_activation_law_drift(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["temporal_representation"].__setitem__("activation_formula", "g(t)=4*(t-0.25)"),
    )


def test_rejects_erasing_chain_rule_term(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["time_derivative_semantics"].__setitem__(
            "extra_transform_time_derivative_terms_present_in_general", False
        ),
    )


def test_rejects_static_momentum_receipt_transfer(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["time_derivative_semantics"].__setitem__(
            "parent_or_static_momentum_residual_receipt_transfers", True
        ),
    )


def test_rejects_reference_energy_as_all_time_acceptance(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["energy_semantics"].__setitem__(
            "reference_time_identity_implies_all_validation_time_energy_acceptance", True
        ),
    )


def test_rejects_order64_as_cr001_energy_ladder(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["energy_semantics"].__setitem__(
            "order64_descriptive_energy_replaces_cr001_24_48_96_ladder", True
        ),
    )


def test_rejects_sampled_divergence_as_cr001_acceptance(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["pr587_evidence_scope"].__setitem__(
            "sampled_divergence_is_cr001_heldout_acceptance", True
        ),
    )


def test_rejects_static_path_receipt_transfer(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["static_evidence_transfer"].__setitem__(
            "static_pr559_full_interval_material_paths_transfer_to_temporal_child", True
        ),
    )


def test_rejects_static_render_wholesale_transfer(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["static_evidence_transfer"].__setitem__(
            "static_pr583_three_time_render_receipt_transfers_wholesale", True
        ),
    )


def test_rejects_temporal_child_export_promotion(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["temporal_representation"].__setitem__("velocity_export_ready", True),
    )


def test_rejects_endpoints_alone_as_candidate_identity(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["delivery_identity_requirements"].__setitem__(
            "same_endpoints_alone_are_sufficient_identity", True
        ),
    )


def test_rejects_pde_failure_blocking_velocity_delivery(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["delivery_identity_requirements"].__setitem__(
            "pde_failure_or_pending_blocks_callable_velocity_delivery", True
        ),
    )


def test_rejects_visual_correspondence_promotion(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["pr587_evidence_scope"].__setitem__(
            "counts_as_visual_correspondence_verified", True
        ),
    )


def test_rejects_invalid_source_class(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["source_classification"].__setitem__(
            "linear_activation_law_g_of_t", "public_source_fact"
        ),
    )


def test_rejects_cr001_threshold_drift(tmp_path: Path) -> None:
    _reject(
        tmp_path,
        lambda c: c["canonical_cr001"].__setitem__("pde_residual_max", 0.01),
    )
