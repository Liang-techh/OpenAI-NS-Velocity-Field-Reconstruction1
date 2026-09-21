from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_partial_composite_osc_runtime_identity_scope import (
    SCOPE_REL,
    _repo_root,
    assert_scope,
    audit_scope,
    load_scope,
)


def _scope() -> dict:
    return load_scope()


def _errors(payload: dict) -> set[str]:
    return set(audit_scope(payload))


def test_registered_scope_passes_fail_closed_audit() -> None:
    payload = _scope()
    assert_scope(payload)
    assert payload["machine_locked_state"]["partial_composite_callable_through_current_Xh"] is True
    assert payload["machine_locked_state"]["full_oscillatory_runtime_realization_identity_bound_into_composite_semantic_sha"] is False
    assert payload["machine_locked_state"]["direct_composite_save_load_available"] is False
    assert payload["machine_locked_state"]["composite_semantic_sha_sufficient_for_saved_velocity_identity"] is False


def test_rejects_adapter_identity_laundered_as_full_runtime_identity() -> None:
    payload = _scope()
    payload["observed_upstream"]["partial_composite"][
        "semantic_payload_binds_full_oscillatory_runtime_payload_or_configuration_sha256"
    ] = True
    assert "full_osc_runtime_binding_laundered" in _errors(payload)


def test_rejects_adapter_semantic_laundered_as_full_public_z_payload() -> None:
    payload = _scope()
    payload["observed_upstream"]["oscillatory_batch_adapter"][
        "adapter_semantic_binds_full_public_z_field_to_payload_digest"
    ] = True
    assert "adapter_full_runtime_binding_laundered" in _errors(payload)


def test_rejects_public_z_runtime_digest_promotion_without_binding() -> None:
    payload = _scope()
    payload["observed_upstream"]["public_z_field"][
        "runtime_payload_digest_is_bound_by_pr970_composite_semantic_identity"
    ] = True
    assert "public_z_runtime_digest_laundered" in _errors(payload)


def test_rejects_direct_save_load_promotion() -> None:
    payload = _scope()
    payload["observed_upstream"]["partial_composite"]["composite_direct_save_load_available"] = True
    payload["machine_locked_state"]["direct_composite_save_load_available"] = True
    errors = _errors(payload)
    assert "direct_composite_save_load_laundered" in errors
    assert "premature_promotion:direct_composite_save_load_available" in errors


def test_rejects_semantic_sha_as_saved_velocity_identity() -> None:
    payload = _scope()
    payload["machine_locked_state"]["composite_semantic_sha_sufficient_for_saved_velocity_identity"] = True
    assert "premature_promotion:composite_semantic_sha_sufficient_for_saved_velocity_identity" in _errors(payload)


@pytest.mark.parametrize(
    "key",
    [
        "velocity_export_ready_for_kokuno_route",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ],
)
def test_rejects_scientific_or_delivery_state_promotion(key: str) -> None:
    payload = _scope()
    payload["machine_locked_state"][key] = True
    assert f"premature_promotion:{key}" in _errors(payload)


def test_mechanics_witness_demonstrates_adapter_hash_collision_only_at_metadata_layer() -> None:
    payload = _scope()
    witness = payload["mechanics_only_witness"]
    adapter = witness["adapter_metadata"]
    a = witness["runtime_a"]
    b = witness["runtime_b"]
    canonical = lambda x: json.dumps(x, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    assert hashlib.sha256(canonical(adapter)).digest() == hashlib.sha256(canonical(adapter)).digest()
    assert a["velocity_at_probe"] != b["velocity_at_probe"]
    assert hashlib.sha256(canonical({"adapter": adapter, "runtime": a})).digest() != hashlib.sha256(
        canonical({"adapter": adapter, "runtime": b})
    ).digest()


def test_rejects_mechanics_witness_promoted_to_public_source_fact() -> None:
    payload = _scope()
    payload["mechanics_only_witness"]["classification"] = "public_source_fact"
    payload["mechanics_only_witness"]["not_source_evidence"] = False
    assert "mechanics_witness_source_laundering" in _errors(payload)


def test_requires_all_four_provenance_classes() -> None:
    payload = _scope()
    del payload["provenance_classes"]["pending_unknown"]
    assert "four_way_provenance_missing" in _errors(payload)


def test_rejects_missing_future_runtime_binding_requirement() -> None:
    payload = _scope()
    payload["future_promotion_requirements"] = [
        item for item in payload["future_promotion_requirements"] if "complete concrete oscillatory runtime" not in item
    ]
    assert "future_identity_requirement_missing:complete concrete oscillatory runtime" in _errors(payload)


def test_rejects_cr001_threshold_relaxation() -> None:
    payload = _scope()
    payload["cr001_binding"]["momentum_max"] = 2.0e-3
    assert "cr001_binding_drift" in _errors(payload)


def test_rejects_cr001_free_force_or_collapse_escape() -> None:
    payload = _scope()
    payload["cr001_binding"]["residual_defined_free_forcing_forbidden"] = False
    assert "cr001_binding_drift" in _errors(payload)
    payload = _scope()
    payload["cr001_binding"]["candidate_collapse_forbidden"] = False
    assert "cr001_binding_drift" in _errors(payload)


def test_contract_is_bound_to_canonical_constraints_blob() -> None:
    root = _repo_root()
    constraints_path = root / "configs" / "constraints.json"
    raw = constraints_path.read_bytes()
    blob = hashlib.sha1(f"blob {len(raw)}\0".encode("ascii") + raw).hexdigest()
    payload = _scope()
    assert blob == payload["cr001_binding"]["constraints_blob_sha1"]
    assert payload["cr001_binding"]["derivative_steps"] == [0.02, 0.01, 0.005]
    assert payload["cr001_binding"]["quadrature_orders_per_axis"] == [24, 48, 96]
    assert payload["cr001_binding"]["momentum_max"] == 1.0e-3
    assert payload["cr001_binding"]["divergence_max"] == 1.0e-5


def test_scope_path_is_unique_and_machine_readable() -> None:
    root = _repo_root()
    path = root / SCOPE_REL
    assert path.is_file()
    assert json.loads(path.read_text(encoding="utf-8"))["schema"].startswith("cr002-")
