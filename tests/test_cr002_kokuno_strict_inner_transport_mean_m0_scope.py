from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_strict_inner_transport_mean_m0_scope import (
    AuditError,
    _manufactured_witness,
    audit_repository,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "configs" / "kokuno_strict_inner_transport_mean_m0_scope.json"


def _contract() -> dict[str, object]:
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def test_cr002_transport_m0_scope_audit_passes() -> None:
    report = audit_repository(ROOT)
    assert report["status"] == "PASS"
    assert report["target_pr"] == 850
    assert report["pde_validated"] is False
    assert report["callable_velocity_delivery_blocked"] is False


def test_manufactured_m2_mode_is_in_m0_nullspace_but_not_zero_full_ring() -> None:
    witness = _manufactured_witness(_contract())
    assert abs(witness["computed_m0_mean_z"]) <= 2e-15
    assert witness["computed_full_ring_rms"] > 0.4
    assert witness["computed_mean_to_full_ratio"] <= 5e-15


@pytest.mark.parametrize(
    ("mutator", "match"),
    [
        (
            lambda c: c["projection_boundary"].__setitem__(
                "zero_m0_mean_implies_zero_full_ring_transport", True
            ),
            "zero m0 mean",
        ),
        (
            lambda c: c["projection_boundary"].__setitem__(
                "small_m0_mean_implies_small_full_ring_transport", True
            ),
            "small m0 mean",
        ),
        (
            lambda c: c["norm_scope"].__setitem__(
                "pr850_full_ring_rms_is_whole_domain_volume_weighted_L2", True
            ),
            "norm-scope promotion",
        ),
        (
            lambda c: c["truth_boundary"].__setitem__("pde_validated", True),
            "truth promotion forbidden",
        ),
        (
            lambda c: c["truth_boundary"].__setitem__(
                "pde_pending_blocks_callable_velocity_delivery", True
            ),
            "truth promotion forbidden",
        ),
        (
            lambda c: c["cr001_snapshot"].__setitem__("pde_residual_L2", 0.002),
            "pde_residual_L2 drifted",
        ),
        (
            lambda c: c["cr001_snapshot"].__setitem__(
                "residual_defined_pointwise_free_force_allowed", True
            ),
            "residual_defined_pointwise_free_force_allowed must remain false",
        ),
        (
            lambda c: c["cr001_snapshot"].__setitem__("amplitude_collapse_allowed", True),
            "amplitude_collapse_allowed must remain false",
        ),
        (
            lambda c: c["forbidden_promotions"].__setitem__(
                "small_m0_mean_to_small_full_transport", True
            ),
            "forbidden promotion enabled",
        ),
        (
            lambda c: c["target"].__setitem__("head_sha", "0" * 40),
            "target head_sha drifted",
        ),
    ],
)
def test_fail_closed_mutations(mutator, match: str) -> None:
    mutated = copy.deepcopy(_contract())
    mutator(mutated)
    with pytest.raises(AuditError, match=match):
        audit_repository(ROOT, contract_override=mutated)


def test_autonomous_angular_protocol_cannot_be_laundered_into_public_source() -> None:
    mutated = copy.deepcopy(_contract())
    entry = mutated["classification"]["autonomous_design"].pop(1)
    mutated["classification"]["public_source_fact"].append(entry)
    with pytest.raises(AuditError, match="angular ladder must remain autonomous_design"):
        audit_repository(ROOT, contract_override=mutated)


def test_pr853_boundary_stays_distinct_from_m0_nullspace_scope() -> None:
    contract = _contract()
    relation = contract["relationship_to_sibling_governance"]
    assert relation["duplicate"] is False
    assert "nonlinear" in relation["pr853"].lower()
    assert "nullspace" in relation["this_increment"].lower()
