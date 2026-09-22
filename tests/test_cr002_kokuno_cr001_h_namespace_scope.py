from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_cr001_h_namespace_scope import (
    CONTRACT_PATH,
    audit_contract,
    parameter_namespace_key,
)


ROOT = Path(__file__).resolve().parents[1]


def _contract():
    return json.loads((ROOT / CONTRACT_PATH).read_text())


def test_current_cross_lineage_h_namespace_contract_passes():
    report = audit_contract(ROOT)
    assert report["canonical_cr001_h"] == pytest.approx(0.005, rel=0.0, abs=0.0)
    assert report["kokuno_a1_h_current"] == pytest.approx(0.005, rel=0.0, abs=2e-15)
    assert report["same_numeric_value_observed"] is True
    assert report["shared_parameter_identity"] is False
    assert report["canonical_eq45_velocity_export_ready"] is True
    assert report["kokuno_global_velocity_export_ready"] is False
    assert report["pde_validated"] is False


def test_same_scalar_is_not_same_semantic_parameter_namespace():
    payload = _contract()
    left = payload["parameter_namespaces"]["canonical_cr001_h"]
    right = payload["parameter_namespaces"]["kokuno_a1_h_current"]
    assert left["current_value"] == right["current_value"] == 0.005
    assert parameter_namespace_key(left) != parameter_namespace_key(right)


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("identity_firewall", "shared_parameter_identity"), True),
        (("identity_firewall", "shared_source_provenance"), True),
        (("identity_firewall", "automatic_value_synchronization"), True),
        (("identity_firewall", "canonical_cr001_evidence_transfers_to_kokuno"), True),
        (("identity_firewall", "kokuno_evidence_transfers_to_canonical_cr001"), True),
        (
            ("identity_firewall", "numeric_equality_can_promote_velocity_visual_pde_or_exactness_state"),
            True,
        ),
        (("parameter_namespaces", "kokuno_a1_h_current", "classification"), "public_source_fact"),
        (("parameter_namespaces", "kokuno_a1_h_current", "source_exact_numeric_value"), True),
        (("parameter_namespaces", "canonical_cr001_h", "source_exact_numeric_value"), True),
        (("independent_state_boundary", "kokuno_a1_lminus1_hold", "pde_validated"), True),
        (
            ("independent_state_boundary", "kokuno_a1_lminus1_hold",
             "unified_global_cartesian_velocity_export_ready"),
            True,
        ),
        (("canonical_cr001_freeze", "pde_residual_max"), 0.002),
        (("canonical_cr001_freeze", "divergence_max"), 2e-5),
    ],
)
def test_mutations_fail_closed(path, value):
    payload = copy.deepcopy(_contract())
    node = payload
    for key in path[:-1]:
        node = node[key]
    node[path[-1]] = value
    with pytest.raises(ValueError):
        audit_contract(ROOT, contract_override=payload)


def test_future_policy_does_not_require_cross_lineage_synchronization():
    payload = _contract()
    policy = payload["future_change_policy"]
    assert "does not force a Kokuno h_current change" in policy["canonical_h_change"]
    assert "does not force a canonical CR001 h change" in policy["kokuno_h_change"]
    assert "insufficient for semantic identity" in policy["same_value_after_change"]
