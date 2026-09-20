from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_inner_convective_consistency_norm_scope import (
    CONTRACT_PATH,
    _REQUIRED_UPSTREAM_CONSTANTS,
    _validate_contract,
    audit,
)


ROOT = Path(__file__).resolve().parents[1]


def _payloads():
    contract = json.loads((ROOT / CONTRACT_PATH).read_text(encoding="utf-8"))
    cr001 = json.loads((ROOT / "configs/constraints.json").read_text(encoding="utf-8"))
    upstream = dict(_REQUIRED_UPSTREAM_CONSTANTS)
    return contract, cr001, upstream


def test_live_scope_audit_passes_without_scientific_promotion():
    report = audit()
    assert report["status"] == "pass"
    assert report["pde_validated"] is False
    assert report["scope"] == (
        "strict_inner_convective_relative_consistency_not_complete_cr001_momentum_residual"
    )


def test_rejects_relative_rms_laundering_as_cr001_momentum_norm():
    contract, cr001, upstream = _payloads()
    contract["upstream_consistency_metric"]["canonical_cr001_momentum_norm"] = True
    with pytest.raises(AssertionError, match="metric-scope laundering"):
        _validate_contract(contract, cr001, upstream)


def test_rejects_convective_subterm_as_complete_momentum_residual():
    contract, cr001, upstream = _payloads()
    contract["upstream_consistency_metric"]["complete_momentum_residual"] = True
    with pytest.raises(AssertionError, match="metric-scope laundering"):
        _validate_contract(contract, cr001, upstream)


def test_rejects_st006_comparability_promotion():
    contract, cr001, upstream = _payloads()
    contract["states"]["st006_comparable"] = True
    with pytest.raises(AssertionError, match="truth-state promotion"):
        _validate_contract(contract, cr001, upstream)


def test_rejects_pde_promotion_from_subterm_consistency():
    contract, cr001, upstream = _payloads()
    contract["states"]["pde_validated"] = True
    with pytest.raises(AssertionError, match="truth-state promotion"):
        _validate_contract(contract, cr001, upstream)


def test_rejects_preclaiming_queued_upstream_ci():
    contract, cr001, upstream = _payloads()
    contract["states"]["strict_inner_convective_consistency_ci_passed"] = True
    with pytest.raises(AssertionError, match="pre-promoted"):
        _validate_contract(contract, cr001, upstream)


def test_rejects_momentum_threshold_relaxation():
    contract, cr001, upstream = _payloads()
    mutated = copy.deepcopy(cr001)
    mutated["validation"]["thresholds"]["pde_residual_max"] = 2.0e-3
    with pytest.raises(AssertionError, match="momentum max threshold drift"):
        _validate_contract(contract, mutated, upstream)


def test_rejects_free_residual_defined_force_permission():
    contract, cr001, upstream = _payloads()
    contract["canonical_cr001"]["free_residual_defined_force_allowed"] = True
    with pytest.raises(AssertionError, match="free residual force enabled"):
        _validate_contract(contract, cr001, upstream)


def test_rejects_amplitude_collapse_permission():
    contract, cr001, upstream = _payloads()
    contract["canonical_cr001"]["reject_amplitude_collapse"] = False
    with pytest.raises(AssertionError, match="amplitude collapse allowed"):
        _validate_contract(contract, cr001, upstream)


def test_rejects_public_source_laundering_by_losing_four_way_buckets():
    contract, cr001, upstream = _payloads()
    contract["classification"]["public_source_fact"].extend(
        contract["classification"].pop("autonomous_design")
    )
    with pytest.raises(AssertionError, match="four-way provenance"):
        _validate_contract(contract, cr001, upstream)


def test_rejects_delivery_blockage_from_scoped_audit_failure():
    contract, cr001, upstream = _payloads()
    contract["promotion_boundary"][
        "upstream_failure_invalidates_existing_callable_velocity_delivery"
    ] = True
    with pytest.raises(AssertionError, match="incorrectly blocks"):
        _validate_contract(contract, cr001, upstream)


def test_rejects_upstream_relative_gate_drift():
    contract, cr001, upstream = _payloads()
    mutated = dict(upstream)
    mutated["FINE_RELATIVE_RMS_GATE"] = 1.0e-3
    with pytest.raises(AssertionError, match="upstream constant drift"):
        _validate_contract(contract, cr001, mutated)
