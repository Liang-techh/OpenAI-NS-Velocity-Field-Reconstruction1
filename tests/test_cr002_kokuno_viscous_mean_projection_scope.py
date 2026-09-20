from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_viscous_mean_projection_scope import (
    CONTRACT_PATH,
    UPSTREAM_PATH,
    _REQUIRED_UPSTREAM_CONSTANTS,
    _validate_contract,
    audit,
    manufactured_vector_laplacian_witness,
)


ROOT = Path(__file__).resolve().parents[1]


def _payloads():
    contract = json.loads((ROOT / CONTRACT_PATH).read_text(encoding="utf-8"))
    cr001 = json.loads((ROOT / "configs/constraints.json").read_text(encoding="utf-8"))
    upstream = dict(_REQUIRED_UPSTREAM_CONSTANTS)
    source = (ROOT / UPSTREAM_PATH).read_text(encoding="utf-8")
    return contract, cr001, upstream, source


def test_live_scope_audit_passes_without_scientific_promotion():
    report = audit()
    assert report["status"] == "pass"
    assert report["pde_validated"] is False
    assert report["scope"] == (
        "cartesian_vector_viscous_mean_projection_not_naive_scalar_laplacian_of_cylindrical_mean"
    )
    witness = report["manufactured_witness"]
    assert abs(witness["vector_radial_laplacian"]) <= 1.0e-15
    assert witness["naive_scalar_radial_laplacian"] > 1.0


def test_manufactured_witness_requires_geometric_basis_term():
    witness = manufactured_vector_laplacian_witness(radius=0.7)
    assert witness["naive_scalar_radial_laplacian"] == pytest.approx(1.0 / 0.7)
    assert witness["geometric_basis_correction"] == pytest.approx(-1.0 / 0.7)
    assert witness["vector_radial_laplacian"] == pytest.approx(0.0, abs=1.0e-15)
    assert witness["cartesian_vector_laplacian_radial_projection"] == 0.0


def test_rejects_naive_scalar_laplacian_promotion():
    contract, cr001, upstream, source = _payloads()
    contract["upstream_operator_semantics"][
        "naive_componentwise_scalar_laplacian_of_mean_components"
    ] = True
    with pytest.raises(AssertionError, match="operator-scope laundering"):
        _validate_contract(contract, cr001, upstream, source)


def test_rejects_operator_on_mean_equivalence_promotion():
    contract, cr001, upstream, source = _payloads()
    contract["upstream_operator_semantics"]["operator_on_mean_equivalence_established_here"] = True
    with pytest.raises(AssertionError, match="operator-scope laundering"):
        _validate_contract(contract, cr001, upstream, source)


def test_rejects_losing_radial_geometric_term():
    contract, cr001, upstream, source = _payloads()
    contract["cylindrical_vector_laplacian_boundary"]["axisymmetric_m0_radial_component"] = (
        "(d_rr + r^-1 d_r + d_zz) u_r"
    )
    with pytest.raises(AssertionError, match="radial geometric term lost"):
        _validate_contract(contract, cr001, upstream, source)


def test_rejects_truth_state_equivalence_promotion():
    contract, cr001, upstream, source = _payloads()
    contract["states"]["operator_on_mean_equivalence_established"] = True
    with pytest.raises(AssertionError, match="truth-state promotion"):
        _validate_contract(contract, cr001, upstream, source)


def test_rejects_preclaiming_queued_upstream_ci():
    contract, cr001, upstream, source = _payloads()
    contract["states"]["upstream_841_ci_passed"] = True
    with pytest.raises(AssertionError, match="pre-promoted"):
        _validate_contract(contract, cr001, upstream, source)


def test_rejects_momentum_threshold_relaxation():
    contract, cr001, upstream, source = _payloads()
    mutated = copy.deepcopy(cr001)
    mutated["validation"]["thresholds"]["pde_residual_max"] = 2.0e-3
    with pytest.raises(AssertionError, match="momentum max threshold drift"):
        _validate_contract(contract, mutated, upstream, source)


def test_rejects_free_residual_defined_force_permission():
    contract, cr001, upstream, source = _payloads()
    contract["canonical_cr001"]["free_residual_defined_force_allowed"] = True
    with pytest.raises(AssertionError, match="free residual force enabled"):
        _validate_contract(contract, cr001, upstream, source)


def test_rejects_amplitude_collapse_permission():
    contract, cr001, upstream, source = _payloads()
    contract["canonical_cr001"]["reject_amplitude_collapse"] = False
    with pytest.raises(AssertionError, match="amplitude collapse allowed"):
        _validate_contract(contract, cr001, upstream, source)


def test_rejects_four_way_provenance_laundering():
    contract, cr001, upstream, source = _payloads()
    contract["classification"]["public_source_fact"].extend(
        contract["classification"].pop("autonomous_design")
    )
    with pytest.raises(AssertionError, match="four-way provenance"):
        _validate_contract(contract, cr001, upstream, source)


def test_rejects_delivery_blockage_from_scoped_upstream_failure():
    contract, cr001, upstream, source = _payloads()
    contract["promotion_boundary"][
        "upstream_failure_invalidates_existing_callable_velocity_delivery"
    ] = True
    with pytest.raises(AssertionError, match="incorrectly blocks"):
        _validate_contract(contract, cr001, upstream, source)


def test_rejects_upstream_viscosity_drift():
    contract, cr001, upstream, source = _payloads()
    mutated = dict(upstream)
    mutated["REPOSITORY_VISCOSITY"] = 0.02
    with pytest.raises(AssertionError, match="upstream constant drift"):
        _validate_contract(contract, cr001, mutated, source)
