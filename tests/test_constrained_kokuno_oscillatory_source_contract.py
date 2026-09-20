import hashlib
import inspect
import json
from dataclasses import fields, replace

import pytest

from openai_ns_reconstruction.kokuno_oscillatory_source_contract import (
    PARENT_EXACT_HEAD,
    SOURCE_DATE,
    SOURCE_PATH,
    SOURCE_VERSION,
    KokunoOscillatorySourceContract,
    default_kokuno_oscillatory_source_contract,
)


def _digest(payload):
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def test_corrected_source_and_parent_are_pinned():
    contract = default_kokuno_oscillatory_source_contract()
    assert contract.source_date == SOURCE_DATE == "2026-09-09"
    assert contract.source_path == SOURCE_PATH == "KokunoYumeto updated/output.pdf"
    assert contract.source_version == SOURCE_VERSION == "corrected-208-page-reconstruction"
    assert contract.parent_exact_head == PARENT_EXACT_HEAD
    assert contract.parent_exact_head == "946b4d13f63ae6761f02e9e66d4c3c7ed485ffa6"


def test_complete_curl_requires_cutoff_gradient_term():
    contract = default_kokuno_oscillatory_source_contract()
    assert contract.complete_curl_terms == (
        "chi_k curl(B^(k))",
        "grad(chi_k) x B^(k)",
    )
    assert contract.cutoff_gradient_term_required is True

    with pytest.raises(ValueError, match="complete curl"):
        replace(contract, cutoff_gradient_term_required=False).validate()

    with pytest.raises(ValueError, match="cutoff-gradient"):
        replace(contract, complete_curl_terms=("chi_k curl(B^(k))",)).validate()


def test_source_and_repository_realization_are_explicitly_partitioned():
    contract = default_kokuno_oscillatory_source_contract()
    payload = contract.semantic_payload()
    source = payload["source_bound"]
    realization = payload["repository_realization"]

    assert "source_phase" in source
    assert "source_structures" in source
    assert "localized_base_field_components" in source
    assert "repository_public_z_pullback" not in source

    assert realization["repository_public_z_pullback"] == "Z_beta = eps_beta z"
    assert realization["paper_exact"] is False
    assert realization["pde_validated"] is False
    assert "not recovered hidden source data" in realization["repository_parameter_status"]
    assert "numerical approximations" in realization["repository_derivative_status"]


def test_contract_has_no_residual_pressure_forcing_or_tuning_surface():
    field_names = {field.name for field in fields(KokunoOscillatorySourceContract)}
    prohibited = {
        "residual",
        "target",
        "pressure",
        "forcing",
        "viscosity",
        "nu",
        "gain",
        "amplitude",
        "phase_offset",
        "spatial_step",
        "scientific_threshold",
    }
    assert field_names.isdisjoint(prohibited)

    signature_names = set(inspect.signature(KokunoOscillatorySourceContract).parameters)
    assert signature_names.isdisjoint(prohibited)


def test_manifest_roundtrip_and_identity_are_deterministic(tmp_path):
    first = default_kokuno_oscillatory_source_contract()
    second = default_kokuno_oscillatory_source_contract()
    assert first.contract_sha256 == second.contract_sha256

    path = first.save_manifest(tmp_path / "source_contract.json")
    loaded = KokunoOscillatorySourceContract.load_manifest(path)
    assert loaded == first
    assert loaded.contract_sha256 == first.contract_sha256


def test_manifest_rejects_checksum_tamper(tmp_path):
    contract = default_kokuno_oscillatory_source_contract()
    path = contract.save_manifest(tmp_path / "source_contract.json")
    raw = json.loads(path.read_text())
    raw["payload"]["source_bound"]["source_phase"] = "tampered"
    path.write_text(json.dumps(raw))

    with pytest.raises(ValueError, match="checksum mismatch"):
        KokunoOscillatorySourceContract.load_manifest(path)


def test_manifest_rejects_semantic_tamper_even_with_recomputed_checksum(tmp_path):
    contract = default_kokuno_oscillatory_source_contract()
    path = contract.save_manifest(tmp_path / "source_contract.json")
    raw = json.loads(path.read_text())
    raw["payload"]["source_bound"]["cutoff_gradient_term_required"] = False
    raw["contract_sha256"] = _digest(raw["payload"])
    path.write_text(json.dumps(raw))

    with pytest.raises(ValueError, match="complete curl"):
        KokunoOscillatorySourceContract.load_manifest(path)


def test_truth_boundary_cannot_be_promoted_by_contract():
    contract = default_kokuno_oscillatory_source_contract()
    with pytest.raises(ValueError, match="paper-exact"):
        replace(contract, paper_exact=True).validate()
    with pytest.raises(ValueError, match="PDE validation"):
        replace(contract, pde_validated=True).validate()
