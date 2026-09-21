from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.kokuno_corrected_oscillation_source_ledger import (
    OSCILLATION_COMPONENT_SHA256,
    PARENT_A2_HEAD,
    PARENT_A2_SOURCE_BLOB,
    SOURCE_BUNDLE_BLOB,
    SOURCE_BUNDLE_SHA256,
    SOURCE_CHECKS_BLOB,
    SOURCE_COMMIT,
    SOURCE_PDF_SHA256,
    SOURCE_READER_PAGES,
    SOURCE_REPOSITORY,
    SOURCE_WORKBENCH_BLOB,
    KokunoCorrectedOscillationSourceLedger,
    _sha256,
    default_kokuno_corrected_oscillation_source_ledger,
)
from openai_ns_reconstruction.kokuno_oscillatory_source_contract import (
    default_kokuno_oscillatory_source_contract,
)


def test_exact_corrected_upstream_identity_and_source_bundle_hashes() -> None:
    ledger = default_kokuno_corrected_oscillation_source_ledger()
    source = ledger.upstream_identity_payload()

    assert source["repository"] == SOURCE_REPOSITORY
    assert source["commit"] == SOURCE_COMMIT
    assert source["reader_pages"] == SOURCE_READER_PAGES == 208
    assert source["workbench_blob"] == SOURCE_WORKBENCH_BLOB
    assert source["checks_blob"] == SOURCE_CHECKS_BLOB
    assert source["source_bundle_blob"] == SOURCE_BUNDLE_BLOB
    assert source["source_bundle_sha256"] == SOURCE_BUNDLE_SHA256
    assert source["pdf_sha256"] == SOURCE_PDF_SHA256
    assert source["oscillation_component_sha256"] == OSCILLATION_COMPONENT_SHA256
    assert len(OSCILLATION_COMPONENT_SHA256) == 10
    assert OSCILLATION_COMPONENT_SHA256[
        "proof_sources/oscillations/03_covariance_curl_and_tails.md"
    ] == "71432f2fffa382cfe4b8c325c51ee4a8caa3c2578cef4fc308f0db2497ca3474"
    assert ledger.parent_a2_head == PARENT_A2_HEAD
    assert ledger.parent_a2_source_blob == PARENT_A2_SOURCE_BLOB


def test_corrected_full_phase_and_required_background_are_explicit() -> None:
    ledger = default_kokuno_corrected_oscillation_source_ledger()

    assert ledger.phase == "Phi = p theta + p_z Z/epsilon + x_0 R - v H_Phi"
    assert ledger.phase_hamiltonian == "H_Phi = p F + p_z G"
    assert ledger.phase_vector == (
        "(n_Phi)_r = x_0 - v (H_Phi)_R",
        "(n_Phi)_theta = p/R",
        "(n_Phi)_z = p_z - epsilon v (H_Phi)_Z",
    )
    required = set(ledger.required_background_inputs)
    assert {"V", "G", "F = V/R", "F_R", "G_R", "g = (R F_R, G_R)"} <= required
    assert "representative leading shear g_0" in required
    assert "representative frame vectors N and K" in required
    assert "b=O(epsilon)" in ledger.imported_background_requirement
    assert "O(epsilon^2)" in ledger.imported_background_requirement


def test_constrained_amplitude_and_source_pressure_relation_are_recorded_without_pressure_promotion() -> None:
    ledger = default_kokuno_corrected_oscillation_source_ledger()

    assert ledger.damping == "d = epsilon k^2 |n_Phi|^2"
    assert ledger.constrained_amplitude_equation == (
        "t_m' + K t_m + m^2 d t_m + i k m n_Phi pi_m = -f_m"
    )
    assert ledger.constrained_amplitude_transversality == "n_Phi dot t_m = 0"
    assert ledger.source_pressure_relation.startswith("pi_m = (i/(k m))")
    assert ledger.matched_pressure_materialized is False


def test_complete_cylindrical_curl_retains_remainder_and_longitudinal_component() -> None:
    ledger = default_kokuno_corrected_oscillation_source_ledger()

    assert ledger.normalized_curl == (
        "(curl_* A)_r = R^(-1) partial_theta A_z - D_z A_theta",
        "(curl_* A)_theta = D_z A_r - D_r A_z",
        "(curl_* A)_z = (D_r + R^(-1)) A_theta - R^(-1) partial_theta A_r",
    )
    assert ledger.vector_potential_coefficient == (
        "C_m = i (n_Phi cross t_m) / (k m |n_Phi|^2)"
    )
    assert ledger.harmonic_vector_potential == "A_m = C_m exp(i k m Phi)"
    assert ledger.exact_curl_remainder == (
        "(r_m)_r = -D_z (C_m)_theta",
        "(r_m)_theta = D_z (C_m)_r - D_r (C_m)_z",
        "(r_m)_z = (D_r + R^(-1)) (C_m)_theta",
    )
    assert ledger.complete_amplitude == "a_m = t_m + r_m"
    assert ledger.exact_curl_divergence_identity == "div_*(curl_* A) = 0"
    assert ledger.complete_amplitude_asserted_transverse is False
    assert ledger.retained_longitudinal_source == "n_Phi dot a_m = n_Phi dot r_m"
    assert "i k m" in ledger.retained_longitudinal_identity
    assert "D_r + R^(-1)" in ledger.retained_longitudinal_identity


def test_source_support_and_cutoff_are_not_relabelled_as_repository_support() -> None:
    ledger = default_kokuno_corrected_oscillation_source_ledger()

    assert "zeta = exp" in ledger.source_support_envelope
    assert "X_a < X < X_b" in ledger.source_support_envelope
    assert "sqrt(zeta)" in ledger.source_wave_support
    assert ledger.cutoff_definition == "hat t_m = psi t_m; hat pi_m = psi pi_m"
    assert "(1-psi) f_m + psi' t_m" in ledger.exact_cutoff_residual
    assert ledger.cutoff_tail_location == "|v - L_s/2| >= L_s/5"
    assert ledger.current_runtime_support_source_exact is False


def test_current_repository_runtime_is_not_promoted_to_source_exact() -> None:
    ledger = default_kokuno_corrected_oscillation_source_ledger()
    truth = ledger.truth_boundary_payload()

    assert truth["source_formula_status"] == (
        "confirmed_in_corrected_208p_workbench_and_source_bundle"
    )
    assert truth["current_runtime_mapping_status"] == "not_verified"
    for key in (
        "current_runtime_phase_source_exact",
        "current_runtime_support_source_exact",
        "current_runtime_complete_curl_source_equivalence_verified",
        "source_to_runtime_parameter_map_complete",
        "current_runtime_paper_exact",
        "matched_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_ns_residual_assessed",
        "pde_validated",
    ):
        assert truth[key] is False


def test_earlier_source_contract_is_preserved_as_a_distinct_notation_capture() -> None:
    earlier = default_kokuno_oscillatory_source_contract()
    ledger = default_kokuno_corrected_oscillation_source_ledger()

    assert earlier.source_phase.startswith("Theta_beta =")
    assert ledger.phase.startswith("Phi =")
    assert earlier.source_phase != ledger.phase
    assert ledger.current_runtime_phase_source_exact is False
    assert ledger.current_runtime_complete_curl_source_equivalence_verified is False


def test_manifest_round_trip_is_deterministic(tmp_path: Path) -> None:
    ledger = default_kokuno_corrected_oscillation_source_ledger()
    path = tmp_path / "ledger.json"
    ledger.save_manifest(path)
    loaded = KokunoCorrectedOscillationSourceLedger.load_manifest(path)

    assert loaded.semantic_payload() == ledger.semantic_payload()
    assert loaded.semantic_sha256 == ledger.semantic_sha256
    assert json.loads(path.read_text(encoding="utf-8"))["semantic_sha256"] == ledger.semantic_sha256


def test_manifest_checksum_tamper_fails_closed(tmp_path: Path) -> None:
    ledger = default_kokuno_corrected_oscillation_source_ledger()
    path = tmp_path / "ledger.json"
    ledger.save_manifest(path)
    raw = json.loads(path.read_text(encoding="utf-8"))
    raw["payload"]["upstream_source"]["date"] = "2099-01-01"
    path.write_text(json.dumps(raw), encoding="utf-8")

    with pytest.raises(ValueError, match="checksum mismatch"):
        KokunoCorrectedOscillationSourceLedger.load_manifest(path)


def test_rehash_consistent_source_or_truth_mutation_fails_closed(tmp_path: Path) -> None:
    ledger = default_kokuno_corrected_oscillation_source_ledger()
    raw = ledger.manifest()

    mutated_source = copy.deepcopy(raw)
    mutated_source["payload"]["upstream_source"]["oscillation_component_sha256"][
        "proof_sources/oscillations/03_covariance_curl_and_tails.md"
    ] = "0" * 64
    mutated_source["semantic_sha256"] = _sha256(mutated_source["payload"])
    path = tmp_path / "mutated-source.json"
    path.write_text(json.dumps(mutated_source), encoding="utf-8")
    with pytest.raises(ValueError, match="manifest semantics"):
        KokunoCorrectedOscillationSourceLedger.load_manifest(path)

    promoted = copy.deepcopy(raw)
    promoted["payload"]["truth_boundary"]["current_runtime_phase_source_exact"] = True
    promoted["semantic_sha256"] = _sha256(promoted["payload"])
    path = tmp_path / "promoted.json"
    path.write_text(json.dumps(promoted), encoding="utf-8")
    with pytest.raises(ValueError, match="cannot promote"):
        KokunoCorrectedOscillationSourceLedger.load_manifest(path)


def test_module_scope_is_provenance_only() -> None:
    ledger = default_kokuno_corrected_oscillation_source_ledger()
    payload = ledger.semantic_payload()
    serialized = json.dumps(payload, sort_keys=True)

    assert "provenance/formula ledger only" in serialized
    assert not hasattr(ledger, "velocity")
    assert not hasattr(ledger, "pressure")
    assert not hasattr(ledger, "forcing")
    assert not hasattr(ledger, "optimize")
    assert ledger.pde_validated is False
