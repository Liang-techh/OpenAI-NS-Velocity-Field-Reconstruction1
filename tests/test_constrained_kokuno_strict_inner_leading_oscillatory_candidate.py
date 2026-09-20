from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_z_pullback_velocity import velocity_osc
from openai_ns_reconstruction.kokuno_strict_inner_leading_oscillatory_candidate import (
    AGENT1_REFERENCE_BLOB,
    AGENT1_REFERENCE_HEAD,
    AGENT1_REFERENCE_PR,
    OSCILLATORY_VELOCITY_BLOB,
    PARENT_AGENT2_HEAD,
    KokunoStrictInnerLeadingOscillatoryCandidate,
    public_contract,
)


class LinearStrictInnerBackend:
    def __init__(self, field_sha256: str = "a" * 64) -> None:
        self.field_sha256 = field_sha256

    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        if np.any(x > 0.9):
            raise ValueError("strict-inner domain exceeded")
        return np.stack((x + t, y - 2.0 * t, z + 0.5 * t), axis=-1)


class WrongShapeBackend(LinearStrictInnerBackend):
    def velocity(self, x, y, z, t):
        x = np.asarray(x, dtype=float)
        return np.zeros(x.shape + (2,), dtype=float)


def sample_points():
    return (
        np.asarray([0.31, 0.42, 0.53]),
        np.asarray([0.27, -0.34, 0.21]),
        np.asarray([0.18, -0.22, 0.11]),
        np.asarray([0.46, 0.50, 0.54]),
    )


def test_candidate_velocity_is_exact_additive_composition():
    backend = LinearStrictInnerBackend()
    candidate = KokunoStrictInnerLeadingOscillatoryCandidate(backend)
    x, y, z, t = sample_points()

    evaluated = candidate.evaluate(x, y, z, t)
    expected_inner = backend.velocity(x, y, z, t)
    expected_osc = velocity_osc(x, y, z, t)

    np.testing.assert_array_equal(evaluated.inner_leading_velocity, expected_inner)
    np.testing.assert_array_equal(evaluated.oscillatory_velocity, expected_osc)
    np.testing.assert_array_equal(evaluated.velocity, expected_inner + expected_osc)
    np.testing.assert_array_equal(candidate.velocity(x, y, z, t), evaluated.velocity)
    assert evaluated.additive_closure_abs_max <= 1.0e-14
    assert np.all(np.isfinite(evaluated.velocity))


def test_candidate_manifest_roundtrip_binds_both_field_identities(tmp_path):
    backend = LinearStrictInnerBackend()
    candidate = KokunoStrictInnerLeadingOscillatoryCandidate(backend)
    path = tmp_path / "candidate.json"
    manifest = candidate.save_manifest(path)

    reloaded = KokunoStrictInnerLeadingOscillatoryCandidate.load_manifest(backend, path)
    assert reloaded.candidate_sha256 == candidate.candidate_sha256
    assert reloaded.field_configuration() == candidate.field_configuration()
    assert manifest["field_configuration"]["inner_leading"]["field_sha256"] == "a" * 64
    assert manifest["field_configuration"]["oscillatory"]["field_sha256"] == candidate.oscillatory_field_sha256
    assert len(manifest["candidate_sha256"]) == 64
    assert len(manifest["manifest_sha256"]) == 64


def test_manifest_rejects_backend_drift_and_payload_tamper(tmp_path):
    backend = LinearStrictInnerBackend()
    candidate = KokunoStrictInnerLeadingOscillatoryCandidate(backend)
    path = tmp_path / "candidate.json"
    candidate.save_manifest(path)

    with pytest.raises(ValueError, match="field configuration/provenance changed"):
        KokunoStrictInnerLeadingOscillatoryCandidate.load_manifest(
            LinearStrictInnerBackend("b" * 64), path
        )

    payload = json.loads(path.read_text())
    payload["truth_boundary"]["pde_validated"] = True
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="manifest sha256 mismatch"):
        KokunoStrictInnerLeadingOscillatoryCandidate.load_manifest(backend, path)


def test_candidate_preserves_inner_domain_fail_closed():
    candidate = KokunoStrictInnerLeadingOscillatoryCandidate(LinearStrictInnerBackend())
    with pytest.raises(ValueError, match="strict-inner domain exceeded"):
        candidate.velocity(1.0, 0.0, 0.0, 0.5)


def test_candidate_rejects_bad_backend_identity_and_shape():
    with pytest.raises(ValueError, match="lowercase hexadecimal sha256"):
        KokunoStrictInnerLeadingOscillatoryCandidate(
            LinearStrictInnerBackend("not-a-sha")
        )

    candidate = KokunoStrictInnerLeadingOscillatoryCandidate(WrongShapeBackend())
    with pytest.raises(ValueError, match="inner leading velocity must have shape"):
        candidate.velocity(0.3, 0.2, 0.1, 0.5)


def test_public_contract_is_pinned_and_fail_closed():
    contract = public_contract()
    assert contract["parent_agent2_head"] == PARENT_AGENT2_HEAD
    assert contract["agent1_reference_pr"] == AGENT1_REFERENCE_PR == 848
    assert contract["agent1_reference_head"] == AGENT1_REFERENCE_HEAD
    assert contract["agent1_reference_blob"] == AGENT1_REFERENCE_BLOB
    assert contract["oscillatory_velocity_blob"] == OSCILLATORY_VELOCITY_BLOB
    assert contract["forbidden_velocity_inputs_present"] == []

    truth = contract["truth_boundary"]
    assert truth["strict_inner_composite_velocity_executable"] is True
    assert truth["oscillatory_complete_curl_path_inherited"] is True
    for key in (
        "inner_leading_is_final_corrected_fixed_point",
        "global_leading_velocity_materialized",
        "outer_join_materialized",
        "matched_pressure_included",
        "restricted_forcing_included",
        "agent3_correction_velocity_included",
        "complete_kokuno_candidate_assembled",
        "complete_ns_momentum_residual_formed",
        "heldout_ns_residual_assessed",
        "st006_same_protocol_comparison_valid",
        "residual_reduction_claimed",
        "paper_exact",
        "pde_validated",
    ):
        assert truth[key] is False
