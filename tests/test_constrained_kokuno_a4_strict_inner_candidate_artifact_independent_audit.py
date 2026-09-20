from __future__ import annotations

import inspect
import json

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_a4_strict_inner_candidate_artifact_independent_audit as audit
from openai_ns_reconstruction.kokuno_strict_inner_leading_oscillatory_candidate import (
    KokunoStrictInnerLeadingOscillatoryCandidate,
)


def test_richardson_derivative_convention_on_cubic_vector_field() -> None:
    def velocity(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(x, y, z, t)
        return np.stack((x**3 + y*t, y**3 + z*t, z**3 + x*t), axis=-1)

    x = np.array([0.17, -0.23, 0.31])
    y = np.array([-0.11, 0.29, 0.07])
    z = np.array([0.21, -0.13, -0.27])
    t = np.array([0.41, 0.53, 0.47])
    level = audit.richardson_divergence_level(velocity, x, y, z, t, step=8.0e-4)

    expected = 3.0 * (x*x + y*y + z*z)
    np.testing.assert_allclose(level.divergence, expected, rtol=1.0e-10, atol=1.0e-11)
    np.testing.assert_allclose(level.jacobian[:, 0, 1], t, rtol=1.0e-10, atol=1.0e-11)
    np.testing.assert_allclose(level.jacobian[:, 1, 2], t, rtol=1.0e-10, atol=1.0e-11)
    np.testing.assert_allclose(level.jacobian[:, 2, 0], t, rtol=1.0e-10, atol=1.0e-11)


def test_frozen_divergence_gates_accept_linear_divergence_free_field() -> None:
    def velocity(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(x, y, z, t)
        return np.stack((y + 0.1*t, -x + z, -y - 0.1*t), axis=-1)

    x = np.array([0.13, -0.19, 0.27, -0.31])
    y = np.array([-0.23, 0.17, -0.09, 0.21])
    z = np.array([0.11, 0.29, -0.25, -0.15])
    t = np.array([0.45, 0.49, 0.53, 0.55])
    levels = tuple(
        audit.richardson_divergence_level(velocity, x, y, z, t, step=h)
        for h in audit.RICHARDSON_STEPS
    )
    assert audit.frozen_gate_failures(levels) == []


def test_divergence_mutation_is_caught_by_fixed_gate() -> None:
    def mutated(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(x, y, z, t)
        return np.stack((y + 1.0e-3*x, -x, np.zeros_like(z)), axis=-1)

    x = np.array([0.14, -0.18, 0.26])
    y = np.array([-0.22, 0.16, 0.08])
    z = np.array([0.12, 0.20, -0.24])
    t = np.array([0.46, 0.50, 0.54])
    levels = tuple(
        audit.richardson_divergence_level(mutated, x, y, z, t, step=h)
        for h in audit.RICHARDSON_STEPS
    )
    failures = audit.frozen_gate_failures(levels)
    assert "divergence_sampled_max" in failures
    assert "divergence_sampled_rms" in failures


class _FailClosedToyBackend:
    field_sha256 = "0" * 64

    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float), np.asarray(y, dtype=float),
            np.asarray(z, dtype=float), np.asarray(t, dtype=float),
        )
        if np.any(np.abs(x) > 1.0):
            raise ValueError("outside toy strict-inner domain")
        return np.stack((y, -x, np.zeros_like(z)), axis=-1)


def test_candidate_manifest_roundtrip_and_tamper_rejection(tmp_path) -> None:
    backend = _FailClosedToyBackend()
    candidate = KokunoStrictInnerLeadingOscillatoryCandidate(backend)
    path = tmp_path / "candidate.json"
    manifest = candidate.save_manifest(path)
    rebound = KokunoStrictInnerLeadingOscillatoryCandidate.load_manifest(backend, path)
    assert rebound.candidate_sha256 == candidate.candidate_sha256

    tampered = dict(manifest)
    tampered["candidate_sha256"] = "f" * 64
    path.write_text(json.dumps(tampered))
    with pytest.raises(ValueError, match="manifest sha256 mismatch"):
        KokunoStrictInnerLeadingOscillatoryCandidate.load_manifest(backend, path)


def test_candidate_propagates_backend_domain_failure() -> None:
    candidate = KokunoStrictInnerLeadingOscillatoryCandidate(_FailClosedToyBackend())
    with pytest.raises(ValueError, match="outside toy strict-inner domain"):
        candidate.velocity(2.0, 0.0, 0.0, 0.5)


def test_a4_source_does_not_use_candidate_derivative_or_component_paths() -> None:
    source = inspect.getsource(audit)
    forbidden = (
        ".velocity_dt(", ".velocity_jacobian(", ".velocity_laplacian(",
        ".self_advection(", ".components(", ".evaluate(",
    )
    for token in forbidden:
        assert token not in source


def test_public_contract_keeps_full_pde_truth_boundary_false() -> None:
    contract = audit.public_contract()
    assert contract["agent2_candidate_head"] == audit.AGENT2_CANDIDATE_HEAD
    assert contract["frozen_gates"]["final_momentum_normalized_max_l2"] == 1.0e-3
    assert contract["frozen_gates"]["final_divergence_max_l2"] == 1.0e-5
    truth = contract["truth_boundary"]
    assert truth["strict_inner_candidate_artifact_independently_consumed"] is True
    assert truth["strict_inner_divergence_independently_assessed"] is True
    for key in (
        "global_leading_velocity_materialized", "outer_join_materialized",
        "matched_pressure_included", "restricted_forcing_included",
        "correction_velocity_included", "complete_kokuno_candidate_assembled",
        "leading_only_ns_residual_assessed", "leading_plus_oscillatory_ns_residual_assessed",
        "after_correction_ns_residual_assessed", "heldout_ns_momentum_residual_assessed",
        "formal_full_domain_pde_gate_assessed", "whole_domain_boundary_support_gate_assessed",
        "pde_validated", "paper_exact",
    ):
        assert truth[key] is False
