from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_oscillatory_time_derivative import (
    velocity_osc_dt,
)
from openai_ns_reconstruction.kokuno_strict_inner_leading_oscillatory_candidate import (
    KokunoStrictInnerLeadingOscillatoryCandidate,
)
from openai_ns_reconstruction.kokuno_strict_inner_leading_oscillatory_differentiable_candidate import (
    A2_TIME_COMPOSITION_BLOB,
    AGENT1_REFERENCE_BLOB,
    AGENT1_REFERENCE_HEAD,
    AGENT1_REFERENCE_PR,
    OSCILLATORY_TIME_DERIVATIVE_BLOB,
    PARENT_AGENT2_HEAD,
    KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate,
    public_contract,
)


class LinearDifferentiableBackend:
    def __init__(
        self,
        field_sha256: str = "a" * 64,
        temporal_derivative_sha256: str = "b" * 64,
    ) -> None:
        self.field_sha256 = field_sha256
        self.temporal_derivative_sha256 = temporal_derivative_sha256

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

    def velocity_dt(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        if np.any(x > 0.9):
            raise ValueError("strict-inner domain exceeded")
        return np.stack(
            (np.ones_like(t), -2.0 * np.ones_like(t), 0.5 * np.ones_like(t)),
            axis=-1,
        )


class MissingTimeBackend:
    field_sha256 = "a" * 64
    temporal_derivative_sha256 = "b" * 64

    def velocity(self, x, y, z, t):
        x = np.asarray(x, dtype=float)
        return np.zeros(x.shape + (3,), dtype=float)


class WrongTimeShapeBackend(LinearDifferentiableBackend):
    def velocity_dt(self, x, y, z, t):
        x = np.asarray(x, dtype=float)
        return np.zeros(x.shape + (2,), dtype=float)


def sample_points():
    return (
        np.asarray([0.31, 0.42, 0.53]),
        np.asarray([0.27, -0.34, 0.21]),
        np.asarray([0.18, -0.22, 0.11]),
        np.asarray([0.46, 0.50, 0.54]),
    )


def test_differentiable_candidate_preserves_base_velocity_and_adds_time_surface():
    backend = LinearDifferentiableBackend()
    base = KokunoStrictInnerLeadingOscillatoryCandidate(backend)
    candidate = KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate(backend)
    x, y, z, t = sample_points()

    evaluated = candidate.evaluate(x, y, z, t)
    np.testing.assert_array_equal(evaluated.velocity, base.velocity(x, y, z, t))
    np.testing.assert_array_equal(candidate.velocity(x, y, z, t), base.velocity(x, y, z, t))
    np.testing.assert_array_equal(
        evaluated.inner_leading_velocity_dt,
        backend.velocity_dt(x, y, z, t),
    )
    np.testing.assert_array_equal(
        evaluated.oscillatory_velocity_dt,
        velocity_osc_dt(x, y, z, t),
    )
    np.testing.assert_array_equal(
        evaluated.velocity_dt,
        backend.velocity_dt(x, y, z, t) + velocity_osc_dt(x, y, z, t),
    )
    np.testing.assert_array_equal(candidate.velocity_dt(x, y, z, t), evaluated.velocity_dt)
    assert candidate.velocity_candidate_sha256 == base.candidate_sha256
    assert evaluated.velocity_replay_abs_max <= 1.0e-14
    assert evaluated.velocity_dt_additive_closure_abs_max <= 1.0e-13


def test_differentiable_manifest_roundtrip_binds_time_identity(tmp_path):
    backend = LinearDifferentiableBackend()
    candidate = KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate(backend)
    path = tmp_path / "candidate-dt.json"
    manifest = candidate.save_manifest(path)

    rebound = KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate.load_manifest(
        backend, path
    )
    assert rebound.velocity_candidate_sha256 == candidate.velocity_candidate_sha256
    assert rebound.differentiable_sha256 == candidate.differentiable_sha256
    assert (
        manifest["differentiable_configuration"]["inner_leading_time_derivative"][
            "temporal_derivative_sha256"
        ]
        == "b" * 64
    )
    assert len(manifest["differentiable_sha256"]) == 64
    assert len(manifest["manifest_sha256"]) == 64


def test_manifest_rejects_temporal_identity_drift_and_tamper(tmp_path):
    backend = LinearDifferentiableBackend()
    candidate = KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate(backend)
    path = tmp_path / "candidate-dt.json"
    candidate.save_manifest(path)

    with pytest.raises(ValueError, match="configuration/provenance changed"):
        KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate.load_manifest(
            LinearDifferentiableBackend(temporal_derivative_sha256="c" * 64), path
        )

    payload = json.loads(path.read_text())
    payload["truth_boundary"]["pde_validated"] = True
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="manifest sha256 mismatch"):
        KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate.load_manifest(
            backend, path
        )


def test_candidate_preserves_strict_inner_fail_closed_for_velocity_and_dt():
    candidate = KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate(
        LinearDifferentiableBackend()
    )
    with pytest.raises(ValueError, match="strict-inner domain exceeded"):
        candidate.velocity(1.0, 0.0, 0.0, 0.5)
    with pytest.raises(ValueError, match="strict-inner domain exceeded"):
        candidate.velocity_dt(1.0, 0.0, 0.0, 0.5)


def test_candidate_rejects_missing_or_malformed_temporal_backend():
    with pytest.raises(TypeError, match="velocity_dt"):
        KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate(MissingTimeBackend())

    with pytest.raises(ValueError, match="temporal_derivative_sha256"):
        KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate(
            LinearDifferentiableBackend(temporal_derivative_sha256="not-a-sha")
        )

    candidate = KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate(
        WrongTimeShapeBackend()
    )
    with pytest.raises(ValueError, match="velocity_dt must have shape"):
        candidate.velocity_dt(0.3, 0.2, 0.1, 0.5)


def test_public_contract_is_pinned_and_fail_closed():
    contract = public_contract()
    assert contract["parent_agent2_head"] == PARENT_AGENT2_HEAD
    assert contract["agent1_reference_pr"] == AGENT1_REFERENCE_PR == 848
    assert contract["agent1_reference_head"] == AGENT1_REFERENCE_HEAD
    assert contract["agent1_reference_blob"] == AGENT1_REFERENCE_BLOB
    assert contract["a2_time_composition_blob"] == A2_TIME_COMPOSITION_BLOB
    assert (
        contract["oscillatory_time_derivative_blob"]
        == OSCILLATORY_TIME_DERIVATIVE_BLOB
    )
    assert contract["forbidden_velocity_inputs_present"] == []
    assert contract["forbidden_velocity_dt_inputs_present"] == []

    truth = contract["truth_boundary"]
    assert truth["strict_inner_composite_velocity_executable"] is True
    assert truth["strict_inner_composite_velocity_dt_executable"] is True
    assert truth["base_velocity_candidate_identity_preserved"] is True
    assert truth["production_velocity_dt_uses_finite_difference"] is False
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
