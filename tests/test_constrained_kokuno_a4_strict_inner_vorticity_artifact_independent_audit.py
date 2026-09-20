from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_strict_inner_vorticity_artifact_independent_audit import (
    FD6_SPATIAL_STEPS,
    FINAL_DIVERGENCE_GATE,
    FINAL_MOMENTUM_NORMALIZED_GATE,
    FINE_RELATIVE_RMS_GATE,
    FINE_RELATIVE_SAMPLED_MAX_GATE,
    HELDOUT_SEED,
    NONTRIVIAL_VELOCITY_RMS_FLOOR,
    component_error_metrics,
    fd6_level,
    fd6_vorticity,
    frozen_gate_failures,
    public_contract,
    relative_rms_error,
    relative_sampled_max_error,
    vector_rms,
)
from openai_ns_reconstruction.kokuno_strict_inner_leading_oscillatory_vorticity_candidate import (
    KokunoStrictInnerLeadingOscillatoryVorticityCandidate,
)


def _canonical_sha(payload: dict[str, object]) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _manufactured_velocity(x, y, z, t):
    x, y, z, t = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    return np.stack((y**3 + z + 0.2 * t, z**3 + x, x**3 + y - 0.1 * t), axis=-1)


def _manufactured_vorticity(x, y, z, t):
    x, y, z, _ = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    return np.stack((1.0 - 3.0 * z**2, 1.0 - 3.0 * x**2, 1.0 - 3.0 * y**2), axis=-1)


class _ToyBackend:
    field_sha256 = "1" * 64
    temporal_derivative_sha256 = "2" * 64
    spatial_derivative_sha256 = "3" * 64

    @staticmethod
    def _arrays(x, y, z, t):
        return np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )

    def velocity(self, x, y, z, t):
        x, y, z, t = self._arrays(x, y, z, t)
        return np.stack((-y, x, 0.25 * z + 0.1 * t), axis=-1)

    def velocity_dt(self, x, y, z, t):
        x, y, z, t = self._arrays(x, y, z, t)
        out = np.zeros(x.shape + (3,), dtype=float)
        out[..., 2] = 0.1
        return out

    def velocity_jacobian(self, x, y, z, t):
        x, y, z, t = self._arrays(x, y, z, t)
        out = np.zeros(x.shape + (3, 3), dtype=float)
        out[..., 0, 1] = -1.0
        out[..., 1, 0] = 1.0
        out[..., 2, 2] = 0.25
        return out


def test_fd6_curl_is_exact_on_degree_three_manufactured_field():
    x = np.asarray([-0.31, -0.07, 0.19, 0.43])
    y = np.asarray([0.22, -0.28, 0.37, -0.11])
    z = np.asarray([-0.26, 0.14, 0.33, -0.41])
    t = np.asarray([0.45, 0.48, 0.52, 0.55])
    expected = _manufactured_vorticity(x, y, z, t)
    observed = fd6_vorticity(_manufactured_velocity, x, y, z, t, step=8.0e-4)
    np.testing.assert_allclose(observed, expected, rtol=0.0, atol=2.0e-12)


def test_fd6_three_level_gate_and_mutations():
    rng = np.random.default_rng(HELDOUT_SEED)
    x = rng.uniform(-0.45, 0.45, 64)
    y = rng.uniform(-0.45, 0.45, 64)
    z = rng.uniform(-0.45, 0.45, 64)
    t = rng.uniform(0.44, 0.56, 64)
    production = _manufactured_vorticity(x, y, z, t)
    levels = tuple(
        fd6_level(_manufactured_velocity, x, y, z, t, step=step)
        for step in FD6_SPATIAL_STEPS
    )
    velocity_rms = vector_rms(_manufactured_velocity(x, y, z, t))
    assert velocity_rms > NONTRIVIAL_VELOCITY_RMS_FLOOR
    assert frozen_gate_failures(production, levels, velocity_rms=velocity_rms) == []

    scaled = production * 0.90
    scaled_failures = frozen_gate_failures(scaled, levels, velocity_rms=velocity_rms)
    assert "fine_relative_rms" in scaled_failures
    assert "fine_relative_sampled_max" in scaled_failures

    flipped = np.array(production, copy=True)
    flipped[..., 2] *= -1.0
    flipped_failures = frozen_gate_failures(flipped, levels, velocity_rms=velocity_rms)
    assert "fine_relative_rms" in flipped_failures


def test_error_metrics_are_componentwise_and_fixed():
    ref = np.asarray([[1.0, 2.0, 3.0], [2.0, -1.0, 4.0]])
    obs = ref + np.asarray([[0.01, 0.0, -0.02], [-0.01, 0.03, 0.0]])
    assert relative_rms_error(ref, obs) > 0.0
    assert relative_sampled_max_error(ref, obs) > 0.0
    metrics = component_error_metrics(ref, obs)
    assert [entry["component"] for entry in metrics] == ["omega_x", "omega_y", "omega_z"]
    assert all(float(entry["absolute_rms"]) >= 0.0 for entry in metrics)


def test_public_contract_locks_independent_path_and_final_gates():
    contract = public_contract()
    assert contract["heldout_seed"] == HELDOUT_SEED
    assert contract["fd6_spatial_steps"] == list(FD6_SPATIAL_STEPS)
    independent = contract["independent_reference"]
    assert independent["artifact_must_be_saved_and_reloaded"] is True
    assert independent["production_surface"] == "reloaded public vorticity(x,y,z,t)"
    assert independent["numerical_reference_uses"] == "reloaded public velocity(x,y,z,t) only"
    assert independent["agent1_analytic_spatial_derivatives_used"] is False
    assert independent["agent2_production_jacobian_used"] is False
    assert independent["agent2_fd4_verifier_used"] is False
    assert contract["frozen_gates"]["fine_relative_rms"] == FINE_RELATIVE_RMS_GATE
    assert contract["frozen_gates"]["fine_relative_sampled_max"] == FINE_RELATIVE_SAMPLED_MAX_GATE
    assert contract["frozen_gates"]["final_momentum_normalized"] == FINAL_MOMENTUM_NORMALIZED_GATE == 1.0e-3
    assert contract["frozen_gates"]["final_divergence"] == FINAL_DIVERGENCE_GATE == 1.0e-5
    truth = contract["truth_boundary"]
    assert truth["leading_only_ns_residual_assessed"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False


def test_vorticity_artifact_manifest_rejects_checksum_valid_semantic_and_gate_mutations(tmp_path: Path):
    backend = _ToyBackend()
    candidate = KokunoStrictInnerLeadingOscillatoryVorticityCandidate(backend)
    path = tmp_path / "candidate.json"
    manifest = candidate.save_manifest(path)
    rebound = KokunoStrictInnerLeadingOscillatoryVorticityCandidate.load_manifest(backend, path)
    assert rebound.vorticity_sha256 == candidate.vorticity_sha256

    semantic = json.loads(json.dumps(manifest))
    semantic["vorticity_sha256"] = "f" * 64
    unhashed = dict(semantic)
    unhashed.pop("manifest_sha256", None)
    semantic["manifest_sha256"] = _canonical_sha(unhashed)
    semantic_path = tmp_path / "semantic.json"
    semantic_path.write_text(json.dumps(semantic))
    with pytest.raises(ValueError, match="vorticity_sha256 changed"):
        KokunoStrictInnerLeadingOscillatoryVorticityCandidate.load_manifest(backend, semantic_path)

    gate = json.loads(json.dumps(manifest))
    gate["scientific_gates"]["momentum_max_l2"] = 1.001e-3
    unhashed = dict(gate)
    unhashed.pop("manifest_sha256", None)
    gate["manifest_sha256"] = _canonical_sha(unhashed)
    gate_path = tmp_path / "gate.json"
    gate_path.write_text(json.dumps(gate))
    with pytest.raises(ValueError, match="scientific_gates changed"):
        KokunoStrictInnerLeadingOscillatoryVorticityCandidate.load_manifest(backend, gate_path)
