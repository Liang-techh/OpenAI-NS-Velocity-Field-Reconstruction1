from __future__ import annotations

import hashlib
import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_oscillatory_nonlinear_cancellation import (
    KokunoOscillatoryNonlinearCancellation,
    public_contract,
)


class PolynomialInnerBackend:
    field_sha256 = "3" * 64

    @staticmethod
    def _broadcast(x, y, z, t):
        return np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )

    def velocity(self, x, y, z, t):
        x, y, z, _ = self._broadcast(x, y, z, t)
        return np.stack((0.19 * x + 0.08, -0.15 * y + 0.03, 0.12 * z - 0.02), axis=-1)

    def velocity_dt(self, x, y, z, t):
        x, _, _, _ = self._broadcast(x, y, z, t)
        return np.zeros(x.shape + (3,), dtype=float)

    def self_advection(self, x, y, z, t):
        u = self.velocity(x, y, z, t)
        return np.stack((0.19 * u[..., 0], -0.15 * u[..., 1], 0.12 * u[..., 2]), axis=-1)

    def velocity_laplacian(self, x, y, z, t):
        return self.velocity_dt(x, y, z, t)


def _probes():
    return (
        np.asarray([0.43, 0.56, 0.69, 0.82, 0.97, 1.09]),
        np.asarray([0.14, -0.19, 0.22, -0.16, 0.11, -0.08]),
        np.asarray([-0.46, -0.24, -0.01, 0.23, 0.48, 0.71]),
        np.asarray([0.474, 0.484, 0.494, 0.504, 0.514, 0.524]),
    )


def _canonical_sha(payload):
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def test_cancellation_metrics_replay_raw_piece_geometry_and_parent_attribution():
    diagnostic = KokunoOscillatoryNonlinearCancellation(PolynomialInnerBackend())
    result = diagnostic.evaluate(*_probes())

    cosines = np.asarray([
        result.inner_to_osc_vs_osc_to_inner_cosine,
        result.inner_to_osc_vs_self_cosine,
        result.osc_to_inner_vs_self_cosine,
        result.mixed_vs_self_cosine,
    ])
    assert np.all(np.isfinite(cosines))
    assert np.all(cosines >= -1.0)
    assert np.all(cosines <= 1.0)
    assert 0.0 <= result.mixed_cancellation_ratio <= 1.0
    assert 0.0 <= result.aggregate_cancellation_ratio <= 1.0
    assert result.gram_reconstruction_relative_error <= 1.0e-12
    assert result.parent_attribution_replay_relative_error <= 1.0e-14


def test_opposite_vectors_show_cancellation_without_being_a_residual_fraction():
    # Lock the mathematical interpretation independently of the production handoff:
    # cancellation ratio ||a+b||/(||a||+||b||) can be zero even with nonzero pieces.
    a = np.asarray([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
    b = -a
    norm = lambda v: float(np.sqrt(np.mean(np.sum(v * v, axis=-1))))
    ratio = norm(a + b) / (norm(a) + norm(b))
    cosine = float(np.mean(np.sum(a * b, axis=-1)) / (norm(a) * norm(b)))
    assert ratio == pytest.approx(0.0)
    assert cosine == pytest.approx(-1.0)

    metrics = KokunoOscillatoryNonlinearCancellation(PolynomialInnerBackend()).semantic_payload()[
        "metric_contract"
    ]
    assert metrics["ratios_are_raw_cartesian_cancellation_diagnostics"] is True
    assert metrics["ratios_are_m0_mean_defect_fractions"] is False
    assert metrics["ratios_are_ns_residual_fractions"] is False
    assert metrics["ratios_are_acceptance_gates"] is False


def test_contract_stays_out_of_agent3_and_full_residual_lanes():
    contract = public_contract()
    assert contract["forbidden_evaluate_inputs_present"] == []
    assert contract["mean_projection_performed"] is False
    assert contract["correction_velocity_constructed"] is False
    assert contract["cancellation_metrics_are_acceptance_gates"] is False
    assert contract["complete_ns_residual"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False

    payload = KokunoOscillatoryNonlinearCancellation(PolynomialInnerBackend()).semantic_payload()
    boundary = payload["scientific_boundary"]
    assert boundary["mean_projection_performed"] is False
    assert boundary["radial_inverse_performed"] is False
    assert boundary["correction_velocity_constructed"] is False
    assert boundary["same_protocol_st006_comparison_valid"] is False
    assert boundary["residual_reduction_claimed"] is False


def test_manifest_round_trip_and_checksum_valid_semantic_tamper_fail_closed(tmp_path):
    backend = PolynomialInnerBackend()
    diagnostic = KokunoOscillatoryNonlinearCancellation(backend)
    path = tmp_path / "nonlinear_cancellation.json"
    diagnostic.save_manifest(path)
    rebound = KokunoOscillatoryNonlinearCancellation.load_manifest(backend, path)
    assert rebound.cancellation_sha256 == diagnostic.cancellation_sha256

    raw = json.loads(path.read_text())
    raw["payload"]["metric_contract"]["ratios_are_acceptance_gates"] = True
    raw["cancellation_sha256"] = _canonical_sha(raw["payload"])
    path.write_text(json.dumps(raw))
    with pytest.raises(ValueError, match="semantics mismatch"):
        KokunoOscillatoryNonlinearCancellation.load_manifest(backend, path)
