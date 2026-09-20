from __future__ import annotations

import json

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_oscillatory_transport_handoff as handoff_mod
from openai_ns_reconstruction.kokuno_public_inner_leading_oscillatory_transport_delta import (
    InnerLeadingOscillatoryTransportDeltaResult,
)


class FakeBackend:
    field_sha256 = "a" * 64

    @staticmethod
    def _zeros(x, y, z, t):
        shape = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )[0].shape
        return np.zeros(shape + (3,), dtype=float)

    def velocity(self, x, y, z, t):
        return self._zeros(x, y, z, t)

    def velocity_dt(self, x, y, z, t):
        return self._zeros(x, y, z, t)

    def self_advection(self, x, y, z, t):
        return self._zeros(x, y, z, t)

    def velocity_laplacian(self, x, y, z, t):
        return self._zeros(x, y, z, t)


def synthetic_parent_result() -> InnerLeadingOscillatoryTransportDeltaResult:
    time = np.asarray([[1.0, 2.0, 3.0], [0.5, -0.5, 1.0]])
    nonlinear = np.asarray([[0.2, -0.1, 0.3], [0.4, 0.2, -0.2]])
    viscous = np.asarray([[-0.05, 0.01, 0.02], [0.03, -0.04, 0.01]])
    delta = time + nonlinear + viscous
    inner = np.asarray([[2.0, 0.0, 1.0], [1.0, -1.0, 0.5]])
    return InnerLeadingOscillatoryTransportDeltaResult(
        inner_leading_transport=inner,
        inner_plus_oscillatory_transport=inner + delta,
        oscillatory_transport_increment=delta,
        oscillatory_time_increment=time,
        oscillatory_nonlinear_increment=nonlinear,
        oscillatory_viscous_increment=viscous,
        inner_leading_velocity=np.ones((2, 3)),
        oscillatory_velocity=np.full((2, 3), 0.25),
        inner_plus_oscillatory_velocity=np.full((2, 3), 1.25),
        viscosity=0.01,
    )


def test_public_contract_is_raw_agent3_handoff_and_not_residual():
    contract = handoff_mod.public_contract()
    assert contract["schema"] == handoff_mod.SCHEMA
    assert contract["forbidden_evaluate_inputs_present"] == []
    assert contract["handoff_quantity"] == "raw Cartesian oscillatory transport increment"
    assert contract["consumer"] == "Kokuno Agent 3"
    assert contract["mean_projection_performed"] is False
    assert contract["correction_velocity_constructed"] is False
    assert contract["complete_ns_residual"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False


def test_evaluate_preserves_parent_delta_and_decomposition(monkeypatch):
    parent = synthetic_parent_result()
    monkeypatch.setattr(
        handoff_mod,
        "evaluate_inner_leading_oscillatory_transport_delta",
        lambda backend, x, y, z, t: parent,
    )
    handoff = handoff_mod.KokunoOscillatoryTransportHandoff(FakeBackend())
    out = handoff.evaluate([0.1, 0.2], [0.0, 0.0], [0.1, 0.2], [0.5, 0.5])
    np.testing.assert_array_equal(
        out.oscillatory_transport_increment,
        parent.oscillatory_transport_increment,
    )
    np.testing.assert_array_equal(
        out.oscillatory_nonlinear_increment,
        parent.oscillatory_nonlinear_increment,
    )
    np.testing.assert_array_equal(
        handoff.transport_delta([0.1, 0.2], 0.0, [0.1, 0.2], 0.5),
        parent.oscillatory_transport_increment,
    )
    np.testing.assert_array_equal(
        handoff.nonlinear_increment([0.1, 0.2], 0.0, [0.1, 0.2], 0.5),
        parent.oscillatory_nonlinear_increment,
    )
    assert out.decomposition_closure_abs_max == 0.0
    assert out.viscosity == 0.01


def test_backend_contract_fails_closed():
    class MissingFieldSha(FakeBackend):
        field_sha256 = "not-a-sha"

    with pytest.raises(TypeError, match="field_sha256"):
        handoff_mod.KokunoOscillatoryTransportHandoff(MissingFieldSha())

    class MissingLaplacian:
        field_sha256 = "b" * 64
        velocity = FakeBackend.velocity
        velocity_dt = FakeBackend.velocity_dt
        self_advection = FakeBackend.self_advection

    with pytest.raises(TypeError, match="velocity_laplacian"):
        handoff_mod.KokunoOscillatoryTransportHandoff(MissingLaplacian())


def test_manifest_roundtrip_and_semantic_tamper_rejection(tmp_path):
    backend = FakeBackend()
    handoff = handoff_mod.KokunoOscillatoryTransportHandoff(backend)
    path = tmp_path / "handoff.json"
    manifest = handoff.save_manifest(path)
    rebound = handoff_mod.KokunoOscillatoryTransportHandoff.load_manifest(backend, path)
    assert rebound.handoff_sha256 == handoff.handoff_sha256
    assert manifest["handoff_sha256"] == handoff.handoff_sha256
    assert manifest["payload"]["source_contract_sha256"] == handoff.source_contract_sha256

    raw = json.loads(path.read_text())
    raw["payload"]["quantity_contract"]["complete_ns_residual"] = True
    raw["handoff_sha256"] = handoff_mod._sha256(raw["payload"])
    path.write_text(json.dumps(raw))
    with pytest.raises(ValueError, match="semantics mismatch"):
        handoff_mod.KokunoOscillatoryTransportHandoff.load_manifest(backend, path)


def test_decomposition_mismatch_is_rejected():
    parent = synthetic_parent_result()
    bad = InnerLeadingOscillatoryTransportDeltaResult(
        inner_leading_transport=parent.inner_leading_transport,
        inner_plus_oscillatory_transport=parent.inner_plus_oscillatory_transport,
        oscillatory_transport_increment=parent.oscillatory_transport_increment + 1.0e-6,
        oscillatory_time_increment=parent.oscillatory_time_increment,
        oscillatory_nonlinear_increment=parent.oscillatory_nonlinear_increment,
        oscillatory_viscous_increment=parent.oscillatory_viscous_increment,
        inner_leading_velocity=parent.inner_leading_velocity,
        oscillatory_velocity=parent.oscillatory_velocity,
        inner_plus_oscillatory_velocity=parent.inner_plus_oscillatory_velocity,
        viscosity=parent.viscosity,
    )
    with pytest.raises(RuntimeError, match="decomposition failed closure"):
        handoff_mod.KokunoOscillatoryTransportHandoff._from_parent_result(bad)


def test_semantic_payload_keeps_agent3_lane_unimplemented():
    handoff = handoff_mod.KokunoOscillatoryTransportHandoff(FakeBackend())
    payload = handoff.semantic_payload()
    quantity = payload["quantity_contract"]
    truth = payload["scientific_boundary"]
    assert quantity["mean_projection_performed"] is False
    assert quantity["radial_inverse_performed"] is False
    assert quantity["correction_velocity_constructed"] is False
    assert truth["strict_inner_only"] is True
    assert truth["global_leading_join_included"] is False
    assert truth["matched_pressure_included"] is False
    assert truth["restricted_forcing_included"] is False
    assert truth["same_protocol_st006_comparison_valid"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["paper_exact"] is False
    assert truth["pde_validated"] is False
