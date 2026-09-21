from __future__ import annotations

import hashlib
import inspect
import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_oscillatory_support_boundary_flux import (
    BOUNDARY_RADIUS_COUNT,
    BOUNDARY_THETA_COUNT,
    BOUNDARY_TIMES,
    BOUNDARY_Z_COUNT,
    INWARD_COLLAR_FRACTIONS,
    KokunoOscillatorySupportBoundaryFlux,
    public_contract,
)


def _canonical_json(payload: dict[str, object]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def test_support_boundary_trace_and_energy_flux_are_exactly_zero() -> None:
    diagnostic = KokunoOscillatorySupportBoundaryFlux()
    result = diagnostic.materialize()

    expected = 2 * len(BOUNDARY_TIMES) * BOUNDARY_THETA_COUNT * (
        BOUNDARY_Z_COUNT + BOUNDARY_RADIUS_COUNT
    )
    assert result.boundary_sample_count == expected
    assert result.outside_sample_count == expected
    assert result.radial_support == pytest.approx((0.15, 1.35), abs=1e-15)
    assert result.axial_support == pytest.approx((-2.0, 2.0), abs=1e-15)
    assert result.radial_support[0] > 0.0

    assert result.pointwise_boundary_trace_zero
    assert result.pointwise_outside_zero
    assert result.boundary_velocity_max_abs == 0.0
    assert result.outside_velocity_max_abs == 0.0
    assert result.boundary_normal_energy_flux_max_abs == 0.0
    assert result.support_boundary_energy_flux_integral == 0.0


def test_manifest_round_trip_and_semantic_tamper_fail_closed(tmp_path) -> None:
    diagnostic = KokunoOscillatorySupportBoundaryFlux()
    path = tmp_path / "support_flux.json"
    saved = diagnostic.save_manifest(path)
    rebound = KokunoOscillatorySupportBoundaryFlux.load_manifest(path)
    assert rebound.support_flux_sha256 == diagnostic.support_flux_sha256
    assert rebound.semantic_payload() == diagnostic.semantic_payload()

    raw = json.loads(path.read_text(encoding="utf-8"))
    raw["payload"]["scientific_boundary"]["pde_validated"] = True
    raw["support_flux_sha256"] = hashlib.sha256(
        _canonical_json(raw["payload"]).encode("utf-8")
    ).hexdigest()
    path.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="semantics mismatch"):
        KokunoOscillatorySupportBoundaryFlux.load_manifest(path)

    # The independent collar protocol is frozen into artifact identity and cannot
    # be post-hoc retuned by the public materializer.
    assert diagnostic.semantic_payload()["frozen_diagnostic_protocol"][
        "independent_inward_collar_fractions"
    ] == list(INWARD_COLLAR_FRACTIONS)


def test_public_contract_has_no_tuning_escape_hatch() -> None:
    contract = public_contract()
    assert contract["forbidden_materialize_inputs_present"] == []
    assert list(inspect.signature(KokunoOscillatorySupportBoundaryFlux.materialize).parameters) == ["self"]
    assert contract["support_boundary_trace_assessed"] is True
    assert contract["support_boundary_energy_flux_assessed"] is True
    assert contract["whole_domain_energy_balance_assessed"] is False
    assert contract["mean_projection_performed"] is False
    assert contract["correction_velocity_constructed"] is False
    assert contract["complete_ns_residual"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False


def test_summary_payload_is_finite_and_truthful() -> None:
    summary = KokunoOscillatorySupportBoundaryFlux().materialize().summary_payload()
    for key in (
        "boundary_velocity_max_abs",
        "outside_velocity_max_abs",
        "boundary_normal_energy_flux_max_abs",
        "support_boundary_energy_flux_integral",
    ):
        assert np.isfinite(summary[key])
        assert summary[key] == 0.0
    assert summary["pointwise_boundary_trace_zero"] is True
    assert summary["pointwise_outside_zero"] is True
