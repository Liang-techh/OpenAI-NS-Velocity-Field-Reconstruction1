from __future__ import annotations

import inspect
import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_oscillatory_self_advection_energy_flux import (
    GRADIENT_FD6_STEP,
    KokunoOscillatorySelfAdvectionEnergyFlux,
    public_contract,
)
from openai_ns_reconstruction.kokuno_public_oscillatory_self_advection import (
    evaluate_oscillatory_self_advection_fd6,
)
from openai_ns_reconstruction.kokuno_public_z_pullback_velocity import velocity_osc


def _points() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    radius = np.asarray([0.47, 0.59, 0.71, 0.83, 0.95, 1.07], dtype=float)
    theta = np.asarray([0.31, 0.89, 1.47, 2.11, -2.53, -1.37], dtype=float)
    z = np.asarray([-0.91, -0.57, -0.18, 0.24, 0.61, 0.97], dtype=float)
    t = np.asarray([0.43, 0.46, 0.49, 0.52, 0.55, 0.58], dtype=float)
    return radius * np.cos(theta), radius * np.sin(theta), z, t


def test_public_contract_has_no_tuning_or_agent3_escape_hatch() -> None:
    contract = public_contract()
    assert contract["forbidden_evaluate_inputs_present"] == []
    assert contract["public_velocity_only_independent_reaudit_required"] is True
    assert contract["production_energy_flux_divergence_independently_differentiated"] is False
    assert contract["mean_projection_performed"] is False
    assert contract["correction_velocity_constructed"] is False
    assert contract["whole_domain_energy_balance_assessed"] is False
    assert contract["complete_ns_residual"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False

    signature = inspect.signature(KokunoOscillatorySelfAdvectionEnergyFlux.evaluate)
    assert list(signature.parameters) == ["self", "x", "y", "z", "t"]


def test_energy_flux_result_replays_frozen_velocity_and_self_advection() -> None:
    x, y, z, t = _points()
    diagnostic = KokunoOscillatorySelfAdvectionEnergyFlux()
    result = diagnostic.evaluate(x, y, z, t)
    inherited = evaluate_oscillatory_self_advection_fd6(
        x, y, z, t, spatial_step=GRADIENT_FD6_STEP
    )

    assert np.array_equal(result.velocity, np.asarray(velocity_osc(x, y, z, t), dtype=float))
    assert np.array_equal(
        result.self_advection, np.asarray(inherited["self_advection"], dtype=float)
    )
    assert np.array_equal(result.divergence, np.asarray(inherited["divergence"], dtype=float))
    assert result.inherited_spatial_step == GRADIENT_FD6_STEP
    assert np.all(np.isfinite(result.kinetic_energy_density))
    assert np.all(result.kinetic_energy_density >= 0.0)
    assert float(np.sqrt(np.mean(result.self_advection_power_density**2))) > 1.0e-12


def test_local_conservative_identity_is_algebraically_closed() -> None:
    x, y, z, t = _points()
    result = KokunoOscillatorySelfAdvectionEnergyFlux().evaluate(x, y, z, t)

    expected_energy = 0.5 * np.sum(result.velocity * result.velocity, axis=-1)
    expected_power = np.sum(result.velocity * result.self_advection, axis=-1)
    expected_remainder = expected_energy * result.divergence

    assert np.array_equal(result.kinetic_energy_density, expected_energy)
    assert np.array_equal(result.self_advection_power_density, expected_power)
    assert np.array_equal(result.compressibility_remainder, expected_remainder)
    assert np.array_equal(
        result.conservative_flux_divergence_equivalent,
        expected_power + expected_remainder,
    )

    summary = result.summary_payload()
    assert summary["kinetic_energy_density_rms"] > 0.0
    assert summary["self_advection_power_rms"] > 0.0
    assert summary["compressibility_remainder_rms"] >= 0.0
    assert np.isfinite(summary["compressibility_to_power_rms_ratio"])


def test_manifest_round_trip_and_semantic_tamper_fail_closed(tmp_path) -> None:
    diagnostic = KokunoOscillatorySelfAdvectionEnergyFlux()
    path = tmp_path / "energy_flux.json"
    manifest = diagnostic.save_manifest(path)
    loaded = KokunoOscillatorySelfAdvectionEnergyFlux.load_manifest(path)
    assert loaded.energy_flux_sha256 == diagnostic.energy_flux_sha256
    assert manifest["energy_flux_sha256"] == diagnostic.energy_flux_sha256

    raw = json.loads(path.read_text(encoding="utf-8"))
    raw["payload"]["scientific_boundary"]["paper_exact"] = True
    path.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="checksum"):
        KokunoOscillatorySelfAdvectionEnergyFlux.load_manifest(path)


def test_semantic_contract_keeps_source_and_repository_realization_separate() -> None:
    payload = KokunoOscillatorySelfAdvectionEnergyFlux().semantic_payload()
    assert payload["source_provenance"]["corrected_reader_date"] == "2026-09-09"
    assert payload["quantity_contract"]["production_derivative"] == "inherits centered Cartesian FD6"
    assert payload["scientific_boundary"]["oscillatory_velocity_changed"] is False
    assert payload["scientific_boundary"]["new_oscillatory_parameters_added"] is False
    assert payload["scientific_boundary"]["whole_domain_energy_balance_assessed"] is False
    assert payload["scientific_boundary"]["paper_exact"] is False
    assert payload["scientific_boundary"]["pde_validated"] is False
