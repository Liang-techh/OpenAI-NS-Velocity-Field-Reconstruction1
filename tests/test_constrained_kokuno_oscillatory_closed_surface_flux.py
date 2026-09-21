from __future__ import annotations

import inspect
import json

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_oscillatory_closed_surface_flux as fluxmod
from openai_ns_reconstruction.kokuno_oscillatory_closed_surface_flux import (
    FACE_FLUX_DENOMINATOR_FLOOR,
    KokunoOscillatoryClosedSurfaceFlux,
    public_contract,
)


def _analytic_divergence_free_velocity(x, y, z, t):
    """Linear divergence-free field with nonzero cylindrical face fluxes."""
    x, y, z, t = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    del t
    return np.stack((x, y, -2.0 * z), axis=-1)


def test_closed_surface_quadrature_recovers_nontrivial_divergence_free_flux(monkeypatch):
    monkeypatch.setattr(fluxmod, "velocity_osc", _analytic_divergence_free_velocity)
    result = KokunoOscillatoryClosedSurfaceFlux().materialize()

    assert len(result.control_volumes) == 2
    assert len(result.levels) == 3
    for level in result.levels:
        assert level["min_absolute_face_flux_sum"] > FACE_FLUX_DENOMINATOR_FLOOR
        assert level["max_relative_closure"] < 2.0e-13
        for case in level["cases"]:
            assert case["absolute_face_flux_sum"] > FACE_FLUX_DENOMINATOR_FLOOR
            assert case["max_absolute_face_flux"] > 1.0e-8
            assert abs(case["radial_inner"]) > 1.0e-8
            assert abs(case["radial_outer"]) > 1.0e-8
            assert abs(case["axial_lower"]) > 1.0e-8
            assert abs(case["axial_upper"]) > 1.0e-8
            assert case["relative_closure"] < 2.0e-13

    shifted = result.shifted_fine
    assert shifted["min_absolute_face_flux_sum"] > FACE_FLUX_DENOMINATOR_FLOOR
    assert shifted["max_relative_closure"] < 2.0e-13


def test_control_volumes_are_strictly_inside_public_support():
    result = KokunoOscillatoryClosedSurfaceFlux().materialize()
    field = fluxmod.default_field()
    for volume in result.control_volumes:
        assert float(field.radial_inner) < volume["r_inner"] < volume["r_outer"] < float(field.radial_outer)
        assert float(field.axial_lower) < volume["z_lower"] < volume["z_upper"] < float(field.axial_upper)


def test_nonfinite_public_velocity_fails_closed(monkeypatch):
    def bad_velocity(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(x, y, z, t)
        return np.full(x.shape + (3,), np.nan, dtype=float)

    monkeypatch.setattr(fluxmod, "velocity_osc", bad_velocity)
    with pytest.raises(RuntimeError, match="invalid Cartesian array"):
        KokunoOscillatoryClosedSurfaceFlux().materialize()


def test_manifest_round_trip_and_semantic_tamper_rejection(tmp_path):
    obj = KokunoOscillatoryClosedSurfaceFlux()
    path = tmp_path / "closed_surface_flux.json"
    manifest = obj.save_manifest(path)
    rebound = KokunoOscillatoryClosedSurfaceFlux.load_manifest(path)
    assert rebound.closed_surface_flux_sha256 == obj.closed_surface_flux_sha256
    assert manifest["closed_surface_flux_sha256"] == obj.closed_surface_flux_sha256

    raw = json.loads(path.read_text(encoding="utf-8"))
    raw["payload"]["fine_theta_shift_cells"] = 0.5
    # Even an attacker who recomputes the outer checksum cannot alter frozen
    # scientific semantics: load must rebuild the canonical contract.
    raw["closed_surface_flux_sha256"] = fluxmod._sha256(raw["payload"])
    path.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="semantics mismatch"):
        KokunoOscillatoryClosedSurfaceFlux.load_manifest(path)


def test_public_contract_has_no_scientific_tuning_or_lane_promotion():
    contract = public_contract()
    assert contract["forbidden_materialize_inputs_present"] == []
    assert contract["generic_interior_closed_surface_flux_assessed"] is True
    assert contract["individual_face_fluxes_required_nontrivial"] is True
    assert contract["uses_only_public_velocity_values"] is True
    assert contract["volume_divergence_integral_assessed"] is False
    assert contract["whole_domain_divergence_volume_l2_assessed"] is False
    assert contract["mean_projection_performed"] is False
    assert contract["correction_velocity_constructed"] is False
    assert contract["complete_ns_residual"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False

    params = inspect.signature(KokunoOscillatoryClosedSurfaceFlux.materialize).parameters
    assert list(params) == ["self"]
