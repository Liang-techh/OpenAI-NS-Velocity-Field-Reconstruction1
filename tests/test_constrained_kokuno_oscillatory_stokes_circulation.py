from __future__ import annotations

import inspect
import json

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_oscillatory_stokes_circulation as stokesmod
from openai_ns_reconstruction.kokuno_oscillatory_stokes_circulation import (
    KokunoOscillatoryStokesCirculation,
    STOKES_SIGNAL_FLOOR,
    public_contract,
)


def _solid_rotation_velocity(x, y, z, t):
    """Divergence-free linear field with curl=(1.2,-0.8,0.6)."""
    x, y, z, t = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    del t
    return np.stack(
        (
            -0.4 * z - 0.3 * y,
            0.3 * x - 0.6 * z,
            0.6 * y + 0.4 * x,
        ),
        axis=-1,
    )


def test_stokes_quadrature_and_fd4_recover_manufactured_curl(monkeypatch):
    monkeypatch.setattr(stokesmod, "velocity_osc", _solid_rotation_velocity)
    result = KokunoOscillatoryStokesCirculation().materialize()

    expected_curl = np.asarray((1.2, -0.8, 0.6), dtype=float)
    assert len(result.surfaces) == 3
    assert len(result.levels) == 3
    for level in result.levels:
        assert level["max_relative_mismatch"] < 2.0e-12
        assert level["min_signal_scale"] > STOKES_SIGNAL_FLOOR
        for case in level["cases"]:
            spec = result.surfaces[case["surface_index"]]
            area = 4.0 * spec["halfwidth_a"] * spec["halfwidth_b"]
            expected = area * expected_curl[spec["normal_axis"]]
            assert case["boundary_circulation"] == pytest.approx(expected, rel=0.0, abs=2.0e-13)
            assert case["surface_curl_flux"] == pytest.approx(expected, rel=0.0, abs=2.0e-12)

    assert result.shifted_fine["max_relative_mismatch"] < 2.0e-12
    assert result.shifted_fine["min_signal_scale"] > STOKES_SIGNAL_FLOOR


def test_registered_surfaces_and_fd_collars_remain_strictly_inside_support():
    field = stokesmod.default_field()
    for raw in stokesmod.SURFACES:
        spec = stokesmod._surface_spec(raw)
        stokesmod._assert_surface_inside_support(spec, max(stokesmod.FD_STEPS))
        corners = stokesmod._surface_corners(spec)
        radius = np.hypot(corners[:, 0], corners[:, 1])
        assert np.min(radius) > float(field.radial_inner)
        assert np.max(radius) < float(field.radial_outer)
        assert np.min(corners[:, 2]) > float(field.axial_lower)
        assert np.max(corners[:, 2]) < float(field.axial_upper)


def test_nonfinite_public_velocity_fails_closed(monkeypatch):
    def bad_velocity(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(x, y, z, t)
        return np.full(x.shape + (3,), np.nan, dtype=float)

    monkeypatch.setattr(stokesmod, "velocity_osc", bad_velocity)
    with pytest.raises(RuntimeError, match="invalid Cartesian array"):
        KokunoOscillatoryStokesCirculation().materialize()


def test_manifest_round_trip_and_semantic_tamper_rejection(tmp_path):
    obj = KokunoOscillatoryStokesCirculation()
    path = tmp_path / "stokes.json"
    manifest = obj.save_manifest(path)
    rebound = KokunoOscillatoryStokesCirculation.load_manifest(path)
    assert rebound.stokes_circulation_sha256 == obj.stokes_circulation_sha256
    assert manifest["stokes_circulation_sha256"] == obj.stokes_circulation_sha256

    raw = json.loads(path.read_text(encoding="utf-8"))
    raw["payload"]["fd_steps"][-1] = 0.0005
    raw["stokes_circulation_sha256"] = stokesmod._sha256(raw["payload"])
    path.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="semantics mismatch"):
        KokunoOscillatoryStokesCirculation.load_manifest(path)


def test_public_contract_has_no_tuning_or_agent3_lane_promotion():
    contract = public_contract()
    assert contract["forbidden_materialize_inputs_present"] == []
    assert contract["uses_only_public_velocity_values"] is True
    assert contract["production_jacobian_used"] is False
    assert contract["production_divergence_used"] is False
    assert contract["production_vorticity_used"] is False
    assert contract["stokes_circulation_consistency_assessed"] is True
    assert contract["whole_domain_vorticity_assessed"] is False
    assert contract["mean_projection_performed"] is False
    assert contract["radial_inverse_performed"] is False
    assert contract["correction_velocity_constructed"] is False
    assert contract["complete_ns_residual"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False

    params = inspect.signature(KokunoOscillatoryStokesCirculation.materialize).parameters
    assert list(params) == ["self"]
