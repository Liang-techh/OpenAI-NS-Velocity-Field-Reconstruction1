from __future__ import annotations

import inspect
from types import SimpleNamespace

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_current_rf40_power_law_nonlinear_radial_force as mod
from openai_ns_reconstruction.kokuno_strict_inner_transport_radial_stress import (
    StrictInnerTransportRadialGeometry,
)


def _geometry() -> StrictInnerTransportRadialGeometry:
    return StrictInnerTransportRadialGeometry(
        time=0.39,
        axial_z=-0.08,
        radii=(0.04, 0.12, 0.20, 0.28, 0.36, 0.44),
        bump_center=0.30,
        bump_halfwidth=0.10,
    )


def test_public_signature_and_truth_boundary_are_fail_closed() -> None:
    signature = inspect.signature(
        mod.materialize_current_rf40_power_law_nonlinear_radial_force
    )
    assert tuple(signature.parameters) == ("backend", "geometry")
    boundary = mod.truth_boundary()
    assert boundary["parent_agent3_pr"] == 1037
    assert boundary["parent_agent3_head"] == (
        "9a5cdcdf78b7862d5ebe171efb9ca56d53664612"
    )
    assert boundary["parent_agent3_source_blob"] == (
        "44dc6d61bc965cef08b119aa5df3026175abd5de"
    )
    assert boundary["source_radial_force_formula"] == "(div T)_r = partial_z sigma_1"
    assert boundary["current_RF40_power_law_axial_e1_stress_consumed"] is True
    assert (
        boundary[
            "current_RF40_power_law_radial_force_from_dz_sigma1_materialized"
        ]
        is True
    )
    assert boundary["independent_third_radial_moment_inverse_introduced"] is False
    assert boundary["recorded_radial_mean_equated_to_radial_force"] is False
    assert boundary["forbidden_public_parameters_absent"] is True
    assert boundary["partial_domain_through_RF40_power_law_only"] is True
    assert boundary["velocity_after_RF40_power_law_materialized"] is False
    assert boundary["newer_agent1_modulated_lineage_consumed"] is False
    assert boundary["complete_ns_defect"] is False
    assert (
        boundary[
            "scoped_current_RF40_power_law_radial_force_authorized_as_correction_target"
        ]
        is False
    )
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["pde_validated"] is False
    assert boundary["final_normalized_momentum_gate"] == pytest.approx(1e-3)
    assert boundary["final_normalized_divergence_gate"] == pytest.approx(1e-5)


def test_fixed_z_ladder_recovers_known_derivative_and_piece_closure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class FakeBackend:
        pass

    geometry = _geometry()
    radii = np.asarray(geometry.radii, dtype=float)

    monkeypatch.setattr(mod, "ExactCurrentRF40PowerLawNonlinearBackend", FakeBackend)
    monkeypatch.setattr(
        mod,
        "_authenticate_parent_source",
        lambda: mod.PARENT_AGENT3_SOURCE_BLOB,
    )

    def fake_stress(_backend: object, shifted: StrictInnerTransportRadialGeometry):
        return SimpleNamespace(z=float(shifted.axial_z))

    def fake_pieces(witness: SimpleNamespace) -> dict[str, np.ndarray]:
        z = float(witness.z)
        quadratic = radii + 1.5 * z + 0.2 * z * z
        mixed = 0.3 * radii - 0.4 * z + 0.05 * z * z * z
        return {
            "quadratic": quadratic,
            "mixed": mixed,
            "aggregate": quadratic + mixed,
        }

    monkeypatch.setattr(
        mod,
        "materialize_current_rf40_power_law_nonlinear_radial_stress",
        fake_stress,
    )
    monkeypatch.setattr(mod, "_axial_stress_pieces", fake_pieces)

    witness = mod.materialize_current_rf40_power_law_nonlinear_radial_force(
        FakeBackend(), geometry
    )
    z0 = geometry.axial_z
    expected_q = np.full_like(radii, 1.5 + 0.4 * z0)
    expected_m = np.full_like(radii, -0.4 + 0.15 * z0 * z0)
    np.testing.assert_allclose(
        witness.finest_quadratic_radial_force, expected_q, rtol=0, atol=1e-11
    )
    np.testing.assert_allclose(
        witness.finest_mixed_radial_force, expected_m, rtol=0, atol=1e-11
    )
    np.testing.assert_allclose(
        witness.finest_aggregate_radial_force,
        expected_q + expected_m,
        rtol=0,
        atol=1e-11,
    )
    assert witness.z_steps == tuple(float(v) for v in mod.Z_DERIVATIVE_STEP_LADDER)
    assert max(witness.piece_closure_relative_max_levels) <= (
        mod.RADIAL_FORCE_PIECE_CLOSURE_RELATIVE_GATE
    )
    assert witness.aggregate_derivative_stability_preflight_passed is True


def test_wrong_backend_fails_before_parent_evaluation() -> None:
    with pytest.raises(TypeError, match="backend must be"):
        mod.materialize_current_rf40_power_law_nonlinear_radial_force(
            object(), _geometry()
        )


def test_parent_source_authentication_rejects_blob_drift(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    fake_source = tmp_path / "fake_parent.py"
    fake_source.write_text("# drifted parent\n", encoding="utf-8")
    monkeypatch.setattr(inspect, "getsourcefile", lambda _fn: str(fake_source))
    with pytest.raises(RuntimeError, match="source blob drifted"):
        mod._authenticate_parent_source()
