from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
    deterministic_receipt,
    evaluate_disjoint_heldin_heldout,
    evaluate_same_cycle_momentum_defect,
)


def _providers(identity: CycleIdentity):
    velocity = VectorFieldProvider(
        identity, "u", lambda x, y, z, t: (y, 0.0, 0.0)
    )
    velocity_dt = VectorFieldProvider(
        identity, "u_t", lambda x, y, z, t: (0.0, 0.0, 0.0)
    )
    pressure = ScalarFieldProvider(
        identity, "p", lambda x, y, z, t: 2.0 * x + 3.0 * z
    )
    forcing = RestrictedForcingProvider(
        identity,
        "f",
        "registered zero forcing; no residual dependence",
        lambda x, y, z, t: (0.0, 0.0, 0.0),
    )
    return velocity, velocity_dt, pressure, forcing


def test_computes_nonzero_actual_momentum_defect_without_residual_input() -> None:
    identity = CycleIdentity("cycle-a", 0, "state-0")
    providers = _providers(identity)
    result = evaluate_same_cycle_momentum_defect(
        *providers,
        [(0.17, -0.23, 0.31, 0.41), (-0.29, 0.12, -0.19, 0.57)],
        viscosity=0.01,
        spatial_step=1.0e-3,
    )
    expected = np.asarray([[2.0, 0.0, 3.0], [2.0, 0.0, 3.0]])
    np.testing.assert_allclose(
        result.residual_values, expected, rtol=0.0, atol=2.0e-11
    )
    assert result.vector_rms == pytest.approx(
        np.sqrt(13.0), rel=0.0, abs=2.0e-11
    )
    assert result.vector_rms > 0.0


def test_rejects_mixed_cycle_providers() -> None:
    identity = CycleIdentity("cycle-a", 0, "state-0")
    velocity, velocity_dt, pressure, forcing = _providers(identity)
    pressure = ScalarFieldProvider(
        CycleIdentity("cycle-b", 0, "state-0"),
        "wrong-pressure",
        pressure.evaluator,
    )
    with pytest.raises(ValueError, match="same-cycle identity"):
        evaluate_same_cycle_momentum_defect(
            velocity,
            velocity_dt,
            pressure,
            forcing,
            [(0.1, 0.2, 0.3, 0.4)],
            viscosity=0.01,
            spatial_step=1.0e-3,
        )


def test_rejects_heldin_heldout_reuse() -> None:
    identity = CycleIdentity("cycle-a", 0, "state-0")
    providers = _providers(identity)
    shared = (0.1, 0.2, 0.3, 0.4)
    with pytest.raises(ValueError, match="must be disjoint"):
        evaluate_disjoint_heldin_heldout(
            *providers,
            held_in_points=[shared, (-0.2, 0.1, -0.3, 0.5)],
            held_out_points=[shared, (0.3, -0.1, 0.2, 0.6)],
            viscosity=0.01,
            spatial_step=1.0e-3,
        )


def test_public_evaluator_has_no_surrogate_injection_parameters() -> None:
    names = set(inspect.signature(evaluate_same_cycle_momentum_defect).parameters)
    forbidden = {
        "residual",
        "defect",
        "stress",
        "target",
        "pressure_gradient",
        "gain",
        "normalized_residual",
    }
    assert names.isdisjoint(forbidden)
    assert "restricted_forcing" in names


def test_receipt_keeps_full_candidate_and_pde_claims_false() -> None:
    receipt = deterministic_receipt()
    regression = receipt["formula_regression"]
    assert regression["maximum_absolute_error"] <= 2.0e-11
    assert regression["held_in"]["vector_rms"] == pytest.approx(
        np.sqrt(13.0), abs=2.0e-11
    )
    assert regression["held_out"]["vector_rms"] == pytest.approx(
        np.sqrt(13.0), abs=2.0e-11
    )
    boundary = receipt["truth_boundary"]
    assert boundary["actual_momentum_defect_formula_executable"] is True
    assert boundary["caller_supplied_residual_allowed"] is False
    assert boundary["heldin_heldout_disjointness_enforced"] is True
    assert boundary["restricted_forcing_semantics_independently_validated_here"] is False
    assert boundary["real_full_candidate_defect_consumed"] is False
    assert boundary["finite_correction_cycle_run"] is False
    assert boundary["heldout_ns_momentum_residual_assessed"] is False
    assert boundary["pde_validated"] is False
