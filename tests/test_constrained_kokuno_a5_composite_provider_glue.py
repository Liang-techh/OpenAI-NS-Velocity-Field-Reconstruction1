from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a5_composite_provider_glue import (
    FINAL_DIVERGENCE_GATE,
    FINAL_MOMENTUM_GATE,
    assemble_composite_provider_set,
    deterministic_receipt,
    stage_from_bundle_evaluator,
    stage_from_vector_evaluators,
    truth_boundary,
)
from openai_ns_reconstruction.kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
)


def _zero_pressure(identity: CycleIdentity) -> ScalarFieldProvider:
    return ScalarFieldProvider(identity, "test:p=0", lambda x, y, z, t: 0.0)


def _zero_force(identity: CycleIdentity) -> RestrictedForcingProvider:
    return RestrictedForcingProvider(
        identity,
        "test:f=0",
        "fixed-zero-force-independent-of-residual",
        lambda x, y, z, t: (0.0, 0.0, 0.0),
    )


def test_deterministic_receipt_is_manufactured_and_fail_closed() -> None:
    receipt = deterministic_receipt()
    q = receipt["manufactured_regression"]
    assert q["assembly_max_abs_error"] == 0.0
    assert q["same_cycle_defect_max_abs_error"] <= 2.0e-11
    assert q["held_in_vector_rms"] > 0.0
    assert q["held_out_vector_rms"] > 0.0
    assert q["restricted_forcing_receipt"] == "manufactured-zero-force-fixed-before-residual"

    boundary = receipt["truth_boundary"]
    assert boundary["typed_composite_velocity_provider_executable"] is True
    assert boundary["leading_plus_oscillatory_assembly_executable"] is True
    assert boundary["optional_correction_slot_executable"] is True
    assert boundary["pressure_forwarded_without_fit"] is True
    assert boundary["restricted_forcing_forwarded_without_redefinition"] is True
    assert boundary["real_agent1_global_leading_handoff_available_here"] is False
    assert boundary["real_agent1_matched_pressure_handoff_available_here"] is False
    assert boundary["real_full_candidate_instantiated"] is False
    assert boundary["correction_ready"] is False
    assert boundary["velocity_export_ready"] is False
    assert boundary["pde_validated"] is False
    assert boundary["final_normalized_momentum_gate"] == FINAL_MOMENTUM_GATE == 1.0e-3
    assert boundary["final_normalized_divergence_gate"] == FINAL_DIVERGENCE_GATE == 1.0e-5
    assert "manufactured" in receipt["scientific_scope"]


def test_composite_sum_and_optional_correction_are_exact() -> None:
    identity = CycleIdentity("sum-test", 2, "state")
    leading = stage_from_vector_evaluators(
        identity=identity,
        stage="leading",
        source_ref="leading",
        velocity_evaluator=lambda x, y, z, t: (x, y, z),
        velocity_dt_evaluator=lambda x, y, z, t: (1.0, 0.0, 0.0),
    )
    oscillatory = stage_from_bundle_evaluator(
        identity=identity,
        stage="oscillatory",
        source_ref="osc",
        bundle_evaluator=lambda x, y, z, t: {
            "velocity": (t, -t, 2.0 * t),
            "velocity_dt": (1.0, -1.0, 2.0),
        },
    )
    correction = stage_from_vector_evaluators(
        identity=identity,
        stage="correction",
        source_ref="corr",
        velocity_evaluator=lambda x, y, z, t: (0.5, -0.25, 0.125),
        velocity_dt_evaluator=lambda x, y, z, t: (0.0, 0.0, 0.0),
    )
    pressure = _zero_pressure(identity)
    force = _zero_force(identity)
    providers = assemble_composite_provider_set(
        identity=identity,
        leading=leading,
        oscillatory=oscillatory,
        correction=correction,
        pressure=pressure,
        restricted_forcing=force,
    )

    point = (0.2, -0.3, 0.4, 0.6)
    np.testing.assert_allclose(
        providers.velocity.evaluator(*point),
        np.asarray((1.3, -1.15, 1.725)),
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        providers.velocity_dt.evaluator(*point),
        np.asarray((2.0, -1.0, 2.0)),
        rtol=0.0,
        atol=0.0,
    )
    assert providers.pressure is pressure
    assert providers.restricted_forcing is force
    assert providers.correction_included is True
    assert providers.stage_source_refs == ("leading", "osc", "corr")


def test_identity_and_stage_role_mismatches_fail_closed() -> None:
    identity = CycleIdentity("cycle", 0, "a")
    other = CycleIdentity("cycle", 0, "b")
    leading = stage_from_vector_evaluators(
        identity=identity,
        stage="leading",
        source_ref="leading",
        velocity_evaluator=lambda x, y, z, t: (0.0, 0.0, 0.0),
        velocity_dt_evaluator=lambda x, y, z, t: (0.0, 0.0, 0.0),
    )
    bad_osc = stage_from_vector_evaluators(
        identity=other,
        stage="oscillatory",
        source_ref="osc",
        velocity_evaluator=lambda x, y, z, t: (0.0, 0.0, 0.0),
        velocity_dt_evaluator=lambda x, y, z, t: (0.0, 0.0, 0.0),
    )
    with pytest.raises(ValueError, match="oscillatory does not share"):
        assemble_composite_provider_set(
            identity=identity,
            leading=leading,
            oscillatory=bad_osc,
            pressure=_zero_pressure(identity),
            restricted_forcing=_zero_force(identity),
        )

    with pytest.raises(ValueError, match="stage must"):
        stage_from_vector_evaluators(
            identity=identity,
            stage="pressure",
            source_ref="bad",
            velocity_evaluator=lambda x, y, z, t: (0.0, 0.0, 0.0),
            velocity_dt_evaluator=lambda x, y, z, t: (0.0, 0.0, 0.0),
        )


def test_bundle_adapter_rejects_missing_or_malformed_fields() -> None:
    identity = CycleIdentity("bundle", 0, "fixed")
    missing_dt = stage_from_bundle_evaluator(
        identity=identity,
        stage="oscillatory",
        source_ref="missing-dt",
        bundle_evaluator=lambda x, y, z, t: {"velocity": (1.0, 2.0, 3.0)},
    )
    with pytest.raises(ValueError, match="omitted velocity_dt"):
        missing_dt.velocity_dt_evaluator(0.0, 0.0, 0.0, 0.0)

    malformed = stage_from_bundle_evaluator(
        identity=identity,
        stage="oscillatory",
        source_ref="malformed",
        bundle_evaluator=lambda x, y, z, t: {
            "velocity": (1.0, 2.0),
            "velocity_dt": (0.0, 0.0, 0.0),
        },
    )
    with pytest.raises(ValueError, match="finite Cartesian 3-vector"):
        malformed.velocity_evaluator(0.0, 0.0, 0.0, 0.0)


def test_public_api_has_no_residual_or_target_injection() -> None:
    forbidden = {"residual", "defect", "stress", "target", "gain", "normalized_score"}
    for func in (
        assemble_composite_provider_set,
        stage_from_bundle_evaluator,
        stage_from_vector_evaluators,
    ):
        assert forbidden.isdisjoint(inspect.signature(func).parameters)

    boundary = truth_boundary()
    assert boundary["caller_supplied_residual_allowed"] is False
    assert boundary["caller_supplied_defect_allowed"] is False
    assert boundary["caller_supplied_stress_allowed"] is False
    assert boundary["caller_supplied_target_allowed"] is False
