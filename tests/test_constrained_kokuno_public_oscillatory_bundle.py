from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_oscillatory_bundle import (
    IDENTITY,
    PARENT_AGENT2_HEAD,
    SCHEMA,
    UPSTREAM,
    evaluate_oscillatory_bundle,
    verification_receipt,
)
from openai_ns_reconstruction.kokuno_public_oscillatory_time_derivative import (
    velocity_osc_dt,
)
from openai_ns_reconstruction.kokuno_public_oscillatory_vector_potential import (
    vector_potential_osc,
)
from openai_ns_reconstruction.kokuno_public_oscillatory_vector_potential_time_derivative import (
    vector_potential_osc_dt,
)
from openai_ns_reconstruction.kokuno_public_z_pullback_velocity import velocity_osc


def _cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray((0.53, -0.44, 0.88, -0.76), dtype=float)
    y = np.asarray((0.17, 0.61, -0.25, -0.50), dtype=float)
    z = np.asarray((-0.41, 0.39, -0.33, 0.37), dtype=float)
    t = np.asarray((0.37, 0.49, 0.57, 0.63), dtype=float)
    return x, y, z, t


def test_bundle_replays_all_four_frozen_public_interfaces_exactly() -> None:
    x, y, z, t = _cloud()
    bundle = evaluate_oscillatory_bundle(x, y, z, t)

    assert np.array_equal(bundle["vector_potential"], vector_potential_osc(x, y, z, t))
    assert np.array_equal(bundle["velocity"], velocity_osc(x, y, z, t))
    assert np.array_equal(
        bundle["vector_potential_dt"], vector_potential_osc_dt(x, y, z, t)
    )
    assert np.array_equal(bundle["velocity_dt"], velocity_osc_dt(x, y, z, t))
    assert bundle["vector_potential"].shape == (4, 3)
    assert bundle["support_mask"].shape == (4,)
    assert np.all(bundle["support_mask"])


def test_bundle_preserves_broadcast_and_registered_support_zero() -> None:
    bundle = evaluate_oscillatory_bundle(
        np.asarray([[0.0], [1.85], [0.0]]),
        np.zeros((3, 1)),
        np.asarray([[0.0, 0.4, 2.15]]),
        0.5,
    )
    assert bundle["velocity"].shape == (3, 3, 3)
    assert bundle["support_mask"].shape == (3, 3)

    # Every broadcast row is outside through the axis or radial exterior; the
    # third column additionally exercises the axial exterior.  All four frozen
    # public fields must remain exact zero.
    assert not np.any(bundle["support_mask"])
    for key in ("vector_potential", "velocity", "vector_potential_dt", "velocity_dt"):
        assert np.max(np.abs(bundle[key])) == 0.0


def test_bundle_identity_is_fail_closed_and_truth_boundary_stays_narrow() -> None:
    payload = IDENTITY.payload()
    digest = IDENTITY.sha256()

    assert payload["schema"] == SCHEMA
    assert payload["parent_agent2_head"] == PARENT_AGENT2_HEAD
    assert len(digest) == 64
    int(digest, 16)
    assert UPSTREAM["velocity"]["pr"] == 561
    assert UPSTREAM["velocity_dt"]["pr"] == 579
    assert UPSTREAM["vector_potential"]["pr"] == 616
    assert UPSTREAM["vector_potential_dt"]["pr"] == 625
    assert "independently_admitted" in UPSTREAM["velocity"]["status"]
    assert "independently_admitted" in UPSTREAM["velocity_dt"]["status"]
    assert "pending_independent_agent4_audit" in UPSTREAM["vector_potential"]["status"]
    assert "pending_independent_agent4_audit" in UPSTREAM["vector_potential_dt"]["status"]

    result = evaluate_oscillatory_bundle(0.62, 0.13, 0.21, 0.5)
    truth = result["truth_boundary"]
    assert truth["velocity_candidate_changed"] is False
    assert truth["vector_potential_candidate_changed"] is False
    assert truth["source_formula_changed"] is False
    assert truth["integration_adapter_only"] is True
    assert truth["independent_agent4_vector_potential_audit_required"] is True
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_bundle_api_does_not_accept_residual_pressure_force_or_target_inputs() -> None:
    names = set(inspect.signature(evaluate_oscillatory_bundle).parameters)
    assert names == {"x", "y", "z", "t"}
    forbidden = {
        "residual",
        "defect",
        "stress",
        "target",
        "pressure",
        "forcing",
        "gain",
        "correction",
    }
    assert names.isdisjoint(forbidden)


def test_verification_receipt_is_exact_replay_only() -> None:
    receipt = verification_receipt()
    assert receipt["schema"] == SCHEMA
    assert receipt["parent_agent2_head"] == PARENT_AGENT2_HEAD
    assert receipt["failed_guards"] == []
    assert receipt["support_interior_count"] > 0
    assert receipt["support_exterior_absolute_max"] == 0.0
    assert all(value == 0.0 for value in receipt["exact_public_api_replay_max_abs"].values())
    assert "no new curl/divergence/NS residual admission" in receipt["scientific_scope"]


def test_bundle_rejects_nonfinite_input_through_frozen_provider_contract() -> None:
    with pytest.raises(ValueError):
        evaluate_oscillatory_bundle(np.nan, 0.0, 0.0, 0.5)
