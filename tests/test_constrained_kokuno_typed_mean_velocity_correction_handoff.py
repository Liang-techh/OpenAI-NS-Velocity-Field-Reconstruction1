from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_same_cycle_defect_contract import CycleIdentity
from openai_ns_reconstruction.kokuno_typed_mean_velocity_correction_handoff import (
    CompleteCurlCorrectionAdapter,
    _analytic_regression_handoff,
    deterministic_receipt,
    materialize_typed_mean_velocity_correction_handoff,
    truth_boundary,
)


def _zero_adapter(*, identity: CycleIdentity, radii: tuple[float, ...], certified: bool, kind: str):
    return CompleteCurlCorrectionAdapter(
        identity=identity,
        radii=radii,
        producer_kind=kind,
        provenance="unit-test adapter",
        velocity_evaluator=lambda amplitude, x, y, z, t: (0.0, 0.0, 0.0),
        velocity_dt_evaluator=lambda amplitude, x, y, z, t: (0.0, 0.0, 0.0),
        source_agent2_complete_curl_certified=certified,
    )


def test_deterministic_handoff_is_typed_nontrivial_plumbing_only() -> None:
    receipt = deterministic_receipt()
    assert receipt["schema"] == "kokuno-a3-typed-velocity-correction-handoff-v1"
    q = receipt["analytic_regression"]
    assert q["forcing_is_zero_and_fixed_independently_of_residual"] is True
    assert q["adapter_is_source_agent2_certified"] is False
    assert q["adapter_is_only_typed_plumbing_probe"] is True
    assert q["correction_nontrivial"] is True
    assert q["correction_vector_rms"] > 0.0
    assert q["correction_vector_max"] >= q["correction_vector_rms"]
    assert q["correction_dt_vector_max"] == 0.0
    assert q["to_identity"]["cycle_index"] == q["from_identity"]["cycle_index"] + 1

    boundary = receipt["truth_boundary"]
    assert boundary["typed_delta_a_recomputed_internally"] is True
    assert boundary["external_complete_curl_adapter_required"] is True
    assert boundary["velocity_and_time_derivative_required_from_adapter"] is True
    assert boundary["agent3_reimplements_agent2_complete_curl"] is False
    assert boundary["source_agent2_complete_curl_adapter_bound_for_full_candidate"] is False
    assert boundary["public_velocity_correction_materialized_for_real_candidate"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["pde_validated"] is False


def test_regression_handoff_returns_finite_cycle_correction_provider() -> None:
    handoff = _analytic_regression_handoff()
    provider = handoff.correction_provider
    assert provider.from_identity == handoff.identity
    assert provider.to_identity == handoff.to_identity
    value = np.asarray(provider.velocity_evaluator(0.21, -0.17, 0.03, 0.0), dtype=float)
    value_dt = np.asarray(provider.velocity_dt_evaluator(0.21, -0.17, 0.03, 0.0), dtype=float)
    assert value.shape == (3,)
    assert value_dt.shape == (3,)
    assert np.all(np.isfinite(value))
    assert np.all(np.isfinite(value_dt))
    assert np.linalg.norm(value) > 0.0
    assert np.linalg.norm(value_dt) == 0.0


def test_source_certified_adapter_must_name_agent2_complete_curl() -> None:
    identity = CycleIdentity("adapter-certification-test", 0, "s0")
    with pytest.raises(ValueError, match="Agent 2 complete-curl provenance"):
        _zero_adapter(
            identity=identity,
            radii=(0.2, 0.3),
            certified=True,
            kind="unrelated-provider",
        )
    adapter = _zero_adapter(
        identity=identity,
        radii=(0.2, 0.3),
        certified=True,
        kind="agent2-public-complete-curl-adapter",
    )
    assert adapter.source_agent2_complete_curl_certified is True


def test_public_handoff_signature_has_no_surrogate_defect_or_target_inputs() -> None:
    parameters = set(inspect.signature(materialize_typed_mean_velocity_correction_handoff).parameters)
    forbidden = {
        "residual",
        "defect",
        "mean_source",
        "source",
        "stress",
        "requested_stress",
        "debt",
        "delta_c",
        "inverse_rhs",
        "rhs",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "delta_y",
        "delta_a",
    }
    assert parameters.isdisjoint(forbidden)
    boundary = truth_boundary()
    assert boundary["forbidden_public_parameters_absent"] is True
    assert boundary["caller_supplied_residual_allowed"] is False
    assert boundary["caller_supplied_defect_allowed"] is False
    assert boundary["caller_supplied_delta_y_allowed"] is False
    assert boundary["caller_supplied_delta_a_allowed"] is False
    assert boundary["caller_supplied_gain_allowed"] is False
    assert boundary["caller_supplied_scientific_threshold_allowed"] is False
