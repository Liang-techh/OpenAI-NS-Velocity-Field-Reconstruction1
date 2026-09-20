import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_same_cycle_defect_contract import CycleIdentity
from openai_ns_reconstruction.kokuno_typed_mean_amplitude_differential import (
    MeanReferenceBaseAmplitudes,
    deterministic_receipt,
    materialize_typed_mean_amplitude_differential,
    truth_boundary,
)


def _base_amplitudes(**overrides):
    kwargs = {
        "identity": CycleIdentity("unit-amplitude-differential", 0, "t0"),
        "radii": (0.2, 0.3),
        "a_plus": (2.0, 2.5),
        "a_minus": (3.0, 3.5),
        "producer_kind": "unit-test",
        "provenance": "fixed unit-test positive amplitudes",
        "source_reference_base_amplitudes_certified": False,
    }
    kwargs.update(overrides)
    return MeanReferenceBaseAmplitudes(**kwargs)


def test_positive_base_amplitude_witness_preserves_signed_ordering() -> None:
    amplitudes = _base_amplitudes()
    np.testing.assert_allclose(
        amplitudes.stacked(),
        np.array([[2.0, 3.0], [2.5, 3.5]]),
        rtol=0.0,
        atol=0.0,
    )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("a_plus", (0.0, 2.5)),
        ("a_plus", (-1.0, 2.5)),
        ("a_minus", (3.0, 0.0)),
        ("a_minus", (3.0, float("nan"))),
    ],
)
def test_base_amplitude_witness_rejects_nonpositive_or_nonfinite_profiles(
    field: str, value: tuple[float, float]
) -> None:
    with pytest.raises(ValueError):
        _base_amplitudes(**{field: value})


def test_base_amplitude_witness_rejects_grid_or_shape_drift() -> None:
    with pytest.raises(ValueError, match="strictly increasing"):
        _base_amplitudes(radii=(0.3, 0.2))
    with pytest.raises(ValueError, match="exact radial grid"):
        _base_amplitudes(a_plus=(2.0,))


def test_deterministic_receipt_materializes_only_uncertified_amplitude_differential() -> None:
    receipt = deterministic_receipt()
    assert receipt["schema"] == "kokuno-a3-typed-mean-amplitude-differential-v1"
    assert receipt["provenance"]["parent_agent3_pr"] == 703
    assert (
        receipt["provenance"]["parent_agent3_head"]
        == "56957cd51c7e6638608f28f8fa21df590194dc84"
    )
    assert (
        receipt["provenance"]["source_corrected_reader_commit"]
        == "143f6773feb424ad9ed3a8d116653200f20346b7"
    )

    evaluation = receipt["analytic_regression"]["evaluation"]
    assert (
        evaluation["base_amplitudes"]["source_reference_base_amplitudes_certified"]
        is False
    )
    assert (
        evaluation["reference_inverse"]["reference_operator"]
        ["source_reference_operator_certified"]
        is False
    )
    assert (
        evaluation["reference_inverse"]["inverse_rhs"]["reference_scale"]
        ["source_reference_scale_certified"]
        is False
    )

    amplitude = evaluation["amplitude_differential"]
    base = np.stack(
        (
            np.asarray(evaluation["base_amplitudes"]["a_plus"], dtype=float),
            np.asarray(evaluation["base_amplitudes"]["a_minus"], dtype=float),
        ),
        axis=-1,
    )
    delta_y = np.asarray(amplitude["delta_y"], dtype=float)
    delta_a = np.asarray(amplitude["delta_a"], dtype=float)
    np.testing.assert_allclose(2.0 * base * delta_a, delta_y, rtol=0.0, atol=1e-15)
    assert amplitude["delta_a_vector_rms"] > 0.0
    assert amplitude["multiplication_back_closure_max_abs"] <= 1e-14
    assert amplitude["maximum_relative_differential"] > 0.0
    assert amplitude["full_step_positivity_is_diagnostic_not_admission_gate"] is True


def test_public_api_has_no_precomputed_success_inputs() -> None:
    signature = inspect.signature(materialize_typed_mean_amplitude_differential)
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
        "correction",
    }
    assert forbidden.isdisjoint(signature.parameters)


def test_truth_boundary_remains_fail_closed_for_real_candidate_and_velocity() -> None:
    boundary = truth_boundary()
    assert boundary["raw_same_cycle_providers_required"] is True
    assert boundary["caller_supplied_delta_y_allowed"] is False
    assert boundary["caller_supplied_delta_a_allowed"] is False
    assert boundary["provenance_bearing_positive_base_amplitudes_required"] is True
    assert boundary["typed_amplitude_differential_executable"] is True
    assert boundary["source_epsilon_materialized_for_full_candidate"] is False
    assert boundary["source_h_ref_materialized_for_full_candidate"] is False
    assert boundary["source_base_amplitudes_materialized_for_full_candidate"] is False
    assert boundary["real_candidate_delta_amplitudes_materialized"] is False
    assert boundary["public_velocity_correction_materialized"] is False
    assert boundary["agent2_complete_curl_reimplemented"] is False
    assert boundary["signed_mean_inverse_input_ready_for_real_candidate"] is False
    assert boundary["real_full_candidate_defect_consumed"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["pde_validated"] is False
