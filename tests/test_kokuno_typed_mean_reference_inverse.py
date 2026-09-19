import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_same_cycle_defect_contract import CycleIdentity
from openai_ns_reconstruction.kokuno_typed_mean_reference_inverse import (
    MeanReferenceInverseOperator,
    deterministic_receipt,
    materialize_typed_mean_reference_inverse,
    truth_boundary,
)


def _operator(**overrides):
    kwargs = {
        "identity": CycleIdentity("unit-reference-inverse", 0, "t0"),
        "radii": (0.2, 0.3),
        "A_c": (2.0, 2.0),
        "u_star": (3.0, 3.0),
        "h_plus": (1.0, 1.0),
        "h_minus": (1.0, 1.0),
        "producer_kind": "unit-test",
        "provenance": "fixed unit-test operator",
        "source_reference_operator_certified": False,
    }
    kwargs.update(overrides)
    return MeanReferenceInverseOperator(**kwargs)


def test_displayed_reference_matrix_and_conditioning() -> None:
    operator = _operator()
    expected = np.array(
        [
            [[-2.0, -2.0], [-3.0, 3.0]],
            [[-2.0, -2.0], [-3.0, 3.0]],
        ]
    )
    np.testing.assert_allclose(operator.matrices(), expected, rtol=0.0, atol=0.0)
    determinants = np.linalg.det(operator.matrices())
    np.testing.assert_allclose(determinants, -12.0, rtol=0.0, atol=2e-15)
    singular_values = np.linalg.svd(operator.matrices(), compute_uv=False)
    np.testing.assert_allclose(
        singular_values[:, 0] / singular_values[:, 1], 1.5, rtol=0.0, atol=2e-15
    )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("A_c", (0.0, 2.0)),
        ("u_star", (3.0, -1.0)),
        ("h_plus", (1.0, 0.0)),
        ("h_minus", (1.0, float("nan"))),
    ],
)
def test_reference_operator_rejects_nonpositive_or_nonfinite_profiles(
    field: str, value: tuple[float, float]
) -> None:
    with pytest.raises(ValueError):
        _operator(**{field: value})


def test_reference_operator_rejects_grid_or_shape_drift() -> None:
    with pytest.raises(ValueError, match="strictly increasing"):
        _operator(radii=(0.3, 0.2))
    with pytest.raises(ValueError, match="exact radial grid"):
        _operator(h_plus=(1.0,))


def test_deterministic_receipt_solves_only_uncertified_regression_operator() -> None:
    receipt = deterministic_receipt()
    assert receipt["schema"] == "kokuno-a3-typed-mean-reference-inverse-v1"
    assert receipt["provenance"]["parent_agent3_pr"] == 696
    assert (
        receipt["provenance"]["parent_agent3_head"]
        == "a8b87637710df0469333e1c3fec03f3713fe4e0e"
    )

    evaluation = receipt["analytic_regression"]["evaluation"]
    assert evaluation["reference_operator"]["source_reference_operator_certified"] is False
    assert (
        evaluation["inverse_rhs"]["reference_scale"]["source_reference_scale_certified"]
        is False
    )
    inverse = evaluation["signed_reference_inverse"]
    assert inverse["channel_ordering"] == ["theta_e2", "axial_e1"]
    assert inverse["delta_y_ordering"] == ["sigma_plus", "sigma_minus"]
    assert inverse["delta_y_vector_rms"] > 0.0
    assert inverse["minimum_abs_determinant"] == pytest.approx(12.0, abs=2e-14)
    assert inverse["maximum_condition_number"] == pytest.approx(1.5, abs=2e-14)
    assert inverse["reconstruction_closure_max_abs"] <= 1e-14


def test_public_api_has_no_precomputed_success_inputs() -> None:
    signature = inspect.signature(materialize_typed_mean_reference_inverse)
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


def test_truth_boundary_remains_fail_closed_for_real_candidate() -> None:
    boundary = truth_boundary()
    assert boundary["raw_same_cycle_providers_required"] is True
    assert boundary["caller_supplied_inverse_rhs_allowed"] is False
    assert boundary["caller_supplied_signed_coordinates_allowed"] is False
    assert boundary["provenance_bearing_reference_operator_required"] is True
    assert boundary["typed_reference_inverse_algebra_executable"] is True
    assert boundary["source_h_ref_materialized_for_full_candidate"] is False
    assert boundary["source_epsilon_materialized_for_full_candidate"] is False
    assert boundary["base_positive_amplitudes_materialized"] is False
    assert boundary["delta_amplitudes_materialized"] is False
    assert boundary["signed_mean_inverse_input_ready"] is False
    assert boundary["public_velocity_correction_materialized"] is False
    assert boundary["real_full_candidate_defect_consumed"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["pde_validated"] is False
