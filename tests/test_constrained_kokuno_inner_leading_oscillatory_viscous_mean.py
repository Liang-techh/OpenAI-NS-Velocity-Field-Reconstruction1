from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_inner_leading_oscillatory_viscous_mean import (
    AGENT1_HEAD,
    AGENT1_SOURCE_BLOB_SHA,
    AGENT2_HEAD,
    AGENT2_SOURCE_BLOB_SHA,
    ANGULAR_ORDERS,
    REPOSITORY_VISCOSITY,
    _cylindrical_field,
    _materialize_from_laplacian_providers,
    build_mechanics_report,
    materialize_inner_leading_oscillatory_viscous_mean,
    truth_boundary,
)


def _mechanics_witness():
    inner = _cylindrical_field(1.0, -2.0, 3.0, 0.7)
    osc = _cylindrical_field(0.5, 1.0, -1.0, -0.4)
    return _materialize_from_laplacian_providers(
        inner,
        osc,
        np.array([0.7, 1.1]),
        np.array([-0.2, 0.3]),
        np.array([0.48, 0.52]),
        inner_backend=None,
        oscillatory_backend=None,
        backend_kind="manufactured-mechanics-only",
    )


def test_rotating_frame_recovers_registered_viscous_means_and_closure():
    witness = _mechanics_witness()
    expected_inner = np.array([-0.01, 0.02, -0.03])
    expected_osc = np.array([-0.005, -0.01, 0.01])
    expected_total = np.array([-0.015, 0.01, -0.02])

    np.testing.assert_allclose(
        witness.mean_inner_viscous_cylindrical,
        np.broadcast_to(expected_inner, witness.mean_inner_viscous_cylindrical.shape),
        atol=2.0e-14,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        witness.mean_oscillatory_viscous_cylindrical,
        np.broadcast_to(expected_osc, witness.mean_oscillatory_viscous_cylindrical.shape),
        atol=2.0e-14,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        witness.mean_total_viscous_cylindrical,
        np.broadcast_to(expected_total, witness.mean_total_viscous_cylindrical.shape),
        atol=2.0e-14,
        rtol=0.0,
    )
    assert witness.angular_orders == ANGULAR_ORDERS == (32, 64, 128)
    assert max(witness.successive_total_viscous_mean_relative_differences) < 2.0e-14
    assert witness.projected_additive_closure_absolute_max < 2.0e-14
    assert witness.total_viscous_mean_rms > 0.0
    assert witness.full_ring_total_viscous_rms >= witness.total_viscous_mean_rms
    assert witness.total_viscous_mean_to_full_rms_ratio <= 1.0 + 1.0e-12
    assert witness.viscosity == REPOSITORY_VISCOSITY == 0.01


def test_malformed_or_nonfinite_laplacian_payloads_fail_closed():
    good = _cylindrical_field(1.0, 0.0, 0.0, 0.1)

    def wrong_shape(x, y, z, t):
        xb = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )[0]
        return np.zeros(xb.shape + (2,))

    with pytest.raises(ValueError, match="unexpected"):
        _materialize_from_laplacian_providers(
            wrong_shape,
            good,
            [0.8],
            [0.0],
            [0.5],
            inner_backend=None,
            oscillatory_backend=None,
            backend_kind="negative-control",
        )

    def nonfinite(x, y, z, t):
        xb = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )[0]
        out = np.zeros(xb.shape + (3,))
        out[..., 1] = np.nan
        return out

    with pytest.raises(ValueError, match="non-finite"):
        _materialize_from_laplacian_providers(
            good,
            nonfinite,
            [0.8],
            [0.0],
            [0.5],
            inner_backend=None,
            oscillatory_backend=None,
            backend_kind="negative-control",
        )


def test_public_materializer_exposes_no_retuning_or_residual_escape_hatches():
    parameters = inspect.signature(
        materialize_inner_leading_oscillatory_viscous_mean
    ).parameters
    assert list(parameters) == [
        "inner_backend",
        "oscillatory_backend",
        "radius",
        "z",
        "t",
    ]
    forbidden = {
        "residual",
        "defect",
        "mean",
        "stress",
        "inverse",
        "pressure",
        "forcing",
        "target",
        "gain",
        "normalized_score",
        "alpha",
        "damping",
        "angular_order",
        "spatial_step",
        "viscosity",
        "nu",
        "delta_y",
        "delta_a",
        "scientific_threshold",
    }
    assert forbidden.isdisjoint(parameters)


def test_truth_boundary_and_exact_upstream_pins_remain_fail_closed():
    boundary = truth_boundary()
    assert AGENT1_HEAD == "e96ee90144976a992b62d76d4361a94eb16bd91e"
    assert AGENT1_SOURCE_BLOB_SHA == "74b0e185e8e0bbb08695d493e0a091c305a03006"
    assert AGENT2_HEAD == "9f8bae37c4d3b2559bee6db050655702e4c058e9"
    assert AGENT2_SOURCE_BLOB_SHA == "7dfd3dabd28147dffac6ff77ec5457c17689bfe1"
    assert boundary["inner_plus_oscillatory_base_viscous_term_included"] is True
    assert boundary["complete_ns_defect"] is False
    assert boundary["pressure_gradient_included"] is False
    assert boundary["restricted_forcing_included"] is False
    assert boundary["global_corrected_leading_join_materialized"] is False
    assert boundary["real_agent3_delta_a_bound"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["pde_validated"] is False
    assert boundary["final_normalized_momentum_gate"] == 1.0e-3
    assert boundary["final_normalized_divergence_gate"] == 1.0e-5


def test_mechanics_report_is_explicitly_not_real_candidate_evidence():
    report = build_mechanics_report()
    assert report["mechanics_only"] is True
    assert report["backend_kind"] == "manufactured-mechanics-only"
    np.testing.assert_allclose(
        report["expected_total_viscous_mean"], [-0.015, 0.01, -0.02]
    )
    assert report["truth_boundary"]["complete_ns_defect"] is False
