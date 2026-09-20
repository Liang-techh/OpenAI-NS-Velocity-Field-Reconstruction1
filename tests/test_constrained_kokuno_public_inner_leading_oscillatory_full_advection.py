from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_inner_leading_oscillatory_full_advection import (
    AGENT1_SELF_ADVECTION_BLOB,
    AGENT1_SELF_ADVECTION_HEAD,
    AGENT1_SELF_ADVECTION_PR,
    PARENT_AGENT2_HEAD,
    PARENT_AGENT2_PR,
    directional_fd4_inner_plus_oscillatory_self_advection,
    evaluate_inner_leading_oscillatory_full_advection,
    public_contract,
)


class ManufacturedInnerSelfAdvectionBackend:
    _J = np.asarray(
        [
            [1.0, 0.2, 0.0],
            [-0.3, 1.0, 0.1],
            [0.0, 0.4, 1.0],
        ],
        dtype=float,
    )

    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        del t
        xyz = np.stack((x, y, z), axis=-1)
        return np.einsum("ij,...j->...i", self._J, xyz)

    def self_advection(self, x, y, z, t):
        u = self.velocity(x, y, z, t)
        return np.einsum("ij,...j->...i", self._J, u)


def test_full_convective_expansion_closes_termwise():
    backend = ManufacturedInnerSelfAdvectionBackend()
    x = np.asarray((0.47, 0.58, 0.69, 0.81), dtype=float)
    y = np.asarray((0.19, -0.33, 0.41, -0.27), dtype=float)
    z = np.asarray((-0.82, -0.31, 0.26, 0.74), dtype=float)
    t = np.asarray((0.44, 0.48, 0.52, 0.56), dtype=float)

    result = evaluate_inner_leading_oscillatory_full_advection(
        backend, x, y, z, t
    )

    np.testing.assert_array_equal(
        result.inner_plus_oscillatory_velocity,
        result.inner_leading_velocity + result.oscillatory_velocity,
    )
    np.testing.assert_array_equal(
        result.mixed_cross_advection,
        result.inner_advects_oscillation + result.oscillation_advects_inner,
    )
    np.testing.assert_array_equal(
        result.inner_plus_oscillatory_self_advection,
        result.inner_leading_self_advection
        + result.mixed_cross_advection
        + result.oscillatory_self_advection,
    )
    assert result.inner_plus_oscillatory_self_advection.shape == (4, 3)
    assert np.all(np.isfinite(result.inner_plus_oscillatory_self_advection))


def test_directional_verifier_does_not_call_backend_self_advection():
    class VelocityOnlyForVerifier(ManufacturedInnerSelfAdvectionBackend):
        def self_advection(self, x, y, z, t):
            raise AssertionError("independent directional verifier called self_advection")

    backend = VelocityOnlyForVerifier()
    value = directional_fd4_inner_plus_oscillatory_self_advection(
        backend,
        np.asarray((0.51, 0.72)),
        np.asarray((0.28, -0.34)),
        np.asarray((-0.55, 0.63)),
        np.asarray((0.47, 0.53)),
        directional_step=1.0e-3,
    )
    assert value.shape == (2, 3)
    assert np.all(np.isfinite(value))


def test_backend_and_resolution_contracts_fail_closed():
    class MissingSelf:
        def velocity(self, x, y, z, t):
            shape = np.broadcast(
                np.asarray(x), np.asarray(y), np.asarray(z), np.asarray(t)
            ).shape
            return np.zeros(shape + (3,))

    with pytest.raises(TypeError, match="self_advection"):
        evaluate_inner_leading_oscillatory_full_advection(
            MissingSelf(), 0.4, 0.2, 0.1, 0.5
        )

    class WrongShape(ManufacturedInnerSelfAdvectionBackend):
        def self_advection(self, x, y, z, t):
            return np.zeros((2,), dtype=float)

    with pytest.raises(ValueError, match="shape"):
        evaluate_inner_leading_oscillatory_full_advection(
            WrongShape(), 0.4, 0.2, 0.1, 0.5
        )

    backend = ManufacturedInnerSelfAdvectionBackend()
    with pytest.raises(ValueError, match="inner_leading_spatial_step"):
        evaluate_inner_leading_oscillatory_full_advection(
            backend, 0.4, 0.2, 0.1, 0.5, inner_leading_spatial_step=0.0
        )
    with pytest.raises(ValueError, match="directional_step"):
        directional_fd4_inner_plus_oscillatory_self_advection(
            backend, 0.4, 0.2, 0.1, 0.5, directional_step=np.nan
        )


def test_public_truth_boundary_remains_inner_only_and_nonretuning():
    contract = public_contract()
    assert contract["parent_agent2_pr"] == PARENT_AGENT2_PR == 821
    assert contract["parent_agent2_head"] == PARENT_AGENT2_HEAD
    assert contract["agent1_self_advection_pr"] == AGENT1_SELF_ADVECTION_PR == 829
    assert contract["agent1_self_advection_head"] == AGENT1_SELF_ADVECTION_HEAD
    assert contract["agent1_self_advection_blob"] == AGENT1_SELF_ADVECTION_BLOB
    assert contract["forbidden_public_inputs_present"] == []

    signature = inspect.signature(evaluate_inner_leading_oscillatory_full_advection)
    assert list(signature.parameters) == [
        "inner_leading_backend", "x", "y", "z", "t",
        "inner_leading_spatial_step", "oscillatory_spatial_step",
    ]

    truth = contract["truth_boundary"]
    assert truth["inner_plus_oscillatory_full_advection_executable"] is True
    assert truth["production_reuses_agent1_analytic_inner_self_advection"] is True
    assert truth["production_reuses_a2_complete_curl_oscillatory_field"] is True
    assert truth["direct_verifier_uses_only_public_velocities"] is True
    for key in (
        "agent1_exact_head_ci_assumed_passed",
        "agent4_inner_field_independent_admission_assumed",
        "agent4_a2_oscillatory_independent_admission_assumed",
        "inner_leading_is_final_corrected_fixed_point",
        "global_leading_velocity_materialized",
        "outer_join_materialized",
        "matched_pressure_included",
        "restricted_forcing_included",
        "agent3_mean_radial_chain_reimplemented",
        "correction_velocity_included",
        "complete_kokuno_candidate_assembled",
        "heldout_ns_residual_assessed",
        "residual_reduction_claimed",
        "paper_exact",
        "pde_validated",
    ):
        assert truth[key] is False
