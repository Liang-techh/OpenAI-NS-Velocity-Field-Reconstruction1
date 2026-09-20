from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_correction_transport_increment import (
    VISCOSITY,
    evaluate_correction_transport_increment,
)
from openai_ns_reconstruction.kokuno_public_signed_amplitude_complete_curl import (
    profile_from_delta_a,
)


def _profile(scale: float = 1.0):
    radii = np.linspace(0.32, 1.18, 31)
    s = (radii - radii[0]) / (radii[-1] - radii[0])
    bump = np.sin(np.pi * s) ** 10
    delta = np.stack(
        (
            scale * 0.0077 * bump * (1.0 + 0.05 * np.cos(2.0 * np.pi * s)),
            -scale * 0.0054 * bump * (1.0 - 0.04 * np.sin(2.0 * np.pi * s)),
        ),
        axis=-1,
    )
    delta[[0, -1], :] = 0.0
    return profile_from_delta_a(
        radii,
        delta,
        reference_time=0.50,
        producer_kind="focused-regression",
        provenance="test-only A2 correction transport increment profile",
        source_mean_amplitude_differential_certified=False,
    )


def test_transport_increment_api_has_no_scientific_or_viscosity_escape_hatch() -> None:
    signature = inspect.signature(evaluate_correction_transport_increment)
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
        "alpha",
        "damping",
        "delta_y",
        "delta_a",
        "normalized_score",
        "scientific_threshold",
        "viscosity",
        "nu",
    }
    assert forbidden.isdisjoint(signature.parameters)
    assert VISCOSITY == 0.01


def test_zero_correction_maps_to_exact_zero_transport_increment() -> None:
    correction = _profile(scale=0.0)
    x = np.asarray((0.46, 0.71, 0.93))
    y = np.asarray((0.09, -0.13, 0.16))
    z = np.asarray((-0.37, 0.08, 0.41))
    t = np.asarray((0.45, 0.50, 0.55))
    result = evaluate_correction_transport_increment(correction, x, y, z, t)
    for key in (
        "correction_velocity",
        "correction_velocity_dt",
        "oscillation_advects_correction",
        "correction_advects_oscillation",
        "cross_advection",
        "correction_self_advection",
        "correction_laplacian",
        "correction_viscous_term",
        "correction_transport_increment",
    ):
        assert np.array_equal(np.asarray(result[key]), np.zeros((3, 3)))


def test_transport_increment_is_exact_registered_termwise_assembly() -> None:
    correction = _profile()
    x = np.asarray((0.48, 0.69, 0.88))
    y = np.asarray((0.12, -0.17, 0.19))
    z = np.asarray((-0.43, 0.06, 0.47))
    t = np.asarray((0.45, 0.50, 0.55))
    result = evaluate_correction_transport_increment(correction, x, y, z, t)

    assembled = (
        np.asarray(result["correction_velocity_dt"])
        + np.asarray(result["cross_advection"])
        + np.asarray(result["correction_self_advection"])
        + np.asarray(result["correction_viscous_term"])
    )
    assert np.array_equal(np.asarray(result["correction_transport_increment"]), assembled)
    assert np.array_equal(
        np.asarray(result["correction_viscous_term"]),
        -VISCOSITY * np.asarray(result["correction_laplacian"]),
    )
    assert np.all(np.isfinite(assembled))
    assert np.sqrt(np.mean(np.sum(assembled * assembled, axis=-1))) > 0.0

    handoff = result["handoff_contract"]
    assert handoff["consumer_lane"] == "Kokuno Agent 3 mean/radial bookkeeping"
    assert handoff["complete_ns_defect"] is False
    assert handoff["includes_correction_time_derivative"] is True
    assert handoff["includes_oscillation_correction_cross_advection"] is True
    assert handoff["includes_correction_self_advection"] is True
    assert handoff["includes_correction_viscosity"] is True
    assert handoff["includes_leading_cross_terms"] is False
    assert handoff["includes_pressure_gradient"] is False
    assert handoff["includes_restricted_forcing"] is False


def test_transport_increment_rejects_invalid_resolution_settings() -> None:
    correction = _profile()
    with pytest.raises(ValueError):
        evaluate_correction_transport_increment(
            correction, 0.61, 0.09, 0.13, 0.50, oscillatory_spatial_step=0.0
        )
    with pytest.raises(ValueError):
        evaluate_correction_transport_increment(
            correction, 0.61, 0.09, 0.13, 0.50, correction_spatial_step=0.2
        )
    with pytest.raises(ValueError):
        evaluate_correction_transport_increment(
            correction, 0.61, 0.09, 0.13, 0.50, correction_laplacian_step=np.nan
        )
