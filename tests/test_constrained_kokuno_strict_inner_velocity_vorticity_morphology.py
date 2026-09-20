from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_strict_inner_velocity_vorticity_morphology import (
    ANGULAR_RESOLUTIONS,
    SOURCE_RINGS,
    frozen_protocol,
    materialize_strict_inner_morphology,
    public_contract,
)


class _AxisymmetricCandidate:
    velocity_candidate_sha256 = "a" * 64
    vorticity_sha256 = "b" * 64

    @staticmethod
    def velocity(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        return np.stack(
            (-0.2 * x - y, -0.2 * y + x, 0.3 * z + 0.05 * t),
            axis=-1,
        )

    @staticmethod
    def vorticity(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        zero = np.zeros_like(x)
        return np.stack((zero, zero, 2.0 + 0.0 * t), axis=-1)


def _mapper(X, eta, t, theta):
    X, eta, t, theta = np.broadcast_arrays(
        np.asarray(X, dtype=float),
        np.asarray(eta, dtype=float),
        np.asarray(t, dtype=float),
        np.asarray(theta, dtype=float),
    )
    return {
        "x": X * np.cos(theta),
        "y": X * np.sin(theta),
        "z": eta,
        "t": t,
    }


def test_frozen_protocol_is_three_resolution_and_strict_inner_only():
    protocol = frozen_protocol()
    assert protocol["angular_resolutions"] == list(ANGULAR_RESOLUTIONS)
    assert len(protocol["source_rings"]) == len(SOURCE_RINGS) == 6
    assert len(protocol["protocol_sha256"]) == 64
    truth = protocol["truth_boundary"]
    assert truth["strict_inner_three_resolution_morphology_protocol_materialized"] is True
    assert truth["renderer_or_camera_used"] is False
    assert truth["source_numeric_targets_used"] is False
    assert truth["whole_domain_velocity_morphology_verified"] is False
    assert truth["whole_domain_vorticity_morphology_verified"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False


def test_axisymmetric_field_has_resolution_and_offgrid_stable_morphology():
    receipt = materialize_strict_inner_morphology(_AxisymmetricCandidate(), _mapper)
    assert receipt["strict_inner_three_resolution_morphology_gate_passed"] is True
    assert set(receipt["levels"]) == {str(n) for n in ANGULAR_RESOLUTIONS}
    assert receipt["candidate_identity"]["velocity_candidate_sha256"] == "a" * 64
    assert receipt["candidate_identity"]["vorticity_sha256"] == "b" * 64

    for key, value in receipt["gate_values"].items():
        assert np.isfinite(value), key
    assert receipt["gate_values"]["velocity_medium_fine_relative_rms"] < 1.0e-12
    assert receipt["gate_values"]["vorticity_medium_fine_relative_rms"] < 1.0e-12
    assert receipt["gate_values"]["velocity_offgrid_fine_relative_rms"] < 1.0e-12
    assert receipt["gate_values"]["vorticity_offgrid_fine_relative_rms"] < 1.0e-12
    assert receipt["gate_values"]["fine_velocity_rms"] > 0.0
    assert receipt["gate_values"]["fine_vorticity_rms"] > 0.0
    assert len(receipt["receipt_sha256"]) == 64


def test_report_records_cylindrical_velocity_and_vorticity_low_modes():
    receipt = materialize_strict_inner_morphology(_AxisymmetricCandidate(), _mapper)
    fine = receipt["levels"][str(ANGULAR_RESOLUTIONS[-1])]
    for field in ("velocity", "vorticity"):
        report = fine[field]
        assert len(report["component_mean"]) == len(SOURCE_RINGS)
        assert len(report["component_rms"]) == len(SOURCE_RINGS)
        assert len(report["vector_norm_rms"]) == len(SOURCE_RINGS)
        assert set(report["harmonic_magnitudes"]) == {"m1_abs", "m2_abs"}


def test_materializer_exposes_no_resolution_or_threshold_knobs():
    params = set(inspect.signature(materialize_strict_inner_morphology).parameters)
    assert params == {"candidate", "coordinate_mapper"}
    contract = public_contract()
    assert contract["scientific_knobs_exposed"] is False
    assert contract["protocol"]["angular_resolutions"] == [16, 32, 64]


def test_nonfinite_candidate_values_fail_closed():
    class _Bad(_AxisymmetricCandidate):
        @staticmethod
        def vorticity(x, y, z, t):
            out = _AxisymmetricCandidate.vorticity(x, y, z, t)
            out[..., 0] = np.nan
            return out

    with pytest.raises(ValueError, match="must be finite"):
        materialize_strict_inner_morphology(_Bad(), _mapper)


def test_bad_coordinate_mapper_fails_closed():
    def bad_mapper(X, eta, t, theta):
        raw = _mapper(X, eta, t, theta)
        raw["x"] = np.asarray(raw["x"])[..., :-1]
        return raw

    with pytest.raises(ValueError, match="coordinate_mapper x must have shape"):
        materialize_strict_inner_morphology(_AxisymmetricCandidate(), bad_mapper)


def test_candidate_identity_must_be_checksum_bound():
    class _BadIdentity(_AxisymmetricCandidate):
        velocity_candidate_sha256 = "not-a-sha"

    with pytest.raises(ValueError, match="velocity_candidate_sha256"):
        materialize_strict_inner_morphology(_BadIdentity(), _mapper)
