from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction.candidate_cylindrical_morphology import (
    SCHEMA,
    fingerprint_identified_velocity_field,
    validate_candidate_morphology_receipt,
)
from openai_ns_reconstruction.st054_cylindrical_morphology import (
    TRUTH_BOUNDARY,
    fingerprint_velocity_field,
)


class ManufacturedVelocityField:
    def velocity(self, x, y, z, t):
        x, y, z = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
        )
        inward = 0.2 + 0.05 * float(t)
        omega = 0.7 + 0.2 * np.hypot(x, y)
        return np.stack(
            (
                -inward * x - omega * y,
                -inward * y + omega * x,
                0.3 * z,
            ),
            axis=-1,
        )


KWARGS = dict(
    radii=(0.5, 1.0),
    z_magnitudes=(0.4,),
    times=(0.25, 0.75),
    azimuth_count=8,
    core_radii=(0.25, 0.5, 0.75, 1.0),
)


def test_identity_bound_wrapper_reuses_measurements_without_st054_identity_laundering():
    field = ManufacturedVelocityField()
    direct = fingerprint_velocity_field(field, model_id="synthetic", **KWARGS)
    receipt = fingerprint_identified_velocity_field(
        field,
        candidate_id="synthetic-frozen-v1",
        candidate_sha256="a" * 64,
        velocity_identity_sha256="b" * 64,
        **KWARGS,
    )

    assert receipt["schema"] == SCHEMA
    assert receipt["candidate_identity"] == {
        "candidate_id": "synthetic-frozen-v1",
        "candidate_sha256": "a" * 64,
        "velocity_identity_sha256": "b" * 64,
    }
    assert receipt["candidate_identity_bound"] is True
    assert receipt["measurements"] == direct["measurements"]
    assert receipt["measurement_sha256"] == direct["measurement_sha256"]
    assert len(receipt["candidate_measurement_binding_sha256"]) == 64
    assert "source_commit" not in receipt
    assert "source_mat_sha256" not in receipt
    for key, expected in TRUTH_BOUNDARY.items():
        assert receipt[key] is expected


def test_identity_is_required_before_field_evaluation():
    class ExplodingField:
        def velocity(self, x, y, z, t):
            raise AssertionError("field must not be evaluated when identity is malformed")

    with pytest.raises(ValueError, match="candidate_id"):
        fingerprint_identified_velocity_field(
            ExplodingField(),
            candidate_id=" ",
            candidate_sha256="a" * 64,
            velocity_identity_sha256="b" * 64,
            **KWARGS,
        )
    with pytest.raises(ValueError, match="candidate_sha256"):
        fingerprint_identified_velocity_field(
            ExplodingField(),
            candidate_id="candidate",
            candidate_sha256="A" * 64,
            velocity_identity_sha256="b" * 64,
            **KWARGS,
        )
    with pytest.raises(ValueError, match="velocity_identity_sha256"):
        fingerprint_identified_velocity_field(
            ExplodingField(),
            candidate_id="candidate",
            candidate_sha256="a" * 64,
            velocity_identity_sha256="not-a-digest",
            **KWARGS,
        )


def test_candidate_receipt_verifier_rejects_identity_tampering_and_truth_promotion():
    receipt = fingerprint_identified_velocity_field(
        ManufacturedVelocityField(),
        candidate_id="synthetic-frozen-v1",
        candidate_sha256="a" * 64,
        velocity_identity_sha256="b" * 64,
        **KWARGS,
    )
    validate_candidate_morphology_receipt(receipt)

    unbound = copy.deepcopy(receipt)
    unbound["candidate_identity_bound"] = False
    with pytest.raises(ValueError, match="identity must be explicitly bound"):
        validate_candidate_morphology_receipt(unbound)

    malformed_identity = copy.deepcopy(receipt)
    malformed_identity["candidate_identity"]["candidate_sha256"] = "c" * 63
    with pytest.raises(ValueError, match="candidate_sha256"):
        validate_candidate_morphology_receipt(malformed_identity)

    valid_shape_identity_tamper = copy.deepcopy(receipt)
    valid_shape_identity_tamper["candidate_identity"]["candidate_sha256"] = "c" * 64
    with pytest.raises(ValueError, match="binding digest mismatch"):
        validate_candidate_morphology_receipt(valid_shape_identity_tamper)

    valid_shape_velocity_identity_tamper = copy.deepcopy(receipt)
    valid_shape_velocity_identity_tamper["candidate_identity"]["velocity_identity_sha256"] = "d" * 64
    with pytest.raises(ValueError, match="binding digest mismatch"):
        validate_candidate_morphology_receipt(valid_shape_velocity_identity_tamper)

    promoted = copy.deepcopy(receipt)
    promoted["visualization_ready"] = True
    with pytest.raises(ValueError, match="truth boundary promoted"):
        validate_candidate_morphology_receipt(promoted)

    invented_target = copy.deepcopy(receipt)
    invented_target["protocol"]["source_numeric_targets_used"] = True
    with pytest.raises(ValueError, match="source_numeric_targets_used"):
        validate_candidate_morphology_receipt(invented_target)

    protocol_tampered = copy.deepcopy(receipt)
    protocol_tampered["protocol"]["radii"][0] += 0.01
    with pytest.raises(ValueError, match="binding digest mismatch"):
        validate_candidate_morphology_receipt(protocol_tampered)

    tampered = copy.deepcopy(receipt)
    tampered["measurements"]["ring_rows"][0]["mean_u_r"] += 1.0
    with pytest.raises(ValueError, match="digest mismatch"):
        validate_candidate_morphology_receipt(tampered)
