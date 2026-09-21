"""Identity-bound wrapper for renderer-independent cylindrical morphology receipts.

The measurement engine lives in :mod:`st054_cylindrical_morphology` and is already
implementation-generic: it consumes only ``velocity(x,y,z,t)``.  This module adds
an explicit immutable candidate identity envelope so the same measurements can be
reused for another frozen velocity candidate without relabelling an ST054 receipt
or treating a free-form model name as provenance.

This is delivery/governance plumbing only.  It does not compare against OpenAI
pixels, infer hidden numerical targets, or perform Navier--Stokes acceptance.
"""
from __future__ import annotations

import hashlib
import json
import re
from typing import Any, Sequence

from .st054_cylindrical_morphology import (
    CORE_RADII,
    DEFAULT_AZIMUTHS,
    DEFAULT_RADII,
    DEFAULT_TIMES,
    DEFAULT_Z_MAGNITUDES,
    EXTERNAL_METHOD_SCREEN,
    TRUTH_BOUNDARY,
    VelocityField,
    fingerprint_velocity_field,
)

SCHEMA = "candidate_cylindrical_morphology_fingerprint_v1"
IMPLEMENTATION = "openai_ns_reconstruction.st054_cylindrical_morphology:fingerprint_velocity_field"
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def _require_nonempty_string(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise ValueError(f"{name} must be a nonempty, whitespace-stable string")
    return value


def _require_sha256(value: str, name: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        raise ValueError(f"{name} must be a lowercase 64-hex SHA-256 digest")
    return value


def _measurement_sha256(measurements: dict[str, Any]) -> str:
    return hashlib.sha256(
        json.dumps(measurements, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def validate_candidate_morphology_receipt(receipt: dict[str, Any]) -> None:
    """Fail closed on missing identity, receipt tampering, or truth promotion."""
    if receipt.get("schema") != SCHEMA:
        raise ValueError("unexpected candidate morphology schema")
    if receipt.get("candidate_identity_bound") is not True:
        raise ValueError("candidate identity must be explicitly bound")
    if receipt.get("cylindrical_morphology_diagnostic_ready") is not True:
        raise ValueError("diagnostic readiness must be explicitly true")

    identity = receipt.get("candidate_identity")
    if not isinstance(identity, dict):
        raise ValueError("missing candidate identity")
    _require_nonempty_string(identity.get("candidate_id"), "candidate_id")
    _require_sha256(identity.get("candidate_sha256"), "candidate_sha256")
    _require_sha256(identity.get("velocity_identity_sha256"), "velocity_identity_sha256")

    if receipt.get("diagnostic_implementation") != IMPLEMENTATION:
        raise ValueError("diagnostic implementation drift")

    protocol = receipt.get("protocol")
    if not isinstance(protocol, dict):
        raise ValueError("missing morphology protocol")
    for key in ("source_numeric_targets_used", "renderer_or_camera_used", "pixel_loss_used"):
        if protocol.get(key) is not False:
            raise ValueError(f"{key} must remain false")

    for key, expected in TRUTH_BOUNDARY.items():
        if receipt.get(key) is not expected:
            raise ValueError(f"truth boundary promoted: {key}")

    screen = receipt.get("external_method_screen")
    if not isinstance(screen, list) or screen != list(EXTERNAL_METHOD_SCREEN):
        raise ValueError("external method screen drift")

    measurements = receipt.get("measurements")
    if not isinstance(measurements, dict):
        raise ValueError("missing morphology measurements")
    if receipt.get("measurement_sha256") != _measurement_sha256(measurements):
        raise ValueError("morphology measurement digest mismatch")


def fingerprint_identified_velocity_field(
    field: VelocityField,
    *,
    candidate_id: str,
    candidate_sha256: str,
    velocity_identity_sha256: str,
    radii: Sequence[float] = DEFAULT_RADII,
    z_magnitudes: Sequence[float] = DEFAULT_Z_MAGNITUDES,
    times: Sequence[float] = DEFAULT_TIMES,
    azimuth_count: int = DEFAULT_AZIMUTHS,
    core_radii: Sequence[float] = CORE_RADII,
) -> dict[str, Any]:
    """Fingerprint a velocity field and bind the receipt to explicit frozen identity.

    The three identity values are provenance inputs, not scientific targets.  They
    are validated before evaluating the field and are preserved verbatim in the
    receipt.  The cylindrical measurements and sampling protocol are delegated to
    the already-merged renderer-independent engine without retuning.
    """
    identity = {
        "candidate_id": _require_nonempty_string(candidate_id, "candidate_id"),
        "candidate_sha256": _require_sha256(candidate_sha256, "candidate_sha256"),
        "velocity_identity_sha256": _require_sha256(
            velocity_identity_sha256, "velocity_identity_sha256"
        ),
    }

    base = fingerprint_velocity_field(
        field,
        model_id=candidate_id,
        radii=radii,
        z_magnitudes=z_magnitudes,
        times=times,
        azimuth_count=azimuth_count,
        core_radii=core_radii,
    )

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "candidate_identity": identity,
        "candidate_identity_bound": True,
        "diagnostic_implementation": IMPLEMENTATION,
        "protocol": base["protocol"],
        "measurements": base["measurements"],
        "measurement_sha256": base["measurement_sha256"],
        "external_method_screen": list(EXTERNAL_METHOD_SCREEN),
        "cylindrical_morphology_diagnostic_ready": True,
        **TRUTH_BOUNDARY,
        "scope": (
            "identity_bound_renderer_independent_velocity_kinematics_"
            "not_source_or_pde_acceptance"
        ),
    }
    validate_candidate_morphology_receipt(receipt)
    return receipt
