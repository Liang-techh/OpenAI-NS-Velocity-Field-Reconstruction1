"""Named public delivery wrapper for the energy-normalized bipolar Eq. (4.5) candidate.

The wrapped field is an explicit experimental child: it uses the existing odd
``Phi(0,1)`` central-flow seed and the already-recorded positive common scale
that gives ``E(0.25)=1`` under the repository's axisymmetric quadrature.  This
module changes no coefficient relative to the committed normalized artifact and
makes no visualization, PDE, paper-exact, or OpenAI-field promotion.
"""
from __future__ import annotations

from dataclasses import replace
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

from .constrained_eq45_bipolar_seed import bipolar_seed
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .eq45_supported_delivery import Eq45SupportedDeliveryField


SOURCE_BIPOLAR_SHA256 = "8027075948fc4cbbf0c79cf080729b543c7a7d00b011a4b35b582b95deab5970"
NORMALIZED_BIPOLAR_SHA256 = "c0e27269adfb6f305c2c0b5a2483f9eddd54f79691f69cacd7bb735ba592a702"
ENERGY_NORMALIZATION_SCALE = 1.8097870818686452


def normalized_bipolar_candidate() -> Eq45SupportedVelocityCandidate:
    """Reconstruct the checked normalized bipolar candidate and fail on drift."""
    source = bipolar_seed(Eq45SupportedDeliveryField().candidate, phi01=1.0)
    if source.sha256 != SOURCE_BIPOLAR_SHA256:
        raise RuntimeError("source bipolar candidate identity drift")

    basis = source.parent.profile_basis
    candidate = replace(
        source,
        parent=replace(
            source.parent,
            profile_basis=replace(
                basis,
                phi_coefficients=tuple(
                    ENERGY_NORMALIZATION_SCALE * value
                    for value in basis.phi_coefficients
                ),
                swirl_coefficients=tuple(
                    ENERGY_NORMALIZATION_SCALE * value
                    for value in basis.swirl_coefficients
                ),
            ),
        ),
    )
    if candidate.sha256 != NORMALIZED_BIPOLAR_SHA256:
        raise RuntimeError("normalized bipolar candidate identity drift")
    return candidate


class Eq45NormalizedBipolarDeliveryField(Eq45SupportedDeliveryField):
    """User-facing wrapper around the frozen normalized bipolar candidate."""

    def __init__(
        self,
        candidate: Eq45SupportedVelocityCandidate | None = None,
        *,
        candidate_path: str | Path | None = None,
    ) -> None:
        if candidate is None and candidate_path is None:
            candidate = normalized_bipolar_candidate()
        super().__init__(candidate=candidate, candidate_path=candidate_path)

    def metadata(self) -> dict[str, Any]:
        metadata = super().metadata()
        metadata.update(
            {
                "entrypoint": (
                    "openai_ns_reconstruction.eq45_normalized_bipolar_delivery:velocity"
                ),
                "candidate_role": "experimental_source_aligned_energy_normalized",
                "source_bipolar_sha256": SOURCE_BIPOLAR_SHA256,
                "energy_normalization_scale": ENERGY_NORMALIZATION_SCALE,
                "selection_status": "not_canonical_not_visualization_selected",
            }
        )
        return metadata


@lru_cache(maxsize=1)
def default_field() -> Eq45NormalizedBipolarDeliveryField:
    return Eq45NormalizedBipolarDeliveryField()


def velocity(x, y, z, t) -> np.ndarray:
    """Return the frozen normalized bipolar field in Cartesian ``[...,3]`` order."""
    return default_field().velocity(x, y, z, t)


def components(x, y, z, t):
    return default_field().components(x, y, z, t)


def u(x, y, z, t):
    return default_field().u(x, y, z, t)


def v(x, y, z, t):
    return default_field().v(x, y, z, t)


def w(x, y, z, t):
    return default_field().w(x, y, z, t)
