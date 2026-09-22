"""Bounded azimuthal frame orientation for source-facing complete-curl waves.

K2-OSC-107 exposes a provider-driven physical Cartesian velocity built from
localized complete curls with spatial and pulse-coordinate support certificates.
The corrected Kokuno reconstruction also carries explicit phase/frame geometry,
including representative frame vectors, but this repository has not recovered a
unique source-exact numerical frame orientation.

This adapter therefore adds exactly one *repository-autonomous* bounded azimuth
``alpha``.  It acts as a rigid physical rotation about the z axis at the level
of both the complex vector potential and its matched complete-curl velocity:

    x_parent = Q(-alpha) x_lab,
    A_alpha(x) = Q(alpha) A_parent(x_parent),
    u_alpha(x) = Q(alpha) u_parent(x_parent).

Because Q is constant and orthogonal, ``curl A`` and ``div u`` are covariant
under this transformation.  No Cartesian velocity component is fitted
independently, no scale parameter is introduced, and no support/cutoff is
multiplied onto an already nonzero velocity.  The exact K2-OSC-107 provider
support logic is evaluated in the rotated parent frame.

The numerical ``alpha`` is candidate-design metadata, not a recovered Kokuno or
OpenAI coefficient.  The background/forcing/remaining partition providers are
still external, so this is not a self-contained or paper-exact field.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
import math
from typing import Any

import numpy as np

from .kokuno_corrected_oscillation_source_ledger import (
    SOURCE_COMMIT,
    SOURCE_REPOSITORY,
    SOURCE_WORKBENCH_BLOB,
    default_kokuno_corrected_oscillation_source_ledger,
)
from .kokuno_source_pulse_coordinate_support_multiharmonic_velocity import (
    KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity,
    pulse_coordinate_support_multiharmonic_contract,
)

TASK = "K2-OSC-108"
PARENT_A2_PR = 1164
PARENT_A2_HEAD = "0cbf74972af6589f7eb9a2edfe21c957d7512a5c"
PARENT_A2_SOURCE_BLOB = "d95265098a9a9eb95ee673fecb3e7b115f2866e6"
SCHEMA = "kokuno-a2-azimuthal-frame-oriented-bounded-multiharmonic-v1"


def _canonical_json(payload: Any) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Any) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _finite_array(value: Any, *, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be finite")
    return out


def _validated_alpha(value: float) -> float:
    alpha = float(value)
    if not math.isfinite(alpha) or alpha < -math.pi or alpha > math.pi:
        raise ValueError("frame_azimuth must be finite and lie in [-pi, pi]")
    return alpha


def _cylindrical_to_cartesian(theta: np.ndarray, vector: np.ndarray) -> np.ndarray:
    """Rotate real or complex cylindrical components into Cartesian components."""
    theta = np.asarray(theta, dtype=float)
    vector = np.asarray(vector)
    if vector.ndim == 0 or vector.shape[-1] != 3:
        raise ValueError("vector must end in three cylindrical components")
    shape = np.broadcast_shapes(theta.shape, vector.shape[:-1])
    theta = np.broadcast_to(theta, shape)
    vector = np.broadcast_to(vector, shape + (3,))
    c = np.cos(theta)
    s = np.sin(theta)
    vr = vector[..., 0]
    vtheta = vector[..., 1]
    vz = vector[..., 2]
    return np.stack((vr * c - vtheta * s, vr * s + vtheta * c, vz), axis=-1)


def _rotate_cartesian_z(vector: np.ndarray, alpha: float) -> np.ndarray:
    vector = np.asarray(vector)
    if vector.ndim == 0 or vector.shape[-1] != 3:
        raise ValueError("vector must end in three Cartesian components")
    c = math.cos(alpha)
    s = math.sin(alpha)
    vx = vector[..., 0]
    vy = vector[..., 1]
    vz = vector[..., 2]
    return np.stack((c * vx - s * vy, s * vx + c * vy, vz), axis=-1)


@dataclass(frozen=True)
class AzimuthalFrameOrientedEvaluation:
    """Physical Cartesian potential/velocity after one rigid frame rotation."""

    complex_vector_potential_cartesian_physical: np.ndarray
    complex_velocity_cartesian_physical: np.ndarray
    real_pair_velocity_cartesian_physical: np.ndarray
    active_term_count: np.ndarray


class KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity:
    """Rigid-z-frame wrapper around exact K2-OSC-107 complete-curl evaluation."""

    def __init__(
        self,
        parent: KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity,
        *,
        frame_azimuth: float,
    ) -> None:
        if not isinstance(parent, KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity):
            raise ValueError(
                "parent must be KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity"
            )
        self._parent = parent
        self._alpha = _validated_alpha(frame_azimuth)

    @property
    def parent(self) -> KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity:
        return self._parent

    @property
    def frame_azimuth(self) -> float:
        return self._alpha

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent": {
                "pr": PARENT_A2_PR,
                "head": PARENT_A2_HEAD,
                "source_blob_sha1": PARENT_A2_SOURCE_BLOB,
                "semantic_sha256": self.parent.semantic_sha256,
                "configuration": self.parent.configuration(),
            },
            "source_provenance": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "workbench_blob_sha1": SOURCE_WORKBENCH_BLOB,
                "source_frame_role": "structural_provenance_only",
            },
            "frame_orientation": {
                "axis": "physical_z",
                "azimuth_radians": self.frame_azimuth,
                "bound": [-math.pi, math.pi],
                "classification": "repository_autonomous_candidate_design",
                "source_exact_claimed": False,
            },
            "truth_boundary": azimuthal_frame_orientation_contract(),
        }

    @property
    def semantic_sha256(self) -> str:
        return _sha256(self.configuration())

    def evaluate(self, x: Any, y: Any, z: Any, t: Any) -> AzimuthalFrameOrientedEvaluation:
        x_a = _finite_array(x, name="x")
        y_a = _finite_array(y, name="y")
        z_a = _finite_array(z, name="z")
        t_a = _finite_array(t, name="t")
        x_a, y_a, z_a, t_a = np.broadcast_arrays(x_a, y_a, z_a, t_a)

        c = math.cos(self.frame_azimuth)
        s = math.sin(self.frame_azimuth)
        # Parent-frame coordinates are Q(-alpha) x_lab.
        xp = c * x_a + s * y_a
        yp = -s * x_a + c * y_a
        theta_parent = np.arctan2(yp, xp)

        parent = self.parent.evaluate(xp, yp, z_a, t_a)
        parent_A_cart = _cylindrical_to_cartesian(
            theta_parent,
            np.asarray(parent.complex_vector_potential_cylindrical_physical, dtype=np.complex128),
        )
        parent_u_cart = _cylindrical_to_cartesian(
            theta_parent,
            np.asarray(parent.complex_velocity_cylindrical_physical, dtype=np.complex128),
        )
        A_cart = _rotate_cartesian_z(parent_A_cart, self.frame_azimuth)
        u_cart = _rotate_cartesian_z(parent_u_cart, self.frame_azimuth)
        real_cart = 2.0 * np.real(u_cart)
        active = np.asarray(parent.active_term_count, dtype=np.int16)

        if not (
            np.all(np.isfinite(A_cart))
            and np.all(np.isfinite(u_cart))
            and np.all(np.isfinite(real_cart))
        ):
            raise ValueError("azimuthal-frame complete-curl evaluation became nonfinite")

        return AzimuthalFrameOrientedEvaluation(
            complex_vector_potential_cartesian_physical=A_cart,
            complex_velocity_cartesian_physical=u_cart,
            real_pair_velocity_cartesian_physical=real_cart,
            active_term_count=active,
        )

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return oriented physical Cartesian velocity; no scientific tuning inputs."""
        return self.evaluate(x, y, z, t).real_pair_velocity_cartesian_physical

    def report(self) -> dict[str, Any]:
        return {
            "task": TASK,
            "schema": SCHEMA,
            "semantic_sha256": self.semantic_sha256,
            "frame_azimuth": self.frame_azimuth,
            "configuration": self.configuration(),
        }


def azimuthal_frame_orientation_contract() -> dict[str, Any]:
    parent = pulse_coordinate_support_multiharmonic_contract()
    if parent["provider_certified_pulse_coordinate_support_materialized"] is not True:
        raise RuntimeError("K2-OSC-107 pulse-coordinate support prerequisite missing")
    if parent["self_contained_velocity_xyzt_provider"] is not False:
        raise RuntimeError("parent self-contained truth boundary drifted")

    ledger = default_kokuno_corrected_oscillation_source_ledger()
    if "representative frame vectors N and K" not in ledger.required_background_inputs:
        raise RuntimeError("corrected-source frame provenance drifted")

    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_a2_pr": PARENT_A2_PR,
        "parent_a2_head": PARENT_A2_HEAD,
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": SOURCE_COMMIT,
        "source_workbench_blob": SOURCE_WORKBENCH_BLOB,
        "parent_complete_curl_and_support_contract_inherited": True,
        "bounded_azimuthal_frame_orientation_materialized": True,
        "orientation_axis": "physical_z",
        "orientation_bound_radians": [-math.pi, math.pi],
        "orientation_classification": "repository_autonomous_candidate_design",
        "same_rigid_rotation_applied_to_vector_potential_and_complete_curl": True,
        "rigid_rotation_preserves_divergence_free_contract_by_covariance": True,
        "axis_and_radial_support_invariant_under_orientation": True,
        "pulse_coordinate_provider_evaluated_in_rotated_parent_frame": True,
        "independent_cartesian_component_fit_used": False,
        "independent_scale_parameter_materialized": False,
        "source_frame_vectors_structural_provenance_consumed": True,
        "source_exact_frame_vectors_recovered": False,
        "source_exact_orientation_recovered": False,
        "source_exact_orientation_claimed": False,
        "project_domain_source_input_provider_materialized": False,
        "self_contained_velocity_xyzt_provider": False,
        "complete_ns_residual_assessed": False,
        "oscillation_before_after_ns_residual_compared": False,
        "residual_reduction_claimed": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def public_contract() -> dict[str, Any]:
    params = set(
        inspect.signature(
            KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity.velocity
        ).parameters
    )
    forbidden = {
        "residual",
        "forcing",
        "pressure",
        "viscosity",
        "gain",
        "threshold",
        "target",
        "tolerance",
        "rtol",
        "atol",
        "panels",
        "steps",
        "optimizer",
        "frame_azimuth",
    }
    return {
        "velocity_inputs": list(
            inspect.signature(
                KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity.velocity
            ).parameters
        ),
        "forbidden_velocity_inputs_present": sorted(params & forbidden),
        "orientation_is_construction_time_bounded_candidate_parameter": True,
        "paper_exact": False,
        "pde_validated": False,
    }


__all__ = [
    "AzimuthalFrameOrientedEvaluation",
    "KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity",
    "azimuthal_frame_orientation_contract",
    "public_contract",
]
