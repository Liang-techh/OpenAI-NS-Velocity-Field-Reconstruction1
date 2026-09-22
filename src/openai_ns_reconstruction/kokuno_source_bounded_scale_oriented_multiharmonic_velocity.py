"""Bounded isotropic physical scaling for oriented Kokuno complete-curl waves.

The parent K2-OSC-108 field already applies bounded amplitude/phase coefficients,
provider-certified spatial/pulse support, corrected physical Q scaling, exact
complete curls, and one repository-autonomous azimuthal frame orientation.
K2-OSC-109 independently audits that public field at three resolutions.

This adapter adds exactly one remaining low-dimensional candidate coordinate: a
bounded *isotropic physical spatial scale* ``s``.  It acts on the vector
potential before curl evaluation by the constant dilation

    x_parent = x_lab / s,          t_parent = t_lab,
    A_s(x,t) = s A_parent(x/s,t),
    u_s(x,t) = u_parent(x/s,t).

Hence ``curl_x A_s = u_s`` by the chain rule, and
``div_x u_s = s^(-1) div_parent u_parent``.  The complete-curl/divergence-free
contract is therefore preserved without fitting Cartesian velocity components
or multiplying a nonzero velocity by a support mask.  Existing provider support
and pulse-coordinate logic are evaluated in the dilated parent coordinates.

The numerical scale bound and value are repository-autonomous candidate design,
not recovered Kokuno/OpenAI source data.  The source background/forcing/remaining
partition providers are still external, so this is not a self-contained or
paper-exact field and no complete Navier--Stokes residual is assessed here.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from . import kokuno_source_azimuthal_frame_oriented_multiharmonic_velocity as _oriented_module
from .kokuno_source_azimuthal_frame_oriented_multiharmonic_velocity import (
    KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity,
    azimuthal_frame_orientation_contract,
)

TASK = "K2-OSC-110"
SCHEMA = "kokuno-a2-bounded-isotropic-scale-oriented-multiharmonic-v1"
STACK_PARENT_A2_PR = 1180
STACK_PARENT_A2_HEAD = "3bdf0d6a33a629c4262875b146e1378d550023d0"
STACK_PARENT_DIAGNOSTIC_SOURCE_BLOB = "16efe4ab3521f10c4d15f166fa054ab9bef2bae8"
ORIENTED_PARENT_A2_PR = 1170
ORIENTED_PARENT_A2_HEAD = "59d3e6d64fcdfeaf145cbee67e68b8a25e89489c"
ORIENTED_PARENT_SOURCE_BLOB = "ca8acae2c3a53b79528cb55f45502443f6ffb2b6"

SPATIAL_SCALE_MIN = 0.5
SPATIAL_SCALE_MAX = 2.0


def _canonical_json(payload: Any) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Any) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _git_blob_sha1(path: str | Path) -> str:
    data = Path(path).read_bytes()
    prefix = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(prefix + data).hexdigest()


def _finite_array(value: Any, *, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be finite")
    return out


def _validated_scale(value: float) -> float:
    scale = float(value)
    if (
        not math.isfinite(scale)
        or scale < SPATIAL_SCALE_MIN
        or scale > SPATIAL_SCALE_MAX
    ):
        raise ValueError(
            f"spatial_scale must be finite and lie in "
            f"[{SPATIAL_SCALE_MIN}, {SPATIAL_SCALE_MAX}]"
        )
    return scale


def _parent_identity() -> None:
    path = getattr(_oriented_module, "__file__", None)
    if not path or _git_blob_sha1(path) != ORIENTED_PARENT_SOURCE_BLOB:
        raise RuntimeError("exact K2-OSC-108 oriented parent source blob drifted")
    contract = azimuthal_frame_orientation_contract()
    if contract.get("bounded_azimuthal_frame_orientation_materialized") is not True:
        raise RuntimeError("K2-OSC-108 orientation prerequisite is unavailable")
    if contract.get("independent_scale_parameter_materialized") is not False:
        raise RuntimeError("K2-OSC-108 scale truth boundary drifted")
    if contract.get("self_contained_velocity_xyzt_provider") is not False:
        raise RuntimeError("K2-OSC-108 provider truth boundary drifted")


@dataclass(frozen=True)
class BoundedIsotropicScaleEvaluation:
    """Physical Cartesian complete-curl data after one isotropic dilation."""

    complex_vector_potential_cartesian_physical: np.ndarray
    complex_velocity_cartesian_physical: np.ndarray
    real_pair_velocity_cartesian_physical: np.ndarray
    active_term_count: np.ndarray


class KokunoBoundedIsotropicScaleOrientedMultiHarmonicPhysicalVelocity:
    """Bounded spatial-dilation wrapper around exact K2-OSC-108 evaluation."""

    def __init__(
        self,
        parent: KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity,
        *,
        spatial_scale: float,
    ) -> None:
        _parent_identity()
        if not isinstance(parent, KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity):
            raise ValueError(
                "parent must be "
                "KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity"
            )
        self._parent = parent
        self._scale = _validated_scale(spatial_scale)

    @property
    def parent(self) -> KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity:
        return self._parent

    @property
    def spatial_scale(self) -> float:
        return self._scale

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "stack_parent": {
                "pr": STACK_PARENT_A2_PR,
                "head": STACK_PARENT_A2_HEAD,
                "diagnostic_source_blob_sha1": STACK_PARENT_DIAGNOSTIC_SOURCE_BLOB,
            },
            "oriented_parent": {
                "pr": ORIENTED_PARENT_A2_PR,
                "head": ORIENTED_PARENT_A2_HEAD,
                "source_blob_sha1": ORIENTED_PARENT_SOURCE_BLOB,
                "semantic_sha256": self.parent.semantic_sha256,
                "configuration": self.parent.configuration(),
            },
            "spatial_scale": {
                "value": self.spatial_scale,
                "bound": [SPATIAL_SCALE_MIN, SPATIAL_SCALE_MAX],
                "coordinate_map": "x_parent=(x,y,z)_lab/s; physical time unchanged",
                "vector_potential_map": "A_s(x,t)=s*A_parent(x/s,t)",
                "velocity_map": "u_s(x,t)=u_parent(x/s,t)",
                "classification": "repository_autonomous_candidate_design",
                "source_exact_claimed": False,
            },
            "truth_boundary": bounded_isotropic_scale_contract(),
        }

    @property
    def semantic_sha256(self) -> str:
        return _sha256(self.configuration())

    def evaluate(self, x: Any, y: Any, z: Any, t: Any) -> BoundedIsotropicScaleEvaluation:
        x_a = _finite_array(x, name="x")
        y_a = _finite_array(y, name="y")
        z_a = _finite_array(z, name="z")
        t_a = _finite_array(t, name="t")
        x_a, y_a, z_a, t_a = np.broadcast_arrays(x_a, y_a, z_a, t_a)

        inv = 1.0 / self.spatial_scale
        parent_eval = self.parent.evaluate(
            x_a * inv,
            y_a * inv,
            z_a * inv,
            t_a,
        )
        A_parent = np.asarray(
            parent_eval.complex_vector_potential_cartesian_physical,
            dtype=np.complex128,
        )
        u_parent = np.asarray(
            parent_eval.complex_velocity_cartesian_physical,
            dtype=np.complex128,
        )
        real_parent = np.asarray(
            parent_eval.real_pair_velocity_cartesian_physical,
            dtype=float,
        )
        A_scaled = self.spatial_scale * A_parent
        active = np.asarray(parent_eval.active_term_count, dtype=np.int16)

        if not (
            np.all(np.isfinite(A_scaled))
            and np.all(np.isfinite(u_parent))
            and np.all(np.isfinite(real_parent))
        ):
            raise ValueError("bounded-scale complete-curl evaluation became nonfinite")
        return BoundedIsotropicScaleEvaluation(
            complex_vector_potential_cartesian_physical=A_scaled,
            complex_velocity_cartesian_physical=u_parent,
            real_pair_velocity_cartesian_physical=real_parent,
            active_term_count=active,
        )

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return physical Cartesian velocity; scale is frozen at construction time."""
        return self.evaluate(x, y, z, t).real_pair_velocity_cartesian_physical

    def report(self) -> dict[str, Any]:
        return {
            "task": TASK,
            "schema": SCHEMA,
            "semantic_sha256": self.semantic_sha256,
            "spatial_scale": self.spatial_scale,
            "configuration": self.configuration(),
        }


def bounded_isotropic_scale_contract() -> dict[str, Any]:
    _parent_identity()
    parent = azimuthal_frame_orientation_contract()
    return {
        "schema": SCHEMA,
        "task": TASK,
        "stack_parent_a2_pr": STACK_PARENT_A2_PR,
        "stack_parent_a2_head": STACK_PARENT_A2_HEAD,
        "oriented_parent_a2_pr": ORIENTED_PARENT_A2_PR,
        "oriented_parent_a2_head": ORIENTED_PARENT_A2_HEAD,
        "parent_complete_curl_orientation_contract_inherited": True,
        "bounded_isotropic_physical_scale_materialized": True,
        "independent_scale_parameter_materialized": True,
        "spatial_scale_bound": [SPATIAL_SCALE_MIN, SPATIAL_SCALE_MAX],
        "scale_classification": "repository_autonomous_candidate_design",
        "all_three_physical_spatial_coordinates_scaled_together": True,
        "physical_time_coordinate_scaled": False,
        "vector_potential_scaling_rule_materialized": True,
        "vector_potential_amplitude_factor_equals_spatial_scale": True,
        "velocity_component_amplitudes_not_independently_fitted": True,
        "velocity_amplitude_gain_parameter_added": False,
        "curl_contract_preserved_by_constant_dilation_chain_rule": True,
        "divergence_free_contract_preserved_by_constant_dilation": True,
        "provider_support_evaluated_in_dilated_parent_coordinates": True,
        "pulse_coordinate_provider_evaluated_in_dilated_parent_coordinates": True,
        "post_velocity_hard_support_mask_used": False,
        "source_exact_scale_recovered": False,
        "source_exact_scale_claimed": False,
        "source_exact_orientation_recovered": parent["source_exact_orientation_recovered"],
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
            KokunoBoundedIsotropicScaleOrientedMultiHarmonicPhysicalVelocity.velocity
        ).parameters
    )
    forbidden = {
        "spatial_scale",
        "scale",
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
                KokunoBoundedIsotropicScaleOrientedMultiHarmonicPhysicalVelocity.velocity
            ).parameters
        ),
        "forbidden_velocity_inputs_present": sorted(params & forbidden),
        "scale_is_construction_time_bounded_candidate_parameter": True,
        "paper_exact": False,
        "pde_validated": False,
    }


__all__ = [
    "BoundedIsotropicScaleEvaluation",
    "KokunoBoundedIsotropicScaleOrientedMultiHarmonicPhysicalVelocity",
    "SPATIAL_SCALE_MIN",
    "SPATIAL_SCALE_MAX",
    "bounded_isotropic_scale_contract",
    "public_contract",
]
