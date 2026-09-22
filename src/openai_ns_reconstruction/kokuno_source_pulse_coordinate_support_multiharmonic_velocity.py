"""Pulse-coordinate support certificates for bounded Kokuno complete curls.

K2-OSC-106 added conservative *spatial* execution certificates to the already
localized, bounded multi-harmonic complete-curl family.  The corrected source
also states that wave coefficients carry a pulse envelope ``P(v)`` with smooth
zero extension, but the repository does not have a source-exact numerical
support interval or a self-contained map from project coordinates to the source
pulse coordinate ``v``.

This increment keeps that distinction explicit.  Each harmonic receives an
identity-bound provider for

    v = v(R, theta, Z, T)

and a conservative interval ``v_min <= v <= v_max`` on which the already
localized harmonic may be nonzero.  The provider must certify exact zero and a
smooth zero extension outside that interval.  Outside it we skip the harmonic
provider entirely.  This is an execution optimization for a field certified to
be zero there; it is *not* multiplication of a nonzero velocity by a hard time
mask, so no new cutoff derivative is silently introduced.

The numerical interval and the map to ``v`` are repository/provider metadata,
not recovered Kokuno constants.  In particular source pulse coordinate ``v``
is not identified with the chart time ``T`` or physical time ``t`` here.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
import math
import re
from typing import Any, Callable, Iterable

import numpy as np

from .kokuno_corrected_oscillation_source_ledger import (
    SOURCE_COMMIT,
    SOURCE_REPOSITORY,
    SOURCE_WORKBENCH_BLOB,
    default_kokuno_corrected_oscillation_source_ledger,
)
from .kokuno_source_bounded_multiharmonic_velocity import (
    BoundedMultiHarmonicEvaluation,
)
from .kokuno_source_physicalized_localized_harmonic import (
    SourcePhysicalScaling,
    physical_to_source_chart_rzt,
)
from .kokuno_source_support_certified_multiharmonic_velocity import (
    KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity,
    SupportCertifiedHarmonicTerm,
    support_certified_multiharmonic_contract,
)

TASK = "K2-OSC-107"
PARENT_A2_PR = 1158
PARENT_A2_HEAD = "f2b7ddb909b8fc2bcf575942ed63ec4a1fbc98bf"
PARENT_A2_SOURCE_BLOB = "5b2bdab4e7fff24ecec8f8cb0fbd8d118a862381"
SCHEMA = "kokuno-a2-pulse-coordinate-support-bounded-multiharmonic-v1"
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")

PulseCoordinateEvaluator = Callable[
    [np.ndarray, np.ndarray, np.ndarray, np.ndarray, SourcePhysicalScaling],
    np.ndarray,
]


@dataclass(frozen=True)
class PulseCoordinateSupportCertifiedHarmonicTerm:
    """One spatially certified harmonic plus a certified pulse-coordinate envelope.

    ``pulse_coordinate_provider_id`` / SHA bind the numerical map used to obtain
    source pulse coordinate ``v`` from chart coordinates.  Neither that map nor
    the numerical interval is claimed to be source-exact.
    """

    spatial: SupportCertifiedHarmonicTerm
    pulse_coordinate_provider_id: str
    pulse_coordinate_provider_semantic_sha256: str
    pulse_v_min: float
    pulse_v_max: float
    pulse_zero_outside_certified: bool
    smooth_pulse_zero_extension_certified: bool
    evaluate_pulse_coordinate: PulseCoordinateEvaluator


def _canonical_json(payload: Any) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Any) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _finite_array(value: Any, *, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _validated_terms(
    terms: Iterable[PulseCoordinateSupportCertifiedHarmonicTerm],
    *,
    ell: int,
    h: float,
) -> tuple[
    tuple[PulseCoordinateSupportCertifiedHarmonicTerm, ...],
    KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity,
]:
    out = tuple(terms)
    if not out:
        raise ValueError("pulse-coordinate support family must contain at least one term")

    spatial_terms: list[SupportCertifiedHarmonicTerm] = []
    for index, term in enumerate(out):
        if not isinstance(term, PulseCoordinateSupportCertifiedHarmonicTerm):
            raise ValueError(
                f"term {index} must be PulseCoordinateSupportCertifiedHarmonicTerm"
            )
        if not isinstance(term.spatial, SupportCertifiedHarmonicTerm):
            raise ValueError(f"term {index}.spatial must be SupportCertifiedHarmonicTerm")
        provider_id = term.pulse_coordinate_provider_id
        provider_sha = term.pulse_coordinate_provider_semantic_sha256
        if not isinstance(provider_id, str) or not provider_id.strip():
            raise ValueError("pulse_coordinate_provider_id must be a nonempty string")
        if not isinstance(provider_sha, str) or not _SHA256_RE.fullmatch(provider_sha):
            raise ValueError(
                "pulse_coordinate_provider_semantic_sha256 must be a lowercase 64-hex digest"
            )
        vmin = float(term.pulse_v_min)
        vmax = float(term.pulse_v_max)
        if not (math.isfinite(vmin) and math.isfinite(vmax) and vmin < vmax):
            raise ValueError("pulse support bounds must be finite with v_min < v_max")
        if term.pulse_zero_outside_certified is not True:
            raise ValueError("provider must certify exact zero outside pulse-coordinate support")
        if term.smooth_pulse_zero_extension_certified is not True:
            raise ValueError("provider must certify smooth zero extension outside pulse support")
        if not callable(term.evaluate_pulse_coordinate):
            raise ValueError("evaluate_pulse_coordinate must be callable")
        spatial_terms.append(term.spatial)

    # Delegate all term-count, bounded-amplitude, phase, duplicate-harmonic-provider,
    # spatial-support and scaling validation to exact K2-OSC-106.
    parent = KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity(
        tuple(spatial_terms), ell=ell, h=h
    )
    return out, parent


def _pulse_values(
    term: PulseCoordinateSupportCertifiedHarmonicTerm,
    R: np.ndarray,
    theta: np.ndarray,
    Z: np.ndarray,
    T: np.ndarray,
    scaling: SourcePhysicalScaling,
) -> np.ndarray:
    raw = term.evaluate_pulse_coordinate(R, theta, Z, T, scaling)
    values = np.asarray(raw, dtype=float)
    expected = (R.size,)
    if values.shape == () and R.size == 1:
        values = values.reshape(1)
    if values.shape != expected:
        raise ValueError(
            "pulse-coordinate provider batch shape does not match requested points: "
            f"expected {expected}, got {values.shape}"
        )
    if not np.all(np.isfinite(values)):
        raise ValueError("pulse-coordinate provider must return finite values")
    return values


def _real_cartesian(theta: np.ndarray, cylindrical: np.ndarray) -> np.ndarray:
    theta = np.asarray(theta, dtype=float)
    cylindrical = np.asarray(cylindrical, dtype=float)
    c = np.cos(theta)
    s = np.sin(theta)
    ur = cylindrical[..., 0]
    uth = cylindrical[..., 1]
    uz = cylindrical[..., 2]
    return np.stack((ur * c - uth * s, ur * s + uth * c, uz), axis=-1)


class KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity:
    """Bounded complete-curl family with spatial and pulse-coordinate zero skipping."""

    def __init__(
        self,
        terms: Iterable[PulseCoordinateSupportCertifiedHarmonicTerm],
        *,
        ell: int,
        h: float,
    ) -> None:
        self._terms, self._spatial_parent = _validated_terms(terms, ell=ell, h=h)
        self._single_term_fields = tuple(
            KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity(
                (term.spatial,), ell=ell, h=h
            )
            for term in self._terms
        )

    @property
    def terms(self) -> tuple[PulseCoordinateSupportCertifiedHarmonicTerm, ...]:
        return self._terms

    @property
    def scaling(self) -> SourcePhysicalScaling:
        return self._spatial_parent.scaling

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent": {
                "pr": PARENT_A2_PR,
                "head": PARENT_A2_HEAD,
                "source_blob_sha1": PARENT_A2_SOURCE_BLOB,
            },
            "source_provenance": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "workbench_blob_sha1": SOURCE_WORKBENCH_BLOB,
                "classification": "structural_support_provenance_only",
            },
            "spatial_parent_configuration": self._spatial_parent.configuration(),
            "pulse_terms": [
                {
                    "harmonic_provider_id": term.spatial.harmonic.provider.provider_id,
                    "pulse_coordinate_provider_id": term.pulse_coordinate_provider_id,
                    "pulse_coordinate_provider_semantic_sha256": (
                        term.pulse_coordinate_provider_semantic_sha256
                    ),
                    "pulse_support_certificate": {
                        "v_min": float(term.pulse_v_min),
                        "v_max": float(term.pulse_v_max),
                        "pulse_zero_outside_certified": True,
                        "smooth_pulse_zero_extension_certified": True,
                        "classification": "repository_provider_execution_metadata",
                        "source_exact_interval_claimed": False,
                        "source_exact_coordinate_map_claimed": False,
                    },
                }
                for term in self._terms
            ],
            "truth_boundary": pulse_coordinate_support_multiharmonic_contract(),
        }

    @property
    def semantic_sha256(self) -> str:
        return _sha256(self.configuration())

    def evaluate(self, x: Any, y: Any, z: Any, t: Any) -> BoundedMultiHarmonicEvaluation:
        x_a = _finite_array(x, name="x")
        y_a = _finite_array(y, name="y")
        z_a = _finite_array(z, name="z")
        t_a = _finite_array(t, name="t")
        x_a, y_a, z_a, t_a = np.broadcast_arrays(x_a, y_a, z_a, t_a)
        shape = x_a.shape

        radius_physical = np.hypot(x_a, y_a)
        theta = np.arctan2(y_a, x_a)
        R, Z, T = physical_to_source_chart_rzt(
            radius_physical,
            z_a,
            t_a,
            ell=self.scaling.ell,
            h=self.scaling.h,
        )

        flat_x = x_a.reshape(-1)
        flat_y = y_a.reshape(-1)
        flat_z = z_a.reshape(-1)
        flat_t = t_a.reshape(-1)
        flat_R = R.reshape(-1)
        flat_theta = theta.reshape(-1)
        flat_Z = Z.reshape(-1)
        flat_T = T.reshape(-1)

        complex_potential = np.zeros(shape + (3,), dtype=np.complex128)
        complex_velocity = np.zeros(shape + (3,), dtype=np.complex128)
        active_count = np.zeros(shape, dtype=np.int16)
        flat_potential = complex_potential.reshape(-1, 3)
        flat_velocity = complex_velocity.reshape(-1, 3)
        flat_count = active_count.reshape(-1)

        for term, one_term_field in zip(self._terms, self._single_term_fields):
            spatial = term.spatial
            harmonic_provider = spatial.harmonic.provider
            spatial_active = (
                (flat_R > float(harmonic_provider.axis_zero_radius_chart))
                & (flat_R <= float(spatial.radial_support_max_chart))
                & (flat_Z >= float(spatial.z_support_min_chart))
                & (flat_Z <= float(spatial.z_support_max_chart))
            )
            spatial_idx = np.flatnonzero(spatial_active)
            if spatial_idx.size == 0:
                continue

            pulse_v = _pulse_values(
                term,
                flat_R[spatial_idx],
                flat_theta[spatial_idx],
                flat_Z[spatial_idx],
                flat_T[spatial_idx],
                self.scaling,
            )
            pulse_active = (
                (pulse_v >= float(term.pulse_v_min))
                & (pulse_v <= float(term.pulse_v_max))
            )
            idx = spatial_idx[pulse_active]
            if idx.size == 0:
                continue

            piece = one_term_field.evaluate(
                flat_x[idx], flat_y[idx], flat_z[idx], flat_t[idx]
            )
            piece_potential = np.asarray(
                piece.complex_vector_potential_cylindrical_physical,
                dtype=np.complex128,
            ).reshape(idx.size, 3)
            piece_velocity = np.asarray(
                piece.complex_velocity_cylindrical_physical,
                dtype=np.complex128,
            ).reshape(idx.size, 3)
            piece_count = np.asarray(piece.active_term_count, dtype=np.int16).reshape(idx.size)
            if not np.array_equal(piece_count, np.ones(idx.size, dtype=np.int16)):
                raise RuntimeError("spatial parent unexpectedly deactivated pulse-interior points")
            flat_potential[idx] += piece_potential
            flat_velocity[idx] += piece_velocity
            flat_count[idx] += piece_count

        real_cylindrical = 2.0 * np.real(complex_velocity)
        real_cartesian = _real_cartesian(theta, real_cylindrical)
        if not (
            np.all(np.isfinite(complex_potential))
            and np.all(np.isfinite(complex_velocity))
            and np.all(np.isfinite(real_cartesian))
        ):
            raise ValueError("pulse-support complete-curl evaluation became nonfinite")

        return BoundedMultiHarmonicEvaluation(
            complex_vector_potential_cylindrical_physical=complex_potential,
            complex_velocity_cylindrical_physical=complex_velocity,
            real_pair_velocity_cylindrical_physical=real_cylindrical,
            real_pair_velocity_cartesian_physical=real_cartesian,
            active_term_count=active_count,
        )

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return physical Cartesian velocity with certified pulse-coordinate skipping."""
        return self.evaluate(x, y, z, t).real_pair_velocity_cartesian_physical

    def report(self) -> dict[str, Any]:
        return {
            "task": TASK,
            "schema": SCHEMA,
            "semantic_sha256": self.semantic_sha256,
            "term_count": len(self._terms),
            "configuration": self.configuration(),
        }


def pulse_coordinate_support_multiharmonic_contract() -> dict[str, Any]:
    parent = support_certified_multiharmonic_contract()
    if parent["provider_certified_outer_spatial_support_materialized"] is not True:
        raise RuntimeError("K2-OSC-106 spatial support prerequisite missing")
    if parent["pulse_time_support_certificate_materialized"] is not False:
        raise RuntimeError("parent pulse-support truth boundary drifted")

    ledger = default_kokuno_corrected_oscillation_source_ledger()
    if "pulse envelope P(v)" not in ledger.source_wave_support:
        raise RuntimeError("corrected-source pulse-envelope provenance drifted")
    if "smooth support extension" not in ledger.cutoff_tail_status:
        raise RuntimeError("corrected-source smooth support provenance drifted")

    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_a2_pr": PARENT_A2_PR,
        "parent_a2_head": PARENT_A2_HEAD,
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": SOURCE_COMMIT,
        "source_workbench_blob": SOURCE_WORKBENCH_BLOB,
        "provider_certified_outer_spatial_support_inherited": True,
        "provider_certified_pulse_coordinate_support_materialized": True,
        "pulse_coordinate_mapping_identity_bound": True,
        "pulse_coordinate_exact_zero_skips_harmonic_provider_evaluation": True,
        "pulse_coordinate_support_applied_as_execution_skip_not_velocity_hard_mask": True,
        "smooth_pulse_zero_extension_certificate_required": True,
        "same_complete_curl_harmonic_and_bounded_coefficients_preserved": True,
        "source_public_pulse_envelope_structure_consumed": True,
        "source_exact_pulse_coordinate_mapping_recovered": False,
        "source_exact_pulse_support_interval_recovered": False,
        "source_exact_pulse_support_numbers_claimed": False,
        "source_pulse_coordinate_identified_with_chart_T": False,
        "source_pulse_coordinate_identified_with_physical_t": False,
        "physical_time_support_certificate_materialized": False,
        "pulse_time_support_certificate_materialized": False,
        "caller_supplies_corrected_background_and_forcing": True,
        "caller_supplies_remaining_partition_and_pulse_cutoff": True,
        "source_exact_duhamel_frame_B_materialized": False,
        "source_exact_duhamel_propagator_Vm_materialized": False,
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
            KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity.velocity
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
    }
    return {
        "velocity_inputs": list(
            inspect.signature(
                KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity.velocity
            ).parameters
        ),
        "forbidden_velocity_inputs_present": sorted(params & forbidden),
        "paper_exact": False,
        "pde_validated": False,
    }


__all__ = [
    "PulseCoordinateSupportCertifiedHarmonicTerm",
    "KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity",
    "pulse_coordinate_support_multiharmonic_contract",
    "public_contract",
]
