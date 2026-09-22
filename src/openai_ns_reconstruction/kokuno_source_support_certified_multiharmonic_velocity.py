"""Support-certified bounded Kokuno multi-harmonic complete-curl family.

K2-OSC-105 deliberately left the outer-support boundary unassessed because the
axis-safe provider contract exposes only an inner certified zero core.  This
increment adds a *provider execution certificate* for a conservative chart-space
spatial envelope.  It does not infer a source-exact numerical support radius.

For each already-bounded harmonic term the provider certifies that its localized
harmonic has a smooth zero extension outside

    axis_zero_radius_chart < R <= radial_support_max_chart,
    z_support_min_chart <= Z <= z_support_max_chart.

The evaluator merely skips provider/cylindrical evaluation where the provider
has certified exact zero.  It does not multiply a nonzero velocity by a hard
mask, so it does not create a new cutoff layer or alter the complete-curl
algebra.  On active points the same autonomous complex coefficient still
multiplies both vector potential and matched curl velocity before summation.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
import math
from typing import Any, Iterable

import numpy as np

from .kokuno_source_bounded_multiharmonic_velocity import (
    BoundedHarmonicTerm,
    BoundedMultiHarmonicEvaluation,
    MAX_L1_AMPLITUDE,
    MAX_TERMS,
    _bounded_terms,
    _finite_array,
    _normalize_vector_batch,
    _real_cartesian,
    bounded_multiharmonic_contract,
)
from .kokuno_source_physicalized_localized_harmonic import (
    SourcePhysicalScaling,
    physical_to_source_chart_rzt,
    physicalize_source_localized_harmonic,
    source_physical_scaling,
)
from .kokuno_source_zero_data_localized_harmonic import (
    SourceLocalizedZeroDataHarmonic,
)

TASK = "K2-OSC-106"
PARENT_A2_PR = 1149
PARENT_A2_HEAD = "db7a9a4efb4faed2fa46d321fdf5a66b8686de1f"
PARENT_A2_DIAGNOSTIC_SOURCE_BLOB = "1a6d07535471d164ae53e22cac560091fe5b7f85"
BOUNDED_FAMILY_SOURCE_BLOB = "06bc24722b4219419e91cd2d5b08cb1294665f24"
SCHEMA = "kokuno-a2-support-certified-bounded-multiharmonic-v1"


@dataclass(frozen=True)
class SupportCertifiedHarmonicTerm:
    """One bounded harmonic plus a conservative certified spatial zero envelope.

    The numerical bounds are repository/provider execution metadata.  They are
    not claimed to be the unique Kokuno support geometry.

    ``spatial_zero_outside_certified`` means the already-localized harmonic is
    exactly zero outside the declared chart envelope.  The stronger
    ``smooth_zero_extension_certified`` prevents this adapter from laundering a
    discontinuous hard mask into a complete-curl field.
    """

    harmonic: BoundedHarmonicTerm
    radial_support_max_chart: float
    z_support_min_chart: float
    z_support_max_chart: float
    spatial_zero_outside_certified: bool
    smooth_zero_extension_certified: bool


def _canonical_json(payload: Any) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Any) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _validated_support_terms(
    terms: Iterable[SupportCertifiedHarmonicTerm],
) -> tuple[SupportCertifiedHarmonicTerm, ...]:
    out = tuple(terms)
    if not out:
        raise ValueError("support-certified family must contain at least one term")
    if len(out) > MAX_TERMS:
        raise ValueError(f"support-certified family may contain at most {MAX_TERMS} terms")

    base_terms: list[BoundedHarmonicTerm] = []
    for index, term in enumerate(out):
        if not isinstance(term, SupportCertifiedHarmonicTerm):
            raise ValueError(f"term {index} must be SupportCertifiedHarmonicTerm")
        if not isinstance(term.harmonic, BoundedHarmonicTerm):
            raise ValueError(f"term {index}.harmonic must be BoundedHarmonicTerm")
        provider = term.harmonic.provider

        rmax = float(term.radial_support_max_chart)
        zmin = float(term.z_support_min_chart)
        zmax = float(term.z_support_max_chart)
        if not math.isfinite(rmax) or rmax <= float(provider.axis_zero_radius_chart):
            raise ValueError(
                "radial_support_max_chart must be finite and strictly exceed the "
                "provider axis-zero radius"
            )
        if not (math.isfinite(zmin) and math.isfinite(zmax) and zmin < zmax):
            raise ValueError("z support bounds must be finite with z_min < z_max")
        if term.spatial_zero_outside_certified is not True:
            raise ValueError("provider must certify exact zero outside the spatial envelope")
        if term.smooth_zero_extension_certified is not True:
            raise ValueError("provider must certify a smooth zero extension outside support")
        base_terms.append(term.harmonic)

    # Reuse the parent's amplitude/phase/provider-identity bounds exactly.
    _bounded_terms(base_terms)
    return out


class KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity:
    """Bounded complete-curl family with certified outer spatial zero skipping."""

    def __init__(
        self,
        terms: Iterable[SupportCertifiedHarmonicTerm],
        *,
        ell: int,
        h: float,
    ) -> None:
        self._terms = _validated_support_terms(terms)
        self._scaling = source_physical_scaling(ell=ell, h=h)

    @property
    def terms(self) -> tuple[SupportCertifiedHarmonicTerm, ...]:
        return self._terms

    @property
    def scaling(self) -> SourcePhysicalScaling:
        return self._scaling

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent": {
                "pr": PARENT_A2_PR,
                "head": PARENT_A2_HEAD,
                "diagnostic_source_blob_sha1": PARENT_A2_DIAGNOSTIC_SOURCE_BLOB,
                "bounded_family_source_blob_sha1": BOUNDED_FAMILY_SOURCE_BLOB,
            },
            "scaling": {"ell": self._scaling.ell, "h": self._scaling.h},
            "terms": [
                {
                    "provider_id": term.harmonic.provider.provider_id,
                    "provider_semantic_sha256": term.harmonic.provider.provider_semantic_sha256,
                    "axis_zero_radius_chart": float(
                        term.harmonic.provider.axis_zero_radius_chart
                    ),
                    "amplitude": float(term.harmonic.amplitude),
                    "phase_offset": float(term.harmonic.phase_offset),
                    "support_certificate": {
                        "radial_support_max_chart": float(term.radial_support_max_chart),
                        "z_support_min_chart": float(term.z_support_min_chart),
                        "z_support_max_chart": float(term.z_support_max_chart),
                        "spatial_zero_outside_certified": True,
                        "smooth_zero_extension_certified": True,
                        "classification": "repository_provider_execution_metadata",
                        "may_be_conservative_envelope": True,
                    },
                }
                for term in self._terms
            ],
            "truth_boundary": support_certified_multiharmonic_contract(),
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
            ell=self._scaling.ell,
            h=self._scaling.h,
        )

        complex_potential = np.zeros(shape + (3,), dtype=np.complex128)
        complex_velocity = np.zeros(shape + (3,), dtype=np.complex128)
        active_count = np.zeros(shape, dtype=np.int16)

        flat_R = R.reshape(-1)
        flat_theta = theta.reshape(-1)
        flat_Z = Z.reshape(-1)
        flat_T = T.reshape(-1)
        flat_potential = complex_potential.reshape(-1, 3)
        flat_velocity = complex_velocity.reshape(-1, 3)
        flat_count = active_count.reshape(-1)

        for support_term in self._terms:
            term = support_term.harmonic
            provider = term.provider
            active = (
                (flat_R > float(provider.axis_zero_radius_chart))
                & (flat_R <= float(support_term.radial_support_max_chart))
                & (flat_Z >= float(support_term.z_support_min_chart))
                & (flat_Z <= float(support_term.z_support_max_chart))
            )
            if not np.any(active):
                continue

            idx = np.flatnonzero(active)
            harmonic = provider.evaluate_chart(
                flat_R[idx],
                flat_theta[idx],
                flat_Z[idx],
                flat_T[idx],
                self._scaling,
            )
            if not isinstance(harmonic, SourceLocalizedZeroDataHarmonic):
                raise ValueError("provider must return SourceLocalizedZeroDataHarmonic")
            physicalized = physicalize_source_localized_harmonic(
                harmonic,
                ell=self._scaling.ell,
                h=self._scaling.h,
            )
            potential = _normalize_vector_batch(
                physicalized.complex_vector_potential_physical,
                size=idx.size,
                name="complex physical vector potential",
                complex_=True,
            )
            velocity = _normalize_vector_batch(
                physicalized.complex_velocity_cylindrical_physical,
                size=idx.size,
                name="complex physical curl velocity",
                complex_=True,
            )
            coefficient = float(term.amplitude) * np.exp(1j * float(term.phase_offset))
            flat_potential[idx] += coefficient * potential
            flat_velocity[idx] += coefficient * velocity
            flat_count[idx] += 1

        real_cylindrical = 2.0 * np.real(complex_velocity)
        real_cartesian = _real_cartesian(theta, real_cylindrical)
        if not (
            np.all(np.isfinite(complex_potential))
            and np.all(np.isfinite(complex_velocity))
            and np.all(np.isfinite(real_cartesian))
        ):
            raise ValueError("support-certified multi-harmonic evaluation became nonfinite")
        return BoundedMultiHarmonicEvaluation(
            complex_vector_potential_cylindrical_physical=complex_potential,
            complex_velocity_cylindrical_physical=complex_velocity,
            real_pair_velocity_cylindrical_physical=real_cylindrical,
            real_pair_velocity_cartesian_physical=real_cartesian,
            active_term_count=active_count,
        )

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return physical Cartesian velocity with exact certified zero skipping."""
        return self.evaluate(x, y, z, t).real_pair_velocity_cartesian_physical

    def report(self) -> dict[str, Any]:
        return {
            "task": TASK,
            "schema": SCHEMA,
            "semantic_sha256": self.semantic_sha256,
            "term_count": len(self._terms),
            "amplitude_l1": float(sum(term.harmonic.amplitude for term in self._terms)),
            "configuration": self.configuration(),
        }


def support_certified_multiharmonic_contract() -> dict[str, Any]:
    parent = bounded_multiharmonic_contract()
    if parent["finite_complete_curl_harmonic_family_materialized"] is not True:
        raise RuntimeError("bounded complete-curl family prerequisite missing")
    if parent["self_contained_velocity_xyzt_provider"] is not False:
        raise RuntimeError("parent self-contained truth boundary drifted")
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_a2_pr": PARENT_A2_PR,
        "parent_a2_head": PARENT_A2_HEAD,
        "bounded_parent_max_terms": MAX_TERMS,
        "bounded_parent_max_l1_amplitude": MAX_L1_AMPLITUDE,
        "provider_certified_outer_spatial_support_materialized": True,
        "outer_support_exact_zero_skips_provider_evaluation": True,
        "axis_zero_core_still_skips_provider_evaluation": True,
        "support_applied_as_execution_skip_not_velocity_hard_mask": True,
        "smooth_zero_extension_certificate_required": True,
        "support_certificate_identity_bound": True,
        "support_envelope_may_be_conservative": True,
        "source_exact_outer_support_geometry_recovered": False,
        "source_exact_support_numbers_claimed": False,
        "pulse_time_support_certificate_materialized": False,
        "same_coefficient_applied_to_vector_potential_and_curl_velocity": True,
        "complete_curl_linearity_preserved_on_active_support": True,
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
            KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity.velocity
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
                KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity.velocity
            ).parameters
        ),
        "forbidden_velocity_inputs_present": sorted(params & forbidden),
        "paper_exact": False,
        "pde_validated": False,
    }


__all__ = [
    "SupportCertifiedHarmonicTerm",
    "KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity",
    "support_certified_multiharmonic_contract",
    "public_contract",
]
