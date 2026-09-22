"""Bounded finite Kokuno oscillatory family from axis-safe complete-curl harmonics.

This Kokuno Agent-2 increment composes a *small fixed number* of already
localized, provider-driven source harmonics from K2-OSC-103.  Each provider
still owns the missing corrected background / forcing / support realization.
The new step is only the finite harmonic-family algebra.

For one positive-harmonic representative the parent source route returns a
complex vector potential ``A_m`` and its matched complete-curl velocity
``u_m = curl A_m``.  The real source wave is represented by the conjugate pair
``2 Re(u_m)``.  Here a bounded autonomous coefficient

    c_j = a_j exp(i phi_j),   0 <= a_j <= 1,

is applied to *both* the complex vector potential and its matched complex curl
velocity.  Summing at most four such terms gives

    A_family = sum_j c_j A_j,
    u_family = sum_j c_j u_j,
    u_real   = 2 Re(u_family).

Because curl is linear, this preserves the complete-curl/divergence-free
contract of every admitted term.  The amplitudes/phases below are repository
candidate parameters, not values recovered from Kokuno's corrected source.
They are deliberately bounded: at most four terms and L1 amplitude <= 1.
No residual, forcing, pressure, optimizer, gain, tolerance or target enters the
surface.

The evaluator remains provider-driven rather than self-contained.  It also
retains the K2-OSC-103 certified-axis rule *per provider*: a term is skipped on
its certified zero core and no cylindrical 1/R evaluation is attempted there.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
import math
from typing import Any, Iterable

import numpy as np

from .kokuno_source_axis_safe_physical_velocity import (
    AxisSafeSourceHarmonicProvider,
    PARENT_A2_HEAD as AXIS_SAFE_PARENT_HEAD,
    TASK as AXIS_SAFE_PARENT_TASK,
    _validated_provider,
    source_axis_safe_physical_velocity_contract,
)
from .kokuno_source_physicalized_localized_harmonic import (
    SourcePhysicalScaling,
    physical_to_source_chart_rzt,
    physicalize_source_localized_harmonic,
    source_physical_scaling,
)
from .kokuno_source_zero_data_localized_harmonic import (
    SourceLocalizedZeroDataHarmonic,
    source_zero_data_localized_harmonic_contract,
)

TASK = "K2-OSC-104"
PARENT_A2_PR = 1132
PARENT_A2_HEAD = "1693c6b0e0e337e7fab3914bbcfcab8ad54bbf3a"
PARENT_A2_SOURCE_BLOB = "9067cb68dc49e08a3da7aaa429366f11d2073cdd"
SCHEMA = "kokuno-a2-bounded-multiharmonic-physical-velocity-v1"
MAX_TERMS = 4
MAX_ABS_AMPLITUDE = 1.0
MAX_L1_AMPLITUDE = 1.0
NONTRIVIALITY_FLOOR = 1.0e-12


@dataclass(frozen=True)
class BoundedHarmonicTerm:
    """One bounded autonomous coefficient multiplying one source provider.

    ``amplitude`` and ``phase_offset`` are autonomous repository candidate
    parameters.  They are not corrected-source recovered data and they never
    depend on a residual in this module.
    """

    provider: AxisSafeSourceHarmonicProvider
    amplitude: float
    phase_offset: float


@dataclass(frozen=True)
class BoundedMultiHarmonicEvaluation:
    """Identity-preserving family evaluation in physical coordinates."""

    complex_vector_potential_cylindrical_physical: np.ndarray
    complex_velocity_cylindrical_physical: np.ndarray
    real_pair_velocity_cylindrical_physical: np.ndarray
    real_pair_velocity_cartesian_physical: np.ndarray
    active_term_count: np.ndarray


def _finite_array(value: Any, *, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be finite")
    return out


def _canonical_json(payload: Any) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Any) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _bounded_terms(terms: Iterable[BoundedHarmonicTerm]) -> tuple[BoundedHarmonicTerm, ...]:
    out = tuple(terms)
    if not (1 <= len(out) <= MAX_TERMS):
        raise ValueError(f"harmonic family must contain between 1 and {MAX_TERMS} terms")

    seen: set[tuple[str, str]] = set()
    l1 = 0.0
    nontrivial = False
    for index, term in enumerate(out):
        if not isinstance(term, BoundedHarmonicTerm):
            raise ValueError(f"term {index} must be BoundedHarmonicTerm")
        provider = _validated_provider(term.provider)
        key = (provider.provider_id, provider.provider_semantic_sha256)
        if key in seen:
            raise ValueError("duplicate provider identity would silently stack one harmonic twice")
        seen.add(key)

        amplitude = float(term.amplitude)
        phase = float(term.phase_offset)
        if not math.isfinite(amplitude) or not (0.0 <= amplitude <= MAX_ABS_AMPLITUDE):
            raise ValueError("term amplitude must be finite and lie in [0,1]")
        if not math.isfinite(phase) or not (-math.pi <= phase <= math.pi):
            raise ValueError("term phase_offset must be finite and lie in [-pi,pi]")
        l1 += amplitude
        nontrivial = nontrivial or amplitude >= NONTRIVIALITY_FLOOR

    if l1 > MAX_L1_AMPLITUDE + 32.0 * np.finfo(float).eps:
        raise ValueError("family L1 amplitude exceeds the frozen bound 1")
    if not nontrivial:
        raise ValueError("family amplitudes collapse below the nontriviality floor")
    return out


def _normalize_vector_batch(array: Any, *, size: int, name: str, complex_: bool) -> np.ndarray:
    dtype = np.complex128 if complex_ else float
    out = np.asarray(array, dtype=dtype)
    if out.shape == (3,) and size == 1:
        out = out.reshape(1, 3)
    expected = (size, 3)
    if out.shape != expected:
        raise ValueError(f"{name} must have shape {expected}, got {out.shape}")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be finite")
    return out


def _real_cartesian(theta: np.ndarray, cylindrical: np.ndarray) -> np.ndarray:
    theta = np.asarray(theta, dtype=float)
    cylindrical = np.asarray(cylindrical, dtype=float)
    if cylindrical.shape != theta.shape + (3,):
        raise ValueError("cylindrical batch is incompatible with theta")
    c = np.cos(theta)
    s = np.sin(theta)
    ur = cylindrical[..., 0]
    uth = cylindrical[..., 1]
    uz = cylindrical[..., 2]
    return np.stack((ur * c - uth * s, ur * s + uth * c, uz), axis=-1)


class KokunoBoundedMultiHarmonicPhysicalVelocity:
    """Low-dimensional bounded physical velocity from complete-curl harmonics."""

    def __init__(
        self,
        terms: Iterable[BoundedHarmonicTerm],
        *,
        ell: int,
        h: float,
    ) -> None:
        self._terms = _bounded_terms(terms)
        self._scaling = source_physical_scaling(ell=ell, h=h)

    @property
    def terms(self) -> tuple[BoundedHarmonicTerm, ...]:
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
                "source_blob_sha1": PARENT_A2_SOURCE_BLOB,
                "parent_task": AXIS_SAFE_PARENT_TASK,
                "parent_declared_head": AXIS_SAFE_PARENT_HEAD,
            },
            "scaling": {"ell": self._scaling.ell, "h": self._scaling.h},
            "bounds": {
                "max_terms": MAX_TERMS,
                "max_abs_amplitude": MAX_ABS_AMPLITUDE,
                "max_l1_amplitude": MAX_L1_AMPLITUDE,
                "phase_interval": [-math.pi, math.pi],
                "nontriviality_floor": NONTRIVIALITY_FLOOR,
                "classification": "repository_autonomous_candidate_bounds",
            },
            "terms": [
                {
                    "provider_id": term.provider.provider_id,
                    "provider_semantic_sha256": term.provider.provider_semantic_sha256,
                    "axis_zero_radius_chart": float(term.provider.axis_zero_radius_chart),
                    "amplitude": float(term.amplitude),
                    "phase_offset": float(term.phase_offset),
                }
                for term in self._terms
            ],
            "truth_boundary": bounded_multiharmonic_contract(),
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

        for term in self._terms:
            provider = term.provider
            active = flat_R > float(provider.axis_zero_radius_chart)
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
            raise ValueError("bounded multi-harmonic evaluation became nonfinite")
        return BoundedMultiHarmonicEvaluation(
            complex_vector_potential_cylindrical_physical=complex_potential,
            complex_velocity_cylindrical_physical=complex_velocity,
            real_pair_velocity_cylindrical_physical=real_cylindrical,
            real_pair_velocity_cartesian_physical=real_cartesian,
            active_term_count=active_count,
        )

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return the bounded real multi-harmonic Cartesian velocity."""
        return self.evaluate(x, y, z, t).real_pair_velocity_cartesian_physical

    def report(self) -> dict[str, Any]:
        cfg = self.configuration()
        return {
            "task": TASK,
            "schema": SCHEMA,
            "semantic_sha256": self.semantic_sha256,
            "term_count": len(self._terms),
            "amplitude_l1": float(sum(term.amplitude for term in self._terms)),
            "configuration": cfg,
        }


def bounded_multiharmonic_contract() -> dict[str, Any]:
    """Return the source/autonomous/pending boundary for the family."""
    parent = source_axis_safe_physical_velocity_contract()
    localized = source_zero_data_localized_harmonic_contract()
    if parent["provider_driven_velocity_xyzt_materialized"] is not True:
        raise RuntimeError("K2-OSC-103 provider-driven velocity prerequisite missing")
    if localized["real_conjugate_pair_velocity_materialized_in_source_chart"] is not True:
        raise RuntimeError("corrected-source conjugate-pair prerequisite missing")
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_a2_pr": PARENT_A2_PR,
        "parent_a2_head": PARENT_A2_HEAD,
        "parent_a2_source_blob": PARENT_A2_SOURCE_BLOB,
        "finite_complete_curl_harmonic_family_materialized": True,
        "same_coefficient_applied_to_vector_potential_and_curl_velocity": True,
        "real_conjugate_pair_sum_materialized": True,
        "linearity_preserves_complete_curl_divergence_free_contract": True,
        "per_provider_axis_zero_core_enforced": True,
        "bounded_term_count": True,
        "bounded_amplitude_l1": True,
        "autonomous_amplitude_phase_parameters": True,
        "source_exact_amplitude_phase_recovered": False,
        "orientation_parameter_materialized": False,
        "independent_scale_parameter_materialized": False,
        "caller_supplies_corrected_background_and_forcing": True,
        "caller_supplies_remaining_partition_and_pulse_cutoff": True,
        "source_exact_duhamel_frame_B_materialized": False,
        "source_exact_duhamel_propagator_Vm_materialized": False,
        "project_domain_source_input_provider_materialized": False,
        "self_contained_velocity_xyzt_provider": False,
        "complete_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def public_contract() -> dict[str, Any]:
    params = set(inspect.signature(KokunoBoundedMultiHarmonicPhysicalVelocity.velocity).parameters)
    forbidden = {
        "residual", "forcing", "pressure", "viscosity", "gain", "threshold",
        "target", "tolerance", "rtol", "atol", "panels", "steps", "optimizer",
    }
    return {
        "velocity_inputs": list(inspect.signature(KokunoBoundedMultiHarmonicPhysicalVelocity.velocity).parameters),
        "forbidden_velocity_inputs_present": sorted(params & forbidden),
        "max_terms": MAX_TERMS,
        "max_l1_amplitude": MAX_L1_AMPLITUDE,
        "paper_exact": False,
        "pde_validated": False,
    }


__all__ = [
    "BoundedHarmonicTerm",
    "BoundedMultiHarmonicEvaluation",
    "KokunoBoundedMultiHarmonicPhysicalVelocity",
    "bounded_multiharmonic_contract",
    "public_contract",
]
