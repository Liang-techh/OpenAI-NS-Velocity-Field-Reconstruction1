"""Executable source-shaped reference for Kokuno's normalized complete curl.

This module implements only algebra that is explicit in the corrected 2026-09-09
Kokuno reconstruction and checksum-bound by
``kokuno_corrected_oscillation_source_ledger``.  It is intentionally *not* a
mapping from the repository's autonomous oscillatory runtime to the corrected
source variables.

The public corrected source uses, on its annular chart,

    curl_* A = (
        R^-1 partial_theta A_z - D_z A_theta,
        D_z A_r - D_r A_z,
        (D_r + R^-1) A_theta - R^-1 partial_theta A_r,
    )

and for a transverse harmonic coefficient ``t_m``

    C_m = i (n_Phi x t_m) / (k m |n_Phi|^2),
    r_m = (-D_z C_theta,
           D_z C_r - D_r C_z,
           (D_r + R^-1) C_theta),
    a_m = t_m + r_m.

``D_r`` and ``D_z`` are the source normalized commuting directional
derivatives.  This reference consumes their already-evaluated jets rather than
pretending that the current repository runtime supplies the source auxiliary
coordinates/background map.  In particular, no source-to-runtime equivalence,
paper exactness, pressure, forcing, complete NS residual, or PDE validation is
claimed here.

The source waves live on annular support.  Therefore the normalized cylindrical
reference fails closed at ``R <= 0`` instead of inventing an axis extension.
The repository's separate physical Cartesian oscillatory runtime retains its own
axis-safe support handling.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from openai_ns_reconstruction.kokuno_corrected_oscillation_source_ledger import (
    SOURCE_COMMIT,
    SOURCE_REPOSITORY,
    SOURCE_WORKBENCH_BLOB,
    default_kokuno_corrected_oscillation_source_ledger,
)


PARENT_A2_HEAD = "bea6517b0f5a827bb2b35f05736e66a54177ed23"
PARENT_SOURCE_LEDGER_BLOB = "241ee6a9210aa2f4c8e3e9678e7b2175731c2b4a"
TRANSVERSALITY_ATOL = 5.0e-12
TRANSVERSALITY_RTOL = 5.0e-12
LONGITUDINAL_CLOSURE_ATOL = 2.0e-12
LONGITUDINAL_CLOSURE_RTOL = 2.0e-12


@dataclass(frozen=True)
class NormalizedPotentialJet:
    """First normalized directional jet needed by the source ``curl_*``.

    Every field may be scalar or an array.  NumPy broadcasting is applied by
    :func:`normalized_curl_from_potential_jet`.
    """

    a_r: Any
    a_theta: Any
    a_z: Any
    dtheta_a_r: Any
    dtheta_a_z: Any
    dr_a_z: Any
    dr_a_theta: Any
    dz_a_r: Any
    dz_a_theta: Any


@dataclass(frozen=True)
class CompleteHarmonicAmplitude:
    """Executable algebraic pieces of one corrected-source harmonic."""

    coefficient: np.ndarray
    remainder: np.ndarray
    amplitude: np.ndarray
    longitudinal_amplitude: np.ndarray
    longitudinal_remainder: np.ndarray



def _finite_array(value: Any, *, name: str, dtype: Any | None = None) -> np.ndarray:
    array = np.asarray(value, dtype=dtype)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array



def _positive_radius(radius: Any) -> np.ndarray:
    r = _finite_array(radius, name="radius", dtype=float)
    if np.any(r <= 0.0):
        raise ValueError("corrected-source normalized curl is annular and requires radius > 0")
    return r



def _vector_last_dim(value: Any, *, name: str, dtype: Any) -> np.ndarray:
    array = _finite_array(value, name=name, dtype=dtype)
    if array.ndim == 0 or array.shape[-1] != 3:
        raise ValueError(f"{name} must have final dimension 3")
    return array



def normalized_curl_from_potential_jet(
    radius: Any,
    jet: NormalizedPotentialJet,
) -> np.ndarray:
    """Evaluate the corrected-source normalized cylindrical ``curl_*``.

    ``dr_*`` and ``dz_*`` are values of the full source operators ``D_r`` and
    ``D_z``.  They are not silently replaced by repository Cartesian
    derivatives.  The returned shape is the broadcast scalar shape plus a final
    component dimension of length three.
    """

    r = _positive_radius(radius)
    values = [
        _finite_array(jet.a_r, name="a_r"),
        _finite_array(jet.a_theta, name="a_theta"),
        _finite_array(jet.a_z, name="a_z"),
        _finite_array(jet.dtheta_a_r, name="dtheta_a_r"),
        _finite_array(jet.dtheta_a_z, name="dtheta_a_z"),
        _finite_array(jet.dr_a_z, name="dr_a_z"),
        _finite_array(jet.dr_a_theta, name="dr_a_theta"),
        _finite_array(jet.dz_a_r, name="dz_a_r"),
        _finite_array(jet.dz_a_theta, name="dz_a_theta"),
    ]
    (
        r,
        a_r,
        a_theta,
        a_z,
        dtheta_a_r,
        dtheta_a_z,
        dr_a_z,
        dr_a_theta,
        dz_a_r,
        dz_a_theta,
    ) = np.broadcast_arrays(r, *values)

    del a_r, a_z  # values are part of the typed potential jet, but not used directly below.
    curl_r = dtheta_a_z / r - dz_a_theta
    curl_theta = dz_a_r - dr_a_z
    curl_z = dr_a_theta + a_theta / r - dtheta_a_r / r
    return np.stack((curl_r, curl_theta, curl_z), axis=-1)



def normalized_divergence_from_vector_jet(
    radius: Any,
    radial_component: Any,
    dr_radial_component: Any,
    dtheta_theta_component: Any,
    dz_axial_component: Any,
) -> np.ndarray:
    """Evaluate the source normalized divergence from a first vector jet."""

    r = _positive_radius(radius)
    radial = _finite_array(radial_component, name="radial_component")
    dr_radial = _finite_array(dr_radial_component, name="dr_radial_component")
    dtheta_theta = _finite_array(
        dtheta_theta_component, name="dtheta_theta_component"
    )
    dz_axial = _finite_array(dz_axial_component, name="dz_axial_component")
    r, radial, dr_radial, dtheta_theta, dz_axial = np.broadcast_arrays(
        r, radial, dr_radial, dtheta_theta, dz_axial
    )
    return dr_radial + radial / r + dtheta_theta / r + dz_axial



def harmonic_vector_potential_coefficient(
    n_phi: Any,
    t_m: Any,
    *,
    k: float,
    m: int,
) -> np.ndarray:
    """Return ``C_m = i(n_Phi x t_m)/(k m |n_Phi|^2)``.

    The corrected source requires ``n_Phi dot t_m = 0``.  This executable
    reference enforces that relation with fixed numerical tolerances; the
    tolerance is deliberately not exposed as a fitting parameter.
    """

    if not np.isfinite(k) or k <= 0.0:
        raise ValueError("k must be a finite positive carrier")
    if isinstance(m, bool) or not isinstance(m, (int, np.integer)) or int(m) == 0:
        raise ValueError("m must be a nonzero integer harmonic")
    m = int(m)

    n = _vector_last_dim(n_phi, name="n_phi", dtype=float)
    t = _vector_last_dim(t_m, name="t_m", dtype=np.complex128)
    n, t = np.broadcast_arrays(n, t)

    norm_sq = np.sum(n * n, axis=-1)
    if np.any(norm_sq <= 0.0) or not np.all(np.isfinite(norm_sq)):
        raise ValueError("n_phi must be nonzero")

    transverse = np.sum(n * t, axis=-1)
    scale = np.maximum(
        1.0,
        np.linalg.norm(n, axis=-1) * np.linalg.norm(t, axis=-1),
    )
    allowed = TRANSVERSALITY_ATOL + TRANSVERSALITY_RTOL * scale
    if np.any(np.abs(transverse) > allowed):
        raise ValueError("corrected-source harmonic requires n_phi dot t_m = 0")

    denominator = float(k * m) * norm_sq
    return 1j * np.cross(n, t) / denominator[..., None]



def complete_harmonic_amplitude_from_coefficient_jet(
    radius: Any,
    n_phi: Any,
    t_m: Any,
    coefficient_dr: Any,
    coefficient_dz: Any,
    *,
    k: float,
    m: int,
) -> CompleteHarmonicAmplitude:
    """Materialize the corrected source's complete ``a_m=t_m+r_m`` algebra.

    ``coefficient_dr`` and ``coefficient_dz`` are the full source normalized
    directional derivatives ``D_r C_m`` and ``D_z C_m``.  They must be supplied
    by a future source background/phase realization; this function does not
    infer them from the current repository oscillatory runtime.
    """

    r = _positive_radius(radius)
    n = _vector_last_dim(n_phi, name="n_phi", dtype=float)
    t = _vector_last_dim(t_m, name="t_m", dtype=np.complex128)
    dr_c = _vector_last_dim(
        coefficient_dr, name="coefficient_dr", dtype=np.complex128
    )
    dz_c = _vector_last_dim(
        coefficient_dz, name="coefficient_dz", dtype=np.complex128
    )
    c = harmonic_vector_potential_coefficient(n, t, k=k, m=m)

    scalar_shape = np.broadcast_shapes(
        np.shape(r),
        np.shape(n)[:-1],
        np.shape(t)[:-1],
        np.shape(dr_c)[:-1],
        np.shape(dz_c)[:-1],
        np.shape(c)[:-1],
    )
    r = np.broadcast_to(r, scalar_shape)
    n = np.broadcast_to(n, scalar_shape + (3,))
    t = np.broadcast_to(t, scalar_shape + (3,))
    dr_c = np.broadcast_to(dr_c, scalar_shape + (3,))
    dz_c = np.broadcast_to(dz_c, scalar_shape + (3,))
    c = np.broadcast_to(c, scalar_shape + (3,))

    remainder = np.stack(
        (
            -dz_c[..., 1],
            dz_c[..., 0] - dr_c[..., 2],
            dr_c[..., 1] + c[..., 1] / r,
        ),
        axis=-1,
    )
    amplitude = t + remainder
    longitudinal_amplitude = np.sum(n * amplitude, axis=-1)
    longitudinal_remainder = np.sum(n * remainder, axis=-1)

    scale = np.maximum(
        1.0,
        np.maximum(np.abs(longitudinal_amplitude), np.abs(longitudinal_remainder)),
    )
    allowed = LONGITUDINAL_CLOSURE_ATOL + LONGITUDINAL_CLOSURE_RTOL * scale
    if np.any(np.abs(longitudinal_amplitude - longitudinal_remainder) > allowed):
        raise RuntimeError("retained longitudinal identity failed")

    return CompleteHarmonicAmplitude(
        coefficient=c,
        remainder=remainder,
        amplitude=amplitude,
        longitudinal_amplitude=longitudinal_amplitude,
        longitudinal_remainder=longitudinal_remainder,
    )



def source_reference_contract() -> dict[str, Any]:
    """Return a fail-closed provenance/truth contract for this reference."""

    ledger = default_kokuno_corrected_oscillation_source_ledger()
    if ledger.normalized_curl != (
        "(curl_* A)_r = R^(-1) partial_theta A_z - D_z A_theta",
        "(curl_* A)_theta = D_z A_r - D_r A_z",
        "(curl_* A)_z = (D_r + R^(-1)) A_theta - R^(-1) partial_theta A_r",
    ):
        raise RuntimeError("corrected-source curl formula drifted from reference")
    if ledger.exact_curl_remainder != (
        "(r_m)_r = -D_z (C_m)_theta",
        "(r_m)_theta = D_z (C_m)_r - D_r (C_m)_z",
        "(r_m)_z = (D_r + R^(-1)) (C_m)_theta",
    ):
        raise RuntimeError("corrected-source complete-curl remainder drifted from reference")

    return {
        "schema": "kokuno-source-normalized-complete-curl-reference-v1",
        "parent_a2_head": PARENT_A2_HEAD,
        "parent_source_ledger_blob": PARENT_SOURCE_LEDGER_BLOB,
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": SOURCE_COMMIT,
        "source_workbench_blob": SOURCE_WORKBENCH_BLOB,
        "source_normalized_curl_formula_executable": True,
        "source_complete_harmonic_amplitude_algebra_executable": True,
        "source_annular_domain_requires_positive_radius": True,
        "source_directional_derivatives_supplied_by_caller": True,
        "current_runtime_phase_source_exact": False,
        "current_runtime_support_source_exact": False,
        "current_runtime_complete_curl_source_equivalence_verified": False,
        "source_to_runtime_parameter_map_complete": False,
        "paper_exact": False,
        "matched_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "complete_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
    }
