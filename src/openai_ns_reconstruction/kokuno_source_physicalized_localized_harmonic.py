"""Physical scaling for the corrected-source localized complete-curl harmonic.

This Kokuno Agent-2 increment consumes the already-materialized source-chart
localized complete curl and applies only the corrected public band scaling.
It does not choose background/forcing/support data and is not a self-contained
project-domain velocity provider.

Corrected-source formulas used here are

    A = 1/2 + h,  D = 1/2 - h,  Q = 2**(-ell),  epsilon = Q**h,
    R = r / Q**(1/2),  Z = z / Q**D,  T = (1-t) / Q,
    A_phys = Q**(1/2-A) A_*,
    u_* = Q**A u_phys.

Hence u_phys = Q**(-A) u_*.  The equality between the physical curl of
A_phys and Q**(-A) times the normalized complete curl is a repository
algebraic consequence of the source scaling: both physical r and z curl
derivatives contribute Q**(-1/2), because D_z = epsilon partial_Z and
Q**(-D) = Q**(-1/2) epsilon.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any
import math

import numpy as np

from .kokuno_corrected_oscillation_source_ledger import (
    SOURCE_COMMIT,
    SOURCE_REPOSITORY,
    SOURCE_WORKBENCH_BLOB,
    default_kokuno_corrected_oscillation_source_ledger,
)
from .kokuno_source_zero_data_localized_harmonic import SourceLocalizedZeroDataHarmonic


TASK = "K2-OSC-102"
PARENT_A2_PR = 1117
PARENT_A2_HEAD = "27741d9c0a27262f7fabf61eebaa2fbd507e9f03"
SOURCE_VELOCITY_SCALING = "u_* = Q^A u_phys"
SOURCE_VECTOR_POTENTIAL_SCALING = "A_phys = Q^(1/2-A) A_*"


@dataclass(frozen=True)
class SourcePhysicalScaling:
    ell: int
    h: float
    Q: float
    A: float
    D: float
    epsilon: float
    radial_scale: float
    axial_scale: float
    time_scale: float
    vector_potential_scale: float
    velocity_scale: float
    curl_derivative_scale: float


@dataclass(frozen=True)
class SourcePhysicalizedLocalizedHarmonic:
    scaling: SourcePhysicalScaling
    source_harmonic: SourceLocalizedZeroDataHarmonic
    complex_vector_potential_physical: np.ndarray
    complex_velocity_cylindrical_physical: np.ndarray
    real_pair_velocity_cylindrical_physical: np.ndarray
    real_pair_velocity_cartesian_physical: np.ndarray


def _finite_array(value: Any, *, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _positive_integer(value: Any, *, name: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be a positive integer")
    out = int(value)
    if out <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return out


def source_physical_scaling(*, ell: int, h: float) -> SourcePhysicalScaling:
    """Return the binary64 realization of the corrected band scaling.

    The public source requires ``0<h<1/100`` and ``Q=2**(-ell)``.  We fail
    closed when that exact dyadic Q or one of the required physical scaling
    factors is not representable as a finite positive binary64 number.
    """

    ell = _positive_integer(ell, name="ell")
    h = float(h)
    if not math.isfinite(h) or not (0.0 < h < 0.01):
        raise ValueError("h must satisfy the corrected-source range 0 < h < 1/100")
    Q = math.ldexp(1.0, -ell)
    if Q <= 0.0 or not math.isfinite(Q):
        raise ValueError("Q=2^(-ell) is outside this binary64 realization")

    A = 0.5 + h
    D = 0.5 - h
    epsilon = Q**h
    radial_scale = math.sqrt(Q)
    axial_scale = Q**D
    time_scale = Q
    vector_potential_scale = Q ** (0.5 - A)
    velocity_scale = Q ** (-A)
    curl_derivative_scale = Q ** (-0.5)
    values = (
        epsilon,
        radial_scale,
        axial_scale,
        time_scale,
        vector_potential_scale,
        velocity_scale,
        curl_derivative_scale,
    )
    if not all(math.isfinite(x) and x > 0.0 for x in values):
        raise ValueError("corrected-source physical scaling is not finite in binary64")
    if not math.isclose(
        vector_potential_scale * curl_derivative_scale,
        velocity_scale,
        rel_tol=4.0e-15,
        abs_tol=0.0,
    ):
        raise RuntimeError("physical vector-potential/curl scaling identity drifted")
    if not math.isclose(Q ** (-D), curl_derivative_scale * epsilon, rel_tol=4.0e-15):
        raise RuntimeError("normalized D_z / physical z derivative scaling drifted")

    return SourcePhysicalScaling(
        ell=ell,
        h=h,
        Q=Q,
        A=A,
        D=D,
        epsilon=epsilon,
        radial_scale=radial_scale,
        axial_scale=axial_scale,
        time_scale=time_scale,
        vector_potential_scale=vector_potential_scale,
        velocity_scale=velocity_scale,
        curl_derivative_scale=curl_derivative_scale,
    )


def source_chart_to_physical_rzt(
    radius_chart: Any,
    z_chart: Any,
    T_chart: Any,
    *,
    ell: int,
    h: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map corrected source chart ``(R,Z,T)`` to physical ``(r,z,t)``."""

    scaling = source_physical_scaling(ell=ell, h=h)
    R = _finite_array(radius_chart, name="radius_chart")
    Z = _finite_array(z_chart, name="z_chart")
    T = _finite_array(T_chart, name="T_chart")
    R, Z, T = np.broadcast_arrays(R, Z, T)
    r = scaling.radial_scale * R
    z = scaling.axial_scale * Z
    t = 1.0 - scaling.time_scale * T
    if not (np.all(np.isfinite(r)) and np.all(np.isfinite(z)) and np.all(np.isfinite(t))):
        raise ValueError("source-chart to physical map produced nonfinite output")
    return r, z, t


def physical_to_source_chart_rzt(
    radius_physical: Any,
    z_physical: Any,
    t_physical: Any,
    *,
    ell: int,
    h: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Invert the corrected source chart map on finite binary64 inputs."""

    scaling = source_physical_scaling(ell=ell, h=h)
    r = _finite_array(radius_physical, name="radius_physical")
    z = _finite_array(z_physical, name="z_physical")
    t = _finite_array(t_physical, name="t_physical")
    r, z, t = np.broadcast_arrays(r, z, t)
    R = r / scaling.radial_scale
    Z = z / scaling.axial_scale
    T = (1.0 - t) / scaling.time_scale
    if not (np.all(np.isfinite(R)) and np.all(np.isfinite(Z)) and np.all(np.isfinite(T))):
        raise ValueError("physical to source-chart map produced nonfinite output")
    return R, Z, T


def physicalize_source_localized_harmonic(
    harmonic: SourceLocalizedZeroDataHarmonic,
    *,
    ell: int,
    h: float,
) -> SourcePhysicalizedLocalizedHarmonic:
    """Apply corrected physical A/u scaling to one localized chart harmonic."""

    if not isinstance(harmonic, SourceLocalizedZeroDataHarmonic):
        raise ValueError("harmonic must be SourceLocalizedZeroDataHarmonic")
    scaling = source_physical_scaling(ell=ell, h=h)

    potential = np.asarray(harmonic.complex_vector_potential, dtype=np.complex128)
    complex_velocity = np.asarray(harmonic.complex_velocity_cylindrical, dtype=np.complex128)
    real_cylindrical = np.asarray(harmonic.real_pair_velocity_cylindrical, dtype=float)
    real_cartesian = np.asarray(harmonic.real_pair_velocity_chart_cartesian, dtype=float)
    for name, array in (
        ("complex_vector_potential", potential),
        ("complex_velocity_cylindrical", complex_velocity),
        ("real_pair_velocity_cylindrical", real_cylindrical),
        ("real_pair_velocity_chart_cartesian", real_cartesian),
    ):
        if array.ndim == 0 or array.shape[-1] != 3 or not np.all(np.isfinite(array)):
            raise ValueError(f"harmonic.{name} must be a finite 3-vector field")

    physical_potential = scaling.vector_potential_scale * potential
    physical_complex_velocity = scaling.velocity_scale * complex_velocity
    physical_real_cylindrical = scaling.velocity_scale * real_cylindrical
    physical_real_cartesian = scaling.velocity_scale * real_cartesian
    arrays = (
        physical_potential,
        physical_complex_velocity,
        physical_real_cylindrical,
        physical_real_cartesian,
    )
    if not all(np.all(np.isfinite(array)) for array in arrays):
        raise ValueError("physicalized harmonic is not finite in binary64")

    return SourcePhysicalizedLocalizedHarmonic(
        scaling=scaling,
        source_harmonic=harmonic,
        complex_vector_potential_physical=physical_potential,
        complex_velocity_cylindrical_physical=physical_complex_velocity,
        real_pair_velocity_cylindrical_physical=physical_real_cylindrical,
        real_pair_velocity_cartesian_physical=physical_real_cartesian,
    )


def source_physicalized_localized_harmonic_contract() -> dict[str, Any]:
    """Return the exact-source / repository-derivation / pending boundary."""

    ledger = default_kokuno_corrected_oscillation_source_ledger()
    if ledger.physical_vector_potential_scaling != "A_phys = Q^(1/2-A) A_*":
        raise RuntimeError("corrected-source physical vector-potential scaling drifted")
    required = {
        "A = 1/2 + h",
        "D = 1/2 - h",
        "Q = 2^(-ell)",
        "epsilon = Q^h",
        "R = r / Q^(1/2)",
        "Z = z / Q^D",
        "T = (1-t) / Q",
    }
    if not required.issubset(set(ledger.chart_scaling)):
        raise RuntimeError("corrected-source chart scaling ledger drifted")

    return {
        "schema": "kokuno-source-physicalized-localized-harmonic-v1",
        "task": TASK,
        "parent_a2_pr": PARENT_A2_PR,
        "parent_a2_head": PARENT_A2_HEAD,
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": SOURCE_COMMIT,
        "source_workbench_blob": SOURCE_WORKBENCH_BLOB,
        "source_formula_vector_potential_scaling": SOURCE_VECTOR_POTENTIAL_SCALING,
        "source_formula_velocity_scaling": SOURCE_VELOCITY_SCALING,
        "source_chart_rzt_map_materialized": True,
        "source_physical_Q_scaling_applied": True,
        "repository_curl_scaling_identity_verified": True,
        "physical_real_pair_velocity_materialized": True,
        "caller_supplies_mode_forcing_directional_jet": True,
        "caller_supplies_corrected_background": True,
        "caller_supplies_remaining_partition_and_pulse_cutoff": True,
        "source_exact_duhamel_frame_B_materialized": False,
        "source_exact_duhamel_propagator_Vm_materialized": False,
        "project_domain_source_input_provider_materialized": False,
        "global_axis_safe_source_velocity_materialized": False,
        "self_contained_velocity_xyzt_provider": False,
        "complete_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "paper_exact": False,
        "pde_validated": False,
    }
