"""Typed corrected-source background-jet adapter for the Kokuno oscillation lane.

The corrected source imports a background angular component ``V`` and axial
component ``G`` and then sets

    F = V/R,
    H_Phi = p F + p_z G.

The phase/coefficient-jet implementation introduced by A2 #1021 consumes
``F, G`` and their ordinary first/second ``R/Z`` derivatives.  This module
removes that caller-supplied algebra layer: callers provide the physical
background ``V,G`` jet, while ``F=V/R`` and every derivative required by
``SourceBackgroundPhaseJet`` are derived analytically.

This is deliberately *not* a corrected-background provider.  The actual
source background construction is still upstream and unresolved here.  In
particular, Agent-1 leading profiles are not substituted for ``V,G`` and this
adapter does not claim paper-exact phase/support/runtime equivalence.
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
from openai_ns_reconstruction.kokuno_source_phase_coefficient_jet import (
    SourceBackgroundPhaseJet,
    SourcePhaseVectorJet,
    source_phase_vector_jet,
)


PARENT_A2_HEAD = "0f4006bfd159281d96a7dc90a3ab50f0566f13e1"
PARENT_PHASE_COEFFICIENT_JET_BLOB = "d3c9684d47cec5dfabf713db01ffab45ea9461c4"


@dataclass(frozen=True)
class SourceActualBackgroundVGJet:
    """Actual-background ``V,G`` values and ordinary ``R/Z`` derivatives.

    ``V`` is the corrected source's angular/tangential background component
    used through ``F=V/R``.  These values remain output-determining source
    realization inputs until an upstream corrected-background provider is
    materialized.
    """

    v: Any
    g: Any
    v_r: Any
    g_r: Any
    v_z: Any
    g_z: Any
    v_rr: Any
    g_rr: Any
    v_rz: Any
    g_rz: Any
    v_zz: Any
    g_zz: Any


@dataclass(frozen=True)
class SourceBackgroundPhaseAdapterResult:
    """Derived ``F,G`` phase input plus the executable #1021 phase-vector jet."""

    background_phase_jet: SourceBackgroundPhaseJet
    phase_vector_jet: SourcePhaseVectorJet


def _finite(value: Any, *, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _positive_radius(radius: Any) -> np.ndarray:
    array = _finite(radius, name="radius")
    if np.any(array <= 0.0):
        raise ValueError("corrected-source background chart requires radius > 0")
    return array


def source_background_phase_jet_from_vg(
    radius: Any,
    actual_background: SourceActualBackgroundVGJet,
) -> SourceBackgroundPhaseJet:
    """Derive the #1021 ``F,G`` background jet from source ``V,G`` data.

    All derivatives here are ordinary ``partial_R`` / ``partial_Z`` derivatives.
    With ``F=V/R``,

    ``F_R  = V_R/R - V/R^2``
    ``F_Z  = V_Z/R``
    ``F_RR = V_RR/R - 2 V_R/R^2 + 2 V/R^3``
    ``F_RZ = V_RZ/R - V_Z/R^2``
    ``F_ZZ = V_ZZ/R``.

    No finite-difference step, fitting coefficient, or free ``H_Phi`` jet is
    introduced.
    """

    radius = _positive_radius(radius)
    values = [
        _finite(actual_background.v, name="actual_background.v"),
        _finite(actual_background.g, name="actual_background.g"),
        _finite(actual_background.v_r, name="actual_background.v_r"),
        _finite(actual_background.g_r, name="actual_background.g_r"),
        _finite(actual_background.v_z, name="actual_background.v_z"),
        _finite(actual_background.g_z, name="actual_background.g_z"),
        _finite(actual_background.v_rr, name="actual_background.v_rr"),
        _finite(actual_background.g_rr, name="actual_background.g_rr"),
        _finite(actual_background.v_rz, name="actual_background.v_rz"),
        _finite(actual_background.g_rz, name="actual_background.g_rz"),
        _finite(actual_background.v_zz, name="actual_background.v_zz"),
        _finite(actual_background.g_zz, name="actual_background.g_zz"),
    ]
    (
        radius,
        v,
        g,
        v_r,
        g_r,
        v_z,
        g_z,
        v_rr,
        g_rr,
        v_rz,
        g_rz,
        v_zz,
        g_zz,
    ) = np.broadcast_arrays(radius, *values)

    inv_r = 1.0 / radius
    inv_r2 = inv_r * inv_r
    inv_r3 = inv_r2 * inv_r

    f = v * inv_r
    f_r = v_r * inv_r - v * inv_r2
    f_z = v_z * inv_r
    f_rr = v_rr * inv_r - 2.0 * v_r * inv_r2 + 2.0 * v * inv_r3
    f_rz = v_rz * inv_r - v_z * inv_r2
    f_zz = v_zz * inv_r

    return SourceBackgroundPhaseJet(
        f=f,
        g=g,
        f_r=f_r,
        g_r=g_r,
        f_z=f_z,
        g_z=g_z,
        f_rr=f_rr,
        g_rr=g_rr,
        f_rz=f_rz,
        g_rz=g_rz,
        f_zz=f_zz,
        g_zz=g_zz,
    )


def source_phase_vector_jet_from_vg(
    radius: Any,
    theta: Any,
    z_normalized: Any,
    pulse_v: Any,
    actual_background: SourceActualBackgroundVGJet,
    *,
    epsilon: float,
    p: float,
    p_z: float,
    x_0: float,
) -> SourceBackgroundPhaseAdapterResult:
    """Derive source ``F=V/R`` inputs and execute the #1021 phase-vector jet."""

    background_phase_jet = source_background_phase_jet_from_vg(
        radius, actual_background
    )
    phase_vector_jet = source_phase_vector_jet(
        radius,
        theta,
        z_normalized,
        pulse_v,
        background_phase_jet,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
    )
    return SourceBackgroundPhaseAdapterResult(
        background_phase_jet=background_phase_jet,
        phase_vector_jet=phase_vector_jet,
    )


def source_background_phase_jet_contract() -> dict[str, Any]:
    """Fail-closed provenance and truth boundary for this algebra adapter."""

    ledger = default_kokuno_corrected_oscillation_source_ledger()
    if "F = V/R" not in ledger.required_background_inputs:
        raise RuntimeError("corrected-source F=V/R requirement drifted")
    if ledger.phase_hamiltonian != "H_Phi = p F + p_z G":
        raise RuntimeError("corrected-source H_Phi formula drifted")

    return {
        "schema": "kokuno-source-background-phase-jet-v1",
        "parent_a2_head": PARENT_A2_HEAD,
        "parent_phase_coefficient_jet_blob": PARENT_PHASE_COEFFICIENT_JET_BLOB,
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": SOURCE_COMMIT,
        "source_workbench_blob": SOURCE_WORKBENCH_BLOB,
        "source_f_equals_v_over_r_executed": True,
        "source_f_first_second_derivatives_derived_analytically": True,
        "source_h_phi_inputs_feed_parent_phase_vector_jet": True,
        "caller_supplies_background_phase_derivative_jet": False,
        "caller_supplies_actual_background_vg_derivative_jet": True,
        "actual_corrected_background_provider_materialized": False,
        "agent1_leading_substituted_for_corrected_background": False,
        "source_amplitude_ode_materialized": False,
        "source_support_cutoff_runtime_materialized": False,
        "self_contained_velocity_xyzt_provider": False,
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
