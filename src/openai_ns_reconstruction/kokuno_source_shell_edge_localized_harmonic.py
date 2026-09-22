"""Materialize the corrected-source shell edge inside A2 complete-curl localization.

This Kokuno Agent-2 increment consumes one explicit support factor from the
corrected 2026-09-09 reconstruction without inventing the still-unresolved
partition or pulse cutoff.  The source defines, on one shell,

    d_a = log(X/X_a),
    d_b = log(X_b/X),
    zeta = exp(-a_a/d_a**2 - a_b/d_b**2),  a_a,a_b > 0,

with smooth zero extension outside ``X_a < X < X_b``.  The oscillatory source
bounds carry the edge factor ``sqrt(zeta)``.  This module evaluates that factor
and its analytic X derivative, pulls the derivative back through caller-supplied
normalized ``D_r X`` / ``D_z X``, and multiplies it into the remaining source
partition jet before the existing #1089 potential-level complete-curl path.

The derivative is

    d_X sqrt(zeta)
      = sqrt(zeta) / X * (a_a/d_a**3 - a_b/d_b**3).

At and outside the shell endpoints both the value and directional jets are
returned exactly as zero, matching the source's smooth zero extension.  The
implementation also avoids ``0 * inf`` at machine-underflow distances from the
edge by evaluating the derivative only where the exponential remains
representable.

Important truth boundary: ``X`` and its normalized directional jets are still
upstream source-chart data.  The concrete squared partition
``chi_l(q) chi_{l,a}``, the pulse cutoff ``psi(v)``, corrected background and
forcing providers, physical Q scaling, and the project-domain coordinate map
are not reconstructed here.  This is therefore an executable source-formula
adapter, not a self-contained ``velocity(x,y,z,t)`` provider and not a
paper-exact/OpenAI-field claim.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .kokuno_corrected_oscillation_source_ledger import (
    SOURCE_COMMIT,
    SOURCE_REPOSITORY,
    SOURCE_WORKBENCH_BLOB,
    default_kokuno_corrected_oscillation_source_ledger,
)
from .kokuno_source_phase_coefficient_jet import SourceBackgroundPhaseJet
from .kokuno_source_zero_data_amplitude_sensitivity import ModeForcingDirectionalJet
from .kokuno_source_zero_data_localized_harmonic import (
    SourceLocalizationJet,
    SourceLocalizedZeroDataHarmonic,
    materialize_source_zero_data_localized_harmonic,
    source_zero_data_localized_harmonic_contract,
)


TASK = "K2-OSC-099"
PARENT_A2_HEAD = "4bc9408d10b2f3c636abd4f6ce33cefab391093d"
PARENT_LOCALIZED_HARMONIC_BLOB = "8087e619bb1812df6dcb616e1a51f1ab4dd728d4"
SOURCE_SHELL_FORMULA = (
    "d_a=log(X/X_a); d_b=log(X_b/X); "
    "zeta=exp(-a_a/d_a^2-a_b/d_b^2); use sqrt(zeta) edge weight"
)


@dataclass(frozen=True)
class SourceShellCoordinateJet:
    """Source shell coordinate and its normalized spatial directional jets."""

    x: Any
    dr_x: Any
    dz_x: Any


@dataclass(frozen=True)
class SourceShellEdgeJet:
    """Executable ``sqrt(zeta)`` shell edge and first normalized jets."""

    weight: np.ndarray
    dx_weight: np.ndarray
    dr_weight: np.ndarray
    dz_weight: np.ndarray
    inside_shell: np.ndarray


@dataclass(frozen=True)
class SourceShellLocalizedZeroDataHarmonic:
    """#1089 localized harmonic with the explicit source shell edge consumed."""

    shell_edge: SourceShellEdgeJet
    effective_localization: SourceLocalizationJet
    harmonic: SourceLocalizedZeroDataHarmonic


def _finite_real(value: Any, *, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _positive_scalar(value: float, *, name: str) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def source_shell_edge_jet(
    coordinate: SourceShellCoordinateJet,
    *,
    x_a: float,
    x_b: float,
    a_a: float,
    a_b: float,
) -> SourceShellEdgeJet:
    """Evaluate the corrected-source ``sqrt(zeta)`` edge and analytic jets.

    ``D_r X`` and ``D_z X`` are deliberately explicit inputs: this function
    does not invent the still-unresolved source-chart -> project-domain map.
    """

    if not isinstance(coordinate, SourceShellCoordinateJet):
        raise ValueError("coordinate must be SourceShellCoordinateJet")
    x_a = _positive_scalar(x_a, name="x_a")
    x_b = _positive_scalar(x_b, name="x_b")
    a_a = _positive_scalar(a_a, name="a_a")
    a_b = _positive_scalar(a_b, name="a_b")
    if not x_b > x_a:
        raise ValueError("x_b must be strictly greater than x_a")

    x = _finite_real(coordinate.x, name="coordinate.x")
    dr_x = _finite_real(coordinate.dr_x, name="coordinate.dr_x")
    dz_x = _finite_real(coordinate.dz_x, name="coordinate.dz_x")
    x, dr_x, dz_x = np.broadcast_arrays(x, dr_x, dz_x)
    if np.any(x <= 0.0):
        raise ValueError("source shell coordinate X must be positive")

    weight = np.zeros(x.shape, dtype=float)
    dx_weight = np.zeros(x.shape, dtype=float)
    inside = (x > x_a) & (x < x_b)

    if np.any(inside):
        x_i = x[inside]
        d_a = np.log(x_i / x_a)
        d_b = np.log(x_b / x_i)
        log_weight = -0.5 * (a_a / d_a**2 + a_b / d_b**2)

        # exp(log_weight) is smoothly negligible before the subnormal floor.
        # Avoid forming an enormous logarithmic derivative when the value has
        # already underflowed to exact zero.
        log_min_subnormal = float(np.log(np.nextafter(0.0, 1.0)))
        representable = log_weight > log_min_subnormal
        local_weight = np.zeros_like(x_i)
        local_dx = np.zeros_like(x_i)
        if np.any(representable):
            d_a_r = d_a[representable]
            d_b_r = d_b[representable]
            x_r = x_i[representable]
            w_r = np.exp(log_weight[representable])
            dlog_dx = (a_a / d_a_r**3 - a_b / d_b_r**3) / x_r
            local_weight[representable] = w_r
            local_dx[representable] = w_r * dlog_dx
        weight[inside] = local_weight
        dx_weight[inside] = local_dx

    dr_weight = dx_weight * dr_x
    dz_weight = dx_weight * dz_x
    if not (
        np.all(np.isfinite(weight))
        and np.all(np.isfinite(dx_weight))
        and np.all(np.isfinite(dr_weight))
        and np.all(np.isfinite(dz_weight))
    ):
        raise RuntimeError("source shell edge evaluation produced nonfinite output")

    return SourceShellEdgeJet(
        weight=weight,
        dx_weight=dx_weight,
        dr_weight=dr_weight,
        dz_weight=dz_weight,
        inside_shell=inside,
    )


def apply_source_shell_edge_to_localization(
    localization: SourceLocalizationJet,
    shell_edge: SourceShellEdgeJet,
) -> SourceLocalizationJet:
    """Multiply ``sqrt(zeta)`` into the remaining source partition jet."""

    if not isinstance(localization, SourceLocalizationJet):
        raise ValueError("localization must be SourceLocalizationJet")
    if not isinstance(shell_edge, SourceShellEdgeJet):
        raise ValueError("shell_edge must be SourceShellEdgeJet")

    eta = _finite_real(localization.eta, name="localization.eta")
    dr_eta = _finite_real(localization.dr_eta, name="localization.dr_eta")
    dz_eta = _finite_real(localization.dz_eta, name="localization.dz_eta")
    shell = _finite_real(shell_edge.weight, name="shell_edge.weight")
    dr_shell = _finite_real(shell_edge.dr_weight, name="shell_edge.dr_weight")
    dz_shell = _finite_real(shell_edge.dz_weight, name="shell_edge.dz_weight")
    eta, dr_eta, dz_eta, shell, dr_shell, dz_shell = np.broadcast_arrays(
        eta, dr_eta, dz_eta, shell, dr_shell, dz_shell
    )

    return SourceLocalizationJet(
        eta=eta * shell,
        dr_eta=dr_eta * shell + eta * dr_shell,
        dz_eta=dz_eta * shell + eta * dz_shell,
        pulse_cutoff=localization.pulse_cutoff,
    )


def materialize_source_zero_data_shell_localized_harmonic(
    pulse_v: float,
    radius: Any,
    theta: Any,
    z_normalized: Any,
    background: SourceBackgroundPhaseJet,
    forcing_jet_m: ModeForcingDirectionalJet,
    localization_without_shell_edge: SourceLocalizationJet,
    shell_coordinate: SourceShellCoordinateJet,
    *,
    x_a: float,
    x_b: float,
    a_a: float,
    a_b: float,
    epsilon: float,
    p: float,
    p_z: float,
    x_0: float,
    k: float,
    m: int,
) -> SourceShellLocalizedZeroDataHarmonic:
    """Consume the public shell edge and execute the existing complete curl."""

    shell = source_shell_edge_jet(
        shell_coordinate,
        x_a=x_a,
        x_b=x_b,
        a_a=a_a,
        a_b=a_b,
    )
    effective = apply_source_shell_edge_to_localization(
        localization_without_shell_edge,
        shell,
    )
    harmonic = materialize_source_zero_data_localized_harmonic(
        pulse_v,
        radius,
        theta,
        z_normalized,
        background,
        forcing_jet_m,
        effective,
        epsilon=epsilon,
        p=p,
        p_z=p_z,
        x_0=x_0,
        k=k,
        m=m,
    )
    return SourceShellLocalizedZeroDataHarmonic(
        shell_edge=shell,
        effective_localization=effective,
        harmonic=harmonic,
    )


def source_shell_edge_localized_harmonic_contract() -> dict[str, Any]:
    """Return a machine-readable source/realization boundary for this increment."""

    ledger = default_kokuno_corrected_oscillation_source_ledger()
    parent = source_zero_data_localized_harmonic_contract()
    if "zeta = exp(-a_a/log^2(X/X_a)" not in ledger.source_support_envelope:
        raise RuntimeError("corrected-source shell envelope statement drifted")
    if "sqrt(zeta) edge weight" not in ledger.source_wave_support:
        raise RuntimeError("corrected-source sqrt(zeta) support statement drifted")
    if parent["potential_level_product_rule_localization_executable"] is not True:
        raise RuntimeError("required parent potential-level localization is unavailable")

    return {
        "schema": "kokuno-source-shell-edge-localized-harmonic-v1",
        "task": TASK,
        "parent_a2_head": PARENT_A2_HEAD,
        "parent_localized_harmonic_blob": PARENT_LOCALIZED_HARMONIC_BLOB,
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": SOURCE_COMMIT,
        "source_workbench_blob": SOURCE_WORKBENCH_BLOB,
        "source_shell_formula": SOURCE_SHELL_FORMULA,
        "source_sqrt_zeta_shell_edge_materialized": True,
        "source_shell_analytic_x_derivative_materialized": True,
        "source_shell_normalized_directional_pullback_materialized": True,
        "source_shell_smooth_zero_extension_materialized": True,
        "source_shell_edge_consumed_before_complete_curl": True,
        "caller_supplies_shell_coordinate_directional_jet": True,
        "caller_supplies_remaining_partition_jet": True,
        "caller_supplies_pulse_cutoff_value": True,
        "concrete_source_chi_partition_bumps_reconstructed": False,
        "source_pulse_cutoff_provider_materialized": False,
        "source_forcing_provider_materialized": False,
        "actual_corrected_background_provider_materialized": False,
        "source_physical_Q_scaling_applied": False,
        "project_domain_coordinate_map_applied": False,
        "global_axis_safe_source_velocity_materialized": False,
        "self_contained_velocity_xyzt_provider": False,
        "complete_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "paper_exact": False,
        "pde_validated": False,
    }
