"""Axis-safe physical Cartesian evaluation of a source-localized Kokuno harmonic.

This Kokuno Agent-2 increment adds only the physical evaluation seam on top of
K2-OSC-102. The actual corrected background, mode forcing, slow partition and
pulse cutoff remain upstream data. A typed provider supplies the already
localized source-chart harmonic at requested off-axis chart points.

The cylindrical complete-curl formulas contain ``1/R`` terms, so the exact
axis must not be reached by inventing a numerical angle or by clipping R. The
provider therefore carries an explicit *repository execution certificate*:
its localized harmonic is identically zero for

    0 <= R <= axis_zero_radius_chart.

Points in that certified zero core are returned as exact Cartesian zero and
the provider is not called there. This radius is not recovered from Kokuno's
paper/source and is not promoted to a source-exact support statement.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable
import math
import re

import numpy as np

from .kokuno_corrected_oscillation_source_ledger import (
    SOURCE_COMMIT,
    SOURCE_REPOSITORY,
    SOURCE_WORKBENCH_BLOB,
)
from .kokuno_source_physicalized_localized_harmonic import (
    SourcePhysicalScaling,
    physical_to_source_chart_rzt,
    physicalize_source_localized_harmonic,
    source_physical_scaling,
    source_physicalized_localized_harmonic_contract,
)
from .kokuno_source_zero_data_localized_harmonic import (
    SourceLocalizedZeroDataHarmonic,
)


TASK = "K2-OSC-103"
PARENT_A2_PR = 1125
PARENT_A2_HEAD = "bab607c1a4f93a586dc8a1a9fa3c0d8b5badeb58"
PARENT_A2_SOURCE_BLOB = "b19805fc2e78d6a51b6ead6c141125a662d3ec4c"
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


ChartHarmonicEvaluator = Callable[
    [np.ndarray, np.ndarray, np.ndarray, np.ndarray, SourcePhysicalScaling],
    SourceLocalizedZeroDataHarmonic,
]


@dataclass(frozen=True)
class AxisSafeSourceHarmonicProvider:
    """Identity-bound caller provider for one localized source harmonic.

    ``axis_zero_radius_chart`` and ``axis_zero_certified`` are a repository
    execution contract, not a public-source numerical datum. The evaluator is
    called only for chart radii strictly larger than the certified zero core.

    The semantic SHA is recorded to make the caller-selected source-input
    realization explicit; this adapter does not claim that SHA identifies a
    paper-exact source field.
    """

    provider_id: str
    provider_semantic_sha256: str
    axis_zero_radius_chart: float
    axis_zero_certified: bool
    evaluate_chart: ChartHarmonicEvaluator


def _finite_array(value: Any, *, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _validated_provider(provider: AxisSafeSourceHarmonicProvider) -> AxisSafeSourceHarmonicProvider:
    if not isinstance(provider, AxisSafeSourceHarmonicProvider):
        raise ValueError("provider must be AxisSafeSourceHarmonicProvider")
    if not isinstance(provider.provider_id, str) or not provider.provider_id.strip():
        raise ValueError("provider_id must be a nonempty string")
    if not isinstance(provider.provider_semantic_sha256, str) or not _SHA256_RE.fullmatch(
        provider.provider_semantic_sha256
    ):
        raise ValueError("provider_semantic_sha256 must be a lowercase 64-hex digest")
    radius = float(provider.axis_zero_radius_chart)
    if not math.isfinite(radius) or radius <= 0.0:
        raise ValueError("axis_zero_radius_chart must be finite and positive")
    if provider.axis_zero_certified is not True:
        raise ValueError("provider must certify exact zero on its declared axis core")
    if not callable(provider.evaluate_chart):
        raise ValueError("provider.evaluate_chart must be callable")
    return provider


class KokunoProviderDrivenAxisSafePhysicalVelocity:
    """Provider-driven ``velocity(x,y,z,t)`` for one real harmonic pair.

    The returned field is physical Cartesian velocity. This object is not a
    self-contained source reconstruction because the provider owns the missing
    corrected background/forcing/support realization.
    """

    def __init__(
        self,
        provider: AxisSafeSourceHarmonicProvider,
        *,
        ell: int,
        h: float,
    ) -> None:
        self._provider = _validated_provider(provider)
        self._scaling = source_physical_scaling(ell=ell, h=h)

    @property
    def provider_id(self) -> str:
        return self._provider.provider_id

    @property
    def provider_semantic_sha256(self) -> str:
        return self._provider.provider_semantic_sha256

    @property
    def scaling(self) -> SourcePhysicalScaling:
        return self._scaling

    @property
    def axis_zero_radius_chart(self) -> float:
        return float(self._provider.axis_zero_radius_chart)

    @property
    def axis_zero_radius_physical(self) -> float:
        radius = self._scaling.radial_scale * self.axis_zero_radius_chart
        if not math.isfinite(radius) or radius <= 0.0:
            raise RuntimeError("physical axis-zero radius is not finite and positive")
        return radius

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Evaluate physical Cartesian velocity with exact certified-axis zero.

        Inputs are broadcast to one common shape and the output has trailing
        dimension three. Provider evaluation is skipped entirely for points
        inside the declared zero core, avoiding arbitrary ``theta`` or clipped
        ``R`` in the cylindrical complete-curl formulas.
        """

        x_a = _finite_array(x, name="x")
        y_a = _finite_array(y, name="y")
        z_a = _finite_array(z, name="z")
        t_a = _finite_array(t, name="t")
        x_a, y_a, z_a, t_a = np.broadcast_arrays(x_a, y_a, z_a, t_a)

        radius_physical = np.hypot(x_a, y_a)
        theta = np.arctan2(y_a, x_a)
        R, Z, T = physical_to_source_chart_rzt(
            radius_physical,
            z_a,
            t_a,
            ell=self._scaling.ell,
            h=self._scaling.h,
        )
        active = R > self.axis_zero_radius_chart

        output = np.zeros(R.shape + (3,), dtype=float)
        if not np.any(active):
            return output

        flat_active = active.reshape(-1)
        R_active = R.reshape(-1)[flat_active]
        theta_active = theta.reshape(-1)[flat_active]
        Z_active = Z.reshape(-1)[flat_active]
        T_active = T.reshape(-1)[flat_active]

        harmonic = self._provider.evaluate_chart(
            R_active,
            theta_active,
            Z_active,
            T_active,
            self._scaling,
        )
        if not isinstance(harmonic, SourceLocalizedZeroDataHarmonic):
            raise ValueError("provider must return SourceLocalizedZeroDataHarmonic")

        physicalized = physicalize_source_localized_harmonic(
            harmonic,
            ell=self._scaling.ell,
            h=self._scaling.h,
        )
        velocity_active = np.asarray(
            physicalized.real_pair_velocity_cartesian_physical,
            dtype=float,
        )
        expected = (R_active.size, 3)
        if velocity_active.shape == (3,) and R_active.size == 1:
            velocity_active = velocity_active.reshape(1, 3)
        if velocity_active.shape != expected:
            raise ValueError(
                "provider harmonic batch shape does not match requested off-axis points: "
                f"expected {expected}, got {velocity_active.shape}"
            )
        if not np.all(np.isfinite(velocity_active)):
            raise ValueError("provider physicalized velocity must be finite")

        output.reshape(-1, 3)[flat_active] = velocity_active
        return output

    def report(self) -> dict[str, Any]:
        return {
            "task": TASK,
            "provider_id": self.provider_id,
            "provider_semantic_sha256": self.provider_semantic_sha256,
            "ell": self._scaling.ell,
            "h": self._scaling.h,
            "Q": self._scaling.Q,
            "axis_zero_radius_chart": self.axis_zero_radius_chart,
            "axis_zero_radius_physical": self.axis_zero_radius_physical,
            "contract": source_axis_safe_physical_velocity_contract(),
        }


def source_axis_safe_physical_velocity_contract() -> dict[str, Any]:
    """Return the source/autonomous/pending boundary for this evaluator."""

    parent = source_physicalized_localized_harmonic_contract()
    if parent["source_physical_Q_scaling_applied"] is not True:
        raise RuntimeError("required K2-OSC-102 physical scaling is unavailable")
    if parent["global_axis_safe_source_velocity_materialized"] is not False:
        raise RuntimeError("parent axis-safe truth boundary drifted")

    return {
        "schema": "kokuno-provider-driven-axis-safe-physical-velocity-v1",
        "task": TASK,
        "parent_a2_pr": PARENT_A2_PR,
        "parent_a2_head": PARENT_A2_HEAD,
        "parent_a2_source_blob": PARENT_A2_SOURCE_BLOB,
        "source_repository": SOURCE_REPOSITORY,
        "source_commit": SOURCE_COMMIT,
        "source_workbench_blob": SOURCE_WORKBENCH_BLOB,
        "source_chart_to_physical_scaling_consumed": True,
        "provider_driven_velocity_xyzt_materialized": True,
        "batch_evaluation_materialized": True,
        "certified_axis_zero_core_enforced": True,
        "axis_core_skips_singular_cylindrical_evaluation": True,
        "axis_zero_core_radius_is_repository_provider_contract": True,
        "source_exact_axis_support_radius_recovered": False,
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
