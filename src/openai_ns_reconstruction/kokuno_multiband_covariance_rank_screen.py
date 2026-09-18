"""Agent-3 mean-covariance rank screen for Agent-2 supplied multi-band families.

This module is downstream of ``KokunoSourceMultiBandRealPairFamily``. It does
not reconstruct pulses, cutoffs, complete curls, or source band schedules.
Instead it consumes the already-Q-scaled physical by-beta velocity family and
the bounded two-coordinate contract introduced for Agent 3.

For base cylindrical velocity ``W`` and a coefficient tangent ``S_j`` the
phase-mean covariance response is repository product-rule algebra

    d_j C_theta = <S_{j,r} W_theta + W_r S_{j,theta}>,
    d_j C_z     = <S_{j,r} W_z     + W_r S_{j,z}>.

The two response columns are then tested for local rank two after averaging
only over caller-declared sample axes. A local PASS is intentionally weaker
than source readiness: caller-supplied/source-compatible mode data can prove
that the public interface has independent covariance capacity, but cannot be
renamed recovered Kokuno data. The current source family is not bound, so this
module cannot promote genuine second-column readiness.

This is a structural preflight, not a Navier--Stokes residual and not a finite
correction cycle.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np

from .kokuno_bounded_family_coefficients import KokunoBoundedFamilyCoefficientCoordinates


@dataclass(frozen=True)
class KokunoMultiBandCovarianceRankScreen:
    """Measure rank of common/band mean-covariance tangent responses."""

    rank_rtol: float = 1.0e-8
    consistency_atol: float = 2.0e-12

    def __post_init__(self) -> None:
        rtol = float(self.rank_rtol)
        atol = float(self.consistency_atol)
        if not np.isfinite(rtol) or not (0.0 < rtol <= 1.0e-3):
            raise ValueError("rank_rtol must lie in (0,1e-3]")
        if not np.isfinite(atol) or not (0.0 < atol <= 1.0e-8):
            raise ValueError("consistency_atol must lie in (0,1e-8]")
        object.__setattr__(self, "rank_rtol", rtol)
        object.__setattr__(self, "consistency_atol", atol)

    @staticmethod
    def _normalize_axes(axes: Sequence[int], sample_ndim: int) -> tuple[int, ...]:
        if sample_ndim < 1:
            raise ValueError("at least one sample axis is required for phase averaging")
        out: list[int] = []
        for axis in axes:
            if isinstance(axis, bool) or not isinstance(axis, (int, np.integer)):
                raise ValueError("averaging_axes must contain integers")
            a = int(axis)
            if a < 0:
                a += sample_ndim
            if not (0 <= a < sample_ndim):
                raise ValueError("averaging axis is outside the sample axes")
            if a in out:
                raise ValueError("averaging_axes must not contain duplicates")
            out.append(a)
        if not out:
            raise ValueError("averaging_axes must be nonempty")
        return tuple(sorted(out))

    def evaluate(
        self,
        physical_family_result: dict[str, Any],
        *,
        averaging_axes: Sequence[int],
        coefficient_max_l1_update: float = 0.125,
    ) -> dict[str, Any]:
        """Consume Agent-2 physical family output and measure covariance rank.

        ``averaging_axes`` index only the sample axes of the cylindrical total
        (all axes before the final component axis). Typical usage is a phase
        or phase/angle quadrature axis, leaving radial/time cells explicit.
        """
        required = {
            "beta_labels",
            "active_ell_bands",
            "velocity_physical_cylindrical_by_beta",
            "velocity_physical_cylindrical_total",
            "deterministic_permutation_invariant_aggregation",
        }
        missing = sorted(required - set(physical_family_result))
        if missing:
            raise ValueError(f"Agent-2 physical family is missing required keys: {missing}")
        active_bands = tuple(int(v) for v in physical_family_result["active_ell_bands"])
        if len(set(active_bands)) < 2:
            raise ValueError("covariance rank screen requires at least two distinct dyadic bands")
        if physical_family_result["deterministic_permutation_invariant_aggregation"] is not True:
            raise ValueError("Agent-2 deterministic aggregation preflight must be true")

        coordinates = KokunoBoundedFamilyCoefficientCoordinates(
            max_l1_update=coefficient_max_l1_update,
            consistency_atol=self.consistency_atol,
        ).consume_physical_family(physical_family_result)
        names = tuple(coordinates["parameter_names"])
        if names != ("delta_common", "delta_band"):
            raise RuntimeError("unexpected bounded-coordinate ordering")
        if coordinates["band_contrast_active"] is not True:
            raise RuntimeError("multi-band family unexpectedly has inactive band contrast")

        W = np.asarray(coordinates["velocity_physical_cylindrical_total_base"], dtype=float)
        S = np.asarray(coordinates["velocity_parameter_tangent_cylindrical"], dtype=float)
        if W.ndim < 2 or W.shape[-1] != 3:
            raise ValueError("base velocity must have shape sample_shape+(3,)")
        if S.shape != W.shape[:-1] + (2, 3):
            raise ValueError("bounded tangent must have shape sample_shape+(2,3)")
        if not np.all(np.isfinite(W)) or not np.all(np.isfinite(S)):
            raise ValueError("base velocity and tangents must be finite")

        sample_ndim = W.ndim - 1
        axes = self._normalize_axes(averaging_axes, sample_ndim)
        Wr = W[..., 0][..., None]
        Wtheta = W[..., 1][..., None]
        Wz = W[..., 2][..., None]
        Sr = S[..., 0]
        Stheta = S[..., 1]
        Sz = S[..., 2]

        response_theta = np.mean(Sr * Wtheta + Wr * Stheta, axis=axes)
        response_z = np.mean(Sr * Wz + Wr * Sz, axis=axes)
        jacobian = np.stack((response_theta, response_z), axis=-1)
        if jacobian.shape[-2:] != (2, 2):
            raise RuntimeError("internal covariance Jacobian must end in (parameter,channel)=(2,2)")

        singular_values = np.linalg.svd(jacobian, compute_uv=False)
        smax = singular_values[..., 0]
        smin = singular_values[..., 1]
        rank_two = (smax > 0.0) & (smin > self.rank_rtol * smax)

        common = jacobian[..., 0, :]
        band = jacobian[..., 1, :]
        common_norm = np.linalg.norm(common, axis=-1)
        band_norm = np.linalg.norm(band, axis=-1)
        dot = np.sum(common * band, axis=-1)
        denom = common_norm * band_norm
        cosine = np.divide(dot, denom, out=np.zeros_like(dot), where=denom > 0.0)
        cosine = np.clip(cosine, -1.0, 1.0)
        novelty = np.sqrt(np.maximum(0.0, 1.0 - cosine * cosine))

        flat_rank = np.asarray(rank_two, dtype=bool).reshape(-1)
        total_cells = int(flat_rank.size)
        rank_two_cells = int(np.count_nonzero(flat_rank))
        local_pass = bool(total_cells > 0 and rank_two_cells == total_cells)
        ratio = np.divide(smin, smax, out=np.zeros_like(smin), where=smax > 0.0)

        return {
            "active_ell_bands": active_bands,
            "parameter_names": names,
            "parameter_units": tuple(coordinates["parameter_units"]),
            "repository_coefficient_unit_mapping_available": bool(
                coordinates["repository_coefficient_unit_mapping_available"]
            ),
            "source_coefficient_unit_mapping_available": bool(
                coordinates["source_coefficient_unit_mapping_available"]
            ),
            "averaging_axes": axes,
            "covariance_response_jacobian": jacobian,
            "singular_values": singular_values,
            "rank_two_mask": rank_two,
            "rank_two_cells": rank_two_cells,
            "total_cells": total_cells,
            "rank_two_fraction": float(rank_two_cells / total_cells) if total_cells else 0.0,
            "minimum_smallest_singular_value": float(np.min(smin)) if total_cells else 0.0,
            "minimum_singular_value_ratio": float(np.min(ratio)) if total_cells else 0.0,
            "band_response_novelty": novelty,
            "minimum_band_response_novelty": float(np.min(novelty)) if total_cells else 0.0,
            "local_supplied_family_covariance_rank_two": local_pass,
            "actual_source_mode_family_bound": False,
            "genuinely_independent_second_covariance_column_ready": False,
            "finite_correction_cycle_rerun_allowed": False,
            "full_ns_residual_assessed": False,
            "residual_reduction_claimed": False,
            "pde_validated": False,
            "scope": "mean covariance rank of bounded common/band tangents; no stress inverse or NS residual",
        }
