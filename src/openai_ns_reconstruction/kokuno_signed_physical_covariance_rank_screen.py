"""Agent-3 phase-mean covariance-rank screen for signed physical complete curls.

This module is deliberately downstream of
``KokunoSourceSignedCompleteCurlFamily``.  Agent 2 owns the signed auxiliary
rectangle amplitudes, slow derivatives, support localization, complete curl,
Fourier reality pair, and Q scaling.  Agent 3 only asks whether the resulting
*physical* ``sigma=+/-`` columns survive the declared phase/angle average as
independent mean-covariance directions.

Let ``W`` be the assembled cylindrical physical velocity and let ``S_sigma``
be the already-Q-scaled physical velocity contributed by rectangle sign
``sigma`` after summing all beta labels.  A dimensionless fractional multiplier
of one sign has tangent ``S_sigma``.  Repository product-rule algebra gives

    d_sigma C_theta = <S_sigma,r W_theta + W_r S_sigma,theta>,
    d_sigma C_z     = <S_sigma,r W_z     + W_r S_sigma,z>.

Using the *assembled* W retains self, cross-sign, and cross-beta interactions.
As a consistency identity, summing the two sign tangents must reproduce W, so
summing the two covariance-response columns must reproduce twice the assembled
covariance.  The screen checks this fail-closed before any SVD rank verdict.

A PASS here is only a supplied/source-compatible structural result.  Agent 2's
current physical family still contains caller-supplied unreleased source data,
so this module cannot promote source binding, correction readiness, finite-cycle
readiness, or a Navier--Stokes residual claim.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np


@dataclass(frozen=True)
class KokunoSignedPhysicalCovarianceRankScreen:
    """Measure phase-mean covariance rank of physical ``sigma=+/-`` columns."""

    rank_rtol: float = 1.0e-8
    consistency_atol: float = 2.0e-12

    def __post_init__(self) -> None:
        rank_rtol = float(self.rank_rtol)
        consistency_atol = float(self.consistency_atol)
        if not np.isfinite(rank_rtol) or not (0.0 < rank_rtol <= 1.0e-3):
            raise ValueError("rank_rtol must lie in (0,1e-3]")
        if not np.isfinite(consistency_atol) or not (0.0 < consistency_atol <= 1.0e-8):
            raise ValueError("consistency_atol must lie in (0,1e-8]")
        object.__setattr__(self, "rank_rtol", rank_rtol)
        object.__setattr__(self, "consistency_atol", consistency_atol)

    @staticmethod
    def _normalize_axes(axes: Sequence[int], sample_ndim: int) -> tuple[int, ...]:
        if sample_ndim < 1:
            raise ValueError("at least one sample axis is required for covariance averaging")
        normalized: list[int] = []
        for axis in axes:
            if isinstance(axis, bool) or not isinstance(axis, (int, np.integer)):
                raise ValueError("averaging_axes must contain integers")
            value = int(axis)
            if value < 0:
                value += sample_ndim
            if not (0 <= value < sample_ndim):
                raise ValueError("averaging axis is outside the physical sample axes")
            if value in normalized:
                raise ValueError("averaging_axes must not contain duplicates")
            normalized.append(value)
        if not normalized:
            raise ValueError("averaging_axes must be nonempty")
        return tuple(sorted(normalized))

    def evaluate(
        self,
        physical_family_result: dict[str, Any],
        *,
        averaging_axes: Sequence[int],
    ) -> dict[str, Any]:
        required = {
            "beta_labels",
            "sigma_labels",
            "velocity_physical_cylindrical_by_beta_sign",
            "velocity_physical_cylindrical_total",
            "reference_covariance_rank_two",
            "physical_complete_curl_covariance_rank_two_assessed",
            "genuinely_independent_second_covariance_column_ready",
        }
        missing = sorted(required - set(physical_family_result))
        if missing:
            raise ValueError(f"Agent-2 signed physical family is missing required keys: {missing}")
        if bool(physical_family_result["reference_covariance_rank_two"]) is not True:
            raise ValueError("source displayed signed reference inverse must be rank two")
        if tuple(physical_family_result["sigma_labels"]) != ("+", "-"):
            raise ValueError("signed physical family must preserve sigma ordering (+,-)")

        by_beta_sign = np.asarray(
            physical_family_result["velocity_physical_cylindrical_by_beta_sign"], dtype=float
        )
        total = np.asarray(physical_family_result["velocity_physical_cylindrical_total"], dtype=float)
        if by_beta_sign.ndim < 4 or by_beta_sign.shape[-2:] != (2, 3):
            raise ValueError(
                "signed physical columns must have shape sample_shape+(beta,2,3) with sample axes"
            )
        if total.shape != by_beta_sign.shape[:-3] + (3,):
            raise ValueError("physical total shape is inconsistent with sign-resolved columns")
        if not np.all(np.isfinite(by_beta_sign)) or not np.all(np.isfinite(total)):
            raise ValueError("physical signed columns and total must be finite")

        by_beta = np.sum(by_beta_sign, axis=-2)
        rebuilt_total = np.sum(by_beta, axis=-2)
        total_scale = max(1.0, float(np.max(np.abs(total), initial=0.0)))
        total_error = float(np.max(np.abs(rebuilt_total - total), initial=0.0))
        if total_error > self.consistency_atol * total_scale:
            raise RuntimeError("sign/beta columns do not reproduce the assembled physical total")

        # Sum over beta, retaining the source rectangle sign as the two tangent axes.
        sign_tangents = np.sum(by_beta_sign, axis=-3)
        tangent_sum = np.sum(sign_tangents, axis=-2)
        tangent_error = float(np.max(np.abs(tangent_sum - total), initial=0.0))
        if tangent_error > self.consistency_atol * total_scale:
            raise RuntimeError("fractional sign tangents do not sum to the assembled physical total")

        sample_ndim = total.ndim - 1
        axes = self._normalize_axes(averaging_axes, sample_ndim)

        Wr = total[..., 0][..., None]
        Wtheta = total[..., 1][..., None]
        Wz = total[..., 2][..., None]
        Sr = sign_tangents[..., 0]
        Stheta = sign_tangents[..., 1]
        Sz = sign_tangents[..., 2]

        response_theta = np.mean(Sr * Wtheta + Wr * Stheta, axis=axes)
        response_z = np.mean(Sr * Wz + Wr * Sz, axis=axes)
        jacobian = np.stack((response_theta, response_z), axis=-1)
        if jacobian.shape[-2:] != (2, 2):
            raise RuntimeError("covariance Jacobian must end in (sigma,channel)=(2,2)")

        covariance_theta = np.mean(total[..., 0] * total[..., 1], axis=axes)
        covariance_z = np.mean(total[..., 0] * total[..., 2], axis=axes)
        covariance = np.stack((covariance_theta, covariance_z), axis=-1)
        homogeneous_lhs = np.sum(jacobian, axis=-2)
        homogeneous_rhs = 2.0 * covariance
        homogeneous_scale = np.maximum(1.0, np.abs(homogeneous_rhs))
        homogeneous_relative_error = np.abs(homogeneous_lhs - homogeneous_rhs) / homogeneous_scale
        maximum_homogeneous_relative_error = float(
            np.max(homogeneous_relative_error, initial=0.0)
        )
        if maximum_homogeneous_relative_error > self.consistency_atol:
            raise RuntimeError("signed covariance response dropped cross interactions")

        singular_values = np.linalg.svd(jacobian, compute_uv=False)
        smax = singular_values[..., 0]
        smin = singular_values[..., 1]
        rank_two = (smax > 0.0) & (smin > self.rank_rtol * smax)
        ratio = np.divide(smin, smax, out=np.zeros_like(smin), where=smax > 0.0)

        plus = jacobian[..., 0, :]
        minus = jacobian[..., 1, :]
        plus_norm = np.linalg.norm(plus, axis=-1)
        minus_norm = np.linalg.norm(minus, axis=-1)
        dot = np.sum(plus * minus, axis=-1)
        denom = plus_norm * minus_norm
        cosine = np.divide(dot, denom, out=np.zeros_like(dot), where=denom > 0.0)
        cosine = np.clip(cosine, -1.0, 1.0)
        novelty = np.sqrt(np.maximum(0.0, 1.0 - cosine * cosine))

        flat_rank = np.asarray(rank_two, dtype=bool).reshape(-1)
        total_cells = int(flat_rank.size)
        rank_two_cells = int(np.count_nonzero(flat_rank))
        local_pass = bool(total_cells > 0 and rank_two_cells == total_cells)

        return {
            "beta_labels": tuple(physical_family_result["beta_labels"]),
            "sigma_labels": ("+", "-"),
            "parameter_names": ("delta_sigma_plus", "delta_sigma_minus"),
            "parameter_units": (
                "dimensionless fractional multiplier of Q-scaled physical sigma+ velocity",
                "dimensionless fractional multiplier of Q-scaled physical sigma- velocity",
            ),
            "averaging_axes": axes,
            "covariance": covariance,
            "covariance_response_jacobian": jacobian,
            "singular_values": singular_values,
            "rank_two_mask": rank_two,
            "rank_two_cells": rank_two_cells,
            "total_cells": total_cells,
            "rank_two_fraction": float(rank_two_cells / total_cells) if total_cells else 0.0,
            "minimum_smallest_singular_value": float(np.min(smin)) if total_cells else 0.0,
            "minimum_singular_value_ratio": float(np.min(ratio)) if total_cells else 0.0,
            "signed_response_novelty": novelty,
            "minimum_signed_response_novelty": float(np.min(novelty)) if total_cells else 0.0,
            "maximum_total_reconstruction_abs_error": total_error,
            "maximum_sign_tangent_sum_abs_error": tangent_error,
            "maximum_quadratic_homogeneity_relative_error": maximum_homogeneous_relative_error,
            "cross_sign_and_cross_beta_terms_retained": True,
            "local_supplied_signed_physical_covariance_rank_two": local_pass,
            "reference_covariance_rank_two": True,
            "physical_complete_curl_covariance_rank_two_assessed": True,
            "actual_positive_order_background_bound": False,
            "actual_source_h_sigma_pulse_integrals_bound": False,
            "actual_signed_auxiliary_rectangles_bound": False,
            "actual_auxiliary_torus_mode_family_bound": False,
            "source_actual_partition_labels_instantiated": False,
            "genuinely_independent_second_covariance_column_ready": False,
            "real_candidate_defect_consumed": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_rerun_allowed": False,
            "heldout_ns_residual_assessed": False,
            "residual_reduction_claimed": False,
            "pde_validated": False,
            "scope": (
                "phase-mean covariance rank of supplied sign-resolved physical complete curls; "
                "no real-defect bounded inverse, materialized correction, or NS residual"
            ),
        }
