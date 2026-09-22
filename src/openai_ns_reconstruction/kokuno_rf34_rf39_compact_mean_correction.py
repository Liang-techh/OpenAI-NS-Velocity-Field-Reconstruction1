"""Source-specific RF34--RF39 compact mean-correction realization.

This module sits directly on the typed RF30/RF31 system from Agent-3 #1126.
The corrected 2026-09-09 Kokuno reader specifies, on the reserved mean patch,

    G_q = 0,
    V_q = a(eta) x^{-1-2 lambda},
    a(eta) = 2^{1/2+lambda} c_patch / (1+eta^2),

with lambda > 0.  It then chooses three multiplicatively separated compact
bumps for the azimuthal-profile correction ``u`` and two for the axial-profile
correction ``d``.  RF36 and RF38 give the five coefficients explicitly; RF34
maps the profile correction back to the fixed-Q chart.

The source leaves the concrete pointwise bump realization existential.  This
repository therefore freezes one deterministic autonomous C-infinity bump
family from the authenticated mean-patch interval *before* defect evaluation.
That choice is provenance, not recovered source data.

The only scientific entry point accepts a typed backend and candidate.  It
does not accept caller-entered defects, P/J targets, coefficients, bump
locations, correction arrays, held-out data, gains, forcing, pressure or
scientific thresholds.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Protocol

import numpy as np

from .kokuno_rf30_rf31_typed_mean_correction import (
    CandidateIdentity,
    RF30PreMeanState,
    RF30RF31Backend,
    RF30RF31ContractError,
    RF30RF31SystemReceipt,
    RF31CorrectionBasis,
    materialize_rf30_rf31_mean_correction_system,
)

TASK = "KOKUNO-A3-RF34-RF39-COMPACT-MEAN-CORRECTION-123"
SCHEMA = "kokuno-a3-rf34-rf39-compact-mean-correction-v1"

PARENT_AGENT3_PR = 1126
PARENT_AGENT3_HEAD = "edcf3a3a2e87f335c268460ca0d4b3b23f4a9847"
PARENT_SOURCE_BLOB = "fa31fdf12adceaa5b295db502cf0ccaaad778307"

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_READER_TEX_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_TEX_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_FORMULAS = ("RF29", "RF33", "RF34", "RF35", "RF36", "RF37", "RF38", "RF39", "RF43")

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5

# Engineering consistency gates for this source-specific operator, not NS gates.
BASIS_MATCH_ATOL = 2.0e-11
BASE_G_ATOL = 2.0e-11
BASE_V_REL_TOL = 2.0e-8
FIVE_ROW_REL_TOL = 2.0e-7
COEFFICIENT_SOLVE_REL_TOL = 2.0e-7
MOMENT_QUADRATURE_ORDER = 128


class RF34RF39ContractError(RF30RF31ContractError):
    """Raised when the source-specific compact correction contract fails."""


@dataclass(frozen=True)
class RF34RF39PatchContext:
    """Identity-bound source reserved-patch parameters.

    ``mean_patch_x_min/max`` describe I_m in profile radius x=r/sqrt(q).
    The concrete bump is *not* supplied by the caller/backend: this module
    deterministically freezes an autonomous one inside that interval.
    """

    source_candidate_id: str
    source_candidate_sha256: str
    source_chart_id: str
    source_chart_sha256: str
    mean_patch_x_min: float
    mean_patch_x_max: float
    scale_s_q_over_Q: float
    similarity_A: float
    eta: float
    lambda_value: float
    c_patch: float
    reserved_mean_patch_verified: bool
    context_frozen_before_defect_evaluation: bool
    heldout_samples_used_to_choose_context: bool = False
    source_exact_bump_claimed: bool = False


class RF34RF39Backend(RF30RF31Backend, Protocol):
    def rf34_rf39_patch_context(
        self, candidate: Any, state: RF30PreMeanState
    ) -> RF34RF39PatchContext: ...


@dataclass(frozen=True)
class RF34RF39CorrectionReceipt:
    candidate: CandidateIdentity
    parent_system_sha256: str
    source_chart_id: str
    source_chart_sha256: str
    patch_context_sha256: str
    autonomous_bump_family_sha256: str
    bump_log_centers: tuple[float, float, float]
    bump_log_halfwidth: float
    profile_moments: tuple[tuple[str, float], ...]
    profile_targets_Pq_Jthetaq_Jzq: tuple[float, float, float]
    coefficients_u0_u1_u2_s0_s1: tuple[float, float, float, float, float]
    direct_linear_solve_coefficients: tuple[float, float, float, float, float]
    delta_v_chart: tuple[float, ...]
    gamma_d_chart: tuple[float, ...]
    basis_match_max_abs: float
    reserved_base_G_max_abs: float
    reserved_base_V_max_relative_error: float
    five_row_closure_max_relative: float
    coefficient_direct_solve_max_relative: float
    correction_chart_l2: float
    correction_chart_max: float
    repository_candidate_correction_evidence: bool
    compact_support_preserved: bool
    two_zero_moments_preserved: bool
    source_rf34_rf39_correction_materialized: bool
    correction_applied_to_candidate: bool
    rf44_rf49_nonlinear_remainder_recomputed: bool
    cartesian_correction_velocity_materialized: bool
    complete_ns_defect: bool
    heldout_ns_residual_assessed: bool
    pde_validated: bool

    def to_dict(self) -> dict[str, object]:
        return {**asdict(self), "truth_boundary": truth_boundary()}


def _canonical_sha(payload: dict[str, object]) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _validate_context(
    identity: CandidateIdentity,
    state: RF30PreMeanState,
    context: RF34RF39PatchContext,
) -> None:
    if context.source_candidate_id != identity.candidate_id:
        raise RF34RF39ContractError("RF34/RF39 patch candidate id mismatch")
    if context.source_candidate_sha256 != identity.candidate_sha256:
        raise RF34RF39ContractError("RF34/RF39 patch candidate hash mismatch")
    if context.source_chart_id != state.source_chart_id:
        raise RF34RF39ContractError("RF34/RF39 patch chart id mismatch")
    if context.source_chart_sha256 != state.source_chart_sha256:
        raise RF34RF39ContractError("RF34/RF39 patch chart hash mismatch")
    values = (
        context.mean_patch_x_min,
        context.mean_patch_x_max,
        context.scale_s_q_over_Q,
        context.similarity_A,
        context.eta,
        context.lambda_value,
        context.c_patch,
    )
    if not all(math.isfinite(v) for v in values):
        raise RF34RF39ContractError("RF34/RF39 patch context is non-finite")
    if not (0.0 < context.mean_patch_x_min < context.mean_patch_x_max):
        raise RF34RF39ContractError("mean patch interval must be positive and ordered")
    if not (0.5 <= context.scale_s_q_over_Q <= 2.0):
        raise RF34RF39ContractError("source fixed-Q band requires 1/2 <= q/Q <= 2")
    if not (0.5 < context.similarity_A < 0.51):
        raise RF34RF39ContractError("source similarity A must equal 1/2+h with 0<h<1/100")
    if abs(context.eta) >= 1.0:
        raise RF34RF39ContractError("source similarity eta must satisfy |eta|<1")
    if context.lambda_value <= 0.0:
        raise RF34RF39ContractError("RF37/RF38 require fixed lambda > 0")
    if context.c_patch <= 0.0:
        raise RF34RF39ContractError("source reserved-patch c_patch must be positive")
    if not context.reserved_mean_patch_verified:
        raise RF34RF39ContractError("source reserved mean patch is unverified")
    if not context.context_frozen_before_defect_evaluation:
        raise RF34RF39ContractError("patch context must be frozen before defect evaluation")
    if context.heldout_samples_used_to_choose_context:
        raise RF34RF39ContractError("held-out samples cannot choose correction patch context")
    if context.source_exact_bump_claimed:
        raise RF34RF39ContractError(
            "the concrete repository bump is autonomous, not recovered source-exact data"
        )


def _autonomous_bump_geometry(
    context: RF34RF39PatchContext,
) -> tuple[float, float, float, tuple[float, float, float]]:
    """Freeze x0,d0,width from I_m without looking at the defect.

    Centers occupy 20%, 50%, 80% of log(I_m); halfwidth is 7.5% of
    log-span, leaving strict support gaps and strict containment.
    """

    lo = math.log(context.mean_patch_x_min)
    hi = math.log(context.mean_patch_x_max)
    span = hi - lo
    if not (span > 0.0):
        raise RF34RF39ContractError("mean-patch log span must be positive")
    log_x0 = lo + 0.20 * span
    d0 = 0.30 * span
    halfwidth = 0.075 * span
    centers = (log_x0, log_x0 + d0, log_x0 + 2.0 * d0)
    if not (lo < centers[0] - halfwidth and centers[2] + halfwidth < hi):
        raise RF34RF39ContractError("autonomous bump family does not fit inside I_m")
    if not (
        centers[0] + halfwidth < centers[1] - halfwidth
        and centers[1] + halfwidth < centers[2] - halfwidth
    ):
        raise RF34RF39ContractError("autonomous bump supports are not disjoint")
    return math.exp(log_x0), d0, halfwidth, centers


def _base_eta0(x: np.ndarray, x0: float, log_halfwidth: float) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    valid = x > 0.0
    y = np.zeros_like(x)
    y[valid] = np.log(x[valid] / x0) / log_halfwidth
    interior = valid & (np.abs(y) < 1.0)
    yy = y[interior]
    out[interior] = np.exp(-1.0 / (1.0 - yy * yy))
    return out


def _eta_j(
    x: np.ndarray, j: int, x0: float, d0: float, log_halfwidth: float
) -> np.ndarray:
    scale = math.exp(-j * d0)
    return scale * _base_eta0(scale * np.asarray(x, dtype=float), x0, log_halfwidth)


def _mu_p(p: float, x0: float, log_halfwidth: float) -> float:
    nodes, weights = np.polynomial.legendre.leggauss(MOMENT_QUADRATURE_ORDER)
    left = x0 * math.exp(-log_halfwidth)
    right = x0 * math.exp(log_halfwidth)
    x = 0.5 * (right - left) * nodes + 0.5 * (right + left)
    vals = np.power(x, p) * _base_eta0(x, x0, log_halfwidth)
    value = 0.5 * (right - left) * float(np.dot(weights, vals))
    if not (math.isfinite(value) and value > 0.0):
        raise RF34RF39ContractError(f"autonomous bump moment mu_{p} is invalid")
    return value


def _canonical_chart_basis(
    context: RF34RF39PatchContext,
    radius_R: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, dict[str, float], dict[str, object]]:
    x0, d0, halfwidth, centers = _autonomous_bump_geometry(context)
    s = context.scale_s_q_over_Q
    A = context.similarity_A
    x = np.asarray(radius_R, dtype=float) / math.sqrt(s)
    profile = np.vstack([_eta_j(x, j, x0, d0, halfwidth) for j in range(3)])
    prefactor = s ** (-A)
    dv = np.zeros((5, x.size), dtype=float)
    gd = np.zeros((5, x.size), dtype=float)
    dv[:3, :] = prefactor * profile
    gd[3, :] = prefactor * profile[0]
    gd[4, :] = prefactor * profile[1]

    lam = context.lambda_value
    required_p = (2.0, -2.0 - 2.0 * lam, -2.0 * lam, 1.0, 1.0 - 2.0 * lam)
    moments = {f"{p:.17g}": _mu_p(p, x0, halfwidth) for p in required_p}
    bump_payload: dict[str, object] = {
        "classification": "autonomous_repository_Cinfinity_bump_realization",
        "chosen_before_defect_evaluation": True,
        "selection_rule": "log-centers at 20/50/80 percent of authenticated I_m; halfwidth 7.5 percent",
        "x0": x0,
        "d0": d0,
        "log_halfwidth": halfwidth,
        "log_centers": centers,
        "moment_quadrature_order": MOMENT_QUADRATURE_ORDER,
        "source_exact_claim": False,
    }
    return dv, gd, moments, bump_payload


def _relative_max(actual: np.ndarray, expected: np.ndarray, floor: float = 1.0e-14) -> float:
    scale = max(float(np.max(np.abs(expected))), floor)
    return float(np.max(np.abs(actual - expected)) / scale)


def _source_coefficients(
    system: RF30RF31SystemReceipt,
    context: RF34RF39PatchContext,
    moments: dict[str, float],
    bump_payload: dict[str, object],
) -> tuple[np.ndarray, tuple[float, float, float]]:
    s = context.scale_s_q_over_Q
    A = context.similarity_A
    lam = context.lambda_value
    eta = context.eta
    c_patch = context.c_patch
    P = float(system.defect.P)
    Jt = float(system.defect.J_theta)
    Jz = float(system.defect.J_z)
    Pq = s ** (2.0 * A) * P
    Jtq = s ** (2.0 * A - 1.5) * Jt
    Jzq = s ** (2.0 * A - 1.0) * Jz

    a = 2.0 ** (0.5 + lam) * c_patch / (1.0 + eta * eta)
    d0 = float(bump_payload["d0"])
    p0 = 2.0
    p1 = -2.0 - 2.0 * lam
    p2 = -2.0 * lam
    t0, t1, t2 = (math.exp(d0 * p) for p in (p0, p1, p2))
    if min(abs(t1 - t0), abs(t2 - t0), abs(t2 - t1)) <= 1.0e-14:
        raise RF34RF39ContractError("RF37 angular determinant is numerically degenerate")

    mu_p1 = moments[f"{p1:.17g}"]
    mu_p2 = moments[f"{p2:.17g}"]
    b1 = -Pq / (2.0 * a * mu_p1)
    b2 = Jzq / (a * mu_p2)

    den1 = (t1 - t0) * (t1 - t2)
    den2 = (t2 - t0) * (t2 - t1)
    poly = np.zeros(3, dtype=float)
    poly += (b1 / den1) * np.array((t0 * t2, -(t0 + t2), 1.0))
    poly += (b2 / den2) * np.array((t0 * t1, -(t0 + t1), 1.0))

    tb = math.exp(d0)
    ta = math.exp(d0 * (1.0 - 2.0 * lam))
    mu_ax = moments[f"{1.0 - 2.0 * lam:.17g}"]
    denom_ax = a * mu_ax * (ta - tb)
    if abs(denom_ax) <= 1.0e-14:
        raise RF34RF39ContractError("RF38 axial determinant is numerically degenerate")
    s1 = -Jtq / denom_ax
    s0 = -tb * s1

    coeff = np.array((poly[0], poly[1], poly[2], s0, s1), dtype=float)
    if not np.all(np.isfinite(coeff)):
        raise RF34RF39ContractError("RF36/RF38 correction coefficients are non-finite")
    return coeff, (Pq, Jtq, Jzq)


def materialize_rf34_rf39_compact_mean_correction(
    backend: RF34RF39Backend,
    candidate: Any,
) -> RF34RF39CorrectionReceipt:
    """Materialize the source-specific compact RF34--RF39 chart correction."""

    context_method = getattr(backend, "rf34_rf39_patch_context", None)
    if not callable(context_method):
        raise RF34RF39ContractError("backend is missing rf34_rf39_patch_context")

    system = materialize_rf30_rf31_mean_correction_system(backend, candidate)
    identity = backend.candidate_identity(candidate)
    state = backend.rf30_pre_mean_state(candidate)
    parent_basis = backend.rf31_correction_basis(candidate, state)
    if not isinstance(identity, CandidateIdentity):
        raise RF34RF39ContractError("backend candidate identity type drifted")
    if not isinstance(state, RF30PreMeanState):
        raise RF34RF39ContractError("backend RF30 state type drifted")
    if not isinstance(parent_basis, RF31CorrectionBasis):
        raise RF34RF39ContractError("backend RF31 basis type drifted")

    context = context_method(candidate, state)
    if not isinstance(context, RF34RF39PatchContext):
        raise RF34RF39ContractError("backend must return RF34RF39PatchContext")
    _validate_context(identity, state, context)

    R = np.asarray(state.radius_R, dtype=float)
    dv_expected, gd_expected, moments, bump_payload = _canonical_chart_basis(context, R)
    dv_parent = np.asarray(parent_basis.delta_v_columns, dtype=float)
    gd_parent = np.asarray(parent_basis.gamma_d_columns, dtype=float)
    basis_match = max(
        float(np.max(np.abs(dv_parent - dv_expected))),
        float(np.max(np.abs(gd_parent - gd_expected))),
    )
    if basis_match > BASIS_MATCH_ATOL:
        raise RF34RF39ContractError(
            f"RF31 basis is not the frozen RF34--RF39 compact source realization: {basis_match}"
        )

    s = context.scale_s_q_over_Q
    A = context.similarity_A
    h = A - 0.5
    lam = context.lambda_value
    a = 2.0 ** (0.5 + lam) * context.c_patch / (1.0 + context.eta * context.eta)
    V_expected = a * s ** (lam - h) * np.power(R, -1.0 - 2.0 * lam)
    V = np.asarray(state.base_V, dtype=float)
    G = np.asarray(state.base_G, dtype=float)
    active = np.sum(np.abs(dv_expected) + np.abs(gd_expected), axis=0) > 1.0e-15
    if not np.any(active):
        raise RF34RF39ContractError("compact correction basis is zero on the RF30 grid")
    base_G_error = float(np.max(np.abs(G[active])))
    base_V_rel = _relative_max(V[active], V_expected[active])
    if base_G_error > BASE_G_ATOL:
        raise RF34RF39ContractError("reserved-patch G=0 law failed")
    if base_V_rel > BASE_V_REL_TOL:
        raise RF34RF39ContractError("reserved-patch RF43 V law failed")

    coeff, profile_targets = _source_coefficients(system, context, moments, bump_payload)
    B = np.asarray(system.matrix_B, dtype=float)
    target = np.asarray(system.defect.target_rows, dtype=float)
    if int(system.matrix_rank) != 5:
        raise RF34RF39ContractError("RF31 five-row matrix is rank deficient")
    direct = np.linalg.solve(B, target)
    direct_rel = _relative_max(coeff, direct)
    if direct_rel > COEFFICIENT_SOLVE_REL_TOL:
        raise RF34RF39ContractError(
            f"RF36/RF38 coefficients disagree with direct five-row solve: {direct_rel}"
        )
    row_result = B @ coeff
    row_rel = _relative_max(row_result, target)
    if row_rel > FIVE_ROW_REL_TOL:
        raise RF34RF39ContractError(
            f"RF34--RF39 correction does not close the five RF31 rows: {row_rel}"
        )

    delta_v = coeff @ dv_expected
    gamma_d = coeff @ gd_expected
    correction_norm = np.sqrt(delta_v * delta_v + gamma_d * gamma_d)
    context_sha = _canonical_sha(asdict(context))
    bump_sha = _canonical_sha(bump_payload)
    repository_evidence = bool(system.repository_candidate_defect_evidence)

    return RF34RF39CorrectionReceipt(
        candidate=identity,
        parent_system_sha256=system.system_sha256,
        source_chart_id=system.source_chart_id,
        source_chart_sha256=system.source_chart_sha256,
        patch_context_sha256=context_sha,
        autonomous_bump_family_sha256=bump_sha,
        bump_log_centers=tuple(float(v) for v in bump_payload["log_centers"]),
        bump_log_halfwidth=float(bump_payload["log_halfwidth"]),
        profile_moments=tuple(sorted((k, float(v)) for k, v in moments.items())),
        profile_targets_Pq_Jthetaq_Jzq=tuple(float(v) for v in profile_targets),
        coefficients_u0_u1_u2_s0_s1=tuple(float(v) for v in coeff),
        direct_linear_solve_coefficients=tuple(float(v) for v in direct),
        delta_v_chart=tuple(float(v) for v in delta_v),
        gamma_d_chart=tuple(float(v) for v in gamma_d),
        basis_match_max_abs=basis_match,
        reserved_base_G_max_abs=base_G_error,
        reserved_base_V_max_relative_error=base_V_rel,
        five_row_closure_max_relative=row_rel,
        coefficient_direct_solve_max_relative=direct_rel,
        correction_chart_l2=float(np.sqrt(np.mean(correction_norm * correction_norm))),
        correction_chart_max=float(np.max(correction_norm)),
        repository_candidate_correction_evidence=repository_evidence,
        compact_support_preserved=True,
        two_zero_moments_preserved=True,
        source_rf34_rf39_correction_materialized=True,
        correction_applied_to_candidate=False,
        rf44_rf49_nonlinear_remainder_recomputed=False,
        cartesian_correction_velocity_materialized=False,
        complete_ns_defect=False,
        heldout_ns_residual_assessed=False,
        pde_validated=False,
    )


def truth_boundary() -> dict[str, object]:
    return {
        "task": TASK,
        "schema": SCHEMA,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_date": SOURCE_READER_DATE,
        "source_reader_tex_path": SOURCE_READER_TEX_PATH,
        "source_reader_tex_blob": SOURCE_READER_TEX_BLOB,
        "source_formulas": SOURCE_FORMULAS,
        "rf34_rf39_compact_correction_realization_materialized": True,
        "concrete_bump_realization_classification": "autonomous_repository_choice",
        "source_exact_bump_recovered": False,
        "current_i4_source_chart_backend_materialized": False,
        "current_i4_rf30_defect_materialized": False,
        "current_i4_rf34_rf39_correction_materialized": False,
        "correction_applied_to_candidate": False,
        "rf44_rf49_nonlinear_remainder_recomputed": False,
        "cartesian_correction_velocity_materialized": False,
        "complete_ns_defect": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        "residual_defined_free_forcing_allowed": False,
        "finite_stage_small_residual_is_blowup_proof": False,
    }


@dataclass(frozen=True)
class _MechanicsCandidate:
    candidate_id: str = "rf34-rf39-mechanics"
    candidate_sha256: str = "4" * 64


class _MechanicsBackend:
    def __init__(self, *, mutate_basis: bool = False, lambda_value: float = 0.08) -> None:
        self.mutate_basis = bool(mutate_basis)
        self.lambda_value = float(lambda_value)
        self._identity = CandidateIdentity(
            candidate_id="rf34-rf39-mechanics",
            candidate_sha256="4" * 64,
            evidence_kind="mechanics-only",
        )
        self._R = np.exp(np.linspace(math.log(0.55), math.log(2.20), 24001))
        self._context = RF34RF39PatchContext(
            source_candidate_id=self._identity.candidate_id,
            source_candidate_sha256=self._identity.candidate_sha256,
            source_chart_id="mechanics-fixed-Q",
            source_chart_sha256="5" * 64,
            mean_patch_x_min=0.60,
            mean_patch_x_max=2.00,
            scale_s_q_over_Q=1.0,
            similarity_A=0.505,
            eta=0.23,
            lambda_value=self.lambda_value,
            c_patch=0.8,
            reserved_mean_patch_verified=True,
            context_frozen_before_defect_evaluation=True,
        )

    def candidate_identity(self, candidate: Any) -> CandidateIdentity:
        return self._identity

    def rf30_pre_mean_state(self, candidate: Any) -> RF30PreMeanState:
        R = self._R
        c = self._context
        h = c.similarity_A - 0.5
        a = 2.0 ** (0.5 + c.lambda_value) * c.c_patch / (1.0 + c.eta * c.eta)
        V = a * c.scale_s_q_over_Q ** (c.lambda_value - h) * R ** (-1.0 - 2.0 * c.lambda_value)
        bump = np.exp(-((R - 1.15) / 0.22) ** 2)
        Wtt = 0.04 * R * bump
        Wzt = -0.025 * bump
        Wzz = 0.03 * bump
        zeros = np.zeros_like(R)
        return RF30PreMeanState(
            source_candidate_id=self._identity.candidate_id,
            source_candidate_sha256=self._identity.candidate_sha256,
            source_chart_id=c.source_chart_id,
            source_chart_sha256=c.source_chart_sha256,
            radius_R=tuple(float(v) for v in R),
            base_V=tuple(float(v) for v in V),
            base_G=tuple(float(v) for v in zeros),
            mean_W_rr=tuple(float(v) for v in zeros),
            mean_W_zr=tuple(float(v) for v in zeros),
            mean_W_thetatheta=tuple(float(v) for v in Wtt),
            mean_W_ztheta=tuple(float(v) for v in Wzt),
            mean_W_zz=tuple(float(v) for v in Wzz),
            dR_mean_W_rr=tuple(float(v) for v in zeros),
            dZ_mean_W_zr=tuple(float(v) for v in zeros),
            actual_candidate_recomputed=False,
            oscillatory_covariance_recomputed=False,
            normalized_haar_mean_used=False,
            source_fixed_q_chart_used=True,
        )

    def rf31_correction_basis(
        self, candidate: Any, state: RF30PreMeanState
    ) -> RF31CorrectionBasis:
        R = np.asarray(state.radius_R, dtype=float)
        dv, gd, _, _ = _canonical_chart_basis(self._context, R)
        if self.mutate_basis:
            dv = dv.copy()
            dv[0, len(R) // 2] += 1.0e-4
        return RF31CorrectionBasis(
            source_candidate_id=self._identity.candidate_id,
            source_candidate_sha256=self._identity.candidate_sha256,
            source_chart_id=self._context.source_chart_id,
            source_chart_sha256=self._context.source_chart_sha256,
            radius_R=tuple(float(v) for v in R),
            delta_v_columns=tuple(tuple(float(v) for v in row) for row in dv),
            gamma_d_columns=tuple(tuple(float(v) for v in row) for row in gd),
            basis_selected_before_defect_evaluation=True,
            compact_support_verified=True,
            torus_independent_verified=True,
        )

    def rf34_rf39_patch_context(
        self, candidate: Any, state: RF30PreMeanState
    ) -> RF34RF39PatchContext:
        return self._context


def build_mechanics_report() -> dict[str, object]:
    candidate = _MechanicsCandidate()
    receipt = materialize_rf34_rf39_compact_mean_correction(_MechanicsBackend(), candidate)
    return {
        "schema": SCHEMA,
        "mechanics_only": True,
        "candidate_residual_evidence": False,
        "receipt": receipt.to_dict(),
    }


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = build_mechanics_report()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    _main()
