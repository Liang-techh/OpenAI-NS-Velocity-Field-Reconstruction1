"""Typed Kokuno RF30/RF31 pre-mean defect and five-row correction system.

Pinned structural provenance is the corrected public 2026-09-09 Kokuno reader.
Before the mean update, RF30 sets v=gamma=beta=0 and retains the oscillatory
covariance W. In the source fixed-Q cylindrical chart,

    g_r = -(D_R + 1/R) W^{rr} - D_Z W^{zr} + W^{theta theta}/R,

    P       = int g_r dR,
    J_theta = int R^2 W^{z theta} dR,
    J_z     = int R W^{zz} dR - 1/2 int R^2 g_r dR.

RF31 then prescribes five linear correction rows for compact torus-independent
increments (Delta v, gamma_d):

    int R^2 Delta v dR                                  = 0,
    int R gamma_d dR                                    = 0,
    int (2V/R) Delta v dR                               = -P,
    int R^2 (G Delta v + V gamma_d) dR                  = -J_theta,
    int (2 R G gamma_d - R V Delta v) dR                = -J_z.

This module makes those formulas executable behind a typed backend.  A caller
cannot enter P/J targets, residuals, forcing, pressure, correction coefficients,
or thresholds.  Repository-candidate evidence is admitted only when the source
state and five-column basis are bound to the same candidate/chart identity and
the state declares actual covariance recomputation rather than a surrogate.

The module intentionally does not solve/apply the correction, recompute RF44--
RF49 nonlinear remainders, construct a Cartesian correction velocity, or claim
a complete Navier--Stokes defect / residual improvement.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from typing import Any, Protocol

import numpy as np


TASK = "KOKUNO-A3-RF30-RF31-TYPED-MEAN-CORRECTION-SYSTEM-122"
SCHEMA = "kokuno-a3-rf30-rf31-typed-mean-correction-system-v1"

PARENT_AGENT3_PR = 1118
PARENT_AGENT3_HEAD = "922f7aa10460ded44af212d313eceff33a2ac647"
PARENT_AGENT3_SOURCE_BLOB = "e28eac6f7fb7f3191a024b052bcf677cff8d74d2"

SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_READER_TEX_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_TEX_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_FORMULAS = ("RF30", "RF31")

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5
CORRECTION_DIMENSION = 5


class RF30RF31ContractError(RuntimeError):
    """Raised when the typed source-state correction contract is violated."""


@dataclass(frozen=True)
class CandidateIdentity:
    candidate_id: str
    candidate_sha256: str
    evidence_kind: str

    def __post_init__(self) -> None:
        if not self.candidate_id or not self.candidate_sha256:
            raise ValueError("candidate identity must be nonempty")
        if self.evidence_kind not in {"mechanics-only", "repository-candidate"}:
            raise ValueError("unsupported evidence_kind")


@dataclass(frozen=True)
class RF30PreMeanState:
    """Azimuthally averaged source-chart state needed by RF30.

    Arrays are sampled on one strictly positive, increasing source-chart radial
    grid at fixed slow variables. Derivatives are derivatives of the averaged
    covariance components in that same chart.
    """

    source_candidate_id: str
    source_candidate_sha256: str
    source_chart_id: str
    source_chart_sha256: str
    radius_R: tuple[float, ...]
    base_V: tuple[float, ...]
    base_G: tuple[float, ...]
    mean_W_rr: tuple[float, ...]
    mean_W_zr: tuple[float, ...]
    mean_W_thetatheta: tuple[float, ...]
    mean_W_ztheta: tuple[float, ...]
    mean_W_zz: tuple[float, ...]
    dR_mean_W_rr: tuple[float, ...]
    dZ_mean_W_zr: tuple[float, ...]
    actual_candidate_recomputed: bool
    oscillatory_covariance_recomputed: bool
    normalized_haar_mean_used: bool
    source_fixed_q_chart_used: bool
    surrogate_defect_used: bool = False
    residual_as_forcing_shortcut_used: bool = False
    heldout_samples_used_to_construct_state: bool = False


@dataclass(frozen=True)
class RF31CorrectionBasis:
    """Five compact torus-independent RF31 basis columns on the same R grid."""

    source_candidate_id: str
    source_candidate_sha256: str
    source_chart_id: str
    source_chart_sha256: str
    radius_R: tuple[float, ...]
    delta_v_columns: tuple[tuple[float, ...], ...]
    gamma_d_columns: tuple[tuple[float, ...], ...]
    basis_selected_before_defect_evaluation: bool
    compact_support_verified: bool
    torus_independent_verified: bool
    surrogate_target_used_to_choose_basis: bool = False
    heldout_samples_used_to_choose_basis: bool = False


class RF30RF31Backend(Protocol):
    def candidate_identity(self, candidate: Any) -> CandidateIdentity: ...

    def rf30_pre_mean_state(self, candidate: Any) -> RF30PreMeanState: ...

    def rf31_correction_basis(
        self, candidate: Any, state: RF30PreMeanState
    ) -> RF31CorrectionBasis: ...


@dataclass(frozen=True)
class RF30DefectTuple:
    P: float
    J_theta: float
    J_z: float
    target_rows: tuple[float, ...]
    radial_source_g_r: tuple[float, ...]
    state_sha256: str


@dataclass(frozen=True)
class RF30RF31SystemReceipt:
    candidate: CandidateIdentity
    source_chart_id: str
    source_chart_sha256: str
    defect: RF30DefectTuple
    matrix_B: tuple[tuple[float, ...], ...]
    matrix_rank: int
    matrix_condition_2: float
    system_sha256: str
    repository_candidate_defect_evidence: bool
    linear_five_row_system_materialized: bool
    correction_coefficients_solved: bool
    correction_applied: bool
    nonlinear_remainder_recomputed: bool
    complete_ns_defect: bool
    heldout_ns_residual_assessed: bool
    pde_validated: bool

    def to_dict(self) -> dict[str, object]:
        return {**asdict(self), "truth_boundary": truth_boundary()}


def _require_backend(backend: Any) -> None:
    missing = [
        name
        for name in (
            "candidate_identity",
            "rf30_pre_mean_state",
            "rf31_correction_basis",
        )
        if not callable(getattr(backend, name, None))
    ]
    if missing:
        raise RF30RF31ContractError(
            "backend is missing typed RF30/RF31 methods: " + ", ".join(missing)
        )


def _as_vector(values: tuple[float, ...], n: int, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != (n,):
        raise RF30RF31ContractError(f"{name} must have shape ({n},)")
    if not np.all(np.isfinite(arr)):
        raise RF30RF31ContractError(f"{name} contains non-finite values")
    return arr


def _validate_radius(values: tuple[float, ...]) -> np.ndarray:
    R = np.asarray(values, dtype=float)
    if R.ndim != 1 or R.size < 3:
        raise RF30RF31ContractError("RF30 radius grid must contain at least 3 nodes")
    if not np.all(np.isfinite(R)):
        raise RF30RF31ContractError("RF30 radius grid contains non-finite values")
    if np.any(R <= 0.0) or np.any(np.diff(R) <= 0.0):
        raise RF30RF31ContractError(
            "RF30 source-chart radius must be strictly positive and increasing"
        )
    return R


def _trapz(values: np.ndarray, R: np.ndarray) -> float:
    return float(np.sum(0.5 * (values[:-1] + values[1:]) * np.diff(R)))


def _canonical_sha(payload: dict[str, object]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _validate_identity(
    identity: CandidateIdentity,
    state: RF30PreMeanState,
    basis: RF31CorrectionBasis,
) -> None:
    for label, candidate_id, candidate_sha in (
        ("RF30 state", state.source_candidate_id, state.source_candidate_sha256),
        ("RF31 basis", basis.source_candidate_id, basis.source_candidate_sha256),
    ):
        if candidate_id != identity.candidate_id:
            raise RF30RF31ContractError(f"{label} candidate id mismatch")
        if candidate_sha != identity.candidate_sha256:
            raise RF30RF31ContractError(f"{label} candidate hash mismatch")

    if not state.source_chart_id or not state.source_chart_sha256:
        raise RF30RF31ContractError("RF30 source chart identity must be nonempty")
    if basis.source_chart_id != state.source_chart_id:
        raise RF30RF31ContractError("RF31 basis source chart id mismatch")
    if basis.source_chart_sha256 != state.source_chart_sha256:
        raise RF30RF31ContractError("RF31 basis source chart hash mismatch")
    if tuple(basis.radius_R) != tuple(state.radius_R):
        raise RF30RF31ContractError("RF31 basis radial grid differs from RF30 state")


def _validate_provenance(
    identity: CandidateIdentity,
    state: RF30PreMeanState,
    basis: RF31CorrectionBasis,
) -> bool:
    if state.surrogate_defect_used:
        raise RF30RF31ContractError("surrogate RF30 defect is forbidden")
    if state.residual_as_forcing_shortcut_used:
        raise RF30RF31ContractError("residual-as-forcing shortcut is forbidden")
    if state.heldout_samples_used_to_construct_state:
        raise RF30RF31ContractError("held-out samples cannot construct RF30 fit state")
    if basis.surrogate_target_used_to_choose_basis:
        raise RF30RF31ContractError("surrogate defect cannot choose the RF31 basis")
    if basis.heldout_samples_used_to_choose_basis:
        raise RF30RF31ContractError("held-out samples cannot choose RF31 basis")
    if not basis.basis_selected_before_defect_evaluation:
        raise RF30RF31ContractError("RF31 basis must be frozen before defect evaluation")
    if not basis.compact_support_verified:
        raise RF30RF31ContractError("RF31 basis compact support is unverified")
    if not basis.torus_independent_verified:
        raise RF30RF31ContractError("RF31 basis must be torus independent")

    repository_evidence = identity.evidence_kind == "repository-candidate"
    if repository_evidence:
        if not state.actual_candidate_recomputed:
            raise RF30RF31ContractError(
                "repository-candidate RF30 state was not recomputed from candidate"
            )
        if not state.oscillatory_covariance_recomputed:
            raise RF30RF31ContractError(
                "repository-candidate oscillatory covariance was not recomputed"
            )
        if not state.normalized_haar_mean_used:
            raise RF30RF31ContractError(
                "repository-candidate RF30 state lacks normalized Haar mean"
            )
        if not state.source_fixed_q_chart_used:
            raise RF30RF31ContractError(
                "repository-candidate RF30 state is not in source fixed-Q chart"
            )
    return repository_evidence


def _rf30_defect(state: RF30PreMeanState) -> tuple[RF30DefectTuple, np.ndarray, np.ndarray, np.ndarray]:
    R = _validate_radius(state.radius_R)
    n = int(R.size)
    V = _as_vector(state.base_V, n, "base_V")
    G = _as_vector(state.base_G, n, "base_G")
    Wrr = _as_vector(state.mean_W_rr, n, "mean_W_rr")
    Wzr = _as_vector(state.mean_W_zr, n, "mean_W_zr")
    Wtt = _as_vector(state.mean_W_thetatheta, n, "mean_W_thetatheta")
    Wzt = _as_vector(state.mean_W_ztheta, n, "mean_W_ztheta")
    Wzz = _as_vector(state.mean_W_zz, n, "mean_W_zz")
    dR_Wrr = _as_vector(state.dR_mean_W_rr, n, "dR_mean_W_rr")
    dZ_Wzr = _as_vector(state.dZ_mean_W_zr, n, "dZ_mean_W_zr")

    # RF30 initial mean state: v = gamma = beta = 0.
    g_r = -(dR_Wrr + Wrr / R) - dZ_Wzr + Wtt / R
    P = _trapz(g_r, R)
    J_theta = _trapz((R**2) * Wzt, R)
    J_z = _trapz(R * Wzz, R) - 0.5 * _trapz((R**2) * g_r, R)

    scalars = (P, J_theta, J_z)
    if not all(math.isfinite(value) for value in scalars):
        raise RF30RF31ContractError("RF30 integrated defects are non-finite")

    state_sha = _canonical_sha(
        {
            "schema": SCHEMA,
            "source_formulas": SOURCE_FORMULAS,
            "source_reader_head": SOURCE_READER_HEAD,
            "source_reader_tex_blob": SOURCE_READER_TEX_BLOB,
            "state": asdict(state),
        }
    )
    defect = RF30DefectTuple(
        P=P,
        J_theta=J_theta,
        J_z=J_z,
        target_rows=(0.0, 0.0, -P, -J_theta, -J_z),
        radial_source_g_r=tuple(float(value) for value in g_r),
        state_sha256=state_sha,
    )
    return defect, R, V, G


def _rf31_matrix(
    basis: RF31CorrectionBasis,
    R: np.ndarray,
    V: np.ndarray,
    G: np.ndarray,
) -> np.ndarray:
    dv = np.asarray(basis.delta_v_columns, dtype=float)
    gd = np.asarray(basis.gamma_d_columns, dtype=float)
    expected = (CORRECTION_DIMENSION, int(R.size))
    if dv.shape != expected or gd.shape != expected:
        raise RF30RF31ContractError(
            f"RF31 basis columns must both have shape {expected}"
        )
    if not np.all(np.isfinite(dv)) or not np.all(np.isfinite(gd)):
        raise RF30RF31ContractError("RF31 basis contains non-finite values")

    B = np.empty((CORRECTION_DIMENSION, CORRECTION_DIMENSION), dtype=float)
    for j in range(CORRECTION_DIMENSION):
        delta_v = dv[j]
        gamma_d = gd[j]
        B[0, j] = _trapz((R**2) * delta_v, R)
        B[1, j] = _trapz(R * gamma_d, R)
        B[2, j] = _trapz((2.0 * V / R) * delta_v, R)
        B[3, j] = _trapz((R**2) * (G * delta_v + V * gamma_d), R)
        B[4, j] = _trapz(2.0 * R * G * gamma_d - R * V * delta_v, R)

    if not np.all(np.isfinite(B)):
        raise RF30RF31ContractError("RF31 five-row matrix is non-finite")
    return B


def materialize_rf30_rf31_mean_correction_system(
    backend: RF30RF31Backend,
    candidate: Any,
) -> RF30RF31SystemReceipt:
    """Materialize exact RF30 targets and RF31 five-row matrix from typed data."""

    _require_backend(backend)
    identity = backend.candidate_identity(candidate)
    if not isinstance(identity, CandidateIdentity):
        raise RF30RF31ContractError("backend must return CandidateIdentity")

    state = backend.rf30_pre_mean_state(candidate)
    if not isinstance(state, RF30PreMeanState):
        raise RF30RF31ContractError("backend must return RF30PreMeanState")

    basis = backend.rf31_correction_basis(candidate, state)
    if not isinstance(basis, RF31CorrectionBasis):
        raise RF30RF31ContractError("backend must return RF31CorrectionBasis")

    _validate_identity(identity, state, basis)
    repository_evidence = _validate_provenance(identity, state, basis)
    defect, R, V, G = _rf30_defect(state)
    B = _rf31_matrix(basis, R, V, G)

    rank = int(np.linalg.matrix_rank(B))
    condition = float(np.linalg.cond(B, p=2))
    if not math.isfinite(condition):
        condition = math.inf

    matrix_tuple = tuple(
        tuple(float(value) for value in row)
        for row in B
    )
    system_sha = _canonical_sha(
        {
            "schema": SCHEMA,
            "candidate": asdict(identity),
            "source_chart_id": state.source_chart_id,
            "source_chart_sha256": state.source_chart_sha256,
            "defect": asdict(defect),
            "matrix_B": matrix_tuple,
            "basis": asdict(basis),
            "source_reader_head": SOURCE_READER_HEAD,
            "source_reader_tex_blob": SOURCE_READER_TEX_BLOB,
        }
    )

    return RF30RF31SystemReceipt(
        candidate=identity,
        source_chart_id=state.source_chart_id,
        source_chart_sha256=state.source_chart_sha256,
        defect=defect,
        matrix_B=matrix_tuple,
        matrix_rank=rank,
        matrix_condition_2=condition,
        system_sha256=system_sha,
        repository_candidate_defect_evidence=repository_evidence,
        linear_five_row_system_materialized=True,
        correction_coefficients_solved=False,
        correction_applied=False,
        nonlinear_remainder_recomputed=False,
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
        "rf30_rf31_formula_executor_materialized": True,
        "current_i4_source_chart_backend_materialized": False,
        "current_i4_rf30_defect_materialized": False,
        "current_i4_rf31_five_row_system_materialized": False,
        "correction_coefficients_solved": False,
        "correction_applied": False,
        "rf44_rf49_nonlinear_remainder_recomputed": False,
        "cartesian_correction_velocity_materialized": False,
        "complete_ns_defect": False,
        "finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        "residual_defined_free_forcing_allowed": False,
        "finite_stage_small_residual_is_blowup_proof": False,
    }
