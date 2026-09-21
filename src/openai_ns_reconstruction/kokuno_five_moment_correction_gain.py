"""Typed Kokuno five-moment correction contraction/gain diagnostic.

The corrected 2026-09-09 Kokuno reader records the exact quadratic moment
correction equation

    B c + Q(c,c) = d

and the fixed-point map

    Phi(c) = B^{-1}(d - Q(c,c)).

If ``nu = ||B^{-1}||``, ``q`` is a bilinear norm bound for ``Q`` and
``r = 2 nu ||d||``, the reader gives

    ||Phi(c)|| <= r/2 + nu q r^2,
    Lip(Phi) <= 2 nu q r,

with the sufficient condition

    8 nu^2 q ||d|| <= 1,

which makes the radius-r ball invariant and gives ``Lip(Phi) <= 1/2``.

This module turns exactly that *local five-moment coefficient solve* into a
fail-closed executable diagnostic.  The public entry points accept only a typed
backend plus a candidate object.  They do not accept caller-entered norms,
gains, success booleans, residuals, pressure, forcing, damping factors or
scientific thresholds.  The backend supplies the actual five-dimensional
``B``, ``Q`` and ``d`` payload; this module recomputes ``nu`` and a rigorous
finite-dimensional bilinear upper bound for ``q`` before deciding whether a
bounded fixed-point solve is allowed.

This is not the full Navier--Stokes correction-cycle gain and it does not turn
current strict-inner radial stress into an authorized NS correction target.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import inspect
import json
import math
from pathlib import Path
from typing import Any, Protocol

import numpy as np


TASK = "KOKUNO-A3-MOMENT-CORRECTION-CONTRACTION-GAIN-088"
SCHEMA = "kokuno-a3-five-moment-correction-gain-v1"
PARENT_AGENT3_PR = 920
PARENT_AGENT3_HEAD = "3ca100009e770f2aecb5859c721bdd99ffb8552c"
SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_READER_TEX_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_TEX_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"

MOMENT_DIMENSION = 5
SOURCE_SMALLNESS_LIMIT = 1.0
SOURCE_LIPSCHITZ_LIMIT = 0.5
FIXED_POINT_MAX_ITERATIONS = 64
FIXED_POINT_RELATIVE_TOLERANCE = 1.0e-12
MATRIX_CONDITION_REJECT_AT = 1.0e12
FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5


class MomentCorrectionGainError(RuntimeError):
    """Raised when typed moment data do not admit the registered solve."""


@dataclass(frozen=True)
class CandidateIdentity:
    candidate_id: str
    candidate_sha256: str
    evidence_kind: str

    def __post_init__(self) -> None:
        if not self.candidate_id or not self.candidate_sha256:
            raise ValueError("candidate identity must be nonempty")
        if self.evidence_kind not in {"mechanics-only", "repository-candidate"}:
            raise ValueError("evidence_kind must be mechanics-only or repository-candidate")


@dataclass(frozen=True)
class FiveMomentCorrectionSystem:
    """Checksumable five-dimensional B/Q/d payload supplied by a typed backend."""

    source_candidate_id: str
    source_candidate_sha256: str
    fit_partition_id: str
    fit_sample_ids: tuple[str, ...]
    matrix_B: tuple[tuple[float, ...], ...]
    bilinear_Q: tuple[tuple[tuple[float, ...], ...], ...]
    discrepancy_d: tuple[float, ...]
    discrepancy_kind: str
    discrepancy_from_complete_ns_defect: bool
    residual_as_forcing_shortcut_used: bool = False

    def __post_init__(self) -> None:
        if not self.source_candidate_id or not self.source_candidate_sha256:
            raise ValueError("source candidate identity must be nonempty")
        if not self.fit_partition_id:
            raise ValueError("fit_partition_id must be nonempty")
        if not self.fit_sample_ids or any(not item for item in self.fit_sample_ids):
            raise ValueError("fit_sample_ids must contain nonempty ids")
        if len(set(self.fit_sample_ids)) != len(self.fit_sample_ids):
            raise ValueError("fit_sample_ids must be unique")
        if not self.discrepancy_kind:
            raise ValueError("discrepancy_kind must be nonempty")


class FiveMomentCorrectionBackend(Protocol):
    def candidate_identity(self, candidate: Any) -> CandidateIdentity: ...

    def five_moment_correction_system(
        self, candidate: Any
    ) -> FiveMomentCorrectionSystem: ...


@dataclass(frozen=True)
class MomentCorrectionGainDiagnostic:
    candidate: CandidateIdentity
    system_sha256: str
    inverse_norm_nu: float
    matrix_condition_2: float
    bilinear_norm_upper_q: float
    discrepancy_norm: float
    registered_radius: float
    image_norm_upper: float
    lipschitz_upper: float
    source_smallness_parameter: float
    source_smallness_condition_met: bool
    self_map_certified: bool
    contraction_certified: bool
    correction_needed: bool
    rejection_reasons: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class FiveMomentCorrectionSolveReceipt:
    diagnostic: MomentCorrectionGainDiagnostic
    coefficient_vector: tuple[float, ...]
    iterations: int
    equation_residual_norm: float
    final_fixed_point_increment_norm: float
    coefficient_norm: float
    correction_nontrivial: bool
    local_moment_fixed_point_converged: bool
    local_moment_correction_accepted: bool
    complete_ns_correction_authorized: bool
    analytic_full_ns_cycle_gain_available: bool
    heldout_ns_residual_assessed: bool
    pde_validated: bool

    def to_dict(self) -> dict[str, object]:
        return {
            **asdict(self),
            "truth_boundary": truth_boundary(),
        }


def _require_backend(backend: Any) -> None:
    missing = [
        name
        for name in ("candidate_identity", "five_moment_correction_system")
        if not callable(getattr(backend, name, None))
    ]
    if missing:
        raise MomentCorrectionGainError(
            "backend is missing typed moment-correction methods: "
            + ", ".join(missing)
        )


def _array_payload(
    system: FiveMomentCorrectionSystem,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    B = np.asarray(system.matrix_B, dtype=float)
    Q = np.asarray(system.bilinear_Q, dtype=float)
    d = np.asarray(system.discrepancy_d, dtype=float)
    if B.shape != (MOMENT_DIMENSION, MOMENT_DIMENSION):
        raise MomentCorrectionGainError("B must be exactly 5x5")
    if Q.shape != (MOMENT_DIMENSION, MOMENT_DIMENSION, MOMENT_DIMENSION):
        raise MomentCorrectionGainError("Q must be exactly 5x5x5")
    if d.shape != (MOMENT_DIMENSION,):
        raise MomentCorrectionGainError("d must have exactly five components")
    if not np.all(np.isfinite(B)) or not np.all(np.isfinite(Q)) or not np.all(np.isfinite(d)):
        raise MomentCorrectionGainError("B/Q/d must contain only finite values")
    return B, Q, d


def _canonical_system_sha256(
    identity: CandidateIdentity,
    system: FiveMomentCorrectionSystem,
) -> str:
    payload = {
        "schema": SCHEMA,
        "candidate": asdict(identity),
        "system": asdict(system),
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_tex_blob": SOURCE_READER_TEX_BLOB,
        "moment_dimension": MOMENT_DIMENSION,
        "norm_convention": "Euclidean vectors / spectral matrix inverse / slice-spectral bilinear upper bound",
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _bilinear_norm_upper(Q: np.ndarray) -> float:
    """Rigorous Euclidean bilinear bound from output-slice spectral norms.

    If ``Q_i`` is the matrix for output component i, then
    ``|x^T Q_i y| <= ||Q_i||_2 ||x||_2 ||y||_2``.  Therefore
    ``sqrt(sum_i ||Q_i||_2^2)`` is a valid upper bound for the vector-valued
    bilinear operator norm.  It is conservative but directly computable.
    """

    slice_norms = np.array([np.linalg.norm(Q[i], ord=2) for i in range(Q.shape[0])])
    return float(np.linalg.norm(slice_norms, ord=2))


def _typed_payload(
    backend: FiveMomentCorrectionBackend,
    candidate: Any,
) -> tuple[CandidateIdentity, FiveMomentCorrectionSystem, np.ndarray, np.ndarray, np.ndarray]:
    _require_backend(backend)
    identity = backend.candidate_identity(candidate)
    if not isinstance(identity, CandidateIdentity):
        raise MomentCorrectionGainError("backend must return CandidateIdentity")
    system = backend.five_moment_correction_system(candidate)
    if not isinstance(system, FiveMomentCorrectionSystem):
        raise MomentCorrectionGainError("backend must return FiveMomentCorrectionSystem")
    if system.source_candidate_id != identity.candidate_id:
        raise MomentCorrectionGainError("five-moment system candidate id mismatch")
    if system.source_candidate_sha256 != identity.candidate_sha256:
        raise MomentCorrectionGainError("five-moment system candidate hash mismatch")
    if system.residual_as_forcing_shortcut_used:
        raise MomentCorrectionGainError("forbidden residual-as-forcing shortcut in moment provenance")
    B, Q, d = _array_payload(system)
    return identity, system, B, Q, d


def assess_five_moment_correction_gain(
    backend: FiveMomentCorrectionBackend,
    candidate: Any,
) -> MomentCorrectionGainDiagnostic:
    """Recompute the public-source local contraction bound from typed B/Q/d.

    No gain/norm/smallness scalar can be supplied by the caller.  A conservative
    finite-dimensional bilinear norm is recomputed from the actual Q tensor.
    """

    identity, system, B, Q, d = _typed_payload(backend, candidate)
    try:
        B_inv = np.linalg.inv(B)
    except np.linalg.LinAlgError as exc:
        raise MomentCorrectionGainError("five-moment matrix B is singular") from exc

    nu = float(np.linalg.norm(B_inv, ord=2))
    condition = float(np.linalg.cond(B, p=2))
    q_upper = _bilinear_norm_upper(Q)
    d_norm = float(np.linalg.norm(d, ord=2))
    radius = 2.0 * nu * d_norm
    image_upper = nu * d_norm + nu * q_upper * radius * radius
    lipschitz_upper = 2.0 * nu * q_upper * radius
    source_smallness = 8.0 * nu * nu * q_upper * d_norm

    finite_scalars = (
        nu,
        condition,
        q_upper,
        d_norm,
        radius,
        image_upper,
        lipschitz_upper,
        source_smallness,
    )
    if not all(math.isfinite(value) and value >= 0.0 for value in finite_scalars):
        raise MomentCorrectionGainError("computed contraction diagnostics are non-finite")

    correction_needed = d_norm > 0.0
    source_condition = source_smallness <= SOURCE_SMALLNESS_LIMIT
    if radius == 0.0:
        self_map = image_upper == 0.0
    else:
        self_map = image_upper <= radius * (1.0 + 1.0e-12)
    contraction = lipschitz_upper <= SOURCE_LIPSCHITZ_LIMIT * (1.0 + 1.0e-12)

    rejection_reasons: list[str] = []
    if not correction_needed:
        rejection_reasons.append("zero_moment_discrepancy")
    if condition > MATRIX_CONDITION_REJECT_AT:
        rejection_reasons.append("moment_matrix_conditioning_rejected")
    if not source_condition:
        rejection_reasons.append("kokuno_source_smallness_condition_failed")
    if not self_map:
        rejection_reasons.append("registered_ball_not_certified_self_mapping")
    if not contraction:
        rejection_reasons.append("fixed_point_lipschitz_bound_exceeds_one_half")

    return MomentCorrectionGainDiagnostic(
        candidate=identity,
        system_sha256=_canonical_system_sha256(identity, system),
        inverse_norm_nu=nu,
        matrix_condition_2=condition,
        bilinear_norm_upper_q=q_upper,
        discrepancy_norm=d_norm,
        registered_radius=radius,
        image_norm_upper=image_upper,
        lipschitz_upper=lipschitz_upper,
        source_smallness_parameter=source_smallness,
        source_smallness_condition_met=source_condition,
        self_map_certified=self_map,
        contraction_certified=contraction,
        correction_needed=correction_needed,
        rejection_reasons=tuple(rejection_reasons),
    )


def _quadratic(Q: np.ndarray, c: np.ndarray) -> np.ndarray:
    return np.einsum("ijk,j,k->i", Q, c, c, optimize=True)


def solve_five_moment_correction(
    backend: FiveMomentCorrectionBackend,
    candidate: Any,
) -> FiveMomentCorrectionSolveReceipt:
    """Run the bounded fixed-point iteration only after the gain firewall passes."""

    diagnostic = assess_five_moment_correction_gain(backend, candidate)
    if diagnostic.rejection_reasons:
        raise MomentCorrectionGainError(
            "typed five-moment correction rejected before iteration: "
            + ", ".join(diagnostic.rejection_reasons)
        )

    identity, system, B, Q, d = _typed_payload(backend, candidate)
    if _canonical_system_sha256(identity, system) != diagnostic.system_sha256:
        raise MomentCorrectionGainError("typed five-moment payload changed after gain assessment")

    c = np.zeros(MOMENT_DIMENSION, dtype=float)
    final_increment = math.inf
    equation_residual = math.inf
    converged = False
    iterations = 0
    for iterations in range(1, FIXED_POINT_MAX_ITERATIONS + 1):
        next_c = np.linalg.solve(B, d - _quadratic(Q, c))
        if not np.all(np.isfinite(next_c)):
            raise MomentCorrectionGainError("fixed-point iterate became non-finite")
        if float(np.linalg.norm(next_c, ord=2)) > diagnostic.registered_radius * (1.0 + 1.0e-10):
            raise MomentCorrectionGainError("fixed-point iterate escaped registered contraction ball")

        final_increment = float(np.linalg.norm(next_c - c, ord=2))
        c = next_c
        equation_residual = float(
            np.linalg.norm(B @ c + _quadratic(Q, c) - d, ord=2)
        )
        scale_c = max(1.0, float(np.linalg.norm(c, ord=2)))
        scale_d = max(1.0, diagnostic.discrepancy_norm)
        if (
            final_increment <= FIXED_POINT_RELATIVE_TOLERANCE * scale_c
            and equation_residual <= FIXED_POINT_RELATIVE_TOLERANCE * scale_d
        ):
            converged = True
            break

    if not converged:
        raise MomentCorrectionGainError(
            f"fixed-point iteration did not converge within {FIXED_POINT_MAX_ITERATIONS} iterations"
        )

    coefficient_norm = float(np.linalg.norm(c, ord=2))
    nontrivial = coefficient_norm > 0.0
    if not nontrivial:
        raise MomentCorrectionGainError("converged moment correction is trivial")

    return FiveMomentCorrectionSolveReceipt(
        diagnostic=diagnostic,
        coefficient_vector=tuple(float(value) for value in c),
        iterations=iterations,
        equation_residual_norm=equation_residual,
        final_fixed_point_increment_norm=final_increment,
        coefficient_norm=coefficient_norm,
        correction_nontrivial=nontrivial,
        local_moment_fixed_point_converged=True,
        local_moment_correction_accepted=True,
        complete_ns_correction_authorized=False,
        analytic_full_ns_cycle_gain_available=False,
        heldout_ns_residual_assessed=False,
        pde_validated=False,
    )


def truth_boundary() -> dict[str, object]:
    assess_signature = inspect.signature(assess_five_moment_correction_gain)
    solve_signature = inspect.signature(solve_five_moment_correction)
    forbidden = {
        "B", "Q", "d", "nu", "q", "gain", "radius", "lipschitz",
        "smallness", "success", "residual", "defect", "pressure", "forcing",
        "damping", "alpha", "threshold", "scientific_threshold", "viscosity",
    }
    return {
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_date": SOURCE_READER_DATE,
        "source_reader_tex_path": SOURCE_READER_TEX_PATH,
        "source_reader_tex_blob": SOURCE_READER_TEX_BLOB,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "moment_dimension": MOMENT_DIMENSION,
        "source_equation": "B c + Q(c,c) = d",
        "source_fixed_point": "Phi(c)=B^{-1}(d-Q(c,c))",
        "source_radius_rule": "r=2*nu*||d||",
        "source_image_bound": "||Phi(c)||<=r/2+nu*q*r^2",
        "source_lipschitz_bound": "Lip(Phi)<=2*nu*q*r",
        "source_smallness_rule": "8*nu^2*q*||d||<=1",
        "source_lipschitz_consequence": "Lip(Phi)<=1/2",
        "bilinear_bound_realization": "sqrt(sum_i ||Q_i||_2^2)",
        "caller_supplied_gain_allowed": False,
        "caller_supplied_norm_bound_allowed": False,
        "caller_supplied_success_boolean_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(assess_signature.parameters)
        and forbidden.isdisjoint(solve_signature.parameters),
        "divergent_local_moment_correction_rejected_before_iteration": True,
        "local_five_moment_analytic_contraction_gain_available": True,
        "analytic_full_ns_cycle_gain_available": False,
        "current_strict_inner_radial_stress_authorized_as_ns_correction": False,
        "complete_ns_defect_required_for_final_ns_cycle": True,
        "complete_ns_correction_velocity_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }


@dataclass(frozen=True)
class _MechanicsCandidate:
    candidate_id: str
    candidate_sha256: str


class _MechanicsBackend:
    def __init__(self, *, quadratic_scale: float) -> None:
        self.quadratic_scale = float(quadratic_scale)

    def candidate_identity(self, candidate: _MechanicsCandidate) -> CandidateIdentity:
        return CandidateIdentity(
            candidate_id=candidate.candidate_id,
            candidate_sha256=candidate.candidate_sha256,
            evidence_kind="mechanics-only",
        )

    def five_moment_correction_system(
        self, candidate: _MechanicsCandidate
    ) -> FiveMomentCorrectionSystem:
        B = np.eye(MOMENT_DIMENSION, dtype=float)
        Q = np.zeros((MOMENT_DIMENSION, MOMENT_DIMENSION, MOMENT_DIMENSION), dtype=float)
        for i in range(MOMENT_DIMENSION):
            Q[i, i, i] = self.quadratic_scale
        d = np.array([0.05, -0.04, 0.03, -0.02, 0.01], dtype=float)
        return FiveMomentCorrectionSystem(
            source_candidate_id=candidate.candidate_id,
            source_candidate_sha256=candidate.candidate_sha256,
            fit_partition_id="mechanics-held-in",
            fit_sample_ids=("m0", "m1", "m2", "m3", "m4"),
            matrix_B=tuple(tuple(float(v) for v in row) for row in B),
            bilinear_Q=tuple(
                tuple(tuple(float(v) for v in row) for row in plane) for plane in Q
            ),
            discrepancy_d=tuple(float(v) for v in d),
            discrepancy_kind="mechanics-only-five-moment-discrepancy",
            discrepancy_from_complete_ns_defect=False,
            residual_as_forcing_shortcut_used=False,
        )


def build_mechanics_report() -> dict[str, object]:
    """Deterministic positive/negative mechanics receipt, never NS evidence."""

    candidate = _MechanicsCandidate(
        candidate_id="mechanics-five-moment-candidate",
        candidate_sha256=hashlib.sha256(b"mechanics-five-moment-candidate").hexdigest(),
    )
    accepted = solve_five_moment_correction(
        _MechanicsBackend(quadratic_scale=0.10), candidate
    )
    rejected = assess_five_moment_correction_gain(
        _MechanicsBackend(quadratic_scale=20.0), candidate
    )
    return {
        "schema": SCHEMA,
        "task": TASK,
        "accepted_mechanics": accepted.to_dict(),
        "rejected_mechanics": rejected.to_dict(),
        "rejected_mechanics_would_be_blocked": bool(rejected.rejection_reasons),
        "candidate_residual_evidence": False,
        "truth_boundary": truth_boundary(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    report = build_mechanics_report()
    encoded = json.dumps(report, indent=2, sort_keys=True)
    if args.out is None:
        print(encoded)
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(encoded + "\n", encoding="utf-8")
        print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
