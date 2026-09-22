"""Typed Kokuno RF44--RF49 post-update mean-defect recomputation.

Stacks on A3 #1134.  The corrected 2026-09-09 reader requires the nonlinear
mean defect to be recomputed after the RF31/RF34--RF39 correction; closing the
five base-linear rows is not a zero-defect certificate.  This module accepts
only candidate-bound pre/post fixed-Q states and the exact #1134 correction.
"""
from __future__ import annotations

import argparse, hashlib, json, math
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Protocol
import numpy as np

from .kokuno_rf30_rf31_typed_mean_correction import CandidateIdentity
from .kokuno_rf34_rf39_compact_mean_correction import (
    RF34RF39Backend, RF34RF39ContractError, RF34RF39CorrectionReceipt,
    _MechanicsBackend, materialize_rf34_rf39_compact_mean_correction,
)

TASK = "KOKUNO-A3-RF44-RF49-POSTUPDATE-RECOMPUTE-124"
SCHEMA = "kokuno-a3-rf44-rf49-postupdate-recompute-v1"
PARENT_AGENT3_PR = 1134
PARENT_AGENT3_HEAD = "ecbc7e3702e366a5f355800dd662b11766273a41"
PARENT_SOURCE_BLOB = "392975eafd3787c91e2dc82ad9b735567513271d"
SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_DATE = "2026-09-09"
SOURCE_READER_TEX_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_TEX_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_FORMULAS = ("RF44", "RF45", "RF46", "RF47", "RF48", "RF49")
FINAL_NORMALIZED_MOMENTUM_GATE = 1e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1e-5
STATE_FIXED_ATOL = 2e-11
INCREMENT_MATCH_ATOL = 2e-10
RF47_REL_TOL = 2e-7
DEFECT_NORM_FLOOR = 1e-18


class RF44RF49ContractError(RF34RF39ContractError):
    pass


@dataclass(frozen=True)
class RF44MeanState:
    source_candidate_id: str
    source_candidate_sha256: str
    source_chart_id: str
    source_chart_sha256: str
    radius_R: tuple[float, ...]
    base_b: tuple[float, ...]
    base_G: tuple[float, ...]
    base_V: tuple[float, ...]
    beta: tuple[float, ...]
    gamma: tuple[float, ...]
    v: tuple[float, ...]
    mean_W_rr: tuple[float, ...]
    mean_W_zr: tuple[float, ...]
    mean_W_thetatheta: tuple[float, ...]
    mean_W_ztheta: tuple[float, ...]
    mean_W_zz: tuple[float, ...]
    tstar_beta: tuple[float, ...]
    dR_rr_bundle: tuple[float, ...]
    dZ_zr_bundle: tuple[float, ...]
    epsilon_laplacian_beta: tuple[float, ...]
    actual_candidate_recomputed: bool
    normalized_haar_mean_used: bool
    source_fixed_q_chart_used: bool
    operator_terms_recomputed_from_state: bool
    correction_applied: bool
    heldout_samples_used_to_construct_state: bool = False
    surrogate_defect_used: bool = False
    residual_as_forcing_shortcut_used: bool = False


@dataclass(frozen=True)
class RF44CorrectionIncrement:
    source_candidate_id: str
    source_candidate_sha256: str
    source_chart_id: str
    source_chart_sha256: str
    parent_correction_sha256: str
    delta_beta: tuple[float, ...]
    delta_v: tuple[float, ...]
    delta_gamma: tuple[float, ...]
    induced_beta_reconstructed_from_solenoidal_correction: bool
    source_solenoidal_realization_verified: bool
    compact_support_verified: bool
    heldout_samples_used_to_construct_increment: bool = False
    surrogate_increment_used: bool = False


class RF44RF49Backend(RF34RF39Backend, Protocol):
    def rf44_pre_update_state(self, candidate: Any, correction: RF34RF39CorrectionReceipt) -> RF44MeanState: ...
    def rf44_post_update_state(self, candidate: Any, correction: RF34RF39CorrectionReceipt) -> RF44MeanState: ...
    def rf44_correction_increment(self, candidate: Any, correction: RF34RF39CorrectionReceipt,
                                  pre_state: RF44MeanState, post_state: RF44MeanState) -> RF44CorrectionIncrement: ...


@dataclass(frozen=True)
class RF44DefectTuple:
    P: float
    J_theta: float
    J_z: float
    radial_source_g_r: tuple[float, ...]
    defect_norm_2: float
    state_sha256: str


@dataclass(frozen=True)
class RF44RF49RecomputeReceipt:
    candidate: CandidateIdentity
    parent_correction_sha256: str
    source_chart_id: str
    source_chart_sha256: str
    pre_defect: RF44DefectTuple
    post_defect: RF44DefectTuple
    radial_remainder_Rg: tuple[float, ...]
    rf47_predicted_post_P_Jtheta_Jz: tuple[float, float, float]
    rf47_closure_max_relative: float
    defect_gain_ratio_2: float
    nonlinear_remainder_nonzero: bool
    repository_candidate_recompute_evidence: bool
    source_rf44_rf49_recompute_materialized: bool
    current_i4_candidate_evidence: bool
    cartesian_correction_velocity_materialized: bool
    complete_ns_defect: bool
    heldout_ns_residual_assessed: bool
    pde_validated: bool
    def to_dict(self) -> dict[str, object]:
        return {**asdict(self), "truth_boundary": truth_boundary()}


def _sha(payload: object) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _trapz(y: np.ndarray, x: np.ndarray) -> float:
    return float(np.sum(0.5 * (y[:-1] + y[1:]) * np.diff(x)))


def _arr(v: tuple[float, ...], n: int, name: str) -> np.ndarray:
    a = np.asarray(v, dtype=float)
    if a.shape != (n,) or not np.all(np.isfinite(a)):
        raise RF44RF49ContractError(f"invalid RF44 vector: {name}")
    return a


def _arrays(s: RF44MeanState) -> dict[str, np.ndarray]:
    R = np.asarray(s.radius_R, dtype=float)
    if R.ndim != 1 or R.size < 3 or np.any(~np.isfinite(R)) or np.any(R <= 0) or np.any(np.diff(R) <= 0):
        raise RF44RF49ContractError("RF44 radius grid must be finite, positive and increasing")
    out = {"R": R}
    for name in ("base_b", "base_G", "base_V", "beta", "gamma", "v", "mean_W_rr", "mean_W_zr",
                 "mean_W_thetatheta", "mean_W_ztheta", "mean_W_zz", "tstar_beta", "dR_rr_bundle",
                 "dZ_zr_bundle", "epsilon_laplacian_beta"):
        out[name] = _arr(getattr(s, name), R.size, name)
    return out


def _check_state(s: RF44MeanState, applied: bool) -> None:
    if s.heldout_samples_used_to_construct_state or s.surrogate_defect_used or s.residual_as_forcing_shortcut_used:
        raise RF44RF49ContractError("RF44 state used held-out/surrogate/residual-as-forcing evidence")
    if not all((s.actual_candidate_recomputed, s.normalized_haar_mean_used,
                s.source_fixed_q_chart_used, s.operator_terms_recomputed_from_state)):
        raise RF44RF49ContractError("RF44 operator/state provenance is incomplete")
    if s.correction_applied != applied:
        raise RF44RF49ContractError("RF44 pre/post correction flag mismatch")


def _defect(s: RF44MeanState) -> tuple[RF44DefectTuple, dict[str, np.ndarray]]:
    a = _arrays(s); R = a["R"]
    b, G, V, beta, gamma, v = (a[k] for k in ("base_b", "base_G", "base_V", "beta", "gamma", "v"))
    g = (-a["tstar_beta"] - a["dR_rr_bundle"] - (2*b*beta + beta*beta + a["mean_W_rr"])/R
         - a["dZ_zr_bundle"] + (2*V*v + v*v + a["mean_W_thetatheta"])/R
         + a["epsilon_laplacian_beta"])
    P = _trapz(g, R)
    Jt = _trapz(R**2 * (G*v + V*gamma + v*gamma + a["mean_W_ztheta"]), R)
    Jz = _trapz(R * (2*G*gamma + gamma*gamma + a["mean_W_zz"]), R) - .5*_trapz(R**2*g, R)
    vals = np.asarray((P, Jt, Jz))
    if not np.all(np.isfinite(vals)):
        raise RF44RF49ContractError("RF44 integrated defect is non-finite")
    return RF44DefectTuple(P, Jt, Jz, tuple(float(x) for x in g), float(np.linalg.norm(vals)), _sha(asdict(s))), a


def materialize_rf44_rf49_postupdate_recompute(backend: RF44RF49Backend, candidate: Any) -> RF44RF49RecomputeReceipt:
    """Apply exact #1134 correction identity and recompute RF44/RF47 defect."""
    for name in ("rf44_pre_update_state", "rf44_post_update_state", "rf44_correction_increment"):
        if not callable(getattr(backend, name, None)):
            raise RF44RF49ContractError(f"backend missing {name}")
    correction = materialize_rf34_rf39_compact_mean_correction(backend, candidate)
    pre = backend.rf44_pre_update_state(candidate, correction)
    post = backend.rf44_post_update_state(candidate, correction)
    if not isinstance(pre, RF44MeanState) or not isinstance(post, RF44MeanState):
        raise RF44RF49ContractError("backend must return typed RF44 states")
    _check_state(pre, False); _check_state(post, True)
    inc = backend.rf44_correction_increment(candidate, correction, pre, post)
    if not isinstance(inc, RF44CorrectionIncrement):
        raise RF44RF49ContractError("backend must return typed RF44 correction increment")
    expected = (correction.candidate.candidate_id, correction.candidate.candidate_sha256,
                correction.source_chart_id, correction.source_chart_sha256)
    for obj in (pre, post, inc):
        if (obj.source_candidate_id, obj.source_candidate_sha256, obj.source_chart_id, obj.source_chart_sha256) != expected:
            raise RF44RF49ContractError("RF44 candidate/chart identity mismatch")
    correction_sha = _sha(asdict(correction))
    if inc.parent_correction_sha256 != correction_sha:
        raise RF44RF49ContractError("RF44 increment is not bound to exact RF34--RF39 correction")
    if inc.heldout_samples_used_to_construct_increment or inc.surrogate_increment_used:
        raise RF44RF49ContractError("RF44 increment used held-out/surrogate data")
    if not all((inc.induced_beta_reconstructed_from_solenoidal_correction,
                inc.source_solenoidal_realization_verified, inc.compact_support_verified)):
        raise RF44RF49ContractError("RF44 solenoidal compact increment provenance is incomplete")

    d0, a0 = _defect(pre); d1, a1 = _defect(post); R = a0["R"]
    if tuple(pre.radius_R) != tuple(post.radius_R):
        raise RF44RF49ContractError("RF44 pre/post grids differ")
    for k in ("base_b", "base_G", "base_V", "mean_W_rr", "mean_W_zr", "mean_W_thetatheta", "mean_W_ztheta", "mean_W_zz"):
        if float(np.max(np.abs(a0[k]-a1[k]))) > STATE_FIXED_ATOL:
            raise RF44RF49ContractError(f"RF44 fixed wave/base drift: {k}")
    B = _arr(inc.delta_beta, R.size, "delta_beta"); u = _arr(inc.delta_v, R.size, "delta_v"); d = _arr(inc.delta_gamma, R.size, "delta_gamma")
    pu = _arr(correction.delta_v_chart, R.size, "RF34 delta_v"); pd = _arr(correction.gamma_d_chart, R.size, "RF34 gamma_d")
    if max(float(np.max(np.abs(u-pu))), float(np.max(np.abs(d-pd)))) > INCREMENT_MATCH_ATOL:
        raise RF44RF49ContractError("RF44 increment does not match exact RF34--RF39 correction")
    for k, delta in (("beta", B), ("v", u), ("gamma", d)):
        if float(np.max(np.abs((a1[k]-a0[k])-delta))) > INCREMENT_MATCH_ATOL:
            raise RF44RF49ContractError(f"RF44 post-pre {k} does not match increment")

    g0 = np.asarray(d0.radial_source_g_r); g1 = np.asarray(d1.radial_source_g_r)
    Rg = g1 - g0 - 2*a0["base_V"]*u/R
    pred = (
        _trapz(Rg, R),
        _trapz(R**2*(a0["gamma"]*u + a0["v"]*d + d*u), R),
        _trapz(R*(2*a0["gamma"]*d + d*d), R) - .5*_trapz(R**2*Rg, R),
    )
    actual = (d1.P, d1.J_theta, d1.J_z)
    rel = max(abs(x-y)/max(abs(y), 1e-14) for x, y in zip(actual, pred))
    if rel > RF47_REL_TOL:
        raise RF44RF49ContractError(f"RF47 post-update identity failed: {rel}")
    gain = d1.defect_norm_2 / max(d0.defect_norm_2, DEFECT_NORM_FLOOR)
    repo = bool(correction.repository_candidate_correction_evidence and correction.candidate.evidence_kind == "repository-candidate")
    return RF44RF49RecomputeReceipt(
        correction.candidate, correction_sha, correction.source_chart_id, correction.source_chart_sha256,
        d0, d1, tuple(float(x) for x in Rg), tuple(float(x) for x in pred), float(rel), float(gain),
        bool(d1.defect_norm_2 > DEFECT_NORM_FLOOR), repo, True, False, False, False, False, False,
    )


def truth_boundary() -> dict[str, object]:
    return {
        "task": TASK, "schema": SCHEMA, "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD, "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD, "source_reader_date": SOURCE_READER_DATE,
        "source_reader_tex_path": SOURCE_READER_TEX_PATH, "source_reader_tex_blob": SOURCE_READER_TEX_BLOB,
        "source_formulas": SOURCE_FORMULAS, "rf44_rf49_postupdate_recompute_contract_materialized": True,
        "source_specific_defect_gain_diagnostic_materialized": True, "new_gain_acceptance_threshold_introduced": False,
        "current_i4_source_chart_backend_materialized": False, "current_i4_rf30_defect_materialized": False,
        "current_i4_rf34_rf39_correction_materialized": False, "current_i4_rf44_rf49_remainder_recomputed": False,
        "cartesian_correction_velocity_materialized": False, "complete_ns_defect": False,
        "real_candidate_finite_correction_cycle_run": False, "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False, "same_protocol_comparable_to_st006": False, "pde_validated": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        "residual_defined_free_forcing_allowed": False, "finite_stage_small_residual_is_blowup_proof": False,
    }


class _RF44MechanicsBackend(_MechanicsBackend):
    """Mechanics-only replay; never repository candidate residual evidence."""
    def _pre(self) -> RF44MeanState:
        s = self.rf30_pre_mean_state(None); R = np.asarray(s.radius_R); z = np.zeros_like(R)
        b = .12 + .015*np.exp(-((R-1)/.4)**2)
        return RF44MeanState(s.source_candidate_id, s.source_candidate_sha256, s.source_chart_id, s.source_chart_sha256,
            s.radius_R, tuple(float(x) for x in b), s.base_G, s.base_V,
            tuple(z), tuple(z), tuple(z), s.mean_W_rr, s.mean_W_zr, s.mean_W_thetatheta, s.mean_W_ztheta, s.mean_W_zz,
            tuple(z), s.dR_mean_W_rr, s.dZ_mean_W_zr, tuple(z), True, True, True, True, False)
    def rf44_pre_update_state(self, candidate: Any, correction: RF34RF39CorrectionReceipt) -> RF44MeanState:
        return self._pre()
    def _increments(self, c: RF34RF39CorrectionReceipt) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        R = np.asarray(self._R); u = np.asarray(c.delta_v_chart); d = np.asarray(c.gamma_d_chart)
        return -.04*np.gradient(R*d, R, edge_order=2), u, d
    def rf44_post_update_state(self, candidate: Any, correction: RF34RF39CorrectionReceipt) -> RF44MeanState:
        s = self._pre(); R = np.asarray(s.radius_R); b = np.asarray(s.base_b); Wrr = np.asarray(s.mean_W_rr)
        B, u, d = self._increments(correction); rr = 2*b*B + B*B + Wrr; z = np.zeros_like(R)
        dB = np.gradient(B, R, edge_order=2); epslap = .03*(np.gradient(dB, R, edge_order=2) + dB/R - B/R**2)
        return replace(s, beta=tuple(B), gamma=tuple(d), v=tuple(u), dR_rr_bundle=tuple(np.gradient(rr, R, edge_order=2)),
                       dZ_zr_bundle=tuple(z), epsilon_laplacian_beta=tuple(epslap), correction_applied=True)
    def rf44_correction_increment(self, candidate: Any, correction: RF34RF39CorrectionReceipt,
                                  pre_state: RF44MeanState, post_state: RF44MeanState) -> RF44CorrectionIncrement:
        B, u, d = self._increments(correction)
        return RF44CorrectionIncrement(pre_state.source_candidate_id, pre_state.source_candidate_sha256,
            pre_state.source_chart_id, pre_state.source_chart_sha256, _sha(asdict(correction)), tuple(B), tuple(u), tuple(d),
            True, True, True)


def build_mechanics_report() -> dict[str, object]:
    r = materialize_rf44_rf49_postupdate_recompute(_RF44MechanicsBackend(), None)
    return {"task": TASK, "schema": SCHEMA, "mechanics_only": True, "candidate_residual_evidence": False, "receipt": r.to_dict()}


def _main() -> None:
    p = argparse.ArgumentParser(); p.add_argument("--out", type=Path, required=True); a = p.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True); a.out.write_text(json.dumps(build_mechanics_report(), indent=2, sort_keys=True)+"\n")


if __name__ == "__main__":
    _main()
