"""Bind the current Agent-1 Xi five-moment handoff into Agent-3 correction input.

This module consumes the *current-lineage* candidate-side PA.15 discrepancy
materialized by Agent 1 PR #947.  It performs only the public PA.16 row
transformation

    (M, I, J, S, C_p)
      -> (M, J - 4 eta I, I, S - 8 eta M, C_p)

and records a checksum-bound typed receipt for the Agent-3 correction lane.
It deliberately does not execute the historical PA.16 inverse, choose T_sh,
construct a velocity correction, or authorize the gain-gated full-NS cycle.
The upstream discrepancy is a real current candidate-side source-moment
mismatch, but it is not a discrepancy derived from a complete Navier--Stokes
defect: the current Kokuno lane still lacks the global Cartesian leading join,
matched pressure, and preregistered restricted forcing.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import inspect
import json
import math
from typing import Any, Mapping, Protocol

import numpy as np


TASK = "KOKUNO-A3-CURRENT-XI-MOMENT-DISCREPANCY-BRIDGE-091"
SCHEMA = "kokuno-a3-current-xi-moment-discrepancy-bridge-v1"
PARENT_AGENT3_PR = 942
PARENT_AGENT3_HEAD = "728b58cf834a45960c8190d5a3898ba7f6ead273"
UPSTREAM_AGENT1_PR = 947
UPSTREAM_AGENT1_HEAD = "2797eed8e3b9b3374c676d8f517f05e3f0fc3e6b"
UPSTREAM_AGENT1_SOURCE_PATH = (
    "src/openai_ns_reconstruction/kokuno_pa10_actual_xi_prefix_moments.py"
)
UPSTREAM_AGENT1_SOURCE_BLOB = "92dce65c9a8793e06497a7e7347c832861032238"
UPSTREAM_AGENT1_SCHEMA = "kokuno-pa10-actual-xi-prefix-moments-v1"
SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_TEX_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_READER_TEX_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_READER_DATE = "2026-09-09"
FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5
NONTRIVIAL_DISCREPANCY_FLOOR = 1.0e-14


class CurrentXiMomentBridgeError(RuntimeError):
    """Raised when the typed current-Xi handoff does not satisfy the contract."""


@dataclass(frozen=True)
class UpstreamA1Identity:
    pr_number: int
    exact_head: str
    source_path: str
    source_blob: str
    report_schema: str


class CurrentXiMomentHandoffBackend(Protocol):
    def upstream_identity(self) -> UpstreamA1Identity: ...

    def current_xi_moment_report(self) -> Mapping[str, Any]: ...


@dataclass(frozen=True)
class CurrentXiMomentSample:
    eta: float
    ell_i: float
    G_i: float
    incoming_scaled_discrepancy: tuple[float, float, float, float, float]
    pa16_row_transformed_discrepancy: tuple[float, float, float, float, float]
    raw_l2_norm: float
    pa16_transformed_l2_norm: float
    correction_needed: bool


@dataclass(frozen=True)
class CurrentXiMomentBridgeReceipt:
    upstream: UpstreamA1Identity
    upstream_semantic_sha256: str
    upstream_report_sha256: str
    outer_schedule_sha256: str
    samples: tuple[CurrentXiMomentSample, ...]
    max_raw_l2_norm: float
    max_pa16_transformed_l2_norm: float
    any_correction_needed: bool
    handoff_sha256: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "upstream": asdict(self.upstream),
            "upstream_semantic_sha256": self.upstream_semantic_sha256,
            "upstream_report_sha256": self.upstream_report_sha256,
            "outer_schedule_sha256": self.outer_schedule_sha256,
            "samples": [asdict(sample) for sample in self.samples],
            "max_raw_l2_norm": self.max_raw_l2_norm,
            "max_pa16_transformed_l2_norm": self.max_pa16_transformed_l2_norm,
            "any_correction_needed": self.any_correction_needed,
            "handoff_sha256": self.handoff_sha256,
            "truth_boundary": truth_boundary(),
        }


def _canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256_json(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _require_backend(backend: Any) -> None:
    missing = [
        name
        for name in ("upstream_identity", "current_xi_moment_report")
        if not callable(getattr(backend, name, None))
    ]
    if missing:
        raise CurrentXiMomentBridgeError(
            "backend is missing typed current-Xi methods: " + ", ".join(missing)
        )


def _require_exact_upstream(identity: UpstreamA1Identity) -> None:
    if not isinstance(identity, UpstreamA1Identity):
        raise CurrentXiMomentBridgeError("backend must return UpstreamA1Identity")
    expected = UpstreamA1Identity(
        pr_number=UPSTREAM_AGENT1_PR,
        exact_head=UPSTREAM_AGENT1_HEAD,
        source_path=UPSTREAM_AGENT1_SOURCE_PATH,
        source_blob=UPSTREAM_AGENT1_SOURCE_BLOB,
        report_schema=UPSTREAM_AGENT1_SCHEMA,
    )
    if identity != expected:
        raise CurrentXiMomentBridgeError("current-Xi handoff is not the pinned Agent-1 #947 artifact")


def _require_mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise CurrentXiMomentBridgeError(f"{name} must be a mapping")
    return value


def _finite_vector(value: Any, name: str, *, size: int | None = None) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if size is not None and array.shape != (size,):
        raise CurrentXiMomentBridgeError(f"{name} must have shape ({size},)")
    if np.any(~np.isfinite(array)):
        raise CurrentXiMomentBridgeError(f"{name} must contain only finite values")
    return array


def _validate_semantic_sha(report: Mapping[str, Any]) -> str:
    source = _require_mapping(report.get("source"), "report.source")
    if report.get("schema") != UPSTREAM_AGENT1_SCHEMA:
        raise CurrentXiMomentBridgeError("unexpected upstream report schema")
    if source.get("repository") != SOURCE_READER_REPO:
        raise CurrentXiMomentBridgeError("upstream source repository drifted")
    if source.get("commit") != SOURCE_READER_HEAD:
        raise CurrentXiMomentBridgeError("upstream corrected-reader commit drifted")
    if source.get("path") != SOURCE_READER_TEX_PATH:
        raise CurrentXiMomentBridgeError("upstream corrected-reader path drifted")
    if source.get("release_date") != SOURCE_READER_DATE:
        raise CurrentXiMomentBridgeError("upstream corrected-reader date drifted")

    semantic_payload = {
        "schema": report.get("schema"),
        "source_commit": source.get("commit"),
        "source_formulas": report.get("source_formulas"),
        "numerical_realization": report.get("numerical_realization"),
        "configuration": report.get("configuration"),
    }
    semantic_sha = _sha256_json(semantic_payload)
    if report.get("semantic_sha256") != semantic_sha:
        raise CurrentXiMomentBridgeError("upstream semantic SHA does not replay")
    return semantic_sha


def _validate_truth_boundary(report: Mapping[str, Any]) -> None:
    truth = _require_mapping(report.get("truth_boundary"), "report.truth_boundary")
    required_true = (
        "current_lineage_profile_0_to_Xi_executable",
        "public_five_prefix_moment_map_executable",
        "public_PA15_scaling_executable",
        "candidate_side_upstream_five_moment_discrepancy_at_Xi_materialized",
        "candidate_side_PA16_input_tuple_materialized",
    )
    required_false = (
        "source_prepared_appendixA_pressure_stress_materialized",
        "source_prepared_upstream_five_moment_discrepancy_materialized",
        "source_T_sh_lower_bound_verified",
        "five_moment_repair_applied",
        "inner_to_outer_join_completed",
        "outer_global_leading_velocity_materialized",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "heldout_ns_residual_assessed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    )
    for key in required_true:
        if truth.get(key) is not True:
            raise CurrentXiMomentBridgeError(f"upstream truth boundary lost required true flag: {key}")
    for key in required_false:
        if truth.get(key) is not False:
            raise CurrentXiMomentBridgeError(f"upstream truth boundary attempted promotion: {key}")


def _pa16_row_transform(discrepancy: np.ndarray, eta: float) -> np.ndarray:
    """Apply only the public PA.16 row operations, preserving exact target data."""
    M, I, J, S, C_p = (float(v) for v in discrepancy)
    return np.asarray(
        [M, J - 4.0 * eta * I, I, S - 8.0 * eta * M, C_p],
        dtype=float,
    )


def materialize_current_xi_moment_bridge(
    backend: CurrentXiMomentHandoffBackend,
) -> CurrentXiMomentBridgeReceipt:
    """Consume the exact current A1 Xi report and expose typed PA.16 input rows.

    The only public input is a typed backend.  In particular there is no public
    residual, defect, moment vector, pressure, forcing, gain, damping, T_sh,
    correction, viscosity, stage budget or scientific-threshold argument.
    """
    _require_backend(backend)
    identity = backend.upstream_identity()
    _require_exact_upstream(identity)
    report = _require_mapping(backend.current_xi_moment_report(), "upstream report")
    semantic_sha = _validate_semantic_sha(report)
    _validate_truth_boundary(report)

    eta_probe = _finite_vector(report.get("eta_probe"), "eta_probe")
    raw_rows = np.asarray(report.get("incoming_PA15_discrepancy"), dtype=float)
    pa16_inputs = report.get("pa16_inputs")
    if raw_rows.shape != (eta_probe.size, 5) or np.any(~np.isfinite(raw_rows)):
        raise CurrentXiMomentBridgeError("incoming_PA15_discrepancy must be finite N x 5")
    if not isinstance(pa16_inputs, list) or len(pa16_inputs) != eta_probe.size:
        raise CurrentXiMomentBridgeError("pa16_inputs must align one-for-one with eta_probe")
    if eta_probe.size == 0:
        raise CurrentXiMomentBridgeError("current-Xi report must contain at least one eta probe")

    samples: list[CurrentXiMomentSample] = []
    schedule_sha: str | None = None
    for index, eta_value in enumerate(eta_probe):
        item = _require_mapping(pa16_inputs[index], f"pa16_inputs[{index}]")
        eta = float(eta_value)
        if not math.isclose(float(item.get("eta")), eta, rel_tol=0.0, abs_tol=0.0):
            raise CurrentXiMomentBridgeError("PA.16 eta row does not match eta_probe")
        ell_i = float(item.get("ell_i"))
        G_i = float(item.get("G_i"))
        if not math.isfinite(ell_i) or not math.isfinite(G_i):
            raise CurrentXiMomentBridgeError("ell_i/G_i must be finite")
        row = _finite_vector(
            item.get("incoming_scaled_discrepancy"),
            f"pa16_inputs[{index}].incoming_scaled_discrepancy",
            size=5,
        )
        if not np.array_equal(row, raw_rows[index]):
            raise CurrentXiMomentBridgeError(
                "top-level incoming discrepancy does not match PA.16 handoff row"
            )
        if item.get("source_T_sh_lower_bound_verified") is not False:
            raise CurrentXiMomentBridgeError(
                "upstream handoff may not promote the unverified T_sh lower bound"
            )
        item_schedule_sha = item.get("outer_schedule_sha256")
        if not isinstance(item_schedule_sha, str) or len(item_schedule_sha) != 64:
            raise CurrentXiMomentBridgeError("outer schedule SHA must be a 64-character digest")
        if schedule_sha is None:
            schedule_sha = item_schedule_sha
        elif item_schedule_sha != schedule_sha:
            raise CurrentXiMomentBridgeError("outer schedule identity changed across eta rows")

        transformed = _pa16_row_transform(row, eta)
        raw_norm = float(np.linalg.norm(row, ord=2))
        transformed_norm = float(np.linalg.norm(transformed, ord=2))
        samples.append(
            CurrentXiMomentSample(
                eta=eta,
                ell_i=ell_i,
                G_i=G_i,
                incoming_scaled_discrepancy=tuple(float(v) for v in row),
                pa16_row_transformed_discrepancy=tuple(float(v) for v in transformed),
                raw_l2_norm=raw_norm,
                pa16_transformed_l2_norm=transformed_norm,
                correction_needed=raw_norm > NONTRIVIAL_DISCREPANCY_FLOOR,
            )
        )

    if schedule_sha is None:
        raise CurrentXiMomentBridgeError("outer schedule identity is missing")
    report_sha = _sha256_json(report)
    handoff_payload = {
        "schema": SCHEMA,
        "upstream": asdict(identity),
        "upstream_semantic_sha256": semantic_sha,
        "upstream_report_sha256": report_sha,
        "outer_schedule_sha256": schedule_sha,
        "samples": [asdict(sample) for sample in samples],
    }
    handoff_sha = _sha256_json(handoff_payload)
    max_raw = max(sample.raw_l2_norm for sample in samples)
    max_transformed = max(sample.pa16_transformed_l2_norm for sample in samples)
    return CurrentXiMomentBridgeReceipt(
        upstream=identity,
        upstream_semantic_sha256=semantic_sha,
        upstream_report_sha256=report_sha,
        outer_schedule_sha256=schedule_sha,
        samples=tuple(samples),
        max_raw_l2_norm=max_raw,
        max_pa16_transformed_l2_norm=max_transformed,
        any_correction_needed=any(sample.correction_needed for sample in samples),
        handoff_sha256=handoff_sha,
    )


def truth_boundary() -> dict[str, Any]:
    signature = inspect.signature(materialize_current_xi_moment_bridge)
    forbidden = {
        "residual",
        "defect",
        "discrepancy",
        "moment",
        "mean",
        "stress",
        "pressure",
        "forcing",
        "gain",
        "damping",
        "alpha",
        "T_sh",
        "correction",
        "viscosity",
        "max_stages",
        "threshold",
        "scientific_threshold",
    }
    return {
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "upstream_agent1_pr": UPSTREAM_AGENT1_PR,
        "upstream_agent1_head": UPSTREAM_AGENT1_HEAD,
        "upstream_agent1_source_blob": UPSTREAM_AGENT1_SOURCE_BLOB,
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_tex_blob": SOURCE_READER_TEX_BLOB,
        "source_reader_date": SOURCE_READER_DATE,
        "public_pa16_row_transform": "(M,I,J,S,C_p)->(M,J-4*eta*I,I,S-8*eta*M,C_p)",
        "current_candidate_side_xi_source_moment_discrepancy_ingested": True,
        "current_candidate_side_pa16_rows_materialized": True,
        "historical_pa16_solver_duplicated": False,
        "source_T_sh_lower_bound_verified": False,
        "pa16_repair_applied": False,
        "discrepancy_from_complete_ns_defect": False,
        "authorized_for_gain_gated_ns_stage": False,
        "current_real_ns_defect_five_moment_discrepancy_materialized": False,
        "current_real_ns_correction_velocity_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }
