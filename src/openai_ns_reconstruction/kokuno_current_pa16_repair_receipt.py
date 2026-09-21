"""Bind the current-lineage PA.16 five-moment repair into an Agent-3 receipt.

This module does *not* implement another PA.16 solver.  It consumes the exact
Agent-1 #953 selected ``T_sh`` certificate and the existing
``KokunoInnerJoinExit.solve_at_eta`` result, while independently checking that
its incoming five moments are the Agent-3 #949 current-Xi handoff and that

    pre_repair + repair_increment = exit_discrepancy.

The receipt is useful correction machinery/provenance, not a Navier--Stokes
residual claim.  In particular the current lane still lacks the completed
global Cartesian leading join, matched pressure and preregistered restricted
forcing, and the source analytic ``T_sh`` lower bound is not proved by #953.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from typing import Any, Mapping, Protocol

import numpy as np

from .kokuno_current_xi_moment_discrepancy_bridge import (
    CurrentXiMomentBridgeReceipt,
    SCHEMA as XI_BRIDGE_SCHEMA,
    TASK as XI_BRIDGE_TASK,
)


TASK = "KOKUNO-A3-CURRENT-PA16-REPAIR-RECEIPT-092"
SCHEMA = "kokuno-a3-current-pa16-repair-receipt-v1"
PARENT_AGENT3_PR = 949
PARENT_AGENT3_HEAD = "87c1ef765d6bf3613c7c6c301669ffccf5260117"
PARENT_AGENT3_SOURCE_PATH = (
    "src/openai_ns_reconstruction/kokuno_current_xi_moment_discrepancy_bridge.py"
)
PARENT_AGENT3_SOURCE_BLOB = "223058e2132bb22298c6a642b0f2a55562181e25"
UPSTREAM_AGENT1_PR = 953
UPSTREAM_AGENT1_HEAD = "bf16b88db92d23923b3e049c00e000a5bfa09992"
UPSTREAM_AGENT1_PARENT_PR = 947
UPSTREAM_AGENT1_PARENT_HEAD = "2797eed8e3b9b3374c676d8f517f05e3f0fc3e6b"
UPSTREAM_TSH_SOURCE_PATH = (
    "src/openai_ns_reconstruction/kokuno_pa16_current_tsh_certificate.py"
)
UPSTREAM_TSH_SOURCE_BLOB = "5fd610b417ccaa1c4673ea64985296b94dfe5946"
UPSTREAM_JOIN_SOURCE_PATH = (
    "src/openai_ns_reconstruction/kokuno_inner_join_exit.py"
)
UPSTREAM_JOIN_SOURCE_BLOB = "baa283c9592e19653e20f5c24371230d680704f5"
SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_PATH = "navier-stokes/navier_stokes_workbench.tex"
BOOKKEEPING_ABS_TOL = 5.0e-12
NONTRIVIAL_COEFFICIENT_FLOOR = 1.0e-14


class CurrentPA16RepairReceiptError(RuntimeError):
    """Raised when a current-lineage PA.16 repair handoff fails closed."""


@dataclass(frozen=True)
class UpstreamPA16Identity:
    pr_number: int
    exact_head: str
    parent_pr_number: int
    parent_exact_head: str
    tsh_source_path: str
    tsh_source_blob: str
    join_source_path: str
    join_source_blob: str


class CurrentPA16RepairBackend(Protocol):
    def current_xi_bridge_receipt(self) -> CurrentXiMomentBridgeReceipt: ...

    def upstream_identity(self) -> UpstreamPA16Identity: ...

    def tsh_report(self) -> Mapping[str, Any]: ...

    def pa16_solve_at_eta(self, eta: float) -> Mapping[str, Any]: ...


@dataclass(frozen=True)
class PA16RepairSample:
    eta: float
    incoming_scaled_discrepancy: tuple[float, float, float, float, float]
    pre_repair_scaled_discrepancy: tuple[float, float, float, float, float]
    repair_achieved_scaled_moments: tuple[float, float, float, float, float]
    exit_scaled_discrepancy: tuple[float, float, float, float, float]
    pa16_transformed_pre_discrepancy: tuple[float, float, float, float, float]
    pa16_transformed_exit_discrepancy: tuple[float, float, float, float, float]
    coefficient_vector: tuple[float, float, float, float, float]
    incoming_alignment_abs_max: float
    repair_algebra_abs_max: float
    target_alignment_abs_max: float
    solver_residual_alignment_abs_max: float
    pre_raw_l2_norm: float
    exit_raw_l2_norm: float
    pre_transformed_l2_norm: float
    exit_transformed_l2_norm: float
    raw_contraction_ratio: float
    transformed_contraction_ratio: float
    coefficient_l2_norm: float
    coefficient_max_abs: float
    coefficient_nontrivial: bool
    max_abs_exit_discrepancy: float
    scaled_jacobian_condition: float


@dataclass(frozen=True)
class CurrentPA16RepairReceipt:
    upstream: UpstreamPA16Identity
    xi_handoff_sha256: str
    xi_upstream_semantic_sha256: str
    tsh_semantic_sha256: str
    selected_route_ready: bool
    selected_T_sh: float
    separation_geometry_feasible: bool
    geometry_margin: float
    status: str
    blocked_reason: str | None
    samples: tuple[PA16RepairSample, ...]
    max_incoming_alignment_abs: float
    max_repair_algebra_abs: float
    max_exit_raw_l2_norm: float
    max_exit_transformed_l2_norm: float
    any_nontrivial_correction: bool
    receipt_sha256: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "upstream": asdict(self.upstream),
            "parent_agent3": {
                "pr_number": PARENT_AGENT3_PR,
                "exact_head": PARENT_AGENT3_HEAD,
                "source_path": PARENT_AGENT3_SOURCE_PATH,
                "source_blob": PARENT_AGENT3_SOURCE_BLOB,
            },
            "xi_handoff_sha256": self.xi_handoff_sha256,
            "xi_upstream_semantic_sha256": self.xi_upstream_semantic_sha256,
            "tsh_semantic_sha256": self.tsh_semantic_sha256,
            "selected_route_ready": self.selected_route_ready,
            "selected_T_sh": self.selected_T_sh,
            "separation_geometry_feasible": self.separation_geometry_feasible,
            "geometry_margin": self.geometry_margin,
            "status": self.status,
            "blocked_reason": self.blocked_reason,
            "samples": [asdict(sample) for sample in self.samples],
            "max_incoming_alignment_abs": self.max_incoming_alignment_abs,
            "max_repair_algebra_abs": self.max_repair_algebra_abs,
            "max_exit_raw_l2_norm": self.max_exit_raw_l2_norm,
            "max_exit_transformed_l2_norm": self.max_exit_transformed_l2_norm,
            "any_nontrivial_correction": self.any_nontrivial_correction,
            "receipt_sha256": self.receipt_sha256,
            "truth_boundary": truth_boundary(),
        }


def _canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256_json(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _finite_vector(value: Any, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape != (5,) or np.any(~np.isfinite(array)):
        raise CurrentPA16RepairReceiptError(f"{name} must be a finite five-vector")
    return array


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise CurrentPA16RepairReceiptError(f"{name} must be a mapping")
    return value


def _ratio(after: float, before: float) -> float:
    if before == 0.0:
        return 0.0 if after == 0.0 else math.inf
    return after / before


def _pa16_row_transform(vector: np.ndarray, eta: float) -> np.ndarray:
    M, I, J, S, C_p = (float(v) for v in vector)
    return np.asarray(
        [M, J - 4.0 * eta * I, I, S - 8.0 * eta * M, C_p], dtype=float
    )


def _require_backend(backend: Any) -> None:
    names = (
        "current_xi_bridge_receipt",
        "upstream_identity",
        "tsh_report",
        "pa16_solve_at_eta",
    )
    missing = [name for name in names if not callable(getattr(backend, name, None))]
    if missing:
        raise CurrentPA16RepairReceiptError(
            "backend is missing typed PA.16 methods: " + ", ".join(missing)
        )


def _require_upstream_identity(identity: Any) -> UpstreamPA16Identity:
    expected = UpstreamPA16Identity(
        pr_number=UPSTREAM_AGENT1_PR,
        exact_head=UPSTREAM_AGENT1_HEAD,
        parent_pr_number=UPSTREAM_AGENT1_PARENT_PR,
        parent_exact_head=UPSTREAM_AGENT1_PARENT_HEAD,
        tsh_source_path=UPSTREAM_TSH_SOURCE_PATH,
        tsh_source_blob=UPSTREAM_TSH_SOURCE_BLOB,
        join_source_path=UPSTREAM_JOIN_SOURCE_PATH,
        join_source_blob=UPSTREAM_JOIN_SOURCE_BLOB,
    )
    if not isinstance(identity, UpstreamPA16Identity) or identity != expected:
        raise CurrentPA16RepairReceiptError("repair handoff is not pinned Agent-1 #953")
    return identity


def _validate_bridge(receipt: Any) -> CurrentXiMomentBridgeReceipt:
    if not isinstance(receipt, CurrentXiMomentBridgeReceipt):
        raise CurrentPA16RepairReceiptError(
            "current_xi_bridge_receipt must be a #949 CurrentXiMomentBridgeReceipt"
        )
    payload = receipt.to_dict()
    if payload.get("schema") != XI_BRIDGE_SCHEMA or payload.get("task") != XI_BRIDGE_TASK:
        raise CurrentPA16RepairReceiptError("current-Xi bridge schema/task drifted")
    truth = _mapping(payload.get("truth_boundary"), "current-Xi truth boundary")
    if truth.get("discrepancy_from_complete_ns_defect") is not False:
        raise CurrentPA16RepairReceiptError("current-Xi bridge attempted full-NS promotion")
    if truth.get("authorized_for_gain_gated_ns_stage") is not False:
        raise CurrentPA16RepairReceiptError("current-Xi bridge attempted correction authorization")
    return receipt


def _validate_tsh_report(
    report: Mapping[str, Any], bridge: CurrentXiMomentBridgeReceipt
) -> tuple[str, bool, float, bool, float]:
    source = _mapping(report.get("source"), "T_sh report source")
    if source.get("repository") != SOURCE_READER_REPO:
        raise CurrentPA16RepairReceiptError("T_sh source repository drifted")
    if source.get("commit") != SOURCE_READER_HEAD or source.get("path") != SOURCE_READER_PATH:
        raise CurrentPA16RepairReceiptError("T_sh corrected-reader provenance drifted")
    if report.get("current_parent_semantic_sha256") != bridge.upstream_semantic_sha256:
        raise CurrentPA16RepairReceiptError("T_sh parent semantic SHA does not match #949 input")
    semantic_sha = report.get("semantic_sha256")
    if not isinstance(semantic_sha, str) or len(semantic_sha) != 64:
        raise CurrentPA16RepairReceiptError("T_sh semantic SHA is missing")
    truth = _mapping(report.get("truth_boundary"), "T_sh truth boundary")
    required_true = (
        "current_lineage_Xi_data_bound",
        "PA10_T_sh_formula_executable",
        "selected_T_sh_lower_bound_numerically_instantiated",
    )
    for key in required_true:
        if truth.get(key) is not True:
            raise CurrentPA16RepairReceiptError(f"T_sh report lost required flag: {key}")
    if truth.get("source_T_sh_lower_bound_verified") is not False:
        raise CurrentPA16RepairReceiptError("finite-sample T_sh evidence was promoted to source theorem")
    if truth.get("five_moment_repair_applied") is not False:
        raise CurrentPA16RepairReceiptError("T_sh report may not pre-claim PA.16 repair")

    geometry = _mapping(report.get("geometry"), "T_sh geometry")
    selected_T_sh = float(geometry.get("selected_T_sh"))
    margin = float(geometry.get("geometry_margin"))
    feasible = geometry.get("separation_geometry_feasible")
    ready = report.get("selected_route_ready")
    if not math.isfinite(selected_T_sh) or selected_T_sh <= 0.0 or not math.isfinite(margin):
        raise CurrentPA16RepairReceiptError("T_sh geometry contains nonfinite values")
    if not isinstance(feasible, bool) or not isinstance(ready, bool):
        raise CurrentPA16RepairReceiptError("T_sh route flags must be boolean")
    if ready and not feasible:
        raise CurrentPA16RepairReceiptError("T_sh report marks an infeasible geometry ready")
    return semantic_sha, ready, selected_T_sh, feasible, margin


def _sample_from_solve(
    bridge_sample: Any, solve: Mapping[str, Any]
) -> PA16RepairSample:
    eta = float(bridge_sample.eta)
    if not math.isclose(float(solve.get("eta")), eta, rel_tol=0.0, abs_tol=0.0):
        raise CurrentPA16RepairReceiptError("PA.16 solve eta does not match current-Xi sample")
    incoming = _finite_vector(solve.get("incoming_scaled_discrepancy"), "incoming")
    expected_incoming = np.asarray(bridge_sample.incoming_scaled_discrepancy, dtype=float)
    incoming_alignment = float(np.max(np.abs(incoming - expected_incoming)))
    if incoming_alignment > BOOKKEEPING_ABS_TOL:
        raise CurrentPA16RepairReceiptError("PA.16 incoming discrepancy is not the #949 handoff")

    pre = _finite_vector(solve.get("pre_repair_scaled_discrepancy"), "pre_repair")
    achieved = _finite_vector(solve.get("repair_achieved_scaled_moments"), "repair_achieved")
    exit_error = _finite_vector(solve.get("exit_scaled_discrepancy"), "exit")
    target = _finite_vector(solve.get("repair_target_scaled_moments"), "repair_target")
    solver_residual = _finite_vector(
        solve.get("repair_residual_scaled_moments"), "repair_residual"
    )
    coefficients = _finite_vector(solve.get("coefficients"), "coefficients")

    algebra_error = float(np.max(np.abs(pre + achieved - exit_error)))
    if algebra_error > BOOKKEEPING_ABS_TOL:
        raise CurrentPA16RepairReceiptError("pre + repair increment does not reproduce exit")
    target_error = float(np.max(np.abs(target + pre)))
    if target_error > BOOKKEEPING_ABS_TOL:
        raise CurrentPA16RepairReceiptError("PA.16 repair target is not -pre_repair")
    residual_error = float(np.max(np.abs((achieved - target) - solver_residual)))
    if residual_error > BOOKKEEPING_ABS_TOL:
        raise CurrentPA16RepairReceiptError("PA.16 solver residual receipt does not replay")

    reported_max = float(solve.get("max_abs_exit_discrepancy"))
    condition = float(solve.get("scaled_jacobian_condition"))
    if not math.isfinite(reported_max) or reported_max < 0.0:
        raise CurrentPA16RepairReceiptError("invalid reported exit discrepancy")
    if abs(reported_max - float(np.max(np.abs(exit_error)))) > BOOKKEEPING_ABS_TOL:
        raise CurrentPA16RepairReceiptError("reported max exit discrepancy does not replay")
    if not math.isfinite(condition) or condition <= 0.0:
        raise CurrentPA16RepairReceiptError("scaled Jacobian condition must be positive finite")

    pre_t = _pa16_row_transform(pre, eta)
    exit_t = _pa16_row_transform(exit_error, eta)
    pre_norm = float(np.linalg.norm(pre))
    exit_norm = float(np.linalg.norm(exit_error))
    pre_t_norm = float(np.linalg.norm(pre_t))
    exit_t_norm = float(np.linalg.norm(exit_t))
    coefficient_l2 = float(np.linalg.norm(coefficients))
    coefficient_max = float(np.max(np.abs(coefficients)))

    return PA16RepairSample(
        eta=eta,
        incoming_scaled_discrepancy=tuple(float(v) for v in incoming),
        pre_repair_scaled_discrepancy=tuple(float(v) for v in pre),
        repair_achieved_scaled_moments=tuple(float(v) for v in achieved),
        exit_scaled_discrepancy=tuple(float(v) for v in exit_error),
        pa16_transformed_pre_discrepancy=tuple(float(v) for v in pre_t),
        pa16_transformed_exit_discrepancy=tuple(float(v) for v in exit_t),
        coefficient_vector=tuple(float(v) for v in coefficients),
        incoming_alignment_abs_max=incoming_alignment,
        repair_algebra_abs_max=algebra_error,
        target_alignment_abs_max=target_error,
        solver_residual_alignment_abs_max=residual_error,
        pre_raw_l2_norm=pre_norm,
        exit_raw_l2_norm=exit_norm,
        pre_transformed_l2_norm=pre_t_norm,
        exit_transformed_l2_norm=exit_t_norm,
        raw_contraction_ratio=_ratio(exit_norm, pre_norm),
        transformed_contraction_ratio=_ratio(exit_t_norm, pre_t_norm),
        coefficient_l2_norm=coefficient_l2,
        coefficient_max_abs=coefficient_max,
        coefficient_nontrivial=coefficient_l2 > NONTRIVIAL_COEFFICIENT_FLOOR,
        max_abs_exit_discrepancy=reported_max,
        scaled_jacobian_condition=condition,
    )


def materialize_current_pa16_repair_receipt(
    backend: CurrentPA16RepairBackend,
) -> CurrentPA16RepairReceipt:
    """Verify the exact current PA.16 repair without exposing scientific knobs."""
    _require_backend(backend)
    identity = _require_upstream_identity(backend.upstream_identity())
    bridge = _validate_bridge(backend.current_xi_bridge_receipt())
    tsh_report = _mapping(backend.tsh_report(), "T_sh report")
    tsh_sha, ready, selected_T_sh, feasible, margin = _validate_tsh_report(tsh_report, bridge)

    if not ready:
        samples: tuple[PA16RepairSample, ...] = ()
        status = "blocked_by_current_tsh_geometry"
        blocked_reason = "Agent-1 #953 selected numerical T_sh route is not executable under current geometry"
    else:
        built = [
            _sample_from_solve(sample, _mapping(backend.pa16_solve_at_eta(sample.eta), "PA.16 solve"))
            for sample in bridge.samples
        ]
        samples = tuple(built)
        status = "current_source_moment_repair_verified"
        blocked_reason = None

    max_incoming = max((s.incoming_alignment_abs_max for s in samples), default=0.0)
    max_algebra = max((s.repair_algebra_abs_max for s in samples), default=0.0)
    max_exit = max((s.exit_raw_l2_norm for s in samples), default=0.0)
    max_exit_t = max((s.exit_transformed_l2_norm for s in samples), default=0.0)
    any_nontrivial = any(s.coefficient_nontrivial for s in samples)
    unsigned = {
        "schema": SCHEMA,
        "task": TASK,
        "upstream": asdict(identity),
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_agent3_source_blob": PARENT_AGENT3_SOURCE_BLOB,
        "xi_handoff_sha256": bridge.handoff_sha256,
        "xi_upstream_semantic_sha256": bridge.upstream_semantic_sha256,
        "tsh_semantic_sha256": tsh_sha,
        "selected_route_ready": ready,
        "selected_T_sh": selected_T_sh,
        "separation_geometry_feasible": feasible,
        "geometry_margin": margin,
        "status": status,
        "blocked_reason": blocked_reason,
        "samples": [asdict(sample) for sample in samples],
    }
    receipt_sha = _sha256_json(unsigned)
    return CurrentPA16RepairReceipt(
        upstream=identity,
        xi_handoff_sha256=bridge.handoff_sha256,
        xi_upstream_semantic_sha256=bridge.upstream_semantic_sha256,
        tsh_semantic_sha256=tsh_sha,
        selected_route_ready=ready,
        selected_T_sh=selected_T_sh,
        separation_geometry_feasible=feasible,
        geometry_margin=margin,
        status=status,
        blocked_reason=blocked_reason,
        samples=samples,
        max_incoming_alignment_abs=max_incoming,
        max_repair_algebra_abs=max_algebra,
        max_exit_raw_l2_norm=max_exit,
        max_exit_transformed_l2_norm=max_exit_t,
        any_nontrivial_correction=any_nontrivial,
        receipt_sha256=receipt_sha,
    )


def truth_boundary() -> dict[str, bool]:
    return {
        "public_reconstruction_source": True,
        "current_lineage_xi_discrepancy_bound": True,
        "current_numerical_tsh_certificate_bound": True,
        "existing_agent1_pa16_solver_reused": True,
        "current_pa16_repair_receipt_materialized_when_route_ready": True,
        "source_T_sh_lower_bound_verified": False,
        "global_inner_to_outer_join_completed": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "discrepancy_from_complete_ns_defect": False,
        "authorized_for_gain_gated_ns_stage": False,
        "current_real_ns_correction_velocity_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
