"""Bind the current PA.16 repair receipt to the actual joined-profile realization.

This Agent-3 increment closes one provenance/application seam only.  The parent
CR002 #958 path first checksum-verifies the Agent-3 #949 current-Xi handoff and
then materializes the Agent-3 #957 PA.16 repair receipt.  Agent-1 #959 now
materializes the current source-coordinate joined profile through ``X_h`` and
states that the existing PA.16 coefficients are consumed in the repair annulus.

The present binder checks that both artifacts use the same current #953
certificate identity and, on frozen interior repair-annulus probes, compares the
actual #959 joined profile with an explicit-coefficient replay using the exact
five-vector recorded by #957.  It also checks that replacing that vector by
zero changes a nontrivial repair realization, preventing an implementation that
silently ignores the repair coefficients from passing.

This is source-coordinate correction-application evidence.  It is not a
Cartesian correction velocity, not a complete Navier--Stokes defect, and not
held-out residual evidence.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from typing import Any, Mapping, Protocol

import numpy as np

from .kokuno_current_pa16_bridge_integrity import (
    materialize_checksum_verified_current_pa16_repair_receipt,
)
from .kokuno_current_pa16_repair_receipt import (
    CurrentPA16RepairBackend,
    CurrentPA16RepairReceipt,
)


TASK = "KOKUNO-A3-CURRENT-PA16-PROFILE-APPLICATION-094"
SCHEMA = "kokuno-a3-current-pa16-profile-application-v1"
PARENT_CR002_PR = 958
PARENT_CR002_HEAD = "0c16e1f1a2ff22b2a3950e77041ce1fbb11e19c7"
PARENT_CR002_SOURCE_PATH = (
    "src/openai_ns_reconstruction/kokuno_current_pa16_bridge_integrity.py"
)
PARENT_CR002_SOURCE_BLOB = "a58349dbe7fa510573349ae3dd26735b6aedde94"
UPSTREAM_AGENT1_PR = 959
UPSTREAM_AGENT1_HEAD = "61b8730fba5623c5acfbd3ec42d54f7cc4f45748"
UPSTREAM_AGENT1_SOURCE_PATH = (
    "src/openai_ns_reconstruction/kokuno_pa16_current_joined_profile.py"
)
UPSTREAM_AGENT1_SOURCE_BLOB = "e38b9b79378aaf97a17226e41da75b477f26ae23"
UPSTREAM_AGENT1_SCHEMA = "kokuno-pa16-current-joined-profile-v1"
SOURCE_READER_REPO = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_READER_PATH = "navier-stokes/navier_stokes_workbench.tex"
REPAIR_LOG_X_PROBES = (-5.90, -5.65, -5.35, -5.10)
APPLICATION_ABS_TOL = 5.0e-10
APPLICATION_REL_TOL = 5.0e-10
NONTRIVIAL_EFFECT_REL_FLOOR = 1.0e-15


class CurrentPA16ProfileApplicationError(RuntimeError):
    """Raised when the repair receipt and joined profile fail to bind."""


@dataclass(frozen=True)
class JoinedProfileIdentity:
    pr_number: int
    exact_head: str
    source_path: str
    source_blob: str
    report_schema: str


class CurrentPA16ProfileApplicationBackend(CurrentPA16RepairBackend, Protocol):
    def joined_profile_identity(self) -> JoinedProfileIdentity: ...

    def joined_profile_report(self) -> Mapping[str, Any]: ...

    def joined_profile_values(self, eta: float, X: tuple[float, ...]) -> Mapping[str, Any]: ...

    def joined_profile_replay_with_coefficients(
        self,
        eta: float,
        X: tuple[float, ...],
        coefficients: tuple[float, float, float, float, float],
    ) -> Mapping[str, Any]: ...


@dataclass(frozen=True)
class PA16ProfileApplicationSample:
    eta: float
    log_x_probes: tuple[float, ...]
    X_probes: tuple[float, ...]
    coefficient_vector: tuple[float, float, float, float, float]
    coefficient_sha256: str
    application_abs_max: float
    application_rel_max: float
    zero_coefficient_effect_l2: float
    zero_coefficient_effect_max: float
    zero_coefficient_effect_relative: float
    nontrivial_repair_effect_observed: bool


@dataclass(frozen=True)
class CurrentPA16ProfileApplicationReceipt:
    upstream: JoinedProfileIdentity
    parent_repair_receipt_sha256: str
    parent_tsh_semantic_sha256: str
    joined_profile_semantic_sha256: str
    selected_route_ready: bool
    selected_T_sh: float
    log_X_R: float
    log_X_h: float
    status: str
    blocked_reason: str | None
    samples: tuple[PA16ProfileApplicationSample, ...]
    max_application_abs: float
    max_application_rel: float
    max_zero_coefficient_effect_l2: float
    any_nontrivial_repair_effect: bool
    receipt_sha256: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "upstream": asdict(self.upstream),
            "parent_cr002": {
                "pr_number": PARENT_CR002_PR,
                "exact_head": PARENT_CR002_HEAD,
                "source_path": PARENT_CR002_SOURCE_PATH,
                "source_blob": PARENT_CR002_SOURCE_BLOB,
            },
            "parent_repair_receipt_sha256": self.parent_repair_receipt_sha256,
            "parent_tsh_semantic_sha256": self.parent_tsh_semantic_sha256,
            "joined_profile_semantic_sha256": self.joined_profile_semantic_sha256,
            "selected_route_ready": self.selected_route_ready,
            "selected_T_sh": self.selected_T_sh,
            "log_X_R": self.log_X_R,
            "log_X_h": self.log_X_h,
            "status": self.status,
            "blocked_reason": self.blocked_reason,
            "samples": [asdict(sample) for sample in self.samples],
            "max_application_abs": self.max_application_abs,
            "max_application_rel": self.max_application_rel,
            "max_zero_coefficient_effect_l2": self.max_zero_coefficient_effect_l2,
            "any_nontrivial_repair_effect": self.any_nontrivial_repair_effect,
            "receipt_sha256": self.receipt_sha256,
            "truth_boundary": truth_boundary(),
        }


def _canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256_json(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise CurrentPA16ProfileApplicationError(f"{name} must be a mapping")
    return value


def _require_backend(backend: Any) -> None:
    names = (
        "joined_profile_identity",
        "joined_profile_report",
        "joined_profile_values",
        "joined_profile_replay_with_coefficients",
    )
    missing = [name for name in names if not callable(getattr(backend, name, None))]
    if missing:
        raise CurrentPA16ProfileApplicationError(
            "backend is missing typed joined-profile methods: " + ", ".join(missing)
        )


def _require_identity(identity: Any) -> JoinedProfileIdentity:
    expected = JoinedProfileIdentity(
        pr_number=UPSTREAM_AGENT1_PR,
        exact_head=UPSTREAM_AGENT1_HEAD,
        source_path=UPSTREAM_AGENT1_SOURCE_PATH,
        source_blob=UPSTREAM_AGENT1_SOURCE_BLOB,
        report_schema=UPSTREAM_AGENT1_SCHEMA,
    )
    if not isinstance(identity, JoinedProfileIdentity) or identity != expected:
        raise CurrentPA16ProfileApplicationError(
            "joined profile is not the pinned Agent-1 #959 artifact"
        )
    return identity


def _finite_profile_matrix(value: Any, name: str, n: int) -> np.ndarray:
    report = _mapping(value, name)
    rows = []
    for key in ("F", "U", "E"):
        array = np.asarray(report.get(key), dtype=float)
        if array.shape != (n,) or np.any(~np.isfinite(array)):
            raise CurrentPA16ProfileApplicationError(
                f"{name}.{key} must be a finite vector of length {n}"
            )
        rows.append(array)
    return np.stack(rows, axis=-1)


def _validate_profile_report(
    report: Mapping[str, Any], repair: CurrentPA16RepairReceipt
) -> tuple[str, bool, float, float, float]:
    source = _mapping(report.get("source"), "joined-profile source")
    if source.get("repository") != SOURCE_READER_REPO:
        raise CurrentPA16ProfileApplicationError("joined-profile source repository drifted")
    if source.get("commit") != SOURCE_READER_HEAD or source.get("path") != SOURCE_READER_PATH:
        raise CurrentPA16ProfileApplicationError("joined-profile corrected-reader provenance drifted")

    semantic_sha = report.get("semantic_sha256")
    if not isinstance(semantic_sha, str) or len(semantic_sha) != 64:
        raise CurrentPA16ProfileApplicationError("joined-profile semantic SHA is missing")
    if report.get("parent_certificate_semantic_sha256") != repair.tsh_semantic_sha256:
        raise CurrentPA16ProfileApplicationError(
            "joined profile does not use the same #953 T_sh certificate as the repair receipt"
        )

    geometry = _mapping(report.get("geometry"), "joined-profile geometry")
    route_ready = geometry.get("route_ready")
    selected_T_sh = float(geometry.get("selected_T_sh"))
    log_X_R = float(geometry.get("log_X_R"))
    log_X_h = float(geometry.get("log_X_h"))
    if not isinstance(route_ready, bool):
        raise CurrentPA16ProfileApplicationError("joined-profile route_ready must be boolean")
    if route_ready != repair.selected_route_ready:
        raise CurrentPA16ProfileApplicationError("joined-profile route readiness disagrees with #957")
    if not all(math.isfinite(v) for v in (selected_T_sh, log_X_R, log_X_h)):
        raise CurrentPA16ProfileApplicationError("joined-profile geometry must be finite")
    if selected_T_sh != repair.selected_T_sh:
        raise CurrentPA16ProfileApplicationError("joined-profile selected T_sh disagrees with #957")
    if abs((log_X_R - 5.0) - log_X_h) > 5.0e-13:
        raise CurrentPA16ProfileApplicationError("joined-profile X_h geometry no longer uses log x=-5")

    truth = _mapping(report.get("truth_boundary"), "joined-profile truth boundary")
    required_true = (
        "current_lineage_prefix_0_to_Xi_bound",
        "current_candidate_side_join_Xi_to_Xh_materialized",
        "current_candidate_side_PA16_coefficients_consumed",
        "physical_X_F_U_E_profile_executable_through_Xh",
    )
    required_false = (
        "source_T_sh_lower_bound_verified",
        "source_global_inner_to_outer_join_admitted",
        "outer_global_leading_velocity_materialized",
        "cartesian_velocity_materialized",
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
            raise CurrentPA16ProfileApplicationError(
                f"joined-profile truth boundary lost required flag: {key}"
            )
    for key in required_false:
        if truth.get(key) is not False:
            raise CurrentPA16ProfileApplicationError(
                f"joined-profile truth boundary attempted promotion: {key}"
            )
    return semantic_sha, route_ready, selected_T_sh, log_X_R, log_X_h


def _sample_application(
    backend: CurrentPA16ProfileApplicationBackend,
    *,
    eta: float,
    coefficient_vector: tuple[float, float, float, float, float],
    coefficient_nontrivial: bool,
    log_X_R: float,
) -> PA16ProfileApplicationSample:
    X = tuple(math.exp(log_X_R + log_x) for log_x in REPAIR_LOG_X_PROBES)
    if any(not math.isfinite(value) or value <= 0.0 for value in X):
        raise CurrentPA16ProfileApplicationError("repair-annulus X probes are nonfinite")

    actual = _finite_profile_matrix(
        backend.joined_profile_values(float(eta), X), "joined_profile_values", len(X)
    )
    replay = _finite_profile_matrix(
        backend.joined_profile_replay_with_coefficients(
            float(eta), X, coefficient_vector
        ),
        "joined_profile_replay_with_coefficients",
        len(X),
    )
    zero = _finite_profile_matrix(
        backend.joined_profile_replay_with_coefficients(
            float(eta), X, (0.0, 0.0, 0.0, 0.0, 0.0)
        ),
        "joined_profile_zero_coefficient_replay",
        len(X),
    )

    delta = actual - replay
    abs_max = float(np.max(np.abs(delta)))
    scale = np.maximum(1.0, np.maximum(np.abs(actual), np.abs(replay)))
    rel_max = float(np.max(np.abs(delta) / scale))
    if abs_max > APPLICATION_ABS_TOL and rel_max > APPLICATION_REL_TOL:
        raise CurrentPA16ProfileApplicationError(
            "joined profile does not reproduce the exact #957 PA.16 coefficient application"
        )

    zero_effect = actual - zero
    effect_l2 = float(np.linalg.norm(zero_effect))
    effect_max = float(np.max(np.abs(zero_effect)))
    actual_scale = max(1.0, float(np.linalg.norm(actual)))
    effect_relative = effect_l2 / actual_scale
    observed = effect_relative > NONTRIVIAL_EFFECT_REL_FLOOR
    if coefficient_nontrivial and not observed:
        raise CurrentPA16ProfileApplicationError(
            "nontrivial #957 coefficient vector has no observable effect in the #959 repair annulus"
        )

    coefficient_sha = _sha256_json(list(coefficient_vector))
    return PA16ProfileApplicationSample(
        eta=float(eta),
        log_x_probes=tuple(float(v) for v in REPAIR_LOG_X_PROBES),
        X_probes=tuple(float(v) for v in X),
        coefficient_vector=coefficient_vector,
        coefficient_sha256=coefficient_sha,
        application_abs_max=abs_max,
        application_rel_max=rel_max,
        zero_coefficient_effect_l2=effect_l2,
        zero_coefficient_effect_max=effect_max,
        zero_coefficient_effect_relative=effect_relative,
        nontrivial_repair_effect_observed=observed,
    )


def materialize_current_pa16_profile_application_receipt(
    backend: CurrentPA16ProfileApplicationBackend,
) -> CurrentPA16ProfileApplicationReceipt:
    """Bind checksum-verified #957 repair coefficients to the #959 profile.

    The public call accepts only one typed backend.  It exposes no residual,
    defect, discrepancy, coefficient, T_sh, pressure, forcing, damping,
    viscosity, repair-probe, tolerance, stage-budget, or scientific-threshold
    argument.
    """
    _require_backend(backend)
    identity = _require_identity(backend.joined_profile_identity())
    repair = materialize_checksum_verified_current_pa16_repair_receipt(backend)
    report = _mapping(backend.joined_profile_report(), "joined_profile_report")
    semantic_sha, route_ready, selected_T_sh, log_X_R, log_X_h = _validate_profile_report(
        report, repair
    )

    samples: list[PA16ProfileApplicationSample] = []
    if route_ready:
        if repair.status != "current_source_moment_repair_verified":
            raise CurrentPA16ProfileApplicationError(
                "ready joined profile requires a verified #957 repair receipt"
            )
        for repair_sample in repair.samples:
            samples.append(
                _sample_application(
                    backend,
                    eta=float(repair_sample.eta),
                    coefficient_vector=tuple(
                        float(v) for v in repair_sample.coefficient_vector
                    ),
                    coefficient_nontrivial=bool(repair_sample.coefficient_nontrivial),
                    log_X_R=log_X_R,
                )
            )
        if not samples:
            raise CurrentPA16ProfileApplicationError(
                "ready joined profile has no PA.16 repair samples"
            )
        status = "current_source_moment_repair_applied_in_joined_profile"
        blocked_reason = None
    else:
        if repair.status != "blocked_by_current_tsh_geometry" or repair.samples:
            raise CurrentPA16ProfileApplicationError(
                "blocked joined profile disagrees with the fail-closed #957 repair receipt"
            )
        status = "blocked_by_current_tsh_geometry"
        blocked_reason = "selected numerical T_sh route is not executable"

    max_abs = max((sample.application_abs_max for sample in samples), default=0.0)
    max_rel = max((sample.application_rel_max for sample in samples), default=0.0)
    max_effect = max((sample.zero_coefficient_effect_l2 for sample in samples), default=0.0)
    any_effect = any(sample.nontrivial_repair_effect_observed for sample in samples)
    if route_ready and repair.any_nontrivial_correction and not any_effect:
        raise CurrentPA16ProfileApplicationError(
            "repair receipt is nontrivial but the joined profile shows no repair effect"
        )

    unsigned = {
        "schema": SCHEMA,
        "upstream": asdict(identity),
        "parent_repair_receipt_sha256": repair.receipt_sha256,
        "parent_tsh_semantic_sha256": repair.tsh_semantic_sha256,
        "joined_profile_semantic_sha256": semantic_sha,
        "selected_route_ready": route_ready,
        "selected_T_sh": selected_T_sh,
        "log_X_R": log_X_R,
        "log_X_h": log_X_h,
        "status": status,
        "blocked_reason": blocked_reason,
        "samples": [asdict(sample) for sample in samples],
        "max_application_abs": max_abs,
        "max_application_rel": max_rel,
        "max_zero_coefficient_effect_l2": max_effect,
        "any_nontrivial_repair_effect": any_effect,
    }
    return CurrentPA16ProfileApplicationReceipt(
        upstream=identity,
        parent_repair_receipt_sha256=repair.receipt_sha256,
        parent_tsh_semantic_sha256=repair.tsh_semantic_sha256,
        joined_profile_semantic_sha256=semantic_sha,
        selected_route_ready=route_ready,
        selected_T_sh=selected_T_sh,
        log_X_R=log_X_R,
        log_X_h=log_X_h,
        status=status,
        blocked_reason=blocked_reason,
        samples=tuple(samples),
        max_application_abs=max_abs,
        max_application_rel=max_rel,
        max_zero_coefficient_effect_l2=max_effect,
        any_nontrivial_repair_effect=any_effect,
        receipt_sha256=_sha256_json(unsigned),
    )


def truth_boundary() -> dict[str, bool]:
    return {
        "checksum_verified_current_xi_lineage_consumed": True,
        "current_pa16_repair_receipt_bound_to_joined_profile": True,
        "exact_repair_coefficients_replayed_in_repair_annulus": True,
        "zero_coefficient_negative_control_required_for_nontrivial_repair": True,
        "source_coordinate_repair_application_verified": True,
        "source_T_sh_lower_bound_verified": False,
        "source_global_inner_to_outer_join_admitted": False,
        "global_cartesian_leading_velocity_materialized": False,
        "current_real_ns_correction_velocity_materialized": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "discrepancy_from_complete_ns_defect": False,
        "authorized_for_gain_gated_ns_stage": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
