"""A2 identity-bound full log-X main-pulse leading + frozen complete-curl oscillation.

Consumes exact Agent-1 PR #1107 without altering its leading field and composes it
with the already-frozen Agent-2 complete-curl oscillatory runtime.  Agent-1 #1107
extends the current principal main-pulse realization through the public source
coordinate xi=11 by evaluating the large-X tail in log X.  Its principal pulse
amplitude remains a repository-autonomous approximation, not the corrected-source
exact Amp(eta) root.  This module therefore claims only project-domain additive
composition, not paper exactness, terminal/global completion, pressure/forcing, or
Navier--Stokes validation.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
import importlib
import inspect
import json
import math
import tempfile
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from . import kokuno_current_partial_composite_identity_save_load as _id
from . import kokuno_current_partial_leading_oscillatory_velocity as _comp
from . import kokuno_current_main_pulse_leading_oscillatory_identity as _parent

TASK = "K2-OSC-101"
SCHEMA = "kokuno-a2-current-main-pulse-logx-leading-oscillatory-v1"
SAVED_SCHEMA = SCHEMA + "-saved"
RECEIPT_SCHEMA = SCHEMA + "-receipt"

PARENT_A2_PR = 1108
PARENT_A2_HEAD = "9f221dd57ef4b5e2d2e80e24c8cf531e991a5387"
PARENT_A2_BLOB = "d8bde2336dc4a8f8dd7cf0ac6fecd184ea693390"

A1_PR = 1107
A1_HEAD = "45da043dd2b4cd067f005a72c7e21fd0f2bcf309"
A1_MODULE = "openai_ns_reconstruction.kokuno_pa16_current_cartesian_main_pulse_logx"
A1_CLASS = "KokunoPA16CurrentCartesianMainPulseLogX"
A1_BLOB = "477bda899694e4c75786d59850a2eddada669f07"
SOURCE_MAIN_XI_END = 11.0

SOURCE = {
    "repository": "KokunoYumeto/yang-mills-interacting-workbench",
    "commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
    "path": "navier-stokes/navier_stokes_workbench.tex",
    "blob_sha1": "205a99807302e21a51c5eaf223390c0dfc42bcd0",
    "corrected_release_date": "2026-09-09",
    "classification": "structural_math_provenance_only",
    "paper_exact_claim": False,
}

COMP_ATOL = 1.0e-13
REPLAY_ATOL = 2.0e-12
SIGNAL_FLOOR = 1.0e-12

INNER = np.asarray(
    [
        (0.32, 0.11, -0.20),
        (0.41, -0.17, -0.08),
        (-0.36, 0.24, 0.05),
        (-0.52, -0.16, 0.16),
        (0.58, 0.21, 0.22),
        (0.47, -0.31, -0.14),
    ],
    dtype=float,
)
INNER_T = np.asarray((0.31, 0.39, 0.47, 0.55, 0.63, 0.71), dtype=float)
AXIS = np.asarray(((0.0, 0.0, -0.10), (0.0, 0.0, 0.0), (0.0, 0.0, 0.12)), dtype=float)
AXIS_T = np.asarray((0.37, 0.51, 0.67), dtype=float)
OUTSIDE = np.asarray(((3.5, 0.0, 0.0), (0.0, -3.5, 0.4), (2.8, 2.8, -0.5)), dtype=float)
OUTSIDE_T = np.asarray((0.41, 0.53, 0.69), dtype=float)
TAIL_FRACTIONS = np.asarray((0.18, 0.55, 0.88), dtype=float)
TAIL_T = np.asarray((0.43, 0.55, 0.67), dtype=float)
TAIL_THETA = np.asarray((0.31, 0.97, 1.63), dtype=float)

A1_TRUTH = {
    "overflow_safe_logX_similarity_materialized": True,
    "overflow_safe_logF_cartesian_swirl_materialized": True,
    "full_source_xi_11_current_cartesian_materialized": True,
    "source_exact_amplitude_root_materialized": False,
    "source_pulse_end_MJ_corrections_materialized": False,
    "source_hidden_parameters_recovered": False,
    "source_terminal_tail_schedule_bound_into_current_velocity": False,
    "source_exterior_heat_replacement_materialized": False,
    "outer_global_leading_velocity_materialized": False,
    "unified_global_cartesian_velocity_export_ready": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "heldout_ns_residual_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "pde_validated": False,
    "paper_exact": False,
}


def _truth() -> dict[str, bool]:
    return {
        "current_logX_main_pulse_leading_consumed": True,
        "full_source_xi_11_current_cartesian_leading_materialized": True,
        "current_logX_main_pulse_leading_plus_frozen_complete_curl_oscillation_materialized": True,
        "frozen_complete_curl_oscillation_composed": True,
        "full_concrete_oscillatory_runtime_digest_bound": True,
        "identity_preserving_save_load_available": True,
        "field_summands_changed_by_this_increment": False,
        "source_exact_main_pulse_amplitude_materialized": False,
        "source_pulse_end_MJ_corrections_materialized": False,
        "terminal_global_leading_velocity_materialized": False,
        "agent3_correction_velocity_composed": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "complete_ns_residual_assessed": False,
        "same_protocol_comparable_to_st006": False,
        "residual_reduction_claimed": False,
        "velocity_export_ready": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def _parent_blob() -> str:
    path = getattr(_parent, "__file__", None)
    if not path or _comp._git_blob_sha1(path) != PARENT_A2_BLOB:
        raise RuntimeError("A2 #1108 parent source drifted")
    return PARENT_A2_BLOB


def _a1_identity(backend: Any) -> dict[str, Any]:
    cls = type(backend)
    if cls.__module__ != A1_MODULE or cls.__name__ != A1_CLASS:
        raise RuntimeError("loaded leading backend is not exact A1 #1107 class")
    path = inspect.getsourcefile(cls)
    if not path or _comp._git_blob_sha1(path) != A1_BLOB:
        raise RuntimeError("A1 #1107 source blob drifted")
    truth = getattr(backend, "truth_boundary", None)
    if not isinstance(truth, Mapping):
        raise RuntimeError("A1 #1107 truth boundary missing")
    for key, expected in A1_TRUTH.items():
        if truth.get(key) is not expected:
            raise RuntimeError(f"A1 #1107 truth drifted at {key}")
    cfg = backend.configuration()
    semantic = str(getattr(backend, "semantic_sha256", ""))
    if not isinstance(cfg, Mapping) or len(semantic) != 64:
        raise RuntimeError("A1 #1107 identity malformed")
    finite_xi = float(getattr(backend.parent, "xi_materializable_max", np.nan))
    log_x_p = float(getattr(backend, "log_X_p", np.nan))
    log_x_end = float(getattr(backend, "log_X_source_end", np.nan))
    lam = float(getattr(backend, "lambda_value", np.nan))
    if not (
        np.isfinite(finite_xi)
        and 0.0 < finite_xi < SOURCE_MAIN_XI_END
        and np.isfinite(log_x_p)
        and np.isfinite(log_x_end)
        and log_x_end > log_x_p
        and np.isfinite(lam)
        and 0.0 < lam < 1.0
    ):
        raise RuntimeError("A1 #1107 log-X domain identity invalid")
    endpoint = float(backend.log_X_from_xi(SOURCE_MAIN_XI_END))
    tol = 256.0 * np.finfo(float).eps * max(1.0, abs(log_x_end))
    if abs(endpoint - log_x_end) > tol:
        raise RuntimeError("A1 #1107 xi=11 endpoint identity drifted")
    return {
        "semantic_sha256": semantic,
        "configuration_sha256": _id._sha256(cfg),
        "finite_xi_max": finite_xi,
        "log_X_p": log_x_p,
        "log_X_source_end": log_x_end,
        "lambda": lam,
    }


def _load(cfg: Mapping[str, Any] | None = None) -> Any:
    try:
        module = importlib.import_module(A1_MODULE)
    except ModuleNotFoundError as exc:
        raise RuntimeError("exact A1 #1107 log-X runtime unavailable") from exc
    cls = getattr(module, A1_CLASS, None)
    if cls is None:
        raise RuntimeError("exact A1 #1107 log-X class unavailable")
    backend = cls() if cfg is None else cls.from_configuration(cfg)
    _a1_identity(backend)
    return backend


def _tail_points(field: "CurrentMainPulseLogXLeadingOscillatoryField"):
    xi0 = field.finite_xi_max
    xi = xi0 + TAIL_FRACTIONS * (SOURCE_MAIN_XI_END - xi0)
    log_x = np.asarray([field.leading_backend.log_X_from_xi(float(v)) for v in xi], dtype=float)
    q = 1.0 - TAIL_T  # eta=0 specialization of the current similarity map.
    log_r = 0.5 * (math.log(2.0) + np.log(q) + log_x)
    if np.any(log_r >= math.log(np.finfo(float).max)):
        raise RuntimeError("A1 log-X tail produced non-representable physical radius")
    radius = np.exp(log_r)
    x = radius * np.cos(TAIL_THETA)
    y = radius * np.sin(TAIL_THETA)
    z = np.zeros_like(radius)
    return x, y, z, TAIL_T.copy(), log_x, xi, radius


@dataclass(frozen=True)
class CurrentMainPulseLogXLeadingOscillatoryField:
    leading_backend: Any
    _leading_identity: dict[str, Any] = field(init=False, repr=False, compare=False)
    _osc_payload: dict[str, Any] = field(init=False, repr=False, compare=False)
    _osc_sha: str = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        _parent_blob()
        identity = _a1_identity(self.leading_backend)
        osc_payload = _id._oscillatory_payload()
        object.__setattr__(self, "_leading_identity", identity)
        object.__setattr__(self, "_osc_payload", osc_payload)
        object.__setattr__(self, "_osc_sha", _id._sha256(osc_payload))

    @property
    def finite_xi_max(self) -> float:
        return float(self._leading_identity["finite_xi_max"])

    @property
    def log_X_source_end(self) -> float:
        return float(self._leading_identity["log_X_source_end"])

    @property
    def oscillatory_runtime_sha256(self) -> str:
        return self._osc_sha

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return _truth()

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return _comp._evaluate_with_backend(self.leading_backend, x, y, z, t).velocity

    def configuration(self) -> dict[str, Any]:
        cfg = self.leading_backend.configuration()
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent2": {
                "pr": PARENT_A2_PR,
                "head": PARENT_A2_HEAD,
                "source_blob_sha1": PARENT_A2_BLOB,
            },
            "agent1": {
                "pr": A1_PR,
                "head": A1_HEAD,
                "module": A1_MODULE,
                "class": A1_CLASS,
                "source_blob_sha1": A1_BLOB,
                "configuration": json.loads(_id._canonical_json(cfg)),
                "configuration_sha256": self._leading_identity["configuration_sha256"],
                "semantic_sha256": self._leading_identity["semantic_sha256"],
            },
            "oscillatory_runtime": {
                "payload": self._osc_payload,
                "payload_sha256": self._osc_sha,
                "public_z_source_blob_sha1": _id.PUBLIC_Z_OSCILLATORY_SOURCE_BLOB_SHA1,
            },
            "source_provenance": SOURCE,
            "composition": "u_main_pulse_logX = u_lead_A1_1107 + u_osc_frozen_complete_curl",
            "truth_boundary": _truth(),
        }

    @property
    def semantic_sha256(self) -> str:
        return _id._sha256(self.configuration())

    def save_candidate(self, path: str | Path) -> dict[str, Any]:
        payload = {
            "schema": SAVED_SCHEMA,
            "configuration": self.configuration(),
            "semantic_sha256": self.semantic_sha256,
        }
        Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return payload

    @classmethod
    def from_configuration(cls, payload: Mapping[str, Any]):
        if (
            not isinstance(payload, Mapping)
            or payload.get("schema") != SCHEMA
            or payload.get("source_provenance") != SOURCE
            or payload.get("truth_boundary") != _truth()
        ):
            raise ValueError("current log-X main-pulse configuration/provenance/truth drifted")
        parent = payload.get("parent_agent2", {})
        expected_parent = {
            "pr": PARENT_A2_PR,
            "head": PARENT_A2_HEAD,
            "source_blob_sha1": PARENT_A2_BLOB,
        }
        if parent != expected_parent:
            raise ValueError("A2 parent identity drifted")
        a1 = payload.get("agent1", {})
        if any(
            (
                a1.get("pr") != A1_PR,
                a1.get("head") != A1_HEAD,
                a1.get("module") != A1_MODULE,
                a1.get("class") != A1_CLASS,
                a1.get("source_blob_sha1") != A1_BLOB,
            )
        ):
            raise ValueError("A1 #1107 provenance drifted")
        cfg = a1.get("configuration")
        cfg_sha = a1.get("configuration_sha256")
        semantic = a1.get("semantic_sha256")
        if not isinstance(cfg, Mapping) or _id._sha256(cfg) != cfg_sha:
            raise ValueError("A1 configuration digest mismatch")
        backend = _load(cfg)
        identity = _a1_identity(backend)
        if identity["configuration_sha256"] != cfg_sha or identity["semantic_sha256"] != semantic:
            raise ValueError("A1 semantic identity drifted")
        osc = payload.get("oscillatory_runtime", {})
        osc_payload = osc.get("payload")
        osc_sha = osc.get("payload_sha256")
        current = _id._oscillatory_payload()
        if (
            not isinstance(osc_payload, Mapping)
            or _id._sha256(osc_payload) != osc_sha
            or _id._sha256(current) != osc_sha
            or _id._canonical_json(osc_payload) != _id._canonical_json(current)
        ):
            raise ValueError("oscillatory runtime drifted")
        obj = cls(backend)
        if obj.configuration() != payload:
            raise ValueError("reconstructed configuration drifted")
        return obj

    @classmethod
    def load_candidate(cls, path: str | Path):
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        cfg = raw.get("configuration")
        if (
            raw.get("schema") != SAVED_SCHEMA
            or not isinstance(cfg, Mapping)
            or _id._sha256(cfg) != raw.get("semantic_sha256")
        ):
            raise ValueError("saved current log-X main-pulse identity invalid")
        obj = cls.from_configuration(cfg)
        if obj.semantic_sha256 != raw["semantic_sha256"]:
            raise ValueError("reloaded semantic identity drifted")
        return obj


@lru_cache(maxsize=1)
def default_field() -> CurrentMainPulseLogXLeadingOscillatoryField:
    return CurrentMainPulseLogXLeadingOscillatoryField(_load())


def velocity(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    return default_field().velocity(x, y, z, t)


def public_contract() -> dict[str, Any]:
    forbidden = {
        "amplitude", "phase", "scale", "orientation", "support", "residual",
        "forcing", "pressure", "viscosity", "gain", "threshold", "mean",
        "stress", "correction", "tolerance", "panels", "steps", "rtol", "atol",
    }
    params = set(inspect.signature(velocity).parameters)
    return {
        "public_inputs": list(inspect.signature(velocity).parameters),
        "forbidden_velocity_inputs_present": sorted(params & forbidden),
        "full_source_xi_11_current_cartesian_leading_materialized": True,
        "source_exact_main_pulse_amplitude_materialized": False,
        "terminal_global_leading_velocity_materialized": False,
        "velocity_export_ready": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def materialize_receipt() -> dict[str, Any]:
    field = default_field()
    inner = _comp._evaluate_with_backend(
        field.leading_backend, INNER[:, 0], INNER[:, 1], INNER[:, 2], INNER_T
    )
    axis = _comp._evaluate_with_backend(
        field.leading_backend, AXIS[:, 0], AXIS[:, 1], AXIS[:, 2], AXIS_T
    )
    tail_x, tail_y, tail_z, tail_t, log_x, xi, radius = _tail_points(field)
    tail = _comp._evaluate_with_backend(field.leading_backend, tail_x, tail_y, tail_z, tail_t)
    outside_osc = np.asarray(_comp.velocity_osc_batch(OUTSIDE, OUTSIDE_T), dtype=float)

    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "candidate.json"
        field.save_candidate(path)
        replayed = CurrentMainPulseLogXLeadingOscillatoryField.load_candidate(path)
        replay = max(
            float(np.max(np.abs(replayed.velocity(INNER[:, 0], INNER[:, 1], INNER[:, 2], INNER_T) - inner.velocity))),
            float(np.max(np.abs(replayed.velocity(tail_x, tail_y, tail_z, tail_t) - tail.velocity))),
        )

    endpoint_xi = SOURCE_MAIN_XI_END
    endpoint_log_x = float(field.leading_backend.log_X_from_xi(endpoint_xi))
    endpoint_t = 0.59
    endpoint_q = 1.0 - endpoint_t
    endpoint_r = math.exp(0.5 * (math.log(2.0) + math.log(endpoint_q) + endpoint_log_x))
    endpoint_velocity = np.asarray(field.velocity(endpoint_r, 0.0, 0.0, endpoint_t), dtype=float)
    endpoint_finite = bool(np.all(np.isfinite(endpoint_velocity)))

    beyond_xi = SOURCE_MAIN_XI_END + 0.05
    beyond_log_x = field._leading_identity["log_X_p"] + beyond_xi / field._leading_identity["lambda"]
    beyond_r = math.exp(0.5 * (math.log(2.0) + math.log(endpoint_q) + beyond_log_x))
    beyond_closed = False
    try:
        field.velocity(beyond_r, 0.0, 0.0, endpoint_t)
    except ValueError:
        beyond_closed = True

    return {
        "schema": RECEIPT_SCHEMA,
        "task": TASK,
        "parent_head": PARENT_A2_HEAD,
        "a1_head": A1_HEAD,
        "source_provenance": SOURCE,
        "semantic_sha256": field.semantic_sha256,
        "domain": {
            "finite_parent_xi_max": field.finite_xi_max,
            "source_main_xi_end": SOURCE_MAIN_XI_END,
            "log_X_source_end": field.log_X_source_end,
            "full_source_xi_11_current_cartesian_leading_materialized": True,
            "terminal_global_leading_velocity_materialized": False,
        },
        "probes": {
            "tail_xi": xi.tolist(),
            "tail_log_X": log_x.tolist(),
            "tail_physical_radius": radius.tolist(),
        },
        "checks": {
            "inner_additive_max_abs": float(np.max(np.abs(inner.velocity - inner.leading - inner.oscillatory))),
            "tail_additive_max_abs": float(np.max(np.abs(tail.velocity - tail.leading - tail.oscillatory))),
            "inner_oscillatory_vector_rms": float(np.sqrt(np.mean(np.sum(inner.oscillatory ** 2, axis=-1)))),
            "tail_oscillatory_vector_rms_observation": float(np.sqrt(np.mean(np.sum(tail.oscillatory ** 2, axis=-1)))),
            "axis_oscillatory_max_abs": float(np.max(np.abs(axis.oscillatory))),
            "outside_support_oscillatory_max_abs": float(np.max(np.abs(outside_osc))),
            "save_load_velocity_replay_max_abs": replay,
            "xi_11_endpoint_velocity_finite": endpoint_finite,
            "post_xi_11_boundary_fail_closed": beyond_closed,
        },
        "truth_boundary": _truth(),
    }


def enforce_receipt(receipt: Mapping[str, Any]) -> None:
    if (
        receipt.get("schema") != RECEIPT_SCHEMA
        or receipt.get("task") != TASK
        or receipt.get("source_provenance") != SOURCE
        or receipt.get("truth_boundary") != _truth()
    ):
        raise AssertionError("receipt identity/truth drifted")
    checks = receipt.get("checks", {})
    domain = receipt.get("domain", {})
    if float(checks.get("inner_additive_max_abs", np.inf)) > COMP_ATOL:
        raise AssertionError("inner additive closure failed")
    if float(checks.get("tail_additive_max_abs", np.inf)) > COMP_ATOL:
        raise AssertionError("log-X tail additive closure failed")
    if float(checks.get("inner_oscillatory_vector_rms", 0.0)) < SIGNAL_FLOOR:
        raise AssertionError("inner oscillation is vacuous")
    if float(checks.get("axis_oscillatory_max_abs", np.inf)) != 0.0:
        raise AssertionError("oscillatory axis contract failed")
    if float(checks.get("outside_support_oscillatory_max_abs", np.inf)) != 0.0:
        raise AssertionError("oscillatory exterior-support contract failed")
    if float(checks.get("save_load_velocity_replay_max_abs", np.inf)) > REPLAY_ATOL:
        raise AssertionError("save/load replay failed")
    if checks.get("xi_11_endpoint_velocity_finite") is not True:
        raise AssertionError("xi=11 endpoint was not materialized")
    if checks.get("post_xi_11_boundary_fail_closed") is not True:
        raise AssertionError("post-xi=11 boundary did not fail closed")
    finite_xi = float(domain.get("finite_parent_xi_max", np.nan))
    source_xi = float(domain.get("source_main_xi_end", np.nan))
    if not (np.isfinite(finite_xi) and 0.0 < finite_xi < source_xi == SOURCE_MAIN_XI_END):
        raise AssertionError("log-X continuation domain identity invalid")
    if domain.get("full_source_xi_11_current_cartesian_leading_materialized") is not True:
        raise AssertionError("A1 xi=11 leading identity missing")
    if domain.get("terminal_global_leading_velocity_materialized") is not False:
        raise AssertionError("terminal/global scope was promoted")
