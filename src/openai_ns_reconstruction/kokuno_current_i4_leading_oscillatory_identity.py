"""Identity-bound current Kokuno I4 leading + frozen complete-curl oscillation.

Kokuno Agent 2 owns only the oscillatory composition in this module.  The
leading field is consumed from exact Agent-1 PR #1079 and is not reimplemented
here::

    u_current_I4(x,y,z,t) = u_lead_A1_1079(x,y,z,t) + u_osc_frozen(x,y,z,t).

The oscillatory summand is the already-frozen bounded complete-curl runtime.
The corrected Kokuno 2026-09-09 reconstruction is structural provenance; this
adapter does not claim a paper-exact oscillation, source I3 positive-order
profiles, the I4 mean correction, a global field, or Navier--Stokes validity.
"""
from __future__ import annotations

from dataclasses import dataclass, field as dataclass_field
from functools import lru_cache
import importlib
import inspect
import json
from pathlib import Path
import tempfile
from typing import Any, Mapping

import numpy as np

from . import kokuno_current_partial_composite_identity_save_load as _identity
from . import kokuno_current_partial_leading_oscillatory_velocity as _composition
from . import kokuno_current_i2_leading_oscillatory_identity as _parent_a2

TASK = "K2-OSC-097"
SCHEMA = "kokuno-a2-current-i4-leading-oscillatory-identity-v1"
SERIALIZED_SCHEMA = "kokuno-a2-current-i4-leading-oscillatory-saved-candidate-v1"
RECEIPT_SCHEMA = "kokuno-a2-current-i4-leading-oscillatory-receipt-v1"

PARENT_AGENT2_PR = 1071
PARENT_AGENT2_HEAD = "48d69e37e78e7f7f0e4e9936f28ff7974719288d"
PARENT_AGENT2_MODULE = "openai_ns_reconstruction.kokuno_current_i2_leading_oscillatory_identity"
PARENT_AGENT2_SOURCE_BLOB_SHA1 = "be68248aca97173a15474324f63efff2c9ffce56"

AGENT1_PR = 1079
AGENT1_HEAD = "b06742ca6e189499192ede3cce40f62cdc1e35ca"
AGENT1_MODULE = "openai_ns_reconstruction.kokuno_pa16_current_cartesian_i4_leading_preservation"
AGENT1_CLASS = "KokunoPA16CurrentCartesianI4LeadingPreservation"
AGENT1_SOURCE_BLOB_SHA1 = "6f04ce0a856b44430402576dad88438da90d1ebb"
AGENT1_TEST_BLOB_SHA1 = "3d876bdfdeca12701b4fb2c300a1f9f4199874f7"
AGENT1_WORKFLOW_BLOB_SHA1 = "48a90b0ac36b4c002663b282208d168a27380c46"

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_BLOB_SHA1 = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_CORRECTED_RELEASE_DATE = "2026-09-09"
PUBLIC_Z_OSCILLATORY_SOURCE_BLOB_SHA1 = _identity.PUBLIC_Z_OSCILLATORY_SOURCE_BLOB_SHA1

COMPOSITION_ATOL = 1.0e-13
VELOCITY_REPLAY_ATOL = 2.0e-12
OSCILLATORY_SIGNAL_FLOOR = 1.0e-12

_INNER_POINTS = np.asarray(
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
_INNER_TIMES = np.asarray((0.31, 0.39, 0.47, 0.55, 0.63, 0.71), dtype=float)
_AXIS_POINTS = np.asarray(((0.0, 0.0, -0.10), (0.0, 0.0, 0.0), (0.0, 0.0, 0.12)))
_AXIS_TIMES = np.asarray((0.37, 0.51, 0.67), dtype=float)
_OUTSIDE_OSC_POINTS = np.asarray(((3.5, 0.0, 0.0), (0.0, -3.5, 0.4), (2.8, 2.8, -0.5)))
_OUTSIDE_OSC_TIMES = np.asarray((0.41, 0.53, 0.69), dtype=float)
_I4_STAGE_FRACTIONS = np.asarray((0.18, 0.57, 0.92), dtype=float)
_I4_ETA = np.asarray((-0.31, 0.08, 0.37), dtype=float)
_I4_TIMES = np.asarray((0.43, 0.55, 0.67), dtype=float)
_I4_THETA = np.asarray((0.23, 0.91, 1.47), dtype=float)


def _require_sha256(value: Any, label: str) -> str:
    if not (
        isinstance(value, str)
        and len(value) == 64
        and all(ch in "0123456789abcdef" for ch in value)
    ):
        raise ValueError(f"{label} must be one lowercase SHA-256 hex digest")
    return value


def _source_provenance() -> dict[str, Any]:
    return {
        "repository": SOURCE_REPOSITORY,
        "commit": SOURCE_COMMIT,
        "path": SOURCE_PATH,
        "blob_sha1": SOURCE_BLOB_SHA1,
        "corrected_release_date": SOURCE_CORRECTED_RELEASE_DATE,
        "classification": "structural_math_provenance_only",
        "scope": "localized waves / complete curls plus I3/I4 source-role organization",
        "paper_exact_claim": False,
    }


def _validate_parent_a2_source() -> str:
    source_path = inspect.getsourcefile(_parent_a2.CurrentI2LeadingOscillatoryField)
    if source_path is None:
        raise RuntimeError("cannot locate exact Agent-2 #1071 parent source")
    blob = _composition._git_blob_sha1(source_path)
    if blob != PARENT_AGENT2_SOURCE_BLOB_SHA1:
        raise RuntimeError("Agent-2 #1071 parent source blob drifted")
    return blob


def _required_agent1_truth() -> dict[str, bool]:
    return {
        "public_reconstruction_source": True,
        "exact_current_I3_parent_consumed": True,
        "source_I4_is_reserved_mean_correction_interval": True,
        "source_leading_E0_U0_unchanged_on_I4": True,
        "source_I4_is_not_a_leading_E0_U0_patch": True,
        "current_leading_preserved_through_I4": True,
        "current_incompressibility_memory_carried_through_I4": True,
        "cartesian_velocity_executable_through_I4": True,
        "velocity_interface_vectorized": True,
        "configuration_serializable": True,
        "source_positive_order_I3_profiles_materialized": False,
        "I3_positive_order_correction_materialized": False,
        "source_I4_mean_correction_materialized": False,
        "current_I4_mean_correction_materialized": False,
        "leading_I4_overlay_invented": False,
        "source_terminal_tail_schedule_bound_into_current_velocity": False,
        "outer_global_leading_velocity_materialized": False,
        "unified_global_cartesian_velocity_export_ready": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "heldout_ns_residual_assessed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def _validate_exact_agent1_backend(backend: Any) -> dict[str, Any]:
    cls = type(backend)
    if cls.__module__ != AGENT1_MODULE or cls.__name__ != AGENT1_CLASS:
        raise RuntimeError("loaded leading backend is not exact Agent-1 #1079")
    source_path = inspect.getsourcefile(cls)
    if source_path is None:
        raise RuntimeError("cannot locate Agent-1 #1079 source")
    source_blob = _composition._git_blob_sha1(source_path)
    if source_blob != AGENT1_SOURCE_BLOB_SHA1:
        raise RuntimeError("Agent-1 #1079 source blob drifted")

    truth = getattr(backend, "truth_boundary", None)
    if not isinstance(truth, Mapping):
        raise RuntimeError("Agent-1 #1079 truth boundary is missing")
    for key, expected in _required_agent1_truth().items():
        if truth.get(key) is not expected:
            raise RuntimeError(f"Agent-1 #1079 truth boundary drifted at {key}")

    semantic = _require_sha256(str(getattr(backend, "semantic_sha256", "")), "Agent-1 semantic SHA")
    configuration = backend.configuration()
    if not isinstance(configuration, Mapping):
        raise RuntimeError("Agent-1 #1079 configuration is missing")
    configuration_sha = _identity._sha256(configuration)

    d = float(getattr(backend, "D", np.nan))
    x_i3_end = float(getattr(backend, "X_I3_end", np.nan))
    x_i4_start = float(getattr(backend, "X_I4_start", np.nan))
    x_i4_end = float(getattr(backend, "X_I4_end", np.nan))
    if not all(np.isfinite(value) for value in (d, x_i3_end, x_i4_start, x_i4_end)):
        raise RuntimeError("Agent-1 #1079 I4 geometry is non-finite")
    if not (d > 0.0 and 0.0 < x_i3_end < x_i4_start < x_i4_end):
        raise RuntimeError("Agent-1 #1079 I4 geometry is invalid")

    return {
        "semantic_sha256": semantic,
        "configuration_sha256": configuration_sha,
        "source_blob_sha1": source_blob,
        "D": d,
        "X_I3_end": x_i3_end,
        "X_I4_start": x_i4_start,
        "X_I4_end": x_i4_end,
    }


@lru_cache(maxsize=1)
def _load_exact_agent1_backend() -> Any:
    try:
        module = importlib.import_module(AGENT1_MODULE)
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "exact Agent-1 #1079 runtime is unavailable; current-I4 composition cannot materialize"
        ) from exc
    cls = getattr(module, AGENT1_CLASS, None)
    if cls is None:
        raise RuntimeError("exact Agent-1 #1079 class is unavailable")
    backend = cls()
    _validate_exact_agent1_backend(backend)
    return backend


def _load_exact_agent1_from_configuration(configuration: Mapping[str, Any]) -> Any:
    try:
        module = importlib.import_module(AGENT1_MODULE)
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "exact Agent-1 #1079 runtime is unavailable; current-I4 candidate cannot reload"
        ) from exc
    cls = getattr(module, AGENT1_CLASS, None)
    if cls is None or not callable(getattr(cls, "from_configuration", None)):
        raise RuntimeError("exact Agent-1 #1079 from_configuration is unavailable")
    backend = cls.from_configuration(configuration)
    _validate_exact_agent1_backend(backend)
    return backend


def _truth_boundary() -> dict[str, bool]:
    return {
        "current_leading_through_I4_consumed": True,
        "current_I4_leading_plus_frozen_complete_curl_oscillation_materialized": True,
        "full_concrete_oscillatory_runtime_digest_bound": True,
        "identity_preserving_current_I4_composite_save_load_available": True,
        "reload_recomputes_leading_and_oscillatory_identities": True,
        "a2_semantic_identity_binds_corrected_source_provenance_block": True,
        "field_summands_changed_by_this_increment": False,
        "source_positive_order_I3_profiles_materialized": False,
        "I3_positive_order_correction_materialized": False,
        "source_I4_mean_correction_materialized": False,
        "current_I4_mean_correction_materialized": False,
        "leading_I4_overlay_invented": False,
        "source_terminal_tail_schedule_bound_into_current_velocity": False,
        "outer_global_leading_velocity_materialized": False,
        "global_compact_support_completed": False,
        "agent3_correction_velocity_composed": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "complete_velocity_pressure_forcing_api": False,
        "heldout_ns_residual_assessed": False,
        "same_protocol_comparable_to_st006": False,
        "residual_reduction_claimed": False,
        "velocity_export_ready": False,
        "visual_correspondence_verified": False,
        "source_forcing_provider_materialized": False,
        "source_background_provider_materialized": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "pde_validated": False,
    }


def _cartesian_i4_probe(
    field: "CurrentI4LeadingOscillatoryField",
    fractions: Any,
    eta: Any,
    t: Any,
    theta: Any,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    fractions, eta, t, theta = np.broadcast_arrays(
        np.asarray(fractions, dtype=float),
        np.asarray(eta, dtype=float),
        np.asarray(t, dtype=float),
        np.asarray(theta, dtype=float),
    )
    if not all(np.all(np.isfinite(v)) for v in (fractions, eta, t, theta)):
        raise ValueError("I4 probe inputs must be finite")
    if np.any((fractions <= 0.0) | (fractions >= 1.0)):
        raise ValueError("I4 stage fractions must lie strictly between zero and one")
    if np.any(np.abs(eta) >= 1.0):
        raise ValueError("I4 probe eta must satisfy |eta|<1")
    if np.any((t <= 0.0) | (t >= 1.0)):
        raise ValueError("I4 probe time must lie strictly between zero and one")

    log_start = np.log(field.X_I4_start)
    log_end = np.log(field.X_I4_end)
    X = np.exp(log_start + fractions * (log_end - log_start))
    tau = 1.0 - t
    q = tau / (1.0 - eta * eta)
    z = np.power(q, field.D) * eta
    r = np.sqrt(2.0 * q * X)
    return r * np.cos(theta), r * np.sin(theta), z, t, X


@dataclass(frozen=True)
class CurrentI4LeadingOscillatoryField:
    """Exact A1 #1079 leading through I4 plus frozen complete-curl oscillation."""

    leading_backend: Any
    _leading_identity: dict[str, Any] = dataclass_field(init=False, repr=False, compare=False)
    _oscillatory_runtime_payload: dict[str, Any] = dataclass_field(init=False, repr=False, compare=False)
    _oscillatory_runtime_sha256: str = dataclass_field(init=False, repr=False, compare=False)
    _parent_agent2_source_blob_sha1: str = dataclass_field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        parent_blob = _validate_parent_a2_source()
        leading_identity = _validate_exact_agent1_backend(self.leading_backend)
        oscillatory_payload = _identity._oscillatory_payload()
        object.__setattr__(self, "_parent_agent2_source_blob_sha1", parent_blob)
        object.__setattr__(self, "_leading_identity", dict(leading_identity))
        object.__setattr__(self, "_oscillatory_runtime_payload", oscillatory_payload)
        object.__setattr__(self, "_oscillatory_runtime_sha256", _identity._sha256(oscillatory_payload))

    @property
    def D(self) -> float:
        return float(self._leading_identity["D"])

    @property
    def X_I3_end(self) -> float:
        return float(self._leading_identity["X_I3_end"])

    @property
    def X_I4_start(self) -> float:
        return float(self._leading_identity["X_I4_start"])

    @property
    def X_I4_end(self) -> float:
        return float(self._leading_identity["X_I4_end"])

    @property
    def oscillatory_runtime_sha256(self) -> str:
        return self._oscillatory_runtime_sha256

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return dict(_truth_boundary())

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return _composition._evaluate_with_backend(self.leading_backend, x, y, z, t).velocity

    def configuration(self) -> dict[str, Any]:
        leading_config = self.leading_backend.configuration()
        if not isinstance(leading_config, Mapping):
            raise RuntimeError("Agent-1 #1079 configuration is missing")
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent2": {
                "pr": PARENT_AGENT2_PR,
                "head": PARENT_AGENT2_HEAD,
                "module": PARENT_AGENT2_MODULE,
                "source_blob_sha1": self._parent_agent2_source_blob_sha1,
            },
            "agent1_leading": {
                "pr": AGENT1_PR,
                "head": AGENT1_HEAD,
                "module": AGENT1_MODULE,
                "class": AGENT1_CLASS,
                "source_blob_sha1": AGENT1_SOURCE_BLOB_SHA1,
                "configuration": json.loads(_identity._canonical_json(leading_config)),
                "configuration_sha256": self._leading_identity["configuration_sha256"],
                "semantic_sha256": self._leading_identity["semantic_sha256"],
            },
            "oscillatory_runtime": {
                "public_z_source_blob_sha1": PUBLIC_Z_OSCILLATORY_SOURCE_BLOB_SHA1,
                "payload": self._oscillatory_runtime_payload,
                "payload_sha256": self._oscillatory_runtime_sha256,
                "binding": "complete concrete default_field().to_payload() realization",
            },
            "source_provenance": _source_provenance(),
            "composition": "u_current_I4 = u_lead_A1_1079 + u_osc_frozen_complete_curl",
            "truth_boundary": self.truth_boundary,
        }

    @property
    def semantic_sha256(self) -> str:
        return _identity._sha256(self.configuration())

    def serialized_payload(self) -> dict[str, Any]:
        return {
            "schema": SERIALIZED_SCHEMA,
            "configuration": self.configuration(),
            "semantic_sha256": self.semantic_sha256,
        }

    def save_candidate(self, path: str | Path) -> dict[str, Any]:
        payload = self.serialized_payload()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return payload

    @classmethod
    def from_configuration(cls, payload: Mapping[str, Any]) -> "CurrentI4LeadingOscillatoryField":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current-I4 composite configuration schema")
        parent = payload.get("parent_agent2")
        if not isinstance(parent, Mapping) or (
            parent.get("pr") != PARENT_AGENT2_PR
            or parent.get("head") != PARENT_AGENT2_HEAD
            or parent.get("module") != PARENT_AGENT2_MODULE
            or parent.get("source_blob_sha1") != PARENT_AGENT2_SOURCE_BLOB_SHA1
        ):
            raise ValueError("Agent-2 #1071 parent identity drifted")
        if payload.get("source_provenance") != _source_provenance():
            raise ValueError("corrected-source provenance drifted")
        if payload.get("composition") != "u_current_I4 = u_lead_A1_1079 + u_osc_frozen_complete_curl":
            raise ValueError("current-I4 composition contract drifted")
        if payload.get("truth_boundary") != _truth_boundary():
            raise ValueError("current-I4 truth boundary drifted")

        leading = payload.get("agent1_leading")
        if not isinstance(leading, Mapping) or (
            leading.get("pr") != AGENT1_PR
            or leading.get("head") != AGENT1_HEAD
            or leading.get("module") != AGENT1_MODULE
            or leading.get("class") != AGENT1_CLASS
            or leading.get("source_blob_sha1") != AGENT1_SOURCE_BLOB_SHA1
        ):
            raise ValueError("Agent-1 #1079 provenance drifted")
        leading_config = leading.get("configuration")
        if not isinstance(leading_config, Mapping):
            raise ValueError("serialized Agent-1 #1079 configuration is missing")
        stored_config_sha = _require_sha256(leading.get("configuration_sha256"), "Agent-1 configuration SHA")
        if _identity._sha256(leading_config) != stored_config_sha:
            raise ValueError("serialized Agent-1 configuration digest mismatch")
        stored_semantic = _require_sha256(leading.get("semantic_sha256"), "Agent-1 semantic SHA")
        backend = _load_exact_agent1_from_configuration(leading_config)
        backend_identity = _validate_exact_agent1_backend(backend)
        if backend_identity["configuration_sha256"] != stored_config_sha:
            raise ValueError("reloaded Agent-1 configuration identity drifted")
        if backend_identity["semantic_sha256"] != stored_semantic:
            raise ValueError("reloaded Agent-1 semantic identity drifted")

        oscillatory = payload.get("oscillatory_runtime")
        if not isinstance(oscillatory, Mapping):
            raise ValueError("serialized oscillatory runtime identity is missing")
        if oscillatory.get("public_z_source_blob_sha1") != PUBLIC_Z_OSCILLATORY_SOURCE_BLOB_SHA1:
            raise ValueError("public-z oscillatory source identity drifted")
        stored_osc_payload = oscillatory.get("payload")
        if not isinstance(stored_osc_payload, Mapping):
            raise ValueError("complete oscillatory runtime payload is missing")
        stored_osc_sha = _require_sha256(oscillatory.get("payload_sha256"), "oscillatory runtime SHA")
        if _identity._sha256(stored_osc_payload) != stored_osc_sha:
            raise ValueError("serialized oscillatory runtime digest mismatch")
        current_osc_payload = _identity._oscillatory_payload()
        if _identity._sha256(current_osc_payload) != stored_osc_sha:
            raise ValueError("current concrete oscillatory runtime identity drifted")
        if _identity._canonical_json(current_osc_payload) != _identity._canonical_json(stored_osc_payload):
            raise ValueError("current concrete oscillatory runtime payload drifted")

        obj = cls(backend)
        if obj.configuration() != payload:
            raise ValueError("reconstructed current-I4 composite configuration drifted")
        return obj

    @classmethod
    def load_candidate(cls, path: str | Path) -> "CurrentI4LeadingOscillatoryField":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, Mapping) or raw.get("schema") != SERIALIZED_SCHEMA:
            raise ValueError("unexpected saved current-I4 candidate schema")
        configuration = raw.get("configuration")
        if not isinstance(configuration, Mapping):
            raise ValueError("saved current-I4 candidate configuration is missing")
        stored_semantic = _require_sha256(raw.get("semantic_sha256"), "saved current-I4 semantic SHA")
        if _identity._sha256(configuration) != stored_semantic:
            raise ValueError("saved current-I4 semantic digest mismatch")
        obj = cls.from_configuration(configuration)
        if obj.semantic_sha256 != stored_semantic:
            raise ValueError("reloaded current-I4 semantic identity drifted")
        return obj


@lru_cache(maxsize=1)
def default_field() -> CurrentI4LeadingOscillatoryField:
    return CurrentI4LeadingOscillatoryField(_load_exact_agent1_backend())


def velocity(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    """Return exact current-I4 leading plus the frozen A2 complete-curl oscillation."""
    return default_field().velocity(x, y, z, t)


def public_contract() -> dict[str, Any]:
    params = inspect.signature(velocity).parameters
    forbidden = {
        "amplitude", "phase", "phase_offset", "scale", "orientation", "support",
        "spatial_step", "derivative_step", "residual", "target", "forcing",
        "pressure", "viscosity", "nu", "gain", "damping", "threshold", "mean",
        "stress", "inverse", "correction", "tolerance", "closure_tolerance",
    }
    return {
        "public_inputs": list(params),
        "forbidden_velocity_inputs_present": sorted(forbidden.intersection(params)),
        "current_leading_through_I4_materialized": True,
        "current_I4_leading_plus_oscillatory_materialized": True,
        "full_oscillatory_runtime_payload_bound": True,
        "identity_preserving_save_load": True,
        "source_I4_mean_correction_materialized": False,
        "field_summands_changed": False,
        "new_oscillatory_parameters": False,
        "leading_profile_reimplemented": False,
        "mean_projection_performed": False,
        "radial_inverse_performed": False,
        "correction_velocity_constructed": False,
        "pressure_or_forcing_added": False,
        "complete_ns_residual": False,
        "velocity_export_ready": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def materialize_current_i4_composite_receipt() -> dict[str, Any]:
    field = default_field()
    inner = _composition._evaluate_with_backend(
        field.leading_backend,
        _INNER_POINTS[:, 0], _INNER_POINTS[:, 1], _INNER_POINTS[:, 2], _INNER_TIMES,
    )
    axis = _composition._evaluate_with_backend(
        field.leading_backend,
        _AXIS_POINTS[:, 0], _AXIS_POINTS[:, 1], _AXIS_POINTS[:, 2], _AXIS_TIMES,
    )
    i4_x, i4_y, i4_z, i4_t, i4_X = _cartesian_i4_probe(
        field, _I4_STAGE_FRACTIONS, _I4_ETA, _I4_TIMES, _I4_THETA
    )
    i4 = _composition._evaluate_with_backend(field.leading_backend, i4_x, i4_y, i4_z, i4_t)
    outside_osc = np.asarray(
        _composition.velocity_osc_batch(_OUTSIDE_OSC_POINTS, _OUTSIDE_OSC_TIMES), dtype=float
    )

    inner_closure = float(np.max(np.abs(inner.velocity - (inner.leading + inner.oscillatory))))
    i4_closure = float(np.max(np.abs(i4.velocity - (i4.leading + i4.oscillatory))))
    inner_osc_rms = float(np.sqrt(np.mean(np.sum(inner.oscillatory * inner.oscillatory, axis=-1))))
    axis_osc_max = float(np.max(np.abs(axis.oscillatory)))
    outside_osc_max = float(np.max(np.abs(outside_osc)))
    i4_osc_rms = float(np.sqrt(np.mean(np.sum(i4.oscillatory * i4.oscillatory, axis=-1))))

    with tempfile.TemporaryDirectory() as tmp:
        candidate_path = Path(tmp) / "current_i4_candidate.json"
        saved = field.save_candidate(candidate_path)
        loaded = CurrentI4LeadingOscillatoryField.load_candidate(candidate_path)
        replay_inner = loaded.velocity(
            _INNER_POINTS[:, 0], _INNER_POINTS[:, 1], _INNER_POINTS[:, 2], _INNER_TIMES
        )
        replay_i4 = loaded.velocity(i4_x, i4_y, i4_z, i4_t)
        replay_max = float(max(
            np.max(np.abs(replay_inner - inner.velocity)),
            np.max(np.abs(replay_i4 - i4.velocity)),
        ))

    return {
        "schema": RECEIPT_SCHEMA,
        "semantic_sha256": field.semantic_sha256,
        "saved_semantic_sha256": saved["semantic_sha256"],
        "oscillatory_runtime_sha256": field.oscillatory_runtime_sha256,
        "agent1_semantic_sha256": field._leading_identity["semantic_sha256"],
        "agent1_configuration_sha256": field._leading_identity["configuration_sha256"],
        "parent_agent2_source_blob_sha1": field._parent_agent2_source_blob_sha1,
        "source_provenance_sha256": _identity._sha256(_source_provenance()),
        "public_contract": public_contract(),
        "diagnostic": {
            "inner_probe_count": int(_INNER_POINTS.shape[0]),
            "i4_probe_count": int(_I4_STAGE_FRACTIONS.shape[0]),
            "i4_stage_fraction": [float(v) for v in _I4_STAGE_FRACTIONS],
            "i4_similarity_X": [float(v) for v in i4_X],
            "X_I4_start": field.X_I4_start,
            "X_I4_end": field.X_I4_end,
            "inner_composition_closure_max_abs": inner_closure,
            "i4_composition_closure_max_abs": i4_closure,
            "inner_oscillatory_vector_rms": inner_osc_rms,
            "i4_oscillatory_vector_rms_observed": i4_osc_rms,
            "axis_oscillatory_abs_max": axis_osc_max,
            "outside_oscillatory_abs_max": outside_osc_max,
            "velocity_save_load_replay_max_abs": replay_max,
        },
        "truth_boundary": field.truth_boundary,
        "frozen_mechanical_gates": {
            "inner_composition_closure_max_abs": COMPOSITION_ATOL,
            "i4_composition_closure_max_abs": COMPOSITION_ATOL,
            "inner_oscillatory_vector_rms_min": OSCILLATORY_SIGNAL_FLOOR,
            "axis_oscillatory_abs_max": 0.0,
            "outside_oscillatory_abs_max": 0.0,
            "velocity_save_load_replay_max_abs": VELOCITY_REPLAY_ATOL,
        },
    }


def enforce_receipt(receipt: Mapping[str, Any]) -> None:
    if not isinstance(receipt, Mapping) or receipt.get("schema") != RECEIPT_SCHEMA:
        raise AssertionError("unexpected current-I4 composite receipt schema")
    _require_sha256(receipt.get("semantic_sha256"), "receipt semantic SHA")
    if receipt.get("saved_semantic_sha256") != receipt.get("semantic_sha256"):
        raise AssertionError("saved/reloaded current-I4 semantic identity drifted")
    _require_sha256(receipt.get("oscillatory_runtime_sha256"), "oscillatory runtime SHA")
    _require_sha256(receipt.get("agent1_semantic_sha256"), "Agent-1 semantic SHA")
    _require_sha256(receipt.get("agent1_configuration_sha256"), "Agent-1 configuration SHA")
    if receipt.get("parent_agent2_source_blob_sha1") != PARENT_AGENT2_SOURCE_BLOB_SHA1:
        raise AssertionError("Agent-2 #1071 parent source identity drifted")
    if receipt.get("source_provenance_sha256") != _identity._sha256(_source_provenance()):
        raise AssertionError("corrected-source provenance identity drifted")
    if receipt.get("public_contract") != public_contract():
        raise AssertionError("public current-I4 velocity contract drifted")
    if receipt.get("truth_boundary") != _truth_boundary():
        raise AssertionError("current-I4 truth boundary drifted")

    diagnostic = receipt.get("diagnostic")
    if not isinstance(diagnostic, Mapping):
        raise AssertionError("current-I4 diagnostic is missing")
    if int(diagnostic.get("inner_probe_count", -1)) != int(_INNER_POINTS.shape[0]):
        raise AssertionError("inner probe count drifted")
    if int(diagnostic.get("i4_probe_count", -1)) != int(_I4_STAGE_FRACTIONS.shape[0]):
        raise AssertionError("I4 probe count drifted")
    if diagnostic.get("i4_stage_fraction") != [float(v) for v in _I4_STAGE_FRACTIONS]:
        raise AssertionError("I4 probe fractions drifted")
    x_start = float(diagnostic.get("X_I4_start", np.nan))
    x_end = float(diagnostic.get("X_I4_end", np.nan))
    xs = np.asarray(diagnostic.get("i4_similarity_X", ()), dtype=float)
    if not (np.isfinite(x_start) and np.isfinite(x_end) and 0.0 < x_start < x_end):
        raise AssertionError("I4 interval is invalid")
    if xs.shape != _I4_STAGE_FRACTIONS.shape or not np.all((xs > x_start) & (xs < x_end)):
        raise AssertionError("I4 probes left the authenticated I4 interval")

    if float(diagnostic.get("inner_composition_closure_max_abs", np.inf)) > COMPOSITION_ATOL:
        raise AssertionError("inner leading+oscillatory addition failed")
    if float(diagnostic.get("i4_composition_closure_max_abs", np.inf)) > COMPOSITION_ATOL:
        raise AssertionError("I4 leading+oscillatory addition failed")
    if float(diagnostic.get("inner_oscillatory_vector_rms", -np.inf)) < OSCILLATORY_SIGNAL_FLOOR:
        raise AssertionError("frozen oscillation lost its nonzero inner signal")
    if not np.isfinite(float(diagnostic.get("i4_oscillatory_vector_rms_observed", np.nan))):
        raise AssertionError("I4 oscillatory observation is non-finite")
    if float(diagnostic.get("axis_oscillatory_abs_max", np.inf)) != 0.0:
        raise AssertionError("axis-safe frozen oscillation leaked onto the axis")
    if float(diagnostic.get("outside_oscillatory_abs_max", np.inf)) != 0.0:
        raise AssertionError("strict-support frozen oscillation leaked outside support")
    if float(diagnostic.get("velocity_save_load_replay_max_abs", np.inf)) > VELOCITY_REPLAY_ATOL:
        raise AssertionError("current-I4 save/load velocity replay drifted")

    expected_gates = {
        "inner_composition_closure_max_abs": COMPOSITION_ATOL,
        "i4_composition_closure_max_abs": COMPOSITION_ATOL,
        "inner_oscillatory_vector_rms_min": OSCILLATORY_SIGNAL_FLOOR,
        "axis_oscillatory_abs_max": 0.0,
        "outside_oscillatory_abs_max": 0.0,
        "velocity_save_load_replay_max_abs": VELOCITY_REPLAY_ATOL,
    }
    if receipt.get("frozen_mechanical_gates") != expected_gates:
        raise AssertionError("current-I4 mechanical gates drifted")


def _main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-output", required=True)
    parser.add_argument("--receipt-output", required=True)
    parser.add_argument("--enforce", action="store_true")
    args = parser.parse_args()
    field = default_field()
    field.save_candidate(args.candidate_output)
    receipt = materialize_current_i4_composite_receipt()
    target = Path(args.receipt_output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.enforce:
        enforce_receipt(receipt)


if __name__ == "__main__":
    _main()
