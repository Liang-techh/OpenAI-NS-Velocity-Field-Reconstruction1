"""Identity-bound current Kokuno I1 leading + frozen complete-curl oscillation.

This Agent-2 increment consumes the exact current modulated Agent-1 #1051
Cartesian leading candidate through I1 and composes it with the already-frozen
Agent-2 axis-safe complete-curl oscillatory runtime::

    u_current_I1(x,y,z,t) = u_lead_A1_1051(x,y,z,t) + u_osc_frozen(x,y,z,t).

Neither summand is changed.  Agent-1 #1051 remains the sole owner of the
leading profile/I1 repair and its immutable 5e-7 closure gate.  The oscillatory
summand remains the existing bounded, strict-support repository realization.
The corrected 2026-09-09 Kokuno reconstruction is structural provenance for the
localized-wave/complete-curl organization only; this composition is not
paper-exact, not globally complete, and not Navier--Stokes validation.
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
from . import kokuno_source_zero_data_amplitude_sensitivity as _parent_a2

TASK = "K2-OSC-095"
SCHEMA = "kokuno-a2-current-i1-leading-oscillatory-identity-v1"
SERIALIZED_SCHEMA = "kokuno-a2-current-i1-leading-oscillatory-saved-candidate-v1"
RECEIPT_SCHEMA = "kokuno-a2-current-i1-leading-oscillatory-receipt-v1"

PARENT_AGENT2_PR = 1052
PARENT_AGENT2_HEAD = "d58cee2bf38bdca8da13a796215a079a569ad7d7"
PARENT_AGENT2_MODULE = "openai_ns_reconstruction.kokuno_source_zero_data_amplitude_sensitivity"
PARENT_AGENT2_SOURCE_BLOB_SHA1 = "eba3c00703bd5763b12e6f48a8c34eb918dec82e"

AGENT1_PR = 1051
AGENT1_HEAD = "ea59dc305b4265fcb0bb0f948c2a45f4477d1ae5"
AGENT1_MODULE = "openai_ns_reconstruction.kokuno_pa16_current_cartesian_i1_frozen_gate"
AGENT1_CLASS = "KokunoPA16CurrentCartesianI1FrozenGate"
AGENT1_SOURCE_BLOB_SHA1 = "fe14dab8dccdbca828d722bdfe2f75f287aa9ea8"
AGENT1_TEST_BLOB_SHA1 = "fd671974c52214771d726c8c0c4fc83152c92f6f"
AGENT1_WORKFLOW_BLOB_SHA1 = "9ba955368a4aa7cbe369efa479b21b556c352d48"
I1_CLOSURE_TOLERANCE = 5.0e-7

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
_I1_STAGE_FRACTIONS = np.asarray((0.18, 0.57, 0.92), dtype=float)
_I1_ETA = np.asarray((-0.31, 0.08, 0.37), dtype=float)
_I1_TIMES = np.asarray((0.43, 0.55, 0.67), dtype=float)
_I1_THETA = np.asarray((0.23, 0.91, 1.47), dtype=float)


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
        "scope": "localized waves / complete curls and their phase/support organization",
        "paper_exact_claim": False,
    }


def _validate_parent_a2_source() -> str:
    source_path = inspect.getsourcefile(_parent_a2.solve_source_zero_data_amplitude_sensitivity)
    if source_path is None:
        raise RuntimeError("cannot locate exact Agent-2 #1052 parent source")
    blob = _composition._git_blob_sha1(source_path)
    if blob != PARENT_AGENT2_SOURCE_BLOB_SHA1:
        raise RuntimeError("Agent-2 #1052 parent source blob drifted")
    return blob


def _required_agent1_truth() -> dict[str, bool]:
    return {
        "current_I1_closure_gate_frozen": True,
        "configuration_mutation_of_I1_closure_gate_fails_closed": True,
        "velocity_formula_changed_by_this_increment": False,
        "PA17_coefficients_changed_by_this_increment": False,
        "current_I2_overlay_applied": False,
        "current_I3_overlay_applied": False,
        "current_I4_overlay_applied": False,
        "outer_global_leading_velocity_materialized": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "heldout_ns_residual_assessed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }


def _validate_exact_agent1_backend(backend: Any) -> dict[str, Any]:
    cls = type(backend)
    if cls.__module__ != AGENT1_MODULE or cls.__name__ != AGENT1_CLASS:
        raise RuntimeError("loaded leading backend is not exact Agent-1 #1051")
    source_path = inspect.getsourcefile(cls)
    if source_path is None:
        raise RuntimeError("cannot locate Agent-1 #1051 source")
    source_blob = _composition._git_blob_sha1(source_path)
    if source_blob != AGENT1_SOURCE_BLOB_SHA1:
        raise RuntimeError("Agent-1 #1051 source blob drifted")

    truth = getattr(backend, "truth_boundary", None)
    if not isinstance(truth, Mapping):
        raise RuntimeError("Agent-1 #1051 truth boundary is missing")
    for key, expected in _required_agent1_truth().items():
        if truth.get(key) is not expected:
            raise RuntimeError(f"Agent-1 #1051 truth boundary drifted at {key}")

    gate = float(getattr(backend, "closure_tolerance", np.nan))
    if not np.isfinite(gate) or gate != I1_CLOSURE_TOLERANCE:
        raise RuntimeError("Agent-1 #1051 frozen I1 closure gate drifted")

    semantic = _require_sha256(str(getattr(backend, "semantic_sha256", "")), "Agent-1 semantic SHA")
    configuration = backend.configuration()
    if not isinstance(configuration, Mapping):
        raise RuntimeError("Agent-1 #1051 configuration is missing")
    configuration_sha = _identity._sha256(configuration)

    d = float(getattr(backend, "D", np.nan))
    x_start = float(getattr(backend, "X_I1_start", np.nan))
    x_end = float(getattr(backend, "X_I1_end", np.nan))
    if not all(np.isfinite(value) for value in (d, x_start, x_end)):
        raise RuntimeError("Agent-1 #1051 I1 geometry is non-finite")
    if not (d > 0.0 and 0.0 < x_start < x_end):
        raise RuntimeError("Agent-1 #1051 I1 geometry is invalid")

    return {
        "semantic_sha256": semantic,
        "configuration_sha256": configuration_sha,
        "source_blob_sha1": source_blob,
        "closure_tolerance": gate,
        "D": d,
        "X_I1_start": x_start,
        "X_I1_end": x_end,
    }


@lru_cache(maxsize=1)
def _load_exact_agent1_backend() -> Any:
    try:
        module = importlib.import_module(AGENT1_MODULE)
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "exact Agent-1 #1051 runtime is unavailable; current-I1 composition cannot materialize"
        ) from exc
    cls = getattr(module, AGENT1_CLASS, None)
    if cls is None:
        raise RuntimeError("exact Agent-1 #1051 class is unavailable")
    backend = cls()
    _validate_exact_agent1_backend(backend)
    return backend


def _load_exact_agent1_from_configuration(configuration: Mapping[str, Any]) -> Any:
    try:
        module = importlib.import_module(AGENT1_MODULE)
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "exact Agent-1 #1051 runtime is unavailable; current-I1 candidate cannot reload"
        ) from exc
    cls = getattr(module, AGENT1_CLASS, None)
    if cls is None or not callable(getattr(cls, "from_configuration", None)):
        raise RuntimeError("exact Agent-1 #1051 from_configuration is unavailable")
    backend = cls.from_configuration(configuration)
    _validate_exact_agent1_backend(backend)
    return backend


def _truth_boundary() -> dict[str, bool]:
    return {
        "current_modulated_leading_through_I1_consumed": True,
        "current_I1_leading_plus_frozen_complete_curl_oscillation_materialized": True,
        "full_concrete_oscillatory_runtime_digest_bound": True,
        "identity_preserving_current_I1_composite_save_load_available": True,
        "reload_recomputes_leading_and_oscillatory_identities": True,
        "a2_semantic_identity_binds_corrected_source_provenance_block": True,
        "field_summands_changed_by_this_increment": False,
        "current_I2_overlay_applied": False,
        "current_I3_overlay_applied": False,
        "current_I4_overlay_applied": False,
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


def _cartesian_i1_probe(
    field: "CurrentI1LeadingOscillatoryField",
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
        raise ValueError("I1 probe inputs must be finite")
    if np.any((fractions <= 0.0) | (fractions >= 1.0)):
        raise ValueError("I1 stage fractions must lie strictly between zero and one")
    if np.any(np.abs(eta) >= 1.0):
        raise ValueError("I1 probe eta must satisfy |eta|<1")
    if np.any((t <= 0.0) | (t >= 1.0)):
        raise ValueError("I1 probe time must lie strictly between zero and one")

    log_start = np.log(field.X_I1_start)
    log_end = np.log(field.X_I1_end)
    X = np.exp(log_start + fractions * (log_end - log_start))
    tau = 1.0 - t
    q = tau / (1.0 - eta * eta)
    z = np.power(q, field.D) * eta
    r = np.sqrt(2.0 * q * X)
    return r * np.cos(theta), r * np.sin(theta), z, t, X


@dataclass(frozen=True)
class CurrentI1LeadingOscillatoryField:
    """Exact A1 #1051 leading through I1 plus frozen A2 complete-curl oscillation."""

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
    def X_I1_start(self) -> float:
        return float(self._leading_identity["X_I1_start"])

    @property
    def X_I1_end(self) -> float:
        return float(self._leading_identity["X_I1_end"])

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
            raise RuntimeError("Agent-1 #1051 configuration is missing")
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
                "frozen_i1_closure_tolerance": self._leading_identity["closure_tolerance"],
            },
            "oscillatory_runtime": {
                "public_z_source_blob_sha1": PUBLIC_Z_OSCILLATORY_SOURCE_BLOB_SHA1,
                "payload": self._oscillatory_runtime_payload,
                "payload_sha256": self._oscillatory_runtime_sha256,
                "binding": "complete concrete default_field().to_payload() realization",
            },
            "source_provenance": _source_provenance(),
            "composition": "u_current_I1 = u_lead_A1_1051 + u_osc_frozen_complete_curl",
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
    def from_configuration(cls, payload: Mapping[str, Any]) -> "CurrentI1LeadingOscillatoryField":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current-I1 composite configuration schema")

        parent = payload.get("parent_agent2")
        if not isinstance(parent, Mapping) or (
            parent.get("pr") != PARENT_AGENT2_PR
            or parent.get("head") != PARENT_AGENT2_HEAD
            or parent.get("module") != PARENT_AGENT2_MODULE
            or parent.get("source_blob_sha1") != PARENT_AGENT2_SOURCE_BLOB_SHA1
        ):
            raise ValueError("Agent-2 #1052 parent identity drifted")

        if payload.get("source_provenance") != _source_provenance():
            raise ValueError("corrected-source provenance drifted")
        if payload.get("composition") != "u_current_I1 = u_lead_A1_1051 + u_osc_frozen_complete_curl":
            raise ValueError("current-I1 composition contract drifted")
        if payload.get("truth_boundary") != _truth_boundary():
            raise ValueError("current-I1 truth boundary drifted")

        leading = payload.get("agent1_leading")
        if not isinstance(leading, Mapping) or (
            leading.get("pr") != AGENT1_PR
            or leading.get("head") != AGENT1_HEAD
            or leading.get("module") != AGENT1_MODULE
            or leading.get("class") != AGENT1_CLASS
            or leading.get("source_blob_sha1") != AGENT1_SOURCE_BLOB_SHA1
            or leading.get("frozen_i1_closure_tolerance") != I1_CLOSURE_TOLERANCE
        ):
            raise ValueError("Agent-1 #1051 provenance or frozen gate drifted")
        leading_config = leading.get("configuration")
        if not isinstance(leading_config, Mapping):
            raise ValueError("serialized Agent-1 #1051 configuration is missing")
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
            raise ValueError("reconstructed current-I1 composite configuration drifted")
        return obj

    @classmethod
    def load_candidate(cls, path: str | Path) -> "CurrentI1LeadingOscillatoryField":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, Mapping) or raw.get("schema") != SERIALIZED_SCHEMA:
            raise ValueError("unexpected saved current-I1 candidate schema")
        configuration = raw.get("configuration")
        if not isinstance(configuration, Mapping):
            raise ValueError("saved current-I1 candidate configuration is missing")
        stored_semantic = _require_sha256(raw.get("semantic_sha256"), "saved current-I1 semantic SHA")
        if _identity._sha256(configuration) != stored_semantic:
            raise ValueError("saved current-I1 semantic digest mismatch")
        obj = cls.from_configuration(configuration)
        if obj.semantic_sha256 != stored_semantic:
            raise ValueError("reloaded current-I1 semantic identity drifted")
        return obj


@lru_cache(maxsize=1)
def default_field() -> CurrentI1LeadingOscillatoryField:
    return CurrentI1LeadingOscillatoryField(_load_exact_agent1_backend())


def velocity(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    """Return exact current-I1 leading plus the frozen A2 complete-curl oscillation."""
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
        "current_modulated_leading_through_I1_materialized": True,
        "current_I1_leading_plus_oscillatory_materialized": True,
        "full_oscillatory_runtime_payload_bound": True,
        "identity_preserving_save_load": True,
        "frozen_agent1_i1_closure_tolerance": I1_CLOSURE_TOLERANCE,
        "a2_support_axis_mask_preserved": True,
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


def materialize_current_i1_composite_receipt() -> dict[str, Any]:
    field = default_field()
    inner = _composition._evaluate_with_backend(
        field.leading_backend,
        _INNER_POINTS[:, 0],
        _INNER_POINTS[:, 1],
        _INNER_POINTS[:, 2],
        _INNER_TIMES,
    )
    axis = _composition._evaluate_with_backend(
        field.leading_backend,
        _AXIS_POINTS[:, 0],
        _AXIS_POINTS[:, 1],
        _AXIS_POINTS[:, 2],
        _AXIS_TIMES,
    )
    i1_x, i1_y, i1_z, i1_t, i1_X = _cartesian_i1_probe(
        field, _I1_STAGE_FRACTIONS, _I1_ETA, _I1_TIMES, _I1_THETA
    )
    i1 = _composition._evaluate_with_backend(field.leading_backend, i1_x, i1_y, i1_z, i1_t)

    outside_osc = np.asarray(
        _composition.velocity_osc_batch(_OUTSIDE_OSC_POINTS, _OUTSIDE_OSC_TIMES),
        dtype=float,
    )

    inner_closure = float(np.max(np.abs(inner.velocity - (inner.leading + inner.oscillatory))))
    i1_closure = float(np.max(np.abs(i1.velocity - (i1.leading + i1.oscillatory))))
    inner_osc_rms = float(np.sqrt(np.mean(np.sum(inner.oscillatory * inner.oscillatory, axis=-1))))
    axis_osc_max = float(np.max(np.abs(axis.oscillatory)))
    outside_osc_max = float(np.max(np.abs(outside_osc)))

    with tempfile.TemporaryDirectory() as tmp:
        candidate_path = Path(tmp) / "current_i1_candidate.json"
        saved = field.save_candidate(candidate_path)
        loaded = CurrentI1LeadingOscillatoryField.load_candidate(candidate_path)
        replay_inner = loaded.velocity(
            _INNER_POINTS[:, 0], _INNER_POINTS[:, 1], _INNER_POINTS[:, 2], _INNER_TIMES
        )
        replay_i1 = loaded.velocity(i1_x, i1_y, i1_z, i1_t)
        replay_max = float(
            max(
                np.max(np.abs(replay_inner - inner.velocity)),
                np.max(np.abs(replay_i1 - i1.velocity)),
            )
        )

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
            "i1_probe_count": int(_I1_STAGE_FRACTIONS.shape[0]),
            "i1_stage_fraction": [float(v) for v in _I1_STAGE_FRACTIONS],
            "i1_similarity_X": [float(v) for v in i1_X],
            "X_I1_start": field.X_I1_start,
            "X_I1_end": field.X_I1_end,
            "frozen_i1_closure_tolerance": float(field._leading_identity["closure_tolerance"]),
            "inner_composition_closure_max_abs": inner_closure,
            "i1_composition_closure_max_abs": i1_closure,
            "inner_oscillatory_vector_rms": inner_osc_rms,
            "axis_oscillatory_abs_max": axis_osc_max,
            "outside_oscillatory_abs_max": outside_osc_max,
            "velocity_save_load_replay_max_abs": replay_max,
        },
        "truth_boundary": field.truth_boundary,
        "frozen_mechanical_gates": {
            "inner_composition_closure_max_abs": COMPOSITION_ATOL,
            "i1_composition_closure_max_abs": COMPOSITION_ATOL,
            "inner_oscillatory_vector_rms_min": OSCILLATORY_SIGNAL_FLOOR,
            "axis_oscillatory_abs_max": 0.0,
            "outside_oscillatory_abs_max": 0.0,
            "velocity_save_load_replay_max_abs": VELOCITY_REPLAY_ATOL,
            "frozen_i1_closure_tolerance": I1_CLOSURE_TOLERANCE,
        },
    }


def enforce_receipt(receipt: Mapping[str, Any]) -> None:
    if not isinstance(receipt, Mapping) or receipt.get("schema") != RECEIPT_SCHEMA:
        raise AssertionError("unexpected current-I1 composite receipt schema")
    _require_sha256(receipt.get("semantic_sha256"), "receipt semantic SHA")
    if receipt.get("saved_semantic_sha256") != receipt.get("semantic_sha256"):
        raise AssertionError("saved/reloaded current-I1 semantic identity drifted")
    _require_sha256(receipt.get("oscillatory_runtime_sha256"), "oscillatory runtime SHA")
    _require_sha256(receipt.get("agent1_semantic_sha256"), "Agent-1 semantic SHA")
    _require_sha256(receipt.get("agent1_configuration_sha256"), "Agent-1 configuration SHA")
    if receipt.get("parent_agent2_source_blob_sha1") != PARENT_AGENT2_SOURCE_BLOB_SHA1:
        raise AssertionError("Agent-2 #1052 parent source identity drifted")
    if receipt.get("source_provenance_sha256") != _identity._sha256(_source_provenance()):
        raise AssertionError("corrected-source provenance identity drifted")
    if receipt.get("public_contract") != public_contract():
        raise AssertionError("public current-I1 velocity contract drifted")
    if receipt.get("truth_boundary") != _truth_boundary():
        raise AssertionError("current-I1 truth boundary drifted")

    diagnostic = receipt.get("diagnostic")
    if not isinstance(diagnostic, Mapping):
        raise AssertionError("current-I1 diagnostic is missing")
    if int(diagnostic.get("inner_probe_count", -1)) != int(_INNER_POINTS.shape[0]):
        raise AssertionError("inner probe count drifted")
    if int(diagnostic.get("i1_probe_count", -1)) != int(_I1_STAGE_FRACTIONS.shape[0]):
        raise AssertionError("I1 probe count drifted")
    if diagnostic.get("i1_stage_fraction") != [float(v) for v in _I1_STAGE_FRACTIONS]:
        raise AssertionError("I1 probe fractions drifted")
    x_start = float(diagnostic.get("X_I1_start", np.nan))
    x_end = float(diagnostic.get("X_I1_end", np.nan))
    xs = np.asarray(diagnostic.get("i1_similarity_X", ()), dtype=float)
    if not (np.isfinite(x_start) and np.isfinite(x_end) and 0.0 < x_start < x_end):
        raise AssertionError("I1 interval is invalid")
    if xs.shape != _I1_STAGE_FRACTIONS.shape or not np.all((xs > x_start) & (xs < x_end)):
        raise AssertionError("I1 probes left the authenticated I1 interval")

    if float(diagnostic.get("frozen_i1_closure_tolerance", np.nan)) != I1_CLOSURE_TOLERANCE:
        raise AssertionError("frozen Agent-1 I1 closure tolerance drifted")
    if float(diagnostic.get("inner_composition_closure_max_abs", np.inf)) > COMPOSITION_ATOL:
        raise AssertionError("inner leading+oscillatory addition failed")
    if float(diagnostic.get("i1_composition_closure_max_abs", np.inf)) > COMPOSITION_ATOL:
        raise AssertionError("I1 leading+oscillatory addition failed")
    if float(diagnostic.get("inner_oscillatory_vector_rms", -np.inf)) < OSCILLATORY_SIGNAL_FLOOR:
        raise AssertionError("frozen oscillation lost its nonzero inner signal")
    if float(diagnostic.get("axis_oscillatory_abs_max", np.inf)) != 0.0:
        raise AssertionError("axis-safe frozen oscillation leaked onto the axis")
    if float(diagnostic.get("outside_oscillatory_abs_max", np.inf)) != 0.0:
        raise AssertionError("strict-support frozen oscillation leaked outside support")
    if float(diagnostic.get("velocity_save_load_replay_max_abs", np.inf)) > VELOCITY_REPLAY_ATOL:
        raise AssertionError("current-I1 save/load velocity replay drifted")

    gates = receipt.get("frozen_mechanical_gates")
    expected_gates = {
        "inner_composition_closure_max_abs": COMPOSITION_ATOL,
        "i1_composition_closure_max_abs": COMPOSITION_ATOL,
        "inner_oscillatory_vector_rms_min": OSCILLATORY_SIGNAL_FLOOR,
        "axis_oscillatory_abs_max": 0.0,
        "outside_oscillatory_abs_max": 0.0,
        "velocity_save_load_replay_max_abs": VELOCITY_REPLAY_ATOL,
        "frozen_i1_closure_tolerance": I1_CLOSURE_TOLERANCE,
    }
    if gates != expected_gates:
        raise AssertionError("current-I1 mechanical gates drifted")


def _main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-output", required=True)
    parser.add_argument("--receipt-output", required=True)
    parser.add_argument("--enforce", action="store_true")
    args = parser.parse_args()

    field = default_field()
    field.save_candidate(args.candidate_output)
    receipt = materialize_current_i1_composite_receipt()
    target = Path(args.receipt_output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.enforce:
        enforce_receipt(receipt)


if __name__ == "__main__":
    _main()
