"""Identity-bound current Kokuno leading + complete-curl oscillation through ``X_R``.

Kokuno Agent 1 PR #980 extends the current candidate-side Cartesian leading field
from the already materialized ``X_h`` seam through the public exterior-preservation
interval ``X_h < X <= X_R``.  Agent 2 PR #981, meanwhile, binds the complete
concrete oscillatory runtime and deterministic save/load, but only for the older
Agent-1 #965 leading field through ``X_h``.

This module closes exactly that composition gap:

    u_current_to_XR(x,t) = u_lead,current-through-XR(x,t) + u_osc,frozen-complete-curl(x,t).

It does not change either velocity summand.  The Agent-1 implementation is loaded
from the exact #980 class/blob and the Agent-2 oscillation is the already-frozen
axis-safe public runtime.  The complete oscillatory ``to_payload()`` realization is
included in the semantic identity and in deterministic save/load.

The corrected 2026-09-09 Kokuno reconstruction remains structural provenance.
The current joined/exterior leading realization, public-z oscillatory pullback,
autonomous bounded oscillatory parameters, support mask, numerical bookkeeping,
and this composition are repository realizations.  This field is therefore not
paper-exact, not post-``X_R``/global, not pressure/forcing complete, and not a
Navier--Stokes validation.
"""
from __future__ import annotations

from dataclasses import dataclass, field as dataclass_field
from functools import lru_cache
import importlib
import inspect
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from . import kokuno_current_partial_composite_identity_save_load as _identity
from . import kokuno_current_partial_leading_oscillatory_velocity as _composition

TASK = "K2-OSC-083"
SCHEMA = "kokuno-a2-current-exterior-xr-leading-oscillatory-identity-v1"
SERIALIZED_SCHEMA = "kokuno-a2-current-exterior-xr-leading-oscillatory-saved-candidate-v1"
RECEIPT_SCHEMA = "kokuno-a2-current-exterior-xr-leading-oscillatory-receipt-v1"

PARENT_AGENT2_PR = 981
PARENT_AGENT2_HEAD = "9997fc55455d126f935643da36bf17eaa0491aa4"
AGENT1_PR = 980
AGENT1_HEAD = "d3c971f2c62e272333e124e532212d23cca4908d"
AGENT1_MODULE = "openai_ns_reconstruction.kokuno_pa16_current_cartesian_exterior_to_xr"
AGENT1_CLASS = "KokunoPA16CurrentCartesianExteriorToXR"
AGENT1_SOURCE_BLOB_SHA1 = "2ba28636a56b252fa485719e7b3e8da754d7e588"
AGENT1_TEST_BLOB_SHA1 = "ecbab9822c54f440f5555f25cd5e7e1648d84bd2"
AGENT1_WORKFLOW_BLOB_SHA1 = "5c2d68211e8fea526e06087fc4475f7bab9eb563"

PUBLIC_Z_OSCILLATORY_SOURCE_BLOB_SHA1 = _identity.PUBLIC_Z_OSCILLATORY_SOURCE_BLOB_SHA1
SOURCE_CORRECTED_READER_COMMIT = _composition.SOURCE_CORRECTED_READER_COMMIT
SOURCE_CORRECTED_READER_DATE = _composition.SOURCE_CORRECTED_READER_DATE

VELOCITY_REPLAY_ATOL = 2.0e-12
COMPOSITION_ATOL = 1.0e-13
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


def _require_sha256(value: Any, label: str) -> str:
    if not (
        isinstance(value, str)
        and len(value) == 64
        and all(ch in "0123456789abcdef" for ch in value)
    ):
        raise ValueError(f"{label} must be one lowercase SHA-256 hex digest")
    return value


def _required_agent1_truth() -> dict[str, bool]:
    return {
        "current_cartesian_leading_through_Xh_consumed": True,
        "candidate_side_exterior_preservation_Xh_to_XR_materialized": True,
        "current_cartesian_spacetime_leading_velocity_through_XR_materialized": True,
        "velocity_beyond_Xh_materialized": True,
        "current_incompressibility_memory_carried": True,
        "velocity_interface_vectorized": True,
        "velocity_configuration_serializable": True,
        "source_exact_exterior_preservation_certified": False,
        "post_XR_RF40_current_lineage_materialized": False,
        "outer_global_leading_velocity_materialized": False,
        "unified_global_cartesian_velocity_export_ready": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "heldout_ns_residual_assessed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "paper_exact": False,
    }


def _validate_exact_agent1_backend(backend: Any) -> dict[str, Any]:
    cls = type(backend)
    if cls.__module__ != AGENT1_MODULE or cls.__name__ != AGENT1_CLASS:
        raise RuntimeError("loaded leading backend is not the pinned Agent-1 #980 class")
    source_path = inspect.getsourcefile(cls)
    if source_path is None:
        raise RuntimeError("cannot locate Agent-1 #980 leading implementation source")
    blob = _composition._git_blob_sha1(source_path)
    if blob != AGENT1_SOURCE_BLOB_SHA1:
        raise RuntimeError("Agent-1 #980 leading implementation source blob drifted")

    truth = getattr(backend, "truth_boundary", None)
    if not isinstance(truth, Mapping):
        raise RuntimeError("Agent-1 #980 truth boundary is missing")
    for key, expected in _required_agent1_truth().items():
        if truth.get(key) is not expected:
            raise RuntimeError(f"Agent-1 #980 truth boundary drifted at {key}")

    semantic = str(getattr(backend, "semantic_sha256", ""))
    _require_sha256(semantic, "Agent-1 #980 semantic SHA")
    configuration = backend.configuration()
    if not isinstance(configuration, Mapping):
        raise RuntimeError("Agent-1 #980 configuration is missing")
    config_sha = _identity._sha256(configuration)

    X_h = float(getattr(backend, "X_h", np.nan))
    X_R = float(getattr(backend, "X_R", np.nan))
    D = float(getattr(backend, "D", np.nan))
    if not all(np.isfinite(v) for v in (X_h, X_R, D)) or not (0.0 < X_h < X_R):
        raise RuntimeError("Agent-1 #980 X_h/X_R/D geometry is invalid")
    if not math.isclose(X_h / X_R, math.exp(-5.0), rel_tol=2.0e-14, abs_tol=0.0):
        raise RuntimeError("Agent-1 #980 X_h/X_R identity drifted")

    return {
        "semantic_sha256": semantic,
        "configuration_sha256": config_sha,
        "source_blob_sha1": blob,
        "X_h": X_h,
        "X_R": X_R,
        "D": D,
    }


@lru_cache(maxsize=1)
def _load_exact_agent1_backend() -> Any:
    try:
        module = importlib.import_module(AGENT1_MODULE)
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "exact Agent-1 #980 runtime is unavailable; cannot materialize the through-X_R composite"
        ) from exc
    cls = getattr(module, AGENT1_CLASS, None)
    if cls is None:
        raise RuntimeError("exact Agent-1 #980 class is unavailable")
    backend = cls()
    _validate_exact_agent1_backend(backend)
    return backend


def _load_exact_agent1_from_configuration(configuration: Mapping[str, Any]) -> Any:
    try:
        module = importlib.import_module(AGENT1_MODULE)
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "exact Agent-1 #980 runtime is unavailable; cannot reload the through-X_R composite"
        ) from exc
    cls = getattr(module, AGENT1_CLASS, None)
    if cls is None or not callable(getattr(cls, "from_configuration", None)):
        raise RuntimeError("exact Agent-1 #980 class/from_configuration is unavailable")
    backend = cls.from_configuration(configuration)
    _validate_exact_agent1_backend(backend)
    return backend


def _truth_boundary() -> dict[str, bool]:
    return {
        "current_leading_plus_oscillatory_velocity_through_XR_materialized": True,
        "velocity_beyond_Xh_through_XR_materialized": True,
        "full_concrete_oscillatory_runtime_digest_bound": True,
        "identity_preserving_through_XR_composite_save_load_available": True,
        "reload_recomputes_leading_and_oscillatory_identities": True,
        "post_XR_velocity_materialized": False,
        "post_XR_RF40_current_lineage_materialized": False,
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
        "paper_exact": False,
        "openai_field_identified": False,
        "pde_validated": False,
    }


def _cartesian_probe(
    field: "CurrentExteriorXRLeadingOscillatoryField",
    X: Any,
    eta: Any,
    t: Any,
    theta: Any,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    X_arr, eta_arr, t_arr, theta_arr = np.broadcast_arrays(
        np.asarray(X, dtype=float),
        np.asarray(eta, dtype=float),
        np.asarray(t, dtype=float),
        np.asarray(theta, dtype=float),
    )
    if not all(np.all(np.isfinite(a)) for a in (X_arr, eta_arr, t_arr, theta_arr)):
        raise ValueError("probe coordinates must be finite")
    if np.any((t_arr <= 0.0) | (t_arr >= 1.0)):
        raise ValueError("probe time must lie strictly between zero and one")
    if np.any(np.abs(eta_arr) >= 1.0):
        raise ValueError("probe eta must satisfy |eta|<1")
    tau = 1.0 - t_arr
    q = tau / (1.0 - eta_arr * eta_arr)
    z = np.power(q, field.D) * eta_arr
    r = np.sqrt(2.0 * q * X_arr)
    return r * np.cos(theta_arr), r * np.sin(theta_arr), z, t_arr


@dataclass(frozen=True)
class CurrentExteriorXRLeadingOscillatoryField:
    """Exact #980 current leading through ``X_R`` plus the frozen A2 oscillation."""

    leading_backend: Any
    _leading_identity: dict[str, Any] = dataclass_field(init=False, repr=False, compare=False)
    _oscillatory_runtime_payload: dict[str, Any] = dataclass_field(
        init=False, repr=False, compare=False
    )
    _oscillatory_runtime_sha256: str = dataclass_field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        identity = _validate_exact_agent1_backend(self.leading_backend)
        oscillatory_payload = _identity._oscillatory_payload()
        object.__setattr__(self, "_leading_identity", dict(identity))
        object.__setattr__(self, "_oscillatory_runtime_payload", oscillatory_payload)
        object.__setattr__(
            self, "_oscillatory_runtime_sha256", _identity._sha256(oscillatory_payload)
        )

    @property
    def X_h(self) -> float:
        return float(self._leading_identity["X_h"])

    @property
    def X_R(self) -> float:
        return float(self._leading_identity["X_R"])

    @property
    def D(self) -> float:
        return float(self._leading_identity["D"])

    @property
    def oscillatory_runtime_sha256(self) -> str:
        return self._oscillatory_runtime_sha256

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return _truth_boundary()

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return _composition._evaluate_with_backend(
            self.leading_backend, x, y, z, t
        ).velocity

    def configuration(self) -> dict[str, Any]:
        leading_configuration = self.leading_backend.configuration()
        if not isinstance(leading_configuration, Mapping):
            raise RuntimeError("Agent-1 #980 leading configuration is missing")
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent2": {"pr": PARENT_AGENT2_PR, "head": PARENT_AGENT2_HEAD},
            "agent1_leading": {
                "pr": AGENT1_PR,
                "head": AGENT1_HEAD,
                "module": AGENT1_MODULE,
                "class": AGENT1_CLASS,
                "source_blob_sha1": AGENT1_SOURCE_BLOB_SHA1,
                "configuration": json.loads(_identity._canonical_json(leading_configuration)),
                "configuration_sha256": self._leading_identity["configuration_sha256"],
                "semantic_sha256": self._leading_identity["semantic_sha256"],
                "X_h": self.X_h,
                "X_R": self.X_R,
            },
            "oscillatory_runtime": {
                "public_z_source_blob_sha1": PUBLIC_Z_OSCILLATORY_SOURCE_BLOB_SHA1,
                "payload": self._oscillatory_runtime_payload,
                "payload_sha256": self._oscillatory_runtime_sha256,
                "binding": "complete concrete default_field().to_payload() realization",
            },
            "composition": "u_current_to_XR = u_lead_current_through_XR + u_osc_frozen_complete_curl",
            "source_provenance": {
                "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
                "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
                "scope": "localized waves / complete curls plus public leading-coordinate structure",
            },
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
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "CurrentExteriorXRLeadingOscillatoryField":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected through-X_R composite configuration schema")
        parent = payload.get("parent_agent2")
        if not isinstance(parent, Mapping) or parent.get("pr") != PARENT_AGENT2_PR or parent.get("head") != PARENT_AGENT2_HEAD:
            raise ValueError("through-X_R composite parent identity drifted")

        leading = payload.get("agent1_leading")
        if not isinstance(leading, Mapping):
            raise ValueError("serialized Agent-1 #980 identity is missing")
        if (
            leading.get("pr") != AGENT1_PR
            or leading.get("head") != AGENT1_HEAD
            or leading.get("module") != AGENT1_MODULE
            or leading.get("class") != AGENT1_CLASS
            or leading.get("source_blob_sha1") != AGENT1_SOURCE_BLOB_SHA1
        ):
            raise ValueError("serialized Agent-1 #980 provenance drifted")
        leading_configuration = leading.get("configuration")
        if not isinstance(leading_configuration, Mapping):
            raise ValueError("serialized Agent-1 #980 configuration is missing")
        stored_config_sha = _require_sha256(
            leading.get("configuration_sha256"), "Agent-1 #980 configuration SHA"
        )
        if _identity._sha256(leading_configuration) != stored_config_sha:
            raise ValueError("serialized Agent-1 #980 configuration digest mismatch")
        stored_semantic = _require_sha256(
            leading.get("semantic_sha256"), "Agent-1 #980 semantic SHA"
        )
        backend = _load_exact_agent1_from_configuration(leading_configuration)
        identity = _validate_exact_agent1_backend(backend)
        if identity["configuration_sha256"] != stored_config_sha:
            raise ValueError("reloaded Agent-1 #980 configuration identity drifted")
        if identity["semantic_sha256"] != stored_semantic:
            raise ValueError("reloaded Agent-1 #980 semantic identity drifted")
        if not math.isclose(float(leading.get("X_h", np.nan)), identity["X_h"], rel_tol=0.0, abs_tol=0.0):
            raise ValueError("reloaded Agent-1 #980 X_h identity drifted")
        if not math.isclose(float(leading.get("X_R", np.nan)), identity["X_R"], rel_tol=0.0, abs_tol=0.0):
            raise ValueError("reloaded Agent-1 #980 X_R identity drifted")

        oscillatory = payload.get("oscillatory_runtime")
        if not isinstance(oscillatory, Mapping):
            raise ValueError("serialized oscillatory runtime identity is missing")
        if oscillatory.get("public_z_source_blob_sha1") != PUBLIC_Z_OSCILLATORY_SOURCE_BLOB_SHA1:
            raise ValueError("serialized public-z oscillatory source identity drifted")
        stored_osc_payload = oscillatory.get("payload")
        if not isinstance(stored_osc_payload, Mapping):
            raise ValueError("serialized complete oscillatory runtime payload is missing")
        stored_osc_sha = _require_sha256(
            oscillatory.get("payload_sha256"), "oscillatory runtime payload SHA"
        )
        if _identity._sha256(stored_osc_payload) != stored_osc_sha:
            raise ValueError("serialized oscillatory runtime payload digest mismatch")
        current_osc_payload = _identity._oscillatory_payload()
        if _identity._sha256(current_osc_payload) != stored_osc_sha:
            raise ValueError("current concrete oscillatory runtime identity drifted")
        if _identity._canonical_json(current_osc_payload) != _identity._canonical_json(stored_osc_payload):
            raise ValueError("current concrete oscillatory runtime payload drifted")

        obj = cls(backend)
        if obj.oscillatory_runtime_sha256 != stored_osc_sha:
            raise ValueError("reloaded oscillatory runtime SHA drifted")
        return obj

    @classmethod
    def load_candidate(cls, path: str | Path) -> "CurrentExteriorXRLeadingOscillatoryField":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, Mapping) or raw.get("schema") != SERIALIZED_SCHEMA:
            raise ValueError("unexpected saved through-X_R composite candidate schema")
        configuration = raw.get("configuration")
        if not isinstance(configuration, Mapping):
            raise ValueError("saved through-X_R composite configuration is missing")
        stored_semantic = _require_sha256(
            raw.get("semantic_sha256"), "saved through-X_R composite semantic SHA"
        )
        if _identity._sha256(configuration) != stored_semantic:
            raise ValueError("saved through-X_R composite semantic digest mismatch")
        obj = cls.from_configuration(configuration)
        if obj.semantic_sha256 != stored_semantic:
            raise ValueError("reloaded through-X_R composite semantic identity drifted")
        return obj


@lru_cache(maxsize=1)
def default_field() -> CurrentExteriorXRLeadingOscillatoryField:
    return CurrentExteriorXRLeadingOscillatoryField(_load_exact_agent1_backend())


def velocity(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    """Return the current A1-through-X_R leading field plus frozen A2 oscillation."""
    return default_field().velocity(x, y, z, t)


def public_contract() -> dict[str, Any]:
    params = inspect.signature(velocity).parameters
    forbidden = {
        "amplitude", "phase", "phase_offset", "scale", "orientation", "support",
        "spatial_step", "derivative_step", "residual", "target", "forcing",
        "pressure", "viscosity", "nu", "gain", "damping", "threshold", "mean",
        "stress", "inverse", "correction",
    }
    return {
        "public_inputs": list(params),
        "forbidden_velocity_inputs_present": sorted(forbidden.intersection(params)),
        "exact_agent1_pr": AGENT1_PR,
        "leading_materialized_through_XR": True,
        "full_oscillatory_runtime_payload_bound": True,
        "identity_preserving_save_load": True,
        "post_XR_fail_closed": True,
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


def _exterior_probes(
    field: CurrentExteriorXRLeadingOscillatoryField,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    X = np.asarray(
        [field.X_h * math.exp(0.7), field.X_h * math.exp(2.3), field.X_R], dtype=float
    )
    eta = np.asarray((-0.21, 0.17, 0.33), dtype=float)
    t = np.asarray((0.37, 0.53, 0.67), dtype=float)
    theta = np.asarray((0.29, 0.83, 1.31), dtype=float)
    x, y, z, t = _cartesian_probe(field, X, eta, t, theta)
    return X, x, y, z, t


def materialize_through_xr_composite_receipt() -> dict[str, Any]:
    field = default_field()

    inner = _composition._evaluate_with_backend(
        field.leading_backend,
        _INNER_POINTS[:, 0],
        _INNER_POINTS[:, 1],
        _INNER_POINTS[:, 2],
        _INNER_TIMES,
    )
    X_ext, x_ext, y_ext, z_ext, t_ext = _exterior_probes(field)
    exterior = _composition._evaluate_with_backend(
        field.leading_backend, x_ext, y_ext, z_ext, t_ext
    )

    inner_closure = float(np.max(np.abs(inner.velocity - (inner.leading + inner.oscillatory))))
    exterior_closure = float(
        np.max(np.abs(exterior.velocity - (exterior.leading + exterior.oscillatory)))
    )
    inner_osc_rms = float(
        np.sqrt(np.mean(np.sum(inner.oscillatory * inner.oscillatory, axis=-1)))
    )
    exterior_osc_max = float(np.max(np.abs(exterior.oscillatory)))

    with __import__("tempfile").TemporaryDirectory() as tmp:
        path = Path(tmp) / "candidate.json"
        saved = field.save_candidate(path)
        loaded = CurrentExteriorXRLeadingOscillatoryField.load_candidate(path)
        replay_points = np.concatenate(
            (
                _INNER_POINTS,
                np.stack((x_ext, y_ext, z_ext), axis=-1),
            ),
            axis=0,
        )
        replay_times = np.concatenate((_INNER_TIMES, t_ext), axis=0)
        before = field.velocity(
            replay_points[:, 0], replay_points[:, 1], replay_points[:, 2], replay_times
        )
        after = loaded.velocity(
            replay_points[:, 0], replay_points[:, 1], replay_points[:, 2], replay_times
        )
        replay_max = float(np.max(np.abs(before - after)))
        saved_semantic = str(saved["semantic_sha256"])

    xb, yb, zb, tb = _cartesian_probe(
        field, field.X_R * 1.001, 0.2, 0.5, 0.3
    )
    beyond_fail_closed = False
    try:
        field.velocity(xb, yb, zb, tb)
    except (ValueError, RuntimeError):
        beyond_fail_closed = True

    return {
        "schema": RECEIPT_SCHEMA,
        "semantic_sha256": field.semantic_sha256,
        "saved_semantic_sha256": saved_semantic,
        "oscillatory_runtime_sha256": field.oscillatory_runtime_sha256,
        "agent1_semantic_sha256": field._leading_identity["semantic_sha256"],
        "agent1_configuration_sha256": field._leading_identity["configuration_sha256"],
        "public_contract": public_contract(),
        "diagnostic": {
            "inner_probe_count": int(_INNER_POINTS.shape[0]),
            "exterior_probe_count": int(X_ext.size),
            "exterior_X_over_Xh": [float(v / field.X_h) for v in X_ext],
            "exterior_X_over_XR": [float(v / field.X_R) for v in X_ext],
            "inner_composition_closure_max_abs": inner_closure,
            "exterior_composition_closure_max_abs": exterior_closure,
            "inner_oscillatory_vector_rms": inner_osc_rms,
            "exterior_oscillatory_abs_max": exterior_osc_max,
            "velocity_save_load_replay_max_abs": replay_max,
            "beyond_current_XR_fail_closed": beyond_fail_closed,
            "X_h": field.X_h,
            "X_R": field.X_R,
        },
        "truth_boundary": field.truth_boundary,
        "frozen_mechanical_gates": {
            "inner_composition_closure_max_abs": COMPOSITION_ATOL,
            "exterior_composition_closure_max_abs": COMPOSITION_ATOL,
            "inner_oscillatory_vector_rms_min": OSCILLATORY_SIGNAL_FLOOR,
            "exterior_oscillatory_abs_max": 0.0,
            "velocity_save_load_replay_max_abs": VELOCITY_REPLAY_ATOL,
            "beyond_current_XR_fail_closed": True,
        },
    }


def enforce_receipt(receipt: Mapping[str, Any]) -> None:
    if not isinstance(receipt, Mapping) or receipt.get("schema") != RECEIPT_SCHEMA:
        raise AssertionError("unexpected through-X_R composite receipt schema")
    diagnostic = receipt.get("diagnostic")
    gates = receipt.get("frozen_mechanical_gates")
    if not isinstance(diagnostic, Mapping) or not isinstance(gates, Mapping):
        raise AssertionError("through-X_R composite receipt is incomplete")

    if float(diagnostic.get("inner_composition_closure_max_abs", np.inf)) > COMPOSITION_ATOL:
        raise AssertionError("inner leading+oscillatory composition closure failed")
    if float(diagnostic.get("exterior_composition_closure_max_abs", np.inf)) > COMPOSITION_ATOL:
        raise AssertionError("exterior leading+oscillatory composition closure failed")
    if float(diagnostic.get("inner_oscillatory_vector_rms", 0.0)) < OSCILLATORY_SIGNAL_FLOOR:
        raise AssertionError("inner oscillatory increment is vacuous")
    if float(diagnostic.get("exterior_oscillatory_abs_max", np.inf)) != 0.0:
        raise AssertionError("frozen compact oscillation leaked into the far exterior probes")
    if float(diagnostic.get("velocity_save_load_replay_max_abs", np.inf)) > VELOCITY_REPLAY_ATOL:
        raise AssertionError("through-X_R candidate save/load velocity replay failed")
    if diagnostic.get("beyond_current_XR_fail_closed") is not True:
        raise AssertionError("through-X_R composite no longer fails closed beyond X_R")

    semantic = _require_sha256(receipt.get("semantic_sha256"), "receipt semantic SHA")
    saved_semantic = _require_sha256(
        receipt.get("saved_semantic_sha256"), "saved semantic SHA"
    )
    if semantic != saved_semantic:
        raise AssertionError("through-X_R semantic identity did not survive save/load")
    _require_sha256(receipt.get("oscillatory_runtime_sha256"), "oscillatory runtime SHA")
    _require_sha256(receipt.get("agent1_semantic_sha256"), "Agent-1 semantic SHA")
    _require_sha256(
        receipt.get("agent1_configuration_sha256"), "Agent-1 configuration SHA"
    )

    if gates != {
        "inner_composition_closure_max_abs": COMPOSITION_ATOL,
        "exterior_composition_closure_max_abs": COMPOSITION_ATOL,
        "inner_oscillatory_vector_rms_min": OSCILLATORY_SIGNAL_FLOOR,
        "exterior_oscillatory_abs_max": 0.0,
        "velocity_save_load_replay_max_abs": VELOCITY_REPLAY_ATOL,
        "beyond_current_XR_fail_closed": True,
    }:
        raise AssertionError("frozen through-X_R mechanical gates drifted")

    truth = receipt.get("truth_boundary")
    if truth != _truth_boundary():
        raise AssertionError("through-X_R composite truth boundary drifted or was illegally promoted")
    contract = receipt.get("public_contract")
    if contract != public_contract():
        raise AssertionError("through-X_R composite public contract drifted")
