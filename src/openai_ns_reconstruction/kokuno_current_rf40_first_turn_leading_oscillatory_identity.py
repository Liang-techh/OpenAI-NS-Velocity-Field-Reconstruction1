"""Identity-bound current Kokuno leading + complete-curl oscillation through RF40 first turn.

Kokuno Agent 1 PR #986 advances the current Cartesian leading field one stage past
``X_R`` through the public first RF40 turn ``X_R < X <= X_1=e X_R``.  Agent 2
PR #987 already composes the preceding Agent-1 #980 field through ``X_R`` with
the frozen complete-curl oscillation and binds the full concrete oscillatory
runtime to deterministic save/load.

This module makes exactly one further composition step:

    u_current_to_X1(x,t) = u_lead,current-through-X1(x,t)
                           + u_osc,frozen-complete-curl(x,t).

Neither summand is changed.  The Agent-1 implementation is loaded from the exact
#986 class/blob, while the oscillatory field remains the existing axis-safe,
strict-support A2 runtime.  The complete oscillatory ``to_payload()`` realization
continues to participate in the semantic identity and reload firewall.

The corrected 2026-09-09 Kokuno reconstruction is structural provenance for the
localized-wave/complete-curl organization and the public RF40 first-turn
schedule.  The concrete oscillatory modes/support/public-z pullback, current
candidate-side leading realization, floating-point incompressibility memory and
this composition remain repository realizations.  This object is therefore not
paper-exact, not a complete post-RF40/global field, not pressure/forcing complete,
and not a Navier--Stokes validation.
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

from . import kokuno_current_exterior_xr_leading_oscillatory_identity as _xr_parent
from . import kokuno_current_partial_composite_identity_save_load as _identity
from . import kokuno_current_partial_leading_oscillatory_velocity as _composition

TASK = "K2-OSC-084"
SCHEMA = "kokuno-a2-current-rf40-first-turn-leading-oscillatory-identity-v1"
SERIALIZED_SCHEMA = "kokuno-a2-current-rf40-first-turn-leading-oscillatory-saved-candidate-v1"
RECEIPT_SCHEMA = "kokuno-a2-current-rf40-first-turn-leading-oscillatory-receipt-v1"

PARENT_AGENT2_PR = 987
PARENT_AGENT2_HEAD = "4be2c9ee898c44dd1ad2217e90601b161fe81964"
PARENT_AGENT2_MODULE = (
    "openai_ns_reconstruction.kokuno_current_exterior_xr_leading_oscillatory_identity"
)
PARENT_AGENT2_SOURCE_BLOB_SHA1 = "93a3da5c3aed4be8310cc510edfab6af940b64dc"

AGENT1_PR = 986
AGENT1_HEAD = "23c00b98526e187ff04d432745300a084ec859f2"
AGENT1_MODULE = "openai_ns_reconstruction.kokuno_pa16_current_cartesian_rf40_first_turn"
AGENT1_CLASS = "KokunoPA16CurrentCartesianRF40FirstTurn"
AGENT1_SOURCE_BLOB_SHA1 = "bf71948cbdb827d791b86ab8f00496bd1b231961"
AGENT1_TEST_BLOB_SHA1 = "14703a4151ef86fecc2a3731f7e0bd3f81602ad2"
AGENT1_WORKFLOW_BLOB_SHA1 = "c723e37ec024e5bbde4a4c05e83401b6d0e9c6ab"

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


def _validate_parent_a2_source() -> str:
    cls = _xr_parent.CurrentExteriorXRLeadingOscillatoryField
    if cls.__module__ != PARENT_AGENT2_MODULE:
        raise RuntimeError("Agent-2 #987 parent module identity drifted")
    source_path = inspect.getsourcefile(cls)
    if source_path is None:
        raise RuntimeError("cannot locate Agent-2 #987 parent implementation source")
    blob = _composition._git_blob_sha1(source_path)
    if blob != PARENT_AGENT2_SOURCE_BLOB_SHA1:
        raise RuntimeError("Agent-2 #987 parent implementation source blob drifted")
    return blob


def _required_agent1_truth() -> dict[str, bool]:
    return {
        "public_reconstruction_source": True,
        "current_cartesian_leading_through_XR_consumed": True,
        "public_RF40_first_turn_formula_reused": True,
        "current_lineage_RF40_first_turn_materialized": True,
        "current_cartesian_spacetime_leading_velocity_through_RF40_first_turn_materialized": True,
        "current_incompressibility_memory_carried_through_RF40_first_turn": True,
        "axis_regular_cartesian_formula_reused": True,
        "velocity_interface_vectorized": True,
        "velocity_configuration_serializable": True,
        "source_hidden_numeric_choices_recovered": False,
        "source_exact_current_lineage_certified": False,
        "full_post_XR_RF40_current_lineage_materialized": False,
        "RF40_axial_shutdown_current_lineage_materialized": False,
        "RF40_lambda_turn_current_lineage_materialized": False,
        "RF40_power_law_current_lineage_materialized": False,
        "current_lineage_cone_I1_I2_I3_I4_outer_overlays_completed": False,
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
        raise RuntimeError("loaded leading backend is not the pinned Agent-1 #986 class")
    source_path = inspect.getsourcefile(cls)
    if source_path is None:
        raise RuntimeError("cannot locate Agent-1 #986 leading implementation source")
    blob = _composition._git_blob_sha1(source_path)
    if blob != AGENT1_SOURCE_BLOB_SHA1:
        raise RuntimeError("Agent-1 #986 leading implementation source blob drifted")

    truth = getattr(backend, "truth_boundary", None)
    if not isinstance(truth, Mapping):
        raise RuntimeError("Agent-1 #986 truth boundary is missing")
    for key, expected in _required_agent1_truth().items():
        if truth.get(key) is not expected:
            raise RuntimeError(f"Agent-1 #986 truth boundary drifted at {key}")

    semantic = str(getattr(backend, "semantic_sha256", ""))
    _require_sha256(semantic, "Agent-1 #986 semantic SHA")
    configuration = backend.configuration()
    if not isinstance(configuration, Mapping):
        raise RuntimeError("Agent-1 #986 configuration is missing")
    config_sha = _identity._sha256(configuration)

    X_h = float(getattr(getattr(backend, "current", None), "X_h", np.nan))
    X_R = float(getattr(backend, "X_R", np.nan))
    X_1 = float(getattr(backend, "X_1", np.nan))
    D = float(getattr(backend, "D", np.nan))
    if not all(np.isfinite(v) for v in (X_h, X_R, X_1, D)):
        raise RuntimeError("Agent-1 #986 X_h/X_R/X_1/D geometry is invalid")
    if not (0.0 < X_h < X_R < X_1):
        raise RuntimeError("Agent-1 #986 staged radial geometry is invalid")
    if not math.isclose(X_h / X_R, math.exp(-5.0), rel_tol=2.0e-14, abs_tol=0.0):
        raise RuntimeError("Agent-1 #986 X_h/X_R identity drifted")
    if not math.isclose(X_1 / X_R, math.e, rel_tol=2.0e-14, abs_tol=0.0):
        raise RuntimeError("Agent-1 #986 RF40 first-turn endpoint identity drifted")

    return {
        "semantic_sha256": semantic,
        "configuration_sha256": config_sha,
        "source_blob_sha1": blob,
        "X_h": X_h,
        "X_R": X_R,
        "X_1": X_1,
        "D": D,
    }


@lru_cache(maxsize=1)
def _load_exact_agent1_backend() -> Any:
    try:
        module = importlib.import_module(AGENT1_MODULE)
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "exact Agent-1 #986 runtime is unavailable; cannot materialize the RF40-first-turn composite"
        ) from exc
    cls = getattr(module, AGENT1_CLASS, None)
    if cls is None:
        raise RuntimeError("exact Agent-1 #986 class is unavailable")
    backend = cls()
    _validate_exact_agent1_backend(backend)
    return backend


def _load_exact_agent1_from_configuration(configuration: Mapping[str, Any]) -> Any:
    try:
        module = importlib.import_module(AGENT1_MODULE)
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "exact Agent-1 #986 runtime is unavailable; cannot reload the RF40-first-turn composite"
        ) from exc
    cls = getattr(module, AGENT1_CLASS, None)
    if cls is None or not callable(getattr(cls, "from_configuration", None)):
        raise RuntimeError("exact Agent-1 #986 class/from_configuration is unavailable")
    backend = cls.from_configuration(configuration)
    _validate_exact_agent1_backend(backend)
    return backend


def _truth_boundary() -> dict[str, bool]:
    return {
        "current_leading_plus_oscillatory_velocity_through_RF40_first_turn_materialized": True,
        "post_XR_RF40_first_turn_leading_plus_oscillatory_materialized": True,
        "full_concrete_oscillatory_runtime_digest_bound": True,
        "identity_preserving_RF40_first_turn_composite_save_load_available": True,
        "reload_recomputes_leading_and_oscillatory_identities": True,
        "velocity_after_RF40_first_turn_materialized": False,
        "full_post_XR_RF40_current_lineage_materialized": False,
        "RF40_axial_shutdown_current_lineage_materialized": False,
        "RF40_lambda_turn_current_lineage_materialized": False,
        "RF40_power_law_current_lineage_materialized": False,
        "current_lineage_cone_I1_I2_I3_I4_outer_overlays_completed": False,
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
    field: "CurrentRF40FirstTurnLeadingOscillatoryField",
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
class CurrentRF40FirstTurnLeadingOscillatoryField:
    """Exact #986 current leading through ``X_1=e X_R`` plus frozen A2 oscillation."""

    leading_backend: Any
    _leading_identity: dict[str, Any] = dataclass_field(init=False, repr=False, compare=False)
    _oscillatory_runtime_payload: dict[str, Any] = dataclass_field(
        init=False, repr=False, compare=False
    )
    _oscillatory_runtime_sha256: str = dataclass_field(init=False, repr=False, compare=False)
    _parent_agent2_source_blob_sha1: str = dataclass_field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        parent_blob = _validate_parent_a2_source()
        identity = _validate_exact_agent1_backend(self.leading_backend)
        oscillatory_payload = _identity._oscillatory_payload()
        object.__setattr__(self, "_parent_agent2_source_blob_sha1", parent_blob)
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
    def X_1(self) -> float:
        return float(self._leading_identity["X_1"])

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
            raise RuntimeError("Agent-1 #986 leading configuration is missing")
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
                "configuration": json.loads(_identity._canonical_json(leading_configuration)),
                "configuration_sha256": self._leading_identity["configuration_sha256"],
                "semantic_sha256": self._leading_identity["semantic_sha256"],
                "X_h": self.X_h,
                "X_R": self.X_R,
                "X_1": self.X_1,
            },
            "oscillatory_runtime": {
                "public_z_source_blob_sha1": PUBLIC_Z_OSCILLATORY_SOURCE_BLOB_SHA1,
                "payload": self._oscillatory_runtime_payload,
                "payload_sha256": self._oscillatory_runtime_sha256,
                "binding": "complete concrete default_field().to_payload() realization",
            },
            "composition": (
                "u_current_to_X1 = u_lead_current_through_RF40_first_turn "
                "+ u_osc_frozen_complete_curl"
            ),
            "source_provenance": {
                "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
                "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
                "scope": (
                    "localized waves / complete curls plus public leading-coordinate "
                    "and RF40 first-turn structure"
                ),
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
    ) -> "CurrentRF40FirstTurnLeadingOscillatoryField":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected RF40-first-turn composite configuration schema")

        parent = payload.get("parent_agent2")
        if not isinstance(parent, Mapping):
            raise ValueError("serialized Agent-2 #987 parent identity is missing")
        if (
            parent.get("pr") != PARENT_AGENT2_PR
            or parent.get("head") != PARENT_AGENT2_HEAD
            or parent.get("module") != PARENT_AGENT2_MODULE
            or parent.get("source_blob_sha1") != PARENT_AGENT2_SOURCE_BLOB_SHA1
        ):
            raise ValueError("RF40-first-turn composite parent identity drifted")
        if _validate_parent_a2_source() != PARENT_AGENT2_SOURCE_BLOB_SHA1:
            raise ValueError("current Agent-2 #987 parent source identity drifted")

        leading = payload.get("agent1_leading")
        if not isinstance(leading, Mapping):
            raise ValueError("serialized Agent-1 #986 identity is missing")
        if (
            leading.get("pr") != AGENT1_PR
            or leading.get("head") != AGENT1_HEAD
            or leading.get("module") != AGENT1_MODULE
            or leading.get("class") != AGENT1_CLASS
            or leading.get("source_blob_sha1") != AGENT1_SOURCE_BLOB_SHA1
        ):
            raise ValueError("serialized Agent-1 #986 provenance drifted")
        leading_configuration = leading.get("configuration")
        if not isinstance(leading_configuration, Mapping):
            raise ValueError("serialized Agent-1 #986 configuration is missing")
        stored_config_sha = _require_sha256(
            leading.get("configuration_sha256"), "Agent-1 #986 configuration SHA"
        )
        if _identity._sha256(leading_configuration) != stored_config_sha:
            raise ValueError("serialized Agent-1 #986 configuration digest mismatch")
        stored_semantic = _require_sha256(
            leading.get("semantic_sha256"), "Agent-1 #986 semantic SHA"
        )
        backend = _load_exact_agent1_from_configuration(leading_configuration)
        identity = _validate_exact_agent1_backend(backend)
        if identity["configuration_sha256"] != stored_config_sha:
            raise ValueError("reloaded Agent-1 #986 configuration identity drifted")
        if identity["semantic_sha256"] != stored_semantic:
            raise ValueError("reloaded Agent-1 #986 semantic identity drifted")
        for key in ("X_h", "X_R", "X_1"):
            if not math.isclose(
                float(leading.get(key, np.nan)), float(identity[key]), rel_tol=0.0, abs_tol=0.0
            ):
                raise ValueError(f"reloaded Agent-1 #986 {key} identity drifted")

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
    def load_candidate(cls, path: str | Path) -> "CurrentRF40FirstTurnLeadingOscillatoryField":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, Mapping) or raw.get("schema") != SERIALIZED_SCHEMA:
            raise ValueError("unexpected saved RF40-first-turn composite candidate schema")
        configuration = raw.get("configuration")
        if not isinstance(configuration, Mapping):
            raise ValueError("saved RF40-first-turn composite configuration is missing")
        stored_semantic = _require_sha256(
            raw.get("semantic_sha256"), "saved RF40-first-turn composite semantic SHA"
        )
        if _identity._sha256(configuration) != stored_semantic:
            raise ValueError("saved RF40-first-turn composite semantic digest mismatch")
        obj = cls.from_configuration(configuration)
        if obj.semantic_sha256 != stored_semantic:
            raise ValueError("reloaded RF40-first-turn composite semantic identity drifted")
        return obj


@lru_cache(maxsize=1)
def default_field() -> CurrentRF40FirstTurnLeadingOscillatoryField:
    return CurrentRF40FirstTurnLeadingOscillatoryField(_load_exact_agent1_backend())


def velocity(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    """Return current leading through RF40 first turn plus frozen A2 oscillation."""
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
        "leading_materialized_through_RF40_first_turn": True,
        "full_oscillatory_runtime_payload_bound": True,
        "identity_preserving_save_load": True,
        "after_RF40_first_turn_fail_closed": True,
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


def _post_xr_probes(
    field: CurrentRF40FirstTurnLeadingOscillatoryField,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    X = np.asarray(
        [
            field.X_R * math.exp(0.18),
            field.X_R * math.exp(0.57),
            field.X_1,
        ],
        dtype=float,
    )
    eta = np.asarray((-0.23, 0.19, 0.31), dtype=float)
    t = np.asarray((0.37, 0.53, 0.67), dtype=float)
    theta = np.asarray((0.31, 0.89, 1.37), dtype=float)
    x, y, z, t = _cartesian_probe(field, X, eta, t, theta)
    return X, x, y, z, t


def materialize_rf40_first_turn_composite_receipt() -> dict[str, Any]:
    field = default_field()

    inner = _composition._evaluate_with_backend(
        field.leading_backend,
        _INNER_POINTS[:, 0],
        _INNER_POINTS[:, 1],
        _INNER_POINTS[:, 2],
        _INNER_TIMES,
    )
    X_post, x_post, y_post, z_post, t_post = _post_xr_probes(field)
    post = _composition._evaluate_with_backend(
        field.leading_backend, x_post, y_post, z_post, t_post
    )

    inner_closure = float(np.max(np.abs(inner.velocity - (inner.leading + inner.oscillatory))))
    post_closure = float(np.max(np.abs(post.velocity - (post.leading + post.oscillatory))))
    inner_osc_rms = float(
        np.sqrt(np.mean(np.sum(inner.oscillatory * inner.oscillatory, axis=-1)))
    )
    post_osc_max = float(np.max(np.abs(post.oscillatory)))

    with __import__("tempfile").TemporaryDirectory() as tmp:
        path = Path(tmp) / "candidate.json"
        saved = field.save_candidate(path)
        loaded = CurrentRF40FirstTurnLeadingOscillatoryField.load_candidate(path)
        replay_points = np.concatenate(
            (
                _INNER_POINTS,
                np.stack((x_post, y_post, z_post), axis=-1),
            ),
            axis=0,
        )
        replay_times = np.concatenate((_INNER_TIMES, t_post), axis=0)
        before = field.velocity(
            replay_points[:, 0], replay_points[:, 1], replay_points[:, 2], replay_times
        )
        after = loaded.velocity(
            replay_points[:, 0], replay_points[:, 1], replay_points[:, 2], replay_times
        )
        replay_max = float(np.max(np.abs(before - after)))
        saved_semantic = str(saved["semantic_sha256"])

    xb, yb, zb, tb = _cartesian_probe(field, field.X_1 * 1.001, 0.2, 0.5, 0.3)
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
        "parent_agent2_source_blob_sha1": field._parent_agent2_source_blob_sha1,
        "public_contract": public_contract(),
        "diagnostic": {
            "inner_probe_count": int(_INNER_POINTS.shape[0]),
            "post_XR_probe_count": int(X_post.size),
            "post_XR_X_over_XR": [float(v / field.X_R) for v in X_post],
            "post_XR_X_over_X1": [float(v / field.X_1) for v in X_post],
            "inner_composition_closure_max_abs": inner_closure,
            "post_XR_composition_closure_max_abs": post_closure,
            "inner_oscillatory_vector_rms": inner_osc_rms,
            "post_XR_oscillatory_abs_max": post_osc_max,
            "velocity_save_load_replay_max_abs": replay_max,
            "beyond_current_X1_fail_closed": beyond_fail_closed,
            "X_h": field.X_h,
            "X_R": field.X_R,
            "X_1": field.X_1,
        },
        "truth_boundary": field.truth_boundary,
        "frozen_mechanical_gates": {
            "inner_composition_closure_max_abs": COMPOSITION_ATOL,
            "post_XR_composition_closure_max_abs": COMPOSITION_ATOL,
            "inner_oscillatory_vector_rms_min": OSCILLATORY_SIGNAL_FLOOR,
            "post_XR_oscillatory_abs_max": 0.0,
            "velocity_save_load_replay_max_abs": VELOCITY_REPLAY_ATOL,
            "beyond_current_X1_fail_closed": True,
        },
    }


def enforce_receipt(receipt: Mapping[str, Any]) -> None:
    if not isinstance(receipt, Mapping) or receipt.get("schema") != RECEIPT_SCHEMA:
        raise AssertionError("unexpected RF40-first-turn composite receipt schema")
    diagnostic = receipt.get("diagnostic")
    gates = receipt.get("frozen_mechanical_gates")
    if not isinstance(diagnostic, Mapping) or not isinstance(gates, Mapping):
        raise AssertionError("RF40-first-turn composite receipt is incomplete")

    if float(diagnostic.get("inner_composition_closure_max_abs", np.inf)) > COMPOSITION_ATOL:
        raise AssertionError("inner leading+oscillatory composition closure failed")
    if float(diagnostic.get("post_XR_composition_closure_max_abs", np.inf)) > COMPOSITION_ATOL:
        raise AssertionError("post-X_R leading+oscillatory composition closure failed")
    if float(diagnostic.get("inner_oscillatory_vector_rms", 0.0)) < OSCILLATORY_SIGNAL_FLOOR:
        raise AssertionError("inner oscillatory increment is vacuous")
    if float(diagnostic.get("post_XR_oscillatory_abs_max", np.inf)) != 0.0:
        raise AssertionError("frozen compact oscillation leaked into post-X_R RF40 probes")
    if float(diagnostic.get("velocity_save_load_replay_max_abs", np.inf)) > VELOCITY_REPLAY_ATOL:
        raise AssertionError("RF40-first-turn candidate save/load velocity replay failed")
    if diagnostic.get("beyond_current_X1_fail_closed") is not True:
        raise AssertionError("RF40-first-turn composite no longer fails closed after X_1")

    semantic = _require_sha256(receipt.get("semantic_sha256"), "receipt semantic SHA")
    saved_semantic = _require_sha256(
        receipt.get("saved_semantic_sha256"), "saved semantic SHA"
    )
    if semantic != saved_semantic:
        raise AssertionError("RF40-first-turn semantic identity did not survive save/load")
    _require_sha256(receipt.get("oscillatory_runtime_sha256"), "oscillatory runtime SHA")
    _require_sha256(receipt.get("agent1_semantic_sha256"), "Agent-1 semantic SHA")
    _require_sha256(
        receipt.get("agent1_configuration_sha256"), "Agent-1 configuration SHA"
    )
    if receipt.get("parent_agent2_source_blob_sha1") != PARENT_AGENT2_SOURCE_BLOB_SHA1:
        raise AssertionError("Agent-2 #987 parent source identity drifted")

    if gates != {
        "inner_composition_closure_max_abs": COMPOSITION_ATOL,
        "post_XR_composition_closure_max_abs": COMPOSITION_ATOL,
        "inner_oscillatory_vector_rms_min": OSCILLATORY_SIGNAL_FLOOR,
        "post_XR_oscillatory_abs_max": 0.0,
        "velocity_save_load_replay_max_abs": VELOCITY_REPLAY_ATOL,
        "beyond_current_X1_fail_closed": True,
    }:
        raise AssertionError("frozen RF40-first-turn mechanical gates drifted")

    truth = receipt.get("truth_boundary")
    if truth != _truth_boundary():
        raise AssertionError("RF40-first-turn truth boundary drifted or was illegally promoted")
    contract = receipt.get("public_contract")
    if contract != public_contract():
        raise AssertionError("RF40-first-turn public contract drifted")
