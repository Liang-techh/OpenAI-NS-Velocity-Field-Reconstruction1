"""Identity-preserving save/load for the current partial Kokuno composite.

This Agent-2 increment does not change the velocity field.  It closes the narrow
representation gap identified after PR #970/#975: the partial Cartesian callable

    u_partial = u_lead,current-through-X_h + u_osc,frozen-complete-curl

already has an authenticated Agent-1 leading configuration, but its composite
semantic identity did not bind the *complete concrete* oscillatory runtime
realization exposed by ``KokunoPublicZPullbackCandidateOscillatoryVelocity.to_payload``.

The surface below therefore serializes the exact Agent-1 leading configuration and
semantic identity together with the complete oscillatory runtime payload and its
SHA-256 digest.  Reload reconstructs the leading backend through its existing
``from_configuration`` path and recomputes both leading and oscillatory identities.
Any runtime/payload drift fails closed before the reloaded field is returned.

This is a partial-domain delivery-identity improvement only.  It does not extend
velocity beyond the currently materialized ``X_h``, add a correction velocity,
pressure, forcing, or a complete Navier--Stokes residual.  The corrected 2026-09-09
Kokuno reconstruction remains structural provenance; the concrete oscillatory
runtime, public-z pullback, current joined leading profile and this serialization
contract remain repository realizations, not paper-exact hidden data.
"""
from __future__ import annotations

from dataclasses import dataclass, field as dataclass_field
from functools import lru_cache
import hashlib
import importlib
import inspect
import json
from pathlib import Path
import tempfile
from typing import Any, Mapping

import numpy as np

from . import kokuno_current_partial_leading_oscillatory_velocity as _composition
from .kokuno_public_z_pullback_velocity import default_field as _default_oscillatory_field

TASK = "K2-OSC-082"
SCHEMA = "kokuno-a2-current-partial-composite-identity-save-load-v1"
SERIALIZED_SCHEMA = "kokuno-a2-current-partial-composite-saved-candidate-v1"
RECEIPT_SCHEMA = "kokuno-a2-current-partial-composite-identity-save-load-receipt-v1"

PARENT_AGENT2_PR = 975
PARENT_AGENT2_HEAD = "b716dfe495763439ac51f6e866211337910017cb"
COMPOSITION_AGENT2_PR = 970
COMPOSITION_AGENT2_HEAD = "3a6405bbd3d10b8c3b38078f0c989c45e5d407b4"
AGENT1_PR = _composition.AGENT1_PR
AGENT1_HEAD = _composition.AGENT1_HEAD
AGENT1_MODULE = _composition.AGENT1_MODULE
AGENT1_CLASS = _composition.AGENT1_CLASS
AGENT1_SOURCE_BLOB_SHA1 = _composition.AGENT1_SOURCE_BLOB_SHA1
PUBLIC_Z_OSCILLATORY_SOURCE_BLOB_SHA1 = "4a3abc1a11f7e054f651dd7bace8ada6d5a6b653"
SOURCE_CORRECTED_READER_COMMIT = _composition.SOURCE_CORRECTED_READER_COMMIT
SOURCE_CORRECTED_READER_DATE = _composition.SOURCE_CORRECTED_READER_DATE

VELOCITY_REPLAY_ATOL = 2.0e-12
OSCILLATORY_SIGNAL_FLOOR = 1.0e-12

_RECEIPT_POINTS = np.asarray(
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
_RECEIPT_TIMES = np.asarray((0.31, 0.39, 0.47, 0.55, 0.63, 0.71), dtype=float)


def _canonical_json(payload: Mapping[str, Any]) -> str:
    if not isinstance(payload, Mapping):
        raise TypeError("identity payload must be a mapping")
    try:
        return json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        raise ValueError("identity payload must be finite and JSON serializable") from exc


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(
        c in "0123456789abcdef" for c in value
    )


def _require_sha256(value: Any, label: str) -> str:
    if not _is_sha256(value):
        raise ValueError(f"{label} must be one lowercase SHA-256 hex digest")
    return str(value)


def _oscillatory_payload(field: Any | None = None) -> dict[str, Any]:
    runtime = _default_oscillatory_field() if field is None else field
    to_payload = getattr(runtime, "to_payload", None)
    if not callable(to_payload):
        raise RuntimeError("frozen oscillatory runtime does not expose to_payload()")
    payload = to_payload()
    if not isinstance(payload, Mapping):
        raise RuntimeError("frozen oscillatory runtime payload is not a mapping")
    # Canonicalization is part of the contract: it rejects NaN/non-JSON runtime state.
    canonical = _canonical_json(payload)
    normalized = json.loads(canonical)
    truth = normalized.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise RuntimeError("oscillatory runtime payload lost its truth boundary")
    if truth.get("paper_exact") is not False or truth.get("pde_validated") is not False:
        raise RuntimeError("oscillatory runtime truth boundary was promoted unexpectedly")
    return normalized


def oscillatory_runtime_sha256(field: Any | None = None) -> str:
    """SHA-256 of the complete concrete public-z oscillatory runtime payload."""
    return _sha256(_oscillatory_payload(field))


def _leading_configuration_sha256(backend: Any) -> str:
    configuration = backend.configuration()
    if not isinstance(configuration, Mapping):
        raise RuntimeError("Agent-1 leading configuration is missing")
    return _sha256(configuration)


def _load_exact_agent1_from_configuration(configuration: Mapping[str, Any]) -> Any:
    try:
        module = importlib.import_module(AGENT1_MODULE)
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "exact Agent-1 #965 runtime is unavailable; cannot reload the partial composite"
        ) from exc
    cls = getattr(module, AGENT1_CLASS, None)
    if cls is None or not callable(getattr(cls, "from_configuration", None)):
        raise RuntimeError("exact Agent-1 #965 class/from_configuration is unavailable")
    backend = cls.from_configuration(configuration)
    _composition._validate_exact_agent1_backend(backend)
    return backend


def _truth_boundary() -> dict[str, bool]:
    return {
        "current_partial_leading_plus_oscillatory_velocity_materialized": True,
        "full_concrete_oscillatory_runtime_digest_bound": True,
        "identity_preserving_partial_composite_save_load_available": True,
        "reload_recomputes_leading_and_oscillatory_identities": True,
        "velocity_beyond_Xh_materialized": False,
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


@dataclass(frozen=True)
class CurrentPartialCompositeIdentityField:
    """Current partial ``leading + oscillatory`` field with bound runtime identity."""

    leading_backend: Any
    _parent_field: Any = dataclass_field(init=False, repr=False, compare=False)
    _leading_identity: dict[str, Any] = dataclass_field(init=False, repr=False, compare=False)
    _oscillatory_runtime_payload: dict[str, Any] = dataclass_field(
        init=False, repr=False, compare=False
    )
    _oscillatory_runtime_sha256: str = dataclass_field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        parent = _composition._CurrentPartialLeadingOscillatoryField(self.leading_backend)
        identity = _composition._validate_exact_agent1_backend(self.leading_backend)
        oscillatory_payload = _oscillatory_payload()
        object.__setattr__(self, "_parent_field", parent)
        object.__setattr__(self, "_leading_identity", dict(identity))
        object.__setattr__(self, "_oscillatory_runtime_payload", oscillatory_payload)
        object.__setattr__(self, "_oscillatory_runtime_sha256", _sha256(oscillatory_payload))

    @property
    def X_h(self) -> float:
        return float(self._parent_field.X_h)

    @property
    def oscillatory_runtime_sha256(self) -> str:
        return self._oscillatory_runtime_sha256

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return _truth_boundary()

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self._parent_field.velocity(x, y, z, t)

    def configuration(self) -> dict[str, Any]:
        leading_configuration = self.leading_backend.configuration()
        if not isinstance(leading_configuration, Mapping):
            raise RuntimeError("Agent-1 leading configuration is missing")
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent2": {
                "pr": PARENT_AGENT2_PR,
                "head": PARENT_AGENT2_HEAD,
            },
            "composition_parent": {
                "pr": COMPOSITION_AGENT2_PR,
                "head": COMPOSITION_AGENT2_HEAD,
                "semantic_sha256": self._parent_field.semantic_sha256,
            },
            "agent1_leading": {
                "pr": AGENT1_PR,
                "head": AGENT1_HEAD,
                "module": AGENT1_MODULE,
                "class": AGENT1_CLASS,
                "source_blob_sha1": AGENT1_SOURCE_BLOB_SHA1,
                "configuration": json.loads(_canonical_json(leading_configuration)),
                "configuration_sha256": _leading_configuration_sha256(self.leading_backend),
                "semantic_sha256": str(self._leading_identity["semantic_sha256"]),
            },
            "oscillatory_runtime": {
                "public_z_source_blob_sha1": PUBLIC_Z_OSCILLATORY_SOURCE_BLOB_SHA1,
                "payload": self._oscillatory_runtime_payload,
                "payload_sha256": self._oscillatory_runtime_sha256,
                "binding": "complete concrete default_field().to_payload() realization",
            },
            "source_provenance": {
                "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
                "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
                "scope": "localized waves / complete curls plus public leading-coordinate structure",
            },
            "truth_boundary": self.truth_boundary,
        }

    @property
    def semantic_sha256(self) -> str:
        # The semantic identity hashes the complete serialized configuration, including
        # the full oscillatory runtime payload, not only adapter/support metadata.
        return _sha256(self.configuration())

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
    def from_configuration(cls, payload: Mapping[str, Any]) -> "CurrentPartialCompositeIdentityField":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected partial-composite identity configuration schema")
        parent = payload.get("parent_agent2")
        if not isinstance(parent, Mapping) or parent.get("pr") != PARENT_AGENT2_PR or parent.get("head") != PARENT_AGENT2_HEAD:
            raise ValueError("partial-composite identity parent drifted")
        composition_parent = payload.get("composition_parent")
        if not isinstance(composition_parent, Mapping):
            raise ValueError("composition parent identity is missing")
        if composition_parent.get("pr") != COMPOSITION_AGENT2_PR or composition_parent.get("head") != COMPOSITION_AGENT2_HEAD:
            raise ValueError("composition parent identity drifted")
        stored_parent_semantic = _require_sha256(
            composition_parent.get("semantic_sha256"), "composition parent semantic SHA"
        )

        leading = payload.get("agent1_leading")
        if not isinstance(leading, Mapping):
            raise ValueError("serialized Agent-1 leading identity is missing")
        if (
            leading.get("pr") != AGENT1_PR
            or leading.get("head") != AGENT1_HEAD
            or leading.get("module") != AGENT1_MODULE
            or leading.get("class") != AGENT1_CLASS
            or leading.get("source_blob_sha1") != AGENT1_SOURCE_BLOB_SHA1
        ):
            raise ValueError("serialized Agent-1 leading provenance drifted")
        configuration = leading.get("configuration")
        if not isinstance(configuration, Mapping):
            raise ValueError("serialized Agent-1 leading configuration is missing")
        stored_leading_config_sha = _require_sha256(
            leading.get("configuration_sha256"), "Agent-1 configuration SHA"
        )
        if _sha256(configuration) != stored_leading_config_sha:
            raise ValueError("serialized Agent-1 leading configuration digest mismatch")
        stored_leading_semantic = _require_sha256(
            leading.get("semantic_sha256"), "Agent-1 semantic SHA"
        )
        backend = _load_exact_agent1_from_configuration(configuration)
        if _leading_configuration_sha256(backend) != stored_leading_config_sha:
            raise ValueError("reloaded Agent-1 leading configuration identity drifted")
        if str(getattr(backend, "semantic_sha256", "")) != stored_leading_semantic:
            raise ValueError("reloaded Agent-1 leading semantic identity drifted")

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
        if _sha256(stored_osc_payload) != stored_osc_sha:
            raise ValueError("serialized oscillatory runtime payload digest mismatch")
        current_osc_payload = _oscillatory_payload()
        if _sha256(current_osc_payload) != stored_osc_sha:
            raise ValueError("current concrete oscillatory runtime identity drifted")
        if _canonical_json(current_osc_payload) != _canonical_json(stored_osc_payload):
            raise ValueError("current concrete oscillatory runtime payload drifted")

        obj = cls(backend)
        if obj._parent_field.semantic_sha256 != stored_parent_semantic:
            raise ValueError("reloaded parent partial-composite semantic identity drifted")
        if obj.oscillatory_runtime_sha256 != stored_osc_sha:
            raise ValueError("reloaded oscillatory runtime SHA drifted")
        return obj

    @classmethod
    def load_candidate(cls, path: str | Path) -> "CurrentPartialCompositeIdentityField":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, Mapping) or raw.get("schema") != SERIALIZED_SCHEMA:
            raise ValueError("unexpected saved partial-composite candidate schema")
        configuration = raw.get("configuration")
        if not isinstance(configuration, Mapping):
            raise ValueError("saved partial-composite configuration is missing")
        stored_semantic = _require_sha256(raw.get("semantic_sha256"), "saved semantic SHA")
        obj = cls.from_configuration(configuration)
        if obj.semantic_sha256 != stored_semantic:
            raise ValueError("saved/reloaded partial-composite semantic identity drifted")
        return obj


@lru_cache(maxsize=1)
def default_field() -> CurrentPartialCompositeIdentityField:
    parent = _composition.default_field()
    return CurrentPartialCompositeIdentityField(parent.leading_backend)


def velocity(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    """Return the identity-bound current partial Cartesian composite velocity."""
    return default_field().velocity(x, y, z, t)


def public_contract() -> dict[str, Any]:
    params = inspect.signature(velocity).parameters
    forbidden = {
        "amplitude", "phase", "phase_offset", "scale", "orientation", "support",
        "spatial_step", "derivative_step", "residual", "defect", "target", "forcing",
        "pressure", "viscosity", "nu", "gain", "damping", "threshold", "mean",
        "stress", "inverse", "correction", "angular_order",
    }
    truth = _truth_boundary()
    return {
        "public_velocity_inputs": list(params),
        "forbidden_velocity_inputs_present": sorted(forbidden.intersection(params)),
        "full_oscillatory_runtime_payload_bound": True,
        "leading_configuration_and_semantic_bound": True,
        "identity_preserving_partial_composite_save_load": True,
        "reload_recomputes_runtime_digests": True,
        "field_formula_changed": False,
        "new_oscillatory_parameters": False,
        "mean_projection_performed": False,
        "radial_inverse_performed": False,
        "correction_velocity_constructed": False,
        "pressure_or_forcing_added": False,
        "complete_ns_residual": False,
        "velocity_export_ready": truth["velocity_export_ready"],
        "paper_exact": truth["paper_exact"],
        "pde_validated": truth["pde_validated"],
    }


def materialize_identity_save_load_receipt() -> dict[str, Any]:
    field = default_field()
    before = np.asarray(
        field.velocity(
            _RECEIPT_POINTS[:, 0],
            _RECEIPT_POINTS[:, 1],
            _RECEIPT_POINTS[:, 2],
            _RECEIPT_TIMES,
        ),
        dtype=float,
    )
    leading = np.asarray(
        field.leading_backend.velocity(
            _RECEIPT_POINTS[:, 0],
            _RECEIPT_POINTS[:, 1],
            _RECEIPT_POINTS[:, 2],
            _RECEIPT_TIMES,
        ),
        dtype=float,
    )
    if before.shape != leading.shape or before.shape != (_RECEIPT_POINTS.shape[0], 3):
        raise RuntimeError("identity receipt velocity shape drifted")

    with tempfile.TemporaryDirectory(prefix="kokuno-a2-identity-") as directory:
        path = Path(directory) / "candidate.json"
        saved = field.save_candidate(path)
        file_sha256 = hashlib.sha256(path.read_bytes()).hexdigest()
        loaded = CurrentPartialCompositeIdentityField.load_candidate(path)
        after = np.asarray(
            loaded.velocity(
                _RECEIPT_POINTS[:, 0],
                _RECEIPT_POINTS[:, 1],
                _RECEIPT_POINTS[:, 2],
                _RECEIPT_TIMES,
            ),
            dtype=float,
        )

    velocity_replay = float(np.max(np.abs(after - before), initial=0.0))
    oscillatory_rms = float(
        np.sqrt(np.mean(np.sum((before - leading) * (before - leading), axis=-1)))
    )
    beyond_fail_closed = False
    try:
        loaded.velocity(1.0e6, 0.0, 0.0, 0.5)
    except (ValueError, RuntimeError):
        beyond_fail_closed = True

    saved_config = saved["configuration"]
    return {
        "schema": RECEIPT_SCHEMA,
        "semantic_sha256": field.semantic_sha256,
        "saved_semantic_sha256": saved["semantic_sha256"],
        "reloaded_semantic_sha256": loaded.semantic_sha256,
        "saved_file_sha256": file_sha256,
        "leading_semantic_sha256": str(field._leading_identity["semantic_sha256"]),
        "oscillatory_runtime_sha256": field.oscillatory_runtime_sha256,
        "saved_oscillatory_runtime_sha256": saved_config["oscillatory_runtime"]["payload_sha256"],
        "reloaded_oscillatory_runtime_sha256": loaded.oscillatory_runtime_sha256,
        "parent_composite_semantic_sha256": field._parent_field.semantic_sha256,
        "reloaded_parent_composite_semantic_sha256": loaded._parent_field.semantic_sha256,
        "public_contract": public_contract(),
        "diagnostic": {
            "probe_count": int(_RECEIPT_POINTS.shape[0]),
            "velocity_save_load_replay_max_abs": velocity_replay,
            "oscillatory_increment_vector_rms": oscillatory_rms,
            "beyond_current_Xh_fail_closed_after_reload": beyond_fail_closed,
        },
        "truth_boundary": field.truth_boundary,
        "frozen_mechanical_gates": {
            "velocity_save_load_replay_max_abs": VELOCITY_REPLAY_ATOL,
            "oscillatory_increment_vector_rms_min": OSCILLATORY_SIGNAL_FLOOR,
            "semantic_identity_exact_replay": True,
            "oscillatory_runtime_identity_exact_replay": True,
            "parent_composite_identity_exact_replay": True,
            "beyond_current_Xh_fail_closed_after_reload": True,
        },
    }


def enforce_receipt(receipt: Mapping[str, Any]) -> None:
    if not isinstance(receipt, Mapping) or receipt.get("schema") != RECEIPT_SCHEMA:
        raise ValueError("unexpected partial-composite identity receipt schema")
    gates = receipt.get("frozen_mechanical_gates")
    diagnostic = receipt.get("diagnostic")
    truth = receipt.get("truth_boundary")
    if not isinstance(gates, Mapping) or not isinstance(diagnostic, Mapping) or not isinstance(truth, Mapping):
        raise ValueError("identity receipt is incomplete")
    semantic = _require_sha256(receipt.get("semantic_sha256"), "receipt semantic SHA")
    if receipt.get("saved_semantic_sha256") != semantic or receipt.get("reloaded_semantic_sha256") != semantic:
        raise AssertionError("partial-composite semantic identity did not replay exactly")
    osc = _require_sha256(receipt.get("oscillatory_runtime_sha256"), "receipt oscillatory SHA")
    if receipt.get("saved_oscillatory_runtime_sha256") != osc or receipt.get("reloaded_oscillatory_runtime_sha256") != osc:
        raise AssertionError("complete oscillatory runtime identity did not replay exactly")
    parent = _require_sha256(receipt.get("parent_composite_semantic_sha256"), "receipt parent composite SHA")
    if receipt.get("reloaded_parent_composite_semantic_sha256") != parent:
        raise AssertionError("parent partial-composite identity did not replay exactly")
    if float(diagnostic.get("velocity_save_load_replay_max_abs", np.inf)) > float(gates["velocity_save_load_replay_max_abs"]):
        raise AssertionError("saved/reloaded velocity replay exceeded frozen gate")
    if float(diagnostic.get("oscillatory_increment_vector_rms", 0.0)) < float(gates["oscillatory_increment_vector_rms_min"]):
        raise AssertionError("oscillatory increment became vacuous")
    if diagnostic.get("beyond_current_Xh_fail_closed_after_reload") is not True:
        raise AssertionError("reloaded partial candidate no longer fails closed beyond X_h")
    for key in (
        "full_concrete_oscillatory_runtime_digest_bound",
        "identity_preserving_partial_composite_save_load_available",
        "reload_recomputes_leading_and_oscillatory_identities",
    ):
        if truth.get(key) is not True:
            raise AssertionError(f"identity truth boundary lost required true state: {key}")
    for key in (
        "velocity_beyond_Xh_materialized",
        "outer_global_leading_velocity_materialized",
        "global_compact_support_completed",
        "agent3_correction_velocity_composed",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_velocity_pressure_forcing_api",
        "heldout_ns_residual_assessed",
        "same_protocol_comparable_to_st006",
        "residual_reduction_claimed",
        "velocity_export_ready",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "pde_validated",
    ):
        if truth.get(key) is not False:
            raise AssertionError(f"identity-only increment illegally promoted: {key}")
