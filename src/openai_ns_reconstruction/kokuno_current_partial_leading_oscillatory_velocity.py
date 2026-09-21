"""Partial current-lineage Kokuno leading + frozen complete-curl oscillation.

This Agent-2 increment consumes, but does not reimplement, the first current-lineage
Agent-1 Cartesian leading velocity through the currently materialized ``X_h`` and
adds the already-frozen Agent-2 oscillatory correction:

    u_partial(x,t) = u_lead_current(x,t) + u_osc(x,t).

The corrected 2026-09-09 Kokuno reconstruction is structural provenance for the
localized-wave / complete-curl organization.  The current joined profile, public-z
pullback, bounded oscillatory parameters, support mask, and this composition are
repository realizations.  This module is therefore explicitly partial-domain and
not paper-exact, globally export-ready, pressure/forcing complete, or PDE-validated.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
import importlib
import inspect
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_oscillatory_batch_axis_safety import (
    batch_axis_safety_sha256,
    velocity_osc_batch,
)
from .kokuno_oscillatory_batch_differentials import batch_differentials_sha256

TASK = "K2-OSC-080"
SCHEMA = "kokuno-a2-current-partial-leading-oscillatory-v1"
RECEIPT_SCHEMA = "kokuno-a2-current-partial-leading-oscillatory-receipt-v1"

PARENT_AGENT2_PR = 960
PARENT_AGENT2_HEAD = "6d2fb1f701a34f783dca15a267ae2ce0734ba741"
AGENT1_PR = 965
AGENT1_HEAD = "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1"
AGENT1_MODULE = "openai_ns_reconstruction.kokuno_pa16_current_cartesian_leading_velocity"
AGENT1_CLASS = "KokunoPA16CurrentCartesianLeadingVelocity"
AGENT1_SOURCE_BLOB_SHA1 = "ca8b80b0451be1a8f31deaf620f542d8c2e92c0c"

SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"

COMPOSITION_ATOL = 1.0e-13
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
_AXIS_POINTS = np.asarray(((0.0, 0.0, -0.10), (0.0, 0.0, 0.0), (0.0, 0.0, 0.12)))
_AXIS_TIMES = np.asarray((0.37, 0.51, 0.67), dtype=float)


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _git_blob_sha1(path: str | Path) -> str:
    raw = Path(path).read_bytes()
    header = f"blob {len(raw)}\0".encode("ascii")
    return hashlib.sha1(header + raw).hexdigest()


def _broadcast_coordinates(
    x: Any, y: Any, z: Any, t: Any
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    try:
        xx, yy, zz, tt = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
    except ValueError as exc:
        raise ValueError("x, y, z, t must be mutually broadcastable") from exc
    if not all(np.all(np.isfinite(a)) for a in (xx, yy, zz, tt)):
        raise ValueError("x, y, z, t must be finite")
    return xx, yy, zz, tt


@dataclass(frozen=True)
class _CompositionEvaluation:
    leading: np.ndarray
    oscillatory: np.ndarray
    velocity: np.ndarray


def _evaluate_with_backend(
    backend: Any, x: Any, y: Any, z: Any, t: Any
) -> _CompositionEvaluation:
    """Compose one supplied leading backend with the frozen A2 oscillation.

    This helper is private so tests can verify addition/fail-closed semantics without
    pretending that an arbitrary backend is the authenticated Agent-1 artifact.
    """
    xx, yy, zz, tt = _broadcast_coordinates(x, y, z, t)
    points = np.stack((xx, yy, zz), axis=-1)

    oscillatory = np.asarray(velocity_osc_batch(points, tt), dtype=float)
    expected = xx.shape + (3,)
    if oscillatory.shape != expected or not np.all(np.isfinite(oscillatory)):
        raise RuntimeError("frozen oscillatory backend returned an invalid velocity")

    leading = np.asarray(backend.velocity(xx, yy, zz, tt), dtype=float)
    if leading.shape != expected or not np.all(np.isfinite(leading)):
        raise RuntimeError("current leading backend returned an invalid velocity")

    total = leading + oscillatory
    if not np.all(np.isfinite(total)):
        raise RuntimeError("leading+oscillatory composition became non-finite")
    return _CompositionEvaluation(leading=leading, oscillatory=oscillatory, velocity=total)


def _required_agent1_truth() -> dict[str, bool]:
    return {
        "current_joined_cartesian_spacetime_leading_velocity_through_Xh_materialized": True,
        "outer_global_leading_velocity_materialized": False,
        "velocity_beyond_Xh_materialized": False,
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
        raise RuntimeError("loaded leading backend is not the pinned Agent-1 class")

    source_path = inspect.getsourcefile(cls)
    if source_path is None:
        raise RuntimeError("cannot locate Agent-1 leading implementation source")
    blob = _git_blob_sha1(source_path)
    if blob != AGENT1_SOURCE_BLOB_SHA1:
        raise RuntimeError("Agent-1 leading implementation source blob drifted")

    semantic = str(getattr(backend, "semantic_sha256", ""))
    if len(semantic) != 64 or any(ch not in "0123456789abcdef" for ch in semantic):
        raise RuntimeError("Agent-1 leading semantic identity is malformed")

    truth = getattr(backend, "truth_boundary", None)
    if not isinstance(truth, Mapping):
        raise RuntimeError("Agent-1 leading truth boundary is missing")
    for key, expected in _required_agent1_truth().items():
        if truth.get(key) is not expected:
            raise RuntimeError(f"Agent-1 leading truth boundary drifted at {key}")

    x_h = float(getattr(backend, "X_h", np.nan))
    if not np.isfinite(x_h) or x_h <= 0.0:
        raise RuntimeError("Agent-1 current X_h is invalid")

    config = backend.configuration()
    if not isinstance(config, Mapping):
        raise RuntimeError("Agent-1 leading configuration is missing")
    return {
        "semantic_sha256": semantic,
        "source_blob_sha1": blob,
        "X_h": x_h,
        "configuration_sha256": hashlib.sha256(
            json.dumps(config, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
        ).hexdigest(),
    }


@lru_cache(maxsize=1)
def _load_exact_agent1_backend() -> Any:
    try:
        module = importlib.import_module(AGENT1_MODULE)
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "exact Agent-1 #965 runtime is not present; this stacked A2 artifact "
            "requires the sibling #965 module to materialize the partial composition"
        ) from exc
    cls = getattr(module, AGENT1_CLASS, None)
    if cls is None:
        raise RuntimeError("exact Agent-1 #965 class is unavailable")
    backend = cls()
    _validate_exact_agent1_backend(backend)
    return backend


@dataclass(frozen=True)
class _CurrentPartialLeadingOscillatoryField:
    leading_backend: Any

    def __post_init__(self) -> None:
        identity = _validate_exact_agent1_backend(self.leading_backend)
        object.__setattr__(self, "_leading_identity", identity)

    @property
    def X_h(self) -> float:
        return float(self._leading_identity["X_h"])

    @property
    def semantic_sha256(self) -> str:
        return _sha256(self.semantic_payload())

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return {
            "current_partial_cartesian_leading_materialized_through_Xh": True,
            "frozen_complete_curl_oscillation_composed": True,
            "leading_plus_oscillatory_partial_velocity_materialized": True,
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
            "paper_exact": False,
            "openai_field_identified": False,
            "pde_validated": False,
        }

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return _evaluate_with_backend(self.leading_backend, x, y, z, t).velocity

    def semantic_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent2": {
                "pr": PARENT_AGENT2_PR,
                "head": PARENT_AGENT2_HEAD,
                "batch_axis_safety_sha256": batch_axis_safety_sha256(),
                "batch_differentials_sha256": batch_differentials_sha256(),
            },
            "agent1_leading": {
                "pr": AGENT1_PR,
                "head": AGENT1_HEAD,
                "module": AGENT1_MODULE,
                "class": AGENT1_CLASS,
                "source_blob_sha1": AGENT1_SOURCE_BLOB_SHA1,
                "runtime_semantic_sha256": self._leading_identity["semantic_sha256"],
                "runtime_configuration_sha256": self._leading_identity["configuration_sha256"],
                "X_h": self.X_h,
            },
            "composition": "u_partial = u_lead_current_through_Xh + u_osc_frozen_complete_curl",
            "source_provenance": {
                "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
                "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
                "scope": "localized waves / complete curls plus public leading-coordinate structure",
            },
            "truth_boundary": self.truth_boundary,
        }


@lru_cache(maxsize=1)
def default_field() -> _CurrentPartialLeadingOscillatoryField:
    return _CurrentPartialLeadingOscillatoryField(_load_exact_agent1_backend())


def velocity(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    """Return the current partial-domain ``u_lead + u_osc`` Cartesian velocity."""
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
        "forbidden_inputs_present": sorted(forbidden.intersection(params)),
        "vectorized_velocity_interface": True,
        "agent1_domain_fail_closed": True,
        "a2_support_axis_mask_preserved": True,
        "new_oscillatory_parameters": False,
        "leading_profile_reimplemented": False,
        "mean_projection_performed": False,
        "correction_velocity_constructed": False,
        "pressure_or_forcing_added": False,
        "complete_ns_residual": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def materialize_current_partial_composition_receipt() -> dict[str, Any]:
    field = default_field()
    main = _evaluate_with_backend(
        field.leading_backend,
        _RECEIPT_POINTS[:, 0],
        _RECEIPT_POINTS[:, 1],
        _RECEIPT_POINTS[:, 2],
        _RECEIPT_TIMES,
    )
    axis = _evaluate_with_backend(
        field.leading_backend,
        _AXIS_POINTS[:, 0],
        _AXIS_POINTS[:, 1],
        _AXIS_POINTS[:, 2],
        _AXIS_TIMES,
    )

    closure = float(np.max(np.abs(main.velocity - (main.leading + main.oscillatory))))
    osc_rms = float(np.sqrt(np.mean(np.sum(main.oscillatory * main.oscillatory, axis=-1))))
    delta_rms = float(np.sqrt(np.mean(np.sum((main.velocity - main.leading) ** 2, axis=-1))))
    axis_osc_max = float(np.max(np.abs(axis.oscillatory)))
    axis_total_lead_max = float(np.max(np.abs(axis.velocity - axis.leading)))

    beyond_fail_closed = False
    try:
        field.velocity(1.0e6, 0.0, 0.0, 0.5)
    except (ValueError, RuntimeError):
        beyond_fail_closed = True

    return {
        "schema": RECEIPT_SCHEMA,
        "semantic_sha256": field.semantic_sha256,
        "semantic_payload": field.semantic_payload(),
        "public_contract": public_contract(),
        "diagnostic": {
            "interior_probe_count": int(_RECEIPT_POINTS.shape[0]),
            "axis_probe_count": int(_AXIS_POINTS.shape[0]),
            "composition_closure_max_abs": closure,
            "oscillatory_vector_rms": osc_rms,
            "total_minus_leading_vector_rms": delta_rms,
            "axis_oscillatory_abs_max": axis_osc_max,
            "axis_total_minus_leading_abs_max": axis_total_lead_max,
            "beyond_current_Xh_fail_closed": beyond_fail_closed,
            "current_X_h": field.X_h,
        },
        "frozen_mechanical_gates": {
            "composition_closure_max_abs": COMPOSITION_ATOL,
            "oscillatory_vector_rms_min": OSCILLATORY_SIGNAL_FLOOR,
            "axis_oscillatory_abs_max": 0.0,
            "axis_total_minus_leading_abs_max": COMPOSITION_ATOL,
            "beyond_current_Xh_fail_closed": True,
        },
    }


def enforce_receipt(receipt: Mapping[str, Any]) -> None:
    if receipt.get("schema") != RECEIPT_SCHEMA:
        raise RuntimeError("composition receipt schema drifted")
    contract = receipt.get("public_contract")
    if not isinstance(contract, Mapping) or contract.get("forbidden_inputs_present") != []:
        raise RuntimeError("composition public API exposes a forbidden tuning input")
    payload = receipt.get("semantic_payload")
    if not isinstance(payload, Mapping) or _sha256(payload) != receipt.get("semantic_sha256"):
        raise RuntimeError("composition semantic identity mismatch")
    truth = payload.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise RuntimeError("composition truth boundary is missing")
    required_false = (
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
        "paper_exact",
        "openai_field_identified",
        "pde_validated",
    )
    if any(truth.get(key) is not False for key in required_false):
        raise RuntimeError("partial composition was promoted beyond its scientific scope")

    d = receipt.get("diagnostic")
    if not isinstance(d, Mapping):
        raise RuntimeError("composition diagnostic is missing")
    if float(d["composition_closure_max_abs"]) > COMPOSITION_ATOL:
        raise RuntimeError("leading+oscillatory addition closure failed")
    if float(d["oscillatory_vector_rms"]) < OSCILLATORY_SIGNAL_FLOOR:
        raise RuntimeError("oscillatory contribution is vacuous on frozen probes")
    if abs(float(d["total_minus_leading_vector_rms"]) - float(d["oscillatory_vector_rms"])) > 1.0e-13:
        raise RuntimeError("total-minus-leading does not replay the oscillatory contribution")
    if float(d["axis_oscillatory_abs_max"]) != 0.0:
        raise RuntimeError("A2 axis/support masking was not preserved")
    if float(d["axis_total_minus_leading_abs_max"]) > COMPOSITION_ATOL:
        raise RuntimeError("axis total does not reduce to the Agent-1 leading field")
    if d["beyond_current_Xh_fail_closed"] is not True:
        raise RuntimeError("Agent-1 current X_h domain boundary was bypassed")


def _main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    receipt = materialize_current_partial_composition_receipt()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.check:
        enforce_receipt(receipt)


if __name__ == "__main__":
    _main()
