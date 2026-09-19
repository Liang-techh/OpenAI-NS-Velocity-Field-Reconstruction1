"""Versioned whole-candidate capsule for the frozen ST052-M linear temporal child.

This module composes already-governed assets without adding representation
freedom: the exact-source ST052-M parent replay capsule and the exact PR #587
linear temporal transform spec.  The resulting directory is a checksum-bound
save/load bundle for ``velocity(x,y,z,t)`` when the authenticated exact ST052
source runtime is provided.

Important truth boundary: the bundled temporal child changes ``u_t``.  Parent
pressure/forcing/residual receipts are therefore not inherited and this module
does not claim PDE validation or OpenAI hidden-field identity.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import importlib
import json
from pathlib import Path
import shutil
import sys
from typing import Any

import numpy as np

from .st052_linear_temporal_transform import (
    DEFAULT_TEMPORAL_SPEC,
    St052LinearTemporalAdapter,
    St052LinearTemporalTransformSpec,
)
from .st052_parent_capsule import CANDIDATE_ID, SOURCE_HEAD
from .st052_parent_reference_binding import St052ReferenceBoundParent
from .st052_source_runtime_identity import (
    authenticate_source_runtime,
    source_runtime_identity_sha256,
)

SCHEMA = "st052-linear-temporal-whole-candidate-capsule/v2"
TASK_ID = "CR-A9-070"
PARENT_CANDIDATE_FILENAME = "parent_candidate.json"
PARENT_VALIDATION_FILENAME = "parent_validation.json"
PARENT_MANIFEST_FILENAME = "parent_replay_manifest.json"
TEMPORAL_SPEC_FILENAME = "temporal_transform_spec.json"
WHOLE_MANIFEST_FILENAME = "whole_candidate_manifest.json"


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    ).encode("utf-8")


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError(f"expected JSON object: {path}")
    return obj


def _require_false(mapping: dict[str, Any], key: str, *, where: str) -> None:
    if mapping.get(key) is not False:
        raise ValueError(f"{where}.{key} must be exactly false")


def _identity_payload(
    *,
    parent_replay_identity: str,
    parent_candidate_sha256: str,
    parent_validation_sha256: str,
    parent_manifest_sha256: str,
    temporal_spec_sha256: str,
    static_transform_spec_sha256: str,
    exact_source_runtime_identity_sha256: str,
) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "parent_candidate_id": CANDIDATE_ID,
        "parent_source_head": SOURCE_HEAD,
        "parent_replay_identity_sha256": parent_replay_identity,
        "parent_candidate_file_sha256": parent_candidate_sha256,
        "parent_validation_file_sha256": parent_validation_sha256,
        "parent_manifest_file_sha256": parent_manifest_sha256,
        "temporal_transform_spec_sha256": temporal_spec_sha256,
        "static_transform_spec_sha256": static_transform_spec_sha256,
        "exact_source_runtime_identity_sha256": exact_source_runtime_identity_sha256,
    }


def build_bundle(
    *,
    parent_candidate: str | Path,
    parent_validation: str | Path,
    parent_manifest: str | Path,
    out_dir: str | Path,
    temporal_spec: St052LinearTemporalTransformSpec = DEFAULT_TEMPORAL_SPEC,
) -> dict[str, Any]:
    """Create a deterministic checksum-bound whole-child bundle directory."""
    parent_candidate = Path(parent_candidate)
    parent_validation = Path(parent_validation)
    parent_manifest = Path(parent_manifest)
    out_dir = Path(out_dir)
    for path in (parent_candidate, parent_validation, parent_manifest):
        if not path.is_file():
            raise FileNotFoundError(path)
    if out_dir.exists() and any(out_dir.iterdir()):
        raise ValueError("whole-candidate output directory must be empty")
    out_dir.mkdir(parents=True, exist_ok=True)

    temporal_spec.validate_static_binding()
    pm = _read_json(parent_manifest)
    if pm.get("candidate_id") != CANDIDATE_ID:
        raise ValueError("parent replay manifest candidate identity mismatch")
    source = pm.get("source")
    replay = pm.get("replay_identity")
    materialized = pm.get("materialized_candidate")
    validation_record = pm.get("validation")
    truth = pm.get("truth_boundary")
    if not all(isinstance(v, dict) for v in (source, replay, materialized, validation_record, truth)):
        raise ValueError("malformed parent replay manifest")
    if source.get("head") != SOURCE_HEAD:
        raise ValueError("parent replay manifest source-head mismatch")
    _require_false(truth, "pde_validated", where="parent_manifest.truth_boundary")
    _require_false(truth, "source_correspondence_verified", where="parent_manifest.truth_boundary")

    candidate_sha = _sha256_file(parent_candidate)
    validation_sha = _sha256_file(parent_validation)
    if materialized.get("sha256") != candidate_sha:
        raise ValueError("parent candidate bytes do not match replay manifest")
    if validation_record.get("sha256") != validation_sha:
        raise ValueError("parent validation bytes do not match replay manifest")

    validation_obj = _read_json(parent_validation)
    _require_false(validation_obj, "pde_validated", where="parent_validation")
    _require_false(validation_obj, "all_numeric_gates_pass", where="parent_validation")
    gates = validation_obj.get("gates")
    if not isinstance(gates, dict) or gates.get("momentum_max") is not False or gates.get("momentum_L2") is not False:
        raise ValueError("expected frozen ST052-M momentum rejection is absent")

    dst_candidate = out_dir / PARENT_CANDIDATE_FILENAME
    dst_validation = out_dir / PARENT_VALIDATION_FILENAME
    dst_parent_manifest = out_dir / PARENT_MANIFEST_FILENAME
    dst_temporal_spec = out_dir / TEMPORAL_SPEC_FILENAME
    shutil.copyfile(parent_candidate, dst_candidate)
    shutil.copyfile(parent_validation, dst_validation)
    shutil.copyfile(parent_manifest, dst_parent_manifest)
    temporal_spec.save(dst_temporal_spec)

    parent_manifest_sha = _sha256_file(dst_parent_manifest)
    temporal_spec_sha = temporal_spec.sha256()
    runtime_id = source_runtime_identity_sha256()
    identity_payload = _identity_payload(
        parent_replay_identity=str(replay.get("sha256")),
        parent_candidate_sha256=_sha256_file(dst_candidate),
        parent_validation_sha256=_sha256_file(dst_validation),
        parent_manifest_sha256=parent_manifest_sha,
        temporal_spec_sha256=temporal_spec_sha,
        static_transform_spec_sha256=temporal_spec.static_transform_spec_sha256,
        exact_source_runtime_identity_sha256=runtime_id,
    )
    whole_id = _sha256_bytes(_canonical_bytes(identity_payload))
    manifest = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "candidate_id": "ST052-M-linear-temporal-child-v1",
        "migration_class": "direct_internal_runtime_identity_binding",
        "whole_candidate_identity_sha256": whole_id,
        "identity_payload": identity_payload,
        "files": {
            PARENT_CANDIDATE_FILENAME: _sha256_file(dst_candidate),
            PARENT_VALIDATION_FILENAME: _sha256_file(dst_validation),
            PARENT_MANIFEST_FILENAME: parent_manifest_sha,
            TEMPORAL_SPEC_FILENAME: _sha256_file(dst_temporal_spec),
        },
        "runtime": {
            "callable": "velocity(x,y,z,t)->[...,3]",
            "exact_parent_source_runtime_required": True,
            "exact_source_runtime_identity_sha256": runtime_id,
            "exact_source_runtime_identity_closed": True,
            "standalone_package_parent_runtime_ready": False,
        },
        "truth_boundary": {
            "whole_child_bundle_materialized": True,
            "whole_child_save_load_ready_with_exact_source_runtime": True,
            "exact_source_runtime_identity_closed": True,
            "standalone_package_parent_runtime_ready": False,
            "velocity_export_ready": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "production_candidate_selected": False,
            "pressure_or_force_changed": False,
            "parent_pressure_or_forcing_receipt_transferred": False,
            "held_out_temporal_child_pde_residual_evaluated": False,
            "pde_validated": False,
            "source_correspondence_verified": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    (out_dir / WHOLE_MANIFEST_FILENAME).write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return manifest


def verify_bundle(bundle_dir: str | Path) -> dict[str, Any]:
    """Recompute every file digest and the composite whole-candidate identity."""
    bundle_dir = Path(bundle_dir)
    manifest_path = bundle_dir / WHOLE_MANIFEST_FILENAME
    manifest = _read_json(manifest_path)
    if manifest.get("schema") != SCHEMA or manifest.get("task_id") != TASK_ID:
        raise ValueError("unsupported ST052 whole-candidate capsule schema/task")
    files = manifest.get("files")
    identity = manifest.get("identity_payload")
    if not isinstance(files, dict) or not isinstance(identity, dict):
        raise ValueError("malformed whole-candidate manifest")
    for name in (
        PARENT_CANDIDATE_FILENAME,
        PARENT_VALIDATION_FILENAME,
        PARENT_MANIFEST_FILENAME,
        TEMPORAL_SPEC_FILENAME,
    ):
        path = bundle_dir / name
        if not path.is_file() or files.get(name) != _sha256_file(path):
            raise ValueError(f"whole-candidate bundle file checksum mismatch: {name}")

    temporal_spec = St052LinearTemporalTransformSpec.load(bundle_dir / TEMPORAL_SPEC_FILENAME)
    parent_manifest = _read_json(bundle_dir / PARENT_MANIFEST_FILENAME)
    replay = parent_manifest.get("replay_identity")
    if not isinstance(replay, dict):
        raise ValueError("parent replay identity missing")
    runtime_id = source_runtime_identity_sha256()
    expected_identity = _identity_payload(
        parent_replay_identity=str(replay.get("sha256")),
        parent_candidate_sha256=_sha256_file(bundle_dir / PARENT_CANDIDATE_FILENAME),
        parent_validation_sha256=_sha256_file(bundle_dir / PARENT_VALIDATION_FILENAME),
        parent_manifest_sha256=_sha256_file(bundle_dir / PARENT_MANIFEST_FILENAME),
        temporal_spec_sha256=temporal_spec.sha256(),
        static_transform_spec_sha256=temporal_spec.static_transform_spec_sha256,
        exact_source_runtime_identity_sha256=runtime_id,
    )
    if identity != expected_identity:
        raise ValueError("whole-candidate identity payload mismatch")
    expected_id = _sha256_bytes(_canonical_bytes(expected_identity))
    if manifest.get("whole_candidate_identity_sha256") != expected_id:
        raise ValueError("whole-candidate identity checksum mismatch")
    runtime = manifest.get("runtime")
    if not isinstance(runtime, dict):
        raise ValueError("whole-candidate runtime metadata missing")
    if runtime.get("exact_source_runtime_identity_sha256") != runtime_id:
        raise ValueError("whole-candidate exact source runtime identity mismatch")
    if runtime.get("exact_source_runtime_identity_closed") is not True:
        raise ValueError("whole-candidate runtime identity closure flag missing")

    truth = manifest.get("truth_boundary")
    if not isinstance(truth, dict):
        raise ValueError("whole-candidate truth boundary missing")
    for key in (
        "pde_validated",
        "source_correspondence_verified",
        "visualization_ready",
        "visual_correspondence_verified",
        "production_candidate_selected",
        "parent_pressure_or_forcing_receipt_transferred",
    ):
        _require_false(truth, key, where="whole_manifest.truth_boundary")
    if truth.get("exact_source_runtime_identity_closed") is not True:
        raise ValueError("whole-candidate runtime identity truth flag missing")
    return manifest


def _load_exact_source_parent(source_root: Path, candidate_path: Path):
    """Load bundled ST052-M bytes only after authenticating the exact #508 runtime."""
    source_root = source_root.resolve()
    authenticate_source_runtime(source_root)

    st052_dir = source_root / "experiments/root_st052"
    inserted = str(st052_dir)
    sys.path.insert(0, inserted)
    try:
        # Source identity is authenticated before any historical code executes.
        # Importing replay_st052 installs the predecessor-chain paths used by
        # the exact source replay; Family.load dispatches the frozen basis kind
        # stored in the candidate bytes.
        sys.modules.pop("replay_st052", None)
        replay = importlib.import_module("replay_st052")
        family, raw = replay.Family.load(candidate_path)
    finally:
        if sys.path and sys.path[0] == inserted:
            sys.path.pop(0)

    def base_velocity(points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        times = np.full(len(points), float(time), dtype=float)
        velocity, _ = family.fields(raw, points, times)
        return np.asarray(velocity, dtype=float)

    return St052ReferenceBoundParent(base_velocity)


@dataclass
class St052LinearTemporalWholeCandidate:
    """Loaded whole-child runtime backed by bundled parent bytes + exact source."""

    adapter: St052LinearTemporalAdapter
    manifest: dict[str, Any]

    @property
    def candidate_id(self) -> str:
        return str(self.manifest["candidate_id"])

    @property
    def identity_sha256(self) -> str:
        return str(self.manifest["whole_candidate_identity_sha256"])

    def velocity(self, x, y, z, t) -> np.ndarray:
        return self.adapter.velocity(x, y, z, t)


def load_bundle_runtime(
    bundle_dir: str | Path,
    *,
    exact_source_root: str | Path,
) -> St052LinearTemporalWholeCandidate:
    """Verify + reload a bundle into the final broadcastable velocity callable.

    The parent evaluator still comes from the historical source checkout, but
    the checkout is now authenticated against the exact #508 Git commit/tree
    and clean-worktree contract before source code is imported.  This closes
    runtime substitution while keeping the separate truth that the parent
    evaluator has not yet been ported into the installable package.
    """
    bundle_dir = Path(bundle_dir)
    manifest = verify_bundle(bundle_dir)
    temporal_spec = St052LinearTemporalTransformSpec.load(
        bundle_dir / TEMPORAL_SPEC_FILENAME
    )
    parent = _load_exact_source_parent(
        Path(exact_source_root), bundle_dir / PARENT_CANDIDATE_FILENAME
    )
    adapter = St052LinearTemporalAdapter(
        parent,
        parent_candidate_id=temporal_spec.parent_candidate_id,
        parent_source_head=temporal_spec.parent_source_head,
        spec=temporal_spec,
    )
    return St052LinearTemporalWholeCandidate(adapter=adapter, manifest=manifest)


TRUTH_BOUNDARY = {
    "whole_child_bundle_materialized": True,
    "whole_child_save_load_ready_with_exact_source_runtime": True,
    "exact_source_runtime_identity_closed": True,
    "standalone_package_parent_runtime_ready": False,
    "velocity_export_ready": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "production_candidate_selected": False,
    "pressure_or_force_changed": False,
    "parent_pressure_or_forcing_receipt_transferred": False,
    "held_out_temporal_child_pde_residual_evaluated": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}
