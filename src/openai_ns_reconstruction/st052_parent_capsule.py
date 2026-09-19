"""Fail-closed materialization manifest for the frozen ST052-M parent.

This module does not claim that a regenerated candidate JSON has the same raw
bytes as the historical ST052-M artifact.  Instead it binds a runnable replay
artifact to the exact source commit and source-native frozen recipe, while
recording the regenerated file checksum separately.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

CANDIDATE_ID = "ST052-M"
SOURCE_REPOSITORY = "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1"
SOURCE_PR = 508
SOURCE_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
SOURCE_RECIPE_PATH = "experiments/root_st052/recipe.json"
SOURCE_RECIPE_GIT_BLOB_SHA1 = "e30c769052379f72afeee46ca264482884cc5ac7"
# Historical source artifact checksum is provenance only: exact replay metadata
# can change candidate.json bytes, as replay_st052.py explicitly documents.
SOURCE_RAW_CHILD_SHA256 = "e078e753fab38ebfa0284d28ba64d26cb8538b700a3afc7705668849c43e12da"
SCHEMA = "st052-parent-replay-capsule/v1"


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _git_blob_sha1(path: Path) -> str:
    payload = path.read_bytes()
    header = f"blob {len(payload)}\0".encode("ascii")
    return hashlib.sha1(header + payload).hexdigest()


def _require_false(mapping: dict[str, Any], key: str, *, where: str) -> None:
    if mapping.get(key) is not False:
        raise ValueError(f"{where}.{key} must be exactly false")


def build_manifest(source_root: Path, candidate: Path, validation: Path) -> dict[str, Any]:
    recipe_path = source_root / SOURCE_RECIPE_PATH
    if not recipe_path.is_file():
        raise FileNotFoundError(recipe_path)
    if not candidate.is_file():
        raise FileNotFoundError(candidate)
    if not validation.is_file():
        raise FileNotFoundError(validation)

    blob_sha1 = _git_blob_sha1(recipe_path)
    if blob_sha1 != SOURCE_RECIPE_GIT_BLOB_SHA1:
        raise ValueError(f"source recipe blob mismatch: {blob_sha1}")

    recipe = json.loads(recipe_path.read_text())
    if recipe.get("parent_id") != "ST051-B":
        raise ValueError("unexpected ST052-M parent_id")
    parent_sha = recipe.get("parent_sha256")
    modifier_sha = recipe.get("modifiers_sha256")
    reference = recipe.get("reference")
    if not isinstance(parent_sha, str) or len(parent_sha) != 64:
        raise ValueError("missing/invalid parent_sha256")
    if not isinstance(modifier_sha, str) or len(modifier_sha) != 64:
        raise ValueError("missing/invalid modifiers_sha256")
    if not isinstance(reference, dict):
        raise ValueError("missing source-native reference table")
    for key in ("points", "times", "velocity", "pressure"):
        if key not in reference:
            raise ValueError(f"reference.{key} missing")
    _require_false(recipe, "pde_validated", where="recipe")
    _require_false(recipe, "source_correspondence_verified", where="recipe")

    report = json.loads(validation.read_text())
    _require_false(report, "pde_validated", where="validation")
    _require_false(report, "all_numeric_gates_pass", where="validation")
    gates = report.get("gates")
    if not isinstance(gates, dict):
        raise ValueError("validation.gates missing")
    if gates.get("momentum_max") is not False or gates.get("momentum_L2") is not False:
        raise ValueError("expected ST052-M momentum rejection is absent")

    reference_sha = _sha256_bytes(_canonical_bytes(reference))
    replay_identity_payload = {
        "schema": SCHEMA,
        "candidate_id": CANDIDATE_ID,
        "source_repository": SOURCE_REPOSITORY,
        "source_pr": SOURCE_PR,
        "source_head": SOURCE_HEAD,
        "source_recipe_git_blob_sha1": SOURCE_RECIPE_GIT_BLOB_SHA1,
        "parent_id": recipe["parent_id"],
        "parent_sha256": parent_sha,
        "modifiers_sha256": modifier_sha,
        "reference_sha256": reference_sha,
    }
    replay_identity = _sha256_bytes(_canonical_bytes(replay_identity_payload))

    return {
        "schema": SCHEMA,
        "candidate_id": CANDIDATE_ID,
        "task_id": "CR-A9-068",
        "migration_class": "direct_internal_exact_source_replay",
        "source": {
            "repository": SOURCE_REPOSITORY,
            "pr": SOURCE_PR,
            "head": SOURCE_HEAD,
            "recipe_path": SOURCE_RECIPE_PATH,
            "recipe_git_blob_sha1": blob_sha1,
            "recipe_sha256": _sha256_file(recipe_path),
            "license_scope": "repository-internal source; no third-party implementation copied",
        },
        "replay_identity": {
            "sha256": replay_identity,
            "scope": "Exact source commit + immutable source recipe identity and frozen source-native reference table; not a global mathematical-function hash.",
            "reference_sha256": reference_sha,
            "parent_id": recipe["parent_id"],
            "parent_sha256": parent_sha,
            "modifiers_sha256": modifier_sha,
        },
        "materialized_candidate": {
            "filename": candidate.name,
            "sha256": _sha256_file(candidate),
            "historical_source_raw_child_sha256": SOURCE_RAW_CHILD_SHA256,
            "historical_raw_sha_identity_claimed": False,
            "runnable_source_replay_artifact": True,
            "save_load_exercised_by_validation": True,
        },
        "validation": {
            "filename": validation.name,
            "sha256": _sha256_file(validation),
            "pde_validated": False,
            "all_numeric_gates_pass": False,
            "momentum_max_gate": False,
            "momentum_L2_gate": False,
        },
        "truth_boundary": {
            "production_candidate_selected": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "source_correspondence_verified": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
            "parent_pressure_or_forcing_receipt_transferred": False,
            "whole_child_save_load_ready": False,
        },
    }


def write_manifest(source_root: Path, candidate: Path, validation: Path, out: Path) -> dict[str, Any]:
    manifest = build_manifest(source_root, candidate, validation)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description="Materialize a fail-closed ST052-M parent replay capsule manifest")
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    manifest = write_manifest(args.source_root, args.candidate, args.validation, args.out)
    print(json.dumps({
        "candidate_id": manifest["candidate_id"],
        "replay_identity_sha256": manifest["replay_identity"]["sha256"],
        "candidate_file_sha256": manifest["materialized_candidate"]["sha256"],
        "pde_validated": manifest["truth_boundary"]["pde_validated"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
