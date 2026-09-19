"""Fail-closed identity guard for the frozen ST052-M source runtime.

The whole-candidate capsule still relies on the historical PR #508 evaluator.
A recipe checksum alone is not enough to identify that evaluator: a caller could
point the loader at another checkout or at a modified worktree that retained the
same recipe bytes.  This module authenticates the immutable Git commit/tree and
requires the tracked worktree to be clean before any source module is imported.

This is delivery/reproducibility governance only.  It does not change the
velocity field, pressure, forcing, residual operator, or any scientific gate.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any

SOURCE_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
SOURCE_TREE = "cbaf188a7079e17984b0f78d7c8c38df8e11401b"
SOURCE_RECIPE_PATH = "experiments/root_st052/recipe.json"
SOURCE_RECIPE_GIT_BLOB_SHA1 = "e30c769052379f72afeee46ca264482884cc5ac7"
SOURCE_REPLAY_PATH = "experiments/root_st052/replay_st052.py"
SCHEMA = "st052-exact-source-runtime-identity/v1"
TASK_ID = "CR-A9-070"


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    ).encode("utf-8")


def _git_blob_sha1(path: Path) -> str:
    payload = path.read_bytes()
    header = f"blob {len(payload)}\0".encode("ascii")
    return hashlib.sha1(header + payload).hexdigest()


def source_runtime_identity_payload() -> dict[str, Any]:
    """Return the immutable runtime identity expected by the package."""
    return {
        "schema": SCHEMA,
        "source_head": SOURCE_HEAD,
        "source_tree": SOURCE_TREE,
        "recipe": {
            "path": SOURCE_RECIPE_PATH,
            "git_blob_sha1": SOURCE_RECIPE_GIT_BLOB_SHA1,
        },
        "replay_entrypoint": {
            "path": SOURCE_REPLAY_PATH,
            "covered_by_source_tree": True,
        },
    }


def source_runtime_identity_sha256() -> str:
    return hashlib.sha256(_canonical_bytes(source_runtime_identity_payload())).hexdigest()


def _git(root: Path, *args: str) -> str:
    try:
        proc = subprocess.run(
            ["git", "-C", str(root), *args],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ValueError(f"unable to authenticate exact ST052 source runtime: git {' '.join(args)}") from exc
    return proc.stdout.strip()


def authenticate_source_runtime(source_root: str | Path) -> dict[str, Any]:
    """Authenticate the exact #508 checkout before executing any source code.

    The contract requires the supplied path itself to be the Git worktree root,
    the exact frozen commit/tree to be checked out, no tracked modifications or
    staged changes, and no untracked Python sources (including ignored ones)
    that could shadow source imports.  The historical recipe keeps its existing
    source-native blob receipt; all tracked runtime code is cryptographically
    covered by the exact commit/tree plus the clean-worktree requirement.
    """
    root = Path(source_root).resolve()
    if not root.is_dir():
        raise FileNotFoundError(root)

    top = Path(_git(root, "rev-parse", "--show-toplevel")).resolve()
    if top != root:
        raise ValueError("exact ST052 source root must be the Git worktree root")
    head = _git(root, "rev-parse", "HEAD")
    tree = _git(root, "rev-parse", "HEAD^{tree}")
    if head != SOURCE_HEAD:
        raise ValueError(f"exact ST052 source HEAD mismatch: {head}")
    if tree != SOURCE_TREE:
        raise ValueError(f"exact ST052 source tree mismatch: {tree}")

    tracked = _git(root, "status", "--porcelain=v1", "--untracked-files=no")
    if tracked:
        raise ValueError("exact ST052 source worktree has tracked modifications")
    untracked_py = _git(root, "ls-files", "--others", "--exclude-standard", "--", "*.py")
    if untracked_py:
        raise ValueError("exact ST052 source worktree has untracked Python sources")
    ignored_untracked_py = _git(
        root, "ls-files", "--others", "--ignored", "--exclude-standard", "--", "*.py"
    )
    if ignored_untracked_py:
        raise ValueError("exact ST052 source worktree has ignored untracked Python sources")

    recipe = root / SOURCE_RECIPE_PATH
    replay = root / SOURCE_REPLAY_PATH
    if not recipe.is_file() or _git_blob_sha1(recipe) != SOURCE_RECIPE_GIT_BLOB_SHA1:
        raise ValueError("exact ST052 source recipe identity mismatch")
    if not replay.is_file():
        raise ValueError("exact ST052 replay entrypoint is missing")

    payload = source_runtime_identity_payload()
    return {
        "task_id": TASK_ID,
        "identity_payload": payload,
        "identity_sha256": source_runtime_identity_sha256(),
        "authenticated_head": head,
        "authenticated_tree": tree,
        "tracked_worktree_clean": True,
        "untracked_python_sources_absent": True,
    }


TRUTH_BOUNDARY = {
    "exact_source_runtime_identity_closed": True,
    "standalone_package_parent_runtime_ready": False,
    "velocity_export_ready": False,
    "visualization_ready": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "openai_field_identified": False,
}
