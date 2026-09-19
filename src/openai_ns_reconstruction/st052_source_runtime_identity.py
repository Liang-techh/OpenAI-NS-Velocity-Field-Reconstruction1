"""Fail-closed identity guard for the frozen ST052-M source runtime.

The whole-candidate capsule still relies on the historical PR #508 evaluator.
A recipe checksum alone is not enough to identify that evaluator: a caller could
point the loader at another checkout or at a modified worktree that retained the
same recipe bytes.  This module authenticates the immutable Git commit/tree,
requires the tracked worktree to be clean, and isolates the historical top-level
module graph from any pre-existing ``sys.modules`` entries before source code is
used.

This is delivery/reproducibility governance only.  It does not change the
velocity field, pressure, forcing, residual operator, or any scientific gate.
"""
from __future__ import annotations

import hashlib
import importlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Any, Iterable

SOURCE_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
SOURCE_TREE = "cbaf188a7079e17984b0f78d7c8c38df8e11401b"
SOURCE_RECIPE_PATH = "experiments/root_st052/recipe.json"
SOURCE_RECIPE_GIT_BLOB_SHA1 = "e30c769052379f72afeee46ca264482884cc5ac7"
SOURCE_REPLAY_PATH = "experiments/root_st052/replay_st052.py"
SCHEMA = "st052-exact-source-runtime-identity/v2"
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
        "module_import_guard": {
            "strategy": "evict_then_verify_tracked_experiments_modules",
            "candidate_scope": "tracked experiments/**/*.py basenames",
            "loaded_module_origin_and_blob_verified": True,
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
        raise ValueError(
            f"unable to authenticate exact ST052 source runtime: git {' '.join(args)}"
        ) from exc
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


def source_python_module_names(source_root: str | Path) -> frozenset[str]:
    """Return importable top-level basenames from tracked historical Python files.

    The PR #508 replay chain uses top-level imports while prepending predecessor
    experiment directories to ``sys.path``.  Enumerating tracked experiment
    basenames gives a conservative cache-eviction set without hard-coding the
    evolving predecessor chain.
    """
    root = Path(source_root).resolve()
    tracked = _git(root, "ls-files", "experiments")
    names: set[str] = set()
    for rel in tracked.splitlines():
        path = Path(rel)
        if path.suffix != ".py" or path.name == "__init__.py":
            continue
        stem = path.stem
        if stem.isidentifier():
            names.add(stem)
    return frozenset(names)


def verify_source_module_cache(
    source_root: str | Path,
    module_names: Iterable[str] | None = None,
) -> tuple[dict[str, str], ...]:
    """Verify loaded historical top-level modules originate from the exact tree.

    Any loaded name in the conservative historical-module set must resolve to a
    tracked ``.py`` file under the authenticated worktree and its current bytes
    must equal the blob recorded by ``HEAD``.  This rejects a foreign or
    in-memory-cache-selected module whose name happens to match the historical
    replay graph.
    """
    root = Path(source_root).resolve()
    names = source_python_module_names(root) if module_names is None else frozenset(module_names)
    records: list[dict[str, str]] = []
    for name in sorted(names):
        module = sys.modules.get(name)
        if module is None:
            continue
        origin = getattr(module, "__file__", None)
        if not origin:
            spec = getattr(module, "__spec__", None)
            origin = getattr(spec, "origin", None) if spec is not None else None
        if not origin or origin in {"built-in", "frozen"}:
            raise ValueError(f"historical source module has unverifiable origin: {name}")
        path = Path(origin).resolve()
        try:
            rel = path.relative_to(root)
        except ValueError as exc:
            raise ValueError(
                f"historical source module escaped exact worktree: {name} -> {path}"
            ) from exc
        if path.suffix != ".py" or not path.is_file():
            raise ValueError(f"historical source module is not tracked Python source: {name}")
        rel_text = rel.as_posix()
        expected_blob = _git(root, "rev-parse", f"HEAD:{rel_text}")
        actual_blob = _git_blob_sha1(path)
        if expected_blob != actual_blob:
            raise ValueError(f"historical source module blob mismatch: {name}")
        records.append(
            {
                "module_name": name,
                "relative_path": rel_text,
                "git_blob_sha1": actual_blob,
            }
        )
    return tuple(records)


def import_authenticated_replay(source_root: str | Path):
    """Import exact PR #508 replay code without trusting inherited module cache.

    Before import, every tracked experiment-module basename is evicted from
    ``sys.modules``.  The historical replay is then imported with its exact
    source directory first on ``sys.path``.  Every historical module that was
    loaded is required to originate from a tracked file in the authenticated
    tree with matching Git-blob bytes.  On failure, the pre-import ``sys.path``
    is restored and partial historical modules are removed.

    On success the exact historical modules remain cached deliberately: this
    prevents later lazy imports in the returned historical runtime from falling
    back to a foreign module that was present before the authenticated load.
    """
    root = Path(source_root).resolve()
    pre = authenticate_source_runtime(root)
    candidate_names = source_python_module_names(root)
    if "replay_st052" not in candidate_names:
        raise ValueError("exact ST052 replay module is not tracked in source tree")

    evicted = sorted(name for name in candidate_names if name in sys.modules)
    for name in candidate_names:
        sys.modules.pop(name, None)

    st052_dir = root / "experiments" / "root_st052"
    before_path = list(sys.path)
    sys.path.insert(0, str(st052_dir))
    importlib.invalidate_caches()
    try:
        replay = importlib.import_module("replay_st052")
        records = verify_source_module_cache(root, candidate_names)
        replay_records = [r for r in records if r["module_name"] == "replay_st052"]
        if len(replay_records) != 1 or replay_records[0]["relative_path"] != SOURCE_REPLAY_PATH:
            raise ValueError("authenticated replay module did not resolve to frozen entrypoint")
        post = authenticate_source_runtime(root)
    except Exception:
        for name in candidate_names:
            sys.modules.pop(name, None)
        sys.path[:] = before_path
        raise

    return replay, {
        "identity_sha256": pre["identity_sha256"],
        "authenticated_head": post["authenticated_head"],
        "authenticated_tree": post["authenticated_tree"],
        "evicted_preexisting_module_names": evicted,
        "loaded_source_modules": list(records),
        "module_cache_isolated": True,
    }


TRUTH_BOUNDARY = {
    "exact_source_runtime_identity_closed": True,
    "historical_module_cache_isolated": True,
    "standalone_package_parent_runtime_ready": False,
    "velocity_export_ready": False,
    "visualization_ready": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "openai_field_identified": False,
}
