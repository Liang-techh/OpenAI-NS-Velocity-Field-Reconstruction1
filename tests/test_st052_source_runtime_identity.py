from __future__ import annotations

import hashlib
from pathlib import Path
import subprocess
import sys
import types

import pytest

import openai_ns_reconstruction.st052_source_runtime_identity as identity


def _run(root: Path, *args: str) -> str:
    p = subprocess.run(
        ["git", "-C", str(root), *args],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return p.stdout.strip()


def _blob(path: Path) -> str:
    payload = path.read_bytes()
    return hashlib.sha1(f"blob {len(payload)}\0".encode() + payload).hexdigest()


def _fake_exact_checkout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    root = tmp_path / "source"
    st052 = root / "experiments" / "root_st052"
    st052.mkdir(parents=True)
    recipe = st052 / "recipe.json"
    replay = st052 / "replay_st052.py"
    dependency = st052 / "minimax_exchange.py"
    recipe.write_text('{"candidate":"ST052-M"}\n')
    dependency.write_text("VALUE = 'frozen-dependency'\n")
    replay.write_text(
        "from minimax_exchange import VALUE as DEP_VALUE\n"
        "VALUE = 'frozen:' + DEP_VALUE\n"
    )
    (root / ".gitignore").write_text("ignored_shadow.py\n")
    _run(root, "init")
    _run(root, "add", ".")
    subprocess.run(
        [
            "git", "-C", str(root),
            "-c", "user.name=Agent9",
            "-c", "user.email=agent9@example.invalid",
            "commit", "-m", "frozen source",
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    monkeypatch.setattr(identity, "SOURCE_HEAD", _run(root, "rev-parse", "HEAD"))
    monkeypatch.setattr(identity, "SOURCE_TREE", _run(root, "rev-parse", "HEAD^{tree}"))
    monkeypatch.setattr(identity, "SOURCE_RECIPE_GIT_BLOB_SHA1", _blob(recipe))
    return root


def test_runtime_identity_is_deterministic():
    a = identity.source_runtime_identity_sha256()
    b = identity.source_runtime_identity_sha256()
    assert len(a) == 64
    assert a == b
    payload = identity.source_runtime_identity_payload()
    assert payload["source_head"] == identity.SOURCE_HEAD
    assert payload["source_tree"] == identity.SOURCE_TREE
    assert payload["replay_entrypoint"]["covered_by_source_tree"] is True
    assert payload["module_import_guard"]["loaded_module_origin_and_blob_verified"] is True


def test_authenticator_accepts_frozen_clean_checkout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root = _fake_exact_checkout(tmp_path, monkeypatch)
    receipt = identity.authenticate_source_runtime(root)
    assert receipt["authenticated_head"] == identity.SOURCE_HEAD
    assert receipt["authenticated_tree"] == identity.SOURCE_TREE
    assert receipt["tracked_worktree_clean"] is True
    assert receipt["untracked_python_sources_absent"] is True


def test_authenticator_rejects_tracked_runtime_tamper(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root = _fake_exact_checkout(tmp_path, monkeypatch)
    replay = root / identity.SOURCE_REPLAY_PATH
    replay.write_text(replay.read_text() + "# tamper\n")
    with pytest.raises(ValueError, match="tracked modifications"):
        identity.authenticate_source_runtime(root)


def test_authenticator_rejects_untracked_python_shadow(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root = _fake_exact_checkout(tmp_path, monkeypatch)
    (root / "shadow.py").write_text("VALUE = 'shadow'\n")
    with pytest.raises(ValueError, match="untracked Python"):
        identity.authenticate_source_runtime(root)


def test_authenticator_rejects_ignored_python_shadow(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root = _fake_exact_checkout(tmp_path, monkeypatch)
    (root / "ignored_shadow.py").write_text("VALUE = 'ignored shadow'\n")
    with pytest.raises(ValueError, match="ignored untracked Python"):
        identity.authenticate_source_runtime(root)


def test_authenticated_import_isolates_preloaded_transitive_module(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    root = _fake_exact_checkout(tmp_path, monkeypatch)
    poison = types.ModuleType("minimax_exchange")
    poison.__file__ = str(tmp_path / "foreign" / "minimax_exchange.py")
    poison.VALUE = "poisoned-cache"
    monkeypatch.setitem(sys.modules, "minimax_exchange", poison)

    try:
        replay, receipt = identity.import_authenticated_replay(root)
        assert replay.VALUE == "frozen:frozen-dependency"
        assert sys.modules["minimax_exchange"] is not poison
        assert receipt["module_cache_isolated"] is True
        assert "minimax_exchange" in receipt["evicted_preexisting_module_names"]
        records = {row["module_name"]: row for row in receipt["loaded_source_modules"]}
        assert records["replay_st052"]["relative_path"] == identity.SOURCE_REPLAY_PATH
        assert records["minimax_exchange"]["relative_path"].endswith(
            "experiments/root_st052/minimax_exchange.py"
        )
    finally:
        sys.modules.pop("replay_st052", None)
        sys.modules.pop("minimax_exchange", None)
