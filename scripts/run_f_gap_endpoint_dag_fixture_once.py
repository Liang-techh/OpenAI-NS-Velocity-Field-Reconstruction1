"""Enter the independent F gap-endpoint DAG fixture exactly once.

The runner authenticates the complete three-file raw manifest before the
exclusive entry receipt.  Only after entry does it read the fixed v2 result
and import the new DAG-fixture leaf.  The fixture is therefore checked against
the previous emitted candidate DAG, while no old producer, evaluator, source
graph, or source value is called or traversed.

Success and failure receipts are persisted before their short summaries are
printed.  A successful fixture result must be JSON-safe and must retain all
four non-admission flags as ``False``.
"""

from __future__ import annotations

import hashlib
import importlib.util
import inspect
import json
import sys
import traceback
from pathlib import Path
from typing import Any, Callable, Mapping


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
COMPONENT = "f_gap_endpoint_dag_fixture_v1"
PHASE = "f_gap_endpoint_dag_fixture_v1_pending"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"
LEAF = ROOT / "src/openai_ns_reconstruction/f_gap_endpoint_dag_fixture.py"

PRIOR_RESULT_REL = "tmp/f_gap_endpoint_candidate_v2/result.json"
PRIOR_RESULT = ROOT / PRIOR_RESULT_REL
PRIOR_RESULT_SHA256 = "70210a7f3582d103e017d60f76a3636d8f42f7f0facea2912e4430ccd661171b"

# This component intentionally has no old source graph or old producer in its
# input scope.  The v2 result is an authenticated, already-emitted candidate
# DAG receipt, not a source program to import or execute.
REQUIRED_INPUTS = {
    "src/openai_ns_reconstruction/f_gap_endpoint_dag_fixture.py",
    "scripts/run_f_gap_endpoint_dag_fixture_once.py",
    PRIOR_RESULT_REL,
}


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _json_bytes(value: Any) -> bytes:
    """Return canonical JSON bytes, rejecting non-JSON values and NaN."""

    return (
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def write_once(path: Path, value: Mapping[str, Any]) -> str:
    raw = _json_bytes(value)
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _repository_path(raw_path: Any) -> tuple[str, Path]:
    require(isinstance(raw_path, str) and raw_path, "manifest path is missing")
    relative = Path(raw_path)
    require(not relative.is_absolute() and ".." not in relative.parts,
            "manifest path must be repository-relative")
    path = (ROOT / relative).resolve()
    require(path.is_relative_to(ROOT.resolve()) and path.is_file(),
            "manifest path is missing or escapes repository: " + raw_path)
    return raw_path.replace("\\", "/"), path


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    """Authenticate every raw manifest file before entry or prior JSON reads."""

    raw_entries = source.get("files")
    entries: list[tuple[Any, Any]] = []
    if isinstance(raw_entries, list):
        for item in raw_entries:
            require(isinstance(item, Mapping), "manifest file entry is not an object")
            entries.append((item.get("path"), item.get("sha256")))
    else:
        mapping = source.get("input_hashes", source.get("sha256"))
        require(isinstance(mapping, Mapping) and mapping,
                "manifest has no files or input_hashes")
        entries.extend(mapping.items())

    hashes: dict[str, str] = {}
    for raw_path, expected in entries:
        key, path = _repository_path(raw_path)
        require(key not in hashes, "duplicate manifest path: " + key)
        require(
            isinstance(expected, str)
            and len(expected) == 64
            and all(c in "0123456789abcdefABCDEF" for c in expected),
            "invalid SHA-256 for " + key,
        )
        actual = sha(path)
        require(actual.lower() == expected.lower(), "raw SHA mismatch: " + key)
        hashes[key] = actual

    required = {item.lower() for item in REQUIRED_INPUTS}
    present = {item.lower() for item in hashes}
    require(
        present == required,
        "manifest must contain exactly the new leaf, new runner, and v2 result",
    )
    require(
        hashes.get(PRIOR_RESULT_REL, "").lower() == PRIOR_RESULT_SHA256,
        "v2 result is not pinned to the fixed emitted-DAG receipt",
    )
    return hashes


def _checkpoint_sha(brief: Mapping[str, Any]) -> str:
    checkpoint = brief.get("source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        checkpoint = brief.get(f"{COMPONENT}_source_checkpoint")
    require(isinstance(checkpoint, Mapping), "pending brief source checkpoint missing")
    raw_path = checkpoint.get("path")
    require(isinstance(raw_path, str), "pending brief source checkpoint path missing")
    checkpoint_path = Path(raw_path)
    if not checkpoint_path.is_absolute():
        checkpoint_path = ROOT / checkpoint_path
    require(
        checkpoint_path.resolve() == SOURCE.resolve(),
        "pending brief source checkpoint path mismatch",
    )
    expected = checkpoint.get("sha256")
    require(
        isinstance(expected, str) and expected.lower() == sha(SOURCE),
        "pending brief source checkpoint SHA mismatch",
    )
    return sha(SOURCE)


def preflight() -> tuple[dict[str, Any], dict[str, str], str]:
    require(not ENTRY.exists() and not RESULT.exists(),
            "one-entry component already entered")
    require(BRIEF.is_file(), "pending brief is missing")
    require(SOURCE.is_file(), "source manifest is missing")
    brief = json.loads(BRIEF.read_bytes())
    source = json.loads(SOURCE.read_bytes())
    require(isinstance(brief, Mapping) and isinstance(source, Mapping),
            "brief or source manifest is not an object")
    require(brief.get("component") == COMPONENT and brief.get("phase") == PHASE,
            "pending brief component or phase mismatch")
    require(
        brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("do_not_reenter") is True,
        "pending brief goal flags mismatch",
    )
    source_sha = _checkpoint_sha(brief)
    require(
        source.get("component") == COMPONENT
        and source.get("one_entry_only") is True
        and source.get("no_source_graph_traversal") is True
        and source.get("all_raw_inputs_authenticated_before_decompression") is True,
        "source manifest safeguards missing",
    )
    hashes = _manifest_hashes(source)
    require(LEAF.is_file(), "new DAG-fixture leaf is missing")
    return dict(source), hashes, source_sha


def _load_new_module() -> Any:
    """Import only the new leaf, and only after the exclusive entry receipt."""

    spec = importlib.util.spec_from_file_location(
        "f_gap_endpoint_dag_fixture_new_leaf", LEAF
    )
    require(spec is not None and spec.loader is not None,
            "new DAG-fixture leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _read_prior_result_after_entry() -> tuple[dict[str, Any], Mapping[str, Any]]:
    """Read the fixed v2 receipt only after entry and pin its emitted DAG."""

    require(PRIOR_RESULT.is_file(), "v2 candidate result is missing")
    actual = sha(PRIOR_RESULT)
    require(
        actual.lower() == PRIOR_RESULT_SHA256,
        "v2 candidate result SHA-256 changed",
    )
    value = json.loads(PRIOR_RESULT.read_bytes())
    require(isinstance(value, Mapping), "v2 candidate result is not an object")
    require(
        value.get("component") == "f_gap_endpoint_candidate_v2"
        and value.get("status") == "returned",
        "v2 candidate result status or component mismatch",
    )
    payload = value.get("result")
    require(isinstance(payload, Mapping), "v2 candidate result payload is missing")
    emitted = payload.get("source_role_candidate")
    require(isinstance(emitted, Mapping),
            "v2 result has no emitted source-role candidate DAG")
    return dict(value), emitted


def _find_fixture(module: Any) -> Callable[..., Any]:
    """Select one explicit emitted-DAG fixture entry point, at most once."""

    for name in (
        "verify_emitted_dag_synthetic_fixture",
        "verify_emitted_dag_fixture",
        "verify_synthetic_fixture",
        "synthetic_fixture",
        "run_fixture",
    ):
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate
    raise RuntimeError("new DAG-fixture leaf has no synthetic fixture function")


def _call_fixture(
    function: Callable[..., Any],
    *,
    root: Path,
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    prior_result: Mapping[str, Any],
    emitted_dag: Mapping[str, Any],
) -> Any:
    """Adapt the coordinated fixture signature without retrying invocation."""

    try:
        signature = inspect.signature(function)
    except (TypeError, ValueError) as exc:
        raise RuntimeError("DAG-fixture function signature is unavailable") from exc

    parameters = list(signature.parameters.values())
    if any(item.kind == inspect.Parameter.VAR_POSITIONAL for item in parameters):
        raise RuntimeError("DAG-fixture function may not use variadic positional inputs")

    prior_fixture: Any = None
    payload = prior_result.get("result")
    if isinstance(payload, Mapping):
        prior_fixture = payload.get("synthetic_fixture")

    values: dict[str, Any] = {}
    consumed_emitted_dag = False
    consumed_prior_result = False
    for item in parameters:
        if item.kind in (inspect.Parameter.VAR_KEYWORD,):
            continue
        if item.kind == inspect.Parameter.POSITIONAL_ONLY:
            raise RuntimeError("DAG-fixture function must expose keyword parameters")
        name = item.name.lower()
        if name in {
            "root", "repo_root", "repository_root", "workspace", "directory",
            "source_root", "root_dir",
        }:
            values[item.name] = root
        elif name in {
            "source", "manifest", "checkpoint", "source_manifest", "ledger",
            "source_ledger_options", "ledger_options", "options",
        }:
            values[item.name] = source
        elif name in {"hashes", "input_hashes", "source_hashes"}:
            values[item.name] = hashes
        elif name in {
            "emitted_dag", "candidate_dag", "emitted_candidate_dag", "dag",
            "candidate", "emitted_candidate", "prior_emitted_dag",
        } or ("dag" in name and ("emitted" in name or "candidate" in name)):
            values[item.name] = emitted_dag
            consumed_emitted_dag = True
        elif name in {
            "prior_result", "candidate_result", "previous_result", "v2_result",
            "prior", "prior_receipt",
        }:
            values[item.name] = prior_result
            consumed_prior_result = True
        elif name in {
            "prior_fixture", "prior_synthetic_fixture", "emitted_fixture",
        }:
            require(isinstance(prior_fixture, Mapping),
                    "v2 result synthetic fixture payload is missing")
            values[item.name] = prior_fixture
            consumed_prior_result = True
        elif name in {"prior_result_sha256", "candidate_result_sha256", "result_sha256"}:
            values[item.name] = PRIOR_RESULT_SHA256
            consumed_prior_result = True
        elif item.default is not inspect.Parameter.empty:
            continue
        else:
            raise RuntimeError("unsupported required DAG-fixture argument: " + item.name)

    require(consumed_emitted_dag or consumed_prior_result,
            "synthetic fixture did not declare the previous emitted DAG/result input")
    return function(**values)


def _require_fixture_output(value: Any) -> None:
    require(isinstance(value, Mapping), "emitted-DAG synthetic fixture did not return an object")
    markers = [name for name in ("passed", "checks_passed") if name in value]
    require(bool(markers), "emitted-DAG synthetic fixture omitted a pass marker")
    require(all(value[name] is True for name in markers),
            "emitted-DAG synthetic fixture did not pass")
    for flag in (
        "source_values_computed",
        "source_graph_traversed",
        "certified",
        "numeric_use_admissible",
    ):
        require(flag in value and value[flag] is False,
                "synthetic fixture scope flag is not false: " + flag)
    # Reject non-JSON values before constructing the success receipt.
    _json_bytes(value)


def _recordable_fixture(value: Any) -> Any:
    """Keep a returned fixture in a failure receipt only if it is JSON-safe."""

    if value is None:
        return None
    try:
        _json_bytes(value)
    except Exception:
        return None
    return value


def main() -> int:
    source, hashes, source_sha = preflight()
    entry_sha = write_once(
        ENTRY,
        {
            "schema": f"{COMPONENT}_entered_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "input_hashes": hashes,
            "consumer_calls_started": 0,
            "verification_invocations": 0,
            "ancestor_math_calls": 0,
            "source_values_computed": False,
            "source_graph_traversed": False,
            "certified": False,
            "numeric_use_admissible": False,
            "runtime_model_metadata": None,
        },
    )

    calls: list[str] = []
    invocation_count = 0
    fixture_result: Any = None
    try:
        # This is the first content read of the prior result.  Its raw bytes
        # were authenticated before entry as required by the manifest gate.
        prior_result, emitted_dag = _read_prior_result_after_entry()
        module = _load_new_module()
        fixture = _find_fixture(module)
        calls.append("emitted_dag_synthetic_fixture")
        invocation_count += 1
        fixture_result = _call_fixture(
            fixture,
            root=ROOT,
            source=source,
            hashes=hashes,
            prior_result=prior_result,
            emitted_dag=emitted_dag,
        )
        _require_fixture_output(fixture_result)

        output: dict[str, Any] = {
            "schema": f"{COMPONENT}_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "source_checkpoint_sha256": source_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "consumer_calls": calls,
            "verification_invocations": invocation_count,
            "ancestor_math_calls": 0,
            "source_values_computed": False,
            "source_graph_traversed": False,
            "certified": False,
            "numeric_use_admissible": False,
            "runtime_model_metadata": None,
            "result": {
                "prior_result_sha256": PRIOR_RESULT_SHA256,
                "prior_result_component": prior_result.get("component"),
                "prior_result_status": prior_result.get("status"),
                "emitted_dag": emitted_dag,
                "synthetic_fixture": fixture_result,
            },
            "error": None,
        }
        # A successful result must be serializable before the exclusive write;
        # otherwise the guarded path records a failure receipt instead.
        _json_bytes(output)
    except Exception as exc:
        output = {
            "schema": f"{COMPONENT}_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "source_checkpoint_sha256": source_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "consumer_calls": calls,
            "verification_invocations": invocation_count,
            "ancestor_math_calls": 0,
            "source_values_computed": False,
            "source_graph_traversed": False,
            "certified": False,
            "numeric_use_admissible": False,
            "runtime_model_metadata": None,
            "fixture_result": _recordable_fixture(fixture_result),
            "error": {"type": type(exc).__name__, "message": str(exc)},
            "traceback": traceback.format_exc(),
        }

    result_sha = write_once(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "result_sha256": result_sha,
                "consumer_calls": calls,
                "verification_invocations": invocation_count,
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
