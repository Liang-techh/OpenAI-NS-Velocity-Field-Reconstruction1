"""Enter the F gap-endpoint candidate leaf exactly once.

This runner is deliberately source-only.  It authenticates every raw file in
``tmp/f_gap_endpoint_candidate_v1/source.json`` before the exclusive entry
receipt, imports only the new candidate leaf after entry, and writes a complete
success or failure result before printing a summary.  It never imports an old
producer, opens a saved source graph, or decompresses an archive itself.

The new leaf may expose ``verify_synthetic_fixture`` (or
``synthetic_fixture``/``run_fixture``) and one of
``build_f_gap_endpoint_candidate``, ``build_source_role_candidate``, or
``build_gap_endpoint_candidate``.  The selected functions are called at most
once each.  The small signature adapter is only for the caller-owned source
ledger interface; it does not retry a failed call.
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
COMPONENT = "f_gap_endpoint_candidate_v1"
PHASE = "f_gap_endpoint_candidate_v1_pending"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"
LEAF = ROOT / "src/openai_ns_reconstruction/f_gap_endpoint_candidate.py"

REQUIRED_INPUTS = {
    "src/openai_ns_reconstruction/f_gap_endpoint_candidate.py",
    "scripts/run_f_gap_endpoint_candidate_once.py",
    "tmp/pulse_gap_actual_source_program_v1/result.json",
    "tmp/pf_noncommon_moment_v1/result.json",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions.py",
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


def write_once(path: Path, value: Mapping[str, Any]) -> str:
    raw = (
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
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
    """Authenticate either the normal ``files`` list or a hash mapping."""

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
        require(isinstance(expected, str) and len(expected) == 64
                and all(c in "0123456789abcdefABCDEF" for c in expected),
                "invalid SHA-256 for " + key)
        actual = sha(path)
        require(actual.lower() == expected.lower(), "raw SHA mismatch: " + key)
        hashes[key] = actual

    required = {item.lower() for item in REQUIRED_INPUTS}
    present = {item.lower() for item in hashes}
    missing = sorted(item for item in required if item not in present)
    require(not missing, "required frozen inputs missing: " + ", ".join(missing))
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
    require(checkpoint_path.resolve() == SOURCE.resolve(),
            "pending brief source checkpoint path mismatch")
    expected = checkpoint.get("sha256")
    require(isinstance(expected, str) and expected.lower() == sha(SOURCE),
            "pending brief source checkpoint SHA mismatch")
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
    require(brief.get("goal_active") is True
            and brief.get("full_goal_complete") is False
            and brief.get("do_not_reenter") is True,
            "pending brief goal flags mismatch")
    source_sha = _checkpoint_sha(brief)
    require(source.get("component") == COMPONENT
            and source.get("one_entry_only") is True
            and source.get("no_source_graph_traversal") is True
            and source.get("all_raw_inputs_authenticated_before_decompression") is True,
            "source manifest safeguards missing")
    hashes = _manifest_hashes(source)
    require(LEAF.is_file(), "new candidate leaf is missing")
    return dict(source), hashes, source_sha


def _load_new_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "f_gap_endpoint_candidate_new_leaf", LEAF
    )
    require(spec is not None and spec.loader is not None,
            "new candidate leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _call_with_declared_inputs(
    function: Callable[..., Any],
    *,
    root: Path,
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
) -> Any:
    """Call one new-leaf function without retrying after an invocation."""

    try:
        signature = inspect.signature(function)
    except (TypeError, ValueError) as exc:
        raise RuntimeError("candidate function signature is unavailable") from exc

    parameters = list(signature.parameters.values())
    if any(item.kind == inspect.Parameter.VAR_POSITIONAL for item in parameters):
        raise RuntimeError("candidate function may not use variadic positional inputs")

    values: dict[str, Any] = {}
    for item in parameters:
        if item.kind in (inspect.Parameter.VAR_KEYWORD,):
            continue
        name = item.name.lower()
        if name in {
            "root",
            "repo_root",
            "repository_root",
            "workspace",
            "directory",
            "source_root",
            "root_dir",
        }:
            values[item.name] = root
        elif name in {
            "source",
            "manifest",
            "checkpoint",
            "source_manifest",
            "ledger",
            "source_ledger_options",
            "ledger_options",
            "options",
        }:
            values[item.name] = source
        elif name in {"hashes", "input_hashes", "source_hashes"}:
            values[item.name] = hashes
        elif item.default is not inspect.Parameter.empty:
            continue
        else:
            raise RuntimeError("unsupported required candidate argument: " + item.name)

    return function(**values)


def _find_function(module: Any, names: tuple[str, ...], label: str) -> Callable[..., Any]:
    for name in names:
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate
    raise RuntimeError("new candidate leaf has no " + label + " function")


def _require_fixture_passed(value: Any) -> None:
    require(isinstance(value, Mapping), "synthetic fixture did not return an object")
    markers = [name for name in ("passed", "checks_passed") if name in value]
    require(bool(markers), "synthetic fixture omitted a pass marker")
    require(all(value[name] is True for name in markers),
            "synthetic fixture did not pass")


def _require_source_only(value: Any) -> None:
    require(isinstance(value, Mapping), "source-role candidate did not return an object")
    require(value.get("source_only") is True,
            "source-role candidate is not explicitly source-only")
    require(value.get("source_values_computed") is False,
            "source-role candidate computed source values")
    require(value.get("certified") is False,
            "source-role candidate claims certification")

    graph_markers = [
        name for name in ("source_graph_traversed", "graph_traversed")
        if name in value
    ]
    require(bool(graph_markers),
            "source-role candidate omitted graph scope flag")
    require(all(value[name] is False for name in graph_markers),
            "source-role candidate traversed a source graph")

    numeric_markers = [
        name for name in ("numeric_use_admissible", "admissible_for_numeric_use")
        if name in value
    ]
    require(bool(numeric_markers),
            "source-role candidate omitted numeric-use scope flag")
    require(all(value[name] is False for name in numeric_markers),
            "source-role candidate is admissible for numeric use")

    false_flags = (
        "source_values_computed",
        "source_value_computation",
        "actual_source_value_computed",
        "source_graph_traversed",
        "graph_traversed",
        "numeric_use_admissible",
        "admissible_for_numeric_use",
        "certified",
    )
    for flag in false_flags:
        if flag in value:
            require(value[flag] is False, "source-only flag is true: " + flag)


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
            "ancestor_math_calls": 0,
            "source_values_computed": False,
            "source_graph_traversed": False,
            "runtime_model_metadata": None,
        },
    )

    calls: list[str] = []
    invocation_count = 0
    try:
        module = _load_new_module()
        synthetic = _find_function(
            module,
            ("verify_synthetic_fixture", "synthetic_fixture", "run_fixture"),
            "synthetic fixture",
        )
        calls.append("synthetic_fixture")
        invocation_count += 1
        synthetic_result = _call_with_declared_inputs(
            synthetic, root=ROOT, source=source, hashes=hashes
        )
        _require_fixture_passed(synthetic_result)

        actual = _find_function(
            module,
            (
                "build_f_gap_endpoint_source_candidate",
                "build_f_gap_endpoint_candidate",
                "build_source_role_candidate",
                "build_gap_endpoint_candidate",
            ),
            "source-role candidate",
        )
        calls.append("actual_source_role_candidate")
        invocation_count += 1
        actual_result = _call_with_declared_inputs(
            actual, root=ROOT, source=source, hashes=hashes
        )
        _require_source_only(actual_result)

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
            "runtime_model_metadata": None,
            "result": {
                "synthetic_fixture": synthetic_result,
                "source_role_candidate": actual_result,
            },
            "error": None,
        }
        # Validate the complete success payload while still inside the
        # guarded block.  A non-JSON candidate result must become a persisted
        # failure result rather than escaping at the final exclusive write.
        json.dumps(output, ensure_ascii=False, allow_nan=False)
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
            "runtime_model_metadata": None,
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
