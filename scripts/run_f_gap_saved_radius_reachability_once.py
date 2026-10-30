"""Run the saved-radius reachability and frontier-role consumers once.

This is a persistence wrapper, not a graph producer.  It authenticates the
raw manifest entries before the entry receipt is written, then imports only
the two new consumer leaves.  The reachability leaf owns the bounded walk of
the already-saved fresh namespace; this wrapper never decodes an archive,
imports an old generator, rebuilds an ancestor graph, or computes source
values.

The frontier leaf's four required prior result paths plus its separately
authenticated optional logRp bridge are deliberately taken from the same
frozen manifest.  The manifest therefore carries five distinct ``result.json``
files; the bridge remains optional to the frontier role semantics, and this
runner does not infer equivalence from a shared numeric ID.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path
from typing import Any, Callable, Mapping


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "f_gap_saved_radius_reachability_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

RUNNER = Path("scripts/run_f_gap_saved_radius_reachability_once.py")
REACHABILITY_LEAF = Path(
    "src/openai_ns_reconstruction/f_gap_saved_radius_reachability.py"
)
FRONTIER_LEAF = Path(
    "src/openai_ns_reconstruction/f_gap_radius_frontier_roles.py"
)
FRESH_RESULT = Path("tmp/frame_gap_fresh_substitution_v1/result.json")
FRESH_AUDIT_RESULT = Path("tmp/frame_gap_fresh_substitution_audit_v1/result.json")
OLD_FRESH_PRODUCER = Path(
    "src/openai_ns_reconstruction/frame_gap_fresh_substitution.py"
)

FRONTIER_RESULT_INPUTS = {
    "tmp/f_gap_radius_spine_binding_v1/result.json":
        "3bbf29e3479632f5140c34931b9f752250ae770a9a1c2093819465c5e822854d",
    "tmp/historical_frame_owner_binding_v1/result.json":
        "8c1acc7746ea6044bfcacc1ca8aa29dd97db83076e795f4af68e6779530b880d",
    "tmp/frame_axis_gamma_roles_v1/result.json":
        "2facb6d999cbbbb386cb8e53945cee571a2a7ab13688bd3c87507aa6cbb67a91",
    "tmp/pulse_gap_current_radius_source_v1/result.json":
        "c6d170f6935e0dad3937e155fbbb7e67c6590ee19369291594d10479883fef05",
}
FRONTIER_OPTIONAL_RESULT_INPUTS = {
    # Required in this component manifest and raw-authenticated when present;
    # the frontier leaf treats this fifth packet as an optional conditional
    # witness rather than a required role-equivalence premise.
    "tmp/original_native_logRp_source_bridge_v1/result.json":
        "41df6f56fc4a071de606372261a01425c445cd98962c76dbb46164e9699b621c",
}

BASE_INPUTS = (
    REACHABILITY_LEAF.as_posix(),
    RUNNER.as_posix(),
    FRESH_RESULT.as_posix(),
    FRESH_AUDIT_RESULT.as_posix(),
    OLD_FRESH_PRODUCER.as_posix(),
)
FRONTIER_INPUTS = tuple(FRONTIER_RESULT_INPUTS)
FRONTIER_INPUTS += tuple(FRONTIER_OPTIONAL_RESULT_INPUTS)
REQUIRED_INPUTS = BASE_INPUTS + (FRONTIER_LEAF.as_posix(),) + FRONTIER_INPUTS

EXPECTED_FIXED_SHA256 = {
    FRESH_RESULT.as_posix(): (
        "a37fc375ff06c06d39b4250ffa3d36c6c315182a81e5d783f4bfb1f981ba0234"
    ),
    FRESH_AUDIT_RESULT.as_posix(): (
        "a3de519abfdc42fe18a0fa6830105c9c2c39a7bb4188d27693123117123087a3"
    ),
    OLD_FRESH_PRODUCER.as_posix(): (
        "8b50e086ddfbf09b3ed610663c7577ef9939925b43eb4babfc46d37a545742ba"
    ),
    **FRONTIER_RESULT_INPUTS,
    **FRONTIER_OPTIONAL_RESULT_INPUTS,
}

REACHABILITY_BUILDER_NAMES = (
    # The first name is the canonical interface for this component.  The two
    # aliases keep this wrapper compatible with the leaf currently staged in
    # the shared worktree without hiding which name was actually selected.
    "build_radius_reachability",
    "build_f_gap_saved_radius_reachability",
    "build_saved_radius_reachability",
)
FRONTIER_BUILDER_NAMES = (
    "build_frontier_roles",
    "build_f_gap_radius_frontier_roles",
    "build_saved_frontier_roles",
)

SCOPE_FALSE_KEYS = (
    "function_identity_proved",
    "capture_free_substitution_proved",
    "capture_free_proved",
    "capture_proved",
    "ordinary_Z_transport_proved",
    "ordinary_z_transport_proved",
    "ordinary_Z_proved",
    "m15_closure_proved",
    "closure_proved",
    "source_values_computed",
    "source_values_selected",
    "source_value_selected",
    "actual_source_values_computed",
    "numeric_values_computed",
    "numeric_use_admissible",
    "old_source_graph_indexed",
    "old_source_graph_traversed",
    "old_graph_rebuilt",
    "ancestor_replayed",
    "protected_root_indexed",
    "protected_frontier_dereferenced",
    "protected_frontier_indexed",
    "protected_ref_dereferenced",
    "saved_boundary_dereferenced",
    "used_as_numeric_value",
)


class RunnerInputError(ValueError):
    """Raised when the frozen one-shot contract is not satisfied."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
        + b"\n"
    )


def _write_once(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as handle:
        handle.write(_json_bytes(value))


def _json_safe(value: Any) -> Any:
    """Retain a returned packet even when a leaf used a non-JSON scalar."""

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted((_json_safe(item) for item in value), key=repr)
    if isinstance(value, Path):
        return value.as_posix()
    return repr(value)


def _relative_path(value: Any, label: str) -> str:
    _require(isinstance(value, str) and value, label + " path is missing")
    normalized = value.replace("\\", "/")
    path = Path(normalized)
    _require(not path.is_absolute(), label + " must be relative")
    _require(".." not in path.parts, label + " escapes repository")
    return normalized


def _safe_repository_path(relative: str) -> Path:
    path = (ROOT / relative).resolve()
    _require(path.is_relative_to(ROOT.resolve()), "manifest path escapes repository")
    _require(path.is_file(), "manifest input is missing: " + relative)
    return path


def _sha_value(value: Any, label: str) -> str:
    _require(isinstance(value, str) and len(value) == 64, label + " SHA malformed")
    value = value.lower()
    _require(
        all(character in "0123456789abcdef" for character in value),
        label + " SHA malformed",
    )
    return value


def _manifest_entries(source: Mapping[str, Any]) -> dict[str, str]:
    entries: dict[str, str] = {}
    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(isinstance(item, Mapping), f"manifest files[{index}] malformed")
            name = _relative_path(item.get("path"), f"manifest files[{index}]")
            digest = _sha_value(item.get("sha256"), f"manifest files[{index}]")
            _require(name not in entries, "manifest duplicate path: " + name)
            entries[name] = digest
    for field in ("sha256", "input_hashes", "hashes"):
        table = source.get(field)
        if not isinstance(table, Mapping):
            continue
        for raw_name, raw_digest in table.items():
            name = _relative_path(str(raw_name), f"manifest {field}")
            digest = _sha_value(raw_digest, f"manifest {field}.{name}")
            _require(
                name not in entries or entries[name] == digest,
                "manifest hash conflict: " + name,
            )
            entries[name] = digest
    _require(entries, "manifest contains no raw SHA entries")
    return entries


def _frontier_result_paths(entries: Mapping[str, str]) -> tuple[str, ...]:
    base = set(BASE_INPUTS)
    _require(FRONTIER_LEAF.as_posix() in entries, "frontier leaf missing from manifest")
    ordered = tuple(FRONTIER_INPUTS)
    _require(
        all(path in entries for path in ordered),
        "manifest is missing one of the five frontier prior results",
    )
    _require(
        set(entries) == base | {FRONTIER_LEAF.as_posix()} | set(ordered),
        "manifest has unexpected raw inputs",
    )
    return ordered


def _required_input_paths(source: Mapping[str, Any]) -> tuple[str, ...]:
    entries = _manifest_entries(source)
    frontier_paths = _frontier_result_paths(entries)
    required = tuple(BASE_INPUTS) + (FRONTIER_LEAF.as_posix(),) + frontier_paths
    _require(set(entries) == set(required), "manifest has unexpected raw inputs")
    return required


def _authenticate_raw_inputs(
    source: Mapping[str, Any], required: tuple[str, ...]
) -> tuple[dict[str, str], dict[str, str]]:
    entries = _manifest_entries(source)
    _require(set(entries) == set(required), "manifest input set changed")
    actual: dict[str, str] = {}
    for name in required:
        expected = entries[name]
        observed = _sha256_file(_safe_repository_path(name))
        _require(observed == expected, "raw SHA changed: " + name)
        fixed = EXPECTED_FIXED_SHA256.get(name)
        if fixed is not None:
            _require(expected == fixed, "fixed SHA pin changed: " + name)
        actual[name] = observed
    return dict(entries), actual


def _checkpoint_sha(brief: Mapping[str, Any]) -> None:
    checkpoint = brief.get("source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "brief source_checkpoint missing")
    raw_path = checkpoint.get("path")
    _require(
        isinstance(raw_path, str) and raw_path,
        "brief source_checkpoint path is missing",
    )
    normalized = raw_path.replace("\\", "/")
    candidate = Path(normalized)
    is_windows_absolute = len(normalized) >= 2 and normalized[1] == ":"
    if candidate.is_absolute() or is_windows_absolute:
        resolved = candidate.resolve()
    else:
        relative = _relative_path(raw_path, "brief source_checkpoint")
        resolved = (ROOT / relative).resolve()
    _require(
        resolved == SOURCE.resolve(),
        "brief checkpoint path changed",
    )
    expected = _sha_value(
        checkpoint.get("sha256"), "brief source_checkpoint"
    )
    _require(expected == _sha256_file(SOURCE), "source checkpoint SHA changed")


def _require_guard(
    source: Mapping[str, Any], names: tuple[str, ...], label: str
) -> None:
    present = [name for name in names if name in source]
    _require(present, "manifest guard missing: " + label)
    _require(any(source[name] is True for name in present), "manifest guard is not true: " + label)


def _load_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not valid JSON") from exc
    _require(isinstance(value, dict), label + " is not a JSON object")
    return value


def _preflight() -> tuple[dict[str, Any], tuple[str, ...], dict[str, str], dict[str, str]]:
    _require(not ENTRY.exists(), "entry receipt already exists")
    _require(not RESULT.exists(), "result receipt already exists")
    _require(BRIEF.is_file(), "pending brief is missing")
    _require(SOURCE.is_file(), "source manifest is missing")
    brief = _load_json_object(BRIEF, "pending brief")
    source = _load_json_object(SOURCE, "source manifest")
    _require(brief.get("component") == COMPONENT, "brief component changed")
    _require(brief.get("phase") == PHASE, "brief phase changed")
    _require(brief.get("goal_active") is True, "brief goal is not active")
    _require(brief.get("full_goal_complete") is False, "brief goal is already complete")
    _require(brief.get("do_not_reenter") is True, "brief re-entry guard missing")
    _checkpoint_sha(brief)
    _require(source.get("component") == COMPONENT, "manifest component changed")
    _require(source.get("one_entry_only") is True, "manifest one-entry guard missing")
    _require(
        source.get("all_raw_inputs_authenticated_before_decompression") is True,
        "manifest raw-authentication guard missing",
    )
    _require_guard(
        source,
        ("no_source_graph_traversal", "no_old_source_graph_traversal"),
        "old source graph traversal",
    )
    _require_guard(
        source,
        ("no_archive_decode", "no_ancestor_archive_decompression"),
        "archive/decompression",
    )
    _require_guard(
        source,
        ("no_old_generator_import", "no_old_producer_import", "no_ancestor_source_import"),
        "old generator/ancestor import",
    )
    required = _required_input_paths(source)
    declared, actual = _authenticate_raw_inputs(source, required)
    return source, required, declared, actual


def _import_leaf(relative: Path, module_name: str) -> Any:
    path = (ROOT / relative).resolve()
    _require(path.is_file(), "new leaf is missing: " + relative.as_posix())
    spec = importlib.util.spec_from_file_location(module_name, path)
    _require(spec is not None and spec.loader is not None, "cannot load new leaf")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _callable(module: Any, names: tuple[str, ...], label: str) -> tuple[str, Callable[..., Any]]:
    for name in names:
        candidate = getattr(module, name, None)
        if callable(candidate):
            return name, candidate
    raise RunnerInputError(label + " callable is missing")


def _invoke(
    steps: list[dict[str, Any]],
    name: str,
    function: Callable[..., Any],
    *args: Any,
) -> Any:
    record: dict[str, Any] = {"step": name, "invoked": True}
    try:
        returned = function(*args)
    except BaseException as exc:
        record["exception"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
        steps.append(record)
        raise
    record["returned_packet"] = _json_safe(returned)
    steps.append(record)
    return returned


def _require_mapping(value: Any, label: str) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), label + " did not return a mapping")
    return value


def _status_passed(packet: Mapping[str, Any], label: str) -> None:
    status = str(packet.get("status", "")).lower()
    checks = packet.get("checks")
    checks_pass = isinstance(checks, Mapping) and bool(checks) and all(
        item is True for item in checks.values()
    )
    passed = (
        packet.get("passed") is True
        or packet.get("checks_passed") is True
        or checks_pass
        or "passed" in status
    )
    _require(passed, label + " synthetic fixture did not pass")


def _reject_positive_scope_flags(packet: Mapping[str, Any], label: str) -> None:
    for key in SCOPE_FALSE_KEYS:
        if packet.get(key) is True:
            raise RunnerInputError(label + " set forbidden flag true: " + key)
    for container_name in (
        "scope",
        "scope_flags",
        "source_gates",
        "flags",
        "protected_frontier_17",
        "protected_ref",
        "protected_refs",
    ):
        container = packet.get(container_name)
        if not isinstance(container, Mapping):
            continue
        for key in SCOPE_FALSE_KEYS:
            if container.get(key) is True:
                raise RunnerInputError(
                    label + " set forbidden " + container_name + " flag true: " + key
                )
        for key in (
            "indexed",
            "dereferenced",
            "source_value_selected",
            "used_as_numeric_value",
        ):
            if container.get(key) is True:
                raise RunnerInputError(
                    label + " set forbidden " + container_name + " flag true: " + key
                )


def _validate_fixture(packet: Any, label: str) -> Mapping[str, Any]:
    body = _require_mapping(packet, label)
    _status_passed(body, label)
    _reject_positive_scope_flags(body, label)
    if "actual_source_read" in body:
        _require(body["actual_source_read"] is False, label + " read actual source")
    if "source_graph_indexed" in body:
        _require(body["source_graph_indexed"] is False, label + " indexed old source graph")
    return body


def _validate_common_build_scope(packet: Any, label: str) -> Mapping[str, Any]:
    body = _require_mapping(packet, label)
    _reject_positive_scope_flags(body, label)
    if "actual_source_values_computed" in body:
        _require(
            body["actual_source_values_computed"] is False,
            label + " computed actual source values",
        )
    if "source_graph_traversed" in body:
        _require(body["source_graph_traversed"] is False, label + " traversed old source graph")
    if "old_source_graph_traversed" in body:
        _require(body["old_source_graph_traversed"] is False, label + " traversed old source graph")
    return body


def _reachability_metrics(packet: Mapping[str, Any]) -> dict[str, Any]:
    scope = packet.get("scope")
    _require(isinstance(scope, Mapping), "reachability result scope is missing")
    _require(scope.get("only_saved_fresh_nodes_read") is True, "reachability left fresh scope")
    _require(
        scope.get("only_fresh_namespace_operand_edges_followed") is True,
        "reachability followed a non-fresh edge",
    )
    _require(scope.get("saved_boundaries_never_dereferenced") is True, "saved ref was dereferenced")
    _require(scope.get("protected_refs_never_dereferenced") is True, "protected ref was dereferenced")
    _require(scope.get("old_source_graph_indexed") is False, "old source graph was indexed")
    _require(scope.get("old_graph_rebuilt") is False, "old graph was rebuilt")
    _require(scope.get("source_values_computed") is False, "source values were computed")
    _require(packet.get("root_count") == 20, "reachability root count changed")
    roots = packet.get("roots")
    _require(isinstance(roots, list) and len(roots) == 20, "reachability roots are not twenty rows")
    counts: list[dict[str, Any]] = []
    for index, root in enumerate(roots):
        _require(isinstance(root, Mapping), f"reachability root {index} malformed")
        count = root.get("visited_node_count")
        _require(type(count) is int and count >= 0, f"reachability root {index} visited count missing")
        counts.append(
            {
                "role": root.get("role"),
                "row_order": root.get("row_order"),
                "visited_node_count": count,
            }
        )
    gates = packet.get("source_gates")
    _require(isinstance(gates, Mapping), "reachability source gates are missing")
    for key in (
        "function_identity_proved",
        "capture_free_substitution_proved",
        "ordinary_Z_transport_proved",
        "m15_closure_proved",
    ):
        _require(gates.get(key) is False, "reachability gate changed: " + key)
    _require(packet.get("target_saved_node") == 33982, "saved boundary target changed")
    return {
        "fresh_saved_DAG_analyzed": True,
        "fresh_saved_DAG_traversed": True,
        "fresh_namespace": packet.get("fresh_namespace"),
        "saved_namespace": packet.get("saved_namespace"),
        "root_count": len(roots),
        "visited_node_count_total": sum(item["visited_node_count"] for item in counts),
        "visited_node_counts_by_root": counts,
        "prior_fresh_node_count": packet.get("prior_fresh_node_count"),
        "saved_boundary_terminations": packet.get("saved_boundary_terminations"),
        "protected_terminations": packet.get("protected_terminations"),
        "target_saved_node": packet.get("target_saved_node"),
        "saved_refs_terminated": True,
        "protected_refs_never_dereferenced": True,
    }


def _frontier_metrics(packet: Mapping[str, Any]) -> dict[str, Any]:
    scope = packet.get("scope")
    if not isinstance(scope, Mapping):
        scope = packet.get("scope_flags")
    if isinstance(scope, Mapping):
        _require(scope.get("old_source_graph_indexed") is not True, "frontier indexed old source graph")
        _require(scope.get("old_source_graph_traversed") is not True, "frontier traversed old source graph")
    visited = packet.get("visited_node_count")
    if type(visited) is not int:
        visited = None
    traversed = (
        packet.get("fresh_saved_DAG_traversed") is True
        or packet.get("fresh_saved_DAG_analyzed") is True
        or (isinstance(scope, Mapping) and scope.get("only_saved_fresh_nodes_read") is True)
        or visited is not None
    )
    return {
        "fresh_saved_DAG_analyzed": bool(
            packet.get("fresh_saved_DAG_analyzed") is True or traversed
        ),
        "fresh_saved_DAG_traversed": bool(traversed),
        "visited_node_count": visited,
        "scope": _json_safe(scope) if isinstance(scope, Mapping) else None,
    }


def _base_flags(
    reach_metrics: Mapping[str, Any] | None,
    frontier_metrics: Mapping[str, Any] | None,
) -> dict[str, Any]:
    reach_traversed = bool(reach_metrics and reach_metrics.get("fresh_saved_DAG_traversed"))
    frontier_traversed = bool(
        frontier_metrics and frontier_metrics.get("fresh_saved_DAG_traversed")
    )
    return {
        # These two names intentionally refer to the old/original source graph
        # only.  The separate fresh_saved_DAG_* keys report the saved fresh walk.
        "source_graph_traversed": False,
        "old_source_graph_traversed": False,
        "source_graph_traversed_scope": "old/original source graph only",
        "fresh_saved_DAG_analyzed": reach_traversed,
        "fresh_saved_DAG_traversed": reach_traversed,
        "frontier_saved_DAG_analyzed": bool(
            frontier_metrics and frontier_metrics.get("fresh_saved_DAG_analyzed")
        ),
        "frontier_saved_DAG_traversed": frontier_traversed,
        "actual_source_values_computed": False,
        "source_values_computed": False,
        "numeric_values_computed": False,
        "numeric_use_admissible": False,
        "function_identity_proved": False,
        "capture_proved": False,
        "capture_free_proved": False,
        "capture_free_substitution_proved": False,
        "ordinary_Z_proved": False,
        "ordinary_Z_transport_proved": False,
        "closure_proved": False,
        "m15_closure_proved": False,
        "certified": False,
        "full_goal_complete": False,
    }


def _envelope(
    *,
    status: str,
    source: Mapping[str, Any] | None,
    required: tuple[str, ...],
    declared: Mapping[str, str],
    actual: Mapping[str, str],
    steps: list[dict[str, Any]],
    reach_packet: Any = None,
    frontier_packet: Any = None,
    reach_metrics: Mapping[str, Any] | None = None,
    frontier_metrics: Mapping[str, Any] | None = None,
    error: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    flags = _base_flags(reach_metrics, frontier_metrics)
    reach_safe = _json_safe(reach_packet) if reach_packet is not None else None
    frontier_safe = _json_safe(frontier_packet) if frontier_packet is not None else None
    result: dict[str, Any] = {
        "schema": COMPONENT + "_result_v1",
        "component": COMPONENT,
        "phase": PHASE,
        "status": status,
        "flags": flags,
        **flags,
        "scope": {
            "old_source_graph_traversed": False,
            "source_graph_traversed": False,
            "fresh_saved_DAG_traversed": flags["fresh_saved_DAG_traversed"],
            "fresh_saved_DAG_analyzed": flags["fresh_saved_DAG_analyzed"],
            "reachability": _json_safe(reach_metrics),
            "frontier_roles": _json_safe(frontier_metrics),
        },
        "required_inputs": list(required),
        "input_hashes": dict(declared),
        "observed_input_hashes": dict(actual),
        "raw_sha256_authenticated_before_entry": bool(actual),
        "old_fresh_producer_schema_evidence": {
            "path": OLD_FRESH_PRODUCER.as_posix(),
            "sha256": declared.get(OLD_FRESH_PRODUCER.as_posix()),
            "raw_hashed": OLD_FRESH_PRODUCER.as_posix() in actual,
            "imported": False,
            "executed": False,
            "decoded": False,
        },
        "entry_path": ENTRY.relative_to(ROOT).as_posix() if ENTRY.exists() else None,
        "entry_sha256": _sha256_file(ENTRY) if ENTRY.is_file() else None,
        "steps": steps,
        "reachability_result": reach_safe,
        "frontier_roles_result": frontier_safe,
    }
    if source is not None:
        result["source_manifest"] = _json_safe(source)
    if error is not None:
        result["error"] = _json_safe(error)
    return result


def main() -> int:
    # Preflight stays outside the analysis try/except, matching the reference
    # one-shot runners: an untrusted/missing manifest cannot create a false
    # result receipt before entry.  Once ENTRY exists, every returned packet
    # and every exception is retained in the single RESULT receipt.
    source, required, declared, actual = _preflight()
    steps: list[dict[str, Any]] = []
    reach_packet: Any = None
    frontier_packet: Any = None
    reach_metrics: dict[str, Any] | None = None
    frontier_metrics: dict[str, Any] | None = None
    _write_once(
        ENTRY,
        {
            "schema": COMPONENT + "_entry_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "source": SOURCE.relative_to(ROOT).as_posix(),
            "source_sha256": _sha256_file(SOURCE),
            "required_inputs": list(required),
            "input_hashes": declared,
            "raw_sha256_authenticated_before_entry": True,
            "old_source_graph_traversed": False,
            "archive_decoded": False,
            "old_generator_imported": False,
        },
    )
    try:
        reach_module = _import_leaf(
            REACHABILITY_LEAF, "_f_gap_saved_radius_reachability_new_leaf"
        )
        fixture_name, fixture = _callable(
            reach_module, ("verify_synthetic_fixture",), "reachability fixture"
        )
        reach_fixture = _invoke(steps, "reachability." + fixture_name, fixture)
        _validate_fixture(reach_fixture, "reachability")
        builder_name, builder = _callable(
            reach_module, REACHABILITY_BUILDER_NAMES, "reachability builder"
        )
        reach_packet = _invoke(
            steps,
            "reachability." + builder_name,
            builder,
            ROOT,
            source,
            declared,
        )
        reach_body = _validate_common_build_scope(reach_packet, "reachability")
        reach_metrics = _reachability_metrics(reach_body)

        frontier_module = _import_leaf(
            FRONTIER_LEAF, "_f_gap_radius_frontier_roles_new_leaf"
        )
        frontier_fixture_name, frontier_fixture = _callable(
            frontier_module, ("verify_synthetic_fixture",), "frontier fixture"
        )
        frontier_fixture_packet = _invoke(
            steps,
            "frontier." + frontier_fixture_name,
            frontier_fixture,
        )
        _validate_fixture(frontier_fixture_packet, "frontier")
        frontier_builder_name, frontier_builder = _callable(
            frontier_module, FRONTIER_BUILDER_NAMES, "frontier builder"
        )
        frontier_packet = _invoke(
            steps,
            "frontier." + frontier_builder_name,
            frontier_builder,
            ROOT,
            source,
            declared,
        )
        frontier_body = _validate_common_build_scope(frontier_packet, "frontier")
        frontier_metrics = _frontier_metrics(frontier_body)
        output = _envelope(
            status="returned",
            source=source,
            required=required,
            declared=declared,
            actual=actual,
            steps=steps,
            reach_packet=reach_packet,
            frontier_packet=frontier_packet,
            reach_metrics=reach_metrics,
            frontier_metrics=frontier_metrics,
        )
        _write_once(RESULT, output)
        print(
            json.dumps(
                {
                    "status": output["status"],
                    "result": RESULT.relative_to(ROOT).as_posix(),
                    "fresh_saved_DAG_traversed": output["fresh_saved_DAG_traversed"],
                    "old_source_graph_traversed": output["old_source_graph_traversed"],
                    "reachability_visited_node_count": (
                        reach_metrics or {}
                    ).get("visited_node_count_total"),
                },
                sort_keys=True,
            )
        )
        return 0
    except BaseException as exc:
        error = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
        output = _envelope(
            status="failed",
            source=source,
            required=required,
            declared=declared,
            actual=actual,
            steps=steps,
            reach_packet=reach_packet,
            frontier_packet=frontier_packet,
            reach_metrics=reach_metrics,
            frontier_metrics=frontier_metrics,
            error=error,
        )
        _write_once(RESULT, output)
        print(
            json.dumps(
                {
                    "status": output["status"],
                    "result": RESULT.relative_to(ROOT).as_posix(),
                    "error_type": error["type"],
                    "fresh_saved_DAG_traversed": output["fresh_saved_DAG_traversed"],
                    "old_source_graph_traversed": output["old_source_graph_traversed"],
                },
                sort_keys=True,
            )
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
