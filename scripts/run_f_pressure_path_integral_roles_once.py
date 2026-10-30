"""Run the bounded F pressure-integral-role consumer exactly once.

This wrapper is the persistence boundary for
``f_pressure_path_integral_roles_v1``.  It authenticates every raw input
named by the pending manifest before the exclusive entry receipt and imports
only the new integral-role leaf after that receipt.  The leaf is permitted to
decode the authenticated F archive as one whole JSON object, then perform the
declared bounded direct row reads only: the five fixed F graph rows and at
most one explicitly declared non-protected product sibling.  It must not
search graph edges, run BFS, replay ancestors, evaluate an integral, or read
source values.

The runner deliberately preserves unresolved and duplicate-match packets.
Those are valid unresolved structural outcomes, not failures and not zeros.
Both success and post-entry failure paths persist a JSON-safe receipt,
including any packet already returned by the leaf.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
import traceback
from pathlib import Path
from typing import Any, Callable, Mapping


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "f_pressure_path_integral_roles_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

RUNNER = Path("scripts/run_f_pressure_path_integral_roles_once.py")
LEAF = Path("src/openai_ns_reconstruction/f_pressure_path_integral_roles.py")
F_ARCHIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
PRIOR_RADIUS_RESULT = Path("tmp/f_pressure_radius_path_v1/result.json")
C3_SOURCE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions.py"
)
PULSE_SOURCE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions.py"
)
PRE_RP_SOURCE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.py"
)

REQUIRED_INPUTS = (
    LEAF.as_posix(),
    RUNNER.as_posix(),
    F_ARCHIVE.as_posix(),
    PRIOR_RADIUS_RESULT.as_posix(),
    C3_SOURCE.as_posix(),
    PULSE_SOURCE.as_posix(),
    PRE_RP_SOURCE.as_posix(),
)

F_ARCHIVE_SHA256 = (
    "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"
)
PRIOR_RADIUS_RESULT_SHA256 = (
    "82674aa39957d87cdbb067e9332a31bceea36cae85e397cf6e708c6dc1949847"
)
EXPECTED_FIXED_SHA256 = {
    F_ARCHIVE.as_posix(): F_ARCHIVE_SHA256,
    PRIOR_RADIUS_RESULT.as_posix(): PRIOR_RADIUS_RESULT_SHA256,
}

SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")

# These handles are structural frontiers.  The runner checks that the leaf
# does not claim to index them; it never reads their rows itself.
PROTECTED_NODE_IDS = frozenset(
    {
        0,
        1,
        17,
        33572,
        25220,
        263344,
        263348,
        263349,
        262918,
        205147,
        263346,
    }
)
TARGET_NODE_IDS = frozenset({33982})
DECLARED_GRAPH_ROW_IDS = frozenset({279383, 279382, 279479, 279478, 34721})
MAX_DIRECT_GRAPH_ROWS = 6
MAX_INTEGRAL_ROLE_MAPPINGS = 36

# A leaf may report structural metadata, but these semantic conclusions are
# outside this component's scope and must never be hard-asserted true.
FORBIDDEN_TRUE_KEYS = frozenset(
    {
        "protected_dereferenced",
        "protected_ref_dereferenced",
        "protected_frontier_dereferenced",
        "protected_indexed",
        "protected_root_indexed",
        "protected_nodes_indexed",
        "target_node_indexed",
        "target_nodes_indexed",
        "target_row_indexed",
        "target_node_dereferenced",
        "target_row_dereferenced",
        "radius_row_indexed",
        "saved_target_dereferenced",
        "source_values_computed",
        "actual_source_values",
        "source_value_selected",
        "source_values_selected",
        "numeric_values_computed",
        "numeric_use_admissible",
        "integrals_evaluated",
        "integral_evaluated",
        "quadrature_evaluated",
        "m15_integral_evaluated",
        "graph_traversal",
        "graph_edge_search",
        "source_graph_traversed",
        "old_graph_traversed",
        "old_F_graph_traversed",
        "bfs_used",
        "ancestor_replayed",
        "old_ancestor_replayed",
        "old_leaf_imported",
        "old_producer_called",
        "duplicate_graph_computation",
        "repeated_graph_computation",
        "source_identity_proved",
        "same_source_function_proved",
        "function_identity_proved",
        "actual_source_identity_proved",
        "capture_free",
        "capture_free_proved",
        "capture_free_substitution_proved",
        "ordinary_Z_proved",
        "ordinary_Z_transport_proved",
        "ordinary_z_transport_proved",
        "endpoint_identity_proved",
        "endpoint_transport_proved",
        "moving_endpoint_rule_proved",
        "closure_proved",
        "m15_closure_proved",
        "five_moment_closure_proved",
        "no_dependence",
        "no_dependency",
        "dependency_absent",
        "radius_not_used",
        "mathematical_nonzero_dependence_proved",
    }
)


class RunnerInputError(ValueError):
    """Raised when the frozen one-shot contract is not satisfied."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha256_file(path: Path) -> str:
    """Hash raw bytes only; this does not decode or import any input."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
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


def _json_safe(value: Any) -> Any:
    """Retain a leaf packet without interpreting source rows or values."""

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


def _write_once(path: Path, value: Any) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = _json_bytes(_json_safe(value))
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _load_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not valid JSON") from exc
    _require(isinstance(value, dict), label + " is not a JSON object")
    return value


def _repository_path(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    normalized = raw_path.replace("\\", "/")
    relative = Path(normalized)
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        label + " must be repository-relative",
    )
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(
        path.is_relative_to(ROOT.resolve()) and path.is_file(),
        label + " is missing or escapes repository: " + normalized,
    )
    return key, path


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    """Authenticate every direct raw path/hash record before archive decode."""

    records: list[tuple[Any, Any, str]] = []
    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(
                isinstance(item, Mapping),
                f"source.files[{index}] is not an object",
            )
            records.append(
                (item.get("path"), item.get("sha256"), f"source.files[{index}]")
            )
    for field in (
        "sha256",
        "input_hashes",
        "input_sha256",
        "raw_input_sha256",
        "code_sha256",
    ):
        table = source.get(field)
        if isinstance(table, Mapping):
            for raw_path, expected in table.items():
                records.append((raw_path, expected, f"source.{field}"))
    _require(records, "source manifest contains no direct raw SHA entries")

    hashes: dict[str, str] = {}
    for raw_path, expected, label in records:
        key, path = _repository_path(raw_path, label)
        _require(
            isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
            "invalid SHA-256 for " + key,
        )
        digest = expected.lower()
        _require(key not in hashes, "duplicate manifest path: " + key)
        actual = _sha256_file(path)
        _require(actual.lower() == digest, "raw SHA changed: " + key)
        fixed = EXPECTED_FIXED_SHA256.get(key)
        if fixed is not None:
            _require(digest == fixed, "fixed SHA pin changed: " + key)
        hashes[key] = actual

    missing = sorted(set(REQUIRED_INPUTS) - set(hashes))
    _require(
        not missing,
        "manifest omits required frozen input(s): " + ", ".join(missing),
    )
    return hashes


def _checkpoint_sha(brief: Mapping[str, Any]) -> str:
    checkpoint = brief.get(f"{COMPONENT}_source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        checkpoint = brief.get("source_checkpoint")
    if isinstance(checkpoint, Mapping):
        raw_path = checkpoint.get("path")
        expected = checkpoint.get("sha256")
    else:
        raw_path = brief.get("source_checkpoint_path")
        expected = brief.get("source_checkpoint_sha256")
    _require(isinstance(raw_path, str) and raw_path, "pending brief source checkpoint path missing")
    candidate = Path(raw_path)
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    _require(candidate.resolve() == SOURCE.resolve(), "pending brief checkpoint points elsewhere")
    actual = _sha256_file(SOURCE)
    _require(
        isinstance(expected, str)
        and SHA256_RE.fullmatch(expected) is not None
        and actual.lower() == expected.lower(),
        "source manifest checkpoint SHA changed",
    )
    return actual


def _guard_true(source: Mapping[str, Any], names: tuple[str, ...], label: str) -> None:
    present = [name for name in names if name in source]
    _require(present, "manifest guard missing: " + label)
    _require(any(source[name] is True for name in present), "manifest guard is not true: " + label)


def _guard_false(source: Mapping[str, Any], names: tuple[str, ...], label: str) -> None:
    present = [name for name in names if name in source]
    _require(present, "manifest guard missing: " + label)
    _require(all(source[name] is False for name in present), "manifest guard is not false: " + label)


def _marker_count(marker: Any) -> int | None:
    if isinstance(marker, bool):
        return None
    if isinstance(marker, int):
        return marker
    if isinstance(marker, (list, tuple, set, frozenset, Mapping)):
        return len(marker)
    return None


def _scope_mappings(value: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    """Return the result plus known scope containers, without graph walking."""

    scopes: list[Mapping[str, Any]] = [value]
    for name in (
        "scope",
        "scope_flags",
        "graph_row_scope",
        "archive_scope",
        "bounded_scope",
        "source_scope",
        "integral_role_scope",
    ):
        candidate = value.get(name)
        if isinstance(candidate, Mapping):
            scopes.append(candidate)
    return scopes


def _scope_values(value: Mapping[str, Any], names: tuple[str, ...]) -> list[tuple[str, Any]]:
    found: list[tuple[str, Any]] = []
    for scope in _scope_mappings(value):
        for name in names:
            if name in scope:
                found.append((name, scope[name]))
    return found


def _walk_true_keys(value: Any, path: str = "$") -> list[str]:
    found: list[str] = []
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_text = str(key)
            child_path = path + "." + key_text
            if key_text in FORBIDDEN_TRUE_KEYS and item is True:
                found.append(child_path)
            found.extend(_walk_true_keys(item, child_path))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            found.extend(_walk_true_keys(item, path + "[" + str(index) + "]"))
    return found


def _require_false_group(
    value: Mapping[str, Any],
    names: tuple[str, ...],
    label: str,
    *,
    required: bool = True,
) -> None:
    present = _scope_values(value, names)
    if required:
        _require(bool(present), label + " flag is missing")
    for name, marker in present:
        _require(marker is False, label + " flag is true: " + name)


def _require_scope_false(value: Mapping[str, Any], label: str) -> None:
    forbidden = _walk_true_keys(value)
    _require(not forbidden, label + " set forbidden true scope: " + ", ".join(forbidden))
    _require_false_group(
        value,
        ("source_values_computed", "actual_source_values", "source_value_selected", "source_values_selected"),
        label + " source-value",
    )
    _require_false_group(
        value,
        ("integrals_evaluated", "integral_evaluated", "quadrature_evaluated", "integral_calls"),
        label + " integral",
    )
    _require_false_group(
        value,
        ("numeric_use_admissible", "numeric_values_computed"),
        label + " numeric-use",
    )
    _require_false_group(
        value,
        ("graph_traversal", "graph_edge_search", "source_graph_traversed", "old_graph_traversed", "bfs_used"),
        label + " graph traversal",
    )
    _require_false_group(
        value,
        ("protected_nodes_indexed", "protected_indexed", "protected_dereferenced"),
        label + " protected-node indexing",
    )
    _require_false_group(
        value,
        ("target_nodes_indexed", "target_node_indexed", "target_row_indexed", "target_node_dereferenced"),
        label + " target-node indexing",
    )
    _require_false_group(
        value,
        ("ancestor_replayed", "old_ancestor_replayed", "old_leaf_imported", "old_producer_called"),
        label + " old-code",
        required=False,
    )
    _require_false_group(
        value,
        ("duplicate_graph_computation", "repeated_graph_computation"),
        label + " duplicate computation",
        required=False,
    )
    _require_false_group(
        value,
        (
            "source_identity_proved",
            "same_source_function_proved",
            "function_identity_proved",
            "actual_source_identity_proved",
            "capture_free",
            "capture_free_proved",
            "ordinary_Z_proved",
            "ordinary_Z_transport_proved",
            "endpoint_identity_proved",
            "endpoint_transport_proved",
            "closure_proved",
            "m15_closure_proved",
            "certified",
        ),
        label + " semantic certification",
        required=False,
    )


def _require_fixture(value: Any) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), "synthetic fixture did not return an object")
    checks = value.get("checks")
    passed = (
        value.get("passed") is True
        or value.get("checks_passed") is True
        or value.get("fixture_passed") is True
        or (
            isinstance(checks, Mapping)
            and bool(checks)
            and all(item is True for item in checks.values())
        )
    )
    _require(passed, "synthetic fixture did not pass")
    _require_scope_false(value, "synthetic fixture")
    _json_bytes(_json_safe(value))
    return value


def _require_true_marker(
    value: Mapping[str, Any], names: tuple[str, ...], label: str
) -> str:
    candidates = _scope_values(value, names)
    _require(candidates, label + " marker is missing")
    for name, marker in candidates:
        if marker is True:
            return name
    raise RunnerInputError(label + " marker is not true")


def _row_id_list(marker: Any) -> list[int] | None:
    if isinstance(marker, (tuple, list, set, frozenset)):
        result: list[int] = []
        for item in marker:
            if isinstance(item, bool) or not isinstance(item, int):
                return None
            result.append(item)
        return result
    return None


def _direct_row_evidence(value: Mapping[str, Any]) -> dict[str, Any]:
    row_names = (
        "declared_graph_row_ids",
        "direct_graph_row_ids",
        "indexed_graph_row_ids",
        "graph_row_ids",
        "direct_row_ids",
    )
    candidates: list[tuple[str, list[int]]] = []
    for scope in _scope_mappings(value):
        for name in row_names:
            if name not in scope:
                continue
            ids = _row_id_list(scope[name])
            _require(ids is not None, name + " must be a list of integer row IDs")
            candidates.append((name, ids))
    _require(candidates, "integral-role result has no direct graph-row evidence")
    first = candidates[0][1]
    _require(len(first) == len(set(first)), "direct graph-row IDs contain duplicates")
    _require(
        DECLARED_GRAPH_ROW_IDS.issubset(first),
        "direct graph-row evidence omits a declared F row",
    )
    _require(
        len(first) <= MAX_DIRECT_GRAPH_ROWS,
        "direct graph-row evidence exceeds five rows plus one sibling",
    )
    extras = sorted(set(first) - DECLARED_GRAPH_ROW_IDS)
    _require(len(extras) <= 1, "more than one optional sibling row was indexed")
    _require(
        not (set(first) & (PROTECTED_NODE_IDS | TARGET_NODE_IDS)),
        "protected or target row was individually indexed",
    )
    for name, ids in candidates[1:]:
        _require(ids == first, "inconsistent direct graph-row evidence: " + name)
    return {
        "declared_rows": sorted(DECLARED_GRAPH_ROW_IDS),
        "direct_rows": first,
        "optional_sibling_rows": extras,
        "direct_row_count": len(first),
        "evidence_fields": [name for name, _ in candidates],
    }


def _integral_role_mapping_evidence(value: Mapping[str, Any]) -> dict[str, Any]:
    names = (
        "integral_role_mappings",
        "integral_roles",
        "archive_integral_roles",
        "source_integral_roles",
        "integral_role_mapping_count",
        "archive_integral_role_mapping_count",
        "integral_roles_read",
    )
    candidates = _scope_values(value, names)
    _require(candidates, "integral-role result has no archive mapping evidence")
    counts: list[tuple[str, int]] = []
    for name, marker in candidates:
        count = _marker_count(marker)
        if count is not None:
            counts.append((name, count))
    _require(counts, "archive integral-role mapping evidence has no count")
    _require(
        all(0 <= count <= MAX_INTEGRAL_ROLE_MAPPINGS for _, count in counts),
        "archive integral-role mapping count exceeds 36",
    )
    first = counts[0][1]
    _require(all(count == first for _, count in counts), "inconsistent integral-role mapping counts")
    return {
        "mapping_count": first,
        "max_mapping_count": MAX_INTEGRAL_ROLE_MAPPINGS,
        "evidence_fields": [name for name, _ in counts],
    }


def _validate_integral_roles(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "integral-role leaf did not return an object")
    _require_scope_false(value, "integral-role leaf")
    archive_marker = _require_true_marker(
        value,
        (
            "archive_whole_json_decoded",
            "whole_json_archive_decoded",
            "archive_json_decoded",
            "archive_decode_whole_json",
        ),
        "whole archive JSON decode",
    )
    direct_marker = _require_true_marker(
        value,
        (
            "bounded_graph_direct_reads",
            "bounded_direct_graph_reads",
            "direct_graph_reads_bounded",
        ),
        "bounded direct graph reads",
    )
    unresolved_marker = _require_true_marker(
        value,
        (
            "unresolved_matches_preserved",
            "duplicate_and_unresolved_matches_preserved",
            "nonhit_and_duplicate_matches_preserved",
            "unresolved_statuses_preserved",
        ),
        "unresolved/duplicate match preservation",
    )
    rows = _direct_row_evidence(value)
    mappings = _integral_role_mapping_evidence(value)
    return {
        "archive_whole_json_decode_marker": archive_marker,
        "bounded_direct_graph_reads_marker": direct_marker,
        "unresolved_matches_preserved_marker": unresolved_marker,
        "direct_graph_rows": rows,
        "integral_role_mappings": mappings,
        "structural_only": True,
        "same_source_function_proved": False,
        "capture_free_proved": False,
        "ordinary_Z_transport_proved": False,
        "endpoint_identity_proved": False,
        "closure_proved": False,
        "numeric_use_admissible": False,
    }


def _import_leaf() -> Any:
    """Import only the new leaf, and only after the exclusive entry receipt."""

    path = (ROOT / LEAF).resolve()
    _require(path.is_file(), "new integral-role leaf is missing")
    spec = importlib.util.spec_from_file_location(
        "f_pressure_path_integral_roles_v1_new_leaf", path
    )
    _require(spec is not None and spec.loader is not None, "new leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _callable(module: Any, name: str) -> Callable[..., Any]:
    function = getattr(module, name, None)
    _require(callable(function), "new leaf callable is missing: " + name)
    return function


def _invoke(
    steps: list[dict[str, Any]],
    name: str,
    function: Callable[..., Any],
    *args: Any,
) -> Any:
    record: dict[str, Any] = {"step": name, "invoked": True}
    try:
        value = function(*args)
    except BaseException as exc:
        record["exception"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
        steps.append(record)
        raise
    record["returned_packet"] = _json_safe(value)
    steps.append(record)
    return value


def _preflight() -> tuple[dict[str, Any], dict[str, str], str]:
    """Authenticate pending metadata and all raw bytes before entry."""

    _require(not ENTRY.exists() and not RESULT.exists(), "one-entry component already entered")
    _require(BRIEF.is_file(), "pending brief is missing")
    _require(SOURCE.is_file(), "source manifest is missing")
    brief = _load_json_object(BRIEF, "pending brief")
    source = _load_json_object(SOURCE, "source manifest")
    _require(brief.get("component") == COMPONENT, "brief component changed")
    _require(brief.get("phase") == PHASE, "brief phase changed")
    _require(
        brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("do_not_reenter") is True,
        "pending brief goal flags changed",
    )
    source_sha = _checkpoint_sha(brief)
    _require(source.get("component") == COMPONENT, "manifest component changed")
    _require(source.get("one_entry_only") is True, "manifest one-entry guard missing")
    _guard_true(
        source,
        (
            "all_raw_inputs_authenticated_before_decompression",
            "all_raw_inputs_authenticated_before_decode",
            "raw_sha256_authenticated_before_entry",
        ),
        "raw authentication",
    )
    _guard_true(
        source,
        (
            "archive_whole_json_decode",
            "archive_whole_json_decoded",
            "decode_entire_json_archive",
            "whole_json_archive_decode",
        ),
        "whole archive JSON decode",
    )
    _guard_true(
        source,
        (
            "bounded_graph_direct_reads",
            "bounded_direct_graph_reads",
            "direct_graph_reads_bounded",
        ),
        "bounded direct graph reads",
    )
    _guard_true(
        source,
        ("protected_and_target_nodes_not_individually_indexed", "no_protected_or_target_indexing"),
        "protected/target indexing",
    )
    _guard_true(
        source,
        (
            "unresolved_matches_preserved",
            "duplicate_and_unresolved_matches_preserved",
            "nonhit_and_duplicate_matches_preserved",
        ),
        "unresolved/duplicate preservation",
    )
    _guard_false(
        source,
        ("graph_traversal", "graph_edge_search", "source_graph_traversed"),
        "graph traversal",
    )
    _guard_false(
        source,
        ("actual_source_values", "source_values_computed", "numeric_use_admissible"),
        "source values/numeric use",
    )
    _guard_true(
        source,
        ("no_ancestor_replay", "no_old_ancestor_replay", "no_ancestor_math_replay"),
        "ancestor replay",
    )
    _guard_true(
        source,
        ("no_integrals", "no_integrals_evaluated", "no_quadrature"),
        "integral evaluation",
    )
    _guard_true(
        source,
        ("no_numeric_evaluation", "no_source_values", "no_numeric_values"),
        "numeric evaluation",
    )
    max_mapping = None
    for name in (
        "max_integral_role_mappings",
        "max_integral_role_mapping_count",
        "bounded_integral_role_mappings",
    ):
        marker = source.get(name)
        if isinstance(marker, int) and not isinstance(marker, bool):
            max_mapping = marker
            break
        if isinstance(marker, Mapping):
            candidate = marker.get("max")
            if isinstance(candidate, int) and not isinstance(candidate, bool):
                max_mapping = candidate
                break
    _require(
        isinstance(max_mapping, int) and 0 <= max_mapping <= MAX_INTEGRAL_ROLE_MAPPINGS,
        "manifest integral-role mapping bound must be in 0..36",
    )
    hashes = _manifest_hashes(source)
    _require(LEAF.as_posix() in hashes, "new integral-role leaf is not pinned")
    return dict(source), hashes, source_sha


def _recordable(value: Any) -> Any:
    """Retain a returned packet in a receipt without evaluating it."""

    return _json_safe(value) if value is not None else None


def main() -> int:
    # No failure packet is written before entry: untrusted pending metadata
    # must not consume the one-entry slot.
    source, hashes, source_sha = _preflight()
    entry_sha = _write_once(
        ENTRY,
        {
            "schema": COMPONENT + "_entry_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "input_hashes": hashes,
            "raw_sha256_authenticated_before_entry": True,
            "new_leaf_imported": False,
            "consumer_calls_started": 0,
            "verification_invocations": 0,
            "archive_whole_json_decoded": False,
            "bounded_graph_direct_reads": False,
            "graph_traversal": False,
            "actual_source_values": False,
            "protected_nodes_indexed": False,
            "target_nodes_indexed": False,
            "integral_role_mapping_count": 0,
            "integrals_evaluated": False,
            "numeric_use_admissible": False,
        },
    )

    steps: list[dict[str, Any]] = []
    fixture_result: Any = None
    roles_result: Any = None
    role_metrics: dict[str, Any] | None = None
    try:
        module = _import_leaf()
        fixture = _callable(module, "verify_synthetic_fixture")
        fixture_result = _invoke(steps, "verify_synthetic_fixture", fixture)
        _require_fixture(fixture_result)

        builder = _callable(module, "build_integral_roles")
        roles_result = _invoke(
            steps,
            "build_integral_roles",
            builder,
            ROOT,
            source,
            hashes,
        )
        role_metrics = _validate_integral_roles(roles_result)

        output: dict[str, Any] = {
            "schema": COMPONENT + "_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "source_checkpoint_sha256": source_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "consumer_calls": [item["step"] for item in steps],
            "verification_invocations": len(steps),
            "archive_whole_json_decoded": True,
            "bounded_graph_direct_reads": True,
            "graph_traversal": False,
            "actual_source_values": False,
            "protected_nodes_indexed": False,
            "target_nodes_indexed": False,
            "integrals_evaluated": False,
            "numeric_use_admissible": False,
            "role_metrics": _json_safe(role_metrics),
            "result": {
                "synthetic_fixture": _json_safe(fixture_result),
                "integral_roles": _json_safe(roles_result),
            },
            "source_manifest": _json_safe(source),
            "steps": _json_safe(steps),
            "error": None,
        }
        _json_bytes(output)
    except BaseException as exc:
        output = {
            "schema": COMPONENT + "_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "source_checkpoint_sha256": source_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "consumer_calls": [item["step"] for item in steps],
            "verification_invocations": len(steps),
            "archive_whole_json_decoded": bool(role_metrics),
            "bounded_graph_direct_reads": bool(role_metrics),
            "graph_traversal": False,
            "actual_source_values": False,
            "protected_nodes_indexed": False,
            "target_nodes_indexed": False,
            "integrals_evaluated": False,
            "numeric_use_admissible": False,
            "role_metrics": _recordable(role_metrics),
            "result": {
                "synthetic_fixture": _recordable(fixture_result),
                "integral_roles": _recordable(roles_result),
            },
            "source_manifest": _json_safe(source),
            "steps": _json_safe(steps),
            "error": {
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
        }

    result_sha = _write_once(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "result_sha256": result_sha,
                "consumer_calls": output["consumer_calls"],
                "verification_invocations": output["verification_invocations"],
                "archive_whole_json_decoded": output["archive_whole_json_decoded"],
                "bounded_graph_direct_reads": output["bounded_graph_direct_reads"],
                "graph_traversal": output["graph_traversal"],
                "actual_source_values": output["actual_source_values"],
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
