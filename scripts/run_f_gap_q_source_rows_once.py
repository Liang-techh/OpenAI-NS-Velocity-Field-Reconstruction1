"""Run the bounded F gap Q-source-row leaf exactly once.

This runner is the persistence boundary for ``f_gap_q_source_rows_v1``.  It
authenticates the four raw inputs named by the pending manifest before the
exclusive-entry receipt is written.  Only the new
``f_gap_q_source_rows.py`` leaf is imported after that receipt.  The leaf may
decode the authenticated F gzip as one JSON object, but this runner accepts
only the two direct source rows ``Q0=277042`` and ``Q1=277055``.  Any child
rows are opaque references: they are not direct reads, indexed rows, or graph
edges.

The synthetic fixture and the real source consumer are separate calls.  The
fixture is recorded only as synthetic evidence and is never used to support
the source packet.  No old producer, ancestor, graph walk, numerical
evaluation, network operation, function-identity claim, derivative proof,
error deletion, or closure claim is admitted here.  A preflight failure does
not consume the one-entry slot; after entry, every exception is persisted in a
complete result receipt before its hash is printed.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
import traceback
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "f_gap_q_source_rows_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

RUNNER = Path("scripts/run_f_gap_q_source_rows_once.py")
LEAF = Path("src/openai_ns_reconstruction/f_gap_q_source_rows.py")
F_ARCHIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
PRIOR_QSQUARE_RESULT = Path("tmp/f_gap_qsquare_transport_v1/result.json")

REQUIRED_INPUTS = (
    LEAF.as_posix(),
    RUNNER.as_posix(),
    F_ARCHIVE.as_posix(),
    PRIOR_QSQUARE_RESULT.as_posix(),
)

F_ARCHIVE_SHA256 = (
    "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"
)
PRIOR_QSQUARE_RESULT_SHA256 = (
    "b516e5a72290062d788114f25237ab5019abc217838a68d125de24caa1323d95"
)
EXPECTED_FIXED_SHA256 = {
    F_ARCHIVE.as_posix(): F_ARCHIVE_SHA256,
    PRIOR_QSQUARE_RESULT.as_posix(): PRIOR_QSQUARE_RESULT_SHA256,
}

SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")

Q_SOURCE_ROWS = frozenset({277042, 277055})
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
        33982,
    }
)

# A returned packet may carry these names only with a false value.  The walk
# is recursive so a nested source-gate object cannot promote an observation to
# a proof or an evaluated source result.
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
        "child_indexed",
        "children_indexed",
        "child_row_indexed",
        "child_rows_indexed",
        "child_dereferenced",
        "children_dereferenced",
        "child_rows_dereferenced",
        "child_graph_rows_indexed",
        "source_values_computed",
        "actual_source_values",
        "source_value_selected",
        "source_values_selected",
        "numeric_values_computed",
        "numeric_use_admissible",
        "integrals_evaluated",
        "integral_evaluated",
        "quadrature_evaluated",
        "graph_traversal",
        "graph_edge_search",
        "source_graph_traversed",
        "old_graph_traversed",
        "bfs_used",
        "ancestor_replayed",
        "old_ancestor_replayed",
        "old_leaf_imported",
        "old_producer_imported",
        "old_producer_called",
        "ancestor_imported",
        "ancestor_called",
        "network_used",
        "network_called",
        "source_identity",
        "source_identity_proved",
        "same_source",
        "same_source_proved",
        "same_source_function_proved",
        "function_identity",
        "function_identity_proved",
        "real_function_identity",
        "real_function_identity_proved",
        "actual_source_identity_proved",
        "Q1_is_partial_Z_Q0_proved",
        "derivative_proved",
        "derivatives_proved",
        "derivative_identity_proved",
        "ordinary_Z_proved",
        "ordinary_z_proved",
        "ordinary_Z_derivative_proved",
        "ordinary_Z_transport_proved",
        "ordinary_z_transport_proved",
        "transport_proved",
        "capture_free",
        "capture_free_proved",
        "capture_free_substitution_proved",
        "endpoint_transport_proved",
        "closure",
        "closed",
        "closure_proved",
        "q_square_closure_proved",
        "m15_closure_proved",
        "error_removal",
        "error_transport_proved",
        "error_removed",
        "old_error_removed",
        "errors_removed",
        "formal_missing_correction_deleted",
        "forbidden_nodes_indexed",
        "semantic_identity_proved",
        "quantitative_error_reduced",
    }
)

SEMANTIC_FALSE_KEYS = frozenset(
    {
        "source_identity",
        "source_identity_proved",
        "same_source",
        "same_source_proved",
        "same_source_function_proved",
        "function_identity",
        "function_identity_proved",
        "real_function_identity",
        "real_function_identity_proved",
        "actual_source_identity_proved",
        "Q1_is_partial_Z_Q0_proved",
        "derivative_proved",
        "derivatives_proved",
        "derivative_identity_proved",
        "ordinary_Z_proved",
        "ordinary_z_proved",
        "ordinary_Z_derivative_proved",
        "ordinary_Z_transport_proved",
        "ordinary_z_transport_proved",
        "transport_proved",
        "capture_free",
        "capture_free_proved",
        "capture_free_substitution_proved",
        "endpoint_transport_proved",
        "closure",
        "closed",
        "closure_proved",
        "q_square_closure_proved",
        "m15_closure_proved",
        "source_values_computed",
        "actual_source_values",
        "source_value_selected",
        "source_values_selected",
        "numeric_values_computed",
        "numeric_use_admissible",
        "integrals_evaluated",
        "integral_evaluated",
        "quadrature_evaluated",
        "graph_traversal",
        "graph_edge_search",
        "source_graph_traversed",
        "old_graph_traversed",
        "bfs_used",
        "ancestor_replayed",
        "old_ancestor_replayed",
        "old_leaf_imported",
        "old_producer_imported",
        "old_producer_called",
        "ancestor_imported",
        "ancestor_called",
        "network_used",
        "network_called",
        "error_removal",
        "error_transport_proved",
        "error_removed",
        "old_error_removed",
        "errors_removed",
        "formal_missing_correction_deleted",
        "forbidden_nodes_indexed",
        "semantic_identity_proved",
        "quantitative_error_reduced",
    }
)

DIRECT_ROW_KEYS = frozenset(
    {
        "direct_graph_row_ids",
        "direct_source_row_ids",
        "direct_row_ids",
        "direct_read_row_ids",
        "direct_read_rows",
        "indexed_graph_row_ids",
        "indexed_source_row_ids",
        "indexed_row_ids",
        "indexed_rows",
        "declared_graph_row_ids",
        "source_row_ids",
        "requested_row_ids",
        "requested_rows",
        "requested_Q_rows",
    }
)


class RunnerInputError(ValueError):
    """Raised when the frozen one-shot contract is not satisfied."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha256_file(path: Path) -> str:
    """Hash raw bytes only; no authenticated input is decoded by the runner."""

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


def _manifest_entries(source: Mapping[str, Any]) -> list[tuple[Any, Any, str]]:
    entries: list[tuple[Any, Any, str]] = []
    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(
                isinstance(item, Mapping),
                f"source.files[{index}] is not an object",
            )
            entries.append(
                (item.get("path"), item.get("sha256"), f"source.files[{index}]")
            )
    else:
        for field in ("sha256", "input_hashes", "input_sha256", "raw_input_sha256"):
            table = source.get(field)
            if isinstance(table, Mapping):
                entries.extend(
                    (raw_path, expected, f"source.{field}")
                    for raw_path, expected in table.items()
                )
    _require(entries, "source manifest contains no direct raw SHA entries")
    return entries


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    """Authenticate exactly four raw paths before importing the new leaf."""

    hashes: dict[str, str] = {}
    for raw_path, expected, label in _manifest_entries(source):
        key, path = _repository_path(raw_path, label)
        _require(key not in hashes, "duplicate manifest path: " + key)
        _require(
            isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
            "invalid SHA-256 for " + key,
        )
        digest = expected.lower()
        actual = _sha256_file(path)
        _require(actual.lower() == digest, "raw SHA changed: " + key)
        fixed = EXPECTED_FIXED_SHA256.get(key)
        if fixed is not None:
            _require(digest == fixed, "fixed SHA pin changed: " + key)
        hashes[key] = actual

    required = set(REQUIRED_INPUTS)
    _require(
        set(hashes) == required,
        "manifest must freeze exactly four inputs; missing or extra paths: "
        + ", ".join(sorted(required ^ set(hashes))),
    )
    return hashes


def _checkpoint_sha(brief: Mapping[str, Any]) -> str:
    checkpoint = brief.get(f"{COMPONENT}_source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        checkpoint = brief.get("source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "pending brief source checkpoint missing")
    raw_path = checkpoint.get("path")
    _require(isinstance(raw_path, str) and raw_path, "source checkpoint path missing")
    candidate = Path(raw_path)
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    _require(
        candidate.resolve() == SOURCE.resolve(),
        "pending brief checkpoint points elsewhere",
    )
    expected = checkpoint.get("sha256")
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
    _require(
        any(source[name] is True for name in present),
        "manifest guard is not true: " + label,
    )


def _guard_false(source: Mapping[str, Any], names: tuple[str, ...], label: str) -> None:
    present = [name for name in names if name in source]
    _require(present, "manifest guard missing: " + label)
    _require(
        all(source[name] is False for name in present),
        "manifest guard is not false: " + label,
    )


def _walk_true_keys(value: Any, path: str = "$") -> list[str]:
    found: list[str] = []
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_text = str(key)
            child = path + "." + key_text
            if key_text in FORBIDDEN_TRUE_KEYS and item is True:
                found.append(child)
            found.extend(_walk_true_keys(item, child))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            found.extend(_walk_true_keys(item, path + "[" + str(index) + "]"))
    return found


def _walk_key_values(
    value: Any, names: set[str], path: str = "$"
) -> list[tuple[str, Any]]:
    found: list[tuple[str, Any]] = []
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_text = str(key)
            child = path + "." + key_text
            if key_text in names:
                found.append((child, item))
            found.extend(_walk_key_values(item, names, child))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            found.extend(_walk_key_values(item, names, path + "[" + str(index) + "]"))
    return found


def _validate_false_scope(value: Mapping[str, Any], label: str) -> None:
    forbidden = _walk_true_keys(value)
    _require(
        not forbidden,
        label + " claims forbidden true scope: " + ", ".join(forbidden),
    )
    for path, marker in _walk_key_values(value, set(SEMANTIC_FALSE_KEYS)):
        _require(marker is False, label + " semantic flag is true: " + path)


def _require_fixture(value: Any) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), "synthetic fixture did not return an object")
    passed = value.get("passed") is True or value.get("checks_passed") is True
    if not passed:
        checks = value.get("checks")
        passed = (
            isinstance(checks, Mapping)
            and bool(checks)
            and all(item is True for item in checks.values())
        )
    _require(passed, "synthetic fixture did not pass")
    synthetic_only = value.get("synthetic_only") is True
    scope = value.get("scope")
    synthetic_only = synthetic_only or (
        isinstance(scope, Mapping) and scope.get("synthetic_only") is True
    )
    _require(synthetic_only, "synthetic fixture is not marked synthetic-only")
    _validate_false_scope(value, "synthetic fixture")
    _json_bytes(_json_safe(value))
    return value


def _row_ids(value: Any) -> list[int] | None:
    if isinstance(value, (list, tuple, set, frozenset)):
        if all(type(item) is int for item in value):
            return list(value)
        ids: list[int] = []
        for item in value:
            if not isinstance(item, Mapping):
                return None
            candidate = next(
                (
                    item[name]
                    for name in ("node", "node_id", "row", "row_id", "id")
                    if name in item
                ),
                None,
            )
            if type(candidate) is not int:
                return None
            ids.append(candidate)
        return ids
    return None


def _direct_row_scope(value: Mapping[str, Any]) -> dict[str, Any]:
    candidates: list[tuple[str, list[int]]] = []
    for path, marker in _walk_key_values(value, set(DIRECT_ROW_KEYS)):
        ids = _row_ids(marker)
        if ids is not None:
            candidates.append((path, ids))
    _require(candidates, "Q-source result has no direct-row evidence")
    first = candidates[0][1]
    _require(len(first) == len(set(first)), "direct-row evidence contains duplicates")
    _require(set(first) == Q_SOURCE_ROWS, "direct-row evidence is not exactly Q0/Q1")
    _require(
        not (set(first) & PROTECTED_NODE_IDS),
        "protected or target row was indexed",
    )
    for path, ids in candidates[1:]:
        _require(ids == first, "inconsistent direct-row evidence: " + path)
    return {
        "required_rows": sorted(Q_SOURCE_ROWS),
        "direct_rows": first,
        "direct_row_count": len(first),
        "children_indexed": False,
        "evidence_paths": [path for path, _ in candidates],
    }


def _import_leaf() -> Any:
    path = (ROOT / LEAF).resolve()
    _require(path.is_file(), "new Q-source leaf is missing")
    spec = importlib.util.spec_from_file_location(
        "f_gap_q_source_rows_v1_new_leaf", path
    )
    _require(spec is not None and spec.loader is not None, "new Q-source leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _callable(module: Any, name: str) -> Callable[..., Any]:
    function = getattr(module, name, None)
    _require(callable(function), "new leaf callable is missing: " + name)
    return function


def _validate_source_result(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "Q-source leaf did not return an object")
    _validate_false_scope(value, "Q-source leaf")
    _require(
        value.get("source_only") is True
        or value.get("source_rows_only") is True,
        "Q-source leaf is not source-only",
    )
    row_scope = _direct_row_scope(value)
    # These markers are checked when present.  Missing markers do not turn a
    # structural packet into a proof, but a positive marker is rejected.
    for key in (
        "children_indexed",
        "child_rows_indexed",
        "children_dereferenced",
        "child_rows_dereferenced",
        "graph_edges_followed",
        "graph_traversal",
        "source_graph_traversed",
    ):
        if key in value:
            _require(value[key] is False, "Q-source leaf set forbidden flag true: " + key)
    return {
        "structural_only": True,
        "direct_row_scope": row_scope,
        "source_identity_proved": False,
        "function_identity_proved": False,
        "derivative_proved": False,
        "capture_free_substitution_proved": False,
        "closure_proved": False,
        "actual_source_values": False,
        "numeric_use_admissible": False,
        "error_removal": False,
    }


def _preflight() -> tuple[dict[str, Any], dict[str, str], str]:
    # Nothing before this return path writes the entry or result receipt.
    _require(
        not ENTRY.exists() and not RESULT.exists(),
        "one-entry component already entered",
    )
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
        ("archive_whole_json_decode", "archive_whole_json_decoded", "whole_json_archive_decode"),
        "whole archive JSON decode",
    )
    _guard_true(
        source,
        ("bounded_graph_direct_reads", "bounded_direct_graph_reads", "direct_graph_reads_bounded"),
        "bounded direct reads",
    )
    _guard_true(
        source,
        (
            "protected_and_target_nodes_not_individually_indexed",
            "no_protected_or_target_indexing",
        ),
        "protected/target indexing",
    )
    _guard_true(
        source,
        (
            "children_collected_as_references_only",
            "children_references_only",
            "child_refs_only",
            "children_not_indexed",
            "child_rows_are_references_only",
        ),
        "child references",
    )
    _guard_false(
        source,
        ("graph_traversal", "graph_edge_search", "source_graph_traversed", "bfs_used"),
        "graph traversal",
    )
    _guard_true(
        source,
        (
            "no_old_code_import",
            "no_old_leaf_import",
            "no_old_producer_import",
            "no_ancestor_replay",
            "no_old_ancestor_replay",
        ),
        "legacy/ancestor imports",
    )
    _guard_true(
        source,
        ("no_numeric_evaluation", "no_source_values", "no_numeric_values"),
        "numeric evaluation",
    )
    _guard_true(
        source,
        ("no_integrals", "no_integrals_evaluated", "no_quadrature"),
        "integral evaluation",
    )
    # A manifest may express these prohibitions either as explicit false
    # result flags or by omitting the fields.  Whichever representation it
    # uses, a positive claim is never admitted.
    _validate_false_scope(source, "source manifest")
    hashes = _manifest_hashes(source)
    _require(LEAF.as_posix() in hashes, "new Q-source leaf is not pinned")
    _require(RUNNER.as_posix() in hashes, "new Q-source runner is not pinned")
    return dict(source), hashes, source_sha


def _recordable(value: Any) -> Any:
    return _json_safe(value) if value is not None else None


def main() -> int:
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
            "new_leaf_imported": False,
            "consumer_calls_started": 0,
            "synthetic_fixture_invocations": 0,
            "source_consumer_invocations": 0,
            "archive_whole_json_decoded": False,
            "direct_graph_rows_indexed": [],
            "child_rows_indexed": False,
            "graph_traversal": False,
            "actual_source_values": False,
            "protected_nodes_indexed": False,
            "target_nodes_indexed": False,
            "integrals_evaluated": False,
            "numeric_use_admissible": False,
            "function_identity_proved": False,
            "derivative_proved": False,
            "error_removal": False,
            "closure_proved": False,
            "runtime_model_metadata": None,
        },
    )

    steps: list[dict[str, Any]] = []
    fixture_result: Any = None
    source_result: Any = None
    source_metrics: dict[str, Any] | None = None
    leaf_imported = False
    try:
        module = _import_leaf()
        leaf_imported = True

        fixture = _callable(module, "verify_synthetic_fixture")
        steps.append({"step": "verify_synthetic_fixture", "invoked": True})
        fixture_result = fixture()
        steps[-1]["returned_packet"] = _recordable(fixture_result)
        _require_fixture(fixture_result)

        builder = _callable(module, "build_q_source_rows")
        steps.append({"step": "build_q_source_rows", "invoked": True})
        # This is the sole real source consumer call.  It receives only the
        # repository root, frozen manifest, and already authenticated hashes.
        source_result = builder(ROOT, source, hashes)
        steps[-1]["returned_packet"] = _recordable(source_result)
        source_metrics = _validate_source_result(source_result)

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
            "synthetic_fixture_invocations": 1,
            "source_consumer_invocations": 1,
            "leaf_imported": leaf_imported,
            "archive_whole_json_decoded": True,
            "direct_graph_rows_indexed": source_metrics["direct_row_scope"]["direct_rows"],
            "child_rows_indexed": False,
            "graph_traversal": False,
            "actual_source_values": False,
            "source_values_computed": False,
            "protected_nodes_indexed": False,
            "target_nodes_indexed": False,
            "integrals_evaluated": False,
            "numeric_use_admissible": False,
            "source_identity_proved": False,
            "function_identity_proved": False,
            "derivative_proved": False,
            "ordinary_Z_derivative_proved": False,
            "capture_free_substitution_proved": False,
            "closure_proved": False,
            "error_removal": False,
            "source_metrics": _json_safe(source_metrics),
            "result": {
                "synthetic_fixture": _json_safe(fixture_result),
                "q_source_rows": _json_safe(source_result),
            },
            "source_manifest": _json_safe(source),
            "steps": _json_safe(steps),
            "error": None,
            "runtime_model_metadata": None,
        }
        _validate_false_scope(output, "runner success result")
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
            "synthetic_fixture_invocations": sum(
                item.get("step") == "verify_synthetic_fixture" for item in steps
            ),
            "source_consumer_invocations": sum(
                item.get("step") == "build_q_source_rows" for item in steps
            ),
            "leaf_imported": leaf_imported,
            "archive_whole_json_decoded": bool(source_metrics),
            "direct_graph_rows_indexed": [],
            "child_rows_indexed": False,
            "graph_traversal": False,
            "actual_source_values": False,
            "source_values_computed": False,
            "protected_nodes_indexed": False,
            "target_nodes_indexed": False,
            "integrals_evaluated": False,
            "numeric_use_admissible": False,
            "source_identity_proved": False,
            "function_identity_proved": False,
            "derivative_proved": False,
            "ordinary_Z_derivative_proved": False,
            "capture_free_substitution_proved": False,
            "closure_proved": False,
            "error_removal": False,
            "source_metrics": _recordable(source_metrics),
            "result": {
                "synthetic_fixture": _recordable(fixture_result),
                "q_source_rows": _recordable(source_result),
            },
            "source_manifest": _json_safe(source),
            "steps": _json_safe(steps),
            "error": {
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
            "runtime_model_metadata": None,
        }

    # Persist the complete JSON-safe success or failure packet before printing.
    result_sha = _write_once(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "result_sha256": result_sha,
                "consumer_calls": output["consumer_calls"],
                "synthetic_fixture_invocations": output["synthetic_fixture_invocations"],
                "source_consumer_invocations": output["source_consumer_invocations"],
                "direct_graph_rows_indexed": output["direct_graph_rows_indexed"],
                "child_rows_indexed": output["child_rows_indexed"],
                "graph_traversal": output["graph_traversal"],
                "actual_source_values": output["actual_source_values"],
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
