"""Run the protected-reference metadata correction exactly once.

This runner is the persistence boundary for ``f_gap_protected_ref_correction_v1``.
The frozen manifest contains only this runner, the new correction leaf, and the
already-returned quotient-operands result.  All three byte streams are hashed
before the exclusive entry receipt is written.  The prior result is not parsed
here; the leaf may parse it only after that authentication and must not open an
archive, import an old producer, traverse a graph, or evaluate source values.

The leaf interface is intentionally exact::

    verify_synthetic_fixture() -> mapping
    build_protected_ref_correction(root, source, hashes) -> mapping

Each callable is invoked at most once.  A returned packet is accepted only if
the prior edge targeting node ``1`` is marked protected while retaining its
original ``False`` marker and its non-indexed scope.  Every other returned edge
is checked against the frozen protected-node set.  Source identity, derivative,
positivity, error-removal, and closure claims remain explicitly out of scope.
Both success and failure packets are written before a result SHA is printed.
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
COMPONENT = "f_gap_protected_ref_correction_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

RUNNER = Path("scripts/run_f_gap_protected_ref_correction_once.py")
LEAF = Path("src/openai_ns_reconstruction/f_gap_protected_ref_correction.py")
PRIOR_RESULT = Path("tmp/f_gap_q_quotient_operands_v1/result.json")

# The manifest is deliberately exact: no archive and no other producer input
# may be smuggled into this correction component.
REQUIRED_INPUTS = (
    RUNNER.as_posix(),
    LEAF.as_posix(),
    PRIOR_RESULT.as_posix(),
)
PRIOR_RESULT_SHA256 = (
    "0ae2076688bd01072504f282bf82a2c3346555abdbb7ec4a8d2828a4ab2a80cc"
)
EXPECTED_FIXED_SHA256 = {PRIOR_RESULT.as_posix(): PRIOR_RESULT_SHA256}
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")

# This is the protected set carried by the preceding quotient-operands
# observation.  IDs are compared as edge metadata only; none is indexed.
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
TARGET_NODE = 1
EXPECTED_DIRECT_EDGE_COUNT = 6


# Returned packets are recursively checked for accidental promotion of a
# structural observation into a graph/numerical/provenance claim.  The same
# names are required to remain explicitly false when present.  The correction
# itself is allowed to set ``protected_target=True`` on an edge.
SEMANTIC_FALSE_KEYS = frozenset(
    {
        "archive_read",
        "archive_opened",
        "archive_decoded",
        "archive_decompressed",
        "archive_whole_json_decoded",
        "old_code_imported",
        "old_leaf_imported",
        "old_producer_imported",
        "old_producer_called",
        "ancestor_imported",
        "ancestor_called",
        "ancestor_replayed",
        "graph_traversal",
        "graph_edge_search",
        "graph_edges_followed",
        "edges_followed",
        "source_graph_traversed",
        "bfs_used",
        "children_indexed",
        "child_rows_indexed",
        "child_indexed",
        "child_dereferenced",
        "children_dereferenced",
        "protected_indexed",
        "protected_node_indexed",
        "protected_nodes_indexed",
        "protected_ref_indexed",
        "protected_ref_dereferenced",
        "node_indexed",
        "node_1_indexed",
        "direct_nodes_indexed",
        "target_node_indexed",
        "target_nodes_indexed",
        "target_node_dereferenced",
        "indexed",
        "source_values_computed",
        "actual_source_values",
        "source_value_selected",
        "numeric_values_computed",
        "numeric_evaluation",
        "numeric_use_admissible",
        "integrals_evaluated",
        "integral_evaluated",
        "quadrature_evaluated",
        "source_identity",
        "source_identity_proved",
        "same_source",
        "same_source_proved",
        "function_identity",
        "function_identity_proved",
        "ordinary_Z",
        "ordinary_z",
        "ordinary_Z_proved",
        "ordinary_z_proved",
        "ordinary_Z_transport",
        "ordinary_z_transport",
        "ordinary_Z_transport_proved",
        "ordinary_z_transport_proved",
        "derivative",
        "derivative_proved",
        "ordinary_Z_derivative",
        "ordinary_Z_derivative_proved",
        "positivity",
        "positive",
        "positivity_proved",
        "positivity_claimed",
        "strictly_positive",
        "nonnegative",
        "denominator_positive_proved",
        "error_removal",
        "error_removed",
        "errors_removed",
        "old_error_removed",
        "error_transport_proved",
        "error_bound_proved",
        "closure",
        "closed",
        "full_closure",
        "closure_proved",
    }
)
FORBIDDEN_TRUE_KEYS = SEMANTIC_FALSE_KEYS

ROW_KEYS = ("Q0_numerator", "Q1_numerator", "shared_denominator")


class RunnerInputError(ValueError):
    """Raised when the frozen one-shot contract is not satisfied."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha256_file(path: Path) -> str:
    """Hash raw bytes only; do not decode prior result or any archive here."""

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
                "source.files[" + str(index) + "] is not an object",
            )
            entries.append(
                (item.get("path"), item.get("sha256"), "source.files[" + str(index) + "]")
            )
    else:
        for field in ("sha256", "input_hashes", "input_sha256", "raw_input_sha256"):
            table = source.get(field)
            if isinstance(table, Mapping):
                entries.extend(
                    (raw_path, expected, "source." + field)
                    for raw_path, expected in table.items()
                )
    _require(entries, "source manifest contains no raw SHA entries")
    return entries


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    """Authenticate exactly runner, leaf, and prior result before entry."""

    hashes: dict[str, str] = {}
    for raw_path, expected, label in _manifest_entries(source):
        key, path = _repository_path(raw_path, label)
        _require(key not in hashes, "duplicate manifest path: " + key)
        _require(
            isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
            "invalid SHA-256 for " + key,
        )
        actual = _sha256_file(path)
        _require(actual.lower() == expected.lower(), "raw SHA mismatch: " + key)
        fixed = EXPECTED_FIXED_SHA256.get(key)
        if fixed is not None:
            _require(actual.lower() == fixed, "fixed SHA pin changed: " + key)
        hashes[key] = actual

    required = set(REQUIRED_INPUTS)
    _require(
        set(hashes) == required,
        "manifest must contain exactly runner, leaf, and prior result; "
        + "missing or extra paths: "
        + ", ".join(sorted(required ^ set(hashes))),
    )
    _require(PRIOR_RESULT.as_posix() in hashes, "prior result is not authenticated")
    _require(hashes[PRIOR_RESULT.as_posix()].lower() == PRIOR_RESULT_SHA256,
             "prior result fixed SHA changed")
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
    _require(candidate.resolve() == SOURCE.resolve(),
             "pending brief source checkpoint points elsewhere")
    expected = checkpoint.get("sha256")
    actual = _sha256_file(SOURCE)
    _require(
        isinstance(expected, str)
        and SHA256_RE.fullmatch(expected) is not None
        and actual.lower() == expected.lower(),
        "source checkpoint raw SHA mismatch",
    )
    return actual


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
    value: Any, names: frozenset[str], path: str = "$") -> list[tuple[str, Any]]:
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
    _require(not forbidden, label + " claims forbidden scope: " + ", ".join(forbidden))
    for path, marker in _walk_key_values(value, SEMANTIC_FALSE_KEYS):
        _require(marker is False, label + " semantic scope is true: " + path)


def _guard_true(source: Mapping[str, Any], names: tuple[str, ...], label: str) -> None:
    present = [name for name in names if name in source]
    _require(present, "manifest guard missing: " + label)
    _require(any(source[name] is True for name in present),
             "manifest guard is not true: " + label)


def _guard_false(source: Mapping[str, Any], names: tuple[str, ...], label: str) -> None:
    present = [name for name in names if name in source]
    _require(present, "manifest guard missing: " + label)
    _require(all(source[name] is False for name in present),
             "manifest guard is not false: " + label)


def _guard_not_claimed(
    source: Mapping[str, Any],
    true_names: tuple[str, ...],
    false_names: tuple[str, ...],
    label: str,
) -> None:
    positive = [name for name in true_names if name in source]
    if positive:
        _require(any(source[name] is True for name in positive),
                 "manifest guard is not true: " + label)
        return
    negative = [name for name in false_names if name in source]
    _require(negative, "manifest guard missing: " + label)
    _require(all(source[name] is False for name in negative),
             "manifest guard is not false: " + label)


def _preflight() -> tuple[dict[str, Any], dict[str, str], str]:
    # Nothing before this return path writes entry/result or imports the leaf.
    _require(not ENTRY.exists() and not RESULT.exists(),
             "one-entry component already entered")
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
            "all_raw_inputs_authenticated_before_decode",
            "all_raw_inputs_authenticated_before_decompression",
            "all_raw_inputs_authenticated_before_entry",
            "raw_sha256_authenticated_before_entry",
        ),
        "raw authentication",
    )
    _guard_not_claimed(
        source,
        ("no_archive_read", "archive_not_read", "no_archive"),
        ("archive_read", "archive_opened", "archive_decoded", "archive_whole_json_decoded"),
        "archive exclusion",
    )
    _guard_not_claimed(
        source,
        ("no_old_code_import", "no_old_leaf_import", "no_old_producer_import"),
        ("old_code_imported", "old_leaf_imported", "old_producer_imported"),
        "old-code exclusion",
    )
    _guard_not_claimed(
        source,
        ("no_graph_traversal", "no_source_graph_traversal"),
        ("graph_traversal", "source_graph_traversed", "graph_edge_search"),
        "graph traversal",
    )
    _guard_not_claimed(
        source,
        ("no_numeric_evaluation", "no_source_values", "no_numeric_values"),
        ("numeric_evaluation", "source_values_computed", "actual_source_values"),
        "numeric evaluation",
    )
    _guard_not_claimed(
        source,
        ("protected_nodes_not_indexed", "no_protected_node_indexing"),
        ("protected_nodes_indexed", "protected_node_indexed", "target_node_indexed"),
        "protected-node indexing",
    )
    _guard_not_claimed(
        source,
        ("no_source_identity_claim", "source_identity_not_claimed", "source_identity_unproved"),
        ("source_identity", "source_identity_proved", "same_source_proved"),
        "source identity",
    )
    _guard_not_claimed(
        source,
        ("no_derivative_claim", "derivative_not_claimed", "derivative_unproved"),
        ("derivative", "derivative_proved", "ordinary_Z_derivative_proved"),
        "derivative",
    )
    _guard_not_claimed(
        source,
        ("no_positivity_claim", "no_positive_claim", "positivity_unproved"),
        ("positivity", "positivity_proved", "positive"),
        "positivity",
    )
    _guard_not_claimed(
        source,
        ("no_error_removal", "errors_retained", "error_removal_unproved"),
        ("error_removal", "error_removed", "errors_removed"),
        "error removal",
    )
    _guard_not_claimed(
        source,
        ("no_closure_claim", "closure_unproved", "not_closed"),
        ("closure", "closure_proved", "closed"),
        "closure",
    )
    _validate_false_scope(source, "source manifest")
    hashes = _manifest_hashes(source)
    _require(LEAF.as_posix() in hashes, "new correction leaf is not pinned")
    _require(RUNNER.as_posix() in hashes, "new correction runner is not pinned")
    return dict(source), hashes, source_sha


def _import_leaf() -> Any:
    path = (ROOT / LEAF).resolve()
    _require(path.is_file(), "new correction leaf is missing")
    spec = importlib.util.spec_from_file_location(
        "f_gap_protected_ref_correction_new_leaf", path
    )
    _require(spec is not None and spec.loader is not None,
             "new correction leaf load failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _callable(module: Any, name: str) -> Callable[..., Any]:
    function = getattr(module, name, None)
    _require(callable(function), "new correction leaf callable is missing: " + name)
    return function


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
    if isinstance(scope, Mapping):
        synthetic_only = synthetic_only or scope.get("synthetic_only") is True
    _require(synthetic_only, "synthetic fixture is not marked synthetic-only")
    _validate_false_scope(value, "synthetic fixture")
    _json_bytes(_json_safe(value))
    return value


def _node_from_value(value: Any) -> int | None:
    if type(value) is int and value >= 0:
        return value
    if not isinstance(value, Mapping):
        return None
    # A corrected record may retain the owning row as ``node`` while the
    # actual edge target lives under ``ref``.  Prefer explicit target fields
    # and nested references before falling back to generic node fields.
    for key in ("target_node", "target_node_id", "target"):
        item = value.get(key)
        if type(item) is int and item >= 0:
            return item
    for key in ("ref", "target_ref", "reference"):
        nested = value.get(key)
        node = _node_from_value(nested)
        if node is not None:
            return node
    for key in ("node", "node_id", "source_node_id"):
        item = value.get(key)
        if type(item) is int and item >= 0:
            return item
    return None


def _validate_edges(value: Mapping[str, Any]) -> dict[str, Any]:
    """Cross-check the top-level review, row reviews, and corrected refs."""

    top_reviews = value.get("edges_reviewed")
    _require(isinstance(top_reviews, list), "top-level edges_reviewed is missing")
    _require(
        len(top_reviews) == EXPECTED_DIRECT_EDGE_COUNT,
        "top-level edges_reviewed must contain all six prior direct edges",
    )
    _require(all(isinstance(item, Mapping) for item in top_reviews),
             "top-level edges_reviewed contains a non-object")

    rows = value.get("rows")
    _require(isinstance(rows, Mapping), "corrected rows are missing")
    _require(set(rows) == set(ROW_KEYS), "corrected row set changed")

    row_reviews: list[Mapping[str, Any]] = []
    row_summaries: list[dict[str, Any]] = []
    target_one_count = 0
    protected_count = 0
    top_corrected = value.get("corrected_typed_edge_refs")
    _require(isinstance(top_corrected, Mapping),
             "top-level corrected_typed_edge_refs is missing")
    _require(set(top_corrected) == set(ROW_KEYS),
             "top-level corrected row keys changed")

    for row_key in ROW_KEYS:
        row = rows.get(row_key)
        _require(isinstance(row, Mapping), "corrected row is missing: " + row_key)
        reviews = row.get("edges_reviewed")
        _require(isinstance(reviews, list), "row edges_reviewed is missing: " + row_key)
        _require(reviews, "row edges_reviewed is empty: " + row_key)
        _require(all(isinstance(item, Mapping) for item in reviews),
                 "row edges_reviewed contains a non-object: " + row_key)
        row_reviews.extend(reviews)

        corrected = row.get("typed_edge_refs")
        corrected_duplicate = row.get("corrected_typed_edge_refs")
        prior = row.get("prior_typed_edge_refs")
        _require(isinstance(corrected, list),
                 "corrected typed refs are missing: " + row_key)
        _require(isinstance(corrected_duplicate, list),
                 "corrected typed-ref duplicate is missing: " + row_key)
        _require(isinstance(prior, list),
                 "prior typed refs are missing: " + row_key)
        _require(len(corrected) == len(reviews) == len(prior),
                 "typed ref counts do not match review: " + row_key)
        _require(corrected == corrected_duplicate,
                 "corrected typed-ref views differ: " + row_key)
        _require(top_corrected[row_key] == corrected,
                 "top-level corrected refs differ: " + row_key)

        scope = row.get("typed_edge_scope")
        _require(isinstance(scope, Mapping), "corrected typed scope is missing: " + row_key)
        _require(scope.get("refs") == corrected,
                 "typed scope does not point to corrected refs: " + row_key)
        _require(scope.get("count") == len(corrected),
                 "corrected typed scope count changed: " + row_key)

        for index, (review, edge, old_edge) in enumerate(zip(reviews, corrected, prior)):
            label = row_key + ".edge[" + str(index) + "]"
            _require(isinstance(review, Mapping), "review is not an object: " + label)
            _require(isinstance(edge, Mapping), "corrected edge is not an object: " + label)
            _require(isinstance(old_edge, Mapping), "prior edge is not an object: " + label)
            _require(review.get("row") == row_key, "review row changed: " + label)
            _require(review.get("edge_index") == index, "review index changed: " + label)
            field = review.get("field")
            node = review.get("node")
            _require(isinstance(field, str), "review field is missing: " + label)
            _require(type(node) is int and node >= 0, "review target is missing: " + label)
            _require(review.get("source_node_id") == node,
                     "review source target changed: " + label)
            previous = review.get("previous_protected_target")
            protected = review.get("protected_target")
            _require(isinstance(previous, bool), "review prior flag is missing: " + label)
            _require(isinstance(protected, bool), "review protected flag is missing: " + label)
            _require(review.get("corrected_protected_target") is protected,
                     "review corrected flag differs: " + label)
            _require(protected is (node in PROTECTED_NODE_IDS),
                     "review protection mismatch: " + label)
            _require(review.get("indexed") is False, "review indexed flag is true: " + label)
            _require(review.get("edge_followed") is False,
                     "review edge-follow flag is true: " + label)

            _require(edge.get("field") == field, "corrected field differs: " + label)
            _require(edge.get("typed_edge") is True, "corrected edge is not typed: " + label)
            _require(edge.get("child_indexed") is review.get("indexed"),
                     "corrected indexed flag differs: " + label)
            _require(edge.get("edge_followed") is review.get("edge_followed"),
                     "corrected edge-follow flag differs: " + label)
            _require(edge.get("protected_target") is protected,
                     "corrected protected flag differs: " + label)
            ref = edge.get("ref")
            _require(isinstance(ref, Mapping), "corrected ref is missing: " + label)
            _require(_node_from_value(ref) == node, "corrected target differs: " + label)
            _require(ref.get("source_node_id") == node,
                     "corrected ref source target differs: " + label)
            _require(ref.get("indexed") is False, "corrected ref indexed: " + label)
            _require(ref.get("protected_node_indexed") is False,
                     "corrected protected ref indexed: " + label)
            _require(ref.get("source_value_selected") is False,
                     "corrected ref selected a source value: " + label)

            # The old list is checked separately and is never used as the
            # corrected list.  This preserves the stale False marker while
            # requiring the newly emitted list to carry the relabel.
            _require(old_edge.get("field") == field, "prior field differs: " + label)
            _require(_node_from_value(old_edge.get("ref")) == node,
                     "prior target differs: " + label)
            _require(old_edge.get("protected_target") is previous,
                     "prior protection flag was not retained: " + label)
            _require(old_edge.get("child_indexed") is False,
                     "prior edge indexed a child: " + label)
            _require(old_edge.get("edge_followed") is False,
                     "prior edge was followed: " + label)

            if node == TARGET_NODE:
                target_one_count += 1
                _require(previous is False and protected is True,
                         "target node 1 was not relabeled from false to true")
                protected_count += 1
            row_summaries.append(
                {
                    "row": row_key,
                    "edge_index": index,
                    "target_node": node,
                    "protected_target": protected,
                    "expected_protected": node in PROTECTED_NODE_IDS,
                    "indexed": False,
                    "edge_followed": False,
                }
            )

    _require(row_reviews == top_reviews,
             "top-level edges_reviewed differs from row reviews")
    _require(target_one_count == 1,
             "correction result must contain exactly one prior edge targeting node 1")
    declared = value.get("protected_node_ids")
    if declared is not None:
        _require(
            isinstance(declared, (list, tuple, set, frozenset))
            and {item for item in declared if type(item) is int} == set(PROTECTED_NODE_IDS)
            and all(type(item) is int for item in declared),
            "correction protected-node set changed",
        )
    return {
        "edge_count": len(row_reviews),
        "target_node_1_corrected": True,
        "target_node_1_count": target_one_count,
        "protected_edge_count": protected_count,
        "protected_node_ids": sorted(PROTECTED_NODE_IDS),
        "edge_summaries": row_summaries,
    }


def _validate_source_result(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "correction leaf did not return an object")
    # ``original_result`` is an authenticated, retained prior receipt.  It
    # legitimately carries that receipt's old archive-decode marker; checking
    # it as if it described this correction would turn preservation into a
    # new archive read.  Validate the newly emitted correction scope and keep
    # the retained prior packet opaque.
    correction_scope = dict(value)
    correction_scope.pop("original_result", None)
    _validate_false_scope(correction_scope, "correction leaf")
    _require(
        value.get("source_only") is True
        or value.get("structural_only") is True
        or value.get("correction_only") is True
        or value.get("protected_ref_correction_only") is True,
        "correction leaf is not source/structural-only",
    )
    return _validate_edges(value)


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
            "manifest_input_count": len(hashes),
            "manifest_paths": list(REQUIRED_INPUTS),
            "prior_result_path": PRIOR_RESULT.as_posix(),
            "prior_result_sha256": hashes[PRIOR_RESULT.as_posix()],
            "new_leaf_imported": False,
            "consumer_calls_started": 0,
            "synthetic_fixture_invocations": 0,
            "source_consumer_invocations": 0,
            "prior_result_parsed": False,
            "archive_read": False,
            "old_code_imported": False,
            "graph_traversal": False,
            "numeric_evaluation": False,
            "protected_node_ids": sorted(PROTECTED_NODE_IDS),
            "protected_nodes_indexed": False,
            "target_node_indexed": False,
            "source_identity_proved": False,
            "derivative_proved": False,
            "positivity_proved": False,
            "error_removal": False,
            "closure_proved": False,
            "runtime_model_metadata": None,
        },
    )

    steps: list[dict[str, Any]] = []
    fixture_result: Any = None
    source_result: Any = None
    edge_metrics: dict[str, Any] | None = None
    leaf_imported = False
    try:
        module = _import_leaf()
        leaf_imported = True

        fixture = _callable(module, "verify_synthetic_fixture")
        steps.append({"step": "verify_synthetic_fixture", "invoked": True})
        fixture_result = fixture()
        steps[-1]["returned_packet"] = _recordable(fixture_result)
        _require_fixture(fixture_result)

        builder = _callable(module, "build_protected_ref_correction")
        steps.append({"step": "build_protected_ref_correction", "invoked": True})
        source_result = builder(ROOT, source, hashes)
        steps[-1]["returned_packet"] = _recordable(source_result)
        edge_metrics = _validate_source_result(source_result)

        output: dict[str, Any] = {
            "schema": COMPONENT + "_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "source_checkpoint_sha256": source_sha,
            "entry_sha256": entry_sha,
            "input_hashes": hashes,
            "manifest_input_count": len(hashes),
            "manifest_paths": list(REQUIRED_INPUTS),
            "prior_result_path": PRIOR_RESULT.as_posix(),
            "prior_result_sha256": hashes[PRIOR_RESULT.as_posix()],
            "consumer_calls": [item["step"] for item in steps],
            "verification_invocations": len(steps),
            "synthetic_fixture_invocations": 1,
            "source_consumer_invocations": 1,
            "leaf_imported": leaf_imported,
            "prior_result_parsed": True,
            "archive_read": False,
            "old_code_imported": False,
            "graph_traversal": False,
            "numeric_evaluation": False,
            "protected_node_ids": sorted(PROTECTED_NODE_IDS),
            "protected_nodes_indexed": False,
            "target_node_indexed": False,
            "source_identity_proved": False,
            "derivative_proved": False,
            "positivity_proved": False,
            "error_removal": False,
            "closure_proved": False,
            "edge_metrics": _json_safe(edge_metrics),
            "result": {
                "synthetic_fixture": _recordable(fixture_result),
                "protected_ref_correction": _recordable(source_result),
            },
            "source_manifest": _json_safe(source),
            "steps": _json_safe(steps),
            "error": None,
            "runtime_model_metadata": None,
        }
        scope_output = dict(output)
        # ``steps`` repeats the full returned packet, including the retained
        # prior receipt, so it is intentionally outside this scope walk.
        scope_output.pop("steps", None)
        scope_result = dict(output["result"])
        if isinstance(source_result, Mapping):
            correction_scope = dict(source_result)
            correction_scope.pop("original_result", None)
            scope_result["protected_ref_correction"] = correction_scope
        scope_output["result"] = scope_result
        _validate_false_scope(scope_output, "runner success result")
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
            "manifest_input_count": len(hashes),
            "manifest_paths": list(REQUIRED_INPUTS),
            "prior_result_path": PRIOR_RESULT.as_posix(),
            "prior_result_sha256": hashes[PRIOR_RESULT.as_posix()],
            "consumer_calls": [item["step"] for item in steps],
            "verification_invocations": len(steps),
            "synthetic_fixture_invocations": sum(
                item.get("step") == "verify_synthetic_fixture" for item in steps
            ),
            "source_consumer_invocations": sum(
                item.get("step") == "build_protected_ref_correction" for item in steps
            ),
            "leaf_imported": leaf_imported,
            "prior_result_parsed": any(
                item.get("step") == "build_protected_ref_correction" for item in steps
            ),
            "archive_read": False,
            "old_code_imported": False,
            "graph_traversal": False,
            "numeric_evaluation": False,
            "protected_node_ids": sorted(PROTECTED_NODE_IDS),
            "protected_nodes_indexed": False,
            "target_node_indexed": False,
            "source_identity_proved": False,
            "derivative_proved": False,
            "positivity_proved": False,
            "error_removal": False,
            "closure_proved": False,
            "edge_metrics": _recordable(edge_metrics),
            "result": {
                "synthetic_fixture": _recordable(fixture_result),
                "protected_ref_correction": _recordable(source_result),
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

    # Persist the complete packet before exposing status or result hash.
    result_sha = _write_once(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "result_sha256": result_sha,
                "entry_sha256": entry_sha,
                "source_checkpoint_sha256": source_sha,
                "input_hashes": hashes,
                "manifest_input_count": len(hashes),
                "manifest_paths": list(REQUIRED_INPUTS),
                "prior_result_path": PRIOR_RESULT.as_posix(),
                "prior_result_sha256": hashes[PRIOR_RESULT.as_posix()],
                "consumer_calls": output["consumer_calls"],
                "synthetic_fixture_invocations": output["synthetic_fixture_invocations"],
                "source_consumer_invocations": output["source_consumer_invocations"],
                "target_node_1_corrected": (
                    edge_metrics.get("target_node_1_corrected") is True
                    if isinstance(edge_metrics, Mapping)
                    else False
                ),
                "protected_nodes_indexed": output["protected_nodes_indexed"],
                "archive_read": output["archive_read"],
                "graph_traversal": output["graph_traversal"],
                "numeric_evaluation": output["numeric_evaluation"],
                "source_identity_proved": output["source_identity_proved"],
                "derivative_proved": output["derivative_proved"],
                "positivity_proved": output["positivity_proved"],
                "error_removal": output["error_removal"],
                "closure_proved": output["closure_proved"],
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
