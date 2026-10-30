"""Run the bounded F pressure-radius path consumer exactly once.

This wrapper is a persistence boundary for
``f_pressure_radius_path_v1``.  It authenticates every raw file named by the
pending manifest before the exclusive entry receipt, then imports only the
new leaf.  The leaf may decode the authenticated F archive once and traverse
at most 4096 old-F nodes per declared root.  This is a new bounded path
query; it is not an ancestor replay, an integral, or a source-value
evaluation.

The wrapper deliberately accepts non-hit and unknown root outcomes.  Those
outcomes remain in the returned packet and are never rewritten as
``no_dependence``.  Saved/protected/terminal references remain opaque, and
the semantic gates for F-radius identity, capture, ordinary-Z transport,
closure, and numeric use are required to remain false.
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
COMPONENT = "f_pressure_radius_path_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"

RUNNER = Path("scripts/run_f_pressure_radius_path_once.py")
LEAF = Path("src/openai_ns_reconstruction/f_pressure_radius_path.py")
F_ARCHIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
F_ARCHIVE_CHECK = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source_check.json"
)
F_SOURCE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.py"
)
PRE_RP_SOURCE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.py"
)
PRIOR_RADIUS_PATH_RESULT = Path(
    "tmp/f_gap_saved_radius_reachability_v2/result.json"
)

F_ARCHIVE_SHA256 = (
    "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"
)
F_ARCHIVE_CHECK_SHA256 = (
    "320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15"
)
F_SOURCE_SHA256 = (
    "369e6404d885dd4548e4f6451f2d5a5ed4160c0abb79a4aca2849520ede03541"
)
PRE_RP_SOURCE_SHA256 = (
    "661933739091ce43f055813510aad016756c2094adf6a6ea2e652def7f129a2e"
)
PRIOR_RADIUS_PATH_RESULT_SHA256 = (
    "5d0ba3e33e9b54941da7f154fb638d3b069554572eaaac55d23391a42a0e31b3"
)

REQUIRED_INPUTS = (
    LEAF.as_posix(),
    RUNNER.as_posix(),
    F_ARCHIVE.as_posix(),
    F_ARCHIVE_CHECK.as_posix(),
    F_SOURCE.as_posix(),
    PRE_RP_SOURCE.as_posix(),
    PRIOR_RADIUS_PATH_RESULT.as_posix(),
)

EXPECTED_FIXED_SHA256 = {
    F_ARCHIVE.as_posix(): F_ARCHIVE_SHA256,
    F_ARCHIVE_CHECK.as_posix(): F_ARCHIVE_CHECK_SHA256,
    F_SOURCE.as_posix(): F_SOURCE_SHA256,
    PRE_RP_SOURCE.as_posix(): PRE_RP_SOURCE_SHA256,
    PRIOR_RADIUS_PATH_RESULT.as_posix(): PRIOR_RADIUS_PATH_RESULT_SHA256,
}

SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")

# These keys are safety gates, not semantic inputs.  A returned leaf may use a
# slightly different nesting shape, so the guard walks mappings and lists.
FORBIDDEN_TRUE_KEYS = frozenset(
    {
        "protected_dereferenced",
        "protected_ref_dereferenced",
        "protected_frontier_dereferenced",
        "protected_indexed",
        "protected_root_indexed",
        "protected_nodes_indexed",
        "terminal_radius_indexed",
        "terminal_radius_dereferenced",
        "radius_terminal_indexed",
        "radius_terminal_dereferenced",
        "radius_row_indexed",
        "saved_boundary_dereferenced",
        "saved_refs_dereferenced",
        "saved_reference_dereferenced",
        "protected_reference_dereferenced",
        "target_row_dereferenced",
        "target_node_dereferenced",
        "target_ref_dereferenced",
        "saved_target_dereferenced",
        "operator_row_dereferenced",
        "F_radius_identity_proved",
        "f_radius_identity_proved",
        "field_radius_identity_proved",
        "actual_field_radius_identity_proved",
        "F_radius_source_identity_proved",
        "F_field_radius_owner_bound",
        "f_field_radius_owner_bound",
        "field_radius_owner_bound",
        "F_radius_owner_bound",
        "f_radius_owner_bound",
        "capture_free",
        "capture_free_proved",
        "capture_free_substitution_proved",
        "capture_free_endpoint_substitution_proved",
        "ordinary_Z_proved",
        "ordinary_Z_transport_proved",
        "ordinary_z_transport_proved",
        "ordinary_Z1_transport_proved",
        "m15_closure_proved",
        "five_moment_closure_proved",
        "closure_proved",
        "source_values_computed",
        "source_value_selected",
        "source_values_selected",
        "numeric_values_computed",
        "numeric_use_admissible",
        "integrals_evaluated",
        "m15_integral_evaluated",
        "quadrature_evaluated",
        "old_source_graph_traversed",
        "no_dependence",
        "no_dependency",
        "dependency_absent",
        "radius_not_used",
        "mathematical_nonzero_dependence_proved",
    }
)

REQUIRED_FALSE_KEYS = (
    "F_radius_identity_proved",
    "capture_free",
    "ordinary_Z_transport_proved",
    "m15_closure_proved",
    "numeric_use_admissible",
)


class RunnerInputError(ValueError):
    """Raised when the frozen one-shot contract is not satisfied."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha256_file(path: Path) -> str:
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
    """Keep a returned fixture/path packet JSON-safe without interpreting it."""

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
    """Hash every declared raw byte before the exclusive entry is written."""

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
        "hashes",
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
    _require(not missing, "manifest omits required frozen input(s): " + ", ".join(missing))
    return hashes


def _checkpoint_sha(brief: Mapping[str, Any]) -> str:
    checkpoint = brief.get(f"{COMPONENT}_source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        checkpoint = brief.get("source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "pending brief source checkpoint missing")
    raw_path = checkpoint.get("path")
    expected = checkpoint.get("sha256")
    _require(isinstance(raw_path, str) and raw_path, "brief checkpoint path missing")
    candidate = Path(raw_path)
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    _require(candidate.resolve() == SOURCE.resolve(), "brief checkpoint points elsewhere")
    _require(
        isinstance(expected, str)
        and SHA256_RE.fullmatch(expected) is not None
        and _sha256_file(SOURCE).lower() == expected.lower(),
        "source checkpoint SHA changed",
    )
    return _sha256_file(SOURCE)


def _guard_true(source: Mapping[str, Any], names: tuple[str, ...], label: str) -> None:
    present = [name for name in names if name in source]
    _require(present, "manifest guard missing: " + label)
    _require(any(source[name] is True for name in present), "manifest guard is not true: " + label)


def _preflight() -> tuple[dict[str, Any], dict[str, str], str]:
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
        ("all_raw_inputs_authenticated_before_decompression", "all_raw_inputs_authenticated_before_decode"),
        "raw authentication",
    )
    _guard_true(
        source,
        ("bounded_old_F_graph_traversal_authorized", "bounded_F_graph_traversal_authorized"),
        "bounded F traversal authorization",
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
    max_nodes = source.get("max_nodes_per_root")
    if max_nodes is None:
        max_nodes = source.get("bounded_max_nodes_per_root")
    _require(type(max_nodes) is int and max_nodes <= 4096 and max_nodes > 0,
             "manifest max_nodes_per_root must be in 1..4096")
    hashes = _manifest_hashes(source)
    _require(LEAF.as_posix() in hashes, "new F pressure-radius leaf is not pinned")
    return dict(source), hashes, source_sha


def _import_leaf() -> Any:
    path = (ROOT / LEAF).resolve()
    _require(path.is_file(), "new F pressure-radius leaf is missing")
    spec = importlib.util.spec_from_file_location(
        "f_pressure_radius_path_v1_new_leaf", path
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


def _require_scope_false(value: Mapping[str, Any], label: str) -> None:
    forbidden = _walk_true_keys(value)
    _require(not forbidden, label + " set forbidden true scope: " + ", ".join(forbidden))
    # If the leaf emits these canonical gates, require them explicitly false;
    # absent optional aliases are not fabricated by the runner.
    for key in REQUIRED_FALSE_KEYS:
        if key in value:
            _require(value[key] is False, label + " gate changed: " + key)


def _require_fixture(value: Any) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), "synthetic fixture did not return an object")
    checks = value.get("checks")
    passed = (
        value.get("passed") is True
        or value.get("checks_passed") is True
        or value.get("fixture_passed") is True
        or (isinstance(checks, Mapping) and bool(checks) and all(item is True for item in checks.values()))
        or "passed" in str(value.get("status", "")).lower()
    )
    _require(passed, "synthetic fixture did not pass")
    _require_scope_false(value, "synthetic fixture")
    _json_bytes(_json_safe(value))
    return value


def _root_count(value: Mapping[str, Any], root: Mapping[str, Any], index: int) -> int | None:
    for key in (
        "visited_node_count",
        "traversed_node_count",
        "nodes_traversed",
        "nodes_visited",
        "bounded_node_count",
        "visit_count",
    ):
        candidate = root.get(key)
        if type(candidate) is int:
            return candidate
    return None


def _validate_bounded_roots(value: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the new bounded walk without requiring any root to hit."""

    roots = value.get("roots")
    if roots is None:
        roots = value.get("root_paths")
    if roots is None:
        roots = value.get("per_root")
    _require(roots is not None, "bounded F path result has no per-root packet")
    if isinstance(roots, Mapping):
        root_items = list(roots.values())
    elif isinstance(roots, list):
        root_items = list(roots)
    else:
        raise RunnerInputError("bounded F path roots are not a list or mapping")

    _require(len(root_items) == 2, "bounded F path must return exactly two roots")
    root_scope = [
        (root.get("role"), root.get("row_order"))
        if isinstance(root, Mapping)
        else (None, None)
        for root in root_items
    ]
    _require(
        root_scope == [
            ("F.pressure.m15", 0),
            ("F.pressure.m15", 1),
        ],
        "bounded F path roots must be F.pressure.m15 row orders 0 and 1",
    )

    counts: list[dict[str, Any]] = []
    for index, root in enumerate(root_items):
        _require(isinstance(root, Mapping), f"bounded root {index} is malformed")
        count = _root_count(value, root, index)
        _require(type(count) is int and 0 <= count <= 4096,
                 f"bounded root {index} exceeds 4096 or has no visit count")
        existence = root.get("existence")
        _require(type(existence) is bool,
                 f"bounded root {index} existence must be a boolean")
        unresolved = root.get("unresolved")
        status = str(root.get("status", ""))
        unknown = status.startswith("unresolved") or (
            isinstance(unresolved, list) and bool(unresolved)
        )
        edge_counts = root.get("edge_counts")
        protected_terminations = (
            edge_counts.get("protected_terminations", 0)
            if isinstance(edge_counts, Mapping)
            else root.get("protected_terminations", 0)
        )
        hit_with_protected_boundary = bool(
            existence and (
                status == "unresolved_protected_boundary"
                or (type(protected_terminations) is int and protected_terminations > 0)
            )
        )
        counts.append(
            {
                "index": index,
                "role": root.get("role"),
                "row_order": root.get("row_order"),
                "visited_node_count": count,
                "hit": existence,
                "unknown": unknown,
                "hit_with_protected_boundary": hit_with_protected_boundary,
            }
        )

    max_nodes = value.get("max_nodes_per_root")
    if max_nodes is None:
        max_nodes = value.get("bounded_max_nodes_per_root")
    archive_scope = value.get("archive_scope")
    if max_nodes is None and isinstance(archive_scope, Mapping):
        max_nodes = archive_scope.get("max_nodes_per_root")
    if max_nodes is not None:
        _require(type(max_nodes) is int and 0 < max_nodes <= 4096,
                 "leaf max_nodes_per_root exceeds 4096")

    scope = value.get("scope")
    archive_scope = value.get("archive_scope")
    old_f_traversed = (
        value.get("old_F_graph_traversed") is True
        or (isinstance(scope, Mapping) and scope.get("F_archive_exact_graph_traversed") is True)
        or (isinstance(archive_scope, Mapping) and archive_scope.get("F_archive_exact_graph_traversed") is True)
    )
    _require(old_f_traversed, "leaf did not record bounded old-F graph traversal")
    _require(isinstance(scope, Mapping), "F pressure-radius path scope is missing")
    _require(scope.get("ancestor_replayed") is False, "leaf replayed an ancestor")
    _require(scope.get("old_producer_replayed") is False, "leaf replayed the old producer")
    _require(scope.get("fresh_graph_generated") is False, "leaf generated a fresh graph")
    _require(scope.get("prior_20_root_query_reentered") is False,
             "leaf reentered the prior twenty-root query")
    _require(scope.get("target_row_dereferenced") is False,
             "leaf dereferenced the terminal radius row")
    _require(scope.get("protected_root_indexed") is False,
             "leaf indexed a protected root")
    return {
        "root_count": len(counts),
        "visited_node_counts_by_root": counts,
        "max_nodes_per_root": max_nodes,
        "old_F_graph_traversed": old_f_traversed,
        "nonhit_or_unknown_preserved": True,
    }


def _validate_path_result(value: Any) -> tuple[Mapping[str, Any], dict[str, Any]]:
    _require(isinstance(value, Mapping), "F pressure-radius path did not return a mapping")
    _require_scope_false(value, "F pressure-radius path")
    metrics = _validate_bounded_roots(value)
    # False and unknown roots are permitted, as are positive structural hits.
    # The only forbidden conversion is a negative/unknown witness into an
    # asserted no-dependence result, guarded recursively above.
    return value, metrics


def _base_envelope(
    *,
    status: str,
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    source_sha: str,
    entry_sha: str | None,
    steps: list[dict[str, Any]],
    fixture_result: Any = None,
    path_result: Any = None,
    path_metrics: Mapping[str, Any] | None = None,
    error: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    flags = {
        "source_values_computed": False,
        "source_value_selected": False,
        "numeric_values_computed": False,
        "numeric_use_admissible": False,
        "integrals_evaluated": False,
        "quadrature_evaluated": False,
        "old_source_graph_traversed": False,
        "source_graph_traversed": False,
        "old_F_graph_traversed": bool(path_metrics and path_metrics.get("old_F_graph_traversed")),
        "protected_dereferenced": False,
        "terminal_radius_indexed": False,
        "F_radius_identity_proved": False,
        "F_radius_owner_bound": False,
        "capture_free": False,
        "ordinary_Z_transport_proved": False,
        "m15_closure_proved": False,
        "closure_proved": False,
        "no_dependence": False,
        "full_goal_complete": False,
    }
    output: dict[str, Any] = {
        "schema": COMPONENT + "_result_v1",
        "component": COMPONENT,
        "phase": PHASE,
        "status": status,
        "flags": flags,
        **flags,
        "source_checkpoint_sha256": source_sha,
        "entry_sha256": entry_sha,
        "required_inputs": list(REQUIRED_INPUTS),
        "input_hashes": dict(hashes),
        "raw_sha256_authenticated_before_entry": True,
        "old_F_graph_traversal": {
            "authorized": True,
            "max_nodes_per_root": 4096,
            "new_query": True,
            "ancestor_replay": False,
            "integrals": False,
        },
        "steps": steps,
        "synthetic_fixture": _json_safe(fixture_result),
        "f_pressure_radius_path": _json_safe(path_result),
        "path_metrics": _json_safe(path_metrics),
        "source_manifest": _json_safe(source),
        "error": _json_safe(error),
    }
    return output


def main() -> int:
    # No failure packet is written before entry: malformed/untrusted pending
    # metadata must not consume the one-entry slot.
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
            "old_F_graph_traversed": False,
            "max_nodes_per_root": 4096,
            "protected_dereferenced": False,
            "terminal_radius_indexed": False,
            "source_values_computed": False,
            "integrals_evaluated": False,
            "numeric_use_admissible": False,
        },
    )

    steps: list[dict[str, Any]] = []
    fixture_result: Any = None
    path_result: Any = None
    path_metrics: dict[str, Any] | None = None
    try:
        module = _import_leaf()
        fixture = _callable(module, "verify_synthetic_fixture")
        fixture_result = _invoke(steps, "verify_synthetic_fixture", fixture)
        _require_fixture(fixture_result)

        builder = _callable(module, "build_f_pressure_radius_path")
        path_result = _invoke(
            steps,
            "build_f_pressure_radius_path",
            builder,
            ROOT,
            source,
            hashes,
        )
        path_result, path_metrics = _validate_path_result(path_result)
        output = _base_envelope(
            status="returned",
            source=source,
            hashes=hashes,
            source_sha=source_sha,
            entry_sha=entry_sha,
            steps=steps,
            fixture_result=fixture_result,
            path_result=path_result,
            path_metrics=path_metrics,
        )
        result_sha = _write_once(RESULT, output)
        print(
            json.dumps(
                {
                    "component": COMPONENT,
                    "phase": PHASE,
                    "status": output["status"],
                    "result_sha256": result_sha,
                    "consumer_calls": [item["step"] for item in steps],
                    "old_F_graph_traversed": output["old_F_graph_traversed"],
                    "max_nodes_per_root": 4096,
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
        output = _base_envelope(
            status="failed",
            source=source,
            hashes=hashes,
            source_sha=source_sha,
            entry_sha=entry_sha,
            steps=steps,
            fixture_result=fixture_result,
            path_result=path_result,
            path_metrics=path_metrics,
            error=error,
        )
        result_sha = _write_once(RESULT, output)
        print(
            json.dumps(
                {
                    "component": COMPONENT,
                    "phase": PHASE,
                    "status": output["status"],
                    "result_sha256": result_sha,
                    "error_type": error["type"],
                    "consumer_calls": [item["step"] for item in steps],
                    "old_F_graph_traversed": output["old_F_graph_traversed"],
                    "max_nodes_per_root": 4096,
                },
                sort_keys=True,
            )
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
