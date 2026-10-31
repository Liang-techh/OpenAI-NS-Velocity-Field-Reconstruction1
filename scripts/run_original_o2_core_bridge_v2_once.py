"""Run the frozen original O2-to-core bridge exactly once.

The runner authenticates the pending brief and every manifest byte before it
writes an exclusive entry receipt.  The synthetic fixture is then run against
the new leaf without reading any real report or decoding the previous
receipt.  Only after that fixture succeeds are the two saved JSON reports,
the compressed core-lineage receipt, and the two producer source texts passed
to the leaf.

The complete leaf comparison is retained.  ``proved`` is reserved for the
leaf's own proved status; ``mismatch`` and ``incomplete`` remain evidence and
are never relabeled by this runner.  The runner does not import or execute a
producer source, does not claim live endpoint or native/corrected identity,
and does not promote interval overlap to function identity.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from collections.abc import Mapping
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "original_o2_core_bridge_v2"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entry.json"
RESULT = DIRECTORY / "result.json"

O2_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_O2_source_parameter_frame_check.json"
)
CORE_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json"
)
LINEAGE_RECEIPT = (
    "experiments/root_st073/original_core_parameter_lineage_v1_receipt.json.gz"
)
CORE_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py"
)
BRIDGE_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source.py"
)
LEAF = "src/openai_ns_reconstruction/original_o2_core_bridge_v2.py"
FIXTURE = "src/openai_ns_reconstruction/original_o2_core_bridge_fixture_v2.py"
RUNNER = "scripts/run_original_o2_core_bridge_v2_once.py"

# These are the immutable pre-existing inputs.  The source manifest also
# freezes the new leaf, fixture, and this runner, while this table prevents a
# later run from silently substituting any archived O2/core input.
FIXED = {
    O2_CHECK: "93321a4d7b26bab5ed3b772448b64e5b0a8c207468d5d401c4b4cb2848134c5c",
    CORE_CHECK: "45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead",
    LINEAGE_RECEIPT: "9b86acdd51dce692d7a9fd34687be0b326662ce3add75d71bd356cbd3bf6e6b1",
    CORE_SOURCE: "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59",
    BRIDGE_SOURCE: "df8b204d023c1d8762da16d07190601ebd1e11447077b01a734bba8baa595a6a",
}
REQUIRED = (
    O2_CHECK,
    CORE_CHECK,
    LINEAGE_RECEIPT,
    CORE_SOURCE,
    BRIDGE_SOURCE,
    LEAF,
    FIXTURE,
    RUNNER,
)
NEW_CODE = (LEAF, FIXTURE, RUNNER)
REPORTS = (O2_CHECK, CORE_CHECK)
PRODUCER_SOURCES = (CORE_SOURCE, BRIDGE_SOURCE)


class RunnerInputError(ValueError):
    """Raised when the frozen one-entry contract is not satisfied."""


def require(condition: object, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _json_constant(token: str) -> Any:
    raise ValueError("non-finite JSON constant is not admitted: " + token)


def load_json(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(raw.decode("utf-8"), parse_constant=_json_constant)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise RunnerInputError(label + " is not finite JSON") from exc
    require(isinstance(value, dict), label + " must be a JSON object")
    return value


def save(path: Path, value: object) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
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
    try:
        with path.open("xb") as handle:
            handle.write(raw)
    except FileExistsError as exc:
        raise RunnerInputError("one-entry receipt already exists: " + str(path)) from exc
    return sha(raw)


def repository_path(relative: str) -> Path:
    path = Path(relative)
    require(not path.is_absolute() and ".." not in path.parts, "non-relative input path")
    resolved = (ROOT / path).resolve()
    require(
        resolved.is_relative_to(ROOT.resolve()),
        "input path escapes repository: " + relative,
    )
    require(resolved.is_file(), "missing input: " + relative)
    return resolved


def manifest_hashes(manifest: Mapping[str, Any]) -> dict[str, str]:
    rows = manifest.get("files")
    require(isinstance(rows, list), "source manifest files list is missing")
    allowed = set(REQUIRED)
    hashes: dict[str, str] = {}
    for index, row in enumerate(rows):
        require(isinstance(row, Mapping), f"source.files[{index}] is not an object")
        path = row.get("path")
        expected = row.get("sha256")
        require(isinstance(path, str) and path in allowed, "unexpected manifest input")
        require(path not in hashes, "duplicate manifest input: " + path)
        require(
            isinstance(expected, str) and len(expected) == 64,
            "invalid manifest SHA: " + str(path),
        )
        actual = sha(repository_path(path).read_bytes())
        require(actual == expected.lower(), "manifest SHA changed: " + path)
        if path in FIXED:
            require(actual == FIXED[path], "fixed SHA changed: " + path)
        hashes[path] = actual
    require(set(REQUIRED) == set(hashes), "manifest must pin exactly the frozen inputs")
    require(all(path in hashes for path in NEW_CODE), "new code files are not all frozen")
    return hashes


def preflight() -> tuple[dict[str, Any], dict[str, str], str, str]:
    require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    require(BRIEF.is_file(), "pending brief is missing")
    require(SOURCE.is_file(), "source manifest is missing")
    brief_raw = BRIEF.read_bytes()
    source_raw = SOURCE.read_bytes()
    brief = load_json(brief_raw, "pending brief")
    source = load_json(source_raw, "source manifest")
    require(
        brief.get("component") == COMPONENT
        and brief.get("phase") == PHASE
        and brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("one_entry_only") is True
        and brief.get("do_not_reenter") is True,
        "pending brief guard",
    )
    require(
        source.get("component") == COMPONENT and source.get("one_entry_only") is True,
        "source manifest guard",
    )
    checkpoint = brief.get("source_checkpoint")
    require(isinstance(checkpoint, Mapping), "pending brief source checkpoint is missing")
    require(
        Path(str(checkpoint.get("path"))).resolve() == SOURCE.resolve(),
        "pending brief points to another source manifest",
    )
    require(
        checkpoint.get("sha256") == sha(source_raw),
        "pending brief source checkpoint SHA changed",
    )
    source_checkpoint = source.get("source_checkpoint")
    if source_checkpoint is not None:
        require(isinstance(source_checkpoint, Mapping), "source checkpoint is not an object")
        require(
            Path(str(source_checkpoint.get("path"))).resolve() == SOURCE.resolve()
            and source_checkpoint.get("sha256") == sha(source_raw),
            "source manifest checkpoint changed",
        )
    hashes = manifest_hashes(source)
    return source, hashes, sha(source_raw), sha(brief_raw)


def authenticated_bytes(hashes: Mapping[str, str]) -> dict[str, bytes]:
    raw: dict[str, bytes] = {}
    for path in hashes:
        data = repository_path(path).read_bytes()
        require(sha(data) == hashes[path], "raw input changed after entry: " + path)
        raw[path] = data
    return raw


def static_api(raw: Mapping[str, bytes], path: str, function_name: str) -> dict[str, Any]:
    """Inspect only the new module's top-level API, without evaluating it."""

    source_text = raw[path].decode("utf-8")
    try:
        tree = ast.parse(source_text, filename=path)
    except (UnicodeDecodeError, SyntaxError) as exc:
        raise RunnerInputError("source API is not parseable: " + path) from exc
    functions = {
        node.name
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    require(
        function_name in functions,
        "required API is missing: " + path + ":" + function_name,
    )
    return {"path": path, "function": function_name, "static_only": True}


def module(relative: str, name: str) -> Any:
    path = repository_path(relative)
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, "module load failed: " + relative)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def decode_receipt(raw: bytes) -> dict[str, Any]:
    try:
        decoded = gzip.decompress(raw)
    except (OSError, EOFError) as exc:
        raise RunnerInputError("core-lineage receipt is not valid gzip") from exc
    return load_json(decoded, LINEAGE_RECEIPT + " decoded JSON")


def validate_fixture(result: Mapping[str, Any]) -> None:
    require(result.get("passed") is True, "synthetic O2-core bridge fixture failed")
    require(
        result.get("synthetic_only") is True
        and result.get("real_sources_read") is False
        and result.get("implementation_imported") is False,
        "synthetic fixture scope changed",
    )
    checks = result.get("checks")
    require(
        isinstance(checks, Mapping) and bool(checks),
        "synthetic fixture checks are missing",
    )
    require(
        all(value is True for value in checks.values()),
        "synthetic fixture did not pass every guard",
    )


def validate_leaf(result: Mapping[str, Any]) -> dict[str, Any]:
    require(
        result.get("schema") == COMPONENT,
        "O2-core bridge leaf schema changed",
    )
    status = result.get("status")
    require(
        status in {"proved", "mismatch", "incomplete"},
        "O2-core bridge source status is invalid",
    )
    require(
        result.get("passed") is (status == "proved"),
        "O2-core bridge passed/status relation changed",
    )
    require(
        result.get("source_execution") is False
        and result.get("real_receipt_opened_here") is False
        and result.get("producer_imported") is False
        and result.get("graph_traversed") is False,
        "O2-core bridge leaf execution scope changed",
    )

    # The bridge may expose additional scope labels, but any known identity
    # claim must remain explicitly false.  Status is the only proof signal.
    scope = result.get("provenance_scope")
    require(isinstance(scope, Mapping), "O2-core bridge provenance scope is missing")
    for key in (
        "cross_run_object_identity_claimed",
        "live_endpoint_basis_identity",
        "native_corrected_RQy_identity",
        "source_function_point_value_identity",
        "interval_overlap_promoted_to_identity",
        "rounding_difference_promoted_to_identity",
    ):
        if key in scope:
            require(scope.get(key) is False, "O2-core bridge overclaimed " + key)

    comparisons = result.get("comparisons")
    if isinstance(comparisons, Mapping):
        require(
            comparisons.get("overlap_used_as_identity") is False,
            "interval overlap was used as identity",
        )
        mismatch_count = comparisons.get("mismatch_count")
    else:
        mismatch_count = None

    return {
        "schema": result["schema"],
        "status": status,
        "passed": status == "proved",
        "mismatch_count": mismatch_count,
        "overlap_used_as_identity": False,
        "live_endpoint_basis_identity": False,
        "native_corrected_RQy_identity": False,
        "old_errors_removed": False,
    }


def consume(
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    raw: Mapping[str, bytes],
    stats: dict[str, Any],
) -> dict[str, Any]:
    leaf_static = static_api(raw, LEAF, "verify")
    fixture_static = static_api(raw, FIXTURE, "run_fixture")
    leaf_module = module(LEAF, "original_o2_core_bridge_v2")
    fixture_module = module(FIXTURE, "original_o2_core_bridge_fixture_v2")
    verify = getattr(leaf_module, "verify", None)
    run_fixture = getattr(fixture_module, "run_fixture", None)
    require(
        callable(verify) and callable(run_fixture),
        "O2-core bridge leaf or fixture API missing",
    )

    # The fixture receives only synthetic mappings and runs before any real
    # JSON report or compressed receipt is decoded.
    fixture_result = run_fixture(verify)
    stats["fixture_invocations"] = 1
    require(isinstance(fixture_result, Mapping), "fixture returned a non-object")
    # Preserve the complete fixture even when a guard fails.
    stats["fixture_result"] = dict(fixture_result)
    validate_fixture(fixture_result)

    # Decode the prior lineage receipt only after the synthetic fixture passes.
    stats["receipt_decode_started"] = True
    lineage_receipt = decode_receipt(raw[LINEAGE_RECEIPT])
    o2_check = load_json(raw[O2_CHECK], O2_CHECK)
    core_check = load_json(raw[CORE_CHECK], CORE_CHECK)
    core_source_text = raw[CORE_SOURCE].decode("utf-8")
    bridge_source_text = raw[BRIDGE_SOURCE].decode("utf-8")
    stats["real_source_reads_started"] = True

    stats["leaf_invocations"] = 1
    leaf_result = verify(
        o2_check,
        core_check,
        lineage_receipt,
        core_source_text,
        bridge_source_text,
    )
    require(isinstance(leaf_result, Mapping), "O2-core bridge leaf returned a non-object")
    leaf_summary = validate_leaf(leaf_result)
    source_status = leaf_summary["status"]
    return {
        "schema": COMPONENT + "_output_v1",
        "api_static": {"leaf": leaf_static, "fixture": fixture_static},
        "fixture": dict(fixture_result),
        "leaf": dict(leaf_result),
        "leaf_summary": leaf_summary,
        "source_comparison_status": source_status,
        "source_comparison_status_preserved": source_status
        in {"proved", "mismatch", "incomplete"},
        "source_comparison_mismatch_is_not_identity_proof": source_status != "proved",
        "report_provenance": {
            path: {"path": path, "sha256": hashes[path]} for path in REPORTS
        },
        "producer_source_provenance": {
            path: {"path": path, "sha256": hashes[path]}
            for path in PRODUCER_SOURCES
        },
        "lineage_receipt_provenance": {
            "path": LINEAGE_RECEIPT,
            "sha256": hashes[LINEAGE_RECEIPT],
        },
        "source_provenance": {
            path: {"path": path, "sha256": hashes[path]} for path in REQUIRED
        },
        "raw_sha_verified_before_receipt_decode": True,
        "receipt_decoded_after_fixture": True,
        "actual_live_endpoint_owner_identity": False,
        "native_corrected_RQy_identity": False,
        "interval_equality_promoted_to_function_identity": False,
        "old_errors_removed": False,
        "ledger_edited": False,
    }


def main() -> int:
    source, hashes, source_sha, brief_sha = preflight()
    entry_sha = save(
        ENTRY,
        {
            "schema": COMPONENT + "_entry_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "brief_sha256": brief_sha,
            "input_hashes": hashes,
            "manifest_input_count": len(hashes),
            "receipt_decode_started": False,
            "fixture_invocations": 0,
            "leaf_invocations": 0,
            "old_errors_removed": False,
        },
    )
    stats: dict[str, Any] = {
        "receipt_decode_started": False,
        "real_source_reads_started": False,
        "fixture_invocations": 0,
        "leaf_invocations": 0,
    }
    outcome: dict[str, Any] = {
        "schema": COMPONENT + "_result_v1",
        "component": COMPONENT,
        "phase": PHASE,
        "status": "failed",
        "source_checkpoint_sha256": source_sha,
        "brief_sha256": brief_sha,
        "entry_sha256": entry_sha,
        "input_hashes": hashes,
        "source_manifest": source,
        "result": None,
        "error": None,
        "old_errors_removed": False,
    }
    try:
        raw = authenticated_bytes(hashes)
        outcome["result"] = consume(source, hashes, raw, stats)
        outcome["status"] = "returned"
        outcome["source_comparison_status"] = outcome["result"][
            "source_comparison_status"
        ]
    except BaseException as exc:
        outcome["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
    outcome.update(stats)
    # Persist the complete outcome before printing any status to the caller.
    result_sha = save(RESULT, outcome)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": outcome["status"],
                "source_comparison_status": outcome.get("source_comparison_status"),
                "entry_sha256": entry_sha,
                "result_sha256": result_sha,
            },
            sort_keys=True,
        )
    )
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())

