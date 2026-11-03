"""Run the saved lower-collar direct-sigma consumer exactly once.

The runner has one deliberately narrow actual scope: the direct lower-collar
``sigma(s)`` occurrence on the saved ``[sc/4, sc/2]`` band.  It does not run
the normalized inner cutoff ``sigma(4*s/sc_quarter - 2)`` and it makes no
owner, live-source, physical ``m15``, or five-moment claim.

The source checkpoint is an ordinary JSON file.  Every frozen byte is hashed
before the exclusive entry, the entry is written before any component import,
and the synthetic fixture runs with its own stub.  Only after that fixture
does this runner parse the pinned JSON and call the new consumer once with an
injected evaluator from the published rational sigma leaf.  No old fixture or
ancestor producer is imported or executed, and no gzip input is accepted.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
import sys
import traceback
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"

COMPONENT = "n1_collar_sigma_actual_domain_v1"
PHASE = COMPONENT + "_pending"
CHECKPOINT_DIR = ROOT / "tmp" / COMPONENT
SOURCE = CHECKPOINT_DIR / "source.json"
ENTRY = CHECKPOINT_DIR / "entry.json"
RESULT = CHECKPOINT_DIR / "result.json"

LEAF_REL = "src/openai_ns_reconstruction/n1_collar_sigma_actual_domain_v1.py"
FIXTURE_REL = (
    "src/openai_ns_reconstruction/n1_collar_sigma_actual_domain_v1_fixture.py"
)
SIGMA_REL = "src/openai_ns_reconstruction/n1_sigma_interval_rational_v1.py"
RUNNER_REL = "scripts/run_n1_collar_sigma_actual_domain_v1_once.py"
SAVED_JSON_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_inner_exit_strict_collar.json"
)
# This is pinned as provenance for the distinction from the normalized inner
# cutoff, but it is never imported, decoded, or executed by this runner.
NORMALIZED_SOURCE_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.py"
)
DIRECT_SOURCE_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.py"
)

FROZEN_INPUTS = (
    LEAF_REL,
    FIXTURE_REL,
    SIGMA_REL,
    RUNNER_REL,
    SAVED_JSON_REL,
    DIRECT_SOURCE_REL,
    NORMALIZED_SOURCE_REL,
)
FROZEN_SET = frozenset(FROZEN_INPUTS)
MANIFEST_SCHEMA = "n1_collar_sigma_actual_domain_v1_manifest"
LEAF_SCHEMA = "n1_collar_sigma_actual_domain_v1"
FIXTURE_SCHEMA = "n1_collar_sigma_actual_domain_v1_fixture"
SIGMA_SCHEMA = "n1_sigma_interval_rational_v1"


class RunnerInputError(ValueError):
    """Raised for a fail-closed input, manifest, or scope violation."""


def _require(condition: object, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_safe(value: Any, seen: set[int] | None = None) -> Any:
    """Serialize the exact result without silently converting rationals."""

    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, Fraction):
        return {
            "type": "Fraction",
            "numerator": value.numerator,
            "denominator": value.denominator,
        }
    if isinstance(value, float):
        return value if math.isfinite(value) else repr(value)
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, bytes):
        return {
            "type": "opaque_bytes",
            "length": len(value),
            "sha256": _sha_bytes(value),
        }

    active = seen if seen is not None else set()
    identity = id(value)
    if identity in active:
        return {"type": "opaque_cycle", "repr": repr(value)}
    active.add(identity)
    try:
        if isinstance(value, Mapping):
            return {str(key): _json_safe(item, active) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [_json_safe(item, active) for item in value]
        if isinstance(value, (set, frozenset)):
            return sorted((_json_safe(item, active) for item in value), key=repr)
        return {
            "type": "opaque_value",
            "class": type(value).__name__,
            "repr": repr(value),
        }
    finally:
        active.discard(identity)


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            _json_safe(value),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _write_exclusive(path: Path, value: Any) -> str:
    """Write one immutable JSON record and return its raw SHA-256."""

    path.parent.mkdir(parents=True, exist_ok=True)
    raw = _json_bytes(value)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0)
    descriptor = os.open(str(path), flags, 0o600)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            descriptor = -1
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
    finally:
        if descriptor != -1:
            os.close(descriptor)
    return _sha_bytes(raw)


def _decode_object(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not ordinary UTF-8 JSON") from exc
    _require(isinstance(value, Mapping), label + " is not a JSON object")
    return dict(value)


def _sha_text(value: Any, label: str) -> str:
    _require(
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdefABCDEF" for character in value),
        "invalid SHA-256 for " + label,
    )
    return value.lower()


def _repo_file(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    relative = Path(raw_path.replace("\\", "/"))
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        label + " path must be repository-relative",
    )
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(
        path.is_relative_to(ROOT.resolve()) and path.is_file(),
        label + " path is missing or escapes repository: " + key,
    )
    return key, path


def _validate_manifest(manifest: Mapping[str, Any]) -> dict[str, str]:
    """Validate the exact seven-file manifest and return its frozen hashes."""

    _require(manifest.get("schema") == MANIFEST_SCHEMA, "manifest schema")
    _require(manifest.get("component") == COMPONENT, "manifest component")
    _require(manifest.get("phase") == PHASE, "manifest phase")
    _require(manifest.get("one_entry_only") is True, "manifest one_entry_only")
    _require(manifest.get("do_not_reenter") is True, "manifest do_not_reenter")
    _require(manifest.get("gzip_forbidden") is True, "manifest gzip_forbidden")
    _require(manifest.get("ordinary_json_source_only") is True, "manifest JSON scope")
    _require(
        manifest.get("actual_scope") == "direct_lower_collar_chi_sigma_s_only",
        "manifest actual scope",
    )

    brief = manifest.get("brief")
    _require(
        brief == {"component": COMPONENT, "phase": PHASE},
        "manifest brief guard",
    )

    contract = manifest.get("input_contract")
    _require(isinstance(contract, Mapping), "manifest input_contract")
    expected_contract = {
        "leaf": LEAF_REL,
        "fixture": FIXTURE_REL,
        "published_sigma_leaf": SIGMA_REL,
        "runner": RUNNER_REL,
        "saved_json_source": SAVED_JSON_REL,
        "direct_lower_collar_source": DIRECT_SOURCE_REL,
        "normalized_cutoff_source": NORMALIZED_SOURCE_REL,
        "fixture_api": "run_fixture(module)",
        "actual_api": "evaluate_actual_collar_sigma(saved_mapping, interval_evaluator)",
        "published_api": "evaluate_sigma_interval(lower, upper)",
    }
    for key, expected in expected_contract.items():
        _require(contract.get(key) == expected, "manifest input_contract." + key)
    _require(
        contract.get("normalized_inner_cutoff_internal_call") is True,
        "manifest normalized-inner internal-call contract",
    )
    _require(
        contract.get("normalized_inner_cutoff_runner_scope") is False,
        "manifest normalized-inner runner scope",
    )

    conditions = manifest.get("conditions")
    _require(
        conditions
        == {
            "fixed_s_in_z": True,
            "fixed_s_over_sc_in_z": True,
            "conditional_not_live_identity": True,
        },
        "manifest condition semantics",
    )

    claims = manifest.get("claims")
    _require(isinstance(claims, Mapping), "manifest claims")
    expected_claims = (
        "direct_lower_collar_sigma_only",
        "normalized_inner_cutoff_claimed",
        "physical_m15",
        "live_owner_binding",
        "source_identity_verified",
        "ancestor_reexecuted",
        "old_fixture_called",
        "gzip_decoded",
    )
    _require(
        all(claims.get(key) is False for key in expected_claims[1:]),
        "manifest forbidden claim",
    )
    _require(claims.get("direct_lower_collar_sigma_only") is True, "manifest direct scope claim")

    rows = manifest.get("files")
    _require(isinstance(rows, list) and len(rows) == len(FROZEN_INPUTS), "manifest file count")
    hashes: dict[str, str] = {}
    for index, row in enumerate(rows):
        _require(isinstance(row, Mapping), f"manifest files[{index}]")
        key, path = _repo_file(row.get("path"), f"manifest files[{index}]")
        _require(key in FROZEN_SET and key not in hashes, "manifest file scope: " + key)
        _require(not key.lower().endswith((".gz", ".gzip")), "gzip input is forbidden")
        digest = _sha_file(path)
        _require(digest == _sha_text(row.get("sha256"), key), "changed input: " + key)
        hashes[key] = digest
    _require(set(hashes) == FROZEN_SET, "manifest missing frozen input")
    _require(not SAVED_JSON_REL.lower().endswith((".gz", ".gzip")), "saved source cannot be gzip")
    _require(SAVED_JSON_REL.lower().endswith(".json"), "saved source must be ordinary JSON")
    return hashes


def _validate_brief(brief: Mapping[str, Any], manifest: Mapping[str, Any], source_sha256: str) -> None:
    _require(brief.get("goal_active") is True, "brief goal_active")
    _require(brief.get("full_goal_complete") is False, "brief full_goal_complete")
    _require(brief.get("one_entry_only") is True, "brief one_entry_only")
    _require(brief.get("do_not_reenter") is True, "brief do_not_reenter")
    guard = manifest["brief"]
    _require(brief.get("component") == guard["component"], "brief component mismatch")
    _require(brief.get("phase") == guard["phase"], "brief phase mismatch")
    checkpoint = brief.get("source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "brief source checkpoint")
    checkpoint_path = Path(str(checkpoint.get("path", "")))
    if not checkpoint_path.is_absolute():
        checkpoint_path = (ROOT / checkpoint_path).resolve()
    else:
        checkpoint_path = checkpoint_path.resolve()
    _require(checkpoint_path == SOURCE.resolve(), "brief checkpoint must target this SOURCE")
    _require(checkpoint.get("sha256") == source_sha256, "brief checkpoint must hash this SOURCE")


def _preflight() -> tuple[dict[str, str], str, str, dict[str, Any]]:
    _require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    brief_raw = BRIEF.read_bytes()
    source_raw = SOURCE.read_bytes()
    brief = _decode_object(brief_raw, "brief")
    manifest = _decode_object(source_raw, "source manifest")
    source_sha256 = _sha_bytes(source_raw)
    _validate_brief(brief, manifest, source_sha256)
    hashes = _validate_manifest(manifest)
    return hashes, source_sha256, _sha_bytes(brief_raw), manifest


def _load_component(relative: str, module_name: str) -> Any:
    _require(relative in {LEAF_REL, FIXTURE_REL, SIGMA_REL}, "unexpected component import")
    spec = importlib.util.spec_from_file_location(module_name, ROOT / relative)
    _require(spec is not None and spec.loader is not None, "component loader missing: " + relative)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def _authenticate_frozen_bytes(hashes: Mapping[str, str]) -> dict[str, bytes]:
    raw: dict[str, bytes] = {}
    for relative in FROZEN_INPUTS:
        data = (ROOT / relative).read_bytes()
        _require(_sha_bytes(data) == hashes[relative], "changed after entry: " + relative)
        raw[relative] = data
    return raw


def _validate_fixture_result(value: Any) -> None:
    _require(isinstance(value, Mapping), "fixture result is not a mapping")
    _require(value.get("schema") == FIXTURE_SCHEMA, "fixture schema")
    _require(value.get("passed") is True, "synthetic fixture did not pass")
    checks = value.get("checks")
    _require(isinstance(checks, Mapping) and all(item is True for item in checks.values()), "fixture checks")
    _require(value.get("synthetic_only") is True, "fixture was not synthetic-only")
    _require(value.get("used_saved_real_mapping") is False, "fixture used real mapping")
    _require(value.get("used_previous_interval_leaf") is False, "fixture used previous leaf")


def _validate_actual_result(value: Any, published_calls: list[tuple[Fraction, Fraction]]) -> None:
    _require(isinstance(value, Mapping), "actual result is not a mapping")
    _require(value.get("schema") == LEAF_SCHEMA + "_result", "actual result schema")
    _require(value.get("status") == "bounded", "actual result status")
    _require(value.get("cutoff_scope") == "direct lower-collar chi sigma(s) only", "actual scope")
    # The consumer's one call performs a second, separately labelled internal
    # normalized-window check.  It is authenticated and retained below, but
    # it is explicitly outside this runner's direct-sigma status scope.
    _require(value.get("normalized_inner_cutoff_claimed") is True, "consumer normalized check missing")
    _require(value.get("normalized_inner_cutoff_node_binding") is False, "normalized cutoff node was bound")
    for key in (
        "source_identity_verified",
        "owner_binding",
        "live_owner_binding",
        "graph_binding",
        "physical_m15",
        "five_moment_closure",
    ):
        _require(value.get(key) is False, "forbidden actual claim: " + key)
    _require(value.get("saved_mapping_only") is True, "actual source scope")
    _require(value.get("sc_value_materialized") is False, "selected sc was materialized")
    _require(value.get("sigma_value") is None, "sigma upper/lower selected as value")
    _require(value.get("d_sigma_ds_value") is None, "derivative bound selected as value")
    _require(value.get("fixed_s_in_z") is True, "direct fixed-s ordinary-Z condition")
    _require(value.get("fixed_s_over_sc_in_z") is True, "normalized fixed-s/sc ordinary-Z condition")
    _require(value.get("sigma_Z") == Fraction(0), "direct conditional sigma_Z")
    _require(value.get("sigma_z") == Fraction(0), "direct conditional sigma_z")
    _require(value.get("normalized_sigma_Z") == Fraction(0), "normalized conditional sigma_Z")
    _require(value.get("normalized_sigma_z") == Fraction(0), "normalized conditional sigma_z")
    direct = value.get("direct_lower_collar")
    normalized = value.get("normalized_inner_cutoff")
    _require(isinstance(direct, Mapping), "direct consumer record missing")
    _require(isinstance(normalized, Mapping), "normalized consumer record missing")
    _require(direct.get("scope") == "direct lower-collar chi sigma(s) only", "direct consumer scope")
    _require(direct.get("interval_leaf_call_count") == 1, "direct consumer interval call count")
    _require(direct.get("fixed_s_in_z") is True, "direct fixed-s condition")
    _require(direct.get("sigma_Z") == Fraction(0), "direct conditional sigma_Z")
    _require(direct.get("sigma_z") == Fraction(0), "direct conditional sigma_z")
    _require(normalized.get("scope") == "first-chart normalized inner cutoff only", "normalized consumer scope")
    _require(normalized.get("normalized_inner_cutoff_claimed") is True, "normalized consumer record")
    _require(normalized.get("normalized_inner_cutoff_live_binding_claimed") is False, "normalized live binding claim")
    _require(normalized.get("normalized_inner_cutoff_node_binding") is False, "normalized node binding claim")
    _require(value.get("expected_interval_leaf_call_count") == 2, "consumer expected interval calls")
    _require(value.get("interval_leaf_call_count") == 2, "consumer interval call count")
    _require(len(published_calls) == 2, "published sigma evaluator call count")
    direct_expected = (Fraction(0), Fraction(1, 2**128))
    normalized_expected = (Fraction(0), Fraction(1))
    _require(published_calls[0] == direct_expected, "direct published sigma domain")
    _require(published_calls[1] == normalized_expected, "normalized published sigma domain")
    _require(direct.get("interval_leaf_input") == direct_expected, "direct consumer interval input")
    _require(
        normalized.get("active_window", {}).get("interval_leaf_input")
        == normalized_expected,
        "normalized consumer interval input",
    )


def _consume(hashes: Mapping[str, str], progress: dict[str, Any]) -> dict[str, Any]:
    # Authenticate all frozen bytes first.  The fixture and JSON decode are
    # deliberately ordered after this check; no gzip API appears in this file.
    raw = _authenticate_frozen_bytes(hashes)
    progress["frozen_bytes_authenticated"] = True

    leaf = _load_component(LEAF_REL, COMPONENT + "_leaf")
    fixture = _load_component(FIXTURE_REL, COMPONENT + "_fixture")
    sigma_leaf = _load_component(SIGMA_REL, COMPONENT + "_published_sigma_leaf")
    progress["imports_after_entry"] = True
    _require(getattr(leaf, "SCHEMA", None) == LEAF_SCHEMA, "new leaf schema")
    _require(getattr(fixture, "SCHEMA", None) == FIXTURE_SCHEMA, "new fixture schema")
    _require(getattr(sigma_leaf, "SCHEMA", None) == SIGMA_SCHEMA, "published sigma schema")

    fixture_runner = getattr(fixture, "run_fixture", None)
    _require(callable(fixture_runner), "fixture run_fixture API")
    fixture_result = fixture_runner(leaf)
    progress["fixture_calls"] = 1
    progress["fixture_result"] = _json_safe(fixture_result)
    _validate_fixture_result(fixture_result)

    # The saved source remains undecoded until the synthetic fixture has
    # completed.  It is an ordinary JSON object, never a gzip archive.
    saved_mapping = _decode_object(raw[SAVED_JSON_REL], "saved ordinary JSON source")
    progress["saved_json_decoded_after_fixture"] = True

    evaluator = getattr(sigma_leaf, "evaluate_sigma_interval", None)
    _require(callable(evaluator), "published sigma evaluator API")
    published_calls: list[tuple[Fraction, Fraction]] = []

    def injected_evaluator(lower: Fraction, upper: Fraction) -> Any:
        published_calls.append((lower, upper))
        return evaluator(lower, upper)

    actual_runner = getattr(leaf, "evaluate_actual_collar_sigma", None)
    _require(callable(actual_runner), "actual consumer API")
    actual_result = actual_runner(
        saved_mapping,
        injected_evaluator,
        fixed_s_in_z=True,
        fixed_s_over_sc_in_z=True,
    )
    progress["actual_consumer_calls"] = 1
    progress["published_sigma_evaluator_calls"] = len(published_calls)
    progress["actual_result"] = _json_safe(actual_result)
    _validate_actual_result(actual_result, published_calls)

    selected = actual_result.get("selected_sc_tuple")
    direct = actual_result["direct_lower_collar"]
    normalized = actual_result["normalized_inner_cutoff"]
    return {
        "fixture": _json_safe(fixture_result),
        "actual_result": {
            "scope": "direct_lower_collar_chi_sigma_s_only",
            "direct_lower_collar": _json_safe(direct),
            "normalized_inner_cutoff": {
                "result": _json_safe(normalized),
                "excluded_from_runner_scope": True,
            },
            "consumer_result": _json_safe(actual_result),
        },
        "selected_sc_tuple": _json_safe(selected),
        "saved_json_source_sha256": hashes[SAVED_JSON_REL],
        "normalized_cutoff_source_sha256": hashes[NORMALIZED_SOURCE_REL],
        "normalized_cutoff_source_pinned_only": True,
        "normalized_cutoff_source_used": False,
        "direct_lower_collar_source_sha256": hashes[DIRECT_SOURCE_REL],
        "direct_lower_collar_source_pinned_only": True,
        "direct_lower_collar_source_used": False,
        "normalized_inner_cutoff_computed_inside_consumer": True,
        "normalized_inner_cutoff_excluded_from_runner_status": True,
        "direct_lower_collar_sigma_only": True,
        "ancestor_reexecuted": False,
        "old_fixture_called": False,
        "physical_m15": False,
        "source_identity_verified": False,
        "live_owner_binding": False,
        "published_sigma_evaluator_calls": len(published_calls),
        "published_sigma_evaluator_call_labels": [
            "direct_lower_collar_sigma_s",
            "normalized_first_chart_active_window_excluded",
        ],
        "actual_consumer_calls": 1,
        "fixed_s_in_z_condition": True,
        "fixed_s_over_sc_in_z_condition": True,
        "fixed_s_conditions_are_conditional_not_live_identity": True,
    }


def main() -> int:
    try:
        hashes, source_sha256, brief_sha256, manifest = _preflight()
    except BaseException:
        print("status=preflight_failed")
        return 1

    try:
        entry_sha256 = _write_exclusive(
            ENTRY,
            {
                "schema": COMPONENT + "_entry",
                "component": COMPONENT,
                "phase": PHASE,
                "manifest_schema": manifest.get("schema"),
                "manifest_sha256": source_sha256,
                "brief_sha256": brief_sha256,
                "input_hashes": hashes,
                "gzip_forbidden": True,
                "ordinary_json_source_only": True,
                "fixed_s_in_z_condition": True,
                "fixed_s_over_sc_in_z_condition": True,
                "fixed_s_conditions_are_conditional_not_live_identity": True,
                "entry_written_before_import": True,
                "entered": True,
            },
        )
    except BaseException:
        print("status=entry_failed")
        return 1

    progress: dict[str, Any] = {
        "frozen_bytes_authenticated": False,
        "imports_after_entry": False,
        "fixture_calls": 0,
        "saved_json_decoded_after_fixture": False,
        "actual_consumer_calls": 0,
        "published_sigma_evaluator_calls": 0,
    }
    outcome: dict[str, Any] = {
        "schema": COMPONENT + "_result",
        "component": COMPONENT,
        "phase": PHASE,
        "manifest_sha256": source_sha256,
        "brief_sha256": brief_sha256,
        "entry_sha256": entry_sha256,
        "input_hashes": hashes,
        "status": "failed",
        "progress": progress,
        "fixture": None,
        "actual_result": None,
        "error": None,
        "gzip_forbidden": True,
        "ordinary_json_source_only": True,
        "direct_lower_collar_sigma_only": True,
        "normalized_inner_cutoff_claimed": False,
        "fixed_s_in_z_condition": True,
        "fixed_s_over_sc_in_z_condition": True,
        "fixed_s_conditions_are_conditional_not_live_identity": True,
        "physical_m15": False,
        "source_identity_verified": False,
        "live_owner_binding": False,
        "ancestor_reexecuted": False,
        "old_fixture_called": False,
    }
    try:
        consumed = _consume(hashes, progress)
        outcome.update(consumed)
        outcome["status"] = "returned"
    except BaseException as exc:
        outcome["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }

    try:
        _write_exclusive(RESULT, outcome)
    except BaseException:
        # The entry exists, so a result write failure is itself a failed state;
        # keep stdout status-only as required by the one-entry protocol.
        print("status=result_write_failed")
        return 1
    print("status=" + str(outcome["status"]))
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
