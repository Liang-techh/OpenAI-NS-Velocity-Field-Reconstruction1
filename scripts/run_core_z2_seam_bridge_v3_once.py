"""Enter the saved core Z2 seam bridge exactly once.

All bytes are authenticated before either compact receipt is decompressed.
The independent synthetic fixture runs before saved receipt consumption.  The
new binder receives the saved exact-MPF tuple dictionaries directly as
``build(old, new_profile, new_seam, low_output)``; this runner never materializes an
interval with mpmath and never imports an ancestor implementation.
"""

from __future__ import annotations

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
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
NAME = "core_z2_seam_bridge_v3"
DIR = ROOT / "tmp" / NAME
SOURCE, ENTRY, RESULT = (DIR / name for name in ("source.json", "entered.json", "result.json"))

LEAF = "src/openai_ns_reconstruction/core_z2_seam_bridge_v3.py"
FIXTURE = "src/openai_ns_reconstruction/core_z2_seam_bridge_fixture_v3.py"
RUNNER = "scripts/run_core_z2_seam_bridge_v3_once.py"
OLD = "experiments/root_st073/core_secondword_partial_v2_receipt.json.gz"
NEW = "experiments/root_st073/core_firstword_z2_profile_v2_receipt.json.gz"
LOW = "experiments/root_st073/n1_core_low_word_v3_receipt.json.gz"
INCIDENT = "experiments/root_st073/core_z2_seam_bridge_v1_uncheckpointed_incident.json"
SUPPLEMENT = "experiments/root_st073/core_z2_seam_bridge_v1_uncheckpointed_supplement.md"
FAILURE = "experiments/root_st073/core_z2_seam_bridge_v2_failure.json.gz"
FILES = (LEAF, FIXTURE, RUNNER, OLD, NEW, LOW, INCIDENT, SUPPLEMENT, FAILURE)
REQUIRED = FILES
RECEIPTS = {OLD, NEW, LOW, FAILURE}

PIN = {
    FAILURE: "9194e46466392eb45cc8d4903fa73993b2d6865bb96c36c693ed5cec2bdc5c79",
    OLD: "d7de2dbcd1fbc90bf039d0b8527cc82be6e7826e475c2c0c1fdfa6ed3c8585fe",
    NEW: "bcbd81914ba8294a8da13790640d63ab12907fd41624fceff10f304a5ca12115",
    LOW: "5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6",
}
LOW_SHA256 = "5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6"
Z2 = "experiments/root_st073/n1_core_word1_z2_v2_receipt.json.gz"
Z2_SHA256 = "c53ff32ad3ce5288279953e84bc7417b71477c82b00d94466d0d2880adbc7ac9"

OLD_SCHEMA = "core_secondword_partial_v2"
NEW_PROFILE_SCHEMA = "core_firstword_z2_profile_v2"
NEW_SEAM_SCHEMA = NEW_PROFILE_SCHEMA + ":missing_seam_v1"
COMMON_LOG_KEYS = ("F0_anchor", "Lambda", "Pstar_squared", "hb")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def require(condition: object, message: str) -> None:
    if not condition:
        raise ValueError(message)


def encode(value: object) -> object:
    """Convert mappings/tuples to JSON values without numeric materialization."""

    if isinstance(value, Mapping):
        return {str(key): encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError("unsupported output value " + type(value).__name__)


def save(path: Path, value: object) -> str:
    raw = (
        json.dumps(
            encode(value),
            sort_keys=True,
            ensure_ascii=False,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    with path.open("xb") as handle:
        handle.write(raw)
    return sha(raw)


def module(path: str, name: str) -> Any:
    """Load only the new leaf or its independent fixture."""

    require(path in (LEAF, FIXTURE), "unexpected module import " + path)
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    require(spec is not None and spec.loader is not None, "module spec " + path)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def manifest_hashes(manifest: object, label: str) -> dict[str, str]:
    require(isinstance(manifest, Mapping), label + " source manifest")
    rows = manifest.get("files")
    require(isinstance(rows, list), label + " manifest files")
    result: dict[str, str] = {}
    for row in rows:
        require(isinstance(row, Mapping), label + " manifest row")
        path = row.get("path")
        digest = row.get("sha256")
        require(isinstance(path, str) and isinstance(digest, str), label + " manifest path/hash")
        require(path not in result, label + " duplicate manifest path " + path)
        result[path] = digest
    return result


def preflight() -> tuple[dict[str, str], str, str]:
    require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    require(SOURCE.exists(), "source manifest missing")

    brief_raw = BRIEF.read_bytes()
    brief = json.loads(brief_raw)
    source_raw = SOURCE.read_bytes()
    source = json.loads(source_raw)
    checkpoint = brief.get("source_checkpoint")

    require(brief.get("component") == NAME, "brief component guard")
    require(brief.get("phase") == NAME + "_pending", "brief phase guard")
    require(brief.get("goal_active") is True and brief.get("full_goal_complete") is False, "goal guard")
    require(brief.get("one_entry_only") is True and brief.get("do_not_reenter") is True, "one-entry guard")
    require(isinstance(checkpoint, Mapping), "source checkpoint missing")
    require(Path(checkpoint["path"]).resolve() == SOURCE.resolve(), "source checkpoint path mismatch")
    require(checkpoint["sha256"] == sha(source_raw), "source checkpoint hash mismatch")
    preimage = Path(checkpoint["preimage_path"])
    require(sha(preimage.read_bytes()) == checkpoint["preimage_sha256"], "brief preimage hash mismatch")
    require(source.get("component") == NAME and source.get("one_entry_only") is True, "manifest guard")

    hashes: dict[str, str] = {}
    rows = source.get("files")
    require(isinstance(rows, list), "manifest files missing")
    for row in rows:
        require(isinstance(row, Mapping), "manifest row must be a mapping")
        path = row.get("path")
        require(isinstance(path, str), "manifest path missing")
        require(path in REQUIRED and path not in hashes, "unexpected or duplicate source file: " + path)
        digest = sha((ROOT / path).read_bytes())
        require(digest == row.get("sha256"), "manifest hash mismatch: " + path)
        if path in PIN:
            require(digest == PIN[path], "pinned input changed: " + path)
        if path in RECEIPTS:
            require(row.get("raw_sha256_before_decode") is True, "receipt is not raw-pinned: " + path)
        hashes[path] = digest
    require(set(hashes) == set(REQUIRED), "manifest input set mismatch")
    return hashes, sha(source_raw), sha(brief_raw)


def receipt(raw: Mapping[str, bytes], path: str) -> dict[str, Any]:
    """Decode a compact receipt only after its raw bytes passed preflight."""

    value = json.loads(gzip.decompress(raw[path]))
    require(isinstance(value, dict), "receipt object " + path)
    result = value.get("result")
    require(isinstance(result, Mapping) and result.get("status") == "returned", "input receipt " + path)
    return value


def exact_tuple_record(value: object, label: str) -> dict[str, list[int]]:
    """Validate one saved exact MPF endpoint pair without creating an interval."""

    require(isinstance(value, Mapping), label + " must be a mapping")
    lower = value.get("lower_exact_mpf_tuple")
    upper = value.get("upper_exact_mpf_tuple")
    require(isinstance(lower, (list, tuple)) and isinstance(upper, (list, tuple)), label + " endpoints missing")
    require(len(lower) == len(upper) == 4, label + " endpoints must be four-tuples")
    require(
        all(type(item) is int for item in tuple(lower) + tuple(upper)),
        label + " endpoints must contain integers",
    )
    return {"lower_exact_mpf_tuple": list(lower), "upper_exact_mpf_tuple": list(upper)}


def tuple_paths(value: object, path: tuple[object, ...] = ()) -> list[list[object]]:
    """Return paths of all exact MPF tuple records in a saved packet."""

    found: list[list[object]] = []
    if isinstance(value, Mapping):
        has_lower = "lower_exact_mpf_tuple" in value
        has_upper = "upper_exact_mpf_tuple" in value
        if has_lower or has_upper:
            exact_tuple_record(value, ".".join(map(str, path)) or "packet")
            found.append(list(path))
        for key, item in value.items():
            found.extend(tuple_paths(item, path + (key,)))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            found.extend(tuple_paths(item, path + (index,)))
    return found


def tuple_digest(record: dict[str, list[int]]) -> str:
    raw = json.dumps(record, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return sha(raw)


def compare_source_log_tuples(
    old_output: Mapping[str, Any],
    new_output: Mapping[str, Any],
    low_output: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    """Compare only saved endpoint tuples at corresponding source-log paths."""

    old_logs = old_output.get("source_logs")
    new_logs = new_output.get("source_logs")
    low_logs = low_output.get("source_logs")
    require(
        isinstance(old_logs, Mapping)
        and isinstance(new_logs, Mapping)
        and isinstance(low_logs, Mapping),
        "source logs missing",
    )
    comparisons: dict[str, dict[str, Any]] = {}
    for key in COMMON_LOG_KEYS:
        old_record = exact_tuple_record(old_logs.get(key), "old source_logs." + key)
        new_record = exact_tuple_record(new_logs.get(key), "new source_logs." + key)
        low_record = exact_tuple_record(low_logs.get(key), "low source_logs." + key)
        equal = old_record == new_record == low_record
        require(equal, "source-log exact tuple mismatch: " + key)
        comparisons[key] = {
            "old_path": ["result", "output", "source_logs", key],
            "new_path": ["result", "output", "source_logs", key],
            "low_path": ["result", "output", "source_logs", key],
            "old_tuple_sha256": tuple_digest(old_record),
            "new_tuple_sha256": tuple_digest(new_record),
            "low_tuple_sha256": tuple_digest(low_record),
            "equal": True,
        }
    return comparisons


def lineage_evidence(
    old_doc: Mapping[str, Any],
    new_doc: Mapping[str, Any],
    low_doc: Mapping[str, Any],
    old_hash: str,
    new_hash: str,
    low_hash: str,
) -> dict[str, Any]:
    """Validate the saved Z2 -> LOW relation without declaring an owner."""

    old_manifest = manifest_hashes(old_doc.get("source"), "old receipt")
    new_manifest = manifest_hashes(new_doc.get("source"), "new receipt")
    low_manifest = manifest_hashes(low_doc.get("source"), "LOW receipt")
    require(old_manifest.get(LOW) == LOW_SHA256, "old source manifest does not pin LOW")
    require(new_manifest.get(LOW) == LOW_SHA256, "new source manifest does not pin LOW")
    require(new_manifest.get(Z2) == Z2_SHA256, "new source manifest does not pin Z2")
    require(low_manifest.get(Z2) == Z2_SHA256, "LOW source manifest does not pin Z2")

    old_result = old_doc.get("result")
    new_result = new_doc.get("result")
    low_result = low_doc.get("result")
    require(
        isinstance(old_result, Mapping)
        and isinstance(new_result, Mapping)
        and isinstance(low_result, Mapping),
        "receipt result envelope",
    )
    old_output = old_result.get("output")
    new_output = new_result.get("output")
    low_output = low_result.get("output")
    require(
        isinstance(old_output, Mapping)
        and isinstance(new_output, Mapping)
        and isinstance(low_output, Mapping),
        "receipt output envelope",
    )

    old_saved = old_output.get("input_hashes")
    if isinstance(old_saved, Mapping) and LOW in old_saved:
        require(old_saved[LOW] == LOW_SHA256, "old saved validation does not pin LOW")
    new_saved = new_output.get("input_receipt_sha256")
    require(isinstance(new_saved, Mapping), "new saved receipt validation missing")
    require(new_saved.get(OLD) == old_hash, "new saved validation does not pin old receipt")
    require(new_saved.get(LOW) == LOW_SHA256, "new saved validation does not pin LOW")
    require(new_saved.get(Z2) == Z2_SHA256, "new saved validation does not pin Z2")

    low_saved = low_output.get("input_receipts")
    require(isinstance(low_saved, Mapping), "LOW saved receipt validation missing")
    require(low_saved.get(Z2) == Z2_SHA256, "LOW saved validation does not pin Z2")

    anchor = new_output.get("anchor_receipt")
    require(isinstance(anchor, Mapping), "new saved anchor missing")
    require(anchor.get("path") == OLD and anchor.get("sha256") == old_hash, "new saved anchor mismatch")
    anchor_manifest = manifest_hashes(anchor.get("source_manifest"), "saved anchor")
    require(anchor_manifest.get(LOW) == LOW_SHA256, "saved anchor does not pin LOW")

    return {
        "old_manifest_low": {"path": LOW, "sha256": old_manifest[LOW]},
        "new_manifest_low": {"path": LOW, "sha256": new_manifest[LOW]},
        "new_manifest_z2": {"path": Z2, "sha256": new_manifest[Z2]},
        "low_manifest_z2": {"path": Z2, "sha256": low_manifest[Z2]},
        "saved_new_input_receipt_sha256": {
            OLD: new_saved[OLD],
            LOW: new_saved[LOW],
            Z2: new_saved[Z2],
        },
        "saved_anchor": {"path": anchor["path"], "sha256": anchor["sha256"]},
        "saved_low_input_receipts": {Z2: low_saved[Z2]},
        "low_receipt_sha256": low_hash,
        "physical_owner_proved": False,
    }


def consume(hashes: dict[str, str], progress: dict[str, Any]) -> dict[str, Any]:
    # No receipt is decoded until every frozen file, including all three gzip
    # inputs, has been authenticated against the source manifest.
    raw: dict[str, bytes] = {}
    for path in REQUIRED:
        data = (ROOT / path).read_bytes()
        require(sha(data) == hashes[path], "input changed after entry: " + path)
        raw[path] = data

    leaf = module(LEAF, NAME + "_leaf")
    fixture = module(FIXTURE, NAME + "_fixture")

    # The independent synthetic fixture is deliberately first.  Its complete
    # return or traceback remains in the saved result even when it fails.
    try:
        fixture_runner = getattr(fixture, "run_fixture", None)
        require(callable(fixture_runner), "fixture must export run_fixture")
        fixture_result = fixture_runner(leaf)
        progress["fixture"] = encode(fixture_result)
        require(isinstance(fixture_result, Mapping) and fixture_result.get("passed") is True, "synthetic fixture failed")
    except BaseException as exc:
        progress["fixture_error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
        raise

    old_doc = receipt(raw, OLD)
    new_doc = receipt(raw, NEW)
    low_doc = receipt(raw, LOW)
    old_result = old_doc["result"]
    new_result = new_doc["result"]
    low_result = low_doc["result"]
    old_output = old_result["output"]
    new_output = new_result["output"]
    low_output = low_result["output"]
    require(
        isinstance(old_output, Mapping)
        and isinstance(new_output, Mapping)
        and isinstance(low_output, Mapping),
        "receipt outputs must be mappings",
    )

    old_packet = old_output.get("result")
    new_profile = new_output.get("result")
    new_seam = new_output.get("missing_seam_result")
    require(isinstance(old_packet, Mapping), "old result.output.result missing")
    require(isinstance(new_profile, Mapping), "new result.output.result missing")
    require(isinstance(new_seam, Mapping), "new result.output.missing_seam_result missing")
    require(old_packet.get("schema") == OLD_SCHEMA, "old packet schema mismatch")
    require(new_profile.get("schema") == NEW_PROFILE_SCHEMA, "new profile schema mismatch")
    require(new_seam.get("schema") == NEW_SEAM_SCHEMA, "new seam schema mismatch")

    lineage = lineage_evidence(old_doc, new_doc, low_doc, hashes[OLD], hashes[NEW], hashes[LOW])
    exact_log_comparisons = compare_source_log_tuples(old_output, new_output, low_output)

    old_family = old_output.get("source_family")
    new_family = new_output.get("source_family")
    low_family = low_output.get("source_family")
    if old_family is not None or new_family is not None or low_family is not None:
        require(old_family == new_family == low_family, "saved source-family mismatch")

    progress["input_evidence"] = {
        "old_packet_path": [OLD, "result", "output", "result"],
        "new_profile_packet_path": [NEW, "result", "output", "result"],
        "new_seam_packet_path": [NEW, "result", "output", "missing_seam_result"],
        "low_output_path": [LOW, "result", "output"],
        "input_receipt_sha256": {OLD: hashes[OLD], NEW: hashes[NEW], LOW: hashes[LOW]},
        "lineage": lineage,
        "exact_source_log_tuple_comparisons": exact_log_comparisons,
        "physical_owner_proved": False,
    }
    progress["actual_build_calls"] = 1
    progress["ancestor_called"] = False

    # This is the sole actual invocation of the new binder.  The four inputs
    # are the authenticated JSON dictionaries; no mpmath conversion or
    # ancestor builder is involved.
    result = leaf.build(old_packet, new_profile, new_seam, low_output)
    require(isinstance(result, Mapping), "binder result must be a mapping")
    require(result.get("physical_owner_proved") is not True, "binder made a physical-owner claim")

    return {
        "source_family": old_family,
        "input_evidence": progress["input_evidence"],
        "fixture": progress["fixture"],
        "result": encode(result),
        "actual_build_calls": 1,
        "ancestor_called": False,
        "physical_owner_proved": False,
        "old_records_merged": result.get("old_records_merged") if isinstance(result, Mapping) else None,
    }


def main() -> int:
    hashes, source_sha256, brief_sha256 = preflight()
    entry_sha256 = save(
        ENTRY,
        {
            "component": NAME,
            "source_sha256": source_sha256,
            "brief_sha256": brief_sha256,
            "input_hashes": hashes,
            "raw_gzip_authenticated_before_decode": True,
            "entered": True,
            "ancestor_imported": False,
        },
    )
    out: dict[str, Any] = {
        "component": NAME,
        "source_sha256": source_sha256,
        "entry_sha256": entry_sha256,
        "status": "failed",
        "raw_gzip_authenticated_before_decode": True,
        "ancestor_called": False,
        "physical_owner_proved": False,
        "actual_build_calls": 0,
    }
    try:
        out.update(status="returned", output=consume(hashes, out), error=None)
    except BaseException as exc:
        out["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
    result_sha256 = save(RESULT, out)
    print(
        json.dumps(
            {
                "component": NAME,
                "status": out["status"],
                "result_sha256": result_sha256,
                "error": out.get("error", {}).get("message"),
            }
        )
    )
    return 0 if out["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
