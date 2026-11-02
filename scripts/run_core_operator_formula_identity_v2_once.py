"""Enter the frozen source-AST operator formula identity proof once.

All input bytes are authenticated before the new leaf or fixture is imported.
The fixture runs first, then the leaf receives only the two authenticated
source texts through ``build(high_source_text, low_source_text)``.  The nearby
gzip receipt is a raw-byte provenance anchor and is deliberately not decoded:
this source-only leaf has no receipt input, so decoding could add no evidence.
"""

from __future__ import annotations

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
NAME = "core_operator_formula_identity_v2"
DIR = ROOT / "tmp" / NAME
SOURCE, ENTRY, RESULT = (DIR / name for name in ("source.json", "entered.json", "result.json"))

LEAF = "src/openai_ns_reconstruction/core_operator_formula_identity_v2.py"
FIXTURE = "src/openai_ns_reconstruction/core_operator_formula_identity_fixture_v2.py"
RUNNER = "scripts/run_core_operator_formula_identity_v2_once.py"
HIGH = "src/openai_ns_reconstruction/n1_core_low_source_matrix_v2.py"
LOW = "src/openai_ns_reconstruction/n1_core_low_matrix_v3.py"
RECEIPT = "experiments/root_st073/core_z2_seam_bridge_v3_receipt.json.gz"
INCIDENT = "experiments/root_st073/core_operator_scan_gzip_preflight_incident.json"
REQUIRED = (LEAF, FIXTURE, RUNNER, HIGH, LOW, RECEIPT, INCIDENT)
RECEIPTS = {RECEIPT, FAILURE}
PIN = {
    FAILURE: "3564de1416485887f1c754979368bed7f36671e2ab4bad0741b8974074c716e6",
    HIGH: "2e2a3e7e4367223f53e31e30818f53aceeab96b37c360b4e095189570ce1a111",
    LOW: "3e790126fa7ca41f3c6045004f5a98f0a4975b04681b9aa997ff7e27a69821c0",
    RECEIPT: "625fa2f92cf821301db7d02a6af92088cf9b4f229258bd4ae56bb1b35da6f2a1",
    INCIDENT: "6c99aff5729d195b33eef3ce119e13a8d233de31da4cbe501cee698706eb51f9",
}

RECEIPT_DECODE_SKIPPED_REASON = (
    "The source-AST leaf consumes only high_source_text and low_source_text; "
    "the receipt is provenance-only, so decoding it would add no input to the proof."
)


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def require(condition: object, message: str) -> None:
    if not condition:
        raise ValueError(message)


def encode(value: object) -> object:
    """Convert a leaf/fixture result to JSON without evaluating source values."""

    if isinstance(value, Mapping):
        return {str(key): encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError("unsupported result value " + type(value).__name__)


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
    """Import only the new leaf or its independent synthetic fixture."""

    require(path in (LEAF, FIXTURE), "unexpected module import " + path)
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    require(spec is not None and spec.loader is not None, "module spec " + path)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def _error(exc: BaseException) -> dict[str, str]:
    return {
        "type": type(exc).__name__,
        "message": str(exc),
        "traceback": traceback.format_exc(),
    }


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
    require(
        brief.get("goal_active") is True and brief.get("full_goal_complete") is False,
        "goal guard",
    )
    require(
        brief.get("one_entry_only") is True and brief.get("do_not_reenter") is True,
        "one-entry guard",
    )
    require(isinstance(checkpoint, Mapping), "source checkpoint missing")
    require(Path(checkpoint["path"]).resolve() == SOURCE.resolve(), "source checkpoint path mismatch")
    require(checkpoint["sha256"] == sha(source_raw), "source checkpoint hash mismatch")
    preimage = Path(checkpoint["preimage_path"])
    require(sha(preimage.read_bytes()) == checkpoint["preimage_sha256"], "brief preimage hash mismatch")

    require(source.get("component") == NAME, "manifest component guard")
    require(source.get("one_entry_only") is True, "manifest one-entry guard")
    rows = source.get("files")
    require(isinstance(rows, list), "manifest files missing")

    hashes: dict[str, str] = {}
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

    contract = source.get("input_contract")
    require(isinstance(contract, Mapping), "manifest input contract missing")
    require(contract.get("high_source_text") == HIGH, "high source contract mismatch")
    require(contract.get("low_source_text") == LOW, "low source contract mismatch")
    require(contract.get("leaf_api") == "build(high_source_text, low_source_text)", "leaf API contract mismatch")
    require(contract.get("fixture_api") == "run_fixture(leaf)", "fixture API contract mismatch")

    return hashes, sha(source_raw), sha(brief_raw)


def consume(hashes: dict[str, str], progress: dict[str, Any]) -> dict[str, Any]:
    # Re-read and authenticate every frozen byte after the exclusive entry.
    raw: dict[str, bytes] = {}
    for path in REQUIRED:
        data = (ROOT / path).read_bytes()
        require(sha(data) == hashes[path], "input changed after entry: " + path)
        raw[path] = data

    leaf = module(LEAF, NAME + "_leaf")
    fixture = module(FIXTURE, NAME + "_fixture")

    # The independent synthetic fixture is deliberately first.  Its complete
    # result remains in the saved result even when the later source build fails.
    fixture_runner = getattr(fixture, "run_fixture", None)
    require(callable(fixture_runner), "fixture must export run_fixture")
    try:
        fixture_result = fixture_runner(leaf)
        progress["fixture"] = encode(fixture_result)
        require(
            isinstance(fixture_result, Mapping) and fixture_result.get("passed") is True,
            "synthetic fixture failed",
        )
    except BaseException as exc:
        progress["fixture_error"] = _error(exc)
        raise

    high_source_text = raw[HIGH].decode("utf-8")
    low_source_text = raw[LOW].decode("utf-8")
    progress["source_texts"] = {
        "high_path": HIGH,
        "high_sha256": hashes[HIGH],
        "high_bytes": len(raw[HIGH]),
        "low_path": LOW,
        "low_sha256": hashes[LOW],
        "low_bytes": len(raw[LOW]),
    }
    progress["associated_receipt"] = {
        "path": RECEIPT,
        "raw_sha256_before_decode": hashes[RECEIPT],
        "expected_raw_sha256": PIN[RECEIPT],
        "decoded": False,
        "lineage_status_authenticated": False,
        "decode_skipped_reason": RECEIPT_DECODE_SKIPPED_REASON,
    }
    progress["preflight_incident"] = {
        "path": INCIDENT,
        "sha256": hashes[INCIDENT],
        "is_gzip": False,
        "purpose": "preserved prior gzip-prefix-before-hash incident; not decoded",
    }

    # This is the sole actual invocation of the new source-AST leaf.  It gets
    # text only; neither source module is imported and no receipt is decoded.
    progress["actual_build_calls"] = 1
    progress["ancestor_called"] = False
    result = leaf.build(high_source_text, low_source_text)
    require(isinstance(result, Mapping), "leaf.build result must be a mapping")

    return {
        "input_hashes": hashes,
        "source_text_hashes": {HIGH: hashes[HIGH], LOW: hashes[LOW]},
        "associated_receipt": progress["associated_receipt"],
        "preflight_incident": progress["preflight_incident"],
        "fixture": progress["fixture"],
        "result": encode(result),
        "actual_build_calls": 1,
        "ancestor_called": False,
        "source_modules_imported": False,
        "source_ast_only": True,
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
            "receipt_decode_skipped": True,
            "receipt_decode_skipped_reason": RECEIPT_DECODE_SKIPPED_REASON,
            "entered": True,
            "ancestor_imported": False,
        },
    )
    out: dict[str, Any] = {
        "component": NAME,
        "source_sha256": source_sha256,
        "entry_sha256": entry_sha256,
        "input_hashes": hashes,
        "status": "failed",
        "raw_gzip_authenticated_before_decode": True,
        "receipt_decode_skipped": True,
        "receipt_decode_skipped_reason": RECEIPT_DECODE_SKIPPED_REASON,
        "ancestor_called": False,
        "source_modules_imported": False,
        "actual_build_calls": 0,
        "fixture": None,
    }
    try:
        output = consume(hashes, out)
    except BaseException as exc:
        out["error"] = _error(exc)
    else:
        # Do not install ``error=None``: the previous launcher then called
        # ``.get`` on None while printing an otherwise successful result.
        out["status"] = "returned"
        out["output"] = output

    result_sha256 = save(RESULT, out)
    saved_error = out.get("error")
    error_message = saved_error.get("message") if isinstance(saved_error, Mapping) else None
    print(
        json.dumps(
            {
                "component": NAME,
                "status": out["status"],
                "result_sha256": result_sha256,
                "error": error_message,
            }
        )
    )
    return 0 if out["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
