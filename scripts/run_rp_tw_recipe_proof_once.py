"""One-entry, raw-authenticated AST proof for old/current ``Tw`` recipes.

The runner reads the pending brief and its source checkpoint, authenticates
every frozen input before decoding an archive, and imports the new leaf only
after exclusive creation of ``entered.json``.  It never imports or calls an
old producer and never selects a graph node.
"""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import traceback


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
COMPONENT = "rp_tw_recipe_proof_v1"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE, ENTRY, RESULT = (DIRECTORY / name for name in ("source.json", "entered.json", "result.json"))
PREFIX = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_"
F_ARCHIVE = PREFIX + "Rp_native_C5_pulse_source.json.gz"
F_RECEIPT = PREFIX + "Rp_native_C5_pulse_source_check.json"
CONTROL_ARCHIVE = PREFIX + "whole_Z_all_N_control_connection.json.gz"
C3_ARCHIVE = PREFIX + "C3_limit_controls.json.gz"
CONTROL_SOURCE = PREFIX + "whole_Z_all_N_control_connection.py"
OLD_SOURCE = PREFIX + "C3_Rp_native_frame.py"
OUTER_SOURCE = PREFIX + "true_outer_C6_C5_functions.py"
LEAF = "src/openai_ns_reconstruction/rp_tw_recipe_proof.py"
RUNNER = "scripts/run_rp_tw_recipe_proof_once.py"

FIXED = {
    F_ARCHIVE: "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe",
    F_RECEIPT: "320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15",
    CONTROL_ARCHIVE: "588f411e9122b9a76e6a42efc2b5754137f6af66a9cae42f8d6a7c2c59920de8",
    C3_ARCHIVE: "23e4f08eb1c7b2f0098953c50d28df71039905cba7ec6bfc6f22e1aa356afe34",
    CONTROL_SOURCE: "aed09a89972ef9de0234e72691c4b5d33dfaa06d82d07b8df0781ef0f7a3807d",
    OLD_SOURCE: "cb7bc1bc1801e04c198f40458dd88246434e74318f500260e8a9043c49398674",
    OUTER_SOURCE: "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3",
}
REQUIRED = frozenset(FIXED) | {LEAF, RUNNER}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def exclusive(path: Path, payload: dict) -> str:
    raw = (json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n").encode()
    with path.open("xb") as stream:
        stream.write(raw)
    return sha(raw)


def preflight() -> tuple[dict[str, str], str, str]:
    if ENTRY.exists() or RESULT.exists():
        raise ValueError("component already entered")
    brief_raw = BRIEF.read_bytes()
    brief = json.loads(brief_raw)
    if (brief.get("component") != COMPONENT or brief.get("phase") != COMPONENT + "_pending"
            or brief.get("goal_active") is not True or brief.get("full_goal_complete") is not False
            or brief.get("do_not_reenter") is not True or brief.get("one_entry_only") is not True):
        raise ValueError("pending brief guard failed")
    checkpoint = brief["source_checkpoint"]
    checkpoint_path = Path(checkpoint["path"])
    if not checkpoint_path.is_absolute() or checkpoint_path.resolve() != SOURCE.resolve():
        raise ValueError("source checkpoint path mismatch")
    source_raw = SOURCE.read_bytes()
    if sha(source_raw) != checkpoint["sha256"]:
        raise ValueError("source checkpoint hash mismatch")
    source = json.loads(source_raw)
    if (source.get("component") != COMPONENT or source.get("schema") != COMPONENT + "_checkpoint_v1"
            or source.get("one_entry_only") is not True):
        raise ValueError("source checkpoint scope mismatch")
    entries = source.get("files")
    if not isinstance(entries, list) or len(entries) != len(REQUIRED):
        raise ValueError("exact source input set required")
    hashes: dict[str, str] = {}
    for entry in entries:
        path, expected = entry["path"], entry["sha256"]
        if path not in REQUIRED or path in hashes:
            raise ValueError("unexpected or duplicate frozen input")
        actual = sha((ROOT / path).read_bytes())
        if actual != expected or (path in FIXED and actual != FIXED[path]):
            raise ValueError("frozen raw hash mismatch: " + path)
        hashes[path] = actual
    if set(hashes) != REQUIRED:
        raise ValueError("missing frozen input")
    return hashes, sha(source_raw), sha(brief_raw)


def parse_pinned(path: str, hashes: dict[str, str]) -> dict:
    raw = (ROOT / path).read_bytes()
    if sha(raw) != hashes[path]:
        raise ValueError("source changed after entry: " + path)
    if path.endswith(".gz"):
        raw = gzip.decompress(raw)
    data = json.loads(raw)
    if not isinstance(data, dict):
        raise ValueError("source object required: " + path)
    return data


def verify_selected_f_receipt_pins(receipt: dict, hashes: dict[str, str]) -> dict:
    pins = receipt.get("input_hashes")
    if not isinstance(pins, dict) or not pins:
        raise ValueError("F receipt input_hashes required")
    selected = (F_ARCHIVE, CONTROL_ARCHIVE, C3_ARCHIVE, CONTROL_SOURCE, OLD_SOURCE, OUTER_SOURCE)
    for path in selected:
        if pins.get(Path(path).name) != hashes[path]:
            raise ValueError("selected F receipt input pin mismatch: " + path)
    return {"selected_input_hashes_verified": True, "count": len(selected), "paths": list(selected)}


def _read_source(path: str, hashes: dict[str, str]) -> str:
    raw = (ROOT / path).read_bytes()
    if sha(raw) != hashes[path]:
        raise ValueError("source changed after entry: " + path)
    return raw.decode("utf-8")


def source_reader(hashes: dict[str, str]) -> dict:
    """Authenticate receipts/archives, then make one restricted leaf call."""
    import_path = ROOT / LEAF
    if sha(import_path.read_bytes()) != hashes[LEAF]:
        raise ValueError("new leaf changed after entry")

    # Authenticate the receipt before decoding; inspect only selected pins.
    f_receipt = parse_pinned(F_RECEIPT, hashes)
    if f_receipt.get("all_passed") is not True:
        raise ValueError("F receipt not accepted")
    receipt_pins = verify_selected_f_receipt_pins(f_receipt, hashes)
    f_archive = parse_pinned(F_ARCHIVE, hashes)
    control_archive = parse_pinned(CONTROL_ARCHIVE, hashes)
    c3_archive = parse_pinned(C3_ARCHIVE, hashes)
    family = f_archive.get("source_family")
    families_match = (isinstance(family, dict) and bool(family)
                      and f_receipt.get("source_family") == family
                      and control_archive.get("source_family") == family
                      and c3_archive.get("source_family") == family)
    if not families_match:
        raise ValueError("F/control/C3 source-family receipt mismatch")

    spec = importlib.util.spec_from_file_location("rp_tw_recipe_proof_v1_new_leaf", import_path)
    if spec is None or spec.loader is None:
        raise ValueError("new leaf cannot be loaded")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    proof = module.prove_tw_recipes_and_compare_refs(
        _read_source(CONTROL_SOURCE, hashes), _read_source(OLD_SOURCE, hashes),
        _read_source(OUTER_SOURCE, hashes), c3_archive, control_archive, f_archive)
    proof["source_family_matches"] = True
    proof["F_receipt_input_pins"] = receipt_pins
    proof["full_archive_json_decoded"] = True
    proof["complete_gzip_json_archives_decoded"] = [F_ARCHIVE, CONTROL_ARCHIVE, C3_ARCHIVE]
    proof["source_values_computed"] = False
    proof["exact_graph_nodes_selected_or_traversed"] = False
    proof["protected_refs_used_as_graph_or_mapping_keys"] = False
    return proof


def main() -> int:
    hashes, source_sha, brief_sha = preflight()
    entry_sha = exclusive(ENTRY, {
        "component": COMPONENT,
        "entered": True,
        "source_checkpoint_sha256": source_sha,
        "brief_sha256": brief_sha,
        "input_hashes": hashes,
        "consumer_calls": 0,
    })
    output = {
        "component": COMPONENT,
        "schema": COMPONENT + "_result_v1",
        "entry_sha256": entry_sha,
        "source_checkpoint_sha256": source_sha,
        "source_reader_calls": 0,
        "status": "failed",
    }
    try:
        output["source_reader_calls"] = 1
        output["proof"] = source_reader(hashes)
        output["status"] = "returned"
        output["error"] = None
    except BaseException as exc:
        output["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
    result_sha = exclusive(RESULT, output)
    print(json.dumps({"component": COMPONENT, "status": output["status"],
                      "entry_sha256": entry_sha, "source_sha256": source_sha,
                      "result_sha256": result_sha}))
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
