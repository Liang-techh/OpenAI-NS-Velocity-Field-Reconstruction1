"""One-entry conditional Rc source-recipe proof; no archived graph access."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
COMPONENT = "rp_rc_recipe_proof_v1"
DIR = ROOT / "tmp" / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / name for name in ("source.json", "entered.json", "result.json"))
PREFIX = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_"
LEAF = "src/openai_ns_reconstruction/rp_rc_recipe_proof.py"
RUNNER = "scripts/run_rp_rc_recipe_proof_once.py"
F_RECEIPT = PREFIX + "Rp_native_C5_pulse_source_check.json"
OLD_RECEIPT = PREFIX + "C3_Rp_native_frame_check.json"
GEOMETRY = PREFIX + "whole_Z_all_N_control_connection.py"
FRAME = PREFIX + "C3_Rp_native_frame.py"
PRESCRIBED = PREFIX + "prescribed_R100_Rm_C6_functions.py"
OUTER = PREFIX + "true_outer_C6_C5_functions.py"
FIXED = {
    F_RECEIPT: "320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15",
    OLD_RECEIPT: "d34d0e4f90b131e68bbde16caad33029ff89b75174741ad02d03d0a6e135d9e8",
    GEOMETRY: "aed09a89972ef9de0234e72691c4b5d33dfaa06d82d07b8df0781ef0f7a3807d",
    FRAME: "cb7bc1bc1801e04c198f40458dd88246434e74318f500260e8a9043c49398674",
    PRESCRIBED: "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    OUTER: "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3",
}
REQUIRED = frozenset(FIXED) | {LEAF, RUNNER}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def exclusive(path, payload):
    raw = (json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n").encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write(raw)
    return sha(raw)


def preflight():
    if ENTRY.exists() or RESULT.exists(): raise ValueError("component already entered")
    brief_raw = BRIEF.read_bytes()
    brief = json.loads(brief_raw)
    if (brief.get("component") != COMPONENT or brief.get("phase") != COMPONENT + "_pending"
            or brief.get("goal_active") is not True or brief.get("full_goal_complete") is not False
            or brief.get("do_not_reenter") is not True or brief.get("one_entry_only") is not True):
        raise ValueError("pending brief guard failed")
    cp = brief["source_checkpoint"]
    p = Path(cp["path"])
    if not p.is_absolute() or p.resolve() != SOURCE.resolve(): raise ValueError("checkpoint path mismatch")
    src = SOURCE.read_bytes()
    if sha(src) != cp["sha256"]: raise ValueError("checkpoint hash mismatch")
    manifest = json.loads(src)
    if (manifest.get("component") != COMPONENT or manifest.get("schema") != COMPONENT + "_checkpoint_v1"
            or manifest.get("one_entry_only") is not True): raise ValueError("manifest guard failed")
    entries = manifest.get("files")
    if not isinstance(entries, list) or len(entries) != len(REQUIRED): raise ValueError("frozen input count")
    hashes = {}
    for entry in entries:
        path, expected = entry["path"], entry["sha256"]
        if path not in REQUIRED or path in hashes: raise ValueError("unexpected/duplicate input")
        actual = sha((ROOT / path).read_bytes())
        if actual != expected or (path in FIXED and actual != FIXED[path]): raise ValueError("input hash " + path)
        hashes[path] = actual
    if set(hashes) != REQUIRED: raise ValueError("missing input")
    return hashes, sha(src), sha(brief_raw)


def pinned(path, hashes):
    raw = (ROOT / path).read_bytes()
    if sha(raw) != hashes[path]: raise ValueError("input changed after entry: " + path)
    return raw


def consume(hashes):
    f = json.loads(pinned(F_RECEIPT, hashes)); old = json.loads(pinned(OLD_RECEIPT, hashes))
    if (f.get("all_passed") is not True or old.get("all_passed") is not True
            or not isinstance(f.get("source_family"), dict) or not f["source_family"]
            or old.get("source_family") != f["source_family"]):
        raise ValueError("receipt acceptance or source family mismatch")
    pins = f.get("input_hashes", {})
    for path in (OLD_RECEIPT, GEOMETRY, FRAME, PRESCRIBED, OUTER):
        if pins.get(Path(path).name) != hashes[path]: raise ValueError("F receipt pin mismatch: " + path)
    if old.get("current_original_Rp_native_exact_radius_identity_installed") is not True:
        raise ValueError("old radius receipt gate missing")
    p = old.get("independent_actual_native_radius", {})
    if p.get("original_Rc_is_actual_O3_power_right_offset_and_R2_is_twice_Rc") is not True:
        raise ValueError("old O3 endpoint provenance gate missing")
    raw_leaf = pinned(LEAF, hashes)
    spec = importlib.util.spec_from_file_location("rp_rc_recipe_proof_v1_new_leaf", ROOT / LEAF)
    if spec is None or spec.loader is None: raise ValueError("new leaf import failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    proof = module.prove(*(pinned(path, hashes).decode("utf-8") for path in (GEOMETRY, FRAME, PRESCRIBED, OUTER)))
    proof["two_receipt_source_families_match"] = True
    proof["old_C3_radius_gate_accepted"] = True
    proof["old_receipt_declares_O3_right_offset"] = True
    proof["old_phase_rc_identity_to_extracted_offset_proved"] = False
    proof["proof_leaf_verified_after_entry"] = sha(raw_leaf) == hashes[LEAF]
    proof["archived_graph_or_gzip_decoded"] = False
    proof["protected_refs_used_as_graph_or_mapping_keys"] = False
    return proof


def main():
    hashes, source_sha, brief_sha = preflight()
    entry_sha = exclusive(ENTRY, {"component": COMPONENT, "entered": True,
        "source_checkpoint_sha256": source_sha, "brief_sha256": brief_sha, "input_hashes": hashes, "consumer_calls": 0})
    output = {"component": COMPONENT, "schema": COMPONENT + "_result_v1", "entry_sha256": entry_sha,
              "source_checkpoint_sha256": source_sha, "consumer_calls": 0}
    try:
        output["consumer_calls"] = 1
        output["proof"] = consume(hashes)
        output["status"] = "returned"
        output["error"] = None
    except BaseException as exc:
        output["status"] = "failed"
        output["error"] = {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()}
    result_sha = exclusive(RESULT, output)
    print(json.dumps({"component": COMPONENT, "status": output["status"], "entry_sha256": entry_sha,
                      "source_sha256": source_sha, "result_sha256": result_sha}))
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
