"""One-entry source-role check and conditional cross-owner radius formula.

The two independent calls are a saved-source role reader and one synthetic
arithmetic fixture. No old producer, graph node or source-value evaluator is
called. The complete F gzip JSON is decoded only after raw authentication.
"""
from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
COMPONENT = "rp_cross_owner_signed_gap_v1"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE, ENTRY, RESULT = (DIRECTORY / n for n in ("source.json", "entered.json", "result.json"))
PREFIX = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_"
LEAF = "src/openai_ns_reconstruction/rp_cross_owner_signed_gap.py"
RUNNER = "scripts/run_rp_cross_owner_signed_gap_once.py"
F_ARCHIVE = PREFIX + "Rp_native_C5_pulse_source.json.gz"
F_RECEIPT = PREFIX + "Rp_native_C5_pulse_source_check.json"
OLD_RECEIPT = PREFIX + "C3_Rp_native_frame_check.json"
OLD_SOURCE = PREFIX + "C3_Rp_native_frame.py"
OLD_CHECKER = PREFIX + "C3_Rp_native_frame_check.py"
OUTER_SOURCE = PREFIX + "true_outer_C6_C5_functions.py"
COMPACT_SOURCE = PREFIX + "true_compact_C5_functions.py"
QUIET_SOURCE = PREFIX + "true_quiet_Rp_C5_pressure_C3_functions.py"
FIXED = {
    F_ARCHIVE: "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe",
    F_RECEIPT: "320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15",
    OLD_RECEIPT: "d34d0e4f90b131e68bbde16caad33029ff89b75174741ad02d03d0a6e135d9e8",
    OLD_SOURCE: "cb7bc1bc1801e04c198f40458dd88246434e74318f500260e8a9043c49398674",
    OLD_CHECKER: "5baaac7050500641aa1a7f1481f2fc3a2994af9d4566bcb96956f5bfd8ebd06c",
    OUTER_SOURCE: "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3",
    COMPACT_SOURCE: "f0b0333b19a0e2504d023b97e841033d3cbcc23427c57ad9672264866c01ec7b",
    QUIET_SOURCE: "0db6ea1a97184662c8e80f011844c12d143eeeb11333510f2fbbb8dcc19de068",
}
REQUIRED = frozenset(FIXED) | {LEAF, RUNNER}
NAMES = {
    "quiet": "lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions.json.gz",
    "compact": "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz",
    "outer": "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz",
}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def exclusive(path, payload):
    raw = (json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n").encode()
    with path.open("xb") as stream:
        stream.write(raw)
    return sha(raw)


def preflight():
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
    raw = SOURCE.read_bytes()
    if sha(raw) != checkpoint["sha256"]:
        raise ValueError("source checkpoint hash mismatch")
    source = json.loads(raw)
    if (source.get("component") != COMPONENT or source.get("schema") != COMPONENT + "_checkpoint_v1"
            or source.get("one_entry_only") is not True):
        raise ValueError("source checkpoint scope mismatch")
    entries = source.get("files")
    if not isinstance(entries, list) or len(entries) != len(REQUIRED):
        raise ValueError("exact source input set required")
    hashes = {}
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
    return hashes, sha(raw), sha(brief_raw)


def parse_pinned(path, hashes):
    raw = (ROOT / path).read_bytes()
    if sha(raw) != hashes[path]:
        raise ValueError("source changed after entry: " + path)
    if path.endswith(".gz"):
        raw = gzip.decompress(raw)
    data = json.loads(raw)
    if not isinstance(data, dict):
        raise ValueError("source object required")
    return data


def equal_refs(*values):
    if not all(isinstance(v, int) and not isinstance(v, bool) and v >= 0 for v in values):
        return {"status": "unresolved_noninteger_or_missing_reference"}
    return {"status": "equal_opaque_refs" if all(v == values[0] for v in values[1:]) else "different_opaque_refs",
            "refs": list(values)}


def source_reader(hashes):
    for path in (OLD_SOURCE, OLD_CHECKER, OUTER_SOURCE, COMPACT_SOURCE, QUIET_SOURCE):
        if sha((ROOT / path).read_bytes()) != hashes[path]:
            raise ValueError("source code changed after entry: " + path)
    archive = parse_pinned(F_ARCHIVE, hashes)
    receipt = parse_pinned(F_RECEIPT, hashes)
    old = parse_pinned(OLD_RECEIPT, hashes)
    family = archive.get("source_family")
    if (not isinstance(family, dict) or not family or receipt.get("source_family") != family
            or old.get("source_family") != family or receipt.get("all_passed") is not True
            or old.get("all_passed") is not True):
        raise ValueError("source family or acceptance mismatch")
    pins = receipt.get("input_hashes", {})
    for path in (F_ARCHIVE, OLD_RECEIPT, OLD_SOURCE, OLD_CHECKER, OUTER_SOURCE, COMPACT_SOURCE, QUIET_SOURCE):
        if pins.get(Path(path).name) != hashes[path]:
            raise ValueError("F receipt input pin mismatch: " + path)
    if old.get("current_original_Rp_native_exact_radius_identity_installed") is not True:
        raise ValueError("old C3 radius gate missing")
    old_proof = old.get("independent_actual_native_radius", {})
    required = ("actual_original_Rw_and_Rp_log_radius_identities",
                "exact_same_current_quiet_Rc_and_Tw_defining_expressions",
                "original_Rc_is_actual_O3_power_right_offset_and_R2_is_twice_Rc")
    if not all(old_proof.get(k) is True for k in required):
        raise ValueError("old radius source proof unavailable")
    snap = archive["Rp_native_C5_pulse_source_bindings"]
    inputs = snap["inputs"]
    q, c, o = (inputs[NAMES[name]] for name in ("quiet", "compact", "outer"))
    roles = {
        "quiet_R2_incoming_compact_terminal": equal_refs(q["R2"], q["incoming"]["radius"], c["terminal"]["radius"]),
        "compact_Rc_outer_O3_power": equal_refs(c["Rc"], o["endpoints"]["O3_power"]["radius"]),
        "quiet_outer_Tw": equal_refs(q["Tw"], o["Tw"]),
    }
    return {"status": "source_roles_returned", "source_family_matches": True,
            "old_C3_radius_gate_accepted": True, "old_C3_radius_proof_scope": "historical C3 quiet, not the current true quiet owner",
            "opaque_role_comparisons": roles,
            "roles_resolved": all(v["status"] != "unresolved_noninteger_or_missing_reference" for v in roles.values()),
            "all_equal_opaque_roles": all(v["status"] == "equal_opaque_refs" for v in roles.values()),
            "manual_source_relations": ["old_logRp=old_logRc+old_Tw-2",
                "current_R2=2*current_Rc under compact endpoint substitution",
                "current_Rp=current_Rc*exp(current_Tw-2) under quiet endpoint substitution"],
            "current_old_Rc_identity_proved": False, "current_old_Tw_identity_proved": False,
            "cross_owner_Rp_identity_proved": False, "source_values_computed": False,
            "full_archive_json_decoded": True, "exact_graph_nodes_selected_or_traversed": False,
            "protected_refs_used_as_graph_or_mapping_keys": False}


def fixture(hashes):
    path = ROOT / LEAF
    if sha(path.read_bytes()) != hashes[LEAF]:
        raise ValueError("new leaf changed after entry")
    spec = importlib.util.spec_from_file_location("rp_cross_owner_signed_gap_v1_new_leaf", path)
    if spec is None or spec.loader is None:
        raise ValueError("new leaf cannot be loaded")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    ledger = {"history": "unchanged", "P1_ne_m15": True}
    r = module.build_rp_cross_owner_signed_gap(
        C=Fraction(1), B=Fraction(1), logRc_old=Fraction(0), Tw_old=Fraction(2),
        d=Fraction(1, 3), t=Fraction(1, 2), f=Fraction(-1, 5), g=Fraction(3, 7),
        h=Fraction(2, 3), dZ=Fraction(2, 11), error_ledger=ledger)
    if (r.logscale != 0 or r.signed_mantissa_delta != 1
            or r.signed_mantissa_delta_Z1 != Fraction(1844, 1155)
            or r.position_mantissa != Fraction(1, 2) or r.duration_mantissa != Fraction(1, 2)
            or r.position_mantissa_Z1 != Fraction(832, 1155)
            or r.duration_mantissa_Z1 != Fraction(92, 105) or r.error_ledger is not ledger):
        raise ValueError("independent rational fixture failed")
    return {"kind": "synthetic_only", "logscale": str(r.logscale), "delta": str(r.signed_mantissa_delta),
            "ordinary_Z1": str(r.signed_mantissa_delta_Z1),
            "position_Z1": str(r.position_mantissa_Z1), "duration_Z1": str(r.duration_mantissa_Z1),
            "ledger_preserved_by_identity": True, "leaf_sha256_verified_after_entry": True}


def main():
    hashes, source_sha, brief_sha = preflight()
    entry_sha = exclusive(ENTRY, {"component": COMPONENT, "entered": True, "source_checkpoint_sha256": source_sha,
                                   "brief_sha256": brief_sha, "input_hashes": hashes, "consumer_calls": 0})
    output = {"component": COMPONENT, "schema": COMPONENT + "_result_v1", "entry_sha256": entry_sha,
              "source_checkpoint_sha256": source_sha, "source_reader_calls": 0, "fixture_calls": 0}
    try:
        output["source_reader_calls"] = 1
        output["source_roles"] = source_reader(hashes)
        output["fixture_calls"] = 1
        output["fixture"] = fixture(hashes)
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
