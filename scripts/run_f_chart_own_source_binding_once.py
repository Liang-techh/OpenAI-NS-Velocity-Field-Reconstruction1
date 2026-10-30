"""One-entry direct-role comparison of old own-field and current F snapshots.

The frozen source consumer compares only opaque references. It does not
evaluate a source, index a graph, infer an identity from a node number, or
call an ancestor. Both full gzip JSON documents are parsed after raw hashing.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
COMPONENT = "f_chart_own_source_binding_v1"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE, ENTRY, RESULT = (DIRECTORY / x for x in ("source.json", "entered.json", "result.json"))
BASE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_"
FILES = {
    "scripts/run_f_chart_own_source_binding_once.py": None,
    BASE + "Rp_native_C5_pulse_source.json.gz": "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe",
    BASE + "Rp_native_C5_pulse_source_check.json": "320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15",
    BASE + "n1_pulse_own_fields.json.gz": "2c0f83a9ac36f4ea7970b7938bac31d36af7c19df03e9c1584a6a11aee14713b",
    BASE + "n1_pulse_own_fields_check.json": "0bc8229a1716871d896161caac907c572f7cd0a58540900aaa00869d3a920fb1",
}
F_ARCHIVE, F_RECEIPT, OWN_ARCHIVE, OWN_RECEIPT = list(FILES)[1:]
PULSE5 = "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions.json.gz"
PROTECTED = frozenset({0, 1, 17, 33572, 25220, 263344, 263348, 263349, 262918, 205147, 263346, 33982})
FIELDS = ("U", "U_y", "m", "m_y", "Q", "Q_y")


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def write_new(path, data):
    raw = (json.dumps(data, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode()
    with path.open("xb") as stream:
        stream.write(raw)
    return digest(raw)


def preflight():
    if ENTRY.exists() or RESULT.exists():
        raise ValueError("component already entered")
    brief_raw = BRIEF.read_bytes()
    brief = json.loads(brief_raw)
    if (brief.get("component") != COMPONENT or brief.get("phase") != COMPONENT + "_pending"
            or brief.get("goal_active") is not True or brief.get("full_goal_complete") is not False
            or brief.get("do_not_reenter") is not True or brief.get("one_entry_only") is not True):
        raise ValueError("pending brief guard failed")
    cp = brief["source_checkpoint"]
    checkpoint_path = Path(cp["path"])
    if not checkpoint_path.is_absolute() or checkpoint_path.resolve() != SOURCE.resolve():
        raise ValueError("checkpoint path mismatch")
    source_raw = SOURCE.read_bytes()
    if digest(source_raw) != cp["sha256"]:
        raise ValueError("checkpoint raw hash mismatch")
    source = json.loads(source_raw)
    if (source.get("component") != COMPONENT or source.get("schema") != COMPONENT + "_checkpoint_v1"
            or source.get("one_entry_only") is not True):
        raise ValueError("source checkpoint scope mismatch")
    entries = source.get("files")
    if not isinstance(entries, list) or len(entries) != len(FILES):
        raise ValueError("five frozen raw inputs required")
    hashes = {}
    for item in entries:
        path, expected = item["path"], item["sha256"]
        if path not in FILES or path in hashes:
            raise ValueError("unexpected or duplicate input")
        raw = (ROOT / path).read_bytes()
        actual = digest(raw)
        if actual != expected or (FILES[path] is not None and actual != FILES[path]):
            raise ValueError("raw hash mismatch: " + path)
        hashes[path] = actual
    if set(hashes) != set(FILES):
        raise ValueError("frozen inputs incomplete")
    return source, hashes, digest(source_raw), digest(brief_raw)


def read_pinned(path, hashes):
    raw = (ROOT / path).read_bytes()
    if digest(raw) != hashes[path]:
        raise ValueError("raw input changed after entry: " + path)
    if path.endswith(".gz"):
        raw = gzip.decompress(raw)
    result = json.loads(raw)
    if not isinstance(result, dict):
        raise ValueError("source object required: " + path)
    return result


def compare(left, right):
    valid = lambda x: isinstance(x, int) and not isinstance(x, bool) and x >= 0
    if not valid(left) or not valid(right):
        return {"status": "unresolved_missing_or_noninteger_reference"}
    return {"status": "equal_opaque_integer" if left == right else "different_opaque_integers",
            "old_ref": left, "current_ref": right,
            "old_protected": left in PROTECTED, "current_protected": right in PROTECTED}


def direct(mapping, *path):
    value = mapping
    for key in path:
        if not isinstance(value, dict) or key not in value:
            return None
        value = value[key]
    return value


def source_consumer(hashes):
    fa = read_pinned(F_ARCHIVE, hashes)
    fr = read_pinned(F_RECEIPT, hashes)
    oa = read_pinned(OWN_ARCHIVE, hashes)
    ore = read_pinned(OWN_RECEIPT, hashes)
    family = fa.get("source_family")
    if (not isinstance(family, dict) or not family or fr.get("source_family") != family
            or oa.get("source_family") != family or ore.get("source_family") != family
            or fr.get("all_passed") is not True or ore.get("all_passed") is not True):
        raise ValueError("source-family/acceptance mismatch")
    fpins, opins = fr.get("input_hashes", {}), ore.get("input_hashes", {})
    for path in (F_ARCHIVE, OWN_ARCHIVE, OWN_RECEIPT):
        if fpins.get(Path(path).name) != hashes[path]:
            raise ValueError("F receipt source pin missing: " + path)
    if opins.get(Path(OWN_ARCHIVE).name) != hashes[OWN_ARCHIVE]:
        raise ValueError("own receipt archive pin mismatch")
    pulse5_sha = "bdfcc70f57d8cea981f508074690a7e940f884382c2452805f126a3f96bdb59d"
    archived_pins = fa.get("input_hashes", {})
    if fpins.get(PULSE5) != pulse5_sha or archived_pins.get(PULSE5) != pulse5_sha:
        raise ValueError("F source pulse5 snapshot pin mismatch")
    for path in (OWN_ARCHIVE, OWN_RECEIPT):
        if archived_pins.get(Path(path).name) != hashes[path]:
            raise ValueError("F archive own source pin mismatch: " + path)
    old = oa.get("own_main_gap_end_fields")
    snap = fa.get("Rp_native_C5_pulse_source_bindings")
    if not isinstance(old, dict) or not isinstance(snap, dict):
        raise ValueError("saved source role containers missing")
    results = {}
    for chart in ("main", "gap", "end"):
        coordinates = {key: compare(direct(old, chart, key), direct(snap, "pulse_coordinates", chart, key))
                       for key in ("query", "y", "radius")}
        jets = {}
        for field in FIELDS:
            a = direct(old, chart, field)
            b = direct(snap, "inputs", PULSE5, "leading_fields", chart, field)
            jets[field] = {str(i): compare(a[i] if isinstance(a, list) and len(a) > i else None,
                                           b[i] if isinstance(b, list) and len(b) > i else None)
                           for i in (0, 1)}
        results[chart] = {"coordinate_roles": coordinates, "old_to_pulse5_field_rows": jets}
    return {
        "schema": COMPONENT, "status": "direct_opaque_role_comparisons_returned", "charts": results,
        "source_family_matches": True, "receipt_pins_match": True,
        "both_entire_archive_json_documents_decoded": True,
        "exact_graph_nodes_selected_or_traversed": False,
        "protected_refs_used_as_graph_or_mapping_keys": False,
        "common_graph_namespace_independently_proved": False,
        "F_field_radius_owner_proved": False, "source_function_identity_proved": False,
        "ordinary_Z_semantics_proved": False, "endpoint_identity_proved": False,
        "source_values_computed": False, "numerical_values_computed": False,
        "old_errors_removed": False, "global_closure_proved": False,
    }


def main():
    source, hashes, source_sha, brief_sha = preflight()
    entry_sha = write_new(ENTRY, {"schema": COMPONENT + "_entry_v1", "component": COMPONENT,
        "source_checkpoint_sha256": source_sha, "brief_sha256": brief_sha,
        "input_hashes": hashes, "entered": True, "source_consumer_invocations": 0})
    try:
        result = source_consumer(hashes)
        output = {"schema": COMPONENT + "_result_v1", "component": COMPONENT, "status": "returned",
                  "entry_sha256": entry_sha, "source_checkpoint_sha256": source_sha,
                  "source_consumer_invocations": 1, "result": result, "error": None}
    except BaseException as exc:
        output = {"schema": COMPONENT + "_result_v1", "component": COMPONENT, "status": "failed",
                  "entry_sha256": entry_sha, "source_checkpoint_sha256": source_sha,
                  "source_consumer_invocations": 1, "result": None,
                  "error": {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()}}
    result_sha = write_new(RESULT, output)
    print(json.dumps({"component": COMPONENT, "status": output["status"], "entry_sha256": entry_sha,
                      "source_sha256": source_sha, "result_sha256": result_sha}))
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
