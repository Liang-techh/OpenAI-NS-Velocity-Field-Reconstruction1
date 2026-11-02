"""Run one safe corrected formal inlet Kchi tail consumer.

The historical formal receipt and its invalid j>=3 majorant are retained as
evidence only.  This runner consumes the saved core inlet and whole-Z axis
packet, invokes only the new safe leaf, and keeps the historical j>=3 and j>=2
errors separate.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path
from typing import Any

from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
NAME = "inlet_kchi_tail_safe_v1"
DIR = ROOT / "tmp" / NAME
SOURCE, ENTRY, RESULT = (DIR / item for item in ("source.json", "entry.json", "result.json"))

LEAF = "src/openai_ns_reconstruction/inlet_kchi_tail_safe_v1.py"
FIXTURE = "src/openai_ns_reconstruction/inlet_kchi_tail_safe_fixture_v1.py"
RUNNER = "scripts/run_inlet_kchi_tail_safe_v1_once.py"
FORMAL = "experiments/root_st073/formal_kchi_collar_v1_receipt.json.gz"
BASE = "experiments/root_st073/n1_core_six_inlet_v2_receipt.json.gz"
KCHI = "experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz"
BIND = "experiments/root_st073/n1_collar_axis_potential_binding_v1_receipt.json.gz"
FIRST = "experiments/root_st073/n1_collar_core_inlet_potential_v5_receipt.json.gz"
OLDTAIL = "experiments/root_st073/n1_collar_inlet_j2_tail_v3_receipt.json.gz"
SELECT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json"

RECEIPTS = (FORMAL, BASE, KCHI, BIND, FIRST, OLDTAIL)
REQUIRED = (LEAF, FIXTURE, RUNNER, FORMAL, BASE, KCHI, BIND, FIRST, OLDTAIL, SELECT)

PIN = {
    FORMAL: "80bcc7912f9cd28781061a94d50598224c30cc3f3071a03a9b5812c6e9a14d35",
    BASE: "59c0ff88652d83d7b92ba9c3a564ee90eac5ae794220b8b858c112936b76342d",
    KCHI: "b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0",
    BIND: "5859d5cb4bb82d372d8ac2f3ea1629e07c9787ffbaca149c99d1aba5edb8c6e8",
    FIRST: "1c644fb0f41a8ae407c09f56b00e7aee1720f91d08163a0c7483b4b1b6c61e7a",
    OLDTAIL: "dd28d5be7e090afc720f7b71d48ddbeb0e542e8e4d5097e86a86d78ab8c36891",
    SELECT: "0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e",
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def require(value: object, message: str) -> None:
    if not value:
        raise ValueError(message)


def save(path: Path, value: object) -> str:
    raw = (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()
    with path.open("xb") as handle:
        handle.write(raw)
    return sha(raw)


def encode_intervals(value: Any) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return {
            "lower_exact_mpf_tuple": list(raw[0]),
            "upper_exact_mpf_tuple": list(raw[1]),
            "directed": True,
        }
    if isinstance(value, dict):
        return {str(key): encode_intervals(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode_intervals(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError("cannot serialize result type " + type(value).__name__)


def module(path: str, name: str) -> Any:
    require(path in (LEAF, FIXTURE), "unexpected import path")
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    if spec is None or spec.loader is None:
        raise RuntimeError("new component unavailable: " + path)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def source_pins(receipt: dict[str, Any]) -> dict[str, str]:
    return {row["path"]: row["sha256"] for row in receipt.get("source", {}).get("files", [])}


def preflight() -> tuple[dict[str, str], str, str]:
    require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    brief_raw = BRIEF.read_bytes()
    source_raw = SOURCE.read_bytes()
    brief = json.loads(brief_raw)
    manifest = json.loads(source_raw)
    require(brief["component"] == NAME and brief["phase"] == NAME + "_pending", "brief guard")
    require(brief["goal_active"] is True and brief["full_goal_complete"] is False, "goal guard")
    require(brief["one_entry_only"] is True and brief["do_not_reenter"] is True, "exclusive-entry guard")
    checkpoint = brief["source_checkpoint"]
    require(checkpoint["sha256"] == sha(source_raw), "source checkpoint hash")
    require(Path(checkpoint["path"]).resolve() == SOURCE.resolve(), "source checkpoint path")
    require(sha(Path(checkpoint["preimage_path"]).read_bytes()) == checkpoint["preimage_sha256"], "preimage hash")
    require(manifest["component"] == NAME and manifest["one_entry_only"] is True, "manifest guard")
    require(len(manifest["files"]) == len(REQUIRED), "manifest file count")
    require({row["path"] for row in manifest["files"]} == set(REQUIRED), "manifest file set")
    hashes: dict[str, str] = {}
    for row in manifest["files"]:
        path = row["path"]
        require(path not in hashes, "duplicate manifest path: " + path)
        digest = sha((ROOT / path).read_bytes())
        require(digest == row["sha256"], "changed input: " + path)
        if path in PIN:
            require(digest == PIN[path], "changed pinned input: " + path)
        hashes[path] = digest
    return hashes, sha(source_raw), sha(brief_raw)


def _output(receipt: dict[str, Any], path: str) -> dict[str, Any]:
    result = receipt.get("result")
    require(isinstance(result, dict), "missing result: " + path)
    require(result.get("status") == "returned", "receipt not returned: " + path)
    output = result.get("output")
    require(isinstance(output, dict), "missing receipt output: " + path)
    return output


def _endpoint_pair(value: Any) -> tuple[tuple[Any, ...], tuple[Any, ...]]:
    current = value
    while isinstance(current, dict):
        if "lower_exact_mpf_tuple" in current or "upper_exact_mpf_tuple" in current:
            require(
                "lower_exact_mpf_tuple" in current and "upper_exact_mpf_tuple" in current,
                "incomplete exact interval",
            )
            lower = tuple(current["lower_exact_mpf_tuple"])
            upper = tuple(current["upper_exact_mpf_tuple"])
            require(len(lower) == len(upper) == 4, "interval endpoint tuples")
            return lower, upper
        nested = next((current[key] for key in ("coefficient", "interval", "enclosure", "box", "value") if key in current), None)
        require(nested is not None and nested is not current, "axis interval endpoints unavailable")
        current = nested
    raise TypeError("axis interval endpoints unavailable")


def _list_values(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def build_inlet(base: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    saved = base["baseline_G0g_fixed_exit"]["G0g"]["aggregate"]
    inlet: dict[str, dict[str, list[dict[str, Any]]]] = {"F": {}, "Fx": {}}
    index: dict[str, Any] = {}
    for field, component in (("F1", "F"), ("d_x_F1", "Fx")):
        for order in ("Z0", "Z1"):
            rows: list[dict[str, Any]] = []
            for item_index, item in enumerate(saved[order][field]):
                require(isinstance(item, dict), "saved inlet term must be a mapping")
                path = ["baseline_G0g_fixed_exit", "G0g", "aggregate", order, field, item_index]
                sid = "saved_core_G0g:" + field + ":" + order + ":" + str(item_index)
                require(sid not in index, "duplicate saved inlet source id: " + sid)
                origin = item.get("origin") if isinstance(item.get("origin"), dict) else {}
                parent_ids = _list_values(item.get("source_ids", [])) + _list_values(origin.get("source_ids", []))
                classes = _list_values(item.get("error_classes", [])) + _list_values(origin.get("error_classes", []))
                source_classes = _list_values(item.get("source_error_classes", [])) + _list_values(origin.get("source_error_classes", []))
                powers = item.get("powers", [0, item.get("F0_power", 0), 0, item.get("Lambda_power", 0)])
                require(isinstance(powers, (list, tuple)) and len(powers) == 4, "saved formal powers")
                original_paths = _list_values(item.get("source_paths", []))
                row = {
                        "coefficient": item["coefficient"],
                        "powers": list(powers),
                        "formal_powers": list(powers),
                        "source_ids": parent_ids + [sid],
                        "parent_source_ids": parent_ids,
                        "source_path": path,
                        "source_paths": original_paths + [path],
                        "error_classes": classes,
                        "source_error_classes": source_classes,
                        "uncertainty": True,
                        "original_term": item,
                    }
                for flag in ("error_only", "symmetric_error", "signed_field"):
                    if flag in item:
                        row[flag] = item[flag]
                    elif flag in origin:
                        row[flag] = origin[flag]
                rows.append(row)
                index[sid] = {"receipt": BASE, "path": path, "parent_source_ids": parent_ids}
            require(rows, "saved inlet rows are empty: " + component + "/" + order)
            inlet[component][order] = rows
    return inlet, index


def build_axis(kchi: dict[str, Any], binding: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    axis: list[dict[str, Any]] = []
    index: dict[str, Any] = {
        "axis_chi:whole_Z": {"receipt": KCHI, "binding_receipt": BIND, "scope": "whole Z"},
    }
    saved_binding = binding["saved_chi_Z0_Z1"]
    for order in (0, 1):
        label = "Z" + str(order)
        coefficient = kchi["chi"]["chi"][order]
        bound = saved_binding[str(order)]["coefficient"]
        require(_endpoint_pair(coefficient) == _endpoint_pair(bound), "axis binding mismatch: " + label)
        require(saved_binding[str(order)].get("source_receipt") == KCHI, "axis source receipt mismatch: " + label)
        path = ["chi", "chi", order]
        axis.append(
            {
                "coefficient": coefficient,
                "powers": [0, 0, 0, 0],
                "formal_powers": [0, 0, 0, 0],
                "source_ids": ["axis_chi:whole_Z", "axis_chi:" + label],
                "source_paths": [path],
                "source_path": path,
                "error_classes": ["whole_Z_axis_potential_interval"],
                "uncertainty": True,
                "error_only": False,
                "whole_Z": True,
                "scope": "same source axis Kchi, whole Z",
            }
        )
        index["axis_chi:" + label] = {
            "receipt": KCHI,
            "binding_receipt": BIND,
            "path": path,
        }
    return axis, index


def consume(hashes: dict[str, str]) -> dict[str, Any]:
    # Every byte is authenticated before any receipt is decompressed or any
    # new leaf/fixture is imported.
    raw: dict[str, bytes] = {}
    for path in REQUIRED:
        data = (ROOT / path).read_bytes()
        require(sha(data) == hashes[path], "changed after entry: " + path)
        raw[path] = data

    leaf = module(LEAF, NAME + "_leaf")
    fixture = module(FIXTURE, NAME + "_fixture")
    context = MPIntervalContext()
    context.prec = 1800
    proof = fixture.run_fixture(leaf.build, context)
    require(proof.get("passed") is True, "fixture did not pass")
    if isinstance(proof.get("checks"), dict):
        require(all(proof["checks"].values()), "fixture has failed checks")

    receipts = {path: json.loads(gzip.decompress(raw[path])) for path in RECEIPTS}
    selected = json.loads(raw[SELECT])
    outputs = {path: _output(receipts[path], path) for path in RECEIPTS if path != FORMAL}
    formal = receipts[FORMAL]["result"]
    require(formal["status"] == "returned", "historical formal receipt not returned")
    base = outputs[BASE]
    kchi = outputs[KCHI]
    binding = outputs[BIND]
    first = outputs[FIRST]
    oldtail = outputs[OLDTAIL]

    family = base["source_family"]
    source_logs = base["source_logs"]
    for path in (KCHI, BIND, FIRST, OLDTAIL):
        require(outputs[path]["source_family"] == family, "source family mismatch: " + path)
        require(outputs[path]["source_logs"] == source_logs, "source logs mismatch: " + path)
    require(selected["source_family"] == family, "selected source family mismatch")
    formal_scope = formal["saved_input_scope"]
    require(formal_scope["family"] == family, "formal saved family mismatch")
    require(formal_scope["source_logs"] == source_logs, "formal saved source logs mismatch")
    require(formal_scope["inlet_receipt"] == BASE, "formal inlet receipt mismatch")
    require(formal_scope["axis_receipt"] == KCHI, "formal axis receipt mismatch")
    require(formal_scope["axis_owner_evidence"] == BIND, "formal axis owner receipt mismatch")
    require(formal_scope["previous_j1_receipt"] == FIRST, "formal first-word receipt mismatch")
    require(formal_scope["previous_j2_tail_receipt"] == OLDTAIL, "formal j2 receipt mismatch")
    require(formal_scope["forcing_supplied"] is False, "formal forcing was supplied")
    require(formal_scope["same_live_owner_proved"] is False, "formal live owner was asserted")
    require(formal_scope["old_and_new_tail_combined"] is False, "formal tails were combined")

    require(binding["source_checks"]["one_query_jet_shared_in_integration"] is True, "axis query-jet binding")
    require(binding["phase_shear_cutoff_identified_with_axis_potential"] is False, "phase shear substituted for axis")
    contract = selected["actual_cutoff_complete_support_measure"]
    require(
        contract["passed"] is True
        and contract["measure_is_Z_independent"] is True
        and contract["source_family"] == family
        and "0<hb<1/2" in contract["inequality_proof"],
        "support-measure geometry contract",
    )

    formal_pins = source_pins(receipts[FORMAL])
    for path in (BASE, KCHI, BIND, FIRST, OLDTAIL, SELECT):
        require(formal_pins.get(path) == hashes[path], "e944 source pin mismatch: " + path)
    first_pins = source_pins(receipts[FIRST])
    for path in (BASE, KCHI):
        require(first_pins.get(path) == hashes[path], "first-word source pin mismatch: " + path)

    sc_bits = tuple(contract["source_sc_exact_mpf_tuple"])
    require(len(sc_bits) == 4 and all(type(item) is int for item in sc_bits), "source sc tuple")
    sc = context.make_mpf((sc_bits, sc_bits))
    require(sc.a > 0, "source sc must be positive")

    inlet, input_index = build_inlet(base)
    axis, axis_index = build_axis(kchi, binding)
    result = leaf.build(context, inlet, axis[0], axis[1], sc)
    require(isinstance(result, dict), "safe leaf result must be a mapping")
    if "tail_starts_at_j" in result:
        require(result["tail_starts_at_j"] == 3, "corrected tail must start at j=3")
    if isinstance(result.get("alternative_envelopes"), dict) and "combined" in result["alternative_envelopes"]:
        require(result["alternative_envelopes"]["combined"] is False, "corrected tail was combined")
    if isinstance(result.get("scope"), dict) and "ledger_written" in result["scope"]:
        require(result["scope"]["ledger_written"] is False, "safe tail wrote the ledger")

    return {
        "source_family": family,
        "source_logs": source_logs,
        "fixture": encode_intervals(proof),
        "corrected_tail": encode_intervals(result),
        "input_index": input_index,
        "axis_index": axis_index,
        "input_receipt_sha256": {path: hashes[path] for path in RECEIPTS + (SELECT,)},
        "e944_formal_reference": {
            "receipt": FORMAL,
            "sha256": hashes[FORMAL],
            "invalid_j_ge_3_tail_retained": True,
            "reentered": False,
        },
        "saved_input_scope": {
            "inlet_receipt": BASE,
            "axis_receipt": KCHI,
            "axis_owner_evidence": BIND,
            "previous_j1_receipt": FIRST,
            "previous_j2_tail_receipt": OLDTAIL,
            "historical_formal_receipt": FORMAL,
            "forcing_supplied": False,
            "corrected_tail_only": True,
            "same_live_owner_proved": False,
            "old_formal_j_ge_3_retained": True,
            "old_j2_tail_retained": True,
            "old_and_new_tail_combined": False,
            "old_errors_removed": False,
            "physical_transformation_applied": False,
            "m15_computed": False,
            "ledger_written": False,
            "ancestor_called": False,
            "previous_words_not_reentered": True,
        },
    }


def main() -> int:
    hashes, source_sha, brief_sha = preflight()
    entry_sha = save(
        ENTRY,
        {
            "component": NAME,
            "source_sha256": source_sha,
            "brief_sha256": brief_sha,
            "input_hashes": hashes,
            "entered": True,
            "ancestor_called": False,
            "previous_words_not_reentered": True,
            "old_formal_j_ge_3_retained": True,
            "old_j2_tail_retained": True,
            "old_and_new_tail_combined": False,
            "raw_gzip_authenticated_before_decompression": True,
        },
    )
    outcome: dict[str, Any] = {
        "component": NAME,
        "source_sha256": source_sha,
        "brief_sha256": brief_sha,
        "entry_sha256": entry_sha,
        "status": "failed",
        "raw_gzip_authenticated_before_decompression": True,
        "ancestor_called": False,
        "previous_words_not_reentered": True,
        "old_formal_j_ge_3_retained": True,
        "old_j2_tail_retained": True,
        "old_and_new_tail_combined": False,
        "fixture": None,
        "corrected_tail": None,
    }
    try:
        payload = consume(hashes)
        outcome.update(
            status="returned",
            output=payload,
            fixture=payload.get("fixture"),
            corrected_tail=payload.get("corrected_tail"),
            error=None,
        )
    except BaseException as exc:
        outcome["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
    result_sha = save(RESULT, outcome)
    print(
        json.dumps(
            {
                "component": NAME,
                "status": outcome["status"],
                "entry_sha256": entry_sha,
                "result_sha256": result_sha,
                "error": (outcome.get("error") or {}).get("message"),
            }
        )
    )
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
