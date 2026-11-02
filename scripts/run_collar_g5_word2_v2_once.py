"""One-entry repaired local-g5 signed word-2 consumer.

The failed v1 archive is authenticated and retained as evidence only.  This
runner never imports or re-enters the v1 leaf, fixture, runner, or result.
All new raw inputs are authenticated before gzip decompression or import; the
exclusive entry is saved before the fixture and repaired leaf are called.
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
NAME = "collar_g5_word2_v2"
DIR = ROOT / "tmp" / NAME
SOURCE, ENTRY, RESULT = (DIR / item for item in ("source.json", "entry.json", "result.json"))

LEAF = "src/openai_ns_reconstruction/original_collar_g5_word2_v2.py"
FIXTURE = "src/openai_ns_reconstruction/original_collar_g5_word2_fixture_v2.py"
RUNNER = "scripts/run_collar_g5_word2_v2_once.py"
G5 = "experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz"
KCHI = "experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz"
BIND = "experiments/root_st073/n1_collar_axis_potential_binding_v1_receipt.json.gz"
FIRST = "experiments/root_st073/n1_collar_kchi_source_first_v3_receipt.json.gz"
OLD_TAIL = "experiments/root_st073/n1_collar_local_g5_j2_tail_v1_receipt.json.gz"
PRESS = "experiments/root_st073/n1_first_collar_g4_pressure_v1_receipt.json.gz"
SELECT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json"
FAILURE = "experiments/root_st073/collar_g5_word2_v1_failure.json.gz"
FAILURE_SHA = "c9a53044de6c8623b9ad2e63742c329a47a98403035cf8fe5dfb8791b4e32eaa"

V1_LEAF = "src/openai_ns_reconstruction/original_collar_g5_word2_v1.py"
V1_FIXTURE = "src/openai_ns_reconstruction/original_collar_g5_word2_fixture_v1.py"
V1_RUNNER = "scripts/run_collar_g5_word2_v1_once.py"

REQUIRED = (LEAF, FIXTURE, RUNNER, G5, KCHI, BIND, FIRST, OLD_TAIL, PRESS, SELECT, FAILURE)
PIN = {
    G5: "96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2",
    KCHI: "b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0",
    BIND: "5859d5cb4bb82d372d8ac2f3ea1629e07c9787ffbaca149c99d1aba5edb8c6e8",
    FIRST: "bdfd705e43db41badd69b737a014f9ff288eaa58a18060fe736aebb196f18262",
    OLD_TAIL: "67894bac58f4fbbcd5437bf970dadafeb918fe2d9638425f8e5165d21af4530b",
    PRESS: "c229cf188d306413f6edd0641a729321ad4567c2587ebf4701ada47720853660",
    SELECT: "0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e",
    FAILURE: FAILURE_SHA,
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


def interval(context: MPIntervalContext, value: Any) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return context.make_mpf((raw[0], raw[1]))
    if isinstance(value, dict):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            require("lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value, "incomplete interval")
            lower = tuple(value["lower_exact_mpf_tuple"])
            upper = tuple(value["upper_exact_mpf_tuple"])
            require(len(lower) == len(upper) == 4, "interval endpoint tuples")
            result = context.make_mpf((lower, upper))
            require(result.a <= result.b, "reversed interval")
            return result
        for key in ("coefficient", "coefficient_interval", "interval", "enclosure", "box", "value"):
            if key in value:
                return interval(context, value[key])
    if isinstance(value, (tuple, list)) and len(value) == 2:
        return context.make_mpf((value[0], value[1]))
    return context.mpf(value)


def module(path: str, name: str) -> Any:
    # The v1 module names are deliberately absent from every import path.
    require(path not in (V1_LEAF, V1_FIXTURE, V1_RUNNER), "v1 component re-entry is forbidden")
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
    require(LEAF not in (V1_LEAF, V1_FIXTURE, V1_RUNNER), "v1 leaf path in v2 manifest")
    hashes: dict[str, str] = {}
    for row in manifest["files"]:
        path = row["path"]
        digest = sha((ROOT / path).read_bytes())
        require(digest == row["sha256"], "changed input: " + path)
        if path in PIN:
            require(digest == PIN[path], "changed pinned input: " + path)
        hashes[path] = digest
    require(hashes[FAILURE] == FAILURE_SHA, "v1 failure archive pin")
    return hashes, sha(source_raw), sha(brief_raw)


def _parents(item: dict[str, Any]) -> tuple[dict[str, Any], list[Any], list[Any]]:
    origin = item.get("origin") if isinstance(item.get("origin"), dict) else {}
    origin_ids = origin.get("source_ids", [])
    item_ids = item.get("source_ids", [])
    origin_ids = list(origin_ids) if isinstance(origin_ids, (list, tuple)) else [origin_ids]
    item_ids = list(item_ids) if isinstance(item_ids, (list, tuple)) else [item_ids]
    parent_ids: list[Any] = []
    for value in origin_ids + item_ids:
        if value not in parent_ids:
            parent_ids.append(value)
    origin_classes = origin.get("error_classes", [])
    item_classes = item.get("error_classes", [])
    origin_classes = list(origin_classes) if isinstance(origin_classes, (list, tuple)) else [origin_classes]
    item_classes = list(item_classes) if isinstance(item_classes, (list, tuple)) else [item_classes]
    classes: list[Any] = []
    for value in origin_classes + item_classes:
        if value not in classes:
            classes.append(value)
    return origin, parent_ids, classes


def _g5_rows(context: MPIntervalContext, forcing: dict[str, Any]) -> tuple[dict[str, list[dict[str, Any]]], dict[str, Any]]:
    output = forcing["result"]["output"]
    rows: dict[str, list[dict[str, Any]]] = {"Z0": [], "Z1": []}
    index: dict[str, Any] = {}
    for order in (0, 1):
        label = "Z" + str(order)
        saved = output["g"][str(order)]["5"]
        require(len(saved) == 16, "expected 16 saved g5 sectors for " + label)
        for item_index, item in enumerate(saved):
            require(isinstance(item, dict), "g5 sector must be a mapping")
            origin, parent_ids, classes = _parents(item)
            require(parent_ids, "g5 sector has no parent source ids")
            sid = f"collar_g5:Z{order}:{item_index}"
            path = ["result", "output", "g", str(order), "5", item_index]
            require(sid not in index, "duplicate g5 source id: " + sid)
            powers = item.get("powers")
            require(isinstance(powers, (list, tuple)) and len(powers) == 4, "g5 formal powers")
            require(
                all(isinstance(power, int) and not isinstance(power, bool) for power in powers)
                and list(powers[1:]) == [1, 0, -1],
                "g5 formal power basis",
            )
            original_paths = item.get("source_paths", [])
            if not isinstance(original_paths, list):
                original_paths = [original_paths]
            row = dict(item)
            row.update(
                {
                    "coefficient": interval(context, item["coefficient"]),
                    "powers": list(powers),
                    "formal_powers": list(powers),
                    "source_ids": parent_ids + [sid],
                    "parent_source_ids": parent_ids,
                    "source_path": path,
                    "source_paths": original_paths + [path],
                    "frame": "raw physical g5",
                    "original_origin": origin,
                    "error_classes": classes,
                    "uncertainty": bool(item.get("uncertainty", origin.get("uncertainty", False))),
                    "error_only": bool(item.get("error_only", origin.get("error_only", False))),
                    "symmetric_error": bool(item.get("symmetric_error", origin.get("symmetric_error", False))),
                }
            )
            rows[label].append(row)
            index[sid] = {
                "receipt": G5,
                "path": path,
                "parent_source_ids": parent_ids,
                "original_origin": origin,
                "error_classes": classes,
            }
    return rows, index


def _axis_term(context: MPIntervalContext, kchi: dict[str, Any], binding: dict[str, Any], order: int) -> dict[str, Any]:
    label = "Z" + str(order)
    saved = kchi["result"]["output"]["chi"]["chi"][order]
    bound = binding["result"]["output"]["saved_chi_Z0_Z1"][str(order)]["coefficient"]
    require(interval(context, saved)._mpi_ == interval(context, bound)._mpi_, "axis chi binding mismatch " + label)
    common_id = "axis_chi:whole_Z"
    order_id = "axis_chi:" + label
    path = ["result", "output", "chi", "chi", order]
    return {
        "coefficient": interval(context, saved),
        "powers": [0, 0, 0, 0],
        "formal_powers": [0, 0, 0, 0],
        "source_ids": [common_id, order_id],
        "source_paths": [path],
        "source_path": path,
        "error_classes": ["whole_Z_axis_potential_interval"],
        "uncertainty": True,
        "error_only": False,
        "symmetric_error": False,
        "whole_Z": True,
        "scope": "same source axis Kchi, whole Z",
    }


def consume(hashes: dict[str, str]) -> dict[str, Any]:
    # Authenticate every byte, including the failed v1 archive, before gzip
    # decompression or importing the repaired leaf/fixture.
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

    receipts = {
        path: json.loads(gzip.decompress(raw[path]))
        for path in (G5, KCHI, BIND, FIRST, OLD_TAIL, PRESS)
    }
    v1_failure = json.loads(gzip.decompress(raw[FAILURE]))
    require(isinstance(v1_failure, dict), "v1 failure archive is not a JSON object")
    require(v1_failure["result"]["component"] == "collar_g5_word2_v1", "wrong v1 failure archive component")
    require(v1_failure["result"]["status"] == "failed", "v1 archive is not a failed result")
    require(all(receipt["result"]["status"] == "returned" for receipt in receipts.values()), "ancestor receipt status")
    forcing = receipts[G5]["result"]["output"]
    kchi = receipts[KCHI]["result"]["output"]
    binding = receipts[BIND]["result"]["output"]
    first = receipts[FIRST]
    old_tail = receipts[OLD_TAIL]
    selected = json.loads(raw[SELECT])
    family = forcing["source_family"]
    source_logs = forcing["source_logs"]
    for name, value in (
        ("kchi", kchi),
        ("binding", binding),
        ("first", first["result"]["output"]),
        ("old_tail", old_tail["result"]["output"]),
        ("pressure", receipts[PRESS]["result"]["output"]),
    ):
        require(value["source_family"] == family and value["source_logs"] == source_logs, "source family/log mismatch: " + name)
    require(selected["source_family"] == family, "selected source family mismatch")
    require(binding["phase_shear_cutoff_identified_with_axis_potential"] is False, "phase shear was used as axis chi")
    require(binding["source_checks"]["one_query_jet_shared_in_integration"] is True, "shared axis query-jet proof")

    first_pins = source_pins(first)
    old_tail_pins = source_pins(old_tail)
    for path in (G5, KCHI, BIND, PRESS, SELECT):
        require(first_pins[path] == hashes[path], "first source pin mismatch: " + path)
        require(old_tail_pins[path] == hashes[path], "old-tail source pin mismatch: " + path)
    require(old_tail_pins[FIRST] == hashes[FIRST], "old-tail first-word pin mismatch")
    require(source_pins(receipts[PRESS])[G5] == hashes[G5], "pressure forcing pin mismatch")
    historical = old_tail["result"]["output"]
    require(
        historical["raw_g5_frame_proven_from_pinned_source"] is True
        and historical["F0_derivative_not_reapplied"] is True
        and historical["hb_not_materialized"] is True,
        "historical raw frame contract",
    )

    contract = selected["actual_cutoff_complete_support_measure"]
    require(
        contract["passed"] is True
        and contract["measure_is_Z_independent"] is True
        and contract["source_family"] == family,
        "physical geometry certificate",
    )
    require("0<hb<1/2" in contract["inequality_proof"], "strict hb certificate")
    sc_bits = tuple(contract["source_sc_exact_mpf_tuple"])
    require(len(sc_bits) == 4, "source sc tuple")
    sc = context.make_mpf((sc_bits, sc_bits))
    require(sc.a > 0 and sc._mpi_ == interval(context, receipts[PRESS]["result"]["output"]["sc"])._mpi_, "same exact sc")
    g5, g5_index = _g5_rows(context, receipts[G5])
    chi0 = _axis_term(context, receipts[KCHI], receipts[BIND], 0)
    chi1 = _axis_term(context, receipts[KCHI], receipts[BIND], 1)
    axis_index = {
        "axis_chi:whole_Z": {"receipt": KCHI, "binding_receipt": BIND, "scope": "whole Z"},
        "axis_chi:Z0": {"receipt": KCHI, "binding_receipt": BIND, "path": ["result", "output", "chi", "chi", 0]},
        "axis_chi:Z1": {"receipt": KCHI, "binding_receipt": BIND, "path": ["result", "output", "chi", "chi", 1]},
    }
    result = leaf.build(context, g5, chi0, chi1, sc)
    require(isinstance(result, dict), "leaf result must be a mapping")
    return {
        "source_family": family,
        "source_logs": source_logs,
        "source_index": {**g5_index, **axis_index},
        "fixture": encode_intervals(proof),
        "g5_input_counts": {label: len(rows) for label, rows in g5.items()},
        "axis_source_ids": {
            "common": "axis_chi:whole_Z",
            "Z0": list(chi0["source_ids"]),
            "Z1": list(chi1["source_ids"]),
        },
        "selected_sc_exact_mpf_tuple": list(sc_bits),
        "g5_word2": encode_intervals(result),
        "v1_failure_archive": {
            "path": FAILURE,
            "sha256": hashes[FAILURE],
            "component": v1_failure["result"]["component"],
            "status": v1_failure["result"]["status"],
            "reentered": False,
        },
        "input_receipt_sha256": {path: hashes[path] for path in (G5, KCHI, BIND, FIRST, OLD_TAIL, PRESS, SELECT, FAILURE)},
        "saved_input_scope": {
            "g5_receipt": G5,
            "axis_receipt": KCHI,
            "axis_owner_evidence": BIND,
            "previous_first_receipt": FIRST,
            "previous_j2_tail_receipt": OLD_TAIL,
            "pressure_receipt": PRESS,
            "v1_failure_archive": FAILURE,
            "n2_conditional_source_enclosure_whole_collar_only": True,
            "raw_g5_to_physical_h5_identity_proved": False,
            "same_live_owner_proved": False,
            "j_ge_3_is_alternative_error_tail": True,
            "old_j2_tail_retained": True,
            "old_and_new_tail_combined": False,
            "old_errors_removed": False,
            "physical_transformation_applied": False,
            "m15_computed": False,
            "ledger_written": False,
            "ancestor_called": False,
            "v1_reentered": False,
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
            "v1_reentered": False,
            "v1_failure_archive": FAILURE,
            "raw_gzip_authenticated_before_decompression": True,
        },
    )
    outcome: dict[str, Any] = {
        "component": NAME,
        "source_sha256": source_sha,
        "entry_sha256": entry_sha,
        "status": "failed",
        "raw_gzip_authenticated_before_decompression": True,
        "ancestor_called": False,
        "v1_reentered": False,
        "fixture": None,
        "g5_word2": None,
    }
    try:
        payload = consume(hashes)
        outcome.update(
            status="returned",
            output=payload,
            fixture=payload.get("fixture"),
            g5_word2=payload.get("g5_word2"),
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
