"""One-entry independent first-word ordinary-Z2 core profile consumer.

The runner authenticates every frozen byte before decoding any saved receipt.
It runs the new synthetic fixture first, then adapts only the saved Z2
``outer_cells`` packet into ``build(c, cells)``.  No ancestor leaf or source
builder is imported or called.  Entry and result files are exclusive and a
failed fixture or leaf call is retained with its complete traceback.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import io
import json
import sys
import traceback
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
NAME = "core_firstword_z2_profile_v1"
DIR = ROOT / "tmp" / NAME
SOURCE, ENTRY, RESULT = (DIR / name for name in ("source.json", "entered.json", "result.json"))

LEAF = "src/openai_ns_reconstruction/core_firstword_z2_profile_v1.py"
FIXTURE = "src/openai_ns_reconstruction/core_firstword_z2_profile_fixture_v1.py"
RUNNER = "scripts/run_core_firstword_z2_profile_v1_once.py"
Z2 = "experiments/root_st073/n1_core_word1_z2_v2_receipt.json.gz"
LOW = "experiments/root_st073/n1_core_low_word_v3_receipt.json.gz"
SECOND = "experiments/root_st073/core_secondword_partial_v2_receipt.json.gz"

FILES = (LEAF, FIXTURE, RUNNER, Z2, LOW, SECOND)
REQUIRED = FILES
RECEIPTS = {Z2, LOW, SECOND}
PIN = {
    Z2: "c53ff32ad3ce5288279953e84bc7417b71477c82b00d94466d0d2880adbc7ac9",
    LOW: "5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6",
    SECOND: "d7de2dbcd1fbc90bf039d0b8527cc82be6e7826e475c2c0c1fdfa6ed3c8585fe",
}
ROWS = ("4", "5", "6")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def require(value: object, message: str) -> None:
    if not value:
        raise ValueError(message)


def save(path: Path, value: object) -> str:
    raw = (
        json.dumps(
            value,
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


def encode(value: object) -> object:
    """Convert interval objects to exact directed endpoint records for JSON."""

    if hasattr(value, "_mpi_"):
        return {
            "lower_exact_mpf_tuple": list(value._mpi_[0]),
            "upper_exact_mpf_tuple": list(value._mpi_[1]),
            "directed": True,
        }
    if isinstance(value, dict):
        return {str(key): encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError("unsupported output " + type(value).__name__)


def decode(c: MPIntervalContext, value: object) -> object:
    """Rebuild every directed endpoint pair in a saved packet as an mpmath iv."""

    if isinstance(value, dict):
        has_lower = "lower_exact_mpf_tuple" in value
        has_upper = "upper_exact_mpf_tuple" in value
        if has_lower or has_upper:
            require(has_lower and has_upper, "incomplete directed interval endpoints")
            lower = tuple(value["lower_exact_mpf_tuple"])
            upper = tuple(value["upper_exact_mpf_tuple"])
            require(
                len(lower) == len(upper) == 4
                and all(type(item) is int for item in lower + upper),
                "invalid exact interval endpoints",
            )
            result = c.make_mpf((lower, upper))
            require(result.a <= result.b, "reversed directed interval endpoints")
            return result
        return {key: decode(c, item) for key, item in value.items()}
    if isinstance(value, list):
        return [decode(c, item) for item in value]
    if isinstance(value, tuple):
        return [decode(c, item) for item in value]
    return value


def lookup(mapping: object, key: str) -> object:
    require(isinstance(mapping, dict), "expected row mapping")
    if key in mapping:
        return mapping[key]
    integer = int(key)
    if integer in mapping:
        return mapping[integer]
    raise KeyError(key)


def module(path: str, name: str):
    require(path in (LEAF, FIXTURE), "unexpected module import")
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    require(spec is not None and spec.loader is not None, "module spec")
    obj = importlib.util.module_from_spec(spec)
    sys.modules[name] = obj
    spec.loader.exec_module(obj)
    return obj


def receipt(raw: dict[str, bytes], path: str) -> dict:
    """Decode one receipt only after its raw bytes passed the frozen hash check."""

    value = json.loads(gzip.decompress(raw[path]))
    require(isinstance(value, dict), "receipt object " + path)
    require(value.get("result", {}).get("status") == "returned", "input receipt " + path)
    return value


def embedded_source_manifest(data: bytes, component: str) -> dict:
    """Read the compact receipt's first source-envelope line only."""

    with gzip.GzipFile(fileobj=io.BytesIO(data), mode="rb") as stream:
        line = stream.readline()
    prefix = (
        "{\"component\":"
        + json.dumps(component, ensure_ascii=False)
        + ",\"source\":"
    ).encode("utf-8")
    require(line.startswith(prefix), "compact source envelope " + component)
    manifest = json.loads(line[len(prefix) :])
    require(isinstance(manifest, dict) and isinstance(manifest.get("files"), list), "source manifest")
    return manifest


def manifest_hashes(manifest: dict) -> dict[str, str]:
    rows = manifest.get("files")
    require(isinstance(rows, list), "manifest files")
    result: dict[str, str] = {}
    for row in rows:
        require(isinstance(row, dict), "manifest row")
        path = row.get("path")
        digest = row.get("sha256")
        require(isinstance(path, str) and isinstance(digest, str), "manifest path/hash")
        require(path not in result, "duplicate manifest path " + path)
        result[path] = digest
    return result


def preflight() -> tuple[dict[str, str], str, str]:
    require(not ENTRY.exists() and not RESULT.exists(), "already entered")
    brief_raw = BRIEF.read_bytes()
    brief = json.loads(brief_raw)
    source_raw = SOURCE.read_bytes()
    source_sha256 = sha(source_raw)
    source = json.loads(source_raw)
    require(
        brief.get("component") == NAME
        and brief.get("phase") == NAME + "_pending"
        and brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("one_entry_only") is True
        and brief.get("do_not_reenter") is True,
        "brief guard",
    )
    checkpoint = brief["source_checkpoint"]
    require(Path(checkpoint["path"]).resolve() == SOURCE.resolve(), "source checkpoint path")
    require(checkpoint["sha256"] == sha(source_raw), "source checkpoint hash")
    preimage = Path(checkpoint["preimage_path"])
    require(sha(preimage.read_bytes()) == checkpoint["preimage_sha256"], "preimage hash")
    require(source.get("component") == NAME and source.get("one_entry_only") is True, "manifest guard")

    hashes: dict[str, str] = {}
    for row in source.get("files", []):
        path = row.get("path")
        require(path in REQUIRED and path not in hashes, "unexpected source file")
        digest = sha((ROOT / path).read_bytes())
        require(digest == row.get("sha256"), "changed source file " + path)
        if path in PIN:
            require(digest == PIN[path], "changed pinned input " + path)
        if path in RECEIPTS:
            require(row.get("raw_sha256_before_decode") is True, "receipt is not raw-pinned " + path)
        hashes[path] = digest
    require(set(hashes) == set(REQUIRED), "missing source file")
    return hashes, source_sha256, sha(brief_raw)


def f0_log(c: MPIntervalContext, logs: object, label: str):
    require(isinstance(logs, dict), label + " source_logs")
    if "F0_of_Z" in logs:
        return decode(c, logs["F0_of_Z"]), "F0_of_Z"
    if "F0_anchor" in logs:
        return decode(c, logs["F0_anchor"]), "F0_anchor"
    raise ValueError(label + " has neither F0_of_Z nor F0_anchor")


def build_cells(c: MPIntervalContext, z2_receipt: dict) -> list[dict]:
    output = z2_receipt["result"]["output"]
    saved = output["result"]["outer_cells"]
    require(isinstance(saved, list) and len(saved) == 8, "eight saved Z2 outer cells")
    cells: list[dict] = []
    for index, saved_cell in enumerate(saved):
        require(isinstance(saved_cell, dict), "saved Z2 cell")
        t = decode(c, saved_cell["t"])
        rho = decode(c, saved_cell["rho"])
        require(hasattr(t, "_mpi_") and hasattr(rho, "_mpi_"), "directed cell coordinates")
        rows: dict[str, list[dict]] = {}
        feedback = saved_cell["feedback_Z2"]
        for row in ROWS:
            terms = lookup(feedback, row)
            require(isinstance(terms, list), "feedback_Z2 row " + row)
            converted: list[dict] = []
            for term_index, saved_term in enumerate(terms):
                term = decode(c, saved_term)
                require(isinstance(term, dict), "feedback_Z2 term")
                # Keep every producer field, then attach an unambiguous local
                # identifier and the exact path in the saved receipt.
                term["input_id"] = f"z2:outer_cells:{index}:feedback_Z2:{row}:{term_index}"
                term["receipt_input_path"] = [
                    Z2,
                    "result",
                    "output",
                    "result",
                    "outer_cells",
                    index,
                    "feedback_Z2",
                    row,
                    term_index,
                ]
                converted.append(term)
            rows[row] = converted
        cells.append({"index": index, "t": t, "rho": rho, "feedback_Z2": rows})
    return cells


def build_missing_seam_cells(c: MPIntervalContext, low_receipt: dict) -> list[dict]:
    """Adapt the optional low-word B1 order-zero seam packet.

    The profile API remains ``build(c, cells)``.  A newer leaf may expose one
    separate ``build_missing_seam(c, profile, seam_cells)`` method; when present, the
    runner supplies only the saved B1 columns 5,1/5,2/5/3 and keeps their
    complete term provenance.  No seam method is synthesized here.
    """

    matrix_cells = low_receipt["result"]["output"]["matrix_cells"]
    require(isinstance(matrix_cells, list) and len(matrix_cells) == 8, "eight low matrix cells")
    seam_cells: list[dict] = []
    for index, saved_cell in enumerate(matrix_cells):
        require(isinstance(saved_cell, dict), "low matrix cell")
        cell = {
            "index": index,
            "t": decode(c, saved_cell["t"]),
            "rho": decode(c, saved_cell["rho"]),
            "B1": {"0": {}},
        }
        b1 = lookup(saved_cell, "B1")
        order_zero = lookup(b1, "0")
        for key in ("5,1", "5,2", "5,3"):
            terms = lookup(order_zero, key)
            require(isinstance(terms, list), "low B1 order-zero row " + key)
            converted: list[dict] = []
            for term_index, saved_term in enumerate(terms):
                term = decode(c, saved_term)
                require(isinstance(term, dict), "low B1 seam term")
                term["input_id"] = f"low:matrix_cells:{index}:B1:Z0:{key}:{term_index}"
                term["receipt_input_path"] = [
                    LOW,
                    "result",
                    "output",
                    "matrix_cells",
                    index,
                    "B1",
                    "0",
                    key,
                    term_index,
                ]
                converted.append(term)
            cell["B1"]["0"][key] = converted
        seam_cells.append(cell)
    return seam_cells


def consume(hashes: dict[str, str], out: dict) -> dict:
    # Re-read and authenticate all frozen bytes before importing the new leaf,
    # running its fixture, or decoding any receipt.
    raw: dict[str, bytes] = {}
    for path in REQUIRED:
        data = (ROOT / path).read_bytes()
        require(sha(data) == hashes[path], "changed after entry " + path)
        raw[path] = data

    c = MPIntervalContext()
    c.prec = 1800
    leaf = module(LEAF, NAME + "_leaf")
    fixture = module(FIXTURE, NAME + "_fixture")

    # The fixture is independent and runs before any saved receipt is decoded.
    try:
        fixture_runner = getattr(fixture, "run_fixture", None)
        require(callable(fixture_runner), "fixture must export run_fixture")
        proof = fixture_runner(leaf, c)
        out["fixture"] = encode(proof)
        require(isinstance(proof, dict) and proof.get("passed") is True, "synthetic fixture failed")
    except BaseException as exc:
        out["fixture_error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
        raise

    # The second-word receipt is a compact source-envelope anchor.  Reading
    # its first line authenticates the low-word parent without loading its
    # multi-megabyte result body.
    second_manifest = embedded_source_manifest(raw[SECOND], "core_secondword_partial_v2")
    second_files = manifest_hashes(second_manifest)
    require(second_files.get(LOW) == hashes[LOW], "secondword anchor does not pin low receipt")

    z2 = receipt(raw, Z2)
    low = receipt(raw, LOW)
    low_files = manifest_hashes(low["source"])
    require(low_files.get(Z2) == hashes[Z2], "low source manifest does not pin Z2 receipt")

    z2_output = z2["result"]["output"]
    low_output = low["result"]["output"]
    saved_word = z2_output["result"]
    require(saved_word["ordinary_Z_order"] == 2 and saved_word["first_collar_added"] is False and saved_word["allword_tail_added"] is False and saved_word["m15_computed"] is False, "saved Z2 scope")
    require(z2_output["source_family"] == low_output["source_family"], "source family mismatch")
    z2_f0, z2_f0_key = f0_log(c, z2_output["source_logs"], "Z2")
    low_f0, low_f0_key = f0_log(c, low_output["source_logs"], "low")
    require(z2_f0._mpi_ == low_f0._mpi_, "F0_of_Z/F0_anchor mismatch")
    z2_lambda = decode(c, z2_output["source_logs"]["Lambda"])
    low_lambda = decode(c, low_output["source_logs"]["Lambda"])
    require(z2_lambda._mpi_ == low_lambda._mpi_, "Lambda source-log mismatch")

    cells = build_cells(c, z2)
    out["actual_build_calls"] = 1
    profile_result = leaf.build(c, cells)
    require(isinstance(profile_result, dict), "leaf result must be a mapping")
    require(profile_result.get("schema") == NAME, "result schema")
    require(len(profile_result["target_cells"]) == 8 and profile_result["complete_firstword_Z2"] is False and profile_result["no_cancellation"] is True, "profile scope")

    seam_cells: list[dict] | None = None
    seam_result: object = {"status": "not_exposed"}
    seam_builder = getattr(leaf, "build_missing_seam", None)
    if callable(seam_builder):
        seam_cells = build_missing_seam_cells(c, low)
        out["actual_missing_seam_calls"] = 1
        seam_result = seam_builder(c, profile_result, seam_cells)
        require(isinstance(seam_result, dict), "missing seam result must be a mapping")
        require(seam_result.get("missing_branch_included") is True and seam_result.get("old_records_merged") is False and seam_result.get("no_cancellation") is True and set(seam_result["terms_by_order"]["Z1"]) == {"F1", "d_x_F1"}, "seam scope")
        require(all(seam_result["terms_by_order"]["Z1"][field] for field in ("F1", "d_x_F1")), "seam terms missing")
    else:
        out["actual_missing_seam_calls"] = 0

    return {
        "source_family": low_output["source_family"],
        "source_logs": low_output["source_logs"],
        "z2_source_logs": z2_output["source_logs"],
        "f0_log_keys_compared": {"z2": z2_f0_key, "low": low_f0_key},
        "anchor_receipt": {
            "path": SECOND,
            "sha256": hashes[SECOND],
            "source_manifest": second_manifest,
        },
        "input_receipt_sha256": {path: hashes[path] for path in (Z2, LOW, SECOND)},
        "input_cells": encode(cells),
        "result": encode(profile_result),
        "missing_seam_inputs": encode(seam_cells) if seam_cells is not None else None,
        "missing_seam_result": encode(seam_result),
        "actual_build_calls": 1,
        "actual_missing_seam_calls": out["actual_missing_seam_calls"],
        "z2_recomputed": False,
        "secondword_anchor_reused": True,
        "allword_tail_added": False,
        "m15_computed": False,
        "physical_owner_proved": False,
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
            "entered": True,
        },
    )
    out = {
        "component": NAME,
        "source_sha256": source_sha256,
        "entry_sha256": entry_sha256,
        "status": "failed",
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
    print(json.dumps({"component": NAME, "status": out["status"], "result_sha256": result_sha256}))
    return 0 if out["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())



