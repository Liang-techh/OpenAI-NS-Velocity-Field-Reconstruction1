"""Package and inspect saved loop-density evidence; never execute mathematics."""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
from pathlib import Path

from run_original_Rm_second_radial_identity import ROOT, BASE, PREFIX, require, sha, path_for, read_pinned, read_saved, write_new

SOURCE = BASE / (PREFIX + "Rm_loop_density_identity.json")
SOURCE_HASH = "529efebffc396163285414c2fe27846f13a9f2bba6441864f86dfff4baa14f7b"
SOURCE_BYTES = 239575760
ARCHIVE = SOURCE.with_suffix(".json.gz")
SELF = "scripts/archive_original_Rm_loop_density_report.py"
PARENT = PREFIX + "Rm_generic_native_identity.json"
PARENT_HASH = "71db4b9121e591878b4047025ed6d98e2a9a78d3659d82430ecfc294c0cafbf2"
CHECKPOINTS = ROOT / "tmp/original_Rm_loop_density_identity"
PACK_RECEIPT = CHECKPOINTS / "new_archive_packaging.json"
RESULT_HASH = "08ee51c4e8657017dad693a8bc1ae8c08a503c820957b7db4335b52e9d74fc1d"
PHASE_HASH = "02ff3c2e0678e11d6c913e87e0a4eaba5728a5c2d34a6870285e2b1eb5ce6970"


def pack():
    raw = SOURCE.read_bytes()
    require(len(raw) == SOURCE_BYTES and hashlib.sha256(raw).hexdigest() == SOURCE_HASH, "saved plain report changed")
    require(not ARCHIVE.exists() and not PACK_RECEIPT.exists(), "archive already packaged")
    report = json.loads(raw)
    require(report["actual_Rm_loop_density_identity"]["passed"] is True
            and report["loop_density_checkpoint_hashes"]["new_loop_density_identity.json"] == RESULT_HASH
            and report["loop_density_checkpoint_hashes"]["new_common_phase_source.json"] == PHASE_HASH,
            "actual saved result changed")
    # Only packaging metadata/pins are added. Every mathematical packet is
    # preserved as parsed JSON data; no expression reducer is imported.
    report["implementation_hashes"][SELF] = sha(ROOT / SELF)
    report["archive_packaging_audit"] = dict(source_plain_sha256=SOURCE_HASH, source_plain_bytes=SOURCE_BYTES,
        canonical_result_checkpoint_sha256=RESULT_HASH, phase_checkpoint_sha256=PHASE_HASH,
        lossless_saved_data_packaging=True, mathematical_packets_modified=False,
        new_or_accepted_math_operator_executed=False, reason="expanded evidence exceeds GitHub single-file limit")
    decoded = (json.dumps(report, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode()
    compressed = gzip.compress(decoded, compresslevel=6, mtime=0)
    with ARCHIVE.open("xb") as stream:
        stream.write(compressed)
    receipt = dict(passed=True, archive_sha256=hashlib.sha256(compressed).hexdigest(), archive_bytes=len(compressed),
        decoded_sha256=hashlib.sha256(decoded).hexdigest(), decoded_bytes=len(decoded),
        original_plain_sha256=SOURCE_HASH, original_plain_bytes=SOURCE_BYTES,
        implementation_sha256=sha(ROOT / SELF), mathematical_operator_executed=False)
    write_new(PACK_RECEIPT, receipt)
    return receipt


def inspect():
    receipt = read_saved(PACK_RECEIPT)
    require(receipt["implementation_sha256"] == sha(ROOT / SELF), "archive inspection source changed")
    compressed = ARCHIVE.read_bytes()
    require(hashlib.sha256(compressed).hexdigest() == receipt["archive_sha256"]
            and len(compressed) == receipt["archive_bytes"], "compressed evidence changed before decode")
    decoded = gzip.decompress(compressed)
    require(hashlib.sha256(decoded).hexdigest() == receipt["decoded_sha256"]
            and len(decoded) == receipt["decoded_bytes"], "lossless decoded evidence changed")
    report = json.loads(decoded)
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    phase, result = report["actual_Rm_common_phase_source"], report["actual_Rm_loop_density_identity"]
    require(phase["passed"] is True and result["passed"] is True and len(result["rows"]) == 22
            and all(row["status"] == "proved" for row in result["rows"]) and not result["unresolved_gates"], "source identity incomplete")
    packet = report["actual_Rm_loop_density_inputs"]
    require(report["parent_report"] == PARENT and report["parent_report_sha256"] == PARENT_HASH
            and phase["source_graph_digest"] == packet["graph_digest"]
            and result["nodes_digest"] == packet["graph_digest"], "source graph or parent changed")
    for key in ("actual_Rm_generic_native_identity", "actual_Rm_history_native_source", "actual_Rm_absolute_m15",
                "corrected_own_swirl_provider", "corrected_reserved_point_numeric", "current_Rm_pressure_inlet_identity"):
        require(report[key] == parent[key], "parent packet changed: " + key)
    for name, expected in {**report["input_hashes"], **report["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "pinned source changed: " + name)
    paths = ("new_source.json", "new_phase_operator_entered.json", "new_common_phase_source.json",
             "new_loop_operator_entered.json", "new_loop_density_identity.json")
    require(report["loop_density_checkpoint_hashes"] == {name:sha(CHECKPOINTS / name) for name in paths}, "checkpoint bytes changed")
    require(report["actual_Rm_loop_density_inputs"] == read_saved(CHECKPOINTS / "new_source.json")
            and phase == read_saved(CHECKPOINTS / "new_common_phase_source.json")
            and result == read_saved(CHECKPOINTS / "new_loop_density_identity.json"), "actual saved packet changed")
    require(phase["same_loop_source_sha256"] == "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"
            and phase["source_checkpoint_sha256"] == sha(CHECKPOINTS / "new_source.json")
            and phase["entry_checkpoint_sha256"] == sha(CHECKPOINTS / "new_phase_operator_entered.json")
            and result["source_checkpoint_sha256"] == sha(CHECKPOINTS / "new_source.json")
            and result["entry_checkpoint_sha256"] == sha(CHECKPOINTS / "new_loop_operator_entered.json"), "execution link changed")
    require(report["scope"]["full_goal_complete"] is False
            and report["scope"]["current_Rm_RQy_germ_identity_proved"] is False, "open goal or pressure gate changed")
    return dict(passed=True, source_identity_passed=True, serialized_only=True,
        bytes=len(compressed), sha256=receipt["archive_sha256"],
        decoded_bytes=len(decoded), decoded_sha256=receipt["decoded_sha256"],
        proved_rows=len(result["rows"]), new_operator_replayed=False,
        old_mathematics_executed=False, full_goal_complete=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pack", action="store_true")
    group.add_argument("--inspect", action="store_true")
    args = parser.parse_args()
    print(json.dumps(pack() if args.pack else inspect()))
