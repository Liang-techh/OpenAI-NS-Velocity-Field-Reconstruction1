"""Freeze the independent first-word ordinary-Z2 core profile consumer.

The freezer only authenticates source bytes, captures the published
``core_secondword_partial_v2`` brief as an independent preimage, and advances
the brief to the one-entry pending state.  It never imports a leaf or fixture,
decodes a receipt, or invokes an ancestor consumer.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
PREVIOUS = "core_firstword_z2_profile_v1"
NAME = "core_firstword_z2_profile_v2"
DIR = ROOT / "tmp" / NAME
PREIMAGE = Path("C:/Users/z5242/Documents/Codex/core_firstword_z2_profile_v2_pre.json")

LEAF = "src/openai_ns_reconstruction/core_firstword_z2_profile_v2.py"
FIXTURE = "src/openai_ns_reconstruction/core_firstword_z2_profile_fixture_v2.py"
RUNNER = "scripts/run_core_firstword_z2_profile_v2_once.py"
Z2 = "experiments/root_st073/n1_core_word1_z2_v2_receipt.json.gz"
LOW = "experiments/root_st073/n1_core_low_word_v3_receipt.json.gz"
SECOND = "experiments/root_st073/core_secondword_partial_v2_receipt.json.gz"

PRIOR = "experiments/root_st073/core_firstword_z2_profile_v1_receipt.json.gz"
FILES = (LEAF, FIXTURE, RUNNER, Z2, LOW, SECOND, PRIOR)
RECEIPTS = {Z2, LOW, SECOND, PRIOR}

# These are raw gzip-byte pins.  They are checked before any later runner
# receipt decode and are repeated in the source manifest for publication QA.
PIN = {
    PRIOR: "3ff8d17e8a5a930f1d8a0c1139be85bb61962137294dd412b3a9087cdb206642",
    Z2: "c53ff32ad3ce5288279953e84bc7417b71477c82b00d94466d0d2880adbc7ac9",
    LOW: "5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6",
    SECOND: "d7de2dbcd1fbc90bf039d0b8527cc82be6e7826e475c2c0c1fdfa6ed3c8585fe",
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def encode(value: object) -> bytes:
    return (json.dumps(value, ensure_ascii=False, separators=(",", ":")) + "\n").encode(
        "utf-8"
    )


def main() -> int:
    if not __debug__:
        raise RuntimeError("optimized mode disables required freeze guards")

    brief_raw = BRIEF.read_bytes()
    brief = json.loads(brief_raw)
    assert brief["component"] == PREVIOUS
    assert brief["phase"] == PREVIOUS + "_pending"
    assert brief["goal_active"] is True and brief["full_goal_complete"] is False
    assert brief["one_entry_only"] is True and brief["do_not_reenter"] is True
    assert not DIR.exists() and not PREIMAGE.exists()

    rows = []
    for path in FILES:
        digest = sha((ROOT / path).read_bytes())
        if path in PIN:
            assert digest == PIN[path], f"pinned input changed: {path}"
        rows.append(
            {
                "path": path,
                "sha256": digest,
                "raw_sha256_before_decode": path in RECEIPTS,
            }
        )
    assert len({row["path"] for row in rows}) == len(FILES)

    source = encode(
        {
            "component": NAME,
            "one_entry_only": True,
            "previous_component": PREVIOUS,
            "previous_source_checkpoint": brief["source_checkpoint"],
            "files": rows,
            "input_contract": {
                "z2_receipt_cells": [
                    "result.output.result.outer_cells[i].t",
                    "result.output.result.outer_cells[i].rho",
                    "result.output.result.outer_cells[i].feedback_Z2[4|5|6]",
                ],
                "leaf_api": "build(c, cells)",
                "term_policy": "retain every term field; attach input_id and exact receipt_input_path",
                "interval_policy": "rebuild directed lower/upper exact MPF tuples in the new interval context",
                "lineage_policy": "low source manifest must pin the z2 receipt; secondword receipt is the saved anchor",
                "f0_log_policy": "compare z2 F0_of_Z to low F0_of_Z or its F0_anchor alias by exact endpoints",
            },
            "source_policy": {
                "source_bind_after_freeze": True,
                "authenticate_all_raw_gzip_bytes_before_decode": True,
                "previous_receipts_decode_only_after_entry": True,
                "no_ancestor_import_or_reconstruction": True,
                "fixture_before_saved_receipt_consumption": True,
                "complete_result_and_failures_are_saved": True,
            },
        }
    )

    brief["component"] = NAME
    brief["phase"] = NAME + "_pending"
    brief["one_entry_only"] = True
    brief["do_not_reenter"] = True
    brief["source_checkpoint"] = {
        "path": str(DIR / "source.json"),
        "sha256": sha(source),
        "preimage_path": str(PREIMAGE),
        "preimage_sha256": sha(brief_raw),
    }
    brief["next_task"] = "Run the frozen independent core first-word Z2 profile consumer once."
    pending = encode(brief)
    assert len(pending) <= 97000

    DIR.mkdir(parents=True)
    with PREIMAGE.open("xb") as handle:
        handle.write(brief_raw)
    with (DIR / "source.json").open("xb") as handle:
        handle.write(source)
    BRIEF.write_bytes(pending)
    print(json.dumps({"source_sha256": sha(source), "brief_bytes": len(pending)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
