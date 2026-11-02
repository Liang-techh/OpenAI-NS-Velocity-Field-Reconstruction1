"""Freeze the source-bound Z2 seam bridge.

The freezer records the exact bytes of the new bridge, its synthetic fixture,
the one-entry runner, and the two already-published compact receipts.  It
only advances the resume brief to the pending state.  It does not import the
new binder, decode a receipt, or invoke an ancestor consumer.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
PREVIOUS = "core_firstword_z2_profile_v2"
NAME = "core_z2_seam_bridge_v2"
DIR = ROOT / "tmp" / NAME
PREIMAGE = Path("C:/Users/z5242/Documents/Codex/core_z2_seam_bridge_v2_pre.json")

LEAF = "src/openai_ns_reconstruction/core_z2_seam_bridge_v2.py"
FIXTURE = "src/openai_ns_reconstruction/core_z2_seam_bridge_fixture_v2.py"
RUNNER = "scripts/run_core_z2_seam_bridge_v2_once.py"
OLD = "experiments/root_st073/core_secondword_partial_v2_receipt.json.gz"
NEW = "experiments/root_st073/core_firstword_z2_profile_v2_receipt.json.gz"
LOW = "experiments/root_st073/n1_core_low_word_v3_receipt.json.gz"

INCIDENT = "experiments/root_st073/core_z2_seam_bridge_v1_uncheckpointed_incident.json"
SUPPLEMENT = "experiments/root_st073/core_z2_seam_bridge_v1_uncheckpointed_supplement.md"
FILES = (LEAF, FIXTURE, RUNNER, OLD, NEW, LOW, INCIDENT, SUPPLEMENT)
RECEIPTS = {OLD, NEW, LOW}
PIN = {
    OLD: "d7de2dbcd1fbc90bf039d0b8527cc82be6e7826e475c2c0c1fdfa6ed3c8585fe",
    NEW: "bcbd81914ba8294a8da13790640d63ab12907fd41624fceff10f304a5ca12115",
    LOW: "5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6",
}

LOW_SHA256 = "5a5fe7c57a2efe10e58089ced608e22937c54fcf54c5f224a42abe54e450c6d6"
Z2 = "experiments/root_st073/n1_core_word1_z2_v2_receipt.json.gz"
Z2_SHA256 = "c53ff32ad3ce5288279953e84bc7417b71477c82b00d94466d0d2880adbc7ac9"


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
    assert brief["phase"] == "component_published"
    assert brief["goal_active"] is True and brief["full_goal_complete"] is False
    assert brief["one_entry_only"] is True and brief["do_not_reenter"] is True
    assert not DIR.exists() and not PREIMAGE.exists()

    rows = []
    for path in FILES:
        raw = (ROOT / path).read_bytes()
        digest = sha(raw)
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
            "previous_source_checkpoint": brief.get("source_checkpoint"),
            "files": rows,
            "input_contract": {
                "old_packet": "result.output.result",
                "new_profile_packet": "result.output.result",
                "new_seam_packet": "result.output.missing_seam_result",
                "low_output": "result.output",
                "leaf_api": "build(old, new_profile, new_seam, low_output)",
                "fixture_api": "run_fixture(leaf)",
                "exact_endpoint_policy": (
                    "inputs are already exact MPF tuple dicts; preserve and compare "
                    "lower_exact_mpf_tuple/upper_exact_mpf_tuple without mpmath materialization"
                ),
                "lineage_policy": {
                    "both_source_manifests_pin": {"path": LOW, "sha256": LOW_SHA256},
                    "new_source_manifest_pins": {"path": Z2, "sha256": Z2_SHA256},
                    "saved_new_validation_pins": [OLD, LOW, Z2],
                    "physical_owner_claim": False,
                },
            },
            "source_policy": {
                "source_bind_after_freeze": True,
                "authenticate_all_raw_gzip_bytes_before_decode": True,
                "previous_receipts_decode_only_after_entry": True,
                "no_ancestor_import_or_reconstruction": True,
                "fixture_before_saved_receipt_consumption": True,
                "complete_result_and_failures_are_saved": True,
                "retain_all_user_branches": True,
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
    brief["next_task"] = "Run the frozen core Z2 seam bridge once."
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
