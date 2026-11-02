"""Freeze the one-entry source-AST operator formula identity proof.

The freezer records the new leaf, its synthetic fixture, this runner, the two
source files whose text is compared, and the nearby v3 receipt.  The receipt
is provenance only for this leaf: its raw gzip bytes are pinned before any
future optional decode, while the leaf itself receives source text only.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
PREVIOUS = "core_operator_formula_identity_v2"
NAME = "core_operator_formula_identity_v3"
DIR = ROOT / "tmp" / NAME
PREIMAGE = Path(
    "C:/Users/z5242/Documents/Codex/core_operator_formula_identity_v3_pre.json"
)

LEAF = "src/openai_ns_reconstruction/core_operator_formula_identity_v3.py"
FIXTURE = "src/openai_ns_reconstruction/core_operator_formula_identity_fixture_v3.py"
RUNNER = "scripts/run_core_operator_formula_identity_v3_once.py"
HIGH = "src/openai_ns_reconstruction/n1_core_low_source_matrix_v2.py"
LOW = "src/openai_ns_reconstruction/n1_core_low_matrix_v3.py"
RECEIPT = "experiments/root_st073/core_z2_seam_bridge_v3_receipt.json.gz"
INCIDENT = "experiments/root_st073/core_operator_scan_gzip_preflight_incident.json"
FAILURE = "experiments/root_st073/core_operator_formula_identity_v1_failure.json.gz"
STARTUP_FAILURE = "experiments/root_st073/core_operator_formula_identity_v2_startup_failure.json.gz"
FILES = (LEAF, FIXTURE, RUNNER, HIGH, LOW, RECEIPT, INCIDENT, FAILURE, STARTUP_FAILURE)
RECEIPTS = {RECEIPT, FAILURE, STARTUP_FAILURE}
PIN = {
    STARTUP_FAILURE: "bf411b7078beabe6148c2a1f1484dbfd5a30c1538c8a1cf39458e9b4463dc800",
    FAILURE: "3564de1416485887f1c754979368bed7f36671e2ab4bad0741b8974074c716e6",
    HIGH: "2e2a3e7e4367223f53e31e30818f53aceeab96b37c360b4e095189570ce1a111",
    LOW: "3e790126fa7ca41f3c6045004f5a98f0a4975b04681b9aa997ff7e27a69821c0",
    RECEIPT: "625fa2f92cf821301db7d02a6af92088cf9b4f229258bd4ae56bb1b35da6f2a1",
    INCIDENT: "6c99aff5729d195b33eef3ce119e13a8d233de31da4cbe501cee698706eb51f9",
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def encode(value: object) -> bytes:
    return (
        json.dumps(value, ensure_ascii=False, separators=(",", ":")) + "\n"
    ).encode("utf-8")


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
                "high_source_text": HIGH,
                "low_source_text": LOW,
                "leaf_api": "build(high_source_text, low_source_text)",
                "fixture_api": "run_fixture(leaf)",
                "source_ast_scope": "lexical source text only; no source module import",
                "associated_receipt": {
                    "path": RECEIPT,
                    "raw_sha256": PIN[RECEIPT],
                    "decode_policy": (
                        "raw gzip is authenticated before any optional decode; this leaf "
                        "does not consume receipt data"
                    ),
                },
                "preflight_incident": {
                    "path": INCIDENT,
                    "sha256": PIN[INCIDENT],
                    "is_gzip": False,
                    "purpose": "preserve the prior gzip-prefix-before-hash incident",
                },
            },
            "source_policy": {
                "source_bind_after_freeze": True,
                "authenticate_all_input_bytes_before_import_or_decode": True,
                "authenticate_raw_gzip_before_decode": True,
                "receipt_decode_optional_lineage_status_only": True,
                "receipt_decode_skipped_by_runner": True,
                "no_ancestor_import_or_reconstruction": True,
                "fixture_before_saved_source_consumption": True,
                "actual_leaf_build_calls": 1,
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
    brief["next_task"] = "Run the frozen core operator formula identity proof once."
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
