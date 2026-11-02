"""Freeze the source-only W0 ordinary-Z2 four-cell consumer.

This freezer records the exact new leaf/fixture/runner bytes together with
the four already-authenticated source records.  It only advances the pending
brief; it does not import a leaf, decode a receipt, or invoke a source
consumer.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
PREVIOUS = "cartesian_z_pullback_v1"
NAME = "w0_z2_four_cell_v1"
DIR = ROOT / "tmp" / NAME
PREIMAGE = Path("C:/Users/z5242/Documents/Codex/w0_z2_fourcell_pre.json")

LEAF = "src/openai_ns_reconstruction/w0_z2_four_cell_v1.py"
FIXTURE = "src/openai_ns_reconstruction/w0_z2_four_cell_fixture_v1.py"
RUNNER = "scripts/run_w0_z2_four_cell_v1_once.py"
GREEN = "experiments/root_st073/n1_actual_collar_green_v3_receipt.json.gz"
FORCING = "experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz"
W0_Z0Z1 = "experiments/root_st073/n1_w0_collar_sectors_v1_receipt.json.gz"
SELECTED = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json"

FILES = (LEAF, FIXTURE, RUNNER, GREEN, FORCING, W0_Z0Z1, SELECTED)
PIN = {
    GREEN: "fad35594ccf5ef81392f94f97ae72b3ffb94fcaa554438fd5702ef2c4343c4b5",
    FORCING: "96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2",
    W0_Z0Z1: "4229c29a2a21207f8b8da3ba53d70bbf6f50e09eb67110bd3b7cb30fa6577c06",
    SELECTED: "0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e",
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def encode(value: object) -> bytes:
    return (json.dumps(value, ensure_ascii=False, separators=(",", ":")) + "\n").encode("utf-8")


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
        digest = sha_path(ROOT / path)
        if path in PIN:
            assert digest == PIN[path], f"pinned input changed: {path}"
        rows.append({"path": path, "sha256": digest, "raw_sha256_only": path in {GREEN, FORCING, W0_Z0Z1}})
    assert len({row["path"] for row in rows}) == len(FILES)

    source = encode(
        {
            "component": NAME,
            "one_entry_only": True,
            "previous_component": PREVIOUS,
            "purpose": "source-only W0 ordinary-Z2 four-cell sectors",
            "source_policy": {
                "source_bind_after_freeze": True,
                "authenticate_all_raw_bytes_before_decode": True,
                "previous_receipts_decode_only_after_entry": True,
                "no_ancestor_import_or_reconstruction": True,
                "old_W0_Z0_Z1_and_all_word_errors_remain_separate": True,
            },
            "actual_input_paths": {
                "green_core": ["result", "output", "result", "core_W0_transport", "2"],
                "forcing_g_rows": ["result", "output", "g", "2", "5|6"],
                "saved_W0_Z0_Z1": "result.output",
                "selected_support_measure": "actual_cutoff_complete_support_measure",
            },
            "ordinary_Z_contract": {
                "actual_core_order": 2,
                "actual_forcing_order": 2,
                "old_W0_receipt_orders": [0, 1],
                "new_q_cell_count": 4,
            },
            "files": rows,
        }
    )

    brief["component"] = NAME
    brief["phase"] = NAME + "_pending"
    brief["source_checkpoint"] = {
        "path": str(DIR / "source.json"),
        "sha256": sha(source),
        "preimage_path": str(PREIMAGE),
        "preimage_sha256": sha(brief_raw),
    }
    brief["next_task"] = "Run the frozen W0 ordinary-Z2 four-cell source consumer once."
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
