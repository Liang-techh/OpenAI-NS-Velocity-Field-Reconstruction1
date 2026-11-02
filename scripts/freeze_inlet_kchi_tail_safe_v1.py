"""Freeze the safe corrected formal inlet Kchi tail consumer."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
NAME = "inlet_kchi_tail_safe_v1"
DIR = ROOT / "tmp" / NAME
PREIMAGE = Path("C:/Users/z5242/Documents/Codex/inlet_tail_safe_pre.json")

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

FILES = (LEAF, FIXTURE, RUNNER, FORMAL, BASE, KCHI, BIND, FIRST, OLDTAIL, SELECT)

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


def encode(value: object) -> bytes:
    return (json.dumps(value, ensure_ascii=False, separators=(",", ":")) + "\n").encode("utf-8")


def main() -> int:
    if not __debug__:
        raise RuntimeError("optimized mode disables required preflight checks")

    brief_raw = BRIEF.read_bytes()
    brief = json.loads(brief_raw)
    assert brief["component"] == "collar_g5_word2_v2"
    assert brief["phase"] == "component_published"
    assert brief["goal_active"] is True and brief["full_goal_complete"] is False
    assert brief["one_entry_only"] is True and brief["do_not_reenter"] is True
    assert not DIR.exists() and not PREIMAGE.exists()

    rows = []
    for path in FILES:
        digest = sha((ROOT / path).read_bytes())
        if path in PIN:
            assert digest == PIN[path], f"pinned input changed: {path}"
        rows.append({"path": path, "sha256": digest})
    assert len({row["path"] for row in rows}) == len(FILES)

    source = encode({"component": NAME, "one_entry_only": True, "files": rows})
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
    brief["next_task"] = "Run the safe corrected formal inlet Kchi tail once."
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
