"""Freeze the repaired local-g5 word-2 consumer and its sealed inputs.

This freezer consumes the failed v1 archive as evidence only.  It never
re-enters the v1 leaf, fixture, runner, or result.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
NAME = "collar_g5_word2_v2"
PREVIOUS = "collar_g5_word2_v1"
PREVIOUS_PHASE = "collar_g5_word2_v1_pending"
PREVIOUS_SOURCE_SHA = "8685a61c741feb832fdc0d6a3eece089f82f5eff642488bbbaa59ca0a8705f27"
DIR = ROOT / "tmp" / NAME
PREIMAGE = Path("C:/Users/z5242/Documents/Codex/g5_word2_v2_pre.json")

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

FILES = (LEAF, FIXTURE, RUNNER, G5, KCHI, BIND, FIRST, OLD_TAIL, PRESS, SELECT, FAILURE)
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


def encode(value: object) -> bytes:
    return (json.dumps(value, ensure_ascii=False, separators=(",", ":")) + "\n").encode("utf-8")


def main() -> int:
    if not __debug__:
        raise RuntimeError("optimized mode disables required preflight checks")

    brief_raw = BRIEF.read_bytes()
    brief = json.loads(brief_raw)
    assert brief["component"] == PREVIOUS and brief["phase"] == PREVIOUS_PHASE
    assert brief["source_checkpoint"]["sha256"] == PREVIOUS_SOURCE_SHA
    assert brief["goal_active"] is True and brief["full_goal_complete"] is False
    assert brief["one_entry_only"] is True and brief["do_not_reenter"] is True
    assert not DIR.exists() and not PREIMAGE.exists()
    assert LEAF != V1_LEAF and FIXTURE != V1_FIXTURE and RUNNER != V1_RUNNER

    failure_raw = (ROOT / FAILURE).read_bytes()
    assert sha(failure_raw) == FAILURE_SHA
    rows = []
    for path in FILES:
        digest = sha((ROOT / path).read_bytes())
        if path in PIN:
            assert digest == PIN[path], f"pinned input changed: {path}"
        rows.append({"path": path, "sha256": digest})

    assert len({row["path"] for row in rows}) == len(FILES)
    assert FAILURE in {row["path"] for row in rows}
    source = encode({"component": NAME, "one_entry_only": True, "files": rows})
    brief["component"] = NAME
    brief["phase"] = NAME + "_pending"
    brief["source_checkpoint"] = {
        "path": str(DIR / "source.json"),
        "sha256": sha(source),
        "preimage_path": str(PREIMAGE),
        "preimage_sha256": sha(brief_raw),
    }
    brief["next_task"] = "Run the repaired local-g5 signed word-2 and alternative j>=3 tail once."
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
