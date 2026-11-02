"""Freeze the source-bound row-5 affine forcing adapter and its inputs.

The freezer only records the exact source manifest and advances the brief.
It does not import the new leaf, read a receipt payload, or invoke a
mathematical/source consumer.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
PREVIOUS = "inlet_kchi_tail_safe_v1"
NAME = "row5_affine_forcing_v1"
DIR = ROOT / "tmp" / NAME
PREIMAGE = Path("C:/Users/z5242/Documents/Codex/row5_affine_pre.json")

LEAF = "src/openai_ns_reconstruction/row5_affine_forcing_v1.py"
FIXTURE = "src/openai_ns_reconstruction/row5_affine_forcing_fixture_v1.py"
RUNNER = "scripts/run_row5_affine_forcing_v1_once.py"
G5 = "experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz"
KCHI = "experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz"
BIND = "experiments/root_st073/n1_collar_axis_potential_binding_v1_receipt.json.gz"
SELECT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json"

FILES = (LEAF, FIXTURE, RUNNER, G5, KCHI, BIND, SELECT)
PIN = {
    G5: "96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2",
    KCHI: "b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0",
    BIND: "5859d5cb4bb82d372d8ac2f3ea1629e07c9787ffbaca149c99d1aba5edb8c6e8",
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
    assert brief["component"] == PREVIOUS
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
    brief["source_checkpoint"] = {
        "path": str(DIR / "source.json"),
        "sha256": sha(source),
        "preimage_path": str(PREIMAGE),
        "preimage_sha256": sha(brief_raw),
    }
    brief["next_task"] = "Run the frozen row-5 affine forcing adapter once."
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
