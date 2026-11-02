#!/usr/bin/env python3
"""Freeze the row5 first-word product component's authenticated preimage.

This is intentionally a source-only freezer.  It does not import, execute, or
decompress any receipt or mathematical code.  The one-entry guard is checked
before the preimage is created, and the brief is updated only after the
preimage and its source manifest have been written exclusively.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
PREVIOUS = "row5_affine_forcing_v1"
NAME = "row5_firstword_product_v1"
DIR = ROOT / "tmp" / NAME
PREIMAGE = Path("C:/Users/z5242/Documents/Codex/row5_firstword_pre.json")

LEAF = "src/openai_ns_reconstruction/row5_firstword_product_v1.py"
FIXTURE = "src/openai_ns_reconstruction/row5_firstword_product_fixture_v1.py"
RUNNER = "scripts/run_row5_firstword_product_v1_once.py"
AFFINE = "experiments/root_st073/row5_affine_forcing_v1_receipt.json.gz"
JOIN = "experiments/root_st073/n1_core_collar_low_join_v2_receipt.json.gz"
GREEN = "experiments/root_st073/n1_actual_collar_green_v3_receipt.json.gz"
FORCING = "experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz"

FILES = (LEAF, FIXTURE, RUNNER, AFFINE, JOIN, GREEN, FORCING)
PIN = {
    AFFINE: "b64d14b59b178502e562fe593d79419cc2e21200f80c9bf0006905d1b6e72515",
    JOIN: "a2e25ccec22f97dc99d1eb80609f0d806d4276e7a75e763ca882d90532e7ce67",
    GREEN: "fad35594ccf5ef81392f94f97ae72b3ffb94fcaa554438fd5702ef2c4343c4b5",
    FORCING: "96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2",
}


def sha(raw: bytes) -> str:
    h = hashlib.sha256()
    h.update(raw)
    return h.hexdigest()


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
        if path in PIN and digest != PIN[path]:
            raise AssertionError(f"pinned input changed: {path}")
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
    brief["next_task"] = "Run the frozen row-5 first-word product consumer once."
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
