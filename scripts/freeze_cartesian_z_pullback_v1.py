"""Freeze the Cartesian fixed-Z pullback before its single source entry.

This freezer records source bytes and advances the pending brief only.  It
does not import the new leaf or fixture and it never decodes the large prior
row-5 receipt; that receipt is authenticated by its raw SHA-256 only.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
PREVIOUS = "row5_firstword_product_v4"
NAME = "cartesian_z_pullback_v1"
DIR = ROOT / "tmp" / NAME
PREIMAGE = Path("C:/Users/z5242/Documents/Codex/cartesian_z_pre.json")

LEAF = "src/openai_ns_reconstruction/cartesian_z_pullback_v1.py"
FIXTURE = "src/openai_ns_reconstruction/cartesian_z_pullback_fixture_v1.py"
RUNNER = "scripts/run_cartesian_z_pullback_v1_once.py"
CARTESIAN_SOURCE = "src/openai_ns_reconstruction/original_n1_cartesian_core_source.py"
PHYSICAL_SOURCE = "src/openai_ns_reconstruction/original_n1_physical_source.py"
PREVIOUS_RECEIPT = "experiments/root_st073/row5_firstword_product_v4_receipt.json.gz"

FILES = (
    LEAF,
    FIXTURE,
    RUNNER,
    CARTESIAN_SOURCE,
    PHYSICAL_SOURCE,
    PREVIOUS_RECEIPT,
)

PIN = {
    CARTESIAN_SOURCE: "bc65a55d79df9812c404a856c8afd0d3a250c0ec6511ea87f50dc80f4546199e",
    PHYSICAL_SOURCE: "542178f6a71dd20e7091dc5f76f115c78274717fd4cce5be391e4c2dbd1158bb",
    PREVIOUS_RECEIPT: "1db10d5abf6b619297d1dd0be282e7dad845b9b5fdf0ecdff000b15693fccfc1",
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
        rows.append(
            {
                "path": path,
                "sha256": digest,
                "raw_sha256_only": path == PREVIOUS_RECEIPT,
            }
        )

    source = encode(
        {
            "component": NAME,
            "one_entry_only": True,
            "previous_component": PREVIOUS,
            "coordinate_formula": {
                "lambda_squared": "tau/(1-Z^2)",
                "R": "(X^2+Y^2)/(2*lambda^2)",
                "rho": "R/epsilon",
                "z_phys": "lambda^(1-delta)*Z",
                "row5_current_x": "sqrt(rho), normalized; not Cartesian X",
            },
            "source_policy": {
                "source_bind_after_freeze": True,
                "static_coordinate_labels_only": True,
                "previous_receipt_decode_forbidden": True,
                "old_leaf_or_calculation_forbidden": True,
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
    brief["next_task"] = "Run the frozen Cartesian fixed-Z pullback leaf and fixture once."
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
