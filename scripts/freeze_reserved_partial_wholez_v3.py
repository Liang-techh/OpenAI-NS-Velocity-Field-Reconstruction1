"""Freeze the reserved partial whole-Z consumer and its source contract.

This freezer only binds immutable input bytes and advances the shared brief.
The root agent runs it after the new leaf and fixture have been published.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
PREVIOUS = "reserved_partial_wholez_v2"
NAME = "reserved_partial_wholez_v3"
DIRECTORY = ROOT / "tmp" / NAME
PREIMAGE = Path("C:/Users/z5242/Documents/Codex/reserved_partial_wholez_v3_pre.json")

LEAF = "src/openai_ns_reconstruction/reserved_partial_wholez_v3.py"
FIXTURE = "src/openai_ns_reconstruction/reserved_partial_wholez_fixture_v3.py"
RUNNER = "scripts/run_reserved_partial_wholez_v3_once.py"
REPAIRED = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_repaired_Rp_m15_numeric.json"
)
COEFF = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_"
    "reserved_coefficient_enclosures.json.gz"
)
COEFF_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_"
    "reserved_coefficient_enclosures_check.json"
)
FULL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_"
    "full_five_numeric_bridge.json.gz"
)
FULL_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_"
    "full_five_numeric_bridge_check.json"
)
LINEAR = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_"
    "reserved_linear_operator.json.gz"
)
LINEAR_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_"
    "reserved_linear_operator_check.json"
)
FAILURE = "experiments/root_st073/reserved_partial_wholez_v1_failure.json.gz"
STARTUP_FAILURE = "experiments/root_st073/reserved_partial_wholez_v2_startup_failure.json.gz"
JOINT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_"
    "joint_reserved_m15_numeric.json"
)

FILES = (
    LEAF,
    FIXTURE,
    RUNNER,
    REPAIRED,
    COEFF,
    COEFF_CHECK,
    FULL,
    FULL_CHECK,
    LINEAR,
    LINEAR_CHECK,
    JOINT,
    FAILURE,
    STARTUP_FAILURE,
)
RAW_ARCHIVES = frozenset((COEFF, FULL, LINEAR, FAILURE, STARTUP_FAILURE))
PIN = {
    FAILURE: "a1dc1f5ea4c1e613ac86a2061aea9ecbe1ca8ade8cf9c1996247d8ac179ef2ca",
    STARTUP_FAILURE: "47be62777ad2b4d11317891cd1a3ecccb277622e0015ad372c74c297018b92be",
    REPAIRED: "bd9fb2fc23072c743eeeadd300041ba187aff908a3bff4a3399f7222b72e6871",
    COEFF: "c4f9f43c52431195338cc1b06d6b399585ee47abefcda826806f51fe27600392",
    COEFF_CHECK: "a005a81c70a14e883b6bb365773027339bef09896398883929788cd1096072e2",
    FULL: "3dfebf4f49f7dd73f63c8418e0b18ff3e30c93a5860b4684d84df0620855a842",
    FULL_CHECK: "e00980861fbcfaca3bad0c8c0237ad0e1afb2ecd1df3163e91c6a7f17601dcea",
    LINEAR: "98fd45a48d72656ec91f19a5564577baed0dbb43fe85ce8c842336856dc17a17",
    LINEAR_CHECK: "b8fe16b23af26a51a3601c5c0dfdc888591d1a6303c74a8ba5d53da7c7b7827c",
    JOINT: "a46413fd755aa96c385aa7b81e1919fe0d22838868b8b5a7d4348708dd16dad6",
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def encode(value: object) -> bytes:
    return (
        json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        + "\n"
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
    assert len(brief_raw) <= 97000
    assert not DIRECTORY.exists() and not PREIMAGE.exists()

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
                "raw_sha256_before_decode": path in RAW_ARCHIVES,
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
                "leaf_api": "build(c, *, response, cells, mu, joint, pulse)",
                "fixture_api": "run_fixture(leaf, c)",
                "interval_precision_dps": 130,
                "interval_endpoint_policy": (
                    "exact mpmath IV endpoint tuples; construct with c.make_mpf; "
                    "no floating-point approximations"
                ),
                "partial_row": {
                    "path": REPAIRED,
                    "x": "1.27",
                    "angular_inverse_response_count": 3,
                    "response_has_explicit_Z": False,
                    "mu_independence_is_prior_contract": True,
                },
                "coefficient_whole_Z": {
                    "path": COEFF,
                    "cells_field": "samples.whole_Z.cells",
                    "clipped_Z_required": True,
                    "normalized_theta_inputs_orders": ["Z0", "Z1"],
                    "preserve_each_raw_cell_record": True,
                    "coverage_domain": "(-1,1)",
                    "boundary_equal_is_allowed": True,
                    "no_function_point_materialization": True,
                    "no_live_owner_claim": True,
                },
                "linear_mu": {
                    "path": LINEAR,
                    "field": "actual_normalized_matrix_enclosures.mu",
                },
                "full_pressure": {
                    "path": FULL,
                    "field": "samples.whole_Z.pressure_component_enclosures.m15.pulse",
                    "orders": ["Z0", "Z1"],
                    "plain_list": True,
                    "positive": True,
                    "P1_participates": False,
                },
                "joint": {
                    "path": JOINT,
                    "field": "joint_reserved_minus_target_source.coefficient_C0_C1",
                    "meaning": "J=-N at Z0/Z1",
                    "full_source_hash_must_match": FULL,
                },
                "source_family_and_check_gates": "same gates and pins as prior reserved-m15 runner",
            },
            "source_policy": {
                "source_bind_after_freeze": True,
                "authenticate_all_input_bytes_before_import_or_decode": True,
                "authenticate_raw_gzip_before_decode": True,
                "fixture_before_saved_receipt_decode": True,
                "no_ancestor_import_or_reconstruction": True,
                "exclusive_entry_and_result": True,
                "actual_leaf_build_calls": 1,
                "complete_result_and_failures_are_saved": True,
                "no_duplicate_ledger": True,
            },
        }
    )
    assert len(source) <= 97000

    brief["component"] = NAME
    brief["phase"] = NAME + "_pending"
    brief["one_entry_only"] = True
    brief["do_not_reenter"] = True
    brief["source_checkpoint"] = {
        "path": str(DIRECTORY / "source.json"),
        "sha256": sha(source),
        "preimage_path": str(PREIMAGE),
        "preimage_sha256": sha(brief_raw),
    }
    brief["next_task"] = "Run the frozen reserved partial whole-Z consumer once."
    pending = encode(brief)
    assert len(pending) <= 97000

    DIRECTORY.mkdir(parents=True)
    with PREIMAGE.open("xb") as handle:
        handle.write(brief_raw)
    with (DIRECTORY / "source.json").open("xb") as handle:
        handle.write(source)
    BRIEF.write_bytes(pending)
    print(json.dumps({"source_sha256": sha(source), "brief_bytes": len(pending)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
