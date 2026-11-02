"""Run the frozen first-collar finite-prefix sector grouping exactly once.

Only the published pressure receipt is consumed.  Its raw gzip bytes are
authenticated before decoding; the synthetic fixture runs after the raw
check and before gzip decoding.  The receipt already contains physically
scaled axial, swirl, and pressure channels, so this consumer performs only
ordinary-Z sector grouping and outward interval addition.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
COMPONENT = "first_collar_m15_group_v1"
DIR = ROOT / "tmp" / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / name for name in ("source.json", "entry.json", "result.json"))

LEAF = "src/openai_ns_reconstruction/first_collar_m15_sector_group.py"
FIXTURE = "src/openai_ns_reconstruction/first_collar_m15_group_fixture_v1.py"
RUNNER = "scripts/run_first_collar_m15_group_v1_once.py"
PRIOR_RECEIPT = "experiments/root_st073/n1_first_collar_g4_pressure_v1_receipt.json.gz"
RECEIPT_SHA256 = "c229cf188d306413f6edd0641a729321ad4567c2587ebf4701ada47720853660"
DOCUMENTED_RESULT_SHA256 = "cce81dc122acba8a0d0e47c40f10aa1fa45f73e2bae61db9114575a23fabaff9"
REQUIRED = (LEAF, FIXTURE, RUNNER, PRIOR_RECEIPT)
PIN = {PRIOR_RECEIPT: RECEIPT_SHA256}
CHANNELS = ("axial", "swirl", "pressure")
EXPECTED_COUNTS = {
    0: {"axial": 520, "swirl": 283, "pressure": 43},
    1: {"axial": 520, "swirl": 566, "pressure": 43},
}
FORBIDDEN_BASENAMES = {
    "resume_first_switch_2026_10_17.json",
    "rm_loop_density_identity.json",
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(value: object, message: str) -> None:
    if not value:
        raise ValueError(message)


def save(path: Path, value: object) -> str:
    raw = (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"
    ).encode("utf-8")
    with path.open("xb") as handle:
        handle.write(raw)
    return sha(raw)


def reject_float_tree(value: Any, label: str) -> None:
    if isinstance(value, float):
        raise TypeError(label + " must not contain floating-point values")
    if isinstance(value, Mapping):
        for key, item in value.items():
            reject_float_tree(key, label + ".key")
            reject_float_tree(item, label)
    elif isinstance(value, (list, tuple, set, frozenset)):
        for item in value:
            reject_float_tree(item, label)


def encode(value: Any) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return {
            "lower_exact_mpf_tuple": list(raw[0]),
            "upper_exact_mpf_tuple": list(raw[1]),
            "directed": True,
        }
    if isinstance(value, Mapping):
        return {str(key): encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return [encode(item) for item in sorted(value, key=repr)]
    if isinstance(value, float):
        raise TypeError("cannot serialize floating-point value")
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    raise TypeError("cannot serialize result type " + type(value).__name__)


def interval(context: MPIntervalContext, value: Any, label: str) -> Any:
    if not isinstance(value, Mapping):
        raise TypeError(label + " must be an encoded directed interval")
    reject_float_tree(value, label)
    require(value.get("directed") is True, label + " is not directed")
    lower = value.get("lower_exact_mpf_tuple")
    upper = value.get("upper_exact_mpf_tuple")
    require(
        isinstance(lower, (list, tuple)) and isinstance(upper, (list, tuple)),
        label + " endpoint tuples missing",
    )
    require(len(lower) == len(upper) == 4, label + " endpoint tuple length")
    require(
        all(isinstance(item, int) and not isinstance(item, bool) for item in lower + upper),
        label + " endpoint tuple is not integral",
    )
    result = context.make_mpf((tuple(lower), tuple(upper)))
    require(result.a <= result.b, label + " reversed interval")
    return result


def module(path: str, name: str) -> Any:
    require(path in (LEAF, FIXTURE), "unexpected import path: " + path)
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    require(spec is not None and spec.loader is not None, "new component unavailable: " + path)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def _forbidden(path: Path) -> bool:
    return path.name.lower() in FORBIDDEN_BASENAMES


def preflight() -> tuple[dict[str, str], str, str]:
    require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    brief_raw = BRIEF.read_bytes()
    source_raw = SOURCE.read_bytes()
    brief = json.loads(brief_raw)
    manifest = json.loads(source_raw)
    require(
        brief.get("component") == COMPONENT
        and brief.get("phase") == COMPONENT + "_pending"
        and brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("one_entry_only") is True
        and brief.get("do_not_reenter") is True,
        "brief guard",
    )
    checkpoint = brief.get("source_checkpoint")
    require(isinstance(checkpoint, Mapping), "source checkpoint missing")
    require(Path(checkpoint["path"]).resolve() == SOURCE.resolve(), "source checkpoint path")
    require(checkpoint["sha256"] == sha(source_raw), "source checkpoint hash")
    preimage_path = Path(checkpoint["preimage_path"])
    require(not _forbidden(preimage_path), "forbidden preimage path")
    require(
        sha_path(preimage_path) == checkpoint["preimage_sha256"],
        "pending brief preimage hash",
    )
    require(
        manifest.get("component") == COMPONENT and manifest.get("one_entry_only") is True,
        "manifest guard",
    )
    require(manifest.get("previous_component") == "reserved_m15_variable_z", "predecessor guard")
    rows = manifest.get("files")
    require(isinstance(rows, list) and len(rows) == len(REQUIRED), "manifest file count")
    hashes: dict[str, str] = {}
    for row in rows:
        require(isinstance(row, Mapping), "manifest row")
        path = row.get("path")
        require(isinstance(path, str) and path in REQUIRED and path not in hashes, "unexpected manifest input")
        require(not _forbidden(Path(path)), "forbidden manifest path")
        digest = sha_path(ROOT / path)
        require(digest == row.get("sha256"), "changed input: " + path)
        if path in PIN:
            require(digest == PIN[path], "changed pinned receipt: " + path)
        hashes[path] = digest
    require(set(hashes) == set(REQUIRED), "missing manifest input")
    require(hashes[PRIOR_RECEIPT] == RECEIPT_SHA256, "prior fixed receipt hash")
    return hashes, sha(source_raw), sha(brief_raw)


def _hydrate_record(
    context: MPIntervalContext,
    raw_record: Any,
    label: str,
) -> dict[str, Any]:
    require(isinstance(raw_record, Mapping), label + " must be a record mapping")
    reject_float_tree(raw_record, label)
    origin = raw_record.get("origin")
    require(isinstance(origin, Mapping), label + ".origin missing")
    source_ids = origin.get("source_ids")
    error_classes = origin.get("error_classes")
    require(
        isinstance(source_ids, Sequence) and not isinstance(source_ids, (str, bytes, bytearray)) and source_ids,
        label + ".origin.source_ids missing",
    )
    require(
        isinstance(origin.get("uncertainty"), bool),
        label + ".origin.uncertainty must be bool",
    )
    require(
        isinstance(error_classes, Sequence) and not isinstance(error_classes, (str, bytes, bytearray)),
        label + ".origin.error_classes must be ordered sequence",
    )
    result = dict(raw_record)
    result["coefficient"] = interval(context, raw_record.get("coefficient"), label + ".coefficient")
    return result


def _hydrate_channels(
    context: MPIntervalContext,
    finite_prefix_channels: Any,
) -> tuple[dict[str, dict[int, list[dict[str, Any]]]], dict[str, Any]]:
    require(isinstance(finite_prefix_channels, Mapping), "finite_prefix_channels mapping")
    reject_float_tree(finite_prefix_channels, "finite_prefix_channels")
    require(set(finite_prefix_channels) == {"0", "1"}, "finite_prefix_channels must be order-major Z0/Z1")
    channels: dict[str, dict[int, list[dict[str, Any]]]] = {
        channel: {0: [], 1: []} for channel in CHANNELS
    }
    input_index: dict[str, Any] = {}
    for order in (0, 1):
        row = finite_prefix_channels[str(order)]
        require(isinstance(row, Mapping), "finite prefix order row")
        require(set(row) == set(CHANNELS), "finite prefix channel set")
        for channel in CHANNELS:
            records = row[channel]
            require(isinstance(records, list), "finite prefix records must be lists")
            require(len(records) == EXPECTED_COUNTS[order][channel], "unexpected " + channel + " Z" + str(order) + " count")
            for index, raw_record in enumerate(records):
                label = "finite_prefix_channels[" + str(order) + "][" + channel + "][" + str(index) + "]"
                record = _hydrate_record(context, raw_record, label)
                channels[channel][order].append(record)
                origin = record["origin"]
                path = ["result", "output", "finite_prefix_channels", str(order), channel, index]
                for source_id in origin["source_ids"]:
                    require(isinstance(source_id, str) and source_id, label + " source id must be nonempty string")
                    occurrence = {
                        "receipt": PRIOR_RECEIPT,
                        "path": path,
                        "channel": channel,
                        "ordinary_Z_order": order,
                        "record_index": index,
                        "uncertainty": origin["uncertainty"],
                        "error_classes": list(origin["error_classes"]),
                        "origin": encode(origin),
                    }
                    bucket = input_index.setdefault(
                        source_id,
                        {"source_id": source_id, "occurrences": []},
                    )
                    bucket["occurrences"].append(occurrence)
    return channels, input_index


def _validate_grouped(
    grouped: Mapping[str, Any],
    channels: Mapping[str, Mapping[int, Sequence[Mapping[str, Any]]]],
) -> dict[str, int]:
    groups = grouped.get("groups")
    require(isinstance(groups, Mapping), "leaf groups missing")
    expected: dict[tuple[str, int, int], Mapping[str, Any]] = {}
    for channel, rows in channels.items():
        for order in (0, 1):
            for index, record in enumerate(rows[order]):
                expected[(channel, order, index)] = record
    seen: set[tuple[str, int, int]] = set()
    counts: dict[str, int] = {}
    for order in (0, 1):
        row = groups.get(order, groups.get(str(order)))
        require(isinstance(row, list), "leaf groups missing Z" + str(order))
        counts[str(order)] = len(row)
        for group in row:
            require(isinstance(group, Mapping), "leaf group mapping")
            members = group.get("members")
            origin = group.get("origin")
            require(isinstance(members, list) and isinstance(origin, Mapping), "leaf group provenance")
            require(origin.get("member_count") == len(members), "leaf member count")
            require(isinstance(origin.get("error_classes"), list), "leaf error class order")
            first_classes: list[Any] | None = None
            for member in members:
                require(isinstance(member, Mapping), "leaf member mapping")
                key = (member.get("channel"), member.get("ordinary_Z_order"), member.get("index"))
                require(key in expected and key not in seen, "leaf member lost or duplicated")
                seen.add(key)
                record = expected[key]
                require(member.get("record") is record, "leaf member record not retained")
                provenance = member.get("provenance")
                require(isinstance(provenance, Mapping), "leaf member provenance missing")
                require(provenance.get("source_ids") == record["origin"]["source_ids"], "source ids changed")
                require(provenance.get("uncertainty") is record["origin"]["uncertainty"], "uncertainty changed")
                require(provenance.get("error_classes") == record["origin"]["error_classes"], "error class order changed")
                require(provenance == record["origin"], "member origin changed")
                if first_classes is None:
                    first_classes = list(record["origin"]["error_classes"])
            require(first_classes == origin.get("error_classes"), "group error class order changed")
    require(seen == set(expected), "not every input record retained")
    return counts


def consume(hashes: Mapping[str, str], progress: dict[str, Any]) -> dict[str, Any]:
    # Read and authenticate every byte first.  No gzip decode occurs here.
    raw: dict[str, bytes] = {}
    for path in REQUIRED:
        data = (ROOT / path).read_bytes()
        require(sha(data) == hashes[path], "changed after entry: " + path)
        raw[path] = data
    require(sha(raw[PRIOR_RECEIPT]) == RECEIPT_SHA256, "prior receipt raw hash")

    # The fixture is deliberately before receipt decoding.
    leaf = module(LEAF, COMPONENT + "_leaf")
    fixture = module(FIXTURE, COMPONENT + "_fixture")
    context = MPIntervalContext()
    context.prec = 1800
    proof = fixture.run_fixture(leaf, context)
    progress["fixture"] = encode(proof)
    require(isinstance(proof, Mapping) and proof.get("passed") is True, "fixture did not pass")
    if isinstance(proof.get("checks"), Mapping):
        require(all(proof["checks"].values()), "fixture has failed checks")
    progress["fixture_completed_before_input_json_decode"] = True

    # Only now is the one pinned receipt decoded.
    progress["input_json_decode_started"] = True
    receipt = json.loads(gzip.decompress(raw[PRIOR_RECEIPT]))
    require(isinstance(receipt, Mapping), "prior receipt is not a mapping")
    require(
        all(key in receipt for key in ("source", "entry", "result", "frozen_code_base64")),
        "prior receipt envelope",
    )
    require(receipt.get("component") == "n1_first_collar_g4_pressure_v1", "wrong prior receipt component")
    prior_result = receipt.get("result")
    require(isinstance(prior_result, Mapping) and prior_result.get("status") == "returned", "prior receipt status")
    prior_result_raw = (
        json.dumps(prior_result, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"
    ).encode("utf-8")
    require(sha(prior_result_raw) == DOCUMENTED_RESULT_SHA256, "prior fixed result hash")
    output = prior_result.get("output")
    require(isinstance(output, Mapping), "prior receipt output missing")
    require("result" not in output, "nested output.result is not accepted")
    source_family = output.get("source_family")
    source_logs = output.get("source_logs")
    sc_raw = output.get("sc")
    domain = output.get("domain")
    finite = output.get("finite_prefix_channels")
    prior_input_index = output.get("input_index")
    for key in ("velocity_prefix_is_W0_plus_W1_once", "finite_prefix_is_unmerged_sector_lists", "pressure_uncut_by_eta", "pressure_direct_density_no_endpoint_added", "source_intervals_not_point_values"):
        require(output.get(key) is True, "prior source scope: " + key)
    for key in ("core_m15_inlet_added", "allword_tail_added", "appended_to_cumulative_m15", "complete_m15_computed", "terminal_closure_proved"):
        require(output.get(key) is False, "prior source scope: " + key)
    require(source_family is not None and source_family != "", "source family missing")
    reject_float_tree(source_family, "source_family")
    require(isinstance(source_logs, Mapping), "source logs must be mapping")
    require("domain" in output, "domain missing")
    require(isinstance(prior_input_index, Mapping), "prior input_index missing")
    require(len(prior_input_index) == 1975, "prior input_index count")
    reject_float_tree(source_logs, "source_logs")
    reject_float_tree(domain, "domain")
    reject_float_tree(prior_input_index, "prior input_index")
    sc = interval(context, sc_raw, "receipt sc")
    channels, input_index = _hydrate_channels(context, finite)

    # Passing the saved source-log mapping lets the leaf validate the actual
    # atom-label shape.  It does not establish owner identity for those atoms.
    source_atoms = source_logs
    progress["actual_leaf_combine_entered"] = True
    grouped = leaf.combine(
        channels["axial"],
        channels["swirl"],
        channels["pressure"],
        source_atoms=source_atoms,
    )
    require(isinstance(grouped, Mapping), "leaf combine did not return a mapping")
    progress["actual_leaf_combine_calls"] = 1
    group_counts = _validate_grouped(grouped, channels)
    input_counts = {
        str(order): {channel: len(channels[channel][order]) for channel in CHANNELS}
        for order in (0, 1)
    }
    member_counts = {
        str(order): sum(
            int(group["origin"]["member_count"])
            for group in grouped["groups"][order]
        )
        for order in (0, 1)
    }
    total_inputs = sum(sum(row.values()) for row in input_counts.values())
    total_members = sum(member_counts.values())
    require(total_inputs == total_members == 1975, "input/member total")
    require(grouped.get("physical_factors_applied") is False, "physical factor applied")
    require(grouped.get("cancellation_performed") is False, "cancellation performed")

    return {
        "source_family": encode(source_family),
        "source_logs": encode(source_logs),
        "source_power_basis": encode(grouped.get("source_power_basis")),
        "grouping_key": "(powers, uncertainty, ordered error_classes tuple)",
        "sc": encode(sc),
        "domain": encode(domain),
        "fixture": encode(proof),
        "input_counts": input_counts,
        "group_counts": group_counts,
        "member_counts": member_counts,
        "total_input_records": total_inputs,
        "total_group_members": total_members,
        "input_index": encode(prior_input_index),
        "prior_input_index": encode(prior_input_index),
        "group_input_index": input_index,
        "input_records": encode(channels),
        "groups": encode(grouped["groups"]),
        "input_receipt_sha256": {PRIOR_RECEIPT: hashes[PRIOR_RECEIPT]},
        "published_receipt": {
            "path": PRIOR_RECEIPT,
            "raw_sha256": RECEIPT_SHA256,
            "result_sha256": DOCUMENTED_RESULT_SHA256,
        },
        "source_atom_owner_identity_proved": False,
        "source_equality_proved": False,
        "source_equality_from_overlap": False,
        "source_intervals_are_not_functions": True,
        "source_intervals_not_point_values": True,
        "complete_kernel_moment_computed": False,
        "already_physically_scaled_channels": True,
        "physical_factors_reapplied": False,
        "inherited_error_conditions": True,
        "errors_independently_certified_here": False,
        "all_input_errors_retained": True,
        "old_errors_removed": False,
        "same_packet_grouping_only": True,
        "actual_leaf_combine_calls": 1,
        "fixture_leaf_combine_calls": 1,
        "ancestor_source_reexecuted": False,
        "appended_to_cumulative_m15": False,
        "complete_m15_computed": False,
        "terminal_closure_proved": False,
    }


def main() -> int:
    hashes, source_sha256, brief_sha256 = preflight()
    entry_sha256 = save(
        ENTRY,
        {
            "component": COMPONENT,
            "source_sha256": source_sha256,
            "brief_sha256": brief_sha256,
            "input_hashes": hashes,
            "entered": True,
            "raw_gzip_authenticated_before_decode": True,
            "prior_receipt_reentered": False,
            "ancestor_source_reexecuted": False,
        },
    )
    progress: dict[str, Any] = {
        "fixture_completed_before_input_json_decode": False,
        "input_json_decode_started": False,
        "actual_leaf_combine_entered": False,
        "actual_leaf_combine_calls": 0,
    }
    outcome: dict[str, Any] = {
        "component": COMPONENT,
        "source_sha256": source_sha256,
        "entry_sha256": entry_sha256,
        "status": "failed",
        "raw_gzip_authenticated_before_decode": True,
        "prior_receipt_reentered": False,
        "progress": progress,
        "output": None,
    }
    try:
        outcome["output"] = consume(hashes, progress)
        outcome["status"] = "returned"
        outcome["error"] = None
    except BaseException as exc:
        outcome["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
    result_sha256 = save(RESULT, outcome)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "status": outcome["status"],
                "entry_sha256": entry_sha256,
                "result_sha256": result_sha256,
                "error": (outcome.get("error") or {}).get("message"),
            }
        )
    )
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())

