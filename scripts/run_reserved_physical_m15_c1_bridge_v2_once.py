"""One-entry reader for the accepted reserved physical ``m15`` C1 bridge.

The root agent freezes ``tmp/reserved_physical_m15_c1_bridge_v2/source.json``
before entry.  This runner first authenticates that manifest, every listed
file, and the four fixed archive/check SHA-256 values.  It then writes an
exclusive entry, imports the new consumer, parses the terminal/check payloads,
and calls the accepted bridge once.  No old producer, fixture, checker, or
integral helper is imported or executed.

The three later calls are independent expected-failure checks on in-memory
copies only: one changes the exact physical sign, one changes the value body,
and one changes the ordinary-Z1 body.  The mutated payloads are never written
and are never returned as source data.
"""

from __future__ import annotations

import copy
import gzip
import hashlib
import importlib.util
import json
from dataclasses import asdict
from pathlib import Path
import sys
import traceback
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "reserved_physical_m15_c1_bridge_v2"
DIR = ROOT / "tmp" / COMPONENT
SOURCE = DIR / "source.json"
ENTRY = DIR / "entered.json"
RESULT = DIR / "result.json"

CONSUMER = "src/openai_ns_reconstruction/reserved_physical_m15_c1_bridge_v2.py"
RUNNER = "scripts/run_reserved_physical_m15_c1_bridge_v2_once.py"
TERMINAL = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_terminal_integral_identities.json.gz"
TERMINAL_CHECK = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_terminal_integral_identities_check.json"
OWN = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_reserved_own_fields.json.gz"
FULL = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_full_defect_linear_functions.json.gz"

REQUIRED = frozenset((CONSUMER, RUNNER, TERMINAL, TERMINAL_CHECK, OWN, FULL))
RAW_ARCHIVE_SHA256 = {
    TERMINAL: "58d90a1d34bbbece45798b3fbd7db75ec751a321eb1740b26f6ac434b76c3b81",
    TERMINAL_CHECK: "a5c0ca5789b00fb2a3cdeffd07d3f4dff77bef2d169f5b86a8048151afe5e5fe",
    OWN: "5d43944d68e9e1c85b52493d1e1a4b86a4928209dd6fb97cba9e1d4c63bdad32",
    FULL: "d6b193bf527d4e8cd4afe415b1e941c928f32288b63362783074381ef7480fc1",
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def exclusive(path: Path, value: Any) -> str:
    raw = (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")
    with path.open("xb") as stream:
        stream.write(raw)
    return sha(raw)


def parse_object(path: Path, raw: bytes) -> dict[str, Any]:
    value = json.loads(raw.decode("utf-8"))
    require(isinstance(value, dict), f"{path} is not a JSON object")
    return value


def preflight() -> tuple[dict[str, str], str, dict[str, Any]]:
    require(not ENTRY.exists() and not RESULT.exists(), "one-entry guard: already entered")
    source_raw = SOURCE.read_bytes()
    source = parse_object(SOURCE, source_raw)
    require(source.get("component") == COMPONENT, "source component mismatch")
    require(source.get("one_entry_only") is True, "source one-entry guard missing")
    files = source.get("files")
    require(isinstance(files, list), "source files manifest missing")
    hashes: dict[str, str] = {}
    for item in files:
        require(isinstance(item, dict), "source file entry is not an object")
        relative = item.get("path")
        expected = item.get("sha256")
        require(isinstance(relative, str) and isinstance(expected, str), "invalid source file pin")
        require(relative in REQUIRED and relative not in hashes, f"unexpected or duplicate source file: {relative}")
        actual = sha((ROOT / relative).read_bytes())
        require(actual == expected, f"changed source file: {relative}")
        hashes[relative] = actual
    require(set(hashes) == REQUIRED, "source manifest does not cover exactly six files")
    for relative, expected in RAW_ARCHIVE_SHA256.items():
        require(hashes[relative] == expected, f"fixed raw archive/check SHA mismatch: {relative}")
    if source.get("raw_archive_sha256") is not None:
        require(source["raw_archive_sha256"] == RAW_ARCHIVE_SHA256, "source raw archive map mismatch")
    return hashes, sha(source_raw), source


def import_consumer(hashes: dict[str, str]):
    path = ROOT / CONSUMER
    require(sha(path.read_bytes()) == hashes[CONSUMER], "consumer changed after entry")
    spec = importlib.util.spec_from_file_location(COMPONENT + "_consumer", path)
    require(spec is not None and spec.loader is not None, "consumer import spec unavailable")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def authenticated_archives(module: Any, hashes: dict[str, str]):
    parsed: dict[str, dict[str, Any]] = {}
    for relative in (TERMINAL, TERMINAL_CHECK, OWN, FULL):
        path = ROOT / relative
        raw = path.read_bytes()
        require(sha(raw) == hashes[relative] == RAW_ARCHIVE_SHA256[relative], f"archive changed after entry: {relative}")
        if relative.endswith(".gz"):
            raw = gzip.decompress(raw)
        parsed[relative] = parse_object(path, raw)
    Archive = module.AuthenticatedArchive
    def make(relative: str):
        return Archive(name=Path(relative).name, sha256=hashes[relative], payload=parsed[relative], hash_authenticated=True)
    return make(TERMINAL), make(TERMINAL_CHECK), make(OWN), make(FULL)


def mutated_archive(archive: Any, payload: dict[str, Any]) -> Any:
    return type(archive)(
        name=archive.name,
        sha256=archive.sha256,
        payload=payload,
        hash_authenticated=True,
    )


def expected_negative_cases(module: Any, terminal: Any, check: Any, own: Any, full: Any, bridge: Any) -> list[dict[str, Any]]:
    sign = bridge.physical.sign_handle
    body = bridge.signed_physical_m15_patch.integrand
    per_ds = bridge.signed_physical_m15_patch.per_ds

    def run_case(label: str, mutate) -> dict[str, Any]:
        # Mutation is deliberately outside the rejection handler.  A malformed
        # negative fixture (for example, a missing graph row) is an actual
        # runner failure, not evidence that the bridge rejected the mutation.
        t_payload = copy.deepcopy(terminal.payload)
        o_payload = copy.deepcopy(own.payload)
        f_payload = copy.deepcopy(full.payload)
        mutate(t_payload, o_payload, f_payload)
        try:
            module.bind_reserved_physical_m15_c1(
                terminal_archive=mutated_archive(terminal, t_payload),
                terminal_check=check,
                reserved_own_archive=mutated_archive(own, o_payload),
                full_defect_archive=mutated_archive(full, f_payload),
            )
        except module.M15BridgeError as exc:
            return {
                "case": label,
                "rejected": True,
                "error": {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()},
            }
        return {"case": label, "rejected": False, "error": None}

    def mutate_sign(terminal_payload, own_payload, full_payload):
        def set_sign(payload, label):
            nodes = payload.get("exact_graph_nodes")
            require(isinstance(nodes, list), f"{label} exact_graph_nodes missing in sign negative")
            require(type(sign) is int and 0 <= sign < len(nodes), f"{label} sign handle path is absent")
            row = nodes[sign]
            require(isinstance(row, dict), f"{label} sign handle row is not an object")
            require(row.get("operation") == "exact_rational", f"{label} sign handle is not exact rational")
            require("numerator" in row, f"{label} sign handle has no numerator")
            row["numerator"] = 1

        # own is the source owner; terminal must contain it through the
        # verified prefix.  A shorter full graph is not indexed by this
        # own-only handle, while a containing full prefix must be changed too
        # or the mutation would be rejected merely by prefix mismatch.
        set_sign(terminal_payload, "terminal")
        set_sign(own_payload, "own")
        full_nodes = full_payload.get("exact_graph_nodes")
        require(isinstance(full_nodes, list), "full exact_graph_nodes missing in sign negative")
        if type(sign) is int and sign < len(full_nodes):
            set_sign(full_payload, "full")

    def mutate_value_body(terminal_payload, own_payload, full_payload):
        terminal_payload["exact_graph_nodes"][body.value]["expression"] = per_ds.Z1
        own_payload["exact_graph_nodes"][body.value]["expression"] = per_ds.Z1

    def mutate_Z1_body(terminal_payload, own_payload, full_payload):
        terminal_payload["exact_graph_nodes"][body.Z1]["expression"] = per_ds.value
        own_payload["exact_graph_nodes"][body.Z1]["expression"] = per_ds.value

    return [
        run_case("negative_physical_sign", mutate_sign),
        run_case("negative_value_body", mutate_value_body),
        run_case("negative_Z1_body", mutate_Z1_body),
    ]


def consume(hashes: dict[str, str], partial: dict[str, Any]) -> dict[str, Any]:
    partial["stage"] = "import_consumer"
    module = import_consumer(hashes)
    partial["stage"] = "decompress_and_parse_authenticated_archives"
    terminal, check, own, full = authenticated_archives(module, hashes)
    partial["stage"] = "accepted_bridge_call"
    bridge = module.bind_reserved_physical_m15_c1(
        terminal_archive=terminal,
        terminal_check=check,
        reserved_own_archive=own,
        full_defect_archive=full,
    )
    accepted_bridge = asdict(bridge)
    partial["accepted_bridge"] = accepted_bridge
    partial["stage"] = "independent_expected_failure_cases"
    negatives = expected_negative_cases(module, terminal, check, own, full, bridge)
    partial["negative_cases"] = negatives
    require(all(case["rejected"] is True for case in negatives), "a graph mutation was unexpectedly accepted")
    partial["stage"] = "returned"
    return {
        "accepted_bridge": accepted_bridge,
        "negative_cases": negatives,
        "production_consumer_calls": 1,
        "negative_consumer_calls": len(negatives),
        "old_producer_or_fixture_called": False,
        "integral_recomputed": False,
    }


def main() -> int:
    out: dict[str, Any] = {"component": COMPONENT, "status": "failed", "partial_output": {}}
    source_sha = None
    entry_sha = None
    try:
        hashes, source_sha, source = preflight()
        out["source_sha256"] = source_sha
        out["partial_output"] = {"stage": "preflight_passed", "source": source, "input_hashes": hashes}
        entry_sha = exclusive(ENTRY, {
            "component": COMPONENT,
            "entered": True,
            "source_sha256": source_sha,
            "input_hashes": hashes,
            "production_consumer_calls": 0,
            "negative_consumer_calls": 0,
        })
        out["entry_sha256"] = entry_sha
        out["partial_output"]["stage"] = "exclusive_entry_written"
        out["proof"] = consume(hashes, out["partial_output"])
        out["status"] = "returned"
        out["error"] = None
    except BaseException as exc:
        out["error"] = {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()}
    result_sha = None
    if not RESULT.exists():
        result_sha = exclusive(RESULT, out)
    print(json.dumps({
        "component": COMPONENT,
        "status": out["status"],
        "source_sha256": source_sha,
        "entry_sha256": entry_sha,
        "result_sha256": result_sha,
    }, sort_keys=True))
    return 0 if out["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
