"""Capture the unfinished physical-only consumer's first failed access.

The accepted radial/slow packets are immutable data. This explicitly repeats
only the unaccepted physical attempt once to preserve its missing traceback.
"""
from __future__ import annotations

import argparse
import ast
from collections.abc import Mapping
import copy
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "experiments/root_st073"
PREFIX = "lei_ren_part1_paper_compliant_current_original_"
PARENT = PREFIX + "Rm_modified_physical_values.json"
PARENT_HASH = "5c26b7bb3b70e874a5790a5f8424233678b5a964200b8ffd88dce102b22110f1"
FINITE = PREFIX + "true_Rm_patch_finite_N_C5_functions.json.gz"
FINITE_HASH = "8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad"
MODIFIED = PREFIX + "modified_Q4_pressure_C3_functions.json.gz"
MODIFIED_HASH = "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a"
PORTS = "experiments/root_st073/original_Rm_loop_ports_admission.json"
PORTS_HASH = "0331c5fc37bd3ec6c0d386c38b07d9863ede622e14f1c6dc65bf8cdf675a16b9"
GRAPH_HASH = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
INPUT_HASH = "081b95e2df85e4f9573d31277db95e009e181a08ab6d84f0a02703b407ea1195"
KERNEL = "src/openai_ns_reconstruction/original_Rm_physical_qy0_value_seam.py"
KERNEL_HASH = "3a6cb37c906d660ddf96758ed0612933fdba3d35772f58a5b14e4749ff66095f"
RUNNER = "scripts/diagnose_original_Rm_physical_qy0_access.py"
CHECKPOINTS = ROOT / "tmp/original_Rm_physical_qy0_access_diagnostic"
SOURCE = CHECKPOINTS / "new_source.json"
ENTRY = CHECKPOINTS / "new_operator_entered.json"
RESULT = BASE / "original_Rm_physical_Qy_Z0_traceback_attempt.json"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()).hexdigest()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def path_for(name):
    path = (ROOT / name if name.startswith(("src/", "scripts/", "experiments/", "docs/")) else BASE / name).resolve()
    require(path.is_relative_to(ROOT.resolve()), "input path escapes repository")
    return path


def read_pinned(path, expected):
    raw = path.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == expected, "source pin changed: " + str(path))
    return json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw)


def write_new(path, value):
    require(not path.exists(), "artifact already exists: " + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True)
        stream.write("\n")


def load_inputs():
    for name in (RUNNER, KERNEL):
        raw = (ROOT / name).read_bytes()
        require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "LF/EOF changed")
        ast.parse(raw, filename=name)
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    for name, expected in {**parent["input_hashes"], **parent["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "inherited pin changed: " + name)
    require(sha(ROOT / KERNEL) == KERNEL_HASH, "frozen failed kernel changed")
    require(parent["actual_Rm_modified_radial_value_seam"]["passed"] is True, "accepted radial packet missing")
    require(parent["actual_Rm_physical_Qy_Z0_value_seam"]["passed"] is False, "physical attempt already passed")
    finite = read_pinned(BASE / FINITE, FINITE_HASH)
    modified = read_pinned(BASE / MODIFIED, MODIFIED_HASH)
    require(finite["source_family"] == modified["source_family"] == parent["source_family"], "source family changed")
    nodes = finite["exact_graph_nodes"]
    require(len(nodes) == 205147 and nodes[:183700] == modified["exact_graph_nodes"], "saved prefixes changed")
    for packet, base_key, full_key in (
        (parent["actual_Rm_native_radial_endpoint_source"], "base_digest", "full_nodes_digest"),
        (parent["actual_Rm_second_radial_endpoint_functions"], "base_digest", "full_nodes_digest"),
        (parent["actual_Rm_history_native_source"]["graph_delta"], "base_nodes_sha256", "full_canonical_nodes_sha256"),
        (parent["actual_Rm_history_seam_identity"]["graph_delta"], "base_nodes_sha256", "full_nodes_sha256"),
        (parent["actual_Rm_modified_moment_germ"]["graph_delta"], "base_nodes_sha256", "full_nodes_sha256"),
    ):
        require(digest(nodes) == packet[base_key]
            and len(nodes) == packet.get("base_node_count", packet.get("base_graph_node_count")), "saved delta base changed")
        nodes = nodes + copy.deepcopy(packet["appended_nodes"])
        require(len(nodes) == packet["full_node_count"] and digest(nodes) == packet[full_key], "saved delta changed")
    require(len(nodes) == 205300 and digest(nodes) == GRAPH_HASH, "saved graph changed")
    inputs = dict(bounded_loop_ports_admission=read_pinned(ROOT / PORTS, PORTS_HASH),
        modified_moment_germ=parent["actual_Rm_modified_moment_germ"],
        modified_radial_value_seam=parent["actual_Rm_modified_radial_value_seam"],
        actual_source_bindings=dict(postrestore=modified["modified_Q4_pressure_C3_bindings"]["charts"]["postrestore"],
            patch=finite["Rm_patch_finite_N_C5_bindings"]), expected_nodes_digest=GRAPH_HASH)
    require(digest(inputs) == INPUT_HASH, "diagnostic input differs from preserved failed attempt")
    packet = dict(schema="Rm_unaccepted_physical_Qy_Z0_traceback_inputs_v1", parent_report=PARENT,
        parent_sha256=PARENT_HASH, graph_digest=GRAPH_HASH, graph_node_count=len(nodes), input_digest=INPUT_HASH,
        implementation_hashes={RUNNER:sha(ROOT / RUNNER), KERNEL:KERNEL_HASH},
        accepted_radial_and_slow_consumed_as_data=True, accepted_operators_replayed=False,
        prior_unaccepted_physical_consumer_reentry_for_traceback=True, full_goal_complete=False)
    return nodes, inputs, packet


def prepare():
    *_, packet = load_inputs()
    require(not any(p.exists() for p in (SOURCE, ENTRY, RESULT)), "diagnostic already prepared")
    write_new(SOURCE, packet)
    return dict(passed=True, source_sha256=sha(SOURCE), input_digest=INPUT_HASH)


def capture_frames(error):
    frames = []
    tb = error.__traceback__
    while tb is not None:
        frame = tb.tb_frame
        mappings = {}
        for key, value in frame.f_locals.items():
            if isinstance(value, Mapping):
                mappings[key] = dict(size=len(value), keys=[repr(k) for k in list(value)[:64]])
        frames.append(dict(file=frame.f_code.co_filename, function=frame.f_code.co_name,
            line=tb.tb_lineno, side=frame.f_locals.get("side") if isinstance(frame.f_locals.get("side"), str) else None,
            mapping_key_summaries=mappings))
        tb = tb.tb_next
    return frames


def run():
    nodes, inputs, packet = load_inputs()
    require(json.loads(SOURCE.read_bytes()) == packet, "diagnostic source changed")
    require(not ENTRY.exists() and not RESULT.exists(), "diagnostic already entered")
    spec = importlib.util.spec_from_file_location("traceback_only_unaccepted_Rm_physical_Qy", ROOT / KERNEL)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    write_new(ENTRY, dict(source_sha256=sha(SOURCE), input_digest=INPUT_HASH,
        kernel_sha256=KERNEL_HASH, entered_at_utc=datetime.now(timezone.utc).isoformat()))
    started = time.perf_counter()
    try:
        returned = module.derive(nodes, **inputs)
        result = dict(exception_captured=False, returned_packet=returned)
    except Exception as error:
        result = dict(exception_captured=True, error_type=type(error).__name__, error=str(error),
            traceback="".join(traceback.format_exception(error)), frames=capture_frames(error))
    result.update(schema="Rm_unaccepted_physical_Qy_Z0_traceback_attempt_v1",
        execution_seconds=time.perf_counter()-started, source_checkpoint=packet,
        source_checkpoint_sha256=sha(SOURCE), entry_checkpoint=json.loads(ENTRY.read_bytes()),
        entry_checkpoint_sha256=sha(ENTRY), input_unchanged=digest(inputs) == INPUT_HASH,
        graph_unchanged=digest(nodes) == GRAPH_HASH, accepted_operators_replayed=False,
        prior_unaccepted_physical_consumer_reentered_once_for_traceback=True, full_goal_complete=False)
    write_new(RESULT, result)
    return dict(saved=True, sha256=sha(RESULT), exception_captured=result["exception_captured"],
        error_type=result.get("error_type"), error=result.get("error"), seconds=result["execution_seconds"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run"))
    print(json.dumps(globals()[parser.parse_args().action]()))
