"""Produce/inspect new source proofs for the two long/reference seams.

Accepted reports and source archives are read as data.  No accepted owner,
builder, checker, norm, recurrence or quadrature is invoked.  The stage keeps
the complete uncontracted pressure source until its remaining obligations
have an actual source-bound proof.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib
import json
from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "experiments/root_st073"
PARENT_NAME = "lei_ren_part1_paper_compliant_current_original_long_reference_finite_N_I1.json"
PARENT_SHA256 = "75918293b8c7432c637eb551a29cb9d248aa7746dcb1cd4f8b3615f72514c14d"
NAME = "lei_ren_part1_paper_compliant_current_original_long_reference_seam_prerequisites.json"
PACKETS = (
    ("actual_long_reference_native_seams", "original_long_reference_native_seams"),
    ("actual_long_reference_history_seams", "original_long_reference_history_seams"),
    ("actual_long_reference_loop_seams", "original_long_reference_loop_seams"),
)
CODE = tuple("src/openai_ns_reconstruction/" + module + ".py" for _, module in PACKETS)
RUNNER = "scripts/run_original_long_reference_seam_prerequisites.py"
REVIEWED = {
    CODE[0]: "15efca2a18d89f096e0ff945280276933a150c10f46e2171111b7e5efae130b2",
    CODE[1]: "3757d3f58ce8674ac515cdb62ed1dc393ff74b764253f2f6587de2ace5700a95",
    CODE[2]: "8fada76ded014715dd9650bf6688a7acbd9b37bd7af97be3cdb651d16a7efc79",
}


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path):
    return json.loads(Path(path).read_bytes())


def parent():
    require(sha(BASE / PARENT_NAME) == PARENT_SHA256, "Accepted parent bytes changed")
    accepted = load(BASE / PARENT_NAME)
    require(accepted["checks"]["passed"] is True, "Accepted parent is unchecked")
    return accepted


def check_pins(report):
    count = 0
    for folder, group in ((BASE, report["input_hashes"]),
                          (ROOT, report["implementation_hashes"])):
        for name, expected in group.items():
            require(isinstance(expected, str) and sha(folder / name) == expected,
                    "Pinned source bytes changed: " + name)
            count += 1
    return count


def merge_pins(accepted, packets):
    inputs = dict(accepted["input_hashes"])
    implementations = dict(accepted["implementation_hashes"])
    inputs[PARENT_NAME] = PARENT_SHA256
    for packet in packets:
        require(packet["passed"] is True
                and packet["parent_report_sha256"] == PARENT_SHA256,
                "New proof belongs to another predecessor")
        for key, item in packet["source_pins"].items():
            path = Path(item.get("path", key) if isinstance(item, dict) else key)
            expected = item.get("expected_sha256") if isinstance(item, dict) else item
            require(path.as_posix() == key and not path.is_absolute()
                    and ".." not in path.parts,
                    "Source pin key/path disagree: " + key)
            if isinstance(item, dict):
                require(item.get("passed") is True and item.get("source_executed") is not True,
                        "Source pin is unchecked or was executed: " + key)
            if len(path.parts) == 1:
                actual, name, target = BASE / path, path.name, inputs
            elif path.as_posix().startswith("experiments/root_st073/"):
                require(path.parent.as_posix() == "experiments/root_st073",
                        "Unexpected nested input source: " + key)
                actual, name, target = ROOT / path, path.name, inputs
            else:
                actual, name, target = ROOT / path, path.as_posix(), implementations
            require(sha(actual) == expected and target.get(name, expected) == expected,
                    "New proof source pin mismatch: " + path.as_posix())
            target[name] = expected
    for path in (*CODE, RUNNER):
        implementations[path] = sha(ROOT / path)
    return inputs, implementations


def validate_packet(name, packet):
    require(packet["passed"] is True
            and packet["parent_report_sha256"] == PARENT_SHA256,
            "Wrong source-proof predecessor: " + name)
    require(packet["source_pins"], "Source-proof pins missing: " + name)
    require(packet["source_executed"] is False,
            "Archived source execution is outside this stage: " + name)
    require(packet.get("limits"), "Proof scope must remain explicit: " + name)
    if name == "actual_long_reference_native_seams":
        require(packet["ordinary_Z_orders"] == [0, 1, 2]
                and packet["physical_native_C2"]["passed"] is True
                and packet["physical_native_C2"]["physical_y_orders"] == [0, 1, 2]
                and packet["physical_native_C2"]["ordinary_Z_parameter_function_identity"] is True
                and packet["native_U_substitution_chain"]["passed"] is True,
                "Actual native physical-C2 source proof is incomplete")
        for seam in ("post_to_long_q0", "long_q1_to_reference_q0"):
            row = packet["physical_native_C2"][seam]
            require(row["passed"] is True and row["log_value_residual"] == "0",
                    "Native endpoint formula did not agree")
            for key in ("phi_C2_residuals", "U_C2_residuals"):
                required = ({"phi", "phi_y", "phi_yy"} if key == "phi_C2_residuals"
                            else {"U", "U_y", "U_yy"})
                require(set(row[key]) == required
                        and all(isinstance(orders, dict) and set(orders) == {"0", "1", "2"}
                            and all(value == "0" for value in orders.values())
                            for orders in row[key].values()),
                        "Native C2 ordinary-Z residual is not zero")
        require(packet["modified_Qy_claimed"] is False
                and packet["full_modified_seam_claimed"] is False,
                "Native proof was promoted to a modified seam")
    elif name == "actual_long_reference_history_seams":
        require(packet["native_history_source"]["passed"] is True
                and packet["c5_history_source"]["passed"] is True
                and packet["q4_background_seed_source"]["independent_seed_algebra"]["passed"] is True
                and set(packet["q4_background_seed_source"]["independent_seed_algebra"][
                    "residuals"].values()) == {"0"}
                and packet["newproof_checks"]["C5_rate_zero_p_memory_retained"] is True,
                "Actual Duhamel/seed source theorem is incomplete")
        require(packet["limits"]["q4_modified_Qy_seam_proved"] is False
                and packet["limits"]["selected_values_used"] is False
                and packet["limits"]["endpoint_value_identity_proved"] is False,
                "Conditional history theorem was promoted to a completed seam")
    elif name == "actual_long_reference_loop_seams":
        require(packet["periodic_formula"]["sympy_exact_check"]["passed"] is True
                and packet["history_binding"]["passed"] is True
                and packet["endpoint_identity_receipts"]["passed"] is True
                and packet["periodic_formula"]["positive_phase_controls_proved_here"] is False,
                "Actual loop formula or its explicit hypotheses were lost")
        require(packet["limits"]["modified_Uy_Qy_closed"] is False
                and packet["limits"]["interval_or_cap_values_used"] is False,
                "Loop periodicity was promoted to modified Qy closure")


def validate_report(report, accepted):
    mutable = {"stage", "parent_report", "parent_report_sha256", "input_hashes",
               "implementation_hashes", "contract", "checks", "execution_seconds",
               "execution_audit"}
    for key, value in accepted.items():
        if key not in mutable:
            require(report.get(key) == value, "Accepted mathematical payload changed: " + key)
    require(report["execution_audit"]["inherited_parent_execution_audit"]
            == accepted["execution_audit"]
            and report["execution_audit"]["current_turn_accepted_ancestors_replayed"] is False,
            "Accepted execution audit was reset")
    for key, _ in PACKETS:
        validate_packet(key, report[key])
    require(report["actual_long_reference_pressure_source"]["joined_pressure_source"][
                "joined_contracted"] is None,
            "An unproved source endpoint was removed")
    require(report["contract"]["modified_RQy_seams_proved"] is False
            and report["contract"]["global_five_moment_closure"] is False
            and report["contract"]["full_goal_complete"] is False,
            "Prerequisite proofs were promoted to full construction")
    return dict(passed=True, accepted_math_payload_preserved=True,
                full_uncontracted_graph_preserved=True, source_proof_packets=len(PACKETS),
                modified_RQy_seams_proved=False, full_goal_complete=False)


def run():
    started = time.perf_counter()
    for path, expected in REVIEWED.items():
        require(expected != "REVIEW_PENDING" and sha(ROOT / path) == expected,
                "Source review pending or bytes changed: " + path)
    accepted = parent()
    check_pins(accepted)
    sys.path.insert(0, str(ROOT / "src"))
    packets = []
    for _, module in PACKETS:
        packets.append(importlib.import_module("openai_ns_reconstruction." + module).build(ROOT, accepted))
    require(accepted == parent(), "New proof changed its accepted input")
    inputs, implementations = merge_pins(accepted, packets)
    report = copy.deepcopy(accepted)
    report.update(stage="actual_long_reference_native_history_loop_seam_prerequisites",
                  parent_report=PARENT_NAME, parent_report_sha256=PARENT_SHA256,
                  input_hashes=inputs, implementation_hashes=implementations)
    for (key, _), packet in zip(PACKETS, packets):
        report[key] = packet
    report["execution_audit"] = dict(current_turn_accepted_ancestors_replayed=False,
        new_exact_source_proofs_only=True,
        inherited_parent_execution_audit=copy.deepcopy(accepted["execution_audit"]))
    report["contract"] = dict(source_bound_seam_prerequisites=True,
        modified_RQy_seams_proved=False, internal_endpoints_contracted=False,
        full_pressure_numeric=False, global_five_moment_closure=False,
        full_goal_complete=False)
    report["checks"] = dict(passed=True,
        packet_invariants=validate_report(report, accepted),
        independent_new_proof_review_completed=True,
        independent_runtime_worker_model_metadata_available=False,
        reviewed_implementation_sha256=REVIEWED)
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    json.dumps(report, allow_nan=False)
    return report


def inspect(path):
    report, accepted = load(path), parent()
    require(report["checks"]["passed"] is True
            and report["checks"]["reviewed_implementation_sha256"] == REVIEWED,
            "Saved source reviews changed")
    for path in (*CODE, RUNNER):
        require(report["implementation_hashes"][path] == sha(ROOT / path),
                "Saved implementation changed: " + path)
    return dict(passed=True, pins=check_pins(report),
        packet_invariants=validate_report(report, accepted), serialized_only=True,
        accepted_builders_or_checkers_replayed=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report", action="store_true")
    group.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    path = BASE / NAME
    if args.write_report:
        require(not path.exists(), "Saved proof report already exists; use --inspect-report")
        report = run()
        tmp = path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                       encoding="utf-8", newline="\n")
        tmp.replace(path)
        print(json.dumps(dict(passed=True, execution_seconds=report["execution_seconds"],
                              bytes=path.stat().st_size, sha256=sha(path))))
    else:
        print(json.dumps(inspect(path), indent=2))
