"""Produce/inspect the actual modified long/reference seam proof stage.

Accepted arithmetic is never replayed.  The four new builders consume saved
source archives and the exact prerequisite report.  This producer requires
the complete source-owned seam proof before enabling a pressure contraction.
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
PARENT_NAME = "lei_ren_part1_paper_compliant_current_original_long_reference_seam_prerequisites.json"
PARENT_SHA256 = "0a91144aca1caaf0dd2b11170b1aee4339f676eb5f213eea236ae8c814cea0c4"
NAME = "lei_ren_part1_paper_compliant_current_original_long_reference_modified_seams.json"
PACKETS = (
    ("actual_long_reference_control_seams", "original_long_reference_control_seams"),
    ("actual_physical_y_dag_derivative", "original_physical_y_dag_derivative"),
    ("actual_long_reference_modified_seams", "original_long_reference_modified_seams"),
    ("actual_long_reference_correlated_pressure", "original_long_reference_correlated_pressure"),
)
CODE = tuple("src/openai_ns_reconstruction/" + module + ".py" for _, module in PACKETS)
RUNNER = "scripts/run_original_long_reference_modified_seams.py"
REVIEWED = {
    CODE[0]: "b192aa1805e357f05102c7ce904388ce4ed50536872c033ad85e475db962f5f0",
    CODE[1]: "1fd508af8b3261006649de58a433bb05b20ae0d077852e847d8026f9f153b49a",
    CODE[2]: "f4e17ce2d1386d8f257909dacbad5bfebe2920dd61c21f538e0396396a445f09",
    CODE[3]: "1086a13ac8f468a631d9bd5f943242da5aecdb8c178daf4eea38dced871e1ec6",
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
    require(accepted["checks"]["passed"] is True, "Accepted prerequisite report is unchecked")
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
        require(packet["passed"] is True and packet["parent_report_sha256"] == PARENT_SHA256,
                "Wrong source-proof predecessor")
        for key, item in packet["source_pins"].items():
            path = Path(item["path"])
            expected = item["expected_sha256"]
            require(path.as_posix() == key and not path.is_absolute()
                    and ".." not in path.parts and item["passed"] is True
                    and item.get("source_executed") is not True,
                    "Malformed or unchecked source pin: " + key)
            if path.parent.as_posix() == "experiments/root_st073":
                name, target = path.name, inputs
            else:
                require(not path.as_posix().startswith("experiments/root_st073/"),
                        "Nested input alias is unsupported: " + key)
                name, target = path.as_posix(), implementations
            require(sha(ROOT / path) == expected and target.get(name, expected) == expected,
                    "Source pin conflict: " + key)
            target[name] = expected
    for path in (*CODE, RUNNER):
        implementations[path] = sha(ROOT / path)
    return inputs, implementations


def validate_report(report, accepted):
    require(report["parent_report"] == PARENT_NAME
            and report["parent_report_sha256"] == PARENT_SHA256,
            "Saved top-level predecessor provenance changed")
    require(all(REVIEWED[path] != "REVIEW_PENDING"
                and report["implementation_hashes"][path] == REVIEWED[path]
                for path in CODE), "Saved implementation hashes were not independently reviewed")
    mutable = {"stage", "parent_report", "parent_report_sha256", "input_hashes",
               "implementation_hashes", "contract", "checks", "execution_seconds", "execution_audit"}
    for key, value in accepted.items():
        if key not in mutable:
            require(report.get(key) == value, "Accepted mathematical payload changed: " + key)
    for key, _ in PACKETS:
        packet = report[key]
        require(packet["passed"] is True and packet["parent_report_sha256"] == PARENT_SHA256
                and packet["source_executed"] is False and packet["source_family"] == accepted["source_family"],
                "Source-proof packet provenance failed: " + key)
    dag = report[PACKETS[1][0]]
    require(dag["conditional_passed"] is True
            and dag["counts"]["slow_derivative_roots_total"] == 36
            and dag["counts"]["slow_derivative_roots_passed"] == 36
            and dag["seed_identity_discharged"] is False,
            "Issued derivative roots or their conditional scope changed")
    proof = report[PACKETS[2][0]]
    require(proof["full_modified_seam_proved"] is True
            and not proof["unresolved_gates"]
            and all(value is True for value in proof["gates"].values())
            and proof["axial_scope"]["Z_interval"] == ["37/100", "93/250"]
            and proof["axial_scope"]["ordinary_Z_orders"] == [0, 1]
            and report["contract"]["certified_Z_interval"] == ["37/100", "93/250"],
            "Full modified endpoint congruence is not established")
    phase = proof["actual_source_composition"]["phase_formula_receipt"]
    require(phase["passed"] is True, "Actual mixed phase formula is not bound")
    for chart in ("post_power", "long_reshape", "reference"):
        mixed = phase["charts"][chart]["phase_partial_graph"]
        convolution = mixed["mixed_lambda_source"]
        require(mixed["passed"] is True and mixed["row_count"] == 21
                and len(mixed["rows"]) == 21
                and all(row["passed"] is True for row in mixed["rows"])
                and convolution["passed"] is True
                and len(convolution["raw_lambda_rows"]) == 15
                and all(row["passed"] is True for row in convolution["raw_lambda_rows"])
                and not convolution["failed"],
                "Actual normalized mixed phase recurrence is incomplete: " + chart)
    for name in ("post_power_upper_to_long_lower", "long_upper_to_reference_lower"):
        seam = proof["seams"][name]
        require(seam["finite_N_Qy_identity_proven"] is True
                and seam["physical_Qy_C1_equal"] is True
                and seam["physical_RQy_C1_equal"] is True,
                "Physical ordinary-C1 seam proof is missing: " + name)
    correlation = report[PACKETS[3][0]]
    patch = correlation["graph_patch"]
    require(patch["base_node_count"] == 185690 and patch["base_prefix_unchanged"] is True
            and patch["base_canonical_nodes_sha256"] == "8422242e6c2bce6750323cf2aafd2a5a4bce2bf406f033ddd6a18ad61ae33d28"
            and patch["full_node_count"] == 185690 + len(patch["appended_nodes"])
            and correlation["joined_contracted"] is not None
            and correlation["exact_linearization"]["targets"] == [160287, 160288, 166149, 166150]
            and correlation["outer_reference_RQy_retained"] == [171869, 171870]
            and correlation["exact_linearization"]["passed"] is True
            and correlation["exact_linearization"]["only_sum_negative_expanded"] is True
            and correlation["exact_linearization"]["conservative_integer_dependency_check"] is True
            and correlation["limits"]["standalone_caps_retained"] is True
            and correlation["limits"]["all_old_data_retained"] is True,
            "Pressure correlation/prefix/outer trace changed")
    for row in correlation["exact_linearization"]["rows"]:
        require(all(value == 0 for value in row["target_coefficients"].values())
                and not any(row["surviving_target_dependencies"].values()),
                "A removed internal endpoint survives in the joined expression")
    require(report["execution_audit"]["inherited_parent_execution_audit"] == accepted["execution_audit"]
            and report["execution_audit"]["current_turn_accepted_ancestors_replayed"] is False,
            "Execution audit was reset")
    require(report["contract"]["full_pressure_numeric"] is False
            and report["contract"]["complete_pressure_numeric"] is False
            and report["contract"]["global_closure"] is False
            and report["contract"]["global_five_moment_closure"] is False
            and report["contract"]["full_goal_complete"] is False,
            "An endpoint identity was promoted to global construction")
    return dict(passed=True, source_proof_packets=4, actual_held_phase_roots=36,
                actual_mixed_phase_coefficients=63, actual_lambda_coefficients=45,
                two_modified_RQy_seams_proved=True, accepted_math_payload_preserved=True,
                all_standalone_caps_retained=True, external_reference_RQy_retained=True,
                full_goal_complete=False)


def run():
    started = time.perf_counter()
    for path, expected in REVIEWED.items():
        require(expected != "REVIEW_PENDING" and sha(ROOT / path) == expected,
                "Source review pending or bytes changed: " + path)
    accepted = parent()
    check_pins(accepted)
    sys.path.insert(0, str(ROOT / "src"))
    modules = [importlib.import_module("openai_ns_reconstruction." + module) for _, module in PACKETS]
    control = modules[0].build(ROOT, accepted)
    dag = modules[1].build(ROOT, accepted)
    proof = modules[2].build(ROOT, accepted, control, dag)
    require(proof["full_modified_seam_proved"] is True,
            "Actual seam proof remains open; no contraction/report is permitted")
    correlation = modules[3].build(ROOT, accepted, proof)
    packets = (control, dag, proof, correlation)
    require(accepted == parent(), "New stage changed its accepted input")
    inputs, implementations = merge_pins(accepted, packets)
    report = copy.deepcopy(accepted)
    report.update(stage="actual_long_reference_modified_RQy_seams_and_correlated_pressure",
                  parent_report=PARENT_NAME, parent_report_sha256=PARENT_SHA256,
                  input_hashes=inputs, implementation_hashes=implementations)
    for (key, _), packet in zip(PACKETS, packets):
        report[key] = packet
    report["execution_audit"] = dict(current_turn_accepted_ancestors_replayed=False,
        new_exact_source_proofs_and_graph_algebra_only=True,
        accepted_graph_rows_restored_as_serialized_data=True,
        inherited_parent_execution_audit=copy.deepcopy(accepted["execution_audit"]))
    report["contract"] = dict(actual_modified_RQy_seam_identities=True,
        internal_endpoint_contractions=True, complete_pressure_numeric=False,
        full_pressure_numeric=False, global_closure=False,
        certified_Z_interval=["37/100", "93/250"], ordinary_Z_orders=[0, 1],
        global_five_moment_closure=False, full_goal_complete=False)
    report["checks"] = dict(passed=True, packet_invariants=validate_report(report, accepted),
        independent_new_math_review_completed=True,
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
    return dict(passed=True, pins=check_pins(report), packet_invariants=validate_report(report, accepted),
                serialized_only=True, accepted_builders_or_checkers_replayed=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report", action="store_true")
    group.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    path = BASE / NAME
    if args.write_report:
        require(not path.exists(), "Saved seam report already exists; use --inspect-report")
        report = run()
        tmp = path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                       encoding="utf-8", newline="\n")
        tmp.replace(path)
        print(json.dumps(dict(passed=True, execution_seconds=report["execution_seconds"],
                              bytes=path.stat().st_size, sha256=sha(path))))
    else:
        print(json.dumps(inspect(path), indent=2))
