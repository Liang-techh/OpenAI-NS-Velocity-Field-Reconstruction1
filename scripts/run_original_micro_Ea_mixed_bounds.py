"""Save or inspect actual ordinary-source E/a mixed bounds.

This runner preserves the accepted C0 packet and records the partial mixed
provider.  It deliberately leaves b/t0/q/p2, physical C_y Z3, mean Z3,
QZZ, and the global five-moment closure open.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "experiments/root_st073"
sys.path.insert(0, str(ROOT / "src"))

# Import the current source module as one object so a source-provenance update
# can add fields without making this serialized-packet checker depend on names
# that are not part of the packet contract.
from openai_ns_reconstruction import original_micro_Ea_mixed_bounds as current


build = current.build
PARENT = current.PARENT
PARENT_SHA256 = current.PARENT_SHA256
HELPER = current.HELPER
HELPER_SHA256 = current.HELPER_SHA256
REVIEWED_MODULE_SHA256 = "449e991587529536eddd2b9ec7e90c0fb0e856fadb223f63951d4a667eaa247d"
require = current.require
sha = current.sha

NAME = "lei_ren_part1_paper_compliant_current_original_micro_Ea_mixed_bounds.json"
CODE = (
    "src/openai_ns_reconstruction/original_micro_Ea_mixed_bounds.py",
    "scripts/run_original_micro_Ea_mixed_bounds.py",
    HELPER,
)


def _load(path):
    return json.loads(Path(path).read_bytes())


def parent():
    require(sha(ROOT / PARENT) == PARENT_SHA256, "Accepted micro C0 parent changed")
    report = _load(ROOT / PARENT)
    require(report["checks"]["passed"] is True, "Checked C0 parent required")
    return report


def pins(report):
    count = 0
    for base, key in ((DIRECTORY, "input_hashes"), (ROOT, "implementation_hashes")):
        for name, expected in report[key].items():
            require(sha(base / name) == expected, "Changed source pin: " + name)
            count += 1
    return count


def _passed(value, label):
    require(isinstance(value, dict) and value.get("passed") is True, label)


def _row_orders(rows, orders, label):
    require(isinstance(rows, list), label + " rows missing")
    require([row.get("ordinary_Z_order") for row in rows] == list(orders), label + " orders missing")


def _packet_fields(report, accepted):
    require(
        report["parent_report"] == Path(PARENT).name
        and report["parent_report_sha256"] == PARENT_SHA256,
        "Same signed C0 predecessor required",
    )
    mutable = {
        "stage",
        "parent_report",
        "parent_report_sha256",
        "input_hashes",
        "implementation_hashes",
        "contract",
        "checks",
        "execution_seconds",
    }
    for key, value in accepted.items():
        if key not in mutable:
            require(report[key] == value, "Inherited signed packet or error changed: " + key)

    packet = report["actual_micro_Ea_mixed_bounds"]
    require(packet.get("actual_Ea_root_bound_provider") is True, "Actual E/a root provider required")
    _passed(packet, "Actual ordinary-source E/a packet missing")
    require(packet["source_family"] == accepted["source_family"], "Source family changed")

    provenance = packet.get("source_provenance")
    _passed(provenance, "Current source provenance receipt missing")
    require(provenance.get("source_producers_replayed") is False, "Accepted source producer replayed")

    provider = packet.get("analytic_derivative_provider")
    _passed(provider, "Analytic derivative provider missing")
    fresh = provider.get("new_rows")
    require(isinstance(fresh, list) and len(fresh) == 3, "Exactly three new Xh rows required")
    require(
        {
            (row.get("scaled_radial_order"), row.get("ordinary_Z_order")) for row in fresh
        }
        == {(0, 4), (1, 3), (1, 4)},
        "Unexpected new Xh rows",
    )
    require(
        provider.get("old_rows_consumed_without_recomputation") is True
        and provider.get("ordinary_derivatives_no_factorial_discount") is True,
        "Old ordinary rows must be consumed and factorial weight retained",
    )

    state = packet.get("comparison_source_jet_proof")
    _passed(state, "Comparison source jet proof missing")
    mixed = packet.get("actual_Ea_mixed_bounds")
    _passed(mixed, "E/a mixed provider missing")
    require(mixed.get("ordinary_Z_orders") == [0, 1, 2, 3], "Ordinary Z0..Z3 provider required")
    require(mixed.get("physical_y_orders") == [0, 1], "Physical y0..y1 provider required")
    for key in (
        "E_Z_over_E",
        "E_yZ_over_E",
        "a_Z_over_a",
        "a_Z_absolute_factored_chi_R",
        "a_yZ_log_absolute_upper_first",
        "a_yZ_log_absolute_upper_second",
    ):
        _row_orders(mixed.get(key), range(4), "E/a " + key)
    require(mixed.get("b_t0_q_p2_bounds_complete") is False, "b/t0/q/p2 gate must remain open")
    require(mixed.get("actual_pre_finite_N_background_Ea") is True
            and mixed.get("finite_N_modified_E_N_not_bounded") is True,
            "This provider bounds the actual background roots before finite-N correction")

    # These are intentionally false: this packet is a partial E/a provider.
    for key in (
        "complete_mixed_root_provider",
        "numerical_Cy_Z3",
        "numerical_mean_Z3",
        "numerical_QZZ",
        "global_five_moment_closure",
        "full_goal_complete",
    ):
        require(packet.get(key) is False, "Unproved completion gate changed: " + key)
    require(packet.get("accepted_ancestors_replayed") is False, "Accepted ancestor replayed")

    require(
        mixed.get("relative_exp_base_one_is_ratio_not_selected_E_value") is True
        and mixed.get("log_absolute_caps_not_signed_derivative_values") is True
        and mixed.get("F0_hb_N_not_materialized") is True,
        "E/a bounds must remain caps with nonmaterialized frequency",
    )
    return dict(
        passed=True,
        ordinary_Z0_Z3=True,
        physical_y0_y1=True,
        exactly_three_new_Xh_rows=True,
        source_provenance_passed=True,
        b_t0_q_p2_open=True,
        Cy_Z3_open=True,
        mean_Z3_open=True,
        QZZ_open=True,
        global_closure_open=True,
        inherited_signed_sources_and_errors_unchanged=True,
        accepted_ancestors_replayed=False,
    )


def run():
    started = time.perf_counter()
    require(sha(ROOT / CODE[0]) == REVIEWED_MODULE_SHA256 and sha(ROOT / HELPER) == HELPER_SHA256,
            "Only the independently reviewed module and helper may produce this report")
    accepted = parent()
    pins(accepted)
    packet = build(ROOT, accepted)
    require(accepted == parent(), "New enclosure mutated its accepted parent")

    inputs = dict(accepted["input_hashes"])
    implementation = dict(accepted["implementation_hashes"])
    for row in packet["source_pins"].values():
        path = row["path"]
        is_input = path.startswith("experiments/root_st073/")
        name = Path(path).name if is_input else path
        target = inputs if is_input else implementation
        require(sha(ROOT / path) == row["expected_sha256"], "Changed new input: " + path)
        require(
            target.get(name, row["expected_sha256"]) == row["expected_sha256"],
            "Conflicting source pin: " + path,
        )
        target[name] = row["expected_sha256"]

    implementation.update({path: sha(ROOT / path) for path in CODE})
    report = copy.deepcopy(accepted)
    report.update(
        stage="actual_micro_Ea_mixed_source_jet_bounds",
        parent_report=Path(PARENT).name,
        parent_report_sha256=PARENT_SHA256,
        actual_micro_Ea_mixed_bounds=packet,
        input_hashes=inputs,
        implementation_hashes=implementation,
    )
    report["contract"] = dict(
        actual_ordinary_source_jet_arithmetic=True,
        ordinary_Z0_to_Z3_and_physical_y0_to_y1_Ea_bounds=True,
        exactly_three_new_Xh_rows=True,
        b_t0_q_p2_Cy_Z3_mean_Z3_QZZ_global_closure_open=True,
        old_signed_sources_and_error_caps_retained=True,
        accepted_ancestors_replayed=False,
    )
    report["checks"] = dict(
        passed=True,
        packet_invariants=_packet_fields(report, accepted),
        independent_new_math_review_completed=True,
        independent_review_completed=True,
        independent_review_scope=(
            "Actual ordinary source jet arithmetic and physical E/a partial provider; "
            "b/t0/q/p2, Cy_Z3, mean_Z3, QZZ, and global closure remain open"
        ),
        independently_reviewed_module_sha256=sha(ROOT / CODE[0]),
        independently_reviewed_model_sha256=sha(ROOT / HELPER),
        independently_reviewed_helper_sha256=sha(ROOT / HELPER),
        independent_runtime_worker_model_metadata_available=False,
        accepted_ancestors_replayed=False,
    )
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    return report


def inspect(path):
    report = _load(path)
    accepted = parent()
    checks = report["checks"]
    require(
        checks["passed"]
        and checks.get("independent_new_math_review_completed") is True
        and checks.get("independent_review_completed") is True
        and checks["independently_reviewed_module_sha256"] == report["implementation_hashes"][CODE[0]] == REVIEWED_MODULE_SHA256
        and checks["independently_reviewed_model_sha256"] == report["implementation_hashes"][HELPER],
        "Same independently reviewed implementation required",
    )
    return dict(
        passed=True,
        pins=pins(report),
        packet_invariants=_packet_fields(report, accepted),
        serialized_only=True,
        accepted_owner_norm_recurrence_or_quadrature_replayed=False,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report", action="store_true")
    group.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    path = DIRECTORY / NAME
    if args.write_report:
        report = run()
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(
            json.dumps(report, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        temporary.replace(path)
        print(
            json.dumps(
                dict(
                    passed=True,
                    execution_seconds=report["execution_seconds"],
                    bytes=path.stat().st_size,
                    sha256=sha(path),
                )
            )
        )
    else:
        print(json.dumps(inspect(path), indent=2))
