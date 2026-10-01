"""Advance the separate Lambda120 core, then regenerate dependent receipts.

Repeated bounded invocations resume the same exact state. Downstream errors
and moments are produced only after all 124 radial orders are available.
"""
import argparse
import hashlib
import json
from pathlib import Path

import lei_ren_part1_paper_candidate_gauge_core_Lambda120 as driver
import lei_ren_part1_paper_candidate_combined_core_budget as combined
import lei_ren_part1_paper_candidate_inlet_trace as angular
import lei_ren_part1_paper_candidate_shared_inlet as shared

HERE = Path(__file__).resolve().parent
TAIL = "lei_ren_part1_paper_candidate_core_tail_budget_Lambda120.json"
BUDGET = "lei_ren_part1_paper_candidate_combined_core_budget_Lambda120.json"
TRACE = "lei_ren_part1_paper_candidate_inlet_trace_Lambda120.json"
MOMENTS = "lei_ren_part1_paper_candidate_shared_inlet_Lambda120.json"
SIGN = "lei_ren_part1_paper_weighted_angular_derivative.json"


def run(seconds=180):
    sign_path = HERE/SIGN
    sign = json.loads(sign_path.read_text(encoding="utf-8"))
    if not sign["Lambda120_weighted_sign_certified"]:
        raise AssertionError("Lambda120 analytic weighted sign estimate missing")
    for name, digest in sign["input_hashes"].items():
        path = (HERE.parent.parent/name if name.startswith("work_paper_cache/") else HERE/name)
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise AssertionError("Weighted sign input changed: " + name)
    state = driver.run(seconds)
    completed = state["completed_radial_order"]
    report = dict(Lambda=driver.LAMBDA, state_file=driver.STATE_PATH.name,
                  state_sha256=hashlib.sha256(driver.STATE_PATH.read_bytes()).hexdigest(),
                  completed_radial_order=completed,
                  required_radial_degree=driver.core.RADIAL_DEGREE,
                  finite_core_complete=completed==driver.core.RADIAL_DEGREE,
                  analytic_fixed_data_weighted_sign_passed=True,
                  analytic_sign_receipt_sha256=hashlib.sha256(sign_path.read_bytes()).hexdigest(),
                  whole_axis_finite_core_certified=False,
                  combined_center_error_budget_passed=False,
                  shared_core_moments_and_inlet_regenerated=False,
                  terminal_functional_five_moment_closure=False,
                  coherent_collar_matched=False, temporal_recursion=False,
                  full_corrected_NS_residual_validated=False)
    if report["finite_core_complete"]:
        budget = combined.run(driver.STATE_PATH.name, TAIL, BUDGET)
        report["combined_center_error_budget_passed"] = budget[
            "all_conditional_normalized_error_budgets_pass"]
        if not report["combined_center_error_budget_passed"]:
            raise AssertionError("Lambda120 combined normalized C3 budget failed")
        angular.run("4", driver.STATE_PATH.name, TRACE)
        shared.run(driver.STATE_PATH.name, "4", BUDGET, TRACE, MOMENTS)
        report["shared_core_moments_and_inlet_regenerated"] = True
        report["dependent_artifact_hashes"] = {
            name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in (TAIL, BUDGET, TRACE, MOMENTS)}
    driver.core._atomic_write_json(Path(__file__).with_suffix(".json"), report)
    print("Lambda120 pipeline", json.dumps(report), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds", type=float, default=180)
    args = parser.parse_args()
    run(args.seconds)
