"""Exact endpoint comparison for the explicit-center incremental recurrence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_candidate_general_center_factory import (
    GeneralCenterCoreFactory,
    PARAMETERS,
)
from lei_ren_part1_paper_functional_core_step import (
    _source_hashes,
    advance_one,
    initial_rows,
)

HERE = Path(__file__).resolve().parent
DEGREE = 4
REQUIRED_DEPTH = 3


def _exact_compare(expected, actual):
    coefficient_count = 0
    endpoint_count = 0
    for name in ("A", "Uz", "P"):
        if len(expected[name]) != len(actual[name]):
            raise AssertionError(f"row count differs for {name}")
        for row_expected, row_actual in zip(expected[name], actual[name]):
            if len(row_expected) != len(row_actual):
                raise AssertionError(f"axial row length differs for {name}")
            for value_expected, value_actual in zip(row_expected, row_actual):
                if not hasattr(value_expected, "_mpi_") or not hasattr(
                    value_actual, "_mpi_"
                ):
                    raise AssertionError("comparison requires directed intervals")
                if value_expected._mpi_ != value_actual._mpi_:
                    raise AssertionError(f"exact recurrence mismatch in {name}")
                coefficient_count += 1
                endpoint_count += 2
    return coefficient_count, endpoint_count


def run():
    factory = GeneralCenterCoreFactory(precision=260)
    cases = {}
    for label, center in (
        ("scalar_Z05", ".5"),
        ("interval_Z049_Z051", [".49", ".51"]),
    ):
        # The factory independently recomputes amplitude, pressure, axis jets,
        # and the complete reference recurrence at this same center.
        built = factory.build(
            center, degree=DEGREE, required_depth=REQUIRED_DEPTH
        )
        fixed = built["fixed"]
        rows = initial_rows(
            factory.ctx,
            fixed,
            built["center_Z"],
            DEGREE,
            required_depth=REQUIRED_DEPTH,
        )
        for n in range(DEGREE):
            advance_one(
                factory.ctx,
                fixed,
                rows,
                n,
                built["center_Z"],
                PARAMETERS["delta"],
            )
        coefficients, endpoints = _exact_compare(built["core_rows"], rows)
        cases[label] = {
            "center": center,
            "radial_degree": DEGREE,
            "required_depth": REQUIRED_DEPTH,
            "exact_coefficient_comparison_count": coefficients,
            "exact_endpoint_comparison_count": endpoints,
            "same_accepted_pressure_datum": built[
                "same_accepted_pressure_datum"
            ],
            "requested_center_recomputed": built[
                "requested_center_recomputed"
            ],
            "old_center_tensor_loaded": built["old_center_tensor_loaded"],
        }

    source_hashes = _source_hashes()
    for name in (
        "lei_ren_part1_paper_candidate_general_center_factory.py",
        "lei_ren_part1_paper_candidate_exact_amplitude.py",
        "lei_ren_part1_paper_candidate_pressure_function.py",
        "lei_ren_part1_paper_uniform_axis_jets.py",
        "lei_ren_part1_paper_factored_core_positivity.py",
    ):
        source_hashes[name] = hashlib.sha256((HERE / name).read_bytes()).hexdigest()

    report = {
        "schema_version": 1,
        "all_checks_passed": True,
        "api": "advance_one(ctx, fixed, rows, n, center, delta)",
        "initializer": "initial_rows(ctx, fixed, center, degree, required_depth=3)",
        "precision": factory.precision,
        "parameters": PARAMETERS,
        "cases": cases,
        "source_hashes": source_hashes,
        "accepted_schedule_sha256": built["pressure"][
            "accepted_schedule_sha256"
        ],
        "all_fourteen_accepted_pressure_stages_retained": True,
        "old_center_tensor_used": False,
        "whole_axis_core_generated": False,
        "terminal_functional_closure": False,
        "equation_provenance": {
            "paper": "Lei--Ren Part I, Section 8.2, equations (8.1)--(8.2)",
            "paper_text_lines": "work_paper_cache/lei_ren_part1.txt:12275-12332",
            "full_reference": "lei_ren_part1_paper_functional_core_recursion.py:coupled_rows",
        },
        "scope": (
            "Exact finite-step comparison at fresh scalar and interval centers; "
            "no infinite-radial remainder, whole-axis core, or terminal matching."
        ),
    }
    output = HERE / "lei_ren_part1_paper_functional_core_step_check.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        "Functional explicit-center step exact comparisons passed:",
        sum(
            item["exact_endpoint_comparison_count"]
            for item in cases.values()
        ),
        "directed endpoints",
        flush=True,
    )
    return report


if __name__ == "__main__":
    with mp.workdps(300):
        run()
