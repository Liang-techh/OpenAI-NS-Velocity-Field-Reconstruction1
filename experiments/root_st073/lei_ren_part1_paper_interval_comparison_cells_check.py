"""Bounded structural check for the directed comparison-cell receipt."""

import hashlib
import json
from pathlib import Path

import mpmath as mp


HERE = Path(__file__).resolve().parent
RECEIPT = HERE / "lei_ren_part1_paper_interval_comparison_cells.json"
SOURCE = HERE / "lei_ren_part1_paper_interval_comparison_cells.py"


def interval(value):
    return mp.mpf(value["lower"]), mp.mpf(value["upper"])


def run():
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    cells = receipt["cells"]
    n = int(receipt["cells_per_half"])
    if len(cells) != 2 * n:
        raise AssertionError("unexpected chronological cell count")
    if [row["phase"] for row in cells[:n]] != ["core"] * n:
        raise AssertionError("first half is not direct core boxes")
    if [row["phase"] for row in cells[n:]] != ["switched"] * n:
        raise AssertionError("second half is not switched boxes")
    if receipt["endpoint"]["phase"] != "switched_endpoint":
        raise AssertionError("endpoint packet is not separate")
    previous_lower = None
    field_names = (
        "phi",
        "U",
        "theta",
        "z",
        "theta_z",
        "p",
        "u_squared",
        "weighted_phi_squared",
    )
    for row in cells:
        ylo, yhi = interval(row["y"])
        if ylo > yhi or (previous_lower is not None and ylo < previous_lower):
            raise AssertionError("cell chronology is not directed")
        previous_lower = ylo
        state = row["state"]
        if any(len(state[name]) != 3 for name in field_names):
            raise AssertionError("state field is not axial order 2")
        if len(row["D"]) != 2 or len(row["I_z"]) != 2:
            raise AssertionError("driver packet is not axial order 1")
    if mp.mpf(receipt["minimum_switch_source_phi_lower"]["value"]) <= 0:
        raise AssertionError("invalid source slope lower bound")
    expected_hash = receipt["input_hashes"][SOURCE.name]
    actual_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    if expected_hash != actual_hash:
        raise AssertionError("receipt source hash mismatch")
    report = {
        "passed": True,
        "cells": len(cells),
        "core_cells": n,
        "switched_cells": n,
        "endpoint_separate": True,
        "field_axial_order": 2,
        "driver_axial_order": 1,
        "source_lower_positive": True,
        "receipt_source_hash_checked": True,
        "physical_exit_ODE_solved": receipt["physical_exit_ODE_solved"],
        "stress_cone_certified": receipt["stress_cone_certified"],
    }
    (HERE / "lei_ren_part1_paper_interval_comparison_cells_check.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print("Directed comparison-cell receipt check passed:", len(cells), "cells", flush=True)
    return report


if __name__ == "__main__":
    run()
