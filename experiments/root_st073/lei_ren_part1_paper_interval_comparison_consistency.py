"""Independent consistency checks for a fresh explicit-center core state.

The state JSON is read once, so an atomic checkpoint update occurring later
cannot mix rows from different radial orders.  The checks are deliberately
independent of the interval comparison adapter: they verify the physical
pressure-row identity and build a fresh scalar ``Z=.5`` core from the same
accepted amplitude and pressure datum.

This receipt does not certify an infinite radial tail, an ODE discretization,
whole-axis transition, or terminal moment closure.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import mpmath as mp

from lei_ren_part1_paper_candidate_gauge_core import (
    _unpack_fixed,
    _unpack_rows,
)
from lei_ren_part1_paper_candidate_general_center_factory import (
    GeneralCenterCoreFactory,
    PARAMETERS,
)
from lei_ren_part1_paper_candidate_pressure_function import ACCEPTED_SHA
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


HERE = Path(__file__).resolve().parent
STATE_NAME = "lei_ren_part1_paper_candidate_interval_core_Z049_Z051_state.json"
OUTPUT_NAME = "lei_ren_part1_paper_interval_comparison_consistency.json"


def _path_for(name: str) -> Path:
    path = Path(name)
    return path if path.is_absolute() else HERE / path


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _overlap(left: Any, right: Any) -> bool:
    ll, lh = endpoints(left)
    rl, rh = endpoints(right)
    return max(ll, rl) <= min(lh, rh)


def _contains(family: Any, scalar: Any) -> bool:
    fl, fh = endpoints(family)
    sl, sh = endpoints(scalar)
    return fl <= sl <= sh <= fh


def _coefficient_count(rows: dict[str, list[list[Any]]]) -> int:
    return sum(len(row) for rowset in rows.values() for row in rowset)


def _pressure_expected_row(ctx, squared, a_rows, radial_index, length):
    """Build the directed interval row expected for physical P[radial_index]."""

    if radial_index < 1:
        raise ValueError("generated pressure rows start at radial index one")
    zero = ctx.mpf(0)
    # P_n = S * sum(A_i A_(n-1-i)) / n.  Every operation below is an
    # interval operation; no midpoint or exact-equality assumption is used.
    result = []
    for degree in range(length):
        value = zero
        for s_degree in range(degree + 1):
            convolution_degree = degree - s_degree
            inner = zero
            for i in range(radial_index):
                left = a_rows[i]
                right = a_rows[radial_index - 1 - i]
                for left_degree in range(convolution_degree + 1):
                    right_degree = convolution_degree - left_degree
                    if left_degree < len(left) and right_degree < len(right):
                        inner += left[left_degree] * right[right_degree]
            value += squared[s_degree] * inner
        result.append(value / radial_index)
    return result


def _check_pressure_identity(ctx, fixed, rows):
    pressure = rows["P"]
    p0 = fixed["P0_Z_taylor"]
    relation_count = 0
    overlap_count = 0
    max_gap = ctx.mpf(0)

    if len(pressure[0]) != len(p0):
        raise AssertionError("P[0] and fixed P0 axial lengths differ")
    for stored, expected in zip(pressure[0], p0):
        if not _overlap(stored, expected):
            raise AssertionError("P[0] does not overlap fixed physical P0")
        overlap_count += 1
        relation_count += 1

    squared = fixed["S_Z_taylor"]
    for radial_index in range(1, len(pressure)):
        stored_row = pressure[radial_index]
        expected_row = _pressure_expected_row(
            ctx,
            squared,
            rows["A"],
            radial_index,
            len(stored_row),
        )
        for stored, expected in zip(stored_row, expected_row):
            if not _overlap(stored, expected):
                sl, sh = endpoints(stored)
                el, eh = endpoints(expected)
                gap = max(sl, el) - min(sh, eh)
                if gap > max_gap:
                    max_gap = gap
                raise AssertionError(
                    "physical pressure row does not overlap S*A*A/n: "
                    f"radial index {radial_index}"
                )
            overlap_count += 1
            relation_count += 1
    return {
        "relation_coefficient_count": relation_count,
        "overlap_coefficient_count": overlap_count,
        "maximum_nonoverlap_gap": mp.nstr(max_gap, 40),
        "P0_fixed_identity_checked": True,
        "physical_pressure_row_identity_checked": True,
    }


def _check_scalar_containment(family_rows, scalar_rows):
    counts = {}
    total = 0
    for name in ("A", "Uz", "P"):
        family_set = family_rows[name]
        scalar_set = scalar_rows[name]
        if len(family_set) != len(scalar_set):
            raise AssertionError(f"scalar/family row count differs for {name}")
        count = 0
        for family_row, scalar_row in zip(family_set, scalar_set):
            if len(family_row) != len(scalar_row):
                raise AssertionError(f"scalar/family axial length differs for {name}")
            for family_value, scalar_value in zip(family_row, scalar_row):
                if not _contains(family_value, scalar_value):
                    raise AssertionError(
                        f"scalar Z=.5 core escaped family interval in {name}"
                    )
                count += 1
        counts[name] = count
        total += count
    return counts, total


def run(state_name: str = STATE_NAME, output_name: str = OUTPUT_NAME):
    state_path = _path_for(state_name)
    # Capture bytes exactly once.  Later atomic checkpoint updates are not
    # allowed to change the state being audited.
    state_bytes = state_path.read_bytes()
    state_hash = hashlib.sha256(state_bytes).hexdigest()
    state = json.loads(state_bytes)
    identity = state["identity"]
    degree = int(state["completed_radial_order"])
    target_degree = int(identity["target_radial_degree"])
    initial_length = int(identity["initial_axis_length"])
    if identity["parameters"] != PARAMETERS:
        raise ValueError("captured state parameters differ from Lambda120 datum")
    if identity["accepted_schedule_sha256"] != ACCEPTED_SHA:
        raise ValueError("captured state accepted schedule differs")
    if initial_length != target_degree + 4:
        raise ValueError("captured state does not use target degree + 4 axial length")
    if not 0 <= degree <= target_degree:
        raise ValueError("captured radial order is outside target range")

    input_hashes = state.get("input_hashes", {})
    if not input_hashes:
        raise ValueError("captured state lacks input hashes")
    for name, digest in input_hashes.items():
        path = _path_for(name)
        if _hash(path) != digest:
            raise ValueError("captured state dependency changed: " + name)

    factory = GeneralCenterCoreFactory(precision=int(identity["precision"]))
    ctx = factory.ctx
    with mp.workdps(factory.precision + 40):
        fixed = _unpack_fixed(ctx, state["fixed"])
        rows = {
            name: _unpack_rows(ctx, state["rows"][name])
            for name in ("A", "Uz", "P")
        }
        for name, rowset in rows.items():
            if len(rowset) != degree + 1:
                raise ValueError(f"captured {name} row count is inconsistent")
            for radial_index, row in enumerate(rowset):
                expected_length = initial_length - radial_index
                if len(row) != expected_length:
                    raise ValueError(
                        f"captured {name}[{radial_index}] has wrong axial length"
                    )

        pressure_check = _check_pressure_identity(ctx, fixed, rows)

        # Recompute a scalar core at .5 from the same exact amplitude and
        # accepted pressure datum.  The required depth makes its axis length
        # exactly the captured state's 128 coefficients.
        scalar_required_depth = target_degree + 3 - degree
        if degree + scalar_required_depth + 1 != initial_length:
            raise AssertionError("scalar comparison axis length is not 128")
        scalar = factory.build(
            ".5",
            degree=degree,
            required_depth=scalar_required_depth,
        )
        if scalar["pressure"]["accepted_schedule_sha256"] != ACCEPTED_SHA:
            raise AssertionError("scalar comparison used a different pressure datum")
        if not scalar["same_accepted_pressure_datum"]:
            raise AssertionError("scalar comparison did not preserve accepted pressure")
        containment_counts, containment_total = _check_scalar_containment(
            rows, scalar["core_rows"]
        )

        source_names = (
            Path(__file__).name,
            "lei_ren_part1_paper_candidate_general_center_factory.py",
            "lei_ren_part1_paper_candidate_exact_amplitude.py",
            "lei_ren_part1_paper_candidate_pressure_function.py",
            "lei_ren_part1_paper_uniform_axis_jets.py",
            "lei_ren_part1_paper_factored_core_positivity.py",
            "lei_ren_part1_paper_functional_core_recursion.py",
            "lei_ren_part1_paper_functional_core_step.py",
        )
        source_hashes = {name: _hash(HERE / name) for name in source_names}
        report = {
            "schema_version": 1,
            "all_checks_passed": True,
            "captured_state_file": state_path.name,
            "captured_state_sha256": state_hash,
            "captured_completed_radial_order": degree,
            "target_radial_degree": target_degree,
            "captured_center_family": identity["center"],
            "initial_axis_length": initial_length,
            "scalar_comparison_center": ".5",
            "scalar_comparison_required_depth": scalar_required_depth,
            "scalar_comparison_axis_length": degree + scalar_required_depth + 1,
            "accepted_schedule_sha256": ACCEPTED_SHA,
            "pressure_identity": pressure_check,
            "scalar_family_containment": {
                "counts_by_component": containment_counts,
                "total_containment_count": containment_total,
                "scalar_recomputed_from_same_amplitude_and_pressure": True,
                "old_Z03_tensor_used": False,
            },
            "source_hashes": source_hashes,
            "captured_state_input_hashes": input_hashes,
            "scalar_center_recomputed": True,
            "comparison_driver_checked": False,
            "analytic_tail_target_checked": False,
            "analytic_tail_target_met": False,
            "ODE_discretization_error_enclosed": False,
            "whole_axis_transition_generated": False,
            "terminal_five_moment_closure": False,
            "temporal_recursion": False,
            "scope": (
                "Fresh-state physical pressure identity and scalar .5 family "
                "containment only; no infinite radial, ODE, whole-axis, or "
                "terminal closure claim."
            ),
        }

    output_path = _path_for(output_name)
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        "Interval core consistency passed at captured degree",
        degree,
        ": pressure coefficients",
        pressure_check["overlap_coefficient_count"],
        "; scalar containment",
        containment_total,
        flush=True,
    )
    return report


if __name__ == "__main__":
    run()
