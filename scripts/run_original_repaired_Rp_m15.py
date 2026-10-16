"""Run and inspect the saved original repaired-Rp m15 forward state.

The runner consumes the accepted compact/quiet forward report and the
independent reserved response.  It never replays an ancestor producer and it
keeps signed source functions, numerical enclosures, and unresolved seam
accuracy as separate facts.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import sys
import time

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "experiments" / "root_st073")]
DIRECTORY = ROOT / "experiments" / "root_st073"

from openai_ns_reconstruction.original_core_numeric import PREFIX, _bounds
from openai_ns_reconstruction.original_repaired_Rp_m15 import OriginalRepairedRpM15


NAME = PREFIX + "compliant_current_original_repaired_Rp_m15_numeric.json"
FIXTURE_NAME = "scripts/check_original_reserved_m15_response.py"
RUNNER_NAME = "scripts/run_original_repaired_Rp_m15.py"
SAVED_Z = ("-.3", "0", ".3")
PARTIAL_X = ("1", "1.27", "2")
REDUCTION_PAIRS = ((4, 8), (124, 143), (125, 144), (126, 145), (133, 151), (134, 152))
PARTITION = (
    ("inner+lower", 0, 6, "axis..rminus; core plus cutoff/lower q=0..2"),
    ("active", 6, 207, "rminus..Rc; all original micro/macro/switch/outer charts"),
    ("compact", 207, 211, "Rc..2Rc; dy=dx/x"),
    ("quiet", 211, 212, "2Rc..Rp; dy=L dq"),
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record_is_valid(row):
    return (
        isinstance(row, dict)
        and isinstance(row.get("lower"), str)
        and isinstance(row.get("upper"), str)
        and isinstance(row.get("source_enclosure_exact_mpf_tuples"), list)
        and len(row["source_enclosure_exact_mpf_tuples"]) == 2
    )


def rows_zero(owner, rows):
    return all(_bounds(owner.response.ctx, owner.response.decode(row)) == (0, 0) for row in rows)


def hull(owner, values):
    c, iv = owner.response.ctx, owner.response.iv
    lo = min(_bounds(c, value)[0] for value in values)
    hi = max(_bounds(c, value)[1] for value in values)
    return iv.mpf([lo, hi])


def select_pulse_sample(owner, Z):
    """Mirror target()'s exact point/boundary/whole-Z source selection."""
    response = owner.response
    selected = owner.pulse["samples"]["whole_Z"]
    zl, zh = _bounds(response.ctx, response.iv.mpf(Z))
    for sample in owner.pulse["samples"].values():
        lo, hi = _bounds(response.ctx, response.decode(sample["Z"]))
        if lo == hi and lo <= zl <= zh <= hi:
            selected = sample
            break
    return selected


def check_partition(packet, terms):
    crosswalk = packet["crosswalk"]
    if len(crosswalk) != len(PARTITION):
        raise AssertionError("four source partition rows required")
    all_indices = []
    for got, (name, left, right, scope) in zip(crosswalk, PARTITION):
        expected_indices = list(range(left, right))
        if got["legacy_component"] != name or got["forward_indices"] != expected_indices:
            raise AssertionError(("partition range", name, got["forward_indices"], expected_indices))
        if got["physical_domain"] != scope:
            raise AssertionError(("partition scope", name, got["physical_domain"], scope))
        if got["forward_names"] != [term["name"] for term in terms[left:right]]:
            raise AssertionError(("partition names", name))
        all_indices.extend(got["forward_indices"])
    if all_indices != list(range(212)) or len(set(all_indices)) != 212:
        raise AssertionError("crosswalk does not cover 0..211 exactly once")
    if len({term["name"] for term in terms[:212]}) != 212:
        raise AssertionError("forward sector names are not unique")
    return dict(
        ranges=[[name, left, right] for name, left, right, _ in PARTITION],
        covered_indices=len(all_indices),
        unique_forward_names=len({term["name"] for term in terms[:212]}),
    )


def check_inner_lower_split(packet):
    """Require the source-owned q=1/q=2 domain split beside the joint six rows."""
    split = packet.get("inner_lower_domain_split")
    if not isinstance(split, dict):
        raise AssertionError("inner/lower source domain split missing")
    expected = {
        "common_coordinate": "R(q)=4epsilon*exp(hb*(sc/4)*q)",
        "core_domain": "axis..4epsilon",
        "cutoff_domain": "q=0..1; 4epsilon..Rin",
        "lower_domain": "q=1..2; Rin..rminus",
        "legacy_inner_identity": "core + collar_prefix(q=1)",
        "legacy_lower_identity": "collar_prefix(q=2) - collar_prefix(q=1)",
        "joint_six_sector_identity": "inner+lower = core + collar_prefix(q=2)",
    }
    for key, value in expected.items():
        if split.get(key) != value:
            raise AssertionError(("inner/lower split", key, split.get(key)))
    if not split.get("shared_q1_endpoint_cancels_as_a_source_function"):
        raise AssertionError("q=1 endpoint source cancellation missing")
    if not split.get("four_collar_sectors_do_not_belong_to_lower_alone"):
        raise AssertionError("q=2 collar sectors misassigned to lower")
    if split.get("saved_q1_terms_appended_to_forward"):
        raise AssertionError("saved q=1 evidence was appended to forward terms")
    if not split.get("enclosure_subtraction_not_selected_as_a_signed_source_value"):
        raise AssertionError("q=1 enclosure subtraction was used as a source value")
    bindings = packet.get("source_partition_bindings", {}).get("archived_nonpulse_components", {})
    if split.get("legacy_inner_graph_pair") != bindings.get("inner"):
        raise AssertionError("legacy inner graph pair is not source-owned")
    if split.get("legacy_lower_graph_pair") != bindings.get("lower"):
        raise AssertionError("legacy lower graph pair is not source-owned")
    q1 = split.get("saved_q1_six_sector_enclosures")
    if not isinstance(q1, list) or len(q1) != 6 or any(not isinstance(row, dict) for row in q1):
        raise AssertionError("saved q=1 six-sector evidence missing")
    return dict(
        common_coordinate=split["common_coordinate"],
        core_domain=split["core_domain"],
        cutoff_domain=split["cutoff_domain"],
        lower_domain=split["lower_domain"],
        legacy_inner_identity=split["legacy_inner_identity"],
        legacy_lower_identity=split["legacy_lower_identity"],
        joint_six_sector_identity=split["joint_six_sector_identity"],
        saved_q1_six_sector_count=len(q1),
        shared_q1_endpoint_cancels_as_source_function=True,
        four_collar_sectors_do_not_belong_to_lower_alone=True,
        saved_q1_terms_appended_to_forward=False,
    )


def check_target(owner, packet, Z):
    response = owner.response
    target = packet["independent_target"]
    pulse = select_pulse_sample(owner, Z)
    pulse_rows = pulse["normalized_true_pulse_C1_enclosures"]["m15"]
    expected = [response.record(-response.decode(row)) for row in pulse_rows]
    if target["coefficient_C0_C1"] != expected:
        raise AssertionError(("independent target sign/source", Z))
    if target["actual_function_binding"]["sign"] != -1:
        raise AssertionError("independent target must carry the negative pulse sign")
    if target["selected_source_Z"] != pulse["Z"]:
        raise AssertionError("target source-Z selection changed")
    return dict(
        selected_source_Z=copy.deepcopy(target["selected_source_Z"]),
        negative_pulse_C0_C1=expected,
        target_matches_negative_pulse=True,
        source=target["actual_function_binding"]["source"],
    )


def check_source_correlation(packet):
    comparison = packet["comparison"]["source_correlated_difference"]
    bindings = packet["source_partition_bindings"]
    if comparison["exact_shared_pulse_physical_handles"] != bindings["complete_pulse_physical"]:
        raise AssertionError("source-correlated pulse handles changed")
    if comparison["identity"] != "reserved_once - (-pulse) = -full + pulse = -(inner+lower+preRp)":
        raise AssertionError("source-correlated identity changed")
    if not comparison["independent_forward_212_remain"] or comparison["no_forward_error_sector_removed"] is not True:
        raise AssertionError("source-correlated difference removed forward sectors")
    if comparison["independent_numeric_residual_resolved"] or packet["comparison"]["numerical_seam_residual_resolved"]:
        raise AssertionError("numeric seam was incorrectly declared resolved")
    names = [row["name"] for row in comparison["normalized_subtracted_nonpulse_components"]]
    if names != ["minus_archived_inner", "minus_archived_lower", "minus_archived_preRp"]:
        raise AssertionError(("nonpulse graph crosswalk", names))
    return dict(
        exact_shared_pulse_physical_handles=copy.deepcopy(comparison["exact_shared_pulse_physical_handles"]),
        nonpulse_component_names=names,
        independent_forward_212_remain=True,
        numerical_seam_residual_resolved=False,
        caps_are_bounds_not_actual_values=packet["comparison"]["caps_are_bounds_not_actual_values"],
    )


def check_correlation_width(owner, packet):
    """Compare old cellwise A*beta boxes with the correlated full response."""
    response = owner.response
    full_response = response.evaluate("0", "2")
    integrals = [response.decode(row) for row in full_response["diagnostics"]["partial_integrals"]]
    # The old independent product is the whole-Z cell hull: retain every
    # source cell, rather than conditioning the comparison on Z=0's two
    # adjacent cells.  The new response is the correlated whole-Z enclosure.
    cells = list(response.coefficients["samples"]["whole_Z"]["cells"])
    if not cells:
        raise AssertionError("whole-Z source cells absent")
    old_rows = [[], []]
    for cell in cells:
        a = [response.decode(row) for row in cell["actual_A_rc_C0_C1_coefficient_enclosures"]]
        b = [
            [response.decode(row) for row in cell["coefficients"][name]["normalized_C0_C1_enclosures"]]
            for name in ("beta1", "beta2", "beta3")
        ]
        old_rows[0].append(-sum(a[0] * b[j][0] * integrals[j] for j in range(3)))
        old_rows[1].append(
            -sum((a[1] * b[j][0] + a[0] * b[j][1]) * integrals[j] for j in range(3))
        )
    old_hull = [hull(owner, old_rows[k]) for k in range(2)]
    new_rows = [response.decode(row) for row in full_response["coefficient_C0_C1"]]
    old_widths, new_widths, ratios = [], [], []
    for old, new in zip(old_hull, new_rows):
        olo, ohi = _bounds(response.ctx, old)
        nlo, nhi = _bounds(response.ctx, new)
        old_width = ohi - olo
        new_width = nhi - nlo
        if old_width <= new_width or new_width <= 0:
            raise AssertionError(("correlation width did not improve", old_width, new_width))
        old_widths.append(response.ctx.nstr(old_width, 16))
        new_widths.append(response.ctx.nstr(new_width, 16))
        ratios.append(response.ctx.nstr(old_width / new_width, 16))
    return dict(
        Z="0",
        x="2",
        source_cells=len(cells),
        old_cellwise_box_widths=old_widths,
        new_correlated_response_widths=new_widths,
        width_improvement_ratios=ratios,
        interpretation="bound-width comparison only; no source value or numeric Rp seam closure",
    )


def check_partial_rows(owner):
    """Exercise the actual source response at three compact coordinates."""
    response = owner.response
    reports = []
    for x in PARTIAL_X:
        out = response.evaluate("0", x)
        angular = out["angular_inverse_response"]
        rows = out["coefficient_C0_C1"]
        if len(angular) != 3 or len(rows) != 2 or any(not record_is_valid(row) for row in angular + rows):
            raise AssertionError(("partial C01 serialization", x))
        if x == "1":
            if not out["diagnostics"]["exact_empty"] or out["diagnostics"]["exact_full"] or not rows_zero(owner, rows):
                raise AssertionError("empty compact response is not exact zero")
        elif x == "2":
            expected = [response.record(response.iv.mpf(v)) for v in (0, 0, 1)]
            if angular != expected:
                raise AssertionError("full compact response is not exact e2")
        else:
            if out["diagnostics"]["exact_empty"] or out["diagnostics"]["exact_full"]:
                raise AssertionError("partial compact response was classified as an endpoint")
            if all(response.decode(row) == 0 for row in angular):
                raise AssertionError("partial compact response unexpectedly zero")
        reports.append(dict(
            x=x,
            exact_empty=out["diagnostics"]["exact_empty"],
            exact_full=out["diagnostics"]["exact_full"],
            angular_inverse_response=angular,
            coefficient_C0_C1=rows,
            source_family=out["source_family"],
        ))
    return reports


def check_packet(owner, packet, Z):
    response = owner.response
    terms = packet["m15_C0_C1_terms"]
    saved = owner.saved[Fraction(Z)]["m15_C0_C1_terms"]
    inner = owner.inner[Fraction(Z)]["m15_C0_C1_terms"]
    if len(terms) != 213 or terms[:212] != saved or packet["original_forward_terms"] != 212:
        raise AssertionError(("raw 213/first 212", Z, len(terms)))
    if terms[:6] != inner[:6]:
        raise AssertionError(("raw six inner terms changed", Z))
    if packet["reserved_terms_added"] != 1 or terms[212]["name"] != "reserved_n1_own_swirl_once":
        raise AssertionError(("reserved term count/name", Z))

    partition = check_partition(packet, terms)
    inner_lower_split = check_inner_lower_split(packet)
    reduced = packet["reduced_m15_C0_C1_terms"]
    cancellation_rows = packet["exact_internal_endpoint_cancellations"]
    actual_pairs = tuple(tuple(row["indices"]) for row in cancellation_rows)
    if actual_pairs != REDUCTION_PAIRS:
        raise AssertionError(("reduction pairs", Z, actual_pairs))
    reduction_indices = sorted(
        index for row in cancellation_rows for index in row["indices"]
    )
    if len(reduction_indices) != len(set(reduction_indices)) or not {4, 8}.issubset(reduction_indices):
        raise AssertionError(("reduction indices", Z, reduction_indices))
    expected_reduced = [term for index, term in enumerate(terms) if index not in reduction_indices]
    if reduced != expected_reduced:
        raise AssertionError(("reduced endpoint rows", Z, reduction_indices))
    retained_error_names = [
        term["name"] for term in terms
        if "finite_N" in term["name"] or "width_error" in term["name"]
    ]
    if any(name not in {term["name"] for term in reduced} for name in retained_error_names):
        raise AssertionError(("finite-N/width error sector dropped", Z))

    target = check_target(owner, packet, Z)
    source_correlation = check_source_correlation(packet)
    own = terms[-1]
    full_rows = owner.full["normalized_full_pressure_moment_C1_enclosures"]["m15"]
    expected_own = [response.record(-response.decode(row)) for row in full_rows]
    if packet["reserved_response"]["angular_inverse_response"] != [
        response.record(response.iv.mpf(v)) for v in (0, 0, 1)
    ]:
        raise AssertionError(("full own response is not e2", Z))
    if own["coefficient_C0_C1"] != expected_own:
        raise AssertionError(("own response is not negative full rows", Z))
    for row in own["coefficient_C0_C1"]:
        if not record_is_valid(row):
            raise AssertionError(("own C01 serialization", Z))
    return dict(
        Z=Z,
        raw_terms=len(terms),
        saved_forward_terms=len(saved),
        raw_first_212_exact_saved=True,
        raw_six_exact_inner=True,
        reserved_terms_added=1,
        partition=partition,
        inner_lower_domain_split=inner_lower_split,
        retained_finite_N_width_error_names=retained_error_names,
        endpoint_reduction_pairs=[list(pair) for pair in actual_pairs],
        endpoint_reduction_dropped_indices=reduction_indices,
        reduced_terms=len(reduced),
        independent_target=target,
        source_correlation=source_correlation,
        full_response_e2=True,
        own_equals_negative_full_C01=True,
    )


def run():
    started = time.perf_counter()
    # Run exactly once in production; the result is retained as evidence in
    # the report and is not recomputed by inspect().
    from check_original_reserved_m15_response import run as run_fixture

    fixture = run_fixture()
    if not fixture["passed"]:
        raise AssertionError("independent reserved-response fixture failed")

    owner = OriginalRepairedRpM15(precision=100, cells=128)
    packets = [owner.evaluate(Z) for Z in SAVED_Z]
    packet_checks = [check_packet(owner, packet, Z) for packet, Z in zip(packets, SAVED_Z)]
    partial_rows = check_partial_rows(owner)
    correlation_width = check_correlation_width(owner, packets[1])

    rejected = 0
    for bad in (".2", "nan", "bad"):
        try:
            owner.evaluate(bad)
        except (ValueError, ArithmeticError):
            rejected += 1
        else:
            raise AssertionError(("invalid/absent Z admitted", bad))
    if rejected != 3:
        raise AssertionError("invalid/absent Z rejection count changed")

    hashes = dict(owner.hashes)
    code = dict(owner.code_hashes)
    code[FIXTURE_NAME] = sha(ROOT / FIXTURE_NAME)
    code[RUNNER_NAME] = sha(ROOT / RUNNER_NAME)
    report = dict(
        source_family=owner.family,
        samples=packets,
        partial_rows=partial_rows,
        input_hashes=hashes,
        implementation_hashes=code,
        checks=dict(
            passed=True,
            independent_fixture=fixture,
            saved_Z=SAVED_Z,
            actual_repaired_Rp_packets=len(packets),
            raw_213_first_212_exact_saved=True,
            partition_ranges=[list(item[:3]) for item in PARTITION],
            inner_lower_domain_split_valid=True,
            raw_six_unchanged=True,
            endpoint_reduction_pairs=[list(pair) for pair in REDUCTION_PAIRS],
            endpoint_reduction_dropped_indices=sorted(index for pair in REDUCTION_PAIRS for index in pair),
            reduced_terms_expected=201,
            all_finite_N_and_width_errors_retained=True,
            independent_target_negative_actual_pulse_C01=True,
            source_correlated_graph_handles=True,
            full_own_response_e2_and_negative_full_C01=True,
            actual_partial_row_queries=len(partial_rows),
            rejected_absent_or_malformed_Z=rejected,
            correlation_width_comparison=correlation_width,
        ),
        execution_seconds=round(time.perf_counter() - started, 3),
        limits=dict(
            saved_forward_only=True,
            ancestor_constructors_reexecuted=False,
            numeric_Rp_seam_closed=False,
            caps_equal_source_values=False,
            signed_finite_N_cancellation_resolved=False,
            full_corrected_NS=False,
        ),
    )
    return report


def inspect(path):
    """Verify report hashes and serialized invariants without recalculation."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if not data["checks"]["passed"] or data["source_family"] is None:
        raise AssertionError("saved repaired-Rp report is not passed")
    for base, group in ((DIRECTORY, data["input_hashes"]), (ROOT, data["implementation_hashes"])):
        for name, digest in group.items():
            if sha(base / name) != digest:
                raise AssertionError(("hash changed", name))
    samples = data["samples"]
    if len(samples) != 3 or tuple(sample["Z"] for sample in samples) != SAVED_Z:
        raise AssertionError("saved Z packet set changed")
    for sample in samples:
        terms = sample["m15_C0_C1_terms"]
        reduced = sample["reduced_m15_C0_C1_terms"]
        if len(terms) != 213:
            raise AssertionError(("serialized packet size", sample["Z"]))
        ranges = [tuple(row["forward_indices"]) for row in sample["crosswalk"]]
        if ranges != [tuple(range(left, right)) for _, left, right, _ in PARTITION]:
            raise AssertionError(("serialized partition ranges", sample["Z"]))
        check_inner_lower_split(sample)
        cancellation_rows = sample["exact_internal_endpoint_cancellations"]
        actual_pairs = tuple(tuple(row["indices"]) for row in cancellation_rows)
        if actual_pairs != REDUCTION_PAIRS:
            raise AssertionError(("serialized reduction pairs", sample["Z"], actual_pairs))
        reduction_indices = sorted(index for row in cancellation_rows for index in row["indices"])
        if len(reduction_indices) != len(set(reduction_indices)) or not {4, 8}.issubset(reduction_indices):
            raise AssertionError(("serialized reduction indices", sample["Z"]))
        if len(reduced) != 201:
            raise AssertionError(("serialized reduced packet size", sample["Z"], len(reduced)))
        if [row["name"] for row in reduced] != [
            term["name"] for index, term in enumerate(terms) if index not in reduction_indices
        ]:
            raise AssertionError(("serialized endpoint reduction", sample["Z"]))
        angular = sample["reserved_response"]["angular_inverse_response"]
        zero = [[0, 0, 0, 0], [0, 0, 0, 0]]
        one = [[0, 1, 0, 1], [0, 1, 0, 1]]
        if angular[0]["source_enclosure_exact_mpf_tuples"] != zero:
            raise AssertionError(("serialized full response first row", sample["Z"]))
        if angular[1]["source_enclosure_exact_mpf_tuples"] != zero:
            raise AssertionError(("serialized full response second row", sample["Z"]))
        if angular[2]["source_enclosure_exact_mpf_tuples"] != one:
            raise AssertionError(("serialized full response third row", sample["Z"]))
        if sample["m15_C0_C1_terms"][-1]["coefficient_C0_C1"] != sample["reserved_response"]["coefficient_C0_C1"]:
            raise AssertionError(("serialized own/response C01", sample["Z"]))
        if not sample["contract"]["saved_forward_212_retained_verbatim"]:
            raise AssertionError(("saved forward contract", sample["Z"]))
    if not data["checks"]["full_own_response_e2_and_negative_full_C01"]:
        raise AssertionError("saved own-response invariant missing")
    if not data["limits"]["saved_forward_only"] or data["limits"]["numeric_Rp_seam_closed"]:
        raise AssertionError("saved report overclaims production scope")
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-report", action="store_true")
    parser.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    path = DIRECTORY / NAME
    data = inspect(path) if args.inspect_report else run()
    if args.write_report:
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8", newline="\n")
        temporary.replace(path)
    print(json.dumps(dict(checks=data["checks"], execution_seconds=data["execution_seconds"]), indent=2), flush=True)
