"""Produce the coupled signed macro packet; inspect it without replay."""

from __future__ import annotations

import argparse
import json
from collections.abc import Mapping
from pathlib import Path
import sys
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from openai_ns_reconstruction.original_macro_coupled_m15 import (  # noqa: E402
    PARENT,
    PARENT_SHA256,
    OriginalMacroCoupledM15,
    digest,
    load,
    materialize_coupled_view,
    sha,
)


DIRECTORY = ROOT / "experiments/root_st073"
NAME = "lei_ren_part1_paper_compliant_current_original_macro_coupled_m15_numeric.json"
CODE = (
    "src/openai_ns_reconstruction/original_native_macro_seam.py",
    "src/openai_ns_reconstruction/original_macro_defect_kernel.py",
    "src/openai_ns_reconstruction/original_macro_coupled_m15.py",
    "scripts/check_original_macro_defect_kernel.py",
    "scripts/run_original_macro_coupled_m15.py",
)
EXPECTED_Z = ("-.3", "0", ".3")
EXPECTED_REDUCED_INDICES = [index for index in range(183) if index != 39]
EXPECTED_REPLACEMENTS = [19, 20, 41, 42]


def require(value: Any, message: str) -> None:
    if not value:
        raise AssertionError(message)


def parent() -> Mapping[str, Any]:
    """Load the accepted weighted parent and verify its immutable bytes."""

    require(sha(DIRECTORY / PARENT) == PARENT_SHA256, "Pinned weighted parent bytes changed")
    report = load(DIRECTORY / PARENT)
    require(isinstance(report, Mapping), "Weighted parent mapping required")
    require(report.get("checks", {}).get("passed") is True, "Passed weighted parent required")
    return report


def _assert_unselected(value: Any, path: str = "source") -> None:
    """Reject interval approximations promoted to exact source values."""

    if isinstance(value, Mapping):
        if "approximation_is_exact_source_value" in value:
            require(
                value["approximation_is_exact_source_value"] is False,
                f"Promoted approximation: {path}",
            )
        if "source_enclosure_is_not_a_selected_value" in value:
            require(
                value["source_enclosure_is_not_a_selected_value"] is True,
                f"Selected enclosure: {path}",
            )
        for key, item in value.items():
            _assert_unselected(item, f"{path}.{key}")
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _assert_unselected(item, f"{path}[{index}]")


def _false_if_present(mapping: Mapping[str, Any], *keys: str) -> None:
    for key in keys:
        if key in mapping:
            require(mapping[key] is False, f"Forbidden true witness: {key}")


def _check_fixture(fixture: Mapping[str, Any]) -> dict[str, Any]:
    """Check the new defect fixture without invoking any producer or owner."""

    require(fixture.get("passed") is True, "Independent defect fixture failed")
    _false_if_present(
        fixture,
        "source_loaders_or_owners_called",
        "accepted_owner_or_producer_called",
        "quadrature_in_production",
        "width_materialized",
        "source_values_selected",
    )
    controls = fixture.get("controls")
    if controls is not None:
        require(isinstance(controls, list), "Defect fixture controls must be a list")
        require(all(isinstance(row, Mapping) and row.get("passed") is True for row in controls),
                "Defect fixture control failed")
    return {"passed": True, "control_count": len(controls) if isinstance(controls, list) else None}


def _check_source_seam(seam: Mapping[str, Any]) -> None:
    """Check stable source-proof gates while allowing proof detail to grow."""

    require(seam.get("passed") is True, "Native/macro source seam failed")
    ast_rows = seam.get("AST")
    require(isinstance(ast_rows, Mapping) and bool(ast_rows), "Native/macro AST source proof missing")
    require(all(isinstance(row, Mapping) and row.get("passed") is True for row in ast_rows.values()),
            "Native/macro AST source proof failed")
    for name in ("same_physical_R", "nonzero_epsilon_cancellation_formal",
                 "normalized_ordinary_Z_orders_0_to_2", "derivative_port_defect_identity"):
        require(isinstance(seam.get(name), Mapping) and seam[name].get("passed") is True,
                f"Required source theorem failed: {name}")
    require(seam["normalized_ordinary_Z_orders_0_to_2"].get("ordinary_Z_orders") == [0,1,2],
            "Normalized value seam ordinary orders changed")
    for key in ("input_hashes", "implementation_hashes"):
        require(isinstance(seam.get(key), Mapping), f"Native/macro seam {key} missing")

    _false_if_present(
        seam,
        "source_values_selected",
        "numeric_caps_used",
        "individual_caps_as_values",
        "individual_positive_caps_used_as_signed_values",
        "native_Uy_ports_equated",
        "native_lower_equals_macro_source_Uy",
        "native_lower_equals_macro0_Uy_proved",
        "actual_port_value_equality_proved",
        "accepted_owner_or_producer_called",
        "numerical_quadrature_called",
        "source_executed",
    )
    if "native_Uy_and_macro_source_Uy_are_distinct_ports" in seam:
        require(seam["native_Uy_and_macro_source_Uy_are_distinct_ports"] is True,
                "Native and macro Uy ports were conflated")

    # The seam may carry these gates in theorem sub-packets rather than at
    # the top level.  Preserve their explicit false semantics wherever they
    # are present, without depending on proof-detail field names.
    for name in ("same_physical_R", "nonzero_epsilon_cancellation_formal",
                 "normalized_ordinary_Z_orders_0_to_2", "derivative_port_defect_identity",
                 "derivative_port_defect"):
        row = seam.get(name)
        if isinstance(row, Mapping):
            _false_if_present(
                row,
                "source_values_selected",
                "numeric_caps_used",
                "individual_caps_as_values",
                "native_lower_equals_macro_source_Uy",
                "native_lower_equals_macro0_Uy_proved",
                "actual_port_value_equality_proved",
                "native_Uy_ports_equated",
            )


def _check_binding(binding: Mapping[str, Any], seam: Mapping[str, Any]) -> None:
    require(binding.get("passed") is True, "Signed macro ledger binding failed")
    require(binding.get("old_raw_indices") == [23, 24, 44], "Whole endpoint raw scope changed")
    require(binding.get("old_reduced_indices") == [19, 20, 39], "Whole endpoint reduced scope changed")
    require("whole_signed_endpoint_binding_digest" not in binding or
            isinstance(binding["whole_signed_endpoint_binding_digest"], str),
            "Malformed whole endpoint binding digest")
    _false_if_present(
        binding,
        "individual_positive_caps_used_as_signed_values",
        "individual_caps_as_values",
        "positive_caps_used_as_values",
        "native_macro_derivative_ports_equated",
        "native_Uy_ports_equated",
        "source_values_selected",
        "numeric_caps_used",
    )
    require(binding.get("actual_lower_radius_rho2_over4_retained") is True,
            "Native lower radius factor was dropped")
    require(binding.get("new_signed_function"), "Signed whole endpoint function identity missing")
    require(binding.get("independent_leading_FTC_passed") is True, "Physical signed-density FTC missing")
    require(binding.get("signed_endpoint_orientation") == {"native_lower":1,"macro_upper":-1,"first_switch_lower_E":1},
            "Signed whole endpoint orientation changed")
    require(all(v.get("passed") is True for v in binding["AST"].values()), "Endpoint category source guard failed")
    require(binding.get("actual_minus_proxy_category_reconstructed_from_original_physical_FTC") is True,
            "Unspecified cap used as the source category")
    # The source proof is the seam witness, not a numerical equality of the
    # two derivative ports.
    if "native_macro_derivative_ports_equated" in binding:
        require(binding["native_macro_derivative_ports_equated"] is False,
                "Unproved native/macro derivative ports equated")
    require(isinstance(seam, Mapping), "Source seam mapping required")


def _ledger_samples(parent_report: Mapping[str, Any]) -> tuple[Mapping[str, Any], dict[str, Mapping[str, Any]]]:
    packet = parent_report.get("inherited_m15_packet")
    require(isinstance(packet, Mapping), "Inherited m15 packet missing")
    ledger_name = packet.get("report")
    ledger_sha = packet.get("report_sha256")
    require(isinstance(ledger_name, str) and isinstance(ledger_sha, str),
            "Inherited ledger identity missing")
    ledger_path = DIRECTORY / ledger_name
    require(sha(ledger_path) == ledger_sha, "Pinned original ledger bytes changed")
    require(parent_report.get("input_hashes", {}).get(ledger_name) == ledger_sha,
            "Ledger input pin changed")
    ledger = load(ledger_path)
    samples = ledger.get("samples")
    require(isinstance(samples, list), "Original ledger samples missing")
    by_z = {sample.get("Z"): sample for sample in samples if isinstance(sample, Mapping)}
    require(tuple(by_z) == EXPECTED_Z, "Original ledger sample order changed")
    return packet, by_z


def _raw_indices(value: Mapping[str, Any]) -> set[int]:
    found: set[int] = set()
    for key in ("raw_index", "raw_indices", "source_raw_index", "source_raw_indices",
                "old_raw_index", "old_raw_indices"):
        item = value.get(key)
        if isinstance(item, int):
            found.add(item)
        elif isinstance(item, (list, tuple)):
            found.update(index for index in item if isinstance(index, int))
    return found


def _check_defect_rows(report: Mapping[str, Any], parent_report: Mapping[str, Any]) -> None:
    rows = report.get("native_macro_defect_rows")
    source_rows = parent_report.get("macro_weighted_source_rows")
    require(isinstance(rows, list) and len(rows) == 3, "Three native/macro defect rows required")
    require(isinstance(source_rows, list), "Parent source rows missing")
    source_by_z = {row.get("Z"): row for row in source_rows if isinstance(row, Mapping)}
    require(tuple(row.get("Z") for row in rows) == EXPECTED_Z, "Defect sample order changed")
    for row in rows:
        require(isinstance(row, Mapping), "Malformed defect row")
        z = row.get("Z")
        require(row.get("ordinary_Z_orders") == [0, 1], f"Defect ordinary Z orders changed: {z}")
        for key in ("defect_C0_C1", "macro0_endpoint_C0_C1", "native_endpoint_C0_C1",
                    "correlated_result_C0_C1"):
            require(isinstance(row.get(key), list) and len(row[key]) == 2,
                    f"Defect row {z} missing {key}")
        expected = source_by_z[z]["native_lower_endpoint_C0_C1"]
        require(row["native_endpoint_C0_C1"] == expected, f"Native endpoint changed: {z}")
        require(row["correlated_result_C0_C1"] == expected,
                f"Correlated defect result was recomputed or changed: {z}")
        require(row.get("source_ports_distinct") is True, f"Defect ports conflated: {z}")
        require(row.get("independent_interval_sum_used_as_identity") is False,
                f"Independent interval sum promoted to identity: {z}")
        require(row.get("defect_added_as_an_extra_ledger_term") is False,
                f"Defect appended as an extra ledger term: {z}")
        _assert_unselected(row, f"defect.{z}")


def _check_new_term(
    entry: Mapping[str, Any],
    reduced: list[Mapping[str, Any]],
    source: Mapping[str, Any],
    seam: Mapping[str, Any],
    binding: Mapping[str, Any],
    index: int,
) -> None:
    require(entry.get("operation") == "same_source_signed_leaf_replacement",
            f"Coupled replacement operation missing at {index}")
    require(entry.get("parent_reduced_index") == index, f"Coupled parent index changed at {index}")
    require(entry.get("parent_term_digest") == digest(reduced[index]),
            f"Coupled parent digest changed at {index}")
    term = entry.get("term")
    require(isinstance(term, Mapping), f"Coupled replacement term missing at {index}")
    require(term.get("log_scale_terms") == reduced[index].get("log_scale_terms"),
            f"Physical log-scale atoms changed at {index}")
    expected_coefficients = source["weighted_C0_C1"] if index == 19 else source["native_lower_endpoint_C0_C1"]
    require(term.get("coefficient_C0_C1") == expected_coefficients,
            f"Coupled source coefficients changed at {index}")
    require(entry.get("coupled_raw_indices") == [23, 24, 44],
            f"Coupled raw source scope changed at {index}")
    term_binding = term.get("actual_function_binding")
    require(isinstance(term_binding, Mapping), f"Coupled source binding missing at {index}")
    require(term_binding.get("actual_signed_coupled_macro_source") is True,
            f"Actual signed source binding missing at {index}")
    require(term_binding.get("error_enclosure_only") is False,
            f"Positive error cap remains at {index}")
    require(term_binding.get("individual_cap_selected_as_value") is False,
            f"Individual cap selected as value at {index}")
    require(term_binding.get("native_macro_value_seam_digest") == digest(seam),
            f"Native/macro seam digest changed at {index}")
    require(term_binding.get("whole_signed_endpoint_binding_digest") == digest(binding),
            f"Whole endpoint binding digest changed at {index}")
    require(term_binding.get("old_raw_indices") == [23, 24, 44],
            f"Term whole endpoint raw scope changed at {index}")
    require(term_binding.get("signed_function_source_digest") == digest(source),
            f"Term source identity changed at {index}")
    _assert_unselected(term, f"view.term.{index}")


def _check_view(
    view: Mapping[str, Any],
    parent_view: Mapping[str, Any],
    sample: Mapping[str, Any],
    source: Mapping[str, Any],
    seam: Mapping[str, Any],
    binding: Mapping[str, Any],
) -> list[Mapping[str, Any]]:
    raw = sample.get("raw_m15_C0_C1_terms")
    reduced = sample.get("correlated_reduced_m15_C0_C1_terms")
    require(isinstance(raw, list) and len(raw) == 213, f"Raw archive changed at {view.get('Z')}")
    require(isinstance(reduced, list) and len(reduced) == 183, f"Reduced archive changed at {view.get('Z')}")
    require(view.get("parent_original_reduced183_digest") == parent_view.get("archived_reduced183_digest"),
            "Parent reduced archive digest changed")
    require(view.get("raw213_digest") == parent_view.get("raw213_digest"), "Raw 213 digest changed")
    require(view.get("parent_signed183_view_digest") == parent_view.get("expanded_signed_view_digest"),
            "Parent signed 183 digest changed")
    entries = view.get("entries")
    require(isinstance(entries, list) and len(entries) == 182, "Coupled view term count changed")
    require(view.get("view_terms") == 182 and view.get("unchanged_original_references") == 178,
            "Coupled view counts changed")
    require([entry.get("parent_reduced_index") for entry in entries] == EXPECTED_REDUCED_INDICES,
            "Coupled parent reduced indices changed")
    require(view.get("four_signed_replacement_slots") == EXPECTED_REPLACEMENTS,
            "Coupled replacement slots changed")
    require(view.get("contracted_parent_slot") == 39, "Coupled contraction slot changed")

    expanded = materialize_coupled_view(reduced, entries)
    require(view.get("expanded_signed182_digest") == digest(expanded), "Coupled expansion digest changed")
    require(expanded[:6] == reduced[:6], "Six inner/lower terms changed")
    require(view.get("six_inner_lower_digest") == parent_view.get("six_inner_lower_digest"),
            "Six inner/lower witness changed")

    by_index = {entry.get("parent_reduced_index"): entry for entry in entries}
    require(set(by_index) == set(EXPECTED_REDUCED_INDICES), "Coupled view index set changed")
    replacements = [entry for entry in entries if entry.get("operation") != "unchanged_parent_reference"]
    require([entry.get("parent_reduced_index") for entry in replacements] == EXPECTED_REPLACEMENTS,
            "Coupled replacement scope changed")
    for index in EXPECTED_REDUCED_INDICES:
        if index in (19, 20):
            _check_new_term(by_index[index], reduced, source, seam, binding, index)
        elif index in (41, 42):
            require(by_index[index] == parent_view["entries"][index],
                    f"Prior signed raw 46/47 entry changed at {index}")
        else:
            require(by_index[index] == parent_view["entries"][index],
                    f"Original parent reference changed at {index}")

    provenance = view.get("contracted_source_endpoint_provenance")
    require(isinstance(provenance, list) and len(provenance) == 1,
            "Removed raw 44 provenance missing")
    removed = provenance[0]
    require(isinstance(removed, Mapping), "Removed endpoint provenance malformed")
    require(removed.get("parent_reduced_index", removed.get("previous_reduced_index")) == 39,
            "Removed endpoint reduced index changed")
    require(removed.get("parent_term_digest", digest(reduced[39])) == digest(reduced[39]),
            "Removed endpoint parent digest changed")
    require(44 in _raw_indices(removed), "Removed raw 44 provenance missing")
    for key in ("active", "provenance_only", "no_cap_value_or_midpoint"):
        if key in removed:
            require(removed[key] is (False if key == "active" else True),
                    f"Removed provenance flag changed: {key}")

    require(view.get("retained_first_switch_signed_pair_digest") ==
            digest([by_index[41], by_index[42]]), "Prior signed raw 46/47 digest changed")
    require(view.get("independent_target") == parent_view.get("independent_target"),
            "Independent target changed")
    require(view.get("independent_target_digest") == digest(sample["independent_target"]),
            "Independent target digest changed")
    target = view["independent_target"]
    require(target.get("actual_function_binding", {}).get("sign") == -1,
            "Independent target sign changed")
    require(target.get("actual_function_binding", {}).get("physical_measure") ==
            "-R*pressure_density*dR", "Independent target physical measure changed")
    require(view.get("selected_N_evidence_digest") == digest(sample["selected_N_evidence"]),
            "Selected N witness changed")
    require(view.get("largest_unchanged_caps") == sample["largest_remaining_individual_log_caps"],
            "Largest unchanged caps changed")
    require(view.get("global_terminal_closure") is False, "Global closure claim introduced")
    _assert_unselected(view, f"view.{view.get('Z')}")
    return expanded


def check_packet(report: Mapping[str, Any], parent_report: Mapping[str, Any]) -> dict[str, Any]:
    require(report.get("parent_report") == PARENT and
            report.get("parent_report_sha256") == PARENT_SHA256, "Parent identity changed")
    require(report.get("source_family") == parent_report.get("source_family"), "Source family changed")
    seam = report.get("native_macro_seam")
    binding = report.get("signed_ledger_binding")
    require(isinstance(seam, Mapping), "Native/macro source seam missing")
    require(isinstance(binding, Mapping), "Signed ledger binding missing")
    _check_source_seam(seam)
    _check_binding(binding, seam)

    require(report.get("inherited_m15_packet") == parent_report.get("inherited_m15_packet"),
            "Inherited m15 packet changed")
    require(report.get("inherited_m15_packet_digest") == parent_report.get("inherited_m15_packet_digest"),
            "Inherited packet digest changed")
    require(report.get("inherited_signed183_views_digest") == digest(parent_report["signed_m15_query_views"]),
            "Inherited signed view digest changed")
    source_rows = report.get("macro_weighted_source_rows")
    require(source_rows == parent_report.get("macro_weighted_source_rows"),
            "Accepted macro weighted source rows changed")
    require(report.get("macro_weighted_source_rows_digest") == digest(source_rows),
            "Macro weighted source row digest changed")
    _check_defect_rows(report, parent_report)

    _, samples = _ledger_samples(parent_report)
    views = report.get("signed_m15_query_views")
    parent_views = parent_report.get("signed_m15_query_views")
    require(isinstance(views, list) and len(views) == 3, "Three coupled views required")
    require(isinstance(parent_views, list) and len(parent_views) == 3, "Three parent signed views required")
    source_by_z = {row.get("Z"): row for row in source_rows if isinstance(row, Mapping)}
    for view, parent_view in zip(views, parent_views):
        require(view.get("Z") == parent_view.get("Z") and view.get("Z") in samples,
                "Coupled view sample order changed")
        require(view["Z"] in source_by_z, "Coupled source row sample missing")
        _check_view(view, parent_view, samples[view["Z"]], source_by_z[view["Z"]], seam, binding)

    contract = report.get("contract")
    require(isinstance(contract, Mapping), "Coupled contract missing")
    expected_true = (
        "raw213_reduced183_archives_unchanged",
        "signed_raw46_raw47_unchanged",
        "six_inner_lower_reserved_once_target_and_N_unchanged",
        "raw43_raw45_and_all_other_errors_retained",
    )
    expected_false = (
        "accepted_ancestor_or_fixture_replayed",
        "saved_macro_primitive_recomputed",
        "mass_core_M0_W1_W2_recomputed",
        "positive_caps_midpoints_or_interval_overlap_used_as_values",
        "hbB_hbS_or_absolute_swirl_amplitude_materialized",
    )
    for key in expected_true:
        require(contract.get(key) is True, f"Contract witness failed: {key}")
    for key in expected_false:
        require(contract.get(key) is False, f"Contract changed: {key}")
    limits = report.get("limits")
    require(isinstance(limits, Mapping) and all(value is False for value in limits.values()),
            "Global completion claim introduced")
    return {
        "passed": True,
        "defect_rows": 3,
        "ordinary_defect_Z_orders": [0, 1],
        "signed_views": 3,
        "view_terms": 182,
        "unchanged_references_per_view": 178,
        "new_replacement_slots": [19, 20],
        "retained_signed_slots": [41, 42],
        "contracted_parent_slot": 39,
        "raw213_reduced183_archives_unchanged": True,
        "independent_target_reserved_once_six_inner_lower_and_N_retained": True,
        "global_limits_all_false": True,
    }


def pins(report: Mapping[str, Any]) -> int:
    count = 0
    for base, key in ((DIRECTORY, "input_hashes"), (ROOT, "implementation_hashes")):
        values = report.get(key)
        require(isinstance(values, Mapping), f"{key} mapping missing")
        for name, expected in values.items():
            require(isinstance(name, str) and isinstance(expected, str), f"Malformed {key} pin")
            path = base / name
            require(sha(path) == expected, "Changed pin: " + name)
            count += 1
    return count


def run() -> dict[str, Any]:
    """Run only the new defect fixture and the new coupled owner."""

    started = time.perf_counter()
    from check_original_macro_defect_kernel import run_fixture

    fixture = run_fixture()
    fixture_summary = _check_fixture(fixture)
    owner = OriginalMacroCoupledM15()
    report = owner.build()
    code = dict(getattr(owner, "implementation_hashes", {}))
    code.update({name: sha(ROOT / name) for name in CODE})
    report.update(input_hashes=owner.input_hashes, implementation_hashes=code)
    report["checks"] = {
        "passed": True,
        "independent_new_defect_fixture": fixture,
        "fixture_invariants": fixture_summary,
        "packet_invariants": check_packet(report, owner.parent),
        "no_accepted_ancestor_replay": True,
    }
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    return report


def inspect(path: Path) -> dict[str, Any]:
    """Inspect only serialized bytes; do not construct an owner or fixture."""

    report = load(path)
    require(isinstance(report, Mapping), "Saved coupled report mapping required")
    checks = report.get("checks")
    require(isinstance(checks, Mapping) and checks.get("passed") is True,
            "Saved coupled report failed")
    fixture = checks.get("independent_new_defect_fixture")
    require(isinstance(fixture, Mapping) and fixture.get("passed") is True,
            "Saved defect fixture failed")
    fixture_summary = _check_fixture(fixture)
    saved_summary = checks.get("fixture_invariants")
    if saved_summary is not None:
        require(saved_summary == fixture_summary, "Saved defect fixture invariants changed")
    packet = check_packet(report, parent())
    return {
        "passed": True,
        "pins": pins(report),
        "packet_invariants": packet,
        "serialized_only": True,
        "owner_or_fixture_or_ancestor_run": False,
    }


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
    else:
        print(json.dumps(inspect(path), indent=2, ensure_ascii=False))
