"""Produce or inspect the actual restoration/post-restore background stage.

The seven source packets are built from the pinned restoration/post-restore
archives and the exact accepted modified-seams report.  The runner never
replays an accepted owner, builder, checker, producer, norm, recurrence, or
quadrature calculation.  Inspection reads only the serialized report and
current source bytes; it does not import any stage module.
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
from collections.abc import Mapping


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "experiments/root_st073"
PARENT_NAME = "lei_ren_part1_paper_compliant_current_original_long_reference_modified_seams.json"
PARENT_PATH = "experiments/root_st073/" + PARENT_NAME
PARENT_SHA256 = "204d8b45aff263a3df1eae319dd7534e304162794ac375a1c0ec2cfe9beafc35"
NAME = "lei_ren_part1_paper_compliant_current_original_restoration_background_I1.json"
RUNNER = "scripts/run_original_restoration_background_I1.py"

PACKETS = (
    ("actual_restoration_pressure_source", "original_restoration_pressure_source"),
    ("actual_restoration_background_I1_numeric", "original_restoration_background_I1_numeric"),
    ("actual_restoration_native_endpoints", "original_restoration_native_endpoints"),
    ("actual_restoration_native_seams", "original_restoration_native_seams"),
    ("actual_restoration_physical_y_dag", "original_restoration_physical_y_dag"),
    ("actual_restoration_phase_mixed", "original_restoration_phase_mixed"),
    ("actual_restoration_correlated_pressure", "original_restoration_correlated_pressure"),
)
CODE = tuple("src/openai_ns_reconstruction/" + module + ".py" for _, module in PACKETS)

# The seven source review hashes are fixed by the current handoff.  The
# independent Bacon/Kepler mathematical review is complete for this runner;
# the report records that review as complete while terminal closure remains
# explicitly disabled in the stage contract.
REVIEWED = {
    CODE[0]: "e52159cd36f34936a05711eab513a8bce0ee78889a613eaf5be495d3a446491f",
    CODE[1]: "aa03b74005d2fff69518baf76298412d4c2ccb668a4d6de9d450d6c3caaff0ee",
    CODE[2]: "d2e5bc67f06fb6f2bb2a14732898817da8773090b2a789fb91f4ec92875407a8",
    CODE[3]: "16a44ff233e64a1ff09fff9d5a21bda4e5ee9708bb67d577909a8e5960a4a91e",
    CODE[4]: "c08973a321690ff2029ef8bede7df5d7c915febbf04ca84d092d5ea84ec390a2",
    CODE[5]: "ee941e87297895371200dffa6486227205cd8ed27adc8dfdf0973578f6b441fb",
    CODE[6]: "64f909cd6fc7be559d73b6924e9ac5efd896a4f187127ef843b06caacf75dbcc",
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


def reviewed_hashes():
    return dict(REVIEWED)


def parent():
    path = BASE / PARENT_NAME
    require(sha(path) == PARENT_SHA256, "Accepted modified-seams parent bytes changed")
    result = load(path)
    require(result["checks"]["passed"] is True, "Accepted parent checks failed")
    return result


def check_pins(report):
    count = 0
    for folder, group in ((BASE, report["input_hashes"]),
                          (ROOT, report["implementation_hashes"])):
        for name, expected in group.items():
            require(isinstance(expected, str) and sha(folder / name) == expected,
                    "Pinned source bytes changed: " + name)
            count += 1
    return count


def _relative(value, label):
    require(isinstance(value, str) and value, "Missing relative source pin: " + label)
    value = value.replace("\\", "/")
    path = Path(value)
    require(not path.is_absolute() and ".." not in path.parts,
            "Unsafe source pin path: " + value)
    return path.as_posix()


def merge_pins(accepted, packets):
    """Merge packet source pins without dropping inherited aliases."""

    inputs = dict(accepted["input_hashes"])
    implementations = dict(accepted["implementation_hashes"])
    inputs[PARENT_NAME] = PARENT_SHA256
    for packet in packets:
        require(packet["passed"] is True
                and packet["parent_report_sha256"] == PARENT_SHA256
                and packet["source_family"] == accepted["source_family"]
                and packet.get("source_executed") is False,
                "Wrong or unchecked restoration source packet")
        pins = packet.get("source_pins")
        require(isinstance(pins, Mapping), "Restoration packet source pins are missing")
        for key, item in pins.items():
            require(isinstance(item, Mapping), "Malformed restoration source pin: " + str(key))
            path_text = _relative(item.get("path", key), str(key))
            key_text = _relative(key, str(key))
            require(path_text == key_text, "Source pin key/path alias changed: " + str(key))
            expected = item.get("expected_sha256")
            require(isinstance(expected, str) and item.get("passed") is True
                    and item.get("source_executed") is not True
                    and ("actual_sha256" not in item or item.get("actual_sha256") == expected),
                    "Malformed or unchecked source pin: " + path_text)
            path = Path(path_text)
            if len(path.parts) == 1 or path_text.startswith("experiments/root_st073/"):
                actual_path = ROOT / path if path_text.startswith("experiments/") else BASE / path
                name, target = path.name, inputs
            else:
                actual_path = ROOT / path
                name, target = path_text, implementations
            require(sha(actual_path) == expected and target.get(name, expected) == expected,
                    "Source pin conflict: " + path_text)
            target[name] = expected
    for path in (*CODE, RUNNER):
        implementations[path] = sha(ROOT / path)
    return inputs, implementations


def _rows(rows, label, count=2):
    require(isinstance(rows, list) and len(rows) == count,
            "Incomplete ordinary-Z rows: " + label)
    require([row.get("ordinary_Z_order") for row in rows] == list(range(count)),
            "Ordinary-Z order changed: " + label)


def _directed_records(value):
    """Check serialized directed interval records without importing mpmath."""

    if isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            low = value.get("lower_exact_mpf_tuple")
            high = value.get("upper_exact_mpf_tuple")
            require(isinstance(low, list) and isinstance(high, list)
                    and len(low) == len(high) == 4
                    and all(type(item) is int for item in low + high)
                    and value.get("directed") is True,
                    "Malformed directed interval record")
            return 1
        return sum(_directed_records(item) for item in value.values())
    if isinstance(value, (list, tuple)):
        return sum(_directed_records(item) for item in value)
    return 0


def _packet(report, key, accepted):
    packet = report.get(key)
    require(isinstance(packet, Mapping) and packet.get("passed") is True,
            "Missing or failed restoration packet: " + key)
    expected_parent_reports = {PARENT_PATH}
    # The pure physical-y DAG intentionally serializes its basename constant;
    # all other packets carry the canonical experiments-relative path.
    if key == "actual_restoration_physical_y_dag":
        expected_parent_reports.add(PARENT_NAME)
    require(packet.get("parent_report") in expected_parent_reports
            and packet.get("parent_report_sha256") == PARENT_SHA256
            and packet.get("source_family") == accepted["source_family"]
            and packet.get("source_executed") is False,
            "Restoration packet provenance changed: " + key)
    pins = packet.get("source_pins")
    require(isinstance(pins, Mapping), "Restoration packet source pins are missing: " + key)
    for pin_key, item in pins.items():
        require(isinstance(item, Mapping), "Malformed packet source pin: " + str(pin_key))
        path_text = _relative(item.get("path", pin_key), str(pin_key))
        key_text = _relative(pin_key, str(pin_key))
        expected = item.get("expected_sha256")
        require(path_text == key_text and isinstance(expected, str)
                and item.get("passed") is True
                and item.get("source_executed") is not True
                and ("actual_sha256" not in item or item.get("actual_sha256") == expected),
                "Malformed packet source pin: " + str(pin_key))
        path = Path(path_text)
        actual_path = BASE / path if len(path.parts) == 1 else ROOT / path
        require(sha(actual_path) == expected,
                "Packet source pin bytes changed: " + path_text)
    return packet


def _validate_source(packet, accepted):
    patch = packet["graph_patch"]
    require(patch["base_node_count"] == 185704
            and patch["new_node_start"] == 185704
            and patch["base_prefix_unchanged"] is True
            and patch["full_node_count"] == 186086
            and patch["appended_node_count"] == len(patch["appended_nodes"])
            and patch["full_node_count"] == patch["base_node_count"] + len(patch["appended_nodes"]),
            "Restoration raw source graph prefix or count changed")
    require(isinstance(patch.get("base_canonical_nodes_sha256"), str)
            and isinstance(patch.get("full_canonical_nodes_sha256"), str),
            "Restoration source graph digests are missing")
    require(packet["joined_source_entry"] == [185700, 185703],
            "Parent joined source handles changed")
    stable = packet["stable_chart_handles"]
    require(stable["parent_joined"] == [185700, 185703]
            and stable["parent_joined_contracted"] == [185700, 185703]
            and stable["new_raw_joined"] == stable["new_raw_joined_uncontracted"]
            and len(stable["new_raw_joined"]) == 2,
            "Restoration raw joined handles changed")
    for chart in ("restoration", "postrestore"):
        handles = stable[chart]
        require(len(handles["lower_RQy_wrappers"]) == 2
                and len(handles["upper_RQy_traces"]) == 2
                and len(handles["local_half_S1"]) == 2
                and len(handles["local_S1"]) == 2,
                "Restoration chart source ports are incomplete: " + chart)
    endpoint = patch["endpoint_handles"]
    require(endpoint["parent_joined"] == [185700, 185703]
            and endpoint["parent_external_reference_RQy"] == [171869, 171870]
            and endpoint["external_Rm_RQy"] == [183696, 183697]
            and endpoint["joined_contracted"] is None,
            "Restoration endpoint handle aliases changed")
    require(packet["preservation"]["accepted_185704_prefix_unchanged"] is True
            and packet["preservation"]["old_212_terms_and_errors_unchanged"] is True
            and packet["preservation"]["old_six_signed_terms_retained"] is True
            and packet["preservation"]["all_standalone_caps_retained"] is True
            and packet["preservation"]["all_error_caps_retained"] is True
            and packet["preservation"]["target_handles_retained"] is True
            and packet["preservation"]["new_seam_atoms_uncontracted"] is True,
            "Accepted restoration data was not retained")
    contract = packet["source_contract"]
    require(contract["physical_Pstar_applied_once"] is True
            and contract["physical_R_squared_applied_once"] is True
            and contract["physical_Jacobian_applied_once"] is True
            and contract["nonreset_memories"] is True
            and contract["all_six_signed_I1_sources_defined"] is True
            and contract["D_B_S1_sources_defined"] is True
            and contract["seam_differences_contracted"] is False
            and contract["external_Rm_RQy_retained"] is True,
            "Restoration source scale, memory, or seam scope changed")
    limits = packet["limits"]
    require(limits["all_caps_retained"] is True
            and limits["all_old_data_retained"] is True
            and limits["new_seams_contracted"] is False
            and limits["external_Rm_RQy_retained"] is True
            and limits["mathematical_completion"] is False,
            "Restoration source limits were promoted")
    for chart in ("restoration", "postrestore"):
        data = packet["charts"][chart]
        require(data["physical_R_squared_applied_once"] is True
                and data["physical_Jacobian_applied_once"] is True
                and data["Pstar_applied_once"] is True
                and data["source_values_selected"] is False
                and data["interval_or_cap_values_used"] is False,
                "Restoration chart source formula changed: " + chart)
        six = data["six_I1_source"]["integral_handles"]
        require(set(six) == {"q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z"}
                and all(len(rows) == 2 for rows in six.values()),
                "Restoration chart does not expose all six signed I1 rows: " + chart)
        require(len(data["D_source"]["rows"]) == 2
                and len(data["B_source"]["rows"]) == 2
                and len(data["S1_source"]["rows"]) == 2,
                "Restoration chart D/B/S1 source rows are incomplete: " + chart)
    return dict(passed=True, raw_source_graph=True, source_ports=True,
                six_signed_I1=True, D_B_S1=True, new_seams_uncontracted=True)


def _validate_numeric(packet, accepted):
    fixture = packet["fixture_bound_check"]
    require(fixture["passed"] is True
            and fixture["alpha_orientation"] == "alpha=1-sigma decreases 1 to 0"
            and fixture["gamma_initial"] == "gamma(0)=1"
            and fixture["no_source_executed"] is True,
            "Numeric fixture-bound receipt is missing or changed")
    require(packet["accepted_ancestors_replayed"] is False
            and packet["ordinary_Z_orders"] == [0, 1, 2, 3, 4]
            and packet["I1_ordinary_Z_orders"] == [0, 1],
            "Restoration numeric scope changed")
    scale = packet["physical_scale"]
    require(scale["Pstar_handle"] == 33572 and scale["applied_once"] is True
            and scale["Wrestore"] == "Rrestore^2" and scale["Wpost"] == "Rm^2"
            and scale["radius_materialized"] is False,
            "Restoration physical scale was changed or selected")
    contract = packet["source_contract"]
    require(contract["retained_boundary_Qy"] is True
            and contract["all_six_I1"] is True
            and contract["actual_signed_C1_interval_enclosures"] is True
            and contract["no_modified_mean_substitution"] is True
            and contract["no_radius_materialization"] is True,
            "Restoration numeric source contract changed")
    runs = packet["actual_runs"]
    require(runs["interval_precision_bits"] == 1600
            and runs["saved_full_unit_restoration_kernel_receipt"] is True
            and runs["saved_kernel_cell_count"] == 128
            and runs["missing_endpoint_memory_kernel_cells"] == 512
            and runs["directed_monotone_enclosure"] is True
            and runs["positive_scalar_kernels"] is True
            and runs["point_or_RK_proxy_used_as_value"] is False
            and runs["interval_midpoint_used_as_source"] is False
            and runs["accepted_ancestor_computation_replayed"] is False,
            "Restoration numeric provenance or enclosure mode changed")
    numeric_limits = packet["limits"]
    require(numeric_limits["finite_N_signed_centers_resolved"] is False
            and numeric_limits["global_five_moment_numeric_join"] is False
            and numeric_limits["common_global_N"] is False,
            "Restoration numeric finite-N scope was promoted")
    total_records = 0
    for chart in ("restoration", "post_restore"):
        chart_packet = packet["charts"][chart]
        arithmetic = chart_packet["arithmetic"]
        require(set(arithmetic["six_integrals"]) ==
                {"q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z"},
                "Restoration numeric six-integral set changed: " + chart)
        for name, rows in arithmetic["six_integrals"].items():
            _rows(rows, chart + " " + name)
            total_records += _directed_records(rows)
        for name in ("D", "B_over_W", "S1_over_W"):
            _rows(arithmetic[name], chart + " " + name)
            total_records += _directed_records(arithmetic[name])
        for name, rows in arithmetic["boundary_Qy"].items():
            _rows(rows, chart + " boundary_Qy " + name, count=2)
            total_records += _directed_records(rows)
        scalar = chart_packet["scalar_kernel_receipt"]
        if chart == "restoration":
            require(scalar["used_as_interval_value"] is True
                    and scalar["point_proxy_as_value"] is False,
                    "Saved scalar kernel receipt was promoted to a point value: " + chart)
            scalar_contract = scalar["contract"]
            require(scalar_contract["original_alpha_source"] is True
                    and scalar_contract["gamma_ODE_positive_recurrence"] is True
                    and scalar_contract["whole_cell_monotone_gamma"] is True
                    and scalar_contract["interval_midpoint_not_source"] is True,
                    "Saved universal-alpha scalar kernel contract changed: " + chart)
        else:
            require(scalar["exact_pure_exponential_primitives"] is True
                    and scalar["source_executed"] is False
                    and scalar["normalized_by"] == "exp(-2) because Wpost=Rm^2"
                    and scalar["point_proxy_as_value"] is False,
                    "Post-restore scalar kernel source changed: " + chart)
        require(chart_packet["source_executed"] is False
                and chart_packet["normalization"] in ("Wrestore=Rrestore^2", "Wpost=Rm^2"),
                "Restoration chart packet provenance changed: " + chart)
    restoration_source = packet["charts"]["restoration"]["native_source"]
    post_source = packet["charts"]["post_restore"]["native_source"]
    require(restoration_source["formula"] ==
            "U=4Z+E*(1-sigma(t)); M=4Z+E*gamma(t)+K*exp(-t)"
            and restoration_source["M_y_formula"] == "M_y=U-M"
            and post_source["formula"] == "U=4Z; M=4Z+(M_restore-4Z)*exp(-q)"
            and post_source["M_y_formula"] == "M_y=U-M",
            "Native U/M restoration decomposition changed")
    gamma = packet["charts"]["restoration"]["arithmetic"]["source"]["gamma_proof"]
    require(gamma["cells"] == 512
            and gamma["alpha_nonincreasing"] is True
            and gamma["positive_cell_recurrence"] is True
            and gamma["point_or_RK_proxy_used_as_value"] is False
            and gamma["source_executed"] is False
            and len(gamma["cell_receipts"]) == 512
            and all(cell["positive_kernel"] is True and cell["monotone_alpha"] is True
                    for cell in gamma["cell_receipts"]),
            "Monotone gamma endpoint enclosure changed")
    require(total_records > 0, "Restoration numeric packet has no directed interval rows")
    return dict(passed=True, both_charts=True, six_signed_I1=True,
                D_B_S1_signed_C1=True, directed_records=total_records,
                retained_boundary_Qy=True)


def _validate_endpoint(packet, accepted):
    require(packet["native_phi_U_physical_C2"] is True
            and packet["radius_endpoint_identities"] is True
            and packet["ordinary_Z_orders"] == [0, 1, 2, 3],
            "Native endpoint C2 scope changed")
    symbolic = packet["symbolic_endpoint_identities"]
    require(symbolic["passed"] is True
            and symbolic["function_identity_not_sampling"] is True,
            "Native endpoint identities are not function identities")
    expected = {
        "reference_to_restoration": ("reference", "restoration"),
        "restoration_to_postrestore": ("restoration", "post_restore"),
    }
    seams = symbolic["seams"]
    require(isinstance(seams, list) and len(seams) == 2,
            "Native endpoint seam count changed")
    seen = set()
    for seam in seams:
        name = seam["name"]
        require(name in expected and seam["left_chart"] == expected[name][0]
                and seam["right_chart"] == expected[name][1]
                and seam["passed"] is True,
                "Native endpoint seam changed: " + name)
        require(all(value == "0" for value in seam["residuals"].values())
                and all(isinstance(values, list) and len(values) == 4
                        and all(value == "0" for value in values)
                        for values in seam["ordinary_Z0_Z3_residuals"].values()),
                "Native endpoint residual is nonzero: " + name)
        require(name not in seen, "Native endpoint seam name is duplicated: " + name)
        seen.add(name)
    require(seen == set(expected), "Native endpoint seam names changed")
    require(packet["limits"]["current_N_numeric"] is False,
            "Native endpoint finite-N scope was promoted")
    return dict(passed=True, endpoint_C2=True, radius_identities=True,
                ordinary_Z0_to_Z3=True)


def _validate_native(packet, accepted):
    require(packet["accepted_ancestors_replayed"] is False
            and all(packet[name] is False for name in (
                "accepted_builders_executed", "accepted_checkers_executed",
                "accepted_producers_executed", "owners_executed", "norms_executed",
                "recurrences_executed", "quadrature_executed"))
            and packet["ordinary_Z_orders"] == list(range(7))
            and packet["physical_y_orders"] == [0, 1, 2],
            "Native history execution or scope changed")
    dependency = packet["endpoint_dependency"]
    require(dependency["consumed"] is True
            and dependency["validation"]["passed"] is True,
            "Native seams did not consume the actual endpoint packet")
    require(packet["actual_native_histories"] is True
            and packet["actual_C5_correction_histories"] is True
            and packet["native_history_source"]["passed"] is True
            and packet["c5_history_source"]["passed"] is True
            and packet["q4_background_seed_source"]["passed"] is True,
            "Native nonreset histories or Q4 seed binding changed")
    require(packet["limits"]["current_N_numeric"] is False,
            "Native history finite-N scope was promoted")
    return dict(passed=True, endpoint_consumed=True, nonreset_histories=True,
                native_C5_sources=True)


def _validate_dag(packet, accepted):
    require(packet["conditional_passed"] is True
            and packet["all_passed"] is False
            and packet["seed_identity_discharged"] is False
            and packet["source_executed"] is False
            and packet["pure_generic_rule_engine_new24_only"] is True,
            "Restoration physical-y DAG scope changed")
    counts = packet["counts"]
    require(counts["slow_derivative_roots_total"] == 24
            and counts["slow_derivative_roots_passed"] == 24,
            "Restoration physical-y DAG root count changed")
    audit = packet["generic_rule_engine_audit"]
    require(audit["applied_root_count"] == 24
            and audit["accepted_old_36_root_arithmetic_replayed"] is False
            and audit["accepted_old_builders_or_checkers_replayed"] is False
            and audit["old_generic_receipt_consumed_as_saved_lemma"] is True,
            "Restoration DAG replay audit changed")
    require(packet["source_family_match"] is True
            and packet["source_pins_passed"] is True
            and all(value is True for value in packet["gates"].values()),
            "Restoration physical-y DAG gate failed")
    require(packet["limits"]["no_numeric_values"] is True
            and packet["limits"]["finite_N_modified_fields_used"] is False,
            "Restoration physical-y DAG numeric scope was promoted")
    return dict(passed=True, conditional_new24=True, old36_replayed=False)


def _validate_phase(packet, accepted):
    receipt = packet["phase_coefficients"]
    require(packet["charts"] == ["restoration", "postrestore"]
            and receipt["passed"] is True
            and receipt["total_rows"] == receipt["matched_rows"] == 42,
            "Restoration phase coefficient receipt changed")
    lambda_rows = 0
    for chart in ("restoration", "postrestore"):
        rows = receipt["charts"][chart]
        require(rows["matched_rows"] == 21
                and rows["raw_lambda_rows_matched"] == 15
                and rows["raw_lambda_convolution_bound"] is True
                and all(row["passed"] is True for row in rows["rows"]),
                "Restoration phase rows changed: " + chart)
        lambda_rows += rows["raw_lambda_rows_matched"]
    require(lambda_rows == 30, "Restoration lambda row count changed")
    require(packet["limits"]["global_five_moment_closure"] is False
            and packet["limits"]["full_goal_complete"] is False,
            "Restoration phase scope was promoted")
    return dict(passed=True, mixed_phase_rows=42, lambda_rows=30)


def _validate_derivative_discharge(discharge):
    """Validate the actual native/endpoint premises consumed by correlation."""

    identity = discharge["actual_seed_function_identity"]
    history = identity["actual_new_history_binding"]
    require(history["passed"] is True
            and history["source_FTC"] is True
            and history["residual_count"] == 0,
            "Actual new native history source-FTC receipt changed")
    seeds = identity["actual_new_9x7_seed_bindings"]
    require(set(seeds) == {"restoration", "postrestore"},
            "Actual new 9x7 seed chart names changed")
    for chart, seed in seeds.items():
        require(seed["passed"] is True
                and seed["native_fields"] == 9
                and seed["ordinary_rows"] == 7
                and seed["residual_count"] == 0,
                "Actual 9x7 seed receipt changed: " + chart)
        normalized = seed["normalized_field_rows"]
        require(len(normalized) == 42
                and all(row["residual"] == "0" for row in normalized),
                "Actual normalized field rows changed: " + chart)
        weight = seed["W_source"]
        squares = weight["W_squared_equals_2FF"]
        require(len(squares) == 7
                and all(row["residual"] == "0" for row in squares),
                "Actual W squared source rows changed: " + chart)
        derivative = weight["derivative_source"]
        relation = derivative["radius_identity"]
        require(relation["R_q_equals_JR"] is True
                and all(derivative["phase_independence"].values())
                and all(value == "0" for value in derivative["residuals"].values()),
                "Actual radial derivative source changed: " + chart)
    controls = discharge["smooth_branch_and_positive_IFT"]
    endpoint_control = controls["endpoint_control"]
    radial = endpoint_control["endpoint_radial_source"]
    require(endpoint_control["passed"] is True
            and radial["passed"] is True
            and set(radial["seam_sides"]) == {
                "reference_to_restoration", "restoration_to_postrestore",
            }, "Actual endpoint radial control receipt changed")
    for seam, sides in radial["seam_sides"].items():
        require(sides["left"]["chart_phase"] == "1"
                and sides["right"]["chart_phase"] == "0",
                "Actual endpoint seam side phase changed: " + seam)
        for side_name in ("left", "right"):
            wrappers = sides[side_name]["source_wrappers"]
            require(isinstance(wrappers, Mapping) and wrappers
                    and all(len(rows) == 7 for rows in wrappers.values()),
                    "Actual endpoint wrapper rows changed: " + seam + "." + side_name)
    require(controls["local_IFT"]["actual_L"]["residual"] == "0",
            "Actual Qy denominator residual changed")


def _validate_correlation(packet, accepted):
    require(packet["accepted_ancestors_replayed"] is False
            and packet["full_modified_seam_proved"] is True
            and packet["composition_conditional"] is False
            and packet["unresolved_gates"] == []
            and all(value is True for value in packet["gates"].values())
            and packet["axial_scope"] == {
                "Z_interval": ["37/100", "93/250"],
                "ordinary_Z_orders": [0, 1],
            },
            "Restoration modified seam proof packet changed")
    require(set(packet["seams"]) == {
        "reference_to_restoration", "restoration_to_postrestore",
    }, "Restoration seam names changed")
    discharge = packet["derivative_premise_discharge"]
    require(discharge["issued_DAG_remains_conditional"] is True
            and discharge["issued_DAG_seed_identity_discharged"] is False
            and discharge["scoped_to_two_seams_and_inherited_Z_window"] is True
            and discharge["global_chart_derivative_numeric_claim"] is False,
            "Restoration derivative premise discharge scope changed")
    _validate_derivative_discharge(discharge)
    patch = packet["graph_patch"]
    require(patch["base_node_count"] == 186086
            and patch["base_prefix_unchanged"] is True
            and patch["full_node_count"] == 186100
            and patch["full_node_count"] == patch["base_node_count"] + len(patch["appended_nodes"])
            and isinstance(patch["base_canonical_nodes_sha256"], str)
            and isinstance(patch["full_canonical_nodes_sha256"], str),
            "Restoration correlated graph prefix or count changed")
    require(packet["joined_contracted"] == [186096, 186099]
            and packet["outer_Rm_RQy_retained"] == [183696, 183697],
            "Restoration joined or external Rm handles changed")
    linear = packet["exact_linearization"]
    require(linear["passed"] is True
            and linear["only_sum_negative_expanded"] is True
            and linear["conservative_integer_dependency_check"] is True
            and linear["targets"] == [171869, 171870, 178093, 178094],
            "Restoration signed linearization scope changed")
    for row in linear["rows"]:
        require(all(value == 0 for value in row["target_coefficients"].values())
                and not any(row["surviving_target_dependencies"].values()),
                "A removed restoration endpoint target survives")
    limits = packet["limits"]
    require(limits["current_N_restoration_numeric"] is False,
            "Restoration finite-N scope changed")
    require(limits["complete_pressure_numeric"] is False
            and limits["standalone_caps_retained"] is True
            and limits["all_old_data_retained"] is True
            and limits["global_five_moment_closure"] is False
            and limits["full_goal_complete"] is False,
            "Restoration correlation limits were promoted")
    return dict(passed=True, full_two_new_modified_seams=True,
                raw_graph_186086=True, contracted_graph_186100=True,
                external_Rm_retained=True)


def validate_report(report, accepted):
    mutable = {"stage", "parent_report", "parent_report_sha256", "input_hashes",
               "implementation_hashes", "contract", "checks", "execution_seconds",
               "execution_audit"}
    for key, value in accepted.items():
        if key not in mutable:
            require(report.get(key) == value, "Accepted mathematical payload changed: " + key)
    require(report["parent_report"] == PARENT_NAME
            and report["parent_report_sha256"] == PARENT_SHA256,
            "Wrong accepted predecessor")
    source = _packet(report, PACKETS[0][0], accepted)
    numeric = _packet(report, PACKETS[1][0], accepted)
    endpoint = _packet(report, PACKETS[2][0], accepted)
    native = _packet(report, PACKETS[3][0], accepted)
    dag = _packet(report, PACKETS[4][0], accepted)
    phase = _packet(report, PACKETS[5][0], accepted)
    correlation = _packet(report, PACKETS[6][0], accepted)
    results = {
        PACKETS[0][0]: _validate_source(source, accepted),
        PACKETS[1][0]: _validate_numeric(numeric, accepted),
        PACKETS[2][0]: _validate_endpoint(endpoint, accepted),
        PACKETS[3][0]: _validate_native(native, accepted),
        PACKETS[4][0]: _validate_dag(dag, accepted),
        PACKETS[5][0]: _validate_phase(phase, accepted),
        PACKETS[6][0]: _validate_correlation(correlation, accepted),
    }
    audit = report["execution_audit"]
    require(audit["inherited_parent_execution_audit"] == accepted["execution_audit"]
            and audit["accepted_ancestor_math_replayed"] is False
            and audit["accepted_builders_or_checkers_replayed"] is False
            and audit["accepted_producer_or_owner_replayed"] is False
            and audit["generic_derivative_engine_applied_new24_only"] is True
            and audit["restoration_numeric_fixture_bound_check_once"] is True,
            "Execution audit was reset")
    contract = report["contract"]
    require(contract["actual_restoration_background_I1"] is True
            and contract["actual_restoration_pressure_source"] is True
            and contract["actual_restoration_native_endpoints"] is True
            and contract["actual_restoration_native_seams"] is True
            and contract["actual_restoration_physical_y_dag"] is True
            and contract["actual_restoration_phase_mixed"] is True
            and contract["actual_restoration_correlated_pressure"] is True
            and contract["actual_endpoint_C2"] is True
            and contract["new_24_physical_y_dag_roots"] is True
            and contract["new_42_phase_rows"] is True
            and contract["new_30_lambda_rows"] is True
            and contract["source_graph_185704_to_186086"] is True
            and contract["full_two_new_modified_RQy_seams"] is True
            and contract["current_N_restoration_numeric"] is False
            and contract["complete_pressure_numeric"] is False
            and contract["global_five_moment_closure"] is False
            and contract["full_goal_complete"] is False,
            "Restoration stage contract changed")
    return dict(passed=True, packet_invariants=results,
                accepted_math_payload_preserved=True,
                all_old_data_retained=True, all_six_signed_rows=True,
                ordinary_C1_both_charts=True)


def _reviewed_check(report):
    current = reviewed_hashes()
    saved = report["checks"]["reviewed_implementation_sha256"]
    require(saved == current, "Saved implementation review hashes changed")
    for path, expected in REVIEWED.items():
        require(sha(ROOT / path) == expected, "Reviewed source bytes changed: " + path)
    require(report["implementation_hashes"][CODE[-1]] == sha(ROOT / CODE[-1]),
            "Correlated pressure self-pin changed")
    return current


def run():
    started = time.perf_counter()
    _reviewed = reviewed_hashes()
    for path, expected in REVIEWED.items():
        require(sha(ROOT / path) == expected, "Reviewed source bytes changed: " + path)
    accepted = parent()
    check_pins(accepted)
    sys.path.insert(0, str(ROOT / "src"))
    modules = [importlib.import_module("openai_ns_reconstruction." + module)
               for _, module in PACKETS]
    source = modules[0].build(ROOT, accepted)
    numeric = modules[1].build(ROOT, accepted)
    fixture = modules[1].fixture_bound_check()
    require(fixture["passed"] is True, "Restoration numeric fixture-bound check failed")
    numeric = copy.deepcopy(numeric)
    numeric["fixture_bound_check"] = fixture
    endpoint = modules[2].build(ROOT, accepted)
    native = modules[3].build(ROOT, accepted, endpoint)
    dag = modules[4].build(ROOT, accepted)
    phase = modules[5].build(ROOT, accepted)
    correlation = modules[6].build(ROOT, accepted, source, endpoint, native, dag, phase)
    packets = (source, numeric, endpoint, native, dag, phase, correlation)
    require(accepted == parent(), "New stage changed its accepted input")
    inputs, implementations = merge_pins(accepted, packets)
    report = copy.deepcopy(accepted)
    report.update(
        stage="actual_restoration_postrestore_background_I1_and_pressure_source",
        parent_report=PARENT_NAME, parent_report_sha256=PARENT_SHA256,
        input_hashes=inputs, implementation_hashes=implementations,
    )
    for (key, _), packet in zip(PACKETS, packets):
        report[key] = packet
    report["execution_audit"] = dict(
        accepted_ancestor_math_replayed=False,
        accepted_builders_or_checkers_replayed=False,
        accepted_producer_or_owner_replayed=False,
        accepted_norms_recurrences_quadrature_replayed=False,
        current_turn_accepted_ancestors_replayed=False,
        generic_derivative_engine_applied_new24_only=True,
        new_source_and_background_I1_only=True,
        restoration_numeric_fixture_bound_check_once=True,
        inherited_parent_execution_audit=copy.deepcopy(accepted["execution_audit"]),
    )
    contract = copy.deepcopy(accepted["contract"])
    contract.update(
        actual_restoration_pressure_source=True,
        actual_restoration_background_I1=True,
        actual_restoration_native_endpoints=True,
        actual_restoration_native_seams=True,
        actual_restoration_physical_y_dag=True,
        actual_restoration_phase_mixed=True,
        actual_restoration_correlated_pressure=True,
        actual_native_restoration_background_I1=True,
        actual_endpoint_C2=True,
        new_24_physical_y_dag_roots=True,
        new_42_phase_rows=True,
        new_30_lambda_rows=True,
        source_graph_185704_to_186086=True,
        full_two_new_modified_RQy_seams=True,
        current_N_restoration_numeric=False,
        complete_pressure_numeric=False,
        global_five_moment_closure=False,
        full_goal_complete=False,
    )
    report["contract"] = contract
    report["checks"] = dict(
        passed=True,
        packet_invariants=validate_report(report, accepted),
        independent_new_math_review_completed=True,
        independent_runtime_worker_model_metadata_available=False,
        correlated_pressure_review_pending=False,
        reviewed_implementation_sha256=_reviewed,
    )
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    json.dumps(report, allow_nan=False)
    return report


def inspect(path):
    """Validate a saved report from JSON and source hashes only."""

    report, accepted = load(path), parent()
    require(report["checks"]["passed"] is True
            and report["checks"]["independent_new_math_review_completed"] is True
            and report["checks"]["correlated_pressure_review_pending"] is False,
            "Saved restoration review status is incomplete")
    _reviewed_check(report)
    for path_text in (*CODE, RUNNER):
        require(report["implementation_hashes"][path_text] == sha(ROOT / path_text),
                "Saved implementation bytes changed: " + path_text)
    check_pins(report)
    return dict(passed=True, pins=len(report["input_hashes"]) + len(report["implementation_hashes"]),
                packet_invariants=validate_report(report, accepted),
                serialized_only=True, accepted_builders_or_checkers_replayed=False,
                independent_new_math_review_completed=True,
                correlated_pressure_review_pending=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report", action="store_true")
    group.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    path = BASE / NAME
    if args.write_report:
        require(not path.exists(), "Saved restoration report already exists; use --inspect-report")
        report = run()
        tmp = path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                       encoding="utf-8", newline="\n")
        tmp.replace(path)
        print(json.dumps(dict(passed=True, execution_seconds=report["execution_seconds"],
                              bytes=path.stat().st_size, sha256=sha(path))))
    else:
        print(json.dumps(inspect(path), indent=2))
