"""Write or inspect the saved postrestore/Rm inlet m15 packet.

``--write-report`` constructs the new owner once and runs the independent
fixture once.  ``--inspect-report`` reads only serialized reports and pinned
sources; it never constructs the owner or reruns the fixture.  The two modes
are mutually exclusive so a routine inspection cannot replay production.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
SCRIPTS = ROOT / "scripts"
DIRECTORY = ROOT / "experiments" / "root_st073"
for path in (str(ROOT), str(SRC), str(SCRIPTS)):
    if path not in sys.path:
        sys.path.insert(0, path)

from openai_ns_reconstruction.original_core_numeric import PREFIX  # noqa: E402
from openai_ns_reconstruction.original_Rm_inlet_m15 import (  # noqa: E402
    PARENT,
    PARENT_SHA,
    OriginalRmInletM15,
)
from openai_ns_reconstruction.original_Rm_inlet_numeric_m15 import (  # noqa: E402
    INCOMING_COUNT,
    NEW_EXTERNAL_TERM,
    NEW_REDUCED_COUNT,
    OLD_EXTERNAL_TERM,
    OLD_RM_TERM,
    RM_ARCHIVE,
    RM_ARCHIVE_SHA256,
    RM_C6_SOURCE,
    RM_C6_SOURCE_SHA256,
    RM_M15_SOURCE,
    RM_M15_SOURCE_SHA256,
    RM_PRESSURE_REPORT,
    RM_PRESSURE_REPORT_SHA256,
    RM_REPORT,
    RM_REPORT_SHA256,
    RM_SOURCE,
    RM_SOURCE_SHA256,
    RM_VISCOUS_INDEX,
    RAW_COUNT,
    SAVED_Z,
    UNCHANGED_COUNT,
    PREVIOUS_REDUCED_COUNT,
)


NAME = PREFIX + "compliant_current_original_Rm_inlet_m15_numeric.json"
FIXTURE_NAME = "scripts/check_original_Rm_inlet_m15.py"
RUNNER_NAME = "scripts/run_original_Rm_inlet_m15.py"
MEAN_SOURCE = "src/openai_ns_reconstruction/original_Rm_inlet_mean.py"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def digest(value: object) -> str:
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def require(condition: bool, message: object) -> None:
    if not condition:
        raise AssertionError(message)


def read_json(path: Path):
    return json.loads(path.read_bytes())


def check_hash_group(base: Path, group, label: str) -> None:
    require(isinstance(group, dict) and group, f"{label} hash map missing")
    for name, expected in group.items():
        require(isinstance(name, str) and isinstance(expected, str) and len(expected) == 64,
                (label, name, "malformed digest"))
        target = base / name
        require(target.is_file(), (label, name, "pinned file missing"))
        require(sha(target) == expected, (label, name, "pinned file changed"))


def _saved_path(name: str) -> Path:
    candidate = DIRECTORY / name
    if candidate.is_file():
        return candidate
    return ROOT / "experiments" / "root_st073" / name


def load_reports():
    parent_path = _saved_path(PARENT)
    require(sha(parent_path) == PARENT_SHA, "accepted long/reference parent changed")
    parent = read_json(parent_path)
    require(parent.get("checks", {}).get("passed") is True, "accepted long/reference parent is not passed")
    rm = read_json(_saved_path(RM_REPORT))
    pressure = read_json(_saved_path(RM_PRESSURE_REPORT))
    require(rm.get("checks", {}).get("passed") is True, "accepted Rm m15 packet is not passed")
    require(pressure.get("checks", {}).get("passed") is True, "accepted Rm pressure packet is not passed")
    require(rm.get("source_family") == parent.get("source_family"), "Rm source family changed")
    require(pressure.get("source_family") == parent.get("source_family"), "Rm pressure source family changed")
    require(sha(_saved_path(RM_REPORT)) == RM_REPORT_SHA256, "saved Rm m15 report changed")
    require(sha(_saved_path(RM_PRESSURE_REPORT)) == RM_PRESSURE_REPORT_SHA256,
            "saved Rm pressure report changed")
    return parent, rm, pressure


def _ctx_pair():
    ctx = mp.mp.clone()
    ctx.dps = 150
    iv = MPIntervalContext()
    iv.dps = 180
    return ctx, iv


def _record_tuples(row):
    require(isinstance(row, dict), ("record", row))
    values = row.get("source_enclosure_exact_mpf_tuples")
    if values is None and "lower_exact_mpf_tuple" in row:
        values = [row["lower_exact_mpf_tuple"], row["upper_exact_mpf_tuple"]]
    if values is None and isinstance(row.get("log_absolute_upper"), dict):
        return _record_tuples(row["log_absolute_upper"])
    require(isinstance(values, list) and len(values) == 2, ("record tuples", row))
    require(all(isinstance(value, list) and len(value) == 4 for value in values),
            ("record tuple shape", row))
    return values


def _decode_record(row, ctx, iv):
    tuples = _record_tuples(row)
    return iv.mpf([ctx.make_mpf(tuple(value)) for value in tuples])


def _bounds(ctx, value):
    return tuple(ctx.make_mpf(row) for row in value._mpi_)


def _check_record(row, label):
    require(isinstance(row, dict), (label, "missing record"))
    _record_tuples(row)
    require(row.get("positive_error_bound_atom_only") is True, (label, "not a positive bound"))
    require(row.get("source_value_selected") is False, (label, "cap used as source value"))


def _check_positive_union(values, record, label):
    """Recompute the serialized positive union from exact MPF tuples."""
    _check_record(record, label)
    ctx, iv = _ctx_pair()
    decoded = [_decode_record(value, ctx, iv) for value in values if value is not None]
    require(decoded, (label, "empty positive union"))
    top = max(_bounds(ctx, value)[1] for value in decoded)
    expected = iv.mpf(top) + iv.log(len(decoded))
    expected_tuples = [list(value) for value in expected._mpi_]
    require(record["source_enclosure_exact_mpf_tuples"] == expected_tuples,
            (label, "positive union arithmetic changed"))
    return {"components": len(decoded), "exact_mpf_union": True, "cap_as_value": False}


def _find(rows, **criteria):
    for row in rows or ():
        if isinstance(row, dict) and all(row.get(key) == value for key, value in criteria.items()):
            return row
    raise AssertionError(("missing row", criteria))


def _check_fixture(fixture):
    require(isinstance(fixture, dict) and fixture.get("passed") is True,
            "independent Rm inlet fixture failed")
    expected = {
        "synthetic_parameters_only": True,
        "original_source_values_evaluated": False,
        "finite_N": 17,
        "actual_Poisson_inverse_evaluated": False,
        "variable_J": "1/x",
        "physical_y": "log(x)",
        "physical_second_derivative_keeps_J_x": True,
        "nonzero_incoming_mean": True,
        "ordinary_Z1_retained": True,
        "not_an_NS_residual": True,
    }
    for key, value in expected.items():
        require(fixture.get(key) == value, ("fixture flag", key, fixture.get(key)))
    controls = fixture.get("controls")
    expected_controls = {
        "variable_J_derivative_required",
        "no_extra_N_on_fast_term",
        "incoming_memory_not_zero",
        "lower_minus_upper_sign",
        "physical_R_required",
        "leading_subtraction_retained",
        "phase_not_reset",
        "ordinary_mixed_Z_derivative_nonzero",
        "external_inlet_retained",
        "external_outlet_retained",
    }
    require(isinstance(controls, dict) and set(controls) == expected_controls and
            all(value is True for value in controls.values()), ("fixture controls", controls))
    require(fixture.get("control_count") == 10, "fixture control count changed")
    require(mp.mpf(fixture.get("maximum_supporting_difference")) < mp.mpf("1e-40"),
            "fixture tolerance changed")
    differences = fixture.get("chain_C0_Z1_differences")
    require(isinstance(differences, list) and len(differences) == 2 and
            all(mp.mpf(value) < mp.mpf("1e-40") for value in differences),
            "fixture chain differences changed")
    return {
        "passed": True,
        "finite_N": 17,
        "control_count": len(controls),
        "tolerance": "1e-40",
        "maximum_supporting_difference": fixture.get("maximum_supporting_difference"),
    }


def _check_flat(flat):
    require(isinstance(flat, dict) and flat.get("passed") is True, "Rm flat/history certificate failed")
    require(flat.get("source_only") is True and flat.get("numeric_caps_used") is False and
            flat.get("global_opaque_alias_used") is False, "flat source/value contract changed")
    scope = flat.get("scope", {})
    require(scope.get("postrestore_native_endpoint") == 1 and scope.get("Rm_patch_native_x") == 1,
            "flat endpoint scope changed")
    for key in ("compact_raw_beta", "unique_implicit_control", "incoming_history"):
        require(flat.get(key, {}).get("passed") is True, ("flat premise", key))
    # The raw modified-field comparison is deliberately a failed, unproved
    # diagnostic.  It must remain visible as only 3 matching rows out of the
    # 27-row source comparison while the independent flat certificate itself
    # stays passed.  This is a raw modified-field probe, not a leading-only
    # comparison.
    modified_probe = flat.get("raw_modified_field_probe")
    require(isinstance(modified_probe, dict) and modified_probe.get("passed") is False,
            "raw modified-field probe must remain failed/unproved")
    require(modified_probe.get("actual_source_group") == "modified_fields" and
            modified_probe.get("not_a_leading_only_comparison") is True,
            "raw modified-field probe source group changed")
    modified_rows = modified_probe.get("rows")
    require(isinstance(modified_rows, list) and len(modified_rows) == 27,
            "raw modified-field probe row count changed")
    matching_rows = sum(row.get("passed") is True for row in modified_rows
                        if isinstance(row, dict))
    require(matching_rows == 3, ("raw modified-field probe 3/27 mismatch", matching_rows))
    if "expected_rows" in modified_probe:
        require(modified_probe.get("expected_rows") == 27,
                "raw modified-field expected row count changed")
    if "matched_rows" in modified_probe:
        require(modified_probe.get("matched_rows") == 3,
                "raw modified-field matched row count changed")
    require(modified_probe.get("source_only") is True and
            modified_probe.get("numeric_values_used") is False,
            "raw modified-field probe source/value contract changed")
    limits = flat.get("limits", {})
    require(isinstance(limits, dict) and all(value is False for value in limits.values()),
            "flat limits changed")
    return {"passed": True, "scope": "postrestore1/Rm_patch_x1", "limits_all_false": True}


def _check_native_guards(guards):
    require(isinstance(guards, dict) and guards.get("passed") is True and
            guards.get("source_only") is True and guards.get("ancestor_replayed") is False,
            "Rm native source guard failed")
    require(guards.get("numeric_caps_used") is False and guards.get("source_values_selected") is False,
            "native source/value contract changed")
    scope = guards.get("scope", {})
    require(scope.get("chart") == "Rm_patch" and isinstance(scope.get("native_coordinate"), int) and
            isinstance(scope.get("native_radius"), int) and isinstance(scope.get("native_physical_measure"), int),
            "Rm native scope changed")
    premises = guards.get("premises", {})
    expected = {
        "R_equals_Rm_times_x": True,
        "J_equals_one_over_x": True,
        "d_y_equals_x_d_x": True,
        "native_J_is_varying": True,
        "native_J_is_constant": False,
        "global_offset_log_R_over_Ra_minus_phase_origin": True,
        "phase_uses_common_N": True,
        "native_background_seed_is_J_times_physical_row": True,
        "slow_derivative_divides_by_moving_J_once": True,
        "physical_U_y_has_one_Pstar": True,
        "physical_R_Q_y_has_one_radius": True,
        "gamma_x_derivatives_seeded": True,
    }
    for key, value in expected.items():
        require(premises.get(key) == value, ("native premise", key, premises.get(key)))
    limits = guards.get("limits", {})
    require(isinstance(limits, dict) and all(value is False for value in limits.values()),
            "native guard limits changed")
    return {"passed": True, "scope": "Rm_patch x=1", "moving_J": True, "limits_all_false": True}


def _check_theorem(theorem):
    require(isinstance(theorem, dict) and theorem.get("passed") is True,
            "Rm source derivative theorem failed")
    for key, value in {
        "source_only": True,
        "numeric_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
    }.items():
        if key in theorem:
            require(theorem.get(key) == value, ("source theorem contract", key, theorem.get(key)))
    rows = theorem.get("rows")
    require(isinstance(rows, list) and len(rows) == 2 and
            {row.get("ordinary_Z_order") for row in rows} == {0, 1} and
            all(row.get("passed") is True for row in rows), "Rm theorem Uy rows changed")
    u_rows = theorem.get("U_rows")
    require(isinstance(u_rows, list) and len(u_rows) == 3 and
            {row.get("ordinary_Z_order") for row in u_rows} == {0, 1, 2} and
            all(row.get("passed") is True and isinstance(row.get("left_node"), int) and
                isinstance(row.get("right_node"), int) for row in u_rows),
            "Rm theorem physical U rows changed")
    u_ports = [node for row in u_rows for node in (row.get("left_node"), row.get("right_node"))]
    require(len(u_ports) == 6 and len(set(u_ports)) == 6,
            "six actual scoped physical U ports required")
    return {
        "passed": True,
        "ordinary_Z_rows": 2,
        "four_scoped_Uy_ports": True,
        "three_scoped_U_rows": True,
        "six_scoped_U_ports": True,
    }


def _check_mean_certificate(mean):
    require(isinstance(mean, dict) and mean.get("passed") is True,
            "Rm leading-mean source certificate failed")
    for key, value in {
        "source_only": True,
        "source_values_selected": False,
        "numeric_values_selected": False,
        "numeric_caps_used": False,
        "assigned_numeric_values": False,
    }.items():
        if key in mean:
            require(mean.get(key) == value, ("leading-mean numeric contract", key, mean.get(key)))
    rows = mean.get("rows")
    require(isinstance(rows, list) and len(rows) == 3 and
            {row.get("ordinary_Z_order") for row in rows} == {0, 1, 2} and
            all(row.get("passed") is True and isinstance(row.get("left_node"), int) and
                isinstance(row.get("right_node"), int) for row in rows),
            "Rm leading-mean rows changed")
    left_ports = [row.get("left_node") for row in rows]
    right_ports = [row.get("right_node") for row in rows]
    require(len(set(left_ports)) == 3 and len(set(right_ports)) == 3 and
            set(left_ports).isdisjoint(right_ports),
            "leading-mean source ports must be six distinct actual nodes")
    for row in rows:
        scope = row.get("endpoint_scope")
        if isinstance(scope, dict) and "left" in scope and "right" in scope:
            require(scope == {
                "left": {"chart": "postrestore", "native_endpoint": 1},
                "right": {"chart": "Rm_patch", "native_endpoint": 1},
            }, ("leading-mean endpoint scope", row))
        else:
            left_scope = row.get("left_scope", row.get("left_endpoint_scope"))
            right_scope = row.get("right_scope", row.get("right_endpoint_scope"))
            require(left_scope == {"chart": "postrestore", "native_endpoint": 1} and
                    right_scope == {"chart": "Rm_patch", "native_endpoint": 1},
                    ("leading-mean endpoint scope", row))
    return {"passed": True, "rows": 3, "six_scoped_leading_mean_ports": True}


def _check_recovery(recovery):
    require(isinstance(recovery, dict) and recovery.get("passed") is True,
            "actual Rm trace recovery failed")
    rows = recovery.get("rows")
    require(isinstance(rows, list) and len(rows) == 17, "Rm recovery row count changed")
    expected_counts = {"U": 3, "m": 3, "m_y": 3, "U_y": 2, "Q": 2, "Q_y": 2, "R_Q_y": 2}
    counts = {name: 0 for name in expected_counts}
    for row in rows:
        require(isinstance(row, dict) and row.get("quantity") in counts and
                row.get("ordinary_Z_order") in range(expected_counts[row["quantity"]]) and
                isinstance(row.get("left_node"), int) and isinstance(row.get("right_node"), int) and
                isinstance(row.get("common_function_sha256"), str), ("recovery row", row))
        counts[row["quantity"]] += 1
    require(counts == expected_counts, ("recovery quantities", counts))
    ports = recovery.get("authorized_physical_Uy_ports")
    require(isinstance(ports, list) and len(ports) == 4 and len(set(ports)) == 4,
            "four actual scoped physical Uy ports required")
    factor_rows = recovery.get("physical_U_V_N_factor_rows")
    require(isinstance(factor_rows, list) and len(factor_rows) == 6,
            "six physical U= Pstar*V_N factor rows required")
    factor_orders = {"left": set(), "right": set()}
    factor_u_nodes = []
    factor_vn_nodes = []
    for row in factor_rows:
        require(isinstance(row, dict) and row.get("side") in factor_orders and
                row.get("ordinary_Z_order") in range(3) and
                isinstance(row.get("U_node"), int) and
                isinstance(row.get("V_N_node"), int) and
                isinstance(row.get("Pstar_node"), int) and
                row.get("passed") is True,
                ("physical U= Pstar*V_N factor row", row))
        factor_orders[row["side"]].add(row["ordinary_Z_order"])
        factor_u_nodes.append(row["U_node"])
        factor_vn_nodes.append(row["V_N_node"])
    require(factor_orders == {"left": {0, 1, 2}, "right": {0, 1, 2}} and
            len(set(factor_u_nodes)) == 6 and len(set(factor_vn_nodes)) == 6,
            "physical U= Pstar*V_N factor rows must cover six scoped ports")
    leading_rows = recovery.get("leading_Qy_input_rows")
    require(isinstance(leading_rows, list) and len(leading_rows) == 8,
            "leading Qy input row count changed")
    expected_leading_counts = {"V": 3, "m": 3, "V_y": 2}
    leading_counts = {key: 0 for key in expected_leading_counts}
    for row in leading_rows:
        require(isinstance(row, dict) and row.get("quantity") in leading_counts and
                row.get("ordinary_Z_order") in range(expected_leading_counts[row["quantity"]]) and
                isinstance(row.get("left_node"), int) and isinstance(row.get("right_node"), int),
                ("leading Qy input row", row))
        leading_counts[row["quantity"]] += 1
    require(leading_counts == expected_leading_counts,
            ("leading Qy input quantities", leading_counts))
    require(recovery.get("original_leading_Qy_subtraction_common") is True,
            "original leading Qy subtraction identity changed")
    consumed = recovery.get("all_source_ports_consumed")
    require(isinstance(consumed, list) and len(consumed) == 16 and len(set(consumed)) == 16 and
            all(isinstance(value, int) for value in consumed),
            "all 16 actual scoped source ports must be consumed exactly once")
    for key in ("original_graph_not_mutated", "endpoint_scope_restricted",
                "ordinary_Z2_U_and_mean_retained", "source_theorem_not_numeric_accuracy"):
        require(recovery.get(key) is True, ("recovery contract", key))
    return {
        "passed": True,
        "rows": 17,
        "four_scoped_Uy_ports": True,
        "six_scoped_U_ports": True,
        "six_physical_U_V_N_factor_rows": True,
        "six_scoped_leading_mean_ports": True,
        "leading_Qy_input_rows": 8,
        "all_source_ports_consumed": 16,
    }


def _check_seam(seam):
    require(isinstance(seam, dict) and seam.get("passed") is True, "Rm source seam missing")
    expected = {
        "left": "postrestore",
        "right": "Rm_patch",
        "left_native_endpoint": 1,
        "right_native_endpoint": 1,
        "ordinary_Z_orders": [0, 1],
        "physical_R_once": True,
        "global_fast_phase_not_reset": True,
        "incoming_histories_retained": True,
        "no_caps_used": True,
        "original_leading_Qy_subtraction_common": True,
    }
    for key, value in expected.items():
        require(seam.get(key) == value, ("source seam", key, seam.get(key)))
    inlet_scope = seam.get("inlet_source_seam_scope")
    require(isinstance(inlet_scope, dict) and inlet_scope.get("endpoint_node") == 1 and
            inlet_scope.get("native_x") == "1" and
            inlet_scope.get("endpoint") == "Rm" and
            inlet_scope.get("native_argument") == {"kind": "exact_decimal_Rm_patch", "value": "1"},
            "Rm source proof must stay at inlet x=1 node 1")
    outlet_scope = seam.get("external_outlet_scope")
    require(isinstance(outlet_scope, dict) and outlet_scope.get("endpoint_node") == 2948 and
            outlet_scope.get("native_x") == "Rh" and outlet_scope.get("endpoint") == "Rh" and
            outlet_scope.get("native_argument") == {"kind": "exact_exp1", "value": "exp(1)"},
            "Rm external outlet scope must stay at Rh=exp(1) node 2948")
    return {"passed": True, "scope": "postrestore native1/Rm patch native1", "ordinary_Z": [0, 1]}


def _check_source_binding(binding, parent, implementation_hashes):
    require(isinstance(binding, dict) and binding.get("passed") is True,
            "numeric source binding missing")
    expected = {
        "physical_radius_power": 1,
        "ordinary_Z_orders": [0, 1],
        "ordinary_Z1_product_rule_retained": True,
        "physical_Pstar_already_in_Qy_error": True,
        "rm_patch_physical_measure_once": True,
        "rm_patch_radius_formula": "Rm*x",
        "postrestore_native1_Rm_patch_x1_seam_required": True,
        "full_original_Rm_viscous_aggregate_used_for_Rh": True,
        "source_functions_values_or_caps_not_selected": True,
        "no_error_cap_as_function_value": True,
        "no_ancestor_constructors_or_producers_run": True,
    }
    for key, value in expected.items():
        require(binding.get(key) == value, ("source binding", key, binding.get(key)))
    require(binding.get("function") ==
            "D=R*(Qy_modified-Qy_leading); long native0 minus Rm/Rh upper source",
            "Rm source function changed")

    producers = binding.get("producer_rows")
    require(isinstance(producers, list) and producers, "producer source rows missing")
    for row in producers:
        require(isinstance(row, dict), ("producer row", row))
        source = row.get("source")
        require(source in implementation_hashes and implementation_hashes[source] == row.get("source_sha256"),
                ("producer pin", source))
        require(isinstance(row.get("method"), str), ("producer method", row))

    pins = binding.get("source_pins")
    require(isinstance(pins, dict) and set(pins) == {
        "rm_m15_report", "rm_pressure_rows_report", "rm_finite_N_archive",
        "rm_finite_N_source", "rm_patch_C6_source", "rm_m15_finite_errors",
    }, "Rm source pins changed")
    for key in ("rm_m15_report", "rm_pressure_rows_report", "rm_finite_N_archive",
                "rm_finite_N_source", "rm_patch_C6_source"):
        pin = pins[key]
        require(pin.get("passed") is True and isinstance(pin.get("path"), str), ("source pin", key))
        if key == "rm_m15_report":
            require(pin.get("sha256") == RM_REPORT_SHA256, ("source pin hash", key))
        elif key == "rm_pressure_rows_report":
            require(pin.get("sha256") == RM_PRESSURE_REPORT_SHA256, ("source pin hash", key))
        elif key == "rm_finite_N_archive":
            require(pin.get("sha256") == RM_ARCHIVE_SHA256, ("source pin hash", key))
        elif key == "rm_finite_N_source":
            require(pin.get("sha256") == RM_SOURCE_SHA256, ("source pin hash", key))
        else:
            require(pin.get("sha256") == RM_C6_SOURCE_SHA256, ("source pin hash", key))
    ast_binding = pins["rm_m15_finite_errors"]
    require(ast_binding.get("passed") is True and ast_binding.get("source") == RM_M15_SOURCE and
            ast_binding.get("source_sha256") == RM_M15_SOURCE_SHA256 and
            ast_binding.get("method") == "OriginalRmM15.finite_errors" and
            ast_binding.get("complete_prefix_error_source") is True and
            ast_binding.get("source_not_executed") is True, "Rm finite_errors AST pin changed")

    graph = binding.get("actual_Rm_patch_finite_N_source_binding")
    require(isinstance(graph, dict) and graph.get("passed") is True and
            graph.get("physical_measure") == "dy=dx/x once" and
            graph.get("radius_formula") == "Rm*x" and graph.get("native_domain") == ["1", "e"],
            "actual Rm graph source binding changed")
    rows = graph.get("graph_function_rows")
    require(isinstance(rows, list) and len(rows) == 2 and
            [row.get("Qy_function_node") for row in rows] == [204500, 204504] and
            [row.get("RQy_function_node") for row in rows] == [204659, 204660] and
            all(row.get("physical_R_explicit_once") is True for row in rows),
            "Rm Qy/RQy graph rows changed")
    require(graph.get("chart_physical", {}).get("Q_y") == [204500, 204504, 204510, 204516, 204522],
            "Rm original chart Qy functions changed")
    scopes = graph.get("endpoint_scopes")
    require(isinstance(scopes, dict) and set(scopes) == {"inlet_x1", "outlet_Rh"},
            "Rm endpoint scopes must include inlet and outlet separately")
    require(scopes.get("inlet_x1") == graph.get("inlet_source_seam_scope") and
            scopes.get("inlet_x1", {}).get("coordinate_node") == graph.get("chart", {}).get("coordinate") and
            scopes.get("inlet_x1", {}).get("endpoint_node") == 1 and
            scopes.get("inlet_x1", {}).get("native_x") == "1" and
            scopes.get("inlet_x1", {}).get("endpoint") == "Rm" and
            scopes.get("inlet_x1", {}).get("native_argument") ==
            {"kind": "exact_decimal_Rm_patch", "value": "1"},
            "Rm inlet x=1 source scope changed")
    require(scopes.get("outlet_Rh") == graph.get("external_outlet_scope") and
            scopes.get("outlet_Rh", {}).get("coordinate_node") == graph.get("chart", {}).get("coordinate") and
            scopes.get("outlet_Rh", {}).get("endpoint_node") == 2948 and
            scopes.get("outlet_Rh", {}).get("native_x") == "Rh" and
            scopes.get("outlet_Rh", {}).get("endpoint") == "Rh" and
            scopes.get("outlet_Rh", {}).get("native_argument") ==
            {"kind": "exact_exp1", "value": "exp(1)"},
            "Rm Rh outlet scope changed")
    require(graph.get("boundary_Q_plus_Q_y") == [204655, 204656, 204657, 204658] and
            graph.get("boundary_R_Q_y") == [204659, 204660, 204661, 204662],
            "Rm endpoint graph handles changed")

    parent_source = parent.get("source_binding") or parent.get("numerical_source_binding")
    require(isinstance(parent_source, dict), "parent source binding missing")
    require(binding.get("inherited_parent_source_binding_sha256") == digest(parent_source),
            "parent source provenance digest changed")
    require(binding.get("inherited_parent_report_sha256") == digest(parent),
            "parent report provenance digest changed")
    return {
        "passed": True,
        "producer_rows": len(producers),
        "source_pins": len(pins),
        "physical_measure": "dy=dx/x once",
        "cap_as_value": False,
    }


def _check_rm_saved_packet(rm_report, parent_sample, Z):
    x1 = _find(rm_report.get("samples"), Z=Z, chart="actual_patch", native_x="1")
    rh = _find(rm_report.get("samples"), Z=Z, chart="actual_patch", native_x="Rh")
    terms = x1.get("m15_C0_C1_terms")
    incoming = x1.get("incoming_m15_terms")
    parent_raw = parent_sample.get("raw_m15_C0_C1_terms")
    require(isinstance(terms, list) and len(terms) == 168 and isinstance(incoming, list) and
            len(incoming) == INCOMING_COUNT and isinstance(parent_raw, list) and
            len(parent_raw) == RAW_COUNT and terms[:INCOMING_COUNT] == incoming and
            incoming == parent_raw[:INCOMING_COUNT],
            "Rm x=1 incoming 156 transfer changed")
    require(x1.get("exact_empty_increment") is True, "Rm x=1 empty increment changed")
    expected_tail = [
        "Rm_patch_constant_background_integral",
        "Rm_patch_R2F_upper_t1",
        "Rm_patch_R2F_lower_t1",
        "Rm_patch_G_integral_t1",
        "Rm_patch_R2F_upper_t2",
        "Rm_patch_R2F_lower_t2",
        "Rm_patch_G_integral_t2",
        "Rm_patch_RQy_upper_t1",
        "Rm_patch_RQy_lower_t1",
        "Rm_patch_finite_N_endpoint_error",
        "Rm_patch_finite_N_volume_error",
        OLD_RM_TERM,
    ]
    tail = terms[INCOMING_COUNT:]
    require([term.get("name") for term in tail] == expected_tail and len(tail) == 12,
            "Rm x=1 twelve empty tail names changed")
    for term in tail:
        coefficients = term.get("coefficient_C0_C1")
        require(isinstance(coefficients, list) and len(coefficients) == 2 and
                all(row.get("approximation") in ("0.0", "0") and
                    row.get("absolute_error_upper") in ("0.0", "0") and
                    row.get("raw_directed_error_exact_mpf_tuple") == [0, 0, 0, 0] and
                    row.get("source_enclosure_exact_mpf_tuples") == [[0, 0, 0, 0], [0, 0, 0, 0]]
                    for row in coefficients),
                ("Rm x=1 empty tail", term.get("name")))
    raw = parent_sample.get("raw_m15_C0_C1_terms")
    rh_terms = rh.get("m15_C0_C1_terms")
    require(isinstance(rh_terms, list) and len(rh_terms) == 168 and
            rh_terms[:INCOMING_COUNT] == incoming and
            rh_terms[:INCOMING_COUNT] == parent_raw[:INCOMING_COUNT] and
            rh_terms[RM_VISCOUS_INDEX] == raw[RM_VISCOUS_INDEX],
            "saved Rm Rh packet must retain the raw 156-term prefix and index 167")
    parent_selected = parent_sample.get("selected_N_evidence", {})
    require(isinstance(parent_selected, dict) and
            isinstance(parent_selected.get("postrestore_native1"), dict) and
            x1.get("actual_selected_integer") == parent_selected["postrestore_native1"] and
            rh.get("actual_selected_integer") == parent_selected["postrestore_native1"],
            "saved Rm selected N must equal parent postrestore selected N evidence")
    aggregate = raw[RM_VISCOUS_INDEX]
    require(aggregate.get("name") == OLD_RM_TERM and
            aggregate.get("actual_function_binding", {}).get("chart") == "actual_patch" and
            aggregate.get("actual_function_binding", {}).get("incoming156_terms_retained_once") is True and
            aggregate.get("actual_function_binding", {}).get("error_enclosure_only") is True,
            "raw Rm viscous source binding changed")
    logs = aggregate.get("log_scale_terms")
    require(isinstance(logs, list) and len(logs) == 1 and
            logs[0].get("positive_error_bound_atom_only") is True and
            logs[0].get("source_value_selected") is False,
            "raw Rm viscous error atom changed")
    return {
        "x1_empty_tail_count": 12,
        "incoming156_transfer": True,
        "aggregate_raw_index": RM_VISCOUS_INDEX,
        "aggregate_matches_Rh": True,
        "cap_as_value": False,
    }


def _check_trace_chain(chain):
    require(isinstance(chain, dict), "Rm signed trace chain missing")
    require(chain.get("signed_coefficients") == {"long_native0": 1, "Rh": -1},
            "Rm signed endpoint coefficients changed")
    require(chain.get("external_endpoints") == {"lower": "long native0", "upper": "Rm/Rh"},
            "Rm external endpoint scope changed")
    require(chain.get("function") == "D(long native0)-D(RmRh)", "Rm trace function changed")
    return {"passed": True, "external": "long native0 -> Rm/Rh", "cap_as_value": False}


def _check_budget(sample):
    budget = sample.get("external_endpoint_budget_witness", {})
    require(isinstance(budget, dict), "Rm external budget missing")
    witness = budget.get("own_longreshape0", {})
    require(isinstance(witness, dict), "own long native0 budget missing")
    own = witness.get("radius_weighted_C0_C1_log_bounds")
    require(isinstance(own, list) and len(own) == 2 and all(value is not None for value in own),
            "own long native0 C0/Z1 rows changed")
    aggregate = budget.get("original_Rm_viscous_aggregate_167", {})
    require(isinstance(aggregate, dict) and aggregate.get("aggregate_raw_index") == RM_VISCOUS_INDEX and
            aggregate.get("aggregate_name") == OLD_RM_TERM and
            aggregate.get("aggregate_matches_saved_Rm_Rh") is True and
            aggregate.get("aggregate_dominates_Rh") is True and
            aggregate.get("incoming156_transfer") is True and
            aggregate.get("empty_tail_count") == 12,
            "Rm viscous aggregate budget changed")
    atom = aggregate.get("aggregate_log_atom")
    inputs = budget.get("positive_union_inputs")
    require(isinstance(atom, dict) and inputs == own + [atom], "Rm positive union inputs changed")
    union = _check_positive_union(inputs, budget.get("positive_union"), "long0 plus Rm167 union")
    require(budget.get("postrestore_upper_separate_budget_claimed") is False and
            budget.get("entire_Rm_viscous_aggregate_dominates_Rh") is True and
            budget.get("no_cap_subtraction_value_midpoint_or_overlap") is True,
            "Rm budget source contract changed")
    return union


def _check_ports(ports, source_binding):
    require(isinstance(ports, list) and len(ports) == 1, "Rm physical endpoint ports changed")
    row = ports[0]
    graph = source_binding.get("actual_Rm_patch_finite_N_source_binding", {})
    graph_rows = graph.get("graph_function_rows", [])
    require(row.get("chart") == "actual_patch" and row.get("side") == "upper" and
            row.get("endpoint_node") == 2948 and
            row.get("native_argument") == {"kind": "exact_exp1", "value": "exp(1)"} and
            row.get("definition") == "D=R*(Qy_modified-Qy_leading)" and
            row.get("ordinary_Z_orders") == [0, 1] and
            row.get("physical_R_explicit_once_in_both_Z_rows") is True and
            row.get("modified_Qy_nodes") == [entry.get("Qy_function_node") for entry in graph_rows] and
            row.get("modified_RQy_nodes") == [entry.get("RQy_function_node") for entry in graph_rows],
            "Rm endpoint port scope changed")
    require(row.get("endpoint_scope") == {
                "chart": "actual_patch",
                "coordinate_node": graph.get("chart", {}).get("coordinate"),
                "endpoint_node": 2948,
                "native_argument": {"kind": "exact_exp1", "value": "exp(1)"},
                "native_x": "Rh",
                "endpoint": "Rh",
                "external_port_role": "upper signed endpoint",
            }, "Rm endpoint coordinate changed")
    return {"passed": True, "ports": 1, "ordinary_Z": [0, 1]}


def _check_sample(sample, parent_sample, rm_report, source_binding):
    require(sample.get("source_family") == parent_sample.get("source_family"), "sample source family changed")
    raw = sample.get("raw_m15_C0_C1_terms")
    parent_raw = parent_sample.get("raw_m15_C0_C1_terms")
    require(isinstance(raw, list) and len(raw) == RAW_COUNT and raw == parent_raw,
            "raw 213 endpoint state changed")
    previous = parent_sample.get("correlated_reduced_m15_C0_C1_terms")
    reduced = sample.get("correlated_reduced_m15_C0_C1_terms")
    require(isinstance(previous, list) and len(previous) == PREVIOUS_REDUCED_COUNT and
            isinstance(reduced, list) and len(reduced) == NEW_REDUCED_COUNT,
            "197 to 196 reduced counts changed")
    old_rm = [term for term in previous if term.get("name") == OLD_RM_TERM]
    old_external = [term for term in previous if term.get("name") == OLD_EXTERNAL_TERM]
    new_external = [term for term in reduced if term.get("name") == NEW_EXTERNAL_TERM]
    require(len(old_rm) == 1 and len(old_external) == 1 and len(new_external) == 1,
            "Rm/external term identities changed")
    unchanged_old = [term for term in previous if term.get("name") not in {OLD_RM_TERM, OLD_EXTERNAL_TERM}]
    unchanged_new = [term for term in reduced if term.get("name") != NEW_EXTERNAL_TERM]
    require(len(unchanged_old) == UNCHANGED_COUNT and unchanged_new == unchanged_old,
            "other 195 reduced entries changed")
    require(sample.get("old_budget_terms_retained_as_provenance") == [raw[RM_VISCOUS_INDEX], old_rm[0]],
            "Rm raw/reduced provenance changed")
    require(sample.get("previous_external_replacement_retained_as_provenance") == old_external[0],
            "old long/reference external provenance changed")
    require(sample.get("external_endpoint_budget_witness", {}).get("source_binding") == source_binding,
            "per-sample Rm source binding changed")
    require(sample.get("independent_target") == parent_sample.get("independent_target"),
            "independent target changed")

    packet_check = _check_rm_saved_packet(rm_report, parent_sample, sample.get("Z"))
    budget_check = _check_budget(sample)
    ports_check = _check_ports(sample.get("actual_signed_function_ports"), source_binding)
    trace_check = _check_trace_chain(sample.get("signed_trace_chain"))
    contract = sample.get("contract", {})
    for key, value in {
        "raw_213_verbatim": True,
        "other_195_verbatim": True,
        "independent_target_unchanged": True,
        "reserved_added_this_stage": 0,
        "ordinary_Z1_and_physical_R_once": True,
        "leading_subtraction_retained": True,
        "long_native0_and_RmRh_external_source_functions_retained": True,
        "parent_long_native0_uses_own_saved_budget": True,
        "entire_original_Rm_viscous_aggregate_167_used": True,
        "aggregate_dominates_Rh": True,
        "raw156_incoming_transfer_exact": True,
        "actual_x1_empty_twelve_tails": True,
        "no_error_cap_as_value": True,
        "all_other_volume_width_numerical_errors_retained": True,
        "ancestors_not_reexecuted": True,
        "source_bound_signed_trace_contraction": True,
    }.items():
        require(contract.get(key) == value, ("sample contract", key, contract.get(key)))
    limits = sample.get("limits", {})
    require(limits.get("long_Rm_Rh_internal_trace_resolved") is True and
            all(value is False for key, value in limits.items() if key != "long_Rm_Rh_internal_trace_resolved"),
            "sample limits changed")
    return {
        "Z": sample.get("Z"),
        "raw_213": True,
        "reduced_196": True,
        "unchanged_195": True,
        "Rm_packet": packet_check,
        "budget": budget_check,
        "ports": ports_check,
        "trace_chain": trace_check,
        "target_unchanged": True,
        "cap_as_value": False,
    }


def check_packet(report, parent, rm_report, pressure_report):
    require(report.get("source_family") == parent.get("source_family"), "source family changed")
    require(report.get("parent_report") == PARENT and report.get("parent_report_sha256") == PARENT_SHA,
            "parent identity changed")
    flat_check = _check_flat(report.get("compact_flat_and_history_certificate"))
    native_check = _check_native_guards(report.get("native_source_guards"))
    theorem_check = _check_theorem(report.get("source_derivative_theorem"))
    mean_check = _check_mean_certificate(report.get("leading_mean_certificate"))
    recovery_check = _check_recovery(report.get("actual_trace_recovery"))
    seam_check = _check_seam(report.get("source_seam_proof"))
    implementation_hashes = report.get("implementation_hashes", {})
    require(MEAN_SOURCE in implementation_hashes and
            isinstance(implementation_hashes.get(MEAN_SOURCE), str) and
            len(implementation_hashes[MEAN_SOURCE]) == 64,
            "leading-mean source hash missing")
    source = report.get("numerical_source_binding", {})
    source_check = _check_source_binding(source, parent, implementation_hashes)
    samples = report.get("samples")
    require(isinstance(samples, list) and tuple(row.get("Z") for row in samples) == SAVED_Z,
            "saved Z set changed")
    parent_by_Z = {row.get("Z"): row for row in parent.get("samples", [])}
    sample_checks = []
    for sample in samples:
        require(sample.get("Z") in parent_by_Z, ("unknown saved Z", sample.get("Z")))
        sample_checks.append(_check_sample(sample, parent_by_Z[sample.get("Z")], rm_report, source))
    require(report.get("contract") == {
        "accepted_ancestor_constructors_not_run": True,
        "source_functions_not_cap_values": True,
        "raw_213_and_independent_target_retained": True,
    }, "top contract changed")
    limits = report.get("limits", {})
    require(isinstance(limits, dict) and all(value is False for value in limits.values()),
            "top limits changed")
    require(isinstance(pressure_report.get("checks"), dict) and pressure_report["checks"].get("passed") is True,
            "pressure report checks changed")
    return {
        "flat": flat_check,
        "native": native_check,
        "theorem": theorem_check,
        "leading_mean": mean_check,
        "recovery": recovery_check,
        "source_seam": seam_check,
        "source_binding": source_check,
        "samples": sample_checks,
        "all_open_limits": True,
        "cap_as_value": False,
    }


def run():
    started = time.perf_counter()
    from check_original_Rm_inlet_m15 import check_fixture

    fixture = check_fixture()
    fixture_check = _check_fixture(fixture)
    owner = OriginalRmInletM15(precision=150)
    packet = owner.build()
    code = dict(owner.code_hashes)
    code[FIXTURE_NAME] = sha(ROOT / FIXTURE_NAME)
    code[RUNNER_NAME] = sha(ROOT / RUNNER_NAME)
    report = dict(packet, input_hashes=dict(owner.hashes), implementation_hashes=code)
    parent, rm_report, pressure = load_reports()
    packet_check = check_packet(report, parent, rm_report, pressure)
    report["checks"] = {
        "passed": True,
        "independent_fixture": fixture,
        "fixture_invariants": fixture_check,
        "packet_invariants": packet_check,
        "saved_Z": list(SAVED_Z),
        "raw_213_same_long_reference_parent": True,
        "correlated_reduced_196": True,
        "unchanged_other_sectors_195": True,
        "own_long0_and_entire_Rm167_budget": True,
        "source_functions_not_caps": True,
        "raw_Rm156_incoming_and_x1_empty_tails": True,
        "production_is_new_owner_and_fixture_once": True,
    }
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    return report


def inspect(path: Path):
    """Inspect hashes and serialized packet shape without owner or fixture."""
    data = read_json(path)
    require(data.get("checks", {}).get("passed") is True, "saved Rm inlet report is not passed")
    check_hash_group(ROOT, data.get("implementation_hashes"), "implementation")
    check_hash_group(DIRECTORY, data.get("input_hashes"), "input")
    parent, rm_report, pressure = load_reports()
    packet_check = check_packet(data, parent, rm_report, pressure)
    stored_fixture = data.get("checks", {}).get("independent_fixture")
    fixture_check = _check_fixture(stored_fixture)
    expected_checks = {
        "saved_Z": list(SAVED_Z),
        "raw_213_same_long_reference_parent": True,
        "correlated_reduced_196": True,
        "unchanged_other_sectors_195": True,
        "own_long0_and_entire_Rm167_budget": True,
        "source_functions_not_caps": True,
        "raw_Rm156_incoming_and_x1_empty_tails": True,
    }
    for key, expected in expected_checks.items():
        require(data.get("checks", {}).get(key) == expected, ("saved check", key))
    return {
        "checks": data["checks"],
        "packet_invariants": packet_check,
        "fixture_invariants": fixture_check,
        "execution_seconds": data.get("execution_seconds"),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write-report", action="store_true")
    modes.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    path = DIRECTORY / NAME
    if args.write_report:
        data = run()
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n",
                             encoding="utf-8", newline="\n")
        temporary.replace(path)
    else:
        data = inspect(path)
    print(json.dumps({"checks": data.get("checks", {}),
                      "execution_seconds": data.get("execution_seconds")}, indent=2), flush=True)
