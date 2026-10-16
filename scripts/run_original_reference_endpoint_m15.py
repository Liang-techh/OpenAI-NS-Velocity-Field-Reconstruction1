"""Run and inspect the signed reference/restoration/postrestore endpoint packet.

Production constructs the new owner once, runs the independent three-piece
fixture once, and serializes the three saved axial packets.  Inspection only
reads the saved report, its pinned inputs, and serialized invariants; it does
not construct the owner or rerun the fixture.
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
from openai_ns_reconstruction.original_reference_endpoint_m15 import (  # noqa: E402
    ACTIVE,
    ACTIVE_SHA,
    NATIVE,
    NATIVE_SHA,
    OLD_COMBINED,
    REPLACEMENT,
    OriginalReferenceEndpointM15,
)


NAME = PREFIX + "compliant_current_original_reference_endpoint_m15_numeric.json"
REFERENCE_NAME = PREFIX + "compliant_current_original_reference_m15_numeric.json"
RESTORE_NAME = PREFIX + "compliant_current_original_restore_m15_numeric.json"
FIXTURE_NAME = "scripts/check_original_reference_endpoint_m15.py"
RUNNER_NAME = "scripts/run_original_reference_endpoint_m15.py"
SAVED_Z = ("-.3", "0", ".3")
RAW_COUNT = 213
PREVIOUS_REDUCED_COUNT = 199
NEW_REDUCED_COUNT = 198
UNCHANGED_COUNT = 197


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def digest(value) -> str:
    # Match the source owner's canonical packet digest.
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def require(condition, message):
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


def load_accepted_reports():
    native_path = DIRECTORY / NATIVE
    active_path = DIRECTORY / ACTIVE
    require(sha(native_path) == NATIVE_SHA, "accepted native parent changed")
    require(sha(active_path) == ACTIVE_SHA, "accepted active parent changed")
    native = read_json(native_path)
    active = read_json(active_path)
    require(native.get("checks", {}).get("passed") is True, "native parent is not passed")
    require(active.get("checks", {}).get("passed") is True, "active parent is not passed")
    reference = read_json(DIRECTORY / REFERENCE_NAME)
    restore = read_json(DIRECTORY / RESTORE_NAME)
    require(reference.get("checks", {}).get("passed") is True, "reference source packet is not passed")
    require(restore.get("checks", {}).get("passed") is True, "restore source packet is not passed")
    return native, active, reference, restore


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
    """Recompute the saved positive log union from exact MPF tuples."""
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


def _find_sample(report, Z, chart, coordinate_value):
    expected_kind = "exact_decimal_" + chart
    rows = report.get("samples", [])
    for row in rows:
        coordinate = row.get("coordinate", {})
        if (row.get("Z") == Z and row.get("chart") == chart and
                coordinate == {"kind": expected_kind, "value": coordinate_value}):
            return row
    raise AssertionError(("missing endpoint sample", Z, chart, coordinate_value))


def _find_reference_pressure_row(reference, Z):
    for row in reference.get("pressure_row_samples", []):
        if row.get("Z") == Z and row.get("coordinate") == {
            "kind": "exact_decimal_reference", "value": "0"
        }:
            return row
    raise AssertionError(("missing reference0 pressure packet", Z))


def _check_fixture(fixture):
    require(isinstance(fixture, dict) and fixture.get("passed") is True,
            "independent endpoint fixture failed")
    for key, expected in {
        "synthetic_parameters_only": True,
        "original_source_values_evaluated": False,
        "finite_N_slow_fast": True,
        "unequal_J_naturality": True,
        "common_physical_seams": 2,
        "nonzero_leading_subtraction": True,
        "physical_R_once": True,
        "physical_y_and_ordinary_Z1_independent_differentiation": True,
        "both_external_endpoints_nonzero": True,
    }.items():
        require(fixture.get(key) is expected, ("fixture flag", key, fixture.get(key)))
    controls = fixture.get("controls")
    require(isinstance(controls, dict) and all(value is True for value in controls.values()),
            ("fixture negative controls", controls))
    trace = fixture.get("signed_trace")
    require(isinstance(trace, dict), "fixture signed trace missing")
    require(trace.get("piece_orientation") == "D(lower)-D(upper)", "fixture orientation changed")
    require(trace.get("delta_definition") == "D=R*(Qy_modified-Qy_leading)",
            "fixture delta definition changed")
    require(trace.get("external_remainder") == "D(reference_native0)-D(postrestore_native1)",
            "fixture external trace changed")
    require(trace.get("internal_seams_cancelled") == [
        "reference_native1=restoration_native0",
        "restoration_native1=postrestore_native0",
    ], "fixture seam list changed")
    return {
        "passed": True,
        "synthetic_parameters_only": True,
        "common_physical_seams": 2,
        "both_external_endpoints_nonzero": True,
        "negative_controls": len(controls),
        "signed_trace_proved": True,
    }


def _check_source_binding(binding, implementation_hashes):
    require(isinstance(binding, dict) and binding.get("passed") is True,
            "source binding missing")
    expected = {
        "function": "D=R*(Qy_modified-Qy_original_leading); chart contribution=D(lower)-D(upper)",
        "physical_radius_power": 1,
        "ordinary_Z_orders": [0, 1],
        "positive_sum_domination": "exp(max(log_bounds)+log(count)) bounds the sum of the nonnegative magnitudes",
        "ordinary_Z1_product_rule_retained": True,
        "physical_Pstar_already_in_Qy_error": True,
        "slow_B_y_over_N_and_fast_B_phi_both_retained": True,
        "separate_inlet_uses_own_reference0_pressure_packet": True,
        "old_reference_lower_packet_is_longreshape1_provenance_only": True,
        "no_longreshape_reference_alias_authorized": True,
        "postrestore_entire_original_155_cap_retained": True,
        "inherited_pressure_bound_validity_not_independently_recertified": True,
    }
    for key, value in expected.items():
        require(binding.get(key) == value, ("source binding", key, binding.get(key)))
    producers = binding.get("producer_rows")
    require(isinstance(producers, list) and len(producers) == 5, "producer source rows changed")
    for row in producers:
        require(isinstance(row, dict), ("producer row", row))
        source = row.get("source")
        require(source in implementation_hashes and implementation_hashes[source] == row.get("source_sha256"),
                ("producer pin", source))
        require(isinstance(row.get("method"), str) and isinstance(row.get("exact_source_statements"), list),
                ("producer row shape", row))
    return {"passed": True, "producer_rows": len(producers), "cap_as_value": False}


def _check_source_seams(seams):
    require(isinstance(seams, dict) and seams.get("passed") is True, "source seam packet missing")
    expected = {
        "delta_definition": "R*(Qy_modified-Qy_original_leading)",
        "leading_Qy_recovery": "Pstar*[2Z*Vy-(1-delta)*Z*(V-m)-(1-Z^2)*(Vz-mz)]/(1-delta*Z^2)",
        "leading_inputs_common_through_Z2": True,
        "physical_radius_common_and_Z_independent": True,
        "fast_phase_not_reset": True,
        "source_functions_not_error_caps": True,
    }
    for key, value in expected.items():
        require(seams.get(key) == value, ("source seam", key, seams.get(key)))
    first = seams.get("reference_restoration_modified_RQy_rows")
    require(isinstance(first, list) and len(first) == 2, "reference/restoration RQy proof rows changed")
    for row in first:
        require(isinstance(row, dict) and row.get("quantity") == "R_Q_y" and
                row.get("ordinary_Z_order") in {0, 1}, ("reference seam row", row))
        require(isinstance(row.get("left_node"), int) and isinstance(row.get("right_node"), int),
                ("reference seam node IDs", row))
        require(isinstance(row.get("common_endpoint_sha256"), str) and
                len(row["common_endpoint_sha256"]) == 64, ("reference seam digest", row))
    second = seams.get("restoration_postrestore_certificate")
    require(isinstance(second, dict) and second.get("passed") is True,
            "restoration/postrestore source proof missing")
    require((second.get("left"), second.get("right")) == ("restoration", "postrestore"),
            "restoration/postrestore seam changed")
    for key in ("global_phase_not_reset", "actual_modified_and_original_leading_traces_equal",
                "independent_of_interval_caps", "physical_R_factor_in_both_ordinary_Z_rows"):
        require(second.get(key) is True, ("outgoing seam flag", key))
    require(second.get("actual_signed_centers_evaluated") is False,
            "outgoing cap was used as a signed value")
    require(second.get("native_Jacobians") == [1, 1], "outgoing native measures changed")
    return {"passed": True, "reference_restoration_rows": 2,
            "restoration_postrestore_proof": True, "caps_are_not_functions": True}


def _check_ports(ports):
    require(isinstance(ports, list) and len(ports) == 6, "six endpoint function ports required")
    seen = set()
    for port in ports:
        require(isinstance(port, dict), ("function port", port))
        chart, side = port.get("chart"), port.get("side")
        seen.add((chart, side))
        require(chart in {"reference", "restoration", "postrestore"} and side in {"lower", "upper"},
                ("function port scope", port))
        require(isinstance(port.get("coordinate_node"), int) and isinstance(port.get("endpoint_node"), int),
                ("function port endpoint", port))
        require(port.get("endpoint_scope") == {
            "coordinate_node": port["coordinate_node"],
            "endpoint_node": port["endpoint_node"],
        }, ("function port endpoint scope", port))
        require(isinstance(port.get("physical_radius_node"), int) and
                isinstance(port.get("native_Jacobian_node"), int), ("function port geometry", port))
        require(isinstance(port.get("modified_RQy_nodes"), list) and
                len(port["modified_RQy_nodes"]) == 2, ("function port RQy rows", port))
        require(isinstance(port.get("leading_normalized_V_nodes"), list) and
                len(port["leading_normalized_V_nodes"]) == 3, ("function port V rows", port))
        require(isinstance(port.get("leading_normalized_mean_nodes"), list) and
                len(port["leading_normalized_mean_nodes"]) == 3, ("function port mean rows", port))
        require(isinstance(port.get("leading_normalized_Vy_nodes"), list) and
                len(port["leading_normalized_Vy_nodes"]) == 2, ("function port Vy rows", port))
        require(port.get("definition") == "D=modified_RQy-R*Pstar*leading_normalized_Qy",
                "function port definition changed")
        require(isinstance(port.get("radius_and_measure_Z_free_dependency_nodes"), list),
                "radius dependency proof missing")
        require(port.get("physical_R_explicit_once_in_both_Z_rows") is True,
                "physical R factor proof missing")
        require(port.get("ordinary_Z_orders") == [0, 1], "ordinary Z orders changed")
        require(port.get("signed_coefficient") == (1 if side == "lower" else -1),
                ("endpoint sign", port))
    require(seen == {(chart, side) for chart in ("reference", "restoration", "postrestore")
                     for side in ("lower", "upper")}, "endpoint port coverage changed")
    return {"ports": 6, "all_charts_and_sides_scoped": True, "physical_R_once": True}


def _check_reference0_graph_envelope(binding, native, implementation_hashes=None):
    """Validate the saved reference-0 source graph certificate verbatim."""
    require(isinstance(binding, dict) and binding.get("passed") is True,
            "reference0 graph-envelope certificate missing")
    require(binding.get("method") ==
            "source_bound_linear_Qy_envelope_at_actual_reference_native0",
            "reference0 graph-envelope method changed")
    scope = binding.get("endpoint_scope")
    require(isinstance(scope, dict) and
            isinstance(scope.get("coordinate_node"), int) and
            isinstance(scope.get("endpoint_node"), int),
            "reference0 graph-envelope endpoint scope missing")
    require(scope == {
        "coordinate_node": native["corrected_history_seam"]["endpoint_scope"]["coordinate_handles"][0],
        "endpoint_node": 0,
    }, "actual reference native0 scope changed")
    rows = binding.get("graph_function_rows")
    require(isinstance(rows, list) and len(rows) == 2,
            "reference0 graph-envelope rows changed")
    for order, row in enumerate(rows):
        require(isinstance(row, dict) and row.get("ordinary_Z_order") == order,
                ("reference0 graph-envelope row order", row))
        for key in ("Qy_function_node", "RQy_function_node"):
            require(isinstance(row.get(key), int),
                    ("reference0 graph-envelope node", key, row))
        recovered = native["actual_modified_trace_recovery"]["rows"]
        for key, quantity in (("Qy_function_node", "Q_y"), ("RQy_function_node", "R_Q_y")):
            actual = next(r["left_node"] for r in recovered
                          if r["quantity"] == quantity and r["ordinary_Z_order"] == order)
            require(row[key] == actual, ("actual reference function port changed", key, order))
        for key in ("Qy_at_reference0_sha256", "RQy_at_reference0_sha256"):
            value = row.get(key)
            require(isinstance(value, str) and len(value) == 64,
                    ("reference0 graph-envelope node digest", key, row))
    guards = native.get("native_derivative_source_guards", {}).get("checks", {})
    native_pressure = guards.get("source_ast", {}).get("rows", {}).get("pressure_formula")
    native_radial = guards.get("source_ast", {}).get("rows", {}).get("radial_formula")
    native_slow_fast = guards.get("slow_B_to_total_Vy_to_physical_Uy", {}).get(
        "charts", {}).get("reference")
    require(binding.get("actual_archived_pressure_formula_certificate") == native_pressure,
            "reference0 pressure native guard certificate changed")
    require(binding.get("actual_archived_radial_formula_certificate") == native_radial,
            "reference0 radial native guard certificate changed")
    require(binding.get("actual_reference_slow_fast_and_Pstar_certificate") == native_slow_fast,
            "reference0 slow/fast native guard certificate changed")
    require(binding.get("reference_wrapper_source") ==
            "src/openai_ns_reconstruction/original_reference_pressure_rows.py",
            "reference0 wrapper source changed")
    require(binding.get("numeric_Qy_coefficient_source") ==
            "src/openai_ns_reconstruction/original_first_micro_pressure_rows.py",
            "reference0 Qy coefficient source changed")
    require(binding.get("exact_coefficient_geometry_statements") == [
        "z=J.variable(ctx,ctx.mpf(Z),order)", "L=1-z*z*delta", "d=1-z*z",
    ], "reference Qy coefficient geometry changed")
    require(isinstance(binding.get("reference_wrapper_sha256"), str) and
            len(binding["reference_wrapper_sha256"]) == 64,
            "reference0 wrapper source digest missing")
    if implementation_hashes is not None:
        require(implementation_hashes.get(binding["reference_wrapper_source"]) ==
                binding["reference_wrapper_sha256"],
                "reference0 wrapper source pin changed")
    require(binding.get("exact_reference_owner_statements") == [
        "self.chart,self.derivative_chart=chart,self.OWNERS[chart]",
        "view._covers=lambda Z,q:self.reference._covers(Z,self.chart,q)",
        "self.loop=json.loads(path.read_bytes())['current_actual_loop_jet_log_bounds_by_chart'][self.derivative_chart]",
    ], "reference0 owner AST statements changed")
    require(binding.get("exact_numeric_Qy_coefficients") ==
            "{'U_y':2*z/L,'U':-z*(1-delta)/L,'m':z*(1-delta)/L,'U_Z':-d/L,'m_Z':d/L}",
            "reference0 numeric Qy coefficients changed")
    require(binding.get("ordinary_Z_product_rule") ==
            "ordinary=coef[j]*(math.factorial(j)*math.comb(k,j))",
            "reference0 ordinary-Z rule changed")
    require(binding.get("own_reference_generic_derivative_owner") == "inner_reference",
            "reference0 generic derivative owner changed")
    for key in ("delta_inputs_from_same_saved_reference_pressure_packet",
                "physical_R_applied_after_Qy_error",
                "source_field_values_or_caps_not_selected",
                "source_function_not_identified_with_longreshape1",
                "inherited_generic_and_history_enclosures_not_recomputed"):
        require(binding.get(key) is True, ("reference0 envelope contract", key))
    return {
        "passed": True,
        "endpoint_scope": copy.deepcopy(scope),
        "rows": 2,
        "actual_node_ids": True,
        "native_guard_certificates_same_exact": True,
    }


def _check_trace_chain(chain):
    require(chain == {
        "signed_coefficients": {"reference0": 1, "Rm": -1},
        "cancelled_interfaces": ["reference1", "restoration1"],
    }, ("signed trace chain", chain))
    return {"signed_coefficients": {"reference0": 1, "Rm": -1},
            "cancelled_interfaces": ["reference1", "restoration1"]}


def _check_budget(sample, raw, reference, Z):
    witness = sample.get("external_endpoint_budget_witness")
    require(isinstance(witness, dict), ("budget witness", Z))
    own = witness.get("own_reference0")
    require(isinstance(own, dict), ("own reference0 budget", Z))
    for key, value in {
        "source_report": REFERENCE_NAME,
        "physical_R_once": True,
        "physical_Pstar_already_in_error": True,
        "extra_epsilon": False,
        "source_function_values_evaluated": False,
        "old_longreshape_lower_packet_not_used": True,
    }.items():
        require(own.get(key) == value, ("own reference budget", key, own.get(key)))
    pressure = _find_reference_pressure_row(reference, Z)
    require(own.get("source_packet_sha256") == digest(pressure), "reference0 packet digest changed")
    require(own.get("coordinate") == pressure["coordinate"], "reference0 pressure coordinate changed")
    require(own.get("selected_integer") == pressure.get("selected_integer"),
            "reference0 selected N changed")
    errors = own.get("actual_own_Qy_error_rows")
    radius = own.get("actual_own_log_radius")
    weighted = own.get("radius_weighted_C0_C1_log_bounds")
    require(errors == pressure.get("physical_pressure_rows", {}).get("Q_y", {}).get(
        "finite_N_absolute_error_log_bounds"),
            ("own Qy source rows changed", Z))
    require(radius == pressure.get("geometry", {}).get("log_physical_radius"),
            ("own physical radius source changed", Z))
    require(own.get("source_cover") == pressure.get("source_cover"),
            ("own source cover changed", Z))
    expected_linear_rows = {
        name: pressure.get("physical_pressure_rows", {}).get(name, {}).get(
            "finite_N_absolute_error_log_bounds")
        for name in ("U", "m", "U_y")
    }
    require(own.get("actual_linear_error_input_rows") == expected_linear_rows,
            ("own linear error input rows changed", Z))
    require(own.get("actual_slow_fast_error_components") ==
            pressure.get("finite_N_axial_radial_error_components"),
            ("own slow-fast error components changed", Z))
    require(isinstance(errors, list) and len(errors) == 2 and isinstance(weighted, list) and len(weighted) == 2,
            ("own Qy rows", Z))
    require(isinstance(radius, dict), ("own physical radius", Z))
    ctx, iv = _ctx_pair()
    for error, saved in zip(errors, weighted):
        if error.get("exact_zero"):
            require(saved is None, ("zero own Qy row", Z))
            continue
        require(saved is not None, ("missing weighted own Qy row", Z))
        expected = _decode_record(error.get("log_absolute_upper"), ctx, iv) + _decode_record(radius, ctx, iv)
        require(saved.get("source_enclosure_exact_mpf_tuples") == [list(value) for value in expected._mpi_],
                ("own physical R weighted bound", Z))
        _check_record(saved, "own weighted bound")
    outgoing = raw[155]
    retained = witness.get("postrestore_full_original_155_aggregate")
    require(retained == outgoing, ("postrestore 155 aggregate changed", Z))
    require(witness.get("postrestore_upper_separate_budget_claimed") is False,
            "postrestore upper was claimed separately")
    post_log = outgoing.get("log_scale_terms")
    require(isinstance(post_log, list) and len(post_log) == 1, "postrestore aggregate shape changed")
    _check_record(post_log[0], "postrestore 155 cap")
    union = witness.get("positive_union")
    components = [row for row in weighted if row is not None] + [post_log[0]]
    union_check = _check_positive_union(components, union, "external endpoint union")
    source = witness.get("source_binding")
    require(isinstance(source, dict), "budget source binding missing")
    return {"own_reference0": True, "postrestore_155_retained": True,
            "positive_union": union_check, "no_cap_as_value": True}


def _check_packet(sample, parent_sample, reference, restore, implementation_hashes,
                  top_reference0_envelope):
    Z = sample.get("Z")
    require(Z in SAVED_Z and sample.get("source_family") == parent_sample.get("source_family"),
            ("packet identity", Z))
    raw = sample.get("raw_m15_C0_C1_terms")
    parent_raw = parent_sample.get("raw_m15_C0_C1_terms")
    old = parent_sample.get("correlated_reduced_m15_C0_C1_terms")
    reduced = sample.get("correlated_reduced_m15_C0_C1_terms")
    require(isinstance(raw, list) and len(raw) == RAW_COUNT and raw == parent_raw,
            ("raw 213 parent packet", Z))
    require(isinstance(old, list) and len(old) == PREVIOUS_REDUCED_COUNT, ("active reduced parent", Z))
    require(isinstance(reduced, list) and len(reduced) == NEW_REDUCED_COUNT,
            ("new reduced packet", Z))
    require(sample.get("original_raw_terms") == RAW_COUNT and
            sample.get("previous_reduced_terms") == PREVIOUS_REDUCED_COUNT and
            sample.get("correlated_reduced_terms") == NEW_REDUCED_COUNT and
            sample.get("unchanged_other_sectors") == UNCHANGED_COUNT, ("term counts", Z))
    ref = _find_sample(reference, Z, "reference", "1")
    rest = _find_sample(restore, Z, "restoration", "1")
    post = _find_sample(restore, Z, "postrestore", "1")
    require(raw[:124] == ref.get("m15_C0_C1_terms") and
            raw[:140] == rest.get("m15_C0_C1_terms") and
            raw[:156] == post.get("m15_C0_C1_terms"), ("exact cumulative raw prefixes", Z))
    selected = [ref.get("actual_selected_integer"), rest.get("actual_selected_integer"),
                post.get("actual_selected_integer")]
    require(selected[0] == selected[1] == selected[2], ("selected N", Z))
    pressure = _find_reference_pressure_row(reference, Z)
    require(pressure.get("selected_integer") == selected[0], ("reference pressure selected N", Z))
    for index, chart in ((123, "reference"), (139, "restoration"), (155, "postrestore")):
        require(raw[index].get("name") == chart + "_finite_N_viscous_error",
                ("raw sector identity", index, chart))
        require(len(raw[index].get("log_scale_terms", [])) == 1,
                ("raw viscous aggregate", index))
        _check_record(raw[index]["log_scale_terms"][0], ("raw cap", index))
    old_names = [term.get("name") for term in old]
    new_names = [term.get("name") for term in reduced]
    require(old_names.count(raw[123]["name"]) == 1 and old_names.count(OLD_COMBINED) == 1,
            ("old endpoint sectors", Z))
    expected_names = []
    for term in old:
        if term["name"] == raw[123]["name"]:
            continue
        expected_names.append(REPLACEMENT if term["name"] == OLD_COMBINED else term["name"])
    require(new_names == expected_names, ("reduced endpoint order", Z))
    unchanged_old = [term for term in old if term["name"] not in {raw[123]["name"], OLD_COMBINED}]
    unchanged_new = [term for term in reduced if term.get("name") != REPLACEMENT]
    require(len(unchanged_old) == UNCHANGED_COUNT and len(unchanged_new) == UNCHANGED_COUNT,
            ("unchanged sectors", Z))
    require(all(canonical(left) == canonical(right) for left, right in zip(unchanged_old, unchanged_new)),
            ("non-endpoint sectors changed", Z))
    replacement = sample.get("replacement_outer_endpoints")
    replacements = [term for term in reduced if term.get("name") == REPLACEMENT]
    require(len(replacements) == 1 and replacements[0] == replacement, ("replacement packet", Z))
    binding = replacement.get("actual_function_binding", {})
    for key, value in {
        "error_enclosure_only": True,
        "delta_definition": "R*(Qy_modified-Qy_original_leading)",
        "ordinary_Z1": True,
        "physical_R_once": True,
        "external_endpoints_not_zeroed": True,
        "signed_centers_evaluated": False,
        "no_error_cap_subtraction": True,
        "independent_own_reference0_budget": True,
        "full_original_postrestore_155_budget_retained": True,
    }.items():
        require(binding.get(key) == value, ("replacement binding", key, binding.get(key)))
    require(binding.get("function") == "D(reference native0)-D(postrestore native1=Rm)",
            ("replacement signed function", Z))
    require(sample.get("old_budget_terms_retained_as_provenance") == [raw[i] for i in (123, 139, 155)],
            ("provenance sectors", Z))
    require(sample.get("independent_target") == parent_sample.get("independent_target"),
            ("independent target changed", Z))
    source = _check_source_binding(sample.get("external_endpoint_budget_witness", {}).get("source_binding"),
                                   implementation_hashes)
    seams = _check_source_seams(sample.get("source_seam_proofs"))
    port_rows = sample.get("actual_signed_function_ports")
    ports = _check_ports(port_rows)
    reference_lower_port = next(port for port in port_rows
                                if port.get("chart") == "reference" and
                                port.get("side") == "lower")
    require(top_reference0_envelope.get("endpoint_scope") ==
            reference_lower_port.get("endpoint_scope"),
            ("reference0 graph-envelope endpoint scope", Z))
    witness = sample.get("external_endpoint_budget_witness", {}).get("own_reference0", {})
    require(witness.get("actual_reference0_graph_envelope_binding") ==
            top_reference0_envelope,
            ("reference0 graph-envelope top/witness mismatch", Z))
    chain = _check_trace_chain(sample.get("signed_trace_chain"))
    budget = _check_budget(sample, raw, reference, Z)
    caps = sample.get("largest_remaining_individual_log_caps")
    require(isinstance(caps, list) and len(caps) == 2, ("remaining caps", Z))
    for cap in caps:
        require(cap.get("cap_is_not_an_actual_residual_value") is True, ("cap as value", cap))
    contract = sample.get("contract")
    expected_contract = {
        "raw_213_verbatim": True,
        "other_197_verbatim": True,
        "independent_target_unchanged": True,
        "reserved_added_this_stage": 0,
        "ordinary_Z1_and_physical_R_once": True,
        "leading_subtraction_retained": True,
        "both_external_source_functions_retained": True,
        "no_longreshape_reference_alias": True,
        "all_other_volume_width_numerical_errors_retained": True,
        "ancestors_not_reexecuted": True,
        "source_bound_signed_trace_contraction": True,
        "error_cap_as_value": False,
    }
    require(isinstance(contract, dict), ("sample contract", Z))
    for key, value in expected_contract.items():
        require(contract.get(key) == value, ("sample contract", key, contract.get(key)))
    limits = sample.get("limits")
    require(isinstance(limits, dict), ("sample limits", Z))
    for key, value in limits.items():
        require(value is (key == "reference_restoration_internal_trace_resolved"),
                ("sample limit", key, value))
    return {
        "Z": Z,
        "raw_terms": RAW_COUNT,
        "raw_same_active_parent": True,
        "exact_raw_prefix": True,
        "selected_N_common": True,
        "correlated_reduced_terms": NEW_REDUCED_COUNT,
        "unchanged_other_sectors": UNCHANGED_COUNT,
        "source_seams": seams,
        "function_ports": ports,
        "signed_trace_chain": chain,
        "external_budget": budget,
        "source_binding": source,
        "cap_as_value": False,
        "independent_target_unchanged": True,
        "reference0_graph_envelope": {
            "passed": True,
            "top_and_witness_identical": True,
            "endpoint_scope_explicit": True,
        },
    }


def check_packet(report, native, active, reference, restore):
    require(report.get("source_family") == native.get("source_family"), "source family changed")
    require(report.get("parent_report") == NATIVE and report.get("parent_report_sha256") == NATIVE_SHA,
            "native parent identity changed")
    require(report.get("numerical_parent_report") == ACTIVE and
            report.get("numerical_parent_report_sha256") == ACTIVE_SHA,
            "active numerical parent identity changed")
    source = _check_source_binding(report.get("source_binding"), report.get("implementation_hashes", {}))
    seams = _check_source_seams(report.get("source_seam_proofs"))
    envelope_source = report.get("reference0_graph_envelope_binding")
    envelope = _check_reference0_graph_envelope(
        envelope_source, native, report.get("implementation_hashes", {}))
    samples = report.get("samples")
    require(isinstance(samples, list) and tuple(row.get("Z") for row in samples) == SAVED_Z,
            "saved Z packet set changed")
    active_by_Z = {row.get("Z"): row for row in active.get("samples", [])}
    require(set(active_by_Z) >= set(SAVED_Z), "active parent saved Z packets missing")
    packet_checks = [_check_packet(row, active_by_Z[row.get("Z")], reference, restore,
                                   report.get("implementation_hashes", {}), envelope_source)
                     for row in samples]
    top_contract = report.get("contract")
    require(top_contract == {
        "saved_packets_only": True,
        "accepted_ancestor_constructors_not_run": True,
        "cap_is_not_signed_function_value": True,
    }, "top-level contract changed")
    top_limits = report.get("limits")
    require(top_limits == {
        "numeric_Rp_closure": False,
        "global_all_five_closure": False,
        "full_goal_complete": False,
    }, "top-level open limits changed")
    return {
        "source_binding": source,
        "source_seams": seams,
        "reference0_graph_envelope": envelope,
        "samples": packet_checks,
        "top_contract": True,
        "all_unresolved_limits_open": True,
    }


def run():
    started = time.perf_counter()
    from check_original_reference_endpoint_m15 import run as run_fixture

    fixture = run_fixture()
    fixture_check = _check_fixture(fixture)
    owner = OriginalReferenceEndpointM15(precision=150)
    packet = owner.build()
    code = dict(owner.code_hashes)
    code[FIXTURE_NAME] = sha(ROOT / FIXTURE_NAME)
    code[RUNNER_NAME] = sha(ROOT / RUNNER_NAME)
    report = dict(packet, input_hashes=dict(owner.hashes), implementation_hashes=code)
    native, active, reference, restore = load_accepted_reports()
    packet_check = check_packet(report, native, active, reference, restore)
    report["checks"] = {
        "passed": True,
        "independent_fixture": fixture,
        "fixture_invariants": fixture_check,
        "packet_invariants": packet_check,
        "saved_Z": SAVED_Z,
        "raw_213_same_active_parent": True,
        "exact_raw_prefixes": True,
        "selected_N_common": True,
        "correlated_reduced_198": True,
        "unchanged_other_sectors_197": True,
        "signed_reference0_to_Rm": True,
        "both_internal_seams_proved": True,
        "reference0_and_Rm_nonzero_fixture": True,
        "independent_target_unchanged": True,
        "reserved_added_this_stage": 0,
        "cap_as_function_value": False,
        "production_is_new_owner_and_fixture_once": True,
    }
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    return report


def inspect(path: Path):
    """Inspect a saved report without constructing the owner or running fixture."""
    data = read_json(path)
    require(data.get("checks", {}).get("passed") is True, "saved endpoint report is not passed")
    check_hash_group(DIRECTORY, data.get("input_hashes"), "input")
    check_hash_group(ROOT, data.get("implementation_hashes"), "implementation")
    native, active, reference, restore = load_accepted_reports()
    packet_check = check_packet(data, native, active, reference, restore)
    fixture = data.get("checks", {}).get("independent_fixture")
    fixture_check = _check_fixture(fixture)
    require(data["checks"].get("saved_Z") == list(SAVED_Z) or
            tuple(data["checks"].get("saved_Z", [])) == SAVED_Z, "saved Z evidence changed")
    for key, expected in {
        "raw_213_same_active_parent": True,
        "exact_raw_prefixes": True,
        "selected_N_common": True,
        "correlated_reduced_198": True,
        "unchanged_other_sectors_197": True,
        "signed_reference0_to_Rm": True,
        "both_internal_seams_proved": True,
        "reference0_and_Rm_nonzero_fixture": True,
        "independent_target_unchanged": True,
        "reserved_added_this_stage": 0,
        "cap_as_function_value": False,
    }.items():
        require(data["checks"].get(key) == expected, ("saved check", key, data["checks"].get(key)))
    return {
        "packet_invariants": packet_check,
        "fixture_invariants": fixture_check,
        "checks": data["checks"],
        "execution_seconds": data.get("execution_seconds"),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-report", action="store_true")
    parser.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    if args.write_report and args.inspect_report:
        parser.error("--write-report and --inspect-report are mutually exclusive")
    path = DIRECTORY / NAME
    data = inspect(path) if args.inspect_report else run()
    if args.write_report:
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n",
                             encoding="utf-8", newline="\n")
        temporary.replace(path)
    print(json.dumps({"checks": data.get("checks", {}),
                      "execution_seconds": data.get("execution_seconds")}, indent=2), flush=True)
