"""Run and inspect the saved long-reshape/reference endpoint contraction."""

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
from openai_ns_reconstruction.original_long_reference_m15 import (  # noqa: E402
    PARENT,
    PARENT_SHA,
    OriginalLongReferenceM15,
)
from openai_ns_reconstruction.original_long_reference_numeric_m15 import (  # noqa: E402
    LONG_Q0,
    LONG_Q1,
    LONG_REPORT,
    LONG_PREFIX_COUNT,
    LONG_VISCOUS_INDEX,
    NEW_EXTERNAL_TERM,
    NEW_REDUCED_COUNT,
    OLD_EXTERNAL_TERM,
    OLD_LONG_TERM,
    POSTRESTORE_VISCOUS_INDEX,
    PREVIOUS_REDUCED_COUNT,
    RAW_COUNT,
    REFERENCE_VISCOUS_INDEX,
    RESTORATION_VISCOUS_INDEX,
    SAVED_Z,
    UNCHANGED_COUNT,
)


NAME = PREFIX + "compliant_current_original_long_reference_m15_numeric.json"
FIXTURE_NAME = "scripts/check_original_long_reference_m15.py"
RUNNER_NAME = "scripts/run_original_long_reference_m15.py"
PARENT_NAME = PREFIX + "compliant_current_original_reference_endpoint_m15_numeric.json"
REFERENCE_NAME = PREFIX + "compliant_current_original_reference_m15_numeric.json"
RESTORE_NAME = PREFIX + "compliant_current_original_restore_m15_numeric.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def digest(value) -> str:
    return hashlib.sha256(canonical(value).encode()).hexdigest()


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


def load_reports():
    parent_path = DIRECTORY / PARENT_NAME
    require(sha(parent_path) == PARENT_SHA, "accepted endpoint parent changed")
    parent = read_json(parent_path)
    require(parent.get("checks", {}).get("passed") is True, "accepted endpoint parent is not passed")
    long = read_json(DIRECTORY / LONG_REPORT)
    reference = read_json(DIRECTORY / REFERENCE_NAME)
    restore = read_json(DIRECTORY / RESTORE_NAME)
    for label, report in (("long", long), ("reference", reference), ("restore", restore)):
        require(report.get("checks", {}).get("passed") is True, f"{label} source packet is not passed")
    return parent, long, reference, restore


def _find(rows, **criteria):
    for row in rows or ():
        if isinstance(row, dict) and all(row.get(key) == value for key, value in criteria.items()):
            return row
    raise AssertionError(("missing row", criteria))


def _check_fixture(fixture):
    require(isinstance(fixture, dict) and fixture.get("passed") is True,
            "independent long/reference fixture failed")
    expected_flags = {
        "synthetic_parameters_only": True,
        "original_source_values_evaluated": False,
        "finite_N": 17,
        "physical_E_two_jets_at_seam": True,
        "ordinary_Z1_E_and_Ey_at_seam": True,
        "genuine_p2_yZ": True,
        "unequal_native_J_divisions": True,
        "common_physical_seams": 3,
        "both_external_endpoints_nonzero": True,
    }
    for key, expected in expected_flags.items():
        require(fixture.get(key) == expected, ("fixture flag", key, fixture.get(key)))
    controls = fixture.get("controls")
    expected_controls = {
        "wrong_sign_rejected", "dropped_external_inlet_rejected",
        "missing_physical_R_rejected", "missing_leading_subtraction_rejected",
        "missing_p2_yZ_rejected", "omitted_fast_phase_term_rejected",
        "extraneous_fast_Bphi_over_N_rejected", "wrong_native_J_rejected",
        "phase_endpoint_reset_rejected", "whole_chart_Ey_equals_point_one_E_rejected",
        "caps_not_used_as_source_values",
    }
    require(isinstance(controls, dict) and set(controls) == expected_controls and
            all(value is True for value in controls.values()),
            ("fixture controls", controls))
    require(fixture.get("negative_controls") == controls, "fixture negative controls changed")
    tolerance = mp.mpf("1e-38")
    for key in ("maximum_absolute_difference", "C0_chain_difference", "Z1_chain_difference"):
        require(mp.fabs(mp.mpf(fixture.get(key))) < tolerance,
                ("fixture tolerance", key, fixture.get(key)))
    require(mp.mpf(fixture.get("maximum_p2_yZ_magnitude")) > 0,
            "fixture p2_yZ witness is empty")
    return {"passed": True, "stored_controls": len(controls), "tolerance": "1e-38"}


def _check_germ(germ):
    require(isinstance(germ, dict) and germ.get("passed") is True, "new physical germ failed")
    require(germ.get("source_only") is True and germ.get("source_values_selected") is False,
            "germ source/value contract changed")
    require(germ.get("numeric_caps_used") is False, "germ cap contract changed")
    evidence = germ.get("source_evidence", {})
    groups = evidence.get("groups", {})
    require(evidence.get("passed") is True and isinstance(groups, dict), "germ source evidence missing")
    expected_groups = {
        "long_reference_outer_formula", "flat_sigma_prime", "bridge_fields_and_radial",
        "duhamel_histories", "flat_sigma_source", "generic_root_recipe", "poisson_loop_recipe",
    }
    require(expected_groups <= set(groups), "germ AST evidence groups changed")
    require(all(groups[key].get("passed") is True for key in expected_groups), "germ AST evidence failed")
    pins = germ.get("source_pins", {})
    require(isinstance(pins, dict) and set(pins) == {"outer", "loop", "bridge"} and
            all(pin.get("pinned") is True for pin in pins.values()),
            "germ source pins changed")
    interface = germ.get("common_source_identity_rows", {})
    scope = interface.get("scope", {})
    require(interface.get("passed") is True and interface.get("source_values_selected") is False and
            isinstance(interface.get("rows"), list) and len(interface["rows"]) == 171,
            "common source interface failed")
    require(scope.get("long_reshape", {}).get("endpoint") == 1 and
            scope.get("reference", {}).get("endpoint") == 0,
            "germ endpoint scope changed")
    physical = germ.get("physical_two_jet_certificate", {})
    require(physical.get("passed") is True and physical.get("source_values_selected") is False,
            "physical two-jet certificate failed")
    require(physical.get("interval_or_cap_values_used") is False and
            physical.get("E_y") == "E/10" and physical.get("E_yy") == "E/100" and
            physical.get("V_y") == "0 at long_reshape native1 and reference native0" and
            physical.get("V_yy") == "0 at the seam",
            "physical two-jet identities changed")
    for key in ("E_y_graph_rows", "V_y_graph_rows"):
        rows = physical.get(key)
        require(isinstance(rows, list) and len(rows) == 6 and
                all(row.get("passed") is True for row in rows),
                ("physical graph rows", key))
    require(physical.get("common_Duhamel_second_jets", {}).get("passed") is True,
            "Duhamel two-jet source premise missing")
    return {
        "passed": True,
        "source_groups": len(expected_groups),
        "scope": {"long_reshape_native1": True, "reference_native0": True},
        "physical_graph_rows": len(physical["E_y_graph_rows"]) + len(physical["V_y_graph_rows"]),
    }


def _check_roots(roots):
    require(isinstance(roots, dict) and roots.get("passed") is True, "new root certificate failed")
    require(roots.get("source_only") is True and roots.get("numeric_values_selected") is False and
            roots.get("numeric_caps_used") is False, "root source/value contract changed")
    require(roots.get("row_count") == 18 and roots.get("root_rows") == 16 and roots.get("q_rows") == 2,
            "root certificate row count changed")
    rows = roots.get("rows")
    require(isinstance(rows, list) and len(rows) == 18 and all(row.get("exact_ring_row") is True for row in rows),
            "root rows are not exact")
    require({row.get("chart") for row in rows} == {"long_reshape", "reference"}, "root chart scope changed")
    require(roots.get("root_ring_passed") is True and roots.get("q_active_proof", {}).get("passed") is True,
            "root ring or active q proof failed")
    require(roots.get("q_active_proof", {}).get("active_branch_evaluated_after_root_ring") is True and
            roots["q_active_proof"].get("root_ring_passed_before_q") is True and
            roots["q_active_proof"].get("inherited_positive_eta_before_q") is True,
            "q branch ordering proof changed")
    scope = roots.get("scope", {})
    require(scope.get("long_reshape", {}).get("native_endpoint") == "native1" and
            scope.get("reference", {}).get("native_endpoint") == "native0",
            "root endpoint scope changed")
    return {"passed": True, "rows": 18, "active_q_after_ring_and_eta": True}


def _check_native_guards(guards):
    require(isinstance(guards, dict) and guards.get("passed") is True and
            guards.get("source_only") is True and guards.get("ancestor_constructors_not_run") is True,
            "native source guards missing")
    charts = guards.get("charts")
    require(isinstance(charts, dict) and set(charts) == {"long_reshape", "reference"},
            "native guard chart scope changed")
    require(all(chart.get("passed") is True for chart in charts.values()), "native chart guard failed")
    require(guards.get("no_native_J_equality_claimed") is True, "native Jacobians were equated")
    return {"passed": True, "charts": 2, "native_J_equality_claimed": False}


def _check_long0_graph_binding(source, guards, parent):
    """Check the scoped long native-0 graph certificate and its owner pins."""
    require(isinstance(source, dict) and source.get("passed") is True,
            "numeric source binding missing")
    require(source.get("physical_radius_power") == 1 and
            source.get("ordinary_Z_orders") == [0, 1] and
            source.get("ordinary_Z1_product_rule_retained") is True and
            source.get("physical_Pstar_already_in_Qy_error") is True and
            source.get("source_function_values_or_caps_not_selected") is True,
            "numeric source binding contract changed")
    require(source.get("separate_inlet_uses_own_longreshape0_pressure_packet") is True and
            source.get("longreshape0_source_packet_is_not_post_power1") is True and
            source.get("post_power1_not_reinterpreted") is True and
            source.get("postrestore_entire_original_155_cap_retained") is True and
            source.get("no_unproved_whole_chart_alias") is True and
            source.get("no_assumed_long0_post_power_seam") is True,
            "long endpoint source contract changed")
    wrapper = source.get("long_pressure_rows_owner")
    require(isinstance(wrapper, dict) and wrapper.get("passed") is True and
            wrapper.get("chart") == "long_reshape" and
            wrapper.get("derivative_chart") == "reshape" and
            wrapper.get("loop_bound_index") == "reshape" and
            wrapper.get("actual_long_source_owner") is True and
            wrapper.get("old_post_power_source_not_reinterpreted") is True,
            "actual long pressure owner binding changed")
    native_guard = source.get("long_native_physical_guard")
    require(isinstance(native_guard, dict) and native_guard.get("passed") is True and
            native_guard.get("chart") == "long_reshape" and
            native_guard.get("physical_native_guard_required") is True and
            native_guard.get("physical_RQy_factor_once") is True and
            native_guard.get("physical_Uy_Pstar_factor_once") is True and
            native_guard.get("source_guard_packet") == guards,
            "actual long native guard binding changed")

    graph = source.get("actual_long0_graph_envelope_binding")
    require(isinstance(graph, dict) and graph.get("passed") is True and
            graph.get("method") == "source_bound_linear_Qy_envelope_at_actual_long_native0",
            "actual long0 graph envelope missing")
    scope = graph.get("endpoint_scope")
    require(isinstance(scope, dict) and isinstance(scope.get("coordinate_node"), int) and
            isinstance(scope.get("endpoint_node"), int),
            "actual long0 graph endpoint scope missing")
    rows = graph.get("graph_function_rows")
    require(isinstance(rows, list) and len(rows) == 2,
            "actual long0 graph rows changed")
    for order, row in enumerate(rows):
        require(isinstance(row, dict) and row.get("ordinary_Z_order") == order and
                isinstance(row.get("Qy_function_node"), int) and
                isinstance(row.get("RQy_function_node"), int) and
                isinstance(row.get("Qy_at_long0_sha256"), str) and
                len(row["Qy_at_long0_sha256"]) == 64 and
                isinstance(row.get("RQy_at_long0_sha256"), str) and
                len(row["RQy_at_long0_sha256"]) == 64 and
                row.get("endpoint_scope") == scope and
                "Qy_at_reference0_sha256" not in row and
                "RQy_at_reference0_sha256" not in row,
                ("actual long0 graph row", order, row))
    require(graph.get("long_pressure_rows_owner") == wrapper and
            graph.get("long_native_guard") == native_guard,
            "actual long0 graph owner certificates changed")

    prior = parent.get("reference0_graph_envelope_binding", {})
    require(isinstance(prior, dict) and prior.get("passed") is True,
            "accepted reference graph certificate missing")
    require(graph.get("inherited_reference0_envelope_sha256") == digest(prior),
            "accepted reference graph certificate digest changed")
    for key in (
        "actual_archived_pressure_formula_certificate",
        "actual_archived_radial_formula_certificate",
        "exact_numeric_Qy_coefficients",
        "ordinary_Z_product_rule",
        "exact_coefficient_geometry_statements",
    ):
        require(graph.get(key) == prior.get(key), ("long0 inherited source certificate", key))
    require(isinstance(prior.get("actual_reference_slow_fast_and_Pstar_certificate"), dict),
            "accepted reference slow/fast/Pstar certificate missing")
    return {
        "passed": True,
        "scope": copy.deepcopy(scope),
        "rows": 2,
        "actual_node_ids": True,
        "native_guard_certificates_same_exact": True,
    }


def _check_source_theorem(theorem, response, roots):
    require(theorem.get("passed") is True and len(theorem.get("rows", [])) == 2,
            "source derivative theorem failed")
    require(theorem.get("actual_root_Z1_rows_consumed") == 18 and
            theorem.get("mixed_yZ_zero_rows_consumed") == 5,
            "source theorem premise counts changed")
    require(theorem.get("new_root_premise_sha256") == digest(roots), "root theorem digest changed")
    require(theorem.get("response_definition_sha256") == digest(response), "response theorem digest changed")
    return {"passed": True, "rows": 2, "root_rows": 18, "mixed_yZ_rows": 5}


def _check_trace(recovery, seam):
    require(recovery.get("passed") is True and len(recovery.get("rows", [])) == 17,
            "actual trace recovery changed")
    require(recovery.get("endpoint_scope_restricted") is True and
            recovery.get("full_raw_Uy_fingerprint") is False,
            "trace recovery scope changed")
    require(seam.get("passed") is True and seam.get("left") == "long_reshape" and
            seam.get("right") == "reference" and seam.get("physical_R_once") is True,
            "long/reference seam proof failed")
    require(seam.get("no_caps_used") is True and seam.get("global_fast_phase_not_reset") is True and
            seam.get("nonzero_external_histories_retained") is True and
            "actual_trace_recovery" in seam and "source_derivative_theorem" in seam,
            "long/reference seam contract changed")
    return {"passed": True, "recovery_rows": 17, "seam": True}


def _check_trace_chain(chain):
    require(chain == {
        "signed_coefficients": {"long_native0": 1, "Rm": -1},
        "cancelled_interfaces": ["reference_native0", "restoration_native0", "postrestore_native0"],
        "external_endpoints": {"lower": "long native0", "upper": "postrestore native1 = Rm"},
        "cancelled_interfaces_explicit": [
            "long native1 = reference native0",
            "reference native1 = restoration native0",
            "restoration native1 = postrestore native0",
        ],
        "function": "D(longreshape native0)-D(postrestore native1=Rm)",
    }, "signed trace chain changed")
    return {"signed_coefficients": {"long_native0": 1, "Rm": -1},
            "cancelled_interfaces": 3, "both_external_endpoints_retained": True}


def _check_raw_scope_probe(probe):
    """Keep the exploratory seam probe separate from the new proof."""
    require(isinstance(probe, dict), "raw long/reference scope probe missing")
    identity = probe.get("identity_rows")
    unproved = probe.get("unproved_rows")
    require(isinstance(identity, list) and len(identity) == 19,
            "raw scope probe identity row count changed")
    require(isinstance(unproved, list) and len(unproved) == 6,
            "raw scope probe unproved row count changed")
    require(probe.get("unproved_fingerprints_are_not_a_discontinuity_theorem") is True and
            probe.get("independent_of_interval_caps") is True and
            probe.get("actual_signed_centers_evaluated") is False,
            "raw scope probe contract changed")
    return {"identity_rows": 19, "unproved_rows": 6, "exploratory_only": True}


def _check_ports(ports):
    require(isinstance(ports, list) and len(ports) == 4, "long/postrestore function ports changed")
    seen = set()
    for port in ports:
        require(isinstance(port, dict), "malformed function port")
        chart, side = port.get("chart"), port.get("side")
        seen.add((chart, side))
        require(chart in {"long_reshape", "postrestore"} and side in {"lower", "upper"},
                ("function port scope", port))
        require(port.get("endpoint_scope") == {
            "coordinate_node": port.get("coordinate_node"),
            "endpoint_node": port.get("endpoint_node"),
        }, "function port endpoint scope changed")
        require(port.get("graph_scope_is_explicit") is True and
                port.get("physical_R_explicit_once_in_both_Z_rows") is True and
                port.get("ordinary_Z_orders") == [0, 1] and
                isinstance(port.get("physical_radius_node"), int) and
                isinstance(port.get("native_Jacobian_node"), int) and
                port.get("definition") == "D=R*(Qy_modified-Qy_leading)" and
                isinstance(port.get("radius_and_measure_Z_free_dependency_nodes"), list) and
                port.get("signed_coefficient") == (1 if side == "lower" else -1),
                "function port contract changed")
        require(len(port.get("modified_RQy_nodes", [])) == 2 and
                len(port.get("leading_normalized_V_nodes", [])) == 3 and
                len(port.get("leading_normalized_mean_nodes", [])) == 3 and
                len(port.get("leading_normalized_Vy_nodes", [])) == 2,
                "function port node rows changed")
    require(seen == {(chart, side) for chart in ("long_reshape", "postrestore")
                     for side in ("lower", "upper")}, "function port coverage changed")
    return {"ports": 4, "endpoint_scope_explicit": True, "physical_R_once": True}


def _record_tuples(row):
    require(isinstance(row, dict), ("record", row))
    values = row.get("source_enclosure_exact_mpf_tuples")
    if values is None and "lower_exact_mpf_tuple" in row:
        values = [row["lower_exact_mpf_tuple"], row["upper_exact_mpf_tuple"]]
    if values is None and isinstance(row.get("log_absolute_upper"), dict):
        return _record_tuples(row["log_absolute_upper"])
    require(isinstance(values, list) and len(values) == 2 and
            all(isinstance(value, list) and len(value) == 4 for value in values),
            ("record tuples", row))
    return values


def _decode(row, ctx, iv):
    return iv.mpf([ctx.make_mpf(tuple(value)) for value in _record_tuples(row)])


def _ctx_pair():
    ctx = mp.mp.clone()
    ctx.dps = 150
    iv = MPIntervalContext()
    iv.dps = 180
    return ctx, iv


def _bounds(ctx, value):
    return tuple(ctx.make_mpf(row) for row in value._mpi_)


def _check_record(row):
    require(isinstance(row, dict), "missing positive error record")
    _record_tuples(row)
    require(row.get("positive_error_bound_atom_only") is True and row.get("source_value_selected") is False,
            "positive error record was used as a source value")


def _check_positive_union(values, record, label):
    """Recompute the directed positive log union from serialized MPF tuples."""
    _check_record(record)
    ctx, iv = _ctx_pair()
    decoded = [_decode(value, ctx, iv) for value in values if value is not None]
    require(decoded, (label, "empty positive union"))
    top = max(_bounds(ctx, value)[1] for value in decoded)
    expected = iv.mpf(top) + iv.log(len(decoded))
    require(record.get("source_enclosure_exact_mpf_tuples") == [list(value) for value in expected._mpi_],
            (label, "positive union arithmetic changed"))
    return {"components": len(decoded), "exact_mpf_union": True, "cap_as_value": False}


def _check_budget(sample, long_report, source_binding):
    witness = sample.get("external_endpoint_budget_witness", {})
    own = witness.get("own_longreshape0")
    require(isinstance(own, dict) and "own_long_reference0" not in witness,
            "canonical long q=0 witness changed")
    for key, expected in {
        "physical_R_once": True,
        "physical_Pstar_already_in_error": True,
        "extra_epsilon": False,
        "source_function_values_evaluated": False,
        "post_power1_packet_not_reinterpreted": True,
        "reference0_packet_not_used_as_long_q0_cap": True,
    }.items():
        require(own.get(key) is expected, ("long q=0 witness", key, own.get(key)))
    long_q0 = _find(long_report.get("pressure_row_samples"), Z=sample["Z"], coordinate=LONG_Q0)
    require(own.get("source_report") == LONG_REPORT and own.get("source_packet_sha256") == digest(long_q0),
            "long q=0 source packet changed")
    require(own.get("coordinate") == long_q0.get("coordinate") and
            own.get("selected_integer") == long_q0.get("selected_integer") and
            own.get("source_cover") == long_q0.get("source_cover"), "long q=0 provenance changed")
    physical = long_q0.get("physical_pressure_rows", {})
    qy_errors = physical.get("Q_y", {}).get("finite_N_absolute_error_log_bounds")
    endpoint = long_q0.get("R_Q_y_endpoint", {})
    geometry = long_q0.get("geometry", {})
    require(own.get("actual_own_Qy_error_rows") == qy_errors,
            "long q=0 Q_y error rows changed")
    require(own.get("actual_own_log_radius") == geometry.get("log_physical_radius") ==
            endpoint.get("finite_N_log_physical_factor"),
            "long q=0 physical radius source changed")
    require(own.get("source_cover") == long_q0.get("source_cover"),
            "long q=0 source cover changed")
    expected_linear = {name: physical[name].get("finite_N_absolute_error_log_bounds") for name in ("U", "m", "U_y")}
    require(own.get("actual_linear_error_input_rows") == expected_linear and
            own.get("own_U_m_Uy_error_rows") == expected_linear, "long q=0 linear rows changed")
    require(own.get("actual_slow_fast_error_components") == long_q0.get("finite_N_axial_radial_error_components") and
            own.get("own_slow_fast_components") == long_q0.get("finite_N_axial_radial_error_components"),
            "long q=0 slow/fast rows changed")
    errors = own.get("actual_own_Qy_error_rows")
    weighted = own.get("radius_weighted_C0_C1_log_bounds")
    require(isinstance(errors, list) and len(errors) == 2 and
            isinstance(weighted, list) and len(weighted) == 2,
            "long q=0 Q_y rows missing")
    ctx, iv = _ctx_pair()
    for error, saved in zip(errors, weighted):
        if error.get("exact_zero"):
            require(saved is None, "zero long q=0 Q_y row was weighted")
            continue
        require(saved is not None, "missing weighted long q=0 Q_y row")
        expected = _decode(error.get("log_absolute_upper"), ctx, iv) + _decode(own["actual_own_log_radius"], ctx, iv)
        require(saved.get("source_enclosure_exact_mpf_tuples") == [list(value) for value in expected._mpi_],
                "long q=0 radius-weighted bound changed")
        _check_record(saved)
    require(own.get("actual_long0_graph_envelope_binding") ==
            source_binding.get("actual_long0_graph_envelope_binding"),
            "per-sample long0 graph binding changed")
    outgoing = sample["raw_m15_C0_C1_terms"][POSTRESTORE_VISCOUS_INDEX]
    require(witness.get("postrestore_full_original_155_aggregate") == outgoing and
            witness.get("postrestore_upper_separate_budget_claimed") is False,
            "postrestore 155 aggregate changed")
    union_inputs = witness.get("positive_union_inputs")
    union = witness.get("positive_union")
    require(isinstance(union_inputs, list) and len(union_inputs) == 3, "positive union components changed")
    require(union_inputs[:2] == own.get("radius_weighted_C0_C1_log_bounds") and
            union_inputs[2] == outgoing.get("log_scale_terms", [None])[0],
            "positive union input order changed")
    union_check = _check_positive_union(union_inputs, union, "long/reference external endpoint union")
    return {"own_long_q0": True, "full_155_retained": True, "positive_union_components": 3,
            "positive_union": union_check, "cap_as_value": False}


def _check_sample(sample, parent_sample, long_report, reference, restore, source_binding):
    Z = sample.get("Z")
    require(Z in SAVED_Z, ("sample Z", Z))
    raw = sample.get("raw_m15_C0_C1_terms")
    previous = parent_sample.get("correlated_reduced_m15_C0_C1_terms")
    reduced = sample.get("correlated_reduced_m15_C0_C1_terms")
    require(isinstance(raw, list) and len(raw) == RAW_COUNT and raw == parent_sample.get("raw_m15_C0_C1_terms"),
            ("raw packet", Z))
    require(isinstance(previous, list) and len(previous) == PREVIOUS_REDUCED_COUNT and
            isinstance(reduced, list) and len(reduced) == NEW_REDUCED_COUNT,
            ("reduced packet", Z))
    require(sample.get("original_raw_terms") == RAW_COUNT and
            sample.get("previous_reduced_terms") == PREVIOUS_REDUCED_COUNT and
            sample.get("correlated_reduced_terms") == NEW_REDUCED_COUNT and
            sample.get("unchanged_other_sectors") == UNCHANGED_COUNT,
            ("term counts", Z))
    long_q1 = _find(long_report.get("samples"), Z=Z, chart="long_reshape", coordinate=LONG_Q1)
    require(len(long_q1.get("m15_C0_C1_terms", [])) == LONG_PREFIX_COUNT and
            raw[:LONG_PREFIX_COUNT] == long_q1["m15_C0_C1_terms"], "long prefix changed")
    require(raw[LONG_VISCOUS_INDEX].get("name") == OLD_LONG_TERM and
            raw[REFERENCE_VISCOUS_INDEX].get("name") == "reference_finite_N_viscous_error" and
            raw[RESTORATION_VISCOUS_INDEX].get("name") == "restoration_finite_N_viscous_error" and
            raw[POSTRESTORE_VISCOUS_INDEX].get("name") == "postrestore_finite_N_viscous_error", "raw sectors changed")
    previous_names = [term.get("name") for term in previous]
    new_names = [term.get("name") for term in reduced]
    require(previous_names.count(OLD_LONG_TERM) == 1 and previous_names.count(OLD_EXTERNAL_TERM) == 1,
            "parent endpoint sectors changed")
    expected_names = [NEW_EXTERNAL_TERM if name == OLD_EXTERNAL_TERM else name
                      for name in previous_names if name != OLD_LONG_TERM]
    require(new_names == expected_names, "new reduced order changed")
    unchanged_old = [term for term in previous if term.get("name") not in {OLD_LONG_TERM, OLD_EXTERNAL_TERM}]
    unchanged_new = [term for term in reduced if term.get("name") != NEW_EXTERNAL_TERM]
    require(len(unchanged_old) == UNCHANGED_COUNT and unchanged_new == unchanged_old,
            "other 196 sectors changed")
    selected = sample.get("selected_N_evidence", {})
    ref0 = _find(reference.get("samples"), Z=Z, chart="reference",
                 coordinate={"kind": "exact_decimal_reference", "value": "0"})
    post1 = _find(restore.get("samples"), Z=Z, chart="postrestore",
                  coordinate={"kind": "exact_decimal_postrestore", "value": "1"})
    require(selected.get("long_reshape_native1") == long_q1.get("actual_selected_integer") and
            selected.get("reference_native0") == ref0.get("actual_selected_integer") and
            selected.get("postrestore_native1") == post1.get("actual_selected_integer") and
            selected.get("same_selected_N") is True and
            selected["long_reshape_native1"] == selected["reference_native0"] == selected["postrestore_native1"],
            "selected N evidence changed")
    budget = _check_budget(sample, long_report, source_binding)
    ports = _check_ports(sample.get("actual_signed_function_ports"))
    chain = _check_trace_chain(sample.get("signed_trace_chain"))
    require(sample.get("external_endpoint_budget_witness", {}).get("source_binding") == source_binding,
            "per-sample source binding changed")
    raw_provenance = [copy.deepcopy(raw[index]) for index in (
        LONG_VISCOUS_INDEX, REFERENCE_VISCOUS_INDEX, RESTORATION_VISCOUS_INDEX,
        POSTRESTORE_VISCOUS_INDEX,
    )]
    require(sample.get("old_budget_terms_retained_as_provenance") == raw_provenance,
            "old viscous provenance changed")
    require(sample.get("previous_external_replacement_retained_as_provenance") ==
            next(term for term in previous if term.get("name") == OLD_EXTERNAL_TERM),
            "previous external provenance changed")
    require(sample.get("independent_target") == parent_sample.get("independent_target"),
            "independent target changed")
    require(sample.get("contract", {}).get("raw_213_verbatim") is True and
            sample["contract"].get("other_196_verbatim") is True and
            sample["contract"].get("long_q0_uses_own_pressure_packet") is True and
            sample["contract"].get("post_power1_not_reinterpreted") is True and
            sample["contract"].get("independent_target_unchanged") is True and
            sample["contract"].get("error_cap_as_value") is False,
            "sample contract changed")
    require(sample.get("limits", {}).get("long_reference_internal_trace_resolved") is True and
            all(value is False for key, value in sample["limits"].items()
                if key != "long_reference_internal_trace_resolved"), "sample limits changed")
    return {"Z": Z, "raw_213": True, "prefix_106": True, "reduced_197": True,
            "unchanged_196": True, "selected_N": True, "budget": budget, "ports": ports,
            "trace_chain": chain, "cap_as_value": False}


def check_packet(report, parent, long_report, reference, restore):
    require(report.get("source_family") == parent.get("source_family"), "source family changed")
    require(report.get("parent_report") == PARENT_NAME and report.get("parent_report_sha256") == PARENT_SHA,
            "parent identity changed")
    raw_probe_check = _check_raw_scope_probe(report.get("raw_scope_probe"))
    germ_check = _check_germ(report.get("new_physical_germ"))
    roots_check = _check_roots(report.get("new_actual_root_certificate"))
    guards_check = _check_native_guards(report.get("native_source_guards"))
    response = report.get("common_physical_response", {})
    require(response.get("endpoint_scope") == "long_reshape native1/reference native0; all ordinary Z functions",
            "physical response scope changed")
    require(report.get("source_derivative_theorem", {}).get("original_p2_y_not_zeroed") is True and
            response.get("true_Z2_field_ports_consumed_by_pyz") is True and
            response.get("no_actual_phase_or_source_values_selected") is True and
            response.get("new_germ_certificate_sha256") == digest(report.get("new_physical_germ")),
            "physical response source contract changed")
    theorem_check = _check_source_theorem(report.get("source_derivative_theorem", {}), response,
                                           report.get("new_actual_root_certificate", {}))
    trace_check = _check_trace(report.get("actual_trace_recovery", {}), report.get("source_seam_proof", {}))
    source = report.get("numerical_source_binding", {})
    source_graph_check = _check_long0_graph_binding(source, report.get("native_source_guards", {}), parent)
    samples = report.get("samples")
    require(isinstance(samples, list) and tuple(row.get("Z") for row in samples) == SAVED_Z,
            "saved Z set changed")
    parent_by_Z = {row.get("Z"): row for row in parent.get("samples", [])}
    source_binding = copy.deepcopy(dict(source))
    sample_checks = [_check_sample(row, parent_by_Z[row.get("Z")], long_report, reference, restore, source_binding)
                     for row in samples]
    require(report.get("contract") == {
        "accepted_ancestor_constructors_not_run": True,
        "source_functions_not_cap_values": True,
        "raw_213_and_independent_target_retained": True,
    }, "top contract changed")
    limits = report.get("limits", {})
    require(limits.get("original_signed_finite_N_endpoint_values") is False and
            limits.get("full_goal_complete") is False and
            all(value is False for key, value in limits.items()
                if key not in {"original_signed_finite_N_endpoint_values"}), "top limits changed")
    return {"raw_scope_probe": raw_probe_check, "germ": germ_check, "roots": roots_check,
            "native_guards": guards_check, "long0_graph_binding": source_graph_check,
            "source_theorem": theorem_check, "trace": trace_check, "samples": sample_checks,
            "all_open_limits": True, "cap_as_value": False}


def run():
    started = time.perf_counter()
    from check_original_long_reference_m15 import run_fixture

    fixture = run_fixture()
    fixture_check = _check_fixture(fixture)
    owner = OriginalLongReferenceM15(precision=150)
    packet = owner.build()
    code = dict(owner.code_hashes)
    code[FIXTURE_NAME] = sha(ROOT / FIXTURE_NAME)
    code[RUNNER_NAME] = sha(ROOT / RUNNER_NAME)
    report = dict(packet, input_hashes=dict(owner.hashes), implementation_hashes=code)
    parent, long_report, reference, restore = load_reports()
    packet_check = check_packet(report, parent, long_report, reference, restore)
    report["checks"] = {
        "passed": True,
        "independent_fixture": fixture,
        "fixture_invariants": fixture_check,
        "packet_invariants": packet_check,
        "saved_Z": list(SAVED_Z),
        "raw_213_same_endpoint_parent": True,
        "exact_long_prefix_106": True,
        "selected_N_common": True,
        "correlated_reduced_197": True,
        "unchanged_other_sectors_196": True,
        "independent_long0_and_full155_budget": True,
        "source_functions_not_caps": True,
        "production_is_new_owner_and_fixture_once": True,
    }
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    return report


def inspect(path: Path):
    """Inspect hashes and serialized packet shape without owner or fixture."""
    data = read_json(path)
    require(data.get("checks", {}).get("passed") is True, "saved long/reference report is not passed")
    check_hash_group(DIRECTORY, data.get("input_hashes"), "input")
    check_hash_group(ROOT, data.get("implementation_hashes"), "implementation")
    parent, long_report, reference, restore = load_reports()
    packet_check = check_packet(data, parent, long_report, reference, restore)
    stored_fixture = data.get("checks", {}).get("independent_fixture")
    fixture_check = _check_fixture(stored_fixture)
    for key, expected in {
        "saved_Z": list(SAVED_Z),
        "raw_213_same_endpoint_parent": True,
        "exact_long_prefix_106": True,
        "selected_N_common": True,
        "correlated_reduced_197": True,
        "unchanged_other_sectors_196": True,
        "independent_long0_and_full155_budget": True,
        "source_functions_not_caps": True,
    }.items():
        require(data.get("checks", {}).get(key) == expected, ("saved check", key))
    return {"checks": data["checks"], "packet_invariants": packet_check,
            "fixture_invariants": fixture_check, "execution_seconds": data.get("execution_seconds")}


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
