"""Run and inspect the source-contracted original micro viscous m15 state.

Production consumes the accepted repaired-Rp packet and the source-bound
micro archives through ``OriginalMicroViscousM15``.  It records the raw 213
term view, the 200-term correlated view, and the independent finite-N
physics fixture.  Inspection only verifies hashes and serialized packet
invariants; it never constructs the production owner or reruns calculations.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DIRECTORY = ROOT / "experiments" / "root_st073"

from openai_ns_reconstruction.original_core_numeric import PREFIX
from openai_ns_reconstruction.original_micro_viscous_m15 import (
    EndpointFingerprint,
    OriginalMicroViscousM15,
)


NAME = PREFIX + "compliant_current_original_micro_viscous_m15_numeric.json"
REPAIRED_NAME = PREFIX + "compliant_current_original_repaired_Rp_m15_numeric.json"
FIXTURE_NAME = "scripts/check_original_micro_viscous_m15.py"
RUNNER_NAME = "scripts/run_original_micro_viscous_m15.py"
SAVED_Z = ("-.3", "0", ".3")
MISSING_Z = "0.1"
FIRST_VISCOUS = "first_micro_finite_N_viscous_error"
SECOND_VISCOUS = "second_micro_finite_N_viscous_error"
REPLACEMENT = "micro_shared_finite_N_viscous_outer_endpoint_error"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value):
    """Stable JSON bytes for asserting unchanged serialized term content."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def check_source_identity(trace_identity):
    if not trace_identity.get("passed"):
        raise AssertionError("source identity proof did not pass")
    if trace_identity.get("common_ordinary_Z_independent_radius") is not True:
        raise AssertionError("ordinary Z/radius source identity missing")
    if trace_identity.get("upper_caps_not_used_as_function_values") is not True:
        raise AssertionError("upper caps were treated as source functions")
    rows = trace_identity.get("source_identity_rows")
    if not isinstance(rows, list) or len(rows) != 14:
        raise AssertionError("fourteen direct U/m/Q/RQy source identity rows required")
    keys = set()
    for row in rows:
        if set(row) != {"field", "physical_y_order", "ordinary_Z_order", "left_node", "right_node", "canonical_endpoint_sha256"}:
            raise AssertionError(("source identity row shape", row))
        if row["field"] not in {"U", "m", "Q", "R_Q_y"}:
            raise AssertionError(("source identity field", row))
        if row["physical_y_order"] not in {0, 1} or row["ordinary_Z_order"] not in {0, 1}:
            raise AssertionError(("source identity derivative orders", row))
        if row["field"] == "R_Q_y" and row["physical_y_order"] != 0:
            raise AssertionError(("physical RQy derivative order", row))
        if not isinstance(row["left_node"], (int, str)) or isinstance(row["left_node"], bool):
            raise AssertionError(("source identity node IDs", row))
        if not isinstance(row["right_node"], (int, str)) or isinstance(row["right_node"], bool):
            raise AssertionError(("source identity node IDs", row))
        digest = row["canonical_endpoint_sha256"]
        if not isinstance(digest, str) or len(digest) != 64:
            raise AssertionError(("source identity digest", row))
        keys.add((row["field"], row["physical_y_order"], row["ordinary_Z_order"]))
    expected_keys = {
        (field, y_order, z_order)
        for field in ("U", "m", "Q")
        for y_order in (0, 1)
        for z_order in (0, 1)
    }
    expected_keys.update({("R_Q_y", 0, 0), ("R_Q_y", 0, 1)})
    if keys != expected_keys:
        raise AssertionError("source identity rows are not a complete U/m/Q/RQy derivative grid")
    lower_flat = trace_identity.get("original_lower_flat_source_zero")
    if not isinstance(lower_flat, dict) or set(lower_flat.get("saved_Z", [])) != set(SAVED_Z):
        raise AssertionError(("lower flat saved Z proof", lower_flat))
    for key in (
        "exact_sc_half_packets_verified",
        "all_five_relative_incoming_histories_zero",
        "source_producer_checks_actual_U_Q_Qy_errors",
    ):
        if lower_flat.get(key) is not True:
            raise AssertionError(("lower flat source proof", key, lower_flat))
    return dict(
        rows=len(rows),
        fields=["U", "m", "Q", "R_Q_y"],
        physical_y_orders=[0, 1],
        ordinary_Z_orders=[0, 1],
        direct_Qy_source_IDs=True,
        exact_source_fingerprints=True,
        lower_flat_saved_sc_half_packets=True,
        lower_flat_all_five_histories_zero=True,
        lower_flat_source_producer_U_Q_Qy_errorzero=True,
    )


def adversarial_endpoint_fingerprint_fixture():
    """Exercise exact empty integrals, scope capture, opaque leaves and 0/0."""
    # An empty exact integral must reduce to zero before its body is visited.
    empty_nodes = [
        {"operation": "bound_variable", "name": "t"},
        {"operation": "exact_rational", "numerator": 1, "denominator": 1},
        {"operation": "definite_integral", "variable": "t", "lower": 1, "upper": 1, "integrand": 0},
        {"operation": "exact_rational", "numerator": 1, "denominator": 1},
    ]
    empty = EndpointFingerprint(empty_nodes, {"coordinate": 3})
    empty_value = empty.at(2)
    if empty_value[1] != 0 or 2 not in empty.empty:
        raise AssertionError(("exact empty integral", empty_value, empty.empty))

    # Substitution of the outer bound variable must not capture a distinct
    # inner variable with a different node identity.
    capture_nodes = [
        {"operation": "bound_variable", "name": "outer"},
        {"operation": "bound_variable", "name": "inner"},
        {"operation": "exact_rational", "numerator": 2, "denominator": 1},
        {"operation": "function_substitution", "variable": 0, "value": 2, "expression": 1, "Z_independent_substitution": True},
        {"operation": "function_substitution", "variable": 0, "value": 2, "expression": 0, "Z_independent_substitution": True},
    ]
    capture = EndpointFingerprint(capture_nodes, {"coordinate": 5})
    inner_direct = capture.at(1)
    inner_after_outer_substitution = capture.at(3)
    outer_after_substitution = capture.at(4)
    if inner_after_outer_substitution != inner_direct or outer_after_substitution[1] != 2:
        raise AssertionError(("bound variable capture", inner_direct, inner_after_outer_substitution, outer_after_substitution))

    # An opaque source function must retain the substituted coordinate in its
    # fingerprint, so f(s=0) and f(s=1) cannot alias.
    opaque_nodes = [
        {"operation": "bound_variable", "name": "Z_variable"},
        {"operation": "bound_variable", "name": "coordinate"},
        {"operation": "exact_rational", "numerator": 0, "denominator": 1},
        {"operation": "exact_rational", "numerator": 1, "denominator": 1},
        {"operation": "current_original_leading_function_recipe", "coordinate": 1, "Z_variable": 0},
        {"operation": "function_substitution", "variable": 1, "value": 2, "expression": 4, "Z_independent_substitution": True},
        {"operation": "function_substitution", "variable": 1, "value": 3, "expression": 4, "Z_independent_substitution": True},
    ]
    opaque = EndpointFingerprint(opaque_nodes, {"coordinate": 1})
    opaque_zero = opaque.at(5)
    opaque_one = opaque.at(6)
    if opaque_zero == opaque_one or opaque_zero[1] is not None or opaque_one[1] is not None:
        raise AssertionError(("opaque function substitution alias", opaque_zero, opaque_one))

    quotient_nodes = [
        {"operation": "exact_rational", "numerator": 1, "denominator": 1},
        {"operation": "exact_rational", "numerator": 0, "denominator": 1},
        {"operation": "positive_quotient", "numerator": 0, "denominator": 1},
    ]
    quotient = EndpointFingerprint(quotient_nodes, {"coordinate": 3})
    try:
        quotient.at(2)
    except ValueError:
        quotient_rejected = True
    else:
        quotient_rejected = False
    if not quotient_rejected:
        raise AssertionError("undefined endpoint quotient was admitted")
    return dict(
        passed=True,
        exact_empty_integral=True,
        bound_variable_capture_blocked=True,
        opaque_function_substitution_distinguished=True,
        undefined_quotient_rejected=True,
    )


def check_contract(packet):
    trace = packet.get("signed_trace_chain")
    if not isinstance(trace, dict):
        raise AssertionError("signed trace chain missing")
    if trace.get("signed_coefficients") != {"rminus_exact_zero": 1, "second_s2": -1}:
        raise AssertionError(("signed trace coefficients", trace))
    if trace.get("cancelled_interfaces") != ["first_s1"]:
        raise AssertionError(("cancelled interface", trace))
    contract = packet.get("contract")
    required_contract = {
        "raw_213_verbatim": True,
        "proved_physical_trace_contracted": True,
        "ordinary_Z1_and_radius_preserved": True,
        "only_first_second_viscous_representation_changed": True,
        "second_positive_aggregate_budget_retained_for_upper": True,
        "all_volume_width_and_source_errors_retained": True,
        "no_interval_cap_subtraction": True,
        "ancestors_not_reexecuted": True,
    }
    if not isinstance(contract, dict) or any(contract.get(k) is not v for k, v in required_contract.items()):
        raise AssertionError(("contract flags", contract))
    limits = packet.get("limits")
    if not isinstance(limits, dict) or limits.get("first_second_finite_N_viscous_internal_trace_resolved") is not True:
        raise AssertionError(("resolved first/second viscous limit", limits))
    for key, value in limits.items():
        if key != "first_second_finite_N_viscous_internal_trace_resolved" and value is not False:
            raise AssertionError(("open limit overclaimed", key, value))
    return dict(
        signed_coefficients=copy.deepcopy(trace["signed_coefficients"]),
        cancelled_interfaces=copy.deepcopy(trace["cancelled_interfaces"]),
        contract_flags=copy.deepcopy(required_contract),
        first_second_internal_trace_resolved=True,
        all_other_limits_open=True,
    )


def check_positive_upper_budget_witness(packet):
    """Check source AST evidence and the two serialized positive caps."""
    witness = packet.get("positive_upper_budget_witness")
    if not isinstance(witness, dict):
        raise AssertionError("positive upper-budget witness missing")
    source = witness.get("source")
    if not isinstance(source, dict) or source.get("passed") is not True:
        raise AssertionError(("positive budget source proof", source))
    if source.get("cap_subtraction_used") is not False:
        raise AssertionError(("positive budget cap subtraction", source))
    if source.get("removed_first_chart_cap_need_not_be_dominated") is not True:
        raise AssertionError(("cancelled first cap obligation", source))
    expected_statements = {
        "OriginalFirstMicroM15.finite_N_errors": [
            "viscous=[None if v is None else v+lr for v in errors['Q_y']]",
        ],
        "OriginalSecondMicroM15.evaluate": [
            "upper=self.finite_N_errors(Z,s,packet)",
            "lower=self.finite_N_errors(Z,iv.mpf(1),lower_packet)",
            "error[name]=[self.rows.sum_logs((upper[name][k],lower[name][k])) for k in (0,1)]",
            "bound=self.rows.sum_logs(error[name])",
            "top=iv.mpf(_bounds(self.audit,bound)[1])",
        ],
        "OriginalFirstMicroPressureRows.sum_logs": [
            "values=[v for v in values if v is not None]",
            "high=max(_bounds(self.ctx,v)[1] for v in values)",
            "return self.iv.mpf(high)+self.iv.log(len(values))",
        ],
        "OriginalFirstMicroM15.evaluate": [
            "lower=self.rows.evaluate(Z,'.5',selected_sc_multiple=True)",
            "if not all(v['exact_zero'] for name in ('U','Q','Q_y') for v in lower['physical_pressure_rows'][name]['finite_N_absolute_error_log_bounds']): raise ValueError('Actual flat lower finite-N pressure endpoint must be zero')",
        ],
    }
    source_rows = source.get("source")
    if not isinstance(source_rows, list) or len(source_rows) != 4:
        raise AssertionError(("positive budget source AST rows", source_rows))
    seen_methods = set()
    for row in source_rows:
        method = row.get("method")
        if method not in expected_statements or row.get("exact_source_statements") != expected_statements[method]:
            raise AssertionError(("positive budget source AST", row))
        digest = row.get("source_sha256")
        if not isinstance(digest, str) or len(digest) != 64:
            raise AssertionError(("positive budget source pin", row))
        seen_methods.add(method)
    if seen_methods != set(expected_statements):
        raise AssertionError(("positive budget source methods", seen_methods))
    rows = witness.get("serialized_C0_C1")
    if not isinstance(rows, list) or len(rows) != 2:
        raise AssertionError(("positive budget serialized rows", rows))
    for row in rows:
        if row.get("ordinary_Z_order") not in {0, 1}:
            raise AssertionError(("positive budget Z order", row))
        combined = row.get("combined_endpoint_log_bound")
        if not isinstance(combined, dict) or combined.get("exact_zero") is not False:
            raise AssertionError(("combined endpoint bound", row))
        log_bound = combined.get("log_absolute_upper")
        if not isinstance(log_bound, dict) or any(
            not isinstance(log_bound.get(key), list) or len(log_bound[key]) != 4
            for key in ("lower_exact_mpf_tuple", "upper_exact_mpf_tuple")
        ):
            raise AssertionError(("combined endpoint exact tuple", row))
        retained = row.get("retained_common_log_cap")
        if not isinstance(retained, dict):
            raise AssertionError(("retained common cap", row))
        if retained.get("positive_error_bound_atom_only") is not True or retained.get("source_value_selected") is not False:
            raise AssertionError(("retained cap interpretation", row))
        if row.get("common_cap_dominates_combined_budget") is not True:
            raise AssertionError(("combined cap domination", row))
    if witness.get("interpretation") != "the combined second bound covers D2; no domination of cancelled D1 cap is required":
        raise AssertionError(("positive budget interpretation", witness))
    return dict(
        source_rows=len(source_rows),
        source_ast_passed=True,
        serialized_C0_C1_rows=len(rows),
        combined_caps_dominated=True,
        cap_subtraction_used=False,
        removed_first_cap_need_not_be_dominated=True,
    )


def check_packet(packet, parent_packet, Z):
    if packet.get("Z") != Z or packet.get("source_family") != parent_packet.get("source_family"):
        raise AssertionError(("packet identity", Z))
    raw = packet.get("raw_m15_C0_C1_terms")
    parent_raw = parent_packet.get("m15_C0_C1_terms")
    if not isinstance(raw, list) or len(raw) != 213 or raw != parent_raw:
        raise AssertionError(("raw 213 terms differ from repaired parent", Z))
    if packet.get("original_raw_terms") != 213 or packet.get("previous_reduced_terms") != 201:
        raise AssertionError(("raw/previous reduced counts", Z))
    reduced = packet.get("correlated_reduced_m15_C0_C1_terms")
    old_reduced = parent_packet.get("reduced_m15_C0_C1_terms")
    if not isinstance(reduced, list) or len(reduced) != 200 or len(old_reduced) != 201:
        raise AssertionError(("correlated reduced count", Z, len(reduced), len(old_reduced)))
    if packet.get("correlated_reduced_terms") != 200:
        raise AssertionError(("serialized correlated count", Z))

    old_names = [term["name"] for term in old_reduced]
    new_names = [term["name"] for term in reduced]
    if old_names.count(FIRST_VISCOUS) != 1 or old_names.count(SECOND_VISCOUS) != 1:
        raise AssertionError(("parent viscous sectors", Z))
    if FIRST_VISCOUS in new_names or SECOND_VISCOUS in new_names or new_names.count(REPLACEMENT) != 1:
        raise AssertionError(("correlated viscous sector names", Z, new_names))
    expected_names = [
        REPLACEMENT if term["name"] == SECOND_VISCOUS else term["name"]
        for term in old_reduced
        if term["name"] != FIRST_VISCOUS
    ]
    if new_names != expected_names:
        raise AssertionError(("correlated term order", Z))

    unchanged_old = [term for term in old_reduced if term["name"] not in {FIRST_VISCOUS, SECOND_VISCOUS}]
    unchanged_new = [term for term in reduced if term["name"] != REPLACEMENT]
    if len(unchanged_old) != 199 or len(unchanged_new) != 199:
        raise AssertionError(("unchanged term count", Z))
    if [term["name"] for term in unchanged_old] != [term["name"] for term in unchanged_new]:
        raise AssertionError(("unchanged term names", Z))
    if any(canonical(left) != canonical(right) for left, right in zip(unchanged_old, unchanged_new)):
        raise AssertionError(("non-viscous term bytes changed", Z))
    retained_error_names = [
        term["name"]
        for term in unchanged_old
        if any(token in term["name"] for token in ("volume", "width", "source_error"))
    ]
    if any(name not in [term["name"] for term in unchanged_new] for name in retained_error_names):
        raise AssertionError(("volume/width/source errors not retained", Z))

    old_budget = packet.get("old_budget_terms_retained_as_provenance")
    if old_budget != [raw[12], raw[16]]:
        raise AssertionError(("old viscous budget provenance", Z))
    if [term["name"] for term in old_budget] != [FIRST_VISCOUS, SECOND_VISCOUS]:
        raise AssertionError(("old viscous budget names", Z))
    replacement = packet.get("replacement_outer_endpoint")
    replacement_rows = [term for term in reduced if term["name"] == REPLACEMENT]
    if len(replacement_rows) != 1 or replacement_rows[0] != replacement:
        raise AssertionError(("replacement endpoint serialization", Z))
    binding = replacement.get("actual_function_binding", {})
    for key, value in {
        "error_enclosure_only": True,
        "ordinary_Z1": True,
        "bound_is_conservative_not_actual_signed_value": True,
        "no_other_error_caps_removed": True,
    }.items():
        if binding.get(key) is not value:
            raise AssertionError(("replacement binding", Z, key, binding))
    if binding.get("function") != "-delta(RQy)(second_micro s=2)":
        raise AssertionError(("replacement physical endpoint", Z))

    source_identity = check_source_identity(packet.get("trace_identity", {}))
    positive_budget = check_positive_upper_budget_witness(packet)
    if packet.get("independent_target") != parent_packet.get("independent_target"):
        raise AssertionError(("independent target changed", Z))
    contract = check_contract(packet)
    caps = packet.get("largest_remaining_individual_log_caps")
    if not isinstance(caps, list) or len(caps) != 2:
        raise AssertionError(("remaining cap ledger", Z))
    for cap in caps:
        if cap.get("ordinary_Z_order") not in {0, 1}:
            raise AssertionError(("remaining cap Z order", Z, cap))
        if cap.get("name") != "restoration_finite_N_viscous_error":
            raise AssertionError(("remaining cap source", Z, cap))
        if not isinstance(cap.get("normalized_log_cap_upper_approximation"), str):
            raise AssertionError(("remaining cap approximation", Z, cap))
        exact_tuple = cap.get("normalized_log_cap_upper_exact_mpf_tuple")
        if not isinstance(exact_tuple, list) or len(exact_tuple) != 4:
            raise AssertionError(("remaining cap exact tuple", Z, cap))
        if cap.get("cap_is_not_an_actual_residual_value") is not True:
            raise AssertionError(("remaining cap interpretation", Z, cap))
    return dict(
        Z=Z,
        raw_terms=len(raw),
        raw_same_parent=True,
        previous_reduced_terms=len(old_reduced),
        correlated_reduced_terms=len(reduced),
        dropped_first_viscous=True,
        renamed_second_viscous=True,
        unchanged_view_terms=199,
        retained_volume_width_source_error_terms=len(retained_error_names),
        old_budget_indices=[12, 16],
        replacement_name=REPLACEMENT,
        source_identity=source_identity,
        positive_upper_budget_witness=positive_budget,
        independent_target_unchanged=True,
        contract=contract,
        remaining_cap_rows=len(caps),
        remaining_cap_exact_mpf_tuple_and_approximation=True,
    )


def run():
    started = time.perf_counter()
    from check_original_micro_viscous_m15 import run as run_fixture

    fixture = run_fixture()
    if not fixture["passed"]:
        raise AssertionError("independent micro viscous fixture failed")
    fingerprint_checks = adversarial_endpoint_fingerprint_fixture()
    owner = OriginalMicroViscousM15(precision=100)
    parent_by_Z = {row["Z"]: row for row in owner.parent["samples"]}
    packets = [owner.evaluate(Z) for Z in SAVED_Z]
    packet_checks = [check_packet(packet, parent_by_Z[Z], Z) for packet, Z in zip(packets, SAVED_Z)]
    missing_rejected = 0
    try:
        owner.evaluate(MISSING_Z)
    except (ValueError, ArithmeticError):
        missing_rejected = 1
    else:
        raise AssertionError(("missing saved Z admitted", MISSING_Z))

    hashes = dict(owner.hashes)
    code = dict(owner.code_hashes)
    code[FIXTURE_NAME] = sha(ROOT / FIXTURE_NAME)
    code[RUNNER_NAME] = sha(ROOT / RUNNER_NAME)
    return dict(
        source_family=owner.family,
        samples=packets,
        input_hashes=hashes,
        implementation_hashes=code,
        checks=dict(
            passed=True,
            independent_fixture=fixture,
            saved_Z=SAVED_Z,
            actual_micro_viscous_packets=len(packets),
            raw_213_same_repaired_parent=True,
            old_reduced_201_view_retained=True,
            correlated_reduced_200=True,
            dropped_first_micro_viscous_index=12,
            renamed_second_micro_viscous_index=16,
            unchanged_nonviscous_view_terms=199,
            all_volume_width_source_errors_retained=True,
            exact_trace_coefficients={"rminus_exact_zero": 1, "second_s2": -1},
            exact_cancelled_interface="first_s1",
            source_identity_rows_per_packet=14,
            direct_Qy_and_Z1_source_ids=True,
            lower_actual_flat_proof_checks=True,
            positive_upper_budget_witness_rows=2,
            positive_upper_budget_source_ast_passed=True,
            positive_upper_budget_cap_subtraction_used=False,
            removed_first_cap_need_not_be_dominated=True,
            endpoint_fingerprint_adversarial=fingerprint_checks,
            remaining_cap_exact_mpf_tuple_and_approximation=True,
            independent_target_unchanged=True,
            missing_Z_rejected=missing_rejected,
            production_is_new_calculation_only=True,
        ),
        packet_checks=packet_checks,
        execution_seconds=round(time.perf_counter() - started, 3),
        limits=dict(
            first_second_finite_N_viscous_internal_trace_resolved=True,
            remaining_outer_trace_signed_value_resolved=False,
            other_internal_traces_resolved=False,
            independent_numeric_Rp_seam=False,
            global_five_moment_numeric_join=False,
            exact_heat_exterior=False,
            stress_cone=False,
            flat_remainder=False,
            genuine_temporal_n_ge_2=False,
            full_corrected_NS=False,
        ),
    )


def inspect(path):
    """Verify report hashes and packet invariants without recalculation."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if not data["checks"]["passed"] or data["source_family"] is None:
        raise AssertionError("saved micro viscous report is not passed")
    for base, group in ((DIRECTORY, data["input_hashes"]), (ROOT, data["implementation_hashes"])):
        for name, digest in group.items():
            if sha(base / name) != digest:
                raise AssertionError(("hash changed", name))
    parent_path = DIRECTORY / REPAIRED_NAME
    parent = json.loads(parent_path.read_text(encoding="utf-8"))
    parent_by_Z = {row["Z"]: row for row in parent["samples"]}
    samples = data["samples"]
    if len(samples) != 3 or tuple(sample["Z"] for sample in samples) != SAVED_Z:
        raise AssertionError("saved Z packet set changed")
    for sample in samples:
        check_packet(sample, parent_by_Z[sample["Z"]], sample["Z"])
    adversarial = data["checks"].get("endpoint_fingerprint_adversarial")
    if not isinstance(adversarial, dict) or adversarial.get("passed") is not True:
        raise AssertionError("endpoint fingerprint adversarial evidence missing")
    for key in (
        "exact_empty_integral",
        "bound_variable_capture_blocked",
        "opaque_function_substitution_distinguished",
        "undefined_quotient_rejected",
    ):
        if adversarial.get(key) is not True:
            raise AssertionError(("endpoint fingerprint evidence", key))
    if data["checks"].get("remaining_cap_exact_mpf_tuple_and_approximation") is not True:
        raise AssertionError("remaining cap exact tuple/approximation evidence missing")
    for key, expected in {
        "lower_actual_flat_proof_checks": True,
        "positive_upper_budget_witness_rows": 2,
        "positive_upper_budget_source_ast_passed": True,
        "positive_upper_budget_cap_subtraction_used": False,
        "removed_first_cap_need_not_be_dominated": True,
    }.items():
        if data["checks"].get(key) != expected:
            raise AssertionError(("positive budget/lower flat evidence", key, data["checks"].get(key)))
    if data["checks"]["missing_Z_rejected"] != 1:
        raise AssertionError("missing-Z rejection evidence changed")
    if not data["checks"]["production_is_new_calculation_only"]:
        raise AssertionError("production scope overclaimed")
    if data["limits"].get("first_second_finite_N_viscous_internal_trace_resolved") is not True:
        raise AssertionError("resolved viscous trace limit missing")
    for key, value in data["limits"].items():
        if key != "first_second_finite_N_viscous_internal_trace_resolved" and value is not False:
            raise AssertionError(("saved report overclaims limit", key, value))
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
