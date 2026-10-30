"""Bind saved original selector primitive roles to branch-local map evidence.

This is a direct-field consumer, not a source evaluator. Raw inputs are
authenticated before decoding; saved references are opaque metadata, and
missing or protected mappings never imply source identity or a zero delta.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

PREFIX = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_"
KERNEL = PREFIX + "n1_pulse_kernel_functions.json.gz"
KERNEL_CHECK = PREFIX + "n1_pulse_kernel_functions_check.json"
COMPLETE = PREFIX + "n1_complete_future_functions.json.gz"
COMPLETE_CHECK = PREFIX + "n1_complete_future_functions_check.json"
F_ARCHIVE = PREFIX + "Rp_native_C5_pulse_source.json.gz"
F_CHECK = PREFIX + "Rp_native_C5_pulse_source_check.json"
P_REPORT = PREFIX + "native_physical_pulse_source.json"
PRIOR = "tmp/pf_selector_root_v1/result.json"
BRIDGE = "tmp/pf_selector_bridge_v1/result.json"
PINS = {
    KERNEL: "9584201b517fe3058ea6d1a169ddb7d319776eb97479b36024af4bffd4973b75",
    KERNEL_CHECK: "47343a839590c3cb08c8d726ef225415e15f4fb12a5f7cf888049ea1f20cffc0",
    COMPLETE: "b63c1cdae91fa6b083a0bc7582c28bd5f9664225f223639b28eea8060ab7185b",
    COMPLETE_CHECK: "fd7fbf5ffd6305198c9f7ed4b14f2fa2d6e660d2f414a0debd78fc140812860f",
    F_ARCHIVE: "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe",
    F_CHECK: "320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15",
    P_REPORT: "8b21f8e868025d8cc95149ffab80b54603833f5f299ec9d100602af355cea82b",
    PRIOR: "6305006836e88f6d7b638bad8e937232ea0c92f65380c199380bb5f0cbde3cf0",
    BRIDGE: "f5f7462c99cd165f814edfeece04e39e0319682b953da41762e9f3d66247ee07",
}
PROTECTED = frozenset({0, 1, 17, 33572, 25220, 263344, 263348, 263349, 262918, 205147, 263346})


class PrimitiveSourceError(RuntimeError):
    pass


def need(ok, message):
    if not ok:
        raise PrimitiveSourceError(message)


def direct(value, *path):
    for key in path:
        need(isinstance(value, dict) and key in value, "missing direct field: " + str(key))
        value = value[key]
    return value


def ref(value, label):
    need(type(value) is int, label + " is not an opaque reference")
    return value


def _map(ref_id, branch, mapping, seen, prior):
    if ref_id in PROTECTED:
        return {"status": "protected_reference_not_queried", "mapped_reference": None}
    if ref_id in prior:
        return {"status": "prior_mapping_record_reused", "record": prior[ref_id],
                "underlying_saved_status": prior[ref_id]["status"]}
    if ref_id in seen:
        return {"status": "new_mapping_metadata_reused", "record": seen[ref_id]}
    key = str(ref_id)
    if key not in mapping:
        record = {"status": "no_saved_mapping_unresolved", "mapped_reference": None}
    else:
        target = mapping[key]
        need(type(target) is int, branch + " mapping target not opaque integer")
        record = {"status": "saved_mapping_entry", "mapped_reference": target}
    seen[ref_id] = record
    return record


def verify_mapping_synthetic():
    p, prior, seen = {}, {}, {}
    missing = _map(900001, "P", p, seen, prior)
    repeated = _map(900001, "P", p, seen, prior)
    protected = _map(17, "P", {"17": 99}, seen, prior)
    protected_prior = _map(17, "P", {"17": 99}, seen, {17: {"status": "saved_mapping_entry"}})
    missing_prior = _map(900002, "P", {}, seen,
                         {900002: {"status": "no_saved_mapping_unresolved"}})
    checks = {"missing_not_identity": missing["status"] == "no_saved_mapping_unresolved",
              "repeat_not_requeried": repeated["status"] == "new_mapping_metadata_reused",
              "protected_not_queried": protected["status"] == "protected_reference_not_queried",
              "protected_takes_priority_over_prior": protected_prior["status"] == "protected_reference_not_queried",
              "prior_missing_remains_missing": missing_prior["underlying_saved_status"] == "no_saved_mapping_unresolved"}
    return {"checks": checks, "checks_passed": all(checks.values()),
            "synthetic_only": True, "source_values_computed": False}


def build_primitive_source(root, manifest):
    root = Path(root).resolve()
    declared = {item["path"]: item["sha256"] for item in manifest["files"]}
    for path, digest in PINS.items():
        need(declared.get(path) == digest, "missing canonical pin " + path)
    raw = {}
    for item in manifest["files"]:
        relative = Path(item["path"])
        need(not relative.is_absolute() and ".." not in relative.parts, "invalid manifest path")
        path = (root / relative).resolve()
        need(path.is_relative_to(root), "manifest path escapes root")
        data = path.read_bytes()
        need(hashlib.sha256(data).hexdigest() == item["sha256"], "raw input hash mismatch")
        if item["path"] in PINS:
            raw[item["path"]] = data
    kr = json.loads(raw[KERNEL_CHECK]); fr = json.loads(raw[F_CHECK]); cr = json.loads(raw[COMPLETE_CHECK])
    need(kr.get("all_passed") is True and kr.get("current_original_n1_exact_pulse_kernel_functions_installed") is True,
         "kernel receipt not accepted")
    f_flags = ("current_original_true_Rp_C5_selected_native_pulse_inlet_bridge_installed",
               "current_original_actual_Rp_native_C5_all_eight_source_functions_consumed",
               "current_original_actual_Rp_C5_frame_equals_canonical_native_smooth_functions")
    need(fr.get("all_passed") is True and all(fr.get(flag) is True for flag in f_flags),
         "F receipt not accepted")
    need(cr.get("all_passed") is True and cr.get("current_original_n1_complete_future_energy_C1_function_installed") is True and
         cr.get("current_original_n1_actual_selected_pulse_coefficient_C1_functions_installed") is True,
         "complete future receipt not accepted")
    for receipt, path in ((kr, KERNEL), (fr, F_ARCHIVE), (cr, COMPLETE)):
        need(direct(receipt, "input_hashes", Path(path).name) == PINS[path], "own archive pin mismatch")
    need(direct(fr, "input_hashes", Path(KERNEL).name) == PINS[KERNEL], "F receipt kernel pin mismatch")
    need(direct(cr, "input_hashes", Path(KERNEL).name) == PINS[KERNEL], "complete receipt kernel pin mismatch")
    kernel = json.loads(gzip.decompress(raw[KERNEL])); frame = json.loads(gzip.decompress(raw[F_ARCHIVE])); complete = json.loads(gzip.decompress(raw[COMPLETE]))
    for archive, receipt in ((kernel, kr), (frame, fr), (complete, cr)):
        family = direct(archive, "source_family")
        need(isinstance(family, dict) and bool(family) and family == direct(receipt, "source_family"),
             "source family receipt mismatch")
    need(kernel["source_family"] == frame["source_family"] == complete["source_family"],
         "saved kernel/F/complete source-family metadata differs")
    previous = json.loads(raw[PRIOR]); prior_bridge = json.loads(raw[BRIDGE])
    need(previous.get("status") == "returned" and prior_bridge.get("status") == "returned", "prior receipt not returned")
    source = direct(previous, "result", "source_roles")
    bridge = direct(prior_bridge, "result", "bridge")
    need(source.get("source_values_computed") is False and source.get("source_graph_traversed") is False and
         source.get("difference_provider_admissible_for_numeric_use") is False and
         bridge.get("P_F_function_equality_proved") is False, "prior source scope mismatch")
    p = direct(json.loads(raw[P_REPORT]), "component_result")
    need(p.get("passed") is True and direct(p, "source_prerequisite_audit").get("authenticated:" + Path(F_ARCHIVE).name) is True,
         "physical P source prerequisite missing")
    pm = direct(p, "substitution_receipt", "changed_row_mapping")
    fm = direct(frame, "Rp_native_C5_pulse_source_bindings", "rewritten_nodes")
    need(isinstance(pm, dict) and isinstance(fm, dict), "saved branch map missing")
    previous_maps = {"P": {}, "F": {}}
    def retain_prior(record, branch, role, order):
        namespace = "physical_current_native_report" if branch == "P" else "frame_Rp_native_C5_archive"
        pointer = "component_result.substitution_receipt.changed_row_mapping" if branch == "P" else "Rp_native_C5_pulse_source_bindings.rewritten_nodes"
        need(isinstance(record, dict) and record.get("branch") == branch and record.get("canonical_role") == role and
             record.get("mapping_owner") == (P_REPORT if branch == "P" else F_ARCHIVE) and
             record.get("mapping_namespace") == namespace and record.get("mapping_pointer") == pointer and
             record.get("ordinary_Z_order") == order and
             type(record.get("canonical_reference")) is int and
             ((record.get("status") == "saved_mapping_entry" and type(record.get("mapped_reference")) is int) or
              (record.get("status") == "no_saved_mapping_unresolved" and record.get("mapped_reference") is None)),
             "prior mapping metadata is malformed")
        key = record["canonical_reference"]
        if key in previous_maps[branch]:
            old = previous_maps[branch][key]
            need((old.get("status"), old.get("mapped_reference")) ==
                 (record.get("status"), record.get("mapped_reference")), "conflicting prior mappings")
        previous_maps[branch][key] = record
    for role_name, role in direct(source, "quadratic_coefficient_roles").items():
        for branch in ("P", "F"):
            for order, record in enumerate(role[branch]):
                retain_prior(record, branch, role_name, order)
    for branch in ("P", "F"):
        for order, record in enumerate(direct(source, "selected_root_roles", branch + "_saved_map")):
            retain_prior(record, branch, "ap", order)
    weights = direct(kernel, "exact_pulse_kernel_functions")
    incoming = direct(kernel, "actual_native_incoming_C1_functions")
    roles = {}
    scalar_names = ("mu", "positive_determinant_absolute", "negative_divided_determinant", "fixed_pulse_energy")
    vector_names = ("first_row", "divided_row", "incoming_factors", "scaled_main_rows", "exact_nonzero_end_energy", "original_end_energy_weights")
    for name in scalar_names:
        roles[name] = [{"canonical_ref": ref(direct(weights, name), name), "component": 0,
                        "ordinary_Z_order": 0, "canonical_Z1_unbound": True}]
    for name in vector_names:
        values = direct(weights, name)
        need(isinstance(values, list) and len(values) == 2, name + " must have two components")
        roles[name] = [{"canonical_ref": ref(value, name), "component": j,
                        "ordinary_Z_order": 0, "canonical_Z1_unbound": True}
                       for j, value in enumerate(values)]
    for name in ("m1", "m2", "energy"):
        values = direct(incoming, name)
        need(isinstance(values, list) and len(values) >= 2, name + " C1 roles missing")
        roles["incoming_" + name] = [{"canonical_ref": ref(value, name), "component": 0,
                                      "ordinary_Z_order": j} for j, value in enumerate(values[:2])]
    future = direct(complete, "function_bindings", "energy", "weighted")
    need(isinstance(future, list) and len(future) >= 2, "weighted future roles missing")
    need(future[:2] == direct(source, "weighted_future_role")[:2], "prior weighted future metadata mismatch")
    roles["weighted_future"] = [{"canonical_ref": ref(value, "weighted future"), "component": 0,
                                 "ordinary_Z_order": j} for j, value in enumerate(future[:2])]
    seen = {"P": {}, "F": {}}
    for name, rows in roles.items():
        for row in rows:
            row["canonical_owner"] = KERNEL if name != "weighted_future" else COMPLETE
            row["canonical_role"] = name
            row["reference_is_opaque_metadata"] = True
            for branch, mapping in (("P", pm), ("F", fm)):
                row[branch] = {"owner": P_REPORT if branch == "P" else F_ARCHIVE,
                    **_map(row["canonical_ref"], branch, mapping, seen[branch], previous_maps[branch])}
            row["conditional_delta"] = {"operation": "unevaluated_P_minus_F_source_role",
                                         "left_P_source_record": row["P"],
                                         "right_F_source_record": row["F"],
                                         "source_value_computed": False,
                                         "identity_and_Z_transport_proved": False}
    return {"schema": "pf_selector_primitive_source_v2", "source_only": True,
            "source_values_computed": False, "source_graph_traversed": False,
            "roles": roles, "prior_quadratic_coefficient_roles": source["quadratic_coefficient_roles"],
            "saved_kernel_F_complete_source_family_metadata_equal": True,
            "base_source_definition": {"formula": "(1-exp(-26))/4",
                "source_path": PREFIX + "n1_pulse_kernel_functions.py",
                "P_F_owner_and_Z_transport_proved": False, "source_value_computed": False},
            "determinant_convention": "selector uses positive_determinant_absolute; negative_divided_determinant is a separate source role",
            "source_function_difference_provider_bound": False,
            "numeric_use_admissible": False, "certified": False, "old_errors_removed": False,
            "unresolved": ["P/F input owners and exact rewrite semantics", "weight and mu ordinary-Z jets",
                           "source-owned correlated primitive differences", "positive determinant and root margins",
                           "capture-free substitution and certified errors"]}
