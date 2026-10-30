"""Bind actual selector-equation roles to saved P/F rewrite entries.

Only fixed direct fields and saved mapping entries are read after all frozen
raw files are authenticated. Handles are opaque metadata; this module does
not follow graph edges, evaluate a source function, or prove rewrite semantics.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

PREFIX = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_"
COMPLETE = PREFIX + "n1_complete_future_functions.json.gz"
COMPLETE_CHECK = PREFIX + "n1_complete_future_functions_check.json"
F_ARCHIVE = PREFIX + "Rp_native_C5_pulse_source.json.gz"
F_CHECK = PREFIX + "Rp_native_C5_pulse_source_check.json"
P_REPORT = PREFIX + "native_physical_pulse_source.json"
PRIOR = "tmp/pf_noncommon_moment_v1/result.json"
BRIDGE = "tmp/pf_selector_bridge_v1/result.json"
PINS = {
    COMPLETE: "b63c1cdae91fa6b083a0bc7582c28bd5f9664225f223639b28eea8060ab7185b",
    COMPLETE_CHECK: "fd7fbf5ffd6305198c9f7ed4b14f2fa2d6e660d2f414a0debd78fc140812860f",
    F_ARCHIVE: "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe",
    F_CHECK: "320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15",
    P_REPORT: "8b21f8e868025d8cc95149ffab80b54603833f5f299ec9d100602af355cea82b",
    PRIOR: "fa5bfea33cfb082964772de2f113b54171c5942a55cdcdb75c9cd0eecde8ee20",
    BRIDGE: "f5f7462c99cd165f814edfeece04e39e0319682b953da41762e9f3d66247ee07",
}


class SelectorSourceError(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise SelectorSourceError(message)


def read_path(obj, *keys):
    for key in keys:
        require(isinstance(obj, dict) and key in obj, "missing direct source role: " + str(key))
        obj = obj[key]
    return obj


def saved_mapping(ref, mapping, branch, role, order):
    require(type(ref) is int and isinstance(mapping, dict), "malformed saved map input")
    out = {"branch": branch, "canonical_role": role, "ordinary_Z_order": order,
           "canonical_reference": ref, "reference_is_opaque_metadata": True,
           "mapped_coefficient_output_identity_proved": False,
           "canonical_owner": COMPLETE,
           "mapping_owner": P_REPORT if branch == "P" else F_ARCHIVE,
           "mapping_namespace": "physical_current_native_report" if branch == "P"
                                else "frame_Rp_native_C5_archive",
           "mapping_pointer": "component_result.substitution_receipt.changed_row_mapping"
                              if branch == "P" else "Rp_native_C5_pulse_source_bindings.rewritten_nodes"}
    key = str(ref)
    if key not in mapping:
        out.update({"status": "no_saved_mapping_unresolved", "mapped_reference": None})
    else:
        value = mapping[key]
        require(type(value) is int, "mapped source reference is not an opaque integer")
        out.update({"status": "saved_mapping_entry", "mapped_reference": value})
    return out


def verify_mapping_synthetic():
    missing = saved_mapping(7, {}, "F", "A0", 0)
    mapped = saved_mapping(7, {"7": 9}, "P", "A0", 0)
    checks = {"missing_stays_unresolved": missing["status"] == "no_saved_mapping_unresolved",
              "saved_entry_retained": mapped["mapped_reference"] == 9,
              "reference_not_value": mapped["reference_is_opaque_metadata"] is True}
    return {"checks": checks, "checks_passed": all(checks.values()),
            "synthetic_only": True, "source_values_computed": False}


def build_source_packet(root, manifest):
    root = Path(root).resolve()
    pins = {x["path"]: x["sha256"] for x in manifest["files"]}
    for path, expected in PINS.items():
        require(pins.get(path) == expected, "canonical manifest pin mismatch: " + path)
    raw = {}
    for item in manifest["files"]:
        rel = Path(item["path"])
        require(not rel.is_absolute() and ".." not in rel.parts, "bad frozen path")
        path = (root / rel).resolve()
        require(path.is_relative_to(root), "frozen path escapes repository")
        data = path.read_bytes()
        require(hashlib.sha256(data).hexdigest() == item["sha256"], "frozen source SHA mismatch")
        if item["path"] in PINS:
            raw[item["path"]] = data
    checks = {}
    for group, receipt_path, archive_path, flags in (
        ("complete", COMPLETE_CHECK, COMPLETE,
         ("current_original_n1_complete_future_energy_C1_function_installed",
          "current_original_n1_actual_selected_pulse_coefficient_C1_functions_installed")),
        ("F", F_CHECK, F_ARCHIVE,
         ("current_original_true_Rp_C5_selected_native_pulse_inlet_bridge_installed",
          "current_original_actual_Rp_C5_frame_equals_canonical_native_smooth_functions")),
    ):
        receipt = json.loads(raw[receipt_path])
        require(receipt.get("all_passed") is True and all(receipt.get(k) is True for k in flags),
                group + " source receipt gate failed")
        require(read_path(receipt, "input_hashes", Path(archive_path).name) == PINS[archive_path],
                group + " archive receipt pin mismatch")
        checks[group] = receipt
    complete = json.loads(gzip.decompress(raw[COMPLETE]))
    f_archive = json.loads(gzip.decompress(raw[F_ARCHIVE]))
    for group, archive in (("complete", complete), ("F", f_archive)):
        family = read_path(archive, "source_family")
        require(isinstance(family, dict) and bool(family) and
                family == read_path(checks[group], "source_family"),
                group + " source family mismatch")
    previous = json.loads(raw[PRIOR])
    require(previous.get("status") == "returned", "previous one-shot result not returned")
    previous_bridge = json.loads(raw[BRIDGE])
    require(previous_bridge.get("status") == "returned", "prior selector bridge not returned")
    bridge = read_path(previous_bridge, "result", "bridge")
    require(bridge.get("source_graph_traversed") is False and bridge.get("source_values_computed") is False
            and bridge.get("numeric_use_admissible") is False and bridge.get("certified") is False
            and bridge.get("P_F_function_equality_proved") is False,
            "prior selector bridge scope mismatch")
    prior_pins = read_path(bridge, "saved_receipt_cross_pins")
    require(len(prior_pins) == 4 and all(x.get("present") is True and x.get("matches") is True
                                        for x in prior_pins.values()), "prior receipt links missing")
    for name in ("C1", "C2"):
        record = read_path(bridge, "rows", name)
        require(record.get("canonical_first_two_equal_c3_saved_metadata") is True and
                record.get("c3_first_four_equal_high_saved_metadata") is True and
                record.get("high_first_four_equal_F_embedded_metadata") is True,
                "prior saved selector role chain is incomplete")
        for branch in ("P", "F"):
            entries = record.get(branch + "_saved_mapping")
            require(isinstance(entries, list) and len(entries) == 4 and
                    all(x.get("status") == "recorded_match" for x in entries),
                    "prior " + branch + " selector map incomplete")
    p = read_path(json.loads(raw[P_REPORT]), "component_result")
    require(p.get("passed") is True, "physical P report not accepted")
    p_audit = read_path(p, "source_prerequisite_audit")
    require(p_audit.get("authenticated:" + Path(F_ARCHIVE).name) is True,
            "P report lacks F archive prerequisite")
    bind = read_path(complete, "function_bindings")
    equations = read_path(bind, "pulse_equations")
    selected = read_path(bind, "selected_pulse")
    pb = read_path(p, "physical_pulse_bindings")
    fb = read_path(f_archive, "Rp_native_C5_pulse_source_bindings")
    pm = read_path(p, "substitution_receipt", "changed_row_mapping")
    fm = read_path(fb, "rewritten_nodes")
    coefficients = {}
    for role in ("A2", "A1", "A0", "target", "discriminant"):
        canonical = read_path(equations, role)
        require(isinstance(canonical, list) and len(canonical) >= 2 and
                all(type(x) is int for x in canonical[:2]), role + " C1 roles missing")
        coefficients[role] = {"canonical_C1_roles": canonical[:2],
            "P": [saved_mapping(ref, pm, "P", role, j) for j, ref in enumerate(canonical[:2])],
            "F": [saved_mapping(ref, fm, "F", role, j) for j, ref in enumerate(canonical[:2])]}
    canonical_ap = read_path(selected, "ap")
    apP = read_path(pb, "selected_pulse", "ap")
    apF = read_path(fb, "selected_pulse", "ap")
    require(all(isinstance(v, list) and len(v) >= 2 and all(type(x) is int for x in v[:2])
                for v in (canonical_ap, apP, apF)), "saved ap ordinary-Z roles missing")
    roots = {"canonical_C1_roles": canonical_ap[:2], "P_C1_roles": apP[:2], "F_C1_roles": apF[:2],
             "P_saved_map": [saved_mapping(ref, pm, "P", "ap", j) for j, ref in enumerate(canonical_ap[:2])],
             "F_saved_map": [saved_mapping(ref, fm, "F", "ap", j) for j, ref in enumerate(canonical_ap[:2])]}
    for branch, output in (("P", apP), ("F", apF)):
        for j, record in enumerate(roots[branch + "_saved_map"]):
            if record["status"] == "saved_mapping_entry":
                require(record["mapped_reference"] == output[j], branch + " saved ap map mismatch")
    weighted = read_path(bind, "energy", "weighted")
    delta_roles = {}
    for role in ("A2", "A1", "A0"):
        delta_roles[role] = []
        for j in (0, 1):
            left = coefficients[role]["P"][j]
            right = coefficients[role]["F"][j]
            delta_roles[role].append({"operation": "unevaluated_source_difference",
                "left_P_role": left, "right_F_role": right, "ordinary_Z_order": j,
                "both_saved_map_entries_present": left["status"] == "saved_mapping_entry" and
                                                    right["status"] == "saved_mapping_entry",
                "difference_value_computed": False, "ordinary_Z_semantics_proved": False})
    return {
        "schema": "pf_selector_quadratic_source_v1", "source_only": True,
        "source_values_computed": False, "source_graph_traversed": False,
        "roles_are_opaque_metadata": True, "weighted_future_role": weighted,
        "prior_selector_bridge": {"path": BRIDGE, "sha256": PINS[BRIDGE],
                                  "saved_role_and_receipt_metadata_checked": True,
                                  "function_equality_proved": False},
        "quadratic_coefficient_roles": coefficients, "selected_root_roles": roots,
        "coefficient_difference_roles": delta_roles,
        "source_identity_proved": False, "P_F_function_equality_proved": False,
        "ordinary_Z_transport_proved": False, "capture_free_proved": False,
        "difference_provider_admissible_for_numeric_use": False, "old_errors_removed": False,
        "conditional_formula": "Delta ap=-(Delta A2*apF^2+Delta A1*apF+Delta A0)/(A2P*(apP+apF)+A1P)",
        "coefficient_transport_contract": {
            "formula": "Delta C_j=Delta u_j+ap_P*Delta v_j+Delta ap*v_F,j",
            "ordinary_Z1": "Delta u_j,Z+ap_P,Z*Delta v_j+ap_P*Delta v_j,Z+Delta ap_Z*v_F,j+Delta ap*v_F,j,Z",
            "u_v_source_owners_bound": False, "signed_difference_provider_bound": False,
            "ordinary_Z_semantics_proved": False, "source_values_computed": False,
        },
        "unresolved": ["missing saved coefficient maps", "positive-root and nonzero denominator on both actual branches",
                       "P/F coefficient difference providers", "capture-free ordinary-Z transport",
                       "common physical coordinate and certified errors"],
    }
