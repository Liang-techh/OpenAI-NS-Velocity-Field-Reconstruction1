"""Direct, source-only selector-row bridge; never traverses saved graphs.

All canonical bytes are hashed before gzip decoding. Saved handles and mapping
entries are opaque metadata; equality is only evidence about saved records,
not function equality, derivative transport, numerical values or P/F ownership.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

PREFIX = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_"
STEMS = {
    "complete": PREFIX + "n1_complete_future_functions",
    "c3": PREFIX + "n1_pulse_C3_functions",
    "high": PREFIX + "n1_pulse_pressure_C3_functions",
    "f": PREFIX + "Rp_native_C5_pulse_source",
    "kernel": PREFIX + "n1_pulse_kernel_functions",
}
P_REPORT = PREFIX + "native_physical_pulse_source.json"
PINS = {
    STEMS["complete"] + ".json.gz": "b63c1cdae91fa6b083a0bc7582c28bd5f9664225f223639b28eea8060ab7185b",
    STEMS["complete"] + "_check.json": "fd7fbf5ffd6305198c9f7ed4b14f2fa2d6e660d2f414a0debd78fc140812860f",
    STEMS["c3"] + ".json.gz": "35c683059732fbb2e9e1cdf79ad86533db472285023a4ec339f00f9920ce6143",
    STEMS["c3"] + "_check.json": "adae772b055774a6c360d4a0760d866cdaf5cf1bf80993aeed8b98c849deedf9",
    STEMS["high"] + ".json.gz": "bdfcc70f57d8cea981f508074690a7e940f884382c2452805f126a3f96bdb59d",
    STEMS["high"] + "_check.json": "cbac4685fc9031a2ed5c5a7ef624b438e1a4ce3f32b85b569512aecb1d38feb8",
    STEMS["f"] + ".json.gz": "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe",
    STEMS["f"] + "_check.json": "320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15",
    STEMS["kernel"] + ".json.gz": "9584201b517fe3058ea6d1a169ddb7d319776eb97479b36024af4bffd4973b75",
    STEMS["kernel"] + "_check.json": "47343a839590c3cb08c8d726ef225415e15f4fb12a5f7cf888049ea1f20cffc0",
    P_REPORT: "8b21f8e868025d8cc95149ffab80b54603833f5f299ec9d100602af355cea82b",
}
GATES = {
    "complete": ("current_original_n1_complete_future_energy_C1_function_installed",
                 "current_original_n1_actual_selected_pulse_coefficient_C1_functions_installed"),
    "c3": ("current_original_n1_actual_future_and_pulse_C3_functions_installed",),
    "high": ("current_original_n1_actual_native_and_Gamma_C5_functions_installed",
             "current_original_n1_actual_leading_pulse_U_mean_C5_Q_C4_functions_installed",
             "current_original_n1_true_pulse_pressure_C3_functions_installed"),
    "f": ("current_original_true_Rp_C5_selected_native_pulse_inlet_bridge_installed",
          "current_original_actual_Rp_native_C5_all_eight_source_functions_consumed",
          "current_original_actual_Rp_C5_frame_equals_canonical_native_smooth_functions"),
    "kernel": ("current_original_n1_exact_pulse_kernel_functions_installed",),
}


class BridgeError(RuntimeError):
    pass


def require(test, message):
    if not test:
        raise BridgeError(message)


def direct(obj, *keys):
    for key in keys:
        require(isinstance(obj, dict) and key in obj, "missing direct field " + str(key))
        obj = obj[key]
    return obj


def rows(obj, label, minimum):
    require(isinstance(obj, list) and len(obj) >= minimum and
            all(type(i) is int for i in obj), label + " row metadata malformed")
    return obj


def mapping_evidence(old, output, mapping):
    """Compare saved map entry in its own declared namespace; never guess missing."""
    require(isinstance(mapping, dict), "saved mapping is not an object")
    key = str(old)
    if key not in mapping:
        return {"old": old, "output": output, "status": "unrecorded_unproved"}
    mapped = mapping[key]
    return {"old": old, "output": output, "saved_mapping": mapped,
            "status": "recorded_match" if mapped == output else "recorded_mismatch"}


def synthetic_negative_fixture():
    checks = {
        "mapping_match": mapping_evidence(4, 8, {"4": 8})["status"] == "recorded_match",
        "mapping_mismatch": mapping_evidence(4, 9, {"4": 8})["status"] == "recorded_mismatch",
        "missing_not_inferred": mapping_evidence(4, 4, {})["status"] == "unrecorded_unproved",
    }
    try:
        rows([1, 2], "synthetic", 4)
    except BridgeError:
        checks["short_jet_rejected"] = True
    else:
        checks["short_jet_rejected"] = False
    return {"checks": checks, "checks_passed": all(checks.values()),
            "source_values_computed": False, "certified": False}


def build_bridge(root, manifest):
    root = Path(root).resolve()
    manifest_files = {x["path"]: x["sha256"] for x in manifest["files"]}
    for path, expected in PINS.items():
        require(manifest_files.get(path) == expected, "canonical manifest pin mismatch: " + path)
    raw = {}
    # Authenticate every frozen raw file before any decompression.
    for item in manifest["files"]:
        relative = Path(item["path"])
        require(not relative.is_absolute() and ".." not in relative.parts, "invalid frozen path")
        path = (root / relative).resolve()
        require(path.is_relative_to(root), "frozen path escapes repository")
        data = path.read_bytes()
        require(hashlib.sha256(data).hexdigest() == item["sha256"], "frozen raw SHA mismatch")
        if item["path"] in PINS:
            raw[item["path"]] = data

    receipts = {}
    for group, stem in STEMS.items():
        receipt = json.loads(raw[stem + "_check.json"])
        require(receipt.get("all_passed") is True, group + " receipt not accepted")
        for flag in GATES[group]:
            require(receipt.get(flag) is True, group + " receipt gate missing: " + flag)
        own = Path(stem + ".json.gz").name
        require(direct(receipt, "input_hashes", own) == PINS[stem + ".json.gz"],
                group + " own archive pin mismatch")
        receipts[group] = receipt

    # Values decoded below are only named direct fields, never graph nodes.
    archives = {group: json.loads(gzip.decompress(raw[stem + ".json.gz"]))
                for group, stem in STEMS.items()}
    for group in STEMS:
        require(direct(archives[group], "source_family") == direct(receipts[group], "source_family"),
                group + " receipt/archive source family mismatch")
    pdoc = json.loads(raw[P_REPORT])
    p = direct(pdoc, "component_result")
    require(p.get("passed") is True, "P source report not accepted")
    audit = direct(p, "source_prerequisite_audit")
    for group in ("high", "f"):
        flag = "authenticated:" + Path(STEMS[group] + ".json.gz").name
        require(audit.get(flag) is True, "P report did not authenticate " + group)
    complete = direct(archives["complete"], "function_bindings")
    c3 = direct(archives["c3"], "function_bindings", "functions", "pulse")
    high = direct(archives["high"], "C5_pulse_C3_pressure_bindings")
    f = direct(archives["f"], "Rp_native_C5_pulse_source_bindings")
    high_selected = direct(high, "selected_pulse")
    f_high = direct(f, "inputs", Path(STEMS["high"] + ".json.gz").name, "selected_pulse")
    f_selected = direct(f, "selected_pulse")
    p_selected = direct(p, "physical_pulse_bindings", "selected_pulse")
    p_map = direct(p, "substitution_receipt", "changed_row_mapping")
    f_map = direct(f, "rewritten_nodes")
    result_rows = {}
    for name in ("C1", "C2"):
        canonical = rows(direct(complete, "selected_pulse", name), name + " canonical", 2)
        c = rows(direct(c3, name), name + " C3", 4)
        h = rows(direct(high_selected, name), name + " HIGH", 4)
        fh = rows(direct(f_high, name), name + " F embedded HIGH", 4)
        fp = rows(direct(f_selected, name), name + " F output", 4)
        pp = rows(direct(p_selected, name), name + " P output", 4)
        result_rows[name] = {
            "canonical_C1_role": canonical[:2], "c3_role": c[:4], "high_role": h[:4],
            "F_embedded_high_role": fh[:4], "F_output_role": fp[:4], "P_output_role": pp[:4],
            "canonical_first_two_equal_c3_saved_metadata": canonical[:2] == c[:2],
            "c3_first_four_equal_high_saved_metadata": c[:4] == h[:4],
            "high_first_four_equal_F_embedded_metadata": h[:4] == fh[:4],
            "F_saved_mapping": [mapping_evidence(h[i], fp[i], f_map) for i in range(4)],
            "P_saved_mapping": [mapping_evidence(h[i], pp[i], p_map) for i in range(4)],
        }

    cross_pins = {}
    for child, parent in (("complete", "kernel"), ("c3", "complete"), ("high", "c3"), ("f", "high")):
        name = Path(STEMS[parent] + ".json.gz").name
        saved = receipts[child].get("input_hashes", {}).get(name)
        require(saved is None or saved == PINS[STEMS[parent] + ".json.gz"],
                child + " receipt has contradictory " + parent + " archive pin")
        cross_pins[child + "_pins_" + parent] = {
            "present": saved is not None,
            "matches": saved == PINS[STEMS[parent] + ".json.gz"],
        }
    return {
        "schema": "pf_selector_direct_bridge_v1", "source_only": True,
        "raw_authenticated_before_decompression": True, "source_graph_traversed": False,
        "source_values_computed": False, "certified": False,
        "weighted_future_role": direct(complete, "energy", "weighted"),
        "kernel_first_row_roles": direct(archives["kernel"], "exact_pulse_kernel_functions", "first_row"),
        "kernel_incoming_roles": direct(archives["kernel"], "actual_native_incoming_C1_functions"),
        "roles_are_opaque_metadata": True, "rows": result_rows,
        "saved_receipt_cross_pins": cross_pins,
        "P_F_function_equality_proved": False, "capture_free_Z0_Z3_proved": False,
        "common_incoming_W_future_owner_proved": False,
        "F_historical_radius_and_endpoints_proved": False,
        "physical_scale_and_error_transport_proved": False,
        "old_errors_removed": False, "numeric_use_admissible": False,
    }
