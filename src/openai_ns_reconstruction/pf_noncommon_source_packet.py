"""Source-bound formal P/F moment DAG with explicit non-common jets.

Only direct records from frozen prior results are read. Opaque saved handles
are bound to symbolic names as metadata; no handle is evaluated as a value.
The algebra is conditional on separately proving the listed delta identities,
ordinary-Z semantics, input owners and coordinate/normalization contracts.
"""
from __future__ import annotations

import hashlib
import json
from math import comb
from pathlib import Path

BRIDGE = "tmp/pf_selector_bridge_v1/result.json"
V2 = "tmp/pulse_gap_source_factorized_delta_v2/result.json"
FREE_NAMES = "tmp/frame_replacement_free_names_v3/result.json"
BRIDGE_SHA = "f5f7462c99cd165f814edfeece04e39e0319682b953da41762e9f3d66247ee07"
V2_SHA = "a2f1e5edc3ebe60f237b6aafc57fd7bd36536fca7ab5697a4ea3e37a797e56e6"
FREE_NAMES_SHA = "134ec429c5f0608a63894c20cd712d0e465c633e881a5d4e86fec1c532e74fff"


class PacketError(RuntimeError):
    pass


def need(condition, label):
    if not condition:
        raise PacketError(label)


def row_identity(row, owner, pointer, order):
    namespace = {"P": "physical_current_native_report", "F": "frame_Rp_native_C5_archive"}[owner]
    need(row.get("status") == "direct_source_ref" and row.get("owner") == owner and
         row.get("namespace") == namespace and row.get("source_role") == pointer and
         row.get("ordinary_Z_order") == order and type(row.get("source_node_id")) is int and
         row.get("node_id_is_not_a_numeric_value") is True and row.get("source_value_selected") is False,
         owner + " row owner/path/opaque metadata mismatch")


def row_set(record, owner):
    need(record.get("status") == "direct_source_rows", owner + " u rows missing")
    rows = record.get("all_saved_rows")
    need(isinstance(rows, list) and len(rows) >= 4, owner + " u jet too short")
    base = {"P": "$.component_result.physical_pulse_bindings", "F": "$.Rp_native_C5_pulse_source_bindings"}[owner]
    for i, row in enumerate(rows[:4]):
        row_identity(row, owner, f"{base}.current_frame.u[{i}]", i)
    return rows[:4]


def name_ref(name):
    return {"op": "symbol", "name": name}


def product(*args):
    return {"op": "product", "args": list(args)}


def addition(*args):
    return {"op": "sum", "args": list(args)}


def negative(arg):
    return {"op": "negate", "arg": arg}


def term(n, k, a, b):
    return product({"op": "exact_integer", "value": comb(n, k)}, name_ref(a), name_ref(b))


def build_source_packet(root, manifest):
    root = Path(root).resolve()
    pinned = {r["path"]: r["sha256"] for r in manifest["files"]}
    for path, digest in ((BRIDGE, BRIDGE_SHA), (V2, V2_SHA), (FREE_NAMES, FREE_NAMES_SHA)):
        need(pinned.get(path) == digest, "prior result is not pinned: " + path)
    raw = {}
    # Verify all frozen bytes before parsing either prior result.
    for record in manifest["files"]:
        relative = Path(record["path"])
        need(not relative.is_absolute() and ".." not in relative.parts, "bad manifest path")
        path = (root / relative).resolve()
        need(path.is_relative_to(root), "manifest path escapes workspace")
        data = path.read_bytes()
        need(hashlib.sha256(data).hexdigest() == record["sha256"], "frozen SHA mismatch")
        if record["path"] in (BRIDGE, V2, FREE_NAMES):
            raw[record["path"]] = data
    bridge_result = json.loads(raw[BRIDGE])
    prior_result = json.loads(raw[V2])
    lexical_result = json.loads(raw[FREE_NAMES])
    need(bridge_result.get("status") == "returned" and prior_result.get("status") == "returned",
         "prior result not returned")
    bridge = bridge_result["result"]["bridge"]
    prior = prior_result["result"]
    lexical = lexical_result["result"]
    need(lexical_result.get("status") == "returned" and lexical.get("saved_pair") == [31655, 265772]
         and lexical.get("actual_capture_proved") is False and
         lexical.get("actual_capture_free_proved") is False and
         lexical.get("ordinary_Z_transport_proved") is False,
         "saved lexical evidence scope mismatch")
    analysis = lexical["replacement_analysis"]
    need(analysis["names"]["Z"]["state"] == "unbound_syntactic_witness",
         "saved syntactic Z witness missing")
    need(bridge.get("source_only") is True and bridge.get("source_graph_traversed") is False and
         bridge.get("source_values_computed") is False and bridge.get("numeric_use_admissible") is False and
         bridge.get("P_F_function_equality_proved") is False and
         prior.get("source_only") is True and prior.get("formal_source_expression_defined") is True and
         prior.get("numeric_values_computed") is False and prior.get("source_values_selected") is False and
         prior.get("scope_audit", {}).get("source_graph_traversed") is False and
         prior.get("admissible_for_numeric_use") is False,
         "prior source scope mismatch")
    cpins = bridge["saved_receipt_cross_pins"]
    need(all(v.get("present") is True and v.get("matches") is True for v in cpins.values())
         and len(cpins) == 4, "recorded receipt chain incomplete")
    symbols = {}
    for side in ("P", "F"):
        side_data = prior["sides"][side]
        u = row_set(side_data["current_frame"]["u"], side)
        for j, item in enumerate(u):
            symbols[f"u_{side}_{j}"] = {"kind": "opaque_saved_role", "source": V2,
                                           "role": item, "ordinary_Z_order": j}
        for r in (1, 2):
            label = f"C{r}"
            saved = bridge["rows"][label]
            refs = saved[side + "_output_role"]
            mappings = saved[side + "_saved_mapping"]
            need(len(refs) == 4 and len(mappings) == 4 and
                 all(m.get("status") == "recorded_match" for m in mappings),
                 side + label + " saved mapping incomplete")
            prior_rows = side_data["selected_pulse"][label]["all_saved_rows"]
            need(len(prior_rows) >= 4, side + label + " prior rows incomplete")
            for j, ref in enumerate(refs):
                item = prior_rows[j]
                base = {"P": "$.component_result.physical_pulse_bindings", "F": "$.Rp_native_C5_pulse_source_bindings"}[side]
                row_identity(item, side, f"{base}.selected_pulse.{label}[{j}]", j)
                need(type(ref) is int and item.get("source_node_id") == ref,
                     side + label + " role metadata differs")
                symbols[f"C_{side}_{r}_{j}"] = {"kind": "opaque_saved_role",
                    "source": V2, "role": item, "ordinary_Z_order": j}
        for j in range(4):
            symbols[f"du_{j}"] = {"kind": "formal_unproved_difference",
                                    "claim": f"u_P_{j}-u_F_{j}"}
            for r in (1, 2):
                symbols[f"w_{side}_{r}_{j}"] = {"kind": "unbound_branch_weight_jet",
                    "owner": side, "ordinary_Z_order": j,
                    "kernel_weight_candidate": bridge["kernel_first_row_roles"][r-1] if j == 0 else None,
                    "candidate_owner_binding_proved": False}
                symbols[f"dw_{r}_{j}"] = {"kind": "formal_unproved_difference",
                                             "claim": f"w_P_{r}_{j}-w_F_{r}_{j}"}
                symbols[f"dC_{r}_{j}"] = {"kind": "formal_unproved_difference",
                                             "claim": f"C_P_{r}_{j}-C_F_{r}_{j}"}
    definitions = {}
    for n in range(4):
        definitions[f"cF_{n}"] = addition(*(term(n, k, f"w_F_{r}_{k}", f"C_F_{r}_{n-k}")
                                             for r in (1, 2) for k in range(n+1)))
        definitions[f"dc_{n}"] = addition(*(x for r in (1, 2) for k in range(n+1) for x in (
            term(n, k, f"w_P_{r}_{k}", f"dC_{r}_{n-k}"),
            term(n, k, f"dw_{r}_{k}", f"C_F_{r}_{n-k}"))))
        definitions[f"mF_{n}"] = negative(addition(*(term(n, k, f"u_F_{k}", f"cF_{n-k}")
                                                     for k in range(n+1))))
        definitions[f"dm_{n}"] = negative(addition(*(x for k in range(n+1) for x in (
            term(n, k, f"u_P_{k}", f"dc_{n-k}"), term(n, k, f"du_{k}", f"cF_{n-k}")))))
    return {
        "schema": "pf_noncommon_source_packet_v1", "source_only": True,
        "symbols": symbols, "definitions": definitions,
        "formula": "dc=Σ(wP*dC+dw*CF); dm=−(uP*dc+du*cF), with full ordinary-Z Leibniz rows",
        "source_values_computed": False, "source_graph_traversed": False,
        "saved_capture_evidence": {
            "source": FREE_NAMES, "saved_pair": lexical["saved_pair"],
            "Z_state": analysis["names"]["Z"]["state"],
            "unresolved_scope_states": len(analysis["unknown"]),
            "prior_audit_orders": [0, 1, 2], "Z3_source_rows_uncovered": True,
            "actual_capture_proved": False, "capture_free_proved": False,
            "ordinary_Z_transport_proved": False,
        },
        "numeric_use_admissible": False, "certified": False, "old_errors_removed": False,
        "unresolved": ["difference identities and correlated providers", "P/F weight owner and Z jets",
                       "P/F inlet and selected coefficient derivative semantics", "capture-free substitution",
                       "delta, scale, physical radius and finite endpoint owners", "certified error transport"],
        "scope": "ordinary Z0–Z3; no temporal coefficient recursion",
        "needs_m4_for_normalized_Q_Z3": True,
    }


def verify_packet_synthetic(packet):
    """Check this emitted DAG against independent polynomial derivatives."""
    import sympy as sp
    z = sp.symbols("z")
    def polynomial(seed):
        return sum(sp.Rational(seed + j + 1, 3 + j) * z**j for j in range(4))
    branch = {}
    for side, offset in (("P", 4), ("F", 1)):
        branch[(side, "u")] = polynomial(offset)
        for r in (1, 2):
            branch[(side, "w", r)] = polynomial(offset + r * 2)
            branch[(side, "C", r)] = polynomial(offset + r * 4)
    def jet(expr, n):
        return sp.diff(expr, z, n).subs(z, 0)
    values = {}
    for n in range(4):
        for side in ("P", "F"):
            values[f"u_{side}_{n}"] = jet(branch[(side, "u")], n)
            for r in (1, 2):
                for kind in ("w", "C"):
                    values[f"{kind}_{side}_{r}_{n}"] = jet(branch[(side, kind, r)], n)
        values[f"du_{n}"] = values[f"u_P_{n}"] - values[f"u_F_{n}"]
        for r in (1, 2):
            values[f"dw_{r}_{n}"] = values[f"w_P_{r}_{n}"] - values[f"w_F_{r}_{n}"]
            values[f"dC_{r}_{n}"] = values[f"C_P_{r}_{n}"] - values[f"C_F_{r}_{n}"]
    defs = packet["definitions"]
    def evaluate(node):
        op = node["op"]
        if op == "symbol":
            name = node["name"]
            return values[name] if name in values else evaluate(defs[name])
        if op == "exact_integer": return sp.Integer(node["value"])
        if op == "negate": return -evaluate(node["arg"])
        if op == "sum": return sum((evaluate(x) for x in node["args"]), sp.Integer(0))
        if op == "product": return sp.prod(evaluate(x) for x in node["args"])
        raise PacketError("unknown synthetic DAG operation")
    cP = sum(branch[("P", "w", r)] * branch[("P", "C", r)] for r in (1, 2))
    cF = sum(branch[("F", "w", r)] * branch[("F", "C", r)] for r in (1, 2))
    mP = -branch[("P", "u")] * cP
    mF = -branch[("F", "u")] * cF
    checks = {}
    for n in range(4):
        for name, expected in (("cF", cF), ("dc", cP-cF), ("mF", mF), ("dm", mP-mF)):
            checks[f"{name}_Z{n}_exact"] = sp.simplify(evaluate(defs[f"{name}_{n}"]) - jet(expected, n)) == 0
    checks["nonzero_weight_difference"] = values["dw_1_0"] != 0
    checks["nonzero_inlet_difference"] = values["du_0"] != 0
    wrong_owner = dict(packet["symbols"]["u_P_0"]["role"])
    wrong_owner["namespace"] = "frame_Rp_native_C5_archive"
    try:
        row_identity(wrong_owner, "P", "$.component_result.physical_pulse_bindings.current_frame.u[0]", 0)
    except PacketError:
        checks["wrong_namespace_rejected"] = True
    else:
        checks["wrong_namespace_rejected"] = False
    checks["branch_weights_remain_unbound"] = all(
        packet["symbols"][f"w_{side}_{r}_0"]["candidate_owner_binding_proved"] is False
        for side in ("P", "F") for r in (1, 2)
    )
    return {"checks": checks, "checks_passed": all(checks.values()),
            "synthetic_only": True, "exact_arithmetic": True,
            "source_values_computed": False, "certified": False}
