"""Bind actual modified moment sums to admitted leading and Duhamel germs.

Only ordinary Z0/Z1 values and first physical-native derivatives of the five
moment fields are admitted.  Linearity is applied to the actual sum nodes.
No old reducer, quadrature, source value or error cap is evaluated.
"""
from __future__ import annotations

import copy
import hashlib
import json
from fractions import Fraction
from pathlib import Path

FIELDS = ("m", "h", "k", "e", "p")
ORDERS = (0, 1)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def require(value, message):
    if not value:
        raise ValueError(message)


def endpoint(nodes, root, body, coordinate):
    row = nodes[root]
    require(row["operation"] == "function_substitution" and row["expression"] == body
            and row["variable"] == coordinate and row["value"] == 1
            and row.get("Z_independent_substitution") is True, "leading endpoint wrapper mismatch")


def derive(nodes, left_chart, right_chart, seed_admission, seed_inputs, raw_native, history_seam, *, expected_nodes_digest):
    require(digest(nodes) == expected_nodes_digest, "modified moment source graph changed")
    require(seed_admission["passed"] is True and raw_native["passed"] is True
            and history_seam["passed"] is True, "leading/history premises incomplete")
    require(raw_native["seed_proof_digest"] == digest(seed_admission), "raw leading proof uses another seed")
    raw_graph = raw_native["source_graph"]
    require(digest(nodes[:raw_graph["node_count"]]) == raw_graph["nodes_digest"], "accepted raw source prefix changed")
    delta = history_seam["graph_delta"]
    base_count = delta["base_graph_node_count"]
    require(digest(nodes[:base_count]) == delta["base_nodes_sha256"]
            and nodes[base_count:] == delta["appended_nodes"]
            and digest(nodes) == delta["full_nodes_sha256"], "history source definition delta changed")
    require(left_chart["coordinate"] == 84528 and right_chart["coordinate"] == 185676, "foreign moment chart")
    seed_rows = {(r["field"], r["order"]):r for r in seed_admission["fields"]}
    native_rows = {(r["field"], r["order"]):r for r in raw_native["row_forms"]}
    history_rows = {(r["field"], r["ordinary_Z_order"]):r for r in history_seam["rows"]}
    require(len(seed_rows) == len(seed_admission["fields"]) and len(native_rows) == len(raw_native["row_forms"])
            and len(history_rows) == len(history_seam["rows"]), "duplicate admitted field/order rows")
    require(set(history_rows) == {(f, k) for f in FIELDS for k in ORDERS}, "history seam rows incomplete")
    appended, rows = [], []
    def append(row):
        handle = len(nodes) + len(appended)
        appended.append(row)
        return handle
    for field in FIELDS:
        for order in ORDERS:
            key = field, order
            seed, radial, history = seed_rows[key], native_rows[key], history_rows[key]
            require(radial["status"] == "proved" and radial["left_form"] == radial["right_form"], "raw native row not admitted")
            require(all(isinstance(radial[key], int) and 0 <= radial[key] < raw_graph["node_count"]
                        for key in ("left_node", "right_node")), "raw derivative handle exceeds admitted raw prefix")
            require(radial["left_form_sha256"] == digest(radial["left_form"])
                    and radial["right_form_sha256"] == digest(radial["right_form"]), "raw derivative form digest changed")
            require(history["status"] == "source_identity_proved", "history native seam not admitted")
            pairing = history["native_pairing"]
            require(all(pairing.get(key) is True for key in ("same_density_source", "same_rate", "same_incoming_source",
                        "post_density_handle_binding_exact", "accepted_patch_density_handle_binding_exact")), "paired history derivative premises incomplete")
            incoming = history["post_history_endpoint_source"]
            endpoint(nodes, incoming, history["post_history_body"], 84528)
            endpoint(nodes, history["patch_history_endpoint_function_handle"], history["patch_history_body"], 185676)
            records = {}
            for side, chart, coordinate, closed_key, raw_key, history_key, native_history_key in (
                ("postrestore", left_chart, 84528, "inlet_node", "left_node", "post_history_body", "post_native_y_definition"),
                ("patch", right_chart, 185676, "patch_node", "right_node", "patch_history_body", "accepted_patch_native_function_handle"),
            ):
                modified = chart["modified_fields"][field][order]
                hbody = chart["history"][field][order]
                require(hbody == history[history_key], "modified moment history body changed")
                descriptor = nodes[modified]
                require(descriptor["operation"] == "sum" and len(descriptor["arguments"]) == 2
                        and descriptor["arguments"].count(hbody) == 1, "actual modified moment is not leading + history")
                leading = next(ref for ref in descriptor["arguments"] if ref != hbody)
                closed = seed[closed_key]
                endpoint(nodes, closed, leading, coordinate)
                if side == "postrestore":
                    require(leading == seed_inputs["postrestore"]["fields"][field][order], "postrestore leading body changed")
                else:
                    require(closed == seed_inputs["leading_endpoint"]["fields"][field][order], "patch leading endpoint seed port changed")
                if field == "m":
                    require(seed["status"] == "external_premise", "leading mean source status changed")
                    mean = seed_inputs["saved_mean_certificate"]
                    require(mean["passed"] is True and mean["source_only"] is True
                            and mean["numeric_values_selected"] is False and mean["numeric_caps_used"] is False,
                            "saved mean premise missing")
                    mean_row = next(r for r in mean["rows"] if r["ordinary_Z_order"] == order)
                    require(mean_row["passed"] is True and mean_row[raw_key] == leading
                            and mean_row["left_function_sha256"] == mean_row["right_function_sha256"], "saved mean source port changed")
                else:
                    require(seed["status"] == "proved" and seed["inlet"] == seed["patch"], "leading value row not admitted")
                source_derivative = radial[raw_key]
                h_derivative = history[native_history_key]
                require(isinstance(h_derivative, int) and 0 <= h_derivative < len(nodes), "history derivative handle invalid")
                derivative_node = nodes[h_derivative]
                require(derivative_node["operation"] == "sum" and len(derivative_node["arguments"]) == 2,
                        "history derivative is not the admitted endpoint source sum")
                density_endpoint, incoming_product = derivative_node["arguments"]
                density_body = chart["finite_N"]["five_signed_density_increments"][field][order]
                endpoint(nodes, density_endpoint, density_body, coordinate)
                product = nodes[incoming_product]
                require(product["operation"] == "product" and len(product["arguments"]) == 2
                        and product["arguments"][1] == incoming, "history derivative incoming source changed")
                coefficient = nodes[product["arguments"][0]]
                require(coefficient["operation"] == "exact_rational"
                        and Fraction(coefficient["numerator"], coefficient["denominator"]) == -Fraction(history["rate"]),
                        "history derivative rate changed")
                if side == "postrestore":
                    require(base_count <= h_derivative < len(nodes), "post derivative is not the new history definition")
                else:
                    require(h_derivative < base_count and h_derivative == history["accepted_patch_native_y_handle"],
                            "patch derivative does not reference the admitted first-native root")
                value = append(dict(operation="function_substitution", expression=modified, variable=coordinate, value=1,
                    Z_independent_substitution=True, guarded=True, field=field, ordinary_Z_order=order,
                    source_formula="actual modified moment endpoint", original_chart=side))
                leading_derivative = append(dict(operation="function_substitution", expression=source_derivative,
                    variable=coordinate, value=1, Z_independent_substitution=True, guarded=True,
                    field=field, ordinary_Z_order=order, source_formula="admitted physical D_y leading moment endpoint", original_chart=side))
                native = append(dict(operation="sum", arguments=[leading_derivative, h_derivative],
                    source_formula="D_y(leading+history)=D_y leading+D_y history", field=field,
                    ordinary_Z_order=order, original_chart=side, no_additional_Jacobian=True))
                records[side] = dict(actual_modified_body=modified, leading_body=leading, history_body=hbody,
                    leading_value_endpoint=closed, leading_native_body=source_derivative,
                    history_native_endpoint=h_derivative, modified_value_endpoint_definition=value,
                    modified_native_endpoint_definition=native)
            rows.append(dict(field=field, order=order, ordinary_Z_order=order, status="proved",
                value_identity_proved=True, native_identity_proved=True, **records,
                value_status="proved_from_saved_mean_premise" if field == "m" else "proved_from_admitted_leading_value",
                native_status="proved_from_admitted_leading_and_history_native_sources",
                value_reason="same admitted leading value plus same actual history endpoint",
                native_reason="same admitted leading physical-native germ plus new same-source history derivative",
                leading_value_status_preserved=seed["status"],
                leading_mean_proof_ref=dict(report_sha256=seed_inputs["saved_mean_report_sha256"],
                    certificate_digest=digest(seed_inputs["saved_mean_certificate"])) if field == "m" else None,
                raw_native_row_digest=digest(radial),
                history_seam_row_digest=digest(history)))
    return dict(schema="actual_Rm_modified_moment_value_native_germ_v1", passed=True, source_only=True,
        actual_saved_graph_run=False, actual_source_graph_inspected=True,
        nodes_digest=expected_nodes_digest, module_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_admission_digests=dict(seed=digest(seed_admission), raw_native=digest(raw_native), history_seam=digest(history_seam)),
        ordinary_orders=list(ORDERS), fields=list(FIELDS), rows=rows, paired_rows=len(rows),
        graph_delta=dict(base_graph_node_count=len(nodes), base_nodes_sha256=expected_nodes_digest,
            appended_nodes=copy.deepcopy(appended), appended_node_count=len(appended),
            full_node_count=len(nodes)+len(appended), full_nodes_sha256=digest(nodes+appended)),
        moment_value_native_seam_proved=True, E_V_native_seam_proved=False,
        slow_native_loop_seam_proved=False, physical_RQy_seam_proved=False,
        numeric_values_computed=False, old_mathematical_operator_replayed=False,
        ledger_packet_not_received=True, ledger_unchanged_verification_required_by_runner=True, full_goal_complete=False)
