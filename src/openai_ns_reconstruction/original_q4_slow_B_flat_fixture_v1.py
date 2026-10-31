"""Fresh synthetic fixture for the original Q4 slow-B flat binding.

The fixture uses opaque integer graph handles only.  It does not import a
producer, decode an archive, evaluate a ``FunctionRef`` or walk a real graph.
"""


SCHEMA = "original_q4_slow_B_flat_fixture_v1"

Q4_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
)
LOOP_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
)
RAW_SHA = {
    Q4_ARCHIVE: "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a",
    LOOP_ARCHIVE: "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
}


def _build_case():
    """Build one entirely fresh synthetic Q4/loop/receipt quartet."""

    family = {"fixture": "q4_slow_B_flat_distinct_family_v1"}
    graph = [
        {"operation": "exact_rational", "numerator": 0, "denominator": 1},
        {"operation": "exact_rational", "numerator": 1, "denominator": 1},
    ]

    def add(row):
        graph.append(row)
        return len(graph) - 1

    coordinate = add({"operation": "bound_variable", "name": "first_micro_coordinate"})
    jacobian = add({"operation": "fixture_physical_measure", "name": "first_micro_J"})
    N = add({
        "operation": "exact_positive_integer_expression",
        "name": "common_N",
        "one_common_N_for_all_charts_and_control_denominators": True,
    })
    predicate = add({
        "operation": "flat_predicate",
        "name": "first_micro_flat_predicate",
        "conditional_on": "sc_positive",
    })
    collar = add({"operation": "fixture_original_collar_upper", "name": "collar_upper"})
    pstar = add({"operation": "analytic_unary", "name": "exp", "argument": 1})
    invjac = add({
        "operation": "positive_quotient",
        "numerator": 1,
        "denominator": jacobian,
    })
    invN = add({
        "operation": "positive_quotient",
        "numerator": 1,
        "denominator": N,
    })

    bphi = [
        add({"operation": "fixture_loop_B_over_Pstar_phi", "ordinary_order": k})
        for k in range(3)
    ]

    derivatives = []
    slow_B_y = []
    for k in range(3):
        derivative = add({
            "operation": "original_lazy_flat_branch",
            "flat_value": 0,
            "flat_predicate": predicate,
            "original_collar_coordinate": coordinate,
            "original_collar_upper": collar,
            "same_original_lazy_native_derivative": True,
            "ordinary_slow_Z_derivative_order": k,
            "same_original_flat_extension_through_order5": True,
            "active_body_not_evaluated_on_flat_branch": True,
        })
        derivatives.append(derivative)
        slow_B_y.append(add({
            "operation": "product",
            "arguments": [invjac, derivative],
        }))

    background = [
        add({"operation": "fixture_radial_background_V", "ordinary_order": k})
        for k in range(3)
    ]
    scaled_slow_B = [
        add({"operation": "product", "arguments": [slow_B_y[k], invN]})
        for k in range(3)
    ]
    total_V_y = [
        add({
            "operation": "sum",
            "arguments": [background[k], scaled_slow_B[k], bphi[k]],
        })
        for k in range(3)
    ]
    # These are opaque graph references to physical rows, not evaluated values.
    physical_U_y = [
        add({"operation": "product", "arguments": [total_V_y[k], pstar]})
        for k in range(3)
    ]

    loop_chart = {
        "coordinate": coordinate,
        "physical_measure": jacobian,
        "physical_scale": pstar,
        "loop": {"B_over_Pstar_phi": list(bphi)},
    }
    loop_bindings = {
        "N": N,
        "charts": {"first_micro": loop_chart},
    }
    q4_chart = {
        "coordinate": coordinate,
        "physical_measure": jacobian,
        "physical_scale": pstar,
        "radial_background": {"V": list(background)},
        "radial": {"V_y": list(total_V_y), "slow_B_y": list(slow_B_y)},
        "physical": {"U_y": list(physical_U_y)},
    }
    q4_bindings = {"charts": {"first_micro": q4_chart}}

    q4 = {
        "source_family": family,
        # The two owners deliberately share this exact synthetic graph array.
        "exact_graph_nodes": graph,
        "modified_Q4_pressure_C3_bindings": q4_bindings,
    }
    loop = {
        "source_family": family,
        "exact_graph_nodes": graph,
        "prescribed_loop_C5_bindings": loop_bindings,
    }
    hashes = dict(RAW_SHA)
    flat_receipt = {
        "result": {
            "archive_result": {
                "status": "verified",
                "evidence": {
                    "predicate": predicate,
                    "coordinate": coordinate,
                    "collar": collar,
                    "lazy_source_rows": {"B_over_Pstar_phi": list(bphi)},
                    "common_N_node": N,
                },
                "scope": {
                    "flat_at_sc_over_two_conditional_on_sc_positive": True,
                    "sc_positive_independently_certified": False,
                },
            },
            "input_hashes": dict(hashes),
        }
    }
    physical_receipt = {
        "result": {
            "archive_result": {
                "status": "verified",
                "evidence": {
                    "source_normalized_rows": {"U_y": list(total_V_y)},
                    "physical_rows": {"U_y": list(physical_U_y)},
                    "Pstar_candidate_node": pstar,
                },
            },
            "input_hashes": dict(hashes),
        }
    }
    meta = {
        "coordinate": coordinate,
        "jacobian": jacobian,
        "N": N,
        "predicate": predicate,
        "collar": collar,
        "pstar": pstar,
        "invjac": invjac,
        "invN": invN,
        "bphi": bphi,
        "derivatives": derivatives,
        "slow_B_y": slow_B_y,
        "background": background,
        "scaled_slow_B": scaled_slow_B,
        "total_V_y": total_V_y,
        "physical_U_y": physical_U_y,
    }
    return q4, loop, flat_receipt, physical_receipt, meta


def _diagnostic(case, result):
    if not isinstance(result, dict):
        return {"case": case, "passed": False, "reason": "verifier did not return a mapping"}
    checks = result.get("checks", {})
    failed = []
    if isinstance(checks, dict):
        failed = [name for name, value in checks.items() if value is not True][:8]
    return {
        "case": case,
        "passed": result.get("passed") is True,
        "status": result.get("status"),
        "failed_checks": failed,
    }


def run_fixture(verify):
    """Run one positive and bounded fail-closed synthetic verification suite."""

    checks = {}
    diagnostics = []
    if not callable(verify):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"verify_callable": False},
            "diagnostics": [{"case": "input", "passed": False, "reason": "verify is not callable"}],
            "synthetic_only": True,
        }
    checks["verify_callable"] = True

    q4, loop, flat_receipt, physical_receipt, meta = _build_case()
    try:
        result = verify(q4, loop, flat_receipt, physical_receipt)
        checks["positive_verified"] = (
            isinstance(result, dict)
            and result.get("status") == "verified"
            and result.get("passed") is True
        )
        checks["positive_keeps_three_ordinary_rows"] = (
            isinstance(result, dict)
            and result.get("scope", {}).get("ordinary_Z_orders") == [0, 1, 2]
            and len(result.get("evidence", {}).get("slow_B_y_rows", ())) == 3
            and len(result.get("evidence", {}).get("total_V_y_rows", ())) == 3
        )
        checks["positive_is_conditional_on_sc_positive"] = (
            isinstance(result, dict)
            and result.get("scope", {}).get(
                "Q4_Uy_correction_zero_at_endpoint_conditional_on_sc_positive"
            ) is True
            and result.get("scope", {}).get("sc_positive_independently_certified") is False
        )
        checks["positive_does_not_claim_corrected_identity"] = (
            isinstance(result, dict)
            and result.get("scope", {}).get("corrected_source_identity_proved") is False
            and result.get("scope", {}).get("history_function_equality_claimed") is False
        )
        checks["positive_does_not_numericize_physical_Uy"] = (
            isinstance(result, dict)
            and result.get("scope", {}).get("numeric_values_computed") is False
            and all(isinstance(item, int) for item in meta["physical_U_y"])
        )
        diagnostics.append(_diagnostic("positive", result))
    except Exception as exc:
        checks["positive_verified"] = False
        diagnostics.append({
            "case": "positive",
            "passed": False,
            "exception": type(exc).__name__ + ": " + str(exc),
        })

    # Each case starts from a new graph and new receipt dictionaries.  The
    # mutations target only the synthetic handles/rows named by the case.
    negative_cases = (
        "wrong_predicate",
        "wrong_flat_value",
        "wrong_ordinary_Z_order",
        "wrong_invjac",
        "wrong_invN",
        "missing_derivative",
        "wrong_background",
        "wrong_Bphi",
        "wrong_receipt",
        "wrong_source_row",
        "changed_real_node_body",
    )
    for case in negative_cases:
        q4, loop, flat_receipt, physical_receipt, meta = _build_case()
        qnodes = q4["exact_graph_nodes"]
        qchart = q4["modified_Q4_pressure_C3_bindings"]["charts"]["first_micro"]
        lchart = loop["prescribed_loop_C5_bindings"]["charts"]["first_micro"]
        slow = qchart["radial"]["slow_B_y"]
        total = qchart["radial"]["V_y"]
        if case == "wrong_predicate":
            qnodes[meta["derivatives"][0]]["flat_predicate"] = meta["collar"]
        elif case == "wrong_flat_value":
            qnodes[meta["derivatives"][1]]["flat_value"] = 1
        elif case == "wrong_ordinary_Z_order":
            qnodes[meta["derivatives"][2]]["ordinary_slow_Z_derivative_order"] = 7
        elif case == "wrong_invjac":
            qnodes[slow[0]]["arguments"][0] = meta["invN"]
        elif case == "wrong_invN":
            qnodes[meta["scaled_slow_B"][0]]["arguments"][1] = meta["invjac"]
        elif case == "missing_derivative":
            qnodes[slow[1]]["arguments"][1] = len(qnodes) + 17
        elif case == "wrong_background":
            qnodes[total[0]]["arguments"][0] = meta["bphi"][0]
        elif case == "wrong_Bphi":
            qnodes[total[1]]["arguments"][2] = meta["background"][1]
        elif case == "wrong_receipt":
            flat_receipt["result"]["archive_result"]["status"] = "mismatch"
        elif case == "wrong_source_row":
            physical_receipt["result"]["archive_result"]["evidence"][
                "source_normalized_rows"
            ]["U_y"][2] = meta["physical_U_y"][2]
        elif case == "changed_real_node_body":
            qnodes[meta["invjac"]]["operation"] = "product"

        try:
            result = verify(q4, loop, flat_receipt, physical_receipt)
            checks["reject_" + case] = (
                isinstance(result, dict)
                and result.get("passed") is False
                and result.get("status") in ("mismatch", "unresolved")
            )
            diagnostics.append(_diagnostic(case, result))
        except Exception as exc:
            checks["reject_" + case] = False
            diagnostics.append({
                "case": case,
                "passed": False,
                "exception": type(exc).__name__ + ": " + str(exc),
            })

    return {
        "schema": SCHEMA,
        "passed": bool(checks) and all(value is True for value in checks.values()),
        "checks": checks,
        "diagnostics": diagnostics,
        "synthetic_only": True,
        "source_code_executed": False,
        "source_graph_walked": False,
        "source_values_evaluated": False,
        "corrected_source_identity_proved": False,
    }


__all__ = ["SCHEMA", "run_fixture"]
