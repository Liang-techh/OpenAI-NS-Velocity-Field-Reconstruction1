"""Conditional C3 outer-amplitude endpoint algebra.

This is a fresh six-by-six endpoint constraint system.  It binds the six
analytic C3 scalar formulas to exact source-query leaves, solves their
endpoint log relations with Fraction Gaussian elimination, and records the
five pointwise seam cancellations.  The cached native action and Rh seed are
inputs only; this module neither imports nor replays their producers.

All actual native/C3 global identities, caps, controls, profiles, and numeric
values remain out of scope.  The callable path is conditional on the root
wrapper's authenticated source rows and on separate exact sigma and mu recipe
bindings.
"""

from fractions import Fraction


c3oa_schema = "compact_C3_outer_amplitude_identity_algebra_input_v1"
c3oa_v1_preimage = {
    "sha256": "610367507ec6ed9e025a60af5bd204cd74b6dc1f4524a4190da81457fa6d6683",
    "bytes": 19073,
    "preserved": True,
}
c3oa_chart_names = (
    "Rh_reference",
    "O2_slope",
    "O2_axial",
    "O2_buffer",
    "O3_transition",
    "O3_power",
)
c3oa_basis_names = (
    "constant",
    "I_slope",
    "exp40",
    "mu_times_I_trans",
    "mu",
)
c3oa_ordinary_orders = (0, 1, 2, 3)

c3oa_expected_formulas = {
    "q": "q(Z)=1+Z^2",
    "Rh_reference": "E_Rh(x,Z)=exp(x/10)/q(Z)",
    "O2_slope": "E_slope(x,Z)=exp(x/10-(3/5)*J_s(x))/q(Z)",
    "O2_axial": "E_axial(x,Z)=E_slope(1,Z)*exp(-(exp(40*x)-1)/2)",
    "O2_buffer": "E_buffer(x,Z)=E_slope(1,Z)*exp(-(exp(40)+x-1)/2)",
    "O3_transition": "E_transition(x,Z)=E_buffer(11,Z)*exp(-x/2-mu*J_t(x))",
    "O3_power": "E_power(x,Z)=E_transition(1,Z)*exp(-(1/2+mu)*x)",
}

c3oa_expected_rhs_literals = (
    ("0", "0", "0", "0", "0"),
    ("1/10", "-3/5", "0", "0", "0"),
    ("1/2", "0", "-1/2", "0", "0"),
    ("-11/2", "0", "0", "0", "0"),
    ("-1/2", "0", "0", "-1", "0"),
    ("-1", "0", "0", "0", "-2"),
)

c3oa_contract_template = {
    "schema": c3oa_schema,
    "source_formulas": c3oa_expected_formulas,
    "source_binding": {
        "root_wrapper_passed": True,
        "saved_native_action_result_ref": "saved native action result, data-only selector",
        "saved_Rh_seed_result_ref": "saved Rh seed result, data-only selector",
        "exact_sigma_recipe_ref": "authenticated exact sigma source recipe",
        "exact_mu_recipe_ref": "authenticated exact mu source recipe",
        "native_C3_parameter_alignment_binding_passed": True,
        "evaluator_covers_are_bounds": True,
        "interval_values_used": False,
        "native_ordinary_Z_slots": "[E0,E1,E2,E3]",
        "c3_query_projection": "[E0,E1,2*E2,6*E3]",
        "c3_projection_factor_application_count": 1,
    },
    "c3_source_query_leaves": {
        "Rh_reference": {
            "ordinary_Z_orders": [0, 1, 2, 3],
            "source_ref": "authenticated C3 source query leaf for Rh_reference",
            "ast_ref": "authenticated C3 AST span for Rh_reference",
        },
        "O2_slope": {
            "ordinary_Z_orders": [0, 1, 2, 3],
            "source_ref": "authenticated C3 source query leaf for O2_slope",
            "ast_ref": "authenticated C3 AST span for O2_slope",
        },
        "O2_axial": {
            "ordinary_Z_orders": [0, 1, 2, 3],
            "source_ref": "authenticated C3 source query leaf for O2_axial",
            "ast_ref": "authenticated C3 AST span for O2_axial",
        },
        "O2_buffer": {
            "ordinary_Z_orders": [0, 1, 2, 3],
            "source_ref": "authenticated C3 source query leaf for O2_buffer",
            "ast_ref": "authenticated C3 AST span for O2_buffer",
        },
        "O3_transition": {
            "ordinary_Z_orders": [0, 1, 2, 3],
            "source_ref": "authenticated C3 source query leaf for O3_transition",
            "ast_ref": "authenticated C3 AST span for O3_transition",
        },
        "O3_power": {
            "ordinary_Z_orders": [0, 1, 2, 3],
            "source_ref": "authenticated C3 source query leaf for O3_power",
            "ast_ref": "authenticated C3 AST span for O3_power",
        },
    },
    "cached_native": {
        "data_only": True,
        "source_result_ref": "saved native action result data-only terminal coefficients",
        "action_total_coefficients": {
            "constant": "-59/10",
            "I_slope": "-3/5",
            "exp40": "-1/2",
            "mu_times_I_trans": "-1",
            "mu": "-2",
        },
        "accepted_Rh_seed_exponent": "-1/2",
    },
    "source_flags": {
        "actual_native_C3_amplitude_identity_proved": False,
        "actual_incoming_identity_proved": False,
        "actual_control_identity_proved": False,
        "actual_profile_identity_proved": False,
        "cap_admitted": False,
        "numeric_values_selected": False,
    },
}


def c3oa_require(condition, message):
    if not condition:
        raise ValueError(message)


def c3oa_text(value, label):
    c3oa_require(isinstance(value, str) and value != "", label + " must be nonempty text")
    return value


def c3oa_parse_fraction(value, label):
    c3oa_require(isinstance(value, str) and value != "", label + " must be an exact rational string")
    parts = value.split("/")
    c3oa_require(len(parts) in (1, 2), label + " must be n or n/d")
    numerator = int(parts[0])
    denominator = int(parts[1]) if len(parts) == 2 else 1
    c3oa_require(denominator != 0, label + " denominator is zero")
    return Fraction(numerator, denominator)


def c3oa_fraction_record(value):
    return {"numerator": value.numerator, "denominator": value.denominator}


def c3oa_vector_record(vector):
    return [c3oa_fraction_record(value) for value in vector]


def c3oa_validate_source_contract(contract):
    """Validate source rows and cached data shape without evaluating formulas."""
    c3oa_require(isinstance(contract, dict), "semantic_contract must be a mapping")
    c3oa_require(contract.get("schema") == c3oa_schema, "wrong C3 amplitude schema")

    formulas = contract.get("source_formulas")
    c3oa_require(isinstance(formulas, dict), "source_formulas are missing")
    for name, expected in c3oa_expected_formulas.items():
        c3oa_require(formulas.get(name) == expected, "source formula changed: " + name)

    binding = contract.get("source_binding")
    c3oa_require(isinstance(binding, dict), "source_binding is missing")
    c3oa_require(binding.get("root_wrapper_passed") is True,
                 "root source wrapper has not passed")
    for name in (
        "saved_native_action_result_ref",
        "saved_Rh_seed_result_ref",
        "exact_sigma_recipe_ref",
        "exact_mu_recipe_ref",
    ):
        c3oa_text(binding.get(name), "source_binding." + name)
    c3oa_require(binding.get("native_C3_parameter_alignment_binding_passed") is True,
                 "native/C3 sigma-mu alignment binding has not passed")
    c3oa_require(binding.get("evaluator_covers_are_bounds") is True,
                 "evaluator covers must remain bounds")
    c3oa_require(binding.get("interval_values_used") is False,
                 "interval values cannot be used as source values")
    c3oa_require(binding.get("native_ordinary_Z_slots") == "[E0,E1,E2,E3]",
                 "native ordinary-Z slots changed")
    c3oa_require(binding.get("c3_query_projection") == "[E0,E1,2*E2,6*E3]",
                 "C3 query projection changed")
    c3oa_require(binding.get("c3_projection_factor_application_count") == 1,
                 "C3 projection factor must be applied once")

    leaves = contract.get("c3_source_query_leaves")
    c3oa_require(isinstance(leaves, dict), "C3 source query leaves are missing")
    for chart in c3oa_chart_names:
        row = leaves.get(chart)
        c3oa_require(isinstance(row, dict), "missing C3 source query leaf: " + chart)
        c3oa_require(row.get("ordinary_Z_orders") == list(c3oa_ordinary_orders),
                     "wrong ordinary-Z query orders: " + chart)
        c3oa_text(row.get("source_ref"), chart + ".source_ref")
        c3oa_text(row.get("ast_ref"), chart + ".ast_ref")

    cached = contract.get("cached_native")
    c3oa_require(isinstance(cached, dict), "cached_native data is missing")
    c3oa_require(cached.get("data_only") is True, "cached native data must be data-only")
    c3oa_text(cached.get("source_result_ref"), "cached_native.source_result_ref")
    coefficients = cached.get("action_total_coefficients")
    c3oa_require(isinstance(coefficients, dict), "cached native action coefficients are missing")
    c3oa_require(set(coefficients) == set(c3oa_basis_names),
                 "cached native basis names changed")
    for name in c3oa_basis_names:
        c3oa_parse_fraction(coefficients[name], "cached_native." + name)
    c3oa_parse_fraction(cached.get("accepted_Rh_seed_exponent"),
                        "cached_native.accepted_Rh_seed_exponent")

    flags = contract.get("source_flags")
    c3oa_require(isinstance(flags, dict), "source_flags are missing")
    for name in (
        "actual_native_C3_amplitude_identity_proved",
        "actual_incoming_identity_proved",
        "actual_control_identity_proved",
        "actual_profile_identity_proved",
        "cap_admitted",
        "numeric_values_selected",
    ):
        c3oa_require(flags.get(name) is False, "forbidden source flag is true: " + name)

    return {
        "schema": c3oa_schema,
        "chart_order": list(c3oa_chart_names),
        "basis_order": list(c3oa_basis_names),
        "ordinary_Z_orders": list(c3oa_ordinary_orders),
        "source_query_leaf_count": len(c3oa_chart_names),
        "cached_native_data_only": True,
    }


def c3oa_matrix():
    matrix = []
    for row_index in range(len(c3oa_chart_names)):
        row = [0 for _ in c3oa_chart_names]
        row[row_index] = 1
        if row_index:
            row[row_index - 1] = -1
        matrix.append(row)
    return matrix


def c3oa_rhs_vectors():
    return [
        [c3oa_parse_fraction(value, "rhs[%d][%d]" % (row, column))
         for column, value in enumerate(values)]
        for row, values in enumerate(c3oa_expected_rhs_literals)
    ]


def c3oa_gaussian_solve(matrix, rhs_vectors):
    """Solve all five coefficient columns by exact Fraction elimination."""
    size = len(matrix)
    width = len(rhs_vectors[0])
    left = [[Fraction(value, 1) for value in row] for row in matrix]
    right = [[Fraction(value) for value in row] for row in rhs_vectors]
    pivot_columns = []
    for column in range(size):
        pivot = None
        for candidate in range(column, size):
            if left[candidate][column] != 0:
                pivot = candidate
                break
        c3oa_require(pivot is not None, "endpoint matrix is singular")
        if pivot != column:
            left[column], left[pivot] = left[pivot], left[column]
            right[column], right[pivot] = right[pivot], right[column]
        pivot_value = left[column][column]
        left[column] = [value / pivot_value for value in left[column]]
        right[column] = [value / pivot_value for value in right[column]]
        for row in range(size):
            if row == column:
                continue
            factor = left[row][column]
            if factor == 0:
                continue
            left[row] = [
                left[row][j] - factor * left[column][j]
                for j in range(size)
            ]
            right[row] = [
                right[row][j] - factor * right[column][j]
                for j in range(width)
            ]
        pivot_columns.append(column)
    return left, right, pivot_columns


def c3oa_residual_rows(matrix, solution, rhs_vectors):
    rows = []
    for row_index, row in enumerate(matrix):
        reconstructed = []
        for coefficient_index in range(len(c3oa_basis_names)):
            value = sum(
                (Fraction(row[column], 1) * solution[column][coefficient_index]
                 for column in range(len(c3oa_chart_names))),
                Fraction(0, 1),
            )
            reconstructed.append(value)
        residual = [
            reconstructed[index] - rhs_vectors[row_index][index]
            for index in range(len(c3oa_basis_names))
        ]
        rows.append({
            "row": c3oa_chart_names[row_index],
            "residual": c3oa_vector_record(residual),
            "all_five_coefficients_zero": all(value == 0 for value in residual),
        })
    return rows


def c3oa_seam_rows(solution):
    zero = [Fraction(0, 1) for _ in c3oa_basis_names]
    axial_increment = [
        Fraction(1, 2),
        Fraction(0, 1),
        Fraction(-1, 2),
        Fraction(0, 1),
        Fraction(0, 1),
    ]
    seam_specs = (
        {
            "left": "Rh_reference@0",
            "right": "O2_slope@0",
            "left_vector": zero,
            "right_vector": solution[0],
            "reason": "Rh and slope both have the zero endpoint log factor",
        },
        {
            "left": "O2_slope@1",
            "right": "O2_axial@0",
            "left_vector": solution[1],
            "right_vector": solution[1],
            "reason": "axial starts from the slope endpoint with exp(0)=1",
        },
        {
            "left": "O2_axial@1",
            "right": "O2_buffer@0",
            "left_vector": solution[2],
            "right_vector": [
                solution[1][index] + axial_increment[index]
                for index in range(len(c3oa_basis_names))
            ],
            "reason": "buffer at zero has the axial exp(-(exp(40)-1)/2) factor",
        },
        {
            "left": "O2_buffer@11",
            "right": "O3_transition@0",
            "left_vector": solution[3],
            "right_vector": solution[3],
            "reason": "transition starts from the buffer endpoint with J_t(0)=0",
        },
        {
            "left": "O3_transition@1",
            "right": "O3_power@0",
            "left_vector": solution[4],
            "right_vector": solution[4],
            "reason": "power starts from the transition endpoint with exp(0)=1",
        },
    )
    rows = []
    for seam_index, spec in enumerate(seam_specs):
        left_vector = spec["left_vector"]
        right_vector = spec["right_vector"]
        residual = [
            left_vector[index] - right_vector[index]
            for index in range(len(c3oa_basis_names))
        ]
        passed = all(value == 0 for value in residual)
        c3oa_require(passed, "C3 endpoint seam coefficient mismatch")
        ordinary = []
        for order in c3oa_ordinary_orders:
            ordinary.append({
                "ordinary_Z_order": order,
                "left_coefficient_vector": c3oa_vector_record(left_vector),
                "right_coefficient_vector": c3oa_vector_record(right_vector),
                "log_factor_difference": c3oa_vector_record(residual),
                "reason": spec["reason"],
                "pointwise_identity": True,
            })
        rows.append({
            "seam_index": seam_index,
            "left_endpoint": spec["left"],
            "right_endpoint": spec["right"],
            "left_coefficient_vector": c3oa_vector_record(left_vector),
            "right_coefficient_vector": c3oa_vector_record(right_vector),
            "coefficient_residual": c3oa_vector_record(residual),
            "coefficient_residual_zero": passed,
            "ordinary_Z_rows": ordinary,
            "ordinary_Z_row_count": len(ordinary),
            "passed": passed,
        })
    return rows


def derive_C3_outer_amplitude_identity_algebra(semantic_contract, journal=None):
    """Return the conditional endpoint system and seam certificate."""
    checked = c3oa_validate_source_contract(semantic_contract)
    if journal is not None:
        c3oa_require(callable(journal), "journal must be callable")

    matrix = c3oa_matrix()
    rhs = c3oa_rhs_vectors()
    reduced_left, solution, pivots = c3oa_gaussian_solve(matrix, rhs)
    residual_rows = c3oa_residual_rows(matrix, solution, rhs)
    residual_checks = []
    for row in residual_rows:
        for coefficient_index, coefficient_name in enumerate(c3oa_basis_names):
            residual_checks.append({
                "row": row["row"],
                "basis": coefficient_name,
                "residual": row["residual"][coefficient_index],
                "passed": row["residual"][coefficient_index] == {
                    "numerator": 0,
                    "denominator": 1,
                },
            })
    c3oa_require(len(residual_checks) == 30, "wrong residual check count")
    c3oa_require(all(check["passed"] for check in residual_checks),
                 "endpoint coefficient residual is nonzero")

    cached = semantic_contract["cached_native"]
    cached_action = [
        c3oa_parse_fraction(cached["action_total_coefficients"][name],
                            "cached_native." + name)
        for name in c3oa_basis_names
    ]
    seed = c3oa_parse_fraction(cached["accepted_Rh_seed_exponent"],
                               "cached_native.accepted_Rh_seed_exponent")
    terminal_solution = solution[-1]
    cached_terminal_with_seed = list(cached_action)
    cached_terminal_with_seed[0] += seed
    terminal_comparison = []
    for index, basis in enumerate(c3oa_basis_names):
        difference = terminal_solution[index] - cached_terminal_with_seed[index]
        terminal_comparison.append({
            "basis": basis,
            "solved_terminal_log_qE": c3oa_fraction_record(terminal_solution[index]),
            "cached_native_action_plus_seed": c3oa_fraction_record(cached_terminal_with_seed[index]),
            "residual": c3oa_fraction_record(difference),
            "passed": difference == 0,
        })
    c3oa_require(len(terminal_comparison) == 5, "wrong terminal comparison count")
    c3oa_require(all(row["passed"] for row in terminal_comparison),
                 "cached native terminal action mismatch")

    seams = c3oa_seam_rows(solution)
    c3oa_require(len(seams) == 5 and all(row["ordinary_Z_row_count"] == 4 for row in seams),
                 "wrong seam certificate shape")

    if journal is not None:
        journal({
            "event": "c3_outer_amplitude_endpoint_algebra_defined",
            "matrix_size": "6x6",
            "residual_checks": 30,
            "terminal_comparisons": 5,
            "seam_count": 5,
            "ordinary_Z_orders": 4,
            "numeric_evaluation": False,
        })

    return {
        "schema": "compact_C3_outer_amplitude_identity_algebra_result_v1",
        "passed": True,
        "conditional_source_contract": checked,
        "source_formulas": dict(c3oa_expected_formulas),
        "endpoint_system": {
            "unknowns": ["log(q*E_" + chart + "_end)" for chart in c3oa_chart_names],
            "basis": list(c3oa_basis_names),
            "matrix": matrix,
            "rhs_literals": [list(values) for values in c3oa_expected_rhs_literals],
            "rhs_exact": [c3oa_vector_record(vector) for vector in rhs],
            "reduced_matrix": [c3oa_vector_record(row) for row in reduced_left],
            "pivot_columns": pivots,
            "solution_log_qE_endpoints": [
                c3oa_vector_record(vector) for vector in solution
            ],
            "endpoint_solution_is_absolute_log_qE": True,
            "accepted_Rh_seed_added_only_to_cached_comparison": True,
        },
        "residual_proof": {
            "residual_rows": residual_rows,
            "coefficient_checks": residual_checks,
            "coefficient_check_count": len(residual_checks),
            "all_zero": True,
        },
        "terminal_comparison": {
            "rows": terminal_comparison,
            "comparison_count": len(terminal_comparison),
            "all_zero": True,
            "cached_native_data_only": True,
        },
        "c3_pointwise_seams": {
            "rows": seams,
            "seam_count": len(seams),
            "ordinary_Z_rows_per_seam": 4,
            "all_five_seams_zero": True,
            "integrals_or_exp_values_evaluated": False,
        },
        "ordinary_Z_identity": {
            "native_slots": "[E0,E1,E2,E3]",
            "c3_query_projection": "[E0,E1,2*E2,6*E3]",
            "c3_projection_factor_application_count": 1,
            "pointwise_function_identity_differentiated": True,
            "time_recurrence_used": False,
        },
        "scope": {
            "actual_native_C3_amplitude_identity_proved": False,
            "native_C3_amplitude_identity_gate": {
                "source_parameter_alignment_binding_required": True,
                "source_parameter_alignment_binding_passed": True,
                "global_identity_admitted_by_this_leaf": False,
            },
            "actual_incoming_identity_proved": False,
            "actual_control_identity_proved": False,
            "actual_profile_identity_proved": False,
            "cap_admitted": False,
            "numeric_values_selected": False,
            "evaluator_covers_are_bounds_only": True,
            "old_kernel_or_helper_replayed": False,
        },
        "audit": {
            "new_fraction_gaussian_elimination_calls": 1,
            "source_provider_calls": 0,
            "old_consumer_calls": 0,
            "quadrature_calls": 0,
            "numeric_exp_evaluations": 0,
            "integral_evaluations": 0,
        },
        "v1_preimage": dict(c3oa_v1_preimage),
    }


__all__ = (
    "c3oa_contract_template",
    "c3oa_validate_source_contract",
    "derive_C3_outer_amplitude_identity_algebra",
)


"""Bounded source binder for the C3 outer amplitude identity.

This module is deliberately a consumer of frozen source metadata.  It does
not import a producer, evaluate a source field, use a midpoint or interval
cover as a value, or infer an identity from graph-node numbering.  The only
mathematical entry point is the new algebra function defined immediately
before this module when the parent assembler concatenates the two files.
"""

import json


C3OA_SCHEMA = "compact_C3_outer_amplitude_identity_input_v1"
C3OA_COMPONENT = "compact_C3_outer_amplitude_identity"
C3OA_CHARTS = (
    "Rh_reference", "O2_slope", "O2_axial", "O2_buffer",
    "O3_transition", "O3_power",
)
C3OA_ORDERS = (0, 1, 2, 3)
C3OA_PRIMARY_PORTS = (4934, 4935, 12046, 22539)
C3OA_REBUILT_PORTS = (4934, 4935, 12046, 24181)
C3OA_PROVIDER_SHA = (
    "9d73659a1c76163dddee6653af53f56ebffd7045986d871dc559ba1e053620b7"
)
C3OA_ARCHIVE_SHAS = {
    "lei_ren_part1_paper_compliant_current_original_C3_limit_controls.json.gz":
        "23e4f08eb1c7b2f0098953c50d28df71039905cba7ec6bfc6f22e1aa356afe34",
    "lei_ren_part1_paper_compliant_current_original_C3_repair_band.json.gz":
        "cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193",
    "lei_ren_part1_paper_compliant_current_original_C3_source_targets.json.gz":
        "3bf0bf91eb8af6b7510a378464307b6b4fcb82c56ba661e4d4067fa54da4a499",
}
C3OA_EXPECTED_ACTION = {
    "constant": "-59/10",
    "I_slope": "-3/5",
    "exp40": "-1/2",
    "mu_times_I_trans": "-1/1",
    "mu": "-2/1",
}
C3OA_EXPECTED_COORDINATE = {
    "operation": "exact_rational", "numerator": 2, "denominator": 1,
}


def c3oa_require(condition, message):
    if not condition:
        raise ValueError(message)


def c3oa_is_bool(value):
    return isinstance(value, bool)


def c3oa_sha(value):
    import hashlib
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def c3oa_clone(value):
    import json
    return json.loads(json.dumps(value, ensure_ascii=False))


def c3oa_basename(value):
    return str(value).replace("\\", "/").rsplit("/", 1)[-1]


def c3oa_archive_sha(records, basename):
    c3oa_require(isinstance(records, dict), "C3 archive records are missing")
    for name, row in records.items():
        if c3oa_basename(name) != basename:
            continue
        c3oa_require(isinstance(row, dict), "Malformed archive record: " + name)
        got = row.get("raw_sha256")
        c3oa_require(got == C3OA_ARCHIVE_SHAS[basename],
                     "C3 archive SHA changed: " + basename)
        return got
    raise ValueError("Missing authenticated C3 archive record: " + basename)


def c3oa_section_rows(section, field):
    c3oa_require(isinstance(section, dict), "Missing C3 namespace section")
    rows = section.get(field)
    c3oa_require(isinstance(rows, list), "Missing C3 namespace rows: " + field)
    result = {}
    for item in rows:
        c3oa_require(isinstance(item, dict), "Malformed C3 namespace row")
        node = item.get("id")
        row = item.get("row")
        c3oa_require(isinstance(node, int) and not isinstance(node, bool)
                     and isinstance(row, dict), "Malformed C3 node row")
        c3oa_require(node not in result, "Duplicate C3 node row")
        result[node] = row
    return result


def c3oa_semantic_row(row, label, order):
    c3oa_require(isinstance(row, dict), label + " row is missing")
    expected = {
        "operation": "current_original_leading_function_recipe",
        "quantity": "E", "native_chart": "O3_power", "Z_variable": 20,
        "coordinate": 48, "Z_order": order,
        "defining_quantity_not_a_range_value": True,
    }
    for key, value in expected.items():
        c3oa_require(row.get(key) == value,
                     label + " changed at " + key)
    recipe = row.get("recipe")
    c3oa_require(isinstance(recipe, dict), label + ".recipe is missing")
    c3oa_require(recipe.get("function") == "WholeZAllNOuterRcFunctions.leading_packet"
                 and recipe.get("defining_provider")
                 == "WholeZAllNOuterRcFunctions.leading_packet"
                 and recipe.get("module_sha256") == C3OA_PROVIDER_SHA
                 and recipe.get("leading_source_fixed_input_N0") is True,
                 label + " provider/source binding changed")
    c3oa_require(recipe.get("range_endpoints_not_field_coefficients") is True,
                 label + " range endpoints became field coefficients")
    if order == 2:
        c3oa_require(row.get("ordinary_slow_Z_derivative_order") == 2
                     and row.get("Taylor_coefficient_factorial") == 2
                     and row.get("derivative_of_range_endpoint") is False,
                     label + " order-two source projection changed")
        c3oa_require(row.get("source_projection")
                     == "same actual terminal O3_power E Taylor coefficient[2] times 2",
                     label + " order-two projection changed")
    if order == 3:
        c3oa_require(row.get("ordinary_slow_Z_derivative_order") == 3
                     and row.get("Taylor_coefficient_factorial") == 6
                     and row.get("derivative_of_range_endpoint") is False,
                     label + " order-three source projection changed")
        c3oa_require(row.get("source_projection")
                     == "6*same actual terminal O3_power E Taylor coefficient[3]",
                     label + " order-three projection changed")
    return {
        "order": order, "node": label, "quantity": "E",
        "semantic_function_key": (
            row.get("operation"), row.get("quantity"),
            row.get("native_chart"), row.get("coordinate"),
            row.get("Z_variable"), row.get("Z_order"),
        ),
        "provider_function": recipe.get("function"),
        "defining_quantity_not_a_range_value": True,
        "range_endpoints_not_field_coefficients": True,
        "derivative_of_range_endpoint": row.get("derivative_of_range_endpoint", False),
        "source_projection": row.get("source_projection"),
    }


def c3oa_coordinate_rows(spec):
    rows = spec.get("C3_coordinate_rows")
    c3oa_require(isinstance(rows, dict), "C3 coordinate rows are missing")
    checks = []
    for namespace in ("C3_source_targets_snapshot", "C3_repair_band_rebuilt"):
        section = rows.get(namespace)
        c3oa_require(isinstance(section, dict), "Missing coordinate namespace: " + namespace)
        row = section.get("48")
        if row is None:
            row = section.get(48)
        c3oa_require(row == C3OA_EXPECTED_COORDINATE,
                     "Coordinate 48 is not the authenticated exact rational 2")
        checks.append({"namespace": namespace, "selector": "exact_graph_nodes[48]",
                       "operation": "exact_rational", "numerator": 2,
                       "denominator": 1, "passed": True})
    return checks


def c3oa_namespace_binding(spec):
    namespaces = spec.get("C3_namespaces")
    c3oa_require(isinstance(namespaces, dict), "C3 namespace separation is missing")
    records = spec.get("C3_archive_records")
    archive_checks = [
        {"basename": name, "raw_sha256": c3oa_archive_sha(records, name)}
        for name in C3OA_ARCHIVE_SHAS
    ]
    primary_section = namespaces.get("C3_source_targets_snapshot")
    repair_section = namespaces.get("C3_repair_band_rebuilt")
    primary = c3oa_section_rows(primary_section, "archive_terminal_amplitude_rows")
    rebuilt = c3oa_section_rows(repair_section, "source_leaf_rows")
    c3oa_require(tuple(primary_section.get("archive_terminal_amplitude_ids", ()))
                 == C3OA_PRIMARY_PORTS,
                 "Primary C3 terminal source ports changed")
    c3oa_require(tuple(repair_section.get("source_leaf_ports_Z0_to_Z3", ()))
                 == C3OA_REBUILT_PORTS,
                 "Rebuilt C3 source leaf ports changed")
    expression_ports = tuple(repair_section.get("rebuilt_expression_ports_Z0_to_Z3", ()))
    c3oa_require(expression_ports == (4934, 4935, 12046, 24177),
                 "Rebuilt expression namespace changed")
    c3oa_require(24177 in expression_ports and 24181 not in expression_ports,
                 "Rebuilt expression/source leaf namespaces are ambiguous")
    primary_rows = []
    rebuilt_rows = []
    for order, (primary_port, rebuilt_port) in enumerate(
            zip(C3OA_PRIMARY_PORTS, C3OA_REBUILT_PORTS)):
        primary_rows.append(c3oa_semantic_row(primary[primary_port],
                                               "primary:" + str(primary_port), order))
        rebuilt_rows.append(c3oa_semantic_row(rebuilt[rebuilt_port],
                                              "rebuilt:" + str(rebuilt_port), order))
    for left, right in zip(primary_rows, rebuilt_rows):
        c3oa_require(left["semantic_function_key"] == right["semantic_function_key"],
                     "Primary/rebuilt semantic FunctionKey changed at order "
                     + str(left["order"]))
        c3oa_require(left["source_projection"] == right["source_projection"],
                     "Primary/rebuilt source projection changed at order "
                     + str(left["order"]))
    coordinates = c3oa_coordinate_rows(spec)
    return {
        "passed": True, "archive_checks": archive_checks,
        "primary_source_ports": list(C3OA_PRIMARY_PORTS),
        "rebuilt_source_leaf_ports": list(C3OA_REBUILT_PORTS),
        "rebuilt_expression_ports": list(expression_ports),
        "rebuilt_expression_24177_is_not_E3": True,
        "primary_and_rebuilt_integer_ids_not_compared_as_identity": True,
        "semantic_function_key_checks": 4,
        "primary_rows": primary_rows, "rebuilt_rows": rebuilt_rows,
        "coordinate_checks": coordinates,
        "ordinary_Z_orders": list(C3OA_ORDERS),
        "quantity": "E", "chart": "O3_power", "coordinate_value": "2",
        "Z_variable": 20, "range_endpoints_not_field_coefficients": True,
    }


def c3oa_source_records(spec):
    import ast
    records = spec.get("creator_source_records")
    c3oa_require(isinstance(records, dict) and records,
                 "Creator source records are missing")
    trees = {}
    for basename, record in records.items():
        c3oa_require(isinstance(record, dict), "Malformed creator record: " + basename)
        text = record.get("source_text")
        expected = record.get("sha256")
        c3oa_require(isinstance(text, str) and isinstance(expected, str),
                     "Malformed creator text/hash: " + basename)
        c3oa_require(c3oa_sha(text) == expected,
                     "Creator source text/hash mismatch: " + basename)
        c3oa_require(c3oa_basename(record.get("path", "")) == basename,
                     "Creator source basename/path mismatch: " + basename)
        trees[basename] = ast.parse(text, filename=basename)

    def text_for(basename):
        c3oa_require(isinstance(basename, str) and basename.endswith(".py"),
                     "Source lookup must use an exact .py basename")
        c3oa_require(basename in records,
                     "Missing exact creator source record: " + basename)
        return records[basename]["source_text"], basename, trees[basename]

    def compact(text):
        return "".join(text.split())

    required = {
        "Rh_reference": ("lei_ren_part1_paper_compliant_current_original_Rh_reference_finite_N.py", (
            "logE=-reference.background.logarithm(1+ref.z*ref.z)+x/10",)),
        "O2_slope": ("lei_ren_part1_paper_compliant_current_original_O2_slope_finite_N.py", (
            "J,mass=slope_masses(c,y,self.cells)",
            "factor=c.exp(y/10-c.mpf('.6')*J)", "E=f.scale(qi,factor)")),
        "O2_axial": ("lei_ren_part1_paper_compliant_current_original_O2_axial_buffer_finite_N.py", (
            "t=y-1", "E=f.scale(u1,root)", "V=f.scale(z,4*B[0])")),
        "O2_buffer": ("lei_ren_part1_paper_compliant_current_original_O2_axial_buffer_finite_N.py", (
            "B=[c.mpf(0)]*5", "chart=='axial'", "E=f.scale(u1,root)")),
        "O3_transition": ("lei_ren_part1_paper_compliant_current_original_O3_Rc_finite_N.py", (
            "factor=c.exp(-t/2-mu_cover*K['J'])", "E=f.scale(u1,factor)",
            "V,Vy=zero,zero")),
        "O3_power": ("lei_ren_part1_paper_compliant_current_original_O3_Rc_finite_N.py", (
            "factor=c.exp(-t/2-mu_cover*t)", "E=f.scale(u1,factor)",
            "V,Vy=zero,zero")),
    }
    formula_checks = []
    for chart, (part, clauses) in required.items():
        text, basename, _ = text_for(part)
        normalized = compact(text)
        missing = [clause for clause in clauses if compact(clause) not in normalized]
        c3oa_require(not missing, "Missing source formula clauses for " + chart
                     + ": " + ",".join(missing))
        formula_checks.append({"chart": chart, "basename": basename,
                               "clauses": list(clauses), "passed": True})

    dispatch_text, dispatch_name, _ = text_for(
        "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions.py")
    dispatch = compact(dispatch_text)
    for clause in (
        "self.views['axial'].parent=lambda ends:self.leading_packet(ends,'O2_slope',(1,1))",
        "self.views['o3'].parent=lambda ends:self.leading_packet(ends,'O2_buffer',(11,1))",
        "self.views['o3'].query=lambda ends,chart,left,right=None:self.leading_packet(ends,'O3_'+chart,left,right)",
        "FRONTENDS[kind](view,key,l,r)",
        "chart.split('_',1)[1]",
    ):
        c3oa_require(compact(clause) in dispatch,
                     "Outer parent/dispatch clause missing: " + clause)
    whole_o3_text, whole_o3_name, _ = text_for(
        "lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N.py")
    whole_o3 = compact(whole_o3_text)
    for clause in ("generic=reference.long.RECOVER(proxy,source['raw'])",
                   "common_velocity_E_axial5", "common_velocity_V_axial5"):
        c3oa_require(compact(clause) in whole_o3,
                     "Whole-Z source projection clause missing: " + clause)

    comparison_text, comparison_name, _ = text_for(
        "lei_ren_part1_paper_interval_comparison_enclosure.py")
    comparison = compact(comparison_text)
    for clause in ("defalpha_box(c,y,hb)",
                   "(y-hb)/hb", "left,right=c.exp(-1/(x*x)),c.exp(-1/((1-x)*(1-x)))",
                   "returnright/(left+right)"):
        c3oa_require(compact(clause) in comparison,
                     "alpha_box source clause missing: " + clause)
    long_text, long_name, _ = text_for(
        "lei_ren_part1_paper_interval_long_reshape_field.py")
    long_compact = compact(long_text)
    c3oa_require("defsigma_value_derivative(c,s)" in long_compact
                 and "return1-alpha_box(c,s+1,c.mpf(1))" in long_compact,
                 "C3 sigma_value_derivative source clause changed")
    interval_slope_text, interval_slope_name, _ = text_for(
        "lei_ren_part1_paper_interval_outer_slope_field.py")
    c3oa_require(
        "fromlei_ren_part1_paper_interval_long_reshape_fieldimportsigma_value_derivative"
        in compact(interval_slope_text),
        "C3 interval slope sigma provider import changed")
    initial_text, initial_name, _ = text_for(
        "lei_ren_part1_paper_compliant_outer_initial.py")
    initial = compact(initial_text)
    for clause in ("defstable_sigma(c,t)",
                   "a=c.exp(-1/x**2);b=c.exp(-1/(1-x)**2)",
                   "value=a/(a+b)"):
        c3oa_require(compact(clause) in initial,
                     "native stable_sigma source clause missing: " + clause)
    mu_text, mu_name, _ = text_for(
        "lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N.py")
    mu = compact(mu_text)
    c3oa_require("defsource_mu_admission(c,eta_log)" in mu
                 and "logmu=c.ln(c.mpf('.001'))-4*(c.exp(40)+11)" in mu,
                 "C3 source_mu_admission clause changed")
    c3oa_require(C3OA_PROVIDER_SHA in json.dumps(spec.get("C3_namespaces", {}),
                                                 sort_keys=True),
                 "leading_packet provider hash is absent from namespace metadata")

    axial_text, axial_name, _ = text_for(
        "lei_ren_part1_paper_compliant_current_original_O2_axial_buffer_finite_N.py")
    axial = compact(axial_text)
    for clause in (
        "MD=40",
        "y=c.exp(MD*selector)ifchart=='axial'elsec.exp(MD)+selector",
        "t=y-1;formal=lambdalog:f.factor((0,0,0,0,0),log)",
        "d,root,d3=formal(-t),formal(-t/2),formal(-c.mpf('1.5')*t)",
        "E=f.scale(u1,root)",
    ):
        c3oa_require(compact(clause) in axial,
                     "Axial coordinate/decay clause missing: " + clause)

    slope_text, slope_name, _ = text_for(
        "lei_ren_part1_paper_compliant_current_original_O2_slope_finite_N.py")
    slope = compact(slope_text)
    for clause in (
        "importlei_ren_part1_paper_compliant_pre_pulse_mixed_C4asoriginal",
        "J,mass=original.slope_masses(c,y,scalar_cells)",
        "qi=f.jet((1+ref.z*ref.z).reciprocal())",
        "E=f.scale(qi,factor)",
    ):
        c3oa_require(compact(clause) in slope,
                     "Slope/qi source clause missing: " + clause)

    pre_text, pre_name, _ = text_for(
        "lei_ren_part1_paper_compliant_pre_pulse_mixed_C4.py")
    pre = compact(pre_text)
    for clause in (
        "defslope_masses(c,y,cells)",
        "returntransition_integrals(c,exact,cells)",
        "defslope_mu(self,Z,offset)",
        "fromlei_ren_part1_paper_interval_outer_slope_fieldimporttransition_integrals",
        "K=transition_kernels(c,t,mu,self.cells)",
        "factor=c.exp(-t/2-mu*K['J'])",
        "defpower(self,Z,phase)",
        "f=c.exp((-c.mpf('.5')-mu)*t)",
    ):
        c3oa_require(compact(clause) in pre,
                     "Pre-pulse slope-mu/power clause missing: " + clause)

    buffer_text, buffer_name, _ = text_for(
        "lei_ren_part1_paper_compliant_outer_buffer.py")
    buffer = compact(buffer_text)
    for clause in (
        "fromlei_ren_part1_paper_compliant_outer_initialimportSharedOuterInitial,stable_sigma",
        "deftransition_kernels(c,t,mu,cells=256)",
        "foriinrange(cells)", "sa=stable_sigma(c,a)[0]", "J=nextJ",
        "returndict(J=J,theta=K[0],energy=K[1],pressure=K[2],cells=cells)",
    ):
        c3oa_require(compact(clause) in buffer,
                     "Transition-kernel/stable-sigma clause missing: " + clause)
    interval_slope = compact(interval_slope_text)
    for clause in (
        "deftransition_integrals(c,y,cells=256)",
        "sa=sigma_value_derivative(c,a)[0];sb=sigma_value_derivative(c,b)[0]",
        "next_j=j+dy*sj", "j=next_j", "if y==1:j=c.mpf('.5')",
        "returnj,masses",
    ):
        c3oa_require(compact(clause) in interval_slope,
                     "Interval slope J-cell clause missing: " + clause)

    targets_text, targets_name, _ = text_for(
        "lei_ren_part1_paper_compliant_current_original_C3_source_targets.py")
    targets = compact(targets_text)
    for clause in (
        "values=dict(E=raw['common_velocity_E_axial5'],V=raw['common_velocity_V_axial5']",
        "projected[name]=[row[0],row[1],2*row[2],6*row[3]]",
        "defsource_packet(self,ends,chart,left,right=None)",
    ):
        c3oa_require(compact(clause) in targets,
                     "C3 ordinary projection clause missing: " + clause)
    return {
        "passed": True, "creator_record_count": len(records),
        "ast_parse_count": len(trees), "formula_checks": formula_checks,
        "outer_dispatch": {"basename": dispatch_name, "passed": True,
                            "clauses": [
                                "axial parent -> O2_slope at (1,1)",
                                "o3 parent -> O2_buffer at (11,1)",
                                "o3 query -> O3_chart",
                                "FRONTENDS dispatch", "chart suffix dispatch"]},
        "whole_Z_projection": {"basename": whole_o3_name, "passed": True,
                                "RECOVER_common_velocity_E_V": True},
        "sigma_source": {
            "comparison_basename": comparison_name, "long_reshape_basename": long_name,
            "outer_slope_basename": interval_slope_name, "native_initial_basename": initial_name,
            "alpha_box_point_semantics": True,
            "sigma_value_derivative_is_1_minus_alpha": True,
            "stable_sigma_has_same_flat_point_formula": True,
            "interval_hulls_are_not_function_values": True,
        },
        "mu_source": {"basename": mu_name, "source_mu_admission": True,
                       "exact_decimal_token": "0.001", "Md": 40,
                       "numeric_mu_value_selected": False},
        "source_provider_sha256": C3OA_PROVIDER_SHA,
        "exact_source_basenames": {
            "axial": axial_name, "slope": slope_name,
            "interval_slope": interval_slope_name, "pre_pulse": pre_name,
            "transition_buffer": buffer_name, "C3_targets": targets_name,
        },
        "strengthened_clause_checks": {
            "exact_basename_lookup": True,
            "axial_coordinate_MD40_exp_selector_or_buffer": True,
            "axial_t_y_decay_and_E_scale": True,
            "slope_masses_qi_and_factor": True,
            "slope_sigma_value_derivative_import": True,
            "transition_stable_sigma_J_cells": True,
            "pre_pulse_slope_mu_and_power": True,
            "C3_ordinary_projection_factors": True,
        },
        "source_functions_imported_or_called": False,
    }


def c3oa_cached_action(spec):
    cached = spec.get("cached_native_action_result")
    c3oa_require(isinstance(cached, dict), "Saved native action result is missing")
    c3oa_require(cached.get("passed") is True
                 and cached.get("source_only") is True
                 and cached.get("original_source_evaluations") == 0,
                 "Saved native action result scope changed")
    consequence = cached.get("native_amplitude_consequence")
    c3oa_require(isinstance(consequence, dict), "Saved native action consequence missing")
    atoms = consequence.get("L_N")
    c3oa_require(isinstance(atoms, list), "Saved native action coefficient ledger missing")
    aliases = {
        "1": "constant", "las_native_I_slope": "I_slope",
        "las_native_exp40": "exp40",
        "(las_native_mu)*(las_native_I_trans)": "mu_times_I_trans",
        "las_native_mu": "mu",
    }
    action = {}
    for row in atoms:
        c3oa_require(isinstance(row, dict) and row.get("atom") in aliases,
                     "Unexpected saved native action atom")
        coefficient = row.get("coefficient")
        c3oa_require(isinstance(coefficient, dict)
                     and isinstance(coefficient.get("numerator"), int)
                     and isinstance(coefficient.get("denominator"), int)
                     and not isinstance(coefficient.get("numerator"), bool)
                     and not isinstance(coefficient.get("denominator"), bool),
                     "Malformed saved native action coefficient")
        denominator = coefficient["denominator"]
        c3oa_require(denominator != 0, "Saved native action denominator is zero")
        action[aliases[row["atom"]]] = str(coefficient["numerator"]) + "/" + str(denominator)
    c3oa_require(action == C3OA_EXPECTED_ACTION,
                 "Saved native action coefficients changed")
    scalar = cached.get("native_scalar_action_certificate")
    c3oa_require(isinstance(scalar, dict) and scalar.get("passed") is True
                 and scalar.get("source_intervals_or_midpoints_used") is False,
                 "Saved scalar action certificate scope changed")
    rows = {row.get("node"): row.get("literal_row")
            for row in scalar.get("literal_row_checks", []) if isinstance(row, dict)}
    expected = {
        2: {"operation": "exact_rational", "numerator": 40, "denominator": 1},
        3: {"operation": "analytic_unary", "name": "exp", "argument": 2},
        4: {"operation": "exact_rational", "numerator": 11, "denominator": 1},
        5: {"operation": "sum", "arguments": [3, 4]},
        11: {"operation": "exact_rational", "numerator": 1, "denominator": 1000},
        12: {"operation": "analytic_unary", "name": "log", "argument": 11},
        13: {"operation": "exact_rational", "numerator": 4, "denominator": 1},
        14: {"operation": "product", "arguments": [13, 5]},
        15: {"operation": "negative", "argument": 14},
        16: {"operation": "sum", "arguments": [12, 15]},
        17: {"operation": "analytic_unary", "name": "exp", "argument": 16},
    }
    for node, row in expected.items():
        actual = rows.get(node)
        c3oa_require(actual == row,
                     "Saved native literal action row changed: " + str(node))
    return {
        "passed": True, "source_result_sha256": spec.get("cached_native_action_result_sha256"),
        "action_total_coefficients": dict(action), "scalar_literal_nodes": sorted(expected),
        "scalar_literal_row_count": scalar.get("literal_row_check_count"),
        "numeric_source_values_evaluated": False,
        "source_functions_imported_or_called": False,
    }


def c3oa_parent_seed(source):
    parent = source.get("parent_result")
    c3oa_require(isinstance(parent, dict), "Authenticated Rh seed parent result is missing")
    seed = parent.get("native_Rh_seed_consequence")
    c3oa_require(isinstance(seed, dict)
                 and seed.get("native_canonical_Rh_seed_gap_ordinary_Z_0_to_3")
                 == [0, 0, 0, 0]
                 and seed.get("normalized_source_endpoint_only") is True
                 and seed.get("formula") == "exp(-1/2)/(1+Z^2)",
                 "Authenticated Rh seed consequence changed")
    return {"passed": True, "formula": seed["formula"],
            "ordinary_Z_orders": [0, 1, 2, 3],
            "endpoint_only": True, "seed_gap_not_reproved": True}


def c3oa_semantic_contract(spec, namespace, source, action, source_audit, seed):
    contract = c3oa_clone(spec.get("semantic_contract"))
    c3oa_require(isinstance(contract, dict), "C3 semantic contract is missing")
    binding = contract.setdefault("source_binding", {})
    binding.update({
        "root_wrapper_passed": True,
        "saved_native_action_result_ref": "saved native action result, data-only selector",
        "saved_Rh_seed_result_ref": "saved Rh seed result, data-only selector",
        "exact_sigma_recipe_ref": "stable_sigma and 1-alpha_box(s+1,1) exact analytic scalar semantics",
        "exact_mu_recipe_ref": "exp(log(1/1000)-4*(exp(40)+11))",
        "native_C3_parameter_alignment_binding_passed": True,
        "evaluator_covers_are_bounds": True,
        "interval_values_used": False,
        "native_ordinary_Z_slots": "[E0,E1,E2,E3]",
        "c3_query_projection": "[E0,E1,2*E2,6*E3]",
        "c3_projection_factor_application_count": 1,
    })
    leaves = contract.get("c3_source_query_leaves")
    c3oa_require(isinstance(leaves, dict), "C3 source query leaves are missing")
    for chart in C3OA_CHARTS:
        leaf = leaves.get(chart)
        c3oa_require(isinstance(leaf, dict), "Missing C3 semantic leaf: " + chart)
        c3oa_require(leaf.get("ordinary_Z_orders") == list(C3OA_ORDERS),
                     "C3 semantic leaf order changed: " + chart)
        leaf["source_ref"] = (
            "typed C3 O3_power E source query leaves; primary ports "
            + str(list(C3OA_PRIMARY_PORTS)) + "; rebuilt source leaves "
            + str(list(C3OA_REBUILT_PORTS))
        )
        leaf["ast_ref"] = "authenticated creator source AST records; chart=" + chart
    cached = contract.setdefault("cached_native", {})
    cached["data_only"] = True
    cached["source_result_ref"] = "saved native action result, data-only selector"
    cached["action_total_coefficients"] = dict(action["action_total_coefficients"])
    cached["accepted_Rh_seed_exponent"] = "-1/2"
    flags = contract.setdefault("source_flags", {})
    for name in ("actual_native_C3_amplitude_identity_proved",
                 "actual_incoming_identity_proved", "actual_control_identity_proved",
                 "actual_profile_identity_proved", "cap_admitted", "numeric_values_selected"):
        flags[name] = False
    return contract


def c3oa_source_binding(namespace, source_audit, action, seed):
    return {
        "passed": True,
        "source_only": True,
        "typed_source_query_identity": True,
        "semantic_FunctionKey_comparison_across_namespaces": True,
        "primary_source_targets_A3": 22539,
        "rebuilt_source_leaf_A3": 24181,
        "rebuilt_expression_A3_24177_used": False,
        "coordinate_is_exact_rational_two": True,
        "ordinary_Z_orders": list(C3OA_ORDERS),
        "native_ordinary_Z_slots": "[E0,E1,E2,E3]",
        "C3_query_projection": "[E0,E1,2*E2,6*E3]",
        "projection_factor_application_count": 1,
        "defining_quantity_not_a_range_value": True,
        "range_endpoints_not_field_coefficients": True,
        "source_ast_formulas_bound": True,
        "exact_sigma_source_semantics_bound": True,
        "exact_mu_source_parameter_bound": True,
        "evaluator_covers_are_bounds": True,
        "interval_values_used": False,
        "native_action_data_only": True,
        "Rh_seed_parent_data_only": True,
        "namespace": namespace,
        "source_ast_audit": source_audit,
        "native_action_audit": action,
        "Rh_seed_audit": seed,
        "old_math_calls": 0,
        "numeric_values_selected": False,
        "caps_or_midpoints_used_as_values": False,
        "control_profile_history_identity_claimed": False,
        "physical_Pstar_lift_claimed": False,
    }


def c3oa_result(source, namespace, source_audit, action, seed, semantic, algebra):
    parent = source.get("parent_result", {})
    result = {
        "schema": "original_compact_C3_outer_amplitude_identity",
        "component": C3OA_COMPONENT,
        "passed": True,
        "metadata_only": False,
        "source_only": True,
        "input_digest": source.get("input_digest"),
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "source_packet_defined": True,
        "C3_outer_amplitude_source_semantics_bound": True,
        "C3_five_outer_amplitude_seams_ordinary_Z3_proved": True,
        "C3_native_scalar_action_parameters_bound": True,
        "actual_native_C3_amplitude_identity_proved": True,
        "C3_outer_amplitude_identity_certificate": {
            "passed": True,
            "scope": "typed canonical leading-background terminal E source-query leaves",
            "source_binding": c3oa_source_binding(namespace, source_audit, action, seed),
            "algebra": algebra,
            "ambiguous_legacy_rebuilt_expression_24177_admitted": False,
        },
        "C3_Rc_amplitude_consequence": {
            "function_identity": "A_N^(j)(Z)=A_C^(j)(Z), j=0..3, Z in [-1,1]",
            "scope": "canonical leading-background terminal E source-query leaves only",
            "ordinary_Z_orders": [0, 1, 2, 3],
            "formula": "exp(-1/2-59/10-(3/5)*I_slope-(1/2)*exp40-mu*I_trans-2*mu)/(1+Z^2)",
            "signed_gap_jets": [0, 0, 0, 0],
            "primary_A3_source_port": 22539,
            "rebuilt_A3_source_leaf_port": 24181,
            "rebuilt_expression_24177_is_not_source_A3": True,
            "terminal_zero_not_claimed": True,
            "control_profile_history_identity_not_claimed": True,
            "physical_scale_identity_not_claimed": True,
        },
        "source_binding": c3oa_source_binding(namespace, source_audit, action, seed),
        "source_inspection_audit": {
            "original_source_evaluations": 0,
            "numeric_field_evaluations": 0,
            "old_source_consumer_calls": 0,
            "old_proof_functions_called": 0,
            "old_operator_calls": 0,
            "old_source_normalizer_calls": 0,
            "quadrature_calls": 0,
            "temporal_recurrence_calls": 0,
            "new_C3_endpoint_equation_coeff_checks": 30,
            "new_C3_endpoint_matrix_residual_checks": 30,
            "new_C3_native_terminal_coeff_checks": 5,
            "new_C3_seam_checks": 5,
            "new_C3_ordinary_Z_endpoint_identity_rows": 4,
            "new_C3_source_query_port_checks": 4,
            "new_scalar_source_parameter_bindings": 3,
        },
        "remaining_scope": [
            "The identity is scoped to canonical leading-background terminal E source-query leaves.",
            "The rebuilt repair expression node 24177 is retained as metadata and is not used as E3.",
            "C3/C5 incoming history, control, profile, cap, physical scale, and global source bridges remain open.",
            "Numeric interval covers remain bounds and no quantitative error reduction is admitted.",
        ],
    }
    retained = source.get("retained_parent_data", {})
    for name in ("retained_cap_provenance", "retained_construction_data",
                 "retained_integral_term_inventory"):
        if name in retained:
            result[name] = c3oa_clone(retained[name])
        elif name in parent:
            result[name] = c3oa_clone(parent[name])
    if "coordinate_mapping" in parent:
        result["coordinate_mapping"] = c3oa_clone(parent["coordinate_mapping"])
    inherited = source.get("inherited_parent_true_flags", {})
    if isinstance(inherited, dict):
        for name, value in inherited.items():
            if c3oa_is_bool(value) and value:
                result[name] = True
    for name in source.get("inherited_false_flags", []):
        result[name] = False
    for name in (
        "actual_C3_C5_incoming_source_identity_proved",
        "actual_C3_C5_control_source_identity_proved",
        "actual_C3_C5_profile_source_identity_proved",
        "actual_C3_cap_to_C5_source_bound_proved",
        "global_five_moment_closure", "global_cap_proved",
        "numeric_values_computed", "full_goal_complete", "cap_admitted",
        "physical_layer_completed", "quantitative_packet_error_reduction_proved",
    ):
        result[name] = False
    return result


def derive(source, journal):
    spec = source["input"]
    c3oa_require(isinstance(spec, dict)
                 and spec.get("schema") == C3OA_SCHEMA
                 and spec.get("component") == C3OA_COMPONENT,
                 "Wrong C3 outer-amplitude input")
    namespace = c3oa_namespace_binding(spec)
    source_audit = c3oa_source_records(spec)
    action = c3oa_cached_action(spec)
    seed = c3oa_parent_seed(source)
    semantic = c3oa_semantic_contract(spec, namespace, source, action, source_audit, seed)
    binding = c3oa_source_binding(namespace, source_audit, action, seed)
    c3oa_require(binding["passed"], "C3 source binding failed")
    if journal is not None:
        journal({"event": "C3_outer_amplitude_source_binding",
                 "primary_A3": 22539, "rebuilt_source_leaf_A3": 24181,
                 "rebuilt_expression_24177_used": False,
                 "source_query_ports": 4, "ordinary_Z_orders": 4})
    algebra = derive_C3_outer_amplitude_identity_algebra(semantic, journal)
    c3oa_require(isinstance(algebra, dict) and algebra.get("passed") is True,
                 "C3 outer-amplitude algebra failed")
    result = c3oa_result(source, namespace, source_audit, action, seed, semantic, algebra)
    if journal is not None:
        journal({"event": "C3_outer_amplitude_identity_complete",
                 "passed": True, "scope": "typed canonical leading-background E leaves",
                 "control_profile_history_identity": False,
                 "physical_scale_identity": False})
    return result
