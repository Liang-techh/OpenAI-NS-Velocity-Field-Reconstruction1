"""Unexecuted exact algebra draft for the quiet-power kernel specialization.

This is a new bounded family.  It does not import a source provider, replay the
generic paired-Duhamel family, evaluate an exponential, integrate a field, or
inspect an ancestor artifact.  The only arithmetic is exact rational arithmetic
on integer numerator/denominator pairs.  The module intentionally has no
top-level calls; the entry point is reserved for a later frozen execution.
"""


QUIET_POWER_KERNEL_INPUT_SCHEMA = "compact_quiet_power_kernel_algebra_input_v4"
QUIET_POWER_KERNEL_RESULT_SCHEMA = "compact_quiet_power_kernel_algebra_result_v4"
QUIET_POWER_KERNEL_COMPONENT = "compact_quiet_power_kernel_algebra"
QUIET_POWER_FIELDS = ("m", "h", "k", "e", "p")
QUIET_POWER_ORDERS = (0, 1, 2, 3)
QUIET_POWER_RATE_EXPECTED = {
    "m": (1, 1),
    "h": (3, 2),
    "k": (3, 2),
    "e": (1, 1),
    "p": (0, 1),
}


QUIET_POWER_KERNEL_CONTRACT_TEMPLATE = {
    "schema": QUIET_POWER_KERNEL_INPUT_SCHEMA,
    "component": QUIET_POWER_KERNEL_COMPONENT,
    "source_binding": {
        "passed": True,
        "actual_quiet_power_kernel_source_bound": True,
    },
    "native_source_geometry": {
        "actual_O3_power_affine": True,
        "O3_chart_stage": 3,
        "order_scope": "O3 chart stage; not temporal n=3",
        "log_offset": {
            "expression": "exp(40)+18+x",
            "exp_argument": [40, 1],
            "constant": [18, 1],
            "slope": [1, 1],
            "variable": "x",
        },
        "jacobian": [1, 1],
        "lower": [0, 1],
        "upper": [2, 1],
        "Z_independent": True,
        "function_keys": {
            "m": "exp(-2)",
            "h": "exp(-3)",
            "k": "exp(-3)",
            "e": "exp(-2)",
            "p": "1",
        },
    },
    "C3_source_geometry": {
        "actual_O3_power_affine": True,
        "O3_chart_stage": 3,
        "order_scope": "O3 chart stage; not temporal n=3",
        "log_offset": {
            "expression": "exp(40)+18+x",
            "exp_argument": [40, 1],
            "constant": [18, 1],
            "slope": [1, 1],
            "variable": "x",
        },
        "jacobian": [1, 1],
        "lower": [0, 1],
        "upper": [2, 1],
        "Z_independent": True,
        "function_keys": {
            "m": "exp(-2)",
            "h": "exp(-3)",
            "k": "exp(-3)",
            "e": "exp(-2)",
            "p": "1",
        },
    },
    "rates": {
        "native": {
            "m": [1, 1],
            "h": [3, 2],
            "k": [3, 2],
            "e": [1, 1],
            "p": [0, 1],
        },
        "C3": {
            "m": [1, 1],
            "h": [3, 2],
            "k": [3, 2],
            "e": [1, 1],
            "p": [0, 1],
        },
    },
    "native_flat_scalar_source_gate": {
        "source_assignments_bound": True,
        "Ey_over_E_coefficients_bound": True,
        "Ey_over_E_coefficients": [[-1, 2], [-1, 1]],
        "Vy_source_zero_bound": True,
        "Vy_over_E_coefficients": [[0, 1], [0, 1]],
        "E_nonzero_source_bound": True,
        "E_positive_source_bound": True,
        "mu_positive_source_bound": True,
        "mu_Z_independent": True,
        "cached_2mu_ge_eta_bound": True,
        "actual_lazy_predicate_bound": True,
        "root_gate": True,
    },
    "source_quiet_local_density_zero": {
        "passed": True,
        "data_gate": True,
        "native_rows": 20,
        "C3_rows": 20,
        "recomputed": False,
        "C3_fixed_zero_root_verified": True,
    },
    "cached_generic_paired_transport": {
        "passed": True,
        "data_gate": True,
        "check_count": 20,
        "source": "parent_saved_result",
        "recomputed": False,
    },
    "source_entry": {
        "delta_not_zeroed": True,
        "unknown_entry_gap_retained": True,
        "values_uncomputed": True,
    },
    "scope": {
        "difference_family": "five finite-N correction history channels",
        "fields": ["m", "h", "k", "e", "p"],
        "ordinary_Z_orders": [0, 1, 2, 3],
        "not_m15": True,
        "not_P1": True,
        "not_full_NS": True,
        "physical_Pstar_lift_claimed": False,
        "global_completion_claimed": False,
        "actual_whole_17_kernel_pair_claimed": False,
    },
}


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _is_integer(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _gcd(left, right):
    left = -left if left < 0 else left
    right = -right if right < 0 else right
    while right:
        left, right = right, left % right
    return left or 1


def _rat(value, label):
    if _is_integer(value):
        numerator, denominator = value, 1
    elif isinstance(value, (list, tuple)) and len(value) == 2:
        numerator, denominator = value
    elif isinstance(value, dict):
        numerator = value.get("numerator", value.get("num"))
        denominator = value.get("denominator", value.get("den"))
    else:
        numerator, denominator = None, None
    _require(_is_integer(numerator) and _is_integer(denominator),
             label + " is not an exact integer pair")
    _require(denominator != 0, label + " has a zero denominator")
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    divisor = _gcd(numerator, denominator)
    return numerator // divisor, denominator // divisor


def _rat_add(left, right):
    return _rat((left[0] * right[1] + right[0] * left[1],
                 left[1] * right[1]), "rational sum")


def _rat_sub(left, right):
    return _rat((left[0] * right[1] - right[0] * left[1],
                 left[1] * right[1]), "rational difference")


def _rat_mul(left, right):
    return _rat((left[0] * right[0], left[1] * right[1]),
                "rational product")


def _rat_neg(value):
    return -value[0], value[1]


def _rat_json(value):
    return [value[0], value[1]]


def _rat_text(value):
    if value[1] == 1:
        return str(value[0])
    return str(value[0]) + "/" + str(value[1])


def _exp_key(exponent):
    if exponent == (0, 1):
        return "1"
    return "exp(" + _rat_text(exponent) + ")"


def _scalar_gcd(left, right):
    left = -left if left < 0 else left
    right = -right if right < 0 else right
    while right:
        left, right = right, left % right
    return left or 1


def _scalar_rat(value, label):
    if _is_integer(value):
        numerator, denominator = value, 1
    elif isinstance(value, (list, tuple)) and len(value) == 2:
        numerator, denominator = value
    elif isinstance(value, dict):
        numerator = value.get("numerator", value.get("num"))
        denominator = value.get("denominator", value.get("den"))
    else:
        numerator, denominator = None, None
    _require(_is_integer(numerator) and _is_integer(denominator),
             label + " is not an exact scalar coefficient")
    _require(denominator != 0, label + " has a zero denominator")
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    divisor = _scalar_gcd(numerator, denominator)
    return numerator // divisor, denominator // divisor


def _scalar_vec_add(left, right):
    return (
        _scalar_rat((left[0][0] * right[0][1] +
                     right[0][0] * left[0][1],
                     left[0][1] * right[0][1]), "scalar constant sum"),
        _scalar_rat((left[1][0] * right[1][1] +
                     right[1][0] * left[1][1],
                     left[1][1] * right[1][1]), "scalar mu sum"),
    )


def _scalar_vec_sub(left, right):
    return (
        _scalar_rat((left[0][0] * right[0][1] -
                     right[0][0] * left[0][1],
                     left[0][1] * right[0][1]), "scalar constant difference"),
        _scalar_rat((left[1][0] * right[1][1] -
                     right[1][0] * left[1][1],
                     left[1][1] * right[1][1]), "scalar mu difference"),
    )


def _scalar_vec_from_coefficients(value, label):
    _require(isinstance(value, (list, tuple)) and len(value) == 2,
             label + " must have constant and mu coefficients")
    return (
        _scalar_rat(value[0], label + " constant"),
        _scalar_rat(value[1], label + " mu"),
    )


def _scalar_vec_scale(value, factor):
    factor = _scalar_rat(factor, "scalar vector factor")
    return (
        _scalar_rat((value[0][0] * factor[0], value[0][1] * factor[1]),
                    "scaled scalar constant"),
        _scalar_rat((value[1][0] * factor[0], value[1][1] * factor[1]),
                    "scaled scalar mu"),
    )


def _scalar_vec_json(value):
    return [[value[0][0], value[0][1]], [value[1][0], value[1][1]]]


def _read_native_flat_scalar_gate(contract):
    gate = contract.get("native_flat_scalar_source_gate")
    _require(isinstance(gate, dict), "native flat scalar source gate is absent")
    for key in (
        "source_assignments_bound",
        "Ey_over_E_coefficients_bound",
        "Vy_source_zero_bound",
        "E_nonzero_source_bound",
        "E_positive_source_bound",
        "mu_positive_source_bound",
        "mu_Z_independent",
        "cached_2mu_ge_eta_bound",
        "actual_lazy_predicate_bound",
        "root_gate",
    ):
        _require(gate.get(key) is True,
                 "native flat scalar source gate is missing: " + key)
    return {
        "source_assignments_bound": True,
        "Ey_over_E_coefficients_bound": True,
        "Ey_over_E_coefficients": [[-1, 2], [-1, 1]],
        "Vy_source_zero_bound": True,
        "Vy_over_E_coefficients": [[0, 1], [0, 1]],
        "E_nonzero_source_bound": True,
        "E_positive_source_bound": True,
        "mu_positive_source_bound": True,
        "mu_Z_independent": True,
        "cached_2mu_ge_eta_bound": True,
        "actual_lazy_predicate_bound": True,
        "root_gate": True,
        "recomputed": False,
    }


def _build_native_scalar_affine_checks(native_gate):
    _require(native_gate.get("root_gate") is True,
             "native scalar affine branch lacks root gate")
    _require(native_gate.get("E_positive_source_bound") is True,
             "native scalar affine branch lacks E positivity gate")
    _require(native_gate.get("mu_positive_source_bound") is True,
             "native scalar affine branch lacks mu positivity gate")
    ey = _scalar_vec_from_coefficients(
        native_gate["Ey_over_E_coefficients"], "literal Ey/E")
    vy = _scalar_vec_from_coefficients(
        native_gate["Vy_over_E_coefficients"], "literal Vy/E")
    _require(ey == ((-1, 2), (-1, 1)),
             "literal Ey/E source coefficients changed")
    _require(vy == ((0, 1), (0, 1)),
             "literal Vy/E source coefficients changed")
    one = ((1, 1), (0, 1))
    two = ((2, 1), (0, 1))
    two_ey = _scalar_vec_scale(ey, (2, 1))
    a = _scalar_vec_sub(one, two_ey)
    b = _scalar_vec_scale(vy, (2, 1))
    zero = ((0, 1), (0, 1))
    _require(native_gate.get("Vy_source_zero_bound") is True,
             "native scalar branch lacks Vy zero gate")
    _require(b == zero, "native scalar b is not zero")
    _require(a == ((2, 1), (2, 1)), "native scalar a changed")
    kappa = a
    _require(a[0][0] > 0 and a[1][0] >= 0,
             "native scalar a positivity is not source-bound")
    delta = _scalar_vec_sub(kappa, two)
    _require(delta == ((0, 1), (2, 1)), "native scalar Delta changed")
    return {
        "root_source_assignments": {
            "Ey_over_E": "-(1/2+mu)",
            "Vy": "0",
            "E": "positive",
            "a": "(E-2Ey)/E",
            "b": "2Vy/E",
            "kappa": "a+b^2/a",
            "Delta": "kappa-2",
        },
        "coefficient_basis": ["constant", "mu"],
        "literal_Ey_over_E": "-(1/2+mu)",
        "literal_Ey_over_E_coefficients": _scalar_vec_json(ey),
        "literal_Vy_over_E_coefficients": _scalar_vec_json(vy),
        "primary_a_check": {
            "identity": "a=(E-2Ey)/E=1-2(Ey/E)",
            "one_coefficients": _scalar_vec_json(one),
            "two_Ey_over_E_coefficients": _scalar_vec_json(two_ey),
            "a_coefficients": _scalar_vec_json(a),
            "passed": True,
        },
        "a_expression": "2+2mu",
        "b_expression": "0",
        "Delta_expression": "2mu",
        "a": _scalar_vec_json(a),
        "b": _scalar_vec_json(b),
        "kappa": _scalar_vec_json(kappa),
        "Delta": _scalar_vec_json(delta),
        "b_derived_from_Vy_zero_gate": True,
        "kappa_equals_a_b_zero": True,
        "a_positive_from_source_gates": True,
        "division_by_zero_introduced": False,
        "source_E_evaluated": False,
        "source_mu_evaluated": False,
        "generic_DAG_walk_performed": False,
        "old_proof_replayed": False,
    }


def _build_native_flat_density_zero_consequence(source_gate, native_gate):
    _require(source_gate.get("C3_fixed_zero_root_verified") is True,
             "C3 fixed-zero root gate is absent")
    _require(native_gate.get("root_gate") is True,
             "native flat scalar root gate is absent")
    return {
        "conditional_on_root_source_gates": True,
        "native_source_rows": 20,
        "C3_source_rows": 20,
        "C3_fixed_zero_root_verified": True,
        "native_source_body": {
            "A": "0",
            "B": "0",
            "r": "E*expm1(0)=0",
            "s": "0/N=0",
            "functions_identically_zero": True,
            "ordinary_Z_rows_zero": 20,
        },
        "N_positive_inherited": True,
        "no_divided_by_zero": True,
        "N_computed": False,
        "exponential_evaluated": False,
        "source_values_evaluated": False,
        "actual_terminal_gap_zero_claimed": False,
        "root_gate_proved_here": False,
    }


def _read_geometry(value, label):
    _require(isinstance(value, dict), label + " geometry is not a mapping")
    _require(value.get("actual_O3_power_affine") is True,
             label + " actual O3 power affine binding is absent")
    chart_stage = value.get("O3_chart_stage")
    _require(chart_stage == 3,
             label + " O3 chart stage is not 3")
    _require(value.get("order_scope", "O3 chart stage; not temporal n=3") ==
             "O3 chart stage; not temporal n=3",
             label + " order scope is not chart-scoped")
    offset = value.get("log_offset")
    _require(isinstance(offset, dict), label + " affine log offset is absent")
    _require(offset.get("expression") == "exp(40)+18+x",
             label + " affine log offset expression changed")
    _require(offset.get("variable") == "x",
             label + " affine log offset variable changed")
    exp_argument = _rat(offset.get("exp_argument"), label + " exp argument")
    constant = _rat(offset.get("constant"), label + " additive constant")
    slope = _rat(offset.get("slope"), label + " affine slope")
    _require(exp_argument == (40, 1), label + " exp argument changed")
    _require(constant == (18, 1), label + " additive constant changed")
    _require(slope == (1, 1), label + " affine slope changed")
    jacobian = _rat(value.get("jacobian"), label + " jacobian")
    lower = _rat(value.get("lower"), label + " lower endpoint")
    upper = _rat(value.get("upper"), label + " upper endpoint")
    _require(jacobian == (1, 1), label + " jacobian changed")
    _require(lower == (0, 1) and upper == (2, 1),
             label + " quiet interval changed")
    _require(value.get("Z_independent") is True,
             label + " is not marked Z independent")
    return {
        "label": label,
        "O3_chart_stage": 3,
        "order_scope": "O3 chart stage; not temporal n=3",
        "exp_argument": exp_argument,
        "constant": constant,
        "slope": slope,
        "jacobian": jacobian,
        "lower": lower,
        "upper": upper,
        "Z_independent": True,
        "function_keys": value.get("function_keys", value.get("exp_function_keys")),
    }


def _read_rates(contract):
    rates = contract.get("rates")
    _require(isinstance(rates, dict), "quiet power rates are absent")
    native = rates.get("native", rates.get("lambda_native"))
    c3 = rates.get("C3", rates.get("lambda_C3"))
    _require(isinstance(native, dict) and isinstance(c3, dict),
             "native and C3 quiet power rates are absent")
    result = {}
    for field in QUIET_POWER_FIELDS:
        native_rate = _rat(native.get(field), "native lambda " + field)
        c3_rate = _rat(c3.get(field), "C3 lambda " + field)
        _require(native_rate == QUIET_POWER_RATE_EXPECTED[field],
                 "native lambda changed: " + field)
        _require(c3_rate == QUIET_POWER_RATE_EXPECTED[field],
                 "C3 lambda changed: " + field)
        _require(native_rate == c3_rate,
                 "native/C3 lambda mismatch: " + field)
        result[field] = {"native": native_rate, "C3": c3_rate, "same": True}
    return result


def _read_source_gate(contract):
    gate = contract.get("source_quiet_local_density_zero")
    _require(isinstance(gate, dict), "quiet local density source gate is absent")
    native_rows = gate.get("native_rows", gate.get("native_zero_rows"))
    c3_rows = gate.get("C3_rows", gate.get("C3_zero_rows"))
    _require(gate.get("passed") is True and gate.get("data_gate") is True,
             "quiet local density source gate is not passed data")
    _require(native_rows == 20 and c3_rows == 20,
             "quiet local density source gate row count changed")
    _require(gate.get("C3_fixed_zero_root_verified") is True,
             "C3 fixed-zero root gate is absent")
    return {
        "passed": True,
        "data_gate": True,
        "native_rows": 20,
        "C3_rows": 20,
        "recomputed": False,
        "C3_fixed_zero_root_verified": True,
    }


def _read_generic_gate(contract):
    gate = contract.get("cached_generic_paired_transport")
    _require(isinstance(gate, dict), "cached generic paired transport is absent")
    count = gate.get("check_count", gate.get("rows"))
    _require(gate.get("passed") is True and gate.get("data_gate") is True,
             "cached generic paired transport is not a passed data gate")
    _require(count == 20, "cached generic paired transport row count changed")
    return {
        "passed": True,
        "data_gate": True,
        "check_count": 20,
        "source": gate.get("source", "parent_saved_result"),
        "recomputed": False,
    }


def _read_function_keys(native, c3, expected):
    native_keys = native.get("function_keys")
    c3_keys = c3.get("function_keys")
    _require(isinstance(native_keys, dict) and isinstance(c3_keys, dict),
             "quiet power function keys are absent")
    result = {}
    for field in QUIET_POWER_FIELDS:
        _require(native_keys.get(field) == expected[field],
                 "native exact exponential function key changed: " + field)
        _require(c3_keys.get(field) == expected[field],
                 "C3 exact exponential function key changed: " + field)
        _require(native_keys.get(field) == c3_keys.get(field),
                 "native/C3 exact exponential function key mismatch: " + field)
        result[field] = {
            "native": native_keys[field],
            "C3": c3_keys[field],
            "identical": True,
        }
    return result


def _journal(journal, event, **payload):
    if journal is not None:
        _require(callable(journal), "journal must be callable")
        journal({"event": event, **payload})


def derive_quiet_power_kernel_algebra(contract, journal=None):
    """Derive the bounded quiet-power transfer factors from a literal contract."""
    _require(isinstance(contract, dict), "quiet power contract is not a mapping")
    _require(contract.get("schema") == QUIET_POWER_KERNEL_INPUT_SCHEMA,
             "quiet power input schema changed")
    _require(contract.get("component") == QUIET_POWER_KERNEL_COMPONENT,
             "quiet power component changed")
    binding = contract.get("source_binding")
    binding_passed = binding.get("passed") if isinstance(binding, dict) else (
        contract.get("source_binding_passed"))
    _require(binding_passed is True, "quiet power source binding is not passed")

    native = _read_geometry(contract.get("native_source_geometry"), "native")
    c3 = _read_geometry(contract.get("C3_source_geometry"), "C3")
    rates = _read_rates(contract)
    native_flat_gate = _read_native_flat_scalar_gate(contract)
    native_scalar = _build_native_scalar_affine_checks(native_flat_gate)
    source_gate = _read_source_gate(contract)
    native_zero = _build_native_flat_density_zero_consequence(
        source_gate, native_flat_gate)
    generic_gate = _read_generic_gate(contract)
    entry = contract.get("source_entry")
    _require(isinstance(entry, dict) and
             entry.get("delta_not_zeroed") is True and
             entry.get("unknown_entry_gap_retained") is True,
             "source entry gap must remain unknown and unzeroed")
    _require(entry.get("values_uncomputed") is True,
             "source entry values must remain uncomputed")

    native_width = _rat_sub(_rat_mul(native["slope"], native["upper"]),
                            _rat_mul(native["slope"], native["lower"]))
    c3_width = _rat_sub(_rat_mul(c3["slope"], c3["upper"]),
                        _rat_mul(c3["slope"], c3["lower"]))
    _require(native_width == (2, 1) and c3_width == (2, 1),
             "quiet power width is not exactly 2")

    native_suffix = {
        "constant": _rat_mul(native["slope"], native["upper"]),
        "slope": _rat_neg(native["slope"]),
        "variable": "s",
        "expression": "2-s",
    }
    c3_suffix = {
        "constant": _rat_mul(c3["slope"], c3["upper"]),
        "slope": _rat_neg(c3["slope"]),
        "variable": "s",
        "expression": "2-s",
    }
    _require(native_suffix["constant"] == (2, 1) and
             c3_suffix["constant"] == (2, 1) and
             native_suffix["slope"] == (-1, 1) and
             c3_suffix["slope"] == (-1, 1),
             "quiet power suffix offset changed")

    factors = {}
    expected_keys = {}
    for field in QUIET_POWER_FIELDS:
        exponent = _rat_neg(_rat_mul(native_width, rates[field]["native"]))
        function_key = _exp_key(exponent)
        factors[field] = {
            "lambda_native": _rat_json(rates[field]["native"]),
            "lambda_C3": _rat_json(rates[field]["C3"]),
            "lambda_native_equals_C3": True,
            "exponent": _rat_json(exponent),
            "exact_exp_expression": function_key,
            "function_key": function_key,
            "positive_exact_factor": True,
            "is_one": exponent == (0, 1),
        }
        expected_keys[field] = function_key
    function_keys = _read_function_keys(native, c3, expected_keys)

    common_rows = []
    transfer_rows = []
    for field in QUIET_POWER_FIELDS:
        for order in QUIET_POWER_ORDERS:
            common_rows.append({
                "channel": field,
                "ordinary_Z_order": order,
                "Z_independent": True,
                "native_width": _rat_json(native_width),
                "C3_width": _rat_json(c3_width),
                "native_suffix_offset": {
                    "constant": _rat_json(native_suffix["constant"]),
                    "slope": _rat_json(native_suffix["slope"]),
                    "variable": "s",
                    "expression": "2-s",
                },
                "C3_suffix_offset": {
                    "constant": _rat_json(c3_suffix["constant"]),
                    "slope": _rat_json(c3_suffix["slope"]),
                    "variable": "s",
                    "expression": "2-s",
                },
                "native_function_key": function_keys[field]["native"],
                "C3_function_key": function_keys[field]["C3"],
                "function_key_identical": True,
                "lambda_native_equals_C3": True,
            })
            transfer_rows.append({
                "channel": field,
                "ordinary_Z_order": order,
                "Z_independent": True,
                "specialization": "deltaM=0;deltaI=0",
                "formula": "deltaT_i,Rc^j=exp(-2*lambda_i)*deltaT_i,Rw^j",
                "factor": factors[field]["exact_exp_expression"],
                "function_key": factors[field]["function_key"],
                "positive_exact_factor": True,
                "delta_out_zero_iff_delta_in_zero": True,
                "source_entry_gap_unknown": True,
                "source_entry_delta_not_zeroed": True,
                "delta_out_zero_claimed": False,
            })
    _require(len(common_rows) == 20 and len(transfer_rows) == 20,
             "quiet power ordinary Z row count changed")
    _journal(journal, "quiet_power_kernel_rows_bound",
             common_kernel_checks=20, source_conditional_transfer_rows=20,
             values_uncomputed=True)

    return {
        "schema": QUIET_POWER_KERNEL_RESULT_SCHEMA,
        "component": QUIET_POWER_KERNEL_COMPONENT,
        "passed": True,
        "source_binding_passed": True,
        "source_quiet_local_density_zero_gate": source_gate,
        "native_flat_scalar_source_gate": native_flat_gate,
        "native_scalar_affine_checks": native_scalar,
        "native_flat_density_zero_consequence": native_zero,
        "cached_generic_paired_transport_gate": generic_gate,
        "source_entry_gap_unknown": True,
        "source_entry_delta_not_zeroed": True,
        "source_values_uncomputed": True,
        "source_independent_specialization": True,
        "geometry": {
            "native": {
                "O3_chart_stage": native["O3_chart_stage"],
                "order_scope": "O3 chart stage; not temporal n=3",
                "exp_argument": _rat_json(native["exp_argument"]),
                "constant": _rat_json(native["constant"]),
                "slope": _rat_json(native["slope"]),
                "jacobian": _rat_json(native["jacobian"]),
                "lower": _rat_json(native["lower"]),
                "upper": _rat_json(native["upper"]),
                "Z_independent": True,
            },
            "C3": {
                "O3_chart_stage": c3["O3_chart_stage"],
                "order_scope": "O3 chart stage; not temporal n=3",
                "exp_argument": _rat_json(c3["exp_argument"]),
                "constant": _rat_json(c3["constant"]),
                "slope": _rat_json(c3["slope"]),
                "jacobian": _rat_json(c3["jacobian"]),
                "lower": _rat_json(c3["lower"]),
                "upper": _rat_json(c3["upper"]),
                "Z_independent": True,
            },
            "absolute_Rm_Rc_source_identity_needed": False,
            "additive_origin_cancels_in_width": True,
        },
        "width_checks": {
            "native_offset_2_minus_offset_0": [2, 1],
            "C3_offset_2_minus_offset_0": [2, 1],
            "native_suffix_offset_2_minus_offset_s": "2-s",
            "C3_suffix_offset_2_minus_offset_s": "2-s",
            "native_and_C3_widths_match": True,
            "native_and_C3_suffixes_match": True,
        },
        "channel_factors": factors,
        "function_key_checks": function_keys,
        "common_kernel_checks": common_rows,
        "source_conditional_transfer_rows": transfer_rows,
        "stage_counts": {
            "common_kernel_checks": 20,
            "source_conditional_transfer_rows": 20,
            "ordinary_Z_rows": 20,
        },
        "specialization": {
            "cached_generic_identity_selector": "$.signed_density_algebra_result.paired_history_transport",
            "deltaM_set_to_zero": True,
            "deltaI_set_to_zero": True,
            "generic_polynomials_rerun": False,
            "formula": "deltaT_i,Rc^j=exp(-2*lambda_i)*deltaT_i,Rw^j",
            "m_factor": "exp(-2)",
            "h_factor": "exp(-3)",
            "k_factor": "exp(-3)",
            "e_factor": "exp(-2)",
            "p_factor": "1",
            "positive_exact_factor_equivalence_scoped_quiet": True,
            "delta_out_zero_iff_delta_in_zero": True,
            "source_entry_gap_unknown": True,
            "source_entry_delta_not_zeroed": True,
            "actual_terminal_delta_out_zero_claimed": False,
        },
        "scope": {
            "difference_family": "five finite-N correction history channels",
            "fields": ["m", "h", "k", "e", "p"],
            "ordinary_Z_orders": [0, 1, 2, 3],
            "not_m15": True,
            "not_P1": True,
            "not_full_NS": True,
            "physical_Pstar_lift_claimed": False,
            "global_completion_claimed": False,
            "actual_whole_17_kernel_pair_claimed": False,
            "old_caps_preserved": True,
            "old_m15_terms_deleted": False,
            "old_212_terms_deleted": False,
            "N_computed": False,
            "integrals_computed": False,
            "exponentials_evaluated": False,
            "native_flat_scalar_branch_conditional": True,
            "source_quiet_local_density_zero_is_root_gate": True,
            "C3_fixed_zero_rows_root_verified": True,
        },
        "runtime_model_metadata": None,
    }


__all__ = (
    "QUIET_POWER_KERNEL_CONTRACT_TEMPLATE",
    "derive_quiet_power_kernel_algebra",
)







"""Unexecuted static C3 quiet-source binder draft.

This module intentionally performs only stdlib text, hash, AST, and bounded
literal-row checks. It never imports or executes project source, providers,
proofs, checkers, normalizers, quadratures, recurrences, or math helpers.
"""

import ast
import hashlib
import json


EXPECTED_J = "1358356628656378313"
EXPECTED_REBUILT_ARCHIVE_SHA256 = (
    "cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193"
)
EXPECTED_ARCHIVE_ROW4898_SELECTOR = "exact_graph_nodes[4898]"
EXPECTED_ARCHIVE_QUIET_SELECTOR = (
    "exact_graph_nodes[4898].recipe.quiet_source_proof"
)
EXPECTED_RECEIPT_PROOF_PROVENANCE_SELECTOR = (
    "v2_receipt.rebuilt_archive_selective_literal_slice.archive_quiet_source_proof"
)
EXPECTED_ARCHIVE_SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": (
        "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894"
    ),
    "datum_enclosure_sha256": (
        "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010"
    ),
    "implicit_source_sha256": (
        "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae"
    ),
}
EXPECTED_ARCHIVE_ROW4898_KEYS = (
    "Z_order",
    "Z_variable",
    "coordinate",
    "defining_quantity_not_a_range_value",
    "native_chart",
    "operation",
    "phase_independent_leading_source",
    "quantity",
    "recipe",
    "source_family",
)
EXPECTED_ARCHIVE_RECIPE_KEYS = (
    "AST_sha256",
    "defining_provider",
    "function",
    "leading_source_fixed_input_N0",
    "module",
    "module_sha256",
    "native_chart",
    "original_independent_P0_and_parent_leading_recipes_retained",
    "quantity_paths",
    "quiet_source_proof",
    "range_endpoints_not_field_coefficients",
)
EXPECTED_ARCHIVE_QUANTITY_PATH_KEYS = (
    "Delta",
    "E",
    "V",
    "a",
    "b",
    "p2",
    "t0",
)
EXPECTED_ARCHIVE_PROOF_KEYS = (
    "admission",
    "checked_receipt",
    "checked_receipt_sha256",
    "original_primitive_A_B_zero_but_incoming_histories_retained",
    "source_frontend_binding",
    "source_identity",
)
EXPECTED_ARCHIVE_ADMISSION_KEYS = (
    "AST_sha256",
    "function",
    "module",
    "module_sha256",
)
EXPECTED_ARCHIVE_FRONTEND_KEYS = (
    "exact_unchanged_leading_frontend_AST_sha256",
    "original_parent_leading_source_recipe_not_changed",
    "original_query_module",
    "original_query_sha256",
    "stop_before_all_u_support_and_fixed_N_density",
)
FROZEN_ENTRY_GUARD_KEYS = [
    "C3_receipt",
    "C3_receipt.receipt_version",
    "C3_receipt.corrected_N_and_invN_literals",
    "C3_receipt.corrected_N_and_invN_literals.J_exact_string",
    "C3_receipt.corrected_N_and_invN_literals.N_node_id",
    "C3_receipt.corrected_N_and_invN_literals.invN_node_id",
    "C3_receipt.corrected_N_and_invN_literals.rebuilt_node_18",
    "C3_receipt.corrected_N_and_invN_literals.rebuilt_node_18.J",
    "C3_receipt.corrected_N_and_invN_literals.primary_node_18",
    "C3_receipt.corrected_N_and_invN_literals.primary_node_18.J",
    "C3_receipt.rebuilt_archive_authorized_selection",
    "C3_receipt.rebuilt_archive_authorized_selection.raw_sha256",
    "C3_receipt.rebuilt_archive_authorized_selection.exact_graph_nodes_count",
    "C3_receipt.rebuilt_incoming_history_slots",
    "C3_receipt.rebuilt_incoming_history_slots.incoming_ids_by_field",
    "C3_receipt.rebuilt_incoming_history_slots.incoming_ids_by_field.m",
    "C3_receipt.rebuilt_incoming_history_slots.incoming_ids_by_field.h",
    "C3_receipt.rebuilt_incoming_history_slots.incoming_ids_by_field.k",
    "C3_receipt.rebuilt_incoming_history_slots.incoming_ids_by_field.e",
    "C3_receipt.rebuilt_incoming_history_slots.incoming_ids_by_field.p",
    "C3_receipt.rebuilt_incoming_history_slots.incoming_node_bodies",
    "C3_receipt.rebuilt_incoming_history_slots.incoming_node_bodies.m",
    "C3_receipt.rebuilt_incoming_history_slots.incoming_node_bodies.h",
    "C3_receipt.rebuilt_incoming_history_slots.incoming_node_bodies.k",
    "C3_receipt.rebuilt_incoming_history_slots.incoming_node_bodies.e",
    "C3_receipt.rebuilt_incoming_history_slots.incoming_node_bodies.p",
    "C3_receipt.additional_literal_rows",
    "C3_receipt.additional_literal_rows.source_guard_and_offset_rows",
    "C3_receipt.additional_literal_rows.source_guard_and_offset_rows.4584",
    "C3_receipt.additional_literal_rows.source_guard_and_offset_rows.4585",
    "C3_receipt.additional_literal_rows.source_guard_and_offset_rows.4890",
    "C3_receipt.additional_literal_rows.source_guard_and_offset_rows.4891",
    "C3_receipt.additional_literal_rows.source_guard_and_offset_rows.4892",
    "C3_receipt.additional_literal_rows.source_guard_and_offset_rows.4893",
    "C3_receipt.additional_literal_rows.source_guard_and_offset_rows.4894",
    "C3_receipt.additional_literal_rows.source_guard_and_offset_rows.4914",
    "C3_receipt.additional_literal_rows.source_guard_and_offset_rows.4915",
    "C3_receipt.additional_literal_rows.source_guard_and_offset_rows.4916",
    "C3_receipt.additional_literal_rows.source_guard_and_offset_rows.4917",
]
EXPECTED_FIELDS = ("m", "h", "k", "e", "p")
EXPECTED_ORDERS = (0, 1, 2, 3)

EXPECTED_INCOMING = {
    "m": [4820, 4821, 12010, 24145],
    "h": [4842, 4843, 12018, 24153],
    "k": [4858, 4859, 12026, 24161],
    "e": [4874, 4875, 12034, 24169],
    "p": [4888, 4889, 12041, 24176],
}

EXPECTED_TERMINALS = {
    "m": [4920, 4921, 12042, 24177],
    "h": [4928, 4929, 12043, 24178],
    "k": [4930, 4931, 12044, 24179],
    "e": [4932, 4933, 12045, 24180],
    "p": [4888, 4889, 12041, 24176],
}

EXPECTED_TERMINAL_ROWS = {
    4920: {"operation": "product", "arguments": [4919, 4820]},
    4921: {"operation": "product", "arguments": [4919, 4821]},
    12042: {"operation": "product", "arguments": [4919, 12010]},
    24177: {"operation": "product", "arguments": [4919, 24145]},
    4928: {"operation": "product", "arguments": [4927, 4842]},
    4929: {"operation": "product", "arguments": [4927, 4843]},
    12043: {"operation": "product", "arguments": [4927, 12018]},
    24178: {"operation": "product", "arguments": [4927, 24153]},
    4930: {"operation": "product", "arguments": [4927, 4858]},
    4931: {"operation": "product", "arguments": [4927, 4859]},
    12044: {"operation": "product", "arguments": [4927, 12026]},
    24179: {"operation": "product", "arguments": [4927, 24161]},
    4932: {"operation": "product", "arguments": [4919, 4874]},
    4933: {"operation": "product", "arguments": [4919, 4875]},
    12045: {"operation": "product", "arguments": [4919, 12034]},
    24180: {"operation": "product", "arguments": [4919, 24169]},
    4888: {"operation": "sum", "arguments": [4580, 4881]},
    4889: {"operation": "sum", "arguments": [4581, 4887]},
    12041: {"operation": "sum", "arguments": [11928, 12040]},
    24176: {"operation": "sum", "arguments": [24053, 24175]},
}

EXPECTED_ADDITIONAL_ROWS = {
    "4584": {"operation": "exact_rational", "numerator": 1012, "denominator": 1},
    "4585": {
        "operation": "sum",
        "arguments": [4584, 2084, 2091, 649, 3, 645, 1507],
    },
    "4890": {"operation": "bound_variable", "name": "native_O3_power"},
    "4891": {
        "operation": "sum",
        "arguments": [4584, 4890, 2084, 2091, 649, 3, 645, 1507],
    },
    "4892": {"operation": "exact_rational", "numerator": 1014, "denominator": 1},
    "4893": {
        "operation": "sum",
        "arguments": [4892, 2084, 2091, 649, 3, 645, 1507],
    },
    "4894": {"operation": "negative", "argument": 4585},
    "4914": {"operation": "negative", "argument": 4891},
    "4915": {"operation": "sum", "arguments": [4893, 4914]},
    "4916": {"operation": "negative", "argument": 4915},
    "4917": {"operation": "analytic_unary", "name": "exp", "argument": 4916},
}

EXPECTED_MEMORY_ROWS = {
    "m": {"memory": 4919, "negative_argument": 4918, "width": 4895},
    "e": {"memory": 4919, "negative_argument": 4918, "width": 4895},
    "h": {
        "memory": 4927,
        "negative_argument": 4926,
        "rate_scaled_width": 4925,
        "rate_parameter": 272,
        "width": 4895,
    },
    "k": {
        "memory": 4927,
        "negative_argument": 4926,
        "rate_scaled_width": 4925,
        "rate_parameter": 272,
        "width": 4895,
    },
    "p": {"memory": 1, "no_exp_factor": True},
}

EXPECTED_MEMORY_SUPPORT_ROWS = {
    4895: {"operation": "sum", "arguments": [4893, 4894]},
    4918: {"operation": "negative", "argument": 4895},
    4919: {"operation": "analytic_unary", "name": "exp", "argument": 4918},
    4925: {"operation": "product", "arguments": [272, 4895]},
    4926: {"operation": "negative", "argument": 4925},
    4927: {"operation": "analytic_unary", "name": "exp", "argument": 4926},
}

EXPECTED_MODULE_SHA256 = {
    "lei_ren_part1_paper_compliant_current_original_C3_source_targets.py": (
        "de41858de44ff9cea9b5d2b17aeca55eeca2c114969ca4df383ec1345b6badeb"
    ),
    "lei_ren_part1_paper_compliant_current_original_C3_repair_band.py": (
        "62e2fa8598bf307126be0ffe0b4e23de83152104524db03f61b4944c225132e4"
    ),
    "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection.py": (
        "aed09a89972ef9de0234e72691c4b5d33dfaa06d82d07b8df0781ef0f7a3807d"
    ),
    "lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N.py": (
        "83e5a8bdae75d20a261d24dd22d30cc7cb5e603e8fa2f0133958ca880b5e76db"
    ),
}

EXPECTED_AST_SHA256 = {
    (
        "lei_ren_part1_paper_compliant_current_original_C3_source_targets.py",
        "build_transport",
    ): "b921b16a0fdb86b3ac26460890c68dea984a56769ffd30835f8026f1f4f8f653",
    (
        "lei_ren_part1_paper_compliant_current_original_C3_repair_band.py",
        "build",
    ): "d229e0b601593275231c6d4f671309b43fbe9a2a92b8b648b2380bfb3a5b4f08",
    (
        "lei_ren_part1_paper_compliant_current_original_C3_repair_band.py",
        "project",
    ): "3f35ad3f8047afbbc564041eddc229c8bb76a82751bd6db2175718edcfc5ce5a",
    (
        "lei_ren_part1_paper_compliant_current_original_C3_repair_band.py",
        "endpoint_source",
    ): "83943428def433eabb15decbab9c3d555c3e55aa3d77138fec8c047113520dee",
    (
        "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection.py",
        "build_graph",
    ): "dca9af938dc9710635301a3d15b462209fc894c840a988a678f3281018349023",
    (
        "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection.py",
        "exact_geometry",
    ): "1508034d7714681948595e89f4f99d1f32db7c4c3a03b73772b974bdc1df498e",
    (
        "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection.py",
        "recipe_manifest",
    ): "e3c8627945af7de31faaecc674ee8c7154d0e991b16fa88317f7008e16944444",
    (
        "lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N.py",
        "source_mu_admission",
    ): "8507a12bb9de7e54ee103da586c0d50fa3e8fafa3e0a5d5034b0dde9e8dbb3f4",
}

REQUIRED_SOURCE_SNIPPETS = {
    (
        "lei_ren_part1_paper_compliant_current_original_C3_source_targets.py",
        "build_transport",
    ): (
        "if base['quiet_local_source_zero']:",
        "density={key:alg.fixed(0) for key in current.RATES};coefficients=None;F=None",
        "third=g.zero if q==g.zero else",
        "predecessor_ZZZ_memory_retained=True",
    ),
    (
        "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection.py",
        "build_graph",
    ): (
        "memory=g.unary('exp',g.neg(g.mul(g.constant(rate),width))) if rate else g.one",
        "measure='original dlogR: Jacobian exactly once'",
        "quiet_local_source_zero=chart=='O3_power'",
    ),
    (
        "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection.py",
        "exact_geometry",
    ): (
        "O3_power=Rm+18+sy.exp(40)+x",
        "O3_power=(sy.Integer(0),sy.Integer(2))",
    ),
    (
        "lei_ren_part1_paper_compliant_current_original_C3_repair_band.py",
        "build",
    ): (
        "invN=g.quotient(g.one,built['N'],'same selected exact positive integer N')",
        "actual_incoming_correction_C3=incoming",
        "partial_endpoint='repair_x in[1,2]'",
    ),
    (
        "lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N.py",
        "source_mu_admission",
    ): (
        "margin=logmu+c.ln(2)-eta_log",
        "Same-source O3 quiet power2mu>=eta not proved",
    ),
}

SOURCE_BASENAMES = {
    "c3_source": "lei_ren_part1_paper_compliant_current_original_C3_source_targets.py",
    "c3_repair": "lei_ren_part1_paper_compliant_current_original_C3_repair_band.py",
    "control": "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection.py",
    "o3": "lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N.py",
}


def _as_text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, str):
        return value
    return None


def _sha256_text(value):
    text = _as_text(value)
    if text is None:
        return None
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _basename(value):
    if not isinstance(value, str):
        return ""
    return value.replace("\\", "/").rsplit("/", 1)[-1]


def _record(records, basename):
    if not isinstance(records, dict):
        return None
    if basename in records:
        return records[basename]
    for key, value in records.items():
        if _basename(key) == basename:
            return value
    return None


def _cached_record(records, basename):
    if not isinstance(records, dict):
        return None
    if basename in records:
        return records[basename]
    for key, value in records.items():
        if _basename(key) == basename:
            return value
    return None


def qpc_journal(journal, event):
    if isinstance(journal, list):
        journal.append(event)


def _find_definition(tree, name):
    short_name = name.rsplit(".", 1)[-1]
    for item in tree.body:
        if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if item.name == short_name:
                return item
        if isinstance(item, ast.ClassDef):
            for member in item.body:
                if isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    if member.name == short_name:
                        return member
    return None


def _ast_hash(node):
    module = ast.Module(body=[node], type_ignores=[])
    dumped = ast.dump(module, include_attributes=False)
    return hashlib.sha256(dumped.encode("utf-8")).hexdigest()


def _get_receipt(spec):
    raw = spec.get("C3_receipt") if isinstance(spec, dict) else None
    if isinstance(raw, dict):
        return raw, None
    if isinstance(raw, str):
        try:
            return json.loads(raw), None
        except (TypeError, ValueError) as exc:
            return None, "C3_receipt string is not valid JSON: %s" % (exc,)
    return None, "C3_receipt v3 content is missing"


def _frozen_entry(spec):
    receipt, error = _get_receipt(spec)
    blockers = []
    if error:
        blockers.append(error)
        return None, blockers, False

    if receipt.get("receipt_version") != "v3":
        blockers.append("C3_receipt must be receipt_version v3")
    corrected = receipt.get("corrected_N_and_invN_literals")
    if not isinstance(corrected, dict):
        blockers.append("corrected_N_and_invN_literals is missing")
        corrected = {}
    if corrected.get("J_exact_string") != EXPECTED_J:
        blockers.append("canonical J is not the exact Python literal string")
    if corrected.get("N_node_id") != 18 or corrected.get("invN_node_id") != 19:
        blockers.append("N/invN node IDs are not 18/19")
    for key in ("rebuilt_node_18", "primary_node_18"):
        row = corrected.get(key)
        if not isinstance(row, dict) or row.get("J") != EXPECTED_J:
            blockers.append("%s does not preserve exact J string" % key)
    archive = receipt.get("rebuilt_archive_authorized_selection")
    if not isinstance(archive, dict):
        blockers.append("rebuilt_archive_authorized_selection is missing")
        archive = {}
    if archive.get("raw_sha256", "").lower() != EXPECTED_REBUILT_ARCHIVE_SHA256:
        blockers.append("rebuilt archive raw SHA-256 is not authenticated")
    if archive.get("exact_graph_nodes_count") != 24952:
        blockers.append("rebuilt archive node count is not 24952")

    slots = receipt.get("rebuilt_incoming_history_slots")
    if not isinstance(slots, dict):
        blockers.append("rebuilt_incoming_history_slots is missing")
        slots = {}
    if slots.get("incoming_ids_by_field") != EXPECTED_INCOMING:
        blockers.append("rebuilt incoming slot arrays do not match v3")
    bodies = slots.get("incoming_node_bodies")
    if not isinstance(bodies, dict):
        blockers.append("rebuilt incoming node bodies are missing")
    else:
        for field in EXPECTED_FIELDS:
            actual = bodies.get(field)
            expected_ids = EXPECTED_INCOMING[field]
            if not isinstance(actual, list):
                blockers.append("incoming node body list missing for %s" % field)
                continue
            by_id = {}
            for item in actual:
                if isinstance(item, dict) and isinstance(item.get("node"), int):
                    by_id[item["node"]] = item.get("row")
            for node_id in expected_ids:
                if node_id not in by_id:
                    blockers.append("incoming node %s is missing" % node_id)

    additional_root = receipt.get("additional_literal_rows")
    if not isinstance(additional_root, dict):
        blockers.append("additional_literal_rows is missing")
        additional_root = {}
    additional = additional_root.get("source_guard_and_offset_rows")
    if not isinstance(additional, dict):
        blockers.append("source_guard_and_offset_rows is missing")
        additional = {}
    for node_id, expected_row in EXPECTED_ADDITIONAL_ROWS.items():
        if additional.get(node_id) != expected_row:
            blockers.append("literal row %s does not match v3" % node_id)

    return receipt, blockers, not blockers


def _port_items(value):
    if isinstance(value, int) and not isinstance(value, bool):
        return [value]
    if isinstance(value, str) and value.isdigit():
        return [int(value)]
    if isinstance(value, dict):
        for key in (
            "right_source_cut",
            "right_source_cut_nodes",
            "nodes",
            "ports",
            "ids",
            "Z0_to_Z3",
            "Z3",
        ):
            if key in value:
                return _port_items(value[key])
        for key in ("node", "function_node", "id", "port"):
            if key in value:
                return _port_items(value[key])
        return []
    if isinstance(value, (list, tuple)):
        result = []
        for item in value:
            result.extend(_port_items(item))
        return result
    return []


def _extract_ports_raw(raw):
    if isinstance(raw, list):
        result = {field: [] for field in EXPECTED_FIELDS}
        for item in raw:
            if not isinstance(item, dict) or item.get("field") not in result:
                continue
            value = item.get("right_source_cut", item.get("node", item.get("id")))
            result[item["field"]].extend(_port_items(value))
        return result
    if not isinstance(raw, dict):
        return None
    result = {}
    for field in EXPECTED_FIELDS:
        value = raw.get(field)
        if value is None:
            for container_key in ("right_source_cut", "roles", "fields"):
                container = raw.get(container_key)
                if isinstance(container, dict) and field in container:
                    value = container[field]
                    break
                if isinstance(container, list):
                    values = []
                    for item in container:
                        if isinstance(item, dict) and item.get("field") == field:
                            values.extend(
                                _port_items(
                                    item.get(
                                        "right_source_cut",
                                        item.get("node", item.get("id")),
                                    )
                                )
                            )
                    if values:
                        value = values
                        break
        result[field] = _port_items(value)
    return result


def _extract_ports(parent_result, literal_slice=None):
    if isinstance(parent_result, dict):
        result = _extract_ports_raw(parent_result.get("signed_history_gap_function_nodes"))
        if result is not None:
            return result
    if isinstance(literal_slice, dict):
        for key in (
            "outgoing_ports",
            "terminal_ports",
            "signed_history_gap_function_nodes",
            "terminal_history_ports",
        ):
            result = _extract_ports_raw(literal_slice.get(key))
            if result is not None:
                return result
        result = _extract_ports_raw(literal_slice)
        if result is not None and any(result.get(field) for field in EXPECTED_FIELDS):
            return result
    return None


def _add_node_container(container, table):
    if isinstance(container, dict):
        for key, value in container.items():
            try:
                node_id = int(key)
            except (TypeError, ValueError):
                continue
            row = value.get("row") if isinstance(value, dict) else None
            if isinstance(row, dict):
                table[node_id] = row
            elif isinstance(value, dict) and "operation" in value:
                table[node_id] = value
    elif isinstance(container, list):
        for item in container:
            if not isinstance(item, dict):
                continue
            node_id = item.get("node", item.get("id"))
            if not isinstance(node_id, int):
                continue
            row = item.get("row")
            if isinstance(row, dict):
                table[node_id] = row
            elif "operation" in item:
                table[node_id] = item


def _extract_node_table(parent_result, literal_slice=None):
    table = {}
    if isinstance(parent_result, dict):
        for key in (
            "node_bodies",
            "selected_node_bodies",
            "terminal_node_bodies",
            "memory_node_bodies",
            "support_node_bodies",
            "nodes",
            "graph_nodes",
        ):
            _add_node_container(parent_result.get(key), table)
        raw = parent_result.get("signed_history_gap_function_nodes")
        if isinstance(raw, dict):
            for field in EXPECTED_FIELDS:
                _add_node_container(raw.get(field), table)
        _add_node_container(parent_result.get("right_source_cut_rows"), table)
    if isinstance(literal_slice, dict):
        for key in (
            "node_bodies",
            "selected_node_bodies",
            "terminal_node_bodies",
            "memory_node_bodies",
            "support_node_bodies",
            "geometry_node_bodies",
            "nodes",
        ):
            _add_node_container(literal_slice.get(key), table)
        incoming = literal_slice.get("incoming_node_bodies")
        if isinstance(incoming, dict):
            for field in EXPECTED_FIELDS:
                _add_node_container(incoming.get(field), table)
    return table


def _has_exact_keys(mapping, expected_keys):
    return isinstance(mapping, dict) and set(mapping) == set(expected_keys)


def _context_source_families(spec, parent_result, literal_slice):
    values = []
    if isinstance(literal_slice, dict) and literal_slice.get("source_family") is not None:
        values.append(("C3_literal_slice.source_family", literal_slice.get("source_family")))
    if isinstance(parent_result, dict):
        for key in ("source_family", "query_source_family", "parent_source_family"):
            if parent_result.get(key) is not None:
                values.append(("parent_result.%s" % key, parent_result.get(key)))
    if isinstance(spec, dict):
        for key in ("source_family", "query_source_family"):
            if spec.get(key) is not None:
                values.append(("spec.%s" % key, spec.get(key)))
    return values


def _check_archive_row4898(spec, parent_result, blockers):
    literal_slice = spec.get("C3_literal_slice") if isinstance(spec, dict) else None
    if not isinstance(literal_slice, dict):
        blockers.append("C3_literal_slice is required for direct archive row 4898 binding")
        return None, None
    if str(literal_slice.get("archive_raw_sha256", "")).lower() != EXPECTED_REBUILT_ARCHIVE_SHA256:
        blockers.append("C3_literal_slice archive raw SHA-256 is not authenticated")
    nodes = literal_slice.get("nodes")
    if not isinstance(nodes, dict):
        blockers.append("C3_literal_slice.nodes is missing")
        return None, None
    row = nodes.get("4898")
    if row is None:
        row = nodes.get(4898)
    if not isinstance(row, dict):
        blockers.append("C3_literal_slice.nodes['4898'] full row is missing")
        return None, None

    if not _has_exact_keys(row, EXPECTED_ARCHIVE_ROW4898_KEYS):
        blockers.append("archive node 4898 row keys are not the full authorized row")
    if row.get("Z_order") != 0:
        blockers.append("archive node 4898 Z_order is not 0")
    if row.get("Z_variable") != 20:
        blockers.append("archive node 4898 Z_variable is not 20")
    if row.get("coordinate") != 4890:
        blockers.append("archive node 4898 coordinate is not native O3 coordinate 4890")
    if row.get("defining_quantity_not_a_range_value") is not True:
        blockers.append("archive node 4898 defining quantity marker is missing")
    if row.get("native_chart") != "O3_power":
        blockers.append("archive node 4898 native chart is not O3_power")
    if row.get("operation") != "current_original_leading_function_recipe":
        blockers.append("archive node 4898 operation is not the leading recipe operation")
    if row.get("phase_independent_leading_source") is not True:
        blockers.append("archive node 4898 phase-independent source marker is missing")
    if row.get("quantity") != "E":
        blockers.append("archive node 4898 quantity is not E")

    recipe = row.get("recipe")
    if not isinstance(recipe, dict) or not _has_exact_keys(recipe, EXPECTED_ARCHIVE_RECIPE_KEYS):
        blockers.append("archive node 4898 recipe is not the full authorized recipe")
        recipe = recipe if isinstance(recipe, dict) else {}
    if recipe.get("AST_sha256") != "5c3b201139ff50c8cc56c3a99df9e2a8054292568cdd2a640364a1743eec7e54":
        blockers.append("archive node 4898 recipe AST SHA is not canonical")
    if recipe.get("defining_provider") != "WholeZAllNOuterRcFunctions.leading_packet":
        blockers.append("archive node 4898 recipe provider is not leading_packet")
    if recipe.get("function") != "WholeZAllNOuterRcFunctions.leading_packet":
        blockers.append("archive node 4898 recipe function is not leading_packet")
    if recipe.get("leading_source_fixed_input_N0") is not True:
        blockers.append("archive node 4898 fixed-input N0 marker is missing")
    if recipe.get("module") != "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions.py":
        blockers.append("archive node 4898 recipe module is not the canonical outer module")
    if recipe.get("module_sha256") != "9d73659a1c76163dddee6653af53f56ebffd7045986d871dc559ba1e053620b7":
        blockers.append("archive node 4898 recipe module SHA is not canonical")
    if recipe.get("native_chart") != "O3_power":
        blockers.append("archive node 4898 recipe native chart is not O3_power")
    if recipe.get("original_independent_P0_and_parent_leading_recipes_retained") is not True:
        blockers.append("archive node 4898 parent leading recipes retention marker is missing")
    if recipe.get("range_endpoints_not_field_coefficients") is not True:
        blockers.append("archive node 4898 range endpoint marker is missing")

    quantity_paths = recipe.get("quantity_paths")
    expected_paths = {
        "Delta": "original_full_source_quotients.actual_Delta_axial5[Z_order]",
        "E": "original_generic_source.common_velocity_E_axial5[Z_order]",
        "V": "original_generic_source.common_velocity_V_axial5[Z_order]",
        "a": "original_roots.a[(0,Z_order)]",
        "b": "original_full_source_quotients.actual_b_axial5[Z_order]",
        "p2": "original_roots.p2[(0,Z_order)]",
        "t0": "original_roots.t0[(0,Z_order)]",
    }
    if not _has_exact_keys(quantity_paths, EXPECTED_ARCHIVE_QUANTITY_PATH_KEYS):
        blockers.append("archive node 4898 quantity paths are not complete")
    elif quantity_paths != expected_paths:
        blockers.append("archive node 4898 quantity paths do not match canonical source")

    proof = recipe.get("quiet_source_proof")
    if not isinstance(proof, dict) or not _has_exact_keys(proof, EXPECTED_ARCHIVE_PROOF_KEYS):
        blockers.append("archive node 4898 quiet proof is not the full direct proof row")
        proof = proof if isinstance(proof, dict) else None
    if isinstance(proof, dict):
        admission = proof.get("admission")
        if not _has_exact_keys(admission, EXPECTED_ARCHIVE_ADMISSION_KEYS):
            blockers.append("archive node 4898 admission binding is incomplete")
        frontend = proof.get("source_frontend_binding")
        if not _has_exact_keys(frontend, EXPECTED_ARCHIVE_FRONTEND_KEYS):
            blockers.append("archive node 4898 frontend binding is incomplete")
        if proof.get("source_identity") != "2*mu>=eta implies q=q_Z=0 throughout current power offset [0,2]":
            blockers.append("archive node 4898 quiet proof identity is not canonical")
        if proof.get("original_primitive_A_B_zero_but_incoming_histories_retained") is not True:
            blockers.append("archive node 4898 quiet proof does not retain incoming histories")
        if proof.get("checked_receipt") != "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions_check.json":
            blockers.append("archive node 4898 checked receipt name is not canonical")
        if proof.get("checked_receipt_sha256") != "f3c4d1b230b3835da6f4f270926fa85f74ca488f8331e05fbbb894a7c22a0158":
            blockers.append("archive node 4898 checked receipt SHA is not canonical")
        if isinstance(admission, dict):
            if admission.get("AST_sha256") != EXPECTED_AST_SHA256[(SOURCE_BASENAMES["o3"], "source_mu_admission")]:
                blockers.append("archive node 4898 admission AST SHA is not canonical")
            if admission.get("function") != "source_mu_admission":
                blockers.append("archive node 4898 admission function is not source_mu_admission")
            if admission.get("module") != SOURCE_BASENAMES["o3"]:
                blockers.append("archive node 4898 admission module is not O3 source module")
            if admission.get("module_sha256") != EXPECTED_MODULE_SHA256[SOURCE_BASENAMES["o3"]]:
                blockers.append("archive node 4898 admission module SHA is not canonical")
        if isinstance(frontend, dict):
            if frontend.get("exact_unchanged_leading_frontend_AST_sha256") != "f58af2bc7d3492b5d5cfbc30d0653be5d6cb1f794fc19c336879e788bd2f9e56":
                blockers.append("archive node 4898 frontend AST SHA is not canonical")
            if frontend.get("original_parent_leading_source_recipe_not_changed") is not True:
                blockers.append("archive node 4898 parent leading recipe marker is missing")
            if frontend.get("original_query_module") != SOURCE_BASENAMES["o3"]:
                blockers.append("archive node 4898 frontend module is not O3 source module")
            if frontend.get("original_query_sha256") != EXPECTED_MODULE_SHA256[SOURCE_BASENAMES["o3"]]:
                blockers.append("archive node 4898 frontend module SHA is not canonical")
            if frontend.get("stop_before_all_u_support_and_fixed_N_density") is not True:
                blockers.append("archive node 4898 frontend stop marker is missing")

    if row.get("source_family") != EXPECTED_ARCHIVE_SOURCE_FAMILY:
        blockers.append("archive node 4898 source_family is not the authenticated family")
    context_families = _context_source_families(spec, parent_result, literal_slice)
    for context_name, context_family in context_families:
        if context_family != EXPECTED_ARCHIVE_SOURCE_FAMILY:
            blockers.append("%s does not match archive node 4898 source_family" % context_name)
        if context_name != "C3_literal_slice.source_family" and context_family != row.get("source_family"):
            blockers.append("%s does not match archive node 4898 row source_family" % context_name)
    return row, proof


def _check_memory_rows(node_table, blockers):
    checks = []
    for node_id, expected in EXPECTED_MEMORY_SUPPORT_ROWS.items():
        observed = node_table.get(node_id)
        passed = observed == expected
        checks.append(
            {
                "node": node_id,
                "expected": expected,
                "observed": observed,
                "passed": passed,
            }
        )
        if not passed:
            blockers.append(
                "memory support node %s does not match expected literal row" % node_id
            )
    return checks


def _check_terminal_rows(node_table, ports, blockers):
    checks = []
    if ports is None:
        blockers.append("parent_result.signed_history_gap_function_nodes is missing")
        return checks
    for field in EXPECTED_FIELDS:
        if ports.get(field) != EXPECTED_TERMINALS[field]:
            blockers.append("signed history ports do not match rebuilt C3 %s roles" % field)
        for node_id in EXPECTED_TERMINALS[field]:
            row = node_table.get(node_id)
            expected = EXPECTED_TERMINAL_ROWS[node_id]
            passed = row == expected
            checks.append(
                {
                    "field": field,
                    "node": node_id,
                    "expected": expected,
                    "observed": row,
                    "passed": passed,
                }
            )
            if not passed:
                blockers.append("terminal node %s does not match expected literal row" % node_id)
    return checks


def _check_source(records, basename, function_hashes, snippets, blockers):
    record = _record(records, basename)
    if not isinstance(record, dict):
        blockers.append("creator source record is missing: %s" % basename)
        return {"basename": basename, "passed": False}
    source_text = _as_text(record.get("source_text"))
    if source_text is None:
        blockers.append("creator source text is missing: %s" % basename)
        return {"basename": basename, "passed": False}
    observed_module_sha = _sha256_text(source_text)
    declared_sha = str(record.get("sha256", "")).lower()
    passed = True
    if declared_sha != observed_module_sha:
        blockers.append("source text SHA mismatch: %s" % basename)
        passed = False
    expected_module_sha = EXPECTED_MODULE_SHA256.get(basename)
    if expected_module_sha and observed_module_sha != expected_module_sha:
        blockers.append("source module SHA is not the authenticated C3 source: %s" % basename)
        passed = False
    try:
        tree = ast.parse(source_text, filename=basename)
    except SyntaxError as exc:
        blockers.append("source AST parse failed for %s: %s" % (basename, exc))
        return {"basename": basename, "module_sha256": observed_module_sha, "passed": False}

    functions = {}
    for function_name, expected_ast_sha in function_hashes.items():
        node = _find_definition(tree, function_name)
        if node is None:
            blockers.append("source function is missing: %s:%s" % (basename, function_name))
            passed = False
            continue
        observed_ast_sha = _ast_hash(node)
        functions[function_name] = {
            "span": [node.lineno, node.end_lineno],
            "ast_sha256": observed_ast_sha,
            "expected_ast_sha256": expected_ast_sha,
            "passed": observed_ast_sha == expected_ast_sha,
        }
        if observed_ast_sha != expected_ast_sha:
            blockers.append("source AST hash mismatch: %s:%s" % (basename, function_name))
            passed = False
        for snippet in snippets.get(function_name, ()):
            if snippet not in source_text:
                blockers.append(
                    "required source literal is missing: %s:%s:%s"
                    % (basename, function_name, snippet)
                )
                passed = False
    return {
        "basename": basename,
        "path": record.get("path"),
        "module_sha256": observed_module_sha,
        "functions": functions,
        "passed": passed,
    }


def _static_source_checks(records, blockers):
    by_function = {}
    for basename, names in (
        (
            SOURCE_BASENAMES["c3_source"],
            {"build_transport": EXPECTED_AST_SHA256[(SOURCE_BASENAMES["c3_source"], "build_transport")]},
        ),
        (
            SOURCE_BASENAMES["c3_repair"],
            {
                "build": EXPECTED_AST_SHA256[(SOURCE_BASENAMES["c3_repair"], "build")],
                "project": EXPECTED_AST_SHA256[(SOURCE_BASENAMES["c3_repair"], "project")],
                "endpoint_source": EXPECTED_AST_SHA256[
                    (SOURCE_BASENAMES["c3_repair"], "endpoint_source")
                ],
            },
        ),
        (
            SOURCE_BASENAMES["control"],
            {
                "build_graph": EXPECTED_AST_SHA256[(SOURCE_BASENAMES["control"], "build_graph")],
                "exact_geometry": EXPECTED_AST_SHA256[
                    (SOURCE_BASENAMES["control"], "exact_geometry")
                ],
                "recipe_manifest": EXPECTED_AST_SHA256[
                    (SOURCE_BASENAMES["control"], "recipe_manifest")
                ],
            },
        ),
        (
            SOURCE_BASENAMES["o3"],
            {
                "source_mu_admission": EXPECTED_AST_SHA256[
                    (SOURCE_BASENAMES["o3"], "source_mu_admission")
                ],
            },
        ),
    ):
        snippets = {}
        for (snippet_basename, function_name), values in REQUIRED_SOURCE_SNIPPETS.items():
            if snippet_basename == basename:
                snippets[function_name] = values
        by_function[basename] = _check_source(
            records, basename, names, snippets, blockers
        )

    optional_outer = "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions.py"
    outer_record = _record(records, optional_outer)
    if outer_record is not None:
        by_function[optional_outer] = _check_source(
            records,
            optional_outer,
            {"leading_packet": "5c3b201139ff50c8cc56c3a99df9e2a8054292568cdd2a640364a1743eec7e54"},
            {},
            blockers,
        )
    quiet_constructor_bound = (
        by_function.get(SOURCE_BASENAMES["c3_source"], {}).get("passed", False)
        and by_function.get(SOURCE_BASENAMES["control"], {}).get("passed", False)
    )
    return by_function, quiet_constructor_bound


def _locate_quiet_proof(spec, receipt):
    # The proof must come directly from the authenticated archive row.  Any
    # receipt copy is provenance only and is never used as the bound value.
    literal_slice = spec.get("C3_literal_slice") if isinstance(spec, dict) else None
    nodes = literal_slice.get("nodes") if isinstance(literal_slice, dict) else None
    row = nodes.get("4898") if isinstance(nodes, dict) else None
    if row is None and isinstance(nodes, dict):
        row = nodes.get(4898)
    recipe = row.get("recipe") if isinstance(row, dict) else None
    proof = recipe.get("quiet_source_proof") if isinstance(recipe, dict) else None
    if isinstance(proof, dict):
        return (
            proof,
            EXPECTED_ARCHIVE_QUIET_SELECTOR,
            EXPECTED_RECEIPT_PROOF_PROVENANCE_SELECTOR,
        )
    return None, None, EXPECTED_RECEIPT_PROOF_PROVENANCE_SELECTOR


def _gate_status(content):
    try:
        raw_text = _as_text(content)
    except UnicodeDecodeError:
        return {
            "inspectable": False,
            "passed": False,
            "failures": ["content_not_utf8"],
        }
    if raw_text is None:
        return {
            "inspectable": False,
            "passed": False,
            "failures": ["content_not_raw_utf8_json_text"],
        }
    try:
        content = json.loads(raw_text, parse_int=str)
    except (TypeError, ValueError):
        return {
            "inspectable": False,
            "passed": False,
            "failures": ["raw_content_not_json"],
        }
    if not isinstance(content, dict):
        return {
            "inspectable": False,
            "passed": False,
            "failures": ["cached_content_not_object"],
        }
    failures = []
    if content.get("all_passed") is not True:
        failures.append("all_passed")
    for key in ("family_pass", "source_family_pass"):
        if key in content and content[key] is not True:
            failures.append(key)
    source_family = content.get("source_family")
    if isinstance(source_family, dict):
        for key in ("passed", "pass", "all_passed", "family_pass", "source_family_pass"):
            if key in source_family and source_family[key] is not True:
                failures.append("source_family.%s" % key)
    gates = content.get("gates")
    if isinstance(gates, dict):
        for gate_name, gate_value in gates.items():
            if isinstance(gate_value, bool) and not gate_value:
                failures.append("gates.%s" % gate_name)
            if isinstance(gate_value, dict):
                for key in ("passed", "all_passed", "family_pass", "source_family_pass"):
                    if key in gate_value and gate_value[key] is not True:
                        failures.append("gates.%s.%s" % (gate_name, key))
    elif isinstance(gates, list):
        for index, gate_value in enumerate(gates):
            if isinstance(gate_value, bool) and not gate_value:
                failures.append("gates[%s]" % index)
            if isinstance(gate_value, dict):
                for key in ("passed", "all_passed", "family_pass", "source_family_pass"):
                    if key in gate_value and gate_value[key] is not True:
                        failures.append("gates[%s].%s" % (index, key))
    return {
        "inspectable": True,
        "passed": not failures,
        "failures": failures,
        "parsed_content": content,
    }


def _check_cached_receipt(proof, cached, blockers):
    checked_name = _basename(proof.get("checked_receipt"))
    expected_sha = str(proof.get("checked_receipt_sha256", "")).lower()
    item = _cached_record(cached, checked_name)
    if not isinstance(item, dict):
        blockers.append("cached admission receipt is missing: %s" % checked_name)
        return {"basename": checked_name, "passed": False}
    declared_sha = str(item.get("sha256", "")).lower()
    if not expected_sha:
        blockers.append("cached admission receipt expected SHA is missing: %s" % checked_name)
    if declared_sha != expected_sha:
        blockers.append("cached admission receipt SHA mismatch: %s" % checked_name)
    content = item.get("content")
    raw_sha = _sha256_text(content)
    if raw_sha is None:
        blockers.append("cached admission receipt raw UTF-8 JSON content is missing: %s" % checked_name)
    elif raw_sha != expected_sha:
        blockers.append("cached admission receipt raw content SHA mismatch: %s" % checked_name)
    gates = _gate_status(content)
    if not gates["passed"]:
        blockers.append("cached admission receipt gates failed: %s" % checked_name)
    return {
        "basename": checked_name,
        "declared_sha256": declared_sha,
        "expected_sha256": expected_sha,
        "raw_content_sha256": raw_sha,
        "gates": gates,
        "parsed_content": gates.get("parsed_content"),
        "passed": (
            declared_sha == expected_sha
            and raw_sha == expected_sha
            and gates["passed"]
        ),
    }


def _check_quiet_proof(spec, receipt, source_records, cached_records, source_info, blockers):
    proof, selector_path, receipt_provenance_path = _locate_quiet_proof(spec, receipt)
    if proof is None:
        blockers.append("direct archive row 4898 quiet-source proof is missing before proof checks")
        return None, {
            "checked": False,
            "passed": False,
            "json_selector_path": None,
            "archive_quiet_source_proof_selector_path": None,
            "receipt_provenance_selector_path": receipt_provenance_path,
            "actual_value": None,
        }
    if proof.get("original_primitive_A_B_zero_but_incoming_histories_retained") is not True:
        blockers.append("quiet proof does not retain incoming histories explicitly")
    identity = proof.get("source_identity", "")
    if "2*mu>=eta" not in identity or "q=q_Z=0" not in identity:
        blockers.append("quiet proof source identity is not the O3 quiet identity")
    admission = proof.get("admission")
    if not isinstance(admission, dict):
        blockers.append("quiet proof source_mu_admission binding is missing")
        admission = {}
    admission_name = _basename(admission.get("module"))
    if admission_name != SOURCE_BASENAMES["o3"]:
        blockers.append("quiet proof admission module is not source_mu_admission module")
    expected_admission_sha = EXPECTED_MODULE_SHA256[SOURCE_BASENAMES["o3"]]
    if admission.get("module_sha256") != expected_admission_sha:
        blockers.append("quiet proof admission module SHA is not canonical")
    if admission.get("function") != "source_mu_admission":
        blockers.append("quiet proof admission function is not source_mu_admission")
    expected_admission_ast = EXPECTED_AST_SHA256[
        (SOURCE_BASENAMES["o3"], "source_mu_admission")
    ]
    if admission.get("AST_sha256") != expected_admission_ast:
        blockers.append("quiet proof admission AST SHA is not canonical")

    front = proof.get("source_frontend_binding")
    if isinstance(front, dict):
        front_name = _basename(front.get("original_query_module"))
        if front_name and front_name != SOURCE_BASENAMES["o3"]:
            blockers.append("quiet proof frontend module is not O3 source module")
        front_sha = front.get("original_query_sha256")
        if front_sha and front_sha != expected_admission_sha:
            blockers.append("quiet proof frontend SHA is not canonical")

    row_family = None
    if isinstance(spec, dict):
        literal_slice = spec.get("C3_literal_slice")
        if isinstance(literal_slice, dict):
            row_family = literal_slice.get("source_family")
    proof_family = proof.get("source_family")
    if proof_family is not None and proof_family != row_family:
        blockers.append("quiet proof source_family does not match archive row 4898")
    family = row_family if row_family is not None else proof_family
    if family is None and isinstance(spec, dict):
        family = spec.get("source_family")
    family_check = {"present": isinstance(family, dict), "passed": True}
    if not isinstance(family, dict):
        blockers.append("archive row 4898 source_family is missing from literal slice")
    cached_result = _check_cached_receipt(proof, cached_records, blockers)
    cached_content = cached_result.get("parsed_content")
    if isinstance(family, dict) and isinstance(cached_content, dict):
        cached_family = cached_content.get("source_family")
        family_check["matched"] = cached_family == family
        if cached_family != family:
            blockers.append("cached receipt source_family does not match archive row 4898 family")
    elif isinstance(family, dict):
        family_check["matched"] = False
        blockers.append("cached receipt source_family object is missing")
    return proof, {
        "checked": True,
        "passed": not blockers,
        "json_selector_path": selector_path,
        "archive_quiet_source_proof_selector_path": selector_path,
        "receipt_provenance_selector_path": receipt_provenance_path,
        "actual_value": proof,
        "archive_quiet_source_proof": proof,
        "admission": admission,
        "checked_receipt": cached_result,
        "source_family": family,
        "source_family_check": family_check,
    }


def _make_result(
    spec,
    receipt,
    blockers,
    phase,
    source_info=None,
    quiet_constructor_bound=False,
    terminal_checks=None,
    memory_checks=None,
    outgoing_ports=None,
    proof_info=None,
    source_guard_difference=None,
    journal_events=None,
):
    source_info = source_info or {}
    terminal_checks = terminal_checks or []
    memory_checks = memory_checks or []
    proof_info = proof_info or {"checked": False, "passed": False}
    parent_result = spec.get("parent_result") if isinstance(spec, dict) else None
    if not isinstance(parent_result, dict):
        parent_result = {}
    native_n_flag = parent_result.get("native_N_scaling_applied_once_bound")
    source_gaps = parent_result.get("source_gaps", spec.get("source_gaps", {}) if isinstance(spec, dict) else {})
    receipt_slots = (
        receipt.get("rebuilt_incoming_history_slots", {}) if isinstance(receipt, dict) else {}
    )
    incoming_ports = receipt_slots.get("incoming_ids_by_field")
    literal_slice = spec.get("C3_literal_slice") if isinstance(spec, dict) else None
    archive_nodes = literal_slice.get("nodes") if isinstance(literal_slice, dict) else None
    archive_row_4898 = None
    if isinstance(archive_nodes, dict):
        archive_row_4898 = archive_nodes.get("4898", archive_nodes.get(4898))
    zero_rows = 20 if quiet_constructor_bound else 0
    local_rows = (
        {field: [0, 0, 0, 0] for field in EXPECTED_FIELDS}
        if quiet_constructor_bound
        else {}
    )
    result = {
        "passed": not blockers,
        "status": "passed" if not blockers else "blocked",
        "phase": phase,
        "unexecuted_draft": False,
        "new_source_binder_entered": True,
        "frozen_entry_guard_keys": list(FROZEN_ENTRY_GUARD_KEYS),
        "draft_sha256": spec.get("draft_sha256") if isinstance(spec, dict) else None,
        "zero_row_count": zero_rows,
        "local_density_rows": local_rows,
        "geometry": {
            "chart": "O3_power",
            "domain": [0, 2],
            "offset_source_literal": "Rm+18+sy.exp(40)+x",
            "measure": "original dlogR: Jacobian exactly once",
            "jacobian_applied_once": True,
            "absolute_Rm_Rc_crossarchive_alias": False,
        },
        "rates_normalized": {
            "m": "lambda_m",
            "h": "lambda_h",
            "k": "lambda_k",
            "e": "lambda_e",
            "p": "lambda_p",
        },
        "actual_memory": EXPECTED_MEMORY_ROWS,
        "memory_source_checks": {
            "memory_4919": {
                "negative_argument": 4918,
                "width_row": 4895,
                "fields": ["m", "e"],
            },
            "memory_4927": {
                "negative_argument": 4926,
                "rate_scaled_width": 4925,
                "rate_parameter": 272,
                "width_row": 4895,
                "fields": ["h", "k"],
            },
            "p_no_exp_factor": {"memory_node": 1, "fields": ["p"]},
            "suffix": {
                "native_variable": 4890,
                "right_native_origin_sum": 4891,
                "negative": 4914,
                "suffix": 4915,
                "negative_suffix": 4916,
                "exp": 4917,
            },
        },
        "incoming_ports": incoming_ports,
        "outgoing_ports": outgoing_ports,
        "archive_row_4898_selector": EXPECTED_ARCHIVE_ROW4898_SELECTOR,
        "archive_row_4898": archive_row_4898,
        "terminal_row_checks": terminal_checks,
        "memory_row_checks": memory_checks,
        "source_checks": source_info,
        "quiet_proof": proof_info,
        "archive_quiet_source_proof_selector_path": proof_info.get(
            "json_selector_path"
        ),
        "archive_quiet_source_proof_receipt_provenance_selector_path": proof_info.get(
            "receipt_provenance_selector_path"
        ),
        "archive_quiet_source_proof": proof_info.get("actual_value"),
        "native_N_scaling_applied_once_bound": native_n_flag,
        "native_N_scaling_reproof": False,
        "source_gaps": source_gaps,
        "entry_gap_status": (
            spec.get("entry_gap_status", "unknown")
            if isinstance(spec, dict)
            else "unknown"
        ),
        "entry_gap_zero_or_nonzero_asserted": False,
        "scope_primary_alias_false": True,
        "full17": False,
        "physical": False,
        "allcaps_untouched": True,
        "journal_events": journal_events or [],
        "blockers": list(blockers),
    }
    if source_guard_difference is not None:
        result["source_guard_difference"] = source_guard_difference
    return result


def bind_quiet_C3_source(spec, journal=None):
    """Return a source binding result under the external frozen-entry workflow.

    The function accepts only metadata and source text records. Source text is
    hashed and parsed with stdlib AST after the frozen receipt entry check; it
    is never imported, compiled for execution, or called.
    """
    blockers = []
    journal_events = []
    if not isinstance(spec, dict):
        return _make_result(
            {},
            None,
            ["spec must be a mapping"],
            "input",
            journal_events=journal_events,
        )

    receipt, frozen_blockers, frozen = _frozen_entry(spec)
    blockers.extend(frozen_blockers)
    qpc_journal(journal, {"stage": "frozen_entry", "passed": frozen, "blockers": list(blockers)})
    journal_events.append({"stage": "frozen_entry", "passed": frozen})
    if not frozen:
        return _make_result(
            spec,
            receipt,
            blockers,
            "frozen_entry",
            journal_events=journal_events,
        )

    # This is the only bounded source-guard arithmetic, after frozen entry.
    additional = receipt["additional_literal_rows"]["source_guard_and_offset_rows"]
    source_guard_difference = (
        int(additional["4892"]["numerator"])
        - int(additional["4584"]["numerator"])
    )
    if source_guard_difference != 2:
        blockers.append("source guard literal difference is not 2")
    qpc_journal(
        journal,
        {
            "stage": "bounded_source_guard",
            "passed": source_guard_difference == 2,
            "difference": source_guard_difference,
        },
    )
    journal_events.append(
        {
            "stage": "bounded_source_guard",
            "passed": source_guard_difference == 2,
        }
    )

    parent_result = spec.get("parent_result", {})
    if not isinstance(parent_result, dict):
        parent_result = {}
    literal_slice = spec.get("C3_literal_slice")
    _check_archive_row4898(spec, parent_result, blockers)
    outgoing_ports = _extract_ports(parent_result, literal_slice)
    node_table = _extract_node_table(parent_result, literal_slice)
    memory_checks = _check_memory_rows(node_table, blockers)
    terminal_checks = _check_terminal_rows(node_table, outgoing_ports, blockers)
    if blockers:
        return _make_result(
            spec,
            receipt,
            blockers,
            "bounded_rows",
            terminal_checks=terminal_checks,
            memory_checks=memory_checks,
            outgoing_ports=outgoing_ports,
            source_guard_difference=source_guard_difference,
            journal_events=journal_events,
        )

    records = spec.get("creator_source_records")
    source_info, quiet_constructor_bound = _static_source_checks(records, blockers)
    qpc_journal(
        journal,
        {
            "stage": "source_ast_and_hash",
            "passed": not blockers,
            "source_count": len(source_info),
        },
    )
    journal_events.append(
        {
            "stage": "source_ast_and_hash",
            "passed": not blockers,
        }
    )
    if blockers:
        return _make_result(
            spec,
            receipt,
            blockers,
            "source_ast",
            source_info=source_info,
            quiet_constructor_bound=quiet_constructor_bound,
            terminal_checks=terminal_checks,
            memory_checks=memory_checks,
            outgoing_ports=outgoing_ports,
            source_guard_difference=source_guard_difference,
            journal_events=journal_events,
        )

    cached = spec.get("cached_admission_receipts")
    proof, proof_info = _check_quiet_proof(
        spec,
        receipt,
        records,
        cached,
        source_info,
        blockers,
    )
    qpc_journal(
        journal,
        {
            "stage": "cached_admission_metadata",
            "passed": not blockers,
            "proof_checked": proof is not None,
        },
    )
    journal_events.append(
        {
            "stage": "cached_admission_metadata",
            "passed": not blockers,
            "proof_checked": proof is not None,
        }
    )

    native_flag = parent_result.get("native_N_scaling_applied_once_bound")
    if native_flag is None:
        blockers.append("parent native_N_scaling_applied_once_bound flag is missing")
    elif native_flag is not True:
        blockers.append("parent native_N_scaling_applied_once_bound flag is not true")

    return _make_result(
        spec,
        receipt,
        blockers,
        "complete",
        source_info=source_info,
        quiet_constructor_bound=quiet_constructor_bound,
        terminal_checks=terminal_checks,
        memory_checks=memory_checks,
        outgoing_ports=outgoing_ports,
        proof_info=proof_info,
        source_guard_difference=source_guard_difference,
        journal_events=journal_events,
    )


__all__ = ["bind_quiet_C3_source"]







"""Unexecuted native O3_power binder draft.

This module is deliberately data-only.  It validates supplied source records and
literal native receipts, parses supplied source text only as AST, and returns
conditional source-gate metadata.  It never imports or executes project code.
"""
import ast as qpn_ast
import hashlib as qpn_hashlib
from collections.abc import Mapping as qpn_Mapping

qpn_NATIVE_RECEIPT_SHA256 = "c963f41b44901683c99dd76e6d26e9d58d360dff54d93941ea5722b773bf8cec"
qpn_LITERAL_ADDENDUM_SHA256 = "3978472d568f94fab5ed852580e404e74d32b27a7db874b03f793fe46efe99fc"
qpn_PARENT_RESULT_SHA256 = "b91da155898c38a9ca4b26595983a38aa04b2965ae3db0a908fadcea97b34241"
qpn_NATIVE_ARCHIVE_SHA256 = "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d"
qpn_OUTER_SOURCE_SHA256 = "4464A934D303A95A0DA157E9E289F6ED8A04D5FF0ADCB0D663857F0D9C57DCD3".lower()
qpn_LOOP_SOURCE_SHA256 = "69D6391120000488472B7961B4BFF1901BFF872DEF95D183C4EE056B0948E17D".lower()
qpn_ADMISSION_MODULE = "lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N.py"
qpn_ADMISSION_MODULE_SHA256 = "83e5a8bdae75d20a261d24dd22d30cc7cb5e603e8fa2f0133958ca880b5e76db"
qpn_ADMISSION_AST_SHA256 = "8507a12bb9de7e54ee103da586c0d50fa3e8fafa3e0a5d5034b0dde9e8dbb3f4"
qpn_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}
qpn_RATES = {"m": "1", "h": "3/2", "k": "3/2", "e": "1", "p": "0"}
qpn_INCOMING_FIRST4 = {
    "m": [250004, 250005, 250006, 250007],
    "h": [250010, 250011, 250012, 250013],
    "k": [250016, 250017, 250018, 250019],
    "e": [250022, 250023, 250024, 250025],
    "p": [250028, 250029, 250030, 250031],
}
qpn_MEMORY = {"m": 250297, "h": 250325, "k": 250325, "e": 250297, "p": 1}
qpn_E_ROOT_FIRST4 = [250226, 250227, 250228, 250229]
qpn_A_ROOT_FIRST4 = [251146, 251150, 251156, 251163]
qpn_A_OVER_N_FIRST4 = [252995, 252996, 252997, 252998]
qpn_DELTA_E_ROWS = [253022, 253025, 253029, 253034]
qpn_DELTA_V_FIRST4 = [253048, 253049, 253050, 253051]
qpn_FINITE_HISTORY = {
    "m": [253257, 253258, 253259, 253260],
    "h": [253290, 253291, 253292, 253293],
    "k": [253320, 253321, 253322, 253323],
    "e": [253350, 253351, 253352, 253353],
    "p": [253374, 253375, 253376, 253377],
}
qpn_LOCAL_HISTORY = {
    "m": [253245, 253246, 253247, 253248],
    "h": [253278, 253279, 253280, 253281],
    "k": [253308, 253309, 253310, 253311],
    "e": [253338, 253339, 253340, 253341],
    "p": [253368, 253369, 253370, 253371],
}
qpn_ENDPOINT_ROWS = [
    258646, 258647, 258648, 258649,
    258652, 258653, 258654, 258655,
    258658, 258659, 258660, 258661,
    258664, 258665, 258666, 258667,
    258670, 258671, 258672, 258673,
]
qpn_FLAT_A_ROWS = [252768, 252769, 252770, 252771]
qpn_FLAT_B_ROWS = [252806, 252807, 252808, 252809]
qpn_SOURCE_TOKEN_GROUPS = {
    "outer_geometry": (
        "'O3_power':g.add(c(18),Y,x)",
        "else g.one",
        "'O3_power':(g.zero,c(2))",
        "decay=g.unary('exp',g.mul(g.constant(-rate),width))",
    ),
    "outer_power_context": (
        "E=a.scale(self.background_endpoints['O3_transition']['fields']['E']",
        "g.constant('-1/2')",
        "g.add(g.constant('1/2'),self.mu)",
        "Ey=a.scale(E,g.neg(g.add(g.constant('1/2'),g.mul(self.mu,sigma(g,x)))))",
        "Ey=a.scale(E,g.neg(g.add(g.constant('1/2'),self.mu)))",
        "Vy=a.fixed(0)",
    ),
    "prescribed_loop": (
        "RATES=dict(m=Fraction(1),h=Fraction(3,2),k=Fraction(3,2),e=Fraction(1),p=Fraction(0))",
        "def generic_roots(g,chart,Pstar):",
        "C=a.sub(E,a.scale(a.trim(radial['E_y']),2))",
        "B=a.scale(a.trim(radial['V_y']),2)",
        "return dict(E=E,V=V,a=aa,b=bb,p1=p1,p2=p2,t0=t0,Delta=a.sub(kap,a.fixed(2)),kappa=kap)",
        "invN=g.quotient(g.one,N,'same exact positive integer N on every source chart')",
        "x=a.scale(loop['A'],invN)",
        "exm1=FunctionJet((g.unary('expm1',x.value),)+rows(ex)[1:])",
        "r=a.mul(roots['E'],exm1)",
        "s=a.scale(loop['B_over_Pstar'],invN)",
        "inc=dict(m=s,h=r,k=a.add(a.mul(V,r),a.mul(E,s),a.mul(r,s))",
        "e=a.sub(a.add(a.scale(a.mul(V,s),2),a.mul(s,s)),a.add(a.mul(E,r),a.scale(a.mul(r,r),'1/2')))",
        "p=a.add(a.mul(E,r),a.scale(a.mul(r,r),'1/2'))",
    ),
}


def qpn_is_mapping(qpn_value):
    return isinstance(qpn_value, qpn_Mapping)


def qpn_receipt_envelope(qpn_value):
    qpn_payload = None
    qpn_declared_sha = None
    qpn_shape = "invalid"
    if qpn_is_mapping(qpn_value):
        qpn_declared = None
        for qpn_key in ("sha256", "raw_sha256", "_raw_sha256"):
            qpn_candidate = qpn_value.get(qpn_key)
            if isinstance(qpn_candidate, str):
                qpn_declared = qpn_candidate
                break
        if "payload" in qpn_value:
            qpn_payload = qpn_value.get("payload")
            qpn_shape = "metadata_wrapper"
        else:
            qpn_payload = qpn_value
            qpn_shape = "payload_with_raw_metadata"
        if isinstance(qpn_declared, str):
            qpn_declared_sha = qpn_declared.lower()
    return {
        "payload": qpn_payload if qpn_is_mapping(qpn_payload) else None,
        "declared_sha256": qpn_declared_sha,
        "shape": qpn_shape,
        "has_declared_sha256": qpn_declared_sha is not None,
        "payload_mapping": qpn_is_mapping(qpn_payload),
    }


def qpn_sha_from_value(qpn_value):
    return qpn_receipt_envelope(qpn_value)["declared_sha256"]


def qpn_norm_path(qpn_path):
    if not isinstance(qpn_path, str):
        return ""
    return qpn_path.replace("\\", "/").lower()


def qpn_record_for_suffix(qpn_records, qpn_suffix):
    if not qpn_is_mapping(qpn_records):
        return None
    qpn_suffix_norm = qpn_suffix.lower()
    for qpn_path, qpn_record in qpn_records.items():
        if qpn_norm_path(qpn_path).endswith(qpn_suffix_norm):
            return qpn_record
    return None


def qpn_guard(qpn_errors, qpn_name, qpn_condition, qpn_detail):
    qpn_ok = bool(qpn_condition)
    if not qpn_ok:
        qpn_errors.append({"guard": qpn_name, "detail": qpn_detail})
    return qpn_ok


def qpn_check_source_record(qpn_record, qpn_expected_sha, qpn_required_tokens, qpn_required_names):
    qpn_result = {
        "record_present": qpn_is_mapping(qpn_record),
        "source_text_present": False,
        "provided_sha256": None,
        "computed_sha256": None,
        "sha256_match": False,
        "ast_parse_ok": False,
        "required_top_level_names_present": False,
        "required_static_tokens_present": False,
    }
    if not qpn_is_mapping(qpn_record):
        return qpn_result
    qpn_source_text = qpn_record.get("source_text")
    qpn_provided_sha = qpn_record.get("sha256")
    qpn_result["provided_sha256"] = (
        qpn_provided_sha.lower() if isinstance(qpn_provided_sha, str) else None
    )
    if not isinstance(qpn_source_text, str):
        return qpn_result
    qpn_result["source_text_present"] = True
    qpn_result["computed_sha256"] = qpn_hashlib.sha256(
        qpn_source_text.encode("utf-8")
    ).hexdigest()
    qpn_result["sha256_match"] = (
        qpn_result["computed_sha256"] == qpn_expected_sha
        and qpn_result["provided_sha256"] == qpn_expected_sha
    )
    try:
        qpn_tree = qpn_ast.parse(qpn_source_text, filename="<creator_source_record>")
    except SyntaxError:
        return qpn_result
    qpn_result["ast_parse_ok"] = True
    qpn_top_level_names = {
        qpn_node.name
        for qpn_node in qpn_tree.body
        if isinstance(qpn_node, (qpn_ast.FunctionDef, qpn_ast.AsyncFunctionDef))
    }
    qpn_result["required_top_level_names_present"] = set(qpn_required_names).issubset(
        qpn_top_level_names
    )
    qpn_result["required_static_tokens_present"] = all(
        qpn_token in qpn_source_text for qpn_token in qpn_required_tokens
    )
    return qpn_result


def qpn_row_body(qpn_rows, qpn_node_id):
    if not qpn_is_mapping(qpn_rows):
        return None
    qpn_entry = qpn_rows.get(str(qpn_node_id))
    if not qpn_is_mapping(qpn_entry):
        return None
    qpn_body = qpn_entry.get("body")
    return qpn_body if qpn_is_mapping(qpn_body) else None


def qpn_body_field(qpn_rows, qpn_node_id, qpn_field, qpn_default=None):
    qpn_body = qpn_row_body(qpn_rows, qpn_node_id)
    if qpn_body is None:
        return qpn_default
    return qpn_body.get(qpn_field, qpn_default)


def qpn_unique_ids(qpn_groups):
    qpn_seen = set()
    qpn_out = []
    for qpn_group in qpn_groups:
        for qpn_node_id in qpn_group:
            if qpn_node_id not in qpn_seen:
                qpn_seen.add(qpn_node_id)
                qpn_out.append(qpn_node_id)
    return qpn_out


def qpn_compute_native_scalar_delta(qpn_ey_ratio, qpn_vy_over_e, qpn_positive_e, qpn_source_guards_ok):
    """Return source-gated symbolic scalar relations without running native proof code."""
    if not qpn_source_guards_ok:
        return {
            "status": "blocked",
            "reason": "source_guards_not_complete",
        }
    if not qpn_positive_e:
        return {
            "status": "blocked",
            "reason": "strictly_positive_actual_prescribed_swirl_E_guard_missing",
        }
    if qpn_ey_ratio != "-1/2-mu" or qpn_vy_over_e != "0":
        return {
            "status": "blocked",
            "reason": "native_Ey_over_E_or_Vy_over_E_binding_mismatch",
        }
    return {
        "status": "conditional_source_relation",
        "E_y_over_E": "-1/2-mu",
        "V_y_over_E": "0",
        "a1_minus_2_Ey_over_E": "1-2*(E_y/E)",
        "a_source_substitution": "1-2*(-1/2-mu)",
        "a_literal": "2+2*mu",
        "b_literal": "0",
        "kappa_literal": "2+2*mu",
        "Delta_literal": "2*mu",
        "proof_executed": False,
    }


def qpn_conditional_rho_zero(qpn_all_typed_bound, qpn_delta_result, qpn_flat_gate):
    if not qpn_all_typed_bound or qpn_delta_result.get("status") != "conditional_source_relation":
        return {
            "status": "conditional_under_native_scalar_branch",
            "rho_zero": None,
            "zero_row_count": 20,
            "conditions_met": False,
            "proof_executed": False,
        }
    if not qpn_flat_gate:
        return {
            "status": "conditional_under_native_scalar_branch",
            "rho_zero": None,
            "zero_row_count": 20,
            "conditions_met": False,
            "proof_executed": False,
        }
    return {
        "status": "conditional_under_native_scalar_branch",
        "rho_zero": {"m": 0, "h": 0, "k": 0, "e": 0, "p": 0},
        "zero_row_count": 20,
        "conditions_met": True,
        "proof_executed": False,
    }


def bind_quiet_native_source(spec, journal=None):
    qpn_errors = []
    qpn_spec = spec if qpn_is_mapping(spec) else {}
    qpn_native_envelope = qpn_receipt_envelope(qpn_spec.get("native_receipt"))
    qpn_addendum_envelope = qpn_receipt_envelope(qpn_spec.get("native_literal_addendum"))
    qpn_parent_envelope = qpn_receipt_envelope(qpn_spec.get("parent_result"))
    qpn_native = qpn_native_envelope["payload"]
    qpn_addendum = qpn_addendum_envelope["payload"]
    qpn_parent = qpn_parent_envelope["payload"]
    qpn_records = qpn_spec.get("creator_source_records")
    qpn_cached_admission = qpn_spec.get("cached_admission_receipts")
    qpn_cached_rh = qpn_spec.get("cached_Rh_seed_result")

    qpn_native_shape_ok = qpn_guard(
        qpn_errors,
        "native_receipt_declared_sha_and_payload",
        qpn_native_envelope["has_declared_sha256"] and qpn_native_envelope["payload_mapping"],
        "native_receipt must be {'sha256' or 'raw_sha256': raw_sha, 'payload': actual JSON}",
    )
    qpn_addendum_shape_ok = qpn_guard(
        qpn_errors,
        "native_literal_addendum_declared_sha_and_payload",
        qpn_addendum_envelope["has_declared_sha256"] and qpn_addendum_envelope["payload_mapping"],
        "native_literal_addendum must be {'sha256' or 'raw_sha256': raw_sha, 'payload': actual JSON}",
    )
    qpn_parent_shape_ok = qpn_guard(
        qpn_errors,
        "parent_result_declared_sha_and_payload",
        qpn_parent_envelope["has_declared_sha256"] and qpn_parent_envelope["payload_mapping"],
        "parent_result must be {'sha256' or 'raw_sha256': raw_sha, 'payload': actual JSON}",
    )

    qpn_native_sha_ok = qpn_guard(
        qpn_errors,
        "native_receipt_packet_sha256",
        qpn_native_envelope["declared_sha256"] == qpn_NATIVE_RECEIPT_SHA256,
        "expected native receipt packet SHA 256",
    )
    qpn_addendum_sha_ok = qpn_guard(
        qpn_errors,
        "native_literal_addendum_sha256",
        qpn_addendum_envelope["declared_sha256"] == qpn_LITERAL_ADDENDUM_SHA256,
        "expected literal addendum SHA 256",
    )
    qpn_parent_sha_ok = qpn_guard(
        qpn_errors,
        "parent_result_sha256",
        qpn_parent_envelope["declared_sha256"] == qpn_PARENT_RESULT_SHA256,
        "expected current parent result SHA 256",
    )
    qpn_rh_seed_sha = qpn_sha_from_value(qpn_cached_rh)
    qpn_rh_seed_optional_ok = qpn_cached_rh is None or (
        qpn_rh_seed_sha is not None and qpn_rh_seed_sha.startswith("0c8b3")
    )

    qpn_archive_auth = qpn_addendum.get("archive_authentication", {}) if qpn_is_mapping(qpn_addendum) else {}
    qpn_archive_ok = qpn_guard(
        qpn_errors,
        "native_outer_archive_sha256",
        isinstance(qpn_archive_auth.get("sha256"), str) and qpn_archive_auth.get("sha256").lower() == qpn_NATIVE_ARCHIVE_SHA256,
        "expected true_outer_C6_C5 raw archive SHA 256",
    )
    qpn_native_archives = qpn_native.get("native_archives") if qpn_is_mapping(qpn_native) else None
    if isinstance(qpn_native_archives, list) and qpn_native_archives:
        qpn_archive_ok = qpn_guard(
            qpn_errors,
            "native_receipt_outer_archive_sha256",
            any(
                qpn_is_mapping(qpn_archive_item)
                and isinstance(qpn_archive_item.get("sha256"), str) and qpn_archive_item.get("sha256").lower() == qpn_NATIVE_ARCHIVE_SHA256
                for qpn_archive_item in qpn_native_archives
            ),
            "native receipt does not carry the authenticated outer archive SHA 256",
        ) and qpn_archive_ok

    qpn_rows = qpn_addendum.get("rows", {}) if qpn_is_mapping(qpn_addendum) else {}
    qpn_rows_shape_ok = qpn_guard(
        qpn_errors,
        "literal_rows_body_shape",
        qpn_is_mapping(qpn_rows)
        and all(
            qpn_is_mapping(qpn_entry) and qpn_is_mapping(qpn_entry.get("body"))
            for qpn_entry in qpn_rows.values()
        ),
        "each literal row must be a mapping containing a body mapping",
    )

    qpn_mu = qpn_addendum.get("mu_proof_selected_literal", {}) if qpn_is_mapping(qpn_addendum) else {}
    qpn_family_ok = qpn_guard(
        qpn_errors,
        "same_source_family",
        qpn_is_mapping(qpn_mu.get("source_family"))
        and dict(qpn_mu.get("source_family")) == qpn_FAMILY,
        "cached mu source family mismatch",
    )
    qpn_mu_positive_ok = qpn_guard(
        qpn_errors,
        "mu_positive_and_2mu_ge_eta",
        qpn_mu.get("positive_mu_and_power_2mu_ge_eta_and_Rp_reservation") is True,
        "cached positive-mu and 2mu>=eta admission flag missing",
    )
    qpn_admission = qpn_mu.get("original_admission_binding", {})
    qpn_admission_ok = (
        qpn_is_mapping(qpn_admission)
        and qpn_admission.get("module") == qpn_ADMISSION_MODULE
        and qpn_admission.get("module_sha256") == qpn_ADMISSION_MODULE_SHA256
        and qpn_admission.get("AST_sha256") == qpn_ADMISSION_AST_SHA256
    )
    qpn_admission_ok = qpn_guard(
        qpn_errors,
        "source_mu_admission_binding",
        qpn_admission_ok,
        "cached source_mu_admission module/hash/AST binding mismatch",
    )
    qpn_bounds_ok = qpn_guard(
        qpn_errors,
        "mu_bounds_not_field_values",
        qpn_mu.get("bounds_not_field_values") is True,
        "mu receipt must be tagged as bounds rather than field values",
    )

    qpn_outer_record = qpn_record_for_suffix(
        qpn_records,
        "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py",
    )
    qpn_loop_record = qpn_record_for_suffix(
        qpn_records,
        "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py",
    )
    qpn_outer_source = qpn_check_source_record(
        qpn_outer_record,
        qpn_OUTER_SOURCE_SHA256,
        qpn_SOURCE_TOKEN_GROUPS["outer_geometry"] + qpn_SOURCE_TOKEN_GROUPS["outer_power_context"],
        ("geometry", "transport_geometry"),
    )
    qpn_loop_source = qpn_check_source_record(
        qpn_loop_record,
        qpn_LOOP_SOURCE_SHA256,
        qpn_SOURCE_TOKEN_GROUPS["prescribed_loop"],
        ("loop_functions", "generic_roots", "finite_N_functions", "transport"),
    )
    qpn_sources_ok = qpn_guard(
        qpn_errors,
        "creator_source_records",
        (
            qpn_outer_source["record_present"]
            and qpn_outer_source["source_text_present"]
            and qpn_outer_source["sha256_match"]
            and qpn_outer_source["ast_parse_ok"]
            and qpn_outer_source["required_top_level_names_present"]
            and qpn_outer_source["required_static_tokens_present"]
            and qpn_loop_source["record_present"]
            and qpn_loop_source["source_text_present"]
            and qpn_loop_source["sha256_match"]
            and qpn_loop_source["ast_parse_ok"]
            and qpn_loop_source["required_top_level_names_present"]
            and qpn_loop_source["required_static_tokens_present"]
        ),
        "source records must hash and parse statically with required literal anchors",
    )

    qpn_admission_source_record = qpn_record_for_suffix(qpn_records, qpn_ADMISSION_MODULE)
    qpn_admission_source_checked = False
    qpn_admission_source_info = {"available": False, "checked": False}
    if qpn_admission_source_record is not None:
        qpn_admission_source_info = qpn_check_source_record(
            qpn_admission_source_record,
            qpn_ADMISSION_MODULE_SHA256,
            ("source_mu_admission",),
            ("source_mu_admission",),
        )
        qpn_admission_source_checked = (
            qpn_admission_source_info["sha256_match"]
            and qpn_admission_source_info["ast_parse_ok"]
            and qpn_admission_source_info["required_top_level_names_present"]
            and qpn_admission_source_info["required_static_tokens_present"]
        )
    elif qpn_is_mapping(qpn_cached_admission):
        qpn_admission_candidate = qpn_cached_admission.get("source_mu_admission")
        if qpn_is_mapping(qpn_admission_candidate) and isinstance(
            qpn_admission_candidate.get("source_text"), str
        ):
            qpn_admission_source_info = qpn_check_source_record(
                qpn_admission_candidate,
                qpn_ADMISSION_MODULE_SHA256,
                ("source_mu_admission",),
                ("source_mu_admission",),
            )
            qpn_admission_source_checked = (
                qpn_admission_source_info["sha256_match"]
                and qpn_admission_source_info["ast_parse_ok"]
                and qpn_admission_source_info["required_top_level_names_present"]
                and qpn_admission_source_info["required_static_tokens_present"]
            )
    qpn_admission_source_info["checked"] = qpn_admission_source_checked
    qpn_admission_receipt_as_data = not qpn_admission_source_checked

    qpn_roles = qpn_addendum.get("row_roles", {}) if qpn_is_mapping(qpn_addendum) else {}
    qpn_root_roles = qpn_roles.get("snapshot_root_maps", {}) if qpn_is_mapping(qpn_roles) else {}
    qpn_background_roles = qpn_roles.get("snapshot_background_fields", {}) if qpn_is_mapping(qpn_roles) else {}
    qpn_radial_roles = qpn_roles.get("snapshot_radial_fields", {}) if qpn_is_mapping(qpn_roles) else {}
    qpn_memory_roles = qpn_roles.get("memory_rows", {}) if qpn_is_mapping(qpn_roles) else {}
    qpn_native_rows_ids = qpn_unique_ids(
        [
            [0, 1, 2, 3, 18, 21, 48, 205156, 250215, 250222, 250296, 250297, 250324, 250325, 250226, 251146, 251216, 251222],
            qpn_root_roles.get("E", [])[:4],
            qpn_root_roles.get("a", [])[:1],
            qpn_background_roles.get("E", [])[:4],
            qpn_radial_roles.get("E_y", [])[:4],
            qpn_radial_roles.get("V_y", [])[:4],
            qpn_FLAT_A_ROWS,
            qpn_FLAT_B_ROWS,
            qpn_A_OVER_N_FIRST4,
            qpn_DELTA_E_ROWS,
            qpn_DELTA_V_FIRST4,
            qpn_unique_ids(list(qpn_FINITE_HISTORY.values())),
            qpn_unique_ids(list(qpn_LOCAL_HISTORY.values())),
            qpn_ENDPOINT_ROWS,
        ]
    )
    qpn_required_rows_ok = qpn_guard(
        qpn_errors,
        "required_native_literal_rows",
        qpn_rows_shape_ok and all(qpn_row_body(qpn_rows, qpn_id) is not None for qpn_id in qpn_native_rows_ids),
        "required O3_power source rows are missing",
    )

    qpn_root_map_ok = (
        qpn_root_roles.get("E", [])[:4] == qpn_E_ROOT_FIRST4
        and qpn_root_roles.get("a", [])[:4] == qpn_A_ROOT_FIRST4
        and qpn_root_roles.get("V", [])[:4] == [0, 0, 0, 0]
        and qpn_root_roles.get("b", [])[:4] == [0, 0, 0, 0]
        and qpn_root_roles.get("Delta", [])[:4] == [251216, 251150, 251156, 251163]
    )
    qpn_root_map_ok = qpn_guard(
        qpn_errors,
        "snapshot_root_maps",
        qpn_root_map_ok,
        "O3_power root E/V/a/b/Delta first-four maps mismatch",
    )
    qpn_a0_body = qpn_row_body(qpn_rows, 251146) or {}
    qpn_e_positive_ok = qpn_guard(
        qpn_errors,
        "actual_E_positive_quotient",
        qpn_a0_body.get("operation") == "positive_quotient"
        and qpn_a0_body.get("denominator") == 250226
        and qpn_a0_body.get("source_positive_certificate")
        == "same strictly positive actual prescribed swirl E",
        "E positivity must use node 251146 quotient certificate over root E node 250226",
    )
    qpn_n_body = qpn_row_body(qpn_rows, 19) or {}
    qpn_n_positive_ok = qpn_guard(
        qpn_errors,
        "N_positive_descriptor",
        qpn_n_body.get("operation") == "positive_quotient"
        and qpn_n_body.get("denominator") == 18,
        "node 19 is only the N positivity certificate",
    )
    qpn_n_descriptor = qpn_row_body(qpn_rows, 18) or {}
    qpn_n_descriptor_ok = qpn_guard(
        qpn_errors,
        "N_node_18_descriptor",
        qpn_n_descriptor.get("operation") == "exact_positive_integer_expression"
        and qpn_n_descriptor.get("J") == "1358356628656378313"
        and qpn_n_descriptor.get("exact_positive_integer") is True,
        "node 18 descriptor/J binding mismatch",
    )
    qpn_literal_node_48 = qpn_row_body(qpn_rows, 48) or {}
    qpn_two_literal_ok = qpn_guard(
        qpn_errors,
        "literal_two_node_48",
        qpn_literal_node_48.get("operation") == "exact_rational"
        and qpn_literal_node_48.get("numerator") == 2
        and qpn_literal_node_48.get("denominator") == 1,
        "node 48 must be exact rational 2",
    )
    qpn_offset_body = qpn_row_body(qpn_rows, 250222) or {}
    qpn_offset_18_body = qpn_row_body(qpn_rows, 205156) or {}
    qpn_offset_exp_body = qpn_row_body(qpn_rows, 3) or {}
    qpn_offset_40_body = qpn_row_body(qpn_rows, 2) or {}
    qpn_coordinate_body = qpn_row_body(qpn_rows, 250215) or {}
    qpn_zero_body = qpn_row_body(qpn_rows, 0) or {}
    qpn_one_body = qpn_row_body(qpn_rows, 1) or {}
    qpn_geometry_rows_ok = qpn_guard(
        qpn_errors,
        "geometry_literal_rows",
        qpn_coordinate_body.get("operation") == "bound_variable"
        and qpn_coordinate_body.get("name") == "true_outer_O3_power"
        and qpn_offset_body.get("operation") == "sum"
        and qpn_offset_body.get("arguments") == [205156, 3, 250215]
        and qpn_offset_18_body.get("operation") == "exact_rational"
        and qpn_offset_18_body.get("numerator") == 18
        and qpn_offset_18_body.get("denominator") == 1
        and qpn_offset_exp_body.get("operation") == "analytic_unary"
        and qpn_offset_exp_body.get("name") == "exp"
        and qpn_offset_exp_body.get("argument") == 2
        and qpn_offset_40_body.get("operation") == "exact_rational"
        and qpn_offset_40_body.get("numerator") == 40
        and qpn_offset_40_body.get("denominator") == 1
        and qpn_zero_body.get("numerator") == 0
        and qpn_zero_body.get("denominator") == 1
        and qpn_one_body.get("numerator") == 1
        and qpn_one_body.get("denominator") == 1,
        "native coordinate/offset/domain/Jacobian literal row chain mismatch",
    )
    qpn_mu_recipe_ok = qpn_guard(
        qpn_errors,
        "mu_source_exp_recipe",
        qpn_mu.get("exact_original_mu_recipe") == "mu=1/1000*exp(-4*(exp(40)+11))",
        "exact source mu recipe mismatch",
    )
    qpn_mu_node_17 = qpn_row_body(qpn_rows, 17) or {}
    qpn_mu_node_ok = qpn_guard(
        qpn_errors,
        "mu_source_exp_node_17",
        qpn_mu_node_17.get("operation") == "analytic_unary"
        and qpn_mu_node_17.get("name") == "exp",
        "node 17 exp source tag missing",
    )

    qpn_native_o3 = qpn_native.get("o3_power", {}) if qpn_is_mapping(qpn_native) else {}
    qpn_geometry_expectations = {
        "coordinate_node": 250215,
        "domain_archive_nodes": [0, 48],
        "domain_semantic": [0, 2],
        "physical_measure_node": 1,
        "log_offset_node": 250222,
        "log_offset_literal": "18 + exp(40) + x",
        "jacobian_literal": "1",
    }
    qpn_geometry_ok = True
    for qpn_key, qpn_expected in qpn_geometry_expectations.items():
        qpn_geometry_ok = qpn_guard(
            qpn_errors,
            "native_geometry_" + qpn_key,
            qpn_native_o3.get(qpn_key) == qpn_expected,
            "native O3_power geometry literal mismatch",
        ) and qpn_geometry_ok
    qpn_memory_ok = qpn_guard(
        qpn_errors,
        "memory_rows",
        qpn_memory_roles.get("m") == 250297
        and qpn_memory_roles.get("h") == 250325
        and qpn_memory_roles.get("k") == 250325
        and qpn_memory_roles.get("e") == 250297
        and qpn_memory_roles.get("p") == 1,
        "O3_power memory row IDs mismatch",
    )
    qpn_memory_shape_ok = qpn_guard(
        qpn_errors,
        "memory_exp_negative_rate_width_shape",
        qpn_row_body(qpn_rows, 250297).get("operation") == "analytic_unary"
        and qpn_row_body(qpn_rows, 250297).get("name") == "exp"
        and qpn_row_body(qpn_rows, 250297).get("argument") == 250296
        and qpn_row_body(qpn_rows, 250296) is not None
        and qpn_row_body(qpn_rows, 250325).get("operation") == "analytic_unary"
        and qpn_row_body(qpn_rows, 250325).get("name") == "exp"
        and qpn_row_body(qpn_rows, 250325).get("argument") == 250324
        and qpn_row_body(qpn_rows, 250324) is not None
        and qpn_row_body(qpn_rows, 1).get("operation") == "exact_rational",
        "memory rows must be native exp rows for m/e and h/k, with p=1",
    )

    qpn_predicate_body = qpn_row_body(qpn_rows, 251222) or {}
    qpn_predicate_ok = qpn_guard(
        qpn_errors,
        "flat_predicate",
        qpn_predicate_body.get("operation") == "exact_real_comparison"
        and qpn_predicate_body.get("operator") == ">="
        and qpn_predicate_body.get("left") == 251216
        and qpn_predicate_body.get("right") == 21,
        "flat predicate must be Delta >= eta with left 251216 and right 21",
    )
    qpn_lazy_row_ok = True
    for qpn_row_group in (qpn_FLAT_A_ROWS, qpn_FLAT_B_ROWS):
        for qpn_index, qpn_node_id in enumerate(qpn_row_group):
            qpn_body = qpn_row_body(qpn_rows, qpn_node_id) or {}
            qpn_this_ok = (
                qpn_body.get("operation") == "original_lazy_flat_branch"
                and qpn_body.get("Delta") == 251216
                and qpn_body.get("eta") == 21
                and qpn_body.get("flat_value") == 0
                and qpn_body.get("flat_predicate") == 251222
                and qpn_body.get("active_body_not_evaluated_on_flat_branch") is True
                and qpn_body.get("ordinary_slow_Z_derivative_order") == qpn_index
                and qpn_body.get("same_original_flat_extension_through_order5") is True
                and qpn_body.get("original_collar_coordinate") is None
                and qpn_body.get("original_collar_upper") is None
            )
            qpn_lazy_row_ok = qpn_guard(
                qpn_errors,
                "lazy_flat_row_" + str(qpn_node_id),
                qpn_this_ok,
                "A/B lazy flat row literal attributes mismatch",
            ) and qpn_lazy_row_ok
    qpn_lazy_row_ok = qpn_guard(
        qpn_errors,
        "lazy_flat_rows",
        qpn_lazy_row_ok,
        "one or more A/B lazy flat rows failed",
    )

    qpn_history_ok = True
    qpn_local_integral_ok = True
    qpn_incoming = qpn_native_o3.get("correction_incoming_nodes")
    for qpn_species in ("m", "h", "k", "e", "p"):
        qpn_species_incoming = (
            qpn_incoming.get(qpn_species, [])[:4]
            if qpn_is_mapping(qpn_incoming)
            else []
        )
        qpn_history_ids = qpn_FINITE_HISTORY[qpn_species]
        qpn_local_ids = qpn_LOCAL_HISTORY[qpn_species]
        qpn_memory_id = qpn_MEMORY[qpn_species]
        qpn_incoming_ok = qpn_guard(
            qpn_errors,
            "incoming_first4_" + qpn_species,
            qpn_species_incoming == qpn_INCOMING_FIRST4[qpn_species],
            "raw incoming node IDs mismatch",
        )
        qpn_history_ok = qpn_incoming_ok and qpn_history_ok
        for qpn_index, qpn_history_id in enumerate(qpn_history_ids):
            qpn_history_body = qpn_row_body(qpn_rows, qpn_history_id) or {}
            qpn_history_args = qpn_history_body.get("arguments")
            qpn_product_id = (
                qpn_history_args[0]
                if isinstance(qpn_history_args, list) and len(qpn_history_args) == 2
                else None
            )
            qpn_local_id = (
                qpn_history_args[1]
                if isinstance(qpn_history_args, list) and len(qpn_history_args) == 2
                else None
            )
            qpn_product_body = qpn_row_body(qpn_rows, qpn_product_id) or {}
            qpn_product_args = qpn_product_body.get("arguments")
            qpn_history_this_ok = (
                qpn_history_body.get("operation") == "sum"
                and qpn_local_id == qpn_local_ids[qpn_index]
                and qpn_product_body.get("operation") == "product"
                and qpn_product_args == [qpn_INCOMING_FIRST4[qpn_species][qpn_index], qpn_memory_id]
            )
            qpn_history_ok = qpn_guard(
                qpn_errors,
                "finite_history_" + qpn_species + "_" + str(qpn_index),
                qpn_history_this_ok,
                "finite history must be incoming*memory + local",
            ) and qpn_history_ok
            qpn_local_body = qpn_row_body(qpn_rows, qpn_local_ids[qpn_index]) or {}
            qpn_local_this_ok = (
                qpn_local_body.get("operation") == "definite_integral"
                and qpn_local_body.get("exact_function_integral") is True
                and qpn_local_body.get("endpoint_independent_of_Z") is True
                and qpn_local_body.get("physical_Jacobian_applied_once") is True
                and qpn_local_body.get("physical_Jacobian") == 1
                and qpn_local_body.get("physical_log_offset") == 250222
                and qpn_local_body.get("lower") == 0
                and qpn_local_body.get("upper") == 250215
                and qpn_local_body.get("original_chart") == "O3_power"
            )
            qpn_local_integral_ok = qpn_guard(
                qpn_errors,
                "local_integral_" + qpn_species + "_" + str(qpn_index),
                qpn_local_this_ok,
                "local integral must be exact, Z-independent, and physically Jacobian-once",
            ) and qpn_local_integral_ok
    qpn_history_ok = qpn_guard(qpn_errors, "finite_history_20_rows", qpn_history_ok, "finite history row binding failed")
    qpn_local_integral_ok = qpn_guard(qpn_errors, "local_integral_20_rows", qpn_local_integral_ok, "local integral row binding failed")

    qpn_endpoint_ok = True
    for qpn_endpoint_id in qpn_ENDPOINT_ROWS:
        qpn_endpoint_body = qpn_row_body(qpn_rows, qpn_endpoint_id) or {}
        qpn_endpoint_ok = qpn_guard(
            qpn_errors,
            "endpoint_substitution_" + str(qpn_endpoint_id),
            qpn_endpoint_body.get("operation") == "function_substitution"
            and qpn_endpoint_body.get("variable") == 250215
            and qpn_endpoint_body.get("value") == 48
            and qpn_endpoint_body.get("Z_independent_substitution") is True,
            "endpoint substitution must be x=2 via variable 250215/value 48",
        ) and qpn_endpoint_ok
    qpn_endpoint_ok = qpn_guard(qpn_errors, "endpoint_substitutions_20_rows", qpn_endpoint_ok, "endpoint row binding failed")

    qpn_parent_roles = {}
    if qpn_is_mapping(qpn_parent):
        for qpn_key in ("N_node", "N_positive", "normalization", "N_scaled_incoming", "n_scaled_incoming"):
            if qpn_key in qpn_parent:
                qpn_parent_roles[qpn_key] = qpn_parent.get(qpn_key)
    qpn_parent_roles["reuse_as_data_only"] = True
    qpn_parent_roles["new_old_scaling_proof"] = False

    qpn_guard_values = {
        "native_receipt_declared_shape": qpn_native_shape_ok,
        "native_literal_addendum_declared_shape": qpn_addendum_shape_ok,
        "parent_result_declared_shape": qpn_parent_shape_ok,
        "native_receipt_packet_sha256": qpn_native_sha_ok,
        "native_literal_addendum_sha256": qpn_addendum_sha_ok,
        "parent_result_sha256": qpn_parent_sha_ok,
        "cached_Rh_seed_optional": qpn_rh_seed_optional_ok,
        "native_outer_archive_sha256": qpn_archive_ok,
        "literal_rows_body_shape": qpn_rows_shape_ok,
        "same_source_family": qpn_family_ok,
        "mu_positive_and_2mu_ge_eta": qpn_mu_positive_ok,
        "source_mu_admission_binding": qpn_admission_ok,
        "mu_bounds_not_field_values": qpn_bounds_ok,
        "creator_source_records": qpn_sources_ok,
        "required_native_literal_rows": qpn_required_rows_ok,
        "snapshot_root_maps": qpn_root_map_ok,
        "actual_E_positive_quotient": qpn_e_positive_ok,
        "N_positive_descriptor": qpn_n_positive_ok,
        "N_node_18_descriptor": qpn_n_descriptor_ok,
        "literal_two_node_48": qpn_two_literal_ok,
        "mu_source_exp_recipe": qpn_mu_recipe_ok,
        "mu_source_exp_node_17": qpn_mu_node_ok,
        "native_geometry": qpn_geometry_ok,
        "geometry_literal_rows": qpn_geometry_rows_ok,
        "memory_rows": qpn_memory_ok,
        "memory_exp_negative_rate_width_shape": qpn_memory_shape_ok,
        "flat_predicate": qpn_predicate_ok,
        "lazy_flat_rows": qpn_lazy_row_ok,
        "finite_history_20_rows": qpn_history_ok,
        "local_integral_20_rows": qpn_local_integral_ok,
        "endpoint_substitutions_20_rows": qpn_endpoint_ok,
    }
    qpn_all_typed_bound = all(qpn_guard_values.values()) and qpn_rh_seed_optional_ok

    qpn_delta_result = qpn_compute_native_scalar_delta(
        "-1/2-mu",
        "0",
        qpn_e_positive_ok,
        qpn_all_typed_bound,
    )
    qpn_rho_result = qpn_conditional_rho_zero(
        qpn_all_typed_bound,
        qpn_delta_result,
        qpn_predicate_ok and qpn_lazy_row_ok,
    )

    qpn_needed_bodies = {}
    for qpn_node_id in qpn_native_rows_ids:
        qpn_body = qpn_row_body(qpn_rows, qpn_node_id)
        if qpn_body is not None:
            qpn_needed_bodies[str(qpn_node_id)] = qpn_body

    qpn_source_gate = {
        "name": "native_flat_scalar_source_gate",
        "all_typed_bound": qpn_all_typed_bound,
        "guards": qpn_guard_values,
        "source_family": dict(qpn_FAMILY),
        "archive_sha256": qpn_NATIVE_ARCHIVE_SHA256,
        "geometry": {
            "coordinate_node": 250215,
            "offset_node": 250222,
            "offset_literal": "18+exp(40)+x",
            "domain_archive_nodes": [0, 48],
            "domain_semantic": [0, 2],
            "jacobian_literal": "1",
            "physical_measure_literal": "1",
        },
        "rates": dict(qpn_RATES),
        "predicate": {
            "node": 251222,
            "operator": ">=",
            "left": 251216,
            "right": 21,
            "Delta_literal": "2*mu",
        },
        "ey_coefficients": [[-1, 2], [-1, 1]],
        "vy_over_e_coefficients": [[0, 1], [0, 1]],
        "zero_row_count": 20,
        "conditional_under_native_scalar_branch": True,
        "incoming_outgoing_local_scope": "incoming/outgoing/local O3_power rows only; not m15, P1, or full physical whole-17 scope",
        "unknown_entry_difference": "retained; neither zero nor nonzero asserted",
    }

    qpn_journal_event = {
        "event": "bind_quiet_native_source",
        "all_typed_bound": qpn_all_typed_bound,
        "project_source_executed": False,
        "proof_checker_normalizer_math_helper_calls": 0,
        "errors": len(qpn_errors),
    }
    if isinstance(journal, list):
        journal.append(qpn_journal_event)
        qpn_journal_out = list(journal)
    else:
        qpn_journal_out = [qpn_journal_event]

    return {
        "status": "bound" if qpn_all_typed_bound else "blocked",
        "spec_digest_field_keys": {
            "native_receipt_wrapper": ["sha256", "raw_sha256", "_raw_sha256", "payload"],
            "native_literal_addendum_wrapper": ["sha256", "raw_sha256", "_raw_sha256", "payload"],
            "parent_result_wrapper": ["sha256", "raw_sha256", "_raw_sha256", "payload"],
            "payload_metadata_key": "_raw_sha256",
        },
        "receipt_wrapper_policy": "declared raw SHA is mandatory; payload is the supplied actual JSON and is not mutated",
        "source_producer_refs": {
            "outer_source_sha256": qpn_OUTER_SOURCE_SHA256,
            "prescribed_loop_source_sha256": qpn_LOOP_SOURCE_SHA256,
            "native_power_E_V_Ey_Vy": "outer source lines 180-183",
            "generic_roots": "prescribed loop source lines 213-223",
            "Delta_binding_target": "Delta=2*mu from E_y/E=-1/2-mu and V_y=0",
        },
        "source_gate": qpn_source_gate,
        "pure_scalar_delta": qpn_delta_result,
        "rho_zero_conditional": qpn_rho_result,
        "native_rows": {
            "row_ids": qpn_native_rows_ids,
            "bodies": qpn_needed_bodies,
            "incoming_first4": dict(qpn_INCOMING_FIRST4),
            "memory_rows": dict(qpn_MEMORY),
            "N_node": 18,
            "N_multiplier_node": 18,
            "N_descriptor_J_decimal_string": qpn_n_descriptor.get("J"),
        },
        "admission": {
            "cached_receipt_as_data": qpn_admission_receipt_as_data,
            "source_text_checked_if_available": qpn_admission_source_checked,
            "binding": dict(qpn_admission) if qpn_is_mapping(qpn_admission) else {},
            "mu_recipe": qpn_mu.get("exact_original_mu_recipe"),
        },
        "E_nonzero_data": {
            "positivity_certificate_node": 251146,
            "positivity_certificate_denominator": 250226,
            "positivity_certificate": qpn_a0_body.get("source_positive_certificate"),
            "N_positivity_certificate_node": 19,
            "accepted_Rh_seed_sha256": qpn_rh_seed_sha,
            "new_numeric_value": False,
        },
        "finite_N": {
            "N_descriptor_node": 18,
            "N_positive_inherited_parent": qpn_parent_sha_ok,
            "invN_literal": "g.quotient(g.one,N,'same exact positive integer N on every source chart')",
            "A_over_N_first4": list(qpn_A_OVER_N_FIRST4),
            "delta_E_rows": list(qpn_DELTA_E_ROWS),
            "delta_V_first4": list(qpn_DELTA_V_FIRST4),
            "r_literal": "a.mul(roots['E'],exm1)",
            "s_literal": "a.scale(loop['B_over_Pstar'],invN)",
            "signed_density_increment_literals": {
                "m": "s",
                "h": "r",
                "k": "a.add(a.mul(V,r),a.mul(E,s),a.mul(r,s))",
                "e": "a.sub(a.add(a.scale(a.mul(V,s),2),a.mul(s,s)),a.add(a.mul(E,r),a.scale(a.mul(r,r),'1/2')))",
                "p": "a.add(a.mul(E,r),a.scale(a.mul(r,r),'1/2'))",
            },
            "formula_source_only": True,
            "parent_normalization_reused_as_data": True,
            "new_old_scaling_proof": False,
        },
        "parent_result_data": qpn_parent_roles,
        "errors": qpn_errors,
        "journal": qpn_journal_out,
    }


"""Unexecuted source-bound glue for the compact quiet-power history pair.

This wrapper is concatenated after the exact algebra draft and the two static
source binders.  It consumes only the frozen input packet, calls the C3 binder
first, then the native binder, and calls the new algebra helper once only when
both source gates are bound.  It never imports a producer, replays a generic
proof, evaluates a source value, computes N, evaluates an exponential, or
asserts an incoming-gap equality.  The endpoint transport is metadata only:
``deltaT_out=exp(-2*lambda)*deltaT_in`` at coordinate 2, with predecessor
history retained.
"""

import hashlib as qpw_hashlib
import json as qpw_json


QPW_INPUT_SCHEMA = "compact_quiet_power_history_pair_input_v1"
QPW_RESULT_SCHEMA = "original_compact_quiet_power_history_pair"
QPW_COMPONENT = "compact_quiet_power_history_pair"
QPW_PARENT_RESULT_SHA256 = (
    "b91da155898c38a9ca4b26595983a38aa04b2965ae3db0a908fadcea97b34241"
)
QPW_C3_RECEIPT_SHA256 = (
    "a7ee951608442791e643dc0167d5442c486b6f41cf792bfcc4130550105ccd4d"
)
QPW_NATIVE_RECEIPT_SHA256 = (
    "c963f41b44901683c99dd76e6d26e9d58d360dff54d93941ea5722b773bf8cec"
)
QPW_NATIVE_LITERAL_ADDENDUM_SHA256 = (
    "3978472d568f94fab5ed852580e404e74d32b27a7db874b03f793fe46efe99fc"
)
QPW_RH_SEED_SHA256 = (
    "0c8b3ca3f53f95bef647332eb7c6edb46d29cea230c377ffa0acabdbd07e2253"
)
QPW_C3_ARCHIVE_SHA256 = (
    "cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193"
)
QPW_NATIVE_ARCHIVE_SHA256 = (
    "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d"
)
QPW_PARENT_NATIVE_TERMINAL_ARCHIVE_SHA256 = (
    "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269"
)
QPW_COMPACT_C5_SOURCE_SHA256 = (
    "f0b0333b19a0e2504d023b97e841033d3cbcc23427c57ad9672264866c01ec7b"
)
QPW_COMPACT_C5_ENDPOINT_IDS = {
    "m": [258646, 258647, 258648, 258649],
    "h": [258652, 258653, 258654, 258655],
    "k": [258658, 258659, 258660, 258661],
    "e": [258664, 258665, 258666, 258667],
    "p": [258670, 258671, 258672, 258673],
}
QPW_NATIVE_INCOMING_IDS = {
    "m": [250004, 250005, 250006, 250007],
    "h": [250010, 250011, 250012, 250013],
    "k": [250016, 250017, 250018, 250019],
    "e": [250022, 250023, 250024, 250025],
    "p": [250028, 250029, 250030, 250031],
}
QPW_C3_INCOMING_IDS = {
    "m": [4820, 4821, 12010, 24145],
    "h": [4842, 4843, 12018, 24153],
    "k": [4858, 4859, 12026, 24161],
    "e": [4874, 4875, 12034, 24169],
    "p": [4888, 4889, 12041, 24176],
}
QPW_FINITE_HISTORY_IDS = {
    "m": [253257, 253258, 253259, 253260],
    "h": [253290, 253291, 253292, 253293],
    "k": [253320, 253321, 253322, 253323],
    "e": [253350, 253351, 253352, 253353],
    "p": [253374, 253375, 253376, 253377],
}
QPW_FIELDS = ("m", "h", "k", "e", "p")
QPW_ORDERS = (0, 1, 2, 3)
QPW_NEW_FLAGS = (
    "quiet_power_native_density_zero_source_proved",
    "quiet_power_C3_density_zero_source_proved",
    "quiet_power_common_history_kernel_source_bound",
    "quiet_power_signed_history_gap_transport_proved",
)
QPW_CURRENT_RESULT_FIELDS = frozenset(
    {
        "schema",
        "component",
        "passed",
        "metadata_only",
        "source_only",
        "input_digest",
        "runtime_model_metadata",
        "original_source_evaluations",
        "source_packet_defined",
        "kernel_imported",
        "derive_entry_calls",
        "derive_calls",
        "entry_attempts",
    }
)
QPW_FORBIDDEN_TRUE_FLAGS = (
    "actual_C3_C5_incoming_source_identity_proved",
    "actual_C3_C5_control_source_identity_proved",
    "actual_C3_C5_profile_identity_proved",
    "actual_C3_cap_to_C5_source_bound_proved",
    "compact_volume_cap_same_source_bound_proved",
    "global_five_moment_closure",
    "numeric_values_computed",
    "full_goal_complete",
)
QPW_SOURCE_SCOPE = {
    "field": "quiet_O3_power",
    "ordinary_Z_orders": [0, 1, 2, 3],
    "coordinate_interval": ["0", "2"],
    "ordinary_Z_domain": ["-1", "1"],
    "ordinary_Z_function_identity_not_samples": True,
    "zero_local_densities_only": True,
    "zero_density_whole_interval": True,
    "zero_density_provenance": "source_bound_flat_branch_and_admission_certificate",
    "existing_C3_partial_history_claimed": False,
}
QPW_KERNEL_SCOPE = {
    "same_relative_log_offset": True,
    "same_physical_measure": True,
    "absolute_Rm_Rc_cross_archive_identity_claimed": False,
    "source_integrand_kernel": "exp(-lambda*(2-s))",
    "memory_kernel": "exp(-2*lambda)",
}
QPW_TRANSPORT_SCOPE = {
    "formula": "deltaT_out=exp(-2*lambda)*deltaT_in",
    "endpoint_coordinate": "2",
    "ordinary_Z_orders": [0, 1, 2, 3],
    "unknown_entry_gap_retained": True,
    "predecessor_inlet_retained": True,
    "predecessor_inlet_zeroed": False,
}


def qpw_require(condition, message):
    if not condition:
        raise ValueError(message)


def qpw_clone(value):
    return qpw_json.loads(qpw_json.dumps(value, ensure_ascii=False, allow_nan=False))


def qpw_sha(value):
    if isinstance(value, bytes):
        raw = value
    elif isinstance(value, str):
        raw = value.encode("utf-8")
    else:
        raw = qpw_json.dumps(
            value, ensure_ascii=False, sort_keys=True, allow_nan=False
        ).encode("utf-8")
    return qpw_hashlib.sha256(raw).hexdigest()


def qpw_emit(journal, event, **fields):
    if journal is None:
        return
    payload = {}
    if isinstance(event, dict):
        payload.update(event)
    elif event is not None:
        payload["event"] = event
    payload.update(fields)
    if callable(journal):
        journal(payload)
    elif isinstance(journal, list):
        journal.append(payload)


def qpw_parent(source, spec):
    candidate = spec.get("parent_result") if isinstance(spec, dict) else None
    if not isinstance(candidate, dict):
        candidate = source.get("parent_result") if isinstance(source, dict) else None
    return candidate if isinstance(candidate, dict) else {}


def qpw_inherited_flags(source, parent):
    inherited = source.get("inherited_parent_true_flags", {})
    if not isinstance(inherited, dict):
        inherited = {}
    if not inherited and isinstance(parent, dict):
        false_flags = parent.get("false_flags", [])
        if not isinstance(false_flags, list):
            false_flags = []
        inherited = {
            key: value
            for key, value in parent.items()
            if isinstance(key, str)
            and isinstance(value, bool)
            and value is True
            and key not in false_flags
        }
    return qpw_clone(inherited)


def qpw_apply_parent_metadata(result, source, parent):
    for key, value in qpw_inherited_flags(source, parent).items():
        if value is True and key not in QPW_CURRENT_RESULT_FIELDS:
            result[key] = True
    false_flags = source.get("inherited_false_flags", [])
    if isinstance(false_flags, list):
        for key in false_flags:
            if isinstance(key, str):
                result[key] = False
    for key in QPW_FORBIDDEN_TRUE_FLAGS:
        result[key] = False
    return result


def qpw_parent_data(source, parent):
    raw_sha = parent.get("_raw_sha256") if isinstance(parent, dict) else None
    if not isinstance(raw_sha, str):
        raw_sha = QPW_PARENT_RESULT_SHA256
    return {
        "parent_result_raw_sha256": raw_sha,
        "parent_result_data_only": True,
        "parent_result_retained_exact": qpw_clone(parent),
        "retained_construction_data": qpw_clone(
            parent.get("retained_construction_data")
        ),
        "retained_integral_term_inventory": qpw_clone(
            parent.get("retained_integral_term_inventory")
        ),
        "retained_cap_provenance": qpw_clone(parent.get("retained_cap_provenance")),
        "parent_source_packet_defined": parent.get("source_packet_defined"),
        "parent_metadata_only": parent.get("metadata_only"),
        "parent_passed": parent.get("passed"),
    }


def qpw_scope():
    return {
        "source": qpw_clone(QPW_SOURCE_SCOPE),
        "common_kernel": qpw_clone(QPW_KERNEL_SCOPE),
        "history_transport": qpw_clone(QPW_TRANSPORT_SCOPE),
    }


def qpw_native_scalar_gate(native):
    gate = native.get("source_gate", {}) if isinstance(native, dict) else {}
    guards = gate.get("guards", {}) if isinstance(gate, dict) else {}
    delta = native.get("pure_scalar_delta", {}) if isinstance(native, dict) else {}
    rho = native.get("rho_zero_conditional", {}) if isinstance(native, dict) else {}
    ey = gate.get("ey_coefficients") if isinstance(gate, dict) else None
    vy = gate.get("vy_over_e_coefficients") if isinstance(gate, dict) else None
    typed = gate.get("all_typed_bound") is True
    exact_affine = ey == [[-1, 2], [-1, 1]] and vy == [[0, 1], [0, 1]]
    positive = guards.get("actual_E_positive_quotient") is True
    mu_gate = guards.get("mu_positive_and_2mu_ge_eta") is True
    root_gate = bool(
        native.get("status") == "bound"
        and typed
        and exact_affine
        and positive
        and mu_gate
        and delta.get("status") == "conditional_source_relation"
        and delta.get("proof_executed") is False
        and rho.get("conditions_met") is True
    )
    return {
        "source_assignments_bound": typed and guards.get("snapshot_root_maps") is True,
        "Ey_over_E_coefficients_bound": exact_affine,
        "Ey_over_E_coefficients": qpw_clone(ey),
        "Vy_source_zero_bound": exact_affine and guards.get("flat_predicate") is True,
        "Vy_over_E_coefficients": qpw_clone(vy),
        "E_nonzero_source_bound": positive and guards.get("N_positive_descriptor") is True,
        "E_positive_source_bound": positive,
        "mu_positive_source_bound": mu_gate,
        "mu_Z_independent": typed and guards.get("native_geometry") is True,
        "cached_2mu_ge_eta_bound": mu_gate,
        "actual_lazy_predicate_bound": guards.get("flat_predicate") is True
        and guards.get("lazy_flat_rows") is True,
        "root_gate": root_gate,
        "recomputed": False,
        "source_gate_status": native.get("status"),
        "E_positive_certificate_node": 251146,
        "N_positive_descriptor_node": 19,
    }


def qpw_c3_density_gate(c3):
    proof = c3.get("quiet_proof", {}) if isinstance(c3, dict) else {}
    source_checks = c3.get("source_checks", {}) if isinstance(c3, dict) else {}
    source_checks_passed = bool(
        isinstance(source_checks, dict)
        and all(
            isinstance(value, dict) and value.get("passed") is True
            for value in source_checks.values()
        )
    )
    passed = bool(
        isinstance(c3, dict)
        and c3.get("status") == "passed"
        and c3.get("passed") is True
        and c3.get("zero_row_count") == 20
        and proof.get("passed") is True
        and source_checks_passed
    )
    return {
        "passed": passed,
        "data_gate": passed,
        "native_rows": 20,
        "C3_rows": c3.get("zero_row_count", 0) if isinstance(c3, dict) else 0,
        "C3_fixed_zero_root_verified": passed,
        "source_bound": passed,
        "recomputed": False,
        "zero_scope": "local O3 densities only; whole coordinate interval [0,2]",
    }


def qpw_geometry(native, c3):
    native_gate = native.get("source_gate", {}) if isinstance(native, dict) else {}
    ng = native_gate.get("geometry", {}) if isinstance(native_gate, dict) else {}
    cg = c3.get("geometry", {}) if isinstance(c3, dict) else {}
    native_ok = (
        native_gate.get("all_typed_bound") is True
        and ng.get("coordinate_node") == 250215
        and ng.get("offset_literal") in ("18+exp(40)+x", "18 + exp(40) + x")
        and ng.get("domain_semantic") == [0, 2]
        and ng.get("jacobian_literal") == "1"
        and ng.get("physical_measure_literal") == "1"
    )
    c3_ok = (
        cg.get("chart") == "O3_power"
        and cg.get("domain") == [0, 2]
        and cg.get("offset_source_literal") == "Rm+18+sy.exp(40)+x"
        and cg.get("measure") == "original dlogR: Jacobian exactly once"
        and cg.get("jacobian_applied_once") is True
        and cg.get("absolute_Rm_Rc_crossarchive_alias") is False
    )
    affine = bool(native_ok and c3_ok)
    rates = native_gate.get("rates") if isinstance(native_gate, dict) else None
    expected_rates = {"m": "1", "h": "3/2", "k": "3/2", "e": "1", "p": "0"}
    rates_ok = rates == expected_rates
    function_keys = (
        {"m": "exp(-2)", "h": "exp(-3)", "k": "exp(-3)", "e": "exp(-2)", "p": "1"}
        if rates_ok and affine
        else {}
    )
    common = {
        "actual_O3_power_affine": affine,
        "O3_chart_stage": 3 if affine else None,
        "order_scope": "O3 chart stage; not temporal n=3",
        "log_offset": {
            "expression": "exp(40)+18+x",
            "exp_argument": [40, 1] if affine else None,
            "constant": [18, 1] if affine else None,
            "slope": [1, 1] if affine else None,
            "variable": "x",
        },
        "jacobian": [1, 1] if affine else None,
        "lower": [0, 1] if affine else None,
        "upper": [2, 1] if affine else None,
        "Z_independent": affine,
        "function_keys": function_keys,
        "source_geometry_evidence": {
            "native_source_gate": native_ok,
            "C3_source_geometry": c3_ok,
            "native_archive_sha256": QPW_NATIVE_ARCHIVE_SHA256,
            "C3_archive_sha256": QPW_C3_ARCHIVE_SHA256,
            "derived_from_bound_source_metadata": affine,
            "arbitrary_passed_true_not_used": True,
        },
    }
    return qpw_clone(common), qpw_clone(common)


def qpw_generic_transport(parent):
    signed = parent.get("signed_density_algebra_result", {}) if isinstance(parent, dict) else {}
    pair = signed.get("paired_history_transport", {}) if isinstance(signed, dict) else {}
    all_zero = pair.get("all_zero") is True
    count = pair.get("check_count") == 20
    return {
        "passed": all_zero and count,
        "data_gate": all_zero and count,
        "all_zero": all_zero,
        "check_count": pair.get("check_count"),
        "source": "parent_result.signed_density_algebra_result.paired_history_transport",
        "data_only": True,
        "recomputed": False,
        "generic_family_replayed": False,
        "quiet_local_zero_does_not_reset_incoming": pair.get(
            "quiet_local_zero_does_not_reset_incoming"
        ),
    }


def qpw_source_entry(spec):
    entry = spec.get("source_entry", {}) if isinstance(spec, dict) else {}
    return {
        "delta_not_zeroed": entry.get("delta_not_zeroed") is True,
        "unknown_entry_gap_retained": entry.get("unknown_entry_gap_retained") is True,
        "values_uncomputed": entry.get("values_uncomputed") is True,
        "entry_gap_zero_or_nonzero_asserted": False,
        "predecessor_inlet_retained": True,
        "predecessor_inlet_zeroed": False,
    }


def qpw_compact_c5_route(spec):
    records = spec.get("creator_source_records", {}) if isinstance(spec, dict) else {}
    target = None
    for name, record in records.items() if isinstance(records, dict) else ():
        if str(name).replace("\\", "/").endswith(
            "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.py"
        ):
            target = record
            break
    required = (
        "self.g=g=outer.g",
        "self.incoming_correction=outer.Rc['correction']",
        "self.targets,self.N_scaled_incoming,self.joint=normalized_targets(g,self.incoming_correction,self.A,self.mu,self.N)",
    )
    text = target.get("source_text") if isinstance(target, dict) else None
    observed_sha = qpw_sha(text) if isinstance(text, str) else None
    endpoint = spec.get("compact_c5_endpoint_slice", {}) if isinstance(spec, dict) else {}
    fields = endpoint.get("fields") if isinstance(endpoint, dict) else None
    raw_rows = endpoint.get("compact_c5_raw_endpoint_rows") if isinstance(endpoint, dict) else None
    native_rows = endpoint.get("original_native_endpoint_rows") if isinstance(endpoint, dict) else None
    incoming_rows = endpoint.get("incoming_port_rows") if isinstance(endpoint, dict) else None
    incoming_role_rows = endpoint.get("incoming_port_row_roles") if isinstance(endpoint, dict) else None
    packet_endpoint_nodes = endpoint.get("packet_endpoint_nodes") if isinstance(endpoint, dict) else None
    clauses_ok = isinstance(text, str) and all(clause in text for clause in required)
    rows_ok = bool(
        fields == list(QPW_FIELDS)
        and isinstance(raw_rows, dict)
        and isinstance(native_rows, dict)
        and endpoint.get("compact_c5_archive_sha256") == QPW_PARENT_NATIVE_TERMINAL_ARCHIVE_SHA256
        and endpoint.get("original_native_archive_sha256") == QPW_NATIVE_ARCHIVE_SHA256
    )
    endpoint_row_pairs = {}
    if rows_ok:
        for field in QPW_FIELDS:
            left = raw_rows.get(field)
            right = native_rows.get(field)
            if not isinstance(left, list) or not isinstance(right, list) or len(left) != 4 or len(right) != 4:
                rows_ok = False
                break
            expected_ids = QPW_COMPACT_C5_ENDPOINT_IDS[field]
            endpoint_row_pairs[field] = []
            for order, (compact_row, native_row) in enumerate(zip(left, right)):
                if not isinstance(compact_row, dict) or not isinstance(native_row, dict):
                    rows_ok = False
                    break
                if (
                    compact_row.get("field") != field
                    or native_row.get("field") != field
                    or compact_row.get("ordinary_Z_order") != order
                    or native_row.get("ordinary_Z_order") != order
                    or compact_row.get("role") != "compact_C5_raw_endpoint_leaf"
                    or native_row.get("role") != "original_native_raw_endpoint_leaf"
                    or compact_row.get("node_id") != expected_ids[order]
                    or native_row.get("node_id") != expected_ids[order]
                    or compact_row.get("row") != native_row.get("row")
                ):
                    rows_ok = False
                    break
                for endpoint_row in (compact_row.get("row"), native_row.get("row")):
                    if (
                        not isinstance(endpoint_row, dict)
                        or endpoint_row.get("operation") != "function_substitution"
                        or endpoint_row.get("variable") != 250215
                        or endpoint_row.get("value") != 48
                        or endpoint_row.get("Z_independent_substitution") is not True
                        or endpoint_row.get("expression") != QPW_FINITE_HISTORY_IDS[field][order]
                    ):
                        rows_ok = False
                        break
                if not rows_ok:
                    break
                endpoint_row_pairs[field].append(
                    {"ordinary_Z_order": order, "node_id": expected_ids[order], "full_row_equal": True}
                )
            if not rows_ok:
                break
    incoming_ok = (
        isinstance(incoming_rows, dict)
        and incoming_rows.get("native") == QPW_NATIVE_INCOMING_IDS
        and incoming_rows.get("C3") == QPW_C3_INCOMING_IDS
    )
    incoming_roles_ok = True
    for role_name, expected_map, expected_role in (
        ("native", QPW_NATIVE_INCOMING_IDS, "native_quiet_incoming_raw_source_leaf"),
        ("C3", QPW_C3_INCOMING_IDS, "rebuilt_C3_quiet_incoming_source_leaf"),
    ):
        rows_by_field = incoming_role_rows.get(role_name) if isinstance(incoming_role_rows, dict) else None
        if not isinstance(rows_by_field, dict):
            incoming_roles_ok = False
            break
        for field in QPW_FIELDS:
            rows = rows_by_field.get(field)
            if not isinstance(rows, list) or len(rows) != 4:
                incoming_roles_ok = False
                break
            for order, row in enumerate(rows):
                if (
                    not isinstance(row, dict)
                    or row.get("field") != field
                    or row.get("ordinary_Z_order") != order
                    or row.get("role") != expected_role
                    or row.get("node_id") != expected_map[field][order]
                ):
                    incoming_roles_ok = False
                    break
            if not incoming_roles_ok:
                break
        if not incoming_roles_ok:
            break
    rows_ok = bool(rows_ok and incoming_ok and incoming_roles_ok)
    rows_ok = bool(
        rows_ok
        and packet_endpoint_nodes == QPW_COMPACT_C5_ENDPOINT_IDS
        and endpoint.get("comparison_deferred_to_frozen_kernel") is True
    )
    return {
        "bound": observed_sha == QPW_COMPACT_C5_SOURCE_SHA256 and clauses_ok and rows_ok,
        "source_path": target.get("path") if isinstance(target, dict) else None,
        "source_sha256": observed_sha,
        "expected_source_sha256": QPW_COMPACT_C5_SOURCE_SHA256,
        "fixed_constructor_clauses": list(required),
        "fixed_constructor_clauses_present": clauses_ok,
        "compact_c5_endpoint_rows": qpw_clone(raw_rows),
        "original_native_endpoint_rows": qpw_clone(native_rows),
        "endpoint_row_pairs": endpoint_row_pairs,
        "endpoint_fields": qpw_clone(fields),
        "incoming_port_rows": qpw_clone(incoming_rows),
        "incoming_port_rows_exact": incoming_ok,
        "incoming_port_row_roles": qpw_clone(incoming_role_rows),
        "incoming_port_row_roles_exact": incoming_roles_ok,
        "packet_endpoint_nodes": qpw_clone(packet_endpoint_nodes),
        "comparison_deferred_to_frozen_kernel": endpoint.get(
            "comparison_deferred_to_frozen_kernel"
        ) is True,
        "compact_c5_archive_sha256": endpoint.get("compact_c5_archive_sha256"),
        "original_native_archive_sha256": endpoint.get("original_native_archive_sha256"),
        "same_graph_constructor_route": True,
        "source_leaf_comparison_only": True,
        "archive_alias_claimed": False,
        "numeric_values_selected": False,
        "blockers": [] if (observed_sha == QPW_COMPACT_C5_SOURCE_SHA256 and clauses_ok and rows_ok) else ["same graph constructor route or endpoint leaves are not bound"],
    }


def qpw_entry_gap_bridge(parent, native, c3, compact_route=None):
    """Bind 20 explicit incoming and parent endpoint gap references as data."""
    parent_rows = parent.get("signed_history_gap_function_nodes") if isinstance(parent, dict) else None
    parent_native_ports = parent.get("native_history_source_ports", {}) if isinstance(parent, dict) else {}
    parent_c3_ports = parent.get("C3_rebuilt_history_source_ports", {}) if isinstance(parent, dict) else {}
    parent_native_scaled = parent_native_ports.get("N_scaled_ports") if isinstance(parent_native_ports, dict) else None
    parent_native_raw = parent_native_ports.get("raw_correction_ports") if isinstance(parent_native_ports, dict) else None
    parent_c3_terminal = parent_c3_ports.get("terminal_ports") if isinstance(parent_c3_ports, dict) else None
    native_rows = native.get("native_rows", {}) if isinstance(native, dict) else {}
    native_incoming = native_rows.get("incoming_first4")
    c3_incoming = c3.get("incoming_ports") if isinstance(c3, dict) else None
    if not isinstance(parent_rows, list) or len(parent_rows) != 20:
        return {"bound": False, "rows": [], "row_count": 0, "blockers": ["parent signed_history_gap_function_nodes must contain 20 rows"]}
    if not isinstance(compact_route, dict) or compact_route.get("bound") is not True:
        return {"bound": False, "rows": [], "row_count": 0, "blockers": ["same compact-C5 graph constructor route is not bound"]}
    if not isinstance(parent_native_scaled, dict) or not isinstance(parent_native_raw, dict) or not isinstance(parent_c3_terminal, dict):
        return {"bound": False, "rows": [], "row_count": 0, "blockers": ["parent native/C3 terminal source-port maps are missing"]}
    if parent_native_ports.get("archive_sha256") != QPW_PARENT_NATIVE_TERMINAL_ARCHIVE_SHA256 or parent_c3_ports.get("archive_sha256") != QPW_C3_ARCHIVE_SHA256:
        return {"bound": False, "rows": [], "row_count": 0, "blockers": ["parent native/C3 terminal archive SHA maps are not authenticated"]}
    if not isinstance(native_incoming, dict) or not isinstance(c3_incoming, dict):
        return {"bound": False, "rows": [], "row_count": 0, "blockers": ["native and C3 incoming port maps are missing"]}
    if native_incoming != QPW_NATIVE_INCOMING_IDS or c3_incoming != QPW_C3_INCOMING_IDS:
        return {"bound": False, "rows": [], "row_count": 0, "blockers": ["native/C3 incoming port maps do not exactly match the quiet source roles"]}
    rows = []
    blockers = []
    for parent_row in parent_rows:
        if not isinstance(parent_row, dict):
            blockers.append("parent gap row is not a mapping")
            continue
        field = parent_row.get("field")
        order = parent_row.get("ordinary_Z_order")
        if field not in QPW_FIELDS or order not in QPW_ORDERS:
            blockers.append("parent gap row field/order is outside quiet O3 scope")
            continue
        native_ids = native_incoming.get(field)
        c3_ids = c3_incoming.get(field)
        if not isinstance(native_ids, list) or not isinstance(c3_ids, list) or order >= len(native_ids) or order >= len(c3_ids):
            blockers.append("incoming source-port row is missing for " + str(field) + "/" + str(order))
            continue
        left = parent_row.get("left")
        right = parent_row.get("right")
        if not isinstance(left, dict) or not isinstance(right, dict):
            blockers.append("parent endpoint gap row lacks left/right source refs")
            continue
        if (
            not isinstance(parent_native_scaled.get(field), list)
            or not isinstance(parent_native_raw.get(field), list)
            or not isinstance(parent_c3_terminal.get(field), list)
            or order >= len(parent_native_scaled[field])
            or order >= len(parent_native_raw[field])
            or order >= len(parent_c3_terminal[field])
            or left.get("source_node") != parent_native_scaled[field][order]
            or left.get("raw_correction_node") != parent_native_raw[field][order]
            or right.get("source_node") != parent_c3_terminal[field][order]
        ):
            blockers.append("parent endpoint row disagrees with authenticated source-port maps")
            continue
        rows.append({
            "field": field,
            "ordinary_Z_order": order,
            "delta_in": {
                "formula": "Delta_in=N*native_quiet_incoming_rawID-rebuilt_C3_quiet_incomingID",
                "N_node": 18,
                "N_integer_materialized": False,
                "native_quiet_incoming_raw_node": native_ids[order],
                "native_archive_sha256": QPW_NATIVE_ARCHIVE_SHA256,
                "C3_quiet_incoming_node": c3_ids[order],
                "C3_archive_sha256": QPW_C3_ARCHIVE_SHA256,
                "gap_zero_or_nonzero_asserted": False,
                "values_uncomputed": True,
            },
            "delta_out": {
                "parent_gap_index": parent_row.get("gap_index"),
                "formula": QPW_TRANSPORT_SCOPE["formula"],
                "endpoint_coordinate": QPW_TRANSPORT_SCOPE["endpoint_coordinate"],
                "native_N_scaled_terminal_node": left.get("source_node"),
                "native_raw_correction_node": left.get("raw_correction_node"),
                "native_terminal_archive_sha256": left.get("archive_sha256"),
                "C3_terminal_node": right.get("source_node"),
                "C3_terminal_archive_sha256": right.get("archive_sha256"),
                "parent_gap_zero_claimed": parent_row.get("gap_zero_claimed") is True,
                "parent_values_uncomputed": parent_row.get("values_uncomputed") is True,
                "gap_zero_or_nonzero_asserted": False,
                "predecessor_inlet_retained": True,
            },
            "parent_row_data_only": True,
        })
    expected = {(field, order) for field in QPW_FIELDS for order in QPW_ORDERS}
    observed = {(row["field"], row["ordinary_Z_order"]) for row in rows}
    if observed != expected or len(rows) != 20:
        blockers.append("explicit quiet entry/endpoint bridge is not exactly 20 field/order rows")
    return {
        "bound": not blockers,
        "row_count": len(rows),
        "rows": rows,
        "blockers": blockers,
        "native_N_scaling_applied_once_inherited": parent.get("native_N_scaling_applied_once_bound") is True,
        "native_N_scaling_reproved": False,
        "normalization_reproved": False,
        "endpoint_coordinate": "2",
        "gap_zero_or_nonzero_asserted": False,
        "data_only": True,
        "same_graph_constructor_route_bound": True,
        "source_leaf_comparison_only": compact_route.get("source_leaf_comparison_only") is True,
        "archive_alias_claimed": compact_route.get("archive_alias_claimed") is True,
        "compact_c5_route": qpw_clone(compact_route),
        "parent_native_history_source_ports": {
            "N_scaled_ports": qpw_clone(parent_native_scaled),
            "raw_correction_ports": qpw_clone(parent_native_raw),
            "archive_sha256": parent_native_ports.get("archive_sha256"),
        },
        "parent_C3_rebuilt_history_source_ports": {
            "terminal_ports": qpw_clone(parent_c3_terminal),
            "archive_sha256": parent_c3_ports.get("archive_sha256"),
        },
    }


def qpw_math_contract(spec, native, c3, parent):
    native_geometry, c3_geometry = qpw_geometry(native, c3)
    native_gate = qpw_native_scalar_gate(native)
    c3_gate = qpw_c3_density_gate(c3)
    source_binding_passed = bool(
        native.get("status") == "bound"
        and c3.get("status") == "passed"
        and native.get("source_gate", {}).get("all_typed_bound") is True
        and c3_gate.get("passed") is True
    )
    source_gate = {
        "passed": bool(source_binding_passed and native_gate.get("root_gate") and c3_gate.get("passed")),
        "data_gate": bool(source_binding_passed and c3_gate.get("passed")),
        "native_rows": 20 if native.get("native_rows", {}).get("row_ids") else 0,
        "C3_rows": c3.get("zero_row_count", 0),
        "C3_fixed_zero_root_verified": c3_gate.get("C3_fixed_zero_root_verified") is True,
        "zero_local_densities_only": True,
        "zero_density_provenance": QPW_SOURCE_SCOPE["zero_density_provenance"],
        "source_bound_flat_branch_and_admission_certificate": True,
        "recomputed": False,
    }
    native_rates = native.get("source_gate", {}).get("rates")
    rates_text = {"m": "1", "h": "3/2", "k": "3/2", "e": "1", "p": "0"}
    qpw_require(native_rates == rates_text, "native source rate map changed")
    rates = {"m": [1, 1], "h": [3, 2], "k": [3, 2], "e": [1, 1], "p": [0, 1]}
    return {
        "schema": "compact_quiet_power_kernel_algebra_input_v4",
        "component": "compact_quiet_power_kernel_algebra",
        "source_binding": {
            "passed": source_binding_passed,
            "native_source_bound": native.get("status") == "bound",
            "C3_source_bound": c3.get("status") == "passed",
            "native_archive_sha256": QPW_NATIVE_ARCHIVE_SHA256,
            "C3_archive_sha256": QPW_C3_ARCHIVE_SHA256,
            "data_only": True,
            "generic_pair_replayed": False,
        },
        "native_source_geometry": native_geometry,
        "C3_source_geometry": c3_geometry,
        "rates": {"native": qpw_clone(rates), "C3": qpw_clone(rates)},
        "native_flat_scalar_source_gate": native_gate,
        "source_quiet_local_density_zero": source_gate,
        "cached_generic_paired_transport": qpw_generic_transport(parent),
        "source_entry": qpw_source_entry(spec),
    }


def qpw_certificate(native, c3, algebra, source, parent, bridge):
    source_entry = qpw_source_entry(source.get("input", {}))
    return {
        "schema": "compact_quiet_power_history_pair_certificate_v1",
        "passed": algebra.get("passed") is True,
        "native_source_bound": native.get("status") == "bound",
        "C3_source_bound": c3.get("status") == "passed",
        "native_zero_rows": 20 if native.get("rho_zero_conditional", {}).get("conditions_met") else 0,
        "C3_zero_rows": c3.get("zero_row_count", 0),
        "common_kernel_rows": algebra.get("stage_counts", {}).get("common_kernel_checks", 0),
        "endpoint_transport_rows": algebra.get("stage_counts", {}).get("source_conditional_transfer_rows", 0),
        "entry_gap_transport_rows": bridge.get("row_count", 0),
        "entry_gap_transport_row_count": bridge.get("row_count", 0),
        "entry_gap_endpoint_map": qpw_clone(bridge),
        "source_geometry": qpw_clone(algebra.get("geometry")),
        "source_integrand_kernel": QPW_KERNEL_SCOPE["source_integrand_kernel"],
        "memory_kernel": QPW_KERNEL_SCOPE["memory_kernel"],
        "transport_formula": QPW_TRANSPORT_SCOPE["formula"],
        "endpoint_coordinate": QPW_TRANSPORT_SCOPE["endpoint_coordinate"],
        "predecessor_inlet_retained": True,
        "unknown_entry_gap_retained": True,
        "entry_gap_zero_or_nonzero_asserted": False,
        "parent_generic_pair_data_only": True,
        "parent_result_raw_sha256": parent.get("_raw_sha256", QPW_PARENT_RESULT_SHA256),
        "values_uncomputed": True,
        "numeric_values_computed": False,
        "source_entry": source_entry,
        "actual_C3_partial_history_function_claimed": False,
        "same_graph_constructor_route_bound": bridge.get("same_graph_constructor_route_bound") is True,
        "source_leaf_comparison_only": bridge.get("source_leaf_comparison_only") is True,
    }


def qpw_consequence(native, c3, algebra, bridge):
    native_ok = native.get("status") == "bound"
    c3_ok = c3.get("status") == "passed"
    algebra_ok = algebra.get("passed") is True
    return {
        "schema": "compact_quiet_power_history_pair_consequence_v1",
        "passed": bool(native_ok and c3_ok and algebra_ok),
        "local_density_scope": "zero local O3 densities on coordinate interval [0,2] only",
        "ordinary_Z_function_scope": "identity on ordinary Z domain [-1,1], not samples",
        "source_integrand_kernel": "exp(-lambda*(2-s))",
        "memory_kernel": "exp(-2*lambda)",
        "transport_formula": "deltaT_out=exp(-2*lambda)*deltaT_in",
        "transport_scope": "endpoint coordinate 2 only",
        "predecessor_inlet_retained": True,
        "predecessor_inlet_zeroed": False,
        "unknown_entry_gap_retained": True,
        "entry_gap_zero_or_nonzero_asserted": False,
        "entry_gap_transport_rows": bridge.get("row_count", 0),
        "entry_gap_endpoint_map_bound": bridge.get("bound") is True,
        "same_graph_constructor_route_bound": bridge.get("same_graph_constructor_route_bound") is True,
        "source_leaf_comparison_only": bridge.get("source_leaf_comparison_only") is True,
        "actual_C3_partial_history_function_claimed": False,
        "actual_C3_C5_incoming_equality_claimed": False,
        "numeric_values_computed": False,
        "physical_completion_claimed": False,
        "full_goal_complete": False,
    }


def qpw_failure(source, parent, message, native=None, c3=None, events=None, bridge=None):
    native = native if isinstance(native, dict) else {}
    c3 = c3 if isinstance(c3, dict) else {}
    result = {
        "schema": QPW_RESULT_SCHEMA,
        "component": QPW_COMPONENT,
        "passed": False,
        "metadata_only": False,
        "source_only": True,
        "input_digest": source.get("input_digest") if isinstance(source, dict) else None,
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "source_packet_defined": False,
        "kernel_imported": True,
        "derive_entry_calls": 1,
        "derive_calls": 1,
        "entry_attempts": 1,
        "source_binding": {"passed": False, "values_uncomputed": True},
        "source_inspection_audit": {
            "original_source_evaluations": 0,
            "numeric_field_evaluations": 0,
            "old_source_consumer_calls": 0,
            "old_proof_functions_called": 0,
            "old_operator_calls": 0,
            "old_source_normalizer_calls": 0,
            "quadratures": 0,
            "temporal_recurrence_calls": 0,
            "new_native_quiet_density_rows": 0,
            "new_C3_quiet_density_rows": 0,
            "new_common_kernel_rows": 0,
            "new_signed_quiet_transport_rows": 0,
        },
        "blockers": [str(message)],
        "native_source_binding": qpw_clone(native),
        "C3_source_binding": qpw_clone(c3),
        "entry_gap_endpoint_map": qpw_clone(bridge or {"bound": False, "rows": [], "row_count": 0}),
        "diagnostic_journal": qpw_clone(events or []),
    }
    for key in QPW_NEW_FLAGS:
        result[key] = False
    result.update(qpw_parent_data(source, parent))
    qpw_apply_parent_metadata(result, source, parent)
    return result


def qpw_success(source, parent, native, c3, algebra, events, bridge=None):
    stage = algebra.get("stage_counts", {}) if isinstance(algebra, dict) else {}
    native_ok = bool(
        native.get("status") == "bound"
        and native.get("source_gate", {}).get("all_typed_bound") is True
        and native.get("rho_zero_conditional", {}).get("conditions_met") is True
    )
    c3_ok = bool(
        c3.get("status") == "passed"
        and c3.get("passed") is True
        and c3.get("zero_row_count") == 20
        and c3.get("quiet_proof", {}).get("passed") is True
    )
    algebra_ok = algebra.get("passed") is True
    bridge = bridge if isinstance(bridge, dict) else qpw_entry_gap_bridge(parent, native, c3)
    flags = {
        QPW_NEW_FLAGS[0]: native_ok,
        QPW_NEW_FLAGS[1]: c3_ok,
        QPW_NEW_FLAGS[2]: bool(algebra_ok and stage.get("common_kernel_checks") == 20),
        QPW_NEW_FLAGS[3]: bool(
            algebra_ok
            and stage.get("source_conditional_transfer_rows") == 20
            and bridge.get("bound") is True
        ),
    }
    passed = all(flags.values())
    result = {
        "schema": QPW_RESULT_SCHEMA,
        "component": QPW_COMPONENT,
        "passed": passed,
        "metadata_only": False,
        "source_only": True,
        "input_digest": source.get("input_digest"),
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "source_packet_defined": passed,
        "kernel_imported": True,
        "derive_entry_calls": 1,
        "derive_calls": 1,
        "entry_attempts": 1,
        "source_binding": {
            "passed": bool(native_ok and c3_ok),
            "native_source_bound": native_ok,
            "C3_source_bound": c3_ok,
            "values_uncomputed": True,
            "numeric_field_values_selected": False,
            "actual_C3_partial_history_claimed": False,
            "entry_gap_endpoint_map_bound": bridge.get("bound") is True,
            "entry_gap_transport_rows": bridge.get("row_count", 0),
            "entry_gap_zero_or_nonzero_asserted": False,
            "same_graph_constructor_route_bound": bridge.get("same_graph_constructor_route_bound") is True,
            "source_leaf_comparison_only": bridge.get("source_leaf_comparison_only") is True,
        },
        "quiet_power_scope": qpw_scope(),
        "native_source_binding": qpw_clone(native),
        "C3_source_binding": qpw_clone(c3),
        "quiet_power_history_pair_certificate": qpw_certificate(
            native, c3, algebra, source, parent, bridge
        ),
        "quiet_power_history_pair_consequence": qpw_consequence(
            native, c3, algebra, bridge
        ),
        "signed_history_gap_status": {
            "endpoint_only": True,
            "ordinary_Z_orders": list(QPW_ORDERS),
            "unknown_entry_gap_retained": True,
            "predecessor_inlet_retained": True,
            "predecessor_inlet_zeroed": False,
            "gap_zero_or_nonzero_asserted": False,
            "actual_C3_partial_history_claimed": False,
        },
        "source_inspection_audit": {
            "original_source_evaluations": 0,
            "numeric_field_evaluations": 0,
            "old_source_consumer_calls": 0,
            "old_proof_functions_called": 0,
            "old_operator_calls": 0,
            "old_source_normalizer_calls": 0,
            "quadratures": 0,
            "temporal_recurrence_calls": 0,
            "new_native_quiet_density_rows": 20 if native_ok else 0,
            "new_C3_quiet_density_rows": 20 if c3_ok else 0,
            "new_common_kernel_rows": stage.get("common_kernel_checks", 0),
            "new_signed_quiet_transport_rows": stage.get(
                "source_conditional_transfer_rows", 0
            ),
        },
        "entry_gap_endpoint_map": qpw_clone(bridge),
        "algebra_contract": qpw_clone(source.get("quiet_power_math_contract")),
        "algebra_result": qpw_clone(algebra),
        "diagnostic_journal": qpw_clone(events),
    }
    result.update(flags)
    result.update(qpw_parent_data(source, parent))
    qpw_apply_parent_metadata(result, source, parent)
    return result


def derive(source, journal=None):
    """Bind C3 then native source metadata, then derive the quiet pair once."""
    source = source if isinstance(source, dict) else {}
    spec = source.get("input") if isinstance(source.get("input"), dict) else {}
    parent = qpw_parent(source, spec)
    events = []
    native = {}
    c3 = {}
    try:
        qpw_require(spec.get("schema") == QPW_INPUT_SCHEMA, "quiet-power input schema changed")
        qpw_require(spec.get("component") == QPW_COMPONENT, "quiet-power input component changed")
        declared_parent_sha = source.get("parent_result_sha256")
        if not isinstance(declared_parent_sha, str):
            declared_parent_sha = spec.get("parent_result_metadata", {}).get("raw_sha256")
        qpw_require(declared_parent_sha == QPW_PARENT_RESULT_SHA256, "parent result raw SHA changed")
        qpw_require(parent.get("_raw_sha256") == QPW_PARENT_RESULT_SHA256, "parent result metadata marker is missing")
        qpw_require(spec.get("C3_receipt_raw_sha256") == QPW_C3_RECEIPT_SHA256, "C3 receipt raw SHA changed")
        qpw_require(spec.get("native_receipt", {}).get("sha256") == QPW_NATIVE_RECEIPT_SHA256, "native receipt SHA changed")
        qpw_require(spec.get("native_literal_addendum", {}).get("sha256") == QPW_NATIVE_LITERAL_ADDENDUM_SHA256, "native literal addendum SHA changed")
        qpw_require(spec.get("quiet_power_source_scope") == QPW_SOURCE_SCOPE, "quiet-power source scope changed")
        qpw_require(spec.get("quiet_power_common_kernel_scope") == QPW_KERNEL_SCOPE, "quiet-power kernel scope changed")
        qpw_require(spec.get("quiet_power_transport_scope") == QPW_TRANSPORT_SCOPE, "quiet-power transport scope changed")
        runtime_spec = qpw_clone(spec)
        runtime_parent = qpw_clone(parent)
        runtime_parent["_raw_sha256"] = QPW_PARENT_RESULT_SHA256
        runtime_spec["parent_result"] = runtime_parent
        qpw_emit(journal, "quiet_power_c3_binding_started", source_only=True)
        c3_journal = []
        c3 = bind_quiet_C3_source(runtime_spec, c3_journal)
        events.extend(qpw_clone(c3_journal))
        for event in c3_journal:
            qpw_emit(journal, event)
        qpw_emit(journal, "quiet_power_c3_binding_finished", passed=c3.get("passed") is True, status=c3.get("status"))
        if c3.get("status") != "passed" or c3.get("passed") is not True:
            return qpw_failure(source, parent, "C3 source binder is blocked", native, c3, events)
        qpw_emit(journal, "quiet_power_native_binding_started", source_only=True)
        native_journal = []
        native = bind_quiet_native_source(runtime_spec, native_journal)
        events.extend(qpw_clone(native_journal))
        for event in native_journal:
            qpw_emit(journal, event)
        qpw_emit(journal, "quiet_power_native_binding_finished", passed=native.get("status") == "bound", status=native.get("status"))
        if native.get("status") != "bound":
            return qpw_failure(source, parent, "native source binder is blocked", native, c3, events)
        compact_route = qpw_compact_c5_route(runtime_spec)
        if compact_route.get("bound") is not True:
            return qpw_failure(
                source,
                parent,
                "same compact-C5 graph constructor route is blocked",
                native,
                c3,
                events,
                {"bound": False, "rows": [], "row_count": 0, "compact_c5_route": compact_route},
            )
        bridge = qpw_entry_gap_bridge(parent, native, c3, compact_route)
        if bridge.get("bound") is not True:
            return qpw_failure(source, parent, "explicit quiet entry/endpoint bridge is blocked", native, c3, events, bridge)
        contract = qpw_math_contract(runtime_spec, native, c3, parent)
        source["quiet_power_math_contract"] = qpw_clone(contract)
        qpw_emit(journal, "quiet_power_math_binding_finished", source_bound=contract.get("source_binding", {}).get("passed") is True, values_uncomputed=True)
        if contract.get("source_binding", {}).get("passed") is not True:
            return qpw_failure(source, parent, "combined source binding is not passed", native, c3, events)
        def math_journal(event=None, *args, **fields):
            payload = {}
            if isinstance(event, dict):
                payload.update(event)
            elif event is not None:
                payload["event"] = event
            if args:
                payload["positional_args"] = list(args)
            payload.update(fields)
            events.append(qpw_clone(payload))
            qpw_emit(journal, payload)

        algebra = derive_quiet_power_kernel_algebra(contract, math_journal)
        result = qpw_success(source, parent, native, c3, algebra, events, bridge)
        qpw_emit(journal, "quiet_power_history_pair_defined", passed=result.get("passed") is True, source_values_uncomputed=True)
        return result
    except Exception as error:
        return qpw_failure(source, parent, error, native, c3, events)


__all__ = ("derive",)

