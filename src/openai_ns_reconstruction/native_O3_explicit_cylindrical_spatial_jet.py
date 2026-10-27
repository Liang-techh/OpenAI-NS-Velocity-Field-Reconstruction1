"""A source-only explicit cylindrical spatial-jet program for O3.

The input is the frozen ``native_O3_explicit_cylindrical_source_capsule_v1``
record.  Its branch fields and saved source rows are data, not executable
functions.  This module builds named formula trees, never imports the project,
calls a source provider, evaluates a source value, materializes N, or walks a
native expression graph.  Unknown or missing source rows are rejected.
"""

COMPONENT = "native_O3_explicit_cylindrical_spatial_jet"
SCHEMA = "native_O3_explicit_cylindrical_spatial_jet_kernel_result_v1"

SUCCESS_COUNTS = {
    "new_source_identity_family_calls": 1,
    "new_explicit_formula_family_calls": 1,
    "new_cylindrical_jet_family_calls": 1,
    "new_cartesian_composition_family_calls": 1,
    "new_mathematical_proofs": 1,
    "source_function_port_calls": 0,
    "saved_radial_derivative_port_calls": 0,
    "saved_flux_source_differential_port_calls": 0,
    "numeric_evaluation_calls": 0,
    "N_materialization_calls": 0,
    "old_math_calls": 0,
    "old_provider_calls": 0,
    "failed_family_calls": 0,
}

SUCCESS_FLAGS = {
    "conditional_closed_parameterization_frame": True,
    "conditional_cylindrical_velocity_formula": True,
    "conditional_physical_direction_operators": True,
    "conditional_cartesian_spatial_jet": True,
    "conditional_curl_and_divergence_formula": True,
    "conditional_pressure_gradient_formula": True,
    "conditional_volume_integrands_defined": True,
    "actual_source_Z_jets_admitted": False,
    "conditional_source_Z1_Z2_contract_consumed": True,
    "actual_incoming_x_independence_proved": False,
    "actual_parameter_admission": False,
    "global_axis_regularity_proved": False,
    "global_phase_seam_regularity_proved": False,
    "native_field_derivative_identity_proved": False,
    "physical_layer_completed": False,
    "numeric_evaluation": False,
    "full_goal_complete": False,
}

CALLBACK_EVENT_COUNTS = {
    "family_entered": 1,
    "source_join_summary": 1,
    "closed_parameterization_summary": 1,
    "cylindrical_formula_summary": 1,
    "physical_direction_summary": 1,
    "cartesian_jet_summary": 1,
    "curl_divergence_summary": 1,
    "pressure_volume_summary": 1,
    "family_completed": 1,
}

COUNT_EVENT_MAP = {
    "family_entered": [
        "new_source_identity_family_calls",
        "new_explicit_formula_family_calls",
        "new_cylindrical_jet_family_calls",
        "new_cartesian_composition_family_calls",
    ],
    "source_join_summary": [],
    "closed_parameterization_summary": [],
    "cylindrical_formula_summary": [],
    "physical_direction_summary": [],
    "cartesian_jet_summary": [],
    "curl_divergence_summary": [],
    "pressure_volume_summary": [],
    "family_completed": ["new_mathematical_proofs"],
}

NATIVE_SOURCE_FUNCTION_PORT = "native_source_function_port"
SAVED_RADIAL_DERIVATIVE_PORT = "saved_radial_derivative_port"
SAVED_FLUX_SOURCE_DIFFERENTIAL_PORT = "saved_flux_source_differential_port"


class KernelInputError(ValueError):
    """The frozen capsule does not satisfy the explicit source contract."""


def _mapping(value, label):
    if not isinstance(value, dict):
        raise KernelInputError("%s must be a mapping" % label)
    return value


def _key(value, name, label):
    if name not in value:
        raise KernelInputError("missing %s[%r]" % (label, name))
    return value[name]


def _literal(value, label):
    return {"kind": "opaque_source_atom", "label": label, "value": value}


def _exact_rational(numerator, denominator=1, label=None, source_node=None,
                    provenance=None):
    result = {
        "kind": "exact_rational",
        "num": numerator,
        "den": denominator,
        "numerator": numerator,
        "denominator": denominator,
        "numeric_evaluation": False,
    }
    if label is not None:
        result["label"] = label
    if source_node is not None:
        result["source_node"] = source_node
    if provenance is not None:
        result["provenance"] = provenance
    return result


def _native_coordinate_port(node, variable_name, role, binding=None):
    result = {
        "kind": "native_coordinate_port",
        "native_node": node,
        "variable_name": variable_name,
        "role": role,
        "numeric_evaluation": False,
    }
    if binding is not None:
        result["physical_composition_binding"] = binding
    return result


def _native_parameter_port(node, variable_name, role, raw_record=None,
                           binding=None, opaque=False):
    result = {
        "kind": "opaque_parameter_port" if opaque else "native_parameter_port",
        "native_node": node,
        "variable_name": variable_name,
        "role": role,
        "numeric_evaluation": False,
    }
    if raw_record is not None:
        result["raw_record"] = raw_record
    if binding is not None:
        result["physical_composition_binding"] = binding
    if opaque:
        result["materialized"] = False
    return result


def _free_variable_port(name, role, binding=None):
    result = {
        "kind": "free_physical_variable_port",
        "variable_name": name,
        "role": role,
        "numeric_evaluation": False,
    }
    if binding is not None:
        result["physical_composition_binding"] = binding
    return result


def _saved_expression_port(expression_id, saved_result_sha, role, source):
    return {
        "kind": "saved_expression_function_port",
        "expression_id": expression_id,
        "saved_result_sha": saved_result_sha,
        "role": role,
        "source": source,
        "numeric_evaluation": False,
    }


def _formula(operation, arguments, **attributes):
    """Create a literal formula node; ``operation`` is never a native id."""
    return {
        "kind": "explicit_formula",
        "operation": operation,
        "arguments": list(arguments),
        "attributes": attributes,
        "numeric_evaluation": False,
    }


def _add(a, b):
    return _formula("add", (a, b))


def _sub(a, b):
    return _formula("subtract", (a, b))


def _mul(a, b):
    return _formula("multiply", (a, b))


def _div(a, b):
    return _formula("divide", (a, b))


def _neg(a):
    return _formula("negative", (a,))


def _pow(a, b):
    return _formula("power", (a, b))


def _sqrt(a):
    return _formula("positive_sqrt", (a,))


def _exp(a):
    return _formula("exp", (a,))


def _exp_source(a, source_node):
    return _formula("exp", (a,), source_node=source_node,
                    source_operation="analytic_unary")


def _cos(a):
    return _formula("cos", (a,))


def _sin(a):
    return _formula("sin", (a,))


def _one():
    return _exact_rational(1, 1, "one")


def _two():
    return _exact_rational(2, 1, "two")


def _zero():
    return _exact_rational(0, 1, "typed_zero_source_contract",
                           provenance="typed P0_R source contract")


def _binding_summary(contract):
    native_ids = _mapping(_key(contract, "native_variable_ids",
                                "source_function_binding_contract"),
                          "source_function_binding_contract.native_variable_ids")
    return {
        "physical_arguments": list(_key(contract, "physical_arguments",
                                          "source_function_binding_contract")),
        "native_variable_ids": dict(native_ids),
        "recursive_bindings_apply_inside_all_source_function_references":
            _key(contract, "recursive_bindings_apply_inside_all_source_function_references",
                 "source_function_binding_contract"),
        "saved_positive_physical_inverse": _key(
            contract, "saved_positive_physical_inverse",
            "source_function_binding_contract"),
        "source_parameters_fixed_during_spatial_derivatives": _key(
            contract, "source_parameters_fixed_during_spatial_derivatives",
            "source_function_binding_contract"),
        "free_physical_composition_variables": {
            "tau": "physical_time_remaining",
            "theta": "physical_azimuth_not_native_phase",
        },
    }


def _source_atom(port, branch, field, order, value, binding):
    return {
        "kind": "atomic_source_function_value",
        "port": port,
        "branch": branch,
        "field": field,
        "ordinary_Z_order": order,
        "value_reference": value,
        "physical_composition_binding": binding,
        "ordinary_Z_source_data_only": True,
        "actual_Z_jet_admission": False,
        "numeric_evaluation": False,
    }


def _source_derivative(port, branch, field, direction, order, value,
                       producer, binding, source_operator=None):
    return {
        "kind": "atomic_source_derivative_value",
        "port": port,
        "branch": branch,
        "field": field,
        "direction": direction,
        "source_operator": source_operator if source_operator is not None else direction,
        "ordinary_Z_order": order,
        "value_reference": value,
        "producer_reference": producer,
        "physical_composition_binding": binding,
        "ordinary_Z_source_data_only": True,
        "actual_Z_jet_admission": False,
        "numeric_evaluation": False,
    }


def _emit(callback, event, row):
    if callback is not None:
        callback(event, row)


def _validate_parameter_records(root):
    records = _mapping(_key(root, "parameter_constructor_records", "capsule"),
                       "parameter_constructor_records")
    if len(records) != 46:
        raise KernelInputError("parameter constructor record count is not 46")
    required = ("1", "2", "3", "4", "5", "8", "10", "13", "14", "18",
                "28", "643", "2083", "4194", "35560", "65604", "79968",
                "81513", "81514", "84529", "84533", "85743", "205154",
                "25220")
    for record_id in required:
        _mapping(_key(records, record_id, "parameter_constructor_records"),
                 "parameter_constructor_records[%s]" % record_id)
    for record_id, name in ((8, "hbB"), (10, "sc"), (35560, "true_core_logC")):
        record = records[str(record_id)]
        if _key(record, "operation", "record %s" % record_id) != "current_original_source_parameter":
            raise KernelInputError("source parameter record %s has wrong operation" % record_id)
        if _key(record, "name", "record %s" % record_id) != name:
            raise KernelInputError("source parameter record %s has wrong name" % record_id)
        if _key(record, "Z_independent", "record %s" % record_id) is not True:
            raise KernelInputError("source parameter record %s lost its fixed-Z premise" % record_id)
    # This is a finite literal-pattern join, not a graph walk.
    if _key(records["2"], "numerator", "record 2") != 40:
        raise KernelInputError("raw atom 2 is not the exact 40 record")
    if _key(records["3"], "operation", "record 3") != "analytic_unary":
        raise KernelInputError("raw atom 3 is not exp(40)")
    if _key(records["4"], "numerator", "record 4") != 11:
        raise KernelInputError("raw atom 4 is not the exact 11 record")
    if _key(records["5"], "arguments", "record 5") != [3, 4]:
        raise KernelInputError("raw atom 5 does not join exp40 and 11")
    rational_records = {
        "1": (1, 1), "2": (40, 1), "4": (11, 1), "13": (4, 1),
        "28": (1, 2), "643": (1000, 1), "2083": (10, 1),
        "4194": (30, 1), "79968": (110, 1), "84529": (-7, 1),
        "205154": (17, 1),
    }
    for record_id, (numerator, denominator) in rational_records.items():
        record = records[record_id]
        if _key(record, "operation", "record %s" % record_id) != "exact_rational":
            raise KernelInputError("raw atom %s is not exact_rational" % record_id)
        if (_key(record, "numerator", "record %s" % record_id) != numerator or
                _key(record, "denominator", "record %s" % record_id) != denominator):
            raise KernelInputError("raw atom %s has unexpected rational literal" % record_id)
    def operation(record_id, expected, arguments=None, name=None, argument=None):
        record = records[str(record_id)]
        if _key(record, "operation", "record %s" % record_id) != expected:
            raise KernelInputError("raw constructor record %s has wrong operation" % record_id)
        if arguments is not None and _key(record, "arguments", "record %s" % record_id) != arguments:
            raise KernelInputError("raw constructor record %s has wrong children" % record_id)
        if name is not None and _key(record, "name", "record %s" % record_id) != name:
            raise KernelInputError("raw constructor record %s has wrong unary name" % record_id)
        if argument is not None and _key(record, "argument", "record %s" % record_id) != argument:
            raise KernelInputError("raw constructor record %s has wrong unary child" % record_id)
    operation(3, "analytic_unary", name="exp", argument=2)
    operation(5, "sum", arguments=[3, 4])
    operation(14, "product", arguments=[13, 5])
    operation(25218, "sum", arguments=[14, 4194])
    operation(25219, "negative", argument=25218)
    operation(25220, "analytic_unary", name="exp", argument=25219)
    operation(33572, "analytic_unary", name="exp", argument=5)
    operation(35562, "sum", arguments=[14, 643])
    operation(35564, "negative", argument=35562)
    operation(35565, "analytic_unary", name="exp", argument=35564)
    operation(65604, "product", arguments=[13, 35565])
    operation(81513, "sum", arguments=[35560, 5])
    operation(81514, "product", arguments=[2083, 81513])
    operation(84530, "sum", arguments=[84529, 84528])
    operation(84531, "sum", arguments=[81514, 84530])
    operation(84532, "analytic_unary", name="exp", argument=84531)
    operation(84533, "product", arguments=[79968, 84532])
    operation(241542, "sum", arguments=[205154, 3, 241536])
    operation(241786, "analytic_unary", name="exp", argument=241542)
    operation(241787, "product", arguments=[85743, 241786])
    operation(24954, "product", arguments=[8, 10, 28])
    operation(85743, "function_substitution")
    substitution = records["85743"]
    if (_key(substitution, "expression", "record 85743") != 84533 or
            _key(substitution, "variable", "record 85743") != 84528 or
            _key(substitution, "value", "record 85743") != 1 or
            _key(substitution, "Z_independent_substitution", "record 85743") is not True):
        raise KernelInputError("Rm substitution record is not the exact source join")
    if _key(substitution, "value", "record 85743") != _key(records["1"],
                                                               "numerator", "record 1"):
        raise KernelInputError("Rm substitution value is not raw node 1")
    if _key(records["84529"], "numerator", "record 84529") != -7:
        raise KernelInputError("Rm constructor exponent does not retain -7")
    core = records["35560"]
    if _key(core, "parameter_enclosure_not_selected_value", "record 35560") is not True:
        raise KernelInputError("core_logC must remain an unselected opaque parameter")
    if _key(core, "exact_fixed_source_parameter", "record 35560") is not True:
        raise KernelInputError("core_logC is not an exact fixed source parameter")
    n_record = records["18"]
    if _key(n_record, "exact_positive_integer", "record 18") is not True:
        raise KernelInputError("N record is not an opaque positive integer")
    if _key(n_record, "outer_integer_and_inner_power_not_materialized",
            "record 18") is not True:
        raise KernelInputError("N materialization is not forbidden by the capsule")
    return records


def _validate_capsule(capsule):
    root = _mapping(capsule, "capsule")
    if _key(root, "schema", "capsule") != "native_O3_explicit_cylindrical_source_capsule_v1":
        raise KernelInputError("unexpected explicit cylindrical source capsule")
    required = (
        "branches", "accepted_radial_result", "accepted_flux_source_rows",
        "accepted_exterior_payload", "saved_geometry_row",
        "physical_source_capsule", "parameter_ports",
        "parameter_frame_domain", "parameter_constructor_records",
        "source_function_binding_contract", "source_jet_domain",
        "scope", "preservation", "source_pins",
    )
    for name in required:
        _key(root, name, "capsule")
    branches = _mapping(root["branches"], "branches")
    for branch_name in ("background", "modified"):
        branch = _mapping(_key(branches, branch_name, "branches"),
                          "branches[%s]" % branch_name)
        fields = _mapping(_key(branch, "fields", "branches[%s]" % branch_name),
                          "branches[%s].fields" % branch_name)
        for field in ("E", "V", "m", "p", "P0", "pressure"):
            rows = _key(fields, field, "branches[%s].fields" % branch_name)
            if not isinstance(rows, list) or len(rows) != 5:
                raise KernelInputError("%s.%s must retain five ordinary-Z rows" %
                                       (branch_name, field))
        targets = _mapping(_key(branch, "radial_target_names",
                                "branches[%s]" % branch_name),
                           "branches[%s].radial_target_names" % branch_name)
        expected = {"background": {"E": "E_background", "V": "V_background"},
                    "modified": {"E": "E_N", "V": "V_N"}}[branch_name]
        if targets != expected:
            raise KernelInputError("radial target names are not exact for %s" % branch_name)
        rows = _mapping(_key(branch, "flux_source_rows",
                             "branches[%s]" % branch_name),
                        "branches[%s].flux_source_rows" % branch_name)
        for field, orders in (("m", ("0", "1")), ("p", ("0",)),
                              ("P0", ("0",))):
            field_rows = _mapping(_key(rows, field, "flux rows"),
                                  "flux rows %s.%s" % (branch_name, field))
            for order in orders:
                _key(field_rows, order, "flux rows %s.%s" % (branch_name, field))
        if _key(branch, "source_values_unevaluated",
                "branches[%s]" % branch_name) is not True:
            raise KernelInputError("branch source values were evaluated")
        exterior = _mapping(_key(branch, "exterior_branch",
                                 "branches[%s]" % branch_name),
                            "branches[%s].exterior_branch" % branch_name)
        if _key(exterior, "P0_used_once", "exterior_branch") is not True:
            raise KernelInputError("P0 usage certificate missing")
        if _key(exterior, "incoming_history_reset", "exterior_branch") is not False:
            raise KernelInputError("incoming histories were reset")
    _validate_parameter_records(root)
    ports = _mapping(root["parameter_ports"], "parameter_ports")
    expected_ports = {"N": 18, "Pstar": 33572, "R": 241787, "Ra": 65604,
                      "Rm": 85743, "Z": 20, "core_logC": 35560,
                      "delta": 25220, "hbB": 8, "native_x": 241536,
                      "phase_origin": 24954, "sc": 10}
    if ports != expected_ports:
        raise KernelInputError("parameter port table is not the admitted explicit table")
    domain = _mapping(root["parameter_frame_domain"], "parameter_frame_domain")
    if _key(domain, "closed_parameter_join_executed", "parameter_frame_domain") is not False:
        raise KernelInputError("parameter join was already executed")
    if _key(domain, "core_logC_value_unselected", "parameter_frame_domain") is not True:
        raise KernelInputError("core_logC value was selected")
    if _key(domain, "unconditional_source_parameter_admission",
            "parameter_frame_domain") is not False:
        raise KernelInputError("unconditional parameter admission is not false")
    contract = _mapping(root["source_function_binding_contract"],
                        "source_function_binding_contract")
    binding = _binding_summary(contract)
    if binding["native_variable_ids"] != {"R": 241787, "Z": 20,
                                           "native_x": 241536}:
        raise KernelInputError("native source binding ids are not exact")
    if binding["recursive_bindings_apply_inside_all_source_function_references"] is not True:
        raise KernelInputError("recursive physical composition binding is missing")
    geometry = _mapping(root["saved_geometry_row"], "saved_geometry_row")
    for name in ("r", "z", "J2", "inverse"):
        _key(geometry, name, "saved_geometry_row")
    exterior_payload = _mapping(root["accepted_exterior_payload"],
                                "accepted_exterior_payload")
    volume = _mapping(_key(exterior_payload, "physical_volume",
                           "accepted_exterior_payload"),
                      "accepted_exterior_payload.physical_volume")
    for name in ("profile_R_Z_theta_absolute_measure",
                 "native_x_Z_theta_absolute_measure"):
        _key(volume, name, "accepted exterior physical volume")
    jet_domain = _mapping(root["source_jet_domain"], "source_jet_domain")
    for name in ("tau_positive", "R_positive", "absolute_Z_less_than_one",
                 "off_phase_seams", "conditional_C2_regular_source_functions",
                 "ordinary_Z_at_fixed_source_parameters"):
        if _key(jet_domain, name, "source_jet_domain") is not True:
            raise KernelInputError("source jet premise %s is not retained" % name)
    if _key(jet_domain, "unconditional_ordinary_Z_jet_admission",
            "source_jet_domain") is not False:
        raise KernelInputError("ordinary-Z source rows were unconditionally admitted")
    radial_payload = _mapping(_key(root["accepted_radial_result"], "result",
                                   "accepted_radial_result"),
                             "accepted_radial_result.result")
    radial_payload = _mapping(_key(radial_payload, "payload",
                                   "accepted_radial_result.result"),
                             "accepted_radial_result.result.payload")
    targets = _key(radial_payload, "targets", "accepted radial payload")
    if not isinstance(targets, list):
        raise KernelInputError("accepted radial targets are not a list")
    for name in ("E_background", "V_background", "E_N", "V_N"):
        matches = [target for target in targets
                   if isinstance(target, dict) and "name" in target
                   and target["name"] == name]
        if len(matches) != 1:
            raise KernelInputError("radial target %s is missing or ambiguous" % name)
        _key(matches[0], "tangent_ref", "radial target %s" % name)
    return root, binding


def _target(root, target_name):
    radial = root["accepted_radial_result"]["result"]["payload"]["targets"]
    matches = [target for target in radial
               if isinstance(target, dict) and "name" in target
               and target["name"] == target_name]
    if len(matches) != 1:
        raise KernelInputError("radial target %s is missing or ambiguous" % target_name)
    return matches[0]


def _field(root, branch_name, name, order):
    rows = root["branches"][branch_name]["fields"][name]
    if not isinstance(order, int) or order < 0 or order >= len(rows):
        raise KernelInputError("ordinary-Z row %s.%s[%s] is unavailable" %
                               (branch_name, name, order))
    return rows[order]


def _flux_row(root, branch_name, field, order):
    rows = root["branches"][branch_name]["flux_source_rows"][field]
    if str(order) not in rows:
        raise KernelInputError("flux source row %s.%s[%s] is unavailable" %
                               (branch_name, field, order))
    return rows[str(order)]


def _source_values(root, branch_name, binding):
    source_port = NATIVE_SOURCE_FUNCTION_PORT
    def value(field, order=0):
        return _source_atom(source_port, branch_name, field, order,
                            _field(root, branch_name, field, order), binding)
    m = value("m")
    v = value("V")
    e = value("E")
    p_hist = value("p")
    p0 = value("P0")
    m_z = value("m", 1)
    m_zz = value("m", 2)
    e_z = value("E", 1)
    v_z = value("V", 1)
    p_hist_z = value("p", 1)
    p0_z = value("P0", 1)
    r_ref = _native_coordinate_port(
        root["parameter_ports"]["R"], "R", "positive similarity profile radius", binding)
    target_e = _target(root, root["branches"][branch_name]["radial_target_names"]["E"])
    target_v = _target(root, root["branches"][branch_name]["radial_target_names"]["V"])
    e_dy = _source_derivative(
        SAVED_RADIAL_DERIVATIVE_PORT, branch_name, "E", "native_y", 0,
        _key(target_e, "tangent_ref", "radial target %s" % target_e["name"]),
        target_e, binding, source_operator="D_y")
    v_dy = _source_derivative(
        SAVED_RADIAL_DERIVATIVE_PORT, branch_name, "V", "native_y", 0,
        _key(target_v, "tangent_ref", "radial target %s" % target_v["name"]),
        target_v, binding, source_operator="D_y")
    e_r = _div(e_dy, r_ref)
    v_r = _div(v_dy, r_ref)
    # The accepted FToC row supplies provenance for the explicit mixed rule.
    m_flux = _flux_row(root, branch_name, "m", 1)
    m_zr = {
        "kind": "accepted_flux_FToC_j1_formula",
        "formula": _div(_sub(v_z, m_z), r_ref),
        "flux_row": m_flux,
        "identity": "m_ZR=(V_Z-m_Z)/R",
        "numeric_evaluation": False,
    }
    p_flux = _flux_row(root, branch_name, "p", 0)
    p_hist_r = _div(_source_derivative(SAVED_FLUX_SOURCE_DIFFERENTIAL_PORT,
                                       branch_name, "p", "native_y", 0,
                                       _key(p_flux, "source_expression",
                                            "flux row %s.p.0" % branch_name),
                                       p_flux, binding, source_operator="D_y"), r_ref)
    p0_flux = _flux_row(root, branch_name, "P0", 0)
    p0_r = _div(_source_derivative(SAVED_FLUX_SOURCE_DIFFERENTIAL_PORT,
                              branch_name, "P0", "native_y", 0,
                              _key(p0_flux, "source_expression",
                                   "flux row %s.P0.0" % branch_name),
                              p0_flux, binding, source_operator="D_y"), r_ref)
    incoming = {
        "m": _literal(_field(root, branch_name, "m", 0),
                       "incoming_%s_m_preserved" % branch_name),
        "p": _literal(_field(root, branch_name, "p", 0),
                       "incoming_%s_p_preserved" % branch_name),
    }
    physical_rows = root["physical_source_capsule"][
        "background_rows" if branch_name == "background" else "modified_rows"]
    incoming["e"] = _literal(physical_rows["e"][0],
                              "incoming_%s_e_preserved" % branch_name)
    return {
        "m": m, "V": v, "E": e, "p_hist": p_hist, "P0": p0,
        "m_Z": m_z, "m_ZZ": m_zz, "E_Z": e_z, "V_Z": v_z,
        "p_hist_Z": p_hist_z, "P0_Z": p0_z, "E_R": e_r, "V_R": v_r,
        "m_ZR": m_zr, "p_hist_R": p_hist_r, "P0_R": p0_r,
        "incoming_preserved": incoming,
        "field_rows_preserved": root["branches"][branch_name]["fields"],
        "flux_rows_preserved": root["branches"][branch_name]["flux_source_rows"],
        "exterior_data": root["branches"][branch_name]["exterior_branch"],
    }


def _closed_parameterization(root):
    records = root["parameter_constructor_records"]
    def exact(record_id, label):
        record = _mapping(records[str(record_id)],
                          "parameter_constructor_records[%s]" % record_id)
        return _exact_rational(
            _key(record, "numerator", "record %s" % record_id),
            _key(record, "denominator", "record %s" % record_id),
            label=label, source_node=record_id)

    def parameter(name, record_id, role, opaque=False):
        return _native_parameter_port(
            record_id, name, role, raw_record=records[str(record_id)],
            opaque=opaque)

    exp40 = _exp_source(exact(2, "exp40_argument"), 3)
    H = _formula("add", (exp40, exact(4, "eleven")),
                 source_node=5, source_operation="sum")
    core = parameter("true_core_logC", 35560, "unselected fixed source parameter")
    hb_b = parameter("hbB", 8, "fixed source parameter")
    sc = parameter("sc", 10, "fixed source parameter")
    half = _div(_one(), _two())
    six = _exact_rational(
        6, 1, label="closed_Rm_constant_6",
        provenance={"raw_nodes": [84529, 84528],
                    "operation": "source_substitution_minus7_plus_fraction",
                    "fraction_record": 84528, "fraction_value": 1})
    rm_exponent = _sub(_mul(exact(2083, "ten"), _add(core, H)), six)
    rm = _formula(
        "multiply",
        (exact(79968, "one_hundred_ten"), _exp_source(rm_exponent, 84532)),
        source_node=84533, source_operation="product",
        substitution_source_node=85743)
    ra_exponent = _neg(_add(_mul(exact(13, "four"), H),
                            exact(643, "one_thousand")))
    ra = _formula("multiply", (exact(13, "four"),
                                _exp_source(ra_exponent, 35565)),
                  source_node=65604, source_operation="product")
    # Raw node 5 is H; raw node 33572 is exp(H), not an evaluated constant 5.
    pstar = _exp_source(H, 33572)
    delta_exponent = _neg(_add(_mul(exact(13, "four"), H),
                               exact(4194, "thirty")))
    delta = _exp_source(delta_exponent, 25220)
    phase = _formula("multiply", (hb_b, sc, half),
                     source_node=24954, source_operation="product")
    n_atom = parameter("opaque_N", 18, "opaque positive integer source premise",
                       opaque=True)
    native_x = _native_coordinate_port(
        241536, "native_x", "O3 log-profile-radius chart variable")
    native_R = _formula(
        "multiply",
        (rm, _exp_source(_add(_add(exact(205154, "seventeen"), exp40),
                              native_x), 241786)),
        source_node=241787, source_operation="product")
    formulas = {
        "Rm": {"text": "Rm=110*exp(10*(true_core_logC+exp40+11)-6)",
               "formula": rm},
        "Ra": {"text": "Ra=4*exp(-(4*(exp40+11)+1000))", "formula": ra},
        "phase_origin": {"text": "phase_origin=hbB*sc*(1/2)", "formula": phase},
        "Pstar": {"text": "Pstar=exp(H), H=exp(40)+11 (native node5 is H)", "formula": pstar},
        "delta": {"text": "delta=exp(-(4*(exp40+11)+30))", "formula": delta},
        "N": {"text": "N=2^(2^J)", "formula": n_atom,
              "materialized": False},
    }
    return {
        "kind": "conditional_closed_source_frame",
        "raw_record_count": 46,
        "raw_records": records,
        "fixed_named_atoms": {"hbB": 8, "sc": 10, "true_core_logC": 35560},
        "opaque_N_atom": {"node": 18, "materialized": False},
        "formulas": formulas,
        "binding_certificate": {
            "Rm_raw_pattern": "110*exp(10*(true_core_logC+exp40+11)-7+fraction)",
            "Rm_substitution": {"record": 85743, "fraction_record": 84528,
                                "fraction_source_node": 1, "fraction_value": 1,
                                "closed_constant": -6},
            "Ra_raw_children": [13, 35565],
            "delta_raw_children": [14, 4194],
            "phase_origin_raw_children": [8, 10, 28],
            "Pstar_raw_children": [33572, 5],
            "core_logC_is_opaque_fixed_input": True,
            "N_is_opaque_and_unmaterialized": True,
        },
        "unselected_core_logC": True,
        "unconditional_source_parameter_admission": False,
        "numeric_evaluation": False,
        "native_R": {
            "text": "R=Rm*exp(17+exp40+native_x)",
            "formula": native_R,
        },
    }


def _geometry(root, binding):
    ports = root["parameter_ports"]
    tau = _free_variable_port("tau", "physical_time_remaining", binding)
    theta = _free_variable_port("theta", "physical_azimuth_not_native_phase", binding)
    z = _native_coordinate_port(ports["Z"], "Z", "ordinary axial source variable", binding)
    R = _native_coordinate_port(ports["R"], "R", "positive similarity profile radius", binding)
    delta = _native_parameter_port(ports["delta"], "delta",
                                   "fixed original similarity exponent", binding=binding)
    pstar = _native_parameter_port(ports["Pstar"], "Pstar",
                                   "fixed original Pstar parameter", binding=binding)
    one, two = _one(), _two()
    d = _sub(one, _mul(z, z))
    L = _sub(one, _mul(delta, _mul(z, z)))
    a = _sub(one, delta)
    lam = _sqrt(_div(tau, d))
    r = _mul(lam, _sqrt(_mul(two, R)))
    kappa = _pow(lam, a)
    C = _mul(pstar, _pow(lam, _neg(_add(one, delta))))
    B = _div(_mul(pstar, R), _mul(r, L))
    saved = root["saved_geometry_row"]
    volume = root["accepted_exterior_payload"]["physical_volume"]
    exterior_sha = _key(root["source_pins"]["exterior"], "sha256",
                        "source_pins.exterior")
    return {
        "tau": tau, "theta": theta, "Z": z, "R": R, "delta": delta,
        "Pstar": pstar, "d": d, "L": L, "a": a, "lambda": lam,
        "r": r, "kappa": kappa, "C": C, "B": B,
        "saved_r": _saved_expression_port(
            _key(_mapping(saved["r"], "saved_geometry_row.r"),
                 "value", "saved_geometry_row.r"),
            exterior_sha, "saved_geometry_r",
            "callback_events[event=geometry_dual].row.r.value"),
        "saved_z": _saved_expression_port(
            _key(_mapping(saved["z"], "saved_geometry_row.z"),
                 "value", "saved_geometry_row.z"),
            exterior_sha, "saved_geometry_z",
            "callback_events[event=geometry_dual].row.z.value"),
        "saved_J2": _saved_expression_port(
            saved["J2"], exterior_sha, "saved_geometry_J2",
            "callback_events[event=geometry_dual].row.J2"),
        "saved_inverse": saved["inverse"],
        "profile_measure": _saved_expression_port(
            volume["profile_R_Z_theta_absolute_measure"], exterior_sha,
            "saved_profile_R_Z_theta_measure_201",
            "accepted_exterior_payload.physical_volume.profile_R_Z_theta_absolute_measure"),
        "native_x_measure": _saved_expression_port(
            volume["native_x_Z_theta_absolute_measure"], exterior_sha,
            "saved_native_x_extra_R_measure_202",
            "accepted_exterior_payload.physical_volume.native_x_Z_theta_absolute_measure"),
    }


def _cylindrical(profile, g):
    S = _sub(_add(_mul(g["d"], profile["m_Z"]),
                  _mul(g["a"], _mul(g["Z"], profile["m"]))),
             _mul(_two(), _mul(g["Z"], profile["V"])))
    return {"S": S, "u_r": _neg(_mul(g["B"], S)),
            "u_theta": _mul(g["C"], profile["E"]),
            "u_z": _mul(g["C"], profile["V"])}


def _partials(profile, g, cyl):
    d, L, a, z, R, B, C = (g[k] for k in ("d", "L", "a", "Z", "R", "B", "C"))
    m, mz, V = profile["m"], profile["m_Z"], profile["V"]
    m_R = _div(_sub(V, m), R)
    mz_R = _div(_sub(profile["V_Z"], mz), R)
    S_R = _sub(_add(_mul(d, mz_R), _mul(a, _mul(z, m_R))),
                _mul(_two(), _mul(z, profile["V_R"])))
    S_Z = _sub(_sub(_add(_mul(d, profile["m_ZZ"]), _mul(a, m)),
                    _mul(_add(_one(), g["delta"]), _mul(z, mz))),
                _add(_mul(_two(), V), _mul(_two(), _mul(z, profile["V_Z"]))))
    B_R = _div(B, _mul(_two(), R))
    B_Z = _mul(B, _add(_neg(_div(z, d)),
                       _div(_mul(_two(), _mul(g["delta"], z)), L)))
    C_R = _zero()
    C_Z = _neg(_div(_mul(_add(_one(), g["delta"]), _mul(z, C)), d))
    return {
        "u_r_R": _neg(_add(_mul(B_R, cyl["S"]), _mul(B, S_R))),
        "u_r_Z": _neg(_add(_mul(B_Z, cyl["S"]), _mul(B, S_Z))),
        "u_theta_R": _add(_mul(C_R, profile["E"]), _mul(C, profile["E_R"])),
        "u_theta_Z": _add(_mul(C_Z, profile["E"]), _mul(C, profile["E_Z"])),
        "u_z_R": _add(_mul(C_R, profile["V"]), _mul(C, profile["V_R"])),
        "u_z_Z": _add(_mul(C_Z, profile["V"]), _mul(C, profile["V_Z"])),
        "m_R": m_R, "m_ZR": mz_R, "S_R": S_R, "S_Z": S_Z,
        "B_R": B_R, "B_Z": B_Z, "C_R": C_R, "C_Z": C_Z,
    }


def _physical_directions(partials, g):
    radial = _div(_mul(_two(), g["R"]), g["r"])
    axial_factor = _div(g["d"], _mul(g["kappa"], g["L"]))
    axial_shift = _div(_mul(_two(), _mul(g["R"], g["Z"])), g["d"])
    out = {}
    for name in ("u_r", "u_theta", "u_z"):
        q_R, q_Z = partials["%s_R" % name], partials["%s_Z" % name]
        out["%s_r" % name] = _mul(radial, q_R)
        out["%s_z" % name] = _mul(axial_factor, _sub(q_Z, _mul(axial_shift, q_R)))
    out["operators"] = {
        "D_r": "(2R/r)*D_R",
        "D_z": "d/(kappa*L)*(D_Z-(2RZ/d)*D_R)",
        "at_fixed": "physical time and z for D_r; physical time and r for D_z",
    }
    return out


def _cartesian(cyl, directions, g):
    c, s, r = _cos(g["theta"]), _sin(g["theta"]), g["r"]
    u = _sub(_mul(c, cyl["u_r"]), _mul(s, cyl["u_theta"]))
    v = _add(_mul(s, cyl["u_r"]), _mul(c, cyl["u_theta"]))
    w = cyl["u_z"]
    theta = {"u": _neg(v), "v": u, "w": _zero()}
    radial = {
        "u": _sub(_mul(c, directions["u_r_r"]),
                   _mul(s, directions["u_theta_r"])),
        "v": _add(_mul(s, directions["u_r_r"]),
                   _mul(c, directions["u_theta_r"])),
        "w": directions["u_z_r"],
    }
    axial = {
        "u": _sub(_mul(c, directions["u_r_z"]),
                   _mul(s, directions["u_theta_z"])),
        "v": _add(_mul(s, directions["u_r_z"]),
                   _mul(c, directions["u_theta_z"])),
        "w": directions["u_z_z"],
    }
    gradient = {}
    for component in ("u", "v", "w"):
        q_r, q_z = radial[component], axial[component]
        q_theta = theta[component]
        gradient[component] = {
            "X": _sub(_mul(c, q_r), _div(_mul(s, q_theta), r)),
            "Y": _add(_mul(s, q_r), _div(_mul(c, q_theta), r)),
            "Z": q_z,
        }
    return {
        "components": {"u": u, "v": v, "w": w},
        "theta_derivatives": theta,
        "gradient": gradient,
        "basis_derivative_receipt": {
            "D_theta_u": "-v", "D_theta_v": "u", "D_theta_w": "0",
            "D_X": "c*D_r-s/r*D_theta", "D_Y": "s*D_r+c/r*D_theta",
        },
    }


def _curl_divergence(cyl, directions, g):
    omega_r = _neg(directions["u_theta_z"])
    omega_theta = _sub(directions["u_r_z"], directions["u_z_r"])
    omega_z = _add(directions["u_theta_r"], _div(cyl["u_theta"], g["r"]))
    c, s = _cos(g["theta"]), _sin(g["theta"])
    omega = {
        "r": omega_r, "theta": omega_theta, "z": omega_z,
        "x": _sub(_mul(c, omega_r), _mul(s, omega_theta)),
        "y": _add(_mul(s, omega_r), _mul(c, omega_theta)),
        "z_cartesian": omega_z,
    }
    divergence = _add(_add(directions["u_r_r"],
                           _div(cyl["u_r"], g["r"])), directions["u_z_z"])
    return {
        "cylindrical": omega,
        "divergence": divergence,
        "curl_receipt": ["omega_r=-u_theta_z", "omega_theta=u_r_z-u_z_r",
                          "omega_z=u_theta_r+u_theta/r"],
        "divergence_receipt": "u_r_r+u_r/r+u_z_z",
    }


def _pressure(profile, g, directions):
    base = _add(profile["P0"], profile["p_hist"])
    C2 = _mul(g["C"], g["C"])
    pressure = _mul(C2, base)
    p_R = _mul(C2, profile["p_hist_R"])
    p_Z = _mul(C2, _sub(_add(profile["P0_Z"], profile["p_hist_Z"]),
                        _mul(_div(_mul(_two(), _mul(_add(_one(), g["delta"]),
                                                      g["Z"])), g["d"]), base)))
    radial = _mul(_div(_mul(_two(), g["R"]), g["r"]), p_R)
    axial = _mul(_div(g["d"], _mul(g["kappa"], g["L"])),
                 _sub(p_Z, _mul(_div(_mul(_two(), _mul(g["R"], g["Z"])),
                                    g["d"]), p_R)))
    c, s = _cos(g["theta"]), _sin(g["theta"])
    return {
        "p": pressure, "p_R": p_R, "p_Z": p_Z,
        "gradient": {"X": _mul(c, radial), "Y": _mul(s, radial), "Z": axial},
        "expression": "p=C^2*(P0+p_hist); P0 enters once",
    }


def _volume(cyl, curl, g):
    speed = _add(_add(_mul(cyl["u_r"], cyl["u_r"]),
                      _mul(cyl["u_theta"], cyl["u_theta"])),
                 _mul(cyl["u_z"], cyl["u_z"]))
    omega = curl["cylindrical"]
    omega_sq = _add(_add(_mul(omega["r"], omega["r"]),
                         _mul(omega["theta"], omega["theta"])),
                    _mul(omega["z"], omega["z"]))
    # The accepted exterior payload owns the primary measure.  The fresh
    # r*J2 composition is retained as a nonnumeric comparison witness only.
    primary = g["profile_measure"]
    fresh_rJ2 = _mul(g["saved_r"], g["saved_J2"])
    extra = g["R"]
    return {
        "speed_squared": speed, "vorticity_squared": omega_sq,
        "primary_saved_rJ2": primary,
        "primary_measure_source":
            "accepted_exterior_payload.physical_volume.profile_R_Z_theta_absolute_measure",
        "fresh_saved_rJ2": fresh_rJ2,
        "fresh_rJ2_comparison": "not performed; source identity receipt required",
        "saved_profile_measure": g["profile_measure"],
        "saved_native_x_measure": g["native_x_measure"],
        "native_x_extra_R": extra,
        "native_x_extra_R_receipt": {
            "kind": "native_x_extra_R_receipt",
            "factor": "R",
            "numeric_evaluation": False,
        },
        "speed_squared_primary_integrand": _mul(speed, primary),
        "vorticity_squared_primary_integrand": _mul(omega_sq, primary),
        "speed_squared_native_x_integrand": _mul(_mul(speed, primary), extra),
        "vorticity_squared_native_x_integrand": _mul(_mul(omega_sq, primary), extra),
        "numeric_L2": False, "numeric_Linf": False, "finite_energy": False,
    }


def _physical_compose(body, binding, label):
    inverse = binding["saved_positive_physical_inverse"]
    old = inverse["bindings"]
    return _formula(
        "compose_saved_positive_physical_inverse", (body,),
        label=label,
        physical_arguments=["X_phys", "Y_phys", "z_phys", "t_phys"],
        binding_ref="source_binding_0",
        native_variable_bindings={
            "241787": old["physical_expr:3"],
            "20": old["physical_expr:5"],
            "241536": old["physical_expr:7"],
        },
        free_variable_bindings={
            "tau": old["physical_expr:8"],
            "theta": old["physical_expr:9"],
        },
        saved_inverse_root=inverse["bound_root"],
        saved_inverse_is_data_only=True,
        bindings_apply_inside_all_referenced_native_and_saved_functions=True,
        numerical_inverse_oracle=False,
        circular_definition=False)


def _branch(root, binding, branch_name, g):
    profile = _source_values(root, branch_name, binding)
    cyl = _cylindrical(profile, g)
    parts = _partials(profile, g, cyl)
    directions = _physical_directions(parts, g)
    cart = _cartesian(cyl, directions, g)
    curl = _curl_divergence(cyl, directions, g)
    pressure = _pressure(profile, g, directions)
    volume = _volume(cyl, curl, g)
    physical_functions = {
        "velocity": {name: _physical_compose(body, binding, branch_name + "." + name)
                     for name, body in cart["components"].items()},
        "gradient": {name: {axis: _physical_compose(body, binding, branch_name + ".D_" + axis + "_" + name)
                              for axis, body in axes.items()}
                     for name, axes in cart["gradient"].items()},
        "curl": {axis: _physical_compose(curl["cylindrical"][key], binding, branch_name + ".omega_" + axis)
                 for axis, key in (("X", "x"), ("Y", "y"), ("Z", "z_cartesian"))},
        "pressure": _physical_compose(pressure["p"], binding, branch_name + ".pressure"),
        "pressure_gradient": {axis: _physical_compose(body, binding, branch_name + ".pressure_" + axis)
                              for axis, body in pressure["gradient"].items()},
    }
    return {
        "branch": branch_name,
        "source_only": True,
        "numeric_evaluation": False,
        "profile": profile,
        "cylindrical": cyl,
        "profile_partials": parts,
        "physical_directions": directions,
        "cartesian": cart,
        "curl_and_divergence": curl,
        "pressure": pressure,
        "volume_integrands": volume,
        "physical_Cartesian_functions": physical_functions,
    }


def _pack_scalar_program(original_payload, binding):
    """Pack formula/source objects by object identity, without normalization.

    This is only a representation pass.  It does not inspect operation
    semantics, reorder algebra, differentiate, evaluate, or dispatch on a
    native node.  Nonformula metadata is copied recursively into the root so
    the full named program remains recoverable.
    """
    scalar_kinds = {
        "explicit_formula",
        "exact_rational",
        "native_coordinate_port",
        "native_parameter_port",
        "opaque_parameter_port",
        "free_physical_variable_port",
        "saved_expression_function_port",
        "atomic_source_function_value",
        "atomic_source_derivative_value",
        "opaque_source_atom",
        "closed_parameter_atom",
        "native_x_extra_R",
        "accepted_flux_FToC_j1_formula",
    }
    memo = {}
    nodes = []

    def pack(value):
        if isinstance(value, dict):
            kind = value["kind"] if "kind" in value else None
            if kind in scalar_kinds:
                identity = id(value)
                if identity in memo:
                    return {"$ref": memo[identity]}
                node_id = "scalar_%d" % len(nodes)
                memo[identity] = node_id
                row = {"id": node_id, "kind": kind}
                nodes.append(row)
                for name, child in value.items():
                    if name == "kind":
                        continue
                    if name == "physical_composition_binding":
                        row[name] = {"$binding_ref": "source_binding_0"}
                    else:
                        row[name] = pack(child)
                return {"$ref": node_id}
            return {name: pack(child) for name, child in value.items()}
        if isinstance(value, (list, tuple)):
            return [pack(child) for child in value]
        return value

    roots = pack(original_payload)
    return {
        "representation": "object_identity_scalar_program_v1",
        "numeric_evaluation": False,
        "nodes": nodes,
        "roots": roots,
        "source_bindings": {"source_binding_0": binding},
    }


def derive(capsule, callback=None):
    """Build one conditional formula receipt from the frozen source capsule."""
    _emit(callback, "family_entered", {
        "component": COMPONENT, "branches": ["background", "modified"],
        "source_only": True, "numeric_evaluation": False,
        "family_entries_counted_on_first_emit": list(COUNT_EVENT_MAP["family_entered"]),
    })
    root, binding = _validate_capsule(capsule)
    _emit(callback, "source_join_summary", {
        "ports": [NATIVE_SOURCE_FUNCTION_PORT, SAVED_RADIAL_DERIVATIVE_PORT,
                  SAVED_FLUX_SOURCE_DIFFERENTIAL_PORT],
        "radial_targets": ["E_background", "V_background", "E_N", "V_N"],
        "flux_rows": {"m": [0, 1], "p": [0], "P0": [0]},
        "incoming_m_e_p_preserved": True,
        "unknown_partial_policy": "reject",
    })
    closed = _closed_parameterization(root)
    _emit(callback, "closed_parameterization_summary", {
        "raw_record_count": 46, "fixed_named_atoms": [8, 10, 35560],
        "opaque_N_node": 18, "core_logC_value_selected": False,
        "conditional_only": True,
    })
    g = _geometry(root, binding)
    _emit(callback, "cylindrical_formula_summary", {
        "definitions": [
            "d=1-Z^2", "L=1-delta*Z^2", "a=1-delta",
            "lambda=sqrt(tau/d)", "r=lambda*sqrt(2R)", "kappa=lambda^a",
            "C=Pstar*lambda^(-1-delta)", "B=Pstar*R/(rL)",
            "S=d*m_Z+a*Z*m-2*Z*V", "u_r=-B*S", "u_theta=C*E", "u_z=C*V",
        ], "branches": ["background", "modified"],
    })
    _emit(callback, "physical_direction_summary", {
        "D_r": "(2R/r)*D_R",
        "D_z": "d/(kappa*L)*(D_Z-(2RZ/d)*D_R)",
        "fixed_physical_time": True, "ordinary_Z_source_jets_retained": True,
    })
    branches = {
        name: _branch(root, binding, name, g)
        for name in ("background", "modified")
    }
    _emit(callback, "cartesian_jet_summary", {
        "entries": 9,
        "basis_derivatives": {"D_theta_u": "-v", "D_theta_v": "u", "D_theta_w": "0"},
        "native_bindings_inside_all_source_functions": {"R": 241787, "Z": 20,
                                                          "native_x": 241536,
                                                          "tau": True, "theta": True},
    })
    _emit(callback, "curl_divergence_summary", {
        "curl": ["omega_r=-u_theta_z", "omega_theta=u_r_z-u_z_r",
                  "omega_z=u_theta_r+u_theta/r"],
        "divergence": "u_r_r+u_r/r+u_z_z", "numeric_zero_claim": False,
    })
    _emit(callback, "pressure_volume_summary", {
        "pressure": "p=C^2*(P0+p_hist)", "P0_entries": 1,
        "primary_volume_factor": "saved r*J2",
        "profile_measure": "physical_expr:201", "native_x_extra": "R",
        "numeric_norms": False,
    })
    full_payload = {
        "component": COMPONENT, "source_only": True, "numeric_evaluation": False,
        "conditional_scope": "regular interior O3 source cells at fixed physical time",
        "source_jet_admission": {
            "actual_source_Z_jets_admitted": False,
            "conditional_source_Z1_Z2_contract_consumed": True,
            "E_V_field_rows_normalized": True,
            "actual_Z_normalization_admission": False,
        },
        "closed_parameterization": closed, "geometry": g, "branches": branches,
        "binding_certificate": binding,
        "accepted_inputs_preserved": {
            "source_pins": root["source_pins"],
            "preservation": root["preservation"],
            "physical_source_capsule": root["physical_source_capsule"],
            "accepted_radial_result": root["accepted_radial_result"],
            "accepted_flux_source_rows": root["accepted_flux_source_rows"],
            "accepted_exterior_payload": root["accepted_exterior_payload"],
        },
        "proof_rows": [
            {"row": "source_lineage", "equation": "atomic source ports; no calls"},
            {"row": "closed_parameter_join", "equation": "finite 46-record literal pattern"},
            {"row": "cylindrical_streamfunction", "equation": "u_r=-B*S,u_theta=C*E,u_z=C*V"},
            {"row": "profile_partials", "equation": "six explicit R/Z product derivatives"},
            {"row": "physical_direction_chain", "equation": "D_r,D_z explicit"},
            {"row": "cartesian_basis", "equation": "D_X,D_Y with D_theta basis derivatives"},
            {"row": "curl_divergence", "equation": "cylindrical curl and divergence"},
            {"row": "pressure_gradient", "equation": "p=C^2*(P0+p_hist)"},
            {"row": "volume_integrands", "equation": "speed^2,omega^2 times saved r*J2"},
        ],
        "all_old_caps_and_errors_preserved": True,
        "no_bounds_midpoints_or_finite_differences": True,
        "no_phase_seam_axis_or_global_NS_admission": True,
    }
    scalar_program = _pack_scalar_program(full_payload, binding)
    packed_roots = scalar_program["roots"]
    branch_output_refs = {
        branch_name: {
            section: packed_roots["branches"][branch_name][section]
            for section in (
                "profile", "cylindrical", "profile_partials",
                "physical_directions", "cartesian", "curl_and_divergence",
                "pressure", "volume_integrands", "physical_Cartesian_functions",
            )
        }
        for branch_name in ("background", "modified")
    }
    payload = {
        "component": COMPONENT,
        "source_only": True,
        "numeric_evaluation": False,
        "conditional_scope": full_payload["conditional_scope"],
        "source_jet_admission": full_payload["source_jet_admission"],
        "scalar_program": scalar_program,
        "branch_output_refs": branch_output_refs,
        "geometry_refs": packed_roots["geometry"],
        "shape_summary": {
            "branches": ["background", "modified"],
            "cartesian_gradient_entries_per_branch": 9,
            "cylindrical_velocity_components_per_branch": 3,
            "physical_direction_operators": 2,
            "proof_rows": 9,
        },
        "binding_certificate": binding,
        "proof_rows": full_payload["proof_rows"],
        "accepted_inputs_preserved": {
            "source_pins": root["source_pins"],
            "preservation": root["preservation"],
        },
        "all_old_caps_and_errors_preserved": True,
        "no_bounds_midpoints_or_finite_differences": True,
        "no_phase_seam_axis_or_global_NS_admission": True,
    }
    result = {
        "schema": SCHEMA, "component": COMPONENT, "passed": True,
        "conditional": True, "metadata_only": False, "source_only": True,
        "numeric_evaluation": False, "full_goal_complete": False,
        "actual_counts": dict(SUCCESS_COUNTS), "new_flags": dict(SUCCESS_FLAGS),
        "proof_count": 1, "callback_event_counts": dict(CALLBACK_EVENT_COUNTS),
        "count_event_map": dict(COUNT_EVENT_MAP), "payload": payload,
    }
    _emit(callback, "family_completed", {
        "component": COMPONENT, "proof_count": 1, "formula_proof_rows": 9,
        "conditional": True, "full_goal_complete": False,
        "numeric_evaluation": False,
    })
    return result
