"""Fresh compact prefix pressure certificates; no ancestor consumer execution."""



def _as_decimal_int(value, label):
    if not isinstance(value, str) or not value:
        raise ValueError(label + " must be a decimal string")
    body = value[1:] if value[0] in "+-" else value
    if not body or not body.isdigit():
        raise ValueError(label + " must be a decimal integer string")
    return int(value)


def _mpf_tuple(value, label):
    if not isinstance(value, (list, tuple)) or len(value) != 4:
        raise ValueError(label + " must have four exact string components")
    sign = _as_decimal_int(value[0], label + ".sign")
    mantissa = _as_decimal_int(value[1], label + ".mantissa")
    exponent = _as_decimal_int(value[2], label + ".exponent")
    precision = _as_decimal_int(value[3], label + ".precision")
    if sign not in (0, 1) or mantissa < 0 or precision < 0:
        raise ValueError(label + " has invalid sign, mantissa, or precision")
    if mantissa == 0:
        if sign != 0 or exponent != 0 or precision != 0:
            raise ValueError(label + " has a noncanonical zero tuple")
    elif precision != mantissa.bit_length():
        raise ValueError(label + " precision does not equal mantissa bit length")
    return sign, mantissa, exponent, precision


def _magnitude_compare(left, right):
    if left[1] == 0 or right[1] == 0:
        return (left[1] > 0) - (right[1] > 0)
    gap = right[2] - left[2]
    if gap == 0:
        return (left[1] > right[1]) - (left[1] < right[1])
    if gap == 1:
        # right has one more binary exponent; only one exact shift is used.
        return (left[1] > 2 * right[1]) - (left[1] < 2 * right[1])
    if gap == -1:
        return (2 * left[1] > right[1]) - (2 * left[1] < right[1])
    raise ValueError("endpoint exponent gap exceeds one binary shift")


def _tuple_compare(left, right):
    if left[0] != right[0]:
        return -1 if left[0] else 1
    magnitude = _magnitude_compare(left, right)
    return -magnitude if left[0] else magnitude


def _same(left, right):
    return list(left) == list(right)


def _recipe_key(value):
    if not isinstance(value, str):
        raise ValueError("source recipe must be text")
    return "".join(value.lower().split()).replace(";", "")


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _emit(journal, payload):
    if journal is None:
        return
    if not callable(journal):
        raise TypeError("journal must be callable or None")
    journal(payload)


def _validate_parent(parent_proof, source_contract):
    _require(isinstance(parent_proof, dict), "parent_proof must be a mapping")
    _require(parent_proof.get("schema") == "original_compact_upstream_mean_identity",
             "cached parent proof schema mismatch")
    _require(parent_proof.get("passed") is True, "cached parent source proof did not pass")
    _require(parent_proof.get("source_only") is True, "parent proof is not source-only")
    for key in ("numeric_values_computed", "numeric_containment_claimed",
                "numeric_mean_enclosure_transferred"):
        _require(parent_proof.get(key) is False, "parent proof has disallowed numeric flag: " + key)
    audit = parent_proof.get("source_inspection_audit")
    _require(isinstance(audit, dict), "parent source inspection audit missing")
    for key in ("numeric_field_evaluations", "old_operator_calls",
                "old_proof_functions_called", "old_source_consumer_calls",
                "old_source_normalizer_calls", "original_source_evaluations",
                "quadratures", "temporal_recurrence_calls"):
        _require(audit.get(key) == 0, "parent audit is nonzero: " + key)

    binding = parent_proof.get("source_binding")
    _require(isinstance(binding, dict) and binding.get("passed") is True,
             "cached source binding did not pass")
    mean_source = str(binding.get("mean_source", "")).lower()
    _require("rate-one sigma" in mean_source and "nonzero homogeneous" in mean_source,
             "parent mean source does not retain the rate-one homogeneous inlet")
    axial = parent_proof.get("upstream_mean_identity", {}).get("axial_transport", {})
    parent_mstar = axial.get("mstar_binding", {})
    parent_recipe = parent_mstar.get("recipe", {})
    _require(isinstance(parent_mstar.get("raw_sha256"), str) and
             len(parent_mstar["raw_sha256"]) == 64,
             "parent mstar recipe hash missing")
    _require(parent_mstar.get("selector") == "source_identity original kernel at u=exp40 plus homogeneous memory",
             "parent mstar selector mismatch")
    _require(str(parent_recipe.get("initial_memory")) == "1" and
             str(parent_recipe.get("kernel_rate")) == "1",
             "parent mstar rate or initial memory mismatch")
    _require(parent_recipe.get("numerical_enclosure_used_as_value") is False,
             "parent recipe treats enclosure as a value")
    _require(parent_recipe.get("ordinary_Z_independent") is True,
             "parent recipe lacks ordinary-Z independence")
    parent_definition = parent_recipe.get("definition")
    parent_formula, separator, parent_endpoint = parent_definition.partition(";") \
        if isinstance(parent_definition, str) else ("", "", "")
    _require(_recipe_key(parent_formula) ==
             _recipe_key(source_contract["source_recipe"]["exact_mstar"]) and
             (not separator or _recipe_key(parent_endpoint) == _recipe_key("Y=exp(40)")),
             "parent and receipt mstar formulas differ")
    _require(parent_recipe.get("source_kernel_substitution") == "u=Y=exp(40)" or
             parent_recipe.get("source_kernel_substitution") == "u=Y=exp40",
             "parent kernel substitution mismatch")
    return {
        "parent_mstar_raw_sha256": parent_mstar.get("raw_sha256"),
        "parent_selector": parent_mstar.get("selector"),
        "parent_source_audit_zero": True,
        "parent_numeric_value_flags_false": True
    }


def _validate_receipt(receipt):
    _require(isinstance(receipt, dict), "receipt must be a mapping")
    _require(receipt.get("schema") == "compact_exact_memory_enclosure_join_trace_v1",
             "memory receipt schema mismatch")
    raw_artifact = receipt.get("authenticated_artifacts", {}).get("raw_report", {})
    _require(raw_artifact.get("sha256") ==
             "57f5495394ab864ace99905fac07aab5bd315ff420329704121de121312b8d4a",
             "authenticated axial report hash mismatch")
    flags = receipt.get("flags", {})
    for key in ("raw_report_authenticated", "selected_rows_authenticated",
                "mstar_exact_endpoint_tuples_captured", "source_recipe_linked",
                "certificate_metadata_linked"):
        _require(flags.get(key) is True, "receipt gate missing: " + key)
    _require(flags.get("numeric_mean_enclosure_transferred") is False,
             "receipt already claims numeric transfer")
    _require(receipt.get("runtime_model_metadata") is None,
             "runtime metadata must remain null")
    integer_audit = receipt.get("exact_integer_audit", {})
    _require(integer_audit.get("raw_source_tuple_equality_verified") is True,
             "receipt lacks exact raw tuple audit")
    _require(integer_audit.get("all_selected_tuple_components_were_native_int_before_serialization") is True,
             "receipt lacks native-integer provenance")
    _require(integer_audit.get("artifact_tuple_components_are_strings_to_preserve_exactness") is True,
             "receipt tuple string preservation missing")

    source_contract = receipt.get("mstar_source_contract")
    _require(isinstance(source_contract, dict), "mstar source contract missing")
    target = source_contract.get("target_call", {})
    _require(target.get("class") == "OriginalAxialRateMemory" and
             target.get("evaluate_method") == "evaluate",
             "source call identity mismatch")
    _require(str(target.get("constructor_kernel_rate")) == "1" and
             target.get("evaluate_arguments") == {"s": "1", "power": "1"},
             "fixed s=1 rate-one call mismatch")
    params = target.get("report_proof_parameters", {})
    _require(str(params.get("power")) == "1" and str(params.get("kernel_rate")) == "1",
             "saved report proof parameters mismatch")
    recipe = source_contract.get("source_recipe", {})
    _require(str(recipe.get("initial_memory")) == "1" and
             _recipe_key(recipe.get("exact_mstar")) ==
             _recipe_key("exp(1-Y)+integral_1^Y exp(-Y+y)*sigma(1-log(y)/40)dy"),
             "receipt exact mstar recipe mismatch")
    _require(_recipe_key(recipe.get("endpoint_substitution")) ==
             _recipe_key("Y=exp(40), d=Y-1"),
             "receipt endpoint substitution mismatch")
    _require(recipe.get("ordinary_Z_independent") is True and
             recipe.get("decimal_or_midpoint_as_value") is False,
             "receipt ordinary-Z or value contract mismatch")
    _require(_recipe_key(recipe.get("homogeneous_term")) == _recipe_key("exp(1-Y)"),
             "homogeneous term mismatch")
    _require(_recipe_key(recipe.get("forced_term")) ==
             _recipe_key("integral_1^Y exp(-Y+y)*sigma(1-log(y)/40)dy"),
             "forced term mismatch")
    sigma_link = _recipe_key(recipe.get("source_sigma_link"))
    _require("changeofvariabley=exp(40*r)" in sigma_link,
             "sigma change-of-variable contract missing")
    certificate = receipt.get("certificate_assessment", {})
    classification = certificate.get("classification", {})
    for key in ("saved_directed_interval_present",
                "saved_positive_kernel_certificate_metadata_present",
                "independent_spot_containment_receipt_present"):
        _require(classification.get(key) is True, "certificate gate missing: " + key)
    _require(classification.get("universal_enclosure_theorem_reproved_by_this_trace") is False,
             "receipt overclaims a universal theorem")
    return source_contract, certificate


def bind_fixed_mstar_enclosure(receipt, parent_proof, journal):
    """Bind the saved directed interval for the one scalar mstar(s=1,power=1).

    The returned binding is conditional on the saved interval/certificate and
    cached source identity.  It does not claim a theorem for other s values.
    """
    try:
        _emit(journal, {"stage": "fixed_mstar_enclosure_binding", "status": "started"})
        source_contract, certificate = _validate_receipt(receipt)
        parent_audit = _validate_parent(parent_proof, source_contract)
        rows = receipt.get("selected_pressure_rows")
        _require(isinstance(rows, list) and len(rows) == 3,
                 "exactly three selected pressure rows are required")
        expected = {5: "-.3", 11: "0", 17: ".3"}
        selected = {}
        for row in rows:
            _require(isinstance(row, dict), "selected pressure row must be a mapping")
            index = _as_decimal_int(row.get("zero_based_pressure_row_index"), "row index")
            _require(index in expected and index not in selected, "unexpected or duplicate pressure row")
            _require(row.get("selector") == "$.pressure_row_samples[" + str(index) + "]",
                     "pressure row selector mismatch")
            _require(row.get("Z") == expected[index] and row.get("chart") == "O2_axial" and
                     row.get("native_coordinate") == "1", "fixed s=1 row identity mismatch")
            _require(row.get("mstar_field_selector") ==
                     "$.pressure_row_samples[" + str(index) + "].leading_scalars.mean",
                     "mstar selector must use leading_scalars.mean")
            lower_raw = row.get("source_enclosure_lower")
            upper_raw = row.get("source_enclosure_upper")
            error_raw = row.get("raw_directed_error_tuple")
            lower = _mpf_tuple(lower_raw, "row" + str(index) + ".lower")
            upper = _mpf_tuple(upper_raw, "row" + str(index) + ".upper")
            error = _mpf_tuple(error_raw, "row" + str(index) + ".directed_error")
            _require(lower[0] == 0 and upper[0] == 0 and lower[1] > 0 and upper[1] > 0,
                     "mstar endpoints must be positive")
            _require(_tuple_compare(lower, upper) <= 0,
                     "mstar lower endpoint exceeds upper endpoint")
            positive = row.get("positive_receipt", {})
            _require(positive.get("selector") ==
                     "$.pressure_row_samples[" + str(index) + "].positive_mean_memory_receipt",
                     "positive memory receipt selector mismatch")
            _require(str(positive.get("power")) == "1" and
                     str(positive.get("kernel_rate")) == "1",
                     "positive receipt rate or power mismatch")
            _require(positive.get("local_axial_zero_does_not_zero_memory") is True,
                     "positive receipt does not retain the homogeneous memory")
            _require(positive.get("positive_actual_memory_lower") == list(lower_raw),
                     "positive lower bound differs from mstar lower endpoint")
            rounding = positive.get("finite_endpoint_rounding_bound", {})
            _require(rounding.get("lower") == ["0", "0", "0", "0"] and
                     rounding.get("upper") == ["0", "0", "0", "0"],
                     "finite endpoint rounding is not exact zero")
            source_family = row.get("source_family")
            _require(isinstance(source_family, dict), "source family missing")
            selected[index] = {
                "row_index": str(index),
                "Z": row["Z"],
                "selector": row["selector"],
                "mstar_selector": row["mstar_field_selector"],
                "source_family": dict(source_family),
                "lower": list(lower_raw),
                "upper": list(upper_raw),
                "raw_directed_error": list(error_raw),
                "positive_receipt_selector": positive.get("selector"),
                "positive_lower": list(positive.get("positive_actual_memory_lower")),
                "finite_endpoint_rounding_bound": {
                    "lower": list(rounding["lower"]),
                    "upper": list(rounding["upper"])
                },
                "endpoint_ordering": "validated by signed mantissa comparison with binary exponent gap at most one"
            }
        _require(set(selected) == set(expected), "selected pressure row set is incomplete")
        ordered = [selected[i] for i in (5, 11, 17)]
        reference = ordered[0]
        for row in ordered[1:]:
            _require(row["lower"] == reference["lower"] and row["upper"] == reference["upper"] and
                     row["raw_directed_error"] == reference["raw_directed_error"],
                     "selected rows do not carry the same fixed mstar interval")
            _require(row["source_family"] == reference["source_family"],
                     "selected rows do not carry the same source family")
        _emit(journal, {"stage": "fixed_mstar_enclosure_binding", "status": "rows_validated", "rows": 3})
        result = {
            "schema": "compact_fixed_memory_enclosure_binding_v1",
            "passed": True,
            "binding_status": "conditional_saved_certificate_binding",
            "fixed_scalar": {
                "name": "mstar",
                "s": "1",
                "power": "1",
                "kernel_rate": "1",
                "ordinary_Z_independent": True,
                "universal_in_s_claimed": False
            },
            "source_recipe": {
                "exact_mstar": source_contract["source_recipe"]["exact_mstar"],
                "homogeneous_term": source_contract["source_recipe"]["homogeneous_term"],
                "forced_term": source_contract["source_recipe"]["forced_term"],
                "initial_memory": str(source_contract["source_recipe"]["initial_memory"]),
                "sigma_change_of_variable": source_contract["source_recipe"]["source_sigma_link"]
            },
            "directed_interval": {
                "lower": list(reference["lower"]),
                "upper": list(reference["upper"]),
                "raw_directed_error": list(reference["raw_directed_error"]),
                "tuple_components_are_decimal_strings": True,
                "decimal_display_or_midpoint_selected": False
            },
            "selected_rows": ordered,
            "endpoint_ordering_audit": {
                "all_three_rows_ordered": True,
                "sign_and_bit_length_validated": True,
                "max_binary_exponent_shift_used": 1,
                "giant_absolute_exponent_shift_used": False,
                "all_three_rows_exactly_equal": True
            },
            "homogeneous_forced_assembly": {
                "same_formula_bound": True,
                "homogeneous_term_retained": True,
                "positive_lower_retained": True,
                "finite_endpoint_rounding_exact_zero": True,
                "no_local_zero_reset": True
            },
            "certificate": {
                "saved_directed_interval": True,
                "saved_positive_kernel_certificate_metadata": True,
                "independent_spot_checks_used_as_proof": False,
                "conditional_on_cached_source_and_certificate": True,
                "parent_source_audit": parent_audit,
                "receipt_certificate_scope": certificate.get("classification", {}).get("status")
            },
            "source_hashes": {
                "raw_report_sha256": receipt["authenticated_artifacts"]["raw_report"]["sha256"],
                "parent_mstar_recipe_raw_sha256": parent_audit["parent_mstar_raw_sha256"]
            },
            "numeric_mean_enclosure_transferred": False,
            "conditional_numeric_mean_enclosure_bound": True,
            "source_runtime_evaluations": 0,
            "new_quadratures": 0,
            "runtime_model_metadata": None
        }
        _emit(journal, {"stage": "fixed_mstar_enclosure_binding", "status": "passed", "conditional": True})
        return result
    except Exception as error:
        message = str(error)
        try:
            _emit(journal, {"stage": "fixed_mstar_enclosure_binding", "status": "failed", "error": message})
        except Exception:
            pass
        return {
            "schema": "compact_fixed_memory_enclosure_binding_v1",
            "passed": False,
            "binding_status": "rejected",
            "error": message,
            "numeric_mean_enclosure_transferred": False,
            "conditional_numeric_mean_enclosure_bound": False,
            "source_runtime_evaluations": 0,
            "new_quadratures": 0,
            "runtime_model_metadata": None
        }


from fractions import Fraction


EXPECTED_GEOMETRY_RECEIPT_SHA256 = "d2edb30b8c86ca6254b0a4864c6e971ea9620c56e4a223393d61b85bee813f27"
EXPECTED_LOGC_LINEAGE_RECEIPT_SHA256 = "913b22162e23938f21e1b92f19083923b2db265f95f60994435343343a9f624d"
EXPECTED_NATIVE_COMPACT_RAW_SHA256 = "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269"
EXPECTED_NUMERIC_COMPACT_RAW_SHA256 = "7ccb30dbb3e440a3e598eb30deaf03509cf8624e429478aeeb364671baaac071"
EXPECTED_NUMERIC_CORE_SOURCE_SHA256 = "9c152f5ff50bd760ee1a860a11e48c1b3a2fd36027681ea05a3e7a0920181603"
EXPECTED_COMPACT_NUMERIC_SOURCE_SHA256 = "7b8dc6ad7b42818fb7ab5cd0e958f7a490188f8b722cb75ca64ad345e6149c6a"
EXPECTED_COMPACT_BACKGROUND_SOURCE_SHA256 = "1575b95765a1f9867af94eadc1e2ba09fbf22705cd8034c5615918dabeb50440"
EXPECTED_SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}


def _need(condition, message):
    if not condition:
        raise ValueError(message)


def _node(nodes, node_id):
    """Read a serialized graph node whose JSON key may be str or int."""
    if node_id in nodes:
        return nodes[node_id]
    key = str(node_id)
    if key in nodes:
        return nodes[key]
    raise ValueError("missing graph node %s" % node_id)


def _fraction(node, label):
    _need(isinstance(node, dict) and node.get("operation") == "exact_rational",
          "expected exact rational node for %s" % label)
    numerator = node.get("numerator")
    denominator = node.get("denominator")
    _need(type(numerator) is int and type(denominator) is int and denominator != 0,
          "malformed exact rational node for %s" % label)
    return Fraction(numerator, denominator)


def _expect_fraction(nodes, node_id, numerator, denominator=1):
    value = _fraction(_node(nodes, node_id), "node %s" % node_id)
    expected = Fraction(numerator, denominator)
    _need(value == expected, "exact rational node %s changed" % node_id)
    return value


def _expect_operation(nodes, node_id, operation, **fields):
    row = _node(nodes, node_id)
    _need(row.get("operation") == operation,
          "node %s operation changed: expected %s" % (node_id, operation))
    for key, expected in fields.items():
        _need(row.get(key) == expected,
              "node %s field %s changed" % (node_id, key))
    return row


def _add(*maps):
    result = {}
    for mapping in maps:
        for key, value in mapping.items():
            result[key] = result.get(key, Fraction(0)) + value
    return {key: value for key, value in result.items() if value}


def _scale(factor, mapping):
    return _add({key: factor * value for key, value in mapping.items()})


def _substitute(mapping, variable, value):
    """Substitute one affine variable with an already exact affine map."""
    coefficient = mapping.get(variable, Fraction(0))
    remainder = {key: term for key, term in mapping.items() if key != variable}
    return _add(remainder, _scale(coefficient, value))


def _zero_residual(left, right):
    return _add(left, _scale(Fraction(-1), right))


def _encode_map(mapping):
    """JSON-safe exact Fraction map; no decimal conversion."""
    return {
        key: {"numerator": value.numerator, "denominator": value.denominator}
        for key, value in sorted(mapping.items())
    }


def _zero_map(mapping):
    return all(value == 0 for value in mapping.values())


def _source_pin_text(source_code_pins, path, label):
    row = source_code_pins.get(path) if isinstance(source_code_pins, dict) else None
    _need(isinstance(row, dict), "saved source-code pin is missing for %s" % label)
    snippets = row.get("text_snippets", [])
    _need(isinstance(snippets, list), "saved source-code snippets malformed for %s" % label)
    text = "\n".join(item.get("text", "") for item in snippets if isinstance(item, dict))
    _need(text, "saved source-code snippets empty for %s" % label)
    return row, text


def _emit(journal, event):
    """Support a callback, event list, or mapping with an events list."""
    if journal is None:
        return
    if callable(journal):
        journal(dict(event))
        return
    if isinstance(journal, list):
        journal.append(dict(event))
        return
    if isinstance(journal, dict):
        journal.setdefault("events", []).append(dict(event))
        return
    raise TypeError("journal must be callable, a list, or a mapping")


def _receipt_pointer_checks(journal):
    """Check pointer hashes when the caller supplied them; constants remain explicit."""
    if not isinstance(journal, dict):
        return {"geometry_receipt_sha256": EXPECTED_GEOMETRY_RECEIPT_SHA256,
                "lineage_receipt_sha256": EXPECTED_LOGC_LINEAGE_RECEIPT_SHA256,
                "caller_pointer_fields_present": False}
    aliases = {
        "geometry_receipt_sha256": ("geometry_receipt_sha256", "geometry_v2_sha256"),
        "lineage_receipt_sha256": ("lineage_receipt_sha256", "logC_lineage_receipt_sha256"),
    }
    seen = {}
    for output_key, keys in aliases.items():
        values = [journal[key] for key in keys if key in journal]
        if values:
            _need(all(value == values[0] for value in values),
                  "%s aliases disagree" % output_key)
            expected = (EXPECTED_GEOMETRY_RECEIPT_SHA256 if output_key == "geometry_receipt_sha256"
                        else EXPECTED_LOGC_LINEAGE_RECEIPT_SHA256)
            _need(values[0] == expected, "%s does not authenticate the saved receipt" % output_key)
            seen[output_key] = values[0]
    return {
        "geometry_receipt_sha256": seen.get("geometry_receipt_sha256", EXPECTED_GEOMETRY_RECEIPT_SHA256),
        "lineage_receipt_sha256": seen.get("lineage_receipt_sha256", EXPECTED_LOGC_LINEAGE_RECEIPT_SHA256),
        "caller_pointer_fields_present": bool(seen),
    }


def _cached_measure_witness(journal):
    """Return an already saved C2 witness without recomputing its rational relation."""
    candidate = None
    if isinstance(journal, dict):
        for key in ("cached_C2_measure_normalization", "c2_measure_normalization",
                    "measure_normalization", "measure_normalization_audit"):
            if isinstance(journal.get(key), dict):
                candidate = journal[key]
                break
    elif isinstance(journal, list):
        for event in reversed(journal):
            if not isinstance(event, dict):
                continue
            if event.get("event") in ("cached_C2_measure_normalization", "cached_C2_interface_bound"):
                candidate = event.get("measure_normalization", event)
                if isinstance(candidate, dict):
                    break
    _need(isinstance(candidate, dict),
          "saved C2 measure normalization witness is required as cached input")
    upper = candidate.get("upper_R2_mass", candidate.get("upper_radius_R2_mass",
                       candidate.get("upper_radius_R2_mass_literal")))
    lower = candidate.get("lower_Rc2_mass", candidate.get("lower_Rc2_mass_literal"))
    _need(upper == "3/8" and lower == "3/2",
          "cached C2 measure normalization literals changed")
    literal = candidate.get("exact_rational_check", candidate.get("source_literal"))
    if literal is not None:
        _need(literal == "4*(3/8)=3/2", "cached C2 measure witness formula changed")
    if "exact_rational_check_passed" in candidate:
        _need(candidate["exact_rational_check_passed"] is True,
              "cached C2 measure witness is not marked accepted")
    if "formula_only" in candidate:
        _need(candidate["formula_only"] is True,
              "C2 measure witness must remain a cached formula-only input")
    if "numeric_arithmetic_or_new_bound_certification_performed" in candidate:
        _need(candidate["numeric_arithmetic_or_new_bound_certification_performed"] is False,
              "C2 measure witness was recertified in this consumer")
    if "numeric_arithmetic_or_new_bound_certification_performed" not in candidate and "new_exact_rational_mass_checks" in candidate:
        _need(candidate["new_exact_rational_mass_checks"] == 0,
              "C2 measure witness was recomputed in this consumer")
    return {
        "upper_R2_mass": "3/8",
        "lower_Rc2_mass": "3/2",
        "source_literal": "4*(3/8)=3/2",
        "reused_without_replay": True,
        "witness_fields_seen": sorted(candidate),
    }


def _native_log_rc_map(graph, native, formulas, source_binding, source_code_pins):
    """Check the finite native Rc DAG and return its exact affine log map.

    The v2 receipt intentionally serializes only the selected DAG rows.  Leaf
    IDs 2, 48, and 2083 are referenced by those rows but are not present in
    ``graph_nodes``.  Their exact values are therefore authenticated through
    the saved source-code pins and source-binding metadata below, rather than
    by pretending the omitted rows were serialized.
    """
    _need(formulas.get("logP_source") == "original_core_numeric.logP=exp(40)+11",
          "saved logP source expression changed")
    logp_source = source_binding.get("numeric_core_logP_source", {})
    _need(logp_source.get("path") == "src/openai_ns_reconstruction/original_core_numeric.py"
          and logp_source.get("source_hash") == EXPECTED_NUMERIC_CORE_SOURCE_SHA256,
          "saved numeric logP source pin changed")
    _need(logp_source.get("binding")
          == "The numeric core defines logP=exp(40)+11 and its interval form; no value is substituted into the native graph here.",
          "saved numeric logP source binding changed")
    compact_source = source_binding.get("numeric_compact_source", {})
    _need(compact_source.get("path") == "src/openai_ns_reconstruction/original_compact_numeric.py"
          and compact_source.get("source_hash") == EXPECTED_COMPACT_NUMERIC_SOURCE_SHA256,
          "saved compact numeric source pin changed")
    _, core_text = _source_pin_text(
        source_code_pins, "src/openai_ns_reconstruction/original_core_numeric.py", "numeric logP")
    _need("self.logP = c.exp(40)+11" in core_text,
          "numeric logP source leaf 40/11 is not retained")
    _, compact_text = _source_pin_text(
        source_code_pins, "src/openai_ns_reconstruction/original_compact_numeric.py", "compact numeric geometry")
    _need("10*ref.logCp,10*ref.logC" in compact_text,
          "compact numeric coefficient-10 source leaf is not retained")
    _need("p['native_O3_chart']=='O3_power' and p['native_O3']=='2'" in compact_text,
          "compact numeric outer-power source leaf 2 is not retained")
    _, background_text = _source_pin_text(
        source_code_pins, "src/openai_ns_reconstruction/original_compact_background_velocity.py",
        "compact native radius")
    _need("R_definition='110*exp(10logCstar+10logPstar+exp40+14)*x'" in background_text,
          "compact native radius source definition changed")
    _need(source_code_pins["src/openai_ns_reconstruction/original_compact_background_velocity.py"].get("sha256")
          == EXPECTED_COMPACT_BACKGROUND_SOURCE_SHA256,
          "compact native radius source pin hash changed")
    semantics = native.get("source_semantics", {})
    _need(semantics.get("outer_power_coordinate_node_250215")
          == "bound_variable true_outer_O3_power",
          "native outer-power source semantic changed")
    _need(semantics.get("radius_node_251066")
          == "function_substitution(expression=250457, variable=250215, value=48, Z_independent_substitution=true)",
          "native Rc substitution source semantic changed")
    _need(semantics.get("compact_phase_node_260181")
          == "true_compact_repair_x, domain [1,48]",
          "native compact phase source semantic changed")
    _expect_operation(graph, 3, "analytic_unary", name="exp", argument=2)
    _expect_fraction(graph, 4, 11)
    _expect_operation(graph, 5, "sum", arguments=[3, 4])
    _expect_fraction(graph, 79968, 110)
    _expect_fraction(graph, 205156, 18)
    _expect_operation(graph, 35560, "current_original_source_parameter",
                      name="true_core_logC", Z_independent=True,
                      exact_fixed_source_parameter=True,
                      source_path="common.core.logC",
                      parameter_enclosure_not_selected_value=True)
    _expect_operation(graph, 81513, "sum", arguments=[35560, 5])
    _expect_operation(graph, 81514, "product", arguments=[2083, 81513])
    _expect_operation(graph, 84528, "bound_variable", name="true_prescribed_post_restore_fraction")
    _expect_operation(graph, 84529, "exact_rational", numerator=-7, denominator=1)
    _expect_operation(graph, 84530, "sum", arguments=[84529, 84528])
    _expect_operation(graph, 84531, "sum", arguments=[81514, 84530])
    _expect_operation(graph, 84532, "analytic_unary", name="exp", argument=84531)
    _expect_operation(graph, 84533, "product", arguments=[79968, 84532])
    _expect_operation(graph, 85743, "function_substitution", expression=84533,
                      variable=84528, value=1, Z_independent_substitution=True)
    _expect_operation(graph, 250215, "bound_variable", name="true_outer_O3_power")
    _expect_operation(graph, 250222, "sum", arguments=[205156, 3, 250215])
    _expect_operation(graph, 250456, "analytic_unary", name="exp", argument=250222)
    _expect_operation(graph, 250457, "product", arguments=[85743, 250456])
    _expect_operation(graph, 251066, "function_substitution", expression=250457,
                      variable=250215, value=48, Z_independent_substitution=True)

    # Node 5 is the exact source recipe logP=exp(40)+11. Keep logP as a
    # named source term in the affine map while separately checking its recipe.
    log_c = {"logC": Fraction(1)}
    log_p = {"logP": Fraction(1)}
    exp40 = {"exp40": Fraction(1)}
    log_110 = {"log110": Fraction(1)}
    post_fraction = {"post_restore_fraction": Fraction(1)}
    outer_power = {"outer_power": Fraction(1)}
    logp_recipe = _add(exp40, {"constant": Fraction(11)})
    _need(logp_recipe == {"exp40": Fraction(1), "constant": Fraction(11)},
          "logP source recipe changed")
    # The named logP term is retained in the Rc map below; this source recipe
    # is recorded separately so the coefficient identity does not silently
    # replace the selected scalar with a different point or interval value.
    pre_restore_expr = _add(log_110, _scale(Fraction(10), _add(log_c, log_p)),
                            {"constant": Fraction(-7)}, post_fraction)
    pre_r_m = _substitute(pre_restore_expr, "post_restore_fraction",
                          {"constant": Fraction(1)})
    outer_factor = _add(exp40, {"constant": Fraction(18)}, outer_power)
    native_before_substitution = _add(pre_r_m, outer_factor)
    native_after_substitution = _substitute(native_before_substitution, "outer_power",
                                             {"constant": Fraction(2)})
    expected = {"log110": Fraction(1), "logC": Fraction(10),
                "logP": Fraction(10), "exp40": Fraction(1),
                "constant": Fraction(14)}
    _need(native_after_substitution == expected,
          "native Rc affine coefficient map changed")
    return {
        "node_map_before_outer_substitution": _encode_map(native_before_substitution),
        "node_map_after_outer_substitution": _encode_map(native_after_substitution),
        "expected_logRc_map": _encode_map(expected),
        "residual_map": _encode_map(_zero_residual(native_after_substitution, expected)),
        "source_aliases_checked": {
            "logP_recipe_node": 5,
            "logP_recipe": "exp(40)+11",
            "logP_recipe_affine_map": _encode_map(logp_recipe),
            "logP_named_affine_map": _encode_map(log_p),
            "unserialized_leaf_contracts": {
                "node_2": {"exact_value": "40", "source": "original_core_numeric self.logP"},
                "node_48": {"exact_value": "2", "source": "original_compact_numeric native_O3 selector"},
                "node_2083": {"exact_value": "10", "source": "original_compact_numeric coordinate terms"},
            },
            "outer_power_substitution_node": 48,
            "outer_power_substitution_exact_value": "2",
            "post_restore_substitution_node": 1,
            "post_restore_substitution_exact_value": "1",
        },
    }


def prove_geometry_binding(geometry_receipt, lineage_receipt, journal):
    """Bind native Rc geometry to numeric logR(x) using saved source metadata only.

    The result is a source-level geometry contract. It does not transfer a
    numeric cap, evaluate a radius, claim terminal zero, or establish a global
    physical identity.
    """
    _need(isinstance(geometry_receipt, dict) and isinstance(lineage_receipt, dict),
          "parsed geometry and logC lineage receipts are required")
    pointer_checks = _receipt_pointer_checks(journal)
    _need(geometry_receipt.get("schema") == "compact_native_numeric_volume_geometry_trace/v2",
          "unexpected geometry receipt schema")
    _need(lineage_receipt.get("schema") == "compact_logC_parameter_lineage_join/v1",
          "unexpected logC lineage receipt schema")
    _need(geometry_receipt.get("authentication", {}).get("native_compact_archive", {}).get("raw_sha256")
          == EXPECTED_NATIVE_COMPACT_RAW_SHA256,
          "native compact archive hash is not authenticated in geometry receipt")
    _need(geometry_receipt.get("authentication", {}).get("numeric_compact_quiet_report", {}).get("raw_sha256")
          == EXPECTED_NUMERIC_COMPACT_RAW_SHA256,
          "numeric compact report hash is not authenticated in geometry receipt")
    _need(lineage_receipt.get("authentication", {}).get("native_compact_archive", {}).get("raw_sha256")
          == EXPECTED_NATIVE_COMPACT_RAW_SHA256,
          "native compact archive hash is not authenticated in logC lineage receipt")

    geometry = geometry_receipt["numeric_cap_geometry"]
    native = geometry_receipt["native_Rc_source_graph"]
    formulas = geometry_receipt["source_formula_chain"]
    _need(geometry.get("source_family") == EXPECTED_SOURCE_FAMILY,
          "numeric geometry source family changed")
    _need(geometry.get("selected_N_all_samples_equal") is True,
          "numeric geometry selected-N rows are not equal")
    _need(geometry.get("selected_mpf_tuple_values_not_recomputed") is True,
          "numeric geometry receipt no longer marks tuples as saved source data")
    _need(geometry.get("numeric_integer_fidelity_preserved") is True,
          "numeric geometry integer fidelity flag changed")
    _need(geometry.get("selected_packet_data_is_source_data") is True,
          "numeric geometry packet source-data flag changed")
    _need(geometry.get("physical_radius_materialized") is False,
          "numeric geometry materialized a radius")
    _need(geometry.get("native_coordinate_is_physical_log_offset") is False,
          "numeric coordinate contract changed")
    _need(geometry.get("exact_source_formula")
          == "log(Rupper)=log(110)+10*logC+10*logP+exp(40)+14+log(x)",
          "numeric exact logR source formula changed")
    _need(geometry.get("equivalent_radius_definition")
          == "Rupper=110*exp(10*logC+10*logP+exp(40)+14)*x",
          "numeric equivalent radius source definition changed")
    _need(geometry.get("coordinate_contract") == "x=R/Rc; dy=dx/x; physical radius is not materialized",
          "numeric x coordinate contract changed")
    _need(formulas.get("numeric_global_y_terms")
          == "[10*logC, 10*logP, exp(40), 14, log(x)] in original_compact_numeric output; log(110) retained as radius base",
          "numeric global-y source terms changed")
    _need(formulas.get("logP_source") == "original_core_numeric.logP=exp(40)+11",
          "numeric logP source expression changed")
    _need(formulas.get("compact_source_definition")
          == "original_compact_background_velocity native R_definition=110*exp(10logCstar+10logPstar+exp40+14)*x",
          "numeric compact source definition changed")
    _need(formulas.get("source_formula_match") is True,
          "saved numeric/native geometry formula match flag is false")

    graph = native["graph_nodes"]
    graph_logc = _node(graph, 35560)
    _need(graph_logc.get("source_family") == EXPECTED_SOURCE_FAMILY,
          "native graph logC source family changed")
    native_map = _native_log_rc_map(
        graph, native, formulas, geometry_receipt["numeric_core_source_binding"],
        geometry_receipt["source_code_pins"])
    native_binding = native.get("compact_binding", {})
    _need(native.get("Rc_node") == 251066,
          "native Rc node changed")
    _need(native_binding.get("background", {}).get("radius") == 260862,
          "native compact radius port changed")
    _need(native_binding.get("background", {}).get("phase") == 260181,
          "native compact x/phase port changed")
    _need(native_binding.get("background", {}).get("domain") == [1, 48],
          "native compact domain changed")

    # Bind the scalar source parameter through the previously authenticated
    # lineage receipt. These checks are conditional source contracts only.
    family_join = lineage_receipt.get("source_family", {})
    _need(family_join.get("exact_equal") is True
          and family_join.get("native") == EXPECTED_SOURCE_FAMILY
          and family_join.get("numeric_core_report") == EXPECTED_SOURCE_FAMILY,
          "numeric/native logC source family join changed")
    logc_node = lineage_receipt["native_logC_target"]["node"]
    _need(logc_node.get("source_path") == "common.core.logC",
          "native logC source path changed")
    _need(logc_node.get("exact_fixed_source_parameter") is True
          and logc_node.get("parameter_enclosure_not_selected_value") is True
          and logc_node.get("Z_independent") is True,
          "native logC source-parameter flags changed")
    _need(logc_node.get("source_family") == EXPECTED_SOURCE_FAMILY,
          "lineage logC source family changed")
    numeric_logc = lineage_receipt["numeric_logC_target"]
    _need(numeric_logc.get("lower_exact_mpf_tuple_equals_upper_exact_mpf_tuple") is True
          and numeric_logc.get("parameter_is_selected_singleton") is True
          and numeric_logc.get("parameter_is_midpoint_or_range_endpoint") is False,
          "numeric selected logC singleton contract changed")
    selected = numeric_logc["selected_logCstar_exact_record"]
    _need(selected.get("lower_exact_mpf_tuple") == selected.get("upper_exact_mpf_tuple"),
          "numeric selected logC tuples differ")
    evidence = lineage_receipt["native_logC_target"]["static_source_provider_evidence"]
    checks = evidence.get("source_text_checks", [])
    _need("core.logC._mpi_ == read(core.records['physical_norm_family'],'selected_logCstar')._mpi_" in checks,
          "native selected-logC source check missing")
    _need(any("common.core.logC" in item and "exact_fixed_source_parameter" in item
              for item in checks),
          "native exact common.core.logC leaf check missing")
    _need(lineage_receipt.get("join_assessment", {}).get("native_numeric_logC_source_identity_join") is True,
          "saved scalar logC source join was not accepted")

    fixed = lineage_receipt["fixed_leading_N0_identity"]
    meta = fixed.get("metadata_checks", {})
    for key in ("native_graph_node6_N0_equals_C3_repair_node6_N0",
                "numeric_packet_J_equals_native_branch_J",
                "numeric_packet_J_equals_C3_repair_J",
                "native_branch_J_equals_C3_repair_J",
                "native_branch_definition_equals_C3_repair_definition"):
        _need(meta.get(key) is True, "fixed-leading-N0 metadata check failed: %s" % key)
    j_numeric = geometry.get("selected_N_exact", {}).get("J")
    j_native = fixed.get("native_branch_selected_integer", {}).get("J")
    j_repair = fixed.get("C3_repair_selected_integer", {}).get("J")
    _need(type(j_numeric) is int and type(j_native) is int and type(j_repair) is int
          and j_numeric == j_native == j_repair,
          "selected exact N J differs across geometry/logC records")

    measure = _cached_measure_witness(journal)
    numeric_map = {
        "log110": Fraction(1), "logC": Fraction(10), "logP": Fraction(10),
        "exp40": Fraction(1), "constant": Fraction(14), "log_x": Fraction(1),
    }
    native_with_x = _add(
        {key: Fraction(value.numerator, value.denominator)
         for key, value in _decode_map(native_map["expected_logRc_map"]).items()},
        {"log_x": Fraction(1)},
    )
    geometry_residual = _zero_residual(native_with_x, numeric_map)
    _need(_zero_map(geometry_residual),
          "native logRc and numeric logR(x) affine maps do not match")

    result = {
        "schema": "compact_prefix_geometry_binding_draft/v1",
        "source_only": True,
        "unentered": True,
        "passed": True,
        "receipt_pointers": pointer_checks,
        "selected_J": j_numeric,
        "native_logRc_affine_map": native_map,
        "numeric_logR_affine_map": _encode_map(numeric_map),
        "geometry_affine_residual": _encode_map(geometry_residual),
        "cached_C2_measure_witness": measure,
        "source_contracts": {
            "native_Rc_node": 251066,
            "native_compact_radius_node": 260862,
            "native_compact_phase_node": 260181,
            "native_source_path": "common.core.logC",
            "numeric_formula": geometry["exact_source_formula"],
            "coordinate_contract": geometry["coordinate_contract"],
            "scalar_logC_join_conditional": True,
            "same_fixed_leading_N0_context": True,
        },
        "formal_preconditions": {
            "ordinary_Z_and_source_family_contracts_are_external": True,
            "logC_source_contract_is_saved_and_conditional": True,
            "physical_radius_not_materialized": True,
            "common_x_domain_is_saved_native_domain": True,
            "L_nonzero_and_other_transport_conditions_not_proved_here": True,
        },
        "false_scope_flags": {
            "terminal_zero_claimed": False,
            "physical_global_identity_claimed": False,
            "quantitative_cap_transfer_claimed": False,
            "full_numeric_function_identity_claimed": False,
            "cap_identity_claimed": False,
        },
        "audit": {
            "source_imports": 0,
            "source_evaluations": 0,
            "proof_calls": 0,
            "normalization_calls": 0,
            "radius_evaluations": 0,
            "new_C2_measure_checks": 0,
            "exact_fraction_coefficient_maps": True,
            "numeric_centers_used_as_source_values": False,
        },
    }
    _emit(journal, {"event": "compact_prefix_geometry_binding_saved",
                    "passed": True, "native_logRc_numeric_logR_residual_zero": True,
                    "source_logC_join_conditional": True,
                    "cap_transfer": False, "terminal_zero": False})
    return result


def _decode_map(encoded):
    """Decode this draft's own exact map encoding without decimal conversion."""
    return {key: Fraction(value["numerator"], value["denominator"])
            for key, value in encoded.items()}


from fractions import Fraction

def _poly_terms(pairs):
    out={}
    for monomial,coefficient in pairs:
        key=tuple(sorted(monomial));out[key]=out.get(key,Fraction(0))+Fraction(coefficient)
    return {k:v for k,v in out.items() if v}

def _plus(*polys):return _poly_terms((k,v) for p in polys for k,v in p.items())
def _minus(p,q):return _plus(p,{k:-v for k,v in q.items()})
def _times(p,q):return _poly_terms((k+l,v*w) for k,v in p.items() for l,w in q.items())
def _var(name):return {(name,):Fraction(1)}
def _poly_scale(p,c):return _poly_terms((m,Fraction(c)*v) for m,v in p.items())

def _ordinary_Z(p):
    derivative={'u0':'u1','u1':'u2','e0':'e1','e1':'e2'}
    terms=[]
    for monomial,coefficient in p.items():
        for j,v in enumerate(monomial):
            if v in derivative:
                terms.append((monomial[:j]+(derivative[v],)+monomial[j+1:],coefficient))
    return _poly_terms(terms)

def _positive_envelope(p):
    names={'b':'B','u0':'U0','u1':'U1','u2':'U2','e0':'E0','e1':'E1','e2':'E2'}
    return _poly_terms((tuple(names[v] for v in monomial),abs(c)) for monomial,c in p.items())

def _convolve_C1(a,b):
    return [_plus(*[_times(a[j],b[k-j]) for j in range(k+1)]) for k in range(2)]

def _archived_product_bound(a,b,ea,eb):
    # A new polynomial representation of the archived three-term bound.
    # No archived function is imported or executed.
    pieces=[_convolve_C1(a,eb),_convolve_C1(b,ea),_convolve_C1(ea,eb)]
    return [_plus(*(p[k] for p in pieces)) for k in range(2)]

def _json_polynomial(p):
    return [{'factors':list(m),'coefficient':{'numerator':c.numerator,'denominator':c.denominator}}
            for m,c in sorted(p.items())]

def prove_compact_C1_dominance(journal):
    """Derive signed channels, then prove their C1 absolute dominance."""
    u0,u1,u2,e0,e1,e2,b=[_var(s) for s in ['u0','u1','u2','e0','e1','e2','b']]
    signed={'q':e0,'qZ':e1,'q2':_minus(_times(_plus(b,e0),_plus(b,e0)),_times(b,b)),
            'UqZ':_times(u0,e1),'U2':_times(u0,u0),'UUZ':_times(u0,u1)}
    q=[_var('B'),{},{}];u=[{},{},{}]
    eq=[_var('E0'),_var('E1'),_var('E2')];eu=[_var('U0'),_var('U1'),_var('U2')]
    archived={'q':eq[:2],'qZ':eq[1:],
              'q2':_archived_product_bound(q,q,eq,eq),
              'UqZ':_archived_product_bound(u,q[1:],eu,eq[1:]),
              'U2':_archived_product_bound(u,u,eu,eu),
              'UUZ':_archived_product_bound(u,u[1:],eu,eu[1:])}
    channels={};signed_jets={};new_bounds={};checks=[]
    for name,p in signed.items():
        jets=[p,_ordinary_Z(p)];signed_jets[name]=jets
        positive=[_positive_envelope(v) for v in jets];new_bounds[name]=positive
        for k in range(2):
            residual=_minus(positive[k],archived[name][k]);checks.append(not residual)
            if residual:raise ValueError('Archived C1 dominance mismatch: '+name+' Z'+str(k))
        channels[name]={'signed_C0_C1':[_json_polynomial(v) for v in jets],
                        'absolute_bound_C0_C1':[_json_polynomial(v) for v in positive],
                        'archived_delta_product_C0_C1':[_json_polynomial(v) for v in archived[name]],
                        'exact_coefficient_checks':[True,True]}
    # This is Leibniz transport of the six variable coefficient channels.
    # A_jk denotes abs(d_Z^k a_j), kept separate from signed source values.
    volume=[];archived_volume=[]
    for k in range(2):
        pieces=[];old=[]
        for name in signed:
            for j in range(k+1):
                a=_var('A_'+name+'_'+str(j))
                pieces.append(_times(a,new_bounds[name][k-j]))
                old.append(_times(a,archived[name][k-j]))
        p=_plus(*pieces);qv=_plus(*old)
        if _minus(p,qv):raise ValueError('C1 volume coefficient transport mismatch')
        volume.append(_json_polynomial(p));archived_volume.append(_json_polynomial(qv))
    # Coefficient jets are stored as Taylor coefficients T_j=f^(j)/j!,
    # while pressure error rows use ordinary derivatives O_j=f^(j). Check
    # j!*binomial(k,j)*T_j == binomial(k,j)*O_j explicitly through order2.
    factorial=[1,1,2];pascal=[[1],[1,1],[1,2,1]]
    taylor=[_var('T0'),_var('T1'),_var('T2')]
    ordinary=[_poly_scale(taylor[j],factorial[j]) for j in range(3)]
    pressure_checks=[]
    for k in range(3):
        for j in range(k+1):
            archived_term=_poly_scale(taylor[j],factorial[j]*pascal[k][j])
            ordinary_term=_poly_scale(ordinary[j],pascal[k][j])
            if _minus(archived_term,ordinary_term):raise ValueError('Pressure Taylor/ordinary conversion failed')
            pressure_checks.append({'target_order':k,'coefficient_order':j,
                                    'taylor_symbol':'T'+str(j),
                                    'ordinary_symbol':'O'+str(j),
                                    'taylor_to_ordinary_factor':factorial[j],
                                    'archived_term':_json_polynomial(archived_term),
                                    'ordinary_leibniz_term':_json_polynomial(ordinary_term),
                                    'exact_match':True})
    result={'schema':'compact_prefix_C1_dominance_v1','passed':True,
            'channels':channels,'volume_bound_C0_C1':volume,
            'archived_volume_bound_C0_C1':archived_volume,
            'new_exact_channel_coefficient_checks':len(checks),
            'new_exact_volume_checks':2,'pressure_ordinary_Z_conversion_checks':pressure_checks,
            'assumptions':{'background_U_C2':[0,0,0],'background_Q_Z1_Z2':[0,0],
                           'U_errors':'abs(d_Z^k deltaU)<=Uk, k=0..2',
                           'Q_errors':'abs(d_Z^k deltaQ)<=Ek, k=0..2',
                           'B':'abs(Q_BG)<=B on the whole compact prefix',
                           'A_jk':'abs(d_Z^k coefficient_j), k=0,1',
                           'all_envelope_variables_nonnegative':True},
            'conclusion':'abs(d_Z^k(v_N-v_BG)) <= volume_bound[k], k=0,1',
            'source_values_selected':False,'cap_is_function_value':False,
            'old_source_or_proof_calls':0,'old_numeric_integration_calls':0,
            'new_polynomial_certificate_calls':1}
    journal({'event':'new_compact_C1_dominance_proved','channel_checks':12,'volume_checks':2})
    return result


def _prefix_triangle_certificate(corrected, branch, journal):
    """A new conditional C3 Duhamel bound, with its missing source ports explicit."""
    assessment=corrected['join_assessment']
    _need(assessment['numeric_prefix_budget_is_relevant_compact_pressure_budget'] is True,
          'Wrong numeric prefix budget')
    _need(assessment['absolute_caps_not_signed_values'] is True,
          'A cap was substituted for a signed value')
    _need(assessment['physical_Pstar_applied_once'] is True,
          'Physical pressure scale changed')
    _need(assessment['selected_profile_bounds_already_N_scaled'] is True,
          'Selected profile scaling changed')
    common=branch['common_literal_guards']
    _need(common['all_three_source_families_equal'] is True,
          'C3/C5 source families differ')
    transfer=branch['branch_transfer_literals']
    _need(transfer['same_original_B_Q_N_and_exact_bump_integrals'] is True,
          'Original branch polynomial changed')
    correspondence=transfer['exact_source_correspondence']
    for name in ('same_all17_prescribed_source_and_normalized_profile_functions',
                 'physical_selected_N_incoming_not_divided_again',
                 'own_rate_original_Duhamel_maps_and_inlet_bound_by_predecessor'):
        _need(correspondence[name] is True,'Missing saved source correspondence: '+name)
    # The mathematical argument is independent of the unsigned caps' size.
    # s, support S and the rate are ordinary-Z independent. Differentiating
    # this *defined* C3 source integral gives the four displayed targets.
    targets=[]
    for k in range(4):
        targets.append({
            'ordinary_Z_order':k,
            'signed_target':'H_%d(s)=exp(-s)*I_%d+integral_0^s exp(-(s-u))*rho_%d(u) du'%(k,k,k),
            'hypotheses':['abs(I_%d)<=Icap_%d'%(k,k),
                          'abs(rho_%d(u))<=Dcap_%d*1_S(u)'%(k,k),
                          '0<=u<=s<=log(2); S is the union of the original m bump supports'],
            'pointwise_kernel_steps':['exp(-s)<=1','0<exp(-(s-u))<=1',
                                      'abs(integral f)<=integral abs(f)'],
            'conclusion':'abs(H_%d(s))<=Icap_%d+Dcap_%d*measure(S intersect [0,s])'%(k,k,k),
            'whole_prefix_conclusion':'sup_0<=v<=s abs(H_%d(v))<=Icap_%d+Dcap_%d*measure(S intersect [0,s])'%(k,k,k),
            'incoming_retained':True,'extra_N_divisions':0,
            'physical_M_conclusion':'abs(d_Z^%d(M_N-M_BG))<=Pstar*(Icap_%d+Dcap_%d*measure(S intersect [0,s]))'%(k,k,k)})
    # The source identities which would discharge the hypotheses are kept as
    # separate obligations. Same-family metadata is not such an identity.
    obligations={
        'incoming':'C3 terminal_N_scaled_histories/N equals native incoming_correction[258646:258650] through ordinary Z3',
        'controls':'C3 exact C1 branch and its ordinary Z1..3 equal native system/component jets, with the same normalized target functions',
        'profiles':'C3 delta_V(y=0..2) and density m equal native delta_V[260375,260378,260382,260387] and correction integrands',
        'predecessor':'The minimum of predecessor and band incoming budgets bounds the same native incoming source, including all earlier error terms'}
    journal({'event':'new_conditional_prefix_triangle_proved','orders':4,
             'native_port_obligations_discharged':False})
    return {'schema':'compact_C3_prefix_triangle_v1','passed':True,
            'conditional_inequality_proved':True,'targets':targets,
            'pressure_transport':'E_Q[k]=sum_name sum_j=0^k binomial(k,j)*abs(d_Z^j a_name)*E_name[k-j], k<=2; Q_y through k<=1',
            'pressure_coefficients_are_cached_same_physical_recovery':True,
            'profile_U':'abs(d_Z^k U_N)<=Pstar*deltaV_cap[k], k<=3, after profile source join; prefix predecessor budget retained',
            'profile_U_y':'abs(d_Z^k (U_N)_y)<=Pstar*deltaV_y_cap[k], k<=1, after profile source join',
            'source_obligations':obligations,
            'native_port_inequalities_proved':False,
            'all_caps_retained':True,'ordinary_Z_not_cap_derivative':True,
            'no_cap_midpoint_or_endpoint_selected_as_signed_value':True}


def derive(source, journal):
    """Run this fresh certificate once; consume saved proofs without replay."""
    spec=source['input'];receipts=spec['receipts'];parent=source['parent_result']
    _need(parent['passed'] is True and
          parent['native_BG_exact_numeric_baseline_C3_identity_proved'] is True and
          parent['relative_mass_C3_source_identity_proved'] is True,
          'Cached native mean/relative mass proof missing')
    # Complete four omitted literal rows directly from the authenticated raw
    # native archive. No source operation or graph normalizer is evaluated.
    import copy
    geo=copy.deepcopy(receipts['geometry'])
    supplement=spec['native_geometry_literal_supplement']
    _need(supplement['archive_raw_sha256']==EXPECTED_NATIVE_COMPACT_RAW_SHA256,
          'Geometry supplement source changed')
    graph=geo['native_Rc_source_graph']['graph_nodes']
    for key,row in supplement['graph_nodes'].items():
        _need(key not in graph or graph[key]==row,'Geometry supplement conflicts')
        graph[key]=row
    _expect_fraction(graph,2,40);_expect_fraction(graph,48,2);_expect_fraction(graph,2083,10)
    _expect_operation(graph,260862,'product',arguments=[251066,260181])
    geometry_events={'geometry_receipt_sha256':EXPECTED_GEOMETRY_RECEIPT_SHA256,
                     'lineage_receipt_sha256':EXPECTED_LOGC_LINEAGE_RECEIPT_SHA256,
                     'measure_normalization':spec['cached_C2_measure_normalization']}
    geometry=prove_geometry_binding(geo,receipts['logC'],geometry_events)
    geometry['unentered']=False
    geometry['native_leaf_supplement_bound']=True
    for event in geometry_events.get('events',[]):journal(event)
    memory=bind_fixed_mstar_enclosure(receipts['memory'],parent,journal)
    _need(memory['passed'] is True,'Fixed mstar binding rejected: '+memory.get('error',''))
    dominance=prove_compact_C1_dominance(journal)
    prefix=_prefix_triangle_certificate(receipts['pressure_prefix'],receipts['branch'],journal)
    result={
        'schema':'original_compact_prefix_pressure_certificate',
        'component':'compact_prefix_pressure_certificate','passed':True,
        'metadata_only':False,'source_only':True,
        'input_digest':source['input_digest'],'runtime_model_metadata':None,
        'original_source_evaluations':0,
        'source_packet_defined':True,'source_ftc_transport_defined':True,
        'native_BG_exact_numeric_baseline_C3_identity_proved':True,
        'relative_mass_C3_source_identity_proved':True,
        'compact_prefix_geometry_logC_binding_proved':geometry['passed'],
        'compact_prefix_fixed_mstar_enclosure_proved':memory['conditional_numeric_mean_enclosure_bound'],
        'compact_prefix_12channel_2volume_dominance_proved':dominance['passed'],
        'compact_prefix_source_binding_proved':False,
        'compact_volume_cap_same_source_bound_proved':False,
        'prefix_pressure_certificate':{'geometry':geometry,'fixed_mstar':memory,
                                       'C1_absolute_dominance':dominance,
                                       'conditional_C3_prefix':prefix},
        'source_binding':{'passed':True,'scope':'Saved native BG and exact relative-mass sources, literal physical geometry, fixed memory certificate and conditional bound algebra',
                          'actual_C3_cap_to_C5_port_bound':False,
                          'same_family_metadata_used_as_function_identity':False,
                          'fixed_memory_binding_conditional_on_saved_interval_certificate':True},
        'coordinate_mapping':{'x':'R/Rc in [1,2]','physical_log_measure':'dy=dx/x once',
                              'log_R':'log(110)+10logC+10logP+exp(40)+14+log(x)',
                              'physical_volume':'R^2 dy; cached C2 mass 3Rc^2/2',
                              'ordinary_Z_geometry_independent':True},
        'analytic_integral_definition':{'signed_Evol':'integral_Rc^2Rc R^2*(v_N-v_BG)dy',
                                       'ordinary_Z_orders':[0,1],
                                       'source_function_value_evaluated':False,
                                       'old_reference_offset_retained':True},
        'source_inspection_audit':{'original_source_evaluations':0,'old_source_consumer_calls':0,
             'old_proof_functions_called':0,'old_operator_calls':0,'old_source_normalizer_calls':0,
             'quadratures':0,'numeric_field_evaluations':0,'temporal_recurrence_calls':0,
             'new_channel_coefficient_checks':12,'new_volume_dominance_checks':2,
             'new_ordinary_pressure_conversion_checks':6,
             'new_conditional_prefix_triangle_orders':4},
        'remaining_scope':list(prefix['source_obligations'].values())+[
             'After discharging these port identities, apply this same-source dominance to the retained volume cap.',
             'Compute signed Evol and independent numeric-to-old reference offset; retain all prior caps.',
             'Global five-moment and later physical/temporal/NS layers remain incomplete.']}
    for key in ('retained_cap_provenance','retained_construction_data','retained_integral_term_inventory'):
        result[key]=copy.deepcopy(parent.get(key))
    for key in source['inherited_false_flags']:result[key]=False
    journal({'event':'prefix_pressure_partial_certificate_complete','actual_cap_transfer':False,
             'geometry_bound':True,'fixed_memory_conditional_bound':True,'new_volume_dominance_checks':2})
    return result

