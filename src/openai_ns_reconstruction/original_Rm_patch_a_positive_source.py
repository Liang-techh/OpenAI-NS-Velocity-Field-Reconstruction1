"""Source-only certificate for the actual leading ``Rm`` loop-root ``a``.

The certificate joins the saved patch DAG to the checked C6 source constructors
and to the old correlated C1 whole-domain theorem.  It never evaluates a field,
replays a producer or checker, or treats an endpoint/profile value as a bound.
The only promoted lower bound is the source theorem's exact ``7/10`` bound.
"""
from __future__ import annotations

import ast
from collections.abc import Mapping, Sequence
from copy import deepcopy
import hashlib
import json
from typing import Any


SCHEMA = "original_Rm_patch_a_positive_source_v1"
PREFIX = "lei_ren_part1_paper_compliant_current_original_"
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
EXPECTED_NODE_COUNT = 205300

Z_NODE = 20
CHART_NODE = 185676
RADIUS_NODE = 187214
PHYSICAL_JACOBIAN_NODE = 187215
E_NODE = 194957
E_Y_NODE = 195164
A_NODE = 196071

SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

# Every pin is an input hash in the immutable parent report.  No gzip archive
# is consumed here; the executable source and JSON receipts are sufficient.
RECORD_HASHES = {
    PREFIX + "true_Rm_patch_C6_functions.py": "a6f1cb5deaf17e9482aa964aef7867aadbc815e4d703a415e4945cf665f14002",
    PREFIX + "true_Rm_patch_C6_functions_check.py": "85f42a1476b0af3df91a8fd3f1077f1e0d7109cc9f25ad24fb9205ac891cb1cf",
    PREFIX + "true_Rm_patch_C6_functions_check.json": "0190fcf738fb8add5be85cf0767464a39cc0d18240e6b36ad4ec7930dce38aa7",
    "lei_ren_part1_paper_compliant_current_patch_relaxed_inputs.py": "686f45a65a77e874d347936eb5aaee7aeb1024f9251d514e0a940285975d5f18",
    "lei_ren_part1_paper_compliant_current_patch_relaxed_inputs_check.py": "1e0a7c69fb5de1603e688723092b877c45d6dcd6bafd5c18d062da4315448f8a",
    "lei_ren_part1_paper_compliant_current_patch_relaxed_inputs_check.json": "b96bddd5ac67a0930d7f22d8c5ba8bb04da339e19636ae797de6b8b416acb9a3",
    "lei_ren_part1_paper_compliant_current_patch_relaxed_inputs.json": "b075c27d4faae5cfcc181ebd3e3e03ce86fb7ed74cea53fc4b634003428c2215",
    "lei_ren_part1_paper_compliant_current_generic_shear_source_bounds.py": "a590cbf346ddbf3cc6298e16b4c28fe7ed330d0acdb14e60d280b83e9ddf7adf",
    "lei_ren_part1_paper_compliant_current_generic_shear_source_bounds_check.py": "e59bd81f2f9e96855764b6d2e781ac0914ff0011751174e0888f17aebe7bf759",
    "lei_ren_part1_paper_compliant_current_generic_shear_source_bounds_check.json": "7350b72291cb6278afaafd407685eb006fd17b503abb7c5f5f852aab878b804c",
    "lei_ren_part1_paper_compliant_current_generic_shear_source_bounds.json": "2817a7849a3bd943acb2e0ffaebae53d0f4d2c396b0c15f0f9b4d6effa8364bd",
    "lei_ren_part1_paper_shared_bump_constants.py": "1c986f72e6967db685d8a464377d4b7b8c679a6eea2bf16aa138702f065f5ced",
    "lei_ren_part1_paper_shared_bump_constants.json": "9dd940dab760c8d1976e6b543bd2170e6374400c142b8d1e2d09abdf90abf77f",
    "lei_ren_part1_paper_compliant_five_moment_repair.py": "978de37342cc2a63aa7c1bfda573e512f344e0157aed28b6017c7dae2ff8eba5",
    "lei_ren_part1_paper_compliant_five_moment_repair_check.py": "2878e3698533208ee07433bfa3a8f89093f55c4d37a766b3d249258b73e2fbf0",
    "lei_ren_part1_paper_compliant_five_moment_repair_check.json": "0a2d901c4fbe5afe822752d190d02d361a72ab647e2c009944935851958f0bdf",
    "lei_ren_part1_paper_compliant_five_moment_repair.json": "87f3f6433bd7be901ab7d7b4cfdb529a088d951da05fb6bbeb2b2d2be4a9ab4d",
    "lei_ren_part1_paper_compliant_actual_feedback_moment_patch.py": "ce009ad8b6b2eece7801660440c8ea25c4aa60d53a3c7b83ddfaebefa8e6c000",
    "lei_ren_part1_paper_compliant_actual_feedback_moment_patch_check.py": "ebf9d90820c3c986ff31180657bd1e7086418d0081cd6cb75eb761f7bea19076",
    "lei_ren_part1_paper_compliant_actual_feedback_moment_patch_check.json": "f094974b3870cad479fb4d981d49c34089970eea76036c3d7869d06b83e24018",
    "lei_ren_part1_paper_compliant_actual_feedback_moment_patch.json": "fedfda0471eefc799c162ce0070e3980bdc2aa1d3f5d218d1932e62abb9a9d54",
}


def canonical_digest(value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"),
                     ensure_ascii=True, allow_nan=False).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _row(nodes: Sequence[Any], index: int) -> Mapping[str, Any]:
    _require(type(index) is int and 0 <= index < len(nodes), f"graph handle out of range: {index}")
    row = nodes[index]
    _require(isinstance(row, Mapping), f"graph handle is not a mapping: {index}")
    return row


def _expect_row(nodes: Sequence[Any], index: int, operation: str, **fields: Any) -> dict[str, Any]:
    row = _row(nodes, index)
    _require(row.get("operation") == operation, f"operation changed at graph handle {index}")
    for key, value in fields.items():
        _require(row.get(key) == value, f"source field {key} changed at graph handle {index}")
    return {"handle": index, "operation": operation, "fields": deepcopy(fields)}


def _record(source_records: Mapping[str, Any], name: str, gates: list[dict[str, Any]]) -> str | None:
    expected = RECORD_HASHES[name]
    record = source_records.get(name) if isinstance(source_records, Mapping) else None
    if not isinstance(record, Mapping):
        gates.append({"gate": "source_record_missing", "record": name})
        return None
    expected_path = "experiments/root_st073/" + name
    text = record.get("source_text")
    actual = hashlib.sha256(text.encode("utf-8")).hexdigest() if isinstance(text, str) else None
    if record.get("path") != expected_path:
        gates.append({"gate": "source_record_path_changed", "record": name, "expected": expected_path})
        return None
    if record.get("sha256") != expected or actual != expected:
        gates.append({"gate": "source_record_hash_changed", "record": name, "expected_sha256": expected})
        return None
    return text


def _json(text: str, name: str, gates: list[dict[str, Any]]) -> Mapping[str, Any] | None:
    try:
        value = json.loads(text, parse_int=str, parse_float=str)
    except (TypeError, json.JSONDecodeError) as exc:
        gates.append({"gate": "source_receipt_invalid", "record": name, "error": str(exc)})
        return None
    if not isinstance(value, Mapping):
        gates.append({"gate": "source_receipt_not_mapping", "record": name})
        return None
    return value


def _path(value: Mapping[str, Any], path: tuple[str, ...]) -> Any:
    current: Any = value
    for key in path:
        current = current.get(key) if isinstance(current, Mapping) else None
    return current


def _receipt(text: str, name: str, required: Mapping[tuple[str, ...], Any], gates: list[dict[str, Any]]) -> Mapping[str, Any] | None:
    value = _json(text, name, gates)
    if value is None:
        return None
    ok = value.get("all_passed") is True
    for path, expected in required.items():
        actual = _path(value, path)
        if actual != expected:
            gates.append({"gate": "accepted_source_receipt_flag_missing", "record": name,
                          "path": ".".join(path), "expected": expected, "actual": actual})
            ok = False
    return value if ok else None


def _definitions(text: str, name: str, wanted: tuple[str, ...], gates: list[dict[str, Any]]) -> dict[str, ast.AST]:
    try:
        tree = ast.parse(text, filename=name)
    except SyntaxError as exc:
        gates.append({"gate": "source_ast_invalid", "record": name, "error": str(exc)})
        return {}
    found: dict[str, ast.AST] = {}
    for item in ast.walk(tree):
        if isinstance(item, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) and item.name in wanted:
            found[item.name] = item
    for item in wanted:
        if item not in found:
            gates.append({"gate": "source_declaration_missing", "record": name, "declaration": item})
    return found


def _tokens(text: str, name: str, wanted: tuple[str, ...], gates: list[dict[str, Any]]) -> bool:
    missing = [token for token in wanted if token not in text]
    if missing:
        gates.append({"gate": "source_constructor_token_missing", "record": name, "tokens": missing})
        return False
    return True


def _call(node: ast.AST, attribute: str) -> ast.Call | None:
    if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
        return None
    return node if node.func.attr == attribute else None


def _patch_constructor_shape(text: str, name: str, gates: list[dict[str, Any]]) -> bool:
    try:
        tree = ast.parse(text, filename=name)
    except SyntaxError:
        return False
    patch = next((x for x in ast.walk(tree) if isinstance(x, ast.FunctionDef) and x.name == "patch_fields"), None)
    if patch is None:
        return False
    assignments: dict[str, ast.AST] = {}
    for item in ast.walk(patch):
        if isinstance(item, ast.Assign) and len(item.targets) == 1 and isinstance(item.targets[0], ast.Name):
            assignments[item.targets[0].id] = item.value
    fields = assignments.get("fields")
    radial = assignments.get("radial")
    if not isinstance(fields, ast.Call) or not isinstance(fields.func, ast.Name) or fields.func.id != "dict":
        gates.append({"gate": "C6_patch_fields_constructor_missing", "record": name, "field": "fields=dict"})
        return False
    e_kw = next((k for k in fields.keywords if k.arg == "E"), None)
    if e_kw is None or _call(e_kw.value, "mul") is None:
        gates.append({"gate": "C6_patch_fields_constructor_changed", "record": name, "field": "E=am*H"})
        return False
    e_args = _call(e_kw.value, "mul").args
    if len(e_args) != 2 or not all(isinstance(x, ast.Name) for x in e_args) or [x.id for x in e_args] != ["am", "H"]:
        gates.append({"gate": "C6_patch_fields_constructor_changed", "record": name, "field": "E=am*H"})
        return False
    if not isinstance(radial, ast.Call) or not isinstance(radial.func, ast.Name) or radial.func.id != "dict":
        gates.append({"gate": "C6_patch_radial_constructor_missing", "record": name})
        return False
    ey_kw = next((k for k in radial.keywords if k.arg == "E_y"), None)
    ey = _call(ey_kw.value, "scale") if ey_kw is not None else None
    inner = _call(ey.args[0], "mul") if ey is not None and ey.args else None
    ok = (ey is not None and len(ey.args) == 2 and isinstance(ey.args[1], ast.Name)
          and ey.args[1].id == "x" and inner is not None and len(inner.args) == 2
          and [x.id for x in inner.args if isinstance(x, ast.Name)] == ["am", "Hx"])
    if not ok:
        gates.append({"gate": "C6_patch_radial_constructor_changed", "record": name, "field": "E_y=x*am*Hx"})
    return ok


def _source_witnesses(source_records: Mapping[str, Any], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
    texts = {name: _record(source_records, name, gates) for name in RECORD_HASHES}
    if any(text is None for text in texts.values()):
        return False, {"records": dict(RECORD_HASHES)}

    c6 = PREFIX + "true_Rm_patch_C6_functions.py"
    c6_check = PREFIX + "true_Rm_patch_C6_functions_check.json"
    relaxed = "lei_ren_part1_paper_compliant_current_patch_relaxed_inputs.py"
    relaxed_check = "lei_ren_part1_paper_compliant_current_patch_relaxed_inputs_check.json"
    relaxed_data = "lei_ren_part1_paper_compliant_current_patch_relaxed_inputs.json"
    bounds = "lei_ren_part1_paper_compliant_current_generic_shear_source_bounds.py"
    bounds_check = "lei_ren_part1_paper_compliant_current_generic_shear_source_bounds_check.json"
    bounds_data = "lei_ren_part1_paper_compliant_current_generic_shear_source_bounds.json"
    feedback = "lei_ren_part1_paper_compliant_actual_feedback_moment_patch.py"
    feedback_check = "lei_ren_part1_paper_compliant_actual_feedback_moment_patch_check.json"
    repair = "lei_ren_part1_paper_compliant_five_moment_repair.py"
    repair_check = "lei_ren_part1_paper_compliant_five_moment_repair_check.json"
    repair_data = "lei_ren_part1_paper_compliant_five_moment_repair.json"
    bump_data = "lei_ren_part1_paper_shared_bump_constants.json"

    declarations = {
        c6: _definitions(texts[c6], c6, ("raw_beta", "ExactBumpWeights", "map_conjugacy", "patch_fields",
                                          "CurrentOriginalTrueRmPatchC6Functions"), gates),
        relaxed: _definitions(texts[relaxed], relaxed, ("patch_source_theorem", "CurrentPatchRelaxedInputs"), gates),
        bounds: _definitions(texts[bounds], bounds, ("positive_denominators", "CurrentGenericSourceBounds"), gates),
        feedback: _definitions(texts[feedback], feedback, ("current_patch_source_bindings", "CompliantActualFeedbackMomentPatch"), gates),
        repair: _definitions(texts[repair], repair, ("SharedFiveMomentRepair",), gates),
    }
    token_ok = all((
        _tokens(texts[c6], c6, (
            "g.constant(40)", "g.constant(1600)", "self.normalization",
            "positive exact original beta normalization", "fields=dict(E=a.mul(am,H),",
            "radial=dict(E_y=a.scale(a.mul(am,Hx),x)", "raw_to_original_normalized_conjugacy",
            "G(u)=S^-1 F(Su); J_G=S^-1 J_F S", "exact_original_beta_weights_normalization_centers_radius_and_support=True",
            "original_unique_common_box_branch_same_full_Jacobian_through_Z6=True",
            "actual_E_V_m_h_k_e_p_P0_pressure_and_Ey_Vy_C6=True",
            "H=a.add(a.fixed(power(g,x,Fraction(1,10))),*(a.scale(h[i+2],gamma[i][0]) for i in range(3)))",
            "Hx=a.add(a.fixed(g.mul(g.constant('1/10'),power(g,x,Fraction(-9,10))))",
            "*(a.scale(h[i+2],gamma[i][1]) for i in range(3)))",
            "repair=parent.linear.read(PREFIX+'five_moment_repair.json');t=repair['axial_box_scale']",
            "angular=g.mul(g.constant(10**10),self.t,self.t", "self.units=[self.t,self.t,angular,angular,angular]",
            "exact_source_correspondence=dict(", "actual_leading_N0_equal_by_exact_integer=True",
            "same_actual_L_Q_defect_and_amplitude_source=True",
            "exact_source_and_weight_branch_transfer=self.branch_transfer",
            "source_projection_binding=parent.active.current.ast_binding(admitted.WholeZRmCenteredInverse.owner)",
        ), gates),
        _patch_constructor_shape(texts[c6], c6, gates),
        _tokens(texts[relaxed], relaxed, (
            "'H':\"f+x**c.mpf('.1')\"", "'u':'am*H'", "'uy':\"am*(fx*x+x**c.mpf('.1')/10)\"",
            "zero('original_full_signed_a'", "a-(s.Rational(4,5)+(f/5-2*x*s.diff(f,x))/H)",
            "def bounds", "denominator=lower(1-B*t)", "amin,amax=c.mpf('.7'),c.mpf('.9')",
            "original_angular_a_bounds=[amin,amax]", "all Z[-1,1]", "source_radius_amplitude_and_width_not_materialized=True",
        ), gates),
        _tokens(texts[bounds], bounds, (
            "def positive_denominators", "minimum_a_from_current_original_whole_middle_theorems='.7'",
            "restoration_and_patch_proof='same checked actual original normalized swirl lower bounds on whole source domains'",
            "source_function_positivity_not_inferred_from_saved_box=True",
        ), gates),
        _tokens(texts[feedback], feedback, (
            "current_patch_source_bindings", "current_actual_moment_patch_installed",
            "current_actual_patch_implicit_axial5_recomputed", "current_five_functional_terminal_identities_connected",
        ), gates),
        _tokens(texts[repair], repair, ("class SharedFiveMomentRepair", "actual_five_defect_family_sha256", "self.W=weights"), gates),
    ))

    c6_value = _receipt(texts[c6_check], c6_check, {
        ("checks", "actual_source", "passed"): True,
        ("checks", "actual_source", "true_Rm_parent_and_same_original_map_nodes"): True,
        ("checks", "actual_source", "no_finite_N_history_reset_or_numeric_containment_claim"): True,
        ("checks", "exact_weights", "passed"): True,
        ("checks", "exact_weights", "each_exact_integral_containment_bound_by_its_pinned_directed_recipe_and_interval"): True,
        ("checks", "true_inlet_normalization", "exact_defining_function_identities_for_arbitrary_Z_raw_moments_Pstar_Rm_F0"): True,
        ("checks", "branch_theorem", "same_exact_original_unique_branch_theorem_bound"): True,
        ("checks", "branch_theorem", "exact_source_map_weights_amplitude_N0_and_initial_data_identified"): True,
        ("checks", "branch_theorem", "original_common_box_strict_Neumann_bound_applies_to_same_exact_Jacobian"): True,
        ("checks", "branch_theorem", "all_four_original_angular_weight_covers_contain_actual_exact_unit"): True,
        ("checks", "branch_theorem", "predecessor_original_source_bindings_compared_by_exact_equality"): True,
        ("checks", "branch_theorem", "defining_five_function_zero_identities_not_interval_zero_containment"): True,
        ("checks", "branch_theorem", "determinant_squared_positive_by_transferred_invertibility"): True,
        ("checks", "normalized_conjugacy", "normalized_Neumann_theorem_proves_raw_Jacobian_invertible"): True,
        ("checks", "normalized_conjugacy", "arbitrary_controls_weights_defects_and_positive_units_checked_symbolically"): True,
        ("current_original_true_Rm_leading_implicit_control_C6_functions_installed",): True,
        ("current_original_true_Rm_patch_background_C6_generic_C5_functions_installed",): True,
    }, gates)
    relaxed_value = _receipt(texts[relaxed_check], relaxed_check, {
        ("current_original_Rm_Rh_relaxed_generic_input_certified",): True,
        ("source_ancestor_constructors_called",): False,
        ("source_radius_amplitude_and_width_not_materialized",): True,
        ("p1_p2_whole_path_norm_bounds_certified",): False,
    }, gates)
    bounds_value = _receipt(texts[bounds_check], bounds_check, {
        ("current_original_noncore_generic_quotient_derivative_log_bounds_certified",): True,
        ("original_core_axis_or_loop_denominator_not_promoted",): True,
        ("source_ancestor_constructors_called",): False,
        ("full_O2_buffer_offset_0_11_retained",): True,
    }, gates)
    feedback_value = _receipt(texts[feedback_check], feedback_check, {
        ("current_actual_moment_patch_installed",): True,
        ("current_actual_patch_implicit_axial5_recomputed",): True,
        ("current_five_functional_terminal_identities_connected",): True,
    }, gates)
    # This admission is a structured provenance record: its explanatory
    # fields bind the current C2 source to the correlated C1 class.  It is
    # intentionally not reduced to a scalar truth value.
    feedback_admission = (
        isinstance(feedback_value, Mapping)
        and isinstance(feedback_value.get("original_correlated_source_class_admission"), Mapping)
    )
    feedback_class = feedback_value.get("original_correlated_source_class_admission") if feedback_value else None
    feedback_admission = (
        feedback_admission
        and feedback_class.get("current_C2_bound_controls_original_summed_C1_norm") is True
        and isinstance(feedback_class.get("C2_to_correlated_C1_source_proof"), str)
        and isinstance(feedback_class.get("correlated_functions_not_arbitrary_independent_boxes"), str)
    )
    if not feedback_admission:
        gates.append({"gate": "correlated_source_class_admission_missing", "record": feedback_check})
    repair_value = _json(texts[repair_check], repair_check, gates)
    repair_data_value = _json(texts[repair_data], repair_data, gates)
    bump_value = _json(texts[bump_data], bump_data, gates)
    relaxed_data_value = _json(texts[relaxed_data], relaxed_data, gates)
    bounds_data_value = _json(texts[bounds_data], bounds_data, gates)
    families = [
        c6_value.get("source_family") if c6_value else None,
        relaxed_value.get("source_family") if relaxed_value else None,
        bounds_value.get("source_family") if bounds_value else None,
    ]
    if any(family != SOURCE_FAMILY for family in families):
        gates.append({"gate": "source_family_binding_changed", "expected": SOURCE_FAMILY, "actual": families})
    for value, name in ((feedback_value, feedback_check), (repair_value, repair_check)):
        if value is not None and value.get("actual_five_defect_family_sha256") != SOURCE_FAMILY["actual_five_defect_family_sha256"]:
            gates.append({"gate": "five_defect_family_binding_changed", "record": name})
    old_proof = relaxed_data_value.get("current_whole_Rm_Rh_relaxed_input_and_bounds", {}) if relaxed_data_value else {}
    repair_inverse = repair_data_value.get("actual_whole_axis_coefficient_inverse", {}) if repair_data_value else {}
    theorem = bounds_data_value.get("positive_noncore_source_denominator_theorem", {}) if bounds_data_value else {}
    chart = theorem.get("source_charts", {}).get("actual_patch", {}) if isinstance(theorem, Mapping) else {}
    theorem_ok = all((
        old_proof.get("physical_domain") == "Rm<=R<=Rh=e*Rm=exp(-5)*Rref, x=R/Rm in[1,e], all Z[-1,1]",
        old_proof.get("source_charts") == ["actual_patch"],
        old_proof.get("original_angular_a_bounds", [{}])[0].get("lower") == "0.7",
        old_proof.get("original_angular_a_bounds", [{}, {}])[1].get("upper") == "0.9",
        theorem.get("minimum_a_from_current_original_whole_middle_theorems") == ".7",
        theorem.get("restoration_and_patch_proof") == "same checked actual original normalized swirl lower bounds on whole source domains",
        theorem.get("log_bounds_not_point_profile_definitions") is True,
        chart.get("source_function_positivity_not_inferred_from_saved_box") is True,
        repair_data_value.get("actual_five_defect_family_sha256") == SOURCE_FAMILY["actual_five_defect_family_sha256"] if repair_data_value else False,
        repair_data_value.get("implicit_source_sha256") == SOURCE_FAMILY["implicit_source_sha256"] if repair_data_value else False,
        repair_data_value.get("datum_enclosure_sha256") == SOURCE_FAMILY["datum_enclosure_sha256"] if repair_data_value else False,
        repair_data_value.get("real_axial_domain") == ["-1", "1"] if repair_data_value else False,
        repair_data_value.get("available_interval_x") == ["1", "e"] if repair_data_value else False,
        repair_inverse.get("certified") is True,
        repair_inverse.get("uniform_implicit_C1_family_exists") is True,
        repair_inverse.get("unique_in_initial_box") is True,
    ))
    if not theorem_ok:
        gates.append({"gate": "old_whole_domain_a_theorem_binding_changed", "required": "actual_patch source theorem with minimum .7"})
    if bump_value is None or bump_value.get("fixed_bump_radius") != "1/40" or bump_value.get("centers") != ["5/4", "3/2", "7/4"]:
        gates.append({"gate": "canonical_normalized_bumps_binding_changed"})
    if repair_value is None or repair_value.get("actual_implicit_functional_five_moment_closure_independently_checked") is not True:
        gates.append({"gate": "correlated_five_moment_source_receipt_changed"})
    ok = token_ok and all(value is not None for value in (c6_value, relaxed_value, bounds_value, feedback_value, repair_value,
                                                          repair_data_value, bump_value, relaxed_data_value, bounds_data_value)) and feedback_admission and theorem_ok
    return ok, {
        "records": dict(RECORD_HASHES),
        "source_family": deepcopy(SOURCE_FAMILY),
        "declarations": {name: sorted(found) for name, found in declarations.items()},
        "accepted_receipts": {
            c6_check: {"sha256": RECORD_HASHES[c6_check], "passed": c6_value is not None},
            relaxed_check: {"sha256": RECORD_HASHES[relaxed_check], "passed": relaxed_value is not None},
            bounds_check: {"sha256": RECORD_HASHES[bounds_check], "passed": bounds_value is not None},
            feedback_check: {"sha256": RECORD_HASHES[feedback_check], "passed": feedback_value is not None},
        },
        "theorem_binding": {
            "old_relaxed_source_receipt": relaxed_check,
            "old_relaxed_theorem_data": relaxed_data,
            "generic_source_chart_data": bounds_data,
            "actual_patch_source_chart": "actual_patch",
            "minimum_a": "7/10" if theorem_ok else None,
            "whole_domain": ["1", "exp(1)"],
            "Z_domain": ["-1", "1"],
            "source_function_not_saved_box": theorem_ok,
            "correlated_C2_to_C1_source_class": feedback_admission,
            "same_C6_exact_branch_transfer": c6_value is not None,
            "same_exact_weights_and_normalized_units": c6_value is not None,
            "same_repair_source_domain": repair_data_value is not None and repair_data_value.get("available_interval_x") == ["1", "e"],
        },
        "constructor_binding": {
            "C6_patch_fields": "E=am*H; E_y=x*am*Hx",
            "H": "x^(1/10)+sum_i h[i+2]*gamma_i",
            "Hx": "(1/10)x^(-9/10)+sum_i h[i+2]*gamma_i'",
            "normalized_bumps": "raw beta, radius 1/40, centers 5/4, 3/2, 7/4, one exact normalization",
            "branch": "same unique normalized common-box branch and exact diagonal conjugacy",
        },
    }


def _graph_witnesses(nodes: Sequence[Any], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
    try:
        rows = [
            _expect_row(nodes, Z_NODE, "bound_variable", name="Z"),
            _expect_row(nodes, CHART_NODE, "bound_variable", name="true_Rm_patch_x"),
            _expect_row(nodes, E_NODE, "product", arguments=[183702, 194763]),
            _expect_row(nodes, 194763, "sum", arguments=[186753, 194742, 194749, 194756]),
            _expect_row(nodes, 195164, "product", arguments=[195130, CHART_NODE]),
            _expect_row(nodes, 195130, "product", arguments=[183702, 194812]),
            _expect_row(nodes, 194812, "sum", arguments=[186805, 194791, 194798, 194805]),
            _expect_row(nodes, 48, "exact_rational", numerator=2, denominator=1),
            _expect_row(nodes, 196047, "product", arguments=[E_Y_NODE, 48]),
            _expect_row(nodes, 196053, "negative", argument=196047),
            _expect_row(nodes, 196059, "sum", arguments=[E_NODE, 196053]),
            _expect_row(nodes, A_NODE, "positive_quotient", numerator=196059, denominator=E_NODE,
                        source_positive_certificate="same strictly positive actual prescribed swirl E"),
            _expect_row(nodes, RADIUS_NODE, "product", arguments=[85743, CHART_NODE]),
            _expect_row(nodes, PHYSICAL_JACOBIAN_NODE, "positive_quotient", numerator=1, denominator=CHART_NODE,
                        source_positive_certificate="positive original actual patch x"),
        ]
        return True, {
            "Z": {"handle": Z_NODE, "domain": ["-1", "1"]},
            "chart": {"handle": CHART_NODE, "name": "true_Rm_patch_x", "domain": ["1", "exp(1)"]},
            "ports": {"E": E_NODE, "E_y": E_Y_NODE, "a": A_NODE},
            "constructors": rows,
            "exact_quotient": {
                "numerator": "E - 2 E_y",
                "denominator": "E",
                "graph_factor_handle": 48,
                "graph_factor_exact": "2",
                "formula": "a=(E-2E_y)/E",
            },
            "denominator_source_guard": "same strictly positive actual prescribed swirl E",
            "physical_measure": PHYSICAL_JACOBIAN_NODE,
        }
    except (TypeError, ValueError) as exc:
        gates.append({"gate": "saved_patch_a_constructor_changed", "error": str(exc)})
        return False, {}


def derive(nodes: Sequence[Any], source_records: Mapping[str, Any], *, expected_nodes_digest: str) -> dict[str, Any]:
    """Return a fail-closed source certificate for the actual patch ``a`` root."""
    gates: list[dict[str, Any]] = []
    if expected_nodes_digest != EXPECTED_NODES_DIGEST:
        gates.append({"gate": "unexpected_source_graph_digest_assertion", "expected": EXPECTED_NODES_DIGEST,
                      "actual": expected_nodes_digest})
    try:
        graph_digest = canonical_digest(nodes)
    except (TypeError, ValueError, OverflowError) as exc:
        graph_digest = None
        gates.append({"gate": "source_graph_digest_unavailable", "error": str(exc)})
    if len(nodes) != EXPECTED_NODE_COUNT or graph_digest != EXPECTED_NODES_DIGEST:
        gates.append({"gate": "saved_source_graph_changed", "expected_count": EXPECTED_NODE_COUNT,
                      "actual_count": len(nodes), "expected_digest": EXPECTED_NODES_DIGEST,
                      "actual_digest": graph_digest})
    source_ok, source_witnesses = _source_witnesses(source_records, gates)
    graph_ok, graph_witnesses = _graph_witnesses(nodes, gates) if graph_digest == EXPECTED_NODES_DIGEST else (False, {})
    admitted = bool(source_ok and graph_ok and not gates)
    return {
        "schema": SCHEMA,
        "source_only": True,
        "nodes_digest": EXPECTED_NODES_DIGEST,
        "passed": admitted,
        "actual_target_handle": A_NODE,
        "E_handle": E_NODE,
        "E_y_handle": E_Y_NODE,
        "chart_handle": CHART_NODE,
        "strict_positive_proved": admitted,
        "whole_native_domain": ["1", "exp(1)"],
        "whole_Z_domain": ["-1", "1"],
        "lower_bound_exact": "7/10" if admitted else None,
        "witnesses": {"source": source_witnesses, "saved_graph": graph_witnesses},
        "unresolved_gates": gates,
        "old_operators_replayed": False,
        "operators_reexecuted": False,
        "numeric_values_computed": False,
        "full_goal_complete": False,
        "source_record_hashes": dict(RECORD_HASHES),
    }


__all__ = ["RECORD_HASHES", "EXPECTED_NODES_DIGEST", "canonical_digest", "derive"]
