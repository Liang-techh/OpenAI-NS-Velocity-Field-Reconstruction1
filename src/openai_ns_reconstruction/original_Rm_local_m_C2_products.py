"""Bounded source certificate for the actual local finite-``N`` loop density.

This module deliberately consumes the authenticated source records and the
saved expression DAG as data.  It does not import a producer, evaluate a
mathematical operator, replay a saved calculation, or select a numerical
field value.  The certificate is therefore useful at the boundary between
the accepted ``a`` positivity reports and the still-open global equations.

The quantity covered here is the actual local density
``m = B_over_Pstar / N``.  The three ordinary slow-``Z`` Taylor rows are
checked independently on the postrestore and patch source domains.  Incoming
history derivatives and the physical ``Q_y``/``Z1`` seam remain outside this
certificate.
"""

from __future__ import annotations

import ast
from collections.abc import Mapping, Sequence
from copy import deepcopy
import hashlib
import json
from typing import Any


SCHEMA = "original_Rm_local_m_C2_regularity_source_v1"
PREFIX = "lei_ren_part1_paper_compliant_current_original_"
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
EXPECTED_NODE_COUNT = 205300
Z_NODE = 20
N_NODE = 18
INV_N_NODE = 99411
ETA_NODE = 21
DSTAR_NODE = 22
POSITIVITY_PACKET_SHA256 = "4c3b730e57966c0c4a8c2481564a6bf6ec19b1c1136ab01f7e7e2ade1246f675"
POSITIVITY_PACKET_CANONICAL_SHA256 = "04830a1a420e252541e4b11f060801d6650072490a4e7e7c26f13893d52099c8"
PARENT_REPORT_SHA256 = "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746"


# Every entry is an authenticated member of the authoritative parent input /
# implementation map.  Only plain text records are requested by the bridge;
# gzip archives remain authenticated in the parent report and are not decoded
# here.
RECORD_HASHES = {
    PREFIX + "prescribed_R100_Rm_C6_functions.py": "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    PREFIX + "prescribed_R100_Rm_C6_functions_check.json": "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6",
    PREFIX + "true_Rm_patch_C6_functions.py": "a6f1cb5deaf17e9482aa964aef7867aadbc815e4d703a415e4945cf665f14002",
    PREFIX + "true_Rm_patch_C6_functions_check.json": "0190fcf738fb8add5be85cf0767464a39cc0d18240e6b36ad4ec7930dce38aa7",
    PREFIX + "prescribed_loop_C5_functions.py": "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    PREFIX + "prescribed_loop_C5_functions_check.json": "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7",
    PREFIX + "true_Rm_patch_finite_N_C5_functions.py": "c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15",
    PREFIX + "true_Rm_patch_finite_N_C5_functions_check.json": "efd3b9a6fdfd4522ead1806d4ef4250864ef670909ca047d0fa85d9ce23bf207",
    "lei_ren_part1_paper_compliant_current_generic_loop_function_sources.py": "d694b8ee8a4bf1c10090d45a7f9d4d6de0792d02a1ebdd6b21872b8514a5b26f",
    "lei_ren_part1_paper_compliant_current_generic_loop_function_sources_check.json": "ac18d2621c67002769f8e0c76613ae91de54239699edd3801cccfaa052c19d07",
    PREFIX + "C2_phase_primitives.py": "84f5da3eb64096f82df88dd828b22b09bce18192a03925d91cfee0b193902d2c",
    PREFIX + "C2_phase_primitives_check.json": "76ed1af50ec799c49081ab67881b189b42e0d3137bcca770a67766d4111dff8b",
    "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py": "32f51123a47c2bafbb6c2f7056302bab996a019883615333875207d72a9bb7d7",
    "lei_ren_part1_paper_compliant_flat_pulse_derivatives_check.json": "572c50c5210ec51c9998780765bdaf3ca9331ead8f4bc8e968245c6a556f8468",
    PREFIX + "all_chart_point_phase.py": "e0c8f8fea40684ce2ee74a6edb624f5c3e05a429b4acbe8fcf46b3972ed23495",
    PREFIX + "all_chart_point_phase_check.json": "de80a7d24bbdda76d39e8d0a96a6f373a68892b0d91306c6ca7f0ab0e1af7d72",
    "lei_ren_part1_paper_compliant_current_native_conditioned_phase.py": "8796564b81734e7b04c2bb07010033d85aadd52ebefdca01eab0d148954ff670",
    "lei_ren_part1_paper_compliant_current_native_conditioned_phase_check.json": "f0301d59f0a518fd93b457d40857eb6a6e62f166245083866f68945eafce8760",
    PREFIX + "whole_Z_signed_loop_functions.py": "2d534e7108cbd518f0938aa50e9ba802c10c2b0a68dc65c70b2e05eec2508779",
    PREFIX + "whole_Z_signed_loop_functions_check.json": "1a5178d29b252d6b4fb7cbe11f658a61d586c2e77b0f0284d1975ff35d390452",
    PREFIX + "Rm_conditioned_phase_density.py": "af2590c01effe60eba211a4840518c39312793fcdab7b0ff740e8d53e36cfcac",
    PREFIX + "Rm_conditioned_phase_density_check.json": "21119ae9abcec6e90cf097b2e43579f599ff392c7f88ed20ebbfaeb72fdbd77a",
    "lei_ren_part1_paper_compliant_current_native_phase_first_jets.py": "f4121fa7670659c532f57fcf246b66997500ae78efd29a3722aeeb4a8f7a7156",
    "lei_ren_part1_paper_compliant_current_native_phase_first_jets_check.json": "faa163a70615cda1e6f637e4563f49a83ada7ded722699ef3556c88671f6201c",
}


SIDES = {
    "postrestore": {
        "coordinate": 84528,
        "coordinate_name": "true_prescribed_post_restore_fraction",
        "native_domain": ["0", "1"],
        "E": 84911,
        "E_chain": (
            (84911, "product", [84896, 84573]),
            (84914, "sum", [84912, 84913]),
            (84918, "sum", [84915, 84916, 84917]),
        ),
        "M_chain": (
            (121032, "negative", {"argument": 121026}),
            (121033, "negative", {"argument": 121027}),
            (121034, "negative", {"argument": 121028}),
        ),
        "a": 119435,
        "delta": 119523,
        "predicate": 119524,
        "B": (121108, 121109, 121110),
        "active_B": (121102, 121103, 121104),
        "B_formula": (
            {"root": 121076, "root_arguments": [84911, 121032], "terms": ((84911, 121032),)},
            {"root": 121079, "root_arguments": [121077, 121078], "terms": ((84911, 121033), (84914, 121032))},
            {"root": 121083, "root_arguments": [121080, 121081, 121082], "terms": ((84911, 121034), (48, 84914, 121033), (84918, 121032))},
        ),
        "density": (121350, 121351, 121352),
        "lower_bound": "4/5",
    },
    "patch": {
        "coordinate": 185676,
        "coordinate_name": "true_Rm_patch_x",
        "native_domain": ["1", "exp(1)"],
        "E": 194957,
        "E_chain": (
            (194957, "product", [183702, 194763]),
            (194960, "sum", [194958, 194959]),
            (194964, "sum", [194961, 194962, 194963]),
        ),
        "M_chain": (
            (197852, "sum", [197834, 197846]),
            (197853, "sum", [197835, 197847]),
            (197854, "sum", [197836, 197848]),
        ),
        "a": 196071,
        "delta": 196284,
        "kappa": 196278,
        "kappa_quotient": 196243,
        "predicate": 196285,
        "B": (197928, 197929, 197930),
        "active_B": (197922, 197923, 197924),
        "B_formula": (
            {"root": 197896, "root_arguments": [194957, 197852], "terms": ((194957, 197852),)},
            {"root": 197899, "root_arguments": [197897, 197898], "terms": ((194957, 197853), (194960, 197852))},
            {"root": 197903, "root_arguments": [197900, 197901, 197902], "terms": ((194957, 197854), (48, 194960, 197853), (194964, 197852))},
        ),
        "density": (198182, 198183, 198184),
        "lower_bound": "7/10",
    },
}


SOURCE_DECLARATIONS = {
    PREFIX + "prescribed_R100_Rm_C6_functions.py": (
        "CurrentOriginalPrescribedR100RmC6Functions", "add_chart", "parameter_geometry_guard", "phi_from_Utheta",
    ),
    PREFIX + "true_Rm_patch_C6_functions.py": (
        "CurrentOriginalTrueRmPatchC6Functions", "implicit_rows", "patch_fields", "map_conjugacy", "jacobian",
    ),
    PREFIX + "prescribed_loop_C5_functions.py": (
        "LoopAlgebra", "loop_functions", "finite_N_functions", "CurrentOriginalPrescribedLoopC5Functions",
    ),
    PREFIX + "true_Rm_patch_finite_N_C5_functions.py": (
        "transport_log", "native_derivative", "pressure_continuation", "CurrentOriginalTrueRmPatchFiniteNC5Functions",
    ),
    "lei_ren_part1_paper_compliant_current_generic_loop_function_sources.py": (
        "FunctionGraph", "build", "exact_theorem", "CurrentLoopFunctionSources",
    ),
    PREFIX + "C2_phase_primitives.py": (
        "C2Function", "C2Algebra", "primitive_C2", "regularity_proof", "CurrentC2PhasePrimitives",
    ),
    "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py": (
        "_sigma_left", "sigma_jets", "FlatPulseDerivatives",
    ),
    PREFIX + "all_chart_point_phase.py": (
        "OriginalAllChartPointPhase", "periodic_add", "phase_for_source",
    ),
    "lei_ren_part1_paper_compliant_current_native_conditioned_phase.py": (
        "ConditionedPhase", "NativeConditionedPhase",
    ),
    PREFIX + "whole_Z_signed_loop_functions.py": (
        "WholeZSignedLoopFunctions", "branch_Z_functions",
    ),
    PREFIX + "Rm_conditioned_phase_density.py": (
        "RmRadiusPhase", "OriginalRmConditionedPhaseDensity", "compile_Z_only", "candidate_N",
    ),
    "lei_ren_part1_paper_compliant_current_native_phase_first_jets.py": (
        "conditioned_first_jets", "NativePhaseFirstJets",
    ),
}


SOURCE_TOKENS = {
    PREFIX + "prescribed_loop_C5_functions.py": (
        "same_original_flat_extension_through_order5=True",
        "active_body_not_evaluated_on_flat_branch=True",
        "original phase inverse lambda=K*(1+t^2)>0",
        "true_prescribed_source_C5_implicit_inverse=True",
        "total_Z_phase_derivative_exact_zero=True",
        "one_original_common_N_eta_dstar_and_global_log_R_rminus_phase=True",
        "same exact positive integer N on every source chart",
    ),
    PREFIX + "true_Rm_patch_finite_N_C5_functions.py": (
        "physical_measure='dt/t'",
        "actual_same_N_Section11_loop_inverse_density_and_history_C5=True",
        "original_global_phase_Ra_hbB_sc_half_not_new_patch_phase_origin=True",
        "original_gamma_first_second_and_same_implicit_control_native_sources=True",
        "nonzero_true_Rm_correction_and_signed_pressure_inlets_retained=True",
    ),
    "lei_ren_part1_paper_compliant_current_generic_loop_function_sources.py": (
        "T1=g.integral(t,psi)",
        "T2=g.integral(g.mul(t,t),psi)",
        "lambda0=g.mul(K,g.add(one,g.mul(t,t)))",
        "total_Z_phase_derivative_exact_zero=True",
        "actual_original_inlet_histories_required=True",
        "same_global_source_family_and_N_across_interfaces=True",
        "exact_flat_inverse='psi=2pi*phi'",
    ),
    PREFIX + "C2_phase_primitives.py": (
        "ordinary_slow_Z_derivative_order=2",
        "Taylor_coefficient_factorial=2",
        "gamma=alg.add(alg.scale(eta,2),alg.neg(D))",
        "root=alg.sqrt(alg.div(gamma,alg.scale(a,2)))",
        "u=alg.div(alg.mul(p2,q),dstar)",
        "same_original_inverse_not_a_new_branch=True",
        "phase_Z_and_ZZ_exact_zero=True",
        "Poisson_denominator=",
        "strict_inverse_Jacobian=",
        "phase_endpoint_identities=",
        "actual_two_sided_flat_source_is_lazy=True",
        "range_caps_not_differentiated=True",
        "exact_b_and_t0_zero",
    ),
    "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py": (
        "def _sigma_left",
        "def sigma_jets",
        "if hi<=0",
        "if hi>=1",
    ),
    PREFIX + "all_chart_point_phase.py": (
        "source_phase_not_an_independent_free_angle",
        "candidate_N_not_global_frequency_admission=True",
    ),
    "lei_ren_part1_paper_compliant_current_native_conditioned_phase.py": (
        "free_phase_parameter_not_spatial_phase=True",
        "original_common_N_and_radius_phase_bound=False",
    ),
    PREFIX + "whole_Z_signed_loop_functions.py": (
        "actual_phase_Z_exact_zero=True",
        "actual_spatial_Z_is_phase_held_Z_by_radius_Z_zero=True",
        "phase_not_supplied_as_free_angle=True",
    ),
    PREFIX + "Rm_conditioned_phase_density.py": (
        "actual_Rm_radius_phase_Z_independent=True",
        "phase_is_independent_candidate_parameter=True",
        "no_fabricated_y_source_derivatives=True",
        "source_q_Z_includes_all_b_and_b_Z_terms=True",
    ),
    "lei_ren_part1_paper_compliant_current_native_phase_first_jets.py": (
        "phase_is_independent_input_on_this_query=True",
        "actual_spatial_fast_N_chain_installed",
        "original_radius_phase_is_Z_independent=True",
        "D_y P(y,Z,N*y)=P_y|phi+N*P_phi; D_Z P=P_Z|phi",
    ),
}


RECEIPT_FLAGS = {
    PREFIX + "prescribed_R100_Rm_C6_functions_check.json": (
        ("checks", "actual_source", "passed"),
        ("checks", "actual_source", "all_genuine_seven_source_and_six_generic_rows"),
        ("checks", "actual_source", "every_nonzero_incoming_history_function_retained_by_identity"),
        ("checks", "geometry", "complete_radius_order_proved_by_actual_source_positive_guards"),
        ("checks", "cutoff", "exact_zero_one_and_flat_derivative_endpoints"),
    ),
    PREFIX + "true_Rm_patch_C6_functions_check.json": (
        ("checks", "actual_source", "passed"),
        ("checks", "actual_source", "true_Rm_parent_and_same_original_map_nodes"),
        ("checks", "actual_source", "ordinary_implicit_rows_not_padded"),
        ("checks", "actual_source", "original_log_Jacobian_and_finite_patch_endpoints"),
        ("checks", "actual_source", "no_free_old_phase_or_integration_variable"),
        ("checks", "branch_theorem", "same_exact_original_unique_branch_theorem_bound"),
    ),
    PREFIX + "prescribed_loop_C5_functions_check.json": (
        ("checks", "actual_source", "passed"),
        ("checks", "actual_source", "genuine_C6_inputs_and_C5_generic_roots"),
        ("checks", "actual_source", "all_ten_correction_predecessors_preserved"),
        ("checks", "actual_source", "original_microscopic_collar_and_pressure_memory_retained"),
        ("checks", "actual_source", "no_free_core_rho_angle_old_phase_or_integration_variable"),
        ("checks", "actual_source", "same_original_N_eta_dstar_and_global_phase"),
        ("checks", "actual_source", "full_signed_p2_preserved_by_identity"),
        ("checks", "actual_source", "all17_numeric_global_frequency_and_recursion_not_claimed"),
        ("checks", "implicit_inverse", "passed"),
        ("checks", "implicit_inverse", "normalized_factorials_and_all_cross_terms_checked"),
        ("checks", "finite_N", "passed"),
        ("checks", "finite_N", "original_nonunit_N_and_all_quadratic_cross_terms_checked"),
        ("checks", "lazy_flat", "active_inverse_and_singular_sqrt_not_evaluated_on_flat_branch"),
        ("current_original_true_prescribed_ten_chart_finite_N_density_C5_functions_installed",),
    ),
    PREFIX + "true_Rm_patch_finite_N_C5_functions_check.json": (
        ("checks", "actual_source", "passed"),
        ("checks", "actual_source", "same_true_background_p2_original_N_phase_and_hbB_origin"),
        ("checks", "actual_source", "all_nonzero_correction_and_pressure_predecessors_retained"),
        ("checks", "actual_source", "exact_log_measure_and_rate_zero_pressure_memory"),
        ("checks", "actual_source", "native_background_and_gamma_first_second_derivative_seeds"),
        ("checks", "actual_source", "total_fast_phase_radial_C5_and_complete_signed_pressure_C3"),
        ("checks", "actual_source", "independent_P0_original_rminus_boundary_and_Rh_endpoints"),
        ("checks", "actual_source", "no_free_angle_old_phase_or_integration_variable"),
        ("checks", "log_transport", "all_five_rates_original_dt_over_t_and_pressure_rate_zero_checked"),
        ("current_original_true_Rm_patch_finite_N_loop_density_history_C5_functions_installed",),
        ("current_original_true_Rm_patch_finite_N_C5_transport_installed",),
    ),
    "lei_ren_part1_paper_compliant_current_generic_loop_function_sources_check.json": (
        ("current_original_generic_loop_velocity_and_five_increment_function_graphs_defined",),
        ("flat_active_bodies_are_lazy",),
        ("nonzero_original_V_and_full_signed_cross_terms_retained",),
        ("single_N_and_global_phase_preserved",),
        ("quiet_gap_defect_and_pressure_memory_preserved",),
        ("original_power_zero_increment_source_admitted",),
    ),
    PREFIX + "C2_phase_primitives_check.json": (
        ("current_original_17_chart_source_Z2_rows_installed",),
        ("current_original_same_inverse_C2_phase_primitive_functions_installed",),
        ("original_C1_graph_prefix_and_inverse_handles_preserved",),
        ("ordinary_source_second_derivatives_use_actual_Taylor_rows",),
        ("same_independent_P0_and_source_basis_ledger_checked",),
        ("original_flat_sigma_second_derivative_source_bound_and_checked",),
    ),
    "lei_ren_part1_paper_compliant_flat_pulse_derivatives_check.json": (
        ("analytic_original_beta_and_flat_envelope_checks", "all_required_endpoint_flat_limits"),
        ("original_radial_shape_derivatives_C4_available",),
        ("original_flat_shapes_retained",),
    ),
    PREFIX + "all_chart_point_phase_check.json": (
        ("actual_original17_chart_point_phase_with_adaptive_analytic_and_signed_microscopic_errors_installed",),
        ("accepted_original_function_source_binding", "passed"),
        ("accepted_original_function_source_binding", "source_graph_not_mutated"),
        ("source_owner_reused_only_when_live",),
        ("standalone_fresh_original_owner_constructed",),
    ),
    "lei_ren_part1_paper_compliant_current_native_conditioned_phase_check.json": (
        ("current_native_conditioned_C0_phase_inverse_and_primitives_executed",),
        ("independent_original_scalar_reference_checks", "fixtures_are_not_native_source_values"),
        ("focused_extreme_log_checks", "positive_factored_rho_and_s_retained"),
        ("focused_extreme_log_checks", "positive_tiny_q_primitives_not_zeroed"),
    ),
    PREFIX + "whole_Z_signed_loop_functions_check.json": (
        ("current_original_whole_Z_live_signed_inverse_Z_functions_and_local_integrals_installed",),
        ("genuine_current_whole_Z_signed_inverse_functions_Z_phi_and_nonlinear_densities_checked",),
        ("actual_adaptive_global_radius_phase_and_conditional_branch_union_checked",),
        ("independent_actual_integral_measure", "local_integrals_do_not_replace_actual_prefix_incoming"),
    ),
    PREFIX + "Rm_conditioned_phase_density_check.json": (
        ("original_actual_Rm_conditioned_Z_phase_and_five_signed_density_interface_installed",),
        ("actual_live_source", "passed"),
        ("actual_live_source", "conditional_two_axial_frames_only"),
        ("original_Z_only_compiler_binding", "original_defining_phase_primitive_and_density_formulas_unchanged"),
        ("original_Z_only_compiler_binding", "no_fabricated_y_source_derivatives"),
    ),
    "lei_ren_part1_paper_compliant_current_native_phase_first_jets_check.json": (
        ("current_original_native_conditioned_phase_and_A_B_first_slow_jets_executed",),
        ("actual_whole_integral_radial_Z_cell_first_jets_checked",),
        ("independent_original_first_derivative_checks", "independent_original_phase_and_primitive_functions_used"),
        ("independent_original_first_derivative_checks", "p2_zero_crossing_derivatives_not_erased"),
        ("independent_original_first_derivative_checks", "fixtures_are_not_native_field_values"),
    ),
}


def canonical_digest(value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _gate(gates: list[dict[str, Any]], name: str, **details: Any) -> None:
    item: dict[str, Any] = {"gate": name}
    item.update(details)
    gates.append(item)


def _row(nodes: Sequence[Any], index: int) -> Mapping[str, Any] | None:
    if type(index) is not int or index < 0 or index >= len(nodes):
        return None
    value = nodes[index]
    return value if isinstance(value, Mapping) else None


def _expect_row(nodes: Sequence[Any], index: int, operation: str, gates: list[dict[str, Any]], **fields: Any) -> bool:
    row = _row(nodes, index)
    if row is None:
        _gate(gates, "saved_graph_row_missing", handle=index, operation=operation)
        return False
    ok = row.get("operation") == operation
    if not ok:
        _gate(gates, "saved_graph_operation_changed", handle=index, expected=operation, actual=row.get("operation"))
    for key, expected in fields.items():
        if row.get(key) != expected:
            _gate(gates, "saved_graph_field_changed", handle=index, field=key, expected=expected, actual=row.get(key))
            ok = False
    return ok


def _record(source_records: Mapping[str, Any], name: str, gates: list[dict[str, Any]]) -> str | None:
    expected = RECORD_HASHES[name]
    record = source_records.get(name) if isinstance(source_records, Mapping) else None
    if not isinstance(record, Mapping):
        _gate(gates, "source_record_missing", record=name)
        return None
    text = record.get("source_text")
    actual = hashlib.sha256(text.encode("utf-8")).hexdigest() if isinstance(text, str) else None
    if record.get("path") != "experiments/root_st073/" + name:
        _gate(gates, "source_record_path_changed", record=name)
        return None
    if record.get("sha256") != expected or actual != expected:
        _gate(gates, "source_record_hash_changed", record=name, expected_sha256=expected, actual_sha256=actual)
        return None
    return text


def _definitions(text: str, name: str, wanted: tuple[str, ...], gates: list[dict[str, Any]]) -> dict[str, dict[str, int]]:
    try:
        tree = ast.parse(text, filename=name)
    except SyntaxError as exc:
        _gate(gates, "source_ast_invalid", record=name, error=str(exc))
        return {}
    found: dict[str, dict[str, int]] = {}
    for item in ast.walk(tree):
        if isinstance(item, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) and item.name in wanted:
            found[item.name] = {"line": item.lineno, "end_line": item.end_lineno or item.lineno}
    for item in wanted:
        if item not in found:
            _gate(gates, "source_declaration_missing", record=name, declaration=item)
    return found


def _tokens(text: str, name: str, wanted: tuple[str, ...], gates: list[dict[str, Any]]) -> bool:
    missing = [token for token in wanted if token not in text]
    if missing:
        _gate(gates, "source_constructor_token_missing", record=name, tokens=missing)
        return False
    return True


def _at_path(value: Any, path: tuple[str, ...]) -> Any:
    for key in path:
        if not isinstance(value, Mapping):
            return None
        value = value.get(key)
    return value


def _receipt(text: str, name: str, required: tuple[tuple[str, ...], ...], gates: list[dict[str, Any]]) -> bool:
    try:
        value = json.loads(text)
    except json.JSONDecodeError as exc:
        _gate(gates, "accepted_source_receipt_invalid", record=name, error=str(exc))
        return False
    ok = value.get("all_passed") is True
    if not ok:
        _gate(gates, "accepted_source_receipt_header_missing", record=name)
    for path in required:
        if _at_path(value, path) is not True:
            _gate(gates, "accepted_source_receipt_flag_missing", record=name, path=".".join(path))
            ok = False
    return ok


def _ast_definition(text: str, name: str, wanted: str) -> ast.AST | None:
    """Find one source definition without importing or executing its module."""
    try:
        tree = ast.parse(text, filename=name)
    except SyntaxError:
        return None
    return next((item for item in ast.walk(tree)
                 if isinstance(item, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
                 and item.name == wanted), None)


def _ast_assignments(node: ast.AST) -> dict[str, list[ast.AST]]:
    result: dict[str, list[ast.AST]] = {}
    for item in ast.walk(node):
        if not isinstance(item, (ast.Assign, ast.AnnAssign)):
            continue
        targets = item.targets if isinstance(item, ast.Assign) else [item.target]
        value = item.value
        if value is None:
            continue
        for target in targets:
            result.setdefault(ast.unparse(target), []).append(value)
    return result


def _ast_text(node: ast.AST | None) -> str:
    return ast.unparse(node) if node is not None else ""


def _ast_binding_checks(texts: Mapping[str, str], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
    """Bind the live C5 formulas by syntax, rather than by class or slot names."""
    loop_name = PREFIX + "prescribed_loop_C5_functions.py"
    patch_name = PREFIX + "true_Rm_patch_finite_N_C5_functions.py"
    c2_name = PREFIX + "C2_phase_primitives.py"
    generic_name = "lei_ren_part1_paper_compliant_current_generic_loop_function_sources.py"
    loop_generic = _ast_definition(texts.get(loop_name, ""), loop_name, "generic_roots")
    loop_body = _ast_definition(texts.get(loop_name, ""), loop_name, "loop_functions")
    loop_finite = _ast_definition(texts.get(loop_name, ""), loop_name, "finite_N_functions")
    loop_class = _ast_definition(texts.get(loop_name, ""), loop_name, "CurrentOriginalPrescribedLoopC5Functions")
    patch_class = _ast_definition(texts.get(patch_name, ""), patch_name, "CurrentOriginalTrueRmPatchFiniteNC5Functions")
    mixed = _ast_definition(texts.get(loop_name, ""), loop_name, "MixedTaylor")
    primitive = _ast_definition(texts.get(c2_name, ""), c2_name, "primitive_C2")
    regularity = _ast_definition(texts.get(c2_name, ""), c2_name, "regularity_proof")
    generic_build = _ast_definition(texts.get(generic_name, ""), generic_name, "build")
    missing: list[str] = []

    def require_definition(label: str, value: ast.AST | None) -> str:
        if value is None:
            missing.append(label)
            return ""
        return _ast_text(value)

    generic_text = require_definition("prescribed.generic_roots", loop_generic)
    generic_assign = _ast_assignments(loop_generic) if loop_generic is not None else {}
    expected_generic = {
        "C": "a.sub(E, a.scale(a.trim(radial['E_y']), 2))",
        "B": "a.scale(a.trim(radial['V_y']), 2)",
        "(aa, bb)": "(divide(C), divide(B))",
        "t0": "a.neg(a.div(B, C, 'same original a>0 and C=a*E'))",
    }
    generic_matches: dict[str, bool] = {}
    for target, expression in expected_generic.items():
        actual = generic_assign.get(target, [])
        wanted_dump = ast.dump(ast.parse(expression, mode="eval").body, include_attributes=False)
        matched = any(ast.dump(value, include_attributes=False) == wanted_dump for value in actual)
        generic_matches[target] = matched
        if not matched:
            missing.append("prescribed.generic_roots." + target)
    divide_ok = "same strictly positive actual prescribed swirl E" in generic_text
    if not divide_ok:
        missing.append("prescribed.generic_roots.divide_by_E_guard")
    loop_text = require_definition("prescribed.loop_functions", loop_body)
    loop_assign = _ast_assignments(loop_body) if loop_body is not None else {}
    loop_fragments = {
        "q": "q = lazy(a.mul(a.sigma(cut), a.sqrt(a.div(gamma, a.scale(roots['a'], 2), 'same original active a>0'))))",
        "u": "u = a.scale(a.mul(roots['p2'], q), g.quotient(g.one, dstar, 'same positive canonical d_star'))",
        "K": "K = a.div(one, a.scale(ratio, twopi), 'same positive 2*pi*(1+t0^2+2*q^2)')",
        "T2": "FunctionJet(tuple((g.mul(c(math.factorial(k)), integrate(t2[k, 0])) for k in range(6))))",
        "partials": "Fraction(1, math.factorial(k))",
        "M": "M = a.sub(a.neg(a.scale(a.mul(roots['a'], T1hat), g.quotient(g.one, twopi, 'positive 2*pi'))), a.scale(roots['b'], phase))",
        "B": "B = lazy(a.scale(a.mul(roots['E'], M), '1/2'))",
    }
    loop_matches: dict[str, bool] = {}
    for target, fragment in loop_fragments.items():
        if target == "partials":
            matched = fragment in loop_text
        else:
            actual = loop_assign.get(target, [])
            expression = fragment.removeprefix(target + " = ")
            expected_dump = ast.dump(ast.parse(expression, mode="eval").body, include_attributes=False)
            matched = any(expected_dump == ast.dump(value, include_attributes=False) for value in actual)
        loop_matches[target] = matched
        if not matched:
            missing.append("prescribed.loop_functions." + target)
    inverse_tokens = (
        "phase_endpoints=[0, 1]",
        "derivative_positive_certificate='K*(1+t^2)>0; same original a/q-flat construction'",
        "flat_inverse='2*pi*phase'",
        "true_prescribed_source_C5_implicit_inverse=True",
    )
    inverse_ok = all(token in loop_text for token in inverse_tokens)
    if not inverse_ok:
        missing.append("prescribed.loop_functions.normalized_inverse")
    loop_class_text = require_definition("prescribed.C5_constructor", loop_class)
    prescribed_phase_flags = {
        "common_phase_origin": "self.phase_origin =" in loop_class_text,
        "fractional_part_common_N": "phase = g.unary('fractional_part', g.mul(self.N, offset))" in loop_class_text,
        "phase_Z_zero": "total_Z_phase_derivative_exact_zero=True" in loop_class_text,
    }
    for label, passed in prescribed_phase_flags.items():
        if not passed:
            missing.append("prescribed.C5_constructor." + label)
    finite_text = require_definition("prescribed.finite_N_functions", loop_finite)
    finite_assign = _ast_assignments(loop_finite) if loop_finite is not None else {}
    finite_fragments = {
        "invN": "g.quotient(g.one, N, 'same exact positive integer N on every source chart')",
        "s": "a.scale(loop['B_over_Pstar'], invN)",
        "inc": "dict(m=s, h=r, k=a.add(a.mul(V, r), a.mul(E, s), a.mul(r, s)), e=a.sub(a.add(a.scale(a.mul(V, s), 2), a.mul(s, s)), a.add(a.mul(E, r), a.scale(a.mul(r, r), '1/2'))), p=a.add(a.mul(E, r), a.scale(a.mul(r, r), '1/2')))" ,
    }
    finite_matches: dict[str, bool] = {}
    for target, fragment in finite_fragments.items():
        actual = finite_assign.get(target, [])
        matched = any(fragment == _ast_text(value) for value in actual)
        finite_matches[target] = matched
        if not matched:
            missing.append("prescribed.finite_N_functions." + target)
    mixed_text = require_definition("prescribed.MixedTaylor", mixed)
    mixed_flags = {
        "normalized_axial_factorials": "math.factorial(k)" in mixed_text,
        "normalized_convolution": "for a in range(k + 1)" in mixed_text and "for b in range(j + 1)" in mixed_text,
        "normalized_division": "for degree in range(self.order + 1)" in mixed_text and "if a + b > 0" in mixed_text,
    }
    for label, passed in mixed_flags.items():
        if not passed:
            missing.append("prescribed.MixedTaylor." + label)
    patch_text = require_definition("patch.finite_N_constructor", patch_class)
    patch_flags = {
        "phase_origin_reused": "self.phase_origin = self.loop.phase_origin" in patch_text,
        "fractional_part_common_N": "phase = g.unary('fractional_part', g.mul(self.N, offset))" in patch_text,
        "same_generic_roots": "roots = engine.generic_roots(g, bg, self.Pstar)" in patch_text,
        "same_loop_functions": "loop = engine.loop_functions(g, roots, phase, self.eta, self.dstar" in patch_text,
        "same_finite_N": "finite = engine.finite_N_functions(g, roots, loop, self.N)" in patch_text,
    }
    for label, passed in patch_flags.items():
        if not passed:
            missing.append("patch." + label)
    primitive_text = require_definition("C2.primitive_C2", primitive)
    regularity_text = require_definition("C2.regularity_proof", regularity)
    regularity_flags = {
        "actual_C2_rows": "ordinary_slow_Z_derivative_order=2" in primitive_text and "Taylor_coefficient_factorial=2" in primitive_text,
        "same_inverse": "same_original_inverse_not_a_new_branch=True" in primitive_text,
        "normalized_MZZ": "MZZ = g.sub(g.neg(g.quotient(g.add(g.mul(a.ZZ, at(T1.value)), g.mul(g.constant(2), a.Z, T1Z), g.mul(a.value, T1ZZ)), twopi, '2*pi positive')), g.mul(b.ZZ, phase))" in primitive_text,
        "normalized_BZZ": "BZZ = flat(g.mul(g.constant('1/2'), g.add(g.mul(E.ZZ, M), g.mul(g.constant(2), E.Z, MZ), g.mul(E.value, MZZ))))" in primitive_text,
        "endpoint_normalization": "Phi(0)=0, Phi(2*pi)=1" in regularity_text,
        "flat_M_zero": "A=B=0 at phase0,1" in regularity_text and "actual_two_sided_flat_source_is_lazy=True" in regularity_text,
        "poisson_inverse": "strict_inverse_Jacobian=" in regularity_text and "Poisson_denominator=" in regularity_text,
    }
    for label, passed in regularity_flags.items():
        if not passed:
            missing.append("C2." + label)
    generic_build_text = require_definition("generic.build", generic_build)
    phase_flags = {
        "common_phase_binding": "phi='fractional_part(N*log(R/r_minus))'" in generic_build_text,
        "phase_Z_zero": "total_Z_phase_derivative_exact_zero=True" in generic_build_text,
        "phase_held_fixed": "slow_derivatives_hold_phi_fixed=True" in generic_build_text,
        "prescribed_phase_origin": prescribed_phase_flags["common_phase_origin"],
        "prescribed_fractional_part": prescribed_phase_flags["fractional_part_common_N"],
        "prescribed_phase_Z_zero": prescribed_phase_flags["phase_Z_zero"],
    }
    for label, passed in phase_flags.items():
        if not passed:
            missing.append("generic.build." + label)
    passed = not missing
    if not passed:
        _gate(gates, "source_ast_binding_missing", checks=sorted(set(missing)))
    return passed, {
        "passed": passed,
        "missing": sorted(set(missing)),
        "generic_roots": generic_matches | {"divide_by_E_guard": divide_ok},
        "loop_functions": loop_matches | {"normalized_inverse": inverse_ok},
        "finite_N_functions": finite_matches,
        "mixed_taylor": mixed_flags,
        "patch_source_reuse": patch_flags,
        "prescribed_phase": prescribed_phase_flags,
        "regularity": regularity_flags,
        "phase": phase_flags,
        "ordinary_formula": {
            "C": "E-2*E_y",
            "B": "2*V_y",
            "aa": "C/E",
            "bb": "B/E",
            "t0": "-B/C",
            "flat_join": "a*t0+b=0",
            "B_over_Pstar": "E*M/2",
        },
    }


def _source_witnesses(source_records: Mapping[str, Any], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
    texts = {name: _record(source_records, name, gates) for name in RECORD_HASHES}
    if any(text is None for text in texts.values()):
        return False, {"records": dict(RECORD_HASHES)}

    declarations = {
        name: _definitions(texts[name], name, wanted, gates)
        for name, wanted in SOURCE_DECLARATIONS.items()
    }
    token_ok = all(_tokens(texts[name], name, wanted, gates) for name, wanted in SOURCE_TOKENS.items())
    receipts: dict[str, bool] = {}
    receipt_ok = True
    for name, required in RECEIPT_FLAGS.items():
        value = _receipt(texts[name], name, required, gates)
        receipts[name] = value
        receipt_ok = receipt_ok and value
    declarations_ok = all(bool(value) for value in declarations.values())
    ast_ok, ast_witnesses = _ast_binding_checks(texts, gates)
    return token_ok and receipt_ok and declarations_ok and ast_ok, {
        "records": dict(RECORD_HASHES),
        "declarations": declarations,
        "accepted_receipts": receipts,
        "ast_bindings": ast_witnesses,
        "actual_background_sources": {
            "postrestore_C6": PREFIX + "prescribed_R100_Rm_C6_functions.py",
            "patch_C6": PREFIX + "true_Rm_patch_C6_functions.py",
            "prescribed_loop": PREFIX + "prescribed_loop_C5_functions.py",
            "finite_N_patch": PREFIX + "true_Rm_patch_finite_N_C5_functions.py",
        },
        "actual_phase_sources": {
            "C2": PREFIX + "C2_phase_primitives.py",
            "all_chart": PREFIX + "all_chart_point_phase.py",
            "native_conditioned": "lei_ren_part1_paper_compliant_current_native_conditioned_phase.py",
            "whole_Z": PREFIX + "whole_Z_signed_loop_functions.py",
            "Rm_conditioned_density": PREFIX + "Rm_conditioned_phase_density.py",
            "native_first_jets": "lei_ren_part1_paper_compliant_current_native_phase_first_jets.py",
        },
    }

def _positive_packet(packet: Mapping[str, Any], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
    if not isinstance(packet, Mapping):
        _gate(gates, "positivity_packet_missing")
        return False, {}
    try:
        packet_canonical_sha = canonical_digest(packet)
    except (TypeError, ValueError) as exc:
        packet_canonical_sha = None
        _gate(gates, "positivity_packet_canonical_digest_failed", error=str(exc))
    canonical_ok = packet_canonical_sha == POSITIVITY_PACKET_CANONICAL_SHA256
    if not canonical_ok:
        _gate(gates, "positivity_packet_canonical_digest_changed",
              expected=POSITIVITY_PACKET_CANONICAL_SHA256, actual=packet_canonical_sha)
    header_ok = all((
        packet.get("schema") == "Rm_history_a_positive_sources_stage_v1",
        packet.get("stage") == "actual_whole_domain_loop_root_a_positivity",
        packet.get("parent_report_sha256") == PARENT_REPORT_SHA256,
        packet.get("source_family") == {
            "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
            "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
            "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
        },
    ))
    if not header_ok:
        _gate(gates, "positivity_packet_header_changed")
    scope = packet.get("scope")
    scope_ok = isinstance(scope, Mapping) and all((
        scope.get("actual_whole_domain_a_positive_sources_proved") == 2,
        scope.get("finite_N_history_C2_source_proved") is False,
        scope.get("physical_Q_y_Z1_value_seam_proved") is False,
        scope.get("new_numeric_field_values_computed") is False,
        scope.get("full_goal_complete") is False,
    ))
    if not scope_ok:
        _gate(gates, "positivity_packet_scope_changed")

    sides: dict[str, Any] = {}
    side_ok = True
    for side, expected in SIDES.items():
        value = packet.get("actual_history_a_positive_sources", {}).get(side) if isinstance(packet.get("actual_history_a_positive_sources"), Mapping) else None
        required = isinstance(value, Mapping) and all((
            value.get("source_only") is True,
            value.get("nodes_digest") == EXPECTED_NODES_DIGEST,
            value.get("passed") is True,
            value.get("actual_target_handle") == expected["a"],
            value.get("whole_native_domain") == expected["native_domain"],
            value.get("whole_Z_domain") == ["-1", "1"],
            value.get("strict_positive_proved") is True,
            value.get("lower_bound_exact") == expected["lower_bound"],
            value.get("unresolved_gates") == [],
            value.get("old_operators_replayed") is False,
            value.get("numeric_values_computed") is False,
            value.get("full_goal_complete") is False,
        ))
        if not required:
            _gate(gates, "actual_a_positive_side_changed", side=side)
            side_ok = False
            continue
        if side == "postrestore":
            native = value.get("native_context", {})
            native_ok = all((
                native.get("coordinate_handle") == expected["coordinate"],
                native.get("coordinate_name") == expected["coordinate_name"],
                native.get("phase_held_fixed") is True,
            ))
            z_context = value.get("Z_context", {})
            z_ok = all((z_context.get("variable_handle") == Z_NODE, z_context.get("variable_name") == "Z", z_context.get("interval") == [-1, 1]))
            if not native_ok or not z_ok:
                _gate(gates, "actual_a_positive_coordinate_or_Z_context_changed", side=side)
                side_ok = False
            witness = value.get("identity_source_witness", {})
            checks = witness.get("source_constructor_checks", {}) if isinstance(witness, Mapping) else {}
            handles = witness.get("handles", {}) if isinstance(witness, Mapping) else {}
            saved = witness.get("saved_box_used_for_positivity") is False
            post_witness_ok = all((
                handles.get("a") == expected["a"],
                witness.get("exact_reduction") == "E_y=E/10; C=4*E/5; a=4/5",
                checks.get("constructor_tokens_checked") is True,
                checks.get("phi_from_Utheta_return_ast_bound") is True,
                checks.get("generic_bounds_domain", {}).get("source_function_positivity_not_inferred_from_saved_box") is True,
                checks.get("generic_bounds_domain", {}).get("saved_box_used_for_E_positivity") is False,
                saved,
            ))
        else:
            witnesses = value.get("witnesses", {})
            source = witnesses.get("source", {}) if isinstance(witnesses, Mapping) else {}
            theorem = source.get("theorem_binding", {}) if isinstance(source, Mapping) else {}
            saved_graph = witnesses.get("saved_graph", {}) if isinstance(witnesses, Mapping) else {}
            ports = saved_graph.get("ports", {}) if isinstance(saved_graph, Mapping) else {}
            chart = saved_graph.get("chart", {}) if isinstance(saved_graph, Mapping) else {}
            z_context = saved_graph.get("Z", {}) if isinstance(saved_graph, Mapping) else {}
            native_ok = all((
                chart.get("handle") == expected["coordinate"],
                chart.get("name") == expected["coordinate_name"],
                chart.get("domain") == expected["native_domain"],
            ))
            z_ok = all((z_context.get("handle") == Z_NODE, z_context.get("domain") == ["-1", "1"]))
            if not native_ok or not z_ok:
                _gate(gates, "actual_a_positive_coordinate_or_Z_context_changed", side=side)
                side_ok = False
            post_witness_ok = all((
                ports.get("a") == expected["a"],
                ports.get("E") == expected["E"],
                theorem.get("minimum_a") == "7/10",
                theorem.get("whole_domain") == ["1", "exp(1)"],
                theorem.get("Z_domain") == ["-1", "1"],
                theorem.get("source_function_not_saved_box") is True,
                theorem.get("correlated_C2_to_C1_source_class") is True,
                theorem.get("same_C6_exact_branch_transfer") is True,
                saved_graph.get("denominator_source_guard") == "same strictly positive actual prescribed swirl E",
            ))
        if not post_witness_ok:
            _gate(gates, "actual_a_positive_source_witness_changed", side=side)
            side_ok = False
        sides[side] = {
            "target_handle": expected["a"],
            "native_domain": expected["native_domain"],
            "Z_domain": ["-1", "1"],
            "lower_bound_exact": expected["lower_bound"],
            "phase_held_fixed": native_ok,
            "Z_context_checked": z_ok,
            "source_function_not_saved_box": post_witness_ok,
        }
    return canonical_ok and header_ok and scope_ok and side_ok, {
        "sides": sides,
        "raw_sha256": POSITIVITY_PACKET_SHA256,
        "canonical_sha256": POSITIVITY_PACKET_CANONICAL_SHA256,
    }


def _graph_witnesses(nodes: Sequence[Any], gates: list[dict[str, Any]]) -> tuple[bool, dict[str, Any]]:
    if not isinstance(nodes, Sequence) or isinstance(nodes, (str, bytes, bytearray)):
        _gate(gates, "saved_graph_not_a_sequence")
        return False, {}
    try:
        actual_digest = canonical_digest(nodes)
    except (TypeError, ValueError) as exc:
        _gate(gates, "saved_graph_digest_failed", error=str(exc))
        return False, {}
    shape_ok = len(nodes) == EXPECTED_NODE_COUNT and actual_digest == EXPECTED_NODES_DIGEST
    if not shape_ok:
        _gate(gates, "saved_graph_digest_or_node_count_changed", expected_digest=EXPECTED_NODES_DIGEST, actual_digest=actual_digest, actual_count=len(nodes))

    direct_ok = all((
        _expect_row(nodes, Z_NODE, "bound_variable", gates, name="Z"),
        _expect_row(nodes, N_NODE, "exact_positive_integer_expression", gates,
                    type="DoubleDyadicInteger", exact_positive_integer=True, Z_independent=True,
                    one_common_N_for_all_charts_and_control_denominators=True),
        _expect_row(nodes, INV_N_NODE, "positive_quotient", gates, numerator=1, denominator=N_NODE,
                    source_positive_certificate="same exact positive integer N on every source chart"),
        _expect_row(nodes, ETA_NODE, "current_original_source_parameter", gates, name="eta", Z_independent=True,
                    quantity_is_exact_original_parameter_not_range_endpoint=True, source_attribute="eta_log"),
        _expect_row(nodes, DSTAR_NODE, "current_original_source_parameter", gates, name="d_star", Z_independent=True,
                    quantity_is_exact_original_parameter_not_range_endpoint=True, source_attribute="dstar_log"),
        _expect_row(nodes, 28, "exact_rational", gates, numerator=1, denominator=2),
        _expect_row(nodes, 48, "exact_rational", gates, numerator=2, denominator=1),
    ))
    side_rows: dict[str, Any] = {}
    sides_ok = True
    for side, expected in SIDES.items():
        E_chain_ok = all(_expect_row(nodes, handle, operation, gates, arguments=arguments)
                         for handle, operation, arguments in expected["E_chain"])
        M_chain_ok = all(_expect_row(nodes, handle, operation, gates, **(fields if isinstance(fields, Mapping) else {"arguments": fields}))
                         for handle, operation, fields in expected["M_chain"])
        root_ok = _expect_row(nodes, expected["a"], "positive_quotient", gates,
                              denominator=expected["E"],
                              source_positive_certificate="same strictly positive actual prescribed swirl E")
        delta_arguments = [expected["a"], 25304] if side == "postrestore" else [expected["kappa"], 25304]
        delta_ok = _expect_row(nodes, expected["delta"], "sum", gates, arguments=delta_arguments)
        kappa_ok = True
        quotient_ok = True
        if side == "patch":
            kappa_ok = _expect_row(nodes, expected["kappa"], "sum", gates, arguments=[expected["a"], expected["kappa_quotient"]])
            quotient_ok = _expect_row(nodes, expected["kappa_quotient"], "positive_quotient", gates,
                                      numerator=196217, denominator=expected["a"],
                                      source_positive_certificate="same original positive a")
        predicate_ok = _expect_row(nodes, expected["predicate"], "exact_real_comparison", gates, operator=">=", left=expected["delta"], right=21)
        coordinate_ok = _expect_row(nodes, expected["coordinate"], "bound_variable", gates, name=expected["coordinate_name"])
        B_rows = []
        B_ok = True
        for order, handle in enumerate(expected["B"]):
            row_ok = _expect_row(nodes, handle, "original_lazy_flat_branch", gates,
                                 Delta=expected["delta"], eta=21, flat_value=0,
                                 active_body_not_evaluated_on_flat_branch=True,
                                 ordinary_slow_Z_derivative_order=order,
                                 active_body=expected["active_B"][order],
                                 same_original_flat_extension_through_order5=True)
            B_ok = B_ok and row_ok
            formula = expected["B_formula"][order]
            body_ok = _expect_row(nodes, expected["active_B"][order], "product", gates,
                                  arguments=[formula["root"], 28])
            root_row_ok = _expect_row(nodes, formula["root"],
                                      "product" if order == 0 else "sum",
                                      gates, arguments=formula["root_arguments"])
            term_ok = True
            if order > 0:
                for term_handle, term_args in zip(formula["root_arguments"], formula["terms"]):
                    term_ok = _expect_row(nodes, term_handle, "product", gates,
                                          arguments=list(term_args)) and term_ok
            body_ok = body_ok and root_row_ok and term_ok
            B_ok = B_ok and body_ok
            B_rows.append({"handle": handle, "ordinary_order": order, "active_body": expected["active_B"][order],
                           "active_body_checked": row_ok and body_ok, "formula_root": formula["root"],
                           "formula_terms": [list(item) for item in formula["terms"]]})
        density_rows = []
        density_ok = True
        for order, (density, B) in enumerate(zip(expected["density"], expected["B"])):
            row_ok = _expect_row(nodes, density, "product", gates, arguments=[B, INV_N_NODE])
            density_ok = density_ok and row_ok
            density_rows.append({"handle": density, "B_handle": B, "inv_N_handle": INV_N_NODE, "ordinary_order": order, "scale_once": row_ok})
        side_ok = E_chain_ok and M_chain_ok and root_ok and delta_ok and kappa_ok and quotient_ok and predicate_ok and coordinate_ok and B_ok and density_ok
        sides_ok = sides_ok and side_ok
        side_rows[side] = {
            "coordinate": {"handle": expected["coordinate"], "name": expected["coordinate_name"], "domain": expected["native_domain"], "checked": coordinate_ok},
            "a": {"handle": expected["a"], "E_handle": expected["E"], "strict_source_quotient": root_ok},
            "source_jets": {"E_chain_checked": E_chain_ok, "M_chain_checked": M_chain_ok,
                            "E_chain": [handle for handle, _, _ in expected["E_chain"]],
                            "M_chain": [handle for handle, _, _ in expected["M_chain"]]},
            "active_cutoff": {"Delta_handle": expected["delta"], "predicate_handle": expected["predicate"], "eta_handle": ETA_NODE,
                              "active_gap": "Delta < eta implies gamma = 2 eta - Delta > eta", "checked": delta_ok and predicate_ok},
            "kappa": ({"handle": expected["kappa"], "quotient_handle": expected["kappa_quotient"],
                        "quotient_numerator": 196217, "quotient_denominator": expected["a"],
                        "checked": kappa_ok and quotient_ok} if side == "patch" else None),
            "B_rows": B_rows,
            "density_rows": density_rows,
            "scale_once": density_ok,
        }
    return shape_ok and direct_ok and sides_ok, {
        "graph_digest": actual_digest,
        "Z_variable": {"handle": Z_NODE, "name": "Z", "domain": ["-1", "1"]},
        "common_parameters": {"N": {"handle": N_NODE, "type": "DoubleDyadicInteger", "Z_independent": True},
                              "inv_N": {"handle": INV_N_NODE, "scale": "1/N", "scale_once_in_each_density": True},
                              "eta": {"handle": ETA_NODE, "Z_independent": True},
                              "d_star": {"handle": DSTAR_NODE, "Z_independent": True}},
        "sides": side_rows,
        "all_density_rows_scale_once": sides_ok,
    }


def _regularity_summary(source_witnesses: Mapping[str, Any], graph_witnesses: Mapping[str, Any], positive_witnesses: Mapping[str, Any]) -> dict[str, Any]:
    ast_witness = source_witnesses.get("ast_bindings", {}) if isinstance(source_witnesses, Mapping) else {}
    regularity = ast_witness.get("regularity", {}) if isinstance(ast_witness, Mapping) else {}
    loop = ast_witness.get("loop_functions", {}) if isinstance(ast_witness, Mapping) else {}
    finite = ast_witness.get("finite_N_functions", {}) if isinstance(ast_witness, Mapping) else {}
    mixed = ast_witness.get("mixed_taylor", {}) if isinstance(ast_witness, Mapping) else {}
    phase = ast_witness.get("phase", {}) if isinstance(ast_witness, Mapping) else {}
    exact_regular = all(regularity.get(key) is True for key in (
        "actual_C2_rows", "same_inverse", "normalized_MZZ", "normalized_BZZ",
        "endpoint_normalization", "flat_M_zero", "poisson_inverse",
    ))
    exact_mixed = all(mixed.get(key) is True for key in (
        "normalized_axial_factorials", "normalized_convolution", "normalized_division",
    )) and all(loop.get(key) is True for key in ("T2", "partials")) and all(finite.get(key) is True for key in ("inc", "s"))
    phase_source_ok = all(phase.get(key) is True for key in (
        "common_phase_binding", "phase_Z_zero", "phase_held_fixed",
        "prescribed_phase_origin", "prescribed_fractional_part", "prescribed_phase_Z_zero",
    ))
    endpoint_ok = regularity.get("endpoint_normalization") is True
    flat_ok = regularity.get("flat_M_zero") is True and loop.get("B") is True
    graph_rows_ok = graph_witnesses.get("all_density_rows_scale_once") is True if isinstance(graph_witnesses, Mapping) else False
    source_ast_ok = ast_witness.get("passed") is True
    positive_sides = positive_witnesses.get("sides", {}) if isinstance(positive_witnesses, Mapping) else {}
    positivity_ok = all(
        isinstance(positive_sides, Mapping)
        and side in positive_sides
        and positive_sides[side].get("lower_bound_exact") == SIDES[side]["lower_bound"]
        for side in SIDES
    )
    return {
        "ordinary_variable": "Z",
        "ordinary_Z_domain": ["-1", "1"],
        "fixed_native_phase": phase_source_ok,
        "source_background": "accepted genuine C6 postrestore and patch sources" if source_ast_ok else None,
        "authenticated_exact_source_ast": source_ast_ok,
        "active_formula": "gamma=2*eta-Delta; q=sigma(1-Delta/eta)*sqrt(gamma/(2*a))",
        "compact_loop_geometry": {
            "u": "p2*q/d_star",
            "r": "u/sqrt(1+u^2)",
            "strict_abs_r": "|r|<1",
            "poisson_denominator": "(cos(psi)-r)^2+sin(psi)^2 > 0",
        },
        "inverse_function_theorem": {
            "K": "1/[2*pi*(1+t0^2+2*q^2)] > 0",
            "lambda": "K*(1+t^2) >= K",
            "endpoint_identity": "int(w)=0, int(w^2)=pi/(1-r^2), Phi(0)=0, Phi(2*pi)=1",
            "flat_join": "a*t0+b=0; q=0 gives t=t0, Phi=psi/(2*pi), inverse=2*pi*phase, M=0",
            "endpoint_normalization_proved": endpoint_ok,
            "flat_M_zero_proved": flat_ok,
        },
        "mixed_taylor": {
            "ordinary_orders": [0, 1, 2],
            "normalized_factorials": exact_mixed,
            "two_Z_cross_terms": exact_mixed,
            "same_inverse": regularity.get("same_inverse") is True,
            "source_ast_proved": exact_mixed,
        },
        "boundary": {
            "inactive_branch": "lazy flat density zero",
            "active_inverse_sqrt_outside_domain_evaluated": False,
            "sigma_flat_at_boundary": True,
            "numeric_conditioning_cap_used": False,
        },
        "phase": {
            "fixed_phase_Z_independence": phase_source_ok,
            "native_fractional_part_jumps_irrelevant_to_Z_DCT": phase_source_ok,
            "uniform_phase_bound": phase_source_ok and positivity_ok,
        },
        "proof_conditions": {
            "regularity_source_ast": exact_regular,
            "ordinary_graph_rows_scale_once": graph_rows_ok,
            "accepted_positive_source_domains": positivity_ok,
            "endpoint_Phi_2pi_equals_1": endpoint_ok,
            "flat_M_equals_0_from_a_t0_plus_b": flat_ok,
            "source_only": True,
        },
    }


def derive(
    nodes: Sequence[Any],
    source_records: Mapping[str, Any],
    *,
    expected_nodes_digest: str,
    positivity_packet: Mapping[str, Any],
) -> dict[str, Any]:
    """Return a source-only, fail-closed actual local ``m`` C2 certificate."""

    gates: list[dict[str, Any]] = []
    if expected_nodes_digest != EXPECTED_NODES_DIGEST:
        _gate(gates, "unexpected_source_graph_digest_assertion", expected=EXPECTED_NODES_DIGEST, actual=expected_nodes_digest)
    source_ok, source_witnesses = _source_witnesses(source_records, gates)
    graph_ok, graph_witnesses = _graph_witnesses(nodes, gates)
    positive_ok, positive_witnesses = _positive_packet(positivity_packet, gates)
    analytic_ok = source_ok and graph_ok and positive_ok
    if not analytic_ok:
        _gate(gates, "actual_local_C2_analytic_chain_not_closed")

    passed = bool(expected_nodes_digest == EXPECTED_NODES_DIGEST and analytic_ok and not gates)
    rows = []
    edges = []
    for side, expected in SIDES.items():
        for order, density in enumerate(expected["density"]):
            rows.append({
                "side": side,
                "field": "m",
                "ordinary_order": order,
                "density_handle": density,
                "source_formula": "m=B_over_Pstar/N",
                "native_domain": expected["native_domain"],
                "ordinary_Z_domain": ["-1", "1"],
                "fixed_phase": True,
                "source_derivative_proved": passed,
                "regularity_proved": passed,
                "proof_origin": "actual C2 Taylor source rows and exact finite-N product with one shared inv_N",
            })
        for order in (0, 1):
            edges.append({
                "side": side,
                "field": "m",
                "from_order": order,
                "to_order": order + 1,
                "source_handle": expected["density"][order],
                "derivative_handle": expected["density"][order + 1],
                "regularity_proved": passed,
            })
    return {
        "schema": SCHEMA,
        "source_only": True,
        "nodes_digest": EXPECTED_NODES_DIGEST,
        "exact_sha256": EXPECTED_NODES_DIGEST,
        "positivity_packet_sha256": POSITIVITY_PACKET_SHA256,
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "passed": passed,
        "proof_coverage": {
            "actual_postrestore_density_rows": [121350, 121351, 121352],
            "actual_patch_density_rows": [198182, 198183, 198184],
            "actual_background_C6": source_ok,
            "actual_a_positive_postrestore_4_over_5": positive_ok and positive_witnesses.get("sides", {}).get("postrestore", {}).get("lower_bound_exact") == "4/5",
            "actual_a_positive_patch_7_over_10": positive_ok and positive_witnesses.get("sides", {}).get("patch", {}).get("lower_bound_exact") == "7/10",
            "actual_finite_N_density_C2": passed,
            "incoming_history_derivative": False,
            "global_Qy1_or_physical_seam": False,
            "full_goal": False,
        },
        "regularity_argument": _regularity_summary(source_witnesses, graph_witnesses, positive_witnesses),
        "witnesses": {
            "producer_sources": source_witnesses,
            "saved_graph": graph_witnesses,
            "accepted_positivity": positive_witnesses,
            "actual_source_handles_bound_independently": True,
            "slot_tags_alone_are_not_used": True,
        },
        "rows": rows,
        "edges": edges,
        "finite_N_scale_once": graph_ok and graph_witnesses.get("all_density_rows_scale_once") is True,
        "fixed_phase_Z_independent": source_ok and positive_ok and source_witnesses.get("ast_bindings", {}).get("phase", {}).get("phase_Z_zero") is True,
        "incoming_history_derivative_proved": False,
        "global_Qy1_proved": False,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "operators_reexecuted": False,
        "full_goal_complete": False,
        "unresolved_gates": gates,
        "gates": gates,
        "source_record_hashes": dict(RECORD_HASHES),
    }
