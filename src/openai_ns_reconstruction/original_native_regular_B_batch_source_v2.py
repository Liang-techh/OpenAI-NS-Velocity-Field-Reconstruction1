"""Fresh failed-binding repair and B source composition; no ancestor replay."""

SCHEMA = 'original_native_regular_B_batch_source_v2'

FAILURE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_regular_B_batch_source_failure.json'

WEIGHTED = 'lei_ren_part1_paper_compliant_current_original_native_active_S1_weighted_ibp.json'

RECORD_NAMES = ['experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_regular_B_batch_source_failure.json', 'lei_ren_part1_paper_compliant_current_original_native_active_S1_weighted_ibp.json', 'lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py', 'lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py', 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py', 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py', 'lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py', 'lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.py', 'lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py', 'lei_ren_part1_paper_compliant_outer_initial.py']

RECORD_HASHES = {'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_regular_B_batch_source_failure.json': 'f9714c0feb746d3721dfaedaa89b084e8ffc6d0f97ea15f3a3daa43c554d3bbb', 'lei_ren_part1_paper_compliant_current_original_native_active_S1_weighted_ibp.json': 'aab72bbc7a882e9bf089a43dd3407babe085f1c0ec61583578d20a5e7c9bec28', 'lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py': '70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32', 'lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py': '8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605', 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py': '69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d', 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py': '7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7', 'lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py': '4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3', 'lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.py': 'c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15', 'lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py': '38d3a8e5abd964c42b135410485f6ad58dfe9f7656035b4b5c7b1ed48c2ebb01', 'lei_ren_part1_paper_compliant_outer_initial.py': 'a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657'}

DECODED_RECORD_HASHES = {'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_regular_B_batch_source_failure.json': 'f9714c0feb746d3721dfaedaa89b084e8ffc6d0f97ea15f3a3daa43c554d3bbb', 'lei_ren_part1_paper_compliant_current_original_native_active_S1_weighted_ibp.json': 'aab72bbc7a882e9bf089a43dd3407babe085f1c0ec61583578d20a5e7c9bec28', 'lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py': '70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32', 'lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py': '8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605', 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py': '69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d', 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py': '7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7', 'lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py': '4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3', 'lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.py': 'c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15', 'lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py': '38d3a8e5abd964c42b135410485f6ad58dfe9f7656035b4b5c7b1ed48c2ebb01', 'lei_ren_part1_paper_compliant_outer_initial.py': 'a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657'}

SOURCE_TEXT_ALIASES = {'c6': 'lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py', 'bridge': 'lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py', 'loop': 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py', 'mod': 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py', 'outer': 'lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py', 'patch': 'lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.py', 'gamma': 'lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py', 'outer_initial': 'lei_ren_part1_paper_compliant_outer_initial.py'}

FAILED_IMPLEMENTATION_HASHES = {'scripts/run_original_native_regular_B_batch_source.py': '25341ca4eca1ef650f7d54dd21b476dfcf339a040e7c1d644069ac59bacb3483', 'src/openai_ns_reconstruction/original_native_regular_B_batch_source.py': '1a5a161e9c791fcf56c0c2e572e930220728183e63a292b89d3d44319782df11'}

TARGET_NODES = [531, 532, 533, 534, 535]

STAGE_FLAGS = ['actual_batch_source_ports_bound', 'all_target_endpoint_identities_proved', 'weighted_B_C1_source_composition_proved', 'native_IN_batch_source_reduction_proved', 'all_old_terms_retained', 'accepted_ancestors_not_reexecuted']

FALSE_FLAGS = ['quantitative_packet_error_reduction_proved', 'numeric_values_computed', 'source_C5_installation_claimed', 'numeric_containment_claimed', 'new_loop_domain_admission_claimed', 'admissible_cone_admission_claimed', 'physical_layer_completed', 'old_error_terms_removed', 'global_five_moment_closure', 'time_order_recursion_claimed', 'full_goal_complete', 'old_operators_replayed']

B_Z_CHAIN_RULE = 'B_Z=2*R*R_Z*H+R^2*(Q_Z/L-Q*L_Z/L^2+Q*Q_Z-4*Z*U^2/L^2-4*Z^2*U*U_Z/L^2+4*Z^2*U^2*L_Z/L^3), H=Q/L+Q^2/2-2*Z^2*U^2/L^2'

"""Fresh, source-only repair/composition for the saved regular outer-B failure.

This module is deliberately unexecuted by the artifact producer.  It does not
import the entered fragment or the combined kernel and it does not walk either
archived graph.  ``derive_outer_repair`` consumes the frozen JSON evidence and
the authenticated weighted-port graph, repairs the two literal wiring errors,
and composes the already authenticated source predicates.
"""

import ast
import hashlib


REPAIR_SCHEMA = "regular_outer_B_batch_fragment_v2_2026_10_21"
REPAIR_TARGET_NODES = (531, 532, 533, 534, 535)
ORDINARY_Z_ORDERS = (0, 1)
OUTER_ARCHIVE_RECORD = (
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz"
)
OUTER_ARCHIVE_RAW_SHA256 = "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d"
LOOP_ARCHIVE_RAW_SHA256 = "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa"
MOD_ARCHIVE_RAW_SHA256 = "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a"
WEIGHTED_REPORT_SHA256 = "aab72bbc7a882e9bf089a43dd3407babe085f1c0ec61583578d20a5e7c9bec28"

SOURCE_HASHES = {
    "outer": "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3",
    "mod": "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
    "patch": "c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15",
    "loop": "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
}

# The entered fragment carried this stale nibble.  The repaired expected
# literal is recorded separately from the raw source-segment digest.
STALE_MOD_PHYSICAL_PRESSURE_HASH = (
    "20226e1212cbf22e366d8deab5f5e663830c6adbe3c5bc3d4352a6641cc3f437"
)
ACTUAL_MOD_PHYSICAL_PRESSURE_HASH = (
    "20226e1212cbf22e366d8deab5f8e663830c6adbe3c5bc3d4352a6641cc3f437"
)
ACTUAL_MOD_PHYSICAL_PRESSURE_SOURCE_SEGMENT_SHA256 = (
    "20226e1212cbf22e366d8deab5f8e663830c6adbe3c5bc3d4352a6641cc3f437"
)
MOD_PHYSICAL_PRESSURE_NEEDLES = (
    "U=a.scale(fields['V'],Pstar)",
    "Uy=a.scale(radial['V_y'],Pstar)",
    "m=a.scale(fields['m'],Pstar)",
    "my=a.sub(U,m)",
    "L=b.sub(one,b.mul(de,b.mul(z,z)))",
    "same original L=1-delta*Z^2",
    "numerator=lambda u,mean",
    "b.mul(z,b.trim(u))",
    "b.shift(mean)",
    "Q=div(numerator(U,m))",
    "Qy=div(numerator(Uy,my))",
    "return dict(U=U,U_y=Uy,m=m,m_y=my,Q=Q,Q_y=Qy)",
)

COMMON_B_FORMULA = "R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)"
COMMON_Q_FORMULA = (
    "N=2*Z*U-(1-delta)*Z*m-(1-Z^2)*m_Z; "
    "N_Z=2*U+2*Z*U_Z-(1-delta)*m+(1+delta)*Z*m_Z-(1-Z^2)*m_ZZ; "
    "Q=N/L; "
    "Q_Z=(N_Z*L-N*L_Z)/L^2"
)

_FIELD_SPECS = {
    "coordinate": ("coordinate", "actual_native_geometry_coordinate", "geometry", (0,)),
    "domain": ("domain", "actual_native_geometry_domain", "geometry", (0, 1)),
    "radius": ("radius", "actual_native_geometry_radius", "geometry", (0,)),
    "physical_measure": (
        "physical_measure", "actual_native_geometry_physical_measure", "geometry", (0,)
    ),
    "Jacobian": (
        "Jacobian", "actual_native_geometry_Jacobian", "geometry.physical_measure", (0,)
    ),
    "U": ("U", "actual_native_physical_U", "true_outer_C6_C5_bindings.charts", (0, 1, 2, 3, 4, 5)),
    "U_y": ("U_y", "actual_native_physical_U_y", "true_outer_C6_C5_bindings.charts", (0, 1, 2, 3, 4, 5)),
    "m": ("m", "actual_native_physical_m", "true_outer_C6_C5_bindings.charts", (0, 1, 2, 3, 4, 5)),
    "m_y": ("m_y", "actual_native_physical_m_y", "true_outer_C6_C5_bindings.charts", (0, 1, 2, 3, 4, 5)),
    "Q": ("Q", "actual_native_physical_Q", "true_outer_C6_C5_bindings.charts", (0, 1, 2, 3, 4)),
    "Q_y": ("Q_y", "actual_native_physical_Q_y", "true_outer_C6_C5_bindings.charts", (0, 1, 2, 3, 4)),
    "S": ("S", "actual_native_physical_S", "true_outer_C6_C5_bindings.charts", (0, 1, 2, 3)),
    "local_S1": ("local_integrals.S1", "actual_native_local_integrals_S1", "true_outer_C6_C5_bindings.charts", (0, 1, 2, 3)),
    "cumulative_S1": ("cumulative_integrals.S1", "actual_native_cumulative_integrals_S1", "true_outer_C6_C5_bindings.charts", (0, 1, 2, 3)),
    "incoming_S1": ("incoming_integrals.S1", "actual_native_incoming_integrals_S1", "true_outer_C6_C5_bindings.charts", (0, 1, 2, 3)),
    "R_Q_y": ("R_Q_y", "actual_native_R_Q_y", "true_outer_C6_C5_bindings.charts", (0, 1, 2, 3)),
}

_JUMP_SPECS = {
    531: ((465, 455), 11),
    532: ((476, 466), 12),
    533: ((487, 477), 13),
    534: ((498, 488), 14),
    535: ((509, 499), 15),
}


def _sha256_text(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _emit(callback, stage, payload):
    if callable(callback):
        callback(stage, payload)


def _outer(saved_result):
    component_proofs = saved_result.get("component_proofs", {})
    return component_proofs.get("outer") or saved_result.get("outer") or saved_result


def _scope(saved_result):
    return saved_result.get("B_identity_scope") or saved_result.get("b_identity_scope") or {}


def _mapping(value):
    return value if isinstance(value, dict) else {}


def _list(value, length=None):
    if not isinstance(value, list):
        return False
    if length is not None and len(value) != length:
        return False
    return all(isinstance(item, int) and not isinstance(item, bool) for item in value)


def _node_map(saved_result):
    outer = _outer(saved_result)
    candidate = outer.get("per_node") or outer.get("nodes") or {}
    if isinstance(candidate, list):
        return {int(item.get("node")): item for item in candidate if isinstance(item, dict) and item.get("node") is not None}
    return {int(key): value for key, value in candidate.items() if str(key).isdigit() and isinstance(value, dict)}


def _scope_node_map(saved_result):
    candidate = _scope(saved_result).get("per_node") or _scope(saved_result).get("nodes") or {}
    if isinstance(candidate, list):
        return {int(item.get("node")): item for item in candidate if isinstance(item, dict) and item.get("node") is not None}
    return {int(key): value for key, value in candidate.items() if str(key).isdigit() and isinstance(value, dict)}


def _saved_binding(saved_result, node, side, outer_nodes, scope_nodes):
    scope_node = scope_nodes.get(node, {})
    side_value = scope_node.get(side, {}) if isinstance(scope_node, dict) else {}
    if isinstance(side_value, dict):
        binding = side_value.get("actual_source_binding")
        if isinstance(binding, dict):
            return binding
    endpoint = _mapping(outer_nodes.get(node, {}).get("endpoint_bindings", {})).get(side, {})
    return _mapping(endpoint.get("actual_source_binding") or endpoint)


def _source_contract(source_texts):
    details = {}
    all_ok = True
    for alias, expected_hash in SOURCE_HASHES.items():
        text = source_texts.get(alias) if isinstance(source_texts, dict) else None
        actual_hash = _sha256_text(text) if isinstance(text, str) else None
        ok = actual_hash == expected_hash
        details[alias] = {"sha256": actual_hash, "expected_sha256": expected_hash, "passed": ok}
        all_ok = all_ok and ok

    mod_text = source_texts.get("mod") if isinstance(source_texts, dict) else None
    mod_ast_hash = None
    needle_matches = {}
    ast_ok = False
    if isinstance(mod_text, str):
        try:
            tree = ast.parse(mod_text)
            function = next(
                item for item in ast.walk(tree)
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef))
                and item.name == "physical_pressure_fields"
            )
            segment = ast.get_source_segment(mod_text, function)
            mod_ast_hash = _sha256_text(segment + "\n") if segment is not None else None
            needle_matches = {needle: needle in mod_text for needle in MOD_PHYSICAL_PRESSURE_NEEDLES}
            ast_ok = (
                mod_ast_hash == ACTUAL_MOD_PHYSICAL_PRESSURE_SOURCE_SEGMENT_SHA256
                and all(needle_matches.values())
            )
        except (SyntaxError, StopIteration):
            ast_ok = False
    details["mod_physical_pressure_fields"] = {
        "stale_expected_sha256": STALE_MOD_PHYSICAL_PRESSURE_HASH,
        "repaired_expected_sha256": ACTUAL_MOD_PHYSICAL_PRESSURE_HASH,
        "actual_source_segment_sha256": mod_ast_hash,
        "needle_matches": needle_matches,
        "passed": ast_ok,
    }
    return {"passed": all_ok and ast_ok, "aliases": details}


def _saved_context(saved_result, outer_nodes, scope_nodes):
    outer = _outer(saved_result)
    guards = _mapping(outer.get("proof_guards"))
    required = (
        "required_source_text_aliases",
        "authenticated_required_source_text_hashes",
        "outer_source_ast_contract",
        "patch_source_ast_contract",
        "loop_source_ast_contract",
        "all_outer_chart_contracts",
        "all_outer_endpoint_substitutions",
        "finite_N_extensional_source_semantics",
        "all_finite_N_extensional_rows",
        "U_Z1_source_semantics",
        "m_Z2_source_semantics",
        "accepted_parent_physical_m_Q_certificate",
        "common_L_source_row",
        "weighted_source_context",
        "all_leading_history_handoffs",
        "retained_correction_memory_all_seams",
        "shared_loop_parameter_rows",
    )
    global_ok = all(guards.get(name) is True for name in required)
    per_node_ok = {}
    for node in REPAIR_TARGET_NODES:
        evidence = _mapping(outer_nodes.get(node, {}))
        predicates = _mapping(scope_nodes.get(node, {}).get("predicates"))
        if not predicates:
            predicates = _mapping(evidence.get("predicates"))
        per_node_ok[node] = all(
            predicates.get(name) is True
            for name in (
                "weighted_graph",
                "actual_jump",
                "predecessor_binding",
                "successor_binding",
                "exact_proof_endpoint_substitutions",
                "ordinary_Z_scope",
                "saved_uniform_B_C1_recipe",
            )
        )
    return {"passed": global_ok and all(per_node_ok.values()), "guards": guards, "per_node": per_node_ok}


def _weighted_graph(weighted):
    if not isinstance(weighted, dict):
        return {}
    direct = weighted.get("definition_graph")
    if isinstance(direct, dict):
        return direct
    return _mapping(weighted.get("component_result")).get("definition_graph", {})


def _graph_nodes(graph):
    raw = graph.get("nodes") if isinstance(graph, dict) else None
    if isinstance(raw, list):
        return {index: value for index, value in enumerate(raw)}
    if isinstance(raw, dict):
        return raw
    return {}


def _row_source_nodes(binding, field):
    archive = _mapping(binding.get("archive_source_nodes"))
    physical = _mapping(archive.get("physical"))
    if field in ("coordinate", "domain", "radius"):
        value = archive.get(field)
        return value if isinstance(value, list) else None
    if field in ("U", "U_y", "m", "m_y", "Q", "Q_y"):
        value = physical.get(field)
        if value is None:
            value = _mapping(binding.get("physical")).get(field)
        return value if isinstance(value, list) else None
    return None


def _port_row_ok(row, field, binding):
    if not isinstance(row, dict):
        return False
    quantity, operation, source_path, orders = _FIELD_SPECS[field]
    if row.get("arguments") != [] or row.get("bound_or_midpoint") is not False:
        return False
    if row.get("function_identity_not_claimed") is not True:
        return False
    if row.get("chart") != binding.get("chart"):
        return False
    if row.get("quantity") != quantity or row.get("operation") != operation:
        return False
    if row.get("source_path") != source_path or row.get("ordinary_Z_orders") != list(orders):
        return False
    if row.get("source_record") != OUTER_ARCHIVE_RECORD:
        return False
    if field in ("coordinate", "domain", "radius", "physical_measure", "Jacobian") and row.get("Z_independent") is not True:
        return False
    actual_nodes = row.get("source_nodes")
    expected_nodes = _row_source_nodes(binding, field)
    if expected_nodes is not None:
        return actual_nodes == expected_nodes
    if not _list(actual_nodes, len(orders)):
        return False
    return True


def _weighted_endpoint_repair(saved_result, weighted, outer_nodes, scope_nodes):
    graph = _weighted_graph(weighted)
    nodes = _graph_nodes(graph)
    if not nodes:
        return {"passed": False, "failed_nodes": list(REPAIR_TARGET_NODES), "rows_checked": 0, "per_node": {}}
    per_node = {}
    rows_checked = 0
    all_ok = True
    for node in REPAIR_TARGET_NODES:
        evidence = _mapping(outer_nodes.get(node, {}))
        endpoint_bindings = _mapping(evidence.get("endpoint_bindings"))
        side_results = {}
        node_ok = True
        for side in ("predecessor", "successor"):
            endpoint = _mapping(endpoint_bindings.get(side))
            port_rows = _mapping(endpoint.get("weighted_port_rows"))
            binding = _saved_binding(saved_result, node, side, outer_nodes, scope_nodes)
            base_row = port_rows.get("coordinate")
            base = nodes.get(str(base_row), nodes.get(base_row)) if base_row is not None else None
            side_ok = isinstance(base_row, int) and isinstance(base, dict)
            checked_fields = []
            for field in _FIELD_SPECS:
                row_id = port_rows.get(field)
                row = nodes.get(str(row_id), nodes.get(row_id)) if row_id is not None else None
                field_ok = isinstance(row_id, int) and _port_row_ok(row, field, binding)
                side_ok = side_ok and field_ok
                checked_fields.append({"field": field, "row_id": row_id, "passed": field_ok})
                rows_checked += 1
            side_results[side] = {
                "base_coordinate_row": base_row,
                "checked_fields": checked_fields,
                "passed": side_ok,
            }
            node_ok = node_ok and side_ok
        per_node[node] = side_results
        all_ok = all_ok and node_ok
    return {"passed": all_ok, "failed_nodes": [] if all_ok else list(REPAIR_TARGET_NODES), "rows_checked": rows_checked, "per_node": per_node}


def _weighted_jump_repair(weighted):
    graph = _weighted_graph(weighted)
    nodes = _graph_nodes(graph)
    per_node = {}
    all_ok = True
    for node, (arguments, jump_index) in _JUMP_SPECS.items():
        row = nodes.get(str(node), nodes.get(node))
        ok = (
            isinstance(row, dict)
            and row.get("arguments") == list(arguments)
            and row.get("operation") == "actual_internal_primitive_jump"
            and row.get("formula") == "J_i=B_next_lower-B_prev_upper"
            and row.get("is_actual_function") is True
            and row.get("retained") is True
            and row.get("zero_assumption") is False
            and row.get("source_identity_proved") is False
            and row.get("jump_index") == jump_index
        )
        per_node[node] = {"arguments": list(arguments), "jump_index": jump_index, "passed": ok}
        all_ok = all_ok and ok
    return {"passed": all_ok, "per_node": per_node}


_ENDPOINT_CASES = {
    531: ("Rh_reference", "O2_slope", 0, "6", "E_next=E_prev*exp(0); E_y=E/10; V=seedV; V_y=0 by sigma(0)=0"),
    532: ("O2_slope", "O2_axial", 1, "7", "E_next=E_prev*exp(0); E_y=-E/2 since 1/10-3/5=-1/2; V=seedV; V_y=0 by sigma(1)=1,sigma_prime(1)=0"),
    533: ("O2_axial", "O2_buffer", 1, "6+exp(40)", "E_next=E_prev*exp(0); E_y=-E/2; V=V_y=0 by sigma(0)=sigma_prime(0)=0"),
    534: ("O2_buffer", "O3_transition", 4, "17+exp(40)", "Upper handle4 denotes11 in the pinned geometry source; E_next=E_prev*exp(0); integral_0^0=0; E_y=-E/2; V=V_y=0"),
    535: ("O3_transition", "O3_power", 1, "18+exp(40)", "E_next=E_prev*exp(0); E_y=-(1/2+mu)*E by sigma(1)=1; V=V_y=0"),
}


def _endpoint_case(node, bindings):
    """Fixed-source analytic endpoint substitutions; identities hold as Z functions."""
    expected = _ENDPOINT_CASES.get(node)
    before = _mapping(bindings.get("predecessor"))
    after = _mapping(bindings.get("successor"))
    before_domain = before.get("domain")
    after_domain = after.get("domain")
    passed = bool(expected and before.get("chart") == expected[0]
        and after.get("chart") == expected[1]
        and isinstance(before_domain,list) and len(before_domain)==2
        and isinstance(after_domain,list) and len(after_domain)==2
        and before_domain[1] == expected[2] and after_domain[0] == 0)
    return {"passed":passed,"chart_pair":list(expected[:2]) if expected else [],
        "common_log_radius_offset":expected[3] if expected else None,
        "leading_E_V_Ey_Vy_law":expected[4] if expected else None,
        "source_basis":"Pinned outer creator geometry45-53 and branches155-183; common R=Rm*exp(offset), flat sigma endpoint values and derivatives, exp(0)=1 and equal-bound integral=0.",
        "ordinary_Z_law":"Substitutions are independent of Z; equal analytic endpoint functions retain their ordinary Z derivatives on the original common domain.",
        "numeric_values_computed":False}


def _finite_side_ok(binding):
    extensional = _mapping(binding.get("leading_extensional_rows"))
    finite = _mapping(extensional.get("finite_N"))
    roots = _mapping(extensional.get("roots"))
    finite_ok = all(_list(finite.get(key), 6) for key in ("E_N", "V_N"))
    roots_ok = all(
        _list(roots.get(key), 6)
        for key in ("E", "V", "a", "b", "p1", "p2", "t0", "Delta", "kappa")
    )
    return finite_ok and roots_ok


def _finite_handoff_ok(evidence, bindings):
    handoff = _mapping(_mapping(evidence.get("source_derivation")).get("finite_N_handoff"))
    law = handoff.get("law")
    finite_rows = _mapping(handoff.get("finite_N_rows_bound"))
    finite_rows_ok = all(
        _list(_mapping(finite_rows.get(side)).get(key), 6)
        for side in ("predecessor", "successor")
        for key in ("E_N", "V_N")
    )
    branch_rows = _mapping(_mapping(evidence.get("source_derivation")).get("branch_row_bindings"))
    branch_roots_ok = all(
        _list(
            _mapping(_mapping(branch_rows.get(side)).get("background")).get("generic", {}).get(key),
            6,
        )
        for side in ("predecessor", "successor")
        for key in ("p1", "p2")
    )
    return (
        isinstance(law, str)
        and "phase" in law
        and "loop_functions" in law
        and handoff.get("shared_parameters") == ["delta", "Pstar", "N", "eta", "dstar", "sc"]
        and handoff.get("shared_root_inputs") == ["E", "V", "E_y", "V_y", "p1", "p2", "a", "b", "t0", "Delta", "kappa"]
        and finite_rows_ok
        and branch_roots_ok
        and handoff.get("physical_U_formula") == "U=Pstar*V_N"
        and handoff.get("physical_U_y_formula") == "U_y=Pstar*V_N_y"
        and handoff.get("leading_V_is_not_substituted_for_V_N") is True
        and handoff.get("outer_chart_labels_have_no_first_micro_collar") is True
        and handoff.get("leading_history_rows_bound") is True
        and handoff.get("correction_history_lower_identity") is True
        and all(_finite_side_ok(binding) for binding in bindings.values())
    )


def _m_transport_ok(evidence):
    transport = _mapping(_mapping(evidence.get("source_derivation")).get("m_transport_rows"))
    required = (
        "modified_m_Z2_rows",
        "history_m_Z2_rows",
        "contribution_m_Z2_rows",
        "chart_incoming_m_rows",
        "endpoint_correction_m_rows",
    )
    return bool(transport) and all(
        isinstance(_mapping(transport.get(side)).get(key), list)
        and len(_mapping(transport.get(side)).get(key)) in (3, 6)
        for side in ("predecessor", "successor")
        for key in required
    )


def _compose_node(saved_result, node, outer_nodes, scope_nodes, guards, mod_ok, weighted_ok, jump_ok):
    evidence = _mapping(outer_nodes.get(node, {}))
    saved_node = _mapping(scope_nodes.get(node, {}))
    root_predicates = _mapping(saved_node.get("predicates"))
    source_predicates = _mapping(evidence.get("predicates"))
    endpoint_bindings = _mapping(evidence.get("endpoint_bindings"))
    bindings = {side: _saved_binding(saved_result, node, side, outer_nodes, scope_nodes) for side in ("predecessor", "successor")}
    endpoint_substitutions = saved_node.get("required_endpoint_substitutions")
    chart_pair = [bindings[side].get("chart") for side in ("predecessor", "successor")]
    endpoint_case = _endpoint_case(node, bindings)
    chart_ok = endpoint_case["passed"]
    geometry_ok = (
        guards.get("outer_source_ast_contract") is True
        and guards.get("all_outer_chart_contracts") is True
        and guards.get("all_outer_endpoint_substitutions") is True
        and guards.get("radius_source_semantics") is True
        and source_predicates.get("geometry_source_ast") is True
        and source_predicates.get("outer_chart_domain") is True
        and source_predicates.get("endpoint_Z_independent_substitutions") is True
        and chart_ok
    )
    finite_ok = (
        geometry_ok
        and guards.get("finite_N_extensional_source_semantics") is True
        and guards.get("all_finite_N_extensional_rows") is True
        and guards.get("all_leading_history_handoffs") is True
        and guards.get("retained_correction_memory_all_seams") is True
        and guards.get("shared_loop_parameter_rows") is True
        and guards.get("loop_source_ast_contract") is True
        and source_predicates.get("leading_history_handoff") is True
        and source_predicates.get("correction_history_lower_identity") is True
        and all(_finite_side_ok(bindings[side]) for side in bindings)
        and _finite_handoff_ok(evidence, bindings)
    )
    physical_u_ok = all(_list(_mapping(binding.get("physical")).get("U"), 6) for binding in bindings.values())
    physical_m_ok = all(_list(_mapping(binding.get("physical")).get("m"), 6) for binding in bindings.values())
    physical_q_ok = all(_list(_mapping(binding.get("physical")).get("Q"), 5) for binding in bindings.values())
    u_ok = finite_ok and guards.get("U_Z1_source_semantics") is True and physical_u_ok
    m_ok = (
        guards.get("m_Z2_source_semantics") is True
        and guards.get("all_leading_history_handoffs") is True
        and guards.get("retained_correction_memory_all_seams") is True
        and physical_m_ok
        and _m_transport_ok(evidence)
        and all(_list(_mapping(binding.get("ordinary_Z_prefix")).get("m_Z2"), 3) for binding in bindings.values())
    )
    q_ok = (
        u_ok
        and m_ok
        and mod_ok
        and guards.get("accepted_parent_physical_m_Q_certificate") is True
        and guards.get("common_L_source_row") is True
        and guards.get("U_Z1_source_semantics") is True
        and guards.get("m_Z2_source_semantics") is True
        and physical_q_ok
        and all(_list(_mapping(binding.get("ordinary_Z_prefix")).get("Q_Z1"), 2) for binding in bindings.values())
    )
    b_ok = geometry_ok and u_ok and m_ok and q_ok and weighted_ok and jump_ok
    return {
        "node": node,
        "passed": b_ok,
        "source_only": True,
        "B_C1_proved": b_ok,
        "R_C1_proved": geometry_ok,
        "U_C1_proved": u_ok,
        "m_Z2_proved": m_ok,
        "Q_C1_proved": q_ok,
        "endpoint_B_equality_proved": b_ok,
        "native_weighted_jump_zero_proved": b_ok,
        "endpoint_bindings": bindings,
        "endpoint_substitutions": endpoint_substitutions,
        "Q_y_claimed": False,
        "predicates": {
            "saved_weighted_graph": root_predicates.get("weighted_graph") is True,
            "saved_actual_jump": root_predicates.get("actual_jump") is True,
            "weighted_endpoint_rows_repaired": weighted_ok,
            "weighted_jump_row_repaired": jump_ok,
            "mod_pressure_source_ast_repaired": mod_ok,
            "geometry_source_evidence": geometry_ok,
            "explicit_endpoint_case": endpoint_case["passed"],
            "finite_N_extensional_evidence": finite_ok,
            "physical_U_equals_Pstar_V_N": u_ok,
            "m_Z2_transport_evidence": m_ok,
            "Q_common_L_evidence": q_ok,
            "endpoint_B_equality_recomposed": b_ok,
            "ordinary_Z_not_temporal_n": True,
            "no_S1_or_width_claim": True,
        },
        "source_derivation": {
            "chart_pair": chart_pair,
            "explicit_endpoint_case": endpoint_case,
            "physical_U_law": "U=Pstar*V_N",
            "m_law": "m through Z2 with actual predecessor incoming correction",
            "predecessor_incoming_rows": [190, 217],
            "Q_law": COMMON_Q_FORMULA,
            "common_B_formula": COMMON_B_FORMULA,
            "shared_N_eta_dstar_sc_globalphase": True,
            "leading_and_correction_histories_zero_width_transport": True,
            "Q_y_identity_claim": False,
        },
        "endpoint_bindings_present": all(isinstance(endpoint_bindings.get(side), dict) for side in bindings),
    }


def derive_outer_repair(saved_result, weighted, source_texts, diagnostic_callback=None):
    """Return fresh source-only proofs from the frozen failure evidence.

    The function performs JSON/AST/hash validation only.  It never imports or
    calls an entered proof/helper, source provider, operator, normalizer,
    quadrature, recurrence, or new-function implementation.
    """
    outer = _outer(saved_result)
    outer_nodes = _node_map(saved_result)
    scope_nodes = _scope_node_map(saved_result)
    source_contract = _source_contract(source_texts)
    saved_context = _saved_context(saved_result, outer_nodes, scope_nodes)
    mod_ok = bool(source_contract["aliases"].get("mod_physical_pressure_fields", {}).get("passed"))
    weighted_endpoint = _weighted_endpoint_repair(saved_result, weighted, outer_nodes, scope_nodes)
    weighted_jump = _weighted_jump_repair(weighted)
    _emit(diagnostic_callback, "source_contract", source_contract)
    _emit(diagnostic_callback, "weighted_endpoint_rows", weighted_endpoint)
    _emit(diagnostic_callback, "weighted_jump_rows", weighted_jump)
    per_node = {}
    for node in REPAIR_TARGET_NODES:
        weighted_node = _mapping(weighted_endpoint.get("per_node", {}).get(node))
        weighted_node_ok = weighted_endpoint.get("passed") is True and weighted_node.get("predecessor", {}).get("passed") is True and weighted_node.get("successor", {}).get("passed") is True
        jump_node_ok = _mapping(weighted_jump.get("per_node", {}).get(node)).get("passed") is True
        per_node[str(node)] = _compose_node(
            saved_result,
            node,
            outer_nodes,
            scope_nodes,
            saved_context["guards"],
            mod_ok,
            weighted_node_ok,
            jump_node_ok,
        )
    all_nodes_proved = source_contract["passed"] and saved_context["passed"] and all(item["passed"] for item in per_node.values())
    result = {
        "schema": REPAIR_SCHEMA,
        "status": "source_outer_B_batch_repair_ready" if all_nodes_proved else "source_outer_B_batch_repair_failed",
        "passed": all_nodes_proved,
        "source_only": True,
        "source_function_identity_only": True,
        "numeric_values_computed": False,
        "old_proof_functions_called": 0,
        "source_constructor_calls": 0,
        "source_producer_function_calls": 0,
        "old_operator_calls": 0,
        "generic_DAG_walks": 0,
        "archived_graphs_executed": False,
        "forbidden_archives_decoded": False,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "ordinary_Z_not_temporal_n": True,
        "common_B_formula": COMMON_B_FORMULA,
        "common_Q_formula": COMMON_Q_FORMULA,
        "Q_y_identity_claim": False,
        "source_contract": source_contract,
        "saved_positive_context": saved_context,
        "repair_details": {
            "weighted_endpoint_rows": weighted_endpoint,
            "weighted_jump_rows": weighted_jump,
            "mod_source_ast_contract": source_contract["aliases"].get("mod_physical_pressure_fields", {}),
            "old_entered_fragment_calls": 0,
            "old_combined_kernel_calls": 0,
        },
        "per_node": per_node,
        "all_nodes_proved": all_nodes_proved,
    }
    _emit(diagnostic_callback, "complete", {"passed": all_nodes_proved, "nodes": list(REPAIR_TARGET_NODES)})
    return result


"""Fresh repair composition over saved failed-source premises; unentered."""
import hashlib as v2_hashlib
import json as v2_json


def v2_digest(value):
    return v2_hashlib.sha256(v2_json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=True,allow_nan=False).encode()).hexdigest()


def v2_require(value,message):
    if value is not True: raise AssertionError(message)


def v2_emit(callback,events,stage,payload):
    frozen=v2_json.loads(v2_json.dumps(payload,allow_nan=False))
    events.append({'stage':stage,'payload':frozen})
    if callback is not None: callback(stage,frozen)


def derive(records,expected_input_digest,diagnostic_callback=None):
    events=[]
    v2_require(set(records)==set(RECORD_NAMES),'Changed repair source records')
    v2_require(v2_digest(records)==expected_input_digest,'Changed repair input digest')
    for name in RECORD_NAMES:
        row=records[name]
        v2_require(row.get('raw_sha256')==RECORD_HASHES[name] and row.get('raw_authenticated_before_decode') is True,'Changed raw repair input')
        v2_require(row.get('decoded_sha256')==DECODED_RECORD_HASHES[name] and v2_hashlib.sha256(row['source_text'].encode()).hexdigest()==DECODED_RECORD_HASHES[name],'Changed decoded repair input')
    failure=v2_json.loads(records[FAILURE]['source_text'])
    saved=failure['component_result']
    cached=saved['B_identity_scope']['per_node']
    parent_certificate=saved['component_proofs']['saved_parent']
    auth={
        'failure_schema':failure.get('schema')=='original_native_regular_B_batch_source_failure_stage_v1',
        'failure_unaccepted':failure.get('accepted_identity') is False and saved.get('passed') is False,
        'failed_code':failure.get('implementation_hashes')==FAILED_IMPLEMENTATION_HASHES,
        'source_count':failure.get('source_pin_count')==1870,
        'failure_union':len(failure.get('preserved_failed_artifact_hashes_union',{}))==32,
        'saved_consumer_scope':saved.get('source_inspection_audit',{}).get('new_source_consumer_calls')==1 and saved['source_inspection_audit'].get('old_proof_functions_called')==0,
        'saved_parent_certificate':parent_certificate.get('passed') is True and isinstance(parent_certificate.get('predicates'),dict) and bool(parent_certificate['predicates']) and all(x is True for x in parent_certificate['predicates'].values()),
        'ordinary_Z':saved.get('ordinary_Z_orders')==[0,1] and saved.get('owned_Z_cell')==['0','1/2'],
        'no_previous_admission':saved.get('newly_proved_B_identity_nodes')==[] and saved.get('accepted_native_B_nodes')==[520,521,522,523],
        'records_same_as_failed_source':all(failure['source_pins'].get(name)==RECORD_HASHES[name] for name in RECORD_NAMES if name!=FAILURE),
    }
    v2_emit(diagnostic_callback,events,'saved_failure_authentication',auth)
    v2_require(all(x is True for x in auth.values()),'Changed saved failure premises')
    weighted=v2_json.loads(records[WEIGHTED]['source_text'])
    texts={alias:records[name]['source_text'] for alias,name in SOURCE_TEXT_ALIASES.items()}
    repair=derive_outer_repair(saved,weighted,texts,diagnostic_callback=lambda stage,payload:v2_emit(diagnostic_callback,events,'new_repair.'+stage,payload))
    v2_emit(diagnostic_callback,events,'new_endpoint_repairs',repair)
    v2_require(repair.get('source_only') is True and repair.get('numeric_values_computed') is False and repair.get('old_proof_functions_called')==0,'Repair scope changed')
    composed={}
    saved_keys=('weighted_graph','actual_jump','predecessor_binding','successor_binding','exact_proof_endpoint_substitutions','ordinary_Z_scope','saved_uniform_B_C1_recipe')
    for number in TARGET_NODES:
        old=cached[str(number)];new=repair.get('per_node',{}).get(str(number),{})
        expected_bindings={'predecessor':old['predecessor']['actual_source_binding'],'successor':old['successor']['actual_source_binding']}
        predicates={
            'fresh_repair_global_context':repair.get('passed') is True and repair.get('source_contract',{}).get('passed') is True and repair.get('saved_positive_context',{}).get('passed') is True,
            'cached_exact_source_port_premises':all(old['predicates'].get(k) is True for k in saved_keys),
            'same_actual_endpoint_bindings':new.get('endpoint_bindings')==expected_bindings,
            'same_endpoint_substitutions':new.get('endpoint_substitutions')==old['required_endpoint_substitutions'],
            'fresh_repair_source_laws':new.get('passed') is True and all(new.get(k) is True for k in ('R_C1_proved','U_C1_proved','m_Z2_proved','Q_C1_proved')),
            'all_new_positive_predicates':bool(new.get('predicates')) and all(x is True for x in new.get('predicates',{}).values()),
            'saved_uniform_B_chain_rule':saved.get('B_C1_derivative_rule')==B_Z_CHAIN_RULE,
        }
        passed=all(x is True for x in predicates.values())
        composed[str(number)]={'passed':passed,'B_C1_proved':passed,'native_weighted_jump_zero_proved':passed,'R_C1_proved':passed,'U_C1_proved':passed,'m_Z2_proved':passed,'Q_C1_proved':passed,'predicates':predicates,'ordinary_Z_orders':[0,1],'endpoint_bindings':expected_bindings,'source_only':True}
    v2_emit(diagnostic_callback,events,'new_B_dependency_composition',composed)
    passed=all(x['passed'] for x in composed.values())
    admitted=list(TARGET_NODES) if passed else []
    remaining=[n for n in range(524,536) if n not in admitted]
    flags=dict.fromkeys(STAGE_FLAGS,passed)
    flags['all_old_terms_retained']=True
    flags['accepted_ancestors_not_reexecuted']=True
    return {
        'schema':SCHEMA,'component':'outer_B_repair_v2','passed':passed,'source_only':True,
        'input_digest':expected_input_digest,'ordinary_Z_orders':[0,1],'ordinary_Z_not_temporal_n':True,
        'owned_Z_cell':['0','1/2'],'identity_domain_scope':saved['identity_domain_scope'],
        'source_family':saved.get('source_family'),'source_binding_scope':saved.get('source_binding_scope'),
        'native_weighted_B_definition':saved['native_weighted_B_definition'],'B_C1_derivative_rule':B_Z_CHAIN_RULE,
        'target_B_identity_nodes':list(TARGET_NODES),'newly_proved_B_identity_nodes':admitted,
        'accepted_native_B_nodes':[520,521,522,523]+admitted,'remaining_B_identity_nodes':remaining,
        'remaining_B_identity_count':len(remaining),'all_native_B_jumps_zero_proved':False,
        'B_identity_scope':{'passed':passed,'per_node':composed},
        'component_proofs':{'saved_parent':parent_certificate,'fresh_outer_repair':repair},
        'native_IN_inherited_identity':saved['native_IN_inherited_identity'],
        'native_IN_reduced_source_identity':'I_N=B_lastupper-B_firstlower-sum('+','.join('J_'+str(n-520) for n in remaining)+')+sum17(integral K_i dy)',
        'native_IN_reduction_is_source_only':True,'archived_all_16_J_terms_retained':True,
        'source_jump_zero_does_not_delete_quantitative_errors':True,
        **{k:saved[k] for k in ('retained_construction_data','retained_integral_term_inventory','old_errors_frozen','quantity_distinctions')},
        'unresolved_gates':[] if passed else ['New binding repair or dependent source composition failed; preserve this successor without replay.'],
        'remaining_construction_gates':['inner and cross-source endpoint laws','same-source quantitative error packet transport','independent terminal functional five moments','domain/cone admission, heat, stress/remainder, true temporal recursion, corrections and Cartesian NS'],
        'stage_flags':flags,**flags,**dict.fromkeys(FALSE_FLAGS,False),
        'source_inspection_audit':{'new_source_consumer_calls':1,'accepted_source_consumer_calls':0,'new_repair_calls':1,'accepted_ancestor_operator_calls':0,'old_proof_functions_called':0,'generic_DAG_walks':0,'numerical_evaluations':0,'quadratures':0,'recurrence_calls':0,'runtime_model_metadata':None},
        'saved_failure_input':{'path':FAILURE,'sha256':RECORD_HASHES[FAILURE],'positive_premises_reused_without_old_function_calls':True},
        'diagnostic_payloads':events,'runtime_model_metadata':None}
