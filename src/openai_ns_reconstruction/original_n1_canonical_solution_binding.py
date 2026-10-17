"""Source-only binding for the canonical regular original ``n=1`` solution.

The inner Volterra owner in ``experiments/root_st073`` is the definition of
the six raw fields.  This module binds that definition to the paper's regular
six-row system and records the homogeneous word estimate which makes the
regular analytic solution unique.  It reads the accepted graph and reports;
it never imports an owner, executes a builder, evaluates a recurrence, or
chooses an interval centre as a function value.

The normalized unknown is

``v=(F1/x**2, Uz1/x**2, K1/x**2, P1/x**2, F1x/x, Uz1x/x)``,

with ``x=sqrt(rho)``.  The raw regular state is consequently
``W=(rho*v1, rho*v2, rho*v3, rho*v4, x*v5, x*v6)``.  The accepted complex
source/operator packet supplies the genuine ``Msrc`` and ``beta`` intervals;
the result below exposes the theorem ``||v|| <= 2*Msrc`` under its
``beta <= 1/2`` gate, but deliberately does not recalculate that bound.

The finite forward implementation remains a finite Taylor-prefix producer.
The accepted all-order conjugacy therefore binds every available forward
coefficient to the canonical coefficients while retaining the analytic
remainder and the separate collar bridge as open scope.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_C5_terminal_primitive.json"
)
PARENT_REPORT_SHA256 = (
    "588d221fd082d5f84d94f1812eb6342bb6647847b114dac458c0165324a063b1"
)

VOLTERRA_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
)
VOLTERRA_GRAPH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.json.gz"
)
VOLTERRA_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions_check.json"
)
OPERATOR_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_operator.py"
)
OPERATOR_GRAPH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_operator.json.gz"
)
OPERATOR_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_operator_check.json"
)
FORWARD_SOURCE = "src/openai_ns_reconstruction/original_n1_core_numeric.py"
SOURCE_VECTOR_MODULE = "src/openai_ns_reconstruction/original_n1_source_vector_identity.py"
SOURCE_VECTOR_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_inlet_source_vector.json"
)
RADIAL_CONJUGACY_MODULE = "src/openai_ns_reconstruction/original_n1_radial_conjugacy.py"
GREEN_MAJORANT_MODULE = "src/openai_ns_reconstruction/original_n1_complex_radial_majorant.py"

EXPECTED_SHA256 = {
    PARENT_REPORT: PARENT_REPORT_SHA256,
    VOLTERRA_SOURCE: "ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208",
    VOLTERRA_GRAPH: "1b44be8c6cbfffb692c7a0430ecbc8a590f46001a233c99c59ca8c5e2c44a6fe",
    VOLTERRA_RECEIPT: "3cfdc3ce6a7e8c2a98a31fae985c02a2acb4e28025d313ee35db4cffc6e7b47c",
    OPERATOR_SOURCE: "268217c6b33bc717472742fa64d76a853c6760a956abfad92ac45bcf385df014",
    OPERATOR_GRAPH: "63b181010f1715d3b20589b06c1c2139ef63bcc73d7b0f3c5b7c43d9c5e03f52",
    OPERATOR_RECEIPT: "846cfedcbd0355a4a8e6a5ed445558d367cd55c7459085dd9c5d74821d969bd4",
    FORWARD_SOURCE: "a7a925360f3a6da2348c5f18ee74cfbb885f8e220730769131e43a18bd954c2f",
    SOURCE_VECTOR_MODULE: "878f2b8781215a52a5730ed8527a561ff0df451adb5da247cec4d9fe0198423f",
    SOURCE_VECTOR_REPORT: "2ac1269c3da29466067d1dfd50a92be321e85291b2634ffa6d2eafd679e9941e",
    RADIAL_CONJUGACY_MODULE: "a6b1763e172fc87a84b68ffe9c188c87f62c1aa1ba6d0b6eeeb452e1ee6c8dea",
    GREEN_MAJORANT_MODULE: "a341ba79df0d43f985bc24320a7a486889c0ab59aeff653655570f9fc7c243f8",
}

D = (0, 0, 2, 0, 3, 1)
AXIS_POWERS = (2, 2, 2, 2, 1, 1)
RAW_FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
CANONICAL_FIELDS = ("F1", "Uz1", "K1", "P1", "F1x", "Uz1x")
EXACT_DEFINITION = "W=sum_k>=0 H^k W0; W0=(I+J)Gg; H=(I+J)G(C0+B1*dZ)"
SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}


def _repo_root(root: str | Path | None = None) -> Path:
    starts: list[Path] = []
    if root is not None:
        starts.append(Path(root))
    starts.extend((Path.cwd(), Path(__file__).resolve().parents[2]))
    for start in starts:
        if start.is_file():
            start = start.parent
        for candidate in (start, *start.parents):
            if (candidate / "src").is_dir() and (candidate / "experiments" / "root_st073").is_dir():
                return candidate
    raise FileNotFoundError("st073 repository root with src/ and experiments/root_st073/ is required")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_json(path: Path) -> Mapping[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"), parse_int=str)
    if not isinstance(value, Mapping):
        raise ValueError(f"saved mapping required: {path}")
    return value


def _read_gzip_json(path: Path) -> Mapping[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        value = json.load(stream, parse_int=str)
    if not isinstance(value, Mapping):
        raise ValueError(f"saved gzip mapping required: {path}")
    return value


def _normal(text: str) -> str:
    return "".join(text.split())


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _file_pin(root: Path, relative: str) -> dict[str, Any]:
    path = root / relative
    actual = _sha256(path) if path.is_file() else None
    expected = EXPECTED_SHA256.get(relative)
    return {
        "path": relative,
        "sha256": actual,
        "expected_sha256": expected,
        "passed": bool(actual is not None and expected is not None and actual == expected),
        "source_executed": False,
        "archive_source_executed": False,
    }


def _function_node(path: Path, name: str, class_name: str | None = None) -> ast.AST:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    if class_name is None:
        candidates = tree.body
    else:
        candidates = []
        for node in tree.body:
            if isinstance(node, ast.ClassDef) and node.name == class_name:
                candidates.extend(node.body)
    for node in candidates:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
            return node
    raise ValueError(f"AST function {class_name + '.' if class_name else ''}{name} not found in {path}")


def _ast_pin(
    root: Path,
    relative: str,
    function: str,
    fragments: tuple[str, ...],
    *,
    class_name: str | None = None,
) -> dict[str, Any]:
    path = root / relative
    actual = _sha256(path)
    expected = EXPECTED_SHA256.get(relative)
    node = _function_node(path, function, class_name)
    source = ast.unparse(node)
    compact = _normal(source)
    fragment_results = {fragment: _normal(fragment) in compact for fragment in fragments}
    return {
        "path": relative,
        "function": f"{class_name + '.' if class_name else ''}{function}",
        "sha256": actual,
        "expected_sha256": expected,
        "function_ast_sha256": hashlib.sha256(
            ast.dump(node, include_attributes=False).encode("utf-8")
        ).hexdigest(),
        "fragments": fragment_results,
        "hash_passed": expected is not None and actual == expected,
        "fragments_passed": all(fragment_results.values()),
        "passed": bool(expected is not None and actual == expected and all(fragment_results.values())),
        "source_executed": False,
        "archive_source_executed": False,
    }


def _mapping_digest(value: Mapping[str, Any]) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _source_pins(root: Path) -> dict[str, Any]:
    """Pin actual definitions without importing their modules."""

    return {
        "volterra_source": _ast_pin(
            root,
            VOLTERRA_SOURCE,
            "raw_matrix_expressions",
            (
                "'4,1': eps * 4 * x * f",
                "'5,1': eps * 2 * (q + E * (-2 + de) / L)",
                "'6,1': -eps ** 2 * 8 * z * r * f / L",
                "'6,4': -eps * 4 * z / L",
            ),
        ),
        "raw_system": _ast_pin(
            root,
            VOLTERRA_SOURCE,
            "raw_system",
            (
                "Np = a.add(a.scale('-1/2', S), a.scale(parameters['Lambda'], a.add(a.scale(2, c1(Q[1])), a.scale(r, c1(Q[2])))))",
                "g4 = g.mul(c(2), eps, x, Np.value)",
                "g5 = g.mul(c(2), eps, Nt)",
                "g6 = g.sub(g.mul(c(2), eps, Nz), g.quotient(g.mul(c(4), zv, r, eps, eps, Np.value), Lv, 'positive original L'))",
                "D = D",
            ),
        ),
        "green_apply": _ast_pin(
            root,
            VOLTERRA_SOURCE,
            "green_apply",
            (
                "actual_raw_regular_n1_Green_body=True",
                "ordinary_slow_Z_derivative_order=k",
                "endpoint_independent_of_Z=True",
                "zero_endpoint_integrand_continuation=0",
                "complete_history_not_reset_at_core_seam=True",
                "signed_body_not_magnitude_cap=True",
                "kernels = (K, log, g.neg(K), g.one, g.mul(t, t, t), t)",
            ),
        ),
        "feedback_vectors": _ast_pin(
            root,
            VOLTERRA_SOURCE,
            "feedback_vectors",
            (
                "for label in ('C0', 'B1')",
                "q = derivative.c3(q.Z)",
                "vector[row] = a.add(vector[row], a.mul(B, q))",
            ),
        ),
        "owner_init": _ast_pin(
            root,
            VOLTERRA_SOURCE,
            "__init__",
            (
                "g.node('current_original_unique_regular_n1_component'",
                "zero_axis_regular_power=POWERS[i]",
                EXACT_DEFINITION,
                "all_words_retained_in_analytic_limit=True",
                "unique_regular_analytic_solution=True",
                "defining_Green_provider=parent.active.current.ast_binding(green_apply)",
                "defining_feedback_provider=parent.active.current.ast_binding(feedback_vectors)",
            ),
            class_name="CurrentOriginalN1InnerVolterraFunctions",
        ),
        "operator_matrix": _ast_pin(
            root,
            OPERATOR_SOURCE,
            "matrix",
            (
                "B0[f'{i},{j}'] = algebra.lift(sgn)",
                "B0['4,1'] = eps(f * (4 * x))",
                "B1[f'{i},{j}'] = eps(value)",
                "return dict(B0={k: truncate(v, 3) for k, v in B0.items()},",
                "D=(0, 0, 2, 0, 3, 1)",
            ),
        ),
        "forward_advance": _ast_pin(
            root,
            FORWARD_SOURCE,
            "_advance",
            (
                "k = n + 1",
                "state['f'].append([v / (4 * k * (k + 1)) for v in rf.rows[n]])",
                "state['u'].append([v / (4 * k * k) for v in ru.rows[n]])",
                "state['p'].append([v / k for v in rp.rows[n]])",
                "state['k'].append([-c.mpf(k) * v / (k + 1) for v in state['u'][-1]])",
            ),
            class_name="OriginalN1CoreNumeric",
        ),
        "forward_finite_coefficients": _ast_pin(
            root,
            FORWARD_SOURCE,
            "coefficients",
            (
                "for n in range(degree)",
                "return packet",
            ),
            class_name="OriginalN1CoreNumeric",
        ),
        "forward_finite_evaluation": _ast_pin(
            root,
            FORWARD_SOURCE,
            "evaluate_prefix",
            (
                "finite_radial_prefix_only=True",
                "infinite_radial_remainder_bounded=False",
                "full_original_n1_point_error_available=False",
            ),
            class_name="OriginalN1CoreNumeric",
        ),
    }


def _load_parent(root: Path, parent: Mapping[str, Any] | None) -> tuple[Mapping[str, Any], dict[str, Any]]:
    path = root / PARENT_REPORT
    pin = _file_pin(root, PARENT_REPORT)
    accepted = _read_json(path)
    # The saved file is authoritative.  A caller-supplied mapping is accepted
    # only as a structural witness; its bytes are never used as a replacement.
    supplied_match = parent is None
    if parent is not None:
        supplied_match = (
            _mapping(parent).get("checks", {}).get("passed") is True
            and _mapping(parent).get("n1_all_order_radial_conjugacy", {}).get("passed") is True
            and _mapping(parent).get("n1_complex_actual_source_bounds", {}).get("passed") is True
        )
    pin["supplied_mapping_shape"] = supplied_match
    pin["passed"] = bool(pin["passed"] and supplied_match)
    return accepted, pin


def _report_checks(root: Path, volterra: Mapping[str, Any], operator: Mapping[str, Any]) -> dict[str, bool]:
    vc = _mapping(volterra.get("checks"))
    ac = _mapping(vc.get("actual_source"))
    ro = _mapping(vc.get("raw_operator"))
    gf = _mapping(vc.get("Green_feedback"))
    oc = _mapping(operator.get("checks"))
    guards = _mapping(oc.get("guards_rejected"))
    return {
        "volterra_report_hash": _file_pin(root, VOLTERRA_RECEIPT)["passed"],
        "volterra_all_passed": volterra.get("all_passed") is True,
        "actual_source_report": ac.get("passed") is True,
        "same_source_operator_complete_limit": ac.get("actual_same_source_operator_bound_to_complete_limit") is True,
        "all_word_convergence_consumed": ac.get("full_raw_preconditioned_all_word_convergence_consumed") is True,
        "B1_derivative_words_retained": ac.get("full_B1_derivative_words_not_omitted") is True,
        "complete_fixed_point_rhs": ac.get("complete_fixed_point_rhs_and_separate_finite_words") is True,
        "raw_operator_report": ro.get("passed") is True,
        "raw_matrix_rows": str(ro.get("accepted_raw_matrix_AST_rows")) == "88",
        "epsilon_and_B1_retained": ro.get("all_epsilon_powers_and_B1_entries_retained") is True,
        "exact_rational_rows": ro.get("exact_rational_field_identities_without_sampling") is True,
        "green_feedback_report": gf.get("passed") is True,
        "zero_axis_regular_values": gf.get("zero_axis_regular_values") is True,
        "genuine_B0_B1_feedback": gf.get("genuine_B0_and_B1_dZ_feedback") is True,
        "operator_report_hash": _file_pin(root, OPERATOR_RECEIPT)["passed"],
        "operator_all_passed": operator.get("all_passed") is True,
        "operator_rows": str(oc.get("independent_direct_paper_operator_and_forcing_rows")) == "306",
        "ordinary_Z_derivatives": oc.get("nonconstant_amplitude_ordinary_Z_derivatives_checked") is True,
        "pressure_feedback": oc.get("pressure_gradient_feedback_checked") is True,
        "operator_guards": set(guards) == {
            "core_comparison_extension",
            "outside_common_inner_collar",
            "foreign_chart",
            "foreign_packet",
            "missing_Z_source_derivative",
        }
        and all(value is True for value in guards.values()),
        "mutation_guards": oc.get("amplitude_frame_mutation_rejected") is True
        and oc.get("live_math_packet_mutation_rejected") is True,
    }


def _node_id(value: Any, nodes: list[Any]) -> Mapping[str, Any] | None:
    try:
        index = int(value)
    except (TypeError, ValueError):
        return None
    if not 0 <= index < len(nodes):
        return None
    return _mapping(nodes[index])


def _owner_graph_binding(root: Path, graph: Mapping[str, Any]) -> dict[str, Any]:
    nodes = graph.get("exact_graph_nodes")
    bindings = _mapping(graph.get("inner_Volterra_bindings"))
    definition = _mapping(bindings.get("definition"))
    functions = _mapping(bindings.get("functions"))
    systems = _mapping(bindings.get("systems"))
    graph_family = _mapping(graph.get("source_family"))
    family_match = graph_family == SOURCE_FAMILY
    nodes_ok = isinstance(nodes, list) and len(nodes) > 0
    definition_checks = {
        "D": list(definition.get("D", ())) == [str(value) for value in D],
        "regular_axis_powers": list(definition.get("regular_axis_powers", ())) == [str(value) for value in AXIS_POWERS],
        "kinematic_entries": definition.get("kinematic_entries") == {"1,5": "1", "2,6": "1", "3,6": "-1"},
        "raw_C0_is_B0_minus_Akin": definition.get("raw_C0_is_B0_minus_Akin") is True,
        "one_core_to_collar_path": definition.get("one_core_to_collar_path_no_seam_reset") is True,
        "full_solution_distinct_from_finite_words": definition.get("full_regular_solution_distinct_from_finite_words_and_tail_enclosures") is True,
    }
    field_nodes: dict[str, list[dict[str, Any]]] = {}
    field_refs: dict[str, list[str]] = {}
    node_checks: dict[str, bool] = {}
    if nodes_ok and set(functions) == {"core", "first"}:
        for chart in ("core", "first"):
            chart_functions = _mapping(functions.get(chart))
            chart_ok = set(chart_functions) == set(RAW_FIELDS)
            for field_index, field in enumerate(RAW_FIELDS):
                refs = chart_functions.get(field)
                entries: list[dict[str, Any]] = []
                if not isinstance(refs, list) or len(refs) != 4:
                    chart_ok = False
                    refs = []
                field_refs[f"{chart}:{field}"] = [str(reference) for reference in refs]
                for reference in refs:
                    node = _node_id(reference, nodes) if nodes_ok else None
                    if node is None:
                        chart_ok = False
                        continue
                    entries.append(node)
                    expected_power = str(AXIS_POWERS[field_index])
                    node_ok = (
                        node.get("operation") == "current_original_unique_regular_n1_component"
                        and node.get("quantity") == field
                        and node.get("chart") == chart
                        and str(node.get("zero_axis_regular_power")) == expected_power
                        and node.get("exact_function_definition") == EXACT_DEFINITION
                        and node.get("all_words_retained_in_analytic_limit") is True
                        and node.get("unique_regular_analytic_solution") is True
                        and node.get("raw_units_not_normalized_potential_kernel") is True
                        and node.get("source_parameter_bounds_not_selected_as_values") is True
                        and _mapping(node.get("source_family")) == SOURCE_FAMILY
                    )
                    node_checks[f"{chart}:{field}:{node.get('ordinary_slow_Z_derivative_order')}"] = node_ok
                field_nodes[f"{chart}:{field}"] = entries
                chart_ok = chart_ok and len(entries) == 4 and all(
                    node_checks.get(f"{chart}:{field}:{node.get('ordinary_slow_Z_derivative_order')}", False)
                    for node in entries
                )
            node_checks[f"{chart}:all_fields"] = chart_ok
    else:
        node_checks["chart_field_shape"] = False

    system_checks: dict[str, bool] = {}
    for chart in ("core", "first"):
        system = _mapping(systems.get(chart))
        b0 = _mapping(system.get("B0"))
        b1 = _mapping(system.get("B1"))
        forcing = _mapping(system.get("g"))
        system_checks[f"{chart}:B0"] = set(b0) == {
            "1,5", "2,6", "3,6", "4,1", "5,1", "5,2", "5,3", "5,5",
            "6,1", "6,2", "6,3", "6,4", "6,6",
        }
        system_checks[f"{chart}:B1"] = set(b1) == {"5,1", "5,2", "5,3", "6,2", "6,3", "6,4"}
        system_checks[f"{chart}:g"] = set(forcing) == {"4", "5", "6"}
        system_checks[f"{chart}:D"] = list(system.get("D", ())) == [str(value) for value in D]
        system_checks[f"{chart}:x_rho"] = isinstance(system.get("x"), str) and isinstance(system.get("rho"), str)

    providers = []
    for node in field_nodes.values():
        providers.extend(node)
    provider_checks = {
        "Green": bool(providers)
        and all(
            _mapping(node.get("defining_Green_provider")).get("function") == "green_apply"
            and _mapping(node.get("defining_Green_provider")).get("module_sha256") == EXPECTED_SHA256[VOLTERRA_SOURCE]
            for node in providers
        ),
        "feedback": bool(providers)
        and all(
            _mapping(node.get("defining_feedback_provider")).get("function") == "feedback_vectors"
            and _mapping(node.get("defining_feedback_provider")).get("module_sha256") == EXPECTED_SHA256[VOLTERRA_SOURCE]
            for node in providers
        ),
    }
    checks = {
        "graph_hash": _file_pin(root, VOLTERRA_GRAPH)["passed"],
        "actual_raw_six_row_operator": graph.get("actual_raw_six_row_operator_and_complete_regular_limit") is True,
        "source_family": family_match,
        "definition": all(definition_checks.values()),
        "functions": all(node_checks.values()) if node_checks else False,
        "systems": all(system_checks.values()),
        "providers": all(provider_checks.values()),
    }
    return {
        "checks": {**definition_checks, **node_checks, **system_checks, **provider_checks, "passed": all(checks.values())},
        "summary_checks": {**checks, "passed": all(checks.values())},
        "definition": {
            "D": list(D),
            "regular_axis_powers": list(AXIS_POWERS),
            "kinematic_entries": {"1,5": 1, "2,6": 1, "3,6": -1},
            "raw_C0_is_B0_minus_Akin": True,
            "exact_function_definition": EXACT_DEFINITION,
        },
        "owner_namespace": "CurrentOriginalN1InnerVolterraFunctions",
        "component_operation": "current_original_unique_regular_n1_component",
        "node_bindings": field_refs,
        "field_nodes": {key: len(value) for key, value in field_nodes.items()},
        "source_family": dict(graph_family),
        "source_executed": False,
        "archive_source_executed": False,
    }


def _accepted_radial_binding(radial: Mapping[str, Any], root: Path) -> dict[str, Any]:
    symbolic = _mapping(radial.get("symbolic"))
    matrix = _mapping(symbolic.get("inherited_matrix_conjugacy", symbolic.get("matrix_conjugacy")))
    source = _mapping(symbolic.get("known_source_map", symbolic.get("inherited_source_map")))
    formal = _mapping(radial.get("formal_statement"))
    available = _mapping(radial.get("available_forward_enclosures"))
    matrix_checks = _mapping(matrix.get("checks"))
    source_checks = _mapping(source.get("checks"))
    prefix_checks = _mapping(_mapping(radial.get("accepted_receipts")).get("prefix_report", {}).get("checks"))
    source_pins = _mapping(radial.get("source_pins"))
    source_pin_ok = bool(source_pins) and all(_mapping(value).get("passed") is True for value in source_pins.values())
    checks = {
        "radial_module_hash": _file_pin(root, RADIAL_CONJUGACY_MODULE)["passed"],
        "accepted_radial_passed": radial.get("passed") is True,
        "radial_source_not_executed": radial.get("source_executed") is False,
        "accepted_receipts_passed": _mapping(radial.get("accepted_receipts")).get("passed") is True,
        "no_recurrence_replay": _mapping(radial.get("accepted_receipts")).get("accepted_recurrences_rerun") is False,
        "prefix_matrix_comparisons": str(matrix.get("state_coefficient_comparisons")) == "30",
        "prefix_mixed_comparisons": str(matrix.get("mixed_coefficient_comparisons")) == "36",
        "matrix_checks": matrix_checks.get("passed") is True,
        "source_checks": source_checks.get("all_three_sources") is True,
        "source_pin_receipts": source_pin_ok,
        "source_vector_module_hash": _file_pin(root, SOURCE_VECTOR_MODULE)["passed"],
        "source_vector_report_hash": _file_pin(root, SOURCE_VECTOR_REPORT)["passed"],
        "prefix_report_checks": prefix_checks.get("passed") is True,
        "formal_state": list(formal.get("canonical_state", ())) == list(CANONICAL_FIELDS),
        "formal_axis": formal.get("canonical_D") == [str(value) for value in D]
        and "W(0,Z)=0" in str(formal.get("axis_extension")),
        "formal_coefficients": _mapping(formal.get("canonical_coefficient_exponents")).get("row_5")
        == "4*(n+1)*(n+2)*F1_(n+1) = [x^(2n)](B0 W+B1 W_Z+g)_5"
        and _mapping(formal.get("canonical_coefficient_exponents")).get("row_6")
        == "4*(n+1)^2*Uz1_(n+1) = [x^(2n)](B0 W+B1 W_Z+g)_6"
        and _mapping(formal.get("canonical_coefficient_exponents")).get("row_4")
        == "2*(n+1)*P1_(n+1) = [x^(2n+1)](B0 W+B1 W_Z+g)_4",
        "available_prefix_degree": str(available.get("available_max_radial_degree")) == "16",
        "available_prefix_enclosures": available.get("formal_identity_applies_to_every_available_forward_coefficient") is True
        and available.get("directed_source_enclosures_present") is True
        and available.get("interval_enclosures_are_not_selected_points") is True,
        "finite_prefix_limit_retained": available.get("finite_prefix_only") is True
        and available.get("radial_remainder_certified") is False
        and available.get("canonical_full_function_equality_claimed") is False,
    }
    return {
        "checks": {**checks, "passed": all(checks.values())},
        "accepted_receipts": {
            "state_coefficient_comparisons": matrix.get("state_coefficient_comparisons"),
            "mixed_coefficient_comparisons": matrix.get("mixed_coefficient_comparisons"),
            "prefix_report_sha256": _file_pin(root, "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_numeric_prefix.json")["sha256"],
        },
        "normalization": dict(source.get("normalization", {})),
        "source_assumptions": dict(source.get("assumptions", {})),
        "source_assumptions_kind": source.get("assumptions_kind"),
        "formal_statement": dict(formal),
        "available_forward_enclosures": dict(available),
        "replayed": False,
        "source_executed": False,
        "point_values_selected": False,
    }


def _accepted_complex_bounds(
    bounds: Mapping[str, Any], green: Mapping[str, Any]
) -> dict[str, Any]:
    gates = _mapping(bounds.get("gates"))
    required_gate_names = [name for name in gates if name not in {"solution_norm", "radial_Taylor_remainder"}]
    green_input = _mapping(bounds.get("green_input"))
    green_controls = _mapping(green.get("exact_controls"))
    checks = {
        "source_bounds_packet": bounds.get("passed") is True,
        "actual_source_bounds": bounds.get("actual_source_bounds_proved") is True,
        "actual_weighted_operator_bounds": bounds.get("actual_weighted_operator_bounds_proved") is True,
        "source_only": bounds.get("source_only") is True and bounds.get("source_executed") is False,
        "no_owner_replay": bounds.get("accepted_owner_checker_fixture_recurrence_quadrature_replayed") is False,
        "diagnostic_false": bounds.get("diagnostic_only") is False,
        "missing_contract_empty": bounds.get("missing_input_contract") == [],
        "accepted_domain_gates": bool(required_gate_names) and all(gates.get(name) is True for name in required_gate_names),
        "norm_not_preclaimed": bounds.get("actual_solution_norm_proved") is False
        and bounds.get("actual_radial_Taylor_remainder_certified") is False,
        "green_packet": green.get("passed") is True
        and _mapping(green.get("actual_weighted_matrix_entries")).get("passed") is True
        and green_controls.get("passed") is True,
        "green_axis_weights": green.get("axis_weights") == [str(value) for value in AXIS_POWERS],
        "word_bound": green_input.get("word_bound") == "Msrc*(R/2)^k*(C0+k*C1/sigma)^k/k!"
        and green_input.get("geometric_beta") == "e*R/2*(C0+C1/sigma)",
        "sigma_definition": green_input.get("sigma_definition") == "Gg_Z_radius - solution_Z_radius",
        "beta_gate": green_input.get("beta_le_half") is True and gates.get("beta_le_half") is True,
        "sigma_gate": green_input.get("sigma_gap_positive") is True and gates.get("sigma_gap_positive") is True,
        "radius_contract": _mapping(bounds.get("domain")).get("rho_bound_radius") == "9/2"
        and _mapping(bounds.get("domain")).get("rho_target_radius") == "4"
        and _mapping(bounds.get("domain")).get("actual_collar_used_as_core_extension") is False,
    }
    return {
        "checks": {**checks, "passed": all(checks.values())},
        "source_bounds_digest": _mapping_digest(bounds),
        "green_packet_digest": _mapping_digest(green),
        "Msrc_interval": green_input.get("Msrc"),
        "sigma_interval": green_input.get("sigma"),
        "beta_interval": green_input.get("beta"),
        "rho_bound_radius": _mapping(bounds.get("domain")).get("rho_bound_radius"),
        "rho_target_radius": _mapping(bounds.get("domain")).get("rho_target_radius"),
        "solution_Z_radius": _mapping(bounds.get("domain")).get("solution_Z_radius"),
        "source_executed": False,
        "numeric_norm_recomputed": False,
    }


def _word_uniqueness(
    bounds: Mapping[str, Any],
    green: Mapping[str, Any],
    owner: Mapping[str, Any],
    radial: Mapping[str, Any],
) -> dict[str, Any]:
    """Separate the word majorant from analytic regular uniqueness.

    The Cauchy-loss word estimate starts with a datum in the larger ``Gg``
    tube and ends in the solution tube.  It is therefore not a contraction
    estimate on one normed space.  Uniqueness is instead proved by the radial
    triangular support of the weighted operator.
    """

    packet_ok = bounds.get("passed") is True
    beta_ok = _mapping(bounds.get("gates")).get("beta_le_half") is True
    weighted = _mapping(green.get("actual_weighted_matrix_entries"))
    entries = weighted.get("entries")
    def nonnegative_rho_degree(entry: Mapping[str, Any]) -> bool:
        try:
            return int(entry.get("after_division_by_rho_min_x_degree")) >= 0
        except (TypeError, ValueError):
            return False

    weighted_entries_ok = (
        isinstance(entries, list)
        and len(entries) == 16
        and all(
            nonnegative_rho_degree(entry)
            and entry.get("even_in_x_after_weighting") is True
            for entry in entries
            if isinstance(entry, Mapping)
        )
        and all(isinstance(entry, Mapping) for entry in entries)
    )
    exact_green = (
        green.get("word_bound") == "Msrc*(R/2)^k*(C0+k*C1/sigma)^k/k!"
        and green.get("geometric_beta") == "e*R/2*(C0+C1/sigma)"
        and green.get("factorial_proof")
        == "sum(log(j),j=1..k)>=integral_1^k log(t)dt>=k*log(k)-k"
    )
    radial_formal = _mapping(radial.get("formal_statement"))
    checks = {
        "accepted_actual_bounds": packet_ok,
        "beta_le_half": beta_ok,
        "exact_green_word_bound": exact_green,
        "nonnegative_norm_inputs_are_bound": packet_ok
        and _mapping(bounds.get("green_input")).get("sigma_gap_positive") is True,
        "weighted_A_holomorphic": weighted.get("passed") is True
        and weighted.get("identity") == "x^(p_j-q_i)*B_ij=rho*A_ij; division means analytic axis continuation"
        and weighted_entries_ok,
        "owner_word_definition": owner.get("definition", {}).get("exact_function_definition") == EXACT_DEFINITION,
        "owner_raw_C0_definition": owner.get("definition", {}).get("raw_C0_is_B0_minus_Akin") is True,
        "accepted_six_kernel_preconditioner": green.get("kinematic_preconditioning_already_in_six_kernels") is True
        and green.get("extra_I_plus_J_factor_applied") is False,
        "accepted_Z_loss_allocation": green.get("Z_loss_allocation")
        == "total gap sigma, sigma/k for each of k applications",
        "ordinary_Z_does_not_lower_radial_order": radial_formal.get("ordinary_Z_rule")
        == "D_Z^k commutes with [rho^n] for n,k>=0 when the source has the required extra rows.",
        "B1_words_are_in_owner": owner.get("summary_checks", {}).get("systems") is True,
        "H_has_rho_gain": weighted.get("identity") == "x^(p_j-q_i)*B_ij=rho*A_ij; division means analytic axis continuation",
    }
    return {
        "checks": {**checks, "passed": all(checks.values())},
        "word_bound": "||H^k Gg|| <= Msrc*(R/2)^k*(C0+k*C1/sigma)^k/k!, k>=0",
        "homogeneous_word_bound": "||H^k h||_solution <= ||h||_Gg*(R/2)^k*(C0+k*C1/sigma)^k/k! <= ||h||_Gg*beta^k, k>=1",
        "homogeneous_norms": {
            "input": "||h||_Gg on the larger source tube",
            "output": "||H^k h||_solution on the smaller solution tube",
            "same_space_contraction_claimed": False,
        },
        "coefficient_inequality": "C0+k*C1/sigma <= k*(C0+C1/sigma), k>=1, C0>=0, C1>=0, sigma>0",
        "factorial_inequality": "k^k/k! <= e^k, k>=1",
        "geometric_beta": "beta=e*R/2*(C0+C1/sigma)",
        "normal_convergence": "sum_k ||H^k Gg|| <= Msrc/(1-beta) <= 2*Msrc when beta<=1/2",
        "homogeneous_uniqueness": "h=Hh; H has a factor rho after weighting, so h_m=0 by radial induction; identity theorem in the connected rho disk gives h=0",
        "radial_triangularity": {
            "weighted_operator": "x^(p_j-q_i)*B_ij=rho*A_ij with A holomorphic",
            "order_action": "H sends rho-order m to rho-order m+1",
            "ordinary_Z_action": "D_Z commutes with [rho^m] and does not reduce m",
            "induction": "h=Hh gives h_0=0, then h_m=0 from h_{m-1}=0 for every m>=1",
            "analytic_closure": "for each Z, the zero rho Taylor series is zero on the connected rho disk",
        },
        "coefficient_uniqueness": {
            "formal_solution_unique": all(
                checks[name]
                for name in (
                    "weighted_A_holomorphic",
                    "H_has_rho_gain",
                    "ordinary_Z_does_not_lower_radial_order",
                )
            ),
            "proof": "h=Hh; coefficient induction in rho gives h_m=0 for every m, then the analytic identity theorem closes the disk",
        },
        "owner_word_definition": EXACT_DEFINITION,
        "preconditioner_binding": "The six accepted weighted Green kernels already include (I+J); no extra (I+J) factor is applied.",
        "termwise_holomorphic_derivatives": True,
        "actual_norm_interval_recomputed": False,
        "source_executed": False,
    }


def _axis_transform(owner: Mapping[str, Any], radial: Mapping[str, Any], green: Mapping[str, Any]) -> dict[str, Any]:
    owner_checks = _mapping(owner.get("checks"))
    radial_formal = _mapping(radial.get("formal_statement"))
    checks = {
        "owner_axis_powers": owner_checks.get("regular_axis_powers") is True,
        "green_axis_weights": green.get("axis_weights") == [str(value) for value in AXIS_POWERS],
        "radial_axis_extension": "W(0,Z)=0" in str(radial_formal.get("axis_extension")),
        "radial_state_order": radial_formal.get("canonical_state") == list(CANONICAL_FIELDS),
    }
    return {
        "checks": {**checks, "passed": all(checks.values())},
        "x": "sqrt(rho)",
        "rho": "x^2",
        "normalized_unknown": "v=(F1/x^2,Uz1/x^2,K1/x^2,P1/x^2,F1x/x,Uz1x/x)",
        "raw_state": "W=(rho*v1,rho*v2,rho*v3,rho*v4,x*v5,x*v6)",
        "W_components": {
            "W1": "rho*v1=F1",
            "W2": "rho*v2=Uz1",
            "W3": "rho*v3=K1",
            "W4": "rho*v4=P1",
            "W5": "x*v5=F1x",
            "W6": "x*v6=Uz1x",
        },
        "zero_axis": "W(0,Z)=0; W1..W4=O(x^2), W5,W6=O(x)",
        "axis_division_not_evaluated": True,
        "source_executed": False,
    }


def _raw_equations(owner: Mapping[str, Any], reports: Mapping[str, bool], pins: Mapping[str, Any]) -> dict[str, Any]:
    checks = {
        "raw_volterra_ast": _mapping(pins.get("volterra_source")).get("passed") is True,
        "raw_system_ast": _mapping(pins.get("raw_system")).get("passed") is True,
        "green_ast": _mapping(pins.get("green_apply")).get("passed") is True,
        "feedback_ast": _mapping(pins.get("feedback_vectors")).get("passed") is True,
        "operator_ast": _mapping(pins.get("operator_matrix")).get("passed") is True,
        "owner_raw_operator": _mapping(owner.get("summary_checks")).get("actual_raw_six_row_operator") is True,
        "accepted_volterra_rows": reports.get("raw_operator_report") is True,
        "accepted_operator_rows": reports.get("operator_report_hash") is True and reports.get("operator_rows") is True,
        "B0_B1_feedback": reports.get("genuine_B0_B1_feedback") is True,
        "zero_endpoint_and_seam": _mapping(pins.get("green_apply")).get("passed") is True,
    }
    return {
        "checks": {**checks, "passed": all(checks.values())},
        "equation": "W_x + D*W/x = B0*W + B1*D_Z W + g",
        "x": "sqrt(rho)",
        "D": list(D),
        "source_rows": {
            "g4": "2*epsilon*x*N1p",
            "g5": "2*epsilon*N1theta",
            "g6": "2*epsilon*N1z - 4*epsilon^2*Z*rho*N1p/L",
        },
        "feedback": "C0=B0 with kinematic columns removed; B1 acts on ordinary Z derivatives",
        "source_executed": False,
    }


def _failure(root: Path, reason: str, parent_pin: Mapping[str, Any] | None = None) -> dict[str, Any]:
    false_keys = {
        "passed": False,
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "failure": reason,
        "checks": {"passed": False, "failure": reason},
        "actual_canonical_owner_bound": False,
        "actual_zero_axis_condition_bound": False,
        "actual_raw_equations_bound": False,
        "actual_raw_B0_B1_g_same_source": False,
        "actual_axis_weighted_transform_bound": False,
        "actual_source_normalization_bound": False,
        "analytic_regular_uniqueness_bound": False,
        "actual_homogeneous_word_uniqueness_bound": False,
        "actual_canonical_regular_solution_exists": False,
        "actual_canonical_regular_solution_norm_proved": False,
        "source_executed": False,
        "old_pure_definitions_executed": False,
        "archive_source_executed": False,
        "accepted_owner_builder_fixture_recurrence_quadrature_replayed": False,
        "full_forward_function_identity": False,
        "q_dependent_collar_bridge": False,
    }
    if parent_pin is not None:
        false_keys["parent_report"] = dict(parent_pin)
    return false_keys


def build(
    root: str | Path | None = None,
    parent: Mapping[str, Any] | None = None,
    *,
    source_bounds: Mapping[str, Any] | None = None,
    radial_conjugacy: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return the source-bound canonical regular-solution receipt.

    ``source_bounds`` and ``radial_conjugacy`` are optional accepted mappings.
    When omitted they are read from the pinned terminal report.  Supplying a
    mapping never executes its producer; it is accepted only when the same
    source/report gates below pass.
    """

    repo = _repo_root(root)
    parent_pin: dict[str, Any] | None = None
    try:
        accepted_parent, parent_pin = _load_parent(repo, parent)
        pins = _source_pins(repo)
        volterra_report = _read_json(repo / VOLTERRA_RECEIPT)
        operator_report = _read_json(repo / OPERATOR_RECEIPT)
        volterra_graph = _read_gzip_json(repo / VOLTERRA_GRAPH)
        _operator_graph = _read_gzip_json(repo / OPERATOR_GRAPH)
        report_checks = _report_checks(repo, volterra_report, operator_report)
        owner = _owner_graph_binding(repo, volterra_graph)
        # The serialized owner graph intentionally carries the historical
        # adapter marker as false, while the accepted check report records its
        # own historical marker as true.  Neither marker is used below: the
        # actual function nodes, AST pins, raw systems, and convergence gates
        # are the acceptance evidence.
        owner["historical_graph_adapter_installed"] = volterra_graph.get(
            "current_original_n1_true_inner_Volterra_function_adapter_installed"
        )
        owner["historical_check_adapter_installed"] = volterra_report.get(
            "current_original_n1_true_inner_Volterra_function_adapter_installed"
        )
        owner["historical_adapter_marker_used_as_acceptance"] = False
        owner["actual_function_nodes_used_as_acceptance"] = owner.get("summary_checks", {}).get(
            "functions"
        ) is True
        radial = radial_conjugacy or _mapping(accepted_parent.get("n1_all_order_radial_conjugacy"))
        bounds = source_bounds or _mapping(accepted_parent.get("n1_complex_actual_source_bounds"))
        green = _mapping(accepted_parent.get("n1_weighted_complex_Green_kernel"))
        radial_binding = _accepted_radial_binding(radial, repo)
        bounds_binding = _accepted_complex_bounds(bounds, green)
        equation_binding = _raw_equations(owner, report_checks, pins)
        axis_binding = _axis_transform(owner, radial_binding, green)
        uniqueness = _word_uniqueness(bounds, green, owner, radial_binding)
        parent_checks = {
            "parent_file_hash": parent_pin.get("passed") is True,
            "parent_checks_passed": _mapping(accepted_parent.get("checks")).get("passed") is True,
            "parent_complex_bounds": _mapping(accepted_parent.get("checks", {}).get("packet_invariants")).get("actual_complex_source_operator_bounds") is True,
            "parent_no_ancestor_replay": _mapping(accepted_parent.get("checks", {}).get("packet_invariants")).get("accepted_ancestor_or_quadrature_replayed") is False,
        }
        source_pin_pass = all(_mapping(pin).get("passed") is True for pin in pins.values())
        checks = {
            "parent": all(parent_checks.values()),
            "source_ast_pins": source_pin_pass,
            "accepted_reports": all(report_checks.values()),
            "actual_canonical_owner": owner.get("summary_checks", {}).get("passed") is True,
            "actual_zero_axis": axis_binding.get("checks", {}).get("passed") is True,
            "actual_raw_equations": equation_binding.get("checks", {}).get("passed") is True,
            "actual_axis_transform": axis_binding.get("checks", {}).get("passed") is True,
            "accepted_source_normalization": radial_binding.get("checks", {}).get("source_checks") is True,
            "accepted_complex_bounds": bounds_binding.get("checks", {}).get("passed") is True,
            "homogeneous_word_estimate": uniqueness.get("checks", {}).get("passed") is True,
        }
        all_passed = all(checks.values())
        return {
            "passed": all_passed,
            "temporal_order": 1,
            "parent_report": PARENT_REPORT,
            "parent_report_sha256": PARENT_REPORT_SHA256,
            "parent_report_pin": parent_pin,
            "checks": {**checks, "passed": all_passed},
            "actual_canonical_owner_bound": checks["actual_canonical_owner"],
            "actual_zero_axis_condition_bound": checks["actual_zero_axis"],
            "actual_raw_equations_bound": checks["actual_raw_equations"],
            "actual_raw_B0_B1_g_same_source": checks["actual_raw_equations"],
            "actual_axis_weighted_transform_bound": checks["actual_axis_transform"],
            "actual_source_normalization_bound": checks["accepted_source_normalization"],
            "analytic_regular_uniqueness_bound": checks["homogeneous_word_estimate"],
            "actual_homogeneous_word_uniqueness_bound": checks["homogeneous_word_estimate"],
            "source_pins": pins,
            "accepted_reports": report_checks,
            "canonical_owner": owner,
            "raw_equations": equation_binding,
            "axis_weighted_transform": axis_binding,
            "source_normalization": radial_binding,
            "accepted_complex_source_bounds": bounds_binding,
            "homogeneous_word_estimate": uniqueness,
            "canonical_solution": {
                "exists_by_normal_convergent_word_series": all_passed,
                "definition": EXACT_DEFINITION,
                "W0_definition": "W0=preconditioned_green_apply(g)=(I+J)Gg",
                "H_definition": "H=(I+J)G(C0+B1*dZ), with the kinematic preconditioner already included in the accepted six kernels",
                "series_consumed_by_root": "W=sum_{k>=0} H^k W0",
                "kinematic_preconditioner_already_in_H": True,
                "normal_convergence_requires": "accepted beta<=1/2 interval and holomorphic common tube",
                "norm_formula": "||v|| <= Msrc/(1-beta) <= 2*Msrc",
                "norm_recomputed_here": False,
                "accepted_numeric_bounds_consumed_here": False,
                "canonical_solution_norm_proved_by_this_module": False,
                "coefficient_uniqueness": uniqueness.get("coefficient_uniqueness"),
            },
            "actual_canonical_regular_solution_exists": all_passed,
            "actual_canonical_regular_solution_norm_proved": False,
            "forward_prefix_binding": {
                "formal_identity": radial_binding.get("formal_statement"),
                "available_enclosures": radial_binding.get("available_forward_enclosures"),
                "degree16_weighted_degree15": True,
                "derivative_map": "W5/x=2*d_rho F1; W6/x=2*d_rho Uz1; degree16 induces degree15",
                "coefficient_intervals_are_not_selected_points": True,
                "finite_prefix_only": True,
            },
            "full_forward_function_identity": False,
            "q_dependent_collar_bridge": False,
            "remaining_analytic_remainder": True,
            "forward_missing_binding": [
                "an actual forward infinite-radial analytic function definition on the same complex rho/Z tube",
                "a forward radial remainder bound or forward Banach fixed-point owner",
                "the separate q-dependent collar bridge and physical field rescaling",
            ],
            "source_executed": False,
            "old_pure_definitions_executed": False,
            "archive_source_executed": False,
            "accepted_owner_builder_fixture_recurrence_quadrature_replayed": False,
            "point_values_selected": False,
            "finite_N_numerical_calculation": False,
            "tool_metadata": {
                "route": "fixed luna_max_worker",
                "initial_fork_turns": "none",
                "requested_model": "GPT-5.6 Luna/max",
                "independent_runtime_model_metadata_exposed": False,
            },
            "source_bounds_packet_digest": bounds_binding.get("source_bounds_digest"),
            "green_packet_digest": bounds_binding.get("green_packet_digest"),
        }
    except Exception as exc:  # fail closed; never turn missing evidence into a proof
        return _failure(repo, f"source binding failed closed: {type(exc).__name__}: {exc}", parent_pin)


def bind_canonical_solution(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compatibility alias for callers naming the binding operation."""

    return build(*args, **kwargs)


def prove_canonical_solution_binding(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compatibility alias; this function only reads accepted evidence."""

    return build(*args, **kwargs)


__all__ = [
    "PARENT_REPORT",
    "PARENT_REPORT_SHA256",
    "D",
    "AXIS_POWERS",
    "build",
    "bind_canonical_solution",
    "prove_canonical_solution_binding",
]
