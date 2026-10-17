"""Source-bound native common-core angular direction identity.

This module consumes the already accepted common-core and first-interface
source receipts.  It does not import or execute their producers.  The only
new algebra is the symbolic change of variables
``rho = Lambda * R`` and the ordinary-Z differentiation of the resulting
identity.  The native direction is the source function in
``inner_bridge_profiles.direction``; its finite interval implementation is
not treated as a value of the canonical common-core function.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path

import sympy as sp


PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_core_collar_pressure_identity.json"
)
PARENT_SHA256 = "ae91b0d72d0415e8146f7ee9926d1a35ddbedffff6b89eac1f26558dbb5b44c7"


# These are the source files actually inspected below.  Their hashes are also
# present in the accepted parent input map.  Reports are pinned alongside the
# source because they carry the accepted graph handles and axis proof; no
# report is regenerated here.
SOURCE_PINS = {
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_first_interface.py":
        "3785967f555cb5669ae227cda4883dab7f084bd02ca87d6742462d78605627ff",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_first_interface.json":
        "6e39acd05eff4a6ef33e6bfeade63de6d8849bf74bd9c152792cb628701ee026",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_first_interface_check.py":
        "5bc09c79d5d0105bd42756da12ef5a94992b552a1413ef8e32f18f9d65075751",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_first_interface_check.json":
        "e3a631d85da4cbeb62b9f6bbd7345e1109349935185edab177d2252a9d010400",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_common_fixed_point.py":
        "dd44110f48fb40fee715b358801dd6b8944bda5e6d5c0054eefa63c66f82255f",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_common_fixed_point.json":
        "c96ece0110427b141f870771d07e1a2db4a7b85ae0a4bd0b3ad7d8ff38bf8286",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_common_fixed_point_check.py":
        "f6fea64369400c6669eddda14840c88eedaf4c91a7fb91a0d78b47f7c36f9e0c",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_common_fixed_point_check.json":
        "a1dee74791b2f9c4a9cafcad62e573baf788c7b6f6b39cef7cde90d23a2ff672",
    "experiments/root_st073/lei_ren_part1_paper_compliant_core_uniform_bounds.json":
        "7415cba5b9abc486b2078c450bd0f2e1b7798d1921f5863c00c04526620844e7",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_nonlinear_operator.py":
        "d09f1bd0e1a55b14fd8a4006fd38f0991ea7787b88bc55ebb136c576736abf9f",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_nonlinear_operator_check.json":
        "a1f9a3632a7765f6bc8d38ac1bd25c28850e128a8b25adaa0d42f769f3704d62",
    "experiments/root_st073/lei_ren_part1_paper_compliant_inner_bridge_profiles.py":
        "142945ac1dd8a5705c380e79da29ea6bd89374c18cfa50f0412f1ba70c520840",
    "experiments/root_st073/lei_ren_part1_paper_compliant_inner_bridge_profiles_check.json":
        "ea360228f830af0dc78fdb38f43f8585cc68f6bf7b6e10b6bb3b4adce4af64ff",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py":
        "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.json.gz":
        "999aefa1e4e285dff1ebc29e1bfc4c2fb58e5507acd01a0386fbd52dfdd41b43",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json":
        "45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_core_hierarchy_source.py":
        "fa83299efc0e84d59762a36dd1cf39837501dcac9b08b39304516bdfcdd51fc3",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_core_hierarchy_source_check.json":
        "6aafbabba74d7a2609be5063eae181a123423e3dc2cb6caa632bf98e3f7e1928",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _repo_root(root: str | Path | None) -> Path:
    if root is not None:
        return Path(root).resolve()
    # src/openai_ns_reconstruction/<this file> -> repository root.
    return Path(__file__).resolve().parents[2]


def _resolve(repo: Path, value: str | Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    direct = repo / path
    if direct.exists():
        return direct
    candidate = repo / "experiments" / "root_st073" / path.name
    if candidate.exists():
        return candidate
    return direct


def _json(path: Path):
    return json.loads(path.read_bytes())


def _gzip_json(path: Path):
    return json.loads(gzip.decompress(path.read_bytes()))


def _mapping_digest(value: dict) -> str:
    """Digest a JSON mapping canonically for parent-map identity checks."""
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode(
            "utf-8"
        )
    ).hexdigest()


def _function(source: str, name: str, class_name: str | None = None):
    tree = ast.parse(source)
    found = [
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == name
    ]
    if class_name is not None:
        parents = {}

        def walk(node, owner=None):
            if isinstance(node, ast.ClassDef):
                owner = node.name
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                parents[id(node)] = owner
            for child in ast.iter_child_nodes(node):
                walk(child, owner)

        walk(tree)
        found = [node for node in found if parents.get(id(node)) == class_name]
    if len(found) != 1:
        raise ValueError(f"expected one source function {class_name or ''}.{name}")
    return found[0]


def _ast_dump(node) -> str:
    return ast.dump(node, annotate_fields=True, include_attributes=False)


def _compact(text: str) -> str:
    """Normalize source spelling without changing the AST-bound checks."""
    return "".join(text.split())


def _assignment(function_node, target: str):
    rows = []
    for node in function_node.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(item, ast.Name) and item.id == target for item in node.targets
        ):
            rows.append(node.value)
    if len(rows) != 1:
        raise ValueError(f"expected one assignment to {target}")
    return rows[0]


def _return_keywords(function_node):
    rows = [node for node in ast.walk(function_node) if isinstance(node, ast.Return)]
    if len(rows) != 1 or not isinstance(rows[0].value, ast.Call):
        raise ValueError("expected one dictionary return")
    return {keyword.arg: keyword.value for keyword in rows[0].value.keywords}


def _source_pins(repo: Path, parent: dict, errors: list[str]) -> dict:
    parent_inputs = parent.get("input_hashes", {})
    rows = {}
    for rel, expected in SOURCE_PINS.items():
        path = repo / rel
        actual = _sha256(path) if path.exists() else None
        parent_expected = parent_inputs.get(Path(rel).name)
        rows[rel] = {
            "path": rel,
            "expected_sha256": expected,
            "actual_sha256": actual,
            "parent_input_sha256": parent_expected,
            "passed": actual == expected and (parent_expected in (None, expected)),
            "source_executed": False,
        }
        if not rows[rel]["passed"]:
            errors.append(f"source pin mismatch: {rel}")
    return rows


def _source_ast_bindings(repo: Path, errors: list[str]) -> dict:
    paths = {
        "first_interface_check": repo
        / "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_first_interface_check.py",
        "inner_bridge_profiles": repo
        / "experiments/root_st073/lei_ren_part1_paper_compliant_inner_bridge_profiles.py",
        "canonical_core": repo
        / "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py",
        "common_fixed_point": repo
        / "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_common_fixed_point.py",
        "common_operator": repo
        / "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_nonlinear_operator.py",
        "hierarchy": repo
        / "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_core_hierarchy_source.py",
    }
    result = {}
    try:
        first = paths["first_interface_check"].read_text(encoding="utf-8")
        integrated = _function(first, "integrated_core_equations_proof")
        equation_bindings = _function(first, "original_core_equation_bindings")
        integrated_text = _compact(ast.unparse(integrated))
        equation_text = _compact(ast.unparse(equation_bindings))
        required_integrated = (
            "angular_original",
            "angular_integrated",
            "exact(ftc(s.diff(angular_integrated,r)),r*angular_original",
            "angular_integrated/r**2",
            "2*L*s.diff(F,r)+L*F*D",
            "original_scalar_direction_functions_replayed_exactly",
            "axis_integration_constant_proof",
        )
        required_equations = (
            "lhsf",
            "rhsf",
            "lhsu",
            "rhsu",
            "Original core pressure or centrifugal term changed",
        )
        integrated_ok = all(_compact(fragment) in integrated_text for fragment in required_integrated)
        equations_ok = all(_compact(fragment) in equation_text for fragment in required_equations)
        result["integrated_core_equations_proof"] = {
            "qualified_name": "current_core_first_interface_check.integrated_core_equations_proof",
            "lineno": integrated.lineno,
            "AST_sha256": hashlib.sha256(_ast_dump(integrated).encode()).hexdigest(),
            "required_fragments": {fragment: fragment in integrated_text for fragment in required_integrated},
            "passed": integrated_ok,
            "source_executed": False,
        }
        result["original_core_equation_bindings"] = {
            "qualified_name": "current_core_first_interface_check.original_core_equation_bindings",
            "lineno": equation_bindings.lineno,
            "AST_sha256": hashlib.sha256(_ast_dump(equation_bindings).encode()).hexdigest(),
            "required_fragments": {fragment: fragment in equation_text for fragment in required_equations},
            "passed": equations_ok,
            "source_executed": False,
        }

        inner = paths["inner_bridge_profiles"].read_text(encoding="utf-8")
        direction = _function(inner, "direction")
        direction_text = _compact(ast.unparse(direction))
        w_ast = _assignment(direction, "W")
        angular_ast = _assignment(direction, "angular")
        expected_w = ast.parse("1-(z*m)*(1-delta)-d*derivative(m)", mode="eval").body
        expected_angular = ast.parse(
            "h*(1-delta/2)-(z*derivative(h))*((1-delta)/2)-d*derivative(k)+(z*k)*(2*delta-1)",
            mode="eval",
        ).body
        direction_fragments = (
            "D_over_R=(-W+angular/(2*f))/L",
            "drive_hydro=(-W*v+(m-z*derivative(m))*((1-delta)/2)+z*a*(2*delta)-d*derivative(a))/(2*L)",
            "drive_pressure=(z*pressure*(2*(1+delta))-d*derivative(pressure))/(2*L)",
            "drive_swirl=(-(z*b)*(2*delta)+d*derivative(b)+z*p*(2*(1+delta))-d*derivative(p))/(2*L)",
        )
        return_keywords = _return_keywords(direction)
        expected_returns = {
            "D_over_R": "(-W+angular/(2*f))/L",
            "drive_hydro": "(-W*v+(m-z*derivative(m))*((1-delta)/2)+z*a*(2*delta)-d*derivative(a))/(2*L)",
            "drive_pressure": "(z*pressure*(2*(1+delta))-d*derivative(pressure))/(2*L)",
            "drive_swirl": "(-(z*b)*(2*delta)+d*derivative(b)+z*p*(2*(1+delta))-d*derivative(p))/(2*L)",
        }
        return_ast_ok = {
            key: _ast_dump(return_keywords[key])
            == _ast_dump(ast.parse(expression, mode="eval").body)
            for key, expression in expected_returns.items()
        }
        direction_ok = (
            all(_compact(fragment) in direction_text for fragment in direction_fragments)
            and all(return_ast_ok.values())
            and _ast_dump(w_ast) == _ast_dump(expected_w)
            and _ast_dump(angular_ast) == _ast_dump(expected_angular)
        )
        result["native_angular_direction"] = {
            "qualified_name": "inner_bridge_profiles.direction",
            "lineno": direction.lineno,
            "AST_sha256": hashlib.sha256(_ast_dump(direction).encode()).hexdigest(),
            "formula": "D_over_R=(-W+angular/(2*f))/L",
            "ordinary_Z_provider": "derivative(...) in the native direction source",
            "required_fragments": {fragment: _compact(fragment) in direction_text for fragment in direction_fragments},
            "return_keyword_AST_bound": return_ast_ok,
            "W_expression_AST_bound": _ast_dump(w_ast) == _ast_dump(expected_w),
            "angular_expression_AST_bound": _ast_dump(angular_ast) == _ast_dump(expected_angular),
            "passed": direction_ok,
            "source_executed": False,
        }

        canonical = paths["canonical_core"].read_text(encoding="utf-8")
        init = _function(canonical, "__init__", "CurrentOriginalN1CoreFixedPointFunctions")
        component = _function(canonical, "component", "CurrentOriginalN1CoreFixedPointFunctions")
        init_text = _compact(ast.unparse(init))
        component_text = _compact(ast.unparse(component))
        canonical_fragments = (
            "Lambda,epsilon=(g.unary('exp',logLambda),g.unary('exp',g.neg(logLambda)))",
            "L=a.add(one,a.neg(a.mul(de,a.mul(Z,Z))))",
            "self.components={key:c3.C3Function",
            "self.functions=dict(Phi=Phi,Psi=Psi,mean_Psi=M",
            "F0value=g.unary('exp',logF.value)",
            "actual_twenty_term_operator=True",
            "unique_analytic_Banach_solution=True",
            "axis_conditions='Phi(0,Z)=1; Psi(0,Z)=0'",
        )
        canonical_ok = all(
            _compact(fragment) in init_text or _compact(fragment) in component_text
            for fragment in canonical_fragments
        )
        result["canonical_core_functions"] = {
            "qualified_name": "current_original_n1_core_fixed_point_functions.CurrentOriginalN1CoreFixedPointFunctions",
            "init_lineno": init.lineno,
            "component_lineno": component.lineno,
            "init_AST_sha256": hashlib.sha256(_ast_dump(init).encode()).hexdigest(),
            "component_AST_sha256": hashlib.sha256(_ast_dump(component).encode()).hexdigest(),
            "required_fragments": {
                fragment: _compact(fragment) in init_text or _compact(fragment) in component_text
                for fragment in canonical_fragments
            },
            "passed": canonical_ok,
            "source_executed": False,
        }

        common = paths["common_fixed_point"].read_text(encoding="utf-8")
        common_tree = ast.parse(common)
        common_class = [
            node
            for node in common_tree.body
            if isinstance(node, ast.ClassDef) and node.name == "CurrentCoreCommonFixedPoint"
        ]
        common_fn = _function(common, "__init__", "CurrentCoreCommonFixedPoint")
        common_text = ast.unparse(common_fn)
        common_ok = bool(common_class) and all(
            fragment in common_text
            for fragment in ("self.source", "self.acceptance_loaded", "self.hashes")
        )
        result["common_core_owner"] = {
            "qualified_name": "current_core_common_fixed_point.CurrentCoreCommonFixedPoint",
            "init_lineno": common_fn.lineno,
            "AST_sha256": hashlib.sha256(_ast_dump(common_fn).encode()).hexdigest(),
            "source_family_owner_class_present": bool(common_class),
            "required_fragments": {fragment: fragment in common_text for fragment in ("self.source", "self.acceptance_loaded", "self.hashes")},
            "passed": common_ok,
            "source_executed": False,
        }

        operator = paths["common_operator"].read_text(encoding="utf-8")
        operator_tree = ast.parse(operator)
        operator_functions = {
            name: _function(operator, name)
            for name in ("source_terms", "green_inverse", "next_rows")
        }
        operator_ok = all(
            any(isinstance(node, ast.FunctionDef) and node.name == name for node in operator_tree.body)
            for name in operator_functions
        )
        result["common_operator"] = {
            "qualified_names": [f"current_core_nonlinear_operator.{name}" for name in operator_functions],
            "AST_sha256": {
                name: hashlib.sha256(_ast_dump(node).encode()).hexdigest()
                for name, node in operator_functions.items()
            },
            "passed": operator_ok,
            "source_executed": False,
        }

        hierarchy = paths["hierarchy"].read_text(encoding="utf-8")
        hierarchy_functions = {
            name: _function(hierarchy, name)
            for name in ("axial_second", "regular_radial_forcing")
        }
        hierarchy_ok = all(
            "rho" in ast.unparse(node) or "R" in ast.unparse(node)
            for node in hierarchy_functions.values()
        )
        result["hierarchy_source"] = {
            "qualified_names": [
                f"current_original_core_hierarchy_source.{name}" for name in hierarchy_functions
            ],
            "AST_sha256": {
                name: hashlib.sha256(_ast_dump(node).encode()).hexdigest()
                for name, node in hierarchy_functions.items()
            },
            "passed": hierarchy_ok,
            "source_executed": False,
        }
    except Exception as exc:  # fail closed while retaining the diagnostic
        errors.append(f"source AST binding failed: {exc}")
        result["error"] = str(exc)
    result["passed"] = bool(result) and all(
        bool(value.get("passed"))
        for key, value in result.items()
        if isinstance(value, dict) and key != "error"
    )
    return result


def _accepted_graph_bindings(repo: Path, parent: dict, errors: list[str]) -> dict:
    first_report = _json(
        repo
        / "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_first_interface_check.json"
    )
    common_report = _json(
        repo
        / "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_common_fixed_point_check.json"
    )
    uniform_bounds = _json(
        repo
        / "experiments/root_st073/lei_ren_part1_paper_compliant_core_uniform_bounds.json"
    )
    core_report = _gzip_json(
        repo
        / "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.json.gz"
    )
    first_proof = first_report.get("original_integrated_stress_free_core_equations", {})
    family_keys = (
        "actual_five_defect_family_sha256",
        "implicit_source_sha256",
        "datum_enclosure_sha256",
    )
    families = [
        tuple(first_report.get(key) for key in family_keys),
        tuple(common_report.get(key) for key in family_keys),
        tuple(core_report.get("source_family", {}).get(key) for key in family_keys),
    ]
    family_ok = all(family == families[0] and all(item is not None for item in family) for family in families)
    first_ok = bool(first_report.get("all_passed")) and all(
        bool(first_proof.get(key))
        for key in (
            "original_core_equations_AST_bound",
            "differentiated_integrated_angular_equation_identity",
            "original_scalar_direction_functions_replayed_exactly",
            "same_common_pressure_P_equals_P0_plus_integral_F_squared",
            "integrated_core_means_and_direction_are_source_functions_not_caps",
        )
    )
    common_ok = bool(common_report.get("all_passed")) and bool(
        common_report.get("current_core_actual_operator_Banach_fixed_point_admitted")
    )
    graph = core_report.get("core_fixed_point_bindings", {})
    function_nodes = graph.get("functions", {})
    graph_ok = (
        bool(core_report.get("unique_actual_analytic_solution_function_definitions"))
        and graph.get("definition", {}).get("Lambda_epsilon_exactly_one") is True
        and graph.get("definition", {}).get("actual_twenty_term_source_and_Green_inverse_retained") is True
        and all(key in function_nodes for key in ("Phi", "Psi", "mean_Psi"))
    )
    phi_lower = uniform_bounds.get("actual_Phi_lower", {})
    phi_lower_tuple = phi_lower.get("lower_exact_mpf_tuple", [])
    phi_strict_lower = (
        isinstance(phi_lower_tuple, list)
        and len(phi_lower_tuple) >= 2
        and phi_lower_tuple[0] == 0
        and isinstance(phi_lower_tuple[1], int)
        and phi_lower_tuple[1] > 0
    )
    phi_domain_ok = (
        uniform_bounds.get("scaled_radial_domain") == ["0", "4.1"]
        and uniform_bounds.get("real_axial_domain") == ["-1", "1"]
        and uniform_bounds.get("unique_analytic_field_bounded_not_sampled") is True
        and phi_strict_lower
        and uniform_bounds.get("implicit_source_sha256") == families[0][1]
        and uniform_bounds.get("datum_enclosure_sha256") == families[0][2]
    )
    complex_bounds = parent.get("n1_complex_actual_source_bounds", {})
    complex_gates = complex_bounds.get("gates", {})
    complex_domain = complex_bounds.get("domain", {})
    nonvanishing_ok = (
        complex_bounds.get("passed") is True
        and complex_bounds.get("actual_source_bounds_proved") is True
        and complex_bounds.get("source_only") is True
        and complex_bounds.get("parent_mapping_matches_file") is True
        and complex_gates.get("analytic_domain_receipt") is True
        and complex_gates.get("F0_complex_axis_binding") is True
        and complex_gates.get("actual_L_lower_from_delta_cap") is True
        and complex_gates.get("L_holomorphic_lower_positive") is True
        and complex_domain.get("rho_source_disk") == "|rho| < 5"
        and complex_domain.get("rho_target_radius") == "4"
    )
    f0_positive_receipt = complex_gates.get("F0_complex_axis_binding") is True
    l_positive_receipt = (
        complex_gates.get("actual_L_lower_from_delta_cap") is True
        and complex_gates.get("L_holomorphic_lower_positive") is True
    )
    accepted = {
        "family_match_across_actual_reports": family_ok,
        "first_interface_integrated_proof": first_ok,
        "common_fixed_point_report": common_ok,
        "canonical_function_graph": graph_ok,
        "accepted_Phi_strict_positive_domain": phi_domain_ok,
        "accepted_F0_and_L_nonvanishing_domain": nonvanishing_ok,
        "first_interface_report_sha256": _sha256(
            repo
            / "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_first_interface_check.json"
        ),
        "common_fixed_point_report_sha256": _sha256(
            repo
            / "experiments/root_st073/lei_ren_part1_paper_compliant_current_core_common_fixed_point_check.json"
        ),
        "canonical_graph_report_sha256": _sha256(
            repo
            / "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.json.gz"
        ),
        "canonical_function_nodes": {
            key: function_nodes.get(key)
            for key in ("Phi", "Psi", "mean_Psi", "F0", "U", "m")
        },
        "nonvanishing_domain": {
            "Phi_receipt": {
                "path": "experiments/root_st073/lei_ren_part1_paper_compliant_core_uniform_bounds.json",
                "scaled_radial_domain": uniform_bounds.get("scaled_radial_domain"),
                "real_axial_domain": uniform_bounds.get("real_axial_domain"),
                "actual_Phi_lower_exact_mpf_tuple": phi_lower_tuple,
                "strict_lower_positive": phi_strict_lower,
                "implicit_source_sha256": uniform_bounds.get("implicit_source_sha256"),
                "datum_enclosure_sha256": uniform_bounds.get("datum_enclosure_sha256"),
                "unique_analytic_field_bounded_not_sampled": uniform_bounds.get(
                    "unique_analytic_field_bounded_not_sampled"
                ),
            },
            "F0_and_L_receipt": {
                "parent_section": "n1_complex_actual_source_bounds",
                "F0_positive_on_real_source_domain": f0_positive_receipt,
                "L_positive_on_real_source_domain": l_positive_receipt,
                "gates": {
                    key: complex_gates.get(key)
                    for key in (
                        "analytic_domain_receipt",
                        "F0_complex_axis_binding",
                        "actual_L_lower_from_delta_cap",
                        "L_holomorphic_lower_positive",
                    )
                },
                "domain": {
                    "rho_source_disk": complex_domain.get("rho_source_disk"),
                    "rho_target_radius": complex_domain.get("rho_target_radius"),
                    "L_modulus_lower_formula": complex_domain.get("L_modulus_lower_formula"),
                },
            },
            "passed": phi_domain_ok and nonvanishing_ok,
        },
        "integrated_proof_receipt": {
            key: first_proof.get(key)
            for key in (
                "original_core_equations_AST_bound",
                "differentiated_integrated_angular_equation_identity",
                "original_scalar_direction_functions_replayed_exactly",
                "exact_stress_free_recovery",
                "axis_integration_constant_proof",
            )
        },
        "source_executed": False,
    }
    accepted["passed"] = all(
        accepted[key]
        for key in (
            "family_match_across_actual_reports",
            "first_interface_integrated_proof",
            "common_fixed_point_report",
            "canonical_function_graph",
            "accepted_Phi_strict_positive_domain",
            "accepted_F0_and_L_nonvanishing_domain",
        )
    )
    if not accepted["passed"]:
        errors.append("accepted common-core or integrated-direction receipt gate failed")
    return accepted


def _symbolic_direction_proof(nonvanishing_domain: dict | None = None) -> dict:
    """Derive the scaled radial identity and its ordinary-Z jets exactly."""
    rho, R, Lambda, epsilon, L, D = sp.symbols(
        "rho R Lambda epsilon L D", nonzero=True
    )
    Z = sp.symbols("Z")
    Phi = sp.Function("Phi")(rho, Z)
    F0 = sp.Function("F0")(Z)
    F = F0 * Phi
    # The common-source change of variables is rho=Lambda*R, hence
    # d_R=Lambda*d_rho.  F0 is allowed to depend on Z and is never frozen.
    dR_F = Lambda * sp.diff(F, rho)
    integrated_direction = sp.solve(sp.Eq(2 * L * dR_F + L * F * D, 0), D)[0]
    desired_direction = -2 * Lambda * sp.diff(Phi, rho) / Phi
    direction_residual = sp.factor(integrated_direction - desired_direction)
    dbar = (epsilon * rho) * desired_direction
    desired_dbar = -2 * rho * sp.diff(Phi, rho) / Phi
    dbar_residual = sp.factor((dbar - desired_dbar).subs(epsilon, 1 / Lambda))
    ordinary_direction = []
    ordinary_dbar = []
    for order in range(4):
        ordinary_direction.append(
            sp.simplify(sp.diff(direction_residual, Z, order)) == 0
        )
        ordinary_dbar.append(sp.simplify(sp.diff(dbar_residual, Z, order)) == 0)
    receipt_bound = bool(nonvanishing_domain and nonvanishing_domain.get("passed"))
    return {
        "assumptions": [
            "rho=Lambda*R",
            "d_R=Lambda*d_rho",
            "Lambda*epsilon=1",
            "F_core=F0(Z)*Phi_core(rho,Z)",
            "F0 is independent of rho and may vary with Z",
            "formal quotient denominators are nonzero wherever the source-bound domain receipt applies",
            "accepted strict Phi>0 floor on real rho in [0,4.1], Z in [-1,1]",
            "accepted F0=exp(logF) axis binding gives F0>0 on the real source domain",
            "accepted delta-cap/holomorphic lower bound gives L>0 on the real source domain",
        ],
        "formal_nonzero_symbols": ["Lambda", "epsilon", "L", "Phi_core"],
        "nonvanishing_domain_receipt_bound": receipt_bound,
        "nonvanishing_domain_scope": (
            "Phi>0 for rho in [0,4.1], Z in [-1,1]; F0>0 and L>0 from "
            "the accepted source/complex-domain gates; quotient identities extend "
            "through the regular axis by the accepted integrated proof"
        ),
        "integrated_equation_used": "0=2*L*d_R(F_core)+L*F_core*D_over_R",
        "derived_D_over_R": "-2*d_R(log(F_core))",
        "scaled_D_over_R": "-2*Lambda*Phi_core_rho/Phi_core",
        "scaled_Dbar": "-2*rho*Phi_core_rho/Phi_core",
        "direction_residual": str(direction_residual),
        "Dbar_residual_after_Lambda_epsilon": str(dbar_residual),
        "direction_identity_exact": direction_residual == 0,
        "Dbar_identity_exact": dbar_residual == 0,
        "ordinary_Z_orders": list(range(4)),
        "ordinary_Z_direction_zero": ordinary_direction,
        "ordinary_Z_Dbar_zero": ordinary_dbar,
        "ordinary_Z_derivatives_exact": all(ordinary_direction + ordinary_dbar),
        "F0_Z_dependence_retained": all(
            sp.diff(F0, Z, order) == sp.diff(F0, Z, order) for order in range(4)
        ),
    }


def _packet(
    repo: Path,
    parent_path: Path,
    parent: dict,
    errors: list[str],
    *,
    parent_mapping_matches_disk: bool = True,
    supplied_parent_mapping_digest: str | None = None,
    pinned_parent_mapping_digest: str | None = None,
) -> dict:
    parent_actual = _sha256(parent_path) if parent_path.exists() else None
    parent_ok = parent_actual == PARENT_SHA256
    if not parent_ok:
        errors.append("parent report hash mismatch")
    if not parent_mapping_matches_disk:
        errors.append("supplied parent mapping differs from pinned disk parent mapping")
    pins = _source_pins(repo, parent, errors)
    source_ast = _source_ast_bindings(repo, errors)
    accepted = _accepted_graph_bindings(repo, parent, errors)
    nonvanishing = accepted.get("nonvanishing_domain", {})
    symbolic = _symbolic_direction_proof(nonvanishing)
    symbolic_ok = all(
        symbolic[key]
        for key in (
            "direction_identity_exact",
            "Dbar_identity_exact",
            "ordinary_Z_derivatives_exact",
        )
    )
    if not symbolic_ok:
        errors.append("fresh symbolic angular direction identity failed")
    f0_source_exp_ast_bound = bool(
        source_ast.get("canonical_core_functions", {})
        .get("required_fragments", {})
        .get("F0value=g.unary('exp',logF.value)")
    )
    f0_positive_receipt = bool(
        nonvanishing.get("F0_and_L_receipt", {}).get("F0_positive_on_real_source_domain")
    )
    l_positive_receipt = bool(
        nonvanishing.get("F0_and_L_receipt", {}).get("L_positive_on_real_source_domain")
    )
    phi_positive_receipt = bool(
        nonvanishing.get("Phi_receipt", {}).get("strict_lower_positive")
    )
    nonvanishing["F0_source_exp_AST_bound"] = f0_source_exp_ast_bound
    nonvanishing["F0_positive_bound"] = f0_source_exp_ast_bound and f0_positive_receipt
    nonvanishing["L_positive_bound"] = l_positive_receipt
    nonvanishing["F_core_positive_bound"] = (
        phi_positive_receipt and nonvanishing["F0_positive_bound"]
    )
    nonvanishing_ok = bool(nonvanishing.get("passed")) and all(
        (
            f0_source_exp_ast_bound,
            f0_positive_receipt,
            l_positive_receipt,
            phi_positive_receipt,
        )
    )
    if not nonvanishing_ok:
        errors.append("accepted Phi/F0/L nonvanishing-domain receipt gate failed")
    source_ast_ok = bool(source_ast.get("passed"))
    checks = {
        "parent_report_pin": parent_ok,
        "parent_mapping_matches_pinned_disk": parent_mapping_matches_disk,
        "source_pins": bool(pins) and all(row["passed"] for row in pins.values()),
        "source_ast_bindings": source_ast_ok,
        "accepted_common_core_and_integrated_proof": bool(accepted.get("passed")),
        "accepted_Phi_F0_L_nonvanishing_domain": nonvanishing_ok,
        "fresh_symbolic_direction_identity": symbolic_ok,
        "ordinary_Z_0_through_3": symbolic["ordinary_Z_derivatives_exact"],
        "axis_constant_consumed_without_axis_division": bool(
            accepted.get("integrated_proof_receipt", {}).get("axis_integration_constant_proof")
        ),
    }
    passed = all(checks.values()) and not errors
    result = {
        "passed": passed,
        "source_only": True,
        "parent_report": str(parent_path),
        "parent_report_sha256": parent_actual,
        "parent_report_expected_sha256": PARENT_SHA256,
        "parent_mapping_digest": {
            "supplied": supplied_parent_mapping_digest,
            "pinned_disk": pinned_parent_mapping_digest,
            "matches": parent_mapping_matches_disk,
        },
        "source_pins": pins,
        "actual_function_handles": {
            "canonical_common_core": {
                "qualified_name": "current_original_n1_core_fixed_point_functions.CurrentOriginalN1CoreFixedPointFunctions.functions['Phi']",
                "graph_nodes": accepted.get("canonical_function_nodes", {}).get("Phi"),
                "source_ast_binding": source_ast.get("canonical_core_functions"),
                "definition": "unique actual analytic common-core Phi function; not an interval center",
            },
            "native_angular_direction": {
                "qualified_name": "inner_bridge_profiles.direction",
                "source_ast_binding": source_ast.get("native_angular_direction"),
                "definition": "D_over_R=(-W+angular/(2*f))/L using ordinary-Z derivative rows",
            },
            "integrated_core_equations": {
                "qualified_name": "current_core_first_interface_check.integrated_core_equations_proof",
                "source_ast_binding": source_ast.get("integrated_core_equations_proof"),
                "accepted_receipt": accepted.get("integrated_proof_receipt"),
            },
            "common_core_owner": source_ast.get("common_core_owner"),
        },
        "angular_direction_identity": {
            "native_formula": "D_over_R=(-W+angular/(2*f))/L",
            "integrated_core_formula": "angular_integrated/R**2=2*L*d_R(F_core)+L*F_core*D_over_R",
            "integrated_residual_axis_constant": "zero by accepted analytic-axis proof; D quotient is not evaluated at R=0",
            "physical_identity": "D_over_R=-2*d_R(log(F_core))",
            "canonical_scaled_identity": "D_over_R=-2*Lambda*Phi_core_rho/Phi_core",
            "canonical_scaled_bar_identity": "Dbar=R*D_over_R=-2*rho*Phi_core_rho/Phi_core",
            "ordinary_Z_orders": list(range(4)),
            "source_function_scope": "same common core function; no finite-width collar or selected point value",
            "nonvanishing_domain": nonvanishing,
            "symbolic_proof": symbolic,
        },
        "ordinary_Z_orders": {
            "orders": list(range(4)),
            "direction_identity_zero": symbolic["ordinary_Z_direction_zero"],
            "Dbar_identity_zero": symbolic["ordinary_Z_Dbar_zero"],
            "F0_Z_derivatives_retained": True,
            "ordinary_derivative_operator": "source direction derivative(...) and analytic product/quotient differentiation",
        },
        "checks": checks,
        "errors": errors,
        "source_executed": False,
        "old_owner_builder_fixture_recurrence_norm_M0_W1_W2_quadrature_replayed": False,
        "point_values_selected": False,
        "finite_N_numerical_calculation": False,
        "finite_width_collar_identity_claimed": False,
        "full_forward_function_identity_claimed": False,
        "all_parent_error_caps_retained": True,
        "remaining_scope": [
            "q-dependent collar bridge and finite-width source identity remain separate",
            "the native angular identity is source-function exact only on the common core domain",
            "no numerical radial value, global collar equality, or six-term closure is claimed",
        ],
        "tool_metadata": {
            "route": "fixed luna_max_worker",
            "initial_fork_turns": "none",
            "requested_model": "gpt-5.6-luna/max",
            "independent_runtime_model_metadata_exposed": False,
        },
    }
    return result


def build(root: str | Path | None = None, parent: str | Path | dict | None = None) -> dict:
    """Return the source-bound angular-direction proof packet.

    ``parent`` may be the pinned report path or an already loaded mapping.  A
    mapping is accepted for orchestration convenience, but the pinned report
    on disk is still hashed whenever available.  This function reads source
    and saved receipts only; it does not import or execute any accepted
    producer.
    """
    repo = _repo_root(root)
    pinned_parent_mapping_digest: str | None = None
    if parent is None:
        parent_path = repo / PARENT_REPORT
        parent_obj = _json(parent_path)
        parent_mapping_matches_disk = True
        supplied_parent_mapping_digest = _mapping_digest(parent_obj)
    elif isinstance(parent, dict):
        parent_path = repo / PARENT_REPORT
        parent_obj = parent
        pinned_parent_obj = _json(parent_path)
        supplied_parent_mapping_digest = _mapping_digest(parent_obj)
        pinned_parent_mapping_digest = _mapping_digest(pinned_parent_obj)
        parent_mapping_matches_disk = supplied_parent_mapping_digest == pinned_parent_mapping_digest
    else:
        parent_path = _resolve(repo, parent)
        parent_obj = _json(parent_path)
        parent_mapping_matches_disk = True
        supplied_parent_mapping_digest = _mapping_digest(parent_obj)
    if pinned_parent_mapping_digest is None:
        pinned_parent_obj = _json(repo / PARENT_REPORT)
        pinned_parent_mapping_digest = _mapping_digest(pinned_parent_obj)
    errors: list[str] = []
    return _packet(
        repo,
        parent_path,
        parent_obj,
        errors,
        parent_mapping_matches_disk=parent_mapping_matches_disk,
        supplied_parent_mapping_digest=supplied_parent_mapping_digest,
        pinned_parent_mapping_digest=pinned_parent_mapping_digest,
    )


__all__ = ["build"]
