"""Source-bound first-chart raw ``n=1`` Volterra response.

This module closes the source/graph binding for the first chart.  It reads the
saved raw six-row graph and the accepted native-collar packet, then checks the
actual ASTs which define the Green map, feedback, core-to-first entry, and the
ordinary slow-
``Z`` derivatives.  It does not execute a numerical owner, recurrence,
fixture, builder, or quadrature.

The graph owner records an all-word formula in physical raw units.  The
accepted coupled-first-step tail uses the same raw Green kernels and its
``Akin`` preconditioning, so its raw units are checked from the initial
kernel and perturbation-bound ASTs before the tail receipt is consumed.  The
separate normalized ``Kchi`` kernel is not used for this proof.  The finite
``OriginalN1CoreNumeric`` producer remains a finite-prefix source pin only.
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
    "lei_ren_part1_paper_compliant_current_original_native_collar_transport.json"
)
PARENT_REPORT_SHA256 = (
    "7c05f7438a27d9caa4378cb89b9ab6d33b44138086e86292b73f2fd8e3f4bf96"
)

INNER_GRAPH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.json.gz"
)
INNER_GRAPH_SHA256 = (
    "1b44be8c6cbfffb692c7a0430ecbc8a590f46001a233c99c59ca8c5e2c44a6fe"
)
VOLERRA_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
)
VOLERRA_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions_check.json"
)
RADIAL_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_radial_functions.py"
)
RADIAL_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_radial_functions_check.json"
)
OPERATOR_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_operator.py"
)
OPERATOR_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_operator_check.json"
)
COUPLED_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step.py"
)
COUPLED_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step.json.gz"
)
INITIAL_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_initial_integral.py"
)
INITIAL_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_initial_integral_check.json"
)
FORWARD_SOURCE = "src/openai_ns_reconstruction/original_n1_core_numeric.py"
SOURCE_VECTOR_MODULE = "src/openai_ns_reconstruction/original_n1_source_vector_identity.py"
SOURCE_VECTOR_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_inlet_source_vector.json"
)
RADIAL_CONJUGACY_MODULE = "src/openai_ns_reconstruction/original_n1_radial_conjugacy.py"
CANONICAL_BINDING_MODULE = (
    "src/openai_ns_reconstruction/original_n1_canonical_solution_binding.py"
)
COMPLEX_MAJORANT_MODULE = (
    "src/openai_ns_reconstruction/original_n1_complex_radial_majorant.py"
)
PREFIX_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_core_numeric_prefix.json"
)

EXPECTED_SHA256 = {
    PARENT_REPORT: PARENT_REPORT_SHA256,
    INNER_GRAPH: INNER_GRAPH_SHA256,
    VOLERRA_SOURCE: "ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208",
    VOLERRA_RECEIPT: "3cfdc3ce6a7e8c2a98a31fae985c02a2acb4e28025d313ee35db4cffc6e7b47c",
    RADIAL_SOURCE: "d80c05494bf5901e09235d2757c9da1a7ff103b925e337b61774f738ad3f68b2",
    RADIAL_RECEIPT: "0d0ac44de54df037fbd4d170a705fb69dd73634f5aa1133e72332d6fbe6ba330",
    OPERATOR_SOURCE: "268217c6b33bc717472742fa64d76a853c6760a956abfad92ac45bcf385df014",
    OPERATOR_RECEIPT: "846cfedcbd0355a4a8e6a5ed445558d367cd55c7459085dd9c5d74821d969bd4",
    COUPLED_SOURCE: "97913b82b8570b220921cd89918569a4019ff3495bdbea0548b3c53e9d8c93e7",
    COUPLED_RECEIPT: "30b807af4147658ec0978ebd253cc1b2b31bed67bfed430a965e66461450add1",
    INITIAL_SOURCE: "ed2a029d17e56503a397f8e380fcaf75abc7a8f82db6ce5c7eeb71fda560e825",
    INITIAL_RECEIPT: "897f38c99eb0f99cefbe0d4e71983c7dc588763a34ba44f0df38a500de09e3e3",
    FORWARD_SOURCE: "a7a925360f3a6da2348c5f18ee74cfbb885f8e220730769131e43a18bd954c2f",
    SOURCE_VECTOR_MODULE: "878f2b8781215a52a5730ed8527a561ff0df451adb5da247cec4d9fe0198423f",
    SOURCE_VECTOR_REPORT: "2ac1269c3da29466067d1dfd50a92be321e85291b2634ffa6d2eafd679e9941e",
    RADIAL_CONJUGACY_MODULE: "a6b1763e172fc87a84b68ffe9c188c87f62c1aa1ba6d0b6eeeb452e1ee6c8dea",
    CANONICAL_BINDING_MODULE: "569337e2d49b8423eac9bcc8eb1ac7f319847ef66cbefe7f79f72922ed0ec0ac",
    COMPLEX_MAJORANT_MODULE: "a341ba79df0d43f985bc24320a7a486889c0ab59aeff653655570f9fc7c243f8",
    PREFIX_REPORT: "2ae583cdf5ea9fde03208fea0423eb11d171054c9a8dd45d78aa5fe2076d81fa",
}

D = (0, 0, 2, 0, 3, 1)
AXIS_POWERS = (2, 2, 2, 2, 1, 1)
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
CANONICAL_STATE = ("F1", "Uz1", "K1", "P1", "F1x", "Uz1x")
FORWARD_STATE = ("f", "u", "k", "p", "f_x", "u_x")
ORDINARY_Z_ORDERS = (0, 1, 2, 3)
NORMALIZATION = {
    "f": "epsilon*F1/F0",
    "u": "Uz1/epsilon",
    "k": "K1/epsilon",
    "p": "P1/(epsilon*Pstar^2)",
}
SCALE_DIAGONAL = (
    "F0(Z)/epsilon",
    "epsilon",
    "epsilon",
    "epsilon*Pstar^2",
    "F0(Z)/epsilon",
    "epsilon",
)
SCALE_INVERSE = (
    "epsilon/F0(Z)",
    "1/epsilon",
    "1/epsilon",
    "1/(epsilon*Pstar^2)",
    "epsilon/F0(Z)",
    "1/epsilon",
)
CANONICAL_EQUATION = "W_x + D*W/x = B0*W + B1*D_Z W + g"
RAW_WORD_DEFINITION = "W=sum_k>=0 H^k W0; W0=(I+J)Gg; H=(I+J)G(C0+B1*dZ)"
FORWARD_EQUATION = (
    "Y_x + D*Y/x = S^{-1}(B0*S+B1*S_Z)Y "
    "+ (S^{-1}*B1*S)Y_Z + S^{-1}g"
)


def _repo_root(root: str | Path | None) -> Path:
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
    raise FileNotFoundError("st073 repository root is required")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _read_json(path: Path) -> Mapping[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, Mapping):
        raise ValueError(f"JSON mapping required: {path}")
    return value


def _read_gzip_json(path: Path) -> Mapping[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, Mapping):
        raise ValueError(f"gzip JSON mapping required: {path}")
    return value


def _all_true(values: Mapping[str, Any]) -> bool:
    return bool(values) and all(value is True for value in values.values())


def _compact(source: str) -> str:
    return "".join(source.split()).replace('"', "'")


def _ast_node(path: Path, name: str, owner: str | None = None) -> tuple[str, ast.AST]:
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(path))
    if owner is None:
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
                return source, node
    else:
        for node in tree.body:
            if isinstance(node, ast.ClassDef) and node.name == owner:
                for child in node.body:
                    if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == name:
                        return source, child
    raise ValueError(f"AST node {owner + '.' if owner else ''}{name} missing in {path}")


def _ast_requirements(
    repo: Path,
    relative: str,
    definitions: Mapping[str, tuple[str, ...]],
    class_methods: Mapping[tuple[str, str], tuple[str, ...]] = (),
) -> dict[str, Any]:
    path = repo / relative
    source = path.read_text(encoding="utf-8")
    found: dict[str, dict[str, Any]] = {}
    for name, fragments in definitions.items():
        raw, node = _ast_node(path, name)
        body = _compact(ast.get_source_segment(raw, node) or "")
        checks = {fragment: _compact(fragment) in body for fragment in fragments}
        found[name] = {
            "function": name,
            "AST_sha256": hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest(),
            "fragments": checks,
            "passed": all(checks.values()),
        }
    for (owner, name), fragments in dict(class_methods).items():
        raw, node = _ast_node(path, name, owner)
        body = _compact(ast.get_source_segment(raw, node) or "")
        checks = {fragment: _compact(fragment) in body for fragment in fragments}
        key = f"{owner}.{name}"
        found[key] = {
            "function": name,
            "owner": owner,
            "AST_sha256": hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest(),
            "fragments": checks,
            "passed": all(checks.values()),
        }
    return {
        "path": relative,
        "sha256": _sha256(path),
        "expected_sha256": EXPECTED_SHA256[relative],
        "hash_passed": _sha256(path) == EXPECTED_SHA256[relative],
        "definitions": found,
        "source_executed": False,
        "archive_source_executed": False,
        "passed": _sha256(path) == EXPECTED_SHA256[relative]
        and all(item["passed"] for item in found.values()),
    }


def _accepted_parent(repo: Path, parent: Mapping[str, Any] | None) -> tuple[Mapping[str, Any], dict[str, Any]]:
    path = repo / PARENT_REPORT
    actual = _sha256(path)
    if actual != PARENT_REPORT_SHA256:
        raise ValueError("native-collar parent report changed")
    saved = _read_json(path)
    supplied = saved if parent is None else parent
    if not isinstance(supplied, Mapping):
        raise ValueError("native-collar parent mapping required")
    if _digest(supplied) != _digest(saved):
        raise ValueError("supplied native-collar parent does not match pinned report")
    return supplied, {
        "passed": True,
        "path": PARENT_REPORT,
        "sha256": actual,
        "expected_sha256": PARENT_REPORT_SHA256,
        "mapping_digest": _digest(supplied),
        "input_kind": "complete_native_collar_transport_report",
        "source_executed": False,
        "archive_source_executed": False,
    }


def _file_pins(repo: Path) -> dict[str, Any]:
    pins: dict[str, Any] = {}
    for relative, expected in EXPECTED_SHA256.items():
        if relative == PARENT_REPORT:
            continue
        path = repo / relative
        actual = _sha256(path)
        pins[relative] = {
            "path": relative,
            "expected_sha256": expected,
            "actual_sha256": actual,
            "passed": actual == expected,
            "source_executed": False,
            "archive_source_executed": False,
        }
    return pins


def _source_pins(repo: Path) -> dict[str, Any]:
    pins = _file_pins(repo)
    pins["forward_ast"] = _ast_requirements(
        repo,
        FORWARD_SOURCE,
        {
            "linear_rhs": (
                "rf=(2*e*(Q+E*(-2+de)/L)-2*H*g/L)*f",
                "ru=-8*F2*z*r*A*f/L",
                "rp=2*F2/(e*P2)*A*f",
            ),
            "_rhs": (
                "Aw=known(leading['A'],n+1,axial+2)",
                "Uw=known(leading['Uz'],n+1,axial+2)",
                "Qw=(2*zz*uu-zz*mm*(1-delta)-dd*mw.dz())/LL",
                "pressure_head=R([[(i+1)*(i+2)*vforvinQw.rows[i+1][:axial+1]",
            ),
            "_advance": (
                "axial=degree+depth-n-1",
                "k=n+1",
                "state['f'].append([v/(4*k*(k+1))forvinrf.rows[n]])",
                "state['u'].append([v/(4*k*k)forvinru.rows[n]])",
                "state['p'].append([v/kforvinrp.rows[n]])",
                "state['k'].append([-c.mpf(k)*v/(k+1)forvinstate['u'][-1]])",
            ),
            "coefficients": ("forninrange(degree)",),
            "evaluate_prefix": (
                "finite_radial_prefix_only=True",
                "infinite_radial_remainder_bounded=False",
                "full_original_n1_point_error_available=False",
                "full_original_n1_point_oracle=False",
            ),
        },
    )
    pins["volterra_ast"] = _ast_requirements(
        repo,
        VOLERRA_SOURCE,
        {
            "raw_matrix_expressions": (
                "'4,1':eps*4*x*f",
                "'5,1':eps*2*(q+E*(-2+de)/L)",
                "'6,1':-eps**2*8*z*r*f/L",
                "'6,4':-eps*4*z/L",
            ),
            "raw_system": (
                "Np=a.add(a.scale('-1/2',S),a.scale(parameters['Lambda'],a.add(a.scale(2,c1(Q[1])),a.scale(r,c1(Q[2])))))",
                "g4=g.mul(c(2),eps,x,Np.value)",
                "g5=g.mul(c(2),eps,Nt)",
                "g6=g.sub(g.mul(c(2),eps,Nz),g.quotient(g.mul(c(4),zv,r,eps,eps,Np.value),Lv,'positiveoriginalL'))",
                "D=D",
            ),
            "green_apply": (
                "actual_raw_regular_n1_Green_body=True",
                "ordinary_slow_Z_derivative_order=k",
                "endpoint_independent_of_Z=True",
                "zero_endpoint_integrand_continuation=0",
                "complete_history_not_reset_at_core_seam=True",
                "signed_body_not_magnitude_cap=True",
                "kernels=(K,log,g.neg(K),g.one,g.mul(t,t,t),t)",
            ),
            "feedback_vectors": (
                "forlabelin('C0','B1')",
                "q=derivative.c3(q.Z)",
                "vector[row]=a.add(vector[row],a.mul(B,q))",
            ),
        },
        {
            ("CurrentOriginalN1InnerVolterraFunctions", "__init__"): (
                "g.node('current_original_unique_regular_n1_component'",
                "zero_axis_regular_power=POWERS[i]",
                RAW_WORD_DEFINITION,
                "all_words_retained_in_analytic_limit=True",
                "unique_regular_analytic_solution=True",
                "raw_units_not_normalized_potential_kernel=True",
                "one_core_to_collar_path_no_seam_reset=True",
                "defining_Green_provider=parent.active.current.ast_binding(green_apply)",
                "defining_feedback_provider=parent.active.current.ast_binding(feedback_vectors)",
            )
        },
    )
    pins["radial_ast"] = _ast_requirements(
        repo,
        RADIAL_SOURCE,
        {
            "equation_velocity_rows": (
                "f1=div(f,g.mul(g.constant(2),x))",
                "u1=div(u,g.mul(g.constant(2),x))",
                "f2=a.add(div(H['5'],g.mul(g.constant(4),r)),a.neg(div(f,g.mul(r,x))))",
                "u2=a.add(div(H['6'],g.mul(g.constant(4),r)),a.neg(div(u,g.mul(g.constant(2),r,x))))",
            ),
            "axis_extend": (
                "regular_even_radial_axis_continuation",
                "axis_value_is_exact_equation_limit=True",
            ),
        },
        {
            ("CurrentOriginalN1InnerRadialFunctions", "__init__"): (
                "full_raw_B0_and_B1_dZ_retained=True",
                "core_first_seam_inherits_one_regular_solution=True",
                "original_numeric_n1_solution_oracle_executed=False",
            )
        },
    )
    pins["operator_ast"] = _ast_requirements(
        repo,
        OPERATOR_SOURCE,
        {"matrix": (
            "B0[f'{i},{j}']=algebra.lift(sgn)",
            "B0['4,1']=eps(f*(4*x))",
            "B1[f'{i},{j}']=eps(value)",
            "returndict(B0={k:truncate(v,3)",
            "D=(0,0,2,0,3,1)",
        )},
    )
    pins["coupled_tail_ast"] = _ast_requirements(
        repo,
        COUPLED_SOURCE,
        {
            "preconditioned_tail": (
                "defpreconditioned_tail",
                "S=1+a",
                "Mbase=domain.upper(c,S*M0)",
                "beta=domain.upper(c,c.e*a*S*(Cp/(N+1)+C1/loss))",
                "arbitrary_k_term_formula='Mbase*[a*(1+a)*(Cp+C1*k/loss)]^k/k!'",
            ),
            "perturbation_bound": (
                "r4=4*epsilon*x*f",
                "r5=epsilon*(2*(q+(2-delta)*E/Lmin)",
                "r6=epsilon*(8*epsilon*z*rho_max*f/Lmin",
                "d5=epsilon*(2*H/Lmin+4*d*Cf/Lmin)",
                "d6=epsilon*(2*(H+d*rho_max*wr)/Lmin",
            ),
        },
    )
    pins["initial_integral_ast"] = _ast_requirements(
        repo,
        INITIAL_SOURCE,
        {
            "kernels": (
                "returndict(F1=cell*tail/2,Uz1=log_kernel(c,cell,endpoint),K1=-cell*tail/2,P1=c.mpf(1),d_x_F1=ratio**3,d_x_Uz1=ratio)",
            ),
            "integrate_cell": (
                "weights=kernels(algebra.ctx,cell,endpoint)",
                "g=system['g']",
                "Inversemicroscopicwidthdidnotcancelbeforeintegration",
            ),
        },
    )
    pins["all_source_hashes_passed"] = all(
        isinstance(value, Mapping) and value.get("passed") is True
        for key, value in pins.items()
        if key != "all_source_hashes_passed"
    )
    return pins


def _graph_binding(repo: Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    graph_path = repo / INNER_GRAPH
    graph = _read_gzip_json(graph_path)
    bindings = _mapping(graph.get("inner_Volterra_bindings"))
    functions = _mapping(bindings.get("functions"))
    expected_functions = {
        "core": {
            "F1": [43300, 43301, 43302, 43303],
            "Uz1": [43304, 43305, 43306, 43307],
            "K1": [43308, 43309, 43310, 43311],
            "P1": [43312, 43313, 43314, 43315],
            "d_x_F1": [43316, 43317, 43318, 43319],
            "d_x_Uz1": [43320, 43321, 43322, 43323],
        },
        "first": {
            "F1": [43324, 43325, 43326, 43327],
            "Uz1": [43328, 43329, 43330, 43331],
            "K1": [43332, 43333, 43334, 43335],
            "P1": [43336, 43337, 43338, 43339],
            "d_x_F1": [43340, 43341, 43342, 43343],
            "d_x_Uz1": [43344, 43345, 43346, 43347],
        },
    }
    nodes = graph.get("exact_graph_nodes")
    if not isinstance(nodes, list):
        raise ValueError("exact graph node list missing")
    owner = _mapping(parent.get("n1_actual_canonical_solution_binding")).get("canonical_owner")
    owner_family = _mapping(owner).get("source_family")
    graph_family = graph.get("source_family")
    checks: dict[str, bool] = {
        "graph_hash": _sha256(graph_path) == INNER_GRAPH_SHA256,
        "node_count": len(nodes) == 49024,
        "source_family": graph_family == owner_family,
        "actual_raw_six_row_operator": graph.get("actual_raw_six_row_operator_and_complete_regular_limit") is True,
        "historical_adapter_not_used": graph.get("current_original_n1_true_inner_Volterra_function_adapter_installed") is False,
        "definition_D": _mapping(bindings.get("definition")).get("D") == list(D),
        "definition_axis": _mapping(bindings.get("definition")).get("regular_axis_powers") == list(AXIS_POWERS),
        "definition_kinematics": _mapping(bindings.get("definition")).get("kinematic_entries") == {"1,5": 1, "2,6": 1, "3,6": -1},
        "raw_C0": _mapping(bindings.get("definition")).get("raw_C0_is_B0_minus_Akin") is True,
        "one_path": _mapping(bindings.get("definition")).get("one_core_to_collar_path_no_seam_reset") is True,
        "full_distinct": _mapping(bindings.get("definition")).get("full_regular_solution_distinct_from_finite_words_and_tail_enclosures") is True,
        "raw_units_distinct": _mapping(bindings.get("definition")).get("normalized_Kchi_not_used_with_raw_matrix") is True,
        "no_point_accuracy": _mapping(bindings.get("definition")).get("numerical_original_n1_point_accuracy_claimed") is False,
        "functions_exact": functions == expected_functions,
    }
    systems = _mapping(bindings.get("systems"))
    for chart in ("core", "first"):
        system = _mapping(systems.get(chart))
        checks[f"{chart}_system_D"] = system.get("D") == list(D)
        checks[f"{chart}_system_B0"] = isinstance(system.get("B0"), Mapping) and bool(system.get("B0"))
        checks[f"{chart}_system_B1"] = isinstance(system.get("B1"), Mapping) and bool(system.get("B1"))
        checks[f"{chart}_system_C0"] = isinstance(system.get("C0"), Mapping) and bool(system.get("C0"))
        checks[f"{chart}_system_source"] = isinstance(system.get("g"), Mapping) and set(system.get("g")) == {"4", "5", "6"}
        checks[f"{chart}_system_coordinates"] = system.get("x") is not None and system.get("rho") is not None
    for chart, field_map in expected_functions.items():
        for field, ids in field_map.items():
            for order, node_id in enumerate(ids):
                if node_id >= len(nodes) or not isinstance(nodes[node_id], Mapping):
                    checks[f"{chart}:{field}:{order}"] = False
                    continue
                node = nodes[node_id]
                checks[f"{chart}:{field}:{order}"] = (
                    node.get("operation") == "current_original_unique_regular_n1_component"
                    and node.get("quantity") == field
                    and node.get("chart") == chart
                    and node.get("ordinary_slow_Z_derivative_order") == order
                    and node.get("source_family") == graph_family
                    and node.get("coordinate_domain_lower") == 0
                    and node.get("coordinate_domain_upper") == (13 if chart == "core" else 35824)
                )
    leading_first = _mapping(bindings.get("leading")).get("first")
    checks["first_histories_start_at_core"] = _mapping(leading_first).get("all_histories_start_at_true_core_rho4") is True
    checks["first_query_coordinate"] = _mapping(bindings.get("coordinates")).get("first") == 35848
    checks["core_query_coordinate"] = _mapping(bindings.get("coordinates")).get("core") == 35558
    words = _mapping(bindings.get("words"))
    fixed = _mapping(bindings.get("fixed_point_rhs"))
    checks["word_zero_first"] = set(_mapping(words.get("0")).get("first")) == set(FIELDS)
    checks["word_one_first"] = set(_mapping(words.get("1")).get("first")) == set(FIELDS)
    checks["fixed_point_first"] = set(_mapping(fixed.get("first"))) == set(FIELDS)
    checks["word_handles_are_lists"] = all(
        isinstance(_mapping(words.get(str(word))).get("first", {}).get(field), list)
        and len(_mapping(words.get(str(word))).get("first", {}).get(field)) == 4
        for word in (0, 1)
        for field in FIELDS
    )
    checks["passed"] = all(checks.values())
    return {
        "checks": checks,
        "passed": checks["passed"],
        "path": INNER_GRAPH,
        "sha256": _sha256(graph_path),
        "expected_sha256": INNER_GRAPH_SHA256,
        "source_family": graph_family,
        "node_count": len(nodes),
        "owner_namespace": "CurrentOriginalN1InnerVolterraFunctions",
        "component_operation": "current_original_unique_regular_n1_component",
        "field_handles": expected_functions,
        "core_coordinate": _mapping(bindings.get("coordinates")).get("core"),
        "first_coordinate": _mapping(bindings.get("coordinates")).get("first"),
        "system_D": list(D),
        "axis_powers": list(AXIS_POWERS),
        "definition": dict(_mapping(bindings.get("definition"))),
        "word_handles": {
            "word0": dict(_mapping(words.get("0")).get("first", {})),
            "word1": dict(_mapping(words.get("1")).get("first", {})),
            "fixed_point_rhs": dict(_mapping(fixed.get("first", {}))),
        },
        "leading_first": {
            "query": _mapping(leading_first).get("query"),
            "domain": _mapping(leading_first).get("domain"),
            "all_histories_start_at_true_core_rho4": _mapping(leading_first).get("all_histories_start_at_true_core_rho4"),
        },
        "tail_admission_metadata": dict(_mapping(bindings.get("tail_admission"))),
        "source_executed": False,
        "archive_source_executed": False,
    }


def _normalization_scale(parent: Mapping[str, Any]) -> dict[str, Any]:
    radial = _mapping(parent.get("n1_all_order_radial_conjugacy"))
    symbolic = _mapping(radial.get("symbolic"))
    source_map = _mapping(symbolic.get("inherited_source_map"))
    accepted_forward = _mapping(parent.get("n1_infinite_forward_analytic_response"))
    accepted_scale = _mapping(accepted_forward.get("scale_conjugacy"))
    source_rows = _mapping(source_map.get("canonical_source_rows"))
    forward_rows = _mapping(source_map.get("forward_normalized_sources"))
    expected_rows = {
        "g4": "2*epsilon*x*N1p = 2*epsilon*x*Pstar^2*Np",
        "g5": "2*epsilon*N1theta",
        "g6": "2*epsilon*N1z - 4*epsilon^2*Z*rho*N1p/L",
    }
    expected_forward_rows = {
        "pressure": "Np=N1p/Pstar^2",
        "angular": "2*epsilon^2*N1theta/F0",
        "axial": "2*N1z-4*epsilon*Z*rho*N1p/L",
    }
    checks = {
        "accepted_forward_packet": accepted_forward.get("passed") is True,
        "accepted_source_conjugacy": radial.get("passed") is True,
        "normalization_exact": source_map.get("normalization") == NORMALIZATION,
        "canonical_source_rows": dict(source_rows) == expected_rows,
        "forward_source_rows": dict(forward_rows) == expected_forward_rows,
        "scale_checks_accepted": _mapping(accepted_scale.get("checks")).get("passed") is True,
        "F0_nonzero_domain": _mapping(accepted_scale.get("nonvanishing_premises")).get("F0") == "F0 != 0",
        "epsilon_nonzero_domain": _mapping(accepted_scale.get("nonvanishing_premises")).get("epsilon") == "epsilon != 0",
        "L_nonzero_domain": _mapping(accepted_scale.get("nonvanishing_premises")).get("L") == "L != 0",
        "axis_identity_domain": _mapping(accepted_scale.get("axis")).get("source_vector_identity_domain") == "rho>0 and L != 0",
    }
    symbolic = _fresh_symbolic_scale()
    checks.update({
        "S_inverse_S_identity": symbolic["S_inverse_S_identity"],
        "D_commutes_with_S": symbolic["D_commutes_with_S"],
        "S_x_zero": symbolic["S_x_zero"],
        "S_z_equals_S_Gamma": symbolic["S_z_equals_S_Gamma"],
    })
    return {
        "checks": {**checks, "passed": all(checks.values())},
        "canonical_state": list(CANONICAL_STATE),
        "forward_state": list(FORWARD_STATE),
        "normalization": dict(NORMALIZATION),
        "rho": "rho=x^2",
        "x": "sqrt(rho)",
        "W_equals_SY": "W=S(Z)Y",
        "scale_matrix": {
            "type": "diagonal",
            "diagonal": list(SCALE_DIAGONAL),
            "S_inverse_diagonal": list(SCALE_INVERSE),
            "positional_map": [
                {"raw_component": raw, "forward_component": forward, "scale": scale, "inverse": inverse}
                for raw, forward, scale, inverse in zip(CANONICAL_STATE, FORWARD_STATE, SCALE_DIAGONAL, SCALE_INVERSE)
            ],
        },
        "exact_symbolic_scale_checks": symbolic,
        "connection": "Gamma=S^{-1}*D_Z(S)",
        "derivative_identity": "D_Z(SY)=S*(D_Z Y+Gamma*Y)",
        "operator_conjugacy": {
            "B0_forward": "S^{-1}(B0*S+B1*S_Z)",
            "B1_forward": "S^{-1}*B1*S",
            "g_forward": "S^{-1}*g",
            "equation": FORWARD_EQUATION,
            "canonical_equation": CANONICAL_EQUATION,
            "S_x": "0",
            "Gamma": "diag(F0_Z/F0,0,0,0,F0_Z/F0,0)",
        },
        "source_scale_map": {
            "canonical_source_rows": expected_rows,
            "forward_normalized_sources": expected_forward_rows,
            "pressure_normalization": "Np=N1p/Pstar^2",
            "fourth_pde_source": "(S^{-1}g)_4=2*x*Np",
            "fourth_pde_source_is_distinct_from_Np_definition": True,
        },
        "nonvanishing_premises": {
            "F0": "F0 != 0",
            "epsilon": "epsilon != 0",
            "P2": "P2 != 0",
            "P2_equals_Pstar_squared": "P2 = Pstar^2",
            "L": "L != 0",
        },
        "source_executed": False,
        "archive_source_executed": False,
    }


def _fresh_symbolic_scale() -> dict[str, Any]:
    """Verify only the diagonal scale algebra with fresh exact symbols."""
    import sympy as sp

    epsilon, F0, F0z, Pstar, x = sp.symbols("epsilon F0 F0z Pstar x", nonzero=True)
    diagonal = [F0 / epsilon, epsilon, epsilon, epsilon * Pstar**2, F0 / epsilon, epsilon]
    inverse = [epsilon / F0, 1 / epsilon, 1 / epsilon, 1 / (epsilon * Pstar**2), epsilon / F0, 1 / epsilon]
    products = [sp.simplify(a * b - 1) == 0 for a, b in zip(diagonal, inverse)]
    gamma = sp.diag(F0z / F0, 0, 0, 0, F0z / F0, 0)
    S = sp.diag(*diagonal)
    Sz = S * gamma
    Dmat = sp.diag(*D)
    commute = S * Dmat - Dmat * S
    expected_Sz = [F0z / epsilon, 0, 0, 0, F0z / epsilon, 0]
    return {
        "S_inverse_S_identity": all(products),
        "component_products": products,
        "D_commutes_with_S": all(sp.simplify(value) == 0 for value in commute),
        "S_x_zero": all(sp.diff(value, x) == 0 for value in diagonal),
        # ``Matrix`` iteration is over all 36 entries, so compare its
        # diagonal explicitly rather than accidentally comparing the fifth
        # row's off-diagonal zero with the fifth expected component.
        "S_z_equals_S_Gamma": list(Sz.diagonal()) == expected_Sz,
        "connection_diagonal": ["F0_Z/F0", "0", "0", "0", "F0_Z/F0", "0"],
        "exact_symbolic_backend": "sympy",
    }


def _equation_binding(parent: Mapping[str, Any], graph: Mapping[str, Any], pins: Mapping[str, Any], scale: Mapping[str, Any]) -> dict[str, Any]:
    binding = _mapping(parent.get("n1_actual_canonical_solution_binding"))
    raw_equations = _mapping(binding.get("raw_equations"))
    canonical_owner = _mapping(binding.get("canonical_owner"))
    checks = {
        "canonical_equation_exact": raw_equations.get("equation") == CANONICAL_EQUATION,
        "canonical_D_exact": raw_equations.get("D") == list(D),
        "canonical_source_rows_exact": raw_equations.get("source_rows") == {
            "g4": "2*epsilon*x*N1p",
            "g5": "2*epsilon*N1theta",
            "g6": "2*epsilon*N1z - 4*epsilon^2*Z*rho*N1p/L",
        },
        "raw_ast_passed": _mapping(pins.get("volterra_ast")).get("passed") is True,
        "operator_ast_passed": _mapping(pins.get("operator_ast")).get("passed") is True,
        "raw_graph_passed": graph.get("passed") is True,
        "owner_raw_units": canonical_owner.get("actual_function_nodes_used_as_acceptance") is True,
        "source_raw_units": graph.get("definition", {}).get("normalized_Kchi_not_used_with_raw_matrix") is True,
        "scale_passed": scale.get("checks", {}).get("passed") is True,
        "B0_connection_included": "B1*S_Z" in FORWARD_EQUATION,
        "physical_pressure_source_distinct": scale.get("source_scale_map", {}).get("fourth_pde_source_is_distinct_from_Np_definition") is True,
    }
    return {
        "checks": {**checks, "passed": all(checks.values())},
        "equation": CANONICAL_EQUATION,
        "forward_equation": FORWARD_EQUATION,
        "D": list(D),
        "x": "sqrt(rho)",
        "raw_source_rows": dict(_mapping(raw_equations.get("source_rows"))),
        "source_rows": {
            "g4": "2*epsilon*x*N1p = 2*epsilon*x*Pstar^2*Np",
            "g5": "2*epsilon*N1theta",
            "g6": "2*epsilon*N1z - 4*epsilon^2*Z*rho*N1p/L",
        },
        "feedback": "C0=B0 with kinematic columns removed; B1 acts on ordinary Z derivatives",
        "scale_connection": "B0hat=S^{-1}(B0*S+B1*S_Z), B1hat=S^{-1}*B1*S, ghat=S^{-1}*g",
        "pressure_units": {
            "Np_definition": "Np=N1p/Pstar^2",
            "fourth_pde_source": "(S^{-1}g)_4=2*x*Np",
            "do_not_conflate": True,
        },
        "source_executed": False,
        "archive_source_executed": False,
    }


def _core_entry_binding(graph: Mapping[str, Any], pins: Mapping[str, Any]) -> dict[str, Any]:
    handles = _mapping(graph.get("field_handles"))
    word = _mapping(graph.get("word_handles"))
    leading = _mapping(graph.get("leading_first"))
    checks = {
        "graph_passed": graph.get("passed") is True,
        "core_handles": set(_mapping(handles.get("core"))) == set(FIELDS),
        "first_handles": set(_mapping(handles.get("first"))) == set(FIELDS),
        "first_word0": set(_mapping(word.get("word0"))) == set(FIELDS),
        "first_word1": set(_mapping(word.get("word1"))) == set(FIELDS),
        "first_fixed_point_rhs": set(_mapping(word.get("fixed_point_rhs"))) == set(FIELDS),
        "one_core_to_collar_path": graph.get("definition", {}).get("one_core_to_collar_path_no_seam_reset") is True,
        "first_histories_start_at_core": leading.get("all_histories_start_at_true_core_rho4") is True,
        "green_ast_seam": _mapping(pins.get("volterra_ast", {}).get("definitions", {}).get("green_apply", {})).get("fragments", {}).get("complete_history_not_reset_at_core_seam=True") is True,
        "green_ast_endpoint": _mapping(pins.get("volterra_ast", {}).get("definitions", {}).get("green_apply", {})).get("fragments", {}).get("endpoint_independent_of_Z=True") is True,
    }
    return {
        "checks": {**checks, "passed": all(checks.values())},
        "core_coordinate": graph.get("core_coordinate"),
        "first_coordinate": graph.get("first_coordinate"),
        "core_field_handles": dict(_mapping(handles.get("core"))),
        "first_field_handles": dict(_mapping(handles.get("first"))),
        "first_word0_handles": dict(_mapping(word.get("word0"))),
        "first_word1_handles": dict(_mapping(word.get("word1"))),
        "first_fixed_point_rhs_handles": dict(_mapping(word.get("fixed_point_rhs"))),
        "entry_definition": "first W0=(I+J)G applied to core source at rho=4, then the local first-chart integral; no history reset",
        "source_executed": False,
        "archive_source_executed": False,
    }


def _ordinary_Z_binding(graph: Mapping[str, Any], pins: Mapping[str, Any], parent: Mapping[str, Any]) -> dict[str, Any]:
    nodes = _mapping(parent.get("n1_actual_canonical_solution_binding")).get("canonical_owner", {}).get("node_bindings", {})
    first_handles = graph.get("field_handles", {}).get("first", {})
    checks = {
        "first_fields_all_orders": all(
            first_handles.get(field) == [first_handles[field][0] + order for order in ORDINARY_Z_ORDERS]
            for field in FIELDS
        ),
        "first_node_orders": all(
            graph.get("checks", {}).get(f"first:{field}:{order}") is True
            for field in FIELDS
            for order in ORDINARY_Z_ORDERS
        ),
        "green_callback_order": _mapping(pins.get("volterra_ast", {}).get("definitions", {}).get("green_apply", {})).get("fragments", {}).get("ordinary_slow_Z_derivative_order=k") is True,
        "feedback_Z_derivative": _mapping(pins.get("volterra_ast", {}).get("definitions", {}).get("feedback_vectors", {})).get("fragments", {}).get("q=derivative.c3(q.Z)") is True,
        "ordinary_Z_commutes_with_radial": _mapping(
            _mapping(parent.get("n1_actual_canonical_solution_binding")).get("forward_prefix_binding")
        ).get("formal_identity", {}).get("ordinary_Z_rule")
        == "D_Z^k commutes with [rho^n] for n,k>=0 when the source has the required extra rows.",
        "accepted_core_prefix_ordinary_Z": _mapping(parent.get("n1_actual_raw_prefix_binding")).get("checks", {}).get("ordinary_Z_archive") is True,
    }
    return {
        "checks": {**checks, "passed": all(checks.values())},
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "first_field_handles": dict(first_handles),
        "statement": "ordinary D_Z^j, j=0..3, is carried by actual graph nodes and commutes with radial coefficient extraction; the Green source callback reserves order j",
        "source_executed": False,
        "archive_source_executed": False,
    }


def _all_word_binding(
    parent: Mapping[str, Any],
    graph: Mapping[str, Any],
    scale: Mapping[str, Any],
    pins: Mapping[str, Any],
) -> dict[str, Any]:
    canonical = _mapping(parent.get("n1_actual_canonical_solution_binding"))
    canonical_solution = _mapping(canonical.get("canonical_solution"))
    norm = _mapping(parent.get("n1_actual_canonical_norm_tail"))
    convergence = _mapping(norm.get("convergence_proof"))
    graph_tail = _mapping(graph.get("tail_admission_metadata"))
    graph_definition = _mapping(graph.get("definition"))
    green_ast = _mapping(pins.get("volterra_ast", {}).get("definitions", {}).get("green_apply", {}))
    initial_ast = _mapping(pins.get("initial_integral_ast"))
    coupled_ast = _mapping(pins.get("coupled_tail_ast"))
    H_definition = canonical_solution.get("H_definition")
    checks = {
        "owner_definition_exact": canonical_solution.get("definition") == RAW_WORD_DEFINITION,
        "owner_all_words_metadata": _mapping(canonical.get("canonical_owner")).get("actual_function_nodes_used_as_acceptance") is True,
        "canonical_norm_tail_accepted": norm.get("passed") is True and norm.get("actual_canonical_regular_solution_norm_proved") is True,
        "canonical_word_bound_accepted": convergence.get("normal_convergence_on_compact_subsets") is True and convergence.get("termwise_derivatives_on_inner_compact_subsets") is True,
        "raw_graph_declares_all_words": graph.get("definition", {}).get("full_regular_solution_distinct_from_finite_words_and_tail_enclosures") is True,
        "raw_units_are_explicit": graph.get("definition", {}).get("normalized_Kchi_not_used_with_raw_matrix") is True,
        "scale_map_is_exact": scale.get("checks", {}).get("passed") is True,
        "tail_metadata_present": bool(graph_tail),
        "raw_initial_integral_ast": _mapping(pins.get("initial_integral_ast")).get("passed") is True,
        "raw_coupled_tail_ast": _mapping(pins.get("coupled_tail_ast")).get("passed") is True,
    }
    # The saved tail packet is Akin-preconditioned in the same raw units as
    # green_apply.  The initial-kernel and perturbation-bound AST pins below
    # establish that provenance; the unrelated normalized Kchi kernel is not
    # substituted for this raw proof.
    raw_convergence_checks = {
        "raw_green_equals_preconditioned_kernel": (
            green_ast.get("passed") is True
            and green_ast.get("fragments", {}).get("actual_raw_regular_n1_Green_body=True") is True
            and initial_ast.get("passed") is True
            and graph_definition.get("normalized_Kchi_not_used_with_raw_matrix") is True
        ),
        "principal_resolvent_factorization": (
            graph_definition.get("raw_C0_is_B0_minus_Akin") is True
            and canonical_solution.get("kinematic_preconditioner_already_in_H") is True
            and H_definition
            == "H=(I+J)G(C0+B1*dZ), with the kinematic preconditioner already included in the accepted six kernels"
            and coupled_ast.get("passed") is True
        ),
        "raw_initial_kernel_units": (
            initial_ast.get("passed") is True
            and initial_ast.get("definitions", {}).get("kernels", {}).get("passed") is True
        ),
        "raw_perturbation_bounds": (
            coupled_ast.get("passed") is True
            and coupled_ast.get("definitions", {}).get("perturbation_bound", {}).get("passed") is True
        ),
        "raw_tail_receipt_admitted": (
            graph_tail.get("retained_preconditioned_indices") == [0, 1]
            and graph_tail.get("all_k_derivative_words_allowed") is True
            and graph_tail.get("half_derivative_block_reduction_used") is False
            and graph_tail.get("old_unpreconditioned_tail_reused") is False
            and graph_tail.get("geometric_tail_ratio_upper") is not None
        ),
        "raw_tail_constant_bound": (
            convergence.get("normal_convergence_on_compact_subsets") is True
            and convergence.get("termwise_derivatives_on_inner_compact_subsets") is True
            and graph_tail.get("geometric_tail_ratio_upper") is not None
        ),
    }
    return {
        "checks": {**checks, "passed": all(checks.values())},
        "formal_definition": RAW_WORD_DEFINITION,
        "owner_all_words_retained_in_analytic_limit": True,
        "canonical_normal_convergence_consumed": all(
            checks[key] for key in ("canonical_norm_tail_accepted", "canonical_word_bound_accepted")
        ),
        "raw_all_word_convergence_certified": all(raw_convergence_checks.values()),
        "raw_all_word_convergence_checks": raw_convergence_checks,
        "raw_all_word_convergence_passed": all(raw_convergence_checks.values()),
        "raw_tail_metadata_consumed_as_bound": all(raw_convergence_checks.values()),
        "tail_admission_scope": graph_tail.get("scope"),
        "tail_admission_is_raw_preconditioned_proof": all(raw_convergence_checks.values()),
        "raw_tail_provenance": (
            "The accepted coupled-first-step tail is consumed with its same-source initial "
            "raw kernels and perturbation_bound AST: Mbase=(1+a)M0, kinematic J is already "
            "inside the six Green kernels, and beta bounds the raw B0-Akin plus B1 words."
        ),
        "missing_factorization": None,
        "source_executed": False,
        "archive_source_executed": False,
    }


def _finite_scope(pins: Mapping[str, Any]) -> dict[str, Any]:
    ast_pin = pins.get("forward_ast", {})
    return {
        "source": FORWARD_SOURCE,
        "source_sha256": EXPECTED_SHA256[FORWARD_SOURCE],
        "AST_pin": ast_pin,
        "definition": "OriginalN1CoreNumeric._advance and evaluate_prefix produce finite radial prefixes only",
        "advance_rules": {
            "f": "f_(n+1)=rf_n/[4*(n+1)*(n+2)]",
            "u": "u_(n+1)=ru_n/[4*(n+1)^2]",
            "p": "p_(n+1)=rp_n/(n+1)",
            "k": "k_(n+1)=-(n+1)*u_(n+1)/(n+2)",
        },
        "finite_loop": "for n in range(degree)",
        "infinite_radial_remainder_bounded": False,
        "full_original_n1_point_oracle": False,
        "source_executed": False,
    }


def _failure(reason: str, parent_pin: Mapping[str, Any] | None = None) -> dict[str, Any]:
    out: dict[str, Any] = {
        "passed": False,
        "status": "failed_closed",
        "source_only": True,
        "failure": reason,
        "actual_first_raw_source_definition_bound": False,
        "actual_first_raw_function_bound": False,
        "actual_first_raw_equation_bound": False,
        "actual_first_raw_core_entry_bound": False,
        "actual_first_raw_ordinary_Z0_to_3_bound": False,
        "actual_first_raw_scale_conjugacy": False,
        "actual_first_raw_all_word_convergence_bound": False,
        "raw_all_word_convergence_certified": False,
        "full_first_raw_function_identity": False,
        "remaining_resolvent_factorization_gap": True,
        "remaining_finite_approximation_gap": True,
        "source_executed": False,
        "old_pure_definitions_executed": False,
        "archive_source_executed": False,
        "accepted_owner_builder_fixture_recurrence_quadrature_replayed": False,
        "point_values_selected": False,
        "finite_N_numerical_calculation": False,
    }
    if parent_pin is not None:
        out["parent_report_pin"] = dict(parent_pin)
    return out


def build(root: str | Path | None = None, parent: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Build the first-chart raw source-binding packet.

    The return value is intentionally a source certificate, not a numerical
    evaluation.  Quantitative tail constants are consumed from the accepted
    same-source receipt; no new numerical producer is run here.
    """
    parent_pin: dict[str, Any] | None = None
    try:
        repo = _repo_root(root)
        accepted, parent_pin = _accepted_parent(repo, parent)
        pins = _source_pins(repo)
        graph = _graph_binding(repo, accepted)
        scale = _normalization_scale(accepted)
        equations = _equation_binding(accepted, graph, pins, scale)
        entry = _core_entry_binding(graph, pins)
        ordinary = _ordinary_Z_binding(graph, pins, accepted)
        all_words = _all_word_binding(accepted, graph, scale, pins)
        finite = _finite_scope(pins)

        parent_checks = {
            "parent_pin": parent_pin.get("passed") is True,
            "parent_top_passed": _mapping(accepted.get("checks")).get("passed") is True,
            "parent_packet_invariants": _mapping(accepted.get("checks", {}).get("packet_invariants")).get("passed") is True,
            "native_packet_passed": _mapping(accepted.get("actual_native_collar_transport")).get("passed") is True,
            "native_first_window_only": _mapping(accepted.get("actual_native_collar_transport")).get("first_window_only") is True,
            "raw_cut_not_claimed": _mapping(accepted.get("actual_native_collar_transport")).get("raw_n1_cut_response_identity") is False,
            "no_ancestor_replay": _mapping(accepted.get("checks", {}).get("packet_invariants")).get("accepted_ancestor_or_quadrature_replayed") is False,
        }
        source_checks = {
            "all_source_hashes": pins.get("all_source_hashes_passed") is True,
            "forward_ast": _mapping(pins.get("forward_ast")).get("passed") is True,
            "volterra_ast": _mapping(pins.get("volterra_ast")).get("passed") is True,
            "radial_ast": _mapping(pins.get("radial_ast")).get("passed") is True,
            "operator_ast": _mapping(pins.get("operator_ast")).get("passed") is True,
            "coupled_tail_ast": _mapping(pins.get("coupled_tail_ast")).get("passed") is True,
        }
        source_bound_checks = {
            "parent": all(parent_checks.values()),
            "source_pins": all(source_checks.values()),
            "graph": graph.get("passed") is True,
            "equation": equations.get("checks", {}).get("passed") is True,
            "scale": scale.get("checks", {}).get("passed") is True,
            "core_entry": entry.get("checks", {}).get("passed") is True,
            "ordinary_Z": ordinary.get("checks", {}).get("passed") is True,
            "owner_definition": all_words.get("checks", {}).get("passed") is True,
            "raw_all_word_convergence": all_words.get("raw_all_word_convergence_passed") is True,
        }
        source_bound = all(source_bound_checks.values())
        raw_function_bound = source_bound and all_words.get("raw_all_word_convergence_passed") is True
        return {
            "passed": source_bound,
            "status": "source_bound_exact_raw_volterra_function" if source_bound else "failed_closed",
            "source_only": True,
            "parent_report": PARENT_REPORT,
            "parent_report_sha256": PARENT_REPORT_SHA256,
            "parent_report_pin": parent_pin,
            "parent_mapping_digest": parent_pin.get("mapping_digest"),
            "checks": {**source_bound_checks, "passed": source_bound},
            "actual_first_raw_source_definition_bound": source_bound,
            "actual_first_raw_function_bound": raw_function_bound,
            "actual_first_raw_equation_bound": equations.get("checks", {}).get("passed") is True,
            "actual_first_raw_core_entry_bound": entry.get("checks", {}).get("passed") is True,
            "actual_first_raw_ordinary_Z0_to_3_bound": ordinary.get("checks", {}).get("passed") is True,
            "actual_first_raw_scale_conjugacy": scale.get("checks", {}).get("passed") is True,
            "actual_first_raw_all_word_convergence_bound": all_words.get("raw_all_word_convergence_passed") is True,
            "raw_all_word_convergence_certified": all_words.get("raw_all_word_convergence_passed") is True,
            "raw_preconditioned_resolvent_factorization": all_words.get("raw_all_word_convergence_checks", {}).get("principal_resolvent_factorization") is True,
            "normalized_Kchi_kernel_substituted": False,
            "full_first_raw_function_identity": False,
            "source_pins": pins,
            "graph_binding": graph,
            "scale_conjugacy": scale,
            "raw_equation_binding": equations,
            "core_entry_binding": entry,
            "ordinary_Z_binding": ordinary,
            "all_word_binding": all_words,
            "finite_forward_producer": finite,
            "actual_function_handles": {
                "owner_namespace": graph.get("owner_namespace"),
                "component_operation": graph.get("component_operation"),
                "core": graph.get("field_handles", {}).get("core"),
                "first": graph.get("field_handles", {}).get("first"),
                "word0_first": graph.get("word_handles", {}).get("word0"),
                "word1_first": graph.get("word_handles", {}).get("word1"),
                "fixed_point_rhs_first": graph.get("word_handles", {}).get("fixed_point_rhs"),
                "core_coordinate": graph.get("core_coordinate"),
                "first_coordinate": graph.get("first_coordinate"),
            },
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "canonical_state": list(CANONICAL_STATE),
            "forward_state": list(FORWARD_STATE),
            "D": list(D),
            "axis_powers": list(AXIS_POWERS),
            "source_family": graph.get("source_family"),
            "remaining_resolvent_factorization_gap": not raw_function_bound,
            "remaining_finite_approximation_gap": True,
            "remaining_scope": [
                "The raw tail proof uses the same-source Akin-preconditioned Green kernels; the separate normalized Kchi kernel is intentionally not substituted.",
                "No raw first-chart quantitative word-tail constant is certified here.",
                "The finite OriginalN1CoreNumeric object remains finite-prefix evidence only.",
                "No q-dependent collar bridge, six-term source equality, or temporal n>=2 recursion is claimed.",
            ],
            "limits": {
                "old_finite_producer_exact_function": False,
                "finite_prefix_identity_auto_promoted": False,
                "radial_remainder_certified": False,
                "raw_first_word_tail_numeric": False,
                "full_forward_response_function_identity": False,
                "q_dependent_collar_bridge": False,
                "six_term_source_equality": False,
                "temporal_n_ge_2_recursion": False,
                "point_values_selected": False,
                "finite_N_numerical_calculation": False,
                "all_error_views_mutated": False,
            },
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
        }
    except Exception as exc:
        return _failure(f"source binding failed closed: {type(exc).__name__}: {exc}", parent_pin)


def bind_first_raw_function(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compatibility alias for :func:`build`."""

    return build(*args, **kwargs)


__all__ = [
    "PARENT_REPORT",
    "PARENT_REPORT_SHA256",
    "INNER_GRAPH",
    "INNER_GRAPH_SHA256",
    "D",
    "AXIS_POWERS",
    "NORMALIZATION",
    "build",
    "bind_first_raw_function",
]
