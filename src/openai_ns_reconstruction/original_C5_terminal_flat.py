"""Source bound terminal flat-loop certificate for the actual O3 transition.

The O3 transition endpoint is a source identity, not a numerical phase
sample.  At the endpoint ``x=1`` the actual background source has

    E_y = -E (1/2 + mu sigma(x)),   sigma(1) = 1,
    a = (E - 2 E_y)/E,             Delta = a - 2.

Consequently ``Delta(1,Z) = 2 mu`` for every ordinary Z jet.  The pinned O3
admission receipt proves ``2 mu >= eta`` from the original logarithmic
source inequality.  The original lazy loop then returns zero, including its
ordinary Z derivatives and the native fixed-phase slow-y derivative.  This
module only reads the saved graph, source ASTs, and the accepted admission
receipt.  It never imports or executes an archived source owner and it does
not evaluate C, H, a phase, or a quadrature.
"""
from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


PARENT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_C5_seam_joined_history.json"
)
PARENT_SHA256 = (
    "35fdc828c112e5fb09bf28a3f709b1d8a9817b54cd3a99fc11f76cfd92251d16"
)
LEGACY_SEAM_PARENT_SHA256 = (
    "67945b4e9a8f2b1cdaa9a255b42a3bd9e5377f514cb12409ba41e9f17f957af0"
)

OUTER_GRAPH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz"
)
OUTER_GRAPH_SHA256 = (
    "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d"
)
OUTER_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py"
)
OUTER_SOURCE_SHA256 = (
    "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3"
)
LOOP_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
)
LOOP_SOURCE_SHA256 = (
    "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"
)
Q4_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
)
Q4_SOURCE_SHA256 = (
    "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7"
)
ADMISSION_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N.py"
)
ADMISSION_SOURCE_SHA256 = (
    "83e5a8bdae75d20a261d24dd22d30cc7cb5e603e8fa2f0133958ca880b5e76db"
)
ADMISSION_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N_check.json"
)
ADMISSION_RECEIPT_SHA256 = (
    "30bd8cd0d502a4e2e9a8a45b9ad37fa55054ea3dd3e1dd690c5fc3c7700a3270"
)
ADMISSION_AST_SHA256 = (
    "8507a12bb9de7e54ee103da586c0d50fa3e8fafa3e0a5d5034b0dde9e8dbb3f4"
)

ORDINARY_Z_ORDERS = (0, 1, 2, 3)
O3_COORDINATE = 241536
O3_ENDPOINT_VALUE = 1
O3_PHASE = 242551
O3_FLAT_PREDICATE = 242552
O3_DELTA = 242546
O3_ETA = 21
O3_MU = 17
O3_LOG_MU = 16
O3_SIGMA = 241559
O3_E = (241552, 241553, 241554, 241555)
O3_EY_MOVING = (241563, 241564, 241565, 241566)
O3_EY_ENDPOINT = (242397, 242398, 242399, 242400)
O3_VY_ENDPOINT = 242346
O3_A = (242476, 242480, 242486, 242493)
O3_DELTA_ROWS = (242546, 242480, 242486, 242493)
O3_B = (244136, 244137, 244138, 244139)
O3_A_ROWS = (244098, 244099, 244100, 244101)
O3_BPHI = (244319, 244320, 244321, 244322)
O3_APHI = (244217, 244218, 244219, 244220)
O3_SLOW_B_Y = (248598, 248632, 248685, 248759)
O3_SLOW_JACOBIAN = 172066


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _load(path: Path) -> Any:
    raw = path.read_bytes()
    if path.name.endswith(".gz"):
        raw = gzip.decompress(raw)
    return json.loads(raw)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _repo(root: Any) -> Path:
    path = Path(root)
    if path.is_file():
        path = path.parent
    candidates = (path, path.parent, Path(__file__).resolve().parents[2])
    for candidate in candidates:
        if (candidate / "src").is_dir() and (candidate / "experiments").is_dir():
            return candidate.resolve()
    raise ValueError("repository root with src/ and experiments/ is required")


def _node(nodes: Sequence[Any], node_id: Any) -> Mapping[str, Any]:
    if isinstance(node_id, int) and 0 <= node_id < len(nodes):
        value = nodes[node_id]
        return value if isinstance(value, Mapping) else {}
    return {}


def _exact(node: Mapping[str, Any], numerator: int, denominator: int = 1) -> bool:
    return (
        node.get("operation") == "exact_rational"
        and node.get("numerator") == numerator
        and node.get("denominator") == denominator
    )


def _contains_text(node: Mapping[str, Any], *parts: str) -> bool:
    text = str(node.get("definition", ""))
    return all(part in text for part in parts)


def _pin(repo: Path, relative: str, expected: str) -> dict[str, Any]:
    path = repo / relative
    try:
        actual = _sha(path)
    except (FileNotFoundError, OSError):
        actual = None
    return {
        "path": relative,
        "sha256": actual,
        "expected_sha256": expected,
        "passed": actual == expected,
        "source_executed": False,
    }


def _ast_guard(
    repo: Path,
    relative: str,
    symbol: str,
    required: Sequence[str],
    expected_ast: str | None = None,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": relative,
        "function": symbol,
        "required_text": {item: False for item in required},
        "AST_sha256": None,
        "passed": False,
        "source_executed": False,
    }
    try:
        source = (repo / relative).read_text(encoding="utf-8")
        tree = ast.parse(source)
    except (FileNotFoundError, OSError, SyntaxError) as exc:
        result["error"] = type(exc).__name__ + ": " + str(exc)
        return result
    matches: list[ast.AST] = []
    if "." in symbol:
        class_name, method_name = symbol.split(".", 1)
        for item in ast.walk(tree):
            if isinstance(item, ast.ClassDef) and item.name == class_name:
                matches.extend(
                    child
                    for child in item.body
                    if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef))
                    and child.name == method_name
                )
    else:
        matches = [
            item
            for item in ast.walk(tree)
            if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef))
            and item.name == symbol
        ]
    if len(matches) != 1:
        result["error"] = "expected one AST function, found " + str(len(matches))
        return result
    function = matches[0]
    segment = ast.get_source_segment(source, function) or ""
    compact = "".join(segment.split())
    checks = {item: "".join(item.split()) in compact for item in required}
    ast_sha = hashlib.sha256(
        ast.dump(function, include_attributes=False).encode()
    ).hexdigest()
    parsed_source_sha = hashlib.sha256(
        ast.dump(ast.parse(segment)).encode()
    ).hexdigest()
    result.update(
        {
            "required_text": checks,
            "AST_sha256": ast_sha,
            "parsed_source_AST_sha256": parsed_source_sha,
            "expected_AST_sha256": expected_ast,
            # The historical source binding records ast.parse(inspect.getsource(fn)),
            # whereas the local guards normally record the function node.  Accept
            # the former only when it is the pinned expected binding hash.
            "AST_hash_matches": expected_ast is None
            or ast_sha == expected_ast
            or parsed_source_sha == expected_ast,
        }
    )
    result["passed"] = all(checks.values()) and result["AST_hash_matches"]
    return result


def _parent_pin(repo: Path, parent: Any) -> tuple[Mapping[str, Any], dict[str, Any]]:
    path = repo / PARENT
    saved = _load(path)
    selected = saved if parent is None else parent
    actual_sha = _sha(path)
    mapping_match = _digest(saved) == _digest(selected)
    pin = {
        "path": PARENT,
        "sha256": actual_sha,
        "expected_sha256": PARENT_SHA256,
        "passed": actual_sha == PARENT_SHA256,
        "mapping_match": mapping_match,
        "source_executed": False,
    }
    if not mapping_match:
        raise ValueError("supplied C5 joined-history parent is not the accepted mapping")
    return selected, pin


def _flat_row_checks(
    nodes: Sequence[Any],
    row_ids: Sequence[int],
    *,
    slow: bool = False,
) -> dict[str, Any]:
    rows = [_node(nodes, row_id) for row_id in row_ids]
    expected = list(ORDINARY_Z_ORDERS)
    orders = [row.get("ordinary_slow_Z_derivative_order") for row in rows]
    checks = {
        "row_count": len(rows) == len(expected),
        "ordinary_Z_orders": orders == expected,
        "all_original_lazy_flat_branch": all(
            row.get("operation") == "original_lazy_flat_branch" for row in rows
        ),
        "shared_actual_Delta": all(row.get("Delta") == O3_DELTA for row in rows),
        "shared_actual_eta": all(row.get("eta") == O3_ETA for row in rows),
        "flat_value_exact_zero": all(row.get("flat_value") == 0 for row in rows),
        "shared_flat_predicate": all(
            row.get("flat_predicate") == O3_FLAT_PREDICATE for row in rows
        ),
        "active_body_not_evaluated": all(
            row.get("active_body_not_evaluated_on_flat_branch") is True
            for row in rows
        ),
        "same_flat_extension_through_order5": all(
            row.get("same_original_flat_extension_through_order5") is True
            for row in rows
        ),
    }
    if slow:
        checks["same_original_lazy_native_derivative"] = all(
            row.get("same_original_lazy_native_derivative") is True for row in rows
        )
    checks["passed"] = all(checks.values())
    return checks


def _slow_row_checks(nodes: Sequence[Any]) -> dict[str, Any]:
    rows = [_node(nodes, row_id) for row_id in O3_SLOW_B_Y]
    body_ids = []
    for row in rows:
        arguments = row.get("arguments", ())
        body_ids.append(next((item for item in arguments if item != O3_SLOW_JACOBIAN), None))
    wrapper_checks = {
        "ordinary_product_rows": all(row.get("operation") == "product" for row in rows),
        "shared_unit_physical_jacobian": all(
            O3_SLOW_JACOBIAN in row.get("arguments", ()) for row in rows
        ),
        "jacobian_node_is_one": _node(nodes, O3_SLOW_JACOBIAN).get("operation")
        == "positive_quotient"
        and _node(nodes, O3_SLOW_JACOBIAN).get("numerator") == 1
        and _node(nodes, O3_SLOW_JACOBIAN).get("denominator") == 1,
        "lazy_native_bodies_present": all(isinstance(item, int) for item in body_ids),
    }
    body = _flat_row_checks(nodes, body_ids, slow=True) if all(
        isinstance(item, int) for item in body_ids
    ) else {"passed": False}
    return {"wrapper": wrapper_checks, "body": body, "passed": all(wrapper_checks.values()) and body.get("passed") is True}


def _source_and_admission_checks(repo: Path, graph: Mapping[str, Any]) -> dict[str, Any]:
    nodes = graph["exact_graph_nodes"]
    bindings = _mapping(graph.get("true_outer_C6_C5_bindings"))
    mu_proof = _mapping(bindings.get("mu_proof"))
    admission = _load(repo / ADMISSION_RECEIPT)
    graph_inputs = _mapping(graph.get("input_hashes"))
    checks: dict[str, Any] = {
        "graph_admission_source_hash_pin": graph_inputs.get(Path(ADMISSION_SOURCE).name)
        == ADMISSION_SOURCE_SHA256,
        "graph_admission_receipt_hash_pin": graph_inputs.get(Path(ADMISSION_RECEIPT).name)
        == ADMISSION_RECEIPT_SHA256,
        "admission_receipt_sha": _sha(repo / ADMISSION_RECEIPT) == ADMISSION_RECEIPT_SHA256,
        "admission_all_passed": admission.get("all_passed") is True,
        "admission_source_family": admission.get("source_family") == graph.get("source_family"),
        "admission_transition_identity": _mapping(
            admission.get("independent_O3_source_identities")
        ).get("exact_transition_a_minus2_equals_2mu_sigma_and_power_equals2mu")
        is True,
        "admission_nonnegative_Delta": _mapping(
            admission.get("independent_O3_source_identities")
        ).get("nonnegative_Delta_active_q_squared_le_eta_over_a_and_a_ge2_implies_eta_over2")
        is True,
        "admission_flat_independent_p2": _mapping(
            admission.get("independent_O3_source_identities")
        ).get("original_flat_q_when_Delta_ge_eta_independent_of_p2")
        is True,
        "admission_positive_mu_transition": admission.get(
            "positive_mu_source_Delta_active_transition_and_exact_quiet_power_retained"
        )
        is True,
        "mu_proof_source_binding": _mapping(mu_proof.get("original_admission_binding")).get(
            "module"
        )
        == Path(ADMISSION_SOURCE).name,
        "mu_proof_AST_pin": _mapping(mu_proof.get("original_admission_binding")).get(
            "AST_sha256"
        )
        == ADMISSION_AST_SHA256,
    }
    checks["ast_admission"] = _ast_guard(
        repo,
        ADMISSION_SOURCE,
        "source_mu_admission",
        (
            "logmu=c.ln(c.mpf('.001'))-4*(c.exp(40)+11)",
            "margin=logmu+c.ln(2)-eta_log",
            "ifep(margin)[0]<0:raiseValueError('Same-source O3 quiet power 2mu>=eta not proved')",
        ),
        ADMISSION_AST_SHA256,
    )
    checks["ast_admission_passed"] = checks["ast_admission"].get("passed") is True
    return checks


def build(root: Any, parent: Any = None) -> dict[str, Any]:
    """Build a fail-closed source receipt for the O3 outgoing flat endpoint."""
    repo = _repo(root)
    accepted, parent_pin = _parent_pin(repo, parent)
    graph = _load(repo / OUTER_GRAPH)
    nodes = graph.get("exact_graph_nodes", ())
    bindings = _mapping(graph.get("true_outer_C6_C5_bindings"))
    charts = _mapping(bindings.get("charts"))
    chart = _mapping(charts.get("O3_transition"))
    roots = _mapping(chart.get("roots"))
    loop = _mapping(chart.get("loop"))
    radial = _mapping(chart.get("radial"))
    endpoint = _mapping(_mapping(bindings.get("background_endpoints")).get("O3_transition"))

    source_pins = {
        "outer_graph": {
            "path": OUTER_GRAPH,
            "sha256": _sha(repo / OUTER_GRAPH),
            "expected_sha256": OUTER_GRAPH_SHA256,
            "passed": _sha(repo / OUTER_GRAPH) == OUTER_GRAPH_SHA256,
            "source_executed": False,
        },
        "outer_source": _pin(repo, OUTER_SOURCE, OUTER_SOURCE_SHA256),
        "prescribed_loop_source": _pin(repo, LOOP_SOURCE, LOOP_SOURCE_SHA256),
        "Q4_source": _pin(repo, Q4_SOURCE, Q4_SOURCE_SHA256),
        "admission_source": _pin(repo, ADMISSION_SOURCE, ADMISSION_SOURCE_SHA256),
        "admission_receipt": _pin(repo, ADMISSION_RECEIPT, ADMISSION_RECEIPT_SHA256),
    }

    checks: dict[str, Any] = {
        "accepted_parent_sha": parent_pin.get("passed") is True,
        "accepted_parent_mapping": parent_pin.get("mapping_match") is True,
        "accepted_parent_passed": _mapping(accepted.get("checks")).get("passed") is True,
        "accepted_seam_passed": _mapping(accepted.get("actual_C5_O2_O3_seam")).get(
            "passed"
        )
        is True,
        "accepted_joined_history_passed": _mapping(
            accepted.get("joined_C5_signed_history")
        ).get("passed")
        is True,
        "outer_graph_hash": source_pins["outer_graph"]["passed"],
        "outer_source_hash": source_pins["outer_source"]["passed"],
        "loop_source_hash": source_pins["prescribed_loop_source"]["passed"],
        "Q4_source_hash": source_pins["Q4_source"]["passed"],
        "admission_source_hash": source_pins["admission_source"]["passed"],
        "admission_receipt_hash": source_pins["admission_receipt"]["passed"],
        "source_family_matches_parent": graph.get("source_family")
        == accepted.get("source_family"),
        "actual_chart_name": chart.get("coordinate") == O3_COORDINATE,
        "actual_chart_domain": chart.get("domain") == [0, O3_ENDPOINT_VALUE],
        "actual_chart_phase": chart.get("actual_phase") == O3_PHASE,
        "actual_flat_predicate_anchor": loop.get("flat_predicate") == O3_FLAT_PREDICATE,
        "actual_delta_anchor": roots.get("Delta", [None])[0] == O3_DELTA,
        "actual_eta_anchor": bindings.get("eta") == O3_ETA,
        "actual_mu_anchor": bindings.get("mu") == O3_MU,
        "actual_log_mu_anchor": bindings.get("logmu") == O3_LOG_MU,
        "actual_A_rows": tuple(loop.get("A", ())[:4]) == O3_A_ROWS,
        "actual_B_rows": tuple(loop.get("B_over_Pstar", ())[:4]) == O3_B,
        "actual_A_phi_rows": tuple(loop.get("A_phi", ())[:4]) == O3_APHI,
        "actual_B_phi_rows": tuple(loop.get("B_over_Pstar_phi", ())[:4]) == O3_BPHI,
        "actual_slow_B_y_rows": tuple(radial.get("slow_B_y", ())[:4]) == O3_SLOW_B_Y,
        "actual_root_a_rows": tuple(roots.get("a", ())[:4]) == O3_A,
        "actual_root_Delta_rows": tuple(roots.get("Delta", ())[:4]) == O3_DELTA_ROWS,
    }
    checks["source_pins"] = all(pin.get("passed") is True for pin in source_pins.values())

    checks["AST_outer_O3"] = _ast_guard(
        repo,
        OUTER_SOURCE,
        "CurrentOriginalTrueOuterC6C5Functions.__init__",
        (
            "elifname=='O3_transition':",
            "Ey=a.scale(E,g.neg(g.add(g.constant('1/2'),g.mul(self.mu,sigma(g,x)))))",
            "V=a.fixed(0);Ey=a.scale(E,g.neg(g.add(g.constant('1/2'),g.mul(self.mu,sigma(g,x)))))",
        ),
        "6fa29f53fbecfd2d9dc79d0184032bb434a32a228ec6a318eb06db237c52c593",
    )
    checks["AST_generic_roots"] = _ast_guard(
        repo,
        LOOP_SOURCE,
        "generic_roots",
        (
            "field=chart['fields'];radial=chart['radial'];generic=chart['generic']",
            "C=a.sub(E,a.scale(a.trim(radial['E_y']),2))",
            "B=a.scale(a.trim(radial['V_y']),2)",
            "aa,bb=divide(C),divide(B)",
            "returndict(E=E,V=V,a=aa,b=bb,p1=p1,p2=p2,t0=t0,Delta=a.sub(kap,a.fixed(2)),kappa=kap)",
        ),
        "24570272611c6745b9e785163e80c67374e6ca1044a90ebebed9ff7dfdd2c708",
    )
    checks["AST_loop_flat"] = _ast_guard(
        repo,
        LOOP_SOURCE,
        "loop_functions",
        (
            "shear=g.node('exact_real_comparison',operator='>=',left=roots['Delta'].value.node,right=eta.node)",
            "lazy=lambdajet:a.lazy(jet,predicate,roots['Delta'].value,eta,coordinate,collar)",
            "phase_endpoints=[0,1]",
            "flat_inverse='2*pi*phase'",
            "branch_inverse_is_original_not_support_cap=True",
        ),
        "fffc175d358e7cce71d3046ec10303e22d46b5ebb3145f9a44c0d1fe88796f8c",
    )
    checks["AST_slow_y"] = _ast_guard(
        repo,
        Q4_SOURCE,
        "total_radial_functions",
        (
            "invjac=g.quotient(g.one,jac,'positive original chart log-radius Jacobian')",
            "slow={key:FunctionJet(tuple(g.mul(invjac,derivative.derivative(ref))forrefinrows(chart['loop'][key])))forkeyin('A','B_over_Pstar')}",
            "slow_B_y=slow['B_over_Pstar']",
            "physical_y_equals_native_partial_divided_by_Jacobian=True",
        ),
        "8c71d236f7e0d17ba836ce2d87070dd8616ecea154dec27cd4926ae39202bed5",
    )
    checks["AST"] = all(
        checks[key].get("passed") is True
        for key in ("AST_outer_O3", "AST_generic_roots", "AST_loop_flat", "AST_slow_y")
    )

    # The graph-level endpoint algebra is deliberately checked from actual
    # nodes.  This is the proof of Delta(1,Z)=2*mu, rather than a copied
    # descriptive flag.
    sigma_node = _node(nodes, O3_SIGMA)
    ey0 = _node(nodes, O3_EY_MOVING[0])
    ey_factor = _node(nodes, 241562)
    ey_sum = _node(nodes, 241561)
    ey_mu_product = _node(nodes, 241560)
    a0 = _node(nodes, O3_A[0])
    a_num = _node(nodes, 242470)
    a_neg_twice_ey = _node(nodes, 242464)
    a_twice_ey = _node(nodes, 242458)
    delta0 = _node(nodes, O3_DELTA)
    delta_minus_two = _node(nodes, 25304)
    endpoint_ey0 = _node(nodes, O3_EY_ENDPOINT[0])
    endpoint_sigma = {
        "coordinate_node": O3_COORDINATE,
        "endpoint_value_node": O3_ENDPOINT_VALUE,
        "coordinate_is_bound_variable": _node(nodes, O3_COORDINATE).get("operation")
        == "bound_variable",
        "endpoint_value_is_exact_one": _exact(_node(nodes, O3_ENDPOINT_VALUE), 1),
        "sigma_node": O3_SIGMA,
        "sigma_operation": sigma_node.get("operation") == "exact_original_source_flat_sigma",
        "sigma_argument_coordinate": sigma_node.get("argument") == O3_COORDINATE,
        "sigma_source_has_one_endpoint": _contains_text(
            sigma_node, "0 for x<=0", "1 for x>=1"
        ),
        "endpoint_Ey_is_actual_function_substitution": endpoint_ey0.get(
            "operation"
        )
        == "function_substitution"
        and endpoint_ey0.get("expression") == O3_EY_MOVING[0]
        and endpoint_ey0.get("variable") == O3_COORDINATE
        and endpoint_ey0.get("value") == O3_ENDPOINT_VALUE,
    }
    checks["sigma_endpoint"] = all(endpoint_sigma.values())

    endpoint_algebra = {
        "mu_graph_is_exp_logmu": _node(nodes, O3_MU).get("operation")
        == "analytic_unary"
        and _node(nodes, O3_MU).get("name") == "exp"
        and _node(nodes, O3_MU).get("argument") == O3_LOG_MU,
        "logmu_exact_recipe": _node(nodes, O3_LOG_MU).get("operation") == "sum"
        and set(_node(nodes, O3_LOG_MU).get("arguments", ())) == {12, 15},
        "Ey_product_is_E_times_factor": ey0.get("operation") == "product"
        and ey0.get("arguments") == [O3_E[0], 241562],
        "Ey_factor_is_negative_half_plus_mu_sigma": ey_factor.get("operation")
        == "negative"
        and ey_factor.get("argument") == 241561
        and ey_sum.get("operation") == "sum"
        and set(ey_sum.get("arguments", ())) == {28, 241560}
        and ey_mu_product.get("operation") == "product"
        and set(ey_mu_product.get("arguments", ())) == {O3_MU, O3_SIGMA},
        "a_is_E_minus_2Ey_over_E": a0.get("operation") == "positive_quotient"
        and a0.get("numerator") == 242470
        and a0.get("denominator") == O3_E[0]
        and a_num.get("operation") == "sum"
        and set(a_num.get("arguments", ())) == {O3_E[0], 242464}
        and a_neg_twice_ey.get("operation") == "negative"
        and a_neg_twice_ey.get("argument") == 242458
        and a_twice_ey.get("operation") == "product"
        and set(a_twice_ey.get("arguments", ())) == {O3_EY_MOVING[0], 48},
        "Delta_is_a_minus_two": delta0.get("operation") == "sum"
        and set(delta0.get("arguments", ())) == {O3_A[0], 25304}
        and delta_minus_two.get("operation") == "negative"
        and delta_minus_two.get("argument") == 48
        and _exact(_node(nodes, 48), 2),
        "Delta_Z_rows_are_actual_a_rows": tuple(roots.get("Delta", ())[:4])
        == O3_DELTA_ROWS
        and tuple(roots.get("a", ())[:4]) == O3_A,
        "V_y_root_rows_exact_zero": all(_exact(_node(nodes, item), 0) for item in roots.get("b", ())[:4]),
        "endpoint_Ey_rows_present": tuple(endpoint.get("radial", {}).get("E_y", ())[:4])
        == O3_EY_ENDPOINT,
        "endpoint_Vy_rows_exact_zero": tuple(endpoint.get("radial", {}).get("V_y", ())[:4])
        == (O3_VY_ENDPOINT,) * 4
        and _node(nodes, O3_VY_ENDPOINT).get("operation") == "function_substitution"
        and _node(nodes, O3_VY_ENDPOINT).get("expression") == 0
        and _node(nodes, O3_VY_ENDPOINT).get("value") == O3_ENDPOINT_VALUE,
    }
    checks["endpoint_algebra"] = all(endpoint_algebra.values())

    admission_checks = _source_and_admission_checks(repo, graph)
    checks["mu_admission"] = all(
        value is True
        for key, value in admission_checks.items()
        if key != "ast_admission"
    ) and admission_checks.get("ast_admission_passed") is True

    flat_rows = {
        "A": _flat_row_checks(nodes, O3_A_ROWS),
        "B_over_Pstar": _flat_row_checks(nodes, O3_B),
        "A_phi": _flat_row_checks(nodes, O3_APHI),
        "B_over_Pstar_phi": _flat_row_checks(nodes, O3_BPHI),
    }
    slow_rows = _slow_row_checks(nodes)
    checks["flat_support_all_ordinary_jets"] = all(item.get("passed") is True for item in flat_rows.values())
    checks["flat_support_fixed_phase_slow_B_y"] = slow_rows.get("passed") is True
    checks["endpoint_flat_predicate"] = (
        _node(nodes, O3_FLAT_PREDICATE).get("operation") == "exact_real_comparison"
        and _node(nodes, O3_FLAT_PREDICATE).get("operator") == ">="
        and _node(nodes, O3_FLAT_PREDICATE).get("left") == O3_DELTA
        and _node(nodes, O3_FLAT_PREDICATE).get("right") == O3_ETA
    )

    outgoing_zero = all(
        (
            checks["accepted_parent_sha"],
            checks["accepted_parent_mapping"],
            checks["accepted_parent_passed"],
            checks["accepted_seam_passed"],
            checks["accepted_joined_history_passed"],
            checks["source_pins"],
            checks["source_family_matches_parent"],
            checks["actual_chart_name"],
            checks["actual_chart_domain"],
            checks["actual_chart_phase"],
            checks["actual_flat_predicate_anchor"],
            checks["actual_delta_anchor"],
            checks["actual_A_rows"],
            checks["actual_B_rows"],
            checks["actual_A_phi_rows"],
            checks["actual_B_phi_rows"],
            checks["actual_slow_B_y_rows"],
            checks["actual_root_a_rows"],
            checks["actual_root_Delta_rows"],
            checks["AST"],
            checks["sigma_endpoint"],
            checks["endpoint_algebra"],
            checks["mu_admission"],
            checks["flat_support_all_ordinary_jets"],
            checks["flat_support_fixed_phase_slow_B_y"],
            checks["endpoint_flat_predicate"],
        )
    )
    return {
        "passed": outgoing_zero,
        "source_only": True,
        "source_executed": False,
        "accepted_parent": PARENT,
        "accepted_parent_sha256": PARENT_SHA256,
        "legacy_seam_parent_sha256": LEGACY_SEAM_PARENT_SHA256,
        "source_pins": source_pins,
        "AST": {
            key: checks[key]
            for key in ("AST_outer_O3", "AST_generic_roots", "AST_loop_flat", "AST_slow_y")
        },
        "admission": admission_checks,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "endpoint": {
            "chart": "O3_transition",
            "coordinate_node": O3_COORDINATE,
            "coordinate_value": O3_ENDPOINT_VALUE,
            "phase_node": O3_PHASE,
            "delta_node": O3_DELTA,
            "delta_endpoint_identity": "Delta(1,Z)=2*mu for ordinary Z orders 0..3",
            "sigma_endpoint": endpoint_sigma,
            "algebra": endpoint_algebra,
        },
        "flat_support": flat_rows,
        "fixed_phase_slow_B_y": slow_rows,
        "checks": checks,
        "outgoing_C_Z0_Z3_exact_zero": outgoing_zero,
        "outgoing_C_y_slow_Z0_Z3_exact_zero": outgoing_zero,
        "arbitrary_issued_phase": True,
        "phase_not_materialized": True,
        "actual_N_not_materialized": _node(nodes, 18).get("outer_integer_and_inner_power_not_materialized") is True,
        "numeric_C_values_evaluated": False,
        "numeric_H_values_evaluated": False,
        "slow_integral_values_evaluated": False,
        "cap_contraction_or_improvement": False,
        "only_terminal_C_endpoint_contracted": True,
        "incoming_history_or_slow_integrals_zeroed": False,
        "limitations": [
            "This proves only the outgoing O3 C and fixed-phase C_y endpoint identities.",
            "It does not evaluate the incoming H0, Cin, I1, I2, or the actual huge-N phase.",
            "It does not claim numeric signed C/H values, a cap improvement, full n1, or terminal closure.",
        ],
    }


def build_original_C5_terminal_flat(root: Any, parent: Any = None) -> dict[str, Any]:
    return build(root, parent)


__all__ = [
    "PARENT",
    "PARENT_SHA256",
    "ORDINARY_Z_ORDERS",
    "build",
    "build_original_C5_terminal_flat",
]
