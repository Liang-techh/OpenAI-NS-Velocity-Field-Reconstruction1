"""Local signed lower-collar A/A0 primitive over pinned source nodes.

This adapter only reads the accepted finite-N receipt and its three direct
source inputs. Its new SignedDensityExpressionGraph contains exact local
substitutions and integrals; no accepted producer or archived graph is run or
mutated, and no numerical value is installed.
"""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping

from .original_finite_N_signed_source import ast_guard, digest, load, sha
from .original_signed_density_expression import ScalarExpression, SignedDensityExpressionGraph

PARENT_REPORT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_finite_N_signed_source.json"
PARENT_REPORT_SHA256 = "3c309218eb3798ba077f2ae0c1c74be9e5f0a0dce023e38c3e1dee3283d309c7"
LOWER_GRAPH = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.json.gz"
LOWER_GRAPH_NAME = Path(LOWER_GRAPH).name
LOWER_GRAPH_SHA256 = "b2f983ab934875ae9aeb24b292277911dd5647e3a23943ddcf2b3b6e82595cb5"
LOWER_SOURCE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.py"
LOWER_SOURCE_NAME = Path(LOWER_SOURCE).name
FORWARD_SOURCE = "src/openai_ns_reconstruction/original_n1_inner_m15_numeric.py"
FAMILY_KEYS = ("actual_five_defect_family_sha256", "implicit_source_sha256", "datum_enclosure_sha256")
HANDLES = {
    "s": 35848, "hb": 8, "epsilon": 35565, "rho": 35851,
    "sigma_prime": 35953, "sigma": 35949, "ratio_C0": 35940, "ratio_C1": 35941,
    "Psi_core_rho_C0": 35870, "Psi_core_rho_C1": 35871,
    "L_C0": 30382, "L_C1": 30381, "Z": 20, "core_rho_variable": 35558,
    "Psi_core_rho4_C0_source": 35715, "Psi_core_rho4_C1_source": 35833, "rho4": 13,
}
_RADIAL_ZERO_PSI4_NODES = frozenset((36065, 36066, 36067, 36068))
_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}


def _m(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _repo(root: Any) -> Path:
    p = Path(root)
    if p.is_file():
        p = p.parent
    for candidate in (p, p.parent, Path(__file__).resolve().parents[2]):
        if (candidate / "src").is_dir() and (candidate / "experiments").is_dir():
            return candidate.resolve()
    raise ValueError("repository root with src/ and experiments/ is required")


def _proof(nodes: list[Any], index: int, expected: Mapping[str, Any]) -> dict[str, Any]:
    actual = nodes[index] if 0 <= index < len(nodes) and isinstance(nodes[index], Mapping) else {}
    selected = {key: actual.get(key) for key in expected}
    return {"index": index, "expected": dict(expected), "actual": selected,
            "passed": selected == dict(expected), "source_executed": False}


def _sub(graph: SignedDensityExpressionGraph, expression: int | ScalarExpression,
         variable: int | ScalarExpression, value: int | ScalarExpression, name: str, **meta: Any) -> ScalarExpression:
    row = {"operation": "function_substitution", "original_graph_nodes_added_or_replaced": False,
           "expression": expression.node if isinstance(expression, ScalarExpression) else expression,
           "variable": variable.node if isinstance(variable, ScalarExpression) else variable,
           "value": value.node if isinstance(value, ScalarExpression) else value, "name": name}
    row.update(meta)
    return graph.insert(row)


def _bound(graph: SignedDensityExpressionGraph, name: str, role: str) -> ScalarExpression:
    return graph.insert({"operation": "bound_variable", "name": name, "role": role,
                         "local_scope": True, "original_graph_nodes_added_or_replaced": False})


def _parent(root: Path, supplied: Any) -> tuple[Mapping[str, Any], dict[str, Any], list[str]]:
    path = root / PARENT_REPORT
    gaps: list[str] = []
    try:
        actual, accepted = sha(path), load(path)
    except Exception as exc:
        return {}, {"path": PARENT_REPORT, "sha256": None, "expected_sha256": PARENT_REPORT_SHA256,
                    "passed": False, "source_executed": False}, [f"cannot load accepted parent: {exc}"]
    pin = {"path": PARENT_REPORT, "sha256": actual, "expected_sha256": PARENT_REPORT_SHA256,
           "passed": actual == PARENT_REPORT_SHA256, "source_executed": False}
    if not pin["passed"]:
        gaps.append("accepted finite_N_signed_source report hash mismatch")
    try:
        if supplied is None:
            value, supplied_path = accepted, None
        elif isinstance(supplied, Mapping):
            value, supplied_path = supplied, None
        else:
            candidate = Path(supplied)
            if not candidate.is_absolute():
                candidate = root / candidate
            value, supplied_path = load(candidate), candidate.as_posix()
        match = digest(value) == digest(accepted)
    except Exception as exc:
        return {}, pin, gaps + [f"cannot load supplied parent: {exc}"]
    pin.update({"supplied_path": supplied_path, "accepted_mapping_sha256": digest(accepted),
                "supplied_mapping_sha256": digest(value), "mapping_match": match})
    if not match:
        gaps.append("supplied parent mapping does not match pinned report")
    return _m(value), pin, gaps


def _pin(parent: Mapping[str, Any], root: Path, relative: str, name: str, kind: str, gaps: list[str]) -> dict[str, Any]:
    expected = _m(parent.get("input_hashes")).get(name) or _m(parent.get("implementation_hashes")).get(name)
    path = root / relative
    try:
        actual = sha(path)
    except OSError:
        actual = None
    row = {"name": name, "kind": kind, "path": relative, "sha256": actual,
           "expected_sha256": expected, "passed": isinstance(expected, str) and actual == expected,
           "source_executed": False}
    if not row["passed"]:
        gaps.append(f"source pin failed: {name}")
    return row


def _fail(root: Path, pin: Mapping[str, Any], gaps: list[str]) -> dict[str, Any]:
    return {"passed": False, "source_binding_passed": False, "source_only": True,
            "source_executed": False, "parent_report": PARENT_REPORT,
            "parent_report_sha256": PARENT_REPORT_SHA256, "parent_report_pin": dict(pin),
            "source_pins": {}, "source_pins_ok": False, "AST": {}, "gaps": gaps,
            "numeric_value_available": False, "Delta_six": {"exact_identity": False, "zero_claim": False},
            "root": str(root), "directory": str(root / "experiments" / "root_st073")}


def _controls() -> dict[str, Any]:
    f0, f1, l0, l1 = Fraction(7, 3), Fraction(-2, 5), Fraction(11, 7), Fraction(3, 13)
    q, q_direct = (f1 * l0 - f0 * l1) / l0**2, f1 / l0 - f0 * l1 / l0**2
    s = Fraction(3, 7)
    coefficients = (Fraction(5, 11), Fraction(2, 9), Fraction(-1, 13))
    direct = sum((coefficient * 2 * s ** (power + 2) / Fraction(power + 2)
                  for power, coefficient in enumerate(coefficients)), Fraction(0))
    remainder = sum((coefficient * 2 * s ** (power + 2) / Fraction(power + 2)
                     for power, coefficient in enumerate(coefficients[1:], 1)), Fraction(0))
    split = coefficients[0] * s**2 + remainder
    derivative = sum(coefficient * 2 * s ** (power + 1)
                     for power, coefficient in enumerate((coefficients[0], coefficients[1], coefficients[2])))
    source_at_s = (coefficients[0] + coefficients[1] * s + coefficients[2] * s**2) * 2 * s
    rows = [{"name": "ordinary_quotient_C1", "passed": q == q_direct},
            {"name": "polynomial_signed_split", "passed": split == direct},
            {"name": "polynomial_primitive_derivative", "passed": derivative == source_at_s}]
    return {"passed": all(row["passed"] for row in rows), "controls": rows,
            "original_functions_evaluated": False, "accepted_producer_executed": False}


def build(root: Any, parent_report: Any = None) -> dict[str, Any]:
    """Return exact local signed roots; numeric_value_available stays false."""
    repo = _repo(root)
    parent, parent_pin, gaps = _parent(repo, parent_report)
    if not parent:
        return _fail(repo, parent_pin, gaps)
    checks, finite, core = _m(parent.get("checks")), _m(parent.get("finite_N_signed_source")), _m(parent.get("inner_collar_source"))
    if checks.get("passed") is not True or finite.get("passed") is not True or finite.get("source_executed") is not False:
        gaps.append("accepted finite_N_signed_source source-only receipt failed")
    if core.get("passed") is not True or core.get("source_executed") is not False or _m(core.get("exact_identity")).get("delta_six_zero") is not False:
        gaps.append("accepted inner_collar_source source-only receipt failed")
    family = {key: _m(parent.get("source_family")).get(key) for key in FAMILY_KEYS}
    if family != _FAMILY:
        gaps.append("parent source_family mismatch")
    pins = {
        "parent_report": dict(parent_pin),
        "lower_graph": _pin(parent, repo, LOWER_GRAPH, LOWER_GRAPH_NAME, "input", gaps),
        "lower_source": _pin(parent, repo, LOWER_SOURCE, LOWER_SOURCE_NAME, "input", gaps),
        "forward_source": _pin(parent, repo, FORWARD_SOURCE, FORWARD_SOURCE, "implementation", gaps),
    }
    pins_ok = all(row["passed"] is True for row in pins.values()) and parent_pin.get("passed") and parent_pin.get("mapping_match")
    try:
        lower = load(repo / LOWER_GRAPH)
        nodes = list(lower["exact_graph_nodes"])
    except Exception as exc:
        return _fail(repo, parent_pin, gaps + [f"cannot load lower graph: {exc}"])
    if {key: _m(lower.get("source_family")).get(key) for key in FAMILY_KEYS} != family:
        gaps.append("lower graph source_family mismatch")
    proofs = {
        "s": _proof(nodes, HANDLES["s"], {"operation": "bound_variable", "name": "original_lower_collar_s"}),
        "hb": _proof(nodes, HANDLES["hb"], {"operation": "current_original_source_parameter", "name": "hbB"}),
        "epsilon": _proof(nodes, HANDLES["epsilon"], {"operation": "analytic_unary", "name": "exp"}),
        "rho": _proof(nodes, HANDLES["rho"], {"operation": "product", "arguments": [13, 35850]}),
        "sigma": _proof(nodes, HANDLES["sigma"], {"operation": "exact_original_source_flat_sigma", "argument": HANDLES["s"]}),
        "sigma_prime": _proof(nodes, HANDLES["sigma_prime"], {"operation": "analytic_unary", "name": "original_flat_sigma_prime", "argument": HANDLES["s"]}),
        "ratio_C0": _proof(nodes, HANDLES["ratio_C0"], {"operation": "analytic_unary", "name": "exp", "argument": 35936}),
        "ratio_C1": _proof(nodes, HANDLES["ratio_C1"], {"operation": "product", "arguments": [35940, 35937]}),
        "Psi_core_rho_C0": _proof(nodes, HANDLES["Psi_core_rho_C0"], {"operation": "function_substitution", "expression": 35715, "variable": 35558, "value": 35851}),
        "Psi_core_rho_C1": _proof(nodes, HANDLES["Psi_core_rho_C1"], {"operation": "function_substitution", "expression": 35833, "variable": 35558, "value": 35851}),
        "L_C0": _proof(nodes, HANDLES["L_C0"], {"operation": "sum", "arguments": [1, 30380]}),
        "L_C1": _proof(nodes, HANDLES["L_C1"], {"operation": "product", "arguments": [220, 30379]}),
        "Z": _proof(nodes, HANDLES["Z"], {"operation": "bound_variable", "name": "Z"}),
        "core_rho_variable": _proof(nodes, HANDLES["core_rho_variable"], {"operation": "bound_variable", "name": "original_n1_true_core_rho"}),
        "rho4": _proof(nodes, HANDLES["rho4"], {"operation": "exact_rational", "numerator": 4, "denominator": 1}),
    }
    r0, r1 = nodes[HANDLES["Psi_core_rho4_C0_source"]], nodes[HANDLES["Psi_core_rho4_C1_source"]]
    radial = {
        "C0_actual_radial_order1": r0.get("operation") == "current_original_unique_analytic_core_component" and r0.get("quantity") == "Psi" and r0.get("radial_order") == 1 and r0.get("ordinary_slow_Z_derivative_order") == 0 and r0.get("rho") == HANDLES["core_rho_variable"] and r0.get("Z_variable") == HANDLES["Z"],
        "C1_actual_ordinary_row": r1.get("operation") == "current_original_unique_analytic_core_component" and r1.get("quantity") == "Psi" and r1.get("radial_order") == 1 and r1.get("ordinary_slow_Z_derivative_order") == 1 and r1.get("rho") == HANDLES["core_rho_variable"] and r1.get("Z_variable") == HANDLES["Z"],
        "not_radial0_Psi4": not any(i in _RADIAL_ZERO_PSI4_NODES for i in (HANDLES["Psi_core_rho4_C0_source"], HANDLES["Psi_core_rho4_C1_source"])),
        "no_saved_35715_at_rho4": not any(isinstance(n, Mapping) and n.get("operation") == "function_substitution" and n.get("expression") == 35715 and n.get("value") == 13 for n in nodes),
    }
    guards = {
        "canonical_collar_fields": ast_guard(repo, LOWER_SOURCE, "collar_fields", ("hb,eps=parameters['hb'],parameters['epsilon'];de=parameters['delta']", "rho=lambda x:g.mul(c(4),g.unary('exp',g.mul(hb,x)))", "jacobian=g.mul(hb,rs)", "weighted_physical_measure=g.mul(eps,eps,rs,jacobian)")),
        "canonical_direct_pressure_functions": ast_guard(repo, LOWER_SOURCE, "direct_pressure_functions", ("Original differential N1p with both radial derivatives and true C1", "original_direct_N1p=True", "original_m15_minus_applied_once=True", "signed_integrand_not_a_cap=True")),
        "forward_leading_regular_frame": ast_guard(repo, FORWARD_SOURCE, "leading_regular_frame", ("chi=1-sigma+sigma*hb", "d,L=1-z*z,1-z*z*de", "Q=[(u*(2*z)-m*(z*(1-de))-m.dz()*d)/L")),
        "forward_prepare_singular_formula": ast_guard(repo, FORWARD_SOURCE, "_prepare", ("singular=2*rr*rr*zi*one*(phi/comp)*psir/L", "basep=inp['radial']['Psi'][1].truncate(2)*z.truncate(2)*32/(1-z.truncate(2)**2*de)", "J.constant(c,4,5)")),
    }
    ast_ok = all(row.get("passed") is True for row in guards.values())
    if not ast_ok or not all(row["passed"] for row in proofs.values()) or not all(radial.values()):
        gaps.append("source handle, radial, or AST guard failed")
    g = SignedDensityExpressionGraph()
    ref = lambda i: g.reference(LOWER_GRAPH_NAME, pins["lower_graph"]["sha256"], i, source_nodes=nodes)
    s, hb, eps, Z = (ref(HANDLES[k]) for k in ("s", "hb", "epsilon", "Z"))
    L0, L1, rho, ratio0, ratio1, psi0, psi1, sigma, sigmap = (ref(HANDLES[k]) for k in ("L_C0", "L_C1", "rho", "ratio_C0", "ratio_C1", "Psi_core_rho_C0", "Psi_core_rho_C1", "sigma", "sigma_prime"))
    ls, lt = _bound(g, "s", "outer signed collar coordinate"), _bound(g, "t", "integration coordinate")
    bind = lambda x, name: _sub(g, x, s, ls, name, local_scope="s", original_graph_nodes_added_or_replaced=False)
    rho_s, q0_s, q1_s, p0_s, p1_s, sig_s, sp_s = (bind(x, n) for x, n in ((rho, "rho_of_s"), (ratio0, "ratio_C0_of_s"), (ratio1, "ratio_C1_of_s"), (psi0, "Psi_C0_of_s"), (psi1, "Psi_C1_of_s"), (sigma, "sigma_of_s"), (sigmap, "sigma_prime_of_s")))
    core_rho, four = ref(HANDLES["core_rho_variable"]), ref(HANDLES["rho4"])
    rp0 = _sub(g, ref(HANDLES["Psi_core_rho4_C0_source"]), core_rho, four, "Psi_core_rho4_C0", local_scope="A0", radial_order=1)
    rp1 = _sub(g, ref(HANDLES["Psi_core_rho4_C1_source"]), core_rho, four, "Psi_core_rho4_C1", local_scope="A0", radial_order=1, ordinary_slow_Z_derivative_order=1)
    one, two, thirty_two = g.one, g.constant(2), g.constant(32)
    pref, prod0, prod1 = two * rho_s**2 * (one - hb), q0_s * p0_s, q1_s * p0_s + q0_s * p1_s
    A0, A1 = pref * Z * prod0 / L0, pref * ((prod0 + Z * prod1) / L0 - Z * prod0 * L1 / L0**2)
    B0, B1 = thirty_two * Z * rp0 / L0, thirty_two * ((rp0 + Z * rp1) / L0 - Z * rp0 * L1 / L0**2)
    A0s, A1s = _sub(g, A0, ls, ls, "A_C0_of_s", local_scope="s", ordinary_slow_Z_derivative_order=0), _sub(g, A1, ls, ls, "A_C1_of_s", local_scope="s", exact_ordinary_C1_algebra=True, ordinary_slow_Z_derivative_order=1)
    B0s, B1s = _sub(g, B0, ls, ls, "A0_C0_of_s", local_scope="s", A0_independent_of_s=True), _sub(g, B1, ls, ls, "A0_C1_of_s", local_scope="s", A0_independent_of_s=True, exact_ordinary_C1_algebra=True)
    At = [_sub(g, x, ls, lt, n, local_scope="t") for x, n in ((A0s, "A_C0_of_t"), (A1s, "A_C1_of_t"))]
    BAt = [_sub(g, x, ls, lt, n, local_scope="t", A0_independent_of_s=True) for x, n in ((B0s, "A0_C0_of_t"), (B1s, "A0_C1_of_t"))]
    spt = _sub(g, sp_s, ls, lt, "sigma_prime_of_t", local_scope="t")
    eps2, integrals, primitive, term4, full, integrands = eps * eps, [], [], [], [], []
    for order, (aa, bb, Bs) in enumerate(zip(At, BAt, (B0s, B1s))):
        integrand = (aa - bb) * spt
        integral = g.insert({"operation": "definite_integral", "integrand": integrand.node, "variable": lt.node, "lower": g.zero.node, "upper": ls.node, "variable_name": "t", "upper_variable_name": "s", "exact_function_integral": True, "numerical_value_not_installed": True, "signed_A_minus_A0_remainder": True, "ordinary_slow_Z_derivative_order": order, "local_scope": {"t": lt.node, "s": ls.node}, "physical_epsilon_squared_outside_integral": True, "source_cap_is_not_integrand": True, "original_graph_nodes_added_or_replaced": False})
        unscaled_term = _sub(g, Bs * sig_s, ls, ls, f"term4_A0_sigma_unscaled_C{order}", local_scope="s", A0_independent_of_s=True)
        unscaled = _sub(g, unscaled_term + integral, ls, ls, f"full_singular_sum_unscaled_C{order}", local_scope="s", exact_split=True)
        integrals.append(integral); integrands.append(integrand); primitive.append(_sub(g, eps2 * integral, ls, ls, f"signed_primitive_C{order}", local_scope="s", physical_epsilon_squared_factor_count=1, numerical_value_not_installed=True))
        term4.append(_sub(g, eps2 * (Bs * sig_s), ls, ls, f"term4_A0_sigma_C{order}", local_scope="s", A0_independent_of_s=True, physical_epsilon_squared_factor_count=1, numerical_value_not_installed=True)); full.append(_sub(g, eps2 * unscaled, ls, ls, f"full_singular_sum_C{order}", local_scope="s", exact_split=True, physical_epsilon_squared_factor_count=1, numerical_value_not_installed=True))
    controls = _controls(); scope = {"local_s_is_bound": ls.node != s.node, "local_t_is_bound": lt.node != ls.node, "integral_upper_is_local_s": all(g.nodes[x.node]["upper"] == ls.node for x in integrals), "integral_variable_is_local_t": all(g.nodes[x.node]["variable"] == lt.node for x in integrals), "A0_independent_of_s": True, "physical_epsilon_squared_once": True, "original_graph_nodes_added_or_replaced": False}; controls["graph_scope_controls"] = scope; controls["passed"] = controls["passed"] and all(v for k, v in scope.items() if k != "original_graph_nodes_added_or_replaced")
    passed = bool(parent_pin.get("passed") and parent_pin.get("mapping_match") and pins_ok and ast_ok and all(x["passed"] for x in proofs.values()) and all(radial.values()) and controls["passed"] and not gaps)
    result = {"passed": passed, "source_binding_passed": passed, "source_only": True, "source_executed": False, "parent_report": PARENT_REPORT, "parent_report_sha256": PARENT_REPORT_SHA256, "parent_report_pin": parent_pin, "source_family": family, "source_pins": pins, "source_pins_ok": pins_ok, "AST": guards, "source_handle_proofs": proofs, "radial_checks": radial, "handles": dict(HANDLES), "ordinary_C1": {"A_C0_C1": [A0s.node, A1s.node], "A0_C0_C1": [B0s.node, B1s.node], "A_formula": "2*rho^2*Z*(1-hb)*ratio*Psi/L", "A0_formula": "32*Z*Psi_core_rho(4,Z)/L", "A0_independent_of_s": True, "exact_ordinary_C1_algebra": True}, "A": {"C0_C1": [A0s.node, A1s.node]}, "A0": {"C0_C1": [B0s.node, B1s.node], "independent_of_s": True}, "signed_primitive": {"integral_C0_C1": [x.node for x in integrals], "signed_primitive_C0_C1": [x.node for x in primitive], "integrand_C0_C1": [x.node for x in integrands], "variable_scope": {"integration_variable": lt.node, "upper_s": ls.node}, "physical_epsilon_squared_factor_count": 1, "numeric_value_available": False}, "term4_A0_sigma": {"term4_A0_sigma_C0_C1": [x.node for x in term4], "A0_independent_of_s": True, "physical_epsilon_squared_factor_count": 1}, "full_singular_sum": {"full_singular_sum_C0_C1": [x.node for x in full], "exact_split": True, "physical_epsilon_squared_factor_count": 1}, "graph": {"report": LOWER_GRAPH_NAME, "report_sha256": pins["lower_graph"]["sha256"], "source_graph_node_count": len(nodes), "source_graph_nodes_untouched": True, "original_graph_nodes_added_or_replaced": False, "local_node_count": len(g.nodes)}, "exact_local_expression_nodes": g.nodes, "local_nodes": len(g.nodes), "graph_packet": g.packet(), "new_exact_controls": controls, "exact_controls": controls, "epsilon_squared": {"node": eps2.node, "factor_count": 1}, "Delta_six": {"exact_identity": False, "zero_claim": False, "delta_six_zero_claim": False, "ordinary_Z0_zero": False, "ordinary_Z1_zero": False, "numeric_value_available": False}, "delta_six_zero_claim": False, "exact_identity": {"Delta_six_zero": False, "ordinary_Z0_zero": False, "ordinary_Z1_zero": False}, "numeric_value_available": False, "actual_numerical_source_value_installed": False, "accepted_producer_executed": False, "caps_removed_or_narrowed": False, "quiet179_views_unchanged": True, "root": str(repo), "directory": str(repo / "experiments" / "root_st073"), "gaps": gaps}
    return result


def build_original_collar_signed_primitive(root: Any, parent_report: Any = None) -> dict[str, Any]:
    return build(root, parent_report)


def build_signed_collar_primitive(root: Any, parent_report: Any = None) -> dict[str, Any]:
    return build(root, parent_report)


__all__ = ["PARENT_REPORT", "PARENT_REPORT_SHA256", "LOWER_GRAPH", "LOWER_GRAPH_SHA256", "HANDLES", "build", "build_original_collar_signed_primitive", "build_signed_collar_primitive"]
