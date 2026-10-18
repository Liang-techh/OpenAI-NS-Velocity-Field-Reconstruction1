"""Actual reference/restoration/post-restore native C2 endpoint identities.

Saved sigma endpoint jets are consumed as a lemma.  Only the two new
native endpoint substitutions are calculated here; no accepted chart is
rebuilt or numerically integrated.
"""
from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_long_reference_modified_seams.json"
PARENT_SHA256 = "204d8b45aff263a3df1eae319dd7534e304162794ac375a1c0ec2cfe9beafc35"
C6 = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions"
MODULE = "src/openai_ns_reconstruction/original_restoration_native_endpoints.py"
SIGMA_SOURCE = BASE + "lei_ren_part1_paper_compliant_outer_initial.py"
SIGMA_JETS = BASE + "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"
CHARTS = ("reference", "restoration", "post_restore")
TOKENS = (
    "self.E=E=a.sub(self.endpoints['R110']['actual']['U'],a.scale(z,4))",
    "t=g.symbol('true_prescribed_restoration_t');offset=g.add(c(-8),t)",
    "U=a.add(a.scale(z,4),a.scale(E,g.sub(g.one,sigma(g,t))))",
    "U_y=a.scale(E,g.neg(sigma_prime(g,t)))",
    "self.add_chart('restoration',t,R,g.one,actual,self.endpoints['Rz']['moments'])",
    "p=g.symbol('true_prescribed_post_restore_fraction');offset=g.add(c(-7),p)",
    "actual=dict(phi=phi,U=a.scale(z,4),phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0)",
    "self.add_chart('post_restore',p,R,g.one,actual,self.endpoints['Rrestore']['moments'])",
    "self.logref=logref=g.mul(c(10),g.add(self.core.parameters['logC'],self.core.parameters['logP']))",
    "logUtheta=a.add(a.neg(logq),a.fixed(g.sub(g.mul(c('1/10'),y),self.core.parameters['logC'])))",
    "logUtheta=a.add(a.neg(logq),a.fixed(g.add(self.core.parameters['logP'],g.mul(c('1/10'),offset))))",
    "R=g.mul(c(110),g.unary('exp',g.add(logref,offset)))",
    "Utheta=a.exp(logUtheta);phi=self.phi_from_Utheta(Utheta,R)",
)


def require(value: Any, message: str) -> None:
    if not value:
        raise ValueError(message)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational(nodes: list[dict], node: int, numerator: int, denominator: int = 1) -> bool:
    row = nodes[node]
    return (row.get("operation") == "exact_rational"
            and row.get("numerator") == numerator and row.get("denominator") == denominator)


def product(nodes: list[dict], node: int, left: int, right: int) -> bool:
    row = nodes[node]
    return row.get("operation") == "product" and sorted(row["arguments"]) == sorted((left, right))


def depends(nodes: list[dict], root: int, variable: int) -> bool:
    """Conservative actual-handle reachability, without name substitution."""
    pending, seen = [root], set()
    fields = ("arguments", "argument", "numerator", "denominator", "expression",
              "variable", "value", "integrand", "lower", "upper", "rho", "Z_variable")
    while pending:
        node = pending.pop()
        if node == variable:
            return True
        if node in seen:
            continue
        seen.add(node)
        row = nodes[node]
        for key in fields:
            value = row.get(key)
            values = value if isinstance(value, list) else (value,)
            pending.extend(item for item in values if type(item) is int and 0 <= item < len(nodes))
    return False


def inlet_chain(nodes: list[dict], start: int, target: int) -> list[dict[str, Any]]:
    chain = []
    while start != target:
        row = nodes[start]
        require(row.get("operation") == "function_substitution", "Native constant-U inlet chain changed")
        variable, body, value = row["variable"], row["expression"], row["value"]
        require(rational(nodes, value, 1) and row.get("Z_independent_substitution") is True
                and not depends(nodes, body, variable), "Constant-U substitution has a moving native input")
        chain.append(dict(node=start, expression=body, variable=variable, value=value,
                          native_coordinate_absent=True, passed=True))
        start = body
    return chain


def symbolic_endpoints() -> dict[str, Any]:
    import sympy as s
    z, T, logC, logP, logF0, ref, t, p = s.symbols("Z T logC logP logF0 q t p", real=True)
    U0 = s.Function("U0")(z)
    logref = 10 * (logC + logP)
    length = logref - T - 8
    yref = T + length * ref
    logRref = s.log(110) + yref
    logRrest = s.log(110) + logref - 8 + t
    logRpost = s.log(110) + logref - 7 + p
    logtheta_ref = yref / 10 - logC - s.log(1 + z*z)
    logtheta_rest = logP + (-8 + t) / 10 - s.log(1 + z*z)
    logtheta_post = logP + (-7 + p) / 10 - s.log(1 + z*z)
    logphis = (logtheta_ref - logF0 - (s.log(2) + logRref)/2,
               logtheta_rest - logF0 - (s.log(2) + logRrest)/2,
               logtheta_post - logF0 - (s.log(2) + logRpost)/2)
    phase = (ref, t, p)
    J = (length, s.Integer(1), s.Integer(1))
    rows = []
    sigma = s.Function("sigma")(t)
    Us = (U0, 4*z + (U0-4*z)*(1-sigma), 4*z)

    def U_endpoint(index: int, endpoint: int, order: int):
        expression = s.diff(Us[index], phase[index], order) / J[index]**order
        if index == 1:
            expression = expression.subs({s.diff(sigma, t, 2): 0,
                                           s.diff(sigma, t): 0, sigma: endpoint})
        return s.simplify(expression.subs(phase[index], endpoint))

    for left, right, name in ((0, 1, "reference_to_restoration"),
                              (1, 2, "restoration_to_postrestore")):
        lp, rp = logphis[left], logphis[right]
        phi_ratio_log = s.simplify(lp.subs(phase[left], 1) - rp.subs(phase[right], 0))
        logR = (logRref, logRrest, logRpost)
        radius_difference = s.simplify(logR[left].subs(phase[left], 1)
                                       - logR[right].subs(phase[right], 0))
        slopes = [s.simplify(s.diff(logphis[k], phase[k]) / J[k]) for k in (left, right)]
        residuals = dict(log_phi=phi_ratio_log, log_R=radius_difference,
                         U=s.simplify(U_endpoint(left, 1, 0)-U_endpoint(right, 0, 0)),
                         U_y=s.simplify(U_endpoint(left, 1, 1)-U_endpoint(right, 0, 1)),
                         U_yy=s.simplify(U_endpoint(left, 1, 2)-U_endpoint(right, 0, 2)),
                         phi_y_over_phi=s.simplify(slopes[0]-slopes[1]),
                         phi_yy_over_phi=s.simplify(slopes[0]**2-slopes[1]**2))
        ordinary = {key: [str(s.simplify(s.diff(value, z, n))) for n in range(4)]
                    for key, value in residuals.items()}
        require(all(value == "0" for values in ordinary.values() for value in values),
                "Actual new endpoint function identity failed")
        rows.append(dict(name=name, left_chart=CHARTS[left], right_chart=CHARTS[right],
                         residuals={key: str(value) for key, value in residuals.items()},
                         ordinary_Z0_Z3_residuals=ordinary, passed=True))
    return dict(passed=True, seams=rows, phi_y_over_phi="-2/5", phi_yy_over_phi="4/25",
                E_y_over_E="1/10", E_yy_over_E="1/100", V_y="0", V_yy="0",
                formulas=dict(logphi=[str(x) for x in logphis], logref=str(logref)),
                function_identity_not_sampling=True)


def build(root: str | Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    root = Path(root)
    raw_parent = (root / PARENT).read_bytes()
    require(sha(root / PARENT) == PARENT_SHA256 and json.loads(raw_parent) == parent,
            "Exact modified-seams predecessor required")
    pins = {}
    for path in (PARENT, C6 + ".py", C6 + ".json.gz", SIGMA_SOURCE, SIGMA_JETS, MODULE):
        actual = sha(root / path)
        expected = (PARENT_SHA256 if path == PARENT else actual if path == MODULE
                    else parent["input_hashes"][Path(path).name])
        require(actual == expected, "Changed source: " + path)
        pins[path] = dict(path=path, expected_sha256=expected, actual_sha256=actual,
                          passed=True, source_executed=False)
    source = (root / (C6 + ".py")).read_text(encoding="utf-8")
    tree = ast.parse(source)
    cls = next(node for node in tree.body if isinstance(node, ast.ClassDef)
               and node.name == "CurrentOriginalPrescribedR100RmC6Functions")
    constructor = next(node for node in cls.body if isinstance(node, ast.FunctionDef) and node.name == "__init__")
    text = ast.get_source_segment(source, constructor)
    guards = [dict(token=token, passed=token in text) for token in TOKENS]
    require(all(item["passed"] for item in guards), "Actual restoration source formula changed")
    cache = json.loads(gzip.decompress((root / (C6 + ".json.gz")).read_bytes()))
    require(cache["source_family"] == parent["source_family"], "Wrong C6 source family")
    nodes = cache["exact_graph_nodes"]
    data = cache["prescribed_R100_Rm_C6_bindings"]
    sigma = parent["actual_long_reference_native_seams"]["sigma_endpoint_identities"]
    require(sigma["passed"] is True and sigma["squared_flat_cutoff"] is True
            and sigma["endpoint_rows"] == {"0": ["0", "0", "0"], "1": ["1", "0", "0"]}
            and sigma["derivative_orders"] == [0, 1, 2], "Saved squared-flat sigma lemma changed")
    target_U = data["endpoints"]["R110"]["actual"]["U"]
    reference_U = data["charts"]["reference"]["actual"]["U"]
    chains = [inlet_chain(nodes, start, target) for start, target in zip(reference_U, target_U)]
    require(len(chains) == 7, "Genuine native U ordinary rows missing")
    E0 = data["E"][0]
    rest = data["charts"]["restoration"]
    post = data["charts"]["post_restore"]
    post_U0 = post["actual"]["U"][0]
    erow = nodes[E0]
    require(erow["operation"] == "sum" and target_U[0] in erow["arguments"]
            and any(nodes[n].get("operation") == "negative" and nodes[n]["argument"] == post_U0
                    for n in erow["arguments"]), "E is not actual U_R110 minus 4Z")
    require(rest["physical_measure"] == post["physical_measure"] == 1,
            "Restoration physical Jacobian is not one")
    # Bind the actual U and U_y graph bodies before substituting flat endpoint jets.
    u_rest = nodes[rest["actual"]["U"][0]]
    require(u_rest.get("operation") == "sum" and post_U0 in u_rest["arguments"],
            "Restoration actual U lacks its 4Z term")
    correction = next((n for n in u_rest["arguments"] if n != post_U0), None)
    require(type(correction) is int and nodes[correction].get("operation") == "product"
            and E0 in nodes[correction]["arguments"], "Restoration actual U lacks actual E")
    alpha = next(n for n in nodes[correction]["arguments"] if n != E0)
    alpha_row = nodes[alpha]
    require(alpha_row.get("operation") == "sum" and 1 in alpha_row["arguments"], "Actual alpha lacks one")
    neg_sigma = next(n for n in alpha_row["arguments"] if n != 1)
    require(nodes[neg_sigma].get("operation") == "negative", "Actual alpha lacks negative sigma")
    sigma_row = nodes[nodes[neg_sigma]["argument"]]
    require(sigma_row.get("operation") == "exact_original_source_flat_sigma"
            and sigma_row.get("argument") == rest["phase"]
            and sigma_row.get("exact_function_not_its_interval_cap") is True
            and sigma_row.get("source_binding", {}).get("module_sha256") == pins[SIGMA_SOURCE]["actual_sha256"],
            "Actual U uses a different sigma")
    uy_row = nodes[rest["actual"]["U_y"][0]]
    require(uy_row.get("operation") == "product" and E0 in uy_row["arguments"], "Actual U_y lacks E")
    neg_prime = next(n for n in uy_row["arguments"] if n != E0)
    require(nodes[neg_prime].get("operation") == "negative", "Actual U_y lacks minus sigma prime")
    prime_row = nodes[nodes[neg_prime]["argument"]]
    require(prime_row.get("operation") == "analytic_unary"
            and prime_row.get("name") == "original_flat_sigma_prime"
            and prime_row.get("argument") == rest["phase"]
            and prime_row.get("derivative_source", {}).get("module_sha256") == pins[SIGMA_JETS]["actual_sha256"]
            and prime_row.get("exact_flat_source_derivative_not_a_bound") is True,
            "Actual U_y uses a different derivative source")
    for chart in ("reference", "post_restore"):
        require(all(rational(nodes, node, 0) for node in data["charts"][chart]["actual"]["U_y"]),
                "Constant-U native first y rows are not zero")
    source_rows = {}
    for chart in CHARTS:
        actual = data["charts"][chart]["actual"]
        phi, phi_y = actual["phi"], actual["phi_y"]
        require(len(phi) == len(phi_y) == 7, "Native phi ordinary rows missing")
        checks = []
        for value, derivative in zip(phi, phi_y):
            row = nodes[derivative]
            factors = row.get("arguments", [])
            checks.append(row.get("operation") == "product" and len(factors) == 2
                          and value in factors and any(rational(nodes, n, -2, 5) for n in factors))
        require(all(checks), "Actual native phi_y differs from -2phi/5")
        phi0 = nodes[phi[0]]
        theta = actual["physical_Utheta"][0]
        require(phi0.get("operation") == "positive_quotient" and phi0.get("numerator") == theta,
                "Actual phi does not use its source Utheta")
        denominator = nodes[phi0["denominator"]]
        require(denominator.get("operation") == "product" and data["F0"][0] in denominator["arguments"],
                "Actual phi lacks the common F0")
        sqrt_id = next(n for n in denominator["arguments"] if n != data["F0"][0])
        sqrt = nodes[sqrt_id]
        require(sqrt.get("operation") == "analytic_unary" and sqrt.get("name") == "sqrt",
                "Actual phi lacks sqrt(2R)")
        radicand = nodes[sqrt["argument"]]
        require(radicand.get("operation") == "product" and data["charts"][chart]["radius"] in radicand["arguments"]
                and any(rational(nodes, n, 2) for n in radicand["arguments"]), "Actual phi radius changed")
        source_rows[chart] = dict(phase=data["charts"][chart]["phase"],
                                  radius=data["charts"][chart]["radius"],
                                  J=data["charts"][chart]["physical_measure"],
                                  actual=actual, actual_phi_y_rows_checked=checks)
    symbolic = symbolic_endpoints()
    return dict(passed=True, parent_report=PARENT, parent_report_sha256=PARENT_SHA256,
                source_family=parent["source_family"], source_pins=pins, source_guards=guards,
                sourceguards_tuple=True, source_executed=False, accepted_ancestors_replayed=False,
                sigma_endpoint_lemma=dict(consumed_saved_only=True, **sigma),
                actual_source_rows=source_rows, reference_U_R110_substitution_chains=chains,
                actual_E_handle=E0, actual_U_R110=target_U, actual_4Z_handle=post_U0,
                symbolic_endpoint_identities=symbolic, native_phi_U_physical_C2=True,
                radius_endpoint_identities=True, ordinary_Z_orders=[0, 1, 2, 3],
                axial_scope=dict(Z_interval=["37/100", "93/250"]),
                limits=dict(full_modified_Qy=False, current_N_numeric=False,
                            global_five_moment_closure=False, full_goal_complete=False))
