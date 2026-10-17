"""Source-only Rc profile and radial-continuity join.

This adapter reads the already saved outer and compact source graphs.  It
does not import or execute the archived owners, quadratures, or checkers.
The returned row IDs are the actual function definitions used by the saved
graphs; endpoint zero tests only certify the flat source branches.
"""
from __future__ import annotations

import ast
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from .original_finite_N_signed_source import ast_guard, digest, load, require, sha

PREFIX = "lei_ren_part1_paper_compliant_current_original_"
REPORT = PREFIX + "correlated_scaled_source.json"
REPORT_SHA256 = "9beead66814da85325e86f4ee2d9c091da34c524f9512f8cbf4c77e18bf1d444"
PARENT = PREFIX + "finite_N_signed_source.json"
OUTER = PREFIX + "true_outer_C6_C5_functions.json.gz"
COMPACT = PREFIX + "true_compact_C5_functions.json.gz"
OUTER_SOURCE = "experiments/root_st073/" + PREFIX + "true_outer_C6_C5_functions.py"
COMPACT_SOURCE = "experiments/root_st073/" + PREFIX + "true_compact_C5_functions.py"
LOOP_SOURCE = "experiments/root_st073/" + PREFIX + "prescribed_loop_C5_functions.py"
BETA_SOURCE = "experiments/root_st073/lei_ren_part1_paper_compliant_five_moment_repair.py"
X_OUTER, X_COMPACT = 250215, 260181


def _method_guard(root, relative, cls, name, required):
    source = (root / relative).read_text(encoding="utf-8")
    tree = ast.parse(source)
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == cls:
            found.extend(x for x in node.body if isinstance(x, (ast.FunctionDef, ast.AsyncFunctionDef)) and x.name == name)
    require(len(found) == 1, "Unique source method required: " + cls + "." + name)
    body = "".join((ast.get_source_segment(source, found[0]) or "").split())
    checks = {text: "".join(text.split()) in body for text in required}
    require(all(checks.values()), "Changed source method: " + cls + "." + name)
    return dict(passed=True, path=relative, class_name=cls, function=name,
                required_text=checks,
                AST_sha256=hashlib.sha256(ast.dump(found[0], include_attributes=False).encode()).hexdigest(),
                source_executed=False)


class _Poly:
    """Tiny exact polynomial in L=log(2), only for compact support tests."""
    def __init__(self, terms=()):
        self.terms = {int(k): Fraction(v) for k, v in dict(terms).items() if v}

    @staticmethod
    def c(value):
        return _Poly({0: Fraction(value)})

    def __add__(self, other):
        out = dict(self.terms)
        for k, v in other.terms.items():
            out[k] = out.get(k, 0) + v
        return _Poly(out)

    def __neg__(self):
        return _Poly({k: -v for k, v in self.terms.items()})

    def __mul__(self, other):
        out = {}
        for k, v in self.terms.items():
            for j, w in other.terms.items():
                out[k + j] = out.get(k + j, 0) + v * w
        return _Poly(out)

    def constant(self):
        return self.terms.get(0) if set(self.terms) <= {0} else None

    def zero(self):
        return not self.terms

    def monomial(self):
        return next(iter(self.terms.items())) if len(self.terms) == 1 else None


def _poly(value):
    return value if isinstance(value, _Poly) else _Poly.c(value)


def _divide(num, den):
    num, den = _poly(num), _poly(den)
    mono = den.monomial()
    if mono is None or mono[1] == 0:
        return None
    power, coeff = mono
    return _Poly({k - power: v / coeff for k, v in num.terms.items()})


class _GraphProof:
    """Evaluate only exact constants, support predicates, and zero factors."""
    def __init__(self, nodes, kind):
        self.nodes, self.kind, self.lazy, self.memo = nodes, kind, [], {}

    def ev(self, node, env=None):
        env = {} if env is None else env
        key = (node, tuple(sorted((k, repr(v.terms) if isinstance(v, _Poly) else repr(v)) for k, v in env.items())))
        if key in self.memo:
            return self.memo[key]
        value = self._ev(node, env)
        self.memo[key] = value
        return value

    def _ev(self, node, env):
        if node in env:
            return env[node]
        row = self.nodes[node]
        op = row.get("operation")
        if op == "exact_rational":
            return _Poly.c(Fraction(row["numerator"], row["denominator"]))
        if op in ("bound_variable", "current_original_source_parameter", "original_compact_unique_Banach_limit_component",
                  "original_Rm_unique_leading_implicit_control", "mathematical_pi"):
            return None
        if op == "function_substitution":
            value = self.ev(row["value"], env)
            new = dict(env); new[row["variable"]] = value
            return self.ev(row["expression"], new)
        if op == "negative":
            value = self.ev(row["argument"], env)
            return None if value is None else -value
        if op in ("sum", "product"):
            values = [self.ev(i, env) for i in row["arguments"]]
            if op == "product":
                if any(v is not None and v.zero() for v in values):
                    return _Poly()
                if any(v is None for v in values):
                    return None
                out = _Poly.c(1)
                for value in values: out = out * value
                return out
            if any(v is None for v in values):
                return None
            out = _Poly()
            for value in values: out = out + value
            return out
        if op in ("positive_quotient", "quotient"):
            num = self.ev(row["numerator"], env)
            if num is not None and num.zero():
                return _Poly()
            den = self.ev(row["denominator"], env)
            return None if num is None or den is None else _divide(num, den)
        if op == "analytic_unary":
            value = self.ev(row["argument"], env)
            if value is None: return None
            if row.get("name") == "log":
                q = value.constant()
                if q == 1: return _Poly()
                if q == 2: return _Poly({1: 1})
            if row.get("name") == "exp" and value.zero(): return _Poly.c(1)
            return None
        if op == "compact_raw_beta":
            value = self.ev(row["argument"], env)
            q = None if value is None else value.constant()
            if q is not None and abs(q) >= 1:
                return _Poly()
            return None
        if op == "logical_or":
            values = [self.ev(i, env) for i in row["arguments"]]
            return True if any(v is True for v in values) else (False if all(v is False for v in values) else None)
        if op == "exact_real_comparison":
            if self.kind == "outer" and node == 251222:
                return True
            left, right = self.ev(row["left"], env), self.ev(row["right"], env)
            left, right = (None if left is None else left.constant()), (None if right is None else right.constant())
            if left is None or right is None: return None
            return {"<=": left <= right, ">=": left >= right, "<": left < right, ">": left > right}[row["operator"]]
        if op == "original_lazy_flat_branch":
            self.lazy.append(node)
            require(self.ev(row["flat_predicate"], env) is True, "Flat branch predicate not proved at endpoint")
            return self.ev(row["flat_value"], env)
        return None


def _zeros(proof, nodes, env=None):
    values = [proof.ev(i, env) for i in nodes]
    require(all(v is not None and v.zero() for v in values), "Actual function row did not reduce to exact zero")
    return True


def _wrappers(graph_nodes, wrapper_ids, expressions, variable, value):
    return all(graph_nodes[i].get("operation") == "function_substitution" and graph_nodes[i].get("expression") == expressions[j]
               and graph_nodes[i].get("variable") == variable and graph_nodes[i].get("value") == value for j, i in enumerate(wrapper_ids))


def _pin(root, saved, relative):
    path = root / relative
    expected = saved["input_hashes"].get(Path(relative).name)
    require(expected is not None, "Saved direct input pin missing: " + relative)
    actual = sha(path)
    require(actual == expected, "Pinned direct input changed: " + relative)
    return dict(path=relative, expected_sha256=expected, actual_sha256=actual, passed=True)


def build(root, parent=None):
    """Build a source-only actual Rc profile join for the root pressure adapter."""
    root = Path(root); work = root / "experiments/root_st073"
    report_path = work / REPORT
    require(sha(report_path) == REPORT_SHA256, "Saved current correlated source report changed")
    saved = json.loads(report_path.read_text(encoding="utf-8"))
    pins = {"saved_report": dict(path="experiments/root_st073/" + REPORT, expected_sha256=REPORT_SHA256, actual_sha256=REPORT_SHA256, passed=True),
        "parent_report": _pin(work, {"input_hashes": {PARENT: saved["parent_report_sha256"]}}, PARENT),
        "outer_archive": _pin(root, saved, "experiments/root_st073/" + OUTER), "compact_archive": _pin(root, saved, "experiments/root_st073/" + COMPACT),
        "outer_source": _pin(root, saved, OUTER_SOURCE), "compact_source": _pin(root, saved, COMPACT_SOURCE),
        "prescribed_loop_source": _pin(root, saved, LOOP_SOURCE), "five_moment_beta_source": _pin(root, saved, BETA_SOURCE)}
    identity = saved["compact_inlet_source_identity"]
    require(identity.get("passed") and identity.get("source_only") and not identity.get("owners_replayed"),
            "Accepted 39-row compact identity receipt required")
    require(identity["identity_count"]["total"] == 39, "Accepted compact identity must retain all 39 rows")
    outer, compact = load(work / OUTER), load(work / COMPACT)
    on, cn = outer["exact_graph_nodes"], compact["exact_graph_nodes"]
    ob, cb = outer["true_outer_C6_C5_bindings"], compact["true_compact_C5_bindings"]
    family = saved["source_family"]
    require(outer["source_family"] == family and compact["source_family"] == family, "Source family changed")
    require(cn[:len(on)] == on, "Compact graph does not retain the actual outer source prefix")
    prefix_digest = digest(on)
    require(prefix_digest == identity["source_graph"]["compact_prefix_digest"], "Actual source prefix digest changed")
    require(cn[X_COMPACT]["operation"] == "bound_variable" and on[X_OUTER]["operation"] == "bound_variable",
            "Rc/outer phase variables changed")
    guards = {
        "compact_bumps": ast_guard(root, COMPACT_SOURCE, "bump_rows",
            ("leading.raw_beta(g,argument)", "first,second=leading.beta_first_two(g,argument,raw)",
             "byy=g.mul(scale,g.add")),
        "compact_profiles": ast_guard(root, COMPACT_SOURCE, "profiles",
            ("for j in range(3):", "dV=a.mul(a.trim(A),G)", "V=dV")),
        "compact_join": _method_guard(root, COMPACT_SOURCE, "CurrentOriginalTrueCompactC5Functions", "__init__",
            ("self.modified_fields.update(E=f[0]['E'],V=f[0]['V'])", "self.radial=dict(E_y=f[1]['E'],V_y=f[1]['V']")),
        "outer_geometry": ast_guard(root, OUTER_SOURCE, "geometry",
            ("domains={'Rh_reference':(c(-5),g.zero),'O2_buffer':(g.zero,c(11)),'O3_power':(g.zero,c(2))}",)),
        "outer_join": _method_guard(root, OUTER_SOURCE, "CurrentOriginalTrueOuterC6C5Functions", "__init__",
            ("V=a.fixed(0)", "loop=engine.loop_functions", "finite=engine.finite_N_functions",
             "modified.update(E=finite['E_N'],V=finite['V_N'])", "total=radial.total_radial_functions")),
        "prescribed_loop": ast_guard(root, LOOP_SOURCE, "loop_functions",
            ("shear=g.node('exact_real_comparison',operator='>=',left=roots['Delta'].value.node,right=eta.node)",
             "lazy=lambda jet:a.lazy(jet,predicate,roots['Delta'].value,eta")),
        "finite_N": ast_guard(root, LOOP_SOURCE, "finite_N_functions",
            ("V_N=a.add(V,s)", "delta_V=s")),
    }
    chart = ob["charts"]["O3_power"]; roots = chart["roots"]; loop = chart["loop"]
    require(roots["V"] == [0] * 6 and roots["b"] == [0] * 6 and roots["Delta"][0] == 251216,
            "O3 source roots are not the quiet V=0 branch")
    require(loop["flat_predicate"] == 251222 and on[251222] == {"operation": "exact_real_comparison", "operator": ">=", "left": 251216, "right": ob["eta"]}
        and ob["mu_proof"]["positive_mu_and_power_2mu_ge_eta_and_Rp_reservation"],
            "Saved original O3 flat admission missing")
    # Check the graph factors rather than trusting the descriptive quiet label.
    require(roots["a"][0] == 251146 and on[251146] == {
        "operation": "positive_quotient", "numerator": 251140, "denominator": 250226,
        "source_positive_certificate": "same strictly positive actual prescribed swirl E"},
        "O3 a source factor changed")
    require(on[251140]["operation"] == "sum" and on[251140]["arguments"] == [250226, 251134]
        and on[251134] == {"operation": "negative", "argument": 251128}
        and on[251128] == {"operation": "product", "arguments": [250233, 48]}
        and on[250233] == {"operation": "product", "arguments": [250226, 4994]}
        and on[4994] == {"operation": "negative", "argument": 4993}
        and on[4993] == {"operation": "sum", "arguments": [28, 17]}
        and on[28] == {"operation": "exact_rational", "numerator": 1, "denominator": 2}
        and ob["mu"] == 17,
        "O3 a is not the exact 2+2*mu source factor")
    require(on[roots["Delta"][0]] == {"operation": "sum", "arguments": [251146, 25304]}
        and on[25304] == {"operation": "negative", "argument": 48}
        and on[48] == {"operation": "exact_rational", "numerator": 2, "denominator": 1},
        "O3 Delta is not the exact a-2=2*mu source factor")
    outer_proof = _GraphProof(on, "outer")
    compact_proof = _GraphProof(cn, "compact")
    outer_full = ob["endpoints"]["O3_power"]["modified_fields"]["V"]
    outer_full_y = ob["endpoints"]["O3_power"]["radial"]["V_y"]
    baseline_outer = ob["background_endpoints"]["O3_power"]["fields"]["V"]
    baseline_outer_y = ob["background_endpoints"]["O3_power"]["radial"]["V_y"]
    require(_zeros(outer_proof, outer_full) and _zeros(outer_proof, outer_full_y), "Outer actual V/V_y not flat")
    require(_zeros(outer_proof, baseline_outer) and _zeros(outer_proof, baseline_outer_y), "Outer baseline V/V_y not zero")
    for node in set(outer_proof.lazy):
        row = on[node]
        require(row.get("flat_predicate") == 251222 and row.get("flat_value") == 0
                and row.get("active_body_not_evaluated_on_flat_branch") is True,
                "Outer endpoint reached a changed lazy branch")
    rows = cb["profiles"]["radial_y_rows"]
    bump_values = []
    for row in cb["profiles"]["bumps"]:
        for ref in row:
            require(compact_proof.ev(ref, {X_COMPACT: _Poly.c(1)}).zero(), "Compact bump is not flat at x=1")
        bump_values.append(list(row))
    compact_full = cb["inlet"]["modified_fields"]["V"]
    compact_full_y = cb["inlet"]["radial"]["V_y"]
    compact_baseline = cb["inlet"]["background_fields"]["V"]
    compact_baseline_y = cb["background"]["radial"]["V_y"]
    env1 = {X_COMPACT: _Poly.c(1)}
    raw_support = [{"node": i, "argument_node": cn[i]["argument"], "argument_at_x1": str(compact_proof.ev(cn[i]["argument"], env1).constant()), "exact_definition": cn[i]["exact_definition"], "source_function": cn[i]["original_source"]["function"]} for i in range(len(on), len(cn)) if cn[i].get("operation") == "compact_raw_beta"]
    require(_zeros(compact_proof, compact_full, env1) and _zeros(compact_proof, compact_full_y, env1),
            "Compact actual V/V_y not flat at x=1")
    require(_zeros(compact_proof, compact_baseline, env1) and _zeros(compact_proof, compact_baseline_y, env1),
            "Compact baseline V/V_y not zero at x=1")
    for node in set(compact_proof.lazy):
        row = cn[node]
        require(row.get("flat_value") == 0 and row.get("active_body_not_evaluated_on_flat_branch") is True,
                "Compact endpoint reached a changed lazy branch")
    profile_rows = {str(j): {key: list(rows[j][key]) for key in ("F", "G", "V")} for j in range(3)}
    for j in range(3):
        require(_zeros(compact_proof, rows[j]["F"], env1) and _zeros(compact_proof, rows[j]["G"], env1),
                "Compact F/G profile row is not flat at x=1")
    env2 = {X_COMPACT: _Poly.c(2)}; upper_density = {key: _zeros(compact_proof, refs, env2) for key, refs in cb["density"].items()}
    upper_profiles = {j: {key: _zeros(compact_proof, rows[j][key], env2) for key in ("F", "G", "V")} for j in (0, 1)}
    upper_full_zero = _zeros(compact_proof, cb["terminal"]["modified_fields"]["V"], env2) and _zeros(compact_proof, cb["terminal"]["radial"]["V_y"], env2); upper_base_zero = _zeros(compact_proof, cb["terminal"]["background_fields"]["V"], env2) and _zeros(compact_proof, cb["background"]["radial"]["V_y"], env2)
    terminal_ports = {}
    for key, cid in cb["terminal_certificates"].items():
        cert = cn[cid]; require(cert["operation"] == "proved_original_compact_relative_C5_terminal_identity" and cert["source_system"] == cb["system"] and cert["relative_correction_only"] and cert["absolute_leading_histories_and_pressure_not_zeroed"] and cert["zero_rows"] == [0] * 6 and cert["actual_endpoint_rows"] == cb["terminal"]["correction"][key], "Saved C5 terminal certificate changed")
        terminal_ports[key] = dict(node=cid, actual_endpoint_rows=cert["actual_endpoint_rows"], residual_rows=cert["equivalent_residual_rows"], source_system=cert["source_system"], zero_rows=cert["zero_rows"])
    upper_delta_f = all(upper_density.values()) and all(all(v.values()) for v in upper_profiles.values()) and all(v["zero_rows"] == [0] * 6 for v in terminal_ports.values())
    require(upper_delta_f and upper_full_zero and upper_base_zero, "Compact upper exact zero profile/terminal proof failed")
    normalization = {
        "compact_full_V_is_profile_row0": _wrappers(cn, compact_full, rows[0]["V"], X_COMPACT, 1),
        "compact_full_V_y_is_profile_row1": _wrappers(cn, compact_full_y, rows[1]["V"], X_COMPACT, 1),
        "compact_baseline_V_is_zero_function": _wrappers(cn, compact_baseline, [0] * 7, X_COMPACT, 1),
        "outer_full_V_is_finite_N_V": _wrappers(on, outer_full, chart["finite_N"]["V_N"], X_OUTER, 48),
        "outer_full_V_y_is_total_radial": _wrappers(on, outer_full_y, chart["radial"]["V_y"], X_OUTER, 48),
        "outer_baseline_V_is_background": _wrappers(on, baseline_outer, [0] * 7, X_OUTER, 48),
        "outer_baseline_V_y_is_background": _wrappers(on, baseline_outer_y, [0] * 7, X_OUTER, 48),
    }
    require(all(normalization.values()), "Full/baseline rows are not bound to actual source functions")
    return {
        "passed": True, "source_only": True, "source_executed": False,
        "owners_replayed": False, "quadratures_replayed": False, "numeric_caps_used": False,
        "source_family": family, "source_pins": pins, "parent_report_sha256": REPORT_SHA256, "input_source_prefix_sha256": prefix_digest,
        "accepted_compact_identity_receipt": dict(passed=True, source_only=True, identity_count=identity["identity_count"]),
        "AST": guards,
        "full_V_and_V_y_exact_zero": True, "baseline_V_and_V_y_exact_zero": True,
        "full_profile_radial_continuity": True,
        "whole_O3_power_quiet": True, "compact_upper_delta_F_C0_C1_exact_zero": upper_delta_f, "compact_upper_profile_V_and_V_y_exact_zero": all(upper_profiles[j]["V"] for j in (0, 1)), "compact_upper_full_V_and_V_y_exact_zero": upper_full_zero, "compact_upper_baseline_V_and_V_y_exact_zero": upper_base_zero,
        "actual_functions": {
            "compact_profile_F_G_V": profile_rows,
            "compact_inlet_V": list(compact_full), "compact_inlet_V_y": list(compact_full_y),
            "outer_O3_power_V": list(outer_full), "outer_O3_power_V_y": list(outer_full_y),
        },
        "baseline_functions": {"compact_V": list(compact_baseline), "compact_V_y": list(compact_baseline_y),
                               "outer_V": list(baseline_outer), "outer_V_y": list(baseline_outer_y)},
        "ordinary_orders": {"V": [0, 1, 2], "V_y": [0, 1], "profile_F_G": [0, 1, 2]},
        "normalization_proof": normalization,
        "outer_O3_power_quiet_branch": dict(V_roots=roots["V"], V_y_roots=ob["background"]["O3_power"]["radial"]["V_y"],
            a_root=roots["a"][0], b_roots=roots["b"], Delta_root=roots["Delta"][0],
            flat_predicate=loop["flat_predicate"], lazy_nodes=list(sorted(set(outer_proof.lazy))),
            exact_source_identity="V=V_y=0; a=2+2*mu; b=0; Delta=2*mu; 2*mu>=eta"),
        "compact_flat_bump_nodes": bump_values,
        "compact_raw_beta_support": raw_support, "compact_upper_terminal_c5_ports": terminal_ports,
        "remaining_gaps": ["Pressure-boundary pairing is owned by the root adapter."],
    }
