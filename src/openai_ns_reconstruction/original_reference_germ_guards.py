"""Source-only guards for the reference/restoration physical germ argument.

This module deliberately inspects the accepted source graph and the pinned
source files only.  It does not construct a producer, replay a numerical
ancestor, or turn a range/receipt into a function value.  ``passed`` means
that the listed source premises were found; it is not an inlet seam theorem.
The global periodic lift is reported separately because the prescribed source
encodes the fundamental inverse interval, while its ``Phi(theta+2*pi)`` lift
is supplied by a separate theorem.
"""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path
from typing import Any, Iterable, Mapping


_SOURCE_FILES = {
    "modified": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py",
    "outer": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py",
    "loop": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py",
    "bridge": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py",
    "bridge_integrals": "experiments/root_st073/lei_ren_part1_paper_compliant_actual_bridge_integrals.py",
    "inlet": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py",
    "control": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection.py",
}


def _unwrap(value: Any, key: str) -> Mapping[str, Any]:
    if isinstance(value, Mapping) and key in value and isinstance(value[key], Mapping):
        return value[key]
    return value if isinstance(value, Mapping) else {}


def _handle_index(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    node = getattr(value, "node", None)
    return node if isinstance(node, int) else None


def _node(nodes: Any, index: Any) -> Mapping[str, Any] | None:
    i = _handle_index(index)
    if i is None:
        return None
    try:
        row = nodes[i] if not isinstance(nodes, Mapping) else nodes.get(i, nodes.get(str(i)))
    except (IndexError, KeyError, TypeError):
        return None
    return row if isinstance(row, Mapping) else None


def _repo_root(root: Any) -> Path:
    if root is not None:
        p = Path(root)
        if p.is_file():
            p = p.parent
        if (p / "experiments").is_dir() or (p / "src").is_dir():
            return p
    return Path(__file__).resolve().parents[2]


def _resolve(root: Path, relative: str) -> Path | None:
    candidates = [root / relative, root / Path(relative).name,
                  root / "experiments" / "root_st073" / Path(relative).name]
    for path in candidates:
        if path.is_file():
            return path
    return None


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _pin_maps(code_hashes: Mapping[str, str] | None, archives: Iterable[Any]) -> list[Mapping[str, str]]:
    maps: list[Mapping[str, str]] = []
    if isinstance(code_hashes, Mapping):
        maps.append(code_hashes)
    for archive in archives:
        if not isinstance(archive, Mapping):
            continue
        for key in ("input_hashes", "implementation_hashes"):
            value = archive.get(key)
            if isinstance(value, Mapping):
                maps.append(value)
        receipt = archive.get("_receipt")
        if isinstance(receipt, Mapping) and isinstance(receipt.get("input_hashes"), Mapping):
            maps.append(receipt["input_hashes"])
    return maps


def _expected_pin(relative: str, path: Path, maps: Iterable[Mapping[str, str]]) -> tuple[str | None, str | None]:
    names = (relative, str(Path(relative)), Path(relative).name,
             "../../" + relative, "../" + relative,
             "src/" + Path(relative).name)
    for mapping in maps:
        for name in names:
            value = mapping.get(name)
            if isinstance(value, str):
                return value, name
        # Some accepted reports store experiment files relative to their
        # archive directory, with the prefix omitted.
        for name, value in mapping.items():
            if isinstance(name, str) and Path(name).name == path.name and isinstance(value, str):
                return value, name
    return None, None


def _source_pin(root: Path, relative: str, maps: Iterable[Mapping[str, str]]) -> dict[str, Any]:
    path = _resolve(root, relative)
    if path is None:
        return {"path": relative, "exists": False, "pinned": False}
    actual = _digest(path)
    expected, pin_name = _expected_pin(relative, path, maps)
    return {"path": str(path.relative_to(root)).replace("\\", "/"),
            "sha256": actual, "expected_sha256": expected,
            "pin_name": pin_name, "exists": True,
            "pinned": expected == actual if expected is not None else False}


def _scope(tree: ast.AST, selector: str) -> list[ast.AST]:
    """Return exact class/method scopes, without importing the source."""
    if selector in {"<module>", "module", "RATES"}:
        return [tree]
    parts = selector.split(".")
    functions = [n for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
    if len(parts) == 1:
        return [n for n in functions if n.name == selector]
    if len(parts) != 2:
        return []
    cls, fn = parts
    out: list[ast.AST] = []
    for n in ast.walk(tree):
        if isinstance(n, ast.ClassDef) and n.name == cls:
            out.extend(x for x in n.body if isinstance(x, (ast.FunctionDef, ast.AsyncFunctionDef)) and x.name == fn)
    return out


def _ast_target(text: str) -> list[ast.AST]:
    try:
        body = ast.parse(text).body
    except SyntaxError:
        return []
    result: list[ast.AST] = []
    for node in body:
        # A call snippet is compared with the Call below an Expr wrapper.
        result.append(node.value if isinstance(node, ast.Expr) else node)
    return result


def _ast_dump(node: ast.AST) -> str:
    return ast.dump(node, annotate_fields=True, include_attributes=False)


def _statement(path: Path, scopes: list[ast.AST], text: str) -> dict[str, Any]:
    targets = _ast_target(text)
    wanted = {_ast_dump(x) for x in targets}
    found = any(_ast_dump(candidate) in wanted for scope in scopes for candidate in ast.walk(scope))
    target_hash = hashlib.sha256("\n".join(sorted(wanted)).encode()).hexdigest() if wanted else None
    return {"text": text, "ast_sha256": target_hash, "found": found}


def _function_evidence(path: Path, selector: str, statements: Iterable[str]) -> dict[str, Any]:
    try:
        tree = ast.parse(path.read_text(encoding="utf8"), filename=str(path))
    except (OSError, SyntaxError) as exc:
        return {"scope": selector, "passed": False, "error": repr(exc), "statements": []}
    scopes = _scope(tree, selector)
    rows = [_statement(path, scopes, text) for text in statements]
    return {"scope": selector, "passed": bool(scopes) and all(row["found"] for row in rows),
            "function_ast_sha256": [hashlib.sha256(_ast_dump(scope).encode()).hexdigest() for scope in scopes],
            "statements": rows}


def _text_has(path: Path, text: str) -> bool:
    try:
        return text in path.read_text(encoding="utf8")
    except OSError:
        return False


def _has_keyword(path: Path, selector: str, function_name: str, keyword: str) -> bool:
    try:
        tree = ast.parse(path.read_text(encoding="utf8"), filename=str(path))
    except (OSError, SyntaxError):
        return False
    for scope in _scope(tree, selector):
        for call in ast.walk(scope):
            if isinstance(call, ast.Call) and _call_name(call) == function_name:
                if any(k.arg == keyword for k in call.keywords):
                    return True
    return False


def _call_name(call: ast.Call) -> str:
    fn = call.func
    if isinstance(fn, ast.Name):
        return fn.id
    if isinstance(fn, ast.Attribute):
        return fn.attr
    return ""


def _for_target(path: Path, selector: str, names: tuple[str, ...]) -> bool:
    try:
        tree = ast.parse(path.read_text(encoding="utf8"), filename=str(path))
    except (OSError, SyntaxError):
        return False
    wanted = tuple(names)
    for scope in _scope(tree, selector):
        for row in ast.walk(scope):
            if not isinstance(row, ast.For) or not isinstance(row.target, ast.Tuple):
                continue
            got = tuple(x.id for x in row.target.elts if isinstance(x, ast.Name))
            if got == wanted:
                return True
    return False


def _keyword_evidence(path: Path, selector: str, keyword_names: Iterable[str]) -> dict[str, Any]:
    """Hash and find exact AST keyword nodes in a source call."""
    wanted = tuple(keyword_names)
    try:
        tree = ast.parse(path.read_text(encoding="utf8"), filename=str(path))
    except (OSError, SyntaxError) as exc:
        return {"scope": selector, "passed": False, "error": repr(exc), "keywords": []}
    scopes = _scope(tree, selector)
    rows = []
    for name in wanted:
        matches = []
        for scope in scopes:
            for call in ast.walk(scope):
                if isinstance(call, ast.Call):
                    for keyword in call.keywords:
                        if keyword.arg == name:
                            matches.append(hashlib.sha256(_ast_dump(keyword).encode()).hexdigest())
        rows.append({"name": name, "ast_sha256": matches[0] if matches else None,
                     "found": bool(matches)})
    return {"scope": selector, "passed": bool(scopes) and all(row["found"] for row in rows),
            "keywords": rows}


def _check(name: str, passed: bool, evidence: Any, *, method: str) -> dict[str, Any]:
    return {"passed": bool(passed), "method": method, "evidence": evidence}


def _binding(value: Any, key: str) -> Mapping[str, Any]:
    return _unwrap(value, key)


def _z_free(nodes: Any, index: Any, seen: set[int] | None = None) -> tuple[bool, list[int]]:
    """Check only that a graph handle has no ordinary Z dependency."""
    seen = set() if seen is None else seen
    i = _handle_index(index)
    if i is None or i in seen:
        return (i is not None, sorted(seen))
    seen.add(i)
    row = _node(nodes, i)
    if row is None:
        return False, sorted(seen)
    op = row.get("operation")
    if op == "bound_variable":
        return row.get("name") != "Z", sorted(seen)
    if op in {"exact_rational", "mathematical_pi", "current_original_source_parameter",
              "exact_positive_integer_expression"}:
        return bool(row.get("Z_independent", True)), sorted(seen)
    if op in {"sum", "product"}:
        results = [_z_free(nodes, x, seen)[0] for x in row.get("arguments", ())]
        return all(results), sorted(seen)
    if op == "negative":
        return _z_free(nodes, row.get("argument"), seen)[0], sorted(seen)
    if op in {"analytic_unary", "exact_original_source_flat_sigma", "fractional_part"}:
        return _z_free(nodes, row.get("argument"), seen)[0], sorted(seen)
    if op == "positive_quotient":
        return (_z_free(nodes, row.get("numerator"), seen)[0]
                and _z_free(nodes, row.get("denominator"), seen)[0]), sorted(seen)
    return False, sorted(seen)


def _unpack(args: tuple[Any, ...], owner: Any = None) -> tuple[Any, Mapping[str, Any], Any, Any, Any, Path, Mapping[str, str]]:
    if owner is None and len(args) == 1:
        owner = args[0]
    if owner is not None:
        nodes = getattr(owner, "nodes")
        loop = getattr(owner, "l")
        outer = getattr(owner, "outer")
        modified = getattr(owner, "modified")
        root = getattr(owner, "root", None)
        if root is None:
            root = getattr(owner, "directory", None)
            if root is not None:
                root = Path(root).parents[1]
        return nodes, loop, outer, getattr(owner, "loop", {}), modified, _repo_root(root), getattr(owner, "code_hashes", {})
    if len(args) != 6:
        raise TypeError("build_guards expects owner or nodes,l,outer,loop,modified,root")
    nodes, l, outer, loop, modified, root = args
    return nodes, l, outer, loop, modified, _repo_root(root), {}


def build_guards(*args: Any, owner: Any = None) -> dict[str, Any]:
    """Inspect accepted source premises and return compact, auditable guards.

    Supported calls are ``build_guards(owner)`` and
    ``build_guards(nodes, l, outer, loop, modified, root)``.  The latter is
    useful for a report reader that has already loaded the accepted archives.
    """
    nodes, l, outer, loop, modified, root, code_hashes = _unpack(args, owner)
    ob = _binding(outer, "prescribed_R100_Rm_C6_bindings")
    lb = _binding(loop, "prescribed_loop_C5_bindings")
    mb = _binding(modified, "modified_Q4_pressure_C3_bindings")
    maps = _pin_maps(code_hashes, (outer, loop, modified, ob, lb, mb))
    pins = {key: _source_pin(root, relative, maps) for key, relative in _SOURCE_FILES.items()}

    def req(key: str, selector: str, statements: Iterable[str]) -> dict[str, Any]:
        path = _resolve(root, _SOURCE_FILES[key])
        if path is None:
            return {"scope": selector, "passed": False, "source": pins[key], "statements": []}
        ev = _function_evidence(path, selector, statements)
        ev["source"] = pins[key]
        return ev

    def archive_flags(name: str, binding: Mapping[str, Any], required: Iterable[str]) -> dict[str, Any]:
        definition = binding.get("definition", {})
        rows = {key: definition.get(key) for key in required}
        def valid(key: str, value: Any) -> bool:
            if value is True:
                return True
            if key == "exact_Rref_Rz_Rrestore_Rm_offsets":
                return value == [-8, -7, -6]
            return False
        return _check(name, all(valid(key, value) for key, value in rows.items()), rows,
                      method="accepted_archive_definition_flags")

    checks: dict[str, Any] = {}
    checks["source_pins"] = _check("source_pins", all(v.get("pinned") for v in pins.values()), pins,
                                    method="SHA256_of_current_source_against_inherited_input_pin")
    checks["accepted_archives"] = _check(
        "accepted_archives",
        all(x["passed"] for x in (
            archive_flags("outer", ob, ("all_six_nonzero_incoming_and_every_chart_memory_retained",
                                         "exact_Rref_Rz_Rrestore_Rm_offsets",
                                         "original_axial_restoration_4Z_plus_actual_E_times_one_minus_sigma",
                                         "original_independent_P0_and_physical_Pstar_F0_factors_once")),
            archive_flags("loop", lb, ("one_original_common_N_eta_dstar_and_global_log_R_rminus_phase",
                                        "all_five_signed_finite_N_density_cross_terms_and_expm1_retained",
                                        "complete_original_rminus_to_Rm_Duhamel_correction_histories_C5")),
            archive_flags("modified", mb, ("actual_slow_native_loop_derivatives_from_true_function_bodies",
                                            "moving_integral_limits_and_implicit_inverse_chain_rule_retained",
                                            "total_modified_Ey_Vy_mean_y_C5_with_fast_phase_terms",
                                            "physical_U_mean_Uy_C5_and_Q_Qy_C4_Pstar_once")))),
        {"outer": ob.get("definition", {}), "loop": lb.get("definition", {}),
         "modified": mb.get("definition", {})}, method="accepted_archive_definition_flags")

    checks["native_derivative"] = _check(
        "native_derivative",
        all(x["passed"] for x in (
            req("modified", "NativeDerivative.__init__", ("self.cache={phase.node:g.zero}", "self.core_extensions=[]")),
            req("modified", "NativeDerivative.sub", ("self.g.node('function_substitution',expression=body.node,variable=variable.node,value=value.node,Z_independent_substitution=True,exact_native_chain_rule_substitution=True)",)),
            req("modified", "NativeDerivative.derivative", (
                "g.quotient(g.sub(g.mul(D(u.node),v),g.mul(u,D(v.node))),g.mul(v,v),row['source_positive_certificate'])",
                "g.mul(outer,dx)",
                "g.add(direct,chain)",
                "g.add(value,g.mul(self.sub(body,var,upper),self.derivative(upper)),g.neg(g.mul(self.sub(body,var,lower),self.derivative(lower))))",
                "g.quotient(g.sub(self.derivative(target),at(held)),lam,'original positive phase lambda; exact implicit native derivative')",
                "attrs['same_original_lazy_native_derivative']=True")),
            req("modified", "CurrentOriginalModifiedQ4PressureC3Functions.__init__", ("derivative.seed(ref,g.mul(jac,dr))",)),
        )),
        {"partial_native_chain_rule": req("modified", "NativeDerivative.derivative", (
            "g.add(direct,chain)",
            "g.add(value,g.mul(self.sub(body,var,upper),self.derivative(upper)),g.neg(g.mul(self.sub(body,var,lower),self.derivative(lower))))",
            "g.quotient(g.sub(self.derivative(target),at(held)),lam,'original positive phase lambda; exact implicit native derivative')",
        )), "phase_is_held": _text_has(_resolve(root, _SOURCE_FILES["modified"]) or Path("."), "Its total-y contribution is"),
         "seed_background_with_J_once": req("modified", "CurrentOriginalModifiedQ4PressureC3Functions.__init__", ("derivative.seed(ref,g.mul(jac,dr))",))},
        method="NativeDerivative_exact_source_chain_rules")

    checks["radial_and_pressure"] = _check(
        "radial_and_pressure",
        all(x["passed"] for x in (
            req("modified", "total_radial_functions", (
                "g.quotient(g.one,jac,'positive original chart log-radius Jacobian')",
                "g.quotient(g.one,N,'same original exact positive common N')",
                "a.add(a.scale(slow['A'],invN),chart['loop']['A_phi'])",
                "a.add(a.trim(bg['radial']['V_y']),a.scale(slow['B_over_Pstar'],invN),chart['loop']['B_over_Pstar_phi'])")),
            req("modified", "physical_pressure_fields", (
                "a.scale(fields['V'],Pstar)", "a.scale(radial['V_y'],Pstar)",
                "a.scale(fields['m'],Pstar)", "b.div(q,L,'same original L=1-delta*Z^2')",
                "Q=div(numerator(U,m))", "Qy=div(numerator(Uy,my))")),
            req("modified", "CurrentOriginalModifiedQ4PressureC3Functions.__init__", ("a.scale(Qy,chart['radius'])",)),
            req("bridge", "generic_fields", (
                "a.scale(F0,g.mul(g.unary('sqrt',g.mul(c(2),radius)),inverse))",
                "a.mul(W,actual['phi'])", "a.scale(actual['U'],inverse)",
                "a.scale(actual['U_y'],inverse)",
                "a.mul(W,a.add(a.scale(actual['phi'],'1/2'),actual['phi_y']))")),
        )),
        {"total_radial": req("modified", "total_radial_functions", (
            "g.quotient(g.one,jac,'positive original chart log-radius Jacobian')",
            "g.quotient(g.one,N,'same original exact positive common N')",
            "a.add(a.trim(bg['radial']['V_y']),a.scale(slow['B_over_Pstar'],invN),chart['loop']['B_over_Pstar_phi'])")),
         "physical_pressure": req("modified", "physical_pressure_fields", (
            "a.scale(fields['V'],Pstar)", "a.scale(radial['V_y'],Pstar)",
            "b.div(q,L,'same original L=1-delta*Z^2')", "Q=div(numerator(U,m))", "Qy=div(numerator(Uy,my))")),
         "physical_RQy": req("modified", "CurrentOriginalModifiedQ4PressureC3Functions.__init__", ("a.scale(Qy,chart['radius'])",))},
        method="source_physical_y_and_pressure_formula")

    checks["leading_germs"] = _check(
        "leading_germs",
        all(x["passed"] for x in (
            req("outer", "CurrentOriginalPrescribedR100RmC6Functions.__init__", (
                "g.mul(c(10),g.add(self.core.parameters['logC'],self.core.parameters['logP']))",
                "g.sub(g.sub(logref,T),c(8))", "a.scale(phi,'-2/5')",
                "a.add(a.scale(z,4),a.scale(E,g.sub(g.one,sigma(g,t))))",
                "a.scale(E,g.neg(sigma_prime(g,t)))", "a.scale(z,4)")),
            req("outer", "CurrentOriginalPrescribedR100RmC6Functions.phi_from_Utheta", (
                "a.div(Utheta,a.scale(self.F0,g.unary('sqrt',g.mul(g.constant(2),R))),'same positive sqrt(2R)*F0; F0 is Z-only')",)),
            req("bridge", "generic_fields", ("a.scale(a.mul(F0,F0),g.mul(radius,inverse,inverse))",)),
        )),
        {"reference_and_restoration_formulas": req("outer", "CurrentOriginalPrescribedR100RmC6Functions.__init__", (
            "g.mul(c(10),g.add(self.core.parameters['logC'],self.core.parameters['logP']))",
            "g.sub(g.sub(logref,T),c(8))", "a.scale(phi,'-2/5')",
            "a.add(a.scale(z,4),a.scale(E,g.sub(g.one,sigma(g,t))))", "a.scale(E,g.neg(sigma_prime(g,t)))")),
         "physical_Pstar_R_once": req("bridge", "generic_fields", ("a.scale(a.mul(F0,F0),g.mul(radius,inverse,inverse))",))},
        method="actual_prescribed_reference_restoration_source_formulas")

    checks["histories"] = _check(
        "histories",
        all(x["passed"] for x in (
            req("bridge_integrals", "<module>", ("RATES=dict(H=2,M=1,K=2,A=1,B=2,C=1)",)),
            req("bridge", "densities", ("a.scale(phi,2)", "a.mul(phi,U)", "a.mul(U,U)", "a.mul(phi,phi)")),
            req("bridge", "histories", ("a.scale(initial[name],decay)", "g.mul(jacobian,kernel)",)),
        )) and _for_target(_resolve(root, _SOURCE_FILES["bridge"]) or Path("."), "histories", ("name", "rate")),
        {"rates_and_densities": req("bridge_integrals", "<module>", ("RATES=dict(H=2,M=1,K=2,A=1,B=2,C=1)",)),
         "duhamel": req("bridge", "histories", ("a.scale(initial[name],decay)", "g.mul(jacobian,kernel)")),
         "loop_for_target_name_rate": _for_target(_resolve(root, _SOURCE_FILES["bridge"]) or Path("."), "histories", ("name", "rate")),
         "retention_markers": {key: _text_has(_resolve(root, _SOURCE_FILES["bridge"]) or Path("."), key) for key in
                              ("original_incoming_memory_retained=True", "source_measure_applied_once=True")}},
        method="actual_bridge_six_history_Duhamel_source")

    checks["generic_roots"] = _check(
        "generic_roots",
        all(x["passed"] for x in (
            req("loop", "generic_roots", (
                "a.sub(E,a.scale(a.trim(radial['E_y']),2))",
                "a.scale(a.trim(radial['V_y']),2)",
                "a.div(jet,E,'same strictly positive actual prescribed swirl E')",
                "aa,bb=divide(C),divide(B)")),
            req("loop", "loop_functions", (
                "a.div(gamma,a.scale(roots['a'],2),'same original active a>0')",
                "a.div(one,a.sqrt(a.add(one,a.mul(u,u))),'positive sqrt(1+u^2)')",
                "a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')")),
            req("loop", "MixedTaylor.div", ("g.quotient(g.sub(left.get((k,j),g.zero),rest),right[0,0],'1-2*r*cos(psi)+r^2=(cos(psi)-r)^2+sin(psi)^2>0 since |r|<1')",)),
        )),
        {"root_recipe": req("loop", "generic_roots", (
            "a.sub(E,a.scale(a.trim(radial['E_y']),2))", "a.scale(a.trim(radial['V_y']),2)",
            "a.div(jet,E,'same strictly positive actual prescribed swirl E')", "aa,bb=divide(C),divide(B)")),
         "universal_loop_positive_branches": req("loop", "loop_functions", (
            "a.div(gamma,a.scale(roots['a'],2),'same original active a>0')",
            "a.div(one,a.sqrt(a.add(one,a.mul(u,u))),'positive sqrt(1+u^2)')",
            "a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')")),
         "poisson_denominator": req("loop", "MixedTaylor.div", ("g.quotient(g.sub(left.get((k,j),g.zero),rest),right[0,0],'1-2*r*cos(psi)+r^2=(cos(psi)-r)^2+sin(psi)^2>0 since |r|<1')",))},
        method="same_universal_loop_source_recipe_C_equals_E_minus_2Ey_B_equals_2Vy")

    checks["flat_cutoff"] = _check(
        "flat_cutoff",
        all(x["passed"] for x in (
            _keyword_evidence(_resolve(root, _SOURCE_FILES["loop"]) or Path("."), "LoopAlgebra.lazy",
                              ("active_body_not_evaluated_on_flat_branch", "same_original_flat_extension_through_order5")),
            _keyword_evidence(_resolve(root, _SOURCE_FILES["loop"]) or Path("."), "LoopAlgebra.sigma",
                              ("original_sigma_logistic_ODE",)),
            _keyword_evidence(_resolve(root, _SOURCE_FILES["modified"]) or Path("."), "NativeDerivative.sigma_prime",
                              ("active_body_not_evaluated_on_flat_branch", "original_sigma_exact_logistic_first_derivative")),
            req("modified", "NativeDerivative.derivative", (
                "attrs['active_body']=D(row['active_body']).node",
                "attrs['flat_value']=D(row['flat_value']).node",
                "attrs['same_original_lazy_native_derivative']=True",)),
            _keyword_evidence(_resolve(root, _SOURCE_FILES["outer"]) or Path("."), "sigma_prime",
                              ("exact_flat_source_derivative_not_a_bound",)),
        )),
        {"lazy_extension": _keyword_evidence(_resolve(root, _SOURCE_FILES["loop"]) or Path("."), "LoopAlgebra.lazy",
                                              ("active_body_not_evaluated_on_flat_branch", "same_original_flat_extension_through_order5")),
         "sigma_logistic_ode": _keyword_evidence(_resolve(root, _SOURCE_FILES["loop"]) or Path("."), "LoopAlgebra.sigma",
                                                   ("original_sigma_logistic_ODE",)),
         "native_sigma_prime": _keyword_evidence(_resolve(root, _SOURCE_FILES["modified"]) or Path("."), "NativeDerivative.sigma_prime",
                                                   ("active_body_not_evaluated_on_flat_branch", "original_sigma_exact_logistic_first_derivative")),
         "native_flat_derivative": req("modified", "NativeDerivative.derivative", (
             "attrs['active_body']=D(row['active_body']).node",
             "attrs['flat_value']=D(row['flat_value']).node",
             "attrs['same_original_lazy_native_derivative']=True",)),
         "outer_sigma_prime": _keyword_evidence(_resolve(root, _SOURCE_FILES["outer"]) or Path("."), "sigma_prime",
                                                 ("exact_flat_source_derivative_not_a_bound",))},
        method="source_C_infinity_flat_branch_and_analytic_active_branch")

    # The source contains the fundamental inverse and explicitly freezes the
    # slow Z phase.  It does not state the global periodic extension identity,
    # so this record remains false and no B_y alias is authorized here.
    loop_path = _resolve(root, _SOURCE_FILES["loop"])
    control_path = _resolve(root, _SOURCE_FILES["control"])
    phase_handles = {name: {key: lb.get("charts", {}).get(name, {}).get(key)
                            for key in ("coordinate", "global_offset", "actual_phase", "radius", "physical_measure")}
                    for name in ("reference", "restoration", "postrestore")}
    phase_rows = []
    for name, row in phase_handles.items():
        actual = _node(nodes, row["actual_phase"])
        phase_rows.append({"chart": name, "total_Z_phase_derivative_exact_zero":
                           lb.get("charts", {}).get(name, {}).get("total_Z_phase_derivative_exact_zero") is True,
                           "actual_phase_fractional_part": bool(actual and actual.get("operation") == "analytic_unary"
                                                                and actual.get("name") == "fractional_part"),
                           "radius_Z_free": _z_free(nodes, row["radius"])[0],
                           "physical_measure_Z_free": _z_free(nodes, row["physical_measure"])[0]})
    phase_local = bool(loop_path and control_path and all(row["total_Z_phase_derivative_exact_zero"] and
                                                          row["actual_phase_fractional_part"] and
                                                          row["radius_Z_free"] and row["physical_measure_Z_free"]
                                                          for row in phase_rows))
    phase_ast = []
    if loop_path:
        phase_ast.append(_keyword_evidence(loop_path, "loop_functions",
                                           ("phase_endpoints", "flat_inverse", "branch_inverse_is_original_not_support_cap")))
    if control_path:
        phase_ast.append(_keyword_evidence(control_path, "original_loop",
                                           ("phase_endpoints", "flat_inverse", "branch_inverse_is_original_not_support_cap")))
    checks["phase_global_R_Z0"] = _check("phase_global_R_Z0", phase_local and all(x["passed"] for x in phase_ast),
                                          {"handles": phase_handles, "rows": phase_rows, "source": phase_ast},
                                          method="fractional_part_N_offset_with_exact_Z_independent_R_and_measure")
    checks["global_periodic_lift"] = _check(
        "global_periodic_lift", False,
        {"fundamental_inverse_source": phase_ast,
         "reason": "No exact source AST statement of Phi(theta+2*pi)=Phi(theta)+1; external periodic C1 theorem required."},
        method="explicit_gap_no_unproved_phase_alias")

    # Handles and exact graph leaves are intentionally small: this is source
    # identity evidence, not a graph dump.
    selected_handles = {
        "N": _handle_index(lb.get("N")), "eta": _handle_index(lb.get("eta")),
        "d_star": _handle_index(lb.get("dstar")), "sc": _handle_index(lb.get("sc")),
        "Z": 20, "delta": 11, "Pstar": _handle_index(ob.get("Pstar")),
        "reference_radius": phase_handles["reference"]["radius"],
        "reference_physical_measure": phase_handles["reference"]["physical_measure"],
    }
    selected_rows = {str(i): _node(nodes, i) for i in (10, 11, 18, 20, selected_handles["Pstar"])
                     if i is not None and _node(nodes, i) is not None}
    nrow = _node(nodes, selected_handles["N"])
    zrow = _node(nodes, selected_handles["Z"])
    drow = _node(nodes, selected_handles["delta"])
    prow = _node(nodes, selected_handles["Pstar"])
    eta_row = _node(nodes, selected_handles["eta"])
    ds_row = _node(nodes, selected_handles["d_star"])
    sc_row = _node(nodes, selected_handles["sc"])
    parameter_pass = bool(
        nrow and nrow.get("operation") == "exact_positive_integer_expression" and
        nrow.get("exact_positive_integer") is True and nrow.get("Z_independent") is True and
        nrow.get("one_common_N_for_all_charts_and_control_denominators") is True and
        zrow and zrow.get("operation") == "bound_variable" and zrow.get("name") == "Z" and
        drow and drow.get("operation") == "exact_rational" and drow.get("numerator") == 1 and drow.get("denominator") == 1000 and
        prow and prow.get("operation") == "analytic_unary" and prow.get("name") == "exp" and prow.get("argument") == 5 and
        all(row and row.get("operation") == "current_original_source_parameter" and row.get("Z_independent") is True
            for row in (eta_row, ds_row, sc_row)) and
        isinstance(eta_row.get("canonical_parameter_proof") if eta_row else None, Mapping) and
        isinstance(ds_row.get("canonical_parameter_proof") if ds_row else None, Mapping))
    checks["parameters"] = _check(
        "parameters", parameter_pass,
        {"handles": selected_handles, "nodes": selected_rows,
         "canonical_log_sources": {key: (row or {}).get("source_attribute") for key, row in
                                    (("eta", eta_row), ("d_star", ds_row), ("sc", sc_row))},
         "Pstar_is_exp_of_logP_node_5": bool(prow and prow.get("argument") == 5),
         "delta_exact_1_over_1000": bool(drow and drow.get("numerator") == 1 and drow.get("denominator") == 1000)},
        method="source_graph_exact_N_delta_Pstar_and_canonical_parameter_handles")

    checks["phase_source_recipe"] = _check(
        "phase_source_recipe",
        bool(control_path and _text_has(control_path, "phase_Z_exact_zero=True") and
             _text_has(control_path, "Full-phase original inverse ranges contain the exact fractional_part(N*offset)")),
        {"source": pins["control"], "phase_Z_exact_zero": True,
         "full_phase_range_statement": True},
        method="accepted_control_source_full_phase_statement")

    ode = req("modified", "background_radial_rows", (
        "a.sub(V,f['m'])", "a.sub(E,a.scale(f['h'],'3/2'))",
        "a.sub(a.mul(E,V),a.scale(f['k'],'3/2'))",
        "a.sub(a.sub(a.mul(V,V),a.scale(E2,'1/2')),f['e'])",
        "a.scale(E2,'1/2')", "result['pressure']=result['p']"))
    axial = req("inlet", "generic_functions", (
        "a.add(a.mul(a.sub(one,de),a.mul(z,m)),a.mul(d,mz))",
        "div(a.add(a.neg(V),half(a.mul(a.sub(one,de),a.sub(m,a.mul(z,mz))))))",
        "div(a.add(a.mul(V,transport),a.scale(a.mul(de,a.mul(z,e)),2),a.neg(a.mul(d,ez)),a.scale(a.mul(a.add(one,de),a.mul(z,pressure)),2),a.neg(a.mul(d,pz))))",
        "a.scale(a.div(a.add(axial_linear,a.scale(axial_quadratic,Pstar)),E,'same actual positive core inlet swirl E'),radius)"))
    checks["normalized_ODE_and_axial_recipe"] = _check(
        "normalized_ODE_and_axial_recipe", ode["passed"] and axial["passed"],
        {"physical_moment_ODE": ode, "root_p2_source_formula": axial},
        method="actual_normalized_history_ODE_and_original_axial_generic_source")

    # ``global_periodic_lift`` and the archive numeric flags are deliberate
    # limits.  The source premises can be passed while the theorem remains
    # conditional on the external periodic C1 proof.
    source_premise_names = tuple(key for key in checks if key not in {"global_periodic_lift"})
    passed = all(checks[key].get("passed", False) for key in source_premise_names)
    return {
        "passed": passed,
        "method": "pinned_source_germ_premises_only",
        "source_only": True,
        "numeric_fixtures_used": False,
        "checks": checks,
        "source_pins": pins,
        "parameter_handles": selected_handles,
        "actual_parameter_handles": selected_handles,
        "global_periodic_lift_proved": False,
        "limits": {
            "full_inlet_theorem": False,
            "B_y_alias_authorized": False,
            "manual_Uy_aliases": False,
            "global_periodic_lift_requires_external_C1_theorem": True,
            "numeric_fixture_or_receipt_used_as_source": False,
        },
    }


__all__ = ["build_guards"]
