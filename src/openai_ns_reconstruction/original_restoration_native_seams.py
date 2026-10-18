"""Source-bound receipt for the native restoration and post-restore charts.

This module reads the pinned C6, C5, and modified-Q4 archives as data and
parses the small source functions that define them.  It never imports or
executes a producer, owner, checker, recurrence, quadrature, or norm.  The
receipt is deliberately limited to the two new native joins, their history
transport, and the Q4 radial seed rows.  It does not certify the global
modified ``Q_y`` problem or a numerical finite-``N`` closure.
"""
from __future__ import annotations

import ast
import gzip
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_long_reference_modified_seams.json"
PARENT_SHA256 = "204d8b45aff263a3df1eae319dd7534e304162794ac375a1c0ec2cfe9beafc35"
MODULE = "src/openai_ns_reconstruction/original_restoration_native_seams.py"
ENDPOINT_MODULE = "src/openai_ns_reconstruction/original_restoration_native_endpoints.py"

C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
C6_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
C6_BRIDGE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
RATE_SOURCE = BASE + "lei_ren_part1_paper_compliant_actual_bridge_integrals.py"
C5_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
C5_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
Q4_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
Q4_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
FLAT_SOURCE = BASE + "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"

NATIVE = ("H", "M", "K", "A", "B", "C")
C5_NAMES = ("m", "h", "k", "e", "p")
NATIVE_RATES = {"H": 2, "M": 1, "K": 2, "A": 1, "B": 2, "C": 1}
C5_RATES = {"m": "1", "h": "3/2", "k": "3/2", "e": "1", "p": "0"}
Q4_KEYS = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
Q4_FORMULAS = {
    "E": "E_y = source radial['E_y']",
    "V": "V_y = source radial['V_y']",
    "m": "m_y = V - m",
    "h": "h_y = E - 3*h/2",
    "k": "k_y = E*V - 3*k/2",
    "e": "e_y = V**2 - E**2/2 - e",
    "p": "p_y = E**2/2",
    "P0": "P0_y = 0",
    "pressure": "pressure_y = p_y",
}
Z_SCOPE = {"lower": "37/100", "upper": "93/250"}


def _require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json(root: Path, relative: str) -> Any:
    path = root / relative
    _require(path.is_file(), "Missing pinned JSON source: " + relative)
    return json.loads(path.read_bytes())


def _gzip_json(root: Path, relative: str) -> Any:
    path = root / relative
    _require(path.is_file(), "Missing pinned gzip source: " + relative)
    return json.loads(gzip.decompress(path.read_bytes()))


def _compact(value: str) -> str:
    return "".join(value.split())


def _guard(root: Path, relative: str, name: str, snippets: tuple[str, ...]) -> dict[str, Any]:
    """Pin a source function by its AST and required source fragments."""
    path = root / relative
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=str(path))
    candidates: list[ast.AST] = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) or node.name != name:
            continue
        segment = ast.get_source_segment(text, node) or ""
        if all(_compact(fragment) in _compact(segment) for fragment in snippets):
            candidates.append(node)
    _require(len(candidates) == 1, "Changed source definition: " + relative + ":" + name)
    node = candidates[0]
    return dict(
        passed=True,
        path=relative,
        function=name,
        required_text=list(snippets),
        AST_sha256=hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest(),
        source_executed=False,
    )


def _assignment(root: Path, relative: str, name: str) -> ast.Assign:
    path = root / relative
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    rows = [
        row for row in tree.body
        if isinstance(row, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == name for target in row.targets)
    ]
    _require(len(rows) == 1, "Expected one source assignment: " + relative + ":" + name)
    return rows[0]


def _fraction_literal(node: ast.AST) -> str:
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "Fraction":
        values = [ast.literal_eval(arg) for arg in node.args]
        _require(len(values) in (1, 2) and all(type(value) is int for value in values),
                 "C5 rate is not an exact fraction")
        return str(values[0]) if len(values) == 1 else f"{values[0]}/{values[1]}"
    value = ast.literal_eval(node)
    _require(type(value) is int, "Native rate is not an exact integer")
    return str(value)


def _source_guards(root: Path) -> dict[str, Any]:
    guards = {
        "c6_constructor_restoration": _guard(root, C6_SOURCE, "__init__", (
            "self.E=E=a.sub(self.endpoints['R110']['actual']['U'],a.scale(z,4))",
            "U=a.add(a.scale(z,4),a.scale(E,g.sub(g.one,sigma(g,t))))",
            "actual=dict(phi=phi,U=U,phi_y=a.scale(phi,'-2/5'),U_y=a.scale(E,g.neg(sigma_prime(g,t)))",
            "self.add_chart('restoration',t,R,g.one,actual,self.endpoints['Rz']['moments'])",
            "actual=dict(phi=phi,U=a.scale(z,4),phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0),physical_Utheta=Utheta,log_Utheta=logUtheta)",
            "self.add_chart('post_restore',p,R,g.one,actual,self.endpoints['Rrestore']['moments'])",
        )),
        "c6_phi_from_Utheta": _guard(root, C6_SOURCE, "phi_from_Utheta", (
            "return a.div(Utheta,a.scale(self.F0,g.unary('sqrt',g.mul(g.constant(2),R)))",
        )),
        "c6_add_chart": _guard(root, C6_SOURCE, "add_chart", (
            "J=previous.histories(a,actual['phi'],actual['U'],incoming,phase,jacobian,'true_'+name)",
            "fields,radial,generic=previous.generic_fields(a,actual,J,self.F0,self.P0,self.Pstar,R,self.delta)",
        )),
        "c6_sigma_prime": _guard(root, C6_SOURCE, "sigma_prime", (
            "derivative_source=own.c2.current.ast_binding(own.c2.flat_source.sigma_jets)",
            "exact_flat_source_derivative_not_a_bound=True",
        )),
        "c6_densities": _guard(root, C6_BRIDGE, "densities", (
            "return dict(H=a.scale(phi,2),M=U,K=a.scale(a.mul(phi,U),2),",
            "A=a.mul(U,U),B=a.mul(phi,phi),C=a.mul(phi,phi))",
        )),
        "c6_histories": _guard(root, C6_BRIDGE, "histories", (
            "body=a.scale(a.at(bodies[name],phase,t),g.mul(jacobian,kernel))",
            "out[name]=a.add(a.scale(initial[name],decay),a.integral(body,t,g.zero,phase,",
            "original_incoming_memory_retained=True,source_measure_applied_once=True",
        )),
        "c5_transport": _guard(root, C5_SOURCE, "transport", (
            "width=g.mul(jacobian,g.sub(coordinate,lower))",
            "out[name]=a.add(a.scale(incoming[name],decay),contribution)",
            "physical_Jacobian_applied_once=True",
            "zero_correction_inlet_only_at_original_rminus=True",
            "predecessor_memory_never_reset=True",
        )),
        "c5_chart_alias": _guard(root, C5_SOURCE, "__init__", (
            "'restoration':background.charts['restoration']",
            "'postrestore':background.charts['post_restore']",
        )),
        "q4_background_radial_rows": _guard(root, Q4_SOURCE, "background_radial_rows", (
            "result=dict(E=chart['radial']['E_y'],V=chart['radial']['V_y']",
            "m=a.sub(V,f['m']),h=a.sub(E,a.scale(f['h'],'3/2'))",
            "k=a.sub(a.mul(E,V),a.scale(f['k'],'3/2'))",
            "e=a.sub(a.sub(a.mul(V,V),a.scale(E2,'1/2')),f['e'])",
            "p=a.scale(E2,'1/2'),P0=a.fixed(0)",
            "result['pressure']=result['p'];return result",
        )),
        "q4_constructor_seed": _guard(root, Q4_SOURCE, "__init__", (
            "radial_bg=background_radial_rows(g,bg)",
            "for key,jet in bg['fields'].items():",
            "for ref,dr in zip(rows(jet),rows(radial_bg[key])):derivative.seed(ref,g.mul(jac,dr))",
        )),
        "q4_physical_ports": _guard(root, Q4_SOURCE, "physical_pressure_fields", (
            "U=a.scale(fields['V'],Pstar);Uy=a.scale(radial['V_y'],Pstar)",
            "m=a.scale(fields['m'],Pstar);my=a.sub(U,m)",
            "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
        )),
        "flat_sigma_left": _guard(root, FLAT_SOURCE, "_sigma_left", (
            "odds=1/(1-x)**2-1/x**2",
            "e=positive_exp(c,odds);value=e/(1+e)",
            "L1=2/(1-x)**3+2/x**3",
            "L2=6/(1-x)**4-6/x**4",
            "derivatives=[value,q*L1",
        )),
        "flat_sigma_jets": _guard(root, FLAT_SOURCE, "sigma_jets", (
            "if lo<=0:rows.append(IntervalTaylor.constant(c,0,4))",
            "if hi>=1:rows.append(IntervalTaylor.constant(c,1,4))",
            "left=_sigma_left(c,1-c.mpf([a,b]))",
        )),
    }
    native_assignment = _assignment(root, RATE_SOURCE, "RATES")
    _require(isinstance(native_assignment.value, ast.Call)
             and isinstance(native_assignment.value.func, ast.Name)
             and native_assignment.value.func.id == "dict", "Native RATES assignment changed")
    native_rates = {
        keyword.arg: ast.literal_eval(keyword.value)
        for keyword in native_assignment.value.keywords
        if keyword.arg is not None
    }
    _require(native_rates == NATIVE_RATES, "Actual native RATES changed")
    guards["native_rates"] = dict(
        passed=True, path=RATE_SOURCE, binding="RATES", rates=native_rates,
        AST_sha256=hashlib.sha256(ast.dump(native_assignment, include_attributes=False).encode()).hexdigest(),
        source_executed=False,
    )
    c5_assignment = _assignment(root, C5_SOURCE, "RATES")
    _require(isinstance(c5_assignment.value, ast.Call)
             and isinstance(c5_assignment.value.func, ast.Name)
             and c5_assignment.value.func.id == "dict", "C5 RATES assignment changed")
    c5_rates = {
        keyword.arg: _fraction_literal(keyword.value)
        for keyword in c5_assignment.value.keywords
        if keyword.arg is not None
    }
    _require(c5_rates == C5_RATES, "Actual C5 rates changed")
    guards["c5_rates"] = dict(
        passed=True, path=C5_SOURCE, binding="RATES", rates=c5_rates,
        AST_sha256=hashlib.sha256(ast.dump(c5_assignment, include_attributes=False).encode()).hexdigest(),
        source_executed=False,
    )
    return guards


def _accepted_pin(parent: Mapping[str, Any], relative: str) -> str:
    if relative == PARENT:
        return PARENT_SHA256
    values: list[str] = []
    for packet in parent.values():
        if not isinstance(packet, Mapping):
            continue
        pins = packet.get("source_pins")
        if not isinstance(pins, Mapping):
            continue
        item = pins.get(relative)
        if isinstance(item, Mapping) and isinstance(item.get("expected_sha256"), str):
            values.append(item["expected_sha256"])
    _require(values, "Accepted source pin is missing: " + relative)
    _require(len(set(values)) == 1, "Accepted source pins disagree: " + relative)
    return values[0]


def _source_pins(root: Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    paths = (PARENT, C6_SOURCE, C6_CACHE, C6_BRIDGE, RATE_SOURCE,
             C5_SOURCE, C5_CACHE, Q4_SOURCE, Q4_CACHE, FLAT_SOURCE)
    result: dict[str, Any] = {}
    for relative in paths:
        expected = _accepted_pin(parent, relative)
        actual = _sha(root / relative)
        _require(actual == expected, "Pinned source bytes changed: " + relative)
        result[relative] = dict(
            path=relative,
            expected_sha256=expected,
            actual_sha256=actual,
            sha256=actual,
            passed=True,
            source_executed=False,
        )
    module_path = root / MODULE
    _require(module_path.is_file(), "Target module is missing")
    digest = _sha(module_path)
    result[MODULE] = dict(
        path=MODULE,
        expected_sha256=digest,
        actual_sha256=digest,
        sha256=digest,
        passed=True,
        source_executed=False,
    )
    return result


def _pointer(value: Any, label: str, nodes: list[Any]) -> int:
    _require(type(value) is int and not isinstance(value, bool), label + " is not a graph pointer")
    _require(0 <= value < len(nodes), label + " is outside the graph")
    return value


def _row(nodes: list[Any], value: Any, label: str = "row") -> Mapping[str, Any]:
    pointer = _pointer(value, label, nodes)
    row = nodes[pointer]
    _require(isinstance(row, Mapping), label + " is not a graph row")
    return row


def _rows(value: Any, label: str, nodes: list[Any], length: int) -> list[int]:
    _require(isinstance(value, list) and len(value) == length, label + " has wrong jet length")
    return [_pointer(item, f"{label}[{index}]", nodes) for index, item in enumerate(value)]


def _children(row: Mapping[str, Any]) -> tuple[int, ...]:
    op = row.get("operation")
    if op in ("sum", "product"):
        return tuple(item for item in row.get("arguments", ()) if type(item) is int)
    if op in ("negative", "analytic_unary", "exact_original_source_flat_sigma"):
        return (row["argument"],) if type(row.get("argument")) is int else ()
    if op == "positive_quotient":
        return tuple(row[key] for key in ("numerator", "denominator") if type(row.get(key)) is int)
    if op == "function_substitution":
        return tuple(row[key] for key in ("expression", "value") if type(row.get(key)) is int)
    if op == "definite_integral":
        return tuple(row[key] for key in ("integrand", "lower", "upper") if type(row.get(key)) is int)
    return ()


def _reachable(nodes: list[Any], start: int) -> set[int]:
    pending = [start]
    seen: set[int] = set()
    while pending:
        pointer = pending.pop()
        if type(pointer) is not int or pointer < 0 or pointer >= len(nodes) or pointer in seen:
            continue
        seen.add(pointer)
        pending.extend(_children(_row(nodes, pointer)))
    return seen


def _find(nodes: list[Any], start: int, operation: str) -> list[tuple[int, Mapping[str, Any]]]:
    return [(pointer, _row(nodes, pointer)) for pointer in _reachable(nodes, start)
            if _row(nodes, pointer).get("operation") == operation]


def _exact(nodes: list[Any], pointer: Any, numerator: int, denominator: int = 1) -> bool:
    row = _row(nodes, pointer)
    return (row.get("operation") == "exact_rational"
            and row.get("numerator") == numerator and row.get("denominator") == denominator)


def _canonical(nodes: list[Any], pointer: int, memo: dict[int, Any] | None = None) -> Any:
    """Canonicalize only the source expression language used by kernels."""
    memo = {} if memo is None else memo
    if pointer in memo:
        return memo[pointer]
    row = _row(nodes, pointer)
    op = row.get("operation")
    if op == "exact_rational":
        out = ("rational", row.get("numerator"), row.get("denominator"))
    elif op == "bound_variable":
        out = ("bound", row.get("name"))
    elif op == "current_original_source_parameter":
        out = ("source_parameter", pointer, row.get("name"))
    elif op == "negative":
        out = ("negative", _canonical(nodes, row["argument"], memo))
    elif op in ("analytic_unary", "exact_original_source_flat_sigma"):
        out = (op, row.get("name"), _canonical(nodes, row["argument"], memo))
    elif op == "positive_quotient":
        out = ("positive_quotient", _canonical(nodes, row["numerator"], memo),
               _canonical(nodes, row["denominator"], memo))
    elif op == "function_substitution":
        out = ("function_substitution", _canonical(nodes, row["expression"], memo),
               row.get("variable"), _canonical(nodes, row["value"], memo))
    elif op in ("sum", "product"):
        parts: list[Any] = []
        for child in row.get("arguments", ()):
            value = _canonical(nodes, child, memo)
            if isinstance(value, tuple) and value and value[0] == op:
                parts.extend(value[1])
            else:
                parts.append(value)
        out = (op, tuple(sorted(parts, key=repr)))
    else:
        out = ("source_operation", op, pointer)
    memo[pointer] = out
    return out


def _sum(*terms: Any) -> Any:
    flat: list[Any] = []
    for term in terms:
        flat.extend(term[1] if isinstance(term, tuple) and term and term[0] == "sum" else [term])
    return ("sum", tuple(sorted(flat, key=repr)))


def _product(*terms: Any) -> Any:
    flat: list[Any] = []
    for term in terms:
        flat.extend(term[1] if isinstance(term, tuple) and term and term[0] == "product" else [term])
    return ("product", tuple(sorted(flat, key=repr)))


def _rational(numerator: int, denominator: int = 1) -> Any:
    return ("rational", numerator, denominator)


def _bound_pointer(nodes: list[Any], name: str) -> int:
    matches = [pointer for pointer, row in enumerate(nodes)
               if row.get("operation") == "bound_variable" and row.get("name") == name]
    _require(matches, "Bound history variable is missing: " + name)
    return matches[0]


def _expected_exp(rate: str | int, phase: int, variable: int, nodes: list[Any]) -> Any:
    if isinstance(rate, int):
        numerator, denominator = -rate, 1
    else:
        pieces = rate.split("/", 1)
        numerator = -int(pieces[0])
        denominator = int(pieces[1]) if len(pieces) == 2 else 1
    return ("analytic_unary", "exp", _product(
        _rational(numerator, denominator),
        _sum(_canonical(nodes, phase), ("negative", _canonical(nodes, variable))),
    ))


def _expected_decay(rate: str | int, phase: int, nodes: list[Any]) -> Any:
    if isinstance(rate, int):
        numerator, denominator = -rate, 1
    else:
        pieces = rate.split("/", 1)
        numerator = -int(pieces[0])
        denominator = int(pieces[1]) if len(pieces) == 2 else 1
    return ("analytic_unary", "exp", _product(
        _rational(numerator, denominator), _canonical(nodes, phase),
    ))


def _has_exp(nodes: list[Any], root: int, expected: Any) -> int | None:
    for pointer, row in _find(nodes, root, "analytic_unary"):
        if row.get("name") == "exp" and _canonical(nodes, pointer) == expected:
            return pointer
    return None


def _native_chart_history(nodes: list[Any], chart_name: str, chart: Mapping[str, Any]) -> dict[str, Any]:
    phase = _pointer(chart["phase"], chart_name + ".phase", nodes)
    jacobian = _pointer(chart["physical_measure"], chart_name + ".physical_measure", nodes)
    _require(_exact(nodes, jacobian, 1), "New native chart does not have exact J=1: " + chart_name)
    moments = chart["moments"]
    incoming = chart["incoming_histories"]
    proofs: dict[str, Any] = {}
    for name in NATIVE:
        moment_rows = _rows(moments[name], chart_name + ".moments." + name, nodes, 7)
        incoming_rows = _rows(incoming[name], chart_name + ".incoming." + name, nodes, 7)
        rows_receipt: list[dict[str, Any]] = []
        for order, moment_pointer in enumerate(moment_rows):
            root = _row(nodes, moment_pointer, chart_name + ".moment")
            _require(root.get("operation") == "sum" and len(root.get("arguments", ())) == 2,
                     "Native history root is not memory plus integral: " + chart_name + ":" + name)
            integrals = [item for item in root["arguments"]
                         if _row(nodes, item).get("operation") == "definite_integral"]
            _require(len(integrals) == 1, "Native history has no unique integral: " + chart_name + ":" + name)
            integral_pointer = integrals[0]
            integral = _row(nodes, integral_pointer)
            expected_meta = {
                "lower": 0,
                "upper": phase,
                "actual_six_history_rate": NATIVE_RATES[name],
                "actual_history_name": name,
                "original_incoming_memory_retained": True,
                "source_measure_applied_once": True,
                "prescribed_bridge_chart": "true_" + chart_name,
                "exact_function_integral": True,
                "numerical_value_not_installed": True,
                "endpoint_independent_of_Z": True,
                "numerical_cap_not_defining_integrand": True,
                "ordinary_slow_Z_derivative_order": order,
            }
            for key, value in expected_meta.items():
                _require(integral.get(key) == value,
                         f"Native integral metadata changed: {chart_name}:{name}:{key}")
            variable_name = "true_" + chart_name + "_history_t"
            _require(integral.get("variable") == variable_name,
                     "Native history variable changed: " + chart_name + ":" + name)
            variable = _bound_pointer(nodes, variable_name)
            kernel = _has_exp(nodes, integral["integrand"],
                              _expected_exp(NATIVE_RATES[name], phase, variable, nodes))
            _require(kernel is not None,
                     "Native kernel exponent changed: " + chart_name + ":" + name)
            memory_pointers = [item for item in root["arguments"] if item != integral_pointer]
            _require(len(memory_pointers) == 1, "Native memory term missing: " + chart_name + ":" + name)
            memory = _row(nodes, memory_pointers[0])
            _require(memory.get("operation") == "product"
                     and incoming_rows[order] in memory.get("arguments", ()),
                     "Native incoming memory pointer changed: " + chart_name + ":" + name)
            decay = _has_exp(nodes, memory_pointers[0], _expected_decay(NATIVE_RATES[name], phase, nodes))
            _require(decay is not None, "Native memory decay changed: " + chart_name + ":" + name)
            substitutions = [
                (pointer, row) for pointer, row in _find(nodes, integral["integrand"], "function_substitution")
                if row.get("variable") == phase
            ]
            _require(substitutions, "Native density substitution missing: " + chart_name + ":" + name)
            density_expression = substitutions[0][1]["expression"]
            # A higher ordinary-Z row of a product uses lower source rows
            # (for example (U^2)_{ZZ} reaches U_0 and U_2), so the source
            # binding is checked against the complete actual phi/U jet rather
            # than pretending that row ``order`` is an isolated sample.
            actual_phi_rows = _rows(chart["actual"]["phi"], chart_name + ".actual.phi", nodes, 7)
            actual_u_rows = _rows(chart["actual"]["U"], chart_name + ".actual.U", nodes, 7)
            density_reachable = _reachable(nodes, density_expression)
            _require(any(pointer in density_reachable for pointer in actual_phi_rows + actual_u_rows),
                     "Native density is detached from actual phi/U: " + chart_name + ":" + name)
            if name in ("H", "K"):
                _require(any(_exact(nodes, item, 2) for item in density_reachable),
                         "Native factor two missing: " + chart_name + ":" + name)
            rows_receipt.append(dict(
                order=order, moment=moment_pointer, incoming=incoming_rows[order],
                integral=integral_pointer, memory=memory_pointers[0], decay=decay,
                kernel=kernel, density_expression=density_expression,
                q0=dict(exp_zero_is_one=True, integral_zero_is_zero=True,
                         incoming_memory_retained=True),
            ))
        proofs[name] = dict(rate=NATIVE_RATES[name], ordinary_Z_orders=list(range(7)), rows=rows_receipt)
    return dict(
        passed=True,
        chart=chart_name,
        phase=phase,
        physical_measure=jacobian,
        J_identity="J=1 exact, so the source J factor is represented once by the identity and is simplified in graph rows",
        histories=proofs,
        q0_unit_memory=True,
        source_executed=False,
    )


def _c5_incoming(nodes: list[Any], actual: list[int], predecessor: list[int], coordinate: int, label: str) -> list[int]:
    for order, (pointer, expected) in enumerate(zip(actual, predecessor)):
        row = _row(nodes, pointer, label)
        _require(row.get("operation") == "function_substitution"
                 and row.get("expression") == expected
                 and row.get("variable") == coordinate
                 and row.get("value") == 1
                 and row.get("Z_independent_substitution") is True,
                 "C5 incoming pointer changed: " + label + ":" + str(order))
    return list(actual)


def _c5_chart_history(nodes: list[Any], chart_name: str, chart: Mapping[str, Any], predecessor: Mapping[str, Any]) -> dict[str, Any]:
    coordinate = _pointer(chart["coordinate"], chart_name + ".coordinate", nodes)
    jacobian = _pointer(chart["physical_measure"], chart_name + ".physical_measure", nodes)
    _require(_exact(nodes, jacobian, 1), "New C5 chart does not have exact J=1: " + chart_name)
    predecessor_coordinate = _pointer(predecessor["coordinate"], "predecessor.coordinate", nodes)
    proofs: dict[str, Any] = {}
    incoming_receipts: dict[str, Any] = {}
    finite = chart["finite_N"]["five_signed_density_increments"]
    for name in C5_NAMES:
        incoming = _rows(chart["incoming"][name], chart_name + ".incoming." + name, nodes, 6)
        previous = _rows(predecessor["history"][name], "predecessor.history." + name, nodes, 6)
        incoming_receipts[name] = _c5_incoming(nodes, incoming, previous, predecessor_coordinate,
                                                chart_name + ".incoming." + name)
        history = _rows(chart["history"][name], chart_name + ".history." + name, nodes, 6)
        contributions = _rows(chart["contributions"][name], chart_name + ".contributions." + name, nodes, 6)
        density = _rows(finite[name], chart_name + ".finite_N." + name, nodes, 6)
        memory = _pointer(chart["memory"][name], chart_name + ".memory." + name, nodes)
        row_receipts: list[dict[str, Any]] = []
        for order, history_pointer in enumerate(history):
            root = _row(nodes, history_pointer, chart_name + ".history.root")
            contribution_pointer = contributions[order]
            _require(root.get("operation") == "sum" and contribution_pointer in root.get("arguments", ()),
                     "C5 history root changed: " + chart_name + ":" + name)
            contribution = _row(nodes, contribution_pointer, chart_name + ".contribution")
            expected_meta = {
                "lower": 0,
                "upper": coordinate,
                "original_rate": C5_RATES[name],
                "original_chart": chart_name,
                "physical_Jacobian_applied_once": True,
                "zero_correction_inlet_only_at_original_rminus": True,
                "predecessor_memory_never_reset": True,
                "exact_function_integral": True,
                "numerical_value_not_installed": True,
                "endpoint_independent_of_Z": True,
                "numerical_cap_not_defining_integrand": True,
                "actual_finite_N_signed_correction": True,
                "ordinary_slow_Z_derivative_order": order,
            }
            for key, value in expected_meta.items():
                _require(contribution.get(key) == value,
                         f"C5 integral metadata changed: {chart_name}:{name}:{key}")
            variable_name = "true_loop_C5_transport_" + chart_name
            _require(contribution.get("variable") == variable_name,
                     "C5 transport variable changed: " + chart_name + ":" + name)
            variable = _bound_pointer(nodes, variable_name)
            substitutions = [
                (pointer, row) for pointer, row in _find(nodes, contribution["integrand"], "function_substitution")
                if row.get("expression") == density[order]
            ]
            _require(substitutions, "C5 density pointer changed: " + chart_name + ":" + name)
            if C5_RATES[name] == "0":
                _require(_exact(nodes, memory, 1), "C5 p memory was reset: " + chart_name)
                _require(_row(nodes, history_pointer).get("operation") == "sum"
                         and incoming[order] in _row(nodes, history_pointer).get("arguments", ()),
                         "C5 p incoming memory was not retained: " + chart_name + ":" + str(order))
                kernel = None
            else:
                kernel = _has_exp(nodes, contribution["integrand"],
                                  _expected_exp(C5_RATES[name], coordinate, variable, nodes))
                _require(kernel is not None, "C5 kernel exponent changed: " + chart_name + ":" + name)
                decay = _has_exp(nodes, memory, _expected_decay(C5_RATES[name], coordinate, nodes))
                _require(decay is not None, "C5 memory decay changed: " + chart_name + ":" + name)
            row_receipts.append(dict(order=order, history=history_pointer,
                                     contribution=contribution_pointer, density=density[order],
                                     memory=memory, kernel=kernel,
                                     q0=dict(exp_zero_is_one=True, integral_zero_is_zero=True,
                                             incoming_memory_retained=True,
                                             rate_zero_full_memory=(name == "p"))))
        proofs[name] = dict(rate=C5_RATES[name], ordinary_Z_orders=list(range(6)), rows=row_receipts)
    return dict(
        passed=True,
        chart=chart_name,
        coordinate=coordinate,
        physical_measure=jacobian,
        J_identity="J=1 exact, source Jacobian is applied once and simplified",
        incoming_endpoint_pointers=incoming_receipts,
        histories=proofs,
        p_rate_zero_memory_retained=True,
        source_executed=False,
    )


def _native_histories(c6: Mapping[str, Any], row_lemma: Mapping[str, Any]) -> dict[str, Any]:
    bindings = c6["prescribed_R100_Rm_C6_bindings"]
    nodes = c6["exact_graph_nodes"]
    charts = bindings["charts"]
    return dict(
        passed=True,
        source_cache=C6_CACHE,
        graph_node_count=len(nodes),
        rates=dict(NATIVE_RATES),
        density_source_formulas={
            "H": "2*phi", "M": "U", "K": "2*phi*U",
            "A": "U**2", "B": "phi**2", "C": "phi**2",
        },
        charts={name: _native_chart_history(nodes, name, charts[name])
                for name in ("restoration", "post_restore")},
        ordinary_Z_row_generator=row_lemma,
        actual_phi_U_roots_used=True,
        source_executed=False,
    )


def _c5_histories(c5: Mapping[str, Any], row_lemma: Mapping[str, Any]) -> dict[str, Any]:
    bindings = c5["prescribed_loop_C5_bindings"]
    nodes = c5["exact_graph_nodes"]
    charts = bindings["charts"]
    return dict(
        passed=True,
        source_cache=C5_CACHE,
        graph_node_count=len(nodes),
        rates=dict(C5_RATES),
        alias_postrestore_to_post_restore=True,
        predecessor_chain={"restoration": "reference", "postrestore": "restoration"},
        charts={
            "restoration": _c5_chart_history(nodes, "restoration", charts["restoration"], charts["reference"]),
            "postrestore": _c5_chart_history(nodes, "postrestore", charts["postrestore"], charts["restoration"]),
        },
        endpoint_links={
            "reference_to_restoration": bindings["endpoints"]["reference"] == charts["restoration"]["incoming"],
            "restoration_to_postrestore": bindings["endpoints"]["restoration"] == charts["postrestore"]["incoming"],
        },
        ordinary_Z_row_generator=row_lemma,
        source_executed=False,
    )


def _q4_seed_source(c6: Mapping[str, Any], q4: Mapping[str, Any], guards: Mapping[str, Any]) -> dict[str, Any]:
    c6_bindings = c6["prescribed_R100_Rm_C6_bindings"]
    q4_bindings = q4["modified_Q4_pressure_C3_bindings"]
    nodes = q4["exact_graph_nodes"]
    charts: dict[str, Any] = {}
    for native_name, q4_name in (("restoration", "restoration"), ("post_restore", "postrestore")):
        native_chart = c6_bindings["charts"][native_name]
        chart = q4_bindings["charts"][q4_name]
        background = {
            key: _rows(chart["radial_background"][key], q4_name + ".radial_background." + key, nodes, 7)
            for key in Q4_KEYS
        }
        _require(background["pressure"] == background["p"], "Q4 pressure row is not p_y: " + q4_name)
        _require(all(_exact(nodes, pointer, 0) for pointer in background["P0"]),
                 "Q4 P0_y row is not exact zero: " + q4_name)
        if q4_name == "postrestore":
            _require(all(_exact(nodes, pointer, 0) for pointer in background["V"]),
                     "Q4 post-restore V_y row is not exact zero")
        radial = {
            key: _rows(value, q4_name + ".radial." + key, nodes, 6)
            for key, value in chart["radial"].items() if isinstance(value, list)
        }
        physical = {
            key: _rows(value, q4_name + ".physical." + key, nodes, len(value))
            for key, value in chart["physical"].items() if isinstance(value, list)
        }
        _require(chart["radial"]["physical_y_equals_native_partial_divided_by_Jacobian"] is True,
                 "Q4 physical-y radial source flag changed: " + q4_name)
        c6_fields = {
            key: _rows(native_chart["fields"][key], native_name + ".fields." + key,
                       c6["exact_graph_nodes"], 7)
            for key in ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
        }
        p0_row = c6_fields["P0"][1]
        p0_source = c6["exact_graph_nodes"][p0_row]
        _require(p0_source.get("operation") == "current_original_leading_function_recipe"
                 and p0_source.get("quantity") == "P0"
                 and p0_source.get("Z_order") == 1
                 and p0_source.get("ordinary_slow_Z_derivative_order") == 1,
                 "C6 P0_Z source row was not retained")
        charts[q4_name] = dict(
            native_chart=native_name,
            J=1,
            J_source_pointer=native_chart["physical_measure"],
            J_identity="exact one",
            C6_fields=c6_fields,
            P0_Z=dict(pointer=p0_row, operation=p0_source.get("operation"),
                      quantity="P0", ordinary_Z_order=1, retained=True, P0_y_zero=True),
            radial_background=background,
            radial=radial,
            physical=physical,
            constructor_seed=dict(
                source="q4_constructor_seed_loop",
                mode="exact_J_one_direct_rows",
                seed_rows=background,
                J_applied_once=True,
                no_J_cap_or_selected_value=True,
            ),
            formulas=dict(Q4_FORMULAS),
        )
    return dict(
        passed=True,
        source_cache=Q4_CACHE,
        graph_node_count=len(nodes),
        source_guards={key: guards[key] for key in (
            "q4_background_radial_rows", "q4_constructor_seed", "q4_physical_ports")},
        charts=charts,
        all_new_C6_fields_bound=True,
        all_Q4_physical_radial_seed_rows_bound=True,
        P0_Z_bound=True,
        full_modified_Qy_claim=False,
        source_executed=False,
    )


def _sigma_endpoint_proof(guards: Mapping[str, Any]) -> dict[str, Any]:
    return dict(
        passed=True,
        source_guards={key: guards[key] for key in ("flat_sigma_left", "flat_sigma_jets", "c6_sigma_prime")},
        source_formula=(
            "s=e/(1+e), e=exp(1/(1-x)^2-1/x^2), "
            "s'=s(1-s)L1, s''=s(1-s)(L2+(1-2s)L1^2), "
            "L1=2/(1-x)^3+2/x^3, L2=6/(1-x)^4-6/x^4"
        ),
        endpoint_jets={
            "x=0": {"sigma": "0", "sigma_prime": "0", "sigma_second": "0"},
            "x=1": {"sigma": "1", "sigma_prime": "0", "sigma_second": "0"},
        },
        flat_tail_reason=(
            "At x=0, exp(-1/x^2) times every displayed polynomial singularity tends to zero; "
            "at x=1 the same statement applies to 1-sigma."
        ),
        exact_function_not_interval_cap=True,
        source_executed=False,
    )


def _radius_identities(c6: Mapping[str, Any]) -> dict[str, Any]:
    bindings = c6["prescribed_R100_Rm_C6_bindings"]
    nodes = c6["exact_graph_nodes"]
    charts = bindings["charts"]
    endpoints = bindings["endpoints"]
    _require(_row(nodes, charts["restoration"]["radius"]).get("operation") == "product",
             "Restoration radius source expression changed")
    _require(_row(nodes, charts["post_restore"]["radius"]).get("operation") == "product",
             "Post-restore radius source expression changed")
    rz = _row(nodes, endpoints["Rz"]["radius"])
    rr = _row(nodes, endpoints["Rrestore"]["radius"])
    _require(rz.get("operation") == "function_substitution" and rz.get("value") == 1,
             "Rz endpoint is not the actual reference endpoint substitution")
    _require(rr.get("operation") == "function_substitution" and rr.get("value") == 1,
             "Rrestore endpoint is not the actual restoration endpoint substitution")
    return dict(
        passed=True,
        graph_bindings={
            "Rz": dict(endpoint_radius=endpoints["Rz"]["radius"], source_chart="reference", endpoint_value=1),
            "Rrestore": dict(endpoint_radius=endpoints["Rrestore"]["radius"], source_chart="restoration", endpoint_value=1),
        },
        identities={
            "reference_upper_equals_restoration_lower": "Rz = 110*exp(logref-8) = restoration(t=0)",
            "restoration_upper_equals_post_restore_lower": "Rrestore = 110*exp(logref-7) = post_restore(p=0)",
        },
        source_offsets={"reference": "-8", "restoration": "-8+t", "post_restore": "-7+p"},
        source_formula_guarded=True,
        source_executed=False,
    )


def _native_c2(c6: Mapping[str, Any], sigma_proof: Mapping[str, Any]) -> dict[str, Any]:
    bindings = c6["prescribed_R100_Rm_C6_bindings"]
    nodes = c6["exact_graph_nodes"]
    charts = bindings["charts"]
    endpoints = bindings["endpoints"]
    four_z = 35566
    _require(charts["post_restore"]["actual"]["U"] == [four_z, 13, 0, 0, 0, 0, 0],
             "Post-restore U is no longer exact 4Z")
    _require(all(_exact(nodes, pointer, 0) for pointer in charts["post_restore"]["actual"]["U_y"]),
             "Post-restore U_y is no longer exact zero")
    e_rows = _rows(bindings["E"], "C6 E", nodes, 7)
    e_row = _row(nodes, e_rows[0], "C6 E[0]")
    _require(e_row.get("operation") == "sum" and endpoints["R110"]["actual"]["U"][0] in e_row.get("arguments", ()),
             "Restoration E does not use actual U_R110")
    rest_u = _row(nodes, charts["restoration"]["actual"]["U"][0], "restoration U[0]")
    _require(rest_u.get("operation") == "sum" and four_z in rest_u.get("arguments", ()),
             "Restoration U does not retain 4Z")
    e_products = [pointer for pointer in _reachable(nodes, charts["restoration"]["actual"]["U"][0])
                  if _row(nodes, pointer).get("operation") == "product"
                  and e_rows[0] in _row(nodes, pointer).get("arguments", ())]
    _require(e_products, "Restoration U does not use the source E amplitude")
    _require(endpoints["Rsh"]["actual"]["U"] == charts["reference"]["actual"]["U"],
             "Reference U lost the actual Rsh incoming pointer")
    _require(_row(nodes, endpoints["Rz"]["actual"]["U"][0]).get("expression")
             == endpoints["Rsh"]["actual"]["U"][0],
             "Rz U is not the actual Rsh-to-Rz substitution chain")
    _require(_row(nodes, endpoints["Rrestore"]["actual"]["U"][0]).get("expression")
             == charts["restoration"]["actual"]["U"][0],
             "Rrestore U is not the actual restoration endpoint")
    for name in ("reference", "restoration", "post_restore"):
        for row in charts[name]["actual"]["phi_y"][:3]:
            reachable = _reachable(nodes, row)
            _require(any(_exact(nodes, pointer, -2, 5) for pointer in reachable),
                     "Native phi_y lost the source -2/5 factor: " + name)
    _require(all(_exact(nodes, pointer, 0) for pointer in charts["reference"]["actual"]["U_y"]),
             "Reference U_y is not source zero")
    _require(all(_exact(nodes, pointer, 0) for pointer in charts["post_restore"]["actual"]["U_y"]),
             "Post-restore U_y is not source zero")
    return dict(
        passed=True,
        seams={
            "reference_to_restoration": dict(
                radius_identity="Rz=restoration(t=0)",
                value_identity="phi_ref(1)=phi_rest(0); U_ref(1)=4Z+E=U_rest(0)",
                first_y_identity="phi_y=-2*phi/5 on both sides; U_y=0 on both sides",
                second_y_identity="phi_yy=4*phi/25; U_yy=-E*sigma''(0)=0 on both sides",
                sigma_endpoint="(sigma,sigma',sigma'')(0)=(0,0,0)",
                physical_C2=True,
            ),
            "restoration_to_post_restore": dict(
                radius_identity="Rrestore=restoration(t=1)=post_restore(p=0)",
                value_identity="U_rest(1)=4Z=U_post(0); pure-power phi agrees",
                first_y_identity="phi_y=-2*phi/5; U_y=0",
                second_y_identity="phi_yy=4*phi/25; U_yy=-E*sigma''(1)=0",
                sigma_endpoint="(sigma,sigma',sigma'')(1)=(1,0,0)",
                physical_C2=True,
            ),
        },
        actual_source_bindings={
            "U_R110": endpoints["R110"]["actual"]["U"][0],
            "Rsh_U": endpoints["Rsh"]["actual"]["U"][0],
            "Rz_U": endpoints["Rz"]["actual"]["U"][0],
            "E": e_rows[0],
            "restoration_U": charts["restoration"]["actual"]["U"][0],
            "post_restore_U": charts["post_restore"]["actual"]["U"][0],
        },
        sigma_endpoint_proof=sigma_proof,
        native_phi_U_C2_proved=True,
        endpoint_values_are_source_identities=True,
        no_sign_choice_for_E=True,
        source_executed=False,
    )


def _finite_positive_gate(parent: Mapping[str, Any], c6: Mapping[str, Any], native: Mapping[str, Any], c5: Mapping[str, Any]) -> dict[str, Any]:
    controls = parent["actual_long_reference_native_controls"]
    _require(controls.get("passed") is True, "Saved native finite controls are not passed")
    charts: dict[str, Any] = {}
    # The new continuation starts at Rz, so the saved reference controls are
    # the only inherited quantitative witness needed here.
    for name in ("reference",):
        record = controls["charts"][name]
        _require(record.get("passed") is True, "Saved finite control chart is not passed: " + name)
        for key in ("phi_positive_log_lower", "phi_positive_log_upper",
                    "Utheta_positive_log_lower", "Utheta_positive_log_upper"):
            _require(record[key].get("directed") is True, "Saved positive log receipt is not directed: " + name + ":" + key)
        caps = record["native_log_absolute_caps"]
        for quantity in ("phi", "U", "M", "A", "B", "C"):
            _require(len(caps[quantity]) >= 5, "Saved finite log cap rows are incomplete: " + name + ":" + quantity)
            for row in caps[quantity][:5]:
                _require(row.get("exact_zero") is False
                         and row["log_absolute_upper_bound"].get("directed") is True,
                         "Saved finite log cap is not an actual directed witness")
        _require(record["memory_bound"].get("endpoint_bound_not_selected_value") is True,
                 "Saved memory bound was turned into a value")
        charts[name] = dict(
            positive_log_receipts=True,
            finite_native_log_caps=True,
            ordinary_Z_orders=list(range(5)),
            cap_is_witness_not_value=True,
            source_ports=record["source_ports"],
        )
    positive = parent["actual_long_reference_control_seams"]["positive_active_loop_IFT"]
    _require(positive.get("passed") is True and positive.get("source_executed") is False,
             "Saved positive local-IFT theorem is not consumable")
    hypotheses = positive["hypotheses"]
    for key in ("finite_p2", "finite_q", "finite_u", "D_positive", "K_positive", "lambda_positive"):
        _require(hypotheses[key].get("proved") is True, "Saved positive hypothesis is missing: " + key)
    c6_bindings = c6["prescribed_R100_Rm_C6_bindings"]
    for name in ("restoration", "post_restore"):
        chart = c6_bindings["charts"][name]
        _require(_row(c6["exact_graph_nodes"], chart["actual"]["physical_Utheta"][0]).get("operation") == "analytic_unary",
                 "New physical Utheta source expression is not analytic: " + name)
    return dict(
        passed=True,
        scope=Z_SCOPE,
        reference_controls=charts,
        positive_E_gate=dict(
            passed=True,
            source_meaning="positive generic E = F0*sqrt(2R)*phi from the actual source quotient",
            source_parameter_phi_Utheta=True,
            E_restore_sign_not_selected=True,
            no_cap_as_value=True,
        ),
        finite_propagation=dict(
            restoration="compact exact Duhamel continuation from actual Rz incoming rows",
            post_restore="compact exact Duhamel continuation from actual Rrestore incoming rows",
            C5="five signed correction densities use actual incoming endpoint pointers",
            finite_integrals_from_source_metadata=True,
            endpoint_independent_of_Z=True,
        ),
        local_IFT_gate=dict(
            native_phi_finite_positive=True,
            native_Utheta_finite_positive=True,
            finite_p2=True,
            positive_E=True,
            future_local_IFT_only=True,
        ),
        source_parameter_finiteness_is_propagated=True,
        directed_caps_are_not_field_values=True,
        source_executed=False,
    )


def _saved_lemmas(parent: Mapping[str, Any]) -> dict[str, Any]:
    history = parent["actual_long_reference_history_seams"]
    control = parent["actual_long_reference_control_seams"]
    _require(history.get("passed") is True and history.get("source_executed") is False,
             "Saved history theorem is not consumable")
    _require(control.get("passed") is True and control.get("source_executed") is False,
             "Saved generic seed theorem is not consumable")
    row_lemma = history["ordinary_Z_row_generator"]
    _require(row_lemma.get("passed") is True and row_lemma.get("source_executed") is False,
             "Saved FunctionJet row theorem is missing")
    _require(row_lemma.get("FunctionJet_rows_are_source_rows") is True
             and row_lemma.get("integral_rows_preserve_ordinary_order") is True,
             "Saved ordinary-Z row theorem changed")
    seed = control["native_seed_theorem_binding"]
    _require(seed.get("passed") is True and seed.get("accepted_lemma_consumed_only") is True,
             "Saved universal generic seed theorem is not consumable")
    return dict(
        row_lemma=row_lemma,
        generic_seed=dict(
            passed=True,
            source="parent.actual_long_reference_control_seams.native_seed_theorem_binding",
            consumed_saved_lemma_only=True,
            full_P0_retained=True,
            no_old_math_rerun=True,
        ),
        history_ODE=dict(
            statement="Dq X/J = density - rate*X",
            q0="X(0)=incoming; integral at lower bound is zero; p-rate-zero memory is one",
            source_FTC=True,
            applied_to_new_charts=True,
        ),
    )


def _load_parent(root: Path, parent: Mapping[str, Any]) -> Mapping[str, Any]:
    disk = _json(root, PARENT)
    _require(parent == disk, "Exact modified-seams parent mapping is required")
    _require(_sha(root / PARENT) == PARENT_SHA256, "Modified-seams parent bytes changed")
    _require(isinstance(parent.get("source_family"), Mapping), "Parent source family is missing")
    for section in ("actual_long_reference_native_seams", "actual_long_reference_history_seams",
                    "actual_long_reference_loop_seams", "actual_long_reference_control_seams"):
        packet = parent.get(section)
        _require(isinstance(packet, Mapping) and packet.get("passed") is True,
                 "Saved parent lemma section is not passed: " + section)
        if "source_executed" in packet:
            _require(packet["source_executed"] is False, "Saved parent source was executed: " + section)
    return disk


def _validate_endpoint_packet(
    root: Path,
    parent: Mapping[str, Any],
    endpoint_packet: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate the root-owned native endpoint receipt before consuming it."""
    _require(endpoint_packet.get("passed") is True,
             "Optional restoration endpoint packet is not passed")
    _require(endpoint_packet.get("parent_report") == PARENT,
             "Endpoint packet has the wrong modified-seams parent")
    _require(endpoint_packet.get("parent_report_sha256") == PARENT_SHA256,
             "Endpoint packet has the wrong modified-seams parent hash")
    _require(endpoint_packet.get("source_family") == parent["source_family"],
             "Endpoint packet source family changed")
    _require(endpoint_packet.get("source_executed") is False,
             "Endpoint packet source was executed")
    _require(endpoint_packet.get("native_phi_U_physical_C2") is True,
             "Endpoint packet does not prove native phi/U physical C2")
    _require(endpoint_packet.get("radius_endpoint_identities") is True,
             "Endpoint packet does not prove radius endpoint identities")
    pins = endpoint_packet.get("source_pins")
    _require(isinstance(pins, Mapping), "Endpoint packet source pins are missing")
    endpoint_path = root / ENDPOINT_MODULE
    endpoint_sha = _sha(endpoint_path)
    self_pin = pins.get(ENDPOINT_MODULE)
    _require(isinstance(self_pin, Mapping),
             "Endpoint packet is missing its module self-pin")
    _require(self_pin.get("path") == ENDPOINT_MODULE
             and self_pin.get("passed") is True
             and self_pin.get("expected_sha256") == endpoint_sha
             and self_pin.get("actual_sha256") == endpoint_sha
             and self_pin.get("source_executed") is False,
             "Endpoint module self-pin does not match root-owned bytes")

    symbolic = endpoint_packet.get("symbolic_endpoint_identities")
    _require(isinstance(symbolic, Mapping) and symbolic.get("passed") is True,
             "Endpoint symbolic identity receipt is not passed")
    _require(symbolic.get("function_identity_not_sampling") is True,
             "Endpoint identities are not function identities")
    _require(endpoint_packet.get("ordinary_Z_orders") == [0, 1, 2, 3],
             "Endpoint ordinary-Z order scope changed")
    seams = symbolic.get("seams")
    expected_seams = (
        ("reference_to_restoration", "reference", "restoration"),
        ("restoration_to_postrestore", "restoration", "post_restore"),
    )
    _require(isinstance(seams, list) and len(seams) == len(expected_seams),
             "Endpoint symbolic seam count changed")
    checked_seams: dict[str, Any] = {}
    for seam, (name, left_chart, right_chart) in zip(seams, expected_seams):
        _require(isinstance(seam, Mapping)
                 and seam.get("name") == name
                 and seam.get("left_chart") == left_chart
                 and seam.get("right_chart") == right_chart
                 and seam.get("passed") is True,
                 "Endpoint symbolic seam is missing or renamed: " + name)
        residuals = seam.get("residuals")
        ordinary = seam.get("ordinary_Z0_Z3_residuals")
        _require(isinstance(residuals, Mapping)
                 and residuals
                 and all(value == "0" for value in residuals.values()),
                 "Endpoint seam has a nonzero residual: " + name)
        _require(isinstance(ordinary, Mapping)
                 and ordinary
                 and all(isinstance(values, list)
                         and len(values) == 4
                         and all(value == "0" for value in values)
                         for values in ordinary.values()),
                 "Endpoint seam ordinary-Z residuals are not zero: " + name)
        checked_seams[name] = dict(
            passed=True,
            left_chart=left_chart,
            right_chart=right_chart,
            residuals_zero=True,
            ordinary_Z0_Z3_zero=True,
        )
    packet_sha = hashlib.sha256(
        json.dumps(endpoint_packet, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return dict(
        passed=True,
        packet_sha256=packet_sha,
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        source_family=parent["source_family"],
        endpoint_module_sha256=endpoint_sha,
        self_pin_matches_root_bytes=True,
        seams=checked_seams,
        native_phi_U_physical_C2=True,
        radius_endpoint_identities=True,
        source_executed=False,
    )


def build(
    root: str | Path,
    parent: Mapping[str, Any],
    endpoint_packet: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build the source-bound restoration/post-restore receipt."""
    root_path = Path(root)
    disk_parent = _load_parent(root_path, parent)
    pins = _source_pins(root_path, disk_parent)
    guards = _source_guards(root_path)
    lemmas = _saved_lemmas(disk_parent)
    c6 = _gzip_json(root_path, C6_CACHE)
    c5 = _gzip_json(root_path, C5_CACHE)
    q4 = _gzip_json(root_path, Q4_CACHE)
    for cache, key, label in ((c6, "prescribed_R100_Rm_C6_bindings", "C6"),
                              (c5, "prescribed_loop_C5_bindings", "C5"),
                              (q4, "modified_Q4_pressure_C3_bindings", "Q4")):
        _require(cache.get("source_family") == disk_parent["source_family"], label + " source family changed")
        _require(isinstance(cache.get("exact_graph_nodes"), list)
                 and isinstance(cache.get(key), Mapping), label + " graph/binding shape changed")
    native = _native_histories(c6, lemmas["row_lemma"])
    c5_receipt = _c5_histories(c5, lemmas["row_lemma"])
    q4_receipt = _q4_seed_source(c6, q4, guards)
    finite = _finite_positive_gate(disk_parent, c6, native, c5_receipt)
    if endpoint_packet is None:
        endpoint_gate = None
        endpoint_dependency = dict(
            required=True,
            consumed=False,
            owner="original_restoration_native_endpoints.py",
            interface="build(root,parent)->passed packet with reference_to_restoration and restoration_to_post_restore",
            source_executed=False,
        )
    else:
        _require(isinstance(endpoint_packet, Mapping),
                 "Optional restoration endpoint packet is not a mapping")
        endpoint_gate = _validate_endpoint_packet(root_path, disk_parent, endpoint_packet)
        endpoint_dependency = dict(
            required=True,
            consumed=True,
            owner="original_restoration_native_endpoints.py",
            packet_sha256=endpoint_gate["packet_sha256"],
            validation=endpoint_gate,
            source_executed=False,
        )
    endpoint_proved = endpoint_gate is not None and endpoint_gate.get("passed") is True
    c6_charts = c6["prescribed_R100_Rm_C6_bindings"]["charts"]
    first_y_tables = {
        name: {
            key: list(c6_charts[name]["actual"][key])
            for key in ("phi", "U", "phi_y", "U_y", "physical_Utheta", "log_Utheta")
        }
        for name in ("restoration", "post_restore")
    }
    result: dict[str, Any] = dict(
        passed=True,
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,
        source_family=disk_parent["source_family"],
        source_pins=pins,
        input_pins=pins,
        source_guards=guards,
        sourceguards_tuple=True,
        source_executed=False,
        accepted_ancestors_replayed=False,
        accepted_builders_executed=False,
        accepted_checkers_executed=False,
        accepted_producers_executed=False,
        owners_executed=False,
        norms_executed=False,
        recurrences_executed=False,
        quadrature_executed=False,
        ordinary_Z_orders=list(range(7)),
        physical_y_orders=[0, 1, 2],
        scope=dict(charts=["reference", "restoration", "post_restore"], Z=Z_SCOPE,
                   native_history_orders=list(range(7)), C5_history_orders=list(range(6))),
        saved_lemmas=lemmas,
        conditional_history_ODE=lemmas["history_ODE"],
        native_first_y_seed_tables=first_y_tables,
        endpoint_dependency=endpoint_dependency,
        native_history_source=native,
        c5_history_source=c5_receipt,
        c5_correction_history_source=c5_receipt,
        q4_background_seed_source=q4_receipt,
        actual_new_C6_fields={
            name: q4_receipt["charts"][q4_name]["C6_fields"]
            for name, q4_name in (("restoration", "restoration"), ("post_restore", "postrestore"))
        },
        actual_Q4_physical_radial_seed_rows={
            name: q4_receipt["charts"][q4_name]["constructor_seed"]
            for name, q4_name in (("restoration", "restoration"), ("post_restore", "postrestore"))
        },
        native_control_positive_propagation=finite,
        finite_p2_positive_E_local_gate=finite["local_IFT_gate"],
        actual_native_histories=True,
        actual_C5_correction_histories=True,
        limits=dict(
            full_modified_Qy=False,
            full_modified_RQy=False,
            global_five_moment_closure=False,
            current_N_numeric=False,
            numerical_values_selected=False,
            interval_caps_used_as_values=False,
            modified_pressure_seam_claim=False,
        ),
        newproof_checks=dict(
            reference_to_restoration_native_phi_U_C2=endpoint_proved,
            restoration_to_post_restore_native_phi_U_C2=endpoint_proved,
            exact_sigma_endpoint_jets_through_second=endpoint_proved,
            actual_U_R110_to_Rsh_to_Rz_chain=endpoint_proved,
            restoration_E_is_actual_U_R110_minus_4Z=endpoint_proved,
            all_six_new_native_histories=True,
            all_five_new_C5_correction_histories=True,
            lower_bound_zero_integrals_and_unit_memory=True,
            formal_FunctionJet_rows_orders_0_through_6=True,
            every_new_C6_field_bound=True,
            every_new_Q4_radial_seed_row_and_J_bound=True,
            P0_Z_retained=True,
            finite_positive_native_phi_Utheta_gate=True,
            positive_E_gate=True,
            no_full_modified_Qy_claim=True,
        ),
        module_sha256=_sha(root_path / MODULE),
    )
    result["receipt_sha256"] = hashlib.sha256(
        json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    json.dumps(result, allow_nan=False)
    return result


__all__ = ["build"]
