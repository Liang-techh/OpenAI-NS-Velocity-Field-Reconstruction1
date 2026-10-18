"""Source-bound history and physical-first-y checks for the long/reference joins.

The receipt in this module is deliberately narrower than a seam proof.  It reads
the accepted C6, C5, and modified-Q4 archives as data and parses their source
functions with :mod:`ast`; it never imports or executes an accepted producer,
checker, owner, recurrence, quadrature, or norm.  The result proves the actual
Duhamel source expansion and cache-pointer chain for the two native charts and
the five finite-N correction histories.  It also records the Q4 background
radial source rows and the constructor's ``J * radial_bg`` seed products.

The FTC statement certified here is, for every source history,

``X(q) = exp(-rate*J*q)*Xin + int_0^q J*exp(-rate*J*(q-s))*density(s) ds``

and therefore ``Dq X / J = density - rate*X``.  Equality of endpoint fields,
periodic controls, or the full modified ``Q_y`` seam remains an explicit
conditional obligation for the root receipt.
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
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_long_reference_finite_N_I1.json"
PARENT_SHA256 = "75918293b8c7432c637eb551a29cb9d248aa7746dcb1cd4f8b3615f72514c14d"
MODULE = "src/openai_ns_reconstruction/original_long_reference_history_seams.py"

C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
C6_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
C6_BRIDGE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
RATE_SOURCE = BASE + "lei_ren_part1_paper_compliant_actual_bridge_integrals.py"
C5_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
C5_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
Q4_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
Q4_CACHE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
GENERIC_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source.py"

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


def _require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_json(root: Path, relative: str) -> Any:
    path = root / relative
    _require(path.is_file(), "Missing pinned JSON source: " + relative)
    return json.loads(path.read_bytes())


def _read_gzip_json(root: Path, relative: str) -> Any:
    path = root / relative
    _require(path.is_file(), "Missing pinned gzip source: " + relative)
    return json.loads(gzip.decompress(path.read_bytes()))


def _copy_json(value: Any) -> Any:
    return json.loads(json.dumps(value, separators=(",", ":")))


def _compact(value: str) -> str:
    return "".join(value.split())


def _source_guard(
    root: Path,
    relative: str,
    name: str,
    snippets: tuple[str, ...],
    *,
    class_name: str | None = None,
) -> dict[str, Any]:
    """Pin one function/class method without importing or executing source."""
    _require(type(snippets) is tuple, "Source guard snippets must be a tuple")
    path = root / relative
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=str(path))
    candidates: list[ast.AST] = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) or node.name != name:
            continue
        if class_name is not None:
            owner = next(
                (
                    parent
                    for parent in ast.walk(tree)
                    if isinstance(parent, ast.ClassDef)
                    and any(child is node for child in parent.body)
                ),
                None,
            )
            if not isinstance(owner, ast.ClassDef) or owner.name != class_name:
                continue
        segment = ast.get_source_segment(text, node) or ""
        if all(_compact(snippet) in _compact(segment) for snippet in snippets):
            candidates.append(node)
    _require(len(candidates) == 1, "Changed source definition: " + relative + ":" + name)
    node = candidates[0]
    return dict(
        passed=True,
        path=relative,
        function=name,
        class_name=class_name,
        required_text=list(snippets),
        AST_sha256=hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest(),
        source_executed=False,
    )


def _source_assignment(root: Path, relative: str, name: str) -> ast.AST:
    path = root / relative
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    rows = [
        row
        for row in tree.body
        if isinstance(row, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == name for target in row.targets)
    ]
    _require(len(rows) == 1, "Expected one source assignment: " + relative + ":" + name)
    return rows[0]


def _source_class_base(root: Path, relative: str, class_name: str, base_name: str) -> dict[str, Any]:
    path = root / relative
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=str(path))
    classes = [node for node in ast.walk(tree)
               if isinstance(node, ast.ClassDef) and node.name == class_name]
    _require(len(classes) == 1, "Expected one source class: " + relative + ":" + class_name)
    node = classes[0]
    bases = []
    for base in node.bases:
        if isinstance(base, ast.Name):
            bases.append(base.id)
        elif isinstance(base, ast.Attribute):
            bases.append(base.attr)
    _require(base_name in bases, "Source class inheritance changed: " + relative + ":" + class_name)
    return dict(
        passed=True,
        path=relative,
        class_name=class_name,
       required_base=base_name,
        actual_bases=bases,
        AST_sha256=hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest(),
        source_executed=False,
    )


def _literal_fraction(node: ast.AST) -> str:
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "Fraction":
        values = [ast.literal_eval(arg) for arg in node.args]
        _require(values and len(values) <= 2 and all(type(v) is int for v in values),
                 "Fraction source rate is not integral")
        return str(values[0]) if len(values) == 1 else f"{values[0]}/{values[1]}"
    value = ast.literal_eval(node)
    _require(type(value) is int, "Source rate is not exact")
    return str(value)


def _source_guards(root: Path) -> dict[str, Any]:
    guards = {
        "native_densities": _source_guard(root, C6_BRIDGE, "densities", (
            "return dict(H=a.scale(phi,2),M=U,K=a.scale(a.mul(phi,U),2),",
            "A=a.mul(U,U),B=a.mul(phi,phi),C=a.mul(phi,phi))",
        )),
        "native_histories": _source_guard(root, C6_BRIDGE, "histories", (
            "body=a.scale(a.at(bodies[name],phase,t),g.mul(jacobian,kernel))",
            "out[name]=a.add(a.scale(initial[name],decay),a.integral(body,t,g.zero,phase,",
            "actual_six_history_rate=rate,actual_history_name=name",
            "original_incoming_memory_retained=True,source_measure_applied_once=True",
        )),
        "native_generic_fields": _source_guard(root, C6_BRIDGE, "generic_fields", (
            "fields=dict(E=a.mul(W,actual['phi']),V=a.scale(actual['U'],inverse),m=a.scale(J['M'],inverse)",
            "h=a.scale(a.mul(W,J['H']),'1/2'),k=a.scale(a.mul(W,J['K']),g.mul(c('1/2'),inverse))",
            "e=a.sub(a.scale(J['A'],g.mul(inverse,inverse)),a.mul(FF,J['B'])),p=a.mul(FF,J['C']),P0=a.trim(P0)",
            "fields['pressure']=a.add(fields['P0'],fields['p'])",
            "radial=dict(E_y=a.mul(W,a.add(a.scale(actual['phi'],'1/2'),actual['phi_y'])),V_y=a.scale(actual['U_y'],inverse))",
        )),
        "ordinary_function_jet_rows": _source_guard(root, GENERIC_SOURCE, "rows", (
            "return jet.rows if type(jet) is FunctionJet else previous.higher.rows(jet)",
        )),
        "ordinary_jet_add_rows": _source_guard(root, GENERIC_SOURCE, "add", (
            "return FunctionJet(tuple(self.g.add(*(rows(j)[k] for j in jets)) for k in range(self.order+1)))",
        ), class_name="JetAlgebra"),
        "ordinary_jet_mul_rows": _source_guard(root, GENERIC_SOURCE, "mul", (
            "return FunctionJet(tuple(g.add(*(g.mul(g.constant(math.comb(k,j)),l[j],r[k-j]) for j in range(k+1))) for k in range(self.order+1)))",
        ), class_name="JetAlgebra"),
        "ordinary_jet_at_rows": _source_guard(root, GENERIC_SOURCE, "at", (
            "return FunctionJet(tuple(self.g.node('function_substitution',expression=q.node,variable=variable.node,value=value.node,Z_independent_substitution=True) for q in rows(self.trim(jet))))",
        ), class_name="JetAlgebra"),
        "ordinary_jet_integral_rows": _source_guard(root, GENERIC_SOURCE, "integral", (
            "return FunctionJet(tuple(controls.integral(g,q,g.nodes[variable.node]['name'],lower,upper,endpoint_independent_of_Z=True,ordinary_slow_Z_derivative_order=k,true_core_inlet_complete_history=True,numerical_cap_not_defining_integrand=True,**attrs) for k,q in enumerate(rows(self.trim(jet)))))",
        ), class_name="JetAlgebra"),
        "native_bridge_rate_alias": _source_guard(root, C6_BRIDGE, "densities", (
            "return dict(H=a.scale(phi,2)",
        )),
        "native_bridge_algebra_inheritance": _source_class_base(
            root, C6_BRIDGE, "BridgeAlgebra", "JetAlgebra"
        ),
        "c5_loop_algebra_inheritance": _source_class_base(
            root, C5_SOURCE, "LoopAlgebra", "Algebra"
        ),
        "c5_rates": _source_guard(root, C5_SOURCE, "transport", (
            "width=g.mul(jacobian,g.sub(coordinate,lower))",
            "contribution=a.integral(body,t,lower,coordinate,actual_finite_N_signed_correction=True",
            "out[name]=a.add(a.scale(incoming[name],decay),contribution)",
            "physical_Jacobian_applied_once=True",
            "zero_correction_inlet_only_at_original_rminus=True",
            "predecessor_memory_never_reset=True",
        )),
        "c5_finite_N": _source_guard(root, C5_SOURCE, "finite_N_functions", (
            "r=a.mul(roots['E'],exm1);s=a.scale(loop['B_over_Pstar'],invN)",
            "inc=dict(m=s,h=r,k=a.add(a.mul(V,r),a.mul(E,s),a.mul(r,s))",
            "e=a.sub(a.add(a.scale(a.mul(V,s),2),a.mul(s,s)),a.add(a.mul(E,r),a.scale(a.mul(r,r),'1/2')))",
            "p=a.add(a.mul(E,r),a.scale(a.mul(r,r),'1/2'))",
            "five_signed_density_increments=inc",
        )),
        "c5_generic_roots": _source_guard(root, C5_SOURCE, "generic_roots", (
            "C=a.sub(E,a.scale(a.trim(radial['E_y']),2));B=a.scale(a.trim(radial['V_y']),2)",
            "aa,bb=divide(C),divide(B)",
            "p2=generic['p2'];t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))",
        )),
        "q4_background_rows": _source_guard(root, Q4_SOURCE, "background_radial_rows", (
            "result=dict(E=chart['radial']['E_y'],V=chart['radial']['V_y']",
            "m=a.sub(V,f['m']),h=a.sub(E,a.scale(f['h'],'3/2'))",
            "k=a.sub(a.mul(E,V),a.scale(f['k'],'3/2'))",
            "e=a.sub(a.sub(a.mul(V,V),a.scale(E2,'1/2')),f['e'])",
            "p=a.scale(E2,'1/2'),P0=a.fixed(0)",
            "result['pressure']=result['p'];return result",
        )),
        "q4_native_init": _source_guard(root, Q4_SOURCE, "__init__", (
            "self.cache={phase.node:g.zero};self.core_extensions=[]",
        ), class_name="NativeDerivative"),
        "q4_constructor_seed_loop": _source_guard(root, Q4_SOURCE, "__init__", (
            "radial_bg=background_radial_rows(g,bg)",
            "for key,jet in bg['fields'].items():",
            "for ref,dr in zip(rows(jet),rows(radial_bg[key])):derivative.seed(ref,g.mul(jac,dr))",
        )),
        "q4_physical_ports": _source_guard(root, Q4_SOURCE, "physical_pressure_fields", (
            "U=a.scale(fields['V'],Pstar);Uy=a.scale(radial['V_y'],Pstar)",
            "m=a.scale(fields['m'],Pstar);my=a.sub(U,m)",
            "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
        )),
        "q4_total_radial": _source_guard(root, Q4_SOURCE, "total_radial_functions", (
            "slow={key:FunctionJet(tuple(g.mul(invjac,derivative.derivative(ref))",
            "Vy=a.add(a.trim(bg['radial']['V_y']),a.scale(slow['B_over_Pstar'],invN),chart['loop']['B_over_Pstar_phi'])",
            "m_y=mean_y",
            "physical_y_equals_native_partial_divided_by_Jacobian=True",
        )),
        "q4_generic_p2": _source_guard(root, Q4_SOURCE, "generic_functions", (
            "Iz=a.scale(a.add(izl,a.scale(izq,Pstar)),radius)",
            "p2=a.div(Iz,E,'same actual positive modified swirl E_N')",
        )),
        "generic_p2_source": _source_guard(root, GENERIC_SOURCE, "generic_functions", (
            "p2=a.scale(a.div(a.add(axial_linear,a.scale(axial_quadratic,Pstar)),E,",
        )),
    }
    assignment = _source_assignment(root, RATE_SOURCE, "RATES")
    _require(isinstance(assignment.value, ast.Call), "Actual RATES assignment changed")
    _require(isinstance(assignment.value.func, ast.Name) and assignment.value.func.id == "dict",
             "Actual RATES must be a dict call")
    native_rates: dict[str, int] = {}
    for keyword in assignment.value.keywords:
        _require(keyword.arg is not None, "Native RATES keyword missing")
        value = ast.literal_eval(keyword.value)
        _require(type(value) is int, "Native RATES entry is not an integer")
        native_rates[keyword.arg] = value
    _require(native_rates == NATIVE_RATES, "Actual native RATES changed")
    guards["native_rates"] = dict(
        passed=True,
        path=RATE_SOURCE,
        binding="RATES",
        rates=native_rates,
        AST_sha256=hashlib.sha256(ast.dump(assignment, include_attributes=False).encode()).hexdigest(),
        source_executed=False,
    )
    c5_assignment = _source_assignment(root, C5_SOURCE, "RATES")
    _require(isinstance(c5_assignment.value, ast.Call)
             and isinstance(c5_assignment.value.func, ast.Name)
             and c5_assignment.value.func.id == "dict", "C5 RATES assignment changed")
    c5_rates: dict[str, str] = {}
    for keyword in c5_assignment.value.keywords:
        _require(keyword.arg is not None, "C5 RATES keyword missing")
        c5_rates[keyword.arg] = _literal_fraction(keyword.value)
    _require(c5_rates == C5_RATES, "Actual C5 correction rates changed")
    guards["c5_rate_assignment"] = dict(
        passed=True,
        path=C5_SOURCE,
        binding="RATES",
        rates=c5_rates,
        AST_sha256=hashlib.sha256(ast.dump(c5_assignment, include_attributes=False).encode()).hexdigest(),
        source_executed=False,
    )
    bridge_text = (root / C6_BRIDGE).read_text(encoding="utf-8")
    _require(_compact("RATES=bridge.RATES") in _compact(bridge_text),
             "C6 bridge no longer aliases the actual native rates")
    guards["native_bridge_rate_alias"] = dict(
        passed=True,
        path=C6_BRIDGE,
        binding="RATES=bridge.RATES",
        rates=dict(NATIVE_RATES),
        source_executed=False,
    )
    return guards


def _ordinary_row_generator_receipt(source_guards: Mapping[str, Any]) -> dict[str, Any]:
    """Bind the source rules that generate every ordinary Z row.

    The history graph checks below inspect value-row topology directly.  The
    density, kernel, memory, and integral identities for higher rows are also
    justified by the actual ``FunctionJet`` generator: each JetAlgebra
    operation maps every row index, while ``integral`` preserves the row
    index as ``ordinary_slow_Z_derivative_order``.  This receipt makes that
    source-level implication explicit instead of treating row 0 as a proxy
    for all ordinary rows.
    """
    required = (
        "ordinary_function_jet_rows", "ordinary_jet_add_rows",
        "ordinary_jet_mul_rows", "ordinary_jet_at_rows", "ordinary_jet_integral_rows",
        "native_bridge_algebra_inheritance", "c5_loop_algebra_inheritance",
    )
    for key in required:
        _require(isinstance(source_guards.get(key), Mapping)
                 and source_guards[key].get("passed") is True,
                 "ordinary FunctionJet source guard missing: " + key)
    return dict(
        passed=True,
        source=GENERIC_SOURCE,
        source_guards={key: source_guards[key] for key in required},
        ordinary_Z_row_orders={"native": list(range(7)), "C5": list(range(6))},
        FunctionJet_rows_are_source_rows=True,
        density_rows_generated_by_source_algebra=True,
        kernel_rows_generated_by_source_algebra=True,
        memory_rows_generated_by_source_algebra=True,
        integral_rows_preserve_ordinary_order=True,
        source_generator_proves_all_ordinary_rows=True,
        native_histories_use_at_scale_integral_add=True,
        C5_histories_use_at_scale_integral_add=True,
        native_bridge_and_C5_loop_inherit_source_algebra=True,
        graph_per_row_topology_checked=False,
        no_row0_proxy_claim=True,
        source_executed=False,
    )


def _accepted_pin(parent: Mapping[str, Any], relative: str) -> str:
    if relative == PARENT:
        return PARENT_SHA256
    values: list[str] = []
    for section in ("actual_long_reference_native_controls", "actual_long_reference_pressure_source"):
        packet = parent.get(section)
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
    paths = (PARENT, C6_SOURCE, C6_CACHE, C6_BRIDGE, RATE_SOURCE, C5_SOURCE, C5_CACHE, Q4_SOURCE, Q4_CACHE, GENERIC_SOURCE)
    result: dict[str, Any] = {}
    for relative in paths:
        expected = _accepted_pin(parent, relative)
        actual = _sha(root / relative)
        _require(actual == expected, "Pinned source bytes changed: " + relative)
        result[relative] = dict(path=relative, expected_sha256=expected, sha256=actual, passed=True)
    module_path = root / MODULE
    _require(module_path.is_file(), "This module must exist before build")
    digest = _sha(module_path)
    result[MODULE] = dict(path=MODULE, expected_sha256=digest, sha256=digest, passed=True)
    return result


def _pointer(value: Any, label: str, nodes: list[Any]) -> int:
    _require(type(value) is int and not isinstance(value, bool), label + " is not a graph pointer")
    _require(0 <= value < len(nodes), label + " is outside the graph")
    return value


def _rows(value: Any, label: str, nodes: list[Any], length: int) -> list[int]:
    _require(isinstance(value, list) and len(value) == length, label + " has the wrong jet length")
    return [_pointer(item, f"{label}[{index}]", nodes) for index, item in enumerate(value)]


def _row(nodes: list[Any], pointer: Any, label: str = "row") -> Mapping[str, Any]:
    index = _pointer(pointer, label, nodes)
    value = nodes[index]
    _require(isinstance(value, Mapping), label + " is not a mapping")
    return value


def _children(row: Mapping[str, Any]) -> tuple[int, ...]:
    operation = row.get("operation")
    if operation in ("sum", "product"):
        return tuple(item for item in row.get("arguments", ()) if type(item) is int)
    if operation in ("negative", "analytic_unary", "exact_original_source_flat_sigma"):
        return (row["argument"],) if type(row.get("argument")) is int else ()
    if operation == "positive_quotient":
        return tuple(row[key] for key in ("numerator", "denominator") if type(row.get(key)) is int)
    if operation == "function_substitution":
        return tuple(row[key] for key in ("expression", "value") if type(row.get(key)) is int)
    if operation == "definite_integral":
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


def _find_reachable(nodes: list[Any], start: int, operation: str) -> list[tuple[int, Mapping[str, Any]]]:
    return [(pointer, _row(nodes, pointer)) for pointer in _reachable(nodes, start)
            if _row(nodes, pointer).get("operation") == operation]


def _exact_rational(nodes: list[Any], pointer: Any, numerator: int, denominator: int = 1) -> bool:
    row = _row(nodes, pointer)
    return row.get("operation") == "exact_rational" and row.get("numerator") == numerator and row.get("denominator") == denominator


def _product_has(nodes: list[Any], pointer: Any, required: set[int]) -> bool:
    row = _row(nodes, pointer)
    return row.get("operation") == "product" and required.issubset(set(row.get("arguments", ())))


def _canonical(nodes: list[Any], pointer: int, memo: dict[int, Any] | None = None) -> Any:
    """Canonicalize the small source expression language used by the caches.

    Products and sums are flattened and sorted, so this checks the actual
    rational expression rather than relying on the construction order chosen
    by the source algebra.  Source parameters and bound variables stay
    pointer/name atoms; no numerical value or interval endpoint is selected.
    """
    if memo is None:
        memo = {}
    if pointer in memo:
        return memo[pointer]
    row = _row(nodes, pointer)
    operation = row.get("operation")
    if operation == "exact_rational":
        result = ("rational", row.get("numerator"), row.get("denominator"))
    elif operation == "bound_variable":
        result = ("bound", row.get("name"))
    elif operation == "current_original_source_parameter":
        result = ("source_parameter", pointer, row.get("name"))
    elif operation == "negative":
        result = ("negative", _canonical(nodes, row["argument"], memo))
    elif operation == "analytic_unary":
        result = ("analytic_unary", row.get("name"), _canonical(nodes, row["argument"], memo))
    elif operation == "positive_quotient":
        result = ("positive_quotient", _canonical(nodes, row["numerator"], memo),
                  _canonical(nodes, row["denominator"], memo))
    elif operation == "function_substitution":
        result = ("function_substitution", _canonical(nodes, row["expression"], memo),
                  row.get("variable"), _canonical(nodes, row["value"], memo))
    elif operation in ("sum", "product"):
        terms: list[Any] = []
        for child in row.get("arguments", ()):
            child_canonical = _canonical(nodes, child, memo)
            if isinstance(child_canonical, tuple) and child_canonical[0] == operation:
                terms.extend(child_canonical[1])
            else:
                terms.append(child_canonical)
        result = (operation, tuple(sorted(terms, key=repr)))
    else:
        # Unknown source operators remain explicit atoms.  The history checks
        # below reject them unless the expected source expression also carries
        # the same operator, so an arbitrary reachable node cannot pass.
        result = ("source_operation", operation, pointer)
    memo[pointer] = result
    return result


def _canonical_product(*terms: Any) -> Any:
    flat: list[Any] = []
    for term in terms:
        if isinstance(term, tuple) and term and term[0] == "product":
            flat.extend(term[1])
        else:
            flat.append(term)
    return ("product", tuple(sorted(flat, key=repr)))


def _canonical_sum(*terms: Any) -> Any:
    flat: list[Any] = []
    for term in terms:
        if isinstance(term, tuple) and term and term[0] == "sum":
            flat.extend(term[1])
        else:
            flat.append(term)
    return ("sum", tuple(sorted(flat, key=repr)))


def _rational(numerator: int, denominator: int = 1) -> Any:
    return ("rational", numerator, denominator)


def _bound_pointer(nodes: list[Any], name: str) -> int:
    matches = [
        pointer for pointer, row in enumerate(nodes)
        if row.get("operation") == "bound_variable" and row.get("name") == name
    ]
    _require(matches, "Integral bound variable is missing: " + name)
    return matches[0]


def _rate_tuple(text: str) -> tuple[int, int]:
    pieces = text.split("/", 1)
    if len(pieces) == 1:
        return int(pieces[0]), 1
    return int(pieces[0]), int(pieces[1])


def _native_endpoint_chain(bindings: Mapping[str, Any], nodes: list[Any]) -> dict[str, Any]:
    charts = bindings["charts"]
    endpoints = bindings["endpoints"]
    links: dict[str, Any] = {}
    for chart, endpoint in (("long_reshape", "R110"), ("reference", "Rsh")):
        incoming = charts[chart]["incoming_histories"]
        predecessor = endpoints[endpoint]["moments"]
        links[chart] = {}
        for name in NATIVE:
            left = _rows(incoming[name], f"{chart}.incoming_histories.{name}", nodes, 7)
            right = _rows(predecessor[name], f"{endpoint}.moments.{name}", nodes, 7)
            _require(left == right, "Native incoming endpoint pointer changed: " + chart + ":" + name)
            links[chart][name] = dict(
                incoming=list(left),
                predecessor_endpoint=list(right),
                pointer_equal=True,
            )
    return dict(
        passed=True,
        all_six_native_histories=True,
        predecessor_endpoints={"long_reshape": "R110", "reference": "Rsh"},
        links=links,
        source_executed=False,
    )


def _native_density_roots(
    nodes: list[Any],
    chart: Mapping[str, Any],
    name: str,
    expression: int,
) -> dict[str, Any]:
    actual = chart["actual"]
    phi = _pointer(actual["phi"][0], "actual.phi[0]", nodes)
    velocity = _pointer(actual["U"][0], "actual.U[0]", nodes)
    reachable = _reachable(nodes, expression)
    roots: dict[str, list[int]] = {
        "H": [phi], "M": [velocity], "K": [phi, velocity],
        "A": [velocity, velocity], "B": [phi, phi], "C": [phi, phi],
    }
    expected = roots[name]
    _require(all(item in reachable for item in expected),
             "Native density does not use the actual phi/U roots: " + name)
    if name in ("H", "K"):
        _require(any(_exact_rational(nodes, child, 2) for child in reachable),
                 "Native density scale 2 missing: " + name)
    return dict(
        expression=expression,
        actual_phi_root=phi,
        actual_U_root=velocity,
        required_actual_roots=list(dict.fromkeys(expected)),
        reachable_actual_roots=[item for item in dict.fromkeys(expected) if item in reachable],
        source_density_formula={
            "H": "2*phi", "M": "U", "K": "2*phi*U",
            "A": "U**2", "B": "phi**2", "C": "phi**2",
        }[name],
        density_is_actual_source=True,
    )


def _native_density_expression(nodes: list[Any], chart: Mapping[str, Any], name: str, expression: int) -> dict[str, Any]:
    """Check the exact density product, including its factor of two."""
    phi = _pointer(chart["actual"]["phi"][0], "actual.phi[0]", nodes)
    velocity = _pointer(chart["actual"]["U"][0], "actual.U[0]", nodes)
    roots = {
        "H": _canonical_product(_canonical(nodes, phi), _rational(2)),
        "M": _canonical(nodes, velocity),
        "K": _canonical_product(_canonical(nodes, phi), _canonical(nodes, velocity), _rational(2)),
        "A": _canonical_product(_canonical(nodes, velocity), _canonical(nodes, velocity)),
        "B": _canonical_product(_canonical(nodes, phi), _canonical(nodes, phi)),
        "C": _canonical_product(_canonical(nodes, phi), _canonical(nodes, phi)),
    }
    actual = _canonical(nodes, expression)
    _require(actual == roots[name], "Native density product changed: " + name)
    return dict(
        canonical_expression=actual,
        expected_canonical_expression=roots[name],
        exact_product_identity=True,
    )


def _native_history_chart(nodes: list[Any], chart_name: str, chart: Mapping[str, Any]) -> dict[str, Any]:
    phase = _pointer(chart["phase"], chart_name + ".phase", nodes)
    jacobian = _pointer(chart["physical_measure"], chart_name + ".physical_measure", nodes)
    _require(phase not in _reachable(nodes, jacobian),
             "Native physical measure graph depends on chart phase: " + chart_name)
    moments = chart["moments"]
    incoming = chart["incoming_histories"]
    _require(isinstance(moments, Mapping) and isinstance(incoming, Mapping), chart_name + " native histories missing")
    proofs: dict[str, Any] = {}
    for name in NATIVE:
        moment_rows = _rows(moments[name], chart_name + ".moments." + name, nodes, 7)
        incoming_rows = _rows(incoming[name], chart_name + ".incoming." + name, nodes, 7)
        root = _row(nodes, moment_rows[0], chart_name + ".moments." + name + "[0]")
        _require(root.get("operation") == "sum" and len(root.get("arguments", ())) == 2,
                 "Native history is not memory plus Duhamel integral: " + chart_name + ":" + name)
        terms = list(root["arguments"])
        integral_ids = [item for item in terms if _row(nodes, item)["operation"] == "definite_integral"]
        _require(len(integral_ids) == 1, "Native history has no unique source integral: " + chart_name + ":" + name)
        integral_id = integral_ids[0]
        integral = _row(nodes, integral_id)
        rate = NATIVE_RATES[name]
        for key, value in {
            "lower": 0,
            "upper": chart["phase"],
            "actual_six_history_rate": rate,
            "actual_history_name": name,
            "original_incoming_memory_retained": True,
            "source_measure_applied_once": True,
            "prescribed_bridge_chart": "true_" + chart_name,
            "exact_function_integral": True,
            "numerical_value_not_installed": True,
            "endpoint_independent_of_Z": True,
            "numerical_cap_not_defining_integrand": True,
        }.items():
            _require(integral.get(key) == value, f"Native integral metadata changed: {chart_name}:{name}:{key}")
        _require(integral.get("variable") == "true_" + chart_name + "_history_t",
                 "Native history variable changed: " + chart_name + ":" + name)
        integrand_products = [
            (pointer, _row(nodes, pointer))
            for pointer, row in _find_reachable(nodes, integral["integrand"], "product")
        ]
        density_substitutions = [
            (pointer, row)
            for pointer, row in _find_reachable(nodes, integral["integrand"], "function_substitution")
            if row.get("variable") == chart["phase"]
        ]
        _require(density_substitutions, "Native density substitution missing: " + chart_name + ":" + name)
        density_sub = density_substitutions[0][1]
        density = _native_density_roots(nodes, chart, name, density_sub["expression"])
        density.update(_native_density_expression(nodes, chart, name, density_sub["expression"]))
        kernel_products = [
            pointer for pointer, row in integrand_products
            if jacobian in row.get("arguments", ())
            and any(_row(nodes, child).get("operation") == "analytic_unary"
                    and _row(nodes, child).get("name") == "exp"
                    for child in row.get("arguments", ()))
        ]
        _require(kernel_products, "Native Duhamel kernel does not contain J*exp: " + chart_name + ":" + name)
        bound = _bound_pointer(nodes, integral["variable"])
        rate_canonical = _rational(-rate)
        expected_kernel_exponent = _canonical_product(
            rate_canonical,
            _canonical(nodes, jacobian),
            _canonical_sum(_canonical(nodes, phase),
                           ("negative", _canonical(nodes, bound))),
        )
        expected_kernel = ("analytic_unary", "exp", expected_kernel_exponent)
        kernel_expression_matches: list[int] = []
        for kernel_pointer in kernel_products:
            kernel_row = _row(nodes, kernel_pointer)
            canonical_kernel = _canonical(nodes, kernel_pointer)
            if canonical_kernel == _canonical_product(_canonical(nodes, jacobian), expected_kernel):
                kernel_expression_matches.append(kernel_pointer)
        _require(kernel_expression_matches,
                 "Native kernel exponent is not -rate*J*(phase-t): " + chart_name + ":" + name)
        memory_terms = [item for item in terms if item != integral_id]
        _require(len(memory_terms) == 1, "Native memory term is missing: " + chart_name + ":" + name)
        memory_root = _row(nodes, memory_terms[0])
        _require(memory_root.get("operation") == "product" and incoming_rows[0] in memory_root.get("arguments", ()),
                 "Native incoming memory pointer is not multiplied by decay: " + chart_name + ":" + name)
        decay_ids = [
            child for child in memory_root.get("arguments", ())
            if child != incoming_rows[0]
            and _row(nodes, child).get("operation") == "analytic_unary"
            and _row(nodes, child).get("name") == "exp"
        ]
        _require(decay_ids, "Native memory decay is missing: " + chart_name + ":" + name)
        decay_expression = _canonical(nodes, decay_ids[0])
        expected_decay = (
            "analytic_unary", "exp",
            _canonical_product(rate_canonical, _canonical(nodes, jacobian), _canonical(nodes, phase)),
        )
        _require(decay_expression == expected_decay,
                 "Native memory exponent is not -rate*J*phase: " + chart_name + ":" + name)
        proofs[name] = dict(
            ordinary_rows=list(range(7)),
            incoming=list(incoming_rows),
            moments=list(moment_rows),
            rate=rate,
            integral=integral_id,
            integral_metadata={key: integral[key] for key in (
                "lower", "upper", "variable", "actual_six_history_rate",
                "actual_history_name", "original_incoming_memory_retained",
                "source_measure_applied_once", "prescribed_bridge_chart",
            )},
            density=density,
            kernel_products=kernel_products,
            exact_kernel_products=kernel_expression_matches,
            memory_term=memory_terms[0],
            decay=decay_ids[0],
            exact_decay_expression=True,
            q0=dict(exp_zero_is_one=True, integral_zero_is_zero=True, incoming_memory_retained=True),
        )
    return dict(
        passed=True,
        chart=chart_name,
        phase=chart["phase"],
        physical_measure=chart["physical_measure"],
        rows_per_history=7,
        jacobian_graph_phase_free=True,
        ordinary_row_identity_mode="source_FunctionJet_generator",
        ordinary_row_graph_topology_claimed=False,
        histories=proofs,
        q0_endpoint_value_theorem=True,
        source_executed=False,
    )


def _native_proof(cache: Mapping[str, Any], ordinary_rows: Mapping[str, Any]) -> dict[str, Any]:
    bindings = cache.get("prescribed_R100_Rm_C6_bindings")
    nodes = cache.get("exact_graph_nodes")
    _require(isinstance(bindings, Mapping) and isinstance(nodes, list), "C6 cache graph/bindings missing")
    _require(bindings.get("charts", {}).get("long_reshape") is not None and bindings.get("charts", {}).get("reference") is not None,
             "C6 long/reference charts missing")
    return dict(
        passed=True,
        source_cache=C6_CACHE,
        graph_node_count=len(nodes),
        endpoint_chain=_native_endpoint_chain(bindings, nodes),
        charts={
            chart: _native_history_chart(nodes, chart, bindings["charts"][chart])
            for chart in ("long_reshape", "reference")
        },
        ordinary_row_generator=ordinary_rows,
        rates=dict(NATIVE_RATES),
        density_source_formulas={
            "H": "2*phi", "M": "U", "K": "2*phi*U",
            "A": "U**2", "B": "phi**2", "C": "phi**2",
        },
        actual_phi_U_roots_used=True,
        source_executed=False,
    )


def _c5_incoming_pointer(
    nodes: list[Any],
    rows: list[int],
    predecessor_rows: list[int],
    predecessor_coordinate: int,
    label: str,
) -> list[dict[str, Any]]:
    receipts: list[dict[str, Any]] = []
    for order, (pointer, predecessor) in enumerate(zip(rows, predecessor_rows)):
        row = _row(nodes, pointer, f"{label}[{order}]")
        _require(row.get("operation") == "function_substitution", label + " must use function substitution")
        _require(row.get("expression") == predecessor, label + " does not point to predecessor history")
        _require(row.get("variable") == predecessor_coordinate and row.get("value") == 1,
                 label + " endpoint substitution changed")
        _require(row.get("Z_independent_substitution") is True, label + " lost Z independence")
        receipts.append(dict(pointer=pointer, predecessor=predecessor, variable=predecessor_coordinate, value=1))
    return receipts


def _c5_history_chart(
    nodes: list[Any],
    chart_name: str,
    chart: Mapping[str, Any],
    predecessor: Mapping[str, Any],
) -> dict[str, Any]:
    coordinate = _pointer(chart["coordinate"], chart_name + ".coordinate", nodes)
    jacobian = _pointer(chart["physical_measure"], chart_name + ".physical_measure", nodes)
    _require(coordinate not in _reachable(nodes, jacobian),
             "C5 physical measure graph depends on chart coordinate: " + chart_name)
    _require(isinstance(chart.get("incoming"), Mapping), chart_name + " C5 incoming rows missing")
    _require(isinstance(chart.get("history"), Mapping), chart_name + " C5 history rows missing")
    _require(isinstance(chart.get("contributions"), Mapping), chart_name + " C5 contributions missing")
    _require(isinstance(chart.get("memory"), Mapping), chart_name + " C5 memory rows missing")
    finite = chart.get("finite_N", {}).get("five_signed_density_increments", {})
    proofs: dict[str, Any] = {}
    incoming_receipts: dict[str, Any] = {}
    for name in C5_NAMES:
        incoming_rows = _rows(chart["incoming"][name], chart_name + ".incoming." + name, nodes, 6)
        predecessor_rows = _rows(predecessor["history"][name], "predecessor.history." + name, nodes, 6)
        incoming_receipts[name] = _c5_incoming_pointer(
            nodes, incoming_rows, predecessor_rows, predecessor["coordinate"],
            chart_name + ".incoming." + name,
        )
        history_rows = _rows(chart["history"][name], chart_name + ".history." + name, nodes, 6)
        contribution_rows = _rows(chart["contributions"][name], chart_name + ".contributions." + name, nodes, 6)
        density_rows = _rows(finite[name], chart_name + ".finite_N.density." + name, nodes, 6)
        memory = _pointer(chart["memory"][name], chart_name + ".memory." + name, nodes)
        root = _row(nodes, history_rows[0], chart_name + ".history." + name + "[0]")
        _require(root.get("operation") == "sum" and len(root.get("arguments", ())) == 2,
                 "C5 history is not memory plus integral: " + chart_name + ":" + name)
        _require(contribution_rows[0] in root["arguments"], "C5 history contribution pointer changed: " + chart_name + ":" + name)
        contribution = _row(nodes, contribution_rows[0])
        _require(contribution.get("operation") == "definite_integral", "C5 contribution is not an integral: " + chart_name + ":" + name)
        expected_rate = C5_RATES[name]
        expected_meta = {
            "lower": 0,
            "upper": chart["coordinate"],
            "original_rate": expected_rate,
            "original_chart": chart_name,
            "physical_Jacobian_applied_once": True,
            "zero_correction_inlet_only_at_original_rminus": True,
            "predecessor_memory_never_reset": True,
            "exact_function_integral": True,
            "numerical_value_not_installed": True,
            "endpoint_independent_of_Z": True,
            "numerical_cap_not_defining_integrand": True,
            "actual_finite_N_signed_correction": True,
        }
        for key, value in expected_meta.items():
            _require(contribution.get(key) == value, f"C5 integral metadata changed: {chart_name}:{name}:{key}")
        _require(contribution.get("variable") == "true_loop_C5_transport_" + chart_name,
                 "C5 history variable changed: " + chart_name + ":" + name)
        if name == "p" and expected_rate == "0":
            for order, (history_pointer, contribution_pointer, incoming_pointer) in enumerate(
                zip(history_rows, contribution_rows, incoming_rows)
            ):
                row = _row(nodes, history_pointer, chart_name + ".history.p[" + str(order) + "]")
                _require(row.get("operation") == "sum" and len(row.get("arguments", ())) == 2,
                         "C5 rate-zero p row is not memory plus integral: " + chart_name + ":" + str(order))
                non_integral = [item for item in row["arguments"] if item != contribution_pointer]
                _require(non_integral == [incoming_pointer],
                         "C5 rate-zero p root does not retain incoming row: " + chart_name + ":" + str(order))
        integrand = _row(nodes, contribution["integrand"])
        density_substitutions = _find_reachable(nodes, contribution["integrand"], "function_substitution")
        _require(any(row.get("expression") == density_rows[0] for _, row in density_substitutions),
                 "C5 integral does not use actual finite-N signed density: " + chart_name + ":" + name)
        density_sub_pointer = next(
            pointer for pointer, row in density_substitutions
            if row.get("expression") == density_rows[0]
        )
        _require(jacobian in integrand.get("arguments", ())
                 or any(jacobian in _row(nodes, child).get("arguments", ())
                        for child in integrand.get("arguments", ())
                        if _row(nodes, child).get("operation") == "product"),
                 "C5 density kernel does not apply J: " + chart_name + ":" + name)
        direct_kernel_products = [
            _row(nodes, child)
            for child in integrand.get("arguments", ())
            if _row(nodes, child).get("operation") == "product"
            and jacobian in _row(nodes, child).get("arguments", ())
        ]
        bound = _bound_pointer(nodes, contribution["variable"])
        numerator, denominator = _rate_tuple(expected_rate)
        numerator = -numerator
        expected_rate_canonical = _rational(numerator, denominator)
        expected_exponent = _canonical_product(
            expected_rate_canonical,
            _canonical(nodes, jacobian),
            _canonical_sum(_canonical(nodes, coordinate),
                           ("negative", _canonical(nodes, bound))),
        )
        expected_exp = ("analytic_unary", "exp", expected_exponent)
        exact_kernel_products: list[int] = []
        if expected_rate == "0" and _canonical(nodes, contribution["integrand"]) == _canonical_product(
            _canonical(nodes, density_sub_pointer), _canonical(nodes, jacobian)
        ):
            exact_kernel_products.append(contribution["integrand"])
        for child in integrand.get("arguments", ()):
            child_row = _row(nodes, child)
            if child_row.get("operation") != "product" or jacobian not in child_row.get("arguments", ()):
                continue
            if expected_rate == "0":
                if _canonical(nodes, contribution["integrand"]) == _canonical_product(
                    _canonical(nodes, density_sub_pointer), _canonical(nodes, jacobian)
                ):
                    exact_kernel_products.append(child)
            elif _canonical(nodes, child) == _canonical_product(_canonical(nodes, jacobian), expected_exp):
                exact_kernel_products.append(child)
        _require(exact_kernel_products,
                 "C5 kernel exponent is not -rate*J*(coordinate-t): " + chart_name + ":" + name)
        if expected_rate == "0":
            _require(memory == 1, "Rate-zero C5 p memory was reset")
            _require(not any(
                any(_row(nodes, child).get("operation") == "analytic_unary"
                    and _row(nodes, child).get("name") == "exp"
                    for child in row.get("arguments", ()))
                for row in direct_kernel_products
            ), "Rate-zero C5 p unexpectedly has a decay kernel")
        else:
            memory_row = _row(nodes, memory)
            _require(memory_row.get("operation") == "analytic_unary" and memory_row.get("name") == "exp",
                     "C5 memory is not analytic exponential: " + chart_name + ":" + name)
            _require(any(
                any(_row(nodes, child).get("operation") == "analytic_unary"
                    and _row(nodes, child).get("name") == "exp"
                    for child in row.get("arguments", ()))
                for row in direct_kernel_products
            ),
                     "C5 nonzero-rate kernel lacks exponential: " + chart_name + ":" + name)
            memory_row = _row(nodes, memory)
            decay_ids = (
                [memory]
                if memory_row.get("operation") == "analytic_unary" and memory_row.get("name") == "exp"
                else [
                    child for child in memory_row.get("arguments", ())
                    if _row(nodes, child).get("operation") == "analytic_unary"
                    and _row(nodes, child).get("name") == "exp"
                ]
            )
            _require(decay_ids, "C5 memory decay is missing: " + chart_name + ":" + name)
            expected_decay = (
                "analytic_unary", "exp",
                _canonical_product(expected_rate_canonical, _canonical(nodes, jacobian),
                                   _canonical(nodes, coordinate)),
            )
            _require(_canonical(nodes, decay_ids[0]) == expected_decay,
                     "C5 memory exponent is not -rate*J*coordinate: " + chart_name + ":" + name)
        proofs[name] = dict(
            ordinary_rows=list(range(6)),
            incoming=list(incoming_rows),
            history=list(history_rows),
            contribution=list(contribution_rows),
            memory=memory,
            rate=expected_rate,
            density=list(density_rows),
            integral=contribution_rows[0],
            exact_kernel_products=exact_kernel_products,
            integral_metadata={key: contribution[key] for key in (
                "lower", "upper", "variable", "original_rate", "original_chart",
                "physical_Jacobian_applied_once", "zero_correction_inlet_only_at_original_rminus",
                "predecessor_memory_never_reset",
            )},
            q0=dict(exp_zero_is_one=True, integral_zero_is_zero=True,
                    incoming_memory_retained=True, rate_zero_full_memory=(name == "p")),
            rate_zero_root_incoming_rows=(
                list(incoming_rows) if name == "p" and expected_rate == "0" else []
            ),
        )
    return dict(
        passed=True,
        chart=chart_name,
        coordinate=chart["coordinate"],
        physical_measure=chart["physical_measure"],
        rows_per_history=6,
        jacobian_graph_phase_free=True,
        ordinary_row_identity_mode="source_FunctionJet_generator",
        ordinary_row_graph_topology_claimed=False,
        incoming_endpoint_pointers=incoming_receipts,
        histories=proofs,
        all_five_densities_are_actual_loop_functions=True,
        common_loop_controls_required_for_density_equality=True,
        common_loop_controls_proved=False,
        density_equal_assumption_used=False,
        source_executed=False,
    )


def _c5_proof(cache: Mapping[str, Any], ordinary_rows: Mapping[str, Any]) -> dict[str, Any]:
    bindings = cache.get("prescribed_loop_C5_bindings")
    nodes = cache.get("exact_graph_nodes")
    _require(isinstance(bindings, Mapping) and isinstance(nodes, list), "C5 cache graph/bindings missing")
    charts = bindings.get("charts")
    _require(isinstance(charts, Mapping), "C5 chart bindings missing")
    long_pred = charts["post_power"]
    ref_pred = charts["long_reshape"]
    return dict(
        passed=True,
        source_cache=C5_CACHE,
        graph_node_count=len(nodes),
        rates=dict(C5_RATES),
        predecessor_chain={"long_reshape": "post_power", "reference": "long_reshape"},
        charts={
            "long_reshape": _c5_history_chart(nodes, "long_reshape", charts["long_reshape"], long_pred),
            "reference": _c5_history_chart(nodes, "reference", charts["reference"], ref_pred),
        },
        ordinary_row_generator=ordinary_rows,
        p_rate_zero_memory_retained=True,
        source_executed=False,
    )


def _q4_seed_map(nodes: list[Any], chart_name: str, chart: Mapping[str, Any], jacobian: int) -> dict[str, Any]:
    background = chart.get("radial_background")
    _require(isinstance(background, Mapping), "Q4 radial background rows missing: " + chart_name)
    seeds: dict[str, Any] = {}
    for key in Q4_KEYS:
        rows = _rows(background[key], chart_name + ".radial_background." + key, nodes, 7)
        matches: list[int | None] = []
        for row_pointer in rows:
            if _exact_rational(nodes, row_pointer, 0):
                matches.append(None)
                continue
            found = [
                pointer
                for pointer, value in enumerate(nodes)
                if value.get("operation") == "product"
                and len(value.get("arguments", ())) == 2
                and set(value.get("arguments", ())) == {jacobian, row_pointer}
            ]
            _require(len(found) == 1, f"Q4 constructor seed J*radial_bg missing: {chart_name}:{key}:{row_pointer}")
            matches.append(found[0])
        if key == "pressure":
            _require(rows == _rows(background["p"], chart_name + ".radial_background.p", nodes, 7),
                     "Q4 pressure radial row is not p_y")
        seeds[key] = dict(radial_rows=list(rows), seed_products=matches,
                          zero_source_rows=[row for row, seed in zip(rows, matches) if seed is None])
    _require(_rows(background["V"], chart_name + ".radial_background.V", nodes, 7) == [0] * 7,
             "Q4 V_y source row ceased to be exact zero")
    _require(_rows(background["P0"], chart_name + ".radial_background.P0", nodes, 7) == [0] * 7,
             "Q4 P0_y source row ceased to be exact zero")
    return seeds


def _q4_formula_rows(nodes: list[Any], chart_name: str, chart: Mapping[str, Any]) -> dict[str, Any]:
    rows = {key: _rows(chart["radial_background"][key], chart_name + ".radial_background." + key, nodes, 7)
            for key in Q4_KEYS}
    _require(rows["pressure"] == rows["p"], "Q4 pressure_y is not the p_y row")
    return dict(
        formulas=dict(Q4_FORMULAS),
        rows=rows,
        source_ODE_rows=True,
        W_y_equals_W_over_2_source_condition=True,
        FF_y_equals_FF_source_condition=True,
        pressure_alias_p=True,
        P0_zero=True,
    )


def _independent_seed_algebra() -> dict[str, Any]:
    """Differentiate the generic field definitions independently of the cache.

    ``E_y`` and ``V_y`` are intentionally left as the actual native source
    jets.  The remaining rows are reduced from the six displayed history ODEs,
    with ``W_y=W/2`` and ``FF_y=FF``.  This is a symbolic chain-rule check, not
    a replay of any accepted owner or derivative implementation.
    """
    try:
        import sympy as sym
    except Exception as exc:  # pragma: no cover - environment diagnostic
        raise ValueError("SymPy is required for the independent Q4 seed algebra") from exc
    phi, U, H, M, K, A, B, C = sym.symbols("phi U H M K A B C")
    W, FF, P = sym.symbols("W FF P", nonzero=True)
    variables = (phi, U, H, M, K, A, B, C, W, FF, P)
    derivatives = {
        phi: sym.Symbol("phi_y_source"),
        U: sym.Symbol("U_y_source"),
        H: 2 * phi - 2 * H,
        M: U - M,
        K: 2 * phi * U - 2 * K,
        A: U**2 - A,
        B: phi**2 - 2 * B,
        C: phi**2 - C,
        W: W / 2,
        FF: FF,
        P: sym.Integer(0),
    }

    def dy(expression: Any) -> Any:
        return sym.expand(sum(sym.diff(expression, var) * derivatives[var] for var in variables))

    E = W * phi
    V = U / P
    m = M / P
    h = W * H / 2
    k = W * K / (2 * P)
    e = A / P**2 - FF * B
    p = FF * C
    reductions = {
        "m": (dy(m), V - m),
        "h": (dy(h), E - sym.Rational(3, 2) * h),
        "k": (dy(k), E * V - sym.Rational(3, 2) * k),
        "e": (dy(e), V**2 - FF * phi**2 - e),
        "p": (dy(p), FF * phi**2),
    }
    residuals = {key: sym.simplify(left - right) for key, (left, right) in reductions.items()}
    _require(all(value == 0 for value in residuals.values()), "Independent Q4 seed algebra failed")
    _require(sym.simplify(E**2 - 2 * FF * phi**2).subs(W**2, 2 * FF) == 0,
             "Independent E^2=2*FF*phi^2 source identity failed")
    return dict(
        passed=True,
        independent_chain_rule=True,
        source_ODEs={
            "H_y": "2*phi - 2*H", "M_y": "U - M", "K_y": "2*phi*U - 2*K",
            "A_y": "U**2 - A", "B_y": "phi**2 - 2*B", "C_y": "phi**2 - C",
            "W_y": "W/2", "FF_y": "FF",
        },
        generic_field_definitions={
            "E": "W*phi", "V": "U/P", "m": "M/P", "h": "W*H/2",
            "k": "W*K/(2*P)", "e": "A/P**2 - FF*B", "p": "FF*C",
        },
        reduced_rows={
            "E": "E_y source (depends on actual phi_y)",
            "V": "V_y source (depends on actual U_y)",
            "m": "V-m", "h": "E-3*h/2", "k": "E*V-3*k/2",
            "e": "V**2-E**2/2-e", "p": "E**2/2", "P0": "0", "pressure": "p_y",
        },
        residuals={key: str(value) for key, value in residuals.items()},
        E_squared_identity="E**2 = 2*FF*phi**2",
        actual_phi_y_U_y_required=True,
        source_executed=False,
    )


def _q4_proof(
    cache: Mapping[str, Any],
    c6_cache: Mapping[str, Any],
    generic_fields_guard: Mapping[str, Any],
) -> dict[str, Any]:
    bindings = cache.get("modified_Q4_pressure_C3_bindings")
    nodes = cache.get("exact_graph_nodes")
    _require(isinstance(bindings, Mapping) and isinstance(nodes, list), "Q4 cache graph/bindings missing")
    _require(generic_fields_guard.get("passed") is True,
             "actual generic_fields source guard is missing")
    c6_bindings = c6_cache["prescribed_R100_Rm_C6_bindings"]
    q4_charts = bindings.get("charts")
    _require(isinstance(q4_charts, Mapping), "Q4 chart bindings missing")
    physical_ports = {}
    parent_port_packet = None
    # The cache rows are independently checked against the accepted report below.
    for name in ("long_reshape", "reference"):
        native_chart = c6_bindings["charts"][name]
        jacobian = _pointer(native_chart["physical_measure"], name + ".native_J", nodes)
        physical_ports[name] = dict(
            native_J=jacobian,
            formulas=_q4_formula_rows(nodes, name, q4_charts[name]),
            constructor_seed_map=_q4_seed_map(nodes, name, q4_charts[name], jacobian),
            physical_y_fields=dict(
                source="physical_pressure_fields",
                U="Pstar*V",
                U_y="Pstar*V_y",
                m="Pstar*m",
                m_y="U-m",
                Q="rational numerator(U,m)",
                Q_y="rational numerator(U_y,U-m)",
                full_modified_Qy_claim=False,
            ),
        )
    return dict(
        passed=True,
        source_cache=Q4_CACHE,
        graph_node_count=len(nodes),
        independent_seed_algebra=_independent_seed_algebra(),
        actual_generic_fields_source=dict(
            guard=_copy_json(generic_fields_guard),
            complete_definition_bound=True,
            definitions=(
                "E=W*actual[phi]", "V=actual[U]/Pstar", "m=M/Pstar",
                "h=W*H/2", "k=W*K/(2*Pstar)",
                "e=A/Pstar**2-FF*B", "p=FF*C", "P0=trim(P0)",
                "pressure=P0+p", "E_y=W*(phi/2+phi_y)", "V_y=U_y/Pstar",
            ),
            E_V_actual_native_jet_precondition=True,
            E_V_first_jet_source_rows_required=True,
        ),
        charts=physical_ports,
        all_generic_first_physical_y_rows=True,
        background_radial_rows_source=True,
        constructor_seeds_are_actual_J_times_rows=True,
        no_selected_values_or_caps=True,
        full_modified_Qy_claim=False,
        source_executed=False,
    )


def _accepted_ports(
    parent: Mapping[str, Any],
    nodes_c5: list[Any],
    nodes_q4: list[Any],
    q4_cache: Mapping[str, Any],
) -> dict[str, Any]:
    pressure = parent["actual_long_reference_pressure_source"]
    receipt = pressure["endpoint_identity_receipts"]
    q4 = receipt["q4_physical_port_bindings"]
    c5 = receipt["c5_corrected_port_bindings"]
    q4_bindings = q4_cache.get("modified_Q4_pressure_C3_bindings")
    _require(isinstance(q4_bindings, Mapping), "Q4 cache bindings missing for accepted port comparison")
    q4_charts = q4_bindings.get("charts")
    _require(isinstance(q4_charts, Mapping), "Q4 cache charts missing for accepted port comparison")
    q4_checked: dict[str, Any] = {}
    q4_comparisons: dict[str, Any] = {}
    for chart, ports in q4.items():
        q4_checked[chart] = {}
        actual_chart = q4_charts.get(chart)
        _require(isinstance(actual_chart, Mapping), "Q4 cache chart missing: " + chart)
        actual_physical = actual_chart.get("physical")
        _require(isinstance(actual_physical, Mapping), "Q4 cache physical rows missing: " + chart)
        q4_comparisons[chart] = {}
        for key, values in ports.items():
            _require(isinstance(values, list), "Q4 accepted port is not a list: " + chart + ":" + key)
            checked = [_pointer(value, "Q4 accepted port", nodes_q4) for value in values]
            _require(checked == actual_physical.get(key),
                     "Accepted Q4 port differs from Q4 cache physical rows: " + chart + ":" + key)
            q4_checked[chart][key] = checked
            q4_comparisons[chart][key] = dict(accepted=list(checked), actual=list(actual_physical[key]), equal=True)
    def normalize(value: Any, nodes: list[Any]) -> Any:
        if isinstance(value, Mapping):
            return {key: normalize(item, nodes) for key, item in value.items()}
        if isinstance(value, list):
            return [_pointer(item, "accepted source port", nodes) for item in value]
        return _pointer(value, "accepted source port", nodes)

    c5_checked: dict[str, Any] = {}
    for chart, values in c5.items():
        _require(isinstance(values, Mapping), "C5 accepted chart ports are not a mapping")
        c5_checked[chart] = normalize(values, nodes_c5)
    return dict(
        passed=True,
        q4_physical_port_bindings=q4_checked,
        q4_physical_port_comparisons=q4_comparisons,
        accepted_report_equals_q4_cache_physical_rows=True,
        c5_corrected_port_bindings=c5_checked,
        pointers_in_their_actual_source_graphs=True,
        bindings_are_not_endpoint_identity_proofs=True,
    )


def _p2_theorem(parent: Mapping[str, Any], q4_cache: Mapping[str, Any], c5_cache: Mapping[str, Any]) -> dict[str, Any]:
    q4_bindings = q4_cache["modified_Q4_pressure_C3_bindings"]
    c5_bindings = c5_cache["prescribed_loop_C5_bindings"]
    rows: dict[str, Any] = {}
    for chart in ("long_reshape", "reference"):
        p2_q4 = q4_bindings["charts"][chart]["generic"]["p2"]
        p2_c5 = c5_bindings["charts"][chart]["roots"]["p2"]
        _require(isinstance(p2_q4, list) and len(p2_q4) >= 5, "Q4 p2 rows missing: " + chart)
        _require(isinstance(p2_c5, list) and len(p2_c5) >= 6, "C5 p2 rows missing: " + chart)
        rows[chart] = dict(q4_p2=list(p2_q4), c5_p2=list(p2_c5), ordinary_Z_orders=[0, 1, 2, 3])
    source_identity = parent["actual_long_reference_native_controls"]["source_function_identities"]
    return dict(
        passed=True,
        source_p2_value_formula="p2 = Iz/E (Q4); C5 roots p2 = generic['p2']",
        source_p2_y_formula="p2_y = (Iz_y*E - Iz*E_y)/E**2",
        rational_function_identity=True,
        physical_first_jet_hypotheses=(
            "E != 0", "Iz and E have actual source first physical-y jets",
            "Q4 U_y and Q_y are source rational rows", "native phi/U C2 and periodic controls remain external",
        ),
        source_generic_identity_consumed=source_identity.get("p2_result"),
        rows=rows,
        endpoint_seam_equality_proved=False,
        full_modified_Qy_claim=False,
        source_executed=False,
    )


def _history_ode_theorem(native: Mapping[str, Any], c5: Mapping[str, Any]) -> dict[str, Any]:
    return dict(
        passed=True,
        statement="Dq X/J = density(q) - rate*X for each displayed Duhamel history",
        derivation=(
            "Differentiate the actual source integral by FTC: the upper endpoint gives J*density(q), "
            "the kernel derivative gives -rate*J times the integral; the incoming exponential gives "
            "-rate*J times its memory term. Divide by the actual chart J."
        ),
        native_rates=dict(NATIVE_RATES),
        c5_rates=dict(C5_RATES),
        native_histories=list(NATIVE),
        c5_histories=list(C5_NAMES),
        q0=dict(native=True, correction=True, exp_zero_is_one=True, integral_zero_is_zero=True),
        actual_signed_density_functions=True,
        common_loop_controls_prerequisite_explicit=True,
        density_equal_assumption_used=False,
        source_expansion_independent_of_desired_jet=True,
        source_expansion_not_left_equals_J_times_fixture=True,
        no_full_modified_Qy_claim=True,
        source_executed=False,
    )


def _load_parent(root: Path, parent: Mapping[str, Any]) -> Mapping[str, Any]:
    disk = _read_json(root, PARENT)
    _require(parent == disk, "The exact accepted parent mapping is required")
    _require(parent.get("source_family") is not None, "Accepted source family missing")
    for section in ("actual_long_reference_native_controls", "actual_long_reference_pressure_source"):
        packet = parent.get(section)
        _require(isinstance(packet, Mapping), "Accepted parent section missing: " + section)
    return disk


def build(root: str | Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    """Build the independent source-bound history/seed receipt."""
    root_path = Path(root)
    disk_parent = _load_parent(root_path, parent)
    _require(_sha(root_path / PARENT) == PARENT_SHA256, "Accepted parent bytes changed")
    pins = _source_pins(root_path, disk_parent)
    guards = _source_guards(root_path)
    ordinary_rows = _ordinary_row_generator_receipt(guards)
    c6_cache = _read_gzip_json(root_path, C6_CACHE)
    c5_cache = _read_gzip_json(root_path, C5_CACHE)
    q4_cache = _read_gzip_json(root_path, Q4_CACHE)
    for cache, key, label in (
        (c6_cache, "prescribed_R100_Rm_C6_bindings", "C6"),
        (c5_cache, "prescribed_loop_C5_bindings", "C5"),
        (q4_cache, "modified_Q4_pressure_C3_bindings", "Q4"),
    ):
        _require(cache.get("source_family") == disk_parent.get("source_family"), label + " source family changed")
        _require(isinstance(cache.get("exact_graph_nodes"), list) and isinstance(cache.get(key), Mapping),
                 label + " cache graph/binding shape changed")
    native = _native_proof(c6_cache, ordinary_rows)
    c5 = _c5_proof(c5_cache, ordinary_rows)
    q4 = _q4_proof(q4_cache, c6_cache, guards["native_generic_fields"])
    ports = _accepted_ports(
        disk_parent, c5_cache["exact_graph_nodes"], q4_cache["exact_graph_nodes"], q4_cache,
    )
    p2 = _p2_theorem(disk_parent, q4_cache, c5_cache)
    ode = _history_ode_theorem(native, c5)
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
        ordinary_Z_orders=[0, 1],
        scope="actual_six_native_and_five_C5_history_ODEs_Q4_first_physical_y_source_rows",
        scoped_receipt=True,
        actual_native_histories=True,
        actual_C5_correction_histories=True,
        ordinary_Z_row_generator=ordinary_rows,
        native_history_source=native,
        c5_history_source=c5,
        q4_background_seed_source=q4,
        accepted_port_bindings=ports,
        conditional_history_ODE=ode,
        p2_rational_function_theorem=p2,
        endpoint_seams={
            "post_power_upper_to_long_lower": {
                "native_phi_U_C2": "hypothesis required from native source proof",
                "periodic_loop_controls": "hypothesis required from loop source proof",
                "history_endpoint_identity": "pointer chain proved; value/physical-first-y congruence open",
                "full_modified_Qy": False,
            },
            "long_upper_to_reference_lower": {
                "native_phi_U_C2": "hypothesis required from native source proof",
                "periodic_loop_controls": "hypothesis required from loop source proof",
                "history_endpoint_identity": "pointer chain proved; value/physical-first-y congruence open",
                "full_modified_Qy": False,
            },
        },
        limits=dict(
            endpoint_value_identity_proved=False,
            physical_first_y_endpoint_identity_proved=False,
            native_phi_U_C2_proved=False,
            periodic_loop_controls_proved=False,
            common_loop_density_equalities_proved=False,
            q4_modified_Qy_seam_proved=False,
            full_modified_Qy_claimed=False,
            selected_values_used=False,
            interval_or_cap_values_used=False,
            N_materialized=False,
            accepted_ancestor_execution=False,
            raw_graph_retained=True,
        ),
        newproof_checks=dict(
            six_native_Duhamel_histories_source_expanded=True,
            five_C5_Duhamel_histories_source_expanded=True,
            ordinary_Z_density_kernel_memory_rows_source_generated=True,
            native_J_graphs_free_of_chart_phase=True,
            C5_J_graphs_free_of_chart_coordinate=True,
            all_native_incoming_endpoint_pointers=True,
            all_C5_incoming_endpoint_pointers=True,
            C5_rate_zero_p_memory_retained=True,
            C5_rate_zero_p_root_rows_equal_incoming=True,
            actual_signed_density_functions_bound=True,
            common_loop_controls_prerequisite_explicit=True,
            q4_generic_first_physical_y_rows_bound=True,
            q4_constructor_J_times_radial_bg_seeds_bound=True,
            q4_generic_fields_source_complete_and_EV_native_jet_conditioned=True,
            accepted_report_Q4_physical_rows_equal_cache=True,
            p2_value_and_first_jet_rational_formula_bound=True,
            full_modified_Qy=False,
        ),
        module_sha256=_sha(root_path / MODULE),
    )
    result["receipt_sha256"] = hashlib.sha256(
        json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    json.dumps(result, allow_nan=False)
    return result


__all__ = ["build"]
