"""Bounded hash-consed proof for the issued Q4/C5 physical-y roots.

The saved Q4 packet contains the C5 and C6 graph prefixes.  This module reads
those packets and the source files as data only; it never imports or executes
an ancestor producer, checker, or ``NativeDerivative``.  Expressions are
canonicalized into a small hash-consed DAG.  A derivative receipt is therefore
bounded by the reachable graph and does not materialize the large nested tuple
expressions used by the superseded draft.

The field equations ``Dq(field) = J * background_radial_row`` are emitted as
explicit upstream obligations.  They are intentionally not discharged from
arbitrary flags in a parent report.  Consequently a successful structural
receipt remains conditional and always reports both physical identity gates as
false.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
import time
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Mapping


Q4_STEM = "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions"
C5_STEM = "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions"
C6_STEM = "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions"
PARENT_NAME = "lei_ren_part1_paper_compliant_current_original_long_reference_seam_prerequisites.json"
SCOPE_CHARTS = ("post_power", "long_reshape", "reference")
C6_CHART_FOR_Q4 = {
    "post_power": "post_switch",
    "long_reshape": "long_reshape",
    "reference": "reference",
}
FIELD_NAMES = ("E", "V", "m", "h", "k", "e", "p", "P0")
ORDERS = tuple(range(6))
EXPECTED_PARENT_SHA256 = "0a91144aca1caaf0dd2b11170b1aee4339f676eb5f213eea236ae8c814cea0c4"
EXPECTED_INPUT_SHA256 = {
    f"experiments/root_st073/{C5_STEM}.py": "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    f"experiments/root_st073/{C5_STEM}.json.gz": "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
    f"experiments/root_st073/{C5_STEM}_check.json": "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7",
    f"experiments/root_st073/{C6_STEM}.py": "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    f"experiments/root_st073/{C6_STEM}.json.gz": "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a",
    f"experiments/root_st073/{C6_STEM}_check.json": "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6",
    f"experiments/root_st073/{Q4_STEM}.py": "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
    f"experiments/root_st073/{Q4_STEM}.json.gz": "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a",
    f"experiments/root_st073/{Q4_STEM}_check.json": "189ef4dacc60bc2feda16046a193a75e89f658f635fd6e3e381200c09637bf6f",
}


class GraphProofError(RuntimeError):
    """Malformed saved input or an invalid proof request."""


class OpenRuleError(GraphProofError):
    """A reachable derivative rule is intentionally left open."""


class BudgetExceeded(GraphProofError):
    """A single root check exceeded its bounded time budget."""


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _locate(root: str | Path, name: str) -> Path:
    base = Path(root)
    candidates = (base / name, base / "experiments" / "root_st073" / name)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    if base.is_file() and base.name == name:
        return base
    # Focused review calls sometimes use the experiments directory as root.
    for candidate in base.rglob(name):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(name)


def _read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _read_gzip_json(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise GraphProofError(f"saved packet is not an object: {path}")
    return value


def _canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def _mapping_digest(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _source_method(tree: ast.AST, class_name: str, method_name: str) -> ast.FunctionDef:
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for child in node.body:
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == method_name:
                    return child  # type: ignore[return-value]
    raise GraphProofError(f"missing {class_name}.{method_name}")


def _compact_source(node: ast.AST) -> str:
    return "".join(ast.unparse(node).split())


def _native_rule_binding(source: Path) -> dict[str, Any]:
    """Hash and inspect the issued NativeDerivative AST without executing it."""

    tree = ast.parse(source.read_text(encoding="utf-8"), filename=str(source))
    methods = {
        name: _source_method(tree, "NativeDerivative", name)
        for name in ("__init__", "seed", "sub", "sigma_prime", "derivative")
    }
    compact = {name: _compact_source(node) for name, node in methods.items()}
    ast_hashes = {
        name: hashlib.sha256(ast.dump(node, include_attributes=False).encode("utf-8")).hexdigest()
        for name, node in methods.items()
    }
    needles: dict[str, tuple[str, ...]] = {
        "constructor_phase_cache": ("self.cache={phase.node:g.zero}", "self.core_extensions=[]"),
        "seed_cache_assignment": ("self.cache[left.node]=right",),
        "substitution_builder": (
            "function_substitution",
            "Z_independent_substitution=True",
            "exact_native_chain_rule_substitution=True",
        ),
        "sigma_prime_active_body": (
            "g.mul(sig,g.sub(g.one,sig),first)",
            "g.quotient(c(2),g.mul(x,x,x)",
            "g.quotient(c(2),g.mul(comp,comp,comp)",
            "Delta=x.node",
            "eta=g.one.node",
            "flat_value=g.zero.node",
            "active_body_not_evaluated_on_flat_branch=True",
            "original_sigma_exact_logistic_first_derivative=True",
        ),
        "function_substitution_chain": (
            "self.sub(self.derivative(body),var,arg)",
            "self.sub(self.other(var).derivative(body),var,arg)",
            "self.derivative(arg)",
        ),
        "moving_integral_partial": (
            "self.derivative(body)",
            "exact_native_partial_of_same_integrand",
            "self.sub(body,var,upper)",
            "self.sub(body,var,lower)",
        ),
        "inverse_angle_chain": (
            "substitute_original_inverse_angle",
            "self.other(angle).derivative(body)",
            "self.derivative(inverse)",
        ),
        "implicit_inverse_rule": (
            "self.derivative(target)",
            "at(held)",
            "lam",
            "exactimplicitnativederivative",
        ),
        "lazy_derivative": ("same_original_lazy_native_derivative",),
        "analytic_core_rule": ("current_original_unique_analytic_core_component",),
    }
    method_for_rule = {
        "constructor_phase_cache": "__init__",
        "seed_cache_assignment": "seed",
        "substitution_builder": "sub",
        "sigma_prime_active_body": "sigma_prime",
    }
    checks: dict[str, bool] = {}
    snippets: list[tuple[str, str]] = []
    for name, required in needles.items():
        method = method_for_rule.get(name, "derivative")
        checks[name] = all(token in compact[method] for token in required)
        snippets.extend((f"NativeDerivative.{method}", token) for token in required)
    return {
        "source": source.name,
        "class": "NativeDerivative",
        "methods": {
            name: {"lineno": node.lineno, "AST_sha256": ast_hashes[name]}
            for name, node in methods.items()
        },
        "rule_checks": checks,
        "all_rule_bindings_present": all(checks.values()),
        "source_guard_snippets": snippets,
    }


class _DAG:
    """A compact, deterministic, hash-consed symbolic expression DAG."""

    def __init__(self) -> None:
        self.nodes: list[tuple[str, tuple[int, ...], tuple[Any, ...]]] = []
        self.keys: dict[tuple[str, tuple[int, ...], tuple[Any, ...]], int] = {}
        self.digests: list[str] = []
        self._substitution_cache: dict[tuple[int, str, int], int] = {}
        self.zero = self.rat(0)
        self.one = self.rat(1)

    def _intern(self, op: str, children: Iterable[int] = (), atoms: Iterable[Any] = ()) -> int:
        child_tuple = tuple(int(child) for child in children)
        atom_tuple = tuple(atoms)
        key = (op, child_tuple, atom_tuple)
        old = self.keys.get(key)
        if old is not None:
            return old
        index = len(self.nodes)
        self.keys[key] = index
        self.nodes.append(key)
        digest = hashlib.sha256(
            _canonical_json({
                "op": op,
                "children": [self.digests[child] for child in child_tuple],
                "atoms": atom_tuple,
            }).encode("utf-8")
        ).hexdigest()
        self.digests.append(digest)
        return index

    def node(self, index: int) -> tuple[str, tuple[int, ...], tuple[Any, ...]]:
        return self.nodes[index]

    def digest(self, index: int) -> str:
        return self.digests[index]

    def rat(self, value: Fraction | int) -> int:
        value = Fraction(value)
        return self._intern("rat", atoms=(value.numerator, value.denominator))

    def var(self, name: str) -> int:
        return self._intern("var", atoms=(str(name),))

    def source(self, operation: str, graph_index: int, raw_digest: str | None = None) -> int:
        # The graph index is retained: distinct source functions are never
        # identified merely because their metadata happens to match.  Retain
        # the raw row digest as a typed source binding as well.
        return self._intern(
            "source_leaf",
            atoms=(str(operation), int(graph_index), str(raw_digest) if raw_digest is not None else ""),
        )

    def is_rat(self, index: int) -> bool:
        return self.nodes[index][0] == "rat"

    def rat_value(self, index: int) -> Fraction | None:
        op, _children, atoms = self.nodes[index]
        if op != "rat":
            return None
        return Fraction(int(atoms[0]), int(atoms[1]))

    def add(self, *values: int) -> int:
        flat: list[int] = []
        rational = Fraction(0)
        for value in values:
            op = self.nodes[value][0]
            if op == "add":
                flat.extend(self.nodes[value][1])
            else:
                rat = self.rat_value(value)
                if rat is None:
                    flat.append(value)
                else:
                    rational += rat
        if rational:
            flat.append(self.rat(rational))
        flat = [value for value in flat if value != self.zero]
        if not flat:
            return self.zero
        if len(flat) == 1:
            return flat[0]
        flat.sort(key=self.digest)
        return self._intern("add", flat)

    def mul(self, *values: int) -> int:
        flat: list[int] = []
        rational = Fraction(1)
        for value in values:
            if value == self.zero:
                return self.zero
            op = self.nodes[value][0]
            if op == "mul":
                flat.extend(self.nodes[value][1])
            else:
                rat = self.rat_value(value)
                if rat is None:
                    flat.append(value)
                else:
                    rational *= rat
        if rational == 0:
            return self.zero
        if rational != 1:
            flat.append(self.rat(rational))
        flat = [value for value in flat if value != self.one]
        if not flat:
            return self.one
        if len(flat) == 1:
            return flat[0]
        flat.sort(key=self.digest)
        return self._intern("mul", flat)

    def neg(self, value: int) -> int:
        rat = self.rat_value(value)
        if rat is not None:
            return self.rat(-rat)
        if self.nodes[value][0] == "neg":
            return self.nodes[value][1][0]
        return self._intern("neg", (value,))

    def quotient(self, numerator: int, denominator: int) -> int:
        nrat = self.rat_value(numerator)
        drat = self.rat_value(denominator)
        if nrat is not None and drat is not None:
            if drat == 0:
                raise GraphProofError("zero rational denominator in exact DAG")
            return self.rat(nrat / drat)
        if numerator == self.zero:
            return self.zero
        if denominator == self.one:
            return numerator
        return self._intern("quot", (numerator, denominator))

    def unary(self, name: str, argument: int) -> int:
        return self._intern("unary", (argument,), (str(name),))

    def substitute(self, body: int, variable: str, value: int) -> int:
        return self._intern("substitute", (body, value), (str(variable),))

    def integral(self, body: int, variable: str, lower: int, upper: int) -> int:
        return self._intern("integral", (body, lower, upper), (str(variable),))

    def at_inverse(self, body: int, inverse: int, angle: str) -> int:
        return self._intern("at_inverse", (body, inverse), (str(angle),))

    def phase_inverse(self, phase: int, angle: str, target: int) -> int:
        return self._intern("phase_inverse", (phase, target), (str(angle),))

    def lazy(self, active: int, delta: int, eta: int, flat: int, predicate: int) -> int:
        return self._intern("lazy", (active, delta, eta, flat, predicate))

    def comparison(self, operator: str, left: int, right: int) -> int:
        return self._intern("comparison", (left, right), (str(operator),))

    def logical_or(self, *arguments: int) -> int:
        args = sorted(arguments, key=self.digest)
        return self._intern("or", args)

    def substitute_tree(self, body: int, variable: str, value: int) -> int:
        key = (body, str(variable), value)
        if key in self._substitution_cache:
            return self._substitution_cache[key]
        op, children, atoms = self.nodes[body]
        if op == "var":
            result = value if atoms[0] == variable else body
        elif op == "rat" or op == "source_leaf":
            result = body
        elif op == "integral":
            bound = str(atoms[0])
            if bound == variable:
                result = self.integral(children[0], bound, children[1], children[2])
            else:
                result = self.integral(
                    self.substitute_tree(children[0], variable, value),
                    bound,
                    self.substitute_tree(children[1], variable, value),
                    self.substitute_tree(children[2], variable, value),
                )
        elif op == "substitute":
            bound = str(atoms[0])
            if bound == variable:
                result = self.substitute(children[0], bound, self.substitute_tree(children[1], variable, value))
            else:
                result = self.substitute(
                    self.substitute_tree(children[0], variable, value),
                    bound,
                    self.substitute_tree(children[1], variable, value),
                )
        else:
            replaced = tuple(self.substitute_tree(child, variable, value) for child in children)
            if op == "add":
                result = self.add(*replaced)
            elif op == "mul":
                result = self.mul(*replaced)
            elif op == "neg":
                result = self.neg(replaced[0])
            elif op == "quot":
                result = self.quotient(*replaced)
            elif op == "unary":
                result = self.unary(str(atoms[0]), replaced[0])
            elif op == "integral":
                result = self.integral(replaced[0], str(atoms[0]), replaced[1], replaced[2])
            elif op == "at_inverse":
                result = self.at_inverse(replaced[0], replaced[1], str(atoms[0]))
            elif op == "phase_inverse":
                result = self.phase_inverse(replaced[0], str(atoms[0]), replaced[1])
            elif op == "lazy":
                result = self.lazy(*replaced)
            elif op == "comparison":
                result = self.comparison(str(atoms[0]), replaced[0], replaced[1])
            elif op == "or":
                result = self.logical_or(*replaced)
            else:
                result = self._intern(op, replaced, atoms)
        self._substitution_cache[key] = result
        return result


GRAPH_FIELDS: dict[str, tuple[str, ...]] = {
    "sum": ("arguments",),
    "product": ("arguments",),
    "negative": ("argument",),
    "positive_quotient": ("numerator", "denominator"),
    "analytic_unary": ("argument",),
    "exact_original_source_flat_sigma": ("argument",),
    "function_substitution": ("expression", "variable", "value"),
    "definite_integral": ("integrand", "lower", "upper"),
    "substitute_original_inverse_angle": ("body", "inverse_angle"),
    "original_monotone_phase_inverse": ("phase_function", "target_phase"),
    "original_lazy_flat_branch": ("active_body", "flat_value", "flat_predicate", "Delta", "eta"),
    "current_original_unique_analytic_core_component": ("rho", "Z_variable"),
    "current_original_leading_function_recipe": ("Z_variable",),
    "current_original_axis_pressure_function_derivative": ("Z_variable",),
    "exact_real_comparison": ("left", "right"),
    "logical_or": ("arguments",),
}


def _reachable(nodes: list[dict[str, Any]], roots: Iterable[int]) -> set[int]:
    seen: set[int] = set()
    stack = [int(root) for root in roots]
    while stack:
        index = stack.pop()
        if index in seen:
            continue
        if index < 0 or index >= len(nodes):
            raise GraphProofError(f"invalid reachable node {index!r}")
        seen.add(index)
        row = nodes[index]
        for field in GRAPH_FIELDS.get(str(row.get("operation")), ()):
            values = row.get(field, ())
            if not isinstance(values, list):
                values = (values,)
            stack.extend(value for value in values if isinstance(value, int) and not isinstance(value, bool))
    return seen


def _operation_counts(nodes: list[dict[str, Any]], indices: Iterable[int]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for index in indices:
        operation = str(nodes[index].get("operation"))
        counts[operation] = counts.get(operation, 0) + 1
    return dict(sorted(counts.items()))


def _first_dag_mismatch(dag: _DAG, expected: int, actual: int) -> dict[str, Any] | None:
    """Return the first structural local difference between two DAG roots."""

    seen: set[tuple[int, int]] = set()

    def visit(left: int, right: int) -> dict[str, Any] | None:
        pair = (int(left), int(right))
        if pair in seen or left == right:
            return None
        seen.add(pair)
        left_node = dag.node(left)
        right_node = dag.node(right)
        if left_node[0] != right_node[0] or left_node[2] != right_node[2]:
            return {
                "expected_node": int(left),
                "actual_node": int(right),
                "expected_digest": dag.digest(left),
                "actual_digest": dag.digest(right),
                "expected_operation": left_node[0],
                "actual_operation": right_node[0],
                "expected_atoms": list(left_node[2]),
                "actual_atoms": list(right_node[2]),
            }
        if len(left_node[1]) != len(right_node[1]):
            return {
                "expected_node": int(left),
                "actual_node": int(right),
                "expected_digest": dag.digest(left),
                "actual_digest": dag.digest(right),
                "expected_operation": left_node[0],
                "actual_operation": right_node[0],
                "expected_child_count": len(left_node[1]),
                "actual_child_count": len(right_node[1]),
            }
        for child_left, child_right in zip(left_node[1], right_node[1]):
            mismatch = visit(child_left, child_right)
            if mismatch is not None:
                return mismatch
        return None

    return visit(int(expected), int(actual))


class _Oracle:
    """Canonicalizer and independent derivative rules over the saved graph."""

    EMPTY_OPERATIONS = {
        "exact_rational",
        "exact_positive_integer_expression",
        "mathematical_pi",
        "mathematical_positive_infinity",
    }
    SUPPORTED_ANALYTIC = {
        "exp",
        "expm1",
        "log",
        "log1p",
        "sqrt",
        "positive_sqrt",
        "sin",
        "cos",
        "original_flat_sigma",
        "original_flat_sigma_prime",
    }

    def __init__(self, nodes: list[dict[str, Any]], dag: _DAG) -> None:
        self.nodes = nodes
        self.dag = dag
        self._canonical_cache: dict[int, int] = {}
        self._free_cache: dict[int, frozenset[str]] = {}
        self._derivative_cache: dict[tuple[int, tuple[Any, ...], int | None, str], int] = {}
        self._core_rows: dict[tuple[int, int, int, int], list[int]] = {}
        self.open_rules: set[str] = set()
        self.lazy_metadata_issues: list[str] = []
        self.lazy_lemma_receipts: list[dict[str, Any]] = []
        self.binder_scope_issues: list[str] = []
        self.binder_scope_receipts: list[dict[str, Any]] = []
        self._raw_digest_cache: dict[int, str] = {}
        self._ticks = 0
        self.deadline: float | None = None
        for index, row in enumerate(nodes):
            if row.get("operation") == "current_original_unique_analytic_core_component":
                key = (
                    int(row.get("rho", -1)),
                    int(row.get("Z_variable", -1)),
                    int(row.get("radial_order", -1)),
                    int(row.get("ordinary_slow_Z_derivative_order", -1)),
                )
                self._core_rows.setdefault(key, []).append(index)

    def _tick(self) -> None:
        self._ticks += 1
        if self.deadline is not None and (self._ticks & 1023) == 0 and time.monotonic() > self.deadline:
            raise BudgetExceeded("DAG derivative root exceeded time budget")

    def row(self, index: int) -> dict[str, Any]:
        if not isinstance(index, int) or isinstance(index, bool) or index < 0 or index >= len(self.nodes):
            raise GraphProofError(f"invalid graph node {index!r}")
        return self.nodes[index]

    def name(self, index: int) -> str:
        row = self.row(index)
        if row.get("operation") == "bound_variable":
            return str(row.get("name"))
        return str(index)

    def raw_digest(self, index: int) -> str:
        old = self._raw_digest_cache.get(index)
        if old is not None:
            return old
        value = hashlib.sha256(_canonical_json(self.row(index)).encode("utf-8")).hexdigest()
        self._raw_digest_cache[index] = value
        return value

    def seed_profile(self, seed_map: Mapping[int, int]) -> str:
        """Return the typed digest for one native seed environment."""

        payload = [
            (int(index), self.dag.digest(int(value)))
            for index, value in sorted(seed_map.items(), key=lambda item: int(item[0]))
        ]
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    def _sigma_source_row_guard(self, index: int, row: Mapping[str, Any]) -> dict[str, Any]:
        """Check the issued flat predicate for the five sigma rows in scope."""

        issues: list[str] = []

        def exact_value(handle: Any) -> Fraction | None:
            if not isinstance(handle, int) or isinstance(handle, bool):
                return None
            try:
                candidate = self.row(handle)
            except GraphProofError:
                return None
            if candidate.get("operation") != "exact_rational":
                return None
            numerator = candidate.get("numerator")
            denominator = candidate.get("denominator")
            if (
                not isinstance(numerator, int)
                or isinstance(numerator, bool)
                or not isinstance(denominator, int)
                or isinstance(denominator, bool)
                or denominator == 0
            ):
                return None
            return Fraction(numerator, denominator)

        eta_value = exact_value(row.get("eta"))
        flat_value = exact_value(row.get("flat_value"))
        if eta_value != Fraction(1):
            issues.append(f"node {index}: sigma eta is not exact one")
        if flat_value != Fraction(0):
            issues.append(f"node {index}: sigma flat_value is not exact zero")

        delta = row.get("Delta")
        predicate = row.get("flat_predicate")
        observed: list[tuple[str, Fraction | None, bool]] = []
        if not isinstance(delta, int) or isinstance(delta, bool):
            issues.append(f"node {index}: sigma Delta is not an integer handle")
        if not isinstance(predicate, int) or isinstance(predicate, bool):
            issues.append(f"node {index}: sigma flat_predicate is not an integer handle")
        else:
            try:
                predicate_row = self.row(predicate)
            except GraphProofError:
                predicate_row = {}
                issues.append(f"node {index}: sigma flat_predicate handle is invalid")
            arguments = predicate_row.get("arguments")
            if predicate_row.get("operation") != "logical_or" or not isinstance(arguments, list) or len(arguments) != 2:
                issues.append(f"node {index}: sigma flat_predicate is not a two-way logical_or")
            else:
                for child in arguments:
                    if not isinstance(child, int) or isinstance(child, bool):
                        issues.append(f"node {index}: sigma predicate child is not an integer handle")
                        continue
                    try:
                        comparison = self.row(child)
                    except GraphProofError:
                        comparison = {}
                        issues.append(f"node {index}: sigma predicate child handle is invalid")
                    operator = comparison.get("operator")
                    left = comparison.get("left")
                    right = comparison.get("right")
                    right_value = exact_value(right)
                    left_matches = (
                        isinstance(delta, int)
                        and not isinstance(delta, bool)
                        and isinstance(left, int)
                        and not isinstance(left, bool)
                        and left == delta
                    )
                    if comparison.get("operation") != "exact_real_comparison":
                        issues.append(f"node {index}: sigma predicate child is not exact_real_comparison")
                    if not left_matches:
                        issues.append(f"node {index}: sigma predicate comparison does not use Delta handle")
                    observed.append((str(operator), right_value, left_matches))
        expected = {("<=", Fraction(0)), (">=", Fraction(1))}
        actual = {(operator, value) for operator, value, left_matches in observed if left_matches}
        if actual != expected or len(observed) != 2:
            issues.append(f"node {index}: sigma predicate is not exactly Delta<=0 or Delta>=1")
        predicate_operation = None
        if isinstance(predicate, int) and not isinstance(predicate, bool) and 0 <= predicate < len(self.nodes):
            predicate_operation = self.row(predicate).get("operation")
        return {
            "applied": True,
            "passed": not issues,
            "eta_exact_value": str(eta_value) if eta_value is not None else None,
            "flat_value_exact_value": str(flat_value) if flat_value is not None else None,
            "predicate_operation": predicate_operation,
            "predicate_comparisons": [
                {"operator": operator, "right_exact_value": str(value) if value is not None else None, "left_is_Delta": left_matches}
                for operator, value, left_matches in observed
            ],
            "active_body_not_evaluated_on_flat_branch": row.get("active_body_not_evaluated_on_flat_branch") is True,
            "issues": issues,
        }

    def _check_true_loop_integral_topology(self, index: int, row: Mapping[str, Any]) -> dict[str, Any] | None:
        """Guard only the pinned C5 loop-integral rows used by actual roots."""

        if row.get("true_loop_C5_original_integrand") is not True:
            return None
        issues: list[str] = []
        for field in ("exact_function_integral", "numerical_value_not_installed", "angle_endpoints_independent_of_Z"):
            if row.get(field) is not True:
                issues.append(f"node {index}: true-loop integral missing {field}")
        variable = row.get("variable")
        integrand = row.get("integrand")
        lower = row.get("lower")
        upper = row.get("upper")
        if not isinstance(variable, str) or not variable:
            issues.append(f"node {index}: true-loop integral variable is empty")
            variable = ""
        if not all(isinstance(value, int) and not isinstance(value, bool) for value in (integrand, lower, upper)):
            issues.append(f"node {index}: true-loop integral child handles are malformed")
            return {"applied": True, "passed": not issues, "issues": issues}
        try:
            body_reachable = _reachable(self.nodes, [int(integrand)])
        except GraphProofError as exc:
            issues.append(f"node {index}: true-loop integral body traversal failed: {exc}")
            body_reachable = set()
        body_handles = {
            child
            for child in body_reachable
            if self.row(child).get("operation") == "bound_variable" and self.row(child).get("name") == variable
        }
        if len(body_handles) != 1:
            issues.append(f"node {index}: true-loop body does not expose one bound handle for {variable!r}")
        lower_row = self.row(int(lower))
        lower_value = None
        if lower_row.get("operation") == "exact_rational":
            numerator = lower_row.get("numerator")
            denominator = lower_row.get("denominator")
            if isinstance(numerator, int) and not isinstance(numerator, bool) and isinstance(denominator, int) and not isinstance(denominator, bool) and denominator != 0:
                lower_value = Fraction(numerator, denominator)
        if lower_value != Fraction(0):
            issues.append(f"node {index}: true-loop lower endpoint is not exact zero")
        upper_row = self.row(int(upper))
        upper_matches = upper_row.get("operation") == "bound_variable" and upper_row.get("name") == variable
        if not upper_matches:
            issues.append(f"node {index}: true-loop upper endpoint is not its bound variable")
        elif body_handles and int(upper) != next(iter(body_handles)):
            issues.append(f"node {index}: true-loop upper handle differs from body bound handle")
        return {
            "applied": True,
            "passed": not issues,
            "variable": variable,
            "integrand_handle": int(integrand),
            "lower_handle": int(lower),
            "upper_handle": int(upper),
            "body_bound_handles": sorted(body_handles),
            "lower_exact_value": str(lower_value) if lower_value is not None else None,
            "upper_is_same_bound_handle": bool(upper_matches and body_handles and int(upper) == next(iter(body_handles))),
            "issues": issues,
        }

    def _check_lazy_metadata(self, index: int, row: Mapping[str, Any]) -> None:
        required = ("active_body", "Delta", "eta", "flat_value", "flat_predicate")
        missing = [field for field in required if field not in row]
        bad = [field for field in required if field in row and (not isinstance(row[field], int) or isinstance(row[field], bool))]
        if missing or bad:
            self.lazy_metadata_issues.append(
                f"node {index}: lazy fields missing={missing!r} non_integer={bad!r}"
            )
        if row.get("active_body_not_evaluated_on_flat_branch") is not True:
            self.lazy_metadata_issues.append(f"node {index}: missing active_body_not_evaluated_on_flat_branch")
        sigma_source_guard = None
        if row.get("original_sigma_exact_logistic_first_derivative") is True:
            sigma_source_guard = self._sigma_source_row_guard(index, row)
            if not sigma_source_guard["passed"]:
                self.lazy_metadata_issues.extend(sigma_source_guard["issues"])
        if "flat_predicate" in row and isinstance(row.get("flat_predicate"), int):
            predicate_row = self.row(int(row["flat_predicate"]))
            if predicate_row.get("operation") not in {"logical_or", "exact_real_comparison", "exact_rational", "bound_variable"}:
                self.lazy_metadata_issues.append(
                    f"node {index}: flat_predicate has operation {predicate_row.get('operation')!r}"
                )
        if missing or bad or not all(isinstance(row.get(field), int) and not isinstance(row.get(field), bool) for field in required):
            return

        # Bind the source sigma-prime construction to its actual active DAG.
        # Metadata alone is insufficient: a malicious or stale row could carry
        # the right flags while exposing a different active expression.
        x = self.canonical(int(row["Delta"]))
        sigma = self.dag.unary("original_flat_sigma", x)
        complement = self.dag.add(self.dag.one, self.dag.neg(sigma))
        x_complement = self.dag.add(self.dag.one, self.dag.neg(x))
        first = self.dag.add(
            self.dag.quotient(self.dag.rat(2), self.dag.mul(x, x, x)),
            self.dag.quotient(self.dag.rat(2), self.dag.mul(x_complement, x_complement, x_complement)),
        )
        sigma_prime_active = self.dag.mul(sigma, complement, first)
        actual_active = self.canonical(int(row["active_body"]))
        if row.get("original_sigma_exact_logistic_first_derivative") is True and row.get("same_original_lazy_native_derivative") is not True:
            passed = actual_active == sigma_prime_active
            receipt = {
                "lazy_node": int(index),
                "lemma": "sigma_prime_active_body=sigma(x)*(1-sigma(x))*(2/x^3+2/(1-x)^3)",
                "active_body_handle": int(row["active_body"]),
                "expected_active_dag_digest": self.dag.digest(sigma_prime_active),
                "actual_active_dag_digest": self.dag.digest(actual_active),
                "expected_active_dag_node": sigma_prime_active,
                "actual_active_dag_node": actual_active,
                "expected_active_dag_nodes": len(self.dag.nodes),
                "passed": passed,
            }
            if sigma_source_guard is not None:
                receipt["source_row_guard"] = sigma_source_guard
            self.lazy_lemma_receipts.append(receipt)
            if not passed:
                self.lazy_metadata_issues.append(
                    f"node {index}: sigma-prime active DAG mismatch expected {self.dag.digest(sigma_prime_active)} "
                    f"actual {self.dag.digest(actual_active)}"
                )
        elif (
            row.get("same_original_lazy_native_derivative") is True
            and row.get("original_sigma_exact_logistic_first_derivative") is True
        ):
            # The issued sigma'' row keeps the full quotient-rule DAG.  Check
            # it against the independent construction while retaining all
            # quotient numerator/denominator factors in the receipt.
            expected_lazy = self._sigma_prime_second(int(row["Delta"]), self.dag.one)
            expected_active = self.dag.node(expected_lazy)[1][0]
            passed = actual_active == expected_active
            receipt = {
                "lazy_node": int(index),
                "lemma": "sigma_second_active_body=independent_quotient_rule",
                "active_body_handle": int(row["active_body"]),
                "expected_active_dag_digest": self.dag.digest(expected_active),
                "actual_active_dag_digest": self.dag.digest(actual_active),
                "expected_active_dag_node": expected_active,
                "actual_active_dag_node": actual_active,
                "expected_active_dag_nodes": len(self.dag.nodes),
                "passed": passed,
            }
            if sigma_source_guard is not None:
                receipt["source_row_guard"] = sigma_source_guard
            self.lazy_lemma_receipts.append(receipt)
            if not passed:
                self.lazy_metadata_issues.append(
                    f"node {index}: sigma-second active DAG mismatch expected {self.dag.digest(expected_active)} "
                    f"actual {self.dag.digest(actual_active)}"
                )

    def _check_binder_scope(self, index: int, row: Mapping[str, Any], operation: str) -> None:
        if operation == "function_substitution":
            variable = row.get("variable")
            if not isinstance(variable, int) or isinstance(variable, bool) or self.row(variable).get("operation") != "bound_variable":
                self.binder_scope_issues.append(f"node {index}: substitution variable is not a bound variable")
            if row.get("Z_independent_substitution") is not True:
                self.binder_scope_issues.append(f"node {index}: substitution lacks Z_independent_substitution guard")
            # Upstream function substitutions predate the NativeDerivative
            # builder and therefore do not carry its optional exact-chain
            # marker.  If the marker is present, it must be true; absence is
            # recorded as an upstream binder rather than misclassified.
            if "exact_native_chain_rule_substitution" in row and row.get("exact_native_chain_rule_substitution") is not True:
                self.binder_scope_issues.append(f"node {index}: substitution lacks exact_native_chain_rule_substitution guard")
            self.binder_scope_receipts.append({
                "node": int(index),
                "binder_kind": "function_substitution",
                "variable_handle": int(variable) if isinstance(variable, int) else None,
                "variable_name": self.name(int(variable)) if isinstance(variable, int) else None,
                "exact_chain_guard": row.get("exact_native_chain_rule_substitution") is True,
                "upstream_guard_only": "exact_native_chain_rule_substitution" not in row,
                "capture_check": "typed variable handle; no same-name identification",
                "passed": isinstance(variable, int) and not isinstance(variable, bool) and self.row(variable).get("operation") == "bound_variable",
            })
        elif operation == "definite_integral":
            variable = str(row.get("variable"))
            if not variable:
                self.binder_scope_issues.append(f"node {index}: integral binder has empty variable name")
            topology_guard = self._check_true_loop_integral_topology(index, row)
            receipt = {
                "node": int(index),
                "binder_kind": "definite_integral",
                "variable_name": variable,
                "lower_handle": row.get("lower"),
                "upper_handle": row.get("upper"),
                "capture_check": "endpoint handles retained; endpoint names are not used as binder identity",
                "passed": bool(variable) and all(isinstance(row.get(field), int) and not isinstance(row.get(field), bool) for field in ("lower", "upper")),
            }
            if topology_guard is not None:
                receipt["source_topology_guard"] = topology_guard
                if not topology_guard["passed"]:
                    self.binder_scope_issues.extend(topology_guard["issues"])
            self.binder_scope_receipts.append(receipt)
        elif operation in ("substitute_original_inverse_angle", "original_monotone_phase_inverse"):
            angle = str(row.get("angle_variable"))
            child = row.get("inverse_angle") if operation == "substitute_original_inverse_angle" else row.get("target_phase")
            if not angle:
                self.binder_scope_issues.append(f"node {index}: inverse-angle binder has empty angle name")
            self.binder_scope_receipts.append({
                "node": int(index),
                "binder_kind": operation,
                "angle_name": angle,
                "child_handle": child,
                "capture_check": "typed child handle retained; angle names are not used as variable identity",
                "passed": bool(angle) and isinstance(child, int) and not isinstance(child, bool),
            })

    def canonical(self, index: int) -> int:
        self._tick()
        old = self._canonical_cache.get(index)
        if old is not None:
            return old
        row = self.row(index)
        operation = str(row.get("operation"))
        if operation == "exact_rational":
            value = Fraction(int(row["numerator"]), int(row["denominator"]))
            result = self.dag.rat(value)
        elif operation == "bound_variable":
            result = self.dag.var(str(row["name"]))
        elif operation == "sum":
            result = self.dag.add(*(self.canonical(child) for child in row["arguments"]))
        elif operation == "product":
            result = self.dag.mul(*(self.canonical(child) for child in row["arguments"]))
        elif operation == "negative":
            result = self.dag.neg(self.canonical(int(row["argument"])))
        elif operation == "positive_quotient":
            result = self.dag.quotient(self.canonical(int(row["numerator"])), self.canonical(int(row["denominator"])))
        elif operation in ("analytic_unary", "exact_original_source_flat_sigma"):
            name = "original_flat_sigma" if operation == "exact_original_source_flat_sigma" else str(row["name"])
            result = self.dag.unary(name, self.canonical(int(row["argument"])))
        elif operation == "function_substitution":
            self._check_binder_scope(index, row, operation)
            result = self.dag.substitute(
                self.canonical(int(row["expression"])),
                self.name(int(row["variable"])),
                self.canonical(int(row["value"])),
            )
        elif operation == "definite_integral":
            self._check_binder_scope(index, row, operation)
            result = self.dag.integral(
                self.canonical(int(row["integrand"])),
                str(row["variable"]),
                self.canonical(int(row["lower"])),
                self.canonical(int(row["upper"])),
            )
        elif operation == "substitute_original_inverse_angle":
            self._check_binder_scope(index, row, operation)
            result = self.dag.at_inverse(
                self.canonical(int(row["body"])),
                self.canonical(int(row["inverse_angle"])),
                str(row["angle_variable"]),
            )
        elif operation == "original_monotone_phase_inverse":
            self._check_binder_scope(index, row, operation)
            result = self.dag.phase_inverse(
                self.canonical(int(row["phase_function"])),
                str(row["angle_variable"]),
                self.canonical(int(row["target_phase"])),
            )
        elif operation == "original_lazy_flat_branch":
            self._check_lazy_metadata(index, row)
            if (
                row.get("original_sigma_exact_logistic_first_derivative") is True
                and row.get("same_original_lazy_native_derivative") is not True
            ):
                # NativeDerivative.sigma_prime emits this issued lazy row.
                # Keep one semantic atom for it so a generated sigma-prime
                # row compares with the unary used by the independent rule.
                result = self.dag.unary("sigma_prime", self.canonical(int(row["Delta"])))
            else:
                result = self.dag.lazy(
                    self.canonical(int(row["active_body"])),
                    self.canonical(int(row["Delta"])),
                    self.canonical(int(row["eta"])),
                    self.canonical(int(row["flat_value"])),
                    self.canonical(int(row["flat_predicate"])),
                )
        elif operation == "exact_real_comparison":
            result = self.dag.comparison(
                str(row.get("operator")),
                self.canonical(int(row["left"])),
                self.canonical(int(row["right"])),
            )
        elif operation == "logical_or":
            result = self.dag.logical_or(*(self.canonical(child) for child in row["arguments"]))
        else:
            result = self.dag.source(operation, index, self.raw_digest(index))
        self._canonical_cache[index] = result
        return result

    def free_variables(self, index: int) -> frozenset[str]:
        self._tick()
        old = self._free_cache.get(index)
        if old is not None:
            return old
        row = self.row(index)
        operation = str(row.get("operation"))
        visit = self.free_variables
        if operation == "bound_variable":
            result: set[str] = {str(row["name"])}
        elif operation in self.EMPTY_OPERATIONS:
            result = set()
        elif operation == "exact_real_comparison":
            result = set(visit(int(row["left"]))) | set(visit(int(row["right"])))
        elif operation == "logical_or":
            result = set().union(*(visit(int(child)) for child in row["arguments"]))
        elif operation == "original_lazy_flat_branch":
            result = set()
            for field in GRAPH_FIELDS[operation]:
                values = row.get(field, ())
                if not isinstance(values, list):
                    values = (values,)
                for child in values:
                    result.update(visit(int(child)))
        elif operation == "original_monotone_phase_inverse":
            result = (set(visit(int(row["phase_function"]))) - {str(row["angle_variable"])}) | set(visit(int(row["target_phase"])))
        elif operation == "substitute_original_inverse_angle":
            result = (set(visit(int(row["body"]))) - {str(row["angle_variable"])}) | set(visit(int(row["inverse_angle"])))
        elif operation == "function_substitution":
            body_vars = set(visit(int(row["expression"])))
            name = self.name(int(row["variable"]))
            result = (body_vars - {name}) | set(visit(int(row["value"]))) if name in body_vars else body_vars
        elif operation == "definite_integral":
            result = (
                set(visit(int(row["integrand"]))) - {str(row["variable"])}
            ) | set(visit(int(row["lower"]))) | set(visit(int(row["upper"])))
        elif operation in ("sum", "product"):
            result = set().union(*(visit(int(child)) for child in row["arguments"]))
        elif operation == "positive_quotient":
            result = set(visit(int(row["numerator"]))) | set(visit(int(row["denominator"])))
        elif operation in ("negative", "analytic_unary", "exact_original_source_flat_sigma"):
            result = set(visit(int(row["argument"])))
        elif operation == "current_original_unique_analytic_core_component":
            result = set(visit(int(row["rho"]))) | set(visit(int(row["Z_variable"])))
        elif operation in ("current_original_leading_function_recipe", "current_original_axis_pressure_function_derivative"):
            if row.get("coordinate_independent") is True and row.get("phase_independent_leading_source") is True:
                result = set()
            else:
                result = set(visit(int(row["Z_variable"])))
        elif operation.startswith("current_original_"):
            if row.get("Z_independent") is True or row.get("phase_independent_leading_source") is True:
                result = set()
            else:
                self.open_rules.add("free_variables:" + operation)
                result = {"<unsupported:" + operation + ">"}
        else:
            self.open_rules.add("free_variables:" + operation)
            result = {"<unsupported:" + operation + ">"}
        frozen = frozenset(result)
        self._free_cache[index] = frozen
        return frozen

    @staticmethod
    def _token_name(token: tuple[Any, ...]) -> str:
        return str(token[1])

    def _core_extension(self, index: int, field: str) -> int:
        row = self.row(index)
        key = (
            int(row.get("rho", -1)),
            int(row.get("Z_variable", -1)),
            int(row.get("radial_order", -1)) + (1 if field == "radial_order" else 0),
            int(row.get("ordinary_slow_Z_derivative_order", -1)) + (1 if field == "ordinary_slow_Z_derivative_order" else 0),
        )
        candidates = self._core_rows.get(key, ())
        if not candidates:
            self.open_rules.add("missing_core_extension:" + field)
            raise OpenRuleError(f"missing core extension for node {index} field {field}")
        return self.canonical(candidates[0])

    def _sigma_prime_second(self, argument: int, d_argument: int) -> int:
        dag = self.dag
        x = self.canonical(argument)
        sigma = dag.unary("original_flat_sigma", x)
        sigma_comp = dag.add(dag.one, dag.neg(sigma))
        x_comp = dag.add(dag.one, dag.neg(x))
        first = dag.add(
            dag.quotient(dag.rat(2), dag.mul(x, x, x)),
            dag.quotient(dag.rat(2), dag.mul(x_comp, x_comp, x_comp)),
        )
        d_sigma = dag.mul(dag.unary("sigma_prime", x), d_argument)
        d_comp = dag.neg(d_sigma)
        # Preserve the issued quotient-rule graph.  In particular, do not
        # collapse ``-(2*3*x^2)/(x^3*x^3)`` to ``-6/x^4``: the source receipt
        # binds the numerator and denominator factors separately.
        x2_sum = dag.add(dag.mul(x, x), dag.mul(x, x), dag.mul(x, x))
        x3 = dag.mul(x, x, x)
        first_numerator = dag.neg(dag.mul(dag.rat(2), x2_sum))
        first_denominator = dag.mul(x3, x3)
        comp2_sum = dag.add(
            dag.mul(dag.neg(dag.one), x_comp, x_comp),
            dag.mul(dag.neg(dag.one), x_comp, x_comp),
            dag.mul(dag.neg(dag.one), x_comp, x_comp),
        )
        second_numerator = dag.neg(dag.mul(dag.rat(2), comp2_sum))
        comp3 = dag.mul(x_comp, x_comp, x_comp)
        d_first = dag.mul(
            dag.add(
                dag.quotient(first_numerator, first_denominator),
                dag.quotient(second_numerator, dag.mul(comp3, comp3)),
            ),
            d_argument,
        )
        active = dag.add(
            dag.mul(d_sigma, sigma_comp, first),
            dag.mul(sigma, d_comp, first),
            dag.mul(sigma, sigma_comp, d_first),
        )
        predicate = dag.logical_or(
            dag.comparison("<=", x, dag.zero),
            dag.comparison(">=", x, dag.one),
        )
        return dag.lazy(active, x, dag.one, dag.zero, predicate)

    def derivative(
        self,
        index: int,
        variable_name: str,
        variable_token: tuple[Any, ...],
        phase: int | None,
        seed_map: Mapping[int, int],
        native_token: tuple[Any, ...],
        seed_profile: str | None = None,
    ) -> int:
        self._tick()
        if seed_profile is None:
            seed_profile = self.seed_profile(seed_map)
        cache_key = (int(index), tuple(variable_token), phase, seed_profile)
        old = self._derivative_cache.get(cache_key)
        if old is not None:
            return old
        if phase is not None and index == phase:
            self._derivative_cache[cache_key] = self.dag.zero
            return self.dag.zero
        if variable_token == native_token and index in seed_map:
            result = seed_map[index]
            self._derivative_cache[cache_key] = result
            return result
        if variable_name not in self.free_variables(index):
            self._derivative_cache[cache_key] = self.dag.zero
            return self.dag.zero
        row = self.row(index)
        operation = str(row.get("operation"))
        D = lambda child: self.derivative(
            int(child), variable_name, variable_token, phase, seed_map, native_token, seed_profile
        )
        C = lambda child: self.canonical(int(child))
        dag = self.dag
        if operation == "bound_variable":
            result = dag.one if str(row["name"]) == variable_name else dag.zero
        elif operation == "sum":
            result = dag.add(*(D(child) for child in row["arguments"]))
        elif operation == "product":
            args = [int(child) for child in row["arguments"]]
            result = dag.add(*(
                dag.mul(*(D(child) if position == selected else C(child) for position, child in enumerate(args)))
                for selected in range(len(args))
            ))
        elif operation == "negative":
            result = dag.neg(D(row["argument"]))
        elif operation == "positive_quotient":
            numerator = int(row["numerator"])
            denominator = int(row["denominator"])
            result = dag.quotient(
                dag.add(dag.mul(D(numerator), C(denominator)), dag.neg(dag.mul(C(numerator), D(denominator)))),
                dag.mul(C(denominator), C(denominator)),
            )
        elif operation in ("analytic_unary", "exact_original_source_flat_sigma"):
            argument = int(row["argument"])
            x = C(argument)
            dx = D(argument)
            name = "original_flat_sigma" if operation == "exact_original_source_flat_sigma" else str(row["name"])
            if name in ("exp", "expm1"):
                outer = dag.unary("exp", x)
            elif name == "log":
                outer = dag.quotient(dag.one, x)
            elif name == "log1p":
                outer = dag.quotient(dag.one, dag.add(dag.one, x))
            elif name in ("sqrt", "positive_sqrt"):
                outer = dag.quotient(dag.one, dag.mul(dag.rat(2), dag.unary("positive_sqrt", x)))
            elif name == "cos":
                outer = dag.neg(dag.unary("sin", x))
            elif name == "sin":
                outer = dag.unary("cos", x)
            elif name == "original_flat_sigma":
                outer = dag.unary("sigma_prime", x)
            elif name == "original_flat_sigma_prime":
                result = self._sigma_prime_second(argument, dx)
                self._derivative_cache[cache_key] = result
                return result
            else:
                self.open_rules.add("analytic_unary:" + name)
                raise OpenRuleError("unissued analytic unary in active derivative: " + name)
            result = dag.mul(outer, dx)
        elif operation == "function_substitution":
            body = int(row["expression"])
            variable = int(row["variable"])
            argument = int(row["value"])
            substituted_name = self.name(variable)
            direct = dag.zero if variable_token == ("node", variable) else dag.substitute(D(body), substituted_name, C(argument))
            other_token = ("node", variable)
            other = self.derivative(
                body,
                substituted_name,
                other_token,
                phase,
                {},
                native_token,
                self.seed_profile({}),
            )
            chain = dag.mul(dag.substitute(other, substituted_name, C(argument)), D(argument))
            result = dag.add(direct, chain)
        elif operation == "definite_integral":
            body = int(row["integrand"])
            lower = int(row["lower"])
            upper = int(row["upper"])
            integration_variable = str(row["variable"])
            result = dag.zero
            if integration_variable != variable_name:
                partial = D(body)
                if partial != dag.zero:
                    result = dag.add(result, dag.integral(partial, integration_variable, C(lower), C(upper)))
            result = dag.add(
                result,
                dag.mul(dag.substitute(C(body), integration_variable, C(upper)), D(upper)),
                dag.neg(dag.mul(dag.substitute(C(body), integration_variable, C(lower)), D(lower))),
            )
        elif operation == "substitute_original_inverse_angle":
            body = int(row["body"])
            inverse = int(row["inverse_angle"])
            angle = str(row["angle_variable"])
            direct = dag.zero if variable_name == angle else dag.at_inverse(D(body), C(inverse), angle)
            other_token = ("name", angle)
            other = self.derivative(
                body, angle, other_token, phase, {}, native_token, self.seed_profile({})
            )
            result = dag.add(direct, dag.mul(dag.at_inverse(other, C(inverse), angle), D(inverse)))
        elif operation == "original_monotone_phase_inverse":
            phase_function = int(row["phase_function"])
            target = int(row["target_phase"])
            angle = str(row["angle_variable"])
            inverse = C(index)
            held = dag.zero if variable_name == angle else D(phase_function)
            other_token = ("name", angle)
            other = self.derivative(
                phase_function,
                angle,
                other_token,
                phase,
                {},
                native_token,
                self.seed_profile({}),
            )
            lam = dag.at_inverse(other, inverse, angle)
            result = dag.quotient(dag.add(D(target), dag.neg(dag.at_inverse(held, inverse, angle))), lam)
        elif operation == "original_lazy_flat_branch":
            result = dag.lazy(
                D(row["active_body"]),
                C(row["Delta"]),
                C(row["eta"]),
                D(row["flat_value"]),
                C(row["flat_predicate"]),
            )
        elif operation == "current_original_unique_analytic_core_component":
            if not (
                row.get("source_family") is not None
                and row.get("unique_analytic_Banach_solution") is True
                and row.get("actual_twenty_term_operator") is True
                and row.get("admitted_solution_not_a_model_or_enclosure_center") is True
            ):
                self.open_rules.add("analytic_core_premise")
                raise OpenRuleError("analytic core source premises are not admitted")
            terms: list[int] = []
            dr = D(row["rho"])
            dz = D(row["Z_variable"])
            if dr != dag.zero:
                if int(row.get("radial_order", 0)) >= 2:
                    self.open_rules.add("analytic_core_radial_order>2")
                    raise OpenRuleError("additional radial core source required")
                terms.append(dag.mul(self._core_extension(index, "radial_order"), dr))
            if dz != dag.zero:
                terms.append(dag.mul(self._core_extension(index, "ordinary_slow_Z_derivative_order"), dz))
            result = dag.add(*terms)
        else:
            self.open_rules.add("derivative:" + operation)
            raise OpenRuleError("unissued actual native derivative operation: " + operation)
        self._derivative_cache[cache_key] = result
        return result


def _parent_validation(root: str | Path, parent: Any) -> dict[str, Any]:
    try:
        path = _locate(root, PARENT_NAME)
        actual_hash = _sha256(path)
        disk = _read_json(path)
        disk_digest = _mapping_digest(disk)
    except Exception as exc:  # pragma: no cover - receipt path
        return {
            "path": str(PARENT_NAME),
            "expected_sha256": EXPECTED_PARENT_SHA256,
            "actual_sha256": None,
            "mapping_sha256": None,
            "passed": False,
            "mapping_exact": False,
            "reason": f"parent read failed: {exc}",
        }
    supplied = parent
    supplied_path: str | None = None
    if isinstance(parent, (str, Path)):
        supplied_path = str(parent)
        try:
            supplied = _read_json(Path(parent))
        except Exception as exc:
            supplied = None
            supplied_path = f"{parent} ({exc})"
    mapping_exact = isinstance(supplied, Mapping) and supplied == disk
    return {
        "path": str(path),
        "supplied_path": supplied_path,
        "expected_sha256": EXPECTED_PARENT_SHA256,
        "actual_sha256": actual_hash,
        "mapping_sha256": disk_digest,
        "supplied_mapping_sha256": _mapping_digest(supplied) if isinstance(supplied, Mapping) else None,
        "mapping_exact": mapping_exact,
        "passed": actual_hash == EXPECTED_PARENT_SHA256 and mapping_exact,
        "source_executed": False,
        "reason": None if mapping_exact else "caller parent mapping is not the accepted disk mapping",
    }


def _source_pins(root: str | Path, self_path: Path) -> tuple[dict[str, Any], bool]:
    pins: dict[str, Any] = {}
    passed = True
    for relative, expected in EXPECTED_INPUT_SHA256.items():
        try:
            path = _locate(root, Path(relative).name)
            actual = _sha256(path)
            ok = actual == expected
            path_string = relative
        except Exception as exc:
            actual = None
            ok = False
            path_string = relative
            pins[relative] = {
                "path": path_string,
                "expected_sha256": expected,
                "actual_sha256": actual,
                "passed": False,
                "source_executed": False,
                "reason": str(exc),
            }
            passed = False
            continue
        pins[relative] = {
            "path": path_string,
            "expected_sha256": expected,
            "actual_sha256": actual,
            "passed": ok,
            "source_executed": False,
        }
        passed = passed and ok
    self_relative = "src/openai_ns_reconstruction/original_physical_y_dag_derivative.py"
    self_hash = _sha256(self_path)
    pins[self_relative] = {
        "path": self_relative,
        "expected_sha256": self_hash,
        "actual_sha256": self_hash,
        "passed": True,
        "source_executed": False,
        "self_binding": True,
    }
    return pins, passed


def _seed_bindings(
    oracle: _Oracle,
    dag: _DAG,
    chart: str,
    c6_bindings: Mapping[str, Any],
    q4_bindings: Mapping[str, Any],
    native_alias: str | None = None,
) -> tuple[dict[int, int], dict[str, Any], list[dict[str, Any]]]:
    c6_name = C6_CHART_FOR_Q4[chart]
    c6_chart = c6_bindings["charts"][c6_name]
    q4_chart = q4_bindings["charts"][chart]
    J = int(c6_chart["physical_measure"])
    if native_alias is None:
        native_alias = str(oracle.row(int(q4_chart.get("coordinate", c6_chart["phase"]))).get("name", q4_chart.get("coordinate", c6_chart["phase"])))
    radial = q4_chart["radial_background"]
    seed_map: dict[int, int] = {}
    field_rows: dict[str, list[dict[str, Any]]] = {}
    table: list[dict[str, Any]] = []
    for field in FIELD_NAMES:
        source_rows = c6_chart["fields"].get(field)
        radial_rows = radial.get(field)
        if not isinstance(source_rows, list) or not isinstance(radial_rows, list) or len(source_rows) != 7 or len(radial_rows) != 7:
            raise GraphProofError(f"seed row length mismatch {chart}.{field}")
        receipts: list[dict[str, Any]] = []
        for order, (source, radial_row) in enumerate(zip(source_rows, radial_rows)):
            source = int(source)
            radial_row = int(radial_row)
            expression = dag.mul(oracle.canonical(J), oracle.canonical(radial_row))
            old = seed_map.get(source)
            if old is not None and old != expression:
                raise GraphProofError(f"conflicting seed for source node {source}")
            seed_map[source] = expression
            row = {
                "chart": chart,
                "native_alias": native_alias,
                "J_handle": J,
                "field_name": field,
                "ordinary_order": order,
                "background_field_handle": source,
                "radial_background_handle": radial_row,
                "seed_Dq_formula": "J*background_radial_row",
                "seed_dag_digest": dag.digest(expression),
                "seed_dag_node_count": len(dag.nodes),
            }
            table.append(row)
            receipts.append({
                "ordinary_order": order,
                "source_field_node": source,
                "physical_radial_row_node": radial_row,
                "seed_expression_sha256": dag.digest(expression),
            })
        field_rows[field] = receipts
    receipt = {
        "variable_phase_handle": int(c6_chart["phase"]),
        "native_alias": native_alias,
        "held_phase_handle": int(q4_chart.get("actual_phase", c6_chart["phase"])),
        "c6_chart_name": c6_name,
        "physical_J_node": J,
        "field_rows": field_rows,
        "seed_identity_obligation": "D_native(field)=J*background_radial_row for every source field row",
        "seed_identity_discharged": False,
    }
    return seed_map, receipt, table


def _positive_j_check(oracle: _Oracle, chart: str, c5_chart: Mapping[str, Any], c6_chart: Mapping[str, Any], q4_chart: Mapping[str, Any]) -> dict[str, Any]:
    J = int(c6_chart["physical_measure"])
    roots = list(q4_chart["radial"]["slow_A_y"]) + list(q4_chart["radial"]["slow_B_y"])
    inverses: list[dict[str, Any]] = []
    for root in roots:
        row = oracle.row(int(root))
        args = row.get("arguments", [])
        inv = int(args[0]) if row.get("operation") == "product" and len(args) == 2 else None
        inv_row = oracle.row(inv) if inv is not None else {}
        ok = bool(
            inv is not None
            and inv_row.get("operation") == "positive_quotient"
            and inv_row.get("numerator") == 1
            and inv_row.get("denominator") == J
            and inv_row.get("source_positive_certificate") == "positive original chart log-radius Jacobian"
        )
        inverses.append({"root": int(root), "inverse_J_node": inv, "passed": ok})
    coordinate = int(c5_chart["coordinate"])
    coordinate_name = oracle.name(coordinate)
    j_free = oracle.free_variables(J)
    inverse_binding_passed = bool(inverses) and all(item["passed"] for item in inverses)
    j_independent_of_coordinate = coordinate_name not in j_free
    j_independent_of_z = "Z" not in j_free
    return {
        "passed": inverse_binding_passed and j_independent_of_coordinate and j_independent_of_z,
        "inverse_J_binding_passed": inverse_binding_passed,
        "physical_J_node": J,
        "physical_J_operation": oracle.row(J).get("operation"),
        "physical_J_free_variables": sorted(j_free),
        "J_independent_of_native_coordinate": j_independent_of_coordinate,
        "J_independent_of_Z": j_independent_of_z,
        "J_independence_gate_passed": j_independent_of_coordinate and j_independent_of_z,
        "inverse_J_rows": inverses,
        "positivity_source": "control/composition source certificate; positivity is not discharged by this DAG module",
        "positive_J_not_numerically_evaluated": True,
    }


def _radius_relation(oracle: _Oracle, c5_chart: Mapping[str, Any], dag: _DAG) -> dict[str, Any]:
    coordinate = int(c5_chart["coordinate"])
    radius = int(c5_chart["radius"])
    J = int(c5_chart["physical_measure"])
    token = ("node", coordinate)
    native_name = oracle.name(coordinate)
    previous_deadline = oracle.deadline
    try:
        oracle.deadline = time.monotonic() + 10.0
        derivative = oracle.derivative(radius, native_name, token, None, {}, token)
        expected = dag.mul(oracle.canonical(J), oracle.canonical(radius))
        passed = derivative == expected
        reason = None if passed else "D_native(R_q) did not normalize to J*R"
    except Exception as exc:
        derivative = dag.source("radius_check_failure", radius)
        expected = dag.source("radius_check_expected", radius)
        passed = False
        reason = f"{type(exc).__name__}: {exc}"
    finally:
        oracle.deadline = previous_deadline
    return {
        "coordinate_node": coordinate,
        "radius_node": radius,
        "J_node": J,
        "R_q_equals_JR": passed,
        "R_q_sha256": dag.digest(derivative),
        "JR_sha256": dag.digest(expected),
        "reason": reason,
    }


def _root_check(
    oracle: _Oracle,
    dag: _DAG,
    source: int,
    root: int,
    native_name: str,
    coordinate: int,
    phase: int,
    seed_map: Mapping[int, int],
    budget_seconds: float = 10.0,
) -> dict[str, Any]:
    started = time.monotonic()
    before_nodes = len(dag.nodes)
    actual_row = oracle.row(root)
    expected_id: int | None = None
    actual_id: int | None = None
    source_factor: int | None = None
    seed_profile = oracle.seed_profile(seed_map)
    reason: str | None = None
    passed = False
    oracle.deadline = started + budget_seconds
    try:
        if actual_row.get("operation") != "product" or len(actual_row.get("arguments", ())) != 2:
            raise GraphProofError("root is not product(invJ, derivative)")
        inv = int(actual_row["arguments"][0])
        source_factor = int(actual_row["arguments"][1])
        token = ("node", coordinate)
        derivative = oracle.derivative(source, native_name, token, phase, seed_map, token, seed_profile)
        expected_id = dag.mul(oracle.canonical(inv), derivative)
        actual_id = oracle.canonical(root)
        passed = actual_id == expected_id
        if not passed:
            reason = "canonical expected DAG differs from issued root DAG"
    except Exception as exc:
        reason = f"{type(exc).__name__}: {exc}"
    finally:
        oracle.deadline = None
    elapsed = time.monotonic() - started
    if elapsed > budget_seconds and passed:
        passed = False
        reason = "root check exceeded time budget"
    if expected_id is None:
        expected_id = dag.source("expected_root_failure", root)
    if actual_id is None:
        actual_id = dag.source("actual_root_failure", root)
    local_mismatch = None if passed else _first_dag_mismatch(dag, expected_id, actual_id)
    return {
        "source_handle": int(source),
        "root_handle": int(root),
        "issued_derivative_handle": source_factor,
        "native_coordinate_handle": int(coordinate),
        "held_phase_handle": int(phase),
        "inverse_J_handle": int(actual_row.get("arguments", [None])[0]) if actual_row.get("operation") == "product" else None,
        "expected_dag_digest": dag.digest(expected_id),
        "actual_dag_digest": dag.digest(actual_id),
        "seed_profile": seed_profile,
        "passed": passed,
        "failure_reason": reason,
        "first_local_mismatch": local_mismatch,
        "elapsed_seconds": elapsed,
        "dag_nodes_before": before_nodes,
        "dag_nodes_after": len(dag.nodes),
        "dag_nodes_added": len(dag.nodes) - before_nodes,
        "oracle_steps": oracle._ticks,
    }


def build(root: str | Path, parent: Any = None) -> dict[str, Any]:
    """Build a JSON-friendly conditional receipt from saved graph packets."""

    started = time.monotonic()
    module_path = Path(__file__).resolve()
    parent_receipt = _parent_validation(root, parent)
    source_paths: dict[str, Path] = {}
    archive_paths: dict[str, Path] = {}
    check_paths: dict[str, Path] = {}
    for tag, stem in (("q4", Q4_STEM), ("c5", C5_STEM), ("c6", C6_STEM)):
        source_paths[tag] = _locate(root, stem + ".py")
        archive_paths[tag] = _locate(root, stem + ".json.gz")
        check_paths[tag] = _locate(root, stem + "_check.json")
    pins, pins_ok = _source_pins(root, module_path)
    q4 = _read_gzip_json(archive_paths["q4"])
    c5 = _read_gzip_json(archive_paths["c5"])
    c6 = _read_gzip_json(archive_paths["c6"])
    nodes = q4.get("exact_graph_nodes")
    c5_nodes = c5.get("exact_graph_nodes")
    c6_nodes = c6.get("exact_graph_nodes")
    if not isinstance(nodes, list) or not isinstance(c5_nodes, list) or not isinstance(c6_nodes, list):
        raise GraphProofError("missing exact graph node arrays")
    q4_prefix = nodes[: len(c5_nodes)] == c5_nodes
    c5_prefix = c5_nodes[: len(c6_nodes)] == c6_nodes
    q4_bindings = q4.get("modified_Q4_pressure_C3_bindings")
    c5_bindings = c5.get("prescribed_loop_C5_bindings")
    c6_bindings = c6.get("prescribed_R100_Rm_C6_bindings")
    if not all(isinstance(value, Mapping) for value in (q4_bindings, c5_bindings, c6_bindings)):
        raise GraphProofError("missing saved source bindings")
    q4_charts = q4_bindings.get("charts")
    c5_charts = c5_bindings.get("charts")
    c6_charts = c6_bindings.get("charts")
    if not all(isinstance(value, Mapping) for value in (q4_charts, c5_charts, c6_charts)):
        raise GraphProofError("missing chart binding maps")
    dag = _DAG()
    oracle = _Oracle(nodes, dag)
    rule_binding = _native_rule_binding(source_paths["q4"])
    chart_receipts: dict[str, Any] = {}
    seed_table: list[dict[str, Any]] = []
    all_roots: list[dict[str, Any]] = []
    all_radius = True
    all_j = True
    all_j_independence = True
    all_phase_zero = True
    all_unsupported_clear = True
    all_open_rules_clear = True
    all_rule_checks = True
    all_roots_passed = True
    for chart in SCOPE_CHARTS:
        c5_chart = c5_charts[chart]
        c6_chart = c6_charts[C6_CHART_FOR_Q4[chart]]
        q4_chart = q4_charts[chart]
        if c5_chart["physical_measure"] != c6_chart["physical_measure"]:
            raise GraphProofError(f"C5/C6 physical J mismatch {chart}")
        coordinate = int(c5_chart["coordinate"])
        native_name = oracle.name(coordinate)
        chart_seed_map, seed_receipt, rows = _seed_bindings(
            oracle, dag, chart, c6_bindings, q4_bindings, native_name
        )
        seed_table.extend(rows)
        phase = int(c5_chart["actual_phase"])
        source_rows = list(c5_chart["loop"]["A"]) + list(c5_chart["loop"]["B_over_Pstar"])
        actual_rows = list(q4_chart["radial"]["slow_A_y"]) + list(q4_chart["radial"]["slow_B_y"])
        source_reachable = _reachable(nodes, source_rows)
        actual_reachable = _reachable(nodes, actual_rows)
        row_checks: list[dict[str, Any]] = []
        for field, source_key, actual_key in (
            ("A", "A", "slow_A_y"),
            ("B_over_Pstar", "B_over_Pstar", "slow_B_y"),
        ):
            source_jet = list(c5_chart["loop"][source_key])
            actual_jet = list(q4_chart["radial"][actual_key])
            if len(source_jet) != 6 or len(actual_jet) != 6:
                raise GraphProofError(f"{chart}.{field} must expose six roots")
            for order, (source, actual) in enumerate(zip(source_jet, actual_jet)):
                receipt = _root_check(
                    oracle,
                    dag,
                    int(source),
                    int(actual),
                    native_name,
                    coordinate,
                    phase,
                    chart_seed_map,
                )
                receipt.update({"chart": chart, "field": field, "ordinary_order": order})
                row_checks.append(receipt)
                all_roots.append(receipt)
        radius = _radius_relation(oracle, c5_chart, dag)
        positive_j = _positive_j_check(oracle, chart, c5_chart, c6_chart, q4_chart)
        all_radius = all_radius and bool(radius["R_q_equals_JR"])
        all_j = all_j and bool(positive_j["passed"])
        all_roots_passed = all_roots_passed and all(bool(row["passed"]) for row in row_checks)
        source_ops = _operation_counts(nodes, source_reachable)
        actual_ops = _operation_counts(nodes, actual_reachable)
        unsupported = sorted(
            operation
            for operation in set(source_ops) | set(actual_ops)
            if operation not in {
                "bound_variable", "sum", "product", "negative", "positive_quotient",
                "analytic_unary", "exact_original_source_flat_sigma", "function_substitution",
                "definite_integral", "substitute_original_inverse_angle", "original_monotone_phase_inverse",
                "original_lazy_flat_branch", "current_original_unique_analytic_core_component",
                "current_original_leading_function_recipe", "current_original_axis_pressure_function_derivative",
                "current_original_source_parameter", "exact_real_comparison", "logical_or",
                "exact_rational", "exact_positive_integer_expression", "mathematical_pi",
                "mathematical_positive_infinity",
            }
        )
        phase_zero = oracle.derivative(phase, native_name, ("node", coordinate), phase, chart_seed_map, ("node", coordinate)) == dag.zero
        open_rules_clear = not oracle.open_rules
        unsupported_clear = not unsupported
        rule_checks = {
            "sum": source_ops.get("sum", 0) > 0 and actual_ops.get("sum", 0) > 0,
            "product": source_ops.get("product", 0) > 0 and actual_ops.get("product", 0) > 0,
            "quotient": source_ops.get("positive_quotient", 0) > 0 and actual_ops.get("positive_quotient", 0) > 0,
            "unary": source_ops.get("analytic_unary", 0) > 0 and actual_ops.get("analytic_unary", 0) > 0,
            "squared_flat_sigma": source_ops.get("original_lazy_flat_branch", 0) > 0 and bool(oracle.lazy_lemma_receipts),
            "function_substitution_chain": source_ops.get("function_substitution", 0) > 0 and actual_ops.get("function_substitution", 0) > 0,
            "definite_integral_partial_and_endpoints": source_ops.get("definite_integral", 0) > 0 and actual_ops.get("definite_integral", 0) > 0,
            "inverse_angle_chain": source_ops.get("substitute_original_inverse_angle", 0) > 0 and actual_ops.get("substitute_original_inverse_angle", 0) > 0,
            "implicit_inverse": source_ops.get("original_monotone_phase_inverse", 0) > 0 and actual_ops.get("original_monotone_phase_inverse", 0) > 0,
            "lazy_flat_branch": source_ops.get("original_lazy_flat_branch", 0) > 0 and actual_ops.get("original_lazy_flat_branch", 0) > 0,
        }
        rule_checks_passed = all(rule_checks.values())
        chart_receipts[chart] = {
            "native_coordinate_handle": coordinate,
            "native_alias": native_name,
            "held_actual_phase_handle": phase,
            "physical_J_handle": int(c6_chart["physical_measure"]),
            "seed_profile": oracle.seed_profile(chart_seed_map),
            "source_loop_roots": {
                "A": [int(value) for value in c5_chart["loop"]["A"]],
                "B_over_Pstar": [int(value) for value in c5_chart["loop"]["B_over_Pstar"]],
            },
            "actual_slow_roots": {
                "A": [int(value) for value in q4_chart["radial"]["slow_A_y"]],
                "B_over_Pstar": [int(value) for value in q4_chart["radial"]["slow_B_y"]],
            },
            "source_reachable_node_count": len(source_reachable),
            "actual_reachable_node_count": len(actual_reachable),
            "seed_bindings": seed_receipt,
            "root_receipts": row_checks,
            "radius_relation": radius,
            "positive_J": positive_j,
            "held_phase_cache_zero": phase_zero,
            "J_independence_gate_passed": bool(positive_j["J_independence_gate_passed"]),
            "rule_coverage": {
                "source_operation_counts": source_ops,
                "actual_derivative_operation_counts": actual_ops,
                "required_rules": [
                    "sum", "product", "positive_quotient", "analytic_unary",
                    "squared_flat_sigma", "function_substitution_chain",
                    "definite_integral_partial_and_endpoints", "inverse_angle_chain",
                    "implicit_inverse", "lazy_flat_branch",
                ],
                "rule_checks": rule_checks,
                "rule_checks_passed": rule_checks_passed,
                "unsupported_reachable_operations": unsupported,
                "active_open_rules": sorted(oracle.open_rules),
                "unsupported_reachable_operations_clear": unsupported_clear,
                "active_open_rules_clear": open_rules_clear,
                "lazy_metadata_issues": list(oracle.lazy_metadata_issues),
                "binder_scope_issues": list(oracle.binder_scope_issues),
                "lazy_lemma_receipts": list(oracle.lazy_lemma_receipts),
                "smooth_control_remains_external": True,
            },
        }
        all_j_independence = all_j_independence and bool(positive_j["J_independence_gate_passed"])
        all_phase_zero = all_phase_zero and phase_zero
        all_unsupported_clear = all_unsupported_clear and unsupported_clear
        all_open_rules_clear = all_open_rules_clear and open_rules_clear
        all_rule_checks = all_rule_checks and rule_checks_passed
    conditional = bool(
        len(all_roots) == 36
        and all_roots_passed
        and all_radius
        and all_j
        and all_j_independence
        and all_phase_zero
        and all_unsupported_clear
        and all_open_rules_clear
        and all_rule_checks
        and not oracle.lazy_metadata_issues
        and not oracle.binder_scope_issues
        and q4_prefix
        and c5_prefix
        and pins_ok
        and parent_receipt["passed"]
        and rule_binding["all_rule_bindings_present"]
    )
    return {
        "passed": conditional,
        "conditional_passed": conditional,
        "all_passed": False,
        "parent_report_sha256": EXPECTED_PARENT_SHA256,
        "parent_report": PARENT_NAME,
        "module_sha256": _sha256(module_path),
        "physical_derivative_identity_proved": False,
        "seed_identity_discharged": False,
        "source_executed": False,
        "accepted_ancestor_execution": False,
        "source_family": q4.get("source_family"),
        "scope": "actual Q4/C5 held-fast-phase slow A/B physical-y roots for post_power, long_reshape and reference",
        "parent_validation": parent_receipt,
        "source_pins": pins,
        "source_pins_passed": pins_ok,
        "source_guard": rule_binding,
        "source_guard_snippets": [tuple(item) for item in rule_binding["source_guard_snippets"]],
        "gates": {
            "all_36_actual_roots": len(all_roots) == 36,
            "all_root_structural_checks": all_roots_passed,
            "all_radius_relations": all_radius,
            "all_inverse_J_bindings_and_independence": all_j,
            "all_J_independence": all_j_independence,
            "all_held_phase_cache_zero": all_phase_zero,
            "no_unsupported_reachable_operations": all_unsupported_clear,
            "no_active_open_rules": all_open_rules_clear,
            "all_required_rule_checks": all_rule_checks,
            "lazy_metadata_and_active_body_lemmas": not oracle.lazy_metadata_issues,
            "binder_scope_and_capture": not oracle.binder_scope_issues,
            "seed_and_smooth_external": True,
        },
        "lazy_lemma_receipts": list(oracle.lazy_lemma_receipts),
        "lazy_metadata_issues": list(oracle.lazy_metadata_issues),
        "binder_scope_issues": list(oracle.binder_scope_issues),
        "binder_scope_receipts": list(oracle.binder_scope_receipts),
        "graph_prefix": {
            "Q4_over_C5": q4_prefix,
            "C5_over_C6": c5_prefix,
            "Q4_nodes": len(nodes),
            "C5_nodes": len(c5_nodes),
            "C6_nodes": len(c6_nodes),
        },
        "charts": chart_receipts,
        "seed_binding_table": seed_table,
        "counts": {
            "Q4_graph_nodes": len(nodes),
            "C5_graph_nodes": len(c5_nodes),
            "C6_graph_nodes": len(c6_nodes),
            "scoped_charts": len(SCOPE_CHARTS),
            "slow_derivative_roots_total": len(all_roots),
            "slow_derivative_roots_passed": sum(1 for row in all_roots if row["passed"]),
            "hash_cons_dag_nodes": len(dag.nodes),
            "hash_cons_dag_unique_keys": len(dag.keys),
            "hash_cons_dag_source_leaf_nodes": sum(1 for node in dag.nodes if node[0] == "source_leaf"),
            "derivative_cache_entries": len(oracle._derivative_cache),
            "elapsed_seconds": time.monotonic() - started,
        },
        "limits": {
            "full_modified_Qy": False,
            "no_numeric_values": True,
            "accepted_execution": False,
            "seed_identity_external": True,
            "smooth_lazy_branch_control_external": True,
            "unscoped_Q4_charts": sorted(set(q4_charts) - set(SCOPE_CHARTS)),
        },
        "open_limits": [
            "physical seed equations D_native(field)=J*background_radial_row remain an upstream obligation",
            "smooth lazy-flat branch control remains an external analytic obligation",
            "full modified Q_y is outside this structural packet",
            "no numeric values or accepted ancestor execution were used",
        ],
    }


__all__ = ["GraphProofError", "OpenRuleError", "BudgetExceeded", "build"]
