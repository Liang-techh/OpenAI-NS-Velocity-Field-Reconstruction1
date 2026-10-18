"""Bounded exact source proof for the ``R_m`` pressure inlet seam.

The accepted Q4 and finite-``N`` patch snapshots share a graph prefix.  This
module consumes that already materialised graph and selected binding maps.  It
does not load a snapshot, import a producer, evaluate a numerical value, or
infer a physical germ identity.  The reducer follows only the graph rows
needed by the two ordinary ``Z`` orders and keeps source endpoint handles as
qualified atoms.

The source pressure formula is

``m15 = 1/2 S1 - R_Q_y(upper) + R_Q_y(initial)``.

At the patch lower endpoint the patch integral has equal bounds after the
guarded ``x -> 1`` substitution.  The remaining boundary endpoint is kept as
its own source atom, so the functional defect is recorded as
``current_RQy - patch_RQy``.  The same mechanism exposes the finite patch
increment to ``R_h`` without assigning a numeric value to the local integral.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from dataclasses import dataclass, field
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any


_OPERATIONS = {
    "exact_rational",
    "sum",
    "negative",
    "product",
    "positive_quotient",
    "quotient",
    "function_substitution",
    "definite_integral",
    "bound_variable",
}
_Z_GUARD = "Z_independent_substitution"


def canonical_digest(value: Any) -> str:
    """Return the canonical JSON SHA-256 digest of a JSON value."""

    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _module_sha256() -> str:
    digest = hashlib.sha256()
    with Path(__file__).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _handle(value: Any, path: str, node_count: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{path} must be an integer graph handle")
    if value < 0:
        raise ValueError(f"{path} must be nonnegative")
    if node_count is not None and value >= node_count:
        raise ValueError(f"{path} handle {value} is outside nodes")
    return value


def _fraction_text(value: Fraction) -> str:
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


@dataclass
class _Expr:
    """A constant plus a finite qualified source-atom linear combination."""

    constant: Fraction = Fraction(0)
    terms: dict[str, Fraction] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.constant = Fraction(self.constant)
        self.terms = {
            str(key): Fraction(value)
            for key, value in self.terms.items()
            if Fraction(value) != 0
        }

    @classmethod
    def rational(cls, value: Fraction | int) -> "_Expr":
        return cls(Fraction(value), {})

    @classmethod
    def atom(cls, name: str) -> "_Expr":
        return cls(Fraction(0), {name: Fraction(1)})

    @property
    def is_rational(self) -> bool:
        return not self.terms

    def scale(self, value: Fraction) -> "_Expr":
        value = Fraction(value)
        return _Expr(
            self.constant * value,
            {key: coefficient * value for key, coefficient in self.terms.items()},
        )

    def __add__(self, other: "_Expr") -> "_Expr":
        terms = dict(self.terms)
        for key, coefficient in other.terms.items():
            terms[key] = terms.get(key, Fraction(0)) + coefficient
        return _Expr(self.constant + other.constant, terms)

    def __neg__(self) -> "_Expr":
        return self.scale(Fraction(-1))

    def __sub__(self, other: "_Expr") -> "_Expr":
        return self + (-other)

    def key(self) -> tuple[Any, ...]:
        return (
            self.constant.numerator,
            self.constant.denominator,
            tuple(
                (key, value.numerator, value.denominator)
                for key, value in sorted(self.terms.items())
            ),
        )

    def json(self) -> dict[str, Any]:
        return {
            "constant": _fraction_text(self.constant),
            "terms": {
                key: _fraction_text(self.terms[key])
                for key in sorted(self.terms)
            },
        }


@dataclass(frozen=True)
class _EnvValue:
    value_handle: int | None
    expression: _Expr


@dataclass(frozen=True)
class _Alias:
    atom: str
    source_node: int
    descriptor: dict[str, Any]
    kind: str


def _path_value(value: Mapping[str, Any], paths: Sequence[Sequence[str]], label: str) -> Any:
    for path in paths:
        current: Any = value
        for key in path:
            if not isinstance(current, Mapping) or key not in current:
                break
            current = current[key]
        else:
            return current
    raise ValueError(f"{label} binding is missing")


def _unwrap_bindings(value: Any, names: Sequence[str]) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError("bindings must be a mapping")
    for name in names:
        nested = value.get(name)
        if isinstance(nested, Mapping):
            return nested
    if isinstance(value.get("exact_graph_nodes"), Sequence):
        for name in names:
            nested = value.get(name)
            if isinstance(nested, Mapping):
                return nested
    return value


def _rows(value: Any, label: str) -> list[int]:
    if not _is_sequence(value) or len(value) < 2:
        raise ValueError(f"{label} must provide at least ordinary Z orders 0 and 1")
    result: list[int] = []
    for order, item in enumerate(value[:2]):
        result.append(_handle(item, f"{label}[{order}]"))
    return result


@dataclass
class _Layout:
    modified_s1: list[int]
    modified_rqy: list[int]
    modified_m15: list[int]
    initial_rqy: list[int]
    patch_lower_s1: list[int]
    patch_rqy_rm: list[int]
    patch_m15_rm: list[int]
    patch_local_s1: list[int] | None
    patch_s1_rh: list[int] | None
    patch_rqy_rh: list[int] | None
    patch_m15_rh: list[int] | None

    @property
    def source_nodes(self) -> dict[str, list[int]]:
        result: dict[str, list[int]] = {
            "current_S1": list(self.modified_s1),
            "current_RQy": list(self.modified_rqy),
            "initial_RQy": list(self.initial_rqy),
            "current_m15": list(self.modified_m15),
            "patch_lower_S1": list(self.patch_lower_s1),
            "patch_RQy_Rm": list(self.patch_rqy_rm),
            "patch_m15_Rm": list(self.patch_m15_rm),
        }
        if self.patch_local_s1 is not None:
            result["patch_local_S1"] = list(self.patch_local_s1)
        if self.patch_s1_rh is not None:
            result["patch_S1_Rh"] = list(self.patch_s1_rh)
        if self.patch_rqy_rh is not None:
            result["patch_RQy_Rh"] = list(self.patch_rqy_rh)
        if self.patch_m15_rh is not None:
            result["patch_m15_Rh"] = list(self.patch_m15_rh)
        return result


def _optional_rows(value: Mapping[str, Any], paths: Sequence[Sequence[str]], label: str) -> list[int] | None:
    try:
        candidate = _path_value(value, paths, label)
    except ValueError:
        return None
    return _rows(candidate, label)


def _layout(modified_bindings: Mapping[str, Any], patch_bindings: Mapping[str, Any]) -> _Layout:
    modified = _unwrap_bindings(
        modified_bindings,
        ("modified_Q4_pressure_C3_bindings", "bindings"),
    )
    patch = _unwrap_bindings(
        patch_bindings,
        ("Rm_patch_finite_N_C5_bindings", "bindings"),
    )
    modified_rm = _path_value(
        modified,
        (("Rm",),),
        "modified Rm",
    )
    if not isinstance(modified_rm, Mapping):
        raise ValueError("modified Rm binding must be a mapping")
    patch_rm = _path_value(
        patch,
        (("endpoints", "Rm"), ("Rm",)),
        "patch Rm endpoint",
    )
    if not isinstance(patch_rm, Mapping):
        raise ValueError("patch Rm endpoint binding must be a mapping")
    modified_initial = _path_value(modified, (("initial",),), "modified initial")
    if not isinstance(modified_initial, Mapping):
        raise ValueError("modified initial binding must be a mapping")

    patch_local = _optional_rows(
        patch,
        (("local_integrals", "S1"),),
        "patch local S1",
    )
    patch_rh = patch.get("endpoints") if isinstance(patch.get("endpoints"), Mapping) else None
    patch_rh = patch_rh.get("Rh") if isinstance(patch_rh, Mapping) else None
    if patch_rh is not None and not isinstance(patch_rh, Mapping):
        raise ValueError("patch Rh endpoint binding must be a mapping")

    def required(mapping: Mapping[str, Any], path: Sequence[str], label: str) -> list[int]:
        return _rows(_path_value(mapping, (tuple(path),), label), label)

    return _Layout(
        modified_s1=required(modified_rm, ("integrals", "S1"), "modified Rm integrals.S1"),
        modified_rqy=required(modified_rm, ("boundary", "R_Q_y"), "modified Rm boundary.R_Q_y"),
        modified_m15=required(modified_rm, ("pressure", "m15"), "modified Rm pressure.m15"),
        initial_rqy=required(modified_initial, ("R_Q_y",), "modified initial.R_Q_y"),
        patch_lower_s1=required(patch_rm, ("integrals", "S1"), "patch endpoints.Rm.integrals.S1"),
        patch_rqy_rm=required(patch_rm, ("boundary", "R_Q_y"), "patch endpoints.Rm.boundary.R_Q_y"),
        patch_m15_rm=required(patch_rm, ("pressure", "m15"), "patch endpoints.Rm.pressure.m15"),
        patch_local_s1=patch_local,
        patch_s1_rh=(
            _rows(patch_rh["integrals"]["S1"], "patch endpoints.Rh.integrals.S1")
            if isinstance(patch_rh, Mapping)
            and isinstance(patch_rh.get("integrals"), Mapping)
            and "S1" in patch_rh["integrals"]
            else None
        ),
        patch_rqy_rh=(
            _rows(patch_rh["boundary"]["R_Q_y"], "patch endpoints.Rh.boundary.R_Q_y")
            if isinstance(patch_rh, Mapping)
            and isinstance(patch_rh.get("boundary"), Mapping)
            and "R_Q_y" in patch_rh["boundary"]
            else None
        ),
        patch_m15_rh=(
            _rows(patch_rh["pressure"]["m15"], "patch endpoints.Rh.pressure.m15")
            if isinstance(patch_rh, Mapping)
            and isinstance(patch_rh.get("pressure"), Mapping)
            and "m15" in patch_rh["pressure"]
            else None
        ),
    )


def _operand_handles(row: Mapping[str, Any], index: int) -> list[tuple[str, int]]:
    operation = row.get("operation")
    if operation not in _OPERATIONS:
        raise ValueError(
            f"visited node {index} uses unsupported operation {operation!r}"
        )
    fields: list[tuple[str, Any]] = []
    if operation in {"sum", "product"}:
        arguments = row.get("arguments")
        if not _is_sequence(arguments) or not arguments:
            raise ValueError(f"node {index}.{operation}.arguments must be nonempty")
        fields.extend((f"arguments[{offset}]", item) for offset, item in enumerate(arguments))
    elif operation == "negative":
        fields.append(("argument", row.get("argument")))
    elif operation in {"positive_quotient", "quotient"}:
        fields.extend(
            (("numerator", row.get("numerator")), ("denominator", row.get("denominator")))
        )
    elif operation == "function_substitution":
        fields.append(("expression", row.get("expression")))
        variable = row.get("variable")
        if isinstance(variable, int) and not isinstance(variable, bool):
            fields.append(("variable", variable))
        elif not isinstance(variable, str):
            raise ValueError(f"node {index}.variable must be a handle or string")
        fields.append(("value", row.get("value")))
    elif operation == "definite_integral":
        fields.extend(
            (
                ("integrand", row.get("integrand")),
                ("lower", row.get("lower")),
                ("upper", row.get("upper")),
            )
        )
        variable = row.get("variable")
        if not isinstance(variable, str) and not (
            isinstance(variable, int) and not isinstance(variable, bool)
        ):
            raise ValueError(f"node {index}.variable must be a handle or string")
        if isinstance(variable, int):
            fields.append(("variable", variable))
    elif operation == "exact_rational":
        numerator = row.get("numerator")
        denominator = row.get("denominator")
        if isinstance(numerator, bool) or not isinstance(numerator, int):
            raise ValueError(f"node {index}.numerator must be an integer literal")
        if isinstance(denominator, bool) or not isinstance(denominator, int):
            raise ValueError(f"node {index}.denominator must be an integer literal")
        if denominator == 0:
            raise ValueError(f"node {index} has zero exact-rational denominator")
    elif operation == "bound_variable":
        pass
    result: list[tuple[str, int]] = []
    for path, value in fields:
        result.append((path, _handle(value, f"node {index}.{path}")))
    return result


class _Reducer:
    """Conservative exact reducer for the selected source closure."""

    def __init__(self, nodes: Sequence[Mapping[str, Any]], prefix_count: int, layout: _Layout):
        self.nodes = nodes
        self.prefix_count = prefix_count
        self.layout = layout
        self.memo: dict[tuple[int, tuple[Any, ...]], _Expr] = {}
        self.used_nodes: set[int] = set()
        self.equal_bound_integrals: set[int] = set()
        self.alias_witnesses: list[dict[str, Any]] = []
        self._witness_keys: set[str] = set()
        self.protected: dict[int, _Alias] = {}
        self.aliases: dict[tuple[int, Any, int], _Alias] = {}
        self.local_rh_alias: dict[int, tuple[int, str, int, int]] = {}
        self._register_layout_aliases()

    def _row(self, index: int) -> Mapping[str, Any]:
        _handle(index, "source node", len(self.nodes))
        row = self.nodes[index]
        if not isinstance(row, Mapping):
            raise ValueError(f"nodes[{index}] must be a mapping")
        return row

    def _descriptor(self, index: int, row: Mapping[str, Any]) -> tuple[int, Any, int, dict[str, Any]]:
        if row.get("operation") != "function_substitution":
            raise ValueError(f"source endpoint node {index} is not function_substitution")
        if row.get(_Z_GUARD) is not True:
            raise ValueError(f"source endpoint node {index} lacks {_Z_GUARD}=True")
        expression = _handle(row.get("expression"), f"node {index}.expression", len(self.nodes))
        variable = row.get("variable")
        if isinstance(variable, bool) or not isinstance(variable, (int, str)):
            raise ValueError(f"node {index}.variable must be a handle or string")
        if isinstance(variable, int):
            _handle(variable, f"node {index}.variable", len(self.nodes))
            if variable >= index:
                raise ValueError(f"forward reference from node {index} variable to {variable}")
        value = _handle(row.get("value"), f"node {index}.value", len(self.nodes))
        if expression >= index:
            raise ValueError(f"forward reference from node {index} expression to {expression}")
        if value >= index:
            raise ValueError(f"forward reference from node {index} value to {value}")
        descriptor = {
            "expression": expression,
            "variable": variable,
            "value": value,
            _Z_GUARD: True,
        }
        return expression, variable, value, descriptor

    def _add_alias(
        self,
        root: int,
        atom: str,
        kind: str,
        *,
        protect: bool,
    ) -> None:
        row = self._row(root)
        if row.get("operation") == "function_substitution":
            expression, variable, value, descriptor = self._descriptor(root, row)
            key = (expression, variable, value)
            previous = self.aliases.get(key)
            if previous is not None and previous.atom != atom:
                raise ValueError(
                    "distinct protected source ports share one exact substitution descriptor"
                )
            info = _Alias(atom, root, descriptor, kind)
            self.aliases[key] = info
        else:
            if row.get("operation") not in {"atom", "source_atom", "qualified_source_atom"}:
                raise ValueError(
                    f"source endpoint node {root} must be a guarded substitution or qualified atom"
                )
            descriptor = {
                "node": root,
                "operation": row.get("operation"),
            }
            info = _Alias(atom, root, descriptor, kind)
        if protect:
            self.protected[root] = info

    def _register_layout_aliases(self) -> None:
        for order in range(2):
            self._add_alias(
                self.layout.modified_s1[order],
                f"S1[{order}]",
                "current_S1",
                protect=True,
            )
            self._add_alias(
                self.layout.modified_rqy[order],
                f"RQy[{order}]",
                "current_RQy",
                protect=True,
            )
            self._add_alias(
                self.layout.initial_rqy[order],
                f"initialRQy[{order}]",
                "initial_RQy",
                protect=True,
            )
            self._add_alias(
                self.layout.patch_rqy_rm[order],
                f"patchRQy[Rm,{order}]",
                "patch_RQy_Rm",
                protect=True,
            )
            if self.layout.patch_rqy_rh is not None:
                self._add_alias(
                    self.layout.patch_rqy_rh[order],
                    f"patchRQy[Rh,{order}]",
                    "patch_RQy_Rh",
                    protect=True,
                )

            if self.layout.patch_local_s1 is not None and self.layout.patch_s1_rh is not None:
                local = self.layout.patch_local_s1[order]
                rh_endpoint = self.layout.patch_s1_rh[order]
                endpoint_row = self._row(rh_endpoint)
                if endpoint_row.get("operation") != "function_substitution":
                    raise ValueError(f"patch Rh S1 endpoint {rh_endpoint} is not a substitution")
                _, variable, value, _ = self._descriptor(rh_endpoint, endpoint_row)
                if not isinstance(variable, int):
                    raise ValueError("patch Rh S1 endpoint variable must be a graph handle")
                self.local_rh_alias[local] = (
                    order,
                    f"patchLocalS1[Rh,{order}]",
                    variable,
                    value,
                )

    def _validate_operands(self, index: int, row: Mapping[str, Any]) -> list[tuple[str, int]]:
        operands = _operand_handles(row, index)
        for path, child in operands:
            if child >= index:
                raise ValueError(f"forward reference from node {index} {path} to {child}")
            if index < self.prefix_count and child >= self.prefix_count:
                raise ValueError(
                    f"prefix node {index} depends on post-patch node {child} through {path}"
                )
        return operands

    @staticmethod
    def _env_signature(env: Mapping[Any, _EnvValue]) -> tuple[Any, ...]:
        return tuple(
            (
                repr(key),
                value.value_handle,
                value.expression.key(),
            )
            for key, value in sorted(env.items(), key=lambda item: repr(item[0]))
        )

    def _record_witness(self, info: _Alias, *, via: str, expression_node: int, environment: Mapping[Any, _EnvValue]) -> None:
        payload = {
            "atom": info.atom,
            "kind": info.kind,
            "source_node": info.source_node,
            "expression_node": expression_node,
            "descriptor": deepcopy(info.descriptor),
            "via": via,
            "scoped_environment": {
                str(key): {
                    "value_handle": value.value_handle,
                    "value": value.expression.json(),
                }
                for key, value in sorted(environment.items(), key=lambda item: repr(item[0]))
            },
        }
        key = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        if key not in self._witness_keys:
            self._witness_keys.add(key)
            self.alias_witnesses.append(payload)

    def _scoped_alias(self, ref: int, env: Mapping[Any, _EnvValue]) -> _Alias | None:
        for variable, binding in sorted(env.items(), key=lambda item: repr(item[0])):
            if binding.value_handle is None:
                continue
            info = self.aliases.get((ref, variable, binding.value_handle))
            if info is not None:
                self._record_witness(
                    info,
                    via="scoped_function_substitution",
                    expression_node=ref,
                    environment=env,
                )
                return info
        return None

    def _rewrite_bound(self, ref: int, env: Mapping[Any, _EnvValue]) -> _Expr:
        for key, value in env.items():
            if key == ref:
                return value.expression
        raise ValueError(f"unbound variable source node {ref} was reached")

    def rewrite(self, ref: int, env: Mapping[Any, _EnvValue] | None = None) -> _Expr:
        env = {} if env is None else env
        _handle(ref, "source handle", len(self.nodes))
        key = (ref, self._env_signature(env))
        if key in self.memo:
            return self.memo[key]
        protected = self.protected.get(ref)
        if protected is not None:
            self.used_nodes.add(ref)
            self._record_witness(protected, via="protected_source_port", expression_node=ref, environment=env)
            result = _Expr.atom(protected.atom)
            self.memo[key] = result
            return result
        scoped = self._scoped_alias(ref, env)
        if scoped is not None:
            self.used_nodes.add(ref)
            result = _Expr.atom(scoped.atom)
            self.memo[key] = result
            return result

        row = self._row(ref)
        self.used_nodes.add(ref)
        operation = row.get("operation")
        self._validate_operands(ref, row)

        if operation == "exact_rational":
            result = _Expr.rational(Fraction(row["numerator"], row["denominator"]))
        elif operation == "bound_variable":
            result = self._rewrite_bound(ref, env)
        elif operation == "sum":
            result = _Expr.rational(0)
            for child in row["arguments"]:
                result = result + self.rewrite(child, env)
        elif operation == "negative":
            result = -self.rewrite(row["argument"], env)
        elif operation == "product":
            factors = [self.rewrite(child, env) for child in row["arguments"]]
            if any(factor.is_rational and factor.constant == 0 for factor in factors):
                result = _Expr.rational(0)
            else:
                symbolic = [factor for factor in factors if not factor.is_rational]
                if len(symbolic) > 1:
                    raise ValueError(
                        f"node {ref} product has more than one nonrational factor"
                    )
                coefficient = Fraction(1)
                for factor in factors:
                    if factor.is_rational:
                        coefficient *= factor.constant
                result = symbolic[0].scale(coefficient) if symbolic else _Expr.rational(coefficient)
        elif operation in {"positive_quotient", "quotient"}:
            numerator = self.rewrite(row["numerator"], env)
            denominator = self.rewrite(row["denominator"], env)
            if not denominator.is_rational:
                raise ValueError(f"node {ref} quotient denominator is not rational")
            if denominator.constant == 0:
                raise ValueError(f"node {ref} quotient denominator is zero")
            result = numerator.scale(Fraction(1, 1) / denominator.constant)
        elif operation == "function_substitution":
            if row.get(_Z_GUARD) is not True:
                raise ValueError(f"node {ref} lacks {_Z_GUARD}=True")
            expression = _handle(row["expression"], f"node {ref}.expression", len(self.nodes))
            variable = row["variable"]
            value_handle = _handle(row["value"], f"node {ref}.value", len(self.nodes))
            value_row = self._row(value_handle)
            # Endpoint descriptors may use a qualified exact source such as
            # ``exp(1)`` for R_h.  The descriptor handle is sufficient for
            # scoped alias matching; evaluating that source would cross the
            # bounded reducer's operation set.  Exact rationals and the
            # declared arithmetic rows remain fully reduced for toy graphs.
            if value_row.get("operation") == "analytic_unary":
                value = _Expr.atom(f"substitution_value[{value_handle}]")
            else:
                value = self.rewrite(value_handle, env)
            scoped_env = dict(env)
            scoped_env[variable] = _EnvValue(value_handle, value)
            result = self.rewrite(expression, scoped_env)
        elif operation == "definite_integral":
            lower = self.rewrite(row["lower"], env)
            upper = self.rewrite(row["upper"], env)
            if lower == upper:
                self.equal_bound_integrals.add(ref)
                result = _Expr.rational(0)
            elif ref in self.local_rh_alias:
                order, atom, variable, value_handle = self.local_rh_alias[ref]
                binding = env.get(variable)
                if binding is None or binding.value_handle != value_handle:
                    raise ValueError(
                        f"patch local S1 node {ref} lacks its exact Rh substitution descriptor"
                    )
                info = _Alias(
                    atom,
                    self.layout.patch_s1_rh[order] if self.layout.patch_s1_rh else ref,
                    {
                        "expression_node": ref,
                        "variable": variable,
                        "value": value_handle,
                        _Z_GUARD: True,
                    },
                    "patch_local_S1_Rh",
                )
                self._record_witness(info, via="matching_Rh_local_integral_substitution", expression_node=ref, environment=env)
                result = _Expr.atom(atom)
            else:
                raise ValueError(
                    f"node {ref} definite_integral has unequal bounds and is protected"
                )
        else:  # pragma: no cover - _validate_operands fails first
            raise ValueError(f"visited node {ref} uses unsupported operation {operation!r}")

        self.memo[key] = result
        return result


def _row_hashes(nodes: Sequence[Mapping[str, Any]], handles: Mapping[str, Sequence[int]]) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for label, values in handles.items():
        result[label] = [canonical_digest(nodes[index]) for index in values]
    return result


def _raw_terms(expression: _Expr) -> dict[str, str]:
    return {
        key: _fraction_text(expression.terms[key])
        for key in sorted(expression.terms)
    }


def _abstract_terms(expression: _Expr, mapping: Mapping[str, str]) -> dict[str, str]:
    result: dict[str, Fraction] = {}
    for atom, coefficient in expression.terms.items():
        result[mapping.get(atom, atom)] = result.get(mapping.get(atom, atom), Fraction(0)) + coefficient
    return {key: _fraction_text(result[key]) for key in sorted(result) if result[key]}


def _expression_text(expression: _Expr) -> str:
    terms: list[str] = []
    if expression.constant:
        terms.append(_fraction_text(expression.constant))
    for atom in sorted(expression.terms):
        coefficient = expression.terms[atom]
        magnitude = abs(coefficient)
        factor = atom if magnitude == 1 else f"{_fraction_text(magnitude)}*{atom}"
        if not terms:
            terms.append(f"-{factor}" if coefficient < 0 else factor)
        else:
            terms.append(("-" if coefficient < 0 else "+") + factor)
    return "0" if not terms else " ".join(terms)


def _require_expression(expression: _Expr, expected: Mapping[str, Fraction], label: str) -> None:
    expected_expr = _Expr(Fraction(0), dict(expected))
    if expression != expected_expr:
        raise ValueError(
            f"{label} did not reduce to the declared linear source form: "
            f"got {_expression_text(expression)}"
        )


def _order_receipt(
    order: int,
    *,
    modified_node: int,
    patch_node: int,
    modified: _Expr,
    patch: _Expr,
    modified_atoms: Mapping[str, str],
    patch_atoms: Mapping[str, str],
) -> dict[str, Any]:
    return {
        "ordinary_Z_order": order,
        "modified_m15_node": modified_node,
        "patch_m15_node": patch_node,
        "modified": {
            "expression": modified.json(),
            "expression_text": _expression_text(modified),
            "terms": _abstract_terms(modified, modified_atoms),
            "raw_terms": _raw_terms(modified),
        },
        "patch": {
            "expression": patch.json(),
            "expression_text": _expression_text(patch),
            "terms": _abstract_terms(patch, patch_atoms),
            "raw_terms": _raw_terms(patch),
        },
    }


def prove(
    nodes: Sequence[Mapping[str, Any]],
    modified_bindings: Mapping[str, Any],
    patch_bindings: Mapping[str, Any],
    *,
    expected_nodes_digest: str,
    common_prefix_count: int = 183700,
) -> dict[str, Any]:
    """Prove the bounded exact ``R_m`` pressure inlet relation.

    Only ordinary orders ``Z^0`` and ``Z^1`` are traversed.  The caller owns
    authentication of the saved Q4 and finite-``N`` prefixes; this function
    authenticates the supplied complete graph against ``expected_nodes_digest``
    and validates backward references encountered by the reducer.
    """

    if not _is_sequence(nodes) or not nodes:
        raise ValueError("nodes must be a nonempty sequence")
    graph = [deepcopy(dict(row)) if isinstance(row, Mapping) else row for row in nodes]
    for index, row in enumerate(graph):
        if not isinstance(row, Mapping):
            raise ValueError(f"nodes[{index}] must be a mapping")
    actual_digest = canonical_digest(graph)
    if not isinstance(expected_nodes_digest, str) or expected_nodes_digest != actual_digest:
        raise ValueError("nodes canonical digest does not match expected_nodes_digest")
    if isinstance(common_prefix_count, bool) or not isinstance(common_prefix_count, int):
        raise ValueError("common_prefix_count must be an integer")
    if common_prefix_count < 1 or common_prefix_count > len(graph):
        raise ValueError("common_prefix_count must lie within nodes")

    layout = _layout(modified_bindings, patch_bindings)
    rh_parts = (
        layout.patch_local_s1,
        layout.patch_s1_rh,
        layout.patch_rqy_rh,
        layout.patch_m15_rh,
    )
    rh_present = [part is not None for part in rh_parts]
    if any(rh_present) and not all(rh_present):
        raise ValueError(
            "patch Rh extension requires local_integrals.S1 and all endpoints.Rh "
            "integrals.S1, boundary.R_Q_y, and pressure.m15 bindings"
        )
    for label, handles in layout.source_nodes.items():
        for order, handle in enumerate(handles):
            _handle(handle, f"{label}[{order}]", len(graph))
    for label, handles in {
        "modified_S1": layout.modified_s1,
        "modified_RQy": layout.modified_rqy,
        "modified_m15": layout.modified_m15,
        "initial_RQy": layout.initial_rqy,
    }.items():
        if any(handle >= common_prefix_count for handle in handles):
            raise ValueError(f"{label} must lie in the authenticated common prefix")

    reducer = _Reducer(graph, common_prefix_count, layout)
    modified_forms: list[_Expr] = []
    patch_forms: list[_Expr] = []
    lower_forms: list[_Expr] = []
    rh_s1_forms: list[_Expr] = []
    rh_forms: list[_Expr] = []

    for order in range(2):
        modified = reducer.rewrite(layout.modified_m15[order])
        patch = reducer.rewrite(layout.patch_m15_rm[order])
        lower = reducer.rewrite(layout.patch_lower_s1[order])
        expected_modified = {
            f"S1[{order}]": Fraction(1, 2),
            f"RQy[{order}]": Fraction(-1),
            f"initialRQy[{order}]": Fraction(1),
        }
        expected_patch = {
            f"S1[{order}]": Fraction(1, 2),
            f"patchRQy[Rm,{order}]": Fraction(-1),
            f"initialRQy[{order}]": Fraction(1),
        }
        _require_expression(modified, expected_modified, f"modified m15 order {order}")
        _require_expression(patch, expected_patch, f"patch m15 order {order}")
        _require_expression(lower, {f"S1[{order}]": Fraction(1)}, f"patch lower S1 order {order}")
        modified_forms.append(modified)
        patch_forms.append(patch)
        lower_forms.append(lower)

        if all(rh_present):
            rh = reducer.rewrite(layout.patch_m15_rh[order])
            expected_rh = {
                f"S1[{order}]": Fraction(1, 2),
                f"patchLocalS1[Rh,{order}]": Fraction(1, 2),
                f"patchRQy[Rh,{order}]": Fraction(-1),
                f"initialRQy[{order}]": Fraction(1),
            }
            _require_expression(rh, expected_rh, f"patch Rh m15 order {order}")
            rh_forms.append(rh)
            if all(rh_present):
                rh_s1 = reducer.rewrite(layout.patch_s1_rh[order])
                _require_expression(
                    rh_s1,
                    {
                        f"S1[{order}]": Fraction(1),
                        f"patchLocalS1[Rh,{order}]": Fraction(1),
                    },
                    f"patch Rh S1 order {order}",
                )
                rh_s1_forms.append(rh_s1)

    defect_forms = [patch - modified for patch, modified in zip(patch_forms, modified_forms)]
    increment_forms = [rh - patch for rh, patch in zip(rh_forms, patch_forms)]
    if all(rh_present) and len(increment_forms) != 2:
        raise ValueError("patch Rh extension did not produce both ordinary Z orders")
    for order, defect in enumerate(defect_forms):
        _require_expression(
            defect,
            {f"RQy[{order}]": Fraction(1), f"patchRQy[Rm,{order}]": Fraction(-1)},
            f"Rm pressure defect order {order}",
        )
    for order, increment in enumerate(increment_forms):
        _require_expression(
            increment,
            {
                f"patchLocalS1[Rh,{order}]": Fraction(1, 2),
                f"patchRQy[Rh,{order}]": Fraction(-1),
                f"patchRQy[Rm,{order}]": Fraction(1),
            },
            f"patch Rh minus Rm increment order {order}",
        )

    # Keep the concise formula names alongside fully qualified raw atoms.
    modified_atom_names = {
        atom: ("S1" if atom.startswith("S1[") else "RQy" if atom.startswith("RQy[") else "initialRQy")
        for order in range(2)
        for atom in (
            f"S1[{order}]",
            f"RQy[{order}]",
            f"initialRQy[{order}]",
        )
    }
    patch_atom_names = {
        atom: ("S1" if atom.startswith("S1[") else "patchRQy" if atom.startswith("patchRQy[Rm") else "initialRQy")
        for order in range(2)
        for atom in (
            f"S1[{order}]",
            f"patchRQy[Rm,{order}]",
            f"initialRQy[{order}]",
        )
    }

    rows_used = sorted(reducer.used_nodes)
    certificate: dict[str, Any] = {
        "passed": True,
        "source_only": True,
        "source_executed": False,
        "accepted_ancestors_replayed": False,
        "pressure_component_only": True,
        "pressure_component": "Phi15=1/2*S1-(RQy-RQy_initial)",
        "full_m15_or_velocity_germ_join": False,
        "old_212_and_other_n1_sectors_untouched": True,
        "nodes_digest": actual_digest,
        "expected_nodes_digest": expected_nodes_digest,
        "common_prefix_count": common_prefix_count,
        "common_prefix_digest": canonical_digest(graph[:common_prefix_count]),
        "ordinary_Z_orders": [0, 1],
        "m15": [
            _order_receipt(
                order,
                modified_node=layout.modified_m15[order],
                patch_node=layout.patch_m15_rm[order],
                modified=modified_forms[order],
                patch=patch_forms[order],
                modified_atoms=modified_atom_names,
                patch_atoms=patch_atom_names,
            )
            for order in range(2)
        ],
        "patch_lower_integrals_S1": [
            {
                "ordinary_Z_order": order,
                "endpoint_node": layout.patch_lower_s1[order],
                "reduced": lower_forms[order].json(),
                "reduced_text": _expression_text(lower_forms[order]),
                "same_current_source_S1": True,
            }
            for order in range(2)
        ],
        "functional_defect": [
            {
                "ordinary_Z_order": order,
                "expression": "patch_pressure_at_Rm-current_pressure_at_Rm",
                "reduced": defect_forms[order].json(),
                "reduced_text": _expression_text(defect_forms[order]),
                "terms": {"current_RQy": "1", "patch_RQy": "-1"},
                "raw_terms": _raw_terms(defect_forms[order]),
                "fixed_sign": "current_RQy-patch_RQy",
            }
            for order in range(2)
        ],
        "patch_Rh_extension_proved": all(rh_present),
        "patch_Rh_increment": (
            [
                {
                    "ordinary_Z_order": order,
                    "expression": "patch_pressure_at_Rh-patch_pressure_at_Rm",
                    "reduced": increment_forms[order].json(),
                    "reduced_text": _expression_text(increment_forms[order]),
                    "terms": {
                        "localS1_Rh": "1/2",
                        "patch_RQy_Rh": "-1",
                        "patch_RQy_Rm": "1",
                    },
                    "raw_terms": _raw_terms(increment_forms[order]),
                    "current_S1_and_initial_RQy_cancelled": True,
                }
                for order in range(2)
            ]
            if all(rh_present)
            else []
        ),
        "source_nodes_used": layout.source_nodes,
        "traversed_nodes": rows_used,
        "equal_bound_integral_reductions": [
            {
                "node": node,
                "reduced": _Expr.rational(0).json(),
            }
            for node in sorted(reducer.equal_bound_integrals)
        ],
        "source_node_row_hashes": _row_hashes(
            graph,
            {
                key: values
                for key, values in layout.source_nodes.items()
            },
        ),
        "alias_witnesses": sorted(
            reducer.alias_witnesses,
            key=lambda item: (
                item["atom"],
                item["source_node"],
                item["expression_node"],
                item["via"],
            ),
        ),
        "equal_bound_integral_ids": sorted(reducer.equal_bound_integrals),
        "protected_source_atoms": sorted(
            {
                info.atom
                for info in reducer.protected.values()
            }
        ),
        "graph_immutable": True,
        "metadata_preserved": True,
        "literal_integer_metadata_not_rewritten": True,
        "backward_DAG_validated_for_traversed_operands": True,
        "prefix_cannot_depend_on_patch_coordinate": True,
        "physical_pressure_sign": {
            "formula": "1/2*S1-R_Q_y(upper)+R_Q_y(initial)",
            "upper_RQy_coefficient": "-1",
            "initial_RQy_coefficient": "1",
            "ordinary_Z1": True,
            "factorial_applied": False,
            "ordinary_Z_normalization": "ordinaryZ1",
        },
        "RQy_identity_proved": False,
        "full_germ_join": False,
        "distinct_boundary_atoms_preserved": True,
        "no_RQy_cancellation_across_distinct_handles": True,
        "numeric_caps_or_values_evaluated": False,
        "all212_or_error_fields_touched": False,
        "own_swirl_added": False,
        "implementation_sha256": _module_sha256(),
        "tool_model_metadata": None,
        "tool_model_metadata_present": False,
    }
    return certificate


__all__ = ["canonical_digest", "prove"]
