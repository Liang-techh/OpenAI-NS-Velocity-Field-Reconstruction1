"""Standalone Decimal point leaf for the original flat cutoff, revision 3.

This module is deliberately a leaf.  It computes only the scalar flat
cutoff and its ordinary first derivative at one supplied ``Decimal`` point;
it never evaluates an interval, a graph, or a Navier--Stokes expression.
The binding function consumes a prior v1 compiled receipt as opaque source
records and exposes the two new operations as v4 callables.  Matching a caller
digest to a saved record is provenance bookkeeping, not a claim that this
implementation has been re-derived from, or numerically certified against,
the historical source.  The saved AST hash strings in ``_canonical_guard``
are assertions retained by the receipt; they are not recomputed here.
"""

from __future__ import annotations

from collections.abc import Mapping
from decimal import (
    Decimal,
    DecimalException,
    DivisionByZero,
    InvalidOperation,
    Overflow,
    Subnormal,
    Underflow,
    localcontext,
)
from typing import Any


__all__ = [
    "SCHEMA",
    "COMPONENT",
    "RESULT_SCHEMA",
    "PRIOR_SCHEMA",
    "SIGMA_NODE_ID",
    "SIGMA_PRIME_NODE_ID",
    "ARGUMENT_NODE_ID",
    "SIGMA_DEFINITION",
    "SIGMA_PRIME_ORDINARY_DERIVATIVE",
    "ARGUMENT_VARIABLE_NAME",
    "FlatSigmaError",
    "InvalidPointInput",
    "PrecisionInsufficientError",
    "SourceBindingError",
    "PrecisionError",
    "BindingError",
    "SigmaPair",
    "stable_sigma",
    "stable_sigma_point",
    "evaluate_point",
    "evaluate",
    "flat_sigma",
    "flat_sigma_prime",
    "stable_sigma_prime",
    "bind_flat_sigma_point_kernel",
    "bind_point_kernel",
    "bind_prior_program",
    "bind",
]


SCHEMA = "n1_flat_sigma_point_kernel_v4"
COMPONENT = SCHEMA
RESULT_SCHEMA = "n1_flat_sigma_point_kernel_v4_result"
# The v4 execution leaf consumes the already compiled v1 receipt format.  A
# new leaf/result schema does not authorize changing the archived input.
PRIOR_SCHEMA = "n1_raw_g5_dag_program_v1"

# These are the exact receipt node IDs.  The v4 distinction is the new leaf
# and result component, not a rewritten archived graph namespace.
SIGMA_NODE_ID = 35949
SIGMA_PRIME_NODE_ID = 35953
ARGUMENT_NODE_ID = 35848

SIGMA_DEFINITION = (
    "0 for x<=0; 1 for x>=1; "
    "exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside"
)
SIGMA_PRIME_ORDINARY_DERIVATIVE = "same original sigma Taylor coefficient[1]"
ARGUMENT_VARIABLE_NAME = "original_lower_collar_s"

# Canonical source details are saved receipt assertions.  They are checked
# when a real receipt names one of these modules, but the AST text is not
# recomputed by this leaf and never turns the result into a certificate.
STABLE_SIGMA_MODULE = "lei_ren_part1_paper_compliant_outer_initial.py"
STABLE_SIGMA_FUNCTION = "stable_sigma"
STABLE_SIGMA_MODULE_SHA256 = (
    "a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657"
)
STABLE_SIGMA_AST_SHA256_SAVED_ASSERTION = (
    "4eca34dd2d6bd23b5bff76496682bece54924d3145ac99e95ed9066d9572ab3c"
)
STABLE_SIGMA_AST_SHA256 = STABLE_SIGMA_AST_SHA256_SAVED_ASSERTION
SIGMA_JETS_MODULE = "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"
SIGMA_JETS_FUNCTION = "sigma_jets"
SIGMA_JETS_MODULE_SHA256 = (
    "32f51123a47c2bafbb6c2f7056302bab996a019883615333875207d72a9bb7d7"
)
SIGMA_JETS_AST_SHA256_SAVED_ASSERTION = (
    "a63c47d05049c86007b530f3095c515d724bef53092c5f7867d53e5bbaa06236"
)
SIGMA_JETS_AST_SHA256 = SIGMA_JETS_AST_SHA256_SAVED_ASSERTION
AST_ASSERTION_STATUS = "saved_assertion_not_recomputed"

_MISSING = object()
_MAX_CONTEXT_EXPONENT = 999_999_999


class FlatSigmaError(ValueError):
    """Base class for this fail-closed point leaf."""


class InvalidPointInput(FlatSigmaError):
    """The point or requested Decimal precision is malformed."""


class PrecisionInsufficientError(FlatSigmaError):
    """A finite interior result would be lost at the requested precision."""


class SourceBindingError(FlatSigmaError):
    """The v4 receipt or caller-supplied source digest is inconsistent."""


PrecisionError = PrecisionInsufficientError
BindingError = SourceBindingError


class SigmaPair(tuple):
    """Tuple-compatible ``(sigma, sigma_prime)`` with point metadata."""

    def __new__(
        cls, sigma: Decimal, sigma_prime: Decimal, precision: int | None = None
    ) -> "SigmaPair":
        value = tuple.__new__(cls, (sigma, sigma_prime))
        value.precision = precision
        value.certified = False
        value.point_evaluation_only = True
        return value

    def __getitem__(self, key: Any) -> Any:
        if isinstance(key, str):
            fields = {
                "sigma": tuple.__getitem__(self, 0),
                "value": tuple.__getitem__(self, 0),
                "derivative": tuple.__getitem__(self, 1),
                "sigma_prime": tuple.__getitem__(self, 1),
                "precision": self.precision,
                "certified": False,
                "point_evaluation_only": True,
            }
            if key not in fields:
                raise KeyError(key)
            return fields[key]
        return tuple.__getitem__(self, key)

    def get(self, key: str, default: Any = None) -> Any:
        try:
            return self[key]
        except KeyError:
            return default

    def as_dict(self, precision: int | None = None) -> dict[str, Any]:
        if precision is None:
            precision = self.precision
        return {
            "sigma": self[0],
            "sigma_prime": self[1],
            "value": self[0],
            "derivative": self[1],
            "precision": precision,
            "certified": False,
            "point_evaluation_only": True,
        }


def _decimal(value: Any, label: str) -> Decimal:
    if type(value) is not Decimal:
        raise TypeError(f"{label} must be a decimal.Decimal")
    if not value.is_finite():
        raise InvalidPointInput(f"{label} must be finite")
    return value


def _precision(value: Any) -> int:
    if type(value) is not int or value <= 0:
        raise InvalidPointInput("precision must be a positive integer")
    return value


def _working_precision(point: Decimal, requested: int) -> int:
    # The extra digits cover both reciprocal squares and the final logistic
    # subtraction while retaining every digit supplied in the point.
    supplied = len(point.as_tuple().digits)
    return max(72, requested + 40, requested + supplied + 12)


def _precision_failure(message: str, cause: BaseException | None = None) -> None:
    error = PrecisionInsufficientError(message)
    if cause is None:
        raise error
    raise error from cause


def _sigma_pair(point: Decimal, requested: int) -> SigmaPair:
    zero = Decimal(0)
    one = Decimal(1)
    if point <= zero:
        return SigmaPair(zero, zero, requested)
    if point >= one:
        return SigmaPair(one, zero, requested)

    try:
        with localcontext() as work:
            work.prec = _working_precision(point, requested)
            work.Emax = _MAX_CONTEXT_EXPONENT
            work.Emin = -_MAX_CONTEXT_EXPONENT
            for signal in (
                DivisionByZero,
                InvalidOperation,
                Overflow,
                Underflow,
                Subnormal,
            ):
                work.traps[signal] = True

            complement = +(one - point)
            if complement <= zero:
                _precision_failure("interior complement rounded to zero")

            # Compute on the left half, so the exponential argument is never
            # positive.  Reflecting the result gives the right half exactly
            # as a point value, not as a selected interval endpoint.
            left_half = point <= complement
            small = point if left_half else complement
            large = +(one - small)
            if small <= zero or large <= zero:
                _precision_failure("interior point rounded to an endpoint")

            inv_small = +(one / small)
            inv_large = +(one / large)
            small_sq = +(inv_small * inv_small)
            large_sq = +(inv_large * inv_large)
            log_odds = +(large_sq - small_sq)
            if log_odds > zero:
                _precision_failure("logistic exponent became positive")
            tail = +log_odds.exp()
            if not tail.is_finite() or tail <= zero:
                _precision_failure("interior logistic tail underflowed")
            denominator = +(one + tail)
            if not denominator.is_finite() or denominator <= zero:
                _precision_failure("interior logistic denominator is invalid")

            left_value = +(tail / denominator)
            q = +(tail / (denominator * denominator))
            if left_value <= zero or left_value >= one or q <= zero:
                _precision_failure("interior logistic value lost precision")
            value_raw = left_value if left_half else +(one - left_value)
            factor = +(Decimal(2) * (small_sq * inv_small + large_sq * inv_large))
            derivative_raw = +(q * factor)

            with localcontext() as output:
                output.prec = requested
                output.Emax = _MAX_CONTEXT_EXPONENT
                output.Emin = -_MAX_CONTEXT_EXPONENT
                for signal in (InvalidOperation, Overflow, Underflow, Subnormal):
                    output.traps[signal] = True
                value = +value_raw
                derivative = +derivative_raw

            if not value.is_finite() or not derivative.is_finite():
                _precision_failure("interior result is non-finite")
            if value <= zero or value >= one:
                _precision_failure("interior value rounded to an endpoint")
            if derivative <= zero:
                _precision_failure("interior derivative rounded to zero")
            return SigmaPair(value, derivative, requested)
    except PrecisionInsufficientError:
        raise
    except DecimalException as cause:
        _precision_failure("Decimal overflow, underflow, or invalid operation", cause)
    except (ArithmeticError, ValueError) as cause:
        _precision_failure("Decimal point evaluation failed", cause)


def stable_sigma(x: Decimal, precision: int) -> SigmaPair:
    """Return the original flat cutoff and its ordinary point derivative."""

    return _sigma_pair(_decimal(x, "x"), _precision(precision))


def stable_sigma_point(x: Decimal, precision: int) -> dict[str, Any]:
    point = _decimal(x, "x")
    requested = _precision(precision)
    pair = _sigma_pair(point, requested)
    return pair.as_dict(requested) | {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "status": "point_evaluation",
        "evaluation_kind": "point",
        "x": point,
        "interval": False,
        "interval_value": False,
        "full_Navier_Stokes": False,
        "claims": {
            "certified": False,
            "point_evaluation_only": True,
            "interval": False,
            "full_Navier_Stokes": False,
        },
    }


evaluate_point = stable_sigma_point
evaluate = stable_sigma_point


def flat_sigma(x: Decimal, precision: int) -> Decimal:
    return stable_sigma(x, precision)[0]


def flat_sigma_prime(x: Decimal, precision: int) -> Decimal:
    return stable_sigma(x, precision)[1]


stable_sigma_prime = flat_sigma_prime


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise SourceBindingError(f"{label} must be a mapping")
    return value


def _node_id(value: Any, label: str) -> int:
    if type(value) is not int or value < 0:
        raise SourceBindingError(f"{label} must be a nonnegative integer")
    return value


def _raw_record(row: Any, label: str) -> Mapping[str, Any]:
    item = _mapping(row, label)
    original = item.get("original_record", _MISSING)
    record = item.get("record", _MISSING)
    if original is not _MISSING:
        original = _mapping(original, f"{label}.original_record")
        if record is not _MISSING:
            record = _mapping(record, f"{label}.record")
            if dict(original) != dict(record):
                raise SourceBindingError(f"{label} record and original_record disagree")
        return original
    if record is not _MISSING:
        return _mapping(record, f"{label}.record")
    if "operation" in item:
        return item
    raise SourceBindingError(f"{label} has no exact original source record")


def _program_candidate(prior: Any) -> Mapping[str, Any]:
    supplied = _mapping(prior, "prior compiled receipt")
    if supplied.get("schema") == PRIOR_SCHEMA:
        candidate = supplied
    else:
        nested = supplied.get("program", _MISSING)
        if isinstance(nested, Mapping) and nested.get("schema") == PRIOR_SCHEMA:
            candidate = nested
        else:
            raise SourceBindingError(f"prior mapping must have schema {PRIOR_SCHEMA!r}")

    # Never infer compiled status from a missing field.  If both fields are
    # present they must agree; a receipt that says failed/opaque is rejected
    # even if its other status field happens to say compiled.
    status = candidate.get("status", _MISSING)
    compile_status = candidate.get("compile_status", _MISSING)
    if status is _MISSING and compile_status is _MISSING:
        raise SourceBindingError("prior receipt compiled status is missing")
    if status is not _MISSING and status != "compiled":
        raise SourceBindingError("prior receipt status is not compiled")
    if compile_status is not _MISSING and compile_status != "compiled":
        raise SourceBindingError("prior receipt compile_status is not compiled")
    return candidate


def _collect_records(program: Mapping[str, Any]) -> dict[int, tuple[Any, Mapping[str, Any]]]:
    result: dict[int, tuple[Any, Mapping[str, Any]]] = {}
    wrapper_nodes: set[int] = set()
    target_nodes = (SIGMA_NODE_ID, SIGMA_PRIME_NODE_ID, ARGUMENT_NODE_ID)
    opaque_nodes = (SIGMA_NODE_ID, SIGMA_PRIME_NODE_ID)

    def add(node: int, row: Any, label: str, *, wrapper: bool) -> None:
        # The archived receipt can contain many unrelated wrappers.  Identify
        # their node ID only; do not read their records or operation metadata.
        if node not in target_nodes:
            return
        raw = _raw_record(row, label)
        if wrapper and node in opaque_nodes:
            operation = raw.get("operation", _MISSING)
            if not isinstance(operation, str):
                raise SourceBindingError(f"{label} wrapper operation is missing")
            # Validate every target sigma wrapper before selecting a row.
            _opaque(row, node, operation)
        elif not wrapper and node in opaque_nodes:
            # A target bare record is not a wrapper, but fail closed if it
            # carries wrapper metadata that falsely claims support.
            operation = raw.get("operation", _MISSING)
            if isinstance(operation, str):
                _opaque(row, node, operation)
        if node in result and dict(result[node][1]) != dict(raw):
            raise SourceBindingError(f"node {node} has conflicting records")
        if wrapper:
            if node not in wrapper_nodes:
                result[node] = (row, raw)
                wrapper_nodes.add(node)
        elif node not in result:
            # A bare record is a fallback only; it never replaces a wrapper.
            result[node] = (row, raw)

    instructions = program.get("instructions", _MISSING)
    if instructions is not _MISSING:
        if not isinstance(instructions, (list, tuple)):
            raise SourceBindingError("prior instructions must be a list or tuple")
        for index, row in enumerate(instructions):
            item = _mapping(row, f"instructions[{index}]")
            add(
                _node_id(item.get("node_id"), f"instructions[{index}].node_id"),
                row,
                f"instructions[{index}]",
                wrapper=True,
            )

    nodes = program.get("nodes", _MISSING)
    if nodes is not _MISSING:
        for key, row in _mapping(nodes, "prior nodes").items():
            if type(key) is int:
                node = key
            elif type(key) is str and key.isdigit():
                node = int(key)
            else:
                continue
            add(node, row, f"nodes[{key!r}]", wrapper=True)

    records = program.get("records", _MISSING)
    if records is not _MISSING:
        records_map = _mapping(records, "prior records")
        for key, row in records_map.items():
            if type(key) is int:
                node = key
            elif type(key) is str and key.isdigit():
                node = int(key)
            else:
                continue
            add(node, row, f"records[{key!r}]", wrapper=False)
    if not result:
        raise SourceBindingError("prior receipt has no records")
    return result

def _opaque(row: Any, node: int, operation: str) -> None:
    if not isinstance(row, Mapping):
        return
    if row.get("status", _MISSING) == "supported":
        raise SourceBindingError(f"node {node} falsely claims {operation!r} is supported")
    op = row.get("op", _MISSING)
    if op is not _MISSING and op != "opaque":
        raise SourceBindingError(f"node {node} op marker is not opaque")
    opaque = row.get("opaque_operation", _MISSING)
    if opaque is not _MISSING and opaque != operation:
        raise SourceBindingError(f"node {node} opaque operation disagrees")


def _operation(record: Mapping[str, Any], node: int) -> str:
    operation = record.get("operation", _MISSING)
    if not isinstance(operation, str):
        raise SourceBindingError(f"node {node} operation is missing")
    return operation


def _function(record: Mapping[str, Any], node: int) -> str:
    names = [record[key] for key in ("function", "name") if key in record]
    if not names or any(not isinstance(item, str) for item in names) or len(set(names)) != 1:
        raise SourceBindingError(f"node {node} function/name identity is malformed")
    return names[0]


def _source(record: Mapping[str, Any], node: int) -> Mapping[str, Any]:
    if node == SIGMA_NODE_ID:
        value = record.get("source_binding", _MISSING)
    else:
        derivative = record.get("derivative_source", _MISSING)
        binding = record.get("source_binding", _MISSING)
        if derivative is not _MISSING and binding is not _MISSING:
            derivative_map = _mapping(derivative, f"node {node} derivative source")
            binding_map = _mapping(binding, f"node {node} source binding")
            if dict(derivative_map) != dict(binding_map):
                raise SourceBindingError(
                    f"node {node} derivative_source and source_binding disagree"
                )
            value = derivative_map
        elif derivative is not _MISSING:
            value = derivative
        else:
            value = binding
    return _mapping(value, f"node {node} source binding")


def _field(source: Mapping[str, Any], names: tuple[str, ...], label: str) -> str:
    values: list[str] = []
    for name in names:
        if name not in source:
            continue
        value = source[name]
        if not isinstance(value, str) or not value:
            raise SourceBindingError(f"{label} alias {name!r} is malformed")
        values.append(value)
    if not values:
        raise SourceBindingError(f"{label} is required")
    if len(set(values)) != 1:
        raise SourceBindingError(f"{label} aliases disagree")
    return values[0]


def _identity(source: Mapping[str, Any], label: str) -> dict[str, str]:
    identity = {
        "module": _field(source, ("module", "source_module"), f"{label}.module"),
        "module_sha256": _field(
            source,
            ("module_sha256", "module_hash", "source_sha256"),
            f"{label}.module_sha256",
        ),
        "AST_sha256": _field(
            source,
            ("AST_sha256", "ast_sha256", "AST_hash", "ast_hash"),
            f"{label}.AST_sha256",
        ),
    }
    if "function" in source:
        if not isinstance(source["function"], str) or not source["function"]:
            raise SourceBindingError(f"{label}.function is malformed")
        identity["function"] = source["function"]
    return identity


def _canonical_guard(identity: Mapping[str, str], label: str) -> None:
    """Compare saved canonical strings; do not call a source or recompute AST."""

    module = identity["module"]
    if module == STABLE_SIGMA_MODULE:
        expected = {
            "function": STABLE_SIGMA_FUNCTION,
            "module_sha256": STABLE_SIGMA_MODULE_SHA256,
            "AST_sha256": STABLE_SIGMA_AST_SHA256_SAVED_ASSERTION,
        }
    elif module == SIGMA_JETS_MODULE:
        expected = {
            "function": SIGMA_JETS_FUNCTION,
            "module_sha256": SIGMA_JETS_MODULE_SHA256,
            "AST_sha256": SIGMA_JETS_AST_SHA256_SAVED_ASSERTION,
        }
    else:
        return
    for key, value in expected.items():
        if identity.get(key) != value:
            raise SourceBindingError(f"{label}: saved canonical assertion mismatch for {key}")


def _digest_entry(source_digests: Any, role: str, node: int) -> Mapping[str, Any]:
    if isinstance(source_digests, (list, tuple)):
        matches = [
            item
            for item in source_digests
            if isinstance(item, Mapping)
            and (item.get("node_id") == node or item.get("role") == role)
        ]
        if len(matches) != 1:
            raise SourceBindingError(f"independent {role} source digest is missing")
        return matches[0]

    source = _mapping(source_digests, "independent source digests")
    for wrapper in ("source_digests", "source_digest", "digests", "sources"):
        nested = source.get(wrapper, _MISSING)
        if nested is not _MISSING:
            return _digest_entry(nested, role, node)
    aliases = {
        "sigma": (
            "sigma",
            "sigma_source",
            "stable_sigma",
            "flat_sigma",
            "exact_original_source_flat_sigma",
            str(SIGMA_NODE_ID),
            SIGMA_NODE_ID,
        ),
        "sigma_prime": (
            "sigma_prime",
            "sigma_prime_source",
            "stable_sigma_prime",
            "flat_sigma_prime",
            "derivative",
            "original_flat_sigma_prime",
            "sigma_jets",
            str(SIGMA_PRIME_NODE_ID),
            SIGMA_PRIME_NODE_ID,
        ),
    }
    found = [source[key] for key in aliases[role] if key in source]
    if len(found) == 1:
        return _mapping(found[0], f"independent {role} digest")
    if len(found) > 1:
        first = _mapping(found[0], f"independent {role} digest")
        if any(dict(_mapping(item, f"independent {role} digest")) != dict(first) for item in found[1:]):
            raise SourceBindingError(f"independent {role} digests disagree")
        return first
    if any(key in source for key in ("module", "source_module", "module_sha256", "AST_sha256", "ast_sha256")):
        return source
    raise SourceBindingError(f"independent {role} source digest is missing")


def _verify_digest(expected: Mapping[str, Any], supplied: Mapping[str, Any], label: str) -> dict[str, str]:
    expected_identity = _identity(expected, f"{label}.record")
    supplied_identity = _identity(supplied, f"{label}.caller")
    _canonical_guard(expected_identity, f"{label}.record")
    for key in ("module", "module_sha256", "AST_sha256"):
        if supplied_identity[key] != expected_identity[key]:
            raise SourceBindingError(f"{label}: caller digest mismatch for {key}")
    if "function" in expected_identity and supplied_identity.get("function") != expected_identity["function"]:
        raise SourceBindingError(f"{label}: caller function mismatch")
    # This return is a saved-record identity only.  It is not a verification
    # of the new leaf's arithmetic against the historical module.
    return expected_identity | {"ast_status": AST_ASSERTION_STATUS}


class _LeafCallable:
    def __init__(self, node_id: int, kind: str) -> None:
        self.node_id = node_id
        self.kind = kind
        self.__name__ = kind

    def __call__(self, x: Decimal, precision: int) -> Decimal:
        pair = stable_sigma(x, precision)
        return pair[0] if self.kind == "flat_sigma" else pair[1]

    def __repr__(self) -> str:
        return f"<{self.kind} v4 leaf node {self.node_id}>"


class _EvaluateCallable:
    __name__ = "flat_sigma_point_evaluate_v4"

    def __call__(self, x: Decimal, precision: int) -> dict[str, Any]:
        return stable_sigma_point(x, precision)


class _BoundKernel(dict[str, Any]):
    def __call__(self, x: Decimal, precision: int) -> dict[str, Any]:
        return self["evaluate"](x, precision)


def bind_flat_sigma_point_kernel(
    prior_program: Any = None,
    source_digest: Any = None,
    source_digest_prime: Any = None,
    *,
    source_digests: Any = None,
    sigma_digest: Any = None,
    sigma_prime_digest: Any = None,
    prior_result: Any = None,
) -> dict[str, Any]:
    """Bind exact v4 receipt records to fresh point leaf callables.

    Caller-supplied records and digests are compared strictly.  The result
    intentionally says ``source_identity_verified=False``: a metadata match
    supplies provenance only, and the saved AST assertions are not a fresh
    AST calculation or an arithmetic equivalence certificate.
    """

    if prior_program is None:
        prior_program = prior_result
    elif prior_result is not None:
        raise SourceBindingError("provide one prior compiled receipt")
    if source_digest is not None and source_digests is not None:
        raise SourceBindingError("provide one independent digest mapping")
    if source_digest_prime is not None:
        if source_digest is None or source_digests is not None:
            raise SourceBindingError("source_digest_prime requires source_digest")
        supplied = {"sigma": source_digest, "sigma_prime": source_digest_prime}
    elif sigma_digest is not None or sigma_prime_digest is not None:
        if source_digest is not None or source_digests is not None or sigma_digest is None or sigma_prime_digest is None:
            raise SourceBindingError("sigma_digest and sigma_prime_digest are required together")
        supplied = {"sigma": sigma_digest, "sigma_prime": sigma_prime_digest}
    else:
        supplied = source_digests if source_digests is not None else source_digest
    if supplied is None:
        raise SourceBindingError("independent source digests are required")

    program = _program_candidate(prior_program)
    records = _collect_records(program)
    for node in (SIGMA_NODE_ID, SIGMA_PRIME_NODE_ID, ARGUMENT_NODE_ID):
        if node not in records:
            raise SourceBindingError(f"prior receipt is missing node {node}")
    sigma_instruction, sigma_record = records[SIGMA_NODE_ID]
    prime_instruction, prime_record = records[SIGMA_PRIME_NODE_ID]
    _argument_record = records[ARGUMENT_NODE_ID][1]

    _opaque(sigma_instruction, SIGMA_NODE_ID, "exact_original_source_flat_sigma")
    _opaque(prime_instruction, SIGMA_PRIME_NODE_ID, "analytic_unary")
    if _operation(sigma_record, SIGMA_NODE_ID) != "exact_original_source_flat_sigma":
        raise SourceBindingError("sigma operation is not exact_original_source_flat_sigma")
    if sigma_record.get("definition") != SIGMA_DEFINITION:
        raise SourceBindingError("sigma definition is not the original flat source")
    if sigma_record.get("exact_function_not_its_interval_cap") is not True:
        raise SourceBindingError("sigma record is not marked as the exact function")
    if _operation(prime_record, SIGMA_PRIME_NODE_ID) != "analytic_unary":
        raise SourceBindingError("sigma prime operation is not analytic_unary")
    if _function(prime_record, SIGMA_PRIME_NODE_ID) != "original_flat_sigma_prime":
        raise SourceBindingError("sigma prime function is not original_flat_sigma_prime")
    if prime_record.get("ordinary_derivative") != SIGMA_PRIME_ORDINARY_DERIVATIVE:
        raise SourceBindingError("sigma prime ordinary derivative semantic mismatch")
    if prime_record.get("exact_flat_source_derivative_not_a_bound") is not True:
        raise SourceBindingError("sigma prime is not marked as an exact flat derivative")
    sigma_argument = _node_id(sigma_record.get("argument"), "sigma argument")
    prime_argument = _node_id(prime_record.get("argument"), "sigma prime argument")
    if sigma_argument != ARGUMENT_NODE_ID or prime_argument != ARGUMENT_NODE_ID:
        raise SourceBindingError("sigma leaves do not share the v4 argument node")
    if _operation(_argument_record, ARGUMENT_NODE_ID) != "bound_variable":
        raise SourceBindingError("argument record is not bound_variable")
    if _argument_record.get("name") != ARGUMENT_VARIABLE_NAME:
        raise SourceBindingError("argument variable name is not original_lower_collar_s")

    sigma_source = _source(sigma_record, SIGMA_NODE_ID)
    prime_source = _source(prime_record, SIGMA_PRIME_NODE_ID)
    sigma_identity = _verify_digest(
        sigma_source,
        _digest_entry(supplied, "sigma", SIGMA_NODE_ID),
        "sigma source",
    )
    prime_identity = _verify_digest(
        prime_source,
        _digest_entry(supplied, "sigma_prime", SIGMA_PRIME_NODE_ID),
        "sigma prime source",
    )

    sigma_leaf = _LeafCallable(SIGMA_NODE_ID, "flat_sigma")
    prime_leaf = _LeafCallable(SIGMA_PRIME_NODE_ID, "flat_sigma_prime")
    evaluate_leaf = _EvaluateCallable()
    functions: dict[Any, Any] = {
        SIGMA_NODE_ID: sigma_leaf,
        str(SIGMA_NODE_ID): sigma_leaf,
        "exact_original_source_flat_sigma": sigma_leaf,
        SIGMA_PRIME_NODE_ID: prime_leaf,
        str(SIGMA_PRIME_NODE_ID): prime_leaf,
        "original_flat_sigma_prime": prime_leaf,
    }
    bound_nodes = {
        "sigma": {
            "node_id": SIGMA_NODE_ID,
            "operation": "exact_original_source_flat_sigma",
            "argument": ARGUMENT_NODE_ID,
            "original_record": dict(sigma_record),
            "source_binding": dict(sigma_source),
            "callable": sigma_leaf,
        },
        "sigma_prime": {
            "node_id": SIGMA_PRIME_NODE_ID,
            "operation": "analytic_unary",
            "function": "original_flat_sigma_prime",
            "argument": ARGUMENT_NODE_ID,
            "original_record": dict(prime_record),
            "source_binding": dict(prime_source),
            "derivative_source": dict(prime_source),
            "callable": prime_leaf,
        },
    }
    return _BoundKernel(
        {
            "schema": SCHEMA,
            "component": COMPONENT,
            "status": "bound",
            "prior_schema": PRIOR_SCHEMA,
            "prior_program_schema": PRIOR_SCHEMA,
            "prior_program_status": program.get("status", "compiled"),
            "prior_program_executable": False,
            "prior_program_execution_supported": False,
            "program_executable": False,
            "compiled_program_executable": False,
            "source_identity_verified": False,
            "caller_supplied_source_records_match": True,
            "source_digest_metadata_checked": True,
            "source_identity_status": "caller_record_match_provenance_only",
            "ast_assertion_status": AST_ASSERTION_STATUS,
            "numeric_equivalence_certified": False,
            "certified": False,
            "point_evaluation_only": True,
            "interval": False,
            "full_Navier_Stokes": False,
            "claims": {
                "certified": False,
                "point_evaluation_only": True,
                "interval": False,
                "full_Navier_Stokes": False,
                "prior_program_executable": False,
                "program_executable": False,
                "source_identity_verified": False,
                "caller_supplied_source_records_match": True,
                "ast_recomputed": False,
                "numeric_equivalence_certified": False,
            },
            "bound_argument_node_id": ARGUMENT_NODE_ID,
            "bound_nodes": bound_nodes,
            "bound_nodes_by_id": {
                SIGMA_NODE_ID: bound_nodes["sigma"],
                SIGMA_PRIME_NODE_ID: bound_nodes["sigma_prime"],
            },
            "source_identities": {
                "sigma": sigma_identity,
                "sigma_prime": prime_identity,
            },
            "sigma": sigma_leaf,
            "sigma_prime": prime_leaf,
            "flat_sigma": sigma_leaf,
            "flat_sigma_prime": prime_leaf,
            "evaluate": evaluate_leaf,
            "point_evaluate": evaluate_leaf,
            "callable": evaluate_leaf,
            "functions": functions,
            "leaf_functions": functions,
        }
    )


bind_point_kernel = bind_flat_sigma_point_kernel
bind_prior_program = bind_flat_sigma_point_kernel
bind = bind_flat_sigma_point_kernel





