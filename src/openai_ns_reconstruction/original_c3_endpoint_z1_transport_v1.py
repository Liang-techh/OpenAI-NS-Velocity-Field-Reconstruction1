"""Conditional Z1 transport for a historical C3 endpoint composite.

The old C1 producer made ``radius`` and ``delta`` with ``C1Algebra.fixed``.
Consequently an old saved first row can be a field-row derivative while still
omitting the ordinary-Z derivatives of the radius and of ``delta``.  This
module emits the small, source-only correction needed to carry that row into
the full moving-radius/moving-``delta`` expression.

``build`` receives exactly twelve opaque source references.  A reference is
copied as JSON metadata and is never used as a scalar, graph node, function,
or source value.  The result is an unevaluated JSON-safe expression DAG; it
does not read an archive, graph, gzip member, or historical producer.

The algebra is, with ``L = 1-delta*Z**2``,

``K = A0/(2*L) + B0/4 - Z**2*C0/L**2 - D0``

and the old fixed-geometry candidate

``K1_fixed = A1/(2*L) - A0*Lfixed_Z/(2*L**2) + B1/4``
``             - 2*Z*C0/L**2 - Z**2*C1/L**2``
``             + 2*Z**2*C0*Lfixed_Z/L**3 - D1``

where ``Lfixed_Z = -2*delta*Z``.  The full row is the unsimplified sum
``K1_full = K1_fixed + rho_Z*(...) + delta_Z*(...)`` specified by the
historical transport contract.  No cancellation or bookkeeping operation is
performed here.
"""

from __future__ import annotations

from collections.abc import Mapping
from math import isfinite
from typing import Any


SCHEMA = "original_c3_endpoint_z1_transport_v1"
COMPONENT = SCHEMA

# The order is stable so serialized input maps and diagnostics are deterministic.
REF_NAMES = (
    "A0",
    "A1",
    "B0",
    "B1",
    "C0",
    "C1",
    "D0",
    "D1",
    "Z",
    "delta",
    "rho_Z",
    "delta_Z",
)

# ``PORT_NAMES`` is a descriptive compatibility alias used by other
# source-only leaves for an exact input-port contract.
PORT_NAMES = REF_NAMES


class OriginalC3EndpointZ1TransportInputError(ValueError):
    """Raised when the opaque-reference mapping is not exactly admitted."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise OriginalC3EndpointZ1TransportInputError(message)


def _copy_json(value: Any, path: str) -> Any:
    """Copy a JSON value while refusing non-JSON source payloads.

    The copy is metadata only.  In particular, numeric-looking members inside
    an opaque reference remain inside the opaque object and are never lifted
    into an arithmetic node.
    """

    if value is None or type(value) is bool or type(value) is int or type(value) is str:
        return value
    if type(value) is float:
        _require(isfinite(value), f"{path}: non-finite JSON number is not admitted")
        return value
    if isinstance(value, Mapping):
        copied: dict[str, Any] = {}
        for key, item in value.items():
            _require(type(key) is str, f"{path}: JSON object keys must be strings")
            copied[key] = _copy_json(item, f"{path}.{key}")
        return copied
    if isinstance(value, (list, tuple)):
        return [_copy_json(item, f"{path}[{index}]") for index, item in enumerate(value)]
    raise OriginalC3EndpointZ1TransportInputError(
        f"{path}: opaque reference must be JSON-safe"
    )


def _validate_refs(refs: Any) -> dict[str, Any]:
    """Require the exact twelve names and copy each opaque JSON reference."""

    _require(isinstance(refs, Mapping), "refs must be a mapping of opaque JSON references")
    supplied = set(refs.keys())
    required = set(REF_NAMES)
    missing = sorted(required - supplied)
    extra = sorted(supplied - required, key=str)
    _require(not missing and not extra,
             "refs key set must be exact (missing=%s; extra=%s)"
             % (",".join(missing), ",".join(str(item) for item in extra)))
    return {name: _copy_json(refs[name], f"refs.{name}") for name in REF_NAMES}


def _expr(operator: str, *arguments: Any, formula_text: str | None = None) -> dict[str, Any]:
    """Create one JSON-safe symbolic node with no source-value selection."""

    result: dict[str, Any] = {
        "operator": operator,
        "arguments": list(arguments),
        "source_values_selected": False,
    }
    if formula_text is not None:
        result["formula_text"] = formula_text
    return result


def _literal(text: str) -> dict[str, Any]:
    return {
        "operator": "exact_literal",
        "text": text,
        "source_values_selected": False,
        "node_id_is_not_a_numeric_value": True,
    }


def _opaque_source_ref(name: str, ref: Any) -> dict[str, Any]:
    """Wrap a copied ref as an opaque operand, never as a numeric literal."""

    return {
        "operator": "opaque_source_ref",
        "arguments": [_copy_json(ref, f"refs.{name}")],
        "source_role": name,
        "token_is_opaque": True,
        "node_id_is_not_a_numeric_value": True,
        "source_values_selected": False,
    }


def _sum(*terms: Any, formula_text: str | None = None) -> dict[str, Any]:
    return _expr("sum", *terms, formula_text=formula_text)


def _product(*factors: Any, formula_text: str | None = None) -> dict[str, Any]:
    return _expr("product", *factors, formula_text=formula_text)


def _quotient(numerator: Any, denominator: Any, formula_text: str | None = None) -> dict[str, Any]:
    return _expr("quotient", numerator, denominator, formula_text=formula_text)


def _negative(argument: Any, formula_text: str | None = None) -> dict[str, Any]:
    return _expr("negative", argument, formula_text=formula_text)


def _power(base: Any, exponent: str, formula_text: str | None = None) -> dict[str, Any]:
    return _expr("power", base, _literal(exponent), formula_text=formula_text)


def _scope() -> dict[str, Any]:
    """Return the deliberately conditional scope ledger for one build."""

    return {
        "old_Z1_is_selected_endpoint_ordinary_derivative": False,
        "old_Z1_capture_safe_same_variable_and_physical_scale": False,
        "positive_radius_and_nonzero_L": False,
        "early_C3_endpoint_source_identity": False,
        "rho_Z_and_delta_Z_independently_supplied": True,
        "rho_Z_and_delta_Z_independence_proved": False,
        "historical_source_nonzero_geometry_derivative_observed": False,
        "synthetic_fixture_nonzero_geometry_only": True,
        "source_identity_proved": False,
        "source_values_selected": False,
        "source_values_computed": False,
        "numeric_evaluation_performed": False,
        "cancellation_performed": False,
        "bookkeeping_or_ledger_changed": False,
        "error_certified": False,
        "condition": (
            "Interpret K1_full only when the old field rows are the selected "
            "ordinary-Z1 rows in one capture-safe variable and physical scale, "
            "the radius is positive, L is nonzero, the early C3 endpoint source "
            "identity remains conditional, and rho_Z and delta_Z are independent "
            "caller inputs."
        ),
    }


def build(refs: Mapping[str, Any]) -> dict[str, Any]:
    """Build the source-opaque historical C3 endpoint Z1 transport DAG.

    ``refs`` must contain exactly :data:`REF_NAMES`.  Each value can be any
    JSON object (or other JSON value) used as an opaque source reference; it is
    copied into the DAG and never interpreted as a scalar.
    """

    copied = _validate_refs(refs)
    p = {name: _opaque_source_ref(name, copied[name]) for name in REF_NAMES}

    Z = p["Z"]
    delta = p["delta"]
    A0, A1 = p["A0"], p["A1"]
    B0, B1 = p["B0"], p["B1"]
    C0, C1 = p["C0"], p["C1"]
    D0, D1 = p["D0"], p["D1"]

    Z2 = _power(Z, "2", formula_text="Z^2")
    Z4 = _power(Z, "4", formula_text="Z^4")
    L = _sum(
        _literal("1"),
        _negative(_product(delta, Z2), formula_text="-delta*Z^2"),
        formula_text="L=1-delta*Z^2",
    )
    L2 = _power(L, "2", formula_text="L^2")
    L3 = _power(L, "3", formula_text="L^3")
    Lfixed_Z = _negative(
        _product(_literal("2"), delta, Z),
        formula_text="Lfixed_Z=-2*delta*Z",
    )

    K = _sum(
        _quotient(A0, _product(_literal("2"), L), formula_text="A0/(2*L)"),
        _quotient(B0, _literal("4"), formula_text="B0/4"),
        _negative(_quotient(_product(Z2, C0), L2), formula_text="-Z^2*C0/L^2"),
        _negative(D0, formula_text="-D0"),
        formula_text="K=A0/(2*L)+B0/4-Z^2*C0/L^2-D0",
    )

    K1_fixed = _sum(
        _quotient(A1, _product(_literal("2"), L), formula_text="A1/(2*L)"),
        _negative(
            _quotient(_product(A0, Lfixed_Z), _product(_literal("2"), L2)),
            formula_text="-A0*Lfixed_Z/(2*L^2)",
        ),
        _quotient(B1, _literal("4"), formula_text="B1/4"),
        _negative(
            _quotient(_product(_literal("2"), Z, C0), L2),
            formula_text="-2*Z*C0/L^2",
        ),
        _negative(
            _quotient(_product(Z2, C1), L2),
            formula_text="-Z^2*C1/L^2",
        ),
        _quotient(
            _product(_literal("2"), Z2, C0, Lfixed_Z),
            L3,
            formula_text="2*Z^2*C0*Lfixed_Z/L^3",
        ),
        _negative(D1, formula_text="-D1"),
        formula_text=(
            "K1_fixed=A1/(2*L)-A0*Lfixed_Z/(2*L^2)+B1/4"
            "-2*Z*C0/L^2-Z^2*C1/L^2+2*Z^2*C0*Lfixed_Z/L^3-D1"
        ),
    )

    rho_inner = _sum(
        _quotient(A0, L, formula_text="A0/L"),
        _quotient(B0, _literal("2"), formula_text="B0/2"),
        _negative(
            _quotient(_product(_literal("2"), Z2, C0), L2),
            formula_text="-2*Z^2*C0/L^2",
        ),
        _negative(D0, formula_text="-D0"),
        formula_text="A0/L+B0/2-2*Z^2*C0/L^2-D0",
    )
    rho_correction = _product(
        p["rho_Z"],
        rho_inner,
        formula_text="rho_Z*(A0/L+B0/2-2*Z^2*C0/L^2-D0)",
    )

    delta_inner = _sum(
        _quotient(
            _product(Z2, A0),
            _product(_literal("2"), L2),
            formula_text="Z^2*A0/(2*L^2)",
        ),
        _negative(
            _quotient(_product(_literal("2"), Z4, C0), L3),
            formula_text="-2*Z^4*C0/L^3",
        ),
        formula_text="Z^2*A0/(2*L^2)-2*Z^4*C0/L^3",
    )
    delta_correction = _product(
        p["delta_Z"],
        delta_inner,
        formula_text="delta_Z*(Z^2*A0/(2*L^2)-2*Z^4*C0/L^3)",
    )
    K1_full = _sum(
        K1_fixed,
        rho_correction,
        delta_correction,
        formula_text="K1_full=K1_fixed+rho_correction+delta_correction",
    )

    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "symbolic_dag_constructed",
        "source_values_selected": False,
        "source_values_computed": False,
        "numeric_values_computed": False,
        "source_graph_accessed": False,
        "source_graph_traversed": False,
        "gzip_read": False,
        "inputs": p,
        "intermediates": {
            "L": L,
            "Lfixed_Z": Lfixed_Z,
            "rho_correction": rho_correction,
            "delta_correction": delta_correction,
        },
        "K": K,
        "K1_fixed": K1_fixed,
        "K1_full": K1_full,
        "scope": _scope(),
    }


build_symbolic_dag = build
construct = build


__all__ = [
    "COMPONENT",
    "OriginalC3EndpointZ1TransportInputError",
    "PORT_NAMES",
    "REF_NAMES",
    "SCHEMA",
    "build",
    "build_symbolic_dag",
    "construct",
]
