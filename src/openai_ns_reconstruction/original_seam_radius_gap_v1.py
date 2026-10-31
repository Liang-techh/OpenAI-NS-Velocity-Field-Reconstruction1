"""Conditional symbolic radius-seam gap constructor.

This leaf constructs only a JSON-safe expression DAG.  ``build`` accepts the
twelve source-role ports explicitly required by the seam algebra; each port
is retained as an opaque token descriptor and is never interpreted as a
number, a graph node, or a source value.  No source function, graph, archive,
gzip member, ledger, or error certificate is read here.

The only radius relation used as a condition premise is
``Ro = Rf * exp(D)``.  Here ``D`` may be the full
``log(R_old/Rf)``, including old/current ``y`` and ``mu`` differences; no
``y_old == yF`` shortcut is taken.  The returned gap is the unsimplified signed expression
``delta_K = K_old - K_new`` with independent ``Ho``/``Hf`` and ``qyo``/``qyf``
ports.  Its ordinary-Z1 expression retains every supplied derivative,
including ``D_Z``.  Owner identity, common coordinates and units, derivative
semantics, certified errors, and seam cancellation are deliberately recorded
as unproved obligations rather than inferred by this constructor.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any


SCHEMA = "original_seam_radius_gap_v1"
COMPONENT = SCHEMA

PORT_NAMES = (
    "Rf",
    "Rf_Z",
    "D",
    "D_Z",
    "Ho",
    "Ho_Z",
    "Hf",
    "Hf_Z",
    "qyo",
    "qyo_Z",
    "qyf",
    "qyf_Z",
)


class SeamRadiusGapInputError(ValueError):
    """Raised when the source-port mapping is missing or has extra keys."""


def _expr(operator: str, *arguments: Any, **metadata: Any) -> dict[str, Any]:
    """Create one unsimplified JSON-safe expression node."""

    result: dict[str, Any] = {
        "operator": operator,
        "arguments": list(arguments),
        "source_values_selected": False,
    }
    result.update(metadata)
    return result


def _literal(text: str) -> dict[str, Any]:
    """Create a syntax literal; it is not a selected source value."""

    return {
        "operator": "exact_literal",
        "text": text,
        "source_values_selected": False,
        "node_id_is_not_a_numeric_value": True,
    }


def _token_descriptor(token: Any) -> dict[str, Any]:
    """Describe an input token without copying or interpreting its payload.

    In particular, numeric-looking token payloads (including floats and
    capacity-like objects) are represented only by their Python type.  The
    descriptor therefore cannot be mistaken for a source scalar by a DAG
    consumer.
    """

    if isinstance(token, Mapping):
        return {
            "kind": "opaque_mapping_token",
            "mapping_key_count": len(token),
            "token_payload_read": False,
        }
    return {
        "kind": "opaque_input_token",
        "python_type": type(token).__name__,
        "token_payload_read": False,
    }


def _port_node(name: str, token: Any) -> dict[str, Any]:
    """Wrap one caller token as an opaque source port."""

    return _expr(
        "opaque_input_token",
        port=name,
        token_descriptor=_token_descriptor(token),
        token_is_opaque=True,
        node_id_is_not_a_numeric_value=True,
    )


def _validate_ports(ports: Any) -> dict[str, Any]:
    """Require exactly the twelve named ports and retain no source payload."""

    if not isinstance(ports, Mapping):
        raise SeamRadiusGapInputError("ports must be a mapping of opaque source tokens")
    supplied = set(ports.keys())
    required = set(PORT_NAMES)
    missing = sorted(required - supplied)
    extra = sorted(supplied - required)
    if missing or extra:
        details: list[str] = []
        if missing:
            details.append("missing=" + ",".join(missing))
        if extra:
            details.append("extra=" + ",".join(str(item) for item in extra))
        raise SeamRadiusGapInputError("ports key set must be exact (" + "; ".join(details) + ")")
    return {name: ports[name] for name in PORT_NAMES}


def _premise(equation: str, reason: str) -> dict[str, Any]:
    return {
        "kind": "symbolic_condition_premise",
        "equation": equation,
        "premise_only": True,
        "proved": False,
        "source_values_selected": False,
        "reason": reason,
    }


def build(ports: Mapping[str, Any]) -> dict[str, Any]:
    """Construct the conditional radius-gap expression DAG.

    ``ports`` must contain exactly ``PORT_NAMES``.  Its values are opaque
    caller tokens: this function records only token shape metadata, never the
    token payload.  The result is a plain mapping containing JSON-safe nodes,
    with no numeric source evaluation.
    """

    supplied = _validate_ports(ports)
    p = {name: _port_node(name, supplied[name]) for name in PORT_NAMES}

    # The deliberately unsimplified factorization requested by the seam
    # contract.  exp(D) and exp(2D) occur only in the derivative terms; e and
    # e2 retain the stable expm1-based gap factors in the value expression.
    e = _expr("expm1", p["D"], formula_text="e=expm1(D)")
    e2 = _expr(
        "multiply",
        e,
        _expr("add", _literal("2"), e),
        formula_text="e2=e*(2+e)",
    )
    exp_D = _expr("exp", p["D"], formula_text="exp(D)")
    exp_2D = _expr(
        "exp",
        _expr("multiply", _literal("2"), p["D"]),
        formula_text="exp(2D)",
    )

    G = _expr(
        "add",
        _expr("subtract", p["Ho"], p["Hf"]),
        _expr("multiply", p["Ho"], e2),
        formula_text="G=(Ho-Hf)+Ho*e2",
    )
    T = _expr(
        "add",
        _expr("subtract", p["qyo"], p["qyf"]),
        _expr("multiply", p["qyo"], e),
        formula_text="T=(qyo-qyf)+qyo*e",
    )
    Rf_squared = _expr("power", p["Rf"], _literal("2"), formula_text="Rf^2")
    delta_K = _expr(
        "subtract",
        _expr("multiply", Rf_squared, G),
        _expr("multiply", p["Rf"], T),
        formula_text="delta_K=Rf^2*G-Rf*T",
    )

    G_Z = _expr(
        "add",
        _expr("subtract", p["Ho_Z"], p["Hf_Z"]),
        _expr(
            "multiply",
            p["Ho_Z"],
            e2,
        ),
        _expr(
            "multiply",
            _literal("2"),
            p["Ho"],
            exp_2D,
            p["D_Z"],
        ),
        formula_text="G_Z=(Ho_Z-Hf_Z)+Ho_Z*e2+2*Ho*exp(2D)*D_Z",
    )
    T_Z = _expr(
        "add",
        _expr("subtract", p["qyo_Z"], p["qyf_Z"]),
        _expr("multiply", p["qyo_Z"], e),
        _expr(
            "multiply",
            p["qyo"],
            exp_D,
            p["D_Z"],
        ),
        formula_text="T_Z=(qyo_Z-qyf_Z)+qyo_Z*e+qyo*exp(D)*D_Z",
    )
    delta_K_Z = _expr(
        "add",
        _expr("multiply", _literal("2"), p["Rf"], p["Rf_Z"], G),
        _expr("multiply", Rf_squared, G_Z),
        _expr("negate", _expr("multiply", p["Rf_Z"], T)),
        _expr("negate", _expr("multiply", p["Rf"], T_Z)),
        formula_text="delta_K_Z=2*Rf*Rf_Z*G+Rf^2*G_Z-Rf_Z*T-Rf*T_Z",
    )

    assumptions = {
        "Ro_equals_Rf_exp_D": _premise(
            "Ro=Rf*exp(D)",
            "The old-radius relation is a condition premise only; Ro is not an input or evaluated node. Positive radii, compatible units, and coordinate pairing are unproved.",
        ),
        "D_full_log_radius_ratio": _premise(
            "D=log(R_old/Rf)",
            "D may include old/current y and mu differences; y_old=yF is not assumed.",
        ),
        "positive_radius_and_unit_compatibility": _premise(
            "Ro>0 and Rf>0 with common radius units",
            "The exponential radius relation requires a positive, unit-compatible radius pairing; this leaf does not certify it.",
        ),
        "owner_identity": _premise(
            "owner(source endpoint o)=owner(source endpoint f)",
            "Different-owner endpoint tokens remain independent; owner equality is not proved.",
        ),
        "same_coordinate_semantics": _premise(
            "coordinate(o)=coordinate(f)",
            "Coordinate identity and pullback semantics are outside this symbolic leaf.",
        ),
        "same_units": _premise(
            "units(o)=units(f)",
            "Dimensional compatibility is not inferred from opaque ports.",
        ),
        "ordinary_Z1_derivative_semantics": _premise(
            "ordinary_Z1=the supplied Z derivative of each port",
            "The DAG preserves every *_Z port but does not certify its source meaning.",
        ),
        "certified_errors": _premise(
            "certified_error_transport",
            "No error bound or certification packet is read or produced.",
        ),
        "seam_cancellation": _premise(
            "seam_cancellation",
            "The independent endpoint difference is retained; cancellation is unproved.",
        ),
    }

    result: dict[str, Any] = {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "symbolic_dag_constructed",
        "constructed_means_symbolic_records_only": True,
        "source_values_selected": False,
        "source_values_computed": False,
        "numeric_values_computed": False,
        "source_graph_accessed": False,
        "source_graph_traversed": False,
        "gzip_read": False,
        "ledger_modified": False,
        "ports": p,
        "input_ports": p,
        "intermediates": {
            "e": e,
            "e2": e2,
            "exp_D": exp_D,
            "exp_2D": exp_2D,
            "G": G,
            "T": T,
            "G_Z": G_Z,
            "T_Z": T_Z,
        },
        "delta_K": delta_K,
        "delta_K_Z": delta_K_Z,
        "ordinary_Z1": {
            "delta_K": delta_K,
            "delta_K_Z": delta_K_Z,
            "derivative_order": 1,
            "source_values_selected": False,
        },
        "premises": assumptions,
        "scope": {
            "owner_identity_proved": False,
            "source_identity_proved": False,
            "same_coordinate_semantics_proved": False,
            "coordinate_pairing_proved": False,
            "same_units_proved": False,
            "positive_radius_and_unit_compatibility_proved": False,
            "D_full_log_radius_ratio_semantics_proved": False,
            "ordinary_Z1_semantics_proved": False,
            "certified_errors": False,
            "seam_cancellation_proved": False,
            "ledger_changed": False,
        },
    }
    return result


# Descriptive aliases keep the one-constructor API discoverable without
# introducing a second implementation or a second source of formula truth.
build_symbolic_dag = build
construct = build


__all__ = [
    "COMPONENT",
    "PORT_NAMES",
    "SCHEMA",
    "SeamRadiusGapInputError",
    "build",
    "build_symbolic_dag",
    "construct",
]
