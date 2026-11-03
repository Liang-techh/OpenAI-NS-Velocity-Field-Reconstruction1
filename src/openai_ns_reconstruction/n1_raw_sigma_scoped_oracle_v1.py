"""Scoped interval oracle for the two saved raw sigma leaves.

This is intentionally a leaf oracle, not a raw-g5 evaluator.  The caller
must provide a runner-pinned graph sidecar and a saved interval result whose
authentication is only the runner caller's assertion.  No historical producer, AST walker,
numeric sigma evaluator, or normalized-cutoff implementation is imported.

The oracle is designed for ``n1_scoped_graph_evaluator_v1``.  Its callable
boundary is ``(node_id, record, category, env)``.  Only nodes ``35949`` and
``35953`` are accepted, and both must use bound-variable node ``35848`` with
the original ``original_lower_collar_s`` coordinate.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType
from typing import Any


SCHEMA = "n1_raw_sigma_scoped_oracle_v1"
RESULT_SCHEMA = SCHEMA + "_result"
COMPONENT = "raw_sigma_scoped_interval_oracle"
AUTHENTICATION_BASIS = "runner_caller_assertion"

GRAPH_SHA256 = "1b44be8c6cbfffb692c7a0430ecbc8a590f46001a233c99c59ca8c5e2c44a6fe"
EXPECTED_GRAPH_SHA256 = GRAPH_SHA256
ARGUMENT_NODE_ID = 35848
SIGMA_NODE_ID = 35949
SIGMA_PRIME_NODE_ID = 35953
AUTHORIZED_NODE_IDS = (SIGMA_NODE_ID, SIGMA_PRIME_NODE_ID)
RAW_S_COORDINATE = "original_lower_collar_s"
RAW_S_FORM = "raw_s"
SAFE_ARGUMENT_INTERVAL = (
    Fraction(0),
    Fraction(1, 2**128),
)
SOURCE_DEFINITION = (
    "0 for x<=0; 1 for x>=1; "
    "exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside"
)
SIGMA_SOURCE = MappingProxyType(
    {
        "AST_sha256": "4eca34dd2d6bd23b5bff76496682bece54924d3145ac99e95ed9066d9572ab3c",
        "function": "stable_sigma",
        "module": "lei_ren_part1_paper_compliant_outer_initial.py",
        "module_sha256": "a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657",
    }
)
SIGMA_PRIME_SOURCE = MappingProxyType(
    {
        "AST_sha256": "a63c47d05049c86007b530f3095c515d724bef53092c5f7867d53e5bbaa06236",
        "function": "sigma_jets",
        "module": "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py",
        "module_sha256": "32f51123a47c2bafbb6c2f7056302bab996a019883615333875207d72a9bb7d7",
    }
)


class ScopedSigmaOracleError(ValueError):
    """Base class for fail-closed saved-source and scope errors."""


class AuthenticationError(ScopedSigmaOracleError):
    """The graph or interval packet was not externally authenticated."""


class GraphRecordError(ScopedSigmaOracleError):
    """A node record is not the authenticated raw record."""


class BinderError(ScopedSigmaOracleError):
    """The required lexical argument binder is missing or malformed."""


class IntervalInputError(ScopedSigmaOracleError):
    """The supplied closed Fraction interval is malformed or out of scope."""


_MISSING = object()


def _sha256(value: Any, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64:
        raise AuthenticationError(f"{label} must be a SHA-256 hex string")
    if any(ch not in "0123456789abcdefABCDEF" for ch in value):
        raise AuthenticationError(f"{label} must be a SHA-256 hex string")
    return value.lower()


def _canonical(value: Any, label: str = "value") -> Any:
    """Canonicalize JSON-shaped records without interpreting source payload."""

    if value is None:
        return ("none", None)
    if type(value) is bool:
        return ("bool", value)
    if type(value) is int:
        return ("int", value)
    if type(value) is str:
        return ("str", value)
    if isinstance(value, Mapping):
        items: list[tuple[str, Any]] = []
        for key, item in value.items():
            if not isinstance(key, str):
                raise GraphRecordError(f"{label} has a non-string key")
            items.append((key, _canonical(item, f"{label}.{key}")))
        return ("mapping", tuple(sorted(items)))
    if isinstance(value, list):
        return ("list", tuple(_canonical(item, f"{label}[{i}]") for i, item in enumerate(value)))
    if isinstance(value, tuple):
        return ("tuple", tuple(_canonical(item, f"{label}[{i}]") for i, item in enumerate(value)))
    raise GraphRecordError(f"{label} is not JSON-shaped")


def _alias(mapping: Mapping[str, Any], names: tuple[str, ...], label: str) -> Any:
    found = [(name, mapping[name]) for name in names if name in mapping]
    if not found:
        return _MISSING
    first = found[0][1]
    first_key = found[0][0]
    for name, value in found[1:]:
        try:
            equal = value == first
        except Exception as exc:  # pragma: no cover - hostile caller boundary
            raise AuthenticationError(f"{label} aliases are not comparable") from exc
        if equal is not True:
            raise AuthenticationError(f"conflicting {label} aliases: {first_key}/{name}")
    return first


def _nested_alias(
    outer: Mapping[str, Any],
    inner: Mapping[str, Any],
    names: tuple[str, ...],
    label: str,
) -> Any:
    """Read an aliased field at both result levels and reject disagreement."""

    outer_value = _alias(outer, names, label)
    inner_value = _alias(inner, names, label)
    if outer_value is _MISSING:
        return inner_value
    if inner_value is _MISSING:
        return outer_value
    try:
        equal = outer_value == inner_value
    except Exception as exc:  # pragma: no cover - hostile caller boundary
        raise AuthenticationError(f"{label} levels are not comparable") from exc
    if equal is not True:
        raise AuthenticationError(f"conflicting {label} outer/inner values")
    return outer_value


def _metadata(mapping: Mapping[str, Any], label: str) -> None:
    for key in ("normalized", "is_normalized", "normalized_cutoff"):
        if key in mapping and mapping[key] is not False:
            raise IntervalInputError(f"{label} is not raw-s: {key}")
    for key in ("coordinate", "coordinate_unit", "unit", "variable", "source_form"):
        if key not in mapping:
            continue
        value = mapping[key]
        if key == "source_form":
            if value != RAW_S_FORM:
                raise IntervalInputError(f"{label} has non-raw source form")
        elif value != RAW_S_COORDINATE and value != RAW_S_FORM:
            raise IntervalInputError(f"{label} has non-raw coordinate/unit")


def _fraction_interval(value: Any, label: str) -> "FractionInterval":
    if isinstance(value, FractionInterval):
        return value
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) != 2:
            raise IntervalInputError(f"{label} sequence must contain exactly two endpoints")
        lower, upper = value
    elif isinstance(value, Mapping):
        _metadata(value, label)
        nested = _alias(value, ("interval", "enclosure"), label + ".nested")
        if nested is not _MISSING:
            return _fraction_interval(nested, label)
        lower = value.get("lower", _MISSING)
        upper = value.get("upper", _MISSING)
    else:
        lower = getattr(value, "lower", _MISSING)
        upper = getattr(value, "upper", _MISSING)
        for key in ("normalized", "is_normalized"):
            marker = getattr(value, key, _MISSING)
            if marker is not _MISSING and marker is not False:
                raise IntervalInputError(f"{label} is not raw-s: {key}")
        for key in ("coordinate", "coordinate_unit", "unit", "variable", "source_form"):
            marker = getattr(value, key, _MISSING)
            if marker is _MISSING:
                continue
            if key == "source_form" and marker != RAW_S_FORM:
                raise IntervalInputError(f"{label} has non-raw source form")
            if key != "source_form" and marker not in (RAW_S_COORDINATE, RAW_S_FORM):
                raise IntervalInputError(f"{label} has non-raw coordinate/unit")
    if not isinstance(lower, Fraction) or not isinstance(upper, Fraction):
        raise IntervalInputError(f"{label} endpoints must be Fraction instances")
    return FractionInterval(lower, upper)


def _interval_alias(mapping: Mapping[str, Any], names: tuple[str, ...], label: str) -> Any:
    """Normalize each interval alias before checking that aliases agree."""

    found = [
        (name, _fraction_interval(mapping[name], f"{label}.{name}"))
        for name in names
        if name in mapping
    ]
    if not found:
        return _MISSING
    first_name, first = found[0]
    for name, value in found[1:]:
        if value != first:
            raise AuthenticationError(f"conflicting {label} aliases: {first_name}/{name}")
    return first


@dataclass(frozen=True)
class FractionInterval:
    """A closed interval whose endpoints are exact ``Fraction`` values."""

    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        if not isinstance(self.lower, Fraction) or not isinstance(self.upper, Fraction):
            raise TypeError("FractionInterval endpoints must be Fraction instances")
        if self.lower > self.upper:
            raise ValueError("FractionInterval endpoints must be ordered")

    def contains(self, value: Fraction) -> bool:
        return isinstance(value, Fraction) and self.lower <= value <= self.upper

    def contains_interval(self, other: "FractionInterval") -> bool:
        return self.lower <= other.lower and other.upper <= self.upper


@dataclass(frozen=True)
class AuthenticatedSigmaBounds:
    """Saved bounds accepted only under a runner-caller assertion.

    The graph hash is the runner's pinned target graph; it is not a claim that
    an older interval receipt independently proves association with that graph.
    """

    graph_sha256: str
    argument_interval: FractionInterval
    sigma_interval: FractionInterval
    sigma_s_interval: FractionInterval
    authenticated: bool = False
    source_result_sha256: str | None = None
    receipt_sha256: str | None = None
    authentication_basis: str = AUTHENTICATION_BASIS
    coordinate: str = RAW_S_COORDINATE
    normalized: bool = False
    source_identity_verified: bool = False
    owner_binding: bool = False

    def __post_init__(self) -> None:
        graph_hash = _sha256(self.graph_sha256, "graph_sha256")
        if graph_hash != EXPECTED_GRAPH_SHA256:
            raise AuthenticationError("unexpected raw graph SHA-256")
        object.__setattr__(self, "graph_sha256", graph_hash)
        if self.authenticated is not True:
            raise AuthenticationError("saved bounds require a runner caller assertion")
        if self.authentication_basis != AUTHENTICATION_BASIS:
            raise AuthenticationError("saved bounds authentication is only a runner-caller assertion")
        for label, interval in (
            ("argument_interval", self.argument_interval),
            ("sigma_interval", self.sigma_interval),
            ("sigma_s_interval", self.sigma_s_interval),
        ):
            if not isinstance(interval, FractionInterval):
                raise IntervalInputError(f"{label} must be a closed FractionInterval")
        if self.coordinate != RAW_S_COORDINATE or self.normalized is not False:
            raise AuthenticationError("normalized cutoff bounds are not accepted")
        if self.source_identity_verified is not False or self.owner_binding is not False:
            raise AuthenticationError("source-owner claims must remain false")
        safe = FractionInterval(*SAFE_ARGUMENT_INTERVAL)
        if self.argument_interval.lower != safe.lower or self.argument_interval.upper > safe.upper:
            raise IntervalInputError("saved argument bounds exceed [0, 2^-128]")
        if self.sigma_interval.lower < 0 or self.sigma_interval.upper > 1:
            raise IntervalInputError("saved sigma bounds must lie in [0, 1]")
        if self.sigma_s_interval.lower < 0:
            raise IntervalInputError("saved sigma_s lower bound must be nonnegative")
        if self.source_result_sha256 is not None:
            object.__setattr__(self, "source_result_sha256", _sha256(self.source_result_sha256, "source_result_sha256"))
        if self.receipt_sha256 is not None:
            object.__setattr__(self, "receipt_sha256", _sha256(self.receipt_sha256, "receipt_sha256"))


@dataclass(frozen=True)
class ScopedSigmaInterval:
    """Frozen interval scalar returned to a scoped evaluator backend."""

    lower: Fraction
    upper: Fraction
    node_id: int
    argument_node_id: int
    input_interval: FractionInterval
    quantity: str
    graph_sha256: str
    source_result_sha256: str | None = None
    receipt_sha256: str | None = None
    certified: bool = False
    conditional: bool = True
    source_identity_verified: bool = False
    owner_binding: bool = False
    normalized: bool = False
    full_g5: bool = False
    physical_m15: bool = False
    coordinate: str = RAW_S_COORDINATE
    source_form: str = RAW_S_FORM

    def __post_init__(self) -> None:
        if not isinstance(self.lower, Fraction) or not isinstance(self.upper, Fraction):
            raise TypeError("ScopedSigmaInterval endpoints must be Fraction instances")
        if self.lower > self.upper:
            raise ValueError("ScopedSigmaInterval endpoints must be ordered")
        if self.node_id not in AUTHORIZED_NODE_IDS or self.argument_node_id != ARGUMENT_NODE_ID:
            raise GraphRecordError("interval result has an unauthorized node identity")
        if self.quantity not in ("sigma", "sigma_s"):
            raise GraphRecordError("interval result has an unauthorized quantity")
        if self.certified is not False or self.conditional is not True:
            raise GraphRecordError("interval certification must remain a conditional caller result")
        if self.source_identity_verified or self.owner_binding or self.normalized:
            raise GraphRecordError("interval result contains a forbidden source claim")
        if self.coordinate != RAW_S_COORDINATE or self.source_form != RAW_S_FORM:
            raise GraphRecordError("interval result is not in raw s form")

    @property
    def interval(self) -> FractionInterval:
        return FractionInterval(self.lower, self.upper)

    @property
    def value(self) -> FractionInterval:
        return self.interval

    @property
    def sigma(self) -> FractionInterval | None:
        return self.interval if self.quantity == "sigma" else None

    @property
    def d_sigma_ds(self) -> FractionInterval | None:
        return self.interval if self.quantity == "sigma_s" else None

    @property
    def derivative(self) -> FractionInterval:
        return self.interval

    def is_finite(self) -> bool:
        return True


def _expected_rows() -> dict[int, Any]:
    return {
        ARGUMENT_NODE_ID: {
            "name": RAW_S_COORDINATE,
            "operation": "bound_variable",
        },
        SIGMA_NODE_ID: {
            "argument": ARGUMENT_NODE_ID,
            "definition": SOURCE_DEFINITION,
            "exact_function_not_its_interval_cap": True,
            "operation": "exact_original_source_flat_sigma",
            "source_binding": dict(SIGMA_SOURCE),
        },
        SIGMA_PRIME_NODE_ID: {
            "argument": ARGUMENT_NODE_ID,
            "derivative_source": dict(SIGMA_PRIME_SOURCE),
            "exact_flat_source_derivative_not_a_bound": True,
            "name": "original_flat_sigma_prime",
            "operation": "analytic_unary",
            "ordinary_derivative": "same original sigma Taylor coefficient[1]",
        },
    }


EXPECTED_GRAPH_ROWS = MappingProxyType(_expected_rows())


def authenticated_bounds_from_result(
    result: Mapping[str, Any],
    *,
    authenticated: bool = False,
    graph_sha256: str = EXPECTED_GRAPH_SHA256,
) -> AuthenticatedSigmaBounds:
    """Adapt an already decoded caller-asserted runner result only."""

    if not isinstance(result, Mapping):
        raise AuthenticationError("saved interval result must be a mapping")
    section = result.get("direct_lower_collar", result)
    if not isinstance(section, Mapping):
        raise AuthenticationError("saved direct sigma result must be a mapping")
    marker = _nested_alias(
        result,
        section,
        ("authenticated", "runner_authenticated", "externally_authenticated"),
        "authentication",
    )
    if marker is _MISSING:
        marker = authenticated
    if marker is not True:
        raise AuthenticationError("saved interval result lacks a runner caller assertion")
    result_hash = _nested_alias(result, section, ("graph_sha256", "graph_hash"), "graph hash")
    if result_hash is _MISSING:
        result_hash = graph_sha256
    if _sha256(result_hash, "graph_sha256") != _sha256(graph_sha256, "graph_sha256"):
        raise AuthenticationError("saved result graph hash conflicts with runner hash")
    argument = _interval_alias(
        section,
        ("argument_interval", "evaluated_argument_interval", "interval_leaf_input"),
        "argument interval",
    )
    sigma = _interval_alias(section, ("sigma", "sigma_interval"), "sigma interval")
    sigma_s = _interval_alias(section, ("d_sigma_ds", "sigma_s_interval", "derivative"), "sigma_s interval")
    if argument is _MISSING or sigma is _MISSING or sigma_s is _MISSING:
        raise AuthenticationError("saved result is missing exact Fraction bounds")
    coordinate = _nested_alias(result, section, ("coordinate", "coordinate_unit", "unit"), "coordinate")
    if coordinate is _MISSING:
        coordinate = RAW_S_COORDINATE
    normalized = _nested_alias(
        result,
        section,
        ("normalized", "normalized_inner_cutoff_claimed"),
        "normalized marker",
    )
    if normalized is _MISSING:
        normalized = False
    source_identity = _nested_alias(result, section, ("source_identity_verified",), "source identity")
    owner = _nested_alias(result, section, ("owner_binding", "live_owner_binding"), "owner binding")
    authentication_basis = _nested_alias(
        result,
        section,
        ("authentication_basis", "auth_basis"),
        "authentication basis",
    )
    if authentication_basis is _MISSING:
        authentication_basis = AUTHENTICATION_BASIS
    source_hash = _alias(result, ("source_result_sha256", "result_sha256"), "source result hash")
    receipt_hash = _alias(result, ("receipt_sha256",), "receipt hash")
    return AuthenticatedSigmaBounds(
        graph_sha256=_sha256(graph_sha256, "graph_sha256"),
        argument_interval=_fraction_interval(argument, "saved argument interval"),
        sigma_interval=_fraction_interval(sigma, "saved sigma interval"),
        sigma_s_interval=_fraction_interval(sigma_s, "saved sigma_s interval"),
        authenticated=True,
        source_result_sha256=None if source_hash is _MISSING else source_hash,
        receipt_sha256=None if receipt_hash is _MISSING else receipt_hash,
        authentication_basis=authentication_basis,
        coordinate=coordinate,
        normalized=False if normalized is _MISSING else normalized,
        source_identity_verified=False if source_identity is _MISSING else source_identity,
        owner_binding=False if owner is _MISSING else owner,
    )


def _sidecar_rows(sidecar: Any) -> tuple[str, dict[int, Any]]:
    if not isinstance(sidecar, Mapping):
        raise AuthenticationError("runner-pinned graph sidecar must be a mapping")
    marker = _alias(sidecar, ("authenticated", "runner_authenticated"), "graph authentication")
    if marker is not True:
        raise AuthenticationError("graph rows require a runner caller assertion")
    graph_hash = _alias(sidecar, ("graph_sha256", "graph_hash"), "graph hash")
    if graph_hash is _MISSING:
        raise AuthenticationError("graph sidecar graph hash is required")
    graph_hash = _sha256(graph_hash, "graph_sha256")
    if graph_hash != EXPECTED_GRAPH_SHA256:
        raise AuthenticationError("unexpected raw graph SHA-256")
    rows = _alias(sidecar, ("rows", "exact_graph_nodes", "graph_nodes", "nodes"), "graph rows")
    if rows is _MISSING:
        raise AuthenticationError("graph sidecar rows are required")
    selected: dict[int, Any] = {}
    if isinstance(rows, Mapping):
        for node_id in AUTHORIZED_NODE_IDS + (ARGUMENT_NODE_ID,):
            if node_id not in rows:
                raise GraphRecordError(f"graph sidecar is missing node {node_id}")
            selected[node_id] = rows[node_id]
    elif isinstance(rows, Sequence) and not isinstance(rows, (str, bytes, bytearray)):
        for node_id in AUTHORIZED_NODE_IDS + (ARGUMENT_NODE_ID,):
            if node_id >= len(rows):
                raise GraphRecordError(f"graph sidecar is missing node {node_id}")
            selected[node_id] = rows[node_id]
    else:
        raise AuthenticationError("graph sidecar rows must be a mapping or indexed sequence")
    expected = _expected_rows()
    for node_id, record in selected.items():
        if _canonical(record, f"graph row {node_id}") != _canonical(expected[node_id], f"expected graph row {node_id}"):
            raise GraphRecordError(f"graph row {node_id} does not match authenticated raw graph")
    return graph_hash, expected


class RawSigmaScopedOracle:
    """Callable category oracle for raw sigma and raw sigma-prime only."""

    def __init__(self, bounds: AuthenticatedSigmaBounds, graph_hash: str, rows: Mapping[int, Any]) -> None:
        if bounds.graph_sha256 != graph_hash:
            raise AuthenticationError("saved bounds and graph sidecar hashes differ")
        self._bounds = bounds
        self._graph_hash = graph_hash
        self._rows = {node_id: _canonical(record, f"graph row {node_id}") for node_id, record in rows.items()}

    @property
    def graph_sha256(self) -> str:
        return self._graph_hash

    @property
    def authorized_node_ids(self) -> tuple[int, int]:
        return AUTHORIZED_NODE_IDS

    def __call__(self, node_id: int, record: Any, category: str, env: Any) -> ScopedSigmaInterval:
        if type(node_id) is not int or node_id not in AUTHORIZED_NODE_IDS:
            raise GraphRecordError(f"node {node_id!r} is not authorized")
        if category != "unknown":
            raise GraphRecordError(f"node {node_id} arrived through category {category!r}")
        if _canonical(record, f"node {node_id} record") != self._rows[node_id]:
            raise GraphRecordError(f"node {node_id} record does not match authenticated graph row")
        lookup = getattr(env, "lookup", None)
        if not callable(lookup):
            raise BinderError("scoped oracle requires env.lookup(node_id)")
        try:
            bound = lookup(ARGUMENT_NODE_ID)
        except Exception as exc:
            raise BinderError("argument binder lookup failed") from exc
        input_interval = _fraction_interval(bound, "argument binder interval")
        safe = FractionInterval(*SAFE_ARGUMENT_INTERVAL)
        if not safe.contains_interval(input_interval):
            raise IntervalInputError("argument binder interval exceeds [0, 2^-128]")
        if not self._bounds.argument_interval.contains_interval(input_interval):
            raise IntervalInputError("argument binder interval exceeds saved interval input")
        if node_id == SIGMA_NODE_ID:
            interval = self._bounds.sigma_interval
            quantity = "sigma"
        else:
            interval = self._bounds.sigma_s_interval
            quantity = "sigma_s"
        return ScopedSigmaInterval(
            lower=interval.lower,
            upper=interval.upper,
            node_id=node_id,
            argument_node_id=ARGUMENT_NODE_ID,
            input_interval=input_interval,
            quantity=quantity,
            graph_sha256=self._graph_hash,
            source_result_sha256=self._bounds.source_result_sha256,
            receipt_sha256=self._bounds.receipt_sha256,
        )


def build_oracle(
    saved_result: AuthenticatedSigmaBounds | Mapping[str, Any],
    graph_sidecar: Mapping[str, Any],
    *,
    saved_result_authenticated: bool = False,
) -> RawSigmaScopedOracle:
    """Build the two-leaf oracle from runner-pinned, caller-asserted sidecars."""

    graph_hash, rows = _sidecar_rows(graph_sidecar)
    bounds = (
        saved_result
        if isinstance(saved_result, AuthenticatedSigmaBounds)
        else authenticated_bounds_from_result(
            saved_result,
            authenticated=saved_result_authenticated,
            graph_sha256=graph_hash,
        )
    )
    if bounds.graph_sha256 != graph_hash:
        raise AuthenticationError("saved bounds and graph sidecar hashes differ")
    return RawSigmaScopedOracle(bounds, graph_hash, rows)


make_oracle = build_oracle
scoped_sigma_oracle = RawSigmaScopedOracle


__all__ = [
    "SCHEMA",
    "RESULT_SCHEMA",
    "COMPONENT",
    "AUTHENTICATION_BASIS",
    "GRAPH_SHA256",
    "EXPECTED_GRAPH_SHA256",
    "ARGUMENT_NODE_ID",
    "SIGMA_NODE_ID",
    "SIGMA_PRIME_NODE_ID",
    "AUTHORIZED_NODE_IDS",
    "RAW_S_COORDINATE",
    "RAW_S_FORM",
    "SAFE_ARGUMENT_INTERVAL",
    "SOURCE_DEFINITION",
    "SIGMA_SOURCE",
    "SIGMA_PRIME_SOURCE",
    "ScopedSigmaOracleError",
    "AuthenticationError",
    "GraphRecordError",
    "BinderError",
    "IntervalInputError",
    "FractionInterval",
    "AuthenticatedSigmaBounds",
    "ScopedSigmaInterval",
    "EXPECTED_GRAPH_ROWS",
    "authenticated_bounds_from_result",
    "RawSigmaScopedOracle",
    "build_oracle",
    "make_oracle",
    "scoped_sigma_oracle",
]
