"""A source-owner guarded bridge from the raw to the ``K_chi`` residual.

The bridge is an algebraic functional, not an estimate and not a source
producer.  With one *identity* shared source owner, the caller supplies the
same regular Green operator ``G`` for both residuals, the operator
``A = -2*chi*E51``, and the resolvent ``K = (I-G*A)**-1*G``.  The represented
functionals are

``d_raw = W - G(g + Q*W)``

``d_chi = W - K(g + (Q-A)*W)``

and, as an operator identity,

``(I-G*A)*d_chi = d_raw`` and ``d_chi = (I+K*A)*d_raw``.

At ordinary fixed ``x`` the caller must certify ``G_Z = 0``.  For a supplied
ordinary ``Z1`` row the bridge then returns

``d_chi_Z = (I+K*A)*(d_raw_Z + G*A_Z*d_chi)``

with ``A_Z = -2*chi_Z*E51``.  No derivative of ``G`` is manufactured here.
The returned source-functional record distinguishes
``Q(W_Z) = C0*W_Z + B1*W_ZZ`` from the differentiated source
``d/dZ(QW) = C0_Z*W + C0*W_Z + B1_Z*W_Z + B1*W_ZZ``.  The supplied
``d_raw_Z`` is required to contain all four differentiated terms; in
particular neither ``B1*W_ZZ`` nor the coefficient-derivative terms may be
dropped.

This module has no numerical backend and does not import ``mpmath``.  Values
are opaque caller fields.  The only arithmetic required of a field is the
caller-defined linear-operator application and field addition (or an
``algebra backend.add`` method).  In particular, the module does not coerce a
field to a scalar, choose a midpoint, merge interval errors, reset a boundary,
or turn a norm cap into a field.  A signed residual is retained as a signed
functional and is not assumed to be zero.

The source-owner and resolvent statements are caller premises.  Passing a
certificate records those premises; it does not verify the actual source.
Labels, source IDs, approximate interval equality, and equal-looking values
cannot establish owner identity.  The same owner object must be supplied by
identity (``is``).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any


SCHEMA = "n1_raw_kchi_residual_bridge_v1"
ORDINARY_Z_ORDERS = (0, 1)
SOURCE_FUNCTIONAL = "signed_residual_functional"
A_DEFINITION = "A=-2*chi*E51"
AZ_DEFINITION = "A_Z=-2*chi_Z*E51"
Q_DEFINITION = "Q=C0+B1*d/dZ"
Q_ON_W_Z_DEFINITION = "Q(W_Z)=C0*W_Z+B1*W_ZZ"
DQW_D_Z_DEFINITION = "d/dZ(QW)=C0_Z*W+C0*W_Z+B1_Z*W_Z+B1*W_ZZ"
# Kept as a public spelling for callers that used the shorter pre-v1 name;
# its value now explicitly denotes the full derivative, not Q(W_Z).
QW_Z_DEFINITION = DQW_D_Z_DEFINITION
RESOLVENT_DEFINITION = "(I-G*A)*K=G"

__all__ = [
    "SCHEMA",
    "ORDINARY_Z_ORDERS",
    "SOURCE_FUNCTIONAL",
    "A_DEFINITION",
    "AZ_DEFINITION",
    "Q_DEFINITION",
    "Q_ON_W_Z_DEFINITION",
    "DQW_D_Z_DEFINITION",
    "QW_Z_DEFINITION",
    "RESOLVENT_DEFINITION",
    "ResidualBridgeError",
    "ResidualBridgeContractError",
    "SourceOwnerContract",
    "build",
    "bridge",
    "bridge_residuals",
    "transform",
]


class ResidualBridgeError(ValueError):
    """Base exception for a rejected residual bridge contract."""


class ResidualBridgeContractError(ResidualBridgeError):
    """The caller did not provide a strict source/algebra contract."""


_MISSING = object()


@dataclass(frozen=True)
class SourceOwnerContract:
    """Small explicit contract helper for callers that prefer an object.

    ``owner`` is intentionally opaque.  The bridge compares it only by
    identity, never by a label or an equality operator.  The remaining fields
    are caller premises and therefore are not independently re-proved here.
    ``units``, ``coordinates``, and ``initial_data`` are required shared
    context records (or explicit matching raw/chi pairs); they are compared
    by the bridge with a caller match flag or one shared object.
    """

    owner: Any
    resolvent_verified: bool
    G_Z: Any = 0
    z1_provided: bool = True
    units: Any = _MISSING
    coordinates: Any = _MISSING
    initial_data: Any = _MISSING
    A_definition: str = A_DEFINITION
    Az_definition: str = AZ_DEFINITION
    q_definition: str = Q_DEFINITION
    q_on_w_z_definition: str = Q_ON_W_Z_DEFINITION
    d_qw_d_z_definition: str = DQW_D_Z_DEFINITION
    q_w_z_includes_B1_W_ZZ: bool = True
    d_qw_d_z_includes_C0_Z_W: bool = True
    d_qw_d_z_includes_B1_Z_W_Z: bool = True
    d_raw_Z_includes_full_dQW: bool = True


class _BackendOperator:
    """Callable adapter for a backend exposing ``apply(name, field)``."""

    def __init__(self, backend: Any, name: str) -> None:
        self.backend = backend
        self.name = name
        # An ``apply`` backend is a stronger explicit owner contract only
        # when it exposes one shared owner for all named operators.  Expose
        # that token on each wrapper so the ordinary per-operator check below
        # remains identity based.
        for attr in ("owner", "source_owner", "shared_owner", "shared_source_owner"):
            if hasattr(backend, attr):
                self.owner = getattr(backend, attr)
                break
        definition_names = (
            "definition",
            f"{name}_definition",
            f"{name}_formula",
            "operator_definition",
        )
        for attr in definition_names:
            if hasattr(backend, attr):
                self.definition = getattr(backend, attr)
                break

    def __call__(self, field: Any) -> Any:
        apply = getattr(self.backend, "apply", None)
        if not callable(apply):
            raise ResidualBridgeContractError(
                "algebra backend.apply(name, field) is not callable"
            )
        return apply(self.name, field)


def _contract_value(contract: Any, names: tuple[str, ...], default: Any = _MISSING) -> Any:
    """Read one contract field without treating a label as an owner token."""

    if isinstance(contract, Mapping):
        for name in names:
            if name in contract:
                return contract[name]
        return default
    for name in names:
        if hasattr(contract, name):
            return getattr(contract, name)
    return default


def _backend_value(backend: Any, names: tuple[str, ...], default: Any = _MISSING) -> Any:
    if backend is None:
        return default
    if isinstance(backend, Mapping):
        for name in names:
            if name in backend:
                return backend[name]
        nested = backend.get("operators", _MISSING)
        if isinstance(nested, Mapping):
            for name in names:
                if name in nested:
                    return nested[name]
        return default
    for name in names:
        if hasattr(backend, name):
            return getattr(backend, name)
    nested = getattr(backend, "operators", None)
    if isinstance(nested, Mapping):
        for name in names:
            if name in nested:
                return nested[name]
    return default


def _operator(
    explicit: Any,
    backend: Any,
    name: str,
) -> Any:
    if explicit is not _MISSING and explicit is not None:
        value = explicit
    else:
        value = _backend_value(backend, (name, name.lower()), _MISSING)
        if value is _MISSING and callable(_backend_value(backend, ("apply",), None)):
            value = _BackendOperator(backend, name)
    if value is _MISSING or value is None:
        raise ResidualBridgeContractError(
            f"callable linear operator {name} is required"
        )
    if not callable(value) and not callable(getattr(value, "__matmul__", None)):
        raise ResidualBridgeContractError(
            f"linear operator {name} must be callable (or provide __matmul__)"
        )
    return value


def _apply(operator: Any, field: Any, name: str) -> Any:
    """Apply an opaque caller operator without converting the field."""

    try:
        if callable(operator):
            return operator(field)
        return operator @ field
    except (TypeError, ValueError, ArithmeticError) as exc:
        raise ResidualBridgeError(
            f"operator {name} could not act on the caller field"
        ) from exc


def _backend_add(backend: Any, left: Any, right: Any) -> Any:
    add = _backend_value(backend, ("add", "field_add"), _MISSING)
    if add is not _MISSING:
        if not callable(add):
            raise ResidualBridgeContractError(
                "algebra backend add must be callable"
            )
        try:
            return add(left, right)
        except (TypeError, ValueError, ArithmeticError) as exc:
            raise ResidualBridgeError("algebra backend could not add fields") from exc
    try:
        return left + right
    except (TypeError, ValueError, ArithmeticError) as exc:
        raise ResidualBridgeError(
            "caller fields must support +, or backend.add must be supplied"
        ) from exc


def _same_metadata(left: Any, right: Any) -> bool:
    """Compare context records without requiring them to be numerics."""

    if left is right:
        return True
    try:
        result = left == right
    except Exception:
        return False
    return result is True or (isinstance(result, bool) and result)


def _pair_from_mapping(contract: Any, name: str) -> tuple[Any, Any]:
    """Find raw/chi metadata pair under common explicit spellings."""

    raw = chi = _MISSING
    if isinstance(contract, Mapping):
        for raw_name in (f"raw_{name}", f"raw{name.title()}"):
            if raw_name in contract:
                raw = contract[raw_name]
                break
        for chi_name in (f"chi_{name}", f"chi{name.title()}", f"kchi_{name}"):
            if chi_name in contract:
                chi = contract[chi_name]
                break
        for branch_name in ("raw", "raw_source", "raw_context"):
            branch = contract.get(branch_name)
            if isinstance(branch, Mapping) and name in branch:
                raw = branch[name]
                break
        for branch_name in ("chi", "kchi", "chi_source", "chi_context"):
            branch = contract.get(branch_name)
            if isinstance(branch, Mapping) and name in branch:
                chi = branch[name]
                break
    return raw, chi


def _check_context_metadata(
    contract: Any,
    *,
    units: Any,
    coordinates: Any,
    initial_data: Any,
) -> dict[str, Any]:
    """Require matching raw/chi units, chart, and initial-data records.

    Each context item must be explicitly supplied either as a raw/chi pair or
    as one shared object.  A pair is admissible only when values compare equal
    *and* the caller supplies a true ``same_*``/``*_match`` flag, or when the
    two values are literally the same object.  Thus a copied label or an
    approximate interval cannot manufacture source identity.
    """

    result: dict[str, Any] = {}
    supplied = {
        "units": units,
        "coordinates": coordinates,
        "initial_data": initial_data,
    }
    aliases = {
        "units": ("units", "unit_system", "normalization_units"),
        "coordinates": ("coordinates", "coordinate_system", "chart"),
        "initial_data": ("initial_data", "initial_condition", "initial_state"),
    }
    for name, values in supplied.items():
        shared = values
        if shared is _MISSING:
            shared = _contract_value(contract, aliases[name], _MISSING)
        raw, chi = _pair_from_mapping(contract, name)
        if shared is not _MISSING:
            # A shared field may itself be an explicit raw/chi pair mapping.
            if isinstance(shared, Mapping):
                if raw is _MISSING:
                    raw = shared.get("raw", shared.get("raw_value", _MISSING))
                if chi is _MISSING:
                    chi = shared.get("chi", shared.get("chi_value", _MISSING))
            if raw is _MISSING and chi is _MISSING:
                raw = chi = shared
        flag_names = {
            "units": ("same_units", "units_match"),
            "coordinates": ("same_coordinates", "coordinates_match", "same_chart"),
            "initial_data": ("same_initial_data", "initial_data_match", "same_initial_condition"),
        }[name]
        flag = _contract_value(contract, flag_names, _MISSING)
        if flag is not _MISSING and flag is not True:
            raise ResidualBridgeContractError(
                f"caller did not certify matching {name}"
            )
        if raw is _MISSING and chi is _MISSING:
            raise ResidualBridgeContractError(
                f"explicit raw/chi {name} context is required"
            )
        if raw is _MISSING or chi is _MISSING:
            raise ResidualBridgeContractError(
                f"both raw and chi {name} records are required"
            )
        if not _same_metadata(raw, chi):
            raise ResidualBridgeContractError(
                f"raw and chi {name} records differ"
            )
        if flag is not True and raw is not chi:
            raise ResidualBridgeContractError(
                f"raw and chi {name} need a true match certificate or one shared object"
            )
        result[name] = {
            "status": "matched",
            "raw": raw,
            "chi": chi,
            "explicit_match_certificate": flag is True,
        }
    return result


def _owner_from_contract(contract: Any) -> Any:
    return _contract_value(
        contract,
        ("owner", "source_owner", "shared_owner", "shared_source_owner"),
        _MISSING,
    )


def _require_owner(owner: Any, contract: Any) -> Any:
    contract_owner = _owner_from_contract(contract)
    if owner is _MISSING or owner is None:
        if contract_owner is _MISSING or contract_owner is None:
            raise ResidualBridgeContractError(
                "identity-based shared source owner is required"
            )
        owner = contract_owner
    if contract_owner is not _MISSING and contract_owner is not owner:
        raise ResidualBridgeContractError(
            "contract owner and supplied owner are not the same object"
        )
    if owner is None:
        raise ResidualBridgeContractError(
            "identity-based shared source owner may not be None"
        )
    return owner


def _require_resolvent(contract: Any) -> None:
    value = _contract_value(
        contract,
        (
            "resolvent_verified",
            "resolvent_relation_verified",
            "K_resolvent_verified",
            "resolvent_certificate",
        ),
        _MISSING,
    )
    if value is not True:
        raise ResidualBridgeContractError(
            f"caller certificate must assert {RESOLVENT_DEFINITION}"
        )


def _is_zero_marker(value: Any) -> bool:
    if isinstance(value, bool):
        return False
    if isinstance(value, str):
        return value.strip().lower() in {
            "0",
            "zero",
            "known_zero",
            "g_z=0",
            "g_z_zero",
        }
    # Do not accept float zero, Fraction/interval zero, or an arbitrary
    # user-defined object whose equality happens to return True.  Such a
    # value is not an explicit ordinary fixed-x certificate.
    return type(value) is int and value == 0


def _require_fixed_x_gz(contract: Any, explicit_g_z: Any) -> None:
    value = explicit_g_z
    if value is _MISSING or value is None:
        value = _contract_value(
            contract,
            ("G_Z", "g_z", "GZ", "g_z_status", "G_Z_status"),
            _MISSING,
        )
        if value is _MISSING:
            marker = _contract_value(
                contract,
                ("G_Z_zero", "g_z_zero", "fixed_x_G_Z_zero"),
                _MISSING,
            )
            if marker is True:
                return
            raise ResidualBridgeContractError(
                "ordinary fixed-x contract must certify known G_Z=0"
            )
    if not _is_zero_marker(value):
        raise ResidualBridgeContractError(
            "ordinary fixed-x bridge rejects unknown or nonzero G_Z"
        )


def _require_q_contract(contract: Any) -> None:
    q_definition = _contract_value(
        contract, ("q_definition", "Q_definition", "Q"), _MISSING
    )
    if q_definition != Q_DEFINITION:
        raise ResidualBridgeContractError(
            "Q must be the raw C0+B1*d/dZ differential operator"
        )
    q_on_w_z = _contract_value(
        contract,
        ("q_on_w_z_definition", "Q_on_W_Z_definition", "QW_Z_operator_definition"),
        _MISSING,
    )
    if q_on_w_z != Q_ON_W_Z_DEFINITION:
        raise ResidualBridgeContractError(
            "Q(W_Z) must be explicitly recorded as C0*W_Z+B1*W_ZZ"
        )
    d_qw_d_z = _contract_value(
        contract,
        (
            "d_qw_d_z_definition",
            "dQW_dZ_definition",
            "d_QW_d_Z_definition",
            "QW_Z_definition",
        ),
        _MISSING,
    )
    if d_qw_d_z != DQW_D_Z_DEFINITION:
        raise ResidualBridgeContractError(
            "d/dZ(QW) must record all C0_Z*W, C0*W_Z, B1_Z*W_Z, and B1*W_ZZ terms"
        )
    includes_second = _contract_value(
        contract,
        (
            "q_w_z_includes_B1_W_ZZ",
            "QW_Z_includes_B1_W_ZZ",
            "includes_B1_W_ZZ",
        ),
        _MISSING,
    )
    if includes_second is not True:
        raise ResidualBridgeContractError(
            "Q contract must explicitly retain the B1*W_ZZ term"
        )
    for names, label in (
        (
            ("d_qw_d_z_includes_C0_Z_W", "dQW_dZ_includes_C0_Z_W"),
            "C0_Z*W",
        ),
        (
            ("d_qw_d_z_includes_B1_Z_W_Z", "dQW_dZ_includes_B1_Z_W_Z"),
            "B1_Z*W_Z",
        ),
    ):
        included = _contract_value(contract, names, _MISSING)
        if included is not True:
            raise ResidualBridgeContractError(
                f"d/dZ(QW) contract must retain the {label} term"
            )
    full_raw_z = _contract_value(
        contract,
        ("d_raw_Z_includes_full_dQW", "d_raw_z_includes_full_dqw", "d_raw_Z_full_derivative"),
        _MISSING,
    )
    if full_raw_z is not True:
        raise ResidualBridgeContractError(
            "d_raw_Z must be explicitly certified to include all four d/dZ(QW) terms"
        )


def _require_z1_contract(contract: Any) -> None:
    """Reject an explicit certificate that disclaims the supplied Z1 row."""

    supplied = _contract_value(
        contract,
        ("z1_provided", "ordinary_Z1_provided", "Z1_provided"),
        _MISSING,
    )
    if supplied is not _MISSING and supplied is not True:
        raise ResidualBridgeContractError(
            "ordinary fixed-x source certificate must provide Z1"
        )


def _require_a_definition(contract: Any, operators: Mapping[str, Any]) -> None:
    definition = _contract_value(
        contract,
        ("A_definition", "a_definition", "A_formula", "a_formula"),
        _MISSING,
    )
    if definition != A_DEFINITION:
        raise ResidualBridgeContractError(
            "A must be the signed operator A=-2*chi*E51"
        )
    az_definition = _contract_value(
        contract,
        ("Az_definition", "A_Z_definition", "a_z_definition", "A_Z_formula"),
        _MISSING,
    )
    if az_definition != AZ_DEFINITION:
        raise ResidualBridgeContractError(
            "A_Z must be the signed derivative operator A_Z=-2*chi_Z*E51"
        )
    # The callable's descriptive source contract is required as an explicit
    # local sign check.  The relation itself remains a caller premise rather
    # than a live-source proof.
    required_definitions = {"A": A_DEFINITION, "Az": AZ_DEFINITION}
    for name in ("A", "Az"):
        operator = operators[name]
        operator_definition = _MISSING
        for attr in ("definition", "formula", "operator_definition"):
            if hasattr(operator, attr):
                operator_definition = getattr(operator, attr)
                break
        if operator_definition != required_definitions[name]:
            raise ResidualBridgeContractError(
                f"operator {name} must carry its explicit signed definition"
            )


def _require_operator_owner(owner: Any, backend: Any, operators: Mapping[str, Any]) -> None:
    backend_owner = _backend_value(
        backend, ("owner", "source_owner", "shared_owner", "shared_source_owner"), _MISSING
    )
    if backend_owner is not _MISSING and backend_owner is not owner:
        raise ResidualBridgeContractError(
            "algebra backend is attached to a different source owner"
        )
    for name, operator in operators.items():
        op_owner = _MISSING
        for attr in ("owner", "source_owner", "shared_owner", "shared_source_owner"):
            if hasattr(operator, attr):
                op_owner = getattr(operator, attr)
                break
        if op_owner is _MISSING:
            raise ResidualBridgeContractError(
                f"operator {name} must carry the shared source owner by identity"
            )
        if op_owner is not owner:
            raise ResidualBridgeContractError(
                f"operator {name} is attached to a different source owner"
            )


def _same_owner_aliases(owner: Any, source_owner: Any) -> Any:
    if owner is _MISSING and source_owner is not _MISSING:
        return source_owner
    if source_owner is not _MISSING and source_owner is not owner:
        raise ResidualBridgeContractError(
            "owner and source_owner aliases must name one identical object"
        )
    return owner


def bridge(
    d_raw: Any = None,
    d_raw_Z: Any = None,
    *,
    owner: Any = _MISSING,
    source_owner: Any = _MISSING,
    contract: Any = None,
    certificate: Any = None,
    source_owner_contract: Any = None,
    backend: Any = None,
    G: Any = _MISSING,
    K: Any = _MISSING,
    A: Any = _MISSING,
    Az: Any = _MISSING,
    A_Z: Any = _MISSING,
    A_z: Any = _MISSING,
    G_Z: Any = _MISSING,
    units: Any = _MISSING,
    coordinates: Any = _MISSING,
    initial_data: Any = _MISSING,
    norm_cap: Any = None,
) -> dict[str, Any]:
    """Return ``d_chi`` and its ordinary fixed-x ``Z1`` bridge.

    Parameters are intentionally caller-owned.  Either four explicit
    callables ``G,K,A,Az`` or a backend exposing those callables (or
    ``apply(name, field)``) is accepted.  A contract/certificate with an
    identity shared owner, the certified resolvent relation, and known
    ``G_Z=0`` is mandatory.  ``d_raw_Z`` is mandatory; passing ``None`` is
    rejected rather than treating a missing Z1 as zero.

    ``units``, ``coordinates``, and ``initial_data`` are convenience arguments
    for explicit context records, but all three matching records are required
    in the certificate or here.  Different units, charts, or initial data are
    rejected.  A norm cap is deliberately not an input field and raises
    immediately when supplied.

    No ``W``, ``g``, or ``Q`` value is evaluated here.  The output's
    ``source_functional`` records the exact raw and chi functionals, the
    operator rule ``Q(W_Z)=C0*W_Z+B1*W_ZZ``, and the full derivative rule
    ``d/dZ(QW)=C0_Z*W+C0*W_Z+B1_Z*W_Z+B1*W_ZZ``.  The caller's ``d_raw_Z``
    is required to include all four differentiated terms.
    ``certificate_is_caller_premise`` is always true and
    ``actual_source_verified`` is always false.
    """

    if norm_cap is not None:
        raise ResidualBridgeContractError(
            "norm cap is not a field of the residual bridge"
        )
    if d_raw is None:
        raise ResidualBridgeContractError("d_raw signed residual field is required")
    if d_raw_Z is None:
        raise ResidualBridgeContractError(
            "ordinary fixed-x d_raw_Z (Z1) signed residual field is required"
        )

    if source_owner_contract is not None:
        if contract is None:
            contract = source_owner_contract
        elif contract is not source_owner_contract:
            raise ResidualBridgeContractError(
                "contract and source_owner_contract were supplied twice"
            )
    if certificate is not None and contract is not None and certificate is not contract:
        # The certificate alias may be combined with an explicit contract
        # only when both mappings carry one same identity owner.
        contract_owner = _owner_from_contract(contract)
        certificate_owner = _owner_from_contract(certificate)
        if (
            contract_owner is _MISSING
            or certificate_owner is _MISSING
            or contract_owner is not certificate_owner
        ):
            raise ResidualBridgeContractError(
                "contract and certificate must carry one identical owner"
            )
        if isinstance(contract, Mapping) and isinstance(certificate, Mapping):
            merged = dict(contract)
            merged.update(certificate)
            contract = merged
        else:
            raise ResidualBridgeContractError(
                "supply one contract/certificate object or two mappings"
            )
    elif contract is not None and certificate is None:
        pass
    elif certificate is not None:
        contract = certificate
    if A_Z is not _MISSING and A_z is not _MISSING and A_Z is not A_z:
        raise ResidualBridgeContractError("A_Z and A_z were supplied twice")
    if A_Z is not _MISSING and Az is not _MISSING and A_Z is not Az:
        raise ResidualBridgeContractError("Az and A_Z were supplied twice")
    if A_z is not _MISSING and Az is not _MISSING and A_z is not Az:
        raise ResidualBridgeContractError("Az and A_z were supplied twice")
    if Az is _MISSING:
        if A_Z is not _MISSING:
            Az = A_Z
        elif A_z is not _MISSING:
            Az = A_z
    if contract is None:
        raise ResidualBridgeContractError(
            "caller source-owner/resolvent certificate is required"
        )

    owner = _same_owner_aliases(owner, source_owner)
    owner = _require_owner(owner, contract)
    _require_resolvent(contract)
    _require_fixed_x_gz(contract, G_Z)
    _require_z1_contract(contract)
    _require_q_contract(contract)

    operators = {
        "G": _operator(G, backend, "G"),
        "K": _operator(K, backend, "K"),
        "A": _operator(A, backend, "A"),
        "Az": _operator(Az, backend, "Az"),
    }
    _require_operator_owner(owner, backend, operators)
    _require_a_definition(contract, operators)
    context = _check_context_metadata(
        contract,
        units=units,
        coordinates=coordinates,
        initial_data=initial_data,
    )

    # d_chi = (I + K*A)d_raw.  Every operation remains in the caller's field
    # algebra; no identity matrix or scalar conversion is constructed.
    a_d_raw = _apply(operators["A"], d_raw, "A")
    k_a_d_raw = _apply(operators["K"], a_d_raw, "K")
    d_chi = _backend_add(backend, d_raw, k_a_d_raw)

    # d_chi_Z = (I + K*A)(d_raw_Z + G*A_Z*d_chi), with G_Z=0 certified above.
    az_d_chi = _apply(operators["Az"], d_chi, "Az")
    g_az_d_chi = _apply(operators["G"], az_d_chi, "G")
    z1_source = _backend_add(backend, d_raw_Z, g_az_d_chi)
    a_z1_source = _apply(operators["A"], z1_source, "A")
    k_a_z1_source = _apply(operators["K"], a_z1_source, "K")
    d_chi_Z = _backend_add(backend, z1_source, k_a_z1_source)

    source_functional = {
        "kind": SOURCE_FUNCTIONAL,
        "raw": "d_raw=W-G(g+QW)",
        "chi": "d_chi=W-K(g+(Q-A)W)",
        "Q": Q_DEFINITION,
        "Q_on_W_Z": Q_ON_W_Z_DEFINITION,
        "dQW_dZ": DQW_D_Z_DEFINITION,
        "QW_Z": DQW_D_Z_DEFINITION,
        "Q(W_Z)_includes_B1_W_ZZ": True,
        "dQW_dZ_includes_C0_Z_W": True,
        "dQW_dZ_includes_C0_W_Z": True,
        "dQW_dZ_includes_B1_Z_W_Z": True,
        "dQW_dZ_includes_B1_W_ZZ": True,
        "d_raw_Z_includes_full_dQW": True,
        "A": A_DEFINITION,
        "A_Z": AZ_DEFINITION,
        "K": "K=(I-G*A)^-1*G",
        "ordinary_fixed_x": "G_Z=0",
        "d_chi_identity": "d_chi=(I+K*A)d_raw",
        "d_chi_Z_identity": "d_chi_Z=(I+K*A)(d_raw_Z+G*A_Z*d_chi)",
        "uses_W_ZZ_in_QW_Z": True,
        "uses_W_ZZ_in_Q_on_W_Z": True,
        "uses_W_ZZ_in_dQW_dZ": True,
        "signed_residual": True,
        "is_value": False,
        "is_error_bound": False,
    }
    return {
        "schema": SCHEMA,
        "d_chi": d_chi,
        "d_chi_Z": d_chi_Z,
        "source_functional": source_functional,
        "source_owner": owner,
        "certificate": {
            "shared_source_owner_identity": True,
            "resolvent_relation": RESOLVENT_DEFINITION,
            "resolvent_relation_verified": True,
            "G_Z": 0,
            "certificate_is_caller_premise": True,
            "actual_source_verified": False,
            "context": context,
        },
        "signed_residual_contract": {
            "same_signed_residual": True,
            "zero_assumed": False,
            "boundary_reset": False,
            "interval_errors_merged": False,
            "norm_cap_is_field": False,
        },
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "operators_are_caller_supplied": True,
        "algebra_backend_used": backend is not None,
    }


def build(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Alias for :func:`bridge`, retained as the primary construction name."""

    return bridge(*args, **kwargs)


def bridge_residuals(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Readable alias for :func:`bridge`."""

    return bridge(*args, **kwargs)


def transform(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Readable alias for :func:`bridge`."""

    return bridge(*args, **kwargs)

