"""Same-graph C1 action for the reserved three-channel matrix.

This module is deliberately a small algebraic consumer.  It does not import
the reserved operator producer, read an archive, evaluate an interval, or
choose a value from an enclosure.  The caller supplies the matrix, its
ordinary-``Z`` derivative, a caller-owned inverse, and the three graph
operations used to combine source handles.

For a matrix ``M(Z)`` and a vector ``y(Z)`` the action implemented here is

``beta = M_inverse * y``
``beta_Z = M_inverse * (y_Z - M_Z * beta)``.

The inverse is an explicit caller premise.  This module checks dimensions and
same-graph ownership, but it deliberately does *not* multiply ``M`` by the
supplied inverse and never reports that the inverse identity was verified.
Likewise, the full third-row helper only promotes the fixed ``(0, 0, 1)``
response conditionally, after the caller supplies the same function handles
for the complete m15 row and the third matrix row (including their ``Z``
handles).  It never infers a function identity from equal values.

The physical m15 sign is intentionally outside this module.  The returned
third-row quantity is ``M[2, :] * beta`` with ordinary C1 derivative
``M_Z[2, :] * beta + M[2, :] * beta_Z``.  A caller representing signed m15
may apply the physical minus sign explicitly.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any, Final


SCHEMA: Final[str] = "reserved_matrix_c1_action_v1"
DIMENSION: Final[int] = 3
FULL_RESPONSE_ROW: Final[tuple[int, int, int]] = (0, 0, 1)
FULL_RESPONSE_ROW_Z: Final[tuple[int, int, int]] = (0, 0, 0)


class MatrixC1ActionError(ValueError):
    """Base error for malformed matrix, pair, or source-handle input."""


class SourceHandleError(MatrixC1ActionError):
    """Raised when a value is not a live same-graph source handle."""


class MatrixShapeError(MatrixC1ActionError):
    """Raised when a matrix/vector does not have the fixed three-channel shape."""


@dataclass(frozen=True, slots=True)
class C1Jet:
    """Minimal C1 pair used only when a caller does not provide a pair factory.

    ``value`` and ``Z`` remain the caller's source handles; this class never
    stores a scalar value or evaluates either handle.
    """

    value: Any
    Z: Any


@dataclass(frozen=True, slots=True)
class MatrixC1ActionReport:
    """Non-certifying provenance report for :func:`inverse_c1_action`.

    ``same_graph_handles`` means that the handles passed to this module and
    the handles returned by the supplied operations share one live graph object
    by ``is`` identity.  It is not an authentication of the producer or an
    archive.  ``inverse_identity_verified`` is permanently false in this
    module: the inverse is a caller premise.
    """

    schema: str
    dimension: int
    same_graph_handles: bool
    source_handles_bound: bool
    caller_inverse_assumed: bool
    inverse_identity_verified: bool
    actual_source_owner_authenticated: bool
    status: str
    derivative_formula: str


@dataclass(frozen=True, slots=True)
class MatrixC1Action:
    """C1 inverse-action result with both component rows and pair objects."""

    beta: tuple[Any, Any, Any]
    beta_Z: tuple[Any, Any, Any]
    beta_pairs: tuple[Any, Any, Any]
    report: MatrixC1ActionReport

    def __iter__(self):
        """Allow ``beta, beta_Z = result`` without hiding the report."""

        yield self.beta
        yield self.beta_Z


@dataclass(frozen=True, slots=True)
class FullThirdRowReport:
    """Conditional status for the complete m15 / third-row relationship."""

    schema: str
    third_row_c1_identity_checked: bool
    full_row_identity_conditionally_available: bool
    conditional_response_row: tuple[int, int, int] | None
    conditional_response_row_Z: tuple[int, int, int] | None
    row_identity_uses_function_handles: bool
    physical_m15_negative_sign_is_external: bool
    inverse_identity_verified: bool
    actual_source_owner_authenticated: bool
    status: str


@dataclass(frozen=True, slots=True)
class FullThirdRowC1Action:
    """Third matrix-row dot product and its ordinary C1 derivative."""

    response: Any
    beta: tuple[Any, Any, Any]
    beta_Z: tuple[Any, Any, Any]
    inverse_action: MatrixC1Action
    report: FullThirdRowReport

    @property
    def pair(self) -> Any:
        """Alias matching the ``C1Function`` style used by source callers."""

        return self.response


@dataclass(frozen=True, slots=True)
class _OwnerRef:
    kind: str
    value: Any
    attribute: str

    def matches(self, other: "_OwnerRef") -> bool:
        if self.kind != other.kind:
            return False
        return self.value is other.value


@dataclass(frozen=True, slots=True)
class _PreparedInputs:
    matrix: tuple[tuple[Any, Any, Any], ...]
    matrix_Z: tuple[tuple[Any, Any, Any], ...]
    inverse: tuple[tuple[Any, Any, Any], ...]
    y: tuple[Any, Any, Any]
    y_Z: tuple[Any, Any, Any]
    owner: _OwnerRef
    pair_template: Any | None


_SCALAR_LIKE: Final[tuple[type[Any], ...]] = (
    type(None),
    bool,
    int,
    float,
    complex,
    str,
    bytes,
    bytearray,
    list,
    tuple,
    dict,
    set,
    frozenset,
)
_OWNER_OBJECT_ATTRIBUTES: Final[tuple[str, ...]] = (
    "graph",
    "source_graph",
    "owner",
)
_FUNCTION_ID_ATTRIBUTES: Final[tuple[str, ...]] = (
    "node",
    "function_id",
    "handle_id",
    "source_node",
)


def _stable_equal(left: Any, right: Any) -> bool:
    """Compare an explicit deterministic ID without comparing expression values."""

    if type(left) is not type(right):
        return False
    if isinstance(left, _SCALAR_LIKE):
        return left == right
    return left is right


def _owner_ref(handle: Any, label: str) -> _OwnerRef:
    """Extract only a non-scalar live graph owner by object identity.

    A graph name, graph ID, node number, string token, or other scalar
    declaration is not an owner.  Repository ``FunctionRef`` values satisfy
    this through their non-scalar ``.graph`` object.
    """

    if isinstance(handle, _SCALAR_LIKE):
        raise SourceHandleError(f"{label} must be a source handle, not a scalar/container")
    for attribute in _OWNER_OBJECT_ATTRIBUTES:
        try:
            candidate = getattr(handle, attribute)
        except (AttributeError, TypeError):
            continue
        if candidate is None:
            continue
        if isinstance(candidate, _SCALAR_LIKE):
            # A scalar graph/name/token is deliberately not accepted as a
            # source owner.  Same-name or same-ID claims are not provenance.
            continue
        return _OwnerRef("object", candidate, attribute)
    raise SourceHandleError(
        f"{label} has no non-scalar live graph owner; graph_id, owner_id, "
        "names, numeric values, and intervals cannot stand in for source functions"
    )


def _require_scalar_handle(handle: Any, label: str, owner: _OwnerRef | None = None) -> _OwnerRef:
    """Validate one scalar source handle and optionally match a known owner."""

    if hasattr(handle, "value") and hasattr(handle, "Z"):
        raise SourceHandleError(f"{label} is a C1 pair; pass its value and Z components separately")
    candidate = _owner_ref(handle, label)
    if owner is not None and not owner.matches(candidate):
        raise SourceHandleError(f"{label} is not on the same source graph owner")
    return candidate


def _pair_components(pair: Any, label: str) -> tuple[Any, Any]:
    """Read a C1 pair without coercing either component to a scalar."""

    if not hasattr(pair, "value") or not hasattr(pair, "Z"):
        raise SourceHandleError(f"{label} must expose source handles as .value and .Z")
    value = getattr(pair, "value")
    derivative = getattr(pair, "Z")
    _require_scalar_handle(value, f"{label}.value")
    _require_scalar_handle(derivative, f"{label}.Z")
    return value, derivative


def _as_vector(value: Any, label: str) -> tuple[Any, Any, Any]:
    if isinstance(value, (str, bytes, bytearray)):
        raise MatrixShapeError(f"{label} must contain exactly three source handles")
    try:
        items = tuple(value)
    except TypeError as exc:
        raise MatrixShapeError(f"{label} must be a length-three sequence") from exc
    if len(items) != DIMENSION:
        raise MatrixShapeError(f"{label} must have length three")
    return items  # type: ignore[return-value]


def _as_matrix(value: Any, label: str) -> tuple[tuple[Any, Any, Any], ...]:
    if isinstance(value, (str, bytes, bytearray)):
        raise MatrixShapeError(f"{label} must be a three-by-three source matrix")
    try:
        rows = tuple(value)
    except TypeError as exc:
        raise MatrixShapeError(f"{label} must be a three-by-three sequence") from exc
    if len(rows) != DIMENSION:
        raise MatrixShapeError(f"{label} must have exactly three rows")
    result = []
    for i, row in enumerate(rows):
        items = _as_vector(row, f"{label}[{i}]")
        result.append(items)
    return tuple(result)  # type: ignore[return-value]


def _same_function(left: Any, right: Any, left_label: str, right_label: str) -> bool:
    """Check function identity by object identity or graph-owner/function ID only."""

    left_owner = _owner_ref(left, left_label)
    right_owner = _owner_ref(right, right_label)
    if not left_owner.matches(right_owner):
        return False
    if left is right:
        return True
    for attribute in _FUNCTION_ID_ATTRIBUTES:
        left_has = hasattr(left, attribute)
        right_has = hasattr(right, attribute)
        if not (left_has and right_has):
            continue
        left_id = getattr(left, attribute)
        right_id = getattr(right, attribute)
        if _stable_equal(left_id, right_id):
            return True
    return False


def same_source_function(left: Any, right: Any) -> bool:
    """Public non-numeric source-function identity check.

    Equal values, equal string representations, and equal interval endpoints
    are intentionally never examined.
    """

    return _same_function(left, right, "left", "right")


def _normalise_c1_vector(
    values: Any,
    derivatives: Any | None,
    label: str,
) -> tuple[tuple[Any, Any, Any], tuple[Any, Any, Any], Any | None]:
    """Split a vector of C1 pairs, or a value vector plus a derivative vector."""

    raw = _as_vector(values, label)
    if derivatives is None:
        components = [_pair_components(item, f"{label}[{i}]") for i, item in enumerate(raw)]
        return (
            tuple(item[0] for item in components),
            tuple(item[1] for item in components),
            raw[0],
        )  # type: ignore[return-value]
    raw_Z = _as_vector(derivatives, label + "_Z")
    for i, (value, derivative) in enumerate(zip(raw, raw_Z)):
        _require_scalar_handle(value, f"{label}[{i}]")
        _require_scalar_handle(derivative, f"{label}_Z[{i}]")
    return raw, raw_Z, None  # type: ignore[return-value]


def _normalise_matrix(
    values: Any,
    derivatives: Any | None,
    label: str,
) -> tuple[tuple[tuple[Any, Any, Any], ...], tuple[tuple[Any, Any, Any], ...]]:
    """Split a C1 matrix or a value matrix plus an explicit derivative matrix."""

    raw = _as_matrix(values, label)
    if derivatives is None:
        value_rows = []
        derivative_rows = []
        for i, row in enumerate(raw):
            components = [_pair_components(item, f"{label}[{i}][{j}]") for j, item in enumerate(row)]
            value_rows.append(tuple(item[0] for item in components))
            derivative_rows.append(tuple(item[1] for item in components))
        return tuple(value_rows), tuple(derivative_rows)
    raw_Z = _as_matrix(derivatives, label + "_Z")
    for i, row in enumerate(raw):
        for j, value in enumerate(row):
            _require_scalar_handle(value, f"{label}[{i}][{j}]")
            _require_scalar_handle(raw_Z[i][j], f"{label}_Z[{i}][{j}]")
    return raw, raw_Z


def _validate_callbacks(add: Callable[..., Any], mul: Callable[..., Any], neg: Callable[..., Any]) -> None:
    if not callable(add) or not callable(mul) or not callable(neg):
        raise TypeError("add, mul, and neg must be caller-supplied algebra operations")


def _check_owner_collection(
    owner: _OwnerRef | None,
    values: Sequence[Any],
    label: str,
) -> _OwnerRef:
    current = owner
    for i, value in enumerate(values):
        candidate = _require_scalar_handle(value, f"{label}[{i}]", current)
        if current is None:
            current = candidate
    if current is None:
        raise SourceHandleError(f"{label} is empty")
    return current


def _prepare_inputs(
    matrix: Any,
    inverse: Any,
    y: Any,
    y_Z: Any | None,
    matrix_Z: Any | None,
) -> _PreparedInputs:
    matrix_value, matrix_derivative = _normalise_matrix(matrix, matrix_Z, "matrix")
    inverse_value = _as_matrix(inverse, "inverse")
    y_value, y_derivative, pair_template = _normalise_c1_vector(y, y_Z, "y")

    owner: _OwnerRef | None = None
    for i, row in enumerate(matrix_value):
        owner = _check_owner_collection(owner, row, f"matrix[{i}]")
    for i, row in enumerate(matrix_derivative):
        owner = _check_owner_collection(owner, row, f"matrix_Z[{i}]")
    for i, row in enumerate(inverse_value):
        owner = _check_owner_collection(owner, row, f"inverse[{i}]")
    owner = _check_owner_collection(owner, y_value, "y")
    owner = _check_owner_collection(owner, y_derivative, "y_Z")
    return _PreparedInputs(
        matrix=matrix_value,
        matrix_Z=matrix_derivative,
        inverse=inverse_value,
        y=y_value,
        y_Z=y_derivative,
        owner=owner,
        pair_template=pair_template,
    )


def _ensure_result_owner(value: Any, owner: _OwnerRef, label: str) -> Any:
    _require_scalar_handle(value, label, owner)
    return value


def _add_many(
    terms: Sequence[Any],
    *,
    add: Callable[..., Any],
    owner: _OwnerRef,
    label: str,
) -> Any:
    if not terms:
        raise MatrixC1ActionError(f"{label} has no terms")
    result = terms[0]
    _ensure_result_owner(result, owner, f"{label}[0]")
    for i, term in enumerate(terms[1:], 1):
        _ensure_result_owner(term, owner, f"{label}[{i}]")
        try:
            result = add(result, term)
        except Exception as exc:
            raise MatrixC1ActionError(f"caller add operation failed for {label}") from exc
        _ensure_result_owner(result, owner, f"{label}.add[{i}]")
    return result


def _multiply(
    left: Any,
    right: Any,
    *,
    mul: Callable[..., Any],
    owner: _OwnerRef,
    label: str,
) -> Any:
    _ensure_result_owner(left, owner, f"{label}.left")
    _ensure_result_owner(right, owner, f"{label}.right")
    try:
        result = mul(left, right)
    except Exception as exc:
        raise MatrixC1ActionError(f"caller mul operation failed for {label}") from exc
    return _ensure_result_owner(result, owner, label)


def _negate(
    value: Any,
    *,
    neg: Callable[..., Any],
    owner: _OwnerRef,
    label: str,
) -> Any:
    _ensure_result_owner(value, owner, label + ".input")
    try:
        result = neg(value)
    except Exception as exc:
        raise MatrixC1ActionError(f"caller neg operation failed for {label}") from exc
    return _ensure_result_owner(result, owner, label)


def _dot(
    row: Sequence[Any],
    vector: Sequence[Any],
    *,
    add: Callable[..., Any],
    mul: Callable[..., Any],
    owner: _OwnerRef,
    label: str,
) -> Any:
    if len(row) != DIMENSION or len(vector) != DIMENSION:
        raise MatrixShapeError(f"{label} dot product must have three terms")
    products = [
        _multiply(row[i], vector[i], mul=mul, owner=owner, label=f"{label}.product[{i}]")
        for i in range(DIMENSION)
    ]
    return _add_many(products, add=add, owner=owner, label=label)


def _make_c1_pair(
    template: Any | None,
    value: Any,
    derivative: Any,
    *,
    pair_factory: Callable[[Any, Any], Any] | None,
    owner: _OwnerRef,
    label: str,
) -> Any:
    if pair_factory is not None:
        if not callable(pair_factory):
            raise TypeError("pair_factory must be callable")
        try:
            pair = pair_factory(value, derivative)
        except Exception as exc:
            raise MatrixC1ActionError(f"caller pair_factory failed for {label}") from exc
    elif template is not None:
        try:
            pair = type(template)(value, derivative)
        except Exception as exc:
            raise MatrixC1ActionError(
                f"{label} pair type cannot be reconstructed; pass pair_factory explicitly"
            ) from exc
    else:
        pair = C1Jet(value, derivative)
    pair_value, pair_derivative = _pair_components(pair, label)
    _require_scalar_handle(pair_value, label + ".value", owner)
    _require_scalar_handle(pair_derivative, label + ".Z", owner)
    return pair


def _normalise_row_identity(
    row: Any,
    row_Z: Any | None,
    label: str,
) -> tuple[tuple[Any, Any, Any], tuple[Any, Any, Any]]:
    values, derivatives, _ = _normalise_c1_vector(row, row_Z, label)
    return values, derivatives


def inverse_c1_action(
    matrix: Any,
    inverse: Any,
    y: Any,
    y_Z: Any | None = None,
    *,
    matrix_Z: Any | None = None,
    add: Callable[..., Any],
    mul: Callable[..., Any],
    neg: Callable[..., Any],
    pair_factory: Callable[[Any, Any], Any] | None = None,
) -> MatrixC1Action:
    """Apply a caller-supplied 3 by 3 inverse to a genuine C1 source vector.

    ``matrix`` may be a 3 by 3 matrix of C1 pairs, in which case ``matrix_Z``
    is omitted, or it may be a value matrix accompanied by the explicit
    ``matrix_Z`` matrix.  ``y`` follows the same convention: it may be a
    length-three vector of C1 pairs, or a value vector accompanied by ``y_Z``.

    The returned derivative is exactly
    ``inverse * (y_Z - matrix_Z * beta)``.  In particular, ``matrix_Z`` is
    never replaced by zero, even when the caller's matrix was built from a
    parameter such as ``mu(Z)``.  No numerical evaluation or inverse check is
    performed.
    """

    _validate_callbacks(add, mul, neg)
    prepared = _prepare_inputs(matrix, inverse, y, y_Z, matrix_Z)
    owner = prepared.owner

    beta = tuple(
        _dot(
            prepared.inverse[i],
            prepared.y,
            add=add,
            mul=mul,
            owner=owner,
            label=f"beta[{i}]",
        )
        for i in range(DIMENSION)
    )
    # Build the full right-hand side vector explicitly.  No synthetic zero
    # handle is fabricated: every component remains a caller-owned source
    # expression, including the complete M_Z beta product.
    rhs_components = []
    for j in range(DIMENSION):
        matrix_Z_beta = _dot(
            prepared.matrix_Z[j],
            beta,
            add=add,
            mul=mul,
            owner=owner,
            label=f"matrix_Z_beta[{j}]",
        )
        rhs_components.append(
            _add_many(
                (
                    prepared.y_Z[j],
                    _negate(
                        matrix_Z_beta,
                        neg=neg,
                        owner=owner,
                        label=f"negative_matrix_Z_beta[{j}]",
                    ),
                ),
                add=add,
                owner=owner,
                label=f"beta_Z_rhs[{j}]",
            )
        )
    beta_Z_terms = []
    for i in range(DIMENSION):
        beta_Z_terms.append(
            _dot(
                prepared.inverse[i],
                rhs_components,
                add=add,
                mul=mul,
                owner=owner,
                label=f"beta_Z[{i}]",
            )
        )
    beta_Z = tuple(beta_Z_terms)

    pair_template = prepared.pair_template
    beta_pairs = tuple(
        _make_c1_pair(
            pair_template,
            beta[i],
            beta_Z[i],
            pair_factory=pair_factory,
            owner=owner,
            label=f"beta_pair[{i}]",
        )
        for i in range(DIMENSION)
    )
    report = MatrixC1ActionReport(
        schema=SCHEMA,
        dimension=DIMENSION,
        same_graph_handles=True,
        source_handles_bound=True,
        caller_inverse_assumed=True,
        inverse_identity_verified=False,
        actual_source_owner_authenticated=False,
        status="conditional_same_graph_caller_inverse",
        derivative_formula="beta_Z=M_inverse*(y_Z-M_Z*beta)",
    )
    return MatrixC1Action(
        beta=beta, beta_Z=beta_Z, beta_pairs=beta_pairs, report=report
    )


def full_third_row_c1_action(
    matrix: Any,
    inverse: Any,
    y: Any,
    y_Z: Any | None = None,
    *,
    matrix_Z: Any | None = None,
    complete_m15_row: Any | None = None,
    complete_m15_row_Z: Any | None = None,
    add: Callable[..., Any],
    mul: Callable[..., Any],
    neg: Callable[..., Any],
    pair_factory: Callable[[Any, Any], Any] | None = None,
    require_full_row_identity: bool = False,
) -> FullThirdRowC1Action:
    """Return the same-matrix third-row action and its ordinary C1 derivative.

    The returned pair is

    ``row[2] * beta``
    ``row[2]_Z * beta + row[2] * beta_Z``.

    If ``complete_m15_row`` is supplied, it must be a vector of the *same
    source function handles* as ``matrix[2]``.  Its ``Z`` components must also
    be the same handles as ``matrix_Z[2]``.  The check uses object identity or
    an explicit graph-owner/function-node identity only; it never compares
    values.  When the check passes, the fixed inverse response row
    ``(0, 0, 1)`` and its zero ordinary derivative are reported as conditional
    on the caller's inverse premise.  The inverse identity itself remains
    unverified.

    The physical signed m15 minus sign is not applied here.  A caller must
    negate the returned pair explicitly when its physical convention requires
    ``-m15``.
    """

    action = inverse_c1_action(
        matrix,
        inverse,
        y,
        y_Z,
        matrix_Z=matrix_Z,
        add=add,
        mul=mul,
        neg=neg,
        pair_factory=pair_factory,
    )
    prepared = _prepare_inputs(matrix, inverse, y, y_Z, matrix_Z)
    owner = prepared.owner
    row_value = _dot(
        prepared.matrix[2],
        action.beta,
        add=add,
        mul=mul,
        owner=owner,
        label="third_row_beta",
    )
    row_derivative = _add_many(
        (
            _dot(
                prepared.matrix_Z[2],
                action.beta,
                add=add,
                mul=mul,
                owner=owner,
                label="third_row_Z_beta",
            ),
            _dot(
                prepared.matrix[2],
                action.beta_Z,
                add=add,
                mul=mul,
                owner=owner,
                label="third_row_beta_Z",
            ),
        ),
        add=add,
        owner=owner,
        label="third_row_derivative",
    )

    identity_checked = False
    if complete_m15_row is not None or complete_m15_row_Z is not None:
        if complete_m15_row is None and complete_m15_row_Z is not None:
            raise MatrixShapeError("complete_m15_row is required with complete_m15_row_Z")
        m15_value, m15_Z = _normalise_row_identity(
            complete_m15_row,
            complete_m15_row_Z,
            "complete_m15_row",
        )
        for j in range(DIMENSION):
            _require_scalar_handle(m15_value[j], f"complete_m15_row[{j}]", owner)
            _require_scalar_handle(m15_Z[j], f"complete_m15_row_Z[{j}]", owner)
            if not _same_function(
                m15_value[j],
                prepared.matrix[2][j],
                f"complete_m15_row[{j}]",
                f"matrix[2][{j}]",
            ):
                raise SourceHandleError(
                    "complete m15 row is not the same source function as matrix third row"
                )
            if not _same_function(
                m15_Z[j],
                prepared.matrix_Z[2][j],
                f"complete_m15_row_Z[{j}]",
                f"matrix_Z[2][{j}]",
            ):
                raise SourceHandleError(
                    "complete m15 row Z derivative is not the same source function as matrix third row Z"
                )
        identity_checked = True
    elif require_full_row_identity:
        raise SourceHandleError(
            "complete_m15_row and its Z derivative are required for the full-row identity"
        )

    response = _make_c1_pair(
        prepared.pair_template,
        row_value,
        row_derivative,
        pair_factory=pair_factory,
        owner=owner,
        label="full_third_row_response",
    )
    report = FullThirdRowReport(
        schema=SCHEMA,
        third_row_c1_identity_checked=identity_checked,
        full_row_identity_conditionally_available=identity_checked,
        conditional_response_row=FULL_RESPONSE_ROW if identity_checked else None,
        conditional_response_row_Z=FULL_RESPONSE_ROW_Z if identity_checked else None,
        row_identity_uses_function_handles=True,
        physical_m15_negative_sign_is_external=True,
        inverse_identity_verified=False,
        actual_source_owner_authenticated=False,
        status=(
            "conditional_full_third_row_same_handles_caller_inverse"
            if identity_checked
            else "third_row_action_only_full_row_identity_unavailable"
        ),
    )
    return FullThirdRowC1Action(
        response=response,
        beta=action.beta,
        beta_Z=action.beta_Z,
        inverse_action=action,
        report=report,
    )


# Explicit aliases make the module easy to consume beside the existing
# ``coefficient_action`` and ``inverse_normalized_action`` names.
apply_inverse_c1 = inverse_c1_action
matrix_c1_inverse_action = inverse_c1_action
full_third_row_action = full_third_row_c1_action


__all__ = [
    "C1Jet",
    "DIMENSION",
    "FULL_RESPONSE_ROW",
    "FULL_RESPONSE_ROW_Z",
    "FullThirdRowC1Action",
    "FullThirdRowReport",
    "MatrixC1Action",
    "MatrixC1ActionError",
    "MatrixC1ActionReport",
    "MatrixShapeError",
    "SCHEMA",
    "SourceHandleError",
    "apply_inverse_c1",
    "full_third_row_action",
    "full_third_row_c1_action",
    "inverse_c1_action",
    "matrix_c1_inverse_action",
    "same_source_function",
]
