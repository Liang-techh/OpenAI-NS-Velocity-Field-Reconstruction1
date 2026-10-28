"""Static-only draft for finite native ``n=1`` Neumann words.

This v2 module is restricted to the common-core, fixed-outer-``x`` Green
contract.  The caller must explicitly declare
``green_scope="common_core_fixed_x"``, ``green_kernel_Z_independent=True``
and ``preconditioned_green=True``.  Missing or false declarations are
rejected.  They are caller premises only: this module cannot authenticate
that a callback really has no ``Z``-dependent kernel, no collar extension, or
the stated ``(I+J)G`` preconditioning, and therefore marks them unverified in
the result metadata.

The caller supplies ordinary-``Z`` jet callbacks ``g4_jets``, ``g5_jets`` and
``g6_jets`` and ordinary-``Z`` matrix callbacks ``C0_jets`` and ``B1_jets``.
Every family must contain orders ``0 .. K+d``.  A callback has shape
``callback(s, Z)`` and is never invoked while a plan is constructed.  The
principal Green callback has the planned response-interface shape::

    principal_green_response(
        *, x, Z, g4, g5, g6, integrate
    ) -> six_vector_or_object_with_W0

The callback may close over its logarithm and kernel implementation.  Its
source rows are lazy ``(s, Z)`` callbacks; ``integrate`` remains caller-owned.
The outer evaluation point ``x`` and ordinary ``Z`` are fixed.  Recursive
word evaluations at a Green kernel's radial coordinate ``s`` are internal
Volterra calls, not eta derivatives and not implicit ``Z`` derivatives.

Only words ``0 .. K`` are materialized.  The omitted infinite tail, source
identity, and quadrature error remain uncertified.  Repeating ``evaluate``
repeats the caller's Green/integration work; no cache is implied.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, TypeAlias


COMPONENT_ORDER = ("F1", "Uz1", "K1", "P1", "F1x", "Uz1x")
COMPONENT_MAPPING_ALIASES = {
    1: ("F1",),
    2: ("Uz1",),
    3: ("K1",),
    4: ("P1",),
    5: ("F1x", "d_x_F1"),
    6: ("Uz1x", "d_x_Uz1"),
}
SOURCE_ROWS = (4, 5, 6)
JET_COORDINATE = "ordinary Z"
GREEN_SCOPE = "common_core_fixed_x"
SCHEMA = "native_O3_eta_n1_volterra_words_draft_v2"

FieldJet: TypeAlias = Callable[[Any, Any], Any]
MatrixJet: TypeAlias = Callable[[Any, Any], Any]
IntegralOperator: TypeAlias = Callable[[Callable[[Any], Any], Any, Any], Any]
PrincipalGreen: TypeAlias = Callable[..., Any]


class VolterraWordsInputError(ValueError):
    """Raised before evaluation when a finite contract is absent or invalid."""


@dataclass(frozen=True)
class FiniteNeumannTerm:
    """One materialized finite word and its requested ordinary-Z jets."""

    word_index: int
    jet_values: tuple[tuple[Any, ...], ...]
    metadata: Mapping[str, Any]


@dataclass(frozen=True)
class FiniteNeumannEvaluation:
    """Finite words and finite partial sums; never the complete ``W``."""

    x: Any
    Z: Any
    K: int
    d: int
    required_jet_order: int
    terms: tuple[FiniteNeumannTerm, ...]
    partial_sum: tuple[tuple[Any, ...], ...]
    quadrature_error: Any
    tail_status: str
    metadata: Mapping[str, Any]


@dataclass(frozen=True)
class FiniteNeumannWordPlan:
    """Validated lazy word graph for one fixed outer ``(x, Z)`` point."""

    x: Any
    Z: Any
    K: int
    d: int
    required_jet_order: int
    green_scope: str
    green_kernel_Z_independent: bool
    preconditioned_green: bool
    _word_providers: tuple[tuple[FieldJet, ...], ...]
    _integrate: IntegralOperator
    _principal_green_response: PrincipalGreen
    _quadrature_error: Any

    def evaluate(self) -> FiniteNeumannEvaluation:
        """Materialize exactly words ``0`` through ``K`` at stored ``(x, Z)``.

        This method is intentionally uncached.  A second call can repeat
        every nested source callback and quadrature operation.
        """

        terms = []
        for word_index, providers in enumerate(self._word_providers):
            values = tuple(
                _green_vector(
                    provider(self.x, self.Z),
                    "W[{}][{}]".format(word_index, jet_order),
                )
                for jet_order, provider in enumerate(providers[: self.d + 1])
            )
            terms.append(
                FiniteNeumannTerm(
                    word_index=word_index,
                    jet_values=values,
                    metadata={
                        "word": "W[{}]".format(word_index),
                        "ordinary_Z_jets": True,
                        "eta_derivatives_used": False,
                        "finite_word_only": True,
                        "source_identity_status": "caller_owned_unverified",
                    },
                )
            )

        finite_terms = tuple(terms)
        partial_sum = tuple(
            _sum_vectors(
                (term.jet_values[jet_order] for term in finite_terms),
                "finite partial sum jet {}".format(jet_order),
            )
            for jet_order in range(self.d + 1)
        )
        error_status = (
            "caller_provided_unvalidated"
            if self._quadrature_error is not None
            else "not_provided"
        )
        return FiniteNeumannEvaluation(
            x=self.x,
            Z=self.Z,
            K=self.K,
            d=self.d,
            required_jet_order=self.required_jet_order,
            terms=finite_terms,
            partial_sum=partial_sum,
            quadrature_error=self._quadrature_error,
            tail_status="uncertified_missing_infinite_tail",
            metadata={
                "schema": SCHEMA,
                "recurrence": (
                    "W[k+1][j]=G_tilde sum_i binom(j,i) "
                    "(C0[i] W[k][j-i] + B1[i] W[k][j-i+1])"
                ),
                "finite_term_indices": tuple(range(self.K + 1)),
                "finite_partial_sum_only": True,
                "required_source_and_matrix_jet_order": self.required_jet_order,
                "source_and_matrix_order_checked_before_evaluation": True,
                "infinite_tail_status": "uncertified_missing_infinite_tail",
                "complete_W": False,
                "n_ge_2_temporal_recursion": False,
                "jet_coordinate": JET_COORDINATE,
                "outer_evaluation_x_fixed": True,
                "outer_evaluation_Z_fixed": True,
                "green_scope": self.green_scope,
                "green_kernel_Z_independent_declared": True,
                "green_kernel_and_endpoints_Z_independent_assumed": True,
                "green_kernel_and_endpoints_Z_independence_authenticated": False,
                "preconditioned_green_declared": True,
                "green_assumptions_authenticated": False,
                "collar_green_rejected_by_contract": True,
                "eta_derivatives_used": False,
                "green_internal_coordinate": (
                    "caller-supplied radial integration coordinate s; "
                    "not an eta derivative and not a Z derivative"
                ),
                "quadrature_error_status": error_status,
                "quadrature_error_certified_here": False,
                "source_identity_authenticated": False,
                "repeat_evaluate_reintegrates": True,
            },
        )


def prepare_finite_neumann_words(
    *,
    x: Any,
    Z: Any,
    K: int,
    d: int,
    g4_jets: Any,
    g5_jets: Any,
    g6_jets: Any,
    C0_jets: Any,
    B1_jets: Any,
    principal_green_response: PrincipalGreen,
    integrate: IntegralOperator,
    green_scope: str | None,
    green_kernel_Z_independent: bool | None,
    preconditioned_green: bool | None,
    quadrature_error: Any = None,
) -> FiniteNeumannWordPlan:
    """Validate the common-core Green contract and build a lazy word graph.

    The three Green assumptions are required keyword inputs with no accepting
    default.  They are declarations, not proofs: the callback is not
    inspected and no source identity or kernel-independence certificate is
    manufactured here.
    """

    _require_finite_nonnegative_int(K, "K")
    _require_finite_nonnegative_int(d, "d")
    _validate_green_assumptions(
        green_scope=green_scope,
        green_kernel_Z_independent=green_kernel_Z_independent,
        preconditioned_green=preconditioned_green,
    )
    required_order = K + d
    for name, value in (
        ("principal_green_response", principal_green_response),
        ("integrate", integrate),
    ):
        if not callable(value):
            raise VolterraWordsInputError(name + " must be callable")

    source_families = (
        _validate_field_family(g4_jets, required_order, "g4_jets"),
        _validate_field_family(g5_jets, required_order, "g5_jets"),
        _validate_field_family(g6_jets, required_order, "g6_jets"),
    )
    c0_family = _validate_matrix_family(C0_jets, required_order, "C0_jets")
    b1_family = _validate_matrix_family(B1_jets, required_order, "B1_jets")

    word_providers = []
    for word_index in range(K + 1):
        current_max = required_order - word_index
        current = []
        for jet_order in range(current_max + 1):
            if word_index == 0:
                current.append(
                    _make_source_word_provider(
                        jet_order=jet_order,
                        source_families=source_families,
                        principal_green_response=principal_green_response,
                        integrate=integrate,
                    )
                )
            else:
                current.append(
                    _make_feedback_word_provider(
                        word_index=word_index,
                        jet_order=jet_order,
                        previous=word_providers[word_index - 1],
                        c0_family=c0_family,
                        b1_family=b1_family,
                        principal_green_response=principal_green_response,
                        integrate=integrate,
                    )
                )
        word_providers.append(tuple(current))

    return FiniteNeumannWordPlan(
        x=x,
        Z=Z,
        K=K,
        d=d,
        required_jet_order=required_order,
        green_scope=green_scope,
        green_kernel_Z_independent=True,
        preconditioned_green=True,
        _word_providers=tuple(word_providers),
        _integrate=integrate,
        _principal_green_response=principal_green_response,
        _quadrature_error=quadrature_error,
    )


def compute_finite_neumann_words(
    *,
    x: Any,
    Z: Any,
    K: int,
    d: int,
    g4_jets: Any,
    g5_jets: Any,
    g6_jets: Any,
    C0_jets: Any,
    B1_jets: Any,
    principal_green_response: PrincipalGreen,
    integrate: IntegralOperator,
    green_scope: str | None,
    green_kernel_Z_independent: bool | None,
    preconditioned_green: bool | None,
    quadrature_error: Any = None,
) -> FiniteNeumannEvaluation:
    """Build and materialize finite words at the fixed outer ``(x, Z)``."""

    return prepare_finite_neumann_words(
        x=x,
        Z=Z,
        K=K,
        d=d,
        g4_jets=g4_jets,
        g5_jets=g5_jets,
        g6_jets=g6_jets,
        C0_jets=C0_jets,
        B1_jets=B1_jets,
        principal_green_response=principal_green_response,
        integrate=integrate,
        green_scope=green_scope,
        green_kernel_Z_independent=green_kernel_Z_independent,
        preconditioned_green=preconditioned_green,
        quadrature_error=quadrature_error,
    ).evaluate()


def finite_neumann_words(**kwargs: Any) -> FiniteNeumannEvaluation:
    """Compatibility spelling for :func:`compute_finite_neumann_words`."""

    return compute_finite_neumann_words(**kwargs)


def _validate_green_assumptions(
    *,
    green_scope: str | None,
    green_kernel_Z_independent: bool | None,
    preconditioned_green: bool | None,
) -> None:
    if green_scope != GREEN_SCOPE:
        raise VolterraWordsInputError(
            "green_scope must be explicitly declared as {}".format(GREEN_SCOPE)
        )
    if green_kernel_Z_independent is not True:
        raise VolterraWordsInputError(
            "green_kernel_Z_independent=True is a required caller declaration"
        )
    if preconditioned_green is not True:
        raise VolterraWordsInputError(
            "preconditioned_green=True is a required caller declaration"
        )


def _make_source_word_provider(
    *,
    jet_order: int,
    source_families: tuple[tuple[FieldJet, ...], ...],
    principal_green_response: PrincipalGreen,
    integrate: IntegralOperator,
) -> FieldJet:
    g4 = source_families[0][jet_order]
    g5 = source_families[1][jet_order]
    g6 = source_families[2][jet_order]

    def evaluate(radial_x: Any, ordinary_Z: Any) -> tuple[Any, ...]:
        value = principal_green_response(
            x=radial_x,
            Z=ordinary_Z,
            g4=g4,
            g5=g5,
            g6=g6,
            integrate=integrate,
        )
        return _green_vector(value, "source Green response jet {}".format(jet_order))

    return evaluate


def _make_feedback_word_provider(
    *,
    word_index: int,
    jet_order: int,
    previous: tuple[FieldJet, ...],
    c0_family: tuple[MatrixJet, ...],
    b1_family: tuple[MatrixJet, ...],
    principal_green_response: PrincipalGreen,
    integrate: IntegralOperator,
) -> FieldJet:
    row_callbacks = tuple(
        _make_feedback_row_callback(
            word_index=word_index,
            jet_order=jet_order,
            row=row,
            previous=previous,
            c0_family=c0_family,
            b1_family=b1_family,
        )
        for row in SOURCE_ROWS
    )
    g4, g5, g6 = row_callbacks

    def evaluate(radial_x: Any, ordinary_Z: Any) -> tuple[Any, ...]:
        value = principal_green_response(
            x=radial_x,
            Z=ordinary_Z,
            g4=g4,
            g5=g5,
            g6=g6,
            integrate=integrate,
        )
        return _green_vector(
            value,
            "feedback Green response W[{}][{}]".format(word_index, jet_order),
        )

    return evaluate


def _make_feedback_row_callback(
    *,
    word_index: int,
    jet_order: int,
    row: int,
    previous: tuple[FieldJet, ...],
    c0_family: tuple[MatrixJet, ...],
    b1_family: tuple[MatrixJet, ...],
) -> FieldJet:
    def feedback_row(radial_s: Any, ordinary_Z: Any) -> Any:
        total = None
        for derivative_order in range(jet_order + 1):
            w = _green_vector(
                previous[jet_order - derivative_order](radial_s, ordinary_Z),
                "W[{}][{}]".format(
                    word_index - 1, jet_order - derivative_order
                ),
            )
            w_next = _green_vector(
                previous[jet_order - derivative_order + 1](
                    radial_s, ordinary_Z
                ),
                "W[{}][{}]".format(
                    word_index - 1, jet_order - derivative_order + 1
                ),
            )
            c0_value = c0_family[derivative_order](radial_s, ordinary_Z)
            b1_value = b1_family[derivative_order](radial_s, ordinary_Z)
            c0_part = _matrix_row_product(c0_value, row, w)
            b1_part = _matrix_row_product(b1_value, row, w_next)
            contribution = c0_part + b1_part
            contribution = _binomial(jet_order, derivative_order) * contribution
            total = contribution if total is None else total + contribution
        return 0 if total is None else total

    return feedback_row


def _validate_field_family(
    family: Any,
    required_order: int,
    label: str,
) -> tuple[FieldJet, ...]:
    values = []
    for order in range(required_order + 1):
        value = _family_entry(family, order, label)
        if not callable(value):
            raise VolterraWordsInputError(
                "{}[{}] must be a callback of (s, Z)".format(label, order)
            )
        values.append(value)
    return tuple(values)


def _validate_matrix_family(
    family: Any,
    required_order: int,
    label: str,
) -> tuple[MatrixJet, ...]:
    values = []
    for order in range(required_order + 1):
        value = _family_entry(family, order, label)
        if not callable(value):
            raise VolterraWordsInputError(
                "{}[{}] must be a matrix callback of (s, Z)".format(
                    label, order
                )
            )
        values.append(value)
    return tuple(values)


def _family_entry(family: Any, order: int, label: str) -> Any:
    if isinstance(family, Mapping):
        if order not in family:
            raise VolterraWordsInputError(
                "{} is missing required ordinary-Z derivative order {}".format(
                    label, order
                )
            )
        return family[order]
    if isinstance(family, Sequence) and not isinstance(
        family, (str, bytes, bytearray)
    ):
        if order >= len(family):
            raise VolterraWordsInputError(
                "{} provides no ordinary-Z derivative order {}".format(
                    label, order
                )
            )
        return family[order]
    raise VolterraWordsInputError(
        "{} must be a mapping or sequence indexed by ordinary-Z order".format(
            label
        )
    )


def _require_finite_nonnegative_int(value: Any, label: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise VolterraWordsInputError(
            "{} must be a finite non-negative integer".format(label)
        )


def _green_vector(value: Any, label: str) -> tuple[Any, ...]:
    if isinstance(value, Mapping) and "W0" in value:
        value = value["W0"]
    else:
        missing = object()
        object_value = getattr(value, "W0", missing)
        if object_value is not missing:
            value = object_value
    if isinstance(value, Mapping):
        result = []
        for index in range(1, 7):
            result.append(_mapping_vector_component(value, index, label))
        return tuple(result)
    if isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    ):
        if len(value) != 6:
            raise VolterraWordsInputError(
                "{} must contain six components in {} order".format(
                    label, COMPONENT_ORDER
                )
            )
        return tuple(value)
    raise VolterraWordsInputError(
        "{} must be a six-vector or an object exposing W0".format(label)
    )


def _mapping_vector_component(
    mapping: Mapping[Any, Any], index: int, label: str
) -> Any:
    keys = [index, str(index)] + list(COMPONENT_MAPPING_ALIASES[index])
    present = [key for key in keys if key in mapping]
    if not present:
        raise VolterraWordsInputError(
            "{} is missing six-vector component {}".format(
                label, COMPONENT_ORDER[index - 1]
            )
        )
    first_key = present[0]
    first_value = mapping[first_key]
    for key in present[1:]:
        other_value = mapping[key]
        if other_value is first_value:
            continue
        try:
            equal = other_value == first_value
        except Exception as exc:
            raise VolterraWordsInputError(
                "conflicting aliases for {} cannot be compared".format(
                    COMPONENT_ORDER[index - 1]
                )
            ) from exc
        if not isinstance(equal, bool) or not equal:
            raise VolterraWordsInputError(
                "conflicting mapping aliases for {}".format(
                    COMPONENT_ORDER[index - 1]
                )
            )
    return first_value


def _matrix_row_product(matrix: Any, row: int, vector: tuple[Any, ...]) -> Any:
    """Multiply one row after strict dense/sparse matrix validation."""

    kind = _validate_matrix_representation(matrix)
    if isinstance(matrix, Mapping):
        if kind == "sparse":
            terms = []
            for key, coefficient in matrix.items():
                row_key, column_key = _sparse_key_parts(key)
                if row_key == row:
                    terms.append(coefficient * _vector_component(vector, column_key))
            return _sum_scalars(terms)
        row_value = _mapping_component(matrix, row)
        return _row_product(row_value, vector)

    return _row_product(matrix[row - 1], vector)


def _validate_matrix_representation(matrix: Any) -> str:
    if isinstance(matrix, Mapping):
        if not matrix:
            return "sparse"
        sparse_keys = []
        dense_keys = []
        for key in matrix:
            sparse = _sparse_key_parts(key)
            if sparse is not None:
                if not all(1 <= value <= 6 for value in sparse):
                    raise VolterraWordsInputError(
                        "sparse feedback matrix key must use row,column in 1..6"
                    )
                sparse_keys.append(sparse)
                continue
            dense = _dense_index_key(key)
            if dense is None:
                raise VolterraWordsInputError(
                    "invalid feedback matrix key; expected row,column sparse key "
                    "or dense row key 1..6"
                )
            dense_keys.append(dense)
        if sparse_keys and dense_keys:
            raise VolterraWordsInputError(
                "feedback matrix cannot mix sparse row,column and dense row keys"
            )
        if sparse_keys:
            return "sparse"
        if len(set(dense_keys)) != len(dense_keys):
            raise VolterraWordsInputError(
                "feedback matrix has duplicate dense row aliases"
            )
        for row_value in matrix.values():
            _validate_dense_row(row_value)
        return "dense"

    if isinstance(matrix, Sequence) and not isinstance(
        matrix, (str, bytes, bytearray)
    ):
        if len(matrix) != 6:
            raise VolterraWordsInputError("feedback matrix must have six rows")
        for row_value in matrix:
            _validate_dense_row(row_value)
        return "dense"

    raise VolterraWordsInputError(
        "feedback matrix must be dense six by six or sparse row,column mapping"
    )


def _validate_dense_row(row_value: Any) -> None:
    if isinstance(row_value, Mapping):
        if not row_value:
            return
        seen = set()
        for key in row_value:
            if _sparse_key_parts(key) is not None:
                raise VolterraWordsInputError(
                    "dense feedback row cannot contain sparse row,column keys"
                )
            column = _dense_index_key(key)
            if column is None:
                raise VolterraWordsInputError(
                    "dense feedback row key must be a column in 1..6"
                )
            if column in seen:
                raise VolterraWordsInputError(
                    "dense feedback row has duplicate column aliases"
                )
            seen.add(column)
        return
    if isinstance(row_value, Sequence) and not isinstance(
        row_value, (str, bytes, bytearray)
    ):
        if len(row_value) != 6:
            raise VolterraWordsInputError("feedback matrix row must have six columns")
        return
    raise VolterraWordsInputError("feedback matrix row must be a six-vector mapping")


def _row_product(row_value: Any, vector: tuple[Any, ...]) -> Any:
    if isinstance(row_value, Mapping):
        return _sum_scalars(
            coefficient * _vector_component(vector, _dense_index_key(key))
            for key, coefficient in row_value.items()
        )
    return _sum_scalars(
        coefficient * _vector_component(vector, column)
        for column, coefficient in enumerate(row_value, start=1)
    )


def _sparse_key_parts(key: Any) -> tuple[int, int] | None:
    if isinstance(key, str):
        pieces = key.split(",")
        if len(pieces) != 2:
            return None
        try:
            values = (int(pieces[0]), int(pieces[1]))
        except ValueError:
            return None
        if any(piece != str(value) for piece, value in zip(pieces, values)):
            return None
        return values
    if isinstance(key, tuple) and len(key) == 2:
        if all(isinstance(value, int) and not isinstance(value, bool) for value in key):
            return key[0], key[1]
    return None


def _dense_index_key(key: Any) -> int | None:
    if isinstance(key, int) and not isinstance(key, bool):
        return key if 1 <= key <= 6 else None
    if isinstance(key, str):
        if key not in {"1", "2", "3", "4", "5", "6"}:
            return None
        return int(key)
    return None


def _vector_component(vector: tuple[Any, ...], column: int) -> Any:
    if not 1 <= column <= 6:
        raise VolterraWordsInputError(
            "matrix column {} is outside the six-component state".format(column)
        )
    return vector[column - 1]


def _mapping_component(mapping: Mapping[Any, Any], index: int) -> Any:
    if index in mapping:
        return mapping[index]
    return mapping[str(index)]


def _sum_scalars(values: Any) -> Any:
    iterator = iter(values)
    try:
        total = next(iterator)
    except StopIteration:
        return 0
    for value in iterator:
        total = total + value
    return total


def _sum_vectors(vectors: Any, label: str) -> tuple[Any, ...]:
    iterator = iter(vectors)
    try:
        total = _green_vector(next(iterator), label)
    except StopIteration as exc:
        raise VolterraWordsInputError(label + " has no vectors") from exc
    for vector in iterator:
        current = _green_vector(vector, label)
        total = tuple(left + right for left, right in zip(total, current))
    return total


def _binomial(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    k = min(k, n - k)
    result = 1
    for factor in range(1, k + 1):
        result = result * (n - k + factor) // factor
    return result


__all__ = [
    "COMPONENT_MAPPING_ALIASES",
    "COMPONENT_ORDER",
    "FiniteNeumannEvaluation",
    "FiniteNeumannTerm",
    "FiniteNeumannWordPlan",
    "GREEN_SCOPE",
    "JET_COORDINATE",
    "SCHEMA",
    "SOURCE_ROWS",
    "VolterraWordsInputError",
    "compute_finite_neumann_words",
    "finite_neumann_words",
    "prepare_finite_neumann_words",
]
