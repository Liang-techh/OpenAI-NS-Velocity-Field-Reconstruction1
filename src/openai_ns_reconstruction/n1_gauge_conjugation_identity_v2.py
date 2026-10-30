"""Formal differential-operator conjugation for the n=1 gauge step.

This is a deliberately small algebra adapter.  It consumes caller-owned
matrices and does only the identities which are valid for a diagonal
``S(Z)`` (and, when supplied, a diagonal extra gauge ``T(Z)``):

``W = S*T*Y``

``B0* = S**-1*B0*S + S**-1*B1*S_Z`` and
``B1* = S**-1*B1*S``.

The implementation accepts ordinary Python matrices, SymPy matrices, and
other matrix-like values whose entries support ``+``, ``-`` and ``*``.  No
source producer, saved receipt, numerical owner, or historical builder is
called here.  In particular, algebraic conjugacy is not promoted to a
same-source identity: the returned packet retains caller ownership and
conditions and sets ``same_source_identity_proved`` to ``False``.

``A0`` and ``A1`` are the raw zero-order and ordinary-Z derivative
coefficients.  ``Akin`` and ``Achi`` are caller-owned principal pieces used
to form ``Qres = B0* - Akin - Achi``.  When the first and fifth diagonal
entries of ``S`` are provably equal, ``P = +2*chi*E51`` is also materialized;
both the residual without that compensation and the residual with it are
returned.  If the equality or ``chi`` value is not established, ``P`` is
zero and the packet records why it was not applied.

Only diagonal gauges are accepted.  This is intentional: an arbitrary
matrix gauge would require additional source/derivative ownership which this
leaf does not prove.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "n1_gauge_conjugation_identity_v2"
COMPONENT = "formal_differential_operator_conjugation_source_adapter"
ORDINARY_Z_ORDERS = (0, 1)
E51 = (4, 0)  # one-based E_51
# A single scalar ``G`` can be useful in a small synthetic algebra witness:
# components 1 and 5 are unchanged and the remaining components carry one
# inverse power.  This is *not* the natural-component gauge used by the real
# n=1 source, whose fourth component has an independent ``Pstar_squared``
# factor.  Real consumers must pass the complete six-entry ``T`` and ``T_Z``.
SYNTHETIC_EXTRA_WEIGHTS = (0, -1, -1, -1, 0, -1)
EXTRA_WEIGHTS = SYNTHETIC_EXTRA_WEIGHTS  # compatibility spelling; synthetic only
REAL_EXTRA_FORMULA = (
    "1",
    "Lambda^-1",
    "Lambda^-1",
    "Pstar_squared*Lambda^-1",
    "1",
    "Lambda^-1",
)
REAL_EXTRA_Z_CONSTANT = True


__all__ = [
    "SCHEMA",
    "COMPONENT",
    "ORDINARY_Z_ORDERS",
    "EXTRA_WEIGHTS",
    "SYNTHETIC_EXTRA_WEIGHTS",
    "REAL_EXTRA_FORMULA",
    "REAL_EXTRA_Z_CONSTANT",
    "E51",
    "conjugate_operator",
    "build",
    "adapt",
    "conjugate",
    "formal_conjugation",
]


def _unwrap(value: Any) -> tuple[Any, Mapping[str, Any]]:
    """Unwrap an optional caller packet without dropping its ownership."""

    if isinstance(value, Mapping) and "matrix" in value:
        metadata = value.get(
            "metadata", value.get("provenance", value.get("owner", {}))
        )
        if isinstance(metadata, Mapping):
            return value["matrix"], metadata
        return value["matrix"], {"owner": metadata}
    return value, {}


def _rows(value: Any, label: str) -> list[list[Any]]:
    """Convert a matrix-like value to rows while retaining entry objects."""

    value, _ = _unwrap(value)
    if hasattr(value, "tolist"):
        raw = value.tolist()
    elif hasattr(value, "rows") and hasattr(value, "cols"):
        raw = [[value[i, j] for j in range(int(value.cols))]
               for i in range(int(value.rows))]
    else:
        raw = value
    if isinstance(raw, (str, bytes, bytearray)) or not isinstance(raw, Sequence):
        raise TypeError(f"{label} must be a square matrix-like value")
    result = [list(row) for row in raw]
    if not result or any(isinstance(row, (str, bytes, bytearray))
                         or not isinstance(row, Sequence) for row in result):
        raise TypeError(f"{label} must be a non-empty square matrix")
    size = len(result)
    if any(len(row) != size for row in result):
        raise ValueError(f"{label} must be square")
    return result


def _is_zero(value: Any) -> bool:
    """Return True only when zero is established, never from a midpoint."""

    if value is None:
        return False
    if isinstance(value, bool):
        return value is False
    try:
        result = value == 0
        if isinstance(result, bool):
            return result
        # SymPy and a few symbolic backends return a BooleanAtom here.
        if getattr(result, "is_Boolean", False):
            return bool(result)
    except Exception:
        pass
    equals = getattr(value, "equals", None)
    if callable(equals):
        try:
            result = equals(0)
            if result is True:
                return True
            if result is False:
                return False
        except Exception:
            pass
    # mpmath interval endpoints and small exact interval objects.
    try:
        if hasattr(value, "_mpi_"):
            lo, hi = value._mpi_
            return lo == hi == 0
    except Exception:
        pass
    return False


def _prove_equal(left: Any, right: Any) -> bool:
    if left is right:
        return True
    try:
        return _is_zero(left - right)
    except Exception:
        return False


def _diagonal_entries(value: Any, size: int, label: str) -> list[Any]:
    """Validate a diagonal matrix/vector and return its diagonal entries."""

    value, _ = _unwrap(value)
    if _is_zero(value):
        return [0 for _ in range(size)]
    # A length-n sequence is a convenient diagonal-vector spelling.
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) == size and all(
            not isinstance(item, Sequence) or isinstance(item, (str, bytes, bytearray))
            for item in value
        ):
            return list(value)
    matrix = _rows(value, label)
    if len(matrix) != size:
        raise ValueError(f"{label} has dimension {len(matrix)}; expected {size}")
    for i, row in enumerate(matrix):
        for j, item in enumerate(row):
            if i != j and not _is_zero(item):
                raise ValueError(
                    f"{label} must be diagonal; off-diagonal entry ({i + 1},{j + 1}) "
                    "is not provably zero"
                )
    return [matrix[i][i] for i in range(size)]


def _zero_matrix(size: int) -> list[list[Any]]:
    return [[0 for _ in range(size)] for _ in range(size)]


def _identity_matrix(size: int) -> list[list[Any]]:
    return [[1 if i == j else 0 for j in range(size)] for i in range(size)]


def _inverse_diagonal(entries: Sequence[Any], label: str) -> list[list[Any]]:
    for index, item in enumerate(entries):
        if _is_zero(item):
            raise ValueError(f"{label}[{index + 1}] must be nonzero")
    return [[1 / entries[i] if i == j else 0 for j in range(len(entries))]
            for i in range(len(entries))]


def _matmul(left: Sequence[Sequence[Any]], right: Sequence[Sequence[Any]]) -> list[list[Any]]:
    n, inner, m = len(left), len(right), len(right[0])
    if len(left[0]) != inner:
        raise ValueError("matrix dimensions do not conform")
    output: list[list[Any]] = []
    for i in range(n):
        row: list[Any] = []
        for j in range(m):
            total: Any = 0
            for k in range(inner):
                total = total + left[i][k] * right[k][j]
            row.append(total)
        output.append(row)
    return output


def _matadd(left: Sequence[Sequence[Any]], right: Sequence[Sequence[Any]], sign: int = 1) -> list[list[Any]]:
    if len(left) != len(right) or any(len(a) != len(b) for a, b in zip(left, right)):
        raise ValueError("matrix dimensions do not conform")
    return [[a + sign * b for a, b in zip(row_a, row_b)]
            for row_a, row_b in zip(left, right)]


def _matzero(value: Sequence[Sequence[Any]]) -> bool:
    return all(_is_zero(item) for row in value for item in row)


def _restore_matrix(template: Any, rows: list[list[Any]]) -> Any:
    """Return a familiar backend type where that can be done safely."""

    raw, _ = _unwrap(template)
    if isinstance(raw, list):
        return [list(row) for row in rows]
    if isinstance(raw, tuple):
        return tuple(tuple(row) for row in rows)
    module = type(raw).__module__
    if module.startswith("sympy"):
        try:
            return type(raw)(rows)
        except Exception:
            pass
    if module.startswith("numpy"):
        try:
            import numpy as np  # optional dependency of the repository
            return np.asarray(rows, dtype=object)
        except Exception:
            pass
    try:
        return type(raw)(rows)
    except Exception:
        return [list(row) for row in rows]


def _as_vector(value: Any, size: int, label: str) -> tuple[list[Any], Any]:
    raw, _ = _unwrap(value)
    if hasattr(raw, "tolist"):
        raw = raw.tolist()
    if hasattr(raw, "rows") and hasattr(raw, "cols"):
        if int(raw.cols) == 1 and int(raw.rows) == size:
            values = [raw[i, 0] for i in range(size)]
            return values, value
        if int(raw.rows) == 1 and int(raw.cols) == size:
            values = [raw[0, i] for i in range(size)]
            return values, value
    if isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
        # numpy ``tolist`` may leave a one-row vector.
        if len(raw) == 1 and isinstance(raw[0], Sequence):
            raw = raw[0]
        values = list(raw)
        if len(values) == size:
            return values, value
    raise ValueError(f"{label} must contain exactly {size} entries")


def _restore_vector(template: Any, values: Sequence[Any]) -> Any:
    raw, _ = _unwrap(template)
    if isinstance(raw, tuple):
        return tuple(values)
    if isinstance(raw, list):
        return list(values)
    module = type(raw).__module__
    if module.startswith("sympy"):
        try:
            return type(raw)(list(values))
        except Exception:
            pass
    if module.startswith("numpy"):
        try:
            import numpy as np
            return np.asarray(list(values), dtype=object)
        except Exception:
            pass
    try:
        return type(raw)(values)
    except Exception:
        return list(values)


def _matvec(matrix: Sequence[Sequence[Any]], vector: Sequence[Any]) -> list[Any]:
    return [sum((matrix[i][j] * vector[j] for j in range(len(vector))), 0)
            for i in range(len(matrix))]


def _diagonal_matrix(entries: Sequence[Any]) -> list[list[Any]]:
    return [[entries[i] if i == j else 0 for j in range(len(entries))]
            for i in range(len(entries))]


def _is_matrix_like(value: Any) -> bool:
    if hasattr(value, "rows") and hasattr(value, "cols"):
        return True
    if hasattr(value, "tolist"):
        return True
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        return False
    return bool(value) and isinstance(value[0], Sequence)


def _is_diagonal_vector(value: Any, size: int) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes, bytearray))
        and len(value) == size
        and all(
            not isinstance(item, Sequence)
            or isinstance(item, (str, bytes, bytearray))
            for item in value
        )
    )


def _scalar_extra_entries(
    gauge: Any,
    gauge_z: Any,
    gauge_zz: Any,
    size: int,
    weights: Sequence[int],
) -> tuple[list[Any], list[Any], list[Any]]:
    """Build ``T=diag(G**w)`` and its ordinary Z derivatives.

    The scalar ``G`` is intentionally independent of ``S``.  For the fixed
    six-component weights, no unknown ``T_Z`` is needed: it is the ordinary
    derivative of the displayed powers of ``G``.
    """

    g0 = _jet(gauge, 0, "G")
    gz = _jet(gauge_z, 0, "G_Z", default=0)
    gzz = _jet(gauge_zz, 0, "G_ZZ", default=0)
    values: list[Any] = []
    derivatives: list[Any] = []
    second: list[Any] = []
    for weight in weights:
        if not isinstance(weight, int) or isinstance(weight, bool):
            raise ValueError("extra_weights must contain integers")
        if weight == 0:
            values.append(1)
            derivatives.append(0)
            second.append(0)
            continue
        value = g0 ** weight
        first = weight * (g0 ** (weight - 1)) * gz
        second_value = (
            weight * (weight - 1) * (g0 ** (weight - 2)) * gz * gz
            + weight * (g0 ** (weight - 1)) * gzz
        )
        values.append(value)
        derivatives.append(first)
        second.append(second_value)
    return values, derivatives, second


def _jet(value: Any, order: int, label: str, *, default: Any = None) -> Any:
    """Read scalar Z0/Z1 jets without interpreting a caller's source."""

    if value is None:
        return default
    if isinstance(value, Mapping):
        if order in value:
            return value[order]
        for key in (f"Z{order}", f"z{order}", "value" if order == 0 else "derivative"):
            if key in value:
                return value[key]
        if order == 1:
            return default
        raise ValueError(f"{label} is missing Z0")
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        # A two-entry scalar jet is distinct from a matrix/vector in this API.
        if len(value) >= order + 1:
            return value[order]
    return value if order == 0 else default


def _same_matrix(left: Sequence[Sequence[Any]], right: Sequence[Sequence[Any]]) -> bool:
    return len(left) == len(right) and all(
        len(a) == len(b) and all(_prove_equal(x, y) for x, y in zip(a, b))
        for a, b in zip(left, right)
    )


def _normalise_matrix(value: Any, label: str) -> tuple[list[list[Any]], Any, Mapping[str, Any]]:
    raw, metadata = _unwrap(value)
    matrix = _rows(raw, label)
    return matrix, raw, metadata


def conjugate_operator(
    A0: Any = None,
    A1: Any = None,
    S: Any = None,
    S_Z: Any = None,
    Akin: Any = None,
    Achi: Any = None,
    *,
    # ``B0``/``B1`` are accepted as explicit aliases for callers whose raw
    # operator is already named B.  They are not a second source.
    B0: Any = None,
    B1: Any = None,
    G: Any = None,
    G_Z: Any = None,
    extra_gauge: Any = None,
    extra_gauge_Z: Any = None,
    T: Any = None,
    T_Z: Any = None,
    chi: Any = None,
    chi_Z: Any = None,
    Y: Any = None,
    Y_Z: Any = None,
    state: Any = None,
    A0_Z: Any = None,
    A1_Z: Any = None,
    S_ZZ: Any = None,
    T_ZZ: Any = None,
    G_ZZ: Any = None,
    extra_weights: Any = None,
    extra_descriptor: Any = None,
    Akin_Z: Any = None,
    Achi_Z: Any = None,
    source_owner: Any = None,
    owner: Any = None,
    conditions: Any = None,
    **unexpected: Any,
) -> dict[str, Any]:
    """Build the forward and reverse formal conjugation packet.

    ``T``/``T_Z`` are the preferred explicit extra diagonal gauge.  The scalar
    ``G``/``G_Z`` spelling is retained only for a synthetic generic witness;
    it must not be labelled as the natural-component gauge because that gauge
    has the independent fourth entry ``Pstar_squared*Lambda^-1``.  Supplying
    both spellings is checked for agreement.  The raw names ``A0``/``A1`` and
    aliases ``B0``/``B1`` are likewise interchangeable.  Missing ``T_Z`` for
    an explicit ``T`` is treated only as a recorded constant-gauge default;
    real consumers should pass ``T_Z=0`` explicitly.  Missing ``S_Z`` is
    always an error because the physical connection may not be silently
    discarded.
    """

    if unexpected:
        names = ", ".join(sorted(str(name) for name in unexpected))
        raise TypeError(f"unexpected conjugation arguments: {names}")
    if G is None:
        G = extra_gauge
    elif extra_gauge is not None and extra_gauge is not G:
        raise TypeError("G and extra_gauge were supplied twice")
    if G_Z is None:
        G_Z = extra_gauge_Z
    elif extra_gauge_Z is not None and extra_gauge_Z is not G_Z:
        raise TypeError("G_Z and extra_gauge_Z were supplied twice")
    if state is not None:
        if not isinstance(state, Mapping) or "Y" not in state:
            raise TypeError("state must be a mapping with Y and optional Y_Z")
        if Y is not None and Y is not state["Y"]:
            raise TypeError("Y and state['Y'] were supplied twice")
        Y = state["Y"]
        if Y_Z is not None and state.get("Y_Z") is not None and Y_Z is not state["Y_Z"]:
            raise TypeError("Y_Z and state['Y_Z'] were supplied twice")
        Y_Z = state.get("Y_Z")
    if A0 is None:
        A0 = B0
    elif B0 is not None and B0 is not A0:
        raise TypeError("A0 and B0 were supplied twice")
    if A1 is None:
        A1 = B1
    elif B1 is not None and B1 is not A1:
        raise TypeError("A1 and B1 were supplied twice")
    if A0 is None or A1 is None or S is None:
        raise TypeError("raw A0, A1 and diagonal S are required")
    if S_Z is None:
        raise TypeError("caller-owned S_Z is required; zero may be passed explicitly")
    if source_owner is None:
        source_owner = owner
    if extra_descriptor is not None and not isinstance(extra_descriptor, Mapping):
        raise TypeError("extra_descriptor must be a mapping when supplied")
    if isinstance(extra_descriptor, Mapping) and extra_descriptor.get("kind") == "natural_component":
        formula = extra_descriptor.get("formula")
        if tuple(formula) != REAL_EXTRA_FORMULA:
            raise ValueError(
                "natural_component extra descriptor must retain the Pstar_squared fourth formula"
            )
        if T is None:
            raise ValueError(
                "natural_component extra gauge requires explicit six-entry T"
            )
        if T_Z is None:
            raise ValueError(
                "natural_component extra gauge requires explicit T_Z=0"
            )

    a0, a0_template, a0_meta = _normalise_matrix(A0, "A0")
    a1, a1_template, a1_meta = _normalise_matrix(A1, "A1")
    _, s_meta = _unwrap(S)
    _, sz_meta = _unwrap(S_Z)
    size = len(a0)
    if isinstance(extra_descriptor, Mapping) and extra_descriptor.get("kind") == "natural_component" and size != len(REAL_EXTRA_FORMULA):
        raise ValueError("natural_component extra gauge requires six components")
    if len(a1) != size:
        raise ValueError("A0 and A1 dimensions must agree")
    if Akin is None:
        akin = _zero_matrix(size)
        akin_template = a0_template
        akin_defaulted = True
    else:
        akin, akin_template, _ = _normalise_matrix(Akin, "Akin")
        akin_defaulted = False
    if Achi is None:
        achi = _zero_matrix(size)
        achi_template = a0_template
        achi_defaulted = True
    else:
        achi, achi_template, _ = _normalise_matrix(Achi, "Achi")
        achi_defaulted = False
    _, akin_meta = _unwrap(Akin) if Akin is not None else (None, {})
    _, achi_meta = _unwrap(Achi) if Achi is not None else (None, {})
    _, g_meta = _unwrap(G) if G is not None else (None, {})
    _, gz_meta = _unwrap(G_Z) if G_Z is not None else (None, {})
    for label, matrix in (("Akin", akin), ("Achi", achi)):
        if len(matrix) != size or any(len(row) != size for row in matrix):
            raise ValueError(f"{label} dimensions must agree with A0")

    s_entries = _diagonal_entries(S, size, "S")
    sz_entries = _diagonal_entries(S_Z, size, "S_Z")
    inv_s = _inverse_diagonal(s_entries, "S")
    s_matrix = _diagonal_matrix(s_entries)
    sz_matrix = _diagonal_matrix(sz_entries)
    inv_s_matrix = inv_s

    # Independent extra gauge: scalar G is lifted with the fixed component
    # weights; a diagonal G/T remains accepted for controlled callers that
    # already own the complete extra diagonal.  Both are independent of S.
    if G is not None and T is not None and G is not T:
        if (not _is_matrix_like(G) and not _is_diagonal_vector(G, size)):
            # Scalar G plus explicit T is checked against the weighted
            # construction below; it is a useful spelling for a controlled
            # consumer and is not a second independent source.
            pass
        elif (_is_matrix_like(G) or _is_matrix_like(T)
              or _is_diagonal_vector(G, size)
              or _is_diagonal_vector(T, size)):
            g_entries = _diagonal_entries(G, size, "G")
            t_entries = _diagonal_entries(T, size, "T")
            if not all(_prove_equal(a, b) for a, b in zip(g_entries, t_entries)):
                raise TypeError("G and T describe different extra gauges")
        else:
            raise TypeError("G and T were supplied twice")
    scalar_extra = (
        G is not None
        and not _is_matrix_like(G)
        and not _is_diagonal_vector(G, size)
    )
    if extra_weights is None:
        # The generic scalar path is synthetic only.  An explicit T has no
        # inferred weights, so its Pstar-squared fourth entry cannot be
        # accidentally relabelled as a scalar-G gauge.
        weights = (
            SYNTHETIC_EXTRA_WEIGHTS
            if scalar_extra and size == len(SYNTHETIC_EXTRA_WEIGHTS)
            else None
        )
    else:
        if (isinstance(extra_weights, (str, bytes, bytearray))
                or not isinstance(extra_weights, Sequence)
                or len(extra_weights) != size):
            raise ValueError("extra_weights must have one integer per component")
        weights = tuple(extra_weights)
    extra = T if T is not None else G
    extra_derivative = T_Z if T_Z is not None else G_Z
    extra_second = T_ZZ if T_ZZ is not None else G_ZZ
    extra_defaulted = extra is None
    extra_derivative_defaulted = False
    scalar_extra_derived = False
    if extra is None:
        t_entries = [1 for _ in range(size)]
        tz_entries = [0 for _ in range(size)]
        tzz_entries = [0 for _ in range(size)]
    elif scalar_extra:
        if weights is None:
            raise ValueError("scalar G requires six components or explicit extra_weights")
        t_entries, tz_entries, tzz_entries = _scalar_extra_entries(
            G, G_Z, G_ZZ, size, weights
        )
        scalar_extra_derived = True
        if T is not None:
            # The scalar path is the authoritative definition.  A caller
            # supplied T is accepted only when it exactly matches it.
            supplied = _diagonal_entries(T, size, "T")
            if not all(_prove_equal(a, b) for a, b in zip(supplied, t_entries)):
                raise TypeError("explicit T disagrees with weighted scalar G")
        if T_Z is not None:
            supplied_z = _diagonal_entries(T_Z, size, "T_Z")
            if not all(_prove_equal(a, b) for a, b in zip(supplied_z, tz_entries)):
                raise TypeError("explicit T_Z disagrees with weighted scalar G_Z")
    else:
        t_entries = _diagonal_entries(extra, size, "T/G")
        if extra_derivative is None:
            # Do not infer a source identity; only expose the constant-gauge
            # default and make it visible to a controlled consumer.
            tz_entries = [0 for _ in range(size)]
            extra_derivative_defaulted = True
        else:
            tz_entries = _diagonal_entries(extra_derivative, size, "T_Z/G_Z")
        if extra_second is None:
            tzz_entries = [0 for _ in range(size)]
        else:
            tzz_entries = _diagonal_entries(extra_second, size, "T_ZZ/G_ZZ")
    t_matrix = _diagonal_matrix(t_entries)
    tz_matrix = _diagonal_matrix(tz_entries)
    tzz_matrix = _diagonal_matrix(tzz_entries)
    if (isinstance(extra_descriptor, Mapping)
            and extra_descriptor.get("kind") == "natural_component"
            and any(not _is_zero(value) for value in tz_entries)):
        raise ValueError("natural_component EXTRA is Z-constant; pass T_Z=0")
    inv_t = _inverse_diagonal(t_entries, "T/G")

    # Raw B0/B1 are the caller-owned A0/A1 values.  Keep all matrix products
    # in entry form so exact symbolic objects are never coerced to floats.
    b0_star = _matadd(
        _matmul(_matmul(inv_s_matrix, a0), s_matrix),
        _matmul(_matmul(inv_s_matrix, a1), sz_matrix),
    )
    b1_star = _matmul(_matmul(inv_s_matrix, a1), s_matrix)
    c0_extra = _matadd(
        _matmul(_matmul(inv_t, b0_star), t_matrix),
        _matmul(_matmul(inv_t, b1_star), tz_matrix),
    )
    c1_extra = _matmul(_matmul(inv_t, b1_star), t_matrix)

    # Achi support is checked independently; P is not inferred from an
    # arbitrary matrix entry.  A supplied chi takes precedence and must agree
    # with the E51 coefficient when that coefficient is itself recognizable.
    nonzero_achi = [(i, j) for i, row in enumerate(achi) for j, value in enumerate(row)
                    if not _is_zero(value)]
    achi_e51_support = size > E51[0] and (
        not nonzero_achi or nonzero_achi == [E51]
    )
    chi_value = _jet(chi, 0, "chi") if chi is not None else None
    chi_derivative = _jet(chi_Z, 0, "chi_Z") if chi_Z is not None else None
    if chi_value is None and nonzero_achi == [E51]:
        chi_value = -achi[E51[0]][E51[1]] / 2
        chi_value_source = "derived from Achi[E51]"
    elif chi_value is not None:
        chi_value_source = "caller-owned chi"
    else:
        chi_value_source = "not supplied"
    if chi_derivative is None and isinstance(chi, Mapping):
        chi_derivative = _jet(chi, 1, "chi")
    if chi_derivative is None:
        chi_derivative = 0
        chi_derivative_defaulted = True
    else:
        chi_derivative_defaulted = False
    if chi_value is not None and nonzero_achi == [E51]:
        if not _prove_equal(achi[E51[0]][E51[1]], -2 * chi_value):
            raise ValueError("Achi[E51] and chi disagree")

    s_equal_1_5 = size > E51[0] and _prove_equal(s_entries[0], s_entries[4])
    p_applied = bool(s_equal_1_5 and chi_value is not None and achi_e51_support)
    p_rows = _zero_matrix(size)
    pz_rows = _zero_matrix(size)
    if p_applied:
        p_rows[E51[0]][E51[1]] = 2 * chi_value
        pz_rows[E51[0]][E51[1]] = 2 * chi_derivative
    p_extra = _matmul(_matmul(inv_t, p_rows), t_matrix)
    pz_extra = _matmul(_matmul(inv_t, pz_rows), t_matrix)

    qres = _matadd(_matadd(b0_star, akin, -1), achi, -1)
    qres_with_p = _matadd(qres, p_rows)
    akin_extra = _matmul(_matmul(inv_t, akin), t_matrix)
    achi_extra = _matmul(_matmul(inv_t, achi), t_matrix)
    qres_extra = _matadd(_matadd(c0_extra, akin_extra, -1), achi_extra, -1)
    qres_extra_with_p = _matadd(qres_extra, p_extra)
    t_akin_commutes = _same_matrix(_matmul(t_matrix, akin), _matmul(akin, t_matrix))
    t_achi_commutes = _same_matrix(_matmul(t_matrix, achi), _matmul(achi, t_matrix))

    # Reverse formulas are kept as explicit algebra, rather than merely
    # returning a boolean.  A controlled consumer can inspect every matrix.
    reverse_b0_star = _matadd(
        _matmul(_matmul(t_matrix, c0_extra), inv_t),
        _matmul(_matmul(t_matrix, c1_extra), _matmul(tz_matrix, inv_t)),
        -1,
    )
    reverse_b1_star = _matmul(_matmul(t_matrix, c1_extra), inv_t)
    reverse_a0 = _matadd(
        _matmul(_matmul(s_matrix, b0_star), inv_s_matrix),
        _matmul(_matmul(a1, sz_matrix), inv_s_matrix),
        -1,
    )
    reverse_a1 = _matmul(_matmul(s_matrix, b1_star), inv_s_matrix)
    reverse_a0_extra = _matadd(
        _matmul(_matmul(s_matrix, reverse_b0_star), inv_s_matrix),
        _matmul(_matmul(a1, sz_matrix), inv_s_matrix),
        -1,
    )
    reverse_checks = {
        "B0star_from_extra": _same_matrix(reverse_b0_star, b0_star),
        "B1star_from_extra": _same_matrix(reverse_b1_star, b1_star),
        "A0_from_S": _same_matrix(reverse_a0, a0),
        "A1_from_S": _same_matrix(reverse_a1, a1),
        "A0_from_ST": _same_matrix(reverse_a0_extra, a0),
    }

    output: dict[str, Any] = {
        "schema": SCHEMA,
        "component": COMPONENT,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "dimension": size,
        "A0": _restore_matrix(a0_template, a0),
        "A1": _restore_matrix(a1_template, a1),
        "B0": _restore_matrix(a0_template, a0),
        "B1": _restore_matrix(a1_template, a1),
        "S": _restore_matrix(S, s_matrix),
        "S_Z": _restore_matrix(S_Z, sz_matrix),
        "S_inverse": _restore_matrix(S, inv_s_matrix),
        "T": _restore_matrix(extra if extra is not None else a0_template, t_matrix),
        "T_Z": _restore_matrix(extra_derivative if extra_derivative is not None else a0_template, tz_matrix),
        "T_ZZ": _restore_matrix(extra_second if extra_second is not None else a0_template, tzz_matrix),
        "T_inverse": _restore_matrix(extra if extra is not None else a0_template, inv_t),
        "extra_weights": tuple(weights) if weights is not None else None,
        "extra_T_Z_derived_from_G": scalar_extra_derived,
        "extra_T_Z_formula": (
            "diag(w_i*G**(w_i-1)*G_Z)" if scalar_extra_derived
            else "caller-owned T_Z (or explicit constant-gauge default)"
        ),
        "extra_descriptor": extra_descriptor,
        "B0_star": _restore_matrix(a0_template, b0_star),
        "B1_star": _restore_matrix(a1_template, b1_star),
        "B0star": _restore_matrix(a0_template, b0_star),
        "B1star": _restore_matrix(a1_template, b1_star),
        "B0_extra": _restore_matrix(a0_template, c0_extra),
        "B1_extra": _restore_matrix(a1_template, c1_extra),
        "B0_star_with_extra_gauge": _restore_matrix(a0_template, c0_extra),
        "B1_star_with_extra_gauge": _restore_matrix(a1_template, c1_extra),
        "Akin": _restore_matrix(akin_template, akin),
        "Achi": _restore_matrix(achi_template, achi),
        "Akin_extra": _restore_matrix(a0_template, akin_extra),
        "Achi_extra": _restore_matrix(a0_template, achi_extra),
        "Qres": _restore_matrix(a0_template, qres),
        "Qres_star": _restore_matrix(a0_template, qres),
        "Qres_without_P": _restore_matrix(a0_template, qres),
        "Qres_with_P": _restore_matrix(a0_template, qres_with_p),
        "Qres_extra": _restore_matrix(a0_template, qres_extra),
        "Qres_total": _restore_matrix(a0_template, qres_extra),
        "Qres_extra_with_P": _restore_matrix(a0_template, qres_extra_with_p),
        "P": _restore_matrix(a0_template, p_rows),
        "P_Z": _restore_matrix(a0_template, pz_rows),
        "P_extra": _restore_matrix(a0_template, p_extra),
        "P_extra_Z": _restore_matrix(a0_template, pz_extra),
        "chi": chi_value,
        "chi_Z": chi_derivative,
        "reverse_identity": {
            "B0_star_from_extra": _restore_matrix(a0_template, reverse_b0_star),
            "B1_star_from_extra": _restore_matrix(a1_template, reverse_b1_star),
            "A0_from_S": _restore_matrix(a0_template, reverse_a0),
            "A1_from_S": _restore_matrix(a1_template, reverse_a1),
            "A0_from_ST": _restore_matrix(a0_template, reverse_a0_extra),
            "checks": reverse_checks,
            "passed": all(reverse_checks.values()),
            "formula": {
                "B0_star": "T*B0_extra*T^{-1} - B1_star*T_Z*T^{-1}",
                "B1_star": "T*B1_extra*T^{-1}",
                "A0": "S*B0_star*S^{-1} - A1*S_Z*S^{-1}",
                "A1": "S*B1_star*S^{-1}",
            },
        },
        "P_status": {
            "applied": p_applied,
            "formula": "+2*chi*E51",
            "index_zero_based": E51,
            "requires_S1_equals_S5": True,
            "S1_equals_S5_proved": s_equal_1_5,
            "Achi_E51_support_proved": achi_e51_support,
            "chi_supplied_or_derived": chi_value is not None,
            "chi_source": chi_value_source,
            "chi_Z_defaulted_zero": chi_derivative_defaulted,
            "not_applied_reason": None if p_applied else (
                "S1 != S5" if not s_equal_1_5 else
                "chi or E51 Achi support unavailable"
            ),
        },
        "gauge_scope": {
            "S_diagonal": True,
            "S_nonconstant_allowed": not all(_is_zero(item) for item in sz_entries),
            "S_Z_caller_owned": True,
            "extra_gauge_diagonal": True,
            "extra_scalar_gauge_weighted": scalar_extra_derived,
            "extra_scalar_gauge_scope": (
                "synthetic_generic_only" if scalar_extra_derived
                else "not_scalar_G"
            ),
            "extra_gauge_kind": (
                "caller_declared_natural_component"
                if isinstance(extra_descriptor, Mapping)
                and extra_descriptor.get("kind") == "natural_component"
                else "explicit_diagonal_or_identity"
            ),
            "real_natural_component_extra_requires_explicit_T": True,
            "real_natural_component_extra_formula": REAL_EXTRA_FORMULA,
            "real_natural_component_extra_power_pattern": (0, -1, -1, -1, 0, -1),
            "real_natural_component_fourth_multiplier": "Pstar_squared",
            "real_natural_component_extra_Z_constant": REAL_EXTRA_Z_CONSTANT,
            "real_natural_component_pstar_squared_slot": 3,
            "real_natural_component_formula_caller_declared": bool(
                isinstance(extra_descriptor, Mapping)
                and extra_descriptor.get("kind") == "natural_component"
            ),
            "real_natural_component_formula_verified": False,
            "extra_principal_commutes": {
                "Akin": t_akin_commutes,
                "Achi": t_achi_commutes,
            },
            "extra_gauge_derivative_defaulted": extra_derivative_defaulted,
            "arbitrary_matrix_gauge_supported": False,
        },
        "principal_inputs_defaulted": {
            "Akin": akin_defaulted,
            "Achi": achi_defaulted,
        },
        "ownership": {
            "source_owner": source_owner,
            "input_metadata": {
                "A0": a0_meta,
                "A1": a1_meta,
                "S": s_meta,
                "S_Z": sz_meta,
                "Akin": akin_meta,
                "Achi": achi_meta,
                "G": g_meta,
                "G_Z": gz_meta,
            },
            "conditions": conditions,
            "caller_conditions_retained": conditions is not None,
            "same_source_identity_proved": False,
            "same_source_identity_required": True,
            "caller_must_certify_source_identity": True,
            "source_executed": False,
            "historical_builder_called": False,
            "receipt_or_json_read": False,
        },
        "same_source_identity_proved": False,
        "caller_owned_algebra_only": True,
        "ordinary_Z1_leibniz": {
            "state_product_rule_available": Y is not None,
            "P_product_rule_available": p_applied,
            "formula": "(P*V)_Z=P_Z*V+P*V_Z",
        },
    }

    if Y is not None:
        y_values, y_template = _as_vector(Y, size, "Y")
        if Y_Z is None:
            y_z_values = [0 for _ in range(size)]
            y_z_defaulted = True
        else:
            y_z_values, _ = _as_vector(Y_Z, size, "Y_Z")
            y_z_defaulted = False
        w_values = _matvec(_matmul(s_matrix, t_matrix), y_values)
        wz_matrix = _matadd(
            _matadd(_matmul(sz_matrix, t_matrix), _matmul(s_matrix, tz_matrix)),
            _matmul(s_matrix, t_matrix),
        )
        wz_values = [
            value + derivative
            for value, derivative in zip(
                _matvec(wz_matrix, y_values),
                _matvec(_matmul(s_matrix, t_matrix), y_z_values),
            )
        ]
        p_state = _matvec(p_rows, y_values)
        p_state_z = [
            value + extra_value
            for value, extra_value in zip(
                _matvec(pz_rows, y_values), _matvec(p_rows, y_z_values)
            )
        ]
        output["state"] = {
            "Y": Y,
            "Y_Z": Y_Z,
            "W": _restore_vector(y_template, w_values),
            "W_Z": _restore_vector(y_template, wz_values),
            "P_Y": _restore_vector(y_template, p_state),
            "P_Y_Z": _restore_vector(y_template, p_state_z),
            "Y_Z_defaulted_zero": y_z_defaulted,
            "formula": {
                "W": "S*T*Y",
                "W_Z": "S_Z*T*Y + S*T_Z*Y + S*T*Y_Z",
                "P_Y_Z": "P_Z*Y + P*Y_Z",
            },
        }
        output["W"] = output["state"]["W"]
        output["W_Z"] = output["state"]["W_Z"]
        output["P_Y"] = output["state"]["P_Y"]
        output["P_Y_Z"] = output["state"]["P_Y_Z"]

    # Optional ordinary derivatives of the operator itself.  They are only
    # emitted when caller-owned derivative rows are present; no zero padding
    # is presented as a source identity.
    if A0_Z is not None or A1_Z is not None or S_ZZ is not None:
        if A0_Z is None or A1_Z is None or S_ZZ is None:
            raise TypeError("A0_Z, A1_Z and S_ZZ are all required for operator Z1")
        a0z, _, _ = _normalise_matrix(A0_Z, "A0_Z")
        a1z, _, _ = _normalise_matrix(A1_Z, "A1_Z")
        szz_entries = _diagonal_entries(S_ZZ, size, "S_ZZ")
        szz_matrix = _diagonal_matrix(szz_entries)
        inv_sz = _matmul(_matmul(inv_s_matrix, sz_matrix), inv_s_matrix)
        inv_sz = [[-value for value in row] for row in inv_sz]
        b0z = _matadd(
            _matadd(
                _matadd(_matmul(_matmul(inv_sz, a0), s_matrix),
                        _matmul(_matmul(inv_s_matrix, a0z), s_matrix)),
                _matmul(_matmul(inv_s_matrix, a0), sz_matrix)),
            _matadd(
                _matadd(_matmul(_matmul(inv_sz, a1), sz_matrix),
                        _matmul(_matmul(inv_s_matrix, a1z), sz_matrix)),
                _matmul(_matmul(inv_s_matrix, a1), szz_matrix)),
        )
        b1z = _matadd(
            _matadd(_matmul(_matmul(inv_sz, a1), s_matrix),
                    _matmul(_matmul(inv_s_matrix, a1z), s_matrix)),
            _matmul(_matmul(inv_s_matrix, a1), sz_matrix),
        )
        output["B0_star_Z"] = _restore_matrix(a0_template, b0z)
        output["B1_star_Z"] = _restore_matrix(a1_template, b1z)
        inv_tz = _matmul(_matmul(inv_t, tz_matrix), inv_t)
        inv_tz = [[-value for value in row] for row in inv_tz]
        c0z = _matadd(
            _matadd(
                _matadd(_matmul(_matmul(inv_tz, b0_star), t_matrix),
                        _matmul(_matmul(inv_t, b0z), t_matrix)),
                _matmul(_matmul(inv_t, b0_star), tz_matrix)),
            _matadd(
                _matadd(_matmul(_matmul(inv_tz, b1_star), tz_matrix),
                        _matmul(_matmul(inv_t, b1z), tz_matrix)),
                _matmul(_matmul(inv_t, b1_star), tzz_matrix)),
        )
        c1z = _matadd(
            _matadd(_matmul(_matmul(inv_tz, b1_star), t_matrix),
                    _matmul(_matmul(inv_t, b1z), t_matrix)),
            _matmul(_matmul(inv_t, b1_star), tz_matrix),
        )
        output["B0_extra_Z"] = _restore_matrix(a0_template, c0z)
        output["B1_extra_Z"] = _restore_matrix(a1_template, c1z)
        if Akin_Z is not None and Achi_Z is not None:
            akin_z, _, _ = _normalise_matrix(Akin_Z, "Akin_Z")
            achi_z, _, _ = _normalise_matrix(Achi_Z, "Achi_Z")
            qres_z = _matadd(_matadd(b0z, akin_z, -1), achi_z, -1)
            qres_z_p = _matadd(qres_z, pz_rows)
            qres_extra_z = _matadd(
                _matadd(c0z, _matadd(
                    _matadd(_matmul(_matmul(inv_tz, akin), t_matrix),
                            _matmul(_matmul(inv_t, akin_z), t_matrix)),
                    _matmul(_matmul(inv_t, akin), tz_matrix)), -1),
                _matadd(
                    _matadd(_matmul(_matmul(inv_tz, achi), t_matrix),
                            _matmul(_matmul(inv_t, achi_z), t_matrix)),
                    _matmul(_matmul(inv_t, achi), tz_matrix)),
                -1,
            )
            qres_extra_z_p = _matadd(
                qres_extra_z,
                _matadd(_matmul(_matmul(inv_tz, p_rows), t_matrix),
                        _matadd(_matmul(_matmul(inv_t, pz_rows), t_matrix),
                                _matmul(_matmul(inv_t, p_rows), tz_matrix))),
            )
            output["Qres_Z"] = _restore_matrix(a0_template, qres_z)
            output["Qres_Z_with_P"] = _restore_matrix(a0_template, qres_z_p)
            output["Qres_extra_Z"] = _restore_matrix(a0_template, qres_extra_z)
            output["Qres_extra_Z_with_P"] = _restore_matrix(a0_template, qres_extra_z_p)
            output["ordinary_Z1_leibniz"]["residual_Z1_available"] = True
        else:
            output["ordinary_Z1_leibniz"]["residual_Z1_available"] = False
        output["ordinary_Z1_leibniz"]["operator_Z1_available"] = True
    else:
        output["ordinary_Z1_leibniz"]["operator_Z1_available"] = False

    return output


build = conjugate_operator
adapt = conjugate_operator
conjugate = conjugate_operator
formal_conjugation = conjugate_operator
