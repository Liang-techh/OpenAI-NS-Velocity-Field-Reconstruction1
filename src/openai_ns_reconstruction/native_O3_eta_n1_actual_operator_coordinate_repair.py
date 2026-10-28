"""Coordinate-repaired v2 actual six-row n=1 coefficient operator.

This is a source-only draft.  Its sparse coefficient expressions are copied
from the paper-compliant ``raw_matrix_expressions`` interface used by the
preceding actual-operator draft.  The coordinate repair is deliberately
narrow: exact ``rho=x**2`` authentication is performed only when both values
have the strict built-in ``int`` or ``Fraction`` type.  Decimal, mpmath,
sympy, interval, and other approximate/enclosure values remain caller
preconditions; no tolerance or approximate comparison is promoted to a proof
of source-coordinate identity.  Safe non-negativity checks are retained.

The caller supplies every source value.  This module does not import or
execute ancestor source, Green terms, checkers, numerical evaluators, or any
m15/pressure construction.  Python ``float`` inputs remain forbidden, and
``rho=x**2`` plus ``L != 0`` remain caller-owned domain preconditions.
"""

from collections.abc import Callable, Mapping
from decimal import Decimal
from fractions import Fraction
from numbers import Real
from typing import Any


RAW_MATRIX_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
    ":raw_matrix_expressions"
)

W_COMPONENTS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
FEEDBACK_ROWS = ("4", "5", "6")
KINEMATIC_ENTRIES = {"1,5": 1, "2,6": 1, "3,6": -1}
_STRICT_EXACT_COORDINATE_TYPES = (int, Fraction)


def _reject_float_inputs(values: Mapping[str, Any]) -> None:
    """Reject Python float values before any source arithmetic is performed."""

    offenders = tuple(name for name, value in values.items() if isinstance(value, float))
    if offenders:
        names = ", ".join(offenders)
        raise TypeError("binary float inputs are forbidden: " + names)


def _safe_real(value: Any) -> bool:
    """Recognize scalar values for which a direct sign check is safe.

    Decimal is included explicitly because it is not registered as
    ``numbers.Real`` on all supported Python versions.  mpmath scalars are
    recognized by their type module without importing mpmath; their sign is
    checked only through the value's own comparison and never used for the
    coordinate identity proof.
    """

    if isinstance(value, bool):
        return False
    if isinstance(value, (Real, Decimal)):
        return True
    return type(value).__module__.startswith("mpmath")


def _safe_nonnegative_status(value: Any) -> bool | None:
    """Return a sign result only for a scalar whose comparison is boolean."""

    if not _safe_real(value):
        return None
    try:
        result = value >= 0
    except Exception:
        return None
    return result if isinstance(result, bool) else None


def _strict_exact_coordinate(value: Any) -> bool:
    """Require the exact built-in int/Fraction type, excluding subclasses."""

    return type(value) in _STRICT_EXACT_COORDINATE_TYPES


def _safe_exact_zero_status(value: Any) -> bool | None:
    """Return True/False only where zero equality is exact and safe.

    Intervals and unresolved symbolic expressions return None.  Sympy-like
    objects are accepted only when their explicit ``is_zero`` flag is a bool.
    mpmath values are checked through their own represented value for the
    separate ``L != 0`` guard; this check is not a coordinate identity proof.
    """

    if isinstance(value, bool):
        return value is False
    if isinstance(value, (int, Fraction, Decimal)):
        return value == 0
    is_zero = getattr(value, "is_zero", None)
    if isinstance(is_zero, bool):
        return is_zero
    if type(value).__module__.startswith("mpmath"):
        result = value == 0
        return result if isinstance(result, bool) else None
    return None


def validate_coordinate_precondition(*, x: Any, rho: Any) -> dict[str, Any]:
    """Check signs and authenticate ``rho=x**2`` only in the strict domain.

    Exact identity is checked only for values whose *runtime types* are
    exactly ``int`` or ``Fraction``.  For every other scalar/enclosure type,
    the returned status is ``caller_precondition_unverified``; no equality,
    tolerance, midpoint, or error-bound decision is made.
    """

    x_sign = _safe_nonnegative_status(x)
    rho_sign = _safe_nonnegative_status(rho)
    if x_sign is False:
        raise ValueError("x must be non-negative when supplied as a safe real value")
    if rho_sign is False:
        raise ValueError("rho must be non-negative when supplied as a safe real value")

    identity_checked = _strict_exact_coordinate(x) and _strict_exact_coordinate(rho)
    if identity_checked and rho != x * x:
        raise ValueError("caller precondition rho=x**2 is not exactly satisfied")

    return {
        "x_nonnegative_checked": x_sign is not None,
        "rho_nonnegative_checked": rho_sign is not None,
        "rho_equals_x_squared_checked": identity_checked,
        "rho_equals_x_squared_exact": identity_checked,
        "rho_equals_x_squared_status": (
            "exact_int_fraction" if identity_checked else "caller_precondition_unverified"
        ),
        "identity_proof_scope": "strict runtime int/Fraction only",
        "approximate_identity_used": False,
        "approximate_identity_rejected_as_proof": True,
    }


def _type_compatible_half(value: Any) -> Any:
    """Construct exact one-half without introducing a Python binary float."""

    if isinstance(value, float):
        raise TypeError("binary float input is forbidden for half construction")
    if isinstance(value, bool):
        raise TypeError("bool is not an accepted source coordinate type")
    if isinstance(value, int):
        return Fraction(1, 2)
    # For Fraction, Decimal, mpmath, sympy, and interval values, zero and one
    # are formed in the value's own algebra before division by integer 2.
    one = value * 0 + 1
    return one / 2


def build_actual_six_row_operator(
    *,
    F: Any,
    Fr: Any,
    FZ: Any,
    U: Any,
    Ur: Any,
    UZ: Any,
    Q: Any,
    epsilon: Any,
    delta: Any,
    rho: Any,
    Z: Any,
    x: Any,
    validate_coordinates: bool = True,
) -> dict[str, Any]:
    """Build source-equivalent sparse ``B0``, ``B1``, and kinematic-free C0.

    The caller owns every source value.  No midpoint or selected enclosure
    value is introduced.  ``B0["6,1"]`` retains the explicit epsilon square.
    """

    source_values = {
        "F": F,
        "Fr": Fr,
        "FZ": FZ,
        "U": U,
        "Ur": Ur,
        "UZ": UZ,
        "Q": Q,
        "epsilon": epsilon,
        "delta": delta,
        "rho": rho,
        "Z": Z,
        "x": x,
    }
    _reject_float_inputs(source_values)
    coordinate_check = (
        validate_coordinate_precondition(x=x, rho=rho)
        if validate_coordinates
        else {
            "x_nonnegative_checked": False,
            "rho_nonnegative_checked": False,
            "rho_equals_x_squared_checked": False,
            "rho_equals_x_squared_exact": False,
            "rho_equals_x_squared_status": "caller_precondition_unverified",
            "identity_proof_scope": "strict runtime int/Fraction only",
            "approximate_identity_used": False,
            "approximate_identity_rejected_as_proof": True,
        }
    )

    half = _type_compatible_half(Z)
    d = 1 - Z * Z
    L = 1 - delta * Z * Z
    L_zero_status = _safe_exact_zero_status(L)
    if L_zero_status is True:
        raise ZeroDivisionError("caller precondition L != 0 is violated exactly")

    H = Z * (1 - delta) * half + d * U
    B = Q + (1 - 2 * Z * U) / L
    E = Z * U - half
    Cf = F + rho * Fr
    Zf = ((-2 - delta) * Z * F + d * FZ - 2 * Z * rho * Fr) / L
    Zw = ((-1 - delta) * Z * U + d * UZ - 2 * Z * rho * Ur) / L

    B0 = {
        "1,5": 1,
        "2,6": 1,
        "3,6": -1,
        "4,1": epsilon * 4 * x * F,
        "5,1": epsilon * 2 * (Q + E * (-2 + delta) / L),
        "5,2": epsilon * (2 * Zf - 2 * (-1 + delta) * Z * Cf / L),
        "5,3": -epsilon * 2 * (1 + delta) * Z * Cf / L,
        "5,5": epsilon * x * B,
        # This is intentionally epsilon**2, not epsilon or a rescaled value.
        "6,1": -epsilon**2 * 8 * Z * rho * F / L,
        "6,2": epsilon
        * 2
        * (E * (-1 + delta) / L + Zw - Ur * Z * (-1 + delta) * rho / L),
        "6,3": -epsilon * 2 * (1 + delta) * Z * rho * Ur / L,
        "6,4": -epsilon * 4 * Z / L,
        "6,6": epsilon * x * B,
    }
    B1 = {
        "5,1": epsilon * 2 * H / L,
        "5,2": -epsilon * 2 * d * Cf / L,
        "5,3": -epsilon * 2 * d * Cf / L,
        "6,2": epsilon * 2 * (H - Ur * d * rho) / L,
        "6,3": -epsilon * 2 * d * rho * Ur / L,
        "6,4": epsilon * 2 * d / L,
    }
    C0 = {key: value for key, value in B0.items() if key not in KINEMATIC_ENTRIES}

    return {
        "B0": B0,
        "B1": B1,
        "C0": C0,
        "kinematic_entries_removed_from_C0": dict(KINEMATIC_ENTRIES),
        "feedback_rows": FEEDBACK_ROWS,
        "W_component_order": W_COMPONENTS,
        "coordinate_check": coordinate_check,
        "domain_conditions": {
            "rho_equals_x_squared": (
                "exact_int_fraction_only_or_caller_precondition_unverified"
            ),
            "L_nonzero": "caller_precondition",
            "L_exact_zero_check_status": L_zero_status,
            "zero_axis_has_no_x_inverse": True,
            "binary_float_inputs_rejected": True,
            "interval_midpoint_selected": False,
            "approximate_coordinate_identity_promoted": False,
        },
        "source_authentication": {
            "raw_matrix_source": RAW_MATRIX_SOURCE,
            "source_digest_verified": False,
            "source_receipt_verified": False,
            "source_ast_identity_verified": False,
            "source_authentication_missing": True,
        },
        "coefficient_operator_only": True,
        "does_not_construct": ("g", "P1", "m15"),
    }


def _vector_component(vector: Any, column: int) -> Any:
    if isinstance(vector, Mapping):
        if column in vector:
            return vector[column]
        return vector[str(column)]
    values = tuple(vector)
    if len(values) != 6:
        raise ValueError("W and W_Z must each have six components")
    return values[column - 1]


def _sparse_row_product(matrix: Mapping[str, Any], row: str, vector: Any) -> Any:
    terms = []
    for key, coefficient in matrix.items():
        row_key, column_key = key.split(",")
        if row_key == row:
            terms.append(coefficient * _vector_component(vector, int(column_key)))
    if not terms:
        return 0
    result = terms[0]
    for term in terms[1:]:
        result = result + term
    return result


def apply_six_row_feedback(operator: Mapping[str, Any], W: Any, W_Z: Any) -> dict[str, Any]:
    """Return rows 4--6 of ``C0*W + B1*W_Z`` only."""

    return {
        row: _sparse_row_product(operator["C0"], row, W)
        + _sparse_row_product(operator["B1"], row, W_Z)
        for row in FEEDBACK_ROWS
    }


def apply_six_row_feedback_parts(
    operator: Mapping[str, Any], W: Any, W_Z: Any
) -> dict[str, dict[str, Any]]:
    """Return ``C0*W``, ``B1*W_Z``, and their explicit rows 4--6 sum."""

    C0_part = {
        row: _sparse_row_product(operator["C0"], row, W) for row in FEEDBACK_ROWS
    }
    B1_part = {
        row: _sparse_row_product(operator["B1"], row, W_Z) for row in FEEDBACK_ROWS
    }
    return {
        "C0_W": C0_part,
        "B1_W_Z": B1_part,
        "feedback": {row: C0_part[row] + B1_part[row] for row in FEEDBACK_ROWS},
    }


def make_deferred_green_feedback(
    operator: Mapping[str, Any],
    *,
    W_callback: Callable[..., Any],
    W_Z_callback: Callable[..., Any],
) -> Callable[..., dict[str, Any]]:
    """Return an explicit deferred adapter for caller-owned Green six-vectors."""

    if not callable(W_callback) or not callable(W_Z_callback):
        raise TypeError("explicit W and W_Z callbacks are required")

    def evaluate(*args: Any, **kwargs: Any) -> dict[str, Any]:
        W = W_callback(*args, **kwargs)
        W_Z = W_Z_callback(*args, **kwargs)
        return apply_six_row_feedback(operator, W, W_Z)

    return evaluate


__all__ = [
    "FEEDBACK_ROWS",
    "KINEMATIC_ENTRIES",
    "RAW_MATRIX_SOURCE",
    "W_COMPONENTS",
    "apply_six_row_feedback",
    "apply_six_row_feedback_parts",
    "build_actual_six_row_operator",
    "make_deferred_green_feedback",
    "validate_coordinate_precondition",
]
