"""Bounded consumer for the saved N1 lower-collar cutoff domain.

This module consumes only the saved JSON mapping for the selected first-phase
endpoint and one injected interval evaluator.  Its scope is the direct
lower-collar ``chi(s)`` occurrence from the source definition, namely
``sigma(s)`` on ``s in [sc/4, sc/2]``.  It is *not* the normalized inner
cutoff ``sigma(4*s/sc_quarter - 2)`` and does not bind any normalized-inner
cutoff: the direct result is kept separate from that additional first-chart
record.  Neither result binds any normalized-inner node (including node
``52490``).  It does not construct the
selected ``sc`` as a Python integer or :class:`fractions.Fraction`: the saved
endpoint is an mpmath-style tuple ``[sign, mantissa, exponent, bitcount]`` and
its exponent is intentionally enormous in magnitude.  The only comparison
needed here is the exact bit-length implication

    mantissa < 2**bitcount,
    sc / 2 < 2**(exponent - 1 + bitcount) <= 2**(-128).

Consequently the previous rational interval leaf is called once on
``[Fraction(0), Fraction(1, 2**128)]`` for the direct occurrence and once on
``[Fraction(0), Fraction(1)]`` for the active first-chart normalized
argument.  The direct returned interval is a bound on the actual
lower-collar domain ``s in [sc/4, sc/2]``; neither upper endpoint is selected
as a function value.  The integral variable is represented as ``0 <= t <= s``
and is enclosed by the same safe dyadic domain.  The normalized result keeps
the flat windows ``q <= 1/2`` and ``q >= 3/4`` exact and separate.

Source files, source SHA-256 values, and any runner entry/receipt are outside
this consumer.  A saved mapping match is not a live source-owner, graph, or
physical identity proof.  In particular, this leaf makes no claim about
physical ``m15`` terms or five-moment closure.

The injected interval result must explicitly mark its derivative as an
ordinary-``x`` derivative (``ordinary_x_derivative=True``,
``z_derivative=False``), and must not present source metadata as an owner
proof.  Only its ``d_sigma_dx`` field is consumed; the direct result relabels
that bound as ``d_sigma_ds`` because its argument is the fixed collar ``s``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any, Callable


SCHEMA = "n1_collar_sigma_actual_domain_v1"
COMPONENT = "actual_lower_collar_sigma_interval_consumer"
RESULT_SCHEMA = SCHEMA + "_result"

SELECTED_ENDPOINT_PATH = (
    "explicit_current_inner_exit_strict_collar.selected_first_phase_endpoint"
)
SC_TUPLE_FIELDS = ("sign", "mantissa", "exponent", "bitcount")

# The evaluator input is deliberately small and fixed.  Constructing 2**128
# is harmless and is unrelated to the selected exponent; no power involving
# the selected exponent is ever constructed.
SAFE_ARGUMENT_LOWER = Fraction(0)
SAFE_ARGUMENT_UPPER = Fraction(1, 2**128)
SAFE_ARGUMENT_INTERVAL = (SAFE_ARGUMENT_LOWER, SAFE_ARGUMENT_UPPER)

# The first-chart normalized cutoff has q = s/(sc/4) in [0, 1].  Its active
# source argument is x = 4*q - 2, so x in [0, 1] exactly on q in
# [1/2, 3/4].  This second small interval is separate from the direct
# lower-collar sigma(s) call above.
NORMALIZED_ARGUMENT_LOWER = Fraction(0)
NORMALIZED_ARGUMENT_UPPER = Fraction(1)
NORMALIZED_ARGUMENT_INTERVAL = (
    NORMALIZED_ARGUMENT_LOWER,
    NORMALIZED_ARGUMENT_UPPER,
)

# These are representation/resource guards, not estimates of the source.
MAX_SC_MANTISSA_BITS = 4096
MAX_SC_EXPONENT_BITS = 256
MAX_SC_BITCOUNT = MAX_SC_MANTISSA_BITS


class ActualCollarSigmaError(ValueError):
    """Base class for this fail-closed consumer."""


class SavedMappingError(ActualCollarSigmaError):
    """The saved endpoint mapping is absent or malformed."""


class ScTupleError(SavedMappingError):
    """The selected lower/upper exact tuple is not a valid exact match."""


class IntervalLeafContractError(ActualCollarSigmaError):
    """The injected interval leaf did not return rational ordered bounds."""


# Compatibility spellings make the failure class discoverable to a runner
# without changing the deliberately narrow source-binding scope.
SourceBindingError = SavedMappingError
BindingError = SavedMappingError


def _factorial(value: int) -> int:
    result = 1
    for index in range(2, value + 1):
        result *= index
    return result


def _actual_derivative_cap() -> Fraction:
    """Return the endpoint cap on ``0 <= x <= 2**(-128)``.

    For ``0 < x <= 1/4``, write

        t = x**(-2) - (1-x)**(-2).

    The endpoint proof uses the elementary inequalities

        t >= (8/9) * x**(-2),
        2*x**(-3) + 2*(1-x)**(-3) <= (56/27) * x**(-3),
        exp(-t) <= 8! / t**8.

    Since ``sigma'(x) = q/(1+q)**2`` times the displayed derivative scale,
    with ``q = exp(-t)``, these imply

        sigma'(x)
        <= (8! / t**8) * (56/27) * x**(-3)
        <= 8! * (9/8)**8 * (56/27) * x**13
        = C * x**13.

    Monotonicity of ``x**13`` on this positive interval and
    ``x <= 2**(-128)`` give the cap below.  This is a conditional rational
    mathematical bound for the direct argument, not a live source-identity
    or owner certificate.

    Replacing ``x`` by ``2**(-128)`` gives this exact rational cap; it is a
    certified bound, not a value selected for ``sigma'``.
    """

    return (
        Fraction(_factorial(8), 1)
        * Fraction(9, 8) ** 8
        * Fraction(56, 27)
        * SAFE_ARGUMENT_UPPER**13
    )


ACTUAL_DERIVATIVE_CAP = _actual_derivative_cap()
DERIVATIVE_CAP = ACTUAL_DERIVATIVE_CAP


def _is_int(value: Any) -> bool:
    return type(value) is int


def validate_sc_tuple(raw: Any, label: str = "sc tuple") -> tuple[int, int, int, int]:
    """Validate and canonicalize one exact ``[sign, mantissa, exponent, bitcount]``.

    ``bitcount`` is checked against ``mantissa.bit_length()``.  This is the
    only large-integer operation performed on the selected endpoint.
    """

    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes, bytearray)):
        raise ScTupleError(f"{label} must be a four-entry sequence")
    if len(raw) != 4:
        raise ScTupleError(f"{label} must have four entries")
    sign, mantissa, exponent, bitcount = raw
    if not _is_int(sign) or sign not in (0, 1):
        raise ScTupleError(f"{label} sign must be 0 or 1")
    if not _is_int(mantissa) or mantissa <= 0:
        raise ScTupleError(f"{label} mantissa must be a positive integer")
    if mantissa.bit_length() > MAX_SC_MANTISSA_BITS:
        raise ScTupleError(f"{label} mantissa exceeds the bit budget")
    if not _is_int(exponent):
        raise ScTupleError(f"{label} exponent must be an integer")
    if abs(exponent).bit_length() > MAX_SC_EXPONENT_BITS:
        raise ScTupleError(f"{label} exponent exceeds the bit budget")
    if not _is_int(bitcount) or bitcount <= 0:
        raise ScTupleError(f"{label} bitcount must be a positive integer")
    if bitcount > MAX_SC_BITCOUNT:
        raise ScTupleError(f"{label} bitcount exceeds the bit budget")
    actual_bitcount = mantissa.bit_length()
    if actual_bitcount != bitcount:
        raise ScTupleError(
            f"{label} bitcount {bitcount} disagrees with mantissa bit length "
            f"{actual_bitcount}"
        )
    return (sign, mantissa, exponent, bitcount)


def _selected_sc_tuple(saved_mapping: Mapping[str, Any]) -> tuple[int, int, int, int]:
    if not isinstance(saved_mapping, Mapping):
        raise SavedMappingError("saved JSON endpoint must be a mapping")
    explicit = saved_mapping.get("explicit_current_inner_exit_strict_collar")
    if not isinstance(explicit, Mapping):
        raise SavedMappingError(
            "saved JSON mapping is missing explicit_current_inner_exit_strict_collar"
        )
    selected = explicit.get("selected_first_phase_endpoint")
    if not isinstance(selected, Mapping):
        raise SavedMappingError(
            "saved JSON mapping is missing selected_first_phase_endpoint"
        )
    if "lower_exact_mpf_tuple" not in selected:
        raise SavedMappingError("selected endpoint lower tuple is missing")
    if "upper_exact_mpf_tuple" not in selected:
        raise SavedMappingError("selected endpoint upper tuple is missing")
    lower = validate_sc_tuple(
        selected["lower_exact_mpf_tuple"], "selected lower_exact_mpf_tuple"
    )
    upper = validate_sc_tuple(
        selected["upper_exact_mpf_tuple"], "selected upper_exact_mpf_tuple"
    )
    if lower != upper:
        raise ScTupleError(
            "selected lower_exact_mpf_tuple and upper_exact_mpf_tuple must match "
            "exactly; an interval mismatch is not an exact sc"
        )
    return lower


def prove_half_sc_domain(raw: Any) -> dict[str, Any]:
    """Prove ``0 <= sc/2 <= 2**(-128)`` without materializing ``sc``.

    The returned mapping records the integer comparison used in the proof;
    it intentionally contains no huge denominator or selected numeric value.
    """

    sc_tuple = validate_sc_tuple(raw)
    sign, _mantissa, exponent, bitcount = sc_tuple
    if sign != 0:
        raise ScTupleError("selected sc must be nonnegative (sign must be 0)")
    half_exponent_bound = exponent - 1 + bitcount
    if half_exponent_bound > -128:
        raise ScTupleError(
            "selected sc/2 is not proved inside the 2**(-128) evaluator domain"
        )
    return {
        "sc_tuple": sc_tuple,
        "sign_zero": True,
        "mantissa_bit_length": bitcount,
        "mantissa_lt_two_to_bitcount": True,
        "half_exponent_bound": half_exponent_bound,
        "target_exponent": -128,
        "comparison": (
            "0 <= sc/2 < 2**(exponent-1+bitcount) <= 2**(-128); "
            "mantissa < 2**bitcount"
        ),
        "strict_bit_length_proof": True,
        "sc_value_materialized": False,
    }


def extract_selected_sc_tuple(
    saved_mapping: Mapping[str, Any],
) -> tuple[int, int, int, int]:
    """Return the exact selected tuple after lower/upper equality checking."""

    return _selected_sc_tuple(saved_mapping)


def _default_interval_evaluator() -> Callable[[Fraction, Fraction], Any]:
    """Resolve the previous interval leaf only when a caller did not inject it."""

    try:
        from openai_ns_reconstruction.n1_sigma_interval_rational_v1 import (
            evaluate_sigma_interval,
        )
    except ImportError:
        # A runner may load this file under a temporary module name while the
        # package itself is still importable only through a relative import.
        from .n1_sigma_interval_rational_v1 import evaluate_sigma_interval

    return evaluate_sigma_interval


def _field(value: Any, names: tuple[str, ...], label: str) -> Any:
    if isinstance(value, Mapping):
        for name in names:
            if name in value:
                return value[name]
        raise IntervalLeafContractError(f"{label} is missing")
    for name in names:
        if hasattr(value, name):
            return getattr(value, name)
    raise IntervalLeafContractError(f"{label} is missing")


def _rational_interval(value: Any, label: str) -> dict[str, Fraction]:
    if isinstance(value, Mapping):
        if "lower" in value and "upper" in value:
            lower, upper = value["lower"], value["upper"]
        else:
            nested = value.get("interval", value.get("enclosure"))
            if nested is None or nested is value:
                raise IntervalLeafContractError(f"{label} is not an interval")
            return _rational_interval(nested, label)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) != 2:
            raise IntervalLeafContractError(f"{label} sequence must have two entries")
        lower, upper = value
    elif hasattr(value, "lower") and hasattr(value, "upper"):
        lower, upper = value.lower, value.upper
    else:
        raise IntervalLeafContractError(f"{label} is not an interval")
    if not isinstance(lower, Fraction) or not isinstance(upper, Fraction):
        raise IntervalLeafContractError(f"{label} endpoints must be Fraction")
    if lower > upper:
        raise IntervalLeafContractError(f"{label} endpoints are unordered")
    return {"lower": lower, "upper": upper}


def _optional_input_interval(value: Any) -> dict[str, Fraction] | None:
    if isinstance(value, Mapping):
        candidate = value.get("input_interval")
    elif hasattr(value, "input_interval"):
        candidate = getattr(value, "input_interval")
    else:
        candidate = None
    if candidate is None:
        return None
    return _rational_interval(candidate, "interval leaf input_interval")


def _call_interval_leaf(
    evaluator: Callable[[Fraction, Fraction], Any],
    lower: Fraction,
    upper: Fraction,
    label: str,
) -> dict[str, dict[str, Fraction]]:
    """Make one checked call and extract only its rational interval fields."""

    leaf_result = evaluator(lower, upper)
    ordinary_x_derivative = _field(
        leaf_result,
        ("ordinary_x_derivative",),
        f"{label} interval leaf derivative coordinate",
    )
    z_derivative = _field(
        leaf_result,
        ("z_derivative",),
        f"{label} interval leaf Z derivative marker",
    )
    source_identity_verified = _field(
        leaf_result,
        ("source_identity_verified",),
        f"{label} interval leaf source-identity marker",
    )
    owner_binding = _field(
        leaf_result,
        ("owner_binding",),
        f"{label} interval leaf owner-binding marker",
    )
    if ordinary_x_derivative is not True:
        raise IntervalLeafContractError(
            f"{label} interval leaf derivative is not explicitly ordinary-x"
        )
    if z_derivative is not False:
        raise IntervalLeafContractError(
            f"{label} interval leaf derivative is marked as a Z derivative"
        )
    if source_identity_verified is not False or owner_binding is not False:
        raise IntervalLeafContractError(
            f"{label} interval leaf source metadata must not prove ownership"
        )
    input_interval = _optional_input_interval(leaf_result)
    expected_input = {"lower": lower, "upper": upper}
    if input_interval is not None and input_interval != expected_input:
        raise IntervalLeafContractError(
            f"{label} interval leaf returned a different input interval"
        )
    sigma = _rational_interval(
        _field(leaf_result, ("sigma",), f"{label} interval leaf sigma"),
        f"{label} interval leaf sigma",
    )
    derivative = _rational_interval(
        _field(
            leaf_result,
            ("d_sigma_dx",),
            f"{label} interval leaf ordinary derivative",
        ),
        f"{label} interval leaf ordinary derivative",
    )
    if sigma["lower"] < 0 or sigma["upper"] > 1:
        raise IntervalLeafContractError(
            f"{label} interval leaf sigma must lie in [0, 1]"
        )
    if derivative["lower"] < 0:
        raise IntervalLeafContractError(
            f"{label} interval leaf derivative must be nonnegative"
        )
    return {"sigma": sigma, "derivative": derivative}


def _resolve_fixed_s(
    fixed_s_in_z: bool | None,
    fixed_s: bool | None,
    fixed_s_in_Z: bool | None,
) -> bool:
    choices = [
        value for value in (fixed_s_in_z, fixed_s, fixed_s_in_Z) if value is not None
    ]
    if any(type(value) is not bool for value in choices):
        raise TypeError("fixed-s Z derivative conditions must be boolean")
    if len(set(choices)) > 1:
        raise ValueError("fixed-s Z derivative conditions disagree")
    return choices[0] if choices else False


def evaluate_actual_collar_sigma(
    saved_mapping: Mapping[str, Any],
    interval_evaluator: Callable[[Fraction, Fraction], Any] | None = None,
    *,
    fixed_s_in_z: bool | None = None,
    fixed_s: bool | None = None,
    fixed_s_in_Z: bool | None = None,
    fixed_s_over_sc_in_z: bool | None = None,
    fixed_s_over_sc: bool | None = None,
    fixed_s_over_sc_in_Z: bool | None = None,
) -> dict[str, Any]:
    """Consume the direct and normalized first-chart cutoff domains.

    ``interval_evaluator`` is normally the previous leaf's
    ``evaluate_sigma_interval``.  A runner may inject it explicitly; the
    small fixture does so with a synthetic callable and therefore does not
    re-enter the previous fixture.  The consumer makes exactly two calls:

    * direct lower collar ``sigma(s)`` on ``[0, 2**(-128)]``;
    * first-chart active argument ``sigma(x)`` on ``x in [0, 1]`` for
      ``x = 4*q - 2`` and ``q in [1/2, 3/4]``.

    The two returned records remain separate.  If omitted, the previous
    evaluator is imported lazily and called only at those two inputs.
    """

    sc_tuple = _selected_sc_tuple(saved_mapping)
    proof = prove_half_sc_domain(sc_tuple)
    fixed = _resolve_fixed_s(fixed_s_in_z, fixed_s, fixed_s_in_Z)
    normalized_fixed_choices = [
        value
        for value in (
            fixed_s_over_sc_in_z,
            fixed_s_over_sc,
            fixed_s_over_sc_in_Z,
        )
        if value is not None
    ]
    if any(type(value) is not bool for value in normalized_fixed_choices):
        raise TypeError("fixed-s/sc Z derivative conditions must be boolean")
    if len(set(normalized_fixed_choices)) > 1:
        raise ValueError("fixed-s/sc Z derivative conditions disagree")
    # A fixed direct ``s`` does not imply fixed ``s/sc``: the saved ``sc``
    # may itself vary with Z.  Normalized sigma_Z=0 therefore requires an
    # explicit fixed-ratio condition, even when fixed_s=True was supplied.
    normalized_fixed = normalized_fixed_choices[0] if normalized_fixed_choices else False
    evaluator = (
        _default_interval_evaluator()
        if interval_evaluator is None
        else interval_evaluator
    )
    if not callable(evaluator):
        raise TypeError("interval_evaluator must be callable")

    # The first call is the direct lower-collar occurrence sigma(s).
    direct_leaf = _call_interval_leaf(
        evaluator,
        SAFE_ARGUMENT_LOWER,
        SAFE_ARGUMENT_UPPER,
        "direct lower-collar",
    )
    sigma = direct_leaf["sigma"]
    direct_derivative_leaf = direct_leaf["derivative"]

    # Apply the independently derived narrow endpoint cap without entering
    # the previous leaf again for this direct occurrence.  An empty
    # intersection is a contract failure, never a silently selected endpoint.
    derivative_upper = min(
        direct_derivative_leaf["upper"], ACTUAL_DERIVATIVE_CAP
    )
    if direct_derivative_leaf["lower"] > derivative_upper:
        raise IntervalLeafContractError(
            "direct interval leaf derivative misses the actual rational cap"
        )
    direct_derivative = {
        "lower": direct_derivative_leaf["lower"],
        "upper": derivative_upper,
    }

    # The second call is the normalized first-chart active window.  Here q is
    # a positive order-preserving quotient s/(sc/4); no sc Fraction is built.
    normalized_leaf = _call_interval_leaf(
        evaluator,
        NORMALIZED_ARGUMENT_LOWER,
        NORMALIZED_ARGUMENT_UPPER,
        "normalized first-chart active window",
    )
    normalized_arg_sigma = normalized_leaf["sigma"]
    normalized_arg_derivative = normalized_leaf["derivative"]
    normalized_derivative = {
        "lower": Fraction(0),
        "upper": 4 * normalized_arg_derivative["upper"],
    }
    normalized_sigma = {
        # The strict flat source definition supplies the two plateau values;
        # the active window joins those endpoint values in the full q domain.
        "lower": Fraction(0),
        "upper": Fraction(1),
    }

    z_value: Fraction | None = Fraction(0) if fixed else None
    normalized_z_value: Fraction | None = (
        Fraction(0) if normalized_fixed else None
    )
    direct_result: dict[str, Any] = {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT + ".direct_lower_collar",
        "status": "bounded",
        "scope": "direct lower-collar chi sigma(s) only",
        "normalized_inner_cutoff_claimed": False,
        "normalized_inner_cutoff_node_binding": False,
        "sc_tuple": sc_tuple,
        "s_band": {
            "lower_factor": Fraction(1, 4),
            "upper_factor": Fraction(1, 2),
            "sc_tuple": sc_tuple,
            "representation": "[sc/4, sc/2]; sc is not materialized",
        },
        "sc_half_domain_proof": proof,
        "evaluated_argument_interval": {
            "lower": SAFE_ARGUMENT_LOWER,
            "upper": SAFE_ARGUMENT_UPPER,
        },
        "integral_variable": {
            "name": "t",
            "domain": "0 <= t <= s",
            "enclosure": {
                "lower": SAFE_ARGUMENT_LOWER,
                "upper": SAFE_ARGUMENT_UPPER,
            },
            "s_band_enclosed": True,
        },
        "sigma": sigma,
        "d_sigma_ds": direct_derivative,
        "d_sigma_dx": direct_derivative,
        "derivative": direct_derivative,
        "sigma_value": None,
        "d_sigma_ds_value": None,
        "sigma_upper_is_not_a_value": True,
        "derivative_upper_is_not_a_value": True,
        "derivative_cap": {
            "bound": ACTUAL_DERIVATIVE_CAP,
            "formula": "8!*(9/8)^8*(56/27)*(2^-128)^13",
            "is_rational_bound": True,
            "is_selected_value": False,
            "intersected_with_interval_leaf_upper": True,
            "conditional_math_bound_not_source_identity": True,
        },
        "interval_leaf_call_count": 1,
        "interval_leaf_input": (SAFE_ARGUMENT_LOWER, SAFE_ARGUMENT_UPPER),
        "fixed_s_in_z": fixed,
        "sigma_Z": z_value,
        "sigma_z": z_value,
        "fixed_s_over_sc_in_z": normalized_fixed,
        "normalized_sigma_Z": normalized_z_value,
        "normalized_sigma_z": normalized_z_value,
        "sigma_x_not_Z1": True,
        "z_derivative_scope": (
            "sigma_Z=0 only under the explicit ordinary-Z fixed-s condition"
            if fixed
            else "sigma_Z not asserted; dependence of s on Z is unverified"
        ),
        "ordinary_s_derivative_is_leaf_x_derivative": True,
    }
    normalized_result: dict[str, Any] = {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT + ".normalized_inner_cutoff",
        "status": "bounded",
        "scope": "first-chart normalized inner cutoff only",
        "normalized_inner_cutoff_claimed": True,
        "normalized_inner_cutoff_live_binding_claimed": False,
        "normalized_inner_cutoff_node_binding": False,
        "source_identity_verified": False,
        "owner_binding": False,
        "graph_binding": False,
        "sc_tuple": sc_tuple,
        "positive_sc_quarter": True,
        "q_definition": "q = s/(sc/4)",
        "q_domain": {"lower": Fraction(0), "upper": Fraction(1)},
        "order_preserving_division": True,
        "argument_definition": "x = 4*q - 2",
        "argument_domain": {"lower": Fraction(-2), "upper": Fraction(2)},
        "active_window": {
            "q_domain": {"lower": Fraction(1, 2), "upper": Fraction(3, 4)},
            "argument_domain": {
                "lower": NORMALIZED_ARGUMENT_LOWER,
                "upper": NORMALIZED_ARGUMENT_UPPER,
            },
            "leaf_sigma": normalized_arg_sigma,
            "leaf_d_sigma_dx": normalized_arg_derivative,
            "interval_leaf_call_count": 1,
            "interval_leaf_input": (
                NORMALIZED_ARGUMENT_LOWER,
                NORMALIZED_ARGUMENT_UPPER,
            ),
        },
        "flat_windows": {
            "q_le_half": {
                "q_domain": {"lower": Fraction(0), "upper": Fraction(1, 2)},
                "argument_domain": {"lower": Fraction(-2), "upper": Fraction(0)},
                "sigma": {"lower": Fraction(0), "upper": Fraction(0)},
                "d_sigma_dq": {"lower": Fraction(0), "upper": Fraction(0)},
                "source_definition": "sigma(x)=0 for x<=0",
            },
            "q_ge_three_quarters": {
                "q_domain": {"lower": Fraction(3, 4), "upper": Fraction(1)},
                "argument_domain": {"lower": Fraction(1), "upper": Fraction(2)},
                "sigma": {"lower": Fraction(1), "upper": Fraction(1)},
                "d_sigma_dq": {"lower": Fraction(0), "upper": Fraction(0)},
                "source_definition": "sigma(x)=1 for x>=1",
            },
        },
        "sigma": normalized_sigma,
        "d_sigma_dq": normalized_derivative,
        "sigma_value": None,
        "d_sigma_dq_value": None,
        "sigma_upper_is_not_a_value": True,
        "derivative_upper_is_not_a_value": True,
        "fixed_s_over_sc_in_z": normalized_fixed,
        "sigma_Z": normalized_z_value,
        "sigma_z": normalized_z_value,
        "z_derivative_scope": (
            "sigma_Z=0 only under fixed s/sc in ordinary Z"
            if normalized_fixed
            else "sigma_Z not asserted; dependence of s/sc on Z is unverified"
        ),
        "sigma_x_not_Z1": True,
        "argument_derivative_is_not_Z1": True,
    }
    result: dict[str, Any] = {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "status": "bounded",
        "saved_mapping_only": True,
        "cutoff_scope": "direct lower-collar chi sigma(s) only",
        "direct_lower_collar_sigma_only": True,
        "normalized_inner_cutoff_claimed": True,
        "normalized_inner_cutoff_live_binding_claimed": False,
        "normalized_inner_cutoff_node_binding": False,
        "selected_endpoint_path": SELECTED_ENDPOINT_PATH,
        "sc_tuple": sc_tuple,
        "selected_sc_tuple": sc_tuple,
        "sc_value_materialized": False,
        "sc_half_domain_proof": proof,
        "direct_lower_collar": direct_result,
        "normalized_inner_cutoff": normalized_result,
        "normalized_inner_cutoff_node_binding": False,
        "expected_interval_leaf_call_count": 2,
        "interval_leaf_call_count": 2,
        "sigma": sigma,
        "d_sigma_ds": direct_derivative,
        "d_sigma_dx": direct_derivative,
        "derivative": direct_derivative,
        "sigma_value": None,
        "d_sigma_ds_value": None,
        "sigma_upper_is_not_a_value": True,
        "derivative_upper_is_not_a_value": True,
        "derivative_cap": direct_result["derivative_cap"],
        "interval_leaf_input": (SAFE_ARGUMENT_LOWER, SAFE_ARGUMENT_UPPER),
        "interval_leaf_inputs": (
            (SAFE_ARGUMENT_LOWER, SAFE_ARGUMENT_UPPER),
            (NORMALIZED_ARGUMENT_LOWER, NORMALIZED_ARGUMENT_UPPER),
        ),
        "source_identity_verified": False,
        "owner_binding": False,
        "live_owner_binding": False,
        "graph_binding": False,
        "physical_m15": False,
        "five_moment_closure": False,
        "source_files_and_sha_authentication": "runner_only",
        "source_files_authenticated_here": False,
        "fixed_s_in_z": fixed,
        "sigma_Z": z_value,
        "sigma_z": z_value,
        "fixed_s_over_sc_in_z": normalized_fixed,
        "normalized_sigma_Z": normalized_z_value,
        "normalized_sigma_z": normalized_z_value,
        "sigma_x_not_Z1": True,
        "z_derivative_scope": (
            "sigma_Z=0 only under the explicit ordinary-Z fixed-s condition"
            if fixed
            else "sigma_Z not asserted; dependence of s on Z is unverified"
        ),
        "z_derivative_is_not_generalized": True,
        "certified": False,
        "rational_bound": True,
    }
    return result


# Concise aliases for callers that describe the operation as binding rather
# than evaluating.  They all preserve the same two-call contract.
bind_actual_collar_sigma = evaluate_actual_collar_sigma
evaluate_actual_domain = evaluate_actual_collar_sigma
consume = evaluate_actual_collar_sigma


__all__ = [
    "ACTUAL_DERIVATIVE_CAP",
    "ActualCollarSigmaError",
    "BindingError",
    "COMPONENT",
    "DERIVATIVE_CAP",
    "IntervalLeafContractError",
    "MAX_SC_BITCOUNT",
    "MAX_SC_EXPONENT_BITS",
    "MAX_SC_MANTISSA_BITS",
    "NORMALIZED_ARGUMENT_INTERVAL",
    "NORMALIZED_ARGUMENT_LOWER",
    "NORMALIZED_ARGUMENT_UPPER",
    "RESULT_SCHEMA",
    "SAFE_ARGUMENT_INTERVAL",
    "SAFE_ARGUMENT_LOWER",
    "SAFE_ARGUMENT_UPPER",
    "SCHEMA",
    "SC_TUPLE_FIELDS",
    "SELECTED_ENDPOINT_PATH",
    "SavedMappingError",
    "ScTupleError",
    "SourceBindingError",
    "bind_actual_collar_sigma",
    "consume",
    "evaluate_actual_collar_sigma",
    "evaluate_actual_domain",
    "extract_selected_sc_tuple",
    "prove_half_sc_domain",
    "validate_sc_tuple",
]
