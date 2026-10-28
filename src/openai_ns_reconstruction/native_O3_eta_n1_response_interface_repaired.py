"""Function-level interface for a same-source original ``n=1`` response.

This is an integration seam, not a source owner.  The caller supplies the
actual canonical response ``W`` and ordinary ``Z`` response ``W_Z`` together
with the leading field and a genuine radial integral operator.  No source
module, finite-prefix table, interval endpoint, sample value, or old solver is
called here.

The accepted raw canonical state is ordered as

``W = (F1, Uz1, K1, P1, F1x, Uz1x)``.

The physical ports used by the original moment formula are kept explicit:

``u = Uz1 / epsilon`` and ``p = P1 / (epsilon * Pstar**2)``.

After substituting those physical ports into the original temporal
first-order pressure moment and cancelling the exact ``Pstar**2`` factor, the
raw same-source form is

``m15^(1) = epsilon * (2 * integral(Ulead*Uz1)``
``                 + integral(P1) - rho*P1(rho))``.

The physical, unreduced form remains part of the binding metadata:
``epsilon**2 * (2 * integral(Ulead*u)``
``                 + Pstar**2 * (integral(p) - rho*p(rho)))``.

Its ordinary ``Z1`` response uses the supplied ``W_Z``, ``Ulead_Z`` and the
same integral operator.  The scale packet is required to declare epsilon and
Pstar ordinary-Z independent, as in the source definition; the code does not
finite-difference or silently differentiate a scale.

Axis conditions and error records are caller-owned evidence.  They are
returned unchanged and never promoted to a proof.  In particular, ``P1`` is a
pressure port, not an ``m15`` value, and this draft makes no claim that the
real source has been connected or that a global five-moment closure exists.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, TypeAlias


SCHEMA = "native_O3_eta_n1_response_interface_draft_v1"
TEMPORAL_ORDER = 1
M15_NAME = "m15"
RUNTIME_MODEL_METADATA = None
PHYSICAL_RHO_DOMAIN = (0, 4)
AXIS_RHO = 0
EXPLICIT_DEPENDENCIES = (
    "W",
    "W_Z",
    "Ulead",
    "Ulead_Z",
    "integrate",
    "epsilon",
    "Pstar",
)

RAW_STATE_ORDER = ("F1", "Uz1", "K1", "P1", "F1x", "Uz1x")
RAW_STATE_AXIS_ORDERS = {
    "F1": "O(rho)",
    "Uz1": "O(rho)",
    "K1": "O(rho)",
    "P1": "O(rho)",
    "F1x": "O(sqrt(rho))",
    "Uz1x": "O(sqrt(rho))",
}

M15_FORMULA = (
    "epsilon*(2*integral(Ulead*Uz1)"
    "+integral(P1)-rho*P1(rho))"
)
M15_ORDINARY_Z1_FORMULA = (
    "epsilon*(2*integral(Ulead_Z*Uz1+Ulead*Uz1_Z)"
    "+integral(P1_Z)-rho*P1_Z(rho))"
)
M15_PHYSICAL_FORMULA = (
    "epsilon^2*(2*integral(Ulead*u)"
    "+Pstar^2*(integral(p)-rho*p(rho)))"
)
M15_PHYSICAL_ORDINARY_Z1_FORMULA = (
    "epsilon^2*(2*integral(Ulead_Z*u+Ulead*u_Z)"
    "+Pstar^2*(integral(p_Z)-rho*p_Z(rho)))"
)
PHYSICAL_U_FORMULA = "u=Uz1/epsilon"
PHYSICAL_P_FORMULA = "p=P1/(epsilon*Pstar^2)"

FieldProvider: TypeAlias = Callable[[Any, Any], Any]
ComponentExtractor: TypeAlias = Callable[[Any], Any]
IntegralOperator: TypeAlias = Callable[
    [Callable[[Any], Any], Any, Any], Any
]
LogOperator: TypeAlias = Callable[[Any], Any]
ScaleProvider: TypeAlias = Any


class N1ResponseInterfaceError(ValueError):
    """Raised when required explicit response dependencies are missing."""


def _component(
    value: Any,
    names: Sequence[str],
    index: int,
    label: str,
) -> Any:
    """Read one raw state component without choosing an interval value."""

    if isinstance(value, Mapping):
        for name in names:
            if name in value:
                return value[name]
    if isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    ):
        try:
            return value[index]
        except IndexError as exc:
            raise N1ResponseInterfaceError(
                label + " sequence has no component " + str(index)
            ) from exc
    raise N1ResponseInterfaceError(
        label + " must expose a named or positional raw state component"
    )


def default_uz1_component(value: Any) -> Any:
    """Extract the raw ``Uz1`` component (canonical index 1)."""

    return _component(value, ("Uz1", "uz1", "W2"), 1, "Uz1")


def default_p1_component(value: Any) -> Any:
    """Extract the raw ``P1`` component (canonical index 3)."""

    return _component(value, ("P1", "p1", "W4"), 3, "P1")


def _resolve_scale(value: ScaleProvider, label: str) -> Any:
    """Resolve an exact scalar without silently accepting binary floats."""

    if callable(value):
        resolved = value()
    else:
        resolved = value
    if resolved is None:
        raise N1ResponseInterfaceError(label + " resolved to None")
    if isinstance(resolved, float):
        raise N1ResponseInterfaceError(
            label + " must be an exact/source scalar, not a binary float"
        )
    return resolved


@dataclass(frozen=True)
class SameSourceN1Response:
    """Explicit dependencies for one same-source regular ``n=1`` response.

    ``W`` and ``W_Z`` are raw canonical state providers.  They must come from
    the same source owner, but this interface cannot authenticate that fact;
    ``source_id`` and ``same_source_evidence`` are retained as caller evidence.
    ``integrate`` has the exact call shape ``integrate(f, 0, rho)`` and must be
    a genuine radial integral operator rather than a sample or quadrature
    shortcut selected by this module.
    """

    W: FieldProvider
    W_Z: FieldProvider
    Ulead: FieldProvider
    Ulead_Z: FieldProvider
    integrate: IntegralOperator
    epsilon: ScaleProvider
    Pstar: ScaleProvider
    scales_Z_independent: bool = False
    uz1_component: ComponentExtractor = default_uz1_component
    p1_component: ComponentExtractor = default_p1_component
    source_id: str | None = None
    same_source_evidence: Mapping[str, Any] = field(default_factory=dict)
    axis_conditions: Mapping[str, Any] = field(default_factory=dict)
    error_data: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class PhysicalN1Point:
    """Expose physical ports at one caller-requested ``(rho, Z)`` point."""

    rho: Any
    Z: Any
    raw_Uz1: Any
    raw_Uz1_Z: Any
    raw_P1: Any
    raw_P1_Z: Any
    epsilon: Any
    Pstar: Any
    Pstar_squared: Any
    u: Any
    u_Z: Any
    p: Any
    p_Z: Any
    scale_ratios: Mapping[str, str]
    formulas: Mapping[str, str]
    P1_is_not_m15: bool = True


@dataclass(frozen=True)
class RawN1Point:
    """Raw ``W/W_Z`` endpoint consumed by the reduced m15 functional."""

    rho: Any
    Z: Any
    raw_Uz1: Any
    raw_Uz1_Z: Any
    raw_P1: Any
    raw_P1_Z: Any
    epsilon: Any
    epsilon_ratio: str
    formulas: Mapping[str, str]
    source_callbacks_called: bool = True
    same_source_authenticated: bool = False
    P1_is_not_m15: bool = True


@dataclass(frozen=True)
class PrincipalGreenResponse:
    """The source-only principal Green response ``W0``.

    This is deliberately separate from the full canonical response ``W``:
    no ``B0`` or ``B1`` feedback is applied here.
    """

    x: Any
    Z: Any
    W0: tuple[Any, Any, Any, Any, Any, Any]
    source_rows: Mapping[str, str]
    metadata: Mapping[str, Any]


@dataclass(frozen=True)
class PrincipalGreenResponseWithZ1:
    """One principal Green response and its ordinary ``Z1`` companion."""

    value: PrincipalGreenResponse
    ordinary_Z1: PrincipalGreenResponse
    metadata: Mapping[str, Any]


@dataclass(frozen=True)
class N1M15Evaluation:
    """A formal affine physical ``m15`` evaluation and ordinary ``Z1`` value."""

    rho: Any
    Z: Any
    m15: Any
    ordinary_Z1: Any
    endpoint: RawN1Point
    leading_integral: Any
    leading_integral_Z1: Any
    pressure_integral: Any
    pressure_integral_Z1: Any
    pressure_boundary: Any
    pressure_boundary_Z1: Any
    retained_error_data: Mapping[str, Any]
    axis_conditions: Mapping[str, Any]
    metadata: Mapping[str, Any]

    @property
    def m15_Z1(self) -> Any:
        """Alias using the explicit moment name for ordinary ``Z1``."""

        return self.ordinary_Z1


def bind_same_source_n1_response(
    *,
    W: FieldProvider,
    W_Z: FieldProvider,
    Ulead: FieldProvider,
    Ulead_Z: FieldProvider,
    integrate: IntegralOperator,
    epsilon: ScaleProvider,
    Pstar: ScaleProvider,
    scales_Z_independent: bool,
    uz1_component: ComponentExtractor = default_uz1_component,
    p1_component: ComponentExtractor = default_p1_component,
    source_id: str | None = None,
    same_source_evidence: Mapping[str, Any] | None = None,
    axis_conditions: Mapping[str, Any] | None = None,
    error_data: Mapping[str, Any] | None = None,
) -> SameSourceN1Response:
    """Create a response binding from explicit caller-owned dependencies.

    This function only checks that callable dependencies are callable.  It
    does not execute ``W``, ``W_Z``, ``Ulead``, ``Ulead_Z``, the integral
    operator, or any source builder, and it does not inspect or reduce errors.
    """

    for name, value in (
        ("W", W),
        ("W_Z", W_Z),
        ("Ulead", Ulead),
        ("Ulead_Z", Ulead_Z),
        ("integrate", integrate),
        ("uz1_component", uz1_component),
        ("p1_component", p1_component),
    ):
        if not callable(value):
            raise N1ResponseInterfaceError(name + " must be callable")
    return SameSourceN1Response(
        W=W,
        W_Z=W_Z,
        Ulead=Ulead,
        Ulead_Z=Ulead_Z,
        integrate=integrate,
        epsilon=epsilon,
        Pstar=Pstar,
        scales_Z_independent=scales_Z_independent,
        uz1_component=uz1_component,
        p1_component=p1_component,
        source_id=source_id,
        same_source_evidence=(
            same_source_evidence if same_source_evidence is not None else {}
        ),
        axis_conditions=(
            axis_conditions if axis_conditions is not None else {}
        ),
        error_data=error_data if error_data is not None else {},
    )


def physical_n1_response_at(
    response: SameSourceN1Response,
    rho: Any,
    Z: Any,
) -> PhysicalN1Point:
    """Materialize only the caller-requested physical ``u,p`` ports.

    The raw source state is rescaled by the exact field diagonal from the
    canonical definition.  Ordinary ``Z`` derivatives are read from ``W_Z``;
    no numerical differentiation is attempted.  This function deliberately
    does not evaluate an integral or call an interval endpoint selector.
    """

    if response.scales_Z_independent is not True:
        raise N1ResponseInterfaceError(
            "ordinary Z1 requires scales_Z_independent=True; "
            "epsilon and Pstar derivatives are not inferred"
        )
    epsilon = _resolve_scale(response.epsilon, "epsilon")
    pstar = _resolve_scale(response.Pstar, "Pstar")
    pstar_squared = pstar * pstar
    w = response.W(rho, Z)
    w_z = response.W_Z(rho, Z)
    raw_uz1 = response.uz1_component(w)
    raw_uz1_z = response.uz1_component(w_z)
    raw_p1 = response.p1_component(w)
    raw_p1_z = response.p1_component(w_z)
    u_ratio = 1 / epsilon
    p_ratio = 1 / (epsilon * pstar_squared)
    return PhysicalN1Point(
        rho=rho,
        Z=Z,
        raw_Uz1=raw_uz1,
        raw_Uz1_Z=raw_uz1_z,
        raw_P1=raw_p1,
        raw_P1_Z=raw_p1_z,
        epsilon=epsilon,
        Pstar=pstar,
        Pstar_squared=pstar_squared,
        u=raw_uz1 * u_ratio,
        u_Z=raw_uz1_z * u_ratio,
        p=raw_p1 * p_ratio,
        p_Z=raw_p1_z * p_ratio,
        scale_ratios={
            "u_from_Uz1": "1/epsilon",
            "p_from_P1": "1/(epsilon*Pstar^2)",
        },
        formulas={
            "u": PHYSICAL_U_FORMULA,
            "p": PHYSICAL_P_FORMULA,
            "ordinary_Z_scale_rule": "D_Z(epsilon)=D_Z(Pstar)=0",
        },
    )


def _ordinary_nonnegative_scalar(x: Any, label: str) -> bool:
    """Require an ordinary scalar domain without interval truth coercion."""

    if isinstance(x, (str, bytes, bytearray, Mapping, Sequence)):
        raise N1ResponseInterfaceError(label + " must be an ordinary scalar")
    try:
        is_zero = x == 0
        is_negative = x < 0
    except Exception as exc:
        raise N1ResponseInterfaceError(
            label + " must support ordinary scalar comparisons"
        ) from exc
    if not isinstance(is_zero, bool) or not isinstance(is_negative, bool):
        raise N1ResponseInterfaceError(
            label + " interval/array truth is delegated to the caller"
        )
    if is_negative:
        raise N1ResponseInterfaceError(label + " must satisfy " + label + " >= 0")
    return is_zero


def _continuous_log_kernel(
    x: Any,
    s: Any,
    g6_at: Callable[[Any], Any],
    log: LogOperator,
) -> Any:
    """Use the continuous value zero at ``s=0`` for ``s*log(x/s)*g6``."""

    is_zero = _ordinary_nonnegative_scalar(s, "Green integration coordinate s")
    if is_zero:
        return 0
    return s * log(x / s) * g6_at(s)


def principal_green_response(
    *,
    x: Any,
    Z: Any,
    g4: FieldProvider,
    g5: FieldProvider,
    g6: FieldProvider,
    integrate: IntegralOperator,
    log: LogOperator,
) -> PrincipalGreenResponse:
    """Return the source-only six-vector ``W0`` from a true integral owner.

    The caller supplies continuous source rows ``g4``, ``g5`` and ``g6`` as
    ``g_i(s, Z)``.  The exact principal kernels are

    ``W0_1 = integral(s/2*(1-(s/x)^2)*g5)``
    ``W0_2 = integral(s*log(x/s)*g6)``
    ``W0_3 = -integral(s/2*(1-(s/x)^2)*g6)``
    ``W0_4 = integral(g4)``
    ``W0_5 = integral((s/x)^3*g5)``
    ``W0_6 = integral((s/x)*g6)``.

    At ``x=0`` the regular source row is defined as the exact zero six-vector;
    no source callback or integral callback is invoked.  This function never
    applies the ``B0``/``B1`` feedback needed for the complete ``W``.
    """

    is_axis = _ordinary_nonnegative_scalar(x, "Green coordinate x")
    for name, value in (
        ("g4", g4),
        ("g5", g5),
        ("g6", g6),
        ("integrate", integrate),
        ("log", log),
    ):
        if not callable(value):
            raise N1ResponseInterfaceError(name + " must be callable")
    if is_axis:
        zeros = (0, 0, 0, 0, 0, 0)
        return PrincipalGreenResponse(
            x=x,
            Z=Z,
            W0=zeros,
            source_rows={"g4": "g4", "g5": "g5", "g6": "g6"},
            metadata={
                "schema": SCHEMA,
                "axis_x_zero": True,
                "axis_W0_exact_zero": True,
                "feedback_applied": False,
                "full_W": False,
                "preconditioned_green": True,
                "preconditioner": "(I+J)G is represented by this W0 kernel",
                "preconditioner_must_not_be_reapplied": True,
                "source_coordinate": "x (not rho); rho=x^2 only at caller boundary",
                "ordinary_scalar_domain_only": True,
                "interval_source_handled_by_caller_integrator": True,
                "source_connected_by_this_module": False,
            },
        )

    def g4_at(s: Any) -> Any:
        return g4(s, Z)

    def g5_at(s: Any) -> Any:
        return g5(s, Z)

    def g6_at(s: Any) -> Any:
        return g6(s, Z)

    w0_1 = integrate(
        lambda s: s / 2 * (1 - (s / x) ** 2) * g5_at(s),
        AXIS_RHO,
        x,
    )
    w0_2 = integrate(
        lambda s: _continuous_log_kernel(x, s, g6_at, log),
        AXIS_RHO,
        x,
    )
    w0_3 = -integrate(
        lambda s: s / 2 * (1 - (s / x) ** 2) * g6_at(s),
        AXIS_RHO,
        x,
    )
    w0_4 = integrate(g4_at, AXIS_RHO, x)
    w0_5 = integrate(
        lambda s: (s / x) ** 3 * g5_at(s),
        AXIS_RHO,
        x,
    )
    w0_6 = integrate(
        lambda s: (s / x) * g6_at(s),
        AXIS_RHO,
        x,
    )
    return PrincipalGreenResponse(
        x=x,
        Z=Z,
        W0=(w0_1, w0_2, w0_3, w0_4, w0_5, w0_6),
        source_rows={"g4": "g4", "g5": "g5", "g6": "g6"},
        metadata={
            "schema": SCHEMA,
            "axis_x_zero": False,
            "axis_W0_exact_zero": False,
            "feedback_applied": False,
            "full_W": False,
            "W0_definition": "principal Green source response only",
            "preconditioned_green": True,
            "preconditioner": "(I+J)G is represented by this W0 kernel",
            "preconditioner_must_not_be_reapplied": True,
            "source_coordinate": "x (not rho); rho=x^2 only at caller boundary",
            "ordinary_scalar_domain_only": True,
            "interval_source_handled_by_caller_integrator": True,
            "s_zero_log_kernel_continuous_extension": "0",
            "B0_B1_feedback_required_for_full_W": True,
            "source_connected_by_this_module": False,
            "interval_or_sample_selection": False,
        },
    )


def principal_green_response_with_Z1(
    *,
    x: Any,
    Z: Any,
    g4: FieldProvider,
    g5: FieldProvider,
    g6: FieldProvider,
    g4_Z: FieldProvider,
    g5_Z: FieldProvider,
    g6_Z: FieldProvider,
    integrate: IntegralOperator,
    log: LogOperator,
) -> PrincipalGreenResponseWithZ1:
    """Apply the same finite principal kernel to source rows and ``g_Z``.

    Since the ordinary ``Z`` derivative does not differentiate the radial
    kernel, the companion response is obtained by the same six formulas with
    ``g4_Z``, ``g5_Z`` and ``g6_Z``.  This is only a finite principal operator;
    it does not add ``B0``/``B1`` feedback or assert canonical-source identity.
    """

    value = principal_green_response(
        x=x,
        Z=Z,
        g4=g4,
        g5=g5,
        g6=g6,
        integrate=integrate,
        log=log,
    )
    ordinary_z1 = principal_green_response(
        x=x,
        Z=Z,
        g4=g4_Z,
        g5=g5_Z,
        g6=g6_Z,
        integrate=integrate,
        log=log,
    )
    return PrincipalGreenResponseWithZ1(
        value=value,
        ordinary_Z1=ordinary_z1,
        metadata={
            "schema": SCHEMA,
            "ordinary_Z1": True,
            "same_radial_kernel": True,
            "kernel_Z_derivative": "zero; only g4_Z,g5_Z,g6_Z are substituted",
            "finite_principal_operator_only": True,
            "feedback_applied": False,
            "preconditioned_green": True,
            "preconditioner_must_not_be_reapplied": True,
            "source_connected_by_this_module": False,
            "source_authentication": False,
            "computed_from_caller_callbacks": True,
            "result_is_certified": False,
        },
    )


def _axis_contract_status(
    conditions: Mapping[str, Any],
) -> str:
    """Classify caller evidence without evaluating a field at the axis."""

    required = (
        "W_axis_zero",
        "Uz1_axis_zero",
        "P1_axis_normalized",
    )
    values = [conditions.get(name) for name in required]
    if all(value is True for value in values):
        return "passed_by_caller_evidence"
    if any(value is False for value in values):
        return "failed_by_caller_evidence"
    return "unverified"


def evaluate_n1_m15(
    response: SameSourceN1Response,
    rho: Any,
    Z: Any,
) -> N1M15Evaluation:
    """Evaluate raw same-source ``m15`` and ordinary ``Z1`` from true ports.

    The supplied integral operator receives the actual integrands and exact
    bounds ``0`` and ``rho``.  It is responsible for performing the true
    source-coordinate integral.  No midpoint, interval endpoint, norm cap,
    or finite sample is selected here.  The raw formula is exactly the
    physical formula after substituting ``u`` and ``p`` and cancelling
    ``Pstar**2``.
    """

    if response.scales_Z_independent is not True:
        raise N1ResponseInterfaceError(
            "ordinary Z1 requires scales_Z_independent=True; "
            "epsilon derivatives are not inferred"
        )
    # The reduced same-source formula only needs epsilon.  In particular,
    # do not force Pstar or the physical p conversion when evaluating m15.
    epsilon = _resolve_scale(response.epsilon, "epsilon")
    w = response.W(rho, Z)
    w_z = response.W_Z(rho, Z)
    raw_uz1 = response.uz1_component(w)
    raw_uz1_z = response.uz1_component(w_z)
    raw_p1 = response.p1_component(w)
    raw_p1_z = response.p1_component(w_z)
    endpoint = RawN1Point(
        rho=rho,
        Z=Z,
        raw_Uz1=raw_uz1,
        raw_Uz1_Z=raw_uz1_z,
        raw_P1=raw_p1,
        raw_P1_Z=raw_p1_z,
        epsilon=epsilon,
        epsilon_ratio="epsilon",
        formulas={
            "m15": M15_FORMULA,
            "ordinary_Z1": M15_ORDINARY_Z1_FORMULA,
        },
    )

    def uz1_at(radial: Any) -> Any:
        return response.uz1_component(response.W(radial, Z))

    def uz1_Z_at(radial: Any) -> Any:
        return response.uz1_component(response.W_Z(radial, Z))

    def p1_at(radial: Any) -> Any:
        return response.p1_component(response.W(radial, Z))

    def p1_Z_at(radial: Any) -> Any:
        return response.p1_component(response.W_Z(radial, Z))

    leading_integral = response.integrate(
        lambda radial: response.Ulead(radial, Z) * uz1_at(radial),
        AXIS_RHO,
        rho,
    )
    leading_integral_Z1 = response.integrate(
        lambda radial: (
            response.Ulead_Z(radial, Z) * uz1_at(radial)
            + response.Ulead(radial, Z) * uz1_Z_at(radial)
        ),
        AXIS_RHO,
        rho,
    )
    pressure_integral = response.integrate(p1_at, AXIS_RHO, rho)
    pressure_integral_Z1 = response.integrate(p1_Z_at, AXIS_RHO, rho)
    pressure_boundary = rho * endpoint.raw_P1
    pressure_boundary_Z1 = rho * endpoint.raw_P1_Z
    m15 = epsilon * (
        2 * leading_integral
        + pressure_integral
        - pressure_boundary
    )
    ordinary_z1 = epsilon * (
        2 * leading_integral_Z1
        + pressure_integral_Z1
        - pressure_boundary_Z1
    )
    axis_status = _axis_contract_status(response.axis_conditions)
    metadata = {
        "schema": SCHEMA,
        "runtime_model_metadata": RUNTIME_MODEL_METADATA,
        "temporal_order": TEMPORAL_ORDER,
        "ordinary_Z": True,
        "ordinary_Z1_is_supplied_product_rule": True,
        "formula": M15_FORMULA,
        "ordinary_Z1_formula": M15_ORDINARY_Z1_FORMULA,
        "physical_formula_before_exact_reduction": M15_PHYSICAL_FORMULA,
        "physical_ordinary_Z1_formula_before_exact_reduction": (
            M15_PHYSICAL_ORDINARY_Z1_FORMULA
        ),
        "exact_Pstar_squared_cancellation": (
            "Pstar^2*(P1/(epsilon*Pstar^2))=P1/epsilon"
        ),
        "endpoint_type": "RawN1Point",
        "Pstar_consumed_by_reduced_m15": False,
        "physical_conversion_is_optional": True,
        "physical_u_formula": PHYSICAL_U_FORMULA,
        "physical_p_formula": PHYSICAL_P_FORMULA,
        "explicit_dependencies": EXPLICIT_DEPENDENCIES,
        "source_id": response.source_id,
        "same_source_evidence": response.same_source_evidence,
        "epsilon_Pstar_ratios_exact": True,
        "epsilon_and_Pstar_Z_independent": True,
        "pressure_axis_condition_status": axis_status,
        "pressure_axis_condition_authenticated": False,
        "pressure_axis_condition_required": (
            "W_axis_zero, Uz1_axis_zero, P1_axis_normalized"
        ),
        "physical_m15_admissible_only_if_axis_contract_passed": True,
        "same_source_identity_checked": False,
        "source_connected_by_this_module": False,
        "P1_is_not_m15": True,
        "raw_P1_consumed_only_as_pressure_port": True,
        "retained_error_data_unauthenticated": True,
        "error_data_consumed": False,
        "interval_or_sample_selection": False,
        "integral_measure": "source-coordinate d rho; epsilon factor retained",
        "integrals_are_raw_ports": True,
        "n1_values_computed": True,
        "computed_from_caller_callbacks": True,
        "n1_source_authenticated": False,
        "result_is_certified": False,
        "physical_rho_domain": "[0,4]",
        "global_moment_closure": False,
        "global_NS_claim": False,
    }
    return N1M15Evaluation(
        rho=rho,
        Z=Z,
        m15=m15,
        ordinary_Z1=ordinary_z1,
        endpoint=endpoint,
        leading_integral=leading_integral,
        leading_integral_Z1=leading_integral_Z1,
        pressure_integral=pressure_integral,
        pressure_integral_Z1=pressure_integral_Z1,
        pressure_boundary=pressure_boundary,
        pressure_boundary_Z1=pressure_boundary_Z1,
        retained_error_data=response.error_data,
        axis_conditions=response.axis_conditions,
        metadata=metadata,
    )


__all__ = [
    "AXIS_RHO",
    "FieldProvider",
    "IntegralOperator",
    "M15_FORMULA",
    "M15_NAME",
    "M15_ORDINARY_Z1_FORMULA",
    "M15_PHYSICAL_FORMULA",
    "M15_PHYSICAL_ORDINARY_Z1_FORMULA",
    "N1M15Evaluation",
    "N1ResponseInterfaceError",
    "PHYSICAL_P_FORMULA",
    "PHYSICAL_RHO_DOMAIN",
    "PHYSICAL_U_FORMULA",
    "PhysicalN1Point",
    "RawN1Point",
    "PrincipalGreenResponse",
    "PrincipalGreenResponseWithZ1",
    "RAW_STATE_AXIS_ORDERS",
    "RAW_STATE_ORDER",
    "RUNTIME_MODEL_METADATA",
    "SCHEMA",
    "SameSourceN1Response",
    "TEMPORAL_ORDER",
    "bind_same_source_n1_response",
    "default_p1_component",
    "default_uz1_component",
    "evaluate_n1_m15",
    "physical_n1_response_at",
    "principal_green_response",
    "principal_green_response_with_Z1",
]


