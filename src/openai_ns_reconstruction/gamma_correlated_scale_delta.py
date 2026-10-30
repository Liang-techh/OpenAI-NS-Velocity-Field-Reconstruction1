"""Exact correlated Gamma scale-difference kernels.

This module is a small, source-independent mathematical leaf.  It does not
load the historical Gamma graph and it does not form a difference of two
large P/F source values.  The exact source-facing definition is retained by
``formal_source_expression``.  The point API evaluates only its normalized
correlated kernels on caller-supplied, finite, moderate point inputs.

For ``a = delta/2``, ``d = 1-Z**2``, ``k = 2*d*exp(-t)``, and

    w(v) = exp(-v) * v**(a+1) / Gamma(1+a),

the normalized kernels are

    K0 = -k**2*(a+1) * int_theta int_s int_v
         w*v*s*(1+s*k*v*(SF+theta*DeltaS))**(-a-2),

    K1 = 4*Z*exp(-t)*(a+1)*k * int_theta int_v
         w*v*(1+k*v*(SF+theta*DeltaS))**(-a-2).

The returned source forms are ``DeltaS*K0`` and ``DeltaS*K1`` only as
unevaluated formal products.  No actual source value, common P/F ownership,
same-SF/DeltaS assertion, or future-selector closure is claimed here.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Mapping

import mpmath as mp


SCHEMA = "gamma_correlated_scale_delta_v1"
COMPONENT = SCHEMA
PHASE = "gamma_correlated_scale_delta_pending"
SOURCE_VALUE_FLAGS = {
    "source_values_computed": False,
    "source_value_computation": False,
    "actual_source_claimed": False,
}


class GammaCorrelatedScaleDeltaInputError(ValueError):
    """Raised when a point or gate input is not a finite valid datum."""


def _mp(value: Any, label: str) -> mp.mpf:
    """Convert one caller scalar without silently accepting an interval."""

    if isinstance(value, bool) or isinstance(value, (bytes, bytearray)):
        raise GammaCorrelatedScaleDeltaInputError(f"{label} must be a real point")
    if getattr(value, "_mpi_", None) is not None:
        raise GammaCorrelatedScaleDeltaInputError(f"{label} must not be an interval")
    try:
        if isinstance(value, Fraction):
            result = mp.mpf(value.numerator) / value.denominator
        else:
            result = mp.mpf(value)
    except (TypeError, ValueError, ArithmeticError) as exc:
        raise GammaCorrelatedScaleDeltaInputError(f"{label} must be a real point") from exc
    if not mp.isfinite(result):
        raise GammaCorrelatedScaleDeltaInputError(f"{label} must be finite")
    return result


def _bool(value: Any, label: str) -> bool:
    if not isinstance(value, bool):
        raise GammaCorrelatedScaleDeltaInputError(f"{label} must be an explicit bool")
    return value


def _nstr(value: Any) -> str:
    if isinstance(value, mp.mpf):
        return mp.nstr(value, 30)
    return str(value)


def _relerr(left: mp.mpf, right: mp.mpf) -> mp.mpf:
    scale = max(mp.mpf(1), abs(left), abs(right))
    return abs(left - right) / scale


def _weight(a: mp.mpf, v: mp.mpf) -> mp.mpf:
    return mp.exp(-v) * mp.power(v, a + 1) / mp.gamma(a + 1)


def _f(a: mp.mpf, q: mp.mpf) -> mp.mpf:
    """Stable continuation of f_a(q) at q = 0."""

    if q == 0:
        return mp.mpf(1)
    if a == 0:
        return mp.log1p(q) / q
    return -mp.expm1(-a * mp.log1p(q)) / (a * q)


def _h(a: mp.mpf, q: mp.mpf) -> mp.mpf:
    return mp.power(1 + q, -a - 1)


def _log1prel(x: mp.mpf) -> mp.mpf:
    """Stable log1p(x)/x with its exact zero continuation."""

    return mp.mpf(1) if x == 0 else mp.log1p(x) / x


def _exprel(x: mp.mpf) -> mp.mpf:
    """Stable expm1(x)/x with its exact zero continuation."""

    return mp.mpf(1) if x == 0 else mp.expm1(x) / x


def _theta_average(
    a: mp.mpf,
    k: mp.mpf,
    v: mp.mpf,
    s: mp.mpf,
    sf: mp.mpf,
    delta_s: mp.mpf,
) -> mp.mpf:
    """Average (1+q+theta*b)^(-p) without a P/F difference.

    Here p = a+2, q = s*k*v*SF, and b = s*k*v*DeltaS.  With
    r=b/(1+q), the exact theta integral is

        (1+q)^(-p) log1prel(r) exprel((1-p)*log1p(r)).

    The two relative functions carry their continuous value one at r = 0,
    so DeltaS -> 0 never requests extra precision and never subtracts two
    nearly equal source kernels.
    """

    p = a + 2
    q = s * k * v * sf
    base = 1 + q
    r = s * k * v * delta_s / base
    if r == 0:
        return mp.power(base, -p)
    # If r rounds to -1 for SP=0 and huge v, the true positive endpoint
    # ratio is recovered without subtracting nearly equal values.
    log_ratio = (mp.log1p(s*k*v*(sf+delta_s)) - mp.log1p(q)) if r <= -1 else mp.log1p(r)
    return mp.power(base, -p) * (log_ratio/r) * _exprel((1-p)*log_ratio)


def _quad_weighted(a: mp.mpf, body: Any) -> mp.mpf:
    """One improper v integral with the normalized Gamma weight."""

    return mp.quad(lambda v: _weight(a, v) * body(v), [0, 1, mp.inf])


def _compute_kernels(
    a: mp.mpf,
    z: mp.mpf,
    t: mp.mpf,
    sf: mp.mpf,
    delta_s: mp.mpf,
    dps: int,
) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    """Evaluate K0/K1 from direct correlated kernels.

    Theta is integrated analytically by ``_theta_average``.  K0 then has one
    finite s integral nested inside its v integral; K1 uses s=1 directly.
    No f/h divided difference and no DeltaS-dependent precision expansion is
    used.  At k=0 the two delta kernels are exactly zero, while the separate
    ordinary-Z boundary derivative remains available below.
    """

    with mp.workdps(dps):
        k = 2 * (1 - z * z) * mp.exp(-t)
        if k == 0:
            return mp.mpf(0), mp.mpf(0), k

        def k0_body(v: mp.mpf) -> mp.mpf:
            s_integral = mp.quad(
                lambda s: s * _theta_average(a, k, v, s, sf, delta_s),
                [0, 1],
            )
            return -k * k * (a + 1) * v * s_integral

        def k1_body(v: mp.mpf) -> mp.mpf:
            theta_at_s_one = _theta_average(a, k, v, mp.mpf(1), sf, delta_s)
            return 4 * z * mp.exp(-t) * (a + 1) * k * v * theta_at_s_one

        k0 = _quad_weighted(a, k0_body)
        k1 = _quad_weighted(a, k1_body)
        return +k0, +k1, +k


def _validate_core(
    a: Any,
    z: Any,
    t: Any,
    sf: Any,
    delta_s: Any,
) -> tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    aa, zz, tt = (_mp(a, "a"), _mp(z, "Z"), _mp(t, "t"))
    ss, dd = _mp(sf, "SF"), _mp(delta_s, "DeltaS")
    if aa < 0:
        raise GammaCorrelatedScaleDeltaInputError("a must be nonnegative")
    if zz < -1 or zz > 1:
        raise GammaCorrelatedScaleDeltaInputError("Z must lie in [-1, 1]")
    if tt < 0:
        raise GammaCorrelatedScaleDeltaInputError("t must be nonnegative")
    if ss < 0:
        raise GammaCorrelatedScaleDeltaInputError("SF must be nonnegative")
    sp = ss + dd
    if sp < 0:
        raise GammaCorrelatedScaleDeltaInputError("SP = SF + DeltaS must be nonnegative")
    k = 2 * (1 - zz * zz) * mp.exp(-tt)
    if k < 0 or not mp.isfinite(k):
        raise GammaCorrelatedScaleDeltaInputError("k must be finite and nonnegative")
    # Return k as a validated point too; the compute path recomputes it at its
    # requested precision rather than using a rounded caller-side value.
    return aa, zz, tt, ss, dd, k


def _gate_record(
    *,
    a: mp.mpf,
    z: mp.mpf,
    t: mp.mpf,
    sf: mp.mpf,
    delta_s: mp.mpf,
    dps: int,
    L: Any,
    actual_mu: Any,
    T: Any,
    T_small: Any,
    ordinary_z_independent: bool,
) -> dict[str, Any]:
    sp = sf + delta_s
    result: dict[str, Any] = {
        "a_nonnegative": bool(a >= 0),
        "Z_in_closed_unit_interval": bool(-1 <= z <= 1),
        "t_finite": bool(mp.isfinite(t)),
        "SF_nonnegative_finite": bool(sf >= 0 and mp.isfinite(sf)),
        "DeltaS_finite": bool(mp.isfinite(delta_s)),
        "SP_nonnegative_finite": bool(sp >= 0 and mp.isfinite(sp)),
        "positive_branch_1_plus_kvS": True,
        "precision_explicit": isinstance(dps, int) and dps > 0,
        "ordinary_Z_independence_assumed": ordinary_z_independent,
        "a_Z_independence_assumed": ordinary_z_independent,
        "K1_ordinary_Z_usable": ordinary_z_independent,
        "Z_endpoint_boundary_retained": bool(z in (-1, 1)),
        "L_positive_finite": False,
        "actual_mu_positive_finite": False,
        "T_is_separate_from_actual_mu": True,
        "T_small_declared": False,
    }
    if L is not None:
        try:
            ll = _mp(L, "L")
            result["L_positive_finite"] = bool(ll > 0)
        except GammaCorrelatedScaleDeltaInputError:
            result["L_positive_finite"] = False
    if actual_mu is not None:
        try:
            mm = _mp(actual_mu, "actual_mu")
            result["actual_mu_positive_finite"] = bool(mm > 0)
        except GammaCorrelatedScaleDeltaInputError:
            result["actual_mu_positive_finite"] = False
    if T_small is not None:
        result["T_small_declared"] = _bool(T_small, "T_small")
    if T is not None:
        _mp(T, "T")
    required = (
        "a_nonnegative",
        "Z_in_closed_unit_interval",
        "t_finite",
        "SF_nonnegative_finite",
        "DeltaS_finite",
        "SP_nonnegative_finite",
        "positive_branch_1_plus_kvS",
        "precision_explicit",
        "ordinary_Z_independence_assumed",
        "T_is_separate_from_actual_mu",
    )
    # Endpoint retention, T-small, L, and actual-mu are metadata rather than
    # unrelated numerical inputs required by the normalized kernel algorithm.
    result["application_gates_satisfied"] = all(result[key] for key in required)
    return result


@dataclass(frozen=True)
class KernelBoundEstimates:
    """Analytic inequality bounds with rounded numerical endpoints.

    The inequalities are analytic, but mpmath's ordinary arithmetic does not
    provide directed rounding.  Therefore these endpoints are estimates and
    are deliberately never marked as certified bounds.
    """

    K0_abs_upper: mp.mpf
    K1_abs_upper: mp.mpf
    certified: bool = False
    bound_kind: str = "analytic inequality with rounded bound estimate"
    source_values_computed: bool = False

    def as_dict(self) -> dict[str, Any]:
        return {
            "K0_abs_upper": _nstr(self.K0_abs_upper),
            "K1_abs_upper": _nstr(self.K1_abs_upper),
            "certified": self.certified,
            "bound_kind": self.bound_kind,
            "analytic_inequality": True,
            "source_values_computed": self.source_values_computed,
            "basis": "w>=0, 1+s*k*v*S>=1, 1+k*v*S>=1",
        }


@dataclass(frozen=True)
class GammaKernelPoint:
    """A normalized point evaluation; it is not a source-value receipt."""

    a: mp.mpf
    Z: mp.mpf
    t: mp.mpf
    SF: mp.mpf
    DeltaS: mp.mpf
    SP: mp.mpf
    k: mp.mpf
    K0: mp.mpf
    K1: mp.mpf
    dps: int
    gates: Mapping[str, Any]
    bound_estimates: KernelBoundEstimates
    estimated_error: Mapping[str, Any] | None = None
    numeric_evaluation: bool = True
    certified: bool = False
    normalized_kernel_only: bool = True
    moderate_point_only: bool = True
    not_for_actual_double_exponential: bool = True
    source_values_computed: bool = False
    source_value_computation: bool = False

    @property
    def k0(self) -> mp.mpf:
        return self.K0

    @property
    def k1(self) -> mp.mpf:
        return self.K1

    def __getitem__(self, key: str) -> Any:
        aliases = {"k0": "K0", "k1": "K1", "delta_s": "DeltaS", "sf": "SF", "sp": "SP"}
        return getattr(self, aliases.get(key, key))

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "a": _nstr(self.a),
            "Z": _nstr(self.Z),
            "t": _nstr(self.t),
            "SF": _nstr(self.SF),
            "DeltaS": _nstr(self.DeltaS),
            "SP": _nstr(self.SP),
            "k": _nstr(self.k),
            "K0": _nstr(self.K0),
            "K1": _nstr(self.K1),
            "dps": self.dps,
            "gates": dict(self.gates),
            "bound_estimates": self.bound_estimates.as_dict(),
            "estimated_error": dict(self.estimated_error or {}),
            "numeric_evaluation": self.numeric_evaluation,
            "certified": self.certified,
            "normalized_kernel_only": self.normalized_kernel_only,
            "moderate_point_only": self.moderate_point_only,
            "not_for_actual_double_exponential": self.not_for_actual_double_exponential,
            "source_values_computed": self.source_values_computed,
            "source_value_computation": self.source_value_computation,
        }


def analytic_bound_estimates(a: Any, z: Any, t: Any, *, dps: int = 50) -> KernelBoundEstimates:
    """Return analytic majorants with explicitly non-certified endpoints.

    The inequalities hold symbolically for all nonnegative SF and SP, but the
    ordinary mpmath endpoints are rounded estimates, not directed-rounding
    certificates and not bounds on a materialized DeltaS times a kernel.
    """

    if not isinstance(dps, int) or isinstance(dps, bool) or dps < 15:
        raise GammaCorrelatedScaleDeltaInputError("dps must be an integer >= 15")
    with mp.workdps(dps):
        aa, zz, tt = _mp(a, "a"), _mp(z, "Z"), _mp(t, "t")
        if aa < 0 or zz < -1 or zz > 1 or tt < 0:
            raise GammaCorrelatedScaleDeltaInputError("invalid a, Z, or branch")
        k = 2 * (1 - zz * zz) * mp.exp(-tt)
        if not mp.isfinite(k):
            raise GammaCorrelatedScaleDeltaInputError("k must be finite")
        moment_v = (aa + 1) * (aa + 2)
        k0 = k * k * (aa + 1) * moment_v / 2
        k1 = 4 * abs(zz) * mp.exp(-tt) * (aa + 1) * k * moment_v
        return KernelBoundEstimates(+k0, +k1, certified=False)


def evaluate_kernels(
    a: Any,
    z: Any,
    t: Any,
    SF: Any,
    DeltaS: Any,
    *,
    dps: int = 50,
    L: Any = None,
    precision: Any = None,
    actual_mu: Any = None,
    mu: Any = None,
    T: Any = None,
    T_small: Any = None,
    ordinary_z_independent: bool | None = None,
    ordinary_Z_independent: bool | None = None,
    sf_sp_z_independent: bool | None = None,
    estimate_error: bool = False,
) -> GammaKernelPoint:
    """Evaluate the normalized K0/K1 point kernels at arbitrary precision.

    ``L``, ``actual_mu``, and ``T`` are caller gates/metadata only; the API
    never reconstructs a radius or obtains S from mu.  Set the explicit
    ordinary-Z independence declaration before requesting K1; an omitted
    declaration is rejected rather than silently treating K1 as ordinary-Z.
    ``bound_estimates`` contains analytic inequalities with rounded,
    explicitly non-certified endpoints; the returned point values and
    optional errors are numerical estimates.
    """

    declared_ordinary = ordinary_z_independent
    if ordinary_Z_independent is not None:
        alias_value = _bool(ordinary_Z_independent, "ordinary_Z_independent")
        if declared_ordinary is not None and _bool(declared_ordinary, "ordinary_z_independent") != alias_value:
            raise GammaCorrelatedScaleDeltaInputError(
                "ordinary-Z independence declarations disagree"
            )
        declared_ordinary = alias_value
    if sf_sp_z_independent is not None:
        alias_value = _bool(sf_sp_z_independent, "sf_sp_z_independent")
        if declared_ordinary is not None and _bool(declared_ordinary, "ordinary_z_independent") != alias_value:
            raise GammaCorrelatedScaleDeltaInputError(
                "ordinary-Z independence declarations disagree"
            )
        declared_ordinary = alias_value
    if declared_ordinary is None:
        raise GammaCorrelatedScaleDeltaInputError(
            "ordinary_z_independent must be explicitly declared before evaluating K1"
        )
    ordinary_z_independent = _bool(declared_ordinary, "ordinary_z_independent")
    if ordinary_z_independent is not True:
        raise GammaCorrelatedScaleDeltaInputError("a, SF and SP must be declared ordinary-Z independent")
    if mu is not None:
        if actual_mu is not None and _mp(mu, "mu") != _mp(actual_mu, "actual_mu"):
            raise GammaCorrelatedScaleDeltaInputError("mu and actual_mu disagree")
        actual_mu = mu
    if precision is not None:
        if isinstance(precision, bool) or not isinstance(precision, int) or precision < 15:
            raise GammaCorrelatedScaleDeltaInputError("precision must be an integer >= 15")
        dps = precision
    if isinstance(dps, bool) or not isinstance(dps, int) or dps < 15:
        raise GammaCorrelatedScaleDeltaInputError("dps must be an integer >= 15")

    with mp.workdps(dps):
        # Convert decimal/string caller inputs inside the requested context;
        # converting before workdps would silently round an arbitrary-
        # precision point to the ambient mpmath precision.
        aa, zz, tt, ss, dd, _ = _validate_core(a, z, t, SF, DeltaS)
        gates = _gate_record(
            a=aa,
            z=zz,
            t=tt,
            sf=ss,
            delta_s=dd,
            dps=dps,
            L=L,
            actual_mu=actual_mu,
            T=T,
            T_small=T_small,
            ordinary_z_independent=ordinary_z_independent,
        )
        k0, k1, kk = _compute_kernels(aa, zz, tt, ss, dd, dps)
        bounds = analytic_bound_estimates(aa, zz, tt, dps=dps)

    estimated: dict[str, Any] | None = None
    if estimate_error:
        low_dps = max(15, dps - max(8, dps // 5))
        low = evaluate_kernels(
            aa,
            zz,
            tt,
            ss,
            dd,
            dps=low_dps,
            L=L,
            actual_mu=actual_mu,
            T=T,
            T_small=T_small,
            ordinary_z_independent=ordinary_z_independent,
            estimate_error=False,
        )
        estimated = {
            "K0": _nstr(_relerr(k0, low.K0)),
            "K1": _nstr(_relerr(k1, low.K1)),
            "kind": "precision-doubling-style estimate, not certified",
        }
    return GammaKernelPoint(
        aa,
        zz,
        tt,
        ss,
        dd,
        ss + dd,
        kk,
        k0,
        k1,
        dps,
        gates,
        bounds,
        estimated,
    )


def point_kernels(*args: Any, **kwargs: Any) -> GammaKernelPoint:
    """Alias for the arbitrary-precision normalized point API."""

    return evaluate_kernels(*args, **kwargs)


def solve_kernels(*args: Any, **kwargs: Any) -> GammaKernelPoint:
    """Compatibility alias emphasizing that only K0/K1 are solved."""

    return evaluate_kernels(*args, **kwargs)


def ordinary_z_derivative_kernel(
    a: Any,
    z: Any,
    t: Any,
    S: Any,
    *,
    dps: int = 50,
) -> mp.mpf:
    """Evaluate the ordinary-Z derivative kernel, retaining the k=0 boundary.

    This is a point-valued normalized reference used by the synthetic check;
    it is not a declaration that a historical actual source has been bound.
    """

    aa, zz, tt, ss, _, kk = _validate_core(a, z, t, S, 0)
    if isinstance(dps, bool) or not isinstance(dps, int) or dps < 15:
        raise GammaCorrelatedScaleDeltaInputError("dps must be an integer >= 15")
    with mp.workdps(dps):
        value = -4 * zz * mp.exp(-tt) * _quad_weighted(
            aa, lambda v: _h(aa, kk * v * ss)
        )
        return +value


def formal_source_expression(
    *,
    L: Any = None,
    precision: Any = None,
    actual_mu: Any = None,
    T_small: Any = None,
    sf_owner: str = "F",
    sp_owner: str = "P",
    delta_s_owner: str = "caller_delta",
    ordinary_z_independent: bool = False,
) -> dict[str, Any]:
    """Return the exact unevaluated correlated source definition.

    The symbols ``SF`` and ``DeltaS`` are intentionally not read from caller
    data here.  Owner and application gates are recorded explicitly so this
    expression cannot be mistaken for a proof of common P/F source ownership.
    """

    owner_ok = all(
        isinstance(value, str) and bool(value.strip())
        for value in (sf_owner, sp_owner, delta_s_owner)
    )
    ordinary_gate = _bool(ordinary_z_independent, "ordinary_z_independent")
    if T_small is not None:
        T_small = _bool(T_small, "T_small")
    return {
        "schema": SCHEMA,
        "kind": "exact_unevaluated_correlated_source_definition",
        "definitions": {
            "a": "delta/2",
            "d": "1-Z^2",
            "k": "2*(1-Z^2)*exp(-t)",
            "t_domain": "t >= 0",
            "w(v)": "exp(-v)*v^(a+1)/Gamma(1+a)",
            "S_theta": "SF + theta*DeltaS",
            "SF_domain": "SF >= 0",
            "SP_domain": "SP = SF + DeltaS >= 0",
            "branch": "1 + s*k*v*S_theta > 0 and 1 + k*v*S_theta > 0",
        },
        "normalized_kernel": {
            "K0": "-k^2*(a+1)*Integral_theta[0,1] Integral_s[0,1] Integral_v[0,inf] "
            "w(v)*v*s*(1+s*k*v*(SF+theta*DeltaS))^(-a-2)",
            "K1": "+4*Z*exp(-t)*(a+1)*k*Integral_theta[0,1] Integral_v[0,inf] "
            "w(v)*v*(1+k*v*(SF+theta*DeltaS))^(-a-2)",
            "analytic_reduction": {
                "s_integral": "Integral_s s*(1+s*q)^(-a-2) ds",
                "theta_average": "A^(-p)*log1prel(r)*exprel((1-p)*log1p(r))",
                "theta_variables": "p=a+2, A=1+s*k*v*SF, r=s*k*v*DeltaS/A",
                "quadrature_dimension": "K0: v with inner s; K1: v only",
                "log1prel": "log1p(r)/r, continued by 1 at r=0",
                "exprel": "expm1(y)/y, continued by 1 at y=0",
                "K0_direct_reduced": "-k^2*(a+1)*Integral_v w*v*Integral_s s*theta_average ds",
                "K1_direct_reduced": "4*Z*exp(-t)*(a+1)*k*Integral_v w*v*theta_average(s=1)",
                "DeltaS_zero_continuation": "K0=-k^2*(a+1)*Integral_v w*v*Integral_s "
                "s*(1+s*k*v*SF)^(-a-2) ds; "
                "K1=4*Z*exp(-t)*(a+1)*k*Integral_v w*v*(1+k*v*SF)^(-a-2)",
            },
        },
        "source_forms": {
            "Delta_Dhat": "DeltaS*K0",
            "Delta_ordinary_Z_Dhat": "DeltaS*K1",
        },
        "retained_inputs": {
            "L_supplied": L is not None,
            "precision_supplied": precision is not None,
            "actual_mu_supplied": actual_mu is not None,
            "T_small_supplied": T_small is not None,
            "T_is_not_actual_mu": True,
            "values_materialized": False,
        },
        "owner_gates": {
            "SF_owner_declared": owner_ok,
            "SP_owner_declared": owner_ok,
            "SF_owner": sf_owner,
            "SP_owner": sp_owner,
            "DeltaS_owner": delta_s_owner,
            "SF_SP_same_owner_assumed": False,
            "SF_SP_Z_independence_assumed": ordinary_gate,
        },
        "gates": {
            "positive_finite_S_and_branch_required": True,
            "ordinary_Z_independence_required_for_K1": True,
            "a_Z_independence_required_for_K1": True,
            "ordinary_Z_independence_declared": ordinary_gate,
            "L_precision_mu_T_gates_explicit": True,
            "Z_equals_pm1_boundary_retained": True,
        },
        "claims": {
            "P_F_common_source_proved": False,
            "same_SF_DeltaS_proved": False,
            "actual_source_value_computed": False,
            "future_selector_closure_proved": False,
            "certified_point_value": False,
        },
        **SOURCE_VALUE_FLAGS,
    }


def formal_normalized_kernel(**kwargs: Any) -> dict[str, Any]:
    """Return only the exact normalized-kernel portion of the formal packet."""

    return formal_source_expression(**kwargs)["normalized_kernel"]


def _synthetic_dhat(a: mp.mpf, z: mp.mpf, t: mp.mpf, S: mp.mpf) -> mp.mpf:
    _, _, _, ss, _, k = _validate_core(a, z, t, S, 0)
    return k * _quad_weighted(a, lambda v: _f(a, k * v * ss))


def _synthetic_fixture_once() -> dict[str, Any]:
    """Independent moderate fixture used only by the one-shot runner."""

    dps = 22
    with mp.workdps(dps + 8):
        a = mp.mpf("0.7")
        z = mp.mpf("0.35")
        t = mp.mpf("0.4")
        sf = mp.mpf("0.8")
        delta_s = mp.mpf("0.25")
        point = evaluate_kernels(
            a,
            z,
            t,
            sf,
            delta_s,
            dps=dps,
            L=mp.mpf("1.2"),
            precision=dps,
            actual_mu=mp.mpf("0.4"),
            T=mp.mpf("0.01"),
            T_small=False,
            ordinary_z_independent=True,
        )
        sf_value = _synthetic_dhat(a, z, t, sf)
        sp_value = _synthetic_dhat(a, z, t, sf + delta_s)
        direct_d0 = sp_value - sf_value
        direct_z0 = ordinary_z_derivative_kernel(a, z, t, sf, dps=dps)
        direct_z1 = ordinary_z_derivative_kernel(a, z, t, sf + delta_s, dps=dps)
        direct_d1 = direct_z1 - direct_z0
        h = mp.mpf("1e-6")
        fd_z = (_synthetic_dhat(a, z + h, t, sf) - _synthetic_dhat(a, z - h, t, sf)) / (2 * h)
        z_fd_error = _relerr(fd_z, direct_z0)

        zero = evaluate_kernels(
            a, z, t, sf, 0, dps=dps, ordinary_z_independent=True,
        )
        tiny = evaluate_kernels(
            a, z, t, sf, mp.mpf("1e-20"), dps=dps, ordinary_z_independent=True,
        )
        continuity_error = {
            "K0": _relerr(zero.K0, tiny.K0),
            "K1": _relerr(zero.K1, tiny.K1),
        }

        endpoints: dict[str, Any] = {}
        endpoint_errors: dict[str, mp.mpf] = {}
        for endpoint in (-1, 1):
            ep = evaluate_kernels(
                a, endpoint, t, sf, delta_s, dps=dps,
                ordinary_z_independent=True,
            )
            boundary = ordinary_z_derivative_kernel(a, endpoint, t, sf, dps=dps)
            expected = -4 * endpoint * mp.exp(-t) * (a + 1)
            endpoint_errors[f"Z{endpoint}_K0_zero"] = _relerr(ep.K0, mp.mpf(0))
            endpoint_errors[f"Z{endpoint}_K1_zero"] = _relerr(ep.K1, mp.mpf(0))
            endpoint_errors[f"Z{endpoint}_ordinary_boundary"] = _relerr(boundary, expected)
            endpoints[str(endpoint)] = {
                "K0": _nstr(ep.K0),
                "K1": _nstr(ep.K1),
                "ordinary_Z_boundary": _nstr(boundary),
                "boundary_error_estimate": _nstr(_relerr(boundary, expected)),
                "boundary_retained": True,
            }

        errors = {
            "K0_delta_product": _relerr(direct_d0, delta_s * point.K0),
            "K1_delta_product": _relerr(direct_d1, delta_s * point.K1),
            "ordinary_Z_finite_difference": z_fd_error,
            "dS_to_zero_K0": continuity_error["K0"],
            "dS_to_zero_K1": continuity_error["K1"],
        }
        errors.update(endpoint_errors)
        checks = all(value < mp.mpf("1e-9") for value in errors.values())
        return {
            "schema": SCHEMA,
            "fixture": "moderate_independent_stable_definition_difference",
            "checks_passed": bool(checks),
            "error_kind": "estimated, not certified",
            "synthetic_reference_values_computed": True,
            "historical_source_values_computed": False,
            "errors": {key: _nstr(value) for key, value in errors.items()},
            "endpoint_checks": endpoints,
            "point": point.as_dict(),
            "numeric_evaluation": True,
            "certified": False,
            "bound_estimates_are_separate": True,
            "source_values_computed": False,
            "source_value_computation": False,
            "claims": {
                "P_F_common_source_proved": False,
                "same_SF_DeltaS_proved": False,
                "actual_source_value_computed": False,
                "future_selector_closure_proved": False,
            },
        }


def _actual_gamma_binding(root: Any, hashes: Mapping[str, str]) -> dict[str, Any]:
    """Authenticate source bytes before selecting direct archived role fields."""
    import gzip
    import hashlib
    import json
    from pathlib import Path
    root = Path(root)
    stem = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions"
    paths = [stem + ".py", stem + "_check.json", stem + ".json.gz"]
    raw = {}
    for path in paths:
        data = (root/path).read_bytes()
        if hashlib.sha256(data).hexdigest() != hashes.get(path):
            raise ValueError("frozen Gamma raw SHA mismatch: " + path)
        raw[path] = data
    check = json.loads(raw[paths[1]])
    archive_sha = hashlib.sha256(raw[paths[2]]).hexdigest()
    if check.get("all_passed") is not True or check.get("current_original_n1_exact_Gamma_three_future_C1_functions_installed") is not True:
        raise ValueError("Gamma source receipt is not accepted")
    if check.get("input_hashes", {}).get(Path(paths[2]).name) != archive_sha:
        raise ValueError("Gamma receipt does not pin the archive")
    archive = json.loads(gzip.decompress(raw[paths[2]]))
    if archive.get("source_family") != check.get("source_family"):
        raise ValueError("Gamma archive and receipt source families differ")
    binding = archive["function_bindings"]
    return {"archive": paths[2], "archive_sha256": archive_sha,
            "source_file": paths[0], "source_sha256": hashes[paths[0]],
            "receipt": paths[1], "receipt_sha256": hashes[paths[1]],
            "source_family": check["source_family"],
            "source_namespace": "historical_Gamma_archive",
            "source_roles": {"schedule.S": binding["schedule"]["S"],
                             "functions.Dhat": binding["functions"]["Dhat"],
                             "functions.energy": binding["functions"]["energy"]},
            "raw_inputs_authenticated_before_decompression": True,
            "source_refs_are_not_values": True,
            "P_F_owner_correspondence_proved": False,
            "source_graph_traversed": False}


def verify_synthetic_fixture(root: Any, hashes: Mapping[str, str]) -> dict[str, Any]:
    """One-shot source metadata binding and independent synthetic checks."""
    binding = _actual_gamma_binding(root, hashes)
    result = _synthetic_fixture_once()
    result["actual_gamma_binding"] = binding
    result["formal_source_expression"] = formal_source_expression()
    return result


# Explicit aliases used by callers that name the two kernels directly.
evaluate_K0_K1 = evaluate_kernels
evaluate_correlated_kernels = evaluate_kernels
point_K0_K1 = evaluate_kernels
build_formal_source_expression = formal_source_expression
analytic_kernel_bounds = analytic_bound_estimates


__all__ = [
    "SCHEMA",
    "COMPONENT",
    "PHASE",
    "GammaCorrelatedScaleDeltaInputError",
    "KernelBoundEstimates",
    "GammaKernelPoint",
    "analytic_bound_estimates",
    "analytic_kernel_bounds",
    "evaluate_kernels",
    "point_kernels",
    "solve_kernels",
    "evaluate_K0_K1",
    "evaluate_correlated_kernels",
    "point_K0_K1",
    "ordinary_z_derivative_kernel",
    "formal_source_expression",
    "formal_normalized_kernel",
    "build_formal_source_expression",
    "verify_synthetic_fixture",
]
