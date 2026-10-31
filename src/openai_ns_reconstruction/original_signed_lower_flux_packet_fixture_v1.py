"""Independent exact-rational fixture for the signed lower-flux adapter.

The fixture creates a fresh live source packet, flow, basis, ledger, context,
and full target-axis interval.  It does not import or call an original source
provider.  Its arithmetic is a deliberately small exact ``Fraction`` box,
independent of the production evaluator's operation order.  ``run_fixture``
only invokes the evaluator supplied by the caller and reports JSON-safe
diagnostics; it never evaluates a real source function.
"""

from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction
from typing import Any, Callable


SCHEMA = "original_signed_lower_flux_packet_fixture_v1"
SOURCE_FAMILY = "synthetic_signed_lower_flux_source_family_v1"


class ExactInterval:
    """Closed rational interval with exact directed endpoint arithmetic."""

    __slots__ = ("lo", "hi", "ctx")

    def __init__(self, lo: Any, hi: Any | None = None, ctx: Any = None):
        self.lo = _fraction(lo)
        self.hi = self.lo if hi is None else _fraction(hi)
        if self.lo > self.hi:
            raise ValueError("reversed exact interval")
        self.ctx = ctx

    def _coerce(self, other: Any) -> "ExactInterval":
        if isinstance(other, ExactInterval):
            return other
        return ExactInterval(other, ctx=self.ctx)

    def __add__(self, other: Any) -> "ExactInterval":
        other = self._coerce(other)
        return ExactInterval(self.lo + other.lo, self.hi + other.hi, self.ctx)

    __radd__ = __add__

    def __neg__(self) -> "ExactInterval":
        return ExactInterval(-self.hi, -self.lo, self.ctx)

    def __sub__(self, other: Any) -> "ExactInterval":
        return self + -self._coerce(other)

    def __rsub__(self, other: Any) -> "ExactInterval":
        return self._coerce(other) - self

    def __mul__(self, other: Any) -> "ExactInterval":
        other = self._coerce(other)
        values = (self.lo * other.lo, self.lo * other.hi,
                  self.hi * other.lo, self.hi * other.hi)
        return ExactInterval(min(values), max(values), self.ctx)

    __rmul__ = __mul__

    def __truediv__(self, other: Any) -> "ExactInterval":
        other = self._coerce(other)
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError("interval denominator crosses zero")
        values = (self.lo / other.lo, self.lo / other.hi,
                  self.hi / other.lo, self.hi / other.hi)
        return ExactInterval(min(values), max(values), self.ctx)

    def __pow__(self, power: int) -> "ExactInterval":
        if power < 0:
            return ExactInterval(1, ctx=self.ctx) / (self ** (-power))
        result = ExactInterval(1, ctx=self.ctx)
        for _ in range(power):
            result = result * self
        return result

    def __repr__(self) -> str:
        return f"ExactInterval({self.lo!r}, {self.hi!r})"

    def __eq__(self, other: Any) -> bool:
        return isinstance(other, ExactInterval) and self.lo == other.lo and self.hi == other.hi


def _fraction(value: Any) -> Fraction:
    if isinstance(value, Fraction):
        return value
    if isinstance(value, ExactInterval):
        if value.lo != value.hi:
            raise TypeError("a non-singleton interval is not an exact scalar")
        return value.lo
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, str):
        return Fraction(value)
    return Fraction(value)


class FixtureContext:
    def mpf(self, value: Any) -> ExactInterval:
        if isinstance(value, ExactInterval):
            if value.ctx is not None and value.ctx is not self:
                raise ValueError("context mismatch")
            value.ctx = self
            return value
        if isinstance(value, (tuple, list)) and len(value) == 2:
            return ExactInterval(value[0], value[1], ctx=self)
        return ExactInterval(value, ctx=self)


class FormalScale:
    def __init__(self, bases: Any, powers: tuple[Any, ...] = (), offset: Any = 0):
        self.bases = bases
        self.ctx = bases[0].ctx
        self.powers = tuple(powers)
        self.offset = offset

    def pair(self, other: "FormalScale") -> None:
        if other.bases is not self.bases:
            raise ValueError("different source basis")

    def __add__(self, other: "FormalScale") -> "FormalScale":
        self.pair(other)
        return FormalScale(self.bases, self.powers, self.offset)

    def __neg__(self) -> "FormalScale":
        return FormalScale(self.bases, self.powers, -self.offset)

    def __sub__(self, other: "FormalScale") -> "FormalScale":
        return self + -other


class ScaledEnclosure:
    """Exact synthetic analogue of the live source ``ScaledEnclosure``."""

    def __init__(self, scale: FormalScale, coefficient: Any, ledger: Any):
        self.scale = scale
        self.ctx = scale.ctx
        self.coefficient = coefficient if isinstance(coefficient, ExactInterval) \
            else self.ctx.mpf(coefficient)
        self.ledger = ledger

    def scalar(self, value: Any) -> "ScaledEnclosure":
        return ScaledEnclosure(FormalScale(self.scale.bases), value, self.ledger)

    def coerce(self, other: Any) -> "ScaledEnclosure":
        if not isinstance(other, ScaledEnclosure):
            other = self.scalar(other)
        if other.scale.bases is not self.scale.bases or other.ledger is not self.ledger:
            raise ValueError("one source basis and ledger required")
        return other

    def __neg__(self) -> "ScaledEnclosure":
        return ScaledEnclosure(-self.scale, -self.coefficient, self.ledger)

    def __add__(self, other: Any) -> "ScaledEnclosure":
        other = self.coerce(other)
        return ScaledEnclosure(self.scale, self.coefficient + other.coefficient, self.ledger)

    __radd__ = __add__

    def __sub__(self, other: Any) -> "ScaledEnclosure":
        return self + -self.coerce(other)

    def __rsub__(self, other: Any) -> "ScaledEnclosure":
        return self.coerce(other) - self

    def __mul__(self, other: Any) -> "ScaledEnclosure":
        other = self.coerce(other)
        return ScaledEnclosure(self.scale + other.scale,
                               self.coefficient * other.coefficient, self.ledger)

    __rmul__ = __mul__

    @property
    def zero(self) -> bool:
        return self.coefficient.lo == 0 and self.coefficient.hi == 0


class FixtureFlow:
    def __init__(self, ctx: FixtureContext, logs: Any, ledger: Any):
        self.c = ctx
        self.logs = logs
        self.ledger = ledger
        self.ctx = ctx

    def scalar(self, value: Any) -> ScaledEnclosure:
        return ScaledEnclosure(FormalScale(self.logs), value, self.ledger)


class FixtureOp:
    def __init__(self, ctx: FixtureContext, flow: FixtureFlow, P0: Any, Pstar: Any = None):
        self.c = ctx
        self.ctx = ctx
        self.flow = flow
        self.P0 = P0
        self.Pstar = Pstar


def _row(flow: FixtureFlow, value: Any) -> ScaledEnclosure:
    return flow.scalar(value)


def _build_inputs() -> dict[str, Any]:
    ctx = FixtureContext()
    logs = [ctx.mpf(Fraction(1, 7)), ctx.mpf(Fraction(2, 7))]
    ledger: dict[str, Any] = {"fixture": True}
    flow = FixtureFlow(ctx, logs, ledger)
    op = FixtureOp(ctx, flow, tuple(_row(flow, Fraction(3, 11)) for _ in range(4)))

    V = [_row(flow, value) for value in
         (Fraction(1, 2), Fraction(1, 3), Fraction(1, 5), Fraction(1, 7))]
    # Every m row is nonzero; the history certificate is checked independently
    # of the formula so a zero-history shortcut cannot pass.
    m = [_row(flow, value) for value in
         (Fraction(1, 10), Fraction(1, 11), Fraction(1, 13), Fraction(1, 17))]
    Vy = [_row(flow, value) for value in
          (Fraction(1, 4), Fraction(1, 6), Fraction(1, 8), Fraction(1, 10))]
    R = _row(flow, Fraction(3, 2))
    Pstar = _row(flow, Fraction(5, 4))
    op.Pstar = Pstar
    generic = {
        "common_velocity_V_axial5": V,
        "common_own_five_histories_axial5": {"m": m},
        "original_physical_radius": R,
    }
    initial = {
        "exact_zero_correction_follows_from_defining_Duhamel_lower_bound": True,
        "actual_background_incoming_histories_remain_nonzero": True,
        "q_A_B_and_Z_rows_flat_on_same_source_two_sided_left_collar": True,
        "no_saved_frame_or_N1024_initial_correction_used": True,
    }
    packet = {
        "source_identity": SOURCE_FAMILY,
        "candidate_N": 1024,
        "chart": "first_micro",
        "exact_Z_cell": ["0", ".5"],
        "exact_left": ["sc_half"],
        "exact_right": ["sc_half"],
        "exact_common_P0_axial5": op.P0,
        "actual_phase_Z_exact_zero": True,
        "actual_background_source": {"V_y": Vy},
        "original_generic_source": generic,
        "source_owned_modification_support": "original_flat_left_collar",
        "source_owned_collar_q_C0_Z_exact_zero": True,
        "original_active_inverse_not_evaluated_on_source_owned_flat_collar": True,
        "actual_source_owned_initial_condition": initial,
        "original_background_nonzero_histories_remain_separate": True,
        "errors": {"synthetic_signed_error": "retained"},
    }
    view = {
        "V": V,
        "m": m,
        "R": R,
        "Pstar": Pstar,
        "source_packet": packet,
        "errors": {"synthetic_signed_error": "retained"},
    }
    Z = ctx.mpf((Fraction(1, 10), Fraction(1, 5)))
    delta = ctx.mpf(Fraction(1, 4))
    L = ctx.mpf(1) - delta * Z * Z
    assumptions = {
        "query": {
            "source_packet": packet,
            "Z": Z,
            "source": "source_query(first_micro, ('sc_half',), ('sc_half',))",
            "checks": {
                "packet_query_exact": True,
                "coordinate_exact": True,
                "candidate_N_is_source_owned": True,
                "Z_matches_exact_Z_cell": True,
            },
            "details": {
                "chart": "first_micro",
                "exact_left": ["sc_half"],
                "exact_right": ["sc_half"],
                "candidate_N": 1024,
                "exact_Z_cell": ["0", ".5"],
            },
        },
        "z_independence": {
            "source_packet": packet,
            "Z": Z,
            "delta": delta,
            "source": "source_packet exact phase/radius/scale/delta contract",
            "checks": {
                "R_Z_independent": True,
                "Pstar_Z_independent": True,
                "delta_Z_independent": True,
                "full_target_Z_interval": True,
            },
        },
        "history": {
            "source_packet": packet,
            "source": "actual_source_owned_initial_condition",
            "checks": {
                "packet_m_rows_exact": True,
                "original_history_not_zeroed": True,
                "source_owned_flat_collar": True,
            },
        },
        "L_positive": {
            "source_packet": packet,
            "Z": Z,
            "delta": delta,
            "source": "ctx interval arithmetic on requested target Z",
            "checks": {"L_strictly_positive": True},
            "lower_bound": L.lo,
        },
    }
    return {
        "view": view, "op": op, "ctx": ctx, "Z": Z, "delta": delta,
        "assumptions": assumptions, "flow": flow, "packet": packet,
        "V": V, "m": m, "Vy": Vy, "R": R, "Pstar": Pstar,
    }


def build_synthetic_inputs() -> dict[str, Any]:
    """Return fresh inputs for a caller-owned evaluator invocation."""

    return _build_inputs()


build_fixture_inputs = build_synthetic_inputs


def _expected(inputs: Mapping[str, Any]) -> tuple[ExactInterval, ExactInterval]:
    """Independent exact arithmetic for RQy and its ordinary first derivative."""

    view = inputs["view"]
    P = view["Pstar"].coefficient
    R = view["R"].coefficient
    V = view["V"]
    m = view["m"]
    vy = view["source_packet"]["actual_background_source"]["V_y"]
    Z = inputs["Z"]
    delta = inputs["delta"]
    U0, U1, U2 = P * V[0].coefficient, P * V[1].coefficient, P * V[2].coefficient * 2
    M0, M1, M2 = P * m[0].coefficient, P * m[1].coefficient, P * m[2].coefficient * 2
    Uy0, Uy1 = P * vy[0].coefficient, P * vy[1].coefficient
    A, Az, Azz = U0 - M0, U1 - M1, U2 - M2
    one = ExactInterval(1, ctx=inputs["ctx"])
    Z2 = Z * Z
    L = one - delta * Z2
    Ny = Uy0 * (Z * 2) - A * ((one - delta) * Z) - Az * (one - Z2)
    NyZ = (Uy0 * 2 + Uy1 * (Z * 2)
           - (A + Az * Z) * (one - delta)
           + Az * (Z * 2) - Azz * (one - Z2))
    RQy = R * Ny / L
    dRQy = R * (NyZ * L + Ny * (delta * Z * 2)) / L / L
    return RQy, dRQy


def _copy_with_packet(inputs: Mapping[str, Any], packet: Mapping[str, Any], **view_changes: Any) -> dict[str, Any]:
    view = dict(inputs["view"])
    view.update(view_changes)
    view["source_packet"] = packet
    result = dict(inputs)
    result["view"] = view
    result["packet"] = packet
    return result


def _invoke(evaluate: Callable[..., Any], inputs: Mapping[str, Any]) -> Any:
    return evaluate(inputs["view"], inputs["op"], inputs["ctx"], inputs["Z"], inputs["delta"],
                    assumptions=inputs["assumptions"])


def run_fixture(evaluate: Callable[..., Any]) -> dict[str, Any]:
    """Run exact positive and fail-closed negative controls for one evaluator."""

    if not callable(evaluate):
        raise TypeError("evaluate must be callable")
    inputs = _build_inputs()
    checks: dict[str, bool] = {}
    diagnostics: dict[str, Any] = {"source_family": SOURCE_FAMILY,
                                    "real_source_provider_called": False,
                                    "actual_source_validation_claimed": False}
    try:
        result = _invoke(evaluate, inputs)
        checks["returned_mapping"] = isinstance(result, Mapping)
        checks["conditional_status"] = result.get("status") == "conditional_source_enclosure"
        source = result.get("source_enclosures")
        checks["returned_source_enclosures"] = isinstance(source, Mapping) and all(
            isinstance(source.get(name), ScaledEnclosure)
            for name in ("R_Q_y", "dZ_R_Q_y", "ordinary_Z1"))
        checks["ordinary_Z1_retained"] = source.get("ordinary_Z1") is result.get("ordinary_Z1")
        checks["same_basis_and_ledger"] = all(
            source[name].scale.bases is inputs["flow"].logs and
            source[name].ledger is inputs["flow"].ledger
            for name in ("R_Q_y", "dZ_R_Q_y"))
        expected_RQy, expected_Z1 = _expected(inputs)
        checks["exact_RQy_formula"] = source["R_Q_y"].coefficient == expected_RQy
        checks["exact_ordinary_Z1_formula"] = source["ordinary_Z1"].coefficient == expected_Z1
        checks["nonzero_m_history_preserved"] = (
            all(not row.zero for row in inputs["m"])
            and result.get("proof", {}).get("original_m_history_not_zeroed") is True)
        checks["full_target_Z_interval_used"] = result.get("scope", {}).get(
            "target_Z_is_full_requested_interval") is True
        checks["errors_preserved"] = result.get("scope", {}).get("errors_preserved") is True
        checks["native_C1_not_claimed"] = result.get("proof", {}).get(
            "native_or_corrected_C1_identity_proved") is False
        # An independent finite difference of actual synthetic polynomials
        # checks the ordinary-Z derivative, including the factorial on row 2.
        point = _build_inputs()
        point_z = point["ctx"].mpf(Fraction(3, 20))
        point["Z"] = point_z
        for key in ("query", "z_independence", "L_positive"):
            point["assumptions"][key]["Z"] = point_z
        point_L = point["ctx"].mpf(1) - point["delta"] * point_z * point_z
        point["assumptions"]["L_positive"]["lower_bound"] = point_L.lo
        point_result = _invoke(evaluate, point)
        factor = point["Pstar"].coefficient.lo
        radius = point["R"].coefficient.lo
        d = point["delta"].lo
        u0, u1, u2 = (factor * point["V"][k].coefficient.lo * (2 if k == 2 else 1) for k in range(3))
        m0, m1, m2 = (factor * point["m"][k].coefficient.lo * (2 if k == 2 else 1) for k in range(3))
        y0, y1 = (factor * point["Vy"][k].coefficient.lo for k in range(2))
        z0 = Fraction(3, 20)
        h = Fraction(1, 100000)
        def analytic_flux(offset):
            z = z0 + offset
            u = u0 + u1 * offset + u2 * offset * offset / 2
            m = m0 + m1 * offset + m2 * offset * offset / 2
            uz = u1 + u2 * offset
            mz = m1 + m2 * offset
            uy = y0 + y1 * offset
            numerator = 2*z*uy - (1-d)*z*(u-m) - (1-z*z)*(uz-mz)
            return radius * numerator / (1-d*z*z)
        finite_difference = (analytic_flux(h) - analytic_flux(-h)) / (2*h)
        derivative = point_result["ordinary_Z1"].coefficient
        checks["independent_polynomial_Z1_difference"] = (
            derivative.lo == derivative.hi and
            abs(derivative.lo - finite_difference) < Fraction(1, 1000000)
        )
    except Exception as exc:
        diagnostics["positive_error"] = type(exc).__name__ + ": " + str(exc)
        for name in ("returned_mapping", "conditional_status", "returned_source_enclosures",
                     "ordinary_Z1_retained", "same_basis_and_ledger", "exact_RQy_formula",
                     "exact_ordinary_Z1_formula", "nonzero_m_history_preserved",
                     "full_target_Z_interval_used", "errors_preserved", "native_C1_not_claimed"):
            checks[name] = False

    # Wrong basis: one live row is moved to a distinct source basis.
    wrong_basis = _build_inputs()
    other_logs = [wrong_basis["ctx"].mpf(1)]
    wrong_basis["view"] = dict(wrong_basis["view"])
    wrong_V = list(wrong_basis["V"])
    wrong_V[0] = ScaledEnclosure(FormalScale(other_logs), Fraction(1, 2), wrong_basis["flow"].ledger)
    wrong_basis["view"]["V"] = wrong_V
    try:
        _invoke(evaluate, wrong_basis)
    except Exception:
        checks["wrong_basis_rejected"] = True
    else:
        checks["wrong_basis_rejected"] = False

    # Wrong exact source query coordinate.
    wrong_query = _build_inputs()
    bad_packet = dict(wrong_query["packet"])
    bad_packet["exact_left"] = ["wrong_query"]
    wrong_query = _copy_with_packet(wrong_query, bad_packet)
    try:
        _invoke(evaluate, wrong_query)
    except Exception:
        checks["wrong_query_rejected"] = True
    else:
        checks["wrong_query_rejected"] = False

    # Unproved conditions: removing the bound certificate must fail closed.
    unproved = _build_inputs()
    unproved["assumptions"] = dict(unproved["assumptions"])
    unproved["assumptions"].pop("z_independence")
    try:
        _invoke(evaluate, unproved)
    except Exception:
        checks["unproved_Z_independence_rejected"] = True
    else:
        checks["unproved_Z_independence_rejected"] = False

    unproved_L = _build_inputs()
    unproved_L["assumptions"] = dict(unproved_L["assumptions"])
    unproved_L["assumptions"].pop("L_positive")
    try:
        _invoke(evaluate, unproved_L)
    except Exception:
        checks["unproved_L_positivity_rejected"] = True
    else:
        checks["unproved_L_positivity_rejected"] = False

    # Zero the live m history while keeping every packet pointer coherent.
    zeroed = _build_inputs()
    zero_m = [zeroed["flow"].scalar(0) for _ in zeroed["m"]]
    zero_packet = dict(zeroed["packet"])
    zero_generic = dict(zero_packet["original_generic_source"])
    zero_generic["common_own_five_histories_axial5"] = {"m": zero_m}
    zero_packet["original_generic_source"] = zero_generic
    zero_initial = dict(zero_packet["actual_source_owned_initial_condition"])
    zero_initial["actual_background_incoming_histories_remain_nonzero"] = False
    zero_packet["actual_source_owned_initial_condition"] = zero_initial
    zeroed["view"] = dict(zeroed["view"], m=zero_m, source_packet=zero_packet)
    zeroed["packet"] = zero_packet
    zeroed["assumptions"] = dict(zeroed["assumptions"])
    for key, cert in tuple(zeroed["assumptions"].items()):
        if isinstance(cert, Mapping):
            updated = dict(cert, source_packet=zero_packet)
            if key == "history":
                updated["checks"] = dict(updated["checks"], original_history_not_zeroed=False)
            zeroed["assumptions"][key] = updated
    try:
        _invoke(evaluate, zeroed)
    except Exception:
        checks["zeroed_history_rejected"] = True
    else:
        checks["zeroed_history_rejected"] = False

    diagnostics.update({
        "ordinary_Z_domain": ["1/10", "1/5"],
        "delta": "1/4",
        "candidate_N": 1024,
        "m_history_rows": 4,
        "checks_count": len(checks),
    })
    return {"schema": SCHEMA, "passed": all(checks.values()),
            "checks": checks, "diagnostics": diagnostics,
            "synthetic_only": True,
            "real_source_validation_claimed": False}


__all__ = [
    "SCHEMA", "SOURCE_FAMILY", "ExactInterval", "FormalScale", "ScaledEnclosure",
    "FixtureContext", "FixtureFlow", "FixtureOp", "build_synthetic_inputs",
    "build_fixture_inputs", "run_fixture",
]
