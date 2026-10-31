"""Fresh exact-rational fixture for the v3 signed lower-flux adapter.

Only synthetic objects are built here.  No source provider is imported, no
source graph is traversed, and no enclosure is converted to a point value.
Normalized generic ``V`` and ``m`` receive ``Pstar`` once.  The packet's
``actual_background_source['V_y']`` rows are already raw physical ``U_y`` and
are passed through without another factor.  ``run_fixture`` invokes only the
evaluator supplied by its caller and returns JSON-safe diagnostics.
"""

from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction
from typing import Any, Callable


SCHEMA = "original_signed_lower_flux_packet_fixture_v3"
SOURCE_FAMILY = "synthetic_signed_lower_flux_source_family_v3_fresh_alpha"


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
        return other if isinstance(other, ExactInterval) else ExactInterval(other, ctx=self.ctx)

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
    """Exact synthetic analogue of the live source scaled enclosure."""

    def __init__(self, scale: FormalScale, coefficient: Any, ledger: Any):
        self.scale = scale
        self.ctx = scale.ctx
        self.coefficient = (coefficient if isinstance(coefficient, ExactInterval)
                            else self.ctx.mpf(coefficient))
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
        self.ctx = ctx
        self.logs = logs
        self.ledger = ledger

    def scalar(self, value: Any) -> ScaledEnclosure:
        return ScaledEnclosure(FormalScale(self.logs), value, self.ledger)


class FixtureOp:
    def __init__(self, ctx: FixtureContext, flow: FixtureFlow, P0: Any,
                 Pstar: Any = None):
        self.c = ctx
        self.ctx = ctx
        self.flow = flow
        self.P0 = P0
        self.Pstar = Pstar


def _row(flow: FixtureFlow, value: Any) -> ScaledEnclosure:
    return flow.scalar(value)


def _build_inputs() -> dict[str, Any]:
    """Build one fresh packet with mixed nonzero and exact-zero histories."""

    ctx = FixtureContext()
    logs = [ctx.mpf(Fraction(7, 17)), ctx.mpf(Fraction(12, 17)),
            ctx.mpf(Fraction(15, 17))]
    ledger: dict[str, Any] = {"fixture_v3_ledger": "fresh_source_ledger_alpha"}
    flow = FixtureFlow(ctx, logs, ledger)
    op = FixtureOp(ctx, flow,
                   tuple(_row(flow, Fraction(13, 29)) for _ in range(6)))

    V = [_row(flow, value) for value in
         (Fraction(-4, 11), Fraction(9, 14), Fraction(13, 25),
          Fraction(-7, 27))]
    # Three rows are nonzero history; one exact zero row must not be treated as
    # evidence that the source-owned history was erased.
    m = [_row(flow, value) for value in
         (Fraction(2, 29), Fraction(0), Fraction(-3, 31), Fraction(5, 37))]
    # V_y is raw physical U_y and must not receive a second Pstar factor.
    Vy = [_row(flow, value) for value in
          (Fraction(-11, 16), Fraction(7, 22), Fraction(-13, 28),
           Fraction(17, 32))]
    R = _row(flow, Fraction(14, 9))
    Pstar = _row(flow, Fraction(9, 4))
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
        "candidate_N": 1729,
        "chart": "first_micro",
        "exact_Z_cell": ["-5/16", "7/19"],
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
        "errors": {"v3_packet_lineage_error": {"kind": "opaque", "token": "alpha"}},
    }
    view = {
        "V": V, "m": m, "R": R, "Pstar": Pstar, "source_packet": packet,
        "errors": {"v3_source_rows_lineage_error": {"kind": "opaque", "token": "beta"}},
    }
    Z = ctx.mpf((Fraction(3, 14), Fraction(5, 18)))
    delta = ctx.mpf(Fraction(4, 13))
    L = ctx.mpf(1) - delta * Z * Z
    assumptions = {
        "query": {
            "source_packet": packet, "Z": Z,
            "source": "synthetic_v3_query_receipt:lower_micro/sc_half",
            "checks": {"packet_query_exact": True, "coordinate_exact": True,
                       "candidate_N_is_source_owned": True,
                       "Z_matches_exact_Z_cell": True},
            "details": {"chart": "first_micro", "exact_left": ["sc_half"],
                        "exact_right": ["sc_half"], "candidate_N": 1729,
                        "exact_Z_cell": ["-5/16", "7/19"]},
        },
        "z_independence": {
            "source_packet": packet, "Z": Z, "delta": delta,
            "source": "synthetic_v3_geometry_receipt:caller_conditional",
            "checks": {"R_Z_independent": True, "Pstar_Z_independent": True,
                       "delta_Z_independent": True, "full_target_Z_interval": True},
        },
        "history": {
            "source_packet": packet,
            "source": "synthetic_v3_history_receipt:mixed_rows",
            "checks": {"packet_m_rows_exact": True,
                       "original_history_not_zeroed": True,
                       "source_owned_flat_collar": True},
        },
        "L_positive": {
            "source_packet": packet, "Z": Z, "delta": delta,
            "source": "synthetic_v3_interval_receipt:positive_L",
            "checks": {"L_strictly_positive": True}, "lower_bound": L.lo,
        },
    }
    return {"view": view, "op": op, "ctx": ctx, "Z": Z, "delta": delta,
            "assumptions": assumptions, "flow": flow, "packet": packet,
            "V": V, "m": m, "Vy": Vy, "R": R, "Pstar": Pstar}


def build_synthetic_inputs() -> dict[str, Any]:
    """Return fresh synthetic inputs for a caller-owned evaluator invocation."""

    return _build_inputs()


build_fixture_inputs = build_synthetic_inputs


def _expected(inputs: Mapping[str, Any], *, duplicate_scale: bool = False
              ) -> tuple[ExactInterval, ExactInterval]:
    """Compute independent exact RQy and ordinary first-Z expressions."""

    view = inputs["view"]
    P, R = view["Pstar"].coefficient, view["R"].coefficient
    V, m = view["V"], view["m"]
    vy = view["source_packet"]["actual_background_source"]["V_y"]
    Z, delta = inputs["Z"], inputs["delta"]
    U0, U1, U2 = P * V[0].coefficient, P * V[1].coefficient, P * V[2].coefficient * 2
    M0, M1, M2 = P * m[0].coefficient, P * m[1].coefficient, P * m[2].coefficient * 2
    if duplicate_scale:
        Uy0, Uy1 = P * vy[0].coefficient, P * vy[1].coefficient
    else:
        Uy0, Uy1 = vy[0].coefficient, vy[1].coefficient
    A, Az, Azz = U0 - M0, U1 - M1, U2 - M2
    one = ExactInterval(1, ctx=inputs["ctx"])
    Z2, one_minus_Z2 = Z * Z, one - Z * Z
    one_minus_delta = one - delta
    L, inverse_L = one - delta * Z2, one / (one - delta * Z2)
    two_Z, two_delta_Z = Z * 2, delta * Z * 2
    Ny = Uy0 * two_Z - A * (one_minus_delta * Z) - Az * one_minus_Z2
    NyZ = (Uy0 * 2 + Uy1 * two_Z - (A + Az * Z) * one_minus_delta
           + Az * two_Z - Azz * one_minus_Z2)
    RQy = R * Ny * inverse_L
    dRQy = R * (NyZ * L + Ny * two_delta_Z) * inverse_L * inverse_L
    return RQy, dRQy


def _copy_with_packet(inputs: Mapping[str, Any], packet: Mapping[str, Any],
                      **view_changes: Any) -> dict[str, Any]:
    """Copy inputs while rebinding every certificate to the replacement packet."""

    old_packet = inputs["packet"]
    view = dict(inputs["view"], **view_changes, source_packet=packet)
    assumptions = {}
    for name, cert in inputs["assumptions"].items():
        if isinstance(cert, Mapping) and cert.get("source_packet") is old_packet:
            cert = dict(cert)
            cert["source_packet"] = packet
        assumptions[name] = cert
    result = dict(inputs, view=view, assumptions=assumptions, packet=packet)
    for name in ("V", "m", "R", "Pstar", "Vy"):
        if name in view_changes:
            result[name] = view_changes[name]
    return result


def _all_zero_history_inputs() -> dict[str, Any]:
    inputs = _build_inputs()
    zero_m = [inputs["flow"].scalar(0) for _ in inputs["m"]]
    packet = dict(inputs["packet"])
    generic = dict(packet["original_generic_source"])
    generic["common_own_five_histories_axial5"] = {"m": zero_m}
    packet["original_generic_source"] = generic
    return _copy_with_packet(inputs, packet, m=zero_m)


def _invoke(evaluate: Callable[..., Any], inputs: Mapping[str, Any]) -> Any:
    return evaluate(inputs["view"], inputs["op"], inputs["ctx"], inputs["Z"],
                    inputs["delta"], assumptions=inputs["assumptions"])


def _different_pstar_inputs() -> dict[str, Any]:
    inputs = _build_inputs()
    Pstar = inputs["flow"].scalar(Fraction(13, 5))
    inputs["Pstar"] = Pstar
    inputs["op"].Pstar = Pstar
    inputs["view"] = dict(inputs["view"], Pstar=Pstar)
    return inputs


def run_fixture(evaluate: Callable[..., Any]) -> dict[str, Any]:
    """Run exact positive checks and fail-closed negative controls."""

    if not callable(evaluate):
        raise TypeError("evaluate must be callable")
    inputs = _build_inputs()
    checks: dict[str, bool] = {}
    diagnostics: dict[str, Any] = {
        "schema": SCHEMA, "source_family": SOURCE_FAMILY,
        "synthetic_only": True, "real_source_provider_called": False,
        "actual_source_validation_claimed": False, "candidate_N": 1729,
        "ordinary_Z_domain": ["3/14", "5/18"], "delta": "4/13", "Pstar": "9/4",
    }
    positive_names = (
        "returned_mapping", "conditional_status", "returned_source_enclosures",
        "ordinary_Z1_retained", "same_basis_and_ledger", "exact_RQy_formula",
        "exact_ordinary_Z1_formula", "physical_Uy_not_rescaled",
        "duplicate_scale_not_present", "nonzero_m_history_preserved",
        "zero_m_row_condition_retained", "all_zero_m_condition_retained",
        "full_target_Z_interval_used", "error_layers_reported_truthfully",
        "errors_preserved_when_both_present", "physical_expression_rows_exact",
        "independence_marked_conditional", "different_Pstar_keeps_physical_Uy",
        "independent_polynomial_Z1_difference", "known_duplicate_scale_defect_Z0_Z1",
    )

    try:
        result = _invoke(evaluate, inputs)
        checks["returned_mapping"] = isinstance(result, Mapping)
        checks["conditional_status"] = result.get("status") == "conditional_source_enclosure"
        source = result.get("source_enclosures")
        uy_rows = source.get("U_y") if isinstance(source, Mapping) else None
        # U_y is a pair: validate both elements, not the pair as one enclosure.
        checks["returned_source_enclosures"] = (
            isinstance(source, Mapping)
            and all(isinstance(source.get(name), ScaledEnclosure)
                    for name in ("R_Q_y", "dZ_R_Q_y", "ordinary_Z1"))
            and isinstance(uy_rows, (tuple, list)) and len(uy_rows) == 2
            and all(isinstance(row, ScaledEnclosure) for row in uy_rows)
        )
        checks["ordinary_Z1_retained"] = (
            isinstance(source, Mapping)
            and source.get("ordinary_Z1") is result.get("ordinary_Z1")
        )
        checks["same_basis_and_ledger"] = all(
            row.scale.bases is inputs["flow"].logs
            and row.ledger is inputs["flow"].ledger
            for row in (source["R_Q_y"], source["dZ_R_Q_y"], uy_rows[0], uy_rows[1])
        )
        expected_RQy, expected_Z1 = _expected(inputs)
        wrong_RQy, wrong_Z1 = _expected(inputs, duplicate_scale=True)
        checks["exact_RQy_formula"] = source["R_Q_y"].coefficient == expected_RQy
        checks["exact_ordinary_Z1_formula"] = source["ordinary_Z1"].coefficient == expected_Z1
        checks["physical_Uy_not_rescaled"] = (
            uy_rows[0] is inputs["Vy"][0] and uy_rows[1] is inputs["Vy"][1]
        )
        checks["duplicate_scale_not_present"] = (
            source["R_Q_y"].coefficient != wrong_RQy
            and source["ordinary_Z1"].coefficient != wrong_Z1
        )
        checks["nonzero_m_history_preserved"] = (
            any(not row.zero for row in inputs["m"])
            and result.get("proof", {}).get("original_m_history_not_zeroed") is True
        )
        checks["zero_m_row_condition_retained"] = (
            any(row.zero for row in inputs["m"])
            and result.get("proof", {}).get("zero_m_rows_not_used_as_history_evidence") is True
        )
        zeroed_result = _invoke(evaluate, _all_zero_history_inputs())
        checks["all_zero_m_condition_retained"] = (
            zeroed_result.get("proof", {}).get("original_m_history_not_zeroed") is True
            and zeroed_result.get("proof", {}).get(
                "zero_m_rows_not_used_as_history_evidence") is True
        )
        checks["full_target_Z_interval_used"] = result.get("scope", {}).get(
            "target_Z_is_full_requested_interval") is True
        errors = result.get("scope", {}).get("error_layers", {})
        checks["error_layers_reported_truthfully"] = (
            errors.get("source_rows_present") is True
            and errors.get("source_packet_present") is True
            and errors.get("both_layers_present") is True
            and result.get("scope", {}).get("retained_error_terms", {}).get(
                "source_rows") is inputs["view"]["errors"]
            and result.get("scope", {}).get("retained_error_terms", {}).get(
                "source_packet") is inputs["packet"]["errors"]
        )
        checks["errors_preserved_when_both_present"] = (
            result.get("scope", {}).get("declared_error_layers_retained") is True
            and result.get("scope", {}).get(
                "all_inherited_errors_verified_preserved") is False
        )
        checks["physical_expression_rows_exact"] = (
            uy_rows[0].coefficient == inputs["Vy"][0].coefficient
            and uy_rows[1].coefficient == inputs["Vy"][1].coefficient
        )
        checks["independence_marked_conditional"] = (
            result.get("proof", {}).get("R_Pstar_delta_Z_independence")
            == "conditional_on_caller_certificate"
            and result.get("scope", {}).get("R_Pstar_delta_Z_independence")
            == "conditional_source_certificate"
        )

        other = _different_pstar_inputs()
        other_result = _invoke(evaluate, other)
        checks["different_Pstar_keeps_physical_Uy"] = (
            other["Pstar"].coefficient != inputs["Pstar"].coefficient
            and other_result["source_enclosures"]["U_y"][0].coefficient
            == other["Vy"][0].coefficient
            and other_result["source_enclosures"]["U_y"][1].coefficient
            == other["Vy"][1].coefficient
        )

        # Independent finite difference of the quadratic Taylor polynomials.
        point = _build_inputs()
        point_z = point["ctx"].mpf(Fraction(5, 21))
        point["Z"] = point_z
        for key in ("query", "z_independence", "L_positive"):
            point["assumptions"][key]["Z"] = point_z
        point_L = point["ctx"].mpf(1) - point["delta"] * point_z * point_z
        point["assumptions"]["L_positive"]["lower_bound"] = point_L.lo
        point_result = _invoke(evaluate, point)
        P, radius, d = (point["Pstar"].coefficient.lo,
                        point["R"].coefficient.lo, point["delta"].lo)
        u0, u1, u2 = (P * point["V"][k].coefficient.lo * (2 if k == 2 else 1)
                      for k in range(3))
        m0, m1, m2 = (P * point["m"][k].coefficient.lo * (2 if k == 2 else 1)
                      for k in range(3))
        y0, y1 = point["Vy"][0].coefficient.lo, point["Vy"][1].coefficient.lo
        z0, h = Fraction(5, 21), Fraction(1, 200000)

        def analytic_flux(offset: Fraction) -> Fraction:
            z = z0 + offset
            u = u0 + u1 * offset + u2 * offset * offset / 2
            m_value = m0 + m1 * offset + m2 * offset * offset / 2
            uz, mz, uy = u1 + u2 * offset, m1 + m2 * offset, y0 + y1 * offset
            numerator = (2 * z * uy - (1 - d) * z * (u - m_value)
                         - (1 - z * z) * (uz - mz))
            return radius * numerator / (1 - d * z * z)

        finite_difference = (analytic_flux(h) - analytic_flux(-h)) / (2 * h)
        derivative = point_result["ordinary_Z1"].coefficient
        checks["independent_polynomial_Z1_difference"] = (
            derivative.lo == derivative.hi
            and abs(derivative.lo - finite_difference) < Fraction(1, 7000000)
        )

        # Exact old double-scale defect, including both Z0 and ordinary Z1.
        legacy_flux, legacy_derivative = _expected(point, duplicate_scale=True)
        physical_flux, physical_derivative = point_result["R_Q_y"].coefficient.lo, derivative.lo
        denominator = 1 - d * z0 * z0
        defect0 = radius * 2 * z0 * (P - 1) * y0 / denominator
        defect1 = radius * 2 * (P - 1) * (
            (y0 + z0 * y1) / denominator
            + 2 * d * z0 * z0 * y0 / denominator**2
        )
        checks["known_duplicate_scale_defect_Z0_Z1"] = (
            legacy_flux.lo - physical_flux == defect0
            and legacy_derivative.lo - physical_derivative == defect1
        )
    except Exception as exc:
        diagnostics["positive_error"] = type(exc).__name__ + ": " + str(exc)
        for name in positive_names:
            checks[name] = False

    wrong_basis = _build_inputs()
    other_logs = [wrong_basis["ctx"].mpf(Fraction(4, 7))]
    wrong_basis["view"] = dict(wrong_basis["view"])
    wrong_V = list(wrong_basis["V"])
    wrong_V[0] = ScaledEnclosure(FormalScale(other_logs), Fraction(-4, 11),
                                 wrong_basis["flow"].ledger)
    wrong_basis["view"]["V"] = wrong_V
    try:
        _invoke(evaluate, wrong_basis)
    except Exception:
        checks["wrong_basis_rejected"] = True
    else:
        checks["wrong_basis_rejected"] = False

    # Rebind certificates so this rejection is caused by the malformed query.
    wrong_query = _build_inputs()
    bad_packet = dict(wrong_query["packet"])
    bad_packet["exact_left"] = ["v3_wrong_lower_query"]
    wrong_query = _copy_with_packet(wrong_query, bad_packet)
    try:
        _invoke(evaluate, wrong_query)
    except Exception:
        checks["wrong_query_rejected"] = True
    else:
        checks["wrong_query_rejected"] = False

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

    diagnostics.update({
        "m_history_rows": 4, "m_nonzero_rows": 3, "m_zero_rows": 1,
        "all_zero_variant_rows": 4, "checks_count": len(checks),
        "physical_Uy_expression": "Uy0=V_y[0]; Uy1=V_y[1]",
        "duplicate_scale_expression": "Pstar*V_y[0], Pstar*V_y[1] (rejected)",
    })
    return {"passed": all(checks.values()), "checks": checks,
            "diagnostics": diagnostics}


__all__ = [
    "SCHEMA", "SOURCE_FAMILY", "ExactInterval", "FormalScale", "ScaledEnclosure",
    "FixtureContext", "FixtureFlow", "FixtureOp", "build_synthetic_inputs",
    "build_fixture_inputs", "run_fixture",
]


