"""Pure post-restore ``Rm`` endpoint arithmetic for the mixed physical-y packet.

This adapter is deliberately a byte-authenticated data boundary.  It consumes
the new mixed-native and mixed-loop packets and the two reviewed report byte
strings, but it never imports or executes a historical producer, builder,
owner, checker, recurrence, quadrature, norm, or endpoint helper.

The native ``Qy`` rows are signed source intervals and are copied from the
post-restore chart.  Finite-``N`` and mixed-loop rows are logarithmic absolute
majorants.  A majorant is never used as a signed function value, and the
normalization keeps the full directed ``log_Rm`` interval without materializing
the astronomical radius.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path
from typing import Any

from mpmath import isfinite, make_mpf, nstr
from mpmath.ctx_iv import MPIntervalContext


NATIVE_REPORT_PATH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json"
)
NATIVE_REPORT_BYTES = 42177499
NATIVE_REPORT_SHA256 = (
    "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"
)
FINITE_REPORT_PATH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_directed_finite_N_controls.json"
)
FINITE_REPORT_BYTES = 1368720
FINITE_REPORT_SHA256 = (
    "152d211cd137e7e0498fd157a0ac1467541dc6214ae59c13af62844112003d9a"
)

CHART = "post_restore"
ORDERS = (0, 1, 2, 3)
ENDPOINT_ORDERS = (0, 1)
Q_ORDERS = (0, 1, 2)
MIXED_INDICES = tuple((i, k) for i in (0, 1) for k in ORDERS)
EXTERNAL_RQY_HANDLES = (183696, 183697)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _lo(value: Any):
    return make_mpf(tuple(value._mpi_[0]))


def _hi(value: Any):
    return make_mpf(tuple(value._mpi_[1]))


def _canonical_endpoint(value: Any, label: str) -> tuple[int, int, int, int]:
    _require(isinstance(value, (list, tuple)) and len(value) == 4,
             label + " must be an exact MPF tuple")
    _require(all(type(part) is int for part in value),
             label + " MPF tuple components must be integers")
    sign, mantissa, exponent, bit_count = tuple(value)
    _require(sign in (0, 1) and mantissa >= 0,
             label + " MPF sign or mantissa is malformed")
    _require(bit_count == mantissa.bit_length(),
             label + " MPF bit count is not canonical")
    if mantissa == 0:
        _require((sign, mantissa, exponent, bit_count) == (0, 0, 0, 0),
                 label + " has a noncanonical zero")
    else:
        _require(mantissa & 1, label + " MPF mantissa is not normalized")
    return sign, mantissa, exponent, bit_count


def _record(value: Any, label: str = "interval") -> dict[str, Any]:
    """Serialize both endpoints of a finite directed interval."""

    _require(hasattr(value, "_mpi_"), label + " must be an interval")
    lower = tuple(value._mpi_[0])
    upper = tuple(value._mpi_[1])
    _canonical_endpoint(lower, label + ".lower")
    _canonical_endpoint(upper, label + ".upper")
    lo, hi = make_mpf(lower), make_mpf(upper)
    _require(lo <= hi, label + " has reversed endpoints")
    _require(isfinite(lo) and isfinite(hi), label + " must be finite")
    return {
        "lower_exact_mpf_tuple": list(lower),
        "upper_exact_mpf_tuple": list(upper),
        "lower_display": nstr(lo, 24),
        "upper_display": nstr(hi, 24),
        "directed": True,
    }


def _read_interval(c: MPIntervalContext, value: Any,
                   label: str = "directed interval"):
    """Decode an exact directed record without selecting an endpoint."""

    if hasattr(value, "_mpi_"):
        out = value
    else:
        if isinstance(value, Mapping) and "value" in value \
                and "lower_exact_mpf_tuple" not in value:
            value = value["value"]
        _require(isinstance(value, Mapping), label + " must be a mapping")
        _require(value.get("directed") is True,
                 label + " must carry directed=True")
        lower = _canonical_endpoint(value.get("lower_exact_mpf_tuple"),
                                    label + ".lower")
        upper = _canonical_endpoint(value.get("upper_exact_mpf_tuple"),
                                    label + ".upper")
        _require(make_mpf(lower) <= make_mpf(upper),
                 label + " has reversed endpoints")
        out = c.mpf([make_mpf(lower), make_mpf(upper)])
    _require(hasattr(out, "_mpi_"), label + " must be an interval")
    lo, hi = _lo(out), _hi(out)
    _require(lo <= hi and isfinite(lo) and isfinite(hi),
             label + " must be finite and ordered")
    return c.mpf([lo, hi])


def _same_interval(left: Any, right: Any) -> bool:
    return (tuple(left._mpi_[0]) == tuple(right._mpi_[0]) and
            tuple(left._mpi_[1]) == tuple(right._mpi_[1]))


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _admit_report(raw: bytes, expected_bytes: int, expected_sha256: str,
                  label: str) -> Mapping[str, Any]:
    _require(type(raw) is bytes, label + " must be immutable bytes")
    _require(len(raw) == expected_bytes,
             label + " byte length changed")
    _require(_sha(raw) == expected_sha256,
             label + " SHA256 changed")
    try:
        report = json.loads(raw)
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError(label + " is not valid JSON") from exc
    _require(isinstance(report, Mapping), label + " must be a report mapping")
    return report


def _at(mapping: Mapping[str, Any], *keys: str) -> Any:
    value: Any = mapping
    for key in keys:
        _require(isinstance(value, Mapping) and key in value,
                 "Missing report path: " + ".".join(keys))
        value = value[key]
    return value


def _packet_chart(packet: Mapping[str, Any]) -> Mapping[str, Any]:
    charts = packet.get("charts")
    if isinstance(charts, Mapping):
        for key in (CHART, "postrestore"):
            candidate = charts.get(key)
            if isinstance(candidate, Mapping):
                return candidate
    for key in (CHART, "postrestore", "post_restore_chart"):
        candidate = packet.get(key)
        if isinstance(candidate, Mapping):
            return candidate
    return packet


def _validate_packets(native_packet: Any, loop_packet: Any,
                      native_sha: str) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    _require(isinstance(native_packet, Mapping),
             "mixed_native_packet must be a mapping")
    _require(isinstance(loop_packet, Mapping),
             "mixed_loop_packet must be a mapping")
    for packet, label in ((native_packet, "mixed-native"),
                          (loop_packet, "mixed-loop")):
        _require(packet.get("passed") is True, label + " packet did not pass")
        for key in ("source_executed", "accepted_ancestors_replayed",
                    "accepted_builders_or_owners_replayed",
                    "accepted_checkers_replayed", "recurrences_replayed",
                    "quadratures_replayed", "norms_replayed"):
            if key in packet:
                _require(packet[key] is False, label + " reports forbidden execution: " + key)
    for owner, label in ((native_packet, "mixed-native"),
                         (loop_packet, "mixed-loop")):
        for key in ("source_native_report_sha256", "native_parent_report_sha256",
                    "source_native_parent_report_sha256"):
            if key in owner:
                _require(owner[key] == native_sha,
                         label + " native parent SHA256 changed")
    for packet, label in ((native_packet, "mixed-native"),
                          (loop_packet, "mixed-loop")):
        chart = packet.get("chart")
        if chart is not None:
            _require(chart == CHART, label + " is not the post_restore chart")
    loop_chart = _packet_chart(loop_packet)
    native_chart = _packet_chart(native_packet)
    return native_chart, loop_chart


def _row_value(c: MPIntervalContext, row: Mapping[str, Any], label: str):
    exact = row.get("exact_zero", False)
    _require(type(exact) is bool, label + " exact_zero must be boolean")
    raw = row.get("log_absolute_upper_bound")
    if raw is None:
        raw = row.get("record", row.get("log_cap"))
    if exact:
        _require(raw is None, label + " exact-zero row carries a cap")
        return None
    _require(raw is not None, label + " nonzero row has no directed record")
    return _read_interval(c, raw, label + ".log_absolute_upper_bound")


def _mixed_rows(c: MPIntervalContext, raw: Any, label: str) -> dict[tuple[int, int], Any]:
    """Decode exactly the canonical eight physical-y/Z rows."""

    if isinstance(raw, Mapping):
        for key in ("rows", "mixed_log_caps", "log_caps", "values"):
            if key in raw:
                return _mixed_rows(c, raw[key], label)
        if set(raw) >= {"y0", "y1"}:
            _require(isinstance(raw["y0"], (list, tuple)) and
                     isinstance(raw["y1"], (list, tuple)),
                     label + " y0/y1 rows must be lists")
            raw = list(raw["y0"]) + list(raw["y1"])
        else:
            # An in-memory mapping may use tuple indices.
            try:
                if set(raw) == set(MIXED_INDICES):
                    return {
                        index: (None if raw[index] is None else
                                _read_interval(c, raw[index], label + str(index)))
                        for index in MIXED_INDICES
                    }
            except TypeError:
                pass
    _require(isinstance(raw, (list, tuple)) and len(raw) == 8,
             label + " must contain exactly eight mixed rows")
    out: dict[tuple[int, int], Any] = {}
    for number, row in enumerate(raw):
        _require(isinstance(row, Mapping), label + " row must be a mapping")
        i = row.get("physical_y_order", row.get("y_order",
            row.get("physical_y", row.get("i"))))
        k = row.get("ordinary_Z_order", row.get("Z_order", row.get("k")))
        _require(type(i) is int and type(k) is int and (i, k) in MIXED_INDICES,
                 label + " has an invalid y/Z row index")
        _require((i, k) not in out,
                 label + " repeats y/Z row " + str((i, k)))
        out[i, k] = _row_value(c, row, label + "[" + str(number) + "]")
    _require(set(out) == set(MIXED_INDICES),
             label + " must cover y0/y1 and Z0..Z3")
    return out


def _packet_rows(c: MPIntervalContext, packet: Mapping[str, Any],
                 chart: Mapping[str, Any], name: str) -> dict[tuple[int, int], Any]:
    raw = None
    for owner in (packet, chart):
        if not isinstance(owner, Mapping):
            continue
        if name in owner:
            raw = owner[name]
            break
        for container in ("mixed_log_caps", "mixed_caps", "log_caps"):
            candidate = owner.get(container)
            if isinstance(candidate, Mapping) and name in candidate:
                raw = candidate[name]
                break
        if raw is not None:
            break
    _require(raw is not None, "Mixed-loop packet is missing " + name)
    return _mixed_rows(c, raw, name)


def _packet_scalar(c: MPIntervalContext, packet: Mapping[str, Any],
                  chart: Mapping[str, Any], names: tuple[str, ...]):
    owners = [packet, chart]
    for owner in (packet, chart):
        if isinstance(owner, Mapping):
            for key in ("canonical_parent_scale", "parent_scale", "scale"):
                candidate = owner.get(key)
                if isinstance(candidate, Mapping):
                    owners.append(candidate)
    for owner in owners:
        if not isinstance(owner, Mapping):
            continue
        for name in names:
            if name in owner:
                return _read_interval(c, owner[name], "mixed-loop " + name)
    return None


def _saved_log_rows(c: MPIntervalContext, raw: Any, count: int,
                    label: str) -> dict[int, Any]:
    _require(isinstance(raw, (list, tuple)), label + " must be a row list")
    out: dict[int, Any] = {}
    for number, row in enumerate(raw):
        _require(isinstance(row, Mapping), label + " row must be a mapping")
        k = row.get("ordinary_Z_order", row.get("Z_order", number))
        _require(type(k) is int and k in range(count),
                 label + " has an invalid ordinary Z order")
        _require(k not in out, label + " repeats an ordinary Z order")
        out[k] = _row_value(c, row, label + "[" + str(number) + "]")
    _require(set(out) == set(range(count)),
             label + " must cover ordinary Z orders 0 through " + str(count - 1))
    return out


def _logsum(c: MPIntervalContext, terms: Any):
    terms = tuple(term for term in terms if term is not None)
    if not terms:
        return None
    if len(terms) == 1:
        return terms[0]
    shift = c.ln(c.mpf(len(terms)))
    shifted = [term + shift for term in terms]
    return c.mpf([min(_lo(term) for term in shifted),
                  max(_hi(term) for term in shifted)])


def _term(*parts: Any):
    if any(part is None for part in parts):
        return None
    result = parts[0]
    for part in parts[1:]:
        result = result + part
    return result


def _positive_log(c: MPIntervalContext, value: Any, label: str):
    _require(_lo(value) > 0, label + " must have a positive lower bound")
    return c.ln(value)


def _absolute_log(c: MPIntervalContext, value: Any, label: str):
    """Bound a signed interval's magnitude without asserting a sign."""

    lo, hi = _lo(value), _hi(value)
    if lo == 0 and hi == 0:
        return None
    if lo <= 0 <= hi:
        magnitude = max(-lo, hi)
        return c.ln(c.mpf([magnitude, magnitude]))
    if lo > 0:
        return c.ln(c.mpf([lo, hi]))
    return c.ln(c.mpf([-hi, -lo]))


def _pack_log(k: int, value: Any) -> dict[str, Any]:
    return {
        "ordinary_Z_order": k,
        "exact_zero": value is None,
        "log_absolute_upper_bound": None if value is None else _record(value),
    }


def _pack_log_map(values: Mapping[int, Any], orders: tuple[int, ...]) -> list[dict[str, Any]]:
    return [_pack_log(k, values[k]) for k in orders]


def _signed_log_scale(c: MPIntervalContext, signed_value: Any,
                      log_rm: Any, label: str) -> dict[str, Any]:
    """Represent ``signed_value / Rm`` by sign and a directed log scale."""

    lo, hi = _lo(signed_value), _hi(signed_value)
    absolute = _absolute_log(c, signed_value, label)
    if absolute is None:
        return {
            "exact_zero": True,
            "contains_zero": True,
            "sign": None,
            "log_absolute_upper_bound": None,
            "signed_source_interval": _record(signed_value, label),
        }
    scaled = absolute - log_rm
    if lo <= 0 <= hi:
        sign = None
        contains_zero = True
    elif lo > 0:
        sign = 1
        contains_zero = False
    else:
        sign = -1
        contains_zero = False
    return {
        "exact_zero": False,
        "contains_zero": contains_zero,
        "sign": sign,
        "log_absolute_upper_bound": _record(scaled, label + ".log_scale"),
        "signed_source_interval": _record(signed_value, label),
        "formula": "log_abs(Qy/Rm)=log_abs(Qy)-log_Rm",
    }


def _finite_sources(c: MPIntervalContext, report: Mapping[str, Any]):
    bounds = _at(report, "actual_directed_finite_N_controls", "charts",
                  CHART, "physical_error_log_bounds")
    delta_u = _saved_log_rows(c, bounds["delta_U"], 4, "finite delta_U")
    delta_m = _saved_log_rows(c, bounds["delta_M"], 4, "finite delta_M")
    delta_q = _saved_log_rows(c, bounds["delta_Q"], 3, "finite delta_Q")
    log_n = _read_interval(c, bounds["logN"], "finite logN")
    log_pe = _read_interval(c, bounds["log_Pstar_E_upper"],
                            "finite log_Pstar_E_upper")
    return bounds, delta_u, delta_m, delta_q, log_n, log_pe


def _native_sources(c: MPIntervalContext, report: Mapping[str, Any]):
    repaired = _at(report, "actual_restoration_post_jet_repair",
                    "corrected_background")
    chart = _at(repaired, "charts", CHART)
    raw_qy = _at(chart, "endpoint_fields", "Qy_end")
    _require(isinstance(raw_qy, (list, tuple)) and len(raw_qy) >= 2,
             "Native post_restore Qy_end has fewer than two rows")
    qy: dict[int, Any] = {}
    source_rows = []
    for index, row in enumerate(raw_qy[:2]):
        _require(isinstance(row, Mapping), "Native Qy_end row is not a mapping")
        _require(row.get("ordinary_Z_order") == index,
                 "Native Qy_end[:2] ordinary-Z order changed")
        qy[index] = _read_interval(c, row.get("value"),
                                   "native Qy_end[" + str(index) + "]")
        source_rows.append(copy.deepcopy(row))
    background = _at(report, "actual_restoration_background_I1_numeric")
    source_scale = _read_interval(
        c, _at(background, "physical_scale", "log_Rm"),
        "actual_restoration_background_I1_numeric.physical_scale.log_Rm",
    )
    background_chart_scale = _read_interval(
        c, _at(background, "charts", CHART, "log_R_upper"),
        "actual_restoration_background_I1_numeric.charts.post_restore.log_R_upper",
    )
    corrected_scale = _read_interval(
        c, _at(repaired, "physical_scale", "log_Rm"),
        "corrected_background.physical_scale.log_Rm",
    )
    chart_scale = _read_interval(c, _at(chart, "log_R_upper"),
                                 "corrected_background.charts.post_restore.log_R_upper")
    _require(_same_interval(source_scale, background_chart_scale) and
             _same_interval(source_scale, corrected_scale) and
             _same_interval(source_scale, chart_scale),
             "Native log_Rm source/chart directed records disagree")
    delta = _read_interval(
        c, _at(report, "actual_core_signed_direction", "parameters", "delta"),
        "native actual_core_signed_direction.parameters.delta",
    )
    domain = _at(background, "strict_Z_domain")
    corrected_domain = _at(repaired, "strict_Z_domain")
    finite_domain = _at(report, "actual_restoration_finite_N_I1_numeric",
                         "strict_Z_domain")
    z_lower = _read_interval(c, domain["lower"], "native strict Z lower")
    z_upper = _read_interval(c, domain["upper"], "native strict Z upper")
    for source, source_label in (
        (corrected_domain, "corrected_background.strict_Z_domain"),
        (finite_domain, "actual_restoration_finite_N_I1_numeric.strict_Z_domain"),
    ):
        other_lower = _read_interval(c, source["lower"], source_label + ".lower")
        other_upper = _read_interval(c, source["upper"], source_label + ".upper")
        _require(_same_interval(z_lower, other_lower) and
                 _same_interval(z_upper, other_upper),
                 "Native strict Z domain source records disagree: " + source_label)
    z = c.mpf([_lo(z_lower), _hi(z_upper)])
    _require(_lo(z) > 0, "Strict Z domain must be positive")
    return {
        "qy": qy,
        "source_rows": source_rows,
        "log_rm": source_scale,
        "background_chart_log_rm": background_chart_scale,
        "corrected_scale_log_rm": corrected_scale,
        "chart_log_rm": chart_scale,
        "delta": delta,
        "z": z,
        "z_lower": z_lower,
        "z_upper": z_upper,
    }


def _finite_delta_uy(c: MPIntervalContext, slow: Mapping[tuple[int, int], Any],
                     fast: Mapping[tuple[int, int], Any], log_pe: Any,
                     log_n: Any):
    slow_terms: dict[int, Any] = {}
    fast_terms: dict[int, Any] = {}
    total: dict[int, Any] = {}
    for k in ORDERS:
        slow_terms[k] = _term(log_pe, slow[1, k], -log_n)
        fast_terms[k] = _term(log_pe, fast[0, k])
        total[k] = _logsum(c, (slow_terms[k], fast_terms[k]))
    return total, slow_terms, fast_terms


def _finite_delta_my(c: MPIntervalContext, delta_u: Mapping[int, Any],
                     delta_m: Mapping[int, Any]):
    # The source identity is signed deltaMy = deltaU - deltaM.  The emitted
    # rows are absolute caps, so cancellation cannot be claimed or selected.
    return {k: _logsum(c, (delta_u[k], delta_m[k])) for k in ORDERS}


def _positive_pressure_coefficients(c: MPIntervalContext, delta: Any,
                                    z: Any):
    one = c.mpf(1)
    h = one - delta
    d = one - z * z
    L = one - delta * z * z
    _require(_lo(h) > 0, "h=1-delta has no positive lower bound")
    _require(_lo(d) > 0, "d=1-Z^2 has no positive lower bound")
    _require(_lo(L) > 0, "L=1-delta*Z^2 has no positive lower bound")
    L_floor = c.mpf([_lo(L), _lo(L)])
    return h, d, L, L_floor


def _finite_delta_qy(c: MPIntervalContext, delta_uy: Mapping[int, Any],
                     delta_my: Mapping[int, Any], delta: Any, z: Any):
    _require(_lo(delta) > 0 and _hi(delta) < 1,
             "delta must have a directed enclosure strictly inside (0,1)")
    h, d, L, L_floor = _positive_pressure_coefficients(c, delta, z)
    log_two = c.ln(c.mpf(2))
    log_z = _positive_log(c, z, "Z")
    log_h = _positive_log(c, h, "h")
    log_d = _positive_log(c, d, "d")
    log_one_plus_delta = _positive_log(c, c.mpf(1) + delta,
                                       "1+delta")
    # Reciprocal uses ln of the positive lower endpoint, never abs_log(L).
    log_L_lower = c.ln(L_floor)
    q0_terms = (
        _term(log_two, log_z, delta_uy[0]),
        _term(log_h, log_z, delta_my[0]),
        _term(log_d, delta_my[1]),
    )
    q1_terms = (
        _term(log_two, delta_uy[0]),
        _term(log_two, log_z, delta_uy[1]),
        _term(log_h, delta_my[0]),
        _term(log_one_plus_delta, log_z, delta_my[1]),
        _term(log_d, delta_my[2]),
    )
    n0 = _logsum(c, q0_terms)
    n1 = _logsum(c, q1_terms)
    q0 = None if n0 is None else n0 - log_L_lower
    # q1 is the derivative of N/L: N1/L + N0*(2*delta*Z)/L^2.
    # The second term is retained even though all rows are absolute caps.
    log_delta = _positive_log(c, delta, "delta")
    denominator_correction = (
        None if n0 is None else
        n0 + log_two + log_delta + log_z - 2 * log_L_lower
    )
    q1 = _logsum(
        c,
        (None if n1 is None else n1 - log_L_lower,
         denominator_correction),
    )
    qy = {0: q0, 1: q1}
    receipt = {
        "formula": (
            "deltaQy=(2Z*deltaUy-(1-delta)*Z*deltaMy-"
            "(1-Z^2)*(deltaMy)_Z)/(1-delta*Z^2)"
        ),
        "C1_formula": (
            "C1=2*deltaUy+2*Z*deltaUy_Z-(1-delta)*deltaMy+"
            "(1+delta)*Z*deltaMy_Z-(1-Z^2)*deltaMy_ZZ"
        ),
        "C1_is_numerator_only": True,
        "positive_L_lower_used": True,
        "L": _record(L),
        "L_positive_lower": _record(L_floor),
        "log_L_lower_for_reciprocal": _record(log_L_lower),
        "h": _record(h),
        "d": _record(d),
        "source_signed_values_selected": False,
        "delta_positive_and_less_than_one_checked": True,
        "q1_denominator_derivative_correction_included": True,
        "q1_denominator_correction": (
            None if denominator_correction is None
            else _record(denominator_correction)
        ),
    }
    return qy, receipt, n1


def _fixture_bound_check() -> dict[str, Any]:
    """Exercise the actual log-cap helpers on a fresh rational fixture."""

    import sympy as sp

    F = Fraction
    c = MPIntervalContext()
    c.prec = 260
    delta = F(1, 11)
    z = F(1, 5)
    N = F(13)
    pstar_e0 = F(7)
    slow_y1 = [F(2, 3), F(1, 2), F(1, 3), F(1, 4)]
    slow_y0 = [F(1, 3), F(1, 4), F(1, 5), F(1, 6)]
    fast_phi_y0 = [F(1, 7), F(1, 8), F(1, 9), F(1, 10)]
    fast_phi_y1 = [F(1, 11), F(1, 12), F(1, 13), F(1, 14)]
    delta_u = [F(7, 100), F(1, 20), F(1, 30), F(1, 40)]
    delta_m = [F(1, 50), F(1, 60), F(1, 70), F(1, 80)]
    delta_my_signed = [delta_u[k] - delta_m[k] for k in ORDERS]

    def _scalar(value):
        return c.mpf(value.numerator) / value.denominator

    def _log_abs(value):
        if value == 0:
            return None
        return c.ln(_scalar(abs(value)))

    def _mixed_caps(y0, y1):
        return {
            (i, k): _log_abs((y0 if i == 0 else y1)[k])
            for i in (0, 1) for k in ORDERS
        }

    slow_caps = _mixed_caps(slow_y0, slow_y1)
    fast_caps = _mixed_caps(fast_phi_y0, fast_phi_y1)
    delta_u_caps = {k: _log_abs(delta_u[k]) for k in ORDERS}
    delta_m_caps = {k: _log_abs(delta_m[k]) for k in ORDERS}
    log_pe = _log_abs(pstar_e0)
    log_n = _log_abs(N)
    delta_uy_caps, slow_terms, fast_terms = _finite_delta_uy(
        c, slow_caps, fast_caps, log_pe, log_n
    )
    delta_my_caps = _finite_delta_my(c, delta_u_caps, delta_m_caps)
    delta_qy_caps, pressure_receipt, c1_cap = _finite_delta_qy(
        c, delta_uy_caps, delta_my_caps, _scalar(delta), _scalar(z)
    )

    x = sp.symbols("x")
    x0 = sp.Rational(z.numerator, z.denominator)

    def _jet_expression(rows):
        return sum(
            sp.Rational(rows[k].numerator, rows[k].denominator)
            * (x - x0) ** k / sp.factorial(k)
            for k in ORDERS
        )

    # Independent signed source functions.  Their ordinary derivatives at x0
    # are exactly the rational rows above, so SymPy differentiates the same
    # N/L operator used by the pressure cap implementation.
    du_signed = [
        pstar_e0 * (slow_y1[k] / N + fast_phi_y0[k]) for k in ORDERS
    ]
    my_signed = delta_my_signed
    du_expr = _jet_expression(du_signed)
    my_expr = _jet_expression(my_signed)
    delta_expr = sp.Rational(delta.numerator, delta.denominator)
    h_expr = 1 - delta_expr
    d_expr = 1 - x**2
    L_expr = 1 - delta_expr * x**2
    numerator_expr = sp.expand(
        2 * x * du_expr - h_expr * x * my_expr
        - d_expr * sp.diff(my_expr, x)
    )
    q_expr = sp.cancel(numerator_expr / L_expr)
    direct_q0 = sp.cancel(q_expr.subs(x, x0))
    direct_q1 = sp.cancel(sp.diff(q_expr, x).subs(x, x0))
    numerator_0 = sp.cancel(numerator_expr.subs(x, x0))
    numerator_1 = sp.cancel(sp.diff(numerator_expr, x).subs(x, x0))
    direct_q1_without_denominator_term = sp.cancel(
        numerator_1 / L_expr.subs(x, x0)
    )
    direct_q1_with_split_formula = sp.cancel(
        numerator_1 / L_expr.subs(x, x0)
        + numerator_0 * 2 * delta_expr * x0
        / L_expr.subs(x, x0) ** 2
    )

    def _fraction(value):
        value = sp.Rational(value)
        return F(int(value.p), int(value.q))

    def _contains_log(cap, value):
        expected = _log_abs(value)
        return expected is None or (cap is not None and
                                    _hi(cap) >= _hi(expected))

    direct_q0_f = _fraction(direct_q0)
    direct_q1_f = _fraction(direct_q1)
    numerator_0_f = _fraction(numerator_0)
    numerator_1_f = _fraction(numerator_1)
    direct_q1_without_denominator_term_f = _fraction(
        direct_q1_without_denominator_term
    )
    direct_q1_with_split_formula_f = _fraction(direct_q1_with_split_formula)
    direct_du_logs = [_log_abs(value) for value in du_signed]
    direct_my_logs = [_log_abs(value) for value in my_signed]
    direct_q_logs = [_log_abs(direct_q0_f), _log_abs(direct_q1_f)]
    expected_slow_logs = [
        _log_abs(pstar_e0 * slow_y1[k] / N) for k in ORDERS
    ]
    expected_fast_logs = [_log_abs(pstar_e0 * fast_phi_y0[k])
                          for k in ORDERS]
    slow_term_checks = [
        slow_terms[k] is not None and
        _lo(slow_terms[k]) <= _lo(expected_slow_logs[k]) and
        _hi(slow_terms[k]) >= _hi(expected_slow_logs[k])
        for k in ORDERS
    ]
    fast_term_checks = [
        fast_terms[k] is not None and
        _lo(fast_terms[k]) <= _lo(expected_fast_logs[k]) and
        _hi(fast_terms[k]) >= _hi(expected_fast_logs[k])
        for k in ORDERS
    ]
    q1_base = c1_cap - _read_interval(
        c, pressure_receipt["log_L_lower_for_reciprocal"]
    )
    q1_correction = _read_interval(
        c, pressure_receipt["q1_denominator_correction"]
    )
    q1_cap_from_terms = _logsum(c, (q1_base, q1_correction))
    checks = {
        "fresh_context_prec_260": c.prec == 260,
        "slow_terms_match_rational_ln": all(slow_term_checks),
        "fast_terms_match_rational_ln_without_N": all(fast_term_checks),
        "fast_term_would_change_if_divided_by_N": (
            expected_fast_logs[0].b != _log_abs(pstar_e0 * fast_phi_y0[0] / N).b
        ),
        "deltaUy_logcaps_contain_signed_source": all(
            _contains_log(delta_uy_caps[k], du_signed[k]) for k in ORDERS
        ),
        "deltaMy_logcaps_contain_signed_difference": all(
            _contains_log(delta_my_caps[k], my_signed[k]) for k in ORDERS
        ),
        "deltaQy_logcaps_contain_sympy_Z0_Z1": (
            _contains_log(delta_qy_caps[0], direct_q0_f) and
            _contains_log(delta_qy_caps[1], direct_q1_f)
        ),
        "sympy_direct_diff_Z0_matches_rational": direct_q0_f == _fraction(direct_q0),
        "sympy_direct_diff_Z1_matches_rational": direct_q1_f == _fraction(direct_q1),
        "q1_denominator_correction_nonzero": (
            direct_q1 != direct_q1_without_denominator_term
        ),
        "q1_split_formula_matches_direct_diff": (
            direct_q1 == direct_q1_with_split_formula
        ),
        "q1_cap_contains_split_terms": (
            _hi(delta_qy_caps[1]) >= _hi(q1_cap_from_terms)
        ),
        "q1_denominator_receipt": (
            pressure_receipt["q1_denominator_derivative_correction_included"] is True
        ),
        "positive_L_lower": 1 - delta * z * z > 0,
        "C1_numerator_cap_contains_direct": _contains_log(c1_cap, numerator_1_f),
        "signed_deltaMy_difference_retained": (
            my_signed[0] == delta_u[0] - delta_m[0]
        ),
    }
    module_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return {
        "passed": all(checks.values()),
        "schema": "directed_Rm_endpoint_rational_fixture_v2",
        "module_sha256": module_sha256,
        "context_precision": c.prec,
        "q1_denominator_derivative_correction_included": True,
        "delta": str(delta),
        "Z": str(z),
        "N": str(N),
        "Pstar_E0": str(pstar_e0),
        "slow_y1": [str(value) for value in slow_y1],
        "fast_phi_y0": [str(value) for value in fast_phi_y0],
        "deltaUy": [str(value) for value in du_signed],
        "signed_deltaU": [str(value) for value in delta_u],
        "signed_deltaM": [str(value) for value in delta_m],
        "signed_deltaMy": [str(value) for value in my_signed],
        "deltaQy_direct_sympy": [str(direct_q0_f), str(direct_q1_f)],
        "deltaQy_Z1_naive_n1_over_L": str(
            direct_q1_without_denominator_term_f
        ),
        "deltaQy_Z1_numerator_n1": str(numerator_1_f),
        "deltaQy_Z1_split_formula": str(direct_q1_with_split_formula_f),
        "C1_numerator_Z0": str(numerator_1_f),
        "checks": checks,
    }


def fixture_bound_check() -> dict[str, Any]:
    """Return the standalone rational fixture receipt.

    The receipt hashes this module's final bytes so a runner can distinguish
    the fixture it executed from a stale saved receipt.
    """

    return _fixture_bound_check()


def compute(c: MPIntervalContext, mixed_native_packet: Mapping[str, Any],
            mixed_loop_packet: Mapping[str, Any], finite_report_bytes: bytes,
            native_report_bytes: bytes) -> dict[str, Any]:
    """Build the signed/log-scale post-restore ``Rm`` endpoint packet.

    ``finite_report_bytes`` and ``native_report_bytes`` are authenticated as
    complete reviewed reports before any nested data is read.  All finite
    corrections are absolute logarithmic caps; the native ``Qy`` source rows
    remain signed directed intervals.
    """

    _require(isinstance(c, MPIntervalContext),
             "c must be an MPIntervalContext")
    finite_report = _admit_report(
        finite_report_bytes, FINITE_REPORT_BYTES, FINITE_REPORT_SHA256,
        "Finite report",
    )
    native_report = _admit_report(
        native_report_bytes, NATIVE_REPORT_BYTES, NATIVE_REPORT_SHA256,
        "Native report",
    )
    native_chart, loop_chart = _validate_packets(
        mixed_native_packet, mixed_loop_packet, NATIVE_REPORT_SHA256,
    )
    native = _native_sources(c, native_report)
    bounds, saved_delta_u, saved_delta_m, saved_delta_q, log_n, log_pe = \
        _finite_sources(c, finite_report)

    slow = _packet_rows(c, mixed_loop_packet, loop_chart,
                        "B_over_Pstar_over_frozen_E0")
    fast = _packet_rows(c, mixed_loop_packet, loop_chart,
                        "B_phi_over_Pstar_over_frozen_E0")
    packet_log_n = _packet_scalar(c, mixed_loop_packet, loop_chart,
                                  ("logN", "log_N"))
    packet_log_pe = _packet_scalar(
        c, mixed_loop_packet, loop_chart,
        ("log_Pstar_E0_upper", "log_Pstar_E_upper", "log_Pstar_E0"),
    )
    if packet_log_n is not None:
        _require(_same_interval(packet_log_n, log_n),
                 "Mixed-loop logN differs from finite report logN")
    if packet_log_pe is not None:
        _require(_same_interval(packet_log_pe, log_pe),
                 "Mixed-loop log_Pstar_E0 differs from finite report scale")
    for owner in (mixed_loop_packet, loop_chart):
        if not isinstance(owner, Mapping):
            continue
        for key in ("B_phi_has_no_N_factor", "fast_term_has_no_1_over_N",
                    "fast_phi_has_no_1_over_N"):
            if key in owner:
                _require(owner[key] is True, key + " must be true")

    delta_uy, slow_terms, fast_terms = _finite_delta_uy(
        c, slow, fast, log_pe, log_n,
    )
    delta_my = _finite_delta_my(c, saved_delta_u, saved_delta_m)
    delta_qy, pressure_receipt, c1_numerator = _finite_delta_qy(
        c, delta_uy, delta_my, native["delta"], native["z"],
    )

    log_rm = native["log_rm"]
    native_abs = {k: _absolute_log(c, native["qy"][k],
                                   "native Qy[" + str(k) + "]")
                  for k in (0, 1)}
    native_norm = {k: (None if native_abs[k] is None else native_abs[k] - log_rm)
                   for k in (0, 1)}
    delta_norm = {k: (None if delta_qy[k] is None else delta_qy[k] - log_rm)
                  for k in (0, 1)}
    combined_norm = {k: _logsum(c, (native_norm[k], delta_norm[k]))
                     for k in (0, 1)}
    threshold = -c.mpf(1024) * c.ln(c.mpf(2))
    threshold_passed = all(
        combined_norm[k] is None or _hi(combined_norm[k]) < _lo(threshold)
        for k in (0, 1)
    )
    honest_enclosure = None
    if threshold_passed:
        radius = c.mpf(2) ** (-1024)
        negative_radius = -radius
        honest_enclosure = _record(
            c.mpf([_lo(negative_radius), _hi(radius)]), "2^-1024 enclosure"
        )

    native_signed_log_scale = {
        k: _signed_log_scale(c, native["qy"][k], log_rm,
                             "native Qy[" + str(k) + "] / Rm")
        for k in (0, 1)
    }
    packet_metadata = {
        "mixed_native_schema": mixed_native_packet.get("schema"),
        "mixed_loop_schema": mixed_loop_packet.get("schema"),
        "mixed_native_source_native_report_sha256": mixed_native_packet.get(
            "source_native_report_sha256"
        ),
        "mixed_loop_source_native_parent_report_sha256": mixed_loop_packet.get(
            "source_native_parent_report_sha256"
        ),
        "chart": CHART,
        "model_metadata_policy": "only report or packet metadata is exposed; no model capability is inferred",
    }
    result = {
        "passed": True,
        "schema": "original_directed_Rm_endpoint_kernel_v1",
        "stage": "current_original_directed_Rm_endpoint_kernel",
        "accepted_ancestors_replayed": False,
        "target_status": {
            "computation_passed": True,
            "structure_passed": True,
            "normalized_2_power_minus_1024_target_passed": threshold_passed,
        },
        "chart": CHART,
        "source_reports": {
            "native": {
                "path": NATIVE_REPORT_PATH,
                "bytes": NATIVE_REPORT_BYTES,
                "sha256": NATIVE_REPORT_SHA256,
            },
            "finite": {
                "path": FINITE_REPORT_PATH,
                "bytes": FINITE_REPORT_BYTES,
                "sha256": FINITE_REPORT_SHA256,
            },
        },
        "native_endpoint": {
            "source_path": (
                "actual_restoration_post_jet_repair.corrected_background."
                "charts.post_restore.endpoint_fields.Qy_end[:2]"
            ),
            "source_rows": native["source_rows"],
            "signed_Qy": [
                dict(ordinary_Z_order=k, value=_record(native["qy"][k]))
                for k in (0, 1)
            ],
            "endpoint_recomputed": False,
        },
        "finite_deltaUy": {
            "formula": (
                "deltaUy=Pstar*E0*(Bslow[1,k]/N+Bphi[0,k]); "
                "Bphi fast term has no 1/N"
            ),
            "slow_source": "B_over_Pstar_over_frozen_E0[physical_y_order=1,k]",
            "fast_source": "B_phi_over_Pstar_over_frozen_E0[physical_y_order=0,k]",
            "slow_terms": _pack_log_map(slow_terms, ORDERS),
            "fast_terms": _pack_log_map(fast_terms, ORDERS),
            "log_absolute_upper_bounds": _pack_log_map(delta_uy, ORDERS),
            "source_signed_value_selected": False,
        },
        "finite_deltaMy": {
            "formula": "deltaMy=deltaU-deltaM",
            "cap_rule": "absolute upper bound uses log_sum(deltaU,deltaM); signed cancellation is not claimed",
            "delta_U": _pack_log_map(saved_delta_u, ORDERS),
            "delta_M": _pack_log_map(saved_delta_m, ORDERS),
            "log_absolute_upper_bounds": _pack_log_map(delta_my, ORDERS),
            "source_signed_value_selected": False,
        },
        "finite_deltaQy": {
            "ordinary_Z_orders": list(ENDPOINT_ORDERS),
            "log_absolute_upper_bounds": _pack_log_map(delta_qy, ENDPOINT_ORDERS),
            "operator_receipt": pressure_receipt,
            "source_signed_value_selected": False,
        },
        "finite_C1_numerator": {
            "ordinary_Z_orders": [0],
            "log_absolute_upper_bounds": _pack_log_map(
                {0: c1_numerator}, (0,)
            ),
            "formula": pressure_receipt["C1_formula"],
            "kind": "C1 numerator at ordinary Z order 0; no Z1 derivative row claimed",
            "source_signed_value_selected": False,
        },
        "saved_finite_error_bounds": {
            "delta_U": _pack_log_map(saved_delta_u, ORDERS),
            "delta_M": _pack_log_map(saved_delta_m, ORDERS),
            "delta_Q": _pack_log_map(saved_delta_q, Q_ORDERS),
            "logN": _record(log_n),
            "log_Pstar_E_upper": _record(log_pe),
            "source_path": (
                "actual_directed_finite_N_controls.charts.post_restore."
                "physical_error_log_bounds"
            ),
        },
        "normalization": {
            "log_Rm": _record(log_rm),
            "log_Rm_source_record": _record(log_rm),
            "log_Rm_source_path": (
                "actual_restoration_background_I1_numeric.physical_scale.log_Rm"
            ),
            "log_Rm_background_chart_record": _record(
                native["background_chart_log_rm"]
            ),
            "log_Rm_corrected_scale_record": _record(
                native["corrected_scale_log_rm"]
            ),
            "log_Rm_chart_record": _record(native["chart_log_rm"]),
            "log_Rm_records_equal": True,
            "log_Rm_reciprocal_uses_full_interval": True,
            "Wpost": "Rm^2",
            "log_Wpost": _record(2 * log_rm),
            "RQy_over_Wpost": "Qy/Rm",
            "native_Qy_over_Rm_signed_log_scale": [
                dict(ordinary_Z_order=k, **native_signed_log_scale[k])
                for k in (0, 1)
            ],
            "finite_delta_Qy_over_Rm_log_absolute_bounds": _pack_log_map(
                delta_norm, (0, 1)
            ),
            "combined_Qy_over_Rm_log_absolute_bounds": _pack_log_map(
                combined_norm, (0, 1)
            ),
            "threshold": _record(threshold),
            "threshold_description": "combined normalized log cap < -1024*log(2)",
            "threshold_passed": threshold_passed,
            "honest_symmetric_enclosure": honest_enclosure,
            "symmetric_enclosure_not_actual_zero_or_value": True,
            "radius_materialized": False,
            "astronomical_logRm_materialized": False,
        },
        "pressure_boundary": {
            "quantity": "-Rm*Qy",
            "coefficient": -1,
            "external_source_handles": list(EXTERNAL_RQY_HANDLES),
            "signed_coefficients_by_Z": {str(k): -1 for k in (0, 1)},
            "source_identity": (
                "postrestore pressure join retains -Rm*Qy; the RQy source "
                "handle carries coefficient -1; (Rm*Qy)/Wpost=Qy/Rm"
            ),
        },
        "native_mixed_receipts": {
            "physical_y_identities": copy.deepcopy(
                mixed_native_packet.get("physical_y_identities", {
                    "phi_y": "-2/5*phi", "U_y": "0", "M_y": "U-M",
                    "A_y": "U^2-A", "B_y": "phi^2-2*B", "C_y": "phi^2-C",
                })
            ),
            "p2_y_formula": "p2_y=(9/10)*p2+R*Nnative_y/(Utheta*L)",
            "p2_y_full_product_retained": True,
            "root_mixed_formula_receipt": {
                "u": "p2*q/dstar", "h": "(1+u^2)^(-1/2)",
                "r": "u*h", "alpha": "2*q*h",
                "K": "1/[2*pi*(1+t0^2+2*q^2)]",
            },
            "derivative_orders": {
                "physical_y_orders": [0, 1],
                "ordinary_Z_orders": list(ORDERS),
                "mixed_rows": 8,
                "LogPhase_total_order_max": 4,
            },
            "positive_floors": copy.deepcopy(
                mixed_native_packet.get("root_positive_floors",
                                       mixed_native_packet.get("conditioning"))
            ),
        },
        "cancellation_review": {
            "deltaUy": "slow and fast are summed only as absolute caps; fast has no 1/N",
            "deltaMy": "source definition is signed deltaU-deltaM, but emitted cap uses addition",
            "native_Qy": "signed Qy rows are copied and kept signed through log-scale normalization",
            "Rm_endpoint_vs_volume": "local and endpoint pressure terms may cancel by signed source graph; the surviving external -Rm*Qy term is represented by the RQy source handle with coefficient -1",
            "finite_correction_vs_native_endpoint": "finite deltaQy is an absolute cap and is not assigned a sign or cancelled against native Qy",
        },
        "packet_metadata": packet_metadata,
        "diagnostic_2_power_minus_1024": {
            "threshold_passed": threshold_passed,
            "may_fail": True,
            "report_only": True,
        },
        "contracts": {
            "accepted_producers_called": False,
            "accepted_builders_called": False,
            "owners_called": False,
            "accepted_checkers_called": False,
            "recurrences_called": False,
            "quadratures_called": False,
            "norms_called": False,
            "native_Qy_recomputed": False,
            "source_values_selected": False,
            "caps_used_as_signed_values": False,
            "symmetric_interval_used_as_zero": False,
            "N_materialized": False,
            "Rm_materialized": False,
            "Wpost_applied_once": True,
        },
    }
    json.dumps(result, allow_nan=False, separators=(",", ":"))
    return result


__all__ = [
    "compute",
    "fixture_bound_check",
    "NATIVE_REPORT_BYTES",
    "NATIVE_REPORT_SHA256",
    "FINITE_REPORT_BYTES",
    "FINITE_REPORT_SHA256",
]
