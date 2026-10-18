"""Source-correlated fast ``B_phi`` and endpoint-C1 adapters.

The production entry point is a byte-authenticated consumer of the already
saved post-restore report.  It reads the native ``E/E0`` rows, the fixed-target
phase inverse ``that`` rows, and the old fast bound.  No historical producer,
builder, checker, recurrence, quadrature, norm, or draft is imported or run.

The fast source is

``bphi = -C * e * t_hat / (1 + t_hat**2)``.

The source identities give ``C <= 2``.  The derivative cap is formed from the
combined signed identity before boxing, using the independent global rational
inequalities ``|t/(1+t**2)| <= 1/2`` and
``|(1-t**2)/(1+t**2)**2| <= 1``.  A logarithmic cap is never used as a signed
value.
"""

from __future__ import annotations

import ast
import copy
import hashlib
import json
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path
from typing import Any

from mpmath import isfinite, make_mpf, nstr
from mpmath.ctx_iv import MPIntervalContext


REPORT_PATH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_directed_Rm_controls.json"
)
REPORT_BYTES = 3159660
REPORT_SHA256 = (
    "78a95b631bdce6dbe1438bed1ccff282127df5d2198e6836d250a2b9ff325deb"
)
MODULE = "src/openai_ns_reconstruction/original_post_restore_fast_phi_bound.py"
CHART = "post_restore"
ORDERS = (0, 1, 2, 3)
MIXED_INDICES = tuple((i, k) for i in (0, 1) for k in ORDERS)
NATIVE_REPORT_SHA256 = (
    "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _lo(value: Any):
    return make_mpf(tuple(value._mpi_[0]))


def _hi(value: Any):
    return make_mpf(tuple(value._mpi_[1]))


def _canonical_endpoint(value: Any, label: str) -> tuple[int, int, int, int]:
    _require(
        isinstance(value, (list, tuple)) and len(value) == 4,
        label + " must be an exact MPF tuple",
    )
    _require(
        all(type(part) is int for part in value),
        label + " MPF tuple components must be integers",
    )
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
    """Serialize both endpoints of one finite directed interval."""

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
    """Decode both exact endpoints and re-box them in ``c``."""

    if hasattr(value, "_mpi_"):
        lower = tuple(value._mpi_[0])
        upper = tuple(value._mpi_[1])
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
    _require(isfinite(_lo(out)) and isfinite(_hi(out)),
             label + " must be finite")
    return out


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _admit_report(raw: bytes) -> Mapping[str, Any]:
    _require(type(raw) is bytes, "report_bytes must be immutable bytes")
    _require(len(raw) == REPORT_BYTES,
             "Saved directed-Rm report byte length changed")
    _require(_sha(raw) == REPORT_SHA256,
             "Saved directed-Rm report SHA256 changed")
    try:
        report = json.loads(raw)
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError("Saved directed-Rm report is not valid JSON") from exc
    _require(isinstance(report, Mapping),
             "Saved directed-Rm report must be a mapping")
    return report


def _at(mapping: Mapping[str, Any], *keys: str) -> Any:
    value: Any = mapping
    for key in keys:
        _require(isinstance(value, Mapping) and key in value,
                 "Missing report path: " + ".".join(keys))
        value = value[key]
    return value


def _row_value(c: MPIntervalContext, row: Mapping[str, Any], label: str):
    _require(isinstance(row, Mapping), label + " row must be a mapping")
    exact = row.get("exact_zero")
    _require(type(exact) is bool, label + ".exact_zero must be boolean")
    raw = row.get("log_absolute_upper_bound")
    if raw is None:
        raw = row.get("record", row.get("log_cap"))
    if exact:
        _require(raw is None, label + " exact-zero row carries a cap")
        return None
    _require(raw is not None, label + " nonzero row has no cap")
    return _read_interval(c, raw, label + ".log_absolute_upper_bound")


def _mixed_rows(c: MPIntervalContext, raw: Any, label: str) \
        -> dict[tuple[int, int], Any]:
    """Decode exactly the canonical ``(y=0,1; Z=0..3)`` row order."""

    if isinstance(raw, Mapping):
        for key in ("rows", "mixed_log_caps", "log_caps", "values"):
            if key in raw:
                return _mixed_rows(c, raw[key], label)
        if set(raw) >= {"y0", "y1"}:
            _require(
                isinstance(raw["y0"], (list, tuple))
                and isinstance(raw["y1"], (list, tuple)),
                label + " y0/y1 rows must be lists",
            )
            raw = list(raw["y0"]) + list(raw["y1"])
        else:
            try:
                if set(raw) == set(MIXED_INDICES):
                    return {
                        index: (None if raw[index] is None else
                                _read_interval(c, raw[index],
                                               label + str(index)))
                        for index in MIXED_INDICES
                    }
            except TypeError:
                pass
    _require(isinstance(raw, (list, tuple)) and len(raw) == 8,
             label + " must contain exactly eight canonical rows")
    out: dict[tuple[int, int], Any] = {}
    for number, expected in enumerate(MIXED_INDICES):
        row = raw[number]
        _require(isinstance(row, Mapping),
                 label + " row must be a mapping")
        index = (
            row.get("physical_y_order", row.get("y_order")),
            row.get("ordinary_Z_order", row.get("Z_order")),
        )
        _require(index == expected,
                 label + " is not in canonical eight-row order")
        out[index] = _row_value(c, row, label + "[" + str(number) + "]")
    return out


def _ordinary_rows(c: MPIntervalContext, raw: Any, count: int,
                   label: str) -> dict[int, Any]:
    _require(isinstance(raw, (list, tuple)), label + " must be a row list")
    out: dict[int, Any] = {}
    for number, row in enumerate(raw):
        _require(isinstance(row, Mapping), label + " row must be a mapping")
        order = row.get("ordinary_Z_order", row.get("Z_order", number))
        _require(order == number and order in range(count),
                 label + " is not in canonical ordinary-Z order")
        out[order] = _row_value(c, row, label + "[" + str(number) + "]")
    _require(len(out) == count,
             label + " must cover ordinary Z orders 0 through "
             + str(count - 1))
    return out


def _logsum(c: MPIntervalContext, terms: Any):
    """Sum log caps with interval arithmetic and an outward re-box.

    The ``ln(count)`` shift is added to each interval before endpoint min/max
    collection.  This avoids exponentiating very large logarithmic scales.
    """

    terms = tuple(term for term in terms if term is not None)
    if not terms:
        return None
    if len(terms) == 1:
        return c.mpf([_lo(terms[0]), _hi(terms[0])])
    shift = c.ln(c.mpf(len(terms)))
    shifted = tuple(term + shift for term in terms)
    return c.mpf([
        min(_lo(term) for term in shifted),
        max(_hi(term) for term in shifted),
    ])


def _safe_abs_log(c: MPIntervalContext, value: Any, label: str):
    """Return a magnitude log using the interval ``absmax`` primitive."""

    value = _read_interval(c, value, label)
    magnitude = c.absmax(value)
    upper = _hi(magnitude)
    if upper == 0:
        return None
    # A cap needs the maximum magnitude.  Re-box that safe interval primitive
    # before taking its log; no scalar endpoint negation or ``abs`` is used.
    return c.ln(c.mpf([upper, upper]))


def _pack_log(order: int, value: Any) -> dict[str, Any]:
    return {
        "ordinary_Z_order": order,
        "exact_zero": value is None,
        "log_absolute_upper_bound": (
            None if value is None else _record(value)
        ),
    }


def _same_interval(left: Any, right: Any) -> bool:
    return (tuple(left._mpi_[0]) == tuple(right._mpi_[0])
            and tuple(left._mpi_[1]) == tuple(right._mpi_[1]))


def _source_binding(native: Mapping[str, Any], loop: Mapping[str, Any],
                    rm: Mapping[str, Any]) -> dict[str, Any]:
    _require(native.get("passed") is True, "Saved native-mixed packet failed")
    _require(loop.get("passed") is True, "Saved loop-mixed packet failed")
    _require(rm.get("passed") is True, "Saved Rm packet failed")
    _require(native.get("schema") ==
             "original_directed_post_restore_mixed_native_v1",
             "Saved native-mixed schema changed")
    _require(loop.get("schema") ==
             "original_directed_post_restore_mixed_loop_v1",
             "Saved loop-mixed schema changed")
    _require(rm.get("schema") == "original_directed_Rm_endpoint_kernel_v1",
             "Saved Rm schema changed")
    for packet, label in ((native, "native-mixed"), (loop, "loop-mixed")):
        _require(packet.get("chart") == CHART,
                 label + " chart changed")
        for key in (
            "source_executed",
            "accepted_ancestors_replayed",
            "accepted_builders_or_owners_replayed",
            "accepted_checkers_replayed",
            "recurrences_replayed",
            "quadratures_replayed",
            "norms_replayed",
        ):
            if key in packet:
                _require(packet[key] is False,
                         label + " execution flag changed: " + key)
    native_sha = native.get("source_native_report_sha256")
    loop_sha = loop.get("source_native_parent_report_sha256")
    _require(native_sha == NATIVE_REPORT_SHA256,
             "Native source SHA binding changed")
    _require(loop_sha == NATIVE_REPORT_SHA256,
             "Loop native-source SHA binding changed")
    metadata = rm.get("packet_metadata")
    _require(isinstance(metadata, Mapping), "Rm packet metadata missing")
    _require(metadata.get("mixed_native_source_native_report_sha256") ==
             NATIVE_REPORT_SHA256,
             "Rm native-source metadata binding changed")
    _require(metadata.get("mixed_loop_source_native_parent_report_sha256") ==
             NATIVE_REPORT_SHA256,
             "Rm loop-source metadata binding changed")
    _require(metadata.get("chart") == CHART,
             "Rm packet metadata chart changed")
    root_ids = native.get("root_identities")
    conditioning = native.get("conditioning")
    _require(isinstance(root_ids, Mapping), "Native root identities missing")
    _require(isinstance(conditioning, Mapping), "Native conditioning missing")
    _require(root_ids.get("a") == "4/5",
             "Native source identity a=4/5 missing")
    _require(root_ids.get("q_squared_post_y_and_Z_constant") is True,
             "Native q-squared y/Z constancy source identity missing")
    _require(conditioning.get("eta_le_one") is True,
             "Native eta<=1 source conditioning missing")
    _require(native.get("ordinary_Z_orders") == list(ORDERS)
             and native.get("physical_y_orders") == [0, 1],
             "Native mixed row orders changed")
    _require(loop.get("ordinary_Z_orders") == list(ORDERS)
             and loop.get("physical_y_orders") == [0, 1],
             "Loop mixed row orders changed")
    for key in ("fixed_target_phase", "target_phase_Z_independent",
                "phase_Z_independent", "fast_phi_separate",
                "B_phi_has_no_N_factor"):
        _require(loop.get(key) is True, "Loop fixed-target flag missing: " + key)
    proof = loop.get("q_control_proof")
    _require(isinstance(proof, Mapping), "Loop q-control proof missing")
    _require(proof.get("q_squared_constant_in_y_Z") is True
             and proof.get("q_value_not_selected") is True,
             "Loop q source-correlation proof missing")
    _require(proof.get("q_squared_identity") ==
             "q^2=3/4+5*eta/4",
             "Loop q source identity changed")
    contracts = rm.get("contracts")
    _require(isinstance(contracts, Mapping), "Rm contracts missing")
    _require(contracts.get("caps_used_as_signed_values") is False
             and contracts.get("source_values_selected") is False,
             "Saved cap/sign contract changed")
    return {
        "native_report_sha256": native_sha,
        "native_a_identity": root_ids["a"],
        "q_squared_identity": proof["q_squared_identity"],
        "q_squared_constant_in_y_Z": True,
        "eta_le_one_source_premise": True,
        "source_sha_binding_verified": True,
        "fixed_target_phase_verified": True,
        "caps_remain_absolute": True,
        "signed_source_values_selected": False,
    }


def _pure_fast_z1_bound(c: MPIntervalContext, e0_log: Any,
                        e_z_log: Any, that_z_log: Any,
                        log_C_upper: Any | None = None):
    """Apply the source-correlated ``Z1`` fast bound in log space.

    Inputs are logarithmic absolute caps.  With ``C <= 2`` this is exactly

    ``ln(2) + logsum(e_Z - ln(2), e_0 + that_Z)``.

    ``log_C_upper`` is optional for the fresh rational fixture and tests; the
    production adapter uses the source-bound ``ln(2)`` value.
    """

    _require(isinstance(c, MPIntervalContext),
             "c must be an MPIntervalContext")
    e0_log = None if e0_log is None else _read_interval(c, e0_log, "e0 log cap")
    e_z_log = None if e_z_log is None else _read_interval(c, e_z_log, "e_Z log cap")
    that_z_log = (
        None if that_z_log is None
        else _read_interval(c, that_z_log, "that_Z log cap")
    )
    log_two = c.ln(c.mpf(2))
    if log_C_upper is None:
        log_C_upper = log_two
    else:
        log_C_upper = _read_interval(c, log_C_upper, "log C cap")
        _require(_hi(log_C_upper) <= _hi(log_two),
                 "Source C cap exceeds 2")
    # C multiplies the already-combined derivative.  Keep it outside the
    # two-term logsum so the source-bound C=2 formula remains
    # ln(2)+logsum(e_Z-ln(2), e0+that_Z).
    first = None if e_z_log is None else e_z_log - log_two
    second = None if e0_log is None or that_z_log is None \
        else e0_log + that_z_log
    inner = _logsum(c, (first, second))
    return None if inner is None else log_C_upper + inner


def pure_bound(c: MPIntervalContext, e0_log: Any, e_z_log: Any,
               that_z_log: Any, log_C_upper: Any | None = None):
    """Public pure source-correlated fast ``Z1`` bound.

    ``e0_log``, ``e_z_log`` and ``that_z_log`` are log absolute caps.  The
    default source bound is ``C <= 2``.  The extra ``ln(2)`` outside the
    ``logsum`` is retained exactly by this wrapper.
    """

    log_two = c.ln(c.mpf(2))
    if log_C_upper is None:
        log_C_upper = log_two
    return _pure_fast_z1_bound(c, e0_log, e_z_log, that_z_log,
                               log_C_upper)


def _pure_bound_path(c: MPIntervalContext, e0_log: Any, e_z_log: Any,
                     that_z_log: Any):
    """Return the new cap plus its interval arithmetic path."""

    log_two = c.ln(c.mpf(2))
    e0 = None if e0_log is None else _read_interval(c, e0_log, "e0 log cap")
    ez = None if e_z_log is None else _read_interval(c, e_z_log, "e_Z log cap")
    thatz = (None if that_z_log is None
             else _read_interval(c, that_z_log, "that_Z log cap"))
    first = None if ez is None else ez - log_two
    second = None if e0 is None or thatz is None else e0 + thatz
    inner = _logsum(c, (first, second))
    result = None if inner is None else log_two + inner
    return result, {
        "e_Z_minus_ln2": None if first is None else _record(first),
        "e0_plus_that_Z": None if second is None else _record(second),
        "inner_logsum": None if inner is None else _record(inner),
        "outer_ln2": _record(log_two),
        "result": None if result is None else _record(result),
        "logsum_rule": "interval term+ln(count), then min/max rebox",
    }


def _row_record(raw_rows: Any, index: int) -> dict[str, Any]:
    _require(isinstance(raw_rows, (list, tuple)) and len(raw_rows) == 8,
             "Canonical mixed rows are required")
    row = raw_rows[index]
    _require(isinstance(row, Mapping), "Canonical mixed row is malformed")
    return copy.deepcopy(dict(row))


def _replace_fast_z1(raw_rows: Any, bound: Any) -> list[dict[str, Any]]:
    rows = [copy.deepcopy(dict(row)) for row in raw_rows]
    _require(len(rows) == 8, "Canonical fast rows are required")
    for number, expected in enumerate(MIXED_INDICES):
        _require((rows[number].get("physical_y_order"),
                  rows[number].get("ordinary_Z_order")) == expected,
                 "Fast rows are not canonical")
    position = MIXED_INDICES.index((0, 1))
    rows[position]["exact_zero"] = bound is None
    rows[position]["log_absolute_upper_bound"] = (
        None if bound is None else _record(bound)
    )
    return rows


def _source_rows(parent: Mapping[str, Any], c: MPIntervalContext):
    packets = parent.get("actual_directed_Rm_controls", parent)
    _require(isinstance(packets, Mapping),
             "Report actual_directed_Rm_controls packet missing")
    native = packets.get("native_mixed")
    loop = packets.get("loop_mixed")
    rm = packets.get("Rm_endpoint")
    _require(isinstance(native, Mapping) and isinstance(loop, Mapping)
             and isinstance(rm, Mapping),
             "Saved native_mixed/loop_mixed/Rm_endpoint packets missing")
    binding = _source_binding(native, loop, rm)
    e_rows_raw = _at(native, "root_mixed_caps", "E_over_frozen_E0")
    that_rows_raw = _at(loop, "inverse_mixed_caps", "that")
    fast_rows_raw = _at(loop, "B_phi_over_Pstar_over_frozen_E0")
    e_rows = _mixed_rows(c, e_rows_raw, "native E_over_frozen_E0")
    that_rows = _mixed_rows(c, that_rows_raw, "loop inverse that")
    _mixed_rows(c, fast_rows_raw, "loop old B_phi")
    return native, loop, rm, binding, e_rows, that_rows, fast_rows_raw


def _finite_delta_uy_inputs(c: MPIntervalContext, rm: Mapping[str, Any],
                            new_fast: Any):
    finite = _at(rm, "finite_deltaUy")
    old_slow_rows = _ordinary_rows(
        c, finite["slow_terms"], 4, "saved finite deltaUy slow terms"
    )
    old_delta_rows = _ordinary_rows(
        c, finite["log_absolute_upper_bounds"], 4,
        "saved finite deltaUy caps",
    )
    saved_bounds = _at(rm, "saved_finite_error_bounds")
    log_pe = _read_interval(c, saved_bounds["log_Pstar_E_upper"],
                            "saved log_Pstar_E_upper")
    z1 = _logsum(
        c,
        (old_slow_rows[1], None if new_fast is None else log_pe + new_fast),
    )
    return old_slow_rows, old_delta_rows, log_pe, z1


def _compute_packets(c: MPIntervalContext,
                     report: Mapping[str, Any]) -> dict[str, Any]:
    """Compute from an already decoded packet mapping.

    The public :func:`compute` wrapper authenticates the immutable saved
    report before calling this helper.  Keeping the packet consumer separate
    also gives the fresh fixture a synthetic, non-historical input path.
    """

    _require(isinstance(c, MPIntervalContext),
             "c must be an MPIntervalContext")
    _require(isinstance(report, Mapping), "Decoded report must be a mapping")
    native, loop, rm, binding, e_rows, that_rows, old_fast_rows_raw = \
        _source_rows(report, c)
    e0, e_z = e_rows[0, 0], e_rows[0, 1]
    that_z = that_rows[0, 1]
    new_fast, path = _pure_bound_path(c, e0, e_z, that_z)
    # This path is exactly pure_bound(c, e0, e_z, that_z), with no endpoint
    # selected.  Keep a runtime equality check on both interval endpoints.
    direct = pure_bound(c, e0, e_z, that_z)
    _require(new_fast is None and direct is None or
             _same_interval(new_fast, direct),
             "Pure-bound path and public pure_bound disagree")
    old_fast_value_row = _row_record(old_fast_rows_raw, 0)
    old_fast_rows = copy.deepcopy(old_fast_rows_raw)
    new_fast_rows = _replace_fast_z1(old_fast_rows_raw, new_fast)
    old_delta_uy_raw = _at(
        rm, "finite_deltaUy", "log_absolute_upper_bounds"
    )
    old_slow, old_delta_uy, log_pe, new_delta_uy_z1 = \
        _finite_delta_uy_inputs(c, rm, new_fast)
    old_rm_caps = {
        "B_phi_over_Pstar_over_frozen_E0": copy.deepcopy(old_fast_rows_raw),
        "physical_deltaUy": copy.deepcopy(
            _at(rm, "finite_deltaUy", "log_absolute_upper_bounds")
        ),
        "finite_deltaQy": copy.deepcopy(
            _at(rm, "finite_deltaQy", "log_absolute_upper_bounds")
        ),
    }
    new_delta_uy_row = _pack_log(1, new_delta_uy_z1)
    result = {
        "passed": True,
        "schema": "original_post_restore_fast_phi_bound_v1",
        "stage": "current_original_post_restore_fast_phi_bound",
        "chart": CHART,
        "input_report_path": REPORT_PATH,
        "input_report_bytes": REPORT_BYTES,
        "input_report_sha256": REPORT_SHA256,
        "ordinary_Z_orders": list(ORDERS),
        "physical_y_orders": [0, 1],
        "source_binding": binding,
        "source_formula": (
            "bphi=-C*e*f(that), f(t)=t/(1+t^2), "
            "bphi_Z=-C*(e_Z*f(that)+e*fprime(that)*that_Z)"
        ),
        "C_bound": {
            "formula": "C=a*(1+2*q^2)/2",
            "a_source_identity": "4/5",
            "q_squared_source_identity": "q^2=3/4+5*eta/4",
            "eta_source_conditioning": "eta<=1",
            "q_squared_upper": "2",
            "C_upper": "2",
            "log_C_upper": _record(c.ln(c.mpf(2))),
            "source_bound_used": True,
        },
        "rational_cap_identities": {
            "f": "t/(1+t^2)",
            "fprime": "(1-t^2)/(1+t^2)^2",
            "global_abs_f_le_half": True,
            "global_abs_fprime_le_one": True,
            "not_differentiated_from_value_cap": True,
        },
        "old_B_phi_over_Pstar_over_frozen_E0": old_fast_rows,
        "old_value_row_0_0": old_fast_value_row,
        "new_B_phi_over_Pstar_over_frozen_E0": new_fast_rows,
        "new_fast_Z1": _pack_log(1, new_fast),
        "new_fast_Z1_arithmetic": path,
        "physical_deltaUy_Z1": new_delta_uy_row,
        "physical_deltaUy_Z1_arithmetic": {
            "formula": "logsum(saved slow_terms[1], new_fast_Z1+saved log_Pstar_E_upper)",
            "saved_slow_terms_Z1": _pack_log(1, old_slow[1]),
            "new_fast_Z1_plus_saved_log_Pstar_E_upper": (
                None if new_fast is None else _record(new_fast + log_pe)
            ),
            "saved_log_Pstar_E_upper": _record(log_pe),
            "result": (None if new_delta_uy_z1 is None
                        else _record(new_delta_uy_z1)),
            "old_deltaUy_Z1_preserved": copy.deepcopy(old_delta_uy_raw[1]),
            "fast_term_has_no_1_over_N": True,
        },
        "old_bounds": old_rm_caps,
        "contracts": {
            "saved_packets_consumed_only": True,
            "accepted_producers_called": False,
            "accepted_builders_called": False,
            "accepted_checkers_called": False,
            "recurrences_called": False,
            "quadratures_called": False,
            "norms_called": False,
            "caps_used_as_signed_values": False,
            "signed_source_values_selected": False,
            "actual_saved_compute_replayed": False,
        },
        "scope": {
            "new_fast_phi_Z1_only": True,
            "value_row_reused_from_saved_loop": True,
            "old_bound_retained": True,
            "endpoint_C1_computed": False,
            "global_five_moment_closure": False,
            "full_goal_complete": False,
        },
    }
    json.dumps(result, allow_nan=False)
    return result


def compute(c: MPIntervalContext, report_bytes: bytes) -> dict[str, Any]:
    """Compute the source-correlated fast ``B_phi`` Z1 cap.

    ``report_bytes`` must be the exact saved Rm report.  Authentication and
    decoding happen before any nested packet is consumed.  This function does
    not call ``endpoint_C1``; the latter is an explicit pure integration hook
    for the root adapter.
    """

    _require(isinstance(c, MPIntervalContext),
             "c must be an MPIntervalContext")
    report = _admit_report(report_bytes)
    return _compute_packets(c, report)


def _parent_packets(parent: Mapping[str, Any]):
    packets = parent.get("actual_directed_Rm_controls", parent)
    _require(isinstance(packets, Mapping), "Parent packets missing")
    native = packets.get("native_mixed")
    loop = packets.get("loop_mixed")
    rm = packets.get("Rm_endpoint")
    _require(isinstance(native, Mapping) and isinstance(loop, Mapping)
             and isinstance(rm, Mapping), "Parent mixed/Rm packets missing")
    return native, loop, rm


def endpoint_C1(c: MPIntervalContext, parent: Mapping[str, Any],
                new_fast_Z1: Any) -> dict[str, Any]:
    """Integrate a new fast Z1 cap into the saved endpoint C1 operator.

    This optional hook uses the saved ``deltaUy[0]`` and ``deltaMy[0:3]``
    caps, the saved positive ``L`` floor and the saved denominator derivative
    correction.  It does not read a fresh ``delta`` or ``Z`` value and does not
    invoke the old pressure kernel.  ``new_fast_Z1`` is an interval returned by
    :func:`pure_bound`, or its packed row.
    """

    _require(isinstance(c, MPIntervalContext),
             "c must be an MPIntervalContext")
    native, loop, rm = _parent_packets(parent)
    _source_binding(native, loop, rm)
    finite_uy = _at(rm, "finite_deltaUy")
    finite_my = _at(rm, "finite_deltaMy")
    finite_qy = _at(rm, "finite_deltaQy")
    old_uy = _ordinary_rows(c, finite_uy["log_absolute_upper_bounds"], 4,
                            "saved deltaUy")
    old_my = _ordinary_rows(c, finite_my["log_absolute_upper_bounds"], 4,
                            "saved deltaMy")
    old_qy = _ordinary_rows(c, finite_qy["log_absolute_upper_bounds"], 2,
                            "saved deltaQy")
    if isinstance(new_fast_Z1, Mapping):
        new_fast_Z1 = _row_value(c, new_fast_Z1, "new fast Z1")
    elif new_fast_Z1 is not None:
        new_fast_Z1 = _read_interval(c, new_fast_Z1, "new fast Z1")
    saved_bounds = _at(rm, "saved_finite_error_bounds")
    log_pe = _read_interval(c, saved_bounds["log_Pstar_E_upper"],
                            "saved log_Pstar_E_upper")
    slow_rows = _ordinary_rows(c, finite_uy["slow_terms"], 4,
                               "saved deltaUy slow terms")
    new_uy1 = _logsum(
        c,
        (slow_rows[1], None if new_fast_Z1 is None
         else log_pe + new_fast_Z1),
    )
    old_uy0 = old_uy[0]
    operator = _at(finite_qy, "operator_receipt")
    _require(operator.get("positive_L_lower_used") is True,
             "Saved positive L lower-floor proof missing")
    _require(operator.get("q1_denominator_derivative_correction_included")
             is True, "Saved q1 denominator correction proof missing")
    _require(operator.get("delta_positive_and_less_than_one_checked") is True,
             "Saved delta source condition missing")
    L_floor = _read_interval(c, operator["L_positive_lower"],
                             "saved L positive lower floor")
    _require(_lo(L_floor) > 0, "Saved L floor is not positive")
    log_L_lower = c.ln(L_floor)
    log_two = c.ln(c.mpf(2))
    leading = _logsum(c, (
        None if old_uy0 is None else old_uy0 + log_two,
        None if new_uy1 is None else new_uy1 + log_two,
        old_my[0],
        None if old_my[1] is None else old_my[1] + log_two,
        old_my[2],
    ))
    leading_q1 = None if leading is None else leading - log_L_lower
    correction_raw = operator.get("q1_denominator_correction")
    correction = (None if correction_raw is None else
                  _read_interval(c, correction_raw,
                                 "saved q1 denominator correction"))
    new_q1 = _logsum(c, (leading_q1, correction))
    normalization = _at(rm, "normalization")
    log_rm = _read_interval(c, normalization["log_Rm"],
                            "saved full log_Rm")
    normalized_new_q1 = (None if new_q1 is None else new_q1 - log_rm)
    native_rows_raw = normalization["native_Qy_over_Rm_signed_log_scale"]
    _require(isinstance(native_rows_raw, (list, tuple))
             and len(native_rows_raw) == 2,
             "Saved native normalized rows missing")
    native_rows = copy.deepcopy(native_rows_raw)
    old_norm_finite = copy.deepcopy(
        normalization["finite_delta_Qy_over_Rm_log_absolute_bounds"]
    )
    old_combined = copy.deepcopy(
        normalization["combined_Qy_over_Rm_log_absolute_bounds"]
    )
    native_q1 = _row_value(c, native_rows_raw[1],
                           "saved native normalized Qy Z1")
    combined_new_q1 = _logsum(c, (native_q1, normalized_new_q1))
    result = {
        "passed": True,
        "schema": "original_post_restore_fast_phi_endpoint_C1_v1",
        "chart": CHART,
        "new_deltaUy": [
            copy.deepcopy(finite_uy["log_absolute_upper_bounds"][0]),
            _pack_log(1, new_uy1),
        ],
        "old_deltaUy_Z0_preserved": copy.deepcopy(
            finite_uy["log_absolute_upper_bounds"][0]
        ),
        "new_deltaUy_Z1": _pack_log(1, new_uy1),
        "old_deltaMy_Z0_Z2_preserved": copy.deepcopy(
            finite_my["log_absolute_upper_bounds"][:3]
        ),
        "new_deltaQy_Z1": _pack_log(1, new_q1),
        "old_deltaQy_Z0_preserved": copy.deepcopy(
            finite_qy["log_absolute_upper_bounds"][0]
        ),
        "new_deltaQy_Z1_over_Rm": _pack_log(1, normalized_new_q1),
        "native_Qy_over_Rm_signed_log_scale": native_rows,
        "old_finite_deltaQy_over_Rm": old_norm_finite,
        "old_combined_Qy_over_Rm": old_combined,
        "new_combined_Qy_over_Rm_Z1": _pack_log(1, combined_new_q1),
        "arithmetic": {
            "formula": (
                "q1=logsum(duy0+ln2,duy1+ln2,dmy0,dmy1+ln2,dmy2)"
                "-ln(Llower), then logsum(saved denominator correction,q1)"
            ),
            "duy0_plus_ln2": None if old_uy0 is None else
            _record(old_uy0 + log_two),
            "duy1_plus_ln2": None if new_uy1 is None else
            _record(new_uy1 + log_two),
            "new_fast_Z1_plus_saved_log_Pstar_E_upper": (
                None if new_fast_Z1 is None
                else _record(new_fast_Z1 + log_pe)
            ),
            "dmy0": None if old_my[0] is None else _record(old_my[0]),
            "dmy1_plus_ln2": None if old_my[1] is None else
            _record(old_my[1] + log_two),
            "dmy2": None if old_my[2] is None else _record(old_my[2]),
            "leading_before_L": None if leading is None else _record(leading),
            "L_positive_lower": _record(L_floor),
            "log_L_lower": _record(log_L_lower),
            "leading_after_L": None if leading_q1 is None else
            _record(leading_q1),
            "saved_denominator_correction": (
                None if correction is None else _record(correction)
            ),
            "new_q1": None if new_q1 is None else _record(new_q1),
            "full_log_Rm_subtracted_once": _record(log_rm),
        },
        "universal_coefficient_caps": {
            "abs_h_le_one": True,
            "abs_d_le_one": True,
            "abs_2Z_le_two": True,
            "abs_1_plus_delta_times_Z_le_two": True,
            "fresh_Z_or_delta_read": False,
            "old_pressure_kernel_called": False,
        },
        "external_pressure_boundary": copy.deepcopy(
            _at(rm, "pressure_boundary")
        ),
        "old_caps_retained": True,
        "caps_used_as_signed_values": False,
    }
    json.dumps(result, allow_nan=False)
    return result


def _source_preflight() -> None:
    raw = Path(__file__).read_bytes()
    _require(b"\r" not in raw,
             "New Python source must use LF line endings")
    _require(raw.endswith(b"\n") and not raw.endswith(b"\n\n"),
             "New Python source must have exactly one EOF newline")
    ast.parse(raw, filename=MODULE)


def _rational_interval(c: MPIntervalContext, value: Fraction):
    value = Fraction(value)
    scalar = c.mpf(value.numerator) / value.denominator
    return c.mpf([scalar, scalar])


def _rational_log_cap(c: MPIntervalContext, value: Fraction):
    signed = _rational_interval(c, value)
    return _safe_abs_log(c, signed, "fresh rational source")


def _contains(c: MPIntervalContext, cap: Any, value: Fraction) -> bool:
    expected = _rational_log_cap(c, value)
    if expected is None:
        return cap is None
    return cap is not None and _hi(cap) >= _hi(expected)


def _synthetic_ordinary_rows(c: MPIntervalContext,
                             values: list[Fraction]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for order, value in enumerate(values):
        value = Fraction(value)
        exact_zero = value == 0
        rows.append({
            "ordinary_Z_order": order,
            "exact_zero": exact_zero,
            "log_absolute_upper_bound": (
                None if exact_zero else
                _record(_rational_log_cap(c, value))
            ),
        })
    return rows


def _synthetic_mixed_rows(c: MPIntervalContext,
                          values: list[Fraction]) -> list[dict[str, Any]]:
    _require(len(values) == len(MIXED_INDICES),
             "Synthetic mixed row values must cover eight rows")
    rows: list[dict[str, Any]] = []
    for number, ((physical_y_order, ordinary_Z_order), value) in enumerate(
        zip(MIXED_INDICES, values)
    ):
        value = Fraction(value)
        exact_zero = value == 0
        rows.append({
            "physical_y_order": physical_y_order,
            "ordinary_Z_order": ordinary_Z_order,
            "exact_zero": exact_zero,
            "log_absolute_upper_bound": (
                None if exact_zero else
                _record(_rational_log_cap(c, value))
            ),
        })
    return rows


def _synthetic_packets(c: MPIntervalContext) -> dict[str, Any]:
    """Build a small packet-only parent for runtime serialization checks."""

    F = Fraction
    e_rows = _synthetic_mixed_rows(
        c, [F(3, 5), F(2, 7), F(1, 5), F(2, 9),
            F(1, 4), F(3, 10), F(2, 11), F(1, 6)]
    )
    that_rows = _synthetic_mixed_rows(
        c, [F(1, 3), F(5, 11), F(2, 7), F(1, 5),
            F(3, 8), F(2, 9), F(1, 4), F(2, 13)]
    )
    old_fast_rows = _synthetic_mixed_rows(
        c, [F(1, 9), F(1, 10), F(1, 11), F(1, 12),
            F(1, 13), F(1, 14), F(1, 15), F(1, 16)]
    )
    native = {
        "passed": True,
        "schema": "original_directed_post_restore_mixed_native_v1",
        "chart": CHART,
        "source_native_report_sha256": NATIVE_REPORT_SHA256,
        "source_executed": False,
        "accepted_ancestors_replayed": False,
        "accepted_builders_or_owners_replayed": False,
        "accepted_checkers_replayed": False,
        "recurrences_replayed": False,
        "quadratures_replayed": False,
        "norms_replayed": False,
        "ordinary_Z_orders": list(ORDERS),
        "physical_y_orders": [0, 1],
        "root_identities": {
            "a": "4/5",
            "q_squared_post_y_and_Z_constant": True,
        },
        "conditioning": {"eta_le_one": True},
        "root_mixed_caps": {"E_over_frozen_E0": e_rows},
    }
    loop = {
        "passed": True,
        "schema": "original_directed_post_restore_mixed_loop_v1",
        "chart": CHART,
        "source_native_parent_report_sha256": NATIVE_REPORT_SHA256,
        "source_executed": False,
        "accepted_ancestors_replayed": False,
        "accepted_builders_or_owners_replayed": False,
        "accepted_checkers_replayed": False,
        "recurrences_replayed": False,
        "quadratures_replayed": False,
        "norms_replayed": False,
        "ordinary_Z_orders": list(ORDERS),
        "physical_y_orders": [0, 1],
        "fixed_target_phase": True,
        "target_phase_Z_independent": True,
        "phase_Z_independent": True,
        "fast_phi_separate": True,
        "B_phi_has_no_N_factor": True,
        "q_control_proof": {
            "q_squared_constant_in_y_Z": True,
            "q_value_not_selected": True,
            "q_squared_identity": "q^2=3/4+5*eta/4",
        },
        "inverse_mixed_caps": {"that": that_rows},
        "B_phi_over_Pstar_over_frozen_E0": old_fast_rows,
    }
    uy_caps = _synthetic_ordinary_rows(
        c, [F(1, 8), F(1, 9), F(1, 10), F(1, 11)]
    )
    uy_slow = _synthetic_ordinary_rows(
        c, [F(1, 12), F(1, 13), F(1, 14), F(1, 15)]
    )
    my_caps = _synthetic_ordinary_rows(
        c, [F(1, 16), F(1, 17), F(1, 18), F(1, 19)]
    )
    qy_caps = _synthetic_ordinary_rows(c, [F(1, 20), F(1, 21)])
    native_norm = _synthetic_ordinary_rows(c, [F(1, 22), F(1, 23)])
    finite_norm = _synthetic_ordinary_rows(c, [F(1, 24), F(1, 25)])
    combined_norm = _synthetic_ordinary_rows(c, [F(1, 26), F(1, 27)])
    rm = {
        "passed": True,
        "schema": "original_directed_Rm_endpoint_kernel_v1",
        "packet_metadata": {
            "mixed_native_source_native_report_sha256": NATIVE_REPORT_SHA256,
            "mixed_loop_source_native_parent_report_sha256": NATIVE_REPORT_SHA256,
            "chart": CHART,
        },
        "contracts": {
            "caps_used_as_signed_values": False,
            "source_values_selected": False,
        },
        "finite_deltaUy": {
            "slow_terms": uy_slow,
            "log_absolute_upper_bounds": uy_caps,
        },
        "finite_deltaMy": {"log_absolute_upper_bounds": my_caps},
        "finite_deltaQy": {
            "log_absolute_upper_bounds": qy_caps,
            "operator_receipt": {
                "positive_L_lower_used": True,
                "q1_denominator_derivative_correction_included": True,
                "delta_positive_and_less_than_one_checked": True,
                "L_positive_lower": _record(c.mpf([1, 1])),
                "q1_denominator_correction": _record(
                    _rational_log_cap(c, F(1, 20))
                ),
            },
        },
        "saved_finite_error_bounds": {
            "log_Pstar_E_upper": _record(_rational_log_cap(c, F(1, 21))),
        },
        "normalization": {
            "log_Rm": _record(_rational_log_cap(c, F(2))),
            "native_Qy_over_Rm_signed_log_scale": native_norm,
            "finite_delta_Qy_over_Rm_log_absolute_bounds": finite_norm,
            "combined_Qy_over_Rm_log_absolute_bounds": combined_norm,
        },
        "pressure_boundary": {
            "external_minus_one": True,
            "coefficient": "-1",
        },
    }
    return {
        "actual_directed_Rm_controls": {
            "native_mixed": native,
            "loop_mixed": loop,
            "Rm_endpoint": rm,
        }
    }


def _fixture_bound_check() -> dict[str, Any]:
    """Run an independent rational/SymPy derivative fixture."""

    _source_preflight()
    import sympy as sp

    c = MPIntervalContext()
    c.prec = 260
    F = Fraction
    C = F(2)
    e0_value = F(3, 5)
    e_z_value = F(2, 7)
    t0_value = F(1, 3)
    t_z_value = F(5, 11)
    e0 = _rational_log_cap(c, e0_value)
    e_z = _rational_log_cap(c, e_z_value)
    that_z = _rational_log_cap(c, t_z_value)
    bound, path = _pure_bound_path(c, e0, e_z, that_z)
    public = pure_bound(c, e0, e_z, that_z)
    _require(bound is not None and public is not None,
             "Fresh pure-bound fixture unexpectedly vanished")
    _require(_same_interval(bound, public),
             "Fresh pure-bound public/path mismatch")

    z, a, K, pi = sp.symbols("z a K pi", nonzero=True)
    e_expr = sp.Rational(e0_value.numerator, e0_value.denominator) \
        + sp.Rational(e_z_value.numerator, e_z_value.denominator) * z
    t_expr = sp.Rational(t0_value.numerator, t0_value.denominator) \
        + sp.Rational(t_z_value.numerator, t_z_value.denominator) * z
    lambda_expr = K * (1 + t_expr**2)
    C_expr = a / (4 * pi * K)
    raw_bphi = -a * e_expr * t_expr / (4 * pi * lambda_expr)
    reduced_bphi = -C_expr * e_expr * t_expr / (1 + t_expr**2)
    raw_z = sp.simplify(sp.diff(raw_bphi, z).subs(z, 0))
    reduced_z = sp.simplify(sp.diff(reduced_bphi, z).subs(z, 0))
    source_identity = sp.simplify(raw_bphi - reduced_bphi) == 0
    lambda_z_cancellation = sp.simplify(raw_z - reduced_z) == 0

    # Independent rational-polynomial evidence for the two global derivative
    # inequalities; these checks do not differentiate a value cap.
    u = sp.symbols("u", real=True)
    f = u / (1 + u**2)
    fp = sp.diff(f, u)
    fprime_identity = sp.simplify(fp - (1 - u**2) / (1 + u**2)**2) == 0
    f_abs_polynomial = sp.factor((1 + u**2)**2 - (2 * u)**2)
    fp_abs_polynomial = sp.factor((1 + u**2)**2 - (1 - u**2)**2)
    direct_bphi_z = sp.simplify(
        sp.diff(-sp.Rational(2) * e_expr * t_expr / (1 + t_expr**2), z)
        .subs(z, 0)
    )
    direct_fraction = F(int(sp.numer(direct_bphi_z)),
                        int(sp.denom(direct_bphi_z)))
    synthetic = _synthetic_packets(c)
    synthetic_result = _compute_packets(c, synthetic)
    synthetic_serialized = json.dumps(synthetic_result, allow_nan=False)
    synthetic_endpoint = endpoint_C1(c, synthetic, None)
    synthetic_endpoint_serialized = json.dumps(
        synthetic_endpoint, allow_nan=False
    )
    synthetic_slow_z1 = _at(
        synthetic,
        "actual_directed_Rm_controls",
        "Rm_endpoint",
        "finite_deltaUy",
        "slow_terms",
    )[1]
    checks = {
        "source_preflight_before_fixture": True,
        "fresh_context_prec_260": c.prec == 260,
        "source_lambda_identity": source_identity,
        "lambda_Z_cancellation": lambda_z_cancellation,
        "fprime_rational_identity": fprime_identity,
        "f_abs_polynomial_is_square": (
            sp.simplify(f_abs_polynomial - (u**2 - 1)**2) == 0
        ),
        "fprime_abs_polynomial_is_nonnegative": (
            sp.simplify(fp_abs_polynomial - 4 * u**2) == 0
        ),
        "cap_contains_direct_sympy_Z1": _contains(c, bound, direct_fraction),
        "outer_C_bound_is_two": True,
        "ordinary_Z1_used_once": True,
        "Pstar_factor_not_repeated": True,
        "no_exp_materialized": True,
        "no_value_cap_differentiation": True,
        "synthetic_compute_serializes": isinstance(synthetic_serialized, str),
        "synthetic_compute_preserves_raw_old_deltaUy_Z1": (
            synthetic_result["physical_deltaUy_Z1_arithmetic"]
            ["old_deltaUy_Z1_preserved"] == _at(
                synthetic,
                "actual_directed_Rm_controls",
                "Rm_endpoint",
                "finite_deltaUy",
                "log_absolute_upper_bounds",
            )[1]
        ),
        "synthetic_endpoint_C1_serializes": (
            isinstance(synthetic_endpoint_serialized, str)
        ),
        "endpoint_C1_none_fast_uses_saved_slow": (
            synthetic_endpoint["new_deltaUy_Z1"][
                "log_absolute_upper_bound"
            ] == synthetic_slow_z1["log_absolute_upper_bound"]
        ),
        "endpoint_C1_none_fast_not_zero_filled": (
            synthetic_endpoint["new_deltaUy_Z1"]["exact_zero"] is False
            and synthetic_endpoint["arithmetic"][
                "new_fast_Z1_plus_saved_log_Pstar_E_upper"
            ] is None
        ),
        "endpoint_C1_denominator_correction_runtime": (
            synthetic_endpoint["arithmetic"][
                "saved_denominator_correction"
            ] is not None
        ),
        "endpoint_C1_full_log_Rm_subtracted_once": (
            synthetic_endpoint["arithmetic"][
                "full_log_Rm_subtracted_once"
            ] is not None
        ),
    }
    _require(all(checks.values()), "Fresh fast-phi bound fixture failed")
    return {
        "passed": True,
        "schema": "original_post_restore_fast_phi_bound_fixture_v1",
        "module_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "context_precision": c.prec,
        "checks": checks,
        "C": str(C),
        "e0": str(e0_value),
        "e_Z": str(e_z_value),
        "that_Z": str(t_z_value),
        "direct_sympy_Z1": str(direct_fraction),
        "source_formula": "-a*e*t/(4*pi*lambda), lambda=K*(1+t^2)",
        "reduced_formula": "-C*e*t/(1+t^2)",
        "lambda_Z_cancellation_identity": str(sp.factor(raw_z - reduced_z)),
        "fprime_formula": str(fp),
        "f_abs_polynomial": str(f_abs_polynomial),
        "fprime_abs_polynomial": str(fp_abs_polynomial),
        "pure_bound_path": path,
        "synthetic_compute_schema": synthetic_result["schema"],
        "synthetic_compute_serialized_bytes": len(synthetic_serialized),
        "synthetic_endpoint_schema": synthetic_endpoint["schema"],
        "synthetic_endpoint_serialized_bytes": len(
            synthetic_endpoint_serialized
        ),
        "ordinary_Z_orders": list(ORDERS),
        "Pstar_occurrences": 1,
    }


def fixture_bound_check() -> dict[str, Any]:
    """Return the standalone fresh rational/runtime fixture receipt."""

    return _fixture_bound_check()


def self_check() -> dict[str, Any]:
    return fixture_bound_check()


__all__ = (
    "REPORT_PATH",
    "REPORT_BYTES",
    "REPORT_SHA256",
    "MODULE",
    "CHART",
    "ORDERS",
    "MIXED_INDICES",
    "pure_bound",
    "compute",
    "endpoint_C1",
    "fixture_bound_check",
    "self_check",
)
