"""Saved-data adapter for the post-restore ``log |u|`` endpoint.

This module consumes authenticated saved packets only.  The changed fast
ordinary-``Z`` row uses the signed native ``p2_Z/p2`` record and the uniform
source estimate ``|W| <= 600``.  The old absolute caps remain inputs to the
finite-``N`` endpoint operator.  No producer, builder, checker, recurrence,
quadrature, norm, or accepted mathematical module is imported or executed.

The two source steps are, in logarithmic absolute-cap form,

``new_fast_Z1 = ln(2) + logsum(e_Z-ln(2), e_0+ln(600)+ln|u_Z/u|)``

and

``Q_y,Z1 = logsum(deltaUy,Z0+ln(2), deltaUy,Z1+ln(2),
                  deltaMy,Z0, deltaMy,Z1+ln(2), deltaMy,Z2)
            - ln(L_lower) + denominator_correction``.

The physical ``Pstar E0`` record is inserted once into ``deltaUy,Z1`` and
the full saved ``log Rm`` record is subtracted once during normalization.
Large logarithms stay in interval arithmetic; no exponential materialization
is used.
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


BASE = "experiments/root_st073/"
SIGNED_REPORT_PATH = BASE + (
    "lei_ren_part1_paper_compliant_current_original_signed_native_sources.json"
)
SIGNED_REPORT_BYTES = 566530
SIGNED_REPORT_SHA256 = (
    "9887a353aeab1b24f301a5b6c240f83f58a5d54d001eaac6d16f932b87979561"
)
FAST_REPORT_PATH = BASE + (
    "lei_ren_part1_paper_compliant_current_original_fast_phi_correlated_C1.json"
)
FAST_REPORT_BYTES = 497555
FAST_REPORT_SHA256 = (
    "8b77d380fcfd4837c4bcfff10a1e4fe512b37d47b30340c5093552bb3b1ca108"
)
RM_REPORT_PATH = BASE + (
    "lei_ren_part1_paper_compliant_current_original_directed_Rm_controls.json"
)
RM_REPORT_BYTES = 3159660
RM_REPORT_SHA256 = (
    "78a95b631bdce6dbe1438bed1ccff282127df5d2198e6836d250a2b9ff325deb"
)
NATIVE_REPORT_PATH = BASE + (
    "lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json"
)
NATIVE_REPORT_BYTES = 42177499
NATIVE_REPORT_SHA256 = (
    "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"
)
MODULE = "src/openai_ns_reconstruction/original_post_restore_logu_endpoint.py"
CHART = "post_restore"
ORDERS = (0, 1, 2, 3)
MIXED_INDICES = tuple((i, k) for i in (0, 1) for k in ORDERS)


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


def _record(value: Any) -> dict[str, Any]:
    _require(hasattr(value, "_mpi_"), "interval record requires MP interval")
    lower = tuple(value._mpi_[0])
    upper = tuple(value._mpi_[1])
    _canonical_endpoint(lower, "interval.lower")
    _canonical_endpoint(upper, "interval.upper")
    lo, hi = make_mpf(lower), make_mpf(upper)
    _require(lo <= hi and isfinite(lo) and isfinite(hi),
             "interval must be finite and ordered")
    return {
        "lower_exact_mpf_tuple": list(lower),
        "upper_exact_mpf_tuple": list(upper),
        "lower_display": nstr(lo, 24),
        "upper_display": nstr(hi, 24),
        "directed": True,
    }


def _read_interval(c: MPIntervalContext, value: Any, label: str):
    """Decode a directed record without selecting an endpoint."""

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


def _read_signed_interval(c: MPIntervalContext, value: Any, label: str):
    """Decode an explicit signed source enclosure, rejecting log caps."""

    if isinstance(value, Mapping):
        _require("log_absolute_upper_bound" not in value and
                 "log_cap" not in value and "cap" not in value,
                 label + " is a logarithmic cap, not a signed source")
    if isinstance(value, Mapping) and "value" in value \
            and "lower_exact_mpf_tuple" not in value:
        kind = value.get("kind")
        if kind is not None:
            _require(kind == "signed_function_enclosure_not_selected_value",
                     label + " has an unknown signed-source kind")
        value = value["value"]
    _require(isinstance(value, Mapping), label + " must be a signed mapping")
    _require("log_absolute_upper_bound" not in value and
             "log_cap" not in value and "cap" not in value,
             label + " is a logarithmic cap, not a signed source")
    return _read_interval(c, value, label)


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _admit_report(raw: bytes, expected_bytes: int, expected_sha: str,
                  label: str) -> Mapping[str, Any]:
    _require(type(raw) is bytes, label + " bytes must be immutable bytes")
    _require(len(raw) == expected_bytes,
             label + " byte length changed")
    _require(_sha(raw) == expected_sha,
             label + " SHA256 changed")
    try:
        report = json.loads(raw)
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError(label + " is not valid JSON") from exc
    _require(isinstance(report, Mapping), label + " must be a mapping")
    return report


def _at(mapping: Mapping[str, Any], *keys: str) -> Any:
    value: Any = mapping
    for key in keys:
        _require(isinstance(value, Mapping) and key in value,
                 "Missing saved path: " + ".".join(keys))
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
    if isinstance(raw, Mapping):
        for key in ("rows", "mixed_log_caps", "log_caps", "values"):
            if key in raw:
                return _mixed_rows(c, raw[key], label)
        if set(raw) == set(MIXED_INDICES):
            return {
                index: (None if raw[index] is None else
                        _read_interval(c, raw[index], label + str(index)))
                for index in MIXED_INDICES
            }
        if set(raw) >= {"y0", "y1"}:
            raw = list(raw["y0"]) + list(raw["y1"])
    _require(isinstance(raw, (list, tuple)) and len(raw) == 8,
             label + " must contain exactly eight rows")
    out: dict[tuple[int, int], Any] = {}
    for number, expected in enumerate(MIXED_INDICES):
        row = raw[number]
        _require(isinstance(row, Mapping), label + " row is malformed")
        index = (row.get("physical_y_order", row.get("y_order")),
                 row.get("ordinary_Z_order", row.get("Z_order")))
        _require(index == expected, label + " is not in canonical y/Z order")
        out[index] = _row_value(c, row, label + "[" + str(number) + "]")
    return out


def _ordinary_rows(c: MPIntervalContext, raw: Any, count: int, label: str):
    _require(isinstance(raw, (list, tuple)), label + " must be a row list")
    out: dict[int, Any] = {}
    for number, row in enumerate(raw):
        _require(isinstance(row, Mapping), label + " row is malformed")
        order = row.get("ordinary_Z_order", row.get("Z_order", number))
        _require(order == number and order in range(count),
                 label + " is not in canonical ordinary-Z order")
        _require(order not in out, label + " repeats an ordinary-Z row")
        out[order] = _row_value(c, row, label + "[" + str(number) + "]")
    _require(len(out) == count,
             label + " must cover ordinary Z orders 0 through " + str(count - 1))
    return out


def _logsum(c: MPIntervalContext, terms: Any):
    """Conservative log-sum without exponentiating any term."""

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
    value = _read_interval(c, value, label)
    magnitude = c.absmax(value)
    upper = _hi(magnitude)
    if upper == 0:
        return None
    return c.ln(c.mpf([upper, upper]))


def _signed_abs_log(c: MPIntervalContext, value: Any, label: str):
    value = _read_signed_interval(c, value, label)
    magnitude = c.absmax(value)
    upper = _hi(magnitude)
    if upper == 0:
        return None
    return c.ln(c.mpf([upper, upper]))


def _pack_log(order: int, value: Any) -> dict[str, Any]:
    return {
        "ordinary_Z_order": order,
        "exact_zero": value is None,
        "log_absolute_upper_bound": None if value is None else _record(value),
    }


def _same_interval(left: Any, right: Any) -> bool:
    return (left is None and right is None) or (
        left is not None and right is not None and
        tuple(left._mpi_[0]) == tuple(right._mpi_[0]) and
        tuple(left._mpi_[1]) == tuple(right._mpi_[1])
    )


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
        _require(packet.get("chart") == CHART, label + " chart changed")
        for key in (
            "source_executed", "accepted_ancestors_replayed",
            "accepted_builders_or_owners_replayed", "accepted_checkers_replayed",
            "recurrences_replayed", "quadratures_replayed", "norms_replayed",
        ):
            if key in packet:
                _require(packet[key] is False,
                         label + " execution flag changed: " + key)
    _require(native.get("source_native_report_sha256") == NATIVE_REPORT_SHA256,
             "Native source SHA binding changed")
    _require(loop.get("source_native_parent_report_sha256") == NATIVE_REPORT_SHA256,
             "Loop source SHA binding changed")
    metadata = rm.get("packet_metadata")
    _require(isinstance(metadata, Mapping), "Rm packet metadata missing")
    _require(metadata.get("mixed_native_source_native_report_sha256") ==
             NATIVE_REPORT_SHA256 and
             metadata.get("mixed_loop_source_native_parent_report_sha256") ==
             NATIVE_REPORT_SHA256,
             "Rm native-source metadata binding changed")
    _require(metadata.get("chart") == CHART, "Rm metadata chart changed")
    roots = native.get("root_identities")
    _require(isinstance(roots, Mapping), "Native root identities missing")
    _require(roots.get("a") == "4/5" and
             roots.get("q_squared_post_y_and_Z_constant") is True,
             "Native endpoint identities changed")
    _require(native.get("ordinary_Z_orders") == list(ORDERS) and
             native.get("physical_y_orders") == [0, 1],
             "Native row orders changed")
    _require(loop.get("ordinary_Z_orders") == list(ORDERS) and
             loop.get("physical_y_orders") == [0, 1],
             "Loop row orders changed")
    proof = loop.get("q_control_proof")
    _require(isinstance(proof, Mapping) and
             proof.get("q_squared_constant_in_y_Z") is True and
             proof.get("q_value_not_selected") is True and
             proof.get("q_squared_identity") == "q^2=3/4+5*eta/4",
             "Saved q source-correlation proof changed")
    contracts = rm.get("contracts")
    _require(isinstance(contracts, Mapping) and
             contracts.get("caps_used_as_signed_values") is False and
             contracts.get("source_values_selected") is False,
             "Saved cap/sign contract changed")
    return {
        "native_report_sha256": NATIVE_REPORT_SHA256,
        "native_a_identity": roots["a"],
        "q_squared_identity": proof["q_squared_identity"],
        "q_squared_constant_in_y_Z": True,
        "q_source_value_not_selected": True,
        "source_sha_binding_verified": True,
        "caps_remain_absolute": True,
        "signed_source_values_selected": False,
    }


def _signed_binding(report: Mapping[str, Any]) -> Mapping[str, Any]:
    _require(report.get("schema") == "original_signed_native_sources_v1",
             "Signed-native schema changed")
    _require(report.get("native_report_sha256") == NATIVE_REPORT_SHA256,
             "Signed-native native report SHA changed")
    actual = report.get("actual_signed_native_sources")
    _require(isinstance(actual, Mapping), "Signed-native payload missing")
    for name in ("gamma2", "histories", "p2"):
        packet = actual.get(name)
        _require(isinstance(packet, Mapping), "Signed-native component missing: " + name)
        _require(packet.get("accepted_ancestors_replayed") is False,
                 "Signed-native ancestor replay flag changed: " + name)
        if name != "gamma2":
            _require(packet.get("payload", {}).get("passed") is True,
                     "Signed-native component failed: " + name)
    p2 = actual["p2"]["payload"]
    _require(p2.get("full_native_p2_signed_scaled_enclosure") is True and
             p2.get("physical_R_applied_once") is True and
             p2.get("midpoint_or_logcap_promoted") is False,
             "Signed p2 source contract changed")
    return actual


def _fast_binding(report: Mapping[str, Any]) -> Mapping[str, Any]:
    actual = report.get("actual_fast_phi_correlated_C1")
    _require(isinstance(actual, Mapping), "Fast correlated packet missing")
    fast = actual.get("fast")
    _require(isinstance(fast, Mapping) and fast.get("passed") is True,
             "Saved fast packet failed")
    _require(fast.get("schema") == "original_post_restore_fast_phi_bound_v1",
             "Saved fast packet schema changed")
    _require(fast.get("input_report_sha256") == RM_REPORT_SHA256,
             "Fast packet Rm source SHA changed")
    _require(fast.get("contracts", {}).get("caps_used_as_signed_values") is False,
             "Fast cap/sign contract changed")
    endpoint = actual.get("endpoint_C1")
    _require(isinstance(endpoint, Mapping) and endpoint.get("passed") is True,
             "Saved failed endpoint packet missing")
    return actual


def _native_binding(report: Mapping[str, Any]) -> Mapping[str, Any]:
    actual = report.get("actual_native_directed_controls_v2")
    _require(isinstance(actual, Mapping), "Native-v2 packet missing")
    _require(actual.get("schema") == "original_native_directed_controls_v2_v1" and
             actual.get("passed") is True,
             "Native-v2 packet failed")
    checks = actual.get("checks", {})
    _require(checks.get("accepted_ancestors_replayed") is False and
             checks.get("source_executed") is False,
             "Native-v2 execution contract changed")
    _require(actual.get("scope", {}).get("corrected_native_controls_verified") is True,
             "Native-v2 corrected-control receipt missing")
    return actual


def _source_rows(parent: Mapping[str, Any], c: MPIntervalContext):
    packets = parent.get("actual_directed_Rm_controls", parent)
    _require(isinstance(packets, Mapping), "Saved Rm packet missing")
    native = packets.get("native_mixed")
    loop = packets.get("loop_mixed")
    rm = packets.get("Rm_endpoint")
    _require(isinstance(native, Mapping) and isinstance(loop, Mapping) and
             isinstance(rm, Mapping), "Saved native/loop/Rm packets missing")
    binding = _source_binding(native, loop, rm)
    e_rows = _mixed_rows(c, _at(native, "root_mixed_caps", "E_over_frozen_E0"),
                         "native E_over_frozen_E0")
    old_fast_rows = _at(loop, "B_phi_over_Pstar_over_frozen_E0")
    _mixed_rows(c, old_fast_rows, "saved old B_phi")
    return native, loop, rm, binding, e_rows, old_fast_rows


def _pure_fast_z1_bound(c: MPIntervalContext, e0_log: Any,
                        e_z_log: Any, ratio_log: Any,
                        log_W_upper: Any | None = None):
    e0_log = None if e0_log is None else _read_interval(c, e0_log, "e0 log cap")
    e_z_log = None if e_z_log is None else _read_interval(c, e_z_log, "e_Z log cap")
    ratio_log = None if ratio_log is None else _read_interval(c, ratio_log, "log|u_Z/u|")
    log_two = c.ln(c.mpf(2))
    if log_W_upper is None:
        log_W_upper = c.ln(c.mpf(600))
    else:
        log_W_upper = _read_interval(c, log_W_upper, "log W cap")
    first = None if e_z_log is None else e_z_log - log_two
    second = None if e0_log is None or ratio_log is None \
        else e0_log + log_W_upper + ratio_log
    inner = _logsum(c, (first, second))
    return None if inner is None else log_two + inner


def pure_bound(c: MPIntervalContext, e0_log: Any, e_z_log: Any,
               ratio_log: Any, log_W_upper: Any | None = None):
    """Return the new source-correlated fast ``Z1`` logarithmic cap."""

    _require(isinstance(c, MPIntervalContext), "c must be an MPIntervalContext")
    return _pure_fast_z1_bound(c, e0_log, e_z_log, ratio_log, log_W_upper)


def _pure_bound_path(c: MPIntervalContext, e0_log: Any, e_z_log: Any,
                     ratio_log: Any):
    log_two = c.ln(c.mpf(2))
    log_600 = c.ln(c.mpf(600))
    e0 = None if e0_log is None else _read_interval(c, e0_log, "e0 log cap")
    ez = None if e_z_log is None else _read_interval(c, e_z_log, "e_Z log cap")
    ratio = None if ratio_log is None else _read_interval(c, ratio_log, "log|u_Z/u|")
    first = None if ez is None else ez - log_two
    second = None if e0 is None or ratio is None \
        else e0 + log_600 + ratio
    inner = _logsum(c, (first, second))
    result = None if inner is None else log_two + inner
    return result, {
        "e_Z_minus_ln2": None if first is None else _record(first),
        "e0_plus_ln600_plus_logabs_uZ_over_u": (
            None if second is None else _record(second)
        ),
        "inner_logsum": None if inner is None else _record(inner),
        "outer_ln2": _record(log_two),
        "uniform_W_upper": _record(log_600),
        "result": None if result is None else _record(result),
        "logsum_rule": "interval term+ln(count), then min/max rebox",
    }


def _domain_receipt(c: MPIntervalContext, native_v2: Mapping[str, Any],
                    signed_actual: Mapping[str, Any],
                    rm_native: Mapping[str, Any]):
    raw = _at(native_v2, "raw_interval_tree")
    shared = _at(raw, "shared")
    restoration = _at(raw, "restoration", "charts", CHART)
    scales = _at(restoration, "scales")
    reference = _at(raw, "long_reference", "charts", "reference")
    reference_scale = _at(reference, "scale")
    logeta = _read_interval(c, shared["logeta"], "raw shared logeta")
    logdstar = _read_interval(c, shared["logdstar"], "raw shared logdstar")
    _require(_hi(logeta) < 0, "Saved logeta does not prove eta<1")
    z = _read_interval(c, shared["z"], "raw shared z")
    log_R = _read_interval(c, scales["log_R"], "post_restore log_R source")
    log_P = _read_interval(c, reference_scale["logP"], "reference logP source")
    offset = _read_interval(c, scales["offsets"], "post_restore offset source")
    _require(_lo(offset) == _lo(c.mpf(-7)) and
             _hi(offset) == _hi(c.mpf(-6)),
             "Saved post_restore offset domain changed")
    p2_chart = _at(signed_actual, "p2", "payload", "charts", CHART)
    _require(p2_chart.get("p2_value_sign") == "negative",
             "Signed native p2 sign changed")
    K_rows = p2_chart.get("K")
    _require(isinstance(K_rows, (list, tuple)) and len(K_rows) >= 1,
             "Signed native p2 K row missing")
    K0 = _read_signed_interval(c, K_rows[0], "signed native p2 K Z0")
    _require(_hi(K0) < 0, "Signed native p2 K Z0 is not negative")
    log_abs_K = c.ln(-K0)
    log_one_plus_z2 = c.ln(c.mpf(1) + z * z)
    log_two = c.ln(c.mpf(2))
    log_u_lower = (
        log_R + log_P + offset / c.mpf(10) - log_one_plus_z2
        + log_abs_K - log_two - log_two - logdstar
    )
    _require(_lo(log_u_lower) > 0,
             "Saved source logs do not certify |u|>=1")
    roots = rm_native.get("root_identities", {})
    _require(roots.get("u") == "p2*q/dstar" and
             roots.get("q_squared_post_y_and_Z_constant") is True,
             "Saved u/q source identity changed")
    q2_lower, q2_upper = c.mpf(3) / 4, c.mpf(2)
    q_floor, q_upper = c.mpf(1) / 2, c.mpf(3) / 2
    _require(q2_lower >= q_floor * q_floor and
             q2_upper <= q_upper * q_upper,
             "q source interval does not imply [1/2,3/2]")
    return {
        "certified_abs_u_ge_one": True,
        "log_u_lower": _record(log_u_lower),
        "log_u_lower_endpoint": _record(c.mpf([_lo(log_u_lower), _lo(log_u_lower)])),
        "log_R_source": _record(log_R),
        "log_P_source": _record(log_P),
        "offset_source": _record(offset),
        "z_source": _record(z),
        "log_one_plus_z_squared": _record(log_one_plus_z2),
        "signed_K_source": _record(K0),
        "log_abs_K": _record(log_abs_K),
        "log_dstar_source": _record(logdstar),
        "q_squared_interval": ["3/4", "2"],
        "q_interval": ["1/2", "3/2"],
        "q_floor_is_bound_not_selected_value": True,
        "eta_upper_log_lt_zero": True,
        "formula": (
            "log|u|=logR+logP+offset/10-log(1+Z^2)+log|K|"
            "-ln2-ln2-logdstar"
        ),
        "source_paths": {
            "log_R": "actual_native_directed_controls_v2.raw_interval_tree.restoration.charts.post_restore.scales.log_R",
            "log_P": "actual_native_directed_controls_v2.raw_interval_tree.long_reference.charts.reference.scale.logP",
            "offset": "actual_native_directed_controls_v2.raw_interval_tree.restoration.charts.post_restore.scales.offsets",
            "logeta": "actual_native_directed_controls_v2.raw_interval_tree.shared.logeta",
            "logdstar": "actual_native_directed_controls_v2.raw_interval_tree.shared.logdstar",
            "p2_K": "actual_signed_native_sources.p2.payload.charts.post_restore.K[0].value",
        },
        "stored_logs_not_exponentiated": True,
    }


def _finite_delta_uy_inputs(c: MPIntervalContext, rm: Mapping[str, Any],
                            new_fast: Any):
    finite = _at(rm, "finite_deltaUy")
    old_slow = _ordinary_rows(c, finite["slow_terms"], 4,
                              "saved deltaUy slow terms")
    old_delta = _ordinary_rows(c, finite["log_absolute_upper_bounds"], 4,
                               "saved deltaUy caps")
    log_pe = _read_interval(c,
                            _at(rm, "saved_finite_error_bounds",
                                "log_Pstar_E_upper"),
                            "saved physical Pstar E0")
    new_uy1 = _logsum(c, (old_slow[1], None if new_fast is None
                          else log_pe + new_fast))
    return old_slow, old_delta, log_pe, new_uy1


def _parent_packets(parent: Mapping[str, Any]):
    packets = parent.get("actual_directed_Rm_controls", parent)
    _require(isinstance(packets, Mapping), "Parent Rm packets missing")
    native, loop, rm = packets.get("native_mixed"), packets.get("loop_mixed"), packets.get("Rm_endpoint")
    _require(isinstance(native, Mapping) and isinstance(loop, Mapping) and
             isinstance(rm, Mapping), "Parent native/loop/Rm packets missing")
    _source_binding(native, loop, rm)
    return native, loop, rm


def _endpoint_packets(c: MPIntervalContext, parent: Mapping[str, Any],
                      new_fast_Z1: Any) -> dict[str, Any]:
    native, loop, rm = _parent_packets(parent)
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
    old_slow, old_delta, log_pe, new_uy1 = _finite_delta_uy_inputs(
        c, rm, new_fast_Z1
    )
    operator = _at(finite_qy, "operator_receipt")
    _require(operator.get("positive_L_lower_used") is True and
             operator.get("q1_denominator_derivative_correction_included") is True and
             operator.get("delta_positive_and_less_than_one_checked") is True,
             "Saved denominator receipt changed")
    L_floor = _read_interval(c, operator["L_positive_lower"], "saved L floor")
    _require(_lo(L_floor) > 0, "Saved L lower floor is not positive")
    log_L_lower = c.ln(c.mpf([_lo(L_floor), _lo(L_floor)]))
    log_two = c.ln(c.mpf(2))
    leading = _logsum(c, (
        None if old_uy[0] is None else old_uy[0] + log_two,
        None if new_uy1 is None else new_uy1 + log_two,
        old_my[0],
        None if old_my[1] is None else old_my[1] + log_two,
        old_my[2],
    ))
    leading_q1 = None if leading is None else leading - log_L_lower
    correction = _read_interval(c, operator["q1_denominator_correction"],
                                "saved denominator correction")
    new_q1 = _logsum(c, (leading_q1, correction))
    normalization = _at(rm, "normalization")
    log_rm = _read_interval(c, normalization["log_Rm"], "saved full log_Rm")
    normalized_new_q1 = None if new_q1 is None else new_q1 - log_rm
    native_rows_raw = normalization["native_Qy_over_Rm_signed_log_scale"]
    _require(isinstance(native_rows_raw, (list, tuple)) and len(native_rows_raw) == 2,
             "Saved native normalized Qy rows missing")
    native_q1 = _row_value(c, native_rows_raw[1], "saved native normalized Qy Z1")
    old_norm_finite = copy.deepcopy(
        normalization["finite_delta_Qy_over_Rm_log_absolute_bounds"]
    )
    old_combined = copy.deepcopy(
        normalization["combined_Qy_over_Rm_log_absolute_bounds"]
    )
    combined_new_q1 = _logsum(c, (native_q1, normalized_new_q1))
    old_combined_rows = _ordinary_rows(c, old_combined, 2,
                                       "old combined normalized Qy")
    threshold = -c.mpf(1024) * log_two
    both_rows = {0: old_combined_rows[0], 1: combined_new_q1}
    threshold_passed = all(value is not None and _hi(value) < _lo(threshold)
                           for value in both_rows.values())
    return {
        "passed": True,
        "schema": "original_post_restore_logu_endpoint_C1_v1",
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
        "native_Qy_over_Rm_signed_log_scale": copy.deepcopy(native_rows_raw),
        "native_Qy_over_Rm_Z1": copy.deepcopy(native_rows_raw[1]),
        "old_finite_deltaQy_over_Rm": old_norm_finite,
        "old_combined_Qy_over_Rm": old_combined,
        "new_combined_Qy_over_Rm_Z1": _pack_log(1, combined_new_q1),
        "both_normalized_rows": {
            "old_Z0": _pack_log(0, old_combined_rows[0]),
            "new_Z1": _pack_log(1, combined_new_q1),
        },
        "threshold": _record(threshold),
        "threshold_passed": threshold_passed,
        "arithmetic": {
            "formula": (
                "duy1=logsum(saved slow_terms[1],saved log_Pstar_E_upper+new_fast_Z1); "
                "q1=logsum(duy0+ln2,duy1+ln2,dmy0,dmy1+ln2,dmy2)"
                "-ln(Llower)+saved denominator correction; full log_Rm subtracted once"
            ),
            "duy0_plus_ln2": None if old_uy[0] is None else _record(old_uy[0] + log_two),
            "duy1_plus_ln2": None if new_uy1 is None else _record(new_uy1 + log_two),
            "new_fast_Z1_plus_saved_log_Pstar_E_upper": (
                None if new_fast_Z1 is None else _record(new_fast_Z1 + log_pe)
            ),
            "saved_slow_Z1": _pack_log(1, old_slow[1]),
            "dmy0": None if old_my[0] is None else _record(old_my[0]),
            "dmy1_plus_ln2": None if old_my[1] is None else _record(old_my[1] + log_two),
            "dmy2": None if old_my[2] is None else _record(old_my[2]),
            "leading_before_L": None if leading is None else _record(leading),
            "L_positive_lower": _record(L_floor),
            "log_L_lower": _record(log_L_lower),
            "leading_after_L": None if leading_q1 is None else _record(leading_q1),
            "saved_denominator_correction": _record(correction),
            "new_q1": None if new_q1 is None else _record(new_q1),
            "full_log_Rm_subtracted_once": _record(log_rm),
        },
        "external_pressure_boundary": copy.deepcopy(_at(rm, "pressure_boundary")),
        "old_caps_retained": True,
        "caps_used_as_signed_values": False,
    }


def endpoint_C1(c: MPIntervalContext, parent: Mapping[str, Any],
                new_fast_Z1: Any) -> dict[str, Any]:
    """Integrate a new fast row into the saved Rm endpoint C1 operator."""

    _require(isinstance(c, MPIntervalContext), "c must be an MPIntervalContext")
    return _endpoint_packets(c, parent, new_fast_Z1)


def _compute_packets(c: MPIntervalContext, signed: Mapping[str, Any],
                     fast: Mapping[str, Any], rm: Mapping[str, Any],
                     native_v2: Mapping[str, Any]) -> dict[str, Any]:
    signed_actual = _signed_binding(signed)
    fast_actual = _fast_binding(fast)
    native_actual = _native_binding(native_v2)
    rm_native, rm_loop, rm_endpoint, source_binding, e_rows, old_fast_rows = \
        _source_rows(rm, c)
    p2_chart = _at(signed_actual, "p2", "payload", "charts", CHART)
    ratio_log = _signed_abs_log(
        c, p2_chart["p2_Z_over_p2"], "signed native p2_Z_over_p2"
    )
    e0, e_z = e_rows[0, 0], e_rows[0, 1]
    new_fast, path = _pure_bound_path(c, e0, e_z, ratio_log)
    direct = pure_bound(c, e0, e_z, ratio_log)
    _require(_same_interval(new_fast, direct),
             "Pure changed-bound path and public bound disagree")
    old_fast = copy.deepcopy(
        _at(fast_actual, "fast", "old_B_phi_over_Pstar_over_frozen_E0")
    )
    _require(old_fast == old_fast_rows,
             "Saved old fast cap differs between loop and fast reports")
    domain = _domain_receipt(c, native_actual, signed_actual, rm_native)
    endpoint = _endpoint_packets(c, rm, new_fast)
    histories = _at(signed_actual, "histories", "payload")
    preserved = signed.get("preserved", {})
    result = {
        "passed": True,
        "schema": "original_post_restore_logu_endpoint_v1",
        "stage": "current_original_post_restore_logu_endpoint",
        "chart": CHART,
        "input_reports": {
            "signed_native_sources": {
                "path": SIGNED_REPORT_PATH, "bytes": SIGNED_REPORT_BYTES,
                "sha256": SIGNED_REPORT_SHA256,
            },
            "fast_phi_correlated_C1": {
                "path": FAST_REPORT_PATH, "bytes": FAST_REPORT_BYTES,
                "sha256": FAST_REPORT_SHA256,
            },
            "directed_Rm_controls": {
                "path": RM_REPORT_PATH, "bytes": RM_REPORT_BYTES,
                "sha256": RM_REPORT_SHA256,
            },
            "native_directed_controls_v2": {
                "path": NATIVE_REPORT_PATH, "bytes": NATIVE_REPORT_BYTES,
                "sha256": NATIVE_REPORT_SHA256,
            },
        },
        "ordinary_Z_orders": list(ORDERS),
        "physical_y_orders": [0, 1],
        "source_binding": source_binding,
        "signed_native_binding": {
            "native_report_sha256": signed.get("native_report_sha256"),
            "signed_p2_value_sign": p2_chart.get("p2_value_sign"),
            "signed_p2_Z_over_p2_consumed": True,
            "caps_promoted_to_values": False,
        },
        "source_formula": (
            "bphi_Z=-C*(e_Z*f(that)+e*W*(u_Z/u)); |W|<=600; C<=2"
        ),
        "uniform_derivative_bound": {
            "W_upper": 600,
            "formula": "W=u*f'(t)*(t_u-tpsi*Iu/(1+t^2))",
            "f": "t/(1+t^2)",
            "fprime": "(1-t^2)/(1+t^2)^2",
            "global_abs_f_le_half": True,
            "global_abs_fprime_le_one": True,
            "not_differentiated_from_value_cap": True,
            "domain_abs_u_ge_one_required": True,
        },
        "C_bound": {
            "formula": "C=a*(1+2*q^2)/2",
            "a_source_identity": "4/5",
            "q_squared_source_identity": "q^2=3/4+5*eta/4",
            "q_squared_upper": "2",
            "C_upper": "2",
            "log_C_upper": _record(c.ln(c.mpf(2))),
            "source_bound_used": True,
        },
        "domain_receipt": domain,
        "old_B_phi_over_Pstar_over_frozen_E0": old_fast,
        "old_fast_caps_retained": copy.deepcopy(old_fast),
        "new_B_phi_over_Pstar_over_frozen_E0": (
            _replace_fast_z1(old_fast, new_fast)
        ),
        "old_value_row_0_0": copy.deepcopy(old_fast[0]),
        "new_fast_Z1": _pack_log(1, new_fast),
        "new_fast_Z1_arithmetic": path,
        "physical_deltaUy_Z1": endpoint["new_deltaUy_Z1"],
        "physical_deltaUy_Z1_arithmetic": endpoint["arithmetic"],
        "endpoint_C1": endpoint,
        "old_failed_endpoint_C1": copy.deepcopy(
            _at(fast_actual, "endpoint_C1")
        ),
        "history_error_receipts": copy.deepcopy(histories),
        "inherited_history_errors_preserved": preserved.get("all_inherited_errors") is True,
        "nonreset_M_memory_preserved": histories.get("incoming_A_B_C_retained") is True,
        "external_Rm_handles": copy.deepcopy(preserved.get("external_Rm_handles",
                                                          [183696, 183697])),
        "external_Rm_coefficient": preserved.get("external_coefficient", -1),
        "external_pressure_boundary": copy.deepcopy(endpoint["external_pressure_boundary"]),
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
            "stored_logs_exponentiated": False,
            "physical_Pstar_E0_applied_once": True,
            "full_log_Rm_subtracted_once": True,
        },
        "scope": {
            "new_signed_source_derivative_interface": True,
            "new_fast_phi_Z1_only": True,
            "value_row_reused_from_saved_loop": True,
            "old_bound_retained": True,
            "endpoint_C1_computed": True,
            "global_five_moment_closure": False,
            "full_goal_complete": False,
        },
    }
    json.dumps(result, allow_nan=False)
    return result


def _replace_fast_z1(raw_rows: Any, bound: Any) -> list[dict[str, Any]]:
    rows = [copy.deepcopy(dict(row)) for row in raw_rows]
    _require(len(rows) == 8, "Canonical old fast rows are required")
    for number, expected in enumerate(MIXED_INDICES):
        _require((rows[number].get("physical_y_order"),
                  rows[number].get("ordinary_Z_order")) == expected,
                 "Old fast rows are not canonical")
    position = MIXED_INDICES.index((0, 1))
    rows[position]["exact_zero"] = bound is None
    rows[position]["log_absolute_upper_bound"] = (
        None if bound is None else _record(bound)
    )
    return rows


def compute(c: MPIntervalContext, report_bytes: Any,
            fast_report_bytes: bytes | Mapping[str, Any] | None = None,
            rm_report_bytes: bytes | Mapping[str, Any] | None = None,
            native_report_bytes: bytes | Mapping[str, Any] | None = None
            ) -> dict[str, Any]:
    """Compute the new adapter from four authenticated saved reports.

    For tests and fresh fixtures, mappings with the same packet shape are
    accepted.  Production bytes must be supplied for all four reports; no
    filesystem fallback is performed.
    """

    _require(isinstance(c, MPIntervalContext), "c must be an MPIntervalContext")
    if isinstance(report_bytes, (str, Path)) and fast_report_bytes is None \
            and rm_report_bytes is None and native_report_bytes is None:
        root = Path(report_bytes)
        _require(root.is_dir(), "Saved-report root must be a directory")
        signed_raw = (root / SIGNED_REPORT_PATH).read_bytes()
        fast_raw = (root / FAST_REPORT_PATH).read_bytes()
        rm_raw = (root / RM_REPORT_PATH).read_bytes()
        native_raw = (root / NATIVE_REPORT_PATH).read_bytes()
        signed = _admit_report(signed_raw, SIGNED_REPORT_BYTES,
                               SIGNED_REPORT_SHA256, "signed-native report")
        fast = _admit_report(fast_raw, FAST_REPORT_BYTES,
                             FAST_REPORT_SHA256, "fast report")
        rm = _admit_report(rm_raw, RM_REPORT_BYTES,
                           RM_REPORT_SHA256, "Rm report")
        native = _admit_report(native_raw, NATIVE_REPORT_BYTES,
                               NATIVE_REPORT_SHA256, "native-v2 report")
        return _compute_packets(c, signed, fast, rm, native)
    if isinstance(report_bytes, Mapping) and fast_report_bytes is None \
            and rm_report_bytes is None and native_report_bytes is None:
        bundle = report_bytes
        signed = bundle.get("signed_native_sources", bundle.get("signed"))
        fast = bundle.get("fast_phi_correlated_C1", bundle.get("fast"))
        rm = bundle.get("directed_Rm_controls", bundle.get("rm"))
        native = bundle.get("native_directed_controls_v2", bundle.get("native"))
        _require(all(isinstance(item, Mapping)
                     for item in (signed, fast, rm, native)),
                 "Mapping compute input must contain four report mappings")
        return _compute_packets(c, signed, fast, rm, native)
    _require(all(item is not None for item in
                 (fast_report_bytes, rm_report_bytes, native_report_bytes)),
             "Four authenticated saved reports are required")
    signed = (_admit_report(report_bytes, SIGNED_REPORT_BYTES,
                            SIGNED_REPORT_SHA256, "signed-native report")
              if isinstance(report_bytes, bytes) else report_bytes)
    fast = (_admit_report(fast_report_bytes, FAST_REPORT_BYTES,
                          FAST_REPORT_SHA256, "fast report")
            if isinstance(fast_report_bytes, bytes) else fast_report_bytes)
    rm = (_admit_report(rm_report_bytes, RM_REPORT_BYTES,
                        RM_REPORT_SHA256, "Rm report")
          if isinstance(rm_report_bytes, bytes) else rm_report_bytes)
    native = (_admit_report(native_report_bytes, NATIVE_REPORT_BYTES,
                            NATIVE_REPORT_SHA256, "native-v2 report")
              if isinstance(native_report_bytes, bytes) else native_report_bytes)
    _require(all(isinstance(item, Mapping) for item in (signed, fast, rm, native)),
             "Decoded reports must be mappings")
    return _compute_packets(c, signed, fast, rm, native)


def _source_preflight() -> None:
    raw = Path(__file__).read_bytes()
    _require(b"\r" not in raw, "New Python source must use LF line endings")
    _require(raw.endswith(b"\n") and not raw.endswith(b"\n\n"),
             "New Python source must have exactly one EOF newline")
    ast.parse(raw, filename=MODULE)


def _rational_interval(c: MPIntervalContext, value: Fraction):
    value = Fraction(value)
    scalar = c.mpf(value.numerator) / value.denominator
    return c.mpf([scalar, scalar])


def _rational_log_cap(c: MPIntervalContext, value: Fraction):
    return _safe_abs_log(c, _rational_interval(c, value), "rational source")


def _contains(c: MPIntervalContext, cap: Any, value: Fraction) -> bool:
    expected = _rational_log_cap(c, value)
    if expected is None:
        return cap is None
    return cap is not None and _hi(cap) >= _hi(expected)


def _synthetic_ordinary_rows(c: MPIntervalContext,
                             values: list[Fraction]) -> list[dict[str, Any]]:
    return [
        {
            "ordinary_Z_order": order,
            "exact_zero": Fraction(value) == 0,
            "log_absolute_upper_bound": (
                None if Fraction(value) == 0
                else _record(_rational_log_cap(c, Fraction(value)))
            ),
        }
        for order, value in enumerate(values)
    ]


def _synthetic_mixed_rows(c: MPIntervalContext,
                          values: list[Fraction]) -> list[dict[str, Any]]:
    _require(len(values) == 8, "Synthetic mixed rows require eight values")
    return [
        {
            "physical_y_order": i,
            "ordinary_Z_order": k,
            "exact_zero": Fraction(value) == 0,
            "log_absolute_upper_bound": (
                None if Fraction(value) == 0
                else _record(_rational_log_cap(c, Fraction(value)))
            ),
        }
        for (i, k), value in zip(MIXED_INDICES, values)
    ]


def _signed_record(c: MPIntervalContext, value: Fraction) -> dict[str, Any]:
    return _record(_rational_interval(c, value))


def _synthetic_reports(c: MPIntervalContext) -> dict[str, Any]:
    F = Fraction
    e_rows = _synthetic_mixed_rows(c, [F(3, 5), F(2, 7), F(1, 5), F(2, 9),
                                      F(1, 4), F(3, 10), F(2, 11), F(1, 6)])
    old_fast = _synthetic_mixed_rows(c, [F(1, 9), F(1, 10), F(1, 11), F(1, 12),
                                        F(1, 13), F(1, 14), F(1, 15), F(1, 16)])
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
            "a": "4/5", "q_squared_post_y_and_Z_constant": True,
            "u": "p2*q/dstar",
        },
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
        "ordinary_Z_orders": list(ORDERS), "physical_y_orders": [0, 1],
        "B_phi_over_Pstar_over_frozen_E0": old_fast,
        "q_control_proof": {
            "q_squared_constant_in_y_Z": True,
            "q_value_not_selected": True,
            "q_squared_identity": "q^2=3/4+5*eta/4",
        },
    }
    uy = _synthetic_ordinary_rows(c, [F(1, 8), F(1, 9), F(1, 10), F(1, 11)])
    slow = _synthetic_ordinary_rows(c, [F(1, 12), F(1, 13), F(1, 14), F(1, 15)])
    my = _synthetic_ordinary_rows(c, [F(1, 16), F(1, 17), F(1, 18), F(1, 19)])
    qy = _synthetic_ordinary_rows(c, [F(1, 20), F(1, 21)])
    native_norm = [
        {"ordinary_Z_order": k, "exact_zero": False,
         "contains_zero": True, "sign": None,
         "log_absolute_upper_bound": _record(_rational_log_cap(c, F(1, 22 + k)))}
        for k in (0, 1)
    ]
    rm_endpoint = {
        "passed": True,
        "schema": "original_directed_Rm_endpoint_kernel_v1",
        "chart": CHART,
        "packet_metadata": {
            "mixed_native_source_native_report_sha256": NATIVE_REPORT_SHA256,
            "mixed_loop_source_native_parent_report_sha256": NATIVE_REPORT_SHA256,
            "chart": CHART,
        },
        "contracts": {"caps_used_as_signed_values": False,
                      "source_values_selected": False},
        "finite_deltaUy": {"slow_terms": slow,
                            "log_absolute_upper_bounds": uy},
        "finite_deltaMy": {"log_absolute_upper_bounds": my},
        "finite_deltaQy": {
            "log_absolute_upper_bounds": qy,
            "operator_receipt": {
                "positive_L_lower_used": True,
                "q1_denominator_derivative_correction_included": True,
                "delta_positive_and_less_than_one_checked": True,
                "L_positive_lower": _record(_rational_interval(c, F(9, 10))),
                "q1_denominator_correction": _record(
                    _rational_log_cap(c, F(1, 20))
                ),
            },
        },
        "saved_finite_error_bounds": {
            "log_Pstar_E_upper": _record(_rational_log_cap(c, F(7))),
        },
        "normalization": {
            "log_Rm": _record(_rational_log_cap(c, F(3))),
            "native_Qy_over_Rm_signed_log_scale": native_norm,
            "finite_delta_Qy_over_Rm_log_absolute_bounds": qy,
            "combined_Qy_over_Rm_log_absolute_bounds": qy,
        },
        "pressure_boundary": {
            "quantity": "-Rm*Qy", "coefficient": -1,
            "external_source_handles": [183696, 183697],
            "signed_coefficients_by_Z": {"0": -1, "1": -1},
        },
    }
    rm = {"actual_directed_Rm_controls": {
        "native_mixed": native, "loop_mixed": loop,
        "Rm_endpoint": rm_endpoint,
    }}
    signed_actual = {
        "gamma2": {"accepted_ancestors_replayed": False},
        "histories": {
            "accepted_ancestors_replayed": False,
            "payload": {
                "passed": True, "incoming_A_B_C_retained": True,
                "all_inherited_errors": True,
            },
        },
        "p2": {
            "accepted_ancestors_replayed": False,
            "payload": {
                "passed": True,
                "physical_R_applied_once": True,
                "midpoint_or_logcap_promoted": False,
                "full_native_p2_signed_scaled_enclosure": True,
                "charts": {CHART: {
                    "p2_value_sign": "negative",
                    "K": [{"ordinary_Z_order": 0,
                           "value": _signed_record(c, F(-1, 2)),
                           "kind": "signed_function_enclosure_not_selected_value"}],
                    "p2_Z_over_p2": _signed_record(c, F(5, 11)),
                }},
            },
        },
    }
    signed = {
        "schema": "original_signed_native_sources_v1",
        "native_report_sha256": NATIVE_REPORT_SHA256,
        "actual_signed_native_sources": signed_actual,
        "preserved": {"all_inherited_errors": True,
                      "external_Rm_handles": [183696, 183697],
                      "external_coefficient": -1},
    }
    fast = {"actual_fast_phi_correlated_C1": {
        "fast": {
            "passed": True, "schema": "original_post_restore_fast_phi_bound_v1",
            "input_report_sha256": RM_REPORT_SHA256,
            "old_B_phi_over_Pstar_over_frozen_E0": old_fast,
            "contracts": {"caps_used_as_signed_values": False},
        },
        "endpoint_C1": {"passed": True, "schema": "original_post_restore_fast_phi_endpoint_C1_v1"},
    }}
    native_v2 = {"actual_native_directed_controls_v2": {
        "schema": "original_native_directed_controls_v2_v1", "passed": True,
        "scope": {"corrected_native_controls_verified": True},
        "checks": {"accepted_ancestors_replayed": False, "source_executed": False},
        "raw_interval_tree": {
            "shared": {"z": _record(_rational_interval(c, F(1, 4))),
                        "logeta": _record(_rational_interval(c, F(-1))),
                        "logdstar": _record(_rational_log_cap(c, F(1, 3)))},
            "long_reference": {"charts": {"reference": {"scale": {
                "logP": _record(_rational_log_cap(c, F(10)))
            }}}},
            "restoration": {"charts": {CHART: {"scales": {
                "log_R": _record(_rational_log_cap(c, F(10))),
                "offsets": _record(c.mpf([-7, -6])),
            }}}},
        },
    }}
    return {"signed": signed, "fast": fast, "rm": rm, "native": native_v2}


def _fixture_bound_check() -> dict[str, Any]:
    """Run a fresh rational changed-bound and endpoint fixture."""

    _source_preflight()
    c = MPIntervalContext()
    c.prec = 320
    F = Fraction
    e0, ez, ratio = F(3, 5), F(2, 7), F(5, 11)
    bound, path = _pure_bound_path(
        c, _rational_log_cap(c, e0), _rational_log_cap(c, ez),
        _rational_log_cap(c, ratio)
    )
    public = pure_bound(c, _rational_log_cap(c, e0),
                        _rational_log_cap(c, ez),
                        _rational_log_cap(c, ratio))
    _require(_same_interval(bound, public), "Fresh changed-bound fixture mismatch")
    t = F(1, 3)
    f = t / (1 + t * t)
    fp = (1 - t * t) / (1 + t * t) ** 2
    direct = -F(2) * (ez * f + e0 * fp * ratio)
    _require(_contains(c, bound, direct),
             "Fresh changed-bound cap does not contain rational derivative")
    # A direct rational witness for the uniform W term stays far inside 600.
    u, tpsi, Iu, t_u = F(2), F(1, 2), F(3), F(3, 2)
    # Synthetic independent t_u witness; it is not a saved source value.
    W = u * fp * (t_u - tpsi * Iu / (1 + t * t))
    _require(abs(W) <= 600, "Fresh W<=600 rational fixture failed")
    mixed_signed_cap_rejected = False
    try:
        _read_signed_interval(
            c,
            {
                "value": _signed_record(c, F(-1, 2)),
                "kind": "signed_function_enclosure_not_selected_value",
                "log_absolute_upper_bound": _record(_rational_log_cap(c, F(1, 2))),
            },
            "fresh mixed signed source",
        )
    except ValueError:
        mixed_signed_cap_rejected = True
    _require(mixed_signed_cap_rejected,
             "Fresh mixed signed cap was accepted as a source value")
    synthetic = _synthetic_reports(c)
    result = _compute_packets(c, synthetic["signed"], synthetic["fast"],
                              synthetic["rm"], synthetic["native"])
    _require(result["domain_receipt"]["certified_abs_u_ge_one"] is True,
             "Fresh log|u| domain receipt failed")
    _require(result["contracts"]["physical_Pstar_E0_applied_once"] is True,
             "Fresh endpoint Pstar factor was not applied once")
    _require(result["contracts"]["full_log_Rm_subtracted_once"] is True,
             "Fresh endpoint logRm factor was not subtracted once")
    serialized = json.dumps(result, allow_nan=False)
    checks = {
        "source_preflight_before_fixture_math": True,
        "fresh_context_prec_320": c.prec == 320,
        "changed_bound_contains_direct_rational_Z1": True,
        "uniform_W_rational_witness": True,
        "mixed_signed_cap_rejected": mixed_signed_cap_rejected,
        "domain_abs_u_ge_one_certified": True,
        "physical_Pstar_E0_once": True,
        "full_log_Rm_once": True,
        "old_fast_caps_retained": result["old_fast_caps_retained"] ==
        synthetic["fast"]["actual_fast_phi_correlated_C1"]["fast"]
        ["old_B_phi_over_Pstar_over_frozen_E0"],
        "native_Qy_Z1_retained": result["endpoint_C1"]["native_Qy_over_Rm_Z1"] is not None,
        "external_handles_retained": result["external_Rm_handles"] == [183696, 183697],
        "no_exp_materialized": result["contracts"]["stored_logs_exponentiated"] is False,
        "serializes": isinstance(serialized, str),
    }
    _require(all(checks.values()), "Fresh logu endpoint fixture failed")
    receipt = {
        "passed": True,
        "schema": "original_post_restore_logu_endpoint_fixture_v1",
        "module_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "context_precision": c.prec,
        "checks": checks,
        "direct_rational_Z1": str(direct),
        "rational_W": str(W),
        "source_formula": "W=u*f'(t)*(t_u-tpsi*Iu/(1+t^2)); |W|<=600",
        "synthetic_t_u": str(t_u),
        "new_bound_path": path,
        "endpoint_schema": result["endpoint_C1"]["schema"],
        "old_endpoint_retained": True,
    }
    scratch = Path("tmp/logu_endpoint_fixture.json")
    scratch.parent.mkdir(parents=True, exist_ok=True)
    scratch.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n",
                       encoding="utf-8", newline="\n")
    return receipt


def fixture_bound_check() -> dict[str, Any]:
    return _fixture_bound_check()


def self_check() -> dict[str, Any]:
    return fixture_bound_check()


__all__ = (
    "SIGNED_REPORT_PATH", "SIGNED_REPORT_BYTES", "SIGNED_REPORT_SHA256",
    "FAST_REPORT_PATH", "FAST_REPORT_BYTES", "FAST_REPORT_SHA256",
    "RM_REPORT_PATH", "RM_REPORT_BYTES", "RM_REPORT_SHA256",
    "NATIVE_REPORT_PATH", "NATIVE_REPORT_BYTES", "NATIVE_REPORT_SHA256",
    "MODULE", "CHART", "ORDERS", "MIXED_INDICES", "pure_bound",
    "compute", "endpoint_C1", "fixture_bound_check", "self_check",
)
