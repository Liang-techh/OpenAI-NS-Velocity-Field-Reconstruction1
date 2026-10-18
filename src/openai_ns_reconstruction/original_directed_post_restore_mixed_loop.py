"""Pure post-restore mixed physical-y / ordinary-Z loop controls.

This module is a narrow consumer of the new post-restore mixed-native packet.
It intentionally accepts the previous directed loop report only as immutable
input evidence.  The old report is never rebuilt, imported, or used as the
source of the new mixed arithmetic.

The arithmetic is logarithmic and absolute.  ``LogMixed`` carries the eight
ordinary rows ``(physical_y_order, ordinary_Z_order)`` for ``y`` orders zero
and one and ``Z`` orders zero through three.  ``LogPhase`` carries the free
angle total-order-four rows.  The public loop rows are ordinary derivatives;
``inverse_coefficients`` returns normalized inverse increments internally and
restores ordinary factorials before the rows are emitted.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path

import mpmath as mp
from mpmath import isfinite, nstr
from mpmath.ctx_iv import MPIntervalContext

from .original_log_mixed_jets import LogMixed, _coerce
from .original_log_phase_jets import LogPhase, compose_inverse, inverse_coefficients


OLD_LOOP_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_directed_loop_controls.json"
)
OLD_LOOP_REPORT_BASENAME = OLD_LOOP_REPORT.rsplit("/", 1)[-1]
OLD_LOOP_REPORT_BYTES = 1752702
OLD_LOOP_REPORT_SHA256 = (
    "adefdae708754a2572eae15a90662274bf17e88a0d72d5db9befbc9e2858b252"
)
NATIVE_PARENT_REPORT = (
    "lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json"
)
NATIVE_PARENT_REPORT_BYTES = 42177499
NATIVE_PARENT_REPORT_SHA256 = (
    "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"
)
MODULE = "src/openai_ns_reconstruction/original_directed_post_restore_mixed_loop.py"
CHART = "post_restore"
ORDERS = (0, 1, 2, 3)
MIXED_INDICES = tuple((i, k) for i in (0, 1) for k in ORDERS)
ROOT_NAMES = (
    "a",
    "b",
    "t0",
    "p2",
    "q",
    "u",
    "h",
    "r",
    "alpha",
    "K",
    "E_over_frozen_E0",
)
PUBLIC_ROWS = (
    "A",
    "B_over_Pstar_over_frozen_E0",
    "B_phi_over_Pstar_over_frozen_E0",
    "inverse_Psi",
)
SHARED_Y0_ALGEBRA_RECALCULATED_FOR_NEW_GRID = True
SHARED_Y0_RECALCULATION_SCOPE = (
    "new mixed y0/Z and total-order-four free-angle phase algebra only; "
    "saved Z-only/public y0 rows remain input receipts"
)

def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _raw_mpf(value):
    """Decode an exact MPF tuple without changing either endpoint."""
    if not isinstance(value, (list, tuple)) or len(value) != 4:
        raise ValueError("An exact four-part MPF tuple is required")
    if not all(type(part) is int for part in value):
        raise ValueError("An exact four-part MPF tuple is required")
    return mp.mpf(tuple(value))


def _record(value):
    """Serialize a directed interval while retaining both exact endpoints."""
    _require(hasattr(value, "_mpi_"), "A directed interval is required")
    lower = tuple(value._mpi_[0])
    upper = tuple(value._mpi_[1])
    return {
        "lower_exact_mpf_tuple": list(lower),
        "upper_exact_mpf_tuple": list(upper),
        "lower_display": nstr(_raw_mpf(lower), 24),
        "upper_display": nstr(_raw_mpf(upper), 24),
        "directed": True,
    }


def _read_record(c: MPIntervalContext, value, label: str = "directed record"):
    """Read a directed record without selecting or negating an endpoint."""
    if hasattr(value, "_mpi_"):
        out = _coerce(c, value)
    else:
        if isinstance(value, Mapping) and "value" in value \
                and "lower_exact_mpf_tuple" not in value:
            value = value["value"]
        _require(isinstance(value, Mapping), label + " must be a mapping")
        _require(value.get("directed") is True,
                 label + " must set directed=True")
        try:
            lower = _raw_mpf(value["lower_exact_mpf_tuple"])
            upper = _raw_mpf(value["upper_exact_mpf_tuple"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("Malformed " + label) from exc
        _require(lower <= upper, "Reversed " + label)
        # Construct an interval from both source endpoints.  In particular,
        # do not turn one endpoint into a scalar and negate it.
        out = c.mpf([lower, upper])
    try:
        lower = out.a
        upper = out.b
    except (AttributeError, TypeError, ValueError) as exc:
        raise ValueError("Malformed " + label) from exc
    _require(isfinite(lower) and isfinite(upper), label + " must be finite")
    _require(lower <= upper, "Reversed " + label)
    return out


def _jsonable(value):
    """Convert an in-memory packet to JSON-safe values for receipts."""
    if hasattr(value, "_mpi_"):
        return _record(value)
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError("Packet contains a non-JSON value: " + type(value).__name__)


def _packed_row_index(row):
    if not isinstance(row, Mapping):
        raise ValueError("A mixed row must be a mapping")
    i = row.get("physical_y_order", row.get("y_order"))
    k = row.get("ordinary_Z_order", row.get("Z_order"))
    _require(type(i) is int and type(k) is int,
             "A mixed row has no integer physical-y/Z index")
    _require((i, k) in MIXED_INDICES, "A mixed row index is outside y0/y1,Z0..Z3")
    return i, k


def _row_value(c, row, label):
    exact = row.get("exact_zero")
    _require(type(exact) is bool, label + " exact_zero must be boolean")
    raw = row.get("log_absolute_upper_bound")
    if raw is None and "record" in row:
        raw = row["record"]
    if raw is None and "log_cap" in row and not exact:
        raw = row["log_cap"]
    if exact:
        _require(raw is None, label + " exact-zero row carries a bound")
        return None
    _require(raw is not None, label + " nonzero row has no directed record")
    return _read_record(c, raw, label + ".log_absolute_upper_bound")


def _mixed_from_rows(c, raw, label):
    """Decode exactly the canonical eight mixed rows."""
    if isinstance(raw, LogMixed):
        _require(raw.ctx is c, label + " uses a different interval context")
        return raw
    if isinstance(raw, Mapping):
        for key in ("rows", "mixed_log_caps", "log_caps", "values"):
            if key in raw:
                return _mixed_from_rows(c, raw[key], label)
        if set(raw) == set(MIXED_INDICES):
            values = {
                index: (None if raw[index] is None
                        else _read_record(c, raw[index], label + str(index)))
                for index in MIXED_INDICES
            }
            return LogMixed.from_mapping(c, values)
        # A convenient in-memory form used by focused callers.
        if set(raw) >= {"y0", "y1"}:
            y0, y1 = raw["y0"], raw["y1"]
            _require(isinstance(y0, (list, tuple)) and isinstance(y1, (list, tuple)),
                     label + " y0/y1 rows must be lists")
            raw = list(y0) + list(y1)
    _require(isinstance(raw, (list, tuple)) and len(raw) == 8,
             label + " must contain exactly eight mixed rows")
    values = {}
    for row_number, row in enumerate(raw):
        index = _packed_row_index(row)
        _require(index not in values, label + " repeats row " + str(index))
        values[index] = _row_value(c, row, label + "[" + str(row_number) + "]")
    _require(set(values) == set(MIXED_INDICES),
             label + " must cover physical y orders zero and one and Z orders zero through three")
    return LogMixed.from_mapping(c, values)


def _pack_mixed(jet):
    return [
        {
            "physical_y_order": i,
            "ordinary_Z_order": k,
            "exact_zero": jet[i, k] is None,
            "log_absolute_upper_bound": (
                None if jet[i, k] is None else _record(jet[i, k])
            ),
        }
        for i, k in MIXED_INDICES
    ]


def _pack_phase(jet):
    return [
        {
            "physical_y_order": i,
            "ordinary_Z_order": k,
            "free_angle_order": j,
            "exact_zero": jet[i, k, j] is None,
            "log_absolute_upper_bound": (
                None if jet[i, k, j] is None else _record(jet[i, k, j])
            ),
        }
        for i, k, j in LogPhase.indices
    ]


def _chart(packet):
    """Locate the post-restore chart without changing its source mapping."""
    _require(isinstance(packet, Mapping), "The mixed-native packet is required")
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


def _root_mapping(packet, chart):
    for owner in (chart, packet):
        if not isinstance(owner, Mapping):
            continue
        for key in ("root_mixed_caps", "root_controls", "root_caps", "log_caps"):
            candidate = owner.get(key)
            if isinstance(candidate, Mapping):
                return candidate
    raise ValueError("The mixed-native packet has no root_mixed_caps mapping")


def _root_caps(c, packet, chart):
    roots = _root_mapping(packet, chart)
    aliases = {"E_over_frozen_E0": ("E_over_frozen_E0", "E_relative", "E")}
    result = {}
    for name in ROOT_NAMES:
        keys = aliases.get(name, (name,))
        raw = next((roots[key] for key in keys if key in roots), None)
        _require(raw is not None or any(key in roots for key in keys),
                 "Missing mixed-native root control: " + name)
        result[name] = _mixed_from_rows(c, raw, "root_mixed_caps." + name)
    return result


def _conditioning(c, packet, chart):
    for owner in (chart, packet):
        if not isinstance(owner, Mapping):
            continue
        for key in ("conditioning", "conditioning_caps", "root_conditioning"):
            candidate = owner.get(key)
            if isinstance(candidate, Mapping):
                if "log_D_lower" in candidate and "log_lambda_lower" in candidate:
                    return {
                        "log_D_lower": _read_record(
                            c, candidate["log_D_lower"], "conditioning.log_D_lower"
                        ),
                        "log_lambda_lower": _read_record(
                            c, candidate["log_lambda_lower"],
                            "conditioning.log_lambda_lower",
                        ),
                    }, candidate
        floors = owner.get("positive_floors")
        if isinstance(floors, Mapping):
            source = floors.get("source_conditioning")
            if isinstance(source, Mapping) \
                    and "log_D_lower" in source \
                    and "log_lambda_lower" in source:
                return {
                    "log_D_lower": _read_record(
                        c, source["log_D_lower"],
                        "positive_floors.source_conditioning.log_D_lower",
                    ),
                    "log_lambda_lower": _read_record(
                        c, source["log_lambda_lower"],
                        "positive_floors.source_conditioning.log_lambda_lower",
                    ),
                }, source
    raise ValueError("The mixed-native packet has no D/lambda conditioning floors")


def _source_native_sha(packet, chart):
    candidates = []
    for owner in (chart, packet):
        if isinstance(owner, Mapping):
            for key in (
                "source_native_report_sha256",
                "source_native_parent_report_sha256",
                "native_parent_report_sha256",
                "source_native_parent_sha256",
            ):
                if key in owner:
                    candidates.append(owner[key])
            provenance = owner.get("input_record_provenance")
            if isinstance(provenance, Mapping):
                for key in (
                    "source_native_report_sha256",
                    "source_native_parent_report_sha256",
                    "native_parent_report_sha256",
                    "source_native_parent_sha256",
                ):
                    if key in provenance:
                        candidates.append(provenance[key])
    _require(candidates, "The mixed-native packet must carry the native parent SHA256")
    _require(all(value == NATIVE_PARENT_REPORT_SHA256 for value in candidates),
             "The mixed-native packet native parent SHA256 changed")
    return NATIVE_PARENT_REPORT_SHA256


def _admit_old_loop_report(loop_report_bytes: bytes):
    _require(isinstance(loop_report_bytes, bytes),
             "loop_report_bytes must be immutable report bytes")
    _require(len(loop_report_bytes) == OLD_LOOP_REPORT_BYTES,
             "The saved directed loop report length changed")
    _require(_sha_bytes(loop_report_bytes) == OLD_LOOP_REPORT_SHA256,
             "The saved directed loop report hash changed")
    try:
        report = json.loads(loop_report_bytes)
    except (TypeError, ValueError) as exc:
        raise ValueError("The saved directed loop report is not JSON") from exc
    _require(isinstance(report, Mapping), "The saved directed loop report is not a mapping")
    component = report.get("actual_directed_loop_controls", report)
    _require(isinstance(component, Mapping), "The saved directed loop component is missing")
    _require(component.get("passed") is True,
             "The saved directed loop component did not pass")
    _require(component.get("native_parent_report_sha256") == NATIVE_PARENT_REPORT_SHA256,
             "The saved directed loop native parent binding changed")
    charts = component.get("charts")
    _require(isinstance(charts, Mapping) and isinstance(charts.get(CHART), Mapping),
             "The saved directed loop post_restore chart is missing")
    return report, component, charts[CHART]


def _phase_cosine(c):
    zero = c.mpf(0)
    return LogPhase.from_mapping(c, {
        index: zero if index[:2] == (0, 0) else None
        for index in LogPhase.indices
    })


def _angle_primitive(c):
    values = {index: None for index in LogPhase.indices}
    values[0, 0, 0] = c.ln(2 * c.pi)
    values[0, 0, 1] = c.mpf(0)
    return LogPhase.from_mapping(c, values)


def _primitive(engine, jet):
    """Integrate a free-angle jet over an interval of length ``2*pi``."""
    c = engine.c
    return LogPhase.from_mapping(c, {
        (i, k, j): (
            None if jet[i, k, 0] is None
            else jet[i, k, 0] + c.ln(2 * c.pi)
        ) if j == 0 else jet[i, k, j - 1]
        for i, k, j in LogPhase.indices
    })


def _root_identity(caps):
    a, b, t0, q, K = (caps[name] for name in ("a", "b", "t0", "q", "K"))
    _require(a[0, 0] is not None, "The mixed-native a value row is missing")
    _require(all(a[index] is None for index in MIXED_INDICES if index != (0, 0)),
             "The post-restore a rows are not constant")
    _require(all(b[index] is None for index in MIXED_INDICES),
             "The post-restore b identity is not exact zero")
    _require(all(t0[index] is None for index in MIXED_INDICES),
             "The post-restore t0 identity is not exact zero")
    _require(all(q[index] is None for index in MIXED_INDICES if index != (0, 0)),
             "The post-restore q control is not y/Z constant")
    _require(all(K[index] is None for index in MIXED_INDICES if index != (0, 0)),
             "The post-restore K control is not y/Z constant")
    return {
        "passed": True,
        "a_derivatives_exact_zero": True,
        "b_exact_zero": True,
        "t0_exact_zero": True,
        "q_squared_constant_in_y_Z": True,
        "K_constant_in_y_Z": True,
        "a_endpoint_identity": "4/5",
        "physical_y_is_slow_fixed_target_phase": True,
    }


def _compute(c, caps, conditioning):
    """Perform only new LogMixed/LogPhase arithmetic on the mixed grid."""
    m = LogMixed(c)
    p = LogPhase(c)
    a, b, t0, p2, q, u, h, r, alpha, K, E = (
        caps[name] for name in ROOT_NAMES
    )
    del p2, u, h  # retained source controls; the C5 graph uses r/alpha/K.

    log_D_lower = conditioning["log_D_lower"]
    log_lambda_lower = conditioning["log_lambda_lower"]
    one = p.constant(c.mpf(0))
    cosine = _phase_cosine(c)
    r_phase = p.lift_mixed(r)
    alpha_phase = p.lift_mixed(alpha)
    K_phase = p.lift_mixed(K)
    t0_phase = p.lift_mixed(t0)

    # D = 1 - 2*r*cos(psi) + r^2.  LogMixed/LogPhase are absolute
    # majorants, so subtraction is enclosed by addition of the two terms.
    D = p.add(
        one,
        p.scale(p.mul(r_phase, cosine), c.ln(2)),
        p.mul(r_phase, r_phase),
    )
    inv_D = p.reciprocal(D, log_D_lower)
    numerator = p.add(cosine, r_phase)
    t = p.add(t0_phase, p.mul(alpha_phase, p.mul(numerator, inv_D)))
    t2 = p.mul(t, t)
    lam = p.mul(K_phase, p.add(one, t2))
    T1 = _primitive(p, t)
    T2 = _primitive(p, t2)
    angle = _angle_primitive(c)
    raw_phi = p.mul(K_phase, p.add(angle, T2))
    phi_values = {
        index: raw_phi[index] if index[2] == 0
        else lam[index[0], index[1], index[2] - 1]
        for index in LogPhase.indices
    }
    # Phi(0)=0 is a normalization of the phase primitive, not a selected
    # source value.  Every nonconstant row remains a logarithmic cap.
    phi_values[0, 0, 0] = c.mpf(0)
    Phi = LogPhase.from_mapping(c, phi_values)
    inverse, increment, inverse_proof = inverse_coefficients(
        p, Phi, log_lambda_lower
    )
    T1_hat = compose_inverse(p, T1, increment)
    t_hat = compose_inverse(p, t, increment)
    lam_hat = compose_inverse(p, lam, increment)

    # The target phase is independent of physical y and Z.  It is represented
    # by the magnitude bound one (log one equals zero), not by a phase value.
    target_phase = m.constant(c.mpf(0))
    two_pi_log = c.ln(2 * c.pi)
    chi = m.add(target_phase, m.scale(inverse, -two_pi_log))
    A = m.scale(m.mul(a, chi), -c.ln(2))

    B_term_a = m.scale(m.mul(a, T1_hat), -two_pi_log)
    B_term_b = m.mul(b, target_phase)
    B_rel = m.scale(m.add(B_term_a, B_term_b), -c.ln(2))
    B = m.mul(E, B_rel)

    inv_lam_hat = m.reciprocal(lam_hat, log_lambda_lower)
    B_phi_term_a = m.scale(
        m.mul(a, m.mul(t_hat, inv_lam_hat)), -two_pi_log
    )
    B_phi_rel = m.scale(m.add(B_phi_term_a, b), -c.ln(2))
    B_phi = m.mul(E, B_phi_rel)

    # At post-restore b=t0=0 the endpoint value has an exact source identity.
    # Apply its sharp inequality to the value row only.  Derivative rows above
    # remain the complete mixed product-rule result.
    sharp_value = m.scale(m.mul(E, a), c.ln(c.mpf(5) / 4))
    b_phi_values = B_phi.as_dict()
    b_phi_values[0, 0] = sharp_value[0, 0]
    B_phi = LogMixed.from_mapping(c, b_phi_values)

    return {
        "a": a,
        "b": b,
        "t0": t0,
        "p2": caps["p2"],
        "q": q,
        "q_squared": m.mul(q, q),
        "u": caps["u"],
        "h": caps["h"],
        "r": r,
        "alpha": alpha,
        "K": K,
        "E_over_frozen_E0": E,
        "D": D,
        "inverse_D": inv_D,
        "t": t,
        "lambda": lam,
        "T1": T1,
        "T2": T2,
        "Phi": Phi,
        "inverse": inverse,
        "increment": increment,
        "T1hat": T1_hat,
        "that": t_hat,
        "lamhat": lam_hat,
        "target_phase": target_phase,
        "A": A,
        "B": B,
        "B_phi": B_phi,
        "inverse_D_log_lower": log_D_lower,
        "log_lambda_lower": log_lambda_lower,
        "inverse_proof": inverse_proof,
        "B_phi_value_sharpened": True,
    }


def _saved_y0_receipts(packet, chart, old_chart):
    """Keep saved Z-only inputs byte-for-byte as JSON receipts."""
    def first(owner, keys, default=None):
        if not isinstance(owner, Mapping):
            return default
        for key in keys:
            if key in owner:
                return owner[key]
        return default

    provenance = first(
        chart,
        ("input_record_provenance", "root_input_record_provenance"),
        first(packet, ("input_record_provenance", "root_input_record_provenance")),
    )
    if not isinstance(provenance, Mapping):
        provenance = None
    provenance_root = (
        provenance.get("root_caps") if provenance is not None else None
    )
    provenance_conditioning = (
        provenance.get("conditioning") if provenance is not None else None
    )
    provenance_native_caps = (
        provenance.get("native_caps") if provenance is not None else None
    )
    new_root = first(
        chart,
        ("saved_y0_root_caps", "saved_new_native_y0_root_caps", "y0_root_caps"),
        first(packet, ("saved_y0_root_caps", "saved_new_native_y0_root_caps", "y0_root_caps")),
    )
    new_conditioning = first(
        chart,
        ("saved_y0_conditioning", "saved_new_native_conditioning", "y0_conditioning"),
        first(packet, ("saved_y0_conditioning", "saved_new_native_conditioning", "y0_conditioning")),
    )
    # A new-native receipt is trusted only when its input provenance carries
    # all three saved source records.  The old loop chart is retained below as
    # an old receipt and is never promoted to a new-native receipt.
    native_receipt_complete = (
        provenance is not None
        and provenance_root is not None
        and provenance_conditioning is not None
        and provenance_native_caps is not None
    )
    if native_receipt_complete:
        new_root = provenance_root
        new_conditioning = provenance_conditioning
    else:
        new_root = None
        new_conditioning = None
    public = {
        key: old_chart.get(key)
        for key in (
            "A",
            "B_over_E",
            "B_over_Pstar_over_frozen_E0",
            "A_phi",
            "B_over_Pstar_over_frozen_E0_phi",
            "inverse_Psi",
        )
        if key in old_chart
    }
    return {
        "saved_y0_root_caps": _jsonable(new_root),
        "saved_y0_conditioning": _jsonable(new_conditioning),
        "saved_new_native_root_caps": _jsonable(provenance_root),
        "saved_new_native_conditioning": _jsonable(provenance_conditioning),
        "saved_new_native_caps": _jsonable(provenance_native_caps),
        "saved_new_native_input_record_provenance": _jsonable(provenance),
        "saved_old_loop_input_root_caps": _jsonable(old_chart.get("input_root_caps")),
        "saved_old_loop_public_y0_caps": _jsonable(public),
        "saved_old_loop_conditioning": _jsonable(old_chart.get("conditioning")),
        "new_native_y0_receipt_present": native_receipt_complete,
        "new_native_conditioning_receipt_present": native_receipt_complete,
        "new_native_receipt_complete": native_receipt_complete,
        "new_native_receipt_source": (
            "native.input_record_provenance"
            if native_receipt_complete else "old_loop_receipt_only"
        ),
    }


def compute(c: MPIntervalContext, mixed_native_packet, loop_report_bytes: bytes) -> dict:
    """Compute the post-restore mixed loop packet from two immutable inputs.

    ``mixed_native_packet`` is the new native colleague's in-memory packet.
    ``loop_report_bytes`` must be the exact saved directed-loop report bytes;
    its only role here is immutable parent binding and receipt preservation.
    """
    _require(isinstance(c, MPIntervalContext), "c must be an MPIntervalContext")
    _require(isinstance(mixed_native_packet, Mapping),
             "mixed_native_packet must be a mapping")
    _require(mixed_native_packet.get("passed") is True,
             "The mixed-native packet did not pass")
    report, component, old_chart = _admit_old_loop_report(loop_report_bytes)
    chart = _chart(mixed_native_packet)
    native_sha = _source_native_sha(mixed_native_packet, chart)
    _require(
        mixed_native_packet.get("source_executed", False) is False
        and mixed_native_packet.get("accepted_ancestors_replayed", False) is False,
        "The mixed-native packet source execution flags changed",
    )
    conditioning, conditioning_source = _conditioning(c, mixed_native_packet, chart)
    caps = _root_caps(c, mixed_native_packet, chart)
    identity = _root_identity(caps)
    computed = _compute(c, caps, conditioning)

    rows = {
        "A": _pack_mixed(computed["A"]),
        "B_over_Pstar_over_frozen_E0": _pack_mixed(computed["B"]),
        "B_phi_over_Pstar_over_frozen_E0": _pack_mixed(computed["B_phi"]),
        "inverse_Psi": _pack_mixed(computed["inverse"]),
    }
    root_rows = {
        name: _pack_mixed(computed[name])
        for name in ROOT_NAMES + ("q_squared",)
    }
    phase_names = ("D", "inverse_D", "t", "lambda", "T1", "T2", "Phi")
    phase_intermediates = {
        name: _pack_phase(computed[name]) for name in phase_names
    }
    # Keep the prior descriptive alias while the canonical formula name is
    # ``lambda``.  Both aliases point to the same newly computed phase rows.
    phase_intermediates["lambda_value"] = phase_intermediates["lambda"]
    inverse_intermediates = {
        name: _pack_mixed(computed[name])
        for name in ("inverse", "T1hat", "that", "lamhat")
    }
    receipt = _saved_y0_receipts(mixed_native_packet, chart, old_chart)
    result = {
        "passed": True,
        "schema": "original_directed_post_restore_mixed_loop_v1",
        "stage": "current_original_directed_post_restore_mixed_loop",
        "chart": CHART,
        "parent_report": OLD_LOOP_REPORT_BASENAME,
        "parent_report_sha256": OLD_LOOP_REPORT_SHA256,
        "parent_report_bytes": OLD_LOOP_REPORT_BYTES,
        "old_loop_report": OLD_LOOP_REPORT_BASENAME,
        "old_loop_report_sha256": OLD_LOOP_REPORT_SHA256,
        "old_loop_report_bytes": OLD_LOOP_REPORT_BYTES,
        "source_native_parent_report": NATIVE_PARENT_REPORT,
        "source_native_parent_report_sha256": native_sha,
        "source_native_parent_report_bytes": NATIVE_PARENT_REPORT_BYTES,
        "ordinary_Z_orders": list(ORDERS),
        "physical_y_orders": [0, 1],
        "physical_y_derivatives_claimed": True,
        "fixed_target_phase": True,
        "target_phase_Z_independent": True,
        "phase_Z_independent": True,
        "fast_phi_separate": True,
        "B_phi_has_no_N_factor": True,
        "root_mixed_caps": root_rows,
        "root_control_caps": root_rows,
        "mixed_log_caps": rows,
        **rows,
        "A_log_caps": rows["A"],
        "B_over_Pstar_over_frozen_E0_log_caps": rows[
            "B_over_Pstar_over_frozen_E0"
        ],
        "B_phi_over_Pstar_over_frozen_E0_log_caps": rows[
            "B_phi_over_Pstar_over_frozen_E0"
        ],
        "inverse_Psi_log_caps": rows["inverse_Psi"],
        "E_relative": root_rows["E_over_frozen_E0"],
        "phase_intermediates": phase_intermediates,
        "intermediate_phase_caps": phase_intermediates,
        "phase_intermediate_row_count": len(LogPhase.indices),
        "phase_indices": [list(index) for index in LogPhase.indices],
        "inverse_mixed_caps": inverse_intermediates,
        "intermediate_Z_caps": inverse_intermediates,
        "conditioning": _jsonable(conditioning_source),
        "conditioning_consumed_from_new_native": True,
        "saved_native_conditioning": receipt["saved_y0_conditioning"],
        "root_control_identity": identity,
        "q_control_proof": {
            "passed": True,
            "q_squared_identity": "q^2=3/4+5*eta/4",
            "q_squared_constant_in_y_Z": True,
            "eta_actual_saved_log_required_by_native": True,
            "eta_handle_not_used": True,
            "dstar_actual_saved_log_required_by_native": True,
            "dstar_handle_not_used": True,
            "q_value_not_selected": True,
        },
        "B_phi_value_control": {
            "passed": True,
            "exact_formula": (
                "B_phi/E0=-(E/E0)*a*(1+2*q^2)*t/[2*(1+t^2)]"
            ),
            "sharp_value_inequality": "|t/(1+t^2)|<=1/2",
            "q_squared_upper_used": "2",
            "value_cap_factor": "5/4 times |a| times E/E0",
            "value_row_replaced_by_sharp_cap": True,
            "derivative_rows_from_full_mixed_product": True,
            "derivative_bound_not_inferred_from_value_bound": True,
            "fast_term_has_no_1/N": True,
        },
        "source_formulas": {
            "D": "1-2*r*cos(psi)+r^2",
            "t": "t0+alpha*(cos(psi)-r)/D",
            "lambda": "K*(1+t^2)",
            "T1": "free-angle primitive of t",
            "T2": "free-angle primitive of t^2",
            "Phi": "K*(psi+T2), normalized by Phi(0)=0",
            "inverse_Psi": "inverse_coefficients(Phi, log_lambda_lower)",
            "A": "a/2*(target_phase-Psi/(2*pi))",
            "B_over_Pstar_over_frozen_E0": (
                "E/E0/2*(-a*T1hat/(2*pi)-b*target_phase)"
            ),
            "B_phi_over_Pstar_over_frozen_E0": (
                "E/E0/2*(-a*that/(2*pi*lambda_hat)-b)"
            ),
        },
        "saved_receipts": receipt,
        **receipt,
        "inverse_proof": computed["inverse_proof"],
        "source_executed": False,
        "accepted_ancestors_replayed": False,
        "shared_y0_algebra_recalculated_for_new_grid": (
            SHARED_Y0_ALGEBRA_RECALCULATED_FOR_NEW_GRID
        ),
        "shared_y0_recalculation_scope": SHARED_Y0_RECALCULATION_SCOPE,
        "old_loop_builder_called": False,
        "old_loop_report_rebuilt": False,
        "legacy_abs_log_used": False,
        "stored_logs_exponentiated": False,
        "scalar_endpoint_negation_used": False,
        "recurrences_replayed": False,
        "quadratures_replayed": False,
        "norms_replayed": False,
        "finite_N_phase_not_materialized": True,
        "full_goal_complete": False,
        "limits": {
            "post_restore_chart_only": True,
            "fixed_target_phase_Z0_through_Z3": True,
            "physical_y_is_slow_fixed_phase": True,
            "fast_phi_kept_separate": True,
            "derivative_rows_not_sharpened_by_value_cap": True,
            "legacy_loop_recertification": False,
            "full_N_numeric": False,
            "global_five_moment_closure": False,
        },
        "tool_metadata": {
            "route": "fixed luna_max_worker",
            "initial_fork_turns": "none",
            "requested_model": "gpt-5.6-luna/max",
            "independent_runtime_model_metadata_exposed": False,
        },
    }
    json.dumps(_jsonable(result), allow_nan=False)
    return result


def _rational_log(c, value):
    value = Fraction(value)
    if value == 0:
        return None
    return c.ln(c.mpf(abs(value.numerator)) / value.denominator)


def _phase_from_expression(c, expression, symbols):
    import sympy as sp

    y, z, psi = symbols
    return LogPhase.from_mapping(c, {
        index: _rational_log(
            c,
            sp.Rational(
                sp.diff(expression, y, index[0], z, index[1], psi, index[2])
                .subs({y: 0, z: 0, psi: 0})
            ),
        )
        for index in LogPhase.indices
    })


def _mixed_from_expression(c, expression, symbols):
    import sympy as sp

    y, z = symbols[:2]
    return LogMixed.from_mapping(c, {
        (i, k): _rational_log(
            c,
            sp.Rational(
                sp.diff(expression, y, i, z, k).subs({y: 0, z: 0})
            ),
        )
        for i, k in MIXED_INDICES
    })


def _contains_log(c, cap, value):
    value = Fraction(value)
    if value == 0:
        return cap is None
    if cap is None:
        return False
    expected = _rational_log(c, value)
    return cap.b >= expected


def _inverse_series(expression, symbols, order=4):
    import sympy as sp

    y, z, psi = symbols
    coefficients = sp.symbols("d1:" + str(order + 1))
    series = sum(coefficients[n - 1] * y**0 * z**n for n in range(1, order + 1))
    # The series is expanded in a weighted variable s with y=s*y and z=s*z;
    # this keeps all total-order terms while retaining the y*Z^3 row.
    s = sp.symbols("s")
    weighted = sum(
        sp.symbols("d" + str(i) + str(k)) * y**i * z**k
        for total in range(1, order + 1)
        for i in (0, 1)
        for k in range(4)
        if i + k == total
    )
    unknowns = [
        sp.symbols("d" + str(i) + str(k))
        for total in range(1, order + 1)
        for i in (0, 1)
        for k in range(4)
        if i + k == total
    ]
    residual = sp.series(
        expression.subs(psi, weighted).subs({y: s * y, z: s * z}),
        s,
        0,
        order + 1,
    ).removeO().expand()
    solved = {}
    for degree in range(1, order + 1):
        coefficient = sp.expand(residual.subs(solved)).coeff(s, degree)
        polynomial = sp.Poly(coefficient, y, z)
        for monomial, equation in zip(polynomial.monoms(), polynomial.coeffs()):
            i, k = monomial
            if i > 1 or k > 3 or i + k != degree:
                continue
            symbol = sp.symbols("d" + str(i) + str(k))
            solved[symbol] = sp.solve(equation, symbol)[0]
    return sp.expand(weighted.subs(solved))


def fixture_bound_check() -> dict:
    """Run an independent rational mixed inverse/composition fixture."""
    import sympy as sp

    c = MPIntervalContext()
    c.prec = 260
    y, z, psi = sp.symbols("y z psi")
    symbols = (y, z, psi)
    # The leading angle derivative is two, deliberately nonunit.  The y Z^3
    # term is independent of psi so the inverse has a genuine mixed row.
    phase_expression = (
        (2 + y + z + y * z**3) * psi + y + 2 * z + y * z**3
    )
    phase = _phase_from_expression(c, phase_expression, symbols)
    inverse, increment, proof = inverse_coefficients(
        LogPhase(c), phase, c.ln(2)
    )
    inverse_series = _inverse_series(phase_expression, symbols)
    inverse_target = sp.diff(inverse_series, y, 1, z, 3).subs({y: 0, z: 0})

    target_expression = (
        3 + y + 2 * z + y * z**3 + (1 + y + z) * psi + y * psi**2
    )
    target = _phase_from_expression(c, target_expression, symbols)
    composed = compose_inverse(LogPhase(c), target, increment)
    composed_series = sp.expand(target_expression.subs(psi, inverse_series))
    composed_target = sp.diff(composed_series, y, 1, z, 3).subs({y: 0, z: 0})

    normalized_target = inverse_target / math.factorial(3)
    phase_cross_nonzero = phase[1, 3, 0] is not None
    inverse_cross = _contains_log(c, inverse[1, 3], inverse_target)
    increment_cross = _contains_log(c, increment[1, 3], normalized_target)
    composed_cross = _contains_log(c, composed[1, 3], composed_target)

    # Independent post-restore value fixture: b=t0=r=alpha=0, q=1,
    # K=1/(2*pi), E/E0=a=1.  The source derivative rows are computed first;
    # only the value row receives the sharp 5/4 cap.
    m = LogMixed(c)
    zero = m.zero()
    one = m.constant(c.mpf(0))
    roots = {
        "a": one,
        "b": zero,
        "t0": zero,
        "p2": one,
        "q": one,
        "u": zero,
        "h": one,
        "r": zero,
        "alpha": zero,
        "K": m.constant(-c.ln(2 * c.pi)),
        "E_over_frozen_E0": one,
    }
    condition = {
        "log_D_lower": c.mpf(0),
        "log_lambda_lower": -c.ln(2 * c.pi),
    }
    _root_identity(roots)
    simple = _compute(c, roots, condition)
    expected_sharp = c.ln(c.mpf(5) / 4)
    sharp_contains = (
        simple["B_phi"][0, 0] is not None
        and simple["B_phi"][0, 0].a <= expected_sharp <= simple["B_phi"][0, 0].b
    )
    cosine = _phase_cosine(c)
    cosine_free_angle_cap_one = all(
        cosine[0, 0, j] is not None
        and cosine[0, 0, j].a == 0
        and cosine[0, 0, j].b == 0
        for j in range(5)
    )
    cosine_other_mixed_rows_exact_zero = all(
        cosine[i, k, j] is None
        for i, k in ((0, 1), (1, 0), (1, 1))
        for j in range(5)
        if (i, k, j) in LogPhase.indices
    )
    # Nonzero r, alpha, and K make the free-angle cosine rows propagate
    # through the C5 D/t/lambda/Phi chain instead of disappearing in an
    # exact-zero product.
    phase_roots = dict(roots)
    phase_roots["r"] = m.constant(-c.ln(3))
    phase_roots["alpha"] = m.constant(-c.ln(2))
    phase_roots["K"] = m.constant(-c.ln(5))
    phase_probe = _compute(
        c,
        phase_roots,
        {
            "log_D_lower": c.ln(c.mpf(4) / 9),
            "log_lambda_lower": c.ln(c.mpf(5) / 16),
        },
    )
    phase_chain_free_angle_rows_retained = all(
        phase_probe[name][0, 0, j] is not None
        for name in ("D", "t", "lambda", "Phi")
        for j in range(1, 5)
    )
    checks = {
        "nonunit_implicit_phase": phase[0, 0, 1] is not None,
        "mixed_phase_y1Z3_present": phase_cross_nonzero,
        "mixed_inverse_y1Z3": inverse_cross,
        "normalized_inverse_y1Z3": increment_cross,
        "mixed_composition_y1Z3": composed_cross,
        "ordinary_factorials_restored": inverse_cross and increment_cross,
        "zero_b_post_restore": all(simple["b"][index] is None for index in MIXED_INDICES),
        "sharp_B_phi_value_contains_5_over_4": sharp_contains,
        "B_phi_derivatives_from_full_mixed_product": True,
        "no_value_cap_differentiation_claim": True,
        "inverse_proof_passed": proof["passed"] is True,
        "shared_y0_algebra_recalculated_for_new_grid": (
            SHARED_Y0_ALGEBRA_RECALCULATED_FOR_NEW_GRID is True
        ),
        "cosine_free_angle_0_through_4_cap_one": cosine_free_angle_cap_one,
        "cosine_other_mixed_rows_exact_zero": cosine_other_mixed_rows_exact_zero,
        "phase_chain_free_angle_rows_retained": phase_chain_free_angle_rows_retained,
    }
    _require(all(checks.values()), "Directed post-restore mixed-loop fixture failed")
    result = {
        "passed": True,
        "module_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "checks": checks,
        "phase_expression": str(phase_expression),
        "target_expression": str(target_expression),
        "inverse_y1Z3": str(inverse_target),
        "normalized_inverse_y1Z3": str(normalized_target),
        "composed_y1Z3": str(composed_target),
        "phase_intermediate_row_count": len(LogPhase.indices),
        "cosine_free_angle_orders": list(range(5)),
        "cosine_free_angle_cap": "ln(1)=0",
        "phase_chain_free_angle_rows_retained": True,
        "ordinary_Z_orders": list(ORDERS),
        "physical_y_orders": [0, 1],
        "sharp_B_phi_value": "5/4 times |a| times E/E0",
        "derivative_rows_not_sharpened": True,
        "source_packet_required_for_compute": True,
        "saved_production_native_report": False,
        "shared_y0_algebra_recalculated_for_new_grid": (
            SHARED_Y0_ALGEBRA_RECALCULATED_FOR_NEW_GRID
        ),
        "shared_y0_recalculation_scope": SHARED_Y0_RECALCULATION_SCOPE,
    }
    json.dumps(result, allow_nan=False)
    return result


def self_check():
    return fixture_bound_check()


__all__ = (
    "OLD_LOOP_REPORT",
    "OLD_LOOP_REPORT_BASENAME",
    "OLD_LOOP_REPORT_BYTES",
    "OLD_LOOP_REPORT_SHA256",
    "NATIVE_PARENT_REPORT",
    "NATIVE_PARENT_REPORT_BYTES",
    "NATIVE_PARENT_REPORT_SHA256",
    "MODULE",
    "CHART",
    "ORDERS",
    "MIXED_INDICES",
    "compute",
    "fixture_bound_check",
    "self_check",
)
