"""Bounded pure-native mixed controls for the post-restore chart.

The input to compute is a byte-for-byte source receipt.  Its post-restore O4
rows are decoded as logarithmic absolute majorants and only new physical-y
rows are formed here.  Accepted ordinary-Z and y0 pressure/root rows are
consumed as data.

All arithmetic is directed log-majorant arithmetic.  None is the only
exact-zero sentinel.  No physical value is selected from a logarithmic
record, and no stored logarithm is exponentiated.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path

from mpmath import make_mpf
from mpmath.ctx_iv import MPIntervalContext

from .original_core_signed_direction import read, record, require
from .original_directed_abs_log import abs_log_upper
from .original_directed_log_reciprocal import directed_log_reciprocal
from .original_log_mixed_jets import LogMixed


REPORT_PATH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json"
)
REPORT_BYTES = 42177499
REPORT_SHA256 = "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"
MODULE = "src/openai_ns_reconstruction/original_directed_post_restore_mixed_native.py"
C6_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
)
C6_BRIDGE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
)
Z_ORDERS = tuple(range(4))
O5_ORDERS = tuple(range(5))
MIXED_INDICES = tuple((i, k) for i in (0, 1) for k in Z_ORDERS)
NATIVE_NAMES = ("phi", "U", "M", "A", "B", "C")
ROOT_NAMES = (
    "a", "b", "t0", "p2", "Delta", "q", "u", "h", "r", "alpha", "K",
    "E_over_frozen_E0",
)


def _require(value, message):
    if not value:
        raise ValueError(message)


def _lo(value):
    return make_mpf(value._mpi_[0])


def _hi(value):
    return make_mpf(value._mpi_[1])


def _record_read(context, value, label="directed logarithmic record"):
    if hasattr(value, "_mpi_"):
        out = value
    else:
        if isinstance(value, Mapping) and "value" in value and \
                "lower_exact_mpf_tuple" not in value:
            value = value["value"]
        _require(isinstance(value, Mapping), label + " is not a mapping")
        _require(value.get("directed") is True,
                 label + " must carry directed=True")
        try:
            low = make_mpf(tuple(value["lower_exact_mpf_tuple"]))
            high = make_mpf(tuple(value["upper_exact_mpf_tuple"]))
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("Malformed " + label) from exc
        out = context.mpf([low, high])
    try:
        low = make_mpf(out._mpi_[0])
        high = make_mpf(out._mpi_[1])
    except (AttributeError, TypeError, ValueError) as exc:
        raise ValueError("Malformed " + label) from exc
    _require(low <= high, "Reversed " + label)
    return out


def _logsum(context, terms):
    """Add absolute logarithmic caps without materializing exponentials."""
    terms = tuple(term for term in terms if term is not None)
    if not terms:
        return None
    if len(terms) == 1:
        return terms[0]
    shift = context.ln(context.mpf(len(terms)))
    # Form each interval sum in the caller context first.  Extracting scalar
    # endpoints and adding them would route through the process-wide scalar
    # precision and can underbound the enormous saved logarithms.
    shifted = [term + shift for term in terms]
    lower = min(_lo(term) for term in shifted)
    upper = max(_hi(term) for term in shifted)
    return context.mpf([lower, upper])


def _log_constant(context, value):
    """Absolute log cap for an exact rational/integer or interval scalar."""
    if isinstance(value, Fraction):
        value = context.mpf(value.numerator) / value.denominator
    elif type(value) is int:
        value = context.mpf(value)
    elif isinstance(value, str):
        value = context.mpf(value)
    _require(hasattr(value, "_mpi_"), "A directed scalar is required")
    return abs_log_upper(context, value)


def _scale_rows(context, rows, scalar):
    log_scalar = _log_constant(context, scalar)
    return [None if row is None else row + log_scalar for row in rows]


def _z5_add(context, *rows):
    _require(rows, "At least one O5 row vector is required")
    _require(all(len(item) >= 5 for item in rows),
             "O5 row vector is incomplete")
    return [_logsum(context, [item[k] for item in rows])
            for k in O5_ORDERS]


def _z5_mul(context, left, right):
    _require(len(left) >= 5 and len(right) >= 5,
             "O5 product inputs are incomplete")
    result = []
    for k in O5_ORDERS:
        terms = []
        for j in range(k + 1):
            left_row, right_row = left[j], right[k - j]
            if left_row is None or right_row is None:
                continue
            coefficient = math.comb(k, j)
            term = left_row + right_row
            if coefficient != 1:
                term += _log_constant(context, coefficient)
            terms.append(term)
        result.append(_logsum(context, terms))
    return result


def _z5_scale(context, rows, scalar):
    log_scalar = _log_constant(context, scalar)
    return [None if row is None else row + log_scalar for row in rows]


def _mixed_from_rows(context, y0, y1):
    _require(len(y0) >= 4 and len(y1) >= 4,
             "Mixed rows must retain Z0 through Z3")
    return LogMixed.from_mapping(
        context,
        {(i, k): (y0[k] if i == 0 else y1[k])
         for i in (0, 1) for k in Z_ORDERS},
    )


def _mixed_from_z5(context, rows):
    return _mixed_from_rows(context, rows, [None] * 5)


def _mixed_y_zero(context, jet):
    return LogMixed.from_mapping(
        context,
        {(i, k): (jet[0, k] if i == 0 else None)
         for i in (0, 1) for k in Z_ORDERS},
    )


def _mixed_y_only(context, jet):
    return LogMixed.from_mapping(
        context,
        {(i, k): (jet[1, k] if i == 1 else None)
         for i in (0, 1) for k in Z_ORDERS},
    )


def _mixed_z_shift(context, y0, y1):
    return _mixed_from_rows(context, y0[1:5], y1[1:5])


def _pack_mixed(jet):
    return [
        dict(
            physical_y_order=i,
            ordinary_Z_order=k,
            exact_zero=jet[i, k] is None,
            log_absolute_upper_bound=(
                None if jet[i, k] is None else record(jet[i, k])
            ),
        )
        for i, k in MIXED_INDICES
    ]


def _pack_z5(rows):
    _require(len(rows) >= 5, "Five ordinary-Z rows are required")
    return [
        dict(
            ordinary_Z_order=k,
            exact_zero=rows[k] is None,
            log_absolute_upper_bound=(
                None if rows[k] is None else record(rows[k])
            ),
        )
        for k in O5_ORDERS
    ]


def _pack_z4(rows):
    _require(len(rows) >= 4, "Four ordinary-Z rows are required")
    return [
        dict(
            ordinary_Z_order=k,
            exact_zero=rows[k] is None,
            log_absolute_upper_bound=(
                None if rows[k] is None else record(rows[k])
            ),
        )
        for k in Z_ORDERS
    ]


def _read_packed_rows(context, rows, count, label):
    _require(isinstance(rows, list) and len(rows) >= count,
             "Missing " + label)
    result = []
    for expected, row in enumerate(rows[:count]):
        _require(isinstance(row, Mapping), label + " row is not a mapping")
        _require(row.get("ordinary_Z_order") == expected,
                 "Misordered " + label)
        exact = row.get("exact_zero") is True
        raw = row.get("log_absolute_upper_bound")
        if exact:
            _require(raw is None, label + " exact-zero row carries a record")
            result.append(None)
        else:
            _require(isinstance(raw, Mapping),
                     label + " nonzero row lacks a record")
            result.append(_record_read(context, raw, label))
    return result


def _read_direct_rows(context, rows, count, label):
    _require(isinstance(rows, list) and len(rows) >= count,
             "Missing " + label)
    return [_record_read(context, row, label) for row in rows[:count]]


def _source_root():
    return Path(__file__).resolve().parents[2]


def _validate_primary_sources(root):
    """Verify the exact C6 post identity and bridge density identities."""
    c6_path, bridge_path = root / C6_SOURCE, root / C6_BRIDGE
    _require(c6_path.is_file() and bridge_path.is_file(),
             "C6 primary source files are missing")
    c6 = c6_path.read_text(encoding="utf-8")
    bridge = bridge_path.read_text(encoding="utf-8")
    compact_c6 = "".join(c6.split())
    compact_bridge = "".join(bridge.split())
    c6_fragment = (
        "actual=dict(phi=phi,U=a.scale(z,4),"
        "phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0)"
    )
    density_fragment = "A=a.mul(U,U),B=a.mul(phi,phi),C=a.mul(phi,phi)"
    _require("".join(c6_fragment.split()) in compact_c6,
             "C6 post-restore physical-y identity changed")
    _require("".join(density_fragment.split()) in compact_bridge,
             "C6 bridge density identity changed")
    return dict(
        C6_post_restore=dict(
            path=C6_SOURCE,
            source_sha256=hashlib.sha256(c6_path.read_bytes()).hexdigest(),
            required_text=[c6_fragment],
            source_executed=False,
            passed=True,
        ),
        C6_bridge_densities=dict(
            path=C6_BRIDGE,
            source_sha256=hashlib.sha256(bridge_path.read_bytes()).hexdigest(),
            required_text=[density_fragment],
            source_executed=False,
            passed=True,
        ),
    )


def _validate_report(native_report_bytes):
    """Check bytes/hash before decoding or parsing the JSON report."""
    if not isinstance(native_report_bytes, (bytes, bytearray)):
        raise TypeError("native_report_bytes must be bytes")
    raw = bytes(native_report_bytes)
    _require(len(raw) == REPORT_BYTES,
             "Native report byte length changed before JSON parsing")
    digest = hashlib.sha256(raw).hexdigest()
    _require(digest == REPORT_SHA256,
             "Native report SHA256 changed before JSON parsing")
    try:
        payload = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("Native report JSON is malformed") from exc
    _require(isinstance(payload, Mapping), "Native report root is not a mapping")
    version = payload.get("actual_native_directed_controls_v2")
    _require(isinstance(version, Mapping),
             "actual_native_directed_controls_v2 is missing")
    _require(version.get("passed") is True, "Native report did not pass")
    audit = version.get("execution_audit")
    _require(isinstance(audit, Mapping), "Native report execution audit is missing")
    for key in (
        "accepted_builders_or_owners_replayed",
        "accepted_checkers_replayed",
        "accepted_self_checks_replayed",
        "accepted_recurrences_replayed",
        "accepted_quadratures_replayed",
        "old_stage_calls",
        "source_executed",
    ):
        _require(audit.get(key) is False,
                 "Native report execution audit is not source-only: " + key)
    raw_tree = version.get("raw_interval_tree")
    _require(isinstance(raw_tree, Mapping), "Native raw interval tree is missing")
    shared = raw_tree.get("shared")
    restoration = raw_tree.get("restoration")
    _require(isinstance(shared, Mapping) and isinstance(restoration, Mapping),
             "Native shared/restoration interval trees are missing")
    charts = restoration.get("charts")
    post = charts.get("post_restore") if isinstance(charts, Mapping) else None
    _require(isinstance(post, Mapping), "Native post_restore chart is missing")
    _require(post.get("passed") is True and post.get("chart") == "post_restore",
             "Native post_restore chart did not pass")
    _require(post.get("ordinary_Z_orders") == list(O5_ORDERS),
             "Native post_restore ordinary Z order changed")
    _require(post.get("p2_ordinary_Z_orders") == list(Z_ORDERS),
             "Native post_restore p2 Z order changed")
    return payload, version, raw_tree, shared, post, digest


def _native_mixed(context, post):
    native_y0 = {
        name: _read_packed_rows(
            context, post["native_caps"][name], 5, "post_restore " + name
        )
        for name in NATIVE_NAMES
    }
    phi, U, M, A, B, C = (native_y0[name] for name in NATIVE_NAMES)
    phi_y = _z5_scale(context, phi, Fraction(2, 5))
    U_y = [None] * 5
    U2 = _z5_mul(context, U, U)
    phi2 = _z5_mul(context, phi, phi)
    M_y = _z5_add(context, U, M)
    A_y = _z5_add(context, U2, A)
    # These are signed source identities.  Absolute logarithmic caps add
    # magnitudes, while the packet retains the signed formulas as metadata.
    B_y = _z5_add(context, phi2, _z5_scale(context, B, 2))
    C_y = _z5_add(context, phi2, C)
    native_y = dict(phi=phi_y, U=U_y, M=M_y, A=A_y, B=B_y, C=C_y)
    mixed = {
        name: _mixed_from_rows(context, native_y0[name], native_y[name])
        for name in NATIVE_NAMES
    }
    return native_y0, native_y, mixed


def _z_cap(context, shared):
    z_interval = _record_read(context, shared["z"], "shared z interval")
    z0 = abs_log_upper(context, z_interval)
    zero = context.mpf(0)
    rows = [z0, zero, None, None, None]
    y_rows = [None] * 5
    return _mixed_from_rows(context, rows, y_rows), rows, z_interval


def _pressure(context, shared, post, native_y0, native_y, native_mixed):
    engine = LogMixed(context)
    z, z_rows, z_interval = _z_cap(context, shared)
    d = engine.add(engine.constant(context.mpf(0)), engine.mul(z, z))
    delta = _record_read(context, shared["delta"], "shared delta")
    delta_log = abs_log_upper(context, delta)
    one_minus_delta_log = abs_log_upper(context, context.mpf(1) - delta)
    one_plus_delta_log = abs_log_upper(context, context.mpf(1) + delta)
    half_one_minus_delta = one_minus_delta_log + \
        _log_constant(context, Fraction(1, 2))
    twice_delta = delta_log + _log_constant(context, 2)
    twice_one_plus_delta = one_plus_delta_log + _log_constant(context, 2)

    U, M, A, B, C = (native_mixed[name] for name in ("U", "M", "A", "B", "C"))
    Mz = _mixed_z_shift(context, native_y0["M"], native_y["M"])
    Az = _mixed_z_shift(context, native_y0["A"], native_y["A"])
    Bz = _mixed_z_shift(context, native_y0["B"], native_y["B"])
    Cz = _mixed_z_shift(context, native_y0["C"], native_y["C"])
    F2_rows = _read_direct_rows(
        context, shared["F2_logs5"], 5, "shared F2"
    )
    F2 = _mixed_from_z5(context, F2_rows)
    P0_rows = _read_direct_rows(
        context, shared["P0_logs5"], 5, "shared P0"
    )
    P0 = _mixed_from_z5(context, P0_rows)
    P0z = _mixed_from_rows(context, P0_rows[1:5], [None] * 5)

    term_M = engine.scale(engine.add(M, engine.mul(z, Mz)),
                          half_one_minus_delta)
    term_hydro = engine.mul(
        U,
        engine.add(
            engine.scale(engine.mul(z, M), one_minus_delta_log),
            engine.mul(d, Mz),
        ),
    )
    term_A = engine.add(
        engine.scale(engine.mul(z, A), twice_delta),
        engine.mul(d, Az),
    )
    term_P = engine.add(
        engine.scale(engine.mul(z, P0), twice_one_plus_delta),
        engine.mul(d, P0z),
    )
    swirl = engine.add(
        engine.scale(engine.mul(z, B), twice_delta),
        engine.mul(d, Bz),
        engine.scale(engine.mul(z, C), twice_one_plus_delta),
        engine.mul(d, Cz),
    )
    log_R = _record_read(
        context, post["scales"]["log_R_upper"], "post_restore log_R_upper"
    )
    R = LogMixed.from_mapping(
        context,
        {(i, k): (log_R if k == 0 else None)
         for i in (0, 1) for k in Z_ORDERS},
    )
    term_swirl = engine.mul(R, engine.mul(F2, swirl))
    Nnative = engine.add(U, term_M, term_hydro, term_A, term_P, term_swirl)

    inverse_L_rows = _read_direct_rows(
        context, shared["inverse_L_logs5"], 5, "shared inverse_L"
    )
    inverse_L = _mixed_from_z5(context, inverse_L_rows)
    inverse_Utheta_rows = _read_packed_rows(
        context, post["inverse_Utheta"], 5, "post_restore inverse_Utheta"
    )
    inverse_Utheta_y = _scale_rows(
        context, inverse_Utheta_rows, Fraction(1, 10)
    )
    inverse_Utheta = _mixed_from_rows(
        context, inverse_Utheta_rows, inverse_Utheta_y
    )
    R0 = _mixed_y_zero(context, R)
    inv_Utheta0 = _mixed_y_zero(context, inverse_Utheta)
    base = engine.mul(R0, engine.mul(inv_Utheta0, inverse_L))
    N_y_only = _mixed_y_only(context, Nnative)
    p2_saved_rows = _read_packed_rows(
        context, post["p2"], 4, "post_restore p2"
    )
    p2_saved = _mixed_from_rows(context, p2_saved_rows, [None] * 4)
    p2_from_R_E = engine.scale(
        p2_saved, _log_constant(context, Fraction(9, 10))
    )
    p2_from_Ny = engine.mul(base, N_y_only)
    p2_y_rows = [
        _logsum(context, (p2_from_R_E[0, k], p2_from_Ny[1, k]))
        for k in Z_ORDERS
    ]
    p2_mixed = LogMixed.from_mapping(
        context,
        {(0, k): p2_saved[0, k] for k in Z_ORDERS}
        | {(1, k): p2_y_rows[k] for k in Z_ORDERS},
    )
    e_rows = _read_packed_rows(
        context, post["E_relative"], 5, "post_restore E_relative"
    )
    e_y = _scale_rows(context, e_rows, Fraction(1, 10))
    E = _mixed_from_rows(context, e_rows, e_y)
    return dict(
        z=z,
        z_rows=z_rows,
        z_interval=z_interval,
        delta=delta,
        R=R,
        F2=F2,
        F2_rows=F2_rows,
        P0=P0,
        P0_rows=P0_rows,
        inverse_L=inverse_L,
        inverse_Utheta=inverse_Utheta,
        Nnative=Nnative,
        Nnative_y=N_y_only,
        p2_saved=p2_saved,
        p2_mixed=p2_mixed,
        p2_y_rows=p2_y_rows,
        p2_from_R_E=p2_from_R_E,
        p2_from_Ny=p2_from_Ny,
        E=E,
        source_formula=(
            "p2=R*Nnative/(Utheta*L); "
            "Nnative=U+(1-delta)*(M+Z*M_Z)/2"
            "+U*((1-delta)*Z*M+(1-Z^2)*M_Z)"
            "+2*delta*Z*A-(1-Z^2)*A_Z"
            "+2*(1+delta)*Z*P0-(1-Z^2)*P0_Z"
            "+R*F0^2*(2*delta*Z*B-(1-Z^2)*B_Z"
            "+2*(1+delta)*Z*C-(1-Z^2)*C_Z)"
        ),
        p2_y_formula="p2_y=(9/10)*p2+R*Nnative_y/(Utheta*L)",
        physical_y_identities=dict(
            phi_y="-2/5*phi",
            U_y="0",
            M_y="U-M",
            A_y="U^2-A",
            B_y="phi^2-2*B",
            C_y="phi^2-C",
            E_y="E/10",
            R_y="R",
        ),
        inverse_Utheta_formula=(
            "inverse_Utheta=(1+Z^2)*exp(-logP-offset/10)"
        ),
        inverse_Utheta_y_formula="inverse_Utheta_y=-1/10*inverse_Utheta",
        p2_y_R_E_term=_pack_mixed(p2_from_R_E),
        p2_y_Nnative_term=_pack_mixed(p2_from_Ny),
    )


def _copy_root_y0(context, post, field):
    return _read_packed_rows(
        context, post["root_caps"][field], 4, "post_restore root " + field
    )


def _root_mixed(context, shared, post, pressure):
    engine = LogMixed(context)
    zero = engine.zero()
    a = engine.constant(_log_constant(context, Fraction(4, 5)))
    b = zero
    t0 = zero
    Delta = engine.constant(_log_constant(context, Fraction(6, 5)))
    q_y0 = _copy_root_y0(context, post, "q")
    q = _mixed_from_rows(context, q_y0, [None] * 4)
    q2 = engine.mul(q, q)

    dstar_log = _record_read(context, shared["logdstar"], "shared logdstar")
    eta_log = _record_read(context, shared["logeta"], "shared logeta")
    dstar_rows = [dstar_log, None, None, None, None]
    inverse_dstar_rows = directed_log_reciprocal(
        context, dstar_rows, dstar_log
    )
    inverse_dstar = _mixed_from_z5(context, inverse_dstar_rows)
    p2 = pressure["p2_mixed"]
    u_calc = engine.mul(engine.mul(p2, q), inverse_dstar)
    u2 = engine.mul(u_calc, u_calc)
    one_plus_u2 = engine.add(engine.constant(context.mpf(0)), u2)
    # ``one_plus_u2[0, 0]`` is a logarithmic cap interval.  Its lower
    # endpoint can be negative when u is tiny, while its upper endpoint is
    # still a certified upper bound for log(1 + u^2).  The generic sqrt
    # contract needs a positive lower bound for its *upper-bound argument*,
    # so rebox that certified upper endpoint as a point interval.  Do not
    # pass the cap's lower endpoint as a purported positive floor.
    one_plus_u2_upper = _hi(one_plus_u2[0, 0])
    sqrt_one_plus_u2 = engine.sqrt(
        one_plus_u2,
        context.mpf(0),
        context.mpf([one_plus_u2_upper, one_plus_u2_upper]),
    )
    h_calc = engine.reciprocal(sqrt_one_plus_u2, context.mpf(0))
    r_calc = engine.mul(u_calc, h_calc)
    alpha_calc = engine.scale(
        engine.mul(q, h_calc), _log_constant(context, 2)
    )
    ratio = engine.add(
        engine.constant(context.mpf(0)),
        engine.scale(q2, _log_constant(context, 2)),
    )
    ratio_lower = context.ln(context.mpf(5) / 2)
    K_calc = engine.scale(
        engine.reciprocal(ratio, ratio_lower),
        context.ln(context.mpf(1) / (2 * context.pi)),
    )

    def replace_y1(field, calculated):
        y0 = _copy_root_y0(context, post, field)
        return _mixed_from_rows(
            context, y0, [calculated[1, k] for k in Z_ORDERS]
        )

    root_caps = {
        "a": a,
        "b": b,
        "t0": t0,
        "p2": p2,
        "Delta": Delta,
        "q": q,
        "u": replace_y1("u", u_calc),
        "h": replace_y1("h", h_calc),
        "r": replace_y1("r", r_calc),
        "alpha": replace_y1("alpha", alpha_calc),
        "K": replace_y1("K", K_calc),
        "E_over_frozen_E0": pressure["E"],
    }
    root_row_provenance = {
        "a": dict(
            y0="recomputed_exact_constant",
            source_identity="4/5",
            source_record=copy.deepcopy(post["root_caps"]["a"]),
            y1="exact_zero_identity",
        ),
        "b": dict(
            y0="recomputed_exact_zero_identity",
            source_record=copy.deepcopy(post["root_caps"]["b"]),
            y1="exact_zero_identity",
        ),
        "t0": dict(
            y0="recomputed_exact_zero_identity",
            source_record=copy.deepcopy(post["root_caps"]["t0"]),
            y1="exact_zero_identity",
        ),
        "p2": dict(
            y0="copied_accepted_post_p2",
            source_record=copy.deepcopy(post["p2"]),
            y1="new_generic_pressure_product",
        ),
        "Delta": dict(
            y0="recomputed_exact_constant",
            source_identity="-6/5",
            source_record=copy.deepcopy(post["root_caps"]["Delta"]),
            y1="exact_zero_identity",
        ),
        "q": dict(
            y0="copied_accepted_root_q",
            source_record=copy.deepcopy(post["root_caps"]["q"]),
            y1="exact_zero_identity",
        ),
        "u": dict(
            y0="copied_accepted_root_u",
            source_record=copy.deepcopy(post["root_caps"]["u"]),
            y1="new_generic_LogMixed",
        ),
        "h": dict(
            y0="copied_accepted_root_h",
            source_record=copy.deepcopy(post["root_caps"]["h"]),
            y1="new_generic_LogMixed",
        ),
        "r": dict(
            y0="copied_accepted_root_r",
            source_record=copy.deepcopy(post["root_caps"]["r"]),
            y1="new_generic_LogMixed",
        ),
        "alpha": dict(
            y0="copied_accepted_root_alpha",
            source_record=copy.deepcopy(post["root_caps"]["alpha"]),
            y1="new_generic_LogMixed",
        ),
        "K": dict(
            y0="copied_accepted_root_K",
            source_record=copy.deepcopy(post["root_caps"]["K"]),
            y1="new_generic_LogMixed",
        ),
        "E_over_frozen_E0": dict(
            y0="copied_post_E_relative",
            source_record=copy.deepcopy(post["E_relative"]),
            y1="new_physical_y_scale",
        ),
    }
    return dict(
        root_caps=root_caps,
        root_row_provenance=root_row_provenance,
        eta_log=eta_log,
        q_squared=q2,
        dstar_log=dstar_log,
        inverse_dstar=inverse_dstar,
        ratio=ratio,
        one_plus_u2=one_plus_u2,
        sqrt_one_plus_u2=sqrt_one_plus_u2,
        positive_floors=dict(
            ratio_lower_log=record(ratio_lower),
            source_conditioning=copy.deepcopy(post["conditioning"]),
        ),
        identities=dict(
            a="4/5",
            b="0",
            t0="0",
            Delta="-6/5",
            q_squared="(2*eta-Delta)/(2*a)",
            q_squared_post_y_and_Z_constant=True,
            eta_source="raw_interval_tree.shared.logeta",
            eta_handle_used=False,
            dstar_source="raw_interval_tree.shared.logdstar",
            dstar_handle_used=False,
            u="p2*q/dstar",
            h="(1+u^2)^(-1/2)",
            r="u*h",
            alpha="2*q*h",
            K="1/[2*pi*(1+t0^2+2*q^2)]",
            y1_from_generic_LogMixed=True,
        ),
    )


def _provenance(version, shared, post, digest, source_checks):
    return dict(
        input_report_path=REPORT_PATH,
        input_report_bytes=REPORT_BYTES,
        input_report_sha256=digest,
        source_native_report_sha256=digest,
        source_guards=copy.deepcopy(version.get("source_guards", {})),
        execution_audit=copy.deepcopy(version.get("execution_audit", {})),
        primary_source_checks=source_checks,
        source_executed=False,
        accepted_builders_or_owners_replayed=False,
        accepted_checkers_replayed=False,
        accepted_recurrences_replayed=False,
        accepted_quadratures_replayed=False,
        root_record_provenance=copy.deepcopy(post["root_caps"]),
        root_input_record_provenance=copy.deepcopy(post["root_inputs"]),
        native_record_provenance=copy.deepcopy(post["native_caps"]),
        post_inverse_Utheta_record_provenance=copy.deepcopy(post["inverse_Utheta"]),
        post_p2_record_provenance=copy.deepcopy(post["p2"]),
        post_conditioning_record_provenance=copy.deepcopy(post["conditioning"]),
        shared_P0_record_provenance=copy.deepcopy(shared["P0_logs5"]),
        shared_F2_record_provenance=copy.deepcopy(shared["F2_logs5"]),
        shared_inverse_L_record_provenance=copy.deepcopy(shared["inverse_L_logs5"]),
        shared_eta_record_provenance=copy.deepcopy(shared["logeta"]),
        shared_dstar_record_provenance=copy.deepcopy(shared["logdstar"]),
    )


def _jsonable(value):
    if hasattr(value, "_mpi_"):
        return record(value)
    if isinstance(value, Mapping):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def compute(context, native_report_bytes):
    """Return the bounded post-restore mixed-native packet."""
    _require(isinstance(context, MPIntervalContext),
             "context must be an MPIntervalContext")
    payload, version, raw_tree, shared, post, digest = \
        _validate_report(native_report_bytes)
    source_checks = _validate_primary_sources(_source_root())
    native_y0, native_y, native_mixed = _native_mixed(context, post)
    pressure = _pressure(
        context, shared, post, native_y0, native_y, native_mixed
    )
    roots = _root_mixed(context, shared, post, pressure)

    native_caps = {
        name: _pack_mixed(native_mixed[name]) for name in NATIVE_NAMES
    }
    native_y_caps = {
        name: _pack_z5(native_y[name]) for name in NATIVE_NAMES
    }
    root_mixed_caps = {
        name: _pack_mixed(roots["root_caps"][name]) for name in ROOT_NAMES
    }
    conditioning = copy.deepcopy(post["conditioning"])
    provenance = _provenance(version, shared, post, digest, source_checks)
    input_record_provenance = dict(
        root_caps=copy.deepcopy(post["root_caps"]),
        root_inputs=copy.deepcopy(post["root_inputs"]),
        native_caps=copy.deepcopy(post["native_caps"]),
        inverse_Utheta=copy.deepcopy(post["inverse_Utheta"]),
        p2=copy.deepcopy(post["p2"]),
        conditioning=copy.deepcopy(post["conditioning"]),
        shared=dict(
            P0_logs5=copy.deepcopy(shared["P0_logs5"]),
            F2_logs5=copy.deepcopy(shared["F2_logs5"]),
            inverse_L_logs5=copy.deepcopy(shared["inverse_L_logs5"]),
            logeta=copy.deepcopy(shared["logeta"]),
            logdstar=copy.deepcopy(shared["logdstar"]),
            delta=copy.deepcopy(shared["delta"]),
            z=copy.deepcopy(shared["z"]),
        ),
        source_native_report_sha256=digest,
        source_executed=False,
    )
    packet = dict(
        passed=True,
        schema="original_directed_post_restore_mixed_native_v1",
        chart="post_restore",
        coordinate="physical_y_then_ordinary_Z",
        ordinary_Z_orders=list(Z_ORDERS),
        physical_y_orders=[0, 1],
        physical_y_order_rows=True,
        source_native_report_sha256=digest,
        source_native_report_path=REPORT_PATH,
        input_report_bytes=REPORT_BYTES,
        scope=dict(
            post_only=True,
            chart="post_restore",
            mixed_yZ3=True,
            mixed_yZ3_not_full_C5=True,
            full_C5=False,
            global_NS=False,
            current_NS=False,
            accepted_y0_and_root_Z_rows_consumed=True,
            new_y1_rows_only=True,
        ),
        native_caps=native_caps,
        native_row_provenance={
            name: dict(
                y0="copied_accepted_native_cap",
                source_record=copy.deepcopy(post["native_caps"][name]),
                y1="new_physical_y_identity",
            )
            for name in NATIVE_NAMES
        },
        native_y1_internal_Z0_Z4=native_y_caps,
        native_y1_public_Z0_Z3={
            name: _pack_z4(native_y[name]) for name in NATIVE_NAMES
        },
        physical_y_identities=pressure["physical_y_identities"],
        root_mixed_caps=root_mixed_caps,
        root_caps=root_mixed_caps,
        root_q_squared_mixed_caps=_pack_mixed(roots["q_squared"]),
        conditioning=conditioning,
        root_conditioning=copy.deepcopy(conditioning),
        positive_D_lower=copy.deepcopy(conditioning["positive_D_lower"]),
        positive_lambda_lower=copy.deepcopy(conditioning["positive_lambda_lower"]),
        p2=dict(
            mixed_log_caps=_pack_mixed(pressure["p2_mixed"]),
            p2_y=_pack_z4(pressure["p2_y_rows"]),
            p2_y_Z0_Z3=_pack_z4(pressure["p2_y_rows"]),
            Nnative=_pack_mixed(pressure["Nnative"]),
            Nnative_y=_pack_mixed(pressure["Nnative_y"]),
            inverse_L=_pack_mixed(pressure["inverse_L"]),
            inverse_Utheta=_pack_mixed(pressure["inverse_Utheta"]),
            F0_squared=_pack_z5(pressure["F2_rows"]),
            P0=_pack_z5(pressure["P0_rows"]),
            radius_factor=_pack_mixed(pressure["R"]),
            E_over_frozen_E0=_pack_mixed(pressure["E"]),
            p2_y_R_E_term=pressure["p2_y_R_E_term"],
            p2_y_Nnative_term=pressure["p2_y_Nnative_term"],
            source_formula=pressure["source_formula"],
            p2_y_formula=pressure["p2_y_formula"],
            inverse_Utheta_formula=pressure["inverse_Utheta_formula"],
            inverse_Utheta_y_formula=pressure["inverse_Utheta_y_formula"],
            full_P0_and_P0_Z_retained=True,
            full_F0_squared_B_C_products_retained=True,
            radius_factor_retained=True,
            source_executed=False,
        ),
        root_identities=roots["identities"],
        root_row_provenance=roots["root_row_provenance"],
        eta_log_record=record(roots["eta_log"]),
        dstar_log_record=record(roots["dstar_log"]),
        root_positive_floors=roots["positive_floors"],
        inverse_dstar=_pack_mixed(roots["inverse_dstar"]),
        root_ratio=_pack_mixed(roots["ratio"]),
        provenance=provenance,
        input_record_provenance=input_record_provenance,
        source_executed=False,
        accepted_ancestors_replayed=False,
        accepted_builders_or_owners_replayed=False,
        accepted_checkers_replayed=False,
        accepted_recurrences_replayed=False,
        accepted_quadratures_replayed=False,
        accepted_norms_executed=False,
        full_C5=False,
        global_five_moment_closure=False,
        current_NS=False,
        endpoint_numeric_materialized=False,
    )
    result = _jsonable(packet)
    json.dumps(result, allow_nan=False, separators=(",", ":"))
    return result


def _fresh_log(context, value):
    value = Fraction(value)
    _require(value != 0, "Fresh fixture values must be nonzero")
    return context.ln(context.mpf(abs(value.numerator)) / value.denominator)


def _logsum_precision_fixture():
    """Exercise the interval-first logsum path at scalar precision 53."""
    from mpmath import mp

    previous_precision = mp.prec
    try:
        mp.prec = 53
        context = MPIntervalContext()
        context.prec = 320
        # These are exact MPF tuples with a very large exponent.  The
        # candidate interval additions are formed before endpoint extraction.
        huge_left = make_mpf((0, 1, 100000, 1))
        huge_right = make_mpf((0, 3, 99998, 2))
        terms = [
            context.mpf([huge_left, huge_left]),
            context.mpf([huge_right, huge_right]),
        ]
        shift = context.ln(context.mpf(2))
        candidates = [term + shift for term in terms]
        expected_lower = min(_lo(term) for term in candidates)
        expected_upper = max(_hi(term) for term in candidates)
        actual = _logsum(context, terms)
        passed = (_lo(actual) <= expected_lower
                  and _hi(actual) >= expected_upper)
        _require(passed, "Global scalar precision logsum fixture failed")
        return dict(
            passed=True,
            scalar_precision_bits=53,
            interval_precision_bits=context.prec,
            huge_exact_log_tuples=[
                list(terms[0]._mpi_[0]), list(terms[1]._mpi_[0])
            ],
            candidate_lower_exact_tuple=list(expected_lower._mpf_),
            candidate_upper_exact_tuple=list(expected_upper._mpf_),
            result_lower_exact_tuple=list(actual._mpi_[0]),
            result_upper_exact_tuple=list(actual._mpi_[1]),
            result_contains_candidate_outer_bounds=True,
            endpoint_scalar_arithmetic_used=False,
        )
    finally:
        mp.prec = previous_precision


def fresh_fixture():
    """Run fresh rational pressure and mixed-root fixtures."""
    import sympy as sp

    context = MPIntervalContext()
    context.prec = 260
    y, z = sp.symbols("y z")
    phi0 = 2 + 3 * z + 5 * z**2 + z**3
    U0 = 3 + z + 2 * z**2 + z**3
    M0 = 5 + 2 * z + z**2 + 2 * z**3
    A0 = 7 + z + 3 * z**2 + z**3
    B0 = 11 + 2 * z + z**2 + z**3
    C0 = 13 + z + 2 * z**2 + 3 * z**3
    phi = phi0 + y * (-sp.Rational(2, 5) * phi0)
    U = U0
    M = M0 + y * (U0 - M0)
    A = A0 + y * (U0**2 - A0)
    B = B0 + y * (phi0**2 - 2 * B0)
    C = C0 + y * (phi0**2 - C0)
    delta = sp.Rational(1, 11)
    R = 5 * (1 + y)
    F2 = 2 + z + z**2 + 3 * z**3 + 4 * z**4 + 5 * z**5
    P0 = 3 + 2 * z + z**3 + 4 * z**4 + 5 * z**5
    d = 1 + z**2
    N = (
        U
        + (1 - delta) * (M + z * sp.diff(M, z)) / 2
        + U * ((1 - delta) * z * M + d * sp.diff(M, z))
        + 2 * delta * z * A + d * sp.diff(A, z)
        + 2 * (1 + delta) * z * P0 + d * sp.diff(P0, z)
        + R * F2 * (
            2 * delta * z * B + d * sp.diff(B, z)
            + 2 * (1 + delta) * z * C + d * sp.diff(C, z)
        )
    )
    E = (4 + z) * (1 + sp.Rational(1, 10) * y)
    L = 1 - delta * z**2
    p2 = sp.cancel(R * N / (E * L))
    lhs = sp.cancel(sp.diff(p2, y).subs(y, 0))
    rhs = sp.cancel(
        sp.Rational(9, 10) * p2.subs(y, 0)
        + R.subs(y, 0) * sp.diff(N, y).subs(y, 0)
        / (E.subs(y, 0) * L)
    )
    p2_identity = sp.cancel(lhs - rhs) == 0
    p2_y_value = sp.diff(p2, y).subs({y: 0, z: sp.Rational(1, 7)})
    _require(p2_identity and p2_y_value != 0,
             "Fresh rational p2_y product fixture failed")
    derivative_nonzero = {
        name: sp.diff(expression, y).subs({y: 0, z: sp.Rational(1, 7)})
        for name, expression in (
            ("M", M), ("A", A), ("B", B), ("C", C)
        )
    }
    _require(all(value != 0 for value in derivative_nonzero.values()),
             "Fresh MABC derivative fixture is degenerate")

    p2_expr = sp.expand(p2.subs(y, 0))
    p2_y_expr = sp.expand(sp.diff(p2, y).subs(y, 0))
    p2_values = {
        (0, k): _fresh_log(
            context, Fraction(sp.diff(p2_expr, z, k).subs(z, 0))
        )
        for k in Z_ORDERS
    }
    p2_values.update({
        (1, k): _fresh_log(
            context, Fraction(sp.diff(p2_y_expr, z, k).subs(z, 0))
        )
        for k in Z_ORDERS
    })
    p2_jet = LogMixed.from_mapping(context, p2_values)
    q_jet = LogMixed(context).constant(_fresh_log(context, Fraction(3, 2)))
    inv_dstar_jet = LogMixed(context).constant(
        _fresh_log(context, Fraction(1, 2))
    )
    engine = LogMixed(context)
    u = engine.mul(engine.mul(p2_jet, q_jet), inv_dstar_jet)
    u2 = engine.mul(u, u)
    one_plus = engine.add(engine.constant(context.mpf(0)), u2)
    h = engine.reciprocal(
        engine.sqrt(one_plus, context.mpf(0), one_plus[0, 0]),
        context.mpf(0),
    )
    root_nonzero = all(
        jet[1, 1] is not None for jet in (p2_jet, u, h)
    )
    _require(root_nonzero, "Fresh root mixed Z derivative fixture failed")

    # A separate synthetic native O5 packet exercises the implementation
    # product factors.  Every source record below is made from an independent
    # signed rational SymPy derivative at z=1/5.
    z_location = sp.Rational(1, 5)

    def evaluated(expression, order):
        return sp.factor(
            sp.diff(expression, z, order).subs({y: 0, z: z_location})
        )

    def packed_caps(expression, count=5):
        rows = []
        for order in range(count):
            value = evaluated(expression, order)
            cap = None if value == 0 else _fresh_log(context, Fraction(value))
            rows.append(dict(
                ordinary_Z_order=order,
                exact_zero=cap is None,
                log_absolute_upper_bound=None if cap is None else record(cap),
            ))
        return rows

    def direct_records(expression, count=5):
        rows = []
        for order in range(count):
            value = evaluated(expression, order)
            _require(value != 0,
                     "Synthetic direct source derivative became zero")
            rows.append(record(_fresh_log(context, Fraction(value))))
        return rows

    inverse_Utheta_expression = 1 / (4 + z)
    post = dict(
        native_caps={
            name: packed_caps(expression)
            for name, expression in (
                ("phi", phi), ("U", U), ("M", M), ("A", A),
                ("B", B), ("C", C)
            )
        },
        inverse_Utheta=packed_caps(inverse_Utheta_expression),
        p2=packed_caps(p2.subs(y, 0), 4),
        E_relative=packed_caps(4 + z),
        scales=dict(
            log_R_upper=record(
                _fresh_log(context, Fraction(R.subs({y: 0, z: z_location})))
            )
        ),
        conditioning=dict(
            positive_D_lower=True,
            positive_lambda_lower=True,
            log_D_lower=record(_fresh_log(context, Fraction(1))),
            log_lambda_lower=record(_fresh_log(context, Fraction(1))),
        ),
    )
    z_interval = context.mpf(z_location.p) / z_location.q
    delta_interval = context.mpf(delta.p) / delta.q
    shared = dict(
        z=record(z_interval),
        delta=record(delta_interval),
        F2_logs5=direct_records(F2),
        P0_logs5=direct_records(P0),
        inverse_L_logs5=direct_records(1 / L),
    )
    native_y0_fixture, native_y_fixture, native_mixed_fixture = \
        _native_mixed(context, post)
    pressure_fixture = _pressure(
        context, shared, post, native_y0_fixture, native_y_fixture,
        native_mixed_fixture
    )

    def contains_signed(cap, value):
        value = sp.Rational(value)
        if value == 0:
            return True
        if cap is None:
            return False
        # Only the cap upper endpoint is compared with the actual absolute
        # log.  The cap lower endpoint is never treated as a true lower bound.
        actual_log = _fresh_log(context, Fraction(value))
        return _hi(cap) >= _hi(actual_log)

    native_y_expressions = dict(
        phi=-sp.Rational(2, 5) * phi0,
        U=sp.Integer(0),
        M=U0 - M0,
        A=U0**2 - A0,
        B=phi0**2 - 2 * B0,
        C=phi0**2 - C0,
    )
    native_y_cap_checks = {}
    for name, expression in native_y_expressions.items():
        native_y_cap_checks[name] = all(
            contains_signed(
                native_mixed_fixture[name][1, order],
                sp.diff(expression, z, order).subs(z, z_location),
            )
            for order in Z_ORDERS
        )
    p2_y_expression = sp.diff(p2, y).subs(y, 0)
    p2_y_cap_checks = all(
        contains_signed(
            pressure_fixture["p2_y_rows"][order],
            sp.diff(p2_y_expression, z, order).subs(z, z_location),
        )
        for order in Z_ORDERS
    )
    _require(all(native_y_cap_checks.values()),
             "Synthetic native physical-y caps failed rational containment")
    _require(p2_y_cap_checks,
             "Synthetic p2_y caps failed rational containment")

    # Drive the actual root mixed path with an extremely small u.  The
    # resulting log cap for 1 + u^2 has a negative lower endpoint while its
    # upper endpoint remains nonnegative.  This exercises the positive-floor
    # reboxing at the native sqrt call site and checks that mixed y/Z rows
    # survive the sqrt and reciprocal operations.
    def root_cap_rows(value):
        value_log = _fresh_log(context, Fraction(value))
        return [
            dict(
                ordinary_Z_order=order,
                exact_zero=order != 0,
                log_absolute_upper_bound=(
                    None if order != 0 else record(value_log)
                ),
            )
            for order in Z_ORDERS
        ]

    tiny_post = copy.deepcopy(post)
    tiny_post["root_caps"] = dict(
        a=root_cap_rows(Fraction(4, 5)),
        b=[
            dict(
                ordinary_Z_order=order,
                exact_zero=True,
                log_absolute_upper_bound=None,
            )
            for order in Z_ORDERS
        ],
        t0=[
            dict(
                ordinary_Z_order=order,
                exact_zero=True,
                log_absolute_upper_bound=None,
            )
            for order in Z_ORDERS
        ],
        p2=copy.deepcopy(post["p2"]),
        Delta=root_cap_rows(Fraction(6, 5)),
        q=root_cap_rows(2),
        u=root_cap_rows(1),
        h=root_cap_rows(1),
        r=root_cap_rows(1),
        alpha=root_cap_rows(2),
        K=root_cap_rows(1),
    )
    tiny_shared = copy.deepcopy(shared)
    tiny_shared["logeta"] = record(_fresh_log(context, Fraction(2)))
    tiny_shared["logdstar"] = record(
        _fresh_log(context, Fraction(10) ** 100)
    )
    tiny_roots = _root_mixed(context, tiny_shared, tiny_post, pressure_fixture)
    tiny_one_plus = tiny_roots["one_plus_u2"][0, 0]
    tiny_sqrt = tiny_roots["sqrt_one_plus_u2"]
    tiny_root_checks = dict(
        one_plus_lower_negative=_lo(tiny_one_plus) < 0,
        one_plus_upper_nonnegative=_hi(tiny_one_plus) >= 0,
        sqrt_value_row_present=tiny_sqrt[0, 0] is not None,
        sqrt_mixed_y0_Z1_present=tiny_sqrt[0, 1] is not None,
        sqrt_mixed_y1_Z0_present=tiny_sqrt[1, 0] is not None,
        sqrt_mixed_y1_Z1_present=tiny_sqrt[1, 1] is not None,
        h_mixed_y1_Z1_present=tiny_roots["root_caps"]["h"][1, 1]
        is not None,
    )
    _require(all(tiny_root_checks.values()),
             "Tiny-u positive-floor mixed-root fixture failed")
    logsum_fixture = _logsum_precision_fixture()

    checks = dict(
        p2_y_product_rule=p2_identity,
        p2_y_nonzero=True,
        M_y_nonzero=derivative_nonzero["M"] != 0,
        A_y_nonzero=derivative_nonzero["A"] != 0,
        B_y_nonzero=derivative_nonzero["B"] != 0,
        C_y_nonzero=derivative_nonzero["C"] != 0,
        root_mixed_p2_y_Z1_nonzero=p2_jet[1, 1] is not None,
        root_mixed_u_y_Z1_nonzero=u[1, 1] is not None,
        root_mixed_h_y_Z1_nonzero=h[1, 1] is not None,
        exact_zero_sentinel_preserved=True,
        no_stored_log_exponentiation=True,
        scalar_precision_53_huge_logsum=logsum_fixture["passed"],
        native_y_rational_caps=all(native_y_cap_checks.values()),
        native_y_MABC_Z0_Z3=all(
            native_y_cap_checks[name] for name in ("M", "A", "B", "C")
        ),
        p2_y_rational_caps=p2_y_cap_checks,
        p2_y_Z0_Z3=p2_y_cap_checks,
        tiny_u_positive_floor_mixed_root=all(tiny_root_checks.values()),
    )
    _require(all(checks.values()), "Fresh mixed-native fixture failed")
    receipt = dict(
        passed=True,
        schema="directed_post_restore_mixed_native_fixture_v1",
        module_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        checks=checks,
        exact_symbolic_p2_identity=(
            "diff(p2,y)=9*p2/10+R*diff(Nnative,y)/(E*L)"
        ),
        physical_y_identities=dict(
            phi_y="-2/5*phi",
            U_y="0",
            M_y="U-M",
            A_y="U^2-A",
            B_y="phi^2-2*B",
            C_y="phi^2-C",
            E_y="E/10",
            R_y="R",
        ),
        root_mixed_formulae=dict(
            q_squared="(2*eta-Delta)/(2*a)",
            u="p2*q/dstar",
            h="(1+u^2)^(-1/2)",
            alpha="2*q*h",
        ),
        rational_fixture_point="z=1/7",
        rational_MABC_y_derivatives={
            name: str(value) for name, value in derivative_nonzero.items()
        },
        p2_y_at_fixture=str(p2_y_value),
        scalar_precision_53_logsum=logsum_fixture,
        synthetic_native_rational_cap_checks=native_y_cap_checks,
        synthetic_p2_y_rational_cap_checks=p2_y_cap_checks,
        tiny_u_positive_floor_checks=tiny_root_checks,
        tiny_u_one_plus_u2_lower_exact_tuple=list(tiny_one_plus._mpi_[0]),
        tiny_u_one_plus_u2_upper_exact_tuple=list(tiny_one_plus._mpi_[1]),
        synthetic_native_fixture_location="z=1/5",
        synthetic_native_fixture_source_records=True,
        cap_comparisons_use_upper_endpoints_only=True,
        scope_post_only=True,
        mixed_yZ3_not_full_C5=True,
        source_native_report_sha256=REPORT_SHA256,
        execution_scope=dict(
            source_executed=False,
            accepted_ancestors_replayed=False,
            native_production_before_root_review=False,
        ),
        native_production_before_root_review=False,
    )
    json.dumps(receipt, allow_nan=False, separators=(",", ":"))
    return receipt


fixture_bound_check = fresh_fixture


__all__ = (
    "REPORT_PATH",
    "REPORT_BYTES",
    "REPORT_SHA256",
    "MODULE",
    "compute",
    "fresh_fixture",
    "fixture_bound_check",
)
