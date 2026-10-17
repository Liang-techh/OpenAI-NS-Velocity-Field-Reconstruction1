"""Same-source forward scale jets for the archived ordinary-Z n=1 prefix.

The n=1 coefficient archive is stored in the normalized fields
``f,u,k,p``.  This module binds its four physical multipliers to the accepted
core source.  The source amplitude is deliberately kept as a logarithm:
``F0/epsilon = exp(-logCstar-Lambda*G(Z)+logLambda)``.  The ordinary-Z
derivatives of that multiplier are nevertheless compact directed intervals;
their binary exponents can be of order ``exp(40)`` and are never converted to
decimal strings or selected as point values.

Only the archived centers ``-0.3``, ``0``, and ``0.3`` are exported.  The
saved log-amplitude packet supplies the absolute ``G`` enclosure at each
center.  Fresh interval jet arithmetic evaluates ``G'``, ``G''`` and ``G'''``
from the pinned rational source formula.  No owner, recurrence, quadrature,
or finite-N producer is executed here.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


PARENT = "lei_ren_part1_paper_compliant_current_original_n1_canonical_norm_tail.json"
PARENT_SHA256 = "b71ee21cf0f1dd2f3fe1d4edf4f45037ced997f1d72e2b8e0c9f1c861aae6989"
LOG_AMPLITUDE = "lei_ren_part1_paper_compliant_current_original_core_log_amplitude.json"
LOG_AMPLITUDE_SHA256 = "cf31caf9f6bf2bf480d241257b2ab6797252a2ab682cab46e5cc59ef6caa32cf"
TUBE = "lei_ren_part1_paper_shared_analytic_tube.json"
TUBE_SHA256 = "1e0de95d87056a3841c753877b51a07655d2d04420e8dc0e7a9d42ba59a608d6"
CORE_NUMERIC = "src/openai_ns_reconstruction/original_core_numeric.py"
CORE_NUMERIC_SHA256 = "9c152f5ff50bd760ee1a860a11e48c1b3a2fd36027681ea05a3e7a0920181603"
N1_NUMERIC = "src/openai_ns_reconstruction/original_n1_core_numeric.py"
N1_NUMERIC_SHA256 = "a7a925360f3a6da2348c5f18ee74cfbb885f8e220730769131e43a18bd954c2f"
CORE_FUNCTIONS = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py"
)
CORE_FUNCTIONS_SHA256 = "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59"
LOG_IMPLEMENTATION = "src/openai_ns_reconstruction/original_core_log_amplitude.py"
LOG_IMPLEMENTATION_SHA256 = "cfa8cb7566a4901f00901dc60f756acd5c8d4aa0995abcb8b181f5b11bd71924"
LOG_SCRIPT = "scripts/run_original_core_log_amplitude.py"
LOG_SCRIPT_SHA256 = "07fa1620e37b1249f0bd6b4dbd7fe4584bb37ccf36e862522ad7af059ad814e9"

CENTERS = ("-0.3", "0", "0.3")
ORDERS = (0, 1, 2, 3)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _read(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _digest(value: Any) -> str:
    """Digest mappings independently of whether JSON integers were strings."""

    def normalize(item: Any) -> Any:
        if item is None or isinstance(item, bool) or isinstance(item, str):
            return item
        if isinstance(item, int):
            return str(item)
        if isinstance(item, float):
            return repr(item)
        if isinstance(item, Mapping):
            return {str(k): normalize(item[k]) for k in sorted(item, key=str)}
        if isinstance(item, (list, tuple)):
            return [normalize(v) for v in item]
        return repr(item)

    payload = json.dumps(normalize(value), sort_keys=True, separators=(",", ":"),
                         ensure_ascii=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _endpoints(value: Any) -> tuple[mp.mpf, mp.mpf]:
    return mp.make_mpf(value._mpi_[0]), mp.make_mpf(value._mpi_[1])


def _interval(ctx: MPIntervalContext, value: Any) -> Any:
    if hasattr(value, "_mpi_"):
        return value
    if isinstance(value, Mapping):
        enclosure = value.get("source_enclosure_exact_mpf_tuples")
        if isinstance(enclosure, (list, tuple)) and len(enclosure) == 2:
            lo, hi = enclosure
            return ctx.mpf([mp.make_mpf(tuple(int(v) for v in lo)),
                            mp.make_mpf(tuple(int(v) for v in hi))])
        lo = value.get("lower_exact_mpf_tuple")
        hi = value.get("upper_exact_mpf_tuple")
        if lo is not None and hi is not None:
            return ctx.mpf([mp.make_mpf(tuple(int(v) for v in lo)),
                            mp.make_mpf(tuple(int(v) for v in hi))])
        if "lower" in value and "upper" in value:
            return ctx.mpf([str(value["lower"]), str(value["upper"])])
        if "upper" in value:
            return ctx.mpf(str(value["upper"]))
    return ctx.mpf(str(value))


def _record(value: Any) -> dict[str, Any]:
    lo, hi = _endpoints(value)
    return {
        "lower_exact_mpf_tuple": list(map(int, lo._mpf_)),
        "upper_exact_mpf_tuple": list(map(int, hi._mpf_)),
        "lower_decimal_display": mp.nstr(lo, 28),
        "upper_decimal_display": mp.nstr(hi, 28),
        "directed_interval": True,
    }


def _roundtrip(ctx: MPIntervalContext, value: Any) -> int:
    if isinstance(value, Mapping):
        if value.get("directed_interval") is True:
            restored = _interval(ctx, value)
            expected = tuple(tuple(int(v) for v in value[key])
                             for key in ("lower_exact_mpf_tuple", "upper_exact_mpf_tuple"))
            _require(restored._mpi_ == expected, "Scale interval did not round trip")
            return 1
        return sum(_roundtrip(ctx, v) for v in value.values())
    if isinstance(value, (list, tuple)):
        return sum(_roundtrip(ctx, v) for v in value)
    return 0


def _sign(value: Any) -> int:
    lo, hi = _endpoints(value)
    if lo > 0:
        return 1
    if hi < 0:
        return -1
    if lo == 0 == hi:
        return 0
    raise ValueError("A source derivative interval crosses zero")


def _constant(ctx: MPIntervalContext, value: Any) -> list[Any]:
    return [value, ctx.mpf(0), ctx.mpf(0), ctx.mpf(0)]


def _add(left: list[Any], right: list[Any]) -> list[Any]:
    return [left[i] + right[i] for i in ORDERS]


def _scale(factor: Any, value: list[Any]) -> list[Any]:
    return [factor * item for item in value]


def _mul(left: list[Any], right: list[Any]) -> list[Any]:
    return [sum((math.comb(i, k) * left[k] * right[i - k]
                 for k in range(i + 1)), left[0] * 0)
            for i in ORDERS]


def _inverse(value: list[Any]) -> list[Any]:
    _require(_sign(value[0]) != 0, "Jet denominator contains zero")
    result = [1 / value[0], value[0] * 0, value[0] * 0, value[0] * 0]
    for i in (1, 2, 3):
        result[i] = -sum((math.comb(i, k) * value[k] * result[i - k]
                          for k in range(1, i + 1)), value[0] * 0) / value[0]
    return result


def _divide(left: list[Any], right: list[Any]) -> list[Any]:
    return _mul(left, _inverse(right))


def _derivative_jet(ctx: MPIntervalContext, center: str, j: Any, delta: Any) -> list[Any]:
    """Return ``[G',G'',G''',G'''']`` as ordinary derivative intervals."""

    zero = ctx.mpf(0)
    one = ctx.mpf(1)
    z = [ctx.mpf(center), one, zero, zero]
    z2 = _mul(z, z)
    de = _constant(ctx, delta)
    L = _add(_constant(ctx, one), _scale(-1, _mul(de, z2)))
    U0 = _add(_scale(4, z), _constant(ctx, j))
    H = _add(_scale(ctx.mpf("0.5"), _mul(z, _add(_constant(ctx, one), _scale(-1, de)))),
             _mul(_add(_constant(ctx, one), _scale(-1, z2)), U0))
    sigma = _scale(ctx.mpf(1) / 500, _constant(ctx, j))
    denominator = _add(_mul(H, H), _mul(sigma, sigma))
    gradient = _divide(_mul(L, H), denominator)
    _require(_sign(denominator[0]) == 1, "G source denominator is not positive")
    return gradient


def _source_rows(root: Path, parent: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    report_path = root / "experiments/root_st073" / PARENT
    _require(_sha(report_path) == PARENT_SHA256, "Canonical norm/tail parent hash changed")
    _require(_digest(parent) == _digest(_read(report_path)),
             "Supplied canonical norm/tail mapping differs from saved parent")
    log_path = root / "experiments/root_st073" / LOG_AMPLITUDE
    tube_path = root / "experiments/root_st073" / TUBE
    _require(_sha(log_path) == LOG_AMPLITUDE_SHA256, "Saved log-amplitude source changed")
    _require(_sha(tube_path) == TUBE_SHA256, "Saved analytic-tube source changed")
    log_report = _read(log_path)
    tube = _read(tube_path)
    family = parent.get("source_family")
    _require(log_report.get("source_family") == family, "Log-amplitude source family differs")
    _require(log_report.get("checks", {}).get("passed") is True,
             "Accepted log-amplitude report is not passed")
    _require(log_report["checks"].get("source_midpoint_selection") is False and
             log_report["checks"].get("absolute_F0_materialized") is False,
             "Log-amplitude source selected a midpoint or materialized F0")
    source_identity = parent["n1_source_vector_identity"]
    _require(source_identity.get("passed") is True, "Accepted source-vector identity required")
    pins = source_identity["source_pins"]
    for key, expected in (("forward_common_core_scales", CORE_NUMERIC_SHA256),
                          ("forward_Pstar_squared_scale", N1_NUMERIC_SHA256)):
        pin = pins[key]
        path = root / pin["path"]
        _require(pin["sha256"] == expected and _sha(path) == expected,
                 "Core scale source pin changed: " + key)
    functions_path = root / CORE_FUNCTIONS
    _require(_sha(functions_path) == CORE_FUNCTIONS_SHA256,
             "Canonical G/F0 source changed")
    for rel, expected in ((LOG_IMPLEMENTATION, LOG_IMPLEMENTATION_SHA256),
                          (LOG_SCRIPT, LOG_SCRIPT_SHA256)):
        _require(_sha(root / rel) == expected, "Log-amplitude implementation changed: " + rel)
        _require(log_report.get("implementation_hashes", {}).get(rel) == expected,
                 "Log-amplitude report does not bind implementation: " + rel)
    for name, expected in ((TUBE, TUBE_SHA256),):
        saved = log_report.get("input_hashes", {}).get(name)
        _require(saved == expected, "Log-amplitude report does not bind the analytic tube")
    delta_admission = parent["n1_complex_actual_source_bounds"]["delta_admission"]
    tube_delta = tube.get("delta")
    admitted_delta = delta_admission.get("tube_delta_record")
    _require(isinstance(tube_delta, Mapping) and isinstance(admitted_delta, Mapping),
             "Analytic-tube directed delta record is missing")
    for bound in ("lower_exact_mpf_tuple", "upper_exact_mpf_tuple"):
        _require(tuple(str(v) for v in tube_delta.get(bound, ())) ==
                 tuple(str(v) for v in admitted_delta.get(bound, ())),
                 "Analytic-tube delta record differs from admitted source record")
    _require(tube_delta.get("upper") in ("1e-200", "1.0e-200"),
             "Analytic-tube delta cap changed")
    _require(delta_admission.get("passed") is True and
             delta_admission.get("actual_delta_positive") is True and
             delta_admission.get("actual_delta_le_tube_cap") is True and
             delta_admission.get("actual_delta_exponential_branch_proved") is True,
             "Actual delta admission is not proved")
    return log_report, tube, delta_admission


def _global_rows(log_report: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    levels = [level for level in log_report.get("levels", [])
              if isinstance(level, Mapping) and isinstance(level.get("global_amplitude"), list)]
    _require(levels, "No saved global amplitude rows")
    level = max(levels, key=lambda item: int(item.get("precision", 0)))
    rows = {str(row.get("Z")): row for row in level["global_amplitude"]}
    _require(all(center in rows for center in CENTERS), "Missing archived scale center")
    _require(len(rows) == len(level["global_amplitude"]), "Duplicate saved scale center")
    for center in CENTERS:
        row = rows[center]
        _require("G" in row and "source_enclosure_exact_mpf_tuples" in row["G"],
                 "Saved G row has no directed source enclosure")
        _require(row.get("positive_F0", {}).get("kind") ==
                 "positive_exponential_of_separate_terms",
                 "Saved F0 source is not positive separated exponential")
        _require(row.get("positive_F0", {}).get("relative_log_term") ==
                 row.get("log_F0_relative_to_anchor"),
                 "F0 relative log source mismatch")
    return rows


def _signed_log_record(ctx: MPIntervalContext, value: Any) -> dict[str, Any]:
    sign = _sign(value)
    if sign == 0:
        return {"kind": "exact_zero", "sign": 0, "value": _record(value),
                "ordinary_derivative_ratio": True}
    magnitude = value if sign > 0 else -value
    return {"kind": "signed_log_interval", "sign": sign,
            "log_abs": _record(ctx.log(magnitude)),
            "ordinary_derivative_ratio": True,
            "materialized_interval": _record(value)}


def _center_packet(ctx: MPIntervalContext, center: str, row: Mapping[str, Any],
                   j: Any, delta: Any, log_p: Any, log_lambda: Any,
                   log_c: Any, epsilon: Any, epsilon_p2: Any,
                   lambda_value: Any) -> dict[str, Any]:
    g0 = _interval(ctx, row["G"])
    derivative_jet = _derivative_jet(ctx, center, j, delta)
    g = [g0, derivative_jet[0], derivative_jet[1], derivative_jet[2]]
    for value in g[1:]:
        _require(_sign(value) != 0, "G derivative source interval crosses zero")
    log_f0_relative = _interval(ctx, row["log_F0_relative_to_anchor"])
    log_f0 = -log_c + log_f0_relative
    log_f = log_f0 + log_lambda
    # Ordinary derivative ratios of exp(-Lambda G):
    # r1=a1, r2=a2+a1^2, r3=a3+3*a1*a2+a1^3, a_k=-Lambda G^(k).
    a1, a2, a3 = (-lambda_value * g[1], -lambda_value * g[2], -lambda_value * g[3])
    relative = [ctx.mpf(1), a1, a2 + a1 * a1,
                a3 + 3 * a1 * a2 + a1 * a1 * a1]
    relative_log_jets = [_signed_log_record(ctx, value) for value in relative]
    one = _record(ctx.mpf(1))
    zero = _record(ctx.mpf(0))

    def entry(name: str, log_scale: Any, relative_values: list[Any], formula: str,
              materialized: Any | None, safe: bool) -> dict[str, Any]:
        return {
            "field_scale": name,
            "formula": formula,
            "log_positive_scale": _record(log_scale),
            "relative_ordinary_Z_jets": [_record(value) for value in relative_values],
            "relative_derivative_Z_jets": [_record(value) for value in relative_values],
            "relative_ordinary_Z_jets_log": relative_log_jets if name == "F0_over_epsilon"
            else [{"kind": "log_abs_zero", "sign": 1 if i == 0 else 0,
                   "log_abs": _record(ctx.mpf(0)) if i == 0 else zero,
                   "ordinary_derivative_ratio": True} for i in ORDERS],
            "relative_jet_convention": "D_Z^j(scale)/scale; no factorial divisor",
            "ordinary_derivative_convention": "D_Z^j S / S, no factorial divisor",
            "materialized_scale": _record(materialized) if materialized is not None else None,
            "exponential_materialization_safe": safe,
            "inverse_log_positive_scale": _record(-log_scale),
        }

    f_entry = entry("F0_over_epsilon", log_f, relative,
                    "exp(-selected_logCstar-Lambda*G(Z)+logLambda)", None, False)
    const_rel = [ctx.mpf(1), ctx.mpf(0), ctx.mpf(0), ctx.mpf(0)]
    eps_entry = entry("epsilon", -log_lambda, const_rel, "exp(-logLambda)", epsilon, True)
    p_entry = entry("epsilon_Pstar_squared", -log_lambda + 2 * log_p, const_rel,
                    "exp(-logLambda+2*logPstar)", epsilon_p2, True)
    # The k scale is exactly epsilon and is intentionally a separate source
    # field even though the numerical multiplier is the same.
    k_entry = entry("epsilon", -log_lambda, const_rel, "exp(-logLambda)", epsilon, True)
    return {
        "Z": center,
        "G_derivatives": [_record(value) for value in g],
        "G_derivative_source": {
            "formula": "Gprime=L*H/(H^2+(j/500)^2)",
            "L": "1-delta*Z^2",
            "H": "Z*(1-delta)/2+(1-Z^2)*(4*Z+j)",
            "delta_interval_is_actual_cap": True,
            "ordinary_derivative_orders": [1, 2, 3],
        },
        "log_F0": _record(log_f0),
        "log_F0_relative_to_anchor": _record(log_f0_relative),
        "raw_scales": {"f": f_entry, "u": eps_entry, "k": k_entry, "p": p_entry},
        "scale_aliases": {
            "F1": "f", "F1x_over_x": "f", "Uz1": "u", "Uz1x_over_x": "u",
            "K1": "k", "P1": "p",
        },
        "normalization": {
            "F1": "(F0/epsilon)*f", "Uz1": "epsilon*u", "K1": "epsilon*k",
            "P1": "(epsilon*Pstar^2)*p",
        },
        "zero_order_relative_jet": {"ordinary_derivative_ratio": True,
                                     "value": one},
    }


def _build(root: Path, parent: Mapping[str, Any], centers: tuple[str, ...],
           precision: int) -> dict[str, Any]:
    _require(type(precision) is int and precision >= 150, "At least150 interval digits required")
    log_report, tube, delta_admission = _source_rows(root, parent)
    rows = _global_rows(log_report)
    prefix = parent["n1_actual_raw_prefix_binding"]["actual_forward_prefix_bound"]
    _require(prefix.get("passed") is True and
             prefix["ordinary_Z_domain"]["coefficient_centers"] == list(CENTERS),
             "Forward prefix centers differ from the saved scale centers")
    _require(prefix["physical_scale_handles"]["F0"] ==
             "F0(Z)=exp(-selected_logCstar-Lambda*G(Z))",
             "Forward prefix F0 scale handle changed")
    ctx = MPIntervalContext()
    ctx.dps = precision
    log_p = ctx.exp(ctx.mpf(40)) + 11
    log_lambda = 4 * log_p + 1000
    delta_record = _interval(ctx, delta_admission["tube_delta_record"])
    delta_upper = _endpoints(delta_record)[1]
    delta = ctx.mpf([mp.make_mpf(ctx.mpf(0)._mpi_[0]),
                     mp.make_mpf(delta_upper._mpf_)])
    j = _interval(ctx, tube["j"])
    epsilon = ctx.exp(-log_lambda)
    epsilon_p2 = ctx.exp(-2 * log_p - 1000)
    # Lambda is a compact MPF intermediate for the relative derivative jets;
    # its decimal expansion and every physical F0 amplitude remain unmaterialized.
    lambda_value = ctx.exp(log_lambda)
    first_row = rows[CENTERS[0]]
    anchor_tuple = first_row["positive_F0"]["exact_anchor_log_term"]["negative_of_exact_mpf_tuple"]
    log_c = ctx.mpf([mp.make_mpf(tuple(int(v) for v in anchor_tuple))] * 2)
    center_packets = {
        center: _center_packet(ctx, center, rows[center], j, delta, log_p, log_lambda,
                               log_c, epsilon, epsilon_p2, lambda_value)
        for center in centers
    }
    packet = {
        "passed": True,
        "parent_report": PARENT,
        "parent_report_sha256": PARENT_SHA256,
        "source_report": LOG_AMPLITUDE,
        "source_report_sha256": LOG_AMPLITUDE_SHA256,
        "analytic_tube_report": TUBE,
        "analytic_tube_sha256": TUBE_SHA256,
        "source_family": parent["source_family"],
        "precision_decimal_digits": precision,
        "centers": center_packets,
        "checks": {
            "parent_report_hash_and_mapping": True,
            "source_family_match": True,
            "source_log_amplitude_hash": True,
            "source_log_amplitude_passed": True,
            "source_no_midpoint_or_F0_materialization": True,
            "core_scale_source_ast_pins": True,
            "canonical_G_source_hash": True,
        "analytic_tube_hash_and_delta_cap": True,
            "actual_delta_positive_and_below_cap": True,
            "forward_prefix_center_and_scale_binding": True,
            "saved_G_center_enclosures": True,
            "fresh_G_derivative_jet_arithmetic": True,
            "ordinary_derivative_ratio_convention": True,
            "physical_large_scales_factored": True,
            "no_owner_builder_fixture_recurrence_or_quadrature_replayed": True,
        },
        "source_pins": {
            "core_scale_source": {"path": CORE_NUMERIC, "sha256": CORE_NUMERIC_SHA256},
            "n1_Pstar_squared_source": {"path": N1_NUMERIC, "sha256": N1_NUMERIC_SHA256},
            "canonical_G_F0_source": {"path": CORE_FUNCTIONS, "sha256": CORE_FUNCTIONS_SHA256},
            "log_amplitude_implementation": {"path": LOG_IMPLEMENTATION,
                                               "sha256": LOG_IMPLEMENTATION_SHA256},
            "log_amplitude_source_script": {"path": LOG_SCRIPT, "sha256": LOG_SCRIPT_SHA256},
        },
        "formulae": {
            "logPstar": "exp(40)+11",
            "logLambda": "4*logPstar+1000",
            "epsilon": "exp(-logLambda)",
            "Pstar_squared": "exp(2*logPstar)",
            "F0": "exp(-selected_logCstar-Lambda*G(Z))",
            "F0_over_epsilon": "exp(-selected_logCstar-Lambda*G(Z)+logLambda)",
            "epsilon_Pstar_squared": "exp(-2*logPstar-1000)",
            "relative_derivative_Z_jets": "D_Z^j S/S; ordinary derivatives, no factorial divisor",
            "G_derivative_source": "Gprime=(1-delta*Z^2)*H/(H^2+(j/500)^2)",
        },
        "global_scale_records": {
            "logPstar": _record(log_p),
            "logLambda": _record(log_lambda),
            "epsilon": _record(epsilon),
            "epsilon_Pstar_squared": _record(epsilon_p2),
            "logCstar": _record(log_c),
            "physical_Lambda_materialized": False,
            "physical_Pstar_materialized": False,
            "compact_lambda_intermediate_used_for_relative_jets": True,
        },
        "delta_admission": delta_admission,
        "limitations": [
            "F0 and F0/epsilon are retained as positive log-domain amplitudes; no decimal physical F0 is selected",
            "Only the three exact ordinary-Z archive centers are exported",
            "This packet supplies scale jets and does not certify a physical NS total or time recursion",
        ],
    }
    packet["exact_interval_records_roundtripped"] = _roundtrip(ctx, packet)
    packet["checks"]["exact_interval_records_roundtripped"] = packet["exact_interval_records_roundtripped"] > 0
    _require(packet["checks"]["exact_interval_records_roundtripped"],
             "Scale packet contains no round-tripped directed intervals")
    return packet


def build(root: str | Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    """Build the accepted same-source scale-jet packet."""

    return _build(Path(root), parent, CENTERS, 180)


def evaluate_scale_jets(root: str | Path, parent: Mapping[str, Any],
                        centers: tuple[str, ...] = CENTERS,
                        precision: int = 180) -> dict[str, Any]:
    """Callable helper returning the same-source center scale records."""

    _require(tuple(centers), "At least one archived center is required")
    _require(all(center in CENTERS for center in centers),
             "Scale helper accepts archived centers only")
    packet = _build(Path(root), parent, tuple(centers), precision)
    return packet
