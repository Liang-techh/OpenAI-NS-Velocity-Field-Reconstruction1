"""Finite-``N`` propagation from the saved directed loop/native reports.

The adapter is intentionally a data boundary.  It authenticates the two
complete saved reports before decoding JSON, reads only their serialized
directed rows, and sends logarithmic caps to
``original_directed_finite_N_kernel``.  No historical finite-``N`` producer,
owner, checker, recurrence, quadrature, or norm is imported or executed.

The four charts are propagated in the order ``long_reshape``, ``reference``,
``restoration``, ``post_restore``.  The first chart receives the saved
post-switch ``delta_M`` cap.  Each later chart receives the previous chart's
outgoing cap; there is no memory reset.  Signed centres are copied from the
saved background rows, and a symmetric ``2**-1024`` radius is installed only
when the new local error target passes.
"""

from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any

from mpmath import isfinite, make_mpf, nstr
from mpmath.ctx_iv import MPIntervalContext


LOOP_REPORT_NAME = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_directed_loop_controls.json"
)
LOOP_REPORT_SHA256 = (
    "adefdae708754a2572eae15a90662274bf17e88a0d72d5db9befbc9e2858b252"
)
LOOP_REPORT_BYTES = 1752702

NATIVE_REPORT_NAME = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json"
)
NATIVE_REPORT_SHA256 = (
    "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"
)
NATIVE_REPORT_BYTES = 42177499
NATIVE_REPORT_BASENAME = NATIVE_REPORT_NAME.rsplit("/", 1)[-1]

MODULE = "src/openai_ns_reconstruction/original_directed_finite_N_controls.py"
KERNEL_MODULE = "src/openai_ns_reconstruction/original_directed_finite_N_kernel.py"
SCHEMA = "original_directed_finite_N_controls_v1"

CHARTS = ("long_reshape", "reference", "restoration", "post_restore")
LONG_REFERENCE = ("long_reshape", "reference")
RESTORATION = ("restoration", "post_restore")
ORDERS = (0, 1, 2, 3, 4)
I1_ORDERS = (0, 1)
Q_ORDERS = (0, 1, 2)
I1_NAMES = ("q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z")
COEFFICIENT_NAMES = (
    "minus_inv_L",
    "hZ_over_2L",
    "minus_half",
    "d_over_L",
    "two_h_Z2_over_L2",
    "two_Zd_over_L2",
)
COEFFICIENT_FOR_I1 = dict(zip(I1_NAMES, COEFFICIENT_NAMES))
NORMALIZATIONS = {
    "long_reshape": "Wlong",
    "reference": "Wref",
    "restoration": "Wrestore",
    "post_restore": "Wpost",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _canonical_endpoint(value: Any, label: str) -> tuple[int, int, int, int]:
    _require(isinstance(value, (list, tuple)) and len(value) == 4,
             label + " must be an MPF tuple")
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
        _require(mantissa & 1,
                 label + " MPF mantissa is not normalized")
    return sign, mantissa, exponent, bit_count


def _record(c: MPIntervalContext, value: Any, label: str = "interval") -> dict:
    """Serialize a finite directed interval without selecting a midpoint."""

    _require(hasattr(value, "_mpi_"), label + " must be an interval")
    lower = tuple(value._mpi_[0])
    upper = tuple(value._mpi_[1])
    _canonical_endpoint(lower, label + ".lower")
    _canonical_endpoint(upper, label + ".upper")
    _require(make_mpf(lower) <= make_mpf(upper),
             label + " has reversed endpoints")
    _require(isfinite(make_mpf(lower)) and isfinite(make_mpf(upper)),
             label + " must be finite")
    return {
        "lower_exact_mpf_tuple": list(lower),
        "upper_exact_mpf_tuple": list(upper),
        "lower_display": nstr(make_mpf(lower), 24),
        "upper_display": nstr(make_mpf(upper), 24),
        "directed": True,
    }


def _read_interval(c: MPIntervalContext, value: Any, label: str) -> Any:
    """Decode one exact serialized interval, rejecting scalar value rows."""

    if hasattr(value, "_mpi_"):
        lower = tuple(value._mpi_[0])
        upper = tuple(value._mpi_[1])
    else:
        _require(isinstance(value, Mapping), label + " must be a directed record")
        _require(value.get("directed") is True,
                 label + " must set directed=True")
        lower = _canonical_endpoint(value.get("lower_exact_mpf_tuple"),
                                    label + ".lower")
        upper = _canonical_endpoint(value.get("upper_exact_mpf_tuple"),
                                    label + ".upper")
    lo = make_mpf(lower)
    hi = make_mpf(upper)
    _require(lo <= hi, label + " has reversed endpoints")
    _require(isfinite(lo) and isfinite(hi), label + " must be finite")
    return c.mpf([lo, hi])


def _admit_bytes(raw: bytes, expected_bytes: int, expected_sha256: str,
                 label: str) -> Mapping[str, Any]:
    """Authenticate complete report bytes before JSON decoding."""

    _require(type(raw) is bytes, label + " must be bytes")
    _require(len(raw) == expected_bytes,
             label + " byte length changed")
    _require(hashlib.sha256(raw).hexdigest() == expected_sha256,
             label + " SHA256 changed")
    try:
        value = json.loads(raw)
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError(label + " is not valid JSON") from exc
    _require(isinstance(value, Mapping), label + " must be a full report mapping")
    return value


def _admit_reports(loop_report_bytes: bytes,
                   native_report_bytes: bytes) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    """Admit only the exact full reports, never a nested packet alone."""

    loop = _admit_bytes(loop_report_bytes, LOOP_REPORT_BYTES,
                        LOOP_REPORT_SHA256, "Loop report")
    native = _admit_bytes(native_report_bytes, NATIVE_REPORT_BYTES,
                          NATIVE_REPORT_SHA256, "Native report")
    _require(loop.get("stage") == "current_original_directed_loop_controls",
             "Loop input is not the reviewed full report")
    _require(loop.get("parent_report") == NATIVE_REPORT_BASENAME
             and loop.get("parent_report_sha256") == NATIVE_REPORT_SHA256,
             "Loop report parent binding changed")
    _require(native.get("stage") == "current_original_native_directed_controls_v2",
             "Native input is not the reviewed full report")
    packet = loop.get("actual_directed_loop_controls")
    _require(isinstance(packet, Mapping)
             and packet.get("passed") is True
             and packet.get("native_parent_report_sha256") == NATIVE_REPORT_SHA256
             and packet.get("native_parent_report_bytes") == NATIVE_REPORT_BYTES,
             "Loop packet binding changed")
    native_packet = native.get("actual_native_directed_controls_v2")
    _require(isinstance(native_packet, Mapping)
             and native_packet.get("schema") == "original_native_directed_controls_v2_v1"
             and native_packet.get("passed") is True,
             "Corrected native packet is missing")
    scope = native_packet.get("scope")
    _require(isinstance(scope, Mapping)
             and scope.get("corrected_native_controls_verified") is True
             and scope.get("inherited_loop_revalidated") is False
             and scope.get("current_N_revalidated") is False
             and scope.get("inherited_log_caps_revalidated") is False,
             "Native upstream scope was promoted or changed")
    _require(isinstance(native_packet.get("raw_interval_tree"), Mapping),
             "Native raw interval tree is missing")
    return loop, native


def _ordered_rows(rows: Any, count: int, label: str) -> list[Any]:
    _require(isinstance(rows, Sequence) and not isinstance(rows, (str, bytes)),
             label + " must be a row sequence")
    _require(len(rows) == count,
             label + " must contain exactly " + str(count) + " rows")
    out = []
    for index, row in enumerate(rows):
        if isinstance(row, Mapping) and "ordinary_Z_order" in row:
            _require(row.get("ordinary_Z_order") == index,
                     label + " rows are not ordered")
        out.append(row)
    return out


def _raw_row(row: Any, *, logarithmic: bool, label: str) -> Any:
    if row is None:
        return None
    if isinstance(row, Mapping):
        if row.get("exact_zero") is True:
            if logarithmic:
                _require(row.get("log_absolute_upper_bound") is None,
                         label + " exact-zero row carries a bound")
            return None
        if logarithmic:
            for key in ("log_absolute_upper_bound", "log_cap", "value"):
                if key in row:
                    return row[key]
        else:
            for key in ("value", "signed_interval", "raw"):
                if key in row:
                    return row[key]
    return row


def _read_log_rows(c: MPIntervalContext, rows: Any, count: int,
                   label: str) -> list[Any]:
    values = []
    for index, row in enumerate(_ordered_rows(rows, count, label)):
        raw = _raw_row(row, logarithmic=True, label=label + "[" + str(index) + "]")
        values.append(None if raw is None else
                      _read_interval(c, raw, label + "[" + str(index) + "]"))
    return values


def _read_signed_rows(c: MPIntervalContext, rows: Any, count: int,
                      label: str) -> list[Any]:
    values = []
    for index, row in enumerate(_ordered_rows(rows, count, label)):
        raw = _raw_row(row, logarithmic=False, label=label + "[" + str(index) + "]")
        _require(raw is not None, label + " cannot contain an omitted signed row")
        values.append(_read_interval(c, raw, label + "[" + str(index) + "]"))
    return values


def _pack_log_rows(c: MPIntervalContext, rows: Sequence[Any], label: str) -> list[dict]:
    _require(len(rows) in (2, 3, 4, 5), label + " has an unsupported row count")
    out = []
    for index, value in enumerate(rows):
        out.append({
            "ordinary_Z_order": index,
            "exact_zero": value is None,
            "log_absolute_upper_bound": (
                None if value is None else _record(c, value, label + "[" + str(index) + "]")
            ),
        })
    return out


def _pack_signed_rows(c: MPIntervalContext, rows: Sequence[Any], label: str) -> list[dict]:
    _require(len(rows) == 2, label + " must contain two rows")
    return [
        {"ordinary_Z_order": index,
         "value": _record(c, value, label + "[" + str(index) + "]")}
        for index, value in enumerate(rows)
    ]


def _path(root: Mapping[str, Any], dotted: str, label: str | None = None) -> Any:
    value: Any = root
    for key in dotted.split("."):
        _require(isinstance(value, Mapping) and key in value,
                 (label or dotted) + " is missing")
        value = value[key]
    return value


def _max_upper_point(c: MPIntervalContext, left: Any, right: Any) -> Any:
    """Return the exact upper-endpoint point used for positive memory caps."""

    values = [value for value in (left, right) if value is not None]
    if not values:
        return None
    upper = max(value.b for value in values)
    return c.mpf([upper, upper])


def _symmetric(c: MPIntervalContext, radius: Any) -> Any:
    return c.mpf([-radius.b, radius.b])


def _log_add(kernel: Any, c: MPIntervalContext, left: Any, right: Any) -> Any:
    values = [value for value in (left, right) if value is not None]
    return None if not values else kernel.logsum(c, values)


def _strictly_below(c: MPIntervalContext, value: Any, threshold: Any) -> bool:
    return value is None or value.b < threshold.a


def _source_packet(native: Mapping[str, Any], chart: str) -> Mapping[str, Any]:
    packet = _path(native, "actual_native_directed_controls_v2")
    tree = _path(packet, "raw_interval_tree")
    family = "long_reference" if chart in LONG_REFERENCE else "restoration"
    section = _path(tree, family)
    return _path(section, "charts." + chart)


def _baseline_chart(native: Mapping[str, Any], chart: str) -> tuple[Mapping[str, Any], str]:
    if chart in LONG_REFERENCE:
        source = _path(native, "actual_long_reference_background_I1_numeric")
        return _path(source, "charts." + chart), (
            "nativeparent.actual_long_reference_background_I1_numeric.charts." + chart
        )
    repair = _path(native, "actual_restoration_post_jet_repair")
    corrected = repair.get("corrected_background")
    _require(isinstance(corrected, Mapping),
             "Corrected restoration background packet is missing")
    _require(corrected.get("corrected_numeric_view") is True
             and corrected.get("source_executed") is False,
             "Corrected restoration background identity changed")
    return _path(corrected, "charts." + chart), (
        "nativeparent.actual_restoration_post_jet_repair.corrected_background.charts."
        + chart
    )


def _baseline_data(c: MPIntervalContext, chart: Mapping[str, Any],
                   name: str) -> dict[str, Any]:
    arithmetic = chart.get("arithmetic")
    _require(isinstance(arithmetic, Mapping), name + ".arithmetic is missing")
    six_source = arithmetic.get("six_integrals")
    coefficients_source = arithmetic.get("coefficients")
    _require(isinstance(six_source, Mapping), name + " six-integral source is missing")
    _require(isinstance(coefficients_source, Mapping),
             name + " coefficient source is missing")
    six = {
        key: _read_signed_rows(c, six_source.get(key), 2,
                               name + ".six_integrals." + key)
        for key in I1_NAMES
    }
    coefficients = {
        key: _read_signed_rows(c, coefficients_source.get(key), 2,
                               name + ".coefficients." + key)
        for key in COEFFICIENT_NAMES
    }
    B = _read_signed_rows(c, arithmetic.get("B_over_W"), 2,
                          name + ".B_over_W")
    return {
        "six": six,
        "coefficients": coefficients,
        "B": B,
        "source": {
            "six_integrals": copy.deepcopy(six_source),
            "coefficients": copy.deepcopy(coefficients_source),
            "B_over_W": copy.deepcopy(arithmetic.get("B_over_W")),
            "D": copy.deepcopy(arithmetic.get("D")),
            "S1_over_W": copy.deepcopy(arithmetic.get("S1_over_W")),
        },
    }


def _full_signed_rows(c: MPIntervalContext, base: dict[str, Any], radius: Any,
                      normalization: str, source_path: str) -> dict[str, Any]:
    """Build signed full rows after a passing target, adding radius once."""

    symmetric = _symmetric(c, radius)
    six_full = {
        key: [value + symmetric for value in rows]
        for key, rows in base["six"].items()
    }
    B_full = [value + symmetric for value in base["B"]]
    products: dict[str, list[Any]] = {}
    for key in I1_NAMES:
        coeff = base["coefficients"][COEFFICIENT_FOR_I1[key]]
        values = six_full[key]
        products[key] = [
            coeff[0] * values[0],
            coeff[1] * values[0] + coeff[0] * values[1],
        ]
    D_full = [
        sum((products[key][order] for key in I1_NAMES),
            products[I1_NAMES[0]][order] * 0)
        for order in I1_ORDERS
    ]
    S1_full = [B_full[order] + D_full[order] for order in I1_ORDERS]
    coeff_rows = {
        key: _pack_signed_rows(c, values, "coefficient." + key)
        for key, values in base["coefficients"].items()
    }
    return {
        "I1_full": {
            "passed": True,
            "normalized_by": normalization,
            "six_integrals": {
                key: _pack_signed_rows(c, values, "I1_full." + key)
                for key, values in six_full.items()
            },
            "signed_background_unchanged": True,
            "symmetric_radius_added_once": _record(c, radius, "radius"),
            "error_radius_is_not_a_selected_value": True,
            "native_centre_source": source_path + ".arithmetic.six_integrals",
            "integrated_Q_used": False,
        },
        "B_full": {
            "passed": True,
            "normalized_by": normalization,
            "rows": _pack_signed_rows(c, B_full, "B_full"),
            "symmetric_radius_added_once": _record(c, radius, "radius"),
            "error_radius_is_not_a_selected_value": True,
            "signed_background_unchanged": True,
            "native_centre_source": source_path + ".arithmetic.B_over_W",
            "own_W_formula": "chart-specific normalized W; no additional W factor",
            "integrated_Q_used": False,
        },
        "D_full": {
            "passed": True,
            "normalized_by": normalization,
            "rows": _pack_signed_rows(c, D_full, "D_full"),
            "D_rows": _pack_signed_rows(c, D_full, "D_full"),
            "coefficient_rows": coeff_rows,
            "formula": "D=sum(saved coefficient C1 rows * full six-I1 C1 rows)",
            "computed_from_full_I1_once": True,
            "separate_D_error_not_added": True,
            "no_extra_D_radius": True,
        },
        "S1_full": {
            "passed": True,
            "normalized_by": normalization,
            "rows": _pack_signed_rows(c, S1_full, "S1_full"),
            "formula": "S1_full=B_full+D_full",
            "computed_from_B_and_D_once": True,
            "separate_S1_radius_not_added": True,
            "no_extra_radius": True,
        },
    }


def _baseline_receipt(base: dict[str, Any], source_path: str) -> dict[str, Any]:
    return {
        "source_path": source_path,
        "signed_centres_retained_exactly": True,
        "six_integrals": copy.deepcopy(base["source"]["six_integrals"]),
        "coefficients": copy.deepcopy(base["source"]["coefficients"]),
        "B_over_W": copy.deepcopy(base["source"]["B_over_W"]),
        "D": copy.deepcopy(base["source"]["D"]),
        "S1_over_W": copy.deepcopy(base["source"]["S1_over_W"]),
        "corrected_postrepair_source_used": "corrected_background" in source_path,
    }


def _chart(c: MPIntervalContext, kernel: Any, chart: str,
           loop_chart: Mapping[str, Any], native_chart: Mapping[str, Any],
           baseline_chart: Mapping[str, Any], source_path: str,
           incoming: Sequence[Any], z: Any, delta: Any, inverseL: Sequence[Any],
           logN: Any, mass_log: Any, radius: Any,
           shared_logN_path: str) -> tuple[dict[str, Any], list[Any]]:
    normalization = NORMALIZATIONS[chart]
    _require(isinstance(loop_chart, Mapping), chart + " loop chart is missing")
    B_rows = _read_log_rows(
        c,
        loop_chart.get("B_over_Pstar_over_frozen_E0"),
        4,
        chart + ".canonical B_over_Pstar_over_frozen_E0",
    )
    logPE_source = (
        _path(native_chart, "inverse_Utheta_detail.log_Utheta_upper")
        if chart in LONG_REFERENCE
        else _path(native_chart, "scales.log_Utheta_upper")
    )
    logPE = _read_interval(c, logPE_source, chart + ".log_Utheta_upper")
    native_U = _read_log_rows(c, _path(native_chart, "native_caps.U"), 5,
                              chart + ".native U caps")
    native_M = _read_log_rows(c, _path(native_chart, "native_caps.M"), 5,
                              chart + ".native M caps")
    delta_U = [None if row is None else logPE + row - logN for row in B_rows]
    delta_M = [_max_upper_point(c, incoming[index], delta_U[index])
               for index in range(4)]

    # The native Q and finite-N delta Q both come from the new pointwise
    # operator.  The old saved Q caps are intentionally never read.
    native_Q = kernel.pressure_Q(c, native_U, native_M, z, delta, inverseL)
    delta_Q = kernel.pressure_Q(c, delta_U, delta_M, z, delta, inverseL)
    six_logs = kernel.six_differences(c, delta_U, delta_Q,
                                      native_U, native_Q)
    density_defects = {
        key: [None if value is None else value + mass_log for value in values]
        for key, values in six_logs.items()
    }
    # Boundary endpoint differences consume the pointwise six differences.
    # The normalized mass belongs only to the integrated defects used by D;
    # applying it here would introduce an extra factor of one half in B.
    B_error = kernel.boundary_difference(c, six_logs, z, delta, inverseL)
    D_result = kernel.D_difference(c, density_defects, z, delta, inverseL)
    D_error = D_result["rows"]
    S1_error = [_log_add(kernel, c, B_error[index], D_error[index])
                for index in I1_ORDERS]
    threshold = -c.mpf(1024) * c.ln(c.mpf(2))
    strict_six = {
        key: all(_strictly_below(c, value, threshold) for value in values)
        for key, values in density_defects.items()
    }
    strict_B = all(_strictly_below(c, value, threshold) for value in B_error)
    target = all(strict_six.values()) and strict_B

    base = _baseline_data(c, baseline_chart, chart)
    full = _full_signed_rows(c, base, radius, normalization, source_path) \
        if target else None
    baseline = _baseline_receipt(base, source_path)
    packed_six = {key: _pack_log_rows(c, values, "six." + key)
                  for key, values in density_defects.items()}
    packed_D = _pack_log_rows(c, D_error, "D_difference")
    packed_S1 = _pack_log_rows(c, S1_error, "S1_difference")
    packed_B = _pack_log_rows(c, B_error, "B_difference")
    coefficient_receipt = {
        "coefficient_rows": {
            key: _pack_signed_rows(c, values, "D.coefficient." + key)
            for key, values in D_result["coefficient_rows"].items()
        },
        "coefficient_log_rows": {
            key: _pack_log_rows(c, values, "D.coefficient_log." + key)
            for key, values in D_result["coefficient_log_rows"].items()
        },
        "term_log_bounds": {
            key: _pack_log_rows(c, values, "D.term." + key)
            for key, values in D_result["term_log_bounds"].items()
        },
        "formula": D_result.get("formula"),
        "source_formula": D_result.get("source_formula"),
        "composed_from_six_I1": D_result.get("composed_from_six_I1") is True,
        "signed_coefficient_source": "new kernel geometric receipt; saved coefficients retained for full D",
    }
    log_W_record = baseline_chart.get("log_W")
    _require(isinstance(log_W_record, Mapping),
             chart + ".log_W baseline record is missing")
    chart_packet = {
        "passed": True,
        "chart": chart,
        "normalization": normalization,
        "normalized_by": normalization,
        "ordinary_Z_orders": list(ORDERS[:4]),
        "I1_ordinary_Z_orders": list(I1_ORDERS),
        "physical_y_derivatives_claimed": False,
        "physical_error_log_bounds": {
            "delta_U4": _pack_log_rows(c, delta_U, chart + ".delta_U"),
            "delta_M4": _pack_log_rows(c, delta_M, chart + ".delta_M"),
            "delta_Q3": _pack_log_rows(c, delta_Q, chart + ".delta_Q"),
            "delta_U": _pack_log_rows(c, delta_U, chart + ".delta_U"),
            "delta_M": _pack_log_rows(c, delta_M, chart + ".delta_M"),
            "delta_Q": _pack_log_rows(c, delta_Q, chart + ".delta_Q"),
            "logN": _record(c, logN, chart + ".logN"),
            "log_Pstar_E_upper": _record(c, logPE, chart + ".log_Utheta_upper"),
        },
        "six_I1_difference_error_log_bounds": packed_six,
        "B_difference_error_log_bounds": packed_B,
        "D_difference_error_log_bounds": packed_D,
        "S1_difference_error_log_bounds": packed_S1,
        "target_2_power_minus1024_met": target,
        "strict_error_checks": {
            "six_I1": strict_six,
            "B": strict_B,
            "threshold": _record(c, threshold, chart + ".I1_error_threshold"),
            "own_W": normalization,
        },
        "I1_error_threshold": _record(c, threshold,
                                       chart + ".I1_error_threshold"),
        "baseline_centres": baseline,
        "I1_full": None if full is None else full["I1_full"],
        "B_full": None if full is None else full["B_full"],
        "D_full": None if full is None else full["D_full"],
        "S1_full": None if full is None else full["S1_full"],
        "signed_I1_radius": None if full is None else
            _record(c, radius, chart + ".radius"),
        "coefficient_receipt": coefficient_receipt,
        "normalized_mass": _record(c, c.mpf(1) / 2, chart + ".J0"),
        "geometry": {
            "normalization": normalization,
            "log_W": copy.deepcopy(log_W_record),
            "log_W_source": source_path + ".log_W",
            "normalized_mass_formula": "J0=1/2 safe universal bound",
            "normalized_mass_source": "safe universal 1/2; no physical W re-applied",
            "own_weight_used": True,
            "physical_Pstar_applied_once": True,
            "no_extra_log_W": True,
        },
        "logW_provenance": {
            "source_path": source_path + ".log_W",
            "record": copy.deepcopy(log_W_record),
            "exact_source_record": True,
            "normalization": normalization,
            "stored_log_W_retained_as_provenance": True,
            "log_W_used_in_error_arithmetic": False,
            "physical_W_applied_again": False,
        },
        "source_row_correspondence": {
            "loop_B": "loop_report.actual_directed_loop_controls.charts." + chart
                       + ".B_over_Pstar_over_frozen_E0",
            "native_U_M": "native_report.actual_native_directed_controls_v2.raw_interval_tree"
                          + (".long_reference" if chart in LONG_REFERENCE else ".restoration")
                          + ".charts." + chart + ".native_caps",
            "native_Q": "new_kernel.pressure_Q(native_U5,native_M5,shared_z,shared_delta,shared_inverse_L5)",
            "delta_Q": "new_kernel.pressure_Q(delta_U4,delta_M4,shared_z,shared_delta,shared_inverse_L5)",
            "baseline": source_path + ".arithmetic",
            "corrected_postrepair_identifier": (
                "actual_restoration_post_jet_repair.corrected_background"
                if chart in RESTORATION else None
            ),
        },
        "source_contract": {
            "passed": True,
            "chart": chart,
            "actual_N_handle": "N18 graph handle",
            "actual_N_not_materialized": True,
            "phase_Z_independent": True,
            "physical_delta_U_local_only": True,
            "incoming_delta_U_not_added": True,
            "incoming_delta_M_retained": True,
            "nonreset_positive_kernel_mean": True,
            "mean_bound": "max(incoming deltaM, local deltaU)",
            "Q_source": "new kernel pointwise pressure_Q",
            "integrated_Q_used": False,
            "normalized_mass": "J0=1/2",
            "own_W": normalization,
            "finite_N_radius_added_once": bool(target),
            "signed_error_caps_not_selected_values": True,
            "physical_y_derivatives_claimed": False,
        },
    }
    return chart_packet, delta_M


def compute(c: MPIntervalContext, loop_report_bytes: bytes,
            native_report_bytes: bytes) -> dict[str, Any]:
    """Compute the JSON finite-``N`` packet from authenticated report bytes."""

    _require(isinstance(c, MPIntervalContext), "compute requires an interval context")
    loop, native = _admit_reports(loop_report_bytes, native_report_bytes)
    # The import is deliberately lazy: report admission and all source checks
    # above remain usable even when the pure kernel is being reviewed alone.
    from . import original_directed_finite_N_kernel as kernel

    native_packet = _path(native, "actual_native_directed_controls_v2")
    tree = _path(native_packet, "raw_interval_tree")
    shared = _path(tree, "shared")
    z = _read_interval(c, shared.get("z"), "shared.z")
    delta = _read_interval(c, shared.get("delta"), "shared.delta")
    logN = _read_interval(c, shared.get("logN"), "shared.logN")
    inverseL = _read_log_rows(c, shared.get("inverse_L_logs5"), 5,
                              "shared.inverse_L_logs5")

    incoming_path = (
        "native_report.actual_switch_post_finite_N_I1_numeric.charts.post_switch"
        ".physical_error_log_bounds.delta_M"
    )
    incoming_source = _path(native, incoming_path.replace("native_report.", ""))
    incoming = _read_log_rows(c, incoming_source, 4, incoming_path)
    incoming_copy = copy.deepcopy(incoming_source)
    radius = kernel.safeiv(c, Fraction(1, 2 ** 1024))
    mass_log = c.ln(c.mpf(1) / 2)
    loop_packet = _path(loop, "actual_directed_loop_controls")
    loop_charts = _path(loop_packet, "charts")
    _require(list(loop_charts.keys()) == list(CHARTS),
             "Loop chart order changed")

    chart_packets: dict[str, dict[str, Any]] = {}
    incoming_current = incoming
    for chart in CHARTS:
        native_chart = _source_packet(native, chart)
        baseline, source_path = _baseline_chart(native, chart)
        chart_packet, outgoing = _chart(
            c, kernel, chart, loop_charts[chart], native_chart, baseline,
            source_path, incoming_current, z, delta, inverseL, logN, mass_log,
            radius, "native_report.actual_native_directed_controls_v2.raw_interval_tree.shared.logN",
        )
        chart_packet["incoming_delta_M"] = _pack_log_rows(
            c, incoming_current, chart + ".incoming_delta_M")
        chart_packet["memory_transfer"] = {
            "incoming_rows": 4,
            "outgoing_rows": 4,
            "no_memory_reset": True,
            "outgoing_becomes_next_incoming": chart != CHARTS[-1],
        }
        chart_packets[chart] = chart_packet
        incoming_current = outgoing

    quantitative = all(
        chart_packets[chart]["target_2_power_minus1024_met"] for chart in CHARTS
    )
    packet = {
        "schema": SCHEMA,
        "passed": True,
        "quantitative_passed": quantitative,
        "all_normalized_I1_and_B_errors_below_2_power_minus1024": quantitative,
        "charts": chart_packets,
        "chart_order": list(CHARTS),
        "ordinary_Z_orders": list(ORDERS),
        "I1_ordinary_Z_orders": list(I1_ORDERS),
        "Q_ordinary_Z_orders": list(Q_ORDERS),
        "parent_reports": {
            "loop": {"path": LOOP_REPORT_NAME, "sha256": LOOP_REPORT_SHA256,
                     "bytes": LOOP_REPORT_BYTES},
            "native": {"path": NATIVE_REPORT_NAME, "sha256": NATIVE_REPORT_SHA256,
                       "bytes": NATIVE_REPORT_BYTES},
            "exact_bytes_authenticated_before_json": True,
            "full_reports_only": True,
        },
        "shared_inputs": {
            "z_source": "native_report.actual_native_directed_controls_v2.raw_interval_tree.shared.z",
            "delta_source": "native_report.actual_native_directed_controls_v2.raw_interval_tree.shared.delta",
            "logN_source": "native_report.actual_native_directed_controls_v2.raw_interval_tree.shared.logN",
            "inverse_L_source": "native_report.actual_native_directed_controls_v2.raw_interval_tree.shared.inverse_L_logs5",
            "actual_N": "N18 graph handle",
            "actual_N_not_materialized": True,
            "inverse_L_rows": 5,
            "shared_logN_consumed_once": True,
        },
        "incoming_delta_M_source": {
            "path": incoming_path,
            "rows": incoming_copy,
            "rows_exact_source_copy": True,
            "initial_rows": 4,
            "memory_not_reset": True,
        },
        "radius_policy": {
            "radius": "2^-1024",
            "symmetric_radius_added_once_per_six_I1_C1_and_B_C1": True,
            "D_and_S1_extra_radius": False,
            "failed_target_keeps_baseline_and_error_logs": True,
        },
        "scope": {
            "finite_N_error_propagation_only": True,
            "conditional_on_hash_bound_incoming_cap": True,
            "upstream_incoming_cap_audit_complete": False,
            "current_N_certified": False,
            "full_pressure_revalidated": False,
            "global_five_moment_closure": False,
            "physical_y_bounds": False,
            "Rm_global_claim": False,
            "full_goal_complete": False,
        },
        "source_contract": {
            "saved_JSON_source_hash_AST_reads_only": True,
            "new_kernel_only": True,
            "accepted_builders_or_producers_called": False,
            "accepted_owners_or_checkers_called": False,
            "accepted_selfchecks_or_recurrences_called": False,
            "quadratures_or_norms_called": False,
            "legacy_finite_N_functions_called": False,
            "stored_huge_exponentials_materialized": False,
            "caps_promoted_to_signed_values": False,
            "native_Q_recomputed_from_new_U_M": True,
            "old_Q_caps_used": False,
            "boundary_from_density_not_integrated": True,
            "D_from_integral_defects": True,
            "S1_from_B_plus_D": True,
            "own_W_labels": dict(NORMALIZATIONS),
        },
        "execution_audit": {
            "source_executed": False,
            "accepted_ancestor_math_replayed": False,
            "accepted_builder_producer_owner_checker_replayed": False,
            "recurrence_quadrature_norm_execution": False,
            "new_pure_kernel_on_changed_inputs": True,
            "fresh_fixture_scope_only": True,
        },
        "limits": {
            "upstream_incoming_delta_M_audit": "required before current_N certification",
            "current_N_certified": False,
            "physical_y_derivatives_claimed": False,
            "Rm_external_pressure_evaluated": False,
            "global_five_moment_closure": False,
            "full_goal_complete": False,
        },
    }
    json.dumps(packet, allow_nan=False)
    return packet


def fresh_fixture() -> dict[str, Any]:
    """Independent rational memory/radius/composition smoke fixture."""

    from . import original_directed_finite_N_kernel as kernel

    c = MPIntervalContext()
    c.prec = 256
    incoming = [
        kernel.safeiv(c, Fraction(1, 10)),
        kernel.safeiv(c, Fraction(1, 20)),
    ]
    local = [
        kernel.safeiv(c, Fraction(1, 20)),
        kernel.safeiv(c, Fraction(1, 5)),
    ]
    maximum = [_max_upper_point(c, incoming[k], local[k]) for k in range(2)]
    memory_checks = [maximum[0].a == incoming[0].b,
                     maximum[1].a == local[1].b]
    radius = kernel.safeiv(c, Fraction(1, 2 ** 1024))
    baseline = [
        kernel.safeiv(c, Fraction(1, 3)),
        kernel.safeiv(c, -Fraction(1, 7)),
    ]
    full = [value + _symmetric(c, radius) for value in baseline]
    radius_checks = [full[0].a < baseline[0].a < full[0].b,
                     full[1].a < baseline[1].a < full[1].b]
    # The kernel is exercised only on fresh rational changed inputs.  This
    # confirms the public D receipt shape without loading either saved report.
    z = kernel.safeiv(c, Fraction(1, 4))
    delta = kernel.safeiv(c, Fraction(1, 10))
    inverseL = [c.mpf(0), c.mpf(0), c.mpf(0), c.mpf(0), c.mpf(0)]
    defects = {name: [c.mpf(-80), c.mpf(-80)] for name in I1_NAMES}
    d = kernel.D_difference(c, defects, z, delta, inverseL)
    composition_checks = [
        isinstance(d, Mapping),
        set(d.get("coefficient_rows", {})) == set(COEFFICIENT_NAMES),
        len(d.get("rows", ())) == 2,
    ]
    result = {
        "passed": all(memory_checks + radius_checks + composition_checks),
        "checks": {
            "max_incoming_preserved": memory_checks[0],
            "max_local_selected": memory_checks[1],
            "symmetric_radius_contains_positive_baseline": radius_checks[0],
            "symmetric_radius_contains_negative_baseline": radius_checks[1],
            "D_coefficient_receipt_shape": composition_checks[1],
            "D_rows_are_C1": composition_checks[2],
        },
        "memory_formula": "max(incoming, local) by exact upper endpoint",
        "radius_formula": "baseline + [-2^-1024, 2^-1024]",
        "composition_formula": "D=sum(C_i*deltaI1_i)",
        "radius": _record(c, radius, "fixture.radius"),
        "source_executed": False,
        "saved_reports_loaded": False,
    }
    json.dumps(result, allow_nan=False)
    _require(result["passed"], "Fresh finite-N adapter fixture failed")
    return result


fixture_bound_check = fresh_fixture


__all__ = [
    "MODULE",
    "KERNEL_MODULE",
    "SCHEMA",
    "CHARTS",
    "compute",
    "fresh_fixture",
    "fixture_bound_check",
]
