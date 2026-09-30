"""Finite high-degree regular-core extension experiment.

This is a bounded numerical experiment around the Section 8 axis data.  It
uses the existing exact-equation radial recurrence in
``lei_ren_part1_paper_core_recursion`` and does not replace that recurrence by
a fitted radial profile.  The source pressure is a prototype only: its axis
amplitude is taken from the finite source outer-pressure quadrature at Z = 0
and extended as K / (1 + Z^2)^2.  The receipt records the omitted future-tail
and finite-quadrature uncertainty explicitly.

The public ``build_coefficients`` function is intended for the parent adapter.
It accepts a shared pressure callback when one is available; otherwise it
uses the labelled prototype pressure described above.  F0 Z-jets are generated
from the exact rational gradient

    d_Z log(F0) = -Lambda L H0 / (H0^2 + sigma0^2),

with normalized Taylor arithmetic.  No scale division is applied to the
returned coefficients, so the recurrence sees the physical F0/P0 relative
scaling.  A future caller may rescale the complete coupled equations, but a
partial F-only rescaling would change the pressure terms and is forbidden.

This file makes no analytic-contraction, global-core, cone, matching, or
completed-pressure claim.
"""

from __future__ import annotations

from decimal import Decimal, getcontext
import json
from pathlib import Path
import sys
from typing import Any, Callable

import mpmath as mp


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lei_ren_part1_paper_axis_jets import RegularCoreAxisJets  # noqa: E402
from lei_ren_part1_paper_core_recursion import (  # noqa: E402
    core_coefficients,
    core_equation_defects,
    evaluate_core_jets,
    pressure_taylor,
)
from lei_ren_part1_paper_outer import PaperOuterSchedule  # noqa: E402
from lei_ren_part1_paper_outer_pressure import normalized_outer_pressure  # noqa: E402


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"


def _mp(value: Any) -> mp.mpf:
    if isinstance(value, mp.mpf):
        return mp.mpf(value)
    return mp.mpf(str(value))


def _nstr(value: Any, digits: int = 40) -> str:
    if value is None:
        return "null"
    return mp.nstr(value if isinstance(value, mp.mpf) else _mp(value), digits)


def _signed_log(value: Any, digits: int = 50) -> dict[str, Any]:
    with mp.workdps(max(80, digits + 20)):
        x = _mp(value)
        if x == 0:
            return {"sign": 0, "log_abs": None, "value": "0"}
        return {
            "sign": 1 if x > 0 else -1,
            "log_abs": _nstr(mp.log(abs(x)), digits),
            "value": _nstr(x, digits),
        }


def _relative(value: mp.mpf, reference: mp.mpf) -> mp.mpf:
    return abs(value - reference) / max(abs(reference), mp.mpf("1e-200"))


def _relative_unfloored(value: mp.mpf, reference: mp.mpf) -> mp.mpf:
    """Relative difference that remains meaningful for exp(-10^18) data."""

    if reference == 0:
        return abs(value)
    return abs(value - reference) / abs(reference)


def _series_zero(length: int) -> list[mp.mpf]:
    return [mp.mpf(0)] * length


def _series_add(*items: list[mp.mpf]) -> list[mp.mpf]:
    length = len(items[0])
    return [sum(item[k] for item in items) for k in range(length)]


def _series_scale(item: list[mp.mpf], scalar: Any) -> list[mp.mpf]:
    c = _mp(scalar)
    return [c * value for value in item]


def _series_mul(left: list[mp.mpf], right: list[mp.mpf]) -> list[mp.mpf]:
    length = min(len(left), len(right))
    return [
        sum(left[i] * right[k - i] for i in range(k + 1))
        for k in range(length)
    ]


def _series_inv(item: list[mp.mpf]) -> list[mp.mpf]:
    if item[0] == 0:
        raise ZeroDivisionError("Taylor reciprocal has zero constant term")
    out = [1 / item[0]]
    for k in range(1, len(item)):
        out.append(-sum(item[i] * out[k - i] for i in range(1, k + 1)) / item[0])
    return out


def _series_coordinate(center: mp.mpf, length: int) -> list[mp.mpf]:
    out = _series_zero(length)
    out[0] = center
    if length > 1:
        out[1] = mp.mpf(1)
    return out


def _axis_gradient_series(
    center: mp.mpf,
    *,
    length: int,
    delta: mp.mpf,
    j: mp.mpf,
    sigma0: mp.mpf,
) -> list[mp.mpf]:
    """Return normalized Taylor coefficients of L H/(H^2+sigma^2)."""

    z = _series_coordinate(center, length)
    one = [mp.mpf(1)] + [mp.mpf(0)] * (length - 1)
    z2 = _series_mul(z, z)
    d = _series_add(one, _series_scale(z2, -1))
    L = _series_add(one, _series_scale(z2, -delta))
    u0 = _series_add(_series_scale(z, 4), [j] + [mp.mpf(0)] * (length - 1))
    H = _series_add(
        _series_scale(z, (1 - delta) / 2),
        _series_mul(d, u0),
    )
    denominator = _series_add(_series_mul(H, H), [sigma0 * sigma0] + [mp.mpf(0)] * (length - 1))
    return _series_mul(_series_mul(L, H), _series_inv(denominator))


def _axis_F0_taylor(axis: RegularCoreAxisJets, center: Any, length: int) -> list[mp.mpf]:
    """Build F0^(k)/k! from the exact rational logarithmic derivative."""

    with mp.workdps(axis.precision):
        z = _mp(center)
        data = axis.axis_data(z)
        q = _axis_gradient_series(
            z,
            length=length,
            delta=axis.delta,
            j=axis.j,
            sigma0=axis.sigma0,
        )
        f = _series_zero(length)
        f[0] = data["F0"]
        for n in range(length - 1):
            f[n + 1] = -axis.Lambda / (n + 1) * sum(
                q[k] * f[n - k] for k in range(n + 1)
            )
        return f


def _axis_u0_taylor(center: Any, *, j: mp.mpf, length: int) -> list[mp.mpf]:
    out = _series_zero(length)
    out[0] = 4 * _mp(center) + j
    if length > 1:
        out[1] = mp.mpf(4)
    return out


class PrototypeAxisPressure:
    """Initial source pressure prototype with explicit provenance metadata."""

    def __init__(
        self,
        *,
        logPstar: Any = "14",
        logRref: Any,
        delta: Any = "1e-32",
        Md: Any = ".5",
        c_mu: Any = ".001",
        c_delta: Any = ".001",
        c_epsilon: Any = ".01",
        pressure_order: int = 64,
        precision: int = 180,
    ) -> None:
        self.precision = max(100, int(precision))
        self.schedule = PaperOuterSchedule(
            logPstar=str(logPstar),
            logRref=str(logRref),
            delta=str(delta),
            Md=str(Md),
            c_mu=str(c_mu),
            c_delta=str(c_delta),
            c_epsilon=str(c_epsilon),
        )
        self.pressure_order = int(pressure_order)
        if self.pressure_order < 16:
            raise ValueError("pressure_order must be at least 16")
        coarse = normalized_outer_pressure(self.schedule, 0.0, order=max(16, self.pressure_order // 2))
        fine = normalized_outer_pressure(self.schedule, 0.0, order=self.pressure_order)
        with mp.workdps(self.precision):
            self.logPstar = _mp(logPstar)
            self.logP0_norm = _mp(fine["reference_extension_P0_over_Pstar_squared"])
            self.logP0_norm_coarse = _mp(coarse["reference_extension_P0_over_Pstar_squared"])
            self.Pstar_squared = mp.exp(2 * self.logPstar)
            self.P0_at_zero = self.logP0_norm * self.Pstar_squared
            self.quadrature_difference = abs(self.logP0_norm - self.logP0_norm_coarse) * self.Pstar_squared
            self.omitted_tail_bound = _mp(fine["omitted_outer_pressure_upper_bound_over_Pstar_squared"]) * self.Pstar_squared
        self.metadata = {
            "kind": "initial_source_axis_amplitude_only",
            "pressure_order_coarse": max(16, self.pressure_order // 2),
            "pressure_order_fine": self.pressure_order,
            "P0_over_Pstar_squared_at_Z0": _nstr(self.logP0_norm, 40),
            "P0_at_Z0": _signed_log(self.P0_at_zero),
            "quadrature_difference_absolute": _nstr(self.quadrature_difference, 30),
            "omitted_future_tail_bound_absolute": _nstr(self.omitted_tail_bound, 30),
            "future_tail_bound_included": True,
            "full_source_pressure_exact": False,
            "pressure_formula": "P0(Z)=P0(0)/(1+Z^2)^2",
        }

    def __call__(self, Z: Any) -> mp.mpf:
        with mp.workdps(self.precision):
            z = _mp(Z)
            if abs(z) > 1:
                raise ValueError("|Z| must be <= 1")
            return self.P0_at_zero / (1 + z * z) ** 2

    def taylor(self, Z: Any, *, length: int) -> list[mp.mpf]:
        with mp.workdps(self.precision):
            z = _mp(Z)
            if abs(z) > 1:
                raise ValueError("|Z| must be <= 1")
            coordinate = _series_coordinate(z, length)
            one = [mp.mpf(1)] + [mp.mpf(0)] * (length - 1)
            denominator = _series_mul(_series_add(one, _series_mul(coordinate, coordinate)),
                                      _series_add(one, _series_mul(coordinate, coordinate)))
            return _series_scale(_series_inv(denominator), self.P0_at_zero)


def _logRref_from_source(*, logPstar: Any, logC: mp.mpf) -> mp.mpf:
    with mp.workdps(220):
        return mp.log(110) + 10 * (logC + _mp(logPstar))


def _default_configuration(
    *,
    precision: int = 180,
    Lambda_value: Any = "1e18",
    pressure: PrototypeAxisPressure | Any | None = None,
) -> tuple[RegularCoreAxisJets, PrototypeAxisPressure | Any, dict[str, Any]]:
    with mp.workdps(precision):
        j = _mp(".02")
        Lambda = _mp(Lambda_value)
        logC = 2 * mp.log(Lambda)
        logPstar = _mp("14")
        logRref = _logRref_from_source(logPstar=logPstar, logC=logC)
    pressure_reused = pressure is not None
    if pressure is None:
        pressure = PrototypeAxisPressure(
            logPstar=_nstr(logPstar, precision),
            logRref=_nstr(logRref, precision),
            delta="1e-32",
            precision=precision,
        )
    axis = RegularCoreAxisJets(
        j=_nstr(j, precision),
        Lambda=_nstr(Lambda, precision),
        # Recompute the equality inside RegularCoreAxisJets' own MP context;
        # str(mp.mpf) outside that context can round the lower-bound value.
        logC=None,
        delta="1e-32",
        pressure=pressure,
        precision=precision,
        pressure_step="1e-6",
    )
    with mp.workdps(precision):
        metadata = {
            "j": _nstr(j, precision),
            "Lambda": _nstr(Lambda, precision),
            "logC": _nstr(logC, precision),
            "logPstar": _nstr(logPstar, precision),
            "logRref": _nstr(logRref, precision),
            "delta": "1e-32",
            "R_a": _nstr(4 / Lambda, precision),
            "pressure_reused_from_other_configuration": pressure_reused,
            "s_domain": "0 <= Lambda*R <= 4.1",
        }
    return axis, pressure, metadata


def build_coefficients(
    Z: Any,
    radial_degree: int,
    *,
    axis: RegularCoreAxisJets | None = None,
    pressure: Any | None = None,
    pressure_taylor_coefficients: list[Any] | None = None,
    pressure_spacing: Any = ".01",
    precision: int = 180,
) -> dict[str, Any]:
    """Build one high-degree local radial coefficient set.

    ``pressure`` may be the actual shared axis-pressure callback.  If
    ``pressure_taylor_coefficients`` is supplied it takes precedence and is
    interpreted in normalized Z-Taylor convention.  With neither input, the
    prototype source pressure is used and marked in the returned metadata.
    """

    degree = int(radial_degree)
    if degree < 1 or degree != radial_degree:
        raise ValueError("radial_degree must be a positive integer")
    with mp.workdps(precision):
        if axis is None:
            axis, default_pressure, _ = _default_configuration(precision=precision)
        else:
            default_pressure = getattr(axis, "pressure", None)
        z = _mp(Z)
        length = degree + 2
        f0 = _axis_F0_taylor(axis, z, length)
        u0 = _axis_u0_taylor(z, j=axis.j, length=length)
        if pressure_taylor_coefficients is not None:
            p0 = [_mp(value) for value in pressure_taylor_coefficients[:length]]
            if len(p0) < length:
                raise ValueError("pressure_taylor_coefficients are too short")
            pressure_source = "caller_supplied_normalized_Z_taylor"
        else:
            provider = pressure if pressure is not None else default_pressure
            if provider is None:
                raise ValueError("actual pressure callback or prototype pressure is required")
            prototype = provider if isinstance(provider, PrototypeAxisPressure) else None
            if prototype is not None:
                p0 = prototype.taylor(z, length=length)
                pressure_source = "prototype_P0_zero_amplitude_with_(1+Z^2)^-2_extension"
            else:
                # The existing helper requires an even interpolation degree;
                # request the next even degree and retain the needed prefix.
                fd_degree = length if length % 2 == 0 else length + 1
                p0 = pressure_taylor(
                    lambda point: provider(point),
                    z,
                    degree=fd_degree,
                    spacing=pressure_spacing,
                    precision=precision,
                )[:length]
                pressure_source = "shared_pressure_callback_fourth_or_higher_order_local_Taylor"
        result = core_coefficients(
            z,
            axis.delta,
            F0_Z_taylor=f0,
            U0_Z_taylor=u0,
            P0_Z_taylor=p0,
            radial_degree=degree,
            precision=precision,
        )
        result["axis_metadata"] = axis.metadata()
        result["pressure_source"] = pressure_source
        result["pressure_metadata"] = getattr(default_pressure, "metadata", None)
        result["source"] = SOURCE
        result["source_version"] = SOURCE_VERSION
        result["prototype_pressure_uncertainty"] = not pressure_source.startswith("shared_") and pressure_taylor_coefficients is None
        result["F0_Z_taylor"] = f0
        result["U0_Z_taylor"] = u0
        result["P0_Z_taylor"] = p0
        return result


def _kappa_from_evaluation(evaluation: dict[str, mp.mpf], R: mp.mpf) -> mp.mpf:
    F = evaluation["F"]
    FR = evaluation["F_R"]
    UzR = evaluation["Uz_R"]
    if R <= 0 or F == 0 or FR == 0:
        return mp.nan
    s_theta = 2 * R * FR
    s_z = mp.sqrt(2 * R) * UzR
    return -(s_theta * s_theta + s_z * s_z) / (F * s_theta)


def diagnose_candidate(coefficients: dict[str, Any], s: Any) -> dict[str, Any]:
    """Evaluate one candidate point at s = Lambda R."""

    with mp.workdps(max(100, int(coefficients["precision"]))):
        Lambda = _mp(coefficients["axis_metadata"]["Lambda"])
        R = _mp(s) / Lambda
        evaluation = evaluate_core_jets(coefficients, R)
        defects = core_equation_defects(coefficients, R)
        kappa = _kappa_from_evaluation(evaluation, R)
        scales = {
            "angular": max(abs(defects["angular"]), abs(evaluation["F"]), mp.mpf("1e-200")),
            "axial": max(abs(defects["axial"]), abs(evaluation["Uz"]), mp.mpf("1e-200")),
            "pressure": max(abs(defects["pressure"]), abs(evaluation["P"]), mp.mpf("1e-200")),
        }
        return {
            "s": str(s),
            "R": _nstr(R, 40),
            "F": _signed_log(evaluation["F"]),
            "F_R": _signed_log(evaluation["F_R"]),
            "Uz": _signed_log(evaluation["Uz"]),
            "Uz_R": _signed_log(evaluation["Uz_R"]),
            "P": _signed_log(evaluation["P"]),
            "kappa": _signed_log(kappa),
            "F_positive": bool(evaluation["F"] > 0),
            "F_R_negative": bool(evaluation["F_R"] < 0),
            "defects": {name: _signed_log(value) for name, value in defects.items()},
            "defect_relative": {
                "angular": _nstr(abs(defects["angular"]) / scales["angular"], 30),
                "axial": _nstr(abs(defects["axial"]) / scales["axial"], 30),
                "pressure": _nstr(abs(defects["pressure"]) / scales["pressure"], 30),
            },
        }


def run_experiment() -> dict[str, Any]:
    precision = 180
    axis, pressure, configuration = _default_configuration(precision=precision)
    degrees = (12, 18, 24)
    with mp.workdps(precision):
        root_z = mp.nstr(axis.Z0, precision)
        narrow_z = mp.nstr(axis.Z0 + 4 * axis._root_width, precision)
    z_values = ["0.3", "0.7", root_z, narrow_z]
    s_values = ("0.1", "1", "2", "4.1")
    rows: list[dict[str, Any]] = []
    coefficient_sets: dict[tuple[str, int], dict[str, Any]] = {}
    for z in z_values:
        for degree in degrees:
            coefficients = build_coefficients(
                z,
                degree,
                axis=axis,
                pressure=pressure,
                precision=precision,
            )
            coefficient_sets[(z, degree)] = coefficients
            for s in s_values:
                row = diagnose_candidate(coefficients, s)
                row.update({"Z": z, "radial_degree": degree})
                rows.append(row)
    convergence: list[dict[str, Any]] = []
    for z in z_values:
        for s in s_values:
            raw_values = {}
            values = {}
            for degree in degrees:
                evaluation = evaluate_core_jets(
                    coefficient_sets[(z, degree)],
                    _mp(s) / _mp(configuration["Lambda"]),
                )
                raw_values[degree] = evaluation
                values[str(degree)] = {
                    key: _signed_log(evaluation[key])
                    for key in ("F", "F_R", "Uz", "Uz_R", "P")
                }
            pairwise = {}
            for lower, upper in zip(degrees, degrees[1:]):
                pairwise[f"{lower}_to_{upper}"] = {
                    key: _nstr(
                        _relative_unfloored(raw_values[upper][key], raw_values[lower][key]),
                        30,
                    )
                    for key in ("F", "F_R", "Uz", "Uz_R", "P")
                }
            convergence.append({"Z": z, "s": str(s), "values": values, "relative_differences": pairwise})
    # Keep the physical prototype pressure fixed while changing Lambda.  This
    # isolates the radial scale effect from a simultaneous pressure change.
    lambda_sweep: list[dict[str, Any]] = []
    for lambda_value in ("1e20", "1e22", "1e24"):
        sweep_axis, sweep_pressure, sweep_configuration = _default_configuration(
            precision=precision,
            Lambda_value=lambda_value,
            pressure=pressure,
        )
        sweep_z = root_z
        sweep_rows = []
        for degree in (12, 18):
            sweep_coefficients = build_coefficients(
                sweep_z,
                degree,
                axis=sweep_axis,
                pressure=sweep_pressure,
                precision=precision,
            )
            for s in ("2", "4.1"):
                row = diagnose_candidate(sweep_coefficients, s)
                row.update({"radial_degree": degree})
                sweep_rows.append(row)
        lambda_sweep.append(
            {
                "Lambda": lambda_value,
                "configuration": sweep_configuration,
                "Z": sweep_z,
                "pressure_reused_from_baseline": True,
                "rows": sweep_rows,
            }
        )

    # Preserve the earlier 15-digit literal as an intentional nearby-point
    # probe.  It is not the exact root and must not be discarded as a logging
    # artefact: the narrow layer has width sigma0/|H0_Z|, and this offset is a
    # genuine source-coordinate sample.
    rounded_root_literal = "-0.00444443468970541"
    rounded_root_probe: list[dict[str, Any]] = []
    with mp.workdps(precision):
        rounded_offset = _mp(rounded_root_literal) - axis.Z0
    for lambda_value in ("1e18", "1e20", "1e22", "1e24"):
        if lambda_value == "1e18":
            probe_axis, probe_pressure, probe_configuration = axis, pressure, configuration
        else:
            probe_axis, probe_pressure, probe_configuration = _default_configuration(
                precision=precision,
                Lambda_value=lambda_value,
                pressure=pressure,
            )
        probe_coefficients = build_coefficients(
            rounded_root_literal,
            18,
            axis=probe_axis,
            pressure=probe_pressure,
            precision=precision,
        )
        probe_rows = []
        for s in ("2", "4.1"):
            row = diagnose_candidate(probe_coefficients, s)
            row.update({"radial_degree": 18})
            probe_rows.append(row)
        rounded_root_probe.append(
            {
                "Lambda": lambda_value,
                "Z_literal": rounded_root_literal,
                "offset_from_exact_root": _nstr(rounded_offset, 40),
                "configuration": probe_configuration,
                "pressure_reused_from_baseline": lambda_value != "1e18",
                "rows": probe_rows,
            }
        )

    # Test the much larger Lambda regime against the entire requested narrow
    # layer.  The offset scale is the source diagnostic
    # |U1(Z0)| / (Lambda |H0_Z(Z0)|), not an arbitrary grid spacing.
    extreme_lambda = "1e36"
    extreme_axis, extreme_pressure, extreme_configuration = _default_configuration(
        precision=precision,
        Lambda_value=extreme_lambda,
        pressure=pressure,
    )
    with mp.workdps(precision):
        root_jet = axis.axis_jets(axis.Z0)
        u1_abs = abs(root_jet.Uz_R)
        h0_z_abs = abs(axis.H0_Z(axis.Z0))
        offset_unit = u1_abs / (_mp(extreme_lambda) * h0_z_abs)
        dangerous_points = {
            "exact_root": (mp.nstr(axis.Z0, precision), "0"),
            "rounded_literal": (rounded_root_literal, None),
        }
        for c in (-16, -4, -1, 0, 1, 4, 16):
            point = axis.Z0 + mp.mpf(c) * offset_unit
            dangerous_points[f"offset_c_{c}"] = (mp.nstr(point, precision), str(c))
    extreme_rows: list[dict[str, Any]] = []
    extreme_convergence: list[dict[str, Any]] = []
    for label, (point_z, c_value) in dangerous_points.items():
        point_coefficients: dict[int, dict[str, Any]] = {}
        for degree in (18, 24):
            point_coefficients[degree] = build_coefficients(
                point_z,
                degree,
                axis=extreme_axis,
                pressure=extreme_pressure,
                precision=precision,
            )
            for s in ("2", "4.1"):
                row = diagnose_candidate(point_coefficients[degree], s)
                row.update({"point_label": label, "Z": point_z, "c": c_value, "radial_degree": degree})
                extreme_rows.append(row)
        for s in ("2", "4.1"):
            raw = {
                degree: evaluate_core_jets(
                    point_coefficients[degree],
                    _mp(s) / _mp(extreme_configuration["Lambda"]),
                )
                for degree in (18, 24)
            }
            extreme_convergence.append(
                {
                    "point_label": label,
                    "Z": point_z,
                    "c": c_value,
                    "s": s,
                    "relative_18_to_24": {
                        key: _nstr(_relative_unfloored(raw[24][key], raw[18][key]), 30)
                        for key in ("F", "F_R", "Uz", "Uz_R", "P")
                    },
                }
            )
    extreme_max_relative: dict[str, str] = {}
    for key in ("F", "F_R", "Uz", "Uz_R", "P"):
        extreme_max_relative[key] = _nstr(
            max(
                _mp(item["relative_18_to_24"][key])
                for item in extreme_convergence
            ),
            30,
        )

    return {
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "configuration": configuration,
        "axis_root": _nstr(axis.Z0, 60),
        "axis_root_width": _nstr(axis._root_width, 60),
        "pressure_metadata": pressure.metadata,
        "kappa_definition": "-( (2 R F_R)^2 + (sqrt(2R) Uz_R)^2 ) / (F * 2 R F_R)",
        "rows": rows,
        "degree_convergence": convergence,
        "lambda_sweep_same_physical_P0": lambda_sweep,
        "rounded_root_probe": rounded_root_probe,
        "lambda_1e36_dangerous_offsets": {
            "configuration": extreme_configuration,
            "pressure_reused_from_baseline": True,
            "axis_Uz_R_at_exact_root": _signed_log(root_jet.Uz_R),
            "abs_H0_Z_at_exact_root": _nstr(h0_z_abs, 40),
            "offset_unit_abs_U1_over_Lambda_abs_H0_Z": _nstr(offset_unit, 40),
            "points": extreme_rows,
            "degree_18_to_24_relative_max": extreme_max_relative,
            "scope": "True root, rounded literal, and requested dangerous offsets only; no global claim.",
        },
        "pressure_jet_uncertainty": {
            "source": "P0(0) from finite outer quadrature, extended as K/(1+Z^2)^2",
            "omitted_future_tail_bound": pressure.metadata["omitted_future_tail_bound_absolute"],
            "full_source_pressure_exact": False,
        },
        "scope": "Finite local radial Taylor prefix through s=4.1; no global core, cone, matching, or contraction claim.",
        "global_core_certified": False,
        "cone_certified": False,
        "matching_certified": False,
    }


def _json_value(value: Any) -> Any:
    if isinstance(value, mp.mpf):
        return _nstr(value, 50)
    if isinstance(value, dict):
        return {str(k): _json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(v) for v in value]
    return value


if __name__ == "__main__":
    report = run_experiment()
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(_json_value(report), indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "rows": len(report["rows"]),
                "degrees": [12, 18, 24],
                "global_core_certified": False,
            }
        ),
        flush=True,
    )


__all__ = [
    "PrototypeAxisPressure",
    "build_coefficients",
    "diagnose_candidate",
    "run_experiment",
]
