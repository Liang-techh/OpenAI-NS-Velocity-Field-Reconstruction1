"""Independent resolved check for flat reshape defect integrals.

The check compares the new signed-log saddle evaluator against direct
high-precision quadrature of

    I(k,m,B) = integral exp(-k ell) expm1(m B sigma(ell/T)) d ell,

and checks ``dI/dB`` against the differentiated integrand.  It also records
the observed positive integrand mode and compares its curvature with the
large-``T`` model ``-k ell - T^2/ell^2``.  These are local resolved checks;
they do not enclose a pressure/width remainder or certify a global field.
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from typing import Any, Mapping

import mpmath as mp

try:
    from .lei_ren_part1_paper_axial_primitive import _sigma_mp
except (ImportError, ValueError):
    from lei_ren_part1_paper_axial_primitive import _sigma_mp


PRECISION = 100
T = mp.mpf(400)
B_VALUES = (mp.mpf(-4), mp.mpf(4))
with mp.workdps(320):
    # Construct decimal inputs before the default mpmath context rounds them;
    # in particular, ``mp.mpf(".2")`` at 15 digits is not the requested 0.2.
    CASES = (
        (mp.mpf("1.6"), mp.mpf(1)),
        (mp.mpf(".2"), mp.mpf(2)),
        (mp.mpf("1.2"), mp.mpf(2)),
    )
CUTS = tuple(mp.mpf(value) for value in (0, 10, 25, 50, 100, 200, 300, 400))
# The production evaluator is a finite saddle-window rule.  The resolved
# comparison is still much tighter than the downstream source tolerances, but
# a 1e-14 log/relative threshold avoids pretending that this unenclosed rule
# is an exact quadrature certificate.
ERROR_THRESHOLD = mp.mpf("1e-14")
MODE_RATIO_THRESHOLD = mp.mpf(".35")
CURVATURE_MODEL_DISCREPANCY_THRESHOLD = mp.mpf(".02")


def _quad_split(function: Any) -> mp.mpf:
    return mp.fsum(
        mp.quad(function, [left, right])
        for left, right in zip(CUTS, CUTS[1:])
    )


def _direct(k: mp.mpf, m: mp.mpf, B: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    def sigma(ell: mp.mpf) -> mp.mpf:
        return _sigma_mp(ell / T)

    value = _quad_split(
        lambda ell: mp.exp(-k * ell) * mp.expm1(m * B * sigma(ell))
    )
    derivative = _quad_split(
        lambda ell: mp.exp(-k * ell)
        * (m * sigma(ell))
        * mp.exp(m * B * sigma(ell))
    )
    return value, derivative


def _log_integrand(ell: mp.mpf, k: mp.mpf, m: mp.mpf, B: mp.mpf) -> mp.mpf:
    if ell <= 0:
        return mp.ninf
    sigma = _sigma_mp(ell / T)
    defect = abs(mp.expm1(m * B * sigma))
    if defect == 0:
        return mp.ninf
    return -k * ell + mp.log(defect)


def _observed_mode(k: mp.mpf, m: mp.mpf, B: mp.mpf) -> tuple[mp.mpf, mp.mpf, bool]:
    """Find the positive mode of the absolute integrand on the full interval."""

    grid = [CUTS[0] + (CUTS[-1] - CUTS[0]) * i / 800 for i in range(801)]
    values = [_log_integrand(x, k, m, B) for x in grid]
    index = max(range(1, len(values) - 1), key=lambda i: values[i])
    boundary = index in (1, len(values) - 2)
    left, right = grid[index - 1], grid[index + 1]
    # Golden-section maximization of log magnitude inside the selected cell.
    phi = (mp.sqrt(5) - 1) / 2
    x1 = right - phi * (right - left)
    x2 = left + phi * (right - left)
    f1, f2 = _log_integrand(x1, k, m, B), _log_integrand(x2, k, m, B)
    for _ in range(100):
        if f1 < f2:
            left = x1
            x1, f1 = x2, f2
            x2 = left + phi * (right - left)
            f2 = _log_integrand(x2, k, m, B)
        else:
            right = x2
            x2, f2 = x1, f1
            x1 = right - phi * (right - left)
            f1 = _log_integrand(x1, k, m, B)
    mode = (left + right) / 2
    return mode, _log_integrand(mode, k, m, B), boundary


def _asymptotic_mode(k: mp.mpf) -> mp.mpf:
    return (2 * T * T / k) ** (mp.mpf(1) / 3)


def _extract(result: Any, names: tuple[str, ...]) -> Any:
    if isinstance(result, Mapping):
        for name in names:
            if name in result:
                return result[name]
    for name in names:
        if hasattr(result, name):
            return getattr(result, name)
    raise KeyError(f"missing evaluator field {names}")


def _signed_log(value: Any) -> tuple[int, mp.mpf]:
    """Read a signed-log record or a numeric value."""

    if isinstance(value, Mapping):
        sign = value.get("sign", value.get("derivative_sign"))
        log_abs = value.get(
            "log_abs",
            value.get(
                "derivative_log_abs",
                value.get("logI", value.get("log_abs_value")),
            ),
        )
        if sign is not None and log_abs is not None:
            return int(sign), mp.mpf(str(log_abs))
    sign = getattr(value, "sign", None)
    log_abs = getattr(value, "log_abs", getattr(value, "logI", None))
    if sign is not None and log_abs is not None:
        return int(sign), mp.mpf(str(log_abs))
    scalar = mp.mpf(str(value))
    if scalar == 0:
        return 0, mp.ninf
    return (1 if scalar > 0 else -1), mp.log(abs(scalar))


def _extract_evaluation(module: Any, T_value: mp.mpf, k: mp.mpf, m: mp.mpf, B: mp.mpf) -> Any:
    """Support the worker's evaluator class/function without changing it."""

    for name in ("evaluate_flat_shape_defect", "evaluate_defect", "integrate_defect"):
        function = getattr(module, name, None)
        if callable(function):
            for args in ((k, m, B, T_value), (T_value, k, m, B)):
                try:
                    return function(*args)
                except TypeError:
                    pass
    for name in ("FlatShapeDefect", "FlatShapeDefectIntegral", "FlatShapeDefectEvaluator"):
        cls = getattr(module, name, None)
        if cls is None:
            continue
        for kwargs in (
            {"T": T_value, "k": k, "m": m, "B": B},
            {"T": T_value, "kappa": k, "m": m, "B": B},
        ):
            try:
                instance = cls(**kwargs)
                for method in ("evaluate", "compute", "value"):
                    candidate = getattr(instance, method, None)
                    if callable(candidate):
                        return candidate()
                    if candidate is not None:
                        return candidate
            except TypeError:
                continue
    raise AttributeError("flat-shape defect evaluator API not found")


def _evaluate(module: Any, T_value: mp.mpf, k: mp.mpf, m: mp.mpf, B: mp.mpf) -> tuple[int, mp.mpf, mp.mpf, Any]:
    result = _extract_evaluation(module, T_value, k, m, B)
    if isinstance(result, Mapping) and "sign" in result and "log_abs" in result:
        sign = int(result["sign"])
        log_abs = mp.mpf(str(result["log_abs"]))
    else:
        value = _extract(result, ("I", "integral", "value", "signed_log"))
        sign, log_abs = _signed_log(value)
    derivative_value = _extract(
        result,
        ("dI_dB", "derivative", "d_integral_dB", "dI_DB", "derivative_value"),
    )
    if isinstance(derivative_value, Mapping):
        derivative = _signed_log(derivative_value)[0] * mp.exp(
            _signed_log(derivative_value)[1]
        )
    else:
        derivative = mp.mpf(str(derivative_value))
    return sign, log_abs, derivative, result


def run_fixture(*, precision: int = PRECISION) -> dict[str, Any]:
    precision = max(80, int(precision))
    with mp.workdps(precision):
        module = importlib.import_module("lei_ren_part1_paper_flat_shape_defect")
        rows: list[dict[str, Any]] = []
        maximum_log_error = mp.mpf(0)
        maximum_derivative_error = mp.mpf(0)
        mode_rows: list[dict[str, Any]] = []
        for k, m in CASES:
            asymptotic = _asymptotic_mode(k)
            for B in B_VALUES:
                direct, direct_derivative = _direct(k, m, B)
                expected_sign = 1 if direct > 0 else -1
                expected_log = mp.log(abs(direct))
                sign, log_value, derivative, raw = _evaluate(module, T, k, m, B)
                log_error = abs(log_value - expected_log)
                derivative_error = abs(derivative - direct_derivative) / max(
                    mp.mpf(1), abs(derivative), abs(direct_derivative)
                )
                maximum_log_error = max(maximum_log_error, log_error)
                maximum_derivative_error = max(maximum_derivative_error, derivative_error)
                mode, mode_log, boundary = _observed_mode(k, m, B)
                curvature = mp.diff(
                    lambda ell: _log_integrand(ell, k, m, B), mode, 2
                )
                curvature_model = -6 * T * T / mode**4
                mode_ratio = mode / asymptotic
                mode_rows.append(
                    {
                        "k": mp.nstr(k, 20),
                        "m": mp.nstr(m, 20),
                        "B": mp.nstr(B, 20),
                        "mode": mp.nstr(mode, 30),
                        "asymptotic_mode": mp.nstr(asymptotic, 30),
                        "mode_ratio": mp.nstr(mode_ratio, 30),
                        "mode_log_abs_integrand": mp.nstr(mode_log, 30),
                        "positive_mode": bool(mode > 0),
                        "boundary_mode": bool(boundary),
                        "curvature": mp.nstr(curvature, 30),
                        "curvature_model": mp.nstr(curvature_model, 30),
                        "curvature_model_scaled_discrepancy": mp.nstr(
                            abs(curvature - curvature_model)
                            / max(mp.mpf(1), abs(curvature), abs(curvature_model)),
                            30,
                        ),
                    }
                )
                rows.append(
                    {
                        "k": mp.nstr(k, 20),
                        "m": mp.nstr(m, 20),
                        "B": mp.nstr(B, 20),
                        "expected_sign": expected_sign,
                        "actual_sign": sign,
                        "expected_log_abs": mp.nstr(expected_log, 30),
                        "actual_log_abs": mp.nstr(log_value, 30),
                        "log_error": mp.nstr(log_error, 30),
                        "expected_dI_dB": mp.nstr(direct_derivative, 30),
                        "actual_dI_dB": mp.nstr(derivative, 30),
                        "relative_dI_dB_error": mp.nstr(derivative_error, 30),
                        "nonzero_defect": bool(direct != 0),
                        "nonzero_dI_dB": bool(direct_derivative != 0),
                        "signed_log_record_present": (
                            isinstance(raw, Mapping)
                            and "sign" in raw
                            and "log_abs" in raw
                        ),
                    }
                )

        max_mode_error = max(
            abs(mp.mpf(row["mode_ratio"]) - 1) for row in mode_rows
        )
        if any(row["expected_sign"] != row["actual_sign"] for row in rows):
            raise AssertionError("signed defect sign mismatch")
        if maximum_log_error >= ERROR_THRESHOLD:
            key = max(rows, key=lambda row: mp.mpf(row["log_error"]))
            raise AssertionError(
                f"signed-log integral mismatch {key['k']},{key['m']},{key['B']}: {maximum_log_error}"
            )
        if maximum_derivative_error >= ERROR_THRESHOLD:
            key = max(rows, key=lambda row: mp.mpf(row["relative_dI_dB_error"]))
            raise AssertionError(
                f"dI/dB mismatch {key['k']},{key['m']},{key['B']}: {maximum_derivative_error}"
            )
        if any(not row["positive_mode"] or row["boundary_mode"] for row in mode_rows):
            raise AssertionError("a defect integrand lacks an interior positive mode")
        if max_mode_error >= MODE_RATIO_THRESHOLD:
            raise AssertionError(
                f"observed mode is outside the asymptotic saddle band: {max_mode_error}"
            )
        maximum_curvature_error = max(
            mp.mpf(row["curvature_model_scaled_discrepancy"]) for row in mode_rows
        )
        if maximum_curvature_error >= CURVATURE_MODEL_DISCREPANCY_THRESHOLD:
            raise AssertionError(
                f"finite-T curvature disagrees with the saddle model: {maximum_curvature_error}"
            )
        report = {
            "passed": True,
            "precision": precision,
            "T": mp.nstr(T, 30),
            "B_values": [mp.nstr(value, 20) for value in B_VALUES],
            "cases": [[mp.nstr(k, 20), mp.nstr(m, 20)] for k, m in CASES],
            "cuts": [mp.nstr(value, 20) for value in CUTS],
            "maximum_log_abs_integral_error": mp.nstr(maximum_log_error, 30),
            "maximum_relative_dI_dB_error": mp.nstr(maximum_derivative_error, 30),
            "comparison_error_threshold": mp.nstr(ERROR_THRESHOLD, 20),
            "maximum_mode_ratio_deviation": mp.nstr(max_mode_error, 30),
            "maximum_curvature_model_scaled_discrepancy": mp.nstr(
                maximum_curvature_error, 30
            ),
            "rows": rows,
            "modes": mode_rows,
            "asymptotic_mode_formula": "ell* = (2 T^2 / k)^(1/3)",
            "curvature_formula": "phi''(ell*) = -6 T^2 / ell*^4 for phi=-k ell-T^2/ell^2",
            "independent_quadrature": "mp.quad over cuts 0,10,25,50,100,200,300,400",
            "dI_dB_integrand": "exp(-k ell) m sigma exp(m B sigma)",
            "signed_log_checked": True,
            "positive_modes_checked": True,
            "asymptotic_mode_checked": True,
            "asymptotic_curvature_model_compared": True,
            "quadrature_error_enclosed": False,
            "saddle_remainder_enclosed": False,
            "global_field_installed": False,
        }
        Path(__file__).with_suffix(".json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return report


def fixture(**kwargs: Any) -> dict[str, Any]:
    return run_fixture(**kwargs)


def main() -> None:
    print(json.dumps(run_fixture(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["fixture", "run_fixture"]
