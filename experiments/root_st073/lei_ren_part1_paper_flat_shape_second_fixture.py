"""Resolved second-Z check for the flat-shape pressure/width lift.

The component adapter is checked against an independent scalar evaluation of
the shifted polynomial

``B(Z) = B0 + B_Z(0) Z + B_ZZ(0) Z**2 / 2``.

The finite-Z stencil only checks the constant pressure/width atom.  The full
finite ring is checked separately by replaying the two second-chain pieces
and their retained Taylor terms before summation.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import mpmath as mp

try:
    from .lei_ren_part1_paper_flat_shape_component import compose_flat_shape_defect
    from .lei_ren_part1_paper_flat_shape_defect import evaluate_flat_shape_defect
    from .lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
except (ImportError, ValueError):
    from lei_ren_part1_paper_flat_shape_component import compose_flat_shape_defect
    from lei_ren_part1_paper_flat_shape_defect import evaluate_flat_shape_defect
    from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet


PRECISION = 120
T = mp.mpf(400)
with mp.workdps(240):
    K = mp.mpf("1.2")
    M = mp.mpf(2)
PRESSURE_ORDER = 2
WIDTH_ORDER = 1
STENCIL_STEP = mp.mpf("1e-3")
SECOND_ERROR_THRESHOLD = mp.mpf("1e-7")
ALGEBRA_ERROR_THRESHOLD = mp.mpf("1e-100")
TERM_ZERO_FLOOR = mp.mpf("1e-250")


def _jet(atoms: Mapping[tuple[int, int], Any]) -> PressureWidthJet:
    return PressureWidthJet(
        atoms,
        pressure_order=PRESSURE_ORDER,
        width_order=WIDTH_ORDER,
    )


def _inputs() -> tuple[
    PressureWidthJet,
    PressureWidthJet,
    PressureWidthJet,
]:
    B = _jet(
        {
            (0, 0): mp.mpf("-4"),
            (1, 0): mp.mpf("1e-3"),
            (0, 1): mp.mpf("2e-3"),
            (2, 0): mp.mpf("-3e-4"),
            (1, 1): mp.mpf("4e-4"),
            (2, 1): mp.mpf("2e-4"),
        }
    )
    B_Z = _jet(
        {
            (0, 0): mp.mpf(".7"),
            (1, 0): mp.mpf("-.25"),
            (0, 1): mp.mpf(".18"),
            (2, 0): mp.mpf(".05"),
            (1, 1): mp.mpf("-.03"),
            (2, 1): mp.mpf(".02"),
        }
    )
    B_ZZ = _jet(
        {
            (0, 0): mp.mpf(".3"),
            (1, 0): mp.mpf(".11"),
            (0, 1): mp.mpf("-.09"),
            (2, 0): mp.mpf(".04"),
            (1, 1): mp.mpf(".02"),
            (2, 1): mp.mpf(".01"),
        }
    )
    return B, B_Z, B_ZZ


def _scalar_value(record: Mapping[str, Any]) -> mp.mpf:
    value = record.get("value")
    if value is None:
        raise KeyError("scalar flat-shape record has no value")
    return mp.mpf(value)


def _shifted_scalar_value(
    B: PressureWidthJet,
    B_Z: PressureWidthJet,
    B_ZZ: PressureWidthJet,
    z: mp.mpf,
    *,
    precision: int,
) -> mp.mpf:
    shifted = (
        B.component(0, 0)
        + B_Z.component(0, 0) * z
        + B_ZZ.component(0, 0) * z * z / 2
    )
    return _scalar_value(
        evaluate_flat_shape_defect(
            K,
            M,
            shifted,
            T,
            precision=precision,
            order=32,
            window=24,
        )
    )


def _relative_error(actual: mp.mpf, expected: mp.mpf) -> mp.mpf:
    return abs(actual - expected) / max(abs(actual), abs(expected), TERM_ZERO_FLOOR)


def _max_atom_error(left: PressureWidthJet, right: PressureWidthJet) -> mp.mpf:
    return max(
        (
            abs(left.component(pressure, width) - right.component(pressure, width))
            for pressure in range(PRESSURE_ORDER + 1)
            for width in range(WIDTH_ORDER + 1)
        ),
        default=mp.mpf(0),
    )


def _json_safe(value: Any) -> Any:
    if isinstance(value, PressureWidthJet):
        return value.as_dict()
    if isinstance(value, mp.mpf):
        return mp.nstr(value, 40)
    if isinstance(value, mp.mpc):
        return mp.nstr(value, 40)
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_safe(item) for item in value]
    return value


def run_fixture(*, precision: int = PRECISION) -> dict[str, Any]:
    precision = max(80, int(precision))
    with mp.workdps(precision):
        B, B_Z, B_ZZ = _inputs()
        composed = compose_flat_shape_defect(
            K,
            M,
            B,
            B_Z,
            T,
            B_ZZ=B_ZZ,
            precision=precision,
            order=32,
            window=24,
        )
        baseline = compose_flat_shape_defect(
            K,
            M,
            B,
            B_Z,
            T,
            precision=precision,
            order=32,
            window=24,
        )

        h = STENCIL_STEP
        samples = {
            offset: _shifted_scalar_value(
                B,
                B_Z,
                B_ZZ,
                offset * h,
                precision=precision,
            )
            for offset in (-1, 0, 1)
        }
        finite_z_second = (samples[-1] - 2 * samples[0] + samples[1]) / (h * h)
        fullsecond_constant = composed["fullsecond"].component(0, 0)
        finite_z_relative_error = _relative_error(
            finite_z_second,
            fullsecond_constant,
        )

        derivative_values = composed["derivative_values"]
        analytic_constant = (
            derivative_values[2] * B_Z.component(0, 0) ** 2
            + derivative_values[1] * B_ZZ.component(0, 0)
        )
        analytic_relative_error = _relative_error(
            fullsecond_constant,
            analytic_constant,
        )

        chain_terms = composed["second_chain_terms"]
        curvature_terms = composed["second_curvature_terms"]
        second_terms = composed["second_terms"]
        split_errors = {
            str(n): _max_atom_error(
                second_terms[n], chain_terms[n] + curvature_terms[n]
            )
            for n in second_terms
        }
        second_sum = composed["fullsecond"]
        replayed_sum = PressureWidthJet(
            0,
            pressure_order=PRESSURE_ORDER,
            width_order=WIDTH_ORDER,
        )
        for term in second_terms.values():
            replayed_sum = replayed_sum + term
        sum_replay_error = _max_atom_error(second_sum, replayed_sum)

        expected_cache = list(range(PRESSURE_ORDER + WIDTH_ORDER + 3))
        cache_complete = sorted(composed["derivative_cache"]) == expected_cache
        old_cache_complete = sorted(baseline["derivative_cache"]) == expected_cache[:-1]
        old_api_preserved = (
            "second_terms" not in baseline
            and "fullsecond" not in baseline
            and baseline["derivative_truncation_order"]
            == PRESSURE_ORDER + WIDTH_ORDER + 1
        )

        if finite_z_relative_error >= SECOND_ERROR_THRESHOLD:
            raise AssertionError(
                f"finite-Z second mismatch: {finite_z_relative_error}"
            )
        if analytic_relative_error >= ALGEBRA_ERROR_THRESHOLD:
            raise AssertionError(
                f"constant second chain mismatch: {analytic_relative_error}"
            )
        if max(split_errors.values(), default=mp.mpf(0)) >= ALGEBRA_ERROR_THRESHOLD:
            raise AssertionError("second chain pieces do not replay their sum")
        if sum_replay_error >= ALGEBRA_ERROR_THRESHOLD:
            raise AssertionError("second Taylor terms do not replay fullsecond")
        if not cache_complete or not old_cache_complete:
            raise AssertionError("derivative caches do not retain the expected orders")
        if not old_api_preserved:
            raise AssertionError("B_ZZ=None changed the first-order API")
        if composed["metadata"]["second_remainder_enclosed"] is not False:
            raise AssertionError("second remainder must remain unenclosed")

        report = {
            "passed": True,
            "precision": precision,
            "k": mp.nstr(K, 30),
            "m": mp.nstr(M, 30),
            "T": mp.nstr(T, 30),
            "B": _json_safe(B),
            "B_Z": _json_safe(B_Z),
            "B_ZZ": _json_safe(B_ZZ),
            "B0": mp.nstr(composed["B0"], 30),
            "finite_Z_stencil_step": mp.nstr(h, 20),
            "finite_Z_second": mp.nstr(finite_z_second, 40),
            "fullsecond_constant": mp.nstr(fullsecond_constant, 40),
            "analytic_constant": mp.nstr(analytic_constant, 40),
            "finite_Z_relative_error": mp.nstr(finite_z_relative_error, 30),
            "analytic_relative_error": mp.nstr(analytic_relative_error, 30),
            "second_term_split_errors": {
                key: mp.nstr(value, 30) for key, value in split_errors.items()
            },
            "second_sum_replay_error": mp.nstr(sum_replay_error, 30),
            "derivative_orders_cached": sorted(composed["derivative_cache"]),
            "old_derivative_orders_cached": sorted(baseline["derivative_cache"]),
            "metadata": _json_safe(composed["metadata"]),
            "second_terms": _json_safe(second_terms),
            "second_chain_terms": _json_safe(chain_terms),
            "second_curvature_terms": _json_safe(curvature_terms),
            "fullsecond": _json_safe(composed["fullsecond"]),
            "old_api_preserved": old_api_preserved,
            "quadrature_remainder_enclosed": False,
            "second_remainder_enclosed": False,
            "global_field_installed": False,
            "limitations": (
                "Resolved scalar finite-Z and finite-ring algebra checks only; "
                "saddle tails, quadrature remainder and global field claims remain unenclosed."
            ),
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
