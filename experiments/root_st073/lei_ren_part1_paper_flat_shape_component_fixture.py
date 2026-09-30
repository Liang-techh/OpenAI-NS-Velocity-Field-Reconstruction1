"""Resolved independent check for the pressure/width flat-shape lift.

The production adapter composes scalar B derivatives with a rectangular
pressure/width jet.  This fixture independently integrates the jet-valued
integrands coefficient by coefficient over the resolved ``T = 400`` interval
and checks both the value and the Z chain-rule tangent.  It also verifies that
each Taylor term and each scalar saddle receipt remains available before the
public jet sum.
"""

from __future__ import annotations

import importlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Callable

import mpmath as mp

try:
    from .lei_ren_part1_paper_axial_primitive import _sigma_mp
    from .lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
    from .lei_ren_part1_paper_flat_shape_component import compose_flat_shape_defect
except (ImportError, ValueError):
    from lei_ren_part1_paper_axial_primitive import _sigma_mp
    from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
    from lei_ren_part1_paper_flat_shape_component import compose_flat_shape_defect


PRECISION = 100
T = mp.mpf(400)
with mp.workdps(320):
    K = mp.mpf("1.2")
    M = mp.mpf(2)
PRESSURE_ORDER = 2
WIDTH_ORDER = 1
CUTS = tuple(mp.mpf(value) for value in (0, 10, 25, 50, 100, 200, 300, 400))
COMPONENT_ERROR_THRESHOLD = mp.mpf("1e-12")
TERM_ZERO_FLOOR = mp.mpf("1e-250")


def _jet(atoms: Mapping[tuple[int, int], Any]) -> PressureWidthJet:
    return PressureWidthJet(
        atoms,
        pressure_order=PRESSURE_ORDER,
        width_order=WIDTH_ORDER,
    )


def _inputs() -> tuple[PressureWidthJet, PressureWidthJet]:
    # The scalar B0 is -4.  All retained nilpotent atoms are small but
    # nonzero, so pressure and width powers can be inspected separately.
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
    return B, B_Z


def _quad_jet(function: Callable[[mp.mpf], PressureWidthJet]) -> PressureWidthJet:
    atoms: dict[tuple[int, int], mp.mpf] = {}
    for pressure in range(PRESSURE_ORDER + 1):
        for width in range(WIDTH_ORDER + 1):
            atoms[(pressure, width)] = mp.fsum(
                mp.quad(
                    lambda ell, p=pressure, w=width: function(ell).component(p, w),
                    [left, right],
                )
                for left, right in zip(CUTS, CUTS[1:])
            )
    return _jet(atoms)


def _direct_integrals(
    B: PressureWidthJet,
    B_Z: PressureWidthJet,
) -> tuple[PressureWidthJet, PressureWidthJet]:
    def argument(ell: mp.mpf) -> PressureWidthJet:
        q = _sigma_mp(ell / T)
        return B * (M * q)

    def value_integrand(ell: mp.mpf) -> PressureWidthJet:
        return (argument(ell).exp() - 1) * mp.exp(-K * ell)

    def tangent_integrand(ell: mp.mpf) -> PressureWidthJet:
        q = _sigma_mp(ell / T)
        return (
            argument(ell).exp()
            * (M * q)
            * B_Z
            * mp.exp(-K * ell)
        )

    return _quad_jet(value_integrand), _quad_jet(tangent_integrand)


def _relative_errors(
    actual: PressureWidthJet,
    expected: PressureWidthJet,
) -> tuple[dict[str, str], mp.mpf]:
    errors: dict[str, str] = {}
    maximum = mp.mpf(0)
    for pressure in range(PRESSURE_ORDER + 1):
        for width in range(WIDTH_ORDER + 1):
            key = f"{pressure},{width}"
            left = actual.component(pressure, width)
            right = expected.component(pressure, width)
            scale = max(abs(right), TERM_ZERO_FLOOR)
            error = abs(left - right) / scale
            errors[key] = mp.nstr(error, 30)
            maximum = max(maximum, error)
    return errors, maximum


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
        B, B_Z = _inputs()
        direct_value, direct_tangent = _direct_integrals(B, B_Z)
        composed = compose_flat_shape_defect(
            K,
            M,
            B,
            B_Z,
            T,
            precision=precision,
            order=32,
            window=24,
        )
        value_errors, maximum_value_error = _relative_errors(
            composed["value"], direct_value
        )
        tangent_errors, maximum_tangent_error = _relative_errors(
            composed["tangent"], direct_tangent
        )

        expected_term_count = PRESSURE_ORDER + WIDTH_ORDER + 1
        value_terms = composed["value_terms"]
        tangent_terms = composed["tangent_terms"]
        derivative_cache = composed["derivative_cache"]
        terms_present = (
            sorted(value_terms) == list(range(expected_term_count))
            and sorted(tangent_terms) == list(range(expected_term_count))
            and sorted(derivative_cache)
            == list(range(expected_term_count + 1))
        )
        value_terms_nonzero = {
            str(n): bool(term.atoms) for n, term in value_terms.items()
        }
        tangent_terms_nonzero = {
            str(n): bool(term.atoms) for n, term in tangent_terms.items()
        }
        all_terms_nonzero = all(value_terms_nonzero.values()) and all(
            tangent_terms_nonzero.values()
        )
        mode_receipts = composed["term_receipts"]
        tangent_mode_receipts = composed["tangent_term_receipts"]
        receipts_present = all(
            any(name in mode_receipts[n] for name in ("mode", "derivative_mode"))
            for n in range(expected_term_count)
        ) and all(
            any(
                name in tangent_mode_receipts[n]
                for name in ("mode", "derivative_mode")
            )
            for n in range(expected_term_count)
        )

        if maximum_value_error >= COMPONENT_ERROR_THRESHOLD:
            raise AssertionError(
                f"value jet mismatch: {maximum_value_error}"
            )
        if maximum_tangent_error >= COMPONENT_ERROR_THRESHOLD:
            raise AssertionError(
                f"tangent jet mismatch: {maximum_tangent_error}"
            )
        if not terms_present:
            raise AssertionError("derivative or term receipts are incomplete")
        if not all_terms_nonzero:
            raise AssertionError("a retained Taylor term unexpectedly vanished")
        if not receipts_present:
            raise AssertionError("a retained Taylor term lacks its saddle receipt")

        report = {
            "passed": True,
            "precision": precision,
            "k": mp.nstr(K, 30),
            "m": mp.nstr(M, 30),
            "T": mp.nstr(T, 30),
            "cuts": [mp.nstr(value, 20) for value in CUTS],
            "B": _json_safe(B),
            "B_Z": _json_safe(B_Z),
            "B0": mp.nstr(composed["B0"], 30),
            "metadata": _json_safe(composed["metadata"]),
            "value": _json_safe(composed["value"]),
            "direct_value": _json_safe(direct_value),
            "tangent": _json_safe(composed["tangent"]),
            "direct_tangent": _json_safe(direct_tangent),
            "value_component_relative_errors": value_errors,
            "tangent_component_relative_errors": tangent_errors,
            "maximum_value_component_relative_error": mp.nstr(
                maximum_value_error, 30
            ),
            "maximum_tangent_component_relative_error": mp.nstr(
                maximum_tangent_error, 30
            ),
            "component_error_threshold": mp.nstr(COMPONENT_ERROR_THRESHOLD, 20),
            "value_terms": _json_safe(value_terms),
            "tangent_terms": _json_safe(tangent_terms),
            "value_terms_nonzero": value_terms_nonzero,
            "tangent_terms_nonzero": tangent_terms_nonzero,
            "derivative_orders_cached": sorted(derivative_cache),
            "terms_preserved_before_sum": terms_present
            and composed["terms_preserved_before_sum"],
            "mode_receipts_present": receipts_present,
            "term_receipts": _json_safe(mode_receipts),
            "tangent_term_receipts": _json_safe(tangent_mode_receipts),
            "rectangular_pressure_width_truncation": True,
            "pressure_parameter_distinct_from_temporal": True,
            "quadrature_remainder_enclosed": False,
            "all_five_moment_claim": False,
            "global_field_installed": False,
            "limitations": (
                "Independent coefficientwise MP quadrature checks this resolved "
                "finite ring only; saddle tails and quadrature remainder are unenclosed."
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
