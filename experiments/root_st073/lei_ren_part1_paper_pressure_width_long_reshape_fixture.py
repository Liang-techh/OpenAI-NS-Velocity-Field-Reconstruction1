"""Independent resolved check for the finite pressure-width long reshape.

The production reshape is evaluated on a compact analytic endpoint provider
whose fields and first ``Z`` tangents are carried by ``AxialDual`` over the
pressure/width jet.  The fixture compares its normalized angular, pressure,
and swirl quadratures with independent high-precision scalar quadrature and
checks the first tangents with a separate fourth-order centered stencil.
Endpoint moment and raw quadratic seeds remain visible in the receipt, even
when their sums are dominated by the corresponding reshape increments.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import mpmath as mp

try:
    from .lei_ren_part1_paper_axial_dual import AxialDual
    from .lei_ren_part1_paper_axial_primitive import _sigma_mp
    from .lei_ren_part1_paper_pressure_width_axial_comparison import (
        AxialPressureWidthComparison,
    )
    from .lei_ren_part1_paper_pressure_width_axial_bridge_fixture import (
        _bundle,
        _make_formal_core,
    )
    from .lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
    from .lei_ren_part1_paper_pressure_width_long_reshape import (
        PressureWidthLongReshape,
    )
except (ImportError, ValueError):
    from lei_ren_part1_paper_axial_dual import AxialDual
    from lei_ren_part1_paper_axial_primitive import _sigma_mp
    from lei_ren_part1_paper_pressure_width_axial_comparison import (
        AxialPressureWidthComparison,
    )
    from lei_ren_part1_paper_pressure_width_axial_bridge_fixture import (
        _bundle,
        _make_formal_core,
    )
    from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
    from lei_ren_part1_paper_pressure_width_long_reshape import (
        PressureWidthLongReshape,
    )


PRECISION = 100
PRESSURE_ORDER = 3
WIDTH_ORDER = 2
RESHAPE_A = mp.mpf(1)
RESHAPE_T = 400
LOG_C = mp.mpf(-4)
QUADRATURE_ORDER = 24
TAIL_DIGITS = 50
Z0 = mp.mpf(".3")
PHASES = (mp.mpf(".5"), mp.mpf(1))
FINITE_DIFFERENCE_STEP = mp.mpf("1e-4")
ERROR_THRESHOLD = mp.mpf("1e-8")
SEED_THRESHOLD = mp.mpf("1e-20")

MOMENT_NAMES = ("theta", "z", "theta_z", "z_theta", "p")
NORMALIZED_NAMES = ("theta", "pressure", "swirl")


def _scalar(value: Any) -> mp.mpf:
    """Materialize an axial dual or pressure/width jet at (1, 1)."""

    if isinstance(value, AxialDual):
        value = value.value
    if isinstance(value, PressureWidthJet):
        return value.evaluate(pressure=1, width=1)
    return mp.mpf(str(value))


def _scaled(actual: Any, expected: Any) -> mp.mpf:
    a = _scalar(actual)
    b = mp.mpf(str(expected)) if not isinstance(expected, (AxialDual, PressureWidthJet)) else _scalar(expected)
    return abs(a - b) / max(mp.mpf(1), abs(a), abs(b))


def _jet(
    comparison: AxialPressureWidthComparison,
    value: mp.mpf,
    *,
    pressure_atom: mp.mpf = mp.mpf("1e-8"),
    width_atom: mp.mpf = mp.mpf("2e-8"),
    mixed_atom: mp.mpf | None = None,
) -> PressureWidthJet:
    """Make a small nontrivial resolved pressure/width endpoint jet."""

    if mixed_atom is None:
        mixed_atom = pressure_atom * width_atom
    return PressureWidthJet(
        {
            (0, 0): value,
            (1, 0): pressure_atom,
            (0, 1): width_atom,
            (1, 1): mixed_atom,
        },
        pressure_order=PRESSURE_ORDER,
        width_order=WIDTH_ORDER,
    )


def _linear_dual(
    comparison: AxialPressureWidthComparison,
    z: mp.mpf,
    base: mp.mpf,
    slope: mp.mpf,
    *,
    pressure_atom: mp.mpf = mp.mpf("1e-8"),
    width_atom: mp.mpf = mp.mpf("2e-8"),
    tangent_pressure_atom: mp.mpf = mp.mpf(0),
    tangent_width_atom: mp.mpf = mp.mpf(0),
    tangent_mixed_atom: mp.mpf | None = None,
) -> AxialDual:
    value = _jet(
        comparison,
        base + slope * z,
        pressure_atom=pressure_atom,
        width_atom=width_atom,
    )
    tangent = _jet(
        comparison,
        slope,
        pressure_atom=tangent_pressure_atom,
        width_atom=tangent_width_atom,
        mixed_atom=tangent_mixed_atom,
    )
    return comparison._dual(value, tangent)


class _ResolvedR110Provider:
    """Small analytic endpoint provider implementing the reshape duck type."""

    def __init__(self, comparison: AxialPressureWidthComparison):
        self.comparison = comparison
        self.precision = comparison.precision

    def evaluate_R(self, radius: Any, Z: Any) -> dict[str, Any]:
        with mp.workdps(self.precision):
            if mp.mpf(str(radius)) != 110:
                raise ValueError("resolved fixture provider only exposes R=110")
            z = mp.mpf(str(Z))
            F = _linear_dual(self.comparison, z, mp.mpf(2), mp.mpf(".2"))
            V = _linear_dual(self.comparison, z, mp.mpf(".4"), mp.mpf(".1"))
            P = _linear_dual(self.comparison, z, mp.mpf(3), mp.mpf(".3"))
            moments = {
                name: _linear_dual(
                    self.comparison,
                    z,
                    mp.mpf(base),
                    mp.mpf(slope),
                )
                for name, base, slope in (
                    ("theta", "2", ".2"),
                    ("z", "3", ".3"),
                    ("theta_z", "4", ".4"),
                    ("z_theta", "5", ".5"),
                    ("p", "6", ".6"),
                )
            }
            axial = _linear_dual(
                self.comparison, z, mp.mpf("7") + mp.mpf(".7") * z, mp.mpf(".7")
            )
            # The swirl seed follows the explicit endpoint reference profile
            # d/dZ log(raw_swirl) = -2 Z/(1+Z^2).
            profile = -2 * z / (1 + z * z)
            swirl_value = mp.mpf(8) / (1 + z * z)
            swirl = _linear_dual(
                self.comparison,
                z,
                swirl_value - swirl_value * profile * z,
                swirl_value * profile,
                tangent_pressure_atom=mp.mpf("1e-8") * profile,
                tangent_width_atom=mp.mpf("2e-8") * profile,
                tangent_mixed_atom=mp.mpf("2e-16") * profile,
            )
            raw = {"axial": axial.value, "swirl": swirl.value}
            raw_Z = {"axial": axial.tangent, "swirl": swirl.tangent}
            zero = self.comparison._base_jet(0)
            return {
                "R": self.comparison._base_jet(110),
                "F": F.value,
                "FZ": F.tangent,
                "F_Z": F.tangent,
                "Uz": V.value,
                "UZ": V.tangent,
                "Uz_Z": V.tangent,
                "P0": (P - moments["p"]).value,
                "P0_Z": (P - moments["p"]).tangent,
                "P": P.value,
                "PZ": P.tangent,
                "P_Z": P.tangent,
                "Ur": zero,
                "moments": {name: value.value for name, value in moments.items()},
                "moments_Z": {name: value.tangent for name, value in moments.items()},
                "momentsZ": {name: value.tangent for name, value in moments.items()},
                "raw_quadratic_integrals": raw,
                "raw_quadratic_integrals_Z": raw_Z,
                "F_R": zero,
                "Uz_R": zero,
                "stress": {"U_r": zero},
                "region": "resolved_R110_endpoint",
            }


def _sigma_prime(s: mp.mpf, sigma: mp.mpf) -> mp.mpf:
    if 0 < s < 1:
        return sigma * (1 - sigma) * (2 / s**3 + 2 / (1 - s) ** 3)
    return mp.mpf(0)


def _quad_split(function: Any, endpoint: mp.mpf) -> mp.mpf:
    """Independent adaptive quadrature with fixed resolved subintervals."""

    if endpoint == 0:
        return mp.mpf(0)
    cuts = [mp.mpf(0)]
    for value in (1, 2, 4, 8, 16, 32, 64, 128, 256):
        value = mp.mpf(value)
        if value < endpoint:
            cuts.append(value)
    cuts.append(endpoint)
    return mp.fsum(mp.quad(function, [left, right]) for left, right in zip(cuts, cuts[1:]))


def _scalar_integrals(
    y: mp.mpf,
    z: mp.mpf,
    *,
    T: mp.mpf,
    f0: mp.mpf,
    log_c: mp.mpf,
) -> dict[str, mp.mpf]:
    """Direct scalar mp.quad replay of the Section 9.30 integrands."""

    root220 = mp.sqrt(220)
    u1 = root220 * f0
    B = mp.log(u1) + log_c + mp.log(1 + z * z)
    sigma_y = _sigma_mp(y / T)

    def ratio(ell: mp.mpf) -> mp.mpf:
        sigma_previous = _sigma_mp((y - ell) / T)
        return mp.exp(-ell / 10 - (sigma_previous - sigma_y) * B)

    return {
        "theta": _quad_split(
            lambda ell: mp.exp(-mp.mpf("1.5") * ell) * ratio(ell), y
        ),
        "pressure": _quad_split(lambda ell: ratio(ell) ** 2, y),
        "swirl": _quad_split(
            lambda ell: mp.exp(-ell) * ratio(ell) ** 2, y
        ),
    }


def _fourth_difference(values: dict[int, mp.mpf], step: mp.mpf) -> mp.mpf:
    return (-values[2] + 8 * values[1] - 8 * values[-1] + values[-2]) / (12 * step)


def _build_comparison(precision: int) -> AxialPressureWidthComparison:
    component, _formal = _make_formal_core(z0=Z0, precision=precision)
    return AxialPressureWidthComparison(
        _bundle(component.axis, component),
        h_b=mp.mpf("1e-5"),
        pressure_order=PRESSURE_ORDER,
        width_order=WIDTH_ORDER,
        transition_steps=4,
    )


def run_fixture(
    *,
    precision: int = PRECISION,
    phases: tuple[Any, ...] = PHASES,
    finite_difference_step: Any = FINITE_DIFFERENCE_STEP,
) -> dict[str, Any]:
    precision = max(80, int(precision))
    step = mp.mpf(str(finite_difference_step))
    if step <= 0:
        raise ValueError("finite_difference_step must be positive")
    with mp.workdps(precision):
        comparison = _build_comparison(precision)
        switches = _ResolvedR110Provider(comparison)
        reshape = PressureWidthLongReshape(
            switches,
            A=RESHAPE_A,
            logC=LOG_C,
            order=QUADRATURE_ORDER,
            tail_digits=TAIL_DIGITS,
        )

        endpoint = switches.evaluate_R(110, Z0)
        phase_zero = reshape.evaluate_phase(0, Z0)
        endpoint_errors: dict[str, mp.mpf] = {
            name: _scaled(phase_zero[name], endpoint[name])
            for name in ("F", "Uz", "P")
        }
        endpoint_errors.update(
            {
                f"moment:{name}": _scaled(
                    phase_zero["moments"][name], endpoint["moments"][name]
                )
                for name in MOMENT_NAMES
            }
        )
        endpoint_errors.update(
            {
                f"raw:{name}": _scaled(
                    phase_zero["raw_quadratic_seeds"][name],
                    endpoint["raw_quadratic_integrals"][name],
                )
                for name in ("axial", "swirl")
            }
        )

        raw_swirl = _scalar(phase_zero["raw_quadratic_seeds"]["swirl"])
        raw_swirl_Z = _scalar(phase_zero["raw_quadratic_seeds_Z"]["swirl"])
        swirl_profile = -2 * Z0 / (1 + Z0 * Z0)
        swirl_profile_error = abs(raw_swirl_Z / raw_swirl - swirl_profile)

        integral_errors: dict[str, mp.mpf] = {}
        tangent_errors: dict[str, mp.mpf] = {}
        phase_rows: list[dict[str, Any]] = []
        all_seed_values: list[mp.mpf] = []
        all_increment_values: list[mp.mpf] = []

        for raw_phase in phases:
            phase = mp.mpf(str(raw_phase))
            output = reshape.evaluate_phase(phase, Z0)
            y = phase * mp.mpf(RESHAPE_T)
            scalar_f0 = _scalar(endpoint["F"])
            expected = _scalar_integrals(
                y,
                Z0,
                T=mp.mpf(RESHAPE_T),
                f0=scalar_f0,
                log_c=LOG_C,
            )
            row_errors = {
                name: _scaled(output["normalized_integrals"][name], expected[name])
                for name in NORMALIZED_NAMES
            }
            integral_errors.update(
                {f"{name}@phase={mp.nstr(phase, 4)}": value for name, value in row_errors.items()}
            )

            tangent_expected: dict[str, mp.mpf] = {}
            fd_values: dict[str, dict[int, mp.mpf]] = {
                name: {} for name in NORMALIZED_NAMES
            }
            for index in (-2, -1, 1, 2):
                z_shift = Z0 + index * step
                endpoint_shift = switches.evaluate_R(110, z_shift)
                f_shift = _scalar(endpoint_shift["F"])
                integrals_shift = _scalar_integrals(
                    y,
                    z_shift,
                    T=mp.mpf(RESHAPE_T),
                    f0=f_shift,
                    log_c=LOG_C,
                )
                for name in NORMALIZED_NAMES:
                    fd_values[name][index] = integrals_shift[name]
            for name in NORMALIZED_NAMES:
                tangent_expected[name] = _fourth_difference(fd_values[name], step)
                tangent_actual = _scalar(output["normalized_integrals_Z"][name])
                tangent_errors[f"{name}@phase={mp.nstr(phase, 4)}"] = (
                    abs(tangent_actual - tangent_expected[name])
                    / max(mp.mpf(1), abs(tangent_actual), abs(tangent_expected[name]))
                )

            seeds = output["moment_seeds"]
            increments = output["moment_increments"]
            all_seed_values.extend(_scalar(seeds[name]) for name in MOMENT_NAMES)
            all_increment_values.extend(
                _scalar(increments[name]) for name in MOMENT_NAMES
            )
            all_seed_values.extend(
                _scalar(output["raw_quadratic_seeds"][name])
                for name in ("axial", "swirl")
            )
            all_increment_values.extend(
                _scalar(output["raw_quadratic_increments"][name])
                for name in ("axial", "swirl")
            )
            phase_rows.append(
                {
                    "phase": mp.nstr(phase, 12),
                    "integral_errors": {
                        name: mp.nstr(value, 24)
                        for name, value in row_errors.items()
                    },
                    "tangent_errors": {
                        name: mp.nstr(
                            tangent_errors[f"{name}@phase={mp.nstr(phase, 4)}"], 24
                        )
                        for name in NORMALIZED_NAMES
                    },
                    "normalized_integrals": {
                        name: mp.nstr(_scalar(output["normalized_integrals"][name]), 24)
                        for name in NORMALIZED_NAMES
                    },
                    "normalized_integrals_Z": {
                        name: mp.nstr(_scalar(output["normalized_integrals_Z"][name]), 24)
                        for name in NORMALIZED_NAMES
                    },
                    "axial_field_constant": _scaled(output["Uz"], endpoint["Uz"])
                    < ERROR_THRESHOLD,
                    "nonzero_moment_seeds": all(
                        abs(_scalar(seeds[name])) > SEED_THRESHOLD for name in MOMENT_NAMES
                    ),
                    "nonzero_moment_increments": all(
                        abs(_scalar(increments[name])) > SEED_THRESHOLD
                        for name in MOMENT_NAMES
                    ),
                }
            )

        max_integral = max(integral_errors.values(), default=mp.mpf(0))
        max_tangent = max(tangent_errors.values(), default=mp.mpf(0))
        max_endpoint = max(endpoint_errors.values(), default=mp.mpf(0))
        if max_endpoint >= ERROR_THRESHOLD:
            raise AssertionError(f"phase-zero endpoint mismatch {mp.nstr(max_endpoint, 20)}")
        if swirl_profile_error >= ERROR_THRESHOLD:
            raise AssertionError(
                "endpoint swirl-Z profile mismatch "
                f"{mp.nstr(swirl_profile_error, 20)}"
            )
        if not phase_zero["source_budget_ok"]:
            raise AssertionError("resolved reshape B/BZ source budget exceeded 2 A")
        if max_integral >= ERROR_THRESHOLD:
            raise AssertionError(
                f"normalized quadrature mismatch {mp.nstr(max_integral, 20)}"
            )
        if max_tangent >= ERROR_THRESHOLD:
            raise AssertionError(
                f"normalized Z tangent mismatch {mp.nstr(max_tangent, 20)}"
            )
        if not all(abs(value) > SEED_THRESHOLD for value in all_seed_values):
            raise AssertionError("a prescribed seed vanished")
        if not all(abs(value) > SEED_THRESHOLD for value in all_increment_values):
            raise AssertionError("a reshape increment vanished")

        report = {
            "passed": True,
            "precision": precision,
            "pressure_order": PRESSURE_ORDER,
            "width_order": WIDTH_ORDER,
            "reshape_A": mp.nstr(RESHAPE_A, 30),
            "reshape_T": mp.nstr(mp.mpf(RESHAPE_T), 30),
            "logC": mp.nstr(LOG_C, 30),
            "quadrature_order": QUADRATURE_ORDER,
            "tail_digits": TAIL_DIGITS,
            "phases": [mp.nstr(mp.mpf(str(value)), 12) for value in phases],
            "finite_difference_step": mp.nstr(step, 20),
            "maximum_scaled_normalized_integral_error": mp.nstr(max_integral, 30),
            "maximum_scaled_normalized_Z_tangent_error": mp.nstr(max_tangent, 30),
            "maximum_scaled_phase_zero_error": mp.nstr(max_endpoint, 30),
            "phase_rows": phase_rows,
            "normalized_integral_errors": {
                key: mp.nstr(value, 24)
                for key, value in sorted(integral_errors.items())
            },
            "normalized_Z_tangent_errors": {
                key: mp.nstr(value, 24)
                for key, value in sorted(tangent_errors.items())
            },
            "phase_zero_endpoint_errors": {
                key: mp.nstr(value, 24)
                for key, value in sorted(endpoint_errors.items())
            },
            "endpoint_swirl_Z_profile": mp.nstr(raw_swirl_Z / raw_swirl, 30),
            "endpoint_swirl_Z_profile_expected": mp.nstr(swirl_profile, 30),
            "endpoint_swirl_Z_profile_error": mp.nstr(swirl_profile_error, 24),
            "source_B_budget_ok": bool(phase_zero["source_budget_ok"]),
            "endpoint_seeds_preserved_separately": True,
            "increments_preserved_separately": True,
            "axial_field_constant_checked": True,
            "independent_quadrature": "mp.quad split over resolved ell intervals",
            "independent_Z_tangent": "fourth-order centered finite difference of scalar quadrature",
            "quadrature_error_enclosed": False,
            "pressure_width_truncation_remainder_enclosed": False,
            "C2_claim": False,
            "cone_certified": False,
            "finite_energy_certified": False,
            "temporal_recursion_certified": False,
            "global_field_installed": False,
            "metadata": reshape.metadata(),
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
