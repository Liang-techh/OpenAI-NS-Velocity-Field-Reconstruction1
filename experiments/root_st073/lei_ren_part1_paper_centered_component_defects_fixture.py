"""Independent resolved fixture for the five centered component defects.

The production adapter composes signed log saddle records over the finite
pressure/width ring.  This fixture replays the same rows with ordinary
high-precision scalar quadrature from the resolved R=110 endpoint, including
the centered seed constants, the reference continuation, and the unit
restoration phase.  A separate fourth-order Z stencil checks the automatic
first tangent.  It is a local finite-row check; it does not certify closure,
smallness, or a global field.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import mpmath as mp

try:
    from .lei_ren_part1_paper_axial_dual import AxialDual
    from .lei_ren_part1_paper_axial_primitive import _sigma_mp
    from .lei_ren_part1_paper_centered_component_defects import (
        CenteredComponentDefects,
    )
    from .lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
    from .lei_ren_part1_paper_pressure_width_long_reshape_fixture import (
        _ResolvedR110Provider,
        _build_comparison,
    )
except (ImportError, ValueError):
    from lei_ren_part1_paper_axial_dual import AxialDual
    from lei_ren_part1_paper_axial_primitive import _sigma_mp
    from lei_ren_part1_paper_centered_component_defects import (
        CenteredComponentDefects,
    )
    from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
    from lei_ren_part1_paper_pressure_width_long_reshape_fixture import (
        _ResolvedR110Provider,
        _build_comparison,
    )


PRECISION = 100
PRESSURE_ORDER = 3
WIDTH_ORDER = 2
T = mp.mpf(400)
A = mp.mpf(1)
LOG_C = mp.mpf(-4)
LOG_PSTAR = mp.mpf("45.3")
Z0 = mp.mpf(".3")
FINITE_DIFFERENCE_STEP = mp.mpf("1e-4")
ERROR_THRESHOLD = mp.mpf("1e-8")

MOMENT_NAMES = ("theta", "z", "theta_z", "z_theta", "p")


def _scalar(value: Any) -> mp.mpf:
    """Materialize one production scalar at pressure=width=1."""

    if isinstance(value, AxialDual):
        value = value.value
    if isinstance(value, PressureWidthJet):
        return value.evaluate(pressure=1, width=1)
    if isinstance(value, mp.mpf):
        return value
    return mp.mpf(str(value))


def _scaled(actual: Any, expected: Any) -> mp.mpf:
    a = _scalar(actual)
    b = _scalar(expected)
    return abs(a - b) / max(mp.mpf(1), abs(a), abs(b))


def _relative(actual: Any, expected: Any) -> mp.mpf:
    """Relative error without a unit-scale floor for tiny centered rows."""

    a = _scalar(actual)
    b = _scalar(expected)
    denominator = max(abs(a), abs(b))
    if denominator == 0:
        return mp.mpf(0)
    return abs(a - b) / denominator


def _quad_split(function: Any, endpoint: mp.mpf) -> mp.mpf:
    """Independent fixed split quadrature, retaining both flat endpoints."""

    cuts = {mp.mpf(0), endpoint}
    cut = mp.mpf(1)
    while cut < endpoint:
        cuts.add(cut)
        cuts.add(endpoint - cut)
        cut *= 2
    cuts = sorted(cuts)
    return mp.fsum(
        mp.quad(function, [left, right])
        for left, right in zip(cuts, cuts[1:])
    )


def _quad_interval(function: Any, left: mp.mpf, right: mp.mpf) -> mp.mpf:
    """Split a short continuation interval independently of the reshape."""

    if right <= left:
        return mp.mpf(0)
    midpoint = (left + right) / 2
    return mp.quad(function, [left, midpoint, right])


def _unit_integral(function: Any) -> mp.mpf:
    return mp.quad(function, [0, mp.mpf(".5"), 1])


def _direct_rows(switches: Any, z: mp.mpf) -> dict[int, mp.mpf]:
    """Replay Eq. (10.3) and its centered source rows at scalar Z."""

    source = switches.evaluate_R(110, z)
    f0 = _scalar(source["F"])
    v0 = _scalar(source["Uz"])
    g = v0 - 4 * z

    rm = 110 * mp.exp(10 * (LOG_C + LOG_PSTAR) - 6)
    rho0 = 110 / rm
    alpha = 110 * mp.exp(T) / rm
    am = mp.exp(LOG_PSTAR - mp.mpf(".6")) / (1 + z * z)
    uref110 = mp.exp(-LOG_C) / (1 + z * z)
    b = mp.log(mp.sqrt(220) * f0) + LOG_C + mp.log(1 + z * z)

    moments = {
        name: _scalar(source["moments"][name]) for name in MOMENT_NAMES
    }
    scale = mp.sqrt(2) * rm ** mp.mpf("1.5") * am
    rows = {
        1: moments["z"] / rm - 440 * z / rm,
        2: (moments["theta_z"] - 4 * z * moments["theta"]) / scale,
        3: (
            moments["theta"]
            - mp.mpf(5) / 8 * mp.sqrt(2) * mp.mpf(110) ** mp.mpf("1.5") * uref110
        )
        / scale,
        4: (
            moments["z_theta"]
            - 8 * z * moments["z"]
            + 16 * z * z * 110
            + mp.mpf(5) / 12 * 110 * uref110 * uref110
        )
        / (rm * am * am),
        5: (moments["p"] - mp.mpf("2.5") * uref110 * uref110) / (am * am),
    }

    # Independent source replay on the actual log-radius intervals.  During
    # the reshape, rho = R/Rm = rho0 exp(x), while the reference ratio is
    # rho**.1.  The expm1 forms keep the angular differences separate when
    # they are many orders smaller than the reference field.
    def rho(x: mp.mpf) -> mp.mpf:
        return rho0 * mp.exp(x)

    def reference_ratio(x: mp.mpf) -> mp.mpf:
        return rho(x) ** mp.mpf(".1")

    def angular_ratio_difference(x: mp.mpf) -> mp.mpf:
        q = _sigma_mp(1 - x / T)
        return reference_ratio(x) * mp.expm1(b * q)

    def angular_square_difference(x: mp.mpf) -> mp.mpf:
        q = _sigma_mp(1 - x / T)
        return rho(x) ** mp.mpf(".2") * mp.expm1(2 * b * q)

    # Reshape interval x in [0,T].  These are direct densities, rather than
    # the analytically reduced flat-shape kernels used by the production map.
    rows[1] += _quad_split(lambda x: rho(x) * g, T)
    rows[2] += _quad_split(
        lambda x: rho(x) ** mp.mpf("1.5")
        * (reference_ratio(x) + angular_ratio_difference(x))
        * g,
        T,
    )
    rows[3] += _quad_split(
        lambda x: rho(x) ** mp.mpf("1.5") * angular_ratio_difference(x),
        T,
    )
    rows[4] += _quad_split(
        lambda x: rho(x)
        * (g * g / (am * am) - mp.mpf(".5") * angular_square_difference(x)),
        T,
    )
    rows[5] += _quad_split(
        lambda x: mp.mpf(".5") * angular_square_difference(x),
        T,
    )

    # Reference interval x in [T,T+5], where the angular difference is zero
    # but the constant axial defect continues to transport the rows.
    rows[1] += _quad_interval(lambda x: rho(x) * g, T, T + 5)
    rows[2] += _quad_interval(
        lambda x: rho(x) ** mp.mpf("1.6") * g, T, T + 5
    )
    rows[4] += _quad_interval(
        lambda x: rho(x) * g * g / (am * am), T, T + 5
    )

    # Restoration t in [0,1] has rho=e^(t-2) and
    # g(t)=g*(1-sigma(t)); after t=1 the centered axial source vanishes.
    def restore_rho(t: mp.mpf) -> mp.mpf:
        return mp.exp(t - 2)

    def restore_gate(t: mp.mpf) -> mp.mpf:
        return 1 - _sigma_mp(t)

    rows[1] += _unit_integral(lambda t: restore_rho(t) * g * restore_gate(t))
    rows[2] += _unit_integral(
        lambda t: restore_rho(t) ** mp.mpf("1.6") * g * restore_gate(t)
    )
    rows[4] += _unit_integral(
        lambda t: restore_rho(t)
        * (g * restore_gate(t) / am) ** 2
    )
    return rows


def _fourth_difference(
    values: dict[int, dict[int, mp.mpf]], step: mp.mpf
) -> dict[int, mp.mpf]:
    return {
        row: (-values[row][2] + 8 * values[row][1] - 8 * values[row][-1] + values[row][-2])
        / (12 * step)
        for row in range(1, 6)
    }


def run_fixture(
    *,
    precision: int = PRECISION,
    Z: Any = Z0,
    finite_difference_step: Any = FINITE_DIFFERENCE_STEP,
) -> dict[str, Any]:
    precision = max(80, int(precision))
    z0 = mp.mpf(str(Z))
    step = mp.mpf(str(finite_difference_step))
    if not (-1 < z0 < 1):
        raise ValueError("Z must lie in (-1, 1)")
    if step <= 0 or z0 - 2 * step <= -1 or z0 + 2 * step >= 1:
        raise ValueError("finite difference stencil must remain in |Z| < 1")

    with mp.workdps(precision):
        comparison = _build_comparison(precision)
        switches = _ResolvedR110Provider(comparison)
        # The small resolved fixture intentionally uses the same compact
        # pressure/width rectangle as the independent long-reshape check.
        switches.pressure_order = PRESSURE_ORDER
        switches.width_order = WIDTH_ORDER
        defects = CenteredComponentDefects(
            switches,
            A=A,
            logC=LOG_C,
            logPstar=LOG_PSTAR,
            order=32,
            window=24,
            restore_order=64,
        )
        output = defects.evaluate(z0)
        expected = _direct_rows(switches, z0)
        row_errors = {
            str(row): _scaled(output["d"][row], expected[row])
            for row in range(1, 6)
        }
        row_relative_errors = {
            str(row): _relative(output["d"][row], expected[row])
            for row in range(1, 6)
        }

        fd_rows: dict[int, dict[int, mp.mpf]] = {
            row: {} for row in range(1, 6)
        }
        for index in (-2, -1, 1, 2):
            rows = _direct_rows(switches, z0 + index * step)
            for row in range(1, 6):
                fd_rows[row][index] = rows[row]
        expected_Z = _fourth_difference(fd_rows, step)
        tangent_errors = {
            str(row): _scaled(output["d_Z"][row], expected_Z[row])
            for row in range(1, 6)
        }
        tangent_relative_errors = {
            str(row): _relative(output["d_Z"][row], expected_Z[row])
            for row in range(1, 6)
        }

        endpoint = switches.evaluate_R(110, z0)
        p0_error = _scaled(output["P0"], endpoint["P0"])
        p0_Z_error = _scaled(output["P0_Z"], endpoint["P0_Z"])

        part_rows = {
            str(row): sum(
                1
                for row_parts in output["parts"].values()
                if row in row_parts and _scalar(row_parts[row]) != 0
            )
            for row in range(1, 6)
        }
        parts_nonzero = all(count > 0 for count in part_rows.values())
        max_row_error = max(row_errors.values())
        max_tangent_error = max(tangent_errors.values())
        max_row_relative_error = max(row_relative_errors.values())
        max_tangent_relative_error = max(tangent_relative_errors.values())
        if max_row_error >= ERROR_THRESHOLD:
            raise AssertionError(f"centered row mismatch {mp.nstr(max_row_error, 20)}")
        if max_tangent_error >= ERROR_THRESHOLD:
            raise AssertionError(
                f"centered Z tangent mismatch {mp.nstr(max_tangent_error, 20)}"
            )
        if p0_error >= ERROR_THRESHOLD or p0_Z_error >= ERROR_THRESHOLD:
            raise AssertionError("explicit P0 datum was not preserved")
        if not parts_nonzero:
            raise AssertionError("at least one centered source row has no receipt")

        report = {
            "passed": True,
            "precision": precision,
            "pressure_order": PRESSURE_ORDER,
            "width_order": WIDTH_ORDER,
            "A": mp.nstr(A, 30),
            "T": mp.nstr(T, 30),
            "logC": mp.nstr(LOG_C, 30),
            "logPstar": mp.nstr(LOG_PSTAR, 30),
            "Z": mp.nstr(z0, 30),
            "finite_difference_step": mp.nstr(step, 30),
            "row_errors": {
                row: mp.nstr(value, 30) for row, value in row_errors.items()
            },
            "row_relative_errors": {
                row: mp.nstr(value, 30)
                for row, value in row_relative_errors.items()
            },
            "Z_tangent_errors": {
                row: mp.nstr(value, 30) for row, value in tangent_errors.items()
            },
            "Z_tangent_relative_errors": {
                row: mp.nstr(value, 30)
                for row, value in tangent_relative_errors.items()
            },
            "expected_rows": {
                str(row): mp.nstr(value, 30) for row, value in expected.items()
            },
            "actual_rows": {
                str(row): mp.nstr(_scalar(output["d"][row]), 30)
                for row in range(1, 6)
            },
            "expected_Z_rows": {
                str(row): mp.nstr(value, 30) for row, value in expected_Z.items()
            },
            "actual_Z_rows": {
                str(row): mp.nstr(_scalar(output["d_Z"][row]), 30)
                for row in range(1, 6)
            },
            "maximum_scaled_row_error": mp.nstr(max_row_error, 30),
            "maximum_scaled_Z_tangent_error": mp.nstr(max_tangent_error, 30),
            "maximum_relative_row_error": mp.nstr(max_row_relative_error, 30),
            "maximum_relative_Z_tangent_error": mp.nstr(
                max_tangent_relative_error, 30
            ),
            "P0_error": mp.nstr(p0_error, 30),
            "P0_Z_error": mp.nstr(p0_Z_error, 30),
            "part_row_nonzero_counts": part_rows,
            "all_rows_have_separate_receipts": parts_nonzero,
            "independent_replay": "mp.quad split over centered shape, reference, and restoration rows",
            "independent_Z_tangent": "fourth-order centered finite difference of scalar replay",
            "quadrature_error_enclosed": False,
            "pressure_width_truncation_remainder_enclosed": False,
            "functional_closure": False,
            "smallness_certified": False,
            "cone_certified": False,
            "global_field_installed": False,
            "metadata": output.get("metadata", {}),
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
