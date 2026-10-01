"""Resolved equivalence fixture for the factored core recurrence."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import mpmath as mp


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lei_ren_part1_paper_amplitude_factored_core import (  # noqa: E402
    factored_core_coefficients,
)
from lei_ren_part1_paper_core_recursion import core_coefficients  # noqa: E402


RADIAL_DEGREE = 4
Z_DERIVATIVE_DEPTH = 3
JET_LENGTH = RADIAL_DEGREE + Z_DERIVATIVE_DEPTH + 1
POINTS = ("-.4", "0", ".3")


def _taylor(function, center, length):
    center = mp.mpf(str(center))
    return [
        mp.diff(function, center, order) / mp.factorial(order)
        for order in range(length)
    ]


def _series_product(left, right, length):
    return [
        sum(left[i] * right[k - i] for i in range(k + 1))
        for k in range(length)
    ]


def _scaled_error(value, reference):
    return abs(value - reference) / max(
        mp.mpf(1), abs(value), abs(reference)
    )


def _compare_at(center):
    center = mp.mpf(str(center))
    f0 = _taylor(lambda x: mp.exp(-2 * x * x), center, JET_LENGTH)
    ell = _taylor(lambda x: -4 * x, center, JET_LENGTH)
    squared_f0 = _taylor(lambda x: mp.exp(-4 * x * x), center, JET_LENGTH)
    u0 = _taylor(lambda x: 4 * x + mp.mpf(".01"), center, JET_LENGTH)
    p0 = _taylor(
        lambda x: -3 / (1 + x * x) ** 2,
        center,
        JET_LENGTH,
    )

    original = core_coefficients(
        center,
        ".01",
        F0_Z_taylor=f0,
        U0_Z_taylor=u0,
        P0_Z_taylor=p0,
        radial_degree=RADIAL_DEGREE,
        precision=170,
    )
    factored = factored_core_coefficients(
        center,
        ".01",
        ell_Z_taylor=ell,
        S_Z_taylor=squared_f0,
        U0_Z_taylor=u0,
        P0_Z_taylor=p0,
        radial_degree=RADIAL_DEGREE,
        precision=170,
    )

    errors = []
    pressure_formula_errors = []
    for row_index, (actual_row, amplitude_row) in enumerate(
        zip(original["F"], factored["A"])
    ):
        reconstructed = _series_product(f0, amplitude_row, len(actual_row))
        for value, reference in zip(reconstructed, actual_row):
            errors.append(_scaled_error(value, reference))
    for name in ("Uz", "P"):
        for actual_row, factored_row in zip(original[name], factored[name]):
            for value, reference in zip(factored_row, actual_row):
                errors.append(_scaled_error(value, reference))

    # Check the physical pressure rule independently of the unfactored replay.
    for row_index in range(1, len(factored["P"])):
        amplitude_convolution = [
            mp.mpf(0)
        ] * len(factored["P"][row_index])
        # The full convolution has all A_i A_j terms with i+j=row_index-1.
        for left_index in range(row_index):
            right_index = row_index - 1 - left_index
            product = _series_product(
                factored["A"][left_index],
                factored["A"][right_index],
                len(amplitude_convolution),
            )
            amplitude_convolution = [
                old + new
                for old, new in zip(amplitude_convolution, product)
            ]
        expected = [
            value / row_index
            for value in _series_product(
                squared_f0,
                amplitude_convolution,
                len(factored["P"][row_index]),
            )
        ]
        for value, reference in zip(factored["P"][row_index], expected):
            pressure_formula_errors.append(_scaled_error(value, reference))

    maximum = max(errors + pressure_formula_errors)
    assert maximum < mp.mpf("1e-130"), maximum
    return {
        "Z": mp.nstr(center, 20),
        "max_scaled_equivalence_error": mp.nstr(max(errors), 30),
        "max_scaled_pressure_formula_error": mp.nstr(
            max(pressure_formula_errors), 30
        ),
        "row_lengths": {
            "A": [len(row) for row in factored["A"]],
            "Uz": [len(row) for row in factored["Uz"]],
            "P": [len(row) for row in factored["P"]],
        },
        "axis_amplitude_row_is_one": all(
            value == (1 if index == 0 else 0)
            for index, value in enumerate(factored["A"][0])
        ),
    }


def run():
    with mp.workdps(180):
        records = [_compare_at(point) for point in POINTS]
        result = {
            "all_checks_passed": True,
            "radial_degree": RADIAL_DEGREE,
            "Z_derivative_depth": Z_DERIVATIVE_DEPTH,
            "jet_length": JET_LENGTH,
            "points": records,
            "data": {
                "F0": "exp(-2 Z^2)",
                "ell": "-4 Z",
                "S": "exp(-4 Z^2)",
                "U0": "4 Z + .01",
                "P0": "-3/(1+Z^2)^2",
                "delta": ".01",
            },
            "comparison": "Reconstructed F rows F0*A against core_coefficients F; Uz and physical P rows compared directly.",
            "scope": "Resolved algebraic fixture only; no actual source generation or global core claim.",
        }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    run()
