"""Independent pointwise checks for the finite flatten-stage enclosure."""

import json
from decimal import Decimal
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp

from lei_ren_part1_paper_flatten_preheat_integrals import (
    integrate_flatten,
    integrate_flatten_core,
)
from lei_ren_part1_paper_schedule_endpoint_enclosures import (
    ScheduleEndpointEnclosures,
    endpoints,
)


def sigma_reference(x):
    if x <= 0:
        return mp.mpf(0)
    if x >= 1:
        return mp.mpf(1)
    phase = 1 / (x * x) - 1 / ((1 - x) * (1 - x))
    if phase >= 0:
        tiny = mp.exp(-phase)
        return tiny / (1 + tiny)
    tiny = mp.exp(phase)
    return 1 / (1 + tiny)


def interval_bounds(value):
    lo, hi = endpoints(value)
    return mp.mpf(lo), mp.mpf(hi)


def reference_integrands(u, z, length, mu):
    beta = 2 * (1 - sigma_reference(u))
    ell = (-mp.mpf(".5") - mu) * length * u
    q = 1 + z * z
    density = mp.exp(2 * ell + (beta - 3) * mp.log(2) - beta * mp.log(q))
    first = density * (-2 * beta * z / q)
    second = density * (-2 * beta / q + 4 * beta * (beta + 1) * z * z / (q * q))
    return density, first, second


def main():
    with mp.workdps(110):
        enclosure = ScheduleEndpointEnclosures(SimpleNamespace(decimal_precision=90), precision=100)
        z = mp.mpf(".3")
        length = mp.mpf(1)
        mu = mp.mpf(".01")
        result = integrate_flatten_core(
            enclosure,
            base=enclosure.iv.mpf(0),
            length=enclosure.iv.mpf(1),
            mu=enclosure.iv.mpf(".01"),
            Z=".3",
            panels=32,
            order=8,
        )
        exact = []
        for component in range(3):
            value = mp.quad(
                lambda u: reference_integrands(u, z, length, mu)[component],
                [0, mp.mpf(".25"), mp.mpf(".5"), mp.mpf(".75"), 1],
            )
            exact.append(value)
        intervals = [
            result["integrals"]["value"],
            result["integrals"]["first_Z"],
            result["integrals"]["second_Z"],
        ]
        checks = []
        for name, interval, target in zip(("value", "first_Z", "second_Z"), intervals, exact):
            lower, upper = interval_bounds(interval)
            assert lower <= target <= upper, (name, lower, target, upper)
            checks.append({
                "name": name,
                "lower": mp.nstr(lower, 35),
                "upper": mp.nstr(upper, 35),
                "reference": mp.nstr(target, 35),
            })

        # Exercise the public schedule wrapper with a resolved stored length
        # while keeping the fixture independent of the expensive source build.
        fake_schedule = SimpleNamespace(
            decimal_precision=90,
            y_v=Decimal(0),
            y_f=Decimal(100),
            mu=Decimal(".01"),
            delta=Decimal(".001"),
            _stage_bounds={"z_flatten": (Decimal(0), Decimal(100))},
        )
        wrapped_enclosure = ScheduleEndpointEnclosures(fake_schedule, precision=100)
        wrapped_enclosure.log_amplitude_ratio = lambda _y: {"interval": wrapped_enclosure.iv.mpf(0)}
        wrapped = integrate_flatten(wrapped_enclosure, Z=".3", panels=16, order=6)
        assert wrapped["stored_length_is_100"]
        assert wrapped["pointwise_Z_only"] and not wrapped["uniform_Z_bound"]

        report = {
            "all_checks_passed": True,
            "directed_interval_arithmetic": True,
            "independent_scalar_quadrature": True,
            "pointwise_Z_only": True,
            "source_certificate": False,
            "original_parameter_errors_enclosed": False,
            "core_checks": checks,
            "core_metadata": {
                "panels": result["panels"],
                "Taylor_order": result["Taylor_order"],
                "edge_method": result["edge_method"],
                "edge_precision_limitation": result["edge_precision_limitation"],
                "interval_widths": {
                    key: mp.nstr(value, 25)
                    for key, value in result["interval_widths"].items()
                },
            },
            "wrapper_stored_length": str(wrapped["stored_flatten_length"]),
        }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"all_checks_passed": True, "output": str(output)}))


if __name__ == "__main__":
    main()
