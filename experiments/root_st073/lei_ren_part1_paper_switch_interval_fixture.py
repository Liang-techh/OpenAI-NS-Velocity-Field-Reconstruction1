"""Small independent check for the directed schedule switch enclosures."""
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp

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


def main():
    mp.mp.dps = 110
    schedule = SimpleNamespace(decimal_precision=80)
    enclosure = ScheduleEndpointEnclosures(schedule, precision=110)
    checks = []
    for text in ("0.1", "0.3", "0.5", "0.9"):
        x = mp.mpf(text)
        s = enclosure.sigma_interval(enclosure.scalar(text))
        slo, shi = interval_bounds(s)
        assert slo <= sigma_reference(x) <= shi
        refined = {}
        for n in (64, 128, 256):
            value = enclosure.primitive_refined(enclosure.scalar(text), panels=n)
            lo, hi = interval_bounds(value)
            exact = mp.quad(sigma_reference, [0, x])
            assert lo <= exact <= hi, (text, n, lo, exact, hi)
            refined[str(n)] = {
                "lower": mp.nstr(lo, 35),
                "upper": mp.nstr(hi, 35),
                "width": mp.nstr(hi - lo, 20),
            }
        widths = [mp.mpf(refined[str(n)]["width"]) for n in (64, 128, 256)]
        assert widths[1] < widths[0] and widths[2] < widths[1], (text, widths)
        checks.append({"x": text, "sigma": [mp.nstr(slo, 35), mp.nstr(shi, 35)], "primitive": refined})

    grids = {}
    for n in (64, 128):
        grid = enclosure.primitive_grid_enclosures(n)
        assert len(grid) == n + 1
        assert interval_bounds(grid[0]) == (mp.mpf(0), mp.mpf(0))
        assert interval_bounds(grid[-1]) == (mp.mpf("0.5"), mp.mpf("0.5"))
        entries = []
        for i in (6, 19, 32, 57) if n == 64 else (12, 38, 64, 114):
            x = mp.mpf(i) / n
            value = grid[i]
            lo, hi = interval_bounds(value)
            exact = mp.quad(sigma_reference, [0, x])
            assert lo <= exact <= hi, (n, i, lo, exact, hi)
            entries.append({"i": i, "x": mp.nstr(x, 20), "lower": mp.nstr(lo, 35), "upper": mp.nstr(hi, 35), "width": mp.nstr(hi - lo, 20)})
        grids[str(n)] = entries
    for small, large in zip(grids["64"], grids["128"]):
        assert mp.mpf(large["width"]) < mp.mpf(small["width"]), (small, large)

    result = {
        "all_checks_passed": True,
        "directed_interval_arithmetic": True,
        "relative_to_stored_schedule_parameters": True,
        "whole_source_certificate": False,
        "primitive_quadrature_reference": "independent mp.quad only in fixture",
        "point_checks": checks,
        "grid_checks": grids,
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"all_checks_passed": True, "output": str(output)}))


if __name__ == "__main__":
    main()
