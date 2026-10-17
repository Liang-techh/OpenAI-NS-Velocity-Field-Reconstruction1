"""Bounded independent checks for the new weighted pulse moments."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from openai_ns_reconstruction.original_weighted_pulse_moments import (  # noqa: E402
    build_weighted_moments,
    positive_weighted_rectangle_bounds,
)


def run_fixture() -> dict[str, object]:
    ctx = mp.mp.clone()
    ctx.dps = 90
    iv = MPIntervalContext()
    iv.dps = 140
    controls: list[dict[str, object]] = []

    def check(name: str, condition: bool) -> None:
        if not condition:
            raise AssertionError(name)
        controls.append({"name": name, "passed": True})

    def ends(box: object) -> tuple[object, object]:
        lower, upper = box._mpi_  # type: ignore[attr-defined]
        return ctx.make_mpf(lower), ctx.make_mpf(upper)

    # These two positive, decreasing test functions have independent closed
    # integrals.  Their values are never obtained from the production pulse
    # quadrature, so this exercises the generic rectangle helper itself.
    fixtures = (
        (
            "kappa_one",
            lambda t: 1,
            lambda s, power: s ** (power + 1) / (power + 1),
        ),
        (
            "kappa_one_minus_t",
            lambda t: 1 - t,
            lambda s, power: s ** (power + 1) / (power + 1)
            - s ** (power + 2) / (power + 2),
        ),
    )
    for fixture_name, kappa, exact_integral in fixtures:
        for phase_text in ("0.25", "0.5", "1"):
            phase = ctx.mpf(phase_text)
            for power in (1, 2):
                coarse, coarse_meta = positive_weighted_rectangle_bounds(
                    ctx, iv, phase, power, 128, kappa, return_details=True
                )
                fine, fine_meta = positive_weighted_rectangle_bounds(
                    ctx, iv, phase, power, 256, kappa, return_details=True
                )
                exact = exact_integral(phase, power)
                coarse_lower, coarse_upper = ends(coarse)
                fine_lower, fine_upper = ends(fine)
                check(
                    f"{fixture_name}_{phase_text}_power{power}_coarse_contains_exact",
                    coarse_lower <= exact <= coarse_upper,
                )
                check(
                    f"{fixture_name}_{phase_text}_power{power}_fine_contains_exact",
                    fine_lower <= exact <= fine_upper,
                )
                check(
                    f"{fixture_name}_{phase_text}_power{power}_refines",
                    coarse_lower <= fine_lower <= fine_upper <= coarse_upper,
                )
                check(
                    f"{fixture_name}_{phase_text}_power{power}_midpoint_not_selected",
                    coarse_meta["rectangle_midpoint_used_as_value"] is False
                    and fine_meta["rectangle_midpoint_used_as_value"] is False,
                )

    report = build_weighted_moments(
        ctx,
        iv,
        ("0", "0.125", "0.5", "1"),
        powers=(1, 2),
        cells=(128, 256),
    )
    rows = report["rows"]
    check("pulse_build_passed", report["passed"] is True)
    check("pulse_row_count", len(rows) == 8)
    check("pulse_has_no_M0_rows", all(row["power"] in (1, 2) for row in rows))
    check(
        "pulse_rows_have_exact_serialized_boxes",
        all(
            isinstance(row["source_enclosure_exact_mpf_tuples"], list)
            and len(row["source_enclosure_exact_mpf_tuples"]) == 2
            for row in rows
        ),
    )
    zero_rows = [row for row in rows if row["phase"] == "0"]
    check(
        "pulse_s0_exact_zero",
        bool(zero_rows)
        and all(
            row["approximation"] == "0"
            and row["absolute_error_upper"] == "0"
            and row["approximation_is_exact_source_value"] is True
            for row in zero_rows
        ),
    )
    check(
        "pulse_interiors_are_not_promoted",
        all(
            row["approximation_is_exact_source_value"] is False
            and row["source_value_not_selected"] is True
            for row in rows
            if row["phase"] != "0"
        ),
    )
    check(
        "pulse_refinement_passed",
        all(row["refinement"]["passed"] is True for row in rows),
    )
    check(
        "pulse_uses_independent_quadrature_and_boxes",
        all(
            row["quadrature"]["method"]
            == "independent scalar adaptive quadrature (ctx.quad)"
            and row["boxes"]
            and row["rectangle_midpoint_used_as_value"] is False
            for row in rows
        ),
    )

    source_text = (
        ROOT / "src/openai_ns_reconstruction/original_weighted_pulse_moments.py"
    ).read_text(encoding="utf-8")
    check(
        "old_M0_producer_not_imported",
        "original_interior_mass_source" not in source_text
        and "build_interior_mass_source" not in source_text,
    )
    check(
        "no_closed_W1_W2_symmetry_claim",
        report["source"]["closed_symmetry_values_claimed"] is False,
    )
    return {
        "passed": True,
        "control_count": len(controls),
        "controls": controls,
        "independent_oracle": "analytic kappa=1 and kappa=1-t positive weighted integrals",
        "accepted_M0_mass_recomputed": False,
        "old_checker_or_owner_called": False,
    }


if __name__ == "__main__":
    print(json.dumps(run_fixture(), indent=2))
