"""Consistent numerical replay of the actual Section 7.5 axial tail.

The source equations first solve the two normalized linear rows (7.31), then
choose ``a_p`` from the energy equation (7.34).  A coefficient solve at one
quadrature order followed by a primitive replay at another order leaves a
small apparent exterior axial mean.  This module uses one canonical order
for the bump matrix, the energy weights, and the terminal primitive replay.

The result is a numerical consistency repair, not an assertion that the
continuous source integrals are exactly closed.  Incoming moments and the
heat-tail target remain the actual profile-derived inputs supplied by the
existing receipts.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

import mpmath as mp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SRC = ROOT / "src"
for path in (SRC, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_paper_axial_correction import (  # noqa: E402
    CONSISTENT_QUADRATURE_ORDER,
    from_signed_log,
    independent_end_bump_replay,
    signed_log,
    solve_actual_axial,
)
from lei_ren_part1_paper_axial_energy_tail import build_default_tail  # noqa: E402
from lei_ren_part1_paper_axial_incoming import incoming_axial_moments  # noqa: E402
from lei_ren_part1_paper_axial_pulse_moments import (  # noqa: E402
    normalized_pulse_integral,
)


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_EQUATIONS = "(2.21)-(2.22), (7.29), (7.31), (7.34)"


def _energy_replay(
    receipt: dict[str, Any],
    *,
    mu: Any,
    target: Any,
    precision: int,
) -> dict[str, Any]:
    """Recompute the source quadratic with the same stored weights."""

    with mp.workdps(precision):
        mu_mp = mp.mpf(str(mu))
        target_mp = mp.mpf(str(target))
        a = mp.mpf(str(receipt["a_p"]))
        coefficients = [from_signed_log(row) for row in receipt["c"]]
        kp = mp.mpf(str(receipt["K_p"]))
        kb = [mp.mpf(str(value)) for value in receipt["K_bump"]]
        lhs = kp * a * a + mu_mp * sum(
            weight * coefficient * coefficient
            for weight, coefficient in zip(kb, coefficients)
        )
        defect = lhs - target_mp
        return {
            "lhs": mp.nstr(lhs, precision),
            "target": mp.nstr(target_mp, precision),
            "absolute_defect": signed_log(defect, precision),
            "relative_defect": mp.nstr(abs(defect / target_mp), 40)
            if target_mp
            else None,
        }


def solve_consistent_actual_tail(
    Z: Any,
    *,
    precision: int = 160,
    requested_order: int = 128,
    closure_order: int = CONSISTENT_QUADRATURE_ORDER,
) -> dict[str, Any]:
    """Solve and replay the actual candidate with one shared quadrature order."""

    requested_closure_order = int(closure_order)
    if requested_closure_order != closure_order or requested_closure_order < 16:
        raise ValueError("closure_order must be an integer at least 16")
    if int(requested_order) != requested_order or requested_order < 16:
        raise ValueError("requested_order must be an integer at least 16")
    if precision < 80:
        raise ValueError("precision must be at least 80")
    # Never let a caller request a replay order below the profile's canonical
    # terminal-bump order.  Every input, solve, and replay then shares one
    # effective quadrature problem.
    closure_order = max(requested_closure_order, CONSISTENT_QUADRATURE_ORDER)

    tail = build_default_tail(precision=precision, quadrature_order=closure_order)
    schedule = tail.schedule
    incoming = incoming_axial_moments(schedule, Z, order=64, precision=precision)
    future = tail.evaluate(Z, quadrature_order=closure_order)
    base = [
        incoming["row_normalization"][key]
        for key in ("scaled_base_m1", "scaled_base_m2")
    ]
    pulse = [
        normalized_pulse_integral(
            schedule.mu,
            row,
            order=closure_order,
            precision=precision,
        )
        for row in (1, 2)
    ]
    with mp.workdps(precision):
        prior = from_signed_log(incoming["E_prior_mu_Mztheta_over_RpEp2"])
        target = (
            (1 - mp.exp(-26)) / 4
            - prior
            + mp.mpf(future["energy_target_contribution_nominal"])
        )
        axial = solve_actual_axial(
            schedule.mu,
            base,
            pulse,
            mp.nstr(target, precision),
            precision=precision,
            # Keep the caller's requested order visible.  The owned solver
            # promotes it to CONSISTENT_QUADRATURE_ORDER, so this exercises
            # the same 128-to-192 compatibility path used by the profile.
            order=requested_order,
        )
    replay = independent_end_bump_replay(
        schedule.mu,
        axial,
        base,
        pulse,
        order=closure_order,
        precision=precision,
    )
    energy = _energy_replay(axial, mu=schedule.mu, target=target, precision=precision)
    row_threshold = mp.mpf("1e-30")
    energy_threshold = mp.mpf("1e-40")
    row_ok = all(mp.mpf(value) <= row_threshold for value in replay["row_relative_differences"])
    energy_ok = mp.mpf(energy["relative_defect"]) <= energy_threshold
    return {
        "Z": float(Z),
        "requested_quadrature_order": int(requested_order),
        "requested_closure_quadrature_order": requested_closure_order,
        "closure_quadrature_order": int(closure_order),
        "incoming_moments": incoming,
        "future_angular_energy": future,
        "pulse_rows": pulse,
        "axial_solution": axial,
        "terminal_normalized_mass_replay": replay,
        "energy_replay": energy,
        "checks_passed": bool(row_ok and energy_ok),
        "checks": {
            "linear_rows_threshold": str(row_threshold),
            "energy_threshold": str(energy_threshold),
            "linear_rows_passed": bool(row_ok),
            "energy_passed": bool(energy_ok),
        },
        "scope": (
            "Actual source candidate inputs with a shared finite quadrature order. "
            "The repair removes order-mismatch drift; it does not certify exact "
            "continuous moment closure, pressure, cone, PDE, or global admissibility."
        ),
    }


def run() -> dict[str, Any]:
    rows = [
        solve_consistent_actual_tail(
            z,
            precision=160,
            requested_order=128,
            closure_order=CONSISTENT_QUADRATURE_ORDER,
        )
        for z in (0.0, 0.5, -0.5)
    ]
    report = {
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_equations": SOURCE_EQUATIONS,
        "canonical_quadrature_order": CONSISTENT_QUADRATURE_ORDER,
        "rows": rows,
        "checks_passed": all(row["checks_passed"] for row in rows),
        "mathematical_limitation": (
            "A shared 192-point replay removes the observed 128-versus-192 "
            "quadrature mismatch. It remains a finite quadrature result with "
            "inherited incoming and heat-tail numerical uncertainty; exact source "
            "integral closure is not inferred."
        ),
        "adapter_integration": (
            "Call solve_actual_axial with closure_order=192 and use the returned "
            "quadrature_order for every terminal primitive and energy replay. Do "
            "not overwrite the terminal mass by hand."
        ),
        "full_outer_closed": False,
    }
    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "checks_passed": report["checks_passed"],
                "linear_replay": [
                    row["terminal_normalized_mass_replay"]["row_relative_differences"]
                    for row in rows
                ],
                "energy_replay": [row["energy_replay"]["relative_defect"] for row in rows],
            }
        )
    )
    return report


if __name__ == "__main__":
    run()
