"""Source-compatible scale-parameter demonstration for the shared profile.

The paper's core discussion uses ``j``, ``Lambda`` and ``C_*`` before the
prepared outer profile is attached.  This script chooses the necessary
lower-bound demonstration ``log(C_*) = 2 log(Lambda)`` (the requested
``A_Omega = 0`` lower bound only) and sets

    log(R_ref) = log(110) + 10 (log(C_*) + log(P_*)),

then passes one explicit schedule through the shared waiting, angular
correction, energy-tail and pressure objects.  The resulting values are
candidate numerical inputs; this does not certify the paper's theorem
conditions, core construction, pressure target, PDE, cone, or scale recursion.
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

from lei_ren_part1_paper_axial_energy_tail import build_default_tail  # noqa: E402
from lei_ren_part1_paper_corrected_profile import CorrectedSourceProfile  # noqa: E402
from lei_ren_part1_paper_outer import PaperOuterSchedule  # noqa: E402


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
PRECISION = 160
QUADRATURE_ORDER = 128


def _selected_coefficients(profile: CorrectedSourceProfile, Z: float) -> dict[str, Any]:
    data = profile.coefficients(Z)
    incoming = data["incoming"]
    axial = data["axial"]
    angular = data["angular"]
    pressure = profile.pressure_at_reference(Z)
    return {
        "Z": Z,
        "incoming_rows": {
            "m1_Mz_over_RpEp": incoming["m1_Mz_over_RpEp"],
            "m2_Mtheta_z_over_Rp_sqrt2RpEp2": incoming[
                "m2_Mtheta_z_over_Rp_sqrt2RpEp2"
            ],
            "E_prior_mu_Mztheta_over_RpEp2": incoming[
                "E_prior_mu_Mztheta_over_RpEp2"
            ],
            "row_normalization": incoming["row_normalization"],
        },
        "angular_candidate_inputs": {
            "r_pre": angular["r_pre"],
            "heat": angular["heat"],
            "log_r": angular["log_r"],
            "log_s": angular["log_s"],
            "d1": angular["d1"],
            "d2": angular["d2"],
            "scaled_algebraic_residuals": angular["scaled_algebraic_residuals"],
            "input_status": angular["input_status"],
        },
        "axial_candidate_inputs": {
            "a_p": axial["a_p"],
            "c": axial["c"],
            "energy_target": axial["energy_target"],
            "linear_relative_replay": axial["linear_relative_replay"],
            "energy_relative_replay": axial["energy_relative_replay"],
            "requested_quadrature_order": axial["requested_quadrature_order"],
            "quadrature_order": axial["quadrature_order"],
            "full_outer_closed": axial["full_outer_closed"],
        },
        "same_axis_P0": {
            "baseline_axis_P0_over_Pstar_squared": pressure[
                "baseline_axis_P0_over_Pstar_squared"
            ],
            "corrected_axis_P0_over_Pstar_squared": pressure[
                "corrected_axis_P0_over_Pstar_squared"
            ],
            "backward_pressure_change_over_Pstar_squared": pressure[
                "backward_pressure_change_over_Pstar_squared"
            ],
            "exact_pressure_target": pressure["exact_pressure_target"],
        },
    }


def run() -> dict[str, Any]:
    with mp.workdps(PRECISION):
        j = mp.mpf(".02")
        Lambda = mp.mpf("2500")
        logPstar = mp.mpf("14")
        logCstar = 2 * mp.log(Lambda)
        log_cstar_diff = logCstar - 2 * mp.log(Lambda)
        logRref = mp.log(110) + 10 * (logCstar + logPstar)
        base = PaperOuterSchedule(
            logPstar=mp.nstr(logPstar, PRECISION),
            logRref=mp.nstr(logRref, PRECISION),
            delta="1e-32",
            Md=".5",
            c_mu=".001",
            c_delta=".001",
            c_epsilon=".01",
        )
        # One explicit base schedule is passed through the waiting-root
        # factory, then the same returned tail and correction are injected
        # into CorrectedSourceProfile and its pressure adapter.
        tail = build_default_tail(
            precision=PRECISION,
            quadrature_order=QUADRATURE_ORDER,
            schedule=base,
            match_waiting=True,
        )
        profile = CorrectedSourceProfile(
            precision=PRECISION,
            order=QUADRATURE_ORDER,
            schedule=tail.schedule,
            angular=tail.correction,
            tail=tail,
        )
        relation_error = mp.mpf(str(profile.schedule.logRref)) - logRref
        rows = [_selected_coefficients(profile, z) for z in (0.0, 0.3)]

    report = {
        "source": SOURCE,
        "source_version": SOURCE_VERSION,
        "source_parameters": {
            "j": str(j),
            "Lambda": str(Lambda),
            "logLambda": mp.nstr(mp.log(Lambda), PRECISION),
            "logCstar": mp.nstr(logCstar, PRECISION),
            "logCstar_minus_2logLambda": mp.nstr(log_cstar_diff, PRECISION),
            "logPstar": mp.nstr(logPstar, PRECISION),
            "logRref_requested": mp.nstr(logRref, PRECISION),
            "logRref_relation_error": mp.nstr(relation_error, PRECISION),
            "A_Omega_lower_bound_used": "0",
            "A_Omega_status": "zero lower-bound demonstration only; not certified",
        },
        "schedule_requested": base.metadata(),
        "schedule_active_after_waiting_match": profile.schedule.metadata(),
        "shared_construction": profile.metadata()["construction"],
        "rows": rows,
        "checks": {
            "logCstar_ge_2logLambda": bool(log_cstar_diff >= 0),
            "logRref_formula_reproduced": bool(abs(relation_error) < mp.mpf("1e-120")),
            "shared_tail_schedule_identity": profile.tail.schedule is profile.schedule,
            "shared_angular_schedule_identity": profile.angular.schedule is profile.schedule,
            "pressure_uses_same_angular_object": profile.pressure.correction is profile.angular,
            "source_theorem_full_conditions_certified": False,
            "same_axis_P0_exactly_restored": False,
        },
        "scope": (
            "Necessary source-scaled parameter demonstration and shared candidate "
            "input receipt only. j, Lambda, C_star and A_Omega choices are not a "
            "certified core theorem construction; pressure, five-moment closure, "
            "PDE, cone and scale recursion remain unresolved."
        ),
    }
    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "logCstar_ge_2logLambda": report["checks"]["logCstar_ge_2logLambda"],
                "logRref": report["source_parameters"]["logRref_requested"],
                "same_axis_P0": [
                    {
                        "Z": row["Z"],
                        "baseline": row["same_axis_P0"][
                            "baseline_axis_P0_over_Pstar_squared"
                        ],
                        "corrected": row["same_axis_P0"][
                            "corrected_axis_P0_over_Pstar_squared"
                        ],
                    }
                    for row in rows
                ],
            },
            default=str,
        )
    )
    return report


if __name__ == "__main__":
    run()
