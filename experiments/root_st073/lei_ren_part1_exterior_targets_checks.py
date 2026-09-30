"""Independent heat identity and source collar moment-transfer receipt."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from lei_ren_part1_exterior_targets import exterior_targets
from lei_ren_part1_heat_moments import heat_tail_moments
from lei_ren_part1_pressure_core import load_core
from lei_ren_part1_heat_collar import HeatCollar
from openai_ns_reconstruction.heat_exterior import HeatExterior


def run():
    heat = HeatExterior(h=0.001, c_inf=0.1)
    rb, ell, epsilon = 2048.0, 0.5, 0.002 / 100
    ra = rb * math.exp(-ell)
    core, _ = load_core()
    collar = HeatCollar(R_b=rb, ell=ell, h=heat.h, c_inf=heat.c_inf,
                        epsilon=epsilon)

    def fh(r, z):
        return heat.profile_E(r, z) / math.sqrt(2 * r)

    def fcollar(r, z):
        y = math.log(rb / r)
        factor = 1.0 if y <= 0 else 1 - epsilon * math.exp(-4 / y**2)
        return factor * fh(r, z)

    rows = []
    identity_errors = []
    refine_errors = []
    derivative_errors = []
    handoff_errors = []
    for z in np.linspace(-1, 1, 33):
        exact = heat_tail_moments(ra, float(z), heat)
        transferred_heat = exterior_targets(fh, heat, ra, rb, float(z))
        identity_errors.extend([
            abs(transferred_heat["angular_target"] - exact["angular_target"]),
            abs(transferred_heat["quadratic_target"] + exact["quadratic_tail"]),
            abs(transferred_heat["pressure_at_inner"] + exact["pressure_tail"]),
        ])
        coarse = exterior_targets(fcollar, heat, ra, rb, float(z), n=32)
        fine = exterior_targets(fcollar, heat, ra, rb, float(z), n=64)
        actual = exterior_targets(collar.F, collar.heat, ra, rb, float(z), n=64)
        handoff_errors.extend([
            abs(actual["angular_target"] - fine["angular_target"]),
            abs(actual["quadratic_target"] - fine["quadratic_target"]),
            abs(actual["pressure_at_inner"] - collar.P(ra, float(z))),
        ])
        for key in ("angular_target", "quadratic_target", "pressure_at_inner"):
            refine_errors.append(abs(fine[key] - coarse[key]))
        dr = ra * 1e-4
        lo = exterior_targets(fcollar, heat, ra-dr, rb, float(z))
        hi = exterior_targets(fcollar, heat, ra+dr, rb, float(z))
        f = fcollar(ra, float(z))
        expected = {"angular_target": 2*ra*f,
                    "quadratic_target": -ra*f*f,
                    "pressure_at_inner": f*f}
        for key, slope in expected.items():
            observed = (hi[key] - lo[key]) / (2*dr)
            derivative_errors.append(abs(observed - slope) / max(abs(slope), 1e-30))
        fine["angular_change_from_pure_heat"] = (
            fine["angular_target"] - transferred_heat["angular_target"])
        bound = float(core.F(0.005, float(z))) * ra**2
        fine["core_exit_monotone_angular_bound"] = bound
        fine["necessary_angular_margin"] = bound - fine["angular_target"]
        fine["necessary_endpoint_margin"] = (
            float(core.F(0.005, float(z))) - fcollar(ra, float(z)))
        rows.append(fine)
    report = {
        "source": "https://arxiv.org/html/2609.35406v1",
        "source_conditions": "(1.2), (2.21)-(2.22), inward heat collar",
        "parameters": {"h": heat.h, "c_inf": heat.c_inf,
                       "Rb": rb, "ell": ell, "epsilon": epsilon},
        "R_inner": ra,
        "heat_transport_identity_max_absolute_error": max(identity_errors),
        "collar_n32_n64_max_absolute_difference": max(refine_errors),
        "independent_radial_derivative_max_relative_error": max(derivative_errors),
        "actual_collar_handoff_max_absolute_error": max(handoff_errors),
        "rows": rows,
        "minimum_inner_necessary_angular_margin": min(
            row["necessary_angular_margin"] for row in rows),
        "minimum_inner_necessary_endpoint_margin": min(
            row["necessary_endpoint_margin"] for row in rows),
        "pressure_target_requires_recomputed_axis_pressure": True,
        "interior_closure_validated": False,
        "stress_or_recursion_validated": False,
    }
    report["checks_passed"] = bool(max(identity_errors) < 1e-8
                                  and max(refine_errors) < 1e-8
                                  and max(derivative_errors) < 2e-6
                                  and max(handoff_errors) < 1e-8
                                  and all(row["necessary_angular_margin"] > 0
                                          and row["necessary_endpoint_margin"] > 0
                                          for row in rows))
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in report.items() if k != "rows"}, indent=2))
    if not report["checks_passed"]:
        raise AssertionError("exterior target transfer checks failed")
    return report


if __name__ == "__main__":
    run()
