"""Regular, source-bound Lei--Ren Omega_0/R on the admitted compliant core.

This exposes the known first pressure forcing from the current v2 compliant
leading core. It is an input/source packet only: it does not solve F_1, Uz_1,
P_1 or certify Assumption 14.1 on the full required interval.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_core_physical_field import (
    CompliantCorePhysicalField, gridkey, Q, V,
)
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent
PREFIX = "lei_ren_part1_paper_compliant_"
NAME = PREFIX + "omega0_regular_source.json"
CHECK = PREFIX + "omega0_regular_source_check.json"


def _sha(path):
    return hashlib.sha256((HERE / path).read_bytes()).hexdigest()


def regular_omega0_over_R(core, rho, z, profile=None):
    """Evaluate Omega_0/R using V_0=R*Q, including the regular axis value.

    All radial derivatives are converted from rho=Lambda*R before use. The
    expression has no division by R, so rho=0 is included directly.
    """
    c = core.ctx
    rho, z = c.mpf(rho), c.mpf(z)
    profile = profile or core.profiles(rho, z)
    g = profile["ordinary_mixed_profile_grids"]
    q = g[Q][gridkey(0, 0)]
    qr = core.Lambda * g[Q][gridkey(1, 0)]
    qrr = core.Lambda**2 * g[Q][gridkey(2, 0)]
    qz = g[Q][gridkey(0, 1)]
    uz = g[V][gridkey(0, 0)]
    d = 1 - z**2
    L = 1 - core.delta * z**2
    D = (1 - core.delta) / 2
    rho_qrho = rho * g[Q][gridkey(1, 0)]

    # T_0(RQ)/R, the i=j=0 transport row divided by R,
    # Uz_0*Z_0(RQ)/R, and the radial viscosity row divided by R.
    time = (q + D*z*qz + rho_qrho) / L
    transport = q * (q/2 + rho_qrho)
    axial = uz * (-2*z*q + d*qz - 2*z*rho_qrho) / L
    viscosity = -4*qr - 2*(rho/core.Lambda)*qrr
    return time + transport + axial + viscosity


def run():
    core = CompliantCorePhysicalField()
    rows = []
    for rho in ("0", "1", "2", "4"):
        for z in ("-0.5", "0", "0.5"):
            value = regular_omega0_over_R(core, rho, z)
            lo, hi = endpoints(value)
            if not (mp.isfinite(lo) and mp.isfinite(hi)):
                raise ArithmeticError("nonfinite regular Omega_0/R enclosure")
            rows.append(dict(rho=rho, Z=z, Omega0_over_R=encode(value),
                             known_pressure_source_over_R=encode(-value/2)))
    result = dict(
        source="Lei-Ren Part I v2 equations (13.16)-(13.19)",
        source_sha256=core.source,
        actual_five_defect_family_sha256=core.family,
        normalized_radial_coordinate="rho=Lambda*R",
        common_core_domain_rho=[0, "4.1"], Z_domain=[-1, 1],
        rows=rows,
        regular_factorization="V0=R*Q; evaluate Omega0/R without division by R",
        known_pressure_source="-Omega0/(2R)",
        coefficient_F1_Uz1_P1_solved=False,
        assumption_14_1_full_interval_admitted=False,
        global_admissible_stress_lift_constructed=False,
        flat_remainder_certified=False,
        input_hashes=dict(core.hashes),
    )
    (HERE / NAME).write_text(json.dumps(encode(result), indent=2) + "\n", encoding="utf8")
    print(f"Wrote {len(rows)} source-bound regular Omega_0/R packets")
    return result


if __name__ == "__main__":
    run()
