"""Independent substitution check for the regular Omega_0/R source packet."""
import hashlib
import json
from pathlib import Path

import sympy as sp

from lei_ren_part1_paper_compliant_core_physical_field import (
    CompliantCorePhysicalField, gridkey, Q, V,
)
from lei_ren_part1_paper_compliant_omega0_regular_source import NAME, regular_omega0_over_R
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent
PREFIX = "lei_ren_part1_paper_compliant_"
OUT = PREFIX + "omega0_regular_source_check.json"


def _sha(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


def _symbolic_regular_identity():
    R, z, delta, q, qr, qrr, qz, uz = sp.symbols("R z delta q q_R q_RR q_Z U_z")
    d, L, D = 1-z**2, 1-delta*z**2, (1-delta)/2
    V0 = R*q
    VR = q+R*qr
    VRR = 2*qr+R*qrr
    VZ = R*qz
    T0V_over_R = (D*z*VZ+R*VR)/(L*R)
    Z0V_over_R = (d*VZ-2*z*R*VR)/(L*R)
    transport_over_R = V0*(VR-V0/(2*R))/R + uz*Z0V_over_R
    viscosity_over_R = -2*R*VRR/R
    direct = sp.cancel(T0V_over_R+transport_over_R+viscosity_over_R)
    regular = (q+D*z*qz+R*qr)/L + q*(q/2+R*qr) \
        + uz*(-2*z*q+d*qz-2*z*R*qr)/L - 4*qr-2*R*qrr
    return sp.simplify(direct-regular) == 0


def run():
    record = json.loads((HERE / NAME).read_bytes())
    core = CompliantCorePhysicalField()
    failures = []
    require = lambda ok, why: failures.append(why) if not ok else None
    require(record["source_sha256"] == core.source, "source family hash")
    require(record["actual_five_defect_family_sha256"] == core.family, "moment family hash")
    require(record["regular_factorization"] == "V0=R*Q; evaluate Omega0/R without division by R", "regular-axis declaration")
    require(record["coefficient_F1_Uz1_P1_solved"] is False, "n=1 scope flag")
    require(record["assumption_14_1_full_interval_admitted"] is False, "Assumption 14.1 scope flag")
    require(_symbolic_regular_identity(), "symbolic V0=R*Q substitution")

    for name, digest in record["input_hashes"].items():
        require(_sha(name) == digest, "stale input: " + name)

    checked = []
    for row in record["rows"]:
        c = core.ctx
        rho, z = c.mpf(row["rho"]), c.mpf(row["Z"])
        profile = core.profiles(rho, z)
        grids = profile["ordinary_mixed_profile_grids"]
        # Independently substitute physical-R jets into the unsimplified
        # Omega formula at R>0. At R=0 the source is checked by the proven
        # regular factorization and finite derivative enclosures.
        q = grids[Q][gridkey(0, 0)]
        q_rho = grids[Q][gridkey(1, 0)]
        q_rhorho = grids[Q][gridkey(2, 0)]
        q_z = grids[Q][gridkey(0, 1)]
        uz = grids[V][gridkey(0, 0)]
        value = regular_omega0_over_R(core, rho, z, profile)
        stored = read_interval(c, row["Omega0_over_R"])
        psource = read_interval(c, row["known_pressure_source_over_R"])
        slo, shi = endpoints(stored)
        vlo, vhi = endpoints(value)
        require(not (shi < vlo or vhi < slo), f"regular enclosure {rho},{z}")
        plo, phi = endpoints(psource)
        elo, ehi = endpoints(-value/2)
        require(not (phi < elo or ehi < plo), f"pressure source {rho},{z}")
        if endpoints(rho)[0] > 0:
            R = rho/core.Lambda
            qr = core.Lambda*q_rho
            qrr = core.Lambda**2*q_rhorho
            d = 1-z*z; L = 1-core.delta*z*z; D=(1-core.delta)/2
            V0=R*q; VR=q+R*qr; VRR=2*qr+R*qrr; VZ=R*q_z
            direct=((D*z*VZ+R*VR)/L + V0*(VR-V0/(2*R))
                    + uz*(d*VZ-2*z*R*VR)/L - 2*R*VRR)/R
            dlo,dhi=endpoints(direct)
            require(not (shi < dlo or dhi < slo), f"direct Omega quotient {rho},{z}")
        checked.append(dict(rho=row["rho"],Z=row["Z"],positive_radius_direct_check=endpoints(rho)[0]>0))

    result = dict(all_passed=not failures, failures=failures,
        symbolic_regular_identity=True, actual_source_rows_checked=len(checked),
        positive_radius_direct_rows=sum(r["positive_radius_direct_check"] for r in checked),
        axis_rows_regularized=sum(not r["positive_radius_direct_check"] for r in checked),
        n1_coefficient_solved=False, assumption_14_1_full_interval_admitted=False,
        global_admissible_stress_lift_constructed=False, flat_remainder_certified=False,
        checked_points=checked,
        input_hashes={NAME:_sha(NAME), Path(__file__).name:_sha(Path(__file__).name),
                      "lei_ren_part1_paper_compliant_omega0_regular_source.py":_sha("lei_ren_part1_paper_compliant_omega0_regular_source.py")})
    (HERE / OUT).write_text(json.dumps(encode(result), indent=2)+"\n", encoding="utf8")
    print("Regular Omega_0/R independent check:", "PASS" if result["all_passed"] else "FAIL")
    if failures:
        raise RuntimeError("; ".join(failures))
    return result


if __name__ == "__main__":
    run()
