"""Actual angular closure, derivative and pressure handoff diagnostics."""
import json
from pathlib import Path

import numpy as np

from lei_ren_part1_extended_swirl import ExtendedSwirl


def run():
    profile=ExtendedSwirl()
    rows=[]
    derivative_errors=[]
    positive=True
    negative_shear=True
    for z in np.linspace(-1,1,33):
        z=float(z)
        row=profile.moments(z,n=256)
        samples=np.geomspace(1e-5,profile.collar.R_b,192)
        fs=np.array([profile.F(r,z) for r in samples])
        frs=np.array([profile.F_R(r,z) for r in samples])
        positive=positive and bool(np.all(fs>0))
        negative_shear=negative_shear and bool(np.all(frs<0))
        row["minimum_F_sampled"]=float(fs.min())
        row["maximum_F_R_sampled"]=float(frs.max())
        # Check each layer with independent finite differences.
        rc=profile._cross(z)
        ra=profile.collar.collar_inner_radius
        for r in (.0075, np.sqrt(profile.R_anchor*rc), np.sqrt(rc*ra),
                  np.sqrt(ra*profile.collar.R_b)):
            dr=r*2e-5
            fd=(profile.F(r+dr,z)-profile.F(r-dr,z))/(2*dr)
            exact=profile.F_R(r,z)
            derivative_errors.append(abs(fd-exact)/max(abs(exact),1e-30))
        rows.append(row)
    holdout=max(abs(row["angular_holdout_defect"]) for row in rows)
    report={"source":"https://arxiv.org/html/2609.35406v1",
            "scope":"actual angular-matched swirl with old finite core; common-pressure iteration and axial repair remain open",
            "R_core":profile.R_core,"R_anchor":profile.R_anchor,
            "anchor_power":profile.anchor_power,"collar":profile.collar.metadata(),
            "Z_count":len(rows),"rows":rows,
            "angular_independent_n256_max_defect":holdout,
            "minimum_positive_swirl_passed":positive,
            "sampled_negative_angular_shear_passed":negative_shear,
            "independent_F_R_max_relative_error":max(derivative_errors),
            "common_pressure_core_rebuilt":False,"whole_stress_cone_validated":False,
            "five_moments_closed":False,"recursion_validated":False}
    report["checks_passed"]=bool(positive and negative_shear and holdout<1e-6
                                 and max(derivative_errors)<2e-5)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('rows','collar')},indent=2))
    if not report["checks_passed"]: raise AssertionError("extended angular bridge failed")
    return report


if __name__=='__main__': run()
