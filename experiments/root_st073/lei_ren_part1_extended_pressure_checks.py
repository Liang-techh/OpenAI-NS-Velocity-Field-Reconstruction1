"""Independent radial-pressure checks for the callable extended profile."""
import json
from pathlib import Path

from lei_ren_part1_extended_pressure_core import load_extended_profile


def run():
    profile=load_extended_profile()
    rows=[]
    ra=profile.collar.collar_inner_radius
    for z in (-.5,0.,.5):
        for r in (.002,.0075,.03,1.,500.,ra*.999):
            dr=r*2e-5
            derivative=(profile.pressure(r+dr,z)-profile.pressure(r-dr,z))/(2*dr)
            density=profile.F(r,z)**2
            rows.append({'R':r,'Z':z,'pressure':profile.pressure(r,z),
                         'derivative':derivative,'F_squared':density,
                         'relative_error':abs(derivative-density)/max(abs(density),1e-30)})
    error=max(row['relative_error'] for row in rows)
    report={'source':'https://arxiv.org/html/2609.35406v1','rows':rows,
            'max_radial_pressure_relative_error':error,
            'radial_pressure_identity_passed':bool(error<2e-5),
            'pressure_definition':'Actual extended F squared integral with exact heat-tail normalization',
            'momentum_validated':False}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
    if not report['radial_pressure_identity_passed']: raise AssertionError('radial pressure identity failed')
    return report


if __name__=='__main__': run()
