"""Replay spatial neighborhoods and nearby times of the local cone repair."""
import json
import numpy as np
from adaptive_bridge_recursive_defect import build_fields
from outer_pressure_modes import OuterPressure,WINDOWS
from outer_swirl_slope import OuterSwirlSlope
from separated_moment_modes import SeparatedMomentModes
from midplane_outer_residual_source import outer_cones
from midplane_integrated_moment_balance import evaluate
from midplane_remote_pressure_feasibility import P_WINDOWS,P_BREAKS
from radial_continuation import ROOT


def build_field(report_name='midplane_outer_slope_pressure_repair.json'):
    inner,fields=build_fields();base=fields['two_sided_cone']
    seed=json.loads((ROOT/'midplane_outer_pressure_staged_k11.json').read_text())
    velocity=SeparatedMomentModes(base,seed['amplitudes'],windows=WINDOWS,knots=(11.,15.,19.),axial_powers=(0,1,2))
    current=OuterPressure(velocity,seed['pressure_coefficients'])
    c=np.array(json.loads((ROOT/report_name).read_text())['coefficients'])
    return inner,OuterSwirlSlope(OuterPressure(current,c[:9],windows=P_WINDOWS),c[9:],11,P_WINDOWS)


def run():
    inner,f=build_field()
    locations=[(eta+de,.75+dy) for eta in (-.2,.2) for de in (-.005,0.,.005) for dy in (-.005,0.,.005)]
    rows=outer_cones(f,11,order=64,radial_breaks=P_BREAKS,locations=locations)
    report=dict(spatial=dict(radial_halfwidth=.005,eta_halfwidth=.005,rows=rows,pass_count=sum(r['cone_pass'] for r in rows)),times=[],
        scope='Finite spatial and time samples of one local affine swirl trajectory; not a continuum or interval certificate.',accepted=False,scale_recursion_established=False)
    print(json.dumps(dict(stage='spatial',passes=report['spatial']['pass_count'],total=len(rows),negative_lambda=sum(r['lambda_squared']<=0 for r in rows))),flush=True)
    path=ROOT/'midplane_outer_slope_neighborhood.json'
    path.write_bytes((json.dumps(report,indent=2)+'\n').encode())
    for k in (10.999,11.,11.001):
        cone=outer_cones(f,k,order=64,radial_breaks=P_BREAKS)
        moments=evaluate(f,inner,k,96,.002,radial_breaks=P_BREAKS)
        row=dict(k=k,cone=cone,pass_count=sum(r['cone_pass'] for r in cone),moments=moments.tolist(),moment_max=float(max(abs(moments))))
        report['times'].append(row);path.write_bytes((json.dumps(report,indent=2)+'\n').encode())
        print(json.dumps(dict(k=k,passes=row['pass_count'],moment_max=row['moment_max'])),flush=True)

if __name__=='__main__':run()
