"""Combine repaired knot coefficients and inspect intermediate-time dynamics."""
import json
import numpy as np
from adaptive_bridge_recursive_defect import build_fields
from adaptive_bridge_moment_fit import moment_slices,outer_moments
from affine_momentum import momentum
from midplane_resolved_feasibility import RADIAL_BREAKS,_cone_replay
from separated_moment_modes import SeparatedMomentModes,RADIAL_WINDOWS_THREE
from radial_continuation import ROOT


def run():
    inner,fields=build_fields();base=fields['two_sided_cone']
    a=np.zeros((3,2,3,3))
    for block,k in enumerate((11,15,19)):
        name='midplane_outer_axial_repair.json' if k==11 else f'midplane_outer_axial_repair_k{k}.json'
        r=json.loads((ROOT/name).read_text())
        a[block]=np.array(r['amplitudes']).reshape(3,2,3,3)[block]
    f=SeparatedMomentModes(base,a.ravel(),windows=RADIAL_WINDOWS_THREE,knots=(11.,15.,19.),axial_powers=(0,1,2))
    supports=json.loads((ROOT/'midplane_physical_covariance_pairs.json').read_text())
    report=dict(amplitudes=a.ravel().tolist(),scales=[],accepted=False,scale_recursion_established=False,
        scope='Two intermediate-time samples of one smooth coefficient-interpolated field; includes time derivatives in full momentum. Not a time-interval bound or recursive contraction.')
    for k in (13,17):
        rows=[]
        for order in (48,96):
            d=moment_slices(inner,base,f,orders=(float(k),),n=order,unit_fields=[f],radial_breaks=RADIAL_BREAKS)[0]
            moments=outer_moments(d,np.zeros(1));r=momentum(d['baseline'])
            rows.append(dict(order=order,moments=moments.tolist(),moment_max=float(max(abs(moments))),
                momentum_peak=float(max(np.linalg.norm(r,axis=1))),component_peaks=np.max(abs(r),axis=0).tolist()))
            print(json.dumps(dict(k=k,**rows[-1])),flush=True)
        cone=_cone_replay(inner,f,supports,k,order=96)
        report['scales'].append(dict(k=k,replay=rows,cone=cone))
        (ROOT/'midplane_outer_axial_interscale.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps(dict(cone_passes=[s['cone']['pass_count'] for s in report['scales']])),flush=True)

if __name__=='__main__':
    run()
