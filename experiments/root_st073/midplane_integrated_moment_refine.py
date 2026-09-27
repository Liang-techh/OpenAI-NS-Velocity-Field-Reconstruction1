"""Refine local k=17 value/slope correction with integrated moment identities."""
import json
import numpy as np
from adaptive_bridge_recursive_defect import build_fields
from midplane_integrated_moment_balance import evaluate
from midplane_outer_moment_dae import LocalTimeCorrection
from midplane_resolved_feasibility import _cone_replay
from separated_moment_modes import SeparatedMomentModes,RADIAL_WINDOWS_THREE
from radial_continuation import ROOT


def run():
    inner,fields=build_fields();base=fields['two_sided_cone'];k=17
    a=np.array(json.loads((ROOT/'midplane_outer_axial_interscale.json').read_text())['amplitudes'])
    current=SeparatedMomentModes(base,a,windows=RADIAL_WINDOWS_THREE,knots=(11.,15.,19.),axial_powers=(0,1,2))
    old=json.loads((ROOT/'midplane_outer_moment_dae_k17.json').read_text())
    values=np.array(old['coefficient_values']);slopes=np.array(old['coefficient_slopes'])
    def field(v,s=slopes):return LocalTimeCorrection(current,v,s,k)
    initial=evaluate(field(values),inner,k,96,.002)
    columns=[];step=1e-3
    for index in (6,7):
        e=np.eye(12)[index]*step
        columns.append((evaluate(field(values+e),inner,k,48,.002)-evaluate(field(values-e),inner,k,48,.002))/(2*step))
    J=np.column_stack(columns);delta=np.linalg.solve(J[[1,3]],-initial[[1,3]])
    values[[6,7]]+=delta
    after_value=evaluate(field(values),inner,k,96,.002)
    D=np.array(old['slope_response']);ds=np.linalg.lstsq(D[[0,2],:6],-after_value[[0,2]],rcond=None)[0]
    slopes[:6]+=ds;f=field(values,slopes)
    print(json.dumps(dict(value_delta=delta.tolist(),slope_delta_norm=float(np.linalg.norm(ds)),initial=initial.tolist(),after_value=after_value.tolist())),flush=True)
    replay=[]
    for n,zstep in ((48,.002),(96,.002),(96,.001)):
        m=evaluate(f,inner,k,n,zstep)
        replay.append(dict(order=n,z_step_factor=zstep,moments=m.tolist(),absolute_max=float(max(abs(m)))))
        print(json.dumps(replay[-1]),flush=True)
    supports=json.loads((ROOT/'midplane_physical_covariance_pairs.json').read_text())
    cone=_cone_replay(inner,f,supports,k,order=96)
    report=dict(k=k,initial_moments=initial.tolist(),value_jacobian=J.tolist(),value_delta=delta.tolist(),
        slope_delta=ds.tolist(),coefficient_values=values.tolist(),coefficient_slopes=slopes.tolist(),replay=replay,cone=cone,
        scope='Refined local affine coefficient trajectory using integrated moment identities. No full momentum improvement, continuous DAE evolution, or recursion established.',accepted=False,scale_recursion_established=False)
    (ROOT/'midplane_integrated_moment_refine.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps(dict(cone_pass=cone['pass_count'])),flush=True)

if __name__=='__main__':
    run()
