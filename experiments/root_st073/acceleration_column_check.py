"""Physical FD regression for imaginary acceleration and pressure columns."""
import json
from pathlib import Path
import numpy as np
from affine_momentum import jets,momentum
from endpoint_acceleration_projection import (
    _acceleration_design,_unpack_acceleration_control,AccelerationCorrectionField)

ROOT=Path(__file__).resolve().parent


class ConstantBase:
    nu=.01
    def fields(self,points,tau):
        return np.broadcast_to([1.,-2.,.5],(len(points),3)).copy(),np.zeros(len(points))


def run():
    snapshot=json.loads((ROOT/'full_wave_frozen_cache.json').read_text())
    g=snapshot['inputs']['wave'];center=np.asarray(g['center']);width=np.asarray(g['widths'])
    angles=np.array([.2,.7,1.2,2.1]);radius=center[0]+width[0]*np.array([-.3,.1,.25,-.2])
    points=np.column_stack((radius*np.cos(angles),radius*np.sin(angles),center[1]+width[1]*np.array([.2,-.3,.1,-.1])))
    tau0=snapshot['inputs']['mean']['tau'];tau1=tau0*2**(-1e-6);dt=tau0-tau1
    base=ConstantBase();u,_=base.fields(points,tau1)
    D,_=_acceleration_design(points,center,width,np.asarray(g['carrier']),dt,u,np.zeros((4,3,3)),base.nu)
    rows=[]
    for index,amplitude in ((37,1e17),(91,1e15),(109,1e17),(253,1e17)):
        c=np.zeros(324);c[index]=amplitude
        results=[]
        for sign in (1.,-1.):
            blocks=_unpack_acceleration_control(sign*c)
            field=AccelerationCorrectionField(base,center,width,g['carrier'],
                {m:v[0] for m,v in blocks.items()},{m:v[1] for m,v in blocks.items()},tau0)
            results.append(momentum(jets(field,points,tau1,2e-7,1e-9)))
        observed=.5*(results[0]-results[1])
        expected=(D@c).reshape(-1,3)
        error=float(np.linalg.norm(observed-expected)/np.linalg.norm(expected))
        if error>1e-5:
            raise AssertionError(f'Imaginary column {index} fails physical derivative check: {error}')
        rows.append(dict(column=index,relative_error=error,expected_norm=float(np.linalg.norm(expected))))
    report=dict(status='completed',scope='Central control perturbations cancel quadratic convection; compares corrected imaginary columns with actual finite-difference momentum.',rows=rows)
    (ROOT/'acceleration_column_check.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    run()
