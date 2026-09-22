"""Smooth parameter continuation of five exact moment equations.
Horizontal lift a_eta=-J^+(partial_eta residual), not pointwise branch switching.
The independent audits occur only after a coefficient polynomial is frozen.
"""
from pathlib import Path
import json,time,hashlib
from datetime import datetime,timezone
import numpy as np
from scipy.integrate import solve_ivp
from numpy.polynomial import chebyshev as ch
from annulus import Annulus
ROOT=Path(__file__).parent

def main():
    s=Annulus(4.,.25,n=40);a0,rec=s.fit([-.5]);a0=a0[0]
    reg=dict(stage='ST069-F2',registered_utc=datetime.now(timezone.utc).isoformat(),Xb=4.,heat_c=.25,seam='C-infinity flat blend with entire finite core polynomial, not merely C2 jets',eta_nodes=[17,33,65,129,193],
        method='minimum-norm horizontal continuation ODE; all12 profile parameters, no pointwise minimum-energy branch switches',
        eta_partial_step=1e-4,eta_ode_rtol=2e-10,eta_ode_atol=2e-12,selection='first calibrated moment defect <1e-8',independent_seed=9226995)
    (ROOT/'evidence/interpolation_registration.json').write_text(json.dumps(reg,indent=2)+'\n')
    counter=[0];conditions=[]
    def rhs(eta,a):
        r,J=s.moments(a,eta,jac=True);step=1e-4
        de=(s.moments(a,eta-2*step)[0]-8*s.moments(a,eta-step)[0]+8*s.moments(a,eta+step)[0]-s.moments(a,eta+2*step)[0])/(12*step)
        # Stabilization damps accumulated equation error without changing exact r=0 solutions.
        d=-J.T@np.linalg.solve(J@J.T,de+5*r)
        counter[0]+=1;conditions.append(np.linalg.cond(J));return d
    t=time.monotonic();sol=solve_ivp(rhs,[-.5,.5],a0,method='DOP853',rtol=2e-10,atol=2e-12,dense_output=True,max_step=.025)
    if not sol.success:raise ValueError(sol.message)
    np.savez_compressed(ROOT/'data/continuation_nodes.npz',eta=sol.t,coefficients=sol.y)
    target=None;models=[]
    for N in [17,33,65,129,193]:
        etas=np.sort(.5*np.cos(np.arange(N)*np.pi/(N-1)));coeff=sol.sol(etas).T
        cheb=ch.chebfit(2*etas,coeff,N-1);cal=np.linspace(-.5,.5,2*N+8);errors=[]
        for eta in cal:errors.append(float(abs(s.moments(ch.chebval(2*eta,cheb),eta)[0]).max()))
        models.append(dict(nodes=N,calibration_max=max(errors)))
        print(N,max(errors),'rhs',counter[0],time.monotonic()-t,flush=True)
        path=ROOT/f'data/annulus_N{N}.npz';np.savez_compressed(path,cheb=cheb,Xc=s.Xc,Xb=s.Xb,c=s.c,h=s.h)
        if max(errors)<1e-8:
            target=ROOT/'data/ST069-M.npz';target.write_bytes(path.read_bytes());break
    (ROOT/'evidence/fitting.json').write_text(json.dumps(dict(models=models,ode_evaluations=counter[0],max_condition=max(conditions),elapsed=time.monotonic()-t,success=sol.success),indent=2)+'\n')
    if target is None:raise ValueError('No calibrated coefficient polynomial passed. No model promoted.')
    freeze=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),id='ST069-M',model_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),core_sha256=hashlib.sha256((ROOT/'data/ST068-I.npz').read_bytes()).hexdigest(),pde_validated=False,stress_cone_admitted=False,global_field_ready=False,independent_seed=9226995)
    (ROOT/'evidence/freeze.json').write_text(json.dumps(freeze,indent=2)+'\n')
if __name__=='__main__':main()
