"""Small fixed-velocity pressure feasibility repair, before fresh control audits.
Minimize the norm of a pressure-gradient correction under finite linear sign guards.
"""
from fit_guarded import Guarded
from fit_profile import Matrix,atomic,ROOT
from pathlib import Path
from scipy.optimize import minimize
import numpy as np,json,hashlib,time
from datetime import datetime,timezone

def run():
    out=ROOT/'evidence/pressure_restore';out.mkdir(exist_ok=False)
    reg=dict(utc=datetime.now(timezone.utc).isoformat(),source='ST065-V final iterate',target='ST065-C',velocity_fixed=True,source_hash=hashlib.sha256((ROOT/'candidates/ST065-V/profile.npz').read_bytes()).hexdigest(),core_samples='15 radial values .0...15, 42 heights outside |z|<=.03 up to .25',minimum_signed_axial_pressure_gradient=.002,minimum_radial_pressure_s_derivative=.001,objective='Minimum squared pressure-gradient change in original 96-order metric',maxiter=200,fresh_log_seed=9226595,fresh_point_seed=9226596)
    atomic(out/'registration.json',reg)
    o=Guarded();x=np.load(ROOT/'candidates/ST065-V/variables.npy');co,st=o.physical(x)
    rr,zz=np.meshgrid(np.linspace(0,.15,15),np.r_[np.linspace(-.25,-.03,21),np.linspace(.03,.25,21)],indexing='ij')
    M=Matrix(rr.ravel()**2,zz.ravel(),o.T)
    A=np.r_[np.sign(M.z)[:,None]*M.D['Pz'],M.D['Ps']];b=np.r_[np.full(len(M.z),.002),np.full(len(M.z),.001)]-A@co[2]
    norms=np.maximum(np.linalg.norm(A,axis=1),1e-9);AA=A/norms[:,None];bb=b/norms
    t=time.monotonic();r=minimize(lambda d:(.5*d@d,d),np.zeros(108),jac=True,method='SLSQP',constraints=[dict(type='ineq',fun=lambda d:AA@d-bb,jac=lambda d:AA)],options=dict(maxiter=200,ftol=1e-13,disp=False))
    assert np.min(AA@r.x-bb)>=-1e-8
    y=x.reshape(3,108).copy();y[2]+=r.x;o.fun(y.ravel());o.save(y.ravel(),ROOT/'candidates/ST065-C')
    ap=np.load(ROOT/'candidates/ST065-V/profile.npz');bp=np.load(ROOT/'candidates/ST065-C/profile.npz')
    for k in ['F','G']:np.testing.assert_array_equal(ap[k],bp[k])
    atomic(out/'receipt.json',dict(success=bool(r.success),message=str(r.message),iterations=int(r.nit),seconds=time.monotonic()-t,minimum_normalized_guard=float(np.min(AA@r.x-bb)),pressure_change_metric_norm=float(np.linalg.norm(r.x)),velocity_arrays_identical=True,stats=o.last,pde_validated=False))
    atomic(out/'FROZEN.json',dict(utc=datetime.now(timezone.utc).isoformat(),profile_sha256=hashlib.sha256((ROOT/'candidates/ST065-C/profile.npz').read_bytes()).hexdigest(),fresh_log_seed=9226595,fresh_point_seed=9226596,no_further_fitting=True))
    print(json.dumps(json.loads((out/'receipt.json').read_text()),indent=2))
if __name__=='__main__':run()
