"""ST032: reuse opposite-z-parity basis with exact pressure-projected curvature.
Global axisymmetry is retained. Only the autonomous reflection ansatz is removed.
"""
from __future__ import annotations
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from spacetime import Family
from asymmetric_basis import AsymmetricFamily
import mixed_exact_newton as solver
from harmonic_fit import Objective as HarmonicObjective, harmonic_tensors

class AllParityMomentObjective(HarmonicObjective):
    def __init__(self, family, seed=9172902, count=6144, moment_weight=10.):
        super().__init__(family, seed, count, moment_weight)
        mp, mb, self.harmonic_gram=harmonic_tensors(family, degrees=tuple(range(2,9)))
        self.Hp=self.torch.tensor(mp);self.Hb=self.torch.tensor(mb)


def embed(warm, target, bias=.003):
    old,x=Family.load(warm)
    if old.basis_kind!='hybrid_legendre7_inverse_even_v1' or old.nr!=9 or old.nz!=9 or old.nt!=8:
        raise ValueError('Only frozen 9x9x8 parent is accepted')
    f=AsymmetricFamily(9,12,8);av,bv,pv,fc,_=old.coefficients(x);parts=[]
    for c,O,T in [(av,old.Tp,f.Tp),(bv,old.Tw,f.Tw),(pv,old.Tq,f.Tq)]:
        poly=(O@c.reshape(old.ns,8)).reshape(9,9,8)
        pad=np.zeros((9,12,8));pad[:,:9]=poly
        parts.append(np.linalg.solve(T,pad.reshape(f.ns,8)).ravel())
    raw=np.r_[*parts,fc];base=raw.copy()
    rng=np.random.default_rng(9172903);pts=rng.uniform(-2,2,(512,3));ts=rng.uniform(.25,.75,512)
    u,p=old.fields(x,pts,ts);ub,pb=f.fields(base,pts,ts)
    pol=np.zeros((9,12,8));pol[0,9,0]=bias
    raw[:f.n]+=np.linalg.solve(f.Tp,pol.reshape(f.ns,8)).ravel()
    raw[:2*f.n]/=max(1.,float(np.max(np.abs(raw[:2*f.n])))/3.9)
    f.save(raw,target,{'scope':'changed asymmetric initialization, not inherited PDE success','bias':bias})
    check={'velocity_embedding_max_error':float(np.max(np.abs(u-ub))), 'pressure_embedding_max_error':float(np.max(np.abs(p-pb))), 'bias':bias,'stored_parameters':len(raw),'global_axisymmetry_retained':True,'pde_validated':False}
    Path(target).with_name('embedding.json').write_text(json.dumps(check,indent=2)+'\n')
    return check


def main():
    p=argparse.ArgumentParser();p.add_argument('--warm',required=True);p.add_argument('--out',required=True);p.add_argument('--iterations',type=int,default=12);a=p.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    reg={'experiment':'CR-ROOT-ST032','created_before_fit':True,'initial_sha256':hashlib.sha256(Path(a.warm).read_bytes()).hexdigest(),'basis':'existing Hybrid9+3 opposite axial parity columns; nt=8 unchanged','stored_parameters':2594,'global_axisymmetry':True,'reflection_symmetry':False,'physical_gates':'Original nu, support, force family/bounds, energy, core and 1e-3 max/L2 gates unchanged','bias':.003,'training_seed':9172902,'validation_seed':9172911,'iterations':a.iterations,'training_count':4096,'method':'same exact mixed-collocation projected Hessian optimizer as ST030','harmonic_degrees':list(range(2,9)),'not_claimed':'source profile reconstruction or OpenAI field identity'}
    (out/'scope_registration.json').write_text(json.dumps(reg,indent=2)+'\n')
    print(json.dumps(embed(a.warm,out/'asymmetric_initial.json'),indent=2),flush=True)
    # Keep generic optimizer registration intact; experiment-level scope is above.
    solver.Objective=AllParityMomentObjective
    solver.run(out/'asymmetric_initial.json',out/'fit',a.iterations,seed=9172902,count=4096,weight=10.)

if __name__=='__main__':main()
