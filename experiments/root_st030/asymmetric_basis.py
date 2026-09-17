"""Axisymmetric candidate basis without the extra z-reflection restriction.
The physical NS, compact support, original force and numerical gates are unchanged.
Preserve nine hybrid columns; add three opposite-parity axial columns.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from hybrid_basis import HybridFamily,hybrid_basis
from spacetime import Family


def asymmetric_basis(s,z,nr,nz,odd):
    if nz<9 or nz>12 or nr>9:raise ValueError('Asymmetric family requires nr<=9, 9<=nz<=12')
    a=hybrid_basis(s,z,nr,9,odd)
    if nz==9:return a
    b=hybrid_basis(s,z,nr,nz-9,not odd)
    return {k:np.concatenate((a[k].reshape(-1,nr,9),b[k].reshape(-1,nr,nz-9)),axis=2).reshape(-1,nr*nz) for k in a}

class AsymmetricFamily(HybridFamily):
    basis=staticmethod(asymmetric_basis)
    basis_kind='hybrid9_axial_opposite3_v1'
    reflection_symmetric=False


def embed(warm,out,perturbation=.03):
    old,raw=Family.load(warm)
    if old.basis_kind!='hybrid_legendre7_inverse_even_v1' or old.nr>9 or old.nz>9:raise ValueError('Only exact preserved hybrid subspace allowed')
    f=AsymmetricFamily(9,12,12);aa,bb,qq,fc,_=old.coefficients(raw);coeffs=[]
    for co,To,Tn in [(aa,old.Tp,f.Tp),(bb,old.Tw,f.Tw),(qq,old.Tq,f.Tq)]:
        pol=(To@co.reshape(old.ns,old.nt)).reshape(old.nr,old.nz,old.nt)
        padded=np.zeros((9,12,12));padded[:old.nr,:old.nz,:old.nt]=pol
        coeffs.append(np.linalg.solve(Tn,padded.reshape(f.ns,f.nt)).ravel())
    x=np.r_[*coeffs,fc];base=x.copy()
    # A small even-in-z streamfunction component supplies axial bias. This is
    # a changed initial field, not claimed to inherit the parent's PDE result.
    p=np.zeros((9,12,12));p[0,9,0]=perturbation
    x[:f.n]+=np.linalg.solve(f.Tp,p.reshape(f.ns,f.nt)).ravel()
    x[:2*f.n]/=max(1.,np.max(abs(x[:2*f.n]))/3.9)
    Path(out).parent.mkdir(parents=True,exist_ok=True);f.save(x,out,dict(stage='independent asymmetric seed before fitting',source=str(warm),axial_bias_coefficient=perturbation))
    rng=np.random.default_rng(9172790);pts=rng.uniform(-2,2,(512,3));tt=rng.uniform(.25,.75,512)
    u,p=old.fields(raw,pts,tt);ub,pb=f.fields(base,pts,tt);uv,pv=f.fields(x,pts,tt)
    info=dict(source_basis=old.basis_kind,target_basis=f.basis_kind,embedding_velocity_error=float(np.max(abs(u-ub))),embedding_pressure_error=float(np.max(abs(p-pb))),perturbation_velocity_max=float(np.max(abs(uv-ub))),initial_energy='renormalized by representation',global_symmetry='axisymmetric retained; reflection removed',pde_validated=False)
    Path(out).with_name('embedding.json').write_text(json.dumps(info,indent=2)+'\n');print(json.dumps(info,indent=2));return f,x
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('warm');p.add_argument('out');p.add_argument('--perturbation',type=float,default=.03);a=p.parse_args();embed(a.warm,a.out,a.perturbation)
