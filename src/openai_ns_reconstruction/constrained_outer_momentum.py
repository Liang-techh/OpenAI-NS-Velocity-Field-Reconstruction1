"""Outer swirl correction constrained by prescribed global torque, not residual force."""
from dataclasses import asdict
import json
from pathlib import Path
import numpy as np
from .cutoffs import standard_cutoff
from .quadrature import unit_rule
from .constrained_force import RestrictedForce,compact_bump
from .constrained_tensor_candidate import TensorCandidate
from .constrained_momentum_budget import angular_moment


class AngularMomentumCandidate:
    def __init__(self,base,force,order=96):
        self.base,self.force,self.order=base,force,order
        ts=np.array([.25,.5,.75])
        values=np.array([angular_moment(base.velocity,float(t),order) for t in ts])
        # Fixed geometry and degree-two time tensor imply J=tau^1.49 P2(t).
        self.coefficients=np.polynomial.polynomial.polyfit(ts-.25,values/(1-ts)**1.49,2)
        self.initial_moment=values[0]
        self.basis_moment=angular_moment(self.outer_basis,.5,order)
        self.force_moment=angular_moment(force,.5,order)
        if self.basis_moment<=0:raise ValueError('nonpositive outer moment')

    @staticmethod
    def outer_basis(points,time):
        p=np.asarray(points,dtype=float);x,y,z=np.moveaxis(p,-1,0)
        r2=x*x+y*y
        def cutoff(a):
            return np.fromiter((standard_cutoff(float(v)) for v in a.ravel()),float,count=a.size).reshape(a.shape)
        b=cutoff(r2/4)*cutoff(z*z/4)*(1-cutoff(r2/.25))
        return np.stack((-y*b,x*b,np.zeros_like(b)),axis=-1)

    def amplitude(self,time):
        t=np.asarray(time,dtype=float);nodes,weights=unit_rule(32)
        q=.25+(t[...,None]-.25)*nodes
        g,_=compact_bump((2*q-1)**2)
        integral=(t-.25)*np.sum(g*weights,axis=-1)
        target=self.initial_moment+self.force_moment*integral
        base=(1-t)**1.49*np.polynomial.polynomial.polyval(t-.25,self.coefficients)
        return (target-base)/self.basis_moment

    def velocity(self,points,time):
        u=self.base.velocity(points,time)
        return u+self.amplitude(time)[...,None]*self.outer_basis(points,time)

    def pressure(self,points,time):return self.base.pressure(points,time)

    def energy(self,time=.25,order=96):
        n,w=unit_rule(order);r=2*n;z=4*n-2
        rr,zz=np.meshgrid(r,z,indexing='ij')
        p=np.stack((rr,np.zeros_like(rr),zz),axis=-1)
        return float(8*np.pi*np.sum(rr*np.sum(self.velocity(p,time)**2,axis=-1)*w[:,None]*w[None,:]))

    def save(self,path):
        data={'family':'angular_momentum_outer_v1','schema_version':1,'status':'candidate',
              'base_parameters':asdict(self.base),'force':asdict(self.force),'quadrature_order':self.order}
        Path(path).write_text(json.dumps(data,indent=2)+'\n')

    @classmethod
    def load(cls,path):
        d=json.loads(Path(path).read_text())
        if d['family']!='angular_momentum_outer_v1':raise ValueError('wrong family')
        return cls(TensorCandidate(**d['base_parameters']),RestrictedForce(**d['force']),d['quadrature_order'])


if __name__=='__main__':
    root=Path('artifacts/constrained/tensor_feasible')
    base=TensorCandidate.load(root/'candidate.json')
    f=RestrictedForce(**json.loads((root/'training.json').read_text())['force'])
    c=AngularMomentumCandidate(base,f)
    out=Path('artifacts/constrained/outer_momentum');out.mkdir(exist_ok=True)
    c.save(out/'candidate.json')
    (out/'training.json').write_text(json.dumps({'status':'deterministic_moment_correction_not_optimized','force':asdict(f)},indent=2)+'\n')
    print('outer moment correction saved')
