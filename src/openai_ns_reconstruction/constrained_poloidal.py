"""Joint compact poloidal streamfunction and pressure candidate."""
from dataclasses import dataclass
from typing import ClassVar
import json
from pathlib import Path
import numpy as np
from .constrained_local_pressure import LocalPressureCandidate
from .constrained_inner_swirl import InnerSwirlCandidate
from .constrained_force import compact_bump
from .constrained_temporal_swirl import _broadcast_physical_coordinates,_bernstein_nonzero,_parent_payload,_parent_from_payload
from .quadrature import unit_rule


@dataclass(frozen=True)
class PoloidalCandidate(LocalPressureCandidate):
    coefficients: tuple[float,...]=(0.,)*27
    anchor_core: bool=False
    COEFFICIENT_COUNT: ClassVar[int]=27
    FAMILY_ID: ClassVar[str]="poloidal_joint_v1"

    def __post_init__(self):
        super().__post_init__()
        a=np.asarray(self.coefficients,dtype=float)
        if a.shape!=(self.COEFFICIENT_COUNT,) or not np.all(np.isfinite(a)) or np.any(np.abs(a)>1):raise ValueError(f'{self.COEFFICIENT_COUNT} bounded velocity coefficients required')
        object.__setattr__(self,'coefficients',tuple(float(v) for v in a))

    def correction_basis(self,points,time):
        x,y,z,t=_broadcast_physical_coordinates(points,time);s=x*x+y*y;v=z*z
        br,dr=compact_bump(s/4);bz,dz=compact_bump(v/4)
        if self.anchor_core:
            core=.01*(1-t);distance=s-core;denominator=s+core+.03
            guard=distance**2/denominator**2
            gs=2*distance/denominator**2-2*distance**2/denominator**3
        else:
            safe=np.maximum(s,.01);guard,dguard=compact_bump(.01/safe);gs=dguard*(-.01/safe**2)
        modes=[]
        for center in (.2,.7,1.2):
            for axial in (0.,.6,1.5):
                g=np.exp(-((s-center)/.4)**2-((v-axial)/.6)**2)
                F=br*bz*guard*g
                Fs=bz*g*(dr*guard/4+br*gs-2*(s-center)/.4**2*br*guard)
                Fv=br*guard*g*(dz/4-2*(v-axial)/.6**2*bz)
                radial=-(F+2*v*Fv)
                modes.append(np.stack((x*radial,y*radial,2*z*(F+s*Fs)),axis=-1))
        spatial=np.stack(modes,axis=-1)
        temporal=_bernstein_nonzero(2*(t-.25))
        return (spatial[..., :, :,None]*temporal[...,None,None,:]).reshape(x.shape+(3,27))

    def velocity(self,points,time):
        return self.base.velocity(points,time)+self.correction_basis(points,time)@np.asarray(self.coefficients)

    def energy(self,time=.25,order=96):
        n,w=unit_rule(order);r,z=np.meshgrid(2*n,4*n-2,indexing='ij')
        u=self.velocity(np.stack((r,np.zeros_like(r),z),axis=-1),time)
        return float(8*np.pi*np.sum(r*np.sum(u*u,axis=-1)*w[:,None]*w[None,:]))

    def save(self,path):
        data={'family':self.FAMILY_ID,'status':'candidate','anchor_core':self.anchor_core,'coefficients':self.coefficients,'pressure_coefficients':self.pressure_coefficients,
              'base':{'coefficients':self.base.coefficients,'order':self.base.order,'parent':_parent_payload(self.base.parent)}}
        Path(path).write_text(json.dumps(data,indent=2)+'\n')

    @classmethod
    def load(cls,path):
        d=json.loads(Path(path).read_text());b=d['base']
        if d['family']!=cls.FAMILY_ID:raise ValueError('wrong family')
        return cls(InnerSwirlCandidate(_parent_from_payload(b['parent']),tuple(b['coefficients']),b['order']),tuple(d['pressure_coefficients']),tuple(d['coefficients']),d.get('anchor_core',False))
