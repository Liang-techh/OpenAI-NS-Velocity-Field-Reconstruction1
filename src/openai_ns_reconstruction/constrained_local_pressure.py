"""Compact pressure corrections in the shrinking velocity coordinates."""
from dataclasses import dataclass
import json
from pathlib import Path
import numpy as np
from .constrained_inner_swirl import InnerSwirlCandidate
from .constrained_temporal_swirl import _parent_payload,_parent_from_payload


@dataclass(frozen=True)
class LocalPressureCandidate:
    base: InnerSwirlCandidate
    pressure_coefficients: tuple[float,...]=(0.,)*27

    def __post_init__(self):
        if not isinstance(self.base,InnerSwirlCandidate):raise TypeError("expected inner swirl base")
        a=np.asarray(self.pressure_coefficients,dtype=float)
        if a.shape!=(27,) or not np.all(np.isfinite(a)) or np.any(np.abs(a)>1):raise ValueError('27 finite pressure coefficients in [-1,1] required')
        object.__setattr__(self,'pressure_coefficients',tuple(float(v) for v in a))

    @property
    def force(self):return self.base.force
    def velocity(self,points,time):return self.base.velocity(points,time)
    def energy(self,time=.25,order=96):return self.base.energy(time,order)

    def pressure_basis(self,points,time):
        chart=self.base.parent.base._chart(points,time)
        tau,R2,Z,b=chart[2],chart[5],chart[6],chart[7]
        w=2*(np.asarray(time)-.25)
        temporal=np.stack(((1-w)**2,2*w*(1-w),w*w),axis=-1)
        spatial=np.stack([tau**(-1.01)*b*np.exp(-((R2-r)/1.2)**2-((Z*Z-z)/1.2)**2)
                          for r in (0.,1.5,3.) for z in (0.,1.5,3.)],axis=-1)
        return (spatial[..., :,None]*temporal[...,None,:]).reshape(spatial.shape[:-1]+(27,))

    def pressure(self,points,time):
        return self.base.pressure(points,time)+self.pressure_basis(points,time)@np.asarray(self.pressure_coefficients)

    def save(self,path):
        data={'family':'local_pressure_v1','status':'candidate','pressure_coefficients':self.pressure_coefficients,
              'base':{'coefficients':self.base.coefficients,'order':self.base.order,'parent':_parent_payload(self.base.parent)}}
        Path(path).write_text(json.dumps(data,indent=2)+'\n')

    @classmethod
    def load(cls,path):
        d=json.loads(Path(path).read_text())
        if d['family']!='local_pressure_v1':raise ValueError('wrong family')
        b=d['base'];return cls(InnerSwirlCandidate(_parent_from_payload(b['parent']),tuple(b['coefficients']),b['order']),tuple(d['pressure_coefficients']))
