"""Compact axisymmetric candidate, not an independently validated NS solution.

Streamfunction psi = a*r^2*tau^(-A)*Z*b(R,Z)*q(R,Z).
u_r=-psi_z/r, u_z=psi_r/r, u_theta=a*c*tau^(-A)*R*b.
The Cartesian formulas extend smoothly across the axis without division by r.
"""
from dataclasses import asdict, dataclass, replace
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from .constrained_force import compact_bump


@dataclass(frozen=True)
class CompactCandidate:
    swirl_ratio: float = 1.0
    radial_width: float = 1.0
    axial_width: float = 1.0
    radial_shape: float = 0.0
    axial_shape: float = 0.0
    swirl_radial_shape: float = 0.0
    swirl_axial_shape: float = 0.0
    swirl_radial_time: float = 0.0
    swirl_axial_time: float = 0.0
    amplitude: float = 1.0
    pressure_constant: float = 0.0
    pressure_radial: float = 0.0
    pressure_axial: float = 0.0

    def __post_init__(self):
        bounds = {
            'swirl_ratio': (0.25, 4), 'radial_width': (0.7, 1.1),
            'axial_width': (0.7, 1.1), 'radial_shape': (-0.2, 0.2),
            'axial_shape': (-0.2, 0.2), 'swirl_radial_shape': (-0.8,0.8),
            'swirl_axial_shape': (-0.8,0.8),
            'swirl_radial_time': (-4,4), 'swirl_axial_time': (-4,4), 'amplitude': (1e-4, 100),
            'pressure_constant': (-100, 100), 'pressure_radial': (-100, 100),
            'pressure_axial': (-100, 100),
        }
        for name, (lo, hi) in bounds.items():
            value = getattr(self, name)
            if not np.isfinite(value) or not lo <= value <= hi:
                raise ValueError(f'{name} must be finite and in [{lo},{hi}]')

    def _chart(self, points, time):
        points = np.asarray(points, dtype=float)
        time = np.asarray(time, dtype=float)
        if points.ndim < 1 or points.shape[-1] != 3:
            raise ValueError('points must have shape (...,3)')
        if not np.all(np.isfinite(points)) or not np.all(np.isfinite(time)):
            raise ValueError('coordinates and time must be finite')
        if np.any((time < 0.25) | (time > 0.75)):
            raise ValueError('candidate time domain is [0.25,0.75]')
        x, y, z, t = np.broadcast_arrays(*np.moveaxis(points, -1, 0), time)
        tau = 1-t
        lr = self.radial_width*np.sqrt(tau)
        lz = self.axial_width*tau**0.495
        R2, Z = (x*x+y*y)/lr**2, z/lz
        br, dbr = compact_bump(R2/4)
        bz, dbz = compact_bump(Z*Z/4)
        b = br*bz
        bR_over_R = dbr*bz/2
        bZ = br*dbz*Z/2
        q = 1+self.radial_shape*R2/4+self.axial_shape*Z*Z/4
        return x, y, tau, lr, lz, R2, Z, b, bR_over_R, bZ, q

    def velocity(self, points, time):
        x,y,tau,lr,lz,R2,Z,b,bR,bZ,q = self._chart(points,time)
        v = self.amplitude*tau**(-0.505)
        # Cartesian u_r/r, u_theta/r and u_z; no removable singularity.
        radial = -v/lz*(b*q+Z*(bZ*q+b*self.axial_shape*Z/2))
        swirl = v*self.swirl_ratio*b/lr*(1+self.swirl_radial_shape*R2/4+self.swirl_axial_shape*Z*Z/4
            +(0.75-tau)*(self.swirl_radial_time*(R2-0.01)/4
                         +self.swirl_axial_time*(Z*Z-0.01)/4))
        axial = v*Z*(2*b*q+R2*(bR*q+b*self.radial_shape/2))
        return np.stack((x*radial-y*swirl, y*radial+x*swirl, axial),axis=-1)

    def vector_potential(self, points, time):
        """Axis-regular potential for the legacy curl/correction composition API."""
        x,y,tau,_,_,_,Z,b,_,_,q = self._chart(points,time)
        H=self.amplitude*tau**(-0.505)*Z*b*q
        return np.stack((-y*H,x*H,np.zeros_like(H)),axis=-1)

    def as_legacy_local_field(self):
        """Reuse LocalField's numerical curl, without returning the direct velocity.

        This scalar adapter enables legacy potential corrections and independently
        checks the direct poloidal formula. It is slower than batched evaluation.
        """
        from .local_field import LocalField
        def potential(x,y,z,t):
            return self.vector_potential([x,y,z],t)
        def swirl(x,y,z,t):
            r=np.hypot(x,y)
            if r==0:
                return 0.0
            u=self.velocity([x,y,z],t)
            return float((-y*u[0]+x*u[1])/r)
        return LocalField(potential,swirl)

    def pressure(self, points, time):
        _,_,tau,_,_,R2,Z,b,_,_,_ = self._chart(points,time)
        return tau**(-1.01)*b*(self.pressure_constant
            +self.pressure_radial*R2+self.pressure_axial*Z*Z)

    def energy(self, time=0.25, order=48):
        """Axisymmetric Gauss quadrature, integral over the full compact support."""
        if not isinstance(order, int) or order < 4:
            raise ValueError('quadrature order must be an integer >=4')
        if not np.isfinite(time) or not 0.25 <= time <= 0.75:
            raise ValueError('candidate time domain is [0.25,0.75]')
        nodes, weights = leggauss(order)
        lr = self.radial_width*np.sqrt(1-time)
        lz = self.axial_width*(1-time)**0.495
        r = lr*(nodes+1)
        z = 2*lz*nodes
        rr,zz = np.meshgrid(r,z,indexing='ij')
        points = np.stack((rr,np.zeros_like(rr),zz),axis=-1)
        speed2 = np.sum(self.velocity(points,time)**2,axis=-1)
        return float(np.pi*np.sum(speed2*rr*weights[:,None]*weights[None,:])*lr*2*lz)

    def normalized(self, order=96):
        """Fix E(0.25)=1 numerically; amplitude is derived, not optimized."""
        e = self.energy(order=order)
        if not np.isfinite(e) or e <= 0:
            raise ValueError('cannot normalize a zero/nonfinite candidate')
        return replace(self, amplitude=self.amplitude/np.sqrt(e))

    def save(self, path):
        Path(path).write_text(json.dumps({
            'schema_version': 1, 'status': 'candidate', 'paper_exact': False,
            'family': 'compact_axisymmetric_window_v3', 'parameters': asdict(self)
        },indent=2)+'\n',encoding='utf-8')

    @classmethod
    def load(cls, path):
        data=json.loads(Path(path).read_text(encoding='utf-8'))
        if data.get('schema_version') != 1 or data.get('family') not in ('compact_axisymmetric_window_v1','compact_axisymmetric_window_v2','compact_axisymmetric_window_v3'):
            raise ValueError('unsupported candidate artifact')
        return cls(**data['parameters'])
