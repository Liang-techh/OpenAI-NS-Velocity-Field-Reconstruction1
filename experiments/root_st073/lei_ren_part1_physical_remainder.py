"""Measured physical R_B, D T_B and E_B for one actual extended candidate.

The selected stress envelope is B(z)^2, with the same physical cutoff as
the velocity. No residual-fitted stress or force is introduced. This leading
stress is not claimed admissible or to give a flat remainder. The completed
tensor has zero radial divergence, including in the cutoff region.
"""
import math
import numpy as np

from lei_ren_part1_extended_field import ExtendedPartIBackgroundField
from lei_ren_part1_extended_profile_stress import ExtendedProfileStressAdapter
from openai_ns_reconstruction.coordinates import similarity_coordinates_from_tau


class PhysicalRemainder:
    def __init__(self, field=None):
        self.field=ExtendedPartIBackgroundField() if field is None else field
        self.adapter=ExtendedProfileStressAdapter(self.field.matched)

    def stress_pair(self, point, tau):
        x,y,z=map(float,point)
        B,_=self.field.z_cutoff(z)
        if B==0: return np.zeros(2)
        s=similarity_coordinates_from_tau(math.hypot(x,y)/math.sqrt(self.field.nu),
                                          z/math.sqrt(self.field.nu),tau,self.field.h)
        if s.X<=0: raise ValueError('Stress diagnostic requires r>0; axis extension is separate')
        return self.field.nu*B*B*s.q**(-1-self.field.h)*np.array(self.adapter.stress(s.X,s.eta))

    def evaluate(self, point, tau, *, spatial_step, tau_step):
        p=np.asarray(point,dtype=float)
        if p.shape!=(3,) or not np.all(np.isfinite(p)):
            raise ValueError('point must have three finite coordinates')
        h,dt=float(spatial_step),float(tau_step)
        r=math.hypot(p[0],p[1])
        if not (0<h<r/4 and 0<dt<tau/4):
            raise ValueError('Require 0<spatial_step<r/4 and 0<tau_step<tau/4')
        velocity=lambda x: self.field.velocity_from_tau(*x,tau)
        pressure=lambda x: self.field.pressure_from_tau(*x,tau)
        v=velocity(p)
        gradient=np.empty((3,3)); pressure_gradient=np.empty(3); lap=np.zeros(3)
        for j in range(3):
            e=np.zeros(3); e[j]=h
            minus2,minus,plus,plus2=[velocity(p+k*e) for k in (-2,-1,1,2)]
            gradient[:,j]=(minus2-8*minus+8*plus-plus2)/(12*h)
            lap+=(-plus2+16*plus-30*v+16*minus-minus2)/(12*h*h)
            pressure_gradient[j]=(pressure(p-2*e)-8*pressure(p-e)
                                  +8*pressure(p+e)-pressure(p+2*e))/(12*h)
        samples=[self.field.velocity_from_tau(*p,tau+k*dt) for k in (-2,-1,1,2)]
        # t=T-tau, so d_t=-d_tau.
        time_derivative=-(samples[0]-8*samples[1]+8*samples[2]-samples[3])/(12*dt)
        residual=time_derivative+gradient@v+pressure_gradient-self.field.nu*lap
        er=np.array([p[0]/r,p[1]/r,0.]); et=np.array([-er[1],er[0],0.])
        pair=self.stress_pair(p,tau)
        pair_samples=[self.stress_pair(p+k*h*er,tau) for k in (-2,-1,1,2)]
        pair_r=(pair_samples[0]-8*pair_samples[1]+8*pair_samples[2]-pair_samples[3])/(12*h)
        div_cyl=np.array([0.,pair_r[0]+2*pair[0]/r,pair_r[1]+pair[1]/r])
        residual_cyl=np.array([er@residual,et@residual,residual[2]])
        remainder=residual_cyl+div_cyl
        return {'R_B':residual_cyl.tolist(),'D_T_B':div_cyl.tolist(),
                'E_B':remainder.tolist(),'divergence':float(np.trace(gradient)),
                'radial_remainder_equals_radial_residual':bool(remainder[0]==residual_cyl[0]),
                'spatial_step':h,'tau_step':dt,
                'forcing':'zero diagnostic forcing',
                'stress_envelope':'B(z)^2 times nu q^(-1-h) actual profile T',
                'whole_cone_validated':False,'flatness_established':False,
                'full_corrected_field':False}
