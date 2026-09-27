"""Factorized nonlinear momentum for a moving meridional coefficient state.

At a fixed evaluation time, retain complete base/value jets and the slope and
pressure columns. Updating coefficients then includes every quadratic
advection product without materializing an O(mode_count**2) tensor. A different
evaluation time requires a new cache; this is not a frozen-in-time NS solver.
"""
import json
from pathlib import Path
import time
import numpy as np

from affine_momentum import jets, momentum
from broad_meridional_momentum import MeridionalSlope, CompactPressureDirection
from broad_meridional_momentum import build_compact_baseline, build_modes_with_parameters
from broad_shear_dynamic_control import BroadShearSlope, load_saved_field
from outer_feedback_evolution import SwirlValue
from outer_swirl_slope import OuterSwirlSlope
from outer_pressure_modes import OuterPressure


def value_and_pressure_modes(dynamic, source):
    """Keep the original time-dependent basis, including its scale knots."""
    base=build_compact_baseline(dynamic,source)
    modes,names=build_modes_with_parameters(dynamic,source,source['k'],compact_pressure=base)
    values=[];pressures=[];value_names=[];pressure_names=[]
    for mode,name in zip(modes,names):
        if type(mode) is MeridionalSlope:
            unit=mode.unit
        elif type(mode) is BroadShearSlope:
            unit=mode.base
        elif type(mode) is OuterSwirlSlope:
            unit=SwirlValue(mode.base,mode.a)
            if tuple(unit.windows)!=tuple(mode.windows):
                raise ValueError('Swirl window mismatch')
        elif isinstance(mode,(OuterPressure,CompactPressureDirection)):
            pressures.append(mode);pressure_names.append(name)
            continue
        else:
            raise TypeError(type(mode).__name__)
        values.append(unit);value_names.append(name)
    return base,values,pressures,value_names,pressure_names


class StatefulMean:
    """Actual field for a local state a(k)=a0+(k-kref)*slope."""
    def __init__(self,base,values,pressures,state,slope,pressure,kref):
        self.base,self.values,self.pressures=base,tuple(values),tuple(pressures)
        self.state,self.slope,self.pressure=map(lambda x:np.asarray(x,float),(state,slope,pressure))
        self.kref=float(kref)
        if self.state.shape!=(len(values),) or self.slope.shape!=self.state.shape or self.pressure.shape!=(len(pressures),):
            raise ValueError('State, slope or pressure dimension mismatch')
        for name in ('inner','nu','join_X','ratio'):setattr(self,name,getattr(base,name))

    def fields(self,points,tau):
        points=np.asarray(points,float)
        times=np.asarray(tau,float).ravel()
        if not len(times) or np.any(times!=times[0]):
            raise ValueError('StatefulMean requires one common remaining time')
        tau=float(times[0])
        amplitude=self.state+(-np.log2(2*tau)-self.kref)*self.slope
        u,p=self.base.fields(points,tau)
        u,p=u.copy(),p.copy()
        for coefficient,unit in zip(amplitude,self.values):
            if coefficient: u+=coefficient*unit.fields(points,tau)[0]
        for coefficient,unit in zip(self.pressure,self.pressures):
            if coefficient: p+=coefficient*unit.fields(points,tau)[1]
        return u,p


class StateCache:
    """Full finite-difference polynomial at exactly one k and spatial grid."""
    def __init__(self,base,values,pressures,points,k,hspace=None,htime=None):
        self.k=float(k);self.tau=.5*2.**(-self.k)
        self.points=np.asarray(points,float)
        self.hspace=5e-4*np.sqrt(base.nu*self.tau) if hspace is None else float(hspace)
        self.htime=1e-4*self.tau if htime is None else float(htime)
        args=(self.points,self.tau,self.hspace,self.htime)
        self.base=jets(base,*args)
        value_jets=[jets(unit,*args) for unit in values]
        self.values=tuple(np.stack([jet[i] for jet in value_jets]) for i in range(3))
        # Exact same temporal stencil as jets(), with zero instantaneous
        # velocity/spatial jets for the local (k-kref)*unit direction.
        slope_columns=[]
        for unit in values:
            samples=[]
            for j in (-2,-1,1,2):
                t=self.tau+j*self.htime
                samples.append((-np.log2(2*t)-self.k)*unit.fields(self.points,t)[0])
            um2,um,up,up2=samples
            slope_columns.append(-(um2-8*um+8*up-up2)/(12*self.htime))
        self.slopes=np.stack(slope_columns)
        self.pressure=np.stack([jets(unit,*args)[2] for unit in pressures])

    def residual(self,state,slope,pressure):
        state,slope,pressure=map(lambda x:np.asarray(x,float),(state,slope,pressure))
        u=self.base[0]+np.tensordot(state,self.values[0],axes=1)
        grad=self.base[1]+np.tensordot(state,self.values[1],axes=1)
        linear=self.base[2]+np.tensordot(state,self.values[2],axes=1)
        linear+=np.tensordot(slope,self.slopes,axes=1)+np.tensordot(pressure,self.pressure,axes=1)
        return linear+np.einsum('nij,nj->ni',grad,u)

    def derivative_problem(self,state):
        """Offset and (point,component,column) array for slope/pressure fitting."""
        return (self.residual(state,np.zeros(len(self.slopes)),np.zeros(len(self.pressure))),
                np.concatenate((self.slopes,self.pressure),axis=0).transpose(1,2,0))


def run():
    dynamic,source=load_saved_field()
    base,values,pressures,vnames,pnames=value_and_pressure_modes(dynamic,source)
    k=source['k'];tau=source['tau']
    ys=np.array([.067,.24,.53,.825,.963])
    points=dynamic.inner.from_similarity(dynamic.join_X*(1+15*ys)**2,
        np.array([-.17,.11,.23,-.29,.07]),tau,np.array([.31,.72,1.2,.47,.91]))
    started=time.perf_counter();cache=StateCache(base,values,pressures,points,k)
    assembly_seconds=time.perf_counter()-started
    rng=np.random.default_rng(731)
    rows=[]
    for size in (.02,.2):
        state=size*rng.normal(size=len(values))
        slope=2*rng.normal(size=len(values))
        pressure=.01*rng.normal(size=len(pressures))
        field=StatefulMean(base,values,pressures,state,slope,pressure,k)
        predicted=cache.residual(state,slope,pressure)
        direct=momentum(jets(field,points,tau,cache.hspace,cache.htime))
        error=np.linalg.norm(predicted-direct,axis=1)
        delta_u=np.tensordot(state,cache.values[0],axes=1)
        delta_j=np.tensordot(state,cache.values[1],axes=1)
        nonlinear=np.einsum('nij,nj->ni',delta_j,delta_u)
        offset,columns=cache.derivative_problem(state)
        affine=offset+np.einsum('ncp,p->nc',columns,np.r_[slope,pressure])
        rows.append(dict(state_scale=size,state=state.tolist(),slope=slope.tolist(),pressure=pressure.tolist(),
            residual_max=float(np.linalg.norm(direct,axis=1).max()),
            cache_direct_absolute_error_max=float(error.max()),
            cache_direct_scaled_error_max=float(np.max(error/np.maximum(np.linalg.norm(direct,axis=1),1.))),
            retained_quadratic_term_max=float(np.linalg.norm(nonlinear,axis=1).max()),
            affine_derivative_problem_error=float(np.max(abs(affine-predicted)))))
    result=dict(accepted=False,pde_validated=False,scale_recursion_established=False,
        k=k,tau=tau,points=points.tolist(),velocity_state_count=len(values),
        slope_count=len(values),pressure_count=len(pressures),
        value_names=vnames,pressure_names=pnames,assembly_seconds=assembly_seconds,
        rows=rows,
        scope='Five-point nonlinear state/derivative cache consistency at one k, two nonzero states. Includes all quadratic convection. No optimized candidate, time integration, moment/cone closure, PDE or recursion acceptance. Reassemble for a different evaluation k.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:value for key,value in result.items() if key not in ('value_names','pressure_names','points','rows')}),flush=True)
    print(json.dumps(rows),flush=True)
    return result


if __name__=='__main__':run()
