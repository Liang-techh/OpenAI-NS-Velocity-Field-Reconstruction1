"""Initial-time momentum columns without redundant zero-velocity stencils.

Only the registered pressure-only and zero-at-k0 slope modes are supported.
This evaluates derivatives of the actual mode formulas, not a fitted operator.
It cannot replace the nonlinear velocity jets after a state has evolved.
"""
import json
from pathlib import Path
import time
import numpy as np

from affine_momentum import jets, momentum
from broad_meridional_momentum import MeridionalSlope, CompactPressureDirection
from broad_meridional_momentum import build_compact_baseline, build_modes_with_parameters
from broad_shear_dynamic_control import BroadShearSlope, load_saved_field
from midplane_resolved_feasibility import ZeroBackground
from outer_swirl_slope import OuterSwirlSlope
from outer_pressure_modes import OuterPressure
from outer_feedback_evolution import SwirlValue
from separated_moment_modes import flat_bump
from joined_field import coordinates


def pressure_gradient(mode, points, tau):
    """Exact Cartesian gradient of the registered scalar bump pressure."""
    points=np.asarray(points,float)
    radius=np.hypot(points[:,0],points[:,1])
    safe=np.where(radius>0,radius,1.)
    sn=np.sqrt(mode.nu)
    co=coordinates(radius/sn,points[:,2]/sn,tau,mode.inner.h)
    q,eta=np.asarray(co['q']),np.asarray(co['eta'])
    qz,etaz=np.asarray(co['q_z'])/sn,np.asarray(co['eta_z'])/sn
    ri=np.sqrt(2*mode.nu*q*mode.join_X)
    width=(mode.ratio-1)*ri
    y=(radius-ri)/width
    yz=-(1+(mode.ratio-1)*y)*qz/(2*q*(mode.ratio-1))
    scale=mode.nu*q**(-2*mode.inner.A)
    pr=np.zeros(len(points));pz=np.zeros(len(points))
    for coefficients,(lo,hi) in zip(mode.coefficients,mode.windows):
        bump,by=flat_bump(y,lo,hi)
        axial=np.polynomial.polynomial.polyval(eta/.3,coefficients)
        az=np.polynomial.polynomial.polyval(eta/.3,np.polynomial.polynomial.polyder(coefficients))*etaz/.3
        pr+=scale*by*axial/width
        pz+=scale*((by*yz-2*mode.inner.A*qz/q*bump)*axial+bump*az)
    return np.column_stack((pr*points[:,0]/safe,pr*points[:,1]/safe,pz))


def instantaneous_column(mode, points, tau, *, pressure_components=(0,1,2), hspace=None):
    """Return dR/dcoefficient at k0; omitted pressure components are zeroed.

    Use pressure_components=(2,) only when assembling tangential moment/cone
    integrals on the x-z half-plane. Default returns the complete vector.
    """
    points=np.asarray(points,float)
    if any(axis not in (0,1,2) for axis in pressure_components):
        raise ValueError('Cartesian pressure component must be 0, 1 or 2')
    if isinstance(mode,CompactPressureDirection):
        hs=5e-4*np.sqrt(mode.nu*tau) if hspace is None else float(hspace)
        out=np.zeros_like(points)
        for axis in pressure_components:
            e=np.zeros(3);e[axis]=hs
            pm2,pm,pp,pp2=[mode.primitive.increment(points+j*e,tau) for j in (-2,-1,1,2)]
            out[:,axis]=(pm2-8*pm+8*pp-pp2)/(12*hs)
        return out
    if type(mode) is OuterPressure:
        if not isinstance(mode.base,ZeroBackground):
            raise ValueError('Pressure direction must have a zero velocity background')
        out=pressure_gradient(mode,points,tau)
        out[:,[j for j in range(3) if j not in pressure_components]]=0.
        return out
    if not hasattr(mode,'k0') or abs(-np.log2(2*tau)-mode.k0)>1e-12:
        raise ValueError('Slope columns are valid only at the registered k0')
    if type(mode) is MeridionalSlope:
        if not isinstance(mode.unit.base,ZeroBackground):
            raise ValueError('Meridional unit must have zero background')
        unit=mode.unit.fields(points,tau)[0]
    elif type(mode) is BroadShearSlope:
        if not isinstance(mode.base.base,ZeroBackground):
            raise ValueError('Broad unit must have zero background')
        unit=mode.base.fields(points,tau)[0]
    elif type(mode) is OuterSwirlSlope:
        if not isinstance(mode.base,ZeroBackground):
            raise ValueError('Swirl unit must have zero background')
        value=SwirlValue(mode.base,mode.a)
        if tuple(mode.windows)!=tuple(value.windows):
            raise ValueError('Only the registered outer swirl windows are supported')
        unit=value.fields(points,tau)[0]
    else:
        raise TypeError(f'Unsupported initial-time mode: {type(mode).__name__}')
    return unit/(tau*np.log(2.))


def run():
    field,report=load_saved_field()
    tau=report['tau'];k=report['k']
    baseline=build_compact_baseline(field,report)
    modes,names=build_modes_with_parameters(field,report,k,compact_pressure=baseline)
    points=field.inner.from_similarity(field.join_X*(1+15*np.array([.067,.24,.53,.825]))**2,
        np.array([-.17,.11,.23,-.29]),tau,np.array([.31,.72,1.2,.47]))
    hs=5e-4*np.sqrt(field.nu*tau)
    started=time.perf_counter()
    fast=np.stack([instantaneous_column(mode,points,tau,hspace=hs) for mode in modes])
    fast_seconds=time.perf_counter()-started
    started=time.perf_counter()
    fd=np.stack([momentum(jets(mode,points,tau,hs,1e-4*tau)) for mode in modes])
    fd_seconds=time.perf_counter()-started
    error=np.max(abs(fast-fd),axis=(1,2))
    scale=np.maximum(np.max(abs(fd),axis=(1,2)),1.)
    result=dict(scope='Focused comparison of all 44 instantaneous columns at four off-axis points; not PDE or continuum acceptance.',
        k=k,tau=tau,points=points.tolist(),fast_seconds=fast_seconds,fd_seconds=fd_seconds,
        maximum_absolute_difference=float(error.max()),maximum_scaled_difference=float((error/scale).max()),
        rows=[dict(mode=name,absolute_difference=float(e),scaled_difference=float(e/s))
              for name,e,s in zip(names,error,scale)])
    path=Path(__file__).with_suffix('.json')
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:value for key,value in result.items() if key!='rows'}),flush=True)
    return result


if __name__=='__main__':run()
