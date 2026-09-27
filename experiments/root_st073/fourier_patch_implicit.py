"""Weak viscous Galerkin evolution with implicit time integration.

Uses the same exact -curl coefficient convention as the spatial basis.
Diffusion is assembled as -nu M^+ K with positive-semidefinite gradient
Gram matrix K; full strong momentum is independently replayed afterward.
"""
import json
from types import SimpleNamespace
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline, CubicHermiteSpline
from numpy.polynomial.legendre import leggauss
from fourier_patch_evolution import Patch, PotentialField, pack, metrics
from outer_feedback_evolution import build_current, Trajectory
from affine_momentum import jets, momentum
from radial_continuation import ROOT


class GalerkinPatch(Patch):
    def __init__(self,*args):
        super().__init__(*args,build_solvers=False)
        self.operators=[]
        self.heat_spectra=[]
        wv=np.repeat(self.weights,3)
        wg=np.repeat(self.weights,9)
        for m,(V,G,L,P,Gp) in enumerate(self.cache):
            v=V.reshape(-1,3*self.q)
            g=G.reshape(-1,3*self.q)
            p=Gp.reshape(-1,self.q)
            U,s,Vh=np.linalg.svd(v*np.sqrt(wv)[:,None],full_matrices=False)
            keep=s>max(s)*1e-10
            mass_inverse=(Vh[keep].conj().T/(s[keep]**2))@Vh[keep]
            force=mass_inverse@(v.conj().T*wv)
            stiffness=g.conj().T@(wg[:,None]*g)
            heat=-self.mean.nu*mass_inverse@stiffness
            pu,ps,pvh=np.linalg.svd(p*np.sqrt(wv)[:,None],full_matrices=False)
            pinverse=(pvh.conj().T*(ps/(ps*ps+(max(ps)*1e-8)**2)))@pu.conj().T
            pinverse*=np.sqrt(wv)[None,:]
            self.operators.append((force,heat,pinverse))
            eig=np.linalg.eigvals(heat)
            self.heat_spectra.append(dict(mode=m,rank=int(sum(keep)),
                max_real=float(max(eig.real)),min_real=float(min(eig.real))))

    def control(self,state,background,with_pressure=True):
        residual=self.residual(state,background)
        slopes=[];pressures=[]
        for m,((V,G,L,P,Gp),(force,heat,pinverse)) in enumerate(zip(self.cache,self.operators)):
            source=np.mean(residual*self.phase[None,:,m,None].conj(),axis=1)
            source=source.real if m==0 else 2*source
            viscous=np.einsum('niq,q->ni',L,state[m])
            slope=-force@(source-viscous).reshape(-1)+heat@state[m]
            if m==0:slope=slope.real
            slopes.append(slope)
            if with_pressure:
                remainder=source+np.einsum('niq,q->ni',V,slope)
                pressure=-pinverse@remainder.reshape(-1)
                pressures.append(pressure.real if m==0 else pressure)
        return np.array(slopes),np.array(pressures),residual


def load_saved_field(path=ROOT/'fourier_patch_implicit.json'):
    """Load the experimental field; acceptance is recorded separately."""
    report=json.loads(path.read_text(encoding='utf-8'))
    data=report['state_data']
    decode=lambda a:np.asarray(a)[...,0]+1j*np.asarray(a)[...,1]
    _,_,current=build_current()
    saved=json.loads((ROOT/'outer_feedback_evolution.json').read_text(encoding='utf-8'))
    mean=Trajectory(current,saved['nodes'])
    times=np.asarray(data['physical_times'])
    state=CubicHermiteSpline(times,decode(data['states']),decode(data['slopes']),extrapolate=False)
    pressure=CubicSpline(times,decode(data['pressures']),extrapolate=False)
    patch=SimpleNamespace(mean=mean,center=tuple(report['center']),widths=tuple(report['widths']),
        degree=report['degree'],modes=report['modes'],carriers={int(m):c for m,c in report['carriers'].items()})
    field=PotentialField(patch,state,pressure);field.interval=(float(times[0]),float(times[-1]))
    return field,report


def run():
    previous=json.loads((ROOT/'fourier_patch_evolution.json').read_text(encoding='utf-8'))
    _,_,current=build_current()
    saved=json.loads((ROOT/'outer_feedback_evolution.json').read_text(encoding='utf-8'))
    mean=Trajectory(current,saved['nodes'])
    center,widths=tuple(previous['center']),tuple(previous['widths'])
    carriers={int(m):np.array(c) for m,c in previous['carriers'].items()}
    data=previous['state_data'];degree=data['degree']
    decode=lambda a:np.asarray(a)[...,0]+1j*np.asarray(a)[...,1]
    initial=decode(data['initial_state'])
    t0,_,t1=data['physical_times'];duration=t1-t0
    axis,w=leggauss(16)
    points=np.array([[center[0]+x*widths[0],0.,center[1]+z*widths[1]] for x in axis for z in axis])
    weights=np.outer(w,w).ravel()*points[:,0];weights/=sum(weights)
    patch=GalerkinPatch(mean,center,widths,carriers,degree,points,weights,-t0)
    report=dict(accepted=False,scale_recursion_established=False,method='Weak gradient-Gram diffusion and BDF in normalized physical time',
        degree=degree,center=list(center),widths=list(widths),carriers=previous['carriers'],modes=patch.modes,
        quadrature_order=16,heat_spectra=patch.heat_spectra,replay=[],
        scope='Spatially supported finite Fourier/polynomial Galerkin trajectory. Base jets interpolated at three times during integration; direct replay uses the original mean field. No volume-L2, energy, full-domain, endpoint or recursion acceptance.')
    path=ROOT/'fourier_patch_implicit.json'
    def save():path.write_bytes((json.dumps(report,indent=2)+'\n').encode())
    save();print(json.dumps(dict(stage='weak_heat_ready',heat_spectra=patch.heat_spectra)),flush=True)
    base=[]
    for fraction in (0.,.5,1.):
        base.append(patch.background(t0+duration*fraction))
        print(json.dumps(dict(stage='background_cached',fraction=fraction)),flush=True)
    interp=[CubicSpline([0.,.5,1.],np.array([row[j] for row in base])) for j in range(3)]
    background=lambda s:tuple(spline(s) for spline in interp)
    count=initial.size
    flatten=lambda a:np.r_[a.real.ravel(),a.imag.ravel()]
    unflatten=lambda a:(a[:count]+1j*a[count:]).reshape(initial.shape)
    calls=0
    def rhs(s,x):
        nonlocal calls
        calls+=1
        value=patch.control(unflatten(x),background(s),with_pressure=False)[0]
        if calls in (1,100,500,1000,2000,5000):
            print(json.dumps(dict(stage='implicit_rhs',calls=calls,fraction=float(s))),flush=True)
        return duration*flatten(value)
    solved=solve_ivp(rhs,(0.,1.),flatten(initial),method='BDF',rtol=1e-6,atol=1e-10,max_step=.1,dense_output=True)
    report['solver']=dict(success=bool(solved.success),message=solved.message,rhs_calls=calls,
                          steps=len(solved.t),nfev=solved.nfev,njev=solved.njev,nlu=solved.nlu)
    save();print(json.dumps(report['solver']),flush=True)
    if not solved.success:return
    fractions=np.unique(np.r_[solved.t,(solved.t[:-1]+solved.t[1:])/2])
    states=[];slopes=[];pressures=[]
    for s in fractions:
        state=unflatten(solved.sol(s))
        slope,pressure,_=patch.control(state,background(s))
        states.append(state);slopes.append(slope);pressures.append(pressure)
    times=t0+duration*fractions
    state_path=CubicHermiteSpline(times,np.array(states),np.array(slopes),extrapolate=False)
    pressure_path=CubicSpline(times,np.array(pressures),extrapolate=False)
    report['state_data']=dict(physical_times=times.tolist(),states=pack(states),slopes=pack(slopes),pressures=pack(pressures))
    field=PotentialField(patch,state_path,pressure_path);field.interval=(t0,t1)
    save()
    held_axis=np.linspace(-.8,.8,6)
    offsets=[(x,z) for x in held_axis for z in held_axis]+[(x,z) for x in (-.97,.97) for z in (-.97,.97)]
    angles=(np.arange(40)+.37)*2*np.pi/40
    held=np.array([[(center[0]+x*widths[0])*np.cos(a),(center[0]+x*widths[0])*np.sin(a),center[1]+z*widths[1]]
                   for x,z in offsets for a in angles])
    for fraction in (.25,.5,.75):
        t=t0+duration*fraction;tau=-t
        jet=jets(field,held,tau,.0005*np.sqrt(mean.nu*tau),.000005*tau)
        row=dict(fraction=fraction,k=float(-np.log2(2*tau)),direct_momentum=metrics(momentum(jet)),
                 sampled_divergence=float(max(abs(np.trace(jet[1],axis1=1,axis2=2)))))
        report['replay'].append(row);save();print(json.dumps(row),flush=True)
        if row['direct_momentum']['max']>=previous['projections'][0]['held_before']['max']:
            report['stopped_reason']='Implicit trajectory fails to reduce the original spatial holdout maximum.'
            save();return


if __name__=='__main__':run()
