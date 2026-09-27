"""Radial centrifugal pressure for a broad swirl increment.

Inner datum preserves core pressure. Its nonzero exterior pressure trace
must be matched dynamically; this wrapper is not a whole-space solution.
"""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from broad_annular_shear import BroadAnnularShear
from fourier_shear_feasibility import load_saved_field, ControlledMean
from midplane_resolved_feasibility import ZeroBackground
from joined_field import coordinates
from affine_momentum import jets, momentum
from radial_continuation import ROOT
from separated_moment_modes import flat_transition


class BroadCentrifugalPressure:
    """Assumes modified velocity = reference velocity + specified broad swirl."""
    def __init__(self, modified, reference, amplitude, parameters, radial_breaks,
                 datum='inner', order=16, compact_window=(.72,.93)):
        self.modified,self.reference=modified,reference
        self.amplitude=amplitude
        self.parameters={k:v for k,v in parameters.items() if k!='amplitude'}
        self.radial_breaks=tuple(radial_breaks)
        self.datum,self.order=datum,int(order)
        self.compact_window=tuple(compact_window)
        if datum not in ('inner','outer','compact'):
            raise ValueError('Pressure datum must be inner, outer or compact.')
        if not 0<self.compact_window[0]<self.compact_window[1]<1:
            raise ValueError('Compact primitive cutoff must lie inside the annulus.')
        for name in ('inner','nu','join_X','ratio'):
            setattr(self,name,getattr(modified,name))

    def increment(self,points,tau):
        points=np.asarray(points,float)
        result=np.zeros(len(points))
        radial=np.hypot(points[:,0],points[:,1])
        g,w=leggauss(self.order)
        broad=BroadAnnularShear(ZeroBackground(self.reference),self.amplitude(tau),**self.parameters)
        for z in np.unique(points[:,2]):
            indices=np.flatnonzero(points[:,2]==z)
            q=float(coordinates(0.,z/np.sqrt(self.nu),tau,self.inner.h)['q'])
            ri=np.sqrt(2*self.nu*q*self.join_X);ro=self.ratio*ri
            cuts=np.unique(np.clip(np.r_[ri,ro,radial[indices],
                             ri*(1+(self.ratio-1)*np.asarray(self.radial_breaks))],ri,ro))
            half=np.diff(cuts)/2
            nodes=(cuts[:-1,None]+cuts[1:,None])/2+half[:,None]*g
            sample=np.column_stack((nodes.ravel(),np.zeros(nodes.size),np.full(nodes.size,z)))
            old=self.reference.fields(sample,tau)[0][:,1].reshape(nodes.shape)
            added=broad.fields(sample,tau)[0][:,1].reshape(nodes.shape)
            integral=np.sum(half[:,None]*w*(2*old*added+added*added)/nodes,axis=1)
            primitive=np.r_[0.,np.cumsum(integral)]
            target=np.clip(radial[indices],ri,ro)
            result[indices]=primitive[np.searchsorted(cuts,target)]
            if self.datum=='outer':result[indices]-=primitive[-1]
            if self.datum=='compact':
                y=(radial[indices]/ri-1)/(self.ratio-1)
                lo,hi=self.compact_window
                chi=np.array([flat_transition((value-lo)/(hi-lo)) for value in y])
                result[indices]-=chi*primitive[-1]
        return result

    def fields(self,points,tau):
        u,p=self.modified.fields(points,tau)
        return u,p+self.increment(points,tau)


def run():
    modified,repair=load_saved_field()
    if not isinstance(modified.current,BroadAnnularShear):
        raise ValueError('Expected saved broad-shear mean.')
    reference=ControlledMean(modified.current.base,repair['state'],repair['control'],repair['k'])
    params=repair['mean_increment']['parameters']
    repaired=BroadCentrifugalPressure(modified,reference,lambda tau:params['amplitude'],
                                     params,repair['radial_breaks'])
    tau=.5*2.**-repair['k']
    seed=json.loads((ROOT/'broad_shear_growth.json').read_text())
    r,z=seed['center']
    # Select the sampled centrifugal peak without running a second full-grid residual.
    ys=np.linspace(.01,.99,61)
    grid=np.concatenate([modified.inner.from_similarity(modified.inner.p.X_max*(1+15*ys)**2,
                             np.full(len(ys),eta),tau) for eta in (-.3,-.1,0.,.1,.3)])
    old=reference.fields(grid,tau)[0];new=modified.fields(grid,tau)[0]
    change=(new[:,1]**2-old[:,1]**2)/grid[:,0]
    peak=grid[int(np.argmax(abs(change)))]
    outer=modified.inner.from_similarity([modified.inner.p.X_max*16**2],[.2],tau)[0]
    points=np.array([[r,0,z],peak,outer])
    hs=.0005*np.sqrt(modified.nu*tau);ht=.0001*tau
    before=momentum(jets(modified,points,tau,hs,ht))
    after=momentum(jets(repaired,points,tau,hs,ht))
    baseline=momentum(jets(reference,points,tau,hs,ht))
    values=repaired.increment(points,tau)
    compact=BroadCentrifugalPressure(modified,reference,lambda tau:params['amplitude'],
                                     params,repair['radial_breaks'],datum='compact')
    collar=modified.inner.from_similarity([modified.inner.p.X_max*(1+15*.825)**2],[0.],tau)[0]
    compact_points=np.vstack((points,collar))
    compact_residual=momentum(jets(compact,compact_points,tau,hs,ht))
    report=dict(accepted=False,source='fourier_shear_feasibility.json',k=repair['k'],
                datum='inner',quadrature_order=16,
                scope='Three diagnostic points only. Radial pressure repair changes axial balance and the outer pressure trace; no volume-L2, global forcing, energy or recursion acceptance.',
                labels=['wave_center','sampled_centrifugal_peak','outer_eta_0.2'],
                points=points.tolist(),pressure_increments=values.tolist(),
                momentum_before=before.tolist(),momentum_after=after.tolist(),
                momentum_reference_without_broad=baseline.tolist(),
                norms_before=np.linalg.norm(before,axis=1).tolist(),
                norms_after=np.linalg.norm(after,axis=1).tolist(),
                radial_difference_from_reference=(after[:,0]-baseline[:,0]).tolist(),
                outer_pressure_continuation_requires_matching=True,
                compact_radial_increment_budget=dict(
                    eta=.2,full_pressure_integral=float(values[2]),
                    annulus_width=float(outer[0]*(1-1/modified.ratio)),
                    supremum_lower_bound=float(abs(values[2])/(outer[0]*(1-1/modified.ratio))),
                    scope='For this pure swirl velocity increment and any pressure increment zero at both radial boundaries, integral(delta R_r)=-integral(g). This bounds the added radial residual, not the full residual of a future meridional/wave-corrected field.'),
                compact_trial=dict(cutoff_window=list(compact.compact_window),
                    labels=['wave_center','sampled_centrifugal_peak','outer_eta_0.2','cutoff_midpoint_eta_0'],
                    points=compact_points.tolist(),
                    pressure_increments=compact.increment(compact_points,tau).tolist(),
                    momentum=compact_residual.tolist(),
                    norms=np.linalg.norm(compact_residual,axis=1).tolist(),
                    remainder='Radial cutoff derivative times the full centrifugal pressure integral is retained in the reported momentum.'))
    (ROOT/'broad_shear_pressure.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report),flush=True)
    return report


if __name__=='__main__':run()
