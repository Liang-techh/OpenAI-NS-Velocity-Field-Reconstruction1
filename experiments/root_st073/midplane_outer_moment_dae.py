"""Local value/slope moment restoration for the interscale field."""
import argparse,json
import numpy as np
from scipy.optimize import least_squares
from adaptive_bridge_recursive_defect import build_fields
from adaptive_bridge_moment_fit import moment_slices,outer_moments
from affine_momentum import combine,momentum
from midplane_resolved_feasibility import ZeroBackground,RADIAL_BREAKS,_integrated_moment_coefficients,_evaluate,_jacobian,_cone_replay
from separated_moment_modes import SeparatedMomentModes,RADIAL_WINDOWS_THREE
from radial_continuation import ROOT
FREE=np.array([3,4,5,6,7,8,12,13,14,15,16,17])

class LocalTimeCorrection:
    def __init__(self,base,values,slopes,k0):
        self.base=base;self.values=np.asarray(values);self.slopes=np.asarray(slopes);self.k0=k0
        for name in ('inner','nu','join_X','ratio'):setattr(self,name,getattr(base,name))
    def fields(self,points,tau):
        times=np.asarray(tau).ravel();k=-np.log2(2*float(times[0]))
        a=np.zeros(18);a[FREE]=self.values+(k-self.k0)*self.slopes
        correction=SeparatedMomentModes(ZeroBackground(self.base),np.tile(a,3),windows=RADIAL_WINDOWS_THREE,knots=(11.,15.,19.),axial_powers=(0,1,2))
        u,p=self.base.fields(points,tau);v,q=correction.fields(points,tau)
        return u+v,p+q


def run(k=13):
    inner,fields=build_fields();base=fields['two_sided_cone']
    a=np.array(json.loads((ROOT/'midplane_outer_axial_interscale.json').read_text())['amplitudes'])
    current=SeparatedMomentModes(base,a,windows=RADIAL_WINDOWS_THREE,knots=(11.,15.,19.),axial_powers=(0,1,2))
    units=[]
    for index in FREE:
        b=np.zeros(18);b[index]=1.
        units.append(SeparatedMomentModes(ZeroBackground(base),np.tile(b,3),windows=RADIAL_WINDOWS_THREE,knots=(11.,15.,19.),axial_powers=(0,1,2)))
    data=moment_slices(inner,base,current,orders=(float(k),),n=48,unit_fields=units,unit_fields_are_deltas=True,radial_breaks=RADIAL_BREAKS)[0]
    c=_integrated_moment_coefficients(data);tau=.5*2.**-k
    # Physical time is minus remaining time, hence -dk/dtau=1/(tau ln2).
    response=data['modes'][0]/(tau*np.log(2))
    rows=[]
    for panel in data['panels']:
        sl=panel['slice'];rr=panel['radii'];ww=panel['weights'];ro=panel['outer_radius']
        rows.extend([-np.einsum('n,pn->p',ww*rr**2/ro**2,response[:,sl,1]),-np.einsum('n,pn->p',ww*rr/ro,response[:,sl,2])])
    D=np.array(rows)
    norm=max(1.,max(abs(c[0])))
    def embed(p):return np.r_[np.zeros(6),p]
    fit=least_squares(lambda p:_evaluate(c,embed(p))[[1,3]]/norm,np.zeros(6),
        jac=lambda p:_jacobian(c,embed(p))[[1,3],6:]/norm,max_nfev=200,ftol=1e-12,xtol=1e-12,gtol=1e-12)
    values=embed(fit.x);m=_evaluate(c,values)
    slopes=np.zeros(12);slopes[:6]=np.linalg.lstsq(D[[0,2],:6],-m[[0,2]],rcond=None)[0]
    predicted=m+D@slopes
    corrected=LocalTimeCorrection(current,values,slopes,k)
    print(json.dumps(dict(k=k,initial=c[0].tolist(),predicted=predicted.tolist(),values_norm=float(np.linalg.norm(values)),slopes_norm=float(np.linalg.norm(slopes)))),flush=True)
    replay=[]
    for order in (96,128):
        d=moment_slices(inner,base,corrected,orders=(float(k),),n=order,unit_fields=[corrected],radial_breaks=RADIAL_BREAKS)[0]
        v=outer_moments(d,np.zeros(1));r=momentum(d['baseline'])
        replay.append(dict(order=order,moments=v.tolist(),absolute_max=float(max(abs(v))),momentum_peak=float(max(np.linalg.norm(r,axis=1)))))
        print(json.dumps(replay[-1]),flush=True)
    supports=json.loads((ROOT/'midplane_physical_covariance_pairs.json').read_text())
    cone=_cone_replay(inner,corrected,supports,k,order=96)
    report=dict(k=k,optimizer_success=bool(fit.success),initial_moments=c[0].tolist(),value_only_moments=m.tolist(),
        coefficient_values=values.tolist(),coefficient_slopes=slopes.tolist(),slope_response=D.tolist(),
        slope_singular_values=np.linalg.svd(D,compute_uv=False).tolist(),predicted_moments=predicted.tolist(),replay=replay,cone=cone,
        scope='One local affine-in-log-scale coefficient trajectory. Independent field derivatives included. Not an integrated DAE, interval solution, or recursive construction.',accepted=False,scale_recursion_established=False)
    (ROOT/f'midplane_outer_moment_dae_k{k}.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps(dict(k=k,cone_pass=cone['pass_count'])),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--k',type=int,choices=(13,17),default=13)
    run(parser.parse_args().k)
