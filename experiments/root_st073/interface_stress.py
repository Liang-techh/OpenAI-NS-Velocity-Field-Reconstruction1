"""Actual ST073 inner stress and relative momentum flux for transition matching."""
from pathlib import Path
import sys,json
import numpy as np
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'NS_ST073_Full_Local_Recurrence'))
from full_radial import FullRadialField

def gradient(f,X,e,tau):
    q=tau/(1-e*e);S=2*q;r=np.sqrt(S*X);co=f.coefficients(float(e),float(q))
    value=[np.polynomial.polynomial.polyval(X,c[:,0]) for c in co]
    dz=[np.polynomial.polynomial.polyval(X,c[:,1]) for c in co]
    ds=[np.polynomial.polynomial.polyval(X,np.polynomial.polynomial.polyder(c[:,0]))/S for c in co]
    a,b,c,p=value;az,bz,cz,pz=dz;ar,br,cr,pr=ds
    return np.array([[a+2*S*X*ar,-b,r*az],[b+2*S*X*br,a,r*bz],[2*r*cr,0,cz]],float)

def run():
    f=FullRadialField.load(ROOT/'NS_ST073_Full_Local_Recurrence/data/ST073-V.json');out=ROOT/'moving_interface';out.mkdir(exist_ok=True);rows=[]
    for order in (12,24):
        nodes,weights=leggauss(order)
        for k in (0,3,6):
            tau=.5*2**(-k);X=[];eta=[];area=[];labels=[]
            for e,w in zip(.5*nodes,.5*weights):
                q=tau/(1-e*e);r=np.sqrt(2*f.nu*q/64);rp=r*e/(1-e*e);zp=np.sqrt(f.nu)*q**f.D*(1+2*f.D*e*e/(1-e*e))
                X.append(1/64);eta.append(e);area.append(2*np.pi*r*w*np.array([zp,0,-rp]));labels.append(0)
            for sign in (-1,1):
                e=sign*.5;q=tau/(1-e*e)
                for x,w in zip((nodes+1)/128,weights/128):
                    X.append(x);eta.append(e);area.append([0,0,sign*2*np.pi*f.nu*q*w]);labels.append(sign)
            X=np.array(X);eta=np.array(eta);area=np.array(area);d=f.evaluate_similarity(X,eta,tau)
            xyz=f.from_similarity(X,eta,tau);boundary=xyz*np.array([-.5,0,-f.D])/tau
            jac=np.array([gradient(f,x,e,tau) for x,e in zip(X,eta)])
            stress=-d['pressure'][:,None,None]*np.eye(3)+f.nu*(jac+jac.transpose(0,2,1))
            traction=np.einsum('nij,nj->ni',stress,area)
            advective=d['velocity']*np.sum((d['velocity']-boundary)*area,axis=1)[:,None]
            flux=advective-traction
            omega=np.column_stack((jac[:,2,1]-jac[:,1,2],jac[:,0,2]-jac[:,2,0],jac[:,1,0]-jac[:,0,1]))
            rows.append(dict(order=order,k=k,axial_stress_integral=float(np.sum(traction[:,2])),axial_relative_advective_flux=float(np.sum(advective[:,2])),axial_net_outward_momentum_flux=float(np.sum(flux[:,2])),axial_torque_stress=float(np.sum(xyz[:,0]*traction[:,1])),axial_angular_momentum_outward_flux=float(np.sum(xyz[:,0]*flux[:,1])),gradient_trace_max=float(np.max(np.abs(np.trace(jac,axis1=1,axis2=2)))),vorticity_replay_error=float(np.max(np.abs(omega-d['vorticity'])))))
            if order==24:np.savez_compressed(out/f'stress_k{k}.npz',X=X,eta=eta,xyz=xyz,area_vector=area,region=labels,velocity=d['velocity'],pressure=d['pressure'],gradient=jac,stress=stress,boundary_velocity=boundary,traction_area=traction,relative_momentum_flux_area=flux)
    report=dict(rows=rows,definition='sigma=-p I+nu(grad u+grad u transpose); outward moving momentum flux=u((u-b) dot n)-sigma n',scope='Axisymmetric ring-integrated axial momentum and torque; Cartesian transverse totals cancel by symmetry. Inner boundary data only, no matched exterior or global acceptance.',global_field_ready=False)
    (out/'stress_flux.json').write_bytes((json.dumps(report,indent=2)+'\n').encode());print(json.dumps(rows[-3:],indent=2))
if __name__=='__main__':run()
