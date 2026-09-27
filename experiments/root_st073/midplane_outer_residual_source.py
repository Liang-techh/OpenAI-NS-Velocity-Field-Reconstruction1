"""Locate full-residual sources and screen outer correction-stress geometry."""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from affine_momentum import jets,momentum
from saved_outer_dae_field import SavedOuterDAEField
from midplane_resolved_feasibility import RADIAL_BREAKS
from joined_field import coordinates
from radial_continuation import ROOT


def terms(field,point,tau):
    p=np.asarray(point)[None,:];h=.0005*np.sqrt(field.nu*tau);ht=.0001*tau
    u,_=field.fields(p,tau);grad=np.zeros((3,3));gp=np.zeros(3);lap=np.zeros(3)
    for i in range(3):
        e=np.eye(3)[i]*h
        samples=[field.fields(p+j*e,tau) for j in (-2,-1,0,1,2)]
        grad[:,i]=(samples[0][0]-8*samples[1][0]+8*samples[3][0]-samples[4][0])[0]/(12*h)
        gp[i]=(samples[0][1]-8*samples[1][1]+8*samples[3][1]-samples[4][1])[0]/(12*h)
        lap+=(-samples[4][0]+16*samples[3][0]-30*samples[2][0]+16*samples[1][0]-samples[0][0])[0]/(12*h*h)
    tt=[field.fields(p,tau+j*ht)[0][0] for j in (-2,-1,1,2)]
    out=dict(time=-(tt[0]-8*tt[1]+8*tt[2]-tt[3])/(12*ht),advection=grad@u[0],pressure=gp,viscous=-field.nu*lap)
    out['sum']=sum(out.values())
    return {key:value.tolist() for key,value in out.items()}


def outer_cones(field,k,order=64):
    inner=field.inner;tau=.5*2.**-k;centers=[];labels=[]
    for eta in (-.2,0.,.2):
        for y in (.5,.75):
            centers.append(inner.from_similarity([inner.p.X_max*(1+15*y)**2],[eta],tau)[0]);labels.append((eta,y))
    points=list(centers);panels=[];g,w=leggauss(order)
    for R,_,z in centers:
        q=float(coordinates(0.,z/np.sqrt(inner.nu),tau,inner.h)['q']);ri=np.sqrt(2*inner.nu*q*inner.p.X_max)
        edges=sorted(set(np.clip([0.,ri,R,*[ri*(1+15*b) for b in RADIAL_BREAKS]],0.,R)))
        rr=np.concatenate([(lo+hi)/2+(hi-lo)*g/2 for lo,hi in zip(edges[:-1],edges[1:])]);ww=np.concatenate([(hi-lo)*w/2 for lo,hi in zip(edges[:-1],edges[1:])])
        start=len(points);points.extend(np.column_stack((rr,np.zeros_like(rr),np.full_like(rr,z))));panels.append((slice(start,len(points)),rr,ww,R))
    jet=jets(field,np.array(points),tau,.0005*np.sqrt(inner.nu*tau),.0001*tau);res=momentum(jet);rows=[]
    for i,(sl,rr,ww,R) in enumerate(panels):
        u,J=jet[0][i],jet[1][i];F=u[1]/R;shear=np.array([J[1,0]-F,J[2,0]]);mag=np.linalg.norm(shear);N=shear/max(mag,1e-30);K=np.array([-N[1],N[0]])
        lam=-2*F*N[0]*(2*F*N[0]+mag);target=np.array([-np.dot(ww*rr**2,res[sl,1])/R**2,-np.dot(ww*rr,res[sl,2])/R]);dn,dk=target@N,target@K
        ratio=float(np.sqrt(max(lam,0))*abs(dk)/max(abs(2*F*N[0]*dn),1e-30))
        rows.append(dict(eta=labels[i][0],radial_fraction=labels[i][1],lambda_squared=float(lam),target_dot_N=float(dn),cone_ratio=ratio,cone_pass=bool(lam>0 and dn<0 and ratio<1),momentum_norm=float(np.linalg.norm(res[i]))))
    return rows


def run():
    f=SavedOuterDAEField();k=13.05;tau=.5*2.**-k
    ys=np.linspace(.01,.99,81);etas=(-.3,-.2,-.1,0.,.1,.2,.3)
    points=np.concatenate([f.inner.from_similarity(f.inner.p.X_max*(1+15*ys)**2,np.full(len(ys),eta),tau) for eta in etas])
    residual=momentum(jets(f,points,tau,.0005*np.sqrt(f.nu*tau),.0001*tau));norm=np.linalg.norm(residual,axis=1);i=int(np.argmax(norm))
    decomposition=terms(f,points[i],tau);np.testing.assert_allclose(decomposition['sum'],residual[i],rtol=1e-8,atol=.1)
    peak=dict(point=points[i].tolist(),eta=etas[i//len(ys)],radial_fraction=float(ys[i%len(ys)]),norm=float(norm[i]),terms=decomposition)
    print(json.dumps(dict(stage='peak',**peak)),flush=True)
    rows=outer_cones(f,k)
    report=dict(k=k,sampled_peak=peak,component_peaks=np.max(abs(residual),axis=0).tolist(),outer_cones=rows,
        scope='Mid-interval diagnostic on an axisymmetric grid and six outer-window stress nodes; not a domain supremum or whole-support cone certificate.',accepted=False,scale_recursion_established=False)
    (ROOT/'midplane_outer_residual_source.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps(dict(outer_cone_pass=sum(r['cone_pass'] for r in rows),rows=rows)),flush=True)

if __name__=='__main__':run()
