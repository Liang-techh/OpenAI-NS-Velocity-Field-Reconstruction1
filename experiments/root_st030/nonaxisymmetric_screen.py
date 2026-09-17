"""ST031: exact-energy, full-curl, nonaxisymmetric annulus capacity screen.
Fixed existing force/pressure support. The curvature and trial scores are training
quantities and cannot promote any full NS acceptance flag.
"""
from __future__ import annotations
import argparse, hashlib, json,time
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import cholesky,solve_triangular,qr,eigh
from spacetime import Family,force,NU
from annular_curl_modes import CurlModes

def cylindrical_grid(order,angles=12):
    g,w=leggauss(order);s,z,th=np.meshgrid(2*(g+1),2*g,2*np.pi*np.arange(angles)/angles,indexing='ij')
    wg=np.broadcast_to((4*np.pi/angles*np.outer(w,w))[:,:,None],s.shape)
    return np.column_stack((np.sqrt(s).ravel()*np.cos(th).ravel(),np.sqrt(s).ravel()*np.sin(th).ravel(),z.ravel())),wg.ravel()

def time_basis(t):
    q=4*np.asarray(t)-2
    return np.stack((np.ones_like(q),q,2*q*q-1),axis=-1),np.stack((np.zeros_like(q),np.full_like(q,4),16*q),axis=-1)

def base_jets(f,raw,x,t):
    xx,yy,z=x.T;s=xx*xx+yy*yy;D=f.bundle(s,z,t);aa,bb,qq,fc,_=f.coefficients(raw)
    V={k:D[k]@(bb if k.startswith('B') else qq if k.startswith('Q') else aa) for k in D if k not in ('s','z','t')}
    A,B,C=V['A'],V['B'],V['C'];As,Az,Bs,Bz=V['As'],V['Az'],V['Bs'],V['Bz'];Cs,Cz=V['Cs'],V['Cz']
    u=np.column_stack((xx*A-yy*B,yy*A+xx*B,C))
    ut=np.column_stack((xx*V['At']-yy*V['Bt'],yy*V['At']+xx*V['Bt'],V['Ct']))
    lap=np.column_stack((xx*V['AL']-yy*V['BL'],yy*V['AL']+xx*V['BL'],V['CL']))
    J=np.empty((len(x),3,3));J[:,0,0]=A+2*xx*xx*As-2*xx*yy*Bs;J[:,0,1]=2*xx*yy*As-B-2*yy*yy*Bs;J[:,0,2]=xx*Az-yy*Bz
    J[:,1,0]=2*xx*yy*As+B+2*xx*xx*Bs;J[:,1,1]=A+2*yy*yy*As+2*xx*yy*Bs;J[:,1,2]=yy*Az+xx*Bz
    J[:,2,0]=2*xx*Cs;J[:,2,1]=2*yy*Cs;J[:,2,2]=Cz
    pg=np.column_stack((2*xx*V['Qs'],2*yy*V['Qs'],V['Qz']))
    return u,J,ut-NU*lap,pg,force(x,t,*fc)

class AnnulusScreen:
    def __init__(self,candidate,order=20,angles=12,time_order=5,energy_order=64):
        self.f,self.raw=Family.load(candidate);self.m=CurlModes(ms=(1,2),radial_degrees=(0,1),axial_degrees=(0,1));self.ns=self.m.size;self.nt=3;self.n=3*self.ns
        # Fine energy normalization independent of the coarser training quadrature.
        xe,we=cylindrical_grid(energy_order,8);ue=self.m.jets(xe)[0]
        M=np.einsum('nki,nli,n->kl',ue,ue,we,optimize=True)
        self.transform=solve_triangular(cholesky(M,lower=True).T,np.eye(self.ns),lower=False)
        self.mass=self.transform.T@M@self.transform
        self.mass0=np.kron(self.mass,np.outer([1,-1,1],[1,-1,1]))
        self.x,self.w=cylindrical_grid(order,angles);self.w=self.w/self.w.sum();U,J,L,P=self.m.jets(self.x)
        self.U=np.einsum('nki,kl->nli',U,self.transform);self.Jw=np.einsum('nkij,kl->nlij',J,self.transform);self.L=np.einsum('nki,kl->nli',L,self.transform)
        tg,tw=leggauss(time_order);self.times=np.r_[.25,.5+.25*tg,.75];self.tw=np.r_[.1,.4*tw,.1]
        self.T,self.Td=time_basis(self.times);self.base=[base_jets(self.f,self.raw,self.x,t) for t in self.times]
        N=len(self.x);n=self.n
        self.Lin=np.empty((len(self.times)*N*3,n));Pcols=np.empty((len(self.times)*N*3,P.shape[1]*self.nt));r=[];H2=np.zeros((n,n));S=0.
        for k,(ub,Jb,Lb,pg,fb) in enumerate(self.base):
            wt=self.w*self.tw[k];R=Lb+np.einsum('nij,nj->ni',Jb,ub)+pg-fb;r.append((np.sqrt(wt)[:,None]*R).ravel())
            cross=np.einsum('nkij,nj->nki',self.Jw,ub)+np.einsum('nij,nkj->nki',Jb,self.U)
            A=(self.U[:,:,:,None]*self.Td[k]+(cross-NU*self.L)[:,:,:,None]*self.T[k]).transpose(0,2,1,3).reshape(N*3,n)
            self.Lin[k*N*3:(k+1)*N*3]=A*np.repeat(np.sqrt(wt),3)[:,None]
            Pcols[k*N*3:(k+1)*N*3]=(P[:,:,:,None]*self.T[k]).transpose(0,2,1,3).reshape(N*3,-1)*np.repeat(np.sqrt(wt),3)[:,None]
            K=np.einsum('na,nib,njab,n->ij',R,self.U,self.Jw,wt,optimize=True)
            H2+=np.kron(K+K.T,np.outer(self.T[k],self.T[k]))
            S+=np.einsum('ni,ni,n',R,Lb+2*np.einsum('nij,nj->ni',Jb,ub),wt)
        self.r0=np.concatenate(r);self.P=Pcols;psc=1/np.linalg.norm(Pcols,axis=0);self.Q,self.PR=qr(Pcols*psc,mode='economic',check_finite=False);self.psc=psc
        Lproj=self.Lin-self.Q@(self.Q.T@self.Lin);self.N=Lproj.T@Lproj
        self.H=self.N+H2-.5*S*self.mass0;self.H=(self.H+self.H.T)/2
        self.grad=Lproj.T@(self.r0-self.Q@(self.Q.T@self.r0))
    def residual(self,c,return_pressure=False):
        c=np.asarray(c);a=c.reshape(self.ns,self.nt);lam=(1+.5*c@self.mass0@c)**-.5;rows=[]
        for k,(ub,Jb,Lb,pg,fb) in enumerate(self.base):
            ct=a@self.T[k];ctd=a@self.Td[k]
            v=ub+np.einsum('nki,k->ni',self.U,ct);J=Jb+np.einsum('nkij,k->nij',self.Jw,ct)
            Lv=Lb+np.einsum('nki,k->ni',self.U,ctd)-NU*np.einsum('nki,k->ni',self.L,ct)
            R=lam*Lv+lam**2*np.einsum('nij,nj->ni',J,v)+pg-fb
            rows.append((np.sqrt(self.w*self.tw[k])[:,None]*R).ravel())
        R=np.concatenate(rows);q=-self.psc*solve_triangular(self.PR,self.Q.T@R,lower=False,check_finite=False);rr=R+self.P@q
        return (rr,q,lam) if return_pressure else rr

def run(candidate,out,order=20,energy_order=64):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    reg=dict(experiment='CR-ROOT-ST031',source_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),mode_type='curl of polynomial annular potentials, angular m=1,2 both phases; radial degree0,1 axial degree0,1 time degree2',velocity_parameters=96,radial_support=[.2,1.8],axial_support=[-1.8,1.8],training_quadrature_order=order,angular_nodes=12,time_gauss_order=5,energy_order=energy_order,amplitude_trials=[0,.001,.003,.01,.03,.1,.2],new_mode_coefficient_bound=.5,force='same frozen original a,c',pressure='frozen parent plus bounded compact nonaxisymmetric modes',physical_gates='nu/time/support/force/energy/thresholds retained; GLOBAL AXISYMMETRY REMOVED FOR DIAGNOSTIC ONLY; not eligible for original acceptance',created_before_run=True)
    (out/'registration.json').write_text(json.dumps(reg,indent=2)+'\n');t=time.time();o=AnnulusScreen(candidate,order=order,energy_order=energy_order)
    sc=1/np.maximum(np.sqrt(np.maximum(np.diag(o.N),0)),1e-8);H=o.H*sc[:,None]*sc[None,:];vals,vec=eigh(H);direction=sc*vec[:,0];direction/=np.max(np.abs(direction))
    baseline=float(o.residual(np.zeros(o.n))@o.residual(np.zeros(o.n)));trials=[]
    for eps in reg['amplitude_trials']:
        for sign in ([1] if eps==0 else [1,-1]):
            rr,pc,lam=o.residual(sign*eps*direction,True);trials.append(dict(amplitude=sign*eps,loss=float(rr@rr),ratio=float(rr@rr)/baseline,normalization=lam,pressure_coeff_max=float(np.max(abs(pc)))))
    # Independent scalar second difference checks the assembled exact Hessian.
    rng=np.random.default_rng(9172907);v=rng.normal(size=o.n);v/=np.linalg.norm(v);hh=1e-4
    e0=.5*baseline;ep=.5*np.sum(o.residual(hh*v)**2);em=.5*np.sum(o.residual(-hh*v)**2);second=(ep-2*e0+em)/hh**2;expected=v@o.H@v
    report=dict(**reg,elapsed=time.time()-t,smallest_scaled_curvatures=vals[:10].tolist(),gradient_max=float(np.max(abs(o.grad))),baseline_training_mse=baseline,pressure_projection_orthogonality=float(np.max(abs(o.Q.T@(o.r0-o.Q@(o.Q.T@o.r0))))),mass_orthogonality=float(np.max(abs(o.mass-np.eye(o.ns)))),hessian_second_difference_relative_error=float(abs(second-expected)/max(abs(expected),1e-20)),trials=trials,pde_validated=False,scope='No held-out PDE acceptance; signs of a finite-dimensional Hessian are only local capacity evidence.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');np.savez_compressed(out/'curvature.npz',H=o.H,N=o.N,gradient=o.grad,direction=direction,transform=o.transform,mass0=o.mass0)
    print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);p.add_argument('--order',type=int,default=20);p.add_argument('--energy-order',type=int,default=64);a=p.parse_args();run(a.candidate,a.out,a.order,a.energy_order)
