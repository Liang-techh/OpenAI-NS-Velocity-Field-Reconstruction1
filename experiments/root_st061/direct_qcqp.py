"""Direct convex quadratic constraints in a Hessian-metric search subspace.
This solves only a linearized-residual subproblem, never certifies the PDE.
All returned steps require full nonlinear checks by the calling optimizer.
"""
from __future__ import annotations
import time
import numpy as np
from scipy.linalg import cholesky,solve_triangular,qr
from scipy.optimize import minimize

def build_subspace(H,g,A,v,Ap,p,cap,damping,rank=96):
    Hreg=(H+H.T)/2+damping*np.eye(len(g));L=cholesky(Hreg,lower=True)
    uw=-solve_triangular(L,g,lower=True)
    E=solve_triangular(L,A.T,lower=True).T
    F=solve_triangular(L,Ap.T,lower=True).T
    norms=np.maximum(np.linalg.norm(E,axis=1),1e-14)
    deficit=(-v-E@uw)/norms
    ids=np.argsort(deficit)[-max(2*rank,128):]
    # Include descent, all peak normals and most obstructing structural normals.
    B=np.column_stack((uw,E[ids].T/norms[ids],F.T/np.maximum(np.linalg.norm(F,axis=1),1e-14)))
    B/=np.maximum(np.linalg.norm(B,axis=0),1e-20)
    # Preserve descent as the first column, then pivot the orthogonal complement.
    first=uw/max(np.linalg.norm(uw),1e-20)
    remainder=B-first[:,None]*(first@B)[None,:]
    U,T,piv=qr(remainder,mode='economic',pivoting=True)
    valid=int(np.count_nonzero(abs(np.diag(T))>1e-9))
    U=np.column_stack((first,U[:,:min(rank-1,valid)]))
    D=solve_triangular(L.T,U,lower=False)
    return D,dict(dimension=D.shape[1],full_dimension=len(g),metric_orthogonality_error=float(abs(D.T@Hreg@D-np.eye(D.shape[1])).max()),initial_gradient_norm=float(np.linalg.norm(uw)))

def quadratic_constraints(R,J,weights,target):
    """Encode weighted ||R+Jd||^2 <= target as PSD quadratic forms."""
    nr,nt,_,k=J.shape;Gs=[];bs=[];fs=[];caps=[]
    for a,w in enumerate(weights):
        ids=np.flatnonzero(w>0)
        for t in range(nt):
            A=J[ids,t].reshape(-1,k);rr=R[ids,t].ravel();ww=np.repeat(w[ids],3)
            Gs.append(A.T@(ww[:,None]*A));bs.append(A.T@(ww*rr));fs.append(float(ww@(rr*rr)));caps.append(float(target[a,t]))
    return np.asarray(Gs),np.asarray(bs),np.asarray(fs),np.asarray(caps)

def solve_subproblem(H,g,A,v,R,J,weights,targets,Rp,Jp,cap,D,damping,maxiter=160):
    started=time.monotonic();k=D.shape[1]
    G,b,f,limits=quadratic_constraints(R,J,weights,targets)
    scales=np.maximum(limits,1e-14)
    LA=A@D;lnorm=np.maximum(np.linalg.norm(LA,axis=1),1e-12);sv=v/lnorm;SA=LA/lnorm[:,None]
    h=D.T@((H+H.T)/2+damping*np.eye(len(g)))@D;gg=D.T@g
    # Keep all linear safeguards, but add them lazily after each subproblem solve.
    uncon=-np.linalg.solve(h,gg)
    active=set(np.argsort(sv+SA@uncon)[:min(220,len(v))].tolist())
    y=np.zeros(k);records=[]
    def energy(y):
        Gy=G@y
        val=(limits-f-2*b@y-np.einsum('i,ni->n',y,Gy))/scales
        der=-2*(b+Gy)/scales[:,None]
        return val,der
    def peak(y):
        rr=Rp+Jp@y
        return 1-np.sum(rr*rr,axis=1)/cap**2,-2*np.einsum('ni,nik->nk',rr,Jp)/cap**2
    # Explicit bounded trust ball, in the Hessian metric of the full search.
    radius=max(5*np.linalg.norm(uncon),.05)
    ret=None
    for exchange in range(5):
        ids=np.array(sorted(active))
        def con(y):
            ev,ej=energy(y);pv,pj=peak(y)
            return np.r_[sv[ids]+SA[ids]@y,ev,pv,1-(y@y)/radius**2],np.vstack((SA[ids],ej,pj,-2*y[None,:]/radius**2))
        last=[None,None]
        def cache(y):
            if last[0] is None or not np.array_equal(last[0],y):last[:]=[y.copy(),con(y)]
            return last[1]
        yscale=max(np.linalg.norm(uncon),1e-8)
        ret=minimize(lambda z:(.5*z@h@z+(gg/yscale)@z,h@z+gg/yscale),y/yscale,jac=True,method='SLSQP',
            constraints=[dict(type='ineq',fun=lambda z:cache(yscale*z)[0],jac=lambda z:yscale*cache(yscale*z)[1])],
            options=dict(maxiter=maxiter,ftol=1e-10,disp=False))
        y=yscale*ret.x;bad=np.flatnonzero(sv+SA@y < -2e-8)
        records.append(dict(exchange=exchange,success=bool(ret.success),message=str(ret.message),iterations=int(ret.nit),linear_min=float(np.min(sv+SA@y)),energy_min=float(energy(y)[0].min()),peak_min=float(peak(y)[0].min()),linear_active=len(ids)))
        if not len(bad):break
        active.update(bad[np.argsort((sv+SA@y)[bad])[:128]].tolist())
    return D@y,dict(records=records,elapsed=time.monotonic()-started,trust_radius=radius,reduced_dimension=k,linear_constraints=len(v),quadratic_constraints=len(f),pointwise_norm_constraints=len(Rp),objective=float(.5*y@h@y+gg@y),scipy_success=bool(ret.success),scope='Numerical convex QCQP in selected subspace; no nonlinear/global/PDE certificate')
