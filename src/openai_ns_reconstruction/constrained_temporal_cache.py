"""Cache the exact quadratic parameter dependence of the training FD operator."""
import numpy as np
from .constrained_optimize import training_residual


def cached_residual(candidate, points, times, nu=.01, step=.001):
    """Return coefficient-to-residual map for additive linear velocity corrections.

    Same centered differences as training_residual; this is an optimization cache,
    never an independent validation operator. Candidate must have zero coefficients.
    """
    if np.any(candidate.coefficients):raise ValueError('cache requires zero corrections')
    x=np.asarray(points);t=np.asarray(times);h=step
    base=training_residual(candidate,candidate.force,x,t,nu,h)
    u=candidate.velocity(x,t);B=candidate.correction_basis(x,t)
    linear=(candidate.correction_basis(x,t+h)-candidate.correction_basis(x,t-h))/(2*h)
    derivatives=[]
    for j in range(3):
        d=np.eye(3)[j]*h
        bp=candidate.correction_basis(x+d,t);bm=candidate.correction_basis(x-d,t)
        db=(bp-bm)/(2*h);derivatives.append(db)
        du=(candidate.velocity(x+d,t)-candidate.velocity(x-d,t))/(2*h)
        linear+=u[:,j,None,None]*db+du[:,:,None]*B[:,j,None,:]-nu*(bp-2*B+bm)/h**2
    def evaluate(coefficients):
        a=np.asarray(coefficients)
        delta=B@a
        result=base+linear@a
        for j,db in enumerate(derivatives):result+=delta[:,j,None]*(db@a)
        return result
    def jacobian(coefficients):
        a=np.asarray(coefficients);delta=B@a
        result=linear.copy()
        for j,db in enumerate(derivatives):
            result+=B[:,j,None,:]*(db@a)[:,:,None]+delta[:,j,None,None]*db
        return result
    evaluate.jacobian=jacobian
    return evaluate


def cached_energy(candidate,times,order=24):
    """Same cylindrical energy quadrature, with exact coefficient derivatives."""
    from .quadrature import unit_rule
    if np.any(candidate.coefficients):raise ValueError('cache requires zero corrections')
    n,w=unit_rule(order);rr,zz=np.meshgrid(2*n,4*n-2,indexing='ij')
    p=np.stack((rr,np.zeros_like(rr),zz),axis=-1)
    weight=(8*np.pi*rr*w[:,None]*w[None,:]).ravel()
    base=np.array([candidate.velocity(p,t).reshape(-1,3) for t in times])
    basis=np.array([candidate.correction_basis(p,t).reshape(-1,3,len(candidate.coefficients)) for t in times])
    def evaluate(a):
        velocity=base+basis@np.asarray(a)
        return np.einsum('tnc,tnc,n->t',velocity,velocity,weight)
    def jacobian(a):
        velocity=base+basis@np.asarray(a)
        return 2*np.einsum('tnc,tnck,n->tk',velocity,basis,weight)
    evaluate.jacobian=jacobian
    return evaluate
