"""Polynomial theta-momentum cache, checked against direct FD evaluation."""
import numpy as np
from .constrained_validation import residual
from .constrained_bipolar_theta import theta
from .constrained_force import RestrictedForce

def build_theta_objective(field, size, times):
    points=np.random.default_rng(9172611).uniform(-2,2,(512,3))
    directions=theta(points);eye=np.eye(size)
    def evaluate(coeff):
        f=field(coeff);chunks=[]
        for t in times:
            r=residual(f.at_points,lambda x,t:np.zeros(len(x)),lambda x,t:np.zeros_like(x),points,t,step=.005)
            chunks.append(np.sum(r['momentum']*directions,axis=1))
        return np.concatenate(chunks)
    positive=[evaluate(v) for v in eye];negative=[evaluate(-v) for v in eye]
    linear=np.stack([(p-m)/2 for p,m in zip(positive,negative)],axis=1)
    quad=np.zeros((len(linear),size,size))
    for i in range(size):quad[:,i,i]=(positive[i]+negative[i])/2
    for i in range(size):
        for j in range(i):
            quad[:,i,j]=quad[:,j,i]=(evaluate(eye[i]+eye[j])-positive[i]-positive[j])/2
    force=np.column_stack([np.concatenate([np.sum(RestrictedForce(a=a,c=c)(points,t)*directions,axis=1) for t in times]) for a,c in [(1,0),(0,1)]])
    probe=np.linspace(-.3,.4,size)
    error=float(np.max(np.abs(evaluate(probe)-(linear@probe+np.einsum('nij,i,j->n',quad,probe,probe)))))
    if error>1e-7:raise RuntimeError(f'Polynomial replay failure {error}')
    def loss(a):
        v=a[:-2];r=linear@v+np.einsum('nij,i,j->n',quad,v,v)-force@a[-2:]
        return float(64*np.mean(r*r))
    return loss,error
