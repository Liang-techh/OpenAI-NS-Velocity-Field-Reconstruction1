"""Calibration of the harmonic witness basis, not a PDE acceptance test."""
import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from harmonic_audit import basis

def test_exact_gram_independent_quadrature_and_nested_bounds():
    H,G,_=basis();s,z=sorted(H[0].free_symbols,key=str)
    assert abs(G[0,0]-352*np.pi/3)<1e-10
    fun=sp.lambdify((s,z),[(sp.diff(h,s),sp.diff(h,z)) for h in H],'numpy')
    x,w=leggauss(16);S,Z=np.meshgrid(2*(x+1),2*x,indexing='ij');s,z=S.ravel(),Z.ravel();r=np.sqrt(s);w=(4*np.pi*np.outer(w,w)).ravel()
    values=fun(s,z);pr=np.array([2*r*np.broadcast_to(a,s.shape) for a,b in values]).T;pz=np.array([np.broadcast_to(b,s.shape) for a,b in values]).T
    actual=pr.T@(w[:,None]*pr)+pz.T@(w[:,None]*pz)
    sc=np.sqrt(np.diag(G));N=G/sc[:,None]/sc[None,:]
    np.testing.assert_allclose(actual/sc[:,None]/sc[None,:],N,atol=2e-13,rtol=2e-13)
    assert np.linalg.eigvalsh(N).min()>0
    d=np.arange(1.,len(H)+1)/sc;bounds=[]
    for k in range(1,len(H)+1):
        c=np.linalg.solve(N[:k,:k],d[:k]);bounds.append(np.sqrt(d[:k]@c))
    assert np.min(np.diff(bounds))>=-1e-14
