import numpy as np
from spacetime import Family
from mixed_exact_newton import ProjectedExact

def test_projected_jacobian_and_hessian():
    f=Family(2,2,3);o=ProjectedExact(f,seed=9172900,count=120,weight=10.);raw=f.initial();rng=np.random.default_rng(44)
    y=o.initial(raw);y[:2*f.n]+=rng.normal(0,.02,2*f.n);y[-2:]=[.2,.3]
    r,J=o.linearize(y);J=J.copy();g=J.T@r;H,N=o.hessian(y)
    v=rng.normal(size=len(y));v/=np.linalg.norm(v);eps=2e-6
    rp,Jp=o.linearize(y+eps*v);gp=Jp.T@rp;rp=rp.copy()
    rm,Jm=o.linearize(y-eps*v);gm=Jm.T@rm
    jerr=np.linalg.norm((rp-rm)/(2*eps)-J@v)/np.linalg.norm(J@v)
    herr=np.linalg.norm((gp-gm)/(2*eps)-H@v)/np.linalg.norm(H@v)
    print('directional errors',jerr,herr)
    assert jerr<1e-6;assert herr<1e-6
    o.evaluate(y)
    assert np.max(np.abs(o.Q.T@o.res[:3*o.N]))<1e-10
