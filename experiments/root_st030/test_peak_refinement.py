import numpy as np
from peak_refinement import PeakObjective
from spacetime import Family

def test_quartic_residual_jacobian():
    f=Family(2,2,3);rng=np.random.default_rng(333)
    pts=np.column_stack((rng.uniform(.05,3.5,120),rng.uniform(-1.8,1.8,120),rng.uniform(.25,.75,120)))
    o=PeakObjective(f,pts);x=f.initial();x[:3*f.n]+=rng.normal(0,.001,3*f.n);x[-2:]=[.2,.1]
    r,J=o.transformed(x);d=rng.normal(size=len(x));d/=np.linalg.norm(d);h=2e-6
    rp,_=o.transformed(x+h*d);rm,_=o.transformed(x-h*d)
    err=np.linalg.norm((rp-rm)/(2*h)-J@d)/np.linalg.norm(J@d)
    print('quartic directional Jacobian relative error',err)
    assert err<1e-6
