import json
import numpy as np
from asymmetric_basis import AsymmetricFamily
from asymmetric_exact_run import AllParityMomentObjective
from spacetime import Family, quad
import mixed_exact_newton as solver


def test_asymmetric_representation_preserves_energy_and_axisymmetry(tmp_path):
    f=AsymmetricFamily(2,10,3);x=f.initial();rng=np.random.default_rng(9172991)
    x[:2*f.n]+=rng.normal(0,.002,2*f.n)
    s,z,w=quad(96);points=np.column_stack((np.sqrt(s),np.zeros(len(s)),z))
    u,p=f.fields(x,points,.25)
    assert abs(.5*w@np.sum(u*u,axis=1)-1)<1e-8
    mirrored,_=f.fields(x,points*np.array([1,1,-1]),.25)
    assert np.max(np.abs(mirrored-u*np.array([1,1,-1])))>1e-4
    angle=.37;Q=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
    actual,_=f.fields(x,points[:80]@Q.T,.25)
    np.testing.assert_allclose(actual,u[:80]@Q.T,atol=1e-12,rtol=1e-12)
    fn=tmp_path/'asymmetric.json';f.save(x,fn);g,xx=Family.load(fn)
    np.testing.assert_allclose(g.fields(xx,points[:40],.25)[0],u[:40],atol=0,rtol=0)


def test_all_parity_projected_hessian(monkeypatch):
    monkeypatch.setattr(solver,'Objective',AllParityMomentObjective)
    f=AsymmetricFamily(2,10,3);o=solver.ProjectedExact(f,seed=9172992,count=360,weight=10.)
    rng=np.random.default_rng(9172993);y=o.initial(f.initial());y[:2*f.n]+=rng.normal(0,.003,2*f.n);y[-2:]=[.1,.2]
    r,J=o.linearize(y);H,N=o.hessian(y);v=rng.normal(size=len(y));v/=np.linalg.norm(v);h=1e-6
    rp,Jp=o.linearize(y+h*v);gp=Jp.T@rp
    rm,Jm=o.linearize(y-h*v);gm=Jm.T@rm
    relative=np.linalg.norm((gp-gm)/(2*h)-H@v)/np.linalg.norm(H@v)
    print('all-parity Hessian relative difference',relative)
    assert relative<1e-5
