import json
import numpy as np
import pytest
from spacetime import Family, force, quad
from validate import cartesian_residual

def test_calibrated_manufactured_shear_and_force_mutation():
    rng=np.random.default_rng(900);x=rng.uniform(-1,1,(50,3))
    def field(x,t):
        return np.column_stack((t**5*np.sin(x[:,1]),np.zeros((len(x),2)))),np.zeros(len(x))
    def forcing(x,t):
        return np.column_stack(((5*t**4+.01*t**5)*np.sin(x[:,1]),np.zeros((len(x),2))))
    for t in (.25,.5,.75):
        errors=[]
        for h in (.02,.01,.005):
            r,d=cartesian_residual(field,forcing,x,t,.002,h);errors.append(np.max(np.abs(r)))
            assert np.max(np.abs(d))<1e-12
        assert 12<errors[0]/errors[1]<20
        assert 12<errors[1]/errors[2]<20
    r,_=cartesian_residual(field,lambda x,t:-forcing(x,t),x,.5,.002,.005)
    assert np.max(np.abs(r))>.1

def test_energy_symmetry_support_and_roundtrip(tmp_path):
    f=Family(3,3,4);raw=f.initial();s,z,w=quad(96)
    x=np.column_stack((np.sqrt(s),np.zeros(len(s)),z));u,p=f.fields(raw,x,.25)
    assert abs(.5*w@np.sum(u*u,axis=1)-1)<1e-9
    um,_=f.fields(raw,x*np.array([1,1,-1]),.25)
    np.testing.assert_allclose(um,u*np.array([1,1,-1]),atol=1e-13)
    xb=np.array([[2,0,0],[0,0,2],[3,0,0],[0,0,-3]],float)
    for v in f.fields(raw,xb,.5): assert np.max(np.abs(v))==0
    path=tmp_path/'candidate.json';f.save(raw,path);ff,rr=Family.load(path)
    np.testing.assert_allclose(ff.fields(rr,x[:12],.5)[0],f.fields(raw,x[:12],.5)[0],atol=0,rtol=0)
    p=json.loads(path.read_text());p['pde_validated']=True;path.write_text(json.dumps(p))
    with pytest.raises(ValueError): Family.load(path)

def test_analytic_residual_matches_independent_fd():
    rng=np.random.default_rng(901);f=Family(3,3,4);raw=f.initial()
    raw[:3*f.n]+=rng.normal(0,.002,3*f.n);raw[-2:]=[.13,.21]
    x=rng.uniform(-1.5,1.5,(80,3));x[:3]=[[0,0,0],[0,0,.3],[0,0,-.4]]
    field=lambda x,t:f.fields(raw,x,t);forcing=lambda x,t:force(x,t,*raw[-2:])
    for t in (.25,.5,.75):
        exact=f.analytic_residual(raw,x,t);errors=[];div=[]
        for h in (.02,.01,.005):
            r,d=cartesian_residual(field,forcing,x,t,h,.001)
            errors.append(np.max(np.linalg.norm(r-exact,axis=1)));div.append(np.max(np.abs(d)))
        assert errors[-1]<3e-5
        assert errors[0]/errors[1]>10
        assert div[0]/div[1]>10

def test_force_is_curl_of_preregistered_potential():
    from spacetime import bump_derivatives
    rng=np.random.default_rng(902);x=rng.uniform(-1.5,1.5,(20,3));a,c=.4,.8;t=.42;h=1e-5
    def potential(x):
        xx,yy,z=x.T;r2=xx*xx+yy*yy
        b=bump_derivatives(r2/4)[0]*bump_derivatives(z*z/4)[0]*bump_derivatives((2*t-1)**2)[0]
        return b[:,None]*np.column_stack((-a*z*yy,a*z*xx,-c*r2/2))
    jac=np.empty((len(x),3,3))
    for i in range(3):
        e=np.eye(3)[i]*h;jac[:,:,i]=(potential(x+e)-potential(x-e))/(2*h)
    curl=np.column_stack((jac[:,2,1]-jac[:,1,2],jac[:,0,2]-jac[:,2,0],jac[:,1,0]-jac[:,0,1]))
    np.testing.assert_allclose(force(x,t,a,c),curl,atol=1e-8)

def test_weak_pressure_independent_identity():
    f=Family();raw=f.initial();s,z,w=quad(96);x=np.column_stack((np.sqrt(s),np.zeros(len(s)),z))
    u,p=f.fields(raw,x,.5);R=f.analytic_residual(raw,x,.5)
    weak=w@(.5*np.sqrt(s)*R[:,0]-z*R[:,2])
    anisotropy=w@(u[:,2]**2-.5*(u[:,0]**2+u[:,1]**2))
    assert abs(weak-anisotropy)<1e-9

def test_fail_closed_invalid_inputs():
    f=Family();r=f.initial()
    with pytest.raises(ValueError):f.fields(r,np.zeros((1,3)),.8)
    with pytest.raises(ValueError):f.coefficients(np.zeros_like(r))
    with pytest.raises(ValueError):force(np.zeros((1,3)),.5,-1,0)


def test_gauss_newton_jacobian_direction():
    pytest.importorskip('torch')
    from gauss_newton import Objective
    f=Family(3,3,4);raw=f.initial();rng=np.random.default_rng(9124)
    raw[:3*f.n]+=rng.normal(0,.005,3*f.n);raw[-2:]=[.3,.4]
    obj=Objective(f,count=300);J=obj.jac(raw).copy()
    direction=rng.normal(size=len(raw));direction/=np.linalg.norm(direction)
    h=1e-5
    finite=(obj.fun(raw+h*direction)-obj.fun(raw-h*direction))/(2*h)
    assert np.max(np.abs(finite-J@direction))<1e-7
