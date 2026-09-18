"""Local algebra/numerical checks, not scientific NS acceptance."""
import json, hashlib
from pathlib import Path
import numpy as np
import pytest
from minimax_exchange import EdgeModel,Epigraph,pool_points,select_points,ROOT,PARENT_SHA
from spacetime import force,quad
from validate import cartesian_residual

@pytest.fixture(scope='module')
def model(tmp_path_factory):
    from replay_st052 import parent_field
    f,r=parent_field();p=tmp_path_factory.mktemp('parent')/'candidate.json';f.save(r,p)
    return EdgeModel(p,2,2,2,True)

def test_sampled_epigraph_jacobian(model):
    rng=np.random.default_rng(9175252);points=np.c_[rng.uniform(.02,3,21),rng.uniform(-1.85,1.85,21),rng.uniform(.25,.75,21)]
    e=Epigraph(model,None,points,.003,1.)
    c=np.zeros(model.dim);d=rng.normal(size=model.dim);d/=np.linalg.norm(d);h=1e-6
    v,J=e.constraints(c,1.)
    fd=(e.constraints(c+h*d,1.)[0]-e.constraints(c-h*d,1.)[0])/(2*h)
    assert np.max(abs(fd-J@d))<2e-6
    np.testing.assert_allclose((e.constraints(c,1+h)[0]-e.constraints(c,1-h)[0])/(2*h),1,atol=1e-8)

def test_epigraph_means_full_vector_peak(model):
    p=np.array([[.04,.17,.5],[.16,1.9,.75],[.04,-1.9,.25]])
    e=Epigraph(model,None,p,.003,1.);c=np.zeros(model.dim)
    R,J=model.momentum(c,e.D);q=float(np.max(np.sum(R*R,axis=1))/.003)
    v,_=e.constraints(c,q);assert v.min()>=-1e-14
    assert e.constraints(c,q-.001)[0].min()<0

def test_actual_residual_matches_independent_cartesian(model):
    x=np.array([[.2,0,.1],[0,0,.13],[.4,0,-1.4],[.57,0,1.92]])
    c=np.zeros(model.dim);raw=model.candidate(c)
    R,_=model.momentum(c,model.cache(x[:,0]**2,x[:,2],.5))
    field=lambda p,t:model.f.fields(raw,p,t)
    forcing=lambda p,t:force(p,t,*raw[-2:])
    fd,div=cartesian_residual(field,forcing,x,.5,.000625,.000625)
    assert np.max(np.linalg.norm(fd-R,axis=1))<1e-6

def test_normalization_support_and_rotation(model):
    raw=model.candidate(np.zeros(model.dim));s,z,w=quad(96)
    u,_=model.f.fields(raw,np.c_[np.sqrt(s),np.zeros(len(s)),z],.25)
    assert abs(.5*w@np.sum(u*u,axis=1)-1)<1e-6
    for v in model.f.fields(raw,np.array([[2.,0,0],[0,0,2.],[3,0,-3.]]),.5):assert not np.any(v)
    a=.31;Q=np.array([[np.cos(a),-np.sin(a),0],[np.sin(a),np.cos(a),0],[0,0,1]])
    x=np.array([[.1,.2,.3],[.7,.2,1.2]])
    np.testing.assert_allclose(model.f.fields(raw,x@Q.T,.5)[0],model.f.fields(raw,x,.5)[0]@Q.T,atol=1e-12)

def test_exchange_is_deterministic_and_cumulative():
    p,ng=pool_points();p2,ng2=pool_points();np.testing.assert_array_equal(p,p2)
    v=np.linalg.norm(p,axis=1);a=select_points(p,v,ng);b=select_points(p,-v,ng,a)
    assert len(a)==528 and len(b)==len(set(b)) and set(a).issubset(b)
    assert p[:,0].min()>=0 and p[:,0].max()<4
    assert np.all((p[:,2]>=.25)&(p[:,2]<=.75))

def test_three_symbolic_local_identities():
    import sympy as s
    a,b,rho=s.symbols('a b rho',real=True);R=s.Matrix(s.symbols('R0:3'));q=s.symbols('q',nonnegative=True)
    assert s.expand(R.dot(R)-(R[0]**2+R[1]**2+R[2]**2))==0
    J=s.Matrix([[a,-b,0],[b,a,0],[0,0,-2*a]])
    assert s.expand(s.trace(J*J)-(6*a*a-2*b*b))==0
    x,z=s.symbols('x z');F=s.Function('F')(x,z);A=-s.diff(F,z);C=2*F+2*x*s.diff(F,x)
    assert s.simplify(2*A+2*x*s.diff(A,x)+s.diff(C,z))==0
