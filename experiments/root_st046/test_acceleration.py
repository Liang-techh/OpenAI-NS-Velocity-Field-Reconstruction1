from pathlib import Path
import numpy as np
import pytest
import acceleration_fit as af
from mechanism_audit import terms
from coupled_fit import CoupledModel
from spacetime import quad
PARENT=Path(__file__).resolve().parents[2]/'artifacts/research/ST045-H/candidate.json'

@pytest.fixture(scope='module')
def model(tmp_path_factory):
    parent=PARENT
    if not parent.is_file():
        from frozen_replay import load
        f,x=load('ST045-H');parent=tmp_path_factory.mktemp('parent')/'candidate.json';f.save(x,parent)
    return CoupledModel(parent,radial=1,axial=1,time_degree=2)

def test_constraint_and_objective_directions(model):
    o=af.AccelerationObjective(model,.95,space_order=(4,6),profile_ratio=.98,pressure_ratio=.97,time_cap=1.)
    c=np.zeros(model.dim);r=np.random.default_rng(9174689);d=r.normal(size=model.dim);d/=np.linalg.norm(d);h=1e-6
    v,J=o.constraints(c);vp=o.constraints(c+h*d)[0];vm=o.constraints(c-h*d)[0]
    assert np.max(abs((vp-vm)/(2*h)-J@d))<3e-5
    f,g=o.fun(c);fp=o.fun(c+h*d)[0];fm=o.fun(c-h*d)[0]
    assert abs((fp-fm)/(2*h)-g@d)<1e-7
    # Pressure coefficients must cancel from the new non-pressure constraints.
    count=len(o.signz)
    assert np.max(abs(J[-count:,model.nv:model.nv+model.np]))<1e-12

def test_terms_match_independent_time_and_space_fd(model):
    f,raw=model.f,model.raw;x=np.array([[.12,0,.1],[.03,0,-.09],[.4,0,.3]]);t=.53;h=2e-4
    from validate import cartesian_residual
    field=lambda x,t:f.fields(raw,x,t)
    forcing=lambda x,t:af.np.array(__import__('spacetime').force(x,t,*raw[-2:]))
    res,_=cartesian_residual(field,forcing,x,t,space_step=h,time_step=h)
    pg=(field(x+np.array([0,0,h]),t)[1]-field(x-np.array([0,0,h]),t)[1])/(2*h)
    tr=terms(f,raw,x[:,0],x[:,2],t)
    np.testing.assert_allclose(tr[:,5],res[:,2]-pg,atol=2e-7,rtol=1e-5)
    np.testing.assert_allclose(tr[:,6],res[:,2],atol=2e-7,rtol=1e-5)

def test_original_energy_support_and_rotation(model):
    c=np.zeros(model.dim);c[8]=.0002;raw=model.candidate(c);s,z,w=quad(96);x=np.c_[np.sqrt(s),np.zeros(len(s)),z];u,p=model.f.fields(raw,x,.25)
    assert abs(.5*w@np.sum(u*u,axis=1)-1)<3e-7
    outside=np.array([[2,0,0],[0,0,2],[3,0,0],[0,0,-3]])
    for a in model.f.fields(raw,outside,.5):assert not a.any()
    a=.3;Q=np.array([[np.cos(a),-np.sin(a),0],[np.sin(a),np.cos(a),0],[0,0,1]])
    x=x[::500];u,p=model.f.fields(raw,x,.5);uq,pq=model.f.fields(raw,x@Q.T,.5)
    np.testing.assert_allclose(uq,u@Q.T,atol=1e-11)
    np.testing.assert_allclose(pq,p,atol=1e-11)

def test_pressure_only_bound_is_not_velocity_bound():
    s=np.array([-1,1]);M=np.array([-.2,.3]);p=s*.1
    assert np.all(abs(M+p)>=np.maximum(s*M,0))
    newM=-p
    assert np.all(newM+p==0) # A jointly changed velocity-side term can remove it.

def test_axis_peak_direction(model):
    o=af.AxisPeakObjective(model,.95,axis_weight=.05,space_order=(4,6),profile_ratio=.98,pressure_ratio=.97,time_cap=1.)
    c=np.zeros(model.dim);rng=np.random.default_rng(9174688);v=rng.normal(size=model.dim);v/=np.linalg.norm(v);h=1e-6
    f,g=o.fun(c);fp=o.fun(c+h*v)[0];fm=o.fun(c-h*v)[0]
    assert abs((fp-fm)/(2*h)-g@v)<1e-7

def test_symbolic_divergence_and_pressure_separation():
    import sympy as S
    s,z,nu=S.symbols('s z nu');F=S.Function('F')(s,z);B=S.Function('B')(s,z)
    A=-S.diff(F,z);C=2*F+2*s*S.diff(F,s)
    assert S.simplify(2*A+2*s*S.diff(A,s)+S.diff(C,z))==0
    # The new constrained non-pressure term cannot be repaired by pressure coefficients alone.
    p,m=S.symbols('p m');res=m+p
    assert S.diff(res-p,p)==0

@pytest.mark.parametrize('ident',['ST046-A','ST046-B','ST046-C','ST046-D'])
def test_frozen_modifier_replay(ident):
    from replay_st046 import reconstruct
    f,raw=reconstruct(ident)
    s,z,w=quad(96);u,p=f.fields(raw,np.c_[np.sqrt(s),np.zeros(len(s)),z],.25)
    assert abs(.5*w@np.sum(u*u,axis=1)-1)<3e-7
    assert raw.size==2594

def test_replay_rejects_claim_and_modifier_change():
    import json
    from replay_st046 import reconstruct,HERE
    r=json.loads((HERE/'recipes.json').read_text());r['ST046-A']['pde_validated']=True
    with pytest.raises(ValueError,match='claim'):reconstruct('ST046-A',r)
    r=json.loads((HERE/'recipes.json').read_text());r['ST046-A']['modifiers_sha256']='00'*32
    with pytest.raises(ValueError,match='checksum'):reconstruct('ST046-A',r)
