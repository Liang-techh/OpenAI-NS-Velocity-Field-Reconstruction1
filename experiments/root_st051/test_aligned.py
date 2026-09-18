"""Numerical/symbolic calibration only; these tests do not accept an NS field."""
from pathlib import Path
import sys
import numpy as np
import pytest
from aligned_continuation import EdgeModel,AlignedObjective,Family
ROOT=Path(__file__).resolve().parents[2]

def parents(tmp_path_factory):
    sys.path.insert(0,str(ROOT/'experiments/root_st050r'))
    from replay_recovery import reconstruct,parent_field
    d=tmp_path_factory.mktemp('parents');paths=[]
    for key in ('ST050R-P','ST048-S'):
        p=ROOT/f'artifacts/research/{key}/candidate.json'
        if p.exists():paths.append(p);continue
        f,r=parent_field() if key=='ST048-S' else reconstruct(key)
        p=d/f'{key}.json';f.save(r,p);paths.append(p)
    return paths

@pytest.fixture(scope='module')
def setup(tmp_path_factory):
    p,s=parents(tmp_path_factory);m=EdgeModel(p,2,2,2,True)
    o=AlignedObjective(m,shear_reference=s,poisson_weight=.003,axis_ratio=1.75,pressure_target=0.,profile_ratio=.995,pressure_ratio=1.,acceleration_ratio=1.,space_order=(12,18),anchor_tolerance=.12,time_cap=1.)
    return m,o

def test_edge_embedding_and_bounds(setup):
    m,_=setup;r=m.candidate(np.zeros(m.dim));x=np.array([[.1,.0,.2],[0,0,.1],[.2,0,1.9]])
    u,p=m.f.fields(r,x,.5);a,b=m.f.fields(m.raw,x,.5)
    np.testing.assert_allclose(u,a,atol=2e-12,rtol=2e-12)
    np.testing.assert_allclose(p,b,atol=2e-12,rtol=2e-12)
    assert m.dim==146 and len(m.bounds)==m.dim and len(m.edge_pairs)==10

def test_complete_momentum_jacobian(setup):
    m,_=setup;rng=np.random.default_rng(9175152);c=rng.normal(0,1e-4,m.dim);c[-2:]=0
    x=np.c_[rng.uniform(.02,1.6,25),rng.uniform(-1.6,1.6,25),rng.uniform(.25,.75,25)]
    D=m.cache(x[:,0]**2,x[:,1],x[:,2]);v,J=m.momentum(c,D);d=rng.normal(size=m.dim);d/=np.linalg.norm(d);h=1e-6
    fd=(m.momentum(c+h*d,D)[0]-m.momentum(c-h*d,D)[0])/(2*h)
    np.testing.assert_allclose(fd,J@d,atol=2e-8,rtol=2e-6)

def test_robust_constraints_and_loss_jacobian(setup):
    m,o=setup;c=np.zeros(m.dim);rng=np.random.default_rng(9175153);d=rng.normal(size=m.dim);d/=np.linalg.norm(d);h=1e-6
    v,J=o.constraints(c);g=o.fun(c)[1];fd=(o.constraints(c+h*d)[0]-o.constraints(c-h*d)[0])/(2*h)
    assert np.max(abs(fd-J@d))<2e-5
    df=(o.fun(c+h*d)[0]-o.fun(c-h*d)[0])/(2*h)
    assert abs(df-g@d)<2e-6

def test_new_shear_is_cartesian_derivative(setup):
    m,o=setup;c=np.zeros(m.dim);v,J=o.robust_shear(c);raw=m.candidate(c);D=o.rs
    indices=np.arange(0,len(v),37);r=np.sqrt(D['s'][indices]);t=D['t'][indices]
    pts=np.c_[r,np.zeros(len(r)),np.zeros(len(r))];h=2e-5;e=np.array([h,0,0])
    fd=(m.f.fields(raw,pts+e,t)[0][:,2]-m.f.fields(raw,pts-e,t)[0][:,2])/(2*h)
    np.testing.assert_allclose(fd,v[indices],atol=1e-7,rtol=1e-5)

def test_energy_support_rotation(setup):
    from spacetime import quad
    m,_=setup;rng=np.random.default_rng(9175154);c=rng.normal(0,1e-4,m.dim);c[-2:]=0;raw=m.candidate(c)
    s,z,w=quad(96);x=np.c_[np.sqrt(s),np.zeros(len(s)),z];u,_=m.f.fields(raw,x,.25)
    assert abs(.5*w@np.sum(u*u,axis=1)-1)<1e-6
    for v in m.f.fields(raw,np.array([[2.,0,0],[0,0,2.],[3,0,0],[0,0,-2.]]),.5):assert np.max(abs(v))==0
    x=np.array([[.13,.09,.12],[.4,0,-.8]]);a=.31;Q=np.array([[np.cos(a),-np.sin(a),0],[np.sin(a),np.cos(a),0],[0,0,1]])
    np.testing.assert_allclose(m.f.fields(raw,x@Q.T,.5)[0],m.f.fields(raw,x,.5)[0]@Q.T,atol=2e-12)

def test_symbolic_structural_identities():
    import sympy as s
    r,z,t=s.symbols('r z t',real=True);q=s.symbols('q',real=True);F=s.Function('F')(q,z,t)
    A=-s.diff(F,z);C=2*F+2*q*s.diff(F,q)
    assert s.simplify(2*A+2*q*s.diff(A,q)+s.diff(C,z))==0
    a,b,az,cz=s.symbols('a b az cz');J=s.Matrix([[a,-b,0],[b,a,0],[0,0,-2*a]])
    assert s.expand(s.trace(J*J)-(6*a*a-2*b*b))==0
    lam,M,N=s.symbols('lam M N',nonzero=True)
    assert s.simplify((lam**2*M)/(lam**2*N)-M/N)==0
    M,P=s.symbols('M P',nonnegative=True)
    assert s.simplify((M+P)-M)==P

def test_frozen_replay_and_mutations():
    import json,copy
    from replay_st051 import reconstruct
    folder=Path(__file__).with_name('recipes')
    records={p.stem:json.loads(p.read_text()) for p in folder.glob('*.json')}
    assert records
    for ident in records:
        f,r=reconstruct(ident,records);assert np.isfinite(r).all()
        bad=copy.deepcopy(records);bad[ident]['pde_validated']=True
        with pytest.raises(ValueError,match='scientific claim'):reconstruct(ident,bad)
        bad=copy.deepcopy(records);bad[ident]['modifiers_sha256']='0'*64
        with pytest.raises(ValueError,match='hash mismatch'):reconstruct(ident,bad)
        bad=copy.deepcopy(records);bad[ident]['parent_id']='invented'
        with pytest.raises(ValueError,match='Parent identity'):reconstruct(ident,bad)
