"""New tests calibrate additions, not a claim of NS acceptance."""
from pathlib import Path
import hashlib
import numpy as np
import pytest
from boundary_shear import BoundaryShearObjective, LocalizedModel
ROOT=Path(__file__).resolve().parents[2]
PARENT=ROOT/'artifacts/research/ST047-E/candidate.json'

@pytest.fixture(scope='module')
def obj(tmp_path_factory):
    from replay_st048 import parent_field
    f,raw=parent_field();path=tmp_path_factory.mktemp('parent')/'candidate.json';f.save(raw,path)
    m=LocalizedModel(path,2,2,2)
    return BoundaryShearObjective(m,acceleration_ratio=.98,axis_weight=.02,profile_ratio=.995,pressure_ratio=.99,space_order=(6,8),peak_weight=0.,softmax_weight=.02,anchor_tolerance=.08,time_cap=1.)

def test_shear_jacobian(obj):
    rng=np.random.default_rng(9174891);c=rng.normal(0,1e-5,obj.m.dim);d=rng.normal(size=obj.m.dim);d/=np.linalg.norm(d);h=2e-6
    v,J=obj.shear(c);fd=(obj.shear(c+h*d)[0]-obj.shear(c-h*d)[0])/(2*h)
    np.testing.assert_allclose(fd,J@d,atol=2e-10,rtol=2e-6)

def test_shear_matches_independent_field_difference(obj):
    m=obj.m;c=np.zeros(m.dim);s=obj.shear_cache['s'];t=obj.shear_cache['t'];x=np.c_[np.sqrt(s),np.zeros(len(s)),np.zeros(len(s))];h=2e-5;e=np.array([h,0.,0.])
    raw=m.candidate(c);fd=(m.f.fields(raw,x+e,t)[0][:,2]-m.f.fields(raw,x-e,t)[0][:,2])/(2*h)
    np.testing.assert_allclose(fd,obj.shear(c)[0],atol=2e-9,rtol=2e-5)

def test_shear_independent_of_pressure_and_force(obj):
    m=obj.m;c=np.zeros(m.dim);v,_=obj.shear(c);c[m.nv:-2]=.001;c[-2:]=.002
    vv,J=obj.shear(c);np.testing.assert_array_equal(v,vv)
    assert np.max(abs(J[:,m.nv:]))==0.

def test_augmented_objective_and_constraints_jacobian(obj):
    rng=np.random.default_rng(9174892);c=rng.normal(0,1e-6,obj.m.dim);d=rng.normal(size=obj.m.dim);d/=np.linalg.norm(d);h=1e-6
    val,g=obj.fun(c);v,J=obj.constraints(c)
    fd=(obj.fun(c+h*d)[0]-obj.fun(c-h*d)[0])/(2*h)
    kd=(obj.constraints(c+h*d)[0]-obj.constraints(c-h*d)[0])/(2*h)
    assert abs(fd-g@d)<1e-7
    np.testing.assert_allclose(kd,J@d,atol=3e-6,rtol=5e-5)

def test_edge_grid_includes_endpoints_and_axis(obj):
    s,z,t=obj.edge_points.T
    assert np.min(s)==0 and t.min()==.25 and t.max()==.75
    assert np.max(abs(z))==1.96 and np.sqrt(s).max()==1.96
    assert np.all((s>=0)&(s<4)&(abs(z)<2))

def test_parent_satisfies_new_signed_shear_safeguards(obj):
    c=np.zeros(obj.m.dim);v,_=obj.constraints(c);n=len(obj.shear_parent)
    assert np.min(v[-2*n:])>=0
    assert np.min(abs(obj.shear_parent))>0

def test_actual_field_support_and_original_initial_energy(obj):
    from spacetime import quad
    m=obj.m;raw=m.candidate(np.zeros(m.dim));s,z,w=quad(96);x=np.c_[np.sqrt(s),np.zeros(len(s)),z]
    u,p=m.f.fields(raw,x,.25);assert abs(.5*w@np.sum(u*u,axis=1)-1)<1e-7
    x=np.array([[2.,0,0],[0,0,2.],[3,2,4]])
    for a in m.f.fields(raw,x,.5):assert not a.any()

def test_symbolic_structure():
    from structural_identities import check
    assert check()['checks_passed']==3

def test_random_time_cartesian_operator_known_linear_flow():
    from extra_audit import spacetime_fd
    class Linear:
        def fields(self,raw,x,t):
            return np.asarray(t)[:,None]*x*np.array([.2,.3,-.5]),np.zeros(len(x))
    rng=np.random.default_rng(9174893);x=rng.uniform(-1,1,(20,3));t=rng.uniform(.3,.7,20);d=np.array([.2,.3,-.5])
    R,div=spacetime_fd(Linear(),np.zeros(2),x,t)
    expected=x*d+t[:,None]**2*x*d*d
    np.testing.assert_allclose(R,expected,atol=1e-10,rtol=1e-9)
    assert np.max(abs(div))<1e-11 and np.max(abs(R))>.1


def test_random_time_operator_rejects_boundary_crossing(obj):
    from extra_audit import spacetime_fd
    with pytest.raises(ValueError,match='central time'):
        spacetime_fd(obj.m.f,obj.m.raw,np.zeros((1,3)),np.array([.25]))
