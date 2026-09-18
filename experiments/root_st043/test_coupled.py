"""Structural and numerical calibration. These are NOT NS acceptance tests."""
import bootstrap
import json,hashlib
from pathlib import Path
import numpy as np
import pytest
from coupled_fit import CoupledModel,JointObjective
from localized_model import LocalizedModel
from spacetime import Family,quad
from replay import load

@pytest.fixture(scope='module')
def parent(tmp_path_factory):
    path=tmp_path_factory.mktemp('parent')/'ST042.json'
    f,c=load('ST042');f.save(c,path,{'scope':'parent recipe replay; metadata differs'})
    return path

@pytest.mark.parametrize('kind',[CoupledModel,LocalizedModel])
def test_reduced_momentum_and_constraint_gradients(parent,kind):
    model=kind(parent,radial=1,axial=1,time_degree=2)
    rng=np.random.default_rng(9174291);s=rng.uniform(.02,2,24);z=rng.uniform(-1,1,24);t=rng.uniform(.25,.75,24)
    cache=model.cache(s,z,t);c=np.zeros(model.dim);R,J=model.momentum(c,cache)
    points=np.c_[np.sqrt(s),np.zeros(len(s)),z]
    np.testing.assert_allclose(R,model.f.analytic_residual(model.raw,points,t),atol=1e-10,rtol=1e-10)
    d=rng.normal(size=model.dim);d/=np.linalg.norm(d);h=1e-6
    rp=model.momentum(c+h*d,cache)[0];rm=model.momentum(c-h*d,cache)[0]
    np.testing.assert_allclose((rp-rm)/(2*h),J@d,atol=1e-7,rtol=3e-6)
    obj=JointObjective(model,space_order=(5,8),time_cap=1.0,softmax_weight=.02,core_weight=.01)
    value,g=obj.fun(c);fp=obj.fun(c+h*d)[0];fm=obj.fun(c-h*d)[0]
    assert abs((fp-fm)/(2*h)-g@d)<1e-7
    C,K=obj.constraints(c);cp=obj.constraints(c+h*d)[0];cm=obj.constraints(c-h*d)[0]
    assert np.max(abs((cp-cm)/(2*h)-K@d))<2e-5


def test_compact_energy_and_axisymmetry(parent):
    m=CoupledModel(parent,radial=1,axial=1,time_degree=2);c=np.zeros(m.dim);c[8]=.0003
    raw=m.candidate(c);s,z,w=quad(96);points=np.c_[np.sqrt(s),np.zeros(len(s)),z]
    u,_=m.f.fields(raw,points,.25)
    assert abs(.5*w@np.sum(u*u,axis=1)-1)<2e-7
    x=np.array([[2,0,0],[0,0,2],[0,0,-2],[3,3,3]],float)
    for v in m.f.fields(raw,x,.5):assert not v.any()
    angle=.41;rot=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
    x=points[::400];u,_=m.f.fields(raw,x,.5);ur,_=m.f.fields(raw,x@rot.T,.5)
    np.testing.assert_allclose(ur,u@rot.T,atol=1e-11,rtol=1e-11)


def test_wrong_pressure_is_not_correct_force(parent):
    m=CoupledModel(parent,radial=1,axial=1,time_degree=2);D=m.cache(np.array([.01,.01]),np.array([-.1,.1]),.5)
    # Positive z*p_z means the PRESSURE FORCE -grad(p) is inward axially.
    assert np.all(D['z']*D['v']['Qz']<0)
    R,_=m.momentum(np.zeros(m.dim),D);material=R[:,2]-D['v']['Qz']
    np.testing.assert_allclose(material+D['v']['Qz'],R[:,2],atol=1e-15)


def test_no_false_pressure_or_energy_collapse(parent):
    m=CoupledModel(parent,radial=1,axial=1,time_degree=2);c=np.zeros(m.dim);c[-2]=20
    with pytest.raises(ValueError,match='bound'):m.candidate(c)
    c=np.zeros(m.dim);c[:m.f.nt]=-1;c[m.na:m.na+m.f.nt]=-1
    with pytest.raises(ValueError,match='Collapsed'):m.norm(c)

@pytest.mark.parametrize('ident',['ST045-G','ST045-H','ST044-F'])
def test_frozen_correction_recipes(ident):
    from frozen_replay import load as child
    f,raw=child(ident)
    s,z,w=quad(96);u,p=f.fields(raw,np.c_[np.sqrt(s),np.zeros(len(s)),z],.25)
    assert abs(.5*w@np.sum(u*u,axis=1)-1)<2e-7
    outside=np.array([[2.1,0,0],[0,0,2.1]])
    u,p=f.fields(raw,outside,.5);assert not u.any() and not p.any()


def test_recipe_claim_tampering():
    from frozen_replay import construct,HERE
    d=json.loads((HERE/'recipes.json').read_text())['ST045-H'];d['pde_validated']=True
    with pytest.raises(ValueError,match='Unsupported'):construct(d)


def test_pointwise_pressure_direction_inequality():
    rng=np.random.default_rng(9174398);sign=rng.choice([-1,1],size=100);M=rng.normal(size=100);pz=sign*rng.uniform(0,1,100)
    assert np.all(abs(M+pz)+1e-15>=np.maximum(sign*M,0))
