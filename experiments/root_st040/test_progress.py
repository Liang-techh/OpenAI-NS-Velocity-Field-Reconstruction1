import bootstrap
from replay import load
import json
from pathlib import Path
import numpy as np
import pytest
from spacetime import Family
from controlled_fit import Model
from temporal_fit import TemporalModel,TemporalBiasModel
from structure_audit import jets
ROOT=Path(__file__).resolve().parent.parent
PARENT=bootstrap.PARENT

@pytest.fixture(scope='module')
def parent():return Family.load(PARENT)

def test_vorticity_and_pressure_independent_cartesian(parent):
 f,c=parent;rng=np.random.default_rng(9174090);r=rng.uniform(.05,.7,30);z=rng.uniform(-.7,.7,30);x=np.column_stack((r,np.zeros(30),z));t=.4375;v=jets(f,c,r,z,t)
 h=.0005;J=np.empty((len(x),3,3));pg=np.empty((len(x),3))
 for k in range(3):
  e=np.eye(3)[k]*h;up,pp=f.fields(c,x+e,t);um,pm=f.fields(c,x-e,t);J[:,:,k]=(up-um)/(2*h);pg[:,k]=(pp-pm)/(2*h)
 curl=np.column_stack((J[:,2,1]-J[:,1,2],J[:,0,2]-J[:,2,0],J[:,1,0]-J[:,0,1]))
 np.testing.assert_allclose(v[:,3:6],curl,atol=6e-6,rtol=2e-5)
 np.testing.assert_allclose(v[:,6:],np.column_stack((pg[:,0],pg[:,2],J[:,2,0])),atol=5e-6,rtol=2e-5)

def test_exact_reflection_missing_bias(parent):
 f,c=parent;v=jets(f,c,np.linspace(0,.6,20),0,.5)
 assert np.max(abs(v[:,2]))==0 and np.max(abs(v[:,8]))==0

def test_pressure_direction_is_force_not_gradient(parent):
 f,c=parent;z=np.array([-.1,.1]);v=jets(f,c,.1,z,.5)
 assert np.all(z*v[:,7]<0) # Fz=-p_z points away, not toward, z=0.
 assert np.all(v[:,6]>0) # radial F=-p_r is correctly inward.

@pytest.mark.parametrize('kind',[TemporalModel,TemporalBiasModel])
def test_reduced_momentum_jacobian_and_embedding(kind,parent):
 m=kind(PARENT);rng=np.random.default_rng(9174091);s=rng.uniform(.01,2,25);z=rng.uniform(-1,1,25);t=rng.uniform(.25,.75,25);K=m.cache(s,z,t);c=np.zeros(m.dim)
 R,J=m.momentum(c,K);f,raw=parent;x=np.column_stack((np.sqrt(s),np.zeros(len(s)),z))
 np.testing.assert_allclose(R,f.analytic_residual(raw,x,t),atol=1e-10,rtol=1e-10)
 d=rng.normal(size=m.dim);d[-2:]=0;d/=np.linalg.norm(d);h=1e-6
 rp=m.momentum(c+h*d,K)[0];rm=m.momentum(c-h*d,K)[0]
 rel=np.linalg.norm((rp-rm)/(2*h)-np.einsum('nij,j->ni',J,d))/np.linalg.norm(np.einsum('nij,j->ni',J,d))
 assert rel<2e-6

def test_children_are_not_falsely_validated():
 for name in ['ST041','ST042']:
  f,c=load(name);outside=np.array([[2.,0,0],[0,0,2],[3,0,3]])
  u,p=f.fields(c,outside,.5);assert not u.any() and not p.any()

def test_source_bias_is_real_not_changed_coordinate():
 f,c=load('ST042');r=np.linspace(.04,.4,25);v=jets(f,c,r,0,.5)
 assert np.min(v[:,2])>0 and np.max(abs(v[:,8]))>1e-3
 # Retains axisymmetry: rotation changes the Cartesian vector, not cylindrical profiles.
 points=np.column_stack((r,np.zeros(len(r)),np.full(len(r),.1)));a=.57;Q=np.array([[np.cos(a),-np.sin(a),0],[np.sin(a),np.cos(a),0],[0,0,1]])
 u,_=f.fields(c,points,.5);uq,_=f.fields(c,points@Q.T,.5)
 np.testing.assert_allclose(uq,u@Q.T,atol=1e-12,rtol=1e-12)

def test_symbolic_parity_obstruction():
 import sympy as s
 x,y,z,t=s.symbols('x y z t');S=s.symbols('S');H=s.Function('H')(S,z*z,t);F=z*H;C=2*F+2*S*s.diff(F,S)
 assert s.simplify(C.subs(z,0))==0
 assert s.simplify(s.diff(C,S).subs(z,0))==0
 A=-s.diff(F,z)
 assert s.simplify(2*A+2*S*s.diff(A,S)+s.diff(C,z))==0


def test_recipe_rejects_unknown_id():
 with pytest.raises(ValueError):load('made_up')

def test_recipe_fields_keep_energy_and_flags():
 from spacetime import quad
 for ident in ['ST041','ST042']:
  f,c=load(ident);s,z,w=quad(96);x=np.column_stack((np.sqrt(s),np.zeros(len(s)),z));u,p=f.fields(c,x,.25)
  assert abs(.5*w@np.sum(u*u,axis=1)-1)<1e-6
