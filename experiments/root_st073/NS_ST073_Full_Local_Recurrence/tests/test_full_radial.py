import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'upstream'))
import numpy as np
import pytest
from full_radial import *
from source_coordinates import coordinates
from local_field import independent_fd

@pytest.mark.parametrize('eta',[-.49,0,.37])
def test_lagrange_axis_values(eta):
 for b in [-2.,-1.01,-1.,-.505]:
  v=q_power_jet(b,eta,(6,4))
  assert abs(v[0,0]-1)<2e-16
  L=1-.01*eta**2;Qz=2*eta/L;Qt=1/L
  assert abs(v[1,0]-b*Qz)<2e-14
  assert abs(v[0,1]-b*Qt)<2e-14

def test_implicit_q_polynomial_identity():
 # Q equation verified by independent polynomial composition to low total degree.
 shape=(8,5);e=LD(.41);q=q_power_jet(1,e,shape);t=np.zeros(shape,dtype=LD);t[0,0]=1-e*e;t[0,1]=1
 z=np.zeros(shape,dtype=LD);z[0,0]=e;z[1,0]=1
 power=q_power_jet(.01,e,shape)
 d=q-multiply(multiply(z,z,shape),power,shape)-t
 assert np.max(abs(d))<1e-9

@pytest.mark.parametrize('N',[2,5])
def test_exact_straining_rotating_solution(N):
 shape=(2*N+6,N+5);z=np.zeros(shape,dtype=LD);z[0,0]=.13;z[1,0]=1
 alpha=LD(.6);omega=LD(.7);B=np.zeros(shape,dtype=LD)
 B[0]=omega*np.array([(-alpha)**j/np.math.factorial(j) for j in range(shape[1])]) if hasattr(np,'math') else omega*np.array([(-alpha)**j/__import__('math').factorial(j) for j in range(shape[1])])
 C=alpha*z;P=-alpha**2/2*multiply(z,z,shape)
 A,bb,cc,pp=recur(B,C,P,N,LD(1),LD(1),LD(1))
 for j in range(1,N+1):
  assert np.max(abs(cc[j][:3,:2]))<1e-15
  assert np.max(abs(bb[j][:3,:2]))<1e-15
  assert np.max(abs(A[j][:3,:2]))<1e-15
 # p_1=(B²-a²)/2 in s=r² coordinates
 expected=(multiply(B,B,shape)-multiply(A[0],A[0],shape))/2
 np.testing.assert_allclose(pp[1][:3,:3],expected[:3,:3],atol=1e-15,rtol=0)

@pytest.mark.parametrize('eta',[-.4,.2])
def test_axis_identity_and_normalization_not_changed(eta):
 f=FullRadialField(Parameters(order=4));tau=.0625;q=tau/(1-eta**2)
 d=f.evaluate_similarity(0,eta,tau)
 np.testing.assert_allclose(d['velocity'],[0,0,np.sqrt(.01)*q**(-f.A)*(4*eta+.02)],atol=1e-13)
 assert abs(d['pressure']-.01*q**(-2*f.A)*(-1+eta**2/2))<1e-12
 assert abs(d['divergence'])<1e-12

def test_divergence_rotation_support_guard():
 f=FullRadialField(Parameters(order=6));d=f.evaluate_similarity([.005,.01],[-.3,.2],.03)
 assert max(abs(d['divergence']))<1e-12
 r=f.evaluate_similarity([.005,.01],[-.3,.2],.03,angle=.72)
 Q=np.array([[np.cos(.72),-np.sin(.72),0],[np.sin(.72),np.cos(.72),0],[0,0,1]])
 np.testing.assert_allclose(r['velocity'],d['velocity']@Q.T,atol=1e-12)
 with pytest.raises(ValueError):f.evaluate_similarity(.016,0,.5)
 with pytest.raises(ValueError):f.evaluate_similarity(.01,.6,.5)
 with pytest.raises(ValueError):f.evaluate_similarity(.01,.1,.001)

def test_more_axis_terms_and_jets():
 a=FullRadialField(Parameters(order=6));b=FullRadialField(Parameters(order=6,axis_terms=160,jet_extra=5))
 for e in [-.5,.43]:
  x=a.evaluate_similarity(.015625,e,.0078125);y=b.evaluate_similarity(.015625,e,.0078125)
  np.testing.assert_allclose(x['velocity'],y['velocity'],rtol=1e-13,atol=1e-13)
  np.testing.assert_allclose(x['residual'],y['residual'],rtol=1e-7,atol=5e-9)

def test_viscosity_change_formula():
 f=FullRadialField(Parameters(order=4,nu=.01));g=FullRadialField(Parameters(order=4,nu=1))
 d=f.evaluate_similarity(.012,.23,.2);e=g.evaluate_similarity(.012,.23,.2)
 for key in ['velocity','residual']:np.testing.assert_allclose(d[key],.1*e[key],atol=1e-13)
 np.testing.assert_allclose(d['pressure'],.01*e['pressure'],atol=1e-13)

def test_independent_full_cartesian_operator():
 f=FullRadialField(Parameters(order=6));tau=.09
 pts=f.from_similarity(np.array([.006,.009]),np.array([-.2,.27]),tau,angle=[.33,.85]);d=f.evaluate(pts,tau)
 r,div=independent_fd(f,pts,tau,1e-5,tau*.001)
 np.testing.assert_allclose(r,d['residual'],rtol=0,atol=1e-5)
 assert max(abs(div))<1e-7

def test_correction_decreases_full_residual():
 vals=[]
 for n in [2,4,6,8]:
  f=FullRadialField(Parameters(order=n));r=f.evaluate_similarity(.015625,.25,.0078125)['residual'];vals.append(float(np.linalg.norm(r)))
 assert all(b<a for a,b in zip(vals,vals[1:])),vals
 assert vals[-1]<1e-5

def test_full_first_angular_coefficient():
 f=FullRadialField(Parameters(order=2,axis_swirl=4));tau=.03125
 co=f.coefficients(0.,tau);first=co[1,1,0]*tau**(1+f.h)/4
 expected=(f.h-3+2*(1+f.h)*tau**(2*f.h))/4
 assert abs(first-expected)<1e-13

def test_frozen_model_identities():
 import hashlib,json
 for fn in ['freeze.json','phase2_freeze.json']:
  d=json.loads((ROOT/'evidence'/fn).read_text());assert hashlib.sha256((ROOT/'full_radial.py').read_bytes()).hexdigest()==d['source_sha256']
  assert hashlib.sha256((ROOT/'data'/(d['id']+'.json')).read_bytes()).hexdigest()==d['model_sha256']

def test_selected_swirl_control_directions():
 f=FullRadialField.load(ROOT/'data/ST073-V.json');x=np.array([.003,.008,.013]);e=np.array([-.31,.08,.37]);d=f.evaluate_similarity(x,e,.07)
 assert np.all(d['velocity'][:,0]<0)
 assert np.all(d['velocity'][:,1]>0)
 assert np.all(d['velocity'][:,2]*e>0)
 assert np.all(d['pressure_gradient'][:,0]>0)
 assert np.all(d['pressure_gradient'][:,2]*e>0)
 assert max(np.linalg.norm(d['residual'],axis=1))<1e-5

def test_independent_audit_scope_and_rejected_original_pressure():
 import json
 for seed in [9237395,9237396]:
  rows=[json.loads(p.read_text()) for p in (ROOT/'evidence').glob(f'V_holdout_{seed}_*.json')]
  assert len(rows)==8 and sum(r['points'] for r in rows)==4096
  assert all(not r['global_pde_validated'] for r in rows)
  assert max(r['full_sampled_max'] for r in rows)<.001
 rows=json.loads((ROOT/'evidence/interface_audit.json').read_text())['rows']
 assert all(r['v_at_midplane']<2 for r in rows)
 assert all(r['inward_pressure_fraction']==0 for r in rows if r['id']=='ST073-F')

def test_evaluator_independent_inputs_and_roundtrip(tmp_path):
 f=FullRadialField.load(ROOT/'data/ST073-V.json');p=tmp_path/'model.json';f.save(p,'ST073-V');g=FullRadialField.load(p)
 d=f.evaluate_similarity(.01,.2,.1);e=g.evaluate_similarity(.01,.2,.1)
 np.testing.assert_array_equal(d['residual'],e['residual'])
 with pytest.raises(FileExistsError):f.save(p)
