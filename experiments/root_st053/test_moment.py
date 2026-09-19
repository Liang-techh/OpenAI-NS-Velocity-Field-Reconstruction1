"""Local mathematical and numerical checks, not scientific candidate acceptance."""
from pathlib import Path
import sys
import numpy as np
import pytest
from moment_step import moment_matrices,weak_moment,training_pool,NORM2,ROOT
from aligned_continuation import EdgeModel
from spacetime import Family,force

@pytest.fixture(scope='module')
def model(tmp_path_factory):
    from replay_st053 import parent_field
    f,r=parent_field();p=tmp_path_factory.mktemp('parent')/'candidate.json';f.save(r,p)
    return EdgeModel(p,2,2,2,False)

def test_derivative(model):
    m=model;G=moment_matrices(m,[.25,.5,.75],(24,32));c=np.zeros(m.dim)
    v=np.random.default_rng(9175357).normal(size=m.dim);v/=np.linalg.norm(v);h=1e-6
    D,J=weak_moment(m,c,G)
    fd=(weak_moment(m,c+h*v,G)[0]-weak_moment(m,c-h*v,G)[0])/(2*h)
    np.testing.assert_allclose(fd,J@v,rtol=1e-7,atol=1e-8)

def test_pressure_force_modifiers_do_not_change_moment(model):
    m=model;G=moment_matrices(m,[.5],(16,24));c=np.zeros(m.dim)
    c[m.nv:]=np.random.default_rng(9175358).normal(0,.001,m.dim-m.nv)
    np.testing.assert_array_equal(weak_moment(m,c,G)[0],weak_moment(m,np.zeros(m.dim),G)[0])
    assert np.max(abs(weak_moment(m,c,G)[1][:,m.nv:]))==0

def integrate(f,raw,nr,nz,t):
    from numpy.polynomial.legendre import leggauss
    x,wx=leggauss(nr);z,wz=leggauss(nz);S,Z=np.meshgrid(2*(x+1),2*z,indexing='ij')
    s,z=S.ravel(),Z.ravel();w=(4*np.pi*np.outer(wx,wz)).ravel();pts=np.c_[np.sqrt(s),np.zeros(len(s)),z]
    U=[];R=[]
    for i in range(0,len(s),256):
        U.append(f.fields(raw,pts[i:i+256],t)[0]);R.append(f.analytic_residual(raw,pts[i:i+256],t))
    u=np.concatenate(U);r=np.concatenate(R);phi=pts*np.array([.5,.5,-1.])
    D=w@(u[:,2]**2-.5*(u[:,0]**2+u[:,1]**2))
    lhs=w@np.sum(phi*r,axis=1);F=force(pts,t,*raw[-2:]);force_moment=w@np.sum(phi*F,axis=1)
    E=.5*w@np.sum(u*u,axis=1);alpha=(w@(u[:,2]**2))/(2*E)
    return np.array([D,lhs,force_moment,E,alpha,np.sqrt(w@np.sum(r*r,axis=1))])

def test_weak_identity_and_quadrature(model):
    a=integrate(model.f,model.raw,48,72,.5);b=integrate(model.f,model.raw,64,96,.5);c=integrate(model.f,model.raw,80,120,.5)
    assert abs(c[0]-c[1])<1e-7,(a,b,c)
    assert abs(c[2])<1e-9
    assert abs(c[0]-c[3]*(3*c[4]-1))<1e-12
    assert abs(b[0]-c[0])<1e-7
    assert abs(c[0])/np.sqrt(NORM2) <= c[5]+1e-10

def test_original_initial_energy_and_support(model):
    a=integrate(model.f,model.raw,64,96,.25);assert abs(a[3]-1)<1e-6
    for v in model.f.fields(model.raw,np.array([[2.,0,0],[0,0,2],[3,0,0]]),.5):assert np.max(abs(v))==0

def test_symbolic_three_identities():
    import sympy as sp
    s,z=sp.symbols('s z');F=sp.Function('F')(s,z);A=-sp.diff(F,z);C=2*F+2*s*sp.diff(F,s)
    assert sp.simplify(2*A+2*s*sp.diff(A,s)+sp.diff(C,z))==0
    norm=sp.pi*sp.integrate(sp.integrate(s/4+z*z,(s,0,4)),(z,-2,2))
    assert sp.simplify(norm-sp.Rational(88,3)*sp.pi)==0
    E,a=sp.symbols('E a');T=2*E
    assert sp.expand(a*T-(T-a*T)/2-E*(3*a-1))==0

def test_training_pool_deterministic():
    a=training_pool();b=training_pool();np.testing.assert_array_equal(a,b)
    assert np.all((a[:,0]>=0)&(a[:,0]<=4))
    assert np.all((a[:,2]>=.25)&(a[:,2]<=.75))


def test_recipe_and_mutation():
    import json,copy
    from replay_st053 import reconstruct
    path=Path(__file__).with_name('recipe.json')
    if not path.exists():pytest.skip('No final candidate frozen yet')
    rec=json.loads(path.read_text());f,r=reconstruct(rec);assert np.isfinite(r).all()
    for key,value in [('pde_validated',True),('parent_id','other'),('modifiers_sha256','0'*64)]:
        bad=copy.deepcopy(rec);bad[key]=value
        with pytest.raises(ValueError):reconstruct(bad)
