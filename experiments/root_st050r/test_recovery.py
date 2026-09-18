"""Independent derivative/identity checks; success is not PDE acceptance."""
from pathlib import Path
import sys
import numpy as np
import pytest
import pressure_morph as pm
from spacetime import Family
from validate import cartesian_residual

ROOT=Path(__file__).resolve().parents[2]
PARENT=ROOT/'artifacts/research/ST048-S/candidate.json'

@pytest.fixture(scope='module')
def model(tmp_path_factory):
    from replay_recovery import parent_field
    f,r=parent_field();p=tmp_path_factory.mktemp('parent')/'candidate.json';f.save(r,p)
    return pm.LocalizedModel(p,2,2,3)

def test_poisson_parameter_derivative(model):
    rng=np.random.default_rng(9175051);s=rng.uniform(.005,3.,41);z=rng.uniform(-1.5,1.5,41);t=rng.uniform(.25,.75,41)
    D=model.cache(s,z,t);lp=pm.pressure_laplacian_cache(model,D);c=np.zeros(model.dim)
    c[:-2]=rng.normal(0,.0001,model.dim-2);direction=rng.normal(size=model.dim);direction/=np.linalg.norm(direction)
    value,J=pm.poisson_residual(model,c,D,lp);h=1e-6
    a=pm.poisson_residual(model,c+h*direction,D,lp)[0];b=pm.poisson_residual(model,c-h*direction,D,lp)[0]
    err=np.linalg.norm((a-b)/(2*h)-J@direction)/np.linalg.norm(J@direction)
    assert err<1e-7,err

def test_poisson_is_cartesian_divergence_of_full_residual(model):
    x=np.array([[.15,.2,.1],[0,0,.13],[.4,.2,-.32],[.8,.1,.65]])
    c=np.zeros(model.dim);D=model.cache(x[:,0]**2+x[:,1]**2,x[:,2],.5);lp=pm.pressure_laplacian_cache(model,D)
    truth=pm.poisson_residual(model,c,D,lp)[0];raw=model.candidate(c);errors=[]
    for h in [.02,.01,.005]:
        div=np.zeros(len(x))
        for k in range(3):
            d=np.eye(3)[k]*h
            fun=lambda pts:model.f.analytic_residual(raw,pts,.5)[:,k]
            div+=(fun(x-2*d)-8*fun(x-d)+8*fun(x+d)-fun(x+2*d))/(12*h)
        errors.append(float(np.max(abs(div-truth))))
    assert errors[-1]<1e-6 and errors[0]/errors[1]>10,errors

def test_independent_cartesian_operator_divergence(model):
    x=np.array([[.2,.1,.17],[0.,0.,-.12]])
    c=np.zeros(model.dim);raw=model.candidate(c)
    from spacetime import force
    field=lambda x,t:model.f.fields(raw,x,t)
    forcing=lambda x,t:force(x,t,*raw[-2:])
    D=model.cache(x[:,0]**2+x[:,1]**2,x[:,2],.5);truth=pm.poisson_residual(model,c,D,pm.pressure_laplacian_cache(model,D))[0]
    errors=[]
    for h in [.01,.005]:
        div=np.zeros(len(x))
        for k in range(3):
            e=np.eye(3)[k]*h
            def F(p):return cartesian_residual(field,forcing,p,.5,h,.001)[0][:,k]
            div+=(F(x-2*e)-8*F(x-e)+8*F(x+e)-F(x+2*e))/(12*h)
        errors.append(float(np.max(abs(div-truth))))
    assert errors[-1]<2e-6,errors

def test_vorticity_gram_derivative(model):
    G=pm.vorticity_grams(model,order=(12,18),times=(.25,.5))
    rng=np.random.default_rng(9175052);c=rng.normal(0,.0001,model.nv);d=np.r_[1,c];v=rng.normal(size=model.nv);v/=np.linalg.norm(v)
    h=1e-6
    fun=lambda cc:np.einsum('i,tkij,j->tk',np.r_[1,cc],G,np.r_[1,cc])
    grad=2*np.einsum('tkij,j->tki',G,d)[:,:,1:]
    np.testing.assert_allclose((fun(c+h*v)-fun(c-h*v))/(2*h),grad@v,atol=1e-7,rtol=1e-7)

def test_symbolic_axis_identity_and_parity():
    import sympy as sp
    A,B,As,Bs,Az,Cs,Cz,s=sp.symbols('A B As Bs Az Cs Cz s',real=True)
    r=sp.symbols('r',real=True)
    J=sp.Matrix([[A+2*r*r*As,-B,r*Az],[B+2*r*r*Bs,A,r*B],[2*r*Cs,0,Cz]])
    tr=sp.expand(sp.trace(J*J));expected=(A+2*r*r*As)**2+A*A+Cz*Cz-2*B*(B+2*r*r*Bs)+4*r*r*Az*Cs
    assert sp.expand(tr-expected)==0
    assert sp.simplify(tr.subs({r:0,Cz:-2*A})-(6*A*A-2*B*B))==0
    s,z=sp.symbols('s z');F=sp.Function('F')(s,z)
    a=-sp.diff(F,z);cc=2*F+2*s*sp.diff(F,s)
    assert sp.simplify(2*a+2*s*sp.diff(a,s)+sp.diff(cc,z))==0

def test_original_energy_and_support(model):
    from spacetime import quad
    s,z,w=quad(96);raw=model.candidate(np.zeros(model.dim));x=np.c_[np.sqrt(s),np.zeros(len(s)),z]
    u,_=model.f.fields(raw,x,.25)
    assert abs(.5*w@np.sum(u*u,axis=1)-1)<1e-6
    x=np.array([[2.,0,0],[0,0,2.],[3.,0,0],[0,0,-2.]])
    for a in model.f.fields(raw,x,.5):assert np.max(abs(a))==0


def test_reconstruction_and_claim_rejection():
    import copy,json
    from replay_recovery import reconstruct
    rec=json.loads((Path(__file__).with_name('recipes.json')).read_text())
    for key in rec:
        f,r=reconstruct(key,rec);assert np.isfinite(r).all()
        bad=copy.deepcopy(rec);bad[key]['pde_validated']=True
        with pytest.raises(ValueError,match='unsupported claim'):reconstruct(key,bad)
        bad=copy.deepcopy(rec);bad[key]['modifiers_sha256']='0'*64
        with pytest.raises(ValueError,match='checksum'):reconstruct(key,bad)
