"""Actual derivative/identity/trajectory checks; test success is not PDE acceptance."""
from pathlib import Path
import base64,copy,json
import numpy as np
import pytest
import continuation as cont
from replay_st047 import reconstruct
from structure_and_particles import integrate
from spacetime import quad,force
from validate import cartesian_residual

@pytest.fixture(scope='module')
def objective():
    root=Path(__file__).resolve().parents[2]
    p=root/'artifacts/research/ST046-A/candidate.json'
    if not p.exists():
        # Same checked mathematical parent on source-only GitHub checkout.
        import sys,tempfile
        sys.path.insert(0,str(root/'experiments/root_st046'))
        from replay_st046 import reconstruct as prior
        f,r=prior('ST046-A');tmp=tempfile.TemporaryDirectory();p=Path(tmp.name)/'parent.json';f.save(r,p)
    else:tmp=None
    m=cont.LocalizedModel(p,1,1,2)
    yield cont.Objective(m,space_order=(8,10),profile_ratio=.995,pressure_ratio=.99,time_cap=1.,softmax_weight=.02)
    if tmp:tmp.cleanup()

def test_conditioned_derivatives(objective):
    o=objective;m=o.m;c=np.zeros(m.dim);T,rs,info=cont.precondition(o,c)
    rng=np.random.default_rng(9174793);d=rng.normal(size=m.dim);d/=np.linalg.norm(d);d=T@d
    h=1e-6;fv,g=o.fun(c);v,j=o.constraints(c)
    assert abs((o.fun(h*d)[0]-o.fun(-h*d)[0])/(2*h)-g@d)<1e-7
    assert np.max(abs((o.constraints(h*d)[0]-o.constraints(-h*d)[0])/(2*h)-j@d))<1e-6
    assert np.isfinite(T).all() and np.isfinite(rs).all()

def test_pressure_cancels_nonpressure_acceleration(objective):
    o=objective;m=o.m;z=np.zeros(m.dim);a=o.constraints(z)[0][-len(o.parent_mz):].copy()
    z[m.nv:m.nv+m.np]=np.random.default_rng(17).normal(0,1e-4,m.np)
    b,j=o.constraints(z)
    np.testing.assert_allclose(b[-len(a):],a,atol=2e-15,rtol=0)
    assert np.max(abs(j[-len(a):,m.nv:m.nv+m.np]))<1e-12

def test_reduced_momentum_against_cartesian_fd(objective):
    m=objective.m;g=np.random.default_rng(71);x=g.uniform(-1.3,1.3,(40,3));x[:3]=[[0,0,0],[0,0,.2],[0,0,-.2]]
    c=g.normal(0,1e-5,m.dim);raw=m.candidate(c);t=.4375
    field=lambda p,t:m.f.fields(raw,p,t);forcing=lambda p,t:force(p,t,*raw[-2:])
    actual,div=cartesian_residual(field,forcing,x,t,.00125,.000625)
    direct=m.f.analytic_residual(raw,x,t)
    assert np.max(abs(actual-direct))<1e-6
    assert np.max(abs(div))<1e-6

@pytest.mark.parametrize('ident',['ST047-W','ST047-E'])
def test_frozen_field_energy_support_rotation(ident):
    f,r=reconstruct(ident);s,z,w=quad(96);pts=np.c_[np.sqrt(s),np.zeros(len(s)),z]
    # Chunk only to bound memory; no change to quadrature.
    u=np.concatenate([f.fields(r,pts[k:k+512],.25)[0] for k in range(0,len(pts),512)])
    assert abs(.5*w@np.sum(u*u,axis=1)-1)<1e-6
    b=np.array([[2,0,0],[0,0,2],[0,0,-2],[3,3,3]],float)
    for value in f.fields(r,b,.5):assert np.max(abs(value))==0
    x=np.random.default_rng(72).uniform(-1,1,(40,3));angle=.32;Q=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
    np.testing.assert_allclose(f.fields(r,x@Q.T,.5)[0],f.fields(r,x,.5)[0]@Q.T,atol=1e-11,rtol=1e-11)

@pytest.mark.parametrize('mutation',['pde_validated','parent_id','modifiers_npy_b64'])
def test_recipe_rejects_mutation(mutation):
    records=json.loads(Path(__file__).with_name('recipes.json').read_text())
    records['ST047-E'][mutation]=True if mutation=='pde_validated' else 'ST006' if mutation=='parent_id' else base64.b64encode(b'invalid').decode()
    with pytest.raises((ValueError,TypeError)):reconstruct('ST047-E',records)

def test_nonautonomous_particle_integrator():
    initial=np.array([[.1,.2,.1],[.3,-.1,-.2],[.1,0,0.]])
    def field(p,t):return np.c_[-.2*p[:,0]-.3*p[:,1],.3*p[:,0]-.2*p[:,1],.4*p[:,2]+.002]
    times,path,_=integrate(field,initial,rtol=1e-10,max_step=.01)
    h=.5;Q=np.array([[np.cos(.3*h),-np.sin(.3*h)],[np.sin(.3*h),np.cos(.3*h)]])
    expected=initial.copy();expected[:,:2]=np.exp(-.2*h)*initial[:,:2]@Q.T;expected[:,2]=(initial[:,2]+.005)*np.exp(.4*h)-.005
    np.testing.assert_allclose(path[-1],expected,atol=1e-10,rtol=1e-10)

def test_symbolic_structure_and_pressure_bound():
    import sympy as s
    v,z=s.symbols('v z');F=s.Function('F')(v,z);A=-s.diff(F,z);C=2*F+2*v*s.diff(F,v)
    assert s.simplify(2*A+2*v*s.diff(A,v)+s.diff(C,z))==0
    N=s.symbols('N',positive=True)
    assert s.simplify(s.sqrt(2/N)**2*N/2-1)==0
    material,pressure=s.symbols('material pressure',real=True)
    assert s.diff((material+pressure)-pressure,pressure)==0
    aligned_pressure=s.symbols('aligned_pressure',nonnegative=True)
    assert ((material+aligned_pressure)-material).is_nonnegative
    for M in [-2.,-.1,0.,.1,2.]:
        for pg in [0.,.01,.5,4.]:assert abs(M+pg)>=max(M,0)
