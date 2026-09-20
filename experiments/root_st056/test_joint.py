"""Focused structure/operator tests: not numerical NS acceptance."""
from pathlib import Path
import numpy as np
import pytest
import joint_step as js
PARENT=js.ROOT/'artifacts/research/ST054-Q2/candidate.json'
@pytest.fixture(scope='module')
def model(tmp_path_factory):
    from replay_st056 import parent_field
    f,r=parent_field();p=tmp_path_factory.mktemp('parent')/'candidate.json';f.save(r,p)
    return js.JointModel(p)

def test_exact_original_evaluator_parity(model):
    rng=np.random.default_rng(9205652);p=np.c_[rng.uniform(0,3.5,27),rng.uniform(-1.8,1.8,27),rng.uniform(.25,.75,27)]
    D=model.cache(p);x=np.c_[np.sqrt(p[:,0]),np.zeros(len(p)),p[:,1]]
    exact=model.f.analytic_residual(model.raw,x,p[:,2])
    np.testing.assert_allclose(D['R'],exact,atol=2e-11,rtol=2e-10)
    np.testing.assert_allclose(D['u'],model.f.fields(model.raw,x,p[:,2])[0],atol=1e-12)

def test_exact_quadratic_updated_field(model):
    rng=np.random.default_rng(9205653);c=rng.normal(0,1e-6,model.dim);p=np.c_[rng.uniform(.001,3.,23),rng.uniform(-1.7,1.7,23),rng.uniform(.25,.75,23)]
    D=model.cache(p);raw=model.candidate(c);x=np.c_[np.sqrt(p[:,0]),np.zeros(len(p)),p[:,1]]
    R,_=model.residual(c,D)
    np.testing.assert_allclose(R,model.f.analytic_residual(raw,x,p[:,2]),atol=1e-11,rtol=1e-8)

def test_joint_pullback(model):
    rng=np.random.default_rng(9205654);c=rng.normal(0,1e-6,model.dim);d=rng.normal(size=model.dim);d/=np.linalg.norm(d)
    p=np.c_[rng.uniform(.001,3,19),rng.uniform(-1.6,1.6,19),rng.uniform(.25,.75,19)];D=model.cache(p);w=rng.normal(size=(len(p),3));h=1e-6
    fun=lambda v:float(np.sum(w*model.residual(v,D)[0]))
    np.testing.assert_allclose((fun(c+h*d)-fun(c-h*d))/(2*h),model.pullback(c,D,w)@d,atol=5e-8,rtol=1e-6)

def test_initial_velocity_poloidal_and_force_preservation(model):
    rng=np.random.default_rng(9205655);c=rng.normal(0,1e-6,model.dim);raw=model.candidate(c);f=model.f
    x=rng.uniform(-1.5,1.5,(43,3));u0=f.fields(model.raw,x,.25)[0];u=f.fields(raw,x,.25)[0]
    np.testing.assert_allclose(u,u0,atol=5e-13,rtol=5e-13)
    np.testing.assert_array_equal(raw[:f.n],model.raw[:f.n]);np.testing.assert_array_equal(raw[-2:],model.raw[-2:])
    for t in [.25,.517,.75]:
        u=f.fields(raw,x,t)[0];v=f.fields(model.raw,x,t)[0]
        np.testing.assert_allclose(u[:,2],v[:,2],atol=2e-13)
        np.testing.assert_allclose(x[:,0]*u[:,0]+x[:,1]*u[:,1],x[:,0]*v[:,0]+x[:,1]*v[:,1],atol=2e-13)

def test_independent_cartesian_operator(model):
    from validate import cartesian_residual
    rng=np.random.default_rng(9205656);c=rng.normal(0,1e-6,model.dim);raw=model.candidate(c)
    x=np.array([[0.,0.,.12],[.11,.04,.17],[.37,.13,-.21],[.27,0,1.8]])
    fun=lambda x,t:model.f.fields(raw,x,t);ff=lambda x,t:js.force(x,t,*raw[-2:])
    R,div=cartesian_residual(fun,ff,x,.537,.00125,.000625)
    np.testing.assert_allclose(R,model.f.analytic_residual(raw,x,.537),atol=2e-7,rtol=2e-5)
    assert abs(div).max()<2e-7

def test_energy_support(model):
    from spacetime import quad
    rng=np.random.default_rng(9205657);raw=model.candidate(rng.normal(0,1e-6,model.dim))
    s,z,w=quad(96);u=model.f.fields(raw,np.c_[np.sqrt(s),np.zeros(len(s)),z],.25)[0]
    assert abs(.5*w@np.sum(u*u,axis=1)-1)<1e-6
    for t in [.25,.5,.75]:
        for a in model.f.fields(raw,np.array([[2,0,0],[0,0,2],[0,0,-2],[3,0,0]]),t):assert np.count_nonzero(a)==0

def test_four_symbolic_identities():
    import sympy as sp
    s,z,r,B,d,t=sp.symbols('s z r B d t');F=sp.Function('F')(s,z)
    A=-sp.diff(F,z);C=2*F+2*s*sp.diff(F,s)
    assert sp.simplify(2*A+2*s*sp.diff(A,s)+sp.diff(C,z))==0
    assert sp.expand(-r*(B+d)**2+r*B**2+r*(2*B*d+d**2))==0
    th=sp.symbols('theta');p=sp.Function('p')(s,z,t)
    assert sp.diff(p,th)==0
    for k in range(1,8):assert sp.chebyshevt(k,-1)-(-1)**k==0

def test_frozen_replay_and_rejections():
    import copy,json
    from pathlib import Path
    from replay_st056 import reconstruct
    paths=sorted(Path(__file__).parent.glob('ST056-*.json'))
    assert any(p.stem=='ST056-J2' for p in paths)
    for p in paths:
        ident=p.stem;record=json.loads(p.read_text())
        f,r=reconstruct(ident,record);assert np.isfinite(r).all()
        bad=copy.deepcopy(record);bad['pde_validated']=True
        with pytest.raises(ValueError,match='scientific claim'):reconstruct(ident,bad)
        bad=copy.deepcopy(record);bad['delta_npy_sha256']='0'*64
        with pytest.raises(ValueError,match='checksum'):reconstruct(ident,bad)
        bad=copy.deepcopy(record);bad['parent']='invented'
        with pytest.raises(ValueError,match='identity'):reconstruct(ident,bad)

def test_harmonic_weak_moment_constants():
    import sympy as s
    x,y,z,r,th=s.symbols('x y z r th',real=True)
    H=(x*x+y*y)/4-z*z/2
    phi=s.Matrix([s.diff(H,v) for v in (x,y,z)])
    assert s.diff(phi[0],x)+s.diff(phi[1],y)+s.diff(phi[2],z)==0
    a,b,c=s.symbols('a b c',real=True);u=s.Matrix([a,b,c])
    assert s.expand(-(u.T*phi.jacobian([x,y,z])*u)[0]-(c*c-(a*a+b*b)/2))==0
    integral=s.integrate((r*r/4+z*z)*r,(r,0,2),(th,0,2*s.pi),(z,-2,2))
    assert s.simplify(integral-88*s.pi/3)==0
    # This algebra does not interval-certify the numerical velocity integrals.
