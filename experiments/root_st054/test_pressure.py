"""Calibration and scope tests; a pass is not NS numerical acceptance."""
from pathlib import Path
import numpy as np
import pytest
from pressure_completion import ROOT,Family,PressureFit,pressure_basis,residual
PARENT=ROOT/'artifacts/research/ST052-M/candidate.json'

@pytest.fixture(scope='module')
def family():
    from replay_st054 import parent_field
    return parent_field('ST052-M')

def test_pressure_matrix_is_finite_difference_gradient(family):
    f,b=family; rng=np.random.default_rng(9175451);s=rng.uniform(.005,3.2,29);z=rng.uniform(-1.8,1.8,29);t=rng.uniform(.25,.75,29)
    p0,p1=pressure_basis(f,s,z,t);dq=rng.normal(0,1e-4,f.n);c=b.copy();c[2*f.n:3*f.n]+=dq
    X=np.c_[np.sqrt(s),np.zeros(len(s)),z];h=1e-5
    diff=lambda x:f.fields(c,x,t)[1]-f.fields(b,x,t)[1]
    for k,truth in [(0,p0@dq),(2,p1@dq)]:
        e=np.eye(3)[k]*h;actual=(diff(X+e)-diff(X-e))/(2*h)
        np.testing.assert_allclose(actual,truth,atol=2e-8,rtol=1e-5)

def test_full_residual_affine_and_azimuthal_invariant(family):
    f,b=family;rng=np.random.default_rng(9175452);p=rng.uniform(-1.6,1.6,(37,3));t=rng.uniform(.25,.75,37);c=b.copy();dq=rng.normal(0,1e-4,f.n);c[2*f.n:3*f.n]+=dq
    rb=residual(f,b,p,t);rc=residual(f,c,p,t);rad=np.linalg.norm(p[:,:2],axis=1);pr,pz=pressure_basis(f,rad**2,p[:,2],t)
    expected=np.c_[p[:,0]/rad*(pr@dq),p[:,1]/rad*(pr@dq),pz@dq]
    np.testing.assert_allclose(rc-rb,expected,atol=5e-13,rtol=1e-10)
    theta=np.c_[-p[:,1]/rad,p[:,0]/rad,np.zeros(len(rad))]
    assert np.max(abs(np.sum(theta*(rc-rb),axis=1)))<1e-14
    np.testing.assert_array_equal(f.fields(b,p,t)[0],f.fields(c,p,t)[0])

def test_independent_cartesian_residual(family):
    from validate import cartesian_residual
    from spacetime import force
    f,b=family;c=b.copy();rng=np.random.default_rng(9175453);c[2*f.n:3*f.n]+=rng.normal(0,1e-4,f.n);x=np.array([[.2,.1,.1],[.35,-.3,.5],[.5,.3,-.7]])
    fd,_=cartesian_residual(lambda x,t:f.fields(c,x,t),lambda x,t:force(x,t,*c[-2:]),x,.5,.00125,.000625)
    np.testing.assert_allclose(fd,residual(f,c,x,.5),atol=2e-7,rtol=1e-5)

def test_pressure_gram_stationarity_and_quadratic(family,tmp_path):
    f,r=family;p=tmp_path/'parent.json';f.save(r,p);o=PressureFit(p,(18,26),11);n=o.f.n
    assert np.max(abs(o.G@o.unconstrained+o.g))<1e-12
    assert np.linalg.eigvalsh(o.G)[0]>0
    d=.2*o.unconstrained
    expected=o.cost+2*o.g@d+d@o.G@d
    assert abs(o.metrics(d)['weighted_mse']-expected)<1e-13
    assert o.metrics(o.unconstrained)['weighted_mse']<=o.cost

def test_symbolic_pressure_invariances():
    import sympy as S
    x,y,z,t=S.symbols('x y z t');s=S.symbols('s');p=S.Function('p')(x,y,z,t)
    assert S.simplify(S.diff(p,x,y)-S.diff(p,y,x))==0
    P=S.Function('P')(s,z,t);F=P.subs(s,x*x+y*y)
    assert S.simplify(-y*S.diff(F,x)+x*S.diff(F,y))==0
    phi=S.Matrix([x/2,y/2,-z]);assert sum(S.diff(phi[i],v) for i,v in enumerate([x,y,z]))==0
    psi=S.Function('psi')(s,z);A=-S.diff(psi,z);C=2*psi+2*s*S.diff(psi,s)
    assert S.simplify(2*A+2*s*S.diff(A,s)+S.diff(C,z))==0

def test_final_velocity_force_and_support(family):
    from spacetime import quad
    for label,par in [('M3','ST052-M'),('Q2','ST053-Q')]:
        from replay_st054 import parent_field,reconstruct
        f,b=parent_field(par);_,c=reconstruct('ST054-'+label)
        np.testing.assert_array_equal(c[:2*f.n],b[:2*f.n]);np.testing.assert_array_equal(c[-2:],b[-2:])
        for v in f.fields(c,np.array([[2.,0.,0.],[0.,0.,2.],[2.1,0,0],[0,0,-2.1]]),.5):assert np.max(abs(v))==0
        s,z,w=quad(96);u=f.fields(c,np.c_[np.sqrt(s),np.zeros(len(s)),z],.25)[0]
        assert abs(.5*w@np.sum(u*u,axis=1)-1)<1e-6


def test_replay_claim_hash_and_parent_mutation():
    import json,copy
    from replay_st054 import HERE,reconstruct
    records=json.loads((HERE/'recipes.json').read_text())
    for ident in records:
        reconstruct(ident,records)
        for key,value,match in [('pde_validated',True,'scientific claim'),('delta_sha256','0'*64,'checksum'),('parent_id','fake','parent identity')]:
            bad=copy.deepcopy(records);bad[ident][key]=value
            with pytest.raises(ValueError,match=match):reconstruct(ident,bad)


def test_multidimensional_midplane_audit_component_axes(tmp_path):
    from replay_st054 import parent_field,reconstruct
    from audit_pressure import audit
    f,b=parent_field('ST052-M');_,c=reconstruct('ST054-M3')
    p=tmp_path/'parent.json';q=tmp_path/'child.json';f.save(b,p);f.save(c,q)
    report=audit(p,q,tmp_path/'audit.json')
    assert report['shear_probe_count']==1025
    assert report['shear_max_difference']==0.
    assert report['midplane_velocity_max_difference']==0.
    assert report['velocity_coefficients_identical'] is True
