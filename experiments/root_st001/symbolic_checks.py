"""Small symbolic identities only; no claim of a complete NS/Lean proof."""
import json
from pathlib import Path
import sympy as S

def run(out=None):
    x,y,z,s,r,a,c=S.symbols('x y z s r a c',real=True)
    F=S.Function('F')(s,z);B=S.Function('B')(s,z)
    A=-S.diff(F,z);C=2*F+2*s*S.diff(F,s)
    # Axisymmetric divergence, written using s=r^2 (regular at r=0).
    div=S.simplify(2*A+2*s*S.diff(A,s)+S.diff(C,z))
    assert div==0
    H=(x*x+y*y)/4-z*z/2
    phi=S.Matrix([S.diff(H,q) for q in (x,y,z)])
    assert sum(S.diff(H,q,2) for q in (x,y,z))==0
    ux,uy,uz=S.symbols('ux uy uz')
    U=S.Matrix([ux,uy,uz])
    contraction=S.expand(-(U.T*S.hessian(H,(x,y,z))*U)[0])
    assert S.simplify(contraction-(uz**2-(ux**2+uy**2)/2))==0
    phi_norm_sq=S.integrate(S.integrate(2*S.pi*r*(r*r/4+z*z),(r,0,2)),(z,-2,2))
    assert S.simplify(phi_norm_sq-88*S.pi/3)==0
    torque_norm_sq=S.integrate(S.integrate(2*S.pi*r**3,(r,0,2)),(z,-2,2))
    assert S.simplify(torque_norm_sq-32*S.pi)==0
    b=S.Function('b')(r)
    # Integral of r*f_theta times cylindrical Jacobian reduces to -c*int r^3 b.
    torque_integrand=r**3*b+r**4*S.diff(b,r)/2
    assert S.simplify(torque_integrand-(S.diff(r**4*b,r)/2-r**3*b))==0
    for j in range(7):
        assert S.expand(S.legendre(2*j+1,-z)+S.legendre(2*j+1,z))==0
        assert S.expand(S.legendre(2*j,-z)-S.legendre(2*j,z))==0
    energy=S.symbols('energy',positive=True)
    assert S.simplify((S.sqrt(2/energy))**2*energy/2)==1
    result={'passed':5,'checks':['divergence identity for arbitrary differentiable F,B','harmonic pressure-independent anisotropy contraction and test-field L2 norm','restricted-force angular torque integration-by-parts identity','odd/even axial polynomial parity for all used orders','initial energy normalization algebra'],'phi_L2_squared':str(phi_norm_sq),'torque_test_L2_squared':str(torque_norm_sq),'scope':'Symbolic local/structural identities, with compact-support integration-by-parts hypotheses stated in README. No NS existence theorem, no continuum residual enclosure, and no Lean build.'}
    if out:Path(out).write_text(json.dumps(result,indent=2)+'\n')
    return result
if __name__=='__main__':
    print(json.dumps(run(Path(__file__).with_name('symbolic_results.json')),indent=2))
