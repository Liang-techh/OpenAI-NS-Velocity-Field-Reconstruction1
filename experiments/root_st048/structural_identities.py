"""Small exact identities only, not full NS existence or source correspondence."""
import json
from pathlib import Path
import sympy as S

def check():
    s,z,r=S.symbols('s z r', real=True)
    F=S.Function('F')(s,z)
    A=-S.diff(F,z);C=2*F+2*s*S.diff(F,s)
    divergence=S.simplify(2*A+2*s*S.diff(A,s)+S.diff(C,z))
    assert divergence==0
    shear=S.diff(C.subs(s,r*r),r)
    assert S.simplify(shear-2*r*S.diff(C,s).subs(s,r*r))==0
    E=S.symbols('E',positive=True)
    assert S.simplify(S.sqrt(2/E)**2*E/2)==1
    return dict(checks_passed=3,identities=['axisymmetric Cartesian divergence for arbitrary smooth F','radial axial shear derivative with s=r^2','nonzero initial energy normalization'],scope='Symbolic identities under smoothness and positive normalization. Not a complete NS solution, pressure-sign certificate or continuum residual bound.',pde_validated=False)
if __name__=='__main__':print(json.dumps(check(),indent=2))
