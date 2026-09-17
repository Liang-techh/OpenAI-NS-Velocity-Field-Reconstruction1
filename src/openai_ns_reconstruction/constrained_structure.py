"""Exact symbolic identities for the candidate ansatz (not an NS proof)."""
import json
from pathlib import Path
import sympy as s


def run():
    x,y,z=s.symbols('x y z',real=True)
    rho=x*x+y*y
    H=s.Function('H')(rho,z)
    S=s.Function('S')(rho,z)
    A=s.Matrix([-y*H,x*H,0])
    curl=s.Matrix([s.diff(A[2],y)-s.diff(A[1],z),
                   s.diff(A[0],z)-s.diff(A[2],x),
                   s.diff(A[1],x)-s.diff(A[0],y)])
    swirl=s.Matrix([-y*S,x*S,0])
    div=lambda v:s.simplify(sum(s.diff(v[i],vbl) for i,vbl in enumerate((x,y,z))))
    poloidal=div(curl);azimuthal=div(swirl)
    assert poloidal==0 and azimuthal==0
    p=s.Function('p')(rho,z)
    # r*(e_theta dot grad p): division by r is unnecessary for identity.
    pressure=s.simplify(-y*s.diff(p,x)+x*s.diff(p,y))
    assert pressure==0
    # General axisymmetric radial/swirl vector is equivariant under rotations.
    angle=s.symbols('angle',real=True)
    rotation=s.Matrix([[s.cos(angle),-s.sin(angle),0],
                       [s.sin(angle),s.cos(angle),0],[0,0,1]])
    radial,spin,axial=s.symbols('radial spin axial')
    def vector(a,b):return s.Matrix([a*radial-b*spin,b*radial+a*spin,axial])
    xy=rotation*s.Matrix([x,y,z])
    norm_identity=s.trigsimp(xy[0]**2+xy[1]**2-rho)
    equivariance=(vector(xy[0],xy[1])-rotation*vector(x,y)).applyfunc(s.trigsimp)
    assert norm_identity==0 and equivariance==s.zeros(3,1)
    tau=s.symbols('tau',positive=True)
    h=s.Rational(1,200);a=s.Rational(1,2)+h;d=s.Rational(1,2)-h
    # At fixed R,Z, radial prefactor=lr/lz*tau^-A; other components=tau^-A.
    scale=s.simplify(tau**s.Rational(1,2)/tau**d*tau**(-a)*tau**s.Rational(1,2)-1)
    assert scale==0
    result={'status':'five_symbolic_identities_verified','engine':'SymPy '+s.__version__,
      'identities':[
        {'id':'curl_divergence','result':str(poloidal),'assumptions':'H(rho,z) C2; commuting mixed partials'},
        {'id':'axisymmetric_swirl_divergence','result':str(azimuthal),'assumptions':'S(rho,z) C1'},
        {'id':'axisymmetric_pressure_no_torque','result':str(pressure),'assumptions':'p(rho,z) C1; azimuthal interpretation r>0'},
        {'id':'rotation_equivariance','result':str(equivariance),'assumptions':'radial, spin and axial coefficients depend only on rho,z,t'},
        {'id':'radial_similarity_exponent','result':str(scale),'assumptions':'tau>0; A=101/200, D=99/200; fixed scaled position and time-independent parameters'}],
      'limitations':'Ansatz identities only. Smoothness of bump extension is assumed, not proved by SymPy here. No momentum equation, interval error bound, singularity theorem or floating-point implementation proof.'}
    out=Path('artifacts/constrained/structure_identities.json')
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])

if __name__=='__main__':run()
