"""Independent exact symbolic sector checks for first-collar forcing jets."""
from __future__ import annotations


def run(leaf,c):
    import sympy as S
    from fractions import Fraction
    r,z=S.symbols('r z',real=True);r0=S.Rational(40401,10000);x0=S.Rational(201,100)
    z0=S.Rational(1,4);delta=S.Rational(1,5);lam=S.Integer(3);j=S.Integer(2)
    d=1-z*z;L=1-delta*z*z;H=z*(1-delta)/2+d*(4*z+j)
    ell=-lam*L*H/(H*H+(j/500)**2)
    lder=[S.diff(ell,z,k).subs(z,z0) for k in range(5)];rel=[S.Integer(1)]
    for n in range(5):rel.append(sum(S.binomial(n,k)*lder[k]*rel[n-k] for k in range(n+1)))
    F0=sum(rel[k]*(z-z0)**k/S.factorial(k) for k in range(6))
    Phi=1+r*z+r*r*z**3/10
    mean=1+r*z**2/2+r*r*z**3/12
    U=mean+r*S.diff(mean,r);F=F0*Phi
    Q=(2*z*U-(1-delta)*z*mean-d*S.diff(mean,z))/L
    qr=S.diff(Q,r);qrr=S.diff(Q,r,2)
    force=(Q+(1-delta)*z*S.diff(Q,z)/2+r*qr)/L+Q*(Q/2+r*qr)
    force+=U*(-2*z*Q+d*S.diff(Q,z)-2*z*r*qr)/L
    P=2*qr+r*qrr
    def axial(f,w):return (w*z*f+d*S.diff(f,z)-2*z*r*S.diff(f,r))/L
    nt=-axial(axial(F,-2-delta),-3)
    nz=-axial(axial(U,-1-delta),-2)
    expected={'4':{(0,0):2*x0*P,(0,-1):-x0*force},
      '5':{(1,-1):2*nt},'6':{(0,-1):2*nz-4*r*z*P/L,(0,-2):2*r*z*force/L}}
    def cv(q):q=S.Rational(q);return c.mpf(int(q.p))/int(q.q)
    def grid(expr):
        return {(i,k):cv(S.diff(expr,r,i,z,k).subs({r:r0,z:z0}))
                for i in range(3) for k in range(6-i)}
    result=leaf.build_first_collar_forcing_jets(c,cv(r0),cv(z0),cv(delta),cv(lam),cv(j),
      Phi_derivatives=grid(Phi),Uz_derivatives=grid(U),radial_average_Uz_derivatives=grid(mean),
      log_F0=c.mpf(0),x=cv(x0))
    def fraction(t):
        sign,m,e,_=t;v=Fraction(int(m))*(2**e if e>=0 else Fraction(1,2**(-e)))
        return -v if sign else v
    checks={}
    for row,sectors in expected.items():
        for order in range(4):
            terms=result['g'+row+'_Z'+str(order)];by={(t['F0_power'],t['Lambda_power']):t['coefficient'] for t in terms}
            checks['sector_shape_'+row+'_'+str(order)]=set(by)==set(sectors) and len(terms)==len(sectors)
            for powers,expr in sectors.items():
                exact=S.cancel(S.diff(expr.subs(r,r0),z,order).subs(z,z0));assert exact.is_Rational
                coefficient=by[powers];target=Fraction(int(exact.p),int(exact.q))
                checks['value_'+row+'_'+str(order)+'_'+str(powers)]=fraction(coefficient._mpi_[0])<=target<=fraction(coefficient._mpi_[1])
    missing=grid(Phi);del missing[(0,5)]
    try:
        leaf.build_first_collar_forcing_jets(c,cv(r0),cv(z0),cv(delta),cv(lam),cv(j),
          Phi_derivatives=missing,Uz_derivatives=grid(U),radial_average_Uz_derivatives=grid(mean),log_F0=c.mpf(0),x=cv(x0))
    except (KeyError,ValueError,TypeError):checks['reject_missing_Z5']=True
    else:checks['reject_missing_Z5']=False
    try:
        leaf.build_first_collar_forcing_jets(c,cv(r0),cv(z0),cv(delta),cv(lam),cv(j),
          Phi_derivatives=grid(Phi),Uz_derivatives=grid(U),radial_average_Uz_derivatives=grid(mean),log_F0=c.mpf(0),x=c.mpf(3))
    except ValueError:checks['reject_wrong_x']=True
    else:checks['reject_wrong_x']=False
    checks['source_not_certified']=result['proof']['same_source_identity_proved'] is False
    return {'passed':all(checks.values()),'check_count':len(checks),'checks':checks,
      'synthetic_only':True,'reference':'independent exact rational symbolic sector differentiation'}
