"""Independent symbolic test of negative-width source factors in forcing."""
from __future__ import annotations


def run(leaf,c):
    import sympy as S
    from fractions import Fraction
    r,z=S.symbols('r z',real=True);r0=S.Rational(40401,10000);x0=S.Rational(201,100)
    z0=S.Rational(1,5);delta=S.Integer(0);lam=S.Integer(3);j=S.Integer(2);hb=S.Rational(1,3);fanchor=S.Rational(7,3)
    d=1-z*z;L=1-delta*z*z;H=z*(1-delta)/2+d*(4*z+j)
    ell=-lam*L*H/(H*H+(j/500)**2)
    lder=[S.diff(ell,z,k).subs(z,z0) for k in range(5)];rel=[S.Integer(1)]
    for n in range(5):rel.append(sum(S.binomial(n,k)*lder[k]*rel[n-k] for k in range(n+1)))
    F0=fanchor*sum(rel[k]*(z-z0)**k/S.factorial(k) for k in range(6))
    p0=1+r*z;p1=r*r*z*z/100
    m0=1+r*z/3;m1=r*r*z*z/120
    u0=m0+r*S.diff(m0,r);u1=m1+r*S.diff(m1,r)
    phi=p0+p1/hb;mean=m0+m1/hb;U=u0+u1/hb;F=F0*phi
    Q=(2*z*U-(1-delta)*z*mean-d*S.diff(mean,z))/L;qr=S.diff(Q,r);qrr=S.diff(Q,r,2)
    P=2*qr+r*qrr
    Sflow=(Q+(1-delta)*z*S.diff(Q,z)/2+r*qr)/L+Q*(Q/2+r*qr)
    Sflow+=U*(-2*z*Q+d*S.diff(Q,z)-2*z*r*qr)/L
    def axial(f,w):return (w*z*f+d*S.diff(f,z)-2*z*r*S.diff(f,r))/L
    Nt=-axial(axial(F,-2-delta),-3);Nz=-axial(axial(U,-1-delta),-2)
    expected={'4':2*x0*(lam*P-Sflow/2)/lam,'5':2*Nt/lam,
      '6':2*Nz/lam-4*z*r*(lam*P-Sflow/2)/(lam**2*L)}
    def cv(q):q=S.Rational(q);return c.mpf(int(q.p))/int(q.q)
    def grid(first,second):
        out={}
        for i in range(3):
            for k in range(6-i):
                a=S.diff(first,r,i,z,k).subs({r:r0,z:z0});b=S.diff(second,r,i,z,k).subs({r:r0,z:z0})
                out[(i,k)]=[{'powers':(0,0,0,0),'coefficient':cv(a),'origin':('part0',i,k)},
                  {'powers':(-1,0,0,0),'coefficient':cv(b),'origin':('hb_negative',i,k)}]
        return out
    args=dict(Phi_derivatives=grid(p0,p1),Uz_derivatives=grid(u0,u1),radial_average_Uz_derivatives=grid(m0,m1),log_F0=c.mpf(0),x=cv(x0))
    out=leaf.build_first_collar_laurent_forcing_jets(c,cv(r0),cv(z0),cv(delta),cv(lam),cv(j),**args)
    bases=[cv(hb),cv(fanchor),c.mpf(3)/2,cv(lam)];checks={}
    def exactfrac(bits):
        sign,m,e,_=bits;v=Fraction(int(m))*(2**e if e>=0 else Fraction(1,2**(-e)))
        return -v if sign else v
    for k in range(4):
        for row,expr in expected.items():
            terms=out['sectors'][str(k)][row];total=c.mpf(0)
            for term in terms:
                val=term['coefficient']
                for base,p in zip(bases,term['powers']):val*=base**p
                total+=val
            target=S.cancel(S.diff(expr.subs(r,r0),z,k).subs(z,z0));assert target.is_Rational
            exact=Fraction(int(target.p),int(target.q))
            checks['g'+row+'_Z'+str(k)]=exactfrac(total._mpi_[0])<=exact<=exactfrac(total._mpi_[1])
            checks['provenance_'+row+'_'+str(k)]=all('origin' in term for term in terms)
    checks['negative_hb_preserved']=any(t['powers'][0]<0 for k in range(4) for t in out['sectors'][str(k)]['5'])
    checks['quadratic_negative_hb_preserved']=any(t['powers'][0]<=-2 for k in range(4) for t in out['sectors'][str(k)]['4'])
    checks['source_not_certified']=out['proof'].get('same_source_identity_proved',False) is False
    return {'passed':all(checks.values()),'check_count':len(checks),'checks':checks,'synthetic_only':True,
      'reference':'independent symbolic physical forcing with explicit nonzero hb^-1 source fields',
      'sector_counts':{k:{row:len(terms) for row,terms in v.items()} for k,v in out['sectors'].items()}}
