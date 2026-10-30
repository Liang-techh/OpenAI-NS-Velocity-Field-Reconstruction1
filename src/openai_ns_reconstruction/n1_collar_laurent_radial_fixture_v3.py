"""Independent symbolic tests of collar radial sector identities."""
from __future__ import annotations


def run(leaf,c):
    import sympy as S
    from fractions import Fraction
    r,z=S.symbols('r z',real=True);r0=S.Rational(40401,10000);z0=S.Rational(1,5)
    hb=S.Rational(1,3);lam=S.Integer(5);sigma=S.Rational(1,4);sp=S.Rational(1,7)
    core_phi=2+r+z+r*r*z*z/50;core_psi=1+r*z+r*r*z*z/30
    phi=3+z+z*z;U=2+z*z+z**3;mean=1+z**3
    A=S.diff(core_phi,r)/core_phi;Ar=S.diff(A,r);chi=1+(hb-1)*sigma;chir=(hb-1)*sp/(hb*r)
    ratio=phi/core_phi
    pr=phi*chi*A;prr=phi*(chir*A+chi**2*A**2+chi*Ar)
    ur=ratio*chi*S.diff(core_psi,r)/lam
    urr=ratio*(chir*S.diff(core_psi,r)+chi*(chi-1)*A*S.diff(core_psi,r)+chi*S.diff(core_psi,r,2))/lam
    mr=(U-mean)/r;mrr=(ur-2*mr)/r
    expected={'Phi':(phi,pr,prr),'U':(U,ur,urr),'mean':(mean,mr,mrr)}
    def cv(q):q=S.Rational(q);return c.mpf(int(q.p))/int(q.q)
    def grid(f):return {(i,k):cv(S.diff(f,r,i,z,k).subs({r:r0,z:z0})) for i in range(3) for k in range(6-i)}
    def actual(f):return {(0,k):[{'powers':(0,0,0,0),'coefficient':cv(S.diff(f,z,k).subs(z,z0)),'origin':('actual',k)}] for k in range(6)}
    cert={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True,'sigma_Z_independent':True}
    atoms={k:'synthetic_'+k for k in ('hb','F0_anchor','Pstar_squared','Lambda')}
    out=leaf.build_collar_laurent_radial(c,cv(r0),cv(z0),cv(hb),actual_Phi=actual(phi),actual_U=actual(U),own_mean=actual(mean),
      comparison_Phi=grid(core_phi),comparison_Psi=grid(core_psi),sigma=[cv(sigma)]+[c.mpf(0)]*5,
      sigma_prime=[cv(sp)]+[c.mpf(0)]*5,source_family='synthetic',source_atoms=atoms,geometry_certificate=cert)
    bases=(cv(hb),c.mpf(7)/3,c.mpf(3)/2,cv(lam));checks={}
    def fraction(bits):
        sign,m,e,_=bits;v=Fraction(int(m))*(2**e if e>=0 else Fraction(1,2**(-e)))
        return -v if sign else v
    for name,rows in expected.items():
        for i,expr in enumerate(rows):
            for k in range(6-i):
                exact=S.cancel(S.diff(expr,z,k).subs({r:r0,z:z0}));assert exact.is_Rational
                total=c.mpf(0)
                for term in out[name][(i,k)]:
                    v=term['coefficient']
                    for base,p in zip(bases,term['powers']):v*=base**p
                    total+=v
                t=Fraction(int(exact.p),int(exact.q));checks[name+'_'+str(i)+'_'+str(k)]=fraction(total._mpi_[0])<=t<=fraction(total._mpi_[1])
    checks['negative_hb_phi2']=any(t['powers'][0]<0 for t in out['Phi'][(2,0)])
    checks['negative_hb_U2']=any(t['powers'][0]<0 for t in out['U'][(2,0)])
    checks['epsilon_factor_Ur']=all(t['powers'][3]==-1 for t in out['U'][(1,0)])
    checks['actual_source_unproved']=out['metadata']['actual_source_identity_proved'] is False
    return {'passed':all(checks.values()),'check_count':len(checks),'checks':checks,'synthetic_only':True,
      'reference':'independent exact symbolic collar radial identities including epsilon and hb^-1'}
