"""Independent synthetic own-mean history and radial-identity checks."""
from __future__ import annotations


def run(leaf,c):
    import mpmath as mp
    from fractions import Fraction
    import math
    z0=Fraction(1,5)
    def jet(coeff,k):
        return sum(Fraction(a*math.factorial(n),math.factorial(n-k))*z0**(n-k)
                   for n,a in enumerate(coeff) if n>=k)
    def I(x):return c.mpf(x.numerator)/x.denominator
    U=[1,2,3,1,2,1];M=[2,-1,1,2,-1,1];A=[1,-2,1,2,1,1]
    def sector(q,powers,origin):return {'powers':powers,'coefficient':I(q),'origin':origin,'source_record':'synthetic_'+origin}
    m0={k:[sector(jet(M,k),(0,1,0,-1),'m0')] for k in range(6)}
    u={k:[sector(jet(U,k),(0,0,0,0),'u')] for k in range(6)}
    slope={k:sector(jet(A,k),(0,0,0,-1),'slope') for k in range(6)}
    def actual_U(v):
        return {k:u[k]+[{**slope[k],'coefficient':slope[k]['coefficient']*v}] for k in range(6)}
    atoms={k:'synthetic_'+k for k in ('hb','F0_anchor','Pstar_squared','Lambda')}
    cert={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True}
    sc=c.mpf(1)/2048;s=sc/4
    inv_rho=c.exp(-(c.mpf(2)/5)*s)/4
    ur={k:[{**slope[k],'powers':(-1,0,0,-1),'coefficient':slope[k]['coefficient']*inv_rho}] for k in range(4)}
    out=leaf.build(c,sc=sc,phase=s,inlet_mean_Z=m0,U_Z=actual_U,U_r_Z=ur,cells=4,source_family='synthetic',source_atoms=atoms,geometry_certificate=cert)
    inlet=leaf.build(c,sc=sc,phase=c.mpf(0),inlet_mean_Z=m0,U_Z=actual_U,U_r_Z=ur,cells=1,source_family='synthetic',source_atoms=atoms,geometry_certificate=cert)
    bases=(c.mpf(2)/5,c.mpf(7)/3,c.mpf(3)/2,c.mpf(3))
    def eval_interval(terms):
        total=c.mpf(0)
        for item in terms:
            v=item['coefficient']
            for b,k in zip(bases,item['powers']):v*=b**k
            total+=v
        return total
    mp.mp.dps=80;h=mp.mpf(2)/5;phase=mp.mpf(1)/8192;rho=4*mp.exp(h*phase)
    checks={}
    for k in range(6):
        mk=mp.mpf(jet(M,k).numerator)/jet(M,k).denominator*mp.mpf(7)/9
        uk=mp.mpf(jet(U,k).numerator)/jet(U,k).denominator
        ak=mp.mpf(jet(A,k).numerator)/jet(A,k).denominator/3
        decay=mp.exp(-h*phase)
        expected=decay*mk+(1-decay)*uk+ak*(phase-(1-decay)/h)
        uk=uk+ak*phase
        actual=eval_interval(out['mean'][k]);ref=c.mpf(mp.nstr(expected,65))
        checks['mean_'+str(k)]=bool(actual.a<=ref.a and actual.b>=ref.b)
        if k<5:
            ref1=c.mpf(mp.nstr((uk-expected)/rho,65));got=eval_interval(out['mean_r'][k])
            checks['mean_r_'+str(k)]=bool(got.a<=ref1.a and got.b>=ref1.b)
        if k<4:
            ref2=c.mpf(mp.nstr((ak/(h*rho)-2*(uk-expected)/rho)/rho,65));got=eval_interval(out['mean_rr'][k])
            checks['mean_rr_'+str(k)]=bool(got.a<=ref2.a and got.b>=ref2.b)
        got=eval_interval(inlet['mean'][k]);ref=eval_interval(m0[k])
        checks['inlet_'+str(k)]=bool(got.a<=ref.a and got.b>=ref.b)
        checks['local_hb_'+str(k)]=all(t['powers'][0]>=1 for t in out['mean'][k][len(m0[k]):])
    checks['not_certified_source']=out['same_source_identity_proved'] is False
    checks['negative_hb_in_Ur_preserved']=any(t['powers'][0]==-1 for t in out['mean_rr'][0])
    try:
        leaf.build(c,sc=sc,phase=s,inlet_mean_Z=m0,U_Z=actual_U,U_r_Z=ur,cells=1,source_family='synthetic',source_atoms=atoms,geometry_certificate={})
    except ValueError:checks['reject_uncertified_geometry']=True
    else:checks['reject_uncertified_geometry']=False
    try:
        leaf.build(c,sc=sc,phase=s,inlet_mean_Z=m0,U_Z=lambda phase:{k:v for k,v in actual_U(phase).items() if k!=5},U_r_Z=ur,cells=1,source_family='synthetic',source_atoms=atoms,geometry_certificate=cert)
    except ValueError:checks['reject_missing_U_Z5']=True
    else:checks['reject_missing_U_Z5']=False
    return {'passed':all(checks.values()),'check_count':len(checks),'checks':checks,
      'synthetic_only':True,'reference':'closed form affine-phase U, own inlet mean and exact radial identities'}
