"""Independent Cauchy/FTC test with actual synthetic analytic functions."""
from __future__ import annotations


def run(leaf,c):
    import sympy as S
    from fractions import Fraction
    z=S.symbols('z');F=1+z+z**2+z**5;U=2-z+z**2-z**5
    f=2+z+z**2;u=1-z+z**2
    def cv(q):q=S.Rational(q);return c.mpf(int(q.p))/int(q.q)
    def fraction(bits):
        sign,m,e,_=bits;v=Fraction(int(m))*(2**e if e>=0 else Fraction(1,2**(-e)))
        return -v if sign else v
    sc=c.mpf(1)/2048;phase=sc/4;hb=S.Rational(1,3);sp=S.Rational(1,8192)
    atoms={k:'synthetic_'+k for k in ('hb','F0_anchor','Pstar_squared','Lambda')}
    contract={k:True for k in ('source_family_checked','core_enclosure_contains_rho4','complex_radial_bound_covers_collar','analytic_Z_disks_cover_real_domain','hb_Z_independent','rho_upper_covers_collar')}
    checks={}
    for q in (S.Rational(-1),S.Rational(1,4),S.Rational(1)):
        core={name:{k:cv(S.diff(expr,z,k).subs(z,q)) for k in range(6)} for name,expr in (('Phi',F),('U',U))}
        out=leaf.build(c,sc=sc,phase=phase,core_inlet_Z=core,radial_tube_bounds={'Phi':{0:c.mpf(2)},'U':{0:c.mpf(2)}},
          cauchy_radius=c.mpf(1),rho_upper=cv(S.Rational(41,10)),source_family='synthetic',source_atoms=atoms,admitted_source_contract=contract)
        for name,expr,change in (('Phi',F,f),('U',U,u)):
            for k in range(6):
                terms=out['fields'][name][k];value=terms[0]['coefficient']+cv(hb)*terms[1]['coefficient']
                exact=S.diff(expr+hb*sp*change,z,k).subs(z,q);target=Fraction(int(exact.p),int(exact.q))
                checks[name+'_'+str(q)+'_Z'+str(k)]=fraction(value._mpi_[0])<=target<=fraction(value._mpi_[1])
                checks['error_not_value_'+name+'_'+str(q)+'_'+str(k)]=terms[1]['semantics']=='symmetric_certified_error_not_field_value' and terms[1]['powers'][0]==1
    checks['not_self_certified']=out['same_source_contract_independently_proved'] is False
    return {'passed':all(checks.values()),'check_count':len(checks),'checks':checks,
      'synthetic_only':True,'reference':'exact analytic phase evolution with explicit complex-disk radial bounds'}
