"""Independent mean-source test using nonconstant actual-phase U."""
from __future__ import annotations


def run(leaf,c):
    import sympy as S
    import mpmath as mp
    z=S.symbols('z');U=2+z-z*z+z**5;M=1+2*z+z**3-z**5;G=1-z+2*z*z
    def cv(q):q=S.Rational(q);return c.mpf(int(q.p))/int(q.q)
    atoms={k:'synthetic_'+k for k in ('hb','F0_anchor','Pstar_squared','Lambda')}
    contract={k:True for k in ('same_source_checked','core_mean_inlet_actual','U_whole_phase_enclosure','hb_Z_independent','rho4_inlet_history')}
    sc=c.mpf(1)/2048;phase=sc/4;checks={};mp.mp.dps=80;h=mp.mpf(1)/3;s=mp.mpf(1)/8192
    for q in (S.Rational(-1),S.Rational(1,3),S.Rational(1)):
        core={k:cv(S.diff(M,z,k).subs(z,q)) for k in range(6)};inputs={}
        for k in range(6):
            u=cv(S.diff(U,z,k).subs(z,q));g=cv(S.diff(G,z,k).subs(z,q));cap=(abs(g)*phase).b
            err=c.make_mpf(((-cap)._mpi_[0],cap._mpi_[1]))
            inputs[k]=[{'powers':(0,0,0,0),'coefficient':u,'origin':('core_U',k),'semantics':'actual_core_interval_not_selected_value'},
              {'powers':(1,0,0,0),'coefficient':err,'origin':('phase_error',k),'semantics':'symmetric_certified_error_not_field_value'}]
        out=leaf.build(c,sc=sc,phase=phase,core_mean_Z=core,U_Z=inputs,source_family='synthetic',source_atoms=atoms,admitted_source_contract=contract)
        inlet=leaf.build(c,sc=sc,phase=c.mpf(0),core_mean_Z=core,U_Z=inputs,source_family='synthetic',source_atoms=atoms,admitted_source_contract=contract)
        for k in range(6):
            def mpq(expr):
                v=S.Rational(S.diff(expr,z,k).subs(z,q));return mp.mpf(int(v.p))/int(v.q)
            e=mp.exp(-h*s);actual=e*mpq(M)+(1-e)*mpq(U)+(h*s-(1-e))*mpq(G)
            total=c.mpf(0)
            for t in out['mean_Z'][k]['all_terms']:total+=t['coefficient']*(c.mpf(1)/3)**t['powers'][0]
            ref=c.mpf(mp.nstr(actual,65));checks['mean_'+str(q)+'_'+str(k)]=bool(total.a<=ref.a and total.b>=ref.b)
            zterms=inlet['mean_Z'][k]['all_terms'];v=c.mpf(0)
            for t in zterms:v+=t['coefficient']*(c.mpf(1)/3)**t['powers'][0]
            checks['inlet_'+str(q)+'_'+str(k)]=bool(v.a<=core[k].a and v.b>=core[k].b)
            errors=out['mean_Z'][k]['phase_uncertainty_terms']
            checks['separate_error_'+str(q)+'_'+str(k)]=len(errors)==1 and errors[0]['powers'][0]==2 and errors[0]['semantics']=='symmetric_certified_error_not_field_value'
    checks['not_self_certified']=out['same_source_contract_independently_proved'] is False
    return {'passed':all(checks.values()),'check_count':len(checks),'checks':checks,'synthetic_only':True,
      'reference':'nonconstant phase U and exact integrated own-mean ODE'}
