"""Independent direct-integral checks for formal-sector collar Green words."""
from __future__ import annotations


def run(leaf,c):
    import mpmath as mp
    from fractions import Fraction
    c0,c1=leaf.C0_KEYS,leaf.B1_KEYS
    def raw(k,i):
        # Nonzero signed synthetic coefficients and independently varying atoms.
        numerator=(k+2)*(i+1)*(-1 if (k+i)%2 else 1)
        return (numerator, (0,(k+i)%2,(k+i)%2,-((k+2*i)%2)))
    def sec(k,i):
        n,p=raw(k,i)
        return {'powers':p,'coefficient':c.mpf(n)/31,'origin':('synthetic',k,i)}
    def eval_iv(terms):
        total=c.mpf(0);bases=(c.mpf(1)/4,c.mpf(7)/3,c.mpf(3)/2,c.mpf(3))
        for term in terms:
            q=term['coefficient']
            for base,p in zip(bases,term['powers']):q*=base**p
            total+=q
        return total
    def val(k,i):
        n,p=raw(k,i);bases=(mp.mpf(1)/4,mp.mpf(7)/3,mp.mpf(3)/2,mp.mpf(3));q=mp.mpf(n)/31
        for b,e in zip(bases,p):q*=b**e
        return q
    core0={j:tuple([sec(j,i)] for i in range(6)) for j in range(4)}
    core1=tuple([sec(5,i)] for i in range(6))
    g={j:{str(r):[sec(j+6,r)] for r in (4,5,6)} for j in range(4)}
    matrix={name:{j:{key:[sec(j+(10 if name=='C0' else 13),i)]
      for i,key in enumerate(sorted(keys))} for j in range(3)} for name,keys in (('C0',c0),('B1',c1))}
    sc=c.mpf(1)/2048;phase=sc/4
    result=leaf.enclose(c,sc=sc,phase=phase,core_W0_Z=core0,core_W1_Z2=core1,
      forcing_Z=lambda v:g,matrix_Z=lambda v:matrix,inner_cells=2,outer_cells=2,source_family='synthetic',source_atoms={k:'synthetic_'+k for k in ('hb','F0_anchor','Pstar_squared','Lambda')},geometry_certificate={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True})
    inlet=leaf.enclose(c,sc=sc,phase=c.mpf(0),core_W0_Z=core0,core_W1_Z2=core1,
      forcing_Z=lambda v:g,matrix_Z=lambda v:matrix,inner_cells=1,outer_cells=1,source_family='synthetic',source_atoms={k:'synthetic_'+k for k in ('hb','F0_anchor','Pstar_squared','Lambda')},geometry_certificate={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True})
    mp.mp.dps=80;hb=mp.mpf(1)/4;end=mp.mpf(1)/8192
    def T(v,x):
        d=hb*v;e=mp.exp(-d);f,u,k,p,fx,ux=x
        return (f+(1-e)*fx,u+d*ux,e*k-(1-e)*ux,p,mp.exp(-3*d/2)*fx,mp.exp(-d/2)*ux)
    def W(v,j):
        X=2*mp.exp(hb*v/2)
        a=(X**2-4)/4-(X**4-16)/(8*X**2)
        b=X**2/4-2*mp.log(X/2)-1
        masses=(a,b,-a,X-2,(X**4-16)/(4*X**3),(X**2-4)/(2*X))
        t=T(v,[val(j,i) for i in range(6)])
        return [t[i]+masses[i]*val(j+6,int(row)) for i,row in enumerate(leaf.ROWS)]
    specs=(('C0',2,0,1),('C0',1,1,2),('C0',0,2,1),('B1',2,1,1),('B1',1,2,2),('B1',0,3,1))
    def source(v):
        word={j:W(v,j) for j in range(4)};out={str(r):mp.mpf(0) for r in (4,5,6)}
        for name,mo,wo,factor in specs:
            keys=c0 if name=='C0' else c1
            for i,key in enumerate(sorted(keys)):
                row,col=key.split(',');out[row]+=factor*val(mo+(10 if name=='C0' else 13),i)*word[wo][int(col)-1]
        return out
    def integrand(v,i):
        u=2*mp.exp(hb*v/2);d=hb*(end-v);ratio=mp.exp(-d/2)
        K=u*(1-mp.exp(-d))/2;kernel=(K,u*d/2,-K,mp.mpf(1),ratio**3,ratio)
        return source(v)[leaf.ROWS[i]]*kernel[i]*hb*u/2
    checks={};core=T(end,[val(5,i) for i in range(6)])
    for i,name in enumerate(leaf.FIELDS):
        reference=core[i]+mp.quad(lambda v:integrand(v,i),[0,end]);got=eval_iv(result['finite_W1_Z2'][name]);ref=c.mpf(mp.nstr(reference,65))
        checks['direct_integral_'+name]=bool(got.a<=ref.a and got.b>=ref.b)
        at=eval_iv(inlet['finite_W1_Z2'][name]);expected=eval_iv([sec(5,i)])
        checks['core_inlet_'+name]=bool(at.a<=expected.a and at.b>=expected.b)
        expected_power=2 if name in ('F1','Uz1','K1') else 1
        checks['local_hb_factor_'+name]=all(t['powers'][0]>=expected_power for t in result['local_signed_sectors'][name])
        checks['source_provenance_'+name]=bool(result['local_signed_sectors'][name]) and all('origin' in t for t in result['local_signed_sectors'][name])
    checks['no_materialized_source']=result['source_power_basis']==['hb','F0_anchor','Pstar_squared','Lambda']
    checks['not_same_source_claim']=result['source_callback_identity_proved'] is False
    return {'passed':all(checks.values()),'check_count':len(checks),'checks':checks,'synthetic_only':True,
       'reference':'independent high-precision integral with exact constant g Green masses',
       'local_sector_counts':{k:len(v) for k,v in result['local_signed_sectors'].items()}}
