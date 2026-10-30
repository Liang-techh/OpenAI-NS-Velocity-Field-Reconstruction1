"""Independent direct-integral checks for formal-sector collar Green words."""
from __future__ import annotations


def run(leaf,c):
    import mpmath as mp
    from fractions import Fraction
    c0,c1=leaf.C0_KEYS,leaf.B1_KEYS
    expected_fields=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1');expected_rows=('5','6','6','4','5','6')
    expected_c0={'4,1','5,1','5,2','5,3','5,5','6,1','6,2','6,3','6,4','6,6'}
    expected_b1={'5,1','5,2','5,3','6,2','6,3','6,4'}
    def raw(k,i):
        # Nonzero signed synthetic coefficients and independently varying atoms.
        numerator=(k+2)*(i+1)*(-1 if (k+i)%2 else 1)
        return (numerator, (0,(k+i)%2,(k+i)%2,-((k+2*i)%2)))
    def sec(k,i):
        n,p=raw(k,i)
        return {'powers':p,'coefficient':c.mpf(n)/31,'origin':{'source_ids':['s:'+str(k)+':'+str(i)],'uncertainty':k%2==0,'error_classes':(['fixture_error'] if k%2==0 else [])}}
    def eval_iv(terms):
        total=c.mpf(0);bases=(c.mpf(1)/3,c.mpf(5)/4,c.mpf(3)/2,c.mpf(4))
        for term in terms:
            q=term['coefficient']
            for base,p in zip(bases,term['powers']):q*=base**p
            total+=q
        return total
    def val(k,i):
        n,p=raw(k,i);bases=(mp.mpf(1)/3,mp.mpf(5)/4,mp.mpf(3)/2,mp.mpf(4));q=mp.mpf(n)/31
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
      forcing_Z=lambda v:(_ for _ in ()).throw(ValueError('zero phase callback')),matrix_Z=lambda v:(_ for _ in ()).throw(ValueError('zero phase callback')),inner_cells=1,outer_cells=1,source_family='synthetic',source_atoms={k:'synthetic_'+k for k in ('hb','F0_anchor','Pstar_squared','Lambda')},geometry_certificate={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True})
    full_phase=c.make_mpf((c.mpf(0)._mpi_[0],phase.b._mpi_[1]))
    full=leaf.enclose(c,sc=sc,phase=full_phase,core_W0_Z=core0,core_W1_Z2=core1,
      forcing_Z=lambda v:g,matrix_Z=lambda v:matrix,inner_cells=1,outer_cells=1,source_family='synthetic',source_atoms={k:'synthetic_'+k for k in ('hb','F0_anchor','Pstar_squared','Lambda')},geometry_certificate={'hb_strictly_between_zero_and_half':True,'sc_phase_Z_independent':True,'hb_epsilon_Z_independent':True})
    mp.mp.dps=80;hb=mp.mpf(1)/3;end=mp.mpf(1)/8192
    def T(v,x):
        d=hb*v;e=mp.exp(-d);f,u,k,p,fx,ux=x
        return (f+(1-e)*fx,u+d*ux,e*k-(1-e)*ux,p,mp.exp(-3*d/2)*fx,mp.exp(-d/2)*ux)
    def W(v,j):
        X=2*mp.exp(hb*v/2)
        a=(X**2-4)/4-(X**4-16)/(8*X**2)
        b=X**2/4-2*mp.log(X/2)-1
        masses=(a,b,-a,X-2,(X**4-16)/(4*X**3),(X**2-4)/(2*X))
        t=T(v,[val(j,i) for i in range(6)])
        return [t[i]+masses[i]*val(j+6,int(row)) for i,row in enumerate(expected_rows)]
    def source(v):
        word={j:W(v,j) for j in range(4)};out={str(r):mp.mpf(0) for r in (4,5,6)}
        # Coefficient convolution of compatible Taylor polynomials; multiplying
        # by 2! differentiates the product, independently of the leaf specs.
        for name in ('C0','B1'):
            keys=expected_c0 if name=='C0' else expected_b1
            for i,key in enumerate(sorted(keys)):
                row,col=key.split(',');col=int(col)-1;offset=10 if name=='C0' else 13
                aa=[val(offset,i),val(offset+1,i),val(offset+2,i)/2]
                bb=([word[0][col],word[1][col],word[2][col]/2] if name=='C0' else [word[1][col],word[2][col],word[3][col]/2])
                out[row]+=2*sum(aa[j]*bb[2-j] for j in range(3))
        return out
    def integrand(v,i):
        u=2*mp.exp(hb*v/2);d=hb*(end-v);ratio=mp.exp(-d/2)
        K=u*(1-mp.exp(-d))/2;kernel=(K,u*d/2,-K,mp.mpf(1),ratio**3,ratio)
        return source(v)[expected_rows[i]]*kernel[i]*hb*u/2
    checks={};core=T(end,[val(5,i) for i in range(6)])
    for i,name in enumerate(leaf.FIELDS):
        reference=core[i]+mp.quad(lambda v:integrand(v,i),[0,end]);got=eval_iv(result['finite_W1_Z2'][name]);ref=c.mpf(mp.nstr(reference,65))
        checks['direct_integral_'+name]=bool(got.a<=ref.a and got.b>=ref.b)
        whole=eval_iv(full['finite_W1_Z2'][name]);checks['full_phase_'+name]=bool(whole.a<=ref.a and whole.b>=ref.b)
        at=eval_iv(inlet['finite_W1_Z2'][name]);expected=eval_iv([sec(5,i)])
        checks['core_inlet_'+name]=bool(at.a<=expected.a and at.b>=expected.b)
        expected_power=2 if name in ('F1','Uz1','K1') else 1
        checks['local_hb_factor_'+name]=all(t['powers'][0]>=expected_power for t in result['local_signed_sectors'][name])
        checks['source_provenance_'+name]=bool(result['local_signed_sectors'][name]) and all('origin' in t for t in result['local_signed_sectors'][name])
    checks['no_materialized_source']=result['source_power_basis']==['hb','F0_anchor','Pstar_squared','Lambda']
    checks['not_same_source_claim']=result['source_callback_identity_proved'] is False
    checks['error_class_retained']=any('fixture_error' in term['origin'].get('error_classes',[]) for vals in result['local_signed_sectors'].values() for term in vals)
    checks['fixed_field_and_row_maps']=leaf.FIELDS==expected_fields and leaf.ROWS==expected_rows
    checks['fixed_matrix_support']=c0==expected_c0 and c1==expected_b1
    return {'passed':all(checks.values()),'check_count':len(checks),'checks':checks,'synthetic_only':True,
       'reference':'independent high-precision integral with source/error groups and exact g Green masses',
       'local_sector_counts':{k:len(v) for k,v in result['local_signed_sectors'].items()}}
